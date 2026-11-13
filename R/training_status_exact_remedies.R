# Exact training-status holdout remedies ----

# Match failed or flagged holdouts to the reviewed remedy ledger and check the
# resulting corrections before they enter the final LOO table.

# Reviewed retries of failed exact checkpoints; passing cases are immutable.
tst_exact_remedy_plan <- function(paths, ledger_path) {
  ledger <- iad_read_csv(ledger_path)
  keys <- c("family", "outcome", "model_index", "observation")
  required <- c(keys, "stage", "iter_warmup", "iter_sampling", "adapt_delta", "max_treedepth", "seed", "initial_sha256", "rationale")
  if (!all(required %in% names(ledger)) || anyNA(ledger[required]) || anyDuplicated(ledger[keys]) ||
      any(!ledger$family %in% iad_families()) || any(!ledger$outcome %in% iad_outcomes()) ||
      any(!grepl("^remedy_[1-9][0-9]*$", ledger$stage)) ||
      any(!ledger$iter_sampling %in% c(6000L, 12000L)) ||
      any(ledger$iter_warmup < 2000L | ledger$iter_warmup %% 1 != 0) ||
      any(ledger$adapt_delta < 0.999 | ledger$adapt_delta >= 1) ||
      any(ledger$max_treedepth < 15L | ledger$max_treedepth %% 1 != 0) ||
      any(ledger$seed < 1 | ledger$seed %% 1 != 0) ||
      any(!grepl("^[a-f0-9]{64}$", ledger$initial_sha256)) || any(!nzchar(trimws(ledger$rationale))))
    stop("Invalid training exact-remedy ledger.")

  matched <- integer(); seen <- character()
  plans <- lapply(unname(unlist(paths)), function(path) {
    initial <- qs::qread(path)
    if (!isTRUE(initial$compact) || isTRUE(initial$skip)) stop("Expected a compact exact case.")
    record <- initial$interaction_record; d <- initial$diagnostics
    if (nrow(d) != 1L || !is.logical(d$diagnostic_pass) || is.na(d$diagnostic_pass) ||
        d$family != record$family || d$outcome != record$outcome ||
        d$model_index != record$model_index || d$observation != record$observation)
      stop("Exact diagnostic identity mismatch.")
    seen <<- c(seen, paste(record$family, record$outcome, record$model_index, record$observation, sep = "|"))
    index <- which(ledger$family == record$family & ledger$outcome == record$outcome &
      ledger$model_index == record$model_index & ledger$observation == record$observation)
    remedy <- NULL
    if (isTRUE(d$diagnostic_pass)) {
      if (length(index) || is.null(initial$correction)) stop("Passing exact case must retain its score and must not be refitted.")
    } else {
      if (length(index) != 1L) stop("Every failed exact case needs one reviewed remedy.")
      remedy <- as.list(ledger[index, , drop = FALSE])
      iad_check_hash(path, remedy$initial_sha256)
      if (remedy$iter_sampling < d$iter_sampling || remedy$iter_warmup < d$iter_warmup)
        stop("Exact remedy cannot shorten the accepted schedule.")
      matched <<- c(matched, index)
    }
    list(initial_path = path, initial_sha256 = iad_hash(path), record = record, remedy = remedy)
  })

  if (anyDuplicated(seen) || length(matched) != nrow(ledger) || !setequal(matched, seq_len(nrow(ledger))))
    stop("Exact-remedy ledger does not match the failed cases uniquely.")
  plans
}

tst_review_exact <- function(plan) {
  iad_require_phase("exact_remedy")
  iad_check_hash(plan$initial_path, plan$initial_sha256)
  if (is.null(plan$remedy)) return(plan$initial_path)
  r <- plan$record; remedy <- plan$remedy
  path <- iad_path(r$family, r$outcome, "loo", "exact_checkpoints", remedy$stage,
    sprintf("model_%02d_observation_%04d.qs", r$model_index, r$observation), archive = TRUE)
  iad_check_hash(r$fit_path, r$fit_sha256); iad_check_hash(r$loo_path, r$loo_sha256)

  if (file.exists(path)) {
    previous <- qs::qread(path)
    if (!isTRUE(previous$compact) || !identical(previous$interaction_record, r) ||
        !identical(previous$exact_remedy, plan)) stop("Existing exact retry has different settings or inputs.")
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  candidate <- qs::qread(r$fit_path); integrated <- qs::qread(r$loo_path)
  if (!isTRUE(loo::pareto_k_values(integrated$loo)[[r$observation]] > 0.7)) stop("Invalid retry observation.")
  s <- candidate$specification
  if (remedy$adapt_delta < s$adapt_delta || remedy$max_treedepth < s$max_treedepth)
    stop("Exact remedy cannot weaken inherited sampler controls.")
  # Only sampler controls change in the copied recipe. The brms formula and
  # prior remain those of the accepted candidate and update(recompile=FALSE).
  candidate$specification$adapt_delta <- remedy$adapt_delta
  candidate$specification$max_treedepth <- as.integer(remedy$max_treedepth)
  refit <- fit_primary_integrated_exact_loo_case(c(r, list(fitted_candidate = candidate, integrated_result = integrated)),
    chains = 4L, cores = 4L, warmup = as.integer(remedy$iter_warmup),
    iter = as.integer(remedy$iter_warmup + remedy$iter_sampling),
    iter_sampling = as.integer(remedy$iter_sampling), seed = as.integer(remedy$seed))
  refit$interaction_record <- r
  rm(candidate, integrated)
  diagnostics <- iad_exact_diagnostic_summary(refit)
  diagnostics$exact_stage <- remedy$stage
  diagnostics$adapt_delta <- remedy$adapt_delta
  diagnostics$max_treedepth <- remedy$max_treedepth
  diagnostics$initial_checkpoint_sha256 <- plan$initial_sha256
  correction <- if (isTRUE(diagnostics$diagnostic_pass)) iad_build_exact_correction(refit) else NULL

  if (!is.null(correction)) {
    correction$exact_remedy <- plan
    correction$exact_refit_summary$exact_stage <- remedy$stage
    correction$exact_refit_summary$iter_sampling <- remedy$iter_sampling
    correction$exact_refit_summary$iter_warmup <- remedy$iter_warmup
  }

  iad_qsave(list(compact = TRUE, interaction_record = r, exact_remedy = plan,
    diagnostics = diagnostics, correction = correction), path)
}
