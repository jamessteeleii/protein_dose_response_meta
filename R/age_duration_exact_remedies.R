# Reviewed convergence remedies for compact interaction exact-LOO checkpoints.
# Initial checkpoints remain immutable; passing cases keep their existing paths.
iad_exact_remedy_plan <- function(paths, ledger_path, family) {
  ledger <- iad_read_csv(ledger_path)
  keys <- c("family", "outcome", "model_index", "observation")
  required <- c(keys, "stage", "iter_sampling", "initial_sha256", "rationale")
  if (!all(required %in% names(ledger)) || anyNA(ledger[required]) ||
      anyDuplicated(ledger[keys]) || !all(ledger$family %in% iad_families()) ||
      !all(ledger$outcome %in% iad_outcomes()) ||
      !all(ledger$model_index > 0 & ledger$model_index %% 1 == 0) ||
      !all(ledger$observation > 0 & ledger$observation %% 1 == 0) ||
      !all(grepl("^remedy_[1-9][0-9]*$", ledger$stage)) ||
      !all(ledger$iter_sampling == 12000L) ||
      !all(grepl("^[a-f0-9]{64}$", ledger$initial_sha256)) ||
      !all(nzchar(trimws(ledger$rationale)))) stop("Invalid exact-refit remedy ledger.")
  ledger <- ledger[ledger$family == family, , drop = FALSE]
  matched <- integer()
  plans <- lapply(unname(unlist(paths)), function(path) {
    initial <- qs::qread(path)
    if (isTRUE(initial$skip)) {
      if (nrow(ledger)) stop("Remedy ledger has no matching exact cases.")
      return(list(initial_path = path, initial_sha256 = iad_hash(path), skip = TRUE))
    }
    record <- initial$interaction_record
    diagnostic <- iad_exact_diagnostics(path)
    if (!identical(record$family, family) || nrow(diagnostic) != 1L ||
        !identical(as.integer(diagnostic$model_index), as.integer(record$model_index)) ||
        !identical(as.integer(diagnostic$observation), as.integer(record$observation)) ||
        !identical(diagnostic$outcome, record$outcome) ||
        !is.logical(diagnostic$diagnostic_pass) || is.na(diagnostic$diagnostic_pass)) {
      stop("Invalid initial exact checkpoint identity or diagnostics.")
    }
    index <- which(ledger$outcome == record$outcome & ledger$model_index == record$model_index &
      ledger$observation == record$observation)
    remedy <- NULL
    if (isTRUE(diagnostic$diagnostic_pass)) {
      if (length(index)) stop("A passing exact case must not be refitted.")
    } else {
      if (length(index) != 1L || !isTRUE(initial$compact) ||
          !isTRUE(diagnostic$geometry_pass) || !identical(diagnostic$convergence_pass, FALSE)) {
        stop("Failed exact case needs an explicit, applicable convergence remedy.")
      }
      remedy <- as.list(ledger[index, , drop = FALSE])
      iad_check_hash(path, remedy$initial_sha256)
      matched <<- c(matched, index)
    }
    list(initial_path = path, initial_sha256 = iad_hash(path), record = record, remedy = remedy)
  })
  if (length(matched) != nrow(ledger) || !setequal(matched, seq_len(nrow(ledger)))) {
    stop("Exact remedy ledger does not match the failed cases uniquely.")
  }
  plans
}

iad_review_exact <- function(plan) {
  iad_require_phase("exact_remedy")
  iad_check_hash(plan$initial_path, plan$initial_sha256)
  if (isTRUE(plan$skip) || is.null(plan$remedy)) return(plan$initial_path)
  record <- plan$record
  remedy <- plan$remedy
  path <- iad_path(record$family, record$outcome, "loo", "exact_checkpoints", remedy$stage,
    sprintf("model_%02d_observation_%04d.qs", record$model_index, record$observation), archive = TRUE)
  iad_check_hash(record$fit_path, record$fit_sha256)
  iad_check_hash(record$loo_path, record$loo_sha256)
  if (file.exists(path)) {
    previous <- qs::qread(path)
    if (!isTRUE(previous$compact) || !identical(previous$interaction_record, record) ||
        !identical(previous$exact_remedy, plan)) stop("Existing exact remedy has different inputs or settings.")
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  candidate <- qs::qread(record$fit_path)
  integrated <- qs::qread(record$loo_path)
  if (!isTRUE(loo::pareto_k_values(integrated$loo)[[record$observation]] > 0.7)) stop("Invalid exact remedy case.")
  refit <- fit_primary_integrated_exact_loo_case(c(record,
    list(fitted_candidate = candidate, integrated_result = integrated)),
    chains = 4L, cores = 4L, warmup = 2000L, iter = 2000L + as.integer(remedy$iter_sampling),
    seed = 1988L, iter_sampling = as.integer(remedy$iter_sampling))
  refit$interaction_record <- record
  rm(candidate, integrated)
  diagnostics <- iad_exact_diagnostic_summary(refit)
  diagnostics$exact_stage <- remedy$stage
  diagnostics$iter_sampling <- as.integer(remedy$iter_sampling)
  diagnostics$initial_checkpoint_sha256 <- plan$initial_sha256
  correction <- if (isTRUE(diagnostics$diagnostic_pass)) iad_build_exact_correction(refit) else NULL
  if (!is.null(correction)) {
    correction$exact_remedy <- plan
    correction$exact_refit_summary$exact_stage <- remedy$stage
    correction$exact_refit_summary$iter_sampling <- as.integer(remedy$iter_sampling)
    correction$exact_refit_summary$initial_checkpoint_sha256 <- plan$initial_sha256
  }
  iad_qsave(list(compact = TRUE, interaction_record = record, exact_remedy = plan,
    diagnostics = diagnostics, correction = correction), path)
}
