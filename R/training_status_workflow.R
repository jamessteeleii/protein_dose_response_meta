# Source AFTER age_duration_interactions.R, only in _targets_training_status.R.
# Reuse the reviewed diagnostics/exact-LOO machinery in a separate R process.
# These adapters deliberately replace iad_* entry points in that process only.
tst_original_read <- iad_read_checkpoints
tst_original_fingerprint <- iad_fingerprint

# Exact refits inherit the accepted full-fit sampling length. This matters for
# supplementation-strength model 14, whose accepted remedy retains 12000 draws.
tst_original_exact_fit <- fit_primary_integrated_exact_loo_case

fit_primary_integrated_exact_loo_case <- function(job, chains = 4L, cores = 4L,
  warmup = 2000L, iter = NULL, seed = 1988L, iter_sampling = NULL) {
  if (is.null(iter_sampling)) iter_sampling <- job$fitted_candidate$specification$iter_sampling
  if (length(iter_sampling) != 1L || is.na(iter_sampling) || !iter_sampling %in% c(6000L, 12000L))
    stop("Missing or invalid accepted exact-refit sampling schedule.")
  if (is.null(iter)) iter <- warmup + iter_sampling
  result <- tst_original_exact_fit(job, chains = chains, cores = cores, warmup = warmup,
    iter = iter, seed = seed, iter_sampling = as.integer(iter_sampling))
  result$training_exact_schedule <- list(chains = chains, iter_warmup = warmup,
    iter_sampling = as.integer(iter_sampling), seed = seed)
  result
}
tst_original_exact_diagnostics <- iad_exact_diagnostic_summary

iad_exact_diagnostic_summary <- function(refit) {
  result <- tst_original_exact_diagnostics(refit)
  result$iter_sampling <- refit$training_exact_schedule$iter_sampling
  result$iter_warmup <- refit$training_exact_schedule$iter_warmup
  result
}

iad_path <- function(..., archive = FALSE) here::here(
  if (archive) "model_archive" else "results", "update_2026", "training_status", ...)

iad_phase <- function(phase) {
  value <- Sys.getenv("PROTEIN_TRAINING_PHASE", "none")
  if (!value %in% c("none", "pilot", "fit", "remedy", "loo", "exact", "exact_remedy", "exact_marginal", "correction", "combine", "report"))
    stop("Unknown PROTEIN_TRAINING_PHASE.")
  identical(value, phase)
}

# Frozen inputs for the training-status comparison ----

# Check the earlier interaction reporting interface and read the primary
# parent fits for reuse in the separate training-status workflow.
tst_age_duration_gate <- function(root = here::here("results", "update_2026", "age_duration"),
  interface = here::here("manuscript", "generated", "age-duration-interactions", "_results.qmd")) {
  # Twelve completed reporting manifests, plus the final materialised interface.
  for (family in iad_families()) for (outcome in iad_outcomes()) {
    path <- file.path(root, family, outcome, "reporting", "file_manifest.csv")
    if (!file.exists(path)) stop("Finish age/duration reporting before training-status execution: ", path)
    manifest <- iad_read_csv(path)
    if (!nrow(manifest) || anyDuplicated(manifest$path)) stop("Invalid age/duration reporting manifest.")
    for (i in seq_len(nrow(manifest))) iad_check_hash(manifest$path[[i]], manifest$sha256[[i]])
  }
  path <- interface
  if (!file.exists(path) || !any(grepl("All twelve primary-family expansions passed", readLines(path))))
    stop("The final age/duration supplement gate has not finished.")
  invisible(TRUE)
}

iad_require_phase <- function(phase) {
  if (!iad_phase(phase)) stop("Training-status phase disabled: ", phase)
  tst_age_duration_gate()
  iad_require_measurement()
}

iad_fingerprint <- function(data) digest::digest(list(
  original = tst_original_fingerprint(data), status = as.character(data$training_status),
  code = data$training_status_code), algo = "sha256")

iad_read_checkpoints <- function(paths) {
  checkpoints <- tst_original_read(paths)
  for (family in iad_families()) for (outcome in iad_outcomes()) {
    source <- checkpoints[[family]]$data[[outcome]]
    status <- tolower(trimws(as.character(source$trained_or_untrained)))
    status[is.na(status) | !nzchar(status)] <- "not reported"
    if (any(!status %in% c("trained", "untrained", "active", "not reported"))) stop("Unreviewed status label.")
    if (any(source$resistance_exercise_code != ifelse(source$resistance_exercise_yes_no == "yes", 0, 1)))
      stop("RT code disagrees with source label.")
    keep <- source$resistance_exercise_code == 0 & status %in% c("trained", "untrained")
    data <- source[keep, , drop = FALSE]
    data$training_status <- factor(status[keep], levels = c("trained", "untrained"))
    data$training_status_code <- as.integer(data$training_status == "untrained")
    data$study_factor <- factor(as.character(data$study_no))
    if (!setequal(data$training_status_code, 0:1)) stop("Both status groups are required.")
    by_arm <- split(as.character(data$training_status), interaction(data$study_no, data$arm_no, drop = TRUE))
    if (any(vapply(by_arm, function(x) length(unique(x)), integer(1)) != 1L)) stop("Arm status is inconsistent.")
    checkpoints[[family]]$data[[outcome]] <- data
    counts <- checkpoints[[family]]$counts
    at <- counts$outcome == outcome
    counts$source_rows[at] <- nrow(source)
    counts$common_rows[at] <- nrow(data)
    counts$studies[at] <- dplyr::n_distinct(data$study_no)
    counts$arms[at] <- dplyr::n_distinct(data$arm_no)
    counts$data_sha256[at] <- iad_fingerprint(data)
    checkpoints[[family]]$counts <- counts
  }

  checkpoints
}

# Training-status candidate preparation ----

# Build candidate specifications from the accepted parent models and retain
# the indexing needed to compare the parent and interaction scores.
iad_specifications <- function(checkpoints) {
  records <- list()
  for (family in iad_families()) for (outcome in iad_outcomes()) {
    table <- checkpoints[[family]]$specifications
    table <- table[table$outcome == outcome, ]
    table <- table[order(table$model_index), ]
    for (i in seq_len(nrow(table))) {
      s <- as.list(table[i, ])
      s$family <- family; s$parent_index <- as.integer(s$model_index)
      s$model_index <- as.integer(s$model_index)
      s$parent_formula <- s$formula
      s$formula <- gsub("resistance_exercise_factor", "training_status", s$formula, fixed = TRUE)
      s$formula <- gsub("resistance_exercise_code", "training_status_code", s$formula, fixed = TRUE)
      s$model_name <- paste0(s$model_name, " [RT only; reported training status]")
      s$interaction <- if (s$model_index >= 15L) "training" else "none"
      s$stage <- "initial"; s$chains <- 4L; s$parallel_chains <- 4L
      s$iter_warmup <- 2000L; s$iter_sampling <- 6000L
      s$adapt_delta <- max(as.numeric(s$adapt_delta), 0.999)
      s$max_treedepth <- max(as.integer(s$max_treedepth), 15L)
      s$accepted_fit_path <- s$accepted_fit_stage <- NULL
      records[[length(records) + 1L]] <- s
    }
  }

  records
}

iad_prepare_input <- function(s, checkpoints, priors) {
  cp <- checkpoints[[s$family]]
  parent <- as.list(cp$specifications[cp$specifications$outcome == s$outcome &
    cp$specifications$model_index == s$parent_index, ])
  prepared <- if (s$family == "relative") prepare_update_2026_primary_candidate_input(parent, cp$data, priors) else
    prepare_update_2026_total_primary_input(parent, cp$data, priors)
  data <- add_primary_half_unit_protein(add_primary_standardised_moderators(prepared$data))

  if (s$family != "relative") data <- add_absolute_primary_body_mass(data)
  data$Intercept <- 1
  if (!iad_same_observations(data, cp$data[[s$outcome]])) stop("Status preparation dropped observations.")
  # Existing RT priors do not describe training status. Keep other applicable
  # parent priors; give every new status coefficient the stated regularising prior.
  p <- prepared$priors[!grepl("resistance_exercise", prepared$priors$coef), ]
  available <- brms::get_prior(stats::as.formula(s$formula), data = data)
  key <- function(x) paste(x$class, x$coef, x$group, x$resp, x$dpar, x$nlpar, sep = "|")
  p <- p[key(p) %in% key(available), ]
  added <- available$coef[available$class == "b" & grepl("training_status", available$coef)]
  for (coef in unique(added)) p <- c(p, brms::set_prior("student_t(3, 0, 0.5)", class = "b", coef = coef))

  if (any(available$class == "sds") && !any(p$class == "sds"))
    p <- c(p, brms::set_prior("normal(0, 0.5)", class = "sds"))
  list(formula = s$formula, data = data, priors = p, specification = s)
}

iad_pilot_specs <- function(specifications) Filter(function(s) s$model_index %in% c(6L, 26L), specifications)
# The shared fitter now calls a family pilot gate. This adapter validates the
# RT/status pilots rather than the age/duration joint-interaction pilot IDs.
iad_check_family_pilot <- function(pilot, family, checkpoints) {
  expected <- dplyr::bind_rows(iad_pilot_specs(iad_specifications(checkpoints)))
  expected <- expected[expected$family == family, c("family", "outcome", "model_index")]
  if (is.null(pilot) || !all(c(names(expected), "diagnostic_pass") %in% names(pilot)))
    stop("Training-status pilot diagnostics are missing.")
  selected <- pilot[pilot$family == family, , drop = FALSE]
  if (nrow(expected) != 8L || nrow(selected) != 8L || anyNA(selected$diagnostic_pass) ||
      !all(selected$diagnostic_pass) || anyDuplicated(selected[names(expected)]) ||
      nrow(dplyr::anti_join(expected, selected, by = names(expected))))
    stop("All eight training-status pilots in the family must pass before remaining fits.")
  invisible(TRUE)
}
# Every model is refitted; no full-sample posterior is a combination candidate.
iad_require_parents <- function(checkpoints) invisible(TRUE)
tst_parent_fit <- iad_fit

iad_fit <- function(specification, checkpoints, priors, preflight, phase, pilot = NULL) {
  tst_parent_fit(specification, checkpoints, priors, preflight, phase, pilot)
}

tst_status_counts <- function(checkpoints) purrr::map_dfr(iad_families(), function(family)
  purrr::map_dfr(iad_outcomes(), function(outcome) {
    data <- checkpoints[[family]]$data[[outcome]]
    data |> dplyr::group_by(training_status) |> dplyr::summarise(
      studies = dplyr::n_distinct(study_no), arms = dplyr::n_distinct(arm_no), observations = dplyr::n(),
      distinct_doses = dplyr::n_distinct(protein_intake_raw), dose_min = min(protein_intake_raw),
      dose_max = max(protein_intake_raw), .groups = "drop") |>
      dplyr::mutate(family = family, outcome = outcome, .before = 1L)
  }))

tst_prepare <- function(checkpoints, specifications, preflight) {
  if (nrow(preflight) != 248L || !all(preflight$valid)) stop("Expected 248 passing inputs.")
  status <- tst_status_counts(checkpoints)
  tables <- list(counts = dplyr::bind_rows(lapply(checkpoints, `[[`, "counts")),
    status_counts = status, specifications = dplyr::bind_rows(specifications), preflight = preflight,
    source_manifest = dplyr::bind_rows(lapply(checkpoints, `[[`, "source_manifest")),
    storage_estimate = iad_storage_estimate(checkpoints, specifications))
  unname(purrr::imap_chr(tables, function(x, name) iad_csv(x, iad_path("preparation", paste0(name, ".csv")))))
}

# Model comparison and prediction ----

# Combine parent and training-status interaction scores for the same
# observations, then use the combined weights in the reporting mixture.
iad_combine <- function(checkpoints, accepted, final_paths) {
  iad_require_phase("combine")
  accepted <- dplyr::bind_rows(accepted)
  results <- lapply(unname(unlist(final_paths)), qs::qread)
  if (nrow(accepted) != 248L || length(results) != 248L || !all(accepted$diagnostic_pass)) stop("Incomplete status family.")
  output <- list()

  for (family in iad_families()) for (outcome in iad_outcomes()) {
    models <- accepted[accepted$family == family & accepted$outcome == outcome, ]
    models <- models[order(models$model_index), ]
    expected <- checkpoints[[family]]$specifications$model_index[checkpoints[[family]]$specifications$outcome == outcome]
    if (anyDuplicated(models$model_index) || !setequal(models$model_index, expected)) stop("Missing status candidate.")
    selected <- Filter(function(x) x$family == family && x$outcome == outcome, results)
    selected <- selected[order(vapply(selected, `[[`, integer(1), "model_index"))]
    if (!identical(vapply(selected, `[[`, integer(1), "model_index"), models$model_index)) stop("LOO identity mismatch.")
    for (i in seq_along(selected)) {
      iad_validate_loo(selected[[i]], checkpoints[[family]]$data[[outcome]])
      iad_check_hash(models$fit_path[[i]], models$sha256[[i]])
      if (selected[[i]]$fit_sha256 != models$sha256[[i]]) stop("LOO belongs to a different fit.")
    }
    objects <- lapply(selected, `[[`, "loo"); names(objects) <- models$model_name
    models$stacking_weight <- iad_stacking_weights(objects)
    models$fit_sha256 <- models$sha256
    output[[paste(family, outcome, sep = "_")]] <- list(family = family, outcome = outcome, models = models,
      loo_summary = purrr::map_dfr(selected, function(r) dplyr::mutate(
        summarise_primary_loo_method(r$loo, "integrated_observation_level_psis_plus_exact"),
        model_index = r$model_index, model_name = r$model_name, .before = 1L)),
      seed = 651988L + 100L * match(family, iad_families()) + match(outcome, iad_outcomes()), gate_passed = TRUE)
  }

  output
}

tst_grid <- function(data, family) {
  # Dose centring/reference are inherited from the frozen exposure definition.
  dose <- sort(unique(c(seq(if (family == "relative") 0.8 else if (family == "total") 50 else 0,
    max(data$protein_intake_raw), by = if (family == "relative") 0.01 else 1),
    if (family == "relative") get_primary_protein_median(data) else unique(data$absolute_exposure_reference_g_day))))
  profiles <- data.frame(profile = c("reference", "typical"), age = c(25,
    median(dplyr::distinct(data, study_no, arm_no, age_centre)$age_centre + 25)),
    weeks = c(12, median(dplyr::distinct(data, study_no, time2)$time2)))
  keys <- expand.grid(dose = dose, training_status = c("trained", "untrained"), profile = profiles$profile,
    stringsAsFactors = FALSE)
  newdata <- data[rep(1L, nrow(keys)), ]; rownames(newdata) <- NULL
  newdata$training_status <- factor(keys$training_status, levels = c("trained", "untrained"))
  newdata$training_status_code <- as.integer(keys$training_status == "untrained")
  reference <- if (family == "relative") get_primary_protein_median(data) else unique(data$absolute_exposure_reference_g_day)

  if (length(reference) != 1L) stop("Ambiguous dose reference.")
  newdata$protein_intake_raw <- keys$dose
  unit <- if (family == "relative") 1 else unique(data$absolute_exposure_unit_g_day)
  if (length(unit) != 1L || !is.finite(unit) || unit <= 0) stop("Invalid exposure scaling unit.")
  newdata$protein_intake <- (keys$dose - reference) / unit
  newdata$protein_intake_05 <- newdata$protein_intake / 0.5
  scaled <- add_primary_standardised_moderators(data)
  centres <- attr(scaled, "moderator_centres"); scales <- attr(scaled, "moderator_scales")
  at <- match(keys$profile, profiles$profile)
  newdata$age_centre <- profiles$age[at] - 25; newdata$duration_centre <- profiles$weeks[at] - 12

  for (v in c("age", "duration")) newdata[[paste0(v, "_outcome_z")]] <-
    (newdata[[paste0(v, "_centre")]] - centres[[paste0(v, "_centre")]]) / scales[[paste0(v, "_centre")]]
  newdata$body_mass_outcome_z <- 0
  list(keys = keys, newdata = newdata, profiles = profiles, reference = reference)
}

tst_population_epred <- function(model, newdata, draw_ids) {
  labels <- attr(stats::terms(model$formula$formula), "term.labels")
  study_smooth <- labels[vapply(labels, function(term) "study_factor" %in% all.vars(str2lang(term)), logical(1))]
  prediction <- as.matrix(brms::posterior_epred(model, newdata = newdata, re_formula = NA, draw_ids = draw_ids))
  if (length(study_smooth)) {
    # brms excludes ordinary group terms via re_formula, but the factor
    # smooth remains an additive spline on the Gaussian identity scale.
    likelihood <- stats::family(model)
    if (length(study_smooth) != 1L || likelihood$family != "gaussian" || likelihood$link != "identity" ||
        !identical(str2lang(study_smooth)[["bs"]], "fs")) stop("Unexpected study smooth in population prediction.")
    label <- gsub("[[:space:]]+", "", study_smooth)
    deviation <- as.matrix(brms::posterior_smooths(model, smooth = label, newdata = newdata, draw_ids = draw_ids))
    if (!identical(dim(prediction), dim(deviation))) stop("Study smooth prediction dimensions differ.")
    prediction <- prediction - deviation
  }

  if (!identical(dim(prediction), c(length(draw_ids), nrow(newdata))) || any(!is.finite(prediction)))
    stop("Invalid population predictions.")
  prediction
}

# Generate and summarize dose curves by exercise training status, including
# contrasts between the status-specific increments.
tst_report <- function(record, checkpoints) {
  iad_require_phase("report")
  if (!isTRUE(record$gate_passed)) stop("Combination gate failed.")
  family <- record$family; outcome <- record$outcome; models <- record$models
  data <- checkpoints[[family]]$data[[outcome]]; grid <- tst_grid(data, family)
  allocation <- iad_model_allocation(models, 4000L, record$seed)
  draws <- matrix(NA_real_, 4000L, nrow(grid$keys))

  for (i in seq_len(nrow(models))) {
    at <- which(allocation == i); if (!length(at)) next
    iad_check_hash(models$fit_path[[i]], models$fit_sha256[[i]])
    fit <- qs::qread(models$fit_path[[i]])
    if (iad_fingerprint(fit$model$data) != iad_fingerprint(data) || fit$model_index != models$model_index[[i]] ||
        fit$outcome != outcome || fit$model_name != models$model_name[[i]]) stop("Prediction identity mismatch.")
    set.seed(record$seed + i)
    ids <- sample.int(posterior::ndraws(fit$model), length(at))
    draws[at, ] <- tst_population_epred(fit$model, grid$newdata, ids)
    rm(fit); gc(verbose = FALSE)
  }

  if (any(!is.finite(draws))) stop("Incomplete posterior mixture.")
  increments <- draws
  for (profile in grid$profiles$profile) for (status in c("trained", "untrained")) {
    cols <- which(grid$keys$profile == profile & grid$keys$training_status == status)
    at <- cols[abs(grid$keys$dose[cols] - grid$reference) < 1e-10]
    if (length(at) != 1L) stop("Missing paired reference dose.")
    increments[, cols] <- sweep(draws[, cols, drop = FALSE], 1L, draws[, at], "-")
  }

  support <- tst_status_counts(checkpoints)
  support <- support[support$family == family & support$outcome == outcome, ]
  doses <- dplyr::bind_cols(grid$keys, iad_summarise_draws(draws))
  changes <- dplyr::bind_cols(grid$keys, iad_summarise_draws(increments))
  moderation <- purrr::map_dfr(grid$profiles$profile, function(profile) {
    trained <- which(grid$keys$profile == profile & grid$keys$training_status == "trained")
    untrained <- which(grid$keys$profile == profile & grid$keys$training_status == "untrained")
    if (!identical(grid$keys$dose[trained], grid$keys$dose[untrained])) stop("Unpaired status grid.")
    dplyr::bind_cols(grid$keys[trained, c("dose", "profile")],
      iad_summarise_draws(increments[, trained, drop = FALSE] - increments[, untrained, drop = FALSE]))
  })
  moderation$current_dose <- grid$reference
  moderation$supported <- all(support$studies >= 3L) & moderation$dose >= max(support$dose_min) &
    moderation$dose <= min(support$dose_max) & grid$reference >= max(support$dose_min) & grid$reference <= min(support$dose_max)
  moderation$interpretation <- if (outcome == "function") "Sparse trained evidence; descriptive and potentially prior-sensitive" else
    "Exploratory study-level association; range screen does not establish joint support"
  public <- here::here("manuscript", "generated", "training-status", family, outcome)
  tables <- list(dose_response = doses, protein_contrasts = changes, training_interaction_contrasts = moderation,
    status_support = support, profiles = grid$profiles,
    stacking_weights = models[, c("family", "outcome", "model_index", "model_name", "interaction", "stacking_weight")])
  paths <- unname(purrr::imap_chr(tables, function(x, name) iad_csv(x, file.path(public, paste0(name, ".csv")))))
  figure <- ggplot2::ggplot(changes, ggplot2::aes(dose, estimate, colour = training_status, fill = training_status)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf_low, ymax = conf_high), alpha = 0.12, colour = NA) +
    ggplot2::geom_line() + ggplot2::facet_wrap(~profile) + ggplot2::theme_bw() +
    ggplot2::labs(x = if (family == "relative") "Protein intake (g/kg/day)" else "Protein dose (g/day)",
      y = "Change relative to reference protein dose (baseline SD)", colour = "Reported status", fill = "Reported status",
      caption = "RT arms; population means. Curves extend beyond shared empirical support; see support tables.")
  figure_path <- file.path(public, "protein_contrasts.png")
  figure_temporary <- paste0(figure_path, ".tmp-", Sys.getpid(), ".png")
  on.exit(unlink(figure_temporary), add = TRUE)
  ggplot2::ggsave(figure_temporary, figure, width = 9, height = 5, dpi = 180)

  if (!file.rename(figure_temporary, figure_path)) stop("Could not finalise training-status figure.")
  paths <- c(paths, figure_path)
  artifacts <- c(tables, list(
    verification_checks = tibble::tibble(check = c("complete family passed gates", "4000 coherent draws", "finite predictions"),
      passed = c(record$gate_passed, length(allocation) == 4000L, all(is.finite(draws)))),
    source_manifest = checkpoints[[family]]$source_manifest,
    specifications = dplyr::bind_rows(iad_specifications(checkpoints)) |>
      dplyr::filter(.data$family == .env$family, .data$outcome == .env$outcome),
    accepted_diagnostics = models[, setdiff(names(models), c("fit_path", "rank"))],
    final_loo_summary = record$loo_summary,
    draw_allocation = tibble::tibble(model_index = models$model_index, draws = tabulate(allocation, nrow(models)))))

  paths <- c(paths, write_verification_bundle(paste0("training_status_", family, "_", outcome, "_update_2026"),
    artifacts, descriptions = stats::setNames(names(artifacts), names(artifacts))))
  manifest <- tibble::tibble(path = paths, bytes = file.info(paths)$size, sha256 = vapply(paths, iad_hash, character(1)))
  c(paths, iad_csv(manifest, iad_path(family, outcome, "reporting", "file_manifest.csv")))
}

tst_supplement <- function(paths, records) {
  iad_require_phase("report")
  manifests <- unname(unlist(paths)); manifests <- manifests[basename(manifests) == "file_manifest.csv"]
  if (length(records) != 12L || length(manifests) != 12L) stop("All twelve scopes must finish.")
  for (path in manifests) {
    m <- iad_read_csv(path)
    for (i in seq_len(nrow(m))) iad_check_hash(m$path[[i]], m$sha256[[i]])
  }

  lines <- c("<!-- Generated by _targets_training_status.R; do not edit. -->", "",
    "All twelve RT-only training-status families passed diagnostic, corrected-LOO and reporting gates.", "",
    "Physical-function findings have only one or two trained studies and are descriptive. Stacking weights are predictive weights, not probabilities of moderation. Contrast thresholds are descriptive for effect modification.", "")
  for (r in records) lines <- c(lines, paste0("### ", r$family, ": ", r$outcome), "",
    paste0("![](../manuscript/generated/training-status/", r$family, "/", r$outcome, "/protein_contrasts.png)"), "",
    paste0("[Paired differences in protein contrasts](../manuscript/generated/training-status/", r$family, "/", r$outcome,
      "/training_interaction_contrasts.csv), [support](../manuscript/generated/training-status/", r$family, "/", r$outcome, "/status_support.csv)."), "")

  path <- here::here("manuscript", "generated", "training-status", "_results.qmd")
  temporary <- paste0(path, ".tmp-", Sys.getpid()); on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary)
  if (!file.rename(temporary, path)) stop("Could not finalise training supplement.")
  path
}
