# Review the recovered primary family before adding interaction candidates.
# New outputs are isolated from the immutable historical primary report.
iad_recovery_reporting_inputs <- function() {
  file.path(here::here("results", "verification", "total_primary_update_2026"),
    c("manifest.csv", "stacking_weights.csv.gz", "dose_response_predictions.csv.gz"))
}

# Reporting records for recovered parents ----

# Join the recovery acceptance mapping to the frozen parent data and weights
# before checking that each recovered ensemble can support reporting.
iad_recovery_reporting_records <- function(accepted_path, checkpoints, source_files) {
  iad_require_phase("recover_report")
  if (!all(file.exists(source_files))) stop("Missing historical reporting inputs.")
  iad_apply_recovery(checkpoints, accepted_path)
  accepted <- iad_read_csv(accepted_path)
  directory <- dirname(source_files[[1L]])
  manifest <- iad_read_csv(file.path(directory, "manifest.csv"))
  historical_weights <- iad_bundle_table(directory, manifest, "stacking_weights")
  historical_dose <- iad_bundle_table(directory, manifest, "dose_response_predictions")
  lapply(iad_outcomes(), function(outcome) {
    models <- accepted[accepted$outcome == outcome, ]
    models <- models[order(models$model_index), ]
    results <- lapply(seq_len(nrow(models)), function(i) {
      iad_check_hash(models$loo_path[[i]], models$loo_sha256[[i]])
      r <- qs::qread(models$loo_path[[i]])
      if (r$outcome != outcome || r$model_index != models$model_index[[i]] ||
          r$fit_path != models$fit_path[[i]] || r$fit_sha256 != models$sha256[[i]]) {
        stop("Recovery reporting fit/LOO identity mismatch.")
      }
      iad_validate_loo(r, checkpoints$total$data[[outcome]])
      r
    })
    if (nrow(models) != 18L || !all(models$diagnostic_pass)) stop("Incomplete recovered ensemble.")
    models$stacking_weight <- iad_stacking_weights(lapply(results, `[[`, "loo"))
    previous <- historical_weights[historical_weights$outcome == outcome, ]
    previous <- previous[match(models$model_index, previous$model_index), ]
    if (anyNA(previous$model_index) || any(previous$model_name != models$model_name)) {
      stop("Historical weight identities differ from the recovery.")
    }
    models$historical_weight <- previous$weight
    list(outcome = outcome, models = models, seed = 941988L + match(outcome, iad_outcomes()),
      historical_dose = historical_dose[historical_dose$outcome == outcome, ],
      accepted_manifest_sha256 = iad_hash(accepted_path))
  })
}

iad_recovery_reference_grid <- function(data) {
  grid <- make_update_2026_total_primary_prediction_grid(data)
  grid$newdata$Intercept <- 1
  grid
}

iad_recovery_model_average <- function(record, grid, draws) {
  structure(list(outcome = record$outcome, draws = draws, grid = grid$keys,
    model_weights = stats::setNames(record$models$stacking_weight, record$models$model_name),
    settings = list(grid_step = grid$grid_step, protein_min = grid$protein_min,
      sesoi = 0.1, conf_level = 0.95, n_draws = nrow(draws), age = 25, weeks = 12)),
    class = "primary_model_average")
}

# Compare the recovered total-intake ensemble with the historical curve on
# the same reference grid and write an audit for this outcome.
iad_recovery_report_scope <- function(record, checkpoints) {
  iad_require_phase("recover_report")
  outcome <- record$outcome; models <- record$models
  data <- checkpoints$total$data[[outcome]]
  grid <- iad_recovery_reference_grid(data)
  n_draws <- 4000L
  allocation <- iad_model_allocation(models, n_draws, record$seed)
  draws <- matrix(NA_real_, n_draws, nrow(grid$newdata))

  for (i in seq_len(nrow(models))) {
    selected <- which(allocation == i)
    if (!length(selected)) next
    iad_check_hash(models$fit_path[[i]], models$sha256[[i]])
    candidate <- qs::qread(models$fit_path[[i]])
    if (candidate$outcome != outcome || candidate$model_index != models$model_index[[i]] ||
        candidate$model_name != models$model_name[[i]] ||
        candidate$data_sha256 != models$data_sha256[[i]]) stop("Recovery prediction identity mismatch.")
    iad_validate_parent_data(candidate, data, "total")
    set.seed(record$seed + models$model_index[[i]])
    ids <- sample.int(posterior::ndraws(candidate$model), length(selected))
    prediction <- as.matrix(brms::posterior_epred(candidate$model,
      newdata = grid$newdata, re_formula = NA, draw_ids = ids))
    if (!identical(dim(prediction), c(length(selected), nrow(grid$newdata))) ||
        any(!is.finite(prediction))) stop("Invalid recovery predictions.")
    draws[selected, ] <- prediction
    rm(candidate, prediction); gc(verbose = FALSE)
  }

  if (anyNA(draws)) stop("Incomplete recovered posterior mixture.")
  average <- iad_recovery_model_average(record, grid, draws)
  dose <- summarise_primary_dose_response(average)
  increments <- draws
  for (rt in c(0, 1)) {
    columns <- which(grid$keys$resistance_exercise_code == rt)
    reference <- which(grid$keys$resistance_exercise_code == rt &
      abs(grid$keys$protein_intake - grid$reference_g_day) < 1e-10)
    if (length(reference) != 1L) stop("Missing unique recovery reference dose.")
    increments[, columns] <- sweep(draws[, columns, drop = FALSE], 1L, draws[, reference], "-")
  }

  contrasts <- compute_primary_heatmap_contrasts(average, data, protein_min = grid$protein_min)
  previous <- record$historical_dose
  keys <- c("resistance_exercise_code", "protein_intake")
  previous <- dplyr::select(previous, dplyr::all_of(c(keys, "estimate", "conf.low", "conf.high")))
  comparison <- dplyr::left_join(dose, previous, by = keys, suffix = c("_recovered", "_historical"))
  if (nrow(comparison) != nrow(dose) || nrow(previous) != nrow(dose) ||
      anyNA(comparison$estimate_historical)) stop("Historical and recovered reporting grids differ.")

  comparison$estimate_difference <- comparison$estimate_recovered - comparison$estimate_historical
  checks <- tibble::tibble(check = c("18 accepted primary candidates", "4000 finite mixture draws",
    "weights sum to one", "age 25 and duration 12", "body mass at outcome mean",
    "historical grid matches", "reference-dose increments exactly zero"),
    passed = c(nrow(models) == 18L && all(models$diagnostic_pass), all(is.finite(draws)),
      abs(sum(models$stacking_weight) - 1) < 1e-8,
      all(grid$newdata$age_centre == 0 & grid$newdata$duration_centre == 0),
      all(grid$newdata$body_mass_outcome_z == 0), nrow(comparison) == nrow(previous),
      all(increments[, abs(grid$keys$protein_intake - grid$reference_g_day) < 1e-10] == 0)))

  if (!all(checks$passed)) stop("Recovery reporting checks failed.")
  allocation_table <- dplyr::mutate(models, draws = tabulate(allocation, nrow(models)),
    weight_difference = stacking_weight - historical_weight)
  settings <- tibble::tibble(outcome = outcome, age = 25, weeks = 12,
    body_mass_outcome_z = 0, reference_g_day = grid$reference_g_day,
    grid_step_g_day = grid$grid_step, n_draws = n_draws, seed = record$seed,
    accepted_manifest_sha256 = record$accepted_manifest_sha256,
    comparison_note = "Fresh draws and weights; numerical agreement requires review, not an automatic tolerance.")
  tables <- list(checks = checks, stacking_weights = allocation_table, settings = settings,
    dose_response = dose,
    reference_contrasts = dplyr::bind_cols(grid$keys, iad_summarise_draws(increments)),
    historical_curve_comparison = comparison, pairwise_contrasts = contrasts$contrasts,
    pairwise_support = contrasts$support,
    local_equivalence = summarise_primary_local_equivalence(contrasts$contrasts, outcome))
  root <- iad_path("recovery", "reporting", outcome)
  paths <- purrr::imap_chr(tables, function(table, name) iad_csv(table,
    file.path(root, paste0(name, if (name %in% c("pairwise_contrasts", "historical_curve_comparison")) ".csv.gz" else ".csv"))))
  file_manifest <- tibble::tibble(path = unname(paths), bytes = file.info(paths)$size,
    sha256 = vapply(paths, iad_hash, character(1)))
  c(unname(paths), iad_csv(file_manifest, file.path(root, "file_manifest.csv")))
}

# Audit and reporting gate ----

# Collect the recovery report results so the interaction workflow can check
# the accepted recovery mapping before proceeding.
iad_recovery_reporting_audit <- function(paths, accepted_path) {
  paths <- unname(unlist(paths))
  manifests <- paths[basename(paths) == "file_manifest.csv"]
  if (length(manifests) != 4L || any(!file.exists(paths))) stop("Incomplete recovery reporting.")
  rows <- lapply(manifests, function(path) {
    m <- iad_read_csv(path)
    for (i in seq_len(nrow(m))) iad_check_hash(m$path[[i]], m$sha256[[i]])
    settings <- iad_read_csv(m$path[basename(m$path) == "settings.csv"])
    checks <- iad_read_csv(m$path[basename(m$path) == "checks.csv"])
    if (!all(checks$passed) || settings$accepted_manifest_sha256 != iad_hash(accepted_path)) {
      stop("Recovery reporting is invalid or belongs to another accepted family.")
    }
    tibble::tibble(outcome = settings$outcome, reporting_checks_pass = TRUE,
      file_manifest = path, sha256 = iad_hash(path),
      accepted_manifest_sha256 = settings$accepted_manifest_sha256)
  })

  rows <- dplyr::bind_rows(rows)
  if (!setequal(rows$outcome, iad_outcomes()) || anyDuplicated(rows$outcome)) stop("Missing recovery outcome.")
  iad_csv(rows, iad_path("recovery", "reporting", "audit.csv"))
}

iad_recovery_reporting_gate <- function(mapping_path) {
  if (!nrow(iad_read_csv(mapping_path))) return(mapping_path)
  path <- iad_path("recovery", "reporting", "audit.csv")
  if (!file.exists(path)) stop("Complete and review recover_report before starting the interaction expansion.")
  audit <- iad_read_csv(path)
  if (nrow(audit) != 4L || !setequal(audit$outcome, iad_outcomes()) ||
      !all(audit$reporting_checks_pass) || any(audit$accepted_manifest_sha256 != iad_hash(mapping_path))) {
    stop("Recovery reporting checkpoint is incomplete or stale.")
  }
  for (i in seq_len(nrow(audit))) iad_check_hash(audit$file_manifest[[i]], audit$sha256[[i]])
  children <- lapply(audit$file_manifest, function(manifest_path) {
    manifest <- iad_read_csv(manifest_path)
    for (i in seq_len(nrow(manifest))) iad_check_hash(manifest$path[[i]], manifest$sha256[[i]])
    manifest$path
  })
  c(path, audit$file_manifest, unname(unlist(children)))
}
