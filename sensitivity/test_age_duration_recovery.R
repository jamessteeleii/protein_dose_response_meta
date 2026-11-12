# Recovery and exact-assembly tests. No pipeline, compilation or sampling.
stopifnot(getRversion() == package_version("4.4.3"))
suppressPackageStartupMessages({ library(tidyverse); library(brms); library(mgcv); library(lspline) })
for (path in c("R/functions.R", "R/primary_model_diagnostics.R", "R/primary_model_loo.R",
  "R/primary_model_averaging.R", "R/incremental_update.R", "R/update_2026_absolute_exposures.R",
  "R/verification_artifacts.R", "R/age_duration_interactions.R", "R/age_duration_recovery.R",
  "R/age_duration_recovery_reporting.R")) source(path)
expect_error <- function(expression) stopifnot(inherits(tryCatch({ force(expression); NULL }, error = identity), "error"))
Sys.setenv(PROTEIN_INTERACTION_PHASE = "none")
for (phase in c("recover_pilot", "recover_fit", "recover_remedy", "recover_loo", "recover_exact", "recover_correction", "recover_accept", "recover_report")) {
  expect_error(iad_require_phase(phase))
}
expect_error(iad_recovery_fit_remedy(list(skip = TRUE), NULL, NULL, NULL))
expect_error(iad_recovery_reporting_records("missing.csv", NULL, character()))
expect_error(iad_recovery_report_scope(NULL, NULL))

# The muscle sampling-length remedy changes metadata, not which pointwise
# predictive scores stacking must use. Temporary normalisation must leave
# both input LOO objects untouched and agree with the pointwise optimiser.
test_stacking_lengths <- function() {
  set.seed(312)
  objects <- lapply(c(400L, 800L), function(n) suppressWarnings(
    loo::loo(matrix(rnorm(n * 5L, -1, 0.2), ncol = 5L))))
  before <- objects
  weights <- iad_stacking_weights(objects)
  expected <- as.numeric(loo::stacking_weights(do.call(cbind,
    lapply(objects, function(x) x$pointwise[, "elpd_loo"]))))
  stopifnot(identical(objects, before), length(weights) == 2L,
    abs(sum(weights) - 1) < 1e-8, isTRUE(all.equal(weights, expected)))
  bad <- objects; bad[[2L]]$pointwise <- bad[[2L]]$pointwise[-1L, , drop = FALSE]
  expect_error(iad_stacking_weights(bad))
}
test_stacking_lengths()

test_recovery_reporting_gate <- function() {
  previous_path <- iad_path
  root <- tempfile("iad-recovery-report-")
  assign("iad_path", function(..., archive = FALSE) file.path(root, ...), envir = .GlobalEnv)
  on.exit(assign("iad_path", previous_path, envir = .GlobalEnv), add = TRUE)
  accepted <- iad_csv(tibble::tibble(stage = "recovery_1"), file.path(root, "accepted.csv"))
  expect_error(iad_recovery_reporting_gate(accepted))
  paths <- unlist(lapply(iad_outcomes(), function(outcome) {
    directory <- file.path(root, outcome)
    files <- c(iad_csv(tibble::tibble(outcome = outcome,
      accepted_manifest_sha256 = iad_hash(accepted)), file.path(directory, "settings.csv")),
      iad_csv(tibble::tibble(check = "fixture", passed = TRUE), file.path(directory, "checks.csv")))
    manifest <- iad_csv(tibble::tibble(path = files, bytes = file.info(files)$size,
      sha256 = vapply(files, iad_hash, character(1))), file.path(directory, "file_manifest.csv"))
    c(files, manifest)
  }))
  iad_recovery_reporting_audit(paths, accepted)
  stopifnot(all(file.exists(iad_recovery_reporting_gate(accepted))))
  changed <- iad_csv(tibble::tibble(stage = "remedy_1"), file.path(root, "changed.csv"))
  expect_error(iad_recovery_reporting_gate(changed))
  iad_csv(tibble::tibble(check = "changed", passed = FALSE), paths[basename(paths) == "checks.csv"][[1L]])
  expect_error(iad_recovery_reporting_gate(accepted))
}
test_recovery_reporting_gate()

# A reviewed remedy replaces only its failed candidate and cannot overwrite
# the initial archive or bypass the diagnostic acceptance gate.
test_recovery_remedy <- function() {
  data <- tibble::tibble(study_no = "1", arm_no = 1, outcome_no = 1,
    time1 = 0, time2 = 12, yi = 0.1, vi = 0.1, protein_intake = 0,
    protein_intake_raw = 100, age_centre = 0, duration_centre = 0,
    resistance_exercise_code = 0)
  original <- list(family = "total", outcome = "strength", model_index = 24L,
    parent_index = 24L, formula = "yi ~ 1", stage = "recovery_1",
    interaction = "none", adapt_delta = 0.995, max_treedepth = 15L,
    chains = 4L, iter_warmup = 2000L, iter_sampling = 6000L)
  other <- original; other$model_index <- other$parent_index <- 2L
  specs <- list(original, other)
  initial <- tibble::tibble(family = "total", outcome = "strength", model_index = c(24L, 2L),
    stage = "recovery_1", diagnostic_pass = c(FALSE, TRUE),
    data_sha256 = iad_fingerprint(data), fit_path = c("failed.qs", "passed.qs"))
  remedies <- iad_remedy_specs(specs, "sensitivity/age-duration-recovery-remedies.csv", initial)
  stopifnot(length(remedies) == 1L, remedies[[1L]]$adapt_delta == 0.999,
    identical(remedies[[1L]][names(original)[!names(original) %in% c("adapt_delta", "stage")]],
      original[!names(original) %in% c("adapt_delta", "stage")]),
    iad_recovery_fit_path(original) != iad_recovery_fit_path(remedies[[1L]]),
    grepl("recovery/remedy_1/strength/model_24.qs$", iad_recovery_fit_path(remedies[[1L]])))
  bad <- original; bad$stage <- "invalid"
  expect_error(iad_recovery_fit_path(bad))
  checkpoint <- list(total = list(data = list(strength = data)))
  expect_error(iad_accepted(initial, tibble::tibble(), specs, checkpoint))
  replacement <- initial[1L, ]; replacement$stage <- "remedy_1"
  expect_error(iad_accepted(initial, replacement, specs, checkpoint))
  replacement$diagnostic_pass <- TRUE; replacement$fit_path <- "remedy.qs"
  accepted <- dplyr::bind_rows(iad_accepted(initial, replacement, specs, checkpoint))
  stopifnot(nrow(accepted) == 2L,
    accepted$fit_path[accepted$model_index == 24L] == "remedy.qs",
    accepted$fit_path[accepted$model_index == 2L] == "passed.qs")
}
test_recovery_remedy()

if (!"--assembly-only" %in% commandArgs(trailingOnly = TRUE)) {
checkpoints <- iad_read_checkpoints(iad_checkpoint_inputs())
specifications <- iad_recovery_specifications(checkpoints, "sensitivity/age-duration-total-recovery.csv")
stopifnot(length(specifications) == 72L, length(iad_recovery_pilot_specs(specifications)) == 8L,
  sum(vapply(specifications, function(s) s$iter_sampling == 12000L, logical(1))) == 1L)
model11 <- Filter(function(s) s$outcome == "hypertrophy" && s$model_index == 11L, specifications)[[1L]]
model20 <- Filter(function(s) s$outcome == "hypertrophy" && s$model_index == 20L, specifications)[[1L]]
stopifnot(model11$iter_sampling == 12000L, model20$adapt_delta == 0.9995)
expect_error(iad_apply_recovery(checkpoints, "sensitivity/age-duration-total-recovery.csv"))
stopifnot(identical(iad_apply_recovery(checkpoints, "sensitivity/age-duration-recovery-status.csv"), checkpoints))
priors <- iad_priors("sensitivity/age-duration-base-priors.csv")
preflight <- iad_preflight(specifications, checkpoints, priors)
stopifnot(nrow(preflight) == 72L, all(preflight$valid))
print(preflight |> dplyr::count(outcome, valid))
print(iad_storage_estimate(checkpoints, iad_specifications(checkpoints)))
historical <- iad_read_csv("results/verification/total_primary_update_2026/dose_response_predictions.csv.gz")
for (outcome in iad_outcomes()) {
  grid <- iad_recovery_reference_grid(checkpoints$total$data[[outcome]])
  old <- historical[historical$outcome == outcome, ]
  formulas <- Filter(function(s) s$outcome == outcome, specifications)
  variables <- unique(unlist(lapply(formulas, function(s) all.vars(stats::as.formula(s$formula)))))
  stopifnot(nrow(grid$keys) == nrow(old), identical(grid$keys$protein_intake, old$protein_intake),
    !length(setdiff(variables, names(grid$newdata))), all(grid$newdata$Intercept == 1),
    all(grid$newdata$age_centre == 0 & grid$newdata$duration_centre == 0),
    all(grid$newdata$body_mass_outcome_z == 0))
}
}

# Zero, multiple, incomplete, duplicated and failed exact corrections.
test_assembly <- function() {
  previous_path <- iad_path
  root <- tempfile("iad-exact-")
  assign("iad_path", function(..., archive = FALSE) file.path(root, ...), envir = .GlobalEnv)
  on.exit(assign("iad_path", previous_path, envir = .GlobalEnv), add = TRUE)
  manifest <- tibble::tibble(study_no = as.character(1:3), arm_no = 1:3,
    outcome_no = 1:3, yi = c(0, 0.1, 0.2), vi = 0.1)
  set.seed(1988)
  base <- suppressWarnings(loo::loo(matrix(rnorm(1200, -1, 0.1), ncol = 3)))
  base$diagnostics$pareto_k[] <- 0.1
  results <- lapply(27:28, function(i) list(family = "relative", outcome = "whole_body",
    model_index = as.integer(i), model_name = paste("Model", i), parent_index = 6L,
    interaction = "both", loo = base, observation_manifest = manifest,
    observation_diagnostics = tibble::tibble(observation = 1:3, integrated_pareto_k = 0.1)))
  make_audit <- function() {
    paths <- vapply(seq_along(results), function(i) iad_qsave(results[[i]],
      file.path(root, paste0("integrated_", i, ".qs"))), character(1))
    cases <- list()
    for (r in results) for (obs in which(loo::pareto_k_values(r$loo) > 0.7)) {
      cases[[length(cases) + 1L]] <- list(family = r$family, outcome = r$outcome,
        model_index = r$model_index, observation = obs)
    }
    list(manifest = tibble::tibble(family = "relative", outcome = "whole_body", model_index = 27:28,
      loo_path = paths, loo_sha256 = vapply(paths, iad_hash, character(1))),
      cases = if (length(cases)) cases else list(list(skip = TRUE)))
  }
  sentinel <- iad_qsave(list(skip = TRUE), file.path(root, "none.qs"))
  audit <- make_audit()
  final <- iad_finalise_loo(audit, sentinel)
  stopifnot(length(final) == 2L)
  results[[1]]$loo$diagnostics$pareto_k[1:2] <- 0.9
  results[[1]]$observation_diagnostics$integrated_pareto_k[1:2] <- 0.9
  results[[2]]$loo$diagnostics$pareto_k[[2]] <- 0.9
  results[[2]]$observation_diagnostics$integrated_pareto_k[[2]] <- 0.9
  audit <- make_audit()
  corrections <- lapply(audit$cases, function(case) {
    r <- results[[match(case$model_index, 27:28)]]
    r$loo$pointwise[case$observation, "elpd_loo"] <- -2 - case$observation
    r$exact_refit_summary <- tibble::tibble(observation = case$observation,
      diagnostic_pass = TRUE, exact_elpd_loo = -2 - case$observation)
    r
  })
  paths <- vapply(seq_along(corrections), function(i) iad_qsave(corrections[[i]],
    file.path(root, paste0("correction_", i, ".qs"))), character(1))
  expect_error(iad_finalise_loo(audit, paths[-1L]))
  expect_error(iad_finalise_loo(audit, c(paths, paths[[1L]])))
  final <- iad_finalise_loo(audit, paths)
  for (path in final) iad_validate_loo(qs::qread(path), manifest)
  stopifnot(nrow(qs::qread(final[[1L]])$exact_refit_summary) == 2L)
  corrections[[1L]]$exact_refit_summary$diagnostic_pass <- FALSE
  iad_qsave(corrections[[1L]], paths[[1L]])
  expect_error(iad_finalise_loo(audit, paths))
}
test_assembly()

# A compact branch can pass its saved correction forward without a refit object.
test_compact_exact_checkpoint <- function() {
  previous_path <- iad_path
  previous_guard <- iad_require_phase
  root <- tempfile("iad-exact-checkpoint-")
  assign("iad_path", function(..., archive = FALSE) file.path(root, ...), envir = .GlobalEnv)
  assign("iad_require_phase", function(phase) stopifnot(phase == "correction"), envir = .GlobalEnv)
  on.exit({
    assign("iad_path", previous_path, envir = .GlobalEnv)
    assign("iad_require_phase", previous_guard, envir = .GlobalEnv)
  })
  fit_path <- iad_qsave(list(fixture = "fit"), file.path(root, "fit.qs"))
  loo_path <- iad_qsave(list(fixture = "loo"), file.path(root, "loo.qs"))
  record <- list(family = "relative", outcome = "strength", model_index = 27L, observation = 2L,
    fit_path = fit_path, fit_sha256 = iad_hash(fit_path), loo_path = loo_path, loo_sha256 = iad_hash(loo_path))
  diagnostic <- tibble::tibble(observation = 2L, diagnostic_pass = TRUE)
  correction <- list(family = "relative", outcome = "strength", model_index = 27L,
    exact_refit_summary = dplyr::mutate(diagnostic, exact_elpd_loo = -2))
  checkpoint <- list(compact = TRUE, interaction_record = record,
    diagnostics = diagnostic, correction = correction)
  path <- iad_qsave(checkpoint, file.path(root, "passed.qs"))
  stopifnot(identical(iad_exact_diagnostics(path), diagnostic),
    identical(qs::qread(iad_score_exact(path)), correction))
  checkpoint$correction$exact_refit_summary$observation <- 3L
  wrong <- iad_qsave(checkpoint, file.path(root, "wrong_case.qs"))
  expect_error(iad_score_exact(wrong))
  checkpoint$diagnostics$diagnostic_pass <- FALSE
  checkpoint$correction <- NULL
  failed <- iad_qsave(checkpoint, file.path(root, "failed.qs"))
  stopifnot(identical(iad_exact_diagnostics(failed)$diagnostic_pass, FALSE))
  expect_error(iad_score_exact(failed))
}
test_compact_exact_checkpoint()

# Leaf-pruning guards must reject every primary-family key.
source("sensitivity/age-duration-storage.R")
for (key in c("primary", "relative_primary", "total_primary", "supplement_primary", "measurement_direct_only")) {
  expect_error(iad_storage_plan(key))
}
cat("Recovery checks completed; exact-LOO assembly and protected-family checks passed. No sampling or pruning.\n")
