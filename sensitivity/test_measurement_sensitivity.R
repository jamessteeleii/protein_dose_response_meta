# Sampling-free checks. Run with the project R/renv environment.
source("R/functions.R")
source("R/primary_model_diagnostics.R")
source("R/primary_model_loo.R")
source("R/primary_model_averaging.R")
source("R/incremental_update.R")
source("R/sensitivity_pilots.R")
source("R/measurement_sensitivity.R")
source("R/measurement_targets.R")

expect_error <- function(expression) {
  stopifnot(inherits(tryCatch({ force(expression); NULL }, error = identity), "error"))
}
Sys.setenv(PROTEIN_MEASUREMENT_PHASES = "none")
stopifnot(!any(vapply(c("pilot", "fit", "remedy", "loo", "exact", "report"),
  measurement_phase_enabled, logical(1))))

checkpoint <- read_measurement_checkpoint(measurement_checkpoint_inputs())
annotations <- read_measurement_annotations("sensitivity/measurement-annotations.csv")
sets <- lapply(measurement_variants(), function(variant) {
  select_measurement_data(checkpoint$data, variant, annotations = annotations)
})
names(sets) <- measurement_variants()
stopifnot(nrow(sets$prereg_methods$data) == 413L, nrow(sets$direct_only$data) == 255L)
confirmed <- sets$prereg_methods$data[sets$prereg_methods$data$measurement_annotation_applied, ]
stopifnot(nrow(confirmed) == 2L, setequal(confirmed$outcome_no, c(66, 72)),
  all(confirmed$measurement_method == "dexa"),
  all(confirmed$measurement_method_reviewed == "dexa - appendicular"),
  !any(sets$direct_only$data$measurement_annotation_applied))
for (selection in sets) {
  flagged <- selection$selection[selection$selection$measurement_annotation_applied, ]
  stopifnot(nrow(flagged) == 2L,
    all(flagged$retained == (selection$variant == "prereg_methods")))
}
bad_annotations <- annotations
bad_annotations$arm_no[[1]] <- "99999"
expect_error(select_measurement_data(checkpoint$data, "prereg_methods", bad_annotations))
bad_annotations <- annotations
bad_annotations$original_measurement_method[[1]] <- "mri"
expect_error(select_measurement_data(checkpoint$data, "prereg_methods", bad_annotations))
expect_error(select_measurement_data(checkpoint$data, "prereg_methods", rbind(annotations, annotations)))
original <- checkpoint$data[checkpoint$data$outcome_group == "direct muscle size/mass", ]
stopifnot(sum(grepl("trunk", original$measurement_method)) == 6L,
  sum(grepl("^skinfold", original$measurement_method)) == 6L)
for (selection in sets) {
  rows <- selection$selection$source_row[selection$selection$retained]
  for (column in c("yi", "vi", "protein_intake", "protein_intake_raw", "protein_median",
                   "study_no", "arm_no", "outcome_no", "time1", "time2")) {
    stopifnot(identical(selection$data[[column]], checkpoint$data[[column]][rows]))
  }
  stopifnot(!any(grepl("trunk|^skinfold", selection$data$measurement_method)))
  copy <- as.data.frame(selection$data)
  copy$study_no <- factor(copy$study_no)
  stopifnot(identical(measurement_data_fingerprint(selection$data), measurement_data_fingerprint(copy)))
  copy$yi[[1L]] <- copy$yi[[1L]] + 0.01
  stopifnot(!identical(measurement_data_fingerprint(selection$data), measurement_data_fingerprint(copy)))
}
# Archived brms objects retain only model-frame columns. Metadata can be absent,
# but row order, response, grouping identities, and all predictors must match.
fit_fixture <- list(specification = list(standardised_moderators = FALSE,
  half_unit_protein = FALSE), model = list(data = as.data.frame(
    sets$direct_only$data[c("yi", "vi", "arm_no", "outcome_no",
      "study_no", "resistance_exercise_code", "protein_intake")])))
fit_fixture$model$data$Intercept <- 1
fit_fixture$model$data[["lspline(protein_intake, 0.5)"]] <-
  lspline::lspline(fit_fixture$model$data$protein_intake, 0.5)
stopifnot(isTRUE(validate_measurement_fit_data(fit_fixture, sets$direct_only)))
bad_fit <- fit_fixture
bad_fit$model$data <- bad_fit$model$data[nrow(bad_fit$model$data):1L, ]
expect_error(validate_measurement_fit_data(bad_fit, sets$direct_only))
bad_fit <- fit_fixture
bad_fit$model$data$protein_intake[[1L]] <- bad_fit$model$data$protein_intake[[1L]] + 0.01
expect_error(validate_measurement_fit_data(bad_fit, sets$direct_only))
bad_fit <- fit_fixture
bad_fit$model$data$vi <- NULL
expect_error(validate_measurement_fit_data(bad_fit, sets$direct_only))
bad_fit <- fit_fixture
bad_fit$model$data[["lspline(protein_intake, 0.5)"]][1L, 1L] <- 99
expect_error(validate_measurement_fit_data(bad_fit, sets$direct_only))

# These fixtures protect the boundary: unfamiliar methods fail rather than
# joining a broad 'not excluded' category; membership is tied to outcome group.
fixture <- checkpoint$data[checkpoint$data$outcome_group == "direct muscle size/mass", ][1:2, ]
fixture$resistance_exercise_code <- c(0, 1)
fixture$measurement_method <- c("mri - quadriceps", "dexa - appendicular")
stopifnot(nrow(select_measurement_data(fixture, "prereg_methods")$data) == 2L)
fixture$measurement_method[[2L]] <- "new unreviewed proxy"
expect_error(select_measurement_data(fixture, "prereg_methods"))
fixture$measurement_method[[2L]] <- "dexa - trunk"
expect_error(select_measurement_data(fixture, "prereg_methods")) # only one RT stratum remains

for (variant in measurement_variants()) {
  specs <- make_measurement_specifications(checkpoint, variant)
  stopifnot(length(specs) == 26L,
    identical(vapply(specs, function(x) x$formula, character(1)),
      stats::setNames(checkpoint$specifications$formula, names(specs))),
    all(vapply(specs, function(x) identical(x$measurement_variant, variant), logical(1))))
  stopifnot(identical(unname(vapply(measurement_pilot_specifications(specs),
    function(x) x$model_index, integer(1))), c(2L, 26L)),
    length(specs[-c(2L, 26L)]) == 24L)
  preflight <- data.frame(valid = rep(TRUE, 26L))
  # A missing opt-in blocks fitting before any model can be opened or sampled.
  expect_error(fit_measurement_candidate(specs[[1L]], sets[[variant]], NULL, preflight))
  expect_error(fit_measurement_exact(list(skip = TRUE, variant = variant)))
  Sys.setenv(PROTEIN_MEASUREMENT_PHASES = "fit")
  expect_error(fit_measurement_candidate(specs[[1L]], sets[[variant]], NULL, preflight))
  Sys.setenv(PROTEIN_MEASUREMENT_PHASES = "none")
}
new_targets <- make_measurement_targets()
target_names <- vapply(new_targets, function(x) x$name, character(1))
stopifnot(length(new_targets) == 64L, !anyDuplicated(target_names))

# Test arbitrary exact-refit counts (including zero and multiple cases for a
# single candidate). All inputs are synthetic; outputs go to an R temp folder.
test_exact_assembly <- function() {
  original_path <- measurement_path
  test_root <- tempfile("measurement-exact-test-")
  assign("measurement_path", function(variant, stage, filename, archive = TRUE) {
    file.path(test_root, variant, stage, filename)
  }, envir = .GlobalEnv)
  on.exit(assign("measurement_path", original_path, envir = .GlobalEnv), add = TRUE)
  set.seed(1988)
  base_loo <- suppressWarnings(loo::loo(matrix(rnorm(1200, -1, 0.1), ncol = 3)))
  observation_manifest <- tibble::tibble(outcome = "hypertrophy", observation = 1:3,
    study_no = as.character(1:3), yi = c(0, 0.1, 0.2), vi = rep(0.01, 3))
  results <- lapply(seq_len(26), function(index) list(
    outcome = "hypertrophy", model_index = as.integer(index), model_name = paste("Model", index),
    measurement_variant = "direct_only", loo = base_loo, observation_manifest = observation_manifest,
    observation_diagnostics = tibble::tibble(observation = 1:3, integrated_pareto_k = rep(0.1, 3))))
  for (i in seq_along(results)) results[[i]]$loo$diagnostics$pareto_k <- rep(0.1, 3)
  write_results <- function() {
    paths <- vapply(seq_along(results), function(i) write_measurement_qs(results[[i]],
      measurement_path("direct_only", "loo/integrated", paste0(i, ".qs"))), character(1))
    tibble::tibble(model_index = seq_len(26), integrated_loo_path = paths)
  }
  manifest <- write_results()
  sentinel <- write_measurement_qs(list(skip = TRUE),
    measurement_path("direct_only", "loo/exact", "none.qs"))
  final_paths <- finalise_measurement_loo(manifest, sentinel, "direct_only")
  stopifnot(length(final_paths) == 26L,
    nrow(qs::qread(final_paths[[1]])$exact_refit_summary) == 0L)
  results[[1]]$loo$diagnostics$pareto_k[1:2] <- 0.9
  results[[1]]$observation_diagnostics$integrated_pareto_k[1:2] <- 0.9
  manifest <- write_results()
  correction_paths <- vapply(1:2, function(observation) {
    correction <- results[[1]]
    correction$loo$pointwise[observation, "elpd_loo"] <- -2 - observation
    correction$exact_refit_summary <- tibble::tibble(model_index = 1L,
      observation = observation, diagnostic_pass = TRUE, exact_elpd_loo = -2 - observation)
    write_measurement_qs(correction, measurement_path("direct_only", "loo/exact",
      paste0(observation, ".qs")))
  }, character(1))
  expect_error(finalise_measurement_loo(manifest, correction_paths[1], "direct_only"))
  final_paths <- finalise_measurement_loo(manifest, correction_paths, "direct_only")
  first <- qs::qread(final_paths[[1]])
  stopifnot(identical(which(is.na(first$loo$diagnostics$pareto_k)), 1:2),
    identical(as.numeric(first$loo$pointwise[1:2, "elpd_loo"]), c(-3, -4)),
    nrow(first$exact_refit_summary) == 2L,
    identical(first$observation_diagnostics$exact_loo_refit, c(TRUE, TRUE, FALSE)))
  weights <- measurement_stacking(final_paths, "direct_only")
  stopifnot(nrow(weights) == 26, abs(sum(weights$stacking_weight) - 1) < 1e-8)
}
test_exact_assembly()

if ("--preflight" %in% commandArgs(trailingOnly = TRUE)) {
  suppressPackageStartupMessages({ library(brms); library(mgcv); library(lspline); library(splines) })
  # Read a small already-frozen prior object only; never build upstream targets.
  priors <- targets::tar_read(steele_priors)
  for (variant in measurement_variants()) {
    cat("Building Stan input matrices (no sampling):", variant, "\n")
    flush.console()
    checks <- measurement_input_preflight(sets[[variant]],
      make_measurement_specifications(checkpoint, variant), priors)
    stopifnot(nrow(checks) == 26L, all(checks$valid))
    cat("Passed 26 input checks:", variant, "\n")
    flush.console()
  }
}
cat("Measurement selection, immutable values, registries, execution guards, and exact-LOO assembly passed.\n")
