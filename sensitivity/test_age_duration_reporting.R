# Reporting regression checks with synthetic predictions. No fits or pipeline.
stopifnot(getRversion() == package_version("4.4.3"))
suppressPackageStartupMessages(library(tidyverse))
for (file in c("R/functions.R", "R/primary_model_averaging.R", "R/verification_artifacts.R",
  "R/update_2026_absolute_exposures.R", "R/age_duration_interactions.R")) source(file)

# This reproduces the tibble data-mask collision: a profile column must not
# shadow the input record used for subsequent age and duration columns.
support_data <- tibble::tibble(study_no = c("1", "1", "2", "3", "4", "5"), arm_no = 1:6,
  resistance_exercise_code = c(0, 0, 0, 1, 1, 1), age_centre = c(0, 0, 30, 0, 5, 35),
  duration_centre = c(0, 0, 12, 0, 2, 12), protein_intake_raw = c(1, 2, 3, 1.5, 2.5, 4))
profile_row <- tibble::tibble(profile = "reference", age = 25, weeks = 12,
  age_neighbourhood_years = 10, duration_neighbourhood_weeks = 4)
support <- iad_profile_support(support_data, profile_row)
stopifnot(nrow(support) == 2L, all(support$profile == "reference"), all(support$age == 25),
  all(support$weeks == 12), identical(support$local_studies, c(1L, 2L)),
  identical(support$local_arms, c(2L, 2L)), identical(support$local_min, c(1, 1.5)),
  identical(support$local_max, c(2, 2.5)), all(support$age_within_range), all(support$duration_within_range))
outside <- profile_row; outside$age <- 100; outside$weeks <- 100
unsupported <- iad_profile_support(support_data, outside)
stopifnot(all(unsupported$local_studies == 0L), all(unsupported$local_arms == 0L),
  all(is.na(unsupported$local_min)), all(is.na(unsupported$local_max)),
  !any(unsupported$age_within_range), !any(unsupported$duration_within_range))

# Exercise the complete reporting formatter with a tiny grid and synthetic S3
# prediction methods. All output, including the static interface, stays in a
# temporary project; the real fitted models and manuscript are never touched.
test_reporting <- function() {
  env <- new.env(parent = .GlobalEnv)
  sys.source("R/age_duration_interactions.R", envir = env)
  root <- tempfile("iad-reporting-fixture-")
  dir.create(root)
  writeLines("# Synthetic reporting fixture", file.path(root, "fixture.R"))
  previous_directory <- getwd()
  on.exit({setwd(previous_directory); suppressMessages(here::i_am("README.md"))}, add = TRUE)
  setwd(root); suppressMessages(here::i_am("fixture.R"))
  env$iad_require_phase <- function(phase) stopifnot(phase == "report")
  env$write_verification_bundle <- write_verification_bundle
  environment(env$write_verification_bundle) <- env
  env$verification_git_commit <- function() NA_character_
  data <- tidyr::expand_grid(study_no = as.character(1:4), resistance_exercise_code = c(0, 1),
    protein_intake_raw = c(1, 1.5, 2)) |>
    dplyr::mutate(arm_no = dplyr::row_number(), outcome_no = arm_no, yi = 0.1, vi = 0.01,
      age_centre = rep(c(0, 10, 20, 30), each = 6L),
      duration_centre = rep(c(-4, 0, 4, 8), each = 6L), protein_intake = protein_intake_raw - 1.5)
  profiles <- tibble::tibble(family = "relative", outcome = "whole_body",
    profile = c("reference", "typical", "age_q25", "age_q75", "duration_q25", "duration_q75"),
    age = c(25, 40, 30, 50, 40, 40), weeks = c(12, 12, 12, 12, 8, 16),
    age_neighbourhood_years = 10, duration_neighbourhood_weeks = 4)
  env$iad_prediction_grid <- function(data, family, profile) {
    keys <- tidyr::expand_grid(resistance_exercise_code = c(0, 1), protein_intake = c(1, 1.5, 2)) |>
      dplyr::mutate(resistance_exercise = ifelse(resistance_exercise_code == 0, "YES", "NO"))
    newdata <- dplyr::mutate(keys, age_centre = profile$age - 25, duration_centre = profile$weeks - 12)
    list(keys = keys, newdata = newdata, grid_step = 0.5, protein_min = 1, reference = 1.5)
  }
  registerS3method("ndraws", "iad_reporting_fixture", function(x) 4000L, envir = asNamespace("posterior"))
  registerS3method("posterior_epred", "iad_reporting_fixture", function(object, newdata, re_formula, draw_ids, ...) {
    mean <- 0.2 * newdata$protein_intake + object$interaction * newdata$protein_intake *
      (0.01 * newdata$age_centre + 0.02 * newdata$duration_centre)
    outer(draw_ids / 100000, rep(1, nrow(newdata))) +
      matrix(rep(mean, each = length(draw_ids)), nrow = length(draw_ids))
  }, envir = asNamespace("rstantools"))
  paths <- vapply(0:1, function(interaction) {
    env$iad_qsave(list(outcome = "whole_body", model_index = interaction + 1L,
      model_name = paste0("fixture_", interaction),
      model = structure(list(data = data, interaction = interaction), class = "iad_reporting_fixture")),
      file.path(root, paste0("model_", interaction, ".qs")))
  }, character(1))
  models <- tibble::tibble(family = "relative", outcome = "whole_body", model_index = 1:2,
    parent_index = 1L, model_name = paste0("fixture_", 0:1), interaction = c("none", "both"),
    fit_path = paths, fit_sha256 = vapply(paths, iad_hash, character(1)), stacking_weight = c(.5, .5),
    primary_only_weight = c(1, 0))
  env$iad_specifications <- function(checkpoints) list(as.list(models[2L, ]))
  checkpoints <- list(relative = list(data = list(whole_body = data),
    source_manifest = tibble::tibble(source = "synthetic fixture")))
  record <- list(family = "relative", outcome = "whole_body", gate_passed = TRUE, models = models,
    seed = 1988L, diagnostics = tibble::tibble(diagnostic_pass = TRUE),
    loo_summary = tibble::tibble(model_index = 1:2))
  files <- env$iad_report_scope(record, checkpoints, profiles)
  stopifnot(all(file.exists(files)))
  directory <- env$iad_path("relative", "whole_body", "reporting")
  contrasts <- iad_read_csv(file.path(directory, "interaction_contrasts.csv.gz"))
  expected <- ifelse(contrasts$moderator == "age", .1, .08) * (contrasts$protein_intake - 1.5)
  stopifnot(nrow(contrasts) == 12L, max(abs(contrasts$estimate - expected)) < 1e-12)
  doses <- iad_read_csv(file.path(directory, "dose_response.csv.gz"))
  stopifnot(nrow(doses) == 72L, all(is.finite(doses$estimate)), dplyr::n_distinct(doses$profile) == 6L)
  manifest <- iad_read_csv(file.path(directory, "file_manifest.csv"))
  for (i in seq_len(nrow(manifest))) iad_check_hash(manifest$path[[i]], manifest$sha256[[i]])
}
test_reporting()
cat("Profile support, empty neighbourhoods, full synthetic reporting and paired moderator contrasts passed. No fitting or pipeline execution.\n")

if ("--review-current" %in% commandArgs(trailingOnly = TRUE)) {
  checkpoints <- targets::tar_read_raw("iad_checkpoints", store = "_targets_age_duration")
  profiles <- targets::tar_read_raw("iad_profile_plan", store = "_targets_age_duration")
  for (i in seq_len(nrow(profiles))) {
    row <- profiles[i, ]
    data <- checkpoints[[row$family]]$data[[row$outcome]]
    support <- iad_profile_support(data, row)
    grid <- iad_prediction_grid(data, row$family, row)
    stopifnot(nrow(support) == 2L, all(support$profile == row$profile),
      all(support$age == row$age), all(support$weeks == row$weeks),
      all(is.finite(grid$newdata$age_outcome_z)), all(is.finite(grid$newdata$duration_outcome_z)))
  }
  cat(nrow(profiles), "saved profile grids and support records passed read-only structural checks.\n")
}
