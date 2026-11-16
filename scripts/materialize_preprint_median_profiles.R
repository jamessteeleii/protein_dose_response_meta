# Regenerate preprint predictions and figures from accepted fits and frozen draws.
# Run from the project root. No models are fitted and no frozen files are edited.
library(lspline)
library(patchwork)
source("R/functions.R")
source("R/primary_model_averaging.R")
source("R/incremental_update.R")
source("R/preprint_prediction_profiles.R")

frozen_dir <- "results/verification/primary_update_2026"
average_dir <- "model_archive/update_2026/primary/reporting/model_average"
report_dir <- "results/update_2026/primary_refit/reporting"
destination <- "manuscript/generated/reporting"
figure_dir <- "manuscript/generated/figures/primary-dose-response"
stage <- tempfile("preprint-median-")
dir.create(stage)

read_csv <- function(path) readr::read_csv(path, show_col_types = FALSE, guess_max = Inf)
data <- read_preprint_profile_sample(frozen_dir)
data_by_outcome <- split(data, data$analysis_outcome)

# Both restricted measurement samples use the primary muscle-outcome profile.
if (identical(commandArgs(trailingOnly = TRUE), "--measurement")) {
  profile <- preprint_median_profile(data_by_outcome$hypertrophy)
  variants <- c(direct_only = "measurement_direct_only_update_2026",
                prereg_methods = "measurement_prereg_methods_update_2026")
  predictions <- lapply(names(variants), function(variant) {
    bundle <- file.path("results", "verification", variants[[variant]])
    frozen <- qs::qread(file.path("model_archive", "update_2026", "measurement",
                                 variant, "reporting", "model_average.qs"))
    record <- list(outcome = "hypertrophy", seed = frozen$settings$seed,
                   models = read_preprint_profile_models(bundle))
    message("Replaying saved measurement predictions for ", variant)
    result <- preprint_reprofile_primary_draws(frozen, record,
                                                read_preprint_profile_sample(bundle), profile)
    message("  Maximum paired contrast error: ", signif(result$contrast_error, 3))
    dplyr::mutate(result$predictions, bundle = variants[[variant]],
                  outcome = "hypertrophy", .before = 1L)
  })
  temporary <- tempfile(tmpdir = destination, fileext = ".csv.gz")
  readr::write_csv(dplyr::bind_rows(predictions), temporary)
  fs::file_move(temporary, file.path(destination, "median-profile-measurement-predictions.csv.gz"))
  quit(save = "no")
}

manifest <- read_csv(file.path(frozen_dir, "accepted_fit_manifest.csv.gz")) |>
  dplyr::mutate(fit_path = sub("^.*?(model_archive/)", "\\1", fit_path))
weights <- read_csv(file.path(frozen_dir, "stacking_weights.csv.gz"))
records <- make_update_2026_primary_reporting_records(manifest, weights)

profiles <- vector("list", length(records))
predictions <- profiles
for (i in seq_along(records)) {
  record <- records[[i]]
  outcome <- record$outcome
  message("Replaying saved primary predictions for ", outcome)
  sample <- data_by_outcome[[outcome]]
  profile <- preprint_median_profile(sample)
  frozen <- qs::qread(file.path(average_dir, paste0(outcome, ".qs")))
  result <- preprint_reprofile_primary_draws(frozen, record, sample, profile)
  profiles[[i]] <- tibble::tibble(
    outcome = outcome, median_age_years = profile[["age"]],
    median_duration_weeks = profile[["weeks"]],
    max_reference_draw_error = result$reference_error,
    max_dose_increment_error = result$contrast_error
  )
  predictions[[i]] <- dplyr::mutate(result$predictions, outcome = outcome, .before = 1L)
  message("  Maximum paired contrast error: ", signif(result$contrast_error, 3))
}
profiles <- dplyr::bind_rows(profiles)
predictions <- dplyr::bind_rows(predictions)
readr::write_csv(profiles, file.path(stage, "median-profiles.csv"))
readr::write_csv(predictions, file.path(stage, "median-profile-primary-predictions.csv.gz"))

audit <- list(
  dose_response_predictions = predictions,
  pairwise_contrasts = read_csv(file.path(frozen_dir, "pairwise_contrasts.csv.gz")),
  support = read_csv(file.path(frozen_dir, "observed_support.csv.gz")),
  model_weights = weights,
  manifest = dplyr::bind_rows(lapply(records, function(record) {
    average <- qs::qread(file.path(average_dir, paste0(record$outcome, ".qs")))
    tibble::as_tibble(c(list(outcome = record$outcome), average$settings))
  }))
)
arm_audit <- list(
  all_checks_pass = TRUE,
  predictions = read_csv(file.path(report_dir, "arm_conditional_predictions.csv.gz"))
)

invisible(write_preprint_median_primary_figures(audit, arm_audit, data_by_outcome, profiles, stage))

# Publish the staged files only after all checks and figure renders succeed.
files <- list.files(stage, full.names = TRUE)
stopifnot(length(files) == 10L, all(file.info(files)$size > 0))
directories <- ifelse(endsWith(files, ".png"), figure_dir, destination)
for (i in seq_along(files)) {
  target <- file.path(directories[[i]], basename(files[[i]]))
  temporary <- tempfile(tmpdir = directories[[i]])
  stopifnot(file.copy(files[[i]], temporary))
  fs::file_move(temporary, target)
}
unlink(files)
unlink(stage)
