#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
  library(qs)
  library(digest)
})

source("R/functions.R")
source("R/primary_model_diagnostics.R")
source("R/primary_model_averaging.R")
source("R/sensitivity_pilots.R")

required_targets <- c(
  "sensitivity_primary_data_by_outcome",
  "sensitivity_priors_by_outcome",
  "sensitivity_final_model_specifications",
  "sensitivity_within_remedy_ensemble_fits",
  "sensitivity_within_remedy_diagnostics"
)
missing_targets <- required_targets[!vapply(
  required_targets,
  targets::tar_exist_objects,
  logical(1)
)]
if (length(missing_targets)) {
  stop(
    "Cannot freeze the within-pilot checkpoint; cached targets are missing: ",
    paste(missing_targets, collapse = ", "),
    call. = FALSE
  )
}

analysis_data <- targets::tar_read_raw(
  "sensitivity_primary_data_by_outcome"
)
priors <- targets::tar_read_raw("sensitivity_priors_by_outcome")
specifications <- targets::tar_read_raw(
  "sensitivity_final_model_specifications"
)
fits <- targets::tar_read_raw("sensitivity_within_remedy_ensemble_fits")
diagnostics <- targets::tar_read_raw(
  "sensitivity_within_remedy_diagnostics"
)

diagnostic_summary <- dplyr::bind_rows(purrr::map(
  diagnostics,
  "summary"
))
parameter_failures <- purrr::map_dfr(diagnostics, function(result) {
  result$parameter_failures |>
    dplyr::mutate(
      outcome = result$summary$outcome[[1L]],
      model_index = result$summary$model_index[[1L]],
      .before = 1L
    )
})
expected_indices <- c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
if (
  length(fits) != 7L || nrow(diagnostic_summary) != 7L ||
    !identical(
      sort(unname(purrr::map_int(fits, "model_index"))),
      expected_indices
    ) ||
    anyNA(diagnostic_summary$diagnostic_pass) ||
    !all(diagnostic_summary$diagnostic_pass)
) {
  stop("The accepted within-pilot ensemble is invalid.", call. = FALSE)
}

fit_paths <- purrr::map_chr(fits, function(candidate) {
  write_sensitivity_fit_file(candidate, stage = "pilot_accepted")
})
fit_manifest <- make_sensitivity_fit_file_manifest(
  fit_paths,
  stage = "pilot_accepted"
) |>
  dplyr::arrange(model_index)

checkpoint <- list(
  checkpoint_version = 1L,
  family = "within",
  generated_utc = format(
    Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"
  ),
  analysis_data_by_outcome = analysis_data,
  priors_by_outcome = priors,
  model_specifications = specifications,
  reused_fit_manifest = fit_manifest,
  pilot_diagnostic_summary = diagnostic_summary,
  pilot_parameter_failures = parameter_failures,
  audit = list(
    accepted_pilot_fits = nrow(fit_manifest),
    all_pilot_diagnostics_pass = all(
      diagnostic_summary$diagnostic_pass
    )
  )
)

output_directory <- file.path("results", "checkpoints")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(
  output_directory, "within_pilot_completed.qs"
)
temporary_path <- paste0(output_path, ".tmp")
on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
qs::qsave(checkpoint, temporary_path, preset = "high")
if (!file.copy(temporary_path, output_path, overwrite = TRUE)) {
  stop("Could not install the within-pilot checkpoint.", call. = FALSE)
}
unlink(temporary_path, force = TRUE)

cat(
  normalizePath(output_path, winslash = "/", mustWork = TRUE),
  "\n",
  sep = ""
)
