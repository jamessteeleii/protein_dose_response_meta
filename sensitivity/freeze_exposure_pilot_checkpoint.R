#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
  library(qs)
  library(digest)
  library(dplyr)
  library(purrr)
})

required_targets <- c(
  "sensitivity_primary_data_by_outcome",
  "sensitivity_priors_by_outcome",
  "sensitivity_final_model_specifications",
  "sensitivity_exposure_comparison_data",
  "sensitivity_achieved_final_ensemble_fits",
  "sensitivity_assignment_remedy_ensemble_fits",
  "sensitivity_achieved_final_diagnostics",
  "sensitivity_assignment_remedy_diagnostics"
)

missing_targets <- required_targets[!vapply(
  required_targets,
  targets::tar_exist_objects,
  logical(1)
)]
if (length(missing_targets)) {
  stop(
    "Cannot freeze the exposure pilot; cached targets are missing: ",
    paste(missing_targets, collapse = ", "),
    call. = FALSE
  )
}

write_candidate <- function(candidate, stage = "pilot_final") {
  if (
    !is.list(candidate) || !inherits(candidate$model, "brmsfit") ||
      !candidate$family %in% c("achieved", "assignment") ||
      !identical(candidate$outcome, "whole_body")
  ) {
    stop("A completed exposure-pilot candidate is invalid.",
         call. = FALSE)
  }
  directory <- file.path(
    "_targets", "user", "sensitivity_fit_files",
    candidate$family,
    candidate$outcome,
    paste0("model_", sprintf("%02d", candidate$model_index))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, paste0(stage, ".qs"))
  temporary_path <- paste0(path, ".tmp")
  on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_path, preset = "high")
  if (!file.copy(temporary_path, path, overwrite = TRUE)) {
    stop("Could not install an exposure-pilot fit file.", call. = FALSE)
  }
  unlink(temporary_path, force = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

make_manifest <- function(fits, family) {
  expected_indices <- c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
  if (
    length(fits) != 7L ||
      !identical(unname(sort(map_int(fits, "model_index"))), expected_indices) ||
      any(map_chr(fits, "family") != family)
  ) {
    stop("The completed exposure-pilot ensemble is invalid.",
         call. = FALSE)
  }
  paths <- map_chr(fits, write_candidate)
  tibble::tibble(
    family = map_chr(fits, "family"),
    outcome = map_chr(fits, "outcome"),
    model_index = map_int(fits, "model_index"),
    model_name = map_chr(fits, "model_name"),
    stage = "pilot_final",
    fit_path = paths,
    bytes = file.info(paths)$size,
    sha256 = map_chr(paths, ~ digest::digest(
      file = .x, algo = "sha256"
    ))
  ) |>
    arrange(model_index)
}

summarise_diagnostics <- function(diagnostics, family) {
  summary <- bind_rows(map(diagnostics, "summary")) |>
    arrange(model_index)
  parameter_failures <- map_dfr(diagnostics, function(result) {
    result$parameter_failures |>
      mutate(
        outcome = result$summary$outcome[[1L]],
        model_index = result$summary$model_index[[1L]],
        .before = 1L
      )
  })
  if (
    nrow(summary) != 7L || anyDuplicated(summary$model_index) ||
      !all(summary$diagnostic_pass)
  ) {
    stop("The final ", family, " pilot diagnostics do not pass.",
         call. = FALSE)
  }
  list(summary = summary, parameter_failures = parameter_failures)
}

analysis_data <- targets::tar_read_raw(
  "sensitivity_primary_data_by_outcome"
)
priors <- targets::tar_read_raw(
  "sensitivity_priors_by_outcome"
)
specifications <- targets::tar_read_raw(
  "sensitivity_final_model_specifications"
)
pilot_comparison <- targets::tar_read_raw(
  "sensitivity_exposure_comparison_data"
)

achieved_fits <- targets::tar_read_raw(
  "sensitivity_achieved_final_ensemble_fits"
)
achieved_manifest <- make_manifest(achieved_fits, "achieved")
rm(achieved_fits)
invisible(gc())

assignment_fits <- targets::tar_read_raw(
  "sensitivity_assignment_remedy_ensemble_fits"
)
assignment_manifest <- make_manifest(assignment_fits, "assignment")
rm(assignment_fits)
invisible(gc())

achieved_diagnostics <- summarise_diagnostics(
  targets::tar_read_raw("sensitivity_achieved_final_diagnostics"),
  "achieved"
)
assignment_diagnostics <- summarise_diagnostics(
  targets::tar_read_raw("sensitivity_assignment_remedy_diagnostics"),
  "assignment"
)

checkpoint <- list(
  checkpoint_version = 1L,
  family = "exposure",
  generated_utc = format(
    Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"
  ),
  analysis_data_by_outcome = analysis_data,
  priors_by_outcome = priors,
  model_specifications = specifications,
  achieved_reused_fit_manifest = achieved_manifest,
  assignment_reused_fit_manifest = assignment_manifest,
  pilot_common_sample_audit = pilot_comparison$audit,
  pilot_diagnostic_summary = bind_rows(
    mutate(
      achieved_diagnostics$summary,
      exposure = "achieved", .before = 1L
    ),
    mutate(
      assignment_diagnostics$summary,
      exposure = "assignment_derived", .before = 1L
    )
  ),
  pilot_parameter_failures = bind_rows(
    mutate(
      achieved_diagnostics$parameter_failures,
      exposure = "achieved", .before = 1L
    ),
    mutate(
      assignment_diagnostics$parameter_failures,
      exposure = "assignment_derived", .before = 1L
    )
  ),
  audit = list(
    reused_fits = 14L,
    achieved_pilot_models = 7L,
    assignment_pilot_models = 7L,
    all_pilot_diagnostics_pass = TRUE
  )
)

output_directory <- file.path("results", "checkpoints")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(
  output_directory, "exposure_pilot_completed.qs"
)
temporary_path <- paste0(output_path, ".tmp")
on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
qs::qsave(checkpoint, temporary_path, preset = "high")
if (!file.copy(temporary_path, output_path, overwrite = TRUE)) {
  stop("Could not install the exposure-pilot checkpoint.", call. = FALSE)
}
unlink(temporary_path, force = TRUE)

cat(
  normalizePath(output_path, winslash = "/", mustWork = TRUE),
  "\n",
  sep = ""
)
