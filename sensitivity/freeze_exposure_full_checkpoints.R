#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
  library(qs)
})

common_targets <- c(
  "sensitivity_exposure_full_data_by_outcome",
  "sensitivity_final_model_specifications"
)
family_targets <- function(family) {
  c(
    paste0("sensitivity_", family, "_full_accepted_fit_manifest"),
    paste0("sensitivity_", family, "_full_final_diagnostic_audit"),
    paste0("sensitivity_", family, "_full_integrated_loo_audit"),
    paste0(
      "sensitivity_", family,
      "_full_exact_loo_diagnostic_audit"
    ),
    paste0("sensitivity_", family, "_full_final_loo_audit"),
    paste0("sensitivity_", family, "_full_stacking_weights")
  )
}

required_targets <- c(
  common_targets,
  family_targets("achieved"),
  family_targets("assignment")
)
missing_targets <- required_targets[!vapply(
  required_targets,
  targets::tar_exist_objects,
  logical(1)
)]
if (length(missing_targets)) {
  stop(
    "Cannot freeze the exposure checkpoints; cached targets are missing: ",
    paste(missing_targets, collapse = ", "),
    call. = FALSE
  )
}

exposure_data <- targets::tar_read_raw(
  "sensitivity_exposure_full_data_by_outcome"
)
specifications <- targets::tar_read_raw(
  "sensitivity_final_model_specifications"
)
output_directory <- file.path("results", "checkpoints")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

write_family_checkpoint <- function(family) {
  target_names <- family_targets(family)
  objects <- setNames(
    lapply(target_names, targets::tar_read_raw),
    c(
      "fit_manifest", "sampling_audit", "integrated_audit",
      "exact_audit", "final_loo_audit", "stacking_weights"
    )
  )
  if (
    !isTRUE(objects$sampling_audit$all_sampling_diagnostics_pass) ||
      !isTRUE(objects$exact_audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(objects$final_loo_audit$stacking_ready) ||
      nrow(objects$fit_manifest) != 104L ||
      nrow(objects$stacking_weights) != 104L ||
      any(objects$fit_manifest$family != family) ||
      any(objects$stacking_weights$family != family) ||
      any(!file.exists(objects$fit_manifest$fit_path))
  ) {
    stop(
      "Cannot freeze the ", family,
      " checkpoint because its final audits fail.",
      call. = FALSE
    )
  }
  checkpoint <- list(
    checkpoint_version = 1L,
    family = family,
    generated_utc = format(
      Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"
    ),
    analysis_data_by_outcome = exposure_data[[family]],
    model_specifications = specifications,
    accepted_fit_manifest = objects$fit_manifest,
    sampling_diagnostic_summary =
      objects$sampling_audit$diagnostic_summary,
    sampling_parameter_failures =
      objects$sampling_audit$parameter_failures,
    integrated_loo_model_summary = objects$integrated_audit$model_summary,
    integrated_loo_outcome_summary =
      objects$integrated_audit$outcome_summary,
    initial_high_pareto_k_cases = objects$integrated_audit$high_k_cases,
    exact_refit_diagnostic_summary = objects$exact_audit$diagnostic_summary,
    exact_loo_correction_summary =
      objects$final_loo_audit$exact_refit_summary,
    final_loo_model_summary = objects$final_loo_audit$model_summary,
    stacking_weights = objects$stacking_weights,
    audit = list(
      models = objects$sampling_audit$models,
      all_sampling_diagnostics_pass =
        objects$sampling_audit$all_sampling_diagnostics_pass,
      exact_refits = objects$exact_audit$exact_refits,
      all_exact_refit_diagnostics_pass =
        objects$exact_audit$all_exact_refit_diagnostics_pass,
      stacking_ready = objects$final_loo_audit$stacking_ready
    )
  )
  output_path <- file.path(
    output_directory,
    paste0(family, "_completed_analysis.qs")
  )
  temporary_path <- paste0(output_path, ".tmp")
  on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
  qs::qsave(checkpoint, temporary_path, preset = "high")
  if (!file.copy(temporary_path, output_path, overwrite = TRUE)) {
    stop("Could not install the ", family, " checkpoint.", call. = FALSE)
  }
  unlink(temporary_path, force = TRUE)
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}

paths <- vapply(
  c("achieved", "assignment"),
  write_family_checkpoint,
  character(1)
)
cat(paste(paths, collapse = "\n"), "\n", sep = "")
