#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
  library(qs)
})

target_names <- c(
  analysis_data = "sensitivity_primary_data_by_outcome",
  specifications = "sensitivity_final_model_specifications",
  fit_manifest = "sensitivity_rescor_full_accepted_fit_manifest",
  sampling_audit =
    "sensitivity_rescor_full_fully_resolved_diagnostic_audit",
  study_level_audit = "sensitivity_rescor_full_integrated_loo_audit",
  integrated_audit =
    "sensitivity_rescor_full_observation_loo_audit",
  exact_audit =
    "sensitivity_rescor_full_observation_exact_fully_resolved_diagnostic_audit",
  final_loo_audit =
    "sensitivity_rescor_full_final_observation_loo_audit",
  stacking_weights = "sensitivity_rescor_full_stacking_weights"
)

missing_targets <- unname(target_names)[!vapply(
  unname(target_names),
  targets::tar_exist_objects,
  logical(1)
)]
if (length(missing_targets)) {
  stop(
    "Cannot freeze the residual-correlation checkpoint; cached targets are missing: ",
    paste(missing_targets, collapse = ", "),
    call. = FALSE
  )
}

objects <- setNames(
  lapply(unname(target_names), targets::tar_read_raw),
  names(target_names)
)
if (
  !isTRUE(objects$sampling_audit$all_sampling_diagnostics_pass) ||
    !isTRUE(objects$exact_audit$all_exact_refit_diagnostics_pass) ||
    !isTRUE(objects$final_loo_audit$stacking_ready) ||
    nrow(objects$fit_manifest) != 104L ||
    nrow(objects$stacking_weights) != 104L ||
    any(objects$fit_manifest$family != "rescor") ||
    any(objects$stacking_weights$family != "rescor") ||
    any(!file.exists(objects$fit_manifest$fit_path))
) {
  stop(
    "Cannot freeze the residual-correlation checkpoint because its final audits fail.",
    call. = FALSE
  )
}

checkpoint <- list(
  checkpoint_version = 1L,
  family = "rescor",
  generated_utc = format(
    Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"
  ),
  analysis_data_by_outcome = objects$analysis_data,
  model_specifications = objects$specifications,
  accepted_fit_manifest = objects$fit_manifest,
  sampling_diagnostic_summary =
    objects$sampling_audit$diagnostic_summary,
  sampling_parameter_failures =
    objects$sampling_audit$parameter_failures,
  integrated_loo_model_summary = objects$integrated_audit$model_summary,
  integrated_loo_outcome_summary = objects$integrated_audit$outcome_summary,
  initial_high_pareto_k_cases = objects$integrated_audit$high_k_cases,
  study_level_loo_model_summary = objects$study_level_audit$model_summary,
  study_level_loo_outcome_summary =
    objects$study_level_audit$outcome_summary,
  study_level_high_pareto_k_cases =
    objects$study_level_audit$high_k_cases,
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
    stacking_ready = objects$final_loo_audit$stacking_ready,
    loo_unit = "observation_conditional_on_study",
    study_level_loo_retained_as_influence_audit = TRUE
  )
)

output_directory <- file.path("results", "checkpoints")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(
  output_directory, "rescor_completed_analysis.qs"
)
temporary_path <- paste0(output_path, ".tmp")
on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
qs::qsave(checkpoint, temporary_path, preset = "high")
if (!file.copy(temporary_path, output_path, overwrite = TRUE)) {
  stop(
    "Could not install the completed residual-correlation checkpoint.",
    call. = FALSE
  )
}
unlink(temporary_path, force = TRUE)
cat(normalizePath(output_path, winslash = "/", mustWork = TRUE), "\n")
