#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
  library(qs)
})

required_targets <- c(
  "sensitivity_strategy_within_full_data_by_outcome",
  "sensitivity_final_model_specifications",
  "sensitivity_strategy_within_full_accepted_fit_manifest",
  "sensitivity_strategy_within_full_final_diagnostic_audit",
  "sensitivity_strategy_within_full_integrated_loo_audit",
  "sensitivity_strategy_within_full_exact_loo_diagnostic_audit",
  "sensitivity_strategy_within_full_final_loo_audit",
  "sensitivity_strategy_within_full_stacking_weights"
)

missing_targets <- required_targets[!vapply(
  required_targets,
  targets::tar_exist_objects,
  logical(1)
)]
if (length(missing_targets)) {
  stop(
    "Cannot freeze the within-strategy checkpoint; cached targets are missing: ",
    paste(missing_targets, collapse = ", "),
    call. = FALSE
  )
}

objects <- setNames(
  lapply(required_targets, targets::tar_read_raw),
  c(
    "analysis_data", "model_specifications", "fit_manifest",
    "sampling_audit", "integrated_audit", "exact_audit",
    "final_loo_audit", "stacking_weights"
  )
)
if (
  !isTRUE(objects$sampling_audit$all_sampling_diagnostics_pass) ||
    !isTRUE(objects$exact_audit$all_exact_refit_diagnostics_pass) ||
    !isTRUE(objects$final_loo_audit$stacking_ready) ||
    nrow(objects$fit_manifest) != 104L ||
    nrow(objects$stacking_weights) != 104L ||
    any(objects$fit_manifest$family != "strategy_within") ||
    any(objects$stacking_weights$family != "strategy_within") ||
    any(!file.exists(objects$fit_manifest$fit_path))
) {
  stop(
    "Cannot freeze the within-strategy checkpoint because a final audit fails.",
    call. = FALSE
  )
}

checkpoint <- list(
  checkpoint_version = 1L,
  family = "strategy_within",
  generated_utc = format(
    Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"
  ),
  analysis_data_by_outcome = objects$analysis_data,
  model_specifications = objects$model_specifications,
  accepted_fit_manifest = objects$fit_manifest,
  sampling_diagnostic_summary =
    objects$sampling_audit$diagnostic_summary,
  sampling_parameter_failures =
    objects$sampling_audit$parameter_failures,
  integrated_loo_model_summary = objects$integrated_audit$model_summary,
  integrated_loo_outcome_summary = objects$integrated_audit$outcome_summary,
  initial_high_pareto_k_cases = objects$integrated_audit$high_k_cases,
  exact_refit_diagnostic_summary = objects$exact_audit$diagnostic_summary,
  exact_loo_correction_summary = objects$final_loo_audit$exact_refit_summary,
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

output_directory <- file.path("results", "checkpoints")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(
  output_directory, "strategy_within_completed_analysis.qs"
)
temporary_path <- paste0(output_path, ".tmp")
on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
qs::qsave(checkpoint, temporary_path, preset = "high")
if (!file.copy(temporary_path, output_path, overwrite = TRUE)) {
  stop("Could not install the within-strategy checkpoint.", call. = FALSE)
}
unlink(temporary_path, force = TRUE)

cat(
  normalizePath(output_path, winslash = "/", mustWork = TRUE),
  "\n",
  sep = ""
)
