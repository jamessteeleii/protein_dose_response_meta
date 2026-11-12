#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
  library(qs)
})

required_targets <- c(
  "sensitivity_within_full_contrast_data_by_outcome",
  "sensitivity_final_model_specifications",
  "sensitivity_within_full_accepted_fit_manifest",
  "sensitivity_within_full_final_diagnostic_audit",
  "sensitivity_within_full_integrated_loo_audit",
  "sensitivity_within_full_exact_loo_diagnostic_audit",
  "sensitivity_within_full_final_loo_audit",
  "sensitivity_within_full_stacking_weights"
)

missing_targets <- required_targets[!vapply(
  required_targets,
  targets::tar_exist_objects,
  logical(1)
)]
if (length(missing_targets)) {
  stop(
    "Cannot freeze the within-study checkpoint; cached targets are missing: ",
    paste(missing_targets, collapse = ", "),
    call. = FALSE
  )
}

analysis_data <- targets::tar_read_raw(
  "sensitivity_within_full_contrast_data_by_outcome"
)
specifications <- targets::tar_read_raw(
  "sensitivity_final_model_specifications"
)
fit_manifest <- targets::tar_read_raw(
  "sensitivity_within_full_accepted_fit_manifest"
)
sampling_audit <- targets::tar_read_raw(
  "sensitivity_within_full_final_diagnostic_audit"
)
integrated_audit <- targets::tar_read_raw(
  "sensitivity_within_full_integrated_loo_audit"
)
exact_audit <- targets::tar_read_raw(
  "sensitivity_within_full_exact_loo_diagnostic_audit"
)
final_loo_audit <- targets::tar_read_raw(
  "sensitivity_within_full_final_loo_audit"
)
stacking_weights <- targets::tar_read_raw(
  "sensitivity_within_full_stacking_weights"
)

if (
  !isTRUE(sampling_audit$all_sampling_diagnostics_pass) ||
    !isTRUE(exact_audit$all_exact_refit_diagnostics_pass) ||
    !isTRUE(final_loo_audit$stacking_ready) ||
    nrow(fit_manifest) != 104L || nrow(stacking_weights) != 104L ||
    any(fit_manifest$family != "within") ||
    any(!file.exists(fit_manifest$fit_path))
) {
  stop(
    "Cannot freeze the within-study checkpoint because its final audits fail.",
    call. = FALSE
  )
}

checkpoint <- list(
  checkpoint_version = 1L,
  family = "within",
  generated_utc = format(
    Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"
  ),
  analysis_data_by_outcome = analysis_data,
  model_specifications = specifications,
  accepted_fit_manifest = fit_manifest,
  sampling_diagnostic_summary = sampling_audit$diagnostic_summary,
  sampling_parameter_failures = sampling_audit$parameter_failures,
  integrated_loo_model_summary = integrated_audit$model_summary,
  integrated_loo_outcome_summary = integrated_audit$outcome_summary,
  initial_high_pareto_k_cases = integrated_audit$high_k_cases,
  exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
  exact_loo_correction_summary = final_loo_audit$exact_refit_summary,
  final_loo_model_summary = final_loo_audit$model_summary,
  stacking_weights = stacking_weights,
  audit = list(
    models = sampling_audit$models,
    all_sampling_diagnostics_pass =
      sampling_audit$all_sampling_diagnostics_pass,
    exact_refits = exact_audit$exact_refits,
    all_exact_refit_diagnostics_pass =
      exact_audit$all_exact_refit_diagnostics_pass,
    stacking_ready = final_loo_audit$stacking_ready
  )
)

output_directory <- file.path("results", "checkpoints")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(
  output_directory, "within_completed_analysis.qs"
)
temporary_path <- paste0(output_path, ".tmp")
on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
qs::qsave(checkpoint, temporary_path, preset = "high")
if (!file.copy(temporary_path, output_path, overwrite = TRUE)) {
  stop("Could not install the within-study checkpoint.", call. = FALSE)
}
unlink(temporary_path, force = TRUE)

cat(
  normalizePath(output_path, winslash = "/", mustWork = TRUE),
  "\n",
  sep = ""
)
