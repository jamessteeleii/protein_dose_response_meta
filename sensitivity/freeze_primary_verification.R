#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(targets)
})

source(file.path("R", "verification_artifacts.R"))

read_target <- function(name) {
  targets::tar_read_raw(name)
}

diagnostic_audit <- read_target("primary_model_diagnostic_audit")
loo_audit <- read_target("primary_psis_loo_final_audit")
stacking_weights <- read_target("primary_stacking_weights_summary")

if (nrow(diagnostic_audit$diagnostic_summary) != 104L) {
  stop("Expected 104 primary-model diagnostic rows.", call. = FALSE)
}
if (!all(diagnostic_audit$diagnostic_summary$diagnostic_pass %in% TRUE)) {
  stop("At least one accepted primary model fails its sampling diagnostics.", call. = FALSE)
}
if (nrow(loo_audit$model_summary) != 104L) {
  stop("Expected 104 final primary LOO model summaries.", call. = FALSE)
}
if (any(loo_audit$outcome_summary$n_models_with_pareto_k_warning > 0L) ||
    any(loo_audit$outcome_summary$total_pareto_k_warnings > 0L) ||
    !all(loo_audit$outcome_summary$stacking_ready %in% TRUE)) {
  stop("The final primary LOO audit still contains high Pareto-k cases.", call. = FALSE)
}
weight_totals <- stacking_weights |>
  group_by(outcome) |>
  summarise(models = n(), weight = sum(stacking_weight), .groups = "drop")
if (any(weight_totals$models != 26L) ||
    any(abs(weight_totals$weight - 1) > 1e-8)) {
  stop("Primary stacking weights do not contain 26 models summing to one.", call. = FALSE)
}

files <- write_verification_bundle(
  family = "primary",
  artifacts = list(
    analysis_data = bind_rows(
      whole_body = read_target("protein_data_complete_whole_body"),
      direct_muscle = read_target("protein_data_complete_hypertrophy"),
      strength = read_target("protein_data_complete_strength"),
      physical_function = read_target("protein_data_complete_function"),
      .id = "analysis_outcome"
    ),
    model_specifications = read_target("primary_model_fit_specification_table"),
    diagnostic_summary = diagnostic_audit$diagnostic_summary,
    diagnostic_classification = diagnostic_audit$definitive_classification,
    diagnostic_parameter_failures = diagnostic_audit$parameter_failures,
    posterior_predictive_observation_flags = diagnostic_audit$ppc_observation_flags,
    posterior_predictive_outcome_summary = diagnostic_audit$ppc_outcome_summary,
    loo_results = loo_audit$results,
    loo_model_summary = loo_audit$model_summary,
    loo_outcome_summary = loo_audit$outcome_summary,
    loo_method_comparison = loo_audit$method_comparison,
    loo_observation_diagnostics = loo_audit$observation_diagnostics,
    loo_exact_refit_summary = loo_audit$exact_refit_summary,
    stacking_weights = stacking_weights,
    dose_response_predictions = bind_rows(
      whole_body = read_target("preds_new_whole_body"),
      direct_muscle = read_target("preds_new_hypertrophy"),
      strength = read_target("preds_new_strength"),
      physical_function = read_target("preds_new_function"),
      .id = "analysis_outcome"
    ),
    pairwise_contrasts = bind_rows(
      whole_body = read_target("primary_heatmap_data_whole_body"),
      direct_muscle = read_target("primary_heatmap_data_hypertrophy"),
      strength = read_target("primary_heatmap_data_strength"),
      physical_function = read_target("primary_heatmap_data_function"),
      .id = "analysis_outcome"
    ),
    local_equivalence = read_target("primary_local_equivalence_summary"),
    reference_values = read_target("primary_reference_values")
  ),
  descriptions = c(
    analysis_data = "Complete-case rows supplied to the primary models.",
    model_specifications = "Final candidate-model registry and remedies.",
    diagnostic_summary = "Model-level convergence and sampling diagnostics.",
    diagnostic_classification = "Definitive diagnostic classification by model.",
    diagnostic_parameter_failures = "Parameters failing diagnostic thresholds.",
    posterior_predictive_observation_flags = "Observation-level PPC flags.",
    posterior_predictive_outcome_summary = "Outcome-level PPC summary.",
    loo_results = "Corrected integrated LOO objects used for stacking.",
    loo_model_summary = "Model-level corrected LOO summary.",
    loo_outcome_summary = "Outcome-level corrected LOO summary.",
    loo_method_comparison = "Conditional and integrated LOO comparison.",
    loo_observation_diagnostics = "Pointwise Pareto-k and LOO diagnostics.",
    loo_exact_refit_summary = "Exact case-deletion refit audit.",
    stacking_weights = "Final outcome-specific stacking weights.",
    dose_response_predictions = "Publication dose-response summaries.",
    pairwise_contrasts = "Publication pairwise contrasts and SESOI probabilities.",
    local_equivalence = "Local practical-equivalence summaries.",
    reference_values = "Raw-unit references for the SESOI."
  )
)

cat("Primary verification bundle written and validated:\n")
cat(paste0("- ", files, collapse = "\n"), "\n")
