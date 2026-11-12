# Isolated opt-in workflow. Never source or modify the main pipeline/store.
library(targets)
workers <- suppressWarnings(as.integer(Sys.getenv("PROTEIN_MODEL_WORKERS", "4")))
if (length(workers) != 1L || is.na(workers) || workers < 1L) stop("Invalid PROTEIN_MODEL_WORKERS.")
future::plan(future::multisession, workers = workers)
tar_option_set(packages = c("here", "qs", "tidyverse", "metafor", "mgcv", "splines",
  "lspline", "brms", "rstan", "posterior", "loo", "matrixStats", "digest"),
  seed = 731988L, format = "qs", memory = "transient", garbage_collection = TRUE,
  storage = "worker", retrieval = "worker")
for (path in c("R/functions.R", "R/primary_model_diagnostics.R", "R/primary_model_loo.R",
  "R/primary_model_averaging.R", "R/incremental_update.R", "R/update_2026_absolute_exposures.R",
  "R/verification_artifacts.R", "R/age_duration_interactions.R", "R/age_duration_recovery.R",
  "R/age_duration_recovery_reporting.R", "R/age_duration_batch.R",
  "R/age_duration_family_targets.R", "R/age_duration_exact_remedies.R")) tar_source(path)

iad_cue <- function(phase) tar_cue(mode = if (iad_phase(phase)) "thorough" else "never")
c(list(
  tar_target(iad_source_files, iad_checkpoint_inputs(), format = "file"),
  tar_target(iad_frozen_checkpoints, iad_read_checkpoints(iad_source_files)),
  tar_target(iad_recovery_mapping_file, iad_recovery_map_input(), cue = tar_cue(mode = "always"), format = "file"),
  tar_target(iad_checkpoints, iad_apply_recovery(iad_frozen_checkpoints, iad_recovery_mapping_file)),
  tar_target(iad_recovery_reporting_gate_file, iad_recovery_reporting_gate(iad_recovery_mapping_file), format = "file"),
  tar_target(iad_prior_file, here::here("sensitivity", "age-duration-base-priors.csv"), format = "file"),
  tar_target(iad_base_priors, iad_priors(iad_prior_file)),
  tar_target(iad_specification, iad_specifications(iad_checkpoints), iteration = "list"),
  tar_target(iad_input_preflight, iad_preflight(iad_specification, iad_checkpoints, iad_base_priors)),
  tar_target(iad_profile_plan, iad_profiles(iad_checkpoints)),
  tar_target(iad_batch_storage_gate_file, iad_batch_storage_gate(iad_specification, iad_checkpoints),
    cue = tar_cue(mode = "always"), format = "file"),
  tar_target(iad_preparation_files, { iad_recovery_reporting_gate_file;
    iad_write_preparation(iad_checkpoints, iad_specification,
      iad_input_preflight, iad_profile_plan) }, format = "file"),
  tar_target(iad_pilot_specification, iad_pilot_specs(iad_specification), iteration = "list"),
  tar_target(iad_remaining_specification, Filter(function(s) !(s$interaction == "both" &&
    s$parent_index %in% c(6L, 26L)), iad_specification), iteration = "list"),
  tar_target(iad_pilot_file, { iad_recovery_reporting_gate_file;
    iad_fit(iad_pilot_specification, iad_checkpoints, iad_base_priors,
      iad_input_preflight, "pilot") }, pattern = map(iad_pilot_specification),
    cue = iad_cue("pilot"), format = "file", iteration = "list", error = "continue"),
  tar_target(iad_pilot_diagnostic, iad_diagnose(iad_pilot_file), pattern = map(iad_pilot_file), iteration = "list"),
  tar_target(iad_pilot_summary, iad_diagnostic_summary(iad_pilot_diagnostic)),
  tar_target(iad_pilot_audit_files, iad_diagnostic_files(iad_pilot_diagnostic, "pilot"), format = "file"),
  tar_target(iad_fit_file, iad_fit(iad_remaining_specification, iad_checkpoints, iad_base_priors,
    iad_input_preflight, "fit", iad_pilot_summary), pattern = map(iad_remaining_specification),
    cue = iad_cue("fit"), format = "file", iteration = "list", error = "continue"),
  tar_target(iad_initial_diagnostic, iad_diagnose(iad_fit_file), pattern = map(iad_fit_file), iteration = "list"),
  tar_target(iad_initial_summary, iad_diagnostic_summary(c(iad_pilot_diagnostic, iad_initial_diagnostic))),
  tar_target(iad_initial_audit_files, iad_diagnostic_files(c(iad_pilot_diagnostic,
    iad_initial_diagnostic), "initial"), format = "file"),
  tar_target(iad_remedy_file, here::here("sensitivity", "age-duration-remedies.csv"), format = "file"),
  tar_target(iad_exact_remedy_ledger_file, here::here("sensitivity", "age-duration-exact-remedies.csv"), format = "file"),
  tar_target(iad_remedy_specification, iad_remedy_specs(iad_specification, iad_remedy_file,
    iad_initial_summary), iteration = "list"),
  tar_target(iad_remedy_fit_file, iad_fit_remedy(iad_remedy_specification, iad_checkpoints,
    iad_base_priors, iad_input_preflight), pattern = map(iad_remedy_specification),
    cue = iad_cue("remedy"), format = "file", iteration = "list", error = "continue"),
  tar_target(iad_remedy_diagnostic, iad_diagnose_remedy(iad_remedy_fit_file),
    pattern = map(iad_remedy_fit_file), iteration = "list"),
  tar_target(iad_accepted_record, iad_accepted(iad_initial_summary,
    iad_diagnostic_summary(iad_remedy_diagnostic), iad_specification, iad_checkpoints), iteration = "list"),
  tar_target(iad_accepted_audit_file, iad_csv(dplyr::bind_rows(iad_accepted_record),
    iad_path("diagnostics", "accepted.csv")), format = "file"),
  tar_target(iad_integrated_loo_file, iad_compute_loo(iad_accepted_record),
    pattern = map(iad_accepted_record), cue = iad_cue("loo"), format = "file", iteration = "list", error = "continue"),
  tar_target(iad_integrated_audit, iad_loo_audit(iad_integrated_loo_file)),
  tar_target(iad_integrated_audit_files, iad_write_loo_audit(iad_integrated_audit), format = "file"),
  tar_target(iad_exact_record, iad_integrated_audit$cases, iteration = "list"),
  tar_target(iad_exact_file, iad_exact(iad_exact_record), pattern = map(iad_exact_record),
    cue = iad_cue("exact"), format = "file", iteration = "list", error = "continue"),
  tar_target(iad_exact_diagnostic, iad_exact_diagnostics(iad_exact_file),
    pattern = map(iad_exact_file), iteration = "list"),
  tar_target(iad_exact_audit_file, iad_csv(dplyr::bind_rows(iad_exact_diagnostic),
    iad_path("loo", "exact_diagnostics.csv")), format = "file"),
  tar_target(iad_exact_correction_file, iad_score_exact(iad_exact_file), pattern = map(iad_exact_file),
    cue = iad_cue("correction"), format = "file", iteration = "list", error = "continue"),
  tar_target(iad_final_loo_files, iad_finalise_loo(iad_integrated_audit, iad_exact_correction_file), format = "file"),
  tar_target(iad_combined_record, { iad_recovery_reporting_gate_file;
    iad_combine(iad_checkpoints, iad_accepted_record, iad_final_loo_files) },
    cue = iad_cue("combine"), iteration = "list"),
  tar_target(iad_combination_audit_files, iad_write_combination(iad_combined_record), format = "file"),
  tar_target(iad_report_files, iad_report_scope(iad_combined_record, iad_checkpoints, iad_profile_plan),
    pattern = map(iad_combined_record), cue = iad_cue("report"), format = "file", iteration = "list", error = "continue"),
  tar_target(iad_supplement_file, iad_write_supplement(iad_report_files, iad_combined_record),
    cue = iad_cue("report"), format = "file"),
  tar_target(iad_batch_supplement_file, iad_write_supplement(
    c(iad_relative_report_files, iad_total_report_files, iad_supplement_report_files),
    c(iad_relative_combined_record, iad_total_combined_record, iad_supplement_combined_record)),
    cue = iad_cue("report"), format = "file")
), iad_recovery_targets(), unlist(lapply(iad_families(), iad_family_targets), recursive = FALSE))
