# One reviewed stage per invocation; default tar_make_future clean R process.
args <- commandArgs(trailingOnly = TRUE)
stages <- c("prepare", "pilot", "fit", "diagnostics", "remedy", "loo", "exact", "correction", "combine", "report",
  "recover_prepare", "recover_pilot", "recover_fit", "recover_diagnostics", "recover_remedy", "recover_loo", "recover_exact", "recover_correction", "recover_accept", "recover_report")
if (length(args) != 1L || !args[[1L]] %in% stages) stop(
  "Usage: Rscript sensitivity/run_age_duration.R {prepare|pilot|fit|diagnostics|remedy|loo|exact|correction|combine|report}")
stopifnot(getRversion() == package_version("4.4.3"))
stage <- args[[1L]]
Sys.setenv(PROTEIN_INTERACTION_PHASE = if (stage %in% c("prepare", "diagnostics", "recover_prepare", "recover_diagnostics")) "none" else stage)
workers <- suppressWarnings(as.integer(Sys.getenv("PROTEIN_MODEL_WORKERS", "4")))
if (length(workers) != 1L || is.na(workers) || workers < 1L) stop("Invalid PROTEIN_MODEL_WORKERS.")
target_names <- switch(stage, prepare = "iad_preparation_files", pilot = "iad_pilot_audit_files",
  fit = "iad_fit_file", diagnostics = "iad_initial_audit_files", remedy = "iad_accepted_audit_file",
  loo = "iad_integrated_audit_files", exact = "iad_exact_audit_file", correction = "iad_final_loo_files",
  combine = "iad_combination_audit_files", report = "iad_supplement_file",
  recover_prepare = "iad_recovery_preparation_files", recover_pilot = "iad_recovery_pilot_audit_files",
  recover_fit = "iad_recovery_fit_file", recover_diagnostics = "iad_recovery_audit_files",
  recover_remedy = "iad_recovery_remedy_audit_files",
  recover_loo = "iad_recovery_loo_audit_files", recover_exact = "iad_recovery_exact_audit_file",
  recover_correction = "iad_recovery_final_loo_files", recover_accept = "iad_recovery_accepted_file",
  recover_report = "iad_recovery_reporting_audit_files")
do.call(targets::tar_make_future, list(names = target_names, workers = workers,
  script = "_targets_age_duration.R", store = "_targets_age_duration", reporter = "balanced"))
