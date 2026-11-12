# Run one reviewed checkpoint in the project's R 4.4.3 console/environment.
args <- commandArgs(trailingOnly = TRUE)
stages <- c(prepare = "tst_preparation_files", pilot = "tst_pilot_audit_files", fit = "tst_fit_file",
  diagnostics = "tst_initial_audit_files", remedy = "tst_accepted_audit_file", loo = "tst_loo_audit_files",
  exact = "tst_exact_audit_file", exact_remedy = "tst_exact_review_audit_file",
  exact_marginal = "tst_marginal_audit_file",
  correction = "tst_final_loo_files", combine = "tst_combination_files", report = "tst_supplement_file")
if (length(args) != 1L || !args[[1L]] %in% names(stages)) stop("Choose one stage: ", paste(names(stages), collapse = " | "))
stopifnot(getRversion() == package_version("4.4.3"))
stage <- args[[1L]]
Sys.setenv(PROTEIN_TRAINING_PHASE = if (stage %in% c("prepare", "diagnostics")) "none" else stage)
workers <- suppressWarnings(as.integer(Sys.getenv("PROTEIN_MODEL_WORKERS", "4")))
if (is.na(workers) || workers < 1L) stop("Invalid worker count.")
do.call(targets::tar_make_future, list(names = unname(stages[[stage]]), script = "_targets_training_status.R",
  store = "_targets_training_status", workers = workers, reporter = "balanced"))
