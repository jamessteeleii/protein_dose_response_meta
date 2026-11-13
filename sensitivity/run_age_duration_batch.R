# Source in the project R console; only an explicit function call runs targets.
run_age_duration_initial_batch <- function(workers = 4L, prune_completed_strategies = FALSE) {
  stopifnot(getRversion() == package_version("4.4.3"))
  if (length(workers) != 1L || is.na(workers) || workers < 1L || workers %% 1 != 0 ||
      length(prune_completed_strategies) != 1L || !is.logical(prune_completed_strategies) ||
      is.na(prune_completed_strategies)) stop("Invalid initial-batch options.")
  if (isTRUE(prune_completed_strategies)) {
    # Author-approved choices only. Source locally so the pruning helpers do
    # not replace anything in the user's analysis session. Verify both plans
    # before deletion; the execution helper rechecks each plan immediately.
    source("sensitivity/age-duration-storage.R", local = TRUE)
    keys <- c("supplement_strategy", "total_strategy")
    plans <- lapply(keys, iad_storage_plan, verify = TRUE)
    print(do.call(rbind, lapply(plans, `[[`, "summary")), row.names = FALSE)
    for (key in keys) iad_prune_completed_family(key, execute = TRUE)
  }
  run_age_duration_checkpoint("initial_fit", families = c("relative", "total", "supplement"),
    workers = workers)
}

run_age_duration_checkpoint <- function(stage = "initial_fit",
    families = c("relative", "total", "supplement"), workers = 4L) {
  stopifnot(getRversion() == package_version("4.4.3"))
  allowed <- c("relative", "total", "supplement")
  stages <- c("prepare", "initial_fit", "diagnostics", "remedy", "loo", "exact", "exact_remedy", "correction", "combine", "report", "supplement")
  if (length(stage) != 1L || !stage %in% stages || !length(families) ||
      anyDuplicated(families) || !all(families %in% allowed) ||
      length(workers) != 1L || is.na(workers) || workers < 1L || workers %% 1 != 0) {
    stop("Choose one documented checkpoint and relative, total and/or supplement dose families.")
  }
  Sys.setenv(PROTEIN_MODEL_WORKERS = as.character(workers),
    PROTEIN_INTERACTION_FAMILIES = paste(families, collapse = ","))
  run <- function(targets, phase) {
    Sys.setenv(PROTEIN_INTERACTION_PHASE = phase)
    targets::tar_make_future(names = tidyselect::all_of(targets),
      script = "_targets_age_duration.R", store = "_targets_age_duration",
      workers = as.integer(workers), reporter = "balanced")
    progress <- targets::tar_progress(names = tidyselect::all_of(targets), store = "_targets_age_duration")
    if (!setequal(progress$name, targets) ||
        any(!progress$progress %in% c("completed", "skipped"))) {
      stop("Checkpoint incomplete: review errored or blocked targets. Successful dose-family fits remain saved; no later stage was started.")
    }
  }
  on.exit(Sys.setenv(PROTEIN_INTERACTION_PHASE = "none"), add = TRUE)
  if (stage %in% c("prepare", "initial_fit")) {
    # Recovery reporting is required only when accepted replacement fits exist.
    recovery_manifest <- file.path(
      "results", "update_2026", "age_duration", "recovery", "accepted_manifest.csv"
    )
    if (file.exists(recovery_manifest)) {
      run("iad_recovery_reporting_audit_files", "recover_report")
    }
    run(c("iad_preparation_files", "iad_batch_storage_gate_file"), "none")
    if (stage == "prepare") return(invisible(TRUE))
  }
  if (stage == "supplement") {
    run("iad_batch_supplement_file", "report")
    return(invisible(TRUE))
  }
  suffix <- switch(stage, initial_fit = "initial_fit_checkpoint_file", diagnostics = "initial_audit_files",
    remedy = "accepted_audit_file", loo = "integrated_audit_files", exact = "exact_audit_file",
    exact_remedy = "exact_review_audit_file", correction = "final_loo_files", combine = "combination_audit_files", report = "report_checkpoint_file")
  phase <- switch(stage, initial_fit = "batch_fit", diagnostics = "none", stage)
  run(paste0("iad_", families, "_", suffix), phase)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || length(args) > 2L) stop(
    "Usage: Rscript sensitivity/run_age_duration_batch.R <checkpoint> [relative,total,supplement]")
  families <- if (length(args) == 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else c("relative", "total", "supplement")
  run_age_duration_checkpoint(args[[1L]], families,
    as.integer(Sys.getenv("PROTEIN_MODEL_WORKERS", "4")))
}
