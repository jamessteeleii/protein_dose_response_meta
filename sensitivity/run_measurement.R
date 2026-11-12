# Reviewed-stage runner: invoke one stage at a time and inspect its checkpoint
# before enabling the next stage. Every tar_make_future() below uses its
# default clean child R process.
args <- commandArgs(trailingOnly = TRUE)
stages <- c("prepare", "pilot", "fit", "diagnostics", "remedy", "loo", "exact", "report")
if (length(args) != 1L || !args[[1L]] %in% stages) {
  stop("Usage: Rscript sensitivity/run_measurement.R {prepare|pilot|fit|diagnostics|remedy|loo|exact|report}",
       call. = FALSE)
}
stage <- args[[1L]]
suffix <- switch(stage,
  prepare = "preparation_files", pilot = "pilot_audit_file", fit = "fit_file",
  diagnostics = "initial_diagnostic_files", remedy = "remedy_file",
  loo = "integrated_audit_files", exact = "final_loo_files",
  report = c("verification_files", "figure_files"))
Sys.setenv(PROTEIN_MEASUREMENT_PHASES = if (stage %in% c("prepare", "diagnostics")) "none" else stage)
workers <- suppressWarnings(as.integer(Sys.getenv("PROTEIN_MODEL_WORKERS", unset = "4")))
if (length(workers) != 1L || is.na(workers) || workers < 1L) stop("Invalid PROTEIN_MODEL_WORKERS.")
target_names <- unlist(lapply(c("prereg_methods", "direct_only"), function(variant) {
  paste0("measurement_", variant, "_", suffix)
}), use.names = FALSE)
do.call(
  targets::tar_make_future,
  list(
    names = target_names,
    workers = workers,
    reporter = "balanced"
  )
)
