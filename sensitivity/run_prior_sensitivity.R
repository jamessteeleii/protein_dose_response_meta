# Launch only the saved-fit prior sensitivity workflow from the project root.
stopifnot(getRversion() == package_version("4.4.3"))
workers <- suppressWarnings(as.integer(Sys.getenv("PROTEIN_MODEL_WORKERS", "4")))
if (is.na(workers) || workers < 1L) stop("Invalid PROTEIN_MODEL_WORKERS.")
targets::tar_make_future(
  script = "_targets_prior_sensitivity.R",
  store = "_targets_prior_sensitivity",
  workers = workers,
  reporter = "balanced"
)
