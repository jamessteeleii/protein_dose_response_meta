library(targets)

prior_workers <- suppressWarnings(as.integer(Sys.getenv("PROTEIN_MODEL_WORKERS", "4")))
if (is.na(prior_workers) || prior_workers < 1L) {
  stop("PROTEIN_MODEL_WORKERS must be a positive integer.")
}
future::plan(future::multisession, workers = prior_workers)

tar_option_set(
  packages = c("brms", "rstan", "lspline", "splines"),
  memory = "transient",
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker"
)
tar_source("R/functions.R")
tar_source("R/primary_model_averaging.R")
tar_source("R/prior_sensitivity.R")

list(
  tar_target(prior_bundle_files, prior_sensitivity_bundle_files(), format = "file"),
  tar_target(prior_inputs, read_prior_sensitivity_inputs(prior_bundle_files)),
  tar_target(prior_registry, make_prior_sensitivity_candidates(prior_inputs)),
  tar_target(prior_candidate,
             lapply(seq_len(nrow(prior_registry)), function(i) prior_registry[i, ]),
             iteration = "list"),
  tar_target(prior_fit_file, prior_sensitivity_fit_path(prior_candidate),
             pattern = map(prior_candidate), format = "file", iteration = "list"),
  tar_target(prior_result,
             analyse_prior_sensitivity_candidate(prior_candidate, prior_fit_file, prior_inputs),
             pattern = map(prior_candidate, prior_fit_file), iteration = "list"),
  tar_target(prior_report_files,
             write_prior_sensitivity_report(prior_result, prior_registry), format = "file")
)
