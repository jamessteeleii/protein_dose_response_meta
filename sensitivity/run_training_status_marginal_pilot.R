# R 4.4.3, in the project root. --check runs only numerical/Stan syntax checks.
stopifnot(getRversion() == package_version("4.4.3"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("--check", "--fit")) stop("Choose --check or --fit explicitly.")
suppressPackageStartupMessages({library(brms); library(mgcv)})
for (file in c("R/primary_model_loo.R", "R/age_duration_interactions.R",
  "R/training_status_workflow.R", "R/training_status_marginal_exact.R")) source(file)
name <- "model_26_observation_0080.qs"
initial <- qs::qread(iad_path("supplement", "strength", "loo", "exact_checkpoints", name, archive = TRUE))
stopifnot(isTRUE(initial$compact), identical(initial$diagnostics$diagnostic_pass, FALSE))
input <- tst_prepare_marginal_exact(initial$interaction_record)
syntax <- rstan::stanc(model_code = input$code, model_name = "training_marginal_exact_26")
stopifnot(isTRUE(syntax$status))
cat("Model 26, observation 80. Maximum held-out log-density discrepancy:", input$score_error, "\n")
if (args == "--fit") {
  tst_age_duration_gate(); iad_require_measurement()
  output <- iad_path("loo", "marginal_pilot", name)
  summary <- iad_path("loo", "marginal_pilot", "diagnostics.csv")
  code_sha256 <- digest::digest(input$code, algo = "sha256")
  if (file.exists(output)) {
    result <- qs::qread(output)
    stopifnot(identical(result$interaction_record, input$record), identical(result$code_sha256, code_sha256))
  } else {
    result <- tst_fit_marginal_exact(input, 6000L, syntax)
    iad_qsave(result, output)
  }
  iad_csv(result$diagnostics, summary)
  print(result$diagnostics)
}
