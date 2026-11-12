# Focused stub tests of schedule forwarding and compact save/resume/failure.
# No actual fitting, targets execution or use of analytical records.
stopifnot(getRversion() == package_version("4.4.3"))
test_training_exact <- function() {
  env <- new.env(parent = .GlobalEnv)
  sys.source("R/age_duration_interactions.R", env)
  env$fit_primary_integrated_exact_loo_case <- function(...) stop("Must be stubbed")
  sys.source("R/training_status_workflow.R", env)
  calls <- 0L
  env$tst_original_exact_fit <- function(job, chains, cores, warmup, iter, seed, iter_sampling) {
    calls <<- calls + 1L
    stopifnot(chains == 4L, cores == 4L, warmup == 2000L, iter == warmup + iter_sampling,
      iter_sampling == job$fitted_candidate$specification$iter_sampling, seed == 1988L)
    list()
  }
  for (n in c(6000L, 12000L)) {
    result <- env$fit_primary_integrated_exact_loo_case(list(fitted_candidate = list(specification = list(iter_sampling = n))))
    stopifnot(result$training_exact_schedule$iter_sampling == n)
  }
  root <- tempfile("training-exact-fixture-")
  env$iad_path <- function(..., archive = FALSE) file.path(root, ...)
  env$iad_require_phase <- function(phase) stopifnot(phase == "exact")
  fit <- env$iad_qsave(list(specification = list(iter_sampling = 12000L)), file.path(root, "parent.qs"))
  set.seed(11)
  loo <- suppressWarnings(loo::loo(matrix(rnorm(1200, -1, .1), ncol = 3)))
  loo$diagnostics$pareto_k[] <- .9
  integrated <- env$iad_qsave(list(loo = loo), file.path(root, "integrated.qs"))
  record <- list(family = "supplement", outcome = "strength", model_index = 14L, observation = 2L,
    fit_path = fit, fit_sha256 = env$iad_hash(fit), loo_path = integrated, loo_sha256 = env$iad_hash(integrated))
  pass <- TRUE
  env$tst_original_exact_diagnostics <- function(refit) data.frame(diagnostic_pass = pass)
  env$iad_build_exact_correction <- function(refit) list(score_fixture = TRUE)
  path <- env$iad_exact(record, compact = TRUE)
  saved <- qs::qread(path)
  stopifnot(calls == 3L, saved$compact, saved$diagnostics$iter_sampling == 12000L,
    saved$diagnostics$iter_warmup == 2000L, saved$correction$score_fixture,
    is.null(saved$refitted_model), identical(env$iad_exact(record, compact = TRUE), path), calls == 3L)
  pass <- FALSE; record$observation <- 3L
  failed <- qs::qread(env$iad_exact(record, compact = TRUE))
  stopifnot(calls == 4L, identical(failed$diagnostics$diagnostic_pass, FALSE), is.null(failed$correction))
  cat("Exact schedule forwarding and compact checkpoint resume/failure tests passed; no sampling.\n")
}
test_training_exact()
