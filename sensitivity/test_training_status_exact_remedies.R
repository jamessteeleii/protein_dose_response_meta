# Non-sampling selection, immutable-resume and sampler-only retry fixtures.
stopifnot(getRversion() == package_version("4.4.3"))
test_training_exact_remedies <- function() {
  env <- new.env(parent = .GlobalEnv)
  sys.source("R/age_duration_interactions.R", env)
  env$fit_primary_integrated_exact_loo_case <- function(...) stop("Stub required")
  sys.source("R/training_status_workflow.R", env)
  sys.source("R/training_status_exact_remedies.R", env)
  expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
  root <- tempfile("training-exact-retry-fixture-")
  env$iad_path <- function(..., archive = FALSE) file.path(root, ...)
  env$iad_require_phase <- function(phase) stopifnot(phase == "exact_remedy")
  fit <- env$iad_qsave(list(specification = list(iter_sampling = 6000L, adapt_delta = .999, max_treedepth = 15L),
    model = list(formula = "unchanged", prior = "unchanged")), file.path(root, "parent.qs"))
  set.seed(11); loo <- suppressWarnings(loo::loo(matrix(rnorm(1200, -1, .1), ncol = 3)))
  loo$diagnostics$pareto_k[] <- .9
  integrated <- env$iad_qsave(list(loo = loo), file.path(root, "integrated.qs"))
  records <- lapply(c(14L, 26L), function(index) list(family = "supplement", outcome = "strength",
    model_index = index, observation = 2L, fit_path = fit, fit_sha256 = env$iad_hash(fit),
    loo_path = integrated, loo_sha256 = env$iad_hash(integrated)))
  paths <- vapply(seq_along(records), function(i) env$iad_qsave(list(compact = TRUE,
    interaction_record = records[[i]], diagnostics = data.frame(family = "supplement", outcome = "strength",
      model_index = records[[i]]$model_index, observation = 2L, diagnostic_pass = i == 2L,
      iter_sampling = 6000L, iter_warmup = 2000L), correction = if(i == 2L) list(passed = TRUE) else NULL),
    file.path(root, paste0("initial_", i, ".qs"))), character(1))
  hashes <- vapply(paths, env$iad_hash, character(1))
  ledger <- data.frame(family = "supplement", outcome = "strength", model_index = 14L, observation = 2L,
    stage = "remedy_1", iter_warmup = 4000L, iter_sampling = 12000L, adapt_delta = .9995,
    max_treedepth = 15L, seed = 1988L, initial_sha256 = hashes[1L], rationale = "Fixture reviewed retry")
  ledger_path <- env$iad_csv(ledger, file.path(root, "ledger.csv"))
  plans <- env$tst_exact_remedy_plan(paths, ledger_path)
  calls <- 0L; pass <- TRUE
  env$tst_original_exact_fit <- function(job, chains, cores, warmup, iter, seed, iter_sampling) {
    calls <<- calls + 1L
    stopifnot(chains == 4L, cores == 4L, warmup == 4000L, iter == 16000L, seed == 1988L,
      iter_sampling == 12000L, job$fitted_candidate$specification$adapt_delta == .9995,
      job$fitted_candidate$specification$max_treedepth == 15L,
      identical(job$fitted_candidate$model, qs::qread(fit)$model))
    list()
  }
  env$tst_original_exact_diagnostics <- function(refit) data.frame(diagnostic_pass = pass)
  env$iad_build_exact_correction <- function(refit) list(exact_refit_summary = data.frame(diagnostic_pass = TRUE))
  stopifnot(identical(env$tst_review_exact(plans[[2L]]), paths[2L]), calls == 0L)
  retry <- env$tst_review_exact(plans[[1L]])
  stopifnot(calls == 1L, !retry %in% paths, qs::qread(retry)$diagnostics$iter_warmup == 4000L,
    identical(env$tst_review_exact(plans[[1L]]), retry), calls == 1L,
    identical(vapply(paths, env$iad_hash, character(1)), hashes))
  changed <- plans[[1L]]; changed$remedy$rationale <- "Different plan"
  expect_error(env$tst_review_exact(changed))
  pass <- FALSE; changed <- plans[[1L]]; changed$remedy$stage <- "remedy_2"
  failed <- qs::qread(env$tst_review_exact(changed))
  stopifnot(calls == 2L, is.null(failed$correction), identical(failed$diagnostics$diagnostic_pass, FALSE))
  expect_error(env$tst_exact_remedy_plan(paths[-1L], ledger_path))
  env$iad_csv(rbind(ledger, ledger), ledger_path)
  expect_error(env$tst_exact_remedy_plan(paths, ledger_path))
  cat("Training exact remedy selection, sampler-only schedule, resume and failure gates passed; no sampling.\n")
}
test_training_exact_remedies()
