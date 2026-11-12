# Synthetic selection/resume tests with a stub fitter. No sampling or pipeline.
stopifnot(getRversion() == package_version("4.4.3"))
source("R/age_duration_interactions.R")
source("R/age_duration_exact_remedies.R")
expect_error <- function(x) stopifnot(inherits(tryCatch({force(x); NULL}, error = identity), "error"))

test_exact_remedies <- function() {
  env <- new.env(parent = .GlobalEnv)
  sys.source("R/age_duration_interactions.R", envir = env)
  sys.source("R/age_duration_exact_remedies.R", envir = env)
  root <- tempfile("iad-exact-remedies-")
  env$iad_path <- function(..., archive = FALSE) file.path(root, ...)
  env$iad_require_phase <- function(phase) stopifnot(phase == "exact_remedy")
  fit <- env$iad_qsave(list(fixture = "parent"), file.path(root, "parent.qs"))
  set.seed(1988)
  loo <- suppressWarnings(loo::loo(matrix(rnorm(1200, -1, 0.1), ncol = 3)))
  loo$diagnostics$pareto_k[] <- 0.9
  integrated <- env$iad_qsave(list(loo = loo), file.path(root, "integrated.qs"))
  records <- lapply(44:45, function(index) list(family = "relative", outcome = "function",
    model_index = as.integer(index), observation = 2L, fit_path = fit, fit_sha256 = env$iad_hash(fit),
    loo_path = integrated, loo_sha256 = env$iad_hash(integrated)))
  checkpoints <- lapply(seq_along(records), function(i) list(compact = TRUE,
    interaction_record = records[[i]], correction = NULL,
    diagnostics = tibble::tibble(family = "relative", outcome = "function", model_index = records[[i]]$model_index,
      observation = 2L, geometry_pass = TRUE, convergence_pass = i == 2L, diagnostic_pass = i == 2L)))
  paths <- vapply(seq_along(checkpoints), function(i) env$iad_qsave(checkpoints[[i]],
    file.path(root, paste0("initial_", i, ".qs"))), character(1))
  initial_hashes <- vapply(paths, env$iad_hash, character(1))
  ledger <- tibble::tibble(family = "relative", outcome = "function", model_index = 44L, observation = 2L,
    stage = "remedy_1", iter_sampling = 12000L, initial_sha256 = initial_hashes[[1L]], rationale = "Convergence retry")
  ledger_path <- env$iad_csv(ledger, file.path(root, "ledger.csv"))
  plans <- env$iad_exact_remedy_plan(paths, ledger_path, "relative")
  stopifnot(length(plans) == 2L, !is.null(plans[[1L]]$remedy), is.null(plans[[2L]]$remedy))
  calls <- 0L
  pass <- TRUE
  env$fit_primary_integrated_exact_loo_case <- function(job, chains, cores, warmup, iter, seed, iter_sampling) {
    calls <<- calls + 1L
    stopifnot(chains == 4L, cores == 4L, warmup == 2000L, iter == 14000L, seed == 1988L,
      iter_sampling == 12000L, identical(job$fitted_candidate, qs::qread(fit)))
    list(job = job)
  }
  env$iad_exact_diagnostic_summary <- function(refit) {
    tibble::tibble(family = "relative", outcome = "function", model_index = 44L,
      observation = 2L, diagnostic_pass = pass)
  }
  env$iad_build_exact_correction <- function(refit) list(exact_refit_summary = tibble::tibble(
    observation = 2L, diagnostic_pass = TRUE, exact_elpd_loo = -2))
  stopifnot(identical(env$iad_review_exact(plans[[2L]]), paths[[2L]]), calls == 0L)
  retry <- env$iad_review_exact(plans[[1L]])
  stopifnot(calls == 1L, !retry %in% paths, isTRUE(qs::qread(retry)$diagnostics$diagnostic_pass),
    qs::qread(retry)$correction$exact_refit_summary$iter_sampling == 12000L,
    identical(vapply(paths, env$iad_hash, character(1)), initial_hashes))
  stopifnot(identical(env$iad_review_exact(plans[[1L]]), retry), calls == 1L)
  changed <- plans[[1L]]; changed$remedy$rationale <- "Changed decision"
  expect_error(env$iad_review_exact(changed))
  failed_plan <- plans[[1L]]; failed_plan$remedy$stage <- "remedy_2"; pass <- FALSE
  failed <- qs::qread(env$iad_review_exact(failed_plan))
  stopifnot(is.null(failed$correction), identical(failed$diagnostics$diagnostic_pass, FALSE), calls == 2L)
  expect_error(env$iad_exact_remedy_plan(paths[-1L], ledger_path, "relative"))
  expect_error(env$iad_exact_remedy_plan(paths, ledger_path, "total"))
  bad <- ledger; bad$model_index <- 45L; bad$initial_sha256 <- initial_hashes[[2L]]
  env$iad_csv(bad, ledger_path)
  expect_error(env$iad_exact_remedy_plan(paths, ledger_path, "relative"))
  env$iad_csv(dplyr::bind_rows(ledger, ledger), ledger_path)
  expect_error(env$iad_exact_remedy_plan(paths, ledger_path, "relative"))
  bad <- ledger; bad$initial_sha256 <- paste(rep("0", 64L), collapse = "")
  env$iad_csv(bad, ledger_path)
  expect_error(env$iad_exact_remedy_plan(paths, ledger_path, "relative"))
  env$iad_csv(ledger, ledger_path)
  checkpoints[[1L]]$diagnostics$geometry_pass <- FALSE
  env$iad_qsave(checkpoints[[1L]], paths[[1L]])
  expect_error(env$iad_exact_remedy_plan(paths, ledger_path, "relative"))
  sentinel <- env$iad_qsave(list(skip = TRUE), file.path(root, "none.qs"))
  no_cases <- env$iad_exact_remedy_plan(sentinel, ledger_path, "total")
  stopifnot(identical(env$iad_review_exact(no_cases[[1L]]), sentinel))
}
test_exact_remedies()
cat("Exact remedy selection, pass-through, immutable originals, resume identity and failed-retry gates passed. No sampling.\n")

if ("--review-current" %in% commandArgs(trailingOnly = TRUE)) {
  for (family in iad_families()) {
    paths <- targets::tar_read_raw(paste0("iad_", family, "_exact_file"), store = "_targets_age_duration")
    plans <- iad_exact_remedy_plan(paths, "sensitivity/age-duration-exact-remedies.csv", family)
    actual <- dplyr::bind_rows(lapply(plans, `[[`, "record"))
    expected <- iad_read_csv(iad_path(family, "loo", "exact_cases.csv"))
    keys <- c("family", "outcome", "model_index", "observation")
    stopifnot(nrow(actual) == nrow(expected), !anyDuplicated(actual[keys]),
      nrow(dplyr::anti_join(actual, expected, by = keys)) == 0L)
    retry <- vapply(plans, function(p) !is.null(p$remedy), logical(1))
    for (plan in plans[!retry]) {
      checkpoint <- qs::qread(plan$initial_path)
      correction <- checkpoint$correction
      stopifnot(isTRUE(checkpoint$diagnostics$diagnostic_pass),
        isTRUE(correction$exact_refit_summary$diagnostic_pass),
        identical(as.integer(correction$exact_refit_summary$observation), as.integer(plan$record$observation)),
        identical(as.integer(correction$model_index), as.integer(plan$record$model_index)),
        identical(correction$family, family), identical(correction$outcome, plan$record$outcome),
        isTRUE(is.finite(correction$exact_refit_summary$exact_elpd_loo)))
    }
    cat(sprintf("%s: %d initial checkpoints, %d preserved, %d reviewed retries.\n",
      family, length(plans), sum(!retry), sum(retry)))
  }
  Sys.setenv(PROTEIN_INTERACTION_PHASE = "exact_remedy", PROTEIN_MODEL_WORKERS = "4")
  targets::tar_validate(script = "_targets_age_duration.R", store = "_targets_age_duration")
  Sys.setenv(PROTEIN_INTERACTION_PHASE = "none")
  cat("Current checkpoint coverage and pipeline validation passed. No pipeline run.\n")
}
