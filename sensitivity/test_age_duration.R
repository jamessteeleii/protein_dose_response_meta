# Non-sampling structural and statistical-invariant checks. No tar_make calls.
stopifnot(getRversion() == package_version("4.4.3"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(brms)
  library(mgcv)
  library(lspline)
})
for (path in c("R/functions.R", "R/primary_model_diagnostics.R", "R/primary_model_loo.R",
  "R/primary_model_averaging.R", "R/incremental_update.R", "R/update_2026_absolute_exposures.R",
  "R/verification_artifacts.R", "R/age_duration_interactions.R")) source(path)
expect_error <- function(expression) {
  stopifnot(inherits(tryCatch({ force(expression); NULL }, error = identity), "error"))
}
Sys.setenv(PROTEIN_INTERACTION_PHASE = "none")
for (phase in c("pilot", "fit", "remedy", "loo", "exact", "combine", "report")) {
  stopifnot(!iad_phase(phase)); expect_error(iad_require_phase(phase))
}
expect_error(iad_fit(list(), NULL, NULL, NULL, "pilot"))
expect_error(iad_compute_loo(list()))
expect_error(iad_exact(list(skip = TRUE)))
expect_error(iad_combine(NULL, NULL, NULL))
expect_error(iad_report_scope(NULL, NULL, NULL))
Sys.setenv(PROTEIN_INTERACTION_PHASE = "invalid")
expect_error(iad_phase("pilot"))
Sys.setenv(PROTEIN_INTERACTION_PHASE = "none")

cat("Reading and checking frozen primary inputs and compact final LOO reductions.\n")
checkpoints <- iad_read_checkpoints(iad_checkpoint_inputs())
print(dplyr::bind_rows(lapply(checkpoints, `[[`, "counts")))
specifications <- iad_specifications(checkpoints)
table <- dplyr::bind_rows(specifications)
stopifnot(nrow(table) == 336L,
  all(table$iter_sampling == 6000L), all(table$iter_warmup == 2000L),
  all(table$chains == 4L), all(table$adapt_delta >= 0.999),
  !anyDuplicated(table[c("family", "outcome", "model_index")]),
  nrow(table[table$family == "relative", ]) == 144L,
  nrow(table[table$family == "total", ]) == 96L,
  nrow(table[table$family == "supplement", ]) == 96L,
  length(iad_pilot_specs(specifications)) == 24L,
  all(grepl("age_outcome_z", table$formula)), all(grepl("duration_outcome_z", table$formula)),
  !any(table$parent_index %in% c(1L, 5L)))
rt <- table$parent_index %in% c(18L, 19L, 20L, 24L, 25L, 26L)
stopifnot(all(grepl("resistance_exercise_code", table$formula[rt])))
for (s in specifications) {
  parent <- checkpoints[[s$family]]$specifications
  parent <- parent[parent$outcome == s$outcome & parent$model_index == s$parent_index, ]
  stopifnot(startsWith(s$formula, replace_primary_standardised_moderators(parent$formula)))
}

profiles <- iad_profiles(checkpoints)
stopifnot(nrow(profiles) == 72L, all(is.finite(profiles$age)), all(is.finite(profiles$weeks)))
for (family in iad_families()) for (outcome in iad_outcomes()) {
  data <- checkpoints[[family]]$data[[outcome]]
  selected <- profiles[profiles$family == family & profiles$outcome == outcome, ]
  reference <- iad_prediction_grid(data, family, selected[selected$profile == "reference", ])
  typical <- iad_prediction_grid(data, family, selected[selected$profile == "typical", ])
  stopifnot(all(reference$newdata$age_centre == 0), all(reference$newdata$duration_centre == 0),
    identical(reference$keys, typical$keys),
    all(reference$newdata$age_outcome_z_rt_no == reference$newdata$age_outcome_z *
      (reference$newdata$resistance_exercise_code == 1)),
    sum(abs(reference$keys$protein_intake - reference$reference) < 1e-10) == 2L)
  copy <- data; copy$age_centre[[1]] <- copy$age_centre[[1]] + 1
  stopifnot(iad_fingerprint(data) != iad_fingerprint(copy))
  copy <- data; copy$duration_centre[[1]] <- copy$duration_centre[[1]] + 1
  stopifnot(iad_fingerprint(data) != iad_fingerprint(copy))
}

# Posterior pairing: additive profile shifts cancel within-dose contrasts,
# whereas an interaction changes the contrast and survives subtraction.
set.seed(1988)
base <- matrix(rnorm(400), nrow = 100)
shift <- rnorm(100)
low <- base; high <- sweep(base, 1L, shift, "+")
increment <- function(x) sweep(x, 1L, x[, 1], "-")
stopifnot(max(abs(increment(high) - increment(low))) < 1e-12)
high <- high + outer(rep(0.2, 100), 0:3)
stopifnot(max(abs((increment(high) - increment(low))[, 4] - 0.6)) < 1e-12)

# Exact-LOO gate tests include reordered observations, missing corrections,
# uncorrected warnings and failed refits. All LOO values are synthetic.
manifest <- tibble::tibble(study_no = as.character(1:3), arm_no = 1:3,
  outcome_no = 1:3, time1 = 0, time2 = 12, yi = 0, vi = 0.1)
set.seed(21)
synthetic <- suppressWarnings(loo::loo(matrix(rnorm(1200, -1, 0.1), ncol = 3)))
synthetic$diagnostics$pareto_k[] <- 0.1
result <- list(loo = synthetic, observation_manifest = manifest)
iad_validate_loo(result, manifest)
expect_error(iad_validate_loo(result, manifest[3:1, ]))
bad <- result; bad$loo$diagnostics$pareto_k[[1]] <- 0.9
expect_error(iad_validate_loo(bad, manifest))
bad <- result; bad$loo$diagnostics$pareto_k[[1]] <- NA_real_
expect_error(iad_validate_loo(bad, manifest))
bad$exact_refit_summary <- tibble::tibble(observation = 1L, diagnostic_pass = FALSE, exact_elpd_loo = -1)
expect_error(iad_validate_loo(bad, manifest))
bad$exact_refit_summary$diagnostic_pass <- TRUE
iad_validate_loo(bad, manifest)

cat("Building all 336 Stan inputs, without compilation or sampling.\n")
priors <- iad_priors(here::here("sensitivity", "age-duration-base-priors.csv"))
preflight <- iad_preflight(specifications, checkpoints, priors)
stopifnot(nrow(preflight) == 336L, all(preflight$valid))
print(preflight |> dplyr::count(family, interaction))
cat("All age/duration non-sampling checks passed. No pipeline or model fits were run.\n")
