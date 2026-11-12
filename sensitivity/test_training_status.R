# Non-sampling tests: actual frozen inputs, all candidate designs, gates and pairing.
stopifnot(getRversion() == package_version("4.4.3"))
suppressPackageStartupMessages({library(tidyverse); library(brms); library(metafor); library(mgcv); library(lspline)})
for (path in c("R/functions.R", "R/primary_model_diagnostics.R", "R/primary_model_loo.R",
  "R/primary_model_averaging.R", "R/incremental_update.R", "R/update_2026_absolute_exposures.R",
  "R/verification_artifacts.R", "R/age_duration_interactions.R", "R/training_status_workflow.R")) source(path)
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
Sys.setenv(PROTEIN_TRAINING_PHASE = "none")
for (phase in c("pilot", "fit", "remedy", "loo", "exact", "correction", "combine", "report")) expect_error(iad_require_phase(phase))
expect_error(tst_age_duration_gate(root = file.path(tempdir(), "absent-age-duration-fixture")))
cp <- iad_read_checkpoints(iad_checkpoint_inputs())
specs <- iad_specifications(cp); table <- dplyr::bind_rows(specs)
stopifnot(nrow(table) == 248L, length(iad_pilot_specs(specs)) == 24L,
  !any(grepl("resistance_exercise", table$formula)),
  all(grepl("training_status_code", table$formula)),
  !anyDuplicated(table[c("family", "outcome", "model_index")]))
pilots <- dplyr::bind_rows(iad_pilot_specs(specs))[, c("family", "outcome", "model_index")]
pilots$diagnostic_pass <- TRUE
for (family in iad_families()) iad_check_family_pilot(pilots, family, cp)
broken <- pilots; broken$diagnostic_pass[1] <- FALSE
expect_error(iad_check_family_pilot(broken, pilots$family[1], cp))
expect_error(iad_check_family_pilot(pilots[-1, ], pilots$family[1], cp))
expect_error(iad_check_family_pilot(dplyr::bind_rows(pilots, pilots[1, ]), pilots$family[1], cp))
for (family in iad_families()) for (outcome in iad_outcomes()) {
  data <- cp[[family]]$data[[outcome]]
  stopifnot(all(data$resistance_exercise_code == 0), setequal(data$training_status_code, 0:1))
  copy <- data; copy$training_status_code[1] <- 1 - copy$training_status_code[1]
  stopifnot(iad_fingerprint(copy) != iad_fingerprint(data))
  grid <- tst_grid(data, family)
  stopifnot(all(grid$newdata$resistance_exercise_code == 0),
    all(grid$newdata$training_status_code == as.integer(grid$keys$training_status == "untrained")),
    sum(abs(grid$keys$dose - grid$reference) < 1e-10) == 4L)
  unit <- if (family == "relative") 1 else unique(data$absolute_exposure_unit_g_day)
  stopifnot(max(abs(grid$newdata$protein_intake - (grid$keys$dose - grid$reference) / unit)) < 1e-12)
}
priors <- iad_priors("sensitivity/age-duration-base-priors.csv")
audit <- iad_preflight(specs, cp, priors)
stopifnot(nrow(audit) == 248L, all(audit$valid))
print(audit |> dplyr::count(family, outcome, valid))
# Coherent status pairing removes additive status differences from moderation.
set.seed(88); x <- matrix(rnorm(400), 100); shift <- rnorm(100)
delta <- function(z) sweep(z, 1, z[, 1], "-")
stopifnot(max(abs(delta(sweep(x, 1, shift, "+")) - delta(x))) < 1e-12)
cat("All 248 RT/status Stan inputs, prediction grids and default-off gates passed; no sampling.\n")
