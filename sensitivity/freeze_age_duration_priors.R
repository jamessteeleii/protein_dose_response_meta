# Freeze only the already accepted population coefficient prior recipe. This
# does not fit a model or rerun the historical prior-calibration meta-analysis.
stopifnot(getRversion() == package_version("4.4.3"))
manifest <- read.csv("results/update_2026/primary_refit/accepted_fit_manifest.csv")
rows <- provenance <- list()
for (outcome in c("whole_body", "hypertrophy", "strength", "function")) {
  path <- manifest$fit_path[manifest$outcome == outcome & manifest$model_index == 1L]
  stopifnot(length(path) == 1L, file.exists(path))
  candidate <- qs::qread(path)
  stopifnot(candidate$outcome == outcome, candidate$model_index == 1L)
  prior <- as.data.frame(candidate$model$prior)
  prior <- prior[prior$class == "b" & prior$coef %in% c("", "Intercept", "resistance_exercise_code"), ]
  stopifnot(nrow(prior) == 3L, all(nzchar(prior$prior)))
  rows[[outcome]] <- data.frame(outcome = outcome, prior = prior$prior,
    class = prior$class, coef = prior$coef, stringsAsFactors = FALSE)
  provenance[[outcome]] <- data.frame(outcome = outcome,
    source_fit = sub(paste0(normalizePath(getwd(), winslash = "/"), "/"), "", path, fixed = TRUE),
    source_fit_sha256 = digest::digest(file = path, algo = "sha256"),
    frozen_on = "2026-11-02", stringsAsFactors = FALSE)
  rm(candidate); gc(verbose = FALSE)
}
write.csv(do.call(rbind, rows), "sensitivity/age-duration-base-priors.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, provenance), "sensitivity/age-duration-prior-provenance.csv", row.names = FALSE, na = "")
cat("Frozen 12 original coefficient-prior records and four source hashes. No fitting performed.\n")
