# Static reporting interface ----

# Locate the pipeline-produced tables and figures that the preprint and
# supplements consume.

# Portable inputs for rendering the preprint from frozen analytical results.
preprint_report_root <- function(root = ".") {
  file.path(root, "manuscript", "generated", "reporting")
}

preprint_report_files <- function(root = ".", section = "") {
  paths <- list.files(file.path(preprint_report_root(root), section),
                      recursive = TRUE, full.names = TRUE)
  if (!length(paths)) stop("Missing preprint report inputs: ", section)
  sort(paths)
}

# Write the portable reporting inputs from the verified analysis artifacts so
# the documents can render from static files.
read_preprint_arm_intake_summary <- function(path) {
  arms <- readr::read_csv(path, col_select = c(study_no, arm_no, protein_intake_raw),
                         show_col_types = FALSE) |>
    dplyr::distinct()
  intake <- arms$protein_intake_raw
  stopifnot(length(intake) > 0L, all(is.finite(intake)),
            !anyDuplicated(arms[c("study_no", "arm_no")]))
  counts <- list(arms = nrow(arms), below_rda = sum(intake < 0.8), above_2_2 = sum(intake > 2.2))
  list(
    `arm-intake-tail-examples` = sprintf(
      "%.1f%% of primary-analysis arms were below 0.8 g/kg/day and %.1f%% were above 2.2 g/kg/day",
      100 * counts$below_rda / counts$arms, 100 * counts$above_2_2 / counts$arms
    ),
    `arm-intake-counts` = counts
  )
}

read_preprint_rt_contrasts <- function(root = ".") {
  bundle <- file.path(root, "results", "verification", "primary_update_2026")
  sesoi <- qs::qread(file.path(bundle, "reference_values.qs"))$rope_half_width
  support <- utils::read.csv(file.path(bundle, "observed_support.csv.gz"))
  doses <- c(0.8, 1.1, 1.6, 2.2)
  average_dir <- file.path(root, "model_archive", "update_2026", "primary",
                           "reporting", "model_average")

  # Saved draws retain the 0.8 g/kg/day physical-function prediction. Additive
  # age and duration terms cancel within each dose contrast.
  dplyr::bind_rows(lapply(c("whole_body", "hypertrophy", "strength", "function"), function(outcome) {
    average <- qs::qread(file.path(average_dir, paste0(outcome, ".qs")))
    observed <- support[support$outcome == outcome & support$resistance_exercise_code == 0, ]
    stopifnot(identical(average$outcome, outcome), nrow(observed) == 1L)
    columns <- vapply(doses, function(dose) {
      index <- which(average$grid$resistance_exercise_code == 0 &
                       abs(average$grid$protein_intake - dose) < 1e-10)
      if (length(index) != 1L) stop("Missing unique RT prediction for ", outcome, " at ", dose)
      index
    }, integer(1))
    lower_doses <- head(doses, -1L)
    higher_doses <- tail(doses, -1L)
    draw_differences <- average$draws[, tail(columns, -1L)] - average$draws[, head(columns, -1L)]
    intervals <- apply(draw_differences, 2L, stats::quantile, probs = c(0.025, 0.975))
    data.frame(
      outcome = outcome, resistance_exercise_code = 0L,
      protein_intake_reference = lower_doses, protein_intake_proposed = higher_doses,
      estimate = colMeans(draw_differences), lower_95 = intervals[1L, ], upper_95 = intervals[2L, ],
      probability_benefit = colMeans(draw_differences > sesoi),
      probability_equivalence = colMeans(abs(draw_differences) <= sesoi),
      probability_detriment = colMeans(draw_differences < -sesoi),
      supported = lower_doses >= observed$observed_min & higher_doses <= observed$observed_max
    )
  })) |>
    dplyr::arrange(protein_intake_reference, outcome)
}

read_preprint_trial_design <- function(representative_contrasts, root = ".", preprocessing_path = file.path(
    root, "results", "update_2026", "primary_refit", "preprocessing.csv")) {
  source(file.path(root, "R", "training_status_power_audit.R"), local = TRUE)
  correlation <- utils::read.csv(preprocessing_path)$pre_post_r_imputation
  if (length(correlation) != 1L || !is.finite(correlation) || abs(correlation) >= 1) {
    stop("Expected the pooled pre-post correlation from the final primary preprocessing.",
         call. = FALSE)
  }
  analysis_data <- utils::read.csv(file.path(
    root, "results", "verification", "primary_update_2026", "analysis_data.csv.gz"
  ))
  dose_pairs <- unique(representative_contrasts[c("protein_intake_reference", "protein_intake_proposed")])
  design <- protein_contrast_trial_design(representative_contrasts, analysis_data,
                                         pre_post_r = correlation, dose_pairs = unname(as.matrix(dose_pairs)))
  transform(design$contrasts, pre_post_r = design$pre_post_r, alpha = design$alpha,
            target_power = design$target_power, trained_min = design$trained_min)
}

write_preprint_report_inputs <- function(root = ".") {
  destination <- preprint_report_root(root)
  verification <- file.path(root, "results", "verification")
  update <- file.path(root, "results", "update_2026")
  primary <- "primary_update_2026"
  measurement <- c("measurement_prereg_methods", "measurement_direct_only")
  relative <- c("primary", "mundlak", "within", "achieved", "assignment",
                "eiv", "rescor", "strategy", "strategy_within")
  supplementation <- c("primary", "mundlak", "within", "rescor",
                        "strategy", "strategy_within")
  bundle_names <- paste0(
    c(relative, paste0("total_", relative),
      paste0("supplement_", supplementation), measurement),
    "_update_2026"
  )
  written <- character()

  # Route each derivative through a shared writer that collects its output
  # path.
  write_output <- function(source, relative_path, write) {
    if (!file.exists(source)) stop("Missing frozen reporting input: ", source)
    path <- file.path(destination, relative_path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    extension <- if (endsWith(path, ".gz")) ".gz" else ""
    temporary <- tempfile(tmpdir = dirname(path), fileext = extension)
    on.exit(unlink(temporary))
    write(source, temporary)
    fs::file_move(temporary, path)
    written <<- c(written, path)
    invisible(path)
  }

  copy_output <- function(source, relative_path) {
    write_output(source, relative_path, function(from, to) {
      if (!file.copy(from, to)) stop("Could not copy ", from)
    })
  }
  write_csv <- function(source, relative_path, transform, columns = NULL) {
    write_output(source, relative_path, function(from, to) {
      data <- if (is.null(columns)) {
        readr::read_csv(from, show_col_types = FALSE)
      } else {
        readr::read_csv(from, show_col_types = FALSE, col_select = dplyr::any_of(columns))
      }
      readr::write_csv(transform(data), to)
    })
  }

  copy_files <- function(source, section, names) {
    for (name in names) copy_output(file.path(source, name), file.path(section, name))
  }

  diagnostic_candidates <- c("final_sampling_diagnostic_summary.csv.gz",
    "sampling_diagnostic_summary.csv.gz", "diagnostic_summary.csv.gz")
  exact_candidates <- c("exact_refit_diagnostic_summary.csv.gz",
    "exact_loo_diagnostic_summary.csv.gz", "final_exact_refit_diagnostic_summary.csv.gz",
    "loo_exact_refit_summary.csv.gz", "accepted_exact_refit_diagnostic_summary.csv.gz")
  for (bundle in bundle_names) {
    source <- file.path(verification, bundle)
    section <- file.path("bundles", bundle)
    columns <- c("study_no", "arm_no", "protein_arm_no", "comparator_arm_no",
                 "analysis_outcome", "condition", "n", "time1", "time2")
    write_csv(file.path(source, "analysis_data.csv.gz"),
              file.path(section, "sample_rows.csv.gz"),
              identity, columns)
    diagnostic <- diagnostic_candidates[file.exists(file.path(source, diagnostic_candidates))]
    exact <- exact_candidates[file.exists(file.path(source, exact_candidates))]
    if (!length(diagnostic)) stop("Missing sampling diagnostics for ", bundle)
    if (bundle == primary) {
      copy_files(source, section, c("diagnostic_summary.csv.gz",
        "diagnostic_classification.csv.gz", "posterior_predictive_outcome_summary.csv.gz",
        "loo_outcome_summary.csv.gz", "loo_model_summary.csv.gz",
        "loo_method_comparison.csv.gz", "loo_observation_diagnostics.csv.gz",
        "loo_exact_refit_summary.csv.gz", "stacking_weights.csv.gz"))
    } else {
      diagnostic <- diagnostic[[1L]]
      for (name in c(diagnostic, head(exact, 1L))) {
        fields <- c("outcome", "accepted_diagnostic_pass", "diagnostic_pass")
        write_csv(file.path(source, name), file.path(section, name), function(data) {
          if (!any(fields[-1L] %in% names(data))) {
            # The source exact-refit summary records accepted refits only.
            data$diagnostic_pass <- rep(TRUE, nrow(data))
          }
          dplyr::select(data, dplyr::any_of(fields))
        })
      }
    }
    copy_files(source, section, c("dose_response_predictions.csv.gz", "local_equivalence.csv.gz"))
    if (bundle %in% paste0(measurement, "_update_2026")) {
      copy_files(source, section, "selection_summary.csv.gz")
      write_csv(file.path(source, "selection_audit.csv.gz"),
                file.path(section, "selection_audit.csv.gz"),
                function(data) dplyr::select(data, -dplyr::any_of(c(
                  "measurement_annotation_source", "measurement_confirmation_date", "source_row"
                ))))
    }
    if (bundle %in% paste0(c("total_primary", "supplement_primary", measurement), "_update_2026")) {
      write_csv(file.path(source, "pairwise_contrasts.csv.gz"),
                file.path(section, "supported_meaningful_contrasts.csv.gz"),
                function(data) dplyr::filter(data,
                  supported %in% TRUE, protein_intake_proposed > protein_intake_current,
                  probability_meaningful_benefit >= 0.95 | probability_meaningful_detriment >= 0.95
                ))
    }
  }

  for (exposure in c("relative", "total", "supplement")) {
    for (outcome in c("whole_body", "hypertrophy", "strength", "function")) {
      bundle <- paste0("age_duration_", exposure, "_", outcome, "_update_2026")
      copy_files(file.path(verification, bundle), file.path("bundles", bundle),
                 "stacking_weights.csv.gz")
    }
  }

  write_output(file.path(verification, primary, "reference_values.qs"),
               "reference-values.yml", function(from, to) {
    values <- lapply(qs::qread(from), function(value) {
      if (inherits(value, "Date")) format(value, "%Y-%m-%d") else value
    })
    writeLines(c("# Generated from frozen primary reference values. Do not edit.",
                 yaml::as.yaml(values, precision = 17)), to, useBytes = TRUE)
  })

  fields <- c("study_no", "doi", "pmid", "arm_no", "outcome_no", "analysis_outcome", "outcome_group",
              "measurement_method", "unit_of_measurement", "condition", "n", "time1", "time2",
              "baseline_mean", "baseline_sd", "baseline_n", "followup_mean", "followup_sd",
              "followup_n", "yi", "vi", "protein_intake", "resistance_exercise_code")
  write_csv(file.path(verification, primary, "analysis_data.csv.gz"),
            "primary-study-effects.csv.gz",
            function(data) dplyr::select(data, dplyr::all_of(fields)), fields)
  copy_files(file.path(update, "sensitivity", "comparisons"), "relative-comparison",
             "reference_increment_overlay_data.csv.gz")

  # The reviewed RoB 2 publication files are already bundled under rob2/.
  # GRADE consumes the study-number keyed judgments in data/GRADE/.
  grade <- file.path(update, "grade")
  copy_files(grade, "grade", c(
    "evidence_profile.csv", "population_effects.csv", "inconsistency_summary.csv",
    "inconsistency_study_summary.csv.gz", "inconsistency_leave_one_out.csv.gz",
    "empirical_support.csv", "rob2_summary.csv", "publication_bias_summary.csv"
  ))
  weights <- "inconsistency_restricted_stacking_weights.csv"
  write_csv(file.path(grade, weights), file.path("grade", weights),
            function(data) dplyr::select(data, -final_loo_path))
  write_csv(file.path(grade, "preflight.csv"), "grade/preflight.csv",
            function(data) dplyr::mutate(data, detail = ifelse(
              grepl("^[A-Za-z]:[/\\\\]|^/", detail), basename(detail), detail
            )))

  bias <- file.path(update, "publication_bias")
  copy_files(bias, "publication-bias", "high_low_outcome_summary.csv")
  write_csv(file.path(bias, "robma", "fit_manifest.csv"),
            "publication-bias/fit_manifest.csv",
            function(data) dplyr::select(data, -fit_path))
  analyses <- c("primary_whole_body", "primary_direct_muscle", "primary_strength",
                "primary_function", "secondary_omnibus")
  reports <- as.vector(outer(analyses, c("-robma_summary.txt", "-robma_models.txt",
                                        "-robma_zplot.txt", "-brma_summary.txt"), paste0))
  copy_files(file.path(bias, "robma", "reporting"), "publication-bias", reports)
  copy_files(file.path(root, "plots", "update_2026", "publication_bias", "robma"),
             "publication-bias", paste0(analyses, "-zplot.png"))

  write_output(file.path(verification, primary, "analysis_data.csv.gz"),
               "arm-intake-summary.yml", function(from, to) {
    writeLines(c("# Generated from frozen primary study arms. Do not edit.",
                 yaml::as.yaml(read_preprint_arm_intake_summary(from))), to, useBytes = TRUE)
  })
  copy_files(file.path(update, "sensitivity", "comparisons"), "relative-comparison",
             "reference_increment_overlay_data.csv.gz")

  representative_contrasts <- read_preprint_rt_contrasts(root)
  trial_design <- read_preprint_trial_design(representative_contrasts, root)
  write_output(file.path(verification, primary, "observed_support.csv.gz"),
               "representative-rt-contrasts.csv", function(from, to) {
    readr::write_csv(representative_contrasts, to)
  })
  write_output(file.path(update, "primary_refit", "preprocessing.csv"),
               "trial-design.csv", function(from, to) {
    readr::write_csv(trial_design, to)
  })
  copy_files(grade, "grade", c(
    "evidence_profile.csv", "population_effects.csv", "inconsistency_summary.csv",
    "inconsistency_study_summary.csv.gz", "inconsistency_leave_one_out.csv.gz",
    "empirical_support.csv", "rob2_summary.csv", "publication_bias_summary.csv"
  ))
  weights <- "inconsistency_restricted_stacking_weights.csv"
  write_csv(file.path(grade, weights), file.path("grade", weights),
            function(data) dplyr::select(data, -final_loo_path))
  write_csv(file.path(grade, "preflight.csv"), "grade/preflight.csv",
            function(data) dplyr::mutate(data, detail = ifelse(
              grepl("^[A-Za-z]:[/\\\\]|^/", detail), basename(detail), detail
            )))

  bias <- file.path(update, "publication_bias")
  copy_files(bias, "publication-bias", "high_low_outcome_summary.csv")
  write_csv(file.path(bias, "robma", "fit_manifest.csv"),
            "publication-bias/fit_manifest.csv",
            function(data) dplyr::select(data, -fit_path))
  analyses <- c("primary_whole_body", "primary_direct_muscle", "primary_strength",
                "primary_function", "secondary_omnibus")
  reports <- as.vector(outer(analyses, c("-robma_summary.txt", "-robma_models.txt",
                                        "-robma_zplot.txt", "-brma_summary.txt"), paste0))
  copy_files(file.path(bias, "robma", "reporting"), "publication-bias", reports)
  copy_files(file.path(root, "plots", "update_2026", "publication_bias", "robma"),
             "publication-bias", paste0(analyses, "-zplot.png"))

  sort(written)
}
