# Power for a two-sample test of change scores. Delta is in baseline SD units.
protein_change_score_power <- function(delta, pre_post_r, n = NULL,
                                       power = NULL, alpha = 0.05) {
  stats::power.t.test(
    n = n, delta = abs(delta), sd = sqrt(2 * (1 - pre_post_r)),
    sig.level = alpha, power = power, type = "two.sample",
    alternative = "two.sided", strict = TRUE, tol = 1e-9
  )
}

# Descriptive audit of saved reporting inputs; does not execute targets or fit models.
training_status_power_audit <- function(
  bundle = "results/verification/primary_update_2026",
  workbooks = c(
    "data/protein_meta_data_all.xlsx",
    "data/update_2026/protein_meta_data_update_2026.xlsx"
  )
) {
  manifest <- utils::read.csv(file.path(bundle, "manifest.csv"))
  read_verified <- function(filename) {
    path <- file.path(bundle, filename)
    expected <- manifest$sha256[manifest$file == filename]
    if (length(expected) != 1L ||
        !identical(digest::digest(file = path, algo = "sha256"), expected)) {
      stop("Missing or mismatched verification checksum: ", path)
    }
    utils::read.csv(gzfile(path), stringsAsFactors = FALSE)
  }
  analysis <- read_verified("analysis_data.csv.gz")
  predictions <- read_verified("dose_response_predictions.csv.gz")
  source_manifest <- read_verified("source_manifest.csv.gz")
  source_checks <- do.call(rbind, lapply(workbooks, function(path) {
    expected <- source_manifest$sha256[
      basename(source_manifest$path) == basename(path)
    ]
    actual <- digest::digest(file = path, algo = "sha256")
    if (length(expected) != 1L || !identical(actual, expected)) {
      stop("Workbook differs from saved primary-analysis input: ", path)
    }
    data.frame(path = path, sha256 = actual)
  }))
  normalise_status <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x[is.na(x) | x == ""] <- "not reported"
    x
  }
  combine_status <- function(x) paste(sort(unique(x)), collapse = ";")
  analysis$status <- normalise_status(analysis$trained_or_untrained)
  arm_groups <- split(analysis, interaction(
    analysis$study_no, analysis$arm_no, drop = TRUE
  ))
  arms <- do.call(rbind, lapply(arm_groups, function(x) {
    stopifnot(length(unique(x$status)) == 1L,
              length(unique(x$resistance_exercise_yes_no)) == 1L)
    data.frame(
      study_no = x$study_no[1], arm_no = x$arm_no[1],
      status = x$status[1], participants = max(x$n),
      rt = x$resistance_exercise_yes_no[1],
      longest_followup_weeks = max(x$time2),
      protein_intake_g_kg_day = x$protein_intake_raw[1]
    )
  }))
  rownames(arms) <- NULL
  studies <- do.call(rbind, lapply(split(arms, arms$study_no), function(x) {
    data.frame(
      study_no = x$study_no[1], status = combine_status(x$status),
      arms = nrow(x), participants = sum(x$participants),
      longest_followup_weeks = max(x$longest_followup_weeks)
    )
  }))
  raw_arms <- do.call(rbind, lapply(workbooks, function(path) {
    x <- readxl::read_excel(path, sheet = "Arm characteristics")
    data.frame(study_no = x[["Study No"]],
               arm_no = x[["Arm No"]],
               status = normalise_status(x[["Trained or Untrained?"]]))
  }))
  raw_arms <- raw_arms[!is.na(raw_arms$study_no) & !is.na(raw_arms$arm_no), ]
  stopifnot(!anyDuplicated(raw_arms[c("study_no", "arm_no")]))
  review_studies <- stats::aggregate(
    status ~ study_no, raw_arms, combine_status
  )
  summarise_status <- function(x, unit, weights = rep(1, nrow(x))) {
    counts <- stats::aggregate(weights, list(status = x$status), sum)
    names(counts)[2] <- "count"
    counts$percent <- 100 * counts$count / sum(counts$count)
    counts$unit <- unit
    counts[c("unit", "status", "count", "percent")]
  }
  status_summary <- rbind(
    summarise_status(review_studies, "All extracted study records"),
    summarise_status(studies, "Primary analysis studies"),
    summarise_status(arms, "Primary analysis arms"),
    summarise_status(arms, "Primary analysis participants", arms$participants)
  )
  rt_arms <- arms[arms$rt == "yes", ]
  rt_studies <- rt_arms[!duplicated(rt_arms$study_no), ]
  status_summary <- rbind(status_summary,
    summarise_status(rt_studies, "Studies contributing an RT arm"))
  outcome_studies <- unique(analysis[c("analysis_outcome", "study_no", "status")])
  outcome_summary <- do.call(rbind, lapply(
    split(outcome_studies, outcome_studies$analysis_outcome),
    function(x) summarise_status(x, x$analysis_outcome[1])
  ))
  design_summary <- do.call(rbind, lapply(sort(unique(studies$status)), function(s) {
    a <- arms[arms$status == s, ]
    w <- studies$longest_followup_weeks[studies$status == s]
    data.frame(status = s, studies = length(w),
      median_n_per_arm = stats::median(a$participants),
      median_weeks = stats::median(w),
      q1_weeks = unname(stats::quantile(w, 0.25)),
      q3_weeks = unname(stats::quantile(w, 0.75)),
      maximum_weeks = max(w), studies_above_12_weeks = sum(w > 12),
      studies_at_least_24_weeks = sum(w >= 24),
      studies_at_least_52_weeks = sum(w >= 52))
  }))
  # These are true hypothetical randomised between-arm differences in change,
  # expressed in baseline SD units, not within-arm changes or observed power.
  # Equal baseline/follow-up variances: SD(change)/SD(baseline) = sqrt(2*(1-r)).
  power <- expand.grid(
    n_per_arm = sort(unique(c(stats::median(arms$participants),
      stats::median(arms$participants[arms$status == "trained"])))),
    pre_post_r = c(0.5, 0.8, 0.9),
    target_difference_baseline_sd = c(0.05, 0.1, 0.2)
  )
  power$alpha <- 0.05
  power$change_sd_over_baseline_sd <- sqrt(2 * (1 - power$pre_post_r))
  power$power <- vapply(seq_len(nrow(power)), function(i) {
    protein_change_score_power(
      delta = power$target_difference_baseline_sd[i],
      pre_post_r = power$pre_post_r[i], n = power$n_per_arm[i],
      alpha = power$alpha[i]
    )$power
  }, numeric(1))
  # Difference of posterior means is the posterior mean contrast. Pointwise
  # intervals must not be subtracted to manufacture a contrast interval.
  mean_contrasts <- do.call(rbind, lapply(
    c("whole_body", "hypertrophy", "strength"), function(outcome) {
      do.call(rbind, lapply(list(c(1.2, 1.6), c(1.6, 2.2)), function(intakes) {
        p <- predictions[predictions$outcome == outcome &
                           predictions$resistance_exercise_code == 0, ]
        estimates <- vapply(intakes, function(intake) {
          value <- p$estimate[abs(p$protein_intake - intake) < 1e-8]
          stopifnot(length(value) == 1L)
          value
        }, numeric(1))
        data.frame(outcome = outcome, current_g_kg_day = intakes[1],
          proposed_g_kg_day = intakes[2],
          posterior_mean_difference_baseline_sd = diff(estimates))
      }))
    }
  ))
  list(status_summary = status_summary, outcome_summary = outcome_summary,
       design_summary = design_summary, arms = arms, studies = studies,
       power = power, illustrative_mean_contrasts = mean_contrasts,
       source_checks = source_checks)
}

write_training_status_power_audit <- function(
  audit, output_dir = "outputs/training_status_power"
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- vapply(names(audit), function(name) {
    path <- file.path(output_dir, paste0(name, ".csv"))
    utils::write.csv(audit[[name]], path, row.names = FALSE, na = "")
    path
  }, character(1))
  md_table <- function(x) {
    cells <- lapply(x, function(v) {
      if (is.numeric(v)) format(round(v, 4), trim = TRUE, scientific = FALSE)
      else as.character(v)
    })
    rows <- do.call(paste, c(cells, list(sep = " | ")))
    c(paste(names(x), collapse = " | "),
      paste(rep("---", ncol(x)), collapse = " | "), rows, "")
  }
  report <- c(
    "<!-- Generated by R/training_status_power_audit.R; do not hand-edit. -->",
    "# Training status and illustrative trial power", "",
    "Inputs: checksum-verified primary_update_2026 analysis_data and dose_response_predictions; both extraction workbooks (arm characteristics).",
    "The extraction set includes review-only/excluded records; primary counts describe only studies contributing to the saved primary analysis.",
    "Active and unreported status remain separate from untrained. Each study/arm is counted once; participant counts use the maximum analysed N per arm across outcomes and follow-ups. These are analysed sample sizes, not recruitment targets. No mixed-status studies occur in the primary set.",
    "These results describe the saved 1,864-row primary checkpoint, before the pending outcome-measurement restriction; they do not summarise a future corrected analysis.", "",
    "## Training status", "", md_table(audit$status_summary),
    "## Studies by outcome", "", md_table(audit$outcome_summary),
    "## Sample sizes and longest eligible follow-up per study", "",
    md_table(audit$design_summary),
    sprintf("Overall median analysed arm N: %g (IQR %g to %g). Median RT-arm N: %g.",
      stats::median(audit$arms$participants),
      stats::quantile(audit$arms$participants, 0.25),
      stats::quantile(audit$arms$participants, 0.75),
      stats::median(audit$arms$participants[audit$arms$rt == "yes"])), "",
    "## Conditional power scenarios", "",
    "Equal allocation, independent arms, normally distributed individual change scores with equal variances, two-sided two-sample t-test of a zero between-arm difference, alpha = 0.05, no multiplicity correction or further attrition. Assumed equal baseline and follow-up SDs imply SD(change)/SD(baseline) = sqrt(2*(1-r)). The pre-post correlations are illustrative assumptions, not estimates from this audit. Power includes both rejection tails (stats::power.t.test, strict = TRUE). Different analysis methods (e.g., ANCOVA), variances, correlations or missingness change power.",
    "Target differences are hypothetical true between-arm differences in mean change in baseline SD units. The 0.10 scenario equals the preregistered threshold; 0.05 is a smaller target and 0.20 a larger comparator. These are conditional sensitivity calculations, not estimates of each trial's actual power, observed/post-hoc power, or power averaged over posterior uncertainty. They test a zero difference; they do not calculate power to establish equivalence or to exceed a 0.10-SD meaningful-effect threshold.", "",
    md_table(audit$power),
    "## Context from the primary RT curves", "",
    "Differences between posterior mean predictions at 12 weeks and age 25. These conditional associations are not identified causal effects or training-status-specific effects. These selected contrasts illustrate scale only, and do not represent a selection rule for effects across the full evidence base. Contrast uncertainty is not calculated here.", "",
    md_table(audit$illustrative_mean_contrasts),
    "## Reproduce", "",
    "From the project R 4.4.3 session with renv restored:", "",
    "```r", "source('R/training_status_power_audit.R')",
    "audit <- training_status_power_audit()",
    "write_training_status_power_audit(audit)", "```", "",
    "No targets pipeline, model fitting, or manuscript rendering is invoked.", "",
    "Method: [R stats power.t.test documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/power.t.test.html).",
    paste("R version used:", getRversion())
  )
  path <- file.path(output_dir, "report.md")
  writeLines(report, path, useBytes = TRUE)
  invisible(c(files, report = path))
}

protein_contrast_trial_design <- function(population_effects, analysis_data,
                                          pre_post_r = 0.8, alpha = 0.05,
                                          target_power = 0.8) {
  outcomes <- c(whole_body = "whole-body lean/fat-free mass",
                hypertrophy = "measures of muscle size or mass",
                strength = "strength")
  dose_pairs <- rbind(c(0.8, 1.6), c(0.8, 2.4), c(1.6, 2.4))
  rt_arms <- subset(analysis_data, resistance_exercise_code == 0)
  trained <- subset(rt_arms, trained_or_untrained == "trained")
  median_arm_n <- function(data) {
    arm_sizes <- stats::aggregate(n ~ study_no + arm_no, data, max)$n
    ceiling(stats::median(arm_sizes))
  }
  n_rt <- median_arm_n(rt_arms)
  n_trained <- median_arm_n(trained)
  effects_rt <- subset(population_effects, resistance_exercise_code == 0)
  results <- lapply(seq_len(nrow(dose_pairs)), function(i) {
    doses <- dose_pairs[i, ]
    effects <- effects_rt[
      effects_rt$outcome %in% names(outcomes) &
        abs(effects_rt$protein_intake_reference - doses[1]) < 1e-8 &
        abs(effects_rt$protein_intake_proposed - doses[2]) < 1e-8, ]
    if (nrow(effects) != length(outcomes) || anyDuplicated(effects$outcome) ||
        !all(effects$supported %in% TRUE) || any(!is.finite(effects$estimate))) {
      stop("Expected one supported estimate per outcome for ",
           paste(doses, collapse = " vs "), " g/kg/day.")
    }
    chosen <- effects[which.min(abs(effects$estimate)), ]
    delta <- chosen$estimate
    power_test <- function(n = NULL, power = NULL) {
      protein_change_score_power(delta, pre_post_r, n, power, alpha)
    }
    required_n <- if (delta == 0) Inf else ceiling(power_test(power = target_power)$n)
    data.frame(
      lower_dose = doses[1], higher_dose = doses[2],
      outcome = unname(outcomes[chosen$outcome]), estimate = delta,
      n_trained = n_trained, power_trained = power_test(n = n_trained)$power,
      n_rt = n_rt, power_rt = power_test(n = n_rt)$power,
      required_n_per_arm = required_n
    )
  })
  list(contrasts = do.call(rbind, results), alpha = alpha,
       pre_post_r = pre_post_r, target_power = target_power,
       trained_min = min(trained$protein_intake_raw[
         trained$analysis_outcome %in% names(outcomes)
       ]))
}
