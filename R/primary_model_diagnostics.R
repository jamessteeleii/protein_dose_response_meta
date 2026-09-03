# Reproducible diagnostics for the final primary candidate-model ensemble.

infer_primary_candidate_source <- function(specification) {
  if (isTRUE(specification$third_remedy)) {
    return("third_remedied")
  }
  if (isTRUE(specification$second_remedy)) {
    return("second_remedied")
  }
  if (isTRUE(specification$requires_refit)) {
    return("remedied")
  }

  "original"
}

diagnose_primary_candidate_set <- function(
    fitted_candidates,
    ppc_draws = 400L
) {
  if (
    length(fitted_candidates) != 26L ||
      !identical(
        sort(unname(purrr::map_int(fitted_candidates, "model_index"))),
        1:26
      )
  ) {
    stop("Expected a complete 26-model candidate set.", call. = FALSE)
  }

  purrr::map(
    fitted_candidates,
    get_primary_candidate_model_diagnostics,
    ppc_draws = ppc_draws
  )
}

primary_sampling_issue <- function(row) {
  reasons <- character()

  if (row[["divergent_transitions"]] > 0L) {
    reasons <- c(
      reasons,
      sprintf(
        "%d divergent transition(s)",
        row[["divergent_transitions"]]
      )
    )
  }
  if (row[["maximum_treedepth_transitions"]] > 0L) {
    reasons <- c(
      reasons,
      sprintf(
        "%d maximum-treedepth transition(s)",
        row[["maximum_treedepth_transitions"]]
      )
    )
  }
  if (
    !is.na(row[["minimum_ebfmi"]]) &&
      row[["minimum_ebfmi"]] < 0.30
  ) {
    reasons <- c(
      reasons,
      sprintf("minimum E-BFMI %.3f", row[["minimum_ebfmi"]])
    )
  }
  if (
    !is.na(row[["maximum_rhat"]]) &&
      row[["maximum_rhat"]] > 1.01
  ) {
    reasons <- c(
      reasons,
      sprintf("maximum R-hat %.3f", row[["maximum_rhat"]])
    )
  }
  if (
    !is.na(row[["minimum_bulk_ess"]]) &&
      row[["minimum_bulk_ess"]] < 400
  ) {
    reasons <- c(
      reasons,
      sprintf("minimum bulk ESS %.0f", row[["minimum_bulk_ess"]])
    )
  }
  if (
    !is.na(row[["minimum_tail_ess"]]) &&
      row[["minimum_tail_ess"]] < 400
  ) {
    reasons <- c(
      reasons,
      sprintf("minimum tail ESS %.0f", row[["minimum_tail_ess"]])
    )
  }
  if (
    !is.na(row[["maximum_mcse_sd_ratio"]]) &&
      row[["maximum_mcse_sd_ratio"]] > 0.10
  ) {
    reasons <- c(
      reasons,
      sprintf(
        "maximum MCSE/SD ratio %.3f",
        row[["maximum_mcse_sd_ratio"]]
      )
    )
  }

  if (length(reasons)) {
    paste(reasons, collapse = "; ")
  } else {
    "None"
  }
}

primary_ppc_warning <- function(summary) {
  summary$maximum_absolute_standardised_residual > 3 |
    summary[["pit_below_0.025"]] > 0L |
    summary[["pit_above_0.975"]] > 0L |
    summary$predictive_90_interval_coverage < 0.80 |
    summary$predictive_90_interval_coverage > 0.98
}

summarise_primary_ppc_outcome <- function(diagnostics, outcome) {
  outcome_diagnostics <- diagnostics[
    vapply(
      diagnostics,
      function(x) identical(x$summary$outcome[[1L]], outcome),
      logical(1L)
    )
  ]
  outcome_diagnostics <- outcome_diagnostics[
    order(vapply(
      outcome_diagnostics,
      function(x) x$summary$model_index[[1L]],
      integer(1L)
    ))
  ]

  if (length(outcome_diagnostics) != 26L) {
    stop("Expected 26 diagnostic results for: ", outcome, call. = FALSE)
  }

  residual_tables <- purrr::map(
    outcome_diagnostics,
    "standardised_residuals"
  )
  observation_counts <- vapply(residual_tables, nrow, integer(1L))

  if (length(unique(observation_counts)) != 1L) {
    stop(
      "Candidate models do not share an observation set for: ",
      outcome,
      call. = FALSE
    )
  }

  first <- residual_tables[[1L]]
  identifiers <- intersect(
    c("study_no", "arm_no", "outcome_no"),
    names(first)
  )
  standardised_residuals <- do.call(
    cbind,
    purrr::map(residual_tables, "standardised_residual")
  )
  pits <- do.call(cbind, purrr::map(residual_tables, "pit"))

  flags <- tibble::as_tibble(first[, c(
    "outcome", "observation", identifiers, "observed"
  )]) |>
    dplyr::mutate(
      models_abs_standardised_residual_gt_3 = rowSums(
        abs(standardised_residuals) > 3
      ),
      models_pit_below_0_025 = rowSums(pits < 0.025),
      models_pit_above_0_975 = rowSums(pits > 0.975),
      median_standardised_residual = apply(
        standardised_residuals,
        1L,
        stats::median
      ),
      maximum_absolute_standardised_residual = apply(
        abs(standardised_residuals),
        1L,
        max
      ),
      median_pit = apply(pits, 1L, stats::median),
      flagged_in_any_model =
        models_abs_standardised_residual_gt_3 > 0L |
        models_pit_below_0_025 > 0L |
        models_pit_above_0_975 > 0L,
      flagged_in_all_models =
        models_abs_standardised_residual_gt_3 == 26L |
        models_pit_below_0_025 == 26L |
        models_pit_above_0_975 == 26L
    )

  outcome_summary <- dplyr::bind_rows(
    purrr::map(outcome_diagnostics, "summary")
  )
  summary <- tibble::tibble(
    outcome = outcome,
    observations = nrow(first),
    models = 26L,
    observations_flagged_in_any_model = sum(flags$flagged_in_any_model),
    observations_flagged_in_all_models = sum(flags$flagged_in_all_models),
    coverage_min = min(
      outcome_summary$predictive_90_interval_coverage
    ),
    coverage_max = max(
      outcome_summary$predictive_90_interval_coverage
    ),
    maximum_residual_min = min(
      outcome_summary$maximum_absolute_standardised_residual
    ),
    maximum_residual_max = max(
      outcome_summary$maximum_absolute_standardised_residual
    ),
    pit_below_min = min(outcome_summary[["pit_below_0.025"]]),
    pit_below_max = max(outcome_summary[["pit_below_0.025"]]),
    pit_above_min = min(outcome_summary[["pit_above_0.975"]]),
    pit_above_max = max(outcome_summary[["pit_above_0.975"]])
  )

  list(
    flags = dplyr::filter(flags, flagged_in_any_model),
    summary = summary
  )
}

compile_primary_model_diagnostic_audit <- function(...) {
  diagnostic_sets <- list(...)

  if (
    length(diagnostic_sets) != 4L ||
      !all(vapply(diagnostic_sets, length, integer(1L)) == 26L)
  ) {
    stop(
      "Expected four complete 26-model diagnostic sets.",
      call. = FALSE
    )
  }

  diagnostics <- unlist(diagnostic_sets, recursive = FALSE)
  summary <- dplyr::bind_rows(purrr::map(diagnostics, "summary")) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    ) |>
    dplyr::mutate(
      ppc_warning = primary_ppc_warning(dplyr::pick(dplyr::everything())),
      status = dplyr::case_when(
        !diagnostic_pass ~ "REQUIRES REMEDIATION",
        ppc_warning ~ "WARNING",
        TRUE ~ "PASS"
      ),
      loo_status = "NOT ASSESSED"
    )

  if (
    nrow(summary) != 104L ||
      anyDuplicated(paste(summary$outcome, summary$model_index))
  ) {
    stop("The primary diagnostic summary is incomplete.", call. = FALSE)
  }

  parameter_failures <- purrr::map_dfr(
    diagnostics,
    function(x) {
      x$parameter_failures |>
        dplyr::mutate(
          outcome = x$summary$outcome[[1L]],
          model_index = x$summary$model_index[[1L]],
          .before = 1L
        )
    }
  )

  classification <- summary |>
    dplyr::rowwise() |>
    dplyr::mutate(
      final_status = dplyr::case_when(
        !diagnostic_pass ~ "REQUIRES REMEDIATION",
        ppc_warning ~ "PASS WITH PPC WARNING",
        TRUE ~ "PASS"
      ),
      sampling_issue = primary_sampling_issue(dplyr::pick(
        dplyr::everything()
      )),
      exclusion_decision = "RETAIN",
      recommended_action = if (
        diagnostic_pass
      ) {
        paste(
          "Eligible for final predictive comparison; retain outcome-wide",
          "posterior-predictive warnings for sensitivity analysis."
        )
      } else {
        paste(
          "Apply the smallest diagnostic-driven remedy, refit, and require",
          "all sampling and convergence rules to pass before stacking."
        )
      }
    ) |>
    dplyr::ungroup()

  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  ppc <- purrr::map(
    outcomes,
    ~ summarise_primary_ppc_outcome(diagnostics, .x)
  )

  list(
    diagnostic_summary = summary,
    definitive_classification = classification,
    parameter_failures = parameter_failures,
    ppc_observation_flags = dplyr::bind_rows(purrr::map(ppc, "flags")),
    ppc_outcome_summary = dplyr::bind_rows(purrr::map(ppc, "summary"))
  )
}

write_primary_model_diagnostic_audit <- function(
    audit,
    output_dir = file.path(
      "_targets", "user", "primary_model_diagnostic_audit"
    )
) {
  if (
    nrow(audit$diagnostic_summary) != 104L ||
      nrow(audit$definitive_classification) != 104L ||
      nrow(audit$ppc_outcome_summary) != 4L
  ) {
    stop("The primary model diagnostic audit is invalid.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    definitive_classification = file.path(
      output_dir,
      "definitive_classification.csv"
    ),
    parameter_failures = file.path(output_dir, "parameter_failures.csv"),
    ppc_observation_flags = file.path(
      output_dir,
      "ppc_observation_flags.csv"
    ),
    ppc_outcome_summary = file.path(
      output_dir,
      "ppc_outcome_summary.csv"
    )
  )

  purrr::walk2(
    audit[names(paths)],
    paths,
    ~ readr::write_csv(.x, .y, na = "")
  )

  unname(paths)
}
