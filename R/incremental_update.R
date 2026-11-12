# Incremental systematic-review update -------------------------------------

update_2026_outcome_names <- function() {
  c(
    whole_body = "whole body lean mass",
    hypertrophy = "direct muscle size/mass",
    strength = "strength",
    `function` = "function"
  )
}

harmonise_workbook_table_types <- function(x, y) {
  common <- intersect(names(x), names(y))

  for (column in common) {
    x_class <- class(x[[column]])[[1]]
    y_class <- class(y[[column]])[[1]]

    if (identical(x_class, y_class)) next

    if (
      is.character(x[[column]]) || is.character(y[[column]]) ||
        is.factor(x[[column]]) || is.factor(y[[column]])
    ) {
      x[[column]] <- as.character(x[[column]])
      y[[column]] <- as.character(y[[column]])
    } else if (is.numeric(x[[column]]) || is.numeric(y[[column]])) {
      x[[column]] <- as.numeric(x[[column]])
      y[[column]] <- as.numeric(y[[column]])
    }
  }

  list(x = x, y = y)
}

combine_protein_workbooks <- function(canonical_data, update_data) {
  required_sheets <- c(
    "study_characteristics",
    "arm_characteristics",
    "arm_based_data_by_timepoint",
    "contrasts_between_timepoints"
  )

  if (
    !setequal(names(canonical_data), required_sheets) ||
      !setequal(names(update_data), required_sheets)
  ) {
    stop(
      "Canonical and update workbooks must contain the same four worksheets.",
      call. = FALSE
    )
  }

  canonical_studies <- canonical_data$study_characteristics$study_no
  update_studies <- update_data$study_characteristics$study_no
  overlap <- intersect(canonical_studies, update_studies)

  if (length(overlap)) {
    stop(
      "Study numbers overlap across source workbooks: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }

  purrr::map(
    required_sheets,
    function(sheet) {
      canonical_columns <- names(canonical_data[[sheet]])
      update_columns <- names(update_data[[sheet]])

      if (!identical(canonical_columns, update_columns)) {
        stop(
          "Workbook schema mismatch on worksheet: ", sheet,
          call. = FALSE
        )
      }

      aligned <- harmonise_workbook_table_types(
        canonical_data[[sheet]],
        update_data[[sheet]]
      )
      dplyr::bind_rows(aligned$x, aligned$y)
    }
  ) |>
    rlang::set_names(required_sheets)
}

apply_combined_review_decisions <- function(prepared, decisions) {
  update_studies <- decisions$study_no

  prepared |>
    dplyr::left_join(
      decisions |>
        dplyr::select(
          study_no,
          quantitative_status,
          include_quantitative,
          review_decision_reason = reason,
          review_decision_source = decision_source
        ),
      by = "study_no"
    ) |>
    dplyr::mutate(
      update_2026_data_source = dplyr::if_else(
        study_no %in% update_studies,
        "updated search",
        "frozen canonical"
      ),
      quantitative_status = dplyr::coalesce(
        quantitative_status,
        "include_frozen_canonical"
      ),
      include_quantitative = dplyr::coalesce(include_quantitative, TRUE),
      review_decision_reason = dplyr::coalesce(
        review_decision_reason,
        "Included in the frozen primary analysis"
      ),
      review_decision_source = dplyr::coalesce(
        review_decision_source,
        "Frozen primary-analysis eligibility decision"
      ),
      update_2026_analysis_ready =
        analysis_eligible_primary & include_quantitative
    )
}

prepare_update_2026_systematic_review_audit_data <- function(reviewed_data) {
  required <- c(
    "include_quantitative", "analysis_eligible_primary",
    "analysis_exclusion_count", "analysis_exclusion_reasons"
  )
  missing <- setdiff(required, names(reviewed_data))
  if (length(missing)) {
    stop(
      "Combined review data are missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  reviewed_data |>
    dplyr::mutate(
      exclusion_author_decision = !.data$include_quantitative,
      analysis_exclusion_count =
        .data$analysis_exclusion_count +
        as.integer(.data$exclusion_author_decision),
      analysis_exclusion_reasons = dplyr::case_when(
        .data$exclusion_author_decision &
          is.na(.data$analysis_exclusion_reasons) ~
          "excluded by author-group decision",
        .data$exclusion_author_decision ~ paste(
          .data$analysis_exclusion_reasons,
          "excluded by author-group decision",
          sep = "; "
        ),
        TRUE ~ .data$analysis_exclusion_reasons
      ),
      analysis_eligible_primary =
        .data$analysis_eligible_primary &
        .data$include_quantitative
    )
}

make_combined_primary_data_by_outcome <- function(reviewed_data) {
  outcome_names <- update_2026_outcome_names()
  analysis_data <- reviewed_data |>
    dplyr::filter(update_2026_analysis_ready)

  result <- purrr::imap(
    outcome_names,
    function(outcome_group_value, outcome_name) {
      outcome_data <- analysis_data |>
        dplyr::filter(outcome_group == outcome_group_value)

      if (identical(outcome_name, "function")) {
        outcome_data <- apply_function_direction(outcome_data)
      }

      add_primary_display_fields(outcome_data) |>
        dplyr::mutate(
          analysis_outcome = dplyr::recode(
            outcome_name,
            hypertrophy = "direct_muscle",
            `function` = "physical_function",
            .default = outcome_name
          )
        )
    }
  )

  names(result) <- names(outcome_names)
  result
}

validate_combined_primary_preparation <- function(
  combined_workbook,
  reviewed_data,
  data_by_outcome,
  frozen_analysis,
  decisions
) {
  model_data <- dplyr::bind_rows(data_by_outcome)
  key_columns <- c(
    "study_no", "arm_no_extracted", "outcome_group",
    "measurement_method", "unit_of_measurement", "time1", "time2"
  )
  duplicate_keys <- reviewed_data |>
    dplyr::count(dplyr::across(dplyr::all_of(key_columns))) |>
    dplyr::filter(n > 1L)

  frozen_model_rows <- model_data |>
    dplyr::filter(update_2026_data_source == "frozen canonical")
  update_model_rows <- model_data |>
    dplyr::filter(update_2026_data_source == "updated search")
  centres <- unique(stats::na.omit(reviewed_data$protein_median))
  imputed_r <- reviewed_data |>
    dplyr::filter(!has_change_sd, !is.na(r)) |>
    dplyr::count(r, sort = TRUE)

  checks <- tibble::tibble(
    check = c(
      "combined workbook contains all source study rows",
      "prepared contrast keys are unique",
      "one protein-intake centre is used",
      "frozen analysis row count is preserved",
      "excluded update studies do not reach model data",
      "all model rows satisfy the primary eligibility gate",
      "function outcome directions are complete",
      "all four outcome datasets are non-empty"
    ),
    passed = c(
      nrow(combined_workbook$study_characteristics) ==
        nrow(unique(combined_workbook$study_characteristics["study_no"])),
      nrow(duplicate_keys) == 0L,
      length(centres) == 1L,
      nrow(frozen_model_rows) == nrow(frozen_analysis),
      !any(
        model_data$study_no %in%
          decisions$study_no[!decisions$include_quantitative]
      ),
      all(model_data$analysis_eligible_primary & model_data$include_quantitative),
      all(!is.na(data_by_outcome$`function`$sign)),
      all(vapply(data_by_outcome, nrow, integer(1L)) > 0L)
    )
  )

  if (any(!checks$passed)) {
    stop(
      "Combined primary preparation failed: ",
      paste(checks$check[!checks$passed], collapse = "; "),
      call. = FALSE
    )
  }

  summary <- list(
    checks = checks,
    preprocessing = tibble::tibble(
      protein_intake_centre = centres[[1]],
      pre_post_r_imputation = if (nrow(imputed_r)) imputed_r$r[[1]] else NA_real_
    ),
    outcome_counts = purrr::imap_dfr(
      data_by_outcome,
      ~ tibble::tibble(
        outcome = .y,
        studies = dplyr::n_distinct(.x$study_no),
        arms = dplyr::n_distinct(.x$arm_no),
        contrasts = nrow(.x),
        frozen_contrasts = sum(
          .x$update_2026_data_source == "frozen canonical"
        ),
        update_contrasts = sum(
          .x$update_2026_data_source == "updated search"
        )
      )
    ),
    update_study_counts = reviewed_data |>
      dplyr::filter(study_no %in% decisions$study_no) |>
      dplyr::group_by(study_no) |>
      dplyr::summarise(
        candidate_contrasts = dplyr::n(),
        eligible_contrasts = sum(analysis_eligible_primary),
        included_contrasts = sum(update_2026_analysis_ready),
        .groups = "drop"
      ) |>
      dplyr::right_join(decisions, by = "study_no") |>
      dplyr::mutate(
        dplyr::across(
          c(candidate_contrasts, eligible_contrasts, included_contrasts),
          ~ tidyr::replace_na(.x, 0L)
        )
      ),
    frozen_effect_size_change = frozen_model_rows |>
      dplyr::inner_join(
        frozen_analysis |>
          dplyr::select(
            study_no, arm_no_extracted, outcome_group,
            measurement_method, unit_of_measurement, time1, time2,
            frozen_yi = yi, frozen_vi = vi
          ),
        by = key_columns
      ) |>
      dplyr::summarise(
        rows = dplyr::n(),
        yi_changed = sum(abs(yi - frozen_yi) > 1e-12),
        vi_changed = sum(abs(vi - frozen_vi) > 1e-12),
        max_abs_yi_change = max(abs(yi - frozen_yi)),
        max_abs_vi_change = max(abs(vi - frozen_vi))
      ),
    total_model_rows = nrow(model_data),
    update_model_rows = nrow(update_model_rows)
  )

  summary
}

write_combined_primary_preparation_audit <- function(
  audit,
  output_dir = here::here("results", "update_2026", "primary_refit")
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  artifacts <- audit[c(
    "checks", "preprocessing", "outcome_counts",
    "update_study_counts", "frozen_effect_size_change"
  )]

  purrr::imap_chr(
    artifacts,
    function(data, name) {
      path <- file.path(output_dir, paste0(name, ".csv"))
      readr::write_csv(data, path)
      path
    }
  ) |>
    unname()
}

make_update_2026_primary_fit_specifications <- function(
  primary_specifications,
  first_remedy_specifications,
  second_remedy_specifications,
  third_remedy_specifications
) {
  materialise_final_primary_specifications(
    primary_specifications,
    first_remedy_specifications,
    second_remedy_specifications,
    third_remedy_specifications
  )
}

prepare_update_2026_primary_candidate_input <- function(
  specification,
  data_by_outcome,
  priors_by_outcome
) {
  outcome <- specification$outcome
  data <- data_by_outcome[[outcome]]
  priors <- priors_by_outcome[[outcome]]

  if (is.null(data) || !nrow(data) || is.null(priors)) {
    stop(
      "Missing data or priors for updated primary outcome: ", outcome,
      call. = FALSE
    )
  }

  priors <- normalize_primary_brms_priors(priors)
  fitting_specification <- specification

  if (isTRUE(specification$third_remedy)) {
    priors <- add_primary_function_tail_priors(
      priors = priors,
      include_interactions =
        specification$targeted_piecewise_interaction_priors
    )
    priors <- c(
      priors,
      brms::set_prior("student_t(3, 0, 1)", class = "sd")
    )
    fitting_specification$regularised_variance_priors <- FALSE
    return(prepare_primary_remedied_candidate(
      specification = fitting_specification,
      priors = priors,
      data = data
    ))
  }

  if (isTRUE(specification$second_remedy)) {
    if (isTRUE(specification$targeted_piecewise_slope_priors)) {
      priors <- add_primary_piecewise_slope_priors(priors)
    }
    if (
      isTRUE(specification$targeted_piecewise_slope_priors) &&
        isTRUE(specification$regularised_variance_priors)
    ) {
      priors <- c(
        priors,
        brms::set_prior("student_t(3, 0, 1)", class = "sd")
      )
      fitting_specification$regularised_variance_priors <- FALSE
    }
    return(prepare_primary_remedied_candidate(
      specification = fitting_specification,
      priors = priors,
      data = data
    ))
  }

  if (isTRUE(specification$requires_refit)) {
    return(prepare_primary_remedied_candidate(
      specification = specification,
      priors = priors,
      data = data
    ))
  }

  if (isTRUE(specification$regularised_variance_priors)) {
    priors <- add_regularised_variance_priors(priors)
  }

  list(
    formula = specification$formula,
    priors = priors,
    data = data
  )
}

validate_update_2026_primary_fit_inputs <- function(
  specifications,
  data_by_outcome,
  priors_by_outcome
) {
  checks <- purrr::map_dfr(
    specifications,
    function(specification) {
      prepared <- prepare_update_2026_primary_candidate_input(
        specification = specification,
        data_by_outcome = data_by_outcome,
        priors_by_outcome = priors_by_outcome
      )
      brms::make_standata(
        formula = stats::as.formula(prepared$formula),
        data = prepared$data,
        prior = prepared$priors
      )

      tibble::tibble(
        outcome = specification$outcome,
        model_index = as.integer(specification$model_index),
        model_name = specification$model_name,
        stage = dplyr::case_when(
          isTRUE(specification$third_remedy) ~ "third_remedy",
          isTRUE(specification$second_remedy) ~ "second_remedy",
          isTRUE(specification$requires_refit) ~ "first_remedy",
          TRUE ~ "original"
        ),
        valid = TRUE
      )
    }
  )

  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Updated primary input preflight is incomplete.", call. = FALSE)
  }

  checks
}

write_update_2026_primary_input_preflight <- function(
  checks,
  output_file = here::here(
    "results", "update_2026", "primary_refit", "input_preflight.csv"
  )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete primary input preflight.", call. = FALSE)
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_primary_candidate_to_file <- function(
  specification,
  data_by_outcome,
  priors_by_outcome,
  output_root = here::here(
    "model_archive", "update_2026", "primary", "candidates"
  )
) {
  outcome <- specification$outcome
  model_index <- as.integer(specification$model_index)
  data <- data_by_outcome[[outcome]]
  priors <- priors_by_outcome[[outcome]]

  if (is.null(data) || !nrow(data) || is.null(priors)) {
    stop(
      "Missing data or priors for updated primary outcome: ", outcome,
      call. = FALSE
    )
  }

  # Frozen prior targets may have been created by an earlier brms release.
  # Rebuild their schema at the refit boundary so current brms can combine
  # them with any remedy-specific priors without column-mismatch failures.
  priors <- normalize_primary_brms_priors(priors)

  candidate <- if (isTRUE(specification$third_remedy)) {
    fit_primary_third_remedied_candidate_model(
      specification = specification,
      priors = priors,
      data = data
    )
  } else if (isTRUE(specification$second_remedy)) {
    fit_primary_second_remedied_candidate_model(
      specification = specification,
      priors = priors,
      data = data
    )
  } else if (isTRUE(specification$requires_refit)) {
    fit_primary_remedied_candidate_model(
      specification = specification,
      priors = priors,
      data = data
    )
  } else {
    fit_primary_candidate_model(
      specification = specification,
      priors = priors,
      data = data
    )
  }

  output_dir <- file.path(output_root, outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", model_index)
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  qs::qsave(candidate, temporary_file, preset = "high")

  if (file.exists(output_file)) file.remove(output_file)
  if (!file.rename(temporary_file, output_file)) {
    stop("Could not finalise updated primary fit file.", call. = FALSE)
  }

  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated primary fit file was not written correctly.", call. = FALSE)
  }

  output_file
}

diagnose_update_2026_primary_fit_file <- function(path, ppc_draws = 400L) {
  candidate <- qs::qread(path)
  get_primary_candidate_model_diagnostics(candidate, ppc_draws = ppc_draws)
}

compile_update_2026_primary_diagnostic_audit <- function(diagnostics) {
  split_diagnostics <- split(
    diagnostics,
    purrr::map_chr(diagnostics, ~ .x$summary$outcome[[1L]])
  )
  required <- c("whole_body", "hypertrophy", "strength", "function")

  if (
    !setequal(names(split_diagnostics), required) ||
      any(vapply(split_diagnostics[required], length, integer(1L)) != 26L)
  ) {
    stop("Updated primary diagnostics are incomplete.", call. = FALSE)
  }

  compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
}

make_update_2026_primary_remedy_specifications <- function(
  specifications,
  diagnostic_audit
) {
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass)

  if (
    nrow(failed) != 1L ||
      failed$outcome[[1L]] != "function" ||
      failed$model_index[[1L]] != 10L ||
      failed$divergent_transitions[[1L]] != 1L ||
      !isTRUE(failed$convergence_pass[[1L]])
  ) {
    stop(
      "The updated primary remedy registry no longer matches its audit.",
      call. = FALSE
    )
  }

  selected <- purrr::keep(
    specifications,
    ~ identical(.x$outcome, "function") && .x$model_index == 10L
  )
  if (length(selected) != 1L) {
    stop("Could not uniquely select function model 10.", call. = FALSE)
  }

  remedied <- selected[[1L]]
  remedied$adapt_delta <- 0.9995
  remedied$update_2026_remedy_stage <- "remedy_1"
  remedied$update_2026_remedy_reason <-
    "One divergent transition at adapt_delta 0.999."

  list(function_model_10 = remedied)
}

fit_update_2026_primary_remedy_to_file <- function(
  specification,
  data_by_outcome,
  priors_by_outcome
) {
  fit_update_2026_primary_candidate_to_file(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    output_root = here::here(
      "model_archive", "update_2026", "primary", "remedies", "remedy_1"
    )
  )
}

compile_update_2026_primary_remedy_diagnostic_audit <- function(diagnostics) {
  if (length(diagnostics) != 1L) {
    stop("Expected one updated primary remedy diagnostic.", call. = FALSE)
  }

  summary <- diagnostics[[1L]]$summary
  if (
    nrow(summary) != 1L ||
      summary$outcome[[1L]] != "function" ||
      summary$model_index[[1L]] != 10L
  ) {
    stop("The updated primary remedy diagnostic is invalid.", call. = FALSE)
  }

  list(
    diagnostic_summary = summary,
    parameter_failures = diagnostics[[1L]]$parameter_failures,
    diagnostic_pass = isTRUE(summary$diagnostic_pass[[1L]])
  )
}

write_update_2026_primary_remedy_diagnostic_audit <- function(
  audit,
  output_dir = here::here(
    "results", "update_2026", "primary_refit", "remedy_1"
  )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    parameter_failures = file.path(output_dir, "parameter_failures.csv")
  )
  readr::write_csv(audit$diagnostic_summary, paths[["diagnostic_summary"]])
  readr::write_csv(audit$parameter_failures, paths[["parameter_failures"]])
  unname(paths)
}

make_update_2026_primary_accepted_fit_manifest <- function(
  initial_fit_paths,
  remedy_fit_paths,
  remedy_diagnostic_audit
) {
  path_manifest <- function(paths, stage) {
    paths <- unlist(paths, use.names = FALSE)
    tibble::tibble(
      outcome = basename(dirname(paths)),
      model_index = as.integer(stringr::str_match(
        basename(paths), "model_([0-9]+)\\.qs$"
      )[, 2]),
      fit_path = paths,
      stage = stage
    )
  }

  initial <- path_manifest(initial_fit_paths, "initial")
  remedy <- path_manifest(remedy_fit_paths, "remedy_1")

  if (
    nrow(initial) != 104L ||
      nrow(remedy) != 1L ||
      remedy$outcome[[1L]] != "function" ||
      remedy$model_index[[1L]] != 10L ||
      !isTRUE(remedy_diagnostic_audit$diagnostic_pass)
  ) {
    stop(
      "Accepted updated primary fits cannot be assembled until the remedy passes.",
      call. = FALSE
    )
  }

  accepted <- initial |>
    dplyr::anti_join(
      remedy |> dplyr::select(outcome, model_index),
      by = c("outcome", "model_index")
    ) |>
    dplyr::bind_rows(remedy) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  if (
    nrow(accepted) != 104L ||
      anyDuplicated(accepted[c("outcome", "model_index")]) ||
      !all(file.exists(accepted$fit_path)) ||
      any(file.info(accepted$fit_path)$size <= 0)
  ) {
    stop("The accepted updated primary fit manifest is invalid.", call. = FALSE)
  }

  accepted
}

write_update_2026_primary_accepted_fit_manifest <- function(
  manifest,
  output_file = here::here(
    "results", "update_2026", "primary_refit", "accepted_fit_manifest.csv"
  )
) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

make_update_2026_primary_loo_records <- function(accepted_fit_manifest) {
  required <- c("outcome", "model_index", "fit_path", "stage")
  if (
    !all(required %in% names(accepted_fit_manifest)) ||
      nrow(accepted_fit_manifest) != 104L ||
      anyDuplicated(accepted_fit_manifest[c("outcome", "model_index")]) ||
      !all(file.exists(accepted_fit_manifest$fit_path))
  ) {
    stop("The updated primary LOO fit manifest is invalid.", call. = FALSE)
  }

  records <- purrr::pmap(
    accepted_fit_manifest[required],
    function(outcome, model_index, fit_path, stage) {
      list(
        outcome = as.character(outcome),
        model_index = as.integer(model_index),
        fit_path = as.character(fit_path),
        fit_stage = as.character(stage)
      )
    }
  )
  names(records) <- paste0(
    accepted_fit_manifest$outcome,
    "_model_",
    stringr::str_pad(accepted_fit_manifest$model_index, 2L, pad = "0")
  )
  records
}

update_2026_primary_loo_file_path <- function(outcome, model_index) {
  directory <- here::here(
    "model_archive", "update_2026", "primary", "loo", "integrated",
    outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf("model_%02d_integrated_loo.qs", as.integer(model_index))
  )
}

compute_update_2026_primary_integrated_loo_to_file <- function(
  record,
  cores = 1L,
  k_threshold = 0.7
) {
  if (
    !file.exists(record$fit_path) ||
      !record$model_index %in% seq_len(26L) ||
      !record$outcome %in% c(
        "whole_body", "hypertrophy", "strength", "function"
      )
  ) {
    stop("An updated primary LOO record is invalid.", call. = FALSE)
  }

  candidate <- qs::qread(record$fit_path)
  if (
    !identical(candidate$outcome, record$outcome) ||
      as.integer(candidate$model_index) != record$model_index
  ) {
    stop("An updated primary fit does not match its LOO record.", call. = FALSE)
  }

  initial <- compute_primary_candidate_psis_loo(
    fitted_candidate = candidate,
    cores = cores,
    save_psis = FALSE
  )
  integrated <- compute_primary_candidate_integrated_psis_loo(
    fitted_candidate = candidate,
    initial_loo_result = initial,
    k_threshold = k_threshold,
    cores = cores
  )
  integrated$fit_path <- normalizePath(
    record$fit_path, winslash = "/", mustWork = TRUE
  )
  integrated$fit_stage <- record$fit_stage

  output_file <- update_2026_primary_loo_file_path(
    record$outcome,
    record$model_index
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(integrated, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated primary LOO reduction.", call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)

  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

audit_update_2026_primary_integrated_loo_files <- function(
  paths,
  k_threshold = 0.7
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (
    length(paths) != 104L || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0)
  ) {
    stop("The updated primary integrated-LOO files are incomplete.",
         call. = FALSE)
  }

  results <- purrr::map(paths, qs::qread)
  valid <- purrr::map_lgl(results, function(result) {
    inherits(result$loo, "loo") &&
      identical(result$remedy, "integrated_observation_level_psis") &&
      result$outcome %in% c(
        "whole_body", "hypertrophy", "strength", "function"
      ) &&
      result$model_index %in% seq_len(26L) &&
      nrow(result$observation_manifest) == nrow(result$loo$pointwise)
  })
  keys <- purrr::map_chr(
    results,
    ~ paste(.x$outcome, .x$model_index, sep = "|")
  )
  expected_keys <- as.vector(outer(
    c("whole_body", "hypertrophy", "strength", "function"),
    seq_len(26L),
    function(outcome, model_index) paste(outcome, model_index, sep = "|")
  ))
  if (!all(valid) || anyDuplicated(keys) || !setequal(keys, expected_keys)) {
    stop("The updated primary integrated-LOO set is invalid.", call. = FALSE)
  }

  method_comparison <- purrr::map_dfr(results, "comparison")
  observation_diagnostics <- purrr::map_dfr(
    results,
    "observation_diagnostics"
  )
  model_summary <- method_comparison |>
    dplyr::filter(method == "integrated_observation_level_psis") |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  high_k_cases <- observation_diagnostics |>
    dplyr::filter(integrated_pareto_k > k_threshold) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  outcome_summary <- model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_high_k = sum(n_pareto_k_gt_threshold > 0L),
      high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_pareto_k = max(maximum_pareto_k),
      compute_hours = sum(compute_seconds) / 3600,
      stacking_ready = all(n_pareto_k_gt_threshold == 0L),
      .groups = "drop"
    )
  reduction_manifest <- purrr::map2_dfr(
    paths,
    results,
    function(path, result) {
      tibble::tibble(
        outcome = result$outcome,
        model_index = as.integer(result$model_index),
        model_name = result$model_name,
        fit_stage = result$fit_stage,
        fit_path = result$fit_path,
        integrated_loo_path = normalizePath(
          path, winslash = "/", mustWork = TRUE
        ),
        bytes = file.info(path)$size,
        sha256 = digest::digest(file = path, algo = "sha256")
      )
    }
  ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  if (
    nrow(model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      nrow(method_comparison) != 208L ||
      nrow(reduction_manifest) != 104L
  ) {
    stop("The updated primary integrated-LOO audit is incomplete.",
         call. = FALSE)
  }

  list(
    reduction_manifest = reduction_manifest,
    model_summary = model_summary,
    outcome_summary = outcome_summary,
    method_comparison = method_comparison,
    observation_diagnostics = observation_diagnostics,
    high_k_cases = high_k_cases,
    exact_refits_required = nrow(high_k_cases),
    stacking_ready = nrow(high_k_cases) == 0L,
    k_threshold = k_threshold
  )
}

write_update_2026_primary_integrated_loo_audit <- function(
  audit,
  output_dir = here::here(
    "results", "update_2026", "primary_refit", "loo_psis"
  )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    method_comparison = file.path(output_dir, "method_comparison.csv"),
    observation_diagnostics = file.path(
      output_dir, "observation_diagnostics.csv.gz"
    ),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(audit$reduction_manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$method_comparison, paths[["method_comparison"]])
  readr::write_csv(
    audit$observation_diagnostics,
    paths[["observation_diagnostics"]]
  )
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

make_update_2026_primary_exact_loo_records <- function(integrated_loo_audit) {
  cases <- integrated_loo_audit$high_k_cases |>
    dplyr::select(outcome, model_index, model_name, observation) |>
    dplyr::distinct() |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    ) |>
    dplyr::mutate(case_id = dplyr::row_number(), .before = 1L) |>
    dplyr::left_join(
      integrated_loo_audit$reduction_manifest |>
        dplyr::select(
          outcome, model_index, fit_path, integrated_loo_path
        ),
      by = c("outcome", "model_index")
    )

  expected <- tibble::tribble(
    ~outcome, ~model_index, ~observation,
    "whole_body", 3L, 3L,
    "whole_body", 19L, 3L,
    "hypertrophy", 13L, 253L
  )
  observed <- cases |>
    dplyr::select(outcome, model_index, observation)
  if (
    !isTRUE(all.equal(observed, expected, check.attributes = FALSE)) ||
      anyNA(cases$fit_path) || anyNA(cases$integrated_loo_path) ||
      any(!file.exists(cases$fit_path)) ||
      any(!file.exists(cases$integrated_loo_path))
  ) {
    stop("The updated primary exact-LOO registry no longer matches its audit.",
         call. = FALSE)
  }

  records <- purrr::pmap(cases, list)
  names(records) <- paste0(
    cases$outcome, "_model_",
    stringr::str_pad(cases$model_index, 2L, pad = "0"),
    "_observation_", cases$observation
  )
  records
}

validate_update_2026_primary_exact_loo_inputs <- function(records) {
  checks <- purrr::map_dfr(records, function(record) {
    candidate <- qs::qread(record$fit_path)
    integrated <- qs::qread(record$integrated_loo_path)
    model_index <- as.integer(record$model_index)
    coefficient_count <- if (model_index %in% c(9L, 12L, 21L, 24L)) {
      2L
    } else if (model_index %in% c(10L, 13L, 22L, 25L)) {
      3L
    } else {
      1L
    }
    expected_design <- if (coefficient_count == 1L) {
      matrix(1, nrow(candidate$model$data), 1L)
    } else if (coefficient_count == 2L) {
      cbind(1, candidate$model$data$protein_intake)
    } else {
      cbind(
        1,
        as.matrix(lspline::lspline(candidate$model$data$protein_intake, 0.5))
      )
    }
    standata <- brms::standata(candidate$model)
    design_names <- paste0("Z_3_", seq_len(coefficient_count))
    if (!all(design_names %in% names(standata))) {
      stop("An exact-LOO study design column is missing.", call. = FALSE)
    }
    observed_design <- do.call(
      cbind,
      lapply(design_names, function(name) as.numeric(standata[[name]]))
    )
    maximum_design_difference <- max(abs(
      observed_design - expected_design
    ))
    pareto_k <- integrated$loo$diagnostics$pareto_k[[record$observation]]
    valid <-
      identical(candidate$outcome, record$outcome) &&
      as.integer(candidate$model_index) == model_index &&
      identical(integrated$outcome, record$outcome) &&
      as.integer(integrated$model_index) == model_index &&
      nrow(expected_design) == nrow(integrated$observation_manifest) &&
      is.finite(maximum_design_difference) &&
      maximum_design_difference < 1e-12 &&
      is.finite(pareto_k) && pareto_k > integrated$k_threshold

    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = model_index,
      observation = as.integer(record$observation),
      study_structure = integrated$study_structure,
      pareto_k = pareto_k,
      maximum_design_difference = maximum_design_difference,
      valid = valid
    )
  })
  if (nrow(checks) != 3L || !all(checks$valid)) {
    stop("The updated primary exact-LOO input preflight failed.",
         call. = FALSE)
  }
  checks
}

write_update_2026_primary_exact_loo_input_preflight <- function(
  checks,
  output_file = here::here(
    "results", "update_2026", "primary_refit", "loo_exact",
    "input_preflight.csv"
  )
) {
  if (nrow(checks) != 3L || !all(checks$valid)) {
    stop("Cannot write an invalid exact-LOO input preflight.", call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

get_update_2026_primary_exact_components <- function(
  refitted_model,
  original_data,
  model_index
) {
  common_variables <- c(
    residual = "sigma",
    arm = "sd_arm_no__Intercept",
    outcome = "sd_outcome_no__Intercept"
  )
  linear_slope_indices <- c(9L, 12L, 21L, 24L)
  spline_slope_indices <- c(10L, 13L, 22L, 25L)
  factor_smooth_indices <- c(11L, 14L, 23L, 26L)
  if (model_index %in% factor_smooth_indices) {
    stop("Updated exact scoring for factor smooths is not required by this audit.",
         call. = FALSE)
  }

  coefficient_labels <- if (model_index %in% linear_slope_indices) {
    c("Intercept", "protein_intake")
  } else if (model_index %in% spline_slope_indices) {
    c(
      "Intercept",
      "lsplineprotein_intake0.51",
      "lsplineprotein_intake0.52"
    )
  } else {
    "Intercept"
  }
  study_design <- if (length(coefficient_labels) == 1L) {
    matrix(1, nrow = nrow(original_data), ncol = 1L)
  } else if (model_index %in% linear_slope_indices) {
    cbind(1, original_data$protein_intake)
  } else {
    cbind(
      1,
      as.matrix(lspline::lspline(original_data$protein_intake, 0.5))
    )
  }
  if (
    nrow(study_design) != nrow(original_data) ||
      ncol(study_design) != length(coefficient_labels) ||
      anyNA(study_design) || any(!is.finite(study_design))
  ) {
    stop("The exact-LOO study design is invalid.", call. = FALSE)
  }

  study_sd_variables <- paste0("sd_study_no__", coefficient_labels)
  model_variables <- posterior::variables(refitted_model)
  correlation_variables <- character()
  correlation_pairs <- list()
  if (length(coefficient_labels) > 1L) {
    for (first in seq_len(length(coefficient_labels) - 1L)) {
      for (second in seq.int(first + 1L, length(coefficient_labels))) {
        forward <- paste0(
          "cor_study_no__", coefficient_labels[[first]], "__",
          coefficient_labels[[second]]
        )
        reverse <- paste0(
          "cor_study_no__", coefficient_labels[[second]], "__",
          coefficient_labels[[first]]
        )
        present <- intersect(c(forward, reverse), model_variables)
        if (length(present) > 1L) {
          stop("An exact-LOO study correlation is duplicated.", call. = FALSE)
        }
        if (length(present) == 1L) {
          correlation_variables <- c(correlation_variables, present)
          correlation_pairs[[length(correlation_pairs) + 1L]] <-
            c(first, second)
        }
      }
    }
  }
  expected_correlations <-
    length(coefficient_labels) * (length(coefficient_labels) - 1L) / 2L
  if (
    length(correlation_variables) != 0L &&
      length(correlation_variables) != expected_correlations
  ) {
    stop("The exact-LOO study correlation matrix is incomplete.",
         call. = FALSE)
  }

  draws <- extract_primary_exact_draws(
    refitted_model,
    c(common_variables, study_sd_variables, correlation_variables)
  )
  n_draws <- nrow(draws)
  n_coefficients <- length(coefficient_labels)
  study_basis <- lapply(
    seq_len(n_coefficients),
    function(x) matrix(NA_real_, n_draws, nrow(original_data))
  )
  if (n_coefficients == 1L) {
    study_basis[[1L]] <- outer(
      draws[[study_sd_variables[[1L]]]],
      study_design[, 1L]
    )
    structure_name <- "study_intercept"
  } else {
    for (draw in seq_len(n_draws)) {
      correlation <- diag(n_coefficients)
      if (length(correlation_variables) > 0L) {
        for (pair_index in seq_along(correlation_variables)) {
          pair <- correlation_pairs[[pair_index]]
          value <- draws[[correlation_variables[[pair_index]]]][[draw]]
          correlation[pair[[1L]], pair[[2L]]] <- value
          correlation[pair[[2L]], pair[[1L]]] <- value
        }
      }
      standard_deviations <- as.numeric(as.matrix(
        draws[draw, study_sd_variables, drop = FALSE]
      ))
      covariance <-
        diag(standard_deviations, n_coefficients) %*%
        correlation %*%
        diag(standard_deviations, n_coefficients)
      cholesky <- tryCatch(t(chol(covariance)), error = function(e) NULL)
      if (is.null(cholesky)) {
        stop("An exact-LOO study covariance is not positive definite.",
             call. = FALSE)
      }
      transformed_design <- study_design %*% cholesky
      for (component in seq_len(n_coefficients)) {
        study_basis[[component]][draw, ] <-
          transformed_design[, component]
      }
    }
    structure_root <- if (model_index %in% linear_slope_indices) {
      "study_intercept_and_linear_slope"
    } else {
      "study_intercept_and_linear_spline_slopes"
    }
    structure_name <- paste0(
      if (length(correlation_variables) > 0L) "correlated_" else "uncorrelated_",
      structure_root
    )
  }

  linear_predictor <- brms::posterior_linpred(
    refitted_model,
    newdata = original_data,
    re_formula = NA,
    transform = FALSE,
    allow_new_levels = TRUE
  )
  if (
    !identical(dim(linear_predictor), c(n_draws, nrow(original_data))) ||
      anyNA(linear_predictor) || any(!is.finite(linear_predictor))
  ) {
    stop("The exact-LOO population predictor is invalid.", call. = FALSE)
  }

  list(
    structure = structure_name,
    linear_predictor = linear_predictor,
    residual_sd = draws[[common_variables[["residual"]]]],
    arm_intercept_sd = draws[[common_variables[["arm"]]]],
    outcome_intercept_sd = draws[[common_variables[["outcome"]]]],
    study_basis = study_basis
  )
}

compute_update_2026_primary_exact_loo_case <- function(refit, components = NULL) {
  integrated_result <- refit$integrated_result
  manifest <- integrated_result$observation_manifest
  observation <- as.integer(refit$observation)
  if (
    observation < 1L || observation > nrow(manifest) ||
      nrow(refit$original_data) != nrow(manifest)
  ) {
    stop("The updated exact-LOO scoring inputs are invalid.", call. = FALSE)
  }

  if (is.null(components)) {
    components <- get_update_2026_primary_exact_components(
      refitted_model = refit$refitted_model,
      original_data = refit$original_data,
      model_index = as.integer(refit$model_index)
    )
  }
  if (!identical(components$structure, integrated_result$study_structure)) {
    stop("The exact-LOO study structure changed during refitting.",
         call. = FALSE)
  }
  integrated_log_lik <- integrate_primary_study_basis_log_lik(
    linear_predictor = components$linear_predictor,
    observed = manifest$yi,
    known_variance = manifest$vi,
    study_id = manifest$study_no,
    arm_id = manifest$arm_no,
    outcome_id = manifest$outcome_no,
    residual_sd = components$residual_sd,
    arm_intercept_sd = components$arm_intercept_sd,
    outcome_intercept_sd = components$outcome_intercept_sd,
    study_basis = components$study_basis
  )
  held_out_log_lik <- integrated_log_lik[, observation]
  maximum_log_lik <- max(held_out_log_lik)
  scaled_likelihood <- exp(held_out_log_lik - maximum_log_lik)
  exact_elpd <- maximum_log_lik + log(mean(scaled_likelihood))
  exact_mcse <- tryCatch(
    posterior::mcse_mean(scaled_likelihood) / mean(scaled_likelihood),
    error = function(e) NA_real_
  )

  patched_loo <- integrated_result$loo
  old_elpd <- patched_loo$pointwise[observation, "elpd_loo"]
  old_p_loo <- patched_loo$pointwise[observation, "p_loo"]
  full_data_lpd <- old_elpd + old_p_loo
  patched_loo$pointwise[observation, "elpd_loo"] <- exact_elpd
  patched_loo$pointwise[observation, "p_loo"] <- full_data_lpd - exact_elpd
  patched_loo$pointwise[observation, "looic"] <- -2 * exact_elpd
  if ("mcse_elpd_loo" %in% colnames(patched_loo$pointwise)) {
    patched_loo$pointwise[observation, "mcse_elpd_loo"] <- exact_mcse
  }
  old_pareto_k <- patched_loo$diagnostics$pareto_k[[observation]]
  patched_loo$diagnostics$pareto_k[[observation]] <- NA_real_
  if (
    !is.null(patched_loo$diagnostics$n_eff) &&
      length(patched_loo$diagnostics$n_eff) >= observation
  ) {
    patched_loo$diagnostics$n_eff[[observation]] <- NA_real_
  }
  patched_loo <- update_primary_exact_loo_estimates(patched_loo)

  initial_summary <- integrated_result$comparison |>
    dplyr::filter(method == "initial_conditional_psis")
  final_summary <- summarise_primary_loo_method(
    patched_loo,
    method = "integrated_observation_level_psis",
    k_threshold = integrated_result$k_threshold
  ) |>
    dplyr::mutate(
      compute_seconds =
        integrated_result$compute_seconds + refit$refit_seconds,
      delta_elpd_vs_initial =
        integrated_result$comparison$delta_elpd_vs_initial[
          integrated_result$comparison$method ==
            "integrated_observation_level_psis"
        ] + exact_elpd - old_elpd,
      se_delta_elpd_vs_initial = NA_real_,
      exact_loo_refits = 1L,
      exact_refit_seconds = refit$refit_seconds
    )
  comparison <- dplyr::bind_rows(
    dplyr::mutate(
      initial_summary,
      exact_loo_refits = 0L,
      exact_refit_seconds = NA_real_
    ),
    final_summary
  ) |>
    dplyr::mutate(
      outcome = integrated_result$outcome,
      model_index = integrated_result$model_index,
      model_name = integrated_result$model_name,
      study_structure = integrated_result$study_structure,
      .before = 1L
    )
  observation_diagnostics <- integrated_result$observation_diagnostics |>
    dplyr::mutate(
      exact_loo_refit = dplyr::row_number() == observation,
      exact_elpd_loo = dplyr::if_else(
        exact_loo_refit, exact_elpd, NA_real_
      ),
      exact_elpd_mcse = dplyr::if_else(
        exact_loo_refit, exact_mcse, NA_real_
      ),
      integrated_pareto_k = dplyr::if_else(
        exact_loo_refit, NA_real_, integrated_pareto_k
      ),
      integrated_resolved_initial_warning = dplyr::if_else(
        exact_loo_refit,
        initial_pareto_k > integrated_result$k_threshold,
        integrated_resolved_initial_warning
      ),
      integrated_new_warning = dplyr::if_else(
        exact_loo_refit, FALSE, integrated_new_warning
      )
    )
  diagnostics <- summarise_primary_exact_loo_refit_diagnostics(refit)
  exact_summary <- diagnostics |>
    dplyr::mutate(
      study_structure = integrated_result$study_structure,
      old_psis_pareto_k = old_pareto_k,
      old_psis_elpd_loo = old_elpd,
      exact_elpd_loo = exact_elpd,
      exact_minus_psis_elpd = exact_elpd - old_elpd,
      exact_elpd_mcse = exact_mcse
    )

  list(
    outcome = integrated_result$outcome,
    model_index = integrated_result$model_index,
    model_name = integrated_result$model_name,
    study_structure = integrated_result$study_structure,
    observation_manifest = manifest,
    loo = patched_loo,
    remedy = "integrated_observation_level_psis_plus_exact",
    k_threshold = integrated_result$k_threshold,
    compute_seconds =
      integrated_result$compute_seconds + refit$refit_seconds,
    comparison = comparison,
    observation_diagnostics = observation_diagnostics,
    exact_refit_summary = exact_summary,
    exact_loo_refits = 1L,
    fit_path = integrated_result$fit_path,
    fit_stage = integrated_result$fit_stage
  )
}

update_2026_primary_exact_loo_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "primary", "loo", "exact",
    record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$observation)
    )
  )
}

fit_and_score_update_2026_primary_exact_loo_to_file <- function(record) {
  candidate <- qs::qread(record$fit_path)
  integrated_result <- qs::qread(record$integrated_loo_path)
  if (
    !identical(candidate$outcome, record$outcome) ||
      as.integer(candidate$model_index) != record$model_index ||
      !identical(integrated_result$outcome, record$outcome) ||
      as.integer(integrated_result$model_index) != record$model_index ||
      integrated_result$loo$diagnostics$pareto_k[[record$observation]] <= 0.7
  ) {
    stop("An updated primary exact-LOO record failed validation.",
         call. = FALSE)
  }

  job <- list(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = as.integer(record$observation),
    fitted_candidate = candidate,
    integrated_result = integrated_result
  )
  refit <- fit_primary_integrated_exact_loo_case(job)
  correction <- compute_update_2026_primary_exact_loo_case(refit)

  output_file <- update_2026_primary_exact_loo_file_path(record)
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(correction, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated primary exact-LOO correction.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

audit_update_2026_primary_exact_loo_files <- function(paths) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (
    length(paths) != 3L || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0)
  ) {
    stop("The updated primary exact-LOO files are incomplete.",
         call. = FALSE)
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(outcome, levels = c("whole_body", "hypertrophy")),
      model_index,
      observation
    )
  correction_manifest <- purrr::map2_dfr(
    paths,
    corrections,
    function(path, correction) {
      tibble::tibble(
        outcome = correction$outcome,
        model_index = as.integer(correction$model_index),
        model_name = correction$model_name,
        observation = correction$exact_refit_summary$observation[[1L]],
        study_structure = correction$study_structure,
        correction_path = normalizePath(
          path, winslash = "/", mustWork = TRUE
        ),
        bytes = file.info(path)$size,
        sha256 = digest::digest(file = path, algo = "sha256")
      )
    }
  )
  if (
    nrow(diagnostic_summary) != 3L ||
      nrow(correction_manifest) != 3L ||
      anyDuplicated(correction_manifest[c(
        "outcome", "model_index", "observation"
      )])
  ) {
    stop("The updated primary exact-LOO audit is invalid.", call. = FALSE)
  }

  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = correction_manifest,
    exact_refits = 3L,
    exact_refits_passing = sum(diagnostic_summary$diagnostic_pass),
    exact_refits_requiring_remediation =
      sum(!diagnostic_summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass =
      all(diagnostic_summary$diagnostic_pass)
  )
}

write_update_2026_primary_exact_loo_audit <- function(
  audit,
  output_dir = here::here(
    "results", "update_2026", "primary_refit", "loo_exact"
  )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(audit$diagnostic_summary, paths[["diagnostic_summary"]])
  readr::write_csv(audit$correction_manifest, paths[["correction_manifest"]])
  unname(paths)
}

make_update_2026_primary_final_loo_manifest <- function(
  integrated_loo_audit,
  exact_loo_audit
) {
  if (!isTRUE(exact_loo_audit$all_exact_refit_diagnostics_pass)) {
    stop("Final LOO assembly requires all exact refits to pass.",
         call. = FALSE)
  }
  exact <- exact_loo_audit$correction_manifest |>
    dplyr::select(
      outcome, model_index, observation, study_structure,
      correction_path, correction_bytes = bytes,
      correction_sha256 = sha256
    )
  manifest <- integrated_loo_audit$reduction_manifest |>
    dplyr::left_join(exact, by = c("outcome", "model_index")) |>
    dplyr::mutate(
      loo_source = dplyr::if_else(
        is.na(correction_path), "integrated_psis", "exact_correction"
      ),
      final_loo_path = dplyr::coalesce(
        correction_path, integrated_loo_path
      ),
      final_bytes = dplyr::coalesce(correction_bytes, bytes),
      final_sha256 = dplyr::coalesce(correction_sha256, sha256)
    ) |>
    dplyr::select(
      outcome, model_index, model_name, fit_stage, fit_path,
      loo_source, observation, study_structure,
      integrated_loo_path, correction_path, final_loo_path,
      final_bytes, final_sha256
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  if (
    nrow(manifest) != 104L ||
      sum(manifest$loo_source == "exact_correction") != 3L ||
      anyDuplicated(manifest[c("outcome", "model_index")]) ||
      any(!file.exists(manifest$final_loo_path)) ||
      any(file.info(manifest$final_loo_path)$size <= 0)
  ) {
    stop("The updated primary final LOO manifest is invalid.",
         call. = FALSE)
  }
  manifest
}

audit_update_2026_primary_final_loo <- function(
  final_loo_manifest,
  exact_loo_audit,
  k_threshold = 0.7
) {
  results <- purrr::map(final_loo_manifest$final_loo_path, qs::qread)
  keys <- purrr::map_chr(
    results,
    ~ paste(.x$outcome, .x$model_index, sep = "|")
  )
  manifest_keys <- paste(
    final_loo_manifest$outcome,
    final_loo_manifest$model_index,
    sep = "|"
  )
  valid <- purrr::map_lgl(results, function(result) {
    inherits(result$loo, "loo") &&
      result$model_index %in% seq_len(26L) &&
      result$outcome %in% c(
        "whole_body", "hypertrophy", "strength", "function"
      ) &&
      nrow(result$loo$pointwise) == nrow(result$observation_manifest)
  })
  if (!all(valid) || !identical(keys, manifest_keys)) {
    stop("A final updated primary LOO object failed validation.",
         call. = FALSE)
  }

  method_comparison <- purrr::map_dfr(results, "comparison")
  model_summary <- method_comparison |>
    dplyr::filter(method == "integrated_observation_level_psis") |>
    dplyr::mutate(
      exact_loo_refits = tidyr::replace_na(exact_loo_refits, 0L),
      exact_refit_seconds = tidyr::replace_na(
        exact_refit_seconds, 0
      )
    )
  observation_diagnostics <- purrr::map_dfr(
    results,
    function(result) {
      diagnostics <- result$observation_diagnostics
      if (!"exact_loo_refit" %in% names(diagnostics)) {
        diagnostics <- diagnostics |>
          dplyr::mutate(
            exact_loo_refit = FALSE,
            exact_elpd_loo = NA_real_,
            exact_elpd_mcse = NA_real_
          )
      }
      diagnostics
    }
  )
  remaining_high_k <- observation_diagnostics |>
    dplyr::filter(integrated_pareto_k > k_threshold)
  exact_refit_summary <- exact_loo_audit$diagnostic_summary
  outcome_summary <- model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_pareto_k_warning = sum(
        n_pareto_k_gt_threshold > 0L
      ),
      remaining_pareto_k_warnings = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k),
      exact_loo_refits = sum(exact_loo_refits),
      total_compute_hours = sum(compute_seconds) / 3600,
      stacking_ready = all(n_pareto_k_gt_threshold == 0L),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      exact_refit_summary |>
        dplyr::group_by(outcome) |>
        dplyr::summarise(
          exact_refit_diagnostics_pass = all(diagnostic_pass),
          .groups = "drop"
        ),
      by = "outcome"
    ) |>
    dplyr::mutate(
      exact_refit_diagnostics_pass = tidyr::replace_na(
        exact_refit_diagnostics_pass, TRUE
      ),
      stacking_ready =
        stacking_ready & exact_refit_diagnostics_pass
    )
  stacking_ready <-
    nrow(remaining_high_k) == 0L &&
    all(exact_refit_summary$diagnostic_pass) &&
    all(outcome_summary$stacking_ready)
  if (
    nrow(model_summary) != 104L ||
      nrow(method_comparison) != 208L ||
      nrow(outcome_summary) != 4L ||
      nrow(exact_refit_summary) != 3L ||
      !stacking_ready
  ) {
    stop("The updated primary final LOO audit is not stacking-ready.",
         call. = FALSE)
  }

  list(
    final_loo_manifest = final_loo_manifest,
    model_summary = model_summary,
    outcome_summary = outcome_summary,
    method_comparison = method_comparison,
    observation_diagnostics = observation_diagnostics,
    exact_refit_summary = exact_refit_summary,
    remaining_high_k = remaining_high_k,
    stacking_ready = stacking_ready
  )
}

write_update_2026_primary_final_loo_audit <- function(
  audit,
  output_dir = here::here(
    "results", "update_2026", "primary_refit", "loo_final"
  )
) {
  if (!isTRUE(audit$stacking_ready)) {
    stop("Cannot write a non-final updated primary LOO audit.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    method_comparison = file.path(output_dir, "method_comparison.csv"),
    observation_diagnostics = file.path(
      output_dir, "observation_diagnostics.csv.gz"
    ),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(audit$final_loo_manifest, paths[["final_loo_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$method_comparison, paths[["method_comparison"]])
  readr::write_csv(
    audit$observation_diagnostics,
    paths[["observation_diagnostics"]]
  )
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

compute_update_2026_primary_stacking_weights <- function(final_loo_audit) {
  if (!isTRUE(final_loo_audit$stacking_ready)) {
    stop("Updated primary LOO is not ready for stacking.", call. = FALSE)
  }
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  final_results <- purrr::map(outcomes, function(outcome) {
    records <- final_loo_audit$final_loo_manifest |>
      dplyr::filter(.data$outcome == .env$outcome) |>
      dplyr::arrange(model_index)
    if (
      nrow(records) != 26L ||
        !identical(as.integer(records$model_index), seq_len(26L))
    ) {
      stop("An updated primary stacking outcome is incomplete.",
           call. = FALSE)
    }
    purrr::map(records$final_loo_path, qs::qread)
  })
  names(final_results) <- outcomes
  compute_primary_stacking_weights(final_results)
}

write_update_2026_primary_stacking_weights <- function(
  summary,
  output_file = here::here(
    "results", "update_2026", "primary_refit", "stacking_weights.csv"
  )
) {
  audit <- summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      weight_sum = sum(stacking_weight),
      nonzero_models = sum(stacking_weight > 1e-12),
      .groups = "drop"
    )
  if (
    nrow(summary) != 104L || nrow(audit) != 4L ||
      any(abs(audit$weight_sum - 1) > 1e-8)
  ) {
    stop("Updated primary stacking weights failed validation.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(summary, output_file)
  output_file
}

make_update_2026_primary_reporting_records <- function(
  accepted_fit_manifest,
  stacking_weights_summary
) {
  records <- accepted_fit_manifest |>
    dplyr::select(outcome, model_index, fit_path, fit_stage = stage) |>
    dplyr::left_join(
      stacking_weights_summary |>
        dplyr::select(
          outcome, model_index, model_name,
          stacking_weight
        ),
      by = c("outcome", "model_index")
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  if (
    nrow(records) != 104L || anyNA(records$stacking_weight) ||
      any(!file.exists(records$fit_path))
  ) {
    stop("Updated primary reporting inputs are incomplete.", call. = FALSE)
  }

  seeds <- c(
    whole_body = 1989L,
    hypertrophy = 1990L,
    strength = 1991L,
    `function` = 1992L
  )
  output <- purrr::map(names(seeds), function(outcome_name) {
    table <- records |>
      dplyr::filter(.data$outcome == .env$outcome_name) |>
      dplyr::arrange(model_index)
    if (
      nrow(table) != 26L ||
        !identical(as.integer(table$model_index), seq_len(26L)) ||
        abs(sum(table$stacking_weight) - 1) > 1e-8
    ) {
      stop("An updated primary reporting outcome is invalid.",
           call. = FALSE)
    }
    list(
      outcome = outcome_name,
      seed = seeds[[outcome_name]],
      models = table
    )
  })
  names(output) <- names(seeds)
  output
}

update_2026_primary_model_average_file_path <- function(outcome) {
  directory <- here::here(
    "model_archive", "update_2026", "primary", "reporting",
    "model_average"
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, paste0(outcome, ".qs"))
}

build_update_2026_primary_model_average_to_file <- function(
  record,
  data_by_outcome,
  protein_min = 0.8,
  grid_step = 0.01,
  n_draws = 4000L,
  output_file = update_2026_primary_model_average_file_path(record$outcome)
) {
  data <- data_by_outcome[[record$outcome]]
  table <- record$models
  if (is.null(data) || !nrow(data) || nrow(table) != 26L) {
    stop("An updated primary model-average input is invalid.", call. = FALSE)
  }
  weights <- table$stacking_weight
  names(weights) <- table$model_name
  prediction_grid <- make_primary_prediction_grid(
    data = data,
    protein_min = protein_min,
    grid_step = grid_step
  )
  counts <- allocate_primary_mixture_models(weights, n_draws)
  set.seed(record$seed)
  allocation <- sample(rep(seq_along(weights), times = counts))
  averaged_draws <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = nrow(prediction_grid$newdata)
  )

  for (model_position in seq_len(nrow(table))) {
    output_rows <- which(allocation == model_position)
    if (length(output_rows) == 0L) next
    candidate <- qs::qread(table$fit_path[[model_position]])
    if (
      !identical(candidate$outcome, record$outcome) ||
        as.integer(candidate$model_index) != table$model_index[[model_position]] ||
        !identical(candidate$model_name, table$model_name[[model_position]])
    ) {
      stop("An archived fit does not match its reporting record.",
           call. = FALSE)
    }
    draw_matrix <- brms::posterior_epred(
      candidate$model,
      newdata = prediction_grid$newdata,
      re_formula = NA,
      ndraws = length(output_rows)
    ) |>
      as.matrix()
    if (
      ncol(draw_matrix) != nrow(prediction_grid$newdata) &&
        nrow(draw_matrix) == nrow(prediction_grid$newdata)
    ) {
      draw_matrix <- t(draw_matrix)
    }
    if (
      !identical(
        dim(draw_matrix),
        c(length(output_rows), nrow(prediction_grid$newdata))
      ) ||
        anyNA(draw_matrix) || any(!is.finite(draw_matrix))
    ) {
      stop("Posterior predictions do not align with the reporting grid.",
           call. = FALSE)
    }
    averaged_draws[output_rows, ] <- draw_matrix
    rm(candidate, draw_matrix)
    gc(verbose = FALSE)
  }
  if (anyNA(averaged_draws)) {
    stop("The updated posterior mixture is incomplete.", call. = FALSE)
  }

  output <- list(
    outcome = record$outcome,
    draws = averaged_draws,
    grid = prediction_grid$keys,
    draw_manifest = tibble::tibble(
      drawid = seq_len(n_draws),
      model_index = table$model_index[allocation],
      model_name = table$model_name[allocation]
    ),
    model_weights = weights,
    settings = list(
      protein_median = prediction_grid$protein_median,
      protein_min = prediction_grid$protein_min,
      protein_max = prediction_grid$protein_max,
      grid_step = prediction_grid$grid_step,
      n_draws = n_draws,
      seed = record$seed
    )
  )
  class(output) <- c("primary_model_average", class(output))

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(output, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated primary model average.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

make_primary_arm_re_formula <- function(model) {
  if (!inherits(model, "brmsfit")) {
    stop("`model` must be a brms fit.", call. = FALSE)
  }
  random_terms <- brms::brmsterms(model$formula)$dpars$mu$re
  if (is.null(random_terms) || !nrow(random_terms)) {
    stop("A primary model has no group-level terms.", call. = FALSE)
  }
  random_terms <- random_terms |>
    dplyr::filter(.data$group %in% c("study_no", "arm_no"))
  if (!nrow(random_terms) || !"arm_no" %in% random_terms$group) {
    stop("A primary model has no arm-level term.", call. = FALSE)
  }

  formula_parts <- vapply(seq_len(nrow(random_terms)), function(index) {
    right_hand_side <- sub(
      "^~",
      "",
      paste(deparse(random_terms$form[[index]]), collapse = " ")
    )
    separator <- if (isTRUE(random_terms$cor[[index]])) "|" else "||"
    sprintf(
      "(%s %s %s)",
      right_hand_side,
      separator,
      random_terms$group[[index]]
    )
  }, character(1L))
  stats::as.formula(
    paste("~", paste(formula_parts, collapse = " + ")),
    env = environment(model$formula)
  )
}

make_primary_arm_prediction_ids <- function(data) {
  key_columns <- c(
    "study_no", "arm_no", "protein_intake_raw",
    "resistance_exercise_code", "duration_centre", "age_centre"
  )
  missing_columns <- setdiff(key_columns, names(data))
  if (length(missing_columns)) {
    stop(
      "Primary arm-prediction data are missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  keys <- lapply(data[key_columns], function(column) {
    if (is.numeric(column)) {
      return(ifelse(is.na(column), "<NA>", sprintf("%.17g", column)))
    }
    ifelse(is.na(column), "<NA>", as.character(column))
  })
  key <- do.call(paste, c(keys, sep = "\r"))
  match(key, unique(key))
}

validate_primary_arm_candidate_data <- function(candidate_data, data) {
  # brms retains only variables required by each candidate formula, so raw
  # exposure, moderator, and study columns are not present in every archived
  # fit. Arm and outcome identifiers are invariant across all 26 candidates
  # and verify the observation order used for posterior prediction; study is
  # checked too whenever the candidate retained it.
  required_key_columns <- c("arm_no", "outcome_no")
  optional_key_columns <- "study_no"
  if (
    nrow(candidate_data) != nrow(data) ||
      !all(required_key_columns %in% names(candidate_data)) ||
      !all(required_key_columns %in% names(data))
  ) {
    stop("An archived fit does not align with its primary data.",
         call. = FALSE)
  }
  key_columns <- c(
    required_key_columns,
    intersect(optional_key_columns, names(candidate_data))
  )
  aligned <- vapply(key_columns, function(variable) {
    identical(
      as.character(candidate_data[[variable]]),
      as.character(data[[variable]])
    )
  }, logical(1L))
  if (!all(aligned)) {
    stop(
      "An archived fit has misaligned primary observations: ",
      paste(key_columns[!aligned], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

update_2026_primary_arm_conditional_file_path <- function(outcome) {
  directory <- here::here(
    "model_archive", "update_2026", "primary", "reporting",
    "arm_conditional"
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, paste0(outcome, ".qs"))
}

build_update_2026_primary_arm_conditional_to_file <- function(
    record,
    data_by_outcome,
    n_draws = 4000L,
    conf_level = 0.95
) {
  data <- data_by_outcome[[record$outcome]]
  table <- record$models
  if (
    is.null(data) || !nrow(data) || nrow(table) != 26L ||
      length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1L ||
      length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1
  ) {
    stop("An updated primary arm-conditional input is invalid.",
         call. = FALSE)
  }
  n_draws <- as.integer(n_draws)
  weights <- table$stacking_weight
  if (any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("Primary arm-conditional stacking weights are invalid.",
         call. = FALSE)
  }
  counts <- allocate_primary_mixture_models(weights, n_draws)
  set.seed(record$seed)
  allocation <- sample(rep(seq_along(weights), times = counts))
  averaged_draws <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = nrow(data)
  )

  for (model_position in seq_len(nrow(table))) {
    output_rows <- which(allocation == model_position)
    if (!length(output_rows)) next
    candidate <- qs::qread(table$fit_path[[model_position]])
    if (
      !identical(candidate$outcome, record$outcome) ||
        as.integer(candidate$model_index) != table$model_index[[model_position]] ||
        !identical(candidate$model_name, table$model_name[[model_position]])
    ) {
      stop("An archived fit does not match its arm-prediction record.",
           call. = FALSE)
    }
    validate_primary_arm_candidate_data(candidate$model$data, data)
    arm_re_formula <- make_primary_arm_re_formula(candidate$model)
    draw_matrix <- brms::posterior_epred(
      candidate$model,
      newdata = candidate$model$data,
      re_formula = arm_re_formula,
      ndraws = length(output_rows)
    ) |>
      as.matrix()
    if (
      ncol(draw_matrix) != nrow(data) &&
        nrow(draw_matrix) == nrow(data)
    ) {
      draw_matrix <- t(draw_matrix)
    }
    if (
      !identical(dim(draw_matrix), c(length(output_rows), nrow(data))) ||
        anyNA(draw_matrix) || any(!is.finite(draw_matrix))
    ) {
      stop(
        "Arm-conditional posterior predictions do not align with the data.",
        call. = FALSE
      )
    }
    averaged_draws[output_rows, ] <- draw_matrix
    rm(candidate, draw_matrix)
    gc(verbose = FALSE)
  }
  if (anyNA(averaged_draws)) {
    stop("The arm-conditional posterior mixture is incomplete.",
         call. = FALSE)
  }

  prediction_id <- make_primary_arm_prediction_ids(data)
  first_profile_row <- match(unique(prediction_id), prediction_id)
  duplicate_rows <- which(duplicated(prediction_id))
  if (length(duplicate_rows)) {
    reference_rows <- first_profile_row[
      match(prediction_id[duplicate_rows], unique(prediction_id))
    ]
    identical_draws <- vapply(seq_along(duplicate_rows), function(index) {
      max(abs(
        averaged_draws[, duplicate_rows[[index]]] -
          averaged_draws[, reference_rows[[index]]]
      )) < 1e-10
    }, logical(1L))
    if (!all(identical_draws)) {
      stop(
        "Observations sharing an arm profile have different predictions.",
        call. = FALSE
      )
    }
  }

  alpha <- (1 - conf_level) / 2
  intervals <- matrixStats::colQuantiles(
    averaged_draws,
    probs = c(alpha, 1 - alpha),
    drop = FALSE
  )
  posterior_sd <- matrixStats::colSds(averaged_draws)
  posterior_precision <- 1 / posterior_sd
  profile_precision <- posterior_precision[first_profile_row]
  precision_range <- range(profile_precision)
  precision_span <- diff(precision_range)
  profile_size <- if (is.finite(precision_span) && precision_span > 0) {
    0.5 + 3 * (profile_precision - precision_range[[1L]]) /
      precision_span
  } else {
    rep(2, length(profile_precision))
  }
  conditional_size <- profile_size[
    match(prediction_id, unique(prediction_id))
  ]

  summary <- tibble::tibble(
    outcome = record$outcome,
    observation = seq_len(nrow(data)),
    prediction_id = prediction_id,
    study_no = data$study_no,
    arm_no = data$arm_no,
    outcome_no = data$outcome_no,
    protein_intake_raw = data$protein_intake_raw,
    resistance_exercise_code = data$resistance_exercise_code,
    duration_centre = data$duration_centre,
    age_centre = data$age_centre,
    estimate = colMeans(averaged_draws),
    conf.low = intervals[, 1L],
    conf.high = intervals[, 2L],
    posterior_sd = posterior_sd,
    posterior_precision = posterior_precision,
    conditional_size = conditional_size
  )
  if (
    any(!is.finite(summary$estimate)) ||
      any(!is.finite(summary$posterior_sd)) ||
      any(summary$posterior_sd <= 0) ||
      any(summary$conf.low > summary$estimate) ||
      any(summary$estimate > summary$conf.high) ||
      any(summary$conditional_size < 0.5 - 1e-10) ||
      any(summary$conditional_size > 3.5 + 1e-10)
  ) {
    stop("Arm-conditional summaries failed validation.", call. = FALSE)
  }

  output <- list(
    outcome = record$outcome,
    predictions = summary,
    draw_allocation = tibble::tibble(
      outcome = record$outcome,
      drawid = seq_len(n_draws),
      model_index = table$model_index[allocation],
      model_name = table$model_name[allocation]
    ),
    settings = list(
      n_draws = n_draws,
      seed = record$seed,
      conf_level = conf_level,
      included_group_effects = c("study", "arm"),
      excluded_group_effects = "outcome"
    )
  )
  output_file <- update_2026_primary_arm_conditional_file_path(record$outcome)
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(output, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise updated arm-conditional predictions.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

audit_update_2026_primary_arm_conditional_files <- function(
    prediction_paths,
    data_by_outcome
) {
  prediction_paths <- unname(unlist(prediction_paths, use.names = FALSE))
  expected_outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  if (
    length(prediction_paths) != length(expected_outcomes) ||
      any(!file.exists(prediction_paths)) ||
      anyDuplicated(normalizePath(prediction_paths, winslash = "/"))
  ) {
    stop("Primary arm-conditional files are incomplete.", call. = FALSE)
  }
  reports <- purrr::map(prediction_paths, qs::qread)
  report_outcomes <- purrr::map_chr(reports, "outcome")
  if (
    !setequal(report_outcomes, expected_outcomes) ||
      anyDuplicated(report_outcomes)
  ) {
    stop("Primary arm-conditional outcomes are invalid.", call. = FALSE)
  }
  reports <- reports[match(expected_outcomes, report_outcomes)]
  predictions <- purrr::map_dfr(reports, "predictions")
  manifest <- purrr::map_dfr(seq_along(reports), function(index) {
    outcome <- expected_outcomes[[index]]
    report <- reports[[index]]
    data <- data_by_outcome[[outcome]]
    outcome_predictions <- report$predictions
    if (
      is.null(data) || nrow(outcome_predictions) != nrow(data) ||
        !identical(
          as.integer(outcome_predictions$observation),
          seq_len(nrow(data))
        )
    ) {
      stop("Primary arm predictions do not align with analysis data.",
           call. = FALSE)
    }
    tibble::tibble(
      outcome = outcome,
      prediction_path = normalizePath(
        prediction_paths[[match(outcome, report_outcomes)]],
        winslash = "/",
        mustWork = TRUE
      ),
      observations = nrow(outcome_predictions),
      arm_profiles = dplyr::n_distinct(outcome_predictions$prediction_id),
      n_draws = as.integer(report$settings$n_draws),
      seed = as.integer(report$settings$seed),
      conf_level = as.numeric(report$settings$conf_level)
    )
  })
  if (
    nrow(predictions) != sum(purrr::map_int(data_by_outcome, nrow)) ||
      any(!is.finite(predictions$estimate)) ||
      any(predictions$posterior_sd <= 0) ||
      any(manifest$n_draws < 1L)
  ) {
    stop("Primary arm-conditional audit failed.", call. = FALSE)
  }
  list(
    predictions = predictions,
    manifest = manifest,
    all_checks_pass = TRUE
  )
}

write_update_2026_primary_arm_conditional_predictions <- function(
    audit,
    output_file = here::here(
      "results", "update_2026", "primary_refit", "reporting",
      "arm_conditional_predictions.csv.gz"
    )
) {
  if (!isTRUE(audit$all_checks_pass)) {
    stop("Cannot write an invalid arm-conditional audit.", call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(audit$predictions, output_file)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

update_2026_primary_reporting_file_path <- function(outcome) {
  directory <- here::here(
    "model_archive", "update_2026", "primary", "reporting", "reports"
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, paste0(outcome, ".qs"))
}

make_update_2026_primary_reporting_file <- function(
  model_average_path,
  data_by_outcome,
  output_file = NULL
) {
  model_average <- qs::qread(model_average_path)
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  dose_response <- summarise_primary_dose_response(model_average)
  contrasts <- compute_primary_heatmap_contrasts(model_average, data)
  local_equivalence <- summarise_primary_local_equivalence(
    contrasts$contrasts,
    outcome
  )
  weight_table <- tibble::tibble(
    outcome = outcome,
    model_index = seq_along(model_average$model_weights),
    model_name = names(model_average$model_weights),
    stacking_weight = as.numeric(model_average$model_weights)
  )
  draw_allocation <- model_average$draw_manifest |>
    dplyr::count(model_index, model_name, name = "draws") |>
    dplyr::mutate(outcome = outcome, .before = 1L)
  report <- list(
    outcome = outcome,
    dose_response = dplyr::mutate(
      dose_response, outcome = outcome, .before = 1L
    ),
    pairwise_contrasts = dplyr::mutate(
      contrasts$contrasts, outcome = outcome, .before = 1L
    ),
    support = dplyr::mutate(
      contrasts$support, outcome = outcome, .before = 1L
    ),
    local_equivalence = local_equivalence,
    model_weights = weight_table,
    draw_allocation = draw_allocation,
    settings = model_average$settings,
    model_average_path = normalizePath(
      model_average_path, winslash = "/", mustWork = TRUE
    )
  )

  if (is.null(output_file)) output_file <- update_2026_primary_reporting_file_path(outcome)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(report, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated primary reporting file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

audit_update_2026_primary_reporting_files <- function(report_paths) {
  report_paths <- unname(unlist(report_paths, use.names = FALSE))
  expected_outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  if (
    length(report_paths) != 4L || any(!file.exists(report_paths)) ||
      anyDuplicated(normalizePath(report_paths, winslash = "/"))
  ) {
    stop("Updated primary reporting files are incomplete.", call. = FALSE)
  }

  reports <- purrr::map(report_paths, qs::qread)
  report_outcomes <- purrr::map_chr(reports, "outcome")
  if (
    !setequal(report_outcomes, expected_outcomes) ||
      anyDuplicated(report_outcomes)
  ) {
    stop("Updated primary reporting outcomes are invalid.", call. = FALSE)
  }
  reports <- reports[match(expected_outcomes, report_outcomes)]

  dose_response <- purrr::map_dfr(reports, "dose_response")
  pairwise_contrasts <- purrr::map_dfr(reports, "pairwise_contrasts")
  local_equivalence <- purrr::map_dfr(reports, "local_equivalence")
  model_weights <- purrr::map_dfr(reports, "model_weights")
  draw_allocation <- purrr::map_dfr(reports, "draw_allocation")
  support <- purrr::map_dfr(reports, "support")
  manifest <- tibble::tibble(
    outcome = expected_outcomes,
    reporting_path = normalizePath(
      report_paths[match(expected_outcomes, report_outcomes)],
      winslash = "/",
      mustWork = TRUE
    ),
    model_average_path = purrr::map_chr(reports, "model_average_path"),
    n_draws = purrr::map_int(reports, ~ as.integer(.x$settings$n_draws)),
    seed = purrr::map_int(reports, ~ as.integer(.x$settings$seed)),
    protein_min = purrr::map_dbl(reports, ~ .x$settings$protein_min),
    protein_max = purrr::map_dbl(reports, ~ .x$settings$protein_max),
    grid_step = purrr::map_dbl(reports, ~ .x$settings$grid_step)
  )

  probability_columns <- c(
    "probability_positive", "probability_negative",
    "probability_meaningful_gain", "probability_negligible",
    "probability_meaningful_loss", "probability_noninferior",
    "probability_meaningful_benefit",
    "probability_practical_equivalence",
    "probability_meaningful_detriment"
  )
  required_contrast_columns <- c(
    "outcome", "protein_intake_current", "protein_intake_proposed",
    "estimate", "conf_low", "conf_high", "supported",
    probability_columns
  )
  missing_contrast_columns <- setdiff(
    required_contrast_columns,
    names(pairwise_contrasts)
  )
  if (length(missing_contrast_columns)) {
    stop(
      "Updated primary contrasts are missing columns: ",
      paste(missing_contrast_columns, collapse = ", "),
      call. = FALSE
    )
  }

  summary <- purrr::map_dfr(expected_outcomes, function(outcome_name) {
    dose <- dose_response |>
      dplyr::filter(.data$outcome == .env$outcome_name)
    contrast <- pairwise_contrasts |>
      dplyr::filter(.data$outcome == .env$outcome_name)
    local <- local_equivalence |>
      dplyr::filter(.data$outcome == .env$outcome_name)
    weights <- model_weights |>
      dplyr::filter(.data$outcome == .env$outcome_name)
    allocation <- draw_allocation |>
      dplyr::filter(.data$outcome == .env$outcome_name)
    observed_support <- support |>
      dplyr::filter(.data$outcome == .env$outcome_name)
    diagonal <- contrast |>
      dplyr::filter(
        abs(
          .data$protein_intake_current -
            .data$protein_intake_proposed
        ) < 1e-10
      )

    probability_matrix <- as.matrix(contrast[, probability_columns])
    partition_sum <-
      contrast$probability_meaningful_gain +
      contrast$probability_negligible +
      contrast$probability_meaningful_loss
    expected_grid_rows <- 2L * (
      round(
        (
          manifest$protein_max[manifest$outcome == outcome_name] -
            manifest$protein_min[manifest$outcome == outcome_name]
        ) /
          manifest$grid_step[manifest$outcome == outcome_name]
      ) + 1L
    )
    tibble::tibble(
      outcome = outcome_name,
      dose_response_rows = nrow(dose),
      pairwise_contrast_rows = nrow(contrast),
      local_equivalence_rows = nrow(local),
      model_count = nrow(weights),
      allocated_models = nrow(allocation),
      allocated_draws = sum(allocation$draws),
      support_rows = nrow(observed_support),
      dose_response_complete =
        nrow(dose) == expected_grid_rows &&
        all(is.finite(dose$estimate)) &&
        all(is.finite(dose$conf.low)) &&
        all(is.finite(dose$conf.high)) &&
        all(dose$conf.low <= dose$estimate) &&
        all(dose$estimate <= dose$conf.high),
      contrast_intervals_valid =
        nrow(contrast) > 0L &&
        all(is.finite(contrast$estimate)) &&
        all(is.finite(contrast$conf_low)) &&
        all(is.finite(contrast$conf_high)) &&
        all(contrast$conf_low <= contrast$estimate) &&
        all(contrast$estimate <= contrast$conf_high),
      contrast_probabilities_valid =
        nrow(probability_matrix) == nrow(contrast) &&
        all(is.finite(probability_matrix)) &&
        all(probability_matrix >= 0) &&
        all(probability_matrix <= 1) &&
        all(abs(partition_sum - 1) < 1e-10),
      contrast_diagonal_valid =
        nrow(diagonal) > 0L &&
        all(abs(diagonal$estimate) < 1e-10) &&
        all(diagonal$probability_negligible == 1),
      weights_valid =
        nrow(weights) == 26L &&
        identical(sort(as.integer(weights$model_index)), seq_len(26L)) &&
        all(is.finite(weights$stacking_weight)) &&
        all(weights$stacking_weight >= 0) &&
        abs(sum(weights$stacking_weight) - 1) < 1e-8,
      allocation_valid =
        nrow(allocation) >= 1L && nrow(allocation) <= 26L &&
        all(allocation$draws > 0L) &&
        sum(allocation$draws) == 4000L &&
        all(allocation$model_index %in% weights$model_index),
      support_valid =
        nrow(observed_support) == 2L &&
        identical(
          sort(as.numeric(observed_support$resistance_exercise_code)),
          c(0, 1)
        ) &&
        all(is.finite(observed_support$observed_min)) &&
        all(is.finite(observed_support$observed_max)) &&
        all(observed_support$observed_min <= observed_support$observed_max) &&
        all(!is.na(contrast$supported)),
      local_equivalence_valid =
        nrow(local) == 2L &&
        all(is.finite(local$equivalent_radius)) &&
        all(is.finite(local$equivalent_radius_q10)) &&
        all(is.finite(local$equivalent_radius_q90)) &&
        all(is.finite(local$median_total_span)) &&
        all(local$supported_current_intakes > 0L)
    )
  }) |>
    dplyr::mutate(
      reporting_checks_pass = dplyr::if_all(
        c(dplyr::ends_with("_valid"), dplyr::ends_with("_complete")),
        identity
      )
    )

  all_reporting_checks_pass <-
    nrow(summary) == 4L && all(summary$reporting_checks_pass)
  if (!isTRUE(all_reporting_checks_pass)) {
    stop("Updated primary reporting audit failed.", call. = FALSE)
  }

  list(
    summary = summary,
    manifest = manifest,
    dose_response_predictions = dose_response,
    pairwise_contrasts = pairwise_contrasts,
    local_equivalence = local_equivalence,
    model_weights = model_weights,
    draw_allocation = draw_allocation,
    support = support,
    all_reporting_checks_pass = all_reporting_checks_pass
  )
}

write_update_2026_primary_reporting_audit <- function(
  audit,
  output_dir = here::here(
    "results", "update_2026", "primary_refit", "reporting"
  )
) {
  if (!isTRUE(audit$all_reporting_checks_pass)) {
    stop("Cannot write a failed updated primary reporting audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    manifest = file.path(output_dir, "manifest.csv"),
    dose_response_predictions = file.path(
      output_dir, "dose_response_predictions.csv.gz"
    ),
    pairwise_contrasts = file.path(
      output_dir, "pairwise_contrasts.csv.gz"
    ),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv"),
    draw_allocation = file.path(output_dir, "draw_allocation.csv"),
    support = file.path(output_dir, "support.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(audit$manifest, paths[["manifest"]])
  readr::write_csv(
    audit$dose_response_predictions,
    paths[["dose_response_predictions"]]
  )
  readr::write_csv(
    audit$pairwise_contrasts,
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    audit$local_equivalence,
    paths[["local_equivalence"]]
  )
  readr::write_csv(audit$model_weights, paths[["model_weights"]])
  readr::write_csv(audit$draw_allocation, paths[["draw_allocation"]])
  readr::write_csv(audit$support, paths[["support"]])
  unname(paths)
}

compile_update_2026_primary_final_diagnostic_audit <- function(
  initial_diagnostics,
  remedy_diagnostics
) {
  if (length(initial_diagnostics) != 104L || length(remedy_diagnostics) != 1L) {
    stop("Updated primary diagnostic inputs are incomplete.", call. = FALSE)
  }
  diagnostic_key <- function(x) {
    paste(x$summary$outcome[[1L]], x$summary$model_index[[1L]], sep = "::")
  }
  initial_keys <- unname(purrr::map_chr(initial_diagnostics, diagnostic_key))
  remedy_keys <- unname(purrr::map_chr(remedy_diagnostics, diagnostic_key))
  if (
    anyDuplicated(initial_keys) || anyDuplicated(remedy_keys) ||
      !identical(remedy_keys, "function::10") ||
      !all(remedy_keys %in% initial_keys)
  ) {
    stop("Updated primary diagnostic replacement is invalid.", call. = FALSE)
  }

  final_diagnostics <- initial_diagnostics
  final_diagnostics[match(remedy_keys, initial_keys)] <- remedy_diagnostics
  split_diagnostics <- split(
    final_diagnostics,
    purrr::map_chr(
      final_diagnostics,
      ~ .x$summary$outcome[[1L]]
    )
  )
  audit <- compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
  remedy_summary <- audit$diagnostic_summary |>
    dplyr::filter(.data$outcome == "function", .data$model_index == 10L)
  if (
    nrow(audit$diagnostic_summary) != 104L ||
      any(!audit$diagnostic_summary$diagnostic_pass) ||
      nrow(remedy_summary) != 1L ||
      remedy_summary$adapt_delta[[1L]] != 0.9995 ||
      remedy_summary$divergent_transitions[[1L]] != 0L
  ) {
    stop("The final updated primary diagnostic audit did not pass.",
         call. = FALSE)
  }
  audit
}

make_update_2026_primary_final_specification_table <- function(
  initial_specifications,
  remedy_specifications,
  accepted_fit_manifest
) {
  initial <- dplyr::bind_rows(initial_specifications)
  remedy <- dplyr::bind_rows(remedy_specifications)
  keys <- c("outcome", "model_index")
  if (
    nrow(initial) != 104L || nrow(remedy) != 1L ||
      remedy$outcome[[1L]] != "function" ||
      remedy$model_index[[1L]] != 10L ||
      remedy$adapt_delta[[1L]] != 0.9995
  ) {
    stop("Updated primary specification replacement is invalid.",
         call. = FALSE)
  }

  final <- dplyr::bind_rows(
    dplyr::anti_join(initial, remedy[keys], by = keys),
    remedy
  ) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::rename(
          accepted_fit_path = fit_path,
          accepted_fit_stage = stage
        ),
      by = keys
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  pair <- paste(final$outcome, final$model_index)
  remedy_row <- final |>
    dplyr::filter(.data$outcome == "function", .data$model_index == 10L)
  if (
    nrow(final) != 104L || anyDuplicated(pair) ||
      anyNA(final$accepted_fit_path) ||
      any(!file.exists(final$accepted_fit_path)) ||
      nrow(remedy_row) != 1L ||
      remedy_row$accepted_fit_stage[[1L]] != "remedy_1" ||
      remedy_row$adapt_delta[[1L]] != 0.9995
  ) {
    stop("The final updated primary specification table is invalid.",
         call. = FALSE)
  }
  final
}

make_update_2026_primary_source_manifest <- function(
  canonical_workbook,
  update_workbook,
  review_decisions
) {
  paths <- c(
    canonical_workbook = canonical_workbook,
    update_workbook = update_workbook,
    review_decisions = review_decisions
  )
  if (any(!file.exists(paths))) {
    stop("An updated primary source file is missing.", call. = FALSE)
  }
  tibble::tibble(
    source = names(paths),
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(paths)$size),
    sha256 = purrr::map_chr(
      paths,
      ~ digest::digest(file = .x, algo = "sha256")
    )
  )
}

audit_update_2026_primary_verification_checkpoint <- function(
  data_by_outcome,
  source_manifest,
  review_decisions,
  preparation_audit,
  final_specifications,
  final_diagnostic_audit,
  accepted_fit_manifest,
  final_loo_audit,
  reporting_audit,
  reference_values
) {
  expected_outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  registry_keys <- function(data) {
    data |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::mutate(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(final_specifications)
  diagnostic_keys <- registry_keys(
    final_diagnostic_audit$diagnostic_summary
  )
  accepted_fit_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  reporting_keys <- registry_keys(reporting_audit$model_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_fit_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, reporting_keys)

  name_alignment <- final_specifications |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      final_diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      reporting_audit$model_weights |>
        dplyr::select(outcome, model_index, reporting_name = model_name),
      by = c("outcome", "model_index")
    )

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "all four combined analysis datasets are present",
      "combined analysis rows are complete for model inputs",
      "the final model registries contain the same 104 keys",
      "model names agree across specifications diagnostics LOO and reporting",
      "all accepted fit paths exist",
      "all final model diagnostics pass",
      "final LOO is corrected and stacking-ready",
      "all model-averaged reporting checks pass",
      "all four stacking-weight vectors contain 26 models and sum to one",
      "the combined protein-intake centre matches the reference value",
      "the combined-data preparation audit remains valid"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      identical(names(data_by_outcome), expected_outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L),
      nrow(analysis_data) > 0L &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) &&
        all(analysis_data$vi > 0) &&
        all(is.finite(analysis_data$protein_intake_raw)),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$specification_name == name_alignment$diagnostic_name) &&
        all(name_alignment$specification_name == name_alignment$loo_name) &&
        all(name_alignment$specification_name == name_alignment$reporting_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(file.exists(accepted_fit_manifest$fit_path)),
      nrow(final_diagnostic_audit$diagnostic_summary) == 104L &&
      all(final_diagnostic_audit$diagnostic_summary$diagnostic_pass),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$remaining_high_k) == 0L &&
        nrow(final_loo_audit$model_summary) == 104L,
      isTRUE(reporting_audit$all_reporting_checks_pass),
      nrow(reporting_audit$model_weights) == 104L &&
        all(
          abs(
            reporting_audit$model_weights |>
              dplyr::summarise(
                weight_sum = sum(stacking_weight),
                .by = outcome
              ) |>
              dplyr::pull(weight_sum) - 1
          ) < 1e-8
        ),
      length(unique(analysis_data$protein_median)) == 1L &&
        isTRUE(all.equal(
          unique(analysis_data$protein_median)[[1L]],
          reference_values$habitual_protein_intake
        )),
      all(preparation_audit$checks$passed)
    )
  )
  failed_checks <- is.na(checks$passed) | !checks$passed
  if (any(failed_checks)) {
    stop(
      "Updated primary verification failed: ",
      paste(checks$check[failed_checks], collapse = "; "),
      call. = FALSE
    )
  }

  list(
    checks = checks,
    analysis_data = analysis_data,
    source_manifest = source_manifest,
    review_decisions = review_decisions,
    preparation_checks = preparation_audit$checks,
    preparation_outcome_counts = preparation_audit$outcome_counts,
    preparation_update_study_counts = preparation_audit$update_study_counts,
    model_name_alignment = name_alignment
  )
}

update_2026_primary_figure_titles <- function() {
  c(
    whole_body = "change in whole-body lean/fat-free mass",
    hypertrophy = "change in measures of muscle size or mass",
    strength = "strength change",
    `function` = "function change"
  )
}

make_update_2026_primary_figure <- function(
  outcome,
  reporting_audit,
  data_by_outcome
) {
  titles <- update_2026_primary_figure_titles()
  if (!outcome %in% names(titles)) {
    stop("Unknown updated primary figure outcome.", call. = FALSE)
  }
  predictions <- reporting_audit$dose_response_predictions |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::select(-outcome)
  contrasts <- reporting_audit$pairwise_contrasts |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::select(-outcome)
  support <- reporting_audit$support |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::select(-outcome)
  weights <- reporting_audit$model_weights |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::arrange(model_index)
  settings <- reporting_audit$manifest |>
    dplyr::filter(.data$outcome == .env$outcome)
  data <- data_by_outcome[[outcome]]
  if (
    !nrow(predictions) || !nrow(contrasts) || nrow(support) != 2L ||
      nrow(weights) != 26L || nrow(settings) != 1L || is.null(data) ||
      !nrow(data)
  ) {
    stop("Updated primary figure inputs are incomplete.", call. = FALSE)
  }

  heatmap <- list(
    contrasts = contrasts,
    support = support,
    model_weights = stats::setNames(
      weights$stacking_weight,
      weights$model_name
    ),
    settings = list(
      protein_min = settings$protein_min[[1L]],
      protein_max = settings$protein_max[[1L]],
      grid_step = settings$grid_step[[1L]],
      n_draws = settings$n_draws[[1L]],
      sesoi = 0.10,
      conf_level = 0.95
    )
  )
  class(heatmap) <- c("primary_heatmap_contrasts", class(heatmap))
  plot_primary_predictions_heatmap(
    data = data,
    predictions = predictions,
    heatmap_contrasts = heatmap,
    title = titles[[outcome]]
  )
}

make_update_2026_primary_arm_figure <- function(
    outcome,
    reporting_audit,
    arm_conditional_audit,
    data_by_outcome,
    include_observed_effects = TRUE
) {
  titles <- update_2026_primary_figure_titles()
  if (
    !outcome %in% names(titles) ||
      length(include_observed_effects) != 1L ||
      !is.logical(include_observed_effects) ||
      is.na(include_observed_effects) ||
      !isTRUE(arm_conditional_audit$all_checks_pass)
  ) {
    stop("Unknown or invalid arm-conditional figure outcome.",
         call. = FALSE)
  }
  predictions <- reporting_audit$dose_response_predictions |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::select(-outcome)
  arm_predictions <- arm_conditional_audit$predictions |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::select(-outcome)
  contrasts <- reporting_audit$pairwise_contrasts |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::select(-outcome)
  support <- reporting_audit$support |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::select(-outcome)
  weights <- reporting_audit$model_weights |>
    dplyr::filter(.data$outcome == .env$outcome) |>
    dplyr::arrange(model_index)
  settings <- reporting_audit$manifest |>
    dplyr::filter(.data$outcome == .env$outcome)
  data <- data_by_outcome[[outcome]]
  if (
    !nrow(predictions) || !nrow(contrasts) || nrow(support) != 2L ||
      nrow(weights) != 26L || nrow(settings) != 1L || is.null(data) ||
      !nrow(data) || nrow(arm_predictions) != nrow(data)
  ) {
    stop("Updated primary arm-figure inputs are incomplete.",
         call. = FALSE)
  }

  heatmap <- list(
    contrasts = contrasts,
    support = support,
    model_weights = stats::setNames(
      weights$stacking_weight,
      weights$model_name
    ),
    settings = list(
      protein_min = settings$protein_min[[1L]],
      protein_max = settings$protein_max[[1L]],
      grid_step = settings$grid_step[[1L]],
      n_draws = settings$n_draws[[1L]],
      sesoi = 0.10,
      conf_level = 0.95
    )
  )
  class(heatmap) <- c("primary_heatmap_contrasts", class(heatmap))
  plot_primary_arm_predictions_heatmap(
    data = data,
    predictions = predictions,
    arm_predictions = arm_predictions,
    heatmap_contrasts = heatmap,
    title = titles[[outcome]],
    include_observed_effects = include_observed_effects
  )
}

write_update_2026_primary_figures <- function(
  reporting_audit,
  data_by_outcome,
  output_dir = here::here("plots", "update_2026", "primary")
) {
  if (!isTRUE(reporting_audit$all_reporting_checks_pass)) {
    stop("Cannot plot an unverified updated primary report.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  figure_paths <- purrr::map(outcomes, function(outcome_name) {
    figure <- make_update_2026_primary_figure(
      outcome = outcome_name,
      reporting_audit = reporting_audit,
      data_by_outcome = data_by_outcome
    )
    tiff_path <- file.path(output_dir, paste0(outcome_name, ".tiff"))
    pdf_path <- file.path(
      output_dir,
      paste0(outcome_name, "-manuscript.pdf")
    )
    ggplot2::ggsave(
      tiff_path,
      figure,
      device = "tiff",
      dpi = 300,
      width = 8,
      height = 12
    )
    ggplot2::ggsave(
      pdf_path,
      figure,
      device = grDevices::cairo_pdf,
      width = 8,
      height = 12
    )
    c(tiff_path, pdf_path)
  }) |>
    unlist(use.names = FALSE)
  if (
    length(figure_paths) != 8L || any(!file.exists(figure_paths)) ||
      any(file.info(figure_paths)$size <= 0)
  ) {
    stop("Updated primary figures were not written correctly.",
         call. = FALSE)
  }

  manifest <- tibble::tibble(
    file = basename(figure_paths),
    path = normalizePath(
      figure_paths,
      winslash = "/",
      mustWork = TRUE
    ),
    bytes = as.numeric(file.info(figure_paths)$size),
    sha256 = purrr::map_chr(
      figure_paths,
      ~ digest::digest(file = .x, algo = "sha256")
    )
  )
  manifest_path <- file.path(output_dir, "manifest.csv")
  readr::write_csv(manifest, manifest_path)
  normalizePath(
    c(figure_paths, manifest_path),
    winslash = "/",
    mustWork = TRUE
  )
}

write_update_2026_primary_arm_figures <- function(
    reporting_audit,
    arm_conditional_audit,
    data_by_outcome,
    include_observed_effects = TRUE,
    output_dir = here::here(
      "plots", "update_2026", "primary_arm_conditional"
    )
) {
  if (
    !isTRUE(reporting_audit$all_reporting_checks_pass) ||
      !isTRUE(arm_conditional_audit$all_checks_pass)
  ) {
    stop("Cannot plot unverified arm-conditional primary reports.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  figure_paths <- purrr::map(outcomes, function(outcome_name) {
    figure <- make_update_2026_primary_arm_figure(
      outcome = outcome_name,
      reporting_audit = reporting_audit,
      arm_conditional_audit = arm_conditional_audit,
      data_by_outcome = data_by_outcome,
      include_observed_effects = include_observed_effects
    )
    tiff_path <- file.path(output_dir, paste0(outcome_name, ".tiff"))
    pdf_path <- file.path(
      output_dir,
      paste0(outcome_name, "-manuscript.pdf")
    )
    png_path <- file.path(output_dir, paste0(outcome_name, ".png"))
    ggplot2::ggsave(
      tiff_path,
      figure,
      device = "tiff",
      dpi = 300,
      width = 8,
      height = 12
    )
    ggplot2::ggsave(
      pdf_path,
      figure,
      device = grDevices::cairo_pdf,
      width = 8,
      height = 12
    )
    ggplot2::ggsave(
      png_path,
      figure,
      device = "png",
      dpi = 180,
      width = 8,
      height = 12
    )
    c(tiff_path, pdf_path, png_path)
  }) |>
    unlist(use.names = FALSE)
  if (
    length(figure_paths) != 12L || any(!file.exists(figure_paths)) ||
      any(file.info(figure_paths)$size <= 0)
  ) {
    stop("Arm-conditional primary figures were not written correctly.",
         call. = FALSE)
  }

  manifest <- tibble::tibble(
    file = basename(figure_paths),
    path = normalizePath(
      figure_paths,
      winslash = "/",
      mustWork = TRUE
    ),
    bytes = as.numeric(file.info(figure_paths)$size),
    sha256 = purrr::map_chr(
      figure_paths,
      ~ digest::digest(file = .x, algo = "sha256")
    )
  )
  manifest_path <- file.path(output_dir, "manifest.csv")
  readr::write_csv(manifest, manifest_path)
  normalizePath(
    c(figure_paths, manifest_path),
    winslash = "/",
    mustWork = TRUE
  )
}

write_update_2026_primary_arm_only_figures <- function(
    reporting_audit,
    arm_conditional_audit,
    data_by_outcome,
    output_dir = here::here(
      "plots", "update_2026", "primary_arm_conditional_only"
    )
) {
  write_update_2026_primary_arm_figures(
    reporting_audit = reporting_audit,
    arm_conditional_audit = arm_conditional_audit,
    data_by_outcome = data_by_outcome,
    include_observed_effects = FALSE,
    output_dir = output_dir
  )
}

# Only the static figures used by the two documents belong in this directory.
write_update_2026_primary_publication_figures <- function(
    arm_only_files,
    arm_observed_files,
    output_dir = here::here(
      "manuscript", "generated", "figures", "primary-dose-response"
    )
) {
  outcomes <- names(update_2026_primary_figure_titles())
  figure_names <- paste0(outcomes, ".png")
  sources <- c(
    arm_only_files[match(figure_names, basename(arm_only_files))],
    arm_observed_files[match(figure_names, basename(arm_observed_files))]
  )
  if (anyNA(sources) || any(!file.exists(sources)) ||
      any(file.info(sources)$size <= 0)) {
    stop("A primary dose-response publication figure is missing.",
         call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  destinations <- file.path(
    output_dir,
    c(paste0(outcomes, "-arm-conditional.png"), paste0(outcomes, "-observed.png"))
  )
  temporary <- tempfile(tmpdir = output_dir)
  on.exit(unlink(temporary), add = TRUE)
  for (index in seq_along(sources)) {
    if (!file.copy(sources[[index]], temporary)) {
      stop("Could not copy figure: ", sources[[index]], call. = FALSE)
    }
    fs::file_move(temporary, destinations[[index]])
  }
  destinations
}

verify_update_2026_manifest_files <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop("Required freeze manifest is missing: ", manifest_path,
         call. = FALSE)
  }
  manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
  required <- c("file", "sha256")
  if (!all(required %in% names(manifest)) || !nrow(manifest)) {
    stop("Required freeze manifest is invalid: ", manifest_path,
         call. = FALSE)
  }
  paths <- file.path(dirname(manifest_path), manifest$file)
  if (any(!file.exists(paths))) {
    stop("A file named by the freeze manifest is missing.", call. = FALSE)
  }
  observed <- purrr::map_chr(
    paths,
    ~ digest::digest(file = .x, algo = "sha256")
  )
  if (any(observed != manifest$sha256)) {
    stop("A freeze-manifest checksum does not match.", call. = FALSE)
  }
  invisible(TRUE)
}

make_update_2026_primary_pruning_plan <- function(
  accepted_fit_manifest,
  verification_files,
  figure_files,
  archive_root = here::here("model_archive", "update_2026", "primary")
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "primary_update_2026"),
    "manifest.csv"
  )
  figure_manifest <- file.path(
    here::here("plots", "update_2026", "primary"),
    "manifest.csv"
  )
  if (
    any(!file.exists(verification_files)) ||
      any(!file.exists(figure_files))
  ) {
    stop("Updated primary freeze outputs are incomplete.", call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)
  verify_update_2026_manifest_files(figure_manifest)

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here("model_archive", "update_2026", "primary"),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated primary archive.",
         call. = FALSE)
  }

  fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  fit_paths <- unlist(
    purrr::map(
      fit_roots,
      ~ list.files(.x, pattern = "\\.qs$", recursive = TRUE,
                   full.names = TRUE)
    ),
    use.names = FALSE
  )
  fit_paths <- normalizePath(fit_paths, winslash = "/", mustWork = TRUE)
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(fit_paths) != 105L || length(accepted_paths) != 104L ||
      anyDuplicated(fit_paths) || anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% fit_paths) ||
      sum(!fit_paths %in% accepted_paths) != 1L
  ) {
    stop("Updated primary fit retention does not match the accepted manifest.",
         call. = FALSE)
  }

  fit_manifest <- tibble::tibble(
    path = fit_paths,
    relative_path = substring(
      fit_paths,
      nchar(resolved_root) + 2L
    ),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", fit_paths, fixed = TRUE) ~ "initial",
      grepl("/remedies/", fit_paths, fixed = TRUE) ~ "remedy",
      TRUE ~ "unknown"
    ),
    accepted_fit = fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(fit_paths)$size),
    conservative_action = ifelse(
      fit_paths %in% accepted_paths,
      "retain",
      "delete superseded fit"
    ),
    lean_action = "delete frozen fit archive"
  ) |>
    dplyr::arrange(accepted_fit, archive_stage, relative_path)

  compact_roots <- file.path(archive_root, c("loo", "reporting"))
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(.x, recursive = TRUE, full.names = TRUE)
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[file.exists(compact_paths)]
  compact_manifest <- tibble::tibble(
    path = normalizePath(
      compact_paths,
      winslash = "/",
      mustWork = TRUE
    ),
    bytes = as.numeric(file.info(compact_paths)$size),
    action = "retain compact LOO/reporting archive"
  )

  summary <- tibble::tibble(
    mode = c("conservative", "lean"),
    fit_files_deleted = c(
      sum(!fit_manifest$accepted_fit),
      nrow(fit_manifest)
    ),
    fit_files_retained = c(
      sum(fit_manifest$accepted_fit),
      0L
    ),
    gib_reclaimed = c(
      sum(fit_manifest$bytes[!fit_manifest$accepted_fit]),
      sum(fit_manifest$bytes)
    ) / 1024^3,
    compact_archive_gib_retained =
      sum(compact_manifest$bytes) / 1024^3,
    restoration = c(
      "All accepted fits remain locally reusable.",
      "Reported results remain frozen; restoring fits requires refitting."
    )
  )

  list(
    summary = summary,
    fit_file_manifest = fit_manifest,
    compact_archive_manifest = compact_manifest,
    verification_manifest = normalizePath(
      verification_manifest, winslash = "/", mustWork = TRUE
    ),
    figure_manifest = normalizePath(
      figure_manifest, winslash = "/", mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_primary_pruning_plan <- function(
  plan,
  output_dir = here::here(
    "results", "update_2026", "primary_refit", "pruning"
  )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    fit_file_manifest = file.path(output_dir, "fit_file_manifest.csv"),
    compact_archive_manifest = file.path(
      output_dir, "compact_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(
    plan$fit_file_manifest,
    paths[["fit_file_manifest"]]
  )
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_primary_pruning <- function(
  plan,
  mode = c("conservative", "lean"),
  confirmation,
  receipt_dir = here::here(
    "results", "update_2026", "primary_refit", "pruning"
  )
) {
  mode <- match.arg(mode)
  required_confirmation <- c(
    conservative = "DELETE SUPERSEDED UPDATED PRIMARY FIT",
    lean = "DELETE ALL UPDATED PRIMARY FIT ARCHIVES"
  )[[mode]]
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)
  verify_update_2026_manifest_files(plan$figure_manifest)

  expected_root <- normalizePath(
    here::here("model_archive", "update_2026", "primary"),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the expected archive root.",
         call. = FALSE)
  }

  selected <- if (identical(mode, "conservative")) {
    plan$fit_file_manifest |>
      dplyr::filter(.data$conservative_action != "retain")
  } else {
    plan$fit_file_manifest
  }
  if (
    (identical(mode, "conservative") &&
       (nrow(selected) != 1L || any(selected$accepted_fit))) ||
      (identical(mode, "lean") && nrow(selected) != 105L)
  ) {
    stop("The selected pruning set is inconsistent with its policy.",
         call. = FALSE)
  }

  receipt_path <- file.path(
    receipt_dir,
    paste0("executed_", mode, ".csv")
  )
  if (file.exists(receipt_path)) {
    stop("This pruning policy already has an execution receipt.",
         call. = FALSE)
  }
  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  if (
    any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the permitted archive subtree.",
         call. = FALSE)
  }
  if (
    identical(mode, "conservative") &&
      !grepl(
        "/candidates/function/model_10\\.qs$",
        selected_paths[[1L]]
      )
  ) {
    stop("The conservative deletion is not the expected superseded fit.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256")
      ),
      pruning_mode = mode,
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop("At least one selected fit archive could not be deleted.",
         call. = FALSE)
  }
  if (
    identical(mode, "conservative") &&
      any(!file.exists(
        plan$fit_file_manifest$path[plan$fit_file_manifest$accepted_fit]
      ))
  ) {
    stop("An accepted updated primary fit is unexpectedly missing.",
         call. = FALSE)
  }

  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, receipt_path)
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

read_frozen_primary_analysis <- function(path) {
  data <- readr::read_csv(
    path,
    show_col_types = FALSE,
    guess_max = 100000
  )

  required <- c(
    "analysis_outcome", "study_no", "arm_no", "outcome_no",
    "protein_median", "protein_intake_raw", "protein_intake", "r",
    "yi", "vi", "outcome_group"
  )
  missing <- setdiff(required, names(data))

  if (length(missing)) {
    stop(
      "Frozen primary analysis is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  data |>
    dplyr::mutate(
      resistance_exercise_factor = factor(
        resistance_exercise_yes_no,
        levels = c("yes", "no")
      ),
      study_factor = factor(study_no)
    )
}

derive_frozen_preprocessing_reference <- function(frozen_analysis) {
  protein_medians <- unique(stats::na.omit(frozen_analysis$protein_median))

  if (length(protein_medians) != 1L) {
    stop(
      "The frozen primary analysis does not contain one protein-intake centre.",
      call. = FALSE
    )
  }

  r_frequency <- frozen_analysis |>
    dplyr::filter(!is.na(r)) |>
    dplyr::count(r, sort = TRUE)

  if (!nrow(r_frequency) || r_frequency$n[[1]] <= nrow(frozen_analysis) / 2) {
    stop(
      "Could not identify the frozen pre-post correlation imputation value.",
      call. = FALSE
    )
  }

  list(
    protein_intake_centre = as.numeric(protein_medians[[1]]),
    pre_post_r_imputation = as.numeric(r_frequency$r[[1]]),
    age_centre_years = 25,
    duration_centre_weeks = 12,
    max_study_no = max(frozen_analysis$study_no, na.rm = TRUE),
    max_arm_no = max(frozen_analysis$arm_no, na.rm = TRUE),
    max_outcome_no = max(frozen_analysis$outcome_no, na.rm = TRUE),
    frozen_rows = nrow(frozen_analysis),
    frozen_studies = dplyr::n_distinct(frozen_analysis$study_no)
  )
}

read_update_review_decisions <- function(path) {
  decisions <- readr::read_csv(
    path,
    col_types = readr::cols(
      study_no = readr::col_double(),
      authors = readr::col_character(),
      quantitative_status = readr::col_character(),
      include_quantitative = readr::col_logical(),
      reason = readr::col_character(),
      decision_source = readr::col_character()
    )
  )

  if (anyDuplicated(decisions$study_no)) {
    stop("Update review decisions contain duplicate study numbers.", call. = FALSE)
  }

  if (any(is.na(decisions$include_quantitative))) {
    stop("Every update study needs an explicit quantitative decision.", call. = FALSE)
  }

  decisions
}

validate_update_workbook <- function(
  update_data,
  frozen_analysis,
  decisions
) {
  required_sheets <- c(
    "study_characteristics",
    "arm_characteristics",
    "arm_based_data_by_timepoint",
    "contrasts_between_timepoints"
  )
  missing_sheets <- setdiff(required_sheets, names(update_data))

  if (length(missing_sheets)) {
    stop(
      "Update workbook is missing worksheets: ",
      paste(missing_sheets, collapse = ", "),
      call. = FALSE
    )
  }

  studies <- update_data$study_characteristics
  update_studies <- sort(unique(studies$study_no))
  frozen_studies <- sort(unique(frozen_analysis$study_no))
  overlapping <- intersect(update_studies, frozen_studies)

  if (length(overlapping)) {
    stop(
      "Update study numbers overlap the frozen analysis: ",
      paste(overlapping, collapse = ", "),
      call. = FALSE
    )
  }

  missing_decisions <- setdiff(update_studies, decisions$study_no)
  extra_decisions <- setdiff(decisions$study_no, update_studies)

  if (length(missing_decisions) || length(extra_decisions)) {
    stop(
      "Review-decision study numbers do not exactly match the update workbook.",
      call. = FALSE
    )
  }

  arm_keys <- update_data$arm_characteristics |>
    dplyr::distinct(study_no, arm_no)
  observation_orphans <- update_data$arm_based_data_by_timepoint |>
    dplyr::anti_join(arm_keys, by = c("study_no", "arm_no"))
  contrast_orphans <- update_data$contrasts_between_timepoints |>
    dplyr::anti_join(arm_keys, by = c("study_no", "arm_no"))

  if (nrow(observation_orphans) || nrow(contrast_orphans)) {
    stop("Update workbook contains outcome rows without a matching arm.", call. = FALSE)
  }

  tibble::tibble(
    check = c(
      "required worksheets present",
      "study numbers disjoint from frozen analysis",
      "one review decision per update study",
      "outcome rows have matching arms"
    ),
    passed = TRUE
  )
}

rebuild_primary_eligibility_flags <- function(data) {
  exclusion_labels <- primary_exclusion_labels()

  data |>
    dplyr::mutate(
      exclusion_missing_mean_change = is.na(mean_change),
      exclusion_missing_baseline_sd = is.na(baseline_sd_pooled),
      exclusion_missing_sample_size = is.na(n),
      exclusion_missing_pre_post_correlation = is.na(r),
      exclusion_effect_size_not_estimable =
        is.na(yi) | is.na(vi) | dplyr::coalesce(vi <= 0, FALSE),
      exclusion_missing_protein_intake = is.na(protein_intake),
      exclusion_missing_resistance_exercise_code =
        is.na(resistance_exercise_code),
      exclusion_missing_mean_age = is.na(age_centre),
      exclusion_missing_duration = is.na(duration_centre)
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      analysis_exclusion_count = sum(
        dplyr::c_across(dplyr::all_of(names(exclusion_labels)))
      ),
      analysis_exclusion_reasons = {
        flags <- dplyr::c_across(dplyr::all_of(names(exclusion_labels)))
        reasons <- unname(exclusion_labels[flags])
        if (length(reasons)) paste(reasons, collapse = "; ") else NA_character_
      },
      analysis_eligible_primary = analysis_exclusion_count == 0
    ) |>
    dplyr::ungroup()
}

prepare_incremental_protein_data <- function(
  update_data,
  preprocessing_reference
) {
  prepared <- withCallingHandlers(
    prepare_protein_data(update_data),
    warning = function(warning_condition) {
      if (grepl(
        "rows with NAs omitted from model fitting",
        conditionMessage(warning_condition),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )

  observed_r <- with(
    prepared,
    ifelse(
      !is.na(baseline_sd) & !is.na(followup_sd) &
        !is.na(standard_deviation),
      (baseline_sd^2 + followup_sd^2 - standard_deviation^2) /
        (2 * baseline_sd * followup_sd),
      NA_real_
    )
  )
  observed_r <- pmin(0.999, pmax(-0.999, observed_r))

  prepared <- prepared |>
    dplyr::mutate(
      update_2026_r_before_freeze_alignment = as.numeric(r),
      update_2026_r_source = dplyr::if_else(
        is.na(observed_r),
        "frozen canonical imputation",
        "update study observation"
      ),
      r = dplyr::coalesce(
        observed_r,
        preprocessing_reference$pre_post_r_imputation
      ),
      protein_median = preprocessing_reference$protein_intake_centre,
      protein_intake = protein_intake_raw - protein_median,
      duration_centre = time2 - preprocessing_reference$duration_centre_weeks,
      age_centre = mean_age_years - preprocessing_reference$age_centre_years,
      arm_no = arm_no + preprocessing_reference$max_arm_no,
      outcome_no = outcome_no + preprocessing_reference$max_outcome_no
    ) |>
    dplyr::select(-yi, -vi)

  prepared <- metafor::escalc(
    measure = "SMCR",
    m1i = mean_change,
    m2i = zero,
    sd1i = baseline_sd_pooled,
    ni = n,
    ri = r,
    data = prepared
  )

  rebuild_primary_eligibility_flags(prepared)
}

apply_update_review_decisions <- function(prepared, decisions) {
  joined <- prepared |>
    dplyr::left_join(
      decisions |>
        dplyr::select(
          study_no,
          quantitative_status,
          include_quantitative,
          review_decision_reason = reason,
          review_decision_source = decision_source
        ),
      by = "study_no"
    )

  if (any(is.na(joined$include_quantitative))) {
    stop("Prepared update rows are missing a review decision.", call. = FALSE)
  }

  joined |>
    dplyr::mutate(
      update_2026_analysis_ready =
        analysis_eligible_primary & include_quantitative
    )
}

apply_function_direction <- function(data) {
  directed <- data |>
    dplyr::mutate(
      sign = dplyr::case_when(
        unit_of_measurement == "score" ~ "increase",
        unit_of_measurement == "seconds" ~ "decrease",
        unit_of_measurement == "m/s" ~ "increase",
        unit_of_measurement == "watts" ~ "increase",
        unit_of_measurement == "number" ~ "increase",
        unit_of_measurement == "metres" ~ "increase",
        unit_of_measurement == "cm" ~ "increase",
        unit_of_measurement == "reps" ~ "increase",
        unit_of_measurement == "duration" ~ "increase",
        unit_of_measurement == "minutes" ~ "decrease"
      )
    )

  unknown_units <- directed |>
    dplyr::filter(is.na(sign)) |>
    dplyr::distinct(unit_of_measurement) |>
    dplyr::pull(unit_of_measurement)

  if (length(unknown_units)) {
    stop(
      "Function outcome units lack a direction rule: ",
      paste(unknown_units, collapse = ", "),
      call. = FALSE
    )
  }

  directed |>
    dplyr::mutate(yi = yi * dplyr::if_else(sign == "increase", 1, -1))
}

add_primary_display_fields <- function(data) {
  weight_range <- range(1 / sqrt(data$vi), na.rm = TRUE)
  weight_span <- diff(weight_range)

  data |>
    dplyr::mutate(
      study_factor = factor(study_no),
      wi = 1 / sqrt(vi),
      size = if (is.finite(weight_span) && weight_span > 0) {
        0.5 + 3 * (wi - weight_range[[1]]) / weight_span
      } else {
        2
      }
    )
}

make_update_primary_data_by_outcome <- function(reviewed_data) {
  outcome_names <- update_2026_outcome_names()
  analysis_data <- reviewed_data |>
    dplyr::filter(update_2026_analysis_ready)

  result <- purrr::imap(
    outcome_names,
    function(outcome_group_value, outcome_name) {
      outcome_data <- analysis_data |>
        dplyr::filter(outcome_group == outcome_group_value)

      if (identical(outcome_name, "function")) {
        outcome_data <- apply_function_direction(outcome_data)
      }

      add_primary_display_fields(outcome_data) |>
        dplyr::mutate(
          analysis_outcome = dplyr::recode(
            outcome_name,
            hypertrophy = "direct_muscle",
            `function` = "physical_function",
            .default = outcome_name
          ),
          update_2026_data_source = "new_update"
        )
    }
  )

  names(result) <- names(outcome_names)
  result
}

make_frozen_primary_data_by_outcome <- function(frozen_analysis) {
  labels <- c(
    whole_body = "whole_body",
    hypertrophy = "direct_muscle",
    strength = "strength",
    `function` = "physical_function"
  )

  purrr::imap(
    labels,
    ~ frozen_analysis |>
      dplyr::filter(analysis_outcome == .x) |>
      dplyr::mutate(update_2026_data_source = "frozen_primary")
  )
}

combine_frozen_and_update_primary_data <- function(
  frozen_by_outcome,
  update_by_outcome
) {
  harmonise_for_bind <- function(frozen, update) {
    common <- intersect(names(frozen), names(update))

    for (column in common) {
      frozen_class <- class(frozen[[column]])[[1]]
      update_class <- class(update[[column]])[[1]]

      if (identical(frozen_class, update_class)) next

      if (
        is.character(frozen[[column]]) ||
          is.character(update[[column]]) ||
          is.factor(frozen[[column]]) ||
          is.factor(update[[column]])
      ) {
        frozen[[column]] <- as.character(frozen[[column]])
        update[[column]] <- as.character(update[[column]])
      } else if (
        is.numeric(frozen[[column]]) || is.numeric(update[[column]])
      ) {
        frozen[[column]] <- as.numeric(frozen[[column]])
        update[[column]] <- as.numeric(update[[column]])
      }
    }

    list(frozen = frozen, update = update)
  }

  purrr::map2(
    frozen_by_outcome,
    update_by_outcome,
    function(frozen, update) {
      aligned <- harmonise_for_bind(frozen, update)

      dplyr::bind_rows(aligned$frozen, aligned$update) |>
        dplyr::mutate(
          resistance_exercise_factor = factor(
            resistance_exercise_yes_no,
            levels = c("yes", "no")
          )
        ) |>
        add_primary_display_fields()
    }
  )
}

make_update_fit_availability <- function(sensitivity_retention) {
  primary <- tibble::tibble(
    family = "primary",
    accepted_fits = 104L,
    full_fit_objects_present = 104L,
    fit_location = "model_archive/primary",
    posterior_update_ready = TRUE,
    next_action = "Implement and validate incremental posterior updating."
  )

  sensitivity <- sensitivity_retention |>
    dplyr::transmute(
      family,
      accepted_fits = as.integer(accepted_fits),
      full_fit_objects_present = as.integer(fit_archives_present),
      fit_location = dplyr::if_else(
        family == "rescor" & fit_archives_present == accepted_fits,
        "_targets/user/sensitivity_fit_files/rescor",
        "not retained locally"
      ),
      posterior_update_ready = fit_archives_present == accepted_fits,
      next_action = dplyr::if_else(
        posterior_update_ready,
        "Archive the retained fits, then implement incremental updating.",
        "Reconstruct the accepted fits before posterior updating."
      )
    )

  dplyr::bind_rows(primary, sensitivity)
}

validate_prepared_update <- function(
  reviewed_data,
  new_by_outcome,
  combined_by_outcome,
  preprocessing_reference,
  frozen_analysis
) {
  key_columns <- c(
    "study_no", "arm_no_extracted", "outcome_group",
    "measurement_method", "unit_of_measurement", "time1", "time2"
  )
  duplicate_keys <- reviewed_data |>
    dplyr::count(dplyr::across(dplyr::all_of(key_columns))) |>
    dplyr::filter(n > 1)

  new_rows <- dplyr::bind_rows(new_by_outcome)
  combined_rows <- dplyr::bind_rows(combined_by_outcome)
  imputed_rows <- reviewed_data |>
    dplyr::filter(update_2026_r_source == "frozen canonical imputation")

  checks <- tibble::tibble(
    check = c(
      "prepared contrast keys are unique",
      "new arm identifiers do not overlap frozen identifiers",
      "new outcome identifiers do not overlap frozen identifiers",
      "frozen protein-intake centre is applied",
      "frozen pre-post correlation imputation is applied",
      "only eligible included rows reach the update interface",
      "outcome interfaces contain every ready row once",
      "combined interface preserves every frozen row",
      "combined interface contains every new row",
      "function outcome directions are complete"
    ),
    passed = c(
      nrow(duplicate_keys) == 0L,
      all(reviewed_data$arm_no > preprocessing_reference$max_arm_no),
      all(reviewed_data$outcome_no > preprocessing_reference$max_outcome_no),
      all(
        reviewed_data$protein_median ==
          preprocessing_reference$protein_intake_centre
      ),
      !nrow(imputed_rows) || all(
        imputed_rows$r == preprocessing_reference$pre_post_r_imputation
      ),
      all(new_rows$analysis_eligible_primary & new_rows$include_quantitative),
      nrow(new_rows) == sum(reviewed_data$update_2026_analysis_ready),
      sum(combined_rows$update_2026_data_source == "frozen_primary") ==
        nrow(frozen_analysis),
      sum(combined_rows$update_2026_data_source == "new_update") ==
        nrow(new_rows),
      all(!is.na(new_by_outcome$`function`$sign))
    )
  )

  if (any(!checks$passed)) {
    stop(
      "Prepared update validation failed: ",
      paste(checks$check[!checks$passed], collapse = "; "),
      call. = FALSE
    )
  }

  checks
}

summarise_update_pipeline_readiness <- function(
  update_data,
  reviewed_data,
  new_by_outcome,
  preprocessing_reference,
  decisions,
  workbook_checks,
  prepared_checks,
  fit_availability
) {
  study_status <- decisions |>
    dplyr::left_join(
      reviewed_data |>
        dplyr::group_by(study_no) |>
        dplyr::summarise(
          prepared_contrasts = dplyr::n(),
          eligible_before_review_decision = sum(analysis_eligible_primary),
          ready_for_update = sum(update_2026_analysis_ready),
          .groups = "drop"
        ),
      by = "study_no"
    ) |>
    dplyr::mutate(
      dplyr::across(
        c(
          prepared_contrasts,
          eligible_before_review_decision,
          ready_for_update
        ),
        ~ tidyr::replace_na(.x, 0L)
      )
    )

  contrast_flow <- tibble::tibble(
    stage = c(
      "prepared candidate contrasts",
      "eligible before review decisions",
      "excluded by review decision",
      "ready for incremental update"
    ),
    contrasts = c(
      nrow(reviewed_data),
      sum(reviewed_data$analysis_eligible_primary),
      sum(
        reviewed_data$analysis_eligible_primary &
          !reviewed_data$include_quantitative
      ),
      sum(reviewed_data$update_2026_analysis_ready)
    )
  )

  outcome_readiness <- purrr::imap_dfr(
    new_by_outcome,
    ~ tibble::tibble(
      outcome = .y,
      studies = dplyr::n_distinct(.x$study_no),
      arms = dplyr::n_distinct(.x$arm_no),
      contrasts = nrow(.x),
      protein_intake_raw_min = min(.x$protein_intake_raw, na.rm = TRUE),
      protein_intake_raw_max = max(.x$protein_intake_raw, na.rm = TRUE)
    )
  )

  reference <- tibble::tibble(
    item = names(preprocessing_reference),
    value = unlist(preprocessing_reference, use.names = FALSE)
  )

  raw_counts <- tibble::tibble(
    table = c("studies", "arms", "arm_timepoints", "reported_contrasts"),
    rows = c(
      nrow(update_data$study_characteristics),
      nrow(update_data$arm_characteristics),
      nrow(update_data$arm_based_data_by_timepoint),
      nrow(update_data$contrasts_between_timepoints)
    )
  )

  list(
    raw_counts = raw_counts,
    workbook_checks = dplyr::bind_rows(workbook_checks, prepared_checks),
    preprocessing_reference = reference,
    study_status = study_status,
    contrast_flow = contrast_flow,
    outcome_readiness = outcome_readiness,
    fit_availability = fit_availability
  )
}

write_update_pipeline_readiness_files <- function(
  readiness,
  output_dir = here::here("results", "update_2026")
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  files <- purrr::imap_chr(
    readiness,
    function(data, name) {
      path <- file.path(output_dir, paste0(name, ".csv"))
      readr::write_csv(data, path)
      path
    }
  )

  unname(files)
}
