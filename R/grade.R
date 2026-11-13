# GRADE evidence-profile and inconsistency-analysis helpers.
#
# The population effects remain those from the complete primary 26-model
# stacking ensemble. Inconsistency is estimated separately from the existing
# study-varying dose-effect candidates (models 9:14 and 21:26), re-stacked
# using their final corrected LOO objects. No function in this file fits a
# model.

# Outcome and contrast definitions ----

# Centralize the dose pairs, exercise groups and decision thresholds so the
# evidence domains are assessed for the same contrasts.
make_grade_primary_config <- function() {
  list(
    fixed_08_proposed_intakes = c(1.2, 1.6, 2.0, 2.4),
    fixed_16_proposed_by_outcome = list(
      whole_body = seq(2.0, 4.4, by = 0.4),
      hypertrophy = c(2.0, 2.4),
      strength = seq(2.0, 3.2, by = 0.4),
      `function` = seq(2.0, 3.2, by = 0.4)
    ),
    no_rt_observed_upper_by_outcome = c(
      whole_body = 1.8,
      hypertrophy = 1.8183864915572232,
      strength = 1.8183864915572232,
      `function` = 1.8183864915572232
    ),
    sesoi = 0.10,
    posterior_draws = 4000L,
    random_dose_model_indices = c(9:14, 21:26),
    di_flag = 50,
    asi_flag = 25,
    endpoint_neighbourhood = 0.10,
    outcomes = c(
      whole_body = "Whole-body lean/fat-free mass",
      hypertrophy = "Measures of muscle size or mass",
      strength = "Strength",
      `function` = "Physical function"
    ),
    rt_labels = c(`0` = "RT", `1` = "No RT"),
    seeds = c(
      whole_body = 78101L,
      hypertrophy = 78102L,
      strength = 78103L,
      `function` = 78104L
    )
  )
}

validate_grade_primary_config <- function(config) {
  required <- c(
    "fixed_08_proposed_intakes", "fixed_16_proposed_by_outcome",
    "no_rt_observed_upper_by_outcome",
    "sesoi", "posterior_draws",
    "random_dose_model_indices", "di_flag", "asi_flag",
    "endpoint_neighbourhood", "outcomes", "rt_labels", "seeds"
  )

  if (!is.list(config) || !all(required %in% names(config))) {
    stop("The GRADE primary configuration is incomplete.", call. = FALSE)
  }
  if (
    !identical(
      as.numeric(config$fixed_08_proposed_intakes),
      c(1.2, 1.6, 2.0, 2.4)
    ) ||
      !identical(
        lapply(config$fixed_16_proposed_by_outcome, as.numeric),
        list(
          whole_body = seq(2.0, 4.4, by = 0.4),
          hypertrophy = c(2.0, 2.4),
          strength = seq(2.0, 3.2, by = 0.4),
          `function` = seq(2.0, 3.2, by = 0.4)
        )
      ) ||
      !identical(
        names(config$no_rt_observed_upper_by_outcome),
        names(config$outcomes)
      ) ||
      any(!is.finite(config$no_rt_observed_upper_by_outcome)) ||
      !isTRUE(all.equal(config$sesoi, 0.10)) ||
      !identical(
        as.integer(config$random_dose_model_indices),
        as.integer(c(9:14, 21:26))
      ) ||
      length(config$outcomes) != 4L ||
      length(config$rt_labels) != 2L ||
      config$posterior_draws < 1000L
  ) {
    stop("The prespecified GRADE primary configuration is invalid.", call. = FALSE)
  }

  invisible(config)
}

grade_primary_contrasts <- function(
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  purrr::imap_dfr(
    config$fixed_16_proposed_by_outcome,
    function(fixed_16_proposed, outcome_name) {
      fixed_08 <- tibble::tibble(
        outcome = outcome_name,
        contrast_family_order = 1L,
        contrast_family = "Fixed 0.8 reference",
        protein_intake_reference = 0.8,
        protein_intake_proposed = config$fixed_08_proposed_intakes
      )
      fixed_16 <- tibble::tibble(
        outcome = outcome_name,
        contrast_family_order = 2L,
        contrast_family = "Fixed 1.6 reference",
        protein_intake_reference = 1.6,
        protein_intake_proposed = fixed_16_proposed
      )
      outcome_contrasts <- dplyr::bind_rows(fixed_08, fixed_16)
      rt_contrasts <- outcome_contrasts |>
        dplyr::mutate(resistance_exercise_code = 0L, .before = 2L)
      no_rt_contrasts <- outcome_contrasts |>
        dplyr::filter(
          .data$protein_intake_proposed <=
            config$no_rt_observed_upper_by_outcome[[outcome_name]] + 1e-10
        ) |>
        dplyr::mutate(resistance_exercise_code = 1L, .before = 2L)
      dplyr::bind_rows(rt_contrasts, no_rt_contrasts) |>
        dplyr::group_by(.data$resistance_exercise_code) |>
        dplyr::mutate(contrast_order = dplyr::row_number()) |>
        dplyr::ungroup()
    }
  ) |>
    dplyr::mutate(
      contrast_label = sprintf(
        "%.1f versus %.1f g/kg/day",
        .data$protein_intake_proposed,
        .data$protein_intake_reference
      )
    )
}

grade_primary_keys <- function(config = make_grade_primary_config()) {
  validate_grade_primary_config(config)
  grade_primary_contrasts(config) |>
    dplyr::mutate(
      outcome_label = unname(config$outcomes[.data$outcome]),
      rt_condition = unname(
        config$rt_labels[as.character(.data$resistance_exercise_code)]
      ),
      contrast_id = sprintf(
        "%s__%s__%.1f_vs_%.1f",
        .data$outcome,
        gsub(" ", "_", tolower(.data$rt_condition)),
        .data$protein_intake_proposed,
        .data$protein_intake_reference
      ),
      contrast_label = .data$contrast_label
    )
}

# Evidence sources and consistency checks ----

# Confirm that accepted fits, holdout scores, reported predictions and
# judgment inputs cover the required outcome and contrast keys.
make_grade_primary_preflight <- function(
    accepted_fit_manifest,
    final_loo_audit,
    reporting_audit,
    data_by_outcome,
    rob2_path,
    publication_bias_estimates,
    publication_bias_high_low,
    judgment_path,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  expected_outcomes <- names(config$outcomes)
  fit_files <- accepted_fit_manifest$fit_path
  fit_subset <- accepted_fit_manifest |>
    dplyr::filter(
      .data$model_index %in% config$random_dose_model_indices
    )
  reporting_checks <- reporting_audit$summary
  judgment_exists <- file.exists(judgment_path)
  checks <- tibble::tribble(
    ~check, ~passed, ~detail,
    "accepted_primary_fits",
    nrow(accepted_fit_manifest) == 104L && all(file.exists(fit_files)),
    sprintf("%d accepted fit records; %d files present", nrow(accepted_fit_manifest), sum(file.exists(fit_files))),
    "study_varying_candidate_fits",
    nrow(fit_subset) == 48L && all(file.exists(fit_subset$fit_path)),
    sprintf("%d of 48 required fit files present", sum(file.exists(fit_subset$fit_path))),
    "final_corrected_loo",
    isTRUE(final_loo_audit$stacking_ready) &&
      nrow(final_loo_audit$final_loo_manifest) == 104L,
    "Final corrected LOO must be complete and stacking-ready",
    "primary_reporting",
    nrow(reporting_checks) == 4L &&
      all(reporting_checks$reporting_checks_pass),
    "Full primary model-average reporting must pass all checks",
    "primary_analysis_data",
    identical(sort(names(data_by_outcome)), sort(expected_outcomes)) &&
      all(vapply(data_by_outcome, nrow, integer(1L)) > 0L),
    "All four primary outcome datasets must be available",
    "rob2_results",
    file.exists(rob2_path),
    normalizePath(rob2_path, winslash = "/", mustWork = FALSE),
    "publication_bias_results",
    is.data.frame(publication_bias_estimates) &&
      nrow(publication_bias_estimates) == 8L &&
      is.list(publication_bias_high_low) &&
      nrow(publication_bias_high_low$outcome_summary) == 4L,
    "Four outcome families with paired BRMA/RoBMA estimates are required",
    "grade_decision_sheet",
    judgment_exists,
    normalizePath(judgment_path, winslash = "/", mustWork = FALSE)
  )

  if (!all(checks$passed)) {
    failed <- checks$check[!checks$passed]
    stop(
      "GRADE primary preflight failed: ",
      paste(failed, collapse = ", "),
      call. = FALSE
    )
  }

  checks
}

write_grade_primary_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "grade", "preflight.csv"
    )
) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

# Study-specific contrasts for inconsistency ----

# Estimate weights for the subset of models that support the study-level
# predictions needed by the inconsistency calculation.
compute_grade_restricted_stacking_weights <- function(
    final_loo_audit,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  if (!isTRUE(final_loo_audit$stacking_ready)) {
    stop("Final primary LOO is not stacking-ready.", call. = FALSE)
  }
  manifest <- final_loo_audit$final_loo_manifest |>
    dplyr::filter(
      .data$model_index %in% config$random_dose_model_indices
    ) |>
    dplyr::arrange(
      factor(.data$outcome, levels = names(config$outcomes)),
      .data$model_index
    )

  if (
    nrow(manifest) != 48L ||
      any(!file.exists(manifest$final_loo_path)) ||
      anyDuplicated(manifest[c("outcome", "model_index")])
  ) {
    stop("The restricted GRADE LOO manifest is incomplete.", call. = FALSE)
  }

  weights <- purrr::map_dfr(names(config$outcomes), function(outcome_name) {
    records <- manifest |>
      dplyr::filter(.data$outcome == .env$outcome_name) |>
      dplyr::arrange(.data$model_index)
    results <- purrr::map(records$final_loo_path, qs::qread)
    loo_objects <- purrr::map(results, "loo")
    model_names <- purrr::map_chr(results, "model_name")
    names(loo_objects) <- model_names
    observed <- loo::loo_model_weights(loo_objects, method = "stacking")
    observed <- pmax(as.numeric(observed), 0)
    observed <- observed / sum(observed)
    tibble::tibble(
      outcome = outcome_name,
      model_index = records$model_index,
      model_name = model_names,
      restricted_stacking_weight = observed,
      loo_source = records$loo_source,
      final_loo_path = records$final_loo_path
    )
  })

  audit <- weights |>
    dplyr::group_by(.data$outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      weight_sum = sum(.data$restricted_stacking_weight),
      .groups = "drop"
    )

  if (
    nrow(weights) != 48L ||
      any(audit$models != 12L) ||
      any(abs(audit$weight_sum - 1) > 1e-8)
  ) {
    stop("Restricted GRADE stacking weights are invalid.", call. = FALSE)
  }

  weights
}

make_grade_inconsistency_records <- function(
    accepted_fit_manifest,
    restricted_weights,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  table <- restricted_weights |>
    dplyr::select(
      outcome, model_index, model_name,
      restricted_stacking_weight
    ) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          fit_stage = stage, fit_path
        ),
      by = c("outcome", "model_index")
    ) |>
    dplyr::arrange(
      factor(.data$outcome, levels = names(config$outcomes)),
      .data$model_index
    )

  if (
    nrow(table) != 48L || anyNA(table$fit_path) ||
      any(!file.exists(table$fit_path))
  ) {
    stop("The GRADE inconsistency fit records are incomplete.", call. = FALSE)
  }

  records <- purrr::map(names(config$outcomes), function(outcome_name) {
    models <- table |>
      dplyr::filter(.data$outcome == .env$outcome_name)
    list(
      outcome = outcome_name,
      seed = unname(config$seeds[[outcome_name]]),
      models = models
    )
  })

  names(records) <- names(config$outcomes)
  records
}

# Build study-specific low- and high-dose rows for the configured contrasts,
# and identify which study-dose combinations have empirical support.
make_grade_study_prediction_grid <- function(
    data,
    outcome,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  if (
    length(outcome) != 1L ||
      !outcome %in% names(config$outcomes)
  ) {
    stop("The GRADE prediction outcome is invalid.", call. = FALSE)
  }

  required <- c(
    "study_no", "arm_no", "outcome_no", "study_factor",
    "resistance_exercise_code", "protein_intake_raw", "protein_intake",
    "protein_median", "vi"
  )
  missing <- setdiff(required, names(data))

  if (length(missing)) {
    stop(
      "GRADE study prediction data are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  protein_median <- get_primary_protein_median(data)
  reference_rows <- data |>
    dplyr::filter(
      !is.na(.data$study_no),
      !is.na(.data$arm_no),
      !is.na(.data$outcome_no),
      !is.na(.data$study_factor),
      .data$resistance_exercise_code %in% c(0, 1)
    ) |>
    dplyr::arrange(
      .data$study_no, .data$resistance_exercise_code,
      .data$arm_no, .data$outcome_no
    ) |>
    dplyr::group_by(
      .data$study_no, .data$resistance_exercise_code
    ) |>
    dplyr::slice(1L) |>
    dplyr::ungroup()

  if (!nrow(reference_rows)) {
    stop("No study rows are available for GRADE predictions.", call. = FALSE)
  }

  contrasts <- grade_primary_contrasts(config) |>
    dplyr::filter(.data$outcome == .env$outcome)
  layout <- tibble::tibble(
    reference_row = seq_len(nrow(reference_rows)),
    resistance_exercise_code =
      reference_rows$resistance_exercise_code
  ) |>
    dplyr::left_join(
      contrasts,
      by = "resistance_exercise_code",
      relationship = "many-to-many"
    ) |>
    tidyr::crossing(endpoint = c("reference", "proposed")) |>
    dplyr::arrange(
      .data$reference_row,
      .data$contrast_order,
      factor(.data$endpoint, levels = c("reference", "proposed"))
    )

  newdata <- reference_rows[layout$reference_row, , drop = FALSE]
  doses <- ifelse(
    layout$endpoint == "reference",
    layout$protein_intake_reference,
    layout$protein_intake_proposed
  )
  newdata$protein_intake_raw <- doses
  newdata$protein_intake <- doses - protein_median
  newdata$protein_intake_05 <- newdata$protein_intake / 0.5
  newdata$vi <- 0

  if ("age_centre" %in% names(newdata)) newdata$age_centre <- 0
  if ("duration_centre" %in% names(newdata)) newdata$duration_centre <- 0

  moderator_data <- add_primary_standardised_moderators(data)
  moderator_centres <- attr(moderator_data, "moderator_centres")
  moderator_scales <- attr(moderator_data, "moderator_scales")
  newdata$age_outcome_z <-
    (newdata$age_centre - moderator_centres[["age_centre"]]) /
      moderator_scales[["age_centre"]]
  newdata$duration_outcome_z <-
    (newdata$duration_centre - moderator_centres[["duration_centre"]]) /
      moderator_scales[["duration_centre"]]

  prediction_keys <- tibble::tibble(
    prediction_column = seq_len(nrow(newdata)),
    study_no = as.integer(newdata$study_no),
    resistance_exercise_code = as.integer(
      newdata$resistance_exercise_code
    ),
    contrast_order = layout$contrast_order,
    contrast_family = layout$contrast_family,
    protein_intake_reference = layout$protein_intake_reference,
    protein_intake_proposed = layout$protein_intake_proposed,
    endpoint = layout$endpoint,
    dose = doses
  ) |>
    dplyr::mutate(
      cell_id = paste(
        .data$study_no,
        .data$resistance_exercise_code,
        sprintf("%.1f", .data$protein_intake_reference),
        sprintf("%.1f", .data$protein_intake_proposed),
        sep = "|"
      )
    )

  cell_keys <- prediction_keys |>
    dplyr::group_by(
      .data$cell_id, .data$study_no, .data$resistance_exercise_code,
      .data$contrast_order, .data$contrast_family,
      .data$protein_intake_reference, .data$protein_intake_proposed
    ) |>
    dplyr::summarise(
      reference_column = .data$prediction_column[
        .data$endpoint == "reference"
      ],
      proposed_column = .data$prediction_column[
        .data$endpoint == "proposed"
      ],
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$resistance_exercise_code,
      .data$contrast_order,
      .data$study_no
    ) |>
    dplyr::mutate(contrast_column = dplyr::row_number())

  if (
    nrow(prediction_keys) != 2L * nrow(cell_keys) ||
      anyDuplicated(cell_keys$cell_id) ||
      any(cell_keys$reference_column == cell_keys$proposed_column)
  ) {
    stop("The GRADE study prediction grid is invalid.", call. = FALSE)
  }

  list(
    newdata = newdata,
    prediction_keys = prediction_keys,
    cell_keys = cell_keys,
    protein_median = protein_median
  )
}

grade_inconsistency_draw_path <- function(outcome) {
  directory <- here::here(
    "model_archive", "update_2026", "grade", "inconsistency"
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, paste0(outcome, ".qs"))
}

# Pair low- and high-dose predictions within each posterior draw, then persist
# the study-level differences and the metadata needed for the inconsistency
# audit.
write_grade_study_contrast_draw_file <- function(
    record,
    data_by_outcome,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  data <- data_by_outcome[[record$outcome]]
  table <- record$models
  if (is.null(data) || !nrow(data) || nrow(table) != 12L) {
    stop("The GRADE inconsistency record is invalid.", call. = FALSE)
  }

  grid <- make_grade_study_prediction_grid(
    data = data,
    outcome = record$outcome,
    config = config
  )
  weights <- table$restricted_stacking_weight
  counts <- allocate_primary_mixture_models(
    weights,
    config$posterior_draws
  )
  set.seed(record$seed)
  allocation <- sample(rep(seq_along(weights), times = counts))
  contrast_draws <- matrix(
    NA_real_,
    nrow = config$posterior_draws,
    ncol = nrow(grid$cell_keys)
  )

  for (model_position in seq_len(nrow(table))) {
    output_rows <- which(allocation == model_position)
    if (!length(output_rows)) next
    candidate <- qs::qread(table$fit_path[[model_position]])
    if (
      !identical(candidate$outcome, record$outcome) ||
        as.integer(candidate$model_index) !=
          as.integer(table$model_index[[model_position]])
    ) {
      stop("An archived GRADE fit does not match its record.", call. = FALSE)
    }
    predictions <- brms::posterior_epred(
      candidate$model,
      newdata = grid$newdata,
      re_formula = NULL,
      ndraws = length(output_rows),
      allow_new_levels = FALSE
    ) |>
      as.matrix()
    if (
      ncol(predictions) != nrow(grid$newdata) &&
        nrow(predictions) == nrow(grid$newdata)
    ) {
      predictions <- t(predictions)
    }
    if (
      !identical(
        dim(predictions),
        c(length(output_rows), nrow(grid$newdata))
      ) ||
        anyNA(predictions) || any(!is.finite(predictions))
    ) {
      stop(
        "Study-specific prediction extraction failed for ",
        record$outcome, " model ", table$model_index[[model_position]],
        ".",
        call. = FALSE
      )
    }
    contrast_draws[output_rows, ] <-
      predictions[, grid$cell_keys$proposed_column, drop = FALSE] -
      predictions[, grid$cell_keys$reference_column, drop = FALSE]
    rm(candidate, predictions)
    gc(verbose = FALSE)
  }

  if (anyNA(contrast_draws) || any(!is.finite(contrast_draws))) {
    stop("The GRADE study-contrast mixture is incomplete.", call. = FALSE)
  }
  colnames(contrast_draws) <- grid$cell_keys$cell_id
  result <- list(
    outcome = record$outcome,
    draws = contrast_draws,
    cell_keys = grid$cell_keys,
    draw_manifest = tibble::tibble(
      drawid = seq_len(config$posterior_draws),
      model_index = table$model_index[allocation],
      model_name = table$model_name[allocation]
    ),
    restricted_weights = table,
    settings = list(
      contrasts = grade_primary_contrasts(config),
      sesoi = config$sesoi,
      posterior_draws = config$posterior_draws,
      protein_median = grid$protein_median,
      seed = record$seed,
      estimand = paste(
        "latent study-specific expected-outcome contrast;",
        "same study/arm/outcome identifiers at both dose endpoints"
      )
    )
  )

  class(result) <- c("grade_study_contrast_draws", class(result))
  output_file <- grade_inconsistency_draw_path(record$outcome)
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(result, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise GRADE study-contrast draws.", call. = FALSE)
  }

  unlink(temporary_file, force = TRUE)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

audit_grade_study_contrast_draw_files <- function(
    paths,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != 4L || any(!file.exists(paths))) {
    stop("The GRADE study-contrast draw files are incomplete.", call. = FALSE)
  }
  objects <- purrr::map(paths, qs::qread)
  outcomes <- purrr::map_chr(objects, "outcome")
  valid <- purrr::map_lgl(objects, function(x) {
    observed_contrasts <- x$cell_keys |>
      dplyr::distinct(
        .data$resistance_exercise_code,
        .data$contrast_family,
        .data$protein_intake_reference,
        .data$protein_intake_proposed
      )
    expected_contrasts <- grade_primary_contrasts(config) |>
      dplyr::filter(.data$outcome == x$outcome) |>
      dplyr::select(
        "resistance_exercise_code",
        "contrast_family",
        "protein_intake_reference",
        "protein_intake_proposed"
      )
    inherits(x, "grade_study_contrast_draws") &&
      nrow(x$draws) == config$posterior_draws &&
      ncol(x$draws) == nrow(x$cell_keys) &&
      all(is.finite(x$draws)) &&
      nrow(observed_contrasts) == nrow(expected_contrasts) &&
      !nrow(dplyr::anti_join(
        observed_contrasts,
        expected_contrasts,
        by = names(expected_contrasts)
      )) &&
      !nrow(dplyr::anti_join(
        expected_contrasts,
        observed_contrasts,
        by = names(expected_contrasts)
      ))
  })

  if (
    !setequal(outcomes, names(config$outcomes)) ||
      anyDuplicated(outcomes) || !all(valid)
  ) {
    stop("A GRADE study-contrast draw file failed audit.", call. = FALSE)
  }

  tibble::tibble(
    outcome = outcomes,
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(paths)$size),
    draws = purrr::map_int(objects, ~ nrow(.x$draws)),
    contrasts = purrr::map_int(objects, ~ ncol(.x$draws)),
    studies = purrr::map_int(
      objects,
      ~ dplyr::n_distinct(.x$cell_keys$study_no)
    )
  ) |>
    dplyr::arrange(factor(.data$outcome, levels = names(config$outcomes)))
}

# Summarize inconsistency and empirical support ----

# Use the study-by-draw matrix and the smallest effect of interest to compute
# the two inconsistency measures used in the evidence profile.
calculate_grade_di_asi <- function(draw_matrix, sesoi = 0.10) {
  draw_matrix <- as.matrix(draw_matrix)
  if (
    !is.numeric(draw_matrix) || !nrow(draw_matrix) || !ncol(draw_matrix) ||
      anyNA(draw_matrix) || any(!is.finite(draw_matrix)) ||
      length(sesoi) != 1L || !is.finite(sesoi) || sesoi <= 0
  ) {
    stop("The DI/ASI draw matrix is invalid.", call. = FALSE)
  }

  categories <- matrix(
    ifelse(
      draw_matrix <= -sesoi,
      "meaningful detriment",
      ifelse(
        draw_matrix > sesoi,
        "meaningful benefit",
        "practical equivalence"
      )
    ),
    nrow = nrow(draw_matrix),
    ncol = ncol(draw_matrix),
    dimnames = dimnames(draw_matrix)
  )

  levels <- c(
    "meaningful detriment", "practical equivalence",
    "meaningful benefit"
  )
  pooled_counts <- table(factor(as.vector(categories), levels = levels))
  pooled_proportions <- as.numeric(pooled_counts) / sum(pooled_counts)
  di <- 100 * (
    1 - 0.5 * sum(abs(1 / 3 - pooled_proportions)) / (2 / 3)
  )

  study_table <- t(vapply(
    seq_len(ncol(categories)),
    function(column) {
      as.numeric(table(factor(categories[, column], levels = levels)))
    },
    numeric(3L)
  ))
  rownames(study_table) <- colnames(draw_matrix)
  colnames(study_table) <- levels
  active <- colSums(study_table) > 0L
  active_table <- study_table[, active, drop = FALSE]
  n_studies <- nrow(study_table)
  n_categories <- ncol(active_table)
  asi <- if (n_studies < 2L) {
    NA_real_
  } else if (n_categories == 1L) {
    0
  } else {
    test <- suppressWarnings(stats::chisq.test(
      active_table,
      correct = FALSE
    ))
    denominator <- sum(active_table) *
      min(nrow(active_table) - 1L, ncol(active_table) - 1L)
    cramers_v <- sqrt(as.numeric(test$statistic) / denominator)
    if (n_categories == 2L && n_studies > 2L) {
      100 * cramers_v / sqrt(2)
    } else {
      100 * cramers_v
    }
  }

  list(
    DI = as.numeric(di),
    ASI = as.numeric(asi),
    class_distribution = tibble::tibble(
      classification = levels,
      samples = as.integer(pooled_counts),
      proportion = pooled_proportions
    ),
    study_classification = study_table,
    categories = categories
  )
}

make_grade_empirical_support <- function(
    data_by_outcome,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  data <- purrr::imap_dfr(data_by_outcome, function(x, outcome_name) {
    dplyr::mutate(x, outcome = outcome_name, .before = 1L)
  })
  contrasts <- grade_primary_contrasts(config)
  study_support <- data |>
    dplyr::filter(
      is.finite(.data$protein_intake_raw),
      .data$resistance_exercise_code %in% c(0, 1)
    ) |>
    dplyr::group_by(
      .data$outcome, .data$resistance_exercise_code, .data$study_no
    ) |>
    dplyr::summarise(
      study_min = min(.data$protein_intake_raw),
      study_max = max(.data$protein_intake_raw),
      observed_arms = dplyr::n_distinct(.data$arm_no),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      contrasts,
      by = c("outcome", "resistance_exercise_code"),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      reference_near_observed =
        .data$study_min <= .data$protein_intake_reference +
          config$endpoint_neighbourhood &
        .data$study_max >= .data$protein_intake_reference -
          config$endpoint_neighbourhood,
      proposed_near_observed =
        .data$study_min <= .data$protein_intake_proposed +
          config$endpoint_neighbourhood &
        .data$study_max >= .data$protein_intake_proposed -
          config$endpoint_neighbourhood,
      spans_contrast =
        .data$study_min <= pmin(
          .data$protein_intake_reference,
          .data$protein_intake_proposed
        ) &
        .data$study_max >= pmax(
          .data$protein_intake_reference,
          .data$protein_intake_proposed
        )
    )

  arm_counts <- data |>
    dplyr::filter(
      .data$resistance_exercise_code %in% c(0, 1),
      is.finite(.data$n), .data$n > 0
    ) |>
    dplyr::group_by(
      .data$outcome, .data$resistance_exercise_code,
      .data$study_no, .data$arm_no
    ) |>
    dplyr::summarise(n = max(.data$n), .groups = "drop") |>
    dplyr::group_by(.data$outcome, .data$resistance_exercise_code) |>
    dplyr::summarise(
      participants = sum(.data$n),
      arms = dplyr::n(),
      .groups = "drop"
    )

  summary <- study_support |>
    dplyr::group_by(
      .data$outcome, .data$resistance_exercise_code,
      .data$contrast_order, .data$contrast_family,
      .data$protein_intake_reference, .data$protein_intake_proposed
    ) |>
    dplyr::summarise(
      studies = dplyr::n_distinct(.data$study_no),
      studies_spanning_contrast = sum(.data$spans_contrast),
      studies_near_reference = sum(.data$reference_near_observed),
      studies_near_proposed = sum(.data$proposed_near_observed),
      observed_min = min(.data$study_min),
      observed_max = max(.data$study_max),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      arm_counts,
      by = c("outcome", "resistance_exercise_code")
    ) |>
    dplyr::mutate(
      reference_in_observed_range =
        .data$protein_intake_reference >= .data$observed_min &
        .data$protein_intake_reference <= .data$observed_max,
      proposed_in_observed_range =
        .data$protein_intake_proposed >= .data$observed_min &
        .data$protein_intake_proposed <= .data$observed_max,
      direct_spanning_fraction =
        .data$studies_spanning_contrast / .data$studies
    )

  list(summary = summary, study_support = study_support)
}

# Population effects and other evidence domains ----

# Match published population contrasts to the GRADE definitions before joining
# risk-of-bias, publication-bias and inconsistency evidence.
make_grade_population_effects <- function(
    reporting_audit,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  manifest <- reporting_audit$manifest |>
    dplyr::arrange(factor(.data$outcome, levels = names(config$outcomes)))
  support <- reporting_audit$support
  if (
    nrow(manifest) != 4L ||
      any(!file.exists(manifest$model_average_path))
  ) {
    stop("Primary model-average files are unavailable for GRADE.", call. = FALSE)
  }

  contrasts <- grade_primary_contrasts(config)
  selected <- purrr::map_dfr(seq_len(nrow(manifest)), function(position) {
    model_average <- qs::qread(manifest$model_average_path[[position]])
    outcome_name <- manifest$outcome[[position]]
    purrr::map_dfr(c(0L, 1L), function(rt_code) {
      outcome_contrasts <- contrasts |>
        dplyr::filter(
          .data$outcome == .env$outcome_name,
          .data$resistance_exercise_code == .env$rt_code
        )
      observed <- support |>
        dplyr::filter(
          .data$outcome == .env$outcome_name,
          .data$resistance_exercise_code == .env$rt_code
        )
      if (nrow(observed) != 1L) {
        stop("Primary support did not align with GRADE.", call. = FALSE)
      }
      purrr::map_dfr(seq_len(nrow(outcome_contrasts)), function(contrast_position) {
        contrast <- outcome_contrasts[contrast_position, ]
        reference <- contrast$protein_intake_reference[[1L]]
        proposed <- contrast$protein_intake_proposed[[1L]]
        reference_column <- which(
          model_average$grid$resistance_exercise_code == rt_code &
            abs(model_average$grid$protein_intake - reference) < 1e-10
        )
        proposed_column <- which(
          model_average$grid$resistance_exercise_code == rt_code &
            abs(model_average$grid$protein_intake - proposed) < 1e-10
        )
        if (
          length(reference_column) != 1L ||
            length(proposed_column) != 1L
        ) {
          stop(
            "The primary mixture lacks a requested GRADE contrast endpoint.",
            call. = FALSE
          )
        }
        draws <-
          model_average$draws[, proposed_column] -
          model_average$draws[, reference_column]
        lower <- stats::quantile(draws, 0.025, names = FALSE)
        upper <- stats::quantile(draws, 0.975, names = FALSE)
        probability_benefit <- mean(draws > config$sesoi)
        probability_equivalence <- mean(
          draws >= -config$sesoi & draws <= config$sesoi
        )
        probability_detriment <- mean(draws < -config$sesoi)
        supported <-
          reference >= observed$observed_min &
          reference <= observed$observed_max &
          proposed >= observed$observed_min &
          proposed <= observed$observed_max
        tibble::tibble(
          outcome = outcome_name,
          resistance_exercise_code = rt_code,
          resistance_exercise = ifelse(rt_code == 0L, "YES", "NO"),
          rt_condition = ifelse(rt_code == 0L, "RT", "No RT"),
          contrast_order = contrast$contrast_order[[1L]],
          contrast_family = contrast$contrast_family[[1L]],
          protein_intake_reference = reference,
          protein_intake_proposed = proposed,
          estimate = mean(draws),
          lower_95 = lower,
          upper_95 = upper,
          probability_benefit = probability_benefit,
          probability_equivalence = probability_equivalence,
          probability_detriment = probability_detriment,
          posterior_conclusion = dplyr::case_when(
            probability_benefit >= 0.95 ~ "meaningful benefit",
            probability_equivalence >= 0.95 ~ "practical equivalence",
            probability_detriment >= 0.95 ~ "meaningful detriment",
            TRUE ~ "uncertain"
          ),
          imprecision_suggestion = dplyr::case_when(
            lower < -config$sesoi & upper > config$sesoi ~ "very serious",
            (lower < -config$sesoi & upper >= -config$sesoi) |
              (lower <= config$sesoi & upper > config$sesoi) ~ "serious",
            TRUE ~ "not serious"
          ),
          supported = supported,
          support_status = ifelse(
            supported,
            "Observed range",
            "Extrapolated beyond outcome/RT observed range"
          ),
          observed_min = observed$observed_min,
          observed_max = observed$observed_max
        )
      })
    })
  })

  expected_rows <- nrow(contrasts)
  if (nrow(selected) != expected_rows) {
    stop(
      "Expected ", expected_rows,
      " population contrasts for the GRADE profile; found ",
      nrow(selected), ".",
      call. = FALSE
    )
  }

  selected
}

make_grade_rob2_summary <- function(
    data_by_outcome,
    rob2_path,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  rob <- readr::read_csv(rob2_path, show_col_types = FALSE) |>
    dplyr::transmute(
      study_no = as.integer(.data$study_no),
      overall_rob2 = .data$overall_final
    )
  studies <- purrr::imap_dfr(data_by_outcome, function(data, outcome_name) {
    data |>
      dplyr::filter(.data$resistance_exercise_code %in% c(0, 1)) |>
      dplyr::distinct(
        outcome = outcome_name,
        .data$resistance_exercise_code,
        .data$study_no
      )
  }) |>
    dplyr::left_join(rob, by = "study_no")

  summary <- studies |>
    dplyr::group_by(.data$outcome, .data$resistance_exercise_code) |>
    dplyr::summarise(
      contributing_studies = dplyr::n_distinct(.data$study_no),
      assessed_studies = sum(!is.na(.data$overall_rob2)),
      low = sum(.data$overall_rob2 == "Low", na.rm = TRUE),
      some_concerns = sum(
        .data$overall_rob2 == "Some concerns", na.rm = TRUE
      ),
      high = sum(.data$overall_rob2 == "High", na.rm = TRUE),
      unassessed_study_numbers = paste(
        sort(unique(.data$study_no[is.na(.data$overall_rob2)])),
        collapse = ";"
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      low_proportion = .data$low / .data$assessed_studies,
      some_concerns_proportion =
        .data$some_concerns / .data$assessed_studies,
      high_proportion = .data$high / .data$assessed_studies
    )

  list(summary = summary, study_judgments = studies)
}

make_grade_publication_bias_summary <- function(
    model_estimates,
    high_low_data,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  outcome_map <- c(
    `whole body lean mass` = "whole_body",
    `direct muscle size/mass` = "hypertrophy",
    strength = "strength",
    `function` = "function"
  )
  estimates <- model_estimates |>
    dplyr::mutate(
      outcome = unname(outcome_map[.data$outcome_group])
    ) |>
    dplyr::filter(!is.na(.data$outcome)) |>
    dplyr::select(
      outcome, model_type, estimate, lower_95, upper_95
    ) |>
    tidyr::pivot_wider(
      names_from = model_type,
      values_from = c(estimate, lower_95, upper_95)
    ) |>
    dplyr::mutate(
      adjusted_minus_unadjusted =
        .data$estimate_robma - .data$estimate_brma,
      absolute_attenuation =
        abs(.data$estimate_brma) - abs(.data$estimate_robma)
    )

  counts <- high_low_data$outcome_summary |>
    dplyr::mutate(
      outcome = unname(outcome_map[.data$outcome_group])
    ) |>
    dplyr::select(
      outcome,
      publication_units = independent_studies,
      selected_contrasts,
      candidate_contrasts
    )

  summary <- estimates |>
    dplyr::left_join(counts, by = "outcome")
  if (nrow(summary) != 4L || anyNA(summary$outcome)) {
    stop("Publication-bias GRADE inputs are incomplete.", call. = FALSE)
  }
  summary
}

summarise_grade_inconsistency <- function(
    draw_manifest,
    empirical_support,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  draw_manifest <- draw_manifest |>
    dplyr::arrange(factor(.data$outcome, levels = names(config$outcomes)))
  objects <- purrr::map(draw_manifest$path, qs::qread)
  summary_rows <- list()
  study_rows <- list()
  loo_rows <- list()
  counter <- 0L
  contrasts <- grade_primary_contrasts(config)

  for (object in objects) {
    for (rt_code in c(0L, 1L)) {
      outcome_contrasts <- contrasts |>
        dplyr::filter(
          .data$outcome == object$outcome,
          .data$resistance_exercise_code == .env$rt_code
        )
      for (contrast_position in seq_len(nrow(outcome_contrasts))) {
        contrast <- outcome_contrasts[contrast_position, ]
        reference <- contrast$protein_intake_reference[[1L]]
        proposed <- contrast$protein_intake_proposed[[1L]]
        counter <- counter + 1L
        keys <- object$cell_keys |>
          dplyr::filter(
            .data$resistance_exercise_code == .env$rt_code,
            abs(.data$protein_intake_reference - .env$reference) < 1e-10,
            abs(.data$protein_intake_proposed - .env$proposed) < 1e-10
          ) |>
          dplyr::arrange(.data$study_no)
        matrix <- object$draws[, keys$contrast_column, drop = FALSE]
        colnames(matrix) <- as.character(keys$study_no)
        index <- calculate_grade_di_asi(matrix, config$sesoi)
        study_sd <- if (ncol(matrix) >= 2L) {
          apply(matrix, 1L, stats::sd)
        } else {
          rep(NA_real_, nrow(matrix))
        }
        support <- empirical_support$study_support |>
          dplyr::filter(
            .data$outcome == object$outcome,
            .data$resistance_exercise_code == .env$rt_code,
            abs(.data$protein_intake_reference - .env$reference) < 1e-10,
            abs(.data$protein_intake_proposed - .env$proposed) < 1e-10
          ) |>
          dplyr::arrange(.data$study_no)
        support <- support[match(keys$study_no, support$study_no), ]
        if (anyNA(support$study_no)) {
          stop("Study support did not align with inconsistency draws.", call. = FALSE)
        }
        direct_columns <- which(support$spans_contrast)
        direct_index <- if (length(direct_columns) >= 2L) {
          calculate_grade_di_asi(
            matrix[, direct_columns, drop = FALSE],
            config$sesoi
          )
        } else {
          list(DI = NA_real_, ASI = NA_real_)
        }
        classes <- index$class_distribution
        summary_rows[[counter]] <- tibble::tibble(
          outcome = object$outcome,
          resistance_exercise_code = rt_code,
          contrast_order = contrast$contrast_order[[1L]],
          contrast_family = contrast$contrast_family[[1L]],
          protein_intake_reference = reference,
          protein_intake_proposed = proposed,
          studies = ncol(matrix),
          directly_spanning_studies = length(direct_columns),
          DI = index$DI,
          ASI = index$ASI,
          direct_support_DI = direct_index$DI,
          direct_support_ASI = direct_index$ASI,
          study_sd_median = stats::median(study_sd, na.rm = TRUE),
          study_sd_lower_95 = stats::quantile(
            study_sd, 0.025, na.rm = TRUE, names = FALSE
          ),
          study_sd_upper_95 = stats::quantile(
            study_sd, 0.975, na.rm = TRUE, names = FALSE
          ),
          probability_detriment = classes$proportion[
            classes$classification == "meaningful detriment"
          ],
          probability_equivalence = classes$proportion[
            classes$classification == "practical equivalence"
          ],
          probability_benefit = classes$proportion[
            classes$classification == "meaningful benefit"
          ],
          DI_flag = index$DI >= config$di_flag,
          ASI_flag = index$ASI >= config$asi_flag,
          inconsistency_suggestion = dplyr::case_when(
            is.na(index$ASI) ~ "unable to assess",
            index$ASI >= config$asi_flag ~ "review for serious inconsistency",
            index$DI >= config$di_flag ~
              "decision dispersion present; primarily review imprecision",
            TRUE ~ "not serious suggested"
          )
        )
        category_matrix <- index$categories
        study_rows[[counter]] <- purrr::map_dfr(
          seq_len(ncol(matrix)),
          function(column) {
            values <- matrix[, column]
            categories <- category_matrix[, column]
            tibble::tibble(
              outcome = object$outcome,
              resistance_exercise_code = rt_code,
              contrast_order = contrast$contrast_order[[1L]],
              contrast_family = contrast$contrast_family[[1L]],
              protein_intake_reference = reference,
              protein_intake_proposed = proposed,
              study_no = as.integer(colnames(matrix)[[column]]),
              estimate = mean(values),
              lower_95 = stats::quantile(
                values, 0.025, names = FALSE
              ),
              upper_95 = stats::quantile(
                values, 0.975, names = FALSE
              ),
              probability_detriment = mean(
                categories == "meaningful detriment"
              ),
              probability_equivalence = mean(
                categories == "practical equivalence"
              ),
              probability_benefit = mean(
                categories == "meaningful benefit"
              ),
              spans_contrast = support$spans_contrast[[column]],
              reference_near_observed =
                support$reference_near_observed[[column]],
              proposed_near_observed =
                support$proposed_near_observed[[column]]
            )
          }
        )
        loo_rows[[counter]] <- if (ncol(matrix) >= 3L) {
          purrr::map_dfr(seq_len(ncol(matrix)), function(column) {
            omitted <- calculate_grade_di_asi(
              matrix[, -column, drop = FALSE],
              config$sesoi
            )
            tibble::tibble(
              outcome = object$outcome,
              resistance_exercise_code = rt_code,
              contrast_order = contrast$contrast_order[[1L]],
              contrast_family = contrast$contrast_family[[1L]],
              protein_intake_reference = reference,
              protein_intake_proposed = proposed,
              omitted_study_no = as.integer(colnames(matrix)[[column]]),
              DI = omitted$DI,
              ASI = omitted$ASI
            )
          })
        } else {
          tibble::tibble(
            outcome = character(),
            resistance_exercise_code = integer(),
            contrast_order = integer(),
            contrast_family = character(),
            protein_intake_reference = numeric(),
            protein_intake_proposed = numeric(),
            omitted_study_no = integer(),
            DI = numeric(),
            ASI = numeric()
          )
        }
      }
    }
  }

  list(
    summary = dplyr::bind_rows(summary_rows),
    study_summary = dplyr::bind_rows(study_rows),
    leave_one_out = dplyr::bind_rows(loo_rows)
  )
}

# Recorded judgments and evidence profiles ----

# Read the review team's domain judgments, check that every required key is
# covered, and retain the rationale in the published evidence profile.
read_grade_domain_judgments <- function(
    path,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  judgments <- readr::read_csv(
    path,
    show_col_types = FALSE,
    na = c("", "NA")
  )
  required <- c(
    "outcome", "rt_condition", "protein_intake_reference",
    "protein_intake_proposed", "starting_certainty",
    "risk_of_bias", "risk_of_bias_downgrade", "risk_of_bias_rationale",
    "inconsistency", "inconsistency_downgrade", "inconsistency_rationale",
    "indirectness", "indirectness_downgrade", "indirectness_rationale",
    "imprecision", "imprecision_downgrade", "imprecision_rationale",
    "publication_bias", "publication_bias_downgrade",
    "publication_bias_rationale", "other_considerations",
    "other_adjustment", "importance"
  )

  missing <- setdiff(required, names(judgments))
  if (length(missing)) {
    stop(
      "The GRADE decision sheet is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  expected <- grade_primary_keys(config) |>
    dplyr::select(
      outcome, rt_condition,
      protein_intake_reference,
      protein_intake_proposed
    )
  observed <- judgments |>
    dplyr::select(dplyr::all_of(names(expected)))

  if (
    nrow(judgments) != nrow(expected) || anyDuplicated(observed) ||
      nrow(dplyr::anti_join(expected, observed, by = names(expected))) ||
      nrow(dplyr::anti_join(observed, expected, by = names(expected)))
  ) {
    stop(
      "The GRADE decision sheet must contain exactly the ", nrow(expected),
      " prespecified ",
      "outcome/RT/contrast rows.",
      call. = FALSE
    )
  }

  judgments
}

grade_certainty_label <- function(score) {
  dplyr::case_when(
    is.na(score) ~ "Pending",
    score >= 4 ~ "High",
    score == 3 ~ "Moderate",
    score == 2 ~ "Low",
    TRUE ~ "Very low"
  )
}

grade_certainty_symbols <- function(score) {
  dplyr::case_when(
    is.na(score) ~ "Pending",
    score >= 4 ~ "⊕⊕⊕⊕",
    score == 3 ~ "⊕⊕⊕◯",
    score == 2 ~ "⊕⊕◯◯",
    TRUE ~ "⊕◯◯◯"
  )
}

assemble_grade_evidence_profile <- function(
    population_effects,
    empirical_support,
    inconsistency,
    rob2,
    publication_bias,
    judgments,
    config = make_grade_primary_config()
) {
  validate_grade_primary_config(config)
  keys <- grade_primary_keys(config)
  profile <- keys |>
    dplyr::left_join(
      population_effects,
      by = c(
        "outcome", "resistance_exercise_code",
        "rt_condition", "contrast_order", "contrast_family",
        "protein_intake_reference",
        "protein_intake_proposed"
      )
    ) |>
    dplyr::left_join(
      empirical_support$summary,
      by = c(
        "outcome", "resistance_exercise_code",
        "contrast_order", "contrast_family",
        "protein_intake_reference", "protein_intake_proposed"
      ),
      suffix = c("", "_empirical")
    ) |>
    dplyr::left_join(
      inconsistency$summary,
      by = c(
        "outcome", "resistance_exercise_code",
        "contrast_order", "contrast_family",
        "protein_intake_reference", "protein_intake_proposed"
      ),
      suffix = c("", "_inconsistency")
    ) |>
    dplyr::left_join(
      rob2$summary,
      by = c("outcome", "resistance_exercise_code")
    ) |>
    dplyr::left_join(publication_bias, by = "outcome") |>
    dplyr::left_join(
      judgments,
      by = c(
        "outcome", "rt_condition", "protein_intake_reference",
        "protein_intake_proposed"
      )
    )

  downgrade_columns <- c(
    "risk_of_bias_downgrade", "inconsistency_downgrade",
    "indirectness_downgrade", "imprecision_downgrade",
    "publication_bias_downgrade", "other_adjustment"
  )
  downgrade_matrix <- as.matrix(profile[, downgrade_columns])
  complete_judgment <-
    is.finite(profile$starting_certainty) &
    apply(downgrade_matrix, 1L, function(x) all(is.finite(x)))
  profile$total_adjustment <- NA_real_
  profile$total_adjustment[complete_judgment] <- rowSums(
    downgrade_matrix[complete_judgment, , drop = FALSE]
  )
  profile$certainty_score <- NA_real_
  profile$certainty_score[complete_judgment] <- pmax(
    1,
    pmin(
      4,
      profile$starting_certainty[complete_judgment] -
        profile$total_adjustment[complete_judgment]
    )
  )

  profile |>
    dplyr::mutate(
      certainty = grade_certainty_label(.data$certainty_score),
      certainty_symbols = grade_certainty_symbols(.data$certainty_score),
      study_design = "Randomised trials; achieved-intake dose-response analysis"
    ) |>
    dplyr::arrange(
      factor(.data$outcome, levels = names(config$outcomes)),
      .data$resistance_exercise_code,
      .data$contrast_order
    )
}

write_grade_analysis_bundle <- function(
    preflight,
    restricted_weights,
    draw_manifest,
    population_effects,
    empirical_support,
    inconsistency,
    rob2,
    publication_bias,
    judgments,
    evidence_profile,
    output_dir = here::here("results", "update_2026", "grade")
) {
  expected_rows <- nrow(grade_primary_keys(make_grade_primary_config()))
  if (!all(preflight$passed) || nrow(evidence_profile) != expected_rows) {
    stop("Cannot write an incomplete GRADE analysis bundle.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    restricted_weights = file.path(
      output_dir, "inconsistency_restricted_stacking_weights.csv"
    ),
    draw_manifest = file.path(
      output_dir, "inconsistency_draw_manifest.csv"
    ),
    population_effects = file.path(output_dir, "population_effects.csv"),
    empirical_support = file.path(output_dir, "empirical_support.csv"),
    study_support = file.path(output_dir, "study_support.csv.gz"),
    inconsistency = file.path(output_dir, "inconsistency_summary.csv"),
    inconsistency_studies = file.path(
      output_dir, "inconsistency_study_summary.csv.gz"
    ),
    inconsistency_loo = file.path(
      output_dir, "inconsistency_leave_one_out.csv.gz"
    ),
    rob2 = file.path(output_dir, "rob2_summary.csv"),
    rob2_studies = file.path(output_dir, "rob2_study_judgments.csv.gz"),
    publication_bias = file.path(output_dir, "publication_bias_summary.csv"),
    judgments = file.path(output_dir, "domain_judgments_snapshot.csv"),
    evidence_profile = file.path(output_dir, "evidence_profile.csv")
  )

  readr::write_csv(restricted_weights, paths[["restricted_weights"]])
  readr::write_csv(draw_manifest, paths[["draw_manifest"]])
  readr::write_csv(population_effects, paths[["population_effects"]])
  readr::write_csv(empirical_support$summary, paths[["empirical_support"]])
  readr::write_csv(empirical_support$study_support, paths[["study_support"]])
  readr::write_csv(inconsistency$summary, paths[["inconsistency"]])
  readr::write_csv(
    inconsistency$study_summary,
    paths[["inconsistency_studies"]]
  )
  readr::write_csv(
    inconsistency$leave_one_out,
    paths[["inconsistency_loo"]]
  )
  readr::write_csv(rob2$summary, paths[["rob2"]])
  readr::write_csv(rob2$study_judgments, paths[["rob2_studies"]])
  readr::write_csv(publication_bias, paths[["publication_bias"]])
  readr::write_csv(judgments, paths[["judgments"]])
  readr::write_csv(evidence_profile, paths[["evidence_profile"]])

  manifest <- tibble::tibble(
    artifact = names(paths),
    path = normalizePath(
      unname(paths), winslash = "/", mustWork = TRUE
    ),
    bytes = as.numeric(file.info(unname(paths))$size),
    sha256 = purrr::map_chr(
      unname(paths),
      ~ digest::digest(file = .x, algo = "sha256")
    )
  )

  manifest_path <- file.path(output_dir, "manifest.csv")
  readr::write_csv(manifest, manifest_path)
  normalizePath(
    c(unname(paths), manifest_path),
    winslash = "/",
    mustWork = TRUE
  )
}
