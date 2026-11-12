# Study-balanced publication-bias models and supplementary reporting ----
#
# The primary inferential data retain one highest-versus-lowest achieved-
# intake contrast per study and outcome. All eligible pairwise contrasts are
# retained for descriptive figures and auditing, but are not treated as
# independent publication units in the RoBMA fits.

make_publication_bias_all_pair_contrasts <- function(data_by_outcome) {
  if (!is.list(data_by_outcome) || !length(data_by_outcome)) {
    stop("Publication-bias input must be a non-empty outcome list.",
         call. = FALSE)
  }

  data <- dplyr::bind_rows(data_by_outcome)
  required <- c(
    "study_no", "arm_no", "arm_no_extracted", "condition",
    "outcome_group", "measurement_method", "unit_of_measurement",
    "time1", "time2", "resistance_exercise_code", "yi", "vi", "n",
    "protein_intake", "protein_intake_raw"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "All-pair contrast input is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"authors" %in% names(data)) {
    data$authors <- paste("Study", data$study_no)
  }
  if (!"update_2026_data_source" %in% names(data)) {
    data$update_2026_data_source <- NA_character_
  }

  keys <- c(
    "study_no", "outcome_group", "measurement_method",
    "unit_of_measurement", "time1", "time2",
    "resistance_exercise_code"
  )
  arms <- data |>
    dplyr::filter(
      .data$condition %in% c("protein", "control", "placebo"),
      is.finite(.data$yi),
      is.finite(.data$vi),
      .data$vi > 0,
      is.finite(.data$n),
      .data$n > 0,
      is.finite(.data$protein_intake_raw)
    ) |>
    dplyr::mutate(publication_bias_row_id = dplyr::row_number())

  select_arm <- function(data, suffix) {
    data |>
      dplyr::transmute(
        dplyr::across(dplyr::all_of(keys)),
        row_id = .data$publication_bias_row_id,
        authors = as.character(.data$authors),
        arm_no = .data$arm_no,
        arm_no_extracted = .data$arm_no_extracted,
        condition = as.character(.data$condition),
        yi = .data$yi,
        vi = .data$vi,
        n = .data$n,
        protein_intake = .data$protein_intake,
        protein_intake_raw = .data$protein_intake_raw,
        update_2026_data_source =
          as.character(.data$update_2026_data_source)
      ) |>
      dplyr::rename_with(
        ~ paste0(.x, "_", suffix),
        -dplyr::all_of(keys)
      )
  }

  arm_a <- select_arm(arms, "a")
  arm_b <- select_arm(arms, "b")
  paired <- dplyr::inner_join(
    arm_a,
    arm_b,
    by = keys,
    relationship = "many-to-many"
  ) |>
    dplyr::filter(
      .data$row_id_a < .data$row_id_b,
      .data$condition_a == "protein" | .data$condition_b == "protein"
    ) |>
    dplyr::mutate(
      # Positive effects always correspond to the arm with the higher
      # achieved intake. Exact intake ties are oriented protein minus
      # comparator, followed by a stable arm-identifier tie-break.
      a_is_high = dplyr::case_when(
        .data$protein_intake_raw_a > .data$protein_intake_raw_b ~ TRUE,
        .data$protein_intake_raw_a < .data$protein_intake_raw_b ~ FALSE,
        .data$condition_a == "protein" &
          .data$condition_b != "protein" ~ TRUE,
        .data$condition_b == "protein" &
          .data$condition_a != "protein" ~ FALSE,
        TRUE ~ as.character(.data$arm_no_extracted_a) >
          as.character(.data$arm_no_extracted_b)
      ),
      high_arm_no = dplyr::if_else(
        .data$a_is_high, .data$arm_no_a, .data$arm_no_b
      ),
      low_arm_no = dplyr::if_else(
        .data$a_is_high, .data$arm_no_b, .data$arm_no_a
      ),
      high_arm_no_extracted = dplyr::if_else(
        .data$a_is_high,
        .data$arm_no_extracted_a,
        .data$arm_no_extracted_b
      ),
      low_arm_no_extracted = dplyr::if_else(
        .data$a_is_high,
        .data$arm_no_extracted_b,
        .data$arm_no_extracted_a
      ),
      high_condition = dplyr::if_else(
        .data$a_is_high, .data$condition_a, .data$condition_b
      ),
      low_condition = dplyr::if_else(
        .data$a_is_high, .data$condition_b, .data$condition_a
      ),
      yi_high = dplyr::if_else(
        .data$a_is_high, .data$yi_a, .data$yi_b
      ),
      yi_low = dplyr::if_else(
        .data$a_is_high, .data$yi_b, .data$yi_a
      ),
      vi_high = dplyr::if_else(
        .data$a_is_high, .data$vi_a, .data$vi_b
      ),
      vi_low = dplyr::if_else(
        .data$a_is_high, .data$vi_b, .data$vi_a
      ),
      n_high = dplyr::if_else(
        .data$a_is_high, .data$n_a, .data$n_b
      ),
      n_low = dplyr::if_else(
        .data$a_is_high, .data$n_b, .data$n_a
      ),
      high_protein_intake = dplyr::if_else(
        .data$a_is_high,
        .data$protein_intake_raw_a,
        .data$protein_intake_raw_b
      ),
      low_protein_intake = dplyr::if_else(
        .data$a_is_high,
        .data$protein_intake_raw_b,
        .data$protein_intake_raw_a
      ),
      yi = .data$yi_high - .data$yi_low,
      vi = .data$vi_high + .data$vi_low,
      sei = sqrt(.data$vi),
      z_statistic = .data$yi / .data$sei,
      p_value_two_sided = 2 * stats::pnorm(
        abs(.data$z_statistic), lower.tail = FALSE
      ),
      total_n = .data$n_high + .data$n_low,
      dose_difference =
        .data$high_protein_intake - .data$low_protein_intake,
      dose_midpoint =
        (.data$high_protein_intake + .data$low_protein_intake) / 2,
      comparison_type = dplyr::if_else(
        .data$condition_a == "protein" &
          .data$condition_b == "protein",
        "Active dose vs active dose",
        "Protein vs comparator"
      ),
      resistance_exercise_label = dplyr::case_when(
        .data$resistance_exercise_code == 0 ~ "With resistance training",
        .data$resistance_exercise_code == 1 ~ "Without resistance training",
        TRUE ~ "Resistance-training status unavailable"
      ),
      authors = dplyr::coalesce(.data$authors_a, .data$authors_b),
      update_2026_data_source = dplyr::coalesce(
        .data$update_2026_data_source_a,
        .data$update_2026_data_source_b
      ),
      effect_id = paste(
        .data$study_no,
        .data$outcome_group,
        .data$measurement_method,
        .data$unit_of_measurement,
        .data$time1,
        .data$time2,
        .data$high_arm_no_extracted,
        .data$low_arm_no_extracted,
        sep = "::"
      )
    ) |>
    dplyr::select(
      dplyr::all_of(c(
        keys,
        "authors",
        "update_2026_data_source",
        "effect_id",
        "comparison_type",
        "high_arm_no",
        "low_arm_no",
        "high_arm_no_extracted",
        "low_arm_no_extracted",
        "high_condition",
        "low_condition",
        "high_protein_intake",
        "low_protein_intake",
        "dose_difference",
        "dose_midpoint",
        "n_high",
        "n_low",
        "total_n",
        "yi_high",
        "yi_low",
        "vi_high",
        "vi_low",
        "yi",
        "vi",
        "sei",
        "z_statistic",
        "p_value_two_sided",
        "resistance_exercise_label"
      ))
    ) |>
    dplyr::arrange(
      factor(
        .data$outcome_group,
        levels = names(publication_bias_outcome_labels())
      ),
      .data$study_no,
      .data$effect_id
    )

  if (
    !nrow(paired) || anyDuplicated(paired$effect_id) ||
      any(!is.finite(paired$yi)) || any(paired$vi <= 0) ||
      any(paired$dose_difference < 0)
  ) {
    stop("All-pair publication-bias contrasts failed validation.",
         call. = FALSE)
  }
  paired
}

select_publication_bias_high_low_contrasts <- function(contrasts) {
  required <- c(
    "study_no", "outcome_group", "dose_difference", "time2", "total_n",
    "measurement_method", "unit_of_measurement", "high_arm_no_extracted",
    "low_arm_no_extracted", "effect_id", "yi", "vi", "sei",
    "p_value_two_sided", "comparison_type", "resistance_exercise_label"
  )
  missing <- setdiff(required, names(contrasts))
  if (length(missing)) {
    stop(
      "High-low selection is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  audit <- contrasts |>
    dplyr::mutate(
      measurement_order = dplyr::coalesce(
        as.character(.data$measurement_method), ""
      ),
      unit_order = dplyr::coalesce(
        as.character(.data$unit_of_measurement), ""
      )
    ) |>
    dplyr::group_by(.data$study_no, .data$outcome_group) |>
    dplyr::arrange(
      dplyr::desc(.data$dose_difference),
      dplyr::desc(.data$time2),
      dplyr::desc(.data$total_n),
      .data$measurement_order,
      .data$unit_order,
      .data$high_arm_no_extracted,
      .data$low_arm_no_extracted,
      .data$effect_id,
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      selection_rank = dplyr::row_number(),
      candidate_count = dplyr::n(),
      selected_index_contrast = .data$selection_rank == 1L,
      selection_hierarchy = paste(
        "largest achieved-intake separation; latest follow-up;",
        "largest combined sample; measurement, unit, arm and effect",
        "identifier tie-breakers"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      -dplyr::all_of(c("measurement_order", "unit_order"))
    )

  selected <- audit |>
    dplyr::filter(.data$selected_index_contrast) |>
    dplyr::mutate(
      outcome_group = factor(
        .data$outcome_group,
        levels = names(publication_bias_outcome_labels())
      )
    ) |>
    dplyr::arrange(.data$outcome_group, .data$study_no) |>
    dplyr::mutate(
      outcome_group = as.character(.data$outcome_group),
      comparison_type = factor(
        .data$comparison_type,
        levels = c(
          "Protein vs comparator",
          "Active dose vs active dose"
        )
      ),
      resistance_exercise_label = factor(
        .data$resistance_exercise_label,
        levels = c(
          "With resistance training",
          "Without resistance training",
          "Resistance-training status unavailable"
        )
      )
    )

  expected <- contrasts |>
    dplyr::distinct(.data$study_no, .data$outcome_group)
  observed <- selected |>
    dplyr::distinct(.data$study_no, .data$outcome_group)
  if (
    nrow(selected) != nrow(expected) ||
      anyDuplicated(paste(observed$study_no, observed$outcome_group)) ||
      !setequal(
        paste(expected$study_no, expected$outcome_group),
        paste(observed$study_no, observed$outcome_group)
      )
  ) {
    stop("High-low selection did not retain one contrast per study-outcome.",
         call. = FALSE)
  }

  summary <- selected |>
    dplyr::group_by(.data$outcome_group) |>
    dplyr::summarise(
      independent_studies = dplyr::n_distinct(.data$study_no),
      selected_contrasts = dplyr::n(),
      candidate_contrasts = sum(.data$candidate_count),
      protein_comparator = sum(
        .data$comparison_type == "Protein vs comparator"
      ),
      active_dose = sum(
        .data$comparison_type == "Active dose vs active dose"
      ),
      median_dose_difference = stats::median(.data$dose_difference),
      median_total_n = stats::median(.data$total_n),
      statistically_significant = sum(.data$p_value_two_sided < 0.05),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      outcome_label = unname(
        publication_bias_outcome_labels()[.data$outcome_group]
      ),
      .before = 1
    )

  list(
    candidate_audit = audit,
    index_contrasts = selected,
    outcome_summary = summary
  )
}

write_publication_bias_extended_data <- function(
    all_pair_contrasts,
    high_low_data,
    nominal_data,
    output_dir = here::here(
      "results", "update_2026", "publication_bias"
    ),
    hosted_dir = here::here(
      "supplement", "data", "publication-bias"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(hosted_dir, recursive = TRUE, showWarnings = FALSE)
  objects <- list(
    all_pair_contrasts = all_pair_contrasts,
    high_low_candidate_audit = high_low_data$candidate_audit,
    high_low_index_contrasts = high_low_data$index_contrasts,
    high_low_outcome_summary = high_low_data$outcome_summary,
    nominal_index_contrasts = nominal_data$index_contrasts,
    nominal_outcome_summary = nominal_data$outcome_summary
  )
  paths <- file.path(output_dir, paste0(names(objects), ".csv"))
  hosted_paths <- file.path(
    hosted_dir, paste0(names(objects), ".csv")
  )
  purrr::walk2(objects, paths, readr::write_csv)
  purrr::walk2(objects, hosted_paths, readr::write_csv)
  all_paths <- c(paths, hosted_paths)
  if (
    any(!file.exists(all_paths)) ||
      any(file.info(all_paths)$size <= 0)
  ) {
    stop("Extended publication-bias data files were not written.",
         call. = FALSE)
  }
  normalizePath(all_paths, winslash = "/", mustWork = TRUE)
}

make_publication_bias_all_contrast_panel <- function(
    all_data,
    selected_data,
    outcome_group,
    model_estimates = NULL
) {
  panel_all <- all_data |>
    dplyr::filter(.data$outcome_group == .env$outcome_group)
  panel_selected <- selected_data |>
    dplyr::filter(.data$outcome_group == .env$outcome_group)
  if (!nrow(panel_all) || !nrow(panel_selected)) {
    stop("Missing all-contrast funnel data for ", outcome_group,
         call. = FALSE)
  }

  se_limit <- max(panel_all$sei) * 1.08
  x_limit <- max(
    abs(panel_all$yi),
    stats::qnorm(0.995) * se_limit,
    0.25
  ) * 1.05
  contour_grid <- make_publication_bias_contour_grid(x_limit, se_limit)
  contour_lines <- tidyr::expand_grid(
    standard_error = seq(0, se_limit, length.out = 200L),
    alpha = c(0.10, 0.05, 0.01),
    side = c(-1, 1)
  ) |>
    dplyr::mutate(
      effect = .data$side *
        stats::qnorm(1 - .data$alpha / 2) *
        .data$standard_error,
      contour = factor(
        .data$alpha,
        levels = c(0.10, 0.05, 0.01),
        labels = c("p = 0.10", "p = 0.05", "p = 0.01")
      ),
      line_group = interaction(.data$alpha, .data$side, drop = TRUE)
    )

  plot <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = contour_grid,
      ggplot2::aes(
        x = .data$effect,
        y = .data$standard_error,
        fill = .data$significance_region
      )
    ) +
    ggplot2::geom_line(
      data = contour_lines,
      ggplot2::aes(
        x = .data$effect,
        y = .data$standard_error,
        group = .data$line_group,
        linetype = .data$contour
      ),
      colour = "grey58",
      linewidth = 0.32
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.45) +
    ggplot2::geom_point(
      data = panel_all,
      ggplot2::aes(
        x = .data$yi,
        y = .data$sei,
        colour = .data$comparison_type
      ),
      size = 1.25,
      alpha = 0.18
    ) +
    ggplot2::geom_point(
      data = panel_selected,
      ggplot2::aes(
        x = .data$yi,
        y = .data$sei,
        fill = .data$comparison_type
      ),
      shape = 21,
      colour = "black",
      stroke = 0.35,
      size = 2.5,
      alpha = 0.95
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "p < 0.01" = "#D5E3F0",
        "0.01 <= p < 0.05" = "#E5ECF3",
        "0.05 <= p < 0.10" = "#F0F3F6",
        "p >= 0.10" = "#FFFFFF",
        "Protein vs comparator" = "#0072B2",
        "Active dose vs active dose" = "#D55E00"
      ),
      name = "Selected contrast / null contour"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Protein vs comparator" = "#0072B2",
        "Active dose vs active dose" = "#D55E00"
      ),
      name = "All eligible contrasts"
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "p = 0.10" = "dotted",
        "p = 0.05" = "dashed",
        "p = 0.01" = "dotdash",
        "Unadjusted brma" = "longdash",
        "RoBMA adjusted" = "solid"
      ),
      name = "Contour / model estimate"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-x_limit, x_limit),
      oob = scales::squish,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_reverse(
      limits = c(se_limit, 0),
      oob = scales::squish,
      expand = ggplot2::expansion(mult = c(0, 0.015))
    ) +
    ggplot2::labs(
      title = unname(publication_bias_outcome_labels()[outcome_group]),
      subtitle = paste0(
        nrow(panel_all), " eligible contrasts from ",
        dplyr::n_distinct(panel_all$study_no), " studies; ",
        nrow(panel_selected), " selected study-level contrasts"
      ),
      x = "Higher-minus-lower achieved-intake effect (SMD)",
      y = "Standard error"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(model_estimates) && nrow(model_estimates)) {
    lines <- model_estimates |>
      dplyr::filter(.data$outcome_group == .env$outcome_group)
    if (nrow(lines)) {
      plot <- plot + ggplot2::geom_vline(
        data = lines,
        ggplot2::aes(
          xintercept = .data$estimate,
          linetype = .data$model_label
        ),
        colour = "black",
        linewidth = 0.65
      )
    }
  }
  plot
}

make_publication_bias_all_contrast_plot <- function(
    all_pair_contrasts,
    high_low_data,
    model_estimates = NULL
) {
  outcomes <- names(publication_bias_outcome_labels())
  plots <- purrr::map(
    outcomes,
    ~ make_publication_bias_all_contrast_panel(
      all_data = all_pair_contrasts,
      selected_data = high_low_data$index_contrasts,
      outcome_group = .x,
      model_estimates = model_estimates
    )
  )
  patchwork::wrap_plots(plots, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = paste(
        "Contour-enhanced funnels of all eligible within-study",
        "dose contrasts"
      ),
      subtitle = paste(
        "Faint points show all derived contrasts; outlined points show the",
        "single highest-versus-lowest achieved-intake contrast selected",
        "per study and outcome."
      )
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left"
    )
}

write_publication_bias_all_contrast_plot <- function(
    plot,
    output_dir = here::here("plots", "update_2026", "publication_bias"),
    filename_stem = "all-contrast-funnel-plots"
) {
  if (
    length(filename_stem) != 1L || is.na(filename_stem) ||
      !grepl("^[A-Za-z0-9_-]+$", filename_stem)
  ) {
    stop("filename_stem must be a safe single file stem.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(
    output_dir,
    paste0(filename_stem, c(".png", ".tiff", ".pdf"))
  )
  ggplot2::ggsave(
    paths[[1L]], plot, device = "png", dpi = 180,
    width = 11, height = 9
  )
  ggplot2::ggsave(
    paths[[2L]], plot, device = "tiff", dpi = 300,
    width = 11, height = 9, compression = "lzw"
  )
  ggplot2::ggsave(
    paths[[3L]], plot, device = grDevices::cairo_pdf,
    width = 11, height = 9
  )
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    stop("All-contrast funnel plots were not written.", call. = FALSE)
  }
  normalizePath(paths, winslash = "/", mustWork = TRUE)
}

check_publication_bias_robma_runtime <- function(
    required_version = "4.0.0"
) {
  installed <- requireNamespace("RoBMA", quietly = TRUE)
  version <- if (installed) {
    as.character(utils::packageVersion("RoBMA"))
  } else {
    NA_character_
  }
  version_ok <- installed &&
    utils::compareVersion(version, required_version) >= 0
  exports_ok <- FALSE
  if (version_ok) {
    required_exports <- c(
      "brma", "RoBMA", "as_zplot", "zplot", "summary_models",
      "pooled_effect"
    )
    exports_ok <- all(
      required_exports %in% getNamespaceExports("RoBMA")
    )
  }
  jags_version <- NA_character_
  if (installed && requireNamespace("rjags", quietly = TRUE)) {
    jags_version <- tryCatch(
      as.character(rjags::jags.version()),
      error = function(error) NA_character_
    )
  }
  jags_ok <- !is.na(jags_version) &&
    utils::compareVersion(jags_version, "4.3.1") >= 0

  tibble::tibble(
    component = c("RoBMA package", "RoBMA 4.0 API", "JAGS"),
    detected = c(
      ifelse(installed, version, "not installed"),
      ifelse(exports_ok, "available", "unavailable"),
      ifelse(is.na(jags_version), "not detected", jags_version)
    ),
    requirement = c(
      paste0(">= ", required_version),
      "brma/RoBMA/as_zplot/zplot/summary_models/pooled_effect",
      ">= 4.3.1"
    ),
    ready = c(version_ok, exports_ok, jags_ok)
  )
}

assert_publication_bias_robma_runtime <- function(preflight) {
  if (
    !is.data.frame(preflight) || !nrow(preflight) ||
      any(!preflight$ready)
  ) {
    stop(
      paste(
        "Publication-bias fitting requires RoBMA >= 4.0.0 and JAGS",
        ">= 4.3.1. Install JAGS, install and snapshot RoBMA in renv,",
        "then rerun the preflight before enabling the fit branch."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

write_publication_bias_robma_preflight <- function(
    preflight,
    path = here::here(
      "results", "update_2026", "publication_bias", "robma_preflight.csv"
    )
) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(preflight, path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

make_publication_bias_robma_config <- function(
    chains = 4L,
    adapt = 5000L,
    burnin = 5000L,
    sample = 10000L,
    seed = 62026L
) {
  values <- c(chains, adapt, burnin, sample, seed)
  if (any(!is.finite(values)) || any(values < 1)) {
    stop("Invalid publication-bias RoBMA configuration.", call. = FALSE)
  }
  list(
    chains = as.integer(chains),
    adapt = as.integer(adapt),
    burnin = as.integer(burnin),
    sample = as.integer(sample),
    seed = as.integer(seed),
    parallel = TRUE,
    required_robma_version = "4.0.0"
  )
}

publication_bias_outcome_slugs <- function() {
  c(
    "whole body lean mass" = "whole_body",
    "direct muscle size/mass" = "direct_muscle",
    "strength" = "strength",
    "function" = "function"
  )
}

make_publication_bias_model_specifications <- function(high_low_data) {
  data <- high_low_data$index_contrasts
  outcomes <- names(publication_bias_outcome_labels())
  missing <- setdiff(outcomes, unique(data$outcome_group))
  if (length(missing)) {
    stop(
      "Publication-bias model data are missing outcomes: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  outcome_specs <- purrr::map_dfr(outcomes, function(outcome) {
    tibble::tibble(
      analysis_id = paste0(
        "primary_", publication_bias_outcome_slugs()[[outcome]]
      ),
      analysis_role = "Primary outcome-specific",
      outcome_group = outcome,
      cluster_by_study = FALSE,
      moderator_formula = NA_character_,
      expected_rows = sum(data$outcome_group == outcome),
      expected_studies = dplyr::n_distinct(
        data$study_no[data$outcome_group == outcome]
      )
    )
  })
  omnibus <- tibble::tibble(
    analysis_id = "secondary_omnibus",
    analysis_role = "Secondary clustered omnibus meta-regression",
    outcome_group = NA_character_,
    cluster_by_study = TRUE,
    moderator_formula = paste(
      "~ outcome_group + dose_difference_z + comparison_type +",
      "resistance_exercise_label"
    ),
    expected_rows = nrow(data),
    expected_studies = dplyr::n_distinct(data$study_no)
  )
  dplyr::bind_rows(outcome_specs, omnibus)
}

make_publication_bias_fit_records <- function(specifications) {
  records <- tidyr::crossing(
    specifications,
    model_type = c("brma", "robma")
  ) |>
    dplyr::mutate(
      fit_id = paste(.data$analysis_id, .data$model_type, sep = "__"),
      seed_offset = dplyr::row_number()
    )
  if (nrow(records) != 10L || anyDuplicated(records$fit_id)) {
    stop("Expected ten unique publication-bias fit records.",
         call. = FALSE)
  }
  split(records, records$fit_id)
}

prepare_publication_bias_fit_data <- function(record, high_low_data) {
  if (is.list(record) && !is.data.frame(record)) {
    record <- tibble::as_tibble(record)
  }
  if (nrow(record) != 1L) {
    stop("A publication-bias fit record must contain one row.",
         call. = FALSE)
  }
  data <- high_low_data$index_contrasts
  if (!is.na(record$outcome_group[[1L]])) {
    data <- data |>
      dplyr::filter(
        .data$outcome_group == record$outcome_group[[1L]]
      )
  }
  dose_sd <- stats::sd(data$dose_difference)
  if (!is.finite(dose_sd) || dose_sd <= 0) {
    dose_sd <- 1
  }
  data <- data |>
    dplyr::mutate(
      outcome_group = factor(
        .data$outcome_group,
        levels = names(publication_bias_outcome_labels())
      ),
      comparison_type = factor(
        .data$comparison_type,
        levels = c(
          "Protein vs comparator",
          "Active dose vs active dose"
        )
      ),
      resistance_exercise_label = factor(
        .data$resistance_exercise_label,
        levels = c(
          "With resistance training",
          "Without resistance training",
          "Resistance-training status unavailable"
        )
      ),
      dose_difference_z =
        (.data$dose_difference - mean(.data$dose_difference)) / dose_sd
    ) |>
    dplyr::mutate(
      dplyr::across(where(is.factor), droplevels)
    )
  if (
    nrow(data) != record$expected_rows[[1L]] ||
      dplyr::n_distinct(data$study_no) !=
        record$expected_studies[[1L]] ||
      any(!is.finite(data$yi)) || any(!is.finite(data$vi)) ||
      any(data$vi <= 0)
  ) {
    stop("Publication-bias fit data failed structural validation.",
         call. = FALSE)
  }
  if (!is.na(record$moderator_formula[[1L]])) {
    moderator_formula <- stats::as.formula(
      record$moderator_formula[[1L]]
    )
    model_frame <- stats::model.frame(
      moderator_formula,
      data = data,
      na.action = stats::na.pass
    )
    model_matrix <- stats::model.matrix(moderator_formula, data = data)
    if (anyNA(model_frame) || qr(model_matrix)$rank < ncol(model_matrix)) {
      stop(
        "Publication-bias moderator matrix contains missing or aliased terms.",
        call. = FALSE
      )
    }
  }
  data
}

fit_publication_bias_model_to_file <- function(
    record,
    high_low_data,
    preflight,
    config,
    enabled = FALSE,
    output_dir = here::here(
      "results", "update_2026", "publication_bias", "robma", "fits"
    )
) {
  if (!isTRUE(enabled)) {
    stop(
      paste(
        "Publication-bias fitting is disabled. Set",
        "PROTEIN_PUBLICATION_BIAS_ROBMA=true only after the active model",
        "family refits have finished."
      ),
      call. = FALSE
    )
  }
  assert_publication_bias_robma_runtime(preflight)
  record <- tibble::as_tibble(record)
  data <- prepare_publication_bias_fit_data(record, high_low_data)
  fit_id <- record$fit_id[[1L]]
  model_type <- record$model_type[[1L]]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(fit_id, ".rds"))
  if (file.exists(path)) {
    existing <- tryCatch(readRDS(path), error = function(error) NULL)
    reusable <- is.list(existing) &&
      all(c("fit", "record", "data", "config", "runtime") %in%
            names(existing)) &&
      identical(existing$record, record) &&
      identical(existing$data, data) &&
      identical(existing$config, config) &&
      identical(existing$runtime, preflight)
    if (isTRUE(reusable)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
  }

  # RoBMA 4.0 resolves its default mean-difference factor contrast by name.
  # Attaching the package within the isolated worker makes contr.meandif
  # visible to model.matrix without attaching RoBMA to the full pipeline.
  robma_was_attached <- "package:RoBMA" %in% search()
  if (!robma_was_attached) {
    suppressPackageStartupMessages(
      library("RoBMA", character.only = TRUE)
    )
    on.exit(
      detach(
        "package:RoBMA",
        unload = FALSE,
        character.only = TRUE
      ),
      add = TRUE
    )
  }

  arguments <- list(
    yi = data$yi,
    vi = data$vi,
    measure = "SMD",
    data = data,
    slab = data$effect_id,
    chains = config$chains,
    parallel = config$parallel,
    adapt = config$adapt,
    burnin = config$burnin,
    sample = config$sample,
    seed = config$seed + record$seed_offset[[1L]],
    silent = TRUE
  )
  if (isTRUE(record$cluster_by_study[[1L]])) {
    arguments$cluster <- data$study_no
  }
  if (!is.na(record$moderator_formula[[1L]])) {
    arguments$mods <- stats::as.formula(
      record$moderator_formula[[1L]]
    )
  }
  fit_function <- if (identical(model_type, "brma")) {
    RoBMA::brma
  } else if (identical(model_type, "robma")) {
    RoBMA::RoBMA
  } else {
    stop("Unknown publication-bias model type: ", model_type,
         call. = FALSE)
  }
  fit <- do.call(fit_function, arguments)

  saveRDS(
    list(
      fit = fit,
      record = record,
      data = data,
      config = config,
      runtime = preflight
    ),
    path,
    compress = "xz"
  )
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

make_publication_bias_fit_manifest_record <- function(record, fit_file) {
  record <- tibble::as_tibble(record)
  valid_fit_file <- length(fit_file) == 1L &&
    !is.na(fit_file) && nzchar(fit_file) && file.exists(fit_file)
  if (nrow(record) != 1L || !isTRUE(valid_fit_file)) {
    stop("Publication-bias fit manifest record is invalid.",
         call. = FALSE)
  }
  record |>
    dplyr::mutate(
      fit_path = normalizePath(
        fit_file, winslash = "/", mustWork = TRUE
      ),
      fit_size_bytes = file.info(fit_file)$size
    )
}

compile_publication_bias_fit_manifest <- function(records) {
  manifest <- dplyr::bind_rows(records) |>
    dplyr::arrange(.data$analysis_id, .data$model_type)
  if (
    nrow(manifest) != 10L || anyDuplicated(manifest$fit_id) ||
      any(!file.exists(manifest$fit_path))
  ) {
    stop("Publication-bias fit manifest is incomplete.", call. = FALSE)
  }
  manifest
}

write_publication_bias_fit_manifest <- function(
    manifest,
    path = here::here(
      "results", "update_2026", "publication_bias", "robma",
      "fit_manifest.csv"
    )
) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

extract_publication_bias_outcome_model_estimates <- function(manifest) {
  required <- c(
    "analysis_id", "outcome_group", "model_type", "fit_path"
  )
  if (!all(required %in% names(manifest))) {
    stop("Publication-bias fit manifest is missing required fields.",
         call. = FALSE)
  }
  rows <- manifest |>
    dplyr::filter(grepl("^primary_", .data$analysis_id))
  if (
    nrow(rows) != 8L ||
      !all(rows$model_type %in% c("brma", "robma")) ||
      any(!file.exists(rows$fit_path))
  ) {
    stop("Outcome-specific publication-bias fits are incomplete.",
         call. = FALSE)
  }

  estimates <- purrr::map_dfr(seq_len(nrow(rows)), function(index) {
    row <- rows[index, , drop = FALSE]
    object <- readRDS(row$fit_path[[1L]])
    if (!is.list(object) || is.null(object$fit)) {
      stop("Publication-bias fit file has an unexpected structure.",
           call. = FALSE)
    }
    samples <- RoBMA::pooled_effect(
      object$fit,
      bias_adjusted = identical(row$model_type[[1L]], "robma")
    )
    draws <- as.numeric(as.matrix(samples))
    if (!length(draws) || any(!is.finite(draws))) {
      stop("Could not extract finite pooled-effect draws.", call. = FALSE)
    }
    tibble::tibble(
      outcome_group = row$outcome_group[[1L]],
      model_type = row$model_type[[1L]],
      model_label = ifelse(
        row$model_type[[1L]] == "robma",
        "RoBMA adjusted",
        "Unadjusted brma"
      ),
      estimate = mean(draws),
      lower_95 = stats::quantile(draws, 0.025, names = FALSE),
      upper_95 = stats::quantile(draws, 0.975, names = FALSE)
    )
  })
  if (
    anyDuplicated(estimates[c("outcome_group", "model_type")]) ||
      any(!is.finite(estimates$estimate))
  ) {
    stop("Extracted publication-bias model estimates are invalid.",
         call. = FALSE)
  }
  estimates |>
    dplyr::arrange(.data$outcome_group, .data$model_type)
}

write_publication_bias_outcome_model_estimates <- function(
    estimates,
    result_path = here::here(
      "results", "update_2026", "publication_bias", "robma",
      "outcome_model_estimates.csv"
    ),
    hosted_path = here::here(
      "supplement", "data", "publication-bias",
      "outcome_model_estimates.csv"
    )
) {
  required <- c(
    "outcome_group", "model_type", "model_label", "estimate",
    "lower_95", "upper_95"
  )
  if (
    !all(required %in% names(estimates)) || nrow(estimates) != 8L ||
      any(!is.finite(estimates$estimate))
  ) {
    stop("Outcome model estimates are incomplete.", call. = FALSE)
  }
  paths <- c(result_path, hosted_path)
  purrr::walk(unique(dirname(paths)), dir.create,
              recursive = TRUE, showWarnings = FALSE)
  purrr::walk(paths, ~ readr::write_csv(estimates, .x))
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    stop("Outcome model estimate files were not written.", call. = FALSE)
  }
  normalizePath(paths, winslash = "/", mustWork = TRUE)
}

write_publication_bias_model_reporting_bundle <- function(
    manifest,
    high_low_data,
    output_dir = here::here(
      "results", "update_2026", "publication_bias", "robma", "reporting"
    ),
    figure_dir = here::here(
      "plots", "update_2026", "publication_bias", "robma"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # Summary and zplot methods for meta-regressions reconstruct the moderator
  # model matrix and therefore need RoBMA's named factor contrasts visible.
  robma_was_attached <- "package:RoBMA" %in% search()
  if (!robma_was_attached) {
    suppressPackageStartupMessages(
      library("RoBMA", character.only = TRUE)
    )
    on.exit(
      detach(
        "package:RoBMA",
        unload = FALSE,
        character.only = TRUE
      ),
      add = TRUE
    )
  }

  analysis_ids <- unique(manifest$analysis_id)
  written <- character()

  for (analysis_id in analysis_ids) {
    rows <- manifest |>
      dplyr::filter(.data$analysis_id == .env$analysis_id)
    if (!setequal(rows$model_type, c("brma", "robma"))) {
      stop("Each publication-bias analysis requires brma and RoBMA fits.",
           call. = FALSE)
    }
    fitted <- stats::setNames(
      lapply(rows$fit_path, readRDS),
      rows$model_type
    )
    fit_brma <- fitted[["brma"]]$fit
    fit_robma <- fitted[["robma"]]$fit

    text_objects <- list(
      brma_summary = capture.output(summary(fit_brma)),
      robma_summary = capture.output(summary(fit_robma)),
      robma_models = capture.output(RoBMA::summary_models(fit_robma)),
      robma_zplot = capture.output(
        summary(RoBMA::as_zplot(fit_robma))
      )
    )
    text_paths <- file.path(
      output_dir,
      paste0(analysis_id, "-", names(text_objects), ".txt")
    )
    purrr::walk2(text_objects, text_paths, writeLines)
    written <- c(written, text_paths)

    observed <- fitted[["robma"]]$data$z_statistic
    lower <- max(-8, floor(min(observed, na.rm = TRUE) * 4) / 4 - 0.25)
    upper <- min(8, ceiling(max(observed, na.rm = TRUE) * 4) / 4 + 0.25)
    if (lower >= upper) {
      lower <- -4
      upper <- 4
    }
    z_brma <- RoBMA::as_zplot(fit_brma)
    z_robma <- RoBMA::as_zplot(fit_robma)
    outcome <- rows$outcome_group[[1L]]
    title <- if (is.na(outcome)) {
      "Secondary Omnibus Analysis"
    } else {
      publication_bias_outcome_labels()[[outcome]]
    }
    png_path <- file.path(figure_dir, paste0(analysis_id, "-zplot.png"))
    pdf_path <- file.path(figure_dir, paste0(analysis_id, "-zplot.pdf"))

    draw_zplot <- function() {
      graphics::hist(
        z_robma,
        from = lower,
        to = upper,
        by = 0.25,
        main = title
      )
      graphics::lines(
        z_brma,
        from = lower,
        to = upper,
        col = "black",
        lty = 2,
        lwd = 2
      )
      graphics::lines(
        z_robma,
        from = lower,
        to = upper,
        col = "#0072B2",
        lty = 1,
        lwd = 2
      )
      graphics::lines(
        z_robma,
        from = lower,
        to = upper,
        extrapolate = TRUE,
        col = "#0072B2",
        lty = 3,
        lwd = 2
      )
      graphics::legend(
        "topright",
        legend = c(
          "Unadjusted brma",
          "RoBMA fitted",
          "RoBMA without selection"
        ),
        col = c("black", "#0072B2", "#0072B2"),
        lty = c(2, 1, 3),
        lwd = 2,
        bty = "n"
      )
    }
    grDevices::png(
      png_path, width = 1800, height = 1200, res = 180
    )
    draw_zplot()
    grDevices::dev.off()
    grDevices::cairo_pdf(pdf_path, width = 10, height = 6.7)
    draw_zplot()
    grDevices::dev.off()
    written <- c(written, png_path, pdf_path)
  }

  if (any(!file.exists(written)) || any(file.info(written)$size <= 0)) {
    stop("Publication-bias reporting bundle is incomplete.",
         call. = FALSE)
  }
  normalizePath(written, winslash = "/", mustWork = TRUE)
}
