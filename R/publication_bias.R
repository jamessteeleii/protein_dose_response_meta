# Publication-bias reporting diagnostics ----
#
# This module deliberately stays upstream of every fitted model. It constructs
# observed protein-arm minus comparator-arm contrasts, selects one result per
# study and outcome using a fully deterministic hierarchy, and draws
# contour-enhanced funnel plots. A fitted-model reference line can be supplied
# later through `model_estimates` without changing the data preparation.

publication_bias_outcome_labels <- function() {
  c(
    "whole body lean mass" = "Whole-body lean/fat-free mass",
    "direct muscle size/mass" = "Measures of muscle size or mass",
    "strength" = "Strength",
    "function" = "Physical function"
  )
}

make_publication_bias_within_study_contrasts <- function(data_by_outcome) {
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
      "Publication-bias input is missing columns: ",
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

  protein <- data |>
    dplyr::filter(.data$condition == "protein") |>
    dplyr::transmute(
      study_no,
      outcome_group,
      measurement_method,
      unit_of_measurement,
      time1,
      time2,
      resistance_exercise_code,
      authors = as.character(.data$authors),
      protein_arm_no = .data$arm_no,
      protein_arm_no_extracted = .data$arm_no_extracted,
      yi_protein = .data$yi,
      vi_protein = .data$vi,
      n_protein = .data$n,
      protein_intake_treatment = .data$protein_intake,
      protein_intake_raw_treatment = .data$protein_intake_raw,
      update_2026_data_source = as.character(.data$update_2026_data_source)
    )

  comparator <- data |>
    dplyr::filter(.data$condition %in% c("control", "placebo")) |>
    dplyr::transmute(
      study_no,
      outcome_group,
      measurement_method,
      unit_of_measurement,
      time1,
      time2,
      resistance_exercise_code,
      comparator_arm_no = .data$arm_no,
      comparator_arm_no_extracted = .data$arm_no_extracted,
      comparator_condition = as.character(.data$condition),
      yi_comparator = .data$yi,
      vi_comparator = .data$vi,
      n_comparator = .data$n,
      protein_intake_comparator = .data$protein_intake,
      protein_intake_raw_comparator = .data$protein_intake_raw
    )

  paired <- dplyr::inner_join(
    protein,
    comparator,
    by = keys,
    relationship = "many-to-many"
  ) |>
    dplyr::filter(
      .data$protein_arm_no != .data$comparator_arm_no,
      is.finite(.data$yi_protein),
      is.finite(.data$vi_protein),
      .data$vi_protein > 0,
      is.finite(.data$yi_comparator),
      is.finite(.data$vi_comparator),
      .data$vi_comparator > 0,
      is.finite(.data$n_protein),
      .data$n_protein > 0,
      is.finite(.data$n_comparator),
      .data$n_comparator > 0
    ) |>
    dplyr::mutate(
      yi = .data$yi_protein - .data$yi_comparator,
      vi = .data$vi_protein + .data$vi_comparator,
      sei = sqrt(.data$vi),
      z_statistic = .data$yi / .data$sei,
      p_value_two_sided = 2 * stats::pnorm(
        abs(.data$z_statistic),
        lower.tail = FALSE
      ),
      total_n = .data$n_protein + .data$n_comparator,
      protein_intake_difference =
        .data$protein_intake_treatment -
        .data$protein_intake_comparator,
      protein_intake_raw_difference =
        .data$protein_intake_raw_treatment -
        .data$protein_intake_raw_comparator,
      resistance_exercise_label = dplyr::case_when(
        .data$resistance_exercise_code == 0 ~ "With resistance training",
        .data$resistance_exercise_code == 1 ~ "Without resistance training",
        TRUE ~ "Resistance-training status unavailable"
      ),
      effect_id = paste(
        .data$study_no,
        .data$outcome_group,
        .data$measurement_method,
        .data$unit_of_measurement,
        .data$time1,
        .data$time2,
        .data$protein_arm_no_extracted,
        .data$comparator_arm_no_extracted,
        sep = "::"
      )
    ) |>
    dplyr::arrange(
      factor(
        .data$outcome_group,
        levels = names(publication_bias_outcome_labels())
      ),
      .data$study_no,
      .data$effect_id
    )

  if (!nrow(paired)) {
    stop("No valid within-study funnel-plot contrasts were formed.",
         call. = FALSE)
  }
  if (
    anyDuplicated(paired$effect_id) ||
      any(!is.finite(paired$yi)) ||
      any(!is.finite(paired$vi)) ||
      any(paired$vi <= 0)
  ) {
    stop("Within-study funnel-plot contrasts failed validation.",
         call. = FALSE)
  }

  paired
}

select_publication_bias_index_contrasts <- function(contrasts) {
  required <- c(
    "study_no", "outcome_group", "time2", "total_n",
    "measurement_method", "unit_of_measurement",
    "protein_arm_no_extracted", "comparator_arm_no_extracted",
    "effect_id", "yi", "vi", "sei", "p_value_two_sided",
    "resistance_exercise_label"
  )
  missing <- setdiff(required, names(contrasts))
  if (length(missing)) {
    stop(
      "Index-contrast selection is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  audit <- contrasts |>
    dplyr::mutate(
      measurement_method_order = dplyr::coalesce(
        as.character(.data$measurement_method),
        ""
      ),
      unit_order = dplyr::coalesce(
        as.character(.data$unit_of_measurement),
        ""
      )
    ) |>
    dplyr::group_by(.data$study_no, .data$outcome_group) |>
    dplyr::arrange(
      dplyr::desc(.data$time2),
      dplyr::desc(.data$total_n),
      .data$measurement_method_order,
      .data$unit_order,
      .data$protein_arm_no_extracted,
      .data$comparator_arm_no_extracted,
      .data$effect_id,
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      selection_rank = dplyr::row_number(),
      candidate_count = dplyr::n(),
      selected_index_contrast = .data$selection_rank == 1L,
      selection_hierarchy = paste(
        "latest follow-up; largest combined sample;",
        "measurement, unit, arm, and effect identifiers"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      -dplyr::all_of(c("measurement_method_order", "unit_order"))
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
    dplyr::mutate(outcome_group = as.character(.data$outcome_group))

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
    stop("Index-contrast selection did not retain one result per study-outcome.",
         call. = FALSE)
  }

  summary <- selected |>
    dplyr::group_by(.data$outcome_group) |>
    dplyr::summarise(
      independent_studies = dplyr::n_distinct(.data$study_no),
      selected_contrasts = dplyr::n(),
      candidate_contrasts = sum(.data$candidate_count),
      median_total_n = stats::median(.data$total_n),
      minimum_total_n = min(.data$total_n),
      maximum_total_n = max(.data$total_n),
      statistically_significant = sum(.data$p_value_two_sided < 0.05),
      with_resistance_training = sum(
        .data$resistance_exercise_label == "With resistance training"
      ),
      without_resistance_training = sum(
        .data$resistance_exercise_label == "Without resistance training"
      ),
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

make_publication_bias_contour_grid <- function(x_limit, se_limit) {
  # Use cell midpoints so geom_raster fills the plotting region without
  # placing cell centres exactly on (and cell edges beyond) scale limits.
  x_breaks <- seq(-x_limit, x_limit, length.out = 321L)
  se_breaks <- seq(0, se_limit, length.out = 241L)
  x <- head(x_breaks, -1L) + diff(x_breaks) / 2
  se <- head(se_breaks, -1L) + diff(se_breaks) / 2
  grid <- tidyr::expand_grid(effect = x, standard_error = se) |>
    dplyr::mutate(
      p_value = 2 * stats::pnorm(
        abs(.data$effect / .data$standard_error),
        lower.tail = FALSE
      ),
      significance_region = cut(
        .data$p_value,
        breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
        labels = c(
          "p < 0.01", "0.01 <= p < 0.05",
          "0.05 <= p < 0.10", "p >= 0.10"
        ),
        right = FALSE
      )
    )
  grid
}

make_publication_bias_funnel_panel <- function(
    data,
    outcome_group,
    model_estimates = NULL
) {
  panel_data <- data |>
    dplyr::filter(.data$outcome_group == .env$outcome_group)
  if (!nrow(panel_data)) {
    stop("No funnel-plot data for outcome: ", outcome_group,
         call. = FALSE)
  }

  se_limit <- max(panel_data$sei) * 1.08
  x_limit <- max(
    abs(panel_data$yi),
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

  outcome_title <- unname(publication_bias_outcome_labels()[outcome_group])
  plot <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = contour_grid,
      ggplot2::aes(
        x = .data$effect,
        y = .data$standard_error,
        fill = .data$significance_region
      ),
      interpolate = FALSE
    ) +
    ggplot2::geom_line(
      data = contour_lines,
      ggplot2::aes(
        x = .data$effect,
        y = .data$standard_error,
        group = .data$line_group,
        linetype = .data$contour
      ),
      colour = "grey55",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = "grey20",
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      data = panel_data,
      ggplot2::aes(
        x = .data$yi,
        y = .data$sei,
        colour = .data$resistance_exercise_label
      ),
      size = 2.4,
      alpha = 0.88
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "p < 0.01" = "#D5E3F0",
        "0.01 <= p < 0.05" = "#E5ECF3",
        "0.05 <= p < 0.10" = "#F0F3F6",
        "p >= 0.10" = "#FFFFFF"
      ),
      drop = FALSE,
      name = "Two-sided null contour"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "With resistance training" = "#0072B2",
        "Without resistance training" = "#D55E00",
        "Resistance-training status unavailable" = "#666666"
      ),
      name = "Trial condition"
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "p = 0.10" = "dotted",
        "p = 0.05" = "dashed",
        "p = 0.01" = "dotdash"
      ),
      name = "Contour boundary"
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
      title = outcome_title,
      subtitle = paste0(
        "One index contrast per study (k = ",
        dplyr::n_distinct(panel_data$study_no),
        ")"
      ),
      x = "Observed protein-minus-comparator effect (SMD)",
      y = "Standard error"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      legend.key.height = grid::unit(0.7, "lines")
    )

  if (!is.null(model_estimates) && nrow(model_estimates)) {
    required_estimates <- c("outcome_group", "estimate")
    if (!all(required_estimates %in% names(model_estimates))) {
      stop(
        "Model estimates require outcome_group and estimate columns.",
        call. = FALSE
      )
    }
    estimate <- model_estimates |>
      dplyr::filter(.data$outcome_group == .env$outcome_group) |>
      dplyr::pull(.data$estimate)
    if (length(estimate) > 1L || any(!is.finite(estimate))) {
      stop("Each funnel panel accepts at most one finite model estimate.",
           call. = FALSE)
    }
    if (length(estimate) == 1L) {
      plot <- plot + ggplot2::geom_vline(
        xintercept = estimate,
        colour = "black",
        linetype = "longdash",
        linewidth = 0.7
      )
    }
  }

  plot
}

make_publication_bias_funnel_plot <- function(
    funnel_data,
    model_estimates = NULL
) {
  if (
    !is.list(funnel_data) ||
      is.null(funnel_data$index_contrasts) ||
      !nrow(funnel_data$index_contrasts)
  ) {
    stop("Funnel plotting requires selected index contrasts.",
         call. = FALSE)
  }

  outcomes <- names(publication_bias_outcome_labels())
  plots <- purrr::map(
    outcomes,
    ~ make_publication_bias_funnel_panel(
      data = funnel_data$index_contrasts,
      outcome_group = .x,
      model_estimates = model_estimates
    )
  )

  patchwork::wrap_plots(plots, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = "Contour-enhanced funnel plots of within-study contrasts",
      subtitle = "Shaded regions mark two-sided significance contours around the null."
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left"
    )
}

write_publication_bias_funnel_data <- function(
    funnel_data,
    output_dir = here::here(
      "results", "update_2026", "publication_bias"
    )
) {
  required <- c("candidate_audit", "index_contrasts", "outcome_summary")
  if (!is.list(funnel_data) || !all(required %in% names(funnel_data))) {
    stop("Publication-bias reporting data are incomplete.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    candidate_audit = file.path(output_dir, "candidate_contrasts.csv"),
    index_contrasts = file.path(output_dir, "index_contrasts.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv")
  )
  purrr::iwalk(paths, function(path, name) {
    readr::write_csv(funnel_data[[name]], path)
  })
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    stop("Publication-bias data files were not written correctly.",
         call. = FALSE)
  }
  normalizePath(unname(paths), winslash = "/", mustWork = TRUE)
}

write_publication_bias_funnel_plot <- function(
    plot,
    output_dir = here::here("plots", "update_2026", "publication_bias")
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    file.path(output_dir, "contour-enhanced-funnel-plots.png"),
    file.path(output_dir, "contour-enhanced-funnel-plots.tiff"),
    file.path(output_dir, "contour-enhanced-funnel-plots.pdf")
  )
  ggplot2::ggsave(
    paths[[1L]],
    plot,
    device = "png",
    dpi = 180,
    width = 11,
    height = 9
  )
  ggplot2::ggsave(
    paths[[2L]],
    plot,
    device = "tiff",
    dpi = 300,
    width = 11,
    height = 9,
    compression = "lzw"
  )
  ggplot2::ggsave(
    paths[[3L]],
    plot,
    device = grDevices::cairo_pdf,
    width = 11,
    height = 9
  )
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    stop("Publication-bias funnel plots were not written correctly.",
         call. = FALSE)
  }
  normalizePath(paths, winslash = "/", mustWork = TRUE)
}
