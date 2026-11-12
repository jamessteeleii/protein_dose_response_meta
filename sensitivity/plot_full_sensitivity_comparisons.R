#!/usr/bin/env Rscript

# Publication-facing comparisons of the primary model average with every
# completed full sensitivity family. The script reads cached primary summaries
# and durable sensitivity verification bundles only; it never fits a model.

suppressPackageStartupMessages({
  library(targets)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

reference_intake <- 0.8
probability_threshold <- 0.95
comparison_root <- file.path(
  "results", "supplementary", "sensitivity_comparisons"
)
figure_root <- file.path("results", "figures", "sensitivity_comparisons")
dir.create(comparison_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

outcome_labels <- c(
  whole_body = "Whole-body lean soft tissue",
  hypertrophy = "Direct muscle mass / size",
  strength = "Strength",
  `function` = "Physical function"
)

family_labels <- c(
  mundlak = "Mundlak decomposition",
  within = "Strict within-study contrast",
  achieved = "Achieved intake",
  assignment = "Assignment-derived intake",
  eiv = "Errors-in-variables",
  rescor = "Residual correlation",
  strategy = "Protein-intake strategy (arm model)",
  strategy_within = "Protein-intake strategy (within-study model)"
)

# The two strategy bundles each contain separate Addition and Substitution
# curves. Treat those as explicit comparison series rather than silently
# stacking two rows under one family label.
comparison_labels <- c(
  family_labels[c("mundlak", "within", "achieved", "assignment", "eiv", "rescor")],
  strategy_addition = "Strategy arm model: Addition",
  strategy_substitution = "Strategy arm model: Substitution",
  strategy_within_addition = "Strategy within-study: Addition",
  strategy_within_substitution = "Strategy within-study: Substitution"
)

analysis_colours <- c(
  Primary = "#000000",
  `Mundlak decomposition` = "#0072B2",
  `Strict within-study contrast` = "#D55E00",
  `Achieved intake` = "#009E73",
  `Assignment-derived intake` = "#CC79A7",
  `Errors-in-variables` = "#E69F00",
  `Residual correlation` = "#56B4E9",
  `Strategy arm model: Addition` = "#332288",
  `Strategy arm model: Substitution` = "#88CCEE",
  `Strategy within-study: Addition` = "#882255",
  `Strategy within-study: Substitution` = "#44AA99"
)

analysis_linetypes <- c(
  Primary = "solid",
  `Mundlak decomposition` = "longdash",
  `Strict within-study contrast` = "dotdash",
  `Achieved intake` = "dashed",
  `Assignment-derived intake` = "twodash",
  `Errors-in-variables` = "dotted",
  `Residual correlation` = "solid",
  `Strategy arm model: Addition` = "longdash",
  `Strategy arm model: Substitution` = "dashed",
  `Strategy within-study: Addition` = "dotdash",
  `Strategy within-study: Substitution` = "twodash"
)

verification_root <- file.path("results", "verification")
verification_directories <- list.dirs(
  verification_root, recursive = FALSE, full.names = FALSE
)
completed_families <- verification_directories[
  file.exists(file.path(
    verification_root, verification_directories,
    "dose_response_predictions.csv.gz"
  )) &
    file.exists(file.path(
      verification_root, verification_directories,
      "pairwise_contrasts.csv.gz"
    ))
]
completed_families <- intersect(names(family_labels), completed_families)
completed_families <- names(family_labels)[
  names(family_labels) %in% completed_families
]
if (length(completed_families) < 1L) {
  stop("No completed sensitivity verification bundles were found.",
       call. = FALSE)
}

add_comparison_identity <- function(data, source_family) {
  if (source_family %in% c("strategy", "strategy_within")) {
    if (!"protein_strategy" %in% names(data)) {
      stop(
        "A protein-strategy table is missing `protein_strategy`.",
        call. = FALSE
      )
    }
    comparison_family <- paste0(
      source_family, "_", tolower(data$protein_strategy)
    )
  } else {
    comparison_family <- rep(source_family, nrow(data))
  }
  data |>
    dplyr::mutate(
      source_family = .env$source_family,
      family = comparison_family,
      analysis = unname(comparison_labels[comparison_family])
    )
}

completed_comparison_ids <- c(
  setdiff(completed_families, c("strategy", "strategy_within")),
  if ("strategy" %in% completed_families) {
    c("strategy_addition", "strategy_substitution")
  },
  if ("strategy_within" %in% completed_families) {
    c("strategy_within_addition", "strategy_within_substitution")
  }
)
completed_comparison_ids <- names(comparison_labels)[
  names(comparison_labels) %in% completed_comparison_ids
]
comparison_registry <- tibble::tibble(
  family = completed_comparison_ids,
  source_family = dplyr::case_when(
    grepl("^strategy_within_", family) ~ "strategy_within",
    grepl("^strategy_", family) ~ "strategy",
    TRUE ~ family
  ),
  analysis = unname(comparison_labels[family])
)

read_verification_table <- function(family, filename) {
  readr::read_csv(
    file.path(verification_root, family, filename),
    show_col_types = FALSE,
    progress = FALSE
  )
}

primary_prediction_targets <- c(
  whole_body = "preds_new_whole_body",
  hypertrophy = "preds_new_hypertrophy",
  strength = "preds_new_strength",
  `function` = "preds_new_function"
)
primary_contrast_targets <- c(
  whole_body = "primary_heatmap_data_whole_body",
  hypertrophy = "primary_heatmap_data_hypertrophy",
  strength = "primary_heatmap_data_strength",
  `function` = "primary_heatmap_data_function"
)

primary_predictions <- purrr::imap_dfr(
  primary_prediction_targets,
  ~ targets::tar_read_raw(.x) |>
    dplyr::mutate(outcome = .y, .before = 1L)
) |>
  dplyr::transmute(
    outcome,
    resistance_exercise_code,
    resistance_exercise,
    protein_intake,
    estimate,
    conf_low = conf.low,
    conf_high = conf.high,
    family = "primary",
    analysis = "Primary"
  )

primary_contrasts <- purrr::imap_dfr(
  primary_contrast_targets,
  ~ targets::tar_read_raw(.x) |>
    dplyr::mutate(outcome = .y, .before = 1L)
) |>
  dplyr::mutate(family = "primary", analysis = "Primary")

primary_support <- primary_contrasts |>
  dplyr::distinct(
    outcome, resistance_exercise_code, observed_min, observed_max
  )
primary_predictions <- primary_predictions |>
  dplyr::left_join(
    primary_support,
    by = c("outcome", "resistance_exercise_code")
  ) |>
  dplyr::mutate(
    supported =
      protein_intake >= observed_min - 1e-10 &
      protein_intake <= observed_max + 1e-10
  ) |>
  dplyr::select(-observed_min, -observed_max)

sensitivity_predictions <- purrr::map_dfr(
  completed_families,
  function(family) {
    data <- read_verification_table(
      family, "dose_response_predictions.csv.gz"
    )
    if ("within_group_observed_support" %in% names(data)) {
      data$supported_for_comparison <- data$within_group_observed_support
    } else if ("within_strategy_observed_support" %in% names(data)) {
      data$supported_for_comparison <- data$within_strategy_observed_support
    } else {
      data$supported_for_comparison <- TRUE
    }
    add_comparison_identity(data, family) |>
      dplyr::transmute(
        outcome,
        resistance_exercise_code,
        resistance_exercise,
        protein_intake,
        estimate,
        conf_low = conf.low,
        conf_high = conf.high,
        source_family,
        family,
        analysis,
        supported = supported_for_comparison
      )
  }
)

all_predictions <- dplyr::bind_rows(
  primary_predictions,
  sensitivity_predictions
)

primary_model_average_targets <- c(
  whole_body = "primary_model_average_whole_body",
  hypertrophy = "primary_model_average_hypertrophy",
  strength = "primary_model_average_strength",
  `function` = "primary_model_average_function"
)

summarise_primary_reference_increments <- function(
    model_average, outcome
) {
  grid <- model_average$grid
  draws <- as.matrix(model_average$draws)
  increment_draws <- draws
  for (exercise_code in c(0, 1)) {
    columns <- which(grid$resistance_exercise_code == exercise_code)
    reference_column <- columns[
      abs(grid$protein_intake[columns] - reference_intake) < 1e-10
    ]
    if (length(reference_column) != 1L) {
      stop("The primary 0.8 g/kg/day reference is not unique.",
           call. = FALSE)
    }
    increment_draws[, columns] <- sweep(
      draws[, columns, drop = FALSE],
      1L,
      draws[, reference_column],
      FUN = "-"
    )
  }
  intervals <- matrixStats::colQuantiles(
    increment_draws,
    probs = c(0.025, 0.975),
    drop = FALSE
  )
  dplyr::bind_cols(
    grid,
    tibble::tibble(
      estimate = colMeans(increment_draws),
      conf_low = intervals[, 1L],
      conf_high = intervals[, 2L]
    )
  ) |>
    dplyr::mutate(outcome = outcome, .before = 1L) |>
    dplyr::left_join(
      primary_support,
      by = c("outcome", "resistance_exercise_code")
    ) |>
    dplyr::transmute(
      outcome,
      resistance_exercise_code,
      resistance_exercise,
      protein_intake,
      estimate,
      conf_low,
      conf_high,
      family = "primary",
      analysis = "Primary",
      supported =
        protein_intake >= observed_min - 1e-10 &
        protein_intake <= observed_max + 1e-10
    )
}

primary_increments <- purrr::imap_dfr(
  primary_model_average_targets,
  ~ summarise_primary_reference_increments(
    targets::tar_read_raw(.x), .y
  )
)

sensitivity_contrast_tables <- purrr::set_names(
  purrr::map(
    completed_families,
    ~ read_verification_table(.x, "pairwise_contrasts.csv.gz")
  ),
  completed_families
)

sensitivity_increments <- purrr::map_dfr(
  completed_families,
  function(family) {
    data <- read_verification_table(
      family, "dose_response_increments.csv.gz"
    )
    if (!"within_group_observed_support" %in% names(data)) {
      data$within_group_observed_support <- TRUE
    }
    add_comparison_identity(data, family) |>
      dplyr::transmute(
        outcome,
        resistance_exercise_code,
        resistance_exercise,
        protein_intake,
        estimate = increment,
        conf_low = increment_conf_low,
        conf_high = increment_conf_high,
        source_family,
        family,
        analysis,
        supported = within_group_observed_support
      )
  }
)
all_increments <- dplyr::bind_rows(
  primary_increments,
  sensitivity_increments
)

analysis_levels <- c(
  "Primary", unname(comparison_labels[completed_comparison_ids])
)
format_plot_data <- function(data) {
  data |>
    dplyr::mutate(
      outcome = factor(
        outcome,
        levels = names(outcome_labels),
        labels = unname(outcome_labels)
      ),
      resistance_exercise = factor(
        resistance_exercise,
        levels = c("NO", "YES"),
        labels = c(
          "Resistance Training: No",
          "Resistance Training: Yes"
        )
      ),
      analysis = factor(analysis, levels = analysis_levels)
    )
}

all_predictions_plot <- format_plot_data(all_predictions)
all_increments_plot <- format_plot_data(all_increments)

primary_increment_comparison <- primary_increments |>
  dplyr::mutate(protein_key = round(protein_intake, 8L)) |>
  dplyr::select(
    outcome,
    resistance_exercise_code,
    protein_key,
    protein_intake,
    primary_estimate = estimate,
    primary_conf_low = conf_low,
    primary_conf_high = conf_high,
    primary_supported = supported
  )

interval_overlap_coefficient <- function(
    low_a, high_a, low_b, high_b
) {
  intersection <- pmax(
    0, pmin(high_a, high_b) - pmax(low_a, low_b)
  )
  smaller_width <- pmin(high_a - low_a, high_b - low_b)
  dplyr::if_else(
    smaller_width <= sqrt(.Machine$double.eps),
    as.numeric(
      abs(low_a - low_b) <= 1e-10 & abs(high_a - high_b) <= 1e-10
    ),
    pmin(1, intersection / smaller_width)
  )
}

increment_comparison <- sensitivity_increments |>
  dplyr::mutate(protein_key = round(protein_intake, 8L)) |>
  dplyr::inner_join(
    primary_increment_comparison,
    by = c("outcome", "resistance_exercise_code", "protein_key"),
    relationship = "many-to-one",
    suffix = c("_sensitivity", "_primary")
  ) |>
  dplyr::mutate(
    common_supported = supported & primary_supported,
    mean_difference = estimate - primary_estimate,
    absolute_mean_difference = abs(mean_difference),
    intervals_overlap =
      pmax(conf_low, primary_conf_low) <=
      pmin(conf_high, primary_conf_high),
    interval_overlap_coefficient = interval_overlap_coefficient(
      conf_low, conf_high, primary_conf_low, primary_conf_high
    )
  )

curve_agreement <- increment_comparison |>
  dplyr::filter(
    common_supported,
    abs(protein_intake_sensitivity - reference_intake) > 1e-10
  ) |>
  dplyr::group_by(family, analysis, outcome, resistance_exercise_code) |>
  dplyr::summarise(
    supported_grid_points = dplyr::n(),
    mean_absolute_difference = mean(absolute_mean_difference),
    root_mean_squared_difference = sqrt(mean(mean_difference^2)),
    maximum_absolute_difference = max(absolute_mean_difference),
    mean_interval_overlap_coefficient =
      mean(interval_overlap_coefficient),
    proportion_with_overlapping_95_intervals = mean(intervals_overlap),
    .groups = "drop"
  )

classify_contrast <- function(gain, negligible, loss) {
  dplyr::case_when(
    gain >= probability_threshold ~ "Meaningful benefit",
    negligible >= probability_threshold ~ "Practical equivalence",
    loss >= probability_threshold ~ "Meaningful detriment",
    TRUE ~ "Uncertain"
  )
}

primary_contrast_comparison <- primary_contrasts |>
  dplyr::mutate(
    current_key = round(protein_intake_current, 8L),
    proposed_key = round(protein_intake_proposed, 8L)
  ) |>
  dplyr::transmute(
    outcome,
    resistance_exercise_code,
    current_key,
    proposed_key,
    protein_intake_current,
    protein_intake_proposed,
    primary_estimate = estimate,
    primary_conf_low = conf_low,
    primary_conf_high = conf_high,
    primary_gain = probability_meaningful_gain,
    primary_negligible = probability_negligible,
    primary_loss = probability_meaningful_loss,
    primary_supported = supported
  )

contrast_comparison <- purrr::imap_dfr(
  sensitivity_contrast_tables,
  function(data, family) {
    add_comparison_identity(data, family) |>
      dplyr::mutate(
        current_key = round(protein_intake_current, 8L),
        proposed_key = round(protein_intake_proposed, 8L)
      ) |>
      dplyr::inner_join(
        primary_contrast_comparison,
        by = c(
          "outcome", "resistance_exercise_code",
          "current_key", "proposed_key"
        ),
        relationship = "many-to-one",
        suffix = c("_sensitivity", "_primary")
      ) |>
      dplyr::transmute(
        source_family,
        family,
        analysis,
        outcome,
        resistance_exercise_code,
        resistance_exercise,
        protein_intake_current = protein_intake_current_sensitivity,
        protein_intake_proposed = protein_intake_proposed_sensitivity,
        sensitivity_estimate = estimate,
        primary_estimate,
        mean_contrast_difference = estimate - primary_estimate,
        sensitivity_conf_low = conf_low,
        sensitivity_conf_high = conf_high,
        primary_conf_low,
        primary_conf_high,
        contrast_intervals_overlap =
          pmax(conf_low, primary_conf_low) <=
          pmin(conf_high, primary_conf_high),
        contrast_interval_overlap_coefficient =
          interval_overlap_coefficient(
            conf_low, conf_high, primary_conf_low, primary_conf_high
          ),
        sensitivity_gain = probability_meaningful_gain,
        sensitivity_negligible = probability_negligible,
        sensitivity_loss = probability_meaningful_loss,
        primary_gain,
        primary_negligible,
        primary_loss,
        sensitivity_supported = supported,
        primary_supported,
        common_supported = supported & primary_supported,
        sensitivity_classification = classify_contrast(
          probability_meaningful_gain,
          probability_negligible,
          probability_meaningful_loss
        ),
        primary_classification = classify_contrast(
          primary_gain, primary_negligible, primary_loss
        ),
        probability_total_variation = 0.5 * (
          abs(probability_meaningful_gain - primary_gain) +
            abs(probability_negligible - primary_negligible) +
            abs(probability_meaningful_loss - primary_loss)
        )
      )
  }
)

contrast_agreement <- contrast_comparison |>
  dplyr::filter(
    common_supported,
    protein_intake_proposed > protein_intake_current + 1e-10
  ) |>
  dplyr::group_by(family, analysis, outcome, resistance_exercise_code) |>
  dplyr::summarise(
    supported_unique_contrasts = dplyr::n(),
    mean_absolute_contrast_difference =
      mean(abs(mean_contrast_difference)),
    root_mean_squared_contrast_difference =
      sqrt(mean(mean_contrast_difference^2)),
    maximum_absolute_contrast_difference =
      max(abs(mean_contrast_difference)),
    mean_interval_overlap_coefficient =
      mean(contrast_interval_overlap_coefficient),
    proportion_with_overlapping_95_intervals =
      mean(contrast_intervals_overlap),
    classification_agreement = mean(
      sensitivity_classification == primary_classification
    ),
    mean_probability_total_variation =
      mean(probability_total_variation),
    maximum_probability_total_variation =
      max(probability_total_variation),
    .groups = "drop"
  )

readr::write_csv(
  all_predictions,
  file.path(comparison_root, "absolute_curve_comparison_data.csv.gz")
)
readr::write_csv(
  all_increments,
  file.path(comparison_root, "reference_increment_comparison_data.csv.gz")
)
readr::write_csv(
  curve_agreement,
  file.path(comparison_root, "curve_agreement_metrics.csv")
)
readr::write_csv(
  contrast_agreement,
  file.path(comparison_root, "contrast_agreement_metrics.csv")
)
readr::write_csv(
  comparison_registry |>
    dplyr::mutate(
    verification_directory = normalizePath(
      file.path(verification_root, source_family),
      winslash = "/",
      mustWork = TRUE
    )
    ),
  file.path(comparison_root, "included_sensitivity_families.csv")
)

comparison_theme <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9.5),
      panel.spacing = grid::unit(0.9, "lines"),
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = grid::unit(1.25, "cm")
    )
}

make_overlay <- function(data, y_label, title, subtitle, caption) {
  supported_data <- dplyr::filter(data, supported)
  primary_supported_data <- dplyr::filter(
    supported_data, analysis == "Primary"
  )
  ggplot2::ggplot(
    supported_data,
    ggplot2::aes(
      x = protein_intake,
      y = estimate,
      colour = analysis,
      fill = analysis,
      linetype = analysis,
      group = analysis
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dotted", colour = "grey55",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = reference_intake, linetype = "dotted",
      colour = "grey55", linewidth = 0.35
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf_low, ymax = conf_high),
      alpha = 0.035, colour = NA, show.legend = FALSE
    ) +
    ggplot2::geom_line(
      linewidth = 0.72, alpha = 0.9
    ) +
    ggplot2::geom_ribbon(
      data = primary_supported_data,
      ggplot2::aes(ymin = conf_low, ymax = conf_high),
      alpha = 0.14, colour = NA, show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = primary_supported_data, linewidth = 1.02,
      show.legend = FALSE
    ) +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(outcome),
      cols = ggplot2::vars(resistance_exercise),
      scales = "free",
      independent = "all"
    ) +
    ggplot2::scale_colour_manual(values = analysis_colours, drop = TRUE) +
    ggplot2::scale_fill_manual(values = analysis_colours, drop = TRUE) +
    ggplot2::scale_linetype_manual(
      values = analysis_linetypes, drop = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_width(0.4),
      expand = ggplot2::expansion(mult = c(0.015, 0.025))
    ) +
    ggplot2::labs(
      x = "Protein intake (g/kg/day)",
      y = y_label,
      title = title,
      subtitle = subtitle,
      caption = caption,
      colour = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
    ) +
    comparison_theme()
}

absolute_plot_data <- all_predictions_plot |>
  dplyr::filter(
    is.na(source_family) |
      !source_family %in% c("within", "strategy_within")
  )
absolute_plot <- make_overlay(
  absolute_plot_data,
  y_label = "Standardised mean change",
  title = "Absolute model-averaged dose-response predictions",
  subtitle = paste(
    "Lines are posterior means; shading gives 95% posterior quantile",
    "intervals within each analysis' observed intake support"
  ),
  caption = paste(
    "The strict within-study estimand is omitted because its absolute level is",
    "an arm contrast rather than an absolute outcome prediction.",
    "Curves are restricted to analysis-specific observed intake support."
  )
)

increment_plot <- make_overlay(
  all_increments_plot,
  y_label = "Change from 0.8 g/kg/day (standardised units)",
  title = "Sensitivity of the model-averaged dose-response shape",
  subtitle = paste(
    "Posterior change from 0.8 g/kg/day; lines are means and shading gives",
    "95% posterior quantile intervals"
  ),
  caption = paste(
    "Reference subtraction removes analysis-specific intercept differences.",
    "Proposed doses are restricted to observed support; in panels where 0.8",
    "g/kg/day lies outside support, the reference prediction is extrapolated."
  )
)

save_plot_pair <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(figure_root, paste0(stem, ".png")),
    plot = plot, width = width, height = height, units = "in",
    dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_root, paste0(stem, ".pdf")),
    plot = plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
}

save_plot_pair(
  absolute_plot, "absolute_dose_response_overlay", 14, 14.5
)
save_plot_pair(
  increment_plot, "reference_0.8_increment_overlay", 14, 14.5
)

# One page per sensitivity family: both marginal 95% intervals are visible
# without the occlusion inherent in an all-family overlay.
paired_pdf <- file.path(
  figure_root, "primary_vs_sensitivity_reference_increment_pages.pdf"
)
grDevices::cairo_pdf(paired_pdf, width = 12, height = 13.5, onefile = TRUE)
for (family in completed_comparison_ids) {
  pair_data <- all_increments_plot |>
    dplyr::filter(
      .data$family %in% c("primary", .env$family), supported
    ) |>
    dplyr::mutate(
      analysis = droplevels(analysis)
    )
  pair_plot <- ggplot2::ggplot(
    pair_data,
    ggplot2::aes(
      x = protein_intake, y = estimate,
      colour = analysis, fill = analysis, group = analysis
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dotted", colour = "grey55",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = reference_intake, linetype = "dotted",
      colour = "grey55", linewidth = 0.35
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf_low, ymax = conf_high),
      alpha = 0.16, colour = NA
    ) +
    ggplot2::geom_line(
      linewidth = 0.85
    ) +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(outcome),
      cols = ggplot2::vars(resistance_exercise),
      scales = "free",
      independent = "all"
    ) +
    ggplot2::scale_colour_manual(values = analysis_colours, drop = TRUE) +
    ggplot2::scale_fill_manual(values = analysis_colours, drop = TRUE) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_width(0.4),
      expand = ggplot2::expansion(mult = c(0.015, 0.025))
    ) +
    ggplot2::labs(
      x = "Protein intake (g/kg/day)",
      y = "Change from 0.8 g/kg/day (standardised units)",
      title = paste(
        "Primary analysis versus", unname(comparison_labels[[family]])
      ),
      subtitle = "Posterior means and marginal 95% posterior quantile intervals",
      caption = stringr::str_wrap(
        paste(
          "Intervals are marginal to each separately fitted analysis; their",
          "overlap is descriptive and is not a joint posterior probability.",
          "Proposed doses are restricted to observed support; the 0.8 reference",
          "may be extrapolated where it lies below the observed range."
        ),
        width = 125
      ),
      colour = NULL,
      fill = NULL
    ) +
    comparison_theme()
  print(pair_plot)
}
grDevices::dev.off()

# Full heatmap coverage: sensitivity-minus-primary posterior mean contrasts on
# every dose pair supported by both analyses. Pages share one colour scale.
heatmap_data <- contrast_comparison |>
  dplyr::filter(common_supported) |>
  dplyr::mutate(
    outcome = factor(
      outcome,
      levels = names(outcome_labels),
      labels = unname(outcome_labels)
    ),
    resistance_exercise = factor(
      resistance_exercise,
      levels = c("NO", "YES"),
      labels = c(
        "Resistance Training: No",
        "Resistance Training: Yes"
      )
    )
  )
heatmap_limit <- max(abs(heatmap_data$mean_contrast_difference), na.rm = TRUE)
if (!is.finite(heatmap_limit) || heatmap_limit <= 0) heatmap_limit <- 0.1
heatmap_pdf <- file.path(
  figure_root, "full_contrast_sensitivity_minus_primary_heatmaps.pdf"
)
grDevices::cairo_pdf(heatmap_pdf, width = 12, height = 13.5, onefile = TRUE)
for (family in completed_comparison_ids) {
  family_heatmap <- dplyr::filter(
    heatmap_data, .data$family == .env$family
  )
  heatmap_plot <- ggplot2::ggplot(
    family_heatmap,
    ggplot2::aes(
      x = protein_intake_current,
      y = protein_intake_proposed,
      fill = mean_contrast_difference
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, colour = "black", linewidth = 0.55
    ) +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(outcome),
      cols = ggplot2::vars(resistance_exercise),
      scales = "free",
      independent = "all"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#D55E00",
      mid = "white",
      high = "#0072B2",
      midpoint = 0,
      limits = c(-heatmap_limit, heatmap_limit),
      oob = scales::squish,
      name = "Sensitivity - primary\nmean contrast"
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_width(0.4),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_width(0.4),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::labs(
      x = "Current protein intake (g/kg/day)",
      y = "Proposed protein intake (g/kg/day)",
      title = paste(
        unname(comparison_labels[[family]]), "versus the primary analysis"
      ),
      subtitle = paste(
        "Difference in posterior mean pairwise contrasts over common",
        "observed dose-pair support"
      ),
      caption = paste(
        "Blue indicates a more positive contrast in the sensitivity analysis;",
        "orange indicates a less positive contrast."
      )
    ) +
    comparison_theme() +
    ggplot2::theme(aspect.ratio = 1)
  print(heatmap_plot)
}
grDevices::dev.off()

agreement_plot_data <- dplyr::bind_rows(
  contrast_agreement |>
    dplyr::transmute(
      family, analysis, outcome, resistance_exercise_code,
      metric = "Full-heatmap mean 95% interval overlap",
      agreement = mean_interval_overlap_coefficient
    ),
  contrast_agreement |>
    dplyr::transmute(
      family, analysis, outcome, resistance_exercise_code,
      metric = "Full-heatmap proportion of intervals overlapping",
      agreement = proportion_with_overlapping_95_intervals
    ),
  contrast_agreement |>
    dplyr::transmute(
      family, analysis, outcome, resistance_exercise_code,
      metric = "Full-heatmap classification agreement",
      agreement = classification_agreement
    ),
  contrast_agreement |>
    dplyr::transmute(
      family, analysis, outcome, resistance_exercise_code,
      metric = "Full-heatmap probability agreement (1 - TV)",
      agreement = 1 - mean_probability_total_variation
    )
) |>
  dplyr::mutate(
    analysis = factor(analysis, levels = analysis_levels[-1L]),
    outcome = factor(
      outcome,
      levels = names(outcome_labels),
      labels = unname(outcome_labels)
    ),
    resistance_exercise = factor(
      resistance_exercise_code,
      levels = c(1, 0),
      labels = c(
        "Resistance Training: No",
        "Resistance Training: Yes"
      )
    ),
    metric = factor(
      metric,
      levels = c(
        "Full-heatmap mean 95% interval overlap",
        "Full-heatmap proportion of intervals overlapping",
        "Full-heatmap classification agreement",
        "Full-heatmap probability agreement (1 - TV)"
      )
    )
  )

agreement_plot <- ggplot2::ggplot(
  agreement_plot_data,
  ggplot2::aes(
    x = analysis,
    y = agreement,
    colour = resistance_exercise,
    shape = resistance_exercise,
    group = resistance_exercise
  )
) +
  ggplot2::geom_hline(
    yintercept = 1, colour = "grey65", linewidth = 0.35
  ) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::geom_line(linewidth = 0.45, alpha = 0.7) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(metric),
    cols = ggplot2::vars(outcome)
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Resistance Training: No" = "#D55E00",
      "Resistance Training: Yes" = "#0072B2"
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Agreement with primary analysis",
    title = "Descriptive agreement of sensitivity analyses with the primary analysis",
    subtitle = paste(
      "All metrics use unique dose increases supported by both analyses",
      "across the full pairwise-contrast heatmap"
    ),
    caption = paste(
      "TV = total-variation distance across posterior probabilities of",
      "meaningful benefit, practical equivalence, and meaningful detriment.",
      "Interval overlap is a marginal descriptive measure, not a joint posterior probability."
    ),
    colour = NULL,
    shape = NULL
  ) +
  comparison_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 40, hjust = 1, vjust = 1, size = 7.5
    ),
    strip.text.y = ggplot2::element_text(angle = 0, hjust = 0)
  )

save_plot_pair(agreement_plot, "sensitivity_agreement_summary", 16, 12)

# Direct draw-wise strategy comparison. This is the scientific contrast that
# distinguishes the two strategy analyses from an ordinary robustness curve:
# positive values favour protein addition over macronutrient substitution at
# the same dose and resistance-training condition.
strategy_difference_data <- purrr::map_dfr(
  intersect(c("strategy", "strategy_within"), completed_families),
  function(family) {
    read_verification_table(
      family, "addition_minus_substitution.csv.gz"
    ) |>
      dplyr::transmute(
        source_family = .env$family,
        analysis = dplyr::if_else(
          .env$family == "strategy",
          "Arm-level strategy model",
          "Strict within-study strategy model"
        ),
        outcome,
        resistance_exercise_code,
        resistance_exercise,
        protein_intake,
        estimate = addition_minus_substitution,
        conf_low = conf.low,
        conf_high = conf.high,
        probability_positive,
        within_common_observed_support,
        within_common_central_90_support
      )
  }
)

readr::write_csv(
  strategy_difference_data,
  file.path(comparison_root, "protein_strategy_difference_data.csv.gz")
)

if (nrow(strategy_difference_data) > 0L) {
  strategy_difference_plot_data <- strategy_difference_data |>
    dplyr::filter(within_common_observed_support) |>
    dplyr::mutate(
      outcome = factor(
        outcome,
        levels = names(outcome_labels),
        labels = unname(outcome_labels)
      ),
      resistance_exercise = factor(
        resistance_exercise,
        levels = c("NO", "YES"),
        labels = c(
          "Resistance Training: No",
          "Resistance Training: Yes"
        )
      )
    )
  strategy_difference_plot <- ggplot2::ggplot(
    strategy_difference_plot_data,
    ggplot2::aes(
      x = protein_intake,
      y = estimate,
      colour = analysis,
      fill = analysis,
      group = analysis
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dotted", colour = "grey50",
      linewidth = 0.4
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf_low, ymax = conf_high),
      alpha = 0.16, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(outcome),
      cols = ggplot2::vars(resistance_exercise),
      scales = "free",
      independent = "all"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Arm-level strategy model" = "#332288",
        "Strict within-study strategy model" = "#882255"
      )
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Arm-level strategy model" = "#332288",
        "Strict within-study strategy model" = "#882255"
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_width(0.4),
      expand = ggplot2::expansion(mult = c(0.015, 0.025))
    ) +
    ggplot2::labs(
      x = "Protein intake (g/kg/day)",
      y = "Addition minus substitution (standardised units)",
      title = "Direct comparison of protein-intake strategies",
      subtitle = paste(
        "Posterior mean difference with marginal 95% posterior quantile",
        "intervals over common observed support"
      ),
      caption = paste(
        "Positive values indicate a more positive predicted outcome or",
        "within-study protein contrast under addition than substitution.",
        "The two lines represent different estimands and are not pooled."
      ),
      colour = NULL,
      fill = NULL
    ) +
    comparison_theme()
  save_plot_pair(
    strategy_difference_plot,
    "protein_strategy_addition_minus_substitution",
    14,
    14.5
  )
}

manifest <- tibble::tibble(
  artifact = c(
    "absolute_curve_data",
    "reference_increment_data",
    "curve_agreement_metrics",
    "contrast_agreement_metrics",
    "included_families",
    "protein_strategy_difference_data",
    "absolute_overlay_png",
    "absolute_overlay_pdf",
    "reference_increment_overlay_png",
    "reference_increment_overlay_pdf",
    "paired_increment_pages_pdf",
    "full_contrast_difference_heatmaps_pdf",
    "agreement_summary_png",
    "agreement_summary_pdf",
    "protein_strategy_difference_png",
    "protein_strategy_difference_pdf"
  ),
  path = c(
    file.path(comparison_root, "absolute_curve_comparison_data.csv.gz"),
    file.path(comparison_root, "reference_increment_comparison_data.csv.gz"),
    file.path(comparison_root, "curve_agreement_metrics.csv"),
    file.path(comparison_root, "contrast_agreement_metrics.csv"),
    file.path(comparison_root, "included_sensitivity_families.csv"),
    file.path(comparison_root, "protein_strategy_difference_data.csv.gz"),
    file.path(figure_root, "absolute_dose_response_overlay.png"),
    file.path(figure_root, "absolute_dose_response_overlay.pdf"),
    file.path(figure_root, "reference_0.8_increment_overlay.png"),
    file.path(figure_root, "reference_0.8_increment_overlay.pdf"),
    paired_pdf,
    heatmap_pdf,
    file.path(figure_root, "sensitivity_agreement_summary.png"),
    file.path(figure_root, "sensitivity_agreement_summary.pdf"),
    file.path(
      figure_root, "protein_strategy_addition_minus_substitution.png"
    ),
    file.path(
      figure_root, "protein_strategy_addition_minus_substitution.pdf"
    )
  )
) |>
  dplyr::mutate(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = vapply(
      path, digest::digest, character(1), file = TRUE, algo = "sha256"
    )
  )
readr::write_csv(
  manifest,
  file.path(comparison_root, "manifest.csv")
)

cat(
  normalizePath(comparison_root, winslash = "/", mustWork = TRUE),
  "\n"
)
