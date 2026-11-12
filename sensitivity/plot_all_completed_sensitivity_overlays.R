#!/usr/bin/env Rscript

# Compare the completed model-averaged primary and sensitivity dose-response
# curves. This script only reads retained targets and verification artifacts;
# it never fits or recompiles a model.

suppressPackageStartupMessages({
  library(targets)
  library(dplyr)
  library(ggplot2)
})

outcome_names <- c(
  whole_body = "Whole-body lean soft tissue",
  hypertrophy = "Direct muscle mass / size",
  strength = "Strength",
  `function` = "Physical function"
)

primary <- bind_rows(
  whole_body = tar_read(preds_new_whole_body),
  hypertrophy = tar_read(preds_new_hypertrophy),
  strength = tar_read(preds_new_strength),
  `function` = tar_read(preds_new_function),
  .id = "outcome"
) |>
  transmute(
    outcome,
    resistance_exercise,
    protein_intake,
    estimate,
    conf.low,
    conf.high,
    analysis = "Primary"
  )

read_sensitivity_predictions <- function(family, label) {
  read.csv(
    gzfile(file.path(
      "results", "verification", family,
      "dose_response_predictions.csv.gz"
    ))
  ) |>
    transmute(
      outcome,
      resistance_exercise,
      protein_intake,
      estimate,
      conf.low,
      conf.high,
      within_group_observed_support,
      within_group_central_90_support,
      analysis = label
    )
}

mundlak <- read_sensitivity_predictions("mundlak", "Mundlak")
primary_support <- mundlak |>
  group_by(outcome, resistance_exercise) |>
  summarise(
    observed_minimum = min(protein_intake[within_group_observed_support]),
    observed_maximum = max(protein_intake[within_group_observed_support]),
    central_90_minimum = min(
      protein_intake[within_group_central_90_support]
    ),
    central_90_maximum = max(
      protein_intake[within_group_central_90_support]
    ),
    .groups = "drop"
  )

primary <- primary |>
  left_join(
    primary_support,
    by = c("outcome", "resistance_exercise")
  ) |>
  mutate(
    within_group_observed_support =
      protein_intake >= observed_minimum &
      protein_intake <= observed_maximum,
    within_group_central_90_support =
      protein_intake >= central_90_minimum &
      protein_intake <= central_90_maximum
  ) |>
  select(-ends_with("minimum"), -ends_with("maximum"))

plot_data <- bind_rows(
  primary,
  mundlak,
  read_sensitivity_predictions("within", "Within-study contrast"),
  read_sensitivity_predictions("achieved", "Achieved intake"),
  read_sensitivity_predictions("assignment", "Assignment-derived intake"),
  read_sensitivity_predictions("eiv", "Errors-in-variables")
) |>
  mutate(
    outcome = factor(
      outcome,
      levels = names(outcome_names),
      labels = unname(outcome_names)
    ),
    resistance_exercise = factor(
      resistance_exercise,
      levels = c("YES", "NO"),
      labels = c("Resistance exercise: Yes", "Resistance exercise: No")
    ),
    analysis = factor(
      analysis,
      levels = c(
        "Primary",
        "Mundlak",
        "Within-study contrast",
        "Achieved intake",
        "Assignment-derived intake",
        "Errors-in-variables"
      )
    )
  )

analysis_colours <- c(
  "Primary" = "#000000",
  "Mundlak" = "#0072B2",
  "Within-study contrast" = "#CC79A7",
  "Achieved intake" = "#009E73",
  "Assignment-derived intake" = "#E69F00",
  "Errors-in-variables" = "#D55E00"
)

analysis_linetypes <- c(
  "Primary" = "solid",
  "Mundlak" = "longdash",
  "Within-study contrast" = "dotdash",
  "Achieved intake" = "dashed",
  "Assignment-derived intake" = "twodash",
  "Errors-in-variables" = "dotted"
)

comparison_plot <- ggplot(
  plot_data,
  aes(
    x = protein_intake,
    y = estimate,
    colour = analysis,
    linetype = analysis,
    group = analysis
  )
) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_vline(
    xintercept = 0.8,
    colour = "grey58",
    linewidth = 0.35,
    linetype = "dotted"
  ) +
  geom_ribbon(
    data = filter(plot_data, analysis == "Primary"),
    aes(ymin = conf.low, ymax = conf.high, fill = analysis),
    inherit.aes = TRUE,
    alpha = 0.025,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    aes(linewidth = analysis == "Primary"),
    alpha = 0.16,
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = filter(
      plot_data,
      analysis == "Primary" & within_group_observed_support
    ),
    aes(ymin = conf.low, ymax = conf.high, fill = analysis),
    inherit.aes = TRUE,
    alpha = 0.045,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = filter(plot_data, within_group_observed_support),
    aes(linewidth = analysis == "Primary"),
    alpha = 0.52,
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = filter(
      plot_data,
      analysis == "Primary" & within_group_central_90_support
    ),
    aes(ymin = conf.low, ymax = conf.high, fill = analysis),
    inherit.aes = TRUE,
    alpha = 0.085,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = filter(plot_data, within_group_central_90_support),
    aes(linewidth = analysis == "Primary"),
    alpha = 0.96
  ) +
  facet_grid(
    rows = vars(outcome),
    cols = vars(resistance_exercise),
    scales = "free_y"
  ) +
  scale_colour_manual(values = analysis_colours) +
  scale_fill_manual(values = analysis_colours) +
  scale_linetype_manual(values = analysis_linetypes) +
  scale_linewidth_manual(values = c(`TRUE` = 1.05, `FALSE` = 0.78), guide = "none") +
  scale_x_continuous(
    breaks = scales::breaks_width(0.4),
    expand = expansion(mult = c(0.015, 0.025))
  ) +
  labs(
    title = "Primary and completed sensitivity dose-response curves",
    subtitle = paste(
      "Posterior means are shown for every analysis; darker sections have",
      "stronger intake support. The shading is the primary 95% credible interval."
    ),
    x = "Protein intake (g/kg/day)",
    y = "Standardised mean change",
    colour = NULL,
    linetype = NULL,
    caption = paste(
      "The dotted vertical reference marks 0.8 g/kg/day.",
      "Dark = central 90% support; medium = observed tails; faint = extrapolation.",
      "Within-study contrast estimates target a distinct contrast-based estimand;",
      "achieved and assignment-derived analyses use their common reduced sample."
    )
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "grey30", margin = margin(b = 8)),
    plot.caption = element_text(colour = "grey35", hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey91", linewidth = 0.3),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.spacing = grid::unit(0.85, "lines"),
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = grid::unit(1.25, "cm"),
    plot.margin = margin(12, 16, 12, 12)
  )

output_directory <- file.path("results", "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(
  output_directory,
  "primary_completed_sensitivity_overlay_data.csv.gz"
)
png_path <- file.path(
  output_directory,
  "primary_completed_sensitivity_dose_response_overlay.png"
)
pdf_path <- file.path(
  output_directory,
  "primary_completed_sensitivity_dose_response_overlay.pdf"
)

write.csv(plot_data, gzfile(data_path), row.names = FALSE)
ggsave(
  filename = png_path,
  plot = comparison_plot,
  width = 14,
  height = 15,
  units = "in",
  dpi = 220,
  bg = "white"
)
ggsave(
  filename = pdf_path,
  plot = comparison_plot,
  width = 14,
  height = 15,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

cat(normalizePath(png_path, winslash = "/", mustWork = TRUE), "\n")
