#!/usr/bin/env Rscript

# Plot finalized model-averaged primary and Mundlak dose-response summaries.
# This script reads existing targets/reporting artifacts and never fits models.

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

mundlak <- read.csv(
  gzfile("results/verification/mundlak/dose_response_predictions.csv.gz")
) |>
  transmute(
    outcome,
    resistance_exercise,
    protein_intake,
    estimate,
    conf.low,
    conf.high,
    analysis = "Mundlak"
  )

plot_data <- bind_rows(primary, mundlak) |>
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
    panel = interaction(
      outcome,
      resistance_exercise,
      sep = "  |  ",
      lex.order = TRUE
    ),
    analysis = factor(analysis, levels = c("Primary", "Mundlak"))
  )

analysis_colours <- c(
  Primary = "#2B6CB0",
  Mundlak = "#D97706"
)

comparison_plot <- ggplot(
  plot_data,
  aes(
    x = protein_intake,
    y = estimate,
    colour = analysis,
    fill = analysis,
    group = analysis
  )
) +
  geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.35) +
  geom_vline(
    xintercept = 0.8,
    colour = "grey55",
    linewidth = 0.4,
    linetype = "dotted"
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.13,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(vars(panel), ncol = 2, scales = "free") +
  scale_colour_manual(values = analysis_colours) +
  scale_fill_manual(values = analysis_colours) +
  scale_x_continuous(
    breaks = scales::breaks_width(0.4),
    expand = expansion(mult = c(0.015, 0.025))
  ) +
  labs(
    title = "Primary and Mundlak model-averaged dose-response curves",
    subtitle = paste(
      "Lines are posterior means; shaded regions are 95% credible intervals.",
      "The dotted line marks 0.8 g/kg/day."
    ),
    x = "Protein intake (g/kg/day)",
    y = "Standardised mean change",
    colour = NULL,
    fill = NULL,
    caption = paste(
      "Primary curves combine within- and between-study information;",
      "Mundlak curves separate study-level mean intake from within-study intake."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "grey30", margin = margin(b = 8)),
    plot.caption = element_text(colour = "grey35", hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = grid::unit(1.3, "cm"),
    plot.margin = margin(12, 16, 12, 12)
  )

output_directory <- file.path("results", "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_directory,
  "primary_mundlak_dose_response_overlay.png"
)

ggsave(
  filename = output_file,
  plot = comparison_plot,
  width = 12,
  height = 14,
  units = "in",
  dpi = 220,
  bg = "white"
)

cat(normalizePath(output_file, winslash = "/", mustWork = TRUE), "\n")
