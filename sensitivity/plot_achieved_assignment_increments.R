#!/usr/bin/env Rscript

# Compare finalized achieved- and assignment-derived protein exposure curves
# as posterior increments from 0.8 g/kg/day. Reads only verified reductions.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

outcome_names <- c(
  whole_body = "Whole-body lean soft tissue",
  hypertrophy = "Direct muscle mass / size",
  strength = "Strength",
  `function` = "Physical function"
)

read_increments <- function(family, label) {
  read.csv(gzfile(file.path(
    "results", "verification", family,
    "dose_response_increments.csv.gz"
  ))) |>
    transmute(
      outcome,
      resistance_exercise_code,
      resistance_exercise,
      protein_intake,
      estimate = increment,
      conf.low = increment_conf_low,
      conf.high = increment_conf_high,
      within_group_observed_support,
      exposure_definition = label
    )
}

plot_data <- bind_rows(
  read_increments("achieved", "Achieved intake"),
  read_increments("assignment", "Assignment-derived intake")
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
    panel = interaction(
      outcome,
      resistance_exercise,
      sep = "  |  ",
      lex.order = TRUE
    ),
    exposure_definition = factor(
      exposure_definition,
      levels = c("Achieved intake", "Assignment-derived intake")
    )
  )

colours <- c(
  `Achieved intake` = "#2B6CB0",
  `Assignment-derived intake` = "#D97706"
)

comparison_plot <- ggplot(
  plot_data,
  aes(
    x = protein_intake,
    y = estimate,
    colour = exposure_definition,
    fill = exposure_definition,
    group = exposure_definition
  )
) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.35) +
  geom_vline(
    xintercept = 0.8,
    colour = "grey50",
    linewidth = 0.4,
    linetype = "dotted"
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 0.95) +
  facet_wrap(vars(panel), ncol = 2, scales = "free") +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(
    breaks = scales::breaks_width(0.4),
    expand = expansion(mult = c(0.015, 0.025))
  ) +
  labs(
    title = "Achieved versus assignment-derived protein dose response",
    subtitle = paste(
      "Posterior change from 0.8 g/kg/day; shaded regions are 95%",
      "credible intervals"
    ),
    x = "Protein intake (g/kg/day)",
    y = "Change in standardised outcome from 0.8 g/kg/day",
    colour = NULL,
    fill = NULL,
    caption = stringr::str_wrap(
      paste(
        "Both analyses use identical common-sample rows, candidate forms,",
        "priors, exact-corrected integrated LOO, and outcome-specific shared",
        "prediction grids. Curve tails should be interpreted with reference",
        "to the exposure-support tables in the verification bundles."
      ),
      width = 135
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
  "achieved_assignment_increment_overlay.png"
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
