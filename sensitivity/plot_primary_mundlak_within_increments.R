#!/usr/bin/env Rscript

# Compare finalized primary, Mundlak, and strict within-study dose-response
# shapes as paired posterior increments from 0.8 g/kg/day. This script reads
# existing posterior mixtures/reporting artifacts and never fits models.

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

summarise_primary_increments <- function(model_average, outcome) {
  grid <- model_average$grid
  draws <- model_average$draws
  increment_draws <- draws
  for (exercise_code in c(0, 1)) {
    group_rows <- which(grid$resistance_exercise_code == exercise_code)
    reference_row <- group_rows[
      abs(grid$protein_intake[group_rows] - 0.8) < 1e-10
    ]
    if (length(reference_row) != 1L) {
      stop("The primary 0.8 g/kg/day reference is not unique.",
           call. = FALSE)
    }
    increment_draws[, group_rows] <- sweep(
      draws[, group_rows, drop = FALSE],
      1L,
      draws[, reference_row],
      FUN = "-"
    )
  }
  intervals <- matrixStats::colQuantiles(
    increment_draws,
    probs = c(0.025, 0.975),
    drop = FALSE
  )
  bind_cols(
    grid,
    tibble::tibble(
      estimate = colMeans(increment_draws),
      conf.low = intervals[, 1L],
      conf.high = intervals[, 2L]
    )
  ) |>
    transmute(
      outcome = outcome,
      resistance_exercise_code,
      resistance_exercise,
      protein_intake,
      estimate,
      conf.low,
      conf.high,
      analysis = "Primary"
    )
}

primary_targets <- c(
  whole_body = "primary_model_average_whole_body",
  hypertrophy = "primary_model_average_hypertrophy",
  strength = "primary_model_average_strength",
  `function` = "primary_model_average_function"
)

primary <- bind_rows(lapply(names(primary_targets), function(outcome) {
  summarise_primary_increments(
    tar_read_raw(primary_targets[[outcome]]),
    outcome
  )
}))

read_sensitivity_increments <- function(family, label) {
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
      analysis = label
    )
}

mundlak <- read_sensitivity_increments("mundlak", "Mundlak")
within <- read_sensitivity_increments("within", "Strict within-study")

within_support <- read.csv(gzfile(file.path(
  "results", "verification", "within",
  "exposure_support_by_group.csv.gz"
))) |>
  select(
    outcome,
    resistance_exercise_code,
    observed_minimum,
    observed_maximum
  )

plot_data <- bind_rows(primary, mundlak, within) |>
  left_join(
    within_support,
    by = c("outcome", "resistance_exercise_code")
  ) |>
  filter(
    protein_intake >= observed_minimum - 1e-10,
    protein_intake <= observed_maximum + 1e-10
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
    analysis = factor(
      analysis,
      levels = c("Primary", "Mundlak", "Strict within-study")
    )
  )

analysis_colours <- c(
  Primary = "#2B6CB0",
  Mundlak = "#D97706",
  `Strict within-study` = "#238B45"
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
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.35) +
  geom_vline(
    xintercept = 0.8,
    colour = "grey50",
    linewidth = 0.4,
    linetype = "dotted"
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.09,
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
    title = "Sensitivity analyses of the protein dose-response shape",
    subtitle = paste(
      "Posterior change from 0.8 g/kg/day; shaded regions are 95%",
      "credible intervals"
    ),
    x = "Protein intake (g/kg/day)",
    y = "Change in standardised outcome from 0.8 g/kg/day",
    colour = NULL,
    fill = NULL,
    caption = paste(
      "Displayed intake ranges are limited to observed treatment-intake",
      "support in the strict within-study data. The strict within-study",
      "analysis compares matched randomized arms and is closest to the",
      "causal dose-contrast estimand."
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
  "primary_mundlak_within_increment_overlay.png"
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
