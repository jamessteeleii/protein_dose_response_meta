#!/usr/bin/env Rscript

# Temporary visual comparison of the available sensitivity pilots with the
# preregistered primary model average. This script reads completed targets only;
# it does not make targets, compile Stan, or sample models.

suppressPackageStartupMessages({
  library(targets)
  library(qs)
  library(tidyverse)
})

reference_intake <- 0.8
plot_min <- 0.8
plot_max <- 2.4

outcome_labels <- c(
  whole_body = "Whole-body lean mass",
  hypertrophy = "Hypertrophy",
  strength = "Strength",
  "function" = "Physical function"
)

primary_targets <- c(
  whole_body = "primary_model_average_whole_body",
  hypertrophy = "primary_model_average_hypertrophy",
  strength = "primary_model_average_strength",
  "function" = "primary_model_average_function"
)

summarise_increment_average <- function(model_average, outcome) {
  grid <- model_average$grid
  draws <- as.matrix(model_average$draws)
  output <- vector("list", 2L)

  for (index in seq_along(c(0, 1))) {
    rt <- c(0, 1)[[index]]
    columns <- which(grid$resistance_exercise_code == rt)
    reference_column <- columns[[which.min(
      abs(grid$protein_intake[columns] - reference_intake)
    )]]
    increment_draws <- sweep(
      draws[, columns, drop = FALSE],
      1L,
      draws[, reference_column],
      FUN = "-"
    )
    intervals <- matrixStats::colQuantiles(
      increment_draws,
      probs = c(0.025, 0.975),
      drop = FALSE
    )
    output[[index]] <- tibble(
      outcome = outcome,
      resistance_exercise_code = rt,
      protein_intake = grid$protein_intake[columns],
      estimate = colMeans(increment_draws),
      conf_low = intervals[, 1L],
      conf_high = intervals[, 2L],
      family = "Primary preregistered"
    )
  }

  bind_rows(output)
}

summary_to_increment_line <- function(summary, outcome, family) {
  summary |>
    transmute(
      outcome = outcome,
      resistance_exercise_code,
      protein_intake,
      estimate
    ) |>
    group_by(outcome, resistance_exercise_code) |>
    arrange(protein_intake, .by_group = TRUE) |>
    mutate(
      reference_estimate = estimate[[which.min(
        abs(protein_intake - reference_intake)
      )]],
      estimate = estimate - reference_estimate,
      conf_low = NA_real_,
      conf_high = NA_real_,
      family = family
    ) |>
    ungroup() |>
    select(-reference_estimate)
}

report_increment_table <- function(report, outcome, family) {
  report$dose_response_increments |>
    transmute(
      outcome = outcome,
      resistance_exercise_code,
      protein_intake,
      estimate = increment,
      conf_low = increment_conf_low,
      conf_high = increment_conf_high,
      family = family
    )
}

pilot_reporting_increment_table <- function(report, outcome, family) {
  report$increments_full_grid |>
    transmute(
      outcome = outcome,
      resistance_exercise_code,
      protein_intake,
      estimate = increment,
      conf_low = increment_conf_low,
      conf_high = increment_conf_high,
      family = family
    )
}

rescor_increment_table <- function(fit_record) {
  protein_grid <- seq(plot_min, plot_max, by = 0.01)
  coefficient <- posterior::as_draws_df(
    fit_record$fit$draws("b_protein")
  )[["b_protein"]]
  increment_draws <- outer(
    coefficient,
    protein_grid - reference_intake
  )
  intervals <- matrixStats::colQuantiles(
    increment_draws,
    probs = c(0.025, 0.975),
    drop = FALSE
  )
  curve <- tibble(
    outcome = "whole_body",
    protein_intake = protein_grid,
    estimate = colMeans(increment_draws),
    conf_low = intervals[, 1L],
    conf_high = intervals[, 2L],
    family = "Residual correlation (model 6)"
  )
  tidyr::crossing(
    resistance_exercise_code = c(0, 1),
    curve
  )
}

primary <- imap_dfr(primary_targets, function(target, outcome) {
  summarise_increment_average(tar_read_raw(target), outcome)
})

mundlak_paths <- unname(unlist(
  tar_read_raw("sensitivity_mundlak_full_reporting_file"),
  use.names = FALSE
))
mundlak <- map_dfr(mundlak_paths, function(path) {
  report <- qs::qread(path)
  report_increment_table(report, report$outcome, "Mundlak (full)")
})

within <- pilot_reporting_increment_table(
  tar_read_raw("sensitivity_within_reporting_pilot"),
  "whole_body",
  "Strict within-study"
)

achieved <- summary_to_increment_line(
  tar_read_raw("sensitivity_achieved_final_summary"),
  "whole_body",
  "Achieved intake (common sample)"
)

assignment <- summary_to_increment_line(
  tar_read_raw("sensitivity_assignment_final_summary"),
  "whole_body",
  "Assignment-derived intake"
)

eiv <- summary_to_increment_line(
  tar_read_raw("sensitivity_eiv_pilot_summary"),
  "whole_body",
  "Errors in variables"
)

rescor <- rescor_increment_table(
  tar_read_raw("sensitivity_rescor_pilot_fit")
)

comparison <- bind_rows(
  primary,
  mundlak,
  within,
  achieved,
  assignment,
  eiv,
  rescor
) |>
  filter(
    protein_intake >= plot_min - 1e-10,
    protein_intake <= plot_max + 1e-10
  ) |>
  mutate(
    outcome = factor(
      outcome,
      levels = names(outcome_labels),
      labels = unname(outcome_labels)
    ),
    resistance_exercise = factor(
      resistance_exercise_code,
      levels = c(0, 1),
      labels = c("Resistance exercise", "No resistance exercise")
    ),
    family = factor(
      family,
      levels = c(
        "Primary preregistered",
        "Mundlak (full)",
        "Strict within-study",
        "Achieved intake (common sample)",
        "Assignment-derived intake",
        "Errors in variables",
        "Residual correlation (model 6)"
      )
    )
  )

palette <- c(
  "Primary preregistered" = "#000000",
  "Mundlak (full)" = "#0072B2",
  "Strict within-study" = "#D55E00",
  "Achieved intake (common sample)" = "#009E73",
  "Assignment-derived intake" = "#CC79A7",
  "Errors in variables" = "#E69F00",
  "Residual correlation (model 6)" = "#56B4E9"
)

line_types <- c(
  "Primary preregistered" = "solid",
  "Mundlak (full)" = "solid",
  "Strict within-study" = "longdash",
  "Achieved intake (common sample)" = "dotdash",
  "Assignment-derived intake" = "twodash",
  "Errors in variables" = "dashed",
  "Residual correlation (model 6)" = "dotted"
)

primary_band <- comparison |>
  filter(family == "Primary preregistered")

plot <- ggplot(
  comparison,
  aes(
    x = protein_intake,
    y = estimate,
    colour = family,
    linetype = family,
    group = family
  )
) +
  geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.35) +
  geom_ribbon(
    data = primary_band,
    aes(x = protein_intake, ymin = conf_low, ymax = conf_high),
    inherit.aes = FALSE,
    fill = "grey55",
    alpha = 0.18
  ) +
  geom_line(linewidth = 0.78, na.rm = TRUE) +
  facet_grid(
    rows = vars(outcome),
    cols = vars(resistance_exercise),
    scales = "free_y"
  ) +
  coord_cartesian(xlim = c(plot_min, plot_max)) +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_linetype_manual(values = line_types, drop = FALSE) +
  scale_x_continuous(breaks = seq(0.8, 2.4, by = 0.4)) +
  labs(
    title = "Sensitivity analyses compared with the preregistered dose-response",
    subtitle = paste0(
      "Posterior mean change from ", reference_intake,
      " g/kg/day; grey ribbon is the primary 95% credible interval"
    ),
    x = "Protein intake (g/kg/day)",
    y = "Change from response at 0.8 g/kg/day (standardized units)",
    colour = NULL,
    linetype = NULL,
    caption = paste(
      "Current availability: all outcomes for the primary and full Mundlak analyses;",
      "whole-body lean mass only for the remaining pilots.",
      "Curves beyond observed group-specific intake support should be interpreted as extrapolation."
    )
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    plot.title.position = "plot",
    plot.caption = element_text(hjust = 0, colour = "grey35")
  )

output_directory <- file.path("results", "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

png_path <- file.path(
  output_directory,
  "pilot_sensitivity_dose_response_overlay.png"
)
pdf_path <- file.path(
  output_directory,
  "pilot_sensitivity_dose_response_overlay.pdf"
)
data_path <- file.path(
  output_directory,
  "pilot_sensitivity_dose_response_overlay_data.csv.gz"
)

ggsave(
  png_path, plot,
  width = 13, height = 12, units = "in", dpi = 220, bg = "white"
)
ggsave(pdf_path, plot, width = 13, height = 12, units = "in", bg = "white")
readr::write_csv(comparison, data_path)

message(normalizePath(png_path, winslash = "/", mustWork = TRUE))
message(normalizePath(pdf_path, winslash = "/", mustWork = TRUE))
message(normalizePath(data_path, winslash = "/", mustWork = TRUE))
