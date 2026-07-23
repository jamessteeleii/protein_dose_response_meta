plot_linear_predictions <- function(
    fit,
    protein_reference,
    lower = 0.6,
    upper = 4.4,
    step = 0.02,
    draws_to_plot = 2000,
    actual_data
) {
  
  draws <- posterior::as_draws_df(
    fit$draws(c("b_intercept", "b_rt", "b_protein"))
  )
  
  draws <- draws |>
    dplyr::slice_sample(n = min(2000, nrow(draws)))
  
  prediction_grid <- tidyr::crossing(
    protein_intake_actual = seq(lower, upper, by = step),
    resistance_exercise_code = c(0, 1)
  ) |>
    dplyr::mutate(
      protein_intake = protein_intake_actual - protein_reference,
      resistance_exercise_yes_no = dplyr::if_else(
        resistance_exercise_code == 0,
        "yes",
        "no"
      )
    )
  
  prediction_draws <- tidyr::crossing(draws, prediction_grid) |>
    dplyr::mutate(
      estimate =
        b_intercept +
        b_rt * resistance_exercise_code +
        b_protein * protein_intake
    )
  
  prediction_summary <- prediction_draws |>
    dplyr::group_by(
      protein_intake_actual,
      resistance_exercise_yes_no
    ) |>
    dplyr::summarise(
      median = median(estimate),
      conf.low = quantile(estimate, 0.025),
      conf.high = quantile(estimate, 0.975),
      .groups = "drop"
    )
  
  ggplot2::ggplot(
    prediction_summary,
    ggplot2::aes(
      x = protein_intake_actual,
      y = median
    )
  ) +
    ggplot2::geom_point(
      data = actual_data |> filter(outcome_group == "whole body lean mass"),
      aes(x = protein_intake, y = yi, color = condition)
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.25
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~ resistance_exercise_yes_no) +
    ggplot2::labs(
      x = "Protein intake (g/kg/day)",
      y = "Expected lean-mass/muscle change (SMD)",
      title = "Linear measurement-error model",
      subtitle = "Population-level posterior predictions"
    ) +
    ggplot2::theme_classic()
}


plot_linear_spline_predictions <- function(
    fit,
    protein_reference,
    lower = 0.4,
    upper = 3.0,
    step = 0.02,
    draws_to_plot = 2000,
    actual_data
) {
  draws <- posterior::as_draws_df(
    fit$draws(c("b_intercept", "b_rt", "b_left", "b_change"))
  )
  
  draws <- draws |>
    dplyr::slice_sample(n = min(2000, nrow(draws)))
  
  prediction_grid <- tidyr::crossing(
    protein_intake_actual = seq(lower, upper, by = step),
    resistance_exercise_code = c(0, 1)
  ) |>
    dplyr::mutate(
      protein_intake = protein_intake_actual - protein_reference,
      resistance_exercise_yes_no = dplyr::if_else(
        resistance_exercise_code == 0,
        "yes",
        "no"
      )
    )
  
  prediction_draws <- tidyr::crossing(draws, prediction_grid) |>
    dplyr::mutate(
      spline_effect =
        b_left * pmin(protein_intake, 1.6 - protein_reference) +
        b_change * pmax(
          protein_intake - (1.6 - protein_reference),
          0
        ),
      
      estimate =
        b_intercept +
        b_rt * resistance_exercise_code +
        spline_effect
    )
  
  prediction_summary <- prediction_draws |>
    dplyr::group_by(
      protein_intake_actual,
      resistance_exercise_yes_no
    ) |>
    dplyr::summarise(
      conf.low = quantile(estimate, 0.025),
      conf.high = quantile(estimate, 0.975),
      estimate = median(estimate),
      .groups = "drop"
    )
  
  ggplot2::ggplot(
    prediction_summary,
    ggplot2::aes(
      x = protein_intake_actual,
      y = estimate
    )
  ) +
    ggplot2::geom_point(
      data = actual_data |> filter(outcome_group == "whole body lean mass"),
      aes(x = protein_intake, y = yi, color = condition)
    ) +
    ggplot2::geom_vline(
      xintercept = 1.6,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.25
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~ resistance_exercise_yes_no) +
    ggplot2::labs(
      x = "Protein intake (g/kg/day)",
      y = "Expected lean-mass/muscle change (SMD)",
      title = "Linear-spline measurement-error model",
      subtitle = "Dashed line: 1.6 g/kg/day threshold"
    ) +
    ggplot2::theme_classic()
}


library(dplyr)
library(ggplot2)
library(posterior)

rcs_orthogonal_basis_r <- function(x, spec) {
  knots <- spec$rcs_knots
  
  d <- vapply(
    knots[1:3],
    function(knot) {
      (
        pmax(x - knot, 0)^3 -
          pmax(x - knots[4], 0)^3
      ) / (knots[4] - knot)
    },
    numeric(length(x))
  )
  
  raw_basis <- cbind(
    d[, 1] - d[, 3],
    d[, 2] - d[, 3]
  )
  
  residual_basis <- raw_basis -
    cbind(1, x) %*% spec$rcs_projection
  
  residual_basis %*%
    spec$rcs_basis_transform
}

plot_rcs_orthogonal_predictions <- function(
    fit,
    spec,
    protein_reference,
    x_limits,
    resistance_exercise_code = 0,
    n_grid = 200,
    interval = 0.95,
    actual_data
) {
  protein_intake <- seq(
    x_limits[1],
    x_limits[2],
    length.out = n_grid
  )
  
  # Stan was fitted using protein intake centred at protein_reference.
  x_centred <- protein_intake - protein_reference
  
  nonlinear_basis <- rcs_orthogonal_basis_r(
    x_centred,
    spec
  )
  
  draws <- fit$draws(
    variables = c(
      "b_intercept",
      "b_rt",
      "b_linear",
      "theta_rcs"
    )
  ) |>
    as_draws_df()
  
  prediction_draws <-
    outer(
      draws$b_intercept +
        draws$b_rt * resistance_exercise_code,
      rep(1, n_grid)
    ) +
    outer(draws$b_linear, x_centred) +
    outer(draws$`theta_rcs[1]`, nonlinear_basis[, 1]) +
    outer(draws$`theta_rcs[2]`, nonlinear_basis[, 2])
  
  alpha <- (1 - interval) / 2
  
  prediction_summary <- tibble(
    protein_intake = protein_intake,
    estimate = apply(prediction_draws, 2, median),
    conf.low = apply(
      prediction_draws,
      2,
      quantile,
      probs = alpha
    ),
    conf.high = apply(
      prediction_draws,
      2,
      quantile,
      probs = 1 - alpha
    )
  )
  
  plot <- ggplot(
    prediction_summary,
    aes(x = protein_intake, y = estimate)
  ) +
    ggplot2::geom_point(
      data = actual_data |> filter(outcome_group == "whole body lean mass"),
      aes(x = protein_intake, y = yi, color = condition)
    ) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      fill = "#56B4E9",
      alpha = 0.25
    ) +
    geom_line(
      colour = "#0072B2",
      linewidth = 1
    ) +
    labs(
      x = "Protein intake (g/kg/day)",
      y = "Population-level predicted effect size",
      subtitle = paste0(
        round(interval * 100),
        "% posterior interval; random intercepts excluded"
      )
    ) +
    theme_minimal()
  
  list(
    summary = prediction_summary,
    draws = prediction_draws,
    plot = plot
  )
}
