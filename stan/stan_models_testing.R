targets::tar_load(c(protein_data_prepared, protein_data_complete, stan_data_whole_body, steele_priors, stan_priors))


stan_data <- c(
  stan_data_whole_body,
  list(
    prior_intercept_location = stan_priors$prior_intercept_location,
    prior_intercept_scale = stan_priors$prior_intercept_scale,
    prior_rt_location = stan_priors$prior_rt_location,
    prior_rt_scale = stan_priors$prior_rt_scale
  )
)

# Intercept only model first ----

intercept_model <- cmdstanr::cmdstan_model(
  "stan/intercept_only.stan"
)

intercept_test_fit <- intercept_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 6000,
  adapt_delta = 0.95,
  seed = 1988
)

intercept_test_fit$summary()
intercept_test_fit$cmdstan_diagnose()

# Linear model ----

linear_model <- cmdstanr::cmdstan_model(
  "stan/linear.stan"
)

linear_fit <- linear_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 6000,
  adapt_delta = 0.95,
  seed = 1988
)

linear_fit$summary(
  variables = c("b_intercept", "b_rt", "b_protein", "x_mean", "x_sd")
)

linear_fit$cmdstan_diagnose()

plot_linear_predictions(
  fit = linear_fit, protein_reference = median(protein_data_prepared$protein_intake, na.rm = TRUE),
  actual_data = protein_data_complete |> filter(outcome_group == "whole body lean mass")
)


# Linear spline model ----

stan_data_linear_spline <- c(
  stan_data,
  list(
    protein_knot = 1.6 - median(protein_data_prepared$protein_intake, na.rm = TRUE)
  )
)

linear_spline_model <- cmdstanr::cmdstan_model(
  "stan/linear_spline.stan"
)

linear_spline_fit <- linear_spline_model$sample(
  data = stan_data_linear_spline,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 6000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 1988
)

linear_spline_fit$cmdstan_diagnose()

draws <- posterior::as_draws_df(
  linear_spline_fit$draws(c("b_left", "b_change"))
) |>
  dplyr::mutate(
    slope_below_1_6 = b_left,
    slope_above_1_6 = b_left + b_change
  )

posterior::summarise_draws(
  posterior::as_draws_matrix(
    posterior::as_draws_df(draws)[
      c("slope_below_1_6", "slope_above_1_6")
    ]
  )
)

linear_spline_fit$summary(
  variables = c("b_intercept", "b_rt", "b_protein", "x_mean", "x_sd")
)

plot_linear_spline_predictions(
  fit = linear_spline_fit, protein_reference = median(protein_data_prepared$protein_intake, na.rm = TRUE),
  actual_data = protein_data_complete |> filter(outcome_group == "whole body lean mass")
)

# # Thin plate spline ----
# 
# make_tps_basis_constants <- function(x_obs, tps_knots) {
#   linear_design <- cbind(
#     intercept = 1,
#     protein_intake = x_obs
#   )
#   
#   radial_basis <- outer(
#     x_obs,
#     tps_knots,
#     FUN = \(x, knot) abs(x - knot)^3
#   )
#   
#   # Remove intercept and linear components.
#   projection <- qr.solve(linear_design, radial_basis)
#   
#   residual_basis <- radial_basis - linear_design %*% projection
#   
#   # Orthonormalise the two remaining nonlinear basis columns.
#   qr_basis <- qr(residual_basis)
#   
#   if (qr_basis$rank < ncol(residual_basis)) {
#     stop("Thin-plate nonlinear basis is not full rank.")
#   }
#   
#   rotation <- solve(qr.R(qr_basis)) * sqrt(nrow(residual_basis))
#   
#   # Verification: should be approximately the identity matrix.
#   orthogonal_basis <- residual_basis %*% rotation
#   print(crossprod(orthogonal_basis) / nrow(orthogonal_basis))
#   
#   list(
#     tps_projection = projection,
#     tps_basis_transform = rotation
#   )
# }
# 
# 
# tps_knots <- unname(
#   quantile(
#     stan_data$x_obs,
#     probs = c(1 / 3, 2 / 3)
#   )
# )
# 
# tps_constants <- make_tps_basis_constants(
#   x_obs = stan_data$x_obs,
#   tps_knots = tps_knots
# )
# 
# stan_data_thin_plate <- c(
#   stan_data,
#   list(
#     N_tps_knots = 2L,
#     tps_knots = tps_knots
#   ),
#   tps_constants
# )
# 
# thin_plate_model <- cmdstanr::cmdstan_model(
#   "stan/thin_plate_spline.stan"
# )
# 
# thin_plate_fit <- thin_plate_model$sample(
#   data = stan_data_thin_plate,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 2000,
#   iter_sampling = 2000,
#   adapt_delta = 0.999,
#   max_treedepth = 15,
#   seed = 1988
# )
# 
# thin_plate_fit$cmdstan_diagnose()

# Restricted cubic spline ----

make_rcs_orthogonal_spec <- function(
    x_obs,
    probs = c(0.05, 0.35, 0.65, 0.95)
) {
  stopifnot(is.numeric(x_obs), all(is.finite(x_obs)))
  
  knots <- as.numeric(
    stats::quantile(x_obs, probs = probs, names = FALSE)
  )
  
  if (any(diff(knots) <= 0)) {
    stop("Restricted cubic-spline knots are not distinct.")
  }
  
  rcs_raw_basis <- function(x) {
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
    
    sweep(
      d[, 1:2, drop = FALSE],
      1,
      d[, 3],
      "-"
    )
  }
  
  H <- rcs_raw_basis(x_obs)
  
  # Remove the intercept and linear protein components from both
  # nonlinear spline columns.
  linear_design <- cbind(1, x_obs)
  rcs_projection <- qr.solve(linear_design, H)
  
  H_residual <- H - linear_design %*% rcs_projection
  
  # Rotate and scale the two residual curvature columns so that they
  # are orthogonal and approximately have SD = 1.
  qr_H <- qr(H_residual)
  
  if (qr_H$rank < 2) {
    stop("The two nonlinear cubic-spline columns are not independent.")
  }
  
  rcs_basis_transform <-
    solve(qr.R(qr_H)) * sqrt(nrow(H_residual))
  
  H_orthogonal <- H_residual %*% rcs_basis_transform
  
  list(
    rcs_knots = knots,
    rcs_projection = rcs_projection,
    rcs_basis_transform = rcs_basis_transform,
    max_linear_correlation = max(
      abs(cor(H_orthogonal, x_obs))
    ),
    basis_crossproduct = crossprod(H_orthogonal) / nrow(H_orthogonal)
  )
}

rcs_orthogonal_spec <- make_rcs_orthogonal_spec(stan_data$x_obs)

rcs_orthogonal_spec$rcs_knots
rcs_orthogonal_spec$max_linear_correlation
rcs_orthogonal_spec$basis_crossproduct

stan_data_rcs_orthogonal <- c(
  stan_data,
  list(
    rcs_knots = rcs_orthogonal_spec$rcs_knots,
    rcs_projection = rcs_orthogonal_spec$rcs_projection,
    rcs_basis_transform =
      rcs_orthogonal_spec$rcs_basis_transform,
    
    # Start with a generous direct curvature prior.
    prior_rcs_scale = 1
  )
)

rcs_model <- cmdstanr::cmdstan_model(
  "stan/restricted_cubic_spline.stan"
)

rcs_fit <- rcs_model$sample(
  data = stan_data_rcs_orthogonal,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  max_treedepth = 15
)

rcs_fit$cmdstan_diagnose()

(plot_rcs_orthogonal_predictions(
  fit = rcs_fit,
  spec = rcs_orthogonal_spec,
  protein_reference = median(protein_data_prepared$protein_intake, na.rm = TRUE),
  x_limits = range(stan_data$x_obs + protein_reference),
  resistance_exercise_code = 1,
  actual_data = protein_data_complete |> filter(outcome_group == "whole body lean mass" &
                                                  resistance_exercise_yes_no == "yes")
)$plot) +

(plot_rcs_orthogonal_predictions(
  fit = rcs_fit,
  spec = rcs_orthogonal_spec,
  protein_reference = median(protein_data_prepared$protein_intake, na.rm = TRUE),
  x_limits = range(stan_data$x_obs + protein_reference),
  resistance_exercise_code = 0,
  actual_data = protein_data_complete |> filter(outcome_group == "whole body lean mass" &
                                                  resistance_exercise_yes_no == "no")
)$plot)

