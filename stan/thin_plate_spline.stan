// functions {
//   // One-dimensional, second-order thin-plate radial basis.
//   real tps_radial(real x, real knot) {
//     return pow(abs(x - knot), 3);
//   }
// }

// functions {
//   real tps_basis(
//     real x,
//     int j,
//     vector knots,
//     matrix projection,
//     vector basis_scale
//   ) {
//     real raw_basis = pow(abs(x - knots[j]), 3);
// 
//     return (
//       raw_basis -
//       projection[1, j] -
//       projection[2, j] * x
//     ) / basis_scale[j];
//   }
// }

functions {
  vector tps_basis_vector(
    real x,
    vector knots,
    matrix projection,
    matrix basis_transform
  ) {
    int K = num_elements(knots);
    vector[K] residual_basis;

    for (j in 1:K) {
      residual_basis[j] =
        pow(abs(x - knots[j]), 3) -
        projection[1, j] -
        projection[2, j] * x;
    }

    // R uses row-vector %*% transform.
    // Stan uses transform' * column-vector.
    return basis_transform' * residual_basis;
  }
}

data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se_y;

  int<lower=1> N_arm;
  array[N] int<lower=1, upper=N_arm> arm;
  vector[N_arm] x_obs;
  vector<lower=0>[N_arm] x_se;

  int<lower=1> N_study;
  array[N] int<lower=1, upper=N_study> study;

  int<lower=1> N_outcome;
  array[N] int<lower=1, upper=N_outcome> outcome;

  vector[N] resistance_exercise_code;

  // k = 4 total basis dimension:
  // intercept + linear term + two penalised radial basis functions.
  int<lower=1> N_tps_knots;
  vector[N_tps_knots] tps_knots;
  
  matrix[2, N_tps_knots] tps_projection;
  // vector<lower=0>[N_tps_knots] tps_basis_scale;
  matrix[N_tps_knots, N_tps_knots] tps_basis_transform;

  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;
  real prior_rt_location;
  real<lower=0> prior_rt_scale;
}

parameters {
  // Non-centred measurement-error representation.
  vector[N_arm] z_x;
  real x_mean;
  real<lower=0> x_sd;

  real b_intercept;
  real b_rt;

  // Unpenalised linear component of the thin-plate spline.
  real b_linear;

  // // Penalised nonlinear part of the thin-plate spline.
  // vector[N_tps_knots] z_tps;
  // real<lower=0> sd_tps;
  vector[N_tps_knots] theta_tps;
  real<lower=0> sd_tps;

  vector[N_study] z_study;
  vector[N_arm] z_arm;
  vector[N_outcome] z_outcome;

  real<lower=0> sd_study;
  real<lower=0> sd_arm;
  real<lower=0> sd_outcome;
}

transformed parameters {
  vector[N_arm] x_true;
  // vector[N_tps_knots] theta_tps;
  vector[N] mu;

  x_true = x_obs + x_se .* z_x;
  // theta_tps = sd_tps * z_tps;

    for (n in 1:N) {
    real smooth_effect;
    vector[N_tps_knots] basis;

    // Unpenalised linear component.
    smooth_effect = b_linear * x_true[arm[n]];

    // Both orthogonalised nonlinear thin-plate basis functions.
    basis = tps_basis_vector(
      x_true[arm[n]],
      tps_knots,
      tps_projection,
      tps_basis_transform
    );

    smooth_effect += dot_product(theta_tps, basis);

    mu[n] =
      b_intercept +
      b_rt * resistance_exercise_code[n] +
      smooth_effect +
      sd_study * z_study[study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
  }
}

model {
  // Latent protein-intake exposure model.
  x_mean ~ normal(0, 2);
  x_sd ~ exponential(1);
  z_x ~ std_normal();
  x_true ~ normal(x_mean, x_sd);

  // Outcome-model priors.
  b_intercept ~ student_t(
    3,
    prior_intercept_location,
    prior_intercept_scale
  );

  b_rt ~ student_t(3, prior_rt_location, prior_rt_scale);

  // Same prior scale as the linear protein slope.
  b_linear ~ student_t(3, 0, 1);

  // Shrinkage of nonlinearity toward the linear model.
  // z_tps ~ std_normal();
  // sd_tps ~ exponential(1);
  theta_tps ~ normal(0, sd_tps);
  sd_tps ~ exponential(2);

  z_study ~ std_normal();
  z_arm ~ std_normal();
  z_outcome ~ std_normal();

  sd_study ~ student_t(3, 0, 2.5);
  sd_arm ~ student_t(3, 0, 2.5);
  sd_outcome ~ student_t(3, 0, 2.5);

  y ~ normal(mu, se_y);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N_study] log_lik_study = rep_vector(0.0, N_study);

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], se_y[n]);
    log_lik_study[study[n]] += log_lik[n];
    y_rep[n] = normal_rng(mu[n], se_y[n]);
  }
}
