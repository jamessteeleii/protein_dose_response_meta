functions {
  // Two nonlinear restricted-cubic-spline basis functions, made
  // orthogonal to the intercept, linear protein effect, and each other.
  vector rcs_orthogonal_basis(
    real x,
    vector rcs_knots,
    matrix rcs_projection,
    matrix rcs_basis_transform
  ) {
    vector[3] d;
    vector[2] raw_basis;
    vector[2] residual_basis;

    // Restricted/natural cubic construction.
    for (j in 1:3) {
      d[j] =
        (
          pow(fmax(x - rcs_knots[j], 0), 3) -
          pow(fmax(x - rcs_knots[4], 0), 3)
        ) /
        (rcs_knots[4] - rcs_knots[j]);
    }

    for (j in 1:2) {
      raw_basis[j] = d[j] - d[3];

      // Remove the parts represented by the intercept and linear x term.
      residual_basis[j] =
        raw_basis[j] -
        rcs_projection[1, j] -
        rcs_projection[2, j] * x;
    }

    // Rotate and scale the residual curvature bases.
    return rcs_basis_transform' * residual_basis;
  }
}

data {
  // Effect sizes and their known sampling standard errors.
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se_y;

  // Arm-level measured protein intake and its measurement error.
  int<lower=1> N_arm;
  array[N] int<lower=1, upper=N_arm> arm;
  vector[N_arm] x_obs;
  vector<lower=0>[N_arm] x_se;

  // Globally unique study indices.
  int<lower=1> N_study;
  array[N] int<lower=1, upper=N_study> study;

  // Globally unique outcome indices.
  int<lower=1> N_outcome;
  array[N] int<lower=1, upper=N_outcome> outcome;

  // Resistance-training indicator/code.
  vector[N] resistance_exercise_code;

  // Fixed restricted-cubic-spline specification created in R.
  vector[4] rcs_knots;
  matrix[2, 2] rcs_projection;
  matrix[2, 2] rcs_basis_transform;

  // Prior information from the existing prior target.
  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;
  real prior_rt_location;
  real<lower=0> prior_rt_scale;

  // Direct prior SD for each standardised nonlinear curvature component.
  real<lower=0> prior_rcs_scale;
}

parameters {
  // Non-centred measurement-error representation of true protein intake.
  vector[N_arm] z_x;
  real x_mean;
  real<lower=0> x_sd;

  // Population-level effects.
  real b_intercept;
  real b_rt;
  real b_linear;
  vector[2] theta_rcs;

  // Non-centred random intercepts.
  vector[N_study] z_study;
  vector[N_arm] z_arm;
  vector[N_outcome] z_outcome;

  real<lower=0> sd_study;
  real<lower=0> sd_arm;
  real<lower=0> sd_outcome;
}

transformed parameters {
  vector[N_arm] x_true;
  vector[N] mu;

  // Measurement-error model for protein intake.
  x_true = x_obs + x_se .* z_x;

  for (n in 1:N) {
    vector[2] nonlinear_basis;
    real protein_effect;

    nonlinear_basis = rcs_orthogonal_basis(
      x_true[arm[n]],
      rcs_knots,
      rcs_projection,
      rcs_basis_transform
    );

    protein_effect =
      b_linear * x_true[arm[n]] +
      dot_product(theta_rcs, nonlinear_basis);

    mu[n] =
      b_intercept +
      b_rt * resistance_exercise_code[n] +
      protein_effect +
      sd_study * z_study[study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
  }
}

model {
  // Latent protein-intake distribution.
  x_mean ~ normal(0, 2);
  x_sd ~ exponential(1);
  z_x ~ std_normal();
  x_true ~ normal(x_mean, x_sd);

  // Population-level priors.
  b_intercept ~ student_t(
    3,
    prior_intercept_location,
    prior_intercept_scale
  );

  b_rt ~ student_t(
    3,
    prior_rt_location,
    prior_rt_scale
  );

  // Corrected common prior scale for a linear protein effect.
  b_linear ~ student_t(3, 0, 1);

  // Direct regularisation of the two orthogonal curvature components.
  theta_rcs ~ normal(0, prior_rcs_scale);

  // Random-intercept priors.
  z_study ~ std_normal();
  z_arm ~ std_normal();
  z_outcome ~ std_normal();

  sd_study ~ student_t(3, 0, 2.5);
  sd_arm ~ student_t(3, 0, 2.5);
  sd_outcome ~ student_t(3, 0, 2.5);

  // Meta-analytic sampling model.
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
