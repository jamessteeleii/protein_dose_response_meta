functions {
  vector rcs_curvature_basis(
    real x,
    vector knots,
    matrix projection,
    matrix basis_transform
  ) {
    vector[3] d;
    vector[2] raw_basis;
    vector[2] residual_basis;
    for (j in 1:3) {
      d[j] = (
        pow(fmax(x - knots[j], 0), 3) -
        pow(fmax(x - knots[4], 0), 3)
      ) / (knots[4] - knots[j]);
    }
    for (j in 1:2) {
      raw_basis[j] = d[j] - d[3];
      residual_basis[j] = raw_basis[j] -
        projection[1, j] - projection[2, j] * x;
    }
    return basis_transform' * residual_basis;
  }

  vector protein_basis(
    real x,
    int basis_type,
    int K_basis,
    real protein_scale,
    real protein_knot,
    vector rcs_knots,
    matrix rcs_projection,
    matrix rcs_basis_transform
  ) {
    vector[K_basis] basis;
    real z = x / protein_scale;
    if (basis_type == 0) {
      basis[1] = 0;
    } else if (basis_type == 1) {
      basis[1] = z;
    } else if (basis_type == 2) {
      real scaled_knot = protein_knot / protein_scale;
      basis[1] = fmin(z, scaled_knot);
      basis[2] = fmax(z - scaled_knot, 0);
    } else {
      vector[2] curvature = rcs_curvature_basis(
        z, rcs_knots, rcs_projection, rcs_basis_transform
      );
      basis[1] = z;
      basis[2:3] = curvature;
    }
    return basis;
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
  vector[N] age_z;
  vector[N] duration_z;
  vector[N] body_mass_z;

  int<lower=0, upper=3> basis_type;
  int<lower=1, upper=3> K_basis;
  int<lower=0, upper=1> include_protein;
  real<lower=0> protein_scale;
  real protein_knot;
  vector[4] rcs_knots;
  matrix[2, 2] rcs_projection;
  matrix[2, 2] rcs_basis_transform;

  int<lower=0, upper=1> include_moderators;
  int<lower=0, upper=1> include_rt_interaction;
  int<lower=0, upper=1> include_random_slope;
  int<lower=1, upper=4> K_study;
  int<lower=0, upper=1> correlated_study_effects;

  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;
  real prior_rt_location;
  real<lower=0> prior_rt_scale;
  vector<lower=0>[K_basis] prior_protein_scale;
  vector<lower=0>[K_study] prior_study_sd_scale;
  real<lower=0> prior_arm_sd_scale;
  real<lower=0> prior_outcome_sd_scale;

  int<lower=1> N_grid;
  vector[N_grid] x_grid;
  vector[N_grid] rt_grid;
  vector[N_grid] age_grid_z;
  vector[N_grid] duration_grid_z;
  vector[N_grid] body_mass_grid_z;
}

transformed data {
  if (basis_type == 0 && K_basis != 1) {
    reject("The null basis requires K_basis = 1.");
  }
  if (basis_type == 1 && K_basis != 1) {
    reject("The linear basis requires K_basis = 1.");
  }
  if (basis_type == 2 && K_basis != 2) {
    reject("The piecewise basis requires K_basis = 2.");
  }
  if (basis_type == 3 && K_basis != 3) {
    reject("The cubic-spline basis requires K_basis = 3.");
  }
  if (K_study != 1 + include_random_slope * K_basis) {
    reject("K_study is inconsistent with the random-slope flag.");
  }
}

parameters {
  vector[N_arm] z_x;
  real x_mean;
  real<lower=0> x_sd;

  real b_intercept;
  real b_rt;
  vector[K_basis] b_protein;
  vector[K_basis] b_protein_rt;
  vector[3] b_moderator;

  matrix[K_study, N_study] z_study;
  vector<lower=0>[K_study] sd_study;
  cholesky_factor_corr[K_study] L_study;
  vector[N_arm] z_arm;
  vector[N_outcome] z_outcome;
  real<lower=0> sd_arm;
  real<lower=0> sd_outcome;
}

transformed parameters {
  vector[N_arm] x_true = x_obs + x_se .* z_x;
  matrix[K_study, N_study] study_effect;
  vector[N] mu;

  if (correlated_study_effects == 1) {
    study_effect = diag_pre_multiply(sd_study, L_study) * z_study;
  } else {
    study_effect = diag_pre_multiply(
      sd_study, diag_matrix(rep_vector(1.0, K_study))
    ) * z_study;
  }

  for (n in 1:N) {
    vector[K_basis] basis = protein_basis(
      x_true[arm[n]], basis_type, K_basis, protein_scale,
      protein_knot, rcs_knots, rcs_projection,
      rcs_basis_transform
    );
    real eta = b_intercept +
      b_rt * resistance_exercise_code[n] +
      study_effect[1, study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
    if (include_protein == 1) {
      eta += dot_product(b_protein, basis);
    }
    if (include_rt_interaction == 1) {
      eta += resistance_exercise_code[n] *
        dot_product(b_protein_rt, basis);
    }
    if (include_moderators == 1) {
      eta += b_moderator[1] * duration_z[n] +
        b_moderator[2] * age_z[n] +
        b_moderator[3] * body_mass_z[n];
    }
    if (include_random_slope == 1) {
      eta += dot_product(
        study_effect[2:K_study, study[n]], basis
      );
    }
    mu[n] = eta;
  }
}

model {
  z_x ~ std_normal();
  x_mean ~ normal(0, 2);
  x_sd ~ exponential(1);
  x_true ~ normal(x_mean, x_sd);

  b_intercept ~ student_t(
    3, prior_intercept_location, prior_intercept_scale
  );
  b_rt ~ student_t(3, prior_rt_location, prior_rt_scale);
  for (k in 1:K_basis) {
    b_protein[k] ~ student_t(3, 0, prior_protein_scale[k]);
    b_protein_rt[k] ~ student_t(3, 0, prior_protein_scale[k]);
  }
  b_moderator ~ student_t(3, 0, 0.5);

  to_vector(z_study) ~ std_normal();
  L_study ~ lkj_corr_cholesky(1);
  for (k in 1:K_study) {
    sd_study[k] ~ student_t(3, 0, prior_study_sd_scale[k]);
  }
  z_arm ~ std_normal();
  z_outcome ~ std_normal();
  sd_arm ~ student_t(3, 0, prior_arm_sd_scale);
  sd_outcome ~ student_t(3, 0, prior_outcome_sd_scale);

  y ~ normal(mu, se_y);
}

generated quantities {
  vector[N] log_lik_y;
  vector[N] y_rep;
  vector[N_grid] mu_grid;
  corr_matrix[K_study] study_correlation;

  if (correlated_study_effects == 1) {
    study_correlation = multiply_lower_tri_self_transpose(L_study);
  } else {
    study_correlation = diag_matrix(rep_vector(1.0, K_study));
  }
  for (n in 1:N) {
    log_lik_y[n] = normal_lpdf(y[n] | mu[n], se_y[n]);
    y_rep[n] = normal_rng(mu[n], se_y[n]);
  }
  for (g in 1:N_grid) {
    vector[K_basis] basis = protein_basis(
      x_grid[g], basis_type, K_basis, protein_scale,
      protein_knot, rcs_knots, rcs_projection,
      rcs_basis_transform
    );
    real eta = b_intercept + b_rt * rt_grid[g];
    if (include_protein == 1) {
      eta += dot_product(b_protein, basis);
    }
    if (include_rt_interaction == 1) {
      eta += rt_grid[g] * dot_product(b_protein_rt, basis);
    }
    if (include_moderators == 1) {
      eta += b_moderator[1] * duration_grid_z[g] +
        b_moderator[2] * age_grid_z[g] +
        b_moderator[3] * body_mass_grid_z[g];
    }
    mu_grid[g] = eta;
  }
}
