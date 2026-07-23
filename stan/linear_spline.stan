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

  // 1.6 g/kg/day expressed on the centred protein scale.
  real protein_knot;

  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;
  real prior_rt_location;
  real<lower=0> prior_rt_scale;
}

parameters {
  // vector[N_arm] x_true;
  vector[N_arm] z_x;
  real x_mean;
  real<lower=0> x_sd;

  real b_intercept;
  real b_rt;

  // Dose-response coefficients.
  real b_left;
  real b_change;

  vector[N_study] z_study;
  vector[N_arm] z_arm;
  vector[N_outcome] z_outcome;

  real<lower=0> sd_study;
  real<lower=0> sd_arm;
  real<lower=0> sd_outcome;
}

transformed parameters {
  // vector[N] mu;
  // 
  // for (n in 1:N) {
  //   mu[n] =
  //     b_intercept +
  //     b_rt * resistance_exercise_code[n] +
  // 
  //     // Continuous linear spline:
  //     // below knot: slope = b_left
  //     // above knot: slope = b_left + b_change
  //     b_left * x_true[arm[n]] +
  //     b_change * fmax(x_true[arm[n]] - protein_knot, 0) +
  // 
  //     sd_study * z_study[study[n]] +
  //     sd_arm * z_arm[arm[n]] +
  //     sd_outcome * z_outcome[outcome[n]];
  // }
  
  vector[N_arm] x_true;
  vector[N] mu;

  x_true = x_obs + x_se .* z_x;

  for (n in 1:N) {
    mu[n] =
      b_intercept +
      b_rt * resistance_exercise_code[n] +
      b_left * x_true[arm[n]] +
      b_change * fmax(x_true[arm[n]] - protein_knot, 0) +
      sd_study * z_study[study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
  }
}

model {
  // Latent protein-intake exposure model.
  x_mean ~ normal(0, 2);
  x_sd ~ exponential(1);
  // x_true ~ normal(x_mean, x_sd);
  z_x ~ std_normal();

  // Known intake measurement error.
  // x_obs ~ normal(x_true, x_se);
  
  // Population distribution of latent true protein intake.
  x_true ~ normal(x_mean, x_sd);

  // Outcome-model priors.
  b_intercept ~ student_t(
    3,
    prior_intercept_location,
    prior_intercept_scale
  );

  b_rt ~ student_t(3, prior_rt_location, prior_rt_scale);
  b_left ~ student_t(3, 0, 1);
  b_change ~ student_t(3, 0, 1);

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
