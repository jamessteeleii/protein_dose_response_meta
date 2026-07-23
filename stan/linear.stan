data {
  // One row per effect size.
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se_y;

  // One observed protein intake per nested study–arm.
  int<lower=1> N_arm;
  array[N] int<lower=1, upper=N_arm> arm;
  vector[N_arm] x_obs;
  vector<lower=0>[N_arm] x_se;

  // Nested random-intercept identifiers.
  int<lower=1> N_study;
  array[N] int<lower=1, upper=N_study> study;

  int<lower=1> N_outcome;
  array[N] int<lower=1, upper=N_outcome> outcome;

  vector[N] resistance_exercise_code;

  // Steele-derived fixed-effect priors.
  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;

  real prior_rt_location;
  real<lower=0> prior_rt_scale;
}

parameters {
  // Latent true protein intake for each study–arm.
  vector[N_arm] x_true;

  // Population distribution of centred true protein intake.
  real x_mean;
  real<lower=0> x_sd;

  // Population-level outcome effects.
  real b_intercept;
  real b_rt;
  real b_protein;

  // Non-centred nested random intercepts.
  vector[N_study] z_study;
  vector[N_arm] z_arm;
  vector[N_outcome] z_outcome;

  real<lower=0> sd_study;
  real<lower=0> sd_arm;
  real<lower=0> sd_outcome;
}

transformed parameters {
  vector[N] mu;

  for (n in 1:N) {
    mu[n] =
      b_intercept +
      b_rt * resistance_exercise_code[n] +
      b_protein * x_true[arm[n]] +
      sd_study * z_study[study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
  }
}

model {
  // Exposure model: protein_intake is already centred.
  x_mean ~ normal(0, 2);
  x_sd ~ exponential(1);
  x_true ~ normal(x_mean, x_sd);

  // Measurement-error model.
  x_obs ~ normal(x_true, x_se);

  // Outcome-model priors.
  b_intercept ~ student_t(
    3,
    prior_intercept_location,
    prior_intercept_scale
  );

  b_rt ~ student_t(3, prior_rt_location, prior_rt_scale);
  b_protein ~ student_t(3, 0, 1);

  // Hierarchical heterogeneity priors.
  z_study ~ std_normal();
  z_arm ~ std_normal();
  z_outcome ~ std_normal();

  sd_study ~ student_t(3, 0, 2.5);
  sd_arm ~ student_t(3, 0, 2.5);
  sd_outcome ~ student_t(3, 0, 2.5);

  // Meta-analytic outcome likelihood.
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
