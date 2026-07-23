data {
  // One row per effect size.
  int<lower=1> N;
  vector[N] y;                       // yi
  vector<lower=0>[N] se_y;          // sqrt(vi)

  // Globally unique nested identifiers.
  int<lower=1> N_study;
  array[N] int<lower=1, upper=N_study> study;

  int<lower=1> N_arm;
  array[N] int<lower=1, upper=N_arm> arm;

  int<lower=1> N_outcome;
  array[N] int<lower=1, upper=N_outcome> outcome;

  // 0 = resistance exercise; 1 = no resistance exercise.
  vector[N] resistance_exercise_code;

  // Numerical priors from the existing Steele-derived prior analysis.
  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;

  real prior_rt_location;
  real<lower=0> prior_rt_scale;
}

parameters {
  // Population-level fixed effects.
  real b_intercept;
  real b_rt;

  // Standard-normal latent variables: non-centred random effects.
  vector[N_study] z_study;
  vector[N_arm] z_arm;
  vector[N_outcome] z_outcome;

  // Estimated between-group heterogeneity.
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
      sd_study * z_study[study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
  }
}

model {
  // Same substantive priors as the brms model.
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

  // Priors for heterogeneity at each nested level.
  z_study ~ std_normal();
  z_arm ~ std_normal();
  z_outcome ~ std_normal();

  sd_study ~ student_t(3, 0, 2.5);
  sd_arm ~ student_t(3, 0, 2.5);
  sd_outcome ~ student_t(3, 0, 2.5);

  // Meta-analytic likelihood: each yi has its known sampling SE.
  y ~ normal(mu, se_y);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N_study] log_lik_study = rep_vector(0.0, N_study);

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], se_y[n]);
    log_lik_study[study[n]] += log_lik[n];

    // Posterior-predictive replicated effect size.
    y_rep[n] = normal_rng(mu[n], se_y[n]);
  }
}
