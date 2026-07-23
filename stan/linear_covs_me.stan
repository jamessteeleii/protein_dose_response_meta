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
  vector[N] duration_centre;
  vector[N] age_centre;

  // Values from your preregistered/Steele-derived priors.
  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;
  real prior_rt_location;
  real<lower=0> prior_rt_scale;
}

parameters {
  // Latent, true protein intake for every study-specific arm.
  vector[N_arm] x_true;
  real x_mean;
  real<lower=0> x_sd;

  real b_intercept;
  real b_rt;
  real b_protein;
  real b_duration;
  real b_age;

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
      b_duration * duration_centre[n] +
      b_age * age_centre[n] +
      sd_study * z_study[study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
  }
}

model {
  // Distribution for latent arm-level protein intakes.
  // protein_intake is already centred, hence x_mean is centred at zero.
  x_mean ~ normal(0, 2);
  x_sd ~ exponential(1);
  x_true ~ normal(x_mean, x_sd);

  // Known sampling error of protein intake.
  x_obs ~ normal(x_true, x_se);

  // Fixed-effect priors.
  b_intercept ~ student_t(
    3,
    prior_intercept_location,
    prior_intercept_scale
  );
  b_rt ~ student_t(3, prior_rt_location, prior_rt_scale);

  b_protein ~ student_t(3, 0, 10);
  b_duration ~ student_t(3, 0, 10);
  b_age ~ student_t(3, 0, 10);

  // Non-centred random-effect parameterisation.
  z_study ~ std_normal();
  z_arm ~ std_normal();
  z_outcome ~ std_normal();

  sd_study ~ student_t(3, 0, 2.5);
  sd_arm ~ student_t(3, 0, 2.5);
  sd_outcome ~ student_t(3, 0, 2.5);

  // Effect-size likelihood: known sampling SE plus hierarchical variation.
  y ~ normal(mu, se_y);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], se_y[n]);
  }
}


