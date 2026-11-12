data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se_y;
  vector[N] protein_intake;
  vector[N] resistance_exercise_code;
  vector[N] age_z;
  vector[N] duration_z;
  int<lower=0, upper=1> include_moderators;
  int<lower=0, upper=1> include_rt_interaction;

  int<lower=1> N_study;
  array[N] int<lower=1, upper=N_study> study;
  int<lower=1> N_arm;
  array[N] int<lower=1, upper=N_arm> arm;
  int<lower=1> N_outcome;
  array[N] int<lower=1, upper=N_outcome> outcome;

  int<lower=1> N_cluster;
  array[N_cluster] int<lower=1, upper=N> cluster_start;
  array[N_cluster] int<lower=1, upper=N> cluster_end;
  real<lower=-1, upper=0> rho_lower;
  real rho_raw_prior_location;

  real prior_intercept_location;
  real<lower=0> prior_intercept_scale;
  real prior_rt_location;
  real<lower=0> prior_rt_scale;
}

parameters {
  real b_intercept;
  real b_rt;
  real b_protein;
  real b_protein_rt;
  vector[2] b_moderator;

  vector[N_study] z_study;
  vector[N_arm] z_arm;
  vector[N_outcome] z_outcome;
  real<lower=0> sd_study;
  real<lower=0> sd_arm;
  real<lower=0> sd_outcome;
  real rho_raw;
}

transformed parameters {
  real<lower=rho_lower, upper=1> rho =
    rho_lower + (1 - rho_lower) * inv_logit(rho_raw);
  vector[N] mu;
  for (n in 1:N) {
    real eta = b_intercept +
      b_rt * resistance_exercise_code[n] +
      b_protein * protein_intake[n] +
      sd_study * z_study[study[n]] +
      sd_arm * z_arm[arm[n]] +
      sd_outcome * z_outcome[outcome[n]];
    if (include_rt_interaction == 1) {
      eta += b_protein_rt * protein_intake[n] *
        resistance_exercise_code[n];
    }
    if (include_moderators == 1) {
      eta += b_moderator[1] * duration_z[n] +
        b_moderator[2] * age_z[n];
    }
    mu[n] = eta;
  }
}

model {
  b_intercept ~ student_t(
    3, prior_intercept_location, prior_intercept_scale
  );
  b_rt ~ student_t(3, prior_rt_location, prior_rt_scale);
  b_protein ~ student_t(3, 0, 1);
  b_protein_rt ~ student_t(3, 0, 1);
  b_moderator ~ student_t(3, 0, 0.5);

  z_study ~ std_normal();
  z_arm ~ std_normal();
  z_outcome ~ std_normal();
  sd_study ~ student_t(3, 0, 2.5);
  sd_arm ~ student_t(3, 0, 2.5);
  sd_outcome ~ student_t(3, 0, 2.5);
  rho_raw ~ normal(rho_raw_prior_location, 1.5);

  for (j in 1:N_cluster) {
    int first = cluster_start[j];
    int last = cluster_end[j];
    int size = last - first + 1;
    vector[size] y_j = segment(y, first, size);
    vector[size] mu_j = segment(mu, first, size);
    vector[size] se_j = segment(se_y, first, size);
    matrix[size, size] covariance = diag_matrix(square(se_j));
    if (size > 1) {
      for (row in 1:(size - 1)) {
        for (column in (row + 1):size) {
          covariance[row, column] = rho * se_j[row] * se_j[column];
          covariance[column, row] = covariance[row, column];
        }
      }
    }
    y_j ~ multi_normal_cholesky(mu_j, cholesky_decompose(covariance));
  }
}

generated quantities {
  vector[N_cluster] log_lik_cluster;
  vector[N] y_rep;
  for (j in 1:N_cluster) {
    int first = cluster_start[j];
    int last = cluster_end[j];
    int size = last - first + 1;
    vector[size] y_j = segment(y, first, size);
    vector[size] mu_j = segment(mu, first, size);
    vector[size] se_j = segment(se_y, first, size);
    matrix[size, size] covariance = diag_matrix(square(se_j));
    matrix[size, size] cholesky;
    if (size > 1) {
      for (row in 1:(size - 1)) {
        for (column in (row + 1):size) {
          covariance[row, column] = rho * se_j[row] * se_j[column];
          covariance[column, row] = covariance[row, column];
        }
      }
    }
    cholesky = cholesky_decompose(covariance);
    log_lik_cluster[j] = multi_normal_cholesky_lpdf(y_j | mu_j, cholesky);
    y_rep[first:last] = multi_normal_cholesky_rng(mu_j, cholesky);
  }
}
