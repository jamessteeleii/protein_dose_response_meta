# Pilot for the two archived supplementation-strength smooth models only.
# Condition on the accepted spline design and marginalise independent Gaussian
# study-smooth, arm and outcome coefficients. Keep every remaining brms prior.
tst_marginal_stan_code <- function(code, model_index) {
  study_spline <- switch(as.character(model_index), `14` = 2L, `26` = 3L,
    stop("Unsupported marginal exact model."))
  population <- if (model_index == 14L) "X * b + Xs * bs + Zs_1_1 * s_1_1" else
    "X * b + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1"
  study_terms <- sprintf("Zs_%d_%d * s_%d_%d", study_spline, 1:3, study_spline, 1:3)
  code <- gsub("\r", "", code, fixed = TRUE)
  replace_once <- function(text, before, after) {
    hits <- gregexpr(before, text, fixed = TRUE)[[1L]]
    if (length(hits) != 1L || hits[[1L]] < 0L) stop("Unexpected saved Stan program layout.")
    sub(before, after, text, fixed = TRUE)
  }
  original_likelihood <- paste(c(
    "  if (!prior_only) {",
    "    // initialize linear predictor term",
    "    vector[N] mu = rep_vector(0.0, N);",
    paste0("    mu += ", paste(c(population, study_terms), collapse = " + "), ";"),
    "    for (n in 1:N) {",
    "      // add more terms to the linear predictor",
    "      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];",
    "    }",
    "    target += normal_lpdf(Y | mu, se);",
    "  }"), collapse = "\n")
  if (!grepl("real sigma = 0;", code, fixed = TRUE)) stop("Expected known sampling variances only.")
  scales <- sprintf("append_row(sds_%d, append_row(sd_1, sd_2))", study_spline)
  likelihood <- paste(c(
    "  if (!prior_only) {",
    sprintf("    vector[N] mu = %s;", population),
    sprintf("    vector[5] scales = %s;", scales),
    "    for (g in 1:tst_studies) {",
    "      int m = tst_train_n[g];",
    "      int start = tst_first[g];",
    "      array[m] int idx = tst_order[start:(start + m - 1)];",
    "      matrix[m, m] L = cholesky_decompose(tst_covariance(tst_kernels, scales, se2, idx));",
    "      target += multi_normal_cholesky_lpdf(Y[idx] | mu[idx], L);",
    "    }",
    "  }"), collapse = "\n")
  code <- replace_once(code, original_likelihood, likelihood)
  # Each removed independent standard-normal latent has its declaration,
  # deterministic scaled assignment, and unit-normal prior in the saved code.
  latent <- c(sprintf("zs_%d_%d", study_spline, 1:3), "z_1", "z_2")
  transformed <- c(sprintf("s_%d_%d", study_spline, 1:3), "r_1_1", "r_2_1")
  lines <- strsplit(code, "\n", fixed = TRUE)[[1L]]
  drop <- rep(FALSE, length(lines))
  for (i in seq_along(latent)) {
    raw_lines <- grepl(paste0("\\b", latent[i], "\\b"), lines)
    scaled_lines <- grepl(paste0("\\b", transformed[i], "\\b"), lines)
    if (sum(raw_lines) != 3L || sum(scaled_lines) != 2L ||
        sum(raw_lines & scaled_lines) != 1L ||
        sum(raw_lines & grepl("target += std_normal_lpdf", lines, fixed = TRUE)) != 1L ||
        any((raw_lines | scaled_lines) & grepl("lprior", lines, fixed = TRUE))) {
      stop("Unexpected dependence on a Gaussian coefficient slated for integration.")
    }
    drop <- drop | raw_lines | scaled_lines
  }
  code <- paste(lines[!drop], collapse = "\n")
  code <- replace_once(code, "functions {\n}", paste(c(
    "functions {",
    "  matrix tst_covariance(array[] matrix kernels, vector scales, vector se2, array[] int idx) {",
    "    int m = size(idx);",
    "    matrix[m, m] V = diag_matrix(se2[idx]);",
    "    for (k in 1:5) V += square(scales[k]) * kernels[k][idx, idx];",
    "    return V;",
    "  }",
    "}"), collapse = "\n"))
  code <- replace_once(code, "  int prior_only;", paste(c(
    "  int<lower=1> tst_studies;",
    "  array[tst_studies] int<lower=1, upper=N> tst_first;",
    "  array[tst_studies] int<lower=1, upper=N> tst_train_n;",
    "  array[N] int<lower=1, upper=N> tst_order;",
    "  int<lower=1, upper=tst_studies> tst_held_study;",
    "  array[5] matrix[N, N] tst_kernels;",
    "  int prior_only;"), collapse = "\n"))
  replace_once(code, "generated quantities {\n}", paste(c(
    "generated quantities {",
    "  real tst_held_log_lik;",
    "  {",
    "    int g = tst_held_study;",
    "    int m = tst_train_n[g] + 1;",
    "    int start = tst_first[g];",
    "    array[m] int idx = tst_order[start:(start + m - 1)];",
    sprintf("    vector[N] mu = %s;", population),
    sprintf("    vector[5] scales = %s;", scales),
    "    matrix[m, m] L = cholesky_decompose(tst_covariance(tst_kernels, scales, se2, idx));",
    "    vector[m] z = mdivide_left_tri_low(L, Y[idx] - mu[idx]);",
    "    tst_held_log_lik = std_normal_lpdf(z[m]) - log(L[m, m]);",
    "  }",
    "}"), collapse = "\n"))
}

tst_prepare_marginal_exact <- function(record) {
  stopifnot(record$family == "supplement", record$outcome == "strength",
    record$model_index %in% c(14L, 26L), record$observation %in% c(80L, 81L))
  iad_check_hash(record$fit_path, record$fit_sha256)
  iad_check_hash(record$loo_path, record$loo_sha256)
  candidate <- qs::qread(record$fit_path)
  integrated <- qs::qread(record$loo_path)
  model <- candidate$model
  manifest <- integrated$observation_manifest
  obs <- as.integer(record$observation)
  stopifnot(candidate$model_index == record$model_index, candidate$outcome == record$outcome,
    candidate$model_name == record$model_name, integrated$fit_sha256 == record$fit_sha256,
    integrated$study_structure == "study_factor_smooth", manifest$study_no[obs] == 57,
    manifest$arm_no[obs] == 85, manifest$outcome_no[obs] == 461,
    loo::pareto_k_values(integrated$loo)[obs] > .7, iad_same_observations(model$data, manifest))
  d <- brms::standata(model)
  group <- if (record$model_index == 14L) 2L else 3L
  design_names <- sprintf("Zs_%d_%d", group, 1:3)
  designs <- d[design_names]
  stopifnot(d$N == nrow(manifest), d$M_1 == 1L, d$M_2 == 1L, d[[paste0("nb_", group)]] == 3L,
    all(d$Z_1_1 == 1), all(d$Z_2_1 == 1), all(d$se > 0), all(is.finite(d$se)),
    isTRUE(all.equal(as.numeric(d$Y), manifest$yi, tolerance = 1e-12)),
    isTRUE(all.equal(as.numeric(d$se)^2, as.numeric(manifest$vi), tolerance = 1e-12)))
  arm_kernel <- outer(as.integer(d$J_1), as.integer(d$J_1), `==`) * 1
  outcome_kernel <- outer(as.integer(d$J_2), as.integer(d$J_2), `==`) * 1
  stopifnot(identical(arm_kernel, outer(manifest$arm_no, manifest$arm_no, `==`) * 1),
    identical(outcome_kernel, outer(manifest$outcome_no, manifest$outcome_no, `==`) * 1))
  kernels <- c(lapply(designs, tcrossprod), list(arm_kernel, outcome_kernel))
  cross_study <- outer(manifest$study_no, manifest$study_no, `!=`)
  if (any(vapply(kernels, function(k) any(k[cross_study] != 0), logical(1))))
    stop("Study blocks are not independent in the saved design.")
  groups <- split(seq_len(d$N), manifest$study_no)
  held <- which(vapply(groups, function(x) obs %in% x, logical(1)))
  groups[[held]] <- c(setdiff(groups[[held]], obs), obs)
  d$tst_studies <- length(groups)
  d$tst_first <- as.integer(cumsum(c(1L, head(lengths(groups), -1L))))
  d$tst_train_n <- lengths(groups)
  d$tst_train_n[[held]] <- d$tst_train_n[[held]] - 1L
  d$tst_order <- unname(unlist(groups))
  d$tst_held_study <- held
  d$tst_kernels <- aperm(array(unlist(kernels), dim = c(d$N, d$N, 5L)), c(3L, 1L, 2L))
  stopifnot(all(d$tst_train_n > 0L), sum(d$tst_train_n) == d$N - 1L,
    d$tst_order[d$tst_first[held] + d$tst_train_n[held]] == obs)
  # Compare the held-out Gaussian score against the production integrated
  # scorer over representative scale settings on the actual archived design.
  bases <- lapply(designs, map_primary_factor_smooth_design, study_factor = model$data$study_factor)
  score_errors <- vapply(c(.05, .3, 1), function(scale) {
    scales <- scale * c(1, 1.3, .7, .8, .5)
    mu <- as.numeric(d$X %*% seq(.1, .5, length.out = d$K))
    basis <- unlist(Map(function(b, s) lapply(b$component_values, function(v) matrix(v * s, nrow = 1L)),
      bases, scales[1:3]), recursive = FALSE)
    reference <- integrate_primary_study_basis_log_lik(matrix(mu, nrow = 1L), d$Y, d$se^2,
      manifest$study_no, manifest$arm_no, manifest$outcome_no, 0, scales[4], scales[5], basis)[1L, obs]
    idx <- groups[[held]]
    V <- diag(d$se[idx]^2) + Reduce(`+`, Map(function(k, s) k[idx, idx] * s^2, kernels, scales))
    L <- t(chol(V)); last <- length(idx)
    z <- forwardsolve(L, d$Y[idx] - mu[idx])
    abs(reference - (stats::dnorm(z[last], log = TRUE) - log(L[last, last])))
  }, numeric(1))
  if (any(!is.finite(score_errors)) || max(score_errors) > 1e-8)
    stop("Marginal score differs from the existing integrated scorer.")
  code <- tst_marginal_stan_code(brms::stancode(model), record$model_index)
  list(code = code, data = d, score_error = max(score_errors), observation = obs,
    record = record, old_k = loo::pareto_k_values(integrated$loo)[obs])
}

tst_marginal_pilot_diagnostics <- function(fit) {
  raw_draws <- as.array(fit)
  draws <- posterior::as_draws_array(raw_draws)
  monitored <- posterior::variables(draws)
  monitored <- setdiff(monitored, c("lp__", "lprior", "sigma"))
  summaries <- posterior::summarise_draws(posterior::subset_draws(draws, variable = monitored),
    "rhat", "ess_bulk", "ess_tail")
  nuts <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(nuts, function(x) sum(x[, "divergent__"]), numeric(1)))
  depth_hits <- sum(vapply(nuts, function(x) sum(x[, "treedepth__"] >= 15L), integer(1)))
  ebfmi <- min(vapply(nuts, function(x) mean(diff(x[, "energy__"])^2) / stats::var(x[, "energy__"]), numeric(1)))
  geometry_pass <- divergences == 0 && depth_hits == 0 && is.finite(ebfmi) && ebfmi >= .30
  convergence_pass <- all(is.finite(as.matrix(summaries[, c("rhat", "ess_bulk", "ess_tail")]))) &&
    max(summaries$rhat) <= 1.01 && min(summaries$ess_bulk) >= 400 && min(summaries$ess_tail) >= 400
  ll <- raw_draws[, , "tst_held_log_lik"]
  peak <- max(ll); likelihood <- exp(ll - peak)
  list(diagnostics = data.frame(divergent_transitions = divergences,
    maximum_treedepth_transitions = depth_hits, minimum_ebfmi = ebfmi,
    maximum_rhat = max(summaries$rhat), minimum_bulk_ess = min(summaries$ess_bulk),
    minimum_tail_ess = min(summaries$ess_tail), geometry_pass = geometry_pass,
    convergence_pass = convergence_pass, diagnostic_pass = geometry_pass && convergence_pass,
    exact_elpd_loo = peak + log(mean(likelihood)),
    exact_elpd_mcse = posterior::mcse_mean(likelihood) / mean(likelihood)), held_log_lik = ll)
}

tst_fit_marginal_exact <- function(input, iter_sampling, syntax = NULL) {
  stopifnot(iter_sampling %in% c(6000L, 12000L))
  if (is.null(syntax)) syntax <- rstan::stanc(model_code = input$code,
    model_name = paste0("training_marginal_exact_", input$record$model_index))
  stopifnot(isTRUE(syntax$status))
  stan_model <- rstan::stan_model(stanc_ret = syntax, auto_write = FALSE)
  started <- proc.time()[["elapsed"]]
  fit <- rstan::sampling(stan_model, data = input$data, chains = 4L, cores = 4L,
    iter = 2000L + iter_sampling, warmup = 2000L, seed = 1988L,
    control = list(adapt_delta = .999, max_treedepth = 15L))
  result <- tst_marginal_pilot_diagnostics(fit)
  result$diagnostics$refit_seconds <- proc.time()[["elapsed"]] - started
  result$diagnostics$model_index <- input$record$model_index
  result$diagnostics$observation <- input$observation
  result$diagnostics$iter_sampling <- iter_sampling
  result$diagnostics$old_psis_pareto_k <- input$old_k
  result$interaction_record <- input$record
  result$code_sha256 <- digest::digest(input$code, algo = "sha256")
  result$score_error <- input$score_error
  result
}

tst_marginal_compact <- function(result, record, integrated) {
  d <- result$diagnostics
  ll <- result$held_log_lik
  # The original pilot used as.matrix(draws_array), which stacks each chain
  # consecutively. Restore its iterations-by-chains layout before computing MCSE.
  if (record$model_index == 26L && record$observation == 80L && d$iter_sampling == 6000L &&
      identical(dim(ll), c(24000L, 1L))) ll <- matrix(ll, nrow = 6000L, ncol = 4L)
  stopifnot(identical(result$interaction_record, record), nrow(d) == 1L,
    d$model_index == record$model_index, d$observation == record$observation,
    identical(dim(ll), c(as.integer(d$iter_sampling), 4L)), all(is.finite(ll)),
    is.finite(result$score_error), result$score_error <= 1e-8,
    integrated$model_index == record$model_index, integrated$fit_sha256 == record$fit_sha256,
    isTRUE(all.equal(d$old_psis_pareto_k, loo::pareto_k_values(integrated$loo)[record$observation])))
  peak <- max(ll)
  likelihood <- exp(ll - peak)
  stopifnot(isTRUE(all.equal(d$exact_elpd_loo, peak + log(mean(likelihood)), tolerance = 1e-12)))
  d$exact_elpd_mcse <- posterior::mcse_mean(likelihood) / mean(likelihood)
  stopifnot(is.finite(d$exact_elpd_mcse))
  d$family <- record$family
  d$outcome <- record$outcome
  d$model_name <- record$model_name
  d$exact_stage <- "marginal_gaussian"
  d$iter_warmup <- 2000L
  d$adapt_delta <- .999
  d$max_treedepth <- 15L
  d$diagnostic_pass <- isTRUE(d$diagnostic_pass) && isTRUE(d$geometry_pass) && isTRUE(d$convergence_pass) &&
    d$divergent_transitions == 0L && d$maximum_treedepth_transitions == 0L && d$minimum_ebfmi >= .30 &&
    d$maximum_rhat <= 1.01 && d$minimum_bulk_ess >= 400 && d$minimum_tail_ess >= 400
  correction <- NULL
  if (isTRUE(d$diagnostic_pass)) {
    obs <- as.integer(record$observation)
    correction <- integrated
    original_elpd <- correction$loo$pointwise[obs, "elpd_loo"]
    d$old_psis_elpd_loo <- original_elpd
    d$exact_minus_psis_elpd <- d$exact_elpd_loo - original_elpd
    correction$loo$pointwise[obs, "elpd_loo"] <- d$exact_elpd_loo
    correction$loo$pointwise[obs, "p_loo"] <- correction$loo$pointwise[obs, "p_loo"] - d$exact_minus_psis_elpd
    correction$loo$pointwise[obs, "looic"] <- -2 * d$exact_elpd_loo
    if ("mcse_elpd_loo" %in% colnames(correction$loo$pointwise))
      correction$loo$pointwise[obs, "mcse_elpd_loo"] <- d$exact_elpd_mcse
    correction$loo$diagnostics$pareto_k[obs] <- NA_real_
    if (!is.null(correction$loo$diagnostics$n_eff)) correction$loo$diagnostics$n_eff[obs] <- NA_real_
    correction$loo <- update_primary_exact_loo_estimates(correction$loo)
    correction$exact_refit_summary <- d
  }
  list(compact = TRUE, interaction_record = record, diagnostics = d, correction = correction,
    code_sha256 = result$code_sha256, score_error = result$score_error)
}

tst_marginal_exact_checkpoint <- function(plan, pilot_file) {
  iad_require_phase("exact_marginal")
  if (is.null(plan$remedy)) stop("Marginal retry requires a reviewed failed case.")
  record <- plan$record
  iad_check_hash(plan$initial_path, plan$initial_sha256)
  name <- sprintf("model_%02d_observation_%04d.qs", record$model_index, record$observation)
  path <- iad_path(record$family, record$outcome, "loo", "exact_checkpoints", "marginal_gaussian", name, archive = TRUE)
  input <- tst_prepare_marginal_exact(record)
  code_sha256 <- digest::digest(input$code, algo = "sha256")
  if (file.exists(path)) {
    saved <- qs::qread(path)
    stopifnot(identical(saved$interaction_record, record), identical(saved$code_sha256, code_sha256),
      identical(saved$initial_checkpoint_sha256, plan$initial_sha256))
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  reuse_pilot <- record$model_index == 26L && record$observation == 80L
  result <- if (reuse_pilot) qs::qread(pilot_file) else
    tst_fit_marginal_exact(input, if (record$model_index == 14L) 12000L else 6000L)
  stopifnot(identical(result$code_sha256, code_sha256))
  saved <- tst_marginal_compact(result, record, qs::qread(record$loo_path))
  if (reuse_pilot && !isTRUE(saved$diagnostics$diagnostic_pass)) stop("The marginal pilot has not passed.")
  saved$initial_checkpoint_sha256 <- plan$initial_sha256
  saved$diagnostics$initial_checkpoint_sha256 <- plan$initial_sha256
  saved$diagnostics$reused_marginal_pilot <- reuse_pilot
  iad_qsave(saved, path)
}

tst_select_marginal_exact_paths <- function(initial_paths, replacement_paths) {
  key <- function(x) {
    r <- x$interaction_record
    paste(r$family, r$outcome, r$model_index, r$observation, sep = "|")
  }
  initial_paths <- unname(unlist(initial_paths)); replacement_paths <- unname(unlist(replacement_paths))
  initial <- lapply(initial_paths, qs::qread); replacements <- lapply(replacement_paths, qs::qread)
  keys <- vapply(initial, key, character(1)); new_keys <- vapply(replacements, key, character(1))
  failed <- !vapply(initial, function(x) isTRUE(x$diagnostics$diagnostic_pass), logical(1))
  if (anyDuplicated(keys) || anyDuplicated(new_keys) || !setequal(keys[failed], new_keys))
    stop("Marginal replacement keys do not match precisely the failed exact cases.")
  for (i in seq_along(replacements)) {
    j <- match(new_keys[i], keys)
    stopifnot(identical(initial[[j]]$interaction_record, replacements[[i]]$interaction_record),
      identical(iad_hash(initial_paths[j]), replacements[[i]]$initial_checkpoint_sha256))
    initial_paths[j] <- replacement_paths[i]
  }
  initial_paths
}
