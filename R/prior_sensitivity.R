# Local prior perturbations for the frozen 2026 relative-intake primary fits.
prior_sensitivity_bundle_files <- function() {
  root <- file.path("results", "verification", "primary_update_2026")
  c(file.path("config", "verification-manifest-index.csv"),
    file.path(root, c("manifest.csv", "accepted_fit_manifest.csv.gz",
                     "stacking_weights.csv.gz", "analysis_data.csv.gz")))
}

read_prior_sensitivity_inputs <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Restore the frozen primary verification bundle. Missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  index <- readr::read_csv(paths[[1L]], show_col_types = FALSE)
  manifest <- readr::read_csv(paths[[2L]], show_col_types = FALSE)
  sha256 <- function(path) digest::digest(file = path, algo = "sha256")
  expected_manifest <- index$manifest_sha256[index$bundle == "primary_update_2026"]
  expected_files <- manifest$sha256[match(basename(paths[-c(1L, 2L)]), manifest$file)]
  if (!identical(sha256(paths[[2L]]), expected_manifest) ||
      !identical(vapply(paths[-c(1L, 2L)], sha256, character(1), USE.NAMES = FALSE),
                 expected_files)) {
    stop("Primary verification inputs do not match the frozen checksums.",
         call. = FALSE)
  }
  list(
    accepted = readr::read_csv(paths[[3L]], show_col_types = FALSE),
    weights = readr::read_csv(paths[[4L]], show_col_types = FALSE),
    analysis_data = readr::read_csv(
      paths[[5L]], show_col_types = FALSE,
      col_select = c(analysis_outcome, study_no, arm_no, outcome_no, study_factor,
                     yi, vi, protein_median, protein_intake_raw, protein_intake,
                     age_centre, duration_centre, resistance_exercise_code,
                     resistance_exercise_yes_no, resistance_exercise_factor)
    ),
    source_commit = unique(manifest$git_commit)
  )
}

make_prior_sensitivity_candidates <- function(inputs) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  keys <- c("outcome", "model_index")
  if (anyDuplicated(inputs$weights[keys]) || anyDuplicated(inputs$accepted[keys]) ||
      length(inputs$source_commit) != 1L) {
    stop("Invalid frozen candidate identities.", call. = FALSE)
  }
  weights <- dplyr::left_join(inputs$weights, inputs$accepted, by = keys)
  selected <- lapply(outcomes, function(outcome_name) {
    candidates <- weights[weights$outcome == outcome_name, ]
    candidates <- candidates[order(candidates$model_index), ]
    if (!identical(as.integer(candidates$model_index), seq_len(26L)) ||
        any(!is.finite(candidates$stacking_weight)) ||
        any(candidates$stacking_weight < 0) ||
        abs(sum(candidates$stacking_weight) - 1) > 1e-8 || anyNA(candidates$fit_path)) {
      stop("Incomplete primary candidate registry for ", outcome_name, call. = FALSE)
    }
    # Keep the archived relative path, regardless of the original drive or clone.
    paths <- gsub("\\", "/", candidates$fit_path, fixed = TRUE)
    paths <- sub("^.*/model_archive/", "model_archive/", paths)
    if (any(!grepl(
      "^model_archive/update_2026/primary/(candidates|remedies/remedy_1)/[a-z_]+/model_[0-9]{2}\\.qs$",
      paths
    ))) {
      stop("Unexpected accepted primary fit archive path.", call. = FALSE)
    }
    candidates$fit_path <- paths
    candidates$source_commit <- inputs$source_commit
    candidates
  })
  dplyr::bind_rows(selected)
}

prior_sensitivity_fit_path <- function(record) {
  path <- record$fit_path
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Restore this accepted primary fit from the research archive: ", path,
         call. = FALSE)
  }
  path
}

prior_sensitivity_contrasts <- function() {
  # Three manuscript/GRADE comparisons plus examples of 0.2 and 0.5 increments.
  tibble::tibble(
    current = c(0.8, 0.8, 0.8, 1.2, 1.6, 1.6, 1.6),
    proposed = c(1.0, 1.6, 2.4, 1.7, 1.8, 2.1, 2.4)
  )
}

prior_sensitivity_prediction_data <- function(data) {
  data$resistance_exercise_factor <- factor(data$resistance_exercise_factor)
  grid <- make_primary_prediction_grid(data)
  pairs <- merge(prior_sensitivity_contrasts(),
                 data.frame(resistance_exercise_code = c(0, 1)))
  pairs <- pairs[order(pairs$resistance_exercise_code, pairs$current, pairs$proposed), ]
  pairs$rt <- ifelse(pairs$resistance_exercise_code == 0, "RT", "No RT")
  support <- split(data$protein_intake_raw, data$resistance_exercise_code)
  pairs$observed_min <- vapply(pairs$resistance_exercise_code, function(code) {
    min(support[[as.character(code)]])
  }, numeric(1))
  pairs$observed_max <- vapply(pairs$resistance_exercise_code, function(code) {
    max(support[[as.character(code)]])
  }, numeric(1))
  pairs$supported <- pairs$current >= pairs$observed_min &
    pairs$proposed <= pairs$observed_max
  pairs$label <- sprintf("%.1f to %.1f", pairs$current, pairs$proposed)

  dose_key <- function(code, intake) sprintf("%d:%.2f", code, intake)
  keys <- dose_key(grid$keys$resistance_exercise_code, grid$keys$protein_intake)
  low <- match(dose_key(pairs$resistance_exercise_code, pairs$current), keys)
  high <- match(dose_key(pairs$resistance_exercise_code, pairs$proposed), keys)
  if (anyNA(c(low[pairs$supported], high[pairs$supported]))) {
    stop("A supported dose contrast is absent from the primary prediction grid.",
         call. = FALSE)
  }
  columns <- sort(unique(c(low[pairs$supported], high[pairs$supported])))
  list(pairs = pairs, newdata = grid$newdata[columns, , drop = FALSE],
       low = match(low, columns), high = match(high, columns))
}

prior_sensitivity_log_densities <- function(model) {
  coefficients <- c("Intercept", "resistance_exercise_code")
  variables <- c("lprior", paste0("b_", coefficients))
  if (!identical(model$backend, "rstan") ||
      !all(variables %in% model$fit@sim$fnames_oi)) {
    stop("The accepted rstan fit is missing required saved prior draws.",
         call. = FALSE)
  }
  # Non-permuted Stan draws are flattened iteration first, then chain, in
  # the same order as brms prediction draw IDs.
  draws <- rstan::extract(model$fit, pars = variables,
                          permuted = FALSE, inc_warmup = FALSE)
  prior_table <- as.data.frame(brms::prior_summary(model))
  anchor <- rep(0, prod(dim(draws)[1:2]))
  for (coefficient in coefficients) {
    prior <- prior_table$prior[prior_table$class == "b" & prior_table$coef == coefficient]
    if (length(prior) != 1L) {
      stop("Expected one explicit historical prior for ", coefficient, call. = FALSE)
    }
    pars <- extract_student_t_prior(prior)
    theta <- as.vector(draws[, , paste0("b_", coefficient)])
    anchor <- anchor + stats::dt((theta - pars[["location"]]) / pars[["scale"]],
                                  df = pars[["df"]], log = TRUE) - log(pars[["scale"]])
  }
  result <- list(
    joint = as.vector(draws[, , "lprior"]),
    empirical = anchor,
    chain_id = rep(seq_len(dim(draws)[2L]), each = dim(draws)[1L])
  )
  if (any(!is.finite(c(result$joint, result$empirical)))) {
    stop("Non-finite log prior density in the accepted fit.", call. = FALSE)
  }
  result
}

prior_sensitivity_weights <- function(log_density, chain_id, power) {
  n <- length(log_density)
  k_limit <- min(0.7, 1 - 1 / log10(n))
  log_ratios <- (power - 1) * log_density
  if (all(log_ratios == log_ratios[[1L]])) {
    return(list(weights = rep(1 / n, n), pareto_k = NA_real_,
                k_limit = k_limit, effective_draws = n, reliable = TRUE))
  }
  r_eff <- loo::relative_eff(exp(log_ratios - max(log_ratios)), chain_id = chain_id)
  # PSIS tail warnings are captured in the per-model diagnostic table.
  ps <- loo::psis(log_ratios, r_eff = r_eff)
  k <- unname(ps$diagnostics$pareto_k[[1L]])
  effective <- unname(ps$diagnostics$n_eff[[1L]])
  list(weights = as.numeric(stats::weights(ps, log = FALSE, normalize = TRUE)),
       pareto_k = k, k_limit = k_limit, effective_draws = effective,
       reliable = is.finite(k) && k < k_limit && is.finite(effective) && effective >= 400)
}

summarise_prior_sensitivity_contrast <- function(draw, weights, sesoi = 0.10) {
  sorted <- order(draw)
  cumulative <- cumsum(weights[sorted])
  interval <- draw[sorted[findInterval(c(0.025, 0.975), cumulative,
                                       left.open = TRUE) + 1L]]
  probabilities <- c(gain = sum(weights[draw > sesoi]),
                     equivalence = sum(weights[abs(draw) <= sesoi]),
                     loss = sum(weights[draw < -sesoi]))
  state <- names(probabilities)[probabilities >= 0.95]
  c(as.list(c(mean = sum(weights * draw), low = interval[[1L]],
              high = interval[[2L]], probabilities)),
    list(conclusion = if (length(state)) state else "uncertain"))
}

write_prior_sensitivity_report <- function(results, candidates,
    directory = file.path("manuscript", "generated", "prior-sensitivity")) {
  tables <- lapply(c("diagnostics", "contrasts", "support"), function(name) {
    dplyr::bind_rows(lapply(results, `[[`, name))
  })
  names(tables) <- c("diagnostics", "contrasts", "support")
  tables$candidates <- candidates[c("outcome", "model_index", "model_name", "stage",
                                    "stacking_weight", "source_commit")]
  contrasts <- tables$contrasts
  baseline <- contrasts[contrasts$prior == "baseline", ]
  keys <- c("outcome", "model_index", "rt", "current", "proposed")
  changes <- dplyr::left_join(contrasts[contrasts$prior != "baseline", ],
    baseline[c(keys, "mean", "gain", "equivalence", "loss", "conclusion")],
    by = keys, suffix = c("", "_baseline"), relationship = "many-to-one")
  changes <- dplyr::mutate(changes,
    mean_change = mean - mean_baseline,
    gain_change = gain - gain_baseline,
    equivalence_change = equivalence - equivalence_baseline,
    loss_change = loss - loss_baseline,
    conclusion_changed = conclusion != conclusion_baseline)
  tables$contrasts <- changes

  labels <- c(whole_body = "Whole-body lean/fat-free mass",
              hypertrophy = "Muscle size or mass", strength = "Strength",
              `function` = "Physical function")
  headline <- changes[changes$current == 0.8 & changes$proposed == 1.6, ]
  headline$outcome <- factor(labels[headline$outcome], levels = unname(labels))
  headline$scenario <- factor(
    paste(ifelse(headline$prior == "joint", "Saved prior terms", "Empirical priors"),
          headline$power),
    levels = rev(c("Saved prior terms 0.8", "Saved prior terms 1.2",
                   "Empirical priors 0.8", "Empirical priors 1.2")))
  figure <- ggplot2::ggplot(headline, ggplot2::aes(
    x = mean_change, y = scenario, colour = prior
  )) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0, height = 0.2, seed = 1988),
      alpha = 0.65, size = 1.8, na.rm = TRUE
    ) +
    ggplot2::facet_wrap(ggplot2::vars(outcome, rt), ncol = 2) +
    ggplot2::scale_colour_manual(values = c(joint = "#2166AC", empirical = "#A63603"),
                                 guide = "none") +
    ggplot2::labs(x = "Shift in the posterior mean contrast (standardised units)", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  save_file <- function(filename, write) {
    path <- file.path(directory, filename)
    tmp <- tempfile(tmpdir = directory, fileext = paste0(".", tools::file_ext(path)))
    on.exit(unlink(tmp))
    write(tmp)
    fs::file_move(tmp, path)
    path
  }
  csvs <- vapply(names(tables), function(name) {
    save_file(paste0(name, ".csv"), function(path) readr::write_csv(tables[[name]], path))
  }, character(1), USE.NAMES = FALSE)
  c(csvs, save_file("headline-contrast.png", function(path) {
    ggplot2::ggsave(path, plot = figure, width = 10, height = 8, dpi = 160, bg = "white")
  }))
}

analyse_prior_sensitivity_candidate <- function(record, path, inputs) {
  fit <- qs::qread(path)
  if (!identical(fit$outcome, record$outcome) ||
      !identical(as.integer(fit$model_index), as.integer(record$model_index)) ||
      !identical(fit$model_name, record$model_name) || !inherits(fit$model, "brmsfit")) {
    stop("Fit identity does not match the frozen reporting record: ", path, call. = FALSE)
  }
  density <- prior_sensitivity_log_densities(fit$model)
  data <- inputs$analysis_data[inputs$analysis_data$analysis_outcome == record$outcome, ]
  if (nrow(data) != stats::nobs(fit$model)) {
    stop("Frozen analysis rows do not match the accepted fit: ", path, call. = FALSE)
  }
  prediction <- prior_sensitivity_prediction_data(data)
  n <- length(density$joint)
  draws <- brms::posterior_epred(fit$model, newdata = prediction$newdata,
                                 re_formula = NA, draw_ids = seq_len(n))
  if (!identical(dim(draws), c(n, nrow(prediction$newdata))) ||
      any(!is.finite(draws))) {
    stop("Posterior predictions do not match the contrast grid: ", path, call. = FALSE)
  }
  supported <- which(prediction$pairs$supported)
  contrasts <- draws[, prediction$high[supported], drop = FALSE] -
    draws[, prediction$low[supported], drop = FALSE]
  scenarios <- data.frame(
    prior = c("baseline", "joint", "joint", "empirical", "empirical"),
    power = c(1, 0.8, 1.2, 0.8, 1.2)
  )
  results <- lapply(seq_len(nrow(scenarios)), function(i) {
    scenario <- scenarios[i, , drop = FALSE]
    log_prior <- if (scenario$prior == "empirical") density$empirical else density$joint
    ps <- prior_sensitivity_weights(log_prior, density$chain_id, scenario$power)
    key <- cbind(record[c("outcome", "model_index")], scenario)
    diagnostics <- cbind(key, n_draws = n, as.data.frame(ps[names(ps) != "weights"]))
    comparison <- dplyr::bind_rows(lapply(seq_along(supported), function(j) {
      cbind(key, prediction$pairs[supported[[j]], c("rt", "current", "proposed", "label")],
            as.data.frame(summarise_prior_sensitivity_contrast(contrasts[, j], ps$weights)))
    }))
    if (!ps$reliable) {
      comparison[c("mean", "low", "high", "gain", "equivalence", "loss", "conclusion")] <- NA
    }
    list(diagnostics = diagnostics, contrasts = comparison)
  })
  list(
    diagnostics = dplyr::bind_rows(lapply(results, `[[`, "diagnostics")),
    contrasts = dplyr::bind_rows(lapply(results, `[[`, "contrasts")),
    support = cbind(record[rep(1L, nrow(prediction$pairs)), c("outcome", "model_index")],
                    prediction$pairs)
  )
}
