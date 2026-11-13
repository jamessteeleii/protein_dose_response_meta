# Final observation-level LOO implementation for the primary candidate ensemble.
# Conditional PSIS supplies the common observation manifest, after which the
# Gaussian multilevel terms are integrated and the few remaining high Pareto-k
# cases are replaced by exact case-deletion refits.

# Align fitted observations with source rows ----

# Define the observation keys and outcome fields used to match LOO scores to
# the same extracted records across candidate models.
primary_loo_manifest_columns <- function(data) {
  study_column <- if ("study_no" %in% names(data)) {
    "study_no"
  } else if ("study_factor" %in% names(data)) {
    "study_factor"
  } else {
    stop(
      "Primary LOO data are missing both study_no and study_factor.",
      call. = FALSE
    )
  }

  required <- c("arm_no", "outcome_no", "yi", "vi")
  missing <- setdiff(required, names(data))

  if (length(missing) > 0L) {
    stop(
      "Primary LOO data are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  optional <- c(
    "outcome_group", "measurement_method", "unit_of_measurement",
    "time1", "time2"
  )

  c(study_column, required, intersect(optional, names(data)))
}

get_primary_loo_observation_manifest <- function(fitted_candidate) {
  required_fields <- c("outcome", "model_index", "model_name", "model")
  missing_fields <- setdiff(required_fields, names(fitted_candidate))

  if (length(missing_fields) > 0L) {
    stop(
      "The fitted primary candidate is missing: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  model <- fitted_candidate$model
  if (!inherits(model, "brmsfit") || is.null(model$data)) {
    stop("The primary candidate must contain a brmsfit with model data.", call. = FALSE)
  }

  data <- model$data
  columns <- primary_loo_manifest_columns(data)
  manifest <- tibble::as_tibble(data[, columns, drop = FALSE])

  # Factors can acquire inconsequential attribute differences across formulae.
  # Store their labels so comparisons concern the observations themselves.
  manifest <- dplyr::mutate(
    manifest,
    dplyr::across(where(is.factor), as.character)
  )
  if (!"study_no" %in% names(manifest)) {
    manifest <- dplyr::rename(manifest, study_no = study_factor)
  }

  manifest <- dplyr::mutate(manifest, study_no = as.character(study_no))

  if (anyNA(manifest$study_no) || any(!nzchar(manifest$study_no))) {
    stop("The primary candidate has invalid study identifiers.", call. = FALSE)
  }

  manifest <- dplyr::mutate(
    manifest,
    outcome = as.character(fitted_candidate$outcome),
    observation = dplyr::row_number(),
    .before = 1L
  )

  if (nrow(manifest) < 1L || anyNA(manifest$yi) || anyNA(manifest$vi)) {
    stop("The primary candidate has an invalid LOO observation manifest.", call. = FALSE)
  }

  manifest
}

# Conditional PSIS-LOO screen ----

# Score the fit's stored log likelihood and attach the observation manifest.
# Integrated scoring below accounts for the held-out observation's random
# effects.
compute_primary_candidate_psis_loo <- function(
    fitted_candidate,
    cores = 1L,
    save_psis = FALSE
) {
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores) || cores < 1L) {
    stop("cores must be one positive integer.", call. = FALSE)
  }

  manifest <- get_primary_loo_observation_manifest(fitted_candidate)
  loo_result <- brms::loo(
    fitted_candidate$model,
    pointwise = TRUE,
    moment_match = FALSE,
    reloo = FALSE,
    save_psis = save_psis,
    cores = as.integer(cores)
  )

  if (
    !inherits(loo_result, "loo") ||
      is.null(loo_result$pointwise) ||
      is.null(loo_result$diagnostics$pareto_k) ||
      nrow(loo_result$pointwise) != nrow(manifest) ||
      length(loo_result$diagnostics$pareto_k) != nrow(manifest)
  ) {
    stop("The pointwise PSIS-LOO result failed validation.", call. = FALSE)
  }

  list(
    outcome = as.character(fitted_candidate$outcome),
    model_index = as.integer(fitted_candidate$model_index),
    model_name = as.character(fitted_candidate$model_name),
    observation_manifest = manifest,
    loo = loo_result,
    remedy = "none_initial_psis"
  )
}

primary_loo_estimate <- function(loo_result, statistic, column) {
  estimates <- loo_result$estimates
  if (
    is.null(estimates) ||
      !statistic %in% rownames(estimates) ||
      !column %in% colnames(estimates)
  ) {
    return(NA_real_)
  }

  unname(estimates[statistic, column])
}












# Integrate Gaussian random effects for holdout scoring ----

# Accept draw-by-observation predictors and residual variance components and
# return one conditional log density per draw and observation. The held-out
# row is conditioned on the remaining rows in its study.
integrate_primary_study_basis_log_lik <- function(
    linear_predictor,
    observed,
    known_variance,
    study_id,
    arm_id,
    outcome_id,
    residual_sd,
    arm_intercept_sd,
    outcome_intercept_sd,
    study_basis
) {
  linear_predictor <- as.matrix(linear_predictor)
  observed <- as.numeric(observed)
  known_variance <- as.numeric(known_variance)
  study_id <- as.character(study_id)
  arm_id <- as.character(arm_id)
  outcome_id <- as.character(outcome_id)
  residual_sd <- as.numeric(residual_sd)
  arm_intercept_sd <- as.numeric(arm_intercept_sd)
  outcome_intercept_sd <- as.numeric(outcome_intercept_sd)
  study_basis <- lapply(study_basis, as.matrix)

  n_draws <- nrow(linear_predictor)
  n_observations <- ncol(linear_predictor)
  observation_vectors <- list(
    observed,
    known_variance,
    study_id,
    arm_id,
    outcome_id
  )
  draw_vectors <- list(
    residual_sd,
    arm_intercept_sd,
    outcome_intercept_sd
  )
  arm_study_counts <- vapply(
    split(study_id, arm_id),
    function(x) length(unique(x)),
    integer(1)
  )
  outcome_arm_counts <- vapply(
    split(arm_id, outcome_id),
    function(x) length(unique(x)),
    integer(1)
  )
  invalid_basis <-
    length(study_basis) < 1L ||
    any(vapply(
      study_basis,
      function(x) {
        !identical(dim(x), c(n_draws, n_observations)) ||
          anyNA(x) || any(!is.finite(x))
      },
      logical(1)
    ))

  if (
    n_draws < 1L || n_observations < 1L ||
      any(vapply(
        observation_vectors,
        length,
        integer(1)
      ) != n_observations) ||
      any(vapply(draw_vectors, length, integer(1)) != n_draws) ||
      anyNA(linear_predictor) || anyNA(observed) ||
      anyNA(known_variance) || any(known_variance <= 0) ||
      anyNA(study_id) || any(!nzchar(study_id)) ||
      anyNA(arm_id) || any(!nzchar(arm_id)) ||
      anyNA(outcome_id) || any(!nzchar(outcome_id)) ||
      any(vapply(draw_vectors, anyNA, logical(1))) ||
      any(vapply(draw_vectors, function(x) any(x < 0), logical(1))) ||
      any(arm_study_counts != 1L) ||
      any(outcome_arm_counts != 1L) ||
      invalid_basis
  ) {
    stop(
      "The inputs to the integrated study-basis likelihood are invalid.",
      call. = FALSE
    )
  }

  arm_variance <- pmax(arm_intercept_sd^2, .Machine$double.eps)
  outcome_variance <- pmax(
    outcome_intercept_sd^2,
    .Machine$double.eps
  )
  log_lik <- matrix(NA_real_, nrow = n_draws, ncol = n_observations)
  study_groups <- split(seq_len(n_observations), study_id)

  for (indices in study_groups) {
    # Begin with sampling and residual variance, then account for the random
    # intercept shared by observations of the same outcome.
    group_variance <- outer(

      residual_sd^2,
      known_variance[indices],
      FUN = "+"
    )
    group_precision <- 1 / group_variance
    group_linear_predictor <- linear_predictor[, indices, drop = FALSE]
    group_residual <- sweep(
      group_linear_predictor,
      2L,
      observed[indices],
      FUN = function(eta, y) y - eta
    )
    n_in_study <- length(indices)

    outcome_precision_one <- matrix(
      NA_real_,
      nrow = n_draws,
      ncol = n_in_study
    )
    outcome_precision_diagonal <- matrix(
      NA_real_,
      nrow = n_draws,
      ncol = n_in_study
    )
    outcomes_within_study <- split(
      seq_len(n_in_study),
      outcome_id[indices]
    )

    apply_outcome_precision <- function(rhs) {
      output <- matrix(
        NA_real_,
        nrow = n_draws,
        ncol = n_in_study
      )
      for (positions in outcomes_within_study) {
        precision <- group_precision[, positions, drop = FALSE]
        total_precision <- rowSums(precision)
        denominator <- 1 / outcome_variance + total_precision
        weighted_rhs <- rowSums(
          rhs[, positions, drop = FALSE] * precision
        )
        for (position in positions) {
          observation_precision <- group_precision[, position]
          output[, position] <-
            observation_precision * rhs[, position] -
            observation_precision * weighted_rhs / denominator
        }
      }
      output
    }

    for (positions in outcomes_within_study) {
      precision <- group_precision[, positions, drop = FALSE]
      total_precision <- rowSums(precision)
      denominator <- 1 / outcome_variance + total_precision
      for (position in positions) {
        observation_precision <- group_precision[, position]
        outcome_precision_one[, position] <-
          observation_precision -
          observation_precision * total_precision / denominator
        outcome_precision_diagonal[, position] <-
          observation_precision -
          observation_precision^2 / denominator
      }
    }

    # Add the arm-level contribution after accounting for the outcome-level
    # covariance.
    arms_within_study <- split(
      seq_len(n_in_study),
      arm_id[indices]
    )
    arm_denominators <- lapply(
      arms_within_study,
      function(positions) {
        1 / arm_variance + rowSums(
          outcome_precision_one[, positions, drop = FALSE]
        )
      }
    )

    apply_arm_precision <- function(outcome_precision_rhs) {
      output <- outcome_precision_rhs
      arm_number <- 0L
      for (positions in arms_within_study) {
        arm_number <- arm_number + 1L
        rhs_total <- rowSums(
          outcome_precision_rhs[, positions, drop = FALSE]
        )
        denominator <- arm_denominators[[arm_number]]
        for (position in positions) {
          output[, position] <-
            outcome_precision_rhs[, position] -
            outcome_precision_one[, position] *
              rhs_total / denominator
        }
      }
      output
    }

    precision_residual <- apply_arm_precision(
      apply_outcome_precision(group_residual)
    )
    precision_diagonal <- outcome_precision_diagonal
    arm_number <- 0L
    for (positions in arms_within_study) {
      arm_number <- arm_number + 1L
      denominator <- arm_denominators[[arm_number]]
      for (position in positions) {
        precision_diagonal[, position] <-
          outcome_precision_diagonal[, position] -
          outcome_precision_one[, position]^2 / denominator
      }
    }

    group_basis <- lapply(
      study_basis,
      function(x) x[, indices, drop = FALSE]
    )
    precision_basis <- lapply(
      group_basis,
      function(x) apply_arm_precision(apply_outcome_precision(x))
    )

    # Add each independent standard-normal study basis component with a
    # vectorised Sherman-Morrison update. The basis matrices already include
    # their draw-specific scales and any Cholesky correlation transform.
    for (basis_index in seq_along(group_basis)) {
      basis <- group_basis[[basis_index]]
      precision_basis_current <- precision_basis[[basis_index]]
      denominator <-
        1 + rowSums(basis * precision_basis_current)
      residual_crossproduct <- rowSums(basis * precision_residual)
      precision_residual <-
        precision_residual -
        precision_basis_current * residual_crossproduct / denominator
      precision_diagonal <-
        precision_diagonal -
        precision_basis_current^2 / denominator

      if (basis_index < length(group_basis)) {
        for (later_index in seq.int(
          basis_index + 1L,
          length(group_basis)
        )) {
          basis_crossproduct <- rowSums(
            basis * precision_basis[[later_index]]
          )
          precision_basis[[later_index]] <-
            precision_basis[[later_index]] -
            precision_basis_current * basis_crossproduct / denominator
        }
      }
    }

    if (
      anyNA(precision_diagonal) ||
        any(!is.finite(precision_diagonal)) ||
        any(precision_diagonal <= 0)
    ) {
      stop(
        "The integrated study-basis conditional precision is invalid.",
        call. = FALSE
      )
    }

    # The diagonal precision gives the conditional variance. The
    # precision-weighted residual supplies the corresponding standardized
    # error.
    log_lik[, indices] <-
      0.5 * (log(precision_diagonal) - log(2 * pi)) -
      0.5 * precision_residual^2 / precision_diagonal
  }

  if (anyNA(log_lik) || any(!is.finite(log_lik))) {
    stop(
      "The integrated study-basis likelihood is non-finite.",
      call. = FALSE
    )
  }

  log_lik
}

summarise_primary_loo_method <- function(
    loo_result,
    method,
    k_threshold = 0.7
) {
  if (
    !inherits(loo_result, "loo") ||
      is.null(loo_result$diagnostics$pareto_k)
  ) {
    stop("Cannot summarise an invalid LOO result.", call. = FALSE)
  }

  pareto_k <- as.numeric(loo_result$diagnostics$pareto_k)
  tibble::tibble(
    method = as.character(method),
    elpd_loo = primary_loo_estimate(loo_result, "elpd_loo", "Estimate"),
    se_elpd_loo = primary_loo_estimate(loo_result, "elpd_loo", "SE"),
    p_loo = primary_loo_estimate(loo_result, "p_loo", "Estimate"),
    looic = primary_loo_estimate(loo_result, "looic", "Estimate"),
    maximum_pareto_k = max(pareto_k, na.rm = TRUE),
    n_pareto_k_gt_threshold = sum(
      pareto_k > k_threshold,
      na.rm = TRUE
    ),
    n_pareto_k_gt_1 = sum(pareto_k > 1, na.rm = TRUE)
  )
}





extract_primary_exact_draws <- function(model, variables) {
  variables <- as.character(variables)
  missing_variables <- setdiff(variables, posterior::variables(model))
  if (length(missing_variables) > 0L) {
    stop(
      "Required posterior variables are missing: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  draws <- posterior::as_draws_df(model, variable = variables)
  if (
    nrow(draws) < 1L ||
      !all(variables %in% names(draws)) ||
      anyNA(draws[, variables, drop = FALSE]) ||
      any(!is.finite(as.matrix(draws[, variables, drop = FALSE])))
  ) {
    stop("The requested posterior draws are invalid.", call. = FALSE)
  }

  draws
}

primary_indexed_posterior_variables <- function(model, prefix) {
  opening <- paste0(prefix, "[")
  candidates <- posterior::variables(model)
  candidates <- candidates[startsWith(candidates, opening)]
  indices <- suppressWarnings(as.integer(sub(
    "\\]$",
    "",
    substring(candidates, nchar(opening) + 1L)
  )))
  valid <- !is.na(indices)
  candidates <- candidates[valid]
  indices <- indices[valid]

  if (
    length(candidates) < 1L || anyDuplicated(indices) ||
      !identical(sort(indices), seq_len(length(indices)))
  ) {
    stop(
      "The indexed posterior variables are incomplete for ",
      prefix,
      ".",
      call. = FALSE
    )
  }

  candidates[order(indices)]
}

# Recover the study-specific basis from the fit ----

# Match the factor-smooth design and indexed posterior coefficients to study
# levels so integration uses the same basis as the fitted model.
map_primary_factor_smooth_design <- function(design, study_factor) {
  design <- as.matrix(design)
  study_factor <- droplevels(as.factor(study_factor))
  study_index <- as.integer(study_factor)
  n_studies <- nlevels(study_factor)
  tolerance <- sqrt(.Machine$double.eps)

  if (
    nrow(design) != length(study_index) || ncol(design) < n_studies ||
      anyNA(design) || any(!is.finite(design))
  ) {
    stop("A factor-smooth design matrix is invalid.", call. = FALSE)
  }

  column_study <- vapply(
    seq_len(ncol(design)),
    function(column) {
      nonzero_rows <- which(abs(design[, column]) > tolerance)
      represented_studies <- unique(study_index[nonzero_rows])
      if (length(represented_studies) != 1L) {
        return(NA_integer_)
      }
      represented_studies[[1L]]
    },
    integer(1)
  )

  columns_by_study <- split(seq_len(ncol(design)), column_study)
  component_counts <- lengths(columns_by_study)
  if (
    anyNA(column_study) || length(columns_by_study) != n_studies ||
      !identical(sort(as.integer(names(columns_by_study))), seq_len(n_studies)) ||
      length(unique(component_counts)) != 1L
  ) {
    stop(
      "The factor-smooth columns do not map evenly to study levels.",
      call. = FALSE
    )
  }

  n_components <- component_counts[[1L]]
  component_columns <- matrix(
    NA_integer_,
    nrow = n_studies,
    ncol = n_components
  )

  for (study in seq_len(n_studies)) {
    component_columns[study, ] <- columns_by_study[[as.character(study)]]
  }

  row <- seq_len(nrow(design))
  list(
    study_index = study_index,
    component_columns = component_columns,
    component_values = lapply(
      seq_len(n_components),
      function(component) {
        columns <- component_columns[, component]
        design[cbind(row, columns[study_index])]
      }
    )
  )
}





find_primary_factor_smooth_design <- function(
    standata,
    coefficient_variables,
    study_factor
) {
  design_names <- grep(
    "^Zs_[0-9]+_[0-9]+$",
    names(standata),
    value = TRUE
  )

  if (length(design_names) < 1L) {
    stop("No smooth design matrices were found.", call. = FALSE)
  }
  smooth_group <- as.integer(sub(
    "^Zs_([0-9]+)_[0-9]+$",
    "\\1",
    design_names
  ))
  smooth_component <- as.integer(sub(
    "^Zs_[0-9]+_([0-9]+)$",
    "\\1",
    design_names
  ))
  grouped_names <- split(design_names, smooth_group)
  grouped_components <- split(smooth_component, smooth_group)
  expected_columns <- lengths(coefficient_variables)

  candidates <- lapply(names(grouped_names), function(group) {
    names_in_group <- grouped_names[[group]]
    components_in_group <- grouped_components[[group]]
    names_in_group <- names_in_group[order(components_in_group)]
    if (
      length(names_in_group) != length(coefficient_variables) ||
        !identical(
          as.integer(vapply(
            names_in_group,
            function(x) ncol(standata[[x]]),
            integer(1)
          )),
          as.integer(expected_columns)
        )
    ) {
      return(NULL)
    }
    mappings <- tryCatch(
      lapply(
        names_in_group,
        function(x) map_primary_factor_smooth_design(
          standata[[x]],
          study_factor
        )
      ),
      error = function(e) NULL
    )
    if (
      is.null(mappings) ||
        !all(vapply(
          mappings,
          function(x) identical(x$study_index, mappings[[1L]]$study_index),
          logical(1)
        ))
    ) {
      return(NULL)
    }
    list(design_names = names_in_group, mappings = mappings)
  })

  candidates <- Filter(Negate(is.null), candidates)
  if (length(candidates) != 1L) {
    stop(
      "The study factor-smooth design could not be identified uniquely.",
      call. = FALSE
    )
  }

  candidates[[1L]]
}

get_primary_integrated_loo_components <- function(fitted_candidate) {
  model_index <- as.integer(fitted_candidate$model_index)
  if (!model_index %in% seq_len(26L)) {
    stop("The primary candidate model index is invalid.", call. = FALSE)
  }

  model <- fitted_candidate$model
  manifest <- get_primary_loo_observation_manifest(fitted_candidate)
  model_family <- stats::family(model)
  if (
    !identical(model_family$family, "gaussian") ||
      !identical(model_family$link, "identity")
  ) {
    stop(
      "Integrated primary LOO requires a Gaussian identity model.",
      call. = FALSE
    )
  }

  standata <- brms::standata(model)
  model_variables <- posterior::variables(model)
  common_variables <- c(
    residual = "sigma",
    arm = "sd_arm_no__Intercept",
    outcome = "sd_outcome_no__Intercept"
  )
  factor_smooth_indices <- c(11L, 14L, 23L, 26L)
  linear_slope_indices <- c(9L, 12L, 21L, 24L)
  spline_slope_indices <- c(10L, 13L, 22L, 25L)

  if (model_index %in% factor_smooth_indices) {
    scale_variables <- model_variables[
      startsWith(model_variables, "sds_") &
        grepl("study_factor_[0-9]+$", model_variables)
    ]
    scale_indices <- suppressWarnings(as.integer(sub(
      ".*_([0-9]+)$",
      "\\1",
      scale_variables
    )))
    scale_variables <- scale_variables[order(scale_indices)]
    scale_indices <- sort(scale_indices)
    if (
      length(scale_variables) < 1L || anyNA(scale_indices) ||
        !identical(scale_indices, seq_along(scale_indices))
    ) {
      stop(
        "The study factor-smooth scale variables are incomplete.",
        call. = FALSE
      )
    }
    coefficient_prefixes <- sub("^sds_", "s_", scale_variables)
    coefficient_variables <- lapply(
      coefficient_prefixes,
      function(x) primary_indexed_posterior_variables(model, x)
    )
    if (
      !"study_factor" %in% names(model$data) ||
        !is.factor(model$data$study_factor)
    ) {
      stop("The study factor used by the smooth is unavailable.", call. = FALSE)
    }
    smooth_design <- find_primary_factor_smooth_design(
      standata = standata,
      coefficient_variables = coefficient_variables,
      study_factor = model$data$study_factor
    )
    draws <- extract_primary_exact_draws(
      model,
      c(common_variables, scale_variables, unlist(coefficient_variables))
    )
    n_draws <- nrow(draws)
    n_observations <- nrow(manifest)
    factor_smooth_contribution <- matrix(
      0,
      nrow = n_draws,
      ncol = n_observations
    )
    study_basis <- list()
    for (block in seq_along(smooth_design$design_names)) {
      mapping <- smooth_design$mappings[[block]]
      coefficient_draws <- as.matrix(
        draws[, coefficient_variables[[block]], drop = FALSE]
      )
      for (component in seq_along(mapping$component_values)) {
        column_for_study <- mapping$component_columns[, component]
        values <- mapping$component_values[[component]]
        selected_columns <- column_for_study[mapping$study_index]
        for (observation in seq_len(n_observations)) {
          factor_smooth_contribution[, observation] <-
            factor_smooth_contribution[, observation] +
            values[[observation]] *
              coefficient_draws[, selected_columns[[observation]]]
        }
        study_basis[[length(study_basis) + 1L]] <- outer(
          draws[[scale_variables[[block]]]],
          values
        )
      }
    }
    conditional_linear_predictor <- brms::posterior_linpred(
      model,
      transform = FALSE,
      re_formula = NA
    )
    linear_predictor <-
      conditional_linear_predictor - factor_smooth_contribution
    structure <- "study_factor_smooth"
  } else {
    coefficient_labels <- if (model_index %in% linear_slope_indices) {
      c("Intercept", "protein_intake")
    } else if (model_index %in% spline_slope_indices) {
      c(
        "Intercept",
        "lsplineprotein_intake0.51",
        "lsplineprotein_intake0.52"
      )
    } else {
      "Intercept"
    }
    study_sd_variables <- paste0("sd_study_no__", coefficient_labels)
    design_names <- paste0("Z_3_", seq_along(coefficient_labels))
    if (!all(design_names %in% names(standata))) {
      stop("The study-level design columns are missing.", call. = FALSE)
    }
    study_design <- do.call(
      cbind,
      lapply(design_names, function(x) as.numeric(standata[[x]]))
    )
    if (
      nrow(study_design) != nrow(manifest) ||
        anyNA(study_design) || any(!is.finite(study_design))
    ) {
      stop("The study-level design matrix is invalid.", call. = FALSE)
    }
    if ("J_3" %in% names(standata)) {
      study_group <- as.integer(standata$J_3)
      group_counts <- vapply(
        split(study_group, manifest$study_no),
        function(x) length(unique(x)),
        integer(1)
      )
      study_counts <- vapply(
        split(manifest$study_no, study_group),
        function(x) length(unique(x)),
        integer(1)
      )
      if (any(group_counts != 1L) || any(study_counts != 1L)) {
        stop("Stan study groups do not match the LOO manifest.", call. = FALSE)
      }
    }

    correlation_variables <- character()
    correlation_pairs <- list()
    if (length(coefficient_labels) > 1L) {
      for (first in seq_len(length(coefficient_labels) - 1L)) {
        for (second in seq.int(first + 1L, length(coefficient_labels))) {
          forward <- paste0(
            "cor_study_no__",
            coefficient_labels[[first]],
            "__",
            coefficient_labels[[second]]
          )
          reverse <- paste0(
            "cor_study_no__",
            coefficient_labels[[second]],
            "__",
            coefficient_labels[[first]]
          )
          present <- intersect(c(forward, reverse), model_variables)
          if (length(present) > 1L) {
            stop("A study correlation is duplicated.", call. = FALSE)
          }
          if (length(present) == 1L) {
            correlation_variables <- c(correlation_variables, present)
            correlation_pairs[[length(correlation_pairs) + 1L]] <-
              c(first, second)
          }
        }
      }
    }
    expected_correlations <-
      length(coefficient_labels) * (length(coefficient_labels) - 1L) / 2L
    if (
      length(correlation_variables) != 0L &&
        length(correlation_variables) != expected_correlations
    ) {
      stop(
        "The study correlation matrix is only partially available.",
        call. = FALSE
      )
    }
    draws <- extract_primary_exact_draws(
      model,
      c(common_variables, study_sd_variables, correlation_variables)
    )
    n_draws <- nrow(draws)
    n_observations <- nrow(manifest)
    n_coefficients <- length(coefficient_labels)
    if (n_coefficients == 1L) {
      study_basis <- list(outer(draws[[study_sd_variables]], study_design[, 1L]))
      structure <- "study_intercept"
    } else {
      study_basis <- lapply(
        seq_len(n_coefficients),
        function(x) matrix(NA_real_, n_draws, n_observations)
      )
      for (draw in seq_len(n_draws)) {
        correlation <- diag(n_coefficients)
        if (length(correlation_variables) > 0L) {
          for (pair_index in seq_along(correlation_variables)) {
            pair <- correlation_pairs[[pair_index]]
            value <- draws[[correlation_variables[[pair_index]]]][[draw]]
            correlation[pair[[1L]], pair[[2L]]] <- value
            correlation[pair[[2L]], pair[[1L]]] <- value
          }
        }
        standard_deviations <- as.numeric(as.matrix(
          draws[draw, study_sd_variables, drop = FALSE]
        ))
        covariance <-
          diag(standard_deviations, n_coefficients) %*%
          correlation %*%
          diag(standard_deviations, n_coefficients)
        cholesky <- tryCatch(
          t(chol(covariance)),
          error = function(e) NULL
        )
        if (is.null(cholesky)) {
          stop(
            "A posterior study covariance matrix is not positive definite.",
            call. = FALSE
          )
        }
        transformed_design <- study_design %*% cholesky
        for (component in seq_len(n_coefficients)) {
          study_basis[[component]][draw, ] <-
            transformed_design[, component]
        }
      }
      structure_root <- if (model_index %in% linear_slope_indices) {
        "study_intercept_and_linear_slope"
      } else {
        "study_intercept_and_linear_spline_slopes"
      }
      structure <- paste0(
        if (length(correlation_variables) > 0L) {
          "correlated_"
        } else {
          "uncorrelated_"
        },
        structure_root
      )
    }
    linear_predictor <- brms::posterior_linpred(
      model,
      transform = FALSE,
      re_formula = NA
    )
  }

  n_draws <- nrow(draws)
  n_observations <- nrow(manifest)
  if (
    !identical(dim(linear_predictor), c(n_draws, n_observations)) ||
      length(study_basis) < 1L ||
      any(vapply(
        study_basis,
        function(x) {
          !identical(dim(x), c(n_draws, n_observations)) ||
            anyNA(x) || any(!is.finite(x))
        },
        logical(1)
      ))
  ) {
    stop(
      "The integrated posterior design does not align with its observations.",
      call. = FALSE
    )
  }

  list(
    structure = structure,
    linear_predictor = linear_predictor,
    residual_sd = draws[[common_variables[["residual"]]]],
    arm_intercept_sd = draws[[common_variables[["arm"]]]],
    outcome_intercept_sd = draws[[common_variables[["outcome"]]]],
    study_basis = study_basis,
    chain_id = as.integer(draws$.chain)
  )
}

# Integrated PSIS-LOO and diagnostics ----

# Rebuild the integrated likelihood for one accepted fit, compute PSIS-LOO,
# and retain the observation keys and Pareto-k values needed to select exact
# refits.
compute_primary_candidate_integrated_psis_loo <- function(
    fitted_candidate,
    initial_loo_result,
    k_threshold = 0.7,
    cores = 1L
) {
  candidate_identity <- c(
    outcome = as.character(fitted_candidate$outcome),
    model_index = as.character(fitted_candidate$model_index),
    model_name = as.character(fitted_candidate$model_name)
  )
  initial_identity <- c(
    outcome = as.character(initial_loo_result$outcome),
    model_index = as.character(initial_loo_result$model_index),
    model_name = as.character(initial_loo_result$model_name)
  )

  if (
    !identical(candidate_identity, initial_identity) ||
      !inherits(initial_loo_result$loo, "loo") ||
      !is.numeric(k_threshold) || length(k_threshold) != 1L ||
      is.na(k_threshold) || k_threshold <= 0 ||
      !is.numeric(cores) || length(cores) != 1L ||
      is.na(cores) || cores < 1L
  ) {
    stop("The integrated primary LOO inputs are invalid.", call. = FALSE)
  }

  manifest <- get_primary_loo_observation_manifest(fitted_candidate)
  if (!isTRUE(all.equal(
    manifest,
    initial_loo_result$observation_manifest,
    tolerance = 0,
    check.attributes = FALSE
  ))) {
    stop("The integrated LOO observation manifest has changed.", call. = FALSE)
  }

  start_time <- proc.time()[["elapsed"]]
  components <- get_primary_integrated_loo_components(fitted_candidate)
  integrated_log_lik <- integrate_primary_study_basis_log_lik(
    linear_predictor = components$linear_predictor,
    observed = manifest$yi,
    known_variance = manifest$vi,
    study_id = manifest$study_no,
    arm_id = manifest$arm_no,
    outcome_id = manifest$outcome_no,
    residual_sd = components$residual_sd,
    arm_intercept_sd = components$arm_intercept_sd,
    outcome_intercept_sd = components$outcome_intercept_sd,
    study_basis = components$study_basis
  )

  relative_efficiency <- loo::relative_eff(
    exp(integrated_log_lik),
    chain_id = components$chain_id,
    cores = as.integer(cores)
  )
  integrated_loo <- loo::loo(
    integrated_log_lik,
    r_eff = relative_efficiency,
    save_psis = FALSE,
    cores = as.integer(cores)
  )
  compute_seconds <- proc.time()[["elapsed"]] - start_time

  if (
    !inherits(integrated_loo, "loo") ||
      is.null(integrated_loo$diagnostics$pareto_k) ||
      length(integrated_loo$diagnostics$pareto_k) != nrow(manifest)
  ) {
    stop("The integrated primary PSIS-LOO result is invalid.", call. = FALSE)
  }

  initial_loo <- initial_loo_result$loo
  pointwise_difference <-
    integrated_loo$pointwise[, "elpd_loo"] -
    initial_loo$pointwise[, "elpd_loo"]
  comparison <- dplyr::bind_rows(
    summarise_primary_loo_method(
      initial_loo,
      method = "initial_conditional_psis",
      k_threshold = k_threshold
    ) |>
      dplyr::mutate(
        compute_seconds = NA_real_,
        delta_elpd_vs_initial = 0,
        se_delta_elpd_vs_initial = 0
      ),
    summarise_primary_loo_method(
      integrated_loo,
      method = "integrated_observation_level_psis",
      k_threshold = k_threshold
    ) |>
      dplyr::mutate(
        compute_seconds = as.numeric(compute_seconds),
        delta_elpd_vs_initial = sum(pointwise_difference),
        se_delta_elpd_vs_initial = sqrt(
          length(pointwise_difference) * stats::var(pointwise_difference)
        )
      )
  ) |>
    dplyr::mutate(
      outcome = candidate_identity[["outcome"]],
      model_index = as.integer(candidate_identity[["model_index"]]),
      model_name = candidate_identity[["model_name"]],
      study_structure = components$structure,
      .before = 1L
    )

  initial_k <- as.numeric(initial_loo$diagnostics$pareto_k)
  integrated_k <- as.numeric(integrated_loo$diagnostics$pareto_k)
  observation_diagnostics <- dplyr::mutate(
    manifest,
    model_index = as.integer(candidate_identity[["model_index"]]),
    model_name = candidate_identity[["model_name"]],
    study_structure = components$structure,
    initial_pareto_k = initial_k,
    integrated_pareto_k = integrated_k,
    integrated_resolved_initial_warning =
      initial_pareto_k > k_threshold & integrated_pareto_k <= k_threshold,
    integrated_new_warning =
      initial_pareto_k <= k_threshold & integrated_pareto_k > k_threshold,
    .after = outcome
  )

  list(
    outcome = candidate_identity[["outcome"]],
    model_index = as.integer(candidate_identity[["model_index"]]),
    model_name = candidate_identity[["model_name"]],
    study_structure = components$structure,
    observation_manifest = manifest,
    loo = integrated_loo,
    remedy = "integrated_observation_level_psis",
    k_threshold = k_threshold,
    compute_seconds = as.numeric(compute_seconds),
    comparison = comparison,
    observation_diagnostics = observation_diagnostics
  )
}

combine_primary_integrated_psis_loo_results <- function(...) {
  outcome_sets <- list(...)
  expected_outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (length(outcome_sets) != length(expected_outcomes)) {
    stop("Expected four integrated primary LOO outcome sets.", call. = FALSE)
  }

  validated_sets <- purrr::map2(
    outcome_sets,
    expected_outcomes,
    function(results, expected_outcome) {
      if (
        length(results) != 26L ||
          !setequal(purrr::map_int(results, "model_index"), seq_len(26L)) ||
          anyDuplicated(purrr::map_int(results, "model_index")) ||
          !all(purrr::map_lgl(results, ~ inherits(.x$loo, "loo"))) ||
          !all(purrr::map_chr(results, "outcome") == expected_outcome) ||
          !all(purrr::map_chr(results, "remedy") ==
            "integrated_observation_level_psis")
      ) {
        stop(
          "The integrated ", expected_outcome, " LOO set is invalid.",
          call. = FALSE
        )
      }
      results[order(purrr::map_int(results, "model_index"))]
    }
  )

  names(validated_sets) <- expected_outcomes

  method_comparison <- purrr::map_dfr(
    unlist(validated_sets, recursive = FALSE),
    "comparison"
  )
  observation_diagnostics <- purrr::map_dfr(
    unlist(validated_sets, recursive = FALSE),
    "observation_diagnostics"
  )
  model_summary <- dplyr::filter(
    method_comparison,
    method == "integrated_observation_level_psis"
  )
  outcome_summary <- model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      n_models = dplyr::n(),
      n_models_with_pareto_k_warning = sum(n_pareto_k_gt_threshold > 0L),
      total_pareto_k_warnings = sum(n_pareto_k_gt_threshold),
      maximum_pareto_k = max(maximum_pareto_k),
      total_compute_seconds = sum(compute_seconds),
      stacking_ready = all(n_pareto_k_gt_threshold == 0L),
      .groups = "drop"
    )

  list(
    results = validated_sets,
    model_summary = model_summary,
    outcome_summary = outcome_summary,
    method_comparison = method_comparison,
    observation_diagnostics = observation_diagnostics,
    stacking_ready = all(outcome_summary$stacking_ready)
  )
}

write_primary_integrated_psis_loo_audit <- function(
    audit,
    output_dir = file.path(
      "_targets", "user", "primary_model_loo_integrated"
    )
) {
  if (
    !is.list(audit$results) || length(audit$results) != 4L ||
      !is.data.frame(audit$model_summary) ||
      nrow(audit$model_summary) != 104L ||
      !is.data.frame(audit$outcome_summary) ||
      nrow(audit$outcome_summary) != 4L ||
      !is.data.frame(audit$method_comparison) ||
      nrow(audit$method_comparison) != 208L ||
      !is.data.frame(audit$observation_diagnostics)
  ) {
    stop("The integrated primary LOO audit is invalid.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    method_comparison = file.path(output_dir, "method_comparison.csv"),
    observations = file.path(output_dir, "observation_diagnostics.csv")
  )
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]], na = "")
  readr::write_csv(audit$model_summary, paths[["model_summary"]], na = "")
  readr::write_csv(
    audit$method_comparison,
    paths[["method_comparison"]],
    na = ""
  )
  readr::write_csv(
    audit$observation_diagnostics,
    paths[["observations"]],
    na = ""
  )
  unname(paths)
}

# Exact holdouts and replacements ----

# Identify the individual observations requiring an exact refit, then replace
# those pointwise scores in the integrated LOO results.
select_primary_integrated_exact_loo_jobs <- function(
    integrated_audit,
    fitted_whole_body,
    fitted_hypertrophy
) {
  expected_cases <- tibble::tribble(
    ~case_id, ~outcome, ~model_index, ~observation,
    1L, "whole_body", 7L, 3L,
    2L, "whole_body", 8L, 3L,
    3L, "hypertrophy", 20L, 252L
  )
  warning_cases <- integrated_audit$observation_diagnostics |>
    dplyr::filter(integrated_pareto_k > 0.7) |>
    dplyr::distinct(outcome, model_index, observation) |>
    dplyr::arrange(outcome, model_index, observation)
  expected_warning_cases <- expected_cases |>
    dplyr::select(outcome, model_index, observation) |>
    dplyr::arrange(outcome, model_index, observation)

  if (
    !isTRUE(all.equal(
      warning_cases,
      expected_warning_cases,
      tolerance = 0,
      check.attributes = FALSE
    )) ||
      length(fitted_whole_body) != 26L ||
      length(fitted_hypertrophy) != 26L
  ) {
    stop(
      "The exact-LOO case registry no longer matches the integrated audit.",
      call. = FALSE
    )
  }

  fitted_sets <- list(
    whole_body = fitted_whole_body,
    hypertrophy = fitted_hypertrophy
  )
  purrr::pmap(
    expected_cases,
    function(case_id, outcome, model_index, observation) {
      fitted_candidates <- fitted_sets[[outcome]]
      fitted_indices <- purrr::map_int(fitted_candidates, "model_index")
      fitted_candidate <- fitted_candidates[[match(model_index, fitted_indices)]]
      integrated_results <- integrated_audit$results[[outcome]]
      integrated_indices <- purrr::map_int(integrated_results, "model_index")
      integrated_result <-
        integrated_results[[match(model_index, integrated_indices)]]
      manifest <- integrated_result$observation_manifest
      if (
        is.null(fitted_candidate) || is.null(integrated_result) ||
          !identical(fitted_candidate$outcome, outcome) ||
          !identical(as.integer(fitted_candidate$model_index), model_index) ||
          !identical(integrated_result$study_structure, "study_intercept") ||
          observation > nrow(manifest) ||
          integrated_result$loo$diagnostics$pareto_k[[observation]] <= 0.7
      ) {
        stop("An exact-LOO job failed identity validation.", call. = FALSE)
      }
      list(
        case_id = case_id,
        outcome = outcome,
        model_index = model_index,
        model_name = fitted_candidate$model_name,
        observation = observation,
        fitted_candidate = fitted_candidate,
        integrated_result = integrated_result
      )
    }
  )
}

fit_primary_integrated_exact_loo_case <- function(
    job,
    chains = 4L,
    cores = 4L,
    warmup = 2000L,
    iter = 8000L,
    seed = 1988L,
    iter_sampling = 6000L
) {
  if (
    !inherits(job$fitted_candidate$model, "brmsfit") ||
      !identical(job$fitted_candidate$outcome, job$outcome) ||
      !identical(
        as.integer(job$fitted_candidate$model_index),
        as.integer(job$model_index)
      ) ||
      length(iter_sampling) != 1L || is.na(iter_sampling) ||
      !iter_sampling %in% c(6000L, 12000L) ||
      iter - warmup != iter_sampling || chains != 4L ||
      cores < 1L || job$observation < 1L
  ) {
    stop("The exact-LOO refit inputs are invalid.", call. = FALSE)
  }

  original_model <- job$fitted_candidate$model
  original_data <- original_model$data
  manifest <- job$integrated_result$observation_manifest
  if (
    nrow(original_data) != nrow(manifest) ||
      job$observation > nrow(original_data)
  ) {
    stop("The exact-LOO data do not match the manifest.", call. = FALSE)
  }

  training_data <- original_data[-job$observation, , drop = FALSE]
  if (nrow(training_data) != nrow(original_data) - 1L) {
    stop("The held-out observation was not removed uniquely.", call. = FALSE)
  }

  specification <- job$fitted_candidate$specification
  # Fits created during the staged remedy work can retain the brms prior-table
  # schema from the package version used at fit time. Normalise its columns
  # before update.brmsfit validates the unchanged priors for the refit.
  sample_prior <- attr(original_model$prior, "sample_prior")
  original_model$prior <- normalize_primary_brms_priors(
    original_model$prior
  )
  attr(original_model$prior, "sample_prior") <- sample_prior
  start_time <- proc.time()[["elapsed"]]
  refitted_model <- stats::update(
    original_model,
    newdata = training_data,
    recompile = FALSE,
    chains = as.integer(chains),
    cores = as.integer(cores),
    warmup = as.integer(warmup),
    iter = as.integer(iter),
    seed = as.integer(seed),
    control = list(
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )

  refit_seconds <- proc.time()[["elapsed"]] - start_time
  if (!inherits(refitted_model, "brmsfit")) {
    stop("The exact-LOO refit did not return a brmsfit.", call. = FALSE)
  }

  list(
    case_id = job$case_id,
    outcome = job$outcome,
    model_index = job$model_index,
    model_name = job$model_name,
    observation = job$observation,
    original_data = original_data,
    integrated_result = job$integrated_result,
    specification = specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(refit_seconds)
  )
}

summarise_primary_exact_loo_refit_diagnostics <- function(refit) {
  model <- refit$refitted_model
  nuts <- tibble::as_tibble(brms::nuts_params(model))
  divergences <- sum(
    nuts$Value[nuts$Parameter == "divergent__"],
    na.rm = TRUE
  )
  treedepth_hits <- sum(
    nuts$Value[nuts$Parameter == "treedepth__"] >=
      refit$specification$max_treedepth,
    na.rm = TRUE
  )
  energy <- nuts |>
    dplyr::filter(Parameter == "energy__") |>
    dplyr::group_by(Chain) |>
    dplyr::summarise(
      ebfmi = mean(diff(Value)^2) / stats::var(Value),
      .groups = "drop"
    )
  monitored_variables <- posterior::variables(model)
  monitored_variables <- monitored_variables[
    !grepl("^(lp__|lprior$|log_lik\\[)", monitored_variables)
  ]
  parameter_summary <- posterior::as_draws_array(
    model,
    variable = monitored_variables
  ) |>
    posterior::summarise_draws("rhat", "ess_bulk", "ess_tail")

  maximum_rhat <- max(parameter_summary$rhat, na.rm = TRUE)
  minimum_bulk_ess <- min(parameter_summary$ess_bulk, na.rm = TRUE)
  minimum_tail_ess <- min(parameter_summary$ess_tail, na.rm = TRUE)
  minimum_ebfmi <- min(energy$ebfmi, na.rm = TRUE)
  tibble::tibble(
    outcome = refit$outcome,
    model_index = refit$model_index,
    model_name = refit$model_name,
    observation = refit$observation,
    divergent_transitions = divergences,
    maximum_treedepth_transitions = treedepth_hits,
    minimum_ebfmi = minimum_ebfmi,
    maximum_rhat = maximum_rhat,
    minimum_bulk_ess = minimum_bulk_ess,
    minimum_tail_ess = minimum_tail_ess,
    geometry_pass = divergences == 0L && treedepth_hits == 0L &&
      minimum_ebfmi >= 0.30,
    convergence_pass = maximum_rhat <= 1.01 &&
      minimum_bulk_ess >= 400 && minimum_tail_ess >= 400,
    diagnostic_pass = geometry_pass && convergence_pass,
    refit_seconds = refit$refit_seconds
  )
}

update_primary_exact_loo_estimates <- function(loo_result) {
  pointwise <- loo_result$pointwise
  required_columns <- c("elpd_loo", "p_loo", "looic")
  if (!all(required_columns %in% colnames(pointwise))) {
    stop("The LOO pointwise matrix is incomplete.", call. = FALSE)
  }

  estimates <- loo_result$estimates
  update_row <- function(statistic, values) {
    if (statistic %in% rownames(estimates)) {
      estimates[statistic, "Estimate"] <<- sum(values)
      estimates[statistic, "SE"] <<-
        sqrt(length(values) * stats::var(values))
    }
  }
  update_row("elpd_loo", pointwise[, "elpd_loo"])
  update_row("p_loo", pointwise[, "p_loo"])
  update_row("looic", pointwise[, "looic"])
  loo_result$estimates <- estimates
  loo_result
}

# Score the excluded observation using the refitted posterior. Keep the
# observation's original key so only its pointwise contribution is replaced.
compute_primary_integrated_exact_loo_case <- function(refit) {
  integrated_result <- refit$integrated_result
  manifest <- integrated_result$observation_manifest
  observation <- as.integer(refit$observation)
  model <- refit$refitted_model
  sd_variables <- c(
    residual = "sigma",
    study = "sd_study_no__Intercept",
    arm = "sd_arm_no__Intercept",
    outcome = "sd_outcome_no__Intercept"
  )

  if (
    !identical(integrated_result$study_structure, "study_intercept") ||
      observation < 1L || observation > nrow(manifest) ||
      nrow(refit$original_data) != nrow(manifest)
  ) {
    stop("The exact integrated scoring inputs are invalid.", call. = FALSE)
  }

  population_linear_predictor <- brms::posterior_linpred(
    model,
    newdata = refit$original_data,
    re_formula = NA,
    transform = FALSE,
    allow_new_levels = TRUE
  )
  sd_draws <- extract_primary_exact_draws(model, sd_variables)

  if (
    nrow(population_linear_predictor) != nrow(sd_draws) ||
      ncol(population_linear_predictor) != nrow(manifest)
  ) {
    stop("Exact-LOO posterior predictions do not align.", call. = FALSE)
  }

  integrated_log_lik <- integrate_primary_study_basis_log_lik(
    linear_predictor = population_linear_predictor,
    observed = manifest$yi,
    known_variance = manifest$vi,
    study_id = manifest$study_no,
    arm_id = manifest$arm_no,
    outcome_id = manifest$outcome_no,
    residual_sd = sd_draws[[sd_variables[["residual"]]]],
    arm_intercept_sd = sd_draws[[sd_variables[["arm"]]]],
    outcome_intercept_sd = sd_draws[[sd_variables[["outcome"]]]],
    study_basis = list(outer(
      sd_draws[[sd_variables[["study"]]]],
      rep(1, nrow(manifest))
    ))
  )

  held_out_log_lik <- integrated_log_lik[, observation]
  maximum_log_lik <- max(held_out_log_lik)
  scaled_likelihood <- exp(held_out_log_lik - maximum_log_lik)
  exact_elpd <- maximum_log_lik + log(mean(scaled_likelihood))
  exact_mcse <- tryCatch(
    posterior::mcse_mean(scaled_likelihood) / mean(scaled_likelihood),
    error = function(e) NA_real_
  )

  patched_loo <- integrated_result$loo
  old_elpd <- patched_loo$pointwise[observation, "elpd_loo"]
  old_p_loo <- patched_loo$pointwise[observation, "p_loo"]
  full_data_lpd <- old_elpd + old_p_loo
  patched_loo$pointwise[observation, "elpd_loo"] <- exact_elpd
  patched_loo$pointwise[observation, "p_loo"] <- full_data_lpd - exact_elpd
  patched_loo$pointwise[observation, "looic"] <- -2 * exact_elpd

  if ("mcse_elpd_loo" %in% colnames(patched_loo$pointwise)) {
    patched_loo$pointwise[observation, "mcse_elpd_loo"] <- exact_mcse
  }
  old_pareto_k <- patched_loo$diagnostics$pareto_k[[observation]]
  patched_loo$diagnostics$pareto_k[[observation]] <- NA_real_
  if (
    !is.null(patched_loo$diagnostics$n_eff) &&
      length(patched_loo$diagnostics$n_eff) >= observation
  ) {
    patched_loo$diagnostics$n_eff[[observation]] <- NA_real_
  }

  patched_loo <- update_primary_exact_loo_estimates(patched_loo)

  initial_loo <- integrated_result$comparison |>
    dplyr::filter(method == "initial_conditional_psis")
  # Recover the initial pointwise object from the stored initial comparison is
  # impossible, so retain its existing summary row and update only the final
  # integrated row below. Pointwise deltas are adjusted by the exact change.
  final_summary <- summarise_primary_loo_method(
    patched_loo,
    method = "integrated_observation_level_psis",
    k_threshold = integrated_result$k_threshold
  ) |>
    dplyr::mutate(
      compute_seconds =
        integrated_result$compute_seconds + refit$refit_seconds,
      delta_elpd_vs_initial =
        integrated_result$comparison$delta_elpd_vs_initial[
          integrated_result$comparison$method ==
            "integrated_observation_level_psis"
        ] + exact_elpd - old_elpd,
      se_delta_elpd_vs_initial = NA_real_,
      exact_loo_refits = 1L,
      exact_refit_seconds = refit$refit_seconds
    )

  comparison <- dplyr::bind_rows(
    dplyr::mutate(
      initial_loo,
      exact_loo_refits = 0L,
      exact_refit_seconds = NA_real_
    ),
    final_summary
  ) |>
    dplyr::mutate(
      outcome = integrated_result$outcome,
      model_index = integrated_result$model_index,
      model_name = integrated_result$model_name,
      study_structure = integrated_result$study_structure,
      .before = 1L
    )

  observation_diagnostics <- integrated_result$observation_diagnostics |>
    dplyr::mutate(
      exact_loo_refit = dplyr::row_number() == observation,
      exact_elpd_loo = dplyr::if_else(
        exact_loo_refit,
        exact_elpd,
        NA_real_
      ),
      exact_elpd_mcse = dplyr::if_else(
        exact_loo_refit,
        exact_mcse,
        NA_real_
      ),
      integrated_pareto_k = dplyr::if_else(
        exact_loo_refit,
        NA_real_,
        integrated_pareto_k
      ),
      integrated_resolved_initial_warning = dplyr::if_else(
        exact_loo_refit,
        initial_pareto_k > integrated_result$k_threshold,
        integrated_resolved_initial_warning
      ),
      integrated_new_warning = dplyr::if_else(
        exact_loo_refit,
        FALSE,
        integrated_new_warning
      )
    )

  diagnostics <- summarise_primary_exact_loo_refit_diagnostics(refit)
  exact_summary <- dplyr::mutate(
    diagnostics,
    old_psis_pareto_k = old_pareto_k,
    old_psis_elpd_loo = old_elpd,
    exact_elpd_loo = exact_elpd,
    exact_minus_psis_elpd = exact_elpd - old_elpd,
    exact_elpd_mcse = exact_mcse
  )

  list(
    outcome = integrated_result$outcome,
    model_index = integrated_result$model_index,
    model_name = integrated_result$model_name,
    study_structure = integrated_result$study_structure,
    observation_manifest = manifest,
    loo = patched_loo,
    remedy = "integrated_observation_level_psis",
    k_threshold = integrated_result$k_threshold,
    compute_seconds =
      integrated_result$compute_seconds + refit$refit_seconds,
    comparison = comparison,
    observation_diagnostics = observation_diagnostics,
    exact_refit_summary = exact_summary,
    exact_loo_refits = 1L
  )
}

select_primary_exact_loo_geometry_remedy_job <- function(exact_jobs) {
  matches <- purrr::keep(
    exact_jobs,
    ~ identical(.x$outcome, "hypertrophy") &&
      identical(as.integer(.x$model_index), 20L) &&
      identical(as.integer(.x$observation), 252L)
  )
  if (length(matches) != 1L) {
    stop(
      "Expected one hypertrophy model-20 exact-LOO remedy job.",
      call. = FALSE
    )
  }
  matches[[1L]]
}

fit_primary_integrated_exact_loo_geometry_remedy <- function(
    job,
    adapt_delta = 0.9999,
    max_treedepth = 15L,
    chains = 4L,
    cores = 4L,
    warmup = 2000L,
    iter = 8000L,
    seed = 1988L
) {
  original_adapt_delta <- job$fitted_candidate$specification$adapt_delta
  original_max_treedepth <-
    job$fitted_candidate$specification$max_treedepth

  if (
    !identical(job$outcome, "hypertrophy") ||
      !identical(as.integer(job$model_index), 20L) ||
      !identical(as.integer(job$observation), 252L) ||
      !is.numeric(adapt_delta) || length(adapt_delta) != 1L ||
      is.na(adapt_delta) || adapt_delta <= original_adapt_delta ||
      adapt_delta >= 1 ||
      !is.numeric(max_treedepth) || length(max_treedepth) != 1L ||
      is.na(max_treedepth) || max_treedepth < original_max_treedepth
  ) {
    stop(
      "The hypertrophy exact-LOO geometry remedy is invalid.",
      call. = FALSE
    )
  }

  job$fitted_candidate$specification$adapt_delta <- adapt_delta
  job$fitted_candidate$specification$max_treedepth <-
    as.integer(max_treedepth)
  remedied <- fit_primary_integrated_exact_loo_case(
    job = job,
    chains = chains,
    cores = cores,
    warmup = warmup,
    iter = iter,
    seed = seed
  )

  remedied$geometry_remedy <- list(
    reason = "one_divergence_in_initial_exact_loo_refit",
    original_adapt_delta = original_adapt_delta,
    remedied_adapt_delta = adapt_delta,
    original_max_treedepth = original_max_treedepth,
    remedied_max_treedepth = as.integer(max_treedepth)
  )
  remedied
}

replace_primary_exact_loo_geometry_remedy <- function(
    exact_results,
    remedied_result
) {
  identities <- tibble::tibble(
    position = seq_along(exact_results),
    outcome = unname(purrr::map_chr(exact_results, "outcome")),
    model_index = unname(purrr::map_int(exact_results, "model_index"))
  )
  replacement_position <- which(
    identities$outcome == "hypertrophy" &
      identities$model_index == 20L
  )

  if (
    length(exact_results) != 3L ||
      length(replacement_position) != 1L ||
      !identical(remedied_result$outcome, "hypertrophy") ||
      !identical(as.integer(remedied_result$model_index), 20L) ||
      !identical(
        as.integer(remedied_result$exact_refit_summary$observation[[1L]]),
        252L
      )
  ) {
    stop("The exact-LOO geometry replacement is invalid.", call. = FALSE)
  }

  exact_results[[replacement_position]] <- remedied_result
  exact_results
}

# Apply the exact pointwise scores to the matching models before recomputing
# LOO summaries for stacking.
apply_primary_integrated_exact_loo_corrections <- function(
    integrated_audit,
    exact_results
) {
  if (
    length(exact_results) != 3L ||
      !identical(
        unname(purrr::map_chr(exact_results, "outcome")),
        c("whole_body", "whole_body", "hypertrophy")
      ) ||
      !identical(
        unname(purrr::map_int(exact_results, "model_index")),
        c(7L, 8L, 20L)
      )
  ) {
    stop("The exact integrated LOO correction set is invalid.", call. = FALSE)
  }

  final_results <- integrated_audit$results
  for (exact_result in exact_results) {
    outcome <- exact_result$outcome
    indices <- purrr::map_int(final_results[[outcome]], "model_index")
    position <- match(exact_result$model_index, indices)
    final_results[[outcome]][[position]] <- exact_result
  }

  flattened <- unlist(final_results, recursive = FALSE)
  method_comparison <- purrr::map_dfr(flattened, "comparison")
  observation_diagnostics <- purrr::map_dfr(
    flattened,
    function(result) {
      diagnostics <- result$observation_diagnostics
      if (!"exact_loo_refit" %in% names(diagnostics)) {
        diagnostics <- dplyr::mutate(
          diagnostics,
          exact_loo_refit = FALSE,
          exact_elpd_loo = NA_real_,
          exact_elpd_mcse = NA_real_
        )
      }
      diagnostics
    }
  )

  model_summary <- dplyr::filter(
    method_comparison,
    method == "integrated_observation_level_psis"
  )
  exact_refit_summary <- purrr::map_dfr(
    exact_results,
    "exact_refit_summary"
  )
  outcome_summary <- model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      n_models = dplyr::n(),
      n_models_with_pareto_k_warning = sum(n_pareto_k_gt_threshold > 0L),
      total_pareto_k_warnings = sum(n_pareto_k_gt_threshold),
      maximum_pareto_k = max(maximum_pareto_k),
      exact_loo_refits = sum(exact_loo_refits, na.rm = TRUE),
      total_compute_seconds = sum(compute_seconds),
      stacking_ready = all(n_pareto_k_gt_threshold == 0L),
      .groups = "drop"
    )

  exact_diagnostics_by_outcome <- exact_refit_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      exact_refit_diagnostics_pass = all(diagnostic_pass),
      .groups = "drop"
    )
  outcome_summary <- outcome_summary |>
    dplyr::left_join(exact_diagnostics_by_outcome, by = "outcome") |>
    dplyr::mutate(
      exact_refit_diagnostics_pass = tidyr::replace_na(
        exact_refit_diagnostics_pass,
        TRUE
      ),
      stacking_ready = stacking_ready & exact_refit_diagnostics_pass
    )

  diagnostic_pass <- all(exact_refit_summary$diagnostic_pass)

  list(
    results = final_results,
    model_summary = model_summary,
    outcome_summary = outcome_summary,
    method_comparison = method_comparison,
    observation_diagnostics = observation_diagnostics,
    exact_refit_summary = exact_refit_summary,
    stacking_ready = all(outcome_summary$stacking_ready) && diagnostic_pass
  )
}

write_primary_final_loo_audit <- function(
    audit,
    output_dir = file.path("_targets", "user", "primary_model_loo_final")
) {
  if (
    !is.list(audit$results) || length(audit$results) != 4L ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$method_comparison) != 208L ||
      nrow(audit$exact_refit_summary) != 3L
  ) {
    stop("The final primary LOO audit is invalid.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    method_comparison = file.path(output_dir, "method_comparison.csv"),
    observations = file.path(output_dir, "observation_diagnostics.csv"),
    exact_refits = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]], na = "")
  readr::write_csv(audit$model_summary, paths[["model_summary"]], na = "")
  readr::write_csv(
    audit$method_comparison,
    paths[["method_comparison"]],
    na = ""
  )
  readr::write_csv(
    audit$observation_diagnostics,
    paths[["observations"]],
    na = ""
  )
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refits"]],
    na = ""
  )
  unname(paths)
}
