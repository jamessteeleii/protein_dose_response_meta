# Pilot infrastructure for supplementary sensitivity analyses.
#
# These functions are deliberately separate from the primary-analysis
# functions. The targets that call a sampler are opt-in (see
# sensitivity/README.md); the structural targets can be built without
# compiling Stan or drawing from a posterior.

sensitivity_model_features <- function(model_index) {
  if (
    length(model_index) != 1L || is.na(model_index) ||
      !model_index %in% seq_len(26L)
  ) {
    stop("`model_index` must be one of 1:26.", call. = FALSE)
  }

  position <- ((as.integer(model_index) - 1L) %% 14L) + 1L
  protein_family <- dplyr::case_when(
    model_index %in% c(1L, 5L) ~ "none",
    model_index %in% c(2L, 6L, 9L, 12L, 15L, 18L, 21L, 24L) ~
      "linear",
    model_index %in% c(3L, 7L, 10L, 13L, 16L, 19L, 22L, 25L) ~
      "piecewise",
    TRUE ~ "spline"
  )

  list(
    protein_family = protein_family,
    includes_protein = protein_family != "none",
    includes_moderators = model_index %in% c(
      5:8, 12:14, 18:20, 24:26
    ),
    includes_random_slope = model_index %in% c(9:14, 21:26),
    includes_rt_interaction = model_index >= 15L,
    position = position
  )
}

materialise_final_primary_specifications <- function(
    primary_specifications,
    first_remedy_specifications,
    second_remedy_specifications,
    third_remedy_specifications
) {
  add_stage <- function(x, stage) {
    dplyr::mutate(dplyr::bind_rows(x), sensitivity_remedy_stage = stage)
  }

  specification_table <- dplyr::bind_rows(
    add_stage(primary_specifications, 0L),
    add_stage(first_remedy_specifications, 1L),
    add_stage(second_remedy_specifications, 2L),
    add_stage(third_remedy_specifications, 3L)
  ) |>
    dplyr::arrange(outcome, model_index, sensitivity_remedy_stage) |>
    dplyr::group_by(outcome, model_index) |>
    dplyr::slice_tail(n = 1L) |>
    dplyr::ungroup()

  logical_fields <- c(
    "uncorrelated_study_effects", "regularised_variance_priors",
    "full_schedule_trial", "standardised_moderators",
    "half_unit_protein", "natural_spline", "targeted_smooth_prior",
    "newly_uncorrelated_study_effects", "requires_refit",
    "targeted_piecewise_slope_priors", "second_remedy",
    "targeted_function_tail_priors",
    "targeted_piecewise_interaction_priors", "third_remedy"
  )
  missing_fields <- setdiff(logical_fields, names(specification_table))
  for (field in missing_fields) specification_table[[field]] <- FALSE
  specification_table <- specification_table |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(logical_fields),
        ~ dplyr::coalesce(.x, FALSE)
      ),
      smooth_prior_scale = dplyr::coalesce(smooth_prior_scale, NA_real_)
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  pair <- paste(specification_table$outcome, specification_table$model_index)
  if (
    nrow(specification_table) != 104L || anyDuplicated(pair) ||
      !identical(
        split(specification_table$model_index, specification_table$outcome),
        stats::setNames(
          rep(list(seq_len(26L)), 4L),
          c("function", "hypertrophy", "strength", "whole_body")
        )[names(split(specification_table$model_index,
                      specification_table$outcome))]
      )
  ) {
    stop("The final sensitivity specification registry is incomplete.",
         call. = FALSE)
  }

  specification_list <- purrr::pmap(specification_table, list)
  names(specification_list) <- paste0(
    specification_table$outcome,
    "_model_",
    stringr::str_pad(specification_table$model_index, 2L, pad = "0")
  )
  specification_list
}

make_sensitivity_pilot_config <- function(
    outcome = "whole_body",
    model_indices = c(2L, 3L, 4L, 13L, 26L),
    residual_correlation_model_index = 6L,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 1000L,
    iter_sampling = 1000L,
    seed = 1988L
) {
  valid_outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (!outcome %in% valid_outcomes) {
    stop("Unknown pilot outcome.", call. = FALSE)
  }
  model_indices <- unique(as.integer(model_indices))
  if (
    length(model_indices) < 1L ||
      any(!model_indices %in% setdiff(seq_len(26L), c(1L, 5L)))
  ) {
    stop("Pilot model indices must identify protein-containing candidates.",
         call. = FALSE)
  }
  if (!residual_correlation_model_index %in% seq_len(26L)) {
    stop("Invalid residual-correlation model index.", call. = FALSE)
  }

  list(
    outcome = outcome,
    model_indices = model_indices,
    null_model_indices = c(1L, 5L),
    residual_correlation_model_index =
      as.integer(residual_correlation_model_index),
    chains = as.integer(chains),
    parallel_chains = as.integer(parallel_chains),
    iter_warmup = as.integer(iter_warmup),
    iter_sampling = as.integer(iter_sampling),
    seed = as.integer(seed),
    adapt_delta_floor = 0.99,
    max_treedepth_floor = 12L,
    prediction_draws = 1000L
  )
}

make_mundlak_remedy_pilot_config <- function(
    base_config,
    remedy_model_indices = c(1L, 2L, 3L, 4L, 13L, 26L),
    regularised_model_indices = 4L,
    iter_warmup = 1500L,
    iter_sampling = 2000L
) {
  remedy_model_indices <- unique(as.integer(remedy_model_indices))
  regularised_model_indices <- unique(as.integer(regularised_model_indices))
  comparison_indices <- sort(unique(c(
    base_config$null_model_indices,
    base_config$model_indices
  )))
  if (
    length(remedy_model_indices) < 1L ||
      !all(remedy_model_indices %in% comparison_indices) ||
      !all(regularised_model_indices %in% remedy_model_indices)
  ) {
    stop("Invalid Mundlak remedy model indices.", call. = FALSE)
  }
  config <- base_config
  config$remedy_model_indices <- remedy_model_indices
  config$regularised_model_indices <- regularised_model_indices
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + 100L)
  config
}

make_mundlak_full_expansion_config <- function(
    base_config,
    iter_warmup = 1500L,
    iter_sampling = 2000L,
    regularised_model_indices = 4L,
    seed_offset = 2000L
) {
  config <- base_config
  config$outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  config$model_indices <- seq_len(26L)
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$regularised_model_indices <- unique(as.integer(
    regularised_model_indices
  ))
  config$seed <- as.integer(base_config$seed + seed_offset)
  if (
    config$iter_warmup < 1L || config$iter_sampling < 1L ||
      any(!config$regularised_model_indices %in% config$model_indices)
  ) {
    stop("Invalid full Mundlak expansion configuration.", call. = FALSE)
  }
  config
}

make_within_remedy_pilot_config <- function(
    base_config,
    remedy_model_indices = 1L,
    iter_warmup = 1500L,
    iter_sampling = 2000L
) {
  remedy_model_indices <- unique(as.integer(remedy_model_indices))
  comparison_indices <- sort(unique(c(
    base_config$null_model_indices,
    base_config$model_indices
  )))
  if (
    length(remedy_model_indices) < 1L ||
      !all(remedy_model_indices %in% comparison_indices)
  ) {
    stop("Invalid within-contrast remedy model indices.", call. = FALSE)
  }
  config <- base_config
  config$remedy_model_indices <- remedy_model_indices
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + 300L)
  config
}

make_within_full_expansion_config <- function(
    base_config,
    iter_warmup = 1500L,
    iter_sampling = 2000L,
    seed_offset = 6000L
) {
  config <- base_config
  config$outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  config$model_indices <- seq_len(26L)
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains < 1L ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("Invalid full within-contrast expansion configuration.",
         call. = FALSE)
  }
  config
}

make_within_full_remedy_config <- function(
    full_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    adapt_delta = 0.999,
    max_treedepth = 15L,
    seed_offset = 1000L
) {
  config <- full_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$adapt_delta_floor <- as.numeric(adapt_delta)
  config$max_treedepth_floor <- as.integer(max_treedepth)
  config$seed <- as.integer(full_config$seed + seed_offset)
  if (
    config$iter_warmup < full_config$iter_warmup ||
      config$iter_sampling <= full_config$iter_sampling ||
      !is.finite(config$adapt_delta_floor) ||
      config$adapt_delta_floor < full_config$adapt_delta_floor ||
      config$adapt_delta_floor >= 1 ||
      config$max_treedepth_floor < full_config$max_treedepth_floor
  ) {
    stop("Invalid full within-contrast remedy configuration.",
         call. = FALSE)
  }
  config
}

make_exposure_full_expansion_config <- function(
    base_config,
    iter_warmup = 1500L,
    iter_sampling = 2000L,
    seed_offset = 12000L
) {
  config <- base_config
  config$outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  config$model_indices <- seq_len(26L)
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains < 1L ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("Invalid full exposure expansion configuration.",
         call. = FALSE)
  }
  config
}

make_exposure_full_remedy_config <- function(
    full_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    adapt_delta = 0.999,
    max_treedepth = 15L,
    seed_offset = 2000L
) {
  config <- full_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$adapt_delta_floor <- as.numeric(adapt_delta)
  config$max_treedepth_floor <- as.integer(max_treedepth)
  config$seed <- as.integer(full_config$seed + seed_offset)
  if (
    config$iter_warmup < full_config$iter_warmup ||
      config$iter_sampling <= full_config$iter_sampling ||
      !is.finite(config$adapt_delta_floor) ||
      config$adapt_delta_floor < full_config$adapt_delta_floor ||
      config$adapt_delta_floor >= 1 ||
      config$max_treedepth_floor < full_config$max_treedepth_floor
  ) {
    stop("Invalid full exposure remedy configuration.", call. = FALSE)
  }
  config
}

make_exposure_remedy_pilot_config <- function(
    base_config,
    achieved_model_indices = c(1L, 2L, 3L, 26L),
    assignment_model_indices = c(2L, 26L),
    shared_null_model_index = 1L,
    iter_warmup = 1500L,
    iter_sampling = 2000L
) {
  comparison_indices <- sort(unique(c(
    base_config$null_model_indices,
    base_config$model_indices
  )))
  achieved_model_indices <- unique(as.integer(achieved_model_indices))
  assignment_model_indices <- unique(as.integer(assignment_model_indices))
  shared_null_model_index <- as.integer(shared_null_model_index)
  if (
    length(shared_null_model_index) != 1L ||
      !shared_null_model_index %in% base_config$null_model_indices ||
      !shared_null_model_index %in% achieved_model_indices ||
      length(achieved_model_indices) < 1L ||
      length(assignment_model_indices) < 1L ||
      !all(c(achieved_model_indices, assignment_model_indices) %in%
        comparison_indices) ||
      shared_null_model_index %in% assignment_model_indices
  ) {
    stop("Invalid exposure remedy model indices.", call. = FALSE)
  }
  config <- base_config
  config$achieved_model_indices <- achieved_model_indices
  config$assignment_model_indices <- assignment_model_indices
  config$shared_null_model_index <- shared_null_model_index
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + 400L)
  config
}

make_exposure_second_remedy_pilot_config <- function(
    base_config,
    model_indices = 3L,
    iter_warmup = 2000L,
    iter_sampling = 4000L
) {
  model_indices <- unique(as.integer(model_indices))
  comparison_indices <- sort(unique(c(
    base_config$null_model_indices,
    base_config$model_indices
  )))
  if (
    length(model_indices) < 1L ||
      !all(model_indices %in% comparison_indices)
  ) {
    stop("Invalid second exposure remedy model indices.", call. = FALSE)
  }
  config <- base_config
  config$model_indices <- model_indices
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + 500L)
  config
}

make_eiv_smoke_remedy_pilot_config <- function(
    base_config,
    iter_warmup = 1500L,
    iter_sampling = 2000L
) {
  config <- base_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + 700L)
  config
}

make_eiv_smoke_second_remedy_pilot_config <- function(
    base_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L
) {
  config <- base_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + 800L)
  config
}

make_eiv_ensemble_remedy_pilot_config <- function(
    base_config,
    remedy_model_indices,
    iter_warmup,
    iter_sampling,
    seed_offset
) {
  eligible <- sort(unique(c(
    base_config$null_model_indices,
    base_config$model_indices
  )))
  remedy_model_indices <- sort(unique(as.integer(remedy_model_indices)))
  if (
    length(remedy_model_indices) < 1L ||
      !all(remedy_model_indices %in% eligible)
  ) {
    stop("Invalid EIV remedy model indices.", call. = FALSE)
  }
  config <- base_config
  config$remedy_model_indices <- remedy_model_indices
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + seed_offset)
  config
}

parse_sensitivity_pilot_selection <- function(selection) {
  families <- c(
    "mundlak", "within", "exposure", "rescor", "eiv", "strategy",
    "strategy_within"
  )
  if (is.null(selection)) selection <- ""
  selection <- tolower(trimws(selection))
  selected <- unique(trimws(unlist(strsplit(selection, ",", fixed = TRUE))))
  selected <- selected[nzchar(selected)]
  unknown <- setdiff(selected, c(families, "all", "none"))
  if (length(unknown) > 0L) {
    stop(
      "Unknown sensitivity pilot family: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if ("all" %in% selected) selected <- families
  if ("none" %in% selected) selected <- character()
  stats::setNames(families %in% selected, families)
}

new_sensitivity_pilot_skip <- function(family) {
  structure(
    list(
      status = "not_run",
      family = family,
      reason = paste0(
        "Set PROTEIN_SENSITIVITY_PILOTS=", family,
        " (or all) before explicitly making this pilot target."
      )
    ),
    class = c("sensitivity_pilot_skip", "list")
  )
}

is_sensitivity_pilot_skip <- function(x) {
  inherits(x, "sensitivity_pilot_skip")
}

select_sensitivity_specifications <- function(
    specifications,
    outcome,
    model_indices
) {
  selected <- purrr::keep(
    specifications,
    ~ identical(.x$outcome, outcome) && .x$model_index %in% model_indices
  )
  selected <- selected[order(purrr::map_int(selected, "model_index"))]
  if (!identical(unname(purrr::map_int(selected, "model_index")),
                 sort(unique(as.integer(model_indices))))) {
    stop("The requested sensitivity specifications were not all found.",
         call. = FALSE)
  }
  selected
}

infer_sensitivity_basis_family <- function(specification) {
  features <- sensitivity_model_features(specification$model_index)
  if (!features$includes_protein) return("none")
  if (features$protein_family != "spline") return(features$protein_family)
  if (isTRUE(specification$natural_spline)) "natural_spline" else "thin_plate"
}

make_sensitivity_basis_spec <- function(data, specification, knot = 0.5) {
  required <- c("study_no", "arm_no", "protein_intake")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Basis data are missing: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  family <- infer_sensitivity_basis_family(specification)
  if (family == "none") {
    stop("A protein basis cannot be made for an intercept-only candidate.",
         call. = FALSE)
  }

  arm_data <- data |>
    dplyr::distinct(study_no, arm_no, protein_intake) |>
    dplyr::filter(is.finite(protein_intake))
  if (nrow(arm_data) < 5L) {
    stop("Too few distinct study-arms to define a protein basis.",
         call. = FALSE)
  }
  x_scale <- if (isTRUE(specification$half_unit_protein)) 0.5 else 1
  x_model <- arm_data$protein_intake / x_scale

  output <- list(
    family = family,
    knot = knot,
    x_scale = x_scale,
    model_index = specification$model_index,
    outcome = specification$outcome
  )
  if (family == "natural_spline") {
    raw <- splines::ns(x_model, df = 3)
    output$knots <- attr(raw, "knots")
    output$boundary_knots <- attr(raw, "Boundary.knots")
    output$penalized_columns <- integer()
  } else if (family == "thin_plate") {
    smooth <- mgcv::smoothCon(
      mgcv::s(x_model, bs = "tp", k = 4),
      data = data.frame(x_model = x_model),
      absorb.cons = TRUE
    )[[1L]]
    penalty <- Reduce(`+`, smooth$S)
    penalty_eigen <- eigen(penalty, symmetric = TRUE)
    tolerance <- max(penalty_eigen$values) * sqrt(.Machine$double.eps)
    penalized <- penalty_eigen$values > tolerance
    basis_transform <- cbind(
      penalty_eigen$vectors[, !penalized, drop = FALSE],
      sweep(
        penalty_eigen$vectors[, penalized, drop = FALSE],
        2L,
        sqrt(penalty_eigen$values[penalized]),
        "/"
      )
    )
    raw <- smooth$X %*% basis_transform
    output$smooth <- smooth
    output$basis_transform <- basis_transform
    output$penalized_columns <- seq.int(
      sum(!penalized) + 1L,
      ncol(basis_transform)
    )
  } else {
    output$penalized_columns <- integer()
  }

  output$basis_dimension <- switch(
    family,
    linear = 1L,
    piecewise = 2L,
    natural_spline = ncol(raw),
    thin_plate = ncol(raw)
  )
  class(output) <- c("sensitivity_basis_spec", class(output))
  output
}

predict_sensitivity_basis <- function(basis_spec, x) {
  if (!inherits(basis_spec, "sensitivity_basis_spec")) {
    stop("Invalid sensitivity basis specification.", call. = FALSE)
  }
  x <- as.numeric(x)
  x_model <- x / basis_spec$x_scale
  basis <- switch(
    basis_spec$family,
    linear = matrix(x_model, ncol = 1L),
    piecewise = cbind(
      pmin(x_model, basis_spec$knot / basis_spec$x_scale),
      pmax(x_model - basis_spec$knot / basis_spec$x_scale, 0)
    ),
    natural_spline = splines::ns(
      x_model,
      knots = basis_spec$knots,
      Boundary.knots = basis_spec$boundary_knots
    ),
    thin_plate = mgcv::PredictMat(
      basis_spec$smooth,
      data.frame(x_model = x_model)
    ) %*% basis_spec$basis_transform
  )
  basis <- as.matrix(basis)
  if (
    nrow(basis) != length(x) ||
      ncol(basis) != basis_spec$basis_dimension ||
      any(!is.finite(basis))
  ) {
    stop("Protein-basis prediction failed validation.", call. = FALSE)
  }
  colnames(basis) <- paste0("basis_", seq_len(ncol(basis)))
  basis
}

add_sensitivity_mundlak_basis <- function(data, basis_spec) {
  arm_data <- data |>
    dplyr::distinct(study_no, arm_no, protein_intake)
  arm_basis <- predict_sensitivity_basis(
    basis_spec,
    arm_data$protein_intake
  )
  basis_names <- paste0("protein_basis_", seq_len(ncol(arm_basis)))
  arm_basis <- tibble::as_tibble(arm_basis, .name_repair = "minimal")
  names(arm_basis) <- basis_names
  arm_data <- dplyr::bind_cols(arm_data, arm_basis)

  between <- arm_data |>
    dplyr::group_by(study_no) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(basis_names), mean),
      .groups = "drop"
    ) |>
    dplyr::rename_with(
      ~ sub("protein_basis_", "protein_between_", .x, fixed = TRUE),
      dplyr::all_of(basis_names)
    )
  arm_data <- arm_data |>
    dplyr::left_join(between, by = "study_no")

  within_names <- sub(
    "protein_basis_", "protein_within_", basis_names, fixed = TRUE
  )
  between_names <- sub(
    "protein_basis_", "protein_between_", basis_names, fixed = TRUE
  )
  for (index in seq_along(basis_names)) {
    arm_data[[within_names[[index]]]] <-
      arm_data[[basis_names[[index]]]] -
      arm_data[[between_names[[index]]]]
  }

  output <- data |>
    dplyr::left_join(
      dplyr::select(
        arm_data,
        study_no, arm_no,
        dplyr::all_of(c(within_names, between_names))
      ),
      by = c("study_no", "arm_no")
    )
  if (anyNA(output[, c(within_names, between_names), drop = FALSE])) {
    stop("Mundlak basis columns contain missing values.", call. = FALSE)
  }
  attr(output, "sensitivity_basis_spec") <- basis_spec
  output
}

make_sensitivity_basis_formula <- function(
    specification,
    within_names = NULL,
    between_names = NULL,
    contrast_names = NULL
) {
  features <- sensitivity_model_features(specification$model_index)
  if (!features$includes_protein) {
    stop("Sensitivity basis formulas require a protein candidate.",
         call. = FALSE)
  }
  mundlak <- !is.null(within_names) && !is.null(between_names)
  contrast <- !is.null(contrast_names)
  if (mundlak == contrast) {
    stop("Specify either Mundlak or contrast basis columns.", call. = FALSE)
  }

  protein_main <- if (mundlak) c(within_names, between_names) else contrast_names
  terms <- c(
    "0 + Intercept",
    "resistance_exercise_code",
    protein_main
  )
  if (features$includes_rt_interaction) {
    terms <- c(
      terms,
      paste0(protein_main, ":resistance_exercise_code")
    )
  }
  if (features$includes_moderators) {
    moderator_names <- if (isTRUE(specification$standardised_moderators)) {
      c("duration_outcome_z", "age_outcome_z")
    } else {
      c("duration_centre", "age_centre")
    }
    terms <- c(terms, moderator_names)
  }

  if (features$includes_random_slope) {
    slope_names <- if (mundlak) within_names else contrast_names
    separator <- if (isTRUE(specification$uncorrelated_study_effects)) {
      " || study_no)"
    } else {
      " | study_no)"
    }
    study_term <- paste0(
      "(1 + ", paste(slope_names, collapse = " + "), separator
    )
  } else {
    study_term <- "(1 | study_no)"
  }

  paste0(
    "yi|se(sqrt(vi)) ~ ",
    paste(c(terms, study_term, "(1 | arm_no)", "(1 | outcome_no)"),
          collapse = " + ")
  )
}

prepare_sensitivity_basis_priors <- function(
    priors,
    specification,
    coefficient_names,
    basis_spec
) {
  priors <- normalize_primary_brms_priors(priors)
  if (isTRUE(specification$standardised_moderators)) {
    priors <- add_primary_standardised_moderator_priors(priors)
  }
  if (isTRUE(specification$regularised_variance_priors)) {
    priors <- c(
      priors,
      brms::set_prior("student_t(3, 0, 1)", class = "sd")
    )
  }

  family <- infer_sensitivity_basis_family(specification)
  if (
    family == "piecewise" &&
      (isTRUE(specification$targeted_piecewise_slope_priors) ||
        isTRUE(specification$targeted_function_tail_priors))
  ) {
    for (coefficient in coefficient_names) {
      priors <- c(
        priors,
        brms::set_prior(
          "student_t(3, 0, 0.5)",
          class = "b",
          coef = coefficient
        )
      )
    }
  }
  if (family == "thin_plate") {
    scale <- if (
      isTRUE(specification$targeted_smooth_prior) &&
        is.finite(specification$smooth_prior_scale)
    ) specification$smooth_prior_scale else 0.5
    coefficient_basis_index <- as.integer(sub(".*_", "", coefficient_names))
    penalized_coefficients <- coefficient_names[
      coefficient_basis_index %in% basis_spec$penalized_columns
    ]
    for (coefficient in penalized_coefficients) {
      priors <- c(
        priors,
        brms::set_prior(
          paste0("normal(0, ", scale, ")"),
          class = "b",
          coef = coefficient
        )
      )
    }
  }
  priors
}

recenter_within_contrast_priors <- function(priors) {
  priors <- normalize_primary_brms_priors(priors)
  prior_table <- as.data.frame(priors)
  for (coefficient in c("Intercept", "resistance_exercise_code")) {
    row <- which(prior_table$class == "b" & prior_table$coef == coefficient)
    if (length(row) != 1L) {
      stop("Expected one prior for contrast coefficient: ", coefficient,
           call. = FALSE)
    }
    parsed <- extract_student_t_prior(prior_table$prior[[row]])
    prior_table$prior[[row]] <- paste0(
      "student_t(", parsed[["df"]], ", 0, ", parsed[["scale"]], ")"
    )
  }
  structure(prior_table, class = c("brmsprior", "data.frame"))
}

prepare_mundlak_pilot_input <- function(data, specification, priors) {
  features <- sensitivity_model_features(specification$model_index)
  if (!features$includes_protein) {
    output <- prepare_sensitivity_final_candidate(
      specification, priors, data
    )
    output$family <- "mundlak"
    output$basis_spec <- NULL
    output$within_names <- NULL
    output$between_names <- NULL
    return(output)
  }
  basis_spec <- make_sensitivity_basis_spec(data, specification)
  data <- add_sensitivity_mundlak_basis(data, basis_spec)
  if (isTRUE(specification$standardised_moderators)) {
    data <- add_primary_standardised_moderators(data)
  }
  within_names <- paste0(
    "protein_within_", seq_len(basis_spec$basis_dimension)
  )
  between_names <- paste0(
    "protein_between_", seq_len(basis_spec$basis_dimension)
  )
  formula <- make_sensitivity_basis_formula(
    specification,
    within_names = within_names,
    between_names = between_names
  )
  priors <- prepare_sensitivity_basis_priors(
    priors,
    specification,
    coefficient_names = c(within_names, between_names),
    basis_spec = basis_spec
  )

  list(
    family = "mundlak",
    specification = specification,
    formula = formula,
    priors = priors,
    data = data,
    basis_spec = basis_spec,
    within_names = within_names,
    between_names = between_names
  )
}

prepare_mundlak_remedy_pilot_input <- function(
    data,
    specification,
    priors,
    remedy_config
) {
  specification$sensitivity_pilot_remedy <- "longer_sampling"
  if (
    specification$model_index %in%
      remedy_config$regularised_model_indices
  ) {
    specification$regularised_variance_priors <- TRUE
    specification$sensitivity_pilot_remedy <-
      "longer_sampling_and_regularised_variance_priors"
  }
  prepare_mundlak_pilot_input(data, specification, priors)
}

select_pending_mundlak_full_specifications <- function(
    specifications,
    reused_fits
) {
  if (is_sensitivity_pilot_skip(reused_fits)) {
    return(new_sensitivity_pilot_skip("mundlak"))
  }
  specification_keys <- paste(
    purrr::map_chr(specifications, "outcome"),
    purrr::map_int(specifications, "model_index")
  )
  reused_keys <- paste(
    purrr::map_chr(reused_fits, "outcome"),
    purrr::map_int(reused_fits, "model_index")
  )
  expected_reused <- paste(
    "whole_body", c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
  )
  if (
    length(specifications) != 104L || anyDuplicated(specification_keys) ||
      !setequal(reused_keys, expected_reused) ||
      !all(reused_keys %in% specification_keys)
  ) {
    stop("The full Mundlak specification/reuse set is invalid.",
         call. = FALSE)
  }
  specifications[!specification_keys %in% reused_keys]
}

prepare_mundlak_full_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  outcome <- specification$outcome
  data <- data_by_outcome[[outcome]]
  priors <- priors_by_outcome[[outcome]]
  if (is.null(data) || is.null(priors) || !outcome %in% config$outcomes) {
    stop("A full Mundlak outcome input is unavailable.", call. = FALSE)
  }
  specification$sensitivity_mundlak_expansion <- TRUE
  specification$sensitivity_mundlak_sampling <- paste0(
    config$iter_warmup, "_warmup_", config$iter_sampling, "_sampling"
  )
  if (specification$model_index %in% config$regularised_model_indices) {
    specification$regularised_variance_priors <- TRUE
    specification$sensitivity_mundlak_pilot_informed_remedy <-
      "regularised_variance_priors"
  } else {
    specification$sensitivity_mundlak_pilot_informed_remedy <-
      "longer_sampling"
  }
  prepare_mundlak_pilot_input(data, specification, priors)
}

audit_mundlak_full_expansion <- function(
    pending_specifications,
    inputs,
    reused_fits,
    config
) {
  if (
    is_sensitivity_pilot_skip(pending_specifications) ||
      is_sensitivity_pilot_skip(inputs) ||
      is_sensitivity_pilot_skip(reused_fits)
  ) {
    return(new_sensitivity_pilot_skip("mundlak"))
  }
  pending_keys <- paste(
    purrr::map_chr(pending_specifications, "outcome"),
    purrr::map_int(pending_specifications, "model_index")
  )
  input_keys <- paste(
    purrr::map_chr(inputs, ~ .x$specification$outcome),
    purrr::map_int(inputs, ~ .x$specification$model_index)
  )
  input_outcomes <- purrr::map_chr(inputs, ~ .x$specification$outcome)
  expected_counts <- c(
    whole_body = 19L,
    hypertrophy = 26L,
    strength = 26L,
    `function` = 26L
  )
  observed_counts <- table(factor(
    input_outcomes,
    levels = names(expected_counts)
  ))
  valid_inputs <- purrr::map_lgl(inputs, function(input) {
    identical(input$family, "mundlak") &&
      is.character(input$formula) && length(input$formula) == 1L &&
      is.data.frame(input$data) && nrow(input$data) > 0L
  })
  regularised_four <- purrr::map_lgl(inputs, function(input) {
    input$specification$model_index != 4L ||
      isTRUE(input$specification$regularised_variance_priors)
  })
  if (
    length(pending_specifications) != 97L || length(inputs) != 97L ||
      length(reused_fits) != 7L || anyDuplicated(pending_keys) ||
      anyDuplicated(input_keys) || !setequal(pending_keys, input_keys) ||
      !identical(as.integer(observed_counts), unname(expected_counts)) ||
      !all(valid_inputs) || !all(regularised_four) ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("The full Mundlak expansion failed its structural audit.",
         call. = FALSE)
  }
  tibble::tibble(
    check = c(
      "complete candidate registry",
      "reused passing whole-body pilots",
      "new independently checkpointed fits",
      "pilot-informed sampling schedule",
      "model-4 regularised variance prior"
    ),
    status = "PASS",
    detail = c(
      "104 candidates: 26 models x 4 outcomes",
      "7 fits: models 1, 2, 3, 4, 5, 13, and 26",
      "97 fits: 19 whole-body plus 26 for each other outcome",
      "1,500 warmup and 2,000 retained iterations per chain",
      "applied prospectively to pending model-4 fits"
    )
  )
}

fit_mundlak_full_candidate <- function(input, enabled, config) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("mundlak"))
  outcome_offset <- match(input$specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the full Mundlak fit.", call. = FALSE)
  }
  local_config <- config
  local_config$seed <- as.integer(
    config$seed + 100L * (outcome_offset - 1L) +
      input$specification$model_index
  )
  set_rstan_options()
  fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "mundlak"
  )
}

collect_mundlak_full_ensemble <- function(
    reused_fits,
    new_fits,
    expected_specifications
) {
  if (
    is_sensitivity_pilot_skip(reused_fits) ||
      is_sensitivity_pilot_skip(new_fits)
  ) {
    return(new_sensitivity_pilot_skip("mundlak"))
  }
  ensemble <- c(reused_fits, new_fits)
  outcome_levels <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  ordering <- order(
    match(purrr::map_chr(ensemble, "outcome"), outcome_levels),
    purrr::map_int(ensemble, "model_index")
  )
  ensemble <- ensemble[ordering]
  observed_keys <- paste(
    purrr::map_chr(ensemble, "outcome"),
    purrr::map_int(ensemble, "model_index")
  )
  expected_keys <- paste(
    purrr::map_chr(expected_specifications, "outcome"),
    purrr::map_int(expected_specifications, "model_index")
  )
  if (
    length(ensemble) != 104L || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !all(purrr::map_lgl(ensemble, ~ identical(.x$family, "mundlak")))
  ) {
    stop("The full Mundlak ensemble is incomplete.", call. = FALSE)
  }
  ensemble
}

diagnose_mundlak_full_candidate <- function(candidate, ppc_draws = 200L) {
  if (is_sensitivity_pilot_skip(candidate)) return(candidate)
  get_primary_candidate_model_diagnostics(
    candidate,
    ppc_draws = ppc_draws
  )
}

compile_mundlak_full_diagnostic_audit <- function(diagnostics) {
  if (is_sensitivity_pilot_skip(diagnostics)) return(diagnostics)
  summary <- dplyr::bind_rows(purrr::map(diagnostics, "summary")) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    ) |>
    dplyr::mutate(
      ppc_warning = primary_ppc_warning(dplyr::pick(dplyr::everything())),
      status = dplyr::case_when(
        !diagnostic_pass ~ "REQUIRES REMEDIATION",
        ppc_warning ~ "PASS WITH PPC WARNING",
        TRUE ~ "PASS"
      )
    )
  keys <- paste(summary$outcome, summary$model_index)
  if (nrow(summary) != 104L || anyDuplicated(keys)) {
    stop("The full Mundlak diagnostic audit is incomplete.", call. = FALSE)
  }
  parameter_failures <- purrr::map_dfr(diagnostics, function(result) {
    result$parameter_failures |>
      dplyr::mutate(
        outcome = result$summary$outcome[[1L]],
        model_index = result$summary$model_index[[1L]],
        .before = 1L
      )
  })
  list(
    diagnostic_summary = summary,
    parameter_failures = parameter_failures,
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

select_mundlak_full_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  failures <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::select(outcome, model_index) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  selected <- purrr::keep(specifications, function(specification) {
    any(
      failures$outcome == specification$outcome &
        failures$model_index == specification$model_index
    )
  })
  selected_keys <- paste(
    purrr::map_chr(selected, "outcome"),
    purrr::map_int(selected, "model_index")
  )
  failure_keys <- paste(failures$outcome, failures$model_index)
  if (
    nrow(failures) < 1L || length(selected) != nrow(failures) ||
      anyDuplicated(selected_keys) || !setequal(selected_keys, failure_keys)
  ) {
    stop("The full Mundlak remedy specification set is invalid.",
         call. = FALSE)
  }
  selected
}

make_mundlak_full_remedy_config <- function(
    full_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    max_treedepth_cases = "strength 6",
    max_treedepth = 15L,
    seed_offset = 1000L
) {
  config <- full_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$max_treedepth_cases <- as.character(max_treedepth_cases)
  config$remedy_max_treedepth <- as.integer(max_treedepth)
  config$seed <- as.integer(full_config$seed + seed_offset)
  if (
    config$iter_warmup < 1L || config$iter_sampling < 1L ||
      config$remedy_max_treedepth < full_config$max_treedepth_floor
  ) {
    stop("Invalid full Mundlak remedy configuration.", call. = FALSE)
  }
  config
}

prepare_mundlak_full_remedy_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  key <- paste(specification$outcome, specification$model_index)
  specification$sensitivity_mundlak_full_remedy <- "longer_sampling"
  if (key %in% config$max_treedepth_cases) {
    specification$max_treedepth <- max(
      specification$max_treedepth,
      config$remedy_max_treedepth
    )
    specification$sensitivity_mundlak_full_remedy <-
      "longer_sampling_and_max_treedepth_15"
  }
  prepare_mundlak_full_input(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config
  )
}

sensitivity_fit_file_path <- function(
    family,
    outcome,
    model_index,
    stage
) {
  safe <- function(x) gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  directory <- here::here(
    "_targets", "user", "sensitivity_fit_files",
    safe(family), safe(outcome),
    paste0("model_", stringr::str_pad(model_index, 2L, pad = "0"))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, paste0(safe(stage), ".qs"))
}

unwrap_sensitivity_fit_candidate <- function(candidate) {
  if (
    is.list(candidate) && !"model" %in% names(candidate) &&
      length(candidate) == 1L && is.list(candidate[[1L]]) &&
      "model" %in% names(candidate[[1L]])
  ) {
    return(candidate[[1L]])
  }
  candidate
}

write_sensitivity_fit_file <- function(candidate, stage) {
  candidate <- unwrap_sensitivity_fit_candidate(candidate)
  if (is_sensitivity_pilot_skip(candidate)) return(candidate)
  if (!inherits(candidate$model, "brmsfit")) {
    stop("Cannot archive an invalid sensitivity fit.", call. = FALSE)
  }
  path <- sensitivity_fit_file_path(
    candidate$family,
    candidate$outcome,
    candidate$model_index,
    stage
  )
  temporary <- paste0(path, ".tmp")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary, preset = "high")
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    stop("Could not install the sensitivity fit file.",
         call. = FALSE)
  }
  unlink(temporary, force = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

fit_mundlak_full_remedy_to_file <- function(
    input,
    enabled,
    config,
    stage = "remedy_1"
) {
  candidate <- fit_mundlak_full_candidate(input, enabled, config)
  if (is_sensitivity_pilot_skip(candidate)) return(candidate)
  write_sensitivity_fit_file(candidate, stage = stage)
}

read_sensitivity_fit_file <- function(path) {
  if (
    length(path) != 1L || is.na(path) || !file.exists(path)
  ) {
    stop("The sensitivity fit file is unavailable.", call. = FALSE)
  }
  candidate <- qs::qread(path)
  if (
    !is.list(candidate) || !inherits(candidate$model, "brmsfit") ||
      is.null(candidate$outcome) || is.null(candidate$model_index)
  ) {
    stop("The sensitivity fit file is invalid.", call. = FALSE)
  }
  family_suffix <- paste0(" [", candidate$family, "]")
  duplicated_suffix <- paste0(family_suffix, family_suffix)
  if (
    length(candidate$model_name) == 1L &&
      stringr::str_ends(
        candidate$model_name,
        stringr::fixed(duplicated_suffix)
      )
  ) {
    candidate$model_name <- substr(
      candidate$model_name,
      1L,
      nchar(candidate$model_name) - nchar(family_suffix)
    )
  }
  candidate
}

diagnose_mundlak_fit_file <- function(path, ppc_draws = 200L) {
  candidate <- read_sensitivity_fit_file(path)
  diagnose_mundlak_full_candidate(candidate, ppc_draws = ppc_draws)
}

combine_mundlak_full_remedy_diagnostic_audit <- function(
    initial_audit,
    remedy_diagnostics
) {
  remedy_summary <- dplyr::bind_rows(purrr::map(
    remedy_diagnostics, "summary"
  ))
  remedy_keys <- paste(remedy_summary$outcome, remedy_summary$model_index)
  expected_keys <- initial_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::transmute(key = paste(outcome, model_index)) |>
    dplyr::pull(key)
  if (
    nrow(remedy_summary) != length(expected_keys) ||
      anyDuplicated(remedy_keys) || !setequal(remedy_keys, expected_keys)
  ) {
    stop("The full Mundlak remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  summary <- initial_audit$diagnostic_summary |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_summary) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    ) |>
    dplyr::mutate(
      ppc_warning = primary_ppc_warning(dplyr::pick(dplyr::everything())),
      status = dplyr::case_when(
        !diagnostic_pass ~ "REQUIRES REMEDIATION",
        ppc_warning ~ "PASS WITH PPC WARNING",
        TRUE ~ "PASS"
      )
    )
  remedy_parameter_failures <- purrr::map_dfr(
    remedy_diagnostics,
    function(result) {
      result$parameter_failures |>
        dplyr::mutate(
          outcome = result$summary$outcome[[1L]],
          model_index = result$summary$model_index[[1L]],
          .before = 1L
        )
    }
  )
  parameter_failures <- initial_audit$parameter_failures |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_parameter_failures)
  if (
    nrow(summary) != 104L ||
      anyDuplicated(paste(summary$outcome, summary$model_index))
  ) {
    stop("The remedied Mundlak diagnostic audit is incomplete.",
         call. = FALSE)
  }
  list(
    diagnostic_summary = summary,
    parameter_failures = parameter_failures,
    remedy_diagnostic_summary = remedy_summary,
    remedy_keys = remedy_keys,
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

make_sensitivity_fit_file_manifest <- function(paths, stage) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) < 1L || anyNA(paths) || any(!file.exists(paths))) {
    stop("Sensitivity fit-file paths are incomplete.", call. = FALSE)
  }
  purrr::map_dfr(paths, function(path) {
    candidate <- read_sensitivity_fit_file(path)
    tibble::tibble(
      family = candidate$family,
      outcome = candidate$outcome,
      model_index = as.integer(candidate$model_index),
      model_name = candidate$model_name,
      stage = stage,
      fit_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  })
}

write_mundlak_full_accepted_fit_file <- function(
    candidate,
    final_diagnostic_audit,
    remedy_manifest,
    stage = "initial"
) {
  candidate <- unwrap_sensitivity_fit_candidate(candidate)
  if (!isTRUE(final_diagnostic_audit$all_sampling_diagnostics_pass)) {
    stop("Mundlak fit files cannot be accepted before diagnostics pass.",
         call. = FALSE)
  }
  key <- paste(candidate$outcome, candidate$model_index)
  remedy_keys <- paste(remedy_manifest$outcome, remedy_manifest$model_index)
  remedy_match <- match(key, remedy_keys)
  if (!is.na(remedy_match)) {
    remedy_path <- remedy_manifest$fit_path[[remedy_match]]
    if (!file.exists(remedy_path)) {
      stop("An accepted Mundlak remedy file is unavailable.", call. = FALSE)
    }
    return(normalizePath(remedy_path, winslash = "/", mustWork = TRUE))
  }
  write_sensitivity_fit_file(candidate, stage = stage)
}

combine_mundlak_accepted_fit_manifest <- function(
    original_manifest,
    remedy_manifest,
    final_diagnostic_audit
) {
  if (!isTRUE(final_diagnostic_audit$all_sampling_diagnostics_pass)) {
    stop("Mundlak fits cannot be accepted before all diagnostics pass.",
         call. = FALSE)
  }
  remedy_keys <- paste(remedy_manifest$outcome, remedy_manifest$model_index)
  accepted <- original_manifest |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_manifest) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  if (
    nrow(accepted) != 104L ||
      anyDuplicated(paste(accepted$outcome, accepted$model_index)) ||
      any(!file.exists(accepted$fit_path))
  ) {
    stop("The accepted Mundlak fit manifest is incomplete.", call. = FALSE)
  }
  accepted
}

make_mundlak_full_fit_records <- function(accepted_manifest) {
  required <- c(
    "family", "outcome", "model_index", "model_name", "stage",
    "fit_path", "bytes", "sha256"
  )
  valid_outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.data.frame(accepted_manifest) ||
      !all(required %in% names(accepted_manifest)) ||
      nrow(accepted_manifest) != 104L ||
      anyDuplicated(paste(
        accepted_manifest$outcome,
        accepted_manifest$model_index
      )) ||
      !setequal(unique(accepted_manifest$outcome), valid_outcomes) ||
      any(!accepted_manifest$family %in% "mundlak") ||
      any(!file.exists(accepted_manifest$fit_path))
  ) {
    stop("The accepted Mundlak fit manifest is invalid.", call. = FALSE)
  }
  ordered <- accepted_manifest |>
    dplyr::arrange(factor(outcome, levels = valid_outcomes), model_index)
  purrr::map(seq_len(nrow(ordered)), function(row) {
    as.list(ordered[row, , drop = FALSE])
  })
}

compute_mundlak_full_file_integrated_loo <- function(
    fit_record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    family = "mundlak",
    slope_description = "mundlak_within_basis_slopes"
) {
  fit_path <- fit_record$fit_path
  if (
    length(fit_path) != 1L || is.na(fit_path) || !file.exists(fit_path)
  ) {
    stop("The Mundlak LOO fit file is unavailable.", call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(fit_path)
  if (
    !identical(unname(as.character(candidate$family)), family) ||
      !identical(candidate$outcome, fit_record$outcome) ||
      !identical(
        as.integer(candidate$model_index),
        as.integer(fit_record$model_index)
      ) ||
      !identical(candidate$model_name, fit_record$model_name) ||
      !inherits(candidate$model, "brmsfit")
  ) {
    stop("The sensitivity fit file does not match its manifest record.",
         call. = FALSE)
  }
  # `brms::loo(pointwise = TRUE)` uses the memory-saving observation-at-a-time
  # evaluator. On Windows it repeatedly creates a PSOCK cluster according to
  # the global `mc.cores` option, which is exceptionally slow for these small
  # 4,000-draw likelihood matrices. Keep the calculation local and construct
  # the full log-likelihood matrix once; even the largest outcome requires
  # only a few million numeric entries.
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  draw_ids <- make_sensitivity_loo_draw_ids(
    brms::ndraws(candidate$model),
    min(as.integer(loo_draws), as.integer(brms::ndraws(candidate$model)))
  )
  conditional_loo <- brms::loo(
    candidate$model,
    draw_ids = draw_ids,
    pointwise = FALSE,
    moment_match = FALSE,
    reloo = FALSE,
    cores = as.integer(cores)
  )
  result <- compute_sensitivity_mundlak_integrated_psis_loo(
    job = list(
      fitted_candidate = candidate,
      conditional_loo = conditional_loo,
      loo_draws = length(draw_ids)
    ),
    k_threshold = k_threshold,
    cores = cores,
    family = family,
    slope_description = slope_description
  )
  result$family <- family
  result$fit_record <- fit_record
  result$loo_draws <- length(draw_ids)
  result
}

compute_exposure_full_file_integrated_loo <- function(
    fit_record,
    family,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L
) {
  if (!family %in% c("achieved", "assignment")) {
    stop("Unknown exposure family for integrated LOO.", call. = FALSE)
  }
  fit_path <- fit_record$fit_path
  if (
    length(fit_path) != 1L || is.na(fit_path) || !file.exists(fit_path)
  ) {
    stop("The exposure LOO fit file is unavailable.", call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(fit_path)
  if (
    !identical(unname(as.character(candidate$family)), family) ||
      !identical(candidate$outcome, fit_record$outcome) ||
      !identical(
        as.integer(candidate$model_index),
        as.integer(fit_record$model_index)
      ) ||
      !identical(candidate$model_name, fit_record$model_name) ||
      !inherits(candidate$model, "brmsfit")
  ) {
    stop("The exposure fit file does not match its manifest record.",
         call. = FALSE)
  }
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  draw_ids <- make_sensitivity_loo_draw_ids(
    brms::ndraws(candidate$model),
    min(as.integer(loo_draws), as.integer(brms::ndraws(candidate$model)))
  )
  conditional_loo <- brms::loo(
    candidate$model,
    draw_ids = draw_ids,
    pointwise = FALSE,
    moment_match = FALSE,
    reloo = FALSE,
    cores = as.integer(cores)
  )
  result <- compute_sensitivity_exposure_integrated_psis_loo(
    job = list(
      fitted_candidate = candidate,
      conditional_loo = conditional_loo,
      loo_draws = length(draw_ids)
    ),
    family = family,
    k_threshold = k_threshold,
    cores = cores
  )
  result$family <- family
  result$fit_record <- fit_record
  result$loo_draws <- length(draw_ids)
  result
}

collect_mundlak_full_integrated_loo <- function(
    results,
    expected_indices = seq_len(26L)
) {
  valid_outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  expected_indices <- as.integer(expected_indices)
  expected_models <- length(valid_outcomes) * length(expected_indices)
  if (
    !is.list(results) || length(results) != expected_models ||
      length(expected_indices) < 1L || anyNA(expected_indices) ||
      anyDuplicated(expected_indices)
  ) {
    stop("The full Mundlak integrated-LOO result set is incomplete.",
         call. = FALSE)
  }
  keys <- purrr::map_chr(
    results,
    ~ paste(.x$outcome, .x$model_index)
  )
  valid <- purrr::map_lgl(results, function(result) {
    inherits(result$loo, "loo") &&
      inherits(result$conditional_loo, "loo") &&
      result$outcome %in% valid_outcomes &&
      result$model_index %in% expected_indices &&
      is.list(result$fit_record) &&
      file.exists(result$fit_record$fit_path)
  })
  if (
    anyDuplicated(keys) || !all(valid) ||
      !setequal(keys, as.vector(outer(
        valid_outcomes,
        expected_indices,
        paste
      )))
  ) {
    stop("The full Mundlak integrated-LOO results failed validation.",
         call. = FALSE)
  }
  results[order(
    match(purrr::map_chr(results, "outcome"), valid_outcomes),
    purrr::map_int(results, "model_index")
  )]
}

audit_mundlak_full_integrated_loo <- function(
    results,
    k_threshold = 0.7,
    expected_indices = seq_len(26L)
) {
  results <- collect_mundlak_full_integrated_loo(
    results,
    expected_indices = expected_indices
  )
  rescor_cluster_loo <- all(purrr::map_lgl(
    results,
    ~ unname(as.character(.x$family)) %in% c(
      "rescor", "total_rescor", "supplement_rescor"
    )
  ))
  conditional_method <- if (rescor_cluster_loo) {
    "conditional_study_cluster_psis"
  } else {
    "conditional_observation_level_psis"
  }
  integrated_method <- if (rescor_cluster_loo) {
    "integrated_study_cluster_psis"
  } else {
    "integrated_observation_level_psis"
  }
  model_summary <- purrr::map_dfr(results, function(result) {
    dplyr::bind_rows(
      summarise_primary_loo_method(
        result$conditional_loo,
        conditional_method,
        k_threshold
      ),
      summarise_primary_loo_method(
        result$loo,
        integrated_method,
        k_threshold
      )
    ) |>
      dplyr::mutate(
        outcome = result$outcome,
        model_index = result$model_index,
        model_name = result$model_name,
        loo_draws = result$loo_draws,
        .before = 1L
      )
  })
  high_k_cases <- purrr::map_dfr(results, function(result) {
    pareto_k <- as.numeric(loo::pareto_k_values(result$loo))
    observations <- which(pareto_k > k_threshold)
    if (length(observations) < 1L) return(tibble::tibble())
    result$observation_manifest[observations, , drop = FALSE] |>
      dplyr::mutate(
        outcome = result$outcome,
        model_index = result$model_index,
        model_name = result$model_name,
        observation = observations,
        pareto_k = pareto_k[observations],
        fit_path = result$fit_record$fit_path,
        .before = 1L
      )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index,
      observation
    )
  outcome_summary <- model_summary |>
    dplyr::filter(method == integrated_method) |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_without_high_k = sum(n_pareto_k_gt_threshold == 0L),
      models_requiring_exact_refits = sum(n_pareto_k_gt_threshold > 0L),
      exact_refits_required = sum(n_pareto_k_gt_threshold),
      maximum_pareto_k = max(maximum_pareto_k),
      .groups = "drop"
    )
  list(
    model_summary = model_summary,
    outcome_summary = outcome_summary,
    high_k_cases = high_k_cases,
    models = length(results),
    exact_refits_required = nrow(high_k_cases),
    integrated_psis_complete = TRUE,
    stacking_ready = nrow(high_k_cases) == 0L,
    loo_unit = if (rescor_cluster_loo) "study" else "observation",
    k_threshold = k_threshold
  )
}

sensitivity_reduction_file_path <- function(
    family,
    outcome,
    model_index,
    stage,
    suffix = NULL
) {
  safe <- function(x) gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  directory <- here::here(
    "_targets", "user", "sensitivity_reductions",
    safe(family), safe(outcome),
    paste0("model_", stringr::str_pad(model_index, 2L, pad = "0"))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- safe(stage)
  if (!is.null(suffix) && nzchar(suffix)) {
    filename <- paste(filename, safe(suffix), sep = "_")
  }
  file.path(directory, paste0(filename, ".qs"))
}

write_atomic_qs <- function(object, path, preset = "high") {
  temporary <- paste0(path, ".tmp")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  qs::qsave(object, temporary, preset = preset)
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    stop("Could not install a sensitivity-analysis file.", call. = FALSE)
  }
  unlink(temporary, force = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_mundlak_full_integrated_loo_file <- function(result) {
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(result$remedy, "integrated_observation_level_psis")
  ) {
    stop("Cannot archive an invalid integrated Mundlak LOO result.",
         call. = FALSE)
  }
  path <- sensitivity_reduction_file_path(
    family = "mundlak",
    outcome = result$outcome,
    model_index = result$model_index,
    stage = "integrated_loo"
  )
  write_atomic_qs(result, path)
}

make_mundlak_full_integrated_loo_file_manifest <- function(paths) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != 104L || anyNA(paths) || any(!file.exists(paths))) {
    stop("The Mundlak integrated-LOO file set is incomplete.", call. = FALSE)
  }
  manifest <- purrr::map_dfr(paths, function(path) {
    result <- qs::qread(path)
    if (!inherits(result$loo, "loo")) {
      stop("A Mundlak integrated-LOO file is invalid.", call. = FALSE)
    }
    tibble::tibble(
      family = "mundlak",
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      integrated_loo_path = normalizePath(
        path, winslash = "/", mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index))
  ) {
    stop("The Mundlak integrated-LOO file manifest is invalid.",
         call. = FALSE)
  }
  manifest
}

write_sensitivity_full_integrated_loo_file <- function(result, family) {
  valid_remedy <- if (identical(family, "eiv")) {
    result$remedy %in% c(
      "outcome_effects_integrated_latent_exposure_psis",
      "outcome_effects_integrated_psis_plus_exact"
    )
  } else if (identical(family, "rescor")) {
    result$remedy %in% c(
      "study_cluster_integrated_psis",
      "study_cluster_integrated_psis_plus_exact"
    )
  } else {
    identical(result$remedy, "integrated_observation_level_psis")
  }
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !isTRUE(valid_remedy) ||
      !identical(unname(as.character(result$family)), family)
  ) {
    stop("Cannot archive an invalid integrated sensitivity LOO result.",
         call. = FALSE)
  }
  path <- sensitivity_reduction_file_path(
    family = family,
    outcome = result$outcome,
    model_index = result$model_index,
    stage = "integrated_loo"
  )
  write_atomic_qs(result, path)
}

make_sensitivity_full_integrated_loo_file_manifest <- function(
    paths,
    family,
    expected_models = 104L,
    expected_indices = seq_len(26L)
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (
    length(paths) != expected_models || anyNA(paths) ||
      any(!file.exists(paths))
  ) {
    stop("The integrated sensitivity LOO file set is incomplete.",
         call. = FALSE)
  }
  manifest <- purrr::map_dfr(paths, function(path) {
    result <- qs::qread(path)
    if (
      !inherits(result$loo, "loo") ||
        !identical(unname(as.character(result$family)), family)
    ) {
      stop("An integrated sensitivity LOO file is invalid.",
           call. = FALSE)
    }
    tibble::tibble(
      family = family,
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      integrated_loo_path = normalizePath(
        path, winslash = "/", mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  if (
    nrow(manifest) != expected_models ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(manifest$family != family) ||
      !setequal(
        paste(manifest$outcome, manifest$model_index, sep = "::"),
        as.vector(outer(
          c("whole_body", "hypertrophy", "strength", "function"),
          as.integer(expected_indices),
          function(outcome, model_index) {
            paste(outcome, model_index, sep = "::")
          }
        ))
      )
  ) {
    stop("The integrated sensitivity LOO manifest is invalid.",
         call. = FALSE)
  }
  manifest
}

audit_sensitivity_full_integrated_loo_files <- function(
    manifest,
    family,
    k_threshold = 0.7,
    expected_models = 104L,
    expected_indices = seq_len(26L)
) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != expected_models ||
      any(manifest$family != family) ||
      any(!file.exists(manifest$integrated_loo_path))
  ) {
    stop("The integrated sensitivity LOO manifest cannot be audited.",
         call. = FALSE)
  }
  results <- purrr::map(manifest$integrated_loo_path, qs::qread)
  valid_family <- purrr::map_lgl(
    results,
    ~ identical(unname(as.character(.x$family)), family)
  )
  if (!all(valid_family)) {
    stop("The integrated sensitivity LOO files have mixed families.",
         call. = FALSE)
  }
  audit <- audit_mundlak_full_integrated_loo(
    results,
    k_threshold = k_threshold,
    expected_indices = expected_indices
  )
  audit$family <- family
  audit
}

make_mundlak_full_exact_loo_records <- function(
    integrated_audit,
    accepted_fit_manifest,
    integrated_loo_manifest
) {
  cases <- integrated_audit$high_k_cases
  expected_cases <- as.integer(integrated_audit$exact_refits_required)
  if (
    !is.data.frame(cases) || length(expected_cases) != 1L ||
      is.na(expected_cases) || expected_cases < 1L ||
      nrow(cases) != expected_cases
  ) {
    stop("The Mundlak exact-LOO case audit has changed.", call. = FALSE)
  }
  fit_lookup <- accepted_fit_manifest |>
    dplyr::select(
      outcome, model_index,
      original_fit_path = fit_path,
      original_fit_sha256 = sha256
    )
  loo_lookup <- integrated_loo_manifest |>
    dplyr::select(
      outcome, model_index, integrated_loo_path,
      integrated_loo_sha256 = sha256
    )
  records <- cases |>
    dplyr::select(outcome, model_index, model_name, observation, pareto_k) |>
    dplyr::left_join(
      fit_lookup,
      by = c("outcome", "model_index"),
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      loo_lookup,
      by = c("outcome", "model_index"),
      relationship = "many-to-one"
    ) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index,
      observation
    ) |>
    dplyr::mutate(case_id = dplyr::row_number(), .before = 1L)
  if (
    nrow(records) != expected_cases ||
      anyNA(records$original_fit_path) ||
      anyNA(records$integrated_loo_path) ||
      any(!file.exists(records$original_fit_path)) ||
      any(!file.exists(records$integrated_loo_path)) ||
      anyDuplicated(paste(
        records$outcome,
        records$model_index,
        records$observation
      ))
  ) {
    stop("The Mundlak exact-LOO records are incomplete.", call. = FALSE)
  }
  purrr::map(seq_len(nrow(records)), function(row) {
    as.list(records[row, , drop = FALSE])
  })
}

make_mundlak_full_exact_loo_config <- function(
    base_config,
    iter_warmup = 2000L,
    iter_sampling = 6000L,
    seed_offset = 4000L
) {
  config <- base_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains < 1L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L
  ) {
    stop("The full Mundlak exact-LOO configuration is invalid.",
         call. = FALSE)
  }
  config
}

mundlak_full_exact_fit_file_path <- function(record) {
  suffix <- paste0(
    "case_", stringr::str_pad(record$case_id, 2L, pad = "0"),
    "_observation_", record$observation
  )
  sensitivity_fit_file_path(
    family = "mundlak_exact",
    outcome = record$outcome,
    model_index = record$model_index,
    stage = suffix
  )
}

fit_mundlak_full_exact_loo_to_file <- function(record, enabled, config) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("mundlak"))
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated_result <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !inherits(candidate$model, "brmsfit") ||
      !identical(candidate$outcome, record$outcome) ||
      !identical(as.integer(candidate$model_index),
                 as.integer(record$model_index)) ||
      !inherits(integrated_result$loo, "loo") ||
      !identical(integrated_result$outcome, record$outcome) ||
      !identical(as.integer(integrated_result$model_index),
                 as.integer(record$model_index)) ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      integrated_result$loo$diagnostics$pareto_k[[observation]] <= 0.7
  ) {
    stop("A full Mundlak exact-LOO record failed validation.",
         call. = FALSE)
  }
  prediction_data <- candidate$model$data
  training_data <- prediction_data[-observation, , drop = FALSE]
  original_model <- candidate$model
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
    chains = config$chains,
    cores = config$parallel_chains,
    warmup = config$iter_warmup,
    iter = config$iter_warmup + config$iter_sampling,
    seed = config$seed + as.integer(record$case_id),
    control = list(
      adapt_delta = candidate$specification$adapt_delta,
      max_treedepth = candidate$specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  archive <- list(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    old_psis_pareto_k = as.numeric(record$pareto_k),
    original_fit_path = record$original_fit_path,
    integrated_loo_path = record$integrated_loo_path,
    prediction_data = prediction_data,
    specification = candidate$specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  write_atomic_qs(archive, mundlak_full_exact_fit_file_path(record))
}

read_mundlak_full_exact_loo_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full Mundlak exact-LOO file is unavailable.", call. = FALSE)
  }
  refit <- qs::qread(path)
  if (
    !is.list(refit) || !inherits(refit$refitted_model, "brmsfit") ||
      is.null(refit$case_id) || is.null(refit$observation) ||
      !file.exists(refit$integrated_loo_path)
  ) {
    stop("A full Mundlak exact-LOO file is invalid.", call. = FALSE)
  }
  refit
}

diagnose_mundlak_full_exact_loo_file <- function(path) {
  refit <- read_mundlak_full_exact_loo_file(path)
  summarise_primary_exact_loo_refit_diagnostics(refit)
}

audit_mundlak_full_exact_loo_diagnostics <- function(
    diagnostics,
    integrated_audit
) {
  summary <- dplyr::bind_rows(diagnostics) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index,
      observation
    )
  if (
    nrow(summary) != integrated_audit$exact_refits_required ||
      anyDuplicated(paste(summary$outcome, summary$model_index,
                          summary$observation))
  ) {
    stop("The full Mundlak exact-refit diagnostic set is incomplete.",
         call. = FALSE)
  }
  list(
    diagnostic_summary = summary,
    exact_refits = nrow(summary),
    exact_refits_passing = sum(summary$diagnostic_pass),
    exact_refits_requiring_remediation = sum(!summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

compute_mundlak_full_exact_loo_file <- function(
    path,
    diagnostic_audit
) {
  if (!isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)) {
    stop("Exact Mundlak corrections require passing refit diagnostics.",
         call. = FALSE)
  }
  refit <- read_mundlak_full_exact_loo_file(path)
  refit$original_data <- refit$prediction_data
  refit$integrated_result <- qs::qread(refit$integrated_loo_path)
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  compute_sensitivity_mundlak_exact_loo_case(refit)
}

sensitivity_full_exact_fit_file_path <- function(record, family) {
  suffix <- paste0(
    "case_", stringr::str_pad(record$case_id, 2L, pad = "0"),
    "_observation_", record$observation
  )
  sensitivity_fit_file_path(
    family = paste0(family, "_exact"),
    outcome = record$outcome,
    model_index = record$model_index,
    stage = suffix
  )
}

fit_sensitivity_full_exact_loo_to_file <- function(
    record,
    enabled,
    config,
    family
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip(family))
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated_result <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !inherits(candidate$model, "brmsfit") ||
      !identical(unname(as.character(candidate$family)), family) ||
      !identical(unname(as.character(candidate$outcome)),
                 unname(as.character(record$outcome))) ||
      !identical(unname(as.integer(candidate$model_index)),
                 unname(as.integer(record$model_index))) ||
      !inherits(integrated_result$loo, "loo") ||
      !identical(unname(as.character(integrated_result$family)), family) ||
      !identical(unname(as.character(integrated_result$outcome)),
                 unname(as.character(record$outcome))) ||
      !identical(unname(as.integer(integrated_result$model_index)),
                 unname(as.integer(record$model_index))) ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      integrated_result$loo$diagnostics$pareto_k[[observation]] <= 0.7
  ) {
    stop("A full sensitivity exact-LOO record failed validation.",
         call. = FALSE)
  }
  prediction_data <- candidate$model$data
  training_data <- prediction_data[-observation, , drop = FALSE]
  original_model <- candidate$model
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
    chains = config$chains,
    cores = config$parallel_chains,
    warmup = config$iter_warmup,
    iter = config$iter_warmup + config$iter_sampling,
    seed = config$seed + as.integer(record$case_id),
    control = list(
      adapt_delta = max(
        candidate$specification$adapt_delta,
        config$adapt_delta_floor
      ),
      max_treedepth = max(
        candidate$specification$max_treedepth,
        config$max_treedepth_floor
      )
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  archive <- list(
    family = family,
    case_id = as.integer(record$case_id),
    outcome = unname(as.character(record$outcome)),
    model_index = unname(as.integer(record$model_index)),
    model_name = unname(as.character(record$model_name)),
    observation = observation,
    old_psis_pareto_k = as.numeric(record$pareto_k),
    original_fit_path = unname(as.character(record$original_fit_path)),
    integrated_loo_path = unname(as.character(record$integrated_loo_path)),
    prediction_data = prediction_data,
    specification = candidate$specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  write_atomic_qs(
    archive,
    sensitivity_full_exact_fit_file_path(record, family)
  )
}

read_sensitivity_full_exact_loo_file <- function(path, family) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full sensitivity exact-LOO file is unavailable.",
         call. = FALSE)
  }
  refit <- qs::qread(path)
  if (
    !is.list(refit) || !inherits(refit$refitted_model, "brmsfit") ||
      !identical(unname(as.character(refit$family)), family) ||
      is.null(refit$case_id) || is.null(refit$observation) ||
      !file.exists(refit$integrated_loo_path)
  ) {
    stop("A full sensitivity exact-LOO file is invalid.",
         call. = FALSE)
  }
  refit
}

diagnose_sensitivity_full_exact_loo_file <- function(path, family) {
  refit <- read_sensitivity_full_exact_loo_file(path, family)
  summarise_primary_exact_loo_refit_diagnostics(refit)
}

compute_sensitivity_full_exact_loo_file <- function(
    path,
    diagnostic_audit,
    family
) {
  if (!isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)) {
    stop("Exact sensitivity corrections require passing refit diagnostics.",
         call. = FALSE)
  }
  refit <- read_sensitivity_full_exact_loo_file(path, family)
  refit$original_data <- refit$prediction_data
  refit$integrated_result <- qs::qread(refit$integrated_loo_path)
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  compute_sensitivity_mundlak_exact_loo_case(refit)
}

apply_sensitivity_full_exact_loo_corrections_from_files <- function(
    integrated_manifest,
    exact_results,
    family,
    expected_models = 104L,
    expected_indices = seq_len(26L)
) {
  if (
    !is.data.frame(integrated_manifest) ||
      nrow(integrated_manifest) != expected_models ||
      any(integrated_manifest$family != family) ||
      any(!file.exists(integrated_manifest$integrated_loo_path))
  ) {
    stop("The integrated sensitivity results are unavailable.",
         call. = FALSE)
  }
  integrated_results <- purrr::map(
    integrated_manifest$integrated_loo_path,
    qs::qread
  )
  audit <- apply_mundlak_full_exact_loo_corrections(
    integrated_results,
    exact_results,
    expected_indices = expected_indices
  )
  audit$family <- family
  audit
}

compute_sensitivity_full_outcome_weights <- function(
    final_loo_audit,
    family,
    expected_models_per_outcome = 26L
) {
  if (!isTRUE(final_loo_audit$stacking_ready)) {
    stop("The full sensitivity LOO audit is not stacking-ready.",
         call. = FALSE)
  }
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  weights <- purrr::map_dfr(outcomes, function(outcome) {
    results <- purrr::keep(
      final_loo_audit$results,
      ~ identical(.x$outcome, outcome)
    )
    if (length(results) != expected_models_per_outcome) {
      stop("An outcome-specific sensitivity LOO ensemble is incomplete.",
           call. = FALSE)
    }
    loo_objects <- purrr::map(results, "loo")
    names(loo_objects) <- purrr::map_chr(results, "model_name")
    values <- loo::loo_model_weights(loo_objects, method = "stacking")
    tibble::tibble(
      family = family,
      outcome = outcome,
      model_index = purrr::map_int(results, "model_index"),
      model_name = purrr::map_chr(results, "model_name"),
      weight = as.numeric(values) / sum(as.numeric(values)),
      method = "stacking_exact_corrected_integrated_loo"
    )
  })
  if (
    nrow(weights) != length(outcomes) * expected_models_per_outcome ||
      any(abs(tapply(weights$weight, weights$outcome, sum) - 1) > 1e-8)
  ) {
    stop("The full sensitivity stacking weights are invalid.",
         call. = FALSE)
  }
  weights
}

collect_mundlak_full_exact_loo_results <- function(results, records) {
  if (length(results) != length(records) || length(results) < 1L) {
    stop("The full Mundlak exact-correction set is incomplete.",
         call. = FALSE)
  }
  results <- results[order(purrr::map_int(results, "case_id"))]
  expected <- sort(purrr::map_int(records, "case_id"))
  if (!identical(unname(purrr::map_int(results, "case_id")), expected)) {
    stop("The full Mundlak exact-correction identities do not match.",
         call. = FALSE)
  }
  results
}

apply_mundlak_full_exact_loo_corrections <- function(
    integrated_results,
    exact_results,
    expected_indices = seq_len(26L)
) {
  integrated_results <- collect_mundlak_full_integrated_loo(
    integrated_results,
    expected_indices = expected_indices
  )
  warning_cases <- purrr::map_dfr(integrated_results, function(result) {
    pareto_k <- loo::pareto_k_values(result$loo)
    observations <- which(pareto_k > 0.7)
    tibble::tibble(
      outcome = result$outcome,
      model_index = result$model_index,
      observation = observations
    )
  }) |>
    dplyr::arrange(outcome, model_index, observation)
  exact_cases <- purrr::map_dfr(exact_results, function(result) {
    tibble::tibble(
      outcome = result$outcome,
      model_index = result$model_index,
      observation = result$exact_refit_summary$observation[[1L]]
    )
  }) |>
    dplyr::arrange(outcome, model_index, observation)
  if (!identical(warning_cases, exact_cases)) {
    stop("Exact Mundlak corrections do not match the full high-k cases.",
         call. = FALSE)
  }
  final_results <- integrated_results
  keys <- unique(paste(exact_cases$outcome, exact_cases$model_index))
  for (key in keys) {
    position <- match(
      key,
      paste(
        purrr::map_chr(final_results, "outcome"),
        purrr::map_int(final_results, "model_index")
      )
    )
    patched <- final_results[[position]]
    corrections <- purrr::keep(
      exact_results,
      ~ identical(paste(.x$outcome, .x$model_index), key)
    )
    for (correction in corrections) {
      observation <- correction$exact_refit_summary$observation[[1L]]
      patched$loo$pointwise[observation, ] <-
        correction$loo$pointwise[observation, ]
      patched$loo$diagnostics$pareto_k[[observation]] <- NA_real_
      if (!is.null(patched$loo$diagnostics$n_eff)) {
        patched$loo$diagnostics$n_eff[[observation]] <- NA_real_
      }
    }
    patched$loo <- update_primary_exact_loo_estimates(patched$loo)
    patched$remedy <- "integrated_observation_level_psis_plus_exact"
    final_results[[position]] <- patched
  }
  exact_summary <- purrr::map_dfr(exact_results, "exact_refit_summary")
  model_summary <- purrr::map_dfr(final_results, function(result) {
    summarise_primary_loo_method(
      result$loo,
      "integrated_observation_level_psis_plus_exact"
    ) |>
      dplyr::mutate(
        outcome = result$outcome,
        model_index = result$model_index,
        model_name = result$model_name,
        .before = 1L
      )
  })
  list(
    results = final_results,
    model_summary = model_summary,
    exact_refit_summary = exact_summary,
    stacking_ready =
      all(model_summary$n_pareto_k_gt_threshold == 0L) &&
      all(exact_summary$diagnostic_pass)
  )
}

compute_mundlak_full_outcome_weights <- function(
    final_loo_audit,
    expected_models_per_outcome = 26L
) {
  if (!isTRUE(final_loo_audit$stacking_ready)) {
    stop("The full Mundlak LOO audit is not stacking-ready.", call. = FALSE)
  }
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  weights <- purrr::map_dfr(outcomes, function(outcome) {
    results <- purrr::keep(
      final_loo_audit$results,
      ~ identical(.x$outcome, outcome)
    )
    if (length(results) != expected_models_per_outcome) {
      stop("An outcome-specific Mundlak LOO ensemble is incomplete.",
           call. = FALSE)
    }
    loo_objects <- purrr::map(results, "loo")
    names(loo_objects) <- purrr::map_chr(results, "model_name")
    values <- loo::loo_model_weights(loo_objects, method = "stacking")
    tibble::tibble(
      outcome = outcome,
      model_index = purrr::map_int(results, "model_index"),
      model_name = purrr::map_chr(results, "model_name"),
      weight = as.numeric(values) / sum(as.numeric(values)),
      method = "stacking_exact_corrected_integrated_loo"
    )
  })
  if (
    nrow(weights) != length(outcomes) * expected_models_per_outcome ||
      any(abs(tapply(weights$weight, weights$outcome, sum) - 1) > 1e-8)
  ) {
    stop("The full Mundlak stacking weights are invalid.", call. = FALSE)
  }
  weights
}

read_mundlak_completed_checkpoint <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(
      "The completed Mundlak checkpoint is unavailable. Run ",
      "sensitivity/freeze_mundlak_checkpoint.R after the final LOO audit.",
      call. = FALSE
    )
  }
  checkpoint <- qs::qread(path)
  required <- c(
    "checkpoint_version", "family", "analysis_data_by_outcome",
    "model_specifications", "accepted_fit_manifest",
    "sampling_diagnostic_summary", "sampling_parameter_failures",
    "integrated_loo_model_summary", "integrated_loo_outcome_summary",
    "initial_high_pareto_k_cases", "exact_refit_diagnostic_summary",
    "exact_loo_correction_summary", "final_loo_model_summary",
    "stacking_weights", "audit"
  )
  if (
    !is.list(checkpoint) || !all(required %in% names(checkpoint)) ||
      !identical(checkpoint$checkpoint_version, 1L) ||
      !identical(checkpoint$family, "mundlak") ||
      !isTRUE(checkpoint$audit$all_sampling_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$stacking_ready) ||
      nrow(checkpoint$accepted_fit_manifest) != 104L ||
      nrow(checkpoint$stacking_weights) != 104L
  ) {
    stop("The completed Mundlak checkpoint is invalid.", call. = FALSE)
  }
  checkpoint
}

make_mundlak_full_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    n_draws = 4000L,
    protein_min = 0.8,
    grid_step = 0.01,
    seed = 1988L
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  required_fit <- c(
    "outcome", "model_index", "model_name", "fit_path", "sha256"
  )
  required_weight <- c(
    "outcome", "model_index", "model_name", "weight", "method"
  )
  n_draws <- as.integer(n_draws)
  if (
    !is.data.frame(accepted_fit_manifest) ||
      !all(required_fit %in% names(accepted_fit_manifest)) ||
      nrow(accepted_fit_manifest) != 104L ||
      !is.data.frame(stacking_weights) ||
      !all(required_weight %in% names(stacking_weights)) ||
      nrow(stacking_weights) != 104L ||
      length(n_draws) != 1L || is.na(n_draws) || n_draws < 1L ||
      length(protein_min) != 1L || !is.finite(protein_min) ||
      length(grid_step) != 1L || !is.finite(grid_step) || grid_step <= 0
  ) {
    stop("The full Mundlak prediction inputs are invalid.", call. = FALSE)
  }
  fit_manifest <- accepted_fit_manifest |>
    dplyr::select(dplyr::all_of(required_fit)) |>
    dplyr::rename(fit_sha256 = sha256)
  joined <- stacking_weights |>
    dplyr::left_join(
      fit_manifest,
      by = c("outcome", "model_index", "model_name"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  if (
    anyNA(joined$fit_path) || any(!file.exists(joined$fit_path)) ||
      anyDuplicated(paste(joined$outcome, joined$model_index)) ||
      !setequal(unique(joined$outcome), outcomes) ||
      any(!is.finite(joined$weight)) || any(joined$weight < 0)
  ) {
    stop("The full Mundlak fit and weight registries do not match.",
         call. = FALSE)
  }
  outcome_seed <- stats::setNames(seed + seq_along(outcomes) * 10000L,
                                  outcomes)
  allocated <- joined |>
    dplyr::group_by(outcome) |>
    dplyr::group_modify(function(table, key) {
      if (nrow(table) != 26L || abs(sum(table$weight) - 1) > 1e-8) {
        stop("An outcome-specific Mundlak prediction ensemble is invalid.",
             call. = FALSE)
      }
      table$draw_count <- allocate_primary_mixture_models(
        table$weight / sum(table$weight), n_draws
      )
      table
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      prediction_seed = as.integer(
        outcome_seed[outcome] + as.integer(model_index)
      ),
      n_draws = n_draws,
      protein_min = protein_min,
      grid_step = grid_step
    )
  allocation_audit <- allocated |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      allocated_draws = sum(draw_count),
      selected_models = sum(draw_count > 0L),
      .groups = "drop"
    )
  if (any(allocation_audit$allocated_draws != n_draws)) {
    stop("The full Mundlak draw allocation is incomplete.", call. = FALSE)
  }
  allocated <- dplyr::filter(allocated, draw_count > 0L)
  purrr::map(seq_len(nrow(allocated)), function(row) {
    purrr::map(
      as.list(allocated[row, , drop = FALSE]),
      unname
    )
  })
}

write_mundlak_full_prediction_draw_file <- function(
    record,
    checkpoint_path
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L
  ) {
    stop("A Mundlak prediction record is invalid.", call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "mundlak") ||
      !identical(
        unname(as.character(candidate$outcome)),
        unname(as.character(record$outcome))
      ) ||
      !identical(
        unname(as.integer(candidate$model_index)),
        unname(as.integer(record$model_index))
      ) ||
      !identical(
        unname(as.character(candidate$model_name)),
        unname(as.character(record$model_name))
      )
  ) {
    stop("A Mundlak prediction fit does not match its record.",
         call. = FALSE)
  }
  checkpoint <- read_mundlak_completed_checkpoint(checkpoint_path)
  prediction_data <- checkpoint$analysis_data_by_outcome[[record$outcome]]
  if (
    is.null(prediction_data) || !is.data.frame(prediction_data) ||
      !all(c("protein_intake_raw", "protein_median") %in%
           names(prediction_data))
  ) {
    stop("The original Mundlak prediction data are unavailable.",
         call. = FALSE)
  }
  grid <- make_mundlak_prediction_newdata(
    fitted_candidate = candidate,
    data = prediction_data,
    protein_min = record$protein_min,
    protein_max = NULL,
    grid_step = record$grid_step
  )
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  set.seed(as.integer(record$prediction_seed))
  draws <- as.matrix(brms::posterior_epred(
    candidate$model,
    newdata = grid$newdata,
    re_formula = NA,
    ndraws = as.integer(record$draw_count)
  ))
  if (
    ncol(draws) != nrow(grid$keys) && nrow(draws) == nrow(grid$keys)
  ) {
    draws <- t(draws)
  }
  if (
    nrow(draws) != as.integer(record$draw_count) ||
      ncol(draws) != nrow(grid$keys) || anyNA(draws)
  ) {
    stop("Mundlak prediction draws do not match their allocation.",
         call. = FALSE)
  }
  reduction <- list(
    family = "mundlak",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    weight = as.numeric(record$weight),
    draw_count = as.integer(record$draw_count),
    prediction_seed = as.integer(record$prediction_seed),
    fit_sha256 = record$fit_sha256,
    grid = grid$keys,
    draws = draws,
    settings = list(
      protein_median = grid$protein_median,
      protein_min = grid$protein_min,
      protein_max = grid$protein_max,
      grid_step = grid$grid_step,
      prediction_type = "mundlak"
    )
  )
  path <- sensitivity_reduction_file_path(
    family = "mundlak",
    outcome = record$outcome,
    model_index = record$model_index,
    stage = "prediction_draws"
  )
  write_atomic_qs(reduction, path)
}

read_mundlak_full_prediction_draw_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full Mundlak prediction-draw file is unavailable.",
         call. = FALSE)
  }
  reduction <- qs::qread(path)
  if (
    !is.list(reduction) || !identical(reduction$family, "mundlak") ||
      !is.matrix(reduction$draws) ||
      nrow(reduction$draws) != reduction$draw_count ||
      ncol(reduction$draws) != nrow(reduction$grid)
  ) {
    stop("A full Mundlak prediction-draw file is invalid.",
         call. = FALSE)
  }
  reduction
}

make_mundlak_full_prediction_file_manifest <- function(paths) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) < 4L || anyNA(paths) || any(!file.exists(paths))) {
    stop("The full Mundlak prediction-file set is incomplete.",
         call. = FALSE)
  }
  manifest <- purrr::map_dfr(paths, function(path) {
    reduction <- read_mundlak_full_prediction_draw_file(path)
    tibble::tibble(
      family = reduction$family,
      outcome = reduction$outcome,
      model_index = reduction$model_index,
      model_name = reduction$model_name,
      weight = reduction$weight,
      draw_count = reduction$draw_count,
      prediction_seed = reduction$prediction_seed,
      fit_sha256 = reduction$fit_sha256,
      prediction_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  audit <- manifest |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      allocated_draws = sum(draw_count),
      selected_models = dplyr::n(),
      .groups = "drop"
    )
  if (
    nrow(audit) != 4L || any(audit$allocated_draws != 4000L) ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index))
  ) {
    stop("The full Mundlak prediction-file manifest is invalid.",
         call. = FALSE)
  }
  manifest
}

make_mundlak_full_model_average_records <- function(
    prediction_manifest,
    stacking_weights,
    seed = 1988L
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.data.frame(prediction_manifest) ||
      !is.data.frame(stacking_weights) || nrow(stacking_weights) != 104L
  ) {
    stop("The full Mundlak averaging inputs are invalid.", call. = FALSE)
  }
  purrr::map(seq_along(outcomes), function(index) {
    outcome <- outcomes[[index]]
    files <- dplyr::filter(
      prediction_manifest,
      .data$outcome == .env$outcome
    )
    weights <- dplyr::filter(
      stacking_weights,
      .data$outcome == .env$outcome
    ) |>
      dplyr::arrange(model_index)
    if (
      nrow(files) < 1L || sum(files$draw_count) != 4000L ||
        nrow(weights) != 26L || abs(sum(weights$weight) - 1) > 1e-8 ||
        any(!file.exists(files$prediction_path))
    ) {
      stop("An outcome-specific Mundlak model-average record is invalid.",
           call. = FALSE)
    }
    list(
      family = "mundlak",
      outcome = outcome,
      prediction_manifest = files,
      stacking_weights = weights,
      n_draws = 4000L,
      averaging_seed = as.integer(seed + index)
    )
  })
}

write_mundlak_full_model_average_file <- function(record) {
  reductions <- purrr::map(
    record$prediction_manifest$prediction_path,
    read_mundlak_full_prediction_draw_file
  )
  key_reference <- reductions[[1L]]$grid
  if (!all(purrr::map_lgl(
    reductions,
    ~ identical(.x$grid, key_reference)
  ))) {
    stop("Full Mundlak prediction grids do not match within outcome.",
         call. = FALSE)
  }
  draws <- do.call(rbind, purrr::map(reductions, "draws"))
  draw_manifest <- purrr::map2_dfr(
    reductions,
    record$prediction_manifest$prediction_path,
    function(reduction, path) {
      tibble::tibble(
        model_index = reduction$model_index,
        model_name = reduction$model_name,
        prediction_file_sha256 = digest::digest(
          file = path, algo = "sha256"
        ),
        within_model_draw = seq_len(reduction$draw_count)
      )
    }
  )
  if (
    nrow(draws) != record$n_draws ||
      nrow(draw_manifest) != record$n_draws || anyNA(draws)
  ) {
    stop("The full Mundlak posterior mixture is incomplete.",
         call. = FALSE)
  }
  set.seed(record$averaging_seed)
  order <- sample(seq_len(record$n_draws))
  draws <- draws[order, , drop = FALSE]
  draw_manifest <- draw_manifest[order, , drop = FALSE] |>
    dplyr::mutate(drawid = dplyr::row_number(), .before = 1L)
  model_weights <- record$stacking_weights$weight
  names(model_weights) <- record$stacking_weights$model_name
  settings <- reductions[[1L]]$settings
  settings$n_draws <- record$n_draws
  settings$seed <- record$averaging_seed
  output <- list(
    outcome = record$outcome,
    draws = draws,
    grid = key_reference,
    draw_manifest = draw_manifest,
    model_weights = model_weights,
    settings = settings
  )
  class(output) <- c("primary_model_average", class(output))
  path <- sensitivity_reduction_file_path(
    family = "mundlak",
    outcome = record$outcome,
    model_index = 0L,
    stage = "model_average"
  )
  write_atomic_qs(output, path)
}

read_mundlak_full_model_average_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full Mundlak model-average file is unavailable.",
         call. = FALSE)
  }
  model_average <- qs::qread(path)
  if (
    !inherits(model_average, "primary_model_average") ||
      !model_average$outcome %in%
        c("whole_body", "hypertrophy", "strength", "function") ||
      nrow(model_average$draws) != 4000L ||
      ncol(model_average$draws) != nrow(model_average$grid)
  ) {
    stop("A full Mundlak model-average file is invalid.", call. = FALSE)
  }
  model_average
}

write_mundlak_full_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_mundlak_full_model_average_file(model_average_path)
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (is.null(data) || !is.data.frame(data)) {
    stop("The Mundlak reporting data are unavailable.", call. = FALSE)
  }
  absolute <- summarise_primary_dose_response(model_average)
  supported <- build_sensitivity_mundlak_reporting_pilot(
    absolute_summary = absolute,
    model_average = model_average,
    data = data,
    reference_intake = reference_intake
  )
  heatmap <- compute_primary_heatmap_contrasts(
    model_average = model_average,
    data = data,
    protein_min = reference_intake,
    sesoi = sesoi
  )
  report <- list(
    family = "mundlak",
    outcome = outcome,
    dose_response = supported$absolute_full_grid,
    dose_response_increments = supported$increments_full_grid,
    exposure_support_by_group = supported$support_by_group,
    common_exposure_support = supported$common_support,
    pairwise_contrasts = heatmap$contrasts,
    pairwise_support = heatmap$support,
    local_equivalence = summarise_primary_local_equivalence(
      heatmap$contrasts,
      outcome = outcome
    ),
    model_weights = tibble::tibble(
      model_name = names(model_average$model_weights),
      weight = as.numeric(model_average$model_weights)
    ),
    settings = utils::modifyList(
      model_average$settings,
      list(reference_intake = reference_intake, sesoi = sesoi)
    )
  )
  path <- sensitivity_reduction_file_path(
    family = "mundlak",
    outcome = outcome,
    model_index = 0L,
    stage = "reporting"
  )
  write_atomic_qs(report, path)
}

read_mundlak_full_reporting_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full Mundlak reporting file is unavailable.", call. = FALSE)
  }
  report <- qs::qread(path)
  required <- c(
    "family", "outcome", "dose_response", "dose_response_increments",
    "exposure_support_by_group", "common_exposure_support",
    "pairwise_contrasts", "pairwise_support", "local_equivalence",
    "model_weights", "settings"
  )
  if (
    !is.list(report) || !all(required %in% names(report)) ||
      !identical(report$family, "mundlak")
  ) {
    stop("A full Mundlak reporting file is invalid.", call. = FALSE)
  }
  report
}

audit_mundlak_full_reporting_files <- function(paths) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != 4L || anyNA(paths) || any(!file.exists(paths))) {
    stop("The full Mundlak reporting-file set is incomplete.",
         call. = FALSE)
  }
  summary <- purrr::map_dfr(paths, function(path) {
    report <- read_mundlak_full_reporting_file(path)
    probability_columns <- c(
      "probability_positive", "probability_negative",
      "probability_meaningful_gain", "probability_negligible",
      "probability_meaningful_loss", "probability_noninferior"
    )
    diagonal <- dplyr::filter(
      report$pairwise_contrasts,
      abs(protein_intake_current - protein_intake_proposed) < 1e-10
    )
    reference <- dplyr::filter(
      report$dose_response_increments,
      abs(protein_intake - report$settings$reference_intake) < 1e-10
    )
    allocated_counts <- allocate_primary_mixture_models(
      report$model_weights$weight,
      report$settings$n_draws
    )
    intervals_valid <- all(
      is.finite(report$dose_response$estimate) &
        is.finite(report$dose_response$conf.low) &
        is.finite(report$dose_response$conf.high) &
        report$dose_response$conf.low <= report$dose_response$estimate &
        report$dose_response$estimate <= report$dose_response$conf.high
    )
    reference_valid <-
      nrow(reference) == 2L &&
      all(abs(reference$increment) < 1e-12) &&
      all(abs(reference$increment_conf_low) < 1e-12) &&
      all(abs(reference$increment_conf_high) < 1e-12)
    probabilities <- as.matrix(
      report$pairwise_contrasts[, probability_columns, drop = FALSE]
    )
    probabilities_valid <-
      all(is.finite(probabilities)) &&
      all(probabilities >= 0 & probabilities <= 1) &&
      all(abs(
        report$pairwise_contrasts$probability_meaningful_gain +
          report$pairwise_contrasts$probability_negligible +
          report$pairwise_contrasts$probability_meaningful_loss - 1
      ) < 1e-12)
    diagonal_valid <-
      nrow(diagonal) >= 2L &&
      all(abs(diagonal$estimate) < 1e-12) &&
      all(abs(diagonal$probability_negligible - 1) < 1e-12)
    tibble::tibble(
      family = report$family,
      outcome = report$outcome,
      dose_response_rows = nrow(report$dose_response),
      increment_rows = nrow(report$dose_response_increments),
      pairwise_contrast_rows = nrow(report$pairwise_contrasts),
      local_equivalence_rows = nrow(report$local_equivalence),
      selected_models = sum(allocated_counts > 0L),
      weight_sum = sum(report$model_weights$weight),
      intervals_valid = intervals_valid,
      reference_increment_valid = reference_valid,
      probabilities_valid = probabilities_valid,
      diagonal_contrasts_valid = diagonal_valid,
      reporting_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  }) |>
    dplyr::arrange(factor(
      outcome,
      levels = c("whole_body", "hypertrophy", "strength", "function")
    ))
  if (
    nrow(summary) != 4L || anyDuplicated(summary$outcome) ||
      any(summary$dose_response_rows < 2L) ||
      any(summary$pairwise_contrast_rows < 2L) ||
      any(abs(summary$weight_sum - 1) > 1e-8) ||
      !all(summary$intervals_valid) ||
      !all(summary$reference_increment_valid) ||
      !all(summary$probabilities_valid) ||
      !all(summary$diagonal_contrasts_valid)
  ) {
    stop("The full Mundlak reporting audit failed.", call. = FALSE)
  }
  list(
    summary = summary,
    outcomes = nrow(summary),
    posterior_mixture_draws_per_outcome = 4000L,
    all_reporting_checks_pass = TRUE
  )
}

write_mundlak_full_verification_bundle <- function(
    reporting_audit,
    checkpoint,
    prediction_manifest
) {
  if (
    !isTRUE(reporting_audit$all_reporting_checks_pass) ||
      !isTRUE(checkpoint$audit$all_sampling_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$stacking_ready)
  ) {
    stop("The Mundlak family is not ready for verification.",
         call. = FALSE)
  }
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_mundlak_full_reporting_file
  )
  bind_report_table <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  analysis_data <- purrr::imap_dfr(
    checkpoint$analysis_data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  write_verification_bundle(
    family = "mundlak",
    artifacts = list(
      analysis_data = analysis_data,
      model_specifications = checkpoint$model_specifications,
      accepted_fit_manifest = checkpoint$accepted_fit_manifest,
      sampling_diagnostic_summary =
        checkpoint$sampling_diagnostic_summary,
      sampling_parameter_failures =
        checkpoint$sampling_parameter_failures,
      initial_integrated_loo_model_summary =
        checkpoint$integrated_loo_model_summary,
      initial_integrated_loo_outcome_summary =
        checkpoint$integrated_loo_outcome_summary,
      initial_high_pareto_k_cases =
        checkpoint$initial_high_pareto_k_cases,
      exact_refit_diagnostic_summary =
        checkpoint$exact_refit_diagnostic_summary,
      exact_loo_correction_summary =
        checkpoint$exact_loo_correction_summary,
      final_loo_model_summary = checkpoint$final_loo_model_summary,
      stacking_weights = checkpoint$stacking_weights,
      prediction_allocation_manifest = prediction_manifest,
      reporting_file_manifest = reporting_audit$summary,
      dose_response_predictions = bind_report_table("dose_response"),
      dose_response_increments =
        bind_report_table("dose_response_increments"),
      exposure_support_by_group =
        bind_report_table("exposure_support_by_group"),
      common_exposure_support =
        bind_report_table("common_exposure_support"),
      pairwise_contrasts = bind_report_table("pairwise_contrasts"),
      pairwise_support = bind_report_table("pairwise_support"),
      local_equivalence = bind_report_table("local_equivalence")
    ),
    descriptions = c(
      analysis_data = "Complete-case rows supplied to the Mundlak models.",
      model_specifications = "Final candidate-model registry inherited from the primary analysis.",
      accepted_fit_manifest = "Accepted one-file-per-model registry and fit checksums.",
      sampling_diagnostic_summary = "Final model-level convergence, geometry, and predictive diagnostics.",
      sampling_parameter_failures = "Any parameter-level diagnostic failures after remediation.",
      initial_integrated_loo_model_summary = "Conditional and integrated PSIS summaries before exact correction.",
      initial_integrated_loo_outcome_summary = "Outcome-level integrated PSIS audit before exact correction.",
      initial_high_pareto_k_cases = "Cases selected prospectively for exact case-deletion correction.",
      exact_refit_diagnostic_summary = "Sampling diagnostics for every exact case-deletion refit.",
      exact_loo_correction_summary = "Exact versus PSIS pointwise correction results.",
      final_loo_model_summary = "Final corrected integrated LOO summaries used for stacking.",
      stacking_weights = "Final outcome-specific stacking weights.",
      prediction_allocation_manifest = "Deterministic allocation of 4,000 mixture curves per outcome.",
      reporting_file_manifest = "Checksums and dimensions of reduced outcome reporting files.",
      dose_response_predictions = "Model-averaged absolute dose-response summaries with support flags.",
      dose_response_increments = "Paired model-averaged increments from 0.8 g/kg/day with support flags.",
      exposure_support_by_group = "Outcome- and exercise-specific observed exposure support.",
      common_exposure_support = "Exposure support shared by both exercise groups.",
      pairwise_contrasts = "Paired dose contrasts and SESOI probabilities.",
      pairwise_support = "Observed support used to flag pairwise contrasts.",
      local_equivalence = "Local practical-equivalence summaries."
    )
  )
}

read_sensitivity_completed_checkpoint <- function(path, family) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(
      "The completed ", family, " checkpoint is unavailable.",
      call. = FALSE
    )
  }
  checkpoint <- qs::qread(path)
  required <- c(
    "checkpoint_version", "family", "analysis_data_by_outcome",
    "model_specifications", "accepted_fit_manifest",
    "sampling_diagnostic_summary", "sampling_parameter_failures",
    "integrated_loo_model_summary", "integrated_loo_outcome_summary",
    "initial_high_pareto_k_cases", "exact_refit_diagnostic_summary",
    "exact_loo_correction_summary", "final_loo_model_summary",
    "stacking_weights", "audit"
  )
  if (
    !is.list(checkpoint) || !all(required %in% names(checkpoint)) ||
      !identical(checkpoint$checkpoint_version, 1L) ||
      !identical(checkpoint$family, family) ||
      !isTRUE(checkpoint$audit$all_sampling_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$stacking_ready) ||
      nrow(checkpoint$accepted_fit_manifest) != 104L ||
      nrow(checkpoint$stacking_weights) != 104L
  ) {
    stop("The completed ", family, " checkpoint is invalid.",
         call. = FALSE)
  }
  checkpoint
}

make_sensitivity_full_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    family,
    prediction_type,
    n_draws = 4000L,
    protein_min = 0.8,
    protein_max_by_outcome = NULL,
    grid_step = 0.01,
    seed = 1988L
) {
  family_prediction_valid <-
    (family %in% c("mundlak", "within") &&
       identical(family, prediction_type)) ||
    (family %in% c("achieved", "assignment") &&
       identical(prediction_type, "exposure"))
  if (
    !family_prediction_valid ||
      any(accepted_fit_manifest$family != family)
  ) {
    stop("The full sensitivity prediction family is invalid.",
         call. = FALSE)
  }
  records <- make_mundlak_full_prediction_records(
    accepted_fit_manifest = accepted_fit_manifest,
    stacking_weights = stacking_weights,
    n_draws = n_draws,
    protein_min = protein_min,
    grid_step = grid_step,
    seed = seed
  )
  purrr::map(records, function(record) {
    if (!is.null(protein_max_by_outcome)) {
      if (
        is.null(names(protein_max_by_outcome)) ||
          !record$outcome %in% names(protein_max_by_outcome) ||
          !is.finite(protein_max_by_outcome[[record$outcome]]) ||
          protein_max_by_outcome[[record$outcome]] <= protein_min
      ) {
        stop("An outcome-specific prediction maximum is invalid.",
             call. = FALSE)
      }
      record$protein_max <- as.numeric(
        protein_max_by_outcome[[record$outcome]]
      )
    }
    record$family <- family
    record$prediction_type <- prediction_type
    record
  })
}

write_sensitivity_full_prediction_draw_file <- function(
    record,
    checkpoint_path
) {
  family <- unname(as.character(record$family))
  prediction_type <- unname(as.character(record$prediction_type))
  family_prediction_valid <-
    (family %in% c("mundlak", "within") &&
       identical(family, prediction_type)) ||
    (family %in% c("achieved", "assignment") &&
       identical(prediction_type, "exposure"))
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !family_prediction_valid
  ) {
    stop("A full sensitivity prediction record is invalid.",
         call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), family) ||
      !identical(
        unname(as.character(candidate$outcome)),
        unname(as.character(record$outcome))
      ) ||
      !identical(
        unname(as.integer(candidate$model_index)),
        unname(as.integer(record$model_index))
      ) ||
      !identical(
        unname(as.character(candidate$model_name)),
        unname(as.character(record$model_name))
      )
  ) {
    stop("A sensitivity prediction fit does not match its record.",
         call. = FALSE)
  }
  checkpoint <- read_sensitivity_completed_checkpoint(
    checkpoint_path, family = family
  )
  prediction_data <- checkpoint$analysis_data_by_outcome[[record$outcome]]
  if (is.null(prediction_data) || !is.data.frame(prediction_data)) {
    stop("The original sensitivity prediction data are unavailable.",
         call. = FALSE)
  }
  grid <- make_sensitivity_prediction_grid(
    fitted_candidate = candidate,
    data = prediction_data,
    prediction_type = prediction_type,
    protein_min = record$protein_min,
    protein_max = if (is.null(record$protein_max)) {
      NULL
    } else {
      as.numeric(record$protein_max)
    },
    grid_step = record$grid_step
  )
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  set.seed(as.integer(record$prediction_seed))
  draws <- as.matrix(brms::posterior_epred(
    candidate$model,
    newdata = grid$newdata,
    re_formula = NA,
    ndraws = as.integer(record$draw_count)
  ))
  if (
    ncol(draws) != nrow(grid$keys) && nrow(draws) == nrow(grid$keys)
  ) {
    draws <- t(draws)
  }
  if (
    nrow(draws) != as.integer(record$draw_count) ||
      ncol(draws) != nrow(grid$keys) || anyNA(draws)
  ) {
    stop("Sensitivity prediction draws do not match their allocation.",
         call. = FALSE)
  }
  reduction <- list(
    family = family,
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    weight = as.numeric(record$weight),
    draw_count = as.integer(record$draw_count),
    prediction_seed = as.integer(record$prediction_seed),
    fit_sha256 = record$fit_sha256,
    grid = grid$keys,
    draws = draws,
    settings = list(
      protein_min = grid$protein_min,
      protein_max = grid$protein_max,
      grid_step = grid$grid_step,
      prediction_type = prediction_type,
      reference_intake = if (prediction_type == "within") 0.8 else NULL
    )
  )
  path <- sensitivity_reduction_file_path(
    family = family,
    outcome = record$outcome,
    model_index = record$model_index,
    stage = "prediction_draws"
  )
  write_atomic_qs(reduction, path)
}

read_sensitivity_full_prediction_draw_file <- function(path, family) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full sensitivity prediction-draw file is unavailable.",
         call. = FALSE)
  }
  reduction <- qs::qread(path)
  if (
    !is.list(reduction) || !identical(reduction$family, family) ||
      !is.matrix(reduction$draws) ||
      nrow(reduction$draws) != reduction$draw_count ||
      ncol(reduction$draws) != nrow(reduction$grid)
  ) {
    stop("A full sensitivity prediction-draw file is invalid.",
         call. = FALSE)
  }
  reduction
}

make_sensitivity_full_prediction_file_manifest <- function(
    paths,
    family,
    expected_draws = 4000L
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) < 4L || anyNA(paths) || any(!file.exists(paths))) {
    stop("The full sensitivity prediction-file set is incomplete.",
         call. = FALSE)
  }
  manifest <- purrr::map_dfr(paths, function(path) {
    reduction <- read_sensitivity_full_prediction_draw_file(path, family)
    tibble::tibble(
      family = reduction$family,
      outcome = reduction$outcome,
      model_index = reduction$model_index,
      model_name = reduction$model_name,
      weight = reduction$weight,
      draw_count = reduction$draw_count,
      prediction_seed = reduction$prediction_seed,
      fit_sha256 = reduction$fit_sha256,
      prediction_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  audit <- manifest |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      allocated_draws = sum(draw_count),
      selected_models = dplyr::n(),
      .groups = "drop"
    )
  if (
    nrow(audit) != 4L || any(audit$allocated_draws != expected_draws) ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(manifest$family != family)
  ) {
    stop("The full sensitivity prediction-file manifest is invalid.",
         call. = FALSE)
  }
  manifest
}

make_sensitivity_full_model_average_records <- function(
    prediction_manifest,
    stacking_weights,
    family,
    n_draws = 4000L,
    seed = 1988L,
    expected_models_per_outcome = 26L
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  expected_models_per_outcome <- as.integer(expected_models_per_outcome)
  if (
    !is.data.frame(prediction_manifest) ||
      !is.data.frame(stacking_weights) ||
      length(expected_models_per_outcome) != 1L ||
      is.na(expected_models_per_outcome) ||
      expected_models_per_outcome < 1L ||
      nrow(stacking_weights) !=
        length(outcomes) * expected_models_per_outcome ||
      any(prediction_manifest$family != family)
  ) {
    stop("The full sensitivity averaging inputs are invalid.",
         call. = FALSE)
  }
  purrr::map(seq_along(outcomes), function(index) {
    outcome <- outcomes[[index]]
    files <- dplyr::filter(
      prediction_manifest, .data$outcome == .env$outcome
    )
    weights <- dplyr::filter(
      stacking_weights, .data$outcome == .env$outcome
    ) |>
      dplyr::arrange(model_index)
    if (
      nrow(files) < 1L || sum(files$draw_count) != n_draws ||
        nrow(weights) != expected_models_per_outcome ||
        abs(sum(weights$weight) - 1) > 1e-8 ||
        any(!file.exists(files$prediction_path))
    ) {
      stop("An outcome-specific sensitivity model-average record is invalid.",
           call. = FALSE)
    }
    list(
      family = family,
      outcome = outcome,
      prediction_manifest = files,
      stacking_weights = weights,
      n_draws = as.integer(n_draws),
      averaging_seed = as.integer(seed + index)
    )
  })
}

write_sensitivity_full_model_average_file <- function(
    record,
    result_path = NULL
) {
  family <- record$family
  reductions <- purrr::map(
    record$prediction_manifest$prediction_path,
    read_sensitivity_full_prediction_draw_file,
    family = family
  )
  key_reference <- reductions[[1L]]$grid
  if (!all(purrr::map_lgl(
    reductions, ~ identical(.x$grid, key_reference)
  ))) {
    stop("Full sensitivity prediction grids do not match within outcome.",
         call. = FALSE)
  }
  draws <- do.call(rbind, purrr::map(reductions, "draws"))
  draw_manifest <- purrr::map2_dfr(
    reductions,
    record$prediction_manifest$prediction_path,
    function(reduction, path) {
      tibble::tibble(
        model_index = reduction$model_index,
        model_name = reduction$model_name,
        prediction_file_sha256 = digest::digest(
          file = path, algo = "sha256"
        ),
        within_model_draw = seq_len(reduction$draw_count)
      )
    }
  )
  if (
    nrow(draws) != record$n_draws ||
      nrow(draw_manifest) != record$n_draws || anyNA(draws)
  ) {
    stop("The full sensitivity posterior mixture is incomplete.",
         call. = FALSE)
  }
  set.seed(record$averaging_seed)
  order <- sample(seq_len(record$n_draws))
  draws <- draws[order, , drop = FALSE]
  draw_manifest <- draw_manifest[order, , drop = FALSE] |>
    dplyr::mutate(drawid = dplyr::row_number(), .before = 1L)
  model_weights <- record$stacking_weights$weight
  names(model_weights) <- record$stacking_weights$model_name
  settings <- reductions[[1L]]$settings
  settings$n_draws <- record$n_draws
  settings$seed <- record$averaging_seed
  output <- list(
    outcome = record$outcome,
    draws = draws,
    grid = key_reference,
    draw_manifest = draw_manifest,
    model_weights = model_weights,
    settings = settings
  )
  class(output) <- c("primary_model_average", class(output))
  path <- if (is.null(result_path)) {
    sensitivity_reduction_file_path(
      family = family,
      outcome = record$outcome,
      model_index = 0L,
      stage = "model_average"
    )
  } else {
    dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
    result_path
  }
  write_atomic_qs(output, path)
}

read_sensitivity_full_model_average_file <- function(path, family) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full sensitivity model-average file is unavailable.",
         call. = FALSE)
  }
  model_average <- qs::qread(path)
  expected_prediction_type <- if (
    family %in% c("achieved", "assignment")
  ) "exposure" else family
  if (
    !inherits(model_average, "primary_model_average") ||
      !model_average$outcome %in%
        c("whole_body", "hypertrophy", "strength", "function") ||
      nrow(model_average$draws) != 4000L ||
      ncol(model_average$draws) != nrow(model_average$grid) ||
      !identical(
        model_average$settings$prediction_type,
        expected_prediction_type
      )
  ) {
    stop("A full sensitivity model-average file is invalid.",
         call. = FALSE)
  }
  model_average
}

write_sensitivity_full_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    family,
    reference_intake = 0.8,
    sesoi = 0.10,
    result_path = NULL
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path, family
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (is.null(data) || !is.data.frame(data)) {
    stop("The sensitivity reporting data are unavailable.", call. = FALSE)
  }
  absolute <- summarise_primary_dose_response(model_average)
  supported <- switch(
    family,
    mundlak = build_sensitivity_mundlak_reporting_pilot(
      absolute, model_average, data, reference_intake
    ),
    within = build_sensitivity_within_reporting_pilot(
      absolute, model_average, data, reference_intake
    ),
    achieved = build_sensitivity_mundlak_reporting_pilot(
      absolute, model_average, data, reference_intake
    ),
    assignment = build_sensitivity_mundlak_reporting_pilot(
      absolute, model_average, data, reference_intake
    ),
    eiv = build_sensitivity_mundlak_reporting_pilot(
      absolute, model_average, data, reference_intake
    ),
    rescor = build_sensitivity_mundlak_reporting_pilot(
      absolute, model_average, data, reference_intake
    ),
    stop("The sensitivity reporting family is unsupported.", call. = FALSE)
  )
  heatmap_data <- if (family == "within") {
    dplyr::mutate(
      data,
      protein_intake_raw = protein_intake_raw_treatment
    )
  } else {
    data
  }
  heatmap <- compute_primary_heatmap_contrasts(
    model_average = model_average,
    data = heatmap_data,
    protein_min = reference_intake,
    sesoi = sesoi
  )
  report <- list(
    family = family,
    outcome = outcome,
    dose_response = supported$absolute_full_grid,
    dose_response_increments = supported$increments_full_grid,
    exposure_support_by_group = supported$support_by_group,
    common_exposure_support = supported$common_support,
    pairwise_contrasts = heatmap$contrasts,
    pairwise_support = heatmap$support,
    local_equivalence = summarise_primary_local_equivalence(
      heatmap$contrasts, outcome = outcome
    ),
    model_weights = tibble::tibble(
      model_name = names(model_average$model_weights),
      weight = as.numeric(model_average$model_weights)
    ),
    settings = utils::modifyList(
      model_average$settings,
      list(reference_intake = reference_intake, sesoi = sesoi)
    )
  )
  path <- if (is.null(result_path)) {
    sensitivity_reduction_file_path(
      family = family,
      outcome = outcome,
      model_index = 0L,
      stage = "reporting"
    )
  } else {
    dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
    result_path
  }
  write_atomic_qs(report, path)
}

read_sensitivity_full_reporting_file <- function(path, family) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A full sensitivity reporting file is unavailable.",
         call. = FALSE)
  }
  report <- qs::qread(path)
  required <- c(
    "family", "outcome", "dose_response", "dose_response_increments",
    "exposure_support_by_group", "common_exposure_support",
    "pairwise_contrasts", "pairwise_support", "local_equivalence",
    "model_weights", "settings"
  )
  if (
    !is.list(report) || !all(required %in% names(report)) ||
      !identical(report$family, family)
  ) {
    stop("A full sensitivity reporting file is invalid.", call. = FALSE)
  }
  report
}

audit_sensitivity_full_reporting_files <- function(paths, family) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != 4L || anyNA(paths) || any(!file.exists(paths))) {
    stop("The full sensitivity reporting-file set is incomplete.",
         call. = FALSE)
  }
  summary <- purrr::map_dfr(paths, function(path) {
    report <- read_sensitivity_full_reporting_file(path, family)
    probability_columns <- c(
      "probability_positive", "probability_negative",
      "probability_meaningful_gain", "probability_negligible",
      "probability_meaningful_loss", "probability_noninferior"
    )
    diagonal <- dplyr::filter(
      report$pairwise_contrasts,
      abs(protein_intake_current - protein_intake_proposed) < 1e-10
    )
    reference <- dplyr::filter(
      report$dose_response_increments,
      abs(protein_intake - report$settings$reference_intake) < 1e-10
    )
    allocated_counts <- allocate_primary_mixture_models(
      report$model_weights$weight,
      report$settings$n_draws
    )
    intervals_valid <- all(
      is.finite(report$dose_response$estimate) &
        is.finite(report$dose_response$conf.low) &
        is.finite(report$dose_response$conf.high) &
        report$dose_response$conf.low <= report$dose_response$estimate &
        report$dose_response$estimate <= report$dose_response$conf.high
    )
    reference_valid <-
      nrow(reference) == 2L &&
      all(abs(reference$increment) < 1e-12) &&
      all(abs(reference$increment_conf_low) < 1e-12) &&
      all(abs(reference$increment_conf_high) < 1e-12)
    probabilities <- as.matrix(
      report$pairwise_contrasts[, probability_columns, drop = FALSE]
    )
    probabilities_valid <-
      all(is.finite(probabilities)) &&
      all(probabilities >= 0 & probabilities <= 1) &&
      all(abs(
        report$pairwise_contrasts$probability_meaningful_gain +
          report$pairwise_contrasts$probability_negligible +
          report$pairwise_contrasts$probability_meaningful_loss - 1
      ) < 1e-12)
    diagonal_valid <-
      nrow(diagonal) >= 2L &&
      all(abs(diagonal$estimate) < 1e-12) &&
      all(abs(diagonal$probability_negligible - 1) < 1e-12)
    tibble::tibble(
      family = report$family,
      outcome = report$outcome,
      dose_response_rows = nrow(report$dose_response),
      increment_rows = nrow(report$dose_response_increments),
      pairwise_contrast_rows = nrow(report$pairwise_contrasts),
      local_equivalence_rows = nrow(report$local_equivalence),
      selected_models = sum(allocated_counts > 0L),
      weight_sum = sum(report$model_weights$weight),
      intervals_valid = intervals_valid,
      reference_increment_valid = reference_valid,
      probabilities_valid = probabilities_valid,
      diagonal_contrasts_valid = diagonal_valid,
      reporting_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  }) |>
    dplyr::arrange(factor(
      outcome,
      levels = c("whole_body", "hypertrophy", "strength", "function")
    ))
  if (
    nrow(summary) != 4L || anyDuplicated(summary$outcome) ||
      any(summary$dose_response_rows < 2L) ||
      any(summary$pairwise_contrast_rows < 2L) ||
      any(abs(summary$weight_sum - 1) > 1e-8) ||
      !all(summary$intervals_valid) ||
      !all(summary$reference_increment_valid) ||
      !all(summary$probabilities_valid) ||
      !all(summary$diagonal_contrasts_valid)
  ) {
    stop("The full sensitivity reporting audit failed.", call. = FALSE)
  }
  list(
    family = family,
    summary = summary,
    outcomes = nrow(summary),
    posterior_mixture_draws_per_outcome = 4000L,
    all_reporting_checks_pass = TRUE
  )
}

make_exposure_full_comparison_records <- function(
    achieved_model_average_paths,
    assignment_model_average_paths,
    achieved_reporting_audit,
    assignment_reporting_audit,
    seed = 23988L
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  make_lookup <- function(paths, family) {
    paths <- unname(unlist(paths, use.names = FALSE))
    if (length(paths) != 4L || any(!file.exists(paths))) {
      stop("An exposure model-average file set is incomplete.",
           call. = FALSE)
    }
    objects <- purrr::map(
      paths,
      read_sensitivity_full_model_average_file,
      family = family
    )
    tibble::tibble(
      outcome = purrr::map_chr(objects, "outcome"),
      model_average_path = normalizePath(
        paths, winslash = "/", mustWork = TRUE
      )
    )
  }
  achieved <- make_lookup(achieved_model_average_paths, "achieved") |>
    dplyr::rename(achieved_model_average_path = model_average_path)
  assignment <- make_lookup(assignment_model_average_paths, "assignment") |>
    dplyr::rename(assignment_model_average_path = model_average_path)
  achieved_reports <- achieved_reporting_audit$summary |>
    dplyr::select(outcome, achieved_reporting_path = reporting_path)
  assignment_reports <- assignment_reporting_audit$summary |>
    dplyr::select(outcome, assignment_reporting_path = reporting_path)
  records <- tibble::tibble(outcome = outcomes) |>
    dplyr::left_join(achieved, by = "outcome", relationship = "one-to-one") |>
    dplyr::left_join(
      assignment, by = "outcome", relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      achieved_reports, by = "outcome", relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      assignment_reports, by = "outcome", relationship = "one-to-one"
    ) |>
    dplyr::mutate(seed = as.integer(seed + dplyr::row_number()))
  path_columns <- grep("_path$", names(records), value = TRUE)
  if (
    nrow(records) != 4L || anyDuplicated(records$outcome) ||
      anyNA(records[, path_columns, drop = FALSE]) ||
      any(!vapply(
        records[, path_columns, drop = FALSE],
        function(x) all(file.exists(x)),
        logical(1)
      ))
  ) {
    stop("The exposure comparison records are invalid.", call. = FALSE)
  }
  purrr::map(seq_len(nrow(records)), function(row) {
    as.list(records[row, , drop = FALSE])
  })
}

make_exposure_increment_draws <- function(draws, grid, reference_intake) {
  output <- draws
  groups <- unique(grid$resistance_exercise_code)
  for (group in groups) {
    columns <- which(grid$resistance_exercise_code == group)
    reference <- columns[[which.min(abs(
      grid$protein_intake[columns] - reference_intake
    ))]]
    if (abs(grid$protein_intake[[reference]] - reference_intake) > 1e-10) {
      stop("The exposure comparison reference is absent.", call. = FALSE)
    }
    output[, columns] <- sweep(
      draws[, columns, drop = FALSE],
      1L,
      draws[, reference],
      "-"
    )
  }
  output
}

summarise_exposure_difference_draws <- function(
    draws,
    grid,
    sesoi = 0.10
) {
  model_average <- list(draws = draws, grid = grid)
  class(model_average) <- c("primary_model_average", class(model_average))
  summary <- summarise_primary_dose_response(model_average)
  summary$probability_assignment_greater <- colMeans(draws > 0)
  summary$probability_assignment_lower <- colMeans(draws < 0)
  summary$probability_meaningfully_greater <- colMeans(draws > sesoi)
  summary$probability_negligible_difference <- colMeans(
    abs(draws) <= sesoi
  )
  summary$probability_meaningfully_lower <- colMeans(draws < -sesoi)
  summary
}

write_exposure_full_comparison_file <- function(
    record,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  achieved <- read_sensitivity_full_model_average_file(
    record$achieved_model_average_path, "achieved"
  )
  assignment <- read_sensitivity_full_model_average_file(
    record$assignment_model_average_path, "assignment"
  )
  achieved_report <- read_sensitivity_full_reporting_file(
    record$achieved_reporting_path, "achieved"
  )
  assignment_report <- read_sensitivity_full_reporting_file(
    record$assignment_reporting_path, "assignment"
  )
  if (
    !identical(achieved$outcome, record$outcome) ||
      !identical(assignment$outcome, record$outcome) ||
      !identical(achieved$grid, assignment$grid) ||
      !identical(dim(achieved$draws), dim(assignment$draws)) ||
      nrow(achieved$draws) != 4000L
  ) {
    stop("The paired exposure model averages are misaligned.",
         call. = FALSE)
  }
  set.seed(as.integer(record$seed))
  assignment_order <- sample(seq_len(nrow(assignment$draws)))
  achieved_absolute <- achieved$draws
  assignment_absolute <- assignment$draws[assignment_order, , drop = FALSE]
  achieved_increment <- make_exposure_increment_draws(
    achieved_absolute, achieved$grid, reference_intake
  )
  assignment_increment <- make_exposure_increment_draws(
    assignment_absolute, assignment$grid, reference_intake
  )
  output <- list(
    family = "exposure_comparison",
    outcome = record$outcome,
    achieved_dose_response = achieved_report$dose_response,
    assignment_dose_response = assignment_report$dose_response,
    assignment_minus_achieved = summarise_exposure_difference_draws(
      assignment_absolute - achieved_absolute,
      achieved$grid,
      sesoi = sesoi
    ),
    achieved_increments = achieved_report$dose_response_increments,
    assignment_increments = assignment_report$dose_response_increments,
    assignment_minus_achieved_increment =
      summarise_exposure_difference_draws(
        assignment_increment - achieved_increment,
        achieved$grid,
        sesoi = sesoi
      ),
    settings = list(
      n_draws = 4000L,
      seed = as.integer(record$seed),
      reference_intake = reference_intake,
      sesoi = sesoi,
      pairing = "independently_shuffled_posterior_mixture_draws"
    )
  )
  path <- sensitivity_reduction_file_path(
    family = "exposure_comparison",
    outcome = record$outcome,
    model_index = 0L,
    stage = "comparison"
  )
  write_atomic_qs(output, path)
}

read_exposure_full_comparison_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("An exposure comparison file is unavailable.", call. = FALSE)
  }
  comparison <- qs::qread(path)
  required <- c(
    "family", "outcome", "achieved_dose_response",
    "assignment_dose_response", "assignment_minus_achieved",
    "achieved_increments", "assignment_increments",
    "assignment_minus_achieved_increment", "settings"
  )
  if (
    !is.list(comparison) || !all(required %in% names(comparison)) ||
      !identical(comparison$family, "exposure_comparison")
  ) {
    stop("An exposure comparison file is invalid.", call. = FALSE)
  }
  comparison
}

audit_exposure_full_comparison_files <- function(paths) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != 4L || any(!file.exists(paths))) {
    stop("The exposure comparison file set is incomplete.", call. = FALSE)
  }
  summary <- purrr::map_dfr(paths, function(path) {
    x <- read_exposure_full_comparison_file(path)
    probability_columns <- c(
      "probability_assignment_greater", "probability_assignment_lower",
      "probability_meaningfully_greater",
      "probability_negligible_difference",
      "probability_meaningfully_lower"
    )
    probabilities <- rbind(
      as.matrix(x$assignment_minus_achieved[, probability_columns]),
      as.matrix(
        x$assignment_minus_achieved_increment[, probability_columns]
      )
    )
    aligned_rows <- length(unique(c(
      nrow(x$achieved_dose_response),
      nrow(x$assignment_dose_response),
      nrow(x$assignment_minus_achieved),
      nrow(x$achieved_increments),
      nrow(x$assignment_increments),
      nrow(x$assignment_minus_achieved_increment)
    ))) == 1L
    tibble::tibble(
      family = x$family,
      outcome = x$outcome,
      grid_rows = nrow(x$assignment_minus_achieved),
      aligned_rows = aligned_rows,
      intervals_valid = all(
        x$assignment_minus_achieved$conf.low <=
          x$assignment_minus_achieved$estimate &
          x$assignment_minus_achieved$estimate <=
            x$assignment_minus_achieved$conf.high &
          x$assignment_minus_achieved_increment$conf.low <=
            x$assignment_minus_achieved_increment$estimate &
          x$assignment_minus_achieved_increment$estimate <=
            x$assignment_minus_achieved_increment$conf.high
      ),
      probabilities_valid = all(
        is.finite(probabilities) & probabilities >= 0 & probabilities <= 1
      ),
      comparison_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256")
    )
  }) |>
    dplyr::arrange(factor(
      outcome,
      levels = c("whole_body", "hypertrophy", "strength", "function")
    ))
  if (
    nrow(summary) != 4L || anyDuplicated(summary$outcome) ||
      !all(summary$aligned_rows) || !all(summary$intervals_valid) ||
      !all(summary$probabilities_valid)
  ) {
    stop("The exposure comparison audit failed.", call. = FALSE)
  }
  list(
    family = "exposure_comparison",
    summary = summary,
    all_comparison_checks_pass = TRUE
  )
}

write_exposure_full_comparison_verification_bundle <- function(
    comparison_audit
) {
  if (!isTRUE(comparison_audit$all_comparison_checks_pass)) {
    stop("The exposure comparison is not ready for verification.",
         call. = FALSE)
  }
  comparisons <- purrr::map(
    comparison_audit$summary$comparison_path,
    read_exposure_full_comparison_file
  )
  bind_component <- function(name) {
    purrr::map_dfr(comparisons, function(x) {
      dplyr::mutate(x[[name]], outcome = x$outcome, .before = 1L)
    })
  }
  write_verification_bundle(
    family = "exposure_comparison",
    artifacts = list(
      comparison_file_manifest = comparison_audit$summary,
      achieved_dose_response = bind_component("achieved_dose_response"),
      assignment_dose_response =
        bind_component("assignment_dose_response"),
      assignment_minus_achieved =
        bind_component("assignment_minus_achieved"),
      achieved_increments = bind_component("achieved_increments"),
      assignment_increments = bind_component("assignment_increments"),
      assignment_minus_achieved_increment =
        bind_component("assignment_minus_achieved_increment")
    ),
    descriptions = c(
      comparison_file_manifest = "Checksums and dimensions of the outcome-specific exposure comparison files.",
      achieved_dose_response = "Model-averaged curves indexed by achieved protein intake.",
      assignment_dose_response = "Model-averaged curves indexed by baseline plus assignment-derived protein intake.",
      assignment_minus_achieved = "Posterior assignment-minus-achieved differences in absolute dose-response curves.",
      achieved_increments = "Achieved-intake increments from 0.8 g/kg/day.",
      assignment_increments = "Assignment-derived increments from 0.8 g/kg/day.",
      assignment_minus_achieved_increment = "Posterior assignment-minus-achieved differences in increments from 0.8 g/kg/day."
    )
  )
}

write_sensitivity_full_verification_bundle <- function(
    reporting_audit,
    checkpoint,
    prediction_manifest,
    family
) {
  if (
    !identical(checkpoint$family, family) ||
      !isTRUE(reporting_audit$all_reporting_checks_pass) ||
      !isTRUE(checkpoint$audit$all_sampling_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$stacking_ready)
  ) {
    stop("The sensitivity family is not ready for verification.",
         call. = FALSE)
  }
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = family
  )
  bind_report_table <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  analysis_data <- purrr::imap_dfr(
    checkpoint$analysis_data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  analysis_description <- switch(
    family,
    within = "Strict matched-arm contrasts supplied to the within-study models.",
    achieved = "Common-sample rows using observed achieved protein intake.",
    assignment = "Common-sample rows using baseline plus assignment-derived protein intake.",
    eiv = "Complete-case rows supplied to the Stan errors-in-variables models, including the arm-level protein-intake standard errors.",
    rescor = "Complete-case rows supplied to the Stan models with estimated within-study sampling-error correlation.",
    mundlak = "Complete-case rows supplied to the Mundlak models.",
    "Outcome-specific rows supplied to the sensitivity models."
  )
  curve_description <- switch(
    family,
    within = "Model-averaged protein-arm minus comparator-arm curves with comparator fixed at 0.8 g/kg/day.",
    achieved = "Model-averaged dose-response curves indexed by observed achieved protein intake.",
    assignment = "Model-averaged dose-response curves indexed by baseline plus assignment-derived protein intake.",
    eiv = "Model-averaged dose-response curves integrating arm-level protein-intake measurement error.",
    rescor = "Model-averaged dose-response curves allowing estimated within-study sampling-error correlation.",
    "Model-averaged sensitivity dose-response curves."
  )
  increment_description <- switch(
    family,
    within = "Paired model-averaged increments from the 0.8 g/kg/day treatment reference.",
    achieved = "Paired achieved-intake model-average increments from 0.8 g/kg/day.",
    assignment = "Paired assignment-derived model-average increments from 0.8 g/kg/day.",
    eiv = "Paired errors-in-variables model-average increments from 0.8 g/kg/day.",
    rescor = "Paired residual-correlation model-average increments from 0.8 g/kg/day.",
    "Paired model-average increments from 0.8 g/kg/day."
  )
  write_verification_bundle(
    family = family,
    artifacts = list(
      analysis_data = analysis_data,
      model_specifications = checkpoint$model_specifications,
      accepted_fit_manifest = checkpoint$accepted_fit_manifest,
      sampling_diagnostic_summary = checkpoint$sampling_diagnostic_summary,
      sampling_parameter_failures = checkpoint$sampling_parameter_failures,
      initial_integrated_loo_model_summary =
        checkpoint$integrated_loo_model_summary,
      initial_integrated_loo_outcome_summary =
        checkpoint$integrated_loo_outcome_summary,
      initial_high_pareto_k_cases = checkpoint$initial_high_pareto_k_cases,
      exact_refit_diagnostic_summary =
        checkpoint$exact_refit_diagnostic_summary,
      exact_loo_correction_summary = checkpoint$exact_loo_correction_summary,
      final_loo_model_summary = checkpoint$final_loo_model_summary,
      stacking_weights = checkpoint$stacking_weights,
      prediction_allocation_manifest = prediction_manifest,
      reporting_file_manifest = reporting_audit$summary,
      dose_response_predictions = bind_report_table("dose_response"),
      dose_response_increments =
        bind_report_table("dose_response_increments"),
      exposure_support_by_group =
        bind_report_table("exposure_support_by_group"),
      common_exposure_support =
        bind_report_table("common_exposure_support"),
      pairwise_contrasts = bind_report_table("pairwise_contrasts"),
      pairwise_support = bind_report_table("pairwise_support"),
      local_equivalence = bind_report_table("local_equivalence")
    ),
    descriptions = c(
      analysis_data = analysis_description,
      model_specifications = "Final candidate-model registry inherited from the primary analysis.",
      accepted_fit_manifest = "Accepted one-file-per-model registry and fit checksums.",
      sampling_diagnostic_summary = "Final model-level convergence, geometry, and predictive diagnostics.",
      sampling_parameter_failures = "Any parameter-level diagnostic failures after remediation.",
      initial_integrated_loo_model_summary = "Conditional and integrated PSIS summaries before exact correction.",
      initial_integrated_loo_outcome_summary = "Outcome-level integrated PSIS audit before exact correction.",
      initial_high_pareto_k_cases = "Cases selected prospectively for exact case-deletion correction.",
      exact_refit_diagnostic_summary = "Sampling diagnostics for every exact case-deletion refit.",
      exact_loo_correction_summary = "Exact versus PSIS pointwise correction results.",
      final_loo_model_summary = "Final corrected integrated LOO summaries used for stacking.",
      stacking_weights = "Final outcome-specific stacking weights.",
      prediction_allocation_manifest = "Deterministic allocation of 4,000 mixture curves per outcome.",
      reporting_file_manifest = "Checksums and dimensions of reduced outcome reporting files.",
      dose_response_predictions = curve_description,
      dose_response_increments = increment_description,
      exposure_support_by_group = "Outcome- and exercise-specific exposure support.",
      common_exposure_support = "Exposure support shared by both exercise groups.",
      pairwise_contrasts = "Paired dose contrasts and SESOI probabilities.",
      pairwise_support = "Observed treatment-intake support used to flag pairwise contrasts.",
      local_equivalence = "Local practical-equivalence summaries."
    )
  )
}

read_within_pilot_completed_checkpoint <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(
      "The completed within-pilot checkpoint is unavailable. Run ",
      "sensitivity/freeze_within_pilot_checkpoint.R first.",
      call. = FALSE
    )
  }
  checkpoint <- qs::qread(path)
  required <- c(
    "checkpoint_version", "family", "analysis_data_by_outcome",
    "priors_by_outcome", "model_specifications",
    "reused_fit_manifest", "pilot_diagnostic_summary",
    "pilot_parameter_failures", "audit"
  )
  if (
    !is.list(checkpoint) || !all(required %in% names(checkpoint)) ||
      !identical(checkpoint$checkpoint_version, 1L) ||
      !identical(checkpoint$family, "within") ||
      !isTRUE(checkpoint$audit$all_pilot_diagnostics_pass) ||
      nrow(checkpoint$reused_fit_manifest) != 7L ||
      any(!file.exists(checkpoint$reused_fit_manifest$fit_path))
  ) {
    stop("The completed within-pilot checkpoint is invalid.",
         call. = FALSE)
  }
  checkpoint
}

make_within_full_contrast_data_by_outcome <- function(
    analysis_data_by_outcome
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.list(analysis_data_by_outcome) ||
      !setequal(names(analysis_data_by_outcome), outcomes)
  ) {
    stop("The within-contrast source data are incomplete.",
         call. = FALSE)
  }
  contrasts <- purrr::map(
    analysis_data_by_outcome[outcomes],
    make_strict_within_study_contrasts
  )
  valid <- purrr::map_lgl(contrasts, function(data) {
    is.data.frame(data) && nrow(data) > 0L &&
      all(c(
        "yi", "vi", "study_no", "arm_no", "outcome_no",
        "protein_intake_raw_treatment",
        "protein_intake_raw_comparator"
      ) %in% names(data)) &&
      all(is.finite(data$yi)) && all(is.finite(data$vi)) &&
      all(data$vi > 0)
  })
  if (!all(valid)) {
    stop("An outcome-specific within-contrast data set is invalid.",
         call. = FALSE)
  }
  contrasts
}

select_pending_within_full_specifications <- function(
    specifications,
    reused_fit_manifest
) {
  specification_keys <- paste(
    purrr::map_chr(specifications, "outcome"),
    purrr::map_int(specifications, "model_index")
  )
  reused_keys <- paste(
    reused_fit_manifest$outcome,
    reused_fit_manifest$model_index
  )
  expected_reused <- paste(
    "whole_body", c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
  )
  if (
    length(specifications) != 104L || anyDuplicated(specification_keys) ||
      nrow(reused_fit_manifest) != 7L ||
      !setequal(reused_keys, expected_reused) ||
      !all(reused_keys %in% specification_keys) ||
      any(!file.exists(reused_fit_manifest$fit_path))
  ) {
    stop("The full within-contrast specification/reuse set is invalid.",
         call. = FALSE)
  }
  specifications[!specification_keys %in% reused_keys]
}

prepare_within_full_input <- function(
    specification,
    contrast_data_by_outcome,
    priors_by_outcome,
    config
) {
  outcome <- specification$outcome
  data <- contrast_data_by_outcome[[outcome]]
  priors <- priors_by_outcome[[outcome]]
  if (is.null(data) || is.null(priors) || !outcome %in% config$outcomes) {
    stop("A full within-contrast outcome input is unavailable.",
         call. = FALSE)
  }
  specification$sensitivity_within_expansion <- TRUE
  specification$sensitivity_within_sampling <- paste0(
    config$iter_warmup, "_warmup_", config$iter_sampling, "_sampling"
  )
  specification$sensitivity_within_pilot_informed_remedy <-
    "longer_sampling"
  prepare_within_contrast_pilot_input(data, specification, priors)
}

audit_within_full_expansion <- function(
    pending_specifications,
    inputs,
    contrast_data_by_outcome,
    reused_fit_manifest,
    config
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  pending_keys <- paste(
    purrr::map_chr(pending_specifications, "outcome"),
    purrr::map_int(pending_specifications, "model_index")
  )
  input_keys <- paste(
    purrr::map_chr(inputs, ~ .x$specification$outcome),
    purrr::map_int(inputs, ~ .x$specification$model_index)
  )
  observed_counts <- table(factor(
    purrr::map_chr(inputs, ~ .x$specification$outcome),
    levels = outcomes
  ))
  expected_counts <- c(
    whole_body = 19L,
    hypertrophy = 26L,
    strength = 26L,
    `function` = 26L
  )
  valid_inputs <- purrr::map_lgl(inputs, function(input) {
    identical(input$family, "within") &&
      is.character(input$formula) && length(input$formula) == 1L &&
      is.data.frame(input$data) && nrow(input$data) > 0L
  })
  contrast_summary <- purrr::imap_dfr(
    contrast_data_by_outcome,
    function(data, outcome) {
      tibble::tibble(
        outcome = outcome,
        contrasts = nrow(data),
        studies = dplyr::n_distinct(data$study_no),
        matched_arm_pairs = dplyr::n_distinct(data$arm_no)
      )
    }
  )
  if (
    length(pending_specifications) != 97L || length(inputs) != 97L ||
      nrow(reused_fit_manifest) != 7L || anyDuplicated(pending_keys) ||
      anyDuplicated(input_keys) || !setequal(pending_keys, input_keys) ||
      !identical(as.integer(observed_counts), unname(expected_counts)) ||
      !all(valid_inputs) || nrow(contrast_summary) != 4L ||
      any(contrast_summary$contrasts < 1L) ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("The full within-contrast expansion failed its structural audit.",
         call. = FALSE)
  }
  list(
    checks = tibble::tibble(
      check = c(
        "complete candidate registry",
        "reused passing whole-body pilots",
        "new independently checkpointed fits",
        "pilot-informed sampling schedule",
        "strict matched-arm contrast data"
      ),
      status = "PASS",
      detail = c(
        "104 candidates: 26 models x 4 outcomes",
        "7 fits: models 1, 2, 3, 4, 5, 13, and 26",
        "97 fits: 19 whole-body plus 26 for each other outcome",
        "1,500 warmup and 2,000 retained iterations per chain",
        "protein arm minus matched comparator arm within study"
      )
    ),
    contrast_summary = contrast_summary,
    pending_fits = length(inputs),
    reused_fits = nrow(reused_fit_manifest),
    structure_pass = TRUE
  )
}

fit_within_full_candidate_to_file <- function(
    input,
    enabled,
    config,
    stage = "initial"
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("within"))
  outcome_offset <- match(input$specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the full within-contrast fit.",
         call. = FALSE)
  }
  local_config <- config
  local_config$seed <- as.integer(
    config$seed + 100L * (outcome_offset - 1L) +
      input$specification$model_index
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "within"
  )
  write_sensitivity_fit_file(candidate, stage = stage)
}

select_within_full_remedy_specifications <- function(
    specifications,
    initial_audit
) {
  failed <- initial_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::transmute(
      outcome = as.character(outcome),
      model_index = as.integer(model_index),
      key = paste(outcome, model_index)
    )
  specification_keys <- paste(
    purrr::map_chr(specifications, "outcome"),
    purrr::map_int(specifications, "model_index")
  )
  selected <- specifications[specification_keys %in% failed$key]
  selected_keys <- paste(
    purrr::map_chr(selected, "outcome"),
    purrr::map_int(selected, "model_index")
  )
  if (
    nrow(failed) < 1L || anyDuplicated(failed$key) ||
      anyDuplicated(specification_keys) ||
      length(selected) != nrow(failed) ||
      !setequal(selected_keys, failed$key)
  ) {
    stop("The full within-contrast remedy specification set is invalid.",
         call. = FALSE)
  }
  selected
}

prepare_within_full_remedy_input <- function(
    specification,
    contrast_data_by_outcome,
    priors_by_outcome,
    config
) {
  specification$adapt_delta <- max(
    specification$adapt_delta,
    config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth,
    config$max_treedepth_floor
  )
  specification$sensitivity_within_full_remedy <-
    "longer_sampling_high_adapt_delta"
  prepare_within_full_input(
    specification,
    contrast_data_by_outcome,
    priors_by_outcome,
    config
  )
}

combine_within_full_remedy_diagnostic_audit <- function(
    initial_audit,
    remedy_diagnostics
) {
  remedy_summary <- dplyr::bind_rows(purrr::map(
    remedy_diagnostics, "summary"
  ))
  remedy_keys <- paste(remedy_summary$outcome, remedy_summary$model_index)
  expected_keys <- initial_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::transmute(key = paste(outcome, model_index)) |>
    dplyr::pull(key)
  if (
    nrow(remedy_summary) != length(expected_keys) ||
      anyDuplicated(remedy_keys) || !setequal(remedy_keys, expected_keys)
  ) {
    stop("The full within-contrast remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  summary <- initial_audit$diagnostic_summary |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_summary) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    ) |>
    dplyr::mutate(
      ppc_warning = primary_ppc_warning(dplyr::pick(dplyr::everything())),
      status = dplyr::case_when(
        !diagnostic_pass ~ "REQUIRES REMEDIATION",
        ppc_warning ~ "PASS WITH PPC WARNING",
        TRUE ~ "PASS"
      )
    )
  remedy_parameter_failures <- purrr::map_dfr(
    remedy_diagnostics,
    function(result) {
      result$parameter_failures |>
        dplyr::mutate(
          outcome = result$summary$outcome[[1L]],
          model_index = result$summary$model_index[[1L]],
          .before = 1L
        )
    }
  )
  parameter_failures <- initial_audit$parameter_failures |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_parameter_failures)
  if (
    nrow(summary) != 104L ||
      anyDuplicated(paste(summary$outcome, summary$model_index))
  ) {
    stop("The remedied within-contrast diagnostic audit is incomplete.",
         call. = FALSE)
  }
  list(
    family = "within",
    diagnostic_summary = summary,
    parameter_failures = parameter_failures,
    remedy_diagnostic_summary = remedy_summary,
    remedy_keys = remedy_keys,
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

combine_within_full_accepted_fit_manifest <- function(
    initial_manifest,
    remedy_manifest,
    final_diagnostic_audit
) {
  if (!isTRUE(final_diagnostic_audit$all_sampling_diagnostics_pass)) {
    stop("Within-contrast fits cannot be accepted before diagnostics pass.",
         call. = FALSE)
  }
  remedy_keys <- paste(remedy_manifest$outcome, remedy_manifest$model_index)
  accepted <- initial_manifest |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_manifest) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  if (
    nrow(accepted) != 104L ||
      anyDuplicated(paste(accepted$outcome, accepted$model_index)) ||
      any(accepted$family != "within") ||
      any(!file.exists(accepted$fit_path))
  ) {
    stop("The accepted within-contrast fit manifest is incomplete.",
         call. = FALSE)
  }
  accepted
}

combine_within_full_fit_manifest <- function(
    reused_manifest,
    new_manifest
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  manifest <- dplyr::bind_rows(reused_manifest, new_manifest) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  counts <- table(factor(manifest$outcome, levels = outcomes))
  if (
    nrow(manifest) != 104L || anyDuplicated(paste(
      manifest$outcome, manifest$model_index
    )) || !identical(as.integer(counts), rep(26L, 4L)) ||
      any(manifest$family != "within") ||
      any(!file.exists(manifest$fit_path))
  ) {
    stop("The full within-contrast fit manifest is incomplete.",
         call. = FALSE)
  }
  manifest
}

make_within_full_fit_records <- function(fit_manifest) {
  if (
    !is.data.frame(fit_manifest) || nrow(fit_manifest) != 104L ||
      any(fit_manifest$family != "within") ||
      any(!file.exists(fit_manifest$fit_path))
  ) {
    stop("The full within-contrast fit manifest is invalid.",
         call. = FALSE)
  }
  purrr::map(seq_len(nrow(fit_manifest)), function(row) {
    purrr::map(
      as.list(fit_manifest[row, , drop = FALSE]),
      unname
    )
  })
}

diagnose_sensitivity_full_fit_file <- function(
    fit_record,
    family,
    ppc_draws = 200L
) {
  candidate <- read_sensitivity_fit_file(fit_record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), family) ||
      !identical(
        unname(as.character(candidate$outcome)),
        unname(as.character(fit_record$outcome))
      ) ||
      !identical(
        unname(as.integer(candidate$model_index)),
        unname(as.integer(fit_record$model_index))
      )
  ) {
    stop("A full sensitivity fit does not match its manifest record.",
         call. = FALSE)
  }
  get_primary_candidate_model_diagnostics(
    candidate,
    ppc_draws = ppc_draws
  )
}

diagnose_sensitivity_fit_path <- function(
    fit_path,
    family,
    ppc_draws = 200L
) {
  candidate <- read_sensitivity_fit_file(fit_path)
  if (!identical(unname(as.character(candidate$family)), family)) {
    stop("A sensitivity fit file has the wrong family.", call. = FALSE)
  }
  get_primary_candidate_model_diagnostics(
    candidate,
    ppc_draws = ppc_draws
  )
}

compile_sensitivity_full_diagnostic_audit <- function(
    diagnostics,
    family,
    expected_models = 104L
) {
  summary <- dplyr::bind_rows(purrr::map(diagnostics, "summary")) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    ) |>
    dplyr::mutate(
      ppc_warning = primary_ppc_warning(dplyr::pick(dplyr::everything())),
      status = dplyr::case_when(
        !diagnostic_pass ~ "REQUIRES REMEDIATION",
        ppc_warning ~ "PASS WITH PPC WARNING",
        TRUE ~ "PASS"
      )
    )
  if (
    nrow(summary) != expected_models ||
      anyDuplicated(paste(summary$outcome, summary$model_index))
  ) {
    stop(
      "The full ", family, " diagnostic audit is incomplete.",
      call. = FALSE
    )
  }
  parameter_failures <- purrr::map_dfr(diagnostics, function(result) {
    result$parameter_failures |>
      dplyr::mutate(
        outcome = result$summary$outcome[[1L]],
        model_index = result$summary$model_index[[1L]],
        .before = 1L
      )
  })
  list(
    family = family,
    diagnostic_summary = summary,
    parameter_failures = parameter_failures,
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

read_exposure_pilot_completed_checkpoint <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(
      "The completed exposure-pilot checkpoint is unavailable. Run ",
      "sensitivity/freeze_exposure_pilot_checkpoint.R first.",
      call. = FALSE
    )
  }
  checkpoint <- qs::qread(path)
  required <- c(
    "checkpoint_version", "family", "analysis_data_by_outcome",
    "priors_by_outcome", "model_specifications",
    "achieved_reused_fit_manifest", "assignment_reused_fit_manifest",
    "pilot_common_sample_audit", "pilot_diagnostic_summary",
    "pilot_parameter_failures", "audit"
  )
  manifests_valid <- function(manifest, family) {
    is.data.frame(manifest) && nrow(manifest) == 7L &&
      all(manifest$family == family) &&
      identical(
        unname(sort(manifest$model_index)),
        c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
      ) &&
      all(file.exists(manifest$fit_path))
  }
  if (
    !is.list(checkpoint) || !all(required %in% names(checkpoint)) ||
      !identical(checkpoint$checkpoint_version, 1L) ||
      !identical(checkpoint$family, "exposure") ||
      !isTRUE(checkpoint$audit$all_pilot_diagnostics_pass) ||
      !manifests_valid(checkpoint$achieved_reused_fit_manifest, "achieved") ||
      !manifests_valid(
        checkpoint$assignment_reused_fit_manifest, "assignment"
      )
  ) {
    stop("The completed exposure-pilot checkpoint is invalid.",
         call. = FALSE)
  }
  checkpoint
}

make_exposure_full_data_by_outcome <- function(
    analysis_data_by_outcome,
    protein_workbook
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.list(analysis_data_by_outcome) ||
      !setequal(names(analysis_data_by_outcome), outcomes)
  ) {
    stop("The full exposure source data are incomplete.", call. = FALSE)
  }
  prepared <- purrr::map(
    analysis_data_by_outcome[outcomes],
    prepare_assignment_exposure_data,
    protein_workbook = protein_workbook,
    common_sample = TRUE
  )
  achieved <- purrr::map(prepared, "achieved")
  assignment <- purrr::map(prepared, "assignment")
  audit <- purrr::imap_dfr(prepared, function(result, outcome) {
    achieved_data <- result$achieved
    assignment_data <- result$assignment
    identity_columns <- intersect(
      c(
        "study_no", "arm_no", "outcome_no", "outcome_group",
        "measurement_method", "unit_of_measurement", "time1", "time2"
      ),
      intersect(names(achieved_data), names(assignment_data))
    )
    identical_sample <- identical(
      achieved_data[, identity_columns, drop = FALSE],
      assignment_data[, identity_columns, drop = FALSE]
    )
    tibble::tibble(
      outcome = outcome,
      rows = nrow(achieved_data),
      studies = dplyr::n_distinct(achieved_data$study_no),
      arms = dplyr::n_distinct(achieved_data$arm_no),
      achieved_minimum = min(achieved_data$protein_intake_raw),
      achieved_maximum = max(achieved_data$protein_intake_raw),
      assignment_minimum = min(assignment_data$protein_intake_raw),
      assignment_maximum = max(assignment_data$protein_intake_raw),
      identical_common_sample = identical_sample
    )
  })
  valid <- purrr::map2_lgl(achieved, assignment, function(x, y) {
    is.data.frame(x) && is.data.frame(y) && nrow(x) == nrow(y) &&
      nrow(x) > 0L && all(is.finite(x$protein_intake_raw)) &&
      all(is.finite(y$protein_intake_raw))
  })
  if (
    !all(valid) || nrow(audit) != 4L ||
      !all(audit$identical_common_sample)
  ) {
    stop("The achieved/assignment common samples are invalid.",
         call. = FALSE)
  }
  list(
    achieved = achieved,
    assignment = assignment,
    audit = audit
  )
}

make_exposure_full_common_prediction_maxima <- function(
    achieved_data_by_outcome,
    assignment_data_by_outcome,
    protein_min = 0.8
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.list(achieved_data_by_outcome) ||
      !is.list(assignment_data_by_outcome) ||
      !setequal(names(achieved_data_by_outcome), outcomes) ||
      !setequal(names(assignment_data_by_outcome), outcomes)
  ) {
    stop("The paired exposure prediction data are incomplete.",
         call. = FALSE)
  }
  maxima <- stats::setNames(
    purrr::map_dbl(outcomes, function(outcome) {
      achieved <- achieved_data_by_outcome[[outcome]]
      assignment <- assignment_data_by_outcome[[outcome]]
      if (
        !is.data.frame(achieved) || !is.data.frame(assignment) ||
          !"protein_intake_raw" %in% names(achieved) ||
          !"protein_intake_raw" %in% names(assignment)
      ) {
        stop("An exposure prediction data set is invalid.", call. = FALSE)
      }
      min(
        max(achieved$protein_intake_raw, na.rm = TRUE),
        max(assignment$protein_intake_raw, na.rm = TRUE)
      )
    }),
    outcomes
  )
  if (any(!is.finite(maxima)) || any(maxima <= protein_min)) {
    stop("The common exposure prediction range is invalid.", call. = FALSE)
  }
  maxima
}

select_pending_exposure_full_specifications <- function(
    specifications,
    reused_fit_manifest,
    family
) {
  specification_keys <- paste(
    purrr::map_chr(specifications, "outcome"),
    purrr::map_int(specifications, "model_index")
  )
  reused_keys <- paste(
    reused_fit_manifest$outcome,
    reused_fit_manifest$model_index
  )
  expected_reused <- paste(
    "whole_body", c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
  )
  if (
    !family %in% c("achieved", "assignment") ||
      length(specifications) != 104L || anyDuplicated(specification_keys) ||
      nrow(reused_fit_manifest) != 7L ||
      any(reused_fit_manifest$family != family) ||
      !setequal(reused_keys, expected_reused) ||
      !all(reused_keys %in% specification_keys) ||
      any(!file.exists(reused_fit_manifest$fit_path))
  ) {
    stop("The full exposure specification/reuse set is invalid.",
         call. = FALSE)
  }
  specifications[!specification_keys %in% reused_keys]
}

prepare_exposure_full_input <- function(
    specification,
    exposure_data_by_outcome,
    priors_by_outcome,
    family,
    config
) {
  outcome <- specification$outcome
  data <- exposure_data_by_outcome[[family]][[outcome]]
  priors <- priors_by_outcome[[outcome]]
  if (
    !family %in% c("achieved", "assignment") ||
      is.null(data) || is.null(priors) || !outcome %in% config$outcomes
  ) {
    stop("A full exposure outcome input is unavailable.", call. = FALSE)
  }
  specification$sensitivity_exposure_expansion <- TRUE
  specification$sensitivity_exposure_definition <- family
  specification$sensitivity_exposure_sampling <- paste0(
    config$iter_warmup, "_warmup_", config$iter_sampling, "_sampling"
  )
  output <- prepare_sensitivity_final_candidate(
    specification, priors, data
  )
  output$family <- family
  output
}

audit_exposure_full_expansion <- function(
    achieved_pending_specifications,
    assignment_pending_specifications,
    achieved_inputs,
    assignment_inputs,
    exposure_data_by_outcome,
    achieved_reused_fit_manifest,
    assignment_reused_fit_manifest,
    config
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  audit_one <- function(specifications, inputs, manifest, family) {
    specification_keys <- paste(
      purrr::map_chr(specifications, "outcome"),
      purrr::map_int(specifications, "model_index")
    )
    input_keys <- paste(
      purrr::map_chr(inputs, ~ .x$specification$outcome),
      purrr::map_int(inputs, ~ .x$specification$model_index)
    )
    counts <- table(factor(
      purrr::map_chr(inputs, ~ .x$specification$outcome),
      levels = outcomes
    ))
    valid_inputs <- purrr::map_lgl(inputs, function(input) {
      identical(input$family, family) &&
        is.character(input$formula) && length(input$formula) == 1L &&
        is.data.frame(input$data) && nrow(input$data) > 0L
    })
    length(specifications) == 97L && length(inputs) == 97L &&
      nrow(manifest) == 7L && all(manifest$family == family) &&
      !anyDuplicated(specification_keys) && !anyDuplicated(input_keys) &&
      setequal(specification_keys, input_keys) &&
      identical(as.integer(counts), c(19L, 26L, 26L, 26L)) &&
      all(valid_inputs)
  }
  achieved_valid <- audit_one(
    achieved_pending_specifications, achieved_inputs,
    achieved_reused_fit_manifest, "achieved"
  )
  assignment_valid <- audit_one(
    assignment_pending_specifications, assignment_inputs,
    assignment_reused_fit_manifest, "assignment"
  )
  if (
    !achieved_valid || !assignment_valid ||
      nrow(exposure_data_by_outcome$audit) != 4L ||
      !all(exposure_data_by_outcome$audit$identical_common_sample) ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("The full exposure expansion failed its structural audit.",
         call. = FALSE)
  }
  list(
    checks = tibble::tibble(
      check = c(
        "complete paired candidate registries",
        "reused passing whole-body pilots",
        "new independently checkpointed fits",
        "pilot-informed sampling schedule",
        "identical achieved/assignment samples"
      ),
      status = "PASS",
      detail = c(
        "208 candidates: 26 models x 4 outcomes x 2 exposure definitions",
        "14 fits: 7 achieved and 7 assignment-derived",
        "194 fits: 97 achieved and 97 assignment-derived",
        "1,500 warmup and 2,000 retained iterations per chain",
        "same rows, studies, arms, outcomes, and follow-up contrasts"
      )
    ),
    common_sample_summary = exposure_data_by_outcome$audit,
    pending_fits = 194L,
    reused_fits = 14L,
    structure_pass = TRUE
  )
}

fit_exposure_full_candidate_to_file <- function(
    input,
    enabled,
    config,
    family,
    stage = "initial"
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip(family))
  if (!family %in% c("achieved", "assignment")) {
    stop("Unknown full exposure family.", call. = FALSE)
  }
  outcome_offset <- match(input$specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the full exposure fit.", call. = FALSE)
  }
  family_offset <- if (family == "achieved") 0L else 1000L
  local_config <- config
  local_config$seed <- as.integer(
    config$seed + family_offset + 100L * (outcome_offset - 1L) +
      input$specification$model_index
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = family
  )
  write_sensitivity_fit_file(candidate, stage = stage)
}

combine_exposure_full_fit_manifest <- function(
    reused_manifest,
    new_manifest,
    family
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  manifest <- dplyr::bind_rows(reused_manifest, new_manifest) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  counts <- table(factor(manifest$outcome, levels = outcomes))
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      any(manifest$family != family) ||
      any(!file.exists(manifest$fit_path))
  ) {
    stop("The full exposure fit manifest is incomplete.", call. = FALSE)
  }
  manifest
}

make_exposure_full_fit_records <- function(fit_manifest, family) {
  if (
    !is.data.frame(fit_manifest) || nrow(fit_manifest) != 104L ||
      any(fit_manifest$family != family) ||
      any(!file.exists(fit_manifest$fit_path))
  ) {
    stop("The full exposure fit manifest is invalid.", call. = FALSE)
  }
  purrr::map(seq_len(nrow(fit_manifest)), function(row) {
    purrr::map(as.list(fit_manifest[row, , drop = FALSE]), unname)
  })
}

combine_exposure_full_diagnostic_audits <- function(
    achieved_audit,
    assignment_audit
) {
  if (
    achieved_audit$models != 104L || assignment_audit$models != 104L ||
      !identical(achieved_audit$family, "achieved") ||
      !identical(assignment_audit$family, "assignment")
  ) {
    stop("The paired full exposure diagnostic audits are incomplete.",
         call. = FALSE)
  }
  summary <- dplyr::bind_rows(
    dplyr::mutate(
      achieved_audit$diagnostic_summary,
      exposure = "achieved", .before = 1L
    ),
    dplyr::mutate(
      assignment_audit$diagnostic_summary,
      exposure = "assignment_derived", .before = 1L
    )
  )
  list(
    family = "exposure",
    diagnostic_summary = summary,
    parameter_failures = dplyr::bind_rows(
      dplyr::mutate(
        achieved_audit$parameter_failures,
        exposure = "achieved", .before = 1L
      ),
      dplyr::mutate(
        assignment_audit$parameter_failures,
        exposure = "assignment_derived", .before = 1L
      )
    ),
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    achieved = achieved_audit,
    assignment = assignment_audit,
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

select_exposure_full_remedy_specifications <- function(
    specifications,
    initial_audit,
    family
) {
  failed <- initial_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::transmute(
      outcome = as.character(outcome),
      model_index = as.integer(model_index),
      key = paste(outcome, model_index)
    )
  specification_keys <- paste(
    purrr::map_chr(specifications, "outcome"),
    purrr::map_int(specifications, "model_index")
  )
  selected <- specifications[specification_keys %in% failed$key]
  selected_keys <- paste(
    purrr::map_chr(selected, "outcome"),
    purrr::map_int(selected, "model_index")
  )
  if (
    !identical(initial_audit$family, family) || nrow(failed) < 1L ||
      anyDuplicated(failed$key) || anyDuplicated(specification_keys) ||
      length(selected) != nrow(failed) ||
      !setequal(selected_keys, failed$key)
  ) {
    stop("The full exposure remedy specification set is invalid.",
         call. = FALSE)
  }
  selected
}

prepare_exposure_full_remedy_input <- function(
    specification,
    exposure_data_by_outcome,
    priors_by_outcome,
    family,
    config
) {
  specification$adapt_delta <- max(
    specification$adapt_delta,
    config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth,
    config$max_treedepth_floor
  )
  specification$sensitivity_exposure_full_remedy <-
    "longer_sampling_high_adapt_delta"
  prepare_exposure_full_input(
    specification = specification,
    exposure_data_by_outcome = exposure_data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    family = family,
    config = config
  )
}

combine_exposure_full_remedy_diagnostic_audit <- function(
    initial_audit,
    remedy_diagnostics,
    family
) {
  remedy_summary <- dplyr::bind_rows(purrr::map(
    remedy_diagnostics, "summary"
  ))
  remedy_keys <- paste(remedy_summary$outcome, remedy_summary$model_index)
  expected_keys <- initial_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::transmute(key = paste(outcome, model_index)) |>
    dplyr::pull(key)
  if (
    !identical(initial_audit$family, family) ||
      nrow(remedy_summary) != length(expected_keys) ||
      anyDuplicated(remedy_keys) || !setequal(remedy_keys, expected_keys)
  ) {
    stop("The full exposure remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  summary <- initial_audit$diagnostic_summary |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_summary) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    ) |>
    dplyr::mutate(
      ppc_warning = primary_ppc_warning(dplyr::pick(dplyr::everything())),
      status = dplyr::case_when(
        !diagnostic_pass ~ "REQUIRES REMEDIATION",
        ppc_warning ~ "PASS WITH PPC WARNING",
        TRUE ~ "PASS"
      )
    )
  remedy_parameter_failures <- purrr::map_dfr(
    remedy_diagnostics,
    function(result) {
      result$parameter_failures |>
        dplyr::mutate(
          outcome = result$summary$outcome[[1L]],
          model_index = result$summary$model_index[[1L]],
          .before = 1L
        )
    }
  )
  parameter_failures <- initial_audit$parameter_failures |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_parameter_failures)
  if (
    nrow(summary) != 104L ||
      anyDuplicated(paste(summary$outcome, summary$model_index))
  ) {
    stop("The remedied exposure diagnostic audit is incomplete.",
         call. = FALSE)
  }
  list(
    family = family,
    diagnostic_summary = summary,
    parameter_failures = parameter_failures,
    remedy_diagnostic_summary = remedy_summary,
    remedy_keys = remedy_keys,
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

combine_exposure_full_accepted_fit_manifest <- function(
    initial_manifest,
    remedy_manifest,
    final_diagnostic_audit,
    family
) {
  if (
    !identical(final_diagnostic_audit$family, family) ||
      !isTRUE(final_diagnostic_audit$all_sampling_diagnostics_pass)
  ) {
    stop("Exposure fits cannot be accepted before diagnostics pass.",
         call. = FALSE)
  }
  remedy_keys <- paste(remedy_manifest$outcome, remedy_manifest$model_index)
  accepted <- initial_manifest |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_manifest) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  if (
    nrow(accepted) != 104L ||
      anyDuplicated(paste(accepted$outcome, accepted$model_index)) ||
      any(accepted$family != family) ||
      any(!file.exists(accepted$fit_path))
  ) {
    stop("The accepted exposure fit manifest is incomplete.",
         call. = FALSE)
  }
  accepted
}

make_strict_within_study_contrasts <- function(data) {
  keys <- c(
    "study_no", "outcome_group", "measurement_method",
    "unit_of_measurement", "time1", "time2",
    "resistance_exercise_code"
  )
  protein <- data |>
    dplyr::filter(condition == "protein")
  comparator <- data |>
    dplyr::filter(condition %in% c("control", "placebo"))
  paired <- dplyr::inner_join(
    protein,
    comparator,
    by = keys,
    suffix = c("_protein", "_comparator"),
    relationship = "many-to-many"
  ) |>
    dplyr::filter(
      arm_no_protein != arm_no_comparator,
      is.finite(yi_protein), is.finite(vi_protein), vi_protein > 0,
      is.finite(yi_comparator), is.finite(vi_comparator), vi_comparator > 0,
      is.finite(protein_intake_protein),
      is.finite(protein_intake_comparator)
    ) |>
    dplyr::transmute(
      study_no,
      outcome_group,
      measurement_method,
      unit_of_measurement,
      time1,
      time2,
      resistance_exercise_code,
      resistance_exercise_factor = resistance_exercise_factor_protein,
      protein_arm_no = arm_no_protein,
      comparator_arm_no = arm_no_comparator,
      protein_arm_no_extracted = arm_no_extracted_protein,
      comparator_arm_no_extracted = arm_no_extracted_comparator,
      yi = yi_protein - yi_comparator,
      vi = vi_protein + vi_comparator,
      protein_intake_treatment = protein_intake_protein,
      protein_intake_comparator = protein_intake_comparator,
      protein_intake_raw_treatment = protein_intake_raw_protein,
      protein_intake_raw_comparator = protein_intake_raw_comparator,
      duration_centre = (duration_centre_protein +
        duration_centre_comparator) / 2,
      age_centre = (age_centre_protein + age_centre_comparator) / 2
    ) |>
    dplyr::mutate(
      arm_no = dplyr::dense_rank(interaction(
        study_no, protein_arm_no, comparator_arm_no, drop = TRUE
      )),
      outcome_no = dplyr::dense_rank(interaction(
        study_no, protein_arm_no, comparator_arm_no, outcome_group,
        measurement_method, unit_of_measurement, drop = TRUE
      )),
      study_factor = factor(study_no),
      resistance_exercise_factor = factor(
        resistance_exercise_factor,
        levels = c("yes", "no")
      )
    )

  if (nrow(paired) < 1L || anyNA(paired$yi) || any(paired$vi <= 0)) {
    stop("No valid strictly within-study contrasts could be formed.",
         call. = FALSE)
  }
  paired
}

prepare_within_contrast_pilot_input <- function(data, specification, priors) {
  priors <- recenter_within_contrast_priors(priors)
  features <- sensitivity_model_features(specification$model_index)
  if (!features$includes_protein) {
    output <- prepare_sensitivity_final_candidate(
      specification, priors, data
    )
    output$family <- "within"
    output$basis_spec <- NULL
    output$contrast_names <- NULL
    return(output)
  }
  pooled_x <- c(data$protein_intake_treatment, data$protein_intake_comparator)
  basis_data <- tibble::tibble(
    study_no = rep(data$study_no, 2L),
    arm_no = seq_along(pooled_x),
    protein_intake = pooled_x
  )
  basis_spec <- make_sensitivity_basis_spec(basis_data, specification)
  treatment_basis <- predict_sensitivity_basis(
    basis_spec, data$protein_intake_treatment
  )
  comparator_basis <- predict_sensitivity_basis(
    basis_spec, data$protein_intake_comparator
  )
  contrast_basis <- treatment_basis - comparator_basis
  contrast_names <- paste0(
    "protein_contrast_", seq_len(ncol(contrast_basis))
  )
  contrast_basis <- tibble::as_tibble(
    contrast_basis,
    .name_repair = "minimal"
  )
  names(contrast_basis) <- contrast_names
  data <- dplyr::bind_cols(data, contrast_basis)
  if (isTRUE(specification$standardised_moderators)) {
    data <- add_primary_standardised_moderators(data)
  }
  formula <- make_sensitivity_basis_formula(
    specification,
    contrast_names = contrast_names
  )
  priors <- prepare_sensitivity_basis_priors(
    priors,
    specification,
    coefficient_names = contrast_names,
    basis_spec = basis_spec
  )

  list(
    family = "within",
    specification = specification,
    formula = formula,
    priors = priors,
    data = data,
    basis_spec = basis_spec,
    contrast_names = contrast_names
  )
}

prepare_assignment_exposure_data <- function(
    data,
    protein_workbook,
    common_sample = TRUE
) {
  arms <- protein_workbook$arm_characteristics |>
    dplyr::transmute(
      study_no,
      arm_no_extracted = arm_no,
      baseline_g_kg = dplyr::coalesce(
        baseline_dietary_protein_intake_g_kg_day,
        baseline_dietary_protein_intake_g_day / mean_body_mass_kg
      ),
      prescribed_protein_g_day =
        absolute_protein_dose_g_day_note_if_placebo_condition_then_check_whether_it_includes_any_protein,
      assignment_body_mass_kg = mean_body_mass_kg
    )
  output <- data |>
    dplyr::left_join(
      arms,
      by = c("study_no", "arm_no_extracted"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      prescribed_protein_g_day = dplyr::case_when(
        condition %in% c("control", "placebo") &
          is.na(prescribed_protein_g_day) ~ 0,
        TRUE ~ prescribed_protein_g_day
      ),
      assigned_increment_g_kg =
        prescribed_protein_g_day / assignment_body_mass_kg,
      assigned_protein_intake_raw = baseline_g_kg + assigned_increment_g_kg,
      achieved_protein_intake_raw = protein_intake_raw,
      assignment_available = is.finite(assigned_protein_intake_raw),
      achieved_available = is.finite(achieved_protein_intake_raw)
    )
  if (isTRUE(common_sample)) {
    output <- dplyr::filter(
      output,
      assignment_available,
      achieved_available
    )
  }

  recode_exposure <- function(input, column, label) {
    comparator_median <- input |>
      dplyr::filter(condition %in% c("control", "placebo")) |>
      dplyr::distinct(study_no, arm_no, .data[[column]]) |>
      dplyr::summarise(value = stats::median(.data[[column]], na.rm = TRUE)) |>
      dplyr::pull(value)
    if (length(comparator_median) != 1L || !is.finite(comparator_median)) {
      stop("Could not centre the ", label, " exposure.", call. = FALSE)
    }
    input |>
      dplyr::mutate(
        protein_intake_raw = .data[[column]],
        protein_intake = protein_intake_raw - comparator_median,
        protein_median = comparator_median,
        sensitivity_exposure = label,
        study_factor = factor(study_no)
      )
  }

  list(
    achieved = recode_exposure(
      output, "achieved_protein_intake_raw", "achieved"
    ),
    assignment = recode_exposure(
      output, "assigned_protein_intake_raw", "assignment_derived"
    ),
    audit = tibble::tibble(
      common_sample = isTRUE(common_sample),
      rows = nrow(output),
      studies = dplyr::n_distinct(output$study_no),
      arms = dplyr::n_distinct(output$arm_no),
      protein_arms = dplyr::n_distinct(
        output$arm_no[output$condition == "protein"]
      ),
      comparator_arms = dplyr::n_distinct(
        output$arm_no[output$condition %in% c("control", "placebo")]
      )
    )
  )
}

prepare_sensitivity_final_candidate <- function(
    specification,
    priors,
    data
) {
  priors <- normalize_primary_brms_priors(priors)
  fitting_specification <- specification

  if (isTRUE(specification$targeted_function_tail_priors)) {
    priors <- add_primary_function_tail_priors(
      priors,
      include_interactions =
        specification$targeted_piecewise_interaction_priors
    )
    priors <- c(
      priors,
      brms::set_prior("student_t(3, 0, 1)", class = "sd")
    )
    fitting_specification$regularised_variance_priors <- FALSE
  } else if (isTRUE(specification$targeted_piecewise_slope_priors)) {
    priors <- add_primary_piecewise_slope_priors(priors)
    if (isTRUE(specification$regularised_variance_priors)) {
      priors <- c(
        priors,
        brms::set_prior("student_t(3, 0, 1)", class = "sd")
      )
      fitting_specification$regularised_variance_priors <- FALSE
    }
  }
  if (isTRUE(fitting_specification$regularised_variance_priors)) {
    priors <- add_regularised_variance_priors(priors)
  }
  if (isTRUE(specification$standardised_moderators)) {
    data <- add_primary_standardised_moderators(data)
    priors <- add_primary_standardised_moderator_priors(priors)
  }
  if (isTRUE(specification$half_unit_protein)) {
    data <- add_primary_half_unit_protein(data)
  }
  if (isTRUE(specification$natural_spline)) {
    priors <- remove_global_primary_prior(priors, "sds")
  }
  if (isTRUE(specification$targeted_smooth_prior)) {
    priors <- add_primary_targeted_smooth_prior(
      priors,
      scale = specification$smooth_prior_scale
    )
  }

  list(
    family = "exposure",
    specification = specification,
    formula = specification$formula,
    priors = priors,
    data = data
  )
}

fit_sensitivity_brms_pilot <- function(
    input,
    enabled,
    config,
    family = input$family
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip(family))
  model <- brms::brm(
    formula = stats::as.formula(input$formula),
    data = input$data,
    prior = input$priors,
    chains = config$chains,
    cores = config$parallel_chains,
    seed = config$seed,
    warmup = config$iter_warmup,
    iter = config$iter_warmup + config$iter_sampling,
    control = list(
      adapt_delta = max(
        input$specification$adapt_delta,
        config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        config$max_treedepth_floor
      )
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  specification <- input$specification
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  family_suffix <- paste0(" [", family, "]")
  model_name <- if (
    stringr::str_ends(
      specification$model_name,
      stringr::fixed(family_suffix)
    )
  ) {
    specification$model_name
  } else {
    paste0(specification$model_name, family_suffix)
  }
  list(
    outcome = specification$outcome,
    model_index = specification$model_index,
    model_name = model_name,
    specification = specification,
    family = family,
    basis_spec = if (is.null(input$basis_spec)) NULL else input$basis_spec,
    within_names = if (is.null(input$within_names)) NULL else input$within_names,
    between_names = if (is.null(input$between_names)) NULL else input$between_names,
    contrast_names = if (is.null(input$contrast_names)) NULL else input$contrast_names,
    model = model
  )
}

collect_completed_sensitivity_fits <- function(fits, family) {
  if (all(purrr::map_lgl(fits, is_sensitivity_pilot_skip))) {
    return(new_sensitivity_pilot_skip(family))
  }
  if (any(purrr::map_lgl(fits, is_sensitivity_pilot_skip))) {
    stop("A sensitivity pilot is only partially fitted.", call. = FALSE)
  }
  fits[order(purrr::map_int(fits, "model_index"))]
}

select_sensitivity_pilot_input <- function(inputs, model_index) {
  selected <- purrr::keep(
    inputs,
    ~ identical(.x$specification$model_index, as.integer(model_index))
  )
  if (length(selected) != 1L) {
    stop("Expected exactly one sensitivity pilot input.", call. = FALSE)
  }
  selected[[1L]]
}

combine_cmdstan_smoke_and_remaining_fits <- function(
    smoke_fit,
    remaining_fits,
    expected_model_indices,
    family
) {
  if (
    is_sensitivity_pilot_skip(smoke_fit) ||
      is_sensitivity_pilot_skip(remaining_fits)
  ) {
    return(new_sensitivity_pilot_skip(family))
  }
  combined <- c(list(smoke_fit), remaining_fits)
  combined <- combined[order(purrr::map_int(combined, "model_index"))]
  if (!identical(
    unname(purrr::map_int(combined, "model_index")),
    sort(unique(as.integer(expected_model_indices)))
  )) {
    stop("The CmdStan smoke/full pilot ensemble is incomplete.",
         call. = FALSE)
  }
  combined
}

validate_cmdstan_pilot_fits <- function(fits, diagnostics, family) {
  if (
    is_sensitivity_pilot_skip(fits) ||
      is_sensitivity_pilot_skip(diagnostics)
  ) {
    return(new_sensitivity_pilot_skip(family))
  }
  if (
    nrow(diagnostics) != length(fits) ||
      anyNA(diagnostics$diagnostic_pass) ||
      !all(diagnostics$diagnostic_pass)
  ) {
    stop(
      paste0("The ", family, " pilot ensemble has unresolved diagnostics."),
      call. = FALSE
    )
  }
  fits
}

combine_mundlak_remedy_fits <- function(
    original_fits,
    remedy_fits,
    expected_model_indices
) {
  if (
    is_sensitivity_pilot_skip(original_fits) ||
      is_sensitivity_pilot_skip(remedy_fits)
  ) {
    return(new_sensitivity_pilot_skip("mundlak"))
  }
  remedy_indices <- purrr::map_int(remedy_fits, "model_index")
  retained <- purrr::keep(
    original_fits,
    ~ !.x$model_index %in% remedy_indices
  )
  combined <- c(retained, remedy_fits)
  combined <- combined[order(purrr::map_int(combined, "model_index"))]
  if (!identical(
    unname(purrr::map_int(combined, "model_index")),
    sort(as.integer(expected_model_indices))
  )) {
    stop("The Mundlak remedy ensemble is incomplete.", call. = FALSE)
  }
  combined
}

combine_sensitivity_remedy_fits <- function(
    original_fits,
    remedy_fits,
    expected_model_indices,
    family
) {
  if (
    is_sensitivity_pilot_skip(original_fits) ||
      is_sensitivity_pilot_skip(remedy_fits)
  ) {
    return(new_sensitivity_pilot_skip(family))
  }
  remedy_indices <- unname(purrr::map_int(remedy_fits, "model_index"))
  retained <- purrr::keep(
    original_fits,
    ~ !.x$model_index %in% remedy_indices
  )
  combined <- c(retained, remedy_fits)
  combined <- combined[order(purrr::map_int(combined, "model_index"))]
  if (!identical(
    unname(purrr::map_int(combined, "model_index")),
    sort(as.integer(expected_model_indices))
  )) {
    stop(
      paste0("The ", family, " remedy ensemble is incomplete."),
      call. = FALSE
    )
  }
  combined
}

relabel_sensitivity_fit <- function(candidate, family) {
  if (is_sensitivity_pilot_skip(candidate)) {
    return(new_sensitivity_pilot_skip(family))
  }
  if (
    length(candidate) < 1L ||
      !is.numeric(candidate$model_index) ||
      length(candidate$model_index) != 1L ||
      !inherits(candidate$model, "brmsfit")
  ) {
    stop("Cannot relabel an invalid sensitivity fit.", call. = FALSE)
  }
  candidate$family <- family
  candidate$model_name <- paste0(
    sub(" \\[[^]]+\\]$", "", candidate$model_name),
    " [", family, "]"
  )
  candidate
}

diagnose_sensitivity_brms_fits <- function(fits, ppc_draws = 200L) {
  if (is_sensitivity_pilot_skip(fits)) return(fits)
  purrr::map(fits, get_primary_candidate_model_diagnostics,
             ppc_draws = ppc_draws)
}

compute_sensitivity_brms_weights <- function(fits) {
  if (is_sensitivity_pilot_skip(fits)) return(fits)
  model_names <- purrr::map_chr(fits, "model_name")
  fit_draws <- purrr::map_int(
    fits,
    ~ as.integer(brms::ndraws(.x$model))
  )
  loo_draws <- min(fit_draws)
  draw_ids <- purrr::map(fit_draws, function(n_draws) {
    if (n_draws == loo_draws) return(seq_len(n_draws))
    as.integer(floor(seq(
      from = 0,
      to = n_draws - 1,
      length.out = loo_draws
    )) + 1L)
  })
  if (
    any(purrr::map_int(draw_ids, length) != loo_draws) ||
      any(purrr::map_int(draw_ids, anyDuplicated) > 0L)
  ) {
    stop("Could not harmonise posterior draws for pilot LOO.", call. = FALSE)
  }
  loo_objects <- purrr::map2(fits, draw_ids, ~ brms::loo(
    .x$model,
    draw_ids = .y,
    pointwise = TRUE,
    moment_match = FALSE,
    reloo = FALSE
  ))
  names(loo_objects) <- model_names
  weights <- loo::loo_model_weights(loo_objects, method = "stacking")
  weights <- pmax(as.numeric(weights), 0)
  weights <- weights / sum(weights)
  names(weights) <- model_names
  names(fit_draws) <- model_names
  list(
    loo = loo_objects,
    weights = weights,
    fit_draws = fit_draws,
    loo_draws = loo_draws
  )
}

make_sensitivity_loo_draw_ids <- function(n_draws, target_draws) {
  n_draws <- as.integer(n_draws)
  target_draws <- as.integer(target_draws)
  if (
    length(n_draws) != 1L || length(target_draws) != 1L ||
      is.na(n_draws) || is.na(target_draws) ||
      target_draws < 1L || target_draws > n_draws
  ) {
    stop("Invalid posterior draw counts for sensitivity LOO.", call. = FALSE)
  }
  if (n_draws == target_draws) return(seq_len(n_draws))
  ids <- as.integer(floor(seq(
    from = 0,
    to = n_draws - 1,
    length.out = target_draws
  )) + 1L)
  if (length(ids) != target_draws || anyDuplicated(ids)) {
    stop("Could not select sensitivity LOO draws.", call. = FALSE)
  }
  ids
}

make_sensitivity_mundlak_integrated_loo_jobs <- function(
    fits,
    conditional_loo,
    family = "mundlak"
) {
  if (
    is_sensitivity_pilot_skip(fits) ||
      is_sensitivity_pilot_skip(conditional_loo)
  ) {
    return(new_sensitivity_pilot_skip(family))
  }
  loo_objects <- conditional_loo$loo
  model_names <- purrr::map_chr(fits, "model_name")
  if (
    length(loo_objects) != length(fits) ||
      !setequal(names(loo_objects), model_names) ||
      !all(purrr::map_lgl(loo_objects, ~ inherits(.x, "loo"))) ||
      !is.numeric(conditional_loo$loo_draws) ||
      length(conditional_loo$loo_draws) != 1L
  ) {
    stop(
      paste0("The conditional ", family, " LOO results are incomplete."),
      call. = FALSE
    )
  }
  purrr::map(fits, function(candidate) {
    list(
      fitted_candidate = candidate,
      conditional_loo = loo_objects[[candidate$model_name]],
      loo_draws = as.integer(conditional_loo$loo_draws)
    )
  })
}

get_sensitivity_mundlak_integrated_components <- function(
    fitted_candidate,
    draw_ids,
    expected_family = "mundlak",
    slope_description = "mundlak_within_basis_slopes"
) {
  model <- fitted_candidate$model
  manifest <- get_primary_loo_observation_manifest(fitted_candidate)
  if (
    !inherits(model, "brmsfit") ||
      !identical(fitted_candidate$family, expected_family)
  ) {
    stop(
      paste0("Integrated ", expected_family, " LOO requires a matching brms fit."),
      call. = FALSE
    )
  }
  total_draws <- as.integer(brms::ndraws(model))
  draw_ids <- as.integer(draw_ids)
  if (
    length(draw_ids) < 1L || anyNA(draw_ids) || anyDuplicated(draw_ids) ||
      min(draw_ids) < 1L || max(draw_ids) > total_draws
  ) {
    stop("Invalid Mundlak LOO draw identifiers.", call. = FALSE)
  }

  standata <- brms::standata(model)
  if (!"J_3" %in% names(standata)) {
    stop("The Mundlak study grouping index is unavailable.", call. = FALSE)
  }
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
    stop("Stan study groups do not match the Mundlak data.", call. = FALSE)
  }

  variance_components <- brms::VarCorr(model, summary = FALSE)
  study_sd <- as.matrix(variance_components$study_no$sd)
  coefficient_names <- colnames(study_sd)
  n_coefficients <- ncol(study_sd)
  design_names <- paste0("Z_3_", seq_len(n_coefficients))
  if (
    nrow(study_sd) != total_draws ||
      length(coefficient_names) != n_coefficients ||
      !identical(coefficient_names[[1L]], "Intercept") ||
      !all(design_names %in% names(standata))
  ) {
    stop("The Mundlak study covariance design is incomplete.", call. = FALSE)
  }
  study_design <- do.call(cbind, lapply(
    design_names,
    function(x) as.numeric(standata[[x]])
  ))
  if (
    !identical(dim(study_design), c(nrow(manifest), n_coefficients)) ||
      anyNA(study_design) || any(!is.finite(study_design)) ||
      any(abs(study_design[, 1L] - 1) > sqrt(.Machine$double.eps))
  ) {
    stop("The Mundlak study design matrix is invalid.", call. = FALSE)
  }

  common_variables <- c(
    residual = "sigma",
    arm = "sd_arm_no__Intercept",
    outcome = "sd_outcome_no__Intercept"
  )
  common_draws <- as.data.frame(extract_primary_exact_draws(
    model,
    common_variables
  ))[draw_ids, , drop = FALSE]
  uncorrelated_study_effects <- isTRUE(
    fitted_candidate$specification$uncorrelated_study_effects
  )
  if (n_coefficients > 1L && !uncorrelated_study_effects) {
    study_covariance <- variance_components$study_no$cov
    if (
      length(dim(study_covariance)) != 3L ||
        dim(study_covariance)[[1L]] != total_draws ||
        dim(study_covariance)[[2L]] != n_coefficients ||
        dim(study_covariance)[[3L]] != n_coefficients
    ) {
      stop("The Mundlak study covariance draws are incomplete.",
           call. = FALSE)
    }
    covariance_draws <- study_covariance[draw_ids, , , drop = FALSE]
  }
  n_draws <- length(draw_ids)
  n_observations <- nrow(manifest)
  if (n_coefficients == 1L || uncorrelated_study_effects) {
    # `brms::VarCorr()` correctly omits `cov` for `(term || study)` effects.
    # Each independent standard-normal component is therefore its design
    # column multiplied directly by its draw-specific standard deviation.
    study_basis <- lapply(seq_len(n_coefficients), function(component) {
      outer(
        study_sd[draw_ids, component],
        study_design[, component]
      )
    })
  } else {
    study_basis <- lapply(
      seq_len(n_coefficients),
      function(x) matrix(NA_real_, n_draws, n_observations)
    )
    for (draw in seq_len(n_draws)) {
      covariance <- matrix(
        covariance_draws[draw, , ],
        nrow = n_coefficients,
        ncol = n_coefficients
      )
      cholesky <- tryCatch(t(chol(covariance)), error = function(e) NULL)
      if (is.null(cholesky)) {
        stop("A Mundlak study covariance draw is not positive definite.",
             call. = FALSE)
      }
      transformed_design <- study_design %*% cholesky
      for (component in seq_len(n_coefficients)) {
        study_basis[[component]][draw, ] <-
          transformed_design[, component]
      }
    }
  }
  linear_predictor <- brms::posterior_linpred(
    model,
    transform = FALSE,
    re_formula = NA,
    draw_ids = draw_ids
  )
  if (!identical(dim(linear_predictor), c(n_draws, n_observations))) {
    stop("The Mundlak integrated linear predictor is misaligned.",
         call. = FALSE)
  }
  list(
    structure = if (n_coefficients == 1L) {
      "study_intercept"
    } else {
      paste0(
        if (isTRUE(fitted_candidate$specification$uncorrelated_study_effects)) {
          "uncorrelated_"
        } else {
          "correlated_"
        },
        "study_intercept_and_",
        n_coefficients - 1L,
        "_", slope_description
      )
    },
    linear_predictor = linear_predictor,
    residual_sd = common_draws[[common_variables[["residual"]]]],
    arm_intercept_sd = common_draws[[common_variables[["arm"]]]],
    outcome_intercept_sd = common_draws[[common_variables[["outcome"]]]],
    study_basis = study_basis,
    chain_id = as.integer(common_draws$.chain),
    coefficient_names = coefficient_names
  )
}

compute_sensitivity_mundlak_integrated_psis_loo <- function(
    job,
    k_threshold = 0.7,
    cores = 1L,
    family = "mundlak",
    slope_description = "mundlak_within_basis_slopes"
) {
  if (is_sensitivity_pilot_skip(job)) return(job)
  candidate <- job$fitted_candidate
  manifest <- get_primary_loo_observation_manifest(candidate)
  conditional_loo <- job$conditional_loo
  if (
    !inherits(conditional_loo, "loo") ||
      nrow(conditional_loo$pointwise) != nrow(manifest)
  ) {
    stop("The Mundlak conditional LOO object is invalid.", call. = FALSE)
  }
  draw_ids <- make_sensitivity_loo_draw_ids(
    brms::ndraws(candidate$model),
    job$loo_draws
  )
  components <- get_sensitivity_mundlak_integrated_components(
    candidate,
    draw_ids,
    expected_family = family,
    slope_description = slope_description
  )
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
  if (
    !inherits(integrated_loo, "loo") ||
      nrow(integrated_loo$pointwise) != nrow(manifest)
  ) {
    stop("The integrated Mundlak LOO result is invalid.", call. = FALSE)
  }
  list(
    outcome = candidate$outcome,
    model_index = candidate$model_index,
    model_name = candidate$model_name,
    study_structure = components$structure,
    observation_manifest = manifest,
    loo = integrated_loo,
    conditional_loo = conditional_loo,
    remedy = "integrated_observation_level_psis",
    k_threshold = k_threshold,
    comparison = dplyr::bind_rows(
      summarise_primary_loo_method(
        conditional_loo,
        "conditional_observation_level_psis",
        k_threshold
      ),
      summarise_primary_loo_method(
        integrated_loo,
        "integrated_observation_level_psis",
        k_threshold
      )
    ) |>
      dplyr::mutate(
        outcome = candidate$outcome,
        model_index = candidate$model_index,
        model_name = candidate$model_name,
        study_structure = components$structure,
        .before = 1L
      )
  )
}

get_sensitivity_exposure_integrated_components <- function(
    fitted_candidate,
    draw_ids,
    expected_family
) {
  if (!identical(fitted_candidate$family, expected_family)) {
    stop("The exposure family does not match its integrated-LOO job.",
         call. = FALSE)
  }
  components <- get_primary_integrated_loo_components(fitted_candidate)
  total_draws <- nrow(components$linear_predictor)
  draw_ids <- as.integer(draw_ids)
  if (
    length(draw_ids) < 1L || anyNA(draw_ids) || anyDuplicated(draw_ids) ||
      min(draw_ids) < 1L || max(draw_ids) > total_draws
  ) {
    stop("Invalid exposure LOO draw identifiers.", call. = FALSE)
  }
  components$linear_predictor <-
    components$linear_predictor[draw_ids, , drop = FALSE]
  components$residual_sd <- components$residual_sd[draw_ids]
  components$arm_intercept_sd <- components$arm_intercept_sd[draw_ids]
  components$outcome_intercept_sd <-
    components$outcome_intercept_sd[draw_ids]
  components$study_basis <- lapply(
    components$study_basis,
    function(x) x[draw_ids, , drop = FALSE]
  )
  components$chain_id <- components$chain_id[draw_ids]
  components
}

compute_sensitivity_exposure_integrated_psis_loo <- function(
    job,
    family,
    k_threshold = 0.7,
    cores = 1L
) {
  if (is_sensitivity_pilot_skip(job)) return(job)
  candidate <- job$fitted_candidate
  manifest <- get_primary_loo_observation_manifest(candidate)
  conditional_loo <- job$conditional_loo
  if (
    !inherits(conditional_loo, "loo") ||
      nrow(conditional_loo$pointwise) != nrow(manifest)
  ) {
    stop("The exposure conditional LOO object is invalid.", call. = FALSE)
  }
  draw_ids <- make_sensitivity_loo_draw_ids(
    brms::ndraws(candidate$model),
    job$loo_draws
  )
  components <- get_sensitivity_exposure_integrated_components(
    candidate,
    draw_ids,
    expected_family = family
  )
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
  if (
    !inherits(integrated_loo, "loo") ||
      nrow(integrated_loo$pointwise) != nrow(manifest)
  ) {
    stop("The integrated exposure LOO result is invalid.", call. = FALSE)
  }
  list(
    outcome = candidate$outcome,
    model_index = candidate$model_index,
    model_name = candidate$model_name,
    exposure = family,
    study_structure = components$structure,
    observation_manifest = manifest,
    loo = integrated_loo,
    conditional_loo = conditional_loo,
    remedy = "integrated_observation_level_psis",
    k_threshold = k_threshold,
    comparison = dplyr::bind_rows(
      summarise_primary_loo_method(
        conditional_loo,
        "conditional_observation_level_psis",
        k_threshold
      ),
      summarise_primary_loo_method(
        integrated_loo,
        "integrated_observation_level_psis",
        k_threshold
      )
    ) |>
      dplyr::mutate(
        exposure = family,
        outcome = candidate$outcome,
        model_index = candidate$model_index,
        model_name = candidate$model_name,
        study_structure = components$structure,
        .before = 1L
      )
  )
}

collect_sensitivity_mundlak_integrated_loo <- function(
    results,
    family = "mundlak"
) {
  if (all(purrr::map_lgl(results, is_sensitivity_pilot_skip))) {
    return(new_sensitivity_pilot_skip(family))
  }
  results <- results[order(purrr::map_int(results, "model_index"))]
  if (
    length(results) < 2L ||
      anyDuplicated(purrr::map_int(results, "model_index")) ||
      !all(purrr::map_lgl(results, ~ inherits(.x$loo, "loo")))
  ) {
    stop(
      paste0("The integrated ", family, " LOO set is incomplete."),
      call. = FALSE
    )
  }
  results
}

summarise_sensitivity_mundlak_integrated_loo <- function(results) {
  if (is_sensitivity_pilot_skip(results)) return(results)
  purrr::map_dfr(results, function(result) {
    dplyr::bind_rows(
      summarise_primary_loo_method(
        result$conditional_loo,
        "conditional_observation_level_psis"
      ),
      summarise_primary_loo_method(
        result$loo,
        "integrated_observation_level_psis"
      )
    ) |>
      dplyr::mutate(
        model_index = result$model_index,
        model_name = result$model_name,
        .before = 1L
      )
  })
}

select_sensitivity_mundlak_exact_loo_jobs <- function(
    integrated_results,
    fitted_candidates,
    k_threshold = 0.7,
    family = "mundlak"
) {
  if (
    is_sensitivity_pilot_skip(integrated_results) ||
      is_sensitivity_pilot_skip(fitted_candidates)
  ) {
    return(new_sensitivity_pilot_skip(family))
  }
  fitted_indices <- purrr::map_int(fitted_candidates, "model_index")
  cases <- purrr::map_dfr(integrated_results, function(result) {
    observations <- which(loo::pareto_k_values(result$loo) > k_threshold)
    tibble::tibble(
      model_index = result$model_index,
      observation = observations,
      pareto_k = loo::pareto_k_values(result$loo)[observations]
    )
  }) |>
    dplyr::arrange(model_index, observation) |>
    dplyr::mutate(case_id = dplyr::row_number(), .before = 1L)
  if (nrow(cases) < 1L) return(list())
  purrr::pmap(cases, function(case_id, model_index, observation, pareto_k) {
    fitted_candidate <- fitted_candidates[[match(
      model_index,
      fitted_indices
    )]]
    integrated_result <- integrated_results[[match(
      model_index,
      purrr::map_int(integrated_results, "model_index")
    )]]
    if (
      is.null(fitted_candidate) || is.null(integrated_result) ||
        observation < 1L ||
        observation > nrow(integrated_result$observation_manifest) ||
        integrated_result$loo$diagnostics$pareto_k[[observation]] <=
          k_threshold
    ) {
      stop(
        paste0("A ", family, " exact-LOO job failed validation."),
        call. = FALSE
      )
    }
    list(
      case_id = as.integer(case_id),
      outcome = fitted_candidate$outcome,
      model_index = as.integer(model_index),
      model_name = fitted_candidate$model_name,
      observation = as.integer(observation),
      pareto_k = as.numeric(pareto_k),
      fitted_candidate = fitted_candidate,
      integrated_result = integrated_result
    )
  })
}

fit_sensitivity_mundlak_exact_loo_case <- function(
    job,
    enabled,
    config,
    family = "mundlak"
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip(family))
  original_model <- job$fitted_candidate$model
  original_data <- original_model$data
  manifest <- job$integrated_result$observation_manifest
  if (
    !inherits(original_model, "brmsfit") ||
      nrow(original_data) != nrow(manifest) ||
      job$observation < 1L || job$observation > nrow(original_data)
  ) {
    stop(
      paste0("Invalid ", family, " exact-refit inputs."),
      call. = FALSE
    )
  }
  training_data <- original_data[-job$observation, , drop = FALSE]
  sample_prior <- attr(original_model$prior, "sample_prior")
  original_model$prior <- normalize_primary_brms_priors(
    original_model$prior
  )
  attr(original_model$prior, "sample_prior") <- sample_prior
  specification <- job$fitted_candidate$specification
  start_time <- proc.time()[["elapsed"]]
  refitted_model <- stats::update(
    original_model,
    newdata = training_data,
    recompile = FALSE,
    chains = config$chains,
    cores = config$parallel_chains,
    warmup = config$iter_warmup,
    iter = config$iter_warmup + config$iter_sampling,
    seed = config$seed + job$case_id,
    control = list(
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
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
    refit_seconds = proc.time()[["elapsed"]] - start_time
  )
}

get_sensitivity_mundlak_exact_components <- function(refit) {
  # Achieved- and assignment-intake candidates retain the primary model's
  # heterogeneous study-effect structures (intercepts, linear/spline slopes,
  # or factor smooths). Reconstruct those structures from the Stan design just
  # as in their integrated PSIS calculation, but substitute the complete
  # prediction data into the case-deletion refit first so the held-out row is
  # represented. The legacy Mundlak/within path below instead reconstructs
  # its study basis directly from named data columns.
  if (
    !is.null(refit$family) &&
      unname(as.character(refit$family)) %in% c(
        "achieved", "assignment", "total_achieved", "total_assignment"
      )
  ) {
    prediction_model <- refit$refitted_model
    prediction_model$data <- refit$original_data
    candidate <- list(
      outcome = unname(as.character(refit$outcome)),
      model_index = unname(as.integer(refit$model_index)),
      model_name = unname(as.character(refit$model_name)),
      model = prediction_model
    )
    return(get_primary_integrated_loo_components(candidate))
  }
  model <- refit$refitted_model
  prediction_data <- refit$original_data
  manifest <- refit$integrated_result$observation_manifest
  total_draws <- as.integer(brms::ndraws(model))
  variance_components <- brms::VarCorr(model, summary = FALSE)
  study_sd <- as.matrix(variance_components$study_no$sd)
  coefficient_names <- colnames(study_sd)
  n_coefficients <- ncol(study_sd)
  slope_names <- coefficient_names[coefficient_names != "Intercept"]
  if (
    nrow(prediction_data) != nrow(manifest) ||
      nrow(study_sd) != total_draws ||
      !identical(coefficient_names[[1L]], "Intercept") ||
      !all(slope_names %in% names(prediction_data))
  ) {
    stop("The exact Mundlak prediction design is incomplete.",
         call. = FALSE)
  }
  study_design <- cbind(
    Intercept = 1,
    as.matrix(prediction_data[, slope_names, drop = FALSE])
  )
  uncorrelated_study_effects <- isTRUE(
    refit$specification$uncorrelated_study_effects
  )
  if (n_coefficients > 1L && !uncorrelated_study_effects) {
    covariance_draws <- variance_components$study_no$cov
  }
  if (
    !identical(dim(study_design), c(nrow(manifest), n_coefficients))
  ) {
    stop("The exact Mundlak study design is invalid.", call. = FALSE)
  }
  if (n_coefficients == 1L || uncorrelated_study_effects) {
    study_basis <- lapply(seq_len(n_coefficients), function(component) {
      outer(study_sd[, component], study_design[, component])
    })
  } else {
    if (
      length(dim(covariance_draws)) != 3L ||
        !identical(
          dim(covariance_draws),
          c(total_draws, n_coefficients, n_coefficients)
        )
    ) {
      stop("The exact Mundlak covariance draws are invalid.", call. = FALSE)
    }
    study_basis <- lapply(
      seq_len(n_coefficients),
      function(x) matrix(NA_real_, total_draws, nrow(manifest))
    )
    for (draw in seq_len(total_draws)) {
      covariance <- matrix(
        covariance_draws[draw, , ],
        nrow = n_coefficients,
        ncol = n_coefficients
      )
      cholesky <- tryCatch(t(chol(covariance)), error = function(e) NULL)
      if (is.null(cholesky)) {
        stop(
          "An exact-refit Mundlak covariance draw is not positive definite.",
          call. = FALSE
        )
      }
      transformed_design <- study_design %*% cholesky
      for (component in seq_len(n_coefficients)) {
        study_basis[[component]][draw, ] <-
          transformed_design[, component]
      }
    }
  }
  common_variables <- c(
    residual = "sigma",
    arm = "sd_arm_no__Intercept",
    outcome = "sd_outcome_no__Intercept"
  )
  common_draws <- as.data.frame(extract_primary_exact_draws(
    model,
    common_variables
  ))
  linear_predictor <- brms::posterior_linpred(
    model,
    newdata = prediction_data,
    re_formula = NA,
    transform = FALSE,
    allow_new_levels = TRUE
  )
  list(
    linear_predictor = linear_predictor,
    residual_sd = common_draws[[common_variables[["residual"]]]],
    arm_intercept_sd = common_draws[[common_variables[["arm"]]]],
    outcome_intercept_sd = common_draws[[common_variables[["outcome"]]]],
    study_basis = study_basis,
    chain_id = as.integer(common_draws$.chain)
  )
}

compute_sensitivity_mundlak_exact_loo_case <- function(refit) {
  if (is_sensitivity_pilot_skip(refit)) return(refit)
  integrated_result <- refit$integrated_result
  manifest <- integrated_result$observation_manifest
  observation <- as.integer(refit$observation)
  components <- get_sensitivity_mundlak_exact_components(refit)
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
  if (!is.null(patched_loo$diagnostics$n_eff)) {
    patched_loo$diagnostics$n_eff[[observation]] <- NA_real_
  }
  patched_loo <- update_primary_exact_loo_estimates(patched_loo)
  exact_diagnostics <- summarise_primary_exact_loo_refit_diagnostics(refit) |>
    dplyr::mutate(
      old_psis_pareto_k = old_pareto_k,
      old_psis_elpd_loo = old_elpd,
      exact_elpd_loo = exact_elpd,
      exact_minus_psis_elpd = exact_elpd - old_elpd,
      exact_elpd_mcse = exact_mcse
    )
  list(
    case_id = refit$case_id,
    outcome = refit$outcome,
    model_index = refit$model_index,
    model_name = refit$model_name,
    study_structure = integrated_result$study_structure,
    observation_manifest = manifest,
    loo = patched_loo,
    conditional_loo = integrated_result$conditional_loo,
    remedy = "integrated_observation_level_psis_plus_exact",
    k_threshold = integrated_result$k_threshold,
    comparison = dplyr::bind_rows(
      summarise_primary_loo_method(
        integrated_result$conditional_loo,
        "conditional_observation_level_psis"
      ),
      summarise_primary_loo_method(
        patched_loo,
        "integrated_observation_level_psis_plus_exact"
      )
    ) |>
      dplyr::mutate(
        outcome = refit$outcome,
        model_index = refit$model_index,
        model_name = refit$model_name,
        study_structure = integrated_result$study_structure,
        .before = 1L
      ),
    exact_refit_summary = exact_diagnostics
  )
}

collect_sensitivity_mundlak_exact_results <- function(
    results,
    jobs,
    family = "mundlak"
) {
  if (is_sensitivity_pilot_skip(jobs)) {
    return(new_sensitivity_pilot_skip(family))
  }
  if (length(results) != length(jobs)) {
    stop(
      paste0("The ", family, " exact-result set is incomplete."),
      call. = FALSE
    )
  }
  results <- results[order(purrr::map_int(results, "case_id"))]
  expected <- sort(purrr::map_int(jobs, "case_id"))
  if (!identical(unname(purrr::map_int(results, "case_id")), expected)) {
    stop(
      paste0("The ", family, " exact-result identities do not match."),
      call. = FALSE
    )
  }
  results
}

apply_sensitivity_mundlak_exact_corrections <- function(
    integrated_results,
    exact_results,
    family = "mundlak"
) {
  if (
    is_sensitivity_pilot_skip(integrated_results) ||
      is_sensitivity_pilot_skip(exact_results)
  ) {
    return(new_sensitivity_pilot_skip(family))
  }
  final_results <- integrated_results
  warning_cases <- purrr::map_dfr(integrated_results, function(result) {
    observations <- which(loo::pareto_k_values(result$loo) > 0.7)
    tibble::tibble(
      model_index = result$model_index,
      observation = observations
    )
  }) |>
    dplyr::arrange(model_index, observation)
  exact_cases <- tibble::tibble(
    model_index = unname(purrr::map_int(exact_results, "model_index")),
    observation = unname(purrr::map_int(exact_results, ~
      .x$exact_refit_summary$observation[[1L]]))
  ) |>
    dplyr::arrange(model_index, observation)
  if (!identical(warning_cases, exact_cases)) {
    stop(
      paste0(
        "Exact ", family,
        " corrections do not match the high-k cases."
      ),
      call. = FALSE
    )
  }
  corrected_models <- unique(unname(purrr::map_int(
    exact_results,
    "model_index"
  )))
  for (model_index in corrected_models) {
    position <- match(
      model_index,
      purrr::map_int(final_results, "model_index")
    )
    patched_result <- final_results[[position]]
    model_corrections <- purrr::keep(
      exact_results,
      ~ identical(as.integer(.x$model_index), as.integer(model_index))
    )
    for (exact_result in model_corrections) {
      observation <- as.integer(
        exact_result$exact_refit_summary$observation[[1L]]
      )
      patched_result$loo$pointwise[observation, ] <-
        exact_result$loo$pointwise[observation, ]
      patched_result$loo$diagnostics$pareto_k[[observation]] <- NA_real_
      if (!is.null(patched_result$loo$diagnostics$n_eff)) {
        patched_result$loo$diagnostics$n_eff[[observation]] <- NA_real_
      }
    }
    patched_result$loo <- update_primary_exact_loo_estimates(
      patched_result$loo
    )
    patched_result$remedy <-
      "integrated_observation_level_psis_plus_exact"
    patched_result$comparison <- dplyr::bind_rows(
      summarise_primary_loo_method(
        patched_result$conditional_loo,
        "conditional_observation_level_psis"
      ),
      summarise_primary_loo_method(
        patched_result$loo,
        "integrated_observation_level_psis_plus_exact"
      )
    ) |>
      dplyr::mutate(
        outcome = patched_result$outcome,
        model_index = patched_result$model_index,
        model_name = patched_result$model_name,
        study_structure = patched_result$study_structure,
        .before = 1L
      )
    final_results[[position]] <- patched_result
  }
  exact_refit_summary <- purrr::map_dfr(
    exact_results,
    "exact_refit_summary"
  )
  model_summary <- summarise_sensitivity_mundlak_integrated_loo(final_results) |>
    dplyr::filter(stringr::str_starts(method, "integrated_"))
  stacking_ready <-
    all(model_summary$n_pareto_k_gt_threshold == 0L) &&
    all(exact_refit_summary$diagnostic_pass)
  list(
    results = final_results,
    model_summary = model_summary,
    exact_refit_summary = exact_refit_summary,
    stacking_ready = stacking_ready
  )
}

compute_sensitivity_integrated_weights <- function(
    results,
    family = "mundlak"
) {
  if (is_sensitivity_pilot_skip(results)) return(results)
  if (!is.null(results$results)) {
    if (!isTRUE(results$stacking_ready)) {
      stop(
        paste0(
          "The final integrated ", family,
          " LOO audit is not stacking-ready."
        ),
        call. = FALSE
      )
    }
    results <- results$results
  }
  high_k <- purrr::map_int(
    results,
    ~ sum(loo::pareto_k_values(.x$loo) > 0.7, na.rm = TRUE)
  )
  if (any(high_k > 0L)) {
    stop(
      "Integrated ", family, " PSIS still has high Pareto-k cases; ",
      "exact corrections are required before stacking.",
      call. = FALSE
    )
  }
  loo_objects <- purrr::map(results, "loo")
  names(loo_objects) <- purrr::map_chr(results, "model_name")
  weights <- loo::loo_model_weights(loo_objects, method = "stacking")
  weights <- pmax(as.numeric(weights), 0)
  weights <- weights / sum(weights)
  names(weights) <- names(loo_objects)
  list(
    loo = loo_objects,
    weights = weights,
    method = "integrated_observation_level_psis"
  )
}

make_mundlak_prediction_newdata <- function(
    fitted_candidate,
    data,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.02,
    study_mean_raw = NULL
) {
  grid <- make_primary_prediction_grid(
    data, protein_min, protein_max, grid_step
  )
  basis_spec <- fitted_candidate$basis_spec
  if (is.null(basis_spec)) return(grid)
  if (is.null(study_mean_raw)) study_mean_raw <- grid$protein_median
  varied <- predict_sensitivity_basis(
    basis_spec,
    grid$newdata$protein_intake
  )
  held <- predict_sensitivity_basis(
    basis_spec,
    study_mean_raw - grid$protein_median
  )
  held <- matrix(held, nrow = nrow(varied), ncol = ncol(varied), byrow = TRUE)
  for (index in seq_len(ncol(varied))) {
    grid$newdata[[fitted_candidate$within_names[[index]]]] <-
      varied[, index] - held[, index]
    grid$newdata[[fitted_candidate$between_names[[index]]]] <- held[, index]
  }
  grid
}

make_within_prediction_newdata <- function(
    fitted_candidate,
    data,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.02,
    reference_raw = NULL
) {
  observed_max <- max(
    c(data$protein_intake_raw_treatment,
      data$protein_intake_raw_comparator),
    na.rm = TRUE
  )
  if (is.null(protein_max)) protein_max <- observed_max
  if (is.null(reference_raw)) reference_raw <- protein_min
  if (
    length(reference_raw) != 1L || !is.finite(reference_raw) ||
      reference_raw < protein_min || reference_raw > protein_max
  ) {
    stop("The within-study prediction reference is outside the grid.",
         call. = FALSE)
  }
  protein_values <- sort(unique(c(
    round(seq(protein_min, protein_max, by = grid_step), 10),
    reference_raw
  )))
  keys <- tidyr::crossing(
    resistance_exercise_code = c(0, 1),
    protein_intake = protein_values
  )
  newdata <- fitted_candidate$model$data[
    rep(1L, nrow(keys)), , drop = FALSE
  ]
  newdata$resistance_exercise_code <- keys$resistance_exercise_code
  if ("resistance_exercise_factor" %in% names(newdata)) {
    labels <- ifelse(keys$resistance_exercise_code == 0, "yes", "no")
    newdata$resistance_exercise_factor <- factor(
      labels,
      levels = levels(fitted_candidate$model$data$resistance_exercise_factor)
    )
  }
  if (!is.null(fitted_candidate$basis_spec)) {
    internal_unit <- get_primary_protein_internal_unit_for_contrasts(data)
    centring_constant <- get_primary_protein_median_for_contrasts(data)
    varied <- predict_sensitivity_basis(
      fitted_candidate$basis_spec,
      (protein_values - centring_constant) / internal_unit
    )
    varied <- varied[
      match(keys$protein_intake, protein_values), , drop = FALSE
    ]
    held <- predict_sensitivity_basis(
      fitted_candidate$basis_spec,
      (reference_raw - centring_constant) / internal_unit
    )
    held <- matrix(
      held, nrow = nrow(varied), ncol = ncol(varied), byrow = TRUE
    )
    for (index in seq_len(ncol(varied))) {
      newdata[[fitted_candidate$contrast_names[[index]]]] <-
        varied[, index] - held[, index]
    }
  }
  if ("age_centre" %in% names(newdata)) newdata$age_centre <- 0
  if ("duration_centre" %in% names(newdata)) newdata$duration_centre <- 0
  if ("age_outcome_z" %in% names(newdata)) {
    moderator_data <- add_primary_standardised_moderators(data)
    centres <- attr(moderator_data, "moderator_centres")
    scales <- attr(moderator_data, "moderator_scales")
    newdata$age_outcome_z <-
      (0 - centres[["age_centre"]]) / scales[["age_centre"]]
    newdata$duration_outcome_z <-
      (0 - centres[["duration_centre"]]) /
        scales[["duration_centre"]]
  }
  list(
    newdata = newdata,
    keys = tibble::tibble(
      resistance_exercise_code = keys$resistance_exercise_code,
      resistance_exercise = ifelse(
        keys$resistance_exercise_code == 0, "YES", "NO"
      ),
      protein_intake = keys$protein_intake,
      reference_protein_intake = reference_raw
    ),
    protein_min = protein_min,
    protein_max = tail(protein_values, 1L),
    grid_step = grid_step
  )
}

get_primary_protein_internal_unit_for_contrasts <- function(data) {
  if (!"absolute_exposure_unit_g_day" %in% names(data)) return(1)
  units <- unique(as.numeric(data$absolute_exposure_unit_g_day))
  units <- units[is.finite(units)]
  if (length(units) != 1L || units[[1L]] <= 0) {
    stop("Within-study contrasts do not share one internal exposure unit.",
         call. = FALSE)
  }
  units[[1L]]
}

get_primary_protein_median_for_contrasts <- function(data) {
  internal_unit <- get_primary_protein_internal_unit_for_contrasts(data)
  values <- unique(c(
    data$protein_intake_raw_treatment -
      data$protein_intake_treatment * internal_unit,
    data$protein_intake_raw_comparator -
      data$protein_intake_comparator * internal_unit
  ))
  values <- values[is.finite(values)]
  centre <- stats::median(values)
  tolerance <- 1e-10 * max(1, abs(centre))
  if (
    length(values) < 1L || !is.finite(centre) ||
      any(abs(values - centre) > tolerance)
  ) {
    stop("Within-study contrasts do not share one centring constant.",
         call. = FALSE)
  }
  if ("absolute_exposure_reference_g_day" %in% names(data)) {
    references <- unique(as.numeric(data$absolute_exposure_reference_g_day))
    references <- references[is.finite(references)]
    if (
      length(references) != 1L ||
        abs(references[[1L]] - centre) > tolerance
    ) {
      stop("The declared within-study centring constant is inconsistent.",
           call. = FALSE)
    }
    return(references[[1L]])
  }
  centre
}

make_sensitivity_prediction_grid <- function(
    fitted_candidate,
    data,
    prediction_type,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.02
) {
  switch(
    prediction_type,
    exposure = make_primary_prediction_grid(
      data, protein_min, protein_max, grid_step
    ),
    mundlak = make_mundlak_prediction_newdata(
      fitted_candidate, data, protein_min, protein_max, grid_step
    ),
    within = make_within_prediction_newdata(
      fitted_candidate, data, protein_min, protein_max, grid_step
    ),
    strategy = make_protein_strategy_prediction_grid(
      fitted_candidate, data, protein_min, protein_max, grid_step
    ),
    strategy_within = make_protein_strategy_within_prediction_grid(
      fitted_candidate, data, protein_min, protein_max, grid_step
    ),
    stop("Unknown sensitivity prediction type.", call. = FALSE)
  )
}

build_sensitivity_brms_model_average <- function(
    fits,
    weight_result,
    data,
    prediction_type,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.02,
    n_draws = 1000L,
    seed = 1988L
) {
  if (is_sensitivity_pilot_skip(fits)) return(fits)
  weights <- weight_result$weights
  models <- purrr::map(fits, "model")
  names(models) <- purrr::map_chr(fits, "model_name")
  counts <- allocate_primary_mixture_models(weights, n_draws)
  set.seed(seed)
  allocation <- sample(rep(seq_along(fits), times = counts))

  grids <- purrr::map(
    fits,
    ~ make_sensitivity_prediction_grid(
      .x,
      data,
      prediction_type,
      protein_min = protein_min,
      protein_max = protein_max,
      grid_step = grid_step
    )
  )
  key_reference <- grids[[1L]]$keys
  if (!all(purrr::map_lgl(grids, ~ identical(.x$keys, key_reference)))) {
    stop("Sensitivity prediction grids do not match across models.",
         call. = FALSE)
  }
  draws <- matrix(NA_real_, nrow = n_draws, ncol = nrow(key_reference))
  for (index in seq_along(fits)) {
    rows <- which(allocation == index)
    if (!length(rows)) next
    model_draws <- brms::posterior_epred(
      fits[[index]]$model,
      newdata = grids[[index]]$newdata,
      re_formula = NA,
      ndraws = length(rows)
    )
    draws[rows, ] <- as.matrix(model_draws)
  }
  output <- list(
    draws = draws,
    grid = key_reference,
    model_weights = weights,
    draw_manifest = tibble::tibble(
      drawid = seq_len(n_draws),
      model_index = purrr::map_int(fits, "model_index")[allocation],
      model_name = names(models)[allocation]
    ),
    settings = list(
      n_draws = n_draws,
      seed = seed,
      prediction_type = prediction_type,
      protein_min = grids[[1L]]$protein_min,
      protein_max = grids[[1L]]$protein_max,
      grid_step = grids[[1L]]$grid_step
    )
  )
  class(output) <- c("primary_model_average", class(output))
  output
}

summarise_sensitivity_model_average <- function(model_average) {
  if (is_sensitivity_pilot_skip(model_average)) return(model_average)
  summarise_primary_dose_response(model_average)
}

summarise_sensitivity_exposure_support <- function(data) {
  raw_intake <- if (
    "protein_intake_raw" %in% names(data)
  ) {
    data$protein_intake_raw
  } else {
    data$protein_intake + data$protein_median
  }
  arms <- data |>
    dplyr::mutate(protein_intake_support = raw_intake) |>
    dplyr::distinct(
      study_no,
      arm_no,
      resistance_exercise_code,
      protein_intake_support
    )
  if (
    nrow(arms) != dplyr::n_distinct(
      interaction(arms$study_no, arms$arm_no, drop = TRUE)
    ) ||
      anyNA(arms$protein_intake_support)
  ) {
    stop("The sensitivity arm-level exposure support is invalid.",
         call. = FALSE)
  }
  by_group <- arms |>
    dplyr::group_by(resistance_exercise_code) |>
    dplyr::summarise(
      arms = dplyr::n(),
      studies = dplyr::n_distinct(study_no),
      observed_minimum = min(protein_intake_support),
      central_90_minimum = as.numeric(stats::quantile(
        protein_intake_support, 0.05
      )),
      median = stats::median(protein_intake_support),
      central_90_maximum = as.numeric(stats::quantile(
        protein_intake_support, 0.95
      )),
      observed_maximum = max(protein_intake_support),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      resistance_exercise = ifelse(
        resistance_exercise_code == 0,
        "YES",
        "NO"
      ),
      .after = resistance_exercise_code
    )
  if (nrow(by_group) != 2L) {
    stop("Expected both resistance-exercise support groups.",
         call. = FALSE)
  }
  common <- tibble::tibble(
    common_observed_minimum = max(by_group$observed_minimum),
    common_observed_maximum = min(by_group$observed_maximum),
    common_central_90_minimum = max(by_group$central_90_minimum),
    common_central_90_maximum = min(by_group$central_90_maximum)
  )
  if (
    common$common_observed_minimum > common$common_observed_maximum ||
      common$common_central_90_minimum > common$common_central_90_maximum
  ) {
    stop("Resistance-exercise exposure supports do not overlap.",
         call. = FALSE)
  }
  list(by_group = by_group, common = common, arms = arms)
}

summarise_sensitivity_within_support <- function(
    data,
    reference_intake = 0.8
) {
  pairs <- data |>
    dplyr::distinct(
      study_no,
      protein_arm_no,
      comparator_arm_no,
      resistance_exercise_code,
      protein_intake_raw_treatment,
      protein_intake_raw_comparator
    ) |>
    dplyr::mutate(
      dose_difference =
        protein_intake_raw_treatment - protein_intake_raw_comparator
    )
  required <- c(
    "protein_intake_raw_treatment",
    "protein_intake_raw_comparator",
    "dose_difference"
  )
  if (nrow(pairs) < 1L || anyNA(pairs[, required, drop = FALSE])) {
    stop("The within-contrast exposure support is invalid.", call. = FALSE)
  }
  by_group <- pairs |>
    dplyr::group_by(resistance_exercise_code) |>
    dplyr::summarise(
      comparisons = dplyr::n(),
      studies = dplyr::n_distinct(study_no),
      observed_minimum = min(protein_intake_raw_treatment),
      central_90_minimum = as.numeric(stats::quantile(
        protein_intake_raw_treatment, 0.05
      )),
      median = stats::median(protein_intake_raw_treatment),
      central_90_maximum = as.numeric(stats::quantile(
        protein_intake_raw_treatment, 0.95
      )),
      observed_maximum = max(protein_intake_raw_treatment),
      comparator_observed_minimum = min(protein_intake_raw_comparator),
      comparator_central_90_minimum = as.numeric(stats::quantile(
        protein_intake_raw_comparator, 0.05
      )),
      comparator_median = stats::median(protein_intake_raw_comparator),
      comparator_central_90_maximum = as.numeric(stats::quantile(
        protein_intake_raw_comparator, 0.95
      )),
      comparator_observed_maximum = max(protein_intake_raw_comparator),
      dose_difference_observed_minimum = min(dose_difference),
      dose_difference_central_90_minimum = as.numeric(stats::quantile(
        dose_difference, 0.05
      )),
      dose_difference_median = stats::median(dose_difference),
      dose_difference_central_90_maximum = as.numeric(stats::quantile(
        dose_difference, 0.95
      )),
      dose_difference_observed_maximum = max(dose_difference),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      resistance_exercise = ifelse(
        resistance_exercise_code == 0, "YES", "NO"
      ),
      reference_intake = reference_intake,
      reference_within_comparator_observed_support =
        reference_intake >= comparator_observed_minimum &
          reference_intake <= comparator_observed_maximum,
      reference_within_comparator_central_90_support =
        reference_intake >= comparator_central_90_minimum &
          reference_intake <= comparator_central_90_maximum,
      .after = resistance_exercise_code
    )
  if (nrow(by_group) != 2L) {
    stop("Expected both within-contrast exercise groups.", call. = FALSE)
  }
  common <- tibble::tibble(
    common_observed_minimum = max(by_group$observed_minimum),
    common_observed_maximum = min(by_group$observed_maximum),
    common_central_90_minimum = max(by_group$central_90_minimum),
    common_central_90_maximum = min(by_group$central_90_maximum)
  )
  if (
    common$common_observed_minimum > common$common_observed_maximum ||
      common$common_central_90_minimum > common$common_central_90_maximum
  ) {
    stop("Within-contrast treatment supports do not overlap.", call. = FALSE)
  }
  list(by_group = by_group, common = common, pairs = pairs)
}

annotate_sensitivity_curve_support <- function(curve, support) {
  common <- support$common
  support_by_group <- support$by_group |>
    dplyr::select(-dplyr::any_of(c(
      "resistance_exercise", "reference_intake"
    )))
  output <- curve |>
    dplyr::left_join(support_by_group, by = "resistance_exercise_code") |>
    dplyr::mutate(
      within_group_observed_support =
        protein_intake >= observed_minimum &
          protein_intake <= observed_maximum,
      within_group_central_90_support =
        protein_intake >= central_90_minimum &
          protein_intake <= central_90_maximum,
      within_common_observed_support =
        protein_intake >= common$common_observed_minimum &
          protein_intake <= common$common_observed_maximum,
      within_common_central_90_support =
        protein_intake >= common$common_central_90_minimum &
          protein_intake <= common$common_central_90_maximum,
      support_region = dplyr::case_when(
        within_common_central_90_support ~ "common_central_90",
        within_common_observed_support ~ "common_observed_tail",
        within_group_central_90_support ~ "group_central_90_only",
        within_group_observed_support ~ "group_observed_tail",
        TRUE ~ "extrapolation"
      )
    )
  if (anyNA(output$support_region)) {
    stop("Sensitivity support annotation failed.", call. = FALSE)
  }
  output
}

summarise_sensitivity_model_average_increments <- function(
    model_average,
    support,
    reference_intake = 0.8
) {
  if (is_sensitivity_pilot_skip(model_average)) return(model_average)
  grid <- model_average$grid
  increment_draws <- model_average$draws
  groups <- sort(unique(grid$resistance_exercise_code))
  for (group in groups) {
    group_rows <- which(grid$resistance_exercise_code == group)
    reference_rows <- group_rows[
      abs(grid$protein_intake[group_rows] - reference_intake) < 1e-10
    ]
    if (length(reference_rows) != 1L) {
      stop("The sensitivity increment reference is not on the grid.",
           call. = FALSE)
    }
    increment_draws[, group_rows] <- sweep(
      increment_draws[, group_rows, drop = FALSE],
      1L,
      model_average$draws[, reference_rows],
      FUN = "-"
    )
  }
  increment_average <- model_average
  increment_average$draws <- increment_draws
  output <- summarise_primary_dose_response(increment_average) |>
    dplyr::rename(
      increment = estimate,
      increment_conf_low = conf.low,
      increment_conf_high = conf.high
    ) |>
    dplyr::mutate(reference_intake = reference_intake)
  annotate_sensitivity_curve_support(output, support)
}

build_sensitivity_mundlak_reporting_pilot <- function(
    absolute_summary,
    model_average,
    data,
    reference_intake = 0.8
) {
  support <- summarise_sensitivity_exposure_support(data)
  absolute <- annotate_sensitivity_curve_support(
    absolute_summary,
    support
  )
  increments <- summarise_sensitivity_model_average_increments(
    model_average,
    support,
    reference_intake = reference_intake
  )
  list(
    support_by_group = support$by_group,
    common_support = support$common,
    absolute_full_grid = absolute,
    increments_full_grid = increments,
    common_support_absolute = dplyr::filter(
      absolute,
      within_common_observed_support
    ),
    common_support_increments = dplyr::filter(
      increments,
      within_common_observed_support
    ),
    group_supported_absolute = dplyr::filter(
      absolute,
      within_group_observed_support
    ),
    group_supported_increments = dplyr::filter(
      increments,
      within_group_observed_support
    ),
    reference_intake = reference_intake
  )
}

build_sensitivity_within_reporting_pilot <- function(
    absolute_summary,
    model_average,
    data,
    reference_intake = 0.8
) {
  support <- summarise_sensitivity_within_support(
    data,
    reference_intake = reference_intake
  )
  absolute <- annotate_sensitivity_curve_support(absolute_summary, support)
  increments <- summarise_sensitivity_model_average_increments(
    model_average,
    support,
    reference_intake = reference_intake
  )
  list(
    support_by_group = support$by_group,
    common_support = support$common,
    absolute_full_grid = absolute,
    increments_full_grid = increments,
    common_support_absolute = dplyr::filter(
      absolute,
      within_common_observed_support
    ),
    common_support_increments = dplyr::filter(
      increments,
      within_common_observed_support
    ),
    group_supported_absolute = dplyr::filter(
      absolute,
      within_group_observed_support
    ),
    group_supported_increments = dplyr::filter(
      increments,
      within_group_observed_support
    ),
    reference_intake = reference_intake
  )
}

compare_sensitivity_exposure_curves <- function(
    achieved_summary,
    assignment_summary
) {
  if (
    is_sensitivity_pilot_skip(achieved_summary) ||
      is_sensitivity_pilot_skip(assignment_summary)
  ) {
    return(new_sensitivity_pilot_skip("exposure"))
  }
  keys <- intersect(
    c("resistance_exercise_code", "resistance_exercise_label",
      "resistance_exercise", "protein_intake"),
    intersect(names(achieved_summary), names(assignment_summary))
  )
  dplyr::full_join(
    dplyr::rename(
      achieved_summary,
      achieved_estimate = estimate,
      achieved_conf_low = conf.low,
      achieved_conf_high = conf.high
    ),
    dplyr::rename(
      assignment_summary,
      assignment_estimate = estimate,
      assignment_conf_low = conf.low,
      assignment_conf_high = conf.high
    ),
    by = keys
  ) |>
    dplyr::mutate(
      assignment_minus_achieved =
        assignment_estimate - achieved_estimate
    )
}

make_rcs_orthogonal_spec <- function(
    x_obs,
    probs = c(0.05, 0.35, 0.65, 0.95),
    distinct_fallback = c("error", "unique_support")
) {
  distinct_fallback <- match.arg(distinct_fallback)
  knots <- as.numeric(stats::quantile(
    x_obs, probs = probs, names = FALSE
  ))
  if (
    all(is.finite(knots)) && any(diff(knots) <= 0) &&
      identical(distinct_fallback, "unique_support")
  ) {
    support <- sort(unique(x_obs[is.finite(x_obs)]))
    if (length(support) < length(probs)) {
      stop(
        "Restricted cubic splines require at least four distinct exposures.",
        call. = FALSE
      )
    }
    target <- knots
    n_knot <- length(target)
    n_support <- length(support)
    scale <- diff(range(support))
    cost <- outer(
      target,
      support,
      function(target_value, support_value) {
        ((target_value - support_value) / scale)^2
      }
    )
    best <- matrix(Inf, nrow = n_knot, ncol = n_support)
    previous <- matrix(NA_integer_, nrow = n_knot, ncol = n_support)
    best[1L, seq_len(n_support - n_knot + 1L)] <-
      cost[1L, seq_len(n_support - n_knot + 1L)]
    for (knot_index in 2:n_knot) {
      last_support <- n_support - n_knot + knot_index
      for (support_index in knot_index:last_support) {
        candidates <- best[knot_index - 1L, seq_len(support_index - 1L)]
        previous_index <- which.min(candidates)
        best[knot_index, support_index] <-
          candidates[[previous_index]] + cost[knot_index, support_index]
        previous[knot_index, support_index] <- previous_index
      }
    }
    selected <- integer(n_knot)
    selected[[n_knot]] <- which.min(best[n_knot, ])
    for (knot_index in seq.int(n_knot, 2L)) {
      selected[[knot_index - 1L]] <-
        previous[knot_index, selected[[knot_index]]]
    }
    knots <- support[selected]
  }
  if (any(!is.finite(knots)) || any(diff(knots) <= 0)) {
    stop("Restricted cubic-spline knots are not distinct.", call. = FALSE)
  }
  raw_basis <- function(x) {
    d <- vapply(
      knots[1:3],
      function(knot) {
        (pmax(x - knot, 0)^3 - pmax(x - knots[4], 0)^3) /
          (knots[4] - knot)
      },
      numeric(length(x))
    )
    sweep(d[, 1:2, drop = FALSE], 1L, d[, 3L], "-")
  }
  h <- raw_basis(x_obs)
  linear_design <- cbind(1, x_obs)
  projection <- qr.solve(linear_design, h)
  residual <- h - linear_design %*% projection
  qr_h <- qr(residual)
  if (qr_h$rank < 2L) {
    stop("Restricted cubic-spline basis is rank deficient.", call. = FALSE)
  }
  transform <- solve(qr.R(qr_h)) * sqrt(nrow(residual))
  list(
    rcs_knots = knots,
    rcs_projection = projection,
    rcs_basis_transform = transform
  )
}

make_eiv_candidate_data <- function(data, specification, priors) {
  features <- sensitivity_model_features(specification$model_index)
  required <- c(
    "yi", "vi", "protein_intake", "protein_intake_se",
    "resistance_exercise_code", "study_no", "arm_no", "outcome_no",
    "age_centre", "duration_centre", "protein_intake_raw",
    "protein_median"
  )
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("EIV data are missing: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  data <- data |>
    dplyr::filter(
      is.finite(yi), is.finite(vi), vi > 0,
      is.finite(protein_intake), is.finite(protein_intake_se),
      protein_intake_se > 0
    )
  moderator_data <- add_primary_standardised_moderators(data)
  moderator_centres <- attr(moderator_data, "moderator_centres")
  moderator_scales <- attr(moderator_data, "moderator_scales")
  data <- moderator_data |>
    dplyr::arrange(study_no, arm_no, outcome_no) |>
    dplyr::mutate(
      study_index = dplyr::dense_rank(study_no),
      arm_index = dplyr::dense_rank(interaction(study_no, arm_no,
                                                drop = TRUE)),
      outcome_index = dplyr::dense_rank(interaction(
        study_no, arm_no, outcome_no, drop = TRUE
      ))
    )
  arm_data <- data |>
    dplyr::distinct(
      study_no, arm_no, arm_index, protein_intake, protein_intake_se
    ) |>
    dplyr::arrange(arm_index)
  arm_check <- arm_data |>
    dplyr::count(study_no, arm_no) |>
    dplyr::filter(n != 1L)
  if (nrow(arm_check) > 0L) {
    stop("EIV study-arms do not have unique exposure summaries.",
         call. = FALSE)
  }

  basis_type <- switch(
    features$protein_family,
    none = 0L,
    linear = 1L,
    piecewise = 2L,
    spline = 3L
  )
  k_basis <- switch(
    as.character(basis_type),
    `0` = 1L, `1` = 1L, `2` = 2L, `3` = 3L
  )
  x_scale <- if (
    features$protein_family == "spline" &&
      isTRUE(specification$half_unit_protein)
  ) 0.5 else 1
  if (basis_type == 3L) {
    rcs_probs <- if (isTRUE(specification$natural_spline)) {
      c(0, 1 / 3, 2 / 3, 1)
    } else {
      c(0.05, 0.35, 0.65, 0.95)
    }
    rcs <- make_rcs_orthogonal_spec(
      arm_data$protein_intake / x_scale,
      probs = rcs_probs
    )
  } else {
    rcs <- list(
      rcs_knots = c(-2, -1, 1, 2),
      rcs_projection = matrix(0, 2, 2),
      rcs_basis_transform = diag(2)
    )
  }

  protein_median <- get_primary_protein_median(data)
  protein_grid_raw <- seq(
    0.8,
    max(data$protein_intake_raw, na.rm = TRUE),
    by = 0.02
  )
  grid <- tidyr::crossing(
    resistance_exercise_code = c(0, 1),
    protein_intake_raw = protein_grid_raw
  ) |>
    dplyr::mutate(protein_intake = protein_intake_raw - protein_median)
  age_grid_z <- (0 - moderator_centres[["age_centre"]]) /
    moderator_scales[["age_centre"]]
  duration_grid_z <- (0 - moderator_centres[["duration_centre"]]) /
    moderator_scales[["duration_centre"]]
  stan_priors <- make_stan_priors(priors)
  k_study <- 1L + as.integer(features$includes_random_slope) * k_basis
  coefficient_scales <- rep(1, k_basis)
  if (basis_type == 3L && !isTRUE(specification$natural_spline)) {
    coefficient_scales[-1L] <- 0.5
  }
  variance_scale <- if (
    isTRUE(specification$regularised_variance_priors)
  ) 1 else 2.5
  study_sd_scales <- rep(variance_scale, k_study)
  if (features$includes_random_slope && basis_type == 3L) {
    study_sd_scales[-1L] <- 0.5
  }
  if (
    features$includes_random_slope && basis_type == 2L &&
      (isTRUE(specification$targeted_piecewise_slope_priors) ||
        isTRUE(specification$targeted_function_tail_priors))
  ) {
    study_sd_scales[-1L] <- 0.5
    if (isTRUE(specification$targeted_function_tail_priors)) {
      study_sd_scales[[k_study]] <- 0.25
    }
  }

  list(
    family = "eiv",
    specification = specification,
    grid = grid,
    stan_data = list(
      N = nrow(data),
      y = data$yi,
      se_y = sqrt(data$vi),
      N_arm = nrow(arm_data),
      arm = data$arm_index,
      x_obs = arm_data$protein_intake,
      x_se = arm_data$protein_intake_se,
      N_study = dplyr::n_distinct(data$study_index),
      study = data$study_index,
      N_outcome = dplyr::n_distinct(data$outcome_index),
      outcome = data$outcome_index,
      resistance_exercise_code = data$resistance_exercise_code,
      age_z = data$age_outcome_z,
      duration_z = data$duration_outcome_z,
      basis_type = basis_type,
      K_basis = k_basis,
      include_protein = as.integer(features$includes_protein),
      protein_scale = x_scale,
      protein_knot = 0.5,
      rcs_knots = rcs$rcs_knots,
      rcs_projection = rcs$rcs_projection,
      rcs_basis_transform = rcs$rcs_basis_transform,
      include_moderators = as.integer(features$includes_moderators),
      include_rt_interaction =
        as.integer(features$includes_rt_interaction),
      include_random_slope =
        as.integer(features$includes_random_slope),
      K_study = k_study,
      correlated_study_effects =
        as.integer(!isTRUE(specification$uncorrelated_study_effects)),
      prior_intercept_location = stan_priors$prior_intercept_location,
      prior_intercept_scale = stan_priors$prior_intercept_scale,
      prior_rt_location = stan_priors$prior_rt_location,
      prior_rt_scale = stan_priors$prior_rt_scale,
      prior_protein_scale = coefficient_scales,
      prior_study_sd_scale = study_sd_scales,
      prior_arm_sd_scale = variance_scale,
      prior_outcome_sd_scale = variance_scale,
      N_grid = nrow(grid),
      x_grid = grid$protein_intake,
      rt_grid = grid$resistance_exercise_code,
      age_grid_z = rep(age_grid_z, nrow(grid)),
      duration_grid_z = rep(duration_grid_z, nrow(grid))
    )
  )
}

compile_sensitivity_stan_model <- function(path, enabled, family) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip(family))
  cmdstanr::cmdstan_model(path)
}

make_sensitivity_cmdstan_output_dir <- function(
    family,
    outcome,
    model_index
) {
  path_component <- function(x) {
    gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  }
  output_dir <- here::here(
    "_targets",
    "user",
    "sensitivity_cmdstan",
    path_component(family),
    path_component(outcome),
    paste0("model_", as.integer(model_index))
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(output_dir, winslash = "/", mustWork = TRUE)
}

fit_eiv_cmdstan_pilot <- function(input, model, enabled, config) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("eiv"))
  if (is_sensitivity_pilot_skip(model)) {
    stop("The EIV Stan model was not compiled.", call. = FALSE)
  }
  specification <- input$specification
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  output_dir <- make_sensitivity_cmdstan_output_dir(
    "eiv",
    input$specification$outcome,
    input$specification$model_index
  )
  fit <- model$sample(
    data = input$stan_data,
    chains = config$chains,
    parallel_chains = config$parallel_chains,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    adapt_delta = specification$adapt_delta,
    max_treedepth = specification$max_treedepth,
    seed = config$seed,
    refresh = 100,
    output_dir = output_dir
  )
  list(
    family = "eiv",
    outcome = input$specification$outcome,
    model_index = input$specification$model_index,
    model_name = paste0(input$specification$model_name, " [EIV Stan]"),
    specification = specification,
    grid = input$grid,
    output_dir = output_dir,
    fit = fit
  )
}

diagnose_cmdstan_pilot_fit <- function(candidate) {
  if (is_sensitivity_pilot_skip(candidate)) return(candidate)
  sampler <- posterior::as_draws_df(candidate$fit$sampler_diagnostics())
  parameter_summary <- candidate$fit$summary() |>
    dplyr::filter(!stringr::str_detect(
      variable,
      paste0(
        "^(log_lik_y|log_lik_cluster|y_rep|mu\\[|mu_grid|x_true|",
        "study_effect|study_correlation)"
      )
    ))
  features <- sensitivity_model_features(candidate$model_index)
  if (!features$includes_protein) {
    parameter_summary <- dplyr::filter(
      parameter_summary,
      !stringr::str_detect(variable, "^b_protein($|\\[)")
    )
  }
  if (!features$includes_rt_interaction) {
    parameter_summary <- dplyr::filter(
      parameter_summary,
      !stringr::str_detect(variable, "^b_protein_rt($|\\[)")
    )
  }
  if (!features$includes_moderators) {
    parameter_summary <- dplyr::filter(
      parameter_summary,
      !stringr::str_detect(variable, "^b_moderator($|\\[)")
    )
  }
  if (
    candidate$family %in% c(
      "eiv", "rescor", "total_rescor", "supplement_rescor"
    ) &&
      isTRUE(candidate$specification$uncorrelated_study_effects)
  ) {
    parameter_summary <- dplyr::filter(
      parameter_summary,
      !stringr::str_detect(variable, "^L_study($|\\[)")
    )
  }
  # Cholesky factors of dimension one and other structural constants have no
  # Monte Carlo error and should not enter convergence minima or maxima.
  parameter_summary <- parameter_summary |>
    dplyr::filter(is.finite(sd), sd > 0)
  if ("mcse_mean" %in% names(parameter_summary)) {
    parameter_summary <- parameter_summary |>
      dplyr::mutate(
        mcse_sd_ratio = dplyr::if_else(
          is.finite(sd) & sd > 0,
          mcse_mean / sd,
          0
        )
      )
  } else {
    # CmdStanR's default summary in some supported versions omits MCSE.
    # For the mean, MCSE / posterior SD is 1 / sqrt(bulk ESS).
    parameter_summary <- parameter_summary |>
      dplyr::mutate(
        mcse_sd_ratio = dplyr::if_else(
          is.finite(ess_bulk) & ess_bulk > 0,
          1 / sqrt(ess_bulk),
          Inf
        )
      )
  }
  energy_by_chain <- sampler |>
    dplyr::group_by(.chain) |>
    dplyr::summarise(
      ebfmi = mean(diff(energy__)^2) / stats::var(energy__),
      .groups = "drop"
    )
  divergent_transitions <- sum(sampler$divergent__)
  maximum_treedepth_transitions <- sum(
    sampler$treedepth__ >= candidate$specification$max_treedepth
  )
  minimum_ebfmi <- min(energy_by_chain$ebfmi)
  maximum_rhat <- max(parameter_summary$rhat, na.rm = TRUE)
  minimum_bulk_ess <- min(parameter_summary$ess_bulk, na.rm = TRUE)
  minimum_tail_ess <- min(parameter_summary$ess_tail, na.rm = TRUE)
  maximum_mcse_sd_ratio <- max(
    parameter_summary$mcse_sd_ratio,
    na.rm = TRUE
  )
  geometry_pass <- divergent_transitions == 0 &&
    maximum_treedepth_transitions == 0 && minimum_ebfmi >= 0.2
  convergence_pass <- maximum_rhat <= 1.01 &&
    minimum_bulk_ess >= 400 && minimum_tail_ess >= 400 &&
    maximum_mcse_sd_ratio <= 0.10

  tibble::tibble(
    family = candidate$family,
    outcome = candidate$outcome,
    model_index = candidate$model_index,
    model_name = candidate$model_name,
    divergent_transitions = divergent_transitions,
    maximum_treedepth_transitions = maximum_treedepth_transitions,
    minimum_ebfmi = minimum_ebfmi,
    maximum_rhat = maximum_rhat,
    minimum_bulk_ess = minimum_bulk_ess,
    minimum_tail_ess = minimum_tail_ess,
    maximum_mcse_sd_ratio = maximum_mcse_sd_ratio,
    geometry_pass = geometry_pass,
    convergence_pass = convergence_pass,
    diagnostic_pass = geometry_pass && convergence_pass
  )
}

summarise_residual_correlation_pilot <- function(candidate) {
  if (is_sensitivity_pilot_skip(candidate)) return(candidate)
  if (!identical(candidate$family, "rescor")) {
    stop("Expected a residual-correlation pilot fit.", call. = FALSE)
  }
  variables <- c(
    "rho",
    "b_intercept",
    "b_rt",
    "b_protein",
    "b_moderator",
    "sd_study",
    "sd_arm",
    "sd_outcome"
  )
  draws <- posterior::as_draws_matrix(candidate$fit$draws(variables))
  purrr::map_dfr(colnames(draws), function(variable) {
    values <- draws[, variable]
    tibble::tibble(
      family = candidate$family,
      outcome = candidate$outcome,
      model_index = candidate$model_index,
      model_name = candidate$model_name,
      parameter = variable,
      estimate = mean(values),
      std_error = stats::sd(values),
      conf_low = unname(stats::quantile(values, 0.025)),
      median = stats::median(values),
      conf_high = unname(stats::quantile(values, 0.975)),
      probability_positive = mean(values > 0)
    )
  })
}

summarise_eiv_candidate_pilot <- function(candidate) {
  if (is_sensitivity_pilot_skip(candidate)) return(candidate)
  if (!identical(candidate$family, "eiv")) {
    stop("Expected an EIV pilot fit.", call. = FALSE)
  }
  variables <- c(
    "x_mean",
    "x_sd",
    "b_intercept",
    "b_rt",
    "b_protein",
    "b_protein_rt",
    "b_moderator",
    "sd_study",
    "sd_arm",
    "sd_outcome"
  )
  draws <- posterior::as_draws_matrix(candidate$fit$draws(variables))
  features <- sensitivity_model_features(candidate$model_index)
  active <- rep(TRUE, ncol(draws))
  if (!features$includes_protein) {
    active <- active & !stringr::str_detect(
      colnames(draws),
      "^b_protein($|\\[)"
    )
  }
  if (!features$includes_rt_interaction) {
    active <- active & !stringr::str_detect(
      colnames(draws),
      "^b_protein_rt($|\\[)"
    )
  }
  if (!features$includes_moderators) {
    active <- active & !stringr::str_detect(
      colnames(draws),
      "^b_moderator($|\\[)"
    )
  }
  draws <- draws[, active, drop = FALSE]
  purrr::map_dfr(colnames(draws), function(variable) {
    values <- draws[, variable]
    tibble::tibble(
      family = candidate$family,
      outcome = candidate$outcome,
      model_index = candidate$model_index,
      model_name = candidate$model_name,
      parameter = variable,
      estimate = mean(values),
      std_error = stats::sd(values),
      conf_low = unname(stats::quantile(values, 0.025)),
      median = stats::median(values),
      conf_high = unname(stats::quantile(values, 0.975)),
      probability_positive = mean(values > 0)
    )
  })
}

compute_cmdstan_pilot_weights <- function(fits) {
  if (is_sensitivity_pilot_skip(fits)) return(fits)
  loo_objects <- purrr::map(fits, function(candidate) {
    log_lik <- posterior::as_draws_matrix(
      candidate$fit$draws("log_lik_y")
    )
    loo::loo(log_lik)
  })
  names(loo_objects) <- purrr::map_chr(fits, "model_name")
  weights <- loo::loo_model_weights(loo_objects, method = "stacking")
  weights <- pmax(as.numeric(weights), 0)
  weights <- weights / sum(weights)
  names(weights) <- names(loo_objects)
  list(loo = loo_objects, weights = weights)
}

# The raw `log_lik_y` emitted by the EIV program conditions on the realised
# study, arm, and outcome effects.  That is useful as a diagnostic, but it is
# not the likelihood used for model comparison in the primary analysis.  The
# helpers below reconstruct the fixed predictor at each latent-intake draw and
# analytically integrate the Gaussian outcome random effects before PSIS-LOO.

extract_cmdstan_indexed_draw_matrix <- function(
    draws,
    prefix,
    expected_length
) {
  draws <- as.data.frame(draws)
  opening <- paste0(prefix, "[")
  variables <- names(draws)[startsWith(names(draws), opening)]
  indices <- suppressWarnings(as.integer(sub(
    "\\]$", "", substring(variables, nchar(opening) + 1L)
  )))
  valid <- !is.na(indices)
  variables <- variables[valid]
  indices <- indices[valid]
  expected_length <- as.integer(expected_length)
  if (
    expected_length < 1L || length(variables) != expected_length ||
      anyDuplicated(indices) ||
      !identical(sort(indices), seq_len(expected_length))
  ) {
    stop(
      "The CmdStan indexed draws are incomplete for `", prefix, "`.",
      call. = FALSE
    )
  }
  output <- as.matrix(draws[, variables[order(indices)], drop = FALSE])
  storage.mode(output) <- "double"
  if (anyNA(output) || any(!is.finite(output))) {
    stop("The CmdStan indexed draws are non-finite.", call. = FALSE)
  }
  output
}

extract_cmdstan_matrix_draw_array <- function(
    draws,
    prefix,
    nrow_matrix,
    ncol_matrix
) {
  draws <- as.data.frame(draws)
  opening <- paste0(prefix, "[")
  variables <- names(draws)[startsWith(names(draws), opening)]
  index_text <- sub("\\]$", "", substring(
    variables, nchar(opening) + 1L
  ))
  split_indices <- strsplit(index_text, ",", fixed = TRUE)
  row_index <- suppressWarnings(as.integer(vapply(
    split_indices, `[`, character(1), 1L
  )))
  column_index <- suppressWarnings(as.integer(vapply(
    split_indices, `[`, character(1), 2L
  )))
  nrow_matrix <- as.integer(nrow_matrix)
  ncol_matrix <- as.integer(ncol_matrix)
  key <- (column_index - 1L) * nrow_matrix + row_index
  if (
    length(variables) != nrow_matrix * ncol_matrix ||
      anyNA(row_index) || anyNA(column_index) || anyDuplicated(key) ||
      !identical(sort(key), seq_len(nrow_matrix * ncol_matrix))
  ) {
    stop(
      "The CmdStan matrix draws are incomplete for `", prefix, "`.",
      call. = FALSE
    )
  }
  ordered <- variables[order(key)]
  values <- as.matrix(draws[, ordered, drop = FALSE])
  storage.mode(values) <- "double"
  if (anyNA(values) || any(!is.finite(values))) {
    stop("The CmdStan matrix draws are non-finite.", call. = FALSE)
  }
  array(
    values,
    dim = c(nrow(values), nrow_matrix, ncol_matrix)
  )
}

evaluate_eiv_stan_basis <- function(x, stan_data) {
  x <- as.matrix(x)
  z <- x / as.numeric(stan_data$protein_scale)
  basis_type <- as.integer(stan_data$basis_type)
  k_basis <- as.integer(stan_data$K_basis)
  if (basis_type == 0L) {
    basis <- list(matrix(0, nrow(x), ncol(x)))
  } else if (basis_type == 1L) {
    basis <- list(z)
  } else if (basis_type == 2L) {
    scaled_knot <- as.numeric(stan_data$protein_knot) /
      as.numeric(stan_data$protein_scale)
    basis <- list(
      pmin(z, scaled_knot),
      pmax(z - scaled_knot, 0)
    )
  } else if (basis_type == 3L) {
    knots <- as.numeric(stan_data$rcs_knots)
    d <- lapply(knots[1:3], function(knot) {
      (pmax(z - knot, 0)^3 - pmax(z - knots[[4L]], 0)^3) /
        (knots[[4L]] - knot)
    })
    raw_one <- d[[1L]] - d[[3L]]
    raw_two <- d[[2L]] - d[[3L]]
    projection <- as.matrix(stan_data$rcs_projection)
    residual_one <- raw_one - projection[1L, 1L] -
      projection[2L, 1L] * z
    residual_two <- raw_two - projection[1L, 2L] -
      projection[2L, 2L] * z
    transform <- as.matrix(stan_data$rcs_basis_transform)
    basis <- list(
      z,
      residual_one * transform[1L, 1L] +
        residual_two * transform[2L, 1L],
      residual_one * transform[1L, 2L] +
        residual_two * transform[2L, 2L]
    )
  } else {
    stop("Unknown EIV protein basis type.", call. = FALSE)
  }
  if (
    length(basis) != k_basis ||
      any(vapply(basis, function(value) {
        !identical(dim(value), dim(x)) || anyNA(value) ||
          any(!is.finite(value))
      }, logical(1)))
  ) {
    stop("The reconstructed EIV basis is invalid.", call. = FALSE)
  }
  basis
}

make_eiv_integrated_loo_jobs <- function(
    fits,
    inputs,
    target_draws = 4000L
) {
  if (
    is_sensitivity_pilot_skip(fits) ||
      is_sensitivity_pilot_skip(inputs) ||
      all(purrr::map_lgl(inputs, is_sensitivity_pilot_skip))
  ) {
    return(new_sensitivity_pilot_skip("eiv"))
  }
  fit_indices <- unname(purrr::map_int(fits, "model_index"))
  input_indices <- unname(purrr::map_int(
    inputs, ~ as.integer(.x$specification$model_index)
  ))
  if (
    length(fits) != length(inputs) || anyDuplicated(fit_indices) ||
      anyDuplicated(input_indices) || !setequal(fit_indices, input_indices)
  ) {
    stop("The EIV fit/input ensemble is misaligned.", call. = FALSE)
  }
  fit_draws <- purrr::map_int(fits, function(candidate) {
    posterior::ndraws(candidate$fit$draws("b_intercept"))
  })
  loo_draws <- min(as.integer(target_draws), fit_draws)
  purrr::map(fits, function(candidate) {
    model_index <- as.integer(candidate$model_index)
    list(
      fitted_candidate = candidate,
      input = inputs[[match(model_index, input_indices)]],
      loo_draws = as.integer(loo_draws),
      fit_draws = as.integer(fit_draws[[match(model_index, fit_indices)]])
    )
  })
}

get_eiv_integrated_loo_components <- function(
    job,
    include_conditional_log_lik = TRUE
) {
  candidate <- job$fitted_candidate
  stan_data <- job$input$stan_data
  if (
    !candidate$family %in% c("eiv", "total_eiv") ||
      candidate$model_index != job$input$specification$model_index ||
      as.integer(stan_data$N) < 1L || as.integer(stan_data$N_arm) < 1L
  ) {
    stop("Invalid EIV integrated-LOO job.", call. = FALSE)
  }
  variables <- c(
    if (isTRUE(include_conditional_log_lik)) "log_lik_y" else character(),
    "x_true", "b_intercept", "b_rt", "b_protein",
    "b_protein_rt", "b_moderator", "sd_study", "L_study",
    "sd_arm", "sd_outcome"
  )
  all_draws <- posterior::as_draws_df(candidate$fit$draws(variables))
  draw_ids <- make_sensitivity_loo_draw_ids(
    nrow(all_draws), job$loo_draws
  )
  draws <- as.data.frame(all_draws)[draw_ids, , drop = FALSE]
  n_draws <- nrow(draws)
  n_observations <- as.integer(stan_data$N)
  k_basis <- as.integer(stan_data$K_basis)
  k_study <- as.integer(stan_data$K_study)
  arm <- as.integer(stan_data$arm)

  x_arm <- extract_cmdstan_indexed_draw_matrix(
    draws, "x_true", stan_data$N_arm
  )
  basis_arm <- evaluate_eiv_stan_basis(x_arm, stan_data)
  basis <- lapply(basis_arm, function(value) {
    value[, arm, drop = FALSE]
  })
  b_intercept <- as.numeric(draws$b_intercept)
  b_rt <- as.numeric(draws$b_rt)
  linear_predictor <- outer(b_intercept, rep(1, n_observations)) +
    outer(b_rt, as.numeric(stan_data$resistance_exercise_code))

  if (as.integer(stan_data$include_protein) == 1L) {
    b_protein <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_protein", k_basis
    )
    for (component in seq_len(k_basis)) {
      linear_predictor <- linear_predictor +
        basis[[component]] * b_protein[, component]
    }
  }
  if (as.integer(stan_data$include_rt_interaction) == 1L) {
    b_interaction <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_protein_rt", k_basis
    )
    rt <- matrix(
      as.numeric(stan_data$resistance_exercise_code),
      nrow = n_draws,
      ncol = n_observations,
      byrow = TRUE
    )
    for (component in seq_len(k_basis)) {
      linear_predictor <- linear_predictor +
        rt * basis[[component]] * b_interaction[, component]
    }
  }
  if (as.integer(stan_data$include_moderators) == 1L) {
    moderator_covariates <- list(
      as.numeric(stan_data$duration_z),
      as.numeric(stan_data$age_z)
    )
    if (identical(candidate$family, "total_eiv")) {
      if (
        !"body_mass_z" %in% names(stan_data) ||
          length(stan_data$body_mass_z) != n_observations ||
          any(!is.finite(stan_data$body_mass_z))
      ) {
        stop("The absolute total-intake EIV body-mass moderator is invalid.",
             call. = FALSE)
      }
      moderator_covariates <- c(
        moderator_covariates,
        list(as.numeric(stan_data$body_mass_z))
      )
    }
    b_moderator <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_moderator", length(moderator_covariates)
    )
    for (moderator in seq_along(moderator_covariates)) {
      linear_predictor <- linear_predictor + outer(
        b_moderator[, moderator],
        moderator_covariates[[moderator]]
      )
    }
  }

  sd_study <- extract_cmdstan_indexed_draw_matrix(
    draws, "sd_study", k_study
  )
  study_design <- c(
    list(matrix(1, n_draws, n_observations)),
    if (as.integer(stan_data$include_random_slope) == 1L) basis else list()
  )
  if (length(study_design) != k_study) {
    stop("The EIV study design is inconsistent.", call. = FALSE)
  }
  if (
    k_study > 1L &&
      as.integer(stan_data$correlated_study_effects) == 1L
  ) {
    l_study <- extract_cmdstan_matrix_draw_array(
      draws, "L_study", k_study, k_study
    )
  } else {
    l_study <- NULL
  }
  study_basis <- lapply(
    seq_len(k_study),
    function(component) matrix(NA_real_, n_draws, n_observations)
  )
  for (draw in seq_len(n_draws)) {
    transform <- if (is.null(l_study)) {
      diag(sd_study[draw, ], nrow = k_study)
    } else {
      diag(sd_study[draw, ], nrow = k_study) %*%
        matrix(l_study[draw, , ], k_study, k_study)
    }
    design <- do.call(cbind, lapply(
      study_design, function(value) value[draw, ]
    ))
    transformed <- design %*% transform
    for (component in seq_len(k_study)) {
      study_basis[[component]][draw, ] <- transformed[, component]
    }
  }

  manifest <- tibble::tibble(
    observation = seq_len(n_observations),
    yi = as.numeric(stan_data$y),
    vi = as.numeric(stan_data$se_y)^2,
    study_no = as.character(stan_data$study),
    arm_no = as.character(stan_data$arm),
    outcome_no = as.character(stan_data$outcome)
  )
  list(
    linear_predictor = linear_predictor,
    residual_sd = rep(0, n_draws),
    arm_intercept_sd = as.numeric(draws$sd_arm),
    outcome_intercept_sd = as.numeric(draws$sd_outcome),
    study_basis = study_basis,
    chain_id = as.integer(draws$.chain),
    conditional_log_lik = if (isTRUE(include_conditional_log_lik)) {
      extract_cmdstan_indexed_draw_matrix(
        draws, "log_lik_y", n_observations
      )
    } else {
      NULL
    },
    observation_manifest = manifest,
    structure = if (k_study == 1L) {
      "study_intercept"
    } else {
      paste0(
        if (as.integer(stan_data$correlated_study_effects) == 1L) {
          "correlated_"
        } else {
          "uncorrelated_"
        },
        "study_intercept_and_", k_study - 1L,
        "_latent_protein_basis_slopes"
      )
    }
  )
}

compute_eiv_integrated_psis_loo <- function(
    job,
    k_threshold = 0.7,
    cores = 1L
) {
  if (is_sensitivity_pilot_skip(job)) return(job)
  candidate <- job$fitted_candidate
  components <- get_eiv_integrated_loo_components(job)
  integrated_log_lik <- integrate_primary_study_basis_log_lik(
    linear_predictor = components$linear_predictor,
    observed = components$observation_manifest$yi,
    known_variance = components$observation_manifest$vi,
    study_id = components$observation_manifest$study_no,
    arm_id = components$observation_manifest$arm_no,
    outcome_id = components$observation_manifest$outcome_no,
    residual_sd = components$residual_sd,
    arm_intercept_sd = components$arm_intercept_sd,
    outcome_intercept_sd = components$outcome_intercept_sd,
    study_basis = components$study_basis
  )
  make_loo <- function(log_lik) {
    relative_efficiency <- loo::relative_eff(
      exp(log_lik),
      chain_id = components$chain_id,
      cores = as.integer(cores)
    )
    loo::loo(
      log_lik,
      r_eff = relative_efficiency,
      save_psis = FALSE,
      cores = as.integer(cores)
    )
  }
  conditional_loo <- make_loo(components$conditional_log_lik)
  integrated_loo <- make_loo(integrated_log_lik)
  comparison <- dplyr::bind_rows(
    summarise_primary_loo_method(
      conditional_loo,
      "conditional_observation_level_psis",
      k_threshold
    ),
    summarise_primary_loo_method(
      integrated_loo,
      "outcome_effects_integrated_latent_exposure_psis",
      k_threshold
    )
  ) |>
    dplyr::mutate(
      outcome = candidate$outcome,
      model_index = candidate$model_index,
      model_name = candidate$model_name,
      study_structure = components$structure,
      loo_draws = as.integer(job$loo_draws),
      .before = 1L
    )
  list(
    family = candidate$family,
    outcome = candidate$outcome,
    model_index = candidate$model_index,
    model_name = candidate$model_name,
    study_structure = components$structure,
    observation_manifest = components$observation_manifest,
    loo = integrated_loo,
    conditional_loo = conditional_loo,
    remedy = "outcome_effects_integrated_latent_exposure_psis",
    k_threshold = k_threshold,
    loo_draws = as.integer(job$loo_draws),
    comparison = comparison
  )
}

summarise_eiv_integrated_loo <- function(results) {
  if (is_sensitivity_pilot_skip(results)) return(results)
  if (
    length(results) < 1L ||
      !all(purrr::map_lgl(results, ~ identical(.x$family, "eiv"))) ||
      !all(purrr::map_lgl(results, ~ is.data.frame(.x$comparison)))
  ) {
    stop("The integrated EIV LOO results are invalid.", call. = FALSE)
  }
  purrr::map_dfr(results, "comparison") |>
    dplyr::arrange(model_index, method)
}

select_eiv_exact_loo_jobs <- function(
    integrated_results,
    fitted_candidates,
    inputs,
    k_threshold = 0.7
) {
  if (
    is_sensitivity_pilot_skip(integrated_results) ||
      is_sensitivity_pilot_skip(fitted_candidates) ||
      is_sensitivity_pilot_skip(inputs)
  ) {
    return(new_sensitivity_pilot_skip("eiv"))
  }
  fitted_indices <- unname(purrr::map_int(
    fitted_candidates, "model_index"
  ))
  input_indices <- unname(purrr::map_int(
    inputs, ~ as.integer(.x$specification$model_index)
  ))
  result_indices <- unname(purrr::map_int(
    integrated_results, "model_index"
  ))
  if (
    anyDuplicated(fitted_indices) || anyDuplicated(input_indices) ||
      anyDuplicated(result_indices) ||
      !setequal(fitted_indices, input_indices) ||
      !setequal(fitted_indices, result_indices)
  ) {
    stop("The EIV exact-LOO inputs are misaligned.", call. = FALSE)
  }
  cases <- purrr::map_dfr(integrated_results, function(result) {
    pareto_k <- loo::pareto_k_values(result$loo)
    observations <- which(pareto_k > k_threshold)
    tibble::tibble(
      model_index = as.integer(result$model_index),
      observation = as.integer(observations),
      pareto_k = as.numeric(pareto_k[observations])
    )
  }) |>
    dplyr::arrange(model_index, observation) |>
    dplyr::mutate(case_id = dplyr::row_number(), .before = 1L)
  if (nrow(cases) < 1L) return(list())
  purrr::pmap(cases, function(
      case_id,
      model_index,
      observation,
      pareto_k
  ) {
    fit_position <- match(model_index, fitted_indices)
    input_position <- match(model_index, input_indices)
    result_position <- match(model_index, result_indices)
    integrated_result <- integrated_results[[result_position]]
    if (
      anyNA(c(fit_position, input_position, result_position)) ||
        observation < 1L ||
        observation > nrow(integrated_result$observation_manifest)
    ) {
      stop("An EIV exact-LOO case failed validation.", call. = FALSE)
    }
    list(
      case_id = as.integer(case_id),
      family = "eiv",
      outcome = fitted_candidates[[fit_position]]$outcome,
      model_index = as.integer(model_index),
      model_name = fitted_candidates[[fit_position]]$model_name,
      observation = as.integer(observation),
      pareto_k = as.numeric(pareto_k),
      fitted_candidate = fitted_candidates[[fit_position]],
      input = inputs[[input_position]],
      integrated_result = integrated_result
    )
  })
}

summarise_eiv_exact_loo_jobs <- function(jobs) {
  if (is_sensitivity_pilot_skip(jobs)) return(jobs)
  if (length(jobs) < 1L) {
    return(tibble::tibble(
      case_id = integer(), model_index = integer(), model_name = character(),
      observation = integer(), pareto_k = double()
    ))
  }
  purrr::map_dfr(jobs, function(job) {
    tibble::tibble(
      case_id = job$case_id,
      model_index = job$model_index,
      model_name = job$model_name,
      observation = job$observation,
      pareto_k = job$pareto_k
    )
  })
}

make_eiv_leave_one_out_stan_data <- function(stan_data, observation) {
  observation <- as.integer(observation)
  n_observations <- as.integer(stan_data$N)
  if (
    length(observation) != 1L || is.na(observation) ||
      observation < 1L || observation > n_observations
  ) {
    stop("Invalid EIV held-out observation.", call. = FALSE)
  }
  row_fields <- c(
    "y", "se_y", "arm", "study", "outcome",
    "resistance_exercise_code", "age_z", "duration_z"
  )
  if ("body_mass_z" %in% names(stan_data)) {
    row_fields <- c(row_fields, "body_mass_z")
  }
  if (
    !all(row_fields %in% names(stan_data)) ||
      any(vapply(stan_data[row_fields], length, integer(1)) != n_observations)
  ) {
    stop("The EIV row-level Stan data are incomplete.", call. = FALSE)
  }
  stan_data$N <- n_observations - 1L
  for (field in row_fields) {
    stan_data[[field]] <- stan_data[[field]][-observation]
  }
  stan_data
}

fit_eiv_exact_loo_case <- function(job, model, enabled, config) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("eiv"))
  if (is_sensitivity_pilot_skip(model)) {
    stop("The EIV Stan model was not compiled.", call. = FALSE)
  }
  training_data <- make_eiv_leave_one_out_stan_data(
    job$input$stan_data,
    job$observation
  )
  specification <- job$fitted_candidate$specification
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  output_dir <- make_sensitivity_cmdstan_output_dir(
    paste0("eiv_exact_case_", job$case_id),
    job$outcome,
    job$model_index
  )
  start_time <- proc.time()[["elapsed"]]
  fit <- model$sample(
    data = training_data,
    chains = config$chains,
    parallel_chains = config$parallel_chains,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    adapt_delta = specification$adapt_delta,
    max_treedepth = specification$max_treedepth,
    seed = config$seed + job$case_id,
    refresh = 100,
    output_dir = output_dir
  )
  candidate <- list(
    family = "eiv",
    outcome = job$outcome,
    model_index = job$model_index,
    model_name = job$model_name,
    specification = specification,
    grid = job$input$grid,
    output_dir = output_dir,
    fit = fit
  )
  list(
    case_id = job$case_id,
    family = "eiv",
    outcome = job$outcome,
    model_index = job$model_index,
    model_name = job$model_name,
    observation = job$observation,
    pareto_k = job$pareto_k,
    fitted_candidate = candidate,
    input = job$input,
    integrated_result = job$integrated_result,
    refit_seconds = proc.time()[["elapsed"]] - start_time
  )
}

summarise_eiv_exact_loo_refit <- function(refit) {
  if (is_sensitivity_pilot_skip(refit)) return(refit)
  diagnose_cmdstan_pilot_fit(refit$fitted_candidate) |>
    dplyr::mutate(
      case_id = refit$case_id,
      observation = refit$observation,
      old_psis_pareto_k = refit$pareto_k,
      refit_seconds = refit$refit_seconds,
      .before = 1L
    )
}

compute_eiv_exact_loo_case <- function(refit, loo_draws = 4000L) {
  if (is_sensitivity_pilot_skip(refit)) return(refit)
  prediction_job <- list(
    fitted_candidate = refit$fitted_candidate,
    input = refit$input,
    loo_draws = min(
      as.integer(loo_draws),
      posterior::ndraws(
        refit$fitted_candidate$fit$draws("b_intercept")
      )
    )
  )
  components <- get_eiv_integrated_loo_components(
    prediction_job,
    include_conditional_log_lik = FALSE
  )
  integrated_log_lik <- integrate_primary_study_basis_log_lik(
    linear_predictor = components$linear_predictor,
    observed = components$observation_manifest$yi,
    known_variance = components$observation_manifest$vi,
    study_id = components$observation_manifest$study_no,
    arm_id = components$observation_manifest$arm_no,
    outcome_id = components$observation_manifest$outcome_no,
    residual_sd = components$residual_sd,
    arm_intercept_sd = components$arm_intercept_sd,
    outcome_intercept_sd = components$outcome_intercept_sd,
    study_basis = components$study_basis
  )
  observation <- as.integer(refit$observation)
  held_out_log_lik <- integrated_log_lik[, observation]
  maximum_log_lik <- max(held_out_log_lik)
  scaled_likelihood <- exp(held_out_log_lik - maximum_log_lik)
  exact_elpd <- maximum_log_lik + log(mean(scaled_likelihood))
  exact_mcse <- tryCatch(
    posterior::mcse_mean(scaled_likelihood) / mean(scaled_likelihood),
    error = function(e) NA_real_
  )
  patched_loo <- refit$integrated_result$loo
  old_elpd <- patched_loo$pointwise[observation, "elpd_loo"]
  old_p_loo <- patched_loo$pointwise[observation, "p_loo"]
  full_data_lpd <- old_elpd + old_p_loo
  patched_loo$pointwise[observation, "elpd_loo"] <- exact_elpd
  patched_loo$pointwise[observation, "p_loo"] <-
    full_data_lpd - exact_elpd
  patched_loo$pointwise[observation, "looic"] <- -2 * exact_elpd
  if ("mcse_elpd_loo" %in% colnames(patched_loo$pointwise)) {
    patched_loo$pointwise[observation, "mcse_elpd_loo"] <- exact_mcse
  }
  patched_loo$diagnostics$pareto_k[[observation]] <- NA_real_
  if (!is.null(patched_loo$diagnostics$n_eff)) {
    patched_loo$diagnostics$n_eff[[observation]] <- NA_real_
  }
  patched_loo <- update_primary_exact_loo_estimates(patched_loo)
  exact_refit_summary <- summarise_eiv_exact_loo_refit(refit) |>
    dplyr::mutate(
      old_psis_elpd_loo = old_elpd,
      exact_elpd_loo = exact_elpd,
      exact_minus_psis_elpd = exact_elpd - old_elpd,
      exact_elpd_mcse = exact_mcse
    )
  list(
    case_id = refit$case_id,
    family = "eiv",
    outcome = refit$outcome,
    model_index = refit$model_index,
    model_name = refit$model_name,
    study_structure = refit$integrated_result$study_structure,
    observation_manifest = components$observation_manifest,
    loo = patched_loo,
    conditional_loo = refit$integrated_result$conditional_loo,
    remedy = "outcome_effects_integrated_psis_plus_exact",
    k_threshold = refit$integrated_result$k_threshold,
    exact_refit_summary = exact_refit_summary
  )
}

build_cmdstan_pilot_model_average <- function(
    fits,
    weight_result,
    n_draws = 1000L,
    seed = 1988L
) {
  if (is_sensitivity_pilot_skip(fits)) return(fits)
  weights <- weight_result$weights
  counts <- allocate_primary_mixture_models(weights, n_draws)
  set.seed(seed)
  allocation <- sample(rep(seq_along(fits), times = counts))
  grid <- fits[[1L]]$grid
  if (!all(purrr::map_lgl(fits, ~ identical(.x$grid, grid)))) {
    stop("CmdStan prediction grids do not match.", call. = FALSE)
  }
  draws <- matrix(NA_real_, n_draws, nrow(grid))
  for (index in seq_along(fits)) {
    rows <- which(allocation == index)
    if (!length(rows)) next
    model_draws <- posterior::as_draws_matrix(
      fits[[index]]$fit$draws("mu_grid")
    )
    selected <- sample(seq_len(nrow(model_draws)), length(rows),
                       replace = nrow(model_draws) < length(rows))
    draws[rows, ] <- model_draws[selected, , drop = FALSE]
  }
  output <- list(
    draws = draws,
    grid = dplyr::transmute(
      grid,
      resistance_exercise_code,
      resistance_exercise = ifelse(
        resistance_exercise_code == 0, "YES", "NO"
      ),
      protein_intake = protein_intake_raw
    ),
    model_weights = weights,
    settings = list(n_draws = n_draws, seed = seed,
                    prediction_type = "eiv")
  )
  class(output) <- c("primary_model_average", class(output))
  output
}

make_residual_correlation_stan_data <- function(
    data,
    specification,
    priors,
    cluster = c("study", "arm")
) {
  cluster <- match.arg(cluster)
  features <- sensitivity_model_features(specification$model_index)
  if (features$protein_family != "linear" ||
      features$includes_random_slope) {
    stop("The residual-correlation pilot uses a fixed-slope linear candidate.",
         call. = FALSE)
  }
  data <- add_primary_standardised_moderators(data)
  cluster_id <- if (cluster == "study") {
    data$study_no
  } else {
    interaction(data$study_no, data$arm_no, drop = TRUE)
  }
  data <- data |>
    dplyr::mutate(
      residual_cluster = as.character(cluster_id),
      study_index = dplyr::dense_rank(study_no),
      arm_index = dplyr::dense_rank(interaction(study_no, arm_no,
                                                drop = TRUE)),
      outcome_index = dplyr::dense_rank(interaction(
        study_no, arm_no, outcome_no, drop = TRUE
      ))
    ) |>
    dplyr::arrange(residual_cluster)
  runs <- rle(data$residual_cluster)
  cluster_end <- cumsum(runs$lengths)
  cluster_start <- c(1L, utils::head(cluster_end, -1L) + 1L)
  max_cluster_size <- max(runs$lengths)
  if (max_cluster_size <= 1L) {
    stop(
      "Residual correlation requires a cluster with repeated observations.",
      call. = FALSE
    )
  }
  rho_lower <- -1 / (max_cluster_size - 1) + 1e-8
  rho_raw_prior_location <- stats::qlogis(
    -rho_lower / (1 - rho_lower)
  )
  stan_priors <- make_stan_priors(priors)

  list(
    family = "rescor",
    specification = specification,
    cluster = cluster,
    stan_data = list(
      N = nrow(data),
      y = data$yi,
      se_y = sqrt(data$vi),
      protein_intake = data$protein_intake,
      resistance_exercise_code = data$resistance_exercise_code,
      age_z = data$age_outcome_z,
      duration_z = data$duration_outcome_z,
      include_moderators = as.integer(features$includes_moderators),
      include_rt_interaction =
        as.integer(features$includes_rt_interaction),
      N_study = dplyr::n_distinct(data$study_index),
      study = data$study_index,
      N_arm = dplyr::n_distinct(data$arm_index),
      arm = data$arm_index,
      N_outcome = dplyr::n_distinct(data$outcome_index),
      outcome = data$outcome_index,
      N_cluster = length(runs$lengths),
      cluster_start = cluster_start,
      cluster_end = cluster_end,
      rho_lower = rho_lower,
      rho_raw_prior_location = rho_raw_prior_location,
      prior_intercept_location = stan_priors$prior_intercept_location,
      prior_intercept_scale = stan_priors$prior_intercept_scale,
      prior_rt_location = stan_priors$prior_rt_location,
      prior_rt_scale = stan_priors$prior_rt_scale
    )
  )
}

fit_residual_correlation_cmdstan_pilot <- function(
    input,
    model,
    enabled,
    config
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("rescor"))
  specification <- input$specification
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  output_dir <- make_sensitivity_cmdstan_output_dir(
    "rescor",
    input$specification$outcome,
    input$specification$model_index
  )
  fit <- model$sample(
    data = input$stan_data,
    chains = config$chains,
    parallel_chains = config$parallel_chains,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    adapt_delta = specification$adapt_delta,
    max_treedepth = specification$max_treedepth,
    seed = config$seed,
    refresh = 100,
    output_dir = output_dir
  )
  list(
    family = "rescor",
    outcome = input$specification$outcome,
    model_index = input$specification$model_index,
    model_name = paste0(
      input$specification$model_name,
      " [residual correlation by ", input$cluster, "]"
    ),
    specification = specification,
    output_dir = output_dir,
    fit = fit
  )
}

audit_sensitivity_pilot_structure <- function(
    final_specifications,
    protein_pilot_specifications,
    comparison_pilot_specifications,
    mundlak_inputs,
    mundlak_remedy_inputs,
    mundlak_remedy_config,
    within_inputs,
    exposure_data,
    eiv_inputs,
    residual_input
) {
  final_table <- dplyr::bind_rows(final_specifications)
  protein_pilot_indices <- unname(purrr::map_int(
    protein_pilot_specifications, "model_index"
  ))
  comparison_pilot_indices <- unname(purrr::map_int(
    comparison_pilot_specifications, "model_index"
  ))
  input_indices <- function(x) unname(purrr::map_int(
    x, ~ .x$specification$model_index
  ))
  remedy_indices <- input_indices(mundlak_remedy_inputs)
  regularised_remedy_indices <- unname(purrr::map_int(
    purrr::keep(
      mundlak_remedy_inputs,
      ~ isTRUE(.x$specification$regularised_variance_priors)
    ),
    ~ .x$specification$model_index
  ))
  residual_data <- residual_input$stan_data
  residual_prior_median <- residual_data$rho_lower +
    (1 - residual_data$rho_lower) *
      stats::plogis(residual_data$rho_raw_prior_location)
  eiv_has_extra_residual_scale <- any(purrr::map_lgl(
    eiv_inputs,
    ~ "prior_residual_sd_scale" %in% names(.x$stan_data)
  ))
  if (
    nrow(final_table) != 104L ||
      !all(protein_pilot_indices %in% comparison_pilot_indices) ||
      !identical(input_indices(mundlak_inputs), comparison_pilot_indices) ||
      !identical(
        remedy_indices,
        mundlak_remedy_config$remedy_model_indices
      ) ||
      !all(
        mundlak_remedy_config$regularised_model_indices %in%
          regularised_remedy_indices
      ) ||
      !identical(input_indices(within_inputs), comparison_pilot_indices) ||
      !identical(input_indices(eiv_inputs), comparison_pilot_indices) ||
      exposure_data$audit$rows < 1L ||
      is.null(residual_data) ||
      abs(residual_prior_median) > sqrt(.Machine$double.eps) ||
      "prior_residual_sd_scale" %in% names(residual_data) ||
      eiv_has_extra_residual_scale
  ) {
    stop("The sensitivity pilot structure failed its dry-run audit.",
         call. = FALSE)
  }
  tibble::tibble(
    check = c(
      "final candidate registry",
      "Mundlak pilot inputs",
      "Mundlak focused-remedy inputs",
      "strict within-study contrast inputs",
      "achieved/assignment common sample",
      "Stan EIV pilot inputs",
      "residual-correlation Stan input"
    ),
    status = "PASS",
    detail = c(
      "104 final outcome-model specifications",
      paste(
        length(mundlak_inputs),
        "candidates (including two unchanged null models)"
      ),
      paste0(
        length(mundlak_remedy_inputs),
        " diagnostic refits; model ",
        paste(mundlak_remedy_config$regularised_model_indices,
              collapse = ", "),
        " gains regularised variance priors"
      ),
      paste(length(within_inputs), "contrast/null candidates"),
      paste(exposure_data$audit$rows, "common-sample rows"),
      paste(
        length(eiv_inputs),
        "Stan candidate data lists (including two null models)"
      ),
      paste0(
        "model ", residual_input$specification$model_index,
        "; correlation prior centred at zero"
      )
    )
  )
}

# Full Stan errors-in-variables expansion -----------------------------------

make_eiv_full_expansion_config <- function(
    base_config,
    iter_warmup = 1500L,
    iter_sampling = 2000L,
    seed_offset = 20000L
) {
  config <- base_config
  config$outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )
  config$model_indices <- seq_len(26L)
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains < 1L ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("Invalid full EIV expansion configuration.", call. = FALSE)
  }
  config
}

make_eiv_full_remedy_config <- function(
    full_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    adapt_delta = 0.999,
    max_treedepth = 15L,
    seed_offset = 1000L
) {
  config <- full_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$adapt_delta_floor <- as.numeric(adapt_delta)
  config$max_treedepth_floor <- as.integer(max_treedepth)
  config$seed <- as.integer(full_config$seed + seed_offset)
  if (
    config$iter_warmup < full_config$iter_warmup ||
      config$iter_sampling <= full_config$iter_sampling ||
      !is.finite(config$adapt_delta_floor) ||
      config$adapt_delta_floor < full_config$adapt_delta_floor ||
      config$adapt_delta_floor >= 1 ||
      config$max_treedepth_floor < full_config$max_treedepth_floor
  ) {
    stop("Invalid full EIV remedy configuration.", call. = FALSE)
  }
  config
}

read_eiv_pilot_completed_checkpoint <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("The completed EIV pilot checkpoint is unavailable.", call. = FALSE)
  }
  checkpoint <- qs::qread(path)
  manifest <- checkpoint$accepted_fit_manifest
  required <- c(
    "family", "outcome", "model_index", "model_name", "fit_path",
    "stage", "bytes", "sha256"
  )
  valid_hash <- function(index) {
    identical(
      digest::digest(manifest$fit_path[[index]], file = TRUE, algo = "sha256"),
      manifest$sha256[[index]]
    )
  }
  if (
    !is.list(checkpoint) || !identical(checkpoint$family, "eiv") ||
      !is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      nrow(manifest) != 7L || any(manifest$family != "eiv") ||
      any(manifest$outcome != "whole_body") ||
      !identical(
        unname(sort(as.integer(manifest$model_index))),
        c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
      ) ||
      any(!file.exists(manifest$fit_path)) ||
      !all(vapply(seq_len(nrow(manifest)), valid_hash, logical(1))) ||
      !is.data.frame(checkpoint$pilot_diagnostics) ||
      !all(checkpoint$pilot_diagnostics$diagnostic_pass)
  ) {
    stop("The completed EIV pilot checkpoint is invalid.", call. = FALSE)
  }
  checkpoint
}

prepare_eiv_full_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome
) {
  outcome <- specification$outcome
  data <- data_by_outcome[[outcome]]
  priors <- priors_by_outcome[[outcome]]
  if (is.null(data) || is.null(priors)) {
    stop("Outcome-specific EIV data or priors are unavailable.", call. = FALSE)
  }
  make_eiv_candidate_data(data, specification, priors)
}

select_pending_eiv_full_inputs <- function(inputs, reused_manifest) {
  reused_keys <- paste(reused_manifest$outcome, reused_manifest$model_index)
  input_keys <- purrr::map_chr(
    inputs,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  pending <- inputs[!input_keys %in% reused_keys]
  pending_keys <- purrr::map_chr(
    pending,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  if (
    length(inputs) != 104L || length(pending) != 97L ||
      anyDuplicated(input_keys) || anyDuplicated(pending_keys)
  ) {
    stop("The pending full EIV input set is invalid.", call. = FALSE)
  }
  pending
}

read_eiv_fit_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("The EIV fit file is unavailable.", call. = FALSE)
  }
  candidate <- qs::qread(path)
  output_files <- tryCatch(
    candidate$fit$output_files(),
    error = function(e) character()
  )
  if (
    !is.list(candidate) || !identical(candidate$family, "eiv") ||
      !inherits(candidate$fit, "CmdStanMCMC") ||
      is.null(candidate$outcome) || is.null(candidate$model_index) ||
      length(output_files) < 1L || any(!file.exists(output_files))
  ) {
    stop("The EIV fit file is invalid or its CmdStan CSV files are missing.",
         call. = FALSE)
  }
  candidate
}

write_eiv_fit_file <- function(candidate, stage) {
  if (
    !is.list(candidate) || !identical(candidate$family, "eiv") ||
      !inherits(candidate$fit, "CmdStanMCMC")
  ) {
    stop("Cannot archive an invalid EIV fit.", call. = FALSE)
  }
  path <- sensitivity_fit_file_path(
    "eiv", candidate$outcome, candidate$model_index, stage
  )
  write_atomic_qs(candidate, path)
}

make_eiv_full_cmdstan_output_dir <- function(stage, outcome, model_index) {
  safe <- function(x) gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  path <- here::here(
    "_targets", "user", "sensitivity_cmdstan", "eiv_full",
    safe(stage), safe(outcome),
    paste0("model_", stringr::str_pad(model_index, 2L, pad = "0"))
  )
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

fit_eiv_full_candidate_to_file <- function(
    input,
    model,
    enabled,
    config,
    stage = "initial"
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("eiv"))
  if (is_sensitivity_pilot_skip(model)) {
    stop("The EIV Stan model was not compiled.", call. = FALSE)
  }
  specification <- input$specification
  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the full EIV fit.", call. = FALSE)
  }
  # A machine restart can occur after the atomic fit archive is installed but
  # before {targets} commits the branch metadata. Reuse a complete matching
  # archive in that case instead of sampling the model a second time.
  existing_path <- sensitivity_fit_file_path(
    "eiv", specification$outcome, specification$model_index, stage
  )
  existing_candidate <- if (file.exists(existing_path)) {
    tryCatch(read_eiv_fit_file(existing_path), error = function(e) NULL)
  } else {
    NULL
  }
  if (
    !is.null(existing_candidate) &&
      identical(existing_candidate$outcome, specification$outcome) &&
      identical(
        as.integer(existing_candidate$model_index),
        as.integer(specification$model_index)
      )
  ) {
    return(normalizePath(existing_path, winslash = "/", mustWork = TRUE))
  }
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  local_seed <- as.integer(
    config$seed + 100L * (outcome_offset - 1L) +
      specification$model_index
  )
  output_dir <- make_eiv_full_cmdstan_output_dir(
    stage, specification$outcome, specification$model_index
  )
  fit <- model$sample(
    data = input$stan_data,
    chains = config$chains,
    parallel_chains = config$parallel_chains,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    adapt_delta = specification$adapt_delta,
    max_treedepth = specification$max_treedepth,
    seed = local_seed,
    refresh = 100,
    output_dir = output_dir
  )
  candidate <- list(
    family = "eiv",
    outcome = specification$outcome,
    model_index = as.integer(specification$model_index),
    model_name = paste0(specification$model_name, " [EIV Stan]"),
    specification = specification,
    grid = input$grid,
    output_dir = output_dir,
    fit = fit
  )
  write_eiv_fit_file(candidate, stage = stage)
}

make_eiv_fit_file_manifest <- function(paths, stage) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) == 0L) {
    return(tibble::tibble(
      family = character(), outcome = character(), model_index = integer(),
      model_name = character(), fit_path = character(), stage = character(),
      bytes = double(), sha256 = character()
    ))
  }
  if (anyNA(paths) || any(!file.exists(paths))) {
    stop("The EIV fit-file set is incomplete.", call. = FALSE)
  }
  manifest <- purrr::map_dfr(paths, function(path) {
    candidate <- read_eiv_fit_file(path)
    tibble::tibble(
      family = "eiv",
      outcome = candidate$outcome,
      model_index = as.integer(candidate$model_index),
      model_name = candidate$model_name,
      fit_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      stage = stage,
      bytes = file.info(path)$size,
      sha256 = digest::digest(path, file = TRUE, algo = "sha256")
    )
  })
  if (anyDuplicated(paste(manifest$outcome, manifest$model_index))) {
    stop("The EIV fit-file manifest contains duplicate models.", call. = FALSE)
  }
  manifest
}

audit_eiv_full_expansion <- function(
    inputs,
    pending_inputs,
    pilot_checkpoint,
    config
) {
  input_keys <- purrr::map_chr(
    inputs,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  pending_keys <- purrr::map_chr(
    pending_inputs,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  expected_keys <- as.vector(outer(config$outcomes, seq_len(26L), paste))
  reused <- pilot_checkpoint$accepted_fit_manifest
  reused_keys <- paste(reused$outcome, reused$model_index)
  pilot_matches <- vapply(seq_len(nrow(reused)), function(index) {
    candidate <- read_eiv_fit_file(reused$fit_path[[index]])
    input <- inputs[[match(reused_keys[[index]], input_keys)]]
    identical(candidate$grid, input$grid) &&
      identical(candidate$model_index, input$specification$model_index)
  }, logical(1))
  valid_inputs <- purrr::map_lgl(inputs, function(input) {
    identical(input$family, "eiv") &&
      is.list(input$stan_data) && input$stan_data$N > 0L &&
      input$stan_data$N_arm > 0L &&
      all(is.finite(input$stan_data$x_se)) &&
      all(input$stan_data$x_se > 0)
  })
  if (
    length(inputs) != 104L || length(pending_inputs) != 97L ||
      nrow(reused) != 7L || anyDuplicated(input_keys) ||
      anyDuplicated(pending_keys) || !setequal(input_keys, expected_keys) ||
      length(intersect(pending_keys, reused_keys)) > 0L ||
      !setequal(c(pending_keys, reused_keys), expected_keys) ||
      !all(valid_inputs) || !all(pilot_matches) ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("The full EIV expansion failed its structural audit.", call. = FALSE)
  }
  tibble::tibble(
    check = c(
      "complete candidate registry",
      "reused passing whole-body pilots",
      "new independently checkpointed fits",
      "arm-level measurement-error inputs",
      "pilot-informed initial schedule"
    ),
    status = "PASS",
    detail = c(
      "104 candidates: 26 models x 4 outcomes",
      "7 accepted fits: models 1, 2, 3, 4, 5, 13, and 26",
      "97 pending fits: 19 whole-body plus 26 for each other outcome",
      "positive protein-intake standard errors in every Stan input",
      "1,500 warmup and 2,000 retained iterations per chain"
    )
  )
}

combine_eiv_initial_fit_manifest <- function(pilot_manifest, initial_manifest) {
  manifest <- dplyr::bind_rows(pilot_manifest, initial_manifest) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  if (
    nrow(manifest) != 104L || any(manifest$family != "eiv") ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$fit_path))
  ) {
    stop("The initial full EIV fit manifest is incomplete.", call. = FALSE)
  }
  manifest
}

make_eiv_full_fit_records <- function(fit_manifest, inputs) {
  input_keys <- purrr::map_chr(
    inputs,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  manifest_keys <- paste(fit_manifest$outcome, fit_manifest$model_index)
  if (
    nrow(fit_manifest) != 104L || length(inputs) != 104L ||
      anyDuplicated(input_keys) || anyDuplicated(manifest_keys) ||
      !setequal(input_keys, manifest_keys)
  ) {
    stop("The full EIV fit/input registry is misaligned.", call. = FALSE)
  }
  purrr::map(seq_len(nrow(fit_manifest)), function(row) {
    record <- purrr::map(as.list(fit_manifest[row, , drop = FALSE]), unname)
    record$input <- inputs[[match(manifest_keys[[row]], input_keys)]]
    record
  })
}

diagnose_eiv_fit_file <- function(path) {
  diagnose_cmdstan_pilot_fit(read_eiv_fit_file(path))
}

compile_eiv_full_diagnostic_audit <- function(
    diagnostics,
    expected_models = 104L
) {
  summary <- dplyr::bind_rows(diagnostics) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    ) |>
    dplyr::mutate(
      status = dplyr::if_else(
        diagnostic_pass, "PASS", "REQUIRES REMEDIATION"
      )
    )
  if (
    nrow(summary) != expected_models || any(summary$family != "eiv") ||
      anyDuplicated(paste(summary$outcome, summary$model_index))
  ) {
    stop("The full EIV diagnostic audit is incomplete.", call. = FALSE)
  }
  list(
    family = "eiv",
    diagnostic_summary = summary,
    parameter_failures = tibble::tibble(),
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

select_eiv_full_remedy_inputs <- function(inputs, initial_audit) {
  failed_keys <- initial_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::transmute(key = paste(outcome, model_index)) |>
    dplyr::pull(key)
  input_keys <- purrr::map_chr(
    inputs,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  selected <- inputs[input_keys %in% failed_keys]
  selected_keys <- purrr::map_chr(
    selected,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  if (
    !identical(initial_audit$family, "eiv") ||
      anyDuplicated(failed_keys) || length(selected) != length(failed_keys) ||
      !setequal(selected_keys, failed_keys)
  ) {
    stop("The full EIV remedy input set is invalid.", call. = FALSE)
  }
  selected
}

combine_eiv_full_remedy_diagnostic_audit <- function(
    initial_audit,
    remedy_diagnostics
) {
  failed <- dplyr::filter(initial_audit$diagnostic_summary, !diagnostic_pass)
  if (nrow(failed) == 0L) return(initial_audit)
  remedy <- dplyr::bind_rows(remedy_diagnostics) |>
    dplyr::mutate(
      status = dplyr::if_else(
        diagnostic_pass, "PASS", "REQUIRES REMEDIATION"
      )
    )
  failed_keys <- paste(failed$outcome, failed$model_index)
  remedy_keys <- paste(remedy$outcome, remedy$model_index)
  if (
    nrow(remedy) != nrow(failed) || anyDuplicated(remedy_keys) ||
      !setequal(failed_keys, remedy_keys)
  ) {
    stop("The full EIV remedy diagnostics are incomplete.", call. = FALSE)
  }
  summary <- initial_audit$diagnostic_summary |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  list(
    family = "eiv",
    diagnostic_summary = summary,
    parameter_failures = tibble::tibble(),
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

combine_eiv_full_accepted_fit_manifest <- function(
    initial_manifest,
    remedy_manifest,
    final_audit
) {
  remedy_keys <- if (nrow(remedy_manifest)) {
    paste(remedy_manifest$outcome, remedy_manifest$model_index)
  } else {
    character()
  }
  accepted <- initial_manifest |>
    dplyr::filter(!paste(outcome, model_index) %in% remedy_keys) |>
    dplyr::bind_rows(remedy_manifest) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  if (
    !isTRUE(final_audit$all_sampling_diagnostics_pass) ||
      nrow(accepted) != 104L || any(accepted$family != "eiv") ||
      anyDuplicated(paste(accepted$outcome, accepted$model_index)) ||
      any(!file.exists(accepted$fit_path))
  ) {
    stop("The accepted full EIV fit manifest is invalid.", call. = FALSE)
  }
  accepted
}

combine_eiv_remedy_fit_manifests <- function(first, second) {
  second_keys <- paste(second$outcome, second$model_index)
  combined <- first |>
    dplyr::filter(!paste(outcome, model_index) %in% second_keys) |>
    dplyr::bind_rows(second) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  if (
    any(combined$family != "eiv") ||
      anyDuplicated(paste(combined$outcome, combined$model_index)) ||
      any(!file.exists(combined$fit_path))
  ) {
    stop("The combined EIV remedy manifest is invalid.", call. = FALSE)
  }
  combined
}

compute_eiv_full_file_integrated_loo <- function(
    record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L
) {
  candidate <- read_eiv_fit_file(record$fit_path)
  input <- record$input
  if (
    candidate$model_index != input$specification$model_index ||
      !identical(candidate$outcome, input$specification$outcome)
  ) {
    stop("The full EIV fit record is misaligned.", call. = FALSE)
  }
  available_draws <- posterior::ndraws(
    candidate$fit$draws("b_intercept")
  )
  job <- list(
    fitted_candidate = candidate,
    input = input,
    loo_draws = min(as.integer(loo_draws), available_draws),
    fit_draws = available_draws
  )
  result <- compute_eiv_integrated_psis_loo(
    job, k_threshold = k_threshold, cores = cores
  )
  result$fit_record <- record[c(
    "family", "outcome", "model_index", "model_name", "fit_path",
    "stage", "bytes", "sha256"
  )]
  result$input <- input
  result
}

eiv_full_exact_fit_file_path <- function(record) {
  directory <- here::here(
    "_targets", "user", "sensitivity_fit_files", "eiv_exact",
    as.character(record$outcome),
    paste0("model_", stringr::str_pad(record$model_index, 2L, pad = "0"))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, paste0("case_", record$case_id, ".qs"))
}

fit_eiv_full_exact_loo_to_file <- function(
    record,
    model,
    enabled,
    config
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("eiv"))
  if (is_sensitivity_pilot_skip(model)) {
    stop("The EIV Stan model was not compiled.", call. = FALSE)
  }
  existing_path <- eiv_full_exact_fit_file_path(record)
  existing_refit <- if (file.exists(existing_path)) {
    tryCatch(
      read_eiv_full_exact_loo_file(existing_path),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (
    !is.null(existing_refit) &&
      identical(as.integer(existing_refit$case_id), as.integer(record$case_id)) &&
      identical(existing_refit$outcome, record$outcome) &&
      identical(
        as.integer(existing_refit$model_index),
        as.integer(record$model_index)
      ) &&
      identical(
        as.integer(existing_refit$observation),
        as.integer(record$observation)
      )
  ) {
    return(normalizePath(existing_path, winslash = "/", mustWork = TRUE))
  }
  candidate <- read_eiv_fit_file(record$original_fit_path)
  integrated_result <- qs::qread(record$integrated_loo_path)
  input <- integrated_result$input
  if (
    !identical(integrated_result$family, "eiv") ||
      candidate$model_index != record$model_index ||
      input$specification$model_index != record$model_index
  ) {
    stop("The full EIV exact-refit record is invalid.", call. = FALSE)
  }
  training_data <- make_eiv_leave_one_out_stan_data(
    input$stan_data, record$observation
  )
  specification <- candidate$specification
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  output_dir <- make_eiv_full_cmdstan_output_dir(
    paste0("exact_case_", record$case_id),
    record$outcome,
    record$model_index
  )
  start_time <- proc.time()[["elapsed"]]
  fit <- model$sample(
    data = training_data,
    chains = config$chains,
    parallel_chains = config$parallel_chains,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    adapt_delta = specification$adapt_delta,
    max_treedepth = specification$max_treedepth,
    seed = as.integer(config$seed + record$case_id),
    refresh = 100,
    output_dir = output_dir
  )
  fitted_candidate <- list(
    family = "eiv",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    specification = specification,
    grid = input$grid,
    output_dir = output_dir,
    fit = fit
  )
  refit <- list(
    case_id = as.integer(record$case_id),
    family = "eiv",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = as.integer(record$observation),
    pareto_k = as.numeric(record$pareto_k),
    fitted_candidate = fitted_candidate,
    input = input,
    integrated_result = integrated_result,
    refit_seconds = proc.time()[["elapsed"]] - start_time
  )
  write_atomic_qs(refit, existing_path)
}

read_eiv_full_exact_loo_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("The full EIV exact-refit file is unavailable.", call. = FALSE)
  }
  refit <- qs::qread(path)
  output_files <- tryCatch(
    refit$fitted_candidate$fit$output_files(),
    error = function(e) character()
  )
  if (
    !is.list(refit) || !identical(refit$family, "eiv") ||
      !inherits(refit$fitted_candidate$fit, "CmdStanMCMC") ||
      length(output_files) < 1L || any(!file.exists(output_files))
  ) {
    stop("The full EIV exact-refit file is invalid.", call. = FALSE)
  }
  refit
}

diagnose_eiv_full_exact_loo_file <- function(path) {
  summarise_eiv_exact_loo_refit(read_eiv_full_exact_loo_file(path))
}

compute_eiv_full_exact_loo_file <- function(path, diagnostic_audit) {
  if (!isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)) {
    stop("The full EIV exact-refit diagnostic gate has not passed.",
         call. = FALSE)
  }
  compute_eiv_exact_loo_case(read_eiv_full_exact_loo_file(path))
}

make_eiv_full_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    n_draws = 4000L,
    seed = 1988L
) {
  records <- make_mundlak_full_prediction_records(
    accepted_fit_manifest,
    stacking_weights,
    n_draws = n_draws,
    protein_min = 0.8,
    grid_step = 0.02,
    seed = seed
  )
  purrr::map(records, function(record) {
    record$family <- "eiv"
    record$prediction_type <- "eiv"
    record
  })
}

write_eiv_full_prediction_draw_file <- function(record, checkpoint_path) {
  if (
    !identical(record$family, "eiv") || record$draw_count < 1L ||
      !file.exists(record$fit_path) ||
      !identical(
        digest::digest(record$fit_path, file = TRUE, algo = "sha256"),
        record$fit_sha256
      )
  ) {
    stop("A full EIV prediction record is invalid.", call. = FALSE)
  }
  checkpoint <- read_sensitivity_completed_checkpoint(
    checkpoint_path, family = "eiv"
  )
  candidate <- read_eiv_fit_file(record$fit_path)
  if (
    candidate$model_index != record$model_index ||
      !identical(candidate$outcome, record$outcome) ||
      !identical(candidate$model_name, record$model_name) ||
      is.null(checkpoint$analysis_data_by_outcome[[record$outcome]])
  ) {
    stop("A full EIV prediction fit does not match its record.", call. = FALSE)
  }
  all_draws <- posterior::as_draws_matrix(
    candidate$fit$draws("mu_grid")
  )
  set.seed(as.integer(record$prediction_seed))
  selected <- sample(
    seq_len(nrow(all_draws)),
    as.integer(record$draw_count),
    replace = nrow(all_draws) < as.integer(record$draw_count)
  )
  draws <- as.matrix(all_draws[selected, , drop = FALSE])
  grid <- candidate$grid |>
    dplyr::transmute(
      resistance_exercise_code,
      resistance_exercise = ifelse(
        resistance_exercise_code == 0, "YES", "NO"
      ),
      protein_intake = protein_intake_raw
    )
  if (
    nrow(draws) != record$draw_count || ncol(draws) != nrow(grid) ||
      anyNA(draws) || any(!is.finite(draws))
  ) {
    stop("The full EIV prediction draws are incomplete.", call. = FALSE)
  }
  reduction <- list(
    family = "eiv",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    weight = as.numeric(record$weight),
    draw_count = as.integer(record$draw_count),
    prediction_seed = as.integer(record$prediction_seed),
    fit_sha256 = record$fit_sha256,
    grid = grid,
    draws = draws,
    settings = list(
      protein_min = min(grid$protein_intake),
      protein_max = max(grid$protein_intake),
      grid_step = 0.02,
      prediction_type = "eiv",
      reference_intake = NULL
    )
  )
  path <- sensitivity_reduction_file_path(
    family = "eiv",
    outcome = record$outcome,
    model_index = record$model_index,
    stage = "prediction_draws"
  )
  write_atomic_qs(reduction, path)
}
