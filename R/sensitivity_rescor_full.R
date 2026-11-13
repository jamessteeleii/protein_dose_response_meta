# Full residual-correlation sensitivity family -----------------------------

# Study-level residual correlation settings ----

# Specify sampling controls for the expansion and subsequent remedies
# separately so accepted fits can be traced to their settings.
make_rescor_full_expansion_config <- function(
    base_config,
    iter_warmup = 1500L,
    iter_sampling = 2000L,
    seed_offset = 30000L
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
    stop("Invalid full residual-correlation configuration.", call. = FALSE)
  }

  config
}

make_rescor_full_remedy_config <- function(
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
    stop("Invalid residual-correlation remedy configuration.",
         call. = FALSE)
  }

  config
}

make_rescor_full_third_remedy_config <- function(
    second_remedy_config,
    iter_warmup = 4000L,
    iter_sampling = 12000L,
    adapt_delta = 0.99999,
    max_treedepth = 17L,
    seed_offset = 1000L
) {
  config <- second_remedy_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$adapt_delta_floor <- as.numeric(adapt_delta)
  config$max_treedepth_floor <- as.integer(max_treedepth)
  config$seed <- as.integer(second_remedy_config$seed + seed_offset)

  if (
    config$iter_warmup <= second_remedy_config$iter_warmup ||
      config$iter_sampling <= second_remedy_config$iter_sampling ||
      !is.finite(config$adapt_delta_floor) ||
      config$adapt_delta_floor <=
        second_remedy_config$adapt_delta_floor ||
      config$adapt_delta_floor >= 1 ||
      config$max_treedepth_floor <
        second_remedy_config$max_treedepth_floor
  ) {
    stop("Invalid residual-correlation third-remedy configuration.",
         call. = FALSE)
  }

  config
}

make_rescor_full_structural_remedy_config <- function(
    second_remedy_config,
    seed_offset = 3000L
) {
  config <- second_remedy_config
  config$seed <- as.integer(second_remedy_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains < 1L ||
      config$iter_warmup != 3000L ||
      config$iter_sampling != 8000L ||
      config$adapt_delta_floor != 0.9999 ||
      config$max_treedepth_floor != 16L
  ) {
    stop("Invalid residual-correlation structural-remedy configuration.",
         call. = FALSE)
  }
  config
}

make_rescor_observation_exact_remedy_config <- function(
    exact_config,
    iter_warmup = 3000L,
    iter_sampling = 8000L,
    adapt_delta = 0.9999,
    max_treedepth = 16L,
    seed_offset = 1000L
) {
  config <- exact_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$adapt_delta_floor <- as.numeric(adapt_delta)
  config$max_treedepth_floor <- as.integer(max_treedepth)
  config$seed <- as.integer(exact_config$seed + seed_offset)

  if (
    config$iter_warmup <= exact_config$iter_warmup ||
      config$iter_sampling <= exact_config$iter_sampling ||
      !is.finite(config$adapt_delta_floor) ||
      config$adapt_delta_floor <= exact_config$adapt_delta_floor ||
      config$adapt_delta_floor >= 1 ||
      config$max_treedepth_floor < exact_config$max_treedepth_floor
  ) {
    stop("Invalid residual-correlation exact-refit remedy configuration.",
         call. = FALSE)
  }

  config
}

make_rescor_observation_exact_second_remedy_config <- function(
    remedy_config,
    iter_warmup = 4000L,
    iter_sampling = 12000L,
    adapt_delta = 0.99999,
    max_treedepth = 17L,
    seed_offset = 1000L
) {
  config <- remedy_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$adapt_delta_floor <- as.numeric(adapt_delta)
  config$max_treedepth_floor <- as.integer(max_treedepth)
  config$seed <- as.integer(remedy_config$seed + seed_offset)

  if (
    config$iter_warmup <= remedy_config$iter_warmup ||
      config$iter_sampling <= remedy_config$iter_sampling ||
      !is.finite(config$adapt_delta_floor) ||
      config$adapt_delta_floor <= remedy_config$adapt_delta_floor ||
      config$adapt_delta_floor >= 1 ||
      config$max_treedepth_floor < remedy_config$max_treedepth_floor
  ) {
    stop(
      "Invalid residual-correlation exact-refit second-remedy configuration.",
      call. = FALSE
    )
  }

  config
}

# Prepare the correlated likelihood inputs ----

# Translate the analysis rows and candidate basis into Stan arrays with the
# study grouping required by the residual-correlation likelihood.
make_rescor_candidate_data <- function(
    data,
    specification,
    priors,
    rcs_distinct_fallback = c("error", "unique_support")
) {
  rcs_distinct_fallback <- match.arg(rcs_distinct_fallback)
  features <- sensitivity_model_features(specification$model_index)
  required <- c(
    "yi", "vi", "protein_intake", "resistance_exercise_code",
    "study_no", "arm_no", "outcome_no", "age_centre",
    "duration_centre", "protein_intake_raw", "protein_median"
  )
  missing <- setdiff(required, names(data))

  if (length(missing) > 0L) {
    stop(
      "Residual-correlation data are missing: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  data <- data |>
    dplyr::filter(
      is.finite(yi), is.finite(vi), vi > 0,
      is.finite(protein_intake),
      is.finite(resistance_exercise_code)
    )
  moderator_data <- add_primary_standardised_moderators(data)
  moderator_centres <- attr(moderator_data, "moderator_centres")
  moderator_scales <- attr(moderator_data, "moderator_scales")
  data <- moderator_data |>
    dplyr::arrange(study_no, arm_no, outcome_no) |>
    dplyr::mutate(
      study_index = dplyr::dense_rank(study_no),
      arm_index = dplyr::dense_rank(interaction(
        study_no, arm_no, drop = TRUE
      )),
      outcome_index = dplyr::dense_rank(interaction(
        study_no, arm_no, outcome_no, drop = TRUE
      )),
      stan_row = dplyr::row_number()
    )

  arm_data <- data |>
    dplyr::distinct(
      study_no, arm_no, arm_index, protein_intake
    ) |>
    dplyr::arrange(arm_index)
  arm_check <- arm_data |>
    dplyr::count(study_no, arm_no) |>
    dplyr::filter(n != 1L)

  if (nrow(arm_check) > 0L) {
    stop("Study arms do not have unique protein-intake values.",
         call. = FALSE)
  }

  runs <- rle(data$study_index)
  cluster_end <- cumsum(runs$lengths)
  cluster_start <- c(1L, utils::head(cluster_end, -1L) + 1L)
  if (
    length(runs$lengths) != dplyr::n_distinct(data$study_index) ||
      max(runs$lengths) <= 1L
  ) {
    stop("Residual-correlation study clusters are invalid.", call. = FALSE)
  }

  cluster_manifest <- tibble::tibble(
    observation = seq_along(runs$lengths),
    study_no = vapply(
      cluster_start,
      function(index) as.character(data$study_no[[index]]),
      character(1)
    ),
    cluster_start = as.integer(cluster_start),
    cluster_end = as.integer(cluster_end),
    cluster_size = as.integer(runs$lengths)
  )

  rho_lower <- -1 / (max(runs$lengths) - 1) + 1e-8
  rho_raw_prior_location <- stats::qlogis(
    -rho_lower / (1 - rho_lower)
  )

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
      probs = rcs_probs,
      distinct_fallback = rcs_distinct_fallback
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
    family = "rescor",
    specification = specification,
    grid = grid,
    cluster_manifest = cluster_manifest,
    stan_data = list(
      N = nrow(data),
      y = data$yi,
      se_y = sqrt(data$vi),
      x = data$protein_intake,
      N_study = dplyr::n_distinct(data$study_index),
      study = data$study_index,
      N_arm = dplyr::n_distinct(data$arm_index),
      arm = data$arm_index,
      N_outcome = dplyr::n_distinct(data$outcome_index),
      outcome = data$outcome_index,
      N_cluster = length(runs$lengths),
      cluster_start = as.integer(cluster_start),
      cluster_end = as.integer(cluster_end),
      rho_lower = rho_lower,
      rho_raw_prior_location = rho_raw_prior_location,
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

prepare_rescor_full_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome
) {
  outcome <- specification$outcome
  data <- data_by_outcome[[outcome]]
  priors <- priors_by_outcome[[outcome]]
  if (is.null(data) || is.null(priors)) {
    stop("Outcome-specific residual-correlation inputs are unavailable.",
         call. = FALSE)
  }
  make_rescor_candidate_data(data, specification, priors)
}

audit_rescor_full_expansion <- function(inputs, config) {
  keys <- purrr::map_chr(
    inputs,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )
  expected <- as.vector(outer(config$outcomes, seq_len(26L), paste))
  valid <- purrr::map_lgl(inputs, function(input) {
    data <- input$stan_data
    prior_median <- data$rho_lower +
      (1 - data$rho_lower) * stats::plogis(data$rho_raw_prior_location)
    identical(input$family, "rescor") && data$N > 0L &&
      data$N_cluster == data$N_study && data$N_cluster > 1L &&
      nrow(input$cluster_manifest) == data$N_cluster &&
      max(input$cluster_manifest$cluster_end) == data$N &&
      abs(prior_median) < 1e-10 &&
      all(is.finite(data$se_y)) && all(data$se_y > 0)
  })

  if (
    length(inputs) != 104L || anyDuplicated(keys) ||
      !setequal(keys, expected) || !all(valid) ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("The full residual-correlation expansion failed its audit.",
         call. = FALSE)
  }

  tibble::tibble(
    check = c(
      "complete candidate registry",
      "study-cluster covariance inputs",
      "correlation prior",
      "candidate structures and priors",
      "restart-safe initial schedule"
    ),
    status = "PASS",
    detail = c(
      "104 candidates: 26 models x 4 outcomes",
      "one exchangeable sampling-error block per study",
      "common rho is estimated with prior median zero",
      "protein bases, moderators, interactions, and random slopes retained",
      "1,500 warmup and 2,000 retained iterations per chain"
    )
  )
}

rescor_fit_file_path <- function(outcome, model_index, stage) {
  sensitivity_fit_file_path("rescor", outcome, model_index, stage)
}

read_rescor_fit_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("The residual-correlation fit file is unavailable.", call. = FALSE)
  }
  candidate <- qs::qread(path)
  output_files <- tryCatch(
    candidate$fit$output_files(),
    error = function(e) character()
  )

  if (
    !is.list(candidate) || !identical(candidate$family, "rescor") ||
      !inherits(candidate$fit, "CmdStanMCMC") ||
      is.null(candidate$outcome) || is.null(candidate$model_index) ||
      length(output_files) < 1L || any(!file.exists(output_files))
  ) {
    stop(
      "The residual-correlation fit or its CmdStan CSV files are invalid.",
      call. = FALSE
    )
  }

  candidate
}

make_rescor_full_cmdstan_output_dir <- function(stage, outcome, model_index) {
  safe <- function(x) gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  path <- here::here(
    "_targets", "user", "sensitivity_cmdstan", "rescor_full",
    safe(stage), safe(outcome),
    paste0("model_", stringr::str_pad(model_index, 2L, pad = "0"))
  )
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# Fit a prepared candidate only when enabled and return its saved artifact for
# the diagnostic and acceptance stages.
fit_rescor_full_candidate_to_file <- function(
    input,
    model,
    enabled,
    config,
    stage = "initial"
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("rescor"))

  if (is_sensitivity_pilot_skip(model)) {
    stop("The full residual-correlation Stan model was not compiled.",
         call. = FALSE)
  }
  specification <- input$specification
  outcome_offset <- match(specification$outcome, config$outcomes)

  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the residual-correlation fit.", call. = FALSE)
  }
  existing_path <- rescor_fit_file_path(
    specification$outcome, specification$model_index, stage
  )
  existing <- if (file.exists(existing_path)) {
    tryCatch(read_rescor_fit_file(existing_path), error = function(e) NULL)
  } else {
    NULL
  }

  if (
    !is.null(existing) &&
      identical(existing$outcome, specification$outcome) &&
      identical(
        as.integer(existing$model_index),
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
  output_dir <- make_rescor_full_cmdstan_output_dir(
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
    family = "rescor",
    outcome = specification$outcome,
    model_index = as.integer(specification$model_index),
    model_name = paste0(
      specification$model_name,
      " [residual correlation by study]"
    ),
    specification = specification,
    grid = input$grid,
    output_dir = output_dir,
    fit = fit
  )

  write_atomic_qs(candidate, existing_path)
}

make_rescor_fit_file_manifest <- function(paths, stage) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) == 0L) {
    return(tibble::tibble(
      family = character(), outcome = character(), model_index = integer(),
      model_name = character(), fit_path = character(), stage = character(),
      bytes = double(), sha256 = character()
    ))
  }

  if (anyNA(paths) || any(!file.exists(paths))) {
    stop("The residual-correlation fit set is incomplete.", call. = FALSE)
  }
  manifest <- purrr::map_dfr(paths, function(path) {
    candidate <- read_rescor_fit_file(path)
    tibble::tibble(
      family = "rescor",
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
    stop("The residual-correlation manifest contains duplicate models.",
         call. = FALSE)
  }
  manifest
}

make_rescor_full_fit_records <- function(fit_manifest, inputs) {
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
    stop("The residual-correlation fit/input registry is misaligned.",
         call. = FALSE)
  }

  purrr::map(seq_len(nrow(fit_manifest)), function(row) {
    record <- purrr::map(as.list(fit_manifest[row, , drop = FALSE]), unname)
    record$input <- inputs[[match(manifest_keys[[row]], input_keys)]]
    record
  })
}

make_rescor_full_accepted_fit_records <- function(fit_manifest, inputs) {
  records <- make_rescor_full_fit_records(fit_manifest, inputs)
  purrr::map(records, function(record) {
    if (identical(record$stage, "remedy_4_uncorrelated")) {
      if (
        !identical(record$outcome, "function") ||
          record$model_index != 23L ||
          isTRUE(
            record$input$specification$uncorrelated_study_effects
          )
      ) {
        stop("The accepted structural-remedy record is invalid.",
             call. = FALSE)
      }
      record$input$specification$uncorrelated_study_effects <- TRUE
      record$input$specification$newly_uncorrelated_study_effects <- TRUE
      record$input$specification$rescor_structural_remedy <-
        "independent_study_random_smooth_coefficients"
      record$input$stan_data$correlated_study_effects <- 0L
    }
    record
  })
}

compile_rescor_full_diagnostic_audit <- function(
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
    nrow(summary) != expected_models || any(summary$family != "rescor") ||
      anyDuplicated(paste(summary$outcome, summary$model_index))
  ) {
    stop("The residual-correlation diagnostic audit is incomplete.",
         call. = FALSE)
  }

  list(
    family = "rescor",
    diagnostic_summary = summary,
    parameter_failures = tibble::tibble(),
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

select_rescor_full_remedy_inputs <- function(inputs, audit) {
  failed_keys <- audit$diagnostic_summary |>
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
    !identical(audit$family, "rescor") || anyDuplicated(failed_keys) ||
      length(selected) != length(failed_keys) ||
      !setequal(selected_keys, failed_keys)
  ) {
    stop("The residual-correlation remedy inputs are invalid.",
         call. = FALSE)
  }

  selected
}

select_rescor_full_uncorrelated_remedy_inputs <- function(inputs, audit) {
  selected <- select_rescor_full_remedy_inputs(inputs, audit)
  keys <- purrr::map_chr(
    selected,
    ~ paste(.x$specification$outcome, .x$specification$model_index)
  )

  if (
    length(selected) != 1L ||
      !identical(unname(keys), "function 23") ||
      isTRUE(selected[[1L]]$specification$uncorrelated_study_effects) ||
      selected[[1L]]$stan_data$include_random_slope != 1L
  ) {
    stop(
      "The targeted residual-correlation structural remedy has changed.",
      call. = FALSE
    )
  }

  purrr::map(selected, function(input) {
    input$specification$uncorrelated_study_effects <- TRUE
    input$specification$newly_uncorrelated_study_effects <- TRUE
    input$specification$rescor_structural_remedy <-
      "independent_study_random_smooth_coefficients"
    input$stan_data$correlated_study_effects <- 0L
    input
  })
}

combine_rescor_full_remedy_diagnostic_audit <- function(
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
    stop("Residual-correlation remedy diagnostics are incomplete.",
         call. = FALSE)
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
    family = "rescor",
    diagnostic_summary = summary,
    parameter_failures = tibble::tibble(),
    models = nrow(summary),
    models_passing_sampling = sum(summary$diagnostic_pass),
    models_requiring_remediation = sum(!summary$diagnostic_pass),
    all_sampling_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

combine_rescor_remedy_fit_manifests <- function(first, second) {
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
    any(combined$family != "rescor") ||
      anyDuplicated(paste(combined$outcome, combined$model_index)) ||
      any(!file.exists(combined$fit_path))
  ) {
    stop("The combined residual-correlation remedy manifest is invalid.",
         call. = FALSE)
  }

  combined
}

combine_rescor_full_accepted_fit_manifest <- function(
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
      nrow(accepted) != 104L || any(accepted$family != "rescor") ||
      anyDuplicated(paste(accepted$outcome, accepted$model_index)) ||
      any(!file.exists(accepted$fit_path))
  ) {
    stop("The accepted residual-correlation manifest is invalid.",
         call. = FALSE)
  }

  accepted
}

# Study-cluster and observation-level scoring ----

# The study-cluster score and the observation-level score have different
# holdout units. Use the observation-level results when comparing against
# candidates scored per observation.
compute_rescor_full_file_integrated_loo <- function(
    record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L
) {
  candidate <- read_rescor_fit_file(record$fit_path)
  input <- record$input

  if (
    candidate$model_index != input$specification$model_index ||
      !identical(candidate$outcome, input$specification$outcome)
  ) {
    stop("The residual-correlation LOO record is misaligned.",
         call. = FALSE)
  }

  available <- posterior::ndraws(candidate$fit$draws("b_intercept"))
  draw_ids <- make_sensitivity_loo_draw_ids(
    available, min(as.integer(loo_draws), available)
  )
  conditional_log_lik <- posterior::as_draws_matrix(
    candidate$fit$draws("log_lik_cluster")
  )[draw_ids, , drop = FALSE]
  integrated_log_lik <- posterior::as_draws_matrix(
    candidate$fit$draws("log_lik_cluster_integrated")
  )[draw_ids, , drop = FALSE]

  if (
    ncol(integrated_log_lik) != nrow(input$cluster_manifest) ||
      anyNA(integrated_log_lik) || any(!is.finite(integrated_log_lik))
  ) {
    stop("Integrated residual-correlation log likelihood is invalid.",
         call. = FALSE)
  }

  conditional_loo <- loo::loo(conditional_log_lik, cores = as.integer(cores))
  integrated_loo <- loo::loo(integrated_log_lik, cores = as.integer(cores))
  list(
    family = "rescor",
    outcome = candidate$outcome,
    model_index = as.integer(candidate$model_index),
    model_name = candidate$model_name,
    loo = integrated_loo,
    conditional_loo = conditional_loo,
    observation_manifest = input$cluster_manifest,
    remedy = "study_cluster_integrated_psis",
    k_threshold = k_threshold,
    loo_draws = length(draw_ids),
    fit_record = record[c(
      "family", "outcome", "model_index", "model_name", "fit_path",
      "stage", "bytes", "sha256"
    )],
    input = input
  )
}

rescor_full_gq_parameter_variables <- function() {
  c(
    "b_intercept", "b_rt", "b_protein", "b_protein_rt",
    "b_moderator", "z_study", "sd_study", "L_study",
    "z_arm", "z_outcome", "sd_arm", "sd_outcome", "rho_raw"
  )
}

select_rescor_full_gq_draws <- function(fit, loo_draws = 4000L) {
  draws <- fit$draws(
    variables = rescor_full_gq_parameter_variables(),
    format = "draws_array"
  )
  dimensions <- dim(draws)

  if (
    length(dimensions) != 3L || dimensions[[1L]] < 1L ||
      dimensions[[2L]] < 1L
  ) {
    stop("Residual-correlation posterior draws are invalid.",
         call. = FALSE)
  }

  requested <- min(as.integer(loo_draws), prod(dimensions[1:2]))
  draws_per_chain <- requested %/% dimensions[[2L]]
  if (draws_per_chain < 1L) {
    stop("Too few draws for the residual-correlation LOO audit.",
         call. = FALSE)
  }

  draw_ids <- make_sensitivity_loo_draw_ids(
    dimensions[[1L]],
    min(draws_per_chain, dimensions[[1L]])
  )
  posterior::as_draws_array(draws[draw_ids, , , drop = FALSE])
}

make_rescor_full_observation_manifest <- function(input) {
  data <- input$stan_data
  clusters <- input$cluster_manifest
  cluster_sizes <- clusters$cluster_end - clusters$cluster_start + 1L
  if (
    length(cluster_sizes) != data$N_cluster ||
      sum(cluster_sizes) != data$N ||
      any(cluster_sizes != clusters$cluster_size)
  ) {
    stop("Residual-correlation observation clusters are misaligned.",
         call. = FALSE)
  }

  cluster <- rep(seq_len(data$N_cluster), cluster_sizes)
  tibble::tibble(
    observation = seq_len(data$N),
    study_cluster = as.integer(cluster),
    study_no = clusters$study_no[cluster],
    cluster_size = as.integer(cluster_sizes[cluster]),
    study_index = as.integer(data$study),
    arm_index = as.integer(data$arm),
    outcome_index = as.integer(data$outcome),
    yi = as.numeric(data$y),
    vi = as.numeric(data$se_y)^2
  )
}

compute_rescor_full_file_observation_loo <- function(
    record,
    generated_quantities_model,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    parallel_chains = 4L,
    seed = 39991L
) {
  if (is_sensitivity_pilot_skip(generated_quantities_model)) {
    stop("The residual-correlation LOO model was not compiled.",
         call. = FALSE)
  }

  candidate <- read_rescor_fit_file(record$fit_path)
  input <- record$input
  if (
    candidate$model_index != input$specification$model_index ||
      !identical(candidate$outcome, input$specification$outcome)
  ) {
    stop("The residual-correlation observation-LOO record is misaligned.",
         call. = FALSE)
  }

  selected_draws <- select_rescor_full_gq_draws(
    candidate$fit, loo_draws = loo_draws
  )
  temporary_directory <- tempfile(
    pattern = paste0(
      "rescor-observation-loo-", candidate$outcome, "-",
      candidate$model_index, "-"
    )
  )
  dir.create(temporary_directory, recursive = TRUE, showWarnings = FALSE)
  on.exit(
    unlink(temporary_directory, recursive = TRUE, force = TRUE),
    add = TRUE
  )
  fitted_parameter_files <- cmdstanr::draws_to_csv(
    selected_draws,
    dir = temporary_directory,
    basename = "fitted_parameters"
  )
  outcome_offset <- match(
    candidate$outcome,
    c("whole_body", "hypertrophy", "strength", "function")
  )
  generated <- generated_quantities_model$generate_quantities(
    fitted_params = fitted_parameter_files,
    data = input$stan_data,
    seed = as.integer(seed + 100L * outcome_offset + candidate$model_index),
    output_dir = temporary_directory,
    output_basename = paste0("model_", candidate$model_index),
    parallel_chains = min(
      as.integer(parallel_chains), dim(selected_draws)[[2L]]
    )
  )

  observation_log_lik <- posterior::as_draws_matrix(
    generated$draws(
      "log_lik_observation_conditional_integrated",
      format = "draws_array"
    )
  )
  observation_manifest <- make_rescor_full_observation_manifest(input)

  if (
    nrow(observation_log_lik) !=
      dim(selected_draws)[[1L]] * dim(selected_draws)[[2L]] ||
      ncol(observation_log_lik) != nrow(observation_manifest) ||
      anyNA(observation_log_lik) ||
      any(!is.finite(observation_log_lik))
  ) {
    stop("Conditional observation-level log likelihood is invalid.",
         call. = FALSE)
  }

  observation_loo <- loo::loo(
    observation_log_lik, cores = as.integer(cores)
  )
  list(
    family = "rescor",
    outcome = candidate$outcome,
    model_index = as.integer(candidate$model_index),
    model_name = candidate$model_name,
    loo = observation_loo,
    observation_manifest = observation_manifest,
    remedy = "correlated_conditional_observation_level_integrated_psis",
    loo_unit = "observation_conditional_on_study",
    k_threshold = as.numeric(k_threshold),
    loo_draws = nrow(observation_log_lik),
    fit_record = record[c(
      "family", "outcome", "model_index", "model_name", "fit_path",
      "stage", "bytes", "sha256"
    )],
    input = input
  )
}

write_rescor_full_observation_loo_file <- function(result) {
  if (
    !inherits(result$loo, "loo") ||
      !identical(result$family, "rescor") ||
      !identical(
        result$remedy,
        "correlated_conditional_observation_level_integrated_psis"
      ) ||
      !identical(result$loo_unit, "observation_conditional_on_study") ||
      nrow(result$observation_manifest) != nrow(result$loo$pointwise)
  ) {
    stop("Cannot archive an invalid observation-level LOO result.",
         call. = FALSE)
  }

  path <- sensitivity_reduction_file_path(
    family = "rescor",
    outcome = result$outcome,
    model_index = result$model_index,
    stage = "observation_conditional_loo"
  )
  write_atomic_qs(result, path)
}

audit_rescor_full_observation_loo_files <- function(
    manifest,
    k_threshold = 0.7
) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 104L ||
      any(manifest$family != "rescor") ||
      any(!file.exists(manifest$integrated_loo_path))
  ) {
    stop("The observation-level residual-correlation LOO files cannot be audited.",
         call. = FALSE)
  }

  results <- purrr::map(manifest$integrated_loo_path, qs::qread)
  valid <- purrr::map_lgl(results, function(result) {
    identical(result$family, "rescor") &&
      identical(result$loo_unit, "observation_conditional_on_study") &&
      inherits(result$loo, "loo") &&
      nrow(result$observation_manifest) == nrow(result$loo$pointwise)
  })
  keys <- purrr::map_chr(
    results, ~ paste(.x$outcome, .x$model_index)
  )
  expected_keys <- as.vector(outer(
    c("whole_body", "hypertrophy", "strength", "function"),
    seq_len(26L),
    paste
  ))

  if (!all(valid) || anyDuplicated(keys) || !setequal(keys, expected_keys)) {
    stop("The observation-level residual-correlation LOO set is invalid.",
         call. = FALSE)
  }
  method <- "correlated_conditional_observation_level_psis"
  model_summary <- purrr::map_dfr(results, function(result) {
    summarise_primary_loo_method(result$loo, method, k_threshold) |>
      dplyr::mutate(
        outcome = result$outcome,
        model_index = result$model_index,
        model_name = result$model_name,
        loo_draws = result$loo_draws,
        .before = 1L
      )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )

  high_k_cases <- purrr::map_dfr(results, function(result) {
    pareto_k <- as.numeric(loo::pareto_k_values(result$loo))
    observations <- which(pareto_k > k_threshold)
    if (length(observations) < 1L) return(tibble::tibble())
    result$observation_manifest[observations, , drop = FALSE] |>
      dplyr::mutate(
        outcome = result$outcome,
        model_index = result$model_index,
        model_name = result$model_name,
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
    models = 104L,
    exact_refits_required = nrow(high_k_cases),
    observation_psis_complete = TRUE,
    stacking_ready = nrow(high_k_cases) == 0L,
    loo_unit = "observation_conditional_on_study",
    study_level_audit_retained = TRUE,
    k_threshold = as.numeric(k_threshold)
  )
}

make_rescor_full_observation_exact_records <- function(
    observation_audit,
    accepted_fit_manifest,
    observation_loo_manifest
) {
  cases <- observation_audit$high_k_cases
  expected <- as.integer(observation_audit$exact_refits_required)

  if (
    !identical(observation_audit$loo_unit,
               "observation_conditional_on_study") ||
      !is.data.frame(cases) || nrow(cases) != expected || expected < 1L
  ) {
    stop("The residual-correlation observation exact cases are invalid.",
         call. = FALSE)
  }

  fit_lookup <- accepted_fit_manifest |>
    dplyr::select(
      outcome, model_index,
      original_fit_path = fit_path,
      original_fit_sha256 = sha256
    )
  loo_lookup <- observation_loo_manifest |>
    dplyr::select(
      outcome, model_index,
      observation_loo_path = integrated_loo_path,
      observation_loo_sha256 = sha256
    )
  records <- cases |>
    dplyr::select(
      outcome, model_index, model_name, observation, pareto_k,
      study_cluster, study_no, cluster_size
    ) |>
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

  valid_hash <- function(path, expected_hash) {
    file.exists(path) && identical(
      digest::digest(path, file = TRUE, algo = "sha256"),
      expected_hash
    )
  }
  hashes_valid <- vapply(
    seq_len(nrow(records)),
    function(index) {
      valid_hash(
        records$original_fit_path[[index]],
        records$original_fit_sha256[[index]]
      ) && valid_hash(
        records$observation_loo_path[[index]],
        records$observation_loo_sha256[[index]]
      )
    },
    logical(1)
  )

  if (
    nrow(records) != expected ||
      anyDuplicated(paste(
        records$outcome, records$model_index, records$observation
      )) ||
      anyNA(records) || !all(hashes_valid)
  ) {
    stop("The residual-correlation exact-record registry is invalid.",
         call. = FALSE)
  }

  purrr::map(seq_len(nrow(records)), function(index) {
    purrr::map(as.list(records[index, , drop = FALSE]), unname)
  })
}

# Exact observation holdouts ----

# Leave out an individual observation while preserving its study context for
# the conditional predictive density.
make_rescor_leave_one_observation_out_data <- function(
    stan_data,
    observation
) {
  observation <- as.integer(observation)
  if (
    length(observation) != 1L || is.na(observation) ||
      observation < 1L || observation > stan_data$N
  ) {
    stop("Invalid held-out residual-correlation observation.",
         call. = FALSE)
  }

  keep <- setdiff(seq_len(stan_data$N), observation)
  row_fields <- c(
    "y", "se_y", "x", "study", "arm", "outcome",
    "resistance_exercise_code", "age_z", "duration_z"
  )
  if (
    !all(row_fields %in% names(stan_data)) ||
      any(vapply(stan_data[row_fields], length, integer(1)) != stan_data$N)
  ) {
    stop("Residual-correlation row data are incomplete.", call. = FALSE)
  }

  for (field in row_fields) {
    stan_data[[field]] <- stan_data[[field]][keep]
  }
  remap <- function(values) match(values, unique(values))
  stan_data$study <- remap(stan_data$study)
  stan_data$arm <- remap(stan_data$arm)
  stan_data$outcome <- remap(stan_data$outcome)
  stan_data$N <- length(keep)
  stan_data$N_study <- max(stan_data$study)
  stan_data$N_arm <- max(stan_data$arm)
  stan_data$N_outcome <- max(stan_data$outcome)
  runs <- rle(stan_data$study)
  stan_data$N_cluster <- length(runs$lengths)
  stan_data$cluster_end <- as.integer(cumsum(runs$lengths))
  stan_data$cluster_start <- as.integer(c(
    1L, utils::head(stan_data$cluster_end, -1L) + 1L
  ))
  # Preserve rho_lower and rho_raw_prior_location so the exact refit uses
  # precisely the same residual-correlation prior as its parent model.
  stan_data
}

rescor_observation_exact_result_file_path <- function(
    record,
    stage = "observation_exact"
) {
  sensitivity_reduction_file_path(
    family = "rescor",
    outcome = record$outcome,
    model_index = record$model_index,
    stage = stage,
    suffix = paste0(
      "case_", record$case_id, "_row_", record$observation
    )
  )
}

rescor_held_out_conditional_integrated_log_lik <- function(
    candidate,
    input,
    observation,
    loo_draws = 4000L
) {
  data <- input$stan_data
  observation <- as.integer(observation)
  cluster <- which(
    observation >= data$cluster_start &
      observation <= data$cluster_end
  )

  if (length(cluster) != 1L) {
    stop("The held-out observation does not identify one study cluster.",
         call. = FALSE)
  }
  rows <- seq.int(
    data$cluster_start[[cluster]], data$cluster_end[[cluster]]
  )
  held_position <- match(observation, rows)
  variables <- c(
    "b_intercept", "b_rt", "b_protein", "b_protein_rt",
    "b_moderator", "sd_study", "L_study", "sd_arm",
    "sd_outcome", "rho"
  )
  all_draws <- posterior::as_draws_df(candidate$fit$draws(variables))
  draw_ids <- make_sensitivity_loo_draw_ids(
    nrow(all_draws), min(as.integer(loo_draws), nrow(all_draws))
  )
  draws <- as.data.frame(all_draws)[draw_ids, , drop = FALSE]
  n_draws <- nrow(draws)
  k_basis <- as.integer(data$K_basis)
  k_study <- as.integer(data$K_study)
  basis_list <- evaluate_eiv_stan_basis(
    matrix(data$x[rows], nrow = 1L), data
  )
  basis <- do.call(cbind, lapply(basis_list, as.numeric))
  fixed_mean <- outer(
    as.numeric(draws$b_intercept), rep(1, length(rows))
  ) + outer(
    as.numeric(draws$b_rt), data$resistance_exercise_code[rows]
  )

  if (data$include_protein == 1L) {
    coefficient <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_protein", k_basis
    )
    fixed_mean <- fixed_mean + coefficient %*% t(basis)
  }

  if (data$include_rt_interaction == 1L) {
    coefficient <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_protein_rt", k_basis
    )
    fixed_mean <- fixed_mean +
      (coefficient %*% t(basis)) *
      matrix(
        data$resistance_exercise_code[rows],
        nrow = n_draws,
        ncol = length(rows),
        byrow = TRUE
      )
  }

  if (data$include_moderators == 1L) {
    moderators <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_moderator", 2L
    )
    fixed_mean <- fixed_mean +
      outer(moderators[, 1L], data$duration_z[rows]) +
      outer(moderators[, 2L], data$age_z[rows])
  }

  study_design <- if (data$include_random_slope == 1L) {
    cbind(1, basis)
  } else {
    matrix(1, nrow = length(rows), ncol = 1L)
  }
  study_sd <- extract_cmdstan_indexed_draw_matrix(
    draws, "sd_study", k_study
  )
  study_cholesky <- if (data$correlated_study_effects == 1L) {
    extract_cmdstan_matrix_draw_array(
      draws, "L_study", k_study, k_study
    )
  } else {
    NULL
  }
  arm_match <- outer(data$arm[rows], data$arm[rows], "==")
  outcome_match <- outer(
    data$outcome[rows], data$outcome[rows], "=="
  )
  residual_scale <- outer(data$se_y[rows], data$se_y[rows])
  observed <- data$y[rows]
  held_out_log_lik <- numeric(n_draws)
  held_basis <- numeric(length(rows))
  held_basis[[held_position]] <- 1

  for (draw in seq_len(n_draws)) {
    study_covariance <- if (data$correlated_study_effects == 1L) {
      scaled_cholesky <- diag(study_sd[draw, ], nrow = k_study) %*%
        study_cholesky[draw, , ]
      tcrossprod(scaled_cholesky)
    } else {
      diag(study_sd[draw, ]^2, nrow = k_study)
    }
    residual_covariance <- residual_scale * as.numeric(draws$rho[[draw]])
    diag(residual_covariance) <- data$se_y[rows]^2
    covariance <- residual_covariance +
      study_design %*% study_covariance %*% t(study_design) +
      arm_match * as.numeric(draws$sd_arm[[draw]])^2 +
      outcome_match * as.numeric(draws$sd_outcome[[draw]])^2
    cholesky <- chol(covariance)
    centred <- observed - fixed_mean[draw, ]
    precision_residual <- backsolve(
      cholesky,
      forwardsolve(t(cholesky), centred)
    )
    precision_column <- backsolve(
      cholesky,
      forwardsolve(t(cholesky), held_basis)
    )
    precision_diagonal <- precision_column[[held_position]]
    conditional_variance <- 1 / precision_diagonal
    conditional_residual <-
      precision_residual[[held_position]] / precision_diagonal
    held_out_log_lik[[draw]] <- stats::dnorm(
      conditional_residual,
      mean = 0,
      sd = sqrt(conditional_variance),
      log = TRUE
    )
  }

  if (
    anyNA(held_out_log_lik) || any(!is.finite(held_out_log_lik))
  ) {
    stop("The exact conditional held-out likelihood is invalid.",
         call. = FALSE)
  }

  held_out_log_lik
}

make_rescor_observation_exact_correction <- function(
    record,
    candidate,
    integrated_result,
    diagnostic,
    refit_seconds,
    loo_draws = 4000L
) {
  held_out_log_lik <- rescor_held_out_conditional_integrated_log_lik(
    candidate,
    integrated_result$input,
    record$observation,
    loo_draws = loo_draws
  )
  maximum_log_lik <- max(held_out_log_lik)
  scaled_likelihood <- exp(held_out_log_lik - maximum_log_lik)
  exact_elpd <- maximum_log_lik + log(mean(scaled_likelihood))
  exact_mcse <- tryCatch(
    posterior::mcse_mean(scaled_likelihood) / mean(scaled_likelihood),
    error = function(e) NA_real_
  )
  patched_loo <- integrated_result$loo
  observation <- as.integer(record$observation)
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
  exact_refit_summary <- diagnostic |>
    dplyr::mutate(
      case_id = as.integer(record$case_id),
      observation = observation,
      old_psis_pareto_k = as.numeric(record$pareto_k),
      refit_seconds = as.numeric(refit_seconds),
      held_out_study = as.character(record$study_no),
      held_out_cluster_size = as.integer(record$cluster_size),
      old_psis_elpd_loo = old_elpd,
      exact_elpd_loo = exact_elpd,
      exact_minus_psis_elpd = exact_elpd - old_elpd,
      exact_elpd_mcse = exact_mcse,
      .before = 1L
    )

  list(
    case_id = as.integer(record$case_id),
    family = "rescor",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation_manifest = integrated_result$observation_manifest,
    loo = patched_loo,
    remedy =
      "correlated_conditional_observation_level_psis_plus_exact",
    loo_unit = "observation_conditional_on_study",
    k_threshold = integrated_result$k_threshold,
    exact_refit_summary = exact_refit_summary
  )
}

fit_and_score_rescor_observation_exact_to_file <- function(
    record,
    model,
    enabled,
    config,
    loo_draws = 4000L,
    stage = "observation_exact",
    result_path = NULL
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("rescor"))

  if (is_sensitivity_pilot_skip(model)) {
    stop("The residual-correlation Stan model was not compiled.",
         call. = FALSE)
  }
  output_path <- if (is.null(result_path)) {
    rescor_observation_exact_result_file_path(record, stage)
  } else {
    if (
      length(result_path) != 1L || is.na(result_path) ||
        !nzchar(result_path)
    ) {
      stop("The observation exact-refit result path is invalid.",
           call. = FALSE)
    }
    dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
    result_path
  }

  config_signature <- list(
    stage = as.character(stage),
    chains = as.integer(config$chains),
    iter_warmup = as.integer(config$iter_warmup),
    iter_sampling = as.integer(config$iter_sampling),
    adapt_delta_floor = as.numeric(config$adapt_delta_floor),
    max_treedepth_floor = as.integer(config$max_treedepth_floor),
    loo_draws = as.integer(loo_draws)
  )

  if (file.exists(output_path)) {
    existing <- tryCatch(qs::qread(output_path), error = function(e) NULL)
    existing_stage <- if (
      is.list(existing) && is.null(existing$stage)
    ) {
      "observation_exact"
    } else if (is.list(existing)) {
      existing$stage
    } else {
      NULL
    }
    signature_fields <- setdiff(names(config_signature), "stage")
    signature_matches <- is.list(existing) &&
      is.list(existing$config_signature) &&
      identical(
        existing$config_signature[signature_fields],
        config_signature[signature_fields]
      )
    if (
      is.list(existing) && identical(existing$family, "rescor") &&
        identical(existing$case_id, as.integer(record$case_id)) &&
        identical(existing$outcome, record$outcome) &&
        identical(existing$model_index, as.integer(record$model_index)) &&
        identical(existing$observation, as.integer(record$observation)) &&
        identical(existing_stage, as.character(stage)) &&
        signature_matches &&
        identical(existing$original_fit_sha256,
                  record$original_fit_sha256) &&
        identical(existing$observation_loo_sha256,
                  record$observation_loo_sha256)
    ) {
      return(normalizePath(output_path, winslash = "/", mustWork = TRUE))
    }
  }

  original <- read_rescor_fit_file(record$original_fit_path)
  integrated_result <- qs::qread(record$observation_loo_path)
  if (
    original$model_index != record$model_index ||
      !identical(original$outcome, record$outcome) ||
      !identical(integrated_result$loo_unit,
                 "observation_conditional_on_study") ||
      integrated_result$model_index != record$model_index
  ) {
    stop("The observation exact-refit record is misaligned.",
         call. = FALSE)
  }

  training_data <- make_rescor_leave_one_observation_out_data(
    integrated_result$input$stan_data,
    record$observation
  )
  specification <- integrated_result$input$specification
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  temporary_root <- here::here(
    "_targets", "user", "rescor_observation_exact_temporary"
  )
  dir.create(temporary_root, recursive = TRUE, showWarnings = FALSE)
  stale_directories <- list.dirs(
    temporary_root, recursive = FALSE, full.names = TRUE
  )
  stale_directories <- stale_directories[startsWith(
    basename(stale_directories), paste0("case_", record$case_id, "_")
  )]

  if (length(stale_directories)) {
    unlink(stale_directories, recursive = TRUE, force = TRUE)
  }
  output_dir <- tempfile(
    pattern = paste0("case_", record$case_id, "_"),
    tmpdir = temporary_root
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
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
    refresh = 500,
    output_dir = output_dir
  )

  refit_seconds <- proc.time()[["elapsed"]] - start_time
  candidate <- list(
    family = "rescor",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    specification = specification,
    grid = integrated_result$input$grid,
    output_dir = output_dir,
    fit = fit
  )

  diagnostic <- diagnose_cmdstan_pilot_fit(candidate)
  correction <- if (isTRUE(diagnostic$diagnostic_pass[[1L]])) {
    make_rescor_observation_exact_correction(
      record,
      candidate,
      integrated_result,
      diagnostic,
      refit_seconds,
      loo_draws = loo_draws
    )
  } else {
    NULL
  }

  diagnostic <- diagnostic |>
    dplyr::mutate(
      case_id = as.integer(record$case_id),
      observation = as.integer(record$observation),
      old_psis_pareto_k = as.numeric(record$pareto_k),
      refit_seconds = as.numeric(refit_seconds),
      .before = 1L
    )

  compact_result <- list(
    case_id = as.integer(record$case_id),
    family = "rescor",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = as.integer(record$observation),
    stage = as.character(stage),
    diagnostic = diagnostic,
    correction = correction,
    config_signature = config_signature,
    original_fit_sha256 = record$original_fit_sha256,
    observation_loo_sha256 = record$observation_loo_sha256
  )

  write_atomic_qs(compact_result, output_path)
}

read_rescor_observation_exact_result_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("A residual-correlation observation exact result is unavailable.",
         call. = FALSE)
  }
  result <- qs::qread(path)
  if (
    !is.list(result) || !identical(result$family, "rescor") ||
      !is.data.frame(result$diagnostic) ||
      nrow(result$diagnostic) != 1L ||
      !identical(result$case_id,
                 as.integer(result$diagnostic$case_id[[1L]]))
  ) {
    stop("A residual-correlation observation exact result is invalid.",
         call. = FALSE)
  }
  result
}

select_rescor_observation_exact_remedy_records <- function(
    records,
    diagnostic_audit
) {
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::select(case_id, outcome, model_index, observation)
  failed_ids <- as.integer(failed$case_id)
  selected <- purrr::keep(
    records,
    ~ as.integer(.x$case_id) %in% failed_ids
  )
  selected_summary <- purrr::map_dfr(selected, function(record) {
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      observation = as.integer(record$observation)
    )
  }) |>
    dplyr::arrange(case_id)

  failed <- dplyr::arrange(failed, case_id)
  if (
    nrow(failed) < 1L || length(selected) != nrow(failed) ||
      !identical(selected_summary, failed)
  ) {
    stop("The observation exact-refit remedy selection is invalid.",
         call. = FALSE)
  }

  selected
}

combine_rescor_observation_exact_result_files <- function(
    initial_paths,
    remedy_paths,
    expected_cases
) {
  initial_paths <- unname(unlist(initial_paths, use.names = FALSE))
  remedy_paths <- unname(unlist(remedy_paths, use.names = FALSE))

  if (
    length(initial_paths) != expected_cases ||
      length(remedy_paths) < 1L ||
      any(!file.exists(c(initial_paths, remedy_paths)))
  ) {
    stop("The observation exact-refit result files are incomplete.",
         call. = FALSE)
  }

  initial <- purrr::map(
    initial_paths, read_rescor_observation_exact_result_file
  )
  remedies <- purrr::map(
    remedy_paths, read_rescor_observation_exact_result_file
  )
  remedy_ids <- purrr::map_int(remedies, "case_id")
  combined <- c(
    purrr::discard(
      initial,
      ~ as.integer(.x$case_id) %in% remedy_ids
    ),
    remedies
  )
  combined <- combined[order(purrr::map_int(combined, "case_id"))]
  case_ids <- purrr::map_int(combined, "case_id")

  if (
    length(combined) != expected_cases || anyDuplicated(case_ids) ||
      !identical(unname(case_ids), seq_len(expected_cases))
  ) {
    stop("The combined observation exact-refit set is invalid.",
         call. = FALSE)
  }

  purrr::map(combined, function(result) {
    normalizePath(
      if (result$case_id %in% remedy_ids) {
        remedy_paths[[match(result$case_id, remedy_ids)]]
      } else {
        initial_paths[[match(
          result$case_id, purrr::map_int(initial, "case_id")
        )]]
      },
      winslash = "/",
      mustWork = TRUE
    )
  })
}

diagnose_rescor_observation_exact_result_file <- function(path) {
  read_rescor_observation_exact_result_file(path)$diagnostic
}

collect_rescor_observation_exact_corrections <- function(
    paths,
    records,
    diagnostic_audit
) {
  if (!isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)) {
    stop("Observation exact corrections require passing diagnostics.",
         call. = FALSE)
  }

  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != length(records) || any(!file.exists(paths))) {
    stop("The observation exact-result file set is incomplete.",
         call. = FALSE)
  }
  results <- purrr::map(paths, read_rescor_observation_exact_result_file)
  results <- results[order(purrr::map_int(results, "case_id"))]
  expected <- sort(purrr::map_int(records, "case_id"))

  if (
    !identical(unname(purrr::map_int(results, "case_id")), expected) ||
      any(!purrr::map_lgl(results, ~ is.list(.x$correction)))
  ) {
    stop("The observation exact corrections are incomplete.",
         call. = FALSE)
  }

  purrr::map(results, "correction")
}

# Replace only the matching pointwise scores and retain the observation
# manifest for downstream comparison.
apply_rescor_observation_exact_corrections <- function(
    observation_loo_manifest,
    exact_results
) {
  initial_results <- purrr::map(
    observation_loo_manifest$integrated_loo_path,
    qs::qread
  )
  initial_keys <- purrr::map_chr(
    initial_results, ~ paste(.x$outcome, .x$model_index)
  )

  if (
    length(initial_results) != 104L || anyDuplicated(initial_keys) ||
      any(!purrr::map_lgl(
        initial_results,
        ~ identical(.x$loo_unit, "observation_conditional_on_study")
      ))
  ) {
    stop("The initial observation-level LOO set is invalid.",
         call. = FALSE)
  }

  warning_cases <- purrr::map_dfr(initial_results, function(result) {
    pareto_k <- as.numeric(loo::pareto_k_values(result$loo))
    tibble::tibble(
      outcome = result$outcome,
      model_index = result$model_index,
      observation = which(pareto_k > result$k_threshold)
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
    stop("Exact observation corrections do not match the high-k cases.",
         call. = FALSE)
  }
  final_results <- initial_results
  keys <- unique(paste(exact_cases$outcome, exact_cases$model_index))

  for (key in keys) {
    position <- match(key, initial_keys)
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
    patched$remedy <-
      "correlated_conditional_observation_level_psis_plus_exact"
    final_results[[position]] <- patched
  }

  exact_summary <- purrr::map_dfr(exact_results, "exact_refit_summary")
  model_summary <- purrr::map_dfr(final_results, function(result) {
    summarise_primary_loo_method(
      result$loo,
      "correlated_conditional_observation_level_psis_plus_exact"
    ) |>
      dplyr::mutate(
        outcome = result$outcome,
        model_index = result$model_index,
        model_name = result$model_name,
        .before = 1L
      )
  })

  list(
    family = "rescor",
    results = final_results,
    model_summary = model_summary,
    exact_refit_summary = exact_summary,
    stacking_ready =
      all(model_summary$n_pareto_k_gt_threshold == 0L) &&
      all(exact_summary$diagnostic_pass),
    loo_unit = "observation_conditional_on_study",
    study_level_audit_retained = TRUE
  )
}

# Exact study holdouts and reporting ----

# Hold out the entire study cluster for the separate study-level check. Later
# reporting reads the accepted fit and checkpoint to generate prediction
# draws.
make_rescor_leave_one_study_out_data <- function(stan_data, cluster) {
  cluster <- as.integer(cluster)
  if (
    length(cluster) != 1L || is.na(cluster) || cluster < 1L ||
      cluster > stan_data$N_cluster
  ) {
    stop("Invalid held-out residual-correlation study.", call. = FALSE)
  }

  first <- stan_data$cluster_start[[cluster]]
  last <- stan_data$cluster_end[[cluster]]
  keep <- setdiff(seq_len(stan_data$N), seq.int(first, last))
  row_fields <- c(
    "y", "se_y", "x", "study", "arm", "outcome",
    "resistance_exercise_code", "age_z", "duration_z"
  )

  if (
    !all(row_fields %in% names(stan_data)) ||
      any(vapply(stan_data[row_fields], length, integer(1)) != stan_data$N)
  ) {
    stop("Residual-correlation row data are incomplete.", call. = FALSE)
  }

  for (field in row_fields) {
    stan_data[[field]] <- stan_data[[field]][keep]
  }
  remap <- function(values) match(values, unique(values))
  stan_data$study <- remap(stan_data$study)
  stan_data$arm <- remap(stan_data$arm)
  stan_data$outcome <- remap(stan_data$outcome)
  stan_data$N <- length(keep)
  stan_data$N_study <- max(stan_data$study)
  stan_data$N_arm <- max(stan_data$arm)
  stan_data$N_outcome <- max(stan_data$outcome)
  runs <- rle(stan_data$study)
  stan_data$N_cluster <- length(runs$lengths)
  stan_data$cluster_end <- as.integer(cumsum(runs$lengths))
  stan_data$cluster_start <- as.integer(c(
    1L, utils::head(stan_data$cluster_end, -1L) + 1L
  ))

  if (max(runs$lengths) <= 1L) {
    stop("The training data cannot identify residual correlation.",
         call. = FALSE)
  }
  stan_data$rho_lower <- -1 / (max(runs$lengths) - 1) + 1e-8
  stan_data$rho_raw_prior_location <- stats::qlogis(
    -stan_data$rho_lower / (1 - stan_data$rho_lower)
  )
  stan_data
}

rescor_full_exact_fit_file_path <- function(record) {
  directory <- here::here(
    "_targets", "user", "sensitivity_fit_files", "rescor_exact",
    as.character(record$outcome),
    paste0("model_", stringr::str_pad(record$model_index, 2L, pad = "0"))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, paste0("case_", record$case_id, ".qs"))
}

read_rescor_full_exact_loo_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("The residual-correlation exact-refit file is unavailable.",
         call. = FALSE)
  }
  refit <- qs::qread(path)
  output_files <- tryCatch(
    refit$fitted_candidate$fit$output_files(),
    error = function(e) character()
  )

  if (
    !is.list(refit) || !identical(refit$family, "rescor") ||
      !inherits(refit$fitted_candidate$fit, "CmdStanMCMC") ||
      length(output_files) < 1L || any(!file.exists(output_files))
  ) {
    stop("The residual-correlation exact-refit file is invalid.",
         call. = FALSE)
  }

  refit
}

fit_rescor_full_exact_loo_to_file <- function(
    record,
    model,
    enabled,
    config
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip("rescor"))

  if (is_sensitivity_pilot_skip(model)) {
    stop("The residual-correlation Stan model was not compiled.",
         call. = FALSE)
  }
  existing_path <- rescor_full_exact_fit_file_path(record)
  existing <- if (file.exists(existing_path)) {
    tryCatch(
      read_rescor_full_exact_loo_file(existing_path),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  if (
    !is.null(existing) &&
      identical(as.integer(existing$case_id), as.integer(record$case_id)) &&
      identical(existing$outcome, record$outcome) &&
      identical(
        as.integer(existing$model_index), as.integer(record$model_index)
      ) &&
      identical(
        as.integer(existing$observation), as.integer(record$observation)
      )
  ) {
    return(normalizePath(existing_path, winslash = "/", mustWork = TRUE))
  }

  integrated_result <- qs::qread(record$integrated_loo_path)
  input <- integrated_result$input
  if (
    !identical(integrated_result$family, "rescor") ||
      input$specification$model_index != record$model_index
  ) {
    stop("The residual-correlation exact-refit record is invalid.",
         call. = FALSE)
  }

  training_data <- make_rescor_leave_one_study_out_data(
    input$stan_data, record$observation
  )
  specification <- input$specification
  specification$adapt_delta <- max(
    specification$adapt_delta, config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth, config$max_treedepth_floor
  )
  output_dir <- make_rescor_full_cmdstan_output_dir(
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
    family = "rescor",
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
    family = "rescor",
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

rescor_held_out_integrated_log_lik <- function(
    candidate,
    input,
    cluster,
    loo_draws = 4000L
) {
  data <- input$stan_data
  cluster <- as.integer(cluster)
  first <- data$cluster_start[[cluster]]
  last <- data$cluster_end[[cluster]]
  rows <- seq.int(first, last)
  variables <- c(
    "b_intercept", "b_rt", "b_protein", "b_protein_rt",
    "b_moderator", "sd_study", "L_study", "sd_arm",
    "sd_outcome", "rho"
  )
  all_draws <- posterior::as_draws_df(candidate$fit$draws(variables))
  draw_ids <- make_sensitivity_loo_draw_ids(
    nrow(all_draws), min(as.integer(loo_draws), nrow(all_draws))
  )
  draws <- as.data.frame(all_draws)[draw_ids, , drop = FALSE]
  n_draws <- nrow(draws)
  k_basis <- as.integer(data$K_basis)
  k_study <- as.integer(data$K_study)
  basis_list <- evaluate_eiv_stan_basis(
    matrix(data$x[rows], nrow = 1L), data
  )
  basis <- do.call(cbind, lapply(basis_list, as.numeric))
  fixed_mean <- outer(
    as.numeric(draws$b_intercept), rep(1, length(rows))
  ) + outer(
    as.numeric(draws$b_rt), data$resistance_exercise_code[rows]
  )

  if (data$include_protein == 1L) {
    coefficient <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_protein", k_basis
    )
    fixed_mean <- fixed_mean + coefficient %*% t(basis)
  }

  if (data$include_rt_interaction == 1L) {
    coefficient <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_protein_rt", k_basis
    )
    fixed_mean <- fixed_mean +
      (coefficient %*% t(basis)) *
      matrix(
        data$resistance_exercise_code[rows],
        nrow = n_draws,
        ncol = length(rows),
        byrow = TRUE
      )
  }

  if (data$include_moderators == 1L) {
    moderators <- extract_cmdstan_indexed_draw_matrix(
      draws, "b_moderator", 2L
    )
    fixed_mean <- fixed_mean +
      outer(moderators[, 1L], data$duration_z[rows]) +
      outer(moderators[, 2L], data$age_z[rows])
  }

  study_design <- if (data$include_random_slope == 1L) {
    cbind(1, basis)
  } else {
    matrix(1, nrow = length(rows), ncol = 1L)
  }
  study_sd <- extract_cmdstan_indexed_draw_matrix(
    draws, "sd_study", k_study
  )
  study_cholesky <- if (data$correlated_study_effects == 1L) {
    extract_cmdstan_matrix_draw_array(
      draws, "L_study", k_study, k_study
    )
  } else {
    NULL
  }
  arm_match <- outer(data$arm[rows], data$arm[rows], "==")
  outcome_match <- outer(
    data$outcome[rows], data$outcome[rows], "=="
  )
  residual_scale <- outer(data$se_y[rows], data$se_y[rows])
  observed <- data$y[rows]
  output <- numeric(n_draws)

  for (draw in seq_len(n_draws)) {
    study_covariance <- if (data$correlated_study_effects == 1L) {
      scaled_cholesky <- diag(study_sd[draw, ], nrow = k_study) %*%
        study_cholesky[draw, , ]
      tcrossprod(scaled_cholesky)
    } else {
      diag(study_sd[draw, ]^2, nrow = k_study)
    }
    residual_covariance <- residual_scale * as.numeric(draws$rho[[draw]])
    diag(residual_covariance) <- data$se_y[rows]^2
    covariance <- residual_covariance +
      study_design %*% study_covariance %*% t(study_design) +
      arm_match * as.numeric(draws$sd_arm[[draw]])^2 +
      outcome_match * as.numeric(draws$sd_outcome[[draw]])^2
    cholesky <- chol(covariance)
    centred <- observed - fixed_mean[draw, ]
    standardised <- forwardsolve(t(cholesky), centred)
    output[[draw]] <- -0.5 * (
      length(rows) * log(2 * pi) +
        2 * sum(log(diag(cholesky))) +
        sum(standardised^2)
    )
  }

  output
}

diagnose_rescor_full_exact_loo_file <- function(path) {
  refit <- read_rescor_full_exact_loo_file(path)
  summarise_eiv_exact_loo_refit(refit)
}

compute_rescor_full_exact_loo_file <- function(path, diagnostic_audit) {
  if (!isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)) {
    stop("Residual-correlation exact-refit diagnostics have not passed.",
         call. = FALSE)
  }
  refit <- read_rescor_full_exact_loo_file(path)
  held_out_log_lik <- rescor_held_out_integrated_log_lik(
    refit$fitted_candidate,
    refit$input,
    refit$observation,
    loo_draws = 4000L
  )
  maximum_log_lik <- max(held_out_log_lik)
  scaled_likelihood <- exp(held_out_log_lik - maximum_log_lik)
  exact_elpd <- maximum_log_lik + log(mean(scaled_likelihood))
  exact_mcse <- tryCatch(
    posterior::mcse_mean(scaled_likelihood) / mean(scaled_likelihood),
    error = function(e) NA_real_
  )
  patched_loo <- refit$integrated_result$loo
  observation <- as.integer(refit$observation)
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
  cluster_row <- refit$input$cluster_manifest[observation, , drop = FALSE]
  exact_refit_summary <- summarise_eiv_exact_loo_refit(refit) |>
    dplyr::mutate(
      held_out_study = cluster_row$study_no[[1L]],
      held_out_cluster_size = cluster_row$cluster_size[[1L]],
      old_psis_elpd_loo = old_elpd,
      exact_elpd_loo = exact_elpd,
      exact_minus_psis_elpd = exact_elpd - old_elpd,
      exact_elpd_mcse = exact_mcse
    )

  list(
    case_id = refit$case_id,
    family = "rescor",
    outcome = refit$outcome,
    model_index = refit$model_index,
    model_name = refit$model_name,
    observation_manifest = refit$input$cluster_manifest,
    loo = patched_loo,
    conditional_loo = refit$integrated_result$conditional_loo,
    remedy = "study_cluster_integrated_psis_plus_exact",
    k_threshold = refit$integrated_result$k_threshold,
    exact_refit_summary = exact_refit_summary
  )
}

apply_rescor_full_exact_loo_corrections <- function(
    integrated_manifest,
    exact_results
) {
  audit <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest,
    exact_results,
    family = "rescor"
  )
  audit$model_summary$method <-
    "integrated_study_cluster_psis_plus_exact"
  audit$results <- purrr::map(audit$results, function(result) {
    result$remedy <- if (anyNA(loo::pareto_k_values(result$loo))) {
      "study_cluster_integrated_psis_plus_exact"
    } else {
      "study_cluster_integrated_psis"
    }
    result
  })

  audit
}

make_rescor_full_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    n_draws = 4000L,
    seed = 39988L
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
    record$family <- "rescor"
    record$prediction_type <- "rescor"
    record
  })
}

write_rescor_full_prediction_draw_file <- function(
    record,
    checkpoint_path = NULL,
    result_path = NULL
) {
  if (
    !identical(record$family, "rescor") || record$draw_count < 1L ||
      !file.exists(record$fit_path) ||
      !identical(
        digest::digest(record$fit_path, file = TRUE, algo = "sha256"),
        record$fit_sha256
      )
  ) {
    stop("A residual-correlation prediction record is invalid.",
         call. = FALSE)
  }

  checkpoint <- if (is.null(checkpoint_path)) {
    NULL
  } else {
    read_sensitivity_completed_checkpoint(
      checkpoint_path, family = "rescor"
    )
  }
  candidate <- read_rescor_fit_file(record$fit_path)

  if (
    candidate$model_index != record$model_index ||
      !identical(candidate$outcome, record$outcome) ||
      !identical(candidate$model_name, record$model_name) ||
      (!is.null(checkpoint) &&
        is.null(checkpoint$analysis_data_by_outcome[[record$outcome]]))
  ) {
    stop("A residual-correlation prediction fit does not match its record.",
         call. = FALSE)
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
  grid <- dplyr::transmute(
    candidate$grid,
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
    stop("Residual-correlation prediction draws are incomplete.",
         call. = FALSE)
  }

  reduction <- list(
    family = "rescor",
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
      prediction_type = "rescor",
      reference_intake = NULL
    )
  )

  path <- if (is.null(result_path)) {
    sensitivity_reduction_file_path(
      family = "rescor",
      outcome = record$outcome,
      model_index = record$model_index,
      stage = "prediction_draws"
    )
  } else {
    dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
    result_path
  }

  write_atomic_qs(reduction, path)
}
