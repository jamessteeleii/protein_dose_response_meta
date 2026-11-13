# Combined-data sensitivity-analysis updates ----

make_update_2026_mundlak_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated Mundlak base specification registry is incomplete.",
         call. = FALSE)
  }

  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  remedy_keys <- unname(purrr::map_chr(primary_remedy_specifications, key))
  if (
    anyDuplicated(primary_keys) || anyDuplicated(remedy_keys) ||
      !identical(remedy_keys, "function::10") ||
      !all(remedy_keys %in% primary_keys)
  ) {
    stop("The updated Mundlak primary remedy replacement is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  specifications <- purrr::map(specifications, function(specification) {
    specification$update_2026_sensitivity_base <- TRUE
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  function_ten <- specifications[[match("function::10", final_keys)]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !identical(as.numeric(function_ten$adapt_delta), 0.9995)
  ) {
    stop("The updated Mundlak specification registry failed validation.",
         call. = FALSE)
  }

  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

validate_update_2026_mundlak_inputs <- function(
    specifications,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  checks <- purrr::map_dfr(specifications, function(specification) {
    input <- prepare_mundlak_full_input(
      specification = specification,
      data_by_outcome = data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      config = config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      family = input$family,
      rows = nrow(input$data),
      basis_dimension = if (is.null(input$basis_spec)) {
        0L
      } else {
        as.integer(input$basis_spec$basis_dimension)
      },
      regularised_variance_priors =
        isTRUE(input$specification$regularised_variance_priors),
      adapt_delta = max(
        input$specification$adapt_delta,
        config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        config$max_treedepth_floor
      ),
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  keys <- paste(checks$outcome, checks$model_index)
  counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  model_four <- checks |>
    dplyr::filter(.data$model_index == 4L)
  function_ten <- checks |>
    dplyr::filter(.data$outcome == "function", .data$model_index == 10L)
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      nrow(model_four) != 4L ||
      !all(model_four$regularised_variance_priors) ||
      nrow(function_ten) != 1L ||
      function_ten$adapt_delta[[1L]] != 0.9995 ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("The updated Mundlak input preflight is incomplete.", call. = FALSE)
  }
  checks
}

write_update_2026_mundlak_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "mundlak",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated Mundlak preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_mundlak_candidate_to_file <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "mundlak",
      "candidates"
    )
) {
  input <- prepare_mundlak_full_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = config
  )
  candidate <- fit_mundlak_full_candidate(
    input = input,
    enabled = TRUE,
    config = config
  )

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated Mundlak fit file.", call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)

  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated Mundlak fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

compile_update_2026_mundlak_diagnostic_audit <- function(diagnostics) {
  split_diagnostics <- split(
    diagnostics,
    purrr::map_chr(diagnostics, ~ .x$summary$outcome[[1L]])
  )
  required <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !setequal(names(split_diagnostics), required) ||
      any(vapply(split_diagnostics[required], length, integer(1L)) != 26L)
  ) {
    stop("Updated Mundlak diagnostics are incomplete.", call. = FALSE)
  }

  compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
}

select_update_2026_mundlak_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  expected_keys <- c(
    "whole_body::3", "whole_body::4", "whole_body::6",
    "whole_body::13", "strength::5", "strength::6",
    "strength::25", "function::5"
  )
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(
      outcome,
      model_index = as.integer(model_index),
      key = paste(outcome, model_index, sep = "::")
    )
  strength_25 <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(.data$outcome == "strength", .data$model_index == 25L)
  if (
    nrow(failed) != 8L || anyDuplicated(failed$key) ||
      !setequal(failed$key, expected_keys) ||
      any(diagnostic_audit$diagnostic_summary$divergent_transitions > 0L) ||
      nrow(strength_25) != 1L ||
      strength_25$maximum_treedepth_transitions[[1L]] != 11L
  ) {
    stop("The updated Mundlak remedy set no longer matches its audit.",
         call. = FALSE)
  }

  selected <- purrr::keep(specifications, function(specification) {
    paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    ) %in% failed$key
  })
  selected_keys <- unname(purrr::map_chr(
    selected,
    ~ paste(.x$outcome, .x$model_index, sep = "::")
  ))
  if (
    length(selected) != 8L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The updated Mundlak remedy specifications are incomplete.",
         call. = FALSE)
  }
  selected
}

validate_update_2026_mundlak_remedy_inputs <- function(
    specifications,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  checks <- purrr::map_dfr(specifications, function(specification) {
    input <- prepare_mundlak_full_remedy_input(
      specification = specification,
      data_by_outcome = data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      config = config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      remedy = input$specification$sensitivity_mundlak_full_remedy,
      adapt_delta = max(
        input$specification$adapt_delta,
        config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        config$max_treedepth_floor
      ),
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  expected_keys <- c(
    "whole_body::3", "whole_body::4", "whole_body::6",
    "whole_body::13", "strength::5", "strength::6",
    "strength::25", "function::5"
  )
  strength_25 <- checks |>
    dplyr::filter(.data$outcome == "strength", .data$model_index == 25L)
  other_cases <- checks |>
    dplyr::filter(!(.data$outcome == "strength" & .data$model_index == 25L))
  if (
    nrow(checks) != 8L || anyDuplicated(keys) ||
      !setequal(keys, expected_keys) || !all(checks$valid) ||
      nrow(strength_25) != 1L || strength_25$max_treedepth[[1L]] != 15L ||
      strength_25$remedy[[1L]] !=
        "longer_sampling_and_max_treedepth_15" ||
      any(other_cases$remedy != "longer_sampling") ||
      config$iter_warmup != 2000L || config$iter_sampling != 4000L
  ) {
    stop("The updated Mundlak remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_mundlak_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "remedy_1",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 8L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated Mundlak remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_mundlak_remedy_to_file <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "mundlak",
      "remedies", "remedy_1"
    )
) {
  input <- prepare_mundlak_full_remedy_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = config
  )
  candidate <- fit_mundlak_full_candidate(
    input = input,
    enabled = TRUE,
    config = config
  )

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated Mundlak remedy fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated Mundlak remedy fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

compile_update_2026_mundlak_remedy_diagnostic_audit <- function(
    initial_diagnostics,
    remedy_diagnostics
) {
  diagnostic_key <- function(x) {
    paste(x$summary$outcome[[1L]], x$summary$model_index[[1L]], sep = "::")
  }
  if (length(initial_diagnostics) != 104L || length(remedy_diagnostics) != 8L) {
    stop("Updated Mundlak remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  initial_keys <- unname(purrr::map_chr(initial_diagnostics, diagnostic_key))
  remedy_keys <- unname(purrr::map_chr(remedy_diagnostics, diagnostic_key))
  expected_keys <- initial_diagnostics |>
    purrr::keep(~ !.x$summary$diagnostic_pass[[1L]]) |>
    purrr::map_chr(diagnostic_key)
  if (
    anyDuplicated(initial_keys) || anyDuplicated(remedy_keys) ||
      !setequal(remedy_keys, expected_keys) ||
      !all(remedy_keys %in% initial_keys)
  ) {
    stop("Updated Mundlak diagnostic replacement is invalid.",
         call. = FALSE)
  }

  final_diagnostics <- initial_diagnostics
  final_diagnostics[match(remedy_keys, initial_keys)] <- remedy_diagnostics
  split_diagnostics <- split(
    final_diagnostics,
    purrr::map_chr(final_diagnostics, ~ .x$summary$outcome[[1L]])
  )
  compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
}

select_update_2026_mundlak_second_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass)
  if (
    nrow(failed) != 1L || failed$outcome[[1L]] != "strength" ||
      failed$model_index[[1L]] != 6L ||
      failed$divergent_transitions[[1L]] != 0L ||
      failed$maximum_treedepth_transitions[[1L]] != 3L ||
      failed$max_treedepth[[1L]] != 12L ||
      isTRUE(failed$geometry_pass[[1L]]) ||
      !isTRUE(failed$convergence_pass[[1L]])
  ) {
    stop("The updated Mundlak second remedy no longer matches its audit.",
         call. = FALSE)
  }

  selected <- purrr::keep(
    specifications,
    ~ identical(.x$outcome, "strength") && .x$model_index == 6L
  )
  if (length(selected) != 1L) {
    stop("Could not uniquely select Mundlak strength model 6.",
         call. = FALSE)
  }
  selected
}

validate_update_2026_mundlak_second_remedy_inputs <- function(
    specifications,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  if (length(specifications) != 1L) {
    stop("Expected one updated Mundlak second-remedy specification.",
         call. = FALSE)
  }
  specification <- specifications[[1L]]
  input <- prepare_mundlak_full_remedy_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = config
  )
  brms::make_standata(
    formula = stats::as.formula(input$formula),
    data = input$data,
    prior = input$priors
  )
  checks <- tibble::tibble(
    outcome = specification$outcome,
    model_index = as.integer(specification$model_index),
    model_name = specification$model_name,
    remedy = "max_treedepth_15_after_first_remedy",
    adapt_delta = max(
      input$specification$adapt_delta,
      config$adapt_delta_floor
    ),
    max_treedepth = max(
      input$specification$max_treedepth,
      config$max_treedepth_floor
    ),
    iter_warmup = as.integer(config$iter_warmup),
    iter_sampling = as.integer(config$iter_sampling),
    valid = TRUE
  )
  if (
    nrow(checks) != 1L || checks$outcome[[1L]] != "strength" ||
      checks$model_index[[1L]] != 6L ||
      checks$max_treedepth[[1L]] != 15L ||
      checks$iter_warmup[[1L]] != 2000L ||
      checks$iter_sampling[[1L]] != 4000L || !checks$valid[[1L]]
  ) {
    stop("The updated Mundlak second-remedy preflight is invalid.",
         call. = FALSE)
  }
  checks
}

write_update_2026_mundlak_second_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "remedy_2",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 1L || !checks$valid[[1L]]) {
    stop("Cannot write an invalid updated Mundlak second-remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

compile_update_2026_mundlak_second_remedy_diagnostic_audit <- function(
    initial_diagnostics,
    first_remedy_diagnostics,
    second_remedy_diagnostics
) {
  diagnostic_key <- function(x) {
    paste(x$summary$outcome[[1L]], x$summary$model_index[[1L]], sep = "::")
  }
  if (
    length(initial_diagnostics) != 104L ||
      length(first_remedy_diagnostics) != 8L ||
      length(second_remedy_diagnostics) != 1L
  ) {
    stop("Updated Mundlak second-remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  initial_keys <- unname(purrr::map_chr(initial_diagnostics, diagnostic_key))
  first_keys <- unname(purrr::map_chr(first_remedy_diagnostics, diagnostic_key))
  second_keys <- unname(purrr::map_chr(second_remedy_diagnostics, diagnostic_key))
  expected_first_keys <- initial_diagnostics |>
    purrr::keep(~ !.x$summary$diagnostic_pass[[1L]]) |>
    purrr::map_chr(diagnostic_key)
  remaining_first_failures <- first_remedy_diagnostics |>
    purrr::keep(~ !.x$summary$diagnostic_pass[[1L]]) |>
    purrr::map_chr(diagnostic_key)
  if (
    anyDuplicated(initial_keys) || anyDuplicated(first_keys) ||
      anyDuplicated(second_keys) ||
      !setequal(first_keys, expected_first_keys) ||
      !identical(unname(remaining_first_failures), "strength::6") ||
      !identical(second_keys, "strength::6")
  ) {
    stop("Updated Mundlak second diagnostic replacement is invalid.",
         call. = FALSE)
  }

  final_diagnostics <- initial_diagnostics
  final_diagnostics[match(first_keys, initial_keys)] <- first_remedy_diagnostics
  final_diagnostics[match(second_keys, initial_keys)] <- second_remedy_diagnostics
  split_diagnostics <- split(
    final_diagnostics,
    purrr::map_chr(final_diagnostics, ~ .x$summary$outcome[[1L]])
  )
  compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
}

make_update_2026_mundlak_accepted_fit_manifest <- function(
    initial_fit_paths,
    first_remedy_fit_paths,
    second_remedy_fit_paths,
    final_diagnostic_audit
) {
  if (
    nrow(final_diagnostic_audit$diagnostic_summary) != 104L ||
      any(!final_diagnostic_audit$diagnostic_summary$diagnostic_pass)
  ) {
    stop("Updated Mundlak fits cannot be accepted before diagnostics pass.",
         call. = FALSE)
  }
  path_table <- function(paths, stage, expected) {
    paths <- unname(unlist(paths, use.names = FALSE))
    table <- tibble::tibble(
      outcome = basename(dirname(paths)),
      model_index = as.integer(stringr::str_match(
        basename(paths), "model_([0-9]+)\\.qs$"
      )[, 2]),
      stage = stage,
      fit_path = normalizePath(paths, winslash = "/", mustWork = TRUE)
    )
    if (
      nrow(table) != expected || anyNA(table$model_index) ||
        anyDuplicated(table[c("outcome", "model_index")]) ||
        any(file.info(table$fit_path)$size <= 0)
    ) {
      stop("An updated Mundlak fit stage is incomplete.", call. = FALSE)
    }
    table
  }

  initial <- path_table(initial_fit_paths, "initial", 104L)
  first <- path_table(first_remedy_fit_paths, "remedy_1", 8L)
  second <- path_table(second_remedy_fit_paths, "remedy_2", 1L)
  keys <- c("outcome", "model_index")
  first_keys <- paste(first$outcome, first$model_index, sep = "::")
  expected_first_keys <- c(
    "whole_body::3", "whole_body::4", "whole_body::6",
    "whole_body::13", "strength::5", "strength::6",
    "strength::25", "function::5"
  )
  second_key <- paste(second$outcome, second$model_index, sep = "::")
  if (
    !setequal(first_keys, expected_first_keys) ||
      !identical(second_key, "strength::6")
  ) {
    stop("Updated Mundlak accepted replacements are invalid.",
         call. = FALSE)
  }

  accepted <- initial |>
    dplyr::anti_join(first[keys], by = keys) |>
    dplyr::bind_rows(
      first |>
        dplyr::anti_join(second[keys], by = keys),
      second
    ) |>
    dplyr::left_join(
      final_diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, model_name = model),
      by = keys,
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      family = "mundlak",
      bytes = file.info(fit_path)$size,
      sha256 = purrr::map_chr(
        fit_path,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      )
    ) |>
    dplyr::select(
      family, outcome, model_index, model_name, stage,
      fit_path, bytes, sha256
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  stage_counts <- table(factor(
    accepted$stage,
    levels = c("initial", "remedy_1", "remedy_2")
  ))
  if (
    nrow(accepted) != 104L ||
      anyDuplicated(accepted[c("outcome", "model_index")]) ||
      anyNA(accepted$model_name) || any(!file.exists(accepted$fit_path)) ||
      !identical(as.integer(stage_counts), c(96L, 7L, 1L)) ||
      any(!nzchar(accepted$sha256))
  ) {
    stop("The updated Mundlak accepted-fit manifest is invalid.",
         call. = FALSE)
  }
  accepted
}

write_update_2026_mundlak_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "mundlak",
      "accepted_fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(manifest[c("outcome", "model_index")]) ||
      any(!file.exists(manifest$fit_path))
  ) {
    stop("Cannot write an invalid updated Mundlak fit manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

compute_update_2026_mundlak_integrated_loo_to_file <- function(
    fit_record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "mundlak", "loo",
      "integrated"
    )
) {
  result <- compute_mundlak_full_file_integrated_loo(
    fit_record = fit_record,
    loo_draws = loo_draws,
    k_threshold = k_threshold,
    cores = cores,
    family = "mundlak",
    slope_description = "mundlak_within_basis_slopes"
  )
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(unname(as.character(result$family)), "mundlak") ||
      !identical(result$remedy, "integrated_observation_level_psis")
  ) {
    stop("Cannot archive an invalid updated Mundlak LOO result.",
         call. = FALSE)
  }
  output_dir <- file.path(output_root, result$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf(
      "model_%02d_integrated_loo.qs",
      as.integer(result$model_index)
    )
  )
  write_atomic_qs(result, output_file)
}

write_update_2026_mundlak_integrated_loo_audit <- function(
    audit,
    manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "loo_psis"
    )
) {
  if (
    !isTRUE(audit$integrated_psis_complete) || audit$models != 104L ||
      nrow(manifest) != 104L ||
      any(!file.exists(manifest$integrated_loo_path))
  ) {
    stop("The updated Mundlak integrated-LOO audit is incomplete.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

validate_update_2026_mundlak_exact_inputs <- function(records, config) {
  if (!is.list(records) || length(records) != 15L) {
    stop("Expected 15 updated Mundlak exact-LOO records.", call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    integrated <- qs::qread(record$integrated_loo_path)
    observation <- as.integer(record$observation)
    if (
      !inherits(integrated$loo, "loo") ||
        !identical(unname(as.character(integrated$family)), "mundlak") ||
        !identical(integrated$outcome, record$outcome) ||
        as.integer(integrated$model_index) != record$model_index ||
        observation < 1L ||
        observation > nrow(integrated$observation_manifest)
    ) {
      stop("An updated Mundlak exact record is invalid.", call. = FALSE)
    }
    held_out <- integrated$observation_manifest[observation, , drop = FALSE]
    observed_k <- as.numeric(
      loo::pareto_k_values(integrated$loo)[[observation]]
    )
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      observation = observation,
      study_no = as.character(held_out$study_no[[1L]]),
      arm_no = as.integer(held_out$arm_no[[1L]]),
      outcome_no = as.integer(held_out$outcome_no[[1L]]),
      pareto_k = observed_k,
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      original_fit_exists = file.exists(record$original_fit_path),
      integrated_loo_exists = file.exists(record$integrated_loo_path),
      valid = is.finite(observed_k) && observed_k > 0.7
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- c(
    "whole_body::4::4", "whole_body::8::3",
    "whole_body::16::3", "whole_body::17::3",
    "whole_body::19::3", "whole_body::20::3",
    "whole_body::20::4", "whole_body::23::4",
    "whole_body::26::3", "hypertrophy::8::252",
    "hypertrophy::8::253", "hypertrophy::8::263",
    "hypertrophy::18::253", "hypertrophy::24::253",
    "hypertrophy::25::252"
  )
  if (
    nrow(checks) != 15L || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(15L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$integrated_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L
  ) {
    stop("The updated Mundlak exact-LOO preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_mundlak_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "loo_exact",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 15L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated Mundlak exact preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_mundlak_exact_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "mundlak", "loo",
    "exact", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%02d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

fit_and_score_update_2026_mundlak_exact_to_file <- function(record, config) {
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !identical(unname(as.character(candidate$family)), "mundlak") ||
      !identical(candidate$outcome, record$outcome) ||
      as.integer(candidate$model_index) != record$model_index ||
      !inherits(integrated$loo, "loo") ||
      !identical(unname(as.character(integrated$family)), "mundlak") ||
      !identical(integrated$outcome, record$outcome) ||
      as.integer(integrated$model_index) != record$model_index ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      loo::pareto_k_values(integrated$loo)[[observation]] <= 0.7
  ) {
    stop("An updated Mundlak exact-LOO record failed validation.",
         call. = FALSE)
  }

  original_data <- candidate$model$data
  training_data <- original_data[-observation, , drop = FALSE]
  original_model <- candidate$model
  sample_prior <- attr(original_model$prior, "sample_prior")
  original_model$prior <- normalize_primary_brms_priors(
    original_model$prior
  )
  attr(original_model$prior, "sample_prior") <- sample_prior
  specification <- candidate$specification
  specification$adapt_delta <- max(
    specification$adapt_delta,
    config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth,
    config$max_treedepth_floor
  )

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
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  refit <- list(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    original_data = original_data,
    integrated_result = integrated,
    specification = specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  correction <- compute_sensitivity_mundlak_exact_loo_case(refit)
  if (
    as.integer(correction$case_id) != record$case_id ||
      !inherits(correction$loo, "loo") ||
      nrow(correction$exact_refit_summary) != 1L
  ) {
    stop("The updated Mundlak exact correction is invalid.",
         call. = FALSE)
  }
  write_atomic_qs(
    correction,
    update_2026_mundlak_exact_file_path(record)
  )
}

audit_update_2026_mundlak_exact_files <- function(paths, integrated_audit) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    length(paths) != expected || expected != 15L || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0)
  ) {
    stop("The updated Mundlak exact-correction files are incomplete.",
         call. = FALSE)
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(outcome, levels = c("whole_body", "hypertrophy")),
      model_index,
      observation
    )
  manifest <- purrr::map2_dfr(paths, corrections, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      observation =
        as.integer(result$exact_refit_summary$observation[[1L]]),
      correction_path = normalizePath(
        path, winslash = "/", mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path, algo = "sha256", serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys)
  ) {
    stop("The updated Mundlak exact audit is invalid.", call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = sum(diagnostic_summary$diagnostic_pass),
    exact_refits_requiring_remediation =
      sum(!diagnostic_summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass =
      all(diagnostic_summary$diagnostic_pass)
  )
}

write_update_2026_mundlak_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "loo_exact"
    )
) {
  if (
    audit$exact_refits != 15L ||
      nrow(audit$diagnostic_summary) != 15L ||
      nrow(audit$correction_manifest) != 15L
  ) {
    stop("Cannot write an incomplete updated Mundlak exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

read_update_2026_mundlak_exact_corrections <- function(manifest) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 15L ||
      any(!file.exists(manifest$correction_path))
  ) {
    stop("The updated Mundlak exact-correction manifest is invalid.",
         call. = FALSE)
  }
  corrections <- purrr::map(manifest$correction_path, qs::qread)
  corrections <- corrections[order(purrr::map_int(corrections, "case_id"))]
  if (
    !identical(
      unname(purrr::map_int(corrections, "case_id")),
      seq_len(15L)
    ) ||
      any(!purrr::map_lgl(corrections, ~ inherits(.x$loo, "loo")))
  ) {
    stop("The updated Mundlak exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections
}

make_update_2026_mundlak_final_loo_audit <- function(
    integrated_manifest,
    exact_corrections,
    exact_audit
) {
  if (!isTRUE(exact_audit$all_exact_refit_diagnostics_pass)) {
    stop("Final Mundlak LOO requires passing exact-refit diagnostics.",
         call. = FALSE)
  }
  final <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest = integrated_manifest,
    exact_results = exact_corrections,
    family = "mundlak"
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- integrated_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "integrated_observation_level_psis_plus_exact",
        "integrated_observation_level_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 15L
  ) {
    stop("The final updated Mundlak LOO audit is not stacking-ready.",
         call. = FALSE)
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 15L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_mundlak_final_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 15L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop("Cannot write an incomplete final Mundlak LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_mundlak_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "mundlak",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L || anyDuplicated(
      paste(weights$outcome, weights$model_index)
    ) ||
      length(outcome_sums) != 4L || any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop("The updated Mundlak stacking weights are invalid.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_mundlak_reporting_preflight <- function(records) {
  if (!is.list(records) || length(records) < 4L) {
    stop("The updated Mundlak prediction allocation is incomplete.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    tibble::tibble(
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      grid_step = as.numeric(record$grid_step),
      fit_path = record$fit_path,
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path)
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop("The updated Mundlak reporting preflight failed.", call. = FALSE)
  }
  checks
}

write_update_2026_mundlak_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "reporting",
      "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(
    checks$allocated_draws,
    checks$outcome,
    sum
  )
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      !all(checks$valid)
  ) {
    stop("Cannot write an invalid updated Mundlak reporting preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_mundlak_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "mundlak",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_mundlak_prediction_draw_file <- function(
    record,
    data_by_outcome
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L
  ) {
    stop("An updated Mundlak prediction record is invalid.",
         call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "mundlak") ||
      !identical(candidate$outcome, unname(as.character(record$outcome))) ||
      as.integer(candidate$model_index) != as.integer(record$model_index) ||
      !identical(candidate$model_name, unname(as.character(record$model_name)))
  ) {
    stop("An updated Mundlak prediction fit does not match its record.",
         call. = FALSE)
  }
  prediction_data <- data_by_outcome[[record$outcome]]
  if (
    is.null(prediction_data) || !is.data.frame(prediction_data) ||
      !all(c("protein_intake_raw", "protein_median") %in%
        names(prediction_data))
  ) {
    stop("The updated Mundlak prediction data are unavailable.",
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
    stop("Updated Mundlak prediction draws are incomplete.",
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
  write_atomic_qs(
    reduction,
    update_2026_mundlak_reporting_path(
      outcome = record$outcome,
      stage = "prediction_draws",
      model_index = record$model_index
    )
  )
}

write_update_2026_mundlak_model_average_file <- function(record) {
  reductions <- purrr::map(
    record$prediction_manifest$prediction_path,
    read_mundlak_full_prediction_draw_file
  )
  key_reference <- reductions[[1L]]$grid
  if (!all(purrr::map_lgl(reductions, ~ identical(.x$grid, key_reference)))) {
    stop("Updated Mundlak prediction grids do not match within outcome.",
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
          file = path, algo = "sha256", serialize = FALSE
        ),
        within_model_draw = seq_len(reduction$draw_count)
      )
    }
  )
  if (
    nrow(draws) != record$n_draws ||
      nrow(draw_manifest) != record$n_draws || anyNA(draws)
  ) {
    stop("The updated Mundlak posterior mixture is incomplete.",
         call. = FALSE)
  }
  set.seed(record$averaging_seed)
  draw_order <- sample(seq_len(record$n_draws))
  draws <- draws[draw_order, , drop = FALSE]
  draw_manifest <- draw_manifest[draw_order, , drop = FALSE] |>
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
  write_atomic_qs(
    output,
    update_2026_mundlak_reporting_path(
      outcome = record$outcome,
      stage = "model_average"
    )
  )
}

write_update_2026_mundlak_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_mundlak_full_model_average_file(model_average_path)
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (is.null(data) || !is.data.frame(data)) {
    stop("The updated Mundlak reporting data are unavailable.",
         call. = FALSE)
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
  write_atomic_qs(
    report,
    update_2026_mundlak_reporting_path(
      outcome = outcome,
      stage = "reporting"
    )
  )
}

write_update_2026_mundlak_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "reporting"
    )
) {
  if (
    !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L
  ) {
    stop("The updated Mundlak reporting audit is incomplete.",
         call. = FALSE)
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_mundlak_full_reporting_file
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir, "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir, "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir, "common_exposure_support.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_mundlak_verification_checkpoint <- function(
    data_by_outcome,
    source_manifest,
    review_decisions,
    preparation_audit,
    specifications,
    accepted_fit_manifest,
    diagnostic_audit,
    integrated_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::mutate(
      mundlak_regularised_variance_priors = model_index == 4L
    ) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_mundlak_full_reporting_file
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(data) {
    data |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_audit$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_mundlak_name <- paste0(name_alignment$specification_name, " [mundlak]")
  weight_sums <- stacking_weights |>
    dplyr::summarise(weight_sum = sum(weight), models = dplyr::n(), .by = outcome)
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      .by = outcome
    )

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "all four combined analysis datasets are present and complete",
      "the combined-data preparation audit remains valid",
      "all five Mundlak registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist and have checksums",
      "all final model diagnostics pass",
      "integrated LOO contains 104 reductions and 15 audited high-k cases",
      "all 15 exact refits pass and compact corrections exist",
      "final corrected LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four model-averaged reporting audits pass"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      identical(names(data_by_outcome), outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L) &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) && all(analysis_data$vi > 0),
      all(preparation_audit$checks$passed),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_mundlak_name) &&
        all(name_alignment$diagnostic_name == expected_mundlak_name) &&
        all(name_alignment$loo_name == expected_mundlak_name) &&
        all(name_alignment$weight_name == expected_mundlak_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      nrow(diagnostic_audit$diagnostic_summary) == 104L &&
        all(diagnostic_audit$diagnostic_summary$diagnostic_pass),
      integrated_loo_audit$models == 104L &&
        integrated_loo_audit$exact_refits_required == 15L &&
        nrow(integrated_loo_audit$high_k_cases) == 15L,
      exact_audit$exact_refits == 15L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        all(file.exists(exact_audit$correction_manifest$correction_path)),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(file.exists(prediction_manifest$prediction_path)),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated Mundlak verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    review_decisions = review_decisions,
    preparation_checks = preparation_audit$checks,
    preparation_outcome_counts = preparation_audit$outcome_counts,
    preparation_update_study_counts = preparation_audit$update_study_counts,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_audit$diagnostic_summary,
    sampling_diagnostic_classification =
      diagnostic_audit$definitive_classification,
    sampling_parameter_failures = diagnostic_audit$parameter_failures,
    posterior_predictive_observation_flags =
      diagnostic_audit$ppc_observation_flags,
    posterior_predictive_outcome_summary =
      diagnostic_audit$ppc_outcome_summary,
    initial_integrated_loo_model_summary = integrated_loo_audit$model_summary,
    initial_integrated_loo_outcome_summary =
      integrated_loo_audit$outcome_summary,
    initial_high_pareto_k_cases = integrated_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_correction_manifest = exact_audit$correction_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final combined-data Mundlak checkpoint checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    review_decisions = "Eligibility decisions for updated-search studies.",
    preparation_checks = "Combined-data preparation invariants.",
    preparation_outcome_counts = "Analysis counts by primary outcome.",
    preparation_update_study_counts = "Updated-study contribution audit.",
    analysis_data = "Combined complete-case rows supplied to Mundlak models.",
    model_specifications = "Final 104-model Mundlak specification registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit path, stage, size, and checksum.",
    sampling_diagnostic_summary = "Final model-level sampling diagnostics.",
    sampling_diagnostic_classification = "Final diagnostic decisions.",
    sampling_parameter_failures = "Parameter-level diagnostic failures.",
    posterior_predictive_observation_flags = "Observation-level PPC flags.",
    posterior_predictive_outcome_summary = "Outcome-level PPC summary.",
    initial_integrated_loo_model_summary = "Conditional and integrated PSIS summaries.",
    initial_integrated_loo_outcome_summary = "Initial LOO audit by outcome.",
    initial_high_pareto_k_cases = "The 15 cases selected for exact refitting.",
    exact_refit_diagnostic_summary = "Diagnostics and exact corrections for 15 refits.",
    exact_correction_manifest = "Compact exact-correction paths and checksums.",
    final_loo_manifest = "Final exact-corrected LOO lineage.",
    final_loo_model_summary = "Final LOO summary for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    stacking_weights = "Final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    dose_response_predictions = "Mundlak model-averaged dose-response summaries.",
    dose_response_increments = "Increments from the 0.8 g/kg/day reference.",
    exposure_support_by_group = "Observed exposure support by RT group.",
    common_exposure_support = "Common exposure support by outcome.",
    pairwise_contrasts = "Model-averaged pairwise contrasts and probabilities.",
    pairwise_support = "Support limits for pairwise contrasts.",
    local_equivalence = "Local practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    figures_deferred = paste(
      "Cross-family primary/Mundlak/within figures are deferred until",
      "the combined-data within-family refit is complete."
    )
  )
}

make_update_2026_mundlak_pruning_plan <- function(
    accepted_fit_manifest,
    verification_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "mundlak"
    )
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "mundlak_update_2026"),
    "manifest.csv"
  )
  if (any(!file.exists(verification_files))) {
    stop("The updated Mundlak verification bundle is incomplete.",
         call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "mundlak"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated Mundlak archive.",
         call. = FALSE)
  }

  fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  if (any(!dir.exists(fit_roots))) {
    stop("An expected updated Mundlak fit directory is missing.",
         call. = FALSE)
  }
  fit_paths <- unlist(
    purrr::map(
      fit_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  )
  fit_paths <- normalizePath(fit_paths, winslash = "/", mustWork = TRUE)
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(fit_paths) != 113L || length(accepted_paths) != 104L ||
      anyDuplicated(fit_paths) || anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% fit_paths) ||
      sum(!fit_paths %in% accepted_paths) != 9L
  ) {
    stop("Updated Mundlak fit retention does not match its audited lineage.",
         call. = FALSE)
  }

  fit_manifest <- tibble::tibble(
    path = fit_paths,
    relative_path = substring(fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", fit_paths, fixed = TRUE) ~ "initial",
      grepl("/remedies/remedy_1/", fit_paths, fixed = TRUE) ~ "remedy_1",
      grepl("/remedies/remedy_2/", fit_paths, fixed = TRUE) ~ "remedy_2",
      TRUE ~ "unknown"
    ),
    accepted_fit = fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(fit_paths)$size),
    lean_action = "delete frozen full fit archive"
  ) |>
    dplyr::arrange(archive_stage, relative_path)
  if (
    any(fit_manifest$archive_stage == "unknown") ||
      !identical(
        as.integer(table(factor(
          fit_manifest$archive_stage,
          levels = c("initial", "remedy_1", "remedy_2")
        ))),
        c(104L, 8L, 1L)
      )
  ) {
    stop("The updated Mundlak archive stages are not as expected.",
         call. = FALSE)
  }

  compact_roots <- file.path(archive_root, c("loo", "reporting"))
  if (any(!dir.exists(compact_roots))) {
    stop("An expected compact Mundlak archive directory is missing.",
         call. = FALSE)
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(.x, recursive = TRUE, full.names = TRUE)
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[file.exists(compact_paths)]
  compact_manifest <- tibble::tibble(
    path = normalizePath(compact_paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(compact_paths)$size),
    action = "retain compact LOO/reporting archive"
  )

  summary <- tibble::tibble(
    mode = "lean",
    fit_files_deleted = nrow(fit_manifest),
    accepted_fit_files_deleted = sum(fit_manifest$accepted_fit),
    superseded_fit_files_deleted = sum(!fit_manifest$accepted_fit),
    fit_files_retained = 0L,
    gib_reclaimed = sum(fit_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained = sum(compact_manifest$bytes) / 1024^3,
    restoration = paste(
      "Reported results remain frozen; restoring full fits requires",
      "refitting the Mundlak family."
    )
  )

  list(
    summary = summary,
    fit_file_manifest = fit_manifest,
    compact_archive_manifest = compact_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_mundlak_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    fit_file_manifest = file.path(output_dir, "fit_file_manifest.csv"),
    compact_archive_manifest = file.path(
      output_dir,
      "compact_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(plan$fit_file_manifest, paths[["fit_file_manifest"]])
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_mundlak_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "mundlak", "pruning"
    )
) {
  required_confirmation <- "DELETE ALL UPDATED MUNDLAK FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "mundlak"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the expected Mundlak archive.",
         call. = FALSE)
  }
  selected <- plan$fit_file_manifest
  if (
    nrow(selected) != 113L || sum(selected$accepted_fit) != 104L ||
      sum(!selected$accepted_fit) != 9L ||
      !identical(
        as.integer(table(factor(
          selected$archive_stage,
          levels = c("initial", "remedy_1", "remedy_2")
        ))),
        c(104L, 8L, 1L)
      )
  ) {
    stop("The selected Mundlak pruning set is inconsistent with its audit.",
         call. = FALSE)
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  if (file.exists(receipt_path)) {
    stop("The updated Mundlak lean-pruning receipt already exists.",
         call. = FALSE)
  }
  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  permitted_subtree <- grepl(
    "/(candidates|remedies)/",
    selected_paths
  )
  if (
    any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted_subtree) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the permitted Mundlak fit subtree.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256")
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop("At least one Mundlak fit archive could not be deleted.",
         call. = FALSE)
  }
  if (
    any(!file.exists(plan$compact_archive_manifest$path)) ||
      !file.exists(plan$verification_manifest)
  ) {
    stop("A retained Mundlak compact artifact is unexpectedly missing.",
         call. = FALSE)
  }

  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, receipt_path)
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

# Combined-data achieved-intake sensitivity refit ----

make_update_2026_achieved_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated achieved-intake base registry is incomplete.",
         call. = FALSE)
  }

  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  remedy_keys <- unname(purrr::map_chr(primary_remedy_specifications, key))
  if (
    anyDuplicated(primary_keys) || anyDuplicated(remedy_keys) ||
      !identical(remedy_keys, "function::10") ||
      !all(remedy_keys %in% primary_keys)
  ) {
    stop("The updated achieved-intake primary replacement is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  historical_remedies <- c(
    "hypertrophy::11", "hypertrophy::14", "hypertrophy::18",
    "hypertrophy::23", "strength::6", "function::14"
  )
  specifications <- purrr::map(specifications, function(specification) {
    specification_key <- key(specification)
    historical_remedy <- specification_key %in% historical_remedies
    if (historical_remedy) {
      specification$adapt_delta <- max(
        as.numeric(specification$adapt_delta),
        0.999
      )
      specification$max_treedepth <- max(
        as.integer(specification$max_treedepth),
        15L
      )
    }
    specification$update_2026_sensitivity_base <- TRUE
    specification$update_2026_achieved_historical_remedy <- historical_remedy
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  historical_keys <- final_keys[purrr::map_lgl(
    specifications,
    ~ isTRUE(.x$update_2026_achieved_historical_remedy)
  )]
  function_ten <- specifications[[match("function::10", final_keys)]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !setequal(historical_keys, historical_remedies) ||
      !identical(as.numeric(function_ten$adapt_delta), 0.9995)
  ) {
    stop("The updated achieved-intake registry failed validation.",
         call. = FALSE)
  }

  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

update_2026_achieved_local_config <- function(specification, config) {
  local_config <- config
  historical_remedy <- isTRUE(
    specification$update_2026_achieved_historical_remedy
  )
  if (historical_remedy) {
    local_config$iter_warmup <- 2000L
    local_config$iter_sampling <- 4000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      15L
    )
  }
  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the updated achieved-intake configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed +
      (if (historical_remedy) 1000L else 0L) +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

validate_update_2026_achieved_inputs <- function(
    specifications,
    exposure_data_by_outcome,
    priors_by_outcome,
    config
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.list(exposure_data_by_outcome) ||
      !all(c("achieved", "assignment", "audit") %in%
        names(exposure_data_by_outcome)) ||
      nrow(exposure_data_by_outcome$audit) != 4L ||
      !all(exposure_data_by_outcome$audit$identical_common_sample)
  ) {
    stop("The updated achieved/assignment common sample is invalid.",
         call. = FALSE)
  }

  checks <- purrr::map_dfr(specifications, function(specification) {
    local_config <- update_2026_achieved_local_config(
      specification,
      config
    )
    input <- prepare_exposure_full_input(
      specification = specification,
      exposure_data_by_outcome = exposure_data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      family = "achieved",
      config = local_config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      rows = nrow(input$data),
      studies = dplyr::n_distinct(input$data$study_no),
      arms = dplyr::n_distinct(input$data$arm_no),
      achieved_minimum = min(input$data$protein_intake_raw),
      achieved_maximum = max(input$data$protein_intake_raw),
      basis_dimension = if (is.null(input$basis_spec)) {
        0L
      } else {
        as.integer(input$basis_spec$basis_dimension)
      },
      historical_remedy_schedule = isTRUE(
        specification$update_2026_achieved_historical_remedy
      ),
      adapt_delta = max(
        input$specification$adapt_delta,
        local_config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        local_config$max_treedepth_floor
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  historical <- checks |>
    dplyr::filter(.data$historical_remedy_schedule)
  standard <- checks |>
    dplyr::filter(!.data$historical_remedy_schedule)
  counts <- table(factor(checks$outcome, levels = outcomes))
  sample_counts <- checks |>
    dplyr::distinct(
      outcome, rows, studies, arms,
      achieved_minimum, achieved_maximum
    )
  expected_historical <- c(
    "hypertrophy::11", "hypertrophy::14", "hypertrophy::18",
    "hypertrophy::23", "strength::6", "function::14"
  )
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      nrow(sample_counts) != 4L || any(sample_counts$rows < 1L) ||
      any(!is.finite(sample_counts$achieved_minimum)) ||
      any(!is.finite(sample_counts$achieved_maximum)) ||
      !setequal(
        paste(historical$outcome, historical$model_index, sep = "::"),
        expected_historical
      ) ||
      any(historical$iter_warmup != 2000L) ||
      any(historical$iter_sampling != 4000L) ||
      any(historical$adapt_delta < 0.999) ||
      any(historical$max_treedepth < 15L) ||
      any(standard$iter_warmup != 1500L) ||
      any(standard$iter_sampling != 2000L) ||
      config$chains != 4L || config$parallel_chains != 4L
  ) {
    stop("The updated achieved-intake preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_achieved_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "achieved",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete achieved-intake preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_achieved_candidate_to_file <- function(
    specification,
    exposure_data_by_outcome,
    priors_by_outcome,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "achieved",
      "candidates"
    )
) {
  local_config <- update_2026_achieved_local_config(specification, config)
  input <- prepare_exposure_full_input(
    specification = specification,
    exposure_data_by_outcome = exposure_data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    family = "achieved",
    config = local_config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "achieved"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated achieved-intake fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated achieved-intake fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

compile_update_2026_achieved_diagnostic_audit <- function(diagnostics) {
  split_diagnostics <- split(
    diagnostics,
    purrr::map_chr(diagnostics, ~ .x$summary$outcome[[1L]])
  )
  required <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !setequal(names(split_diagnostics), required) ||
      any(vapply(split_diagnostics[required], length, integer(1L)) != 26L)
  ) {
    stop("Updated achieved-intake diagnostics are incomplete.",
         call. = FALSE)
  }
  compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
}

make_update_2026_achieved_remedy_config <- function(
    initial_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    extended_key = "hypertrophy::14",
    extended_iter_warmup = 3000L,
    extended_iter_sampling = 6000L,
    seed_offset = 1000L
) {
  config <- initial_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$extended_remedy_key <- as.character(extended_key)
  config$extended_iter_warmup <- as.integer(extended_iter_warmup)
  config$extended_iter_sampling <- as.integer(extended_iter_sampling)
  config$seed <- as.integer(initial_config$seed + seed_offset)
  if (
    config$iter_warmup < initial_config$iter_warmup ||
      config$iter_sampling <= initial_config$iter_sampling ||
      config$extended_remedy_key != "hypertrophy::14" ||
      config$extended_iter_warmup <= config$iter_warmup ||
      config$extended_iter_sampling <= config$iter_sampling ||
      config$adapt_delta_floor != initial_config$adapt_delta_floor ||
      config$max_treedepth_floor != initial_config$max_treedepth_floor
  ) {
    stop("The updated achieved-intake remedy configuration is invalid.",
         call. = FALSE)
  }
  config
}

select_update_2026_achieved_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  expected_keys <- c(
    "hypertrophy::14",
    "hypertrophy::26",
    "strength::14",
    "strength::24"
  )
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(
      outcome,
      model_index = as.integer(model_index),
      key = paste(outcome, model_index, sep = "::"),
      divergent_transitions,
      maximum_treedepth_transitions
    )
  if (
    nrow(failed) != 4L || anyDuplicated(failed$key) ||
      !setequal(failed$key, expected_keys) ||
      any(failed$divergent_transitions != 0L) ||
      any(failed$maximum_treedepth_transitions != 0L)
  ) {
    stop("The updated achieved-intake remedy set no longer matches its audit.",
         call. = FALSE)
  }

  selected <- purrr::keep(specifications, function(specification) {
    paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    ) %in% failed$key
  }) |>
    purrr::map(function(specification) {
      specification_key <- paste(
        specification$outcome,
        specification$model_index,
        sep = "::"
      )
      if (specification_key == "hypertrophy::14") {
        if (!isTRUE(
          specification$update_2026_achieved_historical_remedy
        )) {
          stop("Achieved hypertrophy model 14 lost its historical schedule.",
               call. = FALSE)
        }
        specification$update_2026_achieved_remedy <-
          "longer_sampling_3000_6000_after_historical_schedule"
      } else {
        if (isTRUE(
          specification$update_2026_achieved_historical_remedy
        )) {
          stop("An unexpected historical achieved model requires remedy.",
               call. = FALSE)
        }
        specification$update_2026_achieved_remedy <- "longer_sampling"
      }
      specification$update_2026_achieved_remedy_stage <- "remedy_1"
      specification
    })
  selected_keys <- unname(purrr::map_chr(
    selected,
    ~ paste(.x$outcome, .x$model_index, sep = "::")
  ))
  if (
    length(selected) != 4L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The updated achieved-intake remedy specifications are incomplete.",
         call. = FALSE)
  }
  names(selected) <- paste0(
    purrr::map_chr(selected, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(selected, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  selected
}

update_2026_achieved_remedy_local_config <- function(
    specification,
    config
) {
  local_config <- config
  specification_key <- paste(
    specification$outcome,
    specification$model_index,
    sep = "::"
  )
  if (specification_key == config$extended_remedy_key) {
    local_config$iter_warmup <- config$extended_iter_warmup
    local_config$iter_sampling <- config$extended_iter_sampling
  }
  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in an achieved-intake remedy configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

validate_update_2026_achieved_remedy_inputs <- function(
    specifications,
    exposure_data_by_outcome,
    priors_by_outcome,
    config
) {
  checks <- purrr::map_dfr(specifications, function(specification) {
    local_config <- update_2026_achieved_remedy_local_config(
      specification,
      config
    )
    input <- prepare_exposure_full_input(
      specification = specification,
      exposure_data_by_outcome = exposure_data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      family = "achieved",
      config = local_config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      remedy = specification$update_2026_achieved_remedy,
      adapt_delta = max(
        input$specification$adapt_delta,
        local_config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        local_config$max_treedepth_floor
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  expected_keys <- c(
    "hypertrophy::14",
    "hypertrophy::26",
    "strength::14",
    "strength::24"
  )
  extended <- checks |>
    dplyr::filter(.data$outcome == "hypertrophy", .data$model_index == 14L)
  standard <- checks |>
    dplyr::filter(!(.data$outcome == "hypertrophy" &
      .data$model_index == 14L))
  if (
    nrow(checks) != 4L || anyDuplicated(keys) ||
      !setequal(keys, expected_keys) || !all(checks$valid) ||
      nrow(extended) != 1L ||
      extended$remedy[[1L]] !=
        "longer_sampling_3000_6000_after_historical_schedule" ||
      extended$iter_warmup[[1L]] != 3000L ||
      extended$iter_sampling[[1L]] != 6000L ||
      any(standard$remedy != "longer_sampling") ||
      any(standard$iter_warmup != 2000L) ||
      any(standard$iter_sampling != 4000L)
  ) {
    stop("The updated achieved-intake remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_achieved_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "achieved", "remedy_1",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 4L || !all(checks$valid)) {
    stop("Cannot write an incomplete achieved-intake remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_achieved_remedy_to_file <- function(
    specification,
    exposure_data_by_outcome,
    priors_by_outcome,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "achieved",
      "remedies", "remedy_1"
    )
) {
  local_config <- update_2026_achieved_remedy_local_config(
    specification,
    config
  )
  input <- prepare_exposure_full_input(
    specification = specification,
    exposure_data_by_outcome = exposure_data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    family = "achieved",
    config = local_config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "achieved"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an achieved-intake remedy fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("An achieved-intake remedy fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

compile_update_2026_achieved_remedy_diagnostic_audit <- function(
    initial_diagnostics,
    remedy_diagnostics
) {
  diagnostic_key <- function(x) {
    paste(x$summary$outcome[[1L]], x$summary$model_index[[1L]], sep = "::")
  }
  if (length(initial_diagnostics) != 104L || length(remedy_diagnostics) != 4L) {
    stop("Updated achieved-intake remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  initial_keys <- unname(purrr::map_chr(initial_diagnostics, diagnostic_key))
  remedy_keys <- unname(purrr::map_chr(remedy_diagnostics, diagnostic_key))
  expected_keys <- initial_diagnostics |>
    purrr::keep(~ !.x$summary$diagnostic_pass[[1L]]) |>
    purrr::map_chr(diagnostic_key)
  if (
    anyDuplicated(initial_keys) || anyDuplicated(remedy_keys) ||
      !setequal(remedy_keys, expected_keys) ||
      !all(remedy_keys %in% initial_keys)
  ) {
    stop("Updated achieved-intake diagnostic replacement is invalid.",
         call. = FALSE)
  }
  final_diagnostics <- initial_diagnostics
  final_diagnostics[match(remedy_keys, initial_keys)] <- remedy_diagnostics
  compile_update_2026_achieved_diagnostic_audit(final_diagnostics)
}

make_update_2026_achieved_accepted_fit_manifest <- function(
    initial_fit_paths,
    remedy_fit_paths,
    final_diagnostic_audit
) {
  if (
    nrow(final_diagnostic_audit$diagnostic_summary) != 104L ||
      any(!final_diagnostic_audit$diagnostic_summary$diagnostic_pass)
  ) {
    stop("Updated achieved-intake fits cannot be accepted before diagnostics pass.",
         call. = FALSE)
  }
  path_table <- function(paths, stage, expected) {
    paths <- unname(unlist(paths, use.names = FALSE))
    table <- tibble::tibble(
      outcome = basename(dirname(paths)),
      model_index = as.integer(stringr::str_match(
        basename(paths),
        "model_([0-9]+)\\.qs$"
      )[, 2]),
      stage = stage,
      fit_path = normalizePath(paths, winslash = "/", mustWork = TRUE)
    )
    if (
      nrow(table) != expected || anyNA(table$model_index) ||
        anyDuplicated(table[c("outcome", "model_index")]) ||
        any(file.info(table$fit_path)$size <= 0)
    ) {
      stop("An updated achieved-intake fit stage is incomplete.",
           call. = FALSE)
    }
    table
  }

  initial <- path_table(initial_fit_paths, "initial", 104L)
  remedy <- path_table(remedy_fit_paths, "remedy_1", 4L)
  keys <- c("outcome", "model_index")
  remedy_keys <- paste(remedy$outcome, remedy$model_index, sep = "::")
  expected_remedy_keys <- c(
    "hypertrophy::14",
    "hypertrophy::26",
    "strength::14",
    "strength::24"
  )
  if (!setequal(remedy_keys, expected_remedy_keys)) {
    stop("Updated achieved-intake accepted replacements are invalid.",
         call. = FALSE)
  }

  accepted <- initial |>
    dplyr::anti_join(remedy[keys], by = keys) |>
    dplyr::bind_rows(remedy) |>
    dplyr::left_join(
      final_diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, model_name = model),
      by = keys,
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      family = "achieved",
      bytes = file.info(fit_path)$size,
      sha256 = purrr::map_chr(
        fit_path,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      )
    ) |>
    dplyr::select(
      family, outcome, model_index, model_name, stage,
      fit_path, bytes, sha256
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  stage_counts <- table(factor(
    accepted$stage,
    levels = c("initial", "remedy_1")
  ))
  if (
    nrow(accepted) != 104L ||
      anyDuplicated(accepted[c("outcome", "model_index")]) ||
      anyNA(accepted$model_name) || any(!file.exists(accepted$fit_path)) ||
      !identical(as.integer(stage_counts), c(100L, 4L)) ||
      any(nchar(accepted$sha256) != 64L)
  ) {
    stop("The updated achieved-intake accepted-fit manifest is invalid.",
         call. = FALSE)
  }
  accepted
}

write_update_2026_achieved_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "achieved",
      "accepted_fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(manifest[c("outcome", "model_index")]) ||
      any(!file.exists(manifest$fit_path)) ||
      any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an invalid achieved-intake fit manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

compute_update_2026_achieved_integrated_loo_to_file <- function(
    fit_record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "achieved", "loo",
      "integrated"
    )
) {
  result <- compute_exposure_full_file_integrated_loo(
    fit_record = fit_record,
    family = "achieved",
    loo_draws = loo_draws,
    k_threshold = k_threshold,
    cores = cores
  )
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(unname(as.character(result$family)), "achieved") ||
      !identical(result$remedy, "integrated_observation_level_psis")
  ) {
    stop("Cannot archive an invalid updated achieved-intake LOO result.",
         call. = FALSE)
  }
  output_dir <- file.path(output_root, result$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf(
      "model_%02d_integrated_loo.qs",
      as.integer(result$model_index)
    )
  )
  write_atomic_qs(result, output_file)
}

write_update_2026_achieved_integrated_loo_audit <- function(
    audit,
    manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "loo_psis"
    )
) {
  if (
    !isTRUE(audit$integrated_psis_complete) || audit$models != 104L ||
      nrow(manifest) != 104L ||
      any(!file.exists(manifest$integrated_loo_path))
  ) {
    stop("The updated achieved-intake integrated-LOO audit is incomplete.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

validate_update_2026_achieved_exact_inputs <- function(
    records,
    integrated_audit,
    config
) {
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    !is.list(records) || length(records) != expected || expected != 19L
  ) {
    stop("Expected 19 updated achieved-intake exact-LOO records.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    integrated <- qs::qread(record$integrated_loo_path)
    observation <- as.integer(record$observation)
    if (
      !inherits(integrated$loo, "loo") ||
        !identical(unname(as.character(integrated$family)), "achieved") ||
        !identical(integrated$outcome, record$outcome) ||
        as.integer(integrated$model_index) != record$model_index ||
        observation < 1L ||
        observation > nrow(integrated$observation_manifest)
    ) {
      stop("An updated achieved-intake exact record is invalid.",
           call. = FALSE)
    }
    held_out <- integrated$observation_manifest[observation, , drop = FALSE]
    observed_k <- as.numeric(
      loo::pareto_k_values(integrated$loo)[[observation]]
    )
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      observation = observation,
      study_no = as.character(held_out$study_no[[1L]]),
      arm_no = as.integer(held_out$arm_no[[1L]]),
      outcome_no = as.integer(held_out$outcome_no[[1L]]),
      pareto_k = observed_k,
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      original_fit_exists = file.exists(record$original_fit_path),
      integrated_loo_exists = file.exists(record$integrated_loo_path),
      valid = is.finite(observed_k) && observed_k > 0.7
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  if (
    nrow(checks) != expected || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(expected)) ||
      !identical(as.integer(outcome_counts), c(17L, 2L, 0L, 0L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$integrated_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L
  ) {
    stop("The updated achieved-intake exact-LOO preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_achieved_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "achieved", "loo_exact",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 19L || !all(checks$valid)) {
    stop("Cannot write an incomplete achieved-intake exact preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_achieved_exact_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "achieved", "loo",
    "exact", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

fit_and_score_update_2026_achieved_exact_to_file <- function(
    record,
    config,
    output_file = update_2026_achieved_exact_file_path(record)
) {
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !identical(unname(as.character(candidate$family)), "achieved") ||
      !identical(candidate$outcome, record$outcome) ||
      as.integer(candidate$model_index) != record$model_index ||
      !inherits(integrated$loo, "loo") ||
      !identical(unname(as.character(integrated$family)), "achieved") ||
      !identical(integrated$outcome, record$outcome) ||
      as.integer(integrated$model_index) != record$model_index ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      loo::pareto_k_values(integrated$loo)[[observation]] <= 0.7
  ) {
    stop("An updated achieved-intake exact-LOO record failed validation.",
         call. = FALSE)
  }

  original_data <- candidate$model$data
  training_data <- original_data[-observation, , drop = FALSE]
  original_model <- candidate$model
  sample_prior <- attr(original_model$prior, "sample_prior")
  original_model$prior <- normalize_primary_brms_priors(
    original_model$prior
  )
  attr(original_model$prior, "sample_prior") <- sample_prior
  specification <- candidate$specification
  specification$adapt_delta <- max(
    specification$adapt_delta,
    config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth,
    config$max_treedepth_floor
  )

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
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  refit <- list(
    family = "achieved",
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    original_data = original_data,
    integrated_result = integrated,
    specification = specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  correction <- compute_sensitivity_mundlak_exact_loo_case(refit)
  if (
    as.integer(correction$case_id) != record$case_id ||
      !inherits(correction$loo, "loo") ||
      nrow(correction$exact_refit_summary) != 1L
  ) {
    stop("The updated achieved-intake exact correction is invalid.",
         call. = FALSE)
  }
  write_atomic_qs(correction, output_file)
}

audit_update_2026_achieved_exact_files <- function(
    paths,
    integrated_audit
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    length(paths) != expected || expected != 19L || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0)
  ) {
    stop("The updated achieved-intake exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  manifest <- purrr::map2_dfr(paths, corrections, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      observation =
        as.integer(result$exact_refit_summary$observation[[1L]]),
      correction_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys)
  ) {
    stop("The updated achieved-intake exact audit is invalid.",
         call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = sum(diagnostic_summary$diagnostic_pass),
    exact_refits_requiring_remediation =
      sum(!diagnostic_summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass =
      all(diagnostic_summary$diagnostic_pass)
  )
}

write_update_2026_achieved_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "loo_exact"
    )
) {
  if (
    audit$exact_refits != 19L ||
      nrow(audit$diagnostic_summary) != 19L ||
      nrow(audit$correction_manifest) != 19L
  ) {
    stop("Cannot write an incomplete achieved-intake exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

make_update_2026_achieved_exact_remedy_config <- function(
    exact_config,
    adapt_delta = 0.9995,
    seed_offset = 1000L
) {
  config <- exact_config
  config$adapt_delta_floor <- max(
    as.numeric(exact_config$adapt_delta_floor),
    as.numeric(adapt_delta)
  )
  config$seed <- as.integer(exact_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.9995 ||
      config$max_treedepth_floor != exact_config$max_treedepth_floor
  ) {
    stop("The achieved-intake exact-remedy configuration is invalid.",
         call. = FALSE)
  }
  config
}

select_update_2026_achieved_exact_remedy_record <- function(
    records,
    exact_audit
) {
  failed <- exact_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass)
  if (
    nrow(failed) != 1L || failed$outcome[[1L]] != "whole_body" ||
      failed$model_index[[1L]] != 20L ||
      failed$observation[[1L]] != 62L ||
      failed$divergent_transitions[[1L]] != 1L ||
      failed$maximum_treedepth_transitions[[1L]] != 0L ||
      !isTRUE(failed$convergence_pass[[1L]])
  ) {
    stop("The achieved-intake exact remedy no longer matches its audit.",
         call. = FALSE)
  }
  selected <- purrr::keep(records, function(record) {
    identical(record$outcome, "whole_body") &&
      as.integer(record$model_index) == 20L &&
      as.integer(record$observation) == 62L
  })
  if (
    length(selected) != 1L ||
      as.integer(selected[[1L]]$case_id) != 11L
  ) {
    stop("Could not uniquely select achieved-intake exact case 11.",
         call. = FALSE)
  }
  names(selected) <- "case_011_whole_body_model_20_observation_0062"
  selected
}

validate_update_2026_achieved_exact_remedy_input <- function(
    record,
    config
) {
  if (!is.list(record) || length(record) != 1L) {
    stop("Expected one achieved-intake exact-remedy record.",
         call. = FALSE)
  }
  record <- record[[1L]]
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  observed_k <- as.numeric(
    loo::pareto_k_values(integrated$loo)[[observation]]
  )
  effective_adapt_delta <- max(
    as.numeric(candidate$specification$adapt_delta),
    as.numeric(config$adapt_delta_floor)
  )
  effective_max_treedepth <- max(
    as.integer(candidate$specification$max_treedepth),
    as.integer(config$max_treedepth_floor)
  )
  check <- tibble::tibble(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    remedy = "adapt_delta_0.9995_after_exact_refit",
    adapt_delta = effective_adapt_delta,
    max_treedepth = effective_max_treedepth,
    iter_warmup = as.integer(config$iter_warmup),
    iter_sampling = as.integer(config$iter_sampling),
    pareto_k = observed_k,
    original_fit_exists = file.exists(record$original_fit_path),
    integrated_loo_exists = file.exists(record$integrated_loo_path),
    valid = TRUE
  )
  if (
    check$case_id[[1L]] != 11L ||
      check$outcome[[1L]] != "whole_body" ||
      check$model_index[[1L]] != 20L ||
      check$observation[[1L]] != 62L ||
      check$adapt_delta[[1L]] != 0.9995 ||
      check$max_treedepth[[1L]] != 15L ||
      check$iter_warmup[[1L]] != 2000L ||
      check$iter_sampling[[1L]] != 6000L ||
      !is.finite(observed_k) || observed_k <= 0.7 ||
      !check$original_fit_exists[[1L]] ||
      !check$integrated_loo_exists[[1L]]
  ) {
    stop("The achieved-intake exact-remedy preflight is invalid.",
         call. = FALSE)
  }
  check
}

write_update_2026_achieved_exact_remedy_preflight <- function(
    check,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "achieved", "loo_exact",
      "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(check) != 1L || !isTRUE(check$valid[[1L]])) {
    stop("Cannot write an invalid achieved-intake exact-remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(check, output_file)
  output_file
}

update_2026_achieved_exact_remedy_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "achieved", "loo",
    "exact_remedy", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

audit_update_2026_achieved_exact_remedy_file <- function(path, record) {
  path <- unname(unlist(path, use.names = FALSE))
  if (
    is.list(record) && is.null(record$case_id) && length(record) == 1L
  ) {
    record <- record[[1L]]
  }
  if (
    length(path) != 1L || !file.exists(path) || file.info(path)$size <= 0
  ) {
    stop("The achieved-intake exact-remedy file is unavailable.",
         call. = FALSE)
  }
  correction <- qs::qread(path)
  summary <- correction$exact_refit_summary
  if (
    as.integer(correction$case_id) != 11L || nrow(summary) != 1L ||
      summary$outcome[[1L]] != "whole_body" ||
      summary$model_index[[1L]] != 20L ||
      summary$observation[[1L]] != 62L ||
      as.integer(record$case_id) != 11L
  ) {
    stop("The achieved-intake exact-remedy result is invalid.",
         call. = FALSE)
  }
  manifest <- tibble::tibble(
    case_id = 11L,
    outcome = summary$outcome[[1L]],
    model_index = as.integer(summary$model_index[[1L]]),
    model_name = correction$model_name,
    observation = as.integer(summary$observation[[1L]]),
    correction_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = digest::digest(
      file = path,
      algo = "sha256",
      serialize = FALSE
    )
  )
  list(
    diagnostic_summary = summary,
    correction_manifest = manifest,
    exact_refits = 1L,
    exact_refits_passing = sum(summary$diagnostic_pass),
    exact_refits_requiring_remediation = sum(!summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

write_update_2026_achieved_exact_remedy_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "loo_exact",
      "remedy_1"
    )
) {
  if (
    audit$exact_refits != 1L ||
      nrow(audit$diagnostic_summary) != 1L ||
      nrow(audit$correction_manifest) != 1L
  ) {
    stop("Cannot write an invalid achieved-intake exact-remedy audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

make_update_2026_achieved_final_exact_audit <- function(
    initial_audit,
    remedy_audit
) {
  if (
    initial_audit$exact_refits != 19L ||
      initial_audit$exact_refits_requiring_remediation != 1L ||
      remedy_audit$exact_refits != 1L ||
      !isTRUE(remedy_audit$all_exact_refit_diagnostics_pass)
  ) {
    stop("The achieved-intake exact audits are not ready to combine.",
         call. = FALSE)
  }
  key <- function(data) {
    paste(data$outcome, data$model_index, data$observation, sep = "::")
  }
  remedy_key <- key(remedy_audit$diagnostic_summary)
  expected_key <- "whole_body::20::62"
  if (!identical(remedy_key, expected_key)) {
    stop("The achieved-intake exact remedy has the wrong identity.",
         call. = FALSE)
  }
  diagnostic_summary <- initial_audit$diagnostic_summary |>
    dplyr::filter(key(dplyr::pick(dplyr::everything())) != remedy_key) |>
    dplyr::bind_rows(remedy_audit$diagnostic_summary) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  correction_manifest <- initial_audit$correction_manifest |>
    dplyr::filter(.data$case_id != 11L) |>
    dplyr::bind_rows(remedy_audit$correction_manifest) |>
    dplyr::arrange(case_id)
  if (
    nrow(diagnostic_summary) != 19L ||
      anyDuplicated(key(diagnostic_summary)) ||
      any(!diagnostic_summary$diagnostic_pass) ||
      nrow(correction_manifest) != 19L ||
      anyDuplicated(correction_manifest$case_id) ||
      !identical(correction_manifest$case_id, seq_len(19L)) ||
      any(!file.exists(correction_manifest$correction_path)) ||
      correction_manifest$correction_path[
        correction_manifest$case_id == 11L
      ] != remedy_audit$correction_manifest$correction_path[[1L]]
  ) {
    stop("The final achieved-intake exact audit is invalid.",
         call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = correction_manifest,
    exact_refits = 19L,
    exact_refits_passing = 19L,
    exact_refits_requiring_remediation = 0L,
    all_exact_refit_diagnostics_pass = TRUE,
    remedied_case_id = 11L
  )
}

write_update_2026_achieved_final_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "loo_exact",
      "final"
    )
) {
  if (
    audit$exact_refits != 19L ||
      !isTRUE(audit$all_exact_refit_diagnostics_pass) ||
      nrow(audit$diagnostic_summary) != 19L ||
      nrow(audit$correction_manifest) != 19L
  ) {
    stop("Cannot write an incomplete final achieved-intake exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

read_update_2026_achieved_exact_corrections <- function(manifest) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 19L ||
      any(!file.exists(manifest$correction_path))
  ) {
    stop("The achieved-intake exact-correction manifest is invalid.",
         call. = FALSE)
  }
  corrections <- purrr::map(manifest$correction_path, qs::qread)
  corrections <- corrections[order(purrr::map_int(corrections, "case_id"))]
  if (
    !identical(
      unname(purrr::map_int(corrections, "case_id")),
      seq_len(19L)
    ) ||
      any(!purrr::map_lgl(corrections, ~ inherits(.x$loo, "loo")))
  ) {
    stop("The achieved-intake exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections
}

make_update_2026_achieved_final_loo_audit <- function(
    integrated_manifest,
    exact_corrections,
    exact_audit
) {
  if (!isTRUE(exact_audit$all_exact_refit_diagnostics_pass)) {
    stop("Final achieved-intake LOO requires passing exact diagnostics.",
         call. = FALSE)
  }
  final <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest = integrated_manifest,
    exact_results = exact_corrections,
    family = "achieved"
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- integrated_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "integrated_observation_level_psis_plus_exact",
        "integrated_observation_level_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 19L
  ) {
    stop("The final achieved-intake LOO is not stacking-ready.",
         call. = FALSE)
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 19L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_achieved_final_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 19L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop("Cannot write an incomplete final achieved-intake LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_achieved_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "achieved",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L ||
      anyDuplicated(paste(weights$outcome, weights$model_index)) ||
      length(outcome_sums) != 4L ||
      any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop("The updated achieved-intake stacking weights are invalid.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_achieved_reporting_preflight <- function(records) {
  if (!is.list(records) || length(records) < 4L) {
    stop("The achieved-intake prediction allocation is incomplete.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    tibble::tibble(
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      family = unname(as.character(record$family)),
      prediction_type = unname(as.character(record$prediction_type)),
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      protein_max = as.numeric(record$protein_max),
      grid_step = as.numeric(record$grid_step),
      fit_path = record$fit_path,
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path) &&
        identical(unname(as.character(record$family)), "achieved") &&
        identical(
          unname(as.character(record$prediction_type)),
          "exposure"
        ) &&
        is.finite(record$protein_max) &&
        record$protein_max > record$protein_min
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      allocated_weight = sum(weight),
      common_protein_max = dplyr::first(protein_max),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      any(!is.finite(outcome_audit$allocated_weight)) ||
      any(outcome_audit$allocated_weight <= 0) ||
      any(outcome_audit$allocated_weight > 1 + 1e-8) ||
      any(!is.finite(checks$weight)) || any(checks$weight < 0) ||
      any(!is.finite(outcome_audit$common_protein_max)) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop("The achieved-intake reporting preflight failed.",
         call. = FALSE)
  }
  checks
}

write_update_2026_achieved_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "achieved", "reporting",
      "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(checks$allocated_draws, checks$outcome, sum)
  outcome_weights <- tapply(checks$weight, checks$outcome, sum)
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      length(outcome_weights) != 4L ||
      any(!is.finite(outcome_weights)) || any(outcome_weights <= 0) ||
      any(outcome_weights > 1 + 1e-8) || !all(checks$valid)
  ) {
    stop("Cannot write an invalid achieved-intake reporting preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_achieved_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "achieved",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_achieved_prediction_draw_file <- function(
    record,
    data_by_outcome
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !identical(unname(as.character(record$family)), "achieved") ||
      !identical(unname(as.character(record$prediction_type)), "exposure")
  ) {
    stop("An updated achieved-intake prediction record is invalid.",
         call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "achieved") ||
      !identical(candidate$outcome, unname(as.character(record$outcome))) ||
      as.integer(candidate$model_index) != as.integer(record$model_index) ||
      !identical(candidate$model_name, unname(as.character(record$model_name)))
  ) {
    stop("An updated achieved-intake prediction fit does not match its record.",
         call. = FALSE)
  }
  prediction_data <- data_by_outcome[[record$outcome]]
  if (is.null(prediction_data) || !is.data.frame(prediction_data)) {
    stop("The updated achieved-intake prediction data are unavailable.",
         call. = FALSE)
  }
  grid <- make_sensitivity_prediction_grid(
    fitted_candidate = candidate,
    data = prediction_data,
    prediction_type = "exposure",
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
    stop("Updated achieved-intake prediction draws are incomplete.",
         call. = FALSE)
  }
  reduction <- list(
    family = "achieved",
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
      prediction_type = "exposure",
      reference_intake = NULL
    )
  )
  write_atomic_qs(
    reduction,
    update_2026_achieved_reporting_path(
      outcome = record$outcome,
      stage = "prediction_draws",
      model_index = record$model_index
    )
  )
}

write_update_2026_achieved_model_average_file <- function(record) {
  reductions <- purrr::map(
    record$prediction_manifest$prediction_path,
    read_sensitivity_full_prediction_draw_file,
    family = "achieved"
  )
  key_reference <- reductions[[1L]]$grid
  if (!all(purrr::map_lgl(
    reductions, ~ identical(.x$grid, key_reference)
  ))) {
    stop("Updated achieved-intake prediction grids do not match within outcome.",
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
          file = path, algo = "sha256", serialize = FALSE
        ),
        within_model_draw = seq_len(reduction$draw_count)
      )
    }
  )
  if (
    nrow(draws) != record$n_draws ||
      nrow(draw_manifest) != record$n_draws || anyNA(draws)
  ) {
    stop("The updated achieved-intake posterior mixture is incomplete.",
         call. = FALSE)
  }
  set.seed(record$averaging_seed)
  draw_order <- sample(seq_len(record$n_draws))
  draws <- draws[draw_order, , drop = FALSE]
  draw_manifest <- draw_manifest[draw_order, , drop = FALSE] |>
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
  write_atomic_qs(
    output,
    update_2026_achieved_reporting_path(
      outcome = record$outcome,
      stage = "model_average"
    )
  )
}

write_update_2026_achieved_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path, family = "achieved"
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (is.null(data) || !is.data.frame(data)) {
    stop("The updated achieved-intake reporting data are unavailable.",
         call. = FALSE)
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
    family = "achieved",
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
  write_atomic_qs(
    report,
    update_2026_achieved_reporting_path(
      outcome = outcome,
      stage = "reporting"
    )
  )
}

write_update_2026_achieved_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "reporting"
    )
) {
  if (
    !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L
  ) {
    stop("The updated achieved-intake reporting audit is incomplete.",
         call. = FALSE)
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "achieved"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir, "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir, "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir, "common_exposure_support.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_achieved_verification_checkpoint <- function(
    data_by_outcome,
    exposure_audit,
    source_manifest,
    review_decisions,
    preparation_audit,
    specifications,
    accepted_fit_manifest,
    diagnostic_audit,
    integrated_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "achieved"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(data) {
    data |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_audit$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_achieved_name <- paste0(
    name_alignment$specification_name,
    " [achieved]"
  )
  weight_sums <- stacking_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      .by = outcome
    )

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "all four achieved-intake common-sample datasets are complete",
      "the combined-data preparation and exposure audits remain valid",
      "all five achieved-intake registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist and have checksums",
      "all final model diagnostics pass",
      "integrated LOO contains 104 reductions and 19 audited high-k cases",
      "all 19 final exact refits pass and compact corrections exist",
      "final corrected LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four model-averaged reporting audits pass"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      identical(names(data_by_outcome), outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L) &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) && all(analysis_data$vi > 0) &&
        all(is.finite(analysis_data$protein_intake_raw)),
      all(preparation_audit$checks$passed) &&
        nrow(exposure_audit) == 4L &&
        all(exposure_audit$identical_common_sample) &&
        setequal(exposure_audit$outcome, outcomes),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_achieved_name) &&
        all(name_alignment$diagnostic_name == expected_achieved_name) &&
        all(name_alignment$loo_name == expected_achieved_name) &&
        all(name_alignment$weight_name == expected_achieved_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(accepted_fit_manifest$family == "achieved") &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      nrow(diagnostic_audit$diagnostic_summary) == 104L &&
        all(diagnostic_audit$diagnostic_summary$diagnostic_pass),
      integrated_loo_audit$models == 104L &&
        integrated_loo_audit$exact_refits_required == 19L &&
        nrow(integrated_loo_audit$high_k_cases) == 19L,
      exact_audit$exact_refits == 19L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        all(file.exists(exact_audit$correction_manifest$correction_path)),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(file.exists(prediction_manifest$prediction_path)),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated achieved-intake verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    review_decisions = review_decisions,
    preparation_checks = preparation_audit$checks,
    preparation_outcome_counts = preparation_audit$outcome_counts,
    preparation_update_study_counts = preparation_audit$update_study_counts,
    exposure_data_audit = exposure_audit,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_audit$diagnostic_summary,
    sampling_diagnostic_classification =
      diagnostic_audit$definitive_classification,
    sampling_parameter_failures = diagnostic_audit$parameter_failures,
    posterior_predictive_observation_flags =
      diagnostic_audit$ppc_observation_flags,
    posterior_predictive_outcome_summary =
      diagnostic_audit$ppc_outcome_summary,
    initial_integrated_loo_model_summary = integrated_loo_audit$model_summary,
    initial_integrated_loo_outcome_summary =
      integrated_loo_audit$outcome_summary,
    initial_high_pareto_k_cases = integrated_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_correction_manifest = exact_audit$correction_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final updated achieved-intake checkpoint checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    review_decisions = "Eligibility decisions for updated-search studies.",
    preparation_checks = "Combined-data preparation invariants.",
    preparation_outcome_counts = "Analysis counts by primary outcome.",
    preparation_update_study_counts = "Updated-study contribution audit.",
    exposure_data_audit = "Achieved/assignment common-sample exposure audit.",
    analysis_data = "Combined common-sample rows with achieved intake.",
    model_specifications = "Final 104-model achieved specification registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit path, stage, size, and checksum.",
    sampling_diagnostic_summary = "Final model-level sampling diagnostics.",
    sampling_diagnostic_classification = "Final diagnostic decisions.",
    sampling_parameter_failures = "Parameter-level diagnostic failures.",
    posterior_predictive_observation_flags = "Observation-level PPC flags.",
    posterior_predictive_outcome_summary = "Outcome-level PPC summary.",
    initial_integrated_loo_model_summary = "Initial integrated PSIS summaries.",
    initial_integrated_loo_outcome_summary = "Initial LOO audit by outcome.",
    initial_high_pareto_k_cases = "The 19 cases selected for exact refitting.",
    exact_refit_diagnostic_summary = "Final diagnostics for 19 exact refits.",
    exact_correction_manifest = "Compact exact-correction paths and checksums.",
    final_loo_manifest = "Final exact-corrected LOO lineage.",
    final_loo_model_summary = "Final LOO summary for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    stacking_weights = "Final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    dose_response_predictions = "Achieved-intake model-averaged summaries.",
    dose_response_increments = "Increments from the 0.8 g/kg/day reference.",
    exposure_support_by_group = "Observed achieved-intake support by RT group.",
    common_exposure_support = "Common achieved-intake support by outcome.",
    pairwise_contrasts = "Model-averaged pairwise contrasts and probabilities.",
    pairwise_support = "Support limits for pairwise contrasts.",
    local_equivalence = "Local practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    figures_deferred = paste(
      "Achieved/assignment exposure comparison figures are deferred until",
      "the updated assignment-derived family is complete."
    )
  )
}

make_update_2026_achieved_pruning_plan <- function(
    accepted_fit_manifest,
    verification_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "achieved"
    )
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "achieved_update_2026"),
    "manifest.csv"
  )
  if (
    length(verification_files) != 36L ||
      any(!file.exists(verification_files))
  ) {
    stop("The updated achieved-intake verification bundle is incomplete.",
         call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "achieved"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated achieved archive.",
         call. = FALSE)
  }

  fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  if (any(!dir.exists(fit_roots))) {
    stop("An expected updated achieved fit directory is missing.",
         call. = FALSE)
  }
  fit_paths <- unlist(
    purrr::map(
      fit_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  )
  fit_paths <- normalizePath(fit_paths, winslash = "/", mustWork = TRUE)
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(fit_paths) != 108L || length(accepted_paths) != 104L ||
      anyDuplicated(fit_paths) || anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% fit_paths) ||
      sum(!fit_paths %in% accepted_paths) != 4L
  ) {
    stop("Updated achieved fit retention does not match its lineage.",
         call. = FALSE)
  }

  fit_manifest <- tibble::tibble(
    path = fit_paths,
    relative_path = substring(fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", fit_paths, fixed = TRUE) ~ "initial",
      grepl("/remedies/remedy_1/", fit_paths, fixed = TRUE) ~ "remedy_1",
      TRUE ~ "unknown"
    ),
    accepted_fit = fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(fit_paths)$size),
    lean_action = "delete frozen full fit archive"
  ) |>
    dplyr::arrange(archive_stage, relative_path)
  if (
    any(fit_manifest$archive_stage == "unknown") ||
      !identical(
        as.integer(table(factor(
          fit_manifest$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 4L)
      ) ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "initial") != 100L ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "remedy_1") != 4L
  ) {
    stop("The updated achieved archive stages are not as expected.",
         call. = FALSE)
  }

  compact_roots <- file.path(archive_root, c("loo", "reporting"))
  if (any(!dir.exists(compact_roots))) {
    stop("An expected compact achieved archive directory is missing.",
         call. = FALSE)
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(.x, recursive = TRUE, full.names = TRUE)
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[file.exists(compact_paths)]
  compact_manifest <- tibble::tibble(
    path = normalizePath(compact_paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(compact_paths)$size),
    action = "retain compact LOO/reporting archive"
  )
  if (nrow(compact_manifest) != 170L) {
    stop("The compact achieved-intake archive inventory is incomplete.",
         call. = FALSE)
  }

  summary <- tibble::tibble(
    mode = "lean",
    fit_files_deleted = nrow(fit_manifest),
    accepted_fit_files_deleted = sum(fit_manifest$accepted_fit),
    superseded_fit_files_deleted = sum(!fit_manifest$accepted_fit),
    fit_files_retained = 0L,
    gib_reclaimed = sum(fit_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained = sum(compact_manifest$bytes) / 1024^3,
    restoration = paste(
      "Reported results remain frozen; restoring full fits requires",
      "refitting the achieved-intake family."
    )
  )

  list(
    summary = summary,
    fit_file_manifest = fit_manifest,
    compact_archive_manifest = compact_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_achieved_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    fit_file_manifest = file.path(output_dir, "fit_file_manifest.csv"),
    compact_archive_manifest = file.path(
      output_dir,
      "compact_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(plan$fit_file_manifest, paths[["fit_file_manifest"]])
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_achieved_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "achieved", "pruning"
    )
) {
  required_confirmation <- "DELETE ALL UPDATED ACHIEVED FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "achieved"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the expected achieved archive.",
         call. = FALSE)
  }
  selected <- plan$fit_file_manifest
  if (
    nrow(selected) != 108L || sum(selected$accepted_fit) != 104L ||
      sum(!selected$accepted_fit) != 4L ||
      !identical(
        as.integer(table(factor(
          selected$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 4L)
      )
  ) {
    stop("The achieved pruning set is inconsistent with its audit.",
         call. = FALSE)
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  if (file.exists(receipt_path)) {
    stop("The updated achieved lean-pruning receipt already exists.",
         call. = FALSE)
  }
  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  permitted_subtree <- grepl(
    "/(candidates|remedies)/",
    selected_paths
  )
  if (
    any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted_subtree) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the permitted achieved fit subtree.",
         call. = FALSE)
  }
  if (any(!file.exists(plan$compact_archive_manifest$path))) {
    stop("A retained achieved compact artifact failed revalidation.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop("At least one achieved fit archive could not be deleted.",
         call. = FALSE)
  }
  if (
    any(!file.exists(plan$compact_archive_manifest$path)) ||
      !file.exists(plan$verification_manifest)
  ) {
    stop("A retained achieved artifact is unexpectedly missing.",
         call. = FALSE)
  }

  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, receipt_path)
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

# Combined-data assignment-derived-intake sensitivity refit ----

make_update_2026_assignment_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated assignment-derived base registry is incomplete.",
         call. = FALSE)
  }

  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  remedy_keys <- unname(purrr::map_chr(primary_remedy_specifications, key))
  if (
    anyDuplicated(primary_keys) || anyDuplicated(remedy_keys) ||
      !identical(remedy_keys, "function::10") ||
      !all(remedy_keys %in% primary_keys)
  ) {
    stop("The updated assignment-derived primary replacement is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  historical_remedies <- c(
    "hypertrophy::25", "strength::2", "strength::22"
  )
  specifications <- purrr::map(specifications, function(specification) {
    specification_key <- key(specification)
    historical_remedy <- specification_key %in% historical_remedies
    if (historical_remedy) {
      specification$adapt_delta <- max(
        as.numeric(specification$adapt_delta),
        0.999
      )
      specification$max_treedepth <- max(
        as.integer(specification$max_treedepth),
        15L
      )
    }
    specification$update_2026_sensitivity_base <- TRUE
    specification$update_2026_assignment_historical_remedy <-
      historical_remedy
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  historical_keys <- final_keys[purrr::map_lgl(
    specifications,
    ~ isTRUE(.x$update_2026_assignment_historical_remedy)
  )]
  function_ten <- specifications[[match("function::10", final_keys)]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !setequal(historical_keys, historical_remedies) ||
      !identical(as.numeric(function_ten$adapt_delta), 0.9995)
  ) {
    stop("The updated assignment-derived registry failed validation.",
         call. = FALSE)
  }

  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

update_2026_assignment_local_config <- function(specification, config) {
  local_config <- config
  historical_remedy <- isTRUE(
    specification$update_2026_assignment_historical_remedy
  )
  if (historical_remedy) {
    local_config$iter_warmup <- 2000L
    local_config$iter_sampling <- 4000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      15L
    )
  }
  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the updated assignment configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed +
      (if (historical_remedy) 1000L else 0L) +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

validate_update_2026_assignment_inputs <- function(
    specifications,
    exposure_data_by_outcome,
    priors_by_outcome,
    config
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.list(exposure_data_by_outcome) ||
      !all(c("achieved", "assignment", "audit") %in%
        names(exposure_data_by_outcome)) ||
      nrow(exposure_data_by_outcome$audit) != 4L ||
      !all(exposure_data_by_outcome$audit$identical_common_sample)
  ) {
    stop("The updated achieved/assignment common sample is invalid.",
         call. = FALSE)
  }

  checks <- purrr::map_dfr(specifications, function(specification) {
    local_config <- update_2026_assignment_local_config(
      specification,
      config
    )
    input <- prepare_exposure_full_input(
      specification = specification,
      exposure_data_by_outcome = exposure_data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      family = "assignment",
      config = local_config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      rows = nrow(input$data),
      studies = dplyr::n_distinct(input$data$study_no),
      arms = dplyr::n_distinct(input$data$arm_no),
      assignment_minimum = min(input$data$protein_intake_raw),
      assignment_maximum = max(input$data$protein_intake_raw),
      basis_dimension = if (is.null(input$basis_spec)) {
        0L
      } else {
        as.integer(input$basis_spec$basis_dimension)
      },
      historical_remedy_schedule = isTRUE(
        specification$update_2026_assignment_historical_remedy
      ),
      adapt_delta = max(
        input$specification$adapt_delta,
        local_config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        local_config$max_treedepth_floor
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  historical <- checks |>
    dplyr::filter(.data$historical_remedy_schedule)
  standard <- checks |>
    dplyr::filter(!.data$historical_remedy_schedule)
  counts <- table(factor(checks$outcome, levels = outcomes))
  sample_counts <- checks |>
    dplyr::distinct(
      outcome, rows, studies, arms,
      assignment_minimum, assignment_maximum
    )
  expected_historical <- c(
    "hypertrophy::25", "strength::2", "strength::22"
  )
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      nrow(sample_counts) != 4L || any(sample_counts$rows < 1L) ||
      any(!is.finite(sample_counts$assignment_minimum)) ||
      any(!is.finite(sample_counts$assignment_maximum)) ||
      !setequal(
        paste(historical$outcome, historical$model_index, sep = "::"),
        expected_historical
      ) ||
      any(historical$iter_warmup != 2000L) ||
      any(historical$iter_sampling != 4000L) ||
      any(historical$adapt_delta < 0.999) ||
      any(historical$max_treedepth < 15L) ||
      any(standard$iter_warmup != 1500L) ||
      any(standard$iter_sampling != 2000L) ||
      config$chains != 4L || config$parallel_chains != 4L
  ) {
    stop("The updated assignment-derived preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_assignment_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "assignment",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete assignment-derived preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_assignment_candidate_to_file <- function(
    specification,
    exposure_data_by_outcome,
    priors_by_outcome,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "assignment",
      "candidates"
    )
) {
  local_config <- update_2026_assignment_local_config(specification, config)
  input <- prepare_exposure_full_input(
    specification = specification,
    exposure_data_by_outcome = exposure_data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    family = "assignment",
    config = local_config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "assignment"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated assignment-derived fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated assignment-derived fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

compile_update_2026_assignment_diagnostic_audit <- function(diagnostics) {
  split_diagnostics <- split(
    diagnostics,
    purrr::map_chr(diagnostics, ~ .x$summary$outcome[[1L]])
  )
  required <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !setequal(names(split_diagnostics), required) ||
      any(vapply(split_diagnostics[required], length, integer(1L)) != 26L)
  ) {
    stop("Updated assignment-derived diagnostics are incomplete.",
         call. = FALSE)
  }
  compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
}

make_update_2026_assignment_remedy_config <- function(
    initial_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    seed_offset = 1000L
) {
  config <- initial_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(initial_config$seed + seed_offset)
  if (
    config$iter_warmup < initial_config$iter_warmup ||
      config$iter_sampling <= initial_config$iter_sampling ||
      config$adapt_delta_floor != initial_config$adapt_delta_floor ||
      config$max_treedepth_floor != initial_config$max_treedepth_floor
  ) {
    stop("The updated assignment remedy configuration is invalid.",
         call. = FALSE)
  }
  config
}

select_update_2026_assignment_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  expected_keys <- c(
    "whole_body::23",
    "hypertrophy::10", "hypertrophy::15", "hypertrophy::16",
    "hypertrophy::22",
    "strength::4", "strength::13", "strength::25",
    "function::12"
  )
  divergent_keys <- c("hypertrophy::10", "hypertrophy::22")
  treedepth_keys <- c("strength::4", "strength::13", "strength::25")
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(
      outcome,
      model_index = as.integer(model_index),
      key = paste(outcome, model_index, sep = "::"),
      divergent_transitions,
      maximum_treedepth_transitions
    )
  if (
    nrow(failed) != 9L || anyDuplicated(failed$key) ||
      !setequal(failed$key, expected_keys) ||
      !setequal(
        failed$key[failed$divergent_transitions > 0L],
        divergent_keys
      ) ||
      !setequal(
        failed$key[failed$maximum_treedepth_transitions > 0L],
        treedepth_keys
      ) ||
      any(failed$divergent_transitions[
        !failed$key %in% divergent_keys
      ] != 0L) ||
      any(failed$maximum_treedepth_transitions[
        !failed$key %in% treedepth_keys
      ] != 0L)
  ) {
    stop("The updated assignment remedy set no longer matches its audit.",
         call. = FALSE)
  }

  selected <- purrr::keep(specifications, function(specification) {
    paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    ) %in% failed$key
  }) |>
    purrr::map(function(specification) {
      specification_key <- paste(
        specification$outcome,
        specification$model_index,
        sep = "::"
      )
      if (isTRUE(
        specification$update_2026_assignment_historical_remedy
      )) {
        stop("An unexpected historical assignment model requires remedy.",
             call. = FALSE)
      }
      if (specification_key %in% divergent_keys) {
        specification$adapt_delta <- max(
          as.numeric(specification$adapt_delta),
          0.999
        )
        specification$update_2026_assignment_remedy <-
          "longer_sampling_and_adapt_delta_0.999"
      } else if (specification_key %in% treedepth_keys) {
        specification$max_treedepth <- max(
          as.integer(specification$max_treedepth),
          15L
        )
        specification$update_2026_assignment_remedy <-
          "longer_sampling_and_max_treedepth_15"
      } else {
        specification$update_2026_assignment_remedy <- "longer_sampling"
      }
      specification$update_2026_assignment_remedy_stage <- "remedy_1"
      specification
    })
  selected_keys <- unname(purrr::map_chr(
    selected,
    ~ paste(.x$outcome, .x$model_index, sep = "::")
  ))
  if (
    length(selected) != 9L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The updated assignment remedy specifications are incomplete.",
         call. = FALSE)
  }
  names(selected) <- paste0(
    purrr::map_chr(selected, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(selected, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  selected
}

update_2026_assignment_remedy_local_config <- function(
    specification,
    config
) {
  local_config <- config
  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in an assignment remedy configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

validate_update_2026_assignment_remedy_inputs <- function(
    specifications,
    exposure_data_by_outcome,
    priors_by_outcome,
    config
) {
  checks <- purrr::map_dfr(specifications, function(specification) {
    local_config <- update_2026_assignment_remedy_local_config(
      specification,
      config
    )
    input <- prepare_exposure_full_input(
      specification = specification,
      exposure_data_by_outcome = exposure_data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      family = "assignment",
      config = local_config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      remedy = specification$update_2026_assignment_remedy,
      adapt_delta = max(
        input$specification$adapt_delta,
        local_config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        local_config$max_treedepth_floor
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  expected_keys <- c(
    "whole_body::23",
    "hypertrophy::10", "hypertrophy::15", "hypertrophy::16",
    "hypertrophy::22",
    "strength::4", "strength::13", "strength::25",
    "function::12"
  )
  divergent_keys <- c("hypertrophy::10", "hypertrophy::22")
  treedepth_keys <- c("strength::4", "strength::13", "strength::25")
  divergence_remedies <- dplyr::filter(checks, keys %in% divergent_keys)
  treedepth_remedies <- dplyr::filter(checks, keys %in% treedepth_keys)
  sampling_remedies <- dplyr::filter(
    checks,
    !keys %in% c(divergent_keys, treedepth_keys)
  )
  if (
    nrow(checks) != 9L || anyDuplicated(keys) ||
      !setequal(keys, expected_keys) || !all(checks$valid) ||
      any(checks$iter_warmup != 2000L) ||
      any(checks$iter_sampling != 4000L) ||
      nrow(divergence_remedies) != 2L ||
      any(divergence_remedies$adapt_delta != 0.999) ||
      any(divergence_remedies$remedy !=
        "longer_sampling_and_adapt_delta_0.999") ||
      nrow(treedepth_remedies) != 3L ||
      any(treedepth_remedies$max_treedepth != 15L) ||
      any(treedepth_remedies$remedy !=
        "longer_sampling_and_max_treedepth_15") ||
      nrow(sampling_remedies) != 4L ||
      any(sampling_remedies$remedy != "longer_sampling")
  ) {
    stop("The updated assignment remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_assignment_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "assignment", "remedy_1",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 9L || !all(checks$valid)) {
    stop("Cannot write an incomplete assignment remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_assignment_remedy_to_file <- function(
    specification,
    exposure_data_by_outcome,
    priors_by_outcome,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "assignment",
      "remedies", "remedy_1"
    )
) {
  local_config <- update_2026_assignment_remedy_local_config(
    specification,
    config
  )
  input <- prepare_exposure_full_input(
    specification = specification,
    exposure_data_by_outcome = exposure_data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    family = "assignment",
    config = local_config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "assignment"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an assignment remedy fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("An assignment remedy fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

compile_update_2026_assignment_remedy_diagnostic_audit <- function(
    initial_diagnostics,
    remedy_diagnostics
) {
  diagnostic_key <- function(x) {
    paste(x$summary$outcome[[1L]], x$summary$model_index[[1L]], sep = "::")
  }
  if (length(initial_diagnostics) != 104L || length(remedy_diagnostics) != 9L) {
    stop("Updated assignment remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  initial_keys <- unname(purrr::map_chr(initial_diagnostics, diagnostic_key))
  remedy_keys <- unname(purrr::map_chr(remedy_diagnostics, diagnostic_key))
  expected_keys <- initial_diagnostics |>
    purrr::keep(~ !.x$summary$diagnostic_pass[[1L]]) |>
    purrr::map_chr(diagnostic_key)
  if (
    anyDuplicated(initial_keys) || anyDuplicated(remedy_keys) ||
      !setequal(remedy_keys, expected_keys) ||
      !all(remedy_keys %in% initial_keys)
  ) {
    stop("Updated assignment diagnostic replacement is invalid.",
         call. = FALSE)
  }
  final_diagnostics <- initial_diagnostics
  final_diagnostics[match(remedy_keys, initial_keys)] <- remedy_diagnostics
  compile_update_2026_assignment_diagnostic_audit(final_diagnostics)
}

make_update_2026_assignment_accepted_fit_manifest <- function(
    initial_fit_paths,
    remedy_fit_paths,
    final_diagnostic_audit
) {
  if (
    nrow(final_diagnostic_audit$diagnostic_summary) != 104L ||
      any(!final_diagnostic_audit$diagnostic_summary$diagnostic_pass)
  ) {
    stop("Updated assignment fits cannot be accepted before diagnostics pass.",
         call. = FALSE)
  }
  path_table <- function(paths, stage, expected) {
    paths <- unname(unlist(paths, use.names = FALSE))
    table <- tibble::tibble(
      outcome = basename(dirname(paths)),
      model_index = as.integer(stringr::str_match(
        basename(paths),
        "model_([0-9]+)\\.qs$"
      )[, 2]),
      stage = stage,
      fit_path = normalizePath(paths, winslash = "/", mustWork = TRUE)
    )
    if (
      nrow(table) != expected || anyNA(table$model_index) ||
        anyDuplicated(table[c("outcome", "model_index")]) ||
        any(file.info(table$fit_path)$size <= 0)
    ) {
      stop("An updated assignment fit stage is incomplete.",
           call. = FALSE)
    }
    table
  }

  initial <- path_table(initial_fit_paths, "initial", 104L)
  remedy <- path_table(remedy_fit_paths, "remedy_1", 9L)
  keys <- c("outcome", "model_index")
  remedy_keys <- paste(remedy$outcome, remedy$model_index, sep = "::")
  expected_remedy_keys <- c(
    "whole_body::23",
    "hypertrophy::10", "hypertrophy::15", "hypertrophy::16",
    "hypertrophy::22",
    "strength::4", "strength::13", "strength::25",
    "function::12"
  )
  if (!setequal(remedy_keys, expected_remedy_keys)) {
    stop("Updated assignment accepted replacements are invalid.",
         call. = FALSE)
  }

  accepted <- initial |>
    dplyr::anti_join(remedy[keys], by = keys) |>
    dplyr::bind_rows(remedy) |>
    dplyr::left_join(
      final_diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, model_name = model),
      by = keys,
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      family = "assignment",
      bytes = file.info(fit_path)$size,
      sha256 = purrr::map_chr(
        fit_path,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      )
    ) |>
    dplyr::select(
      family, outcome, model_index, model_name, stage,
      fit_path, bytes, sha256
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  stage_counts <- table(factor(
    accepted$stage,
    levels = c("initial", "remedy_1")
  ))
  if (
    nrow(accepted) != 104L ||
      anyDuplicated(accepted[c("outcome", "model_index")]) ||
      anyNA(accepted$model_name) || any(!file.exists(accepted$fit_path)) ||
      !identical(as.integer(stage_counts), c(95L, 9L)) ||
      any(nchar(accepted$sha256) != 64L)
  ) {
    stop("The updated assignment accepted-fit manifest is invalid.",
         call. = FALSE)
  }
  accepted
}

write_update_2026_assignment_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "assignment",
      "accepted_fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(manifest[c("outcome", "model_index")]) ||
      any(!file.exists(manifest$fit_path)) ||
      any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an invalid assignment accepted-fit manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

compute_update_2026_assignment_integrated_loo_to_file <- function(
    fit_record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "assignment", "loo",
      "integrated"
    )
) {
  result <- compute_exposure_full_file_integrated_loo(
    fit_record = fit_record,
    family = "assignment",
    loo_draws = loo_draws,
    k_threshold = k_threshold,
    cores = cores
  )
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(unname(as.character(result$family)), "assignment") ||
      !identical(result$remedy, "integrated_observation_level_psis")
  ) {
    stop("Cannot archive an invalid updated assignment LOO result.",
         call. = FALSE)
  }
  output_dir <- file.path(output_root, result$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf(
      "model_%02d_integrated_loo.qs",
      as.integer(result$model_index)
    )
  )
  write_atomic_qs(result, output_file)
}

write_update_2026_assignment_integrated_loo_audit <- function(
    audit,
    manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "loo_psis"
    )
) {
  if (
    !isTRUE(audit$integrated_psis_complete) || audit$models != 104L ||
      nrow(manifest) != 104L ||
      any(!file.exists(manifest$integrated_loo_path))
  ) {
    stop("The updated assignment integrated-LOO audit is incomplete.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

validate_update_2026_assignment_exact_inputs <- function(
    records,
    integrated_audit,
    config
) {
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    !is.list(records) || length(records) != expected || expected != 13L
  ) {
    stop("Expected 13 updated assignment exact-LOO records.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    integrated <- qs::qread(record$integrated_loo_path)
    observation <- as.integer(record$observation)
    if (
      !inherits(integrated$loo, "loo") ||
        !identical(unname(as.character(integrated$family)), "assignment") ||
        !identical(integrated$outcome, record$outcome) ||
        as.integer(integrated$model_index) != record$model_index ||
        observation < 1L ||
        observation > nrow(integrated$observation_manifest)
    ) {
      stop("An updated assignment exact record is invalid.",
           call. = FALSE)
    }
    held_out <- integrated$observation_manifest[observation, , drop = FALSE]
    observed_k <- as.numeric(
      loo::pareto_k_values(integrated$loo)[[observation]]
    )
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      observation = observation,
      study_no = as.character(held_out$study_no[[1L]]),
      arm_no = as.integer(held_out$arm_no[[1L]]),
      outcome_no = as.integer(held_out$outcome_no[[1L]]),
      pareto_k = observed_k,
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      original_fit_exists = file.exists(record$original_fit_path),
      integrated_loo_exists = file.exists(record$integrated_loo_path),
      valid = is.finite(observed_k) && observed_k > 0.7
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  if (
    nrow(checks) != expected || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(expected)) ||
      !identical(as.integer(outcome_counts), c(10L, 3L, 0L, 0L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$integrated_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L
  ) {
    stop("The updated assignment exact-LOO preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_assignment_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "assignment", "loo_exact",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 13L || !all(checks$valid)) {
    stop("Cannot write an incomplete assignment exact preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_assignment_exact_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "assignment", "loo",
    "exact", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

fit_and_score_update_2026_assignment_exact_to_file <- function(
    record,
    config,
    output_file = update_2026_assignment_exact_file_path(record)
) {
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !identical(unname(as.character(candidate$family)), "assignment") ||
      !identical(candidate$outcome, record$outcome) ||
      as.integer(candidate$model_index) != record$model_index ||
      !inherits(integrated$loo, "loo") ||
      !identical(unname(as.character(integrated$family)), "assignment") ||
      !identical(integrated$outcome, record$outcome) ||
      as.integer(integrated$model_index) != record$model_index ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      loo::pareto_k_values(integrated$loo)[[observation]] <= 0.7
  ) {
    stop("An updated assignment exact-LOO record failed validation.",
         call. = FALSE)
  }

  original_data <- candidate$model$data
  training_data <- original_data[-observation, , drop = FALSE]
  original_model <- candidate$model
  sample_prior <- attr(original_model$prior, "sample_prior")
  original_model$prior <- normalize_primary_brms_priors(
    original_model$prior
  )
  attr(original_model$prior, "sample_prior") <- sample_prior
  specification <- candidate$specification
  specification$adapt_delta <- max(
    specification$adapt_delta,
    config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth,
    config$max_treedepth_floor
  )

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
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  refit <- list(
    family = "assignment",
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    original_data = original_data,
    integrated_result = integrated,
    specification = specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  correction <- compute_sensitivity_mundlak_exact_loo_case(refit)
  if (
    as.integer(correction$case_id) != record$case_id ||
      !inherits(correction$loo, "loo") ||
      nrow(correction$exact_refit_summary) != 1L
  ) {
    stop("The updated assignment exact correction is invalid.",
         call. = FALSE)
  }
  write_atomic_qs(correction, output_file)
}

audit_update_2026_assignment_exact_files <- function(
    paths,
    integrated_audit
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    length(paths) != expected || expected != 13L || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0)
  ) {
    stop("The updated assignment exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  manifest <- purrr::map2_dfr(paths, corrections, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      observation =
        as.integer(result$exact_refit_summary$observation[[1L]]),
      correction_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys)
  ) {
    stop("The updated assignment exact audit is invalid.",
         call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = sum(diagnostic_summary$diagnostic_pass),
    exact_refits_requiring_remediation =
      sum(!diagnostic_summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass =
      all(diagnostic_summary$diagnostic_pass)
  )
}

write_update_2026_assignment_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "loo_exact"
    )
) {
  if (
    audit$exact_refits != 13L ||
      nrow(audit$diagnostic_summary) != 13L ||
      nrow(audit$correction_manifest) != 13L
  ) {
    stop("Cannot write an incomplete assignment exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

make_update_2026_assignment_exact_remedy_config <- function(
    exact_config,
    adapt_delta = 0.9995,
    seed_offset = 1000L
) {
  config <- exact_config
  config$adapt_delta_floor <- max(
    as.numeric(exact_config$adapt_delta_floor),
    as.numeric(adapt_delta)
  )
  config$seed <- as.integer(exact_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.9995 ||
      config$max_treedepth_floor != exact_config$max_treedepth_floor
  ) {
    stop("The assignment exact-remedy configuration is invalid.",
         call. = FALSE)
  }
  config
}

select_update_2026_assignment_exact_remedy_record <- function(
    records,
    exact_audit
) {
  failed <- exact_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass)
  if (
    nrow(failed) != 1L || failed$outcome[[1L]] != "whole_body" ||
      failed$model_index[[1L]] != 20L ||
      failed$observation[[1L]] != 131L ||
      failed$divergent_transitions[[1L]] != 1L ||
      failed$maximum_treedepth_transitions[[1L]] != 0L ||
      !isTRUE(failed$convergence_pass[[1L]])
  ) {
    stop("The assignment exact remedy no longer matches its audit.",
         call. = FALSE)
  }
  selected <- purrr::keep(records, function(record) {
    identical(record$outcome, "whole_body") &&
      as.integer(record$model_index) == 20L &&
      as.integer(record$observation) == 131L
  })
  if (
    length(selected) != 1L ||
      as.integer(selected[[1L]]$case_id) != 10L
  ) {
    stop("Could not uniquely select assignment exact case 10.",
         call. = FALSE)
  }
  names(selected) <- "case_010_whole_body_model_20_observation_0131"
  selected
}

validate_update_2026_assignment_exact_remedy_input <- function(
    record,
    config
) {
  if (!is.list(record) || length(record) != 1L) {
    stop("Expected one assignment exact-remedy record.", call. = FALSE)
  }
  record <- record[[1L]]
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  observed_k <- as.numeric(
    loo::pareto_k_values(integrated$loo)[[observation]]
  )
  effective_adapt_delta <- max(
    as.numeric(candidate$specification$adapt_delta),
    as.numeric(config$adapt_delta_floor)
  )
  effective_max_treedepth <- max(
    as.integer(candidate$specification$max_treedepth),
    as.integer(config$max_treedepth_floor)
  )
  check <- tibble::tibble(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    remedy = "adapt_delta_0.9995_after_exact_refit",
    adapt_delta = effective_adapt_delta,
    max_treedepth = effective_max_treedepth,
    iter_warmup = as.integer(config$iter_warmup),
    iter_sampling = as.integer(config$iter_sampling),
    pareto_k = observed_k,
    original_fit_exists = file.exists(record$original_fit_path),
    integrated_loo_exists = file.exists(record$integrated_loo_path),
    valid = TRUE
  )
  if (
    check$case_id[[1L]] != 10L ||
      check$outcome[[1L]] != "whole_body" ||
      check$model_index[[1L]] != 20L ||
      check$observation[[1L]] != 131L ||
      check$adapt_delta[[1L]] != 0.9995 ||
      check$max_treedepth[[1L]] != 15L ||
      check$iter_warmup[[1L]] != 2000L ||
      check$iter_sampling[[1L]] != 6000L ||
      !is.finite(observed_k) || observed_k <= 0.7 ||
      !check$original_fit_exists[[1L]] ||
      !check$integrated_loo_exists[[1L]]
  ) {
    stop("The assignment exact-remedy preflight is invalid.",
         call. = FALSE)
  }
  check
}

write_update_2026_assignment_exact_remedy_preflight <- function(
    check,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "assignment", "loo_exact",
      "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(check) != 1L || !isTRUE(check$valid[[1L]])) {
    stop("Cannot write an invalid assignment exact-remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(check, output_file)
  output_file
}

update_2026_assignment_exact_remedy_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "assignment", "loo",
    "exact_remedy", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

audit_update_2026_assignment_exact_remedy_file <- function(path, record) {
  path <- unname(unlist(path, use.names = FALSE))
  if (
    is.list(record) && is.null(record$case_id) && length(record) == 1L
  ) {
    record <- record[[1L]]
  }
  if (
    length(path) != 1L || !file.exists(path) || file.info(path)$size <= 0
  ) {
    stop("The assignment exact-remedy file is unavailable.",
         call. = FALSE)
  }
  correction <- qs::qread(path)
  summary <- correction$exact_refit_summary
  if (
    as.integer(correction$case_id) != 10L || nrow(summary) != 1L ||
      summary$outcome[[1L]] != "whole_body" ||
      summary$model_index[[1L]] != 20L ||
      summary$observation[[1L]] != 131L ||
      as.integer(record$case_id) != 10L
  ) {
    stop("The assignment exact-remedy result is invalid.",
         call. = FALSE)
  }
  manifest <- tibble::tibble(
    case_id = 10L,
    outcome = summary$outcome[[1L]],
    model_index = as.integer(summary$model_index[[1L]]),
    model_name = correction$model_name,
    observation = as.integer(summary$observation[[1L]]),
    correction_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = digest::digest(
      file = path,
      algo = "sha256",
      serialize = FALSE
    )
  )
  list(
    diagnostic_summary = summary,
    correction_manifest = manifest,
    exact_refits = 1L,
    exact_refits_passing = sum(summary$diagnostic_pass),
    exact_refits_requiring_remediation = sum(!summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

write_update_2026_assignment_exact_remedy_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "loo_exact",
      "remedy_1"
    )
) {
  if (
    audit$exact_refits != 1L ||
      nrow(audit$diagnostic_summary) != 1L ||
      nrow(audit$correction_manifest) != 1L
  ) {
    stop("Cannot write an invalid assignment exact-remedy audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

make_update_2026_assignment_final_exact_audit <- function(
    initial_audit,
    remedy_audit
) {
  if (
    initial_audit$exact_refits != 13L ||
      initial_audit$exact_refits_requiring_remediation != 1L ||
      remedy_audit$exact_refits != 1L ||
      !isTRUE(remedy_audit$all_exact_refit_diagnostics_pass)
  ) {
    stop("The assignment exact audits are not ready to combine.",
         call. = FALSE)
  }
  key <- function(data) {
    paste(data$outcome, data$model_index, data$observation, sep = "::")
  }
  remedy_key <- key(remedy_audit$diagnostic_summary)
  expected_key <- "whole_body::20::131"
  if (!identical(remedy_key, expected_key)) {
    stop("The assignment exact remedy has the wrong identity.",
         call. = FALSE)
  }
  diagnostic_summary <- initial_audit$diagnostic_summary |>
    dplyr::filter(key(dplyr::pick(dplyr::everything())) != remedy_key) |>
    dplyr::bind_rows(remedy_audit$diagnostic_summary) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  correction_manifest <- initial_audit$correction_manifest |>
    dplyr::filter(.data$case_id != 10L) |>
    dplyr::bind_rows(remedy_audit$correction_manifest) |>
    dplyr::arrange(case_id)
  if (
    nrow(diagnostic_summary) != 13L ||
      anyDuplicated(key(diagnostic_summary)) ||
      any(!diagnostic_summary$diagnostic_pass) ||
      nrow(correction_manifest) != 13L ||
      anyDuplicated(correction_manifest$case_id) ||
      !identical(correction_manifest$case_id, seq_len(13L)) ||
      any(!file.exists(correction_manifest$correction_path)) ||
      correction_manifest$correction_path[
        correction_manifest$case_id == 10L
      ] != remedy_audit$correction_manifest$correction_path[[1L]]
  ) {
    stop("The final assignment exact audit is invalid.", call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = correction_manifest,
    exact_refits = 13L,
    exact_refits_passing = 13L,
    exact_refits_requiring_remediation = 0L,
    all_exact_refit_diagnostics_pass = TRUE,
    remedied_case_id = 10L
  )
}

write_update_2026_assignment_final_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "loo_exact",
      "final"
    )
) {
  if (
    audit$exact_refits != 13L ||
      !isTRUE(audit$all_exact_refit_diagnostics_pass) ||
      nrow(audit$diagnostic_summary) != 13L ||
      nrow(audit$correction_manifest) != 13L
  ) {
    stop("Cannot write an incomplete final assignment exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

read_update_2026_assignment_exact_corrections <- function(manifest) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 13L ||
      any(!file.exists(manifest$correction_path))
  ) {
    stop("The assignment exact-correction manifest is invalid.",
         call. = FALSE)
  }
  corrections <- purrr::map(manifest$correction_path, qs::qread)
  corrections <- corrections[order(purrr::map_int(corrections, "case_id"))]
  if (
    !identical(
      unname(purrr::map_int(corrections, "case_id")),
      seq_len(13L)
    ) ||
      any(!purrr::map_lgl(corrections, ~ inherits(.x$loo, "loo")))
  ) {
    stop("The assignment exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections
}

make_update_2026_assignment_final_loo_audit <- function(
    integrated_manifest,
    exact_corrections,
    exact_audit
) {
  if (!isTRUE(exact_audit$all_exact_refit_diagnostics_pass)) {
    stop("Final assignment LOO requires passing exact diagnostics.",
         call. = FALSE)
  }
  final <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest = integrated_manifest,
    exact_results = exact_corrections,
    family = "assignment"
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- integrated_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "integrated_observation_level_psis_plus_exact",
        "integrated_observation_level_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 13L
  ) {
    stop("The final assignment LOO is not stacking-ready.",
         call. = FALSE)
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 13L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_assignment_final_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 13L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop("Cannot write an incomplete final assignment LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_assignment_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "assignment",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L ||
      anyDuplicated(paste(weights$outcome, weights$model_index)) ||
      length(outcome_sums) != 4L ||
      any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop("The updated assignment stacking weights are invalid.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_assignment_reporting_preflight <- function(records) {
  if (!is.list(records) || length(records) < 4L) {
    stop("The assignment prediction allocation is incomplete.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    tibble::tibble(
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      family = unname(as.character(record$family)),
      prediction_type = unname(as.character(record$prediction_type)),
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      protein_max = as.numeric(record$protein_max),
      grid_step = as.numeric(record$grid_step),
      fit_path = record$fit_path,
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path) &&
        identical(unname(as.character(record$family)), "assignment") &&
        identical(
          unname(as.character(record$prediction_type)),
          "exposure"
        ) &&
        is.finite(record$protein_max) &&
        record$protein_max > record$protein_min
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      allocated_weight = sum(weight),
      common_protein_max = dplyr::first(protein_max),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      any(!is.finite(outcome_audit$allocated_weight)) ||
      any(outcome_audit$allocated_weight <= 0) ||
      any(outcome_audit$allocated_weight > 1 + 1e-8) ||
      any(!is.finite(checks$weight)) || any(checks$weight < 0) ||
      any(!is.finite(outcome_audit$common_protein_max)) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop("The assignment reporting preflight failed.", call. = FALSE)
  }
  checks
}

write_update_2026_assignment_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "assignment", "reporting",
      "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(checks$allocated_draws, checks$outcome, sum)
  outcome_weights <- tapply(checks$weight, checks$outcome, sum)
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      length(outcome_weights) != 4L ||
      any(!is.finite(outcome_weights)) || any(outcome_weights <= 0) ||
      any(outcome_weights > 1 + 1e-8) || !all(checks$valid)
  ) {
    stop("Cannot write an invalid assignment reporting preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_assignment_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "assignment",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_assignment_prediction_draw_file <- function(
    record,
    data_by_outcome
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !identical(unname(as.character(record$family)), "assignment") ||
      !identical(unname(as.character(record$prediction_type)), "exposure")
  ) {
    stop("An updated assignment prediction record is invalid.",
         call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "assignment") ||
      !identical(candidate$outcome, unname(as.character(record$outcome))) ||
      as.integer(candidate$model_index) != as.integer(record$model_index) ||
      !identical(candidate$model_name, unname(as.character(record$model_name)))
  ) {
    stop("An updated assignment prediction fit does not match its record.",
         call. = FALSE)
  }
  prediction_data <- data_by_outcome[[record$outcome]]
  if (is.null(prediction_data) || !is.data.frame(prediction_data)) {
    stop("The updated assignment prediction data are unavailable.",
         call. = FALSE)
  }
  grid <- make_sensitivity_prediction_grid(
    fitted_candidate = candidate,
    data = prediction_data,
    prediction_type = "exposure",
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
    stop("Updated assignment prediction draws are incomplete.",
         call. = FALSE)
  }
  reduction <- list(
    family = "assignment",
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
      prediction_type = "exposure",
      reference_intake = NULL
    )
  )
  write_atomic_qs(
    reduction,
    update_2026_assignment_reporting_path(
      outcome = record$outcome,
      stage = "prediction_draws",
      model_index = record$model_index
    )
  )
}

write_update_2026_assignment_model_average_file <- function(record) {
  reductions <- purrr::map(
    record$prediction_manifest$prediction_path,
    read_sensitivity_full_prediction_draw_file,
    family = "assignment"
  )
  key_reference <- reductions[[1L]]$grid
  if (!all(purrr::map_lgl(
    reductions, ~ identical(.x$grid, key_reference)
  ))) {
    stop("Updated assignment prediction grids do not match within outcome.",
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
          file = path, algo = "sha256", serialize = FALSE
        ),
        within_model_draw = seq_len(reduction$draw_count)
      )
    }
  )
  if (
    nrow(draws) != record$n_draws ||
      nrow(draw_manifest) != record$n_draws || anyNA(draws)
  ) {
    stop("The updated assignment posterior mixture is incomplete.",
         call. = FALSE)
  }
  set.seed(record$averaging_seed)
  draw_order <- sample(seq_len(record$n_draws))
  draws <- draws[draw_order, , drop = FALSE]
  draw_manifest <- draw_manifest[draw_order, , drop = FALSE] |>
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
  write_atomic_qs(
    output,
    update_2026_assignment_reporting_path(
      outcome = record$outcome,
      stage = "model_average"
    )
  )
}

write_update_2026_assignment_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path, family = "assignment"
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (is.null(data) || !is.data.frame(data)) {
    stop("The updated assignment reporting data are unavailable.",
         call. = FALSE)
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
    family = "assignment",
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
  write_atomic_qs(
    report,
    update_2026_assignment_reporting_path(
      outcome = outcome,
      stage = "reporting"
    )
  )
}

write_update_2026_assignment_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "reporting"
    )
) {
  if (
    !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L
  ) {
    stop("The updated assignment reporting audit is incomplete.",
         call. = FALSE)
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "assignment"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir, "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir, "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir, "common_exposure_support.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_assignment_verification_checkpoint <- function(
    data_by_outcome,
    exposure_audit,
    source_manifest,
    review_decisions,
    preparation_audit,
    specifications,
    accepted_fit_manifest,
    diagnostic_audit,
    integrated_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "assignment"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(data) {
    data |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_audit$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_assignment_name <- paste0(
    name_alignment$specification_name,
    " [assignment]"
  )
  weight_sums <- stacking_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      .by = outcome
    )

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "all four assignment common-sample datasets are complete",
      "the combined-data preparation and exposure audits remain valid",
      "all five assignment registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist and have checksums",
      "all final model diagnostics pass",
      "integrated LOO contains 104 reductions and 13 audited high-k cases",
      "all 13 final exact refits pass and compact corrections exist",
      "final corrected LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four model-averaged reporting audits pass"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      identical(names(data_by_outcome), outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L) &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) && all(analysis_data$vi > 0) &&
        all(is.finite(analysis_data$protein_intake_raw)),
      all(preparation_audit$checks$passed) &&
        nrow(exposure_audit) == 4L &&
        all(exposure_audit$identical_common_sample) &&
        setequal(exposure_audit$outcome, outcomes),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_assignment_name) &&
        all(name_alignment$diagnostic_name == expected_assignment_name) &&
        all(name_alignment$loo_name == expected_assignment_name) &&
        all(name_alignment$weight_name == expected_assignment_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(accepted_fit_manifest$family == "assignment") &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      nrow(diagnostic_audit$diagnostic_summary) == 104L &&
        all(diagnostic_audit$diagnostic_summary$diagnostic_pass),
      integrated_loo_audit$models == 104L &&
        integrated_loo_audit$exact_refits_required == 13L &&
        nrow(integrated_loo_audit$high_k_cases) == 13L,
      exact_audit$exact_refits == 13L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        all(file.exists(exact_audit$correction_manifest$correction_path)),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(file.exists(prediction_manifest$prediction_path)),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated assignment verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    review_decisions = review_decisions,
    preparation_checks = preparation_audit$checks,
    preparation_outcome_counts = preparation_audit$outcome_counts,
    preparation_update_study_counts = preparation_audit$update_study_counts,
    exposure_data_audit = exposure_audit,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_audit$diagnostic_summary,
    sampling_diagnostic_classification =
      diagnostic_audit$definitive_classification,
    sampling_parameter_failures = diagnostic_audit$parameter_failures,
    posterior_predictive_observation_flags =
      diagnostic_audit$ppc_observation_flags,
    posterior_predictive_outcome_summary =
      diagnostic_audit$ppc_outcome_summary,
    initial_integrated_loo_model_summary = integrated_loo_audit$model_summary,
    initial_integrated_loo_outcome_summary =
      integrated_loo_audit$outcome_summary,
    initial_high_pareto_k_cases = integrated_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_correction_manifest = exact_audit$correction_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final updated assignment checkpoint checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    review_decisions = "Eligibility decisions for updated-search studies.",
    preparation_checks = "Combined-data preparation invariants.",
    preparation_outcome_counts = "Analysis counts by primary outcome.",
    preparation_update_study_counts = "Updated-study contribution audit.",
    exposure_data_audit = "Achieved/assignment common-sample exposure audit.",
    analysis_data = "Combined common-sample rows with assigned intake.",
    model_specifications = "Final 104-model assignment specification registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit path, stage, size, and checksum.",
    sampling_diagnostic_summary = "Final model-level sampling diagnostics.",
    sampling_diagnostic_classification = "Final diagnostic decisions.",
    sampling_parameter_failures = "Parameter-level diagnostic failures.",
    posterior_predictive_observation_flags = "Observation-level PPC flags.",
    posterior_predictive_outcome_summary = "Outcome-level PPC summary.",
    initial_integrated_loo_model_summary = "Initial integrated PSIS summaries.",
    initial_integrated_loo_outcome_summary = "Initial LOO audit by outcome.",
    initial_high_pareto_k_cases = "The 13 cases selected for exact refitting.",
    exact_refit_diagnostic_summary = "Final diagnostics for 13 exact refits.",
    exact_correction_manifest = "Compact exact-correction paths and checksums.",
    final_loo_manifest = "Final exact-corrected LOO lineage.",
    final_loo_model_summary = "Final LOO summary for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    stacking_weights = "Final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    dose_response_predictions = "Assignment model-averaged summaries.",
    dose_response_increments = "Increments from the 0.8 g/kg/day reference.",
    exposure_support_by_group = "Observed assignment support by RT group.",
    common_exposure_support = "Common assignment support by outcome.",
    pairwise_contrasts = "Model-averaged pairwise contrasts and probabilities.",
    pairwise_support = "Support limits for pairwise contrasts.",
    local_equivalence = "Local practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    comparison_figures_pending = paste(
      "Achieved/assignment exposure comparison figures are produced as",
      "the next cross-family checkpoint."
    )
  )
}

update_2026_exposure_comparison_path <- function(outcome) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity",
    "exposure_comparison", "comparison", outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(directory, "comparison.qs")
}

write_update_2026_exposure_comparison_file <- function(
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
    stop("The updated paired exposure model averages are misaligned.",
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
  write_atomic_qs(
    output,
    update_2026_exposure_comparison_path(record$outcome)
  )
}

write_update_2026_exposure_comparison_verification_bundle <- function(
    comparison_audit
) {
  if (!isTRUE(comparison_audit$all_comparison_checks_pass)) {
    stop("The updated exposure comparison is not ready for verification.",
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
    family = "exposure_comparison_update_2026",
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
      comparison_file_manifest = paste(
        "Checksums and dimensions of the updated outcome-specific",
        "exposure comparison files."
      ),
      achieved_dose_response = paste(
        "Updated model-averaged curves indexed by achieved protein intake."
      ),
      assignment_dose_response = paste(
        "Updated model-averaged curves indexed by baseline plus",
        "assignment-derived protein intake."
      ),
      assignment_minus_achieved = paste(
        "Posterior assignment-minus-achieved differences in updated",
        "absolute dose-response curves."
      ),
      achieved_increments = paste(
        "Updated achieved-intake increments from 0.8 g/kg/day."
      ),
      assignment_increments = paste(
        "Updated assignment-derived increments from 0.8 g/kg/day."
      ),
      assignment_minus_achieved_increment = paste(
        "Posterior assignment-minus-achieved differences in updated",
        "increments from 0.8 g/kg/day."
      )
    )
  )
}

make_update_2026_exposure_comparison_figure_data <- function(
    comparison_audit,
    reference_intake = 0.8
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !isTRUE(comparison_audit$all_comparison_checks_pass) ||
      nrow(comparison_audit$summary) != 4L
  ) {
    stop("Updated exposure comparison figures require a passing audit.",
         call. = FALSE)
  }
  comparisons <- purrr::map(
    comparison_audit$summary$comparison_path,
    read_exposure_full_comparison_file
  )
  bind_absolute <- function(component, analysis) {
    purrr::map_dfr(comparisons, function(comparison) {
      comparison[[component]] |>
        dplyr::transmute(
          outcome = comparison$outcome,
          resistance_exercise_code,
          resistance_exercise,
          protein_intake,
          estimate,
          conf.low,
          conf.high,
          analysis = analysis
        )
    })
  }
  bind_increment <- function(component, analysis) {
    purrr::map_dfr(comparisons, function(comparison) {
      comparison[[component]] |>
        dplyr::transmute(
          outcome = comparison$outcome,
          resistance_exercise_code,
          resistance_exercise,
          protein_intake,
          estimate = increment,
          conf.low = increment_conf_low,
          conf.high = increment_conf_high,
          analysis = analysis
        )
    })
  }
  bind_difference <- function(component) {
    purrr::map_dfr(comparisons, function(comparison) {
      dplyr::mutate(
        comparison[[component]],
        outcome = comparison$outcome,
        .before = 1L
      )
    })
  }
  absolute <- dplyr::bind_rows(
    bind_absolute("achieved_dose_response", "Achieved intake"),
    bind_absolute("assignment_dose_response", "Assignment-derived intake")
  ) |>
    dplyr::arrange(
      factor(outcome, levels = outcomes),
      resistance_exercise_code,
      factor(
        analysis,
        levels = c("Achieved intake", "Assignment-derived intake")
      ),
      protein_intake
    )
  increments <- dplyr::bind_rows(
    bind_increment("achieved_increments", "Achieved intake"),
    bind_increment(
      "assignment_increments",
      "Assignment-derived intake"
    )
  ) |>
    dplyr::arrange(
      factor(outcome, levels = outcomes),
      resistance_exercise_code,
      factor(
        analysis,
        levels = c("Achieved intake", "Assignment-derived intake")
      ),
      protein_intake
    )
  absolute_difference <- bind_difference("assignment_minus_achieved")
  increment_difference <- bind_difference(
    "assignment_minus_achieved_increment"
  )
  finite_intervals <- function(data) {
    nrow(data) > 0L &&
      all(is.finite(data$estimate)) &&
      all(is.finite(data$conf.low)) &&
      all(is.finite(data$conf.high)) &&
      all(data$conf.low <= data$estimate) &&
      all(data$estimate <= data$conf.high)
  }
  probability_columns <- c(
    "probability_assignment_greater",
    "probability_assignment_lower",
    "probability_meaningfully_greater",
    "probability_negligible_difference",
    "probability_meaningfully_lower"
  )
  probabilities <- rbind(
    as.matrix(absolute_difference[, probability_columns]),
    as.matrix(increment_difference[, probability_columns])
  )
  reference_rows <- increments |>
    dplyr::filter(abs(protein_intake - reference_intake) < 1e-10)
  cells <- increments |>
    dplyr::distinct(outcome, resistance_exercise_code, analysis)
  checks <- tibble::tibble(
    check = c(
      "absolute overlays contain both exposure definitions in eight panels",
      "increment overlays contain both exposure definitions in eight panels",
      "all overlay and direct-difference intervals are finite and ordered",
      "all reference increments are exactly zero",
      "all direct-comparison probabilities are valid",
      "direct-difference and overlay grids have matching dimensions"
    ),
    passed = c(
      setequal(unique(absolute$outcome), outcomes) &&
        nrow(absolute |>
          dplyr::distinct(
            outcome, resistance_exercise_code, analysis
          )) == 16L,
      nrow(cells) == 16L &&
        setequal(
          unique(increments$analysis),
          c("Achieved intake", "Assignment-derived intake")
        ),
      finite_intervals(absolute) && finite_intervals(increments) &&
        finite_intervals(absolute_difference) &&
        finite_intervals(increment_difference),
      nrow(reference_rows) == 16L &&
        all(abs(reference_rows$estimate) < 1e-12) &&
        all(abs(reference_rows$conf.low) < 1e-12) &&
        all(abs(reference_rows$conf.high) < 1e-12),
      all(is.finite(probabilities)) &&
        all(probabilities >= 0 & probabilities <= 1),
      nrow(absolute) == 2L * nrow(absolute_difference) &&
        nrow(increments) == 2L * nrow(increment_difference)
    )
  )
  if (any(is.na(checks$passed) | !checks$passed)) {
    stop(
      "Updated exposure comparison figure data failed validation: ",
      paste(
        checks$check[is.na(checks$passed) | !checks$passed],
        collapse = "; "
      ),
      call. = FALSE
    )
  }
  list(
    absolute = absolute,
    increments = increments,
    absolute_difference = absolute_difference,
    increment_difference = increment_difference,
    checks = checks,
    reference_intake = reference_intake
  )
}

write_update_2026_exposure_comparison_figures <- function(
    figure_data,
    data_dir = here::here(
      "results", "update_2026", "sensitivity", "exposure_comparison"
    ),
    figure_dir = here::here(
      "results", "update_2026", "figures", "sensitivity_comparisons"
    )
) {
  if (!all(figure_data$checks$passed)) {
    stop("Cannot plot failed updated exposure comparison data.",
         call. = FALSE)
  }
  outcome_names <- c(
    whole_body = "Whole-body lean soft tissue",
    hypertrophy = "Direct muscle mass / size",
    strength = "Strength",
    `function` = "Physical function"
  )
  prepare_plot_data <- function(data) {
    data |>
      dplyr::mutate(
        outcome = factor(
          outcome,
          levels = names(outcome_names),
          labels = unname(outcome_names)
        ),
        resistance_exercise = factor(
          resistance_exercise_code,
          levels = c(1, 0),
          labels = c(
            "Resistance exercise: Yes",
            "Resistance exercise: No"
          )
        ),
        panel = interaction(
          outcome,
          resistance_exercise,
          sep = "  |  ",
          lex.order = TRUE
        ),
        analysis = factor(
          analysis,
          levels = c("Achieved intake", "Assignment-derived intake")
        )
      )
  }
  absolute_plot_data <- prepare_plot_data(figure_data$absolute)
  increment_plot_data <- prepare_plot_data(figure_data$increments)
  colours <- c(
    `Achieved intake` = "#2B6CB0",
    `Assignment-derived intake` = "#C05621"
  )
  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(
        colour = "grey30",
        margin = ggplot2::margin(b = 8)
      ),
      plot.caption = ggplot2::element_text(colour = "grey35", hjust = 0),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "grey90", linewidth = 0.3
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey88", linewidth = 0.3
      ),
      legend.position = "top",
      legend.justification = "left",
      legend.key.width = grid::unit(1.3, "cm"),
      plot.margin = ggplot2::margin(12, 16, 12, 12)
    )
  make_plot <- function(data, increment = FALSE) {
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = protein_intake,
        y = estimate,
        colour = analysis,
        fill = analysis,
        group = analysis
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0, colour = "grey68", linewidth = 0.35
      ) +
      ggplot2::geom_vline(
        xintercept = figure_data$reference_intake,
        colour = "grey52",
        linewidth = 0.4,
        linetype = "dotted"
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = conf.low, ymax = conf.high),
        alpha = 0.10,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::facet_wrap(
        ggplot2::vars(panel), ncol = 2, scales = "free"
      ) +
      ggplot2::scale_colour_manual(values = colours) +
      ggplot2::scale_fill_manual(values = colours) +
      ggplot2::scale_x_continuous(
        breaks = scales::breaks_width(0.4),
        expand = ggplot2::expansion(mult = c(0.015, 0.025))
      ) +
      base_theme
    if (isTRUE(increment)) {
      plot + ggplot2::labs(
        title = "Achieved- and assignment-derived intake sensitivity curves",
        subtitle = paste(
          "Updated posterior change from 0.8 g/kg/day; shaded regions are",
          "95% credible intervals"
        ),
        x = "Protein intake (g/kg/day)",
        y = "Change in standardised outcome from 0.8 g/kg/day",
        colour = NULL,
        fill = NULL,
        caption = paste(
          "Both analyses use the same common sample and prediction grid.",
          "The dotted line marks the 0.8 g/kg/day reference."
        )
      )
    } else {
      plot + ggplot2::labs(
        title = "Achieved- and assignment-derived dose-response curves",
        subtitle = paste(
          "Updated common-sample analyses; shaded regions are 95%",
          "credible intervals"
        ),
        x = "Protein intake (g/kg/day)",
        y = "Standardised mean change",
        colour = NULL,
        fill = NULL,
        caption = paste(
          "Achieved intake uses observed end-of-intervention intake;",
          "assignment-derived intake uses baseline plus assigned change."
        )
      )
    }
  }
  absolute_plot <- make_plot(absolute_plot_data)
  increment_plot <- make_plot(increment_plot_data, increment = TRUE)

  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    absolute_overlay_data = file.path(
      data_dir, "absolute_overlay_data.csv.gz"
    ),
    increment_overlay_data = file.path(
      data_dir, "reference_increment_overlay_data.csv.gz"
    ),
    absolute_difference_data = file.path(
      data_dir, "assignment_minus_achieved_absolute.csv.gz"
    ),
    increment_difference_data = file.path(
      data_dir, "assignment_minus_achieved_increment.csv.gz"
    ),
    checks = file.path(data_dir, "checks.csv"),
    absolute_png = file.path(
      figure_dir, "achieved_assignment_absolute_overlay.png"
    ),
    absolute_pdf = file.path(
      figure_dir, "achieved_assignment_absolute_overlay.pdf"
    ),
    increment_png = file.path(
      figure_dir, "achieved_assignment_increment_overlay.png"
    ),
    increment_pdf = file.path(
      figure_dir, "achieved_assignment_increment_overlay.pdf"
    )
  )
  readr::write_csv(
    figure_data$absolute,
    paths[["absolute_overlay_data"]]
  )
  readr::write_csv(
    figure_data$increments,
    paths[["increment_overlay_data"]]
  )
  readr::write_csv(
    figure_data$absolute_difference,
    paths[["absolute_difference_data"]]
  )
  readr::write_csv(
    figure_data$increment_difference,
    paths[["increment_difference_data"]]
  )
  readr::write_csv(figure_data$checks, paths[["checks"]])
  ggplot2::ggsave(
    paths[["absolute_png"]], absolute_plot,
    width = 12, height = 14, units = "in", dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    paths[["absolute_pdf"]], absolute_plot,
    width = 12, height = 14, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    paths[["increment_png"]], increment_plot,
    width = 12, height = 14, units = "in", dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    paths[["increment_pdf"]], increment_plot,
    width = 12, height = 14, units = "in", bg = "white"
  )
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    stop("An updated exposure comparison figure artifact is missing.",
         call. = FALSE)
  }
  manifest <- tibble::tibble(
    artifact = names(paths),
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(paths)$size),
    sha256 = purrr::map_chr(
      paths,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    )
  )
  manifest_path <- file.path(data_dir, "manifest.csv")
  readr::write_csv(manifest, manifest_path)
  normalizePath(
    c(unname(paths), manifest_path),
    winslash = "/",
    mustWork = TRUE
  )
}

make_update_2026_assignment_pruning_plan <- function(
    accepted_fit_manifest,
    verification_files,
    comparison_audit,
    comparison_verification_files,
    figure_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "assignment"
    )
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "assignment_update_2026"),
    "manifest.csv"
  )
  comparison_verification_manifest <- file.path(
    here::here(
      "results", "verification", "exposure_comparison_update_2026"
    ),
    "manifest.csv"
  )
  if (
    length(verification_files) != 36L ||
      any(!file.exists(verification_files))
  ) {
    stop("The updated assignment verification bundle is incomplete.",
         call. = FALSE)
  }
  if (
    length(comparison_verification_files) != 9L ||
      any(!file.exists(comparison_verification_files))
  ) {
    stop("The updated exposure-comparison verification is incomplete.",
         call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)
  verify_update_2026_manifest_files(comparison_verification_manifest)
  figure_files <- unname(unlist(figure_files, use.names = FALSE))
  if (
    length(figure_files) != 10L || anyNA(figure_files) ||
      any(!file.exists(figure_files)) ||
      any(file.info(figure_files)$size <= 0)
  ) {
    stop("The updated exposure-comparison figure archive is incomplete.",
         call. = FALSE)
  }

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "assignment"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated assignment archive.",
         call. = FALSE)
  }

  fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  if (any(!dir.exists(fit_roots))) {
    stop("An expected updated assignment fit directory is missing.",
         call. = FALSE)
  }
  fit_paths <- unlist(
    purrr::map(
      fit_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  )
  fit_paths <- normalizePath(fit_paths, winslash = "/", mustWork = TRUE)
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(fit_paths) != 113L || length(accepted_paths) != 104L ||
      anyDuplicated(fit_paths) || anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% fit_paths) ||
      sum(!fit_paths %in% accepted_paths) != 9L
  ) {
    stop("Updated assignment fit retention does not match its lineage.",
         call. = FALSE)
  }
  fit_manifest <- tibble::tibble(
    path = fit_paths,
    relative_path = substring(fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", fit_paths, fixed = TRUE) ~ "initial",
      grepl("/remedies/remedy_1/", fit_paths, fixed = TRUE) ~ "remedy_1",
      TRUE ~ "unknown"
    ),
    accepted_fit = fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(fit_paths)$size),
    lean_action = "delete frozen full fit archive"
  ) |>
    dplyr::arrange(archive_stage, relative_path)
  if (
    any(fit_manifest$archive_stage == "unknown") ||
      !identical(
        as.integer(table(factor(
          fit_manifest$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 9L)
      ) ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "initial") != 95L ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "remedy_1") != 9L
  ) {
    stop("The updated assignment archive stages are not as expected.",
         call. = FALSE)
  }

  compact_roots <- file.path(archive_root, c("loo", "reporting"))
  if (any(!dir.exists(compact_roots))) {
    stop("An expected compact assignment archive directory is missing.",
         call. = FALSE)
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(.x, recursive = TRUE, full.names = TRUE)
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[file.exists(compact_paths)]
  compact_manifest <- tibble::tibble(
    path = normalizePath(compact_paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(compact_paths)$size),
    action = "retain compact assignment LOO/reporting archive"
  )
  if (nrow(compact_manifest) != 160L) {
    stop("The compact assignment archive inventory is incomplete.",
         call. = FALSE)
  }

  comparison_manifest <- comparison_audit$summary |>
    dplyr::transmute(
      outcome,
      path = normalizePath(
        comparison_path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = as.numeric(bytes),
      sha256,
      action = "retain compact achieved/assignment comparison"
    )
  if (
    !isTRUE(comparison_audit$all_comparison_checks_pass) ||
      nrow(comparison_manifest) != 4L ||
      any(!file.exists(comparison_manifest$path)) ||
      any(nchar(comparison_manifest$sha256) != 64L)
  ) {
    stop("The compact exposure-comparison archive is incomplete.",
         call. = FALSE)
  }
  figure_manifest <- tibble::tibble(
    path = normalizePath(figure_files, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(figure_files)$size),
    sha256 = purrr::map_chr(
      figure_files,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    ),
    action = "retain updated exposure-comparison artifact"
  )

  summary <- tibble::tibble(
    mode = "lean",
    fit_files_deleted = nrow(fit_manifest),
    accepted_fit_files_deleted = sum(fit_manifest$accepted_fit),
    superseded_fit_files_deleted = sum(!fit_manifest$accepted_fit),
    fit_files_retained = 0L,
    gib_reclaimed = sum(fit_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained = sum(compact_manifest$bytes) / 1024^3,
    comparison_archive_files_retained = nrow(comparison_manifest),
    comparison_verification_files_retained =
      length(comparison_verification_files),
    figure_files_retained = nrow(figure_manifest),
    figure_archive_mib_retained = sum(figure_manifest$bytes) / 1024^2,
    restoration = paste(
      "Reported results remain frozen; restoring full fits requires",
      "refitting the assignment-derived family."
    )
  )

  list(
    summary = summary,
    fit_file_manifest = fit_manifest,
    compact_archive_manifest = compact_manifest,
    comparison_archive_manifest = comparison_manifest,
    figure_archive_manifest = figure_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    comparison_verification_manifest = normalizePath(
      comparison_verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_assignment_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    fit_file_manifest = file.path(output_dir, "fit_file_manifest.csv"),
    compact_archive_manifest = file.path(
      output_dir, "compact_archive_manifest.csv"
    ),
    comparison_archive_manifest = file.path(
      output_dir, "comparison_archive_manifest.csv"
    ),
    figure_archive_manifest = file.path(
      output_dir, "figure_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(plan$fit_file_manifest, paths[["fit_file_manifest"]])
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  readr::write_csv(
    plan$comparison_archive_manifest,
    paths[["comparison_archive_manifest"]]
  )
  readr::write_csv(
    plan$figure_archive_manifest,
    paths[["figure_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_assignment_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "assignment", "pruning"
    )
) {
  required_confirmation <- "DELETE ALL UPDATED ASSIGNMENT FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)
  verify_update_2026_manifest_files(
    plan$comparison_verification_manifest
  )

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "assignment"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the assignment archive.",
         call. = FALSE)
  }
  selected <- plan$fit_file_manifest
  if (
    nrow(selected) != 113L || sum(selected$accepted_fit) != 104L ||
      sum(!selected$accepted_fit) != 9L ||
      !identical(
        as.integer(table(factor(
          selected$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 9L)
      )
  ) {
    stop("The assignment pruning set is inconsistent with its audit.",
         call. = FALSE)
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  if (file.exists(receipt_path)) {
    stop("The updated assignment lean-pruning receipt already exists.",
         call. = FALSE)
  }
  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  permitted_subtree <- grepl(
    "/(candidates|remedies)/",
    selected_paths
  )
  if (
    any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted_subtree) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the assignment fit subtree.",
         call. = FALSE)
  }
  hashed_manifests <- list(
    plan$comparison_archive_manifest,
    plan$figure_archive_manifest
  )
  hashes_valid <- purrr::map_lgl(hashed_manifests, function(manifest) {
    all(file.exists(manifest$path)) &&
      all(purrr::map_chr(
        manifest$path,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ) == manifest$sha256)
  })
  if (
    any(!file.exists(plan$compact_archive_manifest$path)) ||
      !all(hashes_valid)
  ) {
    stop("A retained assignment or comparison artifact failed revalidation.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop("At least one assignment fit archive could not be deleted.",
         call. = FALSE)
  }
  if (
    any(!file.exists(plan$compact_archive_manifest$path)) ||
      any(!file.exists(plan$comparison_archive_manifest$path)) ||
      any(!file.exists(plan$figure_archive_manifest$path)) ||
      !file.exists(plan$verification_manifest) ||
      !file.exists(plan$comparison_verification_manifest)
  ) {
    stop("A retained assignment artifact is unexpectedly missing.",
         call. = FALSE)
  }

  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, receipt_path)
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

# Combined-data errors-in-variables sensitivity refit ----

make_update_2026_eiv_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated EIV base specification registry is incomplete.",
         call. = FALSE)
  }

  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  remedy_keys <- unname(purrr::map_chr(primary_remedy_specifications, key))
  if (
    anyDuplicated(primary_keys) || anyDuplicated(remedy_keys) ||
      !identical(remedy_keys, "function::10") ||
      !all(remedy_keys %in% primary_keys)
  ) {
    stop("The updated EIV primary replacement is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  longer_sampling <- c(
    "whole_body::2", "whole_body::5", "whole_body::26"
  )
  first_escalation <- c(
    "whole_body::23", "hypertrophy::23",
    "strength::7", "strength::13", "strength::15"
  )
  second_escalation <- "whole_body::14"
  historical_keys <- c(
    longer_sampling, first_escalation, second_escalation
  )

  specifications <- purrr::map(specifications, function(specification) {
    specification_key <- key(specification)
    schedule <- dplyr::case_when(
      specification_key %in% longer_sampling ~ "longer_sampling",
      specification_key %in% first_escalation ~ "remedy_1",
      specification_key %in% second_escalation ~ "remedy_2",
      TRUE ~ "base"
    )
    specification$update_2026_sensitivity_base <- TRUE
    specification$update_2026_eiv_historical_schedule <- schedule
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  scheduled_keys <- final_keys[purrr::map_chr(
    specifications,
    "update_2026_eiv_historical_schedule"
  ) != "base"]
  function_ten <- specifications[[match("function::10", final_keys)]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !setequal(scheduled_keys, historical_keys) ||
      !identical(as.numeric(function_ten$adapt_delta), 0.9995)
  ) {
    stop("The updated EIV registry failed validation.", call. = FALSE)
  }

  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

update_2026_eiv_local_config <- function(specification, config) {
  local_config <- config
  schedule <- specification$update_2026_eiv_historical_schedule
  if (
    length(schedule) != 1L || is.na(schedule) ||
      !schedule %in% c("base", "longer_sampling", "remedy_1", "remedy_2")
  ) {
    stop("The updated EIV sampling schedule is invalid.", call. = FALSE)
  }

  schedule_offset <- switch(
    schedule,
    base = 0L,
    longer_sampling = 1000L,
    remedy_1 = 2000L,
    remedy_2 = 3000L
  )
  if (schedule == "longer_sampling" || schedule == "remedy_1") {
    local_config$iter_warmup <- 2000L
    local_config$iter_sampling <- 4000L
  } else if (schedule == "remedy_2") {
    local_config$iter_warmup <- 3000L
    local_config$iter_sampling <- 8000L
  }
  if (schedule %in% c("remedy_1", "remedy_2")) {
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      15L
    )
  }

  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the updated EIV configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed + schedule_offset +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

validate_update_2026_eiv_inputs <- function(inputs, config) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    local_config <- update_2026_eiv_local_config(specification, config)
    stan_data <- input$stan_data
    protein_median <- min(input$grid$protein_intake_raw) -
      min(input$grid$protein_intake)
    effective_adapt_delta <- max(
      as.numeric(specification$adapt_delta),
      as.numeric(local_config$adapt_delta_floor)
    )
    effective_max_treedepth <- max(
      as.integer(specification$max_treedepth),
      as.integer(local_config$max_treedepth_floor)
    )
    valid <-
      identical(input$family, "eiv") && is.list(stan_data) &&
      stan_data$N > 0L && stan_data$N_arm > 0L &&
      stan_data$N_study > 0L && stan_data$N_outcome > 0L &&
      length(stan_data$y) == stan_data$N &&
      length(stan_data$arm) == stan_data$N &&
      length(stan_data$x_obs) == stan_data$N_arm &&
      length(stan_data$x_se) == stan_data$N_arm &&
      all(is.finite(stan_data$x_obs)) &&
      all(is.finite(stan_data$x_se)) && all(stan_data$x_se > 0) &&
      all(is.finite(input$grid$protein_intake)) &&
      is.finite(effective_adapt_delta) &&
      effective_adapt_delta > 0 && effective_adapt_delta < 1 &&
      effective_max_treedepth >= 1L

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      rows = as.integer(stan_data$N),
      studies = as.integer(stan_data$N_study),
      arms = as.integer(stan_data$N_arm),
      outcome_rows = as.integer(stan_data$N_outcome),
      protein_minimum = min(stan_data$x_obs) + protein_median,
      protein_maximum = max(stan_data$x_obs) + protein_median,
      minimum_protein_se = min(stan_data$x_se),
      maximum_protein_se = max(stan_data$x_se),
      historical_schedule =
        specification$update_2026_eiv_historical_schedule,
      adapt_delta = effective_adapt_delta,
      max_treedepth = effective_max_treedepth,
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      chains = as.integer(local_config$chains),
      parallel_chains = as.integer(local_config$parallel_chains),
      valid = valid
    )
  }) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  counts <- table(factor(checks$outcome, levels = outcomes))
  schedule_counts <- table(factor(
    checks$historical_schedule,
    levels = c("base", "longer_sampling", "remedy_1", "remedy_2")
  ))
  sample_counts <- checks |>
    dplyr::distinct(
      outcome, rows, studies, arms, outcome_rows,
      protein_minimum, protein_maximum,
      minimum_protein_se, maximum_protein_se
    )
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      !identical(as.integer(schedule_counts), c(95L, 3L, 5L, 1L)) ||
      nrow(sample_counts) != 4L ||
      any(sample_counts$rows < 1L) || any(sample_counts$arms < 1L) ||
      any(sample_counts$minimum_protein_se <= 0) ||
      any(!is.finite(sample_counts$maximum_protein_se)) ||
      any(checks$chains != 4L) || any(checks$parallel_chains != 4L)
  ) {
    stop("The updated EIV input preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_eiv_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "eiv",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated EIV preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_eiv_candidate_to_file <- function(
    input,
    model,
    config,
    stage = "initial",
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "eiv",
      "candidates"
    )
) {
  if (is_sensitivity_pilot_skip(model)) {
    stop("The updated EIV Stan model was not compiled.", call. = FALSE)
  }
  specification <- input$specification
  local_config <- update_2026_eiv_local_config(specification, config)
  effective_adapt_delta <- max(
    as.numeric(specification$adapt_delta),
    as.numeric(local_config$adapt_delta_floor)
  )
  effective_max_treedepth <- max(
    as.integer(specification$max_treedepth),
    as.integer(local_config$max_treedepth_floor)
  )
  specification$adapt_delta <- effective_adapt_delta
  specification$max_treedepth <- effective_max_treedepth
  specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)
  input_sha256 <- digest::digest(
    list(
      specification = specification,
      stan_data = input$stan_data,
      grid = input$grid
    ),
    algo = "sha256"
  )
  stan_file <- model$stan_file()
  if (length(stan_file) != 1L || is.na(stan_file) || !file.exists(stan_file)) {
    stop("The updated EIV Stan source is unavailable.", call. = FALSE)
  }
  stan_sha256 <- digest::digest(
    file = stan_file,
    algo = "sha256",
    serialize = FALSE
  )

  candidate_dir <- file.path(
    output_root,
    specification$outcome,
    sprintf("model_%02d", as.integer(specification$model_index))
  )
  output_file <- file.path(candidate_dir, "fit.qs")
  cmdstan_dir <- file.path(candidate_dir, "cmdstan")
  dir.create(cmdstan_dir, recursive = TRUE, showWarnings = FALSE)

  existing_candidate <- if (file.exists(output_file)) {
    tryCatch(read_eiv_fit_file(output_file), error = function(e) NULL)
  } else {
    NULL
  }
  if (
    !is.null(existing_candidate) &&
      identical(existing_candidate$outcome, specification$outcome) &&
      identical(
        as.integer(existing_candidate$model_index),
        as.integer(specification$model_index)
      ) &&
      identical(
        existing_candidate$specification$update_2026_eiv_historical_schedule,
        specification$update_2026_eiv_historical_schedule
      ) &&
      identical(
        as.integer(
          existing_candidate$specification$update_2026_iter_sampling
        ),
        as.integer(local_config$iter_sampling)
      ) &&
      identical(existing_candidate$input_sha256, input_sha256) &&
      identical(existing_candidate$stan_sha256, stan_sha256)
  ) {
    return(normalizePath(output_file, winslash = "/", mustWork = TRUE))
  }

  fit <- model$sample(
    data = input$stan_data,
    chains = local_config$chains,
    parallel_chains = local_config$parallel_chains,
    iter_warmup = local_config$iter_warmup,
    iter_sampling = local_config$iter_sampling,
    adapt_delta = effective_adapt_delta,
    max_treedepth = effective_max_treedepth,
    seed = local_config$seed,
    refresh = 100,
    output_dir = cmdstan_dir
  )
  candidate <- list(
    family = "eiv",
    archive_stage = stage,
    outcome = specification$outcome,
    model_index = as.integer(specification$model_index),
    model_name = paste0(specification$model_name, " [EIV Stan]"),
    specification = specification,
    grid = input$grid,
    output_dir = normalizePath(
      cmdstan_dir,
      winslash = "/",
      mustWork = TRUE
    ),
    input_sha256 = input_sha256,
    stan_sha256 = stan_sha256,
    fit = fit
  )
  write_atomic_qs(candidate, output_file)
}

make_update_2026_eiv_remedy_config <- function(
    initial_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    seed_offset = 5000L
) {
  config <- initial_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(initial_config$seed + seed_offset)
  if (
    config$iter_warmup < initial_config$iter_warmup ||
      config$iter_sampling <= initial_config$iter_sampling ||
      config$adapt_delta_floor != initial_config$adapt_delta_floor ||
      config$max_treedepth_floor != initial_config$max_treedepth_floor
  ) {
    stop("The updated EIV remedy configuration is invalid.",
         call. = FALSE)
  }
  config
}

select_update_2026_eiv_remedy_inputs <- function(inputs, diagnostic_audit) {
  expected_keys <- c(
    "whole_body::3", "whole_body::9", "whole_body::11",
    "whole_body::13", "whole_body::21", "hypertrophy::6",
    "strength::3", "strength::8", "strength::9", "function::5"
  )
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(
      outcome,
      model_index = as.integer(model_index),
      key = paste(outcome, model_index, sep = "::"),
      divergent_transitions,
      maximum_treedepth_transitions,
      minimum_ebfmi,
      maximum_rhat,
      minimum_bulk_ess,
      minimum_tail_ess,
      maximum_mcse_sd_ratio
    )
  if (
    !identical(diagnostic_audit$family, "eiv") ||
      nrow(failed) != 10L || anyDuplicated(failed$key) ||
      !setequal(failed$key, expected_keys) ||
      any(failed$divergent_transitions != 0L) ||
      any(failed$maximum_treedepth_transitions != 0L) ||
      any(failed$minimum_ebfmi < 0.2) ||
      any(
        failed$maximum_rhat <= 1.01 &
          failed$minimum_bulk_ess >= 400 &
          failed$minimum_tail_ess >= 400 &
          failed$maximum_mcse_sd_ratio <= 0.10
      )
  ) {
    stop("The updated EIV remedy set no longer matches its audit.",
         call. = FALSE)
  }

  selected <- purrr::keep(inputs, function(input) {
    paste(
      input$specification$outcome,
      input$specification$model_index,
      sep = "::"
    ) %in% failed$key
  }) |>
    purrr::map(function(input) {
      if (
        !identical(
          input$specification$update_2026_eiv_historical_schedule,
          "base"
        )
      ) {
        stop("A historically escalated EIV fit unexpectedly needs remedy.",
             call. = FALSE)
      }
      input$specification$update_2026_eiv_remedy <- "longer_sampling"
      input$specification$update_2026_eiv_remedy_stage <- "remedy_1"
      input
    })
  selected_keys <- unname(purrr::map_chr(
    selected,
    ~ paste(
      .x$specification$outcome,
      .x$specification$model_index,
      sep = "::"
    )
  ))
  if (
    length(selected) != 10L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The updated EIV remedy inputs are incomplete.", call. = FALSE)
  }
  names(selected) <- paste0(
    purrr::map_chr(selected, ~ .x$specification$outcome),
    "_model_",
    stringr::str_pad(
      purrr::map_int(selected, ~ .x$specification$model_index),
      width = 2L,
      pad = "0"
    )
  )
  selected
}

validate_update_2026_eiv_remedy_inputs <- function(
    inputs,
    initial_manifest,
    config
) {
  expected_keys <- c(
    "whole_body::3", "whole_body::9", "whole_body::11",
    "whole_body::13", "whole_body::21", "hypertrophy::6",
    "strength::3", "strength::8", "strength::9", "function::5"
  )
  initial_keys <- paste(
    initial_manifest$outcome,
    initial_manifest$model_index,
    sep = "::"
  )
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    key <- paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    )
    local_config <- update_2026_eiv_local_config(specification, config)
    initial_index <- match(key, initial_keys)
    initial_fit_exists <-
      !is.na(initial_index) &&
      file.exists(initial_manifest$fit_path[[initial_index]])
    valid <-
      identical(input$family, "eiv") &&
      identical(specification$update_2026_eiv_remedy, "longer_sampling") &&
      identical(specification$update_2026_eiv_remedy_stage, "remedy_1") &&
      identical(
        specification$update_2026_eiv_historical_schedule,
        "base"
      ) &&
      input$stan_data$N > 0L && input$stan_data$N_arm > 0L &&
      all(is.finite(input$stan_data$x_se)) &&
      all(input$stan_data$x_se > 0) &&
      initial_fit_exists &&
      local_config$iter_warmup == 2000L &&
      local_config$iter_sampling == 4000L

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      remedy = specification$update_2026_eiv_remedy,
      adapt_delta = max(
        as.numeric(specification$adapt_delta),
        as.numeric(local_config$adapt_delta_floor)
      ),
      max_treedepth = max(
        as.integer(specification$max_treedepth),
        as.integer(local_config$max_treedepth_floor)
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      initial_fit_exists = initial_fit_exists,
      valid = valid
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  if (
    nrow(initial_manifest) != 104L || anyDuplicated(initial_keys) ||
      nrow(checks) != 10L || anyDuplicated(keys) ||
      !setequal(keys, expected_keys) || !all(checks$valid) ||
      any(checks$iter_warmup != 2000L) ||
      any(checks$iter_sampling != 4000L)
  ) {
    stop("The updated EIV remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_eiv_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "eiv",
      "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(checks) != 10L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated EIV remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

write_update_2026_eiv_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "eiv",
      "accepted_fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$fit_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an incomplete updated EIV accepted manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

compute_update_2026_eiv_integrated_loo_to_file <- function(
    fit_record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "eiv", "loo",
      "integrated"
    )
) {
  result <- compute_eiv_full_file_integrated_loo(
    record = fit_record,
    loo_draws = loo_draws,
    k_threshold = k_threshold,
    cores = cores
  )
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(unname(as.character(result$family)), "eiv") ||
      !identical(
        result$remedy,
        "outcome_effects_integrated_latent_exposure_psis"
      ) ||
      !is.list(result$input) ||
      result$model_index != fit_record$model_index ||
      !identical(result$outcome, fit_record$outcome)
  ) {
    stop("Cannot archive an invalid updated EIV LOO result.",
         call. = FALSE)
  }
  output_dir <- file.path(output_root, result$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf(
      "model_%02d_integrated_loo.qs",
      as.integer(result$model_index)
    )
  )
  write_atomic_qs(result, output_file)
}

write_update_2026_eiv_integrated_loo_audit <- function(
    audit,
    manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "eiv", "loo_psis"
    )
) {
  if (
    !identical(audit$family, "eiv") ||
      !isTRUE(audit$integrated_psis_complete) ||
      audit$models != 104L || nrow(manifest) != 104L ||
      any(!file.exists(manifest$integrated_loo_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("The updated EIV integrated-LOO audit is incomplete.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

make_update_2026_eiv_exact_config <- function(
    initial_config,
    iter_warmup = 2000L,
    iter_sampling = 6000L,
    adapt_delta = 0.999,
    max_treedepth = 15L,
    seed_offset = 9000L
) {
  config <- initial_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$adapt_delta_floor <- max(
    as.numeric(initial_config$adapt_delta_floor),
    as.numeric(adapt_delta)
  )
  config$max_treedepth_floor <- max(
    as.integer(initial_config$max_treedepth_floor),
    as.integer(max_treedepth)
  )
  config$seed <- as.integer(initial_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.999 ||
      config$max_treedepth_floor != 15L
  ) {
    stop("The updated EIV exact-LOO configuration is invalid.",
         call. = FALSE)
  }
  config
}

validate_update_2026_eiv_exact_inputs <- function(
    records,
    integrated_audit,
    config
) {
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (!is.list(records) || length(records) != expected || expected != 9L) {
    stop("Expected nine updated EIV exact-LOO records.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    integrated_exists <- file.exists(record$integrated_loo_path)
    original_exists <- file.exists(record$original_fit_path)
    integrated <- if (integrated_exists) {
      qs::qread(record$integrated_loo_path)
    } else {
      NULL
    }
    observation <- as.integer(record$observation)
    pareto_values <- if (
      is.list(integrated) && inherits(integrated$loo, "loo")
    ) {
      loo::pareto_k_values(integrated$loo)
    } else {
      numeric()
    }
    input <- if (is.list(integrated)) integrated$input else NULL
    observation_manifest <- if (
      is.list(integrated) && is.data.frame(integrated$observation_manifest) &&
        observation >= 1L &&
        observation <= nrow(integrated$observation_manifest)
    ) {
      integrated$observation_manifest[observation, , drop = FALSE]
    } else {
      NULL
    }
    deletion_valid <- FALSE
    if (
      is.list(input) && is.list(input$stan_data) &&
        observation >= 1L && observation <= input$stan_data$N
    ) {
      training_data <- make_eiv_leave_one_out_stan_data(
        input$stan_data,
        observation
      )
      deletion_valid <-
        training_data$N == input$stan_data$N - 1L &&
        training_data$N_arm == input$stan_data$N_arm &&
        all(is.finite(training_data$x_se)) &&
        all(training_data$x_se > 0)
    }
    valid <-
      original_exists && integrated_exists &&
      nchar(record$original_fit_sha256) == 64L &&
      nchar(record$integrated_loo_sha256) == 64L &&
      is.list(integrated) &&
      identical(unname(as.character(integrated$family)), "eiv") &&
      identical(integrated$outcome, record$outcome) &&
      as.integer(integrated$model_index) == record$model_index &&
      !is.null(observation_manifest) &&
      observation >= 1L && observation <= length(pareto_values) &&
      is.finite(pareto_values[[observation]]) &&
      pareto_values[[observation]] > 0.7 &&
      deletion_valid

    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      observation = observation,
      study_no = if (is.null(observation_manifest)) {
        NA_character_
      } else {
        as.character(observation_manifest$study_no[[1L]])
      },
      arm_no = if (is.null(observation_manifest)) {
        NA_character_
      } else {
        as.character(observation_manifest$arm_no[[1L]])
      },
      outcome_no = if (is.null(observation_manifest)) {
        NA_character_
      } else {
        as.character(observation_manifest$outcome_no[[1L]])
      },
      pareto_k = as.numeric(record$pareto_k),
      adapt_delta = as.numeric(config$adapt_delta_floor),
      max_treedepth = as.integer(config$max_treedepth_floor),
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      original_fit_exists = original_exists,
      integrated_loo_exists = integrated_exists,
      valid = valid
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  if (
    nrow(checks) != expected || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(expected)) ||
      !identical(as.integer(outcome_counts), c(4L, 3L, 0L, 2L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$integrated_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L
  ) {
    stop("The updated EIV exact-LOO preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_eiv_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "eiv", "loo_exact",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 9L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated EIV exact preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_eiv_exact_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "eiv", "loo",
    "exact", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

fit_and_score_update_2026_eiv_exact_to_file <- function(
    record,
    model,
    config,
    loo_draws = 4000L,
    output_file = update_2026_eiv_exact_file_path(record)
) {
  if (is_sensitivity_pilot_skip(model)) {
    stop("The updated EIV Stan model was not compiled.", call. = FALSE)
  }
  stan_file <- model$stan_file()
  if (length(stan_file) != 1L || is.na(stan_file) || !file.exists(stan_file)) {
    stop("The updated EIV Stan source is unavailable.", call. = FALSE)
  }
  stan_sha256 <- digest::digest(
    file = stan_file,
    algo = "sha256",
    serialize = FALSE
  )
  existing <- if (file.exists(output_file)) {
    tryCatch(qs::qread(output_file), error = function(e) NULL)
  } else {
    NULL
  }
  if (
    is.list(existing) && identical(existing$family, "eiv") &&
      as.integer(existing$case_id) == record$case_id &&
      identical(existing$outcome, record$outcome) &&
      as.integer(existing$model_index) == record$model_index &&
      is.data.frame(existing$exact_refit_summary) &&
      nrow(existing$exact_refit_summary) == 1L &&
      existing$exact_refit_summary$observation[[1L]] == record$observation &&
      identical(existing$source_fit_sha256, record$original_fit_sha256) &&
      identical(
        existing$source_integrated_loo_sha256,
        record$integrated_loo_sha256
      ) &&
      identical(existing$stan_sha256, stan_sha256) &&
      identical(
        existing$exact_sampling,
        list(
          iter_warmup = as.integer(config$iter_warmup),
          iter_sampling = as.integer(config$iter_sampling),
          adapt_delta_floor = as.numeric(config$adapt_delta_floor),
          max_treedepth_floor = as.integer(config$max_treedepth_floor),
          loo_draws = as.integer(loo_draws)
        )
      )
  ) {
    return(normalizePath(output_file, winslash = "/", mustWork = TRUE))
  }
  candidate <- read_eiv_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  input <- integrated$input
  observation <- as.integer(record$observation)
  pareto_values <- loo::pareto_k_values(integrated$loo)
  if (
    !identical(candidate$outcome, record$outcome) ||
      as.integer(candidate$model_index) != record$model_index ||
      !identical(unname(as.character(integrated$family)), "eiv") ||
      !identical(integrated$outcome, record$outcome) ||
      as.integer(integrated$model_index) != record$model_index ||
      !is.list(input) || !is.list(input$stan_data) ||
      observation < 1L || observation > input$stan_data$N ||
      pareto_values[[observation]] <= 0.7
  ) {
    stop("An updated EIV exact-LOO record failed validation.",
         call. = FALSE)
  }
  training_data <- make_eiv_leave_one_out_stan_data(
    input$stan_data,
    observation
  )
  specification <- candidate$specification
  specification$adapt_delta <- max(
    as.numeric(specification$adapt_delta),
    as.numeric(config$adapt_delta_floor)
  )
  specification$max_treedepth <- max(
    as.integer(specification$max_treedepth),
    as.integer(config$max_treedepth_floor)
  )
  temporary_output_dir <- tempfile(
    pattern = sprintf("eiv-exact-case-%03d-", record$case_id),
    tmpdir = tempdir()
  )
  dir.create(temporary_output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temporary_output_dir, recursive = TRUE, force = TRUE),
          add = TRUE)
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
    output_dir = temporary_output_dir
  )
  fitted_candidate <- list(
    family = "eiv",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    specification = specification,
    grid = input$grid,
    output_dir = temporary_output_dir,
    fit = fit
  )
  refit <- list(
    case_id = as.integer(record$case_id),
    family = "eiv",
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    pareto_k = as.numeric(record$pareto_k),
    fitted_candidate = fitted_candidate,
    input = input,
    integrated_result = integrated,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  correction <- compute_eiv_exact_loo_case(
    refit,
    loo_draws = loo_draws
  )
  correction$source_fit_sha256 <- record$original_fit_sha256
  correction$source_integrated_loo_sha256 <-
    record$integrated_loo_sha256
  correction$stan_sha256 <- stan_sha256
  correction$exact_sampling <- list(
    iter_warmup = as.integer(config$iter_warmup),
    iter_sampling = as.integer(config$iter_sampling),
    adapt_delta_floor = as.numeric(config$adapt_delta_floor),
    max_treedepth_floor = as.integer(config$max_treedepth_floor),
    loo_draws = as.integer(loo_draws)
  )
  if (
    as.integer(correction$case_id) != record$case_id ||
      !inherits(correction$loo, "loo") ||
      nrow(correction$exact_refit_summary) != 1L ||
      correction$exact_refit_summary$observation[[1L]] != observation
  ) {
    stop("The updated EIV exact correction is invalid.", call. = FALSE)
  }
  write_atomic_qs(correction, output_file)
}

audit_update_2026_eiv_exact_files <- function(paths, integrated_audit) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    length(paths) != expected || expected != 9L || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0)
  ) {
    stop("The updated EIV exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  manifest <- purrr::map2_dfr(paths, corrections, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      observation =
        as.integer(result$exact_refit_summary$observation[[1L]]),
      correction_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys)
  ) {
    stop("The updated EIV exact audit is invalid.", call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = sum(diagnostic_summary$diagnostic_pass),
    exact_refits_requiring_remediation =
      sum(!diagnostic_summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass =
      all(diagnostic_summary$diagnostic_pass)
  )
}

write_update_2026_eiv_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "eiv", "loo_exact"
    )
) {
  if (
    audit$exact_refits != 9L ||
      nrow(audit$diagnostic_summary) != 9L ||
      nrow(audit$correction_manifest) != 9L
  ) {
    stop("Cannot write an incomplete updated EIV exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

read_update_2026_eiv_exact_corrections <- function(manifest) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 9L ||
      any(!file.exists(manifest$correction_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("The updated EIV exact-correction manifest is invalid.",
         call. = FALSE)
  }
  corrections <- purrr::map(manifest$correction_path, qs::qread)
  corrections <- corrections[order(purrr::map_int(corrections, "case_id"))]
  if (
    !identical(
      unname(purrr::map_int(corrections, "case_id")),
      seq_len(9L)
    ) ||
      any(!purrr::map_lgl(corrections, ~ inherits(.x$loo, "loo"))) ||
      any(!purrr::map_lgl(
        corrections,
        ~ identical(unname(as.character(.x$family)), "eiv")
      ))
  ) {
    stop("The updated EIV exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections
}

make_update_2026_eiv_final_loo_audit <- function(
    integrated_manifest,
    exact_corrections,
    exact_audit
) {
  if (
    !isTRUE(exact_audit$all_exact_refit_diagnostics_pass) ||
      exact_audit$exact_refits != 9L
  ) {
    stop("Final EIV LOO requires nine passing exact diagnostics.",
         call. = FALSE)
  }
  final <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest = integrated_manifest,
    exact_results = exact_corrections,
    family = "eiv"
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- integrated_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "outcome_effects_integrated_latent_exposure_psis_plus_exact",
        "outcome_effects_integrated_latent_exposure_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 9L
  ) {
    stop("The final updated EIV LOO is not stacking-ready.",
         call. = FALSE)
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 9L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_eiv_final_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "eiv", "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 9L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop("Cannot write an incomplete final updated EIV LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_eiv_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "eiv",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L ||
      any(weights$family != "eiv") ||
      anyDuplicated(paste(weights$outcome, weights$model_index)) ||
      length(outcome_sums) != 4L ||
      any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop("The updated EIV stacking weights are invalid.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_eiv_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    inputs,
    n_draws = 4000L,
    seed = 31988L
) {
  records <- make_eiv_full_prediction_records(
    accepted_fit_manifest = accepted_fit_manifest,
    stacking_weights = stacking_weights,
    n_draws = n_draws,
    seed = seed
  )
  input_keys <- purrr::map_chr(inputs, function(input) {
    paste(
      input$specification$outcome,
      input$specification$model_index,
      sep = "::"
    )
  })
  records <- purrr::map(records, function(record) {
    key <- paste(record$outcome, record$model_index, sep = "::")
    position <- match(key, input_keys)
    if (is.na(position)) {
      stop("An updated EIV prediction input is unavailable.",
           call. = FALSE)
    }
    record$protein_max <- max(
      inputs[[position]]$grid$protein_intake_raw
    )
    record
  })
  if (
    length(inputs) != 104L || anyDuplicated(input_keys) ||
      length(records) < 4L ||
      sum(purrr::map_int(records, "draw_count")) != 16000L
  ) {
    stop("The updated EIV prediction allocation is incomplete.",
         call. = FALSE)
  }
  records
}

make_update_2026_eiv_reporting_preflight <- function(records) {
  if (!is.list(records) || length(records) < 4L) {
    stop("The updated EIV reporting allocation is incomplete.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    tibble::tibble(
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      family = unname(as.character(record$family)),
      prediction_type = unname(as.character(record$prediction_type)),
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      protein_max = as.numeric(record$protein_max),
      grid_step = as.numeric(record$grid_step),
      fit_path = record$fit_path,
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path) &&
        identical(unname(as.character(record$family)), "eiv") &&
        identical(unname(as.character(record$prediction_type)), "eiv") &&
        nchar(record$fit_sha256) == 64L &&
        is.finite(record$protein_min) &&
        is.finite(record$protein_max) &&
        record$protein_max > record$protein_min
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      retained_weight = sum(weight),
      common_protein_max = dplyr::first(protein_max),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      any(!is.finite(outcome_audit$retained_weight)) ||
      any(outcome_audit$retained_weight <= 0) ||
      any(outcome_audit$retained_weight > 1 + 1e-8) ||
      any(!is.finite(outcome_audit$common_protein_max)) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop("The updated EIV reporting preflight failed.", call. = FALSE)
  }
  checks
}

write_update_2026_eiv_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "eiv", "reporting",
      "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(checks$allocated_draws, checks$outcome, sum)
  outcome_weights <- tapply(checks$weight, checks$outcome, sum)
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      length(outcome_weights) != 4L ||
      any(!is.finite(outcome_weights)) || any(outcome_weights <= 0) ||
      any(outcome_weights > 1 + 1e-8) || !all(checks$valid)
  ) {
    stop("Cannot write an invalid updated EIV reporting preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_eiv_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "eiv",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_eiv_prediction_draw_file <- function(record) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !identical(unname(as.character(record$family)), "eiv") ||
      !identical(unname(as.character(record$prediction_type)), "eiv") ||
      !identical(
        digest::digest(
          file = record$fit_path,
          algo = "sha256",
          serialize = FALSE
        ),
        record$fit_sha256
      )
  ) {
    stop("An updated EIV prediction record is invalid.", call. = FALSE)
  }
  output_file <- update_2026_eiv_reporting_path(
    outcome = record$outcome,
    stage = "prediction_draws",
    model_index = record$model_index
  )
  existing <- if (file.exists(output_file)) {
    tryCatch(
      read_sensitivity_full_prediction_draw_file(output_file, "eiv"),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (
    is.list(existing) && existing$model_index == record$model_index &&
      identical(existing$outcome, record$outcome) &&
      identical(existing$fit_sha256, record$fit_sha256) &&
      existing$draw_count == record$draw_count &&
      existing$prediction_seed == record$prediction_seed
  ) {
    return(normalizePath(output_file, winslash = "/", mustWork = TRUE))
  }
  candidate <- read_eiv_fit_file(record$fit_path)
  if (
    candidate$model_index != record$model_index ||
      !identical(candidate$outcome, record$outcome) ||
      !identical(candidate$model_name, record$model_name)
  ) {
    stop("An updated EIV prediction fit does not match its record.",
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
  grid <- candidate$grid |>
    dplyr::transmute(
      resistance_exercise_code,
      resistance_exercise = ifelse(
        resistance_exercise_code == 0,
        "YES",
        "NO"
      ),
      protein_intake = protein_intake_raw
    )
  if (
    nrow(draws) != record$draw_count || ncol(draws) != nrow(grid) ||
      anyNA(draws) || any(!is.finite(draws)) ||
      abs(max(grid$protein_intake) - record$protein_max) > 1e-10
  ) {
    stop("Updated EIV prediction draws are incomplete.",
         call. = FALSE)
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
  write_atomic_qs(reduction, output_file)
}

write_update_2026_eiv_model_average_file <- function(record) {
  reductions <- purrr::map(
    record$prediction_manifest$prediction_path,
    read_sensitivity_full_prediction_draw_file,
    family = "eiv"
  )
  key_reference <- reductions[[1L]]$grid
  if (!all(purrr::map_lgl(
    reductions,
    ~ identical(.x$grid, key_reference)
  ))) {
    stop("Updated EIV prediction grids do not match within outcome.",
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
          file = path,
          algo = "sha256",
          serialize = FALSE
        ),
        within_model_draw = seq_len(reduction$draw_count)
      )
    }
  )
  if (
    nrow(draws) != record$n_draws ||
      nrow(draw_manifest) != record$n_draws || anyNA(draws)
  ) {
    stop("The updated EIV posterior mixture is incomplete.",
         call. = FALSE)
  }
  set.seed(record$averaging_seed)
  draw_order <- sample(seq_len(record$n_draws))
  draws <- draws[draw_order, , drop = FALSE]
  draw_manifest <- draw_manifest[draw_order, , drop = FALSE] |>
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
  write_atomic_qs(
    output,
    update_2026_eiv_reporting_path(
      outcome = record$outcome,
      stage = "model_average"
    )
  )
}

write_update_2026_eiv_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path,
    family = "eiv"
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (is.null(data) || !is.data.frame(data)) {
    stop("The updated EIV reporting data are unavailable.",
         call. = FALSE)
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
    family = "eiv",
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
  write_atomic_qs(
    report,
    update_2026_eiv_reporting_path(
      outcome = outcome,
      stage = "reporting"
    )
  )
}

write_update_2026_eiv_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "eiv", "reporting"
    )
) {
  if (
    !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L
  ) {
    stop("The updated EIV reporting audit is incomplete.",
         call. = FALSE)
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "eiv"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir,
      "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir,
      "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir,
      "common_exposure_support.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_eiv_verification_checkpoint <- function(
    data_by_outcome,
    input_preflight,
    source_manifest,
    review_decisions,
    preparation_audit,
    specifications,
    accepted_fit_manifest,
    diagnostic_audit,
    integrated_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  raw_analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  # This is the same eligibility filter used by make_eiv_candidate_data().
  # The combined outcome objects intentionally retain rows without usable
  # exposure-error inputs for other analyses, so those rows are not EIV data.
  analysis_data <- raw_analysis_data |>
    dplyr::filter(
      is.finite(yi), is.finite(vi), vi > 0,
      is.finite(protein_intake), is.finite(protein_intake_se),
      protein_intake_se > 0
    )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "eiv"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(data) {
    data |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_audit$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_eiv_name <- paste0(
    name_alignment$specification_name,
    " [EIV Stan]"
  )
  weight_sums <- stacking_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      .by = outcome
    )
  reporting_model_weights <- bind_report("model_weights")
  reporting_weight_sums <- reporting_model_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  sample_preflight <- input_preflight |>
    dplyr::distinct(
      outcome, rows, studies, arms, outcome_rows,
      protein_minimum, protein_maximum,
      minimum_protein_se, maximum_protein_se
    )
  analysis_row_counts <- analysis_data |>
    dplyr::count(analysis_outcome, name = "analysis_rows")
  sample_alignment <- sample_preflight |>
    dplyr::select(outcome, preflight_rows = rows) |>
    dplyr::left_join(
      analysis_row_counts,
      by = c("outcome" = "analysis_outcome")
    )

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "all four EIV datasets contain positive exposure-error inputs",
      "the combined-data preparation audit remains valid",
      "all five EIV registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist and have checksums",
      "all final model diagnostics pass",
      "integrated LOO contains 104 reductions and 9 audited high-k cases",
      "all 9 exact refits pass and compact corrections exist",
      "final corrected LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four reporting audits pass and retain all model weights"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      identical(names(data_by_outcome), outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L) &&
        nrow(sample_alignment) == 4L &&
        all(sample_alignment$preflight_rows ==
              sample_alignment$analysis_rows) &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) && all(analysis_data$vi > 0) &&
        all(is.finite(analysis_data$protein_intake)) &&
        all(is.finite(analysis_data$protein_intake_se)) &&
        all(analysis_data$protein_intake_se > 0) &&
        nrow(input_preflight) == 104L && all(input_preflight$valid) &&
        nrow(sample_preflight) == 4L &&
        setequal(sample_preflight$outcome, outcomes) &&
        all(sample_preflight$minimum_protein_se > 0),
      all(preparation_audit$checks$passed),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_eiv_name) &&
        all(name_alignment$diagnostic_name == expected_eiv_name) &&
        all(name_alignment$loo_name == expected_eiv_name) &&
        all(name_alignment$weight_name == expected_eiv_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(accepted_fit_manifest$family == "eiv") &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(accepted_fit_manifest$bytes > 0) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      nrow(diagnostic_audit$diagnostic_summary) == 104L &&
        isTRUE(diagnostic_audit$all_sampling_diagnostics_pass) &&
        all(diagnostic_audit$diagnostic_summary$diagnostic_pass),
      integrated_loo_audit$models == 104L &&
        integrated_loo_audit$exact_refits_required == 9L &&
        nrow(integrated_loo_audit$high_k_cases) == 9L,
      exact_audit$exact_refits == 9L &&
        exact_audit$exact_refits_passing == 9L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        nrow(exact_audit$correction_manifest) == 9L &&
        all(file.exists(exact_audit$correction_manifest$correction_path)) &&
        all(exact_audit$correction_manifest$bytes > 0) &&
        all(nchar(exact_audit$correction_manifest$sha256) == 64L),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(file.exists(prediction_manifest$prediction_path)) &&
        all(prediction_manifest$bytes > 0) &&
        all(nchar(prediction_manifest$sha256) == 64L),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L &&
        nrow(reporting_model_weights) == 104L &&
        nrow(reporting_weight_sums) == 4L &&
        all(reporting_weight_sums$models == 26L) &&
        all(abs(reporting_weight_sums$weight_sum - 1) < 1e-8)
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated EIV verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    review_decisions = review_decisions,
    preparation_checks = preparation_audit$checks,
    preparation_outcome_counts = preparation_audit$outcome_counts,
    preparation_update_study_counts = preparation_audit$update_study_counts,
    eiv_input_preflight = input_preflight,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_audit$diagnostic_summary,
    sampling_parameter_failures = diagnostic_audit$parameter_failures,
    initial_integrated_loo_model_summary = integrated_loo_audit$model_summary,
    initial_integrated_loo_outcome_summary =
      integrated_loo_audit$outcome_summary,
    initial_high_pareto_k_cases = integrated_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_correction_manifest = exact_audit$correction_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    final_exact_refit_summary = final_loo_audit$exact_refit_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    reporting_model_weights = reporting_model_weights,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final updated EIV checkpoint checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    review_decisions = "Eligibility decisions for updated-search studies.",
    preparation_checks = "Combined-data preparation invariants.",
    preparation_outcome_counts = "Analysis counts by primary outcome.",
    preparation_update_study_counts = "Updated-study contribution audit.",
    eiv_input_preflight = "EIV data, exposure-error, and schedule checks.",
    analysis_data = "Eligible combined outcome rows supplied to EIV models.",
    model_specifications = "Final 104-model EIV specification registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit path, stage, size, and checksum.",
    sampling_diagnostic_summary = "Final model-level sampling diagnostics.",
    sampling_parameter_failures = "Parameter-level diagnostic failures.",
    initial_integrated_loo_model_summary = "Initial integrated PSIS summaries.",
    initial_integrated_loo_outcome_summary = "Initial LOO audit by outcome.",
    initial_high_pareto_k_cases = "The 9 cases selected for exact refitting.",
    exact_refit_diagnostic_summary = "Final diagnostics for 9 exact refits.",
    exact_correction_manifest = "Compact exact-correction paths and checksums.",
    final_loo_manifest = "Final exact-corrected LOO lineage.",
    final_loo_model_summary = "Final LOO summary for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    final_exact_refit_summary = "Final exact-correction audit by case.",
    stacking_weights = "All 104 final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    reporting_model_weights = "All 104 weights retained in reporting objects.",
    dose_response_predictions = "EIV model-averaged summaries.",
    dose_response_increments = "Increments from the 0.8 g/kg/day reference.",
    exposure_support_by_group = "Observed EIV exposure support by RT group.",
    common_exposure_support = "Common EIV exposure support by outcome.",
    pairwise_contrasts = "Model-averaged pairwise contrasts and probabilities.",
    pairwise_support = "Support limits for pairwise contrasts.",
    local_equivalence = "Local practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    next_family = paste(
      "The residual-correlation sensitivity family begins after this EIV",
      "verification checkpoint is frozen and the pruning plan is audited."
    )
  )
}

make_update_2026_eiv_pruning_plan <- function(
    accepted_fit_manifest,
    verification_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "eiv"
    )
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "eiv_update_2026"),
    "manifest.csv"
  )
  if (
    length(verification_files) != 35L ||
      any(!file.exists(verification_files))
  ) {
    stop("The updated EIV verification bundle is incomplete.",
         call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "eiv"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated EIV archive.",
         call. = FALSE)
  }

  fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  if (any(!dir.exists(fit_roots))) {
    stop("An expected updated EIV fit directory is missing.",
         call. = FALSE)
  }
  deletion_paths <- unlist(
    purrr::map(
      fit_roots,
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  deletion_paths <- deletion_paths[
    file.exists(deletion_paths) & !file.info(deletion_paths)$isdir
  ]
  deletion_paths <- normalizePath(
    deletion_paths,
    winslash = "/",
    mustWork = TRUE
  )
  fit_paths <- deletion_paths[basename(deletion_paths) == "fit.qs"]
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(deletion_paths) != 570L || length(fit_paths) != 114L ||
      length(accepted_paths) != 104L ||
      anyDuplicated(deletion_paths) || anyDuplicated(fit_paths) ||
      anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% fit_paths) ||
      sum(!fit_paths %in% accepted_paths) != 10L
  ) {
    stop("Updated EIV fit retention does not match its lineage.",
         call. = FALSE)
  }

  fit_manifest <- tibble::tibble(
    path = fit_paths,
    relative_path = substring(fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", fit_paths, fixed = TRUE) ~ "initial",
      grepl("/remedies/remedy_1/", fit_paths, fixed = TRUE) ~ "remedy_1",
      TRUE ~ "unknown"
    ),
    outcome = basename(dirname(dirname(fit_paths))),
    model_index = as.integer(sub(
      "^model_", "", basename(dirname(fit_paths))
    )),
    accepted_fit = fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(fit_paths)$size),
    lean_action = "delete frozen EIV fitted object"
  ) |>
    dplyr::arrange(archive_stage, outcome, model_index)
  if (
    any(fit_manifest$archive_stage == "unknown") ||
      anyNA(fit_manifest$model_index) ||
      !setequal(fit_manifest$outcome, c(
        "whole_body", "hypertrophy", "strength", "function"
      )) ||
      !identical(
        as.integer(table(factor(
          fit_manifest$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 10L)
      ) ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "initial") != 94L ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "remedy_1") != 10L
  ) {
    stop("The updated EIV fit stages are not as expected.",
         call. = FALSE)
  }

  fit_directories <- dirname(fit_manifest$path)
  file_model_directories <- vapply(
    deletion_paths,
    function(path) {
      if (basename(path) == "fit.qs") {
        dirname(path)
      } else {
        dirname(dirname(path))
      }
    },
    character(1L),
    USE.NAMES = FALSE
  )
  fit_positions <- match(file_model_directories, fit_directories)
  file_types <- dplyr::case_when(
    basename(deletion_paths) == "fit.qs" ~ "fit_object",
    tolower(tools::file_ext(deletion_paths)) == "csv" &
      basename(dirname(deletion_paths)) == "cmdstan" ~
        "cmdstan_chain_csv",
    TRUE ~ "unknown"
  )
  deletion_manifest <- tibble::tibble(
    path = deletion_paths,
    relative_path = substring(
      deletion_paths,
      nchar(resolved_root) + 2L
    ),
    archive_stage = fit_manifest$archive_stage[fit_positions],
    outcome = fit_manifest$outcome[fit_positions],
    model_index = fit_manifest$model_index[fit_positions],
    accepted_fit = fit_manifest$accepted_fit[fit_positions],
    file_type = file_types,
    bytes = as.numeric(file.info(deletion_paths)$size),
    lean_action = "delete refittable EIV sampling archive"
  ) |>
    dplyr::arrange(
      archive_stage, outcome, model_index,
      factor(file_type, levels = c("fit_object", "cmdstan_chain_csv")),
      relative_path
    )
  if (
    anyNA(fit_positions) || any(deletion_manifest$file_type == "unknown") ||
      sum(deletion_manifest$file_type == "fit_object") != 114L ||
      sum(deletion_manifest$file_type == "cmdstan_chain_csv") != 456L ||
      !identical(
        as.integer(table(factor(
          deletion_manifest$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(520L, 50L)
      ) ||
      sum(deletion_manifest$accepted_fit) != 520L ||
      sum(!deletion_manifest$accepted_fit) != 50L ||
      any(deletion_manifest$bytes <= 0)
  ) {
    stop("The updated EIV deletion inventory is incomplete.",
         call. = FALSE)
  }

  compact_roots <- file.path(archive_root, c("loo", "reporting"))
  if (any(!dir.exists(compact_roots))) {
    stop("An expected compact EIV archive directory is missing.",
         call. = FALSE)
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[
    file.exists(compact_paths) & !file.info(compact_paths)$isdir
  ]
  compact_paths <- normalizePath(
    compact_paths,
    winslash = "/",
    mustWork = TRUE
  )
  compact_manifest <- tibble::tibble(
    path = compact_paths,
    bytes = as.numeric(file.info(compact_paths)$size),
    sha256 = purrr::map_chr(
      compact_paths,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    ),
    action = "retain compact EIV LOO/reporting archive"
  )
  if (
    nrow(compact_manifest) != 154L ||
      any(compact_manifest$bytes <= 0) ||
      any(nchar(compact_manifest$sha256) != 64L)
  ) {
    stop("The compact EIV archive inventory is incomplete.",
         call. = FALSE)
  }

  summary <- tibble::tibble(
    mode = "lean",
    model_archives_deleted = nrow(fit_manifest),
    accepted_model_archives_deleted = sum(fit_manifest$accepted_fit),
    superseded_model_archives_deleted = sum(!fit_manifest$accepted_fit),
    archive_files_deleted = nrow(deletion_manifest),
    fit_object_files_deleted = sum(
      deletion_manifest$file_type == "fit_object"
    ),
    cmdstan_csv_files_deleted = sum(
      deletion_manifest$file_type == "cmdstan_chain_csv"
    ),
    fit_archive_files_retained = 0L,
    gib_reclaimed = sum(deletion_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained = sum(compact_manifest$bytes) / 1024^3,
    verification_files_retained = length(verification_files),
    restoration = paste(
      "Reported results and all 104 weights remain frozen; restoring full",
      "fits and chain CSVs requires refitting the EIV family."
    )
  )

  list(
    summary = summary,
    fit_file_manifest = fit_manifest,
    deletion_file_manifest = deletion_manifest,
    compact_archive_manifest = compact_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_eiv_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "eiv", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    fit_file_manifest = file.path(output_dir, "fit_file_manifest.csv"),
    deletion_file_manifest = file.path(
      output_dir, "deletion_file_manifest.csv"
    ),
    compact_archive_manifest = file.path(
      output_dir, "compact_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(plan$fit_file_manifest, paths[["fit_file_manifest"]])
  readr::write_csv(
    plan$deletion_file_manifest,
    paths[["deletion_file_manifest"]]
  )
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_eiv_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "eiv", "pruning"
    )
) {
  required_confirmation <- "DELETE ALL UPDATED EIV FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "eiv"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the EIV archive.",
         call. = FALSE)
  }
  fit_manifest <- plan$fit_file_manifest
  selected <- plan$deletion_file_manifest
  if (
    nrow(fit_manifest) != 114L ||
      sum(fit_manifest$accepted_fit) != 104L ||
      sum(!fit_manifest$accepted_fit) != 10L ||
      nrow(selected) != 570L ||
      sum(selected$file_type == "fit_object") != 114L ||
      sum(selected$file_type == "cmdstan_chain_csv") != 456L ||
      sum(selected$accepted_fit) != 520L ||
      sum(!selected$accepted_fit) != 50L
  ) {
    stop("The EIV pruning set is inconsistent with its audit.",
         call. = FALSE)
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  pending_receipt_path <- file.path(
    receipt_dir,
    "executed_lean.pending.csv"
  )
  if (file.exists(receipt_path) || file.exists(pending_receipt_path)) {
    stop("An updated EIV lean-pruning receipt already exists.",
         call. = FALSE)
  }
  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  current_paths <- unlist(
    purrr::map(
      file.path(plan$archive_root, c("candidates", "remedies")),
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  current_paths <- current_paths[
    file.exists(current_paths) & !file.info(current_paths)$isdir
  ]
  current_paths <- normalizePath(
    current_paths,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  permitted_subtree <- grepl(
    "/(candidates|remedies)/",
    selected_paths
  )
  if (
    !setequal(selected_paths, current_paths) ||
      any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted_subtree) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the audited EIV fit inventory.",
         call. = FALSE)
  }
  compact_paths <- plan$compact_archive_manifest$path
  message("Revalidating 154 compact EIV artifacts...")
  if (any(!file.exists(compact_paths))) {
    stop("A retained compact EIV artifact failed revalidation.",
         call. = FALSE)
  }
  compact_hashes <- purrr::map_chr(
    compact_paths,
    ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
  )
  if (any(compact_hashes != plan$compact_archive_manifest$sha256)) {
    stop("A retained compact EIV artifact failed revalidation.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  message(
    "Hashing 570 EIV sampling files (approximately 29.7 GiB) ",
    "before deletion..."
  )
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, pending_receipt_path)

  message("Deleting the 570 audited EIV sampling files...")
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop(
      "At least one EIV sampling archive could not be deleted; ",
      "the pending receipt records the attempted inventory.",
      call. = FALSE
    )
  }
  if (any(!file.exists(compact_paths))) {
    stop("A retained EIV artifact is unexpectedly missing.",
         call. = FALSE)
  }
  final_compact_hashes <- purrr::map_chr(
    compact_paths,
    ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
  )
  if (
    any(final_compact_hashes != plan$compact_archive_manifest$sha256) ||
      !file.exists(plan$verification_manifest)
  ) {
    stop("A retained EIV artifact is unexpectedly missing.",
         call. = FALSE)
  }

  readr::write_csv(receipt, receipt_path)
  if (!file.remove(pending_receipt_path)) {
    warning("The completed EIV pending receipt could not be removed.")
  }
  message("EIV lean pruning completed; compact reporting artifacts remain.")
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

# Combined-data residual-correlation sensitivity refit ----

make_update_2026_rescor_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated residual-correlation registry is incomplete.",
         call. = FALSE)
  }
  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  remedy_keys <- unname(purrr::map_chr(primary_remedy_specifications, key))
  if (
    anyDuplicated(primary_keys) || anyDuplicated(remedy_keys) ||
      !identical(remedy_keys, "function::10") ||
      !all(remedy_keys %in% primary_keys)
  ) {
    stop("The updated residual-correlation primary replacement is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  remedy_1 <- c(
    "whole_body::1", "whole_body::10", "whole_body::17",
    "hypertrophy::6", "hypertrophy::8", "hypertrophy::9",
    "hypertrophy::21", "hypertrophy::22",
    "strength::3", "strength::17",
    "function::16", "function::18"
  )
  remedy_2 <- c("function::6", "function::10")
  structural <- "function::23"
  historical_keys <- c(remedy_1, remedy_2, structural)

  specifications <- purrr::map(specifications, function(specification) {
    specification_key <- key(specification)
    schedule <- dplyr::case_when(
      specification_key %in% remedy_1 ~ "remedy_1",
      specification_key %in% remedy_2 ~ "remedy_2",
      specification_key %in% structural ~ "structural_uncorrelated",
      TRUE ~ "base"
    )
    specification$update_2026_sensitivity_base <- TRUE
    specification$update_2026_rescor_historical_schedule <- schedule
    if (schedule == "structural_uncorrelated") {
      specification$uncorrelated_study_effects <- TRUE
      specification$newly_uncorrelated_study_effects <- TRUE
      specification$rescor_structural_remedy <-
        "independent_study_random_smooth_coefficients"
    }
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  schedules <- purrr::map_chr(
    specifications,
    "update_2026_rescor_historical_schedule"
  )
  schedule_counts <- table(factor(
    schedules,
    levels = c(
      "base", "remedy_1", "remedy_2", "structural_uncorrelated"
    )
  ))
  scheduled_keys <- final_keys[schedules != "base"]
  structural_specification <- specifications[[
    match(structural, final_keys)
  ]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !identical(as.integer(schedule_counts), c(89L, 12L, 2L, 1L)) ||
      !setequal(scheduled_keys, historical_keys) ||
      !isTRUE(structural_specification$uncorrelated_study_effects) ||
      !identical(
        structural_specification$rescor_structural_remedy,
        "independent_study_random_smooth_coefficients"
      )
  ) {
    stop("The updated residual-correlation registry failed validation.",
         call. = FALSE)
  }
  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

update_2026_rescor_local_config <- function(specification, config) {
  local_config <- config
  schedule <- specification$update_2026_rescor_historical_schedule
  if (
    length(schedule) != 1L || is.na(schedule) ||
      !schedule %in% c(
        "base", "remedy_1", "remedy_2", "structural_uncorrelated"
      )
  ) {
    stop("The updated residual-correlation schedule is invalid.",
         call. = FALSE)
  }
  remedy_level <- specification$update_2026_rescor_remedy_level
  if (is.null(remedy_level)) {
    remedy_level <- NA_character_
  }
  if (
    length(remedy_level) != 1L ||
      (!is.na(remedy_level) &&
         !remedy_level %in% c(
           "remedy_1", "remedy_2", "remedy_3", "remedy_4"
         ))
  ) {
    stop("The updated residual-correlation remedy level is invalid.",
         call. = FALSE)
  }
  active_schedule <- if (is.na(remedy_level)) schedule else remedy_level
  schedule_offset <- switch(
    active_schedule,
    base = 0L,
    remedy_1 = 1000L,
    remedy_2 = 2000L,
    remedy_3 = 3000L,
    remedy_4 = 5000L,
    structural_uncorrelated = 4000L
  )
  if (active_schedule == "remedy_1") {
    local_config$iter_warmup <- 2000L
    local_config$iter_sampling <- 4000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      15L
    )
  } else if (
    active_schedule %in% c("remedy_2", "structural_uncorrelated")
  ) {
    local_config$iter_warmup <- 3000L
    local_config$iter_sampling <- 8000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.9999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      16L
    )
  } else if (active_schedule == "remedy_3") {
    local_config$iter_warmup <- 4000L
    local_config$iter_sampling <- 12000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.99999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      17L
    )
  } else if (active_schedule == "remedy_4") {
    local_config$iter_warmup <- 6000L
    local_config$iter_sampling <- 12000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      18L
    )
  }
  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the residual-correlation configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed + schedule_offset +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

validate_update_2026_rescor_inputs <- function(inputs, config) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    local_config <- update_2026_rescor_local_config(
      specification,
      config
    )
    stan_data <- input$stan_data
    cluster_manifest <- input$cluster_manifest
    schedule <- specification$update_2026_rescor_historical_schedule
    protein_median <- min(input$grid$protein_intake_raw) -
      min(input$grid$protein_intake)
    prior_median <- stan_data$rho_lower +
      (1 - stan_data$rho_lower) *
      stats::plogis(stan_data$rho_raw_prior_location)
    expected_cluster_start <- c(
      1L,
      utils::head(cluster_manifest$cluster_end, -1L) + 1L
    )
    structural_schedule <- schedule == "structural_uncorrelated"
    effective_adapt_delta <- max(
      as.numeric(specification$adapt_delta),
      as.numeric(local_config$adapt_delta_floor)
    )
    effective_max_treedepth <- max(
      as.integer(specification$max_treedepth),
      as.integer(local_config$max_treedepth_floor)
    )
    valid <-
      identical(input$family, "rescor") && is.list(stan_data) &&
      stan_data$N > 0L && stan_data$N_study > 1L &&
      stan_data$N_cluster == stan_data$N_study &&
      stan_data$N_arm > 0L && stan_data$N_outcome > 0L &&
      length(stan_data$y) == stan_data$N &&
      length(stan_data$se_y) == stan_data$N &&
      all(is.finite(stan_data$y)) &&
      all(is.finite(stan_data$se_y)) && all(stan_data$se_y > 0) &&
      all(is.finite(stan_data$x)) &&
      nrow(cluster_manifest) == stan_data$N_cluster &&
      identical(
        as.integer(cluster_manifest$cluster_start),
        as.integer(expected_cluster_start)
      ) &&
      max(cluster_manifest$cluster_end) == stan_data$N &&
      all(cluster_manifest$cluster_size >= 1L) &&
      is.finite(stan_data$rho_lower) &&
      stan_data$rho_lower > -1 && stan_data$rho_lower < 0 &&
      is.finite(prior_median) && abs(prior_median) < 1e-10 &&
      identical(
        as.integer(stan_data$correlated_study_effects),
        as.integer(!isTRUE(specification$uncorrelated_study_effects))
      ) &&
      (!structural_schedule ||
         isTRUE(specification$uncorrelated_study_effects)) &&
      is.finite(effective_adapt_delta) &&
      effective_adapt_delta > 0 && effective_adapt_delta < 1 &&
      effective_max_treedepth >= 1L

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      rows = as.integer(stan_data$N),
      studies = as.integer(stan_data$N_study),
      arms = as.integer(stan_data$N_arm),
      outcome_rows = as.integer(stan_data$N_outcome),
      minimum_cluster_size = min(cluster_manifest$cluster_size),
      maximum_cluster_size = max(cluster_manifest$cluster_size),
      rho_lower = as.numeric(stan_data$rho_lower),
      rho_prior_median = as.numeric(prior_median),
      protein_minimum = min(stan_data$x) + protein_median,
      protein_maximum = max(stan_data$x) + protein_median,
      historical_schedule = schedule,
      uncorrelated_study_effects =
        isTRUE(specification$uncorrelated_study_effects),
      adapt_delta = effective_adapt_delta,
      max_treedepth = effective_max_treedepth,
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      chains = as.integer(local_config$chains),
      parallel_chains = as.integer(local_config$parallel_chains),
      valid = valid
    )
  }) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  counts <- table(factor(checks$outcome, levels = outcomes))
  schedule_counts <- table(factor(
    checks$historical_schedule,
    levels = c(
      "base", "remedy_1", "remedy_2", "structural_uncorrelated"
    )
  ))
  sample_counts <- checks |>
    dplyr::distinct(
      outcome, rows, studies, arms, outcome_rows,
      minimum_cluster_size, maximum_cluster_size,
      rho_lower, rho_prior_median,
      protein_minimum, protein_maximum
    )
  structural_schedule_rows <- checks |>
    dplyr::filter(historical_schedule == "structural_uncorrelated")
  uncorrelated_rows <- checks |>
    dplyr::filter(uncorrelated_study_effects)
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      !identical(as.integer(schedule_counts), c(89L, 12L, 2L, 1L)) ||
      nrow(sample_counts) != 4L ||
      any(sample_counts$rows < 1L) || any(sample_counts$studies < 2L) ||
      any(sample_counts$maximum_cluster_size <= 1L) ||
      any(abs(sample_counts$rho_prior_median) >= 1e-10) ||
      nrow(structural_schedule_rows) != 1L ||
      structural_schedule_rows$outcome[[1L]] != "function" ||
      structural_schedule_rows$model_index[[1L]] != 23L ||
      nrow(uncorrelated_rows) != 12L ||
      any(checks$chains != 4L) || any(checks$parallel_chains != 4L)
  ) {
    stop("The updated residual-correlation input preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_rescor_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete residual-correlation preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_rescor_candidate_to_file <- function(
    input,
    model,
    config,
    stage = "initial",
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "rescor",
      "candidates"
    )
) {
  if (is_sensitivity_pilot_skip(model)) {
    stop("The updated residual-correlation Stan model was not compiled.",
         call. = FALSE)
  }
  specification <- input$specification
  local_config <- update_2026_rescor_local_config(specification, config)
  specification$adapt_delta <- max(
    as.numeric(specification$adapt_delta),
    as.numeric(local_config$adapt_delta_floor)
  )
  specification$max_treedepth <- max(
    as.integer(specification$max_treedepth),
    as.integer(local_config$max_treedepth_floor)
  )
  specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)
  input_sha256 <- digest::digest(
    list(
      specification = specification,
      stan_data = input$stan_data,
      grid = input$grid,
      cluster_manifest = input$cluster_manifest
    ),
    algo = "sha256"
  )
  stan_file <- model$stan_file()
  if (length(stan_file) != 1L || is.na(stan_file) || !file.exists(stan_file)) {
    stop("The updated residual-correlation Stan source is unavailable.",
         call. = FALSE)
  }
  stan_sha256 <- digest::digest(
    file = stan_file,
    algo = "sha256",
    serialize = FALSE
  )

  candidate_dir <- file.path(
    output_root,
    specification$outcome,
    sprintf("model_%02d", as.integer(specification$model_index))
  )
  output_file <- file.path(candidate_dir, "fit.qs")
  cmdstan_dir <- file.path(candidate_dir, "cmdstan")
  dir.create(cmdstan_dir, recursive = TRUE, showWarnings = FALSE)
  existing_candidate <- if (file.exists(output_file)) {
    tryCatch(read_rescor_fit_file(output_file), error = function(e) NULL)
  } else {
    NULL
  }
  if (
    !is.null(existing_candidate) &&
      identical(existing_candidate$outcome, specification$outcome) &&
      identical(
        as.integer(existing_candidate$model_index),
        as.integer(specification$model_index)
      ) &&
      identical(
        existing_candidate$specification$
          update_2026_rescor_historical_schedule,
        specification$update_2026_rescor_historical_schedule
      ) &&
      identical(
        as.integer(
          existing_candidate$specification$update_2026_iter_sampling
        ),
        as.integer(local_config$iter_sampling)
      ) &&
      identical(existing_candidate$input_sha256, input_sha256) &&
      identical(existing_candidate$stan_sha256, stan_sha256)
  ) {
    return(normalizePath(output_file, winslash = "/", mustWork = TRUE))
  }

  fit <- model$sample(
    data = input$stan_data,
    chains = local_config$chains,
    parallel_chains = local_config$parallel_chains,
    iter_warmup = local_config$iter_warmup,
    iter_sampling = local_config$iter_sampling,
    adapt_delta = specification$adapt_delta,
    max_treedepth = specification$max_treedepth,
    seed = local_config$seed,
    refresh = 100,
    output_dir = cmdstan_dir
  )
  candidate <- list(
    family = "rescor",
    archive_stage = stage,
    outcome = specification$outcome,
    model_index = as.integer(specification$model_index),
    model_name = paste0(
      specification$model_name,
      " [residual correlation by study]"
    ),
    specification = specification,
    grid = input$grid,
    output_dir = normalizePath(
      cmdstan_dir,
      winslash = "/",
      mustWork = TRUE
    ),
    input_sha256 = input_sha256,
    stan_sha256 = stan_sha256,
    fit = fit
  )
  write_atomic_qs(candidate, output_file)
}

write_update_2026_rescor_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "initial_fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$fit_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an incomplete residual-correlation fit manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

write_update_2026_rescor_diagnostic_audit <- function(
    audit,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "initial_diagnostic_summary.csv"
    )
) {
  if (
    !identical(audit$family, "rescor") ||
      nrow(audit$diagnostic_summary) != 104L
  ) {
    stop("Cannot write an incomplete residual-correlation diagnostic audit.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(audit$diagnostic_summary, output_file)
  output_file
}

select_update_2026_rescor_remedy_inputs <- function(
    inputs,
    diagnostic_audit
) {
  expected_keys <- c(
    "whole_body::2", "whole_body::9", "whole_body::11",
    "whole_body::13", "whole_body::14", "whole_body::15",
    "whole_body::16", "whole_body::22",
    "hypertrophy::2", "hypertrophy::4",
    "strength::1", "strength::2", "strength::7",
    "strength::9", "strength::15",
    "function::2", "function::6", "function::13",
    "function::16", "function::17", "function::20", "function::23"
  )
  remedy_2_keys <- "function::16"
  remedy_3_keys <- c("function::6", "function::23")
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(
      outcome,
      model_index = as.integer(model_index),
      key = paste(outcome, model_index, sep = "::"),
      divergent_transitions,
      maximum_treedepth_transitions,
      minimum_ebfmi,
      maximum_rhat,
      minimum_bulk_ess,
      minimum_tail_ess,
      maximum_mcse_sd_ratio
    )
  if (
    !identical(diagnostic_audit$family, "rescor") ||
      nrow(failed) != 22L || anyDuplicated(failed$key) ||
      !setequal(failed$key, expected_keys) ||
      sum(failed$divergent_transitions > 0L) != 11L ||
      any(failed$maximum_treedepth_transitions != 0L) ||
      any(failed$minimum_ebfmi < 0.2)
  ) {
    stop(
      "The updated residual-correlation remedy set no longer matches its audit.",
      call. = FALSE
    )
  }

  selected <- purrr::keep(inputs, function(input) {
    paste(
      input$specification$outcome,
      input$specification$model_index,
      sep = "::"
    ) %in% failed$key
  }) |>
    purrr::map(function(input) {
      key <- paste(
        input$specification$outcome,
        input$specification$model_index,
        sep = "::"
      )
      historical_schedule <-
        input$specification$update_2026_rescor_historical_schedule
      remedy_level <- dplyr::case_when(
        key %in% remedy_3_keys ~ "remedy_3",
        key %in% remedy_2_keys ~ "remedy_2",
        TRUE ~ "remedy_1"
      )
      expected_historical_schedule <- dplyr::case_when(
        key == "function::6" ~ "remedy_2",
        key == "function::16" ~ "remedy_1",
        key == "function::23" ~ "structural_uncorrelated",
        TRUE ~ "base"
      )
      if (!identical(historical_schedule, expected_historical_schedule)) {
        stop(
          "A residual-correlation remedy starting schedule has changed.",
          call. = FALSE
        )
      }
      if (
        key == "function::23" &&
          (!isTRUE(input$specification$uncorrelated_study_effects) ||
             input$stan_data$correlated_study_effects != 0L)
      ) {
        stop(
          "The function model 23 structural remedy was not preserved.",
          call. = FALSE
        )
      }
      input$specification$update_2026_rescor_remedy_level <- remedy_level
      input$specification$update_2026_rescor_remedy <- switch(
        remedy_level,
        remedy_1 = "sampling_and_geometry_escalation_1",
        remedy_2 = "sampling_and_geometry_escalation_2",
        remedy_3 = "sampling_and_geometry_escalation_3"
      )
      input$specification$update_2026_rescor_remedy_stage <- "remedy_1"
      input
    })
  selected_keys <- unname(purrr::map_chr(
    selected,
    ~ paste(
      .x$specification$outcome,
      .x$specification$model_index,
      sep = "::"
    )
  ))
  if (
    length(selected) != 22L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The updated residual-correlation remedy inputs are incomplete.",
         call. = FALSE)
  }
  names(selected) <- paste0(
    purrr::map_chr(selected, ~ .x$specification$outcome),
    "_model_",
    stringr::str_pad(
      purrr::map_int(selected, ~ .x$specification$model_index),
      width = 2L,
      pad = "0"
    )
  )
  selected
}

validate_update_2026_rescor_remedy_inputs <- function(
    inputs,
    initial_manifest,
    config
) {
  expected_keys <- c(
    "whole_body::2", "whole_body::9", "whole_body::11",
    "whole_body::13", "whole_body::14", "whole_body::15",
    "whole_body::16", "whole_body::22",
    "hypertrophy::2", "hypertrophy::4",
    "strength::1", "strength::2", "strength::7",
    "strength::9", "strength::15",
    "function::2", "function::6", "function::13",
    "function::16", "function::17", "function::20", "function::23"
  )
  expected_remedy_2_keys <- "function::16"
  expected_remedy_3_keys <- c("function::6", "function::23")
  initial_keys <- paste(
    initial_manifest$outcome,
    initial_manifest$model_index,
    sep = "::"
  )
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    key <- paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    )
    local_config <- update_2026_rescor_local_config(
      specification,
      config
    )
    remedy_level <- specification$update_2026_rescor_remedy_level
    expected_iter_warmup <- switch(
      remedy_level,
      remedy_1 = 2000L,
      remedy_2 = 3000L,
      remedy_3 = 4000L
    )
    expected_iter_sampling <- switch(
      remedy_level,
      remedy_1 = 4000L,
      remedy_2 = 8000L,
      remedy_3 = 12000L
    )
    expected_adapt_delta_floor <- switch(
      remedy_level,
      remedy_1 = 0.999,
      remedy_2 = 0.9999,
      remedy_3 = 0.99999
    )
    expected_max_treedepth_floor <- switch(
      remedy_level,
      remedy_1 = 15L,
      remedy_2 = 16L,
      remedy_3 = 17L
    )
    initial_index <- match(key, initial_keys)
    initial_fit_exists <-
      !is.na(initial_index) &&
      file.exists(initial_manifest$fit_path[[initial_index]])
    structural_preserved <-
      key != "function::23" ||
      (isTRUE(specification$uncorrelated_study_effects) &&
         input$stan_data$correlated_study_effects == 0L)
    valid <-
      identical(input$family, "rescor") &&
      remedy_level %in% c("remedy_1", "remedy_2", "remedy_3") &&
      identical(specification$update_2026_rescor_remedy_stage, "remedy_1") &&
      input$stan_data$N > 0L && input$stan_data$N_study > 1L &&
      initial_fit_exists && structural_preserved &&
      local_config$iter_warmup == expected_iter_warmup &&
      local_config$iter_sampling == expected_iter_sampling &&
      local_config$adapt_delta_floor >= expected_adapt_delta_floor &&
      local_config$max_treedepth_floor >= expected_max_treedepth_floor

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      historical_schedule =
        specification$update_2026_rescor_historical_schedule,
      remedy_level = remedy_level,
      remedy = specification$update_2026_rescor_remedy,
      uncorrelated_study_effects =
        isTRUE(specification$uncorrelated_study_effects),
      adapt_delta = max(
        as.numeric(specification$adapt_delta),
        as.numeric(local_config$adapt_delta_floor)
      ),
      max_treedepth = max(
        as.integer(specification$max_treedepth),
        as.integer(local_config$max_treedepth_floor)
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      initial_fit_exists = initial_fit_exists,
      valid = valid
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  remedy_2_keys <- keys[checks$remedy_level == "remedy_2"]
  remedy_3_keys <- keys[checks$remedy_level == "remedy_3"]
  if (
    nrow(initial_manifest) != 104L || anyDuplicated(initial_keys) ||
      nrow(checks) != 22L || anyDuplicated(keys) ||
      !setequal(keys, expected_keys) || !all(checks$valid) ||
      !setequal(remedy_2_keys, expected_remedy_2_keys) ||
      !setequal(remedy_3_keys, expected_remedy_3_keys) ||
      sum(checks$remedy_level == "remedy_1") != 19L ||
      sum(checks$remedy_level == "remedy_2") != 1L ||
      sum(checks$remedy_level == "remedy_3") != 2L
  ) {
    stop("The updated residual-correlation remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_rescor_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(checks) != 22L || !all(checks$valid)) {
    stop(
      "Cannot write an incomplete residual-correlation remedy preflight.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

write_update_2026_rescor_remedy_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "remedy_1", "fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 22L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$fit_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an incomplete residual-correlation remedy manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

select_update_2026_rescor_second_remedy_inputs <- function(
    inputs,
    diagnostic_audit
) {
  expected_keys <- c("function::6", "function::23")
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(
      outcome,
      model_index = as.integer(model_index),
      key = paste(outcome, model_index, sep = "::"),
      divergent_transitions,
      maximum_treedepth_transitions,
      minimum_ebfmi,
      maximum_rhat,
      minimum_bulk_ess,
      minimum_tail_ess,
      maximum_mcse_sd_ratio
    )
  expected_divergences <- c(
    "function::6" = 1,
    "function::23" = 2
  )
  observed_divergences <- stats::setNames(
    failed$divergent_transitions,
    failed$key
  )
  if (
    !identical(diagnostic_audit$family, "rescor") ||
      nrow(failed) != 2L || anyDuplicated(failed$key) ||
      !setequal(failed$key, expected_keys) ||
      !identical(
        as.numeric(observed_divergences[expected_keys]),
        as.numeric(expected_divergences[expected_keys])
      ) ||
      any(failed$maximum_treedepth_transitions != 0L) ||
      any(failed$minimum_ebfmi < 0.2) ||
      any(failed$maximum_rhat > 1.01) ||
      any(failed$minimum_bulk_ess < 400) ||
      any(failed$minimum_tail_ess < 400) ||
      any(failed$maximum_mcse_sd_ratio > 0.10)
  ) {
    stop(
      "The second residual-correlation remedy set no longer matches its audit.",
      call. = FALSE
    )
  }

  selected <- purrr::keep(inputs, function(input) {
    paste(
      input$specification$outcome,
      input$specification$model_index,
      sep = "::"
    ) %in% expected_keys
  }) |>
    purrr::map(function(input) {
      key <- paste(
        input$specification$outcome,
        input$specification$model_index,
        sep = "::"
      )
      if (
        !identical(
          input$specification$update_2026_rescor_remedy_level,
          "remedy_3"
        ) ||
          (key == "function::23" &&
             (!isTRUE(input$specification$uncorrelated_study_effects) ||
                input$stan_data$correlated_study_effects != 0L))
      ) {
        stop(
          "A second residual-correlation remedy input has changed.",
          call. = FALSE
        )
      }
      input$specification$update_2026_rescor_remedy_level <- "remedy_4"
      input$specification$update_2026_rescor_remedy <-
        "smaller_step_size_after_isolated_divergence"
      input$specification$update_2026_rescor_remedy_stage <- "remedy_2"
      input
    })
  selected_keys <- unname(purrr::map_chr(
    selected,
    ~ paste(
      .x$specification$outcome,
      .x$specification$model_index,
      sep = "::"
    )
  ))
  if (
    length(selected) != 2L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The second residual-correlation remedy inputs are incomplete.",
         call. = FALSE)
  }
  names(selected) <- paste0(
    purrr::map_chr(selected, ~ .x$specification$outcome),
    "_model_",
    stringr::str_pad(
      purrr::map_int(selected, ~ .x$specification$model_index),
      width = 2L,
      pad = "0"
    )
  )
  selected
}

validate_update_2026_rescor_second_remedy_inputs <- function(
    inputs,
    first_remedy_manifest,
    config
) {
  expected_keys <- c("function::6", "function::23")
  first_remedy_keys <- paste(
    first_remedy_manifest$outcome,
    first_remedy_manifest$model_index,
    sep = "::"
  )
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    key <- paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    )
    local_config <- update_2026_rescor_local_config(
      specification,
      config
    )
    previous_index <- match(key, first_remedy_keys)
    previous_fit_exists <-
      !is.na(previous_index) &&
      file.exists(first_remedy_manifest$fit_path[[previous_index]])
    structural_preserved <-
      key != "function::23" ||
      (isTRUE(specification$uncorrelated_study_effects) &&
         input$stan_data$correlated_study_effects == 0L)
    valid <-
      identical(input$family, "rescor") &&
      identical(
        specification$update_2026_rescor_remedy_level,
        "remedy_4"
      ) &&
      identical(
        specification$update_2026_rescor_remedy,
        "smaller_step_size_after_isolated_divergence"
      ) &&
      identical(
        specification$update_2026_rescor_remedy_stage,
        "remedy_2"
      ) &&
      input$stan_data$N > 0L && input$stan_data$N_study > 1L &&
      previous_fit_exists && structural_preserved &&
      local_config$iter_warmup == 6000L &&
      local_config$iter_sampling == 12000L &&
      local_config$adapt_delta_floor >= 0.999999 &&
      local_config$max_treedepth_floor >= 18L

    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      historical_schedule =
        specification$update_2026_rescor_historical_schedule,
      remedy_level = specification$update_2026_rescor_remedy_level,
      remedy = specification$update_2026_rescor_remedy,
      uncorrelated_study_effects =
        isTRUE(specification$uncorrelated_study_effects),
      adapt_delta = max(
        as.numeric(specification$adapt_delta),
        as.numeric(local_config$adapt_delta_floor)
      ),
      max_treedepth = max(
        as.integer(specification$max_treedepth),
        as.integer(local_config$max_treedepth_floor)
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      previous_fit_exists = previous_fit_exists,
      valid = valid
    )
  }) |>
    dplyr::arrange(model_index)
  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  if (
    nrow(first_remedy_manifest) != 22L ||
      anyDuplicated(first_remedy_keys) ||
      nrow(checks) != 2L || anyDuplicated(keys) ||
      !setequal(keys, expected_keys) || !all(checks$valid)
  ) {
    stop(
      "The second residual-correlation remedy preflight is incomplete.",
      call. = FALSE
    )
  }
  checks
}

write_update_2026_rescor_second_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "remedy_2", "input_preflight.csv"
    )
) {
  if (nrow(checks) != 2L || !all(checks$valid)) {
    stop(
      "Cannot write an incomplete second residual-correlation remedy preflight.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

write_update_2026_rescor_second_remedy_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "remedy_2", "fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 2L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$fit_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop(
      "Cannot write an incomplete second residual-correlation remedy manifest.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

write_update_2026_rescor_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "accepted_fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$fit_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an incomplete residual-correlation accepted manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

make_update_2026_rescor_diagnostic_acceptance <- function(
    final_audit,
    second_remedy_preflight,
    second_remedy_manifest
) {
  summary <- final_audit$diagnostic_summary
  keys <- paste(summary$outcome, summary$model_index, sep = "::")
  failed <- summary |>
    dplyr::filter(!.data$diagnostic_pass)
  exception_key <- "function::23"
  exception_row <- failed |>
    dplyr::filter(
      .data$outcome == "function",
      .data$model_index == 23L
    )
  preflight_row <- second_remedy_preflight |>
    dplyr::filter(
      .data$outcome == "function",
      .data$model_index == 23L
    )
  manifest_row <- second_remedy_manifest |>
    dplyr::filter(
      .data$outcome == "function",
      .data$model_index == 23L
    )
  posterior_draws <- 4L * preflight_row$iter_sampling[[1L]]
  exception_valid <-
    nrow(exception_row) == 1L &&
    exception_row$divergent_transitions[[1L]] == 1L &&
    exception_row$maximum_treedepth_transitions[[1L]] == 0L &&
    exception_row$minimum_ebfmi[[1L]] >= 0.2 &&
    exception_row$maximum_rhat[[1L]] <= 1.01 &&
    exception_row$minimum_bulk_ess[[1L]] >= 400 &&
    exception_row$minimum_tail_ess[[1L]] >= 400 &&
    exception_row$maximum_mcse_sd_ratio[[1L]] <= 0.10 &&
    nrow(preflight_row) == 1L && preflight_row$valid[[1L]] &&
    identical(preflight_row$remedy_level[[1L]], "remedy_4") &&
    isTRUE(preflight_row$uncorrelated_study_effects[[1L]]) &&
    preflight_row$adapt_delta[[1L]] >= 0.999999 &&
    preflight_row$max_treedepth[[1L]] >= 18L &&
    preflight_row$iter_warmup[[1L]] == 6000L &&
    preflight_row$iter_sampling[[1L]] == 12000L &&
    nrow(manifest_row) == 1L &&
    identical(manifest_row$stage[[1L]], "remedy_2") &&
    file.exists(manifest_row$fit_path[[1L]]) &&
    manifest_row$bytes[[1L]] > 0 &&
    nchar(manifest_row$sha256[[1L]]) == 64L &&
    posterior_draws == 48000L
  if (
    !identical(final_audit$family, "rescor") ||
      nrow(summary) != 104L || anyDuplicated(keys) ||
      final_audit$models_passing_sampling != 103L ||
      nrow(failed) != 1L || !identical(keys[!summary$diagnostic_pass],
                                      exception_key) ||
      !isTRUE(exception_valid)
  ) {
    stop(
      "The residual-correlation diagnostic exception is not auditable.",
      call. = FALSE
    )
  }

  accepted_summary <- summary |>
    dplyr::mutate(
      raw_diagnostic_pass = .data$diagnostic_pass,
      documented_exception =
        paste(.data$outcome, .data$model_index, sep = "::") ==
          exception_key,
      accepted_diagnostic_pass =
        .data$diagnostic_pass | .data$documented_exception,
      acceptance_status = dplyr::case_when(
        .data$diagnostic_pass ~ "PASS",
        .data$documented_exception ~
          "ACCEPTED WITH DOCUMENTED SINGLE-DIVERGENCE EXCEPTION",
        TRUE ~ "REQUIRES REMEDIATION"
      )
    )
  exception_manifest <- tibble::tibble(
    family = "rescor",
    outcome = "function",
    model_index = 23L,
    model_name = exception_row$model_name[[1L]],
    exception = "single_divergence_after_exhaustive_remediation",
    divergent_transitions = 1L,
    posterior_draws = posterior_draws,
    divergence_rate = 1 / posterior_draws,
    maximum_treedepth_transitions = 0L,
    minimum_ebfmi = exception_row$minimum_ebfmi[[1L]],
    maximum_rhat = exception_row$maximum_rhat[[1L]],
    minimum_bulk_ess = exception_row$minimum_bulk_ess[[1L]],
    minimum_tail_ess = exception_row$minimum_tail_ess[[1L]],
    maximum_mcse_sd_ratio =
      exception_row$maximum_mcse_sd_ratio[[1L]],
    structural_remedy = "independent_study_random_smooth_coefficients",
    adapt_delta = preflight_row$adapt_delta[[1L]],
    max_treedepth = preflight_row$max_treedepth[[1L]],
    iter_warmup = preflight_row$iter_warmup[[1L]],
    iter_sampling = preflight_row$iter_sampling[[1L]],
    fit_path = manifest_row$fit_path[[1L]],
    fit_sha256 = manifest_row$sha256[[1L]],
    rationale = paste(
      "One non-boundary divergence remained among 48,000 posterior draws",
      "after structural decorrelation and two increasingly stringent",
      "sampling remedies; all other geometry and convergence checks passed."
    )
  )
  if (!all(accepted_summary$accepted_diagnostic_pass)) {
    stop("The residual-correlation acceptance audit is incomplete.",
         call. = FALSE)
  }
  list(
    family = "rescor",
    diagnostic_summary = accepted_summary,
    exception_manifest = exception_manifest,
    models = 104L,
    models_raw_passing = sum(accepted_summary$raw_diagnostic_pass),
    models_with_documented_exception =
      sum(accepted_summary$documented_exception),
    models_accepted = sum(accepted_summary$accepted_diagnostic_pass),
    all_raw_sampling_diagnostics_pass =
      all(accepted_summary$raw_diagnostic_pass),
    all_sampling_diagnostics_accepted =
      all(accepted_summary$accepted_diagnostic_pass)
  )
}

write_update_2026_rescor_diagnostic_acceptance <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "acceptance"
    )
) {
  if (
    !identical(audit$family, "rescor") ||
      nrow(audit$diagnostic_summary) != 104L ||
      nrow(audit$exception_manifest) != 1L ||
      !isTRUE(audit$all_sampling_diagnostics_accepted)
  ) {
    stop("Cannot write an invalid residual-correlation acceptance audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(
      output_dir, "diagnostic_acceptance_summary.csv"
    ),
    documented_exceptions = file.path(
      output_dir, "documented_exceptions.csv"
    )
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$exception_manifest,
    paths[["documented_exceptions"]]
  )
  unname(paths)
}

combine_update_2026_rescor_accepted_fit_manifest <- function(
    initial_manifest,
    remedy_manifest,
    acceptance_audit
) {
  remedy_keys <- paste(
    remedy_manifest$outcome,
    remedy_manifest$model_index,
    sep = "::"
  )
  accepted <- initial_manifest |>
    dplyr::filter(
      !paste(.data$outcome, .data$model_index, sep = "::") %in%
        remedy_keys
    ) |>
    dplyr::bind_rows(remedy_manifest) |>
    dplyr::arrange(
      factor(
        .data$outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      .data$model_index
    )
  if (
    !isTRUE(acceptance_audit$all_sampling_diagnostics_accepted) ||
      acceptance_audit$models_raw_passing != 103L ||
      acceptance_audit$models_with_documented_exception != 1L ||
      nrow(accepted) != 104L || any(accepted$family != "rescor") ||
      anyDuplicated(paste(accepted$outcome, accepted$model_index)) ||
      any(!file.exists(accepted$fit_path)) ||
      any(accepted$bytes <= 0) || any(nchar(accepted$sha256) != 64L)
  ) {
    stop("The accepted residual-correlation manifest is invalid.",
         call. = FALSE)
  }
  accepted
}

validate_update_2026_rescor_loo_inputs <- function(
    accepted_manifest,
    inputs,
    acceptance_audit
) {
  manifest_keys <- paste(
    accepted_manifest$outcome,
    accepted_manifest$model_index,
    sep = "::"
  )
  input_keys <- purrr::map_chr(inputs, function(input) {
    paste(
      input$specification$outcome,
      input$specification$model_index,
      sep = "::"
    )
  })
  checks <- accepted_manifest |>
    dplyr::transmute(
      family,
      outcome,
      model_index = as.integer(model_index),
      model_name,
      stage,
      fit_path,
      bytes,
      sha256,
      fit_exists = file.exists(fit_path),
      valid =
        .data$family == "rescor" & .data$fit_exists &
        .data$bytes > 0 & nchar(.data$sha256) == 64L
    ) |>
    dplyr::arrange(
      factor(
        .data$outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      .data$model_index
    )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  stage_counts <- table(factor(
    checks$stage,
    levels = c("initial", "remedy_1", "remedy_2")
  ))
  if (
    !isTRUE(acceptance_audit$all_sampling_diagnostics_accepted) ||
      acceptance_audit$models_accepted != 104L ||
      nrow(checks) != 104L || anyDuplicated(manifest_keys) ||
      length(inputs) != 104L || anyDuplicated(input_keys) ||
      !setequal(manifest_keys, input_keys) || !all(checks$valid) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !identical(as.integer(stage_counts), c(82L, 20L, 2L))
  ) {
    stop("The residual-correlation LOO preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_rescor_loo_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete residual-correlation LOO preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

compute_update_2026_rescor_integrated_loo_to_file <- function(
    fit_record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "rescor", "loo",
      "study_integrated"
    )
) {
  result <- compute_rescor_full_file_integrated_loo(
    record = fit_record,
    loo_draws = loo_draws,
    k_threshold = k_threshold,
    cores = cores
  )
  output_dir <- file.path(
    output_root,
    result$outcome,
    sprintf("model_%02d", as.integer(result$model_index))
  )
  output_file <- file.path(output_dir, "loo.qs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_atomic_qs(result, output_file)
}

write_update_2026_rescor_integrated_loo_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo",
      "study_integrated_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$integrated_loo_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an incomplete residual-correlation LOO manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

write_update_2026_rescor_integrated_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo",
      "study_integrated"
    )
) {
  method_counts <- table(factor(
    audit$model_summary$method,
    levels = c(
      "conditional_study_cluster_psis",
      "integrated_study_cluster_psis"
    )
  ))
  if (
    !identical(audit$family, "rescor") || audit$models != 104L ||
      nrow(audit$model_summary) != 208L ||
      !identical(as.integer(method_counts), c(104L, 104L)) ||
      nrow(audit$outcome_summary) != 4L ||
      audit$exact_refits_required != nrow(audit$high_k_cases)
  ) {
    stop("Cannot write an invalid residual-correlation LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

compute_update_2026_rescor_observation_loo_to_file <- function(
    fit_record,
    generated_quantities_model,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    parallel_chains = 4L,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "rescor", "loo",
      "observation_conditional"
    )
) {
  result <- compute_rescor_full_file_observation_loo(
    record = fit_record,
    generated_quantities_model = generated_quantities_model,
    loo_draws = loo_draws,
    k_threshold = k_threshold,
    cores = cores,
    parallel_chains = parallel_chains
  )
  if (
    !inherits(result$loo, "loo") ||
      !identical(result$family, "rescor") ||
      !identical(
        result$loo_unit,
        "observation_conditional_on_study"
      ) ||
      nrow(result$observation_manifest) != nrow(result$loo$pointwise)
  ) {
    stop("Cannot archive an invalid updated observation-level LOO result.",
         call. = FALSE)
  }
  output_dir <- file.path(
    output_root,
    result$outcome,
    sprintf("model_%02d", as.integer(result$model_index))
  )
  output_file <- file.path(output_dir, "loo.qs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_atomic_qs(result, output_file)
}

write_update_2026_rescor_observation_loo_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo",
      "observation_conditional_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      any(manifest$family != "rescor") ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(!file.exists(manifest$integrated_loo_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an incomplete observation-level LOO manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

write_update_2026_rescor_observation_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo",
      "observation_conditional"
    )
) {
  if (
    audit$models != 104L ||
      !identical(audit$loo_unit, "observation_conditional_on_study") ||
      !isTRUE(audit$observation_psis_complete) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      audit$exact_refits_required != nrow(audit$high_k_cases)
  ) {
    stop("Cannot write an invalid observation-level LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

validate_update_2026_rescor_observation_exact_inputs <- function(
    records,
    observation_audit,
    config
) {
  expected <- as.integer(observation_audit$exact_refits_required)
  if (!is.list(records) || length(records) != expected || expected != 46L) {
    stop("Expected 46 residual-correlation observation exact records.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    original_exists <- file.exists(record$original_fit_path)
    observation_loo_exists <- file.exists(record$observation_loo_path)
    observation_loo <- if (observation_loo_exists) {
      qs::qread(record$observation_loo_path)
    } else {
      NULL
    }
    observation <- as.integer(record$observation)
    pareto_values <- if (
      is.list(observation_loo) && inherits(observation_loo$loo, "loo")
    ) {
      as.numeric(loo::pareto_k_values(observation_loo$loo))
    } else {
      numeric()
    }
    manifest_row <- if (
      is.list(observation_loo) &&
        is.data.frame(observation_loo$observation_manifest) &&
        observation >= 1L &&
        observation <= nrow(observation_loo$observation_manifest)
    ) {
      observation_loo$observation_manifest[observation, , drop = FALSE]
    } else {
      NULL
    }
    deletion_valid <- FALSE
    if (
      is.list(observation_loo) && is.list(observation_loo$input) &&
        is.list(observation_loo$input$stan_data) &&
        observation >= 1L &&
        observation <= observation_loo$input$stan_data$N
    ) {
      training_data <- make_rescor_leave_one_observation_out_data(
        observation_loo$input$stan_data,
        observation
      )
      deletion_valid <-
        training_data$N == observation_loo$input$stan_data$N - 1L &&
        training_data$N_cluster >= 1L &&
        length(training_data$cluster_start) == training_data$N_cluster &&
        length(training_data$cluster_end) == training_data$N_cluster
    }
    specification <- if (
      is.list(observation_loo) && is.list(observation_loo$input)
    ) {
      observation_loo$input$specification
    } else {
      NULL
    }
    actual_adapt_delta <- if (is.list(specification)) {
      max(
        as.numeric(specification$adapt_delta),
        as.numeric(config$adapt_delta_floor)
      )
    } else {
      NA_real_
    }
    actual_max_treedepth <- if (is.list(specification)) {
      max(
        as.integer(specification$max_treedepth),
        as.integer(config$max_treedepth_floor)
      )
    } else {
      NA_integer_
    }
    valid <-
      original_exists && observation_loo_exists &&
      nchar(record$original_fit_sha256) == 64L &&
      nchar(record$observation_loo_sha256) == 64L &&
      is.list(observation_loo) &&
      identical(observation_loo$family, "rescor") &&
      identical(
        observation_loo$loo_unit,
        "observation_conditional_on_study"
      ) &&
      identical(observation_loo$outcome, record$outcome) &&
      as.integer(observation_loo$model_index) == record$model_index &&
      !is.null(manifest_row) &&
      observation >= 1L && observation <= length(pareto_values) &&
      is.finite(pareto_values[[observation]]) &&
      pareto_values[[observation]] > 0.7 &&
      isTRUE(all.equal(
        pareto_values[[observation]],
        as.numeric(record$pareto_k),
        tolerance = 1e-10
      )) &&
      is.finite(actual_adapt_delta) &&
      !is.na(actual_max_treedepth) && deletion_valid

    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      observation = observation,
      study_no = if (is.null(manifest_row)) {
        NA_character_
      } else {
        as.character(manifest_row$study_no[[1L]])
      },
      study_cluster = if (is.null(manifest_row)) {
        NA_integer_
      } else {
        as.integer(manifest_row$study_cluster[[1L]])
      },
      cluster_size = if (is.null(manifest_row)) {
        NA_integer_
      } else {
        as.integer(manifest_row$cluster_size[[1L]])
      },
      pareto_k = as.numeric(record$pareto_k),
      adapt_delta = actual_adapt_delta,
      max_treedepth = actual_max_treedepth,
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      original_fit_exists = original_exists,
      observation_loo_exists = observation_loo_exists,
      valid = valid
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- paste(
    observation_audit$high_k_cases$outcome,
    observation_audit$high_k_cases$model_index,
    observation_audit$high_k_cases$observation,
    sep = "::"
  )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  if (
    nrow(checks) != expected || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(expected)) ||
      !identical(as.integer(outcome_counts), c(1L, 42L, 3L, 0L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$observation_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.999 ||
      config$max_treedepth_floor != 15L
  ) {
    stop("The observation exact-refit preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_rescor_observation_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo_exact",
      "observation_input_preflight.csv"
    )
) {
  if (nrow(checks) != 46L || !all(checks$valid)) {
    stop("Cannot write an incomplete observation exact-refit preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_rescor_observation_exact_result_file_path <- function(
    record,
    stage = "observation_exact"
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "rescor", "loo",
    stage,
    record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

audit_update_2026_rescor_observation_exact_files <- function(
    paths,
    records,
    observation_audit
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(observation_audit$exact_refits_required)
  if (
    length(paths) != expected || length(records) != expected ||
      expected != 46L || anyNA(paths) || any(!file.exists(paths)) ||
      any(file.info(paths)$size <= 0)
  ) {
    stop("The observation exact-refit files are incomplete.",
         call. = FALSE)
  }
  results <- purrr::map(paths, read_rescor_observation_exact_result_file)
  diagnostic_summary <- purrr::map_dfr(results, "diagnostic") |>
    dplyr::arrange(case_id)
  manifest <- purrr::map2_dfr(paths, results, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      observation = as.integer(result$observation),
      stage = result$stage,
      correction_available = is.list(result$correction),
      result_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation,
    sep = "::"
  )
  expected_keys <- paste(
    observation_audit$high_k_cases$outcome,
    observation_audit$high_k_cases$model_index,
    observation_audit$high_k_cases$observation,
    sep = "::"
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(manifest$case_id) ||
      !identical(manifest$case_id, seq_len(expected)) ||
      anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys) ||
      any(nchar(manifest$sha256) != 64L) ||
      any(manifest$correction_available != diagnostic_summary$diagnostic_pass)
  ) {
    stop("The observation exact-refit audit is invalid.", call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    result_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = sum(diagnostic_summary$diagnostic_pass),
    exact_refits_requiring_remediation =
      sum(!diagnostic_summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass =
      all(diagnostic_summary$diagnostic_pass)
  )
}

write_update_2026_rescor_observation_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo_exact",
      "observation_initial"
    )
) {
  if (
    audit$exact_refits != 46L ||
      nrow(audit$diagnostic_summary) != 46L ||
      nrow(audit$result_manifest) != 46L
  ) {
    stop("Cannot write an incomplete observation exact-refit audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    result_manifest = file.path(output_dir, "result_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(audit$result_manifest, paths[["result_manifest"]])
  unname(paths)
}

validate_update_2026_rescor_observation_exact_remedy_inputs <- function(
    records,
    initial_audit,
    config
) {
  failed <- initial_audit$diagnostic_summary |>
    dplyr::filter(!diagnostic_pass) |>
    dplyr::select(case_id, outcome, model_index, observation)
  if (
    !is.list(records) || length(records) != nrow(failed) ||
      nrow(failed) != 1L
  ) {
    stop("Expected one observation exact-refit remedy record.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    initial_row <- initial_audit$result_manifest |>
      dplyr::filter(case_id == as.integer(record$case_id))
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      observation = as.integer(record$observation),
      remedy = "adapt_delta_0.9999_after_observation_exact_refit",
      adapt_delta = as.numeric(config$adapt_delta_floor),
      max_treedepth = as.integer(config$max_treedepth_floor),
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      initial_result_exists =
        nrow(initial_row) == 1L && file.exists(initial_row$result_path),
      original_fit_exists = file.exists(record$original_fit_path),
      observation_loo_exists = file.exists(record$observation_loo_path)
    ) |>
      dplyr::mutate(
        valid =
          .data$initial_result_exists & .data$original_fit_exists &
          .data$observation_loo_exists
      )
  }) |>
    dplyr::arrange(case_id)
  observed <- checks |>
    dplyr::select(case_id, outcome, model_index, observation)
  failed <- dplyr::arrange(failed, case_id)
  if (
    !identical(observed, failed) ||
      checks$case_id[[1L]] != 20L ||
      checks$outcome[[1L]] != "hypertrophy" ||
      checks$model_index[[1L]] != 17L ||
      checks$observation[[1L]] != 252L ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 3000L || config$iter_sampling != 8000L ||
      config$adapt_delta_floor != 0.9999 ||
      config$max_treedepth_floor != 16L || !all(checks$valid)
  ) {
    stop("The observation exact-refit remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_rescor_observation_exact_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo_exact",
      "observation_remedy_preflight.csv"
    )
) {
  if (nrow(checks) != 1L || !all(checks$valid)) {
    stop("Cannot write an incomplete observation exact-remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

make_update_2026_rescor_final_observation_loo_audit <- function(
    observation_loo_manifest,
    exact_corrections,
    exact_audit
) {
  if (
    !isTRUE(exact_audit$all_exact_refit_diagnostics_pass) ||
      exact_audit$exact_refits != 46L || length(exact_corrections) != 46L
  ) {
    stop("Final observation LOO requires 46 passing exact diagnostics.",
         call. = FALSE)
  }
  final <- apply_rescor_observation_exact_corrections(
    observation_loo_manifest = observation_loo_manifest,
    exact_results = exact_corrections
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- observation_loo_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "correlated_conditional_observation_level_psis_plus_exact",
        "correlated_conditional_observation_level_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(final$exact_refit_summary) != 46L ||
      nrow(outcome_summary) != 4L || any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 46L
  ) {
    stop("The final residual-correlation observation LOO is not ready.",
         call. = FALSE)
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 46L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_rescor_final_observation_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      !isTRUE(audit$all_exact_refit_diagnostics_pass) ||
      audit$exact_refits != 46L ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 46L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop("Cannot write an incomplete final observation-level LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_rescor_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L || any(weights$family != "rescor") ||
      anyDuplicated(paste(weights$outcome, weights$model_index)) ||
      length(outcome_sums) != 4L ||
      any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop("The updated residual-correlation stacking weights are invalid.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_rescor_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    inputs,
    n_draws = 4000L,
    seed = 34988L
) {
  records <- make_rescor_full_prediction_records(
    accepted_fit_manifest = accepted_fit_manifest,
    stacking_weights = stacking_weights,
    n_draws = n_draws,
    seed = seed
  )
  input_keys <- purrr::map_chr(inputs, function(input) {
    paste(
      input$specification$outcome,
      input$specification$model_index,
      sep = "::"
    )
  })
  records <- purrr::map(records, function(record) {
    key <- paste(record$outcome, record$model_index, sep = "::")
    position <- match(key, input_keys)
    if (is.na(position)) {
      stop("An updated residual-correlation prediction input is unavailable.",
           call. = FALSE)
    }
    record$protein_max <- max(
      inputs[[position]]$grid$protein_intake_raw
    )
    record
  })
  if (
    length(inputs) != 104L || anyDuplicated(input_keys) ||
      length(records) < 4L ||
      sum(purrr::map_int(records, "draw_count")) != 16000L
  ) {
    stop("The residual-correlation prediction allocation is incomplete.",
         call. = FALSE)
  }
  records
}

make_update_2026_rescor_reporting_preflight <- function(records) {
  if (!is.list(records) || length(records) < 4L) {
    stop("The residual-correlation reporting allocation is incomplete.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    tibble::tibble(
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      family = unname(as.character(record$family)),
      prediction_type = unname(as.character(record$prediction_type)),
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      protein_max = as.numeric(record$protein_max),
      grid_step = as.numeric(record$grid_step),
      fit_path = record$fit_path,
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path) &&
        identical(unname(as.character(record$family)), "rescor") &&
        identical(unname(as.character(record$prediction_type)), "rescor") &&
        nchar(record$fit_sha256) == 64L &&
        is.finite(record$protein_min) && is.finite(record$protein_max) &&
        record$protein_max > record$protein_min
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      retained_weight = sum(weight),
      common_protein_max = dplyr::first(protein_max),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      any(!is.finite(outcome_audit$retained_weight)) ||
      any(outcome_audit$retained_weight <= 0) ||
      any(outcome_audit$retained_weight > 1 + 1e-8) ||
      any(!is.finite(outcome_audit$common_protein_max)) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop("The residual-correlation reporting preflight failed.",
         call. = FALSE)
  }
  checks
}

write_update_2026_rescor_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "rescor", "reporting",
      "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(checks$allocated_draws, checks$outcome, sum)
  outcome_weights <- tapply(checks$weight, checks$outcome, sum)
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      length(outcome_weights) != 4L ||
      any(!is.finite(outcome_weights)) || any(outcome_weights <= 0) ||
      any(outcome_weights > 1 + 1e-8) || !all(checks$valid)
  ) {
    stop("Cannot write an invalid residual-correlation reporting preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_rescor_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "rescor",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_rescor_prediction_draw_file <- function(record) {
  write_rescor_full_prediction_draw_file(
    record = record,
    checkpoint_path = NULL,
    result_path = update_2026_rescor_reporting_path(
      outcome = record$outcome,
      stage = "prediction_draws",
      model_index = record$model_index
    )
  )
}

write_update_2026_rescor_model_average_file <- function(record) {
  write_sensitivity_full_model_average_file(
    record = record,
    result_path = update_2026_rescor_reporting_path(
      outcome = record$outcome,
      stage = "model_average"
    )
  )
}

write_update_2026_rescor_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path,
    family = "rescor"
  )
  write_sensitivity_full_reporting_file(
    model_average_path = model_average_path,
    data_by_outcome = data_by_outcome,
    family = "rescor",
    reference_intake = reference_intake,
    sesoi = sesoi,
    result_path = update_2026_rescor_reporting_path(
      outcome = model_average$outcome,
      stage = "reporting"
    )
  )
}

write_update_2026_rescor_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "reporting"
    )
) {
  if (
    !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L
  ) {
    stop("The residual-correlation reporting audit is incomplete.",
         call. = FALSE)
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "rescor"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir,
      "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir,
      "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir,
      "common_exposure_support.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_rescor_verification_checkpoint <- function(
    data_by_outcome,
    input_preflight,
    source_manifest,
    review_decisions,
    preparation_audit,
    specifications,
    accepted_fit_manifest,
    diagnostic_acceptance,
    study_loo_audit,
    observation_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  ) |>
    dplyr::filter(
      is.finite(yi), is.finite(vi), vi > 0,
      is.finite(protein_intake),
      is.finite(resistance_exercise_code)
    )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "rescor"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(data) {
    data |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_acceptance$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_acceptance$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_name <- paste0(
    name_alignment$specification_name,
    " [residual correlation by study]"
  )
  sample_preflight <- input_preflight |>
    dplyr::distinct(
      outcome, rows, studies, arms, outcome_rows,
      minimum_cluster_size, maximum_cluster_size,
      rho_lower, rho_prior_median,
      protein_minimum, protein_maximum
    )
  analysis_row_counts <- analysis_data |>
    dplyr::count(analysis_outcome, name = "analysis_rows")
  sample_alignment <- sample_preflight |>
    dplyr::select(outcome, preflight_rows = rows) |>
    dplyr::left_join(
      analysis_row_counts,
      by = c("outcome" = "analysis_outcome")
    )
  weight_sums <- stacking_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      .by = outcome
    )
  reporting_model_weights <- bind_report("model_weights")
  reporting_weight_sums <- reporting_model_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "all four residual-correlation datasets and covariance inputs are valid",
      "the combined-data preparation audit remains valid",
      "all five residual-correlation registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist and have checksums",
      "all 104 sampling diagnostics are accepted with one documented exception",
      "study-level LOO retains 104 models and 152 audited high-k cases",
      "conditional observation-level LOO contains 104 models and 46 high-k cases",
      "all 46 final exact refits pass and compact results exist",
      "final corrected observation LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four reporting audits pass and retain all model weights"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      identical(names(data_by_outcome), outcomes) &&
        nrow(input_preflight) == 104L && all(input_preflight$valid) &&
        nrow(sample_preflight) == 4L &&
        setequal(sample_preflight$outcome, outcomes) &&
        nrow(sample_alignment) == 4L &&
        all(sample_alignment$preflight_rows ==
              sample_alignment$analysis_rows) &&
        all(sample_preflight$studies > 1L) &&
        all(sample_preflight$maximum_cluster_size > 1L) &&
        all(sample_preflight$rho_lower > -1) &&
        all(sample_preflight$rho_lower < 0) &&
        all(abs(sample_preflight$rho_prior_median) < 1e-10),
      all(preparation_audit$checks$passed),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_name) &&
        all(name_alignment$diagnostic_name == expected_name) &&
        all(name_alignment$loo_name == expected_name) &&
        all(name_alignment$weight_name == expected_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(accepted_fit_manifest$family == "rescor") &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(accepted_fit_manifest$bytes > 0) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      diagnostic_acceptance$models == 104L &&
        diagnostic_acceptance$models_raw_passing == 103L &&
        diagnostic_acceptance$models_with_documented_exception == 1L &&
        diagnostic_acceptance$models_accepted == 104L &&
        isTRUE(diagnostic_acceptance$all_sampling_diagnostics_accepted) &&
        nrow(diagnostic_acceptance$exception_manifest) == 1L &&
        all(diagnostic_acceptance$diagnostic_summary$accepted_diagnostic_pass),
      study_loo_audit$models == 104L &&
        identical(study_loo_audit$loo_unit, "study") &&
        study_loo_audit$exact_refits_required == 152L &&
        nrow(study_loo_audit$high_k_cases) == 152L &&
        isTRUE(study_loo_audit$integrated_psis_complete),
      observation_loo_audit$models == 104L &&
        identical(
          observation_loo_audit$loo_unit,
          "observation_conditional_on_study"
        ) &&
        observation_loo_audit$exact_refits_required == 46L &&
        nrow(observation_loo_audit$high_k_cases) == 46L &&
        isTRUE(observation_loo_audit$observation_psis_complete),
      exact_audit$exact_refits == 46L &&
        exact_audit$exact_refits_passing == 46L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        nrow(exact_audit$result_manifest) == 46L &&
        all(exact_audit$result_manifest$correction_available) &&
        all(file.exists(exact_audit$result_manifest$result_path)) &&
        all(exact_audit$result_manifest$bytes > 0) &&
        all(nchar(exact_audit$result_manifest$sha256) == 64L),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        final_loo_audit$exact_refits == 46L &&
        isTRUE(final_loo_audit$all_exact_refit_diagnostics_pass) &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(file.exists(prediction_manifest$prediction_path)) &&
        all(prediction_manifest$bytes > 0) &&
        all(nchar(prediction_manifest$sha256) == 64L),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L &&
        nrow(reporting_model_weights) == 104L &&
        nrow(reporting_weight_sums) == 4L &&
        all(reporting_weight_sums$models == 26L) &&
        all(abs(reporting_weight_sums$weight_sum - 1) < 1e-8)
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated residual-correlation verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    review_decisions = review_decisions,
    preparation_checks = preparation_audit$checks,
    preparation_outcome_counts = preparation_audit$outcome_counts,
    preparation_update_study_counts = preparation_audit$update_study_counts,
    rescor_input_preflight = input_preflight,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_acceptance$diagnostic_summary,
    sampling_exception_manifest = diagnostic_acceptance$exception_manifest,
    study_loo_model_summary = study_loo_audit$model_summary,
    study_loo_outcome_summary = study_loo_audit$outcome_summary,
    study_high_pareto_k_cases = study_loo_audit$high_k_cases,
    observation_loo_model_summary = observation_loo_audit$model_summary,
    observation_loo_outcome_summary = observation_loo_audit$outcome_summary,
    observation_high_pareto_k_cases = observation_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_result_manifest = exact_audit$result_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    final_exact_refit_summary = final_loo_audit$exact_refit_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    reporting_model_weights = reporting_model_weights,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final updated residual-correlation checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    review_decisions = "Eligibility decisions for updated-search studies.",
    preparation_checks = "Combined-data preparation invariants.",
    preparation_outcome_counts = "Analysis counts by primary outcome.",
    preparation_update_study_counts = "Updated-study contribution audit.",
    rescor_input_preflight = "Study-block covariance and schedule checks.",
    analysis_data = "Eligible combined rows supplied to rescor models.",
    model_specifications = "Final 104-model rescor specification registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit path, stage, size, and checksum.",
    sampling_diagnostic_summary = "Final accepted sampling diagnostics.",
    sampling_exception_manifest = "Documented one-divergence exception.",
    study_loo_model_summary = "Study-level integrated LOO model summaries.",
    study_loo_outcome_summary = "Study-level influence audit by outcome.",
    study_high_pareto_k_cases = "The 152 study-level high-k influence cases.",
    observation_loo_model_summary = "Conditional observation LOO summaries.",
    observation_loo_outcome_summary = "Initial observation LOO by outcome.",
    observation_high_pareto_k_cases = "The 46 exact-refit observations.",
    exact_refit_diagnostic_summary = "Final diagnostics for 46 exact refits.",
    exact_result_manifest = "Compact exact-result paths and checksums.",
    final_loo_manifest = "Final exact-corrected observation LOO lineage.",
    final_loo_model_summary = "Final observation LOO for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    final_exact_refit_summary = "Final exact-correction audit by case.",
    stacking_weights = "All 104 final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    reporting_model_weights = "All 104 weights retained in reports.",
    dose_response_predictions = "Rescor model-averaged summaries.",
    dose_response_increments = "Increments from the 0.8 g/kg/day reference.",
    exposure_support_by_group = "Observed exposure support by RT group.",
    common_exposure_support = "Common exposure support by outcome.",
    pairwise_contrasts = "Model-averaged contrasts and probabilities.",
    pairwise_support = "Support limits for pairwise contrasts.",
    local_equivalence = "Local practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    next_family = paste(
      "The residual-correlation family is frozen after this verification",
      "checkpoint; its lean-pruning plan can now be audited."
    )
  )
}

make_update_2026_rescor_pruning_plan <- function(
    accepted_fit_manifest,
    verification_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "rescor"
    )
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "rescor_update_2026"),
    "manifest.csv"
  )
  if (
    length(verification_files) != 38L ||
      any(!file.exists(verification_files))
  ) {
    stop("The updated residual-correlation verification bundle is incomplete.",
         call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "rescor"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated rescor archive.",
         call. = FALSE)
  }

  fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  if (any(!dir.exists(fit_roots))) {
    stop("An expected updated rescor fit directory is missing.",
         call. = FALSE)
  }
  deletion_paths <- unlist(
    purrr::map(
      fit_roots,
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  deletion_paths <- deletion_paths[
    file.exists(deletion_paths) & !file.info(deletion_paths)$isdir
  ]
  deletion_paths <- normalizePath(
    deletion_paths,
    winslash = "/",
    mustWork = TRUE
  )
  fit_paths <- deletion_paths[basename(deletion_paths) == "fit.qs"]
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(deletion_paths) != 644L || length(fit_paths) != 128L ||
      length(accepted_paths) != 104L ||
      anyDuplicated(deletion_paths) || anyDuplicated(fit_paths) ||
      anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% fit_paths) ||
      sum(!fit_paths %in% accepted_paths) != 24L
  ) {
    stop("Updated rescor fit retention does not match its lineage.",
         call. = FALSE)
  }

  fit_manifest <- tibble::tibble(
    path = fit_paths,
    relative_path = substring(fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", fit_paths, fixed = TRUE) ~ "initial",
      grepl("/remedies/remedy_1/", fit_paths, fixed = TRUE) ~ "remedy_1",
      grepl("/remedies/remedy_2/", fit_paths, fixed = TRUE) ~ "remedy_2",
      TRUE ~ "unknown"
    ),
    outcome = basename(dirname(dirname(fit_paths))),
    model_index = as.integer(sub(
      "^model_", "", basename(dirname(fit_paths))
    )),
    accepted_fit = fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(fit_paths)$size),
    lean_action = "delete frozen residual-correlation fitted object"
  ) |>
    dplyr::arrange(archive_stage, outcome, model_index)
  stage_levels <- c("initial", "remedy_1", "remedy_2")
  if (
    any(fit_manifest$archive_stage == "unknown") ||
      anyNA(fit_manifest$model_index) ||
      !setequal(fit_manifest$outcome, c(
        "whole_body", "hypertrophy", "strength", "function"
      )) ||
      !identical(
        as.integer(table(factor(
          fit_manifest$archive_stage,
          levels = stage_levels
        ))),
        c(104L, 22L, 2L)
      ) ||
      !identical(
        as.integer(table(factor(
          fit_manifest$archive_stage[fit_manifest$accepted_fit],
          levels = stage_levels
        ))),
        c(82L, 20L, 2L)
      )
  ) {
    stop("The updated rescor fit stages are not as expected.",
         call. = FALSE)
  }

  fit_directories <- dirname(fit_manifest$path)
  file_model_directories <- vapply(
    deletion_paths,
    function(path) {
      if (basename(path) == "fit.qs") dirname(path) else dirname(dirname(path))
    },
    character(1L),
    USE.NAMES = FALSE
  )
  fit_positions <- match(file_model_directories, fit_directories)
  file_types <- dplyr::case_when(
    basename(deletion_paths) == "fit.qs" ~ "fit_object",
    tolower(tools::file_ext(deletion_paths)) == "csv" &
      basename(dirname(deletion_paths)) == "cmdstan" ~
        "cmdstan_chain_csv",
    TRUE ~ "unknown"
  )
  deletion_manifest <- tibble::tibble(
    path = deletion_paths,
    relative_path = substring(
      deletion_paths,
      nchar(resolved_root) + 2L
    ),
    archive_stage = fit_manifest$archive_stage[fit_positions],
    outcome = fit_manifest$outcome[fit_positions],
    model_index = fit_manifest$model_index[fit_positions],
    accepted_fit = fit_manifest$accepted_fit[fit_positions],
    file_type = file_types,
    bytes = as.numeric(file.info(deletion_paths)$size),
    lean_action = "delete refittable residual-correlation sampling archive"
  ) |>
    dplyr::arrange(
      archive_stage, outcome, model_index,
      factor(file_type, levels = c("fit_object", "cmdstan_chain_csv")),
      relative_path
    )
  if (
    anyNA(fit_positions) || any(deletion_manifest$file_type == "unknown") ||
      sum(deletion_manifest$file_type == "fit_object") != 128L ||
      sum(deletion_manifest$file_type == "cmdstan_chain_csv") != 516L ||
      !identical(
        as.integer(table(factor(
          deletion_manifest$archive_stage,
          levels = stage_levels
        ))),
        c(520L, 114L, 10L)
      ) ||
      sum(deletion_manifest$accepted_fit) != 524L ||
      sum(!deletion_manifest$accepted_fit) != 120L ||
      any(deletion_manifest$bytes <= 0)
  ) {
    stop("The updated rescor deletion inventory is incomplete.",
         call. = FALSE)
  }

  compact_roots <- file.path(archive_root, c("loo", "reporting"))
  if (any(!dir.exists(compact_roots))) {
    stop("An expected compact rescor archive directory is missing.",
         call. = FALSE)
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[
    file.exists(compact_paths) & !file.info(compact_paths)$isdir
  ]
  compact_paths <- normalizePath(
    compact_paths,
    winslash = "/",
    mustWork = TRUE
  )
  compact_manifest <- tibble::tibble(
    path = compact_paths,
    bytes = as.numeric(file.info(compact_paths)$size),
    sha256 = purrr::map_chr(
      compact_paths,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    ),
    action = "retain compact residual-correlation LOO/reporting archive"
  )
  if (
    nrow(compact_manifest) != 300L ||
      any(compact_manifest$bytes <= 0) ||
      any(nchar(compact_manifest$sha256) != 64L)
  ) {
    stop("The compact rescor archive inventory is incomplete.",
         call. = FALSE)
  }

  summary <- tibble::tibble(
    mode = "lean",
    model_archives_deleted = nrow(fit_manifest),
    accepted_model_archives_deleted = sum(fit_manifest$accepted_fit),
    superseded_model_archives_deleted = sum(!fit_manifest$accepted_fit),
    archive_files_deleted = nrow(deletion_manifest),
    fit_object_files_deleted = sum(
      deletion_manifest$file_type == "fit_object"
    ),
    cmdstan_csv_files_deleted = sum(
      deletion_manifest$file_type == "cmdstan_chain_csv"
    ),
    fit_archive_files_retained = 0L,
    gib_reclaimed = sum(deletion_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained = sum(compact_manifest$bytes) / 1024^3,
    verification_files_retained = length(verification_files),
    restoration = paste(
      "Reported results, both LOO audits, and all 104 weights remain frozen;",
      "restoring full fits and chain CSVs requires refitting the rescor family."
    )
  )

  list(
    summary = summary,
    fit_file_manifest = fit_manifest,
    deletion_file_manifest = deletion_manifest,
    compact_archive_manifest = compact_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_rescor_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    fit_file_manifest = file.path(output_dir, "fit_file_manifest.csv"),
    deletion_file_manifest = file.path(
      output_dir, "deletion_file_manifest.csv"
    ),
    compact_archive_manifest = file.path(
      output_dir, "compact_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(plan$fit_file_manifest, paths[["fit_file_manifest"]])
  readr::write_csv(
    plan$deletion_file_manifest,
    paths[["deletion_file_manifest"]]
  )
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_rescor_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "rescor", "pruning"
    )
) {
  required_confirmation <- "DELETE ALL UPDATED RESCOR FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "rescor"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the rescor archive.",
         call. = FALSE)
  }
  fit_manifest <- plan$fit_file_manifest
  selected <- plan$deletion_file_manifest
  if (
    nrow(fit_manifest) != 128L ||
      sum(fit_manifest$accepted_fit) != 104L ||
      sum(!fit_manifest$accepted_fit) != 24L ||
      nrow(selected) != 644L ||
      sum(selected$file_type == "fit_object") != 128L ||
      sum(selected$file_type == "cmdstan_chain_csv") != 516L ||
      sum(selected$accepted_fit) != 524L ||
      sum(!selected$accepted_fit) != 120L
  ) {
    stop("The rescor pruning set is inconsistent with its audit.",
         call. = FALSE)
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  pending_receipt_path <- file.path(
    receipt_dir,
    "executed_lean.pending.csv"
  )
  if (file.exists(receipt_path) || file.exists(pending_receipt_path)) {
    stop("An updated rescor lean-pruning receipt already exists.",
         call. = FALSE)
  }
  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  current_paths <- unlist(
    purrr::map(
      file.path(plan$archive_root, c("candidates", "remedies")),
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  current_paths <- current_paths[
    file.exists(current_paths) & !file.info(current_paths)$isdir
  ]
  current_paths <- normalizePath(
    current_paths,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  permitted_subtree <- grepl(
    "/(candidates|remedies)/",
    selected_paths
  )
  if (
    !setequal(selected_paths, current_paths) ||
      any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted_subtree) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the audited rescor fit inventory.",
         call. = FALSE)
  }
  compact_paths <- plan$compact_archive_manifest$path
  message("Revalidating 300 compact residual-correlation artifacts...")
  if (any(!file.exists(compact_paths))) {
    stop("A retained compact rescor artifact failed revalidation.",
         call. = FALSE)
  }
  compact_hashes <- purrr::map_chr(
    compact_paths,
    ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
  )
  if (any(compact_hashes != plan$compact_archive_manifest$sha256)) {
    stop("A retained compact rescor artifact failed revalidation.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  message(
    "Hashing 644 residual-correlation sampling files before deletion..."
  )
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, pending_receipt_path)

  message("Deleting the 644 audited residual-correlation sampling files...")
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop(
      "At least one rescor archive could not be deleted; ",
      "the pending receipt records the attempted inventory.",
      call. = FALSE
    )
  }
  final_compact_hashes <- purrr::map_chr(
    compact_paths,
    ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
  )
  if (
    any(!file.exists(compact_paths)) ||
      any(final_compact_hashes != plan$compact_archive_manifest$sha256) ||
      !file.exists(plan$verification_manifest)
  ) {
    stop("A retained rescor artifact is unexpectedly missing.",
         call. = FALSE)
  }

  readr::write_csv(receipt, receipt_path)
  if (!file.remove(pending_receipt_path)) {
    warning("The completed rescor pending receipt could not be removed.")
  }
  message("Rescor lean pruning completed; compact artifacts remain.")
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

write_update_2026_eiv_diagnostic_audit <- function(
    audit,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "eiv",
      "diagnostics", "diagnostic_summary.csv"
    )
) {
  if (
    !identical(audit$family, "eiv") ||
      nrow(audit$diagnostic_summary) != 104L
  ) {
    stop("Cannot write an incomplete updated EIV diagnostic audit.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(audit$diagnostic_summary, output_file)
  output_file
}

# Combined-data strict within-study sensitivity refit ----

make_update_2026_within_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated within-study base specification registry is incomplete.",
         call. = FALSE)
  }
  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  remedy_keys <- unname(purrr::map_chr(primary_remedy_specifications, key))
  if (
    anyDuplicated(primary_keys) || anyDuplicated(remedy_keys) ||
      !identical(remedy_keys, "function::10") ||
      !all(remedy_keys %in% primary_keys)
  ) {
    stop("The updated within-study primary remedy replacement is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  historical_remedies <- c("hypertrophy::9", "strength::9")
  specifications <- purrr::map(specifications, function(specification) {
    specification_key <- key(specification)
    historical_remedy <- specification_key %in% historical_remedies
    if (historical_remedy) {
      specification$adapt_delta <- max(
        as.numeric(specification$adapt_delta),
        0.999
      )
      specification$max_treedepth <- max(
        as.integer(specification$max_treedepth),
        15L
      )
    }
    specification$update_2026_sensitivity_base <- TRUE
    specification$update_2026_within_historical_remedy <- historical_remedy
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  historical_keys <- final_keys[purrr::map_lgl(
    specifications,
    ~ isTRUE(.x$update_2026_within_historical_remedy)
  )]
  function_ten <- specifications[[match("function::10", final_keys)]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !setequal(historical_keys, historical_remedies) ||
      !identical(as.numeric(function_ten$adapt_delta), 0.9995)
  ) {
    stop("The updated within-study specification registry failed validation.",
         call. = FALSE)
  }
  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

update_2026_within_local_config <- function(specification, config) {
  local_config <- config
  historical_remedy <- isTRUE(
    specification$update_2026_within_historical_remedy
  )
  if (historical_remedy) {
    local_config$iter_warmup <- 2000L
    local_config$iter_sampling <- 4000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      15L
    )
  }
  outcome_offset <- match(specification$outcome, config$outcomes)
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the updated within-study configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed +
      (if (historical_remedy) 1000L else 0L) +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

validate_update_2026_within_inputs <- function(
    specifications,
    contrast_data_by_outcome,
    priors_by_outcome,
    config
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  checks <- purrr::map_dfr(specifications, function(specification) {
    local_config <- update_2026_within_local_config(specification, config)
    input <- prepare_within_full_input(
      specification = specification,
      contrast_data_by_outcome = contrast_data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      config = local_config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      contrasts = nrow(input$data),
      studies = dplyr::n_distinct(input$data$study_no),
      matched_arm_pairs = dplyr::n_distinct(input$data$arm_no),
      basis_dimension = if (is.null(input$basis_spec)) {
        0L
      } else {
        as.integer(input$basis_spec$basis_dimension)
      },
      historical_remedy_schedule = isTRUE(
        specification$update_2026_within_historical_remedy
      ),
      adapt_delta = max(
        input$specification$adapt_delta,
        local_config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        local_config$max_treedepth_floor
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = outcomes),
      model_index
    )

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  historical <- checks |>
    dplyr::filter(.data$historical_remedy_schedule)
  standard <- checks |>
    dplyr::filter(!.data$historical_remedy_schedule)
  counts <- table(factor(checks$outcome, levels = outcomes))
  contrast_counts <- checks |>
    dplyr::distinct(outcome, contrasts, studies, matched_arm_pairs)
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      nrow(contrast_counts) != 4L ||
      any(contrast_counts$contrasts < 1L) ||
      !setequal(
        paste(historical$outcome, historical$model_index, sep = "::"),
        c("hypertrophy::9", "strength::9")
      ) ||
      any(historical$iter_warmup != 2000L) ||
      any(historical$iter_sampling != 4000L) ||
      any(historical$adapt_delta < 0.999) ||
      any(historical$max_treedepth < 15L) ||
      any(standard$iter_warmup != 1500L) ||
      any(standard$iter_sampling != 2000L) ||
      config$chains != 4L || config$parallel_chains != 4L
  ) {
    stop("The updated within-study input preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_within_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "within",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated within-study preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_within_candidate_to_file <- function(
    specification,
    contrast_data_by_outcome,
    priors_by_outcome,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "within",
      "candidates"
    )
) {
  local_config <- update_2026_within_local_config(specification, config)
  input <- prepare_within_full_input(
    specification = specification,
    contrast_data_by_outcome = contrast_data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = local_config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "within"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated within-study fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated within-study fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

compile_update_2026_within_diagnostic_audit <- function(diagnostics) {
  split_diagnostics <- split(
    diagnostics,
    purrr::map_chr(diagnostics, ~ .x$summary$outcome[[1L]])
  )
  required <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !setequal(names(split_diagnostics), required) ||
      any(vapply(split_diagnostics[required], length, integer(1L)) != 26L)
  ) {
    stop("Updated within-study diagnostics are incomplete.", call. = FALSE)
  }
  compile_primary_model_diagnostic_audit(
    split_diagnostics$whole_body,
    split_diagnostics$hypertrophy,
    split_diagnostics$strength,
    split_diagnostics$`function`
  )
}

make_update_2026_within_remedy_config <- function(
    initial_config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    seed_offset = 1000L
) {
  config <- initial_config
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(initial_config$seed + seed_offset)
  if (
    config$iter_warmup < initial_config$iter_warmup ||
      config$iter_sampling <= initial_config$iter_sampling ||
      config$adapt_delta_floor != initial_config$adapt_delta_floor ||
      config$max_treedepth_floor != initial_config$max_treedepth_floor
  ) {
    stop("The updated within-study remedy configuration is invalid.",
         call. = FALSE)
  }
  config
}

select_update_2026_within_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  expected_keys <- c(
    "hypertrophy::10",
    "strength::12",
    "strength::23",
    "strength::25"
  )
  failed <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(
      outcome,
      model_index = as.integer(model_index),
      key = paste(outcome, model_index, sep = "::"),
      divergent_transitions,
      maximum_treedepth_transitions
    )
  hypertrophy_ten <- failed |>
    dplyr::filter(.data$key == "hypertrophy::10")
  if (
    nrow(failed) != 4L || anyDuplicated(failed$key) ||
      !setequal(failed$key, expected_keys) ||
      nrow(hypertrophy_ten) != 1L ||
      hypertrophy_ten$divergent_transitions[[1L]] != 1L ||
      sum(failed$divergent_transitions) != 1L ||
      any(failed$maximum_treedepth_transitions != 0L)
  ) {
    stop("The updated within-study remedy set no longer matches its audit.",
         call. = FALSE)
  }

  selected <- purrr::keep(specifications, function(specification) {
    paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    ) %in% failed$key
  }) |>
    purrr::map(function(specification) {
      specification_key <- paste(
        specification$outcome,
        specification$model_index,
        sep = "::"
      )
      if (specification_key == "hypertrophy::10") {
        specification$adapt_delta <- max(
          as.numeric(specification$adapt_delta),
          0.999
        )
        specification$update_2026_within_remedy <-
          "longer_sampling_and_adapt_delta_0.999"
      } else {
        specification$update_2026_within_remedy <- "longer_sampling"
      }
      specification$update_2026_within_remedy_stage <- "remedy_1"
      specification
    })
  selected_keys <- unname(purrr::map_chr(
    selected,
    ~ paste(.x$outcome, .x$model_index, sep = "::")
  ))
  if (
    length(selected) != 4L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The updated within-study remedy specifications are incomplete.",
         call. = FALSE)
  }
  names(selected) <- paste0(
    purrr::map_chr(selected, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(selected, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  selected
}

validate_update_2026_within_remedy_inputs <- function(
    specifications,
    contrast_data_by_outcome,
    priors_by_outcome,
    config
) {
  checks <- purrr::map_dfr(specifications, function(specification) {
    local_config <- update_2026_within_local_config(specification, config)
    input <- prepare_within_full_input(
      specification = specification,
      contrast_data_by_outcome = contrast_data_by_outcome,
      priors_by_outcome = priors_by_outcome,
      config = local_config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      remedy = specification$update_2026_within_remedy,
      adapt_delta = max(
        input$specification$adapt_delta,
        local_config$adapt_delta_floor
      ),
      max_treedepth = max(
        input$specification$max_treedepth,
        local_config$max_treedepth_floor
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      valid = TRUE
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  expected_keys <- c(
    "hypertrophy::10",
    "strength::12",
    "strength::23",
    "strength::25"
  )
  hypertrophy_ten <- checks |>
    dplyr::filter(.data$outcome == "hypertrophy", .data$model_index == 10L)
  other_cases <- checks |>
    dplyr::filter(!(.data$outcome == "hypertrophy" &
      .data$model_index == 10L))
  if (
    nrow(checks) != 4L || anyDuplicated(keys) ||
      !setequal(keys, expected_keys) || !all(checks$valid) ||
      nrow(hypertrophy_ten) != 1L ||
      hypertrophy_ten$remedy[[1L]] !=
        "longer_sampling_and_adapt_delta_0.999" ||
      hypertrophy_ten$adapt_delta[[1L]] != 0.999 ||
      hypertrophy_ten$max_treedepth[[1L]] != 12L ||
      any(other_cases$remedy != "longer_sampling") ||
      any(checks$iter_warmup != 2000L) ||
      any(checks$iter_sampling != 4000L)
  ) {
    stop("The updated within-study remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_within_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "within", "remedy_1",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 4L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated within-study remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

compile_update_2026_within_remedy_diagnostic_audit <- function(
    initial_diagnostics,
    remedy_diagnostics
) {
  diagnostic_key <- function(x) {
    paste(x$summary$outcome[[1L]], x$summary$model_index[[1L]], sep = "::")
  }
  if (length(initial_diagnostics) != 104L || length(remedy_diagnostics) != 4L) {
    stop("Updated within-study remedy diagnostics are incomplete.",
         call. = FALSE)
  }
  initial_keys <- unname(purrr::map_chr(initial_diagnostics, diagnostic_key))
  remedy_keys <- unname(purrr::map_chr(remedy_diagnostics, diagnostic_key))
  expected_keys <- initial_diagnostics |>
    purrr::keep(~ !.x$summary$diagnostic_pass[[1L]]) |>
    purrr::map_chr(diagnostic_key)
  if (
    anyDuplicated(initial_keys) || anyDuplicated(remedy_keys) ||
      !setequal(remedy_keys, expected_keys) ||
      !all(remedy_keys %in% initial_keys)
  ) {
    stop("Updated within-study diagnostic replacement is invalid.",
         call. = FALSE)
  }
  final_diagnostics <- initial_diagnostics
  final_diagnostics[match(remedy_keys, initial_keys)] <- remedy_diagnostics
  compile_update_2026_within_diagnostic_audit(final_diagnostics)
}

make_update_2026_within_accepted_fit_manifest <- function(
    initial_fit_paths,
    remedy_fit_paths,
    final_diagnostic_audit
) {
  if (
    nrow(final_diagnostic_audit$diagnostic_summary) != 104L ||
      any(!final_diagnostic_audit$diagnostic_summary$diagnostic_pass)
  ) {
    stop("Updated within-study fits cannot be accepted before diagnostics pass.",
         call. = FALSE)
  }
  path_table <- function(paths, stage, expected) {
    paths <- unname(unlist(paths, use.names = FALSE))
    table <- tibble::tibble(
      outcome = basename(dirname(paths)),
      model_index = as.integer(stringr::str_match(
        basename(paths),
        "model_([0-9]+)\\.qs$"
      )[, 2]),
      stage = stage,
      fit_path = normalizePath(paths, winslash = "/", mustWork = TRUE)
    )
    if (
      nrow(table) != expected || anyNA(table$model_index) ||
        anyDuplicated(table[c("outcome", "model_index")]) ||
        any(file.info(table$fit_path)$size <= 0)
    ) {
      stop("An updated within-study fit stage is incomplete.",
           call. = FALSE)
    }
    table
  }

  initial <- path_table(initial_fit_paths, "initial", 104L)
  remedy <- path_table(remedy_fit_paths, "remedy_1", 4L)
  keys <- c("outcome", "model_index")
  remedy_keys <- paste(remedy$outcome, remedy$model_index, sep = "::")
  expected_remedy_keys <- c(
    "hypertrophy::10",
    "strength::12",
    "strength::23",
    "strength::25"
  )
  if (!setequal(remedy_keys, expected_remedy_keys)) {
    stop("Updated within-study accepted replacements are invalid.",
         call. = FALSE)
  }

  accepted <- initial |>
    dplyr::anti_join(remedy[keys], by = keys) |>
    dplyr::bind_rows(remedy) |>
    dplyr::left_join(
      final_diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, model_name = model),
      by = keys,
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      family = "within",
      bytes = file.info(fit_path)$size,
      sha256 = purrr::map_chr(
        fit_path,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      )
    ) |>
    dplyr::select(
      family, outcome, model_index, model_name, stage,
      fit_path, bytes, sha256
    ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )

  stage_counts <- table(factor(
    accepted$stage,
    levels = c("initial", "remedy_1")
  ))
  if (
    nrow(accepted) != 104L ||
      anyDuplicated(accepted[c("outcome", "model_index")]) ||
      anyNA(accepted$model_name) || any(!file.exists(accepted$fit_path)) ||
      !identical(as.integer(stage_counts), c(100L, 4L)) ||
      any(nchar(accepted$sha256) != 64L)
  ) {
    stop("The updated within-study accepted-fit manifest is invalid.",
         call. = FALSE)
  }
  accepted
}

write_update_2026_within_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "within",
      "accepted_fit_manifest.csv"
    )
) {
  if (
    nrow(manifest) != 104L ||
      anyDuplicated(manifest[c("outcome", "model_index")]) ||
      any(!file.exists(manifest$fit_path)) ||
      any(nchar(manifest$sha256) != 64L)
  ) {
    stop("Cannot write an invalid updated within-study fit manifest.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

compute_update_2026_within_integrated_loo_to_file <- function(
    fit_record,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "within", "loo",
      "integrated"
    )
) {
  result <- compute_mundlak_full_file_integrated_loo(
    fit_record = fit_record,
    loo_draws = loo_draws,
    k_threshold = k_threshold,
    cores = cores,
    family = "within",
    slope_description = "within_contrast_basis_slopes"
  )
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(unname(as.character(result$family)), "within") ||
      !identical(result$remedy, "integrated_observation_level_psis")
  ) {
    stop("Cannot archive an invalid updated within-study LOO result.",
         call. = FALSE)
  }
  output_dir <- file.path(output_root, result$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf(
      "model_%02d_integrated_loo.qs",
      as.integer(result$model_index)
    )
  )
  write_atomic_qs(result, output_file)
}

write_update_2026_within_integrated_loo_audit <- function(
    audit,
    manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "loo_psis"
    )
) {
  if (
    !isTRUE(audit$integrated_psis_complete) || audit$models != 104L ||
      nrow(manifest) != 104L ||
      any(!file.exists(manifest$integrated_loo_path))
  ) {
    stop("The updated within-study integrated-LOO audit is incomplete.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

validate_update_2026_within_exact_inputs <- function(
    records,
    integrated_audit,
    config
) {
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    !is.list(records) || length(records) != expected || expected != 59L
  ) {
    stop("Expected 59 updated within-study exact-LOO records.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    integrated <- qs::qread(record$integrated_loo_path)
    observation <- as.integer(record$observation)
    if (
      !inherits(integrated$loo, "loo") ||
        !identical(unname(as.character(integrated$family)), "within") ||
        !identical(integrated$outcome, record$outcome) ||
        as.integer(integrated$model_index) != record$model_index ||
        observation < 1L ||
        observation > nrow(integrated$observation_manifest)
    ) {
      stop("An updated within-study exact record is invalid.",
           call. = FALSE)
    }
    held_out <- integrated$observation_manifest[observation, , drop = FALSE]
    observed_k <- as.numeric(
      loo::pareto_k_values(integrated$loo)[[observation]]
    )
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      observation = observation,
      study_no = as.character(held_out$study_no[[1L]]),
      arm_no = as.integer(held_out$arm_no[[1L]]),
      outcome_no = as.integer(held_out$outcome_no[[1L]]),
      pareto_k = observed_k,
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      original_fit_exists = file.exists(record$original_fit_path),
      integrated_loo_exists = file.exists(record$integrated_loo_path),
      valid = is.finite(observed_k) && observed_k > 0.7
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  if (
    nrow(checks) != expected || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(expected)) ||
      !identical(as.integer(outcome_counts), c(6L, 26L, 1L, 26L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$integrated_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L
  ) {
    stop("The updated within-study exact-LOO preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_within_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "within", "loo_exact",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 59L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated within-study exact preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_within_exact_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "within", "loo",
    "exact", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

fit_and_score_update_2026_within_exact_to_file <- function(
    record,
    config,
    output_file = update_2026_within_exact_file_path(record)
) {
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !identical(unname(as.character(candidate$family)), "within") ||
      !identical(candidate$outcome, record$outcome) ||
      as.integer(candidate$model_index) != record$model_index ||
      !inherits(integrated$loo, "loo") ||
      !identical(unname(as.character(integrated$family)), "within") ||
      !identical(integrated$outcome, record$outcome) ||
      as.integer(integrated$model_index) != record$model_index ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      loo::pareto_k_values(integrated$loo)[[observation]] <= 0.7
  ) {
    stop("An updated within-study exact-LOO record failed validation.",
         call. = FALSE)
  }

  original_data <- candidate$model$data
  training_data <- original_data[-observation, , drop = FALSE]
  original_model <- candidate$model
  sample_prior <- attr(original_model$prior, "sample_prior")
  original_model$prior <- normalize_primary_brms_priors(
    original_model$prior
  )
  attr(original_model$prior, "sample_prior") <- sample_prior
  specification <- candidate$specification
  specification$adapt_delta <- max(
    specification$adapt_delta,
    config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth,
    config$max_treedepth_floor
  )

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
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  refit <- list(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    original_data = original_data,
    integrated_result = integrated,
    specification = specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  correction <- compute_sensitivity_mundlak_exact_loo_case(refit)
  if (
    as.integer(correction$case_id) != record$case_id ||
      !inherits(correction$loo, "loo") ||
      nrow(correction$exact_refit_summary) != 1L
  ) {
    stop("The updated within-study exact correction is invalid.",
         call. = FALSE)
  }
  write_atomic_qs(
    correction,
    output_file
  )
}

audit_update_2026_within_exact_files <- function(paths, integrated_audit) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    length(paths) != expected || expected != 59L || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0)
  ) {
    stop("The updated within-study exact-correction files are incomplete.",
         call. = FALSE)
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  manifest <- purrr::map2_dfr(paths, corrections, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = result$outcome,
      model_index = as.integer(result$model_index),
      model_name = result$model_name,
      observation =
        as.integer(result$exact_refit_summary$observation[[1L]]),
      correction_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys)
  ) {
    stop("The updated within-study exact audit is invalid.",
         call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = sum(diagnostic_summary$diagnostic_pass),
    exact_refits_requiring_remediation =
      sum(!diagnostic_summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass =
      all(diagnostic_summary$diagnostic_pass)
  )
}

write_update_2026_within_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "loo_exact"
    )
) {
  if (
    audit$exact_refits != 59L ||
      nrow(audit$diagnostic_summary) != 59L ||
      nrow(audit$correction_manifest) != 59L
  ) {
    stop("Cannot write an incomplete updated within-study exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

make_update_2026_within_exact_remedy_config <- function(
    exact_config,
    adapt_delta = 0.999,
    seed_offset = 1000L
) {
  config <- exact_config
  config$adapt_delta_floor <- max(
    as.numeric(exact_config$adapt_delta_floor),
    as.numeric(adapt_delta)
  )
  config$seed <- as.integer(exact_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.999 ||
      config$max_treedepth_floor != exact_config$max_treedepth_floor
  ) {
    stop("The updated within-study exact-remedy configuration is invalid.",
         call. = FALSE)
  }
  config
}

select_update_2026_within_exact_remedy_record <- function(
    records,
    exact_audit
) {
  failed <- exact_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass)
  if (
    nrow(failed) != 1L || failed$outcome[[1L]] != "whole_body" ||
      failed$model_index[[1L]] != 21L ||
      failed$observation[[1L]] != 2L ||
      failed$divergent_transitions[[1L]] != 1L ||
      failed$maximum_treedepth_transitions[[1L]] != 0L ||
      !isTRUE(failed$convergence_pass[[1L]])
  ) {
    stop("The within-study exact remedy no longer matches its audit.",
         call. = FALSE)
  }
  selected <- purrr::keep(records, function(record) {
    identical(record$outcome, "whole_body") &&
      as.integer(record$model_index) == 21L &&
      as.integer(record$observation) == 2L
  })
  if (
    length(selected) != 1L ||
      as.integer(selected[[1L]]$case_id) != 6L
  ) {
    stop("Could not uniquely select within-study exact case 6.",
         call. = FALSE)
  }
  names(selected) <- "case_006_whole_body_model_21_observation_0002"
  selected
}

validate_update_2026_within_exact_remedy_input <- function(record, config) {
  if (!is.list(record) || length(record) != 1L) {
    stop("Expected one updated within-study exact-remedy record.",
         call. = FALSE)
  }
  record <- record[[1L]]
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  observed_k <- as.numeric(
    loo::pareto_k_values(integrated$loo)[[observation]]
  )
  check <- tibble::tibble(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    remedy = "adapt_delta_0.999_after_exact_refit",
    adapt_delta = as.numeric(config$adapt_delta_floor),
    max_treedepth = as.integer(config$max_treedepth_floor),
    iter_warmup = as.integer(config$iter_warmup),
    iter_sampling = as.integer(config$iter_sampling),
    pareto_k = observed_k,
    original_fit_exists = file.exists(record$original_fit_path),
    integrated_loo_exists = file.exists(record$integrated_loo_path),
    valid = TRUE
  )
  if (
    check$case_id[[1L]] != 6L || check$outcome[[1L]] != "whole_body" ||
      check$model_index[[1L]] != 21L || check$observation[[1L]] != 2L ||
      check$adapt_delta[[1L]] != 0.999 ||
      check$max_treedepth[[1L]] != 12L ||
      check$iter_warmup[[1L]] != 2000L ||
      check$iter_sampling[[1L]] != 6000L ||
      !is.finite(observed_k) || observed_k <= 0.7 ||
      !check$original_fit_exists[[1L]] ||
      !check$integrated_loo_exists[[1L]]
  ) {
    stop("The updated within-study exact-remedy preflight is invalid.",
         call. = FALSE)
  }
  check
}

write_update_2026_within_exact_remedy_preflight <- function(
    check,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "within", "loo_exact",
      "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(check) != 1L || !isTRUE(check$valid[[1L]])) {
    stop("Cannot write an invalid within-study exact-remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(check, output_file)
  output_file
}

update_2026_within_exact_remedy_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "within", "loo",
    "exact_remedy", record$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

audit_update_2026_within_exact_remedy_file <- function(path, record) {
  path <- unname(unlist(path, use.names = FALSE))
  if (
    is.list(record) && is.null(record$case_id) && length(record) == 1L
  ) {
    record <- record[[1L]]
  }
  if (
    length(path) != 1L || !file.exists(path) || file.info(path)$size <= 0
  ) {
    stop("The updated within-study exact-remedy file is unavailable.",
         call. = FALSE)
  }
  correction <- qs::qread(path)
  summary <- correction$exact_refit_summary
  if (
    as.integer(correction$case_id) != 6L || nrow(summary) != 1L ||
      summary$outcome[[1L]] != "whole_body" ||
      summary$model_index[[1L]] != 21L ||
      summary$observation[[1L]] != 2L ||
      as.integer(record$case_id) != 6L
  ) {
    stop("The updated within-study exact-remedy result is invalid.",
         call. = FALSE)
  }
  manifest <- tibble::tibble(
    case_id = 6L,
    outcome = summary$outcome[[1L]],
    model_index = as.integer(summary$model_index[[1L]]),
    model_name = correction$model_name,
    observation = as.integer(summary$observation[[1L]]),
    correction_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = digest::digest(
      file = path,
      algo = "sha256",
      serialize = FALSE
    )
  )
  list(
    diagnostic_summary = summary,
    correction_manifest = manifest,
    exact_refits = 1L,
    exact_refits_passing = sum(summary$diagnostic_pass),
    exact_refits_requiring_remediation = sum(!summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

write_update_2026_within_exact_remedy_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "loo_exact",
      "remedy_1"
    )
) {
  if (
    audit$exact_refits != 1L ||
      nrow(audit$diagnostic_summary) != 1L ||
      nrow(audit$correction_manifest) != 1L
  ) {
    stop("Cannot write an invalid within-study exact-remedy audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

make_update_2026_within_final_exact_audit <- function(
    initial_audit,
    remedy_audit
) {
  if (
    initial_audit$exact_refits != 59L ||
      initial_audit$exact_refits_requiring_remediation != 1L ||
      remedy_audit$exact_refits != 1L ||
      !isTRUE(remedy_audit$all_exact_refit_diagnostics_pass)
  ) {
    stop("The within-study exact audits are not ready to combine.",
         call. = FALSE)
  }
  key <- function(data) {
    paste(data$outcome, data$model_index, data$observation, sep = "::")
  }
  remedy_key <- key(remedy_audit$diagnostic_summary)
  expected_key <- "whole_body::21::2"
  if (!identical(remedy_key, expected_key)) {
    stop("The within-study exact remedy has the wrong identity.",
         call. = FALSE)
  }
  diagnostic_summary <- initial_audit$diagnostic_summary |>
    dplyr::filter(key(dplyr::pick(dplyr::everything())) != remedy_key) |>
    dplyr::bind_rows(remedy_audit$diagnostic_summary) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  correction_manifest <- initial_audit$correction_manifest |>
    dplyr::filter(.data$case_id != 6L) |>
    dplyr::bind_rows(remedy_audit$correction_manifest) |>
    dplyr::arrange(case_id)
  if (
    nrow(diagnostic_summary) != 59L ||
      anyDuplicated(key(diagnostic_summary)) ||
      any(!diagnostic_summary$diagnostic_pass) ||
      nrow(correction_manifest) != 59L ||
      anyDuplicated(correction_manifest$case_id) ||
      !identical(correction_manifest$case_id, seq_len(59L)) ||
      any(!file.exists(correction_manifest$correction_path)) ||
      correction_manifest$correction_path[
        correction_manifest$case_id == 6L
      ] != remedy_audit$correction_manifest$correction_path[[1L]]
  ) {
    stop("The final within-study exact audit is invalid.", call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = correction_manifest,
    exact_refits = 59L,
    exact_refits_passing = 59L,
    exact_refits_requiring_remediation = 0L,
    all_exact_refit_diagnostics_pass = TRUE,
    remedied_case_id = 6L
  )
}

write_update_2026_within_final_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "loo_exact",
      "final"
    )
) {
  if (
    audit$exact_refits != 59L ||
      !isTRUE(audit$all_exact_refit_diagnostics_pass) ||
      nrow(audit$diagnostic_summary) != 59L ||
      nrow(audit$correction_manifest) != 59L
  ) {
    stop("Cannot write an incomplete final within-study exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

read_update_2026_within_exact_corrections <- function(manifest) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 59L ||
      any(!file.exists(manifest$correction_path))
  ) {
    stop("The updated within-study exact-correction manifest is invalid.",
         call. = FALSE)
  }
  corrections <- purrr::map(manifest$correction_path, qs::qread)
  corrections <- corrections[order(purrr::map_int(corrections, "case_id"))]
  if (
    !identical(
      unname(purrr::map_int(corrections, "case_id")),
      seq_len(59L)
    ) ||
      any(!purrr::map_lgl(corrections, ~ inherits(.x$loo, "loo")))
  ) {
    stop("The updated within-study exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections
}

make_update_2026_within_final_loo_audit <- function(
    integrated_manifest,
    exact_corrections,
    exact_audit
) {
  if (!isTRUE(exact_audit$all_exact_refit_diagnostics_pass)) {
    stop("Final within-study LOO requires passing exact diagnostics.",
         call. = FALSE)
  }
  final <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest = integrated_manifest,
    exact_results = exact_corrections,
    family = "within"
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- integrated_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "integrated_observation_level_psis_plus_exact",
        "integrated_observation_level_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 59L
  ) {
    stop("The final updated within-study LOO is not stacking-ready.",
         call. = FALSE)
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 59L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_within_final_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 59L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop("Cannot write an incomplete final within-study LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_within_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "within",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L ||
      anyDuplicated(paste(weights$outcome, weights$model_index)) ||
      length(outcome_sums) != 4L ||
      any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop("The updated within-study stacking weights are invalid.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_within_reporting_preflight <- function(records) {
  if (!is.list(records) || length(records) < 4L) {
    stop("The updated within-study prediction allocation is incomplete.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    tibble::tibble(
      outcome = record$outcome,
      model_index = as.integer(record$model_index),
      model_name = record$model_name,
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      grid_step = as.numeric(record$grid_step),
      fit_path = record$fit_path,
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path)
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop("The updated within-study reporting preflight failed.",
         call. = FALSE)
  }
  checks
}

write_update_2026_within_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "within", "reporting",
      "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(checks$allocated_draws, checks$outcome, sum)
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      !all(checks$valid)
  ) {
    stop("Cannot write an invalid updated within-study reporting preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_within_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "within",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_within_prediction_draw_file <- function(
    record,
    data_by_outcome
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !identical(unname(as.character(record$family)), "within") ||
      !identical(unname(as.character(record$prediction_type)), "within")
  ) {
    stop("An updated within-study prediction record is invalid.",
         call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "within") ||
      !identical(candidate$outcome, unname(as.character(record$outcome))) ||
      as.integer(candidate$model_index) != as.integer(record$model_index) ||
      !identical(candidate$model_name, unname(as.character(record$model_name)))
  ) {
    stop("An updated within-study prediction fit does not match its record.",
         call. = FALSE)
  }
  prediction_data <- data_by_outcome[[record$outcome]]
  if (is.null(prediction_data) || !is.data.frame(prediction_data)) {
    stop("The updated within-study prediction data are unavailable.",
         call. = FALSE)
  }
  grid <- make_sensitivity_prediction_grid(
    fitted_candidate = candidate,
    data = prediction_data,
    prediction_type = "within",
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
    stop("Updated within-study prediction draws are incomplete.",
         call. = FALSE)
  }
  reduction <- list(
    family = "within",
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
      prediction_type = "within",
      reference_intake = 0.8
    )
  )
  write_atomic_qs(
    reduction,
    update_2026_within_reporting_path(
      outcome = record$outcome,
      stage = "prediction_draws",
      model_index = record$model_index
    )
  )
}

write_update_2026_within_model_average_file <- function(record) {
  reductions <- purrr::map(
    record$prediction_manifest$prediction_path,
    read_sensitivity_full_prediction_draw_file,
    family = "within"
  )
  key_reference <- reductions[[1L]]$grid
  if (!all(purrr::map_lgl(
    reductions, ~ identical(.x$grid, key_reference)
  ))) {
    stop("Updated within-study prediction grids do not match within outcome.",
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
          file = path, algo = "sha256", serialize = FALSE
        ),
        within_model_draw = seq_len(reduction$draw_count)
      )
    }
  )
  if (
    nrow(draws) != record$n_draws ||
      nrow(draw_manifest) != record$n_draws || anyNA(draws)
  ) {
    stop("The updated within-study posterior mixture is incomplete.",
         call. = FALSE)
  }
  set.seed(record$averaging_seed)
  draw_order <- sample(seq_len(record$n_draws))
  draws <- draws[draw_order, , drop = FALSE]
  draw_manifest <- draw_manifest[draw_order, , drop = FALSE] |>
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
  write_atomic_qs(
    output,
    update_2026_within_reporting_path(
      outcome = record$outcome,
      stage = "model_average"
    )
  )
}

write_update_2026_within_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path, family = "within"
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (is.null(data) || !is.data.frame(data)) {
    stop("The updated within-study reporting data are unavailable.",
         call. = FALSE)
  }
  absolute <- summarise_primary_dose_response(model_average)
  supported <- build_sensitivity_within_reporting_pilot(
    absolute_summary = absolute,
    model_average = model_average,
    data = data,
    reference_intake = reference_intake
  )
  heatmap <- compute_primary_heatmap_contrasts(
    model_average = model_average,
    data = dplyr::mutate(
      data,
      protein_intake_raw = protein_intake_raw_treatment
    ),
    protein_min = reference_intake,
    sesoi = sesoi
  )
  report <- list(
    family = "within",
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
  write_atomic_qs(
    report,
    update_2026_within_reporting_path(
      outcome = outcome,
      stage = "reporting"
    )
  )
}

write_update_2026_within_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "reporting"
    )
) {
  if (
    !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L
  ) {
    stop("The updated within-study reporting audit is incomplete.",
         call. = FALSE)
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "within"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir, "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir, "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir, "common_exposure_support.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_within_verification_checkpoint <- function(
    data_by_outcome,
    source_manifest,
    review_decisions,
    preparation_audit,
    specifications,
    accepted_fit_manifest,
    diagnostic_audit,
    integrated_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "within"
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(data) {
    data |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_audit$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_within_name <- paste0(
    name_alignment$specification_name,
    " [within]"
  )
  weight_sums <- stacking_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      .by = outcome
    )

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "all four within-study contrast datasets are present and complete",
      "the combined-data preparation audit remains valid",
      "all five within-study registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist and have checksums",
      "all final model diagnostics pass",
      "integrated LOO contains 104 reductions and 59 audited high-k cases",
      "all 59 final exact refits pass and compact corrections exist",
      "final corrected LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four model-averaged reporting audits pass"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      identical(names(data_by_outcome), outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L) &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) && all(analysis_data$vi > 0),
      all(preparation_audit$checks$passed),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_within_name) &&
        all(name_alignment$diagnostic_name == expected_within_name) &&
        all(name_alignment$loo_name == expected_within_name) &&
        all(name_alignment$weight_name == expected_within_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      nrow(diagnostic_audit$diagnostic_summary) == 104L &&
        all(diagnostic_audit$diagnostic_summary$diagnostic_pass),
      integrated_loo_audit$models == 104L &&
        integrated_loo_audit$exact_refits_required == 59L &&
        nrow(integrated_loo_audit$high_k_cases) == 59L,
      exact_audit$exact_refits == 59L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        all(file.exists(exact_audit$correction_manifest$correction_path)),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(file.exists(prediction_manifest$prediction_path)),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated within-study verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    review_decisions = review_decisions,
    preparation_checks = preparation_audit$checks,
    preparation_outcome_counts = preparation_audit$outcome_counts,
    preparation_update_study_counts = preparation_audit$update_study_counts,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_audit$diagnostic_summary,
    sampling_diagnostic_classification =
      diagnostic_audit$definitive_classification,
    sampling_parameter_failures = diagnostic_audit$parameter_failures,
    posterior_predictive_observation_flags =
      diagnostic_audit$ppc_observation_flags,
    posterior_predictive_outcome_summary =
      diagnostic_audit$ppc_outcome_summary,
    initial_integrated_loo_model_summary = integrated_loo_audit$model_summary,
    initial_integrated_loo_outcome_summary =
      integrated_loo_audit$outcome_summary,
    initial_high_pareto_k_cases = integrated_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_correction_manifest = exact_audit$correction_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final combined-data within-study checkpoint checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    review_decisions = "Eligibility decisions for updated-search studies.",
    preparation_checks = "Combined-data preparation invariants.",
    preparation_outcome_counts = "Analysis counts by primary outcome.",
    preparation_update_study_counts = "Updated-study contribution audit.",
    analysis_data = "Combined within-study contrasts supplied to models.",
    model_specifications = "Final 104-model within-study specification registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit path, stage, size, and checksum.",
    sampling_diagnostic_summary = "Final model-level sampling diagnostics.",
    sampling_diagnostic_classification = "Final diagnostic decisions.",
    sampling_parameter_failures = "Parameter-level diagnostic failures.",
    posterior_predictive_observation_flags = "Observation-level PPC flags.",
    posterior_predictive_outcome_summary = "Outcome-level PPC summary.",
    initial_integrated_loo_model_summary = "Initial integrated PSIS summaries.",
    initial_integrated_loo_outcome_summary = "Initial LOO audit by outcome.",
    initial_high_pareto_k_cases = "The 59 cases selected for exact refitting.",
    exact_refit_diagnostic_summary = "Final diagnostics for 59 exact refits.",
    exact_correction_manifest = "Compact exact-correction paths and checksums.",
    final_loo_manifest = "Final exact-corrected LOO lineage.",
    final_loo_model_summary = "Final LOO summary for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    stacking_weights = "Final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    dose_response_predictions = "Within-study model-averaged summaries.",
    dose_response_increments = "Increments from the 0.8 g/kg/day reference.",
    exposure_support_by_group = "Observed within-study exposure support.",
    common_exposure_support = "Common within-study support by outcome.",
    pairwise_contrasts = "Model-averaged pairwise contrasts and probabilities.",
    pairwise_support = "Support limits for pairwise contrasts.",
    local_equivalence = "Local practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    figures_deferred = paste(
      "Cross-family primary/Mundlak/within figures are now ready to build",
      "after this verification checkpoint is frozen."
    )
  )
}

summarise_update_2026_primary_reference_increments <- function(
    model_average_path,
    outcome,
    reference_intake = 0.8
) {
  if (
    length(model_average_path) != 1L || is.na(model_average_path) ||
      !file.exists(model_average_path)
  ) {
    stop("An updated primary model-average file is unavailable.",
         call. = FALSE)
  }
  model_average <- qs::qread(model_average_path)
  if (
    !inherits(model_average, "primary_model_average") ||
      !identical(model_average$outcome, outcome) ||
      !is.matrix(model_average$draws) ||
      nrow(model_average$draws) != 4000L
  ) {
    stop("An updated primary model-average object is invalid.",
         call. = FALSE)
  }
  grid <- model_average$grid
  increment_draws <- model_average$draws
  groups <- sort(unique(grid$resistance_exercise_code))
  for (group in groups) {
    group_rows <- which(grid$resistance_exercise_code == group)
    reference_rows <- group_rows[
      abs(grid$protein_intake[group_rows] - reference_intake) < 1e-10
    ]
    if (length(reference_rows) != 1L) {
      stop("The updated primary increment reference is not unique.",
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
  summarise_primary_dose_response(increment_average) |>
    dplyr::transmute(
      outcome = outcome,
      resistance_exercise_code,
      resistance_exercise,
      protein_intake,
      estimate,
      conf.low,
      conf.high,
      analysis = "Primary"
    )
}

make_update_2026_primary_mundlak_within_figure_data <- function(
    primary_reporting_audit,
    mundlak_reporting_audit,
    within_reporting_audit,
    reference_intake = 0.8
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !isTRUE(primary_reporting_audit$all_reporting_checks_pass) ||
      !isTRUE(mundlak_reporting_audit$all_reporting_checks_pass) ||
      !isTRUE(within_reporting_audit$all_reporting_checks_pass)
  ) {
    stop("Cross-family figures require three passing reporting audits.",
         call. = FALSE)
  }
  primary_paths <- primary_reporting_audit$manifest |>
    dplyr::arrange(match(outcome, outcomes)) |>
    dplyr::select(outcome, model_average_path)
  if (
    nrow(primary_paths) != 4L ||
      !identical(primary_paths$outcome, outcomes) ||
      any(!file.exists(primary_paths$model_average_path))
  ) {
    stop("The updated primary model-average registry is incomplete.",
         call. = FALSE)
  }
  primary_increments <- purrr::map2_dfr(
    primary_paths$model_average_path,
    primary_paths$outcome,
    summarise_update_2026_primary_reference_increments,
    reference_intake = reference_intake
  )
  primary_absolute <- primary_reporting_audit$dose_response_predictions |>
    dplyr::transmute(
      outcome,
      resistance_exercise_code,
      resistance_exercise,
      protein_intake,
      estimate,
      conf.low,
      conf.high,
      analysis = "Primary"
    )

  mundlak_reports <- purrr::map(
    mundlak_reporting_audit$summary$reporting_path,
    read_mundlak_full_reporting_file
  )
  within_reports <- purrr::map(
    within_reporting_audit$summary$reporting_path,
    read_sensitivity_full_reporting_file,
    family = "within"
  )
  bind_component <- function(reports, component, label, increment = FALSE) {
    purrr::map_dfr(reports, function(report) {
      table <- report[[component]]
      if (isTRUE(increment)) {
        dplyr::transmute(
          table,
          outcome = report$outcome,
          resistance_exercise_code,
          resistance_exercise,
          protein_intake,
          estimate = increment,
          conf.low = increment_conf_low,
          conf.high = increment_conf_high,
          analysis = label
        )
      } else {
        dplyr::transmute(
          table,
          outcome = report$outcome,
          resistance_exercise_code,
          resistance_exercise,
          protein_intake,
          estimate,
          conf.low,
          conf.high,
          analysis = label
        )
      }
    })
  }
  mundlak_absolute <- bind_component(
    mundlak_reports, "dose_response", "Mundlak"
  )
  mundlak_increments <- bind_component(
    mundlak_reports,
    "dose_response_increments",
    "Mundlak",
    increment = TRUE
  )
  within_increments <- bind_component(
    within_reports,
    "dose_response_increments",
    "Strict within-study",
    increment = TRUE
  )
  within_support <- purrr::map_dfr(within_reports, function(report) {
    dplyr::mutate(
      report$exposure_support_by_group,
      outcome = report$outcome,
      .before = 1L
    )
  }) |>
    dplyr::select(
      outcome,
      resistance_exercise_code,
      observed_minimum,
      observed_maximum
    )

  absolute <- dplyr::bind_rows(primary_absolute, mundlak_absolute) |>
    dplyr::arrange(
      factor(outcome, levels = outcomes),
      resistance_exercise_code,
      factor(analysis, levels = c("Primary", "Mundlak")),
      protein_intake
    )
  unfiltered_increments <- dplyr::bind_rows(
    primary_increments,
    mundlak_increments,
    within_increments
  )
  source_reference_rows <- unfiltered_increments |>
    dplyr::filter(abs(protein_intake - reference_intake) < 1e-10)
  increments <- unfiltered_increments |>
    dplyr::left_join(
      within_support,
      by = c("outcome", "resistance_exercise_code"),
      relationship = "many-to-one"
    ) |>
    dplyr::filter(
      protein_intake >= observed_minimum - 1e-10,
      protein_intake <= observed_maximum + 1e-10
    ) |>
    dplyr::arrange(
      factor(outcome, levels = outcomes),
      resistance_exercise_code,
      factor(
        analysis,
        levels = c("Primary", "Mundlak", "Strict within-study")
      ),
      protein_intake
    )

  finite_intervals <- function(data) {
    nrow(data) > 0L &&
      all(is.finite(data$estimate)) &&
      all(is.finite(data$conf.low)) &&
      all(is.finite(data$conf.high)) &&
      all(data$conf.low <= data$estimate) &&
      all(data$estimate <= data$conf.high)
  }
  increment_cells <- increments |>
    dplyr::distinct(outcome, resistance_exercise_code, analysis)
  checks <- tibble::tibble(
    check = c(
      "absolute curves contain primary and Mundlak estimates for four outcomes",
      "increment curves contain all three analyses in eight outcome-RT panels",
      "all plotted intervals are finite and ordered",
      "within-study support bounds are complete and respected",
      "all source reference increments are exactly zero before support restriction"
    ),
    passed = c(
      setequal(unique(absolute$analysis), c("Primary", "Mundlak")) &&
        setequal(unique(absolute$outcome), outcomes),
      nrow(increment_cells) == 24L &&
        setequal(
          unique(increments$analysis),
          c("Primary", "Mundlak", "Strict within-study")
        ),
      finite_intervals(absolute) && finite_intervals(increments),
      nrow(within_support) == 8L &&
        !anyNA(increments$observed_minimum) &&
        !anyNA(increments$observed_maximum) &&
        all(increments$protein_intake >= increments$observed_minimum - 1e-10) &&
        all(increments$protein_intake <= increments$observed_maximum + 1e-10),
      nrow(source_reference_rows) == 24L &&
        all(abs(source_reference_rows$estimate) < 1e-12) &&
        all(abs(source_reference_rows$conf.low) < 1e-12) &&
        all(abs(source_reference_rows$conf.high) < 1e-12)
    )
  )
  if (any(is.na(checks$passed) | !checks$passed)) {
    stop(
      "Updated cross-family figure data failed validation: ",
      paste(checks$check[is.na(checks$passed) | !checks$passed],
            collapse = "; "),
      call. = FALSE
    )
  }
  list(
    absolute = absolute,
    increments = increments,
    within_support = within_support,
    checks = checks,
    reference_intake = reference_intake
  )
}

write_update_2026_primary_mundlak_within_figures <- function(
    figure_data,
    data_dir = here::here(
      "results", "update_2026", "sensitivity", "comparisons"
    ),
    figure_dir = here::here(
      "results", "update_2026", "figures", "sensitivity_comparisons"
    )
) {
  if (!all(figure_data$checks$passed)) {
    stop("Cannot plot failed updated cross-family data.", call. = FALSE)
  }
  outcome_names <- c(
    whole_body = "Whole-body lean soft tissue",
    hypertrophy = "Direct muscle mass / size",
    strength = "Strength",
    `function` = "Physical function"
  )
  prepare_plot_data <- function(data, analyses) {
    data |>
      dplyr::mutate(
        outcome = factor(
          outcome,
          levels = names(outcome_names),
          labels = unname(outcome_names)
        ),
        resistance_exercise = factor(
          resistance_exercise_code,
          levels = c(1, 0),
          labels = c(
            "Resistance exercise: Yes",
            "Resistance exercise: No"
          )
        ),
        panel = interaction(
          outcome,
          resistance_exercise,
          sep = "  |  ",
          lex.order = TRUE
        ),
        analysis = factor(analysis, levels = analyses)
      )
  }
  absolute_plot_data <- prepare_plot_data(
    figure_data$absolute,
    c("Primary", "Mundlak")
  )
  increment_plot_data <- prepare_plot_data(
    figure_data$increments,
    c("Primary", "Mundlak", "Strict within-study")
  )
  absolute_colours <- c(Primary = "#2B6CB0", Mundlak = "#D97706")
  increment_colours <- c(
    Primary = "#2B6CB0",
    Mundlak = "#D97706",
    `Strict within-study` = "#238B45"
  )
  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(
        colour = "grey30",
        margin = ggplot2::margin(b = 8)
      ),
      plot.caption = ggplot2::element_text(colour = "grey35", hjust = 0),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "grey90", linewidth = 0.3
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey88", linewidth = 0.3
      ),
      legend.position = "top",
      legend.justification = "left",
      legend.key.width = grid::unit(1.3, "cm"),
      plot.margin = ggplot2::margin(12, 16, 12, 12)
    )
  absolute_plot <- ggplot2::ggplot(
    absolute_plot_data,
    ggplot2::aes(
      x = protein_intake,
      y = estimate,
      colour = analysis,
      fill = analysis,
      group = analysis
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, colour = "grey72", linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = figure_data$reference_intake,
      colour = "grey55",
      linewidth = 0.4,
      linetype = "dotted"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.13,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(ggplot2::vars(panel), ncol = 2, scales = "free") +
    ggplot2::scale_colour_manual(values = absolute_colours) +
    ggplot2::scale_fill_manual(values = absolute_colours) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_width(0.4),
      expand = ggplot2::expansion(mult = c(0.015, 0.025))
    ) +
    ggplot2::labs(
      title = "Primary and Mundlak model-averaged dose-response curves",
      subtitle = paste(
        "Updated combined dataset; lines are posterior means and shaded",
        "regions are 95% credible intervals"
      ),
      x = "Protein intake (g/kg/day)",
      y = "Standardised mean change",
      colour = NULL,
      fill = NULL,
      caption = paste(
        "The dotted line marks 0.8 g/kg/day. Primary curves combine",
        "within- and between-study information; Mundlak curves separate",
        "study-level mean intake from within-study intake."
      )
    ) +
    base_theme
  increment_plot <- ggplot2::ggplot(
    increment_plot_data,
    ggplot2::aes(
      x = protein_intake,
      y = estimate,
      colour = analysis,
      fill = analysis,
      group = analysis
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, colour = "grey65", linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = figure_data$reference_intake,
      colour = "grey50",
      linewidth = 0.4,
      linetype = "dotted"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.09,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(ggplot2::vars(panel), ncol = 2, scales = "free") +
    ggplot2::scale_colour_manual(values = increment_colours) +
    ggplot2::scale_fill_manual(values = increment_colours) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_width(0.4),
      expand = ggplot2::expansion(mult = c(0.015, 0.025))
    ) +
    ggplot2::labs(
      title = "Sensitivity analyses of the protein dose-response shape",
      subtitle = paste(
        "Updated posterior change from 0.8 g/kg/day; shaded regions are",
        "95% credible intervals"
      ),
      x = "Protein intake (g/kg/day)",
      y = "Change in standardised outcome from 0.8 g/kg/day",
      colour = NULL,
      fill = NULL,
      caption = paste(
        "Displayed ranges are restricted to observed treatment-intake",
        "support in the strict within-study data."
      )
    ) +
    base_theme

  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    absolute_data = file.path(data_dir, "absolute_overlay_data.csv.gz"),
    increment_data = file.path(data_dir, "reference_increment_overlay_data.csv.gz"),
    checks = file.path(data_dir, "checks.csv"),
    absolute_png = file.path(figure_dir, "primary_mundlak_absolute_overlay.png"),
    absolute_pdf = file.path(figure_dir, "primary_mundlak_absolute_overlay.pdf"),
    increment_png = file.path(
      figure_dir, "primary_mundlak_within_increment_overlay.png"
    ),
    increment_pdf = file.path(
      figure_dir, "primary_mundlak_within_increment_overlay.pdf"
    )
  )
  readr::write_csv(figure_data$absolute, paths[["absolute_data"]])
  readr::write_csv(figure_data$increments, paths[["increment_data"]])
  readr::write_csv(figure_data$checks, paths[["checks"]])
  ggplot2::ggsave(
    paths[["absolute_png"]], absolute_plot,
    width = 12, height = 14, units = "in", dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    paths[["absolute_pdf"]], absolute_plot,
    width = 12, height = 14, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    paths[["increment_png"]], increment_plot,
    width = 12, height = 14, units = "in", dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    paths[["increment_pdf"]], increment_plot,
    width = 12, height = 14, units = "in", bg = "white"
  )
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    stop("An updated cross-family figure artifact is missing.",
         call. = FALSE)
  }
  manifest <- tibble::tibble(
    artifact = names(paths),
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(paths)$size),
    sha256 = purrr::map_chr(
      paths,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    )
  )
  manifest_path <- file.path(data_dir, "manifest.csv")
  readr::write_csv(manifest, manifest_path)
  normalizePath(
    c(unname(paths), manifest_path),
    winslash = "/",
    mustWork = TRUE
  )
}

make_update_2026_within_pruning_plan <- function(
    accepted_fit_manifest,
    verification_files,
    figure_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "within"
    )
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "within_update_2026"),
    "manifest.csv"
  )
  if (
    length(verification_files) != 35L ||
      any(!file.exists(verification_files))
  ) {
    stop("The updated within-study verification bundle is incomplete.",
         call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)
  figure_files <- unname(unlist(figure_files, use.names = FALSE))
  if (
    length(figure_files) != 8L || anyNA(figure_files) ||
      any(!file.exists(figure_files)) ||
      any(file.info(figure_files)$size <= 0)
  ) {
    stop("The updated cross-family figure archive is incomplete.",
         call. = FALSE)
  }

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "within"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated within-study archive.",
         call. = FALSE)
  }

  fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  if (any(!dir.exists(fit_roots))) {
    stop("An expected updated within-study fit directory is missing.",
         call. = FALSE)
  }
  fit_paths <- unlist(
    purrr::map(
      fit_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  )
  fit_paths <- normalizePath(fit_paths, winslash = "/", mustWork = TRUE)
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(fit_paths) != 108L || length(accepted_paths) != 104L ||
      anyDuplicated(fit_paths) || anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% fit_paths) ||
      sum(!fit_paths %in% accepted_paths) != 4L
  ) {
    stop("Updated within-study fit retention does not match its lineage.",
         call. = FALSE)
  }

  fit_manifest <- tibble::tibble(
    path = fit_paths,
    relative_path = substring(fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", fit_paths, fixed = TRUE) ~ "initial",
      grepl("/remedies/remedy_1/", fit_paths, fixed = TRUE) ~ "remedy_1",
      TRUE ~ "unknown"
    ),
    accepted_fit = fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(fit_paths)$size),
    lean_action = "delete frozen full fit archive"
  ) |>
    dplyr::arrange(archive_stage, relative_path)
  if (
    any(fit_manifest$archive_stage == "unknown") ||
      !identical(
        as.integer(table(factor(
          fit_manifest$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 4L)
      ) ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "initial") != 100L ||
      sum(fit_manifest$accepted_fit &
            fit_manifest$archive_stage == "remedy_1") != 4L
  ) {
    stop("The updated within-study archive stages are not as expected.",
         call. = FALSE)
  }

  compact_roots <- file.path(archive_root, c("loo", "reporting"))
  if (any(!dir.exists(compact_roots))) {
    stop("An expected compact within-study archive directory is missing.",
         call. = FALSE)
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(.x, recursive = TRUE, full.names = TRUE)
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[file.exists(compact_paths)]
  compact_manifest <- tibble::tibble(
    path = normalizePath(compact_paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(compact_paths)$size),
    action = "retain compact LOO/reporting archive"
  )
  if (nrow(compact_manifest) != 204L) {
    stop("The compact within-study archive inventory is incomplete.",
         call. = FALSE)
  }
  figure_manifest <- tibble::tibble(
    path = normalizePath(figure_files, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(figure_files)$size),
    sha256 = purrr::map_chr(
      figure_files,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    ),
    action = "retain updated cross-family figure artifact"
  )

  summary <- tibble::tibble(
    mode = "lean",
    fit_files_deleted = nrow(fit_manifest),
    accepted_fit_files_deleted = sum(fit_manifest$accepted_fit),
    superseded_fit_files_deleted = sum(!fit_manifest$accepted_fit),
    fit_files_retained = 0L,
    gib_reclaimed = sum(fit_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained = sum(compact_manifest$bytes) / 1024^3,
    figure_files_retained = nrow(figure_manifest),
    figure_archive_mib_retained = sum(figure_manifest$bytes) / 1024^2,
    restoration = paste(
      "Reported results remain frozen; restoring full fits requires",
      "refitting the within-study family."
    )
  )

  list(
    summary = summary,
    fit_file_manifest = fit_manifest,
    compact_archive_manifest = compact_manifest,
    figure_archive_manifest = figure_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_within_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    fit_file_manifest = file.path(output_dir, "fit_file_manifest.csv"),
    compact_archive_manifest = file.path(
      output_dir,
      "compact_archive_manifest.csv"
    ),
    figure_archive_manifest = file.path(
      output_dir,
      "figure_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(plan$fit_file_manifest, paths[["fit_file_manifest"]])
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  readr::write_csv(
    plan$figure_archive_manifest,
    paths[["figure_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_within_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "within", "pruning"
    )
) {
  required_confirmation <- "DELETE ALL UPDATED WITHIN FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "within"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the expected within-study archive.",
         call. = FALSE)
  }
  selected <- plan$fit_file_manifest
  if (
    nrow(selected) != 108L || sum(selected$accepted_fit) != 104L ||
      sum(!selected$accepted_fit) != 4L ||
      !identical(
        as.integer(table(factor(
          selected$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 4L)
      )
  ) {
    stop("The within-study pruning set is inconsistent with its audit.",
         call. = FALSE)
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  if (file.exists(receipt_path)) {
    stop("The updated within-study lean-pruning receipt already exists.",
         call. = FALSE)
  }
  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  permitted_subtree <- grepl(
    "/(candidates|remedies)/",
    selected_paths
  )
  if (
    any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted_subtree) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the permitted within-study subtree.",
         call. = FALSE)
  }
  if (
    any(!file.exists(plan$compact_archive_manifest$path)) ||
      any(!file.exists(plan$figure_archive_manifest$path)) ||
      any(
        purrr::map_chr(
          plan$figure_archive_manifest$path,
          ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
        ) != plan$figure_archive_manifest$sha256
      )
  ) {
    stop("A retained compact or figure artifact failed revalidation.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop("At least one within-study fit archive could not be deleted.",
         call. = FALSE)
  }
  if (
    any(!file.exists(plan$compact_archive_manifest$path)) ||
      any(!file.exists(plan$figure_archive_manifest$path)) ||
      !file.exists(plan$verification_manifest)
  ) {
    stop("A retained within-study artifact is unexpectedly missing.",
         call. = FALSE)
  }

  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, receipt_path)
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

# Combined-data protein-manipulation strategy sensitivity refit ----

make_update_2026_strategy_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated strategy base registries are incomplete.",
         call. = FALSE)
  }
  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  primary_remedy_keys <- unname(purrr::map_chr(
    primary_remedy_specifications,
    key
  ))
  # Accepted pre-update remedy schedule carried into the 2026 fit recipe.
  historical_remedy_keys <- c(
    "whole_body::11",
    "hypertrophy::14", "hypertrophy::19", "hypertrophy::26",
    "strength::3", "strength::10", "strength::15", "strength::22",
    "function::11", "function::14", "function::23", "function::26"
  )
  if (
    anyDuplicated(primary_keys) || anyDuplicated(primary_remedy_keys) ||
      !identical(primary_remedy_keys, "function::10") ||
      !all(primary_remedy_keys %in% primary_keys) ||
      !all(historical_remedy_keys %in% primary_keys)
  ) {
    stop("The historical strategy sampling lineage is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(primary_remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  specifications <- purrr::map(specifications, function(specification) {
    specification_key <- key(specification)
    historical_remedy <- specification_key %in% historical_remedy_keys
    if (historical_remedy) {
      specification$adapt_delta <- max(
        as.numeric(specification$adapt_delta),
        0.999
      )
      specification$max_treedepth <- max(
        as.integer(specification$max_treedepth),
        15L
      )
    }
    specification$update_2026_sensitivity_base <- TRUE
    specification$update_2026_strategy_historical_schedule <- if (
      historical_remedy
    ) {
      "longer_sampling"
    } else {
      "base"
    }
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  scheduled_keys <- final_keys[purrr::map_chr(
    specifications,
    "update_2026_strategy_historical_schedule"
  ) == "longer_sampling"]
  function_ten <- specifications[[match("function::10", final_keys)]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !setequal(scheduled_keys, historical_remedy_keys) ||
      !identical(as.numeric(function_ten$adapt_delta), 0.9995)
  ) {
    stop("The updated strategy specification registry failed validation.",
         call. = FALSE)
  }
  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

update_2026_strategy_local_config <- function(specification, config) {
  schedule <- specification$update_2026_strategy_historical_schedule
  if (
    length(schedule) != 1L || is.na(schedule) ||
      !schedule %in% c("base", "longer_sampling")
  ) {
    stop("The updated strategy sampling schedule is invalid.",
         call. = FALSE)
  }
  local_config <- config
  schedule_offset <- 0L
  if (schedule == "longer_sampling") {
    local_config$iter_warmup <- 2000L
    local_config$iter_sampling <- 4000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      15L
    )
    schedule_offset <- 1000L
  }
  outcome_offset <- match(
    specification$outcome,
    c("whole_body", "hypertrophy", "strength", "function")
  )
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the updated strategy configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed + schedule_offset +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

prepare_update_2026_strategy_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  local_config <- update_2026_strategy_local_config(
    specification,
    config
  )
  prepare_protein_strategy_full_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = local_config,
    family = "strategy"
  )
}

validate_update_2026_strategy_inputs <- function(
    inputs,
    config,
    audit_review
) {
  if (!isTRUE(audit_review$audit_ready)) {
    stop(
      "The updated protein-strategy classifications are not confirmed.",
      call. = FALSE
    )
  }
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    local_config <- update_2026_strategy_local_config(
      specification,
      config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    includes_protein <- !as.integer(specification$model_index) %in%
      c(1L, 5L)
    has_strategy_interaction <- grepl(
      "strategy_code:|by\\s*=\\s*strategy_code",
      input$formula
    )
    valid <-
      identical(input$family, "strategy") &&
      is.data.frame(input$data) && nrow(input$data) > 0L &&
      all(c(
        "study_no", "arm_no", "protein_strategy", "strategy_code"
      ) %in% names(input$data)) &&
      all(input$data$strategy_code %in% c(0, 1)) &&
      setequal(unique(input$data$protein_strategy),
               c("Addition", "Substitution")) &&
      grepl("\\bstrategy_code\\b", input$formula) &&
      identical(has_strategy_interaction, includes_protein) &&
      !grepl(
        paste0(
          "strategy_code:resistance_exercise_code:|",
          "resistance_exercise_code:strategy_code:"
        ),
        gsub("\\s+", "", input$formula)
      )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      rows = nrow(input$data),
      studies = dplyr::n_distinct(input$data$study_no),
      arms = dplyr::n_distinct(input$data$arm_no),
      addition_studies = dplyr::n_distinct(
        input$data$study_no[input$data$protein_strategy == "Addition"]
      ),
      substitution_studies = dplyr::n_distinct(
        input$data$study_no[input$data$protein_strategy == "Substitution"]
      ),
      historical_schedule =
        specification$update_2026_strategy_historical_schedule,
      adapt_delta = max(
        as.numeric(specification$adapt_delta),
        as.numeric(local_config$adapt_delta_floor)
      ),
      max_treedepth = max(
        as.integer(specification$max_treedepth),
        as.integer(local_config$max_treedepth_floor)
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      chains = as.integer(local_config$chains),
      parallel_chains = as.integer(local_config$parallel_chains),
      valid = valid
    )
  }) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  counts <- table(factor(checks$outcome, levels = outcomes))
  schedules <- table(factor(
    checks$historical_schedule,
    levels = c("base", "longer_sampling")
  ))
  samples <- checks |>
    dplyr::distinct(
      outcome, rows, studies, arms,
      addition_studies, substitution_studies
    )
  base <- checks |>
    dplyr::filter(historical_schedule == "base")
  longer <- checks |>
    dplyr::filter(historical_schedule == "longer_sampling")
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      !identical(as.integer(schedules), c(92L, 12L)) ||
      nrow(samples) != 4L ||
      any(samples$rows < 1L) || any(samples$studies < 2L) ||
      any(samples$addition_studies < 1L) ||
      any(samples$substitution_studies < 1L) ||
      any(base$iter_warmup != 1500L) ||
      any(base$iter_sampling != 2000L) ||
      any(longer$iter_warmup != 2000L) ||
      any(longer$iter_sampling != 4000L) ||
      any(longer$adapt_delta < 0.999) ||
      any(longer$max_treedepth < 15L) ||
      config$chains != 4L || config$parallel_chains != 4L
  ) {
    stop("The updated protein-strategy input preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_strategy_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated strategy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_strategy_candidate_to_file <- function(
    input,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "strategy",
      "candidates"
    )
) {
  local_config <- update_2026_strategy_local_config(
    input$specification,
    config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "strategy"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated strategy fit file.", call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated strategy fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

write_update_2026_strategy_diagnostic_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy",
      "diagnostics"
    )
) {
  required_summary_columns <- c(
    "outcome", "model_index", "diagnostic_pass", "status"
  )
  required_failure_columns <- c(
    "outcome", "model_index", "variable", "rhat", "ess_bulk",
    "ess_tail", "mcse_sd_ratio"
  )
  if (
    !identical(audit$family, "strategy") ||
      !is.data.frame(audit$diagnostic_summary) ||
      nrow(audit$diagnostic_summary) != 104L ||
      !all(required_summary_columns %in% names(audit$diagnostic_summary)) ||
      anyDuplicated(paste(
        audit$diagnostic_summary$outcome,
        audit$diagnostic_summary$model_index
      )) ||
      anyNA(audit$diagnostic_summary$diagnostic_pass) ||
      !is.data.frame(audit$parameter_failures) ||
      !all(required_failure_columns %in% names(audit$parameter_failures))
  ) {
    stop("Cannot write an incomplete updated strategy diagnostic audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    parameter_failures = file.path(output_dir, "parameter_failures.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]],
    na = ""
  )
  readr::write_csv(
    audit$parameter_failures,
    paths[["parameter_failures"]],
    na = ""
  )
  unname(paths)
}

make_update_2026_strategy_remedy_config <- function(config) {
  make_protein_strategy_full_remedy_config(
    full_config = config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    adapt_delta = 0.999,
    max_treedepth = 15L,
    seed_offset = 10000L
  )
}

select_update_2026_strategy_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  selected <- select_protein_strategy_full_remedy_specifications(
    specifications = specifications,
    initial_audit = diagnostic_audit,
    family = "strategy"
  )
  expected_keys <- c(
    "whole_body::6", "whole_body::25",
    "hypertrophy::11", "hypertrophy::17",
    "strength::1", "strength::7", "strength::9",
    "strength::17", "strength::20",
    "function::8", "function::11", "function::12"
  )
  selected_keys <- unname(purrr::map_chr(selected, function(specification) {
    paste(specification$outcome, specification$model_index, sep = "::")
  }))
  if (
    length(selected) != 12L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop("The updated strategy remedy set no longer matches its audit.",
         call. = FALSE)
  }
  selected
}

update_2026_strategy_remedy_local_config <- function(
    specification,
    config
) {
  local_config <- config
  if (identical(
    specification$update_2026_strategy_historical_schedule,
    "longer_sampling"
  )) {
    local_config$iter_warmup <- 3000L
    local_config$iter_sampling <- 6000L
  }
  outcome_offset <- match(
    specification$outcome,
    c("whole_body", "hypertrophy", "strength", "function")
  )
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the updated strategy remedy configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed +
      1000L * identical(
        specification$update_2026_strategy_historical_schedule,
        "longer_sampling"
      ) +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

prepare_update_2026_strategy_remedy_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  local_config <- update_2026_strategy_remedy_local_config(
    specification,
    config
  )
  prepare_protein_strategy_full_remedy_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = local_config,
    family = "strategy"
  )
}

validate_update_2026_strategy_remedy_inputs <- function(
    inputs,
    initial_manifest,
    diagnostic_audit,
    config
) {
  expected_keys <- c(
    "whole_body::6", "whole_body::25",
    "hypertrophy::11", "hypertrophy::17",
    "strength::1", "strength::7", "strength::9",
    "strength::17", "strength::20",
    "function::8", "function::11", "function::12"
  )
  failed_keys <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
    dplyr::pull(key)
  manifest_keys <- paste(
    initial_manifest$outcome,
    initial_manifest$model_index,
    sep = "::"
  )
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    local_config <- update_2026_strategy_remedy_local_config(
      specification,
      config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    key <- paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      remedy = if (
        identical(
          specification$update_2026_strategy_historical_schedule,
          "longer_sampling"
        )
      ) {
        "extended_sampling"
      } else {
        "longer_sampling"
      },
      adapt_delta = max(
        as.numeric(specification$adapt_delta),
        as.numeric(local_config$adapt_delta_floor)
      ),
      max_treedepth = max(
        as.integer(specification$max_treedepth),
        as.integer(local_config$max_treedepth_floor)
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      chains = as.integer(local_config$chains),
      parallel_chains = as.integer(local_config$parallel_chains),
      initial_fit_exists = key %in% manifest_keys &&
        file.exists(initial_manifest$fit_path[match(key, manifest_keys)]),
      valid = identical(input$family, "strategy") &&
        identical(
          specification$sensitivity_strategy_full_remedy,
          "longer_sampling_high_adapt_delta"
        ) &&
        is.data.frame(input$data) && nrow(input$data) > 0L
    )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  extended <- checks |>
    dplyr::filter(.data$remedy == "extended_sampling")
  standard <- checks |>
    dplyr::filter(.data$remedy == "longer_sampling")
  if (
    !identical(diagnostic_audit$family, "strategy") ||
      nrow(initial_manifest) != 104L || anyDuplicated(manifest_keys) ||
      length(inputs) != 12L || nrow(checks) != 12L ||
      anyDuplicated(keys) || !setequal(keys, expected_keys) ||
      !setequal(failed_keys, expected_keys) ||
      !all(checks$initial_fit_exists) || !all(checks$valid) ||
      nrow(standard) != 11L ||
      any(standard$iter_warmup != 2000L) ||
      any(standard$iter_sampling != 4000L) ||
      nrow(extended) != 1L ||
      extended$outcome[[1L]] != "function" ||
      extended$model_index[[1L]] != 11L ||
      extended$iter_warmup[[1L]] != 3000L ||
      extended$iter_sampling[[1L]] != 6000L ||
      any(checks$adapt_delta < 0.999) ||
      any(checks$max_treedepth < 15L) ||
      any(checks$chains != 4L) || any(checks$parallel_chains != 4L)
  ) {
    stop("The updated strategy remedy preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_strategy_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy",
      "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(checks) != 12L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated strategy remedy preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_strategy_remedy_to_file <- function(
    input,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "strategy",
      "remedies", "remedy_1"
    )
) {
  local_config <- update_2026_strategy_remedy_local_config(
    input$specification,
    config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "strategy"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated strategy remedy fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated strategy remedy fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

write_update_2026_strategy_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy",
      "accepted_fit_manifest.csv"
    )
) {
  stage_counts <- table(factor(
    manifest$stage,
    levels = c("initial", "remedy_1")
  ))
  if (
    !is.data.frame(manifest) || nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(manifest$family != "strategy") ||
      !identical(as.integer(stage_counts), c(92L, 12L)) ||
      any(!file.exists(manifest$fit_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("The updated strategy accepted-fit manifest is incomplete.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

validate_update_2026_strategy_loo_preflight <- function(
    accepted_manifest,
    final_diagnostic_audit
) {
  stage_counts <- table(factor(
    accepted_manifest$stage,
    levels = c("initial", "remedy_1")
  ))
  checks <- tibble::tibble(
    check = c(
      "all 104 final sampling diagnostics pass",
      "accepted registry contains 26 models per outcome",
      "accepted registry contains 92 initial and 12 remedy fits",
      "all accepted fit files exist and are nonempty",
      "all accepted fit hashes are present"
    ),
    passed = c(
      isTRUE(final_diagnostic_audit$all_sampling_diagnostics_pass) &&
        final_diagnostic_audit$models == 104L,
      identical(
        as.integer(table(factor(
          accepted_manifest$outcome,
          levels = c("whole_body", "hypertrophy", "strength", "function")
        ))),
        rep(26L, 4L)
      ),
      identical(as.integer(stage_counts), c(92L, 12L)),
      nrow(accepted_manifest) == 104L &&
        all(file.exists(accepted_manifest$fit_path)) &&
        all(file.info(accepted_manifest$fit_path)$size > 0),
      nrow(accepted_manifest) == 104L &&
        all(nchar(accepted_manifest$sha256) == 64L)
    )
  )
  if (
    anyDuplicated(paste(
      accepted_manifest$outcome,
      accepted_manifest$model_index
    )) ||
      !all(checks$passed)
  ) {
    stop("The updated strategy LOO preflight is incomplete.",
         call. = FALSE)
  }
  list(
    family = "strategy",
    models = nrow(accepted_manifest),
    checks = checks,
    accepted_fit_manifest = accepted_manifest,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores_per_branch = 1L,
    ready = TRUE
  )
}

write_update_2026_strategy_loo_preflight <- function(
    preflight,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "loo_psis"
    )
) {
  if (
    !identical(preflight$family, "strategy") ||
      preflight$models != 104L || !isTRUE(preflight$ready) ||
      !all(preflight$checks$passed)
  ) {
    stop("Cannot write an invalid updated strategy LOO preflight.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    checks = file.path(output_dir, "preflight_checks.csv"),
    accepted_fit_manifest = file.path(
      output_dir,
      "accepted_fit_manifest.csv"
    )
  )
  readr::write_csv(preflight$checks, paths[["checks"]])
  readr::write_csv(
    preflight$accepted_fit_manifest,
    paths[["accepted_fit_manifest"]]
  )
  unname(paths)
}

write_update_2026_strategy_integrated_loo_file <- function(
    result,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "strategy",
      "loo_psis"
    )
) {
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(result$remedy, "integrated_observation_level_psis") ||
      !identical(unname(as.character(result$family)), "strategy")
  ) {
    stop("Cannot archive an invalid updated strategy LOO result.",
         call. = FALSE)
  }
  output_dir <- file.path(output_root, result$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(result$model_index))
  )
  write_atomic_qs(result, output_file)
}

write_update_2026_strategy_integrated_loo_audit <- function(
    audit,
    manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "loo_psis"
    )
) {
  if (
    !identical(audit$family, "strategy") ||
      !isTRUE(audit$integrated_psis_complete) ||
      audit$models != 104L || nrow(manifest) != 104L ||
      any(!file.exists(manifest$integrated_loo_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("The updated strategy integrated-LOO audit is incomplete.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

validate_update_2026_strategy_exact_inputs <- function(
    records,
    integrated_audit,
    config
) {
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    !is.list(records) || length(records) != expected || expected != 16L
  ) {
    stop("Expected 16 updated strategy exact-LOO records.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    integrated <- qs::qread(record$integrated_loo_path)
    observation <- as.integer(record$observation)
    if (
      !inherits(integrated$loo, "loo") ||
        !identical(unname(as.character(integrated$family)), "strategy") ||
        !identical(
          unname(as.character(integrated$outcome)),
          unname(as.character(record$outcome))
        ) ||
        as.integer(integrated$model_index) !=
          as.integer(record$model_index) ||
        observation < 1L ||
        observation > nrow(integrated$observation_manifest)
    ) {
      stop("An updated strategy exact record is invalid.",
           call. = FALSE)
    }
    held_out <- integrated$observation_manifest[
      observation,
      ,
      drop = FALSE
    ]
    observed_k <- as.numeric(
      loo::pareto_k_values(integrated$loo)[[observation]]
    )
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = unname(as.character(record$outcome)),
      model_index = as.integer(record$model_index),
      model_name = unname(as.character(record$model_name)),
      observation = observation,
      study_no = as.character(held_out$study_no[[1L]]),
      arm_no = as.character(held_out$arm_no[[1L]]),
      outcome_no = as.character(held_out$outcome_no[[1L]]),
      pareto_k = observed_k,
      adapt_delta = as.numeric(config$adapt_delta_floor),
      max_treedepth = as.integer(config$max_treedepth_floor),
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      chains = as.integer(config$chains),
      parallel_chains = as.integer(config$parallel_chains),
      original_fit_exists = file.exists(record$original_fit_path),
      integrated_loo_exists = file.exists(record$integrated_loo_path),
      valid = is.finite(observed_k) && observed_k > 0.7
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  if (
    nrow(checks) != expected || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(expected)) ||
      !identical(as.integer(outcome_counts), c(7L, 7L, 0L, 2L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$integrated_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.999 ||
      config$max_treedepth_floor != 15L
  ) {
    stop("The updated strategy exact-LOO preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_strategy_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy", "loo_exact",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 16L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated strategy exact preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_strategy_exact_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "strategy",
    "loo_exact", "fits", unname(as.character(record$outcome))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

fit_update_2026_strategy_exact_to_file <- function(record, config) {
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !inherits(candidate$model, "brmsfit") ||
      !identical(unname(as.character(candidate$family)), "strategy") ||
      !identical(
        unname(as.character(candidate$outcome)),
        unname(as.character(record$outcome))
      ) ||
      as.integer(candidate$model_index) != as.integer(record$model_index) ||
      !inherits(integrated$loo, "loo") ||
      !identical(unname(as.character(integrated$family)), "strategy") ||
      !identical(
        unname(as.character(integrated$outcome)),
        unname(as.character(record$outcome))
      ) ||
      as.integer(integrated$model_index) != as.integer(record$model_index) ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      loo::pareto_k_values(integrated$loo)[[observation]] <= 0.7
  ) {
    stop("An updated strategy exact-LOO record failed validation.",
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
  exact_specification <- candidate$specification
  exact_specification$adapt_delta <- max(
    candidate$specification$adapt_delta,
    config$adapt_delta_floor
  )
  exact_specification$max_treedepth <- max(
    candidate$specification$max_treedepth,
    config$max_treedepth_floor
  )
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
      adapt_delta = exact_specification$adapt_delta,
      max_treedepth = exact_specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  archive <- list(
    family = "strategy",
    case_id = as.integer(record$case_id),
    outcome = unname(as.character(record$outcome)),
    model_index = as.integer(record$model_index),
    model_name = unname(as.character(record$model_name)),
    observation = observation,
    old_psis_pareto_k = as.numeric(record$pareto_k),
    original_fit_path = unname(as.character(record$original_fit_path)),
    integrated_loo_path = unname(as.character(record$integrated_loo_path)),
    prediction_data = prediction_data,
    specification = exact_specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  write_atomic_qs(
    archive,
    update_2026_strategy_exact_file_path(record)
  )
}

write_update_2026_strategy_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "loo_exact",
      "diagnostics"
    )
) {
  if (
    !is.data.frame(audit$diagnostic_summary) ||
      nrow(audit$diagnostic_summary) != 16L ||
      audit$exact_refits != 16L ||
      audit$exact_refits_passing +
        audit$exact_refits_requiring_remediation != 16L ||
      anyDuplicated(paste(
        audit$diagnostic_summary$outcome,
        audit$diagnostic_summary$model_index,
        audit$diagnostic_summary$observation
      ))
  ) {
    stop("Cannot write an incomplete updated strategy exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, "diagnostic_summary.csv")
  readr::write_csv(audit$diagnostic_summary, output_file, na = "")
  output_file
}

update_2026_strategy_exact_correction_file_path <- function(correction) {
  observation <- as.integer(
    correction$exact_refit_summary$observation[[1L]]
  )
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "strategy",
    "loo_exact", "corrections", correction$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(correction$model_index),
      as.integer(correction$case_id),
      observation
    )
  )
}

write_update_2026_strategy_exact_correction_file <- function(
    path,
    diagnostic_audit
) {
  correction <- compute_protein_strategy_full_exact_loo_file(
    path = path,
    diagnostic_audit = diagnostic_audit,
    family = "strategy"
  )
  if (
    !identical(unname(as.character(correction$family)), "strategy") ||
      !inherits(correction$loo, "loo") ||
      nrow(correction$exact_refit_summary) != 1L ||
      !isTRUE(correction$exact_refit_summary$diagnostic_pass[[1L]])
  ) {
    stop("The updated strategy exact correction is invalid.",
         call. = FALSE)
  }
  write_atomic_qs(
    correction,
    update_2026_strategy_exact_correction_file_path(correction)
  )
}

audit_update_2026_strategy_exact_corrections <- function(
    paths,
    diagnostic_audit,
    integrated_audit
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    expected != 16L || length(paths) != expected || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0) ||
      diagnostic_audit$exact_refits != expected ||
      !isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)
  ) {
    stop("The updated strategy exact-correction files are incomplete.",
         call. = FALSE)
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  manifest <- purrr::map2_dfr(paths, corrections, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = unname(as.character(result$outcome)),
      model_index = as.integer(result$model_index),
      model_name = unname(as.character(result$model_name)),
      observation =
        as.integer(result$exact_refit_summary$observation[[1L]]),
      correction_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  diagnostic_keys <- paste(
    diagnostic_audit$diagnostic_summary$outcome,
    diagnostic_audit$diagnostic_summary$model_index,
    diagnostic_audit$diagnostic_summary$observation,
    sep = "::"
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys) ||
      !setequal(observed_keys, diagnostic_keys) ||
      !identical(manifest$case_id, seq_len(expected)) ||
      any(!diagnostic_summary$diagnostic_pass) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("The updated strategy exact-correction audit is invalid.",
         call. = FALSE)
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = expected,
    exact_refits_requiring_remediation = 0L,
    all_exact_refit_diagnostics_pass = TRUE
  )
}

write_update_2026_strategy_final_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "loo_exact",
      "final"
    )
) {
  if (
    audit$exact_refits != 16L || audit$exact_refits_passing != 16L ||
      audit$exact_refits_requiring_remediation != 0L ||
      !isTRUE(audit$all_exact_refit_diagnostics_pass) ||
      nrow(audit$diagnostic_summary) != 16L ||
      nrow(audit$correction_manifest) != 16L ||
      any(!file.exists(audit$correction_manifest$correction_path))
  ) {
    stop("Cannot write an incomplete final strategy exact audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

read_update_2026_strategy_exact_corrections <- function(manifest) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 16L ||
      any(!file.exists(manifest$correction_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop("The updated strategy exact-correction manifest is invalid.",
         call. = FALSE)
  }
  corrections <- purrr::map(manifest$correction_path, qs::qread)
  corrections <- corrections[order(purrr::map_int(corrections, "case_id"))]
  if (
    !identical(
      unname(purrr::map_int(corrections, "case_id")),
      seq_len(16L)
    ) ||
      any(!purrr::map_lgl(corrections, ~ inherits(.x$loo, "loo"))) ||
      any(!purrr::map_lgl(
        corrections,
        ~ identical(unname(as.character(.x$family)), "strategy")
      ))
  ) {
    stop("The updated strategy exact corrections are incomplete.",
         call. = FALSE)
  }
  corrections
}

make_update_2026_strategy_final_loo_audit <- function(
    integrated_manifest,
    exact_corrections,
    exact_audit
) {
  if (
    exact_audit$exact_refits != 16L ||
      !isTRUE(exact_audit$all_exact_refit_diagnostics_pass)
  ) {
    stop("Final strategy LOO requires passing exact-refit diagnostics.",
         call. = FALSE)
  }
  final <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest = integrated_manifest,
    exact_results = exact_corrections,
    family = "strategy"
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- integrated_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "integrated_observation_level_psis_plus_exact",
        "integrated_observation_level_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 16L
  ) {
    stop("The final updated strategy LOO audit is not stacking-ready.",
         call. = FALSE)
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 16L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_strategy_final_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 16L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop("Cannot write an incomplete final strategy LOO audit.",
         call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_strategy_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L ||
      anyDuplicated(paste(weights$outcome, weights$model_index)) ||
      length(outcome_sums) != 4L || any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop("The updated strategy stacking weights are invalid.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_strategy_reporting_preflight <- function(
    records,
    data_by_outcome
) {
  if (!is.list(records) || length(records) < 4L) {
    stop("The updated strategy prediction allocation is incomplete.",
         call. = FALSE)
  }
  checks <- purrr::map_dfr(records, function(record) {
    prediction_data <- data_by_outcome[[record$outcome]]
    protein_max <- if (
      is.data.frame(prediction_data) &&
        "protein_intake_raw" %in% names(prediction_data)
    ) {
      max(prediction_data$protein_intake_raw, na.rm = TRUE)
    } else {
      NA_real_
    }
    tibble::tibble(
      outcome = unname(as.character(record$outcome)),
      model_index = as.integer(record$model_index),
      model_name = unname(as.character(record$model_name)),
      family = unname(as.character(record$family)),
      prediction_type = unname(as.character(record$prediction_type)),
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      protein_max = as.numeric(protein_max),
      grid_step = as.numeric(record$grid_step),
      fit_path = unname(as.character(record$fit_path)),
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path) &&
        identical(unname(as.character(record$family)), "strategy") &&
        identical(unname(as.character(record$prediction_type)), "strategy") &&
        nchar(record$fit_sha256) == 64L &&
        is.finite(record$protein_min) && is.finite(protein_max) &&
        protein_max > record$protein_min
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      retained_weight = sum(weight),
      protein_max = dplyr::first(protein_max),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      any(!is.finite(outcome_audit$retained_weight)) ||
      any(outcome_audit$retained_weight <= 0) ||
      any(outcome_audit$retained_weight > 1 + 1e-8) ||
      any(!is.finite(outcome_audit$protein_max)) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop("The updated strategy reporting preflight failed.",
         call. = FALSE)
  }
  checks
}

write_update_2026_strategy_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy", "reporting",
      "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(checks$allocated_draws, checks$outcome, sum)
  outcome_weights <- tapply(checks$weight, checks$outcome, sum)
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      length(outcome_weights) != 4L ||
      any(!is.finite(outcome_weights)) || any(outcome_weights <= 0) ||
      any(outcome_weights > 1 + 1e-8) || !all(checks$valid)
  ) {
    stop("Cannot write an invalid updated strategy reporting preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_strategy_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "strategy",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_strategy_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "reporting"
    )
) {
  if (
    !identical(audit$family, "strategy") ||
      !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L ||
      any(!file.exists(prediction_manifest$prediction_path))
  ) {
    stop("The updated strategy reporting audit is incomplete.",
         call. = FALSE)
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_protein_strategy_full_reporting_file
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir, "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir, "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir, "common_exposure_support.csv"
    ),
    addition_minus_substitution = file.path(
      output_dir, "addition_minus_substitution.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("addition_minus_substitution"),
    paths[["addition_minus_substitution"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_strategy_verification_checkpoint <- function(
    data,
    source_manifest,
    preparation_audit,
    strategy_audit_review,
    specifications,
    accepted_fit_manifest,
    diagnostic_audit,
    integrated_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  data_by_outcome <- data$principal_by_outcome
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome, model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_protein_strategy_full_reporting_file
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(table) {
    table |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_audit$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_name <- paste0(name_alignment$specification_name, " [strategy]")
  weight_sums <- stacking_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      retained_weight = sum(weight),
      .by = outcome
    )
  strategy_values <- unique(analysis_data$protein_strategy)

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "the updated strategy audit is complete and confirmed",
      "all four strategy datasets contain valid classified observations",
      "the combined-data preparation audit remains valid",
      "all five strategy registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist and have checksums",
      "all final model diagnostics pass",
      "integrated LOO contains 104 reductions and 16 audited high-k cases",
      "all 16 exact refits pass and compact corrections exist",
      "final corrected LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four model-averaged strategy reporting audits pass"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      isTRUE(strategy_audit_review$audit_ready) &&
        strategy_audit_review$pending_confirmations == 0L &&
        all(strategy_audit_review$structural_checks$passed),
      identical(names(data_by_outcome), outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L) &&
        isTRUE(data$audit$audit_ready) &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) && all(analysis_data$vi > 0) &&
        setequal(strategy_values, c("Addition", "Substitution")),
      all(preparation_audit$checks$passed),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_name) &&
        all(name_alignment$diagnostic_name == expected_name) &&
        all(name_alignment$loo_name == expected_name) &&
        all(name_alignment$weight_name == expected_name),
      nrow(accepted_fit_manifest) == 104L &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(accepted_fit_manifest$bytes > 0) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      diagnostic_audit$models == 104L &&
        diagnostic_audit$models_passing_sampling == 104L &&
        isTRUE(diagnostic_audit$all_sampling_diagnostics_pass),
      integrated_loo_audit$models == 104L &&
        integrated_loo_audit$exact_refits_required == 16L &&
        nrow(integrated_loo_audit$high_k_cases) == 16L,
      exact_audit$exact_refits == 16L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        nrow(exact_audit$correction_manifest) == 16L &&
        all(file.exists(exact_audit$correction_manifest$correction_path)) &&
        all(exact_audit$correction_manifest$bytes > 0) &&
        all(nchar(exact_audit$correction_manifest$sha256) == 64L),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(prediction_sums$retained_weight > 0) &&
        all(prediction_sums$retained_weight <= 1 + 1e-8) &&
        all(file.exists(prediction_manifest$prediction_path)) &&
        all(prediction_manifest$bytes > 0) &&
        all(nchar(prediction_manifest$sha256) == 64L),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L &&
        all(reporting_audit$summary$intervals_valid) &&
        all(reporting_audit$summary$reference_increment_valid) &&
        all(reporting_audit$summary$probabilities_valid) &&
        all(reporting_audit$summary$diagonal_contrasts_valid)
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated strategy verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    strategy_review_rows = strategy_audit_review$review_rows,
    strategy_review_status = strategy_audit_review$status_summary,
    strategy_structural_checks = strategy_audit_review$structural_checks,
    strategy_data_outcome_summary = data$audit$outcome_summary,
    strategy_data_support = data$audit$support,
    preparation_checks = preparation_audit$checks,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_audit$diagnostic_summary,
    sampling_parameter_failures = diagnostic_audit$parameter_failures,
    initial_integrated_loo_model_summary = integrated_loo_audit$model_summary,
    initial_integrated_loo_outcome_summary =
      integrated_loo_audit$outcome_summary,
    initial_high_pareto_k_cases = integrated_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_correction_manifest = exact_audit$correction_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    addition_minus_substitution =
      bind_report("addition_minus_substitution"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final updated strategy checkpoint checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    strategy_review_rows = "Confirmed updated-search strategy classifications.",
    strategy_review_status = "Classification-review status summary.",
    strategy_structural_checks = "Strategy audit structural invariants.",
    strategy_data_outcome_summary = "Strategy analysis counts by outcome.",
    strategy_data_support = "Observed exposure support by strategy and RT group.",
    preparation_checks = "Combined-data preparation invariants.",
    analysis_data = "Complete-case arm rows with confirmed strategy labels.",
    model_specifications = "Final 104-model strategy specification registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit paths, stages, sizes, and checksums.",
    sampling_diagnostic_summary = "Final model-level sampling diagnostics.",
    sampling_parameter_failures = "Parameter-level diagnostic failures.",
    initial_integrated_loo_model_summary = "Initial integrated PSIS summaries.",
    initial_integrated_loo_outcome_summary = "Initial LOO audit by outcome.",
    initial_high_pareto_k_cases = "The 16 cases selected for exact refitting.",
    exact_refit_diagnostic_summary = "Diagnostics for 16 exact refits.",
    exact_correction_manifest = "Compact exact-correction paths and checksums.",
    final_loo_manifest = "Final exact-corrected LOO lineage.",
    final_loo_model_summary = "Final LOO summary for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    stacking_weights = "Final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    dose_response_predictions = "Model-averaged curves by strategy and RT group.",
    dose_response_increments = "Strategy-specific increments from 0.8 g/kg/day.",
    exposure_support_by_group = "Observed support by strategy and RT group.",
    common_exposure_support = "Support shared by both strategies within RT group.",
    addition_minus_substitution = "Addition-minus-substitution posterior contrasts.",
    pairwise_contrasts = "Strategy-specific pairwise contrasts and probabilities.",
    pairwise_support = "Observed support used for pairwise contrasts.",
    local_equivalence = "Strategy-specific practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    next_step = paste(
      "The updated strategy family is frozen after this verification",
      "checkpoint; its lean-pruning plan can now be audited."
    )
  )
}

make_update_2026_strategy_pruning_plan <- function(
    accepted_fit_manifest,
    exact_fit_files,
    verification_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "strategy"
    )
) {
  verification_manifest <- file.path(
    here::here("results", "verification", "strategy_update_2026"),
    "manifest.csv"
  )
  if (
    length(verification_files) != 35L ||
      any(!file.exists(verification_files))
  ) {
    stop("The updated strategy verification bundle is incomplete.",
         call. = FALSE)
  }
  verify_update_2026_manifest_files(verification_manifest)

  resolved_root <- normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "strategy"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop("Refusing to plan pruning outside the updated strategy archive.",
         call. = FALSE)
  }

  model_fit_roots <- file.path(
    archive_root,
    c("candidates", "remedies")
  )
  if (any(!dir.exists(model_fit_roots))) {
    stop("An expected updated strategy fit directory is missing.",
         call. = FALSE)
  }
  model_fit_paths <- unlist(
    purrr::map(
      model_fit_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  )
  model_fit_paths <- normalizePath(
    model_fit_paths,
    winslash = "/",
    mustWork = TRUE
  )
  accepted_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(model_fit_paths) != 116L ||
      length(accepted_paths) != 104L ||
      anyDuplicated(model_fit_paths) ||
      anyDuplicated(accepted_paths) ||
      !all(accepted_paths %in% model_fit_paths) ||
      sum(!model_fit_paths %in% accepted_paths) != 12L
  ) {
    stop("Updated strategy fit retention does not match its lineage.",
         call. = FALSE)
  }

  model_fit_manifest <- tibble::tibble(
    path = model_fit_paths,
    relative_path = substring(
      model_fit_paths,
      nchar(resolved_root) + 2L
    ),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", model_fit_paths, fixed = TRUE) ~ "initial",
      grepl(
        "/remedies/remedy_1/",
        model_fit_paths,
        fixed = TRUE
      ) ~ "remedy_1",
      TRUE ~ "unknown"
    ),
    outcome = basename(dirname(model_fit_paths)),
    model_index = as.integer(sub(
      "^model_([0-9]+)\\.qs$",
      "\\1",
      basename(model_fit_paths)
    )),
    accepted_fit = model_fit_paths %in% accepted_paths,
    bytes = as.numeric(file.info(model_fit_paths)$size),
    lean_action = "delete frozen strategy model fit archive"
  ) |>
    dplyr::arrange(archive_stage, outcome, model_index)
  if (
    any(model_fit_manifest$archive_stage == "unknown") ||
      anyNA(model_fit_manifest$model_index) ||
      !setequal(model_fit_manifest$outcome, c(
        "whole_body", "hypertrophy", "strength", "function"
      )) ||
      !identical(
        as.integer(table(factor(
          model_fit_manifest$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 12L)
      ) ||
      !identical(
        as.integer(table(factor(
          model_fit_manifest$archive_stage[
            model_fit_manifest$accepted_fit
          ],
          levels = c("initial", "remedy_1")
        ))),
        c(92L, 12L)
      ) ||
      any(model_fit_manifest$bytes <= 0)
  ) {
    stop("The updated strategy model-fit stages are not as expected.",
         call. = FALSE)
  }

  exact_fit_root <- file.path(archive_root, "loo_exact", "fits")
  if (!dir.exists(exact_fit_root)) {
    stop("The updated strategy exact-refit directory is missing.",
         call. = FALSE)
  }
  exact_fit_files <- unname(unlist(exact_fit_files, use.names = FALSE))
  exact_fit_paths <- normalizePath(
    exact_fit_files,
    winslash = "/",
    mustWork = TRUE
  )
  archived_exact_paths <- list.files(
    exact_fit_root,
    pattern = "\\.qs$",
    recursive = TRUE,
    full.names = TRUE
  ) |>
    normalizePath(winslash = "/", mustWork = TRUE)
  exact_names <- basename(exact_fit_paths)
  exact_refit_manifest <- tibble::tibble(
    path = exact_fit_paths,
    relative_path = substring(
      exact_fit_paths,
      nchar(resolved_root) + 2L
    ),
    archive_stage = "observation_exact",
    outcome = basename(dirname(exact_fit_paths)),
    model_index = as.integer(sub(
      "^model_([0-9]+)_.*$", "\\1", exact_names
    )),
    case_id = as.integer(sub(
      "^.*_case_([0-9]+)_.*$", "\\1", exact_names
    )),
    observation = as.integer(sub(
      "^.*_observation_([0-9]+)\\.qs$", "\\1", exact_names
    )),
    bytes = as.numeric(file.info(exact_fit_paths)$size),
    lean_action = "delete frozen strategy exact-refit archive"
  ) |>
    dplyr::arrange(case_id)
  if (
    length(exact_fit_paths) != 16L ||
      anyDuplicated(exact_fit_paths) ||
      !setequal(exact_fit_paths, archived_exact_paths) ||
      anyNA(exact_refit_manifest$model_index) ||
      anyNA(exact_refit_manifest$case_id) ||
      anyNA(exact_refit_manifest$observation) ||
      !identical(exact_refit_manifest$case_id, seq_len(16L)) ||
      any(exact_refit_manifest$bytes <= 0)
  ) {
    stop("The updated strategy exact-refit inventory is incomplete.",
         call. = FALSE)
  }

  deletion_manifest <- dplyr::bind_rows(
    model_fit_manifest |>
      dplyr::transmute(
        path,
        relative_path,
        archive_stage,
        artifact_role = dplyr::if_else(
          accepted_fit,
          "accepted_model_fit",
          "superseded_model_fit"
        ),
        outcome,
        model_index,
        case_id = NA_integer_,
        observation = NA_integer_,
        bytes,
        lean_action
      ),
    exact_refit_manifest |>
      dplyr::transmute(
        path,
        relative_path,
        archive_stage,
        artifact_role = "exact_refit",
        outcome,
        model_index,
        case_id,
        observation,
        bytes,
        lean_action
      )
  ) |>
    dplyr::arrange(
      factor(
        artifact_role,
        levels = c(
          "accepted_model_fit",
          "superseded_model_fit",
          "exact_refit"
        )
      ),
      outcome,
      model_index,
      case_id
    )
  if (
    nrow(deletion_manifest) != 132L ||
      sum(deletion_manifest$artifact_role == "accepted_model_fit") != 104L ||
      sum(
        deletion_manifest$artifact_role == "superseded_model_fit"
      ) != 12L ||
      sum(deletion_manifest$artifact_role == "exact_refit") != 16L ||
      anyDuplicated(deletion_manifest$path) ||
      any(deletion_manifest$bytes <= 0)
  ) {
    stop("The updated strategy deletion inventory is incomplete.",
         call. = FALSE)
  }

  compact_roots <- c(
    loo_psis = file.path(archive_root, "loo_psis"),
    exact_corrections = file.path(
      archive_root, "loo_exact", "corrections"
    ),
    reporting = file.path(archive_root, "reporting")
  )
  if (any(!dir.exists(compact_roots))) {
    stop("An expected compact strategy archive directory is missing.",
         call. = FALSE)
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[
    file.exists(compact_paths) & !file.info(compact_paths)$isdir
  ]
  compact_paths <- normalizePath(
    compact_paths,
    winslash = "/",
    mustWork = TRUE
  )
  compact_manifest <- tibble::tibble(
    path = compact_paths,
    relative_path = substring(
      compact_paths,
      nchar(resolved_root) + 2L
    ),
    archive_component = dplyr::case_when(
      grepl("/loo_psis/", compact_paths, fixed = TRUE) ~ "loo_psis",
      grepl(
        "/loo_exact/corrections/",
        compact_paths,
        fixed = TRUE
      ) ~ "exact_corrections",
      grepl("/reporting/", compact_paths, fixed = TRUE) ~ "reporting",
      TRUE ~ "unknown"
    ),
    bytes = as.numeric(file.info(compact_paths)$size),
    sha256 = purrr::map_chr(
      compact_paths,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    ),
    action = "retain compact strategy LOO/correction/reporting archive"
  ) |>
    dplyr::arrange(archive_component, relative_path)
  if (
    nrow(compact_manifest) != 169L ||
      any(compact_manifest$archive_component == "unknown") ||
      !identical(
        as.integer(table(factor(
          compact_manifest$archive_component,
          levels = c("loo_psis", "exact_corrections", "reporting")
        ))),
        c(104L, 16L, 49L)
      ) ||
      any(compact_manifest$bytes <= 0) ||
      any(nchar(compact_manifest$sha256) != 64L)
  ) {
    stop("The compact strategy archive inventory is incomplete.",
         call. = FALSE)
  }

  summary <- tibble::tibble(
    mode = "lean",
    full_fit_archives_deleted = nrow(deletion_manifest),
    accepted_model_fit_archives_deleted = sum(
      deletion_manifest$artifact_role == "accepted_model_fit"
    ),
    superseded_model_fit_archives_deleted = sum(
      deletion_manifest$artifact_role == "superseded_model_fit"
    ),
    exact_refit_archives_deleted = sum(
      deletion_manifest$artifact_role == "exact_refit"
    ),
    full_fit_archives_retained = 0L,
    gib_reclaimed = sum(deletion_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained =
      sum(compact_manifest$bytes) / 1024^3,
    verification_files_retained = length(verification_files),
    restoration = paste(
      "Reported results, corrected LOO, and all 104 weights remain frozen;",
      "restoring model or exact-refit objects requires refitting the",
      "strategy family."
    )
  )

  list(
    summary = summary,
    model_fit_manifest = model_fit_manifest,
    exact_refit_manifest = exact_refit_manifest,
    deletion_file_manifest = deletion_manifest,
    compact_archive_manifest = compact_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_strategy_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    model_fit_manifest = file.path(
      output_dir, "model_fit_manifest.csv"
    ),
    exact_refit_manifest = file.path(
      output_dir, "exact_refit_manifest.csv"
    ),
    deletion_file_manifest = file.path(
      output_dir, "deletion_file_manifest.csv"
    ),
    compact_archive_manifest = file.path(
      output_dir, "compact_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(
    plan$model_fit_manifest,
    paths[["model_fit_manifest"]]
  )
  readr::write_csv(
    plan$exact_refit_manifest,
    paths[["exact_refit_manifest"]]
  )
  readr::write_csv(
    plan$deletion_file_manifest,
    paths[["deletion_file_manifest"]]
  )
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_strategy_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy", "pruning"
    )
) {
  required_confirmation <- "DELETE ALL UPDATED STRATEGY FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "strategy"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop("The pruning plan does not target the strategy archive.",
         call. = FALSE)
  }
  selected <- plan$deletion_file_manifest
  if (
    nrow(plan$model_fit_manifest) != 116L ||
      sum(plan$model_fit_manifest$accepted_fit) != 104L ||
      sum(!plan$model_fit_manifest$accepted_fit) != 12L ||
      nrow(plan$exact_refit_manifest) != 16L ||
      nrow(selected) != 132L ||
      sum(selected$artifact_role == "accepted_model_fit") != 104L ||
      sum(selected$artifact_role == "superseded_model_fit") != 12L ||
      sum(selected$artifact_role == "exact_refit") != 16L
  ) {
    stop("The strategy pruning set is inconsistent with its audit.",
         call. = FALSE)
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  pending_receipt_path <- file.path(
    receipt_dir,
    "executed_lean.pending.csv"
  )
  if (file.exists(receipt_path) || file.exists(pending_receipt_path)) {
    stop("An updated strategy lean-pruning receipt already exists.",
         call. = FALSE)
  }

  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  deletion_roots <- c(
    file.path(plan$archive_root, "candidates"),
    file.path(plan$archive_root, "remedies"),
    file.path(plan$archive_root, "loo_exact", "fits")
  )
  current_paths <- unlist(
    purrr::map(
      deletion_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  current_paths <- current_paths[
    file.exists(current_paths) & !file.info(current_paths)$isdir
  ]
  current_paths <- normalizePath(
    current_paths,
    winslash = "/",
    mustWork = TRUE
  )
  root_prefix <- paste0(tolower(expected_root), "/")
  permitted_subtree <- grepl(
    "/(candidates|remedies|loo_exact/fits)/",
    selected_paths
  )
  if (
    !setequal(selected_paths, current_paths) ||
      any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted_subtree) ||
      any(dirname(selected_paths) == expected_root)
  ) {
    stop("A selected file falls outside the strategy fit subtrees.",
         call. = FALSE)
  }

  retained <- plan$compact_archive_manifest
  if (
    nrow(retained) != 169L ||
      any(!file.exists(retained$path)) ||
      any(purrr::map_chr(
        retained$path,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ) != retained$sha256)
  ) {
    stop("A retained compact strategy artifact failed revalidation.",
         call. = FALSE)
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, pending_receipt_path)
  pending <- readr::read_csv(
    pending_receipt_path,
    show_col_types = FALSE
  )
  if (
    nrow(pending) != 132L ||
      any(nchar(pending$sha256_before_deletion) != 64L)
  ) {
    stop("The pending strategy pruning receipt is incomplete.",
         call. = FALSE)
  }

  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop("At least one strategy fit archive could not be deleted.",
         call. = FALSE)
  }
  if (
    any(!file.exists(retained$path)) ||
      any(purrr::map_chr(
        retained$path,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ) != retained$sha256) ||
      !file.exists(plan$verification_manifest)
  ) {
    stop("A retained strategy artifact is unexpectedly missing or changed.",
         call. = FALSE)
  }
  if (!file.rename(pending_receipt_path, receipt_path)) {
    stop("Could not finalise the strategy lean-pruning receipt.",
         call. = FALSE)
  }
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}

# Combined-data strategy-modified strict within-study sensitivity refit ----

make_update_2026_strategy_within_specifications <- function(
    primary_specifications,
    primary_remedy_specifications
) {
  if (
    length(primary_specifications) != 104L ||
      length(primary_remedy_specifications) != 1L
  ) {
    stop("The updated strategy-within base registries are incomplete.",
         call. = FALSE)
  }
  key <- function(x) paste(x$outcome, x$model_index, sep = "::")
  primary_keys <- unname(purrr::map_chr(primary_specifications, key))
  primary_remedy_keys <- unname(purrr::map_chr(
    primary_remedy_specifications,
    key
  ))
  # Accepted pre-update remedy schedule carried into the 2026 fit recipe.
  historical_remedy_keys <- c(
    "hypertrophy::10", "hypertrophy::22",
    "strength::9", "strength::10", "strength::12", "strength::13"
  )
  if (
    anyDuplicated(primary_keys) || anyDuplicated(primary_remedy_keys) ||
      !identical(primary_remedy_keys, "function::10") ||
      !all(primary_remedy_keys %in% primary_keys) ||
      !all(historical_remedy_keys %in% primary_keys)
  ) {
    stop("The historical strategy-within sampling lineage is invalid.",
         call. = FALSE)
  }

  specifications <- primary_specifications
  specifications[match(primary_remedy_keys, primary_keys)] <-
    primary_remedy_specifications
  specifications <- purrr::map(specifications, function(specification) {
    specification_key <- key(specification)
    historical_remedy <-
      specification_key %in% historical_remedy_keys
    if (historical_remedy) {
      specification$adapt_delta <- max(
        as.numeric(specification$adapt_delta),
        0.999
      )
      specification$max_treedepth <- max(
        as.integer(specification$max_treedepth),
        15L
      )
    }
    specification$update_2026_sensitivity_base <- TRUE
    specification$update_2026_strategy_within_historical_schedule <- if (
      historical_remedy
    ) {
      "longer_sampling"
    } else {
      "base"
    }
    specification
  })

  final_keys <- unname(purrr::map_chr(specifications, key))
  outcome_counts <- table(factor(
    purrr::map_chr(specifications, "outcome"),
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  scheduled_keys <- final_keys[purrr::map_chr(
    specifications,
    "update_2026_strategy_within_historical_schedule"
  ) == "longer_sampling"]
  function_ten <- specifications[[match("function::10", final_keys)]]
  if (
    length(specifications) != 104L || anyDuplicated(final_keys) ||
      !identical(as.integer(outcome_counts), rep(26L, 4L)) ||
      !setequal(scheduled_keys, historical_remedy_keys) ||
      !identical(as.numeric(function_ten$adapt_delta), 0.9995)
  ) {
    stop("The updated strategy-within specification registry is invalid.",
         call. = FALSE)
  }
  names(specifications) <- paste0(
    purrr::map_chr(specifications, "outcome"),
    "_model_",
    stringr::str_pad(
      purrr::map_int(specifications, "model_index"),
      width = 2L,
      pad = "0"
    )
  )
  specifications
}

update_2026_strategy_within_local_config <- function(
    specification,
    config
) {
  schedule <-
    specification$update_2026_strategy_within_historical_schedule
  if (
    length(schedule) != 1L || is.na(schedule) ||
      !schedule %in% c("base", "longer_sampling")
  ) {
    stop("The updated strategy-within sampling schedule is invalid.",
         call. = FALSE)
  }
  local_config <- config
  schedule_offset <- 0L
  if (schedule == "longer_sampling") {
    local_config$iter_warmup <- 2000L
    local_config$iter_sampling <- 4000L
    local_config$adapt_delta_floor <- max(
      as.numeric(local_config$adapt_delta_floor),
      0.999
    )
    local_config$max_treedepth_floor <- max(
      as.integer(local_config$max_treedepth_floor),
      15L
    )
    schedule_offset <- 1000L
  }
  outcome_offset <- match(
    specification$outcome,
    c("whole_body", "hypertrophy", "strength", "function")
  )
  if (is.na(outcome_offset)) {
    stop("Unknown outcome in the updated strategy-within configuration.",
         call. = FALSE)
  }
  local_config$seed <- as.integer(
    config$seed + schedule_offset +
      100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

prepare_update_2026_strategy_within_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  local_config <- update_2026_strategy_within_local_config(
    specification,
    config
  )
  prepare_protein_strategy_full_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = local_config,
    family = "strategy_within"
  )
}

validate_update_2026_strategy_within_inputs <- function(
    inputs,
    config,
    audit_review
) {
  if (!isTRUE(audit_review$audit_ready)) {
    stop(
      "The updated protein-strategy classifications are not confirmed.",
      call. = FALSE
    )
  }
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    local_config <- update_2026_strategy_within_local_config(
      specification,
      config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    includes_protein <- !as.integer(specification$model_index) %in%
      c(1L, 5L)
    contrast_names <- input$contrast_names
    has_each_strategy_interaction <- if (includes_protein) {
      !is.null(contrast_names) && length(contrast_names) > 0L &&
        all(vapply(
          contrast_names,
          function(name) grepl(
            paste0(name, ":strategy_code"),
            input$formula,
            fixed = TRUE
          ),
          logical(1L)
        ))
    } else {
      is.null(contrast_names)
    }
    matched_pairs <- input$data |>
      dplyr::distinct(study_no, protein_arm_no, comparator_arm_no) |>
      nrow()
    valid <-
      identical(input$family, "strategy_within") &&
      is.data.frame(input$data) && nrow(input$data) > 0L &&
      all(c(
        "study_no", "protein_arm_no", "comparator_arm_no",
        "protein_strategy", "strategy_code",
        "protein_intake_raw_treatment",
        "protein_intake_raw_comparator"
      ) %in% names(input$data)) &&
      all(is.finite(input$data$yi)) &&
      all(is.finite(input$data$vi)) && all(input$data$vi > 0) &&
      all(input$data$strategy_code %in% c(0, 1)) &&
      setequal(
        unique(input$data$protein_strategy),
        c("Addition", "Substitution")
      ) &&
      all(input$data$sensitivity_strategy_sample == "principal") &&
      grepl("\\bstrategy_code\\b", input$formula) &&
      has_each_strategy_interaction &&
      !grepl(
        paste0(
          "strategy_code:resistance_exercise_code:|",
          "resistance_exercise_code:strategy_code:"
        ),
        gsub("\\s+", "", input$formula)
      )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      contrasts = nrow(input$data),
      studies = dplyr::n_distinct(input$data$study_no),
      matched_arm_pairs = matched_pairs,
      addition_studies = dplyr::n_distinct(
        input$data$study_no[
          input$data$protein_strategy == "Addition"
        ]
      ),
      substitution_studies = dplyr::n_distinct(
        input$data$study_no[
          input$data$protein_strategy == "Substitution"
        ]
      ),
      contrast_basis_dimension = if (is.null(contrast_names)) {
        0L
      } else {
        length(contrast_names)
      },
      historical_schedule =
        specification$update_2026_strategy_within_historical_schedule,
      adapt_delta = max(
        as.numeric(specification$adapt_delta),
        as.numeric(local_config$adapt_delta_floor)
      ),
      max_treedepth = max(
        as.integer(specification$max_treedepth),
        as.integer(local_config$max_treedepth_floor)
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      chains = as.integer(local_config$chains),
      parallel_chains = as.integer(local_config$parallel_chains),
      valid = valid
    )
  }) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)

  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  counts <- table(factor(checks$outcome, levels = outcomes))
  schedules <- table(factor(
    checks$historical_schedule,
    levels = c("base", "longer_sampling")
  ))
  samples <- checks |>
    dplyr::distinct(
      outcome, contrasts, studies, matched_arm_pairs,
      addition_studies, substitution_studies
    )
  base <- checks |>
    dplyr::filter(historical_schedule == "base")
  longer <- checks |>
    dplyr::filter(historical_schedule == "longer_sampling")
  expected_longer <- c(
    "hypertrophy::10", "hypertrophy::22",
    "strength::9", "strength::10", "strength::12", "strength::13"
  )
  if (
    nrow(checks) != 104L || anyDuplicated(keys) || !all(checks$valid) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      !identical(as.integer(schedules), c(98L, 6L)) ||
      nrow(samples) != 4L ||
      any(samples$contrasts < 1L) || any(samples$studies < 2L) ||
      any(samples$matched_arm_pairs < 1L) ||
      any(samples$addition_studies < 1L) ||
      any(samples$substitution_studies < 1L) ||
      !setequal(
        paste(longer$outcome, longer$model_index, sep = "::"),
        expected_longer
      ) ||
      any(base$iter_warmup != 1500L) ||
      any(base$iter_sampling != 2000L) ||
      any(longer$iter_warmup != 2000L) ||
      any(longer$iter_sampling != 4000L) ||
      any(longer$adapt_delta < 0.999) ||
      any(longer$max_treedepth < 15L) ||
      any(checks$chains != 4L) ||
      any(checks$parallel_chains != 4L)
  ) {
    stop("The updated strategy-within input preflight is incomplete.",
         call. = FALSE)
  }
  checks
}

write_update_2026_strategy_within_input_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "input_preflight.csv"
    )
) {
  if (nrow(checks) != 104L || !all(checks$valid)) {
    stop("Cannot write an incomplete updated strategy-within preflight.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_strategy_within_candidate_to_file <- function(
    input,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity",
      "strategy_within", "candidates"
    )
) {
  local_config <- update_2026_strategy_within_local_config(
    input$specification,
    config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "strategy_within"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop("Could not finalise an updated strategy-within fit file.",
         call. = FALSE)
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop("Updated strategy-within fit file was not written correctly.",
         call. = FALSE)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

write_update_2026_strategy_within_diagnostic_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "diagnostics"
    )
) {
  required_summary_columns <- c(
    "outcome", "model_index", "diagnostic_pass", "status"
  )
  required_failure_columns <- c(
    "outcome", "model_index", "variable", "rhat", "ess_bulk",
    "ess_tail", "mcse_sd_ratio"
  )
  if (
    !identical(audit$family, "strategy_within") ||
      !is.data.frame(audit$diagnostic_summary) ||
      nrow(audit$diagnostic_summary) != 104L ||
      !all(required_summary_columns %in% names(audit$diagnostic_summary)) ||
      anyDuplicated(paste(
        audit$diagnostic_summary$outcome,
        audit$diagnostic_summary$model_index
      )) ||
      anyNA(audit$diagnostic_summary$diagnostic_pass) ||
      !is.data.frame(audit$parameter_failures) ||
      !all(required_failure_columns %in% names(audit$parameter_failures))
  ) {
    stop(
      "Cannot write an incomplete updated strategy-within diagnostic audit.",
      call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    parameter_failures = file.path(output_dir, "parameter_failures.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]],
    na = ""
  )
  readr::write_csv(
    audit$parameter_failures,
    paths[["parameter_failures"]],
    na = ""
  )
  unname(paths)
}

make_update_2026_strategy_within_remedy_config <- function(config) {
  make_protein_strategy_full_remedy_config(
    full_config = config,
    iter_warmup = 2000L,
    iter_sampling = 4000L,
    adapt_delta = 0.999,
    max_treedepth = 15L,
    seed_offset = 10000L
  )
}

select_update_2026_strategy_within_remedy_specifications <- function(
    specifications,
    diagnostic_audit
) {
  selected <- select_protein_strategy_full_remedy_specifications(
    specifications = specifications,
    initial_audit = diagnostic_audit,
    family = "strategy_within"
  )
  expected_keys <- c("strength::14", "strength::21", "strength::22")
  selected_keys <- unname(purrr::map_chr(selected, function(specification) {
    paste(specification$outcome, specification$model_index, sep = "::")
  }))
  if (
    length(selected) != 3L || anyDuplicated(selected_keys) ||
      !setequal(selected_keys, expected_keys)
  ) {
    stop(
      "The updated strategy-within remedy set no longer matches its audit.",
      call. = FALSE
    )
  }
  selected
}

update_2026_strategy_within_remedy_local_config <- function(
    specification,
    config
) {
  outcome_offset <- match(
    specification$outcome,
    c("whole_body", "hypertrophy", "strength", "function")
  )
  if (is.na(outcome_offset)) {
    stop(
      "Unknown outcome in the updated strategy-within remedy configuration.",
      call. = FALSE
    )
  }
  local_config <- config
  local_config$seed <- as.integer(
    config$seed + 100L * (outcome_offset - 1L) +
      as.integer(specification$model_index)
  )
  local_config
}

prepare_update_2026_strategy_within_remedy_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config
) {
  local_config <- update_2026_strategy_within_remedy_local_config(
    specification,
    config
  )
  prepare_protein_strategy_full_remedy_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = local_config,
    family = "strategy_within"
  )
}

validate_update_2026_strategy_within_remedy_inputs <- function(
    inputs,
    initial_manifest,
    diagnostic_audit,
    config
) {
  expected_keys <- c("strength::14", "strength::21", "strength::22")
  failed_keys <- diagnostic_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass) |>
    dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
    dplyr::pull(key)
  manifest_keys <- paste(
    initial_manifest$outcome,
    initial_manifest$model_index,
    sep = "::"
  )
  checks <- purrr::map_dfr(inputs, function(input) {
    specification <- input$specification
    local_config <- update_2026_strategy_within_remedy_local_config(
      specification,
      config
    )
    brms::make_standata(
      formula = stats::as.formula(input$formula),
      data = input$data,
      prior = input$priors
    )
    key <- paste(
      specification$outcome,
      specification$model_index,
      sep = "::"
    )
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      remedy = "longer_sampling",
      adapt_delta = max(
        as.numeric(specification$adapt_delta),
        as.numeric(local_config$adapt_delta_floor)
      ),
      max_treedepth = max(
        as.integer(specification$max_treedepth),
        as.integer(local_config$max_treedepth_floor)
      ),
      iter_warmup = as.integer(local_config$iter_warmup),
      iter_sampling = as.integer(local_config$iter_sampling),
      chains = as.integer(local_config$chains),
      parallel_chains = as.integer(local_config$parallel_chains),
      initial_fit_exists = key %in% manifest_keys &&
        file.exists(initial_manifest$fit_path[match(key, manifest_keys)]),
      valid = identical(input$family, "strategy_within") &&
        identical(
          specification$sensitivity_strategy_full_remedy,
          "longer_sampling_high_adapt_delta"
        ) &&
        is.data.frame(input$data) && nrow(input$data) > 0L
    )
  }) |>
    dplyr::arrange(
      factor(outcome, levels = c(
        "whole_body", "hypertrophy", "strength", "function"
      )),
      model_index
    )
  keys <- paste(checks$outcome, checks$model_index, sep = "::")
  if (
    !identical(diagnostic_audit$family, "strategy_within") ||
      nrow(initial_manifest) != 104L || anyDuplicated(manifest_keys) ||
      length(inputs) != 3L || nrow(checks) != 3L ||
      anyDuplicated(keys) || !setequal(keys, expected_keys) ||
      !setequal(failed_keys, expected_keys) ||
      !all(checks$initial_fit_exists) || !all(checks$valid) ||
      any(checks$iter_warmup != 2000L) ||
      any(checks$iter_sampling != 4000L) ||
      any(checks$adapt_delta < 0.999) ||
      any(checks$max_treedepth < 15L) ||
      any(checks$chains != 4L) || any(checks$parallel_chains != 4L)
  ) {
    stop(
      "The updated strategy-within remedy preflight is incomplete.",
      call. = FALSE
    )
  }
  checks
}

write_update_2026_strategy_within_remedy_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(checks) != 3L || !all(checks$valid)) {
    stop(
      "Cannot write an incomplete updated strategy-within remedy preflight.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

fit_update_2026_strategy_within_remedy_to_file <- function(
    input,
    config,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity",
      "strategy_within", "remedies", "remedy_1"
    )
) {
  local_config <- update_2026_strategy_within_remedy_local_config(
    input$specification,
    config
  )
  set_rstan_options()
  candidate <- fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = local_config,
    family = "strategy_within"
  )
  candidate$specification$update_2026_iter_warmup <-
    as.integer(local_config$iter_warmup)
  candidate$specification$update_2026_iter_sampling <-
    as.integer(local_config$iter_sampling)

  output_dir <- file.path(output_root, candidate$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(candidate$model_index))
  )
  temporary_file <- paste0(output_file, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  qs::qsave(candidate, temporary_file, preset = "high")
  if (!file.copy(temporary_file, output_file, overwrite = TRUE)) {
    stop(
      "Could not finalise an updated strategy-within remedy fit file.",
      call. = FALSE
    )
  }
  unlink(temporary_file, force = TRUE)
  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop(
      "Updated strategy-within remedy fit file was not written correctly.",
      call. = FALSE
    )
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

write_update_2026_strategy_within_accepted_fit_manifest <- function(
    manifest,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "accepted_fit_manifest.csv"
    )
) {
  stage_counts <- table(factor(
    manifest$stage,
    levels = c("initial", "remedy_1")
  ))
  if (
    !is.data.frame(manifest) || nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      any(manifest$family != "strategy_within") ||
      !identical(as.integer(stage_counts), c(101L, 3L)) ||
      any(!file.exists(manifest$fit_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop(
      "The updated strategy-within accepted-fit manifest is incomplete.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, output_file)
  output_file
}

validate_update_2026_strategy_within_loo_preflight <- function(
    accepted_manifest,
    final_diagnostic_audit
) {
  stage_counts <- table(factor(
    accepted_manifest$stage,
    levels = c("initial", "remedy_1")
  ))
  checks <- tibble::tibble(
    check = c(
      "all 104 final sampling diagnostics pass",
      "accepted registry contains 26 models per outcome",
      "accepted registry contains 101 initial and 3 remedy fits",
      "all accepted fit files exist and are nonempty",
      "all accepted fit hashes are present"
    ),
    passed = c(
      isTRUE(final_diagnostic_audit$all_sampling_diagnostics_pass) &&
        final_diagnostic_audit$models == 104L,
      identical(
        as.integer(table(factor(
          accepted_manifest$outcome,
          levels = c("whole_body", "hypertrophy", "strength", "function")
        ))),
        rep(26L, 4L)
      ),
      identical(as.integer(stage_counts), c(101L, 3L)),
      nrow(accepted_manifest) == 104L &&
        all(file.exists(accepted_manifest$fit_path)) &&
        all(file.info(accepted_manifest$fit_path)$size > 0),
      nrow(accepted_manifest) == 104L &&
        all(nchar(accepted_manifest$sha256) == 64L)
    )
  )
  if (
    anyDuplicated(paste(
      accepted_manifest$outcome,
      accepted_manifest$model_index
    )) ||
      !all(checks$passed)
  ) {
    stop(
      "The updated strategy-within LOO preflight is incomplete.",
      call. = FALSE
    )
  }
  list(
    family = "strategy_within",
    models = nrow(accepted_manifest),
    checks = checks,
    accepted_fit_manifest = accepted_manifest,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores_per_branch = 1L,
    ready = TRUE
  )
}

write_update_2026_strategy_within_loo_preflight <- function(
    preflight,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_psis"
    )
) {
  if (
    !identical(preflight$family, "strategy_within") ||
      preflight$models != 104L || !isTRUE(preflight$ready) ||
      !all(preflight$checks$passed)
  ) {
    stop(
      "Cannot write an invalid updated strategy-within LOO preflight.",
      call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    checks = file.path(output_dir, "preflight_checks.csv"),
    accepted_fit_manifest = file.path(
      output_dir,
      "accepted_fit_manifest.csv"
    )
  )
  readr::write_csv(preflight$checks, paths[["checks"]])
  readr::write_csv(
    preflight$accepted_fit_manifest,
    paths[["accepted_fit_manifest"]]
  )
  unname(paths)
}

write_update_2026_strategy_within_integrated_loo_file <- function(
    result,
    output_root = here::here(
      "model_archive", "update_2026", "sensitivity", "strategy_within",
      "loo_psis"
    )
) {
  if (
    !inherits(result$loo, "loo") ||
      !inherits(result$conditional_loo, "loo") ||
      !identical(result$remedy, "integrated_observation_level_psis") ||
      !identical(unname(as.character(result$family)), "strategy_within")
  ) {
    stop(
      "Cannot archive an invalid updated strategy-within LOO result.",
      call. = FALSE
    )
  }
  output_dir <- file.path(output_root, result$outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_dir,
    sprintf("model_%02d.qs", as.integer(result$model_index))
  )
  write_atomic_qs(result, output_file)
}

write_update_2026_strategy_within_integrated_loo_audit <- function(
    audit,
    manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_psis"
    )
) {
  if (
    !identical(audit$family, "strategy_within") ||
      !isTRUE(audit$integrated_psis_complete) ||
      audit$models != 104L || nrow(manifest) != 104L ||
      any(!file.exists(manifest$integrated_loo_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop(
      "The updated strategy-within integrated-LOO audit is incomplete.",
      call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    reduction_manifest = file.path(output_dir, "reduction_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    high_k_cases = file.path(output_dir, "high_k_cases.csv")
  )
  readr::write_csv(manifest, paths[["reduction_manifest"]])
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(audit$high_k_cases, paths[["high_k_cases"]])
  unname(paths)
}

validate_update_2026_strategy_within_exact_inputs <- function(
    records,
    integrated_audit,
    config
) {
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    !is.list(records) || length(records) != expected || expected != 65L
  ) {
    stop(
      "Expected 65 updated strategy-within exact-LOO records.",
      call. = FALSE
    )
  }
  checks <- purrr::map_dfr(records, function(record) {
    integrated <- qs::qread(record$integrated_loo_path)
    observation <- as.integer(record$observation)
    if (
      !inherits(integrated$loo, "loo") ||
        !identical(
          unname(as.character(integrated$family)),
          "strategy_within"
        ) ||
        !identical(
          unname(as.character(integrated$outcome)),
          unname(as.character(record$outcome))
        ) ||
        as.integer(integrated$model_index) !=
          as.integer(record$model_index) ||
        observation < 1L ||
        observation > nrow(integrated$observation_manifest)
    ) {
      stop(
        "An updated strategy-within exact record is invalid.",
        call. = FALSE
      )
    }
    held_out <- integrated$observation_manifest[
      observation,
      ,
      drop = FALSE
    ]
    observed_k <- as.numeric(
      loo::pareto_k_values(integrated$loo)[[observation]]
    )
    tibble::tibble(
      case_id = as.integer(record$case_id),
      outcome = unname(as.character(record$outcome)),
      model_index = as.integer(record$model_index),
      model_name = unname(as.character(record$model_name)),
      observation = observation,
      study_no = as.character(held_out$study_no[[1L]]),
      arm_no = as.character(held_out$arm_no[[1L]]),
      outcome_no = as.character(held_out$outcome_no[[1L]]),
      pareto_k = observed_k,
      adapt_delta = as.numeric(config$adapt_delta_floor),
      max_treedepth = as.integer(config$max_treedepth_floor),
      iter_warmup = as.integer(config$iter_warmup),
      iter_sampling = as.integer(config$iter_sampling),
      chains = as.integer(config$chains),
      parallel_chains = as.integer(config$parallel_chains),
      original_fit_exists = file.exists(record$original_fit_path),
      integrated_loo_exists = file.exists(record$integrated_loo_path),
      valid = is.finite(observed_k) && observed_k > 0.7
    )
  }) |>
    dplyr::arrange(case_id)

  observed_keys <- paste(
    checks$outcome,
    checks$model_index,
    checks$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  outcome_counts <- table(factor(
    checks$outcome,
    levels = c("whole_body", "hypertrophy", "strength", "function")
  ))
  if (
    nrow(checks) != expected || anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !identical(checks$case_id, seq_len(expected)) ||
      !identical(as.integer(outcome_counts), c(10L, 26L, 0L, 29L)) ||
      !all(checks$original_fit_exists) ||
      !all(checks$integrated_loo_exists) || !all(checks$valid) ||
      config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.999 ||
      config$max_treedepth_floor != 15L
  ) {
    stop(
      "The updated strategy-within exact-LOO preflight is incomplete.",
      call. = FALSE
    )
  }
  checks
}

write_update_2026_strategy_within_exact_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_exact", "input_preflight.csv"
    )
) {
  if (nrow(checks) != 65L || !all(checks$valid)) {
    stop(
      "Cannot write an incomplete updated strategy-within exact preflight.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_strategy_within_exact_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "strategy_within",
    "loo_exact", "fits", unname(as.character(record$outcome))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

fit_update_2026_strategy_within_exact_to_file <- function(
    record,
    config,
    output_file = update_2026_strategy_within_exact_file_path(record)
) {
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  if (
    !inherits(candidate$model, "brmsfit") ||
      !identical(
        unname(as.character(candidate$family)),
        "strategy_within"
      ) ||
      !identical(
        unname(as.character(candidate$outcome)),
        unname(as.character(record$outcome))
      ) ||
      as.integer(candidate$model_index) != as.integer(record$model_index) ||
      !inherits(integrated$loo, "loo") ||
      !identical(
        unname(as.character(integrated$family)),
        "strategy_within"
      ) ||
      !identical(
        unname(as.character(integrated$outcome)),
        unname(as.character(record$outcome))
      ) ||
      as.integer(integrated$model_index) != as.integer(record$model_index) ||
      observation < 1L || observation > nrow(candidate$model$data) ||
      loo::pareto_k_values(integrated$loo)[[observation]] <= 0.7
  ) {
    stop(
      "An updated strategy-within exact-LOO record failed validation.",
      call. = FALSE
    )
  }
  prediction_data <- candidate$model$data
  training_data <- prediction_data[-observation, , drop = FALSE]
  original_model <- candidate$model
  sample_prior <- attr(original_model$prior, "sample_prior")
  original_model$prior <- normalize_primary_brms_priors(
    original_model$prior
  )
  attr(original_model$prior, "sample_prior") <- sample_prior
  exact_specification <- candidate$specification
  exact_specification$adapt_delta <- max(
    candidate$specification$adapt_delta,
    config$adapt_delta_floor
  )
  exact_specification$max_treedepth <- max(
    candidate$specification$max_treedepth,
    config$max_treedepth_floor
  )
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
      adapt_delta = exact_specification$adapt_delta,
      max_treedepth = exact_specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )
  archive <- list(
    family = "strategy_within",
    case_id = as.integer(record$case_id),
    outcome = unname(as.character(record$outcome)),
    model_index = as.integer(record$model_index),
    model_name = unname(as.character(record$model_name)),
    observation = observation,
    old_psis_pareto_k = as.numeric(record$pareto_k),
    original_fit_path = unname(as.character(record$original_fit_path)),
    integrated_loo_path = unname(as.character(record$integrated_loo_path)),
    prediction_data = prediction_data,
    specification = exact_specification,
    refitted_model = refitted_model,
    refit_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
  write_atomic_qs(
    archive,
    output_file
  )
}

write_update_2026_strategy_within_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_exact", "diagnostics"
    )
) {
  if (
    !is.data.frame(audit$diagnostic_summary) ||
      nrow(audit$diagnostic_summary) != 65L ||
      audit$exact_refits != 65L ||
      audit$exact_refits_passing +
        audit$exact_refits_requiring_remediation != 65L ||
      anyDuplicated(paste(
        audit$diagnostic_summary$outcome,
        audit$diagnostic_summary$model_index,
        audit$diagnostic_summary$observation
      ))
  ) {
    stop(
      "Cannot write an incomplete updated strategy-within exact audit.",
      call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, "diagnostic_summary.csv")
  readr::write_csv(audit$diagnostic_summary, output_file, na = "")
  output_file
}

make_update_2026_strategy_within_exact_remedy_config <- function(
    exact_config,
    adapt_delta = 0.9995,
    seed_offset = 1000L
) {
  config <- exact_config
  config$adapt_delta_floor <- max(
    as.numeric(exact_config$adapt_delta_floor),
    as.numeric(adapt_delta)
  )
  config$seed <- as.integer(exact_config$seed + seed_offset)
  if (
    config$chains != 4L || config$parallel_chains != 4L ||
      config$iter_warmup != 2000L || config$iter_sampling != 6000L ||
      config$adapt_delta_floor != 0.9995 ||
      config$max_treedepth_floor != 15L
  ) {
    stop(
      "The updated strategy-within exact-remedy configuration is invalid.",
      call. = FALSE
    )
  }
  config
}

select_update_2026_strategy_within_exact_remedy_record <- function(
    records,
    exact_audit
) {
  failed <- exact_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass)
  if (
    nrow(failed) != 1L || failed$outcome[[1L]] != "hypertrophy" ||
      failed$model_index[[1L]] != 22L ||
      failed$observation[[1L]] != 158L ||
      failed$divergent_transitions[[1L]] != 1L ||
      failed$maximum_treedepth_transitions[[1L]] != 0L ||
      !isTRUE(failed$convergence_pass[[1L]])
  ) {
    stop(
      "The strategy-within exact remedy no longer matches its audit.",
      call. = FALSE
    )
  }
  selected <- purrr::keep(records, function(record) {
    identical(record$outcome, "hypertrophy") &&
      as.integer(record$model_index) == 22L &&
      as.integer(record$observation) == 158L
  })
  if (
    length(selected) != 1L ||
      as.integer(selected[[1L]]$case_id) != 30L
  ) {
    stop(
      "Could not uniquely select strategy-within exact case 30.",
      call. = FALSE
    )
  }
  names(selected) <-
    "case_030_hypertrophy_model_22_observation_0158"
  selected
}

validate_update_2026_strategy_within_exact_remedy_input <- function(
    record,
    config
) {
  if (!is.list(record) || length(record) != 1L) {
    stop(
      "Expected one strategy-within exact-remedy record.",
      call. = FALSE
    )
  }
  record <- record[[1L]]
  candidate <- read_sensitivity_fit_file(record$original_fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  observation <- as.integer(record$observation)
  observed_k <- as.numeric(
    loo::pareto_k_values(integrated$loo)[[observation]]
  )
  effective_adapt_delta <- max(
    as.numeric(candidate$specification$adapt_delta),
    as.numeric(config$adapt_delta_floor)
  )
  effective_max_treedepth <- max(
    as.integer(candidate$specification$max_treedepth),
    as.integer(config$max_treedepth_floor)
  )
  check <- tibble::tibble(
    case_id = as.integer(record$case_id),
    outcome = record$outcome,
    model_index = as.integer(record$model_index),
    model_name = record$model_name,
    observation = observation,
    remedy = "adapt_delta_0.9995_after_exact_refit",
    adapt_delta = effective_adapt_delta,
    max_treedepth = effective_max_treedepth,
    iter_warmup = as.integer(config$iter_warmup),
    iter_sampling = as.integer(config$iter_sampling),
    chains = as.integer(config$chains),
    parallel_chains = as.integer(config$parallel_chains),
    pareto_k = observed_k,
    original_fit_exists = file.exists(record$original_fit_path),
    integrated_loo_exists = file.exists(record$integrated_loo_path),
    valid = TRUE
  )
  if (
    check$case_id[[1L]] != 30L ||
      check$outcome[[1L]] != "hypertrophy" ||
      check$model_index[[1L]] != 22L ||
      check$observation[[1L]] != 158L ||
      check$adapt_delta[[1L]] != 0.9995 ||
      check$max_treedepth[[1L]] != 15L ||
      check$iter_warmup[[1L]] != 2000L ||
      check$iter_sampling[[1L]] != 6000L ||
      check$chains[[1L]] != 4L || check$parallel_chains[[1L]] != 4L ||
      !is.finite(observed_k) || observed_k <= 0.7 ||
      !check$original_fit_exists[[1L]] ||
      !check$integrated_loo_exists[[1L]]
  ) {
    stop(
      "The updated strategy-within exact-remedy preflight is invalid.",
      call. = FALSE
    )
  }
  check
}

write_update_2026_strategy_within_exact_remedy_preflight <- function(
    check,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_exact", "remedy_1", "input_preflight.csv"
    )
) {
  if (nrow(check) != 1L || !isTRUE(check$valid[[1L]])) {
    stop(
      "Cannot write an invalid strategy-within exact-remedy preflight.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(check, output_file)
  output_file
}

update_2026_strategy_within_exact_remedy_file_path <- function(record) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "strategy_within",
    "loo_exact", "remedies", "remedy_1",
    unname(as.character(record$outcome))
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(record$model_index),
      as.integer(record$case_id),
      as.integer(record$observation)
    )
  )
}

audit_update_2026_strategy_within_exact_remedy_diagnostic <- function(
    diagnostics,
    initial_audit
) {
  summary <- dplyr::bind_rows(diagnostics)
  initial_failed <- initial_audit$diagnostic_summary |>
    dplyr::filter(!.data$diagnostic_pass)
  if (
    initial_audit$exact_refits != 65L ||
      initial_audit$exact_refits_requiring_remediation != 1L ||
      nrow(initial_failed) != 1L || nrow(summary) != 1L ||
      summary$outcome[[1L]] != "hypertrophy" ||
      summary$model_index[[1L]] != 22L ||
      summary$observation[[1L]] != 158L
  ) {
    stop(
      "The strategy-within exact-remedy diagnostic is incomplete.",
      call. = FALSE
    )
  }
  list(
    diagnostic_summary = summary,
    exact_refits = 1L,
    exact_refits_passing = sum(summary$diagnostic_pass),
    exact_refits_requiring_remediation = sum(!summary$diagnostic_pass),
    all_exact_refit_diagnostics_pass = all(summary$diagnostic_pass)
  )
}

write_update_2026_strategy_within_exact_remedy_audit <- function(
    audit,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_exact", "remedy_1", "diagnostic_summary.csv"
    )
) {
  if (
    audit$exact_refits != 1L ||
      nrow(audit$diagnostic_summary) != 1L
  ) {
    stop(
      "Cannot write an invalid strategy-within exact-remedy audit.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(audit$diagnostic_summary, output_file, na = "")
  output_file
}

make_update_2026_strategy_within_final_exact_diagnostic_audit <- function(
    initial_audit,
    remedy_audit,
    records,
    initial_paths,
    remedy_path
) {
  initial_paths <- unname(unlist(initial_paths, use.names = FALSE))
  remedy_path <- unname(unlist(remedy_path, use.names = FALSE))
  if (
    initial_audit$exact_refits != 65L ||
      initial_audit$exact_refits_requiring_remediation != 1L ||
      remedy_audit$exact_refits != 1L ||
      !isTRUE(remedy_audit$all_exact_refit_diagnostics_pass) ||
      length(records) != 65L || length(initial_paths) != 65L ||
      length(remedy_path) != 1L ||
      any(!file.exists(initial_paths)) || !file.exists(remedy_path)
  ) {
    stop(
      "The strategy-within exact diagnostics are not ready to combine.",
      call. = FALSE
    )
  }
  record_manifest <- purrr::map2_dfr(
    records,
    initial_paths,
    function(record, path) {
      tibble::tibble(
        case_id = as.integer(record$case_id),
        outcome = unname(as.character(record$outcome)),
        model_index = as.integer(record$model_index),
        model_name = unname(as.character(record$model_name)),
        observation = as.integer(record$observation),
        stage = "observation_exact",
        fit_path = normalizePath(path, winslash = "/", mustWork = TRUE)
      )
    }
  ) |>
    dplyr::arrange(case_id)
  remedy_summary <- remedy_audit$diagnostic_summary
  remedy_key <- paste(
    remedy_summary$outcome,
    remedy_summary$model_index,
    remedy_summary$observation,
    sep = "::"
  )
  expected_remedy_key <- "hypertrophy::22::158"
  diagnostic_key <- function(data) {
    paste(data$outcome, data$model_index, data$observation, sep = "::")
  }
  initial_summary <- initial_audit$diagnostic_summary
  final_summary <- initial_summary |>
    dplyr::filter(
      paste(outcome, model_index, observation, sep = "::") !=
        expected_remedy_key
    ) |>
    dplyr::bind_rows(remedy_summary) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  record_keys <- diagnostic_key(record_manifest)
  remedy_row <- match(expected_remedy_key, record_keys)
  if (
    !identical(remedy_key, expected_remedy_key) ||
      is.na(remedy_row) ||
      as.integer(record_manifest$case_id[[remedy_row]]) != 30L
  ) {
    stop(
      "The strategy-within exact remedy has the wrong identity.",
      call. = FALSE
    )
  }
  record_manifest$stage[[remedy_row]] <- "observation_exact_remedy_1"
  record_manifest$fit_path[[remedy_row]] <- normalizePath(
    remedy_path,
    winslash = "/",
    mustWork = TRUE
  )
  record_manifest <- record_manifest |>
    dplyr::mutate(
      bytes = file.info(fit_path)$size
    )
  final_keys <- diagnostic_key(final_summary)
  stage_counts <- table(factor(
    record_manifest$stage,
    levels = c("observation_exact", "observation_exact_remedy_1")
  ))
  if (
    nrow(final_summary) != 65L || anyDuplicated(final_keys) ||
      !setequal(final_keys, record_keys) ||
      any(!final_summary$diagnostic_pass) ||
      nrow(record_manifest) != 65L || anyDuplicated(record_keys) ||
      !identical(record_manifest$case_id, seq_len(65L)) ||
      !identical(as.integer(stage_counts), c(64L, 1L)) ||
      any(!file.exists(record_manifest$fit_path)) ||
      any(record_manifest$bytes <= 0)
  ) {
    stop(
      "The final strategy-within exact diagnostic audit is incomplete.",
      call. = FALSE
    )
  }
  list(
    family = "strategy_within",
    diagnostic_summary = final_summary,
    result_manifest = record_manifest,
    exact_refits = 65L,
    exact_refits_passing = 65L,
    exact_refits_requiring_remediation = 0L,
    all_exact_refit_diagnostics_pass = TRUE
  )
}

make_update_2026_strategy_within_accepted_exact_fit_records <- function(
    audit
) {
  manifest <- audit$result_manifest
  if (
    audit$exact_refits != 65L ||
      !isTRUE(audit$all_exact_refit_diagnostics_pass) ||
      !is.data.frame(manifest) || nrow(manifest) != 65L ||
      any(!file.exists(manifest$fit_path))
  ) {
    stop(
      "The accepted strategy-within exact-fit registry is incomplete.",
      call. = FALSE
    )
  }
  purrr::map(seq_len(nrow(manifest)), function(row) {
    purrr::map(as.list(manifest[row, , drop = FALSE]), unname)
  })
}

update_2026_strategy_within_exact_correction_file_path <- function(
    correction
) {
  observation <- as.integer(
    correction$exact_refit_summary$observation[[1L]]
  )
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "strategy_within",
    "loo_exact", "corrections", correction$outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  file.path(
    directory,
    sprintf(
      "model_%02d_case_%03d_observation_%04d.qs",
      as.integer(correction$model_index),
      as.integer(correction$case_id),
      observation
    )
  )
}

write_update_2026_strategy_within_exact_correction_file <- function(
    record,
    diagnostic_audit
) {
  correction <- compute_protein_strategy_full_exact_loo_file(
    path = record$fit_path,
    diagnostic_audit = diagnostic_audit,
    family = "strategy_within"
  )
  # The shared strategy exact-correction helper restores the family label for
  # arm-based strategy models, but returns the generic Mundlak result directly
  # for strategy-within models. Restore that identifying metadata before the
  # correction is validated and archived.
  correction$family <- "strategy_within"
  if (
    !identical(
      unname(as.character(correction$family)),
      "strategy_within"
    ) ||
      !inherits(correction$loo, "loo") ||
      nrow(correction$exact_refit_summary) != 1L ||
      !isTRUE(correction$exact_refit_summary$diagnostic_pass[[1L]]) ||
      as.integer(correction$case_id) != as.integer(record$case_id) ||
      as.integer(correction$model_index) != as.integer(record$model_index) ||
      as.integer(correction$exact_refit_summary$observation[[1L]]) !=
        as.integer(record$observation)
  ) {
    stop(
      "The updated strategy-within exact correction is invalid.",
      call. = FALSE
    )
  }
  write_atomic_qs(
    correction,
    update_2026_strategy_within_exact_correction_file_path(correction)
  )
}

audit_update_2026_strategy_within_exact_corrections <- function(
    paths,
    diagnostic_audit,
    integrated_audit
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  expected <- as.integer(integrated_audit$exact_refits_required)
  if (
    expected != 65L || length(paths) != expected || anyNA(paths) ||
      any(!file.exists(paths)) || any(file.info(paths)$size <= 0) ||
      diagnostic_audit$exact_refits != expected ||
      !isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)
  ) {
    stop(
      "The updated strategy-within exact-correction files are incomplete.",
      call. = FALSE
    )
  }
  corrections <- purrr::map(paths, qs::qread)
  diagnostic_summary <- purrr::map_dfr(
    corrections,
    "exact_refit_summary"
  ) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index,
      observation
    )
  manifest <- purrr::map2_dfr(paths, corrections, function(path, result) {
    tibble::tibble(
      case_id = as.integer(result$case_id),
      outcome = unname(as.character(result$outcome)),
      model_index = as.integer(result$model_index),
      model_name = unname(as.character(result$model_name)),
      observation =
        as.integer(result$exact_refit_summary$observation[[1L]]),
      correction_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(path)$size,
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }) |>
    dplyr::arrange(case_id)
  observed_keys <- paste(
    manifest$outcome,
    manifest$model_index,
    manifest$observation,
    sep = "::"
  )
  expected_keys <- paste(
    integrated_audit$high_k_cases$outcome,
    integrated_audit$high_k_cases$model_index,
    integrated_audit$high_k_cases$observation,
    sep = "::"
  )
  diagnostic_keys <- paste(
    diagnostic_audit$diagnostic_summary$outcome,
    diagnostic_audit$diagnostic_summary$model_index,
    diagnostic_audit$diagnostic_summary$observation,
    sep = "::"
  )
  if (
    nrow(diagnostic_summary) != expected || nrow(manifest) != expected ||
      anyDuplicated(observed_keys) ||
      !setequal(observed_keys, expected_keys) ||
      !setequal(observed_keys, diagnostic_keys) ||
      !identical(manifest$case_id, seq_len(expected)) ||
      any(!diagnostic_summary$diagnostic_pass) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop(
      "The updated strategy-within exact-correction audit is invalid.",
      call. = FALSE
    )
  }
  list(
    diagnostic_summary = diagnostic_summary,
    correction_manifest = manifest,
    exact_refits = expected,
    exact_refits_passing = expected,
    exact_refits_requiring_remediation = 0L,
    all_exact_refit_diagnostics_pass = TRUE
  )
}

write_update_2026_strategy_within_final_exact_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_exact", "final"
    )
) {
  if (
    audit$exact_refits != 65L || audit$exact_refits_passing != 65L ||
      audit$exact_refits_requiring_remediation != 0L ||
      !isTRUE(audit$all_exact_refit_diagnostics_pass) ||
      nrow(audit$diagnostic_summary) != 65L ||
      nrow(audit$correction_manifest) != 65L ||
      any(!file.exists(audit$correction_manifest$correction_path))
  ) {
    stop(
      "Cannot write an incomplete final strategy-within exact audit.",
      call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    diagnostic_summary = file.path(output_dir, "diagnostic_summary.csv"),
    correction_manifest = file.path(output_dir, "correction_manifest.csv")
  )
  readr::write_csv(
    audit$diagnostic_summary,
    paths[["diagnostic_summary"]]
  )
  readr::write_csv(
    audit$correction_manifest,
    paths[["correction_manifest"]]
  )
  unname(paths)
}

read_update_2026_strategy_within_exact_corrections <- function(manifest) {
  if (
    !is.data.frame(manifest) || nrow(manifest) != 65L ||
      any(!file.exists(manifest$correction_path)) ||
      any(manifest$bytes <= 0) || any(nchar(manifest$sha256) != 64L)
  ) {
    stop(
      "The updated strategy-within exact-correction manifest is invalid.",
      call. = FALSE
    )
  }
  corrections <- purrr::map(manifest$correction_path, qs::qread)
  corrections <- corrections[order(purrr::map_int(corrections, "case_id"))]
  if (
    !identical(
      unname(purrr::map_int(corrections, "case_id")),
      seq_len(65L)
    ) ||
      any(!purrr::map_lgl(corrections, ~ inherits(.x$loo, "loo"))) ||
      any(!purrr::map_lgl(
        corrections,
        ~ identical(
          unname(as.character(.x$family)),
          "strategy_within"
        )
      ))
  ) {
    stop(
      "The updated strategy-within exact corrections are incomplete.",
      call. = FALSE
    )
  }
  corrections
}

make_update_2026_strategy_within_final_loo_audit <- function(
    integrated_manifest,
    exact_corrections,
    exact_audit
) {
  if (
    exact_audit$exact_refits != 65L ||
      !isTRUE(exact_audit$all_exact_refit_diagnostics_pass)
  ) {
    stop(
      "Final strategy-within LOO requires passing exact diagnostics.",
      call. = FALSE
    )
  }
  final <- apply_sensitivity_full_exact_loo_corrections_from_files(
    integrated_manifest = integrated_manifest,
    exact_results = exact_corrections,
    family = "strategy_within"
  )
  outcome_summary <- final$model_summary |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      models = dplyr::n(),
      models_with_unresolved_high_k =
        sum(n_pareto_k_gt_threshold > 0L),
      unresolved_high_k_cases = sum(n_pareto_k_gt_threshold),
      maximum_remaining_pareto_k = max(maximum_pareto_k, na.rm = TRUE),
      .groups = "drop"
    )
  correction_counts <- exact_audit$diagnostic_summary |>
    dplyr::count(outcome, model_index, name = "exact_corrections")
  final_manifest <- integrated_manifest |>
    dplyr::left_join(
      correction_counts,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      exact_corrections = tidyr::replace_na(exact_corrections, 0L),
      final_method = dplyr::if_else(
        exact_corrections > 0L,
        "integrated_observation_level_psis_plus_exact",
        "integrated_observation_level_psis"
      )
    )
  if (
    !isTRUE(final$stacking_ready) || nrow(final$model_summary) != 104L ||
      nrow(outcome_summary) != 4L ||
      any(outcome_summary$models != 26L) ||
      any(outcome_summary$unresolved_high_k_cases != 0L) ||
      nrow(final_manifest) != 104L ||
      sum(final_manifest$exact_corrections) != 65L
  ) {
    stop(
      "The final updated strategy-within LOO audit is not stacking-ready.",
      call. = FALSE
    )
  }
  final$outcome_summary <- outcome_summary
  final$final_loo_manifest <- final_manifest
  final$exact_refits <- 65L
  final$all_exact_refit_diagnostics_pass <- TRUE
  final
}

write_update_2026_strategy_within_final_loo_audit <- function(
    audit,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "loo_final"
    )
) {
  if (
    !isTRUE(audit$stacking_ready) ||
      nrow(audit$model_summary) != 104L ||
      nrow(audit$outcome_summary) != 4L ||
      nrow(audit$exact_refit_summary) != 65L ||
      nrow(audit$final_loo_manifest) != 104L
  ) {
    stop(
      "Cannot write an incomplete final strategy-within LOO audit.",
      call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    final_loo_manifest = file.path(output_dir, "final_loo_manifest.csv"),
    model_summary = file.path(output_dir, "model_summary.csv"),
    outcome_summary = file.path(output_dir, "outcome_summary.csv"),
    exact_refit_summary = file.path(output_dir, "exact_refit_summary.csv")
  )
  readr::write_csv(
    audit$final_loo_manifest,
    paths[["final_loo_manifest"]]
  )
  readr::write_csv(audit$model_summary, paths[["model_summary"]])
  readr::write_csv(audit$outcome_summary, paths[["outcome_summary"]])
  readr::write_csv(
    audit$exact_refit_summary,
    paths[["exact_refit_summary"]]
  )
  unname(paths)
}

write_update_2026_strategy_within_stacking_weights <- function(
    weights,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "stacking_weights.csv"
    )
) {
  outcome_sums <- tapply(weights$weight, weights$outcome, sum)
  if (
    nrow(weights) != 104L ||
      anyDuplicated(paste(weights$outcome, weights$model_index)) ||
      length(outcome_sums) != 4L || any(abs(outcome_sums - 1) > 1e-8)
  ) {
    stop(
      "The updated strategy-within stacking weights are invalid.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(weights, output_file)
  output_file
}

make_update_2026_strategy_within_reporting_preflight <- function(
    records,
    data_by_outcome
) {
  if (!is.list(records) || length(records) < 4L) {
    stop(
      "The updated strategy-within prediction allocation is incomplete.",
      call. = FALSE
    )
  }
  checks <- purrr::map_dfr(records, function(record) {
    prediction_data <- data_by_outcome[[record$outcome]]
    protein_max <- if (
      is.data.frame(prediction_data) &&
        "protein_intake_raw_treatment" %in% names(prediction_data)
    ) {
      max(prediction_data$protein_intake_raw_treatment, na.rm = TRUE)
    } else {
      NA_real_
    }
    tibble::tibble(
      outcome = unname(as.character(record$outcome)),
      model_index = as.integer(record$model_index),
      model_name = unname(as.character(record$model_name)),
      family = unname(as.character(record$family)),
      prediction_type = unname(as.character(record$prediction_type)),
      weight = as.numeric(record$weight),
      allocated_draws = as.integer(record$draw_count),
      prediction_seed = as.integer(record$prediction_seed),
      protein_min = as.numeric(record$protein_min),
      protein_max = as.numeric(protein_max),
      grid_step = as.numeric(record$grid_step),
      fit_path = unname(as.character(record$fit_path)),
      fit_exists = file.exists(record$fit_path),
      valid = record$draw_count > 0L && file.exists(record$fit_path) &&
        identical(unname(as.character(record$family)), "strategy_within") &&
        identical(
          unname(as.character(record$prediction_type)),
          "strategy_within"
        ) &&
        nchar(record$fit_sha256) == 64L &&
        is.finite(record$protein_min) && is.finite(protein_max) &&
        protein_max > record$protein_min
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_audit <- checks |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      selected_models = dplyr::n(),
      allocated_draws = sum(allocated_draws),
      retained_weight = sum(weight),
      protein_max = dplyr::first(protein_max),
      .groups = "drop"
    )
  if (
    nrow(outcome_audit) != 4L ||
      any(outcome_audit$selected_models < 1L) ||
      any(outcome_audit$allocated_draws != 4000L) ||
      any(!is.finite(outcome_audit$retained_weight)) ||
      any(outcome_audit$retained_weight <= 0) ||
      any(outcome_audit$retained_weight > 1 + 1e-8) ||
      any(!is.finite(outcome_audit$protein_max)) ||
      anyDuplicated(checks[c("outcome", "model_index")]) ||
      !all(checks$valid)
  ) {
    stop(
      "The updated strategy-within reporting preflight failed.",
      call. = FALSE
    )
  }
  checks
}

write_update_2026_strategy_within_reporting_preflight <- function(
    checks,
    output_file = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "reporting", "input_preflight.csv"
    )
) {
  outcome_draws <- tapply(checks$allocated_draws, checks$outcome, sum)
  outcome_weights <- tapply(checks$weight, checks$outcome, sum)
  if (
    length(outcome_draws) != 4L || any(outcome_draws != 4000L) ||
      length(outcome_weights) != 4L ||
      any(!is.finite(outcome_weights)) || any(outcome_weights <= 0) ||
      any(outcome_weights > 1 + 1e-8) || !all(checks$valid)
  ) {
    stop(
      "Cannot write an invalid updated strategy-within reporting preflight.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(checks, output_file)
  output_file
}

update_2026_strategy_within_reporting_path <- function(
    outcome,
    stage,
    model_index = NULL
) {
  directory <- here::here(
    "model_archive", "update_2026", "sensitivity", "strategy_within",
    "reporting", stage, outcome
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- if (is.null(model_index)) {
    paste0(stage, ".qs")
  } else {
    sprintf("model_%02d_%s.qs", as.integer(model_index), stage)
  }
  file.path(directory, filename)
}

write_update_2026_strategy_within_prediction_draw_file <- function(
    record,
    data_by_outcome
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !identical(unname(as.character(record$family)), "strategy_within") ||
      !identical(
        unname(as.character(record$prediction_type)),
        "strategy_within"
      )
  ) {
    stop(
      "An updated strategy-within prediction record is invalid.",
      call. = FALSE
    )
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "strategy_within") ||
      !identical(
        unname(as.character(candidate$outcome)),
        unname(as.character(record$outcome))
      ) ||
      as.integer(candidate$model_index) != as.integer(record$model_index) ||
      !identical(
        unname(as.character(candidate$model_name)),
        unname(as.character(record$model_name))
      )
  ) {
    stop(
      "An updated strategy-within prediction fit does not match its record.",
      call. = FALSE
    )
  }
  prediction_data <- data_by_outcome[[record$outcome]]
  if (is.null(prediction_data) || !is.data.frame(prediction_data)) {
    stop(
      "The updated strategy-within prediction data are unavailable.",
      call. = FALSE
    )
  }
  grid <- make_protein_strategy_within_prediction_grid(
    fitted_candidate = candidate,
    data = prediction_data,
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
    stop(
      "Updated strategy-within prediction draws are incomplete.",
      call. = FALSE
    )
  }
  reduction <- list(
    family = "strategy_within",
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
      prediction_type = "strategy_within",
      reference_intake = 0.8
    )
  )
  write_atomic_qs(
    reduction,
    update_2026_strategy_within_reporting_path(
      outcome = record$outcome,
      stage = "prediction_draws",
      model_index = record$model_index
    )
  )
}

write_update_2026_strategy_within_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path,
    family = "strategy_within"
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (
    is.null(data) || !is.data.frame(data) ||
      !all(c(
        "protein_strategy", "protein_intake_raw_treatment"
      ) %in% names(data))
  ) {
    stop(
      "The updated strategy-within reporting data are unavailable.",
      call. = FALSE
    )
  }
  strategy_summary <- build_protein_strategy_within_pilot_reporting(
    model_average = model_average,
    data = data,
    outcome = outcome
  )
  strategy_reports <- purrr::map(
    c("Substitution", "Addition"),
    function(strategy) {
      strategy_average <- subset_protein_strategy_model_average(
        model_average,
        strategy
      )
      strategy_data <- dplyr::filter(
        data,
        protein_strategy == strategy
      )
      absolute <- summarise_primary_dose_response(strategy_average)
      supported <- build_sensitivity_within_reporting_pilot(
        absolute,
        strategy_average,
        strategy_data,
        reference_intake = reference_intake
      )
      heatmap <- compute_primary_heatmap_contrasts(
        model_average = strategy_average,
        data = dplyr::mutate(
          strategy_data,
          protein_intake_raw = protein_intake_raw_treatment
        ),
        protein_min = reference_intake,
        sesoi = sesoi
      )
      annotate <- function(table) {
        dplyr::mutate(table, protein_strategy = strategy, .before = 1L)
      }
      list(
        increments = annotate(supported$increments_full_grid),
        support = annotate(supported$support_by_group),
        common_support = annotate(supported$common_support),
        pairwise = annotate(heatmap$contrasts),
        pairwise_support = annotate(heatmap$support),
        local_equivalence = annotate(
          summarise_primary_local_equivalence(
            heatmap$contrasts,
            outcome = outcome
          )
        )
      )
    }
  )
  bind_component <- function(name) {
    dplyr::bind_rows(purrr::map(strategy_reports, name))
  }
  report <- list(
    family = "strategy_within",
    outcome = outcome,
    dose_response = strategy_summary$dose_response_contrasts,
    dose_response_increments = bind_component("increments"),
    exposure_support_by_group = strategy_summary$support,
    common_exposure_support = strategy_summary$common_support,
    addition_minus_substitution =
      strategy_summary$addition_minus_substitution,
    pairwise_contrasts = bind_component("pairwise"),
    pairwise_support = bind_component("pairwise_support"),
    local_equivalence = bind_component("local_equivalence"),
    model_weights = tibble::tibble(
      model_name = names(model_average$model_weights),
      weight = as.numeric(model_average$model_weights)
    ),
    settings = utils::modifyList(
      model_average$settings,
      list(reference_intake = reference_intake, sesoi = sesoi)
    )
  )
  write_atomic_qs(
    report,
    update_2026_strategy_within_reporting_path(
      outcome = outcome,
      stage = "reporting"
    )
  )
}

write_update_2026_strategy_within_reporting_audit <- function(
    audit,
    prediction_manifest,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within",
      "reporting"
    )
) {
  if (
    !identical(audit$family, "strategy_within") ||
      !isTRUE(audit$all_reporting_checks_pass) ||
      nrow(audit$summary) != 4L ||
      sum(prediction_manifest$draw_count) != 16000L ||
      any(!file.exists(prediction_manifest$prediction_path))
  ) {
    stop(
      "The updated strategy-within reporting audit is incomplete.",
      call. = FALSE
    )
  }
  reports <- purrr::map(
    audit$summary$reporting_path,
    read_protein_strategy_within_full_reporting_file
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  model_weights <- purrr::map_dfr(reports, function(report) {
    dplyr::mutate(
      report$model_weights,
      outcome = report$outcome,
      .before = 1L
    )
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    prediction_manifest = file.path(output_dir, "prediction_manifest.csv"),
    dose_response = file.path(output_dir, "dose_response.csv"),
    dose_response_increments = file.path(
      output_dir,
      "dose_response_increments.csv"
    ),
    exposure_support_by_group = file.path(
      output_dir,
      "exposure_support_by_group.csv"
    ),
    common_exposure_support = file.path(
      output_dir,
      "common_exposure_support.csv"
    ),
    addition_minus_substitution = file.path(
      output_dir,
      "addition_minus_substitution.csv"
    ),
    pairwise_contrasts = file.path(output_dir, "pairwise_contrasts.csv.gz"),
    pairwise_support = file.path(output_dir, "pairwise_support.csv"),
    local_equivalence = file.path(output_dir, "local_equivalence.csv"),
    model_weights = file.path(output_dir, "model_weights.csv")
  )
  readr::write_csv(audit$summary, paths[["summary"]])
  readr::write_csv(prediction_manifest, paths[["prediction_manifest"]])
  readr::write_csv(bind_report("dose_response"), paths[["dose_response"]])
  readr::write_csv(
    bind_report("dose_response_increments"),
    paths[["dose_response_increments"]]
  )
  readr::write_csv(
    bind_report("exposure_support_by_group"),
    paths[["exposure_support_by_group"]]
  )
  readr::write_csv(
    bind_report("common_exposure_support"),
    paths[["common_exposure_support"]]
  )
  readr::write_csv(
    bind_report("addition_minus_substitution"),
    paths[["addition_minus_substitution"]]
  )
  readr::write_csv(
    bind_report("pairwise_contrasts"),
    paths[["pairwise_contrasts"]]
  )
  readr::write_csv(
    bind_report("pairwise_support"),
    paths[["pairwise_support"]]
  )
  readr::write_csv(
    bind_report("local_equivalence"),
    paths[["local_equivalence"]]
  )
  readr::write_csv(model_weights, paths[["model_weights"]])
  unname(paths)
}

make_update_2026_strategy_within_verification_checkpoint <- function(
    data_by_outcome,
    source_manifest,
    preparation_audit,
    strategy_audit_review,
    specifications,
    accepted_fit_manifest,
    diagnostic_audit,
    integrated_loo_audit,
    exact_audit,
    final_loo_audit,
    stacking_weights,
    prediction_manifest,
    reporting_audit
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  analysis_data <- purrr::imap_dfr(
    data_by_outcome,
    ~ dplyr::mutate(.x, analysis_outcome = .y, .before = 1L)
  )
  data_summary <- analysis_data |>
    dplyr::summarise(
      contrasts = dplyr::n(),
      studies = dplyr::n_distinct(study_no),
      matched_arm_pairs = dplyr::n_distinct(
        paste(study_no, protein_arm_no, comparator_arm_no, sep = "::")
      ),
      .by = c(analysis_outcome, protein_strategy)
    ) |>
    dplyr::arrange(
      factor(analysis_outcome, levels = outcomes),
      protein_strategy
    )
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(
          outcome,
          model_index,
          accepted_stage = stage,
          accepted_fit_sha256 = sha256
        ),
      by = c("outcome", "model_index"),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_protein_strategy_within_full_reporting_file
  )
  bind_report <- function(name) {
    purrr::map_dfr(reports, function(report) {
      dplyr::mutate(report[[name]], outcome = report$outcome, .before = 1L)
    })
  }
  registry_keys <- function(table) {
    table |>
      dplyr::distinct(outcome, model_index) |>
      dplyr::arrange(outcome, model_index) |>
      dplyr::transmute(key = paste(outcome, model_index, sep = "::")) |>
      dplyr::pull(key)
  }
  specification_keys <- registry_keys(specification_table)
  diagnostic_keys <- registry_keys(diagnostic_audit$diagnostic_summary)
  accepted_keys <- registry_keys(accepted_fit_manifest)
  loo_keys <- registry_keys(final_loo_audit$model_summary)
  weight_keys <- registry_keys(stacking_weights)
  common_registry <-
    identical(specification_keys, diagnostic_keys) &&
    identical(specification_keys, accepted_keys) &&
    identical(specification_keys, loo_keys) &&
    identical(specification_keys, weight_keys)
  name_alignment <- specification_table |>
    dplyr::select(outcome, model_index, specification_name = model_name) |>
    dplyr::left_join(
      accepted_fit_manifest |>
        dplyr::select(outcome, model_index, accepted_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      diagnostic_audit$diagnostic_summary |>
        dplyr::select(outcome, model_index, diagnostic_name = model),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      final_loo_audit$model_summary |>
        dplyr::select(outcome, model_index, loo_name = model_name),
      by = c("outcome", "model_index")
    ) |>
    dplyr::left_join(
      stacking_weights |>
        dplyr::select(outcome, model_index, weight_name = model_name),
      by = c("outcome", "model_index")
    )
  expected_name <- paste0(
    name_alignment$specification_name,
    " [strategy_within]"
  )
  weight_sums <- stacking_weights |>
    dplyr::summarise(
      weight_sum = sum(weight),
      models = dplyr::n(),
      .by = outcome
    )
  prediction_sums <- prediction_manifest |>
    dplyr::summarise(
      draws = sum(draw_count),
      selected_models = dplyr::n(),
      retained_weight = sum(weight),
      .by = outcome
    )
  accepted_stage_counts <- table(accepted_fit_manifest$stage)

  checks <- tibble::tibble(
    check = c(
      "three immutable source files have hashes",
      "the updated strategy audit is complete and confirmed",
      "all four strategy-within datasets contain valid matched contrasts",
      "the combined-data preparation audit remains valid",
      "all five strategy-within registries contain the same 104 model keys",
      "model names agree across specifications fits diagnostics LOO and weights",
      "all accepted fit files exist with 101 initial and 3 remedy fits",
      "all final model diagnostics pass",
      "integrated LOO contains 104 reductions and 65 audited high-k cases",
      "all 65 exact refits pass and compact corrections exist",
      "final corrected LOO contains no unresolved high-k cases",
      "four stacking vectors contain 26 models and sum to one",
      "prediction allocation contains 4000 draws per outcome",
      "all four model-averaged strategy-within reporting audits pass"
    ),
    passed = c(
      nrow(source_manifest) == 3L &&
        all(nchar(source_manifest$sha256) == 64L),
      isTRUE(strategy_audit_review$audit_ready) &&
        strategy_audit_review$pending_confirmations == 0L &&
        all(strategy_audit_review$structural_checks$passed),
      identical(names(data_by_outcome), outcomes) &&
        all(vapply(data_by_outcome, nrow, integer(1L)) > 0L) &&
        all(is.finite(analysis_data$yi)) &&
        all(is.finite(analysis_data$vi)) && all(analysis_data$vi > 0) &&
        setequal(
          unique(analysis_data$protein_strategy),
          c("Addition", "Substitution")
        ),
      all(preparation_audit$checks$passed),
      common_registry && length(specification_keys) == 104L,
      nrow(name_alignment) == 104L &&
        all(stats::complete.cases(name_alignment)) &&
        all(name_alignment$accepted_name == expected_name) &&
        all(name_alignment$diagnostic_name == expected_name) &&
        all(name_alignment$loo_name == expected_name) &&
        all(name_alignment$weight_name == expected_name),
      nrow(accepted_fit_manifest) == 104L &&
        identical(unname(as.integer(accepted_stage_counts["initial"])), 101L) &&
        identical(unname(as.integer(accepted_stage_counts["remedy_1"])), 3L) &&
        all(file.exists(accepted_fit_manifest$fit_path)) &&
        all(accepted_fit_manifest$bytes > 0) &&
        all(nchar(accepted_fit_manifest$sha256) == 64L),
      diagnostic_audit$models == 104L &&
        diagnostic_audit$models_passing_sampling == 104L &&
        isTRUE(diagnostic_audit$all_sampling_diagnostics_pass),
      integrated_loo_audit$models == 104L &&
        integrated_loo_audit$exact_refits_required == 65L &&
        nrow(integrated_loo_audit$high_k_cases) == 65L,
      exact_audit$exact_refits == 65L &&
        isTRUE(exact_audit$all_exact_refit_diagnostics_pass) &&
        nrow(exact_audit$correction_manifest) == 65L &&
        all(file.exists(exact_audit$correction_manifest$correction_path)) &&
        all(exact_audit$correction_manifest$bytes > 0) &&
        all(nchar(exact_audit$correction_manifest$sha256) == 64L),
      isTRUE(final_loo_audit$stacking_ready) &&
        nrow(final_loo_audit$model_summary) == 104L &&
        all(final_loo_audit$outcome_summary$unresolved_high_k_cases == 0L),
      nrow(weight_sums) == 4L && all(weight_sums$models == 26L) &&
        all(abs(weight_sums$weight_sum - 1) < 1e-8),
      nrow(prediction_sums) == 4L &&
        all(prediction_sums$draws == 4000L) &&
        all(prediction_sums$retained_weight > 0) &&
        all(prediction_sums$retained_weight <= 1 + 1e-8) &&
        all(file.exists(prediction_manifest$prediction_path)) &&
        all(prediction_manifest$bytes > 0) &&
        all(nchar(prediction_manifest$sha256) == 64L),
      isTRUE(reporting_audit$all_reporting_checks_pass) &&
        nrow(reporting_audit$summary) == 4L &&
        all(reporting_audit$summary$intervals_valid) &&
        all(reporting_audit$summary$reference_increment_valid) &&
        all(reporting_audit$summary$probabilities_valid) &&
        all(reporting_audit$summary$diagonal_contrasts_valid)
    )
  )
  failed <- is.na(checks$passed) | !checks$passed
  if (any(failed)) {
    stop(
      "Updated strategy-within verification failed: ",
      paste(checks$check[failed], collapse = "; "),
      call. = FALSE
    )
  }

  artifacts <- list(
    verification_checks = checks,
    source_manifest = source_manifest,
    strategy_review_rows = strategy_audit_review$review_rows,
    strategy_review_status = strategy_audit_review$status_summary,
    strategy_structural_checks = strategy_audit_review$structural_checks,
    preparation_checks = preparation_audit$checks,
    strategy_within_data_summary = data_summary,
    analysis_data = analysis_data,
    model_specifications = specification_table,
    model_name_alignment = name_alignment,
    accepted_fit_manifest = accepted_fit_manifest,
    sampling_diagnostic_summary = diagnostic_audit$diagnostic_summary,
    sampling_parameter_failures = diagnostic_audit$parameter_failures,
    initial_integrated_loo_model_summary = integrated_loo_audit$model_summary,
    initial_integrated_loo_outcome_summary =
      integrated_loo_audit$outcome_summary,
    initial_high_pareto_k_cases = integrated_loo_audit$high_k_cases,
    exact_refit_diagnostic_summary = exact_audit$diagnostic_summary,
    exact_correction_manifest = exact_audit$correction_manifest,
    final_loo_manifest = final_loo_audit$final_loo_manifest,
    final_loo_model_summary = final_loo_audit$model_summary,
    final_loo_outcome_summary = final_loo_audit$outcome_summary,
    stacking_weights = stacking_weights,
    prediction_allocation_manifest = prediction_manifest,
    reporting_file_manifest = reporting_audit$summary,
    dose_response_predictions = bind_report("dose_response"),
    dose_response_increments = bind_report("dose_response_increments"),
    exposure_support_by_group = bind_report("exposure_support_by_group"),
    common_exposure_support = bind_report("common_exposure_support"),
    addition_minus_substitution = bind_report("addition_minus_substitution"),
    pairwise_contrasts = bind_report("pairwise_contrasts"),
    pairwise_support = bind_report("pairwise_support"),
    local_equivalence = bind_report("local_equivalence")
  )
  descriptions <- c(
    verification_checks = "Final updated strategy-within checkpoint checks.",
    source_manifest = "Paths, sizes, and SHA-256 hashes for source data.",
    strategy_review_rows = "Confirmed updated-search strategy classifications.",
    strategy_review_status = "Classification-review status summary.",
    strategy_structural_checks = "Strategy audit structural invariants.",
    preparation_checks = "Combined-data preparation invariants.",
    strategy_within_data_summary = "Matched contrasts by outcome and strategy.",
    analysis_data = "Matched within-study contrasts supplied to models.",
    model_specifications = "Final 104-model strategy-within registry.",
    model_name_alignment = "Cross-stage model-name identity audit.",
    accepted_fit_manifest = "Accepted fit paths, stages, sizes, and checksums.",
    sampling_diagnostic_summary = "Final model-level sampling diagnostics.",
    sampling_parameter_failures = "Parameter-level diagnostic failures.",
    initial_integrated_loo_model_summary = "Initial integrated PSIS summaries.",
    initial_integrated_loo_outcome_summary = "Initial LOO audit by outcome.",
    initial_high_pareto_k_cases = "The 65 cases selected for exact refitting.",
    exact_refit_diagnostic_summary = "Diagnostics for 65 accepted exact refits.",
    exact_correction_manifest = "Compact exact-correction paths and checksums.",
    final_loo_manifest = "Final exact-corrected LOO lineage.",
    final_loo_model_summary = "Final LOO summary for all 104 models.",
    final_loo_outcome_summary = "Final high-k resolution by outcome.",
    stacking_weights = "Final exact-corrected stacking weights.",
    prediction_allocation_manifest = "Allocation of 4,000 draws per outcome.",
    reporting_file_manifest = "Reporting file dimensions and checksums.",
    dose_response_predictions = "Model-averaged curves by strategy and RT group.",
    dose_response_increments = "Strategy-specific increments from 0.8 g/kg/day.",
    exposure_support_by_group = "Observed support by strategy and RT group.",
    common_exposure_support = "Support shared by both strategies within RT group.",
    addition_minus_substitution = "Addition-minus-substitution contrasts.",
    pairwise_contrasts = "Strategy-specific contrasts and probabilities.",
    pairwise_support = "Observed support used for pairwise contrasts.",
    local_equivalence = "Strategy-specific practical-equivalence summaries."
  )
  list(
    checks = checks,
    artifacts = artifacts,
    descriptions = descriptions,
    next_step = paste(
      "The updated strategy-within family is frozen after this verification",
      "checkpoint; its lean-pruning plan can now be audited."
    )
  )
}

make_update_2026_strategy_within_pruning_plan <- function(
    accepted_fit_manifest,
    accepted_exact_fit_manifest,
    verification_files,
    archive_root = here::here(
      "model_archive", "update_2026", "sensitivity", "strategy_within"
    )
) {
  verification_manifest <- file.path(
    here::here(
      "results", "verification", "strategy_within_update_2026"
    ),
    "manifest.csv"
  )
  if (
    length(verification_files) != 34L ||
      any(!file.exists(verification_files))
  ) {
    stop(
      "The updated strategy-within verification bundle is incomplete.",
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(verification_manifest)

  resolved_root <- normalizePath(archive_root, winslash = "/", mustWork = TRUE)
  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "strategy_within"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(resolved_root), tolower(expected_root))) {
    stop(
      "Refusing to plan pruning outside the updated strategy-within archive.",
      call. = FALSE
    )
  }

  model_fit_roots <- file.path(archive_root, c("candidates", "remedies"))
  if (any(!dir.exists(model_fit_roots))) {
    stop(
      "An expected updated strategy-within fit directory is missing.",
      call. = FALSE
    )
  }
  model_fit_paths <- unlist(
    purrr::map(
      model_fit_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  ) |>
    normalizePath(winslash = "/", mustWork = TRUE)
  accepted_model_paths <- normalizePath(
    accepted_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (
    length(model_fit_paths) != 107L ||
      length(accepted_model_paths) != 104L ||
      anyDuplicated(model_fit_paths) || anyDuplicated(accepted_model_paths) ||
      !all(accepted_model_paths %in% model_fit_paths) ||
      sum(!model_fit_paths %in% accepted_model_paths) != 3L
  ) {
    stop(
      "Updated strategy-within model retention does not match its lineage.",
      call. = FALSE
    )
  }
  model_fit_manifest <- tibble::tibble(
    path = model_fit_paths,
    relative_path = substring(model_fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/candidates/", model_fit_paths, fixed = TRUE) ~ "initial",
      grepl(
        "/remedies/remedy_1/",
        model_fit_paths,
        fixed = TRUE
      ) ~ "remedy_1",
      TRUE ~ "unknown"
    ),
    outcome = basename(dirname(model_fit_paths)),
    model_index = as.integer(sub(
      "^model_([0-9]+)\\.qs$",
      "\\1",
      basename(model_fit_paths)
    )),
    accepted_fit = model_fit_paths %in% accepted_model_paths,
    bytes = as.numeric(file.info(model_fit_paths)$size),
    lean_action = "delete frozen strategy-within model fit archive"
  ) |>
    dplyr::arrange(archive_stage, outcome, model_index)
  if (
    any(model_fit_manifest$archive_stage == "unknown") ||
      anyNA(model_fit_manifest$model_index) ||
      !identical(
        as.integer(table(factor(
          model_fit_manifest$archive_stage,
          levels = c("initial", "remedy_1")
        ))),
        c(104L, 3L)
      ) ||
      !identical(
        as.integer(table(factor(
          model_fit_manifest$archive_stage[
            model_fit_manifest$accepted_fit
          ],
          levels = c("initial", "remedy_1")
        ))),
        c(101L, 3L)
      ) ||
      any(model_fit_manifest$bytes <= 0)
  ) {
    stop(
      "The updated strategy-within model-fit stages are not as expected.",
      call. = FALSE
    )
  }

  exact_fit_roots <- file.path(
    archive_root,
    c("loo_exact/fits", "loo_exact/remedies")
  )
  if (any(!dir.exists(exact_fit_roots))) {
    stop(
      "An updated strategy-within exact-refit directory is missing.",
      call. = FALSE
    )
  }
  exact_fit_paths <- unlist(
    purrr::map(
      exact_fit_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  ) |>
    normalizePath(winslash = "/", mustWork = TRUE)
  accepted_exact_paths <- normalizePath(
    accepted_exact_fit_manifest$fit_path,
    winslash = "/",
    mustWork = TRUE
  )
  exact_names <- basename(exact_fit_paths)
  exact_refit_manifest <- tibble::tibble(
    path = exact_fit_paths,
    relative_path = substring(exact_fit_paths, nchar(resolved_root) + 2L),
    archive_stage = dplyr::case_when(
      grepl("/loo_exact/fits/", exact_fit_paths, fixed = TRUE) ~
        "observation_exact",
      grepl(
        "/loo_exact/remedies/remedy_1/",
        exact_fit_paths,
        fixed = TRUE
      ) ~ "observation_exact_remedy_1",
      TRUE ~ "unknown"
    ),
    outcome = basename(dirname(exact_fit_paths)),
    model_index = as.integer(sub(
      "^model_([0-9]+)_.*$",
      "\\1",
      exact_names
    )),
    case_id = as.integer(sub(
      "^.*_case_([0-9]+)_.*$",
      "\\1",
      exact_names
    )),
    observation = as.integer(sub(
      "^.*_observation_([0-9]+)\\.qs$",
      "\\1",
      exact_names
    )),
    accepted_fit = exact_fit_paths %in% accepted_exact_paths,
    bytes = as.numeric(file.info(exact_fit_paths)$size),
    lean_action = "delete frozen strategy-within exact-refit archive"
  ) |>
    dplyr::arrange(archive_stage, case_id)
  if (
    length(exact_fit_paths) != 66L ||
      length(accepted_exact_paths) != 65L ||
      anyDuplicated(exact_fit_paths) || anyDuplicated(accepted_exact_paths) ||
      !all(accepted_exact_paths %in% exact_fit_paths) ||
      sum(!exact_fit_paths %in% accepted_exact_paths) != 1L ||
      any(exact_refit_manifest$archive_stage == "unknown") ||
      anyNA(exact_refit_manifest$model_index) ||
      anyNA(exact_refit_manifest$case_id) ||
      anyNA(exact_refit_manifest$observation) ||
      !identical(
        as.integer(table(factor(
          exact_refit_manifest$archive_stage,
          levels = c(
            "observation_exact",
            "observation_exact_remedy_1"
          )
        ))),
        c(65L, 1L)
      ) ||
      sum(exact_refit_manifest$accepted_fit) != 65L ||
      any(exact_refit_manifest$bytes <= 0)
  ) {
    stop(
      "The updated strategy-within exact-refit inventory is incomplete.",
      call. = FALSE
    )
  }

  deletion_manifest <- dplyr::bind_rows(
    model_fit_manifest |>
      dplyr::transmute(
        path,
        relative_path,
        archive_stage,
        artifact_role = dplyr::if_else(
          accepted_fit,
          "accepted_model_fit",
          "superseded_model_fit"
        ),
        outcome,
        model_index,
        case_id = NA_integer_,
        observation = NA_integer_,
        bytes,
        lean_action
      ),
    exact_refit_manifest |>
      dplyr::transmute(
        path,
        relative_path,
        archive_stage,
        artifact_role = dplyr::if_else(
          accepted_fit,
          "accepted_exact_refit",
          "superseded_exact_refit"
        ),
        outcome,
        model_index,
        case_id,
        observation,
        bytes,
        lean_action
      )
  ) |>
    dplyr::arrange(artifact_role, outcome, model_index, case_id)
  role_counts <- table(deletion_manifest$artifact_role)
  if (
    nrow(deletion_manifest) != 173L ||
      unname(role_counts["accepted_model_fit"]) != 104L ||
      unname(role_counts["superseded_model_fit"]) != 3L ||
      unname(role_counts["accepted_exact_refit"]) != 65L ||
      unname(role_counts["superseded_exact_refit"]) != 1L ||
      anyDuplicated(deletion_manifest$path) ||
      any(deletion_manifest$bytes <= 0)
  ) {
    stop(
      "The updated strategy-within deletion inventory is incomplete.",
      call. = FALSE
    )
  }

  compact_roots <- c(
    loo_psis = file.path(archive_root, "loo_psis"),
    exact_corrections = file.path(
      archive_root,
      "loo_exact",
      "corrections"
    ),
    reporting = file.path(archive_root, "reporting")
  )
  if (any(!dir.exists(compact_roots))) {
    stop(
      "An expected compact strategy-within directory is missing.",
      call. = FALSE
    )
  }
  compact_paths <- unlist(
    purrr::map(
      compact_roots,
      ~ list.files(
        .x,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  compact_paths <- compact_paths[
    file.exists(compact_paths) & !file.info(compact_paths)$isdir
  ] |>
    normalizePath(winslash = "/", mustWork = TRUE)
  compact_manifest <- tibble::tibble(
    path = compact_paths,
    relative_path = substring(compact_paths, nchar(resolved_root) + 2L),
    archive_component = dplyr::case_when(
      grepl("/loo_psis/", compact_paths, fixed = TRUE) ~ "loo_psis",
      grepl(
        "/loo_exact/corrections/",
        compact_paths,
        fixed = TRUE
      ) ~ "exact_corrections",
      grepl("/reporting/", compact_paths, fixed = TRUE) ~ "reporting",
      TRUE ~ "unknown"
    ),
    bytes = as.numeric(file.info(compact_paths)$size),
    sha256 = purrr::map_chr(
      compact_paths,
      ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
    ),
    action = "retain compact strategy-within LOO/correction/reporting archive"
  ) |>
    dplyr::arrange(archive_component, relative_path)
  if (
    nrow(compact_manifest) != 212L ||
      any(compact_manifest$archive_component == "unknown") ||
      !identical(
        as.integer(table(factor(
          compact_manifest$archive_component,
          levels = c("loo_psis", "exact_corrections", "reporting")
        ))),
        c(104L, 65L, 43L)
      ) ||
      any(compact_manifest$bytes <= 0) ||
      any(nchar(compact_manifest$sha256) != 64L)
  ) {
    stop(
      "The compact strategy-within archive inventory is incomplete.",
      call. = FALSE
    )
  }

  summary <- tibble::tibble(
    mode = "lean",
    full_fit_archives_deleted = nrow(deletion_manifest),
    accepted_model_fit_archives_deleted = 104L,
    superseded_model_fit_archives_deleted = 3L,
    accepted_exact_refit_archives_deleted = 65L,
    superseded_exact_refit_archives_deleted = 1L,
    full_fit_archives_retained = 0L,
    gib_reclaimed = sum(deletion_manifest$bytes) / 1024^3,
    compact_archive_files_retained = nrow(compact_manifest),
    compact_archive_gib_retained = sum(compact_manifest$bytes) / 1024^3,
    verification_files_retained = length(verification_files),
    restoration = paste(
      "Reported results, corrected LOO, and all 104 weights remain frozen;",
      "restoring model or exact-refit objects requires refitting the",
      "strategy-within family."
    )
  )
  list(
    summary = summary,
    model_fit_manifest = model_fit_manifest,
    exact_refit_manifest = exact_refit_manifest,
    deletion_file_manifest = deletion_manifest,
    compact_archive_manifest = compact_manifest,
    verification_manifest = normalizePath(
      verification_manifest,
      winslash = "/",
      mustWork = TRUE
    ),
    archive_root = resolved_root
  )
}

write_update_2026_strategy_within_pruning_plan <- function(
    plan,
    output_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within", "pruning"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    summary = file.path(output_dir, "summary.csv"),
    model_fit_manifest = file.path(output_dir, "model_fit_manifest.csv"),
    exact_refit_manifest = file.path(output_dir, "exact_refit_manifest.csv"),
    deletion_file_manifest = file.path(
      output_dir,
      "deletion_file_manifest.csv"
    ),
    compact_archive_manifest = file.path(
      output_dir,
      "compact_archive_manifest.csv"
    )
  )
  readr::write_csv(plan$summary, paths[["summary"]])
  readr::write_csv(plan$model_fit_manifest, paths[["model_fit_manifest"]])
  readr::write_csv(plan$exact_refit_manifest, paths[["exact_refit_manifest"]])
  readr::write_csv(
    plan$deletion_file_manifest,
    paths[["deletion_file_manifest"]]
  )
  readr::write_csv(
    plan$compact_archive_manifest,
    paths[["compact_archive_manifest"]]
  )
  unname(paths)
}

execute_update_2026_strategy_within_pruning <- function(
    plan,
    confirmation,
    receipt_dir = here::here(
      "results", "update_2026", "sensitivity", "strategy_within", "pruning"
    )
) {
  required_confirmation <-
    "DELETE ALL UPDATED STRATEGY-WITHIN FIT ARCHIVES"
  if (
    length(confirmation) != 1L ||
      !identical(confirmation, required_confirmation)
  ) {
    stop(
      "Pruning confirmation must be exactly: ",
      required_confirmation,
      call. = FALSE
    )
  }
  verify_update_2026_manifest_files(plan$verification_manifest)

  expected_root <- normalizePath(
    here::here(
      "model_archive", "update_2026", "sensitivity", "strategy_within"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(tolower(plan$archive_root), tolower(expected_root))) {
    stop(
      "The pruning plan does not target the strategy-within archive.",
      call. = FALSE
    )
  }
  selected <- plan$deletion_file_manifest
  role_counts <- table(selected$artifact_role)
  if (
    nrow(plan$model_fit_manifest) != 107L ||
      sum(plan$model_fit_manifest$accepted_fit) != 104L ||
      sum(!plan$model_fit_manifest$accepted_fit) != 3L ||
      nrow(plan$exact_refit_manifest) != 66L ||
      sum(plan$exact_refit_manifest$accepted_fit) != 65L ||
      sum(!plan$exact_refit_manifest$accepted_fit) != 1L ||
      nrow(selected) != 173L ||
      unname(role_counts["accepted_model_fit"]) != 104L ||
      unname(role_counts["superseded_model_fit"]) != 3L ||
      unname(role_counts["accepted_exact_refit"]) != 65L ||
      unname(role_counts["superseded_exact_refit"]) != 1L
  ) {
    stop(
      "The strategy-within pruning set is inconsistent with its audit.",
      call. = FALSE
    )
  }

  receipt_path <- file.path(receipt_dir, "executed_lean.csv")
  pending_receipt_path <- file.path(
    receipt_dir,
    "executed_lean.pending.csv"
  )
  if (file.exists(receipt_path) || file.exists(pending_receipt_path)) {
    stop(
      "An updated strategy-within lean-pruning receipt already exists.",
      call. = FALSE
    )
  }

  selected_paths <- normalizePath(
    selected$path,
    winslash = "/",
    mustWork = TRUE
  )
  deletion_roots <- file.path(
    plan$archive_root,
    c(
      "candidates",
      "remedies",
      "loo_exact/fits",
      "loo_exact/remedies"
    )
  )
  if (any(!dir.exists(deletion_roots))) {
    stop(
      "An audited strategy-within deletion subtree is missing.",
      call. = FALSE
    )
  }
  current_paths <- unlist(
    purrr::map(
      deletion_roots,
      ~ list.files(
        .x,
        pattern = "\\.qs$",
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE
      )
    ),
    use.names = FALSE
  )
  current_paths <- current_paths[
    file.exists(current_paths) & !file.info(current_paths)$isdir
  ] |>
    normalizePath(winslash = "/", mustWork = TRUE)
  permitted_roots <- normalizePath(
    deletion_roots,
    winslash = "/",
    mustWork = TRUE
  )
  permitted <- vapply(selected_paths, function(path) {
    any(vapply(permitted_roots, function(root) {
      startsWith(tolower(path), paste0(tolower(root), "/"))
    }, logical(1L)))
  }, logical(1L))
  root_prefix <- paste0(tolower(expected_root), "/")
  if (
    !setequal(selected_paths, current_paths) ||
      any(!startsWith(tolower(selected_paths), root_prefix)) ||
      any(!permitted) || any(dirname(selected_paths) == expected_root)
  ) {
    stop(
      "A selected file falls outside the audited strategy-within fit inventory.",
      call. = FALSE
    )
  }

  retained <- plan$compact_archive_manifest
  retained_paths <- retained$path
  message("Revalidating 212 compact strategy-within artifacts...")
  if (
    nrow(retained) != 212L ||
      any(!file.exists(retained_paths))
  ) {
    stop(
      "A retained compact strategy-within artifact failed revalidation.",
      call. = FALSE
    )
  }
  retained_hashes <- purrr::map_chr(
    retained_paths,
    ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
  )
  if (any(retained_hashes != retained$sha256)) {
    stop(
      "A retained compact strategy-within artifact failed revalidation.",
      call. = FALSE
    )
  }

  current_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(current_commit) != 1L || !nzchar(current_commit)) {
    current_commit <- NA_character_
  }
  message("Hashing 173 strategy-within fit archives before deletion...")
  receipt <- selected |>
    dplyr::mutate(
      sha256_before_deletion = purrr::map_chr(
        selected_paths,
        ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
      ),
      pruning_mode = "lean",
      deleted_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = current_commit,
      .before = 1L
    )
  dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(receipt, pending_receipt_path)
  pending <- readr::read_csv(
    pending_receipt_path,
    show_col_types = FALSE
  )
  if (
    nrow(pending) != 173L ||
      any(nchar(pending$sha256_before_deletion) != 64L)
  ) {
    stop(
      "The pending strategy-within pruning receipt is incomplete.",
      call. = FALSE
    )
  }

  message("Deleting the 173 audited strategy-within fit archives...")
  deleted <- file.remove(selected_paths)
  if (!all(deleted) || any(file.exists(selected_paths))) {
    stop(
      "At least one strategy-within archive could not be deleted; ",
      "the pending receipt records the attempted inventory.",
      call. = FALSE
    )
  }
  final_retained_hashes <- purrr::map_chr(
    retained_paths,
    ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
  )
  if (
    any(!file.exists(retained_paths)) ||
      any(final_retained_hashes != retained$sha256) ||
      !file.exists(plan$verification_manifest)
  ) {
    stop(
      "A retained strategy-within artifact is unexpectedly missing or changed.",
      call. = FALSE
    )
  }
  if (!file.rename(pending_receipt_path, receipt_path)) {
    stop(
      "Could not finalise the strategy-within lean-pruning receipt.",
      call. = FALSE
    )
  }
  message("Strategy-within lean pruning completed; compact artifacts remain.")
  normalizePath(receipt_path, winslash = "/", mustWork = TRUE)
}
