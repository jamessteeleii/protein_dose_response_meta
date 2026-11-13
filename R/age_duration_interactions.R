# Supplementary expansion of the frozen primary families only. This module is
# sourced by a separate targets script; it never dispatches the main pipeline.

# Workflow families and archive paths ----

# The age and duration analysis starts from accepted primary candidates, adds
# interaction terms where specified, and reports dose contrasts at supported
# moderator profiles.
iad_families <- function() c("relative", "total", "supplement")

iad_outcomes <- function() c("whole_body", "hypertrophy", "strength", "function")

iad_phase <- function(phase) {
  value <- Sys.getenv("PROTEIN_INTERACTION_PHASE", "none")
  if (!value %in% c("none", "pilot", "fit", "batch_fit", "remedy", "loo", "exact", "exact_remedy", "correction", "combine", "report",
      "recover_pilot", "recover_fit", "recover_remedy", "recover_loo", "recover_exact", "recover_correction", "recover_accept", "recover_report")) {
    stop("Unknown PROTEIN_INTERACTION_PHASE.")
  }
  identical(value, phase) || (identical(value, "batch_fit") && phase %in% c("pilot", "fit"))
}

iad_path <- function(..., archive = FALSE) {
  here::here(if (archive) "model_archive" else "results", "update_2026", "age_duration", ...)
}

iad_qsave <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  qs::qsave(object, temporary, preset = "high")
  if (!file.rename(temporary, path)) stop("Could not atomically finalise interaction archive.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

iad_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  if (grepl("\\.gz$", path)) write_verification_csv(object, temporary) else
    readr::write_csv(object, temporary, na = "")
  if (!file.rename(temporary, path)) stop("Could not atomically write interaction table.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

iad_read_csv <- function(path) {
  data <- readr::read_csv(path, guess_max = Inf, show_col_types = FALSE)
  if (nrow(readr::problems(data))) stop("Lossy checkpoint parsing: ", basename(path))
  data
}

iad_hash <- function(path) digest::digest(file = path, algo = "sha256")

iad_check_hash <- function(path, hash) {
  if (!file.exists(path) || length(hash) != 1L || is.na(hash) ||
      !identical(iad_hash(path), unname(as.character(hash)))) stop("Checksum mismatch: ", path)
  invisible(TRUE)
}

iad_bundle_table <- function(directory, manifest, artifact) {
  row <- manifest[manifest$artifact == artifact, , drop = FALSE]
  if (nrow(row) != 1L || basename(row$file) != row$file) stop("Invalid source artifact: ", artifact)
  path <- file.path(directory, row$file)
  iad_check_hash(path, row$sha256)
  iad_read_csv(path)
}

# Verify frozen parent analysis inputs ----

# Identify the frozen tables and fit artifacts used to reproduce each parent
# ensemble, and match the original observation order before reusing its LOO
# scores.
iad_checkpoint_inputs <- function() {
  folders <- c(relative = "primary_update_2026", total = "total_primary_update_2026",
               supplement = "supplement_primary_update_2026")
  unlist(lapply(folders, function(folder) {
    root <- here::here("results", "verification", folder)
    file.path(root, c("manifest.csv", "analysis_data.csv.gz", "model_specifications.csv.gz",
      "verification_checks.csv.gz", if (folder == "primary_update_2026") "loo_final_manifest.csv.gz" else "final_loo_manifest.csv.gz",
      if (folder == "primary_update_2026") "accepted_fit_manifest.csv.gz" else "final_fit_manifest.csv.gz",
      if (folder == "primary_update_2026") "diagnostic_summary.csv.gz" else "final_sampling_diagnostic_summary.csv.gz"))
  }), use.names = FALSE)
}

iad_identity_columns <- function() c("study_no", "arm_no", "outcome_no", "time1", "time2", "yi", "vi")

iad_same_observations <- function(first, second) {
  # Older brms spline fits keep only formula variables. Their study_factor
  # labels carry the same study IDs; factor integer codes are not those IDs.
  if (!"study_no" %in% names(first) && "study_factor" %in% names(first)) {
    first$study_no <- as.character(first$study_factor)
  }

  if (!"study_no" %in% names(second) && "study_factor" %in% names(second)) {
    second$study_no <- as.character(second$study_factor)
  }
  columns <- c("study_no", "arm_no", "outcome_no", "yi", "vi")
  if (!all(columns %in% names(first)) || !all(columns %in% names(second))) return(FALSE)
  columns <- c(columns, intersect(c("time1", "time2"), intersect(names(first), names(second))))
  # Identifiers may be integer in fits and double/character in CSV checkpoints.
  normalise <- function(x) {
    lapply(columns, function(column) {
      if (column %in% c("study_no", "arm_no", "outcome_no")) as.character(x[[column]]) else
        as.numeric(x[[column]])
    })
  }
  isTRUE(all.equal(normalise(first), normalise(second), tolerance = 1e-13, check.attributes = FALSE))
}

iad_fingerprint <- function(data) {
  columns <- c(iad_identity_columns(), "protein_intake", "protein_intake_raw",
               "age_centre", "duration_centre", "resistance_exercise_code")
  if ("mean_body_mass_kg" %in% names(data)) columns <- c(columns, "mean_body_mass_kg")
  value <- lapply(columns, function(column) {
    if (column %in% c("study_no", "arm_no", "outcome_no")) as.character(data[[column]]) else
      signif(as.numeric(data[[column]]), 13L)
  })
  digest::digest(value, algo = "sha256")
}

iad_validate_parent_data <- function(candidate, reference, family) {
  data <- candidate$model$data
  expected <- add_primary_half_unit_protein(add_primary_standardised_moderators(reference))
  if (family != "relative") expected <- add_absolute_primary_body_mass(expected)
  expected$Intercept <- 1
  if (!iad_same_observations(data, reference)) stop("Parent fit observations differ from the frozen common sample.")
  columns <- intersect(names(data), names(expected))
  for (column in columns) {
    first <- data[[column]]; second <- expected[[column]]
    if (is.factor(first) || is.factor(second)) {
      equal <- identical(as.character(first), as.character(second))
    } else equal <- isTRUE(all.equal(first, second, tolerance = 1e-12, check.attributes = FALSE))
    if (!equal) stop("Parent fit predictor differs from checkpoint: ", column)
  }
  invisible(TRUE)
}

iad_validate_loo <- function(result, expected = NULL) {
  if (!inherits(result$loo, "loo") || is.null(result$observation_manifest) ||
      nrow(result$loo$pointwise) != nrow(result$observation_manifest) ||
      any(!is.finite(result$loo$pointwise[, "elpd_loo"]))) stop("Invalid pointwise LOO reduction.")
  if (!is.null(expected) && !iad_same_observations(result$observation_manifest, expected)) {
    stop("LOO observations, ordering or likelihood inputs differ; combination is blocked.")
  }
  k <- loo::pareto_k_values(result$loo)
  exact <- result$exact_refit_summary
  observations <- if (is.data.frame(exact) && "observation" %in% names(exact)) exact$observation else integer()
  if (any(k > 0.7, na.rm = TRUE) || !setequal(which(is.na(k)), observations) ||
      anyDuplicated(observations) || (length(observations) &&
      (!all(exact$diagnostic_pass) || any(!is.finite(exact$exact_elpd_loo))))) {
    stop("Unresolved or unaudited LOO cases; combination is blocked.")
  }
  invisible(TRUE)
}

# Read and validate the parent manifests, model data and fit paths needed to
# prepare the interaction candidates.
iad_read_checkpoints <- function(paths) {
  if (!all(file.exists(paths))) stop("Frozen primary checkpoint inputs are missing.")
  folders <- c(relative = "primary_update_2026", total = "total_primary_update_2026",
               supplement = "supplement_primary_update_2026")
  output <- lapply(iad_families(), function(family) {
    directory <- here::here("results", "verification", folders[[family]])
    manifest <- iad_read_csv(file.path(directory, "manifest.csv"))
    get <- function(name) iad_bundle_table(directory, manifest, name)
    checks <- get("verification_checks")
    if (!nrow(checks) || !all(checks$passed)) stop("Parent checkpoint checks failed: ", family)
    data <- get("analysis_data")
    specifications <- get("model_specifications")
    fits <- get(if (family == "relative") "accepted_fit_manifest" else "final_fit_manifest")
    loo_manifest <- get(if (family == "relative") "loo_final_manifest" else "final_loo_manifest")
    diagnostics <- get(if (family == "relative") "diagnostic_summary" else "final_sampling_diagnostic_summary")
    if (!all(diagnostics$diagnostic_pass)) stop("Parent sampling diagnostics failed.")
    records <- list(); datasets <- list(); source_rows <- list()
    for (outcome in iad_outcomes()) {
      s <- specifications[specifications$outcome == outcome, , drop = FALSE]
      s <- s[order(s$model_index), , drop = FALSE]
      f <- fits[fits$outcome == outcome, , drop = FALSE]
      l <- loo_manifest[loo_manifest$outcome == outcome, , drop = FALSE]
      l <- l[order(l$model_index), , drop = FALSE]
      n <- if (family == "relative") 26L else 18L
      if (nrow(s) != n || nrow(f) != n || nrow(l) != n ||
          anyDuplicated(s$model_index) || !setequal(s$model_index, f$model_index) ||
          !identical(s$model_index, l$model_index)) stop("Incomplete frozen primary registry.")
      reference <- NULL
      for (i in seq_len(nrow(l))) {
        iad_check_hash(l$final_loo_path[[i]], l$final_sha256[[i]])
        reduction <- qs::qread(l$final_loo_path[[i]])
        iad_validate_loo(reduction, reference)
        if (reduction$model_index != l$model_index[[i]] || reduction$outcome != outcome ||
            reduction$model_name != s$model_name[[i]]) stop("Parent LOO identity mismatch.")
        if (is.null(reference)) reference <- reduction$observation_manifest
      }
      source <- data[data$analysis_outcome == outcome, , drop = FALSE]
      source_row_count <- nrow(source)
      # Historical fits and LOO manifests omit some source fields. Restore the
      # documented common sample, preserve its order and check the likelihood
      # row by row. When available, validate every saved predictor in one
      # accepted moderator fit, including its standardised duration.
      if (family != "relative") source <- source[is.finite(source$mean_body_mass_kg), ]
      if (!iad_same_observations(reference, source)) {
        stop("The frozen source ordering/common sample does not match final LOO: ", family, "/", outcome)
      }
      retained <- source
      anchor_path <- f$fit_path[f$model_index == 6L]
      if (length(anchor_path) != 1L) stop("Missing parent row-order anchor record.")
      if (file.exists(anchor_path)) {
        anchor <- qs::qread(anchor_path)
        if (anchor$model_index != 6L || anchor$outcome != outcome) stop("Parent fit identity differs.")
        iad_validate_parent_data(anchor, retained, family)
        rm(anchor); gc(verbose = FALSE)
      }
      retained$resistance_exercise_factor <- factor(retained$resistance_exercise_factor, c("yes", "no"))
      retained$study_factor <- factor(as.character(retained$study_no))
      if (anyNA(retained$resistance_exercise_factor) || !setequal(retained$resistance_exercise_code, c(0, 1)) ||
          any(retained$resistance_exercise_code != ifelse(retained$resistance_exercise_factor == "yes", 0, 1))) {
        stop("Parent resistance-training coding is invalid.")
      }
      f <- f[match(s$model_index, f$model_index), , drop = FALSE]
      if (any(f$fit_path != l$fit_path)) stop("Parent fit and LOO manifests disagree.")
      records[[outcome]] <- tibble::tibble(family = family, outcome = outcome,
        model_index = as.integer(s$model_index), parent_index = as.integer(s$model_index),
        model_name = s$model_name, interaction = "none", fit_path = f$fit_path,
        stage = f$stage, loo_path = l$final_loo_path, loo_sha256 = l$final_sha256,
        available = file.exists(f$fit_path),
        expected_bytes = if ("bytes" %in% names(f)) f$bytes else file.info(f$fit_path)$size)
      datasets[[outcome]] <- retained
      source_rows[[outcome]] <- tibble::tibble(family = family, outcome = outcome,
        source_rows = source_row_count, common_rows = nrow(retained),
        data_sha256 = iad_fingerprint(retained), studies = dplyr::n_distinct(retained$study_no),
        arms = dplyr::n_distinct(retained$arm_no), required_parent_fits = nrow(f),
        available_parent_fits = sum(file.exists(f$fit_path)))
    }
    list(family = family, data = datasets, specifications = specifications, diagnostics = diagnostics,
      parents = dplyr::bind_rows(records), counts = dplyr::bind_rows(source_rows),
      source_manifest = manifest[manifest$artifact %in% c("analysis_data", "model_specifications",
        "verification_checks", "final_loo_manifest", "loo_final_manifest", "accepted_fit_manifest", "final_fit_manifest",
        "diagnostic_summary", "final_sampling_diagnostic_summary"), ])
  })

  names(output) <- iad_families()
  output
}

# There is no target dependency on measurement fits: only their completed,
# checksummed public reporting checkpoints can release these expensive stages.
iad_require_measurement <- function() {
  for (variant in c("prereg_methods", "direct_only")) {
    root <- here::here("results", "verification", paste0("measurement_", variant, "_update_2026"))
    path <- file.path(root, "manifest.csv")
    if (!file.exists(path)) stop("Wait for both muscle measurement analyses to finish: ", variant)
    manifest <- iad_read_csv(path)
    for (i in seq_len(nrow(manifest))) {
      if (basename(manifest$file[[i]]) != manifest$file[[i]]) stop("Invalid measurement manifest path.")
      iad_check_hash(file.path(root, manifest$file[[i]]), manifest$sha256[[i]])
    }
    checks <- iad_bundle_table(root, manifest, "verification_checks")
    if (!nrow(checks) || !all(checks$passed)) stop("Measurement verification has not passed: ", variant)
  }
  invisible(TRUE)
}

iad_require_phase <- function(phase) {
  if (!iad_phase(phase)) stop("Interaction phase disabled: ", phase)
  iad_require_measurement()
}

iad_require_parents <- function(checkpoints) {
  parents <- dplyr::bind_rows(lapply(checkpoints, `[[`, "parents"))
  if (any(!file.exists(parents$fit_path))) stop(
    "Restore the missing accepted primary fit archives before fitting/combining this extension; see preparation/counts.csv.")
  invisible(TRUE)
}

# Interaction formulas, priors and profiles ----

# Extend a parent dose-response formula for age or duration moderation while
# using the parent model registry to retain candidate identity.
iad_formula <- function(parent, type) {
  if (!type %in% c("age", "duration", "both")) stop("Invalid interaction type.")
  formula <- replace_primary_standardised_moderators(parent$formula)
  moderators <- switch(type, age = "age_outcome_z", duration = "duration_outcome_z",
                       both = c("age_outcome_z", "duration_outcome_z"))
  index <- as.integer(parent$model_index)
  rt <- index %in% c(18L, 19L, 20L, 24L, 25L, 26L)
  smooth <- index %in% c(8L, 14L, 20L, 26L)
  spline <- index %in% c(7L, 13L, 19L, 25L)
  x <- if (grepl("protein_intake_05", formula, fixed = TRUE)) "protein_intake_05" else "protein_intake"
  # Natural-spline remedies are extended using their same fixed basis.
  natural <- grepl("splines::ns", formula, fixed = TRUE)
  basis <- if (natural) "splines::ns(protein_intake, df = 3)" else
    if (spline) "lspline(protein_intake, 0.5)" else x
  terms <- character()
  for (moderator in moderators) {
    if (rt) terms <- c(terms, paste0(moderator, ":resistance_exercise_code"))
    if (smooth && !natural) {
      by <- if (rt) paste0(moderator, c("_rt_yes", "_rt_no")) else moderator
      # A point constraint prevents the numeric-by smooth from duplicating the
      # moderator's additive main effect. Zero is the frozen centred dose.
      terms <- c(terms, sprintf("s(%s, by=%s, bs='tp', k=4, pc=0)", x, by))
    } else {
      terms <- c(terms, paste0(basis, ":", moderator))
      if (rt) terms <- c(terms, paste0(basis, ":", moderator, ":resistance_exercise_code"))
    }
  }

  paste(formula, paste(terms, collapse = " + "), sep = " + ")
}

iad_specifications <- function(checkpoints) {
  records <- list()
  for (family in iad_families()) for (outcome in iad_outcomes()) {
    table <- checkpoints[[family]]$specifications
    parents <- table[table$outcome == outcome & table$model_index %in%
      c(6L, 7L, 8L, 12L, 13L, 14L, 18L, 19L, 20L, 24L, 25L, 26L), , drop = FALSE]
    parents <- parents[order(parents$model_index), , drop = FALSE]
    if (nrow(parents) != if (family == "relative") 12L else 8L) stop("Interaction parent registry is incomplete.")
    next_index <- 26L
    for (type in c("age", "duration", "both")) for (i in seq_len(nrow(parents))) {
      parent <- as.list(parents[i, , drop = FALSE])
      next_index <- next_index + 1L
      s <- parent
      s$family <- family; s$parent_index <- as.integer(parent$model_index)
      s$model_index <- next_index; s$interaction <- type
      s$model_name <- paste0(parent$model_name, " [protein x ", type, "]")
      s$formula <- iad_formula(parent, type)
      s$stage <- "initial"
      s$chains <- 4L; s$parallel_chains <- 4L
      s$iter_warmup <- 2000L; s$iter_sampling <- 6000L
      s$adapt_delta <- max(as.numeric(s$adapt_delta), 0.999)
      s$max_treedepth <- max(as.integer(s$max_treedepth), 15L)
      s$parent_formula <- parent$formula
      # Remove historical paths: these new fits can only write into their own archive.
      s$accepted_fit_path <- s$accepted_fit_stage <- NULL
      records[[length(records) + 1L]] <- s
    }
  }

  names(records) <- vapply(records, function(s) paste(s$family, s$outcome, s$model_index, sep = "_"), character(1))
  records
}

iad_augment_data <- function(data) {
  for (moderator in c("age_outcome_z", "duration_outcome_z")) {
    data[[paste0(moderator, "_rt_yes")]] <- data[[moderator]] * (data$resistance_exercise_code == 0)
    data[[paste0(moderator, "_rt_no")]] <- data[[moderator]] * (data$resistance_exercise_code == 1)
  }
  data
}

iad_priors <- function(path) {
  data <- iad_read_csv(path)
  data$coef[is.na(data$coef)] <- ""
  if (nrow(data) != 12L || anyDuplicated(data[c("outcome", "coef")]) ||
      any(data$class != "b") || !setequal(data$outcome, iad_outcomes())) stop("Invalid frozen coefficient priors.")
  result <- lapply(iad_outcomes(), function(outcome) {
    selected <- data[data$outcome == outcome, ]
    if (!setequal(selected$coef, c("", "Intercept", "resistance_exercise_code"))) stop("Missing base prior coefficient.")
    priors <- lapply(seq_len(nrow(selected)), function(i) brms::set_prior(selected$prior[[i]],
      class = "b", coef = selected$coef[[i]]))
    do.call(c, priors)
  })
  names(result) <- iad_outcomes()
  result
}

iad_prepare_input <- function(specification, checkpoints, priors) {
  s <- specification
  checkpoint <- checkpoints[[s$family]]
  table <- checkpoint$specifications
  parent <- as.list(table[table$outcome == s$outcome & table$model_index == s$parent_index, , drop = FALSE])
  prepared <- if (s$family == "relative") prepare_update_2026_primary_candidate_input(
    parent, checkpoint$data, priors) else prepare_update_2026_total_primary_input(parent, checkpoint$data, priors)

  if (!iad_same_observations(prepared$data, checkpoint$data[[s$outcome]])) stop("Interaction preparation dropped rows.")
  # Reuse all parent priors, including diagnostic remedies, before adding new
  # parameters. Standardisation is identical to the frozen parent's data.
  prepared$data <- add_primary_half_unit_protein(add_primary_standardised_moderators(prepared$data))
  prepared$data <- iad_augment_data(prepared$data)
  base_prior <- brms::get_prior(stats::as.formula(prepared$formula), data = prepared$data)
  expanded_prior <- brms::get_prior(stats::as.formula(s$formula), data = prepared$data)
  added <- expanded_prior$coef[expanded_prior$class == "b" & nzchar(expanded_prior$coef) &
    !expanded_prior$coef %in% base_prior$coef[base_prior$class == "b"]]

  for (coefficient in unique(added)) prepared$priors <- c(prepared$priors,
    brms::set_prior("student_t(3, 0, 0.5)", class = "b", coef = coefficient))
  # Numeric-by smooth deviations use the parent's sds prior, or an explicitly
  # regularising half-normal when the parent has no smooth-SD prior.
  if (any(expanded_prior$class == "sds") && !any(prepared$priors$class == "sds")) {
    prepared$priors <- c(prepared$priors, brms::set_prior("normal(0, 0.5)", class = "sds"))
  }
  prepared$formula <- s$formula
  prepared$specification <- s
  prepared
}

# Check the prepared model data and specifications before the opt-in fits. The
# profile table records where moderator comparisons can be reported.
iad_preflight <- function(specifications, checkpoints, priors) {
  purrr::map_dfr(specifications, function(s) {
    messages <- character()
    prepared <- withCallingHandlers(iad_prepare_input(s, checkpoints, priors), warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    stan_data <- withCallingHandlers(brms::make_standata(stats::as.formula(prepared$formula),
      data = prepared$data, prior = prepared$priors), warning = function(w) {
        messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    valid <- stan_data$N == nrow(checkpoints[[s$family]]$data[[s$outcome]]) &&
      all(is.finite(stan_data$X)) && qr(stan_data$X)$rank == ncol(stan_data$X)
    if (!valid) stop("Invalid or rank-deficient interaction input: ", s$family, "/", s$outcome, "/", s$model_index)
    tibble::tibble(family = s$family, outcome = s$outcome, model_index = s$model_index,
      parent_index = s$parent_index, interaction = s$interaction, observations = stan_data$N,
      fixed_columns = ncol(stan_data$X), valid = valid,
      data_sha256 = iad_fingerprint(prepared$data), warnings = paste(unique(messages), collapse = " | "))
  })
}

iad_profiles <- function(checkpoints) {
  purrr::map_dfr(iad_families(), function(family) purrr::map_dfr(iad_outcomes(), function(outcome) {
    data <- checkpoints[[family]]$data[[outcome]]
    ages <- dplyr::distinct(data, study_no, arm_no, age_centre)$age_centre + 25
    durations <- dplyr::distinct(data, study_no, time2)$time2
    age <- as.numeric(stats::quantile(ages, c(0.25, 0.5, 0.75)))
    weeks <- as.numeric(stats::quantile(durations, c(0.25, 0.5, 0.75)))
    tibble::tibble(family = family, outcome = outcome,
      profile = c("reference", "typical", "age_q25", "age_q75", "duration_q25", "duration_q75"),
      age = c(25, age[[2]], age[[1]], age[[3]], age[[2]], age[[2]]),
      weeks = c(12, weeks[[2]], weeks[[2]], weeks[[2]], weeks[[1]], weeks[[3]]),
      age_neighbourhood_years = 10, duration_neighbourhood_weeks = 4)
  }))
}

iad_write_preparation <- function(checkpoints, specifications, preflight, profiles) {
  if (nrow(preflight) != length(specifications) || !all(preflight$valid)) stop("Incomplete input audit.")
  tables <- list(counts = dplyr::bind_rows(lapply(checkpoints, `[[`, "counts")),
    specifications = dplyr::bind_rows(specifications), preflight = preflight, profiles = profiles,
    source_manifest = dplyr::bind_rows(lapply(checkpoints, `[[`, "source_manifest")),
    protected_parent_files = dplyr::bind_rows(lapply(checkpoints, `[[`, "parents")),
    storage_estimate = iad_storage_estimate(checkpoints, specifications))
  purrr::imap_chr(tables, function(table, name) iad_csv(table,
    iad_path("preparation", paste0(name, ".csv")))) |> unname()
}

iad_storage_estimate <- function(checkpoints, specifications) {
  table <- purrr::map_dfr(specifications, function(s) {
    parents <- checkpoints[[s$family]]$parents
    path <- parents$fit_path[parents$outcome == s$outcome & parents$model_index == s$parent_index]
    tibble::tibble(family = s$family, outcome = s$outcome,
      interaction = s$interaction, parent_bytes = parents$expected_bytes[parents$outcome == s$outcome & parents$model_index == s$parent_index])
  })
  table |> dplyr::group_by(family) |> dplyr::summarise(new_candidates = dplyr::n(),
    same_size_as_parent_gib = sum(parent_bytes) / 1024^3,
    planning_allowance_gib = 1.5 * sum(parent_bytes) / 1024^3,
    note = "Estimate only; remedies, exact refits and reporting require additional space.", .groups = "drop")
}

iad_pilot_specs <- function(specifications) {
  # Linear and most complex study/RT-smooth joint-interaction candidates, per scope.
  Filter(function(s) s$interaction == "both" && s$parent_index %in% c(6L, 26L), specifications)
}

# Fit, diagnose and accept interaction candidates ----

# Fit an individual interaction candidate. The pilot, remedy and acceptance
# helpers below govern which saved fits can enter the combined ensemble.
iad_fit <- function(specification, checkpoints, priors, preflight, phase, pilot = NULL) {
  iad_require_phase(phase)
  iad_require_parents(checkpoints)
  if (!all(preflight$valid) || nrow(preflight) != length(iad_specifications(checkpoints))) stop("Input gate failed.")
  if (phase == "fit") iad_check_family_pilot(pilot, specification$family, checkpoints)
  input <- iad_prepare_input(specification, checkpoints, priors)
  s <- specification
  path <- iad_path(s$family, s$outcome, s$stage, sprintf("model_%02d.qs", s$model_index), archive = TRUE)

  if (file.exists(path)) {
    previous <- qs::qread(path)
    fields <- c("family", "outcome", "model_index", "formula", "stage", "adapt_delta",
      "max_treedepth", "iter_warmup", "iter_sampling")
    if (!identical(previous$specification[fields], s[fields]) ||
        iad_fingerprint(previous$model$data) != iad_fingerprint(input$data)) {
      stop("Existing fit has a different specification; preserve it and use a reviewed new stage.")
    }
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  started <- proc.time()[["elapsed"]]
  model <- brms::brm(stats::as.formula(input$formula), data = input$data, prior = input$priors,
    chains = 4L, cores = 4L, warmup = s$iter_warmup, iter = s$iter_warmup + s$iter_sampling,
    seed = 731988L + 10000L * match(s$family, iad_families()) +
      100L * match(s$outcome, iad_outcomes()) + s$model_index,
    control = list(adapt_delta = s$adapt_delta, max_treedepth = s$max_treedepth),
    save_pars = brms::save_pars(all = TRUE))

  for (column in setdiff(names(input$data), names(model$data))) model$data[[column]] <- input$data[[column]]
  candidate <- list(outcome = s$outcome, model_index = s$model_index, model_name = s$model_name,
    specification = s, model = model, elapsed_seconds = proc.time()[["elapsed"]] - started,
    data_sha256 = iad_fingerprint(input$data))
  iad_qsave(candidate, path)
}

iad_diagnose <- function(path) {
  candidate <- qs::qread(path)
  s <- candidate$specification
  result <- get_primary_candidate_model_diagnostics(candidate, ppc_draws = 400L)
  result$summary$family <- s$family
  result$summary$model_name <- s$model_name
  result$summary$parent_index <- s$parent_index
  result$summary$interaction <- s$interaction
  result$summary$stage <- s$stage
  result$summary$elapsed_seconds <- candidate$elapsed_seconds
  result$summary$fit_path <- path
  result$summary$sha256 <- iad_hash(path)
  result$summary$data_sha256 <- iad_fingerprint(candidate$model$data)
  result$parameter_failures$family <- rep(s$family, nrow(result$parameter_failures))
  result$parameter_failures$outcome <- rep(s$outcome, nrow(result$parameter_failures))
  result$parameter_failures$model_index <- rep(s$model_index, nrow(result$parameter_failures))
  result$standardised_residuals$family <- s$family

  if (!identical(result$summary$data_sha256, candidate$data_sha256)) stop("Archived interaction data changed.")
  result
}

iad_diagnostic_files <- function(diagnostics, stage, family = NULL) {
  summary <- purrr::map_dfr(diagnostics, "summary")
  summary$ppc_warning <- primary_ppc_warning(summary)
  tables <- list(summary = summary, parameter_failures = purrr::map_dfr(diagnostics, "parameter_failures"),
    posterior_predictive_observations = purrr::map_dfr(diagnostics, "standardised_residuals"))
  root <- if (is.null(family)) iad_path("diagnostics", stage) else iad_path(family, "diagnostics", stage)
  unname(purrr::imap_chr(tables, function(table, name) iad_csv(table,
    file.path(root, paste0(name, ".csv")))))
}

iad_diagnostic_summary <- function(diagnostics) purrr::map_dfr(diagnostics, "summary")

iad_remedy_specs <- function(specifications, path, diagnostics, family = NULL) {
  decisions <- iad_read_csv(path)
  if (!is.null(family)) decisions <- decisions[decisions$family == family, , drop = FALSE]
  records <- list()
  if (!nrow(decisions)) return(list(list(skip = TRUE)))
  if (anyNA(decisions) || any(!grepl("^remedy_[1-9][0-9]*$", decisions$stage)) ||
      anyDuplicated(decisions[c("family", "outcome", "model_index", "stage")]) ||
      any(!nzchar(trimws(decisions$rationale)))) stop("Invalid interaction remedy ledger.")

  for (i in seq_len(nrow(decisions))) {
    row <- decisions[i, ]
    found <- Filter(function(s) s$family == row$family && s$outcome == row$outcome &&
      s$model_index == row$model_index, specifications)
    if (length(found) != 1L) stop("Unknown remedy candidate.")
    d <- diagnostics[diagnostics$family == row$family & diagnostics$outcome == row$outcome &
      diagnostics$model_index == row$model_index, ]
    if (nrow(d) != 1L || isTRUE(d$diagnostic_pass)) stop("Remedy must address a failed initial diagnostic.")
    s <- found[[1L]]
    if (row$adapt_delta < s$adapt_delta || row$adapt_delta >= 1 ||
        row$max_treedepth < s$max_treedepth || row$max_treedepth %% 1 != 0) stop("Invalid sampler escalation.")
    s$adapt_delta <- row$adapt_delta; s$max_treedepth <- as.integer(row$max_treedepth)
    if ("iter_sampling" %in% names(row)) {
      sampling <- row$iter_sampling[[1L]]
      if (!is.numeric(sampling) || !is.finite(sampling) || sampling %% 1 != 0 ||
          sampling < s$iter_sampling || sampling > .Machine$integer.max - s$iter_warmup) {
        stop("Invalid retained sampling length in remedy ledger.")
      }
      s$iter_sampling <- as.integer(sampling)
    }
    s$stage <- row$stage; s$remedy_rationale <- row$rationale
    records[[i]] <- s
  }

  records
}

iad_fit_remedy <- function(record, checkpoints, priors, preflight, family = NULL) {
  iad_require_phase("remedy")
  if (isTRUE(record$skip)) return(iad_qsave(record,
    if (is.null(family)) iad_path("remedies", "none.qs", archive = TRUE) else
      iad_path(family, "remedies", "none.qs", archive = TRUE)))
  iad_fit(record, checkpoints, priors, preflight, "remedy")
}

iad_diagnose_remedy <- function(path) {
  if (basename(path) == "none.qs") return(list(summary = tibble::tibble(), parameter_failures = tibble::tibble()))
  iad_diagnose(path)
}

iad_accepted <- function(initial, remedies, specifications, checkpoints) {
  rows <- dplyr::bind_rows(initial, remedies)
  rows$rank <- suppressWarnings(as.integer(sub("remedy_", "", rows$stage)))
  rows$rank[rows$stage == "initial"] <- 0L
  rows$rank[rows$stage == "recovery_1"] <- 0L
  rows <- rows |> dplyr::group_by(family, outcome, model_index) |>
    dplyr::slice_max(rank, n = 1, with_ties = FALSE) |> dplyr::ungroup() |>
    dplyr::arrange(family, outcome, model_index)
  expected <- dplyr::bind_rows(specifications)[c("family", "outcome", "model_index")]

  if (nrow(rows) != nrow(expected) || nrow(dplyr::anti_join(expected, rows,
    by = c("family", "outcome", "model_index"))) || !all(rows$diagnostic_pass)) {
    stop("Every interaction candidate must pass diagnostics before LOO. No automatic exclusions.")
  }
  for (i in seq_len(nrow(rows))) {
    reference <- checkpoints[[rows$family[[i]]]]$data[[rows$outcome[[i]]]]
    if (rows$data_sha256[[i]] != iad_fingerprint(reference)) stop("Accepted fit row/exposure/moderator mismatch.")
  }

  purrr::pmap(rows, list)
}

# Integrated and exact holdout scores ----

# Compute observation-level integrated PSIS-LOO on the parent sample. Use the
# exact refit helpers for observations with problematic importance weights.
iad_compute_loo <- function(record, phase = "loo") {
  iad_require_phase(phase)
  iad_check_hash(record$fit_path, record$sha256)
  candidate <- qs::qread(record$fit_path)
  # The existing integration helper chooses the study random-effect structure
  # by original candidate index. Supply that parent index only to the adapter;
  # the saved fit and all returned interaction identities remain unique.
  candidate$model_index <- as.integer(record$parent_index)
  initial <- compute_primary_candidate_psis_loo(candidate, cores = 1L, save_psis = FALSE)
  result <- compute_primary_candidate_integrated_psis_loo(candidate, initial, cores = 1L, k_threshold = 0.7)
  result <- iad_relabel_loo(result, record)
  result$fit_path <- record$fit_path; result$fit_sha256 <- record$sha256
  result$fit_stage <- record$stage
  iad_qsave(result, iad_path(record$family, record$outcome, "loo", "integrated",
    sprintf("model_%02d.qs", record$model_index), archive = TRUE))
}

iad_relabel_loo <- function(result, record) {
  result$model_index <- as.integer(record$model_index)
  result$model_name <- record$model_name
  result$family <- record$family; result$parent_index <- as.integer(record$parent_index)
  result$interaction <- record$interaction
  for (name in c("comparison", "observation_diagnostics", "exact_refit_summary")) {
    if (is.data.frame(result[[name]]) && nrow(result[[name]])) {
      result[[name]]$model_index <- as.integer(record$model_index)
      result[[name]]$model_name <- record$model_name
      result[[name]]$family <- record$family
    }
  }
  result
}

iad_loo_audit <- function(paths) {
  rows <- list(); cases <- list(); summary <- list()
  for (path in unname(unlist(paths))) {
    r <- qs::qread(path)
    k <- loo::pareto_k_values(r$loo)
    if (anyNA(k) || any(!is.finite(r$loo$pointwise[, "elpd_loo"]))) stop("Invalid integrated interaction LOO.")
    rows[[length(rows) + 1L]] <- tibble::tibble(family = r$family, outcome = r$outcome,
      model_index = r$model_index, model_name = r$model_name, parent_index = r$parent_index,
      interaction = r$interaction, fit_path = r$fit_path, fit_sha256 = r$fit_sha256,
      fit_stage = r$fit_stage, loo_path = path, loo_sha256 = iad_hash(path))
    for (observation in which(k > 0.7)) cases[[length(cases) + 1L]] <- c(as.list(rows[[length(rows)]]),
      list(observation = as.integer(observation), case_id = length(cases) + 1L))
    summary[[length(summary) + 1L]] <- tibble::tibble(family = r$family, outcome = r$outcome,
      model_index = r$model_index, observations = length(k), high_k_cases = sum(k > 0.7), maximum_k = max(k))
  }
  list(manifest = dplyr::bind_rows(rows), summary = dplyr::bind_rows(summary),
       cases = if (length(cases)) cases else list(list(skip = TRUE)))
}

iad_write_loo_audit <- function(audit, family = NULL) {
  root <- if (is.null(family)) iad_path("loo") else iad_path(family, "loo")
  c(iad_csv(audit$manifest, file.path(root, "manifest.csv")),
    iad_csv(audit$summary, file.path(root, "summary.csv")),
    iad_csv(dplyr::bind_rows(audit$cases), file.path(root, "exact_cases.csv")))
}

iad_exact <- function(record, phase = "exact", family = NULL, compact = FALSE) {
  iad_require_phase(phase)
  archive <- if (compact) "exact_checkpoints" else "exact_fits"
  if (isTRUE(record$skip)) return(iad_qsave(record,
    if (is.null(family)) iad_path("loo", archive, "none.qs", archive = TRUE) else
      iad_path(family, "loo", archive, "none.qs", archive = TRUE)))

  iad_check_hash(record$fit_path, record$fit_sha256)
  iad_check_hash(record$loo_path, record$loo_sha256)
  path <- iad_path(record$family, record$outcome, "loo", archive,
    sprintf("model_%02d_observation_%04d.qs", record$model_index, record$observation), archive = TRUE)
  if (file.exists(path)) {
    previous <- qs::qread(path)
    if (!identical(previous$interaction_record, record)) stop("An existing exact refit belongs to different inputs.")
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  candidate <- qs::qread(record$fit_path)
  integrated <- qs::qread(record$loo_path)
  if (!isTRUE(loo::pareto_k_values(integrated$loo)[[record$observation]] > 0.7)) stop("Invalid exact case.")
  refit <- fit_primary_integrated_exact_loo_case(c(record,
    list(fitted_candidate = candidate, integrated_result = integrated)))
  refit$interaction_record <- record

  if (compact) {
    rm(candidate, integrated)
    diagnostics <- iad_exact_diagnostic_summary(refit)
    correction <- if (isTRUE(diagnostics$diagnostic_pass)) iad_build_exact_correction(refit) else NULL
    return(iad_qsave(list(compact = TRUE, interaction_record = record,
      diagnostics = diagnostics, correction = correction), path))
  }

  iad_qsave(refit, path)
}

iad_exact_diagnostic_summary <- function(refit) {
  result <- summarise_primary_exact_loo_refit_diagnostics(refit)
  result$family <- refit$interaction_record$family
  result
}

iad_exact_diagnostics <- function(path) {
  refit <- qs::qread(path)
  if (isTRUE(refit$skip)) return(tibble::tibble())
  if (isTRUE(refit$compact)) return(refit$diagnostics)
  iad_exact_diagnostic_summary(refit)
}

iad_build_exact_correction <- function(refit) {
  record <- refit$interaction_record
  prediction_model <- refit$refitted_model
  prediction_model$data <- refit$original_data
  components <- get_primary_integrated_loo_components(list(outcome = record$outcome,
    model_index = as.integer(record$parent_index), model_name = record$model_name, model = prediction_model))
  correction <- compute_update_2026_primary_exact_loo_case(refit, components = components)
  iad_relabel_loo(correction, record)
}

iad_score_exact <- function(path, phase = "correction", family = NULL) {
  iad_require_phase(phase)
  refit <- qs::qread(path)
  if (isTRUE(refit$skip)) return(iad_qsave(refit,
    if (is.null(family)) iad_path("loo", "exact_none.qs", archive = TRUE) else
      iad_path(family, "loo", "exact_none.qs", archive = TRUE)))

  record <- refit$interaction_record
  iad_check_hash(record$fit_path, record$fit_sha256)
  iad_check_hash(record$loo_path, record$loo_sha256)
  diagnostics <- if (isTRUE(refit$compact)) refit$diagnostics else iad_exact_diagnostic_summary(refit)
  if (!isTRUE(diagnostics$diagnostic_pass)) {
    stop("Exact refit diagnostics require a reviewed remedy before scoring.")
  }

  correction <- if (isTRUE(refit$compact)) refit$correction else iad_build_exact_correction(refit)
  if (!identical(correction$family, record$family) || !identical(correction$outcome, record$outcome) ||
      !identical(correction$model_index, as.integer(record$model_index)) ||
      !identical(as.integer(correction$exact_refit_summary$observation), as.integer(record$observation)) ||
      !isTRUE(correction$exact_refit_summary$diagnostic_pass) ||
      !isTRUE(is.finite(correction$exact_refit_summary$exact_elpd_loo))) {
    stop("Invalid exact correction checkpoint.")
  }

  iad_qsave(correction, iad_path(record$family, record$outcome, "loo", "exact",
    sprintf("model_%02d_observation_%04d.qs", record$model_index, record$observation), archive = TRUE))
}

iad_finalise_loo <- function(audit, exact_paths) {
  corrections <- lapply(unname(unlist(exact_paths)), qs::qread)
  corrections <- Filter(function(x) !isTRUE(x$skip), corrections)
  actual <- purrr::map_dfr(corrections, function(x) tibble::tibble(family = x$family,
    outcome = x$outcome, model_index = x$model_index, observation = x$exact_refit_summary$observation))
  expected <- dplyr::bind_rows(audit$cases)
  columns <- c("family", "outcome", "model_index", "observation")

  if ("skip" %in% names(expected)) expected <- actual[FALSE, , drop = FALSE]
  if (nrow(expected) != nrow(actual) || (nrow(expected) &&
      (anyDuplicated(actual[columns]) || nrow(dplyr::anti_join(expected, actual, by = columns))))) {
    stop("Exact corrections do not cover precisely all interaction high-k cases.")
  }
  paths <- character()

  for (i in seq_len(nrow(audit$manifest))) {
    row <- as.list(audit$manifest[i, ])
    iad_check_hash(row$loo_path, row$loo_sha256)
    r <- qs::qread(row$loo_path)
    selected <- Filter(function(x) x$family == row$family && x$outcome == row$outcome &&
      x$model_index == row$model_index, corrections)
    r$exact_refit_summary <- if (length(selected)) purrr::map_dfr(selected, "exact_refit_summary") else
      tibble::tibble(observation = integer(), diagnostic_pass = logical(), exact_elpd_loo = double())
    for (correction in selected) {
      if (!iad_same_observations(correction$observation_manifest, r$observation_manifest) ||
          !all(correction$exact_refit_summary$diagnostic_pass)) stop("Exact refit identity/diagnostic gate failed.")
      obs <- correction$exact_refit_summary$observation
      if (length(obs) != 1L) stop("An exact archive must correct one observation.")
      r$loo$pointwise[obs, ] <- correction$loo$pointwise[obs, ]
      r$loo$diagnostics$pareto_k[obs] <- NA_real_
      if (!is.null(r$loo$diagnostics$n_eff)) r$loo$diagnostics$n_eff[obs] <- NA_real_
      r$observation_diagnostics$integrated_pareto_k[obs] <- NA_real_
    }
    r$observation_diagnostics$exact_loo_refit <- seq_len(nrow(r$observation_manifest)) %in%
      r$exact_refit_summary$observation
    r$loo <- update_primary_exact_loo_estimates(r$loo)
    iad_validate_loo(r)
    paths <- c(paths, iad_qsave(r, iad_path(row$family, row$outcome, "loo", "final",
      sprintf("model_%02d.qs", row$model_index), archive = TRUE)))
  }

  paths
}

iad_stacking_weights <- function(objects) {
  counts <- vapply(objects, function(x) nrow(x$pointwise), integer(1))
  if (!length(counts) || length(unique(counts)) != 1L || any(counts < 1L) ||
      any(!vapply(objects, function(x) all(is.finite(x$pointwise[, "elpd_loo"])), logical(1)))) {
    stop("Stacking requires finite pointwise scores on matching observations.")
  }
  # Stacking consumes pointwise ELPD. Different accepted sampling lengths
  # (including muscle model 11) need compatible draw-count metadata only on
  # these temporary copies, as in the original absolute-primary workflow.
  objects <- lapply(objects, function(x) { attr(x, "dims") <- c(1L, counts[[1L]]); x })
  weights <- as.numeric(loo::loo_model_weights(objects, method = "stacking", cores = 1L))
  if (length(weights) != length(objects) || any(!is.finite(weights)) ||
      any(weights < -1e-8) || abs(sum(weights) - 1) > 1e-8) stop("Invalid stacking weights.")
  pmax(weights, 0) / sum(pmax(weights, 0))
}

# Compare primary and interaction models ----

# Assemble the accepted interaction candidates and their primary parents on
# the same observations, then estimate weights for reporting.
iad_combine <- function(checkpoints, accepted, final_paths, families = iad_families()) {
  iad_require_phase("combine")
  iad_require_parents(checkpoints)
  accepted <- dplyr::bind_rows(accepted)
  new_loo <- lapply(unname(unlist(final_paths)), qs::qread)
  output <- list()

  for (family in families) for (outcome in iad_outcomes()) {
    parents <- checkpoints[[family]]$parents
    parents <- parents[parents$outcome == outcome, , drop = FALSE]
    new <- accepted[accepted$family == family & accepted$outcome == outcome, ]
    additions <- Filter(function(x) x$family == family && x$outcome == outcome, new_loo)
    additions <- additions[order(vapply(additions, `[[`, integer(1), "model_index"))]
    new <- new[order(new$model_index), ]
    if (length(additions) != (if (family == "relative") 36L else 24L) ||
        !identical(vapply(additions, `[[`, integer(1), "model_index"), new$model_index) ||
        !all(new$diagnostic_pass)) stop("Incomplete interaction expansion.")
    results <- lapply(seq_len(nrow(parents)), function(i) {
      iad_check_hash(parents$loo_path[[i]], parents$loo_sha256[[i]])
      qs::qread(parents$loo_path[[i]])
    })
    reference <- checkpoints[[family]]$data[[outcome]]
    for (r in c(results, additions)) iad_validate_loo(r, reference)
    # Validate identity and data against the actual saved parent fits here,
    # rather than assuming an archived reporting manifest describes the fit.
    for (i in seq_len(nrow(parents))) {
      candidate <- qs::qread(parents$fit_path[[i]])
      identity_checks <- c(model_index = isTRUE(candidate$model_index == parents$model_index[[i]]),
        outcome = isTRUE(candidate$outcome == outcome),
        model_name = isTRUE(candidate$model_name == parents$model_name[[i]]),
        observations = iad_same_observations(candidate$model$data, reference),
        distribution = identical(stats::family(candidate$model)$family, "gaussian"),
        link = identical(stats::family(candidate$model)$link, "identity"))
      if (!all(identity_checks)) stop("Parent fit combination gate failed for ", family, "/", outcome,
        "/model_", parents$model_index[[i]], ": ", paste(names(identity_checks)[!identity_checks], collapse = ", "))
      iad_validate_parent_data(candidate, reference, family)
      specification <- checkpoints[[family]]$specifications
      specification <- specification[specification$outcome == outcome &
        specification$model_index == candidate$model_index, ]
      if (nrow(specification) != 1L || !identical(
          gsub("\\s+", "", candidate$specification$formula),
          gsub("\\s+", "", specification$formula[[1L]]))) stop("Parent fitted formula differs from its frozen specification.")
      if (!"fit_sha256" %in% names(parents)) parents$fit_sha256 <- rep(NA_character_, nrow(parents))
      parents$fit_sha256[[i]] <- iad_hash(parents$fit_path[[i]])
      rm(candidate); gc(verbose = FALSE)
    }
    all_results <- c(results, additions)
    models <- dplyr::bind_rows(parents,
      dplyr::transmute(new, family, outcome, model_index, parent_index, model_name, interaction,
        fit_path, stage, fit_sha256 = sha256))
    if (anyDuplicated(models$model_index) || anyDuplicated(models$model_name)) stop("Duplicate combined candidate.")
    objects <- lapply(all_results, `[[`, "loo"); names(objects) <- models$model_name
    w <- iad_stacking_weights(objects)
    if (any(!is.finite(w)) || any(w < -1e-8) || abs(sum(w) - 1) > 1e-8) stop("Invalid combined stacking weights.")
    models$stacking_weight <- pmax(w, 0) / sum(pmax(w, 0))
    old_w <- iad_stacking_weights(objects[seq_len(nrow(parents))])
    if (any(!is.finite(old_w)) || any(old_w < -1e-8)) stop("Invalid primary-only reference weights.")
    old_w <- pmax(old_w, 0) / sum(pmax(old_w, 0))
    models$primary_only_weight <- c(old_w, rep(0, nrow(new)))
    models$fit_sha256[models$interaction != "none"] <- new$sha256
    output[[paste(family, outcome, sep = "_")]] <- list(family = family, outcome = outcome,
      models = models, seed = 931988L + 100L * match(family, iad_families()) + match(outcome, iad_outcomes()),
      diagnostics = new,
      loo_summary = purrr::map_dfr(all_results, function(r) dplyr::mutate(
        summarise_primary_loo_method(r$loo, "integrated_observation_level_psis_plus_exact"),
        model_index = r$model_index, model_name = r$model_name, .before = 1L)),
      gate_passed = TRUE)
  }

  output
}

iad_write_combination <- function(records, family = NULL) {
  table <- purrr::map_dfr(records, "models")
  root <- if (is.null(family)) iad_path("combination") else iad_path(family, "combination")
  c(iad_csv(table, file.path(root, "stacking_weights.csv")),
    iad_csv(table |> dplyr::group_by(family, outcome, interaction) |>
      dplyr::summarise(candidates = dplyr::n(), stacking_weight = sum(stacking_weight), .groups = "drop"),
      file.path(root, "summary.csv")))
}

# Predict and report at moderator profiles ----

# Use the dose and moderator reference grid to compare profile-specific
# predictions and increments with their empirical support.
iad_prediction_grid <- function(data, family, profile) {
  grid <- if (family == "relative") make_primary_prediction_grid(data) else
    make_update_2026_total_primary_prediction_grid(data,
      exposure_minimum_g_day = if (family == "total") 50 else 0)
  reference <- if (family == "relative") get_primary_protein_median(data) else
    unique(data$absolute_exposure_reference_g_day)
  # Retain the frozen dose reference explicitly for paired dose contrasts.
  if (!any(abs(grid$keys$protein_intake - reference) < 1e-10)) {
    for (rt in c(0, 1)) {
      position <- which(grid$keys$resistance_exercise_code == rt)[[1L]]
      row <- grid$newdata[position, , drop = FALSE]
      key <- grid$keys[position, , drop = FALSE]
      key$protein_intake <- reference
      if ("protein_intake_g_day" %in% names(key)) key$protein_intake_g_day <- reference
      row$protein_intake_raw <- reference
      row$protein_intake <- 0
      grid$newdata <- dplyr::bind_rows(grid$newdata, row)
      grid$keys <- dplyr::bind_rows(grid$keys, key)
    }
    order <- order(grid$keys$resistance_exercise_code, grid$keys$protein_intake)
    grid$newdata <- grid$newdata[order, ]; grid$keys <- grid$keys[order, ]
  }

  grid$newdata$Intercept <- 1
  grid$newdata$age_centre <- profile$age - 25
  grid$newdata$duration_centre <- profile$weeks - 12
  moderator_data <- add_primary_standardised_moderators(data)
  centres <- attr(moderator_data, "moderator_centres")
  scales <- attr(moderator_data, "moderator_scales")
  grid$newdata$age_outcome_z <- (profile$age - 25 - centres[["age_centre"]]) / scales[["age_centre"]]
  grid$newdata$duration_outcome_z <- (profile$weeks - 12 - centres[["duration_centre"]]) / scales[["duration_centre"]]
  grid$newdata$protein_intake_05 <- grid$newdata$protein_intake / 0.5

  if (family != "relative") grid$newdata$body_mass_outcome_z <- 0
  grid$newdata <- iad_augment_data(grid$newdata)
  grid$reference <- reference
  grid
}

iad_profile_support <- function(data, profile_row) {
  purrr::map_dfr(c(0, 1), function(rt) {
    rows <- data[data$resistance_exercise_code == rt, ]
    local <- rows[abs(rows$age_centre + 25 - profile_row$age) <= profile_row$age_neighbourhood_years &
      abs(rows$duration_centre + 12 - profile_row$weeks) <= profile_row$duration_neighbourhood_weeks, ]
    tibble::tibble(resistance_exercise_code = rt, profile = profile_row$profile,
      age = profile_row$age, weeks = profile_row$weeks,
      local_studies = dplyr::n_distinct(local$study_no), local_arms = dplyr::n_distinct(local$arm_no),
      local_min = if (nrow(local)) min(local$protein_intake_raw) else NA_real_,
      local_max = if (nrow(local)) max(local$protein_intake_raw) else NA_real_,
      age_within_range = profile_row$age >= min(rows$age_centre + 25) && profile_row$age <= max(rows$age_centre + 25),
      duration_within_range = profile_row$weeks >= min(rows$duration_centre + 12) && profile_row$weeks <= max(rows$duration_centre + 12))
  })
}

iad_summarise_draws <- function(draws) {
  limits <- matrixStats::colQuantiles(draws, probs = c(0.025, 0.975), drop = FALSE)
  tibble::tibble(estimate = colMeans(draws), conf_low = limits[, 1L], conf_high = limits[, 2L],
    probability_meaningful_gain = colMeans(draws > 0.1),
    probability_practical_equivalence = colMeans(abs(draws) <= 0.1),
    probability_meaningful_loss = colMeans(draws < -0.1))
}

iad_model_allocation <- function(models, n_draws, seed, column = "stacking_weight") {
  weights <- models[[column]]
  if (!length(weights) || any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) stop("Invalid posterior-mixture allocation weights.")
  set.seed(seed)
  sample(rep(seq_len(nrow(models)), times = allocate_primary_mixture_models(weights, n_draws)))
}

# Produce predictions for one outcome and family. Sample from the combined
# ensemble and the primary-only ensemble, then write the dose curves, support
# checks and moderator contrasts.
iad_report_scope <- function(record, checkpoints, profiles) {
  iad_require_phase("report")
  if (!isTRUE(record$gate_passed)) stop("Combined ensemble has not passed its gate.")
  family <- record$family; outcome <- record$outcome; models <- record$models
  data <- checkpoints[[family]]$data[[outcome]]
  profiles <- profiles[profiles$family == family & profiles$outcome == outcome, ]

  # Set up the same prediction profiles for the combined and primary-only
  # ensembles so the resulting increments can be compared.
  grids <- lapply(seq_len(nrow(profiles)), function(i) iad_prediction_grid(data, family, profiles[i, ]))
  n_draws <- 4000L
  allocation <- iad_model_allocation(models, n_draws, record$seed)
  old_allocation <- iad_model_allocation(models, n_draws, record$seed, "primary_only_weight")
  width <- nrow(grids[[1L]]$newdata)
  newdata <- dplyr::bind_rows(lapply(grids, `[[`, "newdata"))
  draws <- old_draws <- matrix(NA_real_, n_draws, nrow(newdata))

  for (i in seq_len(nrow(models))) {
    selected <- which(allocation == i); old_selected <- which(old_allocation == i)
    count <- max(length(selected), length(old_selected))
    if (!count) next
    iad_check_hash(models$fit_path[[i]], models$fit_sha256[[i]])
    candidate <- qs::qread(models$fit_path[[i]])
    if (candidate$model_index != models$model_index[[i]] || candidate$model_name != models$model_name[[i]] ||
        candidate$outcome != outcome || !iad_same_observations(candidate$model$data, data)) {
      stop("Prediction fit identity or data mismatch.")
    }
    iad_validate_parent_data(candidate, data, family)
    set.seed(record$seed + models$model_index[[i]])
    ids <- sample.int(posterior::ndraws(candidate$model), count)
    prediction <- as.matrix(brms::posterior_epred(candidate$model, newdata = newdata,
      re_formula = NA, draw_ids = ids))
    if (!identical(dim(prediction), c(count, nrow(newdata))) || any(!is.finite(prediction))) {
      stop("Interaction posterior prediction shape mismatch.")
    }
    if (length(selected)) draws[selected, ] <- prediction[seq_along(selected), , drop = FALSE]
    if (length(old_selected)) old_draws[old_selected, ] <- prediction[seq_along(old_selected), , drop = FALSE]
    rm(candidate, prediction); gc(verbose = FALSE)
  }

  if (anyNA(draws) || anyNA(old_draws)) stop("Incomplete interaction posterior mixture.")
  paths <- character(); dose_tables <- support_tables <- list(); increment_draws <- list()
  root <- iad_path(family, outcome, "reporting")

  # Summarize each profile on its own observed support before comparing
  # increments across profiles.
  for (i in seq_len(nrow(profiles))) {
    profile <- profiles[i, ]; columns <- (i - 1L) * width + seq_len(width)
    grid <- grids[[i]]
    average <- list(outcome = outcome, draws = draws[, columns, drop = FALSE], grid = grid$keys,
      model_weights = stats::setNames(models$stacking_weight, models$model_name),
      settings = list(grid_step = grid$grid_step, protein_min = grid$protein_min, sesoi = 0.1,
        conf_level = 0.95, n_draws = n_draws, age = profile$age, weeks = profile$weeks))
    class(average) <- "primary_model_average"
    support <- iad_profile_support(data, profile)
    dose <- summarise_primary_dose_response(average) |>
      dplyr::rename(conf_low = conf.low, conf_high = conf.high) |>
      dplyr::left_join(support, by = "resistance_exercise_code")
    dose$locally_supported <- dose$local_studies >= 3L & dose$age_within_range & dose$duration_within_range &
      dose$protein_intake >= dose$local_min & dose$protein_intake <= dose$local_max
    dose$locally_supported[is.na(dose$locally_supported)] <- FALSE
    old_average <- average; old_average$draws <- old_draws[, columns, drop = FALSE]
    old_dose <- summarise_primary_dose_response(old_average) |>
      dplyr::rename(conf_low = conf.low, conf_high = conf.high)
    old_dose$profile <- profile$profile
    old_dose$ensemble <- "Primary only (matched frozen sample)"
    dose$ensemble <- "Primary plus age/duration interactions"
    dose_tables[[i]] <- dplyr::bind_rows(dose, old_dose)
    support_tables[[i]] <- support
    contrasts <- compute_primary_heatmap_contrasts(average, data, protein_min = grid$protein_min)
    table <- contrasts$contrasts |> dplyr::left_join(support, by = "resistance_exercise_code")
    table$supported <- table$supported & table$local_studies >= 3L & table$age_within_range &
      table$duration_within_range & table$protein_intake_current >= table$local_min &
      table$protein_intake_proposed >= table$local_min & table$protein_intake_current <= table$local_max &
      table$protein_intake_proposed <= table$local_max
    table$supported[is.na(table$supported)] <- FALSE
    table$support_status <- ifelse(table$supported, "Within local empirical range", "Limited local support")
    table$profile <- profile$profile
    paths <- c(paths, iad_csv(table, file.path(root, paste0(profile$profile, "_pairwise.csv.gz"))))
    local <- if (any(table$supported)) summarise_primary_local_equivalence(table, outcome) else
      tibble::tibble(outcome = character(), rt_condition = character(), median_total_span = double())
    local$profile <- rep(profile$profile, nrow(local))
    paths <- c(paths, iad_csv(local, file.path(root, paste0(profile$profile, "_local_equivalence.csv.gz"))))
    increments <- average$draws
    for (rt in c(0, 1)) {
      at <- which(grid$keys$resistance_exercise_code == rt &
        abs(grid$keys$protein_intake - grid$reference) < 1e-10)
      cols <- which(grid$keys$resistance_exercise_code == rt)
      if (length(at) != 1L) stop("Missing unique reference-dose prediction.")
      increments[, cols] <- sweep(average$draws[, cols, drop = FALSE], 1L, average$draws[, at], "-")
    }
    increment_draws[[profile$profile]] <- increments
  }
  # Paired differences of protein contrasts across profiles use the same
  # component and posterior draw, preserving their posterior covariance.
  comparisons <- list(age = c("age_q75", "age_q25"), duration = c("duration_q75", "duration_q25"))
  moderation <- purrr::imap_dfr(comparisons, function(pair, moderator) {
    difference <- increment_draws[[pair[[1L]]]] - increment_draws[[pair[[2L]]]]
    dplyr::bind_cols(grids[[1L]]$keys, iad_summarise_draws(difference)) |>
      dplyr::mutate(moderator = moderator, higher_profile = pair[[1L]], lower_profile = pair[[2L]],
        current_intake = grids[[1L]]$reference, .before = 1L)
  })
  doses <- dplyr::bind_rows(dose_tables)
  support <- dplyr::bind_rows(support_tables)
  moderation <- moderation |>
    dplyr::left_join(dplyr::rename(support, higher_profile = profile),
      by = c("resistance_exercise_code", "higher_profile")) |>
    dplyr::rename(higher_local_studies = local_studies, higher_min = local_min, higher_max = local_max) |>
    dplyr::left_join(dplyr::rename(support, lower_profile = profile),
      by = c("resistance_exercise_code", "lower_profile"), suffix = c("_higher", "_lower"))
  moderation$supported <- moderation$higher_local_studies >= 3L & moderation$local_studies >= 3L &
    moderation$current_intake >= moderation$higher_min & moderation$current_intake <= moderation$higher_max &
    moderation$current_intake >= moderation$local_min & moderation$current_intake <= moderation$local_max &
    moderation$protein_intake >= moderation$higher_min & moderation$protein_intake <= moderation$higher_max &
    moderation$protein_intake >= moderation$local_min & moderation$protein_intake <= moderation$local_max &
    moderation$age_within_range_higher & moderation$age_within_range_lower &
    moderation$duration_within_range_higher & moderation$duration_within_range_lower
  moderation$supported[is.na(moderation$supported)] <- FALSE
  allocation <- tibble::tibble(model_index = models$model_index, model_name = models$model_name,
    draws = tabulate(allocation, nrow(models)), primary_only_draws = tabulate(old_allocation, nrow(models)))
  paths <- c(paths, iad_csv(doses, file.path(root, "dose_response.csv.gz")),
    iad_csv(support, file.path(root, "profile_support.csv")),
    iad_csv(moderation, file.path(root, "interaction_contrasts.csv.gz")),
    iad_csv(allocation, file.path(root, "draw_allocation.csv")),
    iad_csv(models, file.path(root, "stacking_weights.csv")))
  # Small publication interface: static tables/figures, no fit or targets access.
  interface <- here::here("manuscript", "generated", "age-duration-interactions", family, outcome)
  paths <- c(paths,
    iad_csv(doses, file.path(interface, "dose_response.csv.gz")),
    iad_csv(support, file.path(interface, "profile_support.csv")),
    iad_csv(moderation, file.path(interface, "interaction_contrasts.csv.gz")))
  figure <- ggplot2::ggplot(doses, ggplot2::aes(x = protein_intake, y = estimate,
    colour = ensemble, fill = ensemble)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf_low, ymax = conf_high), alpha = 0.1, colour = NA) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::facet_grid(profile ~ resistance_exercise) + ggplot2::theme_bw() +
    ggplot2::labs(x = if (family == "relative") "Protein intake (g/kg/day)" else "Protein dose (g/day)",
      y = "Expected standardised baseline-to-follow-up change", title = paste(family, outcome),
      caption = "Population means; random effects excluded. Curves may extend beyond local support; inspect profile_support.csv.") +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())

  figure_path <- file.path(interface, "conditional_predictions.png")
  ggplot2::ggsave(figure_path, figure, width = 10, height = 14, dpi = 180)
  paths <- c(paths, normalizePath(figure_path, winslash = "/", mustWork = TRUE))
  public_weights <- dplyr::select(models, family, outcome, model_index, parent_index,
    model_name, interaction, stacking_weight, primary_only_weight)
  diagnostics <- dplyr::select(record$diagnostics, -dplyr::any_of("fit_path"))
  specifications <- dplyr::bind_rows(iad_specifications(checkpoints))
  specifications <- specifications[specifications$family == family & specifications$outcome == outcome, ]
  artifacts <- list(verification_checks = tibble::tibble(check = c("combined gate passed",
    "4000 mixture draws", "all conditional means finite", "six profiles retained", "posterior draws paired across profiles"),
    passed = c(record$gate_passed, sum(allocation$draws) == 4000L,
      all(is.finite(doses$estimate)), nrow(profiles) == 6L, TRUE)),
    source_manifest = checkpoints[[family]]$source_manifest,
    initial_interaction_specifications = specifications,
    accepted_interaction_diagnostics = diagnostics, final_loo_summary = record$loo_summary,
    stacking_weights = public_weights, profiles = profiles, profile_support = support,
    dose_response_predictions = doses, interaction_contrasts = moderation, draw_allocation = allocation)

  if (!all(artifacts$verification_checks$passed)) stop("Reporting verification failed.")
  paths <- c(paths, write_verification_bundle(paste0("age_duration_", family, "_", outcome, "_update_2026"),
    artifacts, descriptions = stats::setNames(names(artifacts), names(artifacts))))
  manifest <- tibble::tibble(path = paths, bytes = file.info(paths)$size,
    sha256 = vapply(paths, iad_hash, character(1)))
  c(paths, iad_csv(manifest, file.path(root, "file_manifest.csv")))
}

iad_write_supplement <- function(report_paths, records) {
  iad_require_phase("report")
  paths <- unname(unlist(report_paths))
  if (length(records) != 12L || any(!file.exists(paths))) stop("All twelve reporting scopes must finish first.")
  manifests <- paths[basename(paths) == "file_manifest.csv"]
  if (length(manifests) != 12L) stop("Missing supplementary reporting manifests.")

  for (manifest_path in manifests) {
    manifest <- iad_read_csv(manifest_path)
    for (i in seq_len(nrow(manifest))) iad_check_hash(manifest$path[[i]], manifest$sha256[[i]])
  }
  lines <- c("<!-- Generated by the age-duration interaction pipeline; do not edit. -->",
    "", "All twelve primary-family expansions passed the combination and reporting gates.", "",
    "These exploratory analyses were specified after registration. Stacking weights are estimated anew within each exposure and outcome; they are not posterior probabilities of moderator hypotheses.", "",
    "The plots show expected standardised baseline-to-follow-up changes with random effects excluded. Age is an arm/study mean. Profiles use the frozen common sample: age quartiles give each distinct arm one observation, and duration quartiles give each distinct study/follow-up one observation. Body mass is held at its outcome mean for absolute-exposure moderator candidates.", "",
    "Local support is a descriptive screen: at least three studies within 10 years and four weeks of a profile, with both doses inside that neighbourhood's observed range. It is not a guarantee of joint positivity or a causal comparison. Full support counts and paired differences of protein contrasts are provided alongside each plot.", "")

  for (record in records) {
    family <- record$family; outcome <- record$outcome
    weights <- record$models
    lines <- c(lines, paste0("### ", family, ": ", outcome), "",
      sprintf("Combined candidates: %d. Stacking weight on added interaction candidates: %.3f.",
        nrow(weights), sum(weights$stacking_weight[weights$interaction != "none"])), "",
      paste0("![](../manuscript/generated/age-duration-interactions/", family, "/", outcome, "/conditional_predictions.png)"), "",
      paste0("[Conditional predictions](../manuscript/generated/age-duration-interactions/", family, "/", outcome, "/dose_response.csv.gz), ",
        "[profile support](../manuscript/generated/age-duration-interactions/", family, "/", outcome, "/profile_support.csv), ",
        "[paired interaction contrasts](../manuscript/generated/age-duration-interactions/", family, "/", outcome, "/interaction_contrasts.csv.gz)."), "")
  }

  path <- here::here("manuscript", "generated", "age-duration-interactions", "_results.qmd")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(lines, temporary)
  on.exit(unlink(temporary), add = TRUE)
  if (!file.rename(temporary, path)) stop("Could not atomically materialise supplement interface.")

  path
}

