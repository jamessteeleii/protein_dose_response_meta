# Deferred measurement restrictions for the primary muscle-outcome ensemble.
# Selection is made from the checksummed, frozen 2026 primary rows. No original
# workbook, fitted primary model, or existing sensitivity output is modified.

measurement_variants <- function() c("prereg_methods", "direct_only")

measurement_phase_enabled <- function(phase) {
  phases <- strsplit(tolower(Sys.getenv(
    "PROTEIN_MEASUREMENT_PHASES", unset = "none"
  )), "[, ]+")[[1L]]
  allowed <- c("none", "pilot", "fit", "remedy", "loo", "exact", "report")
  if (any(!phases %in% allowed)) {
    stop("Unknown PROTEIN_MEASUREMENT_PHASES value.", call. = FALSE)
  }
  phase %in% phases
}

measurement_path <- function(variant, stage, filename, archive = TRUE) {
  stopifnot(variant %in% measurement_variants())
  root <- if (archive) "model_archive" else "results"
  here::here(root, "update_2026", "measurement", variant, stage, filename)
}

write_measurement_qs <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_atomic_qs(object, path)
}

measurement_checkpoint_inputs <- function() {
  here::here("results", "verification", "primary_update_2026", c(
    "manifest.csv", "analysis_data.csv.gz", "model_specifications.csv.gz"
  ))
}

read_measurement_checkpoint <- function(paths) {
  expected <- c("manifest.csv", "analysis_data.csv.gz", "model_specifications.csv.gz")
  if (!identical(basename(paths), expected) || any(!file.exists(paths))) {
    stop("Measurement analysis requires the frozen primary checkpoint.", call. = FALSE)
  }
  manifest <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE)
  for (path in paths[-1L]) {
    record <- manifest[manifest$file == basename(path), , drop = FALSE]
    if (nrow(record) != 1L || !identical(
      digest::digest(file = path, algo = "sha256"), record$sha256[[1L]]
    )) stop("Primary checkpoint checksum mismatch: ", basename(path), call. = FALSE)
  }
  data <- readr::read_csv(paths[[2L]], show_col_types = FALSE, guess_max = Inf)
  specifications <- readr::read_csv(paths[[3L]], show_col_types = FALSE, guess_max = Inf)
  if (nrow(readr::problems(data)) || nrow(readr::problems(specifications))) {
    stop("Frozen primary checkpoint could not be parsed without data loss.")
  }
  # CSV does not preserve factors. Keep the primary RT coding (yes = 0) and
  # rebuild study-factor levels from the retained rows for each restricted set.
  data$resistance_exercise_factor <- factor(
    data$resistance_exercise_factor, levels = c("yes", "no")
  )
  if (anyNA(data$resistance_exercise_factor) || any(
    data$resistance_exercise_code != ifelse(data$resistance_exercise_factor == "yes", 0, 1)
  )) stop("Frozen primary RT coding is inconsistent.", call. = FALSE)
  specifications <- specifications[specifications$outcome == "hypertrophy", , drop = FALSE]
  specifications <- specifications[order(specifications$model_index), , drop = FALSE]
  if (!identical(as.integer(specifications$model_index), seq_len(26L))) {
    stop("Expected all 26 frozen primary muscle specifications.", call. = FALSE)
  }
  # Historical absolute archive paths are not inputs to these new fits.
  specifications$accepted_fit_path <- NULL
  specifications$accepted_fit_stage <- NULL
  list(data = data, specifications = specifications, source_manifest = manifest[
    manifest$file %in% expected[-1L], , drop = FALSE
  ])
}

read_measurement_annotations <- function(path) {
  annotations <- utils::read.csv(path, colClasses = "character", stringsAsFactors = FALSE)
  required <- c("study_no", "arm_no", "outcome_no", "time1", "time2",
    "original_measurement_method", "reviewed_measurement_method",
    "confirmation_date", "confirmation_source")
  if (!setequal(names(annotations), required) || anyNA(annotations) ||
      any(!nzchar(trimws(as.matrix(annotations))))) {
    stop("Measurement annotations require complete identifiers, labels, and provenance.")
  }
  annotations
}

review_measurement_methods <- function(muscle, annotations = NULL) {
  reviewed <- muscle$measurement_method
  applied <- rep(FALSE, nrow(muscle))
  provenance <- confirmed_on <- rep("", nrow(muscle))
  if (!is.null(annotations) && nrow(annotations)) {
    identifiers <- c("study_no", "arm_no", "outcome_no", "time1", "time2")
    key <- function(data) do.call(paste, c(lapply(data[identifiers], as.character), sep = "|"))
    annotation_keys <- key(annotations)
    source_keys <- key(muscle)
    if (anyDuplicated(annotation_keys)) stop("Duplicate measurement annotation identifiers.")
    positions <- match(annotation_keys, source_keys)
    if (anyNA(positions) || any(vapply(annotation_keys,
      function(x) sum(source_keys == x), integer(1)) != 1L)) {
      stop("Each measurement annotation must identify exactly one frozen muscle row.")
    }
    original <- tolower(trimws(annotations$original_measurement_method))
    replacement <- tolower(trimws(annotations$reviewed_measurement_method))
    if (any(tolower(trimws(reviewed[positions])) != original) ||
        any(!replacement %in% c("dexa - appendicular", "dxa - appendicular", "bia - appendicular")) ||
        any(sub("^dexa", "dxa", original) != sub("^dexa", "dxa", sub(" -.*$", "", replacement)))) {
      stop("Appendicular annotations must match the original method and retain DXA/BIA proxy status.")
    }
    reviewed[positions] <- annotations$reviewed_measurement_method
    applied[positions] <- TRUE
    provenance[positions] <- annotations$confirmation_source
    confirmed_on[positions] <- annotations$confirmation_date
  }
  data.frame(measurement_method_reviewed = reviewed,
    measurement_annotation_applied = applied,
    measurement_annotation_source = provenance,
    measurement_confirmation_date = confirmed_on, stringsAsFactors = FALSE)
}

select_measurement_data <- function(data, variant, annotations = NULL) {
  stopifnot(variant %in% measurement_variants())
  required <- c("outcome_group", "measurement_method", "study_no", "arm_no",
                "outcome_no", "time1", "time2", "yi", "vi")
  if (!all(required %in% names(data))) stop("Measurement data fields are missing.")
  muscle <- data[data$outcome_group == "direct muscle size/mass", , drop = FALSE]
  review <- review_measurement_methods(muscle, annotations)
  method <- tolower(trimws(review$measurement_method_reviewed))
  if (anyNA(method) || any(!nzchar(method))) stop("Missing measurement method.")
  direct <- grepl("^(mri|ct|computed tomography|ultrasound|biopsy)( -|$)", method)
  composition <- grepl("^(dexa|dxa|bia)( -|$)", method)
  trunk <- grepl("(^|[^a-z])trunk([^a-z]|$)", method)
  skinfold <- grepl("^skinfold", method)
  # An allow-list prevents an unrecognised proxy from silently entering either
  # analysis. MRI/CT/US/biopsy eligibility remains within the muscle group.
  known <- direct | composition | skinfold
  if (any(!known)) stop("Unreviewed muscle measurement method: ",
    paste(unique(method[!known]), collapse = ", "), call. = FALSE)
  keep <- if (variant == "direct_only") direct else (direct | composition) & !trunk
  reason <- ifelse(keep, "retained", ifelse(
    trunk, "trunk regional estimate", ifelse(skinfold, "skinfold estimate", "DXA/BIA proxy")
  ))
  audit <- muscle[, required, drop = FALSE]
  audit <- dplyr::bind_cols(audit, review)
  audit$source_row <- which(data$outcome_group == "direct muscle size/mass")
  audit$retained <- keep
  audit$decision <- reason
  retained <- muscle[keep, , drop = FALSE]
  # New reporting columns preserve both the frozen label and author-confirmed
  # interpretation. The source label and every numerical model input are intact.
  retained <- dplyr::bind_cols(retained, review[keep, , drop = FALSE])
  if (!nrow(retained) || length(unique(retained$resistance_exercise_code)) != 2L) {
    stop("Both RT strata require retained observations.", call. = FALSE)
  }
  retained$study_factor <- factor(retained$study_no)
  # Preserve the original exposure centring, effect sizes, variances, and IDs.
  # Recompute only the point-size display field, which is not a model input.
  retained <- add_primary_display_fields(retained)
  list(variant = variant, data = retained, selection = audit,
       annotations = if (is.null(annotations)) data.frame() else annotations,
       summary = tibble::tibble(
         variant = variant, original_rows = nrow(muscle), retained_rows = nrow(retained),
         excluded_rows = sum(!keep), studies = length(unique(retained$study_no)),
         arms = length(unique(retained$arm_no))
       ))
}

make_measurement_specifications <- function(checkpoint, variant) {
  table <- checkpoint$specifications
  records <- lapply(seq_len(nrow(table)), function(i) {
    record <- as.list(table[i, , drop = FALSE])
    record$model_index <- as.integer(record$model_index)
    record$measurement_variant <- variant
    record$measurement_stage <- "initial"
    record
  })
  names(records) <- sprintf("model_%02d", seq_len(26L))
  records
}

measurement_input_preflight <- function(selection, specifications, priors) {
  checks <- purrr::map_dfr(specifications, function(specification) {
    prepared <- prepare_update_2026_primary_candidate_input(
      specification, list(hypertrophy = selection$data), list(hypertrophy = priors)
    )
    warning_messages <- character()
    stan_data <- withCallingHandlers(brms::make_standata(
      stats::as.formula(prepared$formula), data = prepared$data, prior = prepared$priors
    ), warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    tibble::tibble(variant = selection$variant,
      model_index = specification$model_index, observations = stan_data$N, valid = TRUE,
      warnings = paste(unique(warning_messages), collapse = " | "))
  })
  if (nrow(checks) != 26L || !all(checks$valid) ||
      any(checks$observations != nrow(selection$data))) stop("Measurement preflight failed.")
  checks
}

write_measurement_preparation <- function(selection, specifications, preflight) {
  directory <- dirname(measurement_path(selection$variant, "preparation", "x", FALSE))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  tables <- list(selection = selection$selection, annotations = selection$annotations,
                 counts = selection$summary,
                 specifications = dplyr::bind_rows(specifications), preflight = preflight)
  purrr::imap_chr(tables, function(table, name) {
    path <- file.path(directory, paste0(name, ".csv"))
    readr::write_csv(table, path)
    path
  }) |> unname()
}

measurement_fit_paths <- function(paths) {
  paths <- unname(unlist(paths, use.names = FALSE))
  paths[grepl("\\.qs$", paths)]
}

measurement_pilot_specifications <- function(specifications) {
  lapply(specifications[c(2L, 26L)], function(x) {
    x$measurement_stage <- "pilot"
    x
  })
}

fit_measurement_candidate <- function(specification, selection, priors, preflight,
                                      pilot_audit = NULL) {
  stage <- specification$measurement_stage
  phase <- if (stage == "initial") "fit" else if (stage == "pilot") "pilot" else "remedy"
  if (!measurement_phase_enabled(phase)) {
    stop("Measurement ", phase, " is disabled; run only after the extended families finish.")
  }
  if (nrow(preflight) != 26L || !all(preflight$valid)) stop("Preflight must pass before fitting.")
  if (stage == "initial" && !isTRUE(pilot_audit$all_pass)) {
    stop("Time and diagnose candidates 2 and 26 before full measurement fan-out.")
  }
  started <- proc.time()[["elapsed"]]
  path <- fit_update_2026_primary_candidate_to_file(
    specification, list(hypertrophy = selection$data), list(hypertrophy = priors),
    output_root = dirname(measurement_path(selection$variant, stage, "x"))
  )
  timing <- sub("\\.qs$", "-timing.csv", path)
  readr::write_csv(tibble::tibble(variant = selection$variant,
    model_index = specification$model_index, stage = stage,
    elapsed_seconds = proc.time()[["elapsed"]] - started), timing)
  c(path, timing)
}

measurement_pilot_audit <- function(paths, variant) {
  paths <- measurement_fit_paths(paths)
  diagnostics <- lapply(paths, diagnose_update_2026_primary_fit_file)
  summary <- dplyr::bind_rows(purrr::map(diagnostics, "summary")) |>
    dplyr::arrange(model_index)
  if (!identical(as.integer(summary$model_index), c(2L, 26L))) stop("Incomplete measurement pilot.")
  timings <- purrr::map_dfr(sub("\\.qs$", "-timing.csv", paths),
    ~ readr::read_csv(.x, show_col_types = FALSE))
  if (any(timings$variant != variant) || any(!is.finite(timings$elapsed_seconds))) {
    stop("Measurement pilot timings are unavailable.")
  }
  list(variant = variant, summary = dplyr::left_join(summary, timings,
    by = "model_index"), all_pass = all(summary$diagnostic_pass))
}

write_measurement_pilot_audit <- function(audit) {
  path <- measurement_path(audit$variant, "performance_pilot", "audit.csv", FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(audit$summary, path)
  path
}

make_measurement_remedy_records <- function(specifications, decisions, variant,
                                           initial_diagnostics) {
  if (!identical(initial_diagnostics$variant, variant) ||
      nrow(initial_diagnostics$diagnostic_summary) != 26L) {
    stop("Inspect a complete initial diagnostic audit before specifying remedies.")
  }
  if (any(!decisions$variant %in% measurement_variants())) stop("Unknown remedy measurement set.")
  if (!setequal(names(decisions), c("variant", "model_index", "stage",
      "adapt_delta", "max_treedepth", "rationale"))) stop("Unrecognised remedy decision columns.")
  decisions <- decisions[decisions$variant == variant, , drop = FALSE]
  if (!nrow(decisions)) return(list(list(skip = TRUE, measurement_variant = variant)))
  required <- c("model_index", "stage", "adapt_delta", "max_treedepth", "rationale")
  if (!all(required %in% names(decisions)) || anyNA(decisions[required]) ||
      any(!nzchar(trimws(decisions$rationale))) ||
      any(!decisions$model_index %in% seq_len(26L)) ||
      any(!grepl("^remedy_[1-9][0-9]*$", decisions$stage)) ||
      anyDuplicated(decisions[c("model_index", "stage")])) stop("Invalid measurement remedy decisions.")
  purrr::pmap(decisions, function(variant, model_index, stage, adapt_delta,
                                 max_treedepth, rationale, ...) {
    specification <- specifications[[as.integer(model_index)]]
    if (adapt_delta < specification$adapt_delta || adapt_delta >= 1 ||
        max_treedepth < specification$max_treedepth || max_treedepth %% 1 != 0) {
      stop("A measurement remedy must not weaken sampler controls.")
    }
    specification$measurement_stage <- stage
    specification$adapt_delta <- adapt_delta
    specification$max_treedepth <- as.integer(max_treedepth)
    specification$measurement_remedy_reason <- rationale
    specification
  })
}

fit_measurement_remedy <- function(record, selection, priors, preflight) {
  if (isTRUE(record$skip)) {
    if (!measurement_phase_enabled("remedy")) stop("Measurement remedy phase is disabled.")
    return(write_measurement_qs(list(skip = TRUE), measurement_path(
      selection$variant, "remedies", "no_remedies.qs"
    )))
  }
  fit_measurement_candidate(record, selection, priors, preflight)
}

measurement_diagnostic_audit <- function(paths, variant) {
  diagnostics <- lapply(measurement_fit_paths(paths), function(path) {
    result <- diagnose_update_2026_primary_fit_file(path)
    gc(verbose = FALSE)
    result
  })
  summary <- dplyr::bind_rows(purrr::map(diagnostics, "summary"))
  summary <- summary[order(summary$model_index), , drop = FALSE]
  if (!identical(as.integer(summary$model_index), seq_len(26L)) ||
      any(summary$outcome != "hypertrophy")) stop("Incomplete measurement diagnostic set.")
  summary$ppc_warning <- primary_ppc_warning(summary)
  summary$status <- ifelse(!summary$diagnostic_pass, "REQUIRES REMEDIATION",
                          ifelse(summary$ppc_warning, "PASS WITH PPC WARNING", "PASS"))
  ppc <- summarise_primary_ppc_outcome(diagnostics, "hypertrophy")
  list(variant = variant, diagnostic_summary = summary,
       parameter_failures = purrr::map_dfr(diagnostics, function(x) {
         dplyr::mutate(x$parameter_failures, model_index = x$summary$model_index[[1L]])
       }), ppc_observation_flags = ppc$flags, ppc_outcome_summary = ppc$summary)
}

write_measurement_diagnostics <- function(audit, stage) {
  tables <- audit[setdiff(names(audit), "variant")]
  directory <- dirname(measurement_path(audit$variant,
    paste0("diagnostics/", stage), "x", FALSE))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  purrr::imap_chr(tables, function(table, name) {
    path <- file.path(directory, paste0(name, ".csv"))
    readr::write_csv(table, path)
    path
  }) |> unname()
}

measurement_data_fingerprint <- function(data) {
  columns <- c("study_no", "arm_no", "outcome_no", "time1", "time2", "yi", "vi",
    "protein_intake", "protein_intake_raw", "protein_median", "resistance_exercise_code")
  if (!all(columns %in% names(data))) stop("Fit data lack measurement identity fields.")
  ids <- c("study_no", "arm_no", "outcome_no")
  values <- lapply(columns, function(column) {
    x <- data[[column]]
    if (is.factor(x)) x <- as.character(x)
    unname(if (column %in% ids) as.character(x) else as.numeric(x))
  })
  names(values) <- columns
  digest::digest(values, algo = "sha256")
}

validate_measurement_fit_data <- function(candidate, selection) {
  # brms retains a model frame, not the complete extraction table. Reconstruct
  # its predictors from the frozen selection and verify every retained column.
  expected <- selection$data
  specification <- candidate$specification
  if (isTRUE(specification$standardised_moderators)) {
    expected <- add_primary_standardised_moderators(expected)
  }
  if (isTRUE(specification$half_unit_protein)) {
    expected <- add_primary_half_unit_protein(expected)
  }
  actual <- candidate$model$data
  required <- c("yi", "vi", "arm_no", "outcome_no", "resistance_exercise_code")
  if (!all(required %in% names(actual)) ||
      !any(c("study_no", "study_factor") %in% names(actual)) ||
      nrow(actual) != nrow(expected)) {
    stop("Archived measurement fit lacks row identity fields or has a different row count.")
  }
  expected$Intercept <- rep(1, nrow(expected))
  for (column in names(actual)) {
    if (!column %in% names(expected)) {
      # brms also retains evaluated piecewise-linear basis matrices.
      if (!grepl("^lspline\\(protein_intake, [0-9.]+\\)$", column)) {
        stop("Unrecognised archived measurement predictor: ", column)
      }
      expected[[column]] <- eval(parse(text = column),
        envir = as.list(expected), enclos = asNamespace("lspline"))
    }
    normalise <- function(x) {
      if (is.factor(x) || is.character(x)) return(as.character(x))
      if (is.matrix(x)) return(unname(x))
      unname(as.numeric(x))
    }
    if (!isTRUE(all.equal(normalise(actual[[column]]),
                          normalise(expected[[column]]), tolerance = 1e-12))) {
      stop("Archived measurement fit does not match frozen rows/predictor: ", column)
    }
  }
  invisible(TRUE)
}

measurement_accepted_manifest <- function(initial_paths, remedy_paths, variant, selection) {
  expected_fingerprint <- measurement_data_fingerprint(selection$data)
  describe <- function(path) {
    candidate <- qs::qread(path)
    if (isTRUE(candidate$skip)) return(NULL)
    s <- candidate$specification
    if (!identical(s$measurement_variant, variant) || candidate$outcome != "hypertrophy") {
      stop("A fit belongs to a different measurement set.")
    }
    validate_measurement_fit_data(candidate, selection)
    row <- tibble::tibble(outcome = "hypertrophy", model_index = as.integer(candidate$model_index),
      model_name = candidate$model_name, stage = s$measurement_stage,
      fit_path = path, bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256"),
      data_sha256 = expected_fingerprint,
      rank = if (s$measurement_stage %in% c("initial", "pilot")) 0L else
        as.integer(sub("remedy_", "", s$measurement_stage)))
    rm(candidate)
    gc(verbose = FALSE)
    row
  }
  manifest <- purrr::map_dfr(c(measurement_fit_paths(initial_paths),
    measurement_fit_paths(remedy_paths)), describe)
  manifest <- manifest |>
    dplyr::group_by(model_index) |> dplyr::slice_max(rank, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup() |> dplyr::arrange(model_index)
  if (!identical(manifest$model_index, seq_len(26L))) stop("Incomplete accepted measurement fits.")
  manifest
}

measurement_loo_records <- function(manifest, diagnostics, variant) {
  if (!all(diagnostics$diagnostic_summary$diagnostic_pass)) {
    stop("Measurement diagnostics require remediation; no candidate may be dropped before LOO.")
  }
  purrr::pmap(manifest, function(...) c(list(...), list(variant = variant)))
}

compute_measurement_integrated_loo <- function(record) {
  if (!measurement_phase_enabled("loo")) stop("Measurement LOO is disabled.")
  if (!identical(digest::digest(file = record$fit_path, algo = "sha256"), record$sha256)) {
    stop("Accepted measurement fit checksum changed before LOO.")
  }
  candidate <- qs::qread(record$fit_path)
  if (!identical(candidate$specification$measurement_variant, record$variant) ||
      candidate$model_index != record$model_index) stop("LOO fit identity mismatch.")
  initial <- compute_primary_candidate_psis_loo(candidate, cores = 1L, save_psis = FALSE)
  result <- compute_primary_candidate_integrated_psis_loo(
    candidate, initial, k_threshold = 0.7, cores = 1L
  )
  result$fit_path <- record$fit_path
  result$fit_stage <- record$stage
  result$measurement_variant <- record$variant
  write_measurement_qs(result, measurement_path(record$variant, "loo/integrated",
    sprintf("model_%02d.qs", record$model_index)))
}

measurement_integrated_manifest <- function(paths, variant) {
  manifest <- purrr::map_dfr(unname(unlist(paths)), function(path) {
    result <- qs::qread(path)
    if (!inherits(result$loo, "loo") || !identical(result$measurement_variant, variant)) {
      stop("Invalid measurement integrated LOO reduction.")
    }
    if (anyNA(loo::pareto_k_values(result$loo)) ||
        nrow(result$loo$pointwise) != nrow(result$observation_manifest)) {
      stop("Integrated LOO must provide diagnostics for every retained observation.")
    }
    tibble::tibble(outcome = result$outcome, model_index = as.integer(result$model_index),
      model_name = result$model_name, fit_path = result$fit_path, fit_stage = result$fit_stage,
      integrated_loo_path = path, bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256"))
  }) |> dplyr::arrange(model_index)
  if (!identical(manifest$model_index, seq_len(26L))) stop("Incomplete measurement LOO reductions.")
  manifest
}

write_measurement_loo_audit <- function(manifest, variant) {
  results <- lapply(manifest$integrated_loo_path, qs::qread)
  tables <- list(reduction_manifest = manifest,
    method_comparison = purrr::map_dfr(results, "comparison"),
    observation_diagnostics = purrr::map_dfr(results, "observation_diagnostics"))
  tables$high_k_cases <- tables$observation_diagnostics |>
    dplyr::filter(integrated_pareto_k > 0.7)
  directory <- dirname(measurement_path(variant, "loo_psis", "x", FALSE))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  purrr::imap_chr(tables, function(table, name) {
    path <- file.path(directory, paste0(name, ".csv.gz"))
    write_verification_csv(table, path)
    path
  }) |> unname()
}

measurement_exact_records <- function(manifest, variant) {
  records <- list()
  for (i in seq_len(nrow(manifest))) {
    result <- qs::qread(manifest$integrated_loo_path[[i]])
    cases <- which(loo::pareto_k_values(result$loo) > 0.7)
    for (observation in cases) {
      records[[length(records) + 1L]] <- c(as.list(manifest[i, , drop = FALSE]),
        list(observation = as.integer(observation), case_id = length(records) + 1L, variant = variant))
    }
  }
  if (!length(records)) list(list(skip = TRUE, variant = variant)) else records
}

fit_measurement_exact <- function(record) {
  if (!measurement_phase_enabled("exact")) stop("Measurement exact-LOO refitting is disabled.")
  if (isTRUE(record$skip)) return(write_measurement_qs(list(skip = TRUE),
    measurement_path(record$variant, "loo/exact", "no_exact_refits.qs")))
  candidate <- qs::qread(record$fit_path)
  integrated <- qs::qread(record$integrated_loo_path)
  if (!identical(candidate$specification$measurement_variant, record$variant) ||
      !identical(integrated$measurement_variant, record$variant) ||
      candidate$model_index != record$model_index ||
      !isTRUE(loo::pareto_k_values(integrated$loo)[[record$observation]] > 0.7)) {
    stop("Invalid measurement exact-LOO case.")
  }
  job <- c(record, list(fitted_candidate = candidate, integrated_result = integrated))
  refit <- fit_primary_integrated_exact_loo_case(job)
  # Use the same general primary design reconstruction as achieved/assignment
  # exact refits, including random slopes and study factor smooths.
  prediction_model <- refit$refitted_model
  prediction_model$data <- refit$original_data
  components <- get_primary_integrated_loo_components(list(
    outcome = refit$outcome, model_index = refit$model_index,
    model_name = refit$model_name, model = prediction_model
  ))
  correction <- compute_update_2026_primary_exact_loo_case(refit, components = components)
  correction$measurement_variant <- record$variant
  write_measurement_qs(correction, measurement_path(record$variant, "loo/exact",
    sprintf("model_%02d_observation_%04d.qs", record$model_index, record$observation)))
}

finalise_measurement_loo <- function(manifest, exact_paths, variant) {
  corrections <- lapply(unname(unlist(exact_paths)), qs::qread)
  corrections <- Filter(function(x) !isTRUE(x$skip), corrections)
  exact_summary <- if (length(corrections)) {
    dplyr::bind_rows(purrr::map(corrections, "exact_refit_summary"))
  } else tibble::tibble(model_index = integer(), observation = integer(), diagnostic_pass = logical())
  if (length(corrections) && (!all(exact_summary$diagnostic_pass) ||
      any(!is.finite(exact_summary$exact_elpd_loo)))) {
    stop("Measurement exact-LOO diagnostics failed; remediate before stacking.")
  }
  results <- lapply(manifest$integrated_loo_path, qs::qread)
  expected <- purrr::map_dfr(results, function(x) tibble::tibble(
    model_index = x$model_index, observation = which(loo::pareto_k_values(x$loo) > 0.7)))
  actual <- purrr::map_dfr(corrections, function(x) tibble::tibble(
    model_index = x$model_index, observation = x$exact_refit_summary$observation[[1L]]))
  if (!length(corrections)) actual <- expected[FALSE, , drop = FALSE]
  order_cases <- function(x) x[order(x$model_index, x$observation), , drop = FALSE]
  if (!isTRUE(all.equal(order_cases(expected), order_cases(actual), check.attributes = FALSE))) {
    stop("Exact corrections do not cover precisely the new high-k cases.")
  }
  paths <- character(length(results))
  for (i in seq_along(results)) {
    result <- results[[i]]
    selected <- Filter(function(x) x$model_index == result$model_index, corrections)
    for (correction in selected) {
      if (!identical(correction$measurement_variant, variant) ||
          !identical(correction$observation_manifest, result$observation_manifest)) {
        stop("Exact-LOO observation identity mismatch.")
      }
      obs <- correction$exact_refit_summary$observation[[1L]]
      result$loo$pointwise[obs, ] <- correction$loo$pointwise[obs, ]
      result$loo$diagnostics$pareto_k[[obs]] <- NA_real_
      if (!is.null(result$loo$diagnostics$n_eff)) result$loo$diagnostics$n_eff[[obs]] <- NA_real_
      result$observation_diagnostics$integrated_pareto_k[[obs]] <- NA_real_
    }
    result$loo <- update_primary_exact_loo_estimates(result$loo)
    result$exact_refit_summary <- exact_summary[exact_summary$model_index == result$model_index, , drop = FALSE]
    result$observation_diagnostics$exact_loo_refit <-
      result$observation_diagnostics$observation %in% result$exact_refit_summary$observation
    paths[[i]] <- write_measurement_qs(result, measurement_path(variant, "loo/final",
      sprintf("model_%02d.qs", result$model_index)))
  }
  paths
}

measurement_stacking <- function(paths, variant) {
  results <- lapply(unname(unlist(paths)), qs::qread)
  indices <- vapply(results, function(x) as.integer(x$model_index), integer(1))
  if (!identical(indices, seq_len(26L))) stop("Measurement final LOO order is invalid.")
  first_manifest <- results[[1L]]$observation_manifest
  for (result in results) {
    exact_observations <- result$exact_refit_summary$observation
    missing_k <- which(is.na(loo::pareto_k_values(result$loo)))
    if (!identical(result$measurement_variant, variant) ||
        !identical(result$observation_manifest, first_manifest) ||
        !setequal(missing_k, exact_observations) ||
        !all(result$exact_refit_summary$diagnostic_pass) ||
        any(loo::pareto_k_values(result$loo) > 0.7, na.rm = TRUE)) {
      stop("Measurement final LOO set is not stacking-ready.")
    }
  }
  objects <- lapply(results, `[[`, "loo")
  names(objects) <- vapply(results, `[[`, character(1), "model_name")
  weights <- as.numeric(loo::loo_model_weights(objects, method = "stacking"))
  if (any(!is.finite(weights)) || any(weights < -1e-8) || abs(sum(weights) - 1) > 1e-8) {
    stop("Invalid measurement stacking weights.")
  }
  weights <- pmax(weights, 0); weights <- weights / sum(weights)
  tibble::tibble(outcome = "hypertrophy", model_index = seq_len(26L),
    model_name = names(objects), stacking_weight = weights)
}

write_measurement_model_average <- function(manifest, weights, selection) {
  if (!measurement_phase_enabled("report")) stop("Measurement reporting is disabled.")
  models <- dplyr::left_join(manifest, weights,
    by = c("outcome", "model_index", "model_name")) |> dplyr::arrange(model_index)
  if (anyNA(models$stacking_weight)) stop("Missing measurement stacking weights.")
  record <- list(outcome = "hypertrophy", seed = 1990L, models = models)
  build_update_2026_primary_model_average_to_file(
    record, list(hypertrophy = selection$data), output_file = measurement_path(
      selection$variant, "reporting", "model_average.qs"
    )
  )
}

write_measurement_report <- function(average_path, selection) {
  if (!measurement_phase_enabled("report")) stop("Measurement reporting is disabled.")
  make_update_2026_primary_reporting_file(average_path, list(hypertrophy = selection$data),
    output_file = measurement_path(selection$variant, "reporting", "report.qs"))
}

write_measurement_checkpoint <- function(report_path, selection, specifications,
                                        diagnostics, manifest, loo_paths, checkpoint) {
  report <- qs::qread(report_path)
  if (nrow(report$model_weights) != 26L ||
      abs(sum(report$model_weights$stacking_weight) - 1) > 1e-8 ||
      any(!is.finite(report$dose_response$estimate))) stop("Invalid measurement report.")
  reductions <- lapply(unname(unlist(loo_paths)), qs::qread)
  accepted_specifications <- purrr::map_dfr(manifest$fit_path, function(path) {
    candidate <- qs::qread(path)
    specification <- tibble::as_tibble(candidate$specification)
    if (nrow(specification) != 1L) stop("Accepted measurement specification is not a single record.")
    rm(candidate)
    gc(verbose = FALSE)
    specification
  })
  artifacts <- list(
    source_manifest = checkpoint$source_manifest, selection_audit = selection$selection,
    measurement_annotations = selection$annotations,
    selection_summary = selection$summary,
    analysis_data = dplyr::select(selection$data, -dplyr::any_of(c(
      "corresponding_author", "corresponding_author_email"))),
    initial_model_specifications = dplyr::bind_rows(specifications),
    model_specifications = accepted_specifications,
    accepted_fit_manifest = manifest,
    diagnostic_summary = diagnostics$diagnostic_summary,
    diagnostic_parameter_failures = diagnostics$parameter_failures,
    posterior_predictive_observation_flags = diagnostics$ppc_observation_flags,
    posterior_predictive_outcome_summary = diagnostics$ppc_outcome_summary,
    loo_model_summary = purrr::map_dfr(reductions, function(x) dplyr::mutate(
      summarise_primary_loo_method(x$loo, "integrated_observation_level_psis_plus_exact"),
      outcome = x$outcome, model_index = x$model_index, model_name = x$model_name)),
    loo_exact_refit_summary = dplyr::bind_rows(purrr::map(reductions, "exact_refit_summary")),
    loo_observation_diagnostics = purrr::map_dfr(reductions, "observation_diagnostics"),
    stacking_weights = report$model_weights, mixture_draw_allocation = report$draw_allocation,
    dose_response_predictions = report$dose_response,
    pairwise_contrasts = report$pairwise_contrasts, observed_support = report$support,
    local_equivalence = report$local_equivalence
  )
  artifacts$reporting_settings <- tibble::as_tibble(c(report$settings,
    list(sesoi = 0.1, probability_criterion = 0.95, reference_age = 25, reference_weeks = 12)))
  artifacts$final_loo_manifest <- tibble::tibble(model_index = seq_len(26L),
    final_loo_path = unname(unlist(loo_paths)),
    bytes = file.info(unname(unlist(loo_paths)))$size,
    sha256 = vapply(unname(unlist(loo_paths)), function(path) {
      digest::digest(file = path, algo = "sha256")
    }, character(1)))
  artifacts$verification_checks <- tibble::tibble(
    check = c("26 accepted candidates", "sampling diagnostics pass", "stacking sums to one",
      "all new high-k cases corrected", "reported contrasts preserve support"),
    passed = c(nrow(manifest) == 26L, all(diagnostics$diagnostic_summary$diagnostic_pass),
      abs(sum(report$model_weights$stacking_weight) - 1) < 1e-8,
      all(vapply(reductions, function(x) !any(loo::pareto_k_values(x$loo) > 0.7, na.rm = TRUE), logical(1))),
      all(!is.na(report$pairwise_contrasts$supported))))
  if (!all(artifacts$verification_checks$passed)) stop("Measurement checkpoint verification failed.")
  descriptions <- c(
    source_manifest = "Checksums of frozen primary data and candidate specifications.",
    selection_audit = "Every original muscle row with retained/excluded decision and reason.",
    measurement_annotations = "Author-confirmed measurement interpretations, row identifiers, dates, and provenance; frozen source labels unchanged.",
    selection_summary = "Study, arm, and contrast counts after measurement restriction.",
    analysis_data = "Retained analysis rows; no participant-level data or contact fields.",
    initial_model_specifications = "Inherited final-primary specifications before new diagnostics.",
    model_specifications = "Specifications actually used in accepted fits, including new remedies.",
    accepted_fit_manifest = "Single-copy accepted fit paths, stages, byte sizes, and SHA-256 hashes.",
    diagnostic_summary = "Unchanged primary sampling/convergence and PPC diagnostics.",
    diagnostic_parameter_failures = "Parameters that failed the primary diagnostic thresholds.",
    posterior_predictive_observation_flags = "Observation-level PPC warnings across the ensemble.",
    posterior_predictive_outcome_summary = "Outcome-wide PPC summary for the retained measurements.",
    loo_model_summary = "Final integrated observation-level LOO after exact corrections.",
    loo_exact_refit_summary = "All newly required case-deletion refits and their diagnostic results.",
    loo_observation_diagnostics = "Observation diagnostics identifying exact-corrected cases.",
    stacking_weights = "New 26-model stacking weights for this restricted dataset.",
    mixture_draw_allocation = "Deterministic posterior-mixture model draw counts.",
    dose_response_predictions = "Population dose-response predictions at age 25 and week 12.",
    pairwise_contrasts = "Pairwise SESOI probabilities, intervals, and empirical support flags.",
    observed_support = "RT-specific observed intake ranges.",
    local_equivalence = "Outcome-specific local practical-equivalence summaries.",
    reporting_settings = "Prediction grid, draws, seed, reference profile, and decision thresholds.",
    final_loo_manifest = "Final LOO reduction paths, byte sizes, and checksums.",
    verification_checks = "Required final reporting and diagnostic invariants."
  )
  paths <- write_verification_bundle(paste0("measurement_", selection$variant, "_update_2026"),
    artifacts, descriptions)
  bundle_manifest <- utils::read.csv(paths[basename(paths) == "manifest.csv"], stringsAsFactors = FALSE)
  directory <- dirname(paths[[1L]])
  hashes <- vapply(file.path(directory, bundle_manifest$file), function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1))
  if (!identical(unname(hashes), bundle_manifest$sha256)) stop("Measurement bundle checksum audit failed.")
  paths
}

write_measurement_figure <- function(report_path, selection) {
  report <- qs::qread(report_path)
  heatmap <- list(contrasts = report$pairwise_contrasts, support = report$support,
    settings = c(report$settings, list(sesoi = 0.1, conf_level = 0.95)))
  class(heatmap) <- c("primary_heatmap_contrasts", class(heatmap))
  title <- if (selection$variant == "direct_only") {
    "changes in muscle size"
  } else "changes in measures of muscle size or mass"
  figure <- plot_primary_predictions_heatmap(selection$data, report$dose_response,
    heatmap, title = title)
  caption <- if (selection$variant == "direct_only") {
    "Restricted to MRI, CT, ultrasound, and muscle biopsy; whole-muscle and fibre measures differ."
  } else {
    "Trunk and skinfold estimates excluded; limb DXA/BIA proxies retained alongside imaging and biopsy."
  }
  figure <- figure + patchwork::plot_annotation(caption = caption)
  directory <- here::here("plots", "update_2026", "measurement", selection$variant)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(directory, c("muscle-manuscript.pdf", "muscle.tiff"))
  ggplot2::ggsave(paths[[1L]], figure, device = grDevices::cairo_pdf, width = 8, height = 12)
  ggplot2::ggsave(paths[[2L]], figure, device = "tiff", dpi = 300, width = 8, height = 12)
  paths
}
