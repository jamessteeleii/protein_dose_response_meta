# Restore missing parent fits for the interaction workflow ----

# Use the recovery ledger to identify the missing primary candidates and the
# reviewed settings needed to reconstruct them.

# Recovery of the already-pruned total-intake primary family. Historical
# results are immutable; the replacement posterior/LOO files have new paths.
iad_recovery_specifications <- function(checkpoints, ledger_path) {
  checkpoint <- checkpoints$total
  specs <- checkpoint$specifications
  diagnostics <- checkpoint$diagnostics
  decisions <- iad_read_csv(ledger_path)
  if (nrow(specs) != 72L || nrow(diagnostics) != 72L) stop("Incomplete frozen total-primary recipe.")

  records <- lapply(seq_len(nrow(specs)), function(i) {
    s <- as.list(specs[i, ])
    d <- diagnostics[diagnostics$outcome == s$outcome & diagnostics$model_index == s$model_index, ]
    if (nrow(d) != 1L || !d$diagnostic_pass) stop("Missing accepted parent diagnostic recipe.")
    s$adapt_delta <- d$adapt_delta[[1L]]; s$max_treedepth <- as.integer(d$max_treedepth[[1L]])
    s$family <- "total"; s$parent_index <- as.integer(s$model_index)
    s$model_index <- as.integer(s$model_index); s$interaction <- "none"
    s$stage <- "recovery_1"
    s$accepted_fit_path <- s$accepted_fit_stage <- NULL
    override <- decisions[decisions$outcome == s$outcome & decisions$model_index == s$model_index, ]
    if (nrow(override) > 1L) stop("Duplicate recovery schedule override.")
    if (nrow(override)) {
      s$iter_warmup <- as.integer(override$iter_warmup)
      s$iter_sampling <- as.integer(override$iter_sampling)
      s$recovery_rationale <- override$rationale
    }
    if (s$chains != 4L || s$iter_warmup != 2000L || !s$iter_sampling %in% c(6000L, 12000L)) {
      stop("Unexpected recovery schedule; review the frozen recipe.")
    }
    s
  })

  if (nrow(decisions) != 1L || decisions$outcome != "hypertrophy" || decisions$model_index != 11L ||
      decisions$iter_sampling != 12000L) stop("Recovery override must preserve the documented model-11 sampling remedy.")
  records
}

iad_recovery_map_input <- function() {
  path <- iad_path("recovery", "accepted_manifest.csv")
  if (file.exists(path)) path else here::here("sensitivity", "age-duration-recovery-status.csv")
}

iad_apply_recovery <- function(checkpoints, mapping_path) {
  mapping <- iad_read_csv(mapping_path)
  if (!nrow(mapping)) return(checkpoints)
  parents <- checkpoints$total$parents
  key <- function(x) paste(x$outcome, x$model_index, sep = "|")
  if (nrow(mapping) != 72L || anyDuplicated(key(mapping)) || !setequal(key(mapping), key(parents)) ||
      !all(mapping$diagnostic_pass) || any(!file.exists(mapping$fit_path))) stop("Incomplete recovery mapping.")

  source_hash <- iad_hash(here::here("results", "verification", "total_primary_update_2026", "manifest.csv"))
  if (any(mapping$source_manifest_sha256 != source_hash)) stop("Recovery mapping belongs to different frozen inputs.")
  mapping <- mapping[match(key(parents), key(mapping)), ]
  for (i in seq_len(nrow(mapping))) {
    iad_check_hash(mapping$loo_path[[i]], mapping$loo_sha256[[i]])
    result <- qs::qread(mapping$loo_path[[i]])
    iad_validate_loo(result, checkpoints$total$data[[mapping$outcome[[i]]]])
    if (mapping$data_sha256[[i]] != iad_fingerprint(checkpoints$total$data[[mapping$outcome[[i]]]])) {
      stop("Recovery predictor identity differs from the frozen common sample.")
    }
  }

  parents$fit_path <- mapping$fit_path; parents$stage <- mapping$stage
  parents$loo_path <- mapping$loo_path; parents$loo_sha256 <- mapping$loo_sha256
  parents$available <- TRUE; parents$expected_bytes <- mapping$bytes
  checkpoints$total$parents <- parents
  checkpoints$total$counts$available_parent_fits <- 18L
  checkpoints
}

iad_recovery_pilot_specs <- function(specifications) {
  Filter(function(s) s$model_index %in% c(2L, 26L), specifications)
}

iad_recovery_fit_path <- function(s) {
  if (identical(s$stage, "recovery_1")) {
    iad_path("recovery", s$outcome, sprintf("model_%02d.qs", s$model_index), archive = TRUE)
  } else {
    if (!grepl("^remedy_[1-9][0-9]*$", s$stage)) stop("Invalid recovery archive stage.")
    iad_path("recovery", s$stage, s$outcome, sprintf("model_%02d.qs", s$model_index), archive = TRUE)
  }
}

iad_recovery_fit_remedy <- function(record, checkpoints, priors, preflight) {
  iad_require_phase("recover_remedy")
  if (isTRUE(record$skip)) return(iad_qsave(record,
    iad_path("recovery", "remedies", "none.qs", archive = TRUE)))
  iad_recovery_fit(record, checkpoints, priors, preflight, "recover_remedy")
}

# Fit a recovery candidate from the frozen parent inputs and write the
# artifact for the acceptance checks.
iad_recovery_fit <- function(s, checkpoints, priors, preflight, phase, pilot = NULL) {
  iad_require_phase(phase)
  if (nrow(preflight) != 72L || !all(preflight$valid)) stop("Recovery input gate failed.")
  if (phase == "recover_fit" && (nrow(pilot) != 8L || !all(pilot$diagnostic_pass))) stop("Recovery pilot gate failed.")
  input <- iad_prepare_input(s, checkpoints, priors)
  path <- iad_recovery_fit_path(s)

  if (file.exists(path)) {
    previous <- qs::qread(path)
    if (!identical(previous$specification, s) || previous$data_sha256 != iad_fingerprint(input$data)) {
      stop("Recovery file differs; retain it and use a reviewed remedy stage.")
    }
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  started <- proc.time()[["elapsed"]]
  model <- brms::brm(stats::as.formula(input$formula), data = input$data, prior = input$priors,
    chains = 4L, cores = 4L, warmup = s$iter_warmup, iter = s$iter_warmup + s$iter_sampling,
    seed = 831988L + 100L * match(s$outcome, iad_outcomes()) + s$model_index,
    control = list(adapt_delta = s$adapt_delta, max_treedepth = s$max_treedepth),
    save_pars = brms::save_pars(all = TRUE))

  for (column in setdiff(names(input$data), names(model$data))) model$data[[column]] <- input$data[[column]]
  iad_qsave(list(outcome = s$outcome, model_index = s$model_index, model_name = s$model_name,
    model = model, specification = s, data_sha256 = iad_fingerprint(input$data),
    elapsed_seconds = proc.time()[["elapsed"]] - started), path)
}

iad_accept_recovery <- function(records, final_paths, checkpoints) {
  iad_require_phase("recover_accept")
  rows <- dplyr::bind_rows(records)
  if (nrow(rows) != 72L || !all(rows$diagnostic_pass)) stop("Recovery diagnostics incomplete.")
  final_paths <- unname(unlist(final_paths))
  if (length(final_paths) != 72L) stop("Recovery LOO incomplete.")
  results <- lapply(final_paths, qs::qread)
  loo <- purrr::map2_dfr(results, final_paths, function(r, path) {
    iad_validate_loo(r, checkpoints$total$data[[r$outcome]])
    tibble::tibble(outcome = r$outcome, model_index = r$model_index,
      loo_path = path, loo_sha256 = iad_hash(path))
  })
  rows <- dplyr::left_join(rows, loo, by = c("outcome", "model_index"))
  rows$bytes <- file.info(rows$fit_path)$size
  rows$source_manifest_sha256 <- iad_hash(here::here("results", "verification", "total_primary_update_2026", "manifest.csv"))
  if (anyNA(rows$loo_path) || anyDuplicated(rows[c("outcome", "model_index")])) stop("Recovery identity mismatch.")
  iad_csv(rows, iad_path("recovery", "accepted_manifest.csv"))
}

# Recovery pipeline and acceptance ----

# Connect the recovery pilot, fits, diagnostics and accepted-fit mapping
# consumed by the interaction pipeline.
iad_recovery_targets <- function() {
  target <- targets::tar_target_raw
  add <- function(name, command, pattern = NULL, phase = NULL, format = "qs", iteration = "vector") {
    args <- list(name = paste0("iad_recovery_", name), command = command, format = format, iteration = iteration)
    if (!is.null(pattern)) args$pattern <- pattern
    if (!is.null(phase)) { args$cue <- iad_cue(phase); args$error <- "continue" }
    do.call(target, args, quote = TRUE)
  }
  list(
    add("schedule_file", quote(here::here("sensitivity", "age-duration-total-recovery.csv")), format = "file"),
    add("specification", quote(iad_recovery_specifications(iad_frozen_checkpoints, iad_recovery_schedule_file)), iteration = "list"),
    add("preflight", quote(iad_preflight(iad_recovery_specification, iad_frozen_checkpoints, iad_base_priors))),
    add("preparation_files", quote(c(iad_csv(dplyr::bind_rows(iad_recovery_specification),
      iad_path("recovery", "specifications.csv")), iad_csv(iad_recovery_preflight,
      iad_path("recovery", "preflight.csv")))), format = "file"),
    add("pilot_specification", quote(iad_recovery_pilot_specs(iad_recovery_specification)), iteration = "list"),
    add("remaining_specification", quote(Filter(function(s) !s$model_index %in% c(2L, 26L),
      iad_recovery_specification)), iteration = "list"),
    add("pilot_file", quote(iad_recovery_fit(iad_recovery_pilot_specification, iad_frozen_checkpoints,
      iad_base_priors, iad_recovery_preflight, "recover_pilot")), pattern = quote(map(iad_recovery_pilot_specification)),
      phase = "recover_pilot", format = "file", iteration = "list"),
    add("pilot_diagnostic", quote(iad_diagnose(iad_recovery_pilot_file)),
      pattern = quote(map(iad_recovery_pilot_file)), iteration = "list"),
    add("pilot_summary", quote(iad_diagnostic_summary(iad_recovery_pilot_diagnostic))),
    add("pilot_audit_files", quote(iad_diagnostic_files(iad_recovery_pilot_diagnostic, "recovery_pilot")), format = "file"),
    add("fit_file", quote(iad_recovery_fit(iad_recovery_remaining_specification, iad_frozen_checkpoints,
      iad_base_priors, iad_recovery_preflight, "recover_fit", iad_recovery_pilot_summary)),
      pattern = quote(map(iad_recovery_remaining_specification)), phase = "recover_fit", format = "file", iteration = "list"),
    add("fit_diagnostic", quote(iad_diagnose(iad_recovery_fit_file)), pattern = quote(map(iad_recovery_fit_file)), iteration = "list"),
    add("audit_files", quote(iad_diagnostic_files(c(iad_recovery_pilot_diagnostic,
      iad_recovery_fit_diagnostic), "recovery_initial")), format = "file"),
    add("initial_summary", quote(iad_diagnostic_summary(c(iad_recovery_pilot_diagnostic,
      iad_recovery_fit_diagnostic)))),
    add("remedy_ledger_file", quote(here::here("sensitivity", "age-duration-recovery-remedies.csv")), format = "file"),
    add("remedy_specification", quote(iad_remedy_specs(iad_recovery_specification,
      iad_recovery_remedy_ledger_file, iad_recovery_initial_summary)), iteration = "list"),
    add("remedy_fit_file", quote(iad_recovery_fit_remedy(iad_recovery_remedy_specification,
      iad_frozen_checkpoints, iad_base_priors, iad_recovery_preflight)),
      pattern = quote(map(iad_recovery_remedy_specification)), phase = "recover_remedy", format = "file", iteration = "list"),
    add("remedy_diagnostic", quote(iad_diagnose_remedy(iad_recovery_remedy_fit_file)),
      pattern = quote(map(iad_recovery_remedy_fit_file)), iteration = "list"),
    add("remedy_audit_files", quote(iad_diagnostic_files(iad_recovery_remedy_diagnostic,
      "recovery_remedy")), format = "file"),
    add("accepted_record", quote(iad_accepted(iad_recovery_initial_summary,
      iad_diagnostic_summary(iad_recovery_remedy_diagnostic), iad_recovery_specification, iad_frozen_checkpoints)), iteration = "list"),
    add("accepted_audit_file", quote(iad_csv(dplyr::bind_rows(iad_recovery_accepted_record),
      iad_path("diagnostics", "recovery_accepted.csv"))), format = "file"),
    add("loo_file", quote(iad_compute_loo(iad_recovery_accepted_record, "recover_loo")),
      pattern = quote(map(iad_recovery_accepted_record)), phase = "recover_loo", format = "file", iteration = "list"),
    add("loo_audit", quote(iad_loo_audit(iad_recovery_loo_file))),
    add("loo_audit_files", quote(c(iad_csv(iad_recovery_loo_audit$summary, iad_path("recovery", "loo_summary.csv")),
      iad_csv(dplyr::bind_rows(iad_recovery_loo_audit$cases), iad_path("recovery", "exact_cases.csv")))), format = "file"),
    add("exact_record", quote(iad_recovery_loo_audit$cases), iteration = "list"),
    add("exact_file", quote(iad_exact(iad_recovery_exact_record, "recover_exact")),
      pattern = quote(map(iad_recovery_exact_record)), phase = "recover_exact", format = "file", iteration = "list"),
    add("exact_diagnostic", quote(iad_exact_diagnostics(iad_recovery_exact_file)),
      pattern = quote(map(iad_recovery_exact_file)), iteration = "list"),
    add("exact_audit_file", quote(iad_csv(dplyr::bind_rows(iad_recovery_exact_diagnostic),
      iad_path("recovery", "exact_diagnostics.csv"))), format = "file"),
    add("correction_file", quote(iad_score_exact(iad_recovery_exact_file, "recover_correction")),
      pattern = quote(map(iad_recovery_exact_file)), phase = "recover_correction", format = "file", iteration = "list"),
    add("final_loo_files", quote(iad_finalise_loo(iad_recovery_loo_audit, iad_recovery_correction_file)), format = "file"),
    add("accepted_file", quote(iad_accept_recovery(iad_recovery_accepted_record,
      iad_recovery_final_loo_files, iad_frozen_checkpoints)), phase = "recover_accept", format = "file"),
    add("report_source_files", quote(iad_recovery_reporting_inputs()), format = "file"),
    add("report_record", quote(iad_recovery_reporting_records(iad_recovery_accepted_file,
      iad_frozen_checkpoints, iad_recovery_report_source_files)), phase = "recover_report", iteration = "list"),
    add("report_files", quote(iad_recovery_report_scope(iad_recovery_report_record,
      iad_frozen_checkpoints)), pattern = quote(map(iad_recovery_report_record)),
      phase = "recover_report", format = "file", iteration = "list"),
    add("reporting_audit_files", quote(iad_recovery_reporting_audit(iad_recovery_report_files,
      iad_recovery_accepted_file)), format = "file")
  )
}
