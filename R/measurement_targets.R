# Targets for the measurement restrictions ----

# Define the two restriction branches from frozen input selection through the
# accepted fits, holdout scores and manuscript outputs.

# Target definitions are generated here to keep the two measurement restrictions
# on identical workflows without duplicating the analysis implementation.
make_measurement_targets <- function() {
  target <- targets::tar_target_raw
  common <- list(
    target("measurement_checkpoint_files", quote(measurement_checkpoint_inputs()), format = "file"),
    target("measurement_checkpoint", quote(read_measurement_checkpoint(measurement_checkpoint_files))),
    target("measurement_annotation_file", quote(here::here(
      "sensitivity", "measurement-annotations.csv"
    )), format = "file"),
    target("measurement_annotations", quote(read_measurement_annotations(measurement_annotation_file))),
    target("measurement_remedy_decision_file", quote(here::here(
      "sensitivity", "measurement-remedies.csv"
    )), format = "file"),
    target("measurement_remedy_decisions", quote(utils::read.csv(
      measurement_remedy_decision_file, stringsAsFactors = FALSE
    )))
  )

  variants <- lapply(measurement_variants(), function(variant) {
    prefix <- paste0("measurement_", variant, "_")
    symbol <- function(suffix) as.name(paste0(prefix, suffix))
    substitutions <- list(V = variant,
      SEL = symbol("selection"), SPEC = symbol("specification"),
      PRE = symbol("preflight"), FIT = symbol("fit_file"),
      PILOTSPEC = symbol("pilot_specification"), FULLSPEC = symbol("fit_specification"),
      PILOTFIT = symbol("pilot_file"), PILOT = symbol("pilot_audit"),
      INITDIAG = symbol("initial_diagnostic_audit"),
      REMSPEC = symbol("remedy_specification"), REMFIT = symbol("remedy_file"),
      MAN = symbol("accepted_manifest"), DIAG = symbol("diagnostic_audit"),
      LOOREC = symbol("loo_record"), LOO = symbol("integrated_loo_file"),
      LOOMAN = symbol("integrated_manifest"), EXREC = symbol("exact_record"),
      EXACT = symbol("exact_file"), FINAL = symbol("final_loo_files"),
      WEIGHT = symbol("stacking_weights"), AVG = symbol("model_average_file"),
      REPORT = symbol("report_file"))
    command <- function(x) {
      if (is.symbol(x) && as.character(x) %in% names(substitutions)) {
        return(substitutions[[as.character(x)]])
      }
      if (is.call(x)) return(as.call(lapply(as.list(x), command)))
      x
    }
    add <- function(suffix, expression, pattern = NULL, phase = NULL,
                    format = "qs", iteration = "vector", error = "stop") {
      arguments <- list(name = paste0(prefix, suffix), command = command(expression),
        format = format, iteration = iteration, error = error)
      if (!is.null(pattern)) arguments$pattern <- command(pattern)
      if (!is.null(phase)) arguments$cue <- targets::tar_cue(
        mode = if (measurement_phase_enabled(phase)) "thorough" else "never"
      )
      do.call(target, arguments, quote = TRUE)
    }
    list(
      add("selection", quote(select_measurement_data(measurement_checkpoint$data, V,
        annotations = measurement_annotations))),
      add("specification", quote(make_measurement_specifications(measurement_checkpoint, V)),
        iteration = "list"),
      add("preflight", quote(measurement_input_preflight(SEL, SPEC, steele_priors))),
      add("preparation_files", quote(write_measurement_preparation(SEL, SPEC, PRE)), format = "file"),
      add("pilot_specification", quote(measurement_pilot_specifications(SPEC)), iteration = "list"),
      add("fit_specification", quote(SPEC[-c(2L, 26L)]), iteration = "list"),
      add("pilot_file", quote(fit_measurement_candidate(PILOTSPEC, SEL, steele_priors, PRE)),
        pattern = quote(map(PILOTSPEC)), phase = "pilot", format = "file", iteration = "list", error = "continue"),
      add("pilot_audit", quote(measurement_pilot_audit(PILOTFIT, V))),
      add("pilot_audit_file", quote(write_measurement_pilot_audit(PILOT)), format = "file"),
      add("fit_file", quote(fit_measurement_candidate(FULLSPEC, SEL, steele_priors, PRE, PILOT)),
        pattern = quote(map(FULLSPEC)), phase = "fit", format = "file", iteration = "list", error = "continue"),
      add("initial_diagnostic_audit", quote(measurement_diagnostic_audit(c(PILOTFIT, FIT), V))),
      add("initial_diagnostic_files", quote(write_measurement_diagnostics(INITDIAG, "initial")),
        format = "file"),
      add("remedy_specification", quote(make_measurement_remedy_records(
        SPEC, measurement_remedy_decisions, V, initial_diagnostics = INITDIAG
      )), iteration = "list"),
      add("remedy_file", quote(fit_measurement_remedy(REMSPEC, SEL, steele_priors, PRE)),
        pattern = quote(map(REMSPEC)), phase = "remedy", format = "file", iteration = "list", error = "continue"),
      add("accepted_manifest", quote(measurement_accepted_manifest(c(PILOTFIT, FIT), REMFIT, V, SEL))),
      add("diagnostic_audit", quote(measurement_diagnostic_audit(MAN$fit_path, V))),
      add("diagnostic_files", quote(write_measurement_diagnostics(DIAG, "accepted")), format = "file"),
      add("loo_record", quote(measurement_loo_records(MAN, DIAG, V)), iteration = "list"),
      add("integrated_loo_file", quote(compute_measurement_integrated_loo(LOOREC)),
        pattern = quote(map(LOOREC)), phase = "loo", format = "file", iteration = "list", error = "continue"),
      add("integrated_manifest", quote(measurement_integrated_manifest(LOO, V))),
      add("integrated_audit_files", quote(write_measurement_loo_audit(LOOMAN, V)), format = "file"),
      add("exact_record", quote(measurement_exact_records(LOOMAN, V)), iteration = "list"),
      add("exact_file", quote(fit_measurement_exact(EXREC)), pattern = quote(map(EXREC)),
        phase = "exact", format = "file", iteration = "list", error = "continue"),
      add("final_loo_files", quote(finalise_measurement_loo(LOOMAN, EXACT, V)), format = "file"),
      add("stacking_weights", quote(measurement_stacking(FINAL, V))),
      add("model_average_file", quote(write_measurement_model_average(MAN, WEIGHT, SEL)),
        phase = "report", format = "file"),
      add("report_file", quote(write_measurement_report(AVG, SEL)), phase = "report", format = "file"),
      add("verification_files", quote(write_measurement_checkpoint(
        REPORT, SEL, SPEC, DIAG, MAN, FINAL, measurement_checkpoint
      )), phase = "report", format = "file"),
      add("figure_files", quote(write_measurement_figure(REPORT, SEL)), phase = "report", format = "file")
    )
  })

  c(common, unlist(variants, recursive = FALSE))
}
