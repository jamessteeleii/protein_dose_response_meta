# Targets for one age or duration family ----

# Assemble the pilot, fitting, diagnostics, holdout scoring and reporting
# stages for one family. The phase checks control which stages are enabled.

# Each exposure has independent pilot, diagnostic, LOO and reporting gates.
# All families share one future scheduler and one cache, never three runs.
iad_family_targets <- function(family) {
  stopifnot(length(family) == 1L, family %in% iad_families())
  prefix <- paste0("iad_", family, "_")
  ref <- function(suffix) as.name(paste0(prefix, suffix))
  add <- function(suffix, command, pattern = NULL, phase = NULL, format = "qs", iteration = "vector") {
    args <- list(name = paste0(prefix, suffix), command = command, format = format,
      iteration = iteration, error = "continue")
    if (!is.null(pattern)) args$pattern <- pattern
    if (!is.null(phase)) args$cue <- iad_cue(phase)
    do.call(targets::tar_target_raw, args, quote = TRUE)
  }
  list(
    add("specification", bquote(Filter(function(s) s$family == .(family), iad_specification)), iteration = "list"),
    add("pilot_specification", bquote(iad_pilot_specs(.(ref("specification")))), iteration = "list"),
    add("remaining_specification", bquote(Filter(function(s) !(s$interaction == "both" &&
      s$parent_index %in% c(6L, 26L)), .(ref("specification")))), iteration = "list"),
    add("pilot_file", bquote({ iad_preparation_files; iad_batch_storage_gate_file;
      iad_fit(.(ref("pilot_specification")), iad_checkpoints, iad_base_priors, iad_input_preflight, "pilot") }),
      pattern = bquote(map(.(ref("pilot_specification")))), phase = "pilot", format = "file", iteration = "list"),
    add("pilot_diagnostic", bquote(iad_diagnose(.(ref("pilot_file")))),
      pattern = bquote(map(.(ref("pilot_file")))), iteration = "list"),
    add("pilot_summary", bquote(iad_diagnostic_summary(.(ref("pilot_diagnostic"))))),
    add("pilot_audit_files", bquote(iad_diagnostic_files(.(ref("pilot_diagnostic")), "pilot", .(family))), format = "file"),
    add("pilot_gate_file", bquote({ .(ref("pilot_audit_files"));
      iad_family_pilot_gate(.(ref("pilot_summary")), .(family), iad_checkpoints) }), format = "file"),
    add("fit_file", bquote({ .(ref("pilot_gate_file")); iad_batch_storage_gate_file;
      iad_fit(.(ref("remaining_specification")), iad_checkpoints, iad_base_priors,
        iad_input_preflight, "fit", .(ref("pilot_summary"))) }),
      pattern = bquote(map(.(ref("remaining_specification")))), phase = "fit", format = "file", iteration = "list"),
    add("initial_fit_checkpoint_file", bquote(iad_family_fit_checkpoint(
      c(.(ref("pilot_file")), .(ref("fit_file"))), .(ref("specification")), .(family))), format = "file"),
    add("initial_diagnostic", bquote(iad_diagnose(.(ref("fit_file")))),
      pattern = bquote(map(.(ref("fit_file")))), iteration = "list"),
    add("initial_summary", bquote(iad_diagnostic_summary(c(.(ref("pilot_diagnostic")), .(ref("initial_diagnostic")))))),
    add("initial_audit_files", bquote(iad_diagnostic_files(c(.(ref("pilot_diagnostic")),
      .(ref("initial_diagnostic"))), "initial", .(family))), format = "file"),
    add("remedy_specification", bquote(iad_remedy_specs(.(ref("specification")), iad_remedy_file,
      .(ref("initial_summary")), family = .(family))), iteration = "list"),
    add("remedy_fit_file", bquote(iad_fit_remedy(.(ref("remedy_specification")), iad_checkpoints,
      iad_base_priors, iad_input_preflight, family = .(family))),
      pattern = bquote(map(.(ref("remedy_specification")))), phase = "remedy", format = "file", iteration = "list"),
    add("remedy_diagnostic", bquote(iad_diagnose_remedy(.(ref("remedy_fit_file")))),
      pattern = bquote(map(.(ref("remedy_fit_file")))), iteration = "list"),
    add("accepted_record", bquote(iad_accepted(.(ref("initial_summary")),
      iad_diagnostic_summary(.(ref("remedy_diagnostic"))), .(ref("specification")), iad_checkpoints)), iteration = "list"),
    add("accepted_audit_file", bquote(iad_csv(dplyr::bind_rows(.(ref("accepted_record"))),
      iad_path(.(family), "diagnostics", "accepted.csv"))), format = "file"),
    add("integrated_loo_file", bquote(iad_compute_loo(.(ref("accepted_record")))),
      pattern = bquote(map(.(ref("accepted_record")))), phase = "loo", format = "file", iteration = "list"),
    add("integrated_audit", bquote(iad_loo_audit(.(ref("integrated_loo_file"))))),
    add("integrated_audit_files", bquote(iad_write_loo_audit(.(ref("integrated_audit")), .(family))), format = "file"),
    add("exact_record", bquote(.(ref("integrated_audit"))$cases), iteration = "list"),
    add("exact_file", bquote(iad_exact(.(ref("exact_record")), family = .(family), compact = TRUE)),
      pattern = bquote(map(.(ref("exact_record")))), phase = "exact", format = "file", iteration = "list"),
    add("exact_diagnostic", bquote(iad_exact_diagnostics(.(ref("exact_file")))),
      pattern = bquote(map(.(ref("exact_file")))), iteration = "list"),
    add("exact_audit_file", bquote(iad_csv(dplyr::bind_rows(.(ref("exact_diagnostic"))),
      iad_path(.(family), "loo", "exact_diagnostics.csv"))), format = "file"),
    add("exact_remedy_plan", bquote(iad_exact_remedy_plan(.(ref("exact_file")),
      iad_exact_remedy_ledger_file, .(family))), iteration = "list"),
    add("exact_review_file", bquote(iad_review_exact(.(ref("exact_remedy_plan")))),
      pattern = bquote(map(.(ref("exact_remedy_plan")))), phase = "exact_remedy", format = "file", iteration = "list"),
    add("exact_review_diagnostic", bquote(iad_exact_diagnostics(.(ref("exact_review_file")))),
      pattern = bquote(map(.(ref("exact_review_file")))), iteration = "list"),
    add("exact_review_audit_file", bquote(iad_csv(dplyr::bind_rows(.(ref("exact_review_diagnostic"))),
      iad_path(.(family), "loo", "exact_review_diagnostics.csv"))), format = "file"),
    add("exact_correction_file", bquote(iad_score_exact(.(ref("exact_review_file")), family = .(family))),
      pattern = bquote(map(.(ref("exact_review_file")))), phase = "correction", format = "file", iteration = "list"),
    add("final_loo_files", bquote(iad_finalise_loo(.(ref("integrated_audit")), .(ref("exact_correction_file")))), format = "file"),
    add("combined_record", bquote({ iad_recovery_reporting_gate_file;
      iad_combine(iad_checkpoints, .(ref("accepted_record")), .(ref("final_loo_files")), families = .(family)) }),
      phase = "combine", iteration = "list"),
    add("combination_audit_files", bquote(iad_write_combination(.(ref("combined_record")), .(family))), format = "file"),
    add("report_files", bquote(iad_report_scope(.(ref("combined_record")), iad_checkpoints, iad_profile_plan)),
      pattern = bquote(map(.(ref("combined_record")))), phase = "report", format = "file", iteration = "list"),
    add("report_checkpoint_file", bquote(iad_family_reporting_checkpoint(.(ref("report_files")), .(family))), format = "file")
  )
}
