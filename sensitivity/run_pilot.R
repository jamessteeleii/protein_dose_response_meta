args <- commandArgs(trailingOnly = TRUE)

families <- c("mundlak", "within", "exposure", "rescor", "eiv")
pilot_names <- c(families, "mundlak-remedy")
stages <- c("fit", "diagnostics", "summary")

if (length(args) == 1L && identical(args[[1L]], "structure")) {
  Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
  targets::tar_make(sensitivity_pilot_structure_audit)
  quit(save = "no", status = 0L)
}

if (
  length(args) != 2L || !args[[1L]] %in% pilot_names ||
    !args[[2L]] %in% stages
) {
  stop(
    paste0(
      "Usage: Rscript sensitivity/run_pilot.R structure\n",
      "   or: Rscript sensitivity/run_pilot.R ",
      "{mundlak|mundlak-remedy|within|exposure|rescor|eiv} ",
      "{fit|diagnostics|summary}"
    ),
    call. = FALSE
  )
}

family <- args[[1L]]
stage <- args[[2L]]
targets_by_family <- list(
  mundlak = list(
    fit = "sensitivity_mundlak_pilot_fit",
    diagnostics = "sensitivity_mundlak_pilot_diagnostics",
    summary = "sensitivity_mundlak_pilot_summary"
  ),
  `mundlak-remedy` = list(
    fit = "sensitivity_mundlak_remedy_fit",
    diagnostics = "sensitivity_mundlak_remedy_diagnostics",
    summary = "sensitivity_mundlak_remedy_summary"
  ),
  within = list(
    fit = "sensitivity_within_pilot_fit",
    diagnostics = "sensitivity_within_pilot_diagnostics",
    summary = "sensitivity_within_pilot_summary"
  ),
  exposure = list(
    fit = c(
      "sensitivity_achieved_pilot_fit",
      "sensitivity_assignment_pilot_fit"
    ),
    diagnostics = c(
      "sensitivity_achieved_pilot_diagnostics",
      "sensitivity_assignment_pilot_diagnostics"
    ),
    summary = "sensitivity_exposure_pilot_comparison"
  ),
  rescor = list(
    fit = "sensitivity_rescor_pilot_fit",
    diagnostics = "sensitivity_rescor_pilot_diagnostics",
    summary = "sensitivity_rescor_pilot_diagnostics"
  ),
  eiv = list(
    fit = "sensitivity_eiv_pilot_fit",
    diagnostics = "sensitivity_eiv_pilot_diagnostics",
    summary = "sensitivity_eiv_pilot_summary"
  )
)

enabled_family <- if (identical(family, "mundlak-remedy")) {
  "mundlak"
} else {
  family
}
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = enabled_family)
target <- unname(targets_by_family[[family]][[stage]])
targets::tar_make(names = tidyselect::all_of(target))
