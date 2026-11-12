# Shared readers and presentation helpers for publication supplements.
# Reports use durable verification bundles rather than the disposable targets cache.

report_project_root <- function() {
  if (file.exists("_targets.R")) {
    normalizePath(".", winslash = "/", mustWork = TRUE)
  } else {
    normalizePath("..", winslash = "/", mustWork = TRUE)
  }
}

report_input_root <- function() {
  file.path(report_project_root(), "manuscript", "generated", "reporting")
}

report_bundle_root <- function() {
  file.path(report_input_root(), "bundles")
}

bundle_path <- function(bundle, filename) {
  file.path(report_bundle_root(), bundle, filename)
}

bundle_file <- function(bundle, candidates, required = TRUE) {
  paths <- bundle_path(bundle, candidates)
  found <- paths[file.exists(paths)]
  if (!length(found)) {
    if (!required) return(NA_character_)
    stop(
      "Missing required verification artifact for ", bundle, ": ",
      paste(candidates, collapse = ", "), call. = FALSE
    )
  }
  found[[1L]]
}

read_bundle <- function(bundle, candidates, required = TRUE) {
  path <- bundle_file(bundle, candidates, required = required)
  if (is.na(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

report_outcome_labels <- c(
  whole_body = "Whole-body lean/fat-free mass",
  hypertrophy = "Measures of muscle size or mass",
  strength = "Strength",
  `function` = "Physical function"
)

report_outcome_group_labels <- stats::setNames(
  report_outcome_labels,
  c("whole body lean mass", "direct muscle size/mass", "strength", "function")
)

report_outcome_group_label <- function(x) {
  keys <- sub("^whole-body", "whole body", tolower(as.character(x)))
  dplyr::coalesce(unname(report_outcome_group_labels[keys]), as.character(x))
}

report_rt_label <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    x %in% c("YES", "RT", "0") ~ "RT",
    x %in% c("NO", "NO RT", "1") ~ "No RT",
    TRUE ~ as.character(x)
  )
}

report_rt_display <- c(
  "RT" = "Resistance Training: Yes",
  "No RT" = "Resistance Training: No"
)

report_facet_labels <- function() {
  ggplot2::labeller(
    rt = report_rt_display, `RT condition` = report_rt_display,
    outcome = ggplot2::label_wrap_gen(24),
    Outcome = ggplot2::label_wrap_gen(24)
  )
}

report_curve_colours <- c(
  "Primary achieved intake" = "black",
  "Primary absolute total intake" = "black",
  "Primary prescribed supplement" = "black",
  "Primary measurement eligibility" = "black",
  "Mundlak decomposition" = "#0072B2",
  "Strict within-study" = "#D55E00",
  "Achieved intake, common sample" = "#009E73",
  "Assignment-derived, common sample" = "#8E5C8A",
  "Errors-in-variables" = "#8C6500",
  "Residual correlation" = "#537487",
  "Registered measurement methods" = "#0072B2",
  "MRI/CT/ultrasound/biopsy only" = "#D55E00",
  "Addition" = "#0072B2",
  "Substitution" = "#D55E00"
)

report_dose_axes <- c(
  relative = "Protein Intake (g/kg/day)",
  total = "Achieved Total Protein Intake (g/day)",
  supplement = "Prescribed Supplemental Protein (g/day)"
)

normalise_predictions <- function(data, analysis, dose_unit) {
  rt <- if ("rt_condition" %in% names(data)) {
    data$rt_condition
  } else if ("resistance_exercise" %in% names(data)) {
    data$resistance_exercise
  } else {
    data$resistance_exercise_code
  }
  strategy <- if ("protein_strategy" %in% names(data)) {
    data$protein_strategy
  } else {
    rep(NA_character_, nrow(data))
  }
  tibble::tibble(
    analysis = analysis,
    outcome = unname(report_outcome_labels[data$outcome]),
    rt = report_rt_label(rt),
    strategy = strategy,
    dose = data$protein_intake,
    dose_unit = dose_unit,
    estimate = data$estimate,
    conf_low = data$conf.low,
    conf_high = data$conf.high
  )
}

normalise_increments <- function(data, analysis, dose_unit) {
  rt <- if ("rt_condition" %in% names(data)) {
    data$rt_condition
  } else if ("resistance_exercise" %in% names(data)) {
    data$resistance_exercise
  } else {
    data$resistance_exercise_code
  }
  strategy <- if ("protein_strategy" %in% names(data)) {
    data$protein_strategy
  } else {
    rep(NA_character_, nrow(data))
  }
  tibble::tibble(
    analysis = analysis,
    outcome = unname(report_outcome_labels[data$outcome]),
    rt = report_rt_label(rt),
    strategy = strategy,
    dose = data$protein_intake,
    dose_unit = dose_unit,
    estimate = data$increment,
    conf_low = data$increment_conf_low,
    conf_high = data$increment_conf_high,
    reference = data$reference_intake
  )
}

analysis_sample_row <- function(data, analysis, within = FALSE) {
  if (within) {
    source_arms <- unique(c(data$protein_arm_no, data$comparator_arm_no))
    arms <- length(source_arms)
    comparisons <- dplyr::n_distinct(data$arm_no)
    unit <- "Matched within-study contrast"
  } else {
    arms <- dplyr::n_distinct(data$arm_no)
    comparisons <- NA_integer_
    unit <- "Arm-level effect"
  }
  tibble::tibble(
    Analysis = analysis,
    Studies = dplyr::n_distinct(data$study_no),
    `Source arms` = arms,
    `Paired comparisons` = comparisons,
    `Analysis rows` = nrow(data),
    `Analysis unit` = unit
  )
}

diagnostic_row <- function(bundle, analysis) {
  sampling <- read_bundle(
    bundle,
    c(
      "final_sampling_diagnostic_summary.csv.gz",
      "sampling_diagnostic_summary.csv.gz",
      "diagnostic_summary.csv.gz"
    )
  )
  pass_column <- intersect(
    c("accepted_diagnostic_pass", "diagnostic_pass"), names(sampling)
  )[[1L]]
  exact <- read_bundle(
    bundle,
    c(
      "exact_refit_diagnostic_summary.csv.gz",
      "exact_loo_diagnostic_summary.csv.gz",
      "final_exact_refit_diagnostic_summary.csv.gz",
      "loo_exact_refit_summary.csv.gz",
      "accepted_exact_refit_diagnostic_summary.csv.gz"
    ),
    required = FALSE
  )
  exact_n <- if (is.null(exact)) 0L else nrow(exact)
  exact_pass <- if (is.null(exact) || !nrow(exact)) {
    0L
  } else if ("diagnostic_pass" %in% names(exact)) {
    sum(exact$diagnostic_pass %in% TRUE, na.rm = TRUE)
  } else {
    exact_n
  }
  tibble::tibble(
    Analysis = analysis,
    `Accepted fits` = nrow(sampling),
    `Fits passing final diagnostics` = sum(sampling[[pass_column]] %in% TRUE, na.rm = TRUE),
    `Exact LOO refits` = exact_n,
    `Exact refits passing` = exact_pass
  )
}

local_equivalence_summary <- function(bundle, analysis) {
  x <- read_bundle(bundle, "local_equivalence.csv.gz")
  x |>
    dplyr::mutate(
      Analysis = analysis,
      Outcome = unname(report_outcome_labels[outcome]),
      `RT condition` = report_rt_label(rt_condition)
    ) |>
    dplyr::select(
      Analysis, dplyr::any_of("protein_strategy"), Outcome, `RT condition`,
      equivalent_radius, equivalent_radius_q10, equivalent_radius_q90,
      supported_current_intakes
    )
}

report_theme <- function() {
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0),
      panel.spacing = grid::unit(0.85, "lines"),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(hjust = 0)
    )
}

format_count <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}
