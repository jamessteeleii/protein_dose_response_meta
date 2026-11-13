# Protein-manipulation-strategy sensitivity analysis.
#
# The study-level audit is intentionally external to the primary data. It is
# joined only after the four outcome datasets have been prepared, so audit
# edits cannot invalidate the completed primary or earlier sensitivity fits.

read_protein_strategy_audit <- function(path) {
  audit <- readr::read_csv(
    path,
    col_types = readr::cols(
      .default = readr::col_character(),
      study_no = readr::col_character(),
      year = readr::col_character(),
      arm_count = readr::col_integer()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  required <- c(
    "study_no", "pmid", "doi", "authors", "year", "title",
    "arm_count", "arm_mix", "primary_strategy", "strategy_detail",
    "energy_relation", "displaced_energy", "confidence", "coding_reason",
    "reported_intervention", "evidence_excerpt", "page"
  )
  missing <- setdiff(required, names(audit))
  if (length(missing) > 0L) {
    stop(
      "The protein-strategy audit is missing: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  valid_strategy <- c(
    "Addition", "Substitution", "Mixed",
    "Protein-only comparison", "Unclear"
  )
  valid_confidence <- c("High", "Moderate", "Low")
  # The public historical audit leaves verbatim source-text fields blank.
  required_content <- setdiff(
    required, c("pmid", "doi", "page", "reported_intervention", "evidence_excerpt")
  )
  missing_content <- vapply(
    audit[required_content],
    function(value) any(is.na(value) | !nzchar(trimws(as.character(value)))),
    logical(1L)
  )
  clearly_matched <- c(
    "Matched", "Matched/controlled", "Approximately matched",
    "Matched by powder dose", "Probably matched"
  )
  clearly_not_matched <- c("Not matched", "Not constrained")
  strategy_energy_conflict <-
    (audit$primary_strategy == "Addition" &
       audit$energy_relation %in% clearly_matched) |
    (audit$primary_strategy == "Substitution" &
       audit$energy_relation %in% clearly_not_matched)
  if (
    nrow(audit) != 109L || anyNA(audit$study_no) ||
      anyDuplicated(audit$study_no) ||
      any(!audit$primary_strategy %in% valid_strategy) ||
      any(!audit$confidence %in% valid_confidence) ||
      any(missing_content) || any(strategy_energy_conflict)
  ) {
    stop(
      "The final protein-strategy audit failed its schema or value audit. ",
      "Run sensitivity/finalise_strategy_audit_csv.mjs and inspect its report.",
      call. = FALSE
    )
  }
  audit |>
    dplyr::mutate(
      # The lean CSV is the completed, confirmed audit. These working fields
      # are derived at read time so review-only flags cannot become stale.
      review_status = "Confirmed",
      classification_confirmed = TRUE,
      principal_eligible = primary_strategy %in%
        c("Addition", "Substitution"),
      strict_eligible = principal_eligible & confidence == "High"
    )
}

read_update_2026_protein_strategy_audit <- function(
  historical_path,
  additions_path
) {
  historical <- read_protein_strategy_audit(historical_path)
  additions <- readr::read_csv(
    additions_path,
    col_types = readr::cols(
      .default = readr::col_character(),
      study_no = readr::col_character(),
      year = readr::col_character(),
      arm_count = readr::col_integer()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  audit_columns <- c(
    "study_no", "pmid", "doi", "authors", "year", "title",
    "arm_count", "arm_mix", "primary_strategy", "strategy_detail",
    "energy_relation", "displaced_energy", "confidence", "coding_reason",
    "reported_intervention", "evidence_excerpt", "page"
  )
  provenance_columns <- c(
    "zotero_key", "source_filename", "source_sha256", "review_status"
  )
  missing <- setdiff(
    c(audit_columns, provenance_columns),
    names(additions)
  )
  if (length(missing) > 0L) {
    stop(
      "The updated protein-strategy audit additions are missing: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  expected_studies <- as.character(c(163L, 165L, 166L, 167L, 171L, 172L))
  required_content <- setdiff(
    c(audit_columns, provenance_columns),
    "pmid"
  )
  missing_content <- vapply(
    additions[required_content],
    function(value) any(is.na(value) | !nzchar(trimws(as.character(value)))),
    logical(1L)
  )
  clearly_matched <- c(
    "Matched", "Matched/controlled", "Approximately matched",
    "Matched by powder dose", "Probably matched"
  )
  clearly_not_matched <- c("Not matched", "Not constrained")
  strategy_energy_conflict <-
    (additions$primary_strategy == "Addition" &
       additions$energy_relation %in% clearly_matched) |
    (additions$primary_strategy == "Substitution" &
       additions$energy_relation %in% clearly_not_matched)
  if (
    nrow(additions) != 6L ||
      !setequal(additions$study_no, expected_studies) ||
      anyDuplicated(additions$study_no) ||
      any(!additions$primary_strategy %in% c("Addition", "Substitution")) ||
      any(!additions$confidence %in% c("High", "Moderate", "Low")) ||
      any(!additions$review_status %in% c("Provisional", "Confirmed")) ||
      any(missing_content) || any(strategy_energy_conflict)
  ) {
    stop(
      "The updated protein-strategy audit additions failed their schema ",
      "or value audit.",
      call. = FALSE
    )
  }
  historical <- historical |>
    dplyr::mutate(
      zotero_key = NA_character_,
      source_filename = NA_character_,
      source_sha256 = NA_character_
    ) |>
    dplyr::select(
      dplyr::all_of(audit_columns),
      zotero_key,
      source_filename,
      source_sha256,
      review_status,
      classification_confirmed,
      principal_eligible,
      strict_eligible
    )
  additions <- additions |>
    dplyr::mutate(
      classification_confirmed = review_status == "Confirmed",
      principal_eligible = primary_strategy %in%
        c("Addition", "Substitution"),
      strict_eligible = principal_eligible & confidence == "High"
    ) |>
    dplyr::select(dplyr::all_of(names(historical)))
  audit <- dplyr::bind_rows(historical, additions) |>
    dplyr::arrange(as.integer(study_no))
  if (
    nrow(audit) != 115L ||
      anyDuplicated(audit$study_no) ||
      !all(expected_studies %in% audit$study_no)
  ) {
    stop(
      "The combined updated protein-strategy audit is incomplete.",
      call. = FALSE
    )
  }
  audit
}

audit_update_2026_protein_strategy <- function(audit, data_by_outcome) {
  expected_new_studies <- as.character(c(163L, 165L, 166L, 167L, 171L, 172L))
  if (!is.list(data_by_outcome) || length(data_by_outcome) != 4L) {
    stop("The four updated outcome datasets are required.", call. = FALSE)
  }
  analysis_studies <- data_by_outcome |>
    purrr::map(~ as.character(.x$study_no)) |>
    unlist(use.names = FALSE) |>
    unique() |>
    sort()
  audit_studies <- sort(unique(as.character(audit$study_no)))
  missing_studies <- setdiff(analysis_studies, audit_studies)
  review_rows <- audit |>
    dplyr::filter(study_no %in% expected_new_studies) |>
    dplyr::arrange(as.integer(study_no))
  pending <- review_rows |>
    dplyr::filter(!classification_confirmed)
  structural_checks <- tibble::tibble(
    check = c(
      "combined strategy audit contains 115 unique studies",
      "updated analysis contains the expected 95 studies",
      "every updated-analysis study has a strategy classification",
      "the six updated-search studies have audit records",
      "the updated audit retains 20 historical audit-only studies",
      "all six updated records identify their Zotero source PDF"
    ),
    passed = c(
      nrow(audit) == 115L && !anyDuplicated(audit$study_no),
      length(analysis_studies) == 95L,
      length(missing_studies) == 0L,
      setequal(review_rows$study_no, expected_new_studies),
      length(setdiff(audit_studies, analysis_studies)) == 20L,
      nrow(review_rows) == 6L &&
        all(nzchar(review_rows$zotero_key)) &&
        all(nzchar(review_rows$source_filename)) &&
        all(nchar(review_rows$source_sha256) == 64L)
    )
  )
  if (!all(structural_checks$passed)) {
    stop(
      "The updated protein-strategy audit failed structural validation.",
      call. = FALSE
    )
  }
  list(
    review_rows = review_rows,
    status_summary = review_rows |>
      dplyr::count(primary_strategy, confidence, review_status),
    structural_checks = structural_checks,
    analysis_studies = length(analysis_studies),
    audited_studies = nrow(audit),
    pending_confirmations = nrow(pending),
    audit_ready = nrow(pending) == 0L,
    next_step = if (nrow(pending) == 0L) {
      "All updated classifications are confirmed; strategy preflight may proceed."
    } else {
      paste0(
        nrow(pending),
        " updated classifications remain provisional and require author review."
      )
    }
  )
}

write_update_2026_strategy_audit_review <- function(
  audit,
  output_dir = here::here(
    "results", "update_2026", "sensitivity", "strategy", "audit"
  )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    classifications = file.path(output_dir, "classifications_for_review.csv"),
    status = file.path(output_dir, "review_status_summary.csv"),
    checks = file.path(output_dir, "structural_checks.csv")
  )
  readr::write_csv(audit$review_rows, paths[["classifications"]], na = "")
  readr::write_csv(audit$status_summary, paths[["status"]], na = "")
  readr::write_csv(audit$structural_checks, paths[["checks"]], na = "")
  unname(paths)
}

join_protein_strategy_audit <- function(data, audit) {
  if (!is.data.frame(data) || !"study_no" %in% names(data)) {
    stop("Prepared outcome data must contain `study_no`.", call. = FALSE)
  }
  lookup <- audit |>
    dplyr::transmute(
      strategy_study_key = as.character(study_no),
      protein_strategy = primary_strategy,
      strategy_classification_confirmed = classification_confirmed,
      strategy_principal_eligible = principal_eligible,
      strategy_strict_eligible = strict_eligible
    )
  if (anyDuplicated(lookup$strategy_study_key)) {
    stop("The protein-strategy lookup is not one row per study.", call. = FALSE)
  }
  output <- data |>
    dplyr::mutate(strategy_study_key = as.character(study_no)) |>
    dplyr::left_join(lookup, by = "strategy_study_key")
  if (
    nrow(output) != nrow(data) ||
      anyNA(output$strategy_principal_eligible)
  ) {
    missing_studies <- unique(output$strategy_study_key[
      is.na(output$strategy_principal_eligible)
    ])
    stop(
      "The strategy audit did not match every prepared-data study",
      if (length(missing_studies)) {
        paste0(": ", paste(missing_studies, collapse = ", "))
      } else "",
      ".",
      call. = FALSE
    )
  }
  output
}

finalise_strategy_analysis_data <- function(data, analysis) {
  if (!analysis %in% c("principal", "strict")) {
    stop("Unknown protein-strategy analysis sample.", call. = FALSE)
  }
  eligible <- if (analysis == "principal") {
    data$strategy_principal_eligible
  } else {
    data$strategy_strict_eligible
  }
  output <- data[eligible, , drop = FALSE] |>
    dplyr::mutate(
      strategy_code = ifelse(protein_strategy == "Addition", 1, 0),
      strategy_factor = factor(
        protein_strategy,
        levels = c("Substitution", "Addition")
      ),
      study_factor = factor(study_no),
      sensitivity_strategy_sample = analysis
    ) |>
    dplyr::select(
      -strategy_study_key,
      -strategy_classification_confirmed,
      -strategy_principal_eligible,
      -strategy_strict_eligible
    )
  study_check <- output |>
    dplyr::distinct(study_no, protein_strategy) |>
    dplyr::count(study_no, name = "strategies")
  if (
    nrow(output) < 1L || anyNA(output$strategy_code) ||
      any(study_check$strategies != 1L) ||
      !setequal(unique(output$strategy_code), c(0, 1))
  ) {
    stop("The strategy analysis sample failed validation.", call. = FALSE)
  }
  output
}

summarise_strategy_support <- function(data, outcome, analysis) {
  arms <- data |>
    dplyr::distinct(
      study_no, arm_no, protein_strategy, resistance_exercise_code,
      protein_intake_raw
    )
  arms |>
    dplyr::group_by(protein_strategy, resistance_exercise_code) |>
    dplyr::summarise(
      rows = dplyr::n(),
      arms = dplyr::n_distinct(arm_no),
      studies = dplyr::n_distinct(study_no),
      observed_minimum = min(protein_intake_raw),
      central_90_minimum = as.numeric(stats::quantile(
        protein_intake_raw, 0.05
      )),
      median = stats::median(protein_intake_raw),
      central_90_maximum = as.numeric(stats::quantile(
        protein_intake_raw, 0.95
      )),
      observed_maximum = max(protein_intake_raw),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      outcome = outcome,
      analysis = analysis,
      resistance_exercise = ifelse(
        resistance_exercise_code == 0, "YES", "NO"
      ),
      .before = 1L
    )
}

prepare_protein_strategy_data <- function(data_by_outcome, audit) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  if (
    !is.list(data_by_outcome) ||
      !all(outcomes %in% names(data_by_outcome))
  ) {
    stop("The four prepared outcome datasets are required.", call. = FALSE)
  }
  joined <- purrr::map(data_by_outcome[outcomes],
                       join_protein_strategy_audit, audit = audit)
  principal <- purrr::map(joined, finalise_strategy_analysis_data,
                          analysis = "principal")
  strict <- purrr::map(joined, finalise_strategy_analysis_data,
                       analysis = "strict")
  support <- dplyr::bind_rows(
    purrr::imap_dfr(principal, summarise_strategy_support,
                    analysis = "principal"),
    purrr::imap_dfr(strict, summarise_strategy_support,
                    analysis = "strict")
  )
  principal_audit <- audit |>
    dplyr::filter(principal_eligible)
  blocking <- principal_audit |>
    dplyr::filter(!classification_confirmed) |>
    dplyr::select(
      study_no, authors, year, primary_strategy, confidence, review_status
    )
  outcome_summary <- purrr::imap_dfr(principal, function(data, outcome) {
    tibble::tibble(
      outcome = outcome,
      rows = nrow(data),
      arms = dplyr::n_distinct(data$arm_no),
      studies = dplyr::n_distinct(data$study_no),
      addition_studies = dplyr::n_distinct(
        data$study_no[data$protein_strategy == "Addition"]
      ),
      substitution_studies = dplyr::n_distinct(
        data$study_no[data$protein_strategy == "Substitution"]
      )
    )
  })
  list(
    joined_by_outcome = joined,
    principal_by_outcome = principal,
    strict_by_outcome = strict,
    audit = list(
      strategy_counts = audit |>
        dplyr::count(primary_strategy, review_status, confidence),
      outcome_summary = outcome_summary,
      support = support,
      blocking_studies = blocking,
      principal_studies = nrow(principal_audit),
      strict_studies = sum(audit$strict_eligible),
      audit_ready = nrow(blocking) == 0L
    )
  )
}

insert_strategy_formula_term <- function(formula, matched_term, new_term) {
  location <- regexpr(matched_term, formula, fixed = TRUE)
  if (location[[1L]] < 1L) {
    stop("Could not locate a strategy formula insertion point.", call. = FALSE)
  }
  end <- location[[1L]] + attr(location, "match.length") - 1L
  paste0(
    substr(formula, 1L, end),
    " +\n                         ", new_term,
    substr(formula, end + 1L, nchar(formula))
  )
}

first_strategy_dose_term <- function(formula, model_index) {
  if (model_index %in% c(1L, 5L)) return(NULL)
  patterns <- c(
    "splines::ns\\(\\s*protein_intake(?:_05)?\\s*,\\s*df\\s*=\\s*3\\s*\\)",
    paste0(
      "s\\(\\s*protein_intake(?:_05)?\\s*,[^)]*",
      "bs\\s*=\\s*['\"]tp['\"][^)]*\\)"
    ),
    "lspline\\(\\s*protein_intake(?:_05)?\\s*,\\s*0\\.5\\s*\\)",
    "\\bprotein_intake(?:_05)?\\b"
  )
  for (pattern in patterns) {
    location <- regexpr(pattern, formula, perl = TRUE)
    if (location[[1L]] > 0L) return(regmatches(formula, location))
  }
  stop("Could not identify the fixed protein term.", call. = FALSE)
}

make_protein_strategy_formula <- function(specification) {
  formula <- unname(as.character(specification$formula))
  if (length(formula) != 1L || !nzchar(formula)) {
    stop("A candidate formula is required.", call. = FALSE)
  }
  formula <- insert_strategy_formula_term(
    formula, "resistance_exercise_code", "strategy_code"
  )
  dose_term <- first_strategy_dose_term(formula, specification$model_index)
  if (is.null(dose_term)) return(formula)
  if (grepl("^s\\(", dose_term)) {
    variable <- regmatches(
      dose_term,
      regexpr("protein_intake(?:_05)?", dose_term, perl = TRUE)
    )
    interaction <- paste0(
      "s(", variable,
      ", by=strategy_code, bs='tp', k=4)"
    )
  } else {
    interaction <- paste0("strategy_code:", dose_term)
  }
  formula <- insert_strategy_formula_term(formula, dose_term, interaction)
  if (
    grepl("strategy_code.*resistance_exercise_code|resistance_exercise_code.*strategy_code.*:",
          interaction)
  ) {
    stop("A strategy-by-RT three-way term was introduced.", call. = FALSE)
  }
  formula
}

audit_protein_strategy_full_formula_registry <- function(specifications) {
  registry <- purrr::map_dfr(specifications, function(specification) {
    formula <- make_protein_strategy_formula(specification)
    includes_protein <- !specification$model_index %in% c(1L, 5L)
    tibble::tibble(
      outcome = specification$outcome,
      model_index = as.integer(specification$model_index),
      model_name = specification$model_name,
      formula = formula,
      has_strategy_main = grepl("\\bstrategy_code\\b", formula),
      has_strategy_dose_interaction = grepl(
        "strategy_code:|by\\s*=\\s*strategy_code", formula
      ),
      expected_strategy_dose_interaction = includes_protein,
      malformed_strategy_smooth = grepl(
        "s\\([^)]*strategy_code:", formula
      ),
      formula_parses = !inherits(
        try(stats::as.formula(formula), silent = TRUE), "try-error"
      )
    )
  }) |>
    dplyr::arrange(
      factor(
        outcome,
        levels = c("whole_body", "hypertrophy", "strength", "function")
      ),
      model_index
    )
  outcome_counts <- registry |>
    dplyr::count(outcome, name = "models")
  if (
    nrow(registry) != 104L ||
      anyDuplicated(paste(registry$outcome, registry$model_index)) ||
      nrow(outcome_counts) != 4L || any(outcome_counts$models != 26L) ||
      !all(registry$has_strategy_main) ||
      !identical(
        registry$has_strategy_dose_interaction,
        registry$expected_strategy_dose_interaction
      ) ||
      any(registry$malformed_strategy_smooth) ||
      !all(registry$formula_parses)
  ) {
    stop("The 104-model arm-strategy formula registry is invalid.",
         call. = FALSE)
  }
  list(
    models = nrow(registry),
    outcomes = nrow(outcome_counts),
    all_formulas_pass = TRUE,
    outcome_summary = outcome_counts,
    formula_registry = registry
  )
}

prepare_protein_strategy_pilot_input <- function(data, specification, priors) {
  base <- prepare_sensitivity_final_candidate(specification, priors, data)
  formula <- make_protein_strategy_formula(specification)
  original_formula <- specification$formula
  specification$primary_formula <- original_formula
  specification$formula <- formula
  specification$sensitivity_strategy_interaction <- TRUE
  base$family <- "strategy"
  base$formula <- formula
  base$specification <- specification
  base$strategy_reference <- "Substitution"
  base
}

make_protein_strategy_remedy_config <- function(
    base_config,
    remedy_model_indices = 1:4,
    iter_warmup = 1500L,
    iter_sampling = 3000L
) {
  remedy_model_indices <- unique(as.integer(remedy_model_indices))
  comparison_indices <- sort(unique(c(
    base_config$null_model_indices,
    base_config$model_indices
  )))
  if (
    length(remedy_model_indices) < 1L ||
      !all(remedy_model_indices %in% comparison_indices) ||
      iter_warmup < 1L || iter_sampling < 1L
  ) {
    stop("Invalid protein-strategy remedy configuration.", call. = FALSE)
  }
  config <- base_config
  config$remedy_model_indices <- remedy_model_indices
  config$iter_warmup <- as.integer(iter_warmup)
  config$iter_sampling <- as.integer(iter_sampling)
  config$seed <- as.integer(base_config$seed + 700L)
  config
}

audit_protein_strategy_pilot_structure <- function(
    strategy_data,
    inputs,
    expected_indices = c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
) {
  if (!is.list(strategy_data) || !is.list(inputs)) {
    stop("The strategy pilot inputs are incomplete.", call. = FALSE)
  }
  indices <- unname(purrr::map_int(
    inputs,
    ~ .x$specification$model_index
  ))
  formula_summary <- purrr::map_dfr(inputs, function(input) {
    includes_protein <- !input$specification$model_index %in% c(1L, 5L)
    has_strategy_dose <- grepl(
      "strategy_code:|by\\s*=\\s*strategy_code",
      input$formula
    )
    tibble::tibble(
      model_index = input$specification$model_index,
      model_name = input$specification$model_name,
      has_strategy_main = grepl("\\bstrategy_code\\b", input$formula),
      has_strategy_dose_interaction = has_strategy_dose,
      expected_strategy_dose_interaction = includes_protein,
      malformed_strategy_smooth = grepl(
        "s\\([^)]*strategy_code:", input$formula
      ),
      formula_parses = !inherits(
        try(stats::as.formula(input$formula), silent = TRUE),
        "try-error"
      ),
      has_strategy_rt_three_way = grepl(
        "strategy_code:resistance_exercise_code:|resistance_exercise_code:strategy_code:",
        gsub("\\s+", "", input$formula)
      )
    )
  })
  checks <- c(
    candidate_indices = identical(
      sort(indices), sort(as.integer(expected_indices))
    ),
    unique_indices = !anyDuplicated(indices),
    strategy_main = all(formula_summary$has_strategy_main),
    strategy_dose_interactions = identical(
        formula_summary$has_strategy_dose_interaction,
        formula_summary$expected_strategy_dose_interaction
      ),
    smooth_structure = !any(formula_summary$malformed_strategy_smooth),
    formula_parsing = all(formula_summary$formula_parses),
    no_strategy_rt_three_way = !any(
      formula_summary$has_strategy_rt_three_way
    ),
    family_labels = all(purrr::map_lgl(
      inputs, ~ identical(.x$family, "strategy")
    ))
  )
  if (!all(checks)) {
    stop(
      "The protein-strategy pilot failed: ",
      paste(names(checks)[!checks], collapse = ", "), ".",
      call. = FALSE
    )
  }
  list(
    ready_to_sample = isTRUE(strategy_data$audit$audit_ready),
    analysis = "principal resolved addition/substitution studies",
    reference_strategy = "Substitution",
    candidate_models = length(inputs),
    formula_summary = formula_summary,
    strategy_counts = strategy_data$audit$strategy_counts,
    outcome_summary = strategy_data$audit$outcome_summary,
    support = strategy_data$audit$support,
    blocking_studies = strategy_data$audit$blocking_studies,
    note = if (isTRUE(strategy_data$audit$audit_ready)) {
      "The classification gate is open."
    } else {
      paste0(
        nrow(strategy_data$audit$blocking_studies),
        " resolved addition/substitution studies still require confirmation."
      )
    }
  )
}

check_protein_strategy_sampling_gate <- function(
    enabled,
    structure_audit,
    family = "strategy"
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip(family))
  if (!isTRUE(structure_audit$ready_to_sample)) {
    stop(
      "The protein-strategy classification audit is not complete. ",
      "Resolve the studies listed in ",
      "sensitivity_strategy_pilot_structure_audit$blocking_studies first.",
      call. = FALSE
    )
  }
  TRUE
}

fit_protein_strategy_pilot <- function(input, sampling_gate, config) {
  if (is_sensitivity_pilot_skip(sampling_gate)) return(sampling_gate)
  if (!isTRUE(sampling_gate)) {
    stop("The protein-strategy sampling gate is invalid.", call. = FALSE)
  }
  fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = config,
    family = "strategy"
  )
}

new_protein_strategy_no_exact_loo_job <- function(family) {
  structure(
    list(family = family, no_exact_refit_required = TRUE),
    class = c("protein_strategy_no_exact_loo_job", "list")
  )
}

is_protein_strategy_no_exact_loo_job <- function(x) {
  inherits(x, "protein_strategy_no_exact_loo_job") ||
    isTRUE(x$no_exact_refit_required)
}

select_protein_strategy_exact_loo_jobs <- function(
    integrated_results,
    fitted_candidates,
    family
) {
  jobs <- select_sensitivity_mundlak_exact_loo_jobs(
    integrated_results = integrated_results,
    fitted_candidates = fitted_candidates,
    family = family
  )
  if (is_sensitivity_pilot_skip(jobs) || length(jobs) > 0L) return(jobs)
  # {targets} cannot dynamically branch over an empty list. A single inert
  # branch preserves the branched workflow and is removed by the collector.
  list(new_protein_strategy_no_exact_loo_job(family))
}

fit_protein_strategy_exact_loo_case <- function(
    job,
    enabled,
    config,
    family = "strategy"
) {
  if (is_protein_strategy_no_exact_loo_job(job)) return(job)
  refit <- fit_sensitivity_mundlak_exact_loo_case(
    job = job,
    enabled = enabled,
    config = config,
    family = family
  )
  if (is_sensitivity_pilot_skip(refit)) return(refit)
  refit$family <- family
  refit
}

compute_protein_strategy_exact_loo_case <- function(refit) {
  if (is_protein_strategy_no_exact_loo_job(refit)) return(refit)
  if (is_sensitivity_pilot_skip(refit)) return(refit)
  if (!refit$family %in% c("strategy", "strategy_within")) {
    stop("The exact-LOO refit is not a strategy candidate.", call. = FALSE)
  }
  if (identical(refit$family, "strategy_within")) {
    return(compute_sensitivity_mundlak_exact_loo_case(refit))
  }
  # Strategy models retain the primary study-effect design. Reuse the same
  # reconstruction path as the achieved/assignment families without changing
  # the shared exact-LOO implementation (and hence without invalidating their
  # completed targets).
  generic_refit <- refit
  generic_refit$family <- "achieved"
  result <- compute_sensitivity_mundlak_exact_loo_case(generic_refit)
  result$family <- "strategy"
  result
}

collect_protein_strategy_exact_results <- function(results, jobs, family) {
  if (is_sensitivity_pilot_skip(jobs)) {
    return(new_sensitivity_pilot_skip(family))
  }
  no_refit_jobs <- vapply(
    jobs,
    is_protein_strategy_no_exact_loo_job,
    logical(1L)
  )
  if (any(no_refit_jobs)) {
    if (
      length(jobs) != 1L || length(results) != 1L ||
        !isTRUE(no_refit_jobs[[1L]]) ||
        !is_protein_strategy_no_exact_loo_job(results[[1L]])
    ) {
      stop(
        paste0("The ", family, " no-refit placeholder is invalid."),
        call. = FALSE
      )
    }
    return(list())
  }
  collect_sensitivity_mundlak_exact_results(
    results = results,
    jobs = jobs,
    family = family
  )
}

prepare_protein_strategy_within_data <- function(
    contrast_data,
    audit,
    analysis = "principal"
) {
  contrast_data |>
    join_protein_strategy_audit(audit) |>
    finalise_strategy_analysis_data(analysis = analysis)
}

make_protein_strategy_within_formula <- function(
    formula,
    model_index,
    contrast_names
) {
  formula <- insert_strategy_formula_term(
    formula, "resistance_exercise_code", "strategy_code"
  )
  includes_protein <- !model_index %in% c(1L, 5L)
  if (!includes_protein) {
    if (!is.null(contrast_names)) {
      stop("A null within-strategy model has contrast columns.", call. = FALSE)
    }
    return(formula)
  }
  if (is.null(contrast_names) || length(contrast_names) < 1L) {
    stop("A protein within-strategy model lacks contrast columns.", call. = FALSE)
  }
  for (contrast_name in contrast_names) {
    formula <- insert_strategy_formula_term(
      formula,
      contrast_name,
      paste0(contrast_name, ":strategy_code")
    )
  }
  formula
}

prepare_protein_strategy_within_pilot_input <- function(
    data,
    specification,
    priors
) {
  base <- prepare_within_contrast_pilot_input(data, specification, priors)
  formula <- make_protein_strategy_within_formula(
    base$formula,
    specification$model_index,
    base$contrast_names
  )
  specification$within_formula <- base$formula
  specification$formula <- formula
  specification$sensitivity_strategy_interaction <- TRUE
  specification$sensitivity_within_contrast <- TRUE
  base$family <- "strategy_within"
  base$formula <- formula
  base$specification <- specification
  base$strategy_reference <- "Substitution"
  base
}

summarise_strategy_within_support <- function(
    data,
    outcome = "whole_body",
    analysis = "principal"
) {
  contrasts <- data |>
    dplyr::distinct(
      study_no, protein_arm_no, comparator_arm_no, protein_strategy,
      resistance_exercise_code, protein_intake_raw_treatment,
      protein_intake_raw_comparator
    )
  contrasts |>
    dplyr::group_by(protein_strategy, resistance_exercise_code) |>
    dplyr::summarise(
      contrasts = dplyr::n(),
      studies = dplyr::n_distinct(study_no),
      treatment_observed_minimum = min(protein_intake_raw_treatment),
      treatment_central_90_minimum = as.numeric(stats::quantile(
        protein_intake_raw_treatment, 0.05
      )),
      treatment_median = stats::median(protein_intake_raw_treatment),
      treatment_central_90_maximum = as.numeric(stats::quantile(
        protein_intake_raw_treatment, 0.95
      )),
      treatment_observed_maximum = max(protein_intake_raw_treatment),
      comparator_median = stats::median(protein_intake_raw_comparator),
      comparator_observed_minimum = min(protein_intake_raw_comparator),
      comparator_observed_maximum = max(protein_intake_raw_comparator),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      outcome = outcome,
      analysis = analysis,
      resistance_exercise = ifelse(
        resistance_exercise_code == 0, "YES", "NO"
      ),
      .before = 1L
    )
}

audit_protein_strategy_within_pilot_structure <- function(
    strategy_data,
    contrast_data,
    inputs,
    expected_indices = c(1L, 2L, 3L, 4L, 5L, 13L, 26L)
) {
  indices <- unname(purrr::map_int(
    inputs,
    ~ .x$specification$model_index
  ))
  formula_summary <- purrr::map_dfr(inputs, function(input) {
    includes_protein <- !input$specification$model_index %in% c(1L, 5L)
    contrast_names <- input$contrast_names
    has_each_strategy_interaction <- if (includes_protein) {
      length(contrast_names) > 0L && all(vapply(
        contrast_names,
        function(name) grepl(
          paste0(name, ":strategy_code"), input$formula, fixed = TRUE
        ),
        logical(1)
      ))
    } else {
      FALSE
    }
    tibble::tibble(
      model_index = input$specification$model_index,
      model_name = input$specification$model_name,
      contrast_basis_dimension = if (is.null(contrast_names)) {
        0L
      } else {
        length(contrast_names)
      },
      contrast_basis_valid = if (includes_protein) {
        length(contrast_names) > 0L
      } else {
        is.null(contrast_names)
      },
      has_strategy_main = grepl("\\bstrategy_code\\b", input$formula),
      has_strategy_contrast_interaction = has_each_strategy_interaction,
      expected_strategy_contrast_interaction = includes_protein,
      formula_parses = !inherits(
        try(stats::as.formula(input$formula), silent = TRUE),
        "try-error"
      ),
      has_strategy_rt_three_way = grepl(
        "strategy_code:resistance_exercise_code:|resistance_exercise_code:strategy_code:",
        gsub("\\s+", "", input$formula)
      )
    )
  })
  support <- summarise_strategy_within_support(contrast_data)
  study_counts <- contrast_data |>
    dplyr::distinct(study_no, protein_strategy) |>
    dplyr::count(protein_strategy, name = "studies")
  checks <- c(
    candidate_indices = identical(
      sort(indices), sort(as.integer(expected_indices))
    ),
    unique_indices = !anyDuplicated(indices),
    strategy_main = all(formula_summary$has_strategy_main),
    contrast_basis = all(formula_summary$contrast_basis_valid),
    strategy_contrast_interactions = identical(
        formula_summary$has_strategy_contrast_interaction,
        formula_summary$expected_strategy_contrast_interaction
      ),
    formula_parsing = all(formula_summary$formula_parses),
    no_strategy_rt_three_way = !any(
      formula_summary$has_strategy_rt_three_way
    ),
    family_labels = all(purrr::map_lgl(
      inputs, ~ identical(.x$family, "strategy_within")
    )),
    both_strategy_groups = setequal(
      study_counts$protein_strategy,
      c("Addition", "Substitution")
    )
  )
  if (!all(checks)) {
    stop(
      "The within-study strategy pilot failed: ",
      paste(names(checks)[!checks], collapse = ", "), ".",
      call. = FALSE
    )
  }
  list(
    ready_to_sample = isTRUE(strategy_data$audit$audit_ready),
    analysis = paste(
      "strict within-study protein-arm minus comparator-arm contrasts",
      "among resolved addition/substitution studies"
    ),
    reference_strategy = "Substitution",
    candidate_models = length(inputs),
    contrast_rows = nrow(contrast_data),
    study_counts = study_counts,
    support = support,
    formula_summary = formula_summary,
    blocking_studies = strategy_data$audit$blocking_studies,
    note = if (isTRUE(strategy_data$audit$audit_ready)) {
      "The classification gate is open."
    } else {
      paste0(
        nrow(strategy_data$audit$blocking_studies),
        " resolved addition/substitution studies still require confirmation."
      )
    }
  )
}

fit_protein_strategy_within_pilot <- function(
    input,
    sampling_gate,
    config
) {
  if (is_sensitivity_pilot_skip(sampling_gate)) return(sampling_gate)
  if (!isTRUE(sampling_gate)) {
    stop("The within-strategy sampling gate is invalid.", call. = FALSE)
  }
  fit_sensitivity_brms_pilot(
    input = input,
    enabled = TRUE,
    config = config,
    family = "strategy_within"
  )
}

make_protein_strategy_full_config <- function(
    base_config,
    seed_offset,
    iter_warmup = 1500L,
    iter_sampling = 2000L
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
      config$iter_warmup != 1500L || config$iter_sampling != 2000L ||
      length(config$model_indices) != 26L ||
      anyDuplicated(config$model_indices)
  ) {
    stop("Invalid full protein-strategy expansion configuration.",
         call. = FALSE)
  }
  config
}

make_protein_strategy_full_remedy_config <- function(
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
    stop("Invalid full protein-strategy remedy configuration.",
         call. = FALSE)
  }
  config
}

select_pending_protein_strategy_full_specifications <- function(
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
    !family %in% c("strategy", "strategy_within") ||
      length(specifications) != 104L || anyDuplicated(specification_keys) ||
      nrow(reused_fit_manifest) != 7L ||
      any(reused_fit_manifest$family != family) ||
      !setequal(reused_keys, expected_reused) ||
      !all(reused_keys %in% specification_keys) ||
      any(!file.exists(reused_fit_manifest$fit_path))
  ) {
    stop("The full protein-strategy specification/reuse set is invalid.",
         call. = FALSE)
  }
  specifications[!specification_keys %in% reused_keys]
}

prepare_protein_strategy_full_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config,
    family
) {
  outcome <- specification$outcome
  data <- data_by_outcome[[outcome]]
  priors <- priors_by_outcome[[outcome]]
  if (
    !family %in% c("strategy", "strategy_within") ||
      is.null(data) || is.null(priors) || !outcome %in% config$outcomes
  ) {
    stop("A full protein-strategy outcome input is unavailable.",
         call. = FALSE)
  }
  specification$sensitivity_strategy_expansion <- TRUE
  specification$sensitivity_strategy_family <- family
  specification$sensitivity_strategy_sampling <- paste0(
    config$iter_warmup, "_warmup_", config$iter_sampling, "_sampling"
  )
  if (identical(family, "strategy")) {
    prepare_protein_strategy_pilot_input(data, specification, priors)
  } else {
    prepare_protein_strategy_within_pilot_input(
      data, specification, priors
    )
  }
}

select_protein_strategy_full_remedy_specifications <- function(
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
    !family %in% c("strategy", "strategy_within") ||
      !identical(initial_audit$family, family) ||
      nrow(failed) < 1L || anyDuplicated(failed$key) ||
      anyDuplicated(specification_keys) ||
      length(selected) != nrow(failed) ||
      !setequal(selected_keys, failed$key)
  ) {
    stop("The full protein-strategy remedy set is invalid.",
         call. = FALSE)
  }
  selected
}

prepare_protein_strategy_full_remedy_input <- function(
    specification,
    data_by_outcome,
    priors_by_outcome,
    config,
    family
) {
  specification$adapt_delta <- max(
    specification$adapt_delta,
    config$adapt_delta_floor
  )
  specification$max_treedepth <- max(
    specification$max_treedepth,
    config$max_treedepth_floor
  )
  specification$sensitivity_strategy_full_remedy <-
    "longer_sampling_high_adapt_delta"
  prepare_protein_strategy_full_input(
    specification = specification,
    data_by_outcome = data_by_outcome,
    priors_by_outcome = priors_by_outcome,
    config = config,
    family = family
  )
}

make_protein_strategy_within_full_data_by_outcome <- function(
    primary_data_by_outcome,
    audit
) {
  contrast_data <- make_within_full_contrast_data_by_outcome(
    primary_data_by_outcome
  )
  purrr::map(
    contrast_data,
    prepare_protein_strategy_within_data,
    audit = audit,
    analysis = "principal"
  )
}

audit_protein_strategy_full_expansion <- function(
    pending_specifications,
    inputs,
    data_by_outcome,
    reused_fit_manifest,
    config,
    family
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
    model_index <- as.integer(input$specification$model_index)
    includes_protein <- !model_index %in% c(1L, 5L)
    has_strategy_interaction <- if (!includes_protein) {
      TRUE
    } else if (identical(family, "strategy")) {
      grepl("strategy_code:|by\\s*=\\s*strategy_code", input$formula)
    } else {
      !is.null(input$contrast_names) &&
        length(input$contrast_names) > 0L &&
        all(vapply(
          input$contrast_names,
          function(name) grepl(
            paste0(name, ":strategy_code"),
            input$formula,
            fixed = TRUE
          ),
          logical(1L)
        ))
    }
    identical(input$family, family) &&
      is.character(input$formula) && length(input$formula) == 1L &&
      !inherits(
        try(stats::as.formula(input$formula), silent = TRUE),
        "try-error"
      ) &&
      grepl("\\bstrategy_code\\b", input$formula) &&
      has_strategy_interaction &&
      is.data.frame(input$data) && nrow(input$data) > 0L &&
      "strategy_code" %in% names(input$data) &&
      all(input$data$strategy_code %in% c(0, 1))
  })
  data_summary <- purrr::imap_dfr(data_by_outcome, function(data, outcome) {
    tibble::tibble(
      outcome = outcome,
      rows = nrow(data),
      studies = dplyr::n_distinct(data$study_no),
      addition_studies = dplyr::n_distinct(
        data$study_no[data$protein_strategy == "Addition"]
      ),
      substitution_studies = dplyr::n_distinct(
        data$study_no[data$protein_strategy == "Substitution"]
      )
    )
  })
  if (
    !family %in% c("strategy", "strategy_within") ||
      length(pending_specifications) != 97L || length(inputs) != 97L ||
      nrow(reused_fit_manifest) != 7L || anyDuplicated(pending_keys) ||
      anyDuplicated(input_keys) || !setequal(pending_keys, input_keys) ||
      !identical(as.integer(observed_counts), unname(expected_counts)) ||
      !all(valid_inputs) || nrow(data_summary) != 4L ||
      any(data_summary$rows < 1L) ||
      any(data_summary$addition_studies < 1L) ||
      any(data_summary$substitution_studies < 1L) ||
      config$iter_warmup != 1500L || config$iter_sampling != 2000L
  ) {
    stop("The full protein-strategy expansion failed its structural audit.",
         call. = FALSE)
  }
  list(
    family = family,
    checks = tibble::tibble(
      check = c(
        "complete candidate registry",
        "reused passing whole-body pilots",
        "new independently checkpointed fits",
        "pilot-informed sampling schedule",
        "resolved addition/substitution analysis data"
      ),
      status = "PASS",
      detail = c(
        "104 candidates: 26 models x 4 outcomes",
        "7 whole-body pilot fits",
        "97 fits: 19 whole-body plus 26 for each other outcome",
        "1,500 warmup and 2,000 retained iterations per chain",
        if (identical(family, "strategy")) {
          "arm-based strategy-modified candidate set"
        } else {
          "strategy-modified strict within-study contrasts"
        }
      )
    ),
    data_summary = data_summary,
    pending_fits = length(inputs),
    reused_fits = nrow(reused_fit_manifest),
    structure_pass = TRUE
  )
}

fit_protein_strategy_full_candidate_to_file <- function(
    input,
    enabled,
    config,
    family,
    stage = "initial"
) {
  if (!isTRUE(enabled)) return(new_sensitivity_pilot_skip(family))
  if (!family %in% c("strategy", "strategy_within")) {
    stop("Unknown full protein-strategy family.", call. = FALSE)
  }
  outcome_offset <- match(input$specification$outcome, config$outcomes)
  if (is.na(outcome_offset) || !identical(input$family, family)) {
    stop("Unknown outcome or family in the full strategy fit.",
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
    family = family
  )
  write_sensitivity_fit_file(candidate, stage = stage)
}

combine_protein_strategy_full_fit_manifest <- function(
    reused_manifest,
    new_manifest,
    family
) {
  outcomes <- c("whole_body", "hypertrophy", "strength", "function")
  manifest <- dplyr::bind_rows(reused_manifest, new_manifest) |>
    dplyr::arrange(factor(outcome, levels = outcomes), model_index)
  counts <- table(factor(manifest$outcome, levels = outcomes))
  if (
    !family %in% c("strategy", "strategy_within") ||
      nrow(manifest) != 104L ||
      anyDuplicated(paste(manifest$outcome, manifest$model_index)) ||
      !identical(as.integer(counts), rep(26L, 4L)) ||
      any(manifest$family != family) ||
      any(!file.exists(manifest$fit_path))
  ) {
    stop("The full protein-strategy fit manifest is incomplete.",
         call. = FALSE)
  }
  manifest
}

make_protein_strategy_full_fit_records <- function(
    fit_manifest,
    family
) {
  if (
    !family %in% c("strategy", "strategy_within") ||
      !is.data.frame(fit_manifest) || nrow(fit_manifest) != 104L ||
      any(fit_manifest$family != family) ||
      any(!file.exists(fit_manifest$fit_path))
  ) {
    stop("The full protein-strategy fit manifest is invalid.",
         call. = FALSE)
  }
  purrr::map(seq_len(nrow(fit_manifest)), function(row) {
    purrr::map(as.list(fit_manifest[row, , drop = FALSE]), unname)
  })
}

compute_protein_strategy_full_file_integrated_loo <- function(
    fit_record,
    family,
    loo_draws = 4000L,
    k_threshold = 0.7,
    cores = 1L
) {
  if (!family %in% c("strategy", "strategy_within")) {
    stop("Unknown protein-strategy family for integrated LOO.",
         call. = FALSE)
  }
  if (identical(family, "strategy_within")) {
    return(compute_mundlak_full_file_integrated_loo(
      fit_record = fit_record,
      loo_draws = loo_draws,
      k_threshold = k_threshold,
      cores = cores,
      family = family,
      slope_description =
        "strategy_modified_within_contrast_basis_slopes"
    ))
  }
  fit_path <- fit_record$fit_path
  if (
    length(fit_path) != 1L || is.na(fit_path) || !file.exists(fit_path)
  ) {
    stop("The protein-strategy LOO fit file is unavailable.",
         call. = FALSE)
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
    stop("The strategy fit file does not match its manifest record.",
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

compute_protein_strategy_full_exact_loo_file <- function(
    path,
    diagnostic_audit,
    family
) {
  if (!isTRUE(diagnostic_audit$all_exact_refit_diagnostics_pass)) {
    stop("Exact strategy corrections require passing refit diagnostics.",
         call. = FALSE)
  }
  if (!family %in% c("strategy", "strategy_within")) {
    stop("Unknown protein-strategy exact-LOO family.", call. = FALSE)
  }
  refit <- read_sensitivity_full_exact_loo_file(path, family)
  refit$original_data <- refit$prediction_data
  refit$integrated_result <- qs::qread(refit$integrated_loo_path)
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = 1L)
  compute_protein_strategy_exact_loo_case(refit)
}

combine_protein_strategy_full_remedy_diagnostic_audit <- function(
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
    !family %in% c("strategy", "strategy_within") ||
      !identical(initial_audit$family, family) ||
      nrow(remedy_summary) != length(expected_keys) ||
      anyDuplicated(remedy_keys) || !setequal(remedy_keys, expected_keys)
  ) {
    stop("The full protein-strategy remedy diagnostics are incomplete.",
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
    stop("The remedied protein-strategy audit is incomplete.",
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

combine_protein_strategy_full_accepted_fit_manifest <- function(
    initial_manifest,
    remedy_manifest,
    final_diagnostic_audit,
    family
) {
  if (
    !family %in% c("strategy", "strategy_within") ||
      !identical(final_diagnostic_audit$family, family) ||
      !isTRUE(final_diagnostic_audit$all_sampling_diagnostics_pass)
  ) {
    stop("Protein-strategy fits cannot be accepted before diagnostics pass.",
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
    stop("The accepted protein-strategy fit manifest is incomplete.",
         call. = FALSE)
  }
  accepted
}

make_protein_strategy_prediction_grid <- function(
    fitted_candidate,
    data,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.02
) {
  base <- make_primary_prediction_grid(
    data, protein_min = protein_min, protein_max = protein_max,
    grid_step = grid_step
  )
  strategies <- c("Substitution", "Addition")
  newdata <- dplyr::bind_rows(lapply(strategies, function(strategy) {
    result <- base$newdata
    result$protein_strategy <- strategy
    result$strategy_code <- ifelse(strategy == "Addition", 1, 0)
    result$strategy_factor <- factor(
      strategy, levels = c("Substitution", "Addition")
    )
    result
  }))
  keys <- dplyr::bind_rows(lapply(strategies, function(strategy) {
    dplyr::mutate(
      base$keys,
      protein_strategy = strategy,
      strategy_code = ifelse(strategy == "Addition", 1, 0),
      .before = protein_intake
    )
  }))
  list(
    newdata = newdata,
    keys = keys,
    protein_median = base$protein_median,
    protein_min = base$protein_min,
    protein_max = base$protein_max,
    grid_step = base$grid_step
  )
}

make_protein_strategy_within_prediction_grid <- function(
    fitted_candidate,
    data,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.02,
    reference_raw = NULL
) {
  base <- make_within_prediction_newdata(
    fitted_candidate = fitted_candidate,
    data = data,
    protein_min = protein_min,
    protein_max = protein_max,
    grid_step = grid_step,
    reference_raw = reference_raw
  )
  strategies <- c("Substitution", "Addition")
  newdata <- dplyr::bind_rows(lapply(strategies, function(strategy) {
    result <- base$newdata
    result$protein_strategy <- strategy
    result$strategy_code <- ifelse(strategy == "Addition", 1, 0)
    result$strategy_factor <- factor(
      strategy, levels = c("Substitution", "Addition")
    )
    result
  }))
  keys <- dplyr::bind_rows(lapply(strategies, function(strategy) {
    dplyr::mutate(
      base$keys,
      protein_strategy = strategy,
      strategy_code = ifelse(strategy == "Addition", 1, 0),
      .before = protein_intake
    )
  }))
  list(
    newdata = newdata,
    keys = keys,
    protein_min = base$protein_min,
    protein_max = base$protein_max,
    grid_step = base$grid_step
  )
}

build_protein_strategy_pilot_reporting <- function(
    model_average,
    data,
    conf_level = 0.95,
    outcome = "whole_body"
) {
  curve <- summarise_primary_dose_response(model_average)
  support <- summarise_strategy_support(
    data, outcome = outcome, analysis = "principal"
  )
  common_support <- support |>
    dplyr::group_by(resistance_exercise_code) |>
    dplyr::summarise(
      common_observed_minimum = max(observed_minimum),
      common_observed_maximum = min(observed_maximum),
      common_central_90_minimum = max(central_90_minimum),
      common_central_90_maximum = min(central_90_maximum),
      .groups = "drop"
    )
  curve <- curve |>
    dplyr::left_join(
      dplyr::select(
        support, protein_strategy, resistance_exercise_code,
        observed_minimum, observed_maximum,
        central_90_minimum, central_90_maximum
      ),
      by = c("protein_strategy", "resistance_exercise_code")
    ) |>
    dplyr::left_join(common_support, by = "resistance_exercise_code") |>
    dplyr::mutate(
      within_strategy_observed_support =
        protein_intake >= observed_minimum &
        protein_intake <= observed_maximum,
      within_common_observed_support =
        protein_intake >= common_observed_minimum &
        protein_intake <= common_observed_maximum,
      within_common_central_90_support =
        protein_intake >= common_central_90_minimum &
        protein_intake <= common_central_90_maximum
    )
  grid <- model_average$grid
  alpha <- (1 - conf_level) / 2
  comparison_keys <- grid |>
    dplyr::distinct(resistance_exercise_code, resistance_exercise,
                    protein_intake) |>
    dplyr::arrange(resistance_exercise_code, protein_intake)
  differences <- purrr::pmap_dfr(comparison_keys, function(
      resistance_exercise_code, resistance_exercise, protein_intake
  ) {
    addition <- which(
      grid$resistance_exercise_code == resistance_exercise_code &
        grid$protein_strategy == "Addition" &
        abs(grid$protein_intake - protein_intake) < 1e-10
    )
    substitution <- which(
      grid$resistance_exercise_code == resistance_exercise_code &
        grid$protein_strategy == "Substitution" &
        abs(grid$protein_intake - protein_intake) < 1e-10
    )
    if (length(addition) != 1L || length(substitution) != 1L) {
      stop("Strategy prediction columns do not pair uniquely.", call. = FALSE)
    }
    draws <- model_average$draws[, addition] -
      model_average$draws[, substitution]
    interval <- stats::quantile(draws, c(alpha, 1 - alpha))
    tibble::tibble(
      resistance_exercise_code = resistance_exercise_code,
      resistance_exercise = resistance_exercise,
      protein_intake = protein_intake,
      addition_minus_substitution = mean(draws),
      conf.low = interval[[1L]],
      conf.high = interval[[2L]],
      probability_positive = mean(draws > 0)
    )
  }) |>
    dplyr::left_join(common_support, by = "resistance_exercise_code") |>
    dplyr::mutate(
      within_common_observed_support =
        protein_intake >= common_observed_minimum &
        protein_intake <= common_observed_maximum,
      within_common_central_90_support =
        protein_intake >= common_central_90_minimum &
        protein_intake <= common_central_90_maximum
    )
  list(
    dose_response = curve,
    addition_minus_substitution = differences,
    support = support,
    common_support = common_support,
    model_weights = model_average$model_weights,
    settings = model_average$settings
  )
}

build_protein_strategy_within_pilot_reporting <- function(
    model_average,
    data,
    conf_level = 0.95,
    outcome = "whole_body"
) {
  curve <- summarise_primary_dose_response(model_average)
  support <- summarise_strategy_within_support(data, outcome = outcome)
  common_support <- support |>
    dplyr::group_by(resistance_exercise_code) |>
    dplyr::summarise(
      common_observed_minimum = max(treatment_observed_minimum),
      common_observed_maximum = min(treatment_observed_maximum),
      common_central_90_minimum = max(treatment_central_90_minimum),
      common_central_90_maximum = min(treatment_central_90_maximum),
      .groups = "drop"
    )
  curve <- curve |>
    dplyr::left_join(
      dplyr::select(
        support, protein_strategy, resistance_exercise_code,
        treatment_observed_minimum, treatment_observed_maximum,
        treatment_central_90_minimum, treatment_central_90_maximum
      ),
      by = c("protein_strategy", "resistance_exercise_code")
    ) |>
    dplyr::left_join(common_support, by = "resistance_exercise_code") |>
    dplyr::mutate(
      within_strategy_observed_support =
        protein_intake >= treatment_observed_minimum &
        protein_intake <= treatment_observed_maximum,
      within_common_observed_support =
        protein_intake >= common_observed_minimum &
        protein_intake <= common_observed_maximum,
      within_common_central_90_support =
        protein_intake >= common_central_90_minimum &
        protein_intake <= common_central_90_maximum
    )
  grid <- model_average$grid
  alpha <- (1 - conf_level) / 2
  comparison_keys <- grid |>
    dplyr::distinct(
      resistance_exercise_code, resistance_exercise, protein_intake
    ) |>
    dplyr::arrange(resistance_exercise_code, protein_intake)
  differences <- purrr::pmap_dfr(comparison_keys, function(
      resistance_exercise_code, resistance_exercise, protein_intake
  ) {
    addition <- which(
      grid$resistance_exercise_code == resistance_exercise_code &
        grid$protein_strategy == "Addition" &
        abs(grid$protein_intake - protein_intake) < 1e-10
    )
    substitution <- which(
      grid$resistance_exercise_code == resistance_exercise_code &
        grid$protein_strategy == "Substitution" &
        abs(grid$protein_intake - protein_intake) < 1e-10
    )
    if (length(addition) != 1L || length(substitution) != 1L) {
      stop("Within-strategy prediction columns do not pair uniquely.",
           call. = FALSE)
    }
    draws <- model_average$draws[, addition] -
      model_average$draws[, substitution]
    interval <- stats::quantile(draws, c(alpha, 1 - alpha))
    tibble::tibble(
      resistance_exercise_code = resistance_exercise_code,
      resistance_exercise = resistance_exercise,
      protein_intake = protein_intake,
      addition_minus_substitution = mean(draws),
      conf.low = interval[[1L]],
      conf.high = interval[[2L]],
      probability_positive = mean(draws > 0)
    )
  }) |>
    dplyr::left_join(common_support, by = "resistance_exercise_code") |>
    dplyr::mutate(
      within_common_observed_support =
        protein_intake >= common_observed_minimum &
        protein_intake <= common_observed_maximum,
      within_common_central_90_support =
        protein_intake >= common_central_90_minimum &
        protein_intake <= common_central_90_maximum
    )
  list(
    dose_response_contrasts = curve,
    addition_minus_substitution = differences,
    support = support,
    common_support = common_support,
    model_weights = model_average$model_weights,
    settings = model_average$settings
  )
}

make_protein_strategy_full_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    n_draws = 4000L,
    protein_min = 0.8,
    grid_step = 0.01,
    seed = 32988L
) {
  if (
    !is.data.frame(accepted_fit_manifest) ||
      nrow(accepted_fit_manifest) != 104L ||
      any(accepted_fit_manifest$family != "strategy") ||
      !is.data.frame(stacking_weights) ||
      nrow(stacking_weights) != 104L ||
      any(stacking_weights$family != "strategy")
  ) {
    stop("The full protein-strategy prediction inputs are invalid.",
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
    record$family <- "strategy"
    record$prediction_type <- "strategy"
    record
  })
}

read_protein_strategy_completed_checkpoint <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("The completed protein-strategy checkpoint is unavailable.",
         call. = FALSE)
  }
  checkpoint <- readRDS(path)
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
      !identical(checkpoint$family, "strategy") ||
      !isTRUE(checkpoint$audit$all_sampling_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$stacking_ready) ||
      nrow(checkpoint$accepted_fit_manifest) != 104L ||
      nrow(checkpoint$stacking_weights) != 104L
  ) {
    stop("The completed protein-strategy checkpoint is invalid.",
         call. = FALSE)
  }
  checkpoint
}

write_protein_strategy_full_prediction_draw_file <- function(
    record,
    checkpoint_path = NULL,
    data_by_outcome = NULL,
    result_path = NULL
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !identical(unname(as.character(record$family)), "strategy") ||
      !identical(unname(as.character(record$prediction_type)), "strategy")
  ) {
    stop("A full protein-strategy prediction record is invalid.",
         call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "strategy") ||
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
    stop("A protein-strategy prediction fit does not match its record.",
         call. = FALSE)
  }
  if (is.null(data_by_outcome)) {
    checkpoint <- read_protein_strategy_completed_checkpoint(checkpoint_path)
    data_by_outcome <- checkpoint$analysis_data_by_outcome
  }
  prediction_data <- data_by_outcome[[record$outcome]]
  if (is.null(prediction_data) || !is.data.frame(prediction_data)) {
    stop("The protein-strategy prediction data are unavailable.",
         call. = FALSE)
  }
  grid <- make_protein_strategy_prediction_grid(
    fitted_candidate = candidate,
    data = prediction_data,
    protein_min = record$protein_min,
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
    stop("Protein-strategy prediction draws do not match their allocation.",
         call. = FALSE)
  }
  reduction <- list(
    family = "strategy",
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
      prediction_type = "strategy"
    )
  )
  path <- if (is.null(result_path)) {
    sensitivity_reduction_file_path(
      family = "strategy",
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

subset_protein_strategy_model_average <- function(
    model_average,
    protein_strategy
) {
  if (
    !inherits(model_average, "primary_model_average") ||
      !protein_strategy %in% c("Substitution", "Addition") ||
      !"protein_strategy" %in% names(model_average$grid)
  ) {
    stop("The protein-strategy model average cannot be subset.",
         call. = FALSE)
  }
  columns <- which(model_average$grid$protein_strategy == protein_strategy)
  if (!length(columns)) {
    stop("The requested protein strategy is absent from the grid.",
         call. = FALSE)
  }
  output <- model_average
  output$grid <- model_average$grid[columns, , drop = FALSE]
  output$draws <- model_average$draws[, columns, drop = FALSE]
  output$settings$protein_strategy <- protein_strategy
  output
}

write_protein_strategy_full_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10,
    result_path = NULL,
    family = "strategy",
    protein_min = reference_intake
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path, family = family
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (
    is.null(data) || !is.data.frame(data) ||
      !all(c("protein_strategy", "protein_intake_raw") %in% names(data))
  ) {
    stop("The full protein-strategy reporting data are unavailable.",
         call. = FALSE)
  }
  strategy_summary <- build_protein_strategy_pilot_reporting(
    model_average = model_average,
    data = data,
    outcome = outcome
  )
  strategy_reports <- purrr::map(
    c("Substitution", "Addition"),
    function(strategy) {
      strategy_average <- subset_protein_strategy_model_average(
        model_average, strategy
      )
      strategy_data <- dplyr::filter(
        data, protein_strategy == strategy
      )
      absolute <- summarise_primary_dose_response(strategy_average)
      supported <- build_sensitivity_mundlak_reporting_pilot(
        absolute,
        strategy_average,
        strategy_data,
        reference_intake = reference_intake
      )
      heatmap <- compute_primary_heatmap_contrasts(
        model_average = strategy_average,
        data = strategy_data,
        protein_min = protein_min,
        sesoi = sesoi
      )
      annotate <- function(table) {
        dplyr::mutate(
          table,
          protein_strategy = strategy,
          .before = 1L
        )
      }
      list(
        increments = annotate(supported$increments_full_grid),
        support = annotate(supported$support_by_group),
        common_support = annotate(supported$common_support),
        pairwise = annotate(heatmap$contrasts),
        pairwise_support = annotate(heatmap$support),
        local_equivalence = annotate(
          summarise_primary_local_equivalence(
            heatmap$contrasts, outcome = outcome
          )
        )
      )
    }
  )
  bind_component <- function(name) {
    dplyr::bind_rows(purrr::map(strategy_reports, name))
  }
  report <- list(
    family = family,
    outcome = outcome,
    dose_response = strategy_summary$dose_response,
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

read_protein_strategy_full_reporting_file <- function(
    path,
    family = "strategy"
) {
  report <- read_sensitivity_full_reporting_file(path, family = family)
  if (
    !"addition_minus_substitution" %in% names(report) ||
      !is.data.frame(report$addition_minus_substitution) ||
      !all(c("Substitution", "Addition") %in%
           unique(report$dose_response$protein_strategy))
  ) {
    stop("A full protein-strategy reporting file is invalid.",
         call. = FALSE)
  }
  report
}

audit_protein_strategy_full_reporting_files <- function(
    paths,
    family = "strategy"
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != 4L || anyNA(paths) || any(!file.exists(paths))) {
    stop("The protein-strategy reporting-file set is incomplete.",
         call. = FALSE)
  }
  summary <- purrr::map_dfr(paths, function(path) {
    report <- read_protein_strategy_full_reporting_file(path, family)
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
    ) && all(
      is.finite(report$addition_minus_substitution$addition_minus_substitution) &
        report$addition_minus_substitution$conf.low <=
          report$addition_minus_substitution$addition_minus_substitution &
        report$addition_minus_substitution$addition_minus_substitution <=
          report$addition_minus_substitution$conf.high
    )
    reference_valid <-
      nrow(reference) == 4L &&
      all(abs(reference$increment) < 1e-12) &&
      all(abs(reference$increment_conf_low) < 1e-12) &&
      all(abs(reference$increment_conf_high) < 1e-12)
    probabilities <- as.matrix(
      report$pairwise_contrasts[, probability_columns, drop = FALSE]
    )
    strategy_probabilities <-
      report$addition_minus_substitution$probability_positive
    probabilities_valid <-
      all(is.finite(probabilities)) &&
      all(probabilities >= 0 & probabilities <= 1) &&
      all(is.finite(strategy_probabilities)) &&
      all(strategy_probabilities >= 0 & strategy_probabilities <= 1) &&
      all(abs(
        report$pairwise_contrasts$probability_meaningful_gain +
          report$pairwise_contrasts$probability_negligible +
          report$pairwise_contrasts$probability_meaningful_loss - 1
      ) < 1e-12)
    diagonal_valid <-
      nrow(diagonal) >= 4L &&
      all(abs(diagonal$estimate) < 1e-12) &&
      all(abs(diagonal$probability_negligible - 1) < 1e-12)
    tibble::tibble(
      family = report$family,
      outcome = report$outcome,
      dose_response_rows = nrow(report$dose_response),
      increment_rows = nrow(report$dose_response_increments),
      strategy_contrast_rows = nrow(report$addition_minus_substitution),
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
      any(summary$dose_response_rows < 4L) ||
      any(summary$strategy_contrast_rows < 2L) ||
      any(summary$pairwise_contrast_rows < 4L) ||
      any(abs(summary$weight_sum - 1) > 1e-8) ||
      !all(summary$intervals_valid) ||
      !all(summary$reference_increment_valid) ||
      !all(summary$probabilities_valid) ||
      !all(summary$diagonal_contrasts_valid)
  ) {
    stop("The full protein-strategy reporting audit failed.",
         call. = FALSE)
  }
  list(
    family = family,
    summary = summary,
    outcomes = nrow(summary),
    posterior_mixture_draws_per_outcome = 4000L,
    all_reporting_checks_pass = TRUE
  )
}

write_protein_strategy_full_verification_bundle <- function(
    reporting_audit,
    checkpoint,
    prediction_manifest
) {
  if (
    !identical(checkpoint$family, "strategy") ||
      !isTRUE(reporting_audit$all_reporting_checks_pass) ||
      !isTRUE(checkpoint$audit$all_sampling_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$stacking_ready)
  ) {
    stop("The protein-strategy family is not ready for verification.",
         call. = FALSE)
  }
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_protein_strategy_full_reporting_file
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
    family = "strategy",
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
      addition_minus_substitution =
        bind_report_table("addition_minus_substitution"),
      pairwise_contrasts = bind_report_table("pairwise_contrasts"),
      pairwise_support = bind_report_table("pairwise_support"),
      local_equivalence = bind_report_table("local_equivalence")
    ),
    descriptions = c(
      analysis_data = "Outcome-specific arm rows joined to the audited randomized protein-manipulation strategy.",
      model_specifications = "Final candidate-model registry inherited from the primary analysis.",
      accepted_fit_manifest = "Accepted one-file-per-model registry and fit checksums.",
      sampling_diagnostic_summary = "Final model-level convergence, geometry, and predictive diagnostics.",
      sampling_parameter_failures = "Any parameter-level diagnostic failures after remediation.",
      initial_integrated_loo_model_summary = "Integrated observation-level PSIS summaries before exact correction.",
      initial_integrated_loo_outcome_summary = "Outcome-level integrated PSIS audit before exact correction.",
      initial_high_pareto_k_cases = "Cases selected prospectively for exact case-deletion correction.",
      exact_refit_diagnostic_summary = "Sampling diagnostics for every exact case-deletion refit.",
      exact_loo_correction_summary = "Exact versus PSIS pointwise correction results.",
      final_loo_model_summary = "Final corrected integrated LOO summaries used for stacking.",
      stacking_weights = "Final outcome-specific stacking weights.",
      prediction_allocation_manifest = "Deterministic allocation of 4,000 mixture curves per outcome.",
      reporting_file_manifest = "Checksums and dimensions of reduced outcome reporting files.",
      dose_response_predictions = "Model-averaged arm-level dose-response curves under addition and substitution strategies.",
      dose_response_increments = "Strategy-specific paired increments from 0.8 g/kg/day.",
      exposure_support_by_group = "Observed protein-intake support by strategy and resistance-exercise condition.",
      common_exposure_support = "Protein-intake support shared by both strategies within resistance-exercise condition.",
      addition_minus_substitution = "Direct posterior addition-minus-substitution contrast at each dose and exercise condition.",
      pairwise_contrasts = "Strategy-specific paired dose contrasts and SESOI probabilities.",
      pairwise_support = "Observed support used to flag strategy-specific pairwise contrasts.",
      local_equivalence = "Strategy-specific local practical-equivalence summaries."
    )
  )
}

make_protein_strategy_within_full_prediction_records <- function(
    accepted_fit_manifest,
    stacking_weights,
    n_draws = 4000L,
    protein_min = 0.8,
    grid_step = 0.01,
    seed = 34988L
) {
  if (
    !is.data.frame(accepted_fit_manifest) ||
      nrow(accepted_fit_manifest) != 104L ||
      any(accepted_fit_manifest$family != "strategy_within") ||
      !is.data.frame(stacking_weights) ||
      nrow(stacking_weights) != 104L ||
      any(stacking_weights$family != "strategy_within")
  ) {
    stop("The full within-strategy prediction inputs are invalid.",
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
    record$family <- "strategy_within"
    record$prediction_type <- "strategy_within"
    record
  })
}

write_protein_strategy_within_full_prediction_draw_file <- function(
    record,
    checkpoint_path
) {
  if (
    length(record$fit_path) != 1L || is.na(record$fit_path) ||
      !file.exists(record$fit_path) || record$draw_count < 1L ||
      !identical(unname(as.character(record$family)), "strategy_within") ||
      !identical(
        unname(as.character(record$prediction_type)), "strategy_within"
      )
  ) {
    stop("A full within-strategy prediction record is invalid.",
         call. = FALSE)
  }
  candidate <- read_sensitivity_fit_file(record$fit_path)
  if (
    !identical(unname(as.character(candidate$family)), "strategy_within") ||
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
    stop("A within-strategy prediction fit does not match its record.",
         call. = FALSE)
  }
  checkpoint <- read_sensitivity_completed_checkpoint(
    checkpoint_path, family = "strategy_within"
  )
  prediction_data <- checkpoint$analysis_data_by_outcome[[record$outcome]]
  if (is.null(prediction_data) || !is.data.frame(prediction_data)) {
    stop("The within-strategy prediction data are unavailable.",
         call. = FALSE)
  }
  grid <- make_protein_strategy_within_prediction_grid(
    fitted_candidate = candidate,
    data = prediction_data,
    protein_min = record$protein_min,
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
    stop("Within-strategy prediction draws do not match their allocation.",
         call. = FALSE)
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
  path <- sensitivity_reduction_file_path(
    family = "strategy_within",
    outcome = record$outcome,
    model_index = record$model_index,
    stage = "prediction_draws"
  )
  write_atomic_qs(reduction, path)
}

write_protein_strategy_within_full_reporting_file <- function(
    model_average_path,
    data_by_outcome,
    reference_intake = 0.8,
    sesoi = 0.10
) {
  model_average <- read_sensitivity_full_model_average_file(
    model_average_path, family = "strategy_within"
  )
  outcome <- model_average$outcome
  data <- data_by_outcome[[outcome]]
  if (
    is.null(data) || !is.data.frame(data) ||
      !all(c(
        "protein_strategy", "protein_intake_raw_treatment"
      ) %in% names(data))
  ) {
    stop("The full within-strategy reporting data are unavailable.",
         call. = FALSE)
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
        model_average, strategy
      )
      strategy_data <- dplyr::filter(
        data, protein_strategy == strategy
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
        dplyr::mutate(
          table,
          protein_strategy = strategy,
          .before = 1L
        )
      }
      list(
        increments = annotate(supported$increments_full_grid),
        support = annotate(supported$support_by_group),
        common_support = annotate(supported$common_support),
        pairwise = annotate(heatmap$contrasts),
        pairwise_support = annotate(heatmap$support),
        local_equivalence = annotate(
          summarise_primary_local_equivalence(
            heatmap$contrasts, outcome = outcome
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
  path <- sensitivity_reduction_file_path(
    family = "strategy_within",
    outcome = outcome,
    model_index = 0L,
    stage = "reporting"
  )
  write_atomic_qs(report, path)
}

read_protein_strategy_within_full_reporting_file <- function(
    path,
    family = "strategy_within"
) {
  report <- read_sensitivity_full_reporting_file(
    path, family = family
  )
  if (
    !"addition_minus_substitution" %in% names(report) ||
      !is.data.frame(report$addition_minus_substitution) ||
      !all(c("Substitution", "Addition") %in%
           unique(report$dose_response$protein_strategy))
  ) {
    stop("A full within-strategy reporting file is invalid.",
         call. = FALSE)
  }
  report
}

audit_protein_strategy_within_full_reporting_files <- function(
    paths,
    family = "strategy_within"
) {
  paths <- unname(unlist(paths, use.names = FALSE))
  if (length(paths) != 4L || anyNA(paths) || any(!file.exists(paths))) {
    stop("The within-strategy reporting-file set is incomplete.",
         call. = FALSE)
  }
  summary <- purrr::map_dfr(paths, function(path) {
    report <- read_protein_strategy_within_full_reporting_file(path, family)
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
    ) && all(
      is.finite(
        report$addition_minus_substitution$addition_minus_substitution
      ) &
        report$addition_minus_substitution$conf.low <=
          report$addition_minus_substitution$addition_minus_substitution &
        report$addition_minus_substitution$addition_minus_substitution <=
          report$addition_minus_substitution$conf.high
    )
    reference_valid <-
      nrow(reference) == 4L &&
      all(abs(reference$increment) < 1e-12) &&
      all(abs(reference$increment_conf_low) < 1e-12) &&
      all(abs(reference$increment_conf_high) < 1e-12)
    probabilities <- as.matrix(
      report$pairwise_contrasts[, probability_columns, drop = FALSE]
    )
    strategy_probabilities <-
      report$addition_minus_substitution$probability_positive
    probabilities_valid <-
      all(is.finite(probabilities)) &&
      all(probabilities >= 0 & probabilities <= 1) &&
      all(is.finite(strategy_probabilities)) &&
      all(strategy_probabilities >= 0 & strategy_probabilities <= 1) &&
      all(abs(
        report$pairwise_contrasts$probability_meaningful_gain +
          report$pairwise_contrasts$probability_negligible +
          report$pairwise_contrasts$probability_meaningful_loss - 1
      ) < 1e-12)
    diagonal_valid <-
      nrow(diagonal) >= 4L &&
      all(abs(diagonal$estimate) < 1e-12) &&
      all(abs(diagonal$probability_negligible - 1) < 1e-12)
    tibble::tibble(
      family = report$family,
      outcome = report$outcome,
      dose_response_rows = nrow(report$dose_response),
      increment_rows = nrow(report$dose_response_increments),
      strategy_contrast_rows = nrow(report$addition_minus_substitution),
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
      any(summary$dose_response_rows < 4L) ||
      any(summary$strategy_contrast_rows < 2L) ||
      any(summary$pairwise_contrast_rows < 4L) ||
      any(abs(summary$weight_sum - 1) > 1e-8) ||
      !all(summary$intervals_valid) ||
      !all(summary$reference_increment_valid) ||
      !all(summary$probabilities_valid) ||
      !all(summary$diagonal_contrasts_valid)
  ) {
    stop("The full within-strategy reporting audit failed.",
         call. = FALSE)
  }
  list(
    family = family,
    summary = summary,
    outcomes = nrow(summary),
    posterior_mixture_draws_per_outcome = 4000L,
    all_reporting_checks_pass = TRUE
  )
}

write_protein_strategy_within_full_verification_bundle <- function(
    reporting_audit,
    checkpoint,
    prediction_manifest
) {
  if (
    !identical(checkpoint$family, "strategy_within") ||
      !isTRUE(reporting_audit$all_reporting_checks_pass) ||
      !isTRUE(checkpoint$audit$all_sampling_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$all_exact_refit_diagnostics_pass) ||
      !isTRUE(checkpoint$audit$stacking_ready)
  ) {
    stop("The within-strategy family is not ready for verification.",
         call. = FALSE)
  }
  reports <- purrr::map(
    reporting_audit$summary$reporting_path,
    read_protein_strategy_within_full_reporting_file
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
    family = "strategy_within",
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
      addition_minus_substitution =
        bind_report_table("addition_minus_substitution"),
      pairwise_contrasts = bind_report_table("pairwise_contrasts"),
      pairwise_support = bind_report_table("pairwise_support"),
      local_equivalence = bind_report_table("local_equivalence")
    ),
    descriptions = c(
      analysis_data = "Strict matched-arm contrasts joined to the audited randomized protein-manipulation strategy.",
      model_specifications = "Final candidate-model registry inherited from the primary analysis.",
      accepted_fit_manifest = "Accepted one-file-per-model registry and fit checksums.",
      sampling_diagnostic_summary = "Final model-level convergence, geometry, and predictive diagnostics.",
      sampling_parameter_failures = "Any parameter-level diagnostic failures after remediation.",
      initial_integrated_loo_model_summary = "Integrated observation-level PSIS summaries before exact correction.",
      initial_integrated_loo_outcome_summary = "Outcome-level integrated PSIS audit before exact correction.",
      initial_high_pareto_k_cases = "Cases selected prospectively for exact case-deletion correction.",
      exact_refit_diagnostic_summary = "Sampling diagnostics for every exact case-deletion refit.",
      exact_loo_correction_summary = "Exact versus PSIS pointwise correction results.",
      final_loo_model_summary = "Final corrected integrated LOO summaries used for stacking.",
      stacking_weights = "Final outcome-specific stacking weights.",
      prediction_allocation_manifest = "Deterministic allocation of 4,000 mixture curves per outcome.",
      reporting_file_manifest = "Checksums and dimensions of reduced outcome reporting files.",
      dose_response_predictions = "Model-averaged within-study protein-arm minus comparator-arm curves under addition and substitution strategies.",
      dose_response_increments = "Strategy-specific paired within-study increments from 0.8 g/kg/day.",
      exposure_support_by_group = "Observed treatment-intake support by strategy and resistance-exercise condition.",
      common_exposure_support = "Treatment-intake support shared by both strategies within resistance-exercise condition.",
      addition_minus_substitution = "Direct posterior addition-minus-substitution contrast in the within-study estimand at each dose and exercise condition.",
      pairwise_contrasts = "Strategy-specific paired dose contrasts and SESOI probabilities.",
      pairwise_support = "Observed treatment-intake support used to flag strategy-specific pairwise contrasts.",
      local_equivalence = "Strategy-specific local practical-equivalence summaries."
    )
  )
}
