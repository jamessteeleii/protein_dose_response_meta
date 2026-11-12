make_prisma_study_characteristics <- function(workbook, primary) {
  format_range <- function(x, digits = 1L) {
    values <- x[is.finite(x)]
    if (!length(values)) return("&mdash;")
    bounds <- unique(formatC(range(values), digits = digits, format = "f"))
    paste(bounds, collapse = "&ndash;")
  }

  list_values <- function(x) {
    values <- unique(trimws(x[!is.na(x)]))
    values <- values[nzchar(values)]
    if (length(values)) paste(values, collapse = ", ") else "&mdash;"
  }

  format_baseline_n <- function(x) {
    known <- is.finite(x)
    if (!any(known)) return("&mdash;")
    total <- format(sum(x[known]), big.mark = ",", scientific = FALSE, trim = TRUE)
    if (all(known)) return(total)
    sprintf("%s (%d of %d arms)", total, sum(known), length(x))
  }

  make_report_link <- function(author, year, doi, pmid, title) {
    year <- stringr::str_extract(year, "[12][0-9]{3}")
    label <- paste0(author, ifelse(is.na(year), "", paste0(" (", year, ")")))
    dois <- unique(unlist(stringr::str_extract_all(
      doi, "10\\.[0-9]{4,9}/[^ ;,|\\s]+"
    )))
    dois <- dois[!is.na(dois)]
    if (length(dois)) {
      urls <- paste0("https://doi.org/", vapply(
        dois, utils::URLencode, character(1), reserved = FALSE
      ))
      links <- sprintf("[%s](<%s>)", c(label, rep("Related report", length(dois) - 1L)), urls)
      return(paste(links, collapse = " &middot; "))
    }
    if (!is.na(pmid) && grepl("^[0-9]+$", pmid)) {
      return(sprintf("[%s](https://pubmed.ncbi.nlm.nih.gov/%s/)", label, pmid))
    }
    if (is.na(title) || !nzchar(title)) label else paste0(label, ". ", title)
  }

  arms <- workbook$arm_characteristics
  stopifnot(!anyDuplicated(arms[c("study_no", "arm_no")]))
  arm_summary <- arms |>
    dplyr::transmute(
      study_no,
      age = mean_age_years,
      male = 100 * proportion_male_0_1_i_e_all_male_group_1_all_female_group_0,
      n = sample_size_at_baseline,
      rt = resistance_exercise_yes_no,
      condition = intervention_i_e_protein_supplement_or_comparator_i_e_placebo_or_control,
      baseline_intake = baseline_dietary_protein_intake_g_kg_day,
      intervention_intake = intervention_dietary_protein_intake_g_kg_day,
      dose = absolute_protein_dose_g_day_note_if_placebo_condition_then_check_whether_it_includes_any_protein
    ) |>
    dplyr::group_by(study_no) |>
    dplyr::summarise(
      baseline_n = format_baseline_n(n),
      mean_age = format_range(age),
      male_pct = if (dplyr::first(study_no) == 164L) "&mdash;" else format_range(male, 0L),
      arm_types = list_values(condition),
      dose = format_range(dose[which(condition == "protein")]),
      baseline_protein = format_range(baseline_intake, 2L),
      intervention_protein = format_range(intervention_intake, 2L),
      rt = list_values(rt),
      .groups = "drop"
    )

  followup <- workbook$arm_based_data_by_timepoint |>
    dplyr::filter(timepoint_week_note_baseline_0 > 0) |>
    dplyr::group_by(study_no) |>
    dplyr::summarise(
      followup_weeks = format_range(timepoint_week_note_baseline_0),
      .groups = "drop"
    )

  outcome_labels <- c(
    "whole body lean mass" = "Whole-body lean/fat-free mass",
    "direct muscle size/mass" = "Measures of muscle size or mass",
    strength = "Strength",
    "function" = "Physical function"
  )
  stopifnot(all(primary$outcome_group %in% names(outcome_labels)))
  contribution <- primary |>
    dplyr::group_by(study_no) |>
    dplyr::summarise(
      primary_outcomes = list_values(
        unname(outcome_labels[names(outcome_labels) %in% outcome_group])
      ),
      .groups = "drop"
    )

  studies <- workbook$study_characteristics
  stopifnot(
    !anyNA(studies$study_no), !anyDuplicated(studies$study_no),
    setequal(studies$study_no, arm_summary$study_no),
    all(primary$study_no %in% studies$study_no)
  )
  studies |>
    dplyr::mutate(report = purrr::pmap_chr(
      list(
        authors_e_g_smith_or_smith_and_jones_or_smith_et_al,
        year_of_publication, doi, pmid, title_of_article
      ),
      make_report_link
    )) |>
    dplyr::left_join(arm_summary, by = "study_no") |>
    dplyr::left_join(followup, by = "study_no") |>
    dplyr::left_join(contribution, by = "study_no") |>
    dplyr::arrange(study_no) |>
    dplyr::transmute(
      ID = study_no,
      Study = report,
      Arms = how_many_arms_included,
      `Baseline n` = baseline_n,
      `Age (years)` = mean_age,
      `Male (%)` = male_pct,
      `Extracted arm types` = arm_types,
      `Prescribed protein (g/day)` = dose,
      `Baseline total (g/kg/day)` = baseline_protein,
      `Intervention total (g/kg/day)` = intervention_protein,
      RT = rt,
      `Follow-up (weeks)` = dplyr::coalesce(followup_weeks, "&mdash;"),
      `Primary outcomes contributed` = dplyr::coalesce(primary_outcomes, "None")
    )
}
