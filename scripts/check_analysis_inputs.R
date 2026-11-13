# Run from the project root after restoring renv. This check does not sample.
suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
  library(stringr)
  library(tidyr)
})
source("R/functions.R")
source("R/incremental_update.R")
source("R/sensitivity_strategy.R")
source("R/update_2026_absolute_exposures.R")
source("R/grade.R")

canonical <- read_protein_workbook("data/protein_meta_data_all.xlsx")
update <- read_protein_workbook("data/update_2026/protein_meta_data_update_2026.xlsx")
decisions <- read_update_review_decisions("data/update_2026/review_decisions.csv")
reference <- read_frozen_primary_analysis("data/update_2026/frozen_pre_update_primary.csv")
workbook_checks <- validate_update_workbook(update, reference, decisions)

combined <- combine_protein_workbooks(canonical, update)
reviewed <- apply_combined_review_decisions(prepare_protein_data(combined), decisions)
outcomes <- make_combined_primary_data_by_outcome(reviewed)
preparation <- validate_combined_primary_preparation(
  combined, reviewed, outcomes, reference, decisions
)
stopifnot(
  preparation$frozen_effect_size_change$rows == nrow(reference),
  identical(preparation$outcome_counts$contrasts, c(283L, 425L, 794L, 362L))
)

keys <- c("study_no", "arm_no", "outcome_no", "time1", "time2")
fields <- c(keys, "yi", "vi", "protein_intake", "resistance_exercise_code", "n")
published <- readr::read_csv(
  "manuscript/generated/reporting/primary-study-effects.csv.gz",
  col_select = all_of(fields), show_col_types = FALSE
)
prepared <- bind_rows(outcomes) |>
  select(all_of(fields)) |>
  arrange(across(all_of(keys)))
published <- arrange(published, across(all_of(keys)))
stopifnot(
  !anyDuplicated(prepared[keys]),
  isTRUE(all.equal(prepared, published, tolerance = 1e-9, check.attributes = FALSE))
)

strategy <- read_update_2026_protein_strategy_audit(
  "data/protein_manipulation_strategy_audit.csv",
  "data/update_2026/protein_manipulation_strategy_audit_additions.csv"
)
strategy_checks <- audit_update_2026_protein_strategy(strategy, outcomes)
stopifnot(strategy_checks$audit_ready)
supplement_strategy <- read_update_2026_supplement_strategy_audit(
  "data/update_2026/protein_manipulation_strategy_audit_supplement_only.csv",
  strategy
)
stopifnot(nrow(supplement_strategy) == nrow(strategy) + 2L)

rob2 <- make_grade_rob2_summary(outcomes, "data/GRADE/rob2_overall_judgments.csv")
stopifnot(all(rob2$summary$assessed_studies == rob2$summary$contributing_studies))

print(workbook_checks)
print(preparation$checks)
print(preparation$outcome_counts)
print(strategy_checks$structural_checks)
cat("Verified the bundled extraction, strategy and GRADE inputs.\n")
