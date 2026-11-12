args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this script with Rscript.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
data_dir <- file.path(project_root, "supplement", "data")
include_dir <- file.path(project_root, "supplement", "includes")
dir.create(include_dir, recursive = TRUE, showWarnings = FALSE)

read_utf8 <- function(path) {
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    encoding = "UTF-8"
  )
}

escape_md <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("[\r\n]+", " ", x)
  gsub("|", "\\|", x, fixed = TRUE)
}

format_count <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

markdown_table <- function(headers, rows, align = rep("l", length(headers))) {
  separators <- ifelse(
    align == "r", "---:",
    ifelse(align == "c", ":---:", ":---")
  )
  rows <- as.matrix(rows)
  body <- apply(
    rows,
    1L,
    function(row) paste0("| ", paste(escape_md(row), collapse = " | "), " |")
  )
  c(
    paste0("| ", paste(escape_md(headers), collapse = " | "), " |"),
    paste0("| ", paste(separators, collapse = " | "), " |"),
    body
  )
}

write_include <- function(name, caption, headers, rows, align) {
  writeLines(
    c(caption, "", markdown_table(headers, rows, align), ""),
    file.path(include_dir, name),
    useBytes = TRUE
  )
}

search <- read_utf8(file.path(data_dir, "prisma-search-counts.csv"))
flow <- read_utf8(file.path(data_dir, "prisma-flow-counts.csv"))
excluded <- read_utf8(file.path(data_dir, "prisma-full-text-exclusions.csv"))
checklist <- read_utf8(file.path(data_dir, "prisma-checklist-register.csv"))
update <- read_utf8(file.path(data_dir, "update-2026-study-status.csv"))

source(file.path(project_root, "R", "functions.R"))
source(file.path(project_root, "R", "incremental_update.R"))
source(file.path(dirname(script_path), "prisma-study-characteristics.R"))
workbook <- combine_protein_workbooks(
  read_protein_workbook(file.path(project_root, "data", "protein_meta_data_all.xlsx")),
  read_protein_workbook(file.path(
    project_root, "data", "update_2026", "protein_meta_data_update_2026.xlsx"
  ))
)
primary <- readr::read_csv(file.path(
  project_root, "manuscript", "generated", "reporting", "primary-study-effects.csv.gz"
), col_select = c(study_no, outcome_group), show_col_types = FALSE)
characteristics <- make_prisma_study_characteristics(workbook, primary)
stopifnot(nrow(characteristics) == flow$count[flow$stage == "Total studies included in current review"])
write_include(
  "prisma-study-characteristics.md",
  "**Characteristics of the included studies.**",
  names(characteristics), characteristics,
  c("r", "l", "r", "r", "r", "r", "l", "r", "r", "r", "l", "r", "l")
)

for (wave in c("July 2025", "September 2026")) {
  x <- search[search$search_wave == wave, , drop = FALSE]
  include_coverage <- identical(wave, "September 2026")
  rows <- data.frame(
    Database = x$database,
    Protein = format_count(x$protein_search),
    Protein_RE = format_count(x$protein_plus_re_search),
    Total = format_count(x$total),
    stringsAsFactors = FALSE
  )
  if (include_coverage) rows$Coverage <- x$coverage_note
  total <- data.frame(
    Database = "**Total**",
    Protein = paste0("**", format_count(sum(x$protein_search)), "**"),
    Protein_RE = paste0("**", format_count(sum(x$protein_plus_re_search)), "**"),
    Total = paste0("**", format_count(sum(x$total)), "**"),
    stringsAsFactors = FALSE
  )
  if (include_coverage) total$Coverage <- ""
  rows <- rbind(rows, total)

  write_include(
    if (wave == "July 2025") "prisma-search-july-2025.md" else "prisma-search-september-2026.md",
    if (wave == "July 2025") {
      "**Database yields for the 13 July 2025 search.**"
    } else {
      "**Database yields for the 14 September 2026 update.**"
    },
    c("Database", "Protein", "Protein + resistance exercise", "Total",
      if (include_coverage) "Coverage"),
    rows,
    c("l", "r", "r", "r", if (include_coverage) "l")
  )
}

write_include(
  "prisma-flow-ledger.md",
  "**Selection-flow ledger and reconciliation status.**",
  c("Search", "Stage", "Count", "Calculation or note", "Status"),
  data.frame(
    flow$search_wave,
    flow$stage,
    format_count(flow$count),
    flow$calculation_or_note,
    flow$reconciliation_status
  ),
  c("l", "l", "r", "l", "l")
)

for (wave in c("July 2025", "September 2026")) {
  x <- excluded[excluded$search_wave == wave, , drop = FALSE]
  write_include(
    if (wave == "July 2025") "prisma-exclusions-july-2025.md" else "prisma-exclusions-september-2026.md",
    if (wave == "July 2025") {
      "**Reports excluded after full-text assessment in the July 2025 search.**"
    } else {
      "**Reports excluded after eligibility/full-text assessment in the September 2026 update.**"
    },
    c("No.", "Complete reference", "Recorded reason"),
    x[c("number", "citation", "recorded_reason")],
    c("r", "l", "l")
  )
}

write_include(
  "prisma-update-study-status.md",
  "**Study-level status of the September 2026 extraction tranche.**",
  c(
    "Study", "Citation", "DOI", "Arms", "Timepoint rows", "Contrast rows",
    "Analysis status", "Audit note"
  ),
  data.frame(
    update$study_no,
    paste0(update$authors, " (", update$year, ")"),
    update$doi,
    update$arms,
    update$timepoint_rows,
    update$contrast_rows,
    update$analysis_status,
    update$audit_note
  ),
  c("r", "l", "l", "r", "r", "r", "l", "l")
)

for (kind in c("Main", "Abstract")) {
  x <- checklist[checklist$checklist == kind, , drop = FALSE]
  write_include(
    if (kind == "Main") "prisma-checklist-main.md" else "prisma-checklist-abstract.md",
    if (kind == "Main") {
      "**PRISMA 2020 main checklist.**"
    } else {
      "**PRISMA 2020 for Abstracts checklist.**"
    },
    c(
      "Item", "Topic", "Checklist requirement", "Reporting location", "Manuscript page"
    ),
    x[c("item", "topic", "requirement", "reporting", "manuscript_page")],
    c("l", "l", "l", "l", "l")
  )
}

message("PRISMA Markdown includes rebuilt from the editable CSV ledgers.")
