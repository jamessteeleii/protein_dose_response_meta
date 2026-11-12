#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(digest))

root <- file.path("results", "verification")
manifest_paths <- list.files(
  root,
  pattern = "^manifest\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(manifest_paths)) {
  stop("No verification manifests were found.", call. = FALSE)
}

audit_one <- function(manifest_path) {
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  paths <- file.path(dirname(manifest_path), manifest$file)
  exists <- file.exists(paths)
  observed_hash <- rep(NA_character_, length(paths))
  observed_hash[exists] <- vapply(
    paths[exists],
    digest::digest,
    character(1),
    file = TRUE,
    algo = "sha256"
  )
  data.frame(
    family = basename(dirname(manifest_path)),
    manifest = gsub("\\\\", "/", manifest_path),
    artifacts = nrow(manifest),
    missing_files = sum(!exists),
    checksum_failures = sum(exists & observed_hash != manifest$sha256),
    bytes = sum(manifest$bytes, na.rm = TRUE),
    verified_utc = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, lapply(manifest_paths, audit_one))
audit <- audit[order(audit$family), , drop = FALSE]
dir.create(file.path("results", "freeze"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  audit,
  file.path("results", "freeze", "verification-bundle-audit.csv"),
  row.names = FALSE,
  na = ""
)
print(audit, row.names = FALSE)
if (any(audit$missing_files > 0L | audit$checksum_failures > 0L)) {
  stop("At least one verification bundle failed validation.", call. = FALSE)
}

fit_manifest_paths <- list.files(
  root,
  pattern = "^accepted_fit_manifest\\.csv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)
fit_retention <- do.call(rbind, lapply(fit_manifest_paths, function(path) {
  manifest <- utils::read.csv(gzfile(path), stringsAsFactors = FALSE)
  data.frame(
    family = basename(dirname(path)),
    accepted_fits = nrow(manifest),
    fit_archives_present = sum(file.exists(manifest$fit_path)),
    fit_archives_absent = sum(!file.exists(manifest$fit_path)),
    stringsAsFactors = FALSE
  )
}))
fit_retention <- fit_retention[order(fit_retention$family), , drop = FALSE]
utils::write.csv(
  fit_retention,
  file.path("results", "freeze", "sensitivity-fit-retention-audit.csv"),
  row.names = FALSE,
  na = ""
)
print(fit_retention, row.names = FALSE)
