#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(targets)
})

archive_dir <- file.path("model_archive", "primary")
manifest_dir <- file.path("results", "freeze")
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

target_names <- c(
  "model_fit_final_whole_body",
  "model_fit_final_hypertrophy",
  "model_fit_final_strength",
  "model_fit_final_function"
)
meta <- targets::tar_meta(fields = c(name, path, bytes, time))
meta <- meta[match(target_names, meta$name), , drop = FALSE]
if (any(is.na(meta$name))) {
  stop("One or more final primary-fit targets are absent.", call. = FALSE)
}

rows <- vector("list", nrow(meta))
for (index in seq_len(nrow(meta))) {
  source_path <- file.path("_targets", "objects", meta$name[[index]])
  if (!file.exists(source_path)) {
    stop("Missing cached object for ", meta$name[[index]], call. = FALSE)
  }
  outcome <- sub("^model_fit_final_", "", meta$name[[index]])
  destination <- file.path(archive_dir, paste0(outcome, ".qs"))
  source_hash <- digest::digest(file = source_path, algo = "sha256")
  if (!file.exists(destination) ||
      !identical(digest::digest(file = destination, algo = "sha256"), source_hash)) {
    copied <- file.copy(source_path, destination, overwrite = TRUE, copy.date = TRUE)
    if (!copied) stop("Could not archive ", meta$name[[index]], call. = FALSE)
  }
  archive_hash <- digest::digest(file = destination, algo = "sha256")
  if (!identical(source_hash, archive_hash)) {
    stop("Checksum mismatch after archiving ", meta$name[[index]], call. = FALSE)
  }
  rows[[index]] <- data.frame(
    target = meta$name[[index]],
    outcome = outcome,
    archive_file = gsub("\\\\", "/", destination),
    bytes = file.info(destination)$size,
    sha256 = archive_hash,
    target_timestamp = format(meta$time[[index]], tz = "UTC", usetz = TRUE),
    archived_utc = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
    stringsAsFactors = FALSE
  )
  cat("Archived and verified ", meta$name[[index]], "\n", sep = "")
}

manifest <- do.call(rbind, rows)
utils::write.csv(
  manifest,
  file.path(manifest_dir, "primary-fit-archive-manifest.csv"),
  row.names = FALSE,
  na = ""
)
cat("Primary fit archive manifest written.\n")
