# Reviewed leaf-file pruning for completed non-primary absolute families.
# Source this file, then call iad_storage_audit(). Deletion requires an explicit
# call to iad_prune_completed_family(..., execute = TRUE). No targets are run.
iad_storage_families <- function() {
  data.frame(
    key = c(paste0("total_", c("mundlak", "within", "achieved", "assignment", "eiv", "rescor", "strategy", "strategy_within")),
      paste0("supplement_", c("mundlak", "within", "rescor", "strategy", "strategy_within"))),
    exposure = c(rep("total_g_day", 8), rep("supplement_g_day", 5)),
    folder = c("mundlak", "within", "achieved_common_sample", "assignment_common_sample", "eiv", "rescor", "strategy", "strategy_within",
      "mundlak", "within", "rescor", "strategy", "strategy_within"), stringsAsFactors = FALSE)
}
iad_storage_protected <- function() {
  roots <- c("primary_update_2026", "total_primary_update_2026", "supplement_primary_update_2026")
  paths <- character()
  for (root in roots) {
    directory <- file.path("results", "verification", root)
    fit_file <- if (root == "primary_update_2026") "accepted_fit_manifest.csv.gz" else "final_fit_manifest.csv.gz"
    loo_file <- if (root == "primary_update_2026") "loo_final_manifest.csv.gz" else "final_loo_manifest.csv.gz"
    fits <- read.csv(gzfile(file.path(directory, fit_file)), stringsAsFactors = FALSE)
    loo <- read.csv(gzfile(file.path(directory, loo_file)), stringsAsFactors = FALSE)
    paths <- c(paths, fits$fit_path, loo$final_loo_path)
  }
  unique(tolower(gsub("\\", "/", paths, fixed = TRUE)))
}
iad_storage_plan <- function(key, verify = TRUE) {
  families <- iad_storage_families()
  row <- families[families$key == key, ]
  if (nrow(row) != 1L) stop("Only the 13 named completed non-primary absolute families are eligible.")
  root <- normalizePath(file.path("model_archive", "update_2026", "absolute_exposure", row$exposure, row$folder),
    winslash = "/", mustWork = TRUE)
  directory <- file.path("results", "update_2026", "absolute_exposure", row$exposure, row$folder, "pruning")
  selected <- read.csv(file.path(directory, "deletion_file_manifest.csv"), stringsAsFactors = FALSE)
  compact <- read.csv(file.path(directory, "compact_archive_manifest.csv"), stringsAsFactors = FALSE)
  bundle <- file.path("results", "verification", paste0(key, "_update_2026"))
  manifest <- read.csv(file.path(bundle, "manifest.csv"), stringsAsFactors = FALSE)
  check_path <- file.path(bundle, manifest$file[manifest$artifact == "verification_checks"])
  checks <- read.csv(gzfile(check_path), stringsAsFactors = FALSE)
  if (!nrow(checks) || !all(checks$passed)) stop("Final verification has not passed: ", key)
  existing <- selected[file.exists(selected$path), , drop = FALSE]
  resolved <- normalizePath(existing$path, winslash = "/", mustWork = TRUE)
  if (anyDuplicated(tolower(resolved)) || any(!startsWith(tolower(resolved), paste0(tolower(root), "/"))) ||
      any(file.info(resolved)$isdir) || any(tolower(resolved) %in% iad_storage_protected()) ||
      any(tolower(resolved) %in% tolower(gsub("\\", "/", compact$path, fixed = TRUE)))) {
    stop("Deletion plan overlaps a protected artifact or leaves its named family root.")
  }
  if (any(file.info(resolved)$size != existing$bytes)) stop("A proposed deletion changed size; rebuild/review the plan.")
  verify_hash <- function(paths, hashes) {
    if (any(!file.exists(paths)) || anyNA(hashes)) stop("A retained artifact is missing.")
    observed <- vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1))
    if (any(observed != hashes)) stop("A retained artifact failed its SHA-256 check.")
  }
  if (verify) {
    if (any(basename(manifest$file) != manifest$file)) stop("Unsafe verification manifest filename.")
    verify_hash(file.path(bundle, manifest$file), manifest$sha256)
    verify_hash(compact$path, compact$sha256)
  }
  list(key = key, directory = directory, root = root, selected = existing,
    retained_manifest = manifest, compact = compact,
    summary = data.frame(family = key, present_files = nrow(existing),
      reclaimable_gib = sum(file.info(resolved)$size) / 1024^3,
      retained_artifacts_verified = verify, stringsAsFactors = FALSE))
}
iad_storage_audit <- function() {
  rows <- lapply(iad_storage_families()$key, function(key) {
    plan <- iad_storage_plan(key, verify = TRUE)
    plan$summary
  })
  table <- do.call(rbind, rows)
  table <- table[order(table$reclaimable_gib, decreasing = TRUE), ]
  dir.create("results/update_2026/age_duration/storage", recursive = TRUE, showWarnings = FALSE)
  write.csv(table, "results/update_2026/age_duration/storage/pruning_audit.csv", row.names = FALSE)
  print(table, row.names = FALSE)
  cat("Total reclaimable GiB:", sum(table$reclaimable_gib), "\nNo files deleted.\n")
  invisible(table)
}
iad_prune_completed_family <- function(key, execute = FALSE) {
  # Recheck current paths, verification and every retained compact hash in the
  # same process immediately before deletion. Never delete a directory.
  plan <- iad_storage_plan(key, verify = TRUE)
  print(plan$summary, row.names = FALSE)
  if (!isTRUE(execute)) return(invisible(plan))
  if (!nrow(plan$selected)) return(invisible(plan))
  receipt <- plan$selected
  receipt$sha256_before_deletion <- vapply(receipt$path, function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1))
  receipt$deleted_utc <- NA_character_
  receipt_path <- file.path(plan$directory, "executed_for_age_duration.csv")
  if (file.exists(receipt_path) || file.exists(paste0(receipt_path, ".pending"))) {
    stop("An age/duration pruning receipt already exists; review it before retrying.")
  }
  write.csv(receipt, paste0(receipt_path, ".pending"), row.names = FALSE, na = "")
  for (i in seq_len(nrow(receipt))) {
    path <- normalizePath(receipt$path[[i]], winslash = "/", mustWork = TRUE)
    if (!startsWith(tolower(path), paste0(tolower(plan$root), "/")) ||
        tolower(path) %in% iad_storage_protected()) stop("Protected path encountered before deletion.")
    if (unlink(path, recursive = FALSE) != 0L || file.exists(path)) stop("Could not delete planned leaf file: ", path)
    receipt$deleted_utc[[i]] <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    write.csv(receipt, paste0(receipt_path, ".pending"), row.names = FALSE, na = "")
  }
  if (!file.rename(paste0(receipt_path, ".pending"), receipt_path)) stop("Could not finalise pruning receipt.")
  cat("Deleted only planned files for", key, "; retained compact outputs and primary archives.\n")
  invisible(receipt)
}
