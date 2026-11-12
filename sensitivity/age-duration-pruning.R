# Prepare and verify a lean plan for the completed interaction fits.
# Run from the project root. This script writes manifests only.
iad_plan_interaction_pruning <- function() {
  source("R/age_duration_interactions.R", local = TRUE)
  project <- normalizePath(".", winslash = "/", mustWork = TRUE)
  archive <- file.path(project, "model_archive", "update_2026", "age_duration")
  interface <- file.path(project, "manuscript", "generated", "age-duration-interactions")
  relative_path <- function(path) substring(path, nchar(project) + 2L)
  key <- function(rows) paste(rows$family, rows$outcome, rows$model_index, sep = "/")
  case_key <- function(rows) paste(key(rows), rows$observation, sep = "/")

  retain <- function(paths, role, expected = NULL) {
    paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
    info <- file.info(paths)
    if (any(!startsWith(tolower(paths), paste0(tolower(project), "/"))) ||
        anyNA(info$size) || any(info$isdir) || any(info$size <= 0)) {
      stop("Invalid retained file for ", role, ".")
    }
    hashes <- vapply(paths, iad_hash, character(1), USE.NAMES = FALSE)
    if (!is.null(expected) && !identical(hashes, unname(as.character(expected)))) {
      stop("Retained checksum mismatch for ", role, ": ", paths[which(hashes != expected)[1L]])
    }
    data.frame(role = role, path = relative_path(paths), bytes = info$size, sha256 = hashes)
  }

  results <- file.path(interface, "_results.qmd")
  if (sum(grepl("^### ", readLines(results))) != 12L) stop("Finish the supplementary assembly first.")
  retained <- list(retain(c(results, "supplement/age-duration-interactions.qmd"), "supplement source"))
  deletion <- list()

  for (family in iad_families()) {
    paths <- function(...) iad_path(family, ...)
    initial_path <- paths("diagnostics", "initial", "summary.csv")
    initial_manifest_path <- paths("checkpoints", "initial_fits.csv")
    accepted_path <- paths("diagnostics", "accepted.csv")
    integrated_path <- paths("loo", "manifest.csv")
    cases_path <- paths("loo", "exact_cases.csv")
    reviewed_path <- paths("loo", "exact_review_diagnostics.csv")
    report_path <- paths("checkpoints", "reporting.csv")
    initial <- iad_read_csv(initial_path)
    initial_manifest <- iad_read_csv(initial_manifest_path)
    accepted <- iad_read_csv(accepted_path)
    integrated <- iad_read_csv(integrated_path)
    cases <- iad_read_csv(cases_path)
    reviewed <- iad_read_csv(reviewed_path)
    report <- iad_read_csv(report_path)
    models <- if (family == "relative") 27:62 else 27:50
    expected_keys <- key(expand.grid(family = family, outcome = iad_outcomes(), model_index = models))
    for (table in list(initial, initial_manifest, accepted, integrated)) {
      if (anyDuplicated(key(table)) || !setequal(key(table), expected_keys)) {
        stop("Incomplete interaction model registry for ", family, ".")
      }
    }
    integrated <- integrated[match(key(accepted), key(integrated)), ]
    if (!all(accepted$diagnostic_pass) || !identical(accepted$sha256, integrated$fit_sha256) ||
        !identical(accepted$fit_path, integrated$fit_path)) stop("Accepted fit and LOO identity differ for ", family, ".")
    initial_manifest <- initial_manifest[match(key(initial), key(initial_manifest)), ]
    if (!identical(initial$fit_path, initial_manifest$fit_path) ||
        !all(initial$diagnostic_pass == (initial$fit_path %in% accepted$fit_path))) {
      stop("The initial and accepted fit records disagree for ", family, ".")
    }

    columns <- c("family", "outcome", "model_index", "stage", "fit_path", "sha256")
    fits <- rbind(accepted[columns], initial[!initial$fit_path %in% accepted$fit_path, columns])
    if (anyDuplicated(fits$fit_path) || !all(grepl("^(initial|remedy_[1-9][0-9]*)$", fits$stage)) ||
        !all(grepl("^[a-f0-9]{64}$", fits$sha256))) stop("Invalid interaction fit registry for ", family, ".")
    resolved <- normalizePath(fits$fit_path, winslash = "/", mustWork = TRUE)
    expected_paths <- file.path(archive, family, fits$outcome, fits$stage, sprintf("model_%02d.qs", fits$model_index))
    sizes <- file.info(resolved)
    if (!identical(tolower(resolved), tolower(expected_paths)) || any(sizes$isdir) || any(sizes$size <= 0) ||
        !identical(file.info(initial$fit_path)$size, as.numeric(initial_manifest$bytes))) {
      stop("A fit left the permitted archive location or changed size for ", family, ".")
    }
    deletion[[family]] <- data.frame(fits[c("family", "outcome", "model_index", "stage")],
      disposition = ifelse(fits$fit_path %in% accepted$fit_path, "accepted", "superseded"),
      path = relative_path(resolved), bytes = sizes$size, recorded_sha256 = fits$sha256)

    if (nrow(report) != 4L || !all(report$family == family) || !setequal(report$outcome, iad_outcomes())) {
      stop("Reporting is incomplete for ", family, ".")
    }
    retained <- c(retained, list(
      retain(c(initial_path, initial_manifest_path, accepted_path, integrated_path, cases_path,
        reviewed_path, report_path, paths("combination", "stacking_weights.csv")), "interaction audit"),
      retain(report$file_manifest, "report manifests", report$sha256),
      retain(integrated$loo_path, "integrated LOO", integrated$loo_sha256)))
    for (i in seq_len(nrow(report))) {
      manifest <- iad_read_csv(report$file_manifest[[i]])
      outcome <- report$outcome[[i]]
      required <- file.path(interface, family, outcome, c("conditional_predictions.png", "dose_response.csv.gz",
        "profile_support.csv", "interaction_contrasts.csv.gz"))
      if (!all(required %in% manifest$path) || anyDuplicated(manifest$path)) stop("Incomplete report for ", family, "/", outcome, ".")
      retained <- c(retained, list(retain(manifest$path, "report artifacts", manifest$sha256)))
      checks_path <- manifest$path[basename(manifest$path) == "verification_checks.csv.gz"]
      if (length(checks_path) != 1L) stop("Missing report verification checks.")
      checks <- iad_read_csv(checks_path)
      if (!nrow(checks) || !all(checks$passed)) stop("Report verification failed for ", family, "/", outcome, ".")
    }

    if (anyDuplicated(case_key(cases)) || anyDuplicated(case_key(reviewed)) ||
        !setequal(case_key(cases), case_key(reviewed)) || !all(reviewed$diagnostic_pass) ||
        !all(key(cases) %in% expected_keys)) stop("Incomplete exact LOO review for ", family, ".")
    final <- file.path(archive, family, accepted$outcome, "loo", "final", sprintf("model_%02d.qs", accepted$model_index))
    for (i in seq_along(final)) {
      result <- qs::qread(final[[i]])
      iad_validate_loo(result)
      if (!identical(key(result), key(accepted[i, ])) || !identical(result$fit_sha256, accepted$sha256[[i]]) ||
          !setequal(result$exact_refit_summary$observation, cases$observation[key(cases) == key(result)])) {
        stop("Final LOO provenance differs for ", final[[i]])
      }
    }
    exact_name <- sprintf("model_%02d_observation_%04d.qs", cases$model_index, cases$observation)
    exact_root <- file.path(archive, family, cases$outcome, "loo")
    retry <- reviewed[0L, ]
    if ("exact_stage" %in% names(reviewed)) {
      retry <- reviewed[!is.na(reviewed$exact_stage) & nzchar(reviewed$exact_stage), ]
    }
    if (nrow(retry) && !all(grepl("^remedy_[1-9][0-9]*$", retry$exact_stage))) stop("Invalid exact remedy stage.")
    exact_checkpoints <- c(file.path(exact_root, "exact_checkpoints", exact_name),
      if (nrow(retry)) file.path(archive, family, retry$outcome, "loo", "exact_checkpoints", retry$exact_stage,
        sprintf("model_%02d_observation_%04d.qs", retry$model_index, retry$observation)))
    retained <- c(retained, list(retain(final, "final LOO"),
      retain(file.path(exact_root, "exact", exact_name), "exact LOO corrections"),
      retain(exact_checkpoints, "compact exact checkpoints")))
  }

  recovery_path <- iad_path("recovery", "accepted_manifest.csv")
  recovery <- iad_read_csv(recovery_path)
  if (nrow(recovery) != 72L || !all(recovery$diagnostic_pass) || !all(file.exists(recovery$fit_path))) {
    stop("Recovered primary fits are incomplete.")
  }
  retained <- c(retained, list(retain(recovery_path, "recovery audit"),
    retain(recovery$loo_path, "recovered primary final LOO", recovery$loo_sha256)))
  retained <- do.call(rbind, retained)
  deletion <- do.call(rbind, deletion)
  if (anyDuplicated(tolower(deletion$path)) || any(tolower(deletion$path) %in% tolower(retained$path))) {
    stop("The proposed deletion list duplicates or overlaps retained files.")
  }
  summary <- do.call(rbind, lapply(split(deletion, deletion$family), function(rows) {
    data.frame(family = rows$family[[1L]], accepted = sum(rows$disposition == "accepted"),
      superseded = sum(rows$disposition == "superseded"), files = nrow(rows),
      reclaimable_gib = sum(rows$bytes) / 1024^3)
  }))
  directory <- iad_path("storage", "interaction_pruning")
  iad_csv(deletion, file.path(directory, "deletion_file_manifest.csv"))
  iad_csv(retained, file.path(directory, "retained_file_manifest.csv"))
  iad_csv(summary, file.path(directory, "summary.csv"))
  print(summary, row.names = FALSE)
  cat(sprintf("Plan lists %d fit files (%.2f GiB). Verified %d retained file entries. No files deleted.\n",
    nrow(deletion), sum(summary$reclaimable_gib), nrow(retained)))
  invisible(list(deletion = deletion, retained = retained, summary = summary))
}

if (sys.nframe() == 0L) iad_plan_interaction_pruning()
