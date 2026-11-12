# One scheduler for unattended initial fits, with independent dose-family gates.
iad_check_family_pilot <- function(pilot, family, checkpoints,
    specifications = iad_specifications(checkpoints)) {
  expected <- dplyr::bind_rows(iad_pilot_specs(specifications))
  expected <- expected[expected$family == family, ]
  required <- c("family", "outcome", "model_index", "parent_index", "diagnostic_pass", "data_sha256")
  if (!all(required %in% names(pilot))) stop("Pilot diagnostics are incomplete: ", family)
  observed <- pilot[pilot$family == family, ]
  key <- function(x) paste(x$outcome, x$model_index)
  index <- match(key(expected), key(observed))
  if (nrow(expected) != 8L || nrow(observed) != 8L || anyNA(index) ||
      anyDuplicated(key(observed)) || anyNA(observed$diagnostic_pass) ||
      !all(observed$diagnostic_pass) || any(observed$parent_index[index] != expected$parent_index)) {
    stop("Pilot gate failed for ", family, "; other dose families may continue independently.")
  }
  for (i in seq_len(nrow(observed))) {
    if (observed$data_sha256[[i]] != iad_fingerprint(checkpoints[[family]]$data[[observed$outcome[[i]]]])) {
      stop("Pilot data fingerprint differs: ", family)
    }
  }
  invisible(TRUE)
}

iad_family_pilot_gate <- function(pilot, family, checkpoints) {
  iad_check_family_pilot(pilot, family, checkpoints)
  iad_csv(pilot[pilot$family == family, ], iad_path(family, "diagnostics", "pilot", "passed_gate.csv"))
}

iad_disk_free_gib <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (.Platform$OS.type == "windows") {
    command <- paste0("([System.IO.DriveInfo]::new((Get-Item -LiteralPath '",
      gsub("'", "''", path, fixed = TRUE), "').PSDrive.Root)).AvailableFreeSpace")
    value <- suppressWarnings(as.numeric(system2("powershell.exe",
      c("-NoProfile", "-NonInteractive", "-Command", shQuote(command)), stdout = TRUE)))
  } else {
    output <- system2("df", c("-Pk", shQuote(path)), stdout = TRUE)
    fields <- strsplit(trimws(tail(output, 1L)), "[[:space:]]+")[[1L]]
    value <- suppressWarnings(as.numeric(fields[[4L]]) * 1024)
  }
  value <- value[is.finite(value)]
  if (length(value) != 1L || value < 0) stop("Could not determine free space on the project drive.")
  value / 1024^3
}

iad_batch_storage_plan <- function(specifications, checkpoints, families, available_gib,
    multiplier = 1.5, reserve_gib = 20) {
  if (!length(families) || anyDuplicated(families) || !all(families %in% iad_families()) ||
      length(available_gib) != 1L || !is.finite(available_gib) || available_gib < 0 ||
      multiplier < 1.5 || reserve_gib < 20) stop("Invalid batch storage settings.")
  records <- Filter(function(s) s$family %in% families, specifications)
  rows <- purrr::map_dfr(records, function(s) {
    parents <- checkpoints[[s$family]]$parents
    parent <- parents[parents$outcome == s$outcome & parents$model_index == s$parent_index, ]
    if (nrow(parent) != 1L || !is.finite(parent$expected_bytes) || parent$expected_bytes <= 0) {
      stop("Missing parent size for batch storage planning.")
    }
    path <- iad_path(s$family, s$outcome, s$stage, sprintf("model_%02d.qs", s$model_index), archive = TRUE)
    tibble::tibble(family = s$family, already_saved = file.exists(path), parent_bytes = parent$expected_bytes)
  })
  expected <- c(relative = 144L, total = 96L, supplement = 96L)[families]
  if (nrow(rows) != sum(expected) || !setequal(rows$family, families)) stop("Incomplete batch model registry.")
  summary <- rows |> dplyr::group_by(family) |> dplyr::summarise(candidates = dplyr::n(),
    saved = sum(already_saved), remaining = sum(!already_saved),
    estimated_remaining_gib = sum(parent_bytes[!already_saved]) / 1024^3,
    planning_allowance_gib = multiplier * estimated_remaining_gib, .groups = "drop")
  summary$available_gib <- available_gib
  summary$reserve_gib <- reserve_gib
  summary$total_required_gib <- sum(summary$planning_allowance_gib) + reserve_gib
  summary$storage_pass <- available_gib >= summary$total_required_gib
  summary
}

iad_batch_storage_gate <- function(specifications, checkpoints) {
  families <- strsplit(Sys.getenv("PROTEIN_INTERACTION_FAMILIES", "relative,total,supplement"), ",", fixed = TRUE)[[1L]]
  plan <- iad_batch_storage_plan(specifications, checkpoints, families, iad_disk_free_gib())
  path <- iad_csv(plan, iad_path("preparation", "batch_storage_gate.csv"))
  if (!all(plan$storage_pass)) stop(sprintf(
    "Insufficient space for the queued interaction fits: %.1f GiB free; %.1f GiB required including reserve. Review batch_storage_gate.csv and completed-family pruning plans. Nothing was deleted.",
    plan$available_gib[[1L]], plan$total_required_gib[[1L]]))
  path
}

iad_family_fit_checkpoint <- function(paths, specifications, family) {
  paths <- unname(unlist(paths))
  specs <- dplyr::bind_rows(specifications)
  if (any(specs$family != family) || nrow(specs) != if (family == "relative") 144L else 96L) {
    stop("Incomplete initial-fit registry: ", family)
  }
  expected <- vapply(specifications, function(s) iad_path(s$family, s$outcome, s$stage,
    sprintf("model_%02d.qs", s$model_index), archive = TRUE), character(1))
  normalise <- function(x) normalizePath(x, winslash = "/", mustWork = TRUE)
  if (length(paths) != length(expected) || anyDuplicated(normalise(paths)) ||
      !setequal(normalise(paths), normalise(expected)) || any(file.info(paths)$size <= 0)) {
    stop("Initial fits are incomplete: ", family)
  }
  specs$fit_path <- normalise(expected); specs$bytes <- file.info(expected)$size
  specs$status <- "saved; full diagnostic review pending"
  iad_csv(specs, iad_path(family, "checkpoints", "initial_fits.csv"))
}

iad_family_reporting_checkpoint <- function(paths, family) {
  paths <- unname(unlist(paths))
  manifests <- paths[basename(paths) == "file_manifest.csv"]
  outcomes <- basename(dirname(dirname(manifests)))
  if (length(manifests) != 4L || anyDuplicated(outcomes) || !setequal(outcomes, iad_outcomes())) {
    stop("Incomplete family reporting: ", family)
  }
  for (path in manifests) {
    manifest <- iad_read_csv(path)
    for (i in seq_len(nrow(manifest))) iad_check_hash(manifest$path[[i]], manifest$sha256[[i]])
  }
  iad_csv(tibble::tibble(family = family, outcome = outcomes, file_manifest = manifests,
    sha256 = vapply(manifests, iad_hash, character(1))), iad_path(family, "checkpoints", "reporting.csv"))
}
