# Read-only design audit before choosing training-status interaction models.
# No pipeline stages, model fits, prior calibration or workbook edits.
training_interaction_audit <- function(root = ".") {
  folders <- c(relative = "primary_update_2026", total = "total_primary_update_2026",
               supplement = "supplement_primary_update_2026")
  read_verified <- function(directory, artifact) {
    manifest <- read.csv(file.path(directory, "manifest.csv"), stringsAsFactors = FALSE)
    record <- manifest[manifest$artifact == artifact, , drop = FALSE]
    if (nrow(record) != 1L || basename(record$file) != record$file) stop("Invalid audit source: ", artifact)
    path <- file.path(directory, record$file)
    if (!identical(digest::digest(file = path, algo = "sha256"), record$sha256)) stop("Audit source hash changed.")
    read.csv(gzfile(path), stringsAsFactors = FALSE)
  }
  normalise <- function(x) {
    value <- tolower(trimws(as.character(x)))
    value[is.na(value) | !nzchar(value)] <- "not reported"
    if (any(!value %in% c("trained", "untrained", "active", "not reported"))) {
      stop("Unreviewed training-status label: ", paste(setdiff(unique(value),
        c("trained", "untrained", "active", "not reported")), collapse = ", "))
    }
    value
  }
  rows <- lapply(names(folders), function(family) {
    directory <- file.path(root, "results", "verification", folders[[family]])
    data <- read_verified(directory, "analysis_data")
    loo <- read_verified(directory, if (family == "relative") "loo_final_manifest" else "final_loo_manifest")
    if (family != "relative") data <- data[is.finite(data$mean_body_mass_kg), , drop = FALSE]
    for (outcome in unique(data$analysis_outcome)) {
      entry <- loo[loo$outcome == outcome, , drop = FALSE][1L, ]
      if (!identical(digest::digest(file = entry$final_loo_path, algo = "sha256"), entry$final_sha256)) {
        stop("Final LOO hash mismatch.")
      }
      reference <- qs::qread(entry$final_loo_path)$observation_manifest
      current <- data[data$analysis_outcome == outcome, , drop = FALSE]
      columns <- c("study_no", "arm_no", "outcome_no", "yi", "vi")
      for (column in columns) {
        first <- if (column %in% c("yi", "vi")) as.numeric(reference[[column]]) else as.character(reference[[column]])
        second <- if (column %in% c("yi", "vi")) as.numeric(current[[column]]) else as.character(current[[column]])
        if (!isTRUE(all.equal(first, second, tolerance = 1e-13, check.attributes = FALSE))) {
          stop("Common-sample/order mismatch in training audit: ", family, "/", outcome)
        }
      }
    }
    data$family <- family
    data$training_status_audit <- normalise(data$trained_or_untrained)
    data
  })
  data <- dplyr::bind_rows(rows)
  if (any(!is.finite(data$age_centre)) || any(!is.finite(data$duration_centre))) stop("Unexpected missing age/duration.")
  arm_consistency <- data |>
    dplyr::group_by(family, study_no, arm_no) |>
    dplyr::summarise(statuses = dplyr::n_distinct(training_status_audit),
      rt_conditions = dplyr::n_distinct(resistance_exercise_code), .groups = "drop")
  if (any(arm_consistency$statuses != 1L) || any(arm_consistency$rt_conditions != 1L)) stop("Inconsistent arm status/RT coding.")
  scope <- c("family", "analysis_outcome")
  status_counts <- data |>
    dplyr::group_by(family, analysis_outcome, training_status_audit) |>
    dplyr::summarise(studies = dplyr::n_distinct(study_no), arms = dplyr::n_distinct(arm_no),
      contrasts = dplyr::n(), .groups = "drop")
  rt_counts <- data |>
    dplyr::group_by(family, analysis_outcome, training_status_audit, resistance_exercise_code) |>
    dplyr::summarise(studies = dplyr::n_distinct(study_no), arms = dplyr::n_distinct(arm_no),
      contrasts = dplyr::n(), distinct_doses = dplyr::n_distinct(protein_intake_raw),
      dose_min = min(protein_intake_raw), dose_max = max(protein_intake_raw),
      age_min = min(age_centre + 25), age_max = max(age_centre + 25),
      duration_min = min(duration_centre + 12), duration_max = max(duration_centre + 12), .groups = "drop")
  scenarios <- lapply(c("reported_three_level", "binary_only"), function(scenario) {
    included <- if (scenario == "reported_three_level") c("trained", "untrained", "active") else c("trained", "untrained")
    data |>
      dplyr::group_by(family, analysis_outcome) |>
      dplyr::summarise(original_contrasts = dplyr::n(), original_studies = dplyr::n_distinct(study_no),
        original_arms = dplyr::n_distinct(arm_no),
        retained_contrasts = sum(training_status_audit %in% included),
        retained_studies = dplyr::n_distinct(study_no[training_status_audit %in% included]),
        retained_arms = dplyr::n_distinct(arm_no[training_status_audit %in% included]), .groups = "drop") |>
      dplyr::mutate(scenario = scenario, excluded_contrasts = original_contrasts - retained_contrasts,
        original_candidates_to_refit = ifelse(excluded_contrasts > 0L, ifelse(family == "relative", 26L, 18L), 0L))
  })
  scenarios <- dplyr::bind_rows(scenarios)
  # Rank screens catch structural non-identification, not statistical power.
  rank <- lapply(split(data, interaction(data$family, data$analysis_outcome, drop = TRUE)), function(current) {
    dplyr::bind_rows(lapply(c("reported_three_level", "binary_only"), function(scenario) {
      retained <- current[current$training_status_audit %in%
        if (scenario == "reported_three_level") c("trained", "untrained", "active") else c("trained", "untrained"), ]
      retained$status <- factor(retained$training_status_audit)
      retained$rt <- retained$resistance_exercise_code
      retained$p <- retained$protein_intake
      formulas <- list(pooled_status_interaction = ~ age_centre + duration_centre + rt + p * status + p:rt,
        rt_specific_status_interaction = ~ age_centre + duration_centre + p * status * rt)
      dplyr::bind_rows(lapply(names(formulas), function(name) {
        if (nlevels(retained$status) < 2L) return(data.frame(family = current$family[[1L]],
          outcome = current$analysis_outcome[[1L]], scenario = scenario, design = name,
          columns = NA_integer_, rank = NA_integer_, full_rank = FALSE))
        matrix <- stats::model.matrix(formulas[[name]], retained)
        data.frame(family = current$family[[1L]], outcome = current$analysis_outcome[[1L]],
          scenario = scenario, design = name, columns = ncol(matrix), rank = qr(matrix)$rank,
          full_rank = qr(matrix)$rank == ncol(matrix))
      }))
    }))
  })
  unknown <- data[data$training_status_audit == "not reported", ] |>
    dplyr::group_by(family, analysis_outcome, study_no, arm_no, resistance_exercise_code) |>
    dplyr::summarise(contrasts = dplyr::n(), .groups = "drop")
  overall <- data |>
    dplyr::group_by(family, training_status_audit) |>
    dplyr::summarise(studies = dplyr::n_distinct(study_no), arms = dplyr::n_distinct(arm_no),
      contrasts = dplyr::n(), .groups = "drop")
  list(overall_status = overall, status_by_outcome = status_counts, status_by_outcome_rt = rt_counts,
    sample_scenarios = scenarios, design_rank = dplyr::bind_rows(rank), unknown_status_arms = unknown,
    checks = data.frame(check = c("source inputs and final LOO hashes pass", "source order matches final common samples",
      "age and duration complete", "arm status and RT coding consistent"), passed = TRUE))
}

write_training_interaction_audit <- function(audit, output_dir = "outputs/training_status_interactions") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  vapply(names(audit), function(name) {
    path <- file.path(output_dir, paste0(name, ".csv"))
    utils::write.csv(audit[[name]], path, row.names = FALSE, na = "")
    path
  }, character(1))
}
