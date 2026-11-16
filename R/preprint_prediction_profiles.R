preprint_median_profile <- function(data) {
  ages <- dplyr::distinct(data, study_no, arm_no, age_centre)$age_centre + 25
  weeks <- dplyr::distinct(data, study_no, time2)$time2
  c(age = stats::median(ages), weeks = stats::median(weeks))
}

read_preprint_profile_sample <- function(bundle) {
  readr::read_csv(file.path(bundle, "analysis_data.csv.gz"),
                 show_col_types = FALSE, guess_max = Inf) |>
    dplyr::mutate(
      resistance_exercise_factor = factor(resistance_exercise_factor, c("yes", "no")),
      study_factor = factor(as.character(study_no))
    )
}

read_preprint_profile_models <- function(bundle) {
  fits <- readr::read_csv(file.path(bundle, "accepted_fit_manifest.csv.gz"), show_col_types = FALSE)
  weights <- readr::read_csv(file.path(bundle, "stacking_weights.csv.gz"), show_col_types = FALSE)
  fits |>
    dplyr::transmute(outcome, model_index, fit_path = sub("^.*(model_archive/)", "\\1", fit_path)) |>
    dplyr::inner_join(weights, by = c("outcome", "model_index")) |>
    dplyr::arrange(outcome, model_index)
}

preprint_reprofile_primary_draws <- function(frozen, record, data, profile) {
  settings <- frozen$settings
  grid <- make_primary_prediction_grid(
    data, protein_min = settings$protein_min, grid_step = settings$grid_step
  )
  if (!identical(grid$keys, frozen$grid) || record$seed != settings$seed) {
    stop("Frozen grid or seed differs for ", record$outcome)
  }

  # For additive moderators a single offset applies to each draw across doses.
  # Check separated doses in both RT conditions before applying the offset.
  anchors <- unlist(lapply(c(0, 1), function(rt) {
    columns <- which(grid$keys$resistance_exercise_code == rt)
    columns[c(1L, ceiling(length(columns) / 2), length(columns))]
  }), use.names = FALSE)
  old <- grid$newdata[anchors, , drop = FALSE]
  typical <- old
  typical$age_centre <- unname(profile[["age"]] - 25)
  typical$duration_centre <- unname(profile[["weeks"]] - 12)
  moderators <- add_primary_standardised_moderators(data)
  centres <- attr(moderators, "moderator_centres")
  scales <- attr(moderators, "moderator_scales")
  typical$age_outcome_z <- (typical$age_centre - centres[["age_centre"]]) / scales[["age_centre"]]
  typical$duration_outcome_z <- (typical$duration_centre - centres[["duration_centre"]]) / scales[["duration_centre"]]
  paired <- dplyr::bind_rows(old, typical)

  models <- record$models
  set.seed(settings$seed)
  counts <- allocate_primary_mixture_models(models$stacking_weight, settings$n_draws)
  allocation <- sample(rep(seq_len(nrow(models)), times = counts))
  if (!identical(as.integer(models$model_index[allocation]),
                 as.integer(frozen$draw_manifest$model_index))) {
    stop("Frozen mixture allocation differs for ", record$outcome)
  }

  offsets <- numeric(settings$n_draws)
  error <- 0
  for (i in which(counts > 0)) {
    positions <- which(allocation == i)
    fit <- qs::qread(models$fit_path[[i]], nthreads = 4L)
    if (fit$outcome != record$outcome || fit$model_index != models$model_index[[i]] ||
        fit$model_name != models$model_name[[i]]) {
      stop("A saved fit differs from the reporting manifest.")
    }
    newdata <- paired
    for (column in intersect(names(newdata), names(fit$model$data))) {
      if (is.factor(fit$model$data[[column]])) {
        newdata[[column]] <- factor(as.character(newdata[[column]]),
                                    levels = levels(fit$model$data[[column]]))
      }
    }
    prediction <- as.matrix(brms::posterior_epred(
      fit$model, newdata = newdata, re_formula = NA, ndraws = length(positions)
    ))
    if (nrow(prediction) != length(positions)) prediction <- t(prediction)
    reference <- prediction[, seq_along(anchors), drop = FALSE]
    shift <- prediction[, length(anchors) + seq_along(anchors), drop = FALSE] - reference
    error <- max(error, abs(reference - frozen$draws[positions, anchors, drop = FALSE]))
    if (!is.finite(error) || error > 1e-10 || max(abs(shift - shift[, 1L])) > 1e-10) {
      stop("Could not replay the frozen additive draws for ", record$outcome,
           " model ", models$model_index[[i]], " (maximum reference error ", error, ").")
    }
    offsets[positions] <- shift[, 1L]
    rm(fit, prediction)
    gc(verbose = FALSE)
  }

  updated <- frozen
  updated$draws <- sweep(frozen$draws, 1L, offsets, "+")
  reference_increments <- frozen$draws - frozen$draws[, 1L]
  updated_increments <- updated$draws - updated$draws[, 1L]
  contrast_error <- max(abs(reference_increments - updated_increments))
  if (contrast_error > 1e-10) stop("Changing the display profile altered dose contrasts.")
  list(
    predictions = summarise_primary_dose_response(updated),
    reference_error = error,
    contrast_error = contrast_error
  )
}

write_preprint_median_primary_figures <- function(audit, arm_audit, data_by_outcome,
                                                 profiles, directory) {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("The patchwork package is required.")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(tmpdir = directory)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  for (i in seq_len(nrow(profiles))) {
    profile <- profiles[i, ]
    outcome <- profile$outcome
    message("Rendering the paired figures for ", outcome)
    for (show_observed in c(FALSE, TRUE)) {
      figure <- make_update_2026_primary_arm_figure(
        outcome, audit, arm_audit, data_by_outcome, include_observed_effects = show_observed
      ) + patchwork::plot_annotation(subtitle = sprintf(
        "Population predictions at median age %g years and follow-up %g weeks",
        profile$median_age_years, profile$median_duration_weeks
      ))
      name <- paste0(outcome, if (show_observed) "-observed.png" else "-arm-conditional.png")
      ggplot2::ggsave(file.path(stage, name), figure, device = ragg::agg_png,
                      width = 8, height = 12, dpi = 300)
    }
  }
  paths <- file.path(directory, c(paste0(profiles$outcome, "-arm-conditional.png"),
                                  paste0(profiles$outcome, "-observed.png")))
  staged <- file.path(stage, basename(paths))
  stopifnot(all(file.exists(staged)), all(file.info(staged)$size > 0))
  fs::file_move(staged, paths)
  paths
}
