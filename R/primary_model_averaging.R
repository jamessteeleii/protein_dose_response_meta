# Corrected primary-model stacking and posterior-mixture predictions.
#
# The functions in this file deliberately accept completed `loo` objects. They
# never call brms::loo() and do not use the earlier conditional-LOO helpers.

validate_primary_final_loo_results <- function(final_loo_results) {
  expected_outcomes <- c(
    "whole_body", "hypertrophy", "strength", "function"
  )

  if (
    !is.list(final_loo_results) ||
      !identical(names(final_loo_results), expected_outcomes)
  ) {
    stop(
      "`final_loo_results` must be `primary_psis_loo_final_audit$results` ",
      "with outcomes in their audited order.",
      call. = FALSE
    )
  }

  purrr::walk2(
    final_loo_results,
    expected_outcomes,
    function(results, outcome) {
      model_indices <- unname(purrr::map_int(results, "model_index"))
      model_names <- unname(purrr::map_chr(results, "model_name"))

      if (
        length(results) != 26L ||
          !identical(model_indices, seq_len(26L)) ||
          anyDuplicated(model_names) ||
          !all(purrr::map_chr(results, "outcome") == outcome) ||
          !all(purrr::map_lgl(results, ~ inherits(.x$loo, "loo")))
      ) {
        stop(
          "The corrected LOO results for '", outcome,
          "' are incomplete or out of model order.",
          call. = FALSE
        )
      }

      pointwise_rows <- purrr::map_int(
        results,
        ~ nrow(.x$loo$pointwise)
      )
      if (length(unique(pointwise_rows)) != 1L) {
        stop(
          "Corrected LOO objects for '", outcome,
          "' do not use a common observation set.",
          call. = FALSE
        )
      }
    }
  )

  invisible(final_loo_results)
}

compute_primary_stacking_weights <- function(final_loo_results) {
  validate_primary_final_loo_results(final_loo_results)

  purrr::imap(
    final_loo_results,
    function(results, outcome) {
      loo_objects <- purrr::map(results, "loo")
      model_names <- purrr::map_chr(results, "model_name")
      names(loo_objects) <- model_names

      weights <- loo::loo_model_weights(
        loo_objects,
        method = "stacking"
      )
      weights <- as.numeric(weights)
      names(weights) <- model_names

      if (
        length(weights) != 26L ||
          any(!is.finite(weights)) ||
          any(weights < -sqrt(.Machine$double.eps)) ||
          !isTRUE(all.equal(sum(weights), 1, tolerance = 1e-8))
      ) {
        stop(
          "Stacking returned invalid weights for '", outcome, "'.",
          call. = FALSE
        )
      }

      # Optimisation can return negligible negative round-off at the boundary.
      weights <- pmax(weights, 0)
      weights / sum(weights)
    }
  )
}

summarise_primary_stacking_weights <- function(stacking_weights) {
  purrr::imap_dfr(
    stacking_weights,
    ~ tibble::tibble(
      outcome = .y,
      model_index = seq_along(.x),
      model_name = names(.x),
      stacking_weight = as.numeric(.x)
    )
  ) |>
    dplyr::arrange(outcome, dplyr::desc(stacking_weight), model_index)
}

validate_primary_models_and_weights <- function(models, weights) {
  if (is.null(names(models)) || anyNA(names(models)) || any(names(models) == "")) {
    stop("`models` must be a named list.", call. = FALSE)
  }
  if (anyDuplicated(names(models))) {
    stop("Model names must be unique.", call. = FALSE)
  }
  if (!all(purrr::map_lgl(models, inherits, what = "brmsfit"))) {
    stop("Every object in `models` must be a brmsfit.", call. = FALSE)
  }
  if (
    length(weights) != length(models) ||
      is.null(names(weights)) ||
      !identical(names(weights), names(models))
  ) {
    stop(
      "Named stacking weights must match the candidate models exactly and ",
      "in audited model order.",
      call. = FALSE
    )
  }
  if (
    any(!is.finite(weights)) ||
      any(weights < 0) ||
      !isTRUE(all.equal(sum(weights), 1, tolerance = 1e-8))
  ) {
    stop("`weights` must be non-negative and sum to one.", call. = FALSE)
  }

  invisible(TRUE)
}

get_primary_protein_median <- function(data) {
  if (!"protein_median" %in% names(data)) {
    stop("`data` must contain `protein_median`.", call. = FALSE)
  }
  protein_median <- unique(stats::na.omit(data$protein_median))
  if (length(protein_median) != 1L || !is.finite(protein_median)) {
    stop(
      "`data$protein_median` must contain one finite value.",
      call. = FALSE
    )
  }
  protein_median[[1L]]
}

make_primary_prediction_grid <- function(
    data,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.01
) {
  if (!"protein_intake_raw" %in% names(data)) {
    stop("`data` must contain `protein_intake_raw`.", call. = FALSE)
  }
  if (is.null(protein_max)) {
    protein_max <- max(data$protein_intake_raw, na.rm = TRUE)
  }
  if (
    length(protein_min) != 1L || length(protein_max) != 1L ||
      length(grid_step) != 1L ||
      !all(is.finite(c(protein_min, protein_max, grid_step))) ||
      protein_min >= protein_max || grid_step <= 0
  ) {
    stop("The protein prediction-grid limits are invalid.", call. = FALSE)
  }

  protein_median <- get_primary_protein_median(data)
  protein_values <- round(
    seq(protein_min, protein_max, by = grid_step),
    10
  )

  keys <- tidyr::crossing(
    resistance_exercise_code = c(0, 1),
    protein_intake = protein_values
  )
  grouping_variables <- intersect(
    c("study_no", "arm_no", "outcome_no", "study_factor"),
    names(data)
  )
  reference_index <- 1L
  if (length(grouping_variables) > 0L) {
    complete <- stats::complete.cases(
      data[, grouping_variables, drop = FALSE]
    )
    if (!any(complete)) {
      stop("No complete row is available for the prediction grid.", call. = FALSE)
    }
    reference_index <- which(complete)[[1L]]
  }

  newdata <- data[rep(reference_index, nrow(keys)), , drop = FALSE]
  newdata$resistance_exercise_code <- keys$resistance_exercise_code
  newdata$protein_intake <- keys$protein_intake - protein_median
  newdata$protein_intake_raw <- keys$protein_intake

  exercise_label <- ifelse(keys$resistance_exercise_code == 0, "yes", "no")
  if ("resistance_exercise_yes_no" %in% names(newdata)) {
    newdata$resistance_exercise_yes_no <- exercise_label
  }
  if ("resistance_exercise_factor" %in% names(newdata)) {
    newdata$resistance_exercise_factor <- factor(
      exercise_label,
      levels = levels(data$resistance_exercise_factor)
    )
  }
  if ("age_centre" %in% names(newdata)) newdata$age_centre <- 0
  if ("duration_centre" %in% names(newdata)) newdata$duration_centre <- 0
  if ("vi" %in% names(newdata)) newdata$vi <- 0

  # Post-pilot remedies changed the formula-side variable names for some
  # candidates. Recreate those columns with exactly the transformations used
  # at fitting time while retaining age_centre = duration_centre = 0 as the
  # common prediction profile.
  moderator_data <- add_primary_standardised_moderators(data)
  moderator_centres <- attr(moderator_data, "moderator_centres")
  moderator_scales <- attr(moderator_data, "moderator_scales")
  newdata$age_outcome_z <-
    (newdata$age_centre - moderator_centres[["age_centre"]]) /
      moderator_scales[["age_centre"]]
  newdata$duration_outcome_z <-
    (newdata$duration_centre - moderator_centres[["duration_centre"]]) /
      moderator_scales[["duration_centre"]]
  newdata$protein_intake_05 <- newdata$protein_intake / 0.5

  list(
    newdata = newdata,
    keys = tibble::tibble(
      resistance_exercise_code = keys$resistance_exercise_code,
      resistance_exercise_label = "Resistance Exercise?",
      resistance_exercise = ifelse(
        keys$resistance_exercise_code == 0,
        "YES",
        "NO"
      ),
      protein_intake = keys$protein_intake
    ),
    protein_median = protein_median,
    protein_min = protein_min,
    protein_max = tail(protein_values, 1L),
    grid_step = grid_step
  )
}

allocate_primary_mixture_models <- function(weights, n_draws) {
  expected <- weights * n_draws
  counts <- floor(expected)
  remainder <- n_draws - sum(counts)
  if (remainder > 0L) {
    order_remainder <- order(expected - counts, decreasing = TRUE)
    counts[order_remainder[seq_len(remainder)]] <-
      counts[order_remainder[seq_len(remainder)]] + 1L
  }
  counts
}

build_primary_model_average <- function(
    models,
    weights,
    data,
    protein_min = 0.8,
    protein_max = NULL,
    grid_step = 0.01,
    n_draws = 4000L,
    seed = 1988L
) {
  validate_primary_models_and_weights(models, weights)
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1) {
    stop("`n_draws` must be a positive integer.", call. = FALSE)
  }
  n_draws <- as.integer(n_draws)
  prediction_grid <- make_primary_prediction_grid(
    data = data,
    protein_min = protein_min,
    protein_max = protein_max,
    grid_step = grid_step
  )

  keep <- weights > 0
  models <- models[keep]
  weights <- weights[keep]
  weights <- weights / sum(weights)
  counts <- allocate_primary_mixture_models(weights, n_draws)

  set.seed(seed)
  allocation <- sample(rep(seq_along(models), times = counts))
  averaged_draws <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = nrow(prediction_grid$newdata)
  )

  for (model_index in seq_along(models)) {
    output_rows <- which(allocation == model_index)
    if (length(output_rows) == 0L) next

    model_name <- names(models)[[model_index]]
    draw_matrix <- tryCatch(
      brms::posterior_epred(
        models[[model_index]],
        newdata = prediction_grid$newdata,
        re_formula = NA,
        ndraws = length(output_rows)
      ),
      error = function(error) {
        stop(
          "Prediction extraction failed for '", model_name, "': ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
    draw_matrix <- as.matrix(draw_matrix)
    if (
      ncol(draw_matrix) != nrow(prediction_grid$newdata) &&
        nrow(draw_matrix) == nrow(prediction_grid$newdata)
    ) {
      draw_matrix <- t(draw_matrix)
    }
    if (ncol(draw_matrix) != nrow(prediction_grid$newdata)) {
      stop(
        "Posterior predictions for '", model_name,
        "' do not align with the shared prediction grid.",
        call. = FALSE
      )
    }
    if (nrow(draw_matrix) != length(output_rows)) {
      stop(
        "Posterior prediction sampling for '", model_name,
        "' returned the wrong number of draws.",
        call. = FALSE
      )
    }
    averaged_draws[output_rows, ] <- draw_matrix
  }

  if (anyNA(averaged_draws)) {
    stop("Some posterior-mixture draws were not populated.", call. = FALSE)
  }

  output <- list(
    draws = averaged_draws,
    grid = prediction_grid$keys,
    draw_manifest = tibble::tibble(
      drawid = seq_len(n_draws),
      model_index = which(keep)[allocation],
      model_name = names(models)[allocation]
    ),
    model_weights = weights,
    settings = list(
      protein_median = prediction_grid$protein_median,
      protein_min = prediction_grid$protein_min,
      protein_max = prediction_grid$protein_max,
      grid_step = prediction_grid$grid_step,
      n_draws = n_draws,
      seed = seed
    )
  )
  class(output) <- c("primary_model_average", class(output))
  output
}

summarise_primary_dose_response <- function(
    model_average,
    conf_level = 0.95
) {
  if (!inherits(model_average, "primary_model_average")) {
    stop("`model_average` is not a primary model average.", call. = FALSE)
  }
  alpha <- (1 - conf_level) / 2
  intervals <- matrixStats::colQuantiles(
    model_average$draws,
    probs = c(alpha, 1 - alpha),
    drop = FALSE
  )

  dplyr::bind_cols(
    model_average$grid,
    tibble::tibble(
      estimate = colMeans(model_average$draws),
      conf.low = intervals[, 1L],
      conf.high = intervals[, 2L]
    )
  )
}

compute_primary_heatmap_contrasts <- function(
    model_average,
    data,
    protein_min = 0.8,
    sesoi = 0.10,
    conf_level = 0.95
) {
  if (!inherits(model_average, "primary_model_average")) {
    stop("`model_average` is not a primary model average.", call. = FALSE)
  }
  required_data <- c("protein_intake_raw", "resistance_exercise_code")
  missing_data <- setdiff(required_data, names(data))
  if (length(missing_data) > 0L) {
    stop(
      "`data` is missing required columns: ",
      paste(missing_data, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (
    length(sesoi) != 1L || !is.finite(sesoi) || sesoi <= 0 ||
      length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1
  ) {
    stop("`sesoi` and `conf_level` are invalid.", call. = FALSE)
  }

  grid <- model_average$grid
  observed_support <- data |>
    dplyr::filter(
      is.finite(protein_intake_raw),
      resistance_exercise_code %in% c(0, 1)
    ) |>
    dplyr::group_by(resistance_exercise_code) |>
    dplyr::summarise(
      observed_min = min(protein_intake_raw),
      observed_max = max(protein_intake_raw),
      .groups = "drop"
    ) |>
    dplyr::arrange(resistance_exercise_code)
  if (!identical(
    as.numeric(observed_support$resistance_exercise_code),
    c(0, 1)
  )) {
    stop(
      "Both resistance-training conditions require observed protein intakes.",
      call. = FALSE
    )
  }

  protein_values <- sort(unique(grid$protein_intake))
  observed_min <- min(observed_support$observed_min)
  observed_max <- max(observed_support$observed_max)
  protein_values <- protein_values[
    protein_values >= max(protein_min, observed_min) &
      protein_values <= observed_max
  ]
  if (length(protein_values) < 2L) {
    stop(
      "The model-average grid and observed protein-intake range do not overlap.",
      call. = FALSE
    )
  }

  alpha <- (1 - conf_level) / 2
  results <- vector("list", 2L * length(protein_values))
  result_index <- 1L

  for (exercise_code in c(0, 1)) {
    columns <- which(
      grid$resistance_exercise_code == exercise_code &
        grid$protein_intake %in% protein_values
    )
    columns <- columns[order(grid$protein_intake[columns])]
    exercise_draws <- model_average$draws[, columns, drop = FALSE]

    for (current_index in seq_along(protein_values)) {
      proposed_indices <- seq_along(protein_values)
      differences <- sweep(
        exercise_draws[, proposed_indices, drop = FALSE],
        1L,
        exercise_draws[, current_index],
        "-"
      )
      interval <- matrixStats::colQuantiles(
        differences,
        probs = c(alpha, 1 - alpha),
        drop = FALSE
      )
      direction <- ifelse(
        protein_values[proposed_indices] > protein_values[current_index],
        "Increase",
        ifelse(
          protein_values[proposed_indices] < protein_values[current_index],
          "Decrease",
          "No change"
        )
      )
      probability_meaningful_gain <- colMeans(differences > sesoi)
      probability_negligible <- colMeans(
        differences > -sesoi & differences < sesoi
      )
      probability_meaningful_loss <- colMeans(differences < -sesoi)

      results[[result_index]] <- tibble::tibble(
        resistance_exercise_code = exercise_code,
        resistance_exercise = ifelse(exercise_code == 0, "YES", "NO"),
        rt_condition = factor(
          ifelse(exercise_code == 0, "RT", "No RT"),
          levels = c("No RT", "RT")
        ),
        protein_intake_current = protein_values[current_index],
        protein_intake_proposed = protein_values[proposed_indices],
        direction = direction,
        estimate = colMeans(differences),
        conf_low = interval[, 1L],
        conf_high = interval[, 2L],
        probability_positive = colMeans(differences > 0),
        probability_negative = colMeans(differences < 0),
        probability_meaningful_gain = probability_meaningful_gain,
        probability_negligible = probability_negligible,
        probability_meaningful_loss = probability_meaningful_loss,
        probability_noninferior = colMeans(differences > -sesoi)
      )
      result_index <- result_index + 1L
    }
  }

  contrasts <- dplyr::bind_rows(results) |>
    dplyr::left_join(observed_support, by = "resistance_exercise_code") |>
    dplyr::mutate(
      contrast_estimate = estimate,
      probability_meaningful_benefit = probability_meaningful_gain,
      probability_practical_equivalence = probability_negligible,
      probability_meaningful_detriment = probability_meaningful_loss,
      supported =
        protein_intake_current >= observed_min &
        protein_intake_current <= observed_max &
        protein_intake_proposed >= observed_min &
        protein_intake_proposed <= observed_max,
      support_status = factor(
        ifelse(supported, "Observed range", "Outside RT-specific range"),
        levels = c("Observed range", "Outside RT-specific range")
      )
    )
  output <- list(
    contrasts = contrasts,
    support = observed_support |>
      dplyr::mutate(
        resistance_exercise = ifelse(
          resistance_exercise_code == 0,
          "YES",
          "NO"
        ),
        rt_condition = factor(
          ifelse(resistance_exercise_code == 0, "RT", "No RT"),
          levels = c("No RT", "RT")
        )
      ),
    model_weights = model_average$model_weights,
    settings = utils::modifyList(
      model_average$settings,
      list(
        protein_min = protein_min,
        sesoi = sesoi,
        conf_level = conf_level
      )
    )
  )
  class(output) <- c("primary_heatmap_contrasts", class(output))
  output
}

plot_primary_heatmap_contrasts <- function(
    x,
    title,
    show_contours = TRUE
) {
  if (!inherits(x, "primary_heatmap_contrasts")) {
    stop("`x` is not a primary heatmap-contrast object.", call. = FALSE)
  }
  plot_data <- dplyr::filter(x$contrasts, supported) |>
    dplyr::mutate(plotted_estimate = estimate)
  supported_data <- plot_data
  grid_step <- x$settings$grid_step
  sesoi <- x$settings$sesoi
  fill_limit <- max(abs(supported_data$estimate), na.rm = TRUE)
  if (!is.finite(fill_limit) || fill_limit == 0) fill_limit <- sesoi
  extend_nearby_upper_tenth <- function(limits) {
    rounded_upper <- ceiling(
      (limits[[2L]] - sqrt(.Machine$double.eps)) * 10
    ) / 10
    if (rounded_upper - limits[[2L]] <= 0.03) {
      # Keep a negligible tolerance above the labelled endpoint so ggplot2
      # does not censor a break such as 2.8 because of floating-point error.
      limits[[2L]] <- rounded_upper + sqrt(.Machine$double.eps)
    }
    limits
  }
  contour_names <- c(
    "Meaningful benefit",
    "Practical equivalence",
    "Meaningful detriment"
  )
  contour_specs <- tibble::tibble(
    probability = c(
      "probability_meaningful_benefit",
      "probability_practical_equivalence",
      "probability_meaningful_detriment"
    ),
    contour_region = contour_names
  )
  classification_present <- purrr::map_lgl(
    contour_specs$probability,
    ~ any(supported_data[[.x]] >= 0.95, na.rm = TRUE)
  )
  absent_classifications <- tolower(
    contour_specs$contour_region[!classification_present]
  )
  if (length(absent_classifications) == 0L) {
    classification_note <- paste0(
      "Areas outside the displayed classification regions remain uncertain ",
      "relative to the ROPE."
    )
  } else {
    absent_text <- if (length(absent_classifications) == 1L) {
      absent_classifications
    } else {
      paste(
        paste(utils::head(absent_classifications, -1L), collapse = ", "),
        utils::tail(absent_classifications, 1L),
        sep = " or "
      )
    }
    classification_note <- paste0(
      "No cells met the 95% posterior criterion for ", absent_text, ". ",
      "Areas outside the displayed classification regions remain uncertain ",
      "relative to the ROPE."
    )
  }
  caption_text <- stringr::str_wrap(
    paste0(
      "ROPE = Region of Practical Equivalence (-0.1 to +0.1). ",
      "The solid diagonal is zero contrast (proposed = current). ",
      "Each panel is limited to the observed protein-intake range for ",
      "that RT condition. ",
      classification_note
    ),
    width = 110
  )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = protein_intake_current,
      y = protein_intake_proposed
    )
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = plotted_estimate),
      width = grid_step,
      height = grid_step
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      colour = "black",
      linewidth = 1
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(rt_condition),
      nrow = 1,
      scales = "free",
      drop = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#D55E00",
      mid = "white",
      high = "#0072B2",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      oob = scales::squish,
      na.value = "#E6E6E6",
      name = "Estimated contrast\n(proposed - current)",
      guide = ggplot2::guide_colourbar(
        order = 1,
        direction = "horizontal",
        title.position = "top",
        barwidth = grid::unit(10, "lines")
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_width(0.4),
      limits = extend_nearby_upper_tenth,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_width(0.4),
      limits = extend_nearby_upper_tenth,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::labs(
      x = "Current protein intake (g/kg/day)",
      y = "Proposed protein intake (g/kg/day)",
      title = title,
      subtitle = paste0(
        "Fill shows the mean model-averaged contrast; contours enclose ",
        "regions with at least 95% posterior classification"
      ),
      caption = caption_text
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(1, "lines"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 7, hjust = 1),
      axis.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      legend.box = "horizontal",
      legend.box.just = "center"
    )

  if (show_contours) {
    contour_data <- purrr::map(
      contour_specs$probability,
      function(probability) {
        supported_data |>
          dplyr::group_by(rt_condition) |>
          dplyr::filter(
            min(.data[[probability]], na.rm = TRUE) <= 0.95,
            max(.data[[probability]], na.rm = TRUE) >= 0.95
          ) |>
          dplyr::ungroup()
      }
    )

    for (index in seq_len(nrow(contour_specs))) {
      if (nrow(contour_data[[index]]) == 0L) next
      probability <- contour_specs$probability[[index]]
      contour_region <- contour_specs$contour_region[[index]]
      plot <- plot + ggplot2::geom_contour(
        data = contour_data[[index]],
        mapping = ggplot2::aes(
          x = protein_intake_current,
          y = protein_intake_proposed,
          z = !!rlang::sym(probability),
          linetype = !!contour_region
        ),
        breaks = 0.95,
        colour = "black",
        linewidth = 0.6,
        na.rm = TRUE,
        inherit.aes = FALSE
      )
    }

    boundary_present <- purrr::map_lgl(contour_data, ~ nrow(.x) > 0L)
    present_contours <- contour_names[boundary_present]
    contour_labels <- stats::setNames(
      c(
        sprintf("Benefit (> +%.2f)", sesoi),
        paste0(
          "At least 95% posterior within the\n",
          sprintf(
            "Region of Practical Equivalence (ROPE; -%.1f to +%.1f)",
            sesoi,
            sesoi
          )
        ),
        sprintf("Detriment (< -%.2f)", sesoi)
      ),
      contour_names
    )
    if (length(present_contours) > 0L) {
      plot <- plot +
        ggplot2::scale_linetype_manual(
          values = c(
            "Meaningful benefit" = "dashed",
            "Practical equivalence" = "longdash",
            "Meaningful detriment" = "dotted"
          ),
          limits = present_contours,
          breaks = present_contours,
          labels = contour_labels[present_contours],
          name = NULL
        ) +
        ggplot2::guides(
          linetype = ggplot2::guide_legend(
            order = 2,
            nrow = 1,
            byrow = TRUE,
            override.aes = list(
              alpha = 1,
              colour = "black",
              linewidth = 0.6
            )
          )
        )
    }
  }
  plot
}

plot_primary_predictions_panel <- function(data, predictions) {
  max_no_rt <- max(
    data$protein_intake_raw[data$resistance_exercise_code == 1],
    na.rm = TRUE
  )
  prediction_data <- predictions |>
    dplyr::mutate(
      rt_condition = factor(
        ifelse(resistance_exercise_code == 0, "RT", "No RT"),
        levels = c("No RT", "RT")
      )
    ) |>
    dplyr::filter(rt_condition == "RT" | protein_intake <= max_no_rt)
  point_data <- data |>
    dplyr::mutate(
      rt_condition = factor(
        ifelse(resistance_exercise_code == 0, "RT", "No RT"),
        levels = c("No RT", "RT")
      )
    )

  ggplot2::ggplot(
    prediction_data,
    ggplot2::aes(x = protein_intake, y = estimate)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted",
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(
        x = protein_intake_raw,
        y = yi,
        size = size,
        colour = condition
      ),
      alpha = 0.5
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.5,
      colour = "black",
      linewidth = 0.25
    ) +
    ggplot2::geom_line(
      linewidth = 0.75,
      colour = "black"
    ) +
    ggplot2::scale_colour_manual(
      values = c("#0072B2", "darkorange", "darkgreen")
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_width(width = 0.5, offset = 0)
    ) +
    ggh4x::facet_nested(
      . ~ rt_condition,
      scales = "free_x"
    ) +
    ggh4x::facetted_pos_scales(
      x = list(
        rt_condition == "No RT" ~ ggplot2::scale_x_continuous(
          breaks = scales::breaks_width(0.2)
        ),
        rt_condition == "RT" ~ ggplot2::scale_x_continuous(
          breaks = scales::breaks_width(0.4)
        )
      )
    ) +
    ggplot2::labs(
      x = "Protein Intake (g/kg/day)",
      y = "Standardised Mean Change",
      colour = "Condition",
      title = "Predicted effects",
      subtitle = "Global grand mean predictions and 95% posterior quantile intervals"
    ) +
    ggplot2::guides(size = "none") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

plot_primary_predictions_heatmap <- function(
    data,
    predictions,
    heatmap_contrasts,
    title
) {
  prediction_plot <- plot_primary_predictions_panel(data, predictions)
  heatmap_plot <- plot_primary_heatmap_contrasts(
    heatmap_contrasts,
    title = "Pairwise model-averaged contrasts"
  )

  (prediction_plot / heatmap_plot) +
    patchwork::plot_layout(
      heights = c(1, 1.35)
    ) +
    patchwork::plot_annotation(
      title = paste("Dose-response association of protein intake with", title),
      subtitle = "Bayesian model averaging with LOO-CV stacking weights"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}
