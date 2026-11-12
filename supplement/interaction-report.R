interaction_exposures <- c(
  relative = "Relative Total Protein Intake",
  total = "Absolute Total Protein Intake",
  supplement = "Absolute Protein Supplementation"
)

interaction_report_grid <- function(report) {
  tidyr::expand_grid(
    exposure = names(interaction_exposures),
    outcome = names(report_outcome_labels)
  ) |>
    dplyr::mutate(
      relative_path = file.path("..", "manuscript", "generated", report, exposure, outcome),
      path = file.path(report_project_root(), "manuscript", "generated", report, exposure, outcome)
    )
}

interaction_csv <- function(path, file) {
  readr::read_csv(file.path(path, file), show_col_types = FALSE, progress = FALSE)
}

interaction_links <- function(path, files) {
  links <- sprintf("[%s](%s)", names(files), file.path(path, unname(files)))
  cat("\n\n", paste(links, collapse = " · "), "\n\n", sep = "")
}

age_duration_curves <- function(data, exposure) {
  ensemble_labels <- c(
    "Primary only (matched frozen sample)" = "Primary candidates",
    "Primary plus age/duration interactions" = "Primary and interaction candidates"
  )
  profile_labels <- c(
    age_q25 = "Age (25th Percentile)",
    reference_age = "Reference Age",
    age_q75 = "Age (75th Percentile)",
    duration_q25 = "Duration (25th Percentile)",
    reference_duration = "Reference Duration",
    duration_q75 = "Duration (75th Percentile)"
  )
  rt_labels <- unname(report_rt_display[c("RT", "No RT")])
  data <- data |>
    dplyr::mutate(
      ensemble = factor(ensemble_labels[ensemble], levels = unname(ensemble_labels)),
      rt = unname(report_rt_display[report_rt_label(resistance_exercise)])
    )
  reference <- dplyr::filter(data, profile == "reference")
  quartiles <- dplyr::bind_rows(
    dplyr::filter(data, profile %in% names(profile_labels)),
    dplyr::mutate(reference, profile = "reference_age"),
    dplyr::mutate(reference, profile = "reference_duration")
  ) |>
    dplyr::mutate(
      profile_label = unname(profile_labels[profile]),
      panel = factor(
        paste(profile_label, rt, sep = "\n"),
        levels = as.vector(outer(unname(profile_labels), rt_labels, paste, sep = "\n"))
      )
    )
  typical <- dplyr::filter(data, profile == "typical") |>
    dplyr::mutate(panel = factor(paste("Typical Profile", rt, sep = "\n"),
                                levels = paste("Typical Profile", rt_labels, sep = "\n")))

  plot_profiles <- function(x, columns) {
    ggplot2::ggplot(x, ggplot2::aes(protein_intake, estimate, colour = ensemble, fill = ensemble)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.25) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = conf_low, ymax = conf_high), alpha = 0.13, colour = NA) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::facet_wrap(ggplot2::vars(panel), ncol = columns) +
      ggplot2::scale_colour_manual(values = c("black", "#0072B2")) +
      ggplot2::scale_fill_manual(values = c("black", "#0072B2")) +
      ggplot2::labs(x = report_dose_axes[[exposure]], y = "Standardised Mean Change",
                    colour = NULL, fill = NULL) +
      report_theme()
  }

  patchwork::wrap_plots(
    plot_profiles(quartiles, 3),
    plot_profiles(typical, 2),
    ncol = 1, heights = c(4, 1.25), guides = "collect"
  ) & ggplot2::theme(legend.position = "bottom")
}

training_status_curves <- function(data, profiles, exposure, reference_dose) {
  profile_labels <- stats::setNames(
    sprintf("%s Profile\n%s Years, %s Weeks",
            ifelse(profiles$profile == "reference", "Reference", "Typical"),
            round(profiles$age, 1), round(profiles$weeks, 1)),
    profiles$profile
  )
  data <- data |>
    dplyr::mutate(
      training_status = factor(training_status, levels = c("trained", "untrained"),
                               labels = c("Trained", "Untrained")),
      profile = factor(profile, levels = c("reference", "typical"),
                       labels = unname(profile_labels[c("reference", "typical")]))
    )
  ggplot2::ggplot(data, ggplot2::aes(dose, estimate, colour = training_status, fill = training_status)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = reference_dose, colour = "grey65", linewidth = 0.25,
                       linetype = "dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf_low, ymax = conf_high), alpha = 0.13, colour = NA) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(ggplot2::vars(profile), nrow = 1) +
    ggplot2::scale_colour_manual(values = c("black", "#0072B2")) +
    ggplot2::scale_fill_manual(values = c("black", "#0072B2")) +
    ggplot2::labs(
      x = report_dose_axes[[exposure]], y = "Outcome Contrast From Reference Dose (SMD)",
      colour = "Reported Training Status", fill = "Reported Training Status"
    ) +
    report_theme()
}
