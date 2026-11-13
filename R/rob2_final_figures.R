output_dir <- file.path("manuscript", "generated", "reporting", "rob2")
judgments_file <- file.path(output_dir, "rob2_final_judgments.csv")

if (!requireNamespace("robvis", quietly = TRUE)) {
  stop("The robvis package is required. Restore the project renv environment before running this script.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required. Restore the project renv environment before running this script.")
}

judgments <- read.csv(
  judgments_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

judgment_fields <- c("study_label", paste0("D", 1:5, "_final"), "overall_final")
if (!all(judgment_fields %in% names(judgments))) {
  stop("The public RoB 2 judgments must include the study label and final domain judgments.")
}
rob <- judgments[judgment_fields]
names(rob) <- c("Study", paste0("D", 1:5), "Overall")
rob$Weight <- 1

required_columns <- c("Study", "D1", "D2", "D3", "D4", "D5", "Overall", "Weight")
risk_levels <- c("Low", "Some concerns", "High")

if (!identical(names(rob), required_columns)) {
  stop("The robvis input must contain Study, D1-D5, Overall, and Weight in that order.")
}
if (any(!unlist(rob[c("D1", "D2", "D3", "D4", "D5", "Overall")]) %in% risk_levels)) {
  stop("Unexpected RoB 2 judgment found.")
}
if (any(!is.finite(rob$Weight)) || any(rob$Weight <= 0)) {
  stop("All robvis weights must be finite and positive.")
}

robvis_theme <- ggplot2::theme(
  text = ggplot2::element_text(family = "sans", colour = "#222222"),
  plot.title = ggplot2::element_text(face = "bold", size = 13),
  plot.subtitle = ggplot2::element_text(size = 10, colour = "#555555"),
  plot.caption = ggplot2::element_text(size = 8, colour = "#555555", hjust = 0),
  plot.margin = ggplot2::margin(10, 16, 10, 10)
)

traffic_theme <- ggplot2::theme(
  strip.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = 7.5),
  strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 7.5),
  legend.position = "bottom",
  legend.justification = "center",
  legend.direction = "horizontal",
  legend.margin = ggplot2::margin(t = 4, r = 0, b = 0, l = 0),
  plot.caption = ggplot2::element_text(size = 7, colour = "#555555", hjust = 0)
)

# RoB 2 traffic-light figure ----

# Use the validated trial-level judgments and shared theme to display the five
# domains and overall judgment for each index result.
make_traffic <- function(data, subtitle = NULL, point_size = 8) {
  robvis::rob_traffic_light(
    data = data,
    tool = "ROB2",
    colour = "cochrane",
    psize = point_size,
    quiet = FALSE
  ) +
    ggplot2::labs(
      title = "Risk of bias for the prespecified index result",
      subtitle = subtitle,
      caption = sprintf(
        paste0(
          "Domains: D1 randomization process; D2 deviations from intended interventions; ",
          "D3 missing outcome data; D4 outcome measurement; D5 selection of the reported result.\n",
          "Cochrane RoB 2; one prespecified index result per trial; n = %d trials."
        ),
        nrow(rob)
      )
    ) +
    robvis_theme +
    traffic_theme
}

traffic_height <- max(14, 0.25 * nrow(rob) + 3.5)
traffic_plot <- make_traffic(
  rob,
  subtitle = "Effect of assignment to intervention",
  point_size = 8
)

ggplot2::ggsave(
  filename = file.path(output_dir, "rob2_traffic_light.png"),
  plot = traffic_plot,
  width = 9.5,
  height = traffic_height,
  units = "in",
  dpi = 240,
  bg = "white",
  limitsize = FALSE
)

# A paginated PDF is easier to inspect and use as supplementary material.
dir.create("outputs", showWarnings = FALSE)
page_size <- 30L
page_starts <- seq.int(1L, nrow(rob), by = page_size)
grDevices::pdf(
  file.path("outputs", "rob2_traffic_light.pdf"),
  width = 9.5,
  height = 11,
  onefile = TRUE,
  useDingbats = FALSE
)
for (page_index in seq_along(page_starts)) {
  first_row <- page_starts[[page_index]]
  last_row <- min(first_row + page_size - 1L, nrow(rob))
  page_data <- rob[first_row:last_row, , drop = FALSE]
  print(make_traffic(
    page_data,
    subtitle = sprintf(
      "Effect of assignment to intervention — page %d of %d (trials %d-%d)",
      page_index,
      length(page_starts),
      first_row,
      last_row
    ),
    point_size = 11
  ))
}
grDevices::dev.off()

summary_plot <- robvis::rob_summary(
  data = rob,
  tool = "ROB2",
  overall = TRUE,
  weighted = FALSE,
  colour = "cochrane",
  quiet = FALSE
) +
  ggplot2::labs(
    title = "Risk-of-bias judgments across included trials",
    subtitle = "Unweighted proportions for the prespecified index result",
    caption = sprintf("Cochrane RoB 2; n = %d trials.", nrow(rob))
  ) +
  robvis_theme

ggplot2::ggsave(
  filename = file.path(output_dir, "rob2_proportions.png"),
  plot = summary_plot,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  filename = file.path(output_dir, "rob2_proportions.pdf"),
  plot = summary_plot,
  width = 11,
  height = 7,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat(
  "Created robvis/ggplot2 RoB 2 figures in",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)
