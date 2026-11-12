# Refresh portable manuscript inputs from the accepted analytical outputs.
# Run from the project root after restoring the recorded source bundles.
source("R/preprint_reporting.R")
paths <- write_preprint_report_inputs()
cat("Wrote ", length(paths), " portable reporting inputs.\n", sep = "")
