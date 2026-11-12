# Run from the project root, adding --archived after restoring archived inputs.
manifest <- yaml::read_yaml("config/data-manifest.yml")
inputs <- manifest$bundled
if ("--archived" %in% commandArgs(trailingOnly = TRUE)) {
  inputs <- c(inputs, manifest$archived)
}
for (input in inputs) {
  path <- input$path
  if (!file.exists(path)) stop("Missing input: ", path, ". See REPRODUCIBILITY.md.")
  if (isTRUE(input$text_lf)) {
    value <- readChar(path, nchars = file.info(path)$size, useBytes = TRUE)
    value <- gsub("\r\n", "\n", value, fixed = TRUE, useBytes = TRUE)
    sha256 <- digest::digest(charToRaw(value), algo = "sha256", serialize = FALSE)
  } else {
    sha256 <- digest::digest(path, algo = "sha256", file = TRUE)
  }
  if (sha256 != input$sha256) stop("Input checksum differs: ", path)
}
cat("Verified ", length(inputs), " input checksums.\n", sep = "")
