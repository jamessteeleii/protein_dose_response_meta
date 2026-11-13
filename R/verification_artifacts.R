# Verification bundle paths and provenance ----

# Use filesystem-safe artifact names and record the Git commit so an archived
# bundle can be traced to its source code.
verification_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

write_verification_csv <- function(x, path) {
  connection <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(connection), add = TRUE)
  utils::write.csv(x, connection, row.names = FALSE, na = "")
  invisible(path)
}

verification_git_commit <- function() {
  commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(commit) == 1L && nzchar(commit)) commit else NA_character_
}

# Write the named artifacts and a manifest with their sizes and SHA-256
# hashes. The manifest allows downstream readers to verify the exact files
# used.
write_verification_bundle <- function(
  family,
  artifacts,
  descriptions = NULL,
  output_root = here::here("results", "verification")
) {
  stopifnot(
    is.character(family),
    length(family) == 1L,
    nzchar(family),
    is.list(artifacts),
    length(artifacts) > 0L,
    !is.null(names(artifacts)),
    all(nzchar(names(artifacts))),
    !anyDuplicated(names(artifacts))
  )

  if (is.null(descriptions)) {
    descriptions <- stats::setNames(rep("", length(artifacts)), names(artifacts))
  }
  missing_descriptions <- setdiff(names(artifacts), names(descriptions))
  if (length(missing_descriptions) > 0L) {
    descriptions[missing_descriptions] <- ""
  }

  output_dir <- file.path(output_root, verification_slug(family))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  artifact_paths <- character(length(artifacts))
  manifest_rows <- vector("list", length(artifacts))

  for (index in seq_along(artifacts)) {
    artifact_name <- names(artifacts)[[index]]
    artifact <- artifacts[[index]]
    slug <- verification_slug(artifact_name)

    if (is.data.frame(artifact)) {
      path <- file.path(output_dir, paste0(slug, ".csv.gz"))
      write_verification_csv(artifact, path)
      storage_format <- "csv.gz"
      rows <- nrow(artifact)
      columns <- ncol(artifact)
    } else {
      path <- file.path(output_dir, paste0(slug, ".qs"))
      qs::qsave(artifact, path, preset = "high")
      storage_format <- "qs"
      rows <- NA_integer_
      columns <- NA_integer_
    }

    artifact_paths[[index]] <- path
    manifest_rows[[index]] <- data.frame(
      artifact = artifact_name,
      file = basename(path),
      storage_format = storage_format,
      object_class = paste(class(artifact), collapse = "/"),
      rows = rows,
      columns = columns,
      bytes = file.info(path)$size,
      sha256 = digest::digest(file = path, algo = "sha256"),
      description = unname(descriptions[[artifact_name]]),
      stringsAsFactors = FALSE
    )
  }

  manifest <- dplyr::bind_rows(manifest_rows) |>
    dplyr::mutate(
      family = family,
      generated_utc = format(
        Sys.time(),
        tz = "UTC",
        format = "%Y-%m-%dT%H:%M:%SZ"
      ),
      git_commit = verification_git_commit(),
      .before = 1L
    )

  manifest_path <- file.path(output_dir, "manifest.csv")
  utils::write.csv(manifest, manifest_path, row.names = FALSE, na = "")

  readme_path <- file.path(output_dir, "README.md")
  writeLines(
    c(
      paste0("# ", family, " verification bundle"),
      "",
      paste0(
        "Generated from commit `",
        ifelse(is.na(manifest$git_commit[[1L]]), "unknown", manifest$git_commit[[1L]]),
        "` at ", manifest$generated_utc[[1L]], "."
      ),
      "",
      paste0(
        "`manifest.csv` records the type, dimensions, byte size, SHA-256 ",
        "checksum, and purpose of every artifact. CSV tables are gzip-",
        "compressed; non-tabular R objects use the `qs` format."
      ),
      "",
      paste0(
        "These files are sufficient to inspect reported estimates, diagnostics, ",
        "model-comparison results, and plotting data without restoring the ",
        "disposable `_targets` cache. They verify the reporting layer; an ",
        "independent refit is still required to reproduce posterior sampling ",
        "from the raw inputs."
      )
    ),
    readme_path,
    useBytes = TRUE
  )

  normalizePath(
    c(artifact_paths, manifest_path, readme_path),
    winslash = "/",
    mustWork = TRUE
  )
}
