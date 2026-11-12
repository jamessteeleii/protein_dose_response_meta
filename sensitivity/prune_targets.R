#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
family_arg <- grep("^--family=", args, value = TRUE)
family <- if (length(family_arg)) {
  sub("^--family=", "", family_arg[[1L]])
} else {
  "obsolete"
}
execute <- "--execute" %in% args

suppressPackageStartupMessages({
  library(targets)
  library(digest)
})

family_patterns <- list(
  primary = c(
    "^model_fit_new_",
    "^model_fit_remedied_",
    "^model_fit_second_remedied_",
    "^model_fit_third_remedied_",
    "^model_fit_final_",
    "^models_new_",
    "^primary_model_average_",
    "^primary_heatmap_plot",
    "^primary_psis_loo_exact_jobs$",
    "^primary_psis_loo_exact_.*_job$",
    "^primary_psis_loo_exact_.*_refit$"
  ),
  mundlak = c(
    "^sensitivity_mundlak_.*_fit$",
    "^sensitivity_mundlak_.*_fit_file$",
    "^sensitivity_mundlak_.*_fits$",
    "^sensitivity_mundlak_.*_average$",
    "^sensitivity_mundlak_.*_loo_jobs$",
    "^sensitivity_mundlak_full_integrated_loo($|_)",
    "^sensitivity_mundlak_full_exact_loo($|_)",
    "^sensitivity_mundlak_full_final_loo_audit$",
    "^sensitivity_mundlak_full_prediction($|_)",
    "^sensitivity_mundlak_full_model_average($|_)",
    "^sensitivity_mundlak_full_reporting($|_)"
  ),
  within = c(
    "^sensitivity_within_.*_fit$",
    "^sensitivity_within_.*_fit_file$",
    "^sensitivity_within_.*_fits$",
    "^sensitivity_within_.*_average$",
    "^sensitivity_within_.*_loo_jobs$",
    "^sensitivity_within_full_integrated_loo($|_)",
    "^sensitivity_within_full_exact_loo($|_)",
    "^sensitivity_within_full_final_loo_audit$",
    "^sensitivity_within_full_prediction($|_)",
    "^sensitivity_within_full_model_average($|_)",
    "^sensitivity_within_full_reporting($|_)"
  ),
  exposure = c(
    "^sensitivity_achieved_.*_fit$",
    "^sensitivity_achieved_.*_fit_file$",
    "^sensitivity_achieved_.*_fits$",
    "^sensitivity_achieved_.*_average$",
    "^sensitivity_achieved_full_integrated_loo($|_)",
    "^sensitivity_achieved_full_exact_loo($|_)",
    "^sensitivity_achieved_full_final_loo_audit$",
    "^sensitivity_achieved_full_prediction($|_)",
    "^sensitivity_achieved_full_model_average($|_)",
    "^sensitivity_achieved_full_reporting($|_)",
    "^sensitivity_assignment_.*_fit$",
    "^sensitivity_assignment_.*_fit_file$",
    "^sensitivity_assignment_.*_fits$",
    "^sensitivity_assignment_.*_average$",
    "^sensitivity_assignment_full_integrated_loo($|_)",
    "^sensitivity_assignment_full_exact_loo($|_)",
    "^sensitivity_assignment_full_final_loo_audit$",
    "^sensitivity_assignment_full_prediction($|_)",
    "^sensitivity_assignment_full_model_average($|_)",
    "^sensitivity_assignment_full_reporting($|_)",
    "^sensitivity_exposure_.*_loo_jobs$",
    "^sensitivity_exposure_full_comparison($|_)"
  ),
  eiv = c(
    "^sensitivity_eiv_.*_fit$",
    "^sensitivity_eiv_.*_fit_file$",
    "^sensitivity_eiv_.*_fits$",
    "^sensitivity_eiv_.*_average$",
    "^sensitivity_eiv_.*_loo_jobs$",
    "^sensitivity_eiv_stan_model$",
    "^sensitivity_eiv_full_integrated_loo($|_)",
    "^sensitivity_eiv_full_exact_loo($|_)",
    "^sensitivity_eiv_full_final_loo_audit$",
    "^sensitivity_eiv_full_prediction($|_)",
    "^sensitivity_eiv_full_model_average($|_)",
    "^sensitivity_eiv_full_reporting($|_)"
  ),
  rescor = c(
    "^sensitivity_rescor_.*_fit$",
    "^sensitivity_rescor_.*_fit_file$",
    "^sensitivity_rescor_.*_fits$",
    "^sensitivity_rescor_.*_average$",
    "^sensitivity_rescor_.*_loo_jobs$",
    "^sensitivity_rescor(_full)?_stan_model$",
    "^sensitivity_rescor_full_observation_loo_stan_model$",
    "^sensitivity_rescor_full_observation_loo($|_)",
    "^sensitivity_rescor_full_observation_exact($|_)",
    "^sensitivity_rescor_full_final_observation_loo_audit$",
    "^sensitivity_rescor_full_integrated_loo($|_)",
    "^sensitivity_rescor_full_exact_loo($|_)",
    "^sensitivity_rescor_full_final_loo_audit$",
    "^sensitivity_rescor_full_prediction($|_)",
    "^sensitivity_rescor_full_model_average($|_)",
    "^sensitivity_rescor_full_reporting($|_)"
  ),
  strategy = c(
    "^sensitivity_strategy_pilot_fit($|_)",
    "^sensitivity_strategy_pilot_fits$",
    "^sensitivity_strategy_remedy_fit($|_)",
    "^sensitivity_strategy_remedy_refits$",
    "^sensitivity_strategy_remedy_ensemble_fits$",
    "^sensitivity_strategy_pilot_average$",
    "^sensitivity_strategy_pilot_.*loo($|_)",
    "^sensitivity_strategy_full_.*fit($|_)",
    "^sensitivity_strategy_full_integrated_loo($|_)",
    "^sensitivity_strategy_full_exact_loo($|_)",
    "^sensitivity_strategy_full_final_loo_audit$",
    "^sensitivity_strategy_full_prediction($|_)",
    "^sensitivity_strategy_full_model_average($|_)",
    "^sensitivity_strategy_full_reporting($|_)"
  ),
  strategy_within = c(
    "^sensitivity_strategy_within_pilot_fit($|_)",
    "^sensitivity_strategy_within_remedy_fit($|_)",
    "^sensitivity_strategy_within_remedy_ensemble_fits$",
    "^sensitivity_strategy_within_pilot_average$",
    "^sensitivity_strategy_within_pilot_.*loo($|_)",
    "^sensitivity_strategy_within_full_.*fit($|_)",
    "^sensitivity_strategy_within_full_integrated_loo($|_)",
    "^sensitivity_strategy_within_full_exact_loo($|_)",
    "^sensitivity_strategy_within_full_final_loo_audit$",
    "^sensitivity_strategy_within_full_prediction($|_)",
    "^sensitivity_strategy_within_full_model_average($|_)",
    "^sensitivity_strategy_within_full_reporting($|_)"
  )
)

family_user_directories <- list(
  mundlak = c(
    file.path(
      "_targets", "user", "sensitivity_fit_files", "mundlak"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "mundlak_exact"
    )
  ),
  within = c(
    file.path(
      "_targets", "user", "sensitivity_fit_files", "within"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "within_exact"
    )
  ),
  exposure = c(
    file.path(
      "_targets", "user", "sensitivity_fit_files", "achieved"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "achieved_exact"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "assignment"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "assignment_exact"
    )
  ),
  eiv = c(
    file.path(
      "_targets", "user", "sensitivity_fit_files", "eiv"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "eiv_exact"
    ),
    file.path(
      "_targets", "user", "sensitivity_cmdstan", "eiv"
    ),
    file.path(
      "_targets", "user", "sensitivity_cmdstan", "eiv_full"
    )
  ),
  rescor = c(
    file.path(
      "_targets", "user", "sensitivity_fit_files", "rescor"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "rescor_exact"
    ),
    file.path(
      "_targets", "user", "sensitivity_cmdstan", "rescor"
    ),
    file.path(
      "_targets", "user", "sensitivity_cmdstan", "rescor_full"
    )
  ),
  strategy = c(
    file.path(
      "_targets", "user", "sensitivity_fit_files", "strategy"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "strategy_exact"
    )
  ),
  strategy_within = c(
    file.path(
      "_targets", "user", "sensitivity_fit_files", "strategy_within"
    ),
    file.path(
      "_targets", "user", "sensitivity_fit_files", "strategy_within_exact"
    )
  )
)

directory_size <- function(path) {
  files <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )
  if (!length(files)) return(0)
  sum(file.info(files)$size, na.rm = TRUE)
}

validate_user_directory <- function(path) {
  if (!dir.exists(path)) return(NA_character_)
  user_root <- normalizePath(
    file.path("_targets", "user"),
    winslash = "/",
    mustWork = TRUE
  )
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root_prefix <- paste0(tolower(user_root), "/")
  if (
    identical(tolower(resolved), tolower(user_root)) ||
      !startsWith(tolower(resolved), root_prefix)
  ) {
    stop(
      "Refusing to prune a directory outside _targets/user: ",
      resolved,
      call. = FALSE
    )
  }
  resolved
}

verify_bundle <- function(family) {
  manifest_path <- file.path(
    "results", "verification", family, "manifest.csv"
  )
  if (!file.exists(manifest_path)) {
    stop(
      "Refusing to prune: missing verification manifest at ",
      manifest_path,
      call. = FALSE
    )
  }

  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  paths <- file.path(dirname(manifest_path), manifest$file)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Refusing to prune: bundle files are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  observed <- vapply(
    paths,
    digest::digest,
    character(1),
    file = TRUE,
    algo = "sha256"
  )
  bad <- paths[observed != manifest$sha256]
  if (length(bad)) {
    stop(
      "Refusing to prune: checksum mismatch for: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

if (identical(family, "obsolete")) {
  candidates <- targets::tar_prune_list()
  cat("Obsolete cache objects eligible for targets::tar_prune():\n")
  print(candidates)
  if (execute) {
    targets::tar_prune()
  } else {
    cat("\nDry run only. Re-run with --execute to delete these objects.\n")
  }
  quit(save = "no", status = 0L)
}

if (!family %in% names(family_patterns)) {
  stop(
    "Unknown family '", family, "'. Choose obsolete, ",
    paste(names(family_patterns), collapse = ", "),
    ".",
    call. = FALSE
  )
}

verification_families <- if (identical(family, "exposure")) {
  c("achieved", "assignment", "exposure_comparison")
} else {
  family
}
invisible(lapply(verification_families, verify_bundle))
meta <- targets::tar_meta(fields = c(name, bytes, parent))
matches_pattern <- function(name) {
  any(vapply(family_patterns[[family]], grepl, logical(1), x = name))
}
candidates <- meta$name[vapply(meta$name, matches_pattern, logical(1))]
expanded_candidates <- unique(c(
  candidates,
  meta$name[!is.na(meta$parent) & meta$parent %in% candidates]
))
candidate_meta <- meta[meta$name %in% expanded_candidates, , drop = FALSE]
candidate_meta$selected_directly <- candidate_meta$name %in% candidates
candidate_meta <- candidate_meta[order(candidate_meta$bytes, decreasing = TRUE), ]

cat("Verified the ", family, " reporting bundle.\n", sep = "")
cat("Heavy cached targets eligible for deletion:\n")
print(candidate_meta, n = Inf)
cat(
  "Total selected size including dynamic branches: ",
  round(sum(candidate_meta$bytes, na.rm = TRUE) / 1024^3, 2),
  " GiB\n",
  sep = ""
)

user_directories <- family_user_directories[[family]]
if (is.null(user_directories)) user_directories <- character()
user_directories <- vapply(
  user_directories,
  validate_user_directory,
  character(1)
)
user_directories <- stats::na.omit(user_directories)
if (length(user_directories)) {
  directory_meta <- data.frame(
    directory = unname(user_directories),
    bytes = vapply(user_directories, directory_size, numeric(1)),
    stringsAsFactors = FALSE
  )
  cat("Persistent family directories eligible for deletion:\n")
  print(directory_meta, row.names = FALSE)
  cat(
    "Total persistent-directory size: ",
    round(sum(directory_meta$bytes) / 1024^3, 2),
    " GiB\n",
    sep = ""
  )
}

if (!execute) {
  cat("\nDry run only. Re-run with --execute to delete these objects.\n")
  quit(save = "no", status = 0L)
}

if (length(candidates)) {
  rlang::inject(targets::tar_delete(!!candidates))
}
for (directory in user_directories) {
  unlink(directory, recursive = TRUE, force = TRUE)
  if (dir.exists(directory)) {
    stop("Could not remove persistent directory: ", directory,
         call. = FALSE)
  }
}
cat("Pruning complete. Versioned verification artifacts were not removed.\n")
