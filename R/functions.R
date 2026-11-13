# Functions for targets pipeline

# Quarto and sampler setup ----
# These helpers configure local execution and resolve the Quarto executable
# before a render. The renderer checks its declared report inputs and restores
# the caller's environment on exit.
set_rstan_options <- function() {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)
}

find_quarto_cli <- function() {
  candidates <- c(
    Sys.getenv("QUARTO_PATH", unset = ""),
    unname(Sys.which("quarto")),
    Sys.getenv("RSTUDIO_QUARTO", unset = "")
  )

  if (.Platform$OS.type == "windows") {
    candidates <- c(
      candidates,
      file.path(
        Sys.getenv("ProgramFiles", unset = ""),
        "RStudio", "resources", "app", "bin", "quarto", "bin",
        "quarto.exe"
      ),
      file.path(
        Sys.getenv("ProgramFiles", unset = ""),
        "Quarto", "bin", "quarto.exe"
      )
    )
  }

  candidates <- unique(candidates[nzchar(candidates)])
  valid <- candidates[file.exists(candidates)]

  if (!length(valid)) {
    stop(
      paste(
        "Quarto CLI not found. Install Quarto or use an RStudio release",
        "that bundles it before building a document target."
      ),
      call. = FALSE
    )
  }

  normalizePath(valid[[1L]], winslash = "/", mustWork = TRUE)
}

render_quarto_document <- function(input, output, dependencies = character()) {
  missing_dependencies <- dependencies[!file.exists(dependencies)]
  if (length(missing_dependencies)) {
    stop(
      "Required report inputs are missing: ",
      paste(missing_dependencies, collapse = ", "),
      call. = FALSE
    )
  }

  quarto <- find_quarto_cli()
  old_quarto_r <- Sys.getenv("QUARTO_R", unset = NA_character_)
  old_local_app_data <- Sys.getenv("LOCALAPPDATA", unset = NA_character_)
  on.exit(
    {
      if (is.na(old_quarto_r)) {
        Sys.unsetenv("QUARTO_R")
      } else {
        Sys.setenv(QUARTO_R = old_quarto_r)
      }
      if (.Platform$OS.type == "windows") {
        if (is.na(old_local_app_data)) {
          Sys.unsetenv("LOCALAPPDATA")
        } else {
          Sys.setenv(LOCALAPPDATA = old_local_app_data)
        }
      }
    },
    add = TRUE
  )

  Sys.setenv(QUARTO_R = normalizePath(R.home("bin"), winslash = "/"))
  if (.Platform$OS.type == "windows") {
    quarto_cache <- file.path("_targets", "quarto-cache")
    dir.create(quarto_cache, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(
      LOCALAPPDATA = normalizePath(
        quarto_cache,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  status <- system2(
    quarto,
    args = c("render", shQuote(input)),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L) || !file.exists(output)) {
    stop("Quarto failed to render: ", input, call. = FALSE)
  }

  normalizePath(output, winslash = "/", mustWork = TRUE)
}

# Preliminary re-analysis with Nunes et al. (https://doi.org/10.1002/jcsm.12922) lean soft tissue/muscle data ----

# Historical re-analysis and prior calibration ----

# Combine the arm characteristics and lean-tissue outcome extracts from Nunes
# et al. Estimate pooled pre-post correlations and prepare standardized
# changes for the historical re-analysis.
read_prepare_nunes_data <- function(file1, file2) {

  # read in and prepare study characteristics csv for combining with effect sizes
  study_characteristics <- read_csv(file1) |>
    janitor::clean_names() |>
    select(id, ref, sex, age, duration, resistance_exercise, int_group, int_protein_intake, con_protein_intake) |>

    # pivot to arm based format
    pivot_longer(
      8:9,
      names_to = "condition",
      values_to = "protein_intake"
    ) |>

    # convert to NA missing values (i.e, 10000s in original dataset) and relabel condition
    mutate(
      protein_intake = if_else(protein_intake == 10000, NA, protein_intake),
      condition = case_when(
        condition == "int_protein_intake" ~ "int",
        condition == "con_protein_intake" ~ "con"
      )
    ) |>
    rename(intervention = "int_group")

  # read in and prepare lean soft tissue/muscle outcome data to calculate effect sizes
  data <- read_csv(file2) |>
    janitor::clean_names() |>
    mutate(
      # Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
      int_ri = (int_pre_sd^2 + int_post_sd^2 - int_delta_sd^2)/(2 * int_pre_sd * int_post_sd),
      con_ri = (con_pre_sd^2 + con_post_sd^2 - con_delta_sd^2)/(2 * con_pre_sd * con_post_sd)
    ) |>
    mutate(
      # Remove values outside the range of -1 to +1 as they are likely due to misreporting, miscalculations, or unreported dropout in original studies
      int_ri = if_else(between(int_ri,-1,1) == FALSE, NA, int_ri),
      con_ri = if_else(between(con_ri,-1,1) == FALSE, NA, con_ri)
    ) |>

    # calculate pooled SD to use for denominator in effect size calculations
    mutate(pre_sd_pool = sqrt(((int_n - 1) * int_pre_sd ^ 2 +
                                 (con_n - 1) * con_pre_sd ^ 2) /  (int_n + con_n - 2))) |>
    # add effect code
    rowid_to_column("effect") |>
    # add arm code
    unite("arm", c(id,intervention), sep = "_", remove = FALSE) |>
    mutate(arm = dense_rank(arm))

  # Calculate pre-post correlations for SMCR calculations

  data <- escalc(measure = "ZCOR", ri = int_ri, ni = int_n, data = data)

  meta_int_ri <- rma.mv(yi, V=vi, data=data,
                        random = list(~ 1 | id, ~1 | arm, ~1 | effect), method="REML", test="t",
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))

  rob_meta_int_ri <- robust(meta_int_ri, data$id)

  z2r_int <- psych::fisherz2r(rob_meta_int_ri$b[1])

  data$int_ri <- if_else(is.na(data$int_ri), z2r_int, data$int_ri)

  data_con <- data |>
    group_by(id) |>
    slice_head()

  data_con <- escalc(measure = "ZCOR", ri = con_ri, ni = con_n, data = data_con)

  meta_con_ri <- rma.mv(yi, V=vi, data=data_con,
                        random = list(~ 1 | id, ~1 | arm, ~1 | effect), method="REML", test="t",
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))

  rob_meta_con_ri <- robust(meta_con_ri, data_con$id)

  z2r_con <- psych::fisherz2r(rob_meta_con_ri$b[1])

  data$con_ri <- if_else(is.na(data$con_ri), z2r_con, data$con_ri)

  # Pivot longer so in arm based format

  data <- data |>
    select(-yi, -vi) |>
    pivot_longer(
      cols = contains(c("int_", "con_")),
      names_to = "what",
      values_to = "value"
    ) |>
    separate("what", into = c("condition", "what"), sep = "_", extra = "merge") |>
    pivot_wider(names_from = "what",
                values_from = "value")


  # Join with study characteristics
  data_joined <- left_join(data, study_characteristics, by = c("id", "condition", "intervention")) |>
    filter(!is.na(protein_intake)) |>
    group_by(id, condition, intervention) |>
    slice_head(n=1)

  # calculate for pre-post and delta separately
  data_joined_pre_post <- data_joined |>
    filter(!is.na(post_m))

  data_joined_pre_post <- escalc(
    measure = "SMCR",
    m1i = post_m,
    m2i = pre_m,
    sd1i = pre_sd_pool,
    # sd2i = pre_sd,
    ni = n,
    ri = ri,
    data = data_joined_pre_post
  )

  data_joined_delta <- data_joined |>
    filter(is.na(post_m)) |>
    mutate(ref = 0)

  data_joined_delta <- escalc(
    measure = "SMCR",
    m1i = delta_m,
    m2i = ref,
    sd1i = pre_sd_pool,
    # sd2i = pre_sd,
    ni = n,
    ri = ri,
    data = data_joined_delta
  )

  data_joined <- bind_rows(data_joined_pre_post, data_joined_delta)

  data_joined <- data_joined |>

    # add study weights/sizes
    mutate(
      wi = 1/sqrt(vi),
      size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE))) |>

    # centre  predictors
    mutate(
      protein_intake = protein_intake - 1.12,
      duration_centre = duration - 12,
      age_centre = age - 25,
      resistance_exercise_code = case_when(
        resistance_exercise == "YES" ~ 0,
        resistance_exercise == "NO" ~ 1,
      )
    ) |>

  # re-add effect and arm codes now all in long form

  # add effect code
  select(-effect) |>
  rowid_to_column("effect") |>
    # add arm code
    unite("arm", c(id,intervention, condition), sep = "_", remove = FALSE) |>
    mutate(arm = dense_rank(arm))

}

# Pool standardized changes in lean tissue and muscle from the historical
# resistance-exercise studies. Use the estimates and standard errors to set
# priors for the primary intercept and exercise coefficient.

# Note these are for hypertrophy (lean mass/muscle mass/size) only as Nunes reanalysis was only with these outcomes
get_steele_priors <- function(data) {

  # uses data from Steele et al. (https://doi.org/10.1080/02640414.2023.2286748) for lean soft tissue/muscle data (i.e., "hypertrophy" outcomes)

  steele_data <- data |>
    mutate(
      # calculate pre-post SDs from SEs
      RT_pre_sd = if_else(is.na(RT_pre_se), RT_pre_sd, RT_pre_se * sqrt(RT_n)),
      CON_pre_sd = if_else(is.na(CON_pre_se), CON_pre_sd, CON_pre_se * sqrt(CON_n)),
      RT_post_sd = if_else(is.na(RT_post_se), RT_post_sd, RT_post_se * sqrt(RT_n)),
      CON_post_sd = if_else(is.na(CON_post_se), CON_post_sd, CON_post_se * sqrt(CON_n)),

      # convert p to t (Change scores)
      RT_delta_t_value = if_else(is.na(RT_delta_t_value), qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE), RT_delta_t_value),
      CON_delta_t_value = if_else(is.na(CON_delta_t_value), qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE), CON_delta_t_value),

      # convert t to SE (Change scores)
      RT_delta_se = if_else(is.na(RT_delta_se), if_else(is.na(RT_delta_m), (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value), RT_delta_se),
      CON_delta_se = if_else(is.na(CON_delta_se), if_else(is.na(CON_delta_m), (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value), CON_delta_se),

      # make positive
      RT_delta_se = if_else(RT_delta_se < 0, RT_delta_se * -1, RT_delta_se),
      CON_delta_se = if_else(CON_delta_se < 0, CON_delta_se * -1, CON_delta_se),

      # convert CI to SE (Change scores)
      RT_delta_se = if_else(is.na(RT_delta_se), (RT_delta_CI_upper - RT_delta_CI_lower)/3.92, RT_delta_se),
      CON_delta_se = if_else(is.na(CON_delta_se), (CON_delta_CI_upper - CON_delta_CI_lower)/3.92, CON_delta_se),

      # convert SE to SD (Change scores)
      RT_delta_sd = if_else(is.na(RT_delta_sd), RT_delta_se * sqrt(RT_n), RT_delta_sd),
      CON_delta_sd = if_else(is.na(CON_delta_sd), CON_delta_se * sqrt(CON_n), CON_delta_sd),

      # calculate pre-post correlation coefficient for those with pre, post, and delta SDs
      RT_ri = (RT_pre_sd^2 + RT_post_sd^2 - RT_delta_sd^2)/(2 * RT_pre_sd * RT_post_sd),
      CON_ri = (CON_pre_sd^2 + CON_post_sd^2 - CON_delta_sd^2)/(2 * CON_pre_sd * CON_post_sd)
    ) |>
    mutate(
      # Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
      RT_ri = if_else(between(RT_ri,-1,1) == FALSE, NA, RT_ri),
      CON_ri = if_else(between(CON_ri,-1,1) == FALSE, NA, CON_ri)
    ) |>
    mutate(pre_sd_pool = sqrt(((RT_n - 1) * RT_pre_sd ^ 2 +
                                 (CON_n - 1) * CON_pre_sd ^ 2) /  (RT_n + CON_n - 2)))

  # Calculate pre-post correlations for SMCR calculations

  steele_data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = steele_data)

  meta_RT_ri <- rma.mv(yi, V=vi, data=steele_data,
                        random = list(~ 1 | study, ~1 | arm, ~1 | es), method="REML", test="t",
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))

  rob_meta_RT_ri <- robust(meta_RT_ri, steele_data$study)

  z2r_RT <- psych::fisherz2r(rob_meta_RT_ri$b[1])

  steele_data$RT_ri <- if_else(is.na(steele_data$RT_ri), z2r_RT, steele_data$RT_ri)

  steele_data_CON <- steele_data |>
    select(study, arm, es, CON_ri, CON_n) |>
    distinct()

  steele_data_CON <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = steele_data_CON)

  meta_CON_ri <- rma.mv(yi, V=vi, data=steele_data_CON,
                        random = list(~ 1 | study, ~1 | arm, ~1 | es), method="REML", test="t",
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))

  rob_meta_CON_ri <- robust(meta_CON_ri, steele_data_CON$study)

  z2r_CON <- psych::fisherz2r(rob_meta_CON_ri$b[1])

  steele_data$CON_ri <- if_else(is.na(steele_data$CON_ri), z2r_CON, steele_data$CON_ri)


  # calculate for pre-post and delta separately

  # Pivot longer so in arm based format

  steele_data <- steele_data |>
    select(-yi, -vi, -RT_only) |>
    pivot_longer(
      cols = contains(c("RT_", "CON_")),
      names_to = "what",
      values_to = "value"
    ) |>
    separate("what", into = c("condition", "what"), sep = "_", extra = "merge") |>
    pivot_wider(names_from = "what",
                values_from = "value")

  # calculate for pre-post and delta separately
  steele_data_pre_post <- steele_data |>
    filter(!is.na(post_m))

  steele_data_pre_post <- escalc(
    measure = "SMCR",
    m1i = post_m,
    m2i = pre_m,
    sd1i = pre_sd_pool,
    # sd2i = pre_sd,
    ni = n,
    ri = ri,
    data = steele_data_pre_post
  )

  steele_data_joined_delta <- steele_data |>
    filter(is.na(post_m)) |>
    mutate(ref = 0)

  data_joined_delta <- escalc(
    measure = "SMCR",
    m1i = delta_m,
    m2i = ref,
    sd1i = pre_sd_pool,
    # sd2i = pre_sd,
    ni = n,
    ri = ri,
    data = steele_data_joined_delta
  )

  steele_data_joined <- bind_rows(steele_data_pre_post, steele_data_joined_delta) |>
    filter(outcome == "hypertrophy") |>
    unite("arm", c(outcome,study,arm,condition), sep = "_", remove = FALSE) |>
    mutate(arm = dense_rank(arm)) |>
    rowid_to_column("effect") |>
    mutate(condition_code = case_when(
      condition == "RT" ~ 0,
      condition == "CON" ~ 1,
    ))

  # frequentist meta-analysis estimates
  steele_meta <- rma.mv(yi, vi,
                        random = list(~ 1 | study, ~1 | arm, ~1 | effect),
                        mods = ~ condition_code,
                        data = steele_data_joined,
                        method="REML", test="t",
                        # control=list(optimizer="optim", optmethod="Nelder-Mead")
  )

  rob_steele_meta <- robust(steele_meta, cluster = steele_data_joined$study, clubSandwich = TRUE)

  # set priors based on estimates (assumes habitual protein intake to get effects of RT and CON at that level)

  priors <-
    c(
      # # We set priors based on the overall estimates for both control and RT conditions in an arm based model using Steele et al. data
      set_prior(paste("student_t(3,", rob_steele_meta$b[1],",", rob_steele_meta$se[1],")"), class = "b", coef = "Intercept"),
      set_prior(paste("student_t(3,", rob_steele_meta$b[2],",", rob_steele_meta$se[2],")"), class = "b", coef = "resistance_exercise_code"),
      # All other b parameters are kept as uninformative but weakly regularising
      set_prior("student_t(3,0,1)", class = "b")
      # All other priors for variance parameters are kept as default weakly regularising
    )

  priors

}

# Determining baseline standard deviations to provide reference values for smallest effect size of interest of 0.1

# Reference values for interpreting effects ----

# Put the outcome-specific reference standard deviations and regions of
# practical equivalence on a consistent scale for the reporting functions.
estimate_lean_mass_sd <- function(data) {

  benito_data <-  data |>
    separate(study_arm, into = c("study", "arm"), sep = "_")|>
    group_by(study) |>
    mutate(
      study = cur_group_id()
    ) |>
    group_by(study, arm) |>
    mutate(
      arm = cur_group_id()
    ) |>
    mutate(
      duration_centre = duration/12
    ) |>
    filter(!is.na(duration)) |>
    rowid_to_column("effect")

  benito_data <- escalc(
    measure = "SDLN",
    sdi = pre_sd,
    ni = pre_n,
    data = benito_data
  )

  benito_meta_sd <- rma.mv(yi, vi,
                           random = list(~ 1 | study, ~ 1 | arm, ~ 1 | effect),
                           mods = ~ 0 + outcome,
                           data = benito_data,
                           method="REML", test="t")

  rob_benito_meta_sd <- robust(benito_meta_sd, benito_data$study)


  benito_sds <- broom::tidy(rob_benito_meta_sd) |>
    mutate(across(c(estimate, std.error), exp))

  return(benito_sds)
}

estimate_muscle_thickness_sd <- function(data) {

  buckner_data <- data |>
    janitor::clean_names() |>
    filter(retracted == "NO") |>
    mutate(centred_pre_mean = pre_mean - 5)

  buckner_data <- escalc(
    measure = "SDLN",
    sdi = pre_sd,
    ni = n,
    data = buckner_data
  )

  buckner_meta_sd <- rma.mv(yi, vi,
                            random = list(~ 1 | study, ~ 1 | group, ~ 1 | effect),
                            mods = ~ log(pre_mean),
                           data = buckner_data,
                           method="REML", test="t")

  rob_buckner_meta_sd <- robust(buckner_meta_sd, buckner_data$study)


  reference_mean <- median(buckner_data$pre_mean)
  buckner_sd <- as_tibble(predict(
    rob_buckner_meta_sd,
    newmods = log(reference_mean)
  )) |>
    mutate(
      across(where(is.double), exp),
      reference_mean = reference_mean
    )


  return(buckner_sd)

}

prepare_open_powerlifting_strength_data <- function(data) {

  data |>
    filter(Equipment == "Raw" & Tested == "Yes" & Division == "Open") |>
    mutate(Best3SquatKg = if_else(Best3SquatKg > 0, Best3SquatKg, NA), # We also remove any failed lifts as it impacts the totals
           Best3BenchKg = if_else(Best3BenchKg > 0, Best3BenchKg, NA),
           Best3DeadliftKg = if_else(Best3DeadliftKg > 0, Best3DeadliftKg, NA)) |>
    group_by(Name) |>
    arrange(Date) |>
    slice_head(n=1) |>
    ungroup()

}

read_open_powerlifting_strength_snapshot <- function(path) {
  snapshot <- readr::read_csv(
    path,
    col_types = readr::cols(
      snapshot_date = readr::col_date(),
      accessed_date = readr::col_date(),
      source_url = readr::col_character(),
      source_file = readr::col_character(),
      source_sha256 = readr::col_character(),
      source_rows = readr::col_double(),
      selected_lifters = readr::col_double(),
      squat_n = readr::col_double(),
      bench_n = readr::col_double(),
      deadlift_n = readr::col_double(),
      squat_sd = readr::col_double(),
      bench_sd = readr::col_double(),
      deadlift_sd = readr::col_double()
    ),
    show_col_types = FALSE
  )

  required_columns <- c(
    "snapshot_date", "accessed_date", "source_url", "source_file",
    "source_sha256", "source_rows", "selected_lifters", "squat_n",
    "bench_n", "deadlift_n", "squat_sd", "bench_sd", "deadlift_sd"
  )

  if (
    nrow(snapshot) != 1L ||
      !identical(names(snapshot), required_columns) ||
      !grepl("^[A-F0-9]{64}$", snapshot$source_sha256) ||
      any(!is.finite(unlist(snapshot[c(
        "source_rows", "selected_lifters", "squat_n", "bench_n",
        "deadlift_n", "squat_sd", "bench_sd", "deadlift_sd"
      )], use.names = FALSE)))
  ) {
    stop("The pinned OpenPowerlifting snapshot is invalid.", call. = FALSE)
  }

  snapshot
}

# Return the reference scales and equivalence thresholds for whole-body lean
# mass, muscle thickness, strength and function. Read the pinned strength
# snapshot when one is supplied.
make_primary_reference_values <- function(
    protein_data,
    lean_mass_sd,
    muscle_thickness_sd,
    open_powerlifting_snapshot,
    rope_half_width = 0.10
) {
  protein_median <- unique(protein_data$protein_median)
  lean_soft_tissue_sd <- lean_mass_sd$estimate[
    lean_mass_sd$term == "outcomelmm"
  ]
  fat_free_mass_sd <- lean_mass_sd$estimate[
    lean_mass_sd$term == "outcomeffm"
  ]

  if (
    length(protein_median) != 1L ||
      length(lean_soft_tissue_sd) != 1L ||
      length(fat_free_mass_sd) != 1L ||
      nrow(muscle_thickness_sd) != 1L ||
      nrow(open_powerlifting_snapshot) != 1L ||
      !is.finite(rope_half_width) ||
      rope_half_width <= 0
  ) {
    stop("Could not assemble the manuscript reference values.", call. = FALSE)
  }

  list(
    habitual_protein_intake = protein_median,
    rope_half_width = rope_half_width,
    lean_soft_tissue_mass_change_kg =
      rope_half_width * lean_soft_tissue_sd,
    fat_free_mass_change_kg = rope_half_width * fat_free_mass_sd,
    muscle_thickness_reference_mean_cm =
      muscle_thickness_sd$reference_mean[[1L]],
    muscle_thickness_change_cm =
      rope_half_width * muscle_thickness_sd$pred[[1L]],
    squat_change_kg =
      rope_half_width * open_powerlifting_snapshot$squat_sd[[1L]],
    bench_change_kg =
      rope_half_width * open_powerlifting_snapshot$bench_sd[[1L]],
    deadlift_change_kg =
      rope_half_width * open_powerlifting_snapshot$deadlift_sd[[1L]],
    openpowerlifting_snapshot_date =
      open_powerlifting_snapshot$snapshot_date[[1L]],
    openpowerlifting_accessed_date =
      open_powerlifting_snapshot$accessed_date[[1L]],
    openpowerlifting_source_file =
      open_powerlifting_snapshot$source_file[[1L]],
    openpowerlifting_source_sha256 =
      open_powerlifting_snapshot$source_sha256[[1L]]
  )
}

download_open_powerlifting_strength_sd <- function(url) {

  archive <- tempfile(fileext = ".zip")
  extract_dir <- tempfile(pattern = "openpowerlifting-")
  dir.create(extract_dir)

  on.exit(
    unlink(c(archive, extract_dir), recursive = TRUE, force = TRUE),
    add = TRUE
  )

  utils::download.file(
    url = url,
    destfile = archive,
    mode = "wb",
    quiet = FALSE
  )

  archive_listing <- utils::unzip(archive, list = TRUE)
  csv_member <- archive_listing$Name[
    grepl(
      "^openpowerlifting-.*\\.csv$",
      basename(archive_listing$Name),
      ignore.case = TRUE
    )
  ]

  if (length(csv_member) != 1) {
    stop(
      "Expected exactly one OpenPowerlifting CSV in the downloaded archive; found ",
      length(csv_member),
      ".",
      call. = FALSE
    )
  }

  utils::unzip(
    zipfile = archive,
    files = csv_member,
    exdir = extract_dir,
    junkpaths = TRUE
  )

  extracted_csv <- file.path(extract_dir, basename(csv_member))
  required_columns <- c(
    "Name",
    "Equipment",
    "Tested",
    "Division",
    "Date",
    "Best3SquatKg",
    "Best3BenchKg",
    "Best3DeadliftKg"
  )

  source_data <- readr::read_csv(
    extracted_csv,
    col_select = tidyselect::all_of(required_columns),
    show_col_types = FALSE,
    progress = FALSE
  )

  selected_data <- prepare_open_powerlifting_strength_data(source_data)

  selected_data |>
    summarise(
      source_url = url,
      source_file = basename(csv_member),
      source_sha256 = toupper(
        digest::digest(
          file = archive,
          algo = "sha256",
          serialize = FALSE
        )
      ),
      source_rows = nrow(source_data),
      selected_lifters = n(),
      squat_n = sum(!is.na(Best3SquatKg)),
      bench_n = sum(!is.na(Best3BenchKg)),
      deadlift_n = sum(!is.na(Best3DeadliftKg)),
      squat_sd = sd(Best3SquatKg, na.rm = TRUE),
      bench_sd = sd(Best3BenchKg, na.rm = TRUE),
      deadlift_sd = sd(Best3DeadliftKg, na.rm = TRUE)
    )

}

# Function to fit all candidate models
fit_candidate_models <- function(formula, priors, data) {
  formula <- as.formula(paste0(formula))
  # Stored preliminary priors may come from an older brms release whose
  # `brmsprior` table included columns (for example `tag`) that are absent from
  # the current schema. Normalise before brms combines them with default priors.
  priors <- normalize_primary_brms_priors(priors)
  model <- brm(
    formula,
    data = data,
    prior = priors,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 6000, # just for pre-reg so it doesn't take so long to run
    # iter = 42000, # so as to have 40k iterations for BF calculations
    control = list(adapt_delta = 0.99),
    save_pars = save_pars(all = TRUE)
  )

  return(model)
}

# Functions to define, fit, and diagnose the primary candidate models ----

# Candidate formula and prior adjustments ----

# These helpers modify individual model terms and priors. The specification
# builders below record which adjustments apply to each candidate and each
# refit stage.
make_uncorrelated_study_formula <- function(formula) {
  updated_formula <- stringr::str_replace_all(
    formula,
    stringr::fixed(" | study_no)"),
    " || study_no)"
  )

  if (identical(updated_formula, formula)) {
    stop(
      "Could not find a correlated study-level term to make uncorrelated.",
      call. = FALSE
    )
  }

  updated_formula
}

add_regularised_variance_priors <- function(priors) {
  c(
    priors,
    brms::set_prior("student_t(3, 0, 1)", class = "sd"),
    brms::set_prior("student_t(3, 0, 1)", class = "sds")
  )
}

remove_global_primary_prior <- function(priors, class) {
  priors[
    !(
      priors$class == class &
        priors$coef == "" &
        priors$group == ""
    ),
  ]
}

add_primary_standardised_moderators <- function(data) {
  moderator_centres <- c(
    age_centre = mean(data$age_centre),
    duration_centre = mean(data$duration_centre)
  )
  moderator_scales <- c(
    age_centre = stats::sd(data$age_centre),
    duration_centre = stats::sd(data$duration_centre)
  )

  if (
    any(!is.finite(moderator_centres)) ||
      any(!is.finite(moderator_scales)) ||
      any(moderator_scales <= 0)
  ) {
    stop(
      "Age and duration could not be standardised for this outcome.",
      call. = FALSE
    )
  }

  data |>
    dplyr::mutate(
      age_outcome_z =
        (age_centre - moderator_centres[["age_centre"]]) /
          moderator_scales[["age_centre"]],
      duration_outcome_z =
        (duration_centre - moderator_centres[["duration_centre"]]) /
          moderator_scales[["duration_centre"]]
    ) |>
    structure(
      moderator_centres = moderator_centres,
      moderator_scales = moderator_scales
    )
}

add_primary_half_unit_protein <- function(data) {
  dplyr::mutate(
    data,
    protein_intake_05 = protein_intake / 0.5
  )
}

replace_primary_standardised_moderators <- function(formula) {
  formula |>
    stringr::str_replace_all(
      "\\bage_centre\\b",
      "age_outcome_z"
    ) |>
    stringr::str_replace_all(
      "\\bduration_centre\\b",
      "duration_outcome_z"
    )
}

replace_primary_half_unit_protein <- function(formula) {
  stringr::str_replace_all(
    formula,
    "\\bprotein_intake\\b",
    "protein_intake_05"
  )
}

replace_primary_tps_with_natural_spline <- function(formula) {
  updated_formula <- stringr::str_replace(
    formula,
    paste0(
      "s\\(\\s*protein_intake\\s*,\\s*",
      "bs\\s*=\\s*['\"]tp['\"]\\s*,\\s*",
      "k\\s*=\\s*4\\s*\\)"
    ),
    "splines::ns(protein_intake, df = 3)"
  )

  if (identical(updated_formula, formula)) {
    stop(
      "Could not find the simple TPS term to replace.",
      call. = FALSE
    )
  }

  updated_formula
}

add_primary_standardised_moderator_priors <- function(priors) {
  c(
    priors,
    brms::set_prior(
      "student_t(3, 0, 0.5)",
      class = "b",
      coef = "age_outcome_z"
    ),
    brms::set_prior(
      "student_t(3, 0, 0.5)",
      class = "b",
      coef = "duration_outcome_z"
    )
  )
}

add_primary_targeted_smooth_prior <- function(
    priors,
    scale = 0.5
) {
  priors <- remove_global_primary_prior(priors, "sds")

  c(
    priors,
    brms::set_prior(
      paste0("normal(0, ", scale, ")"),
      class = "sds"
    )
  )
}

normalize_primary_brms_priors <- function(priors) {
  template <- as.data.frame(
    brms::set_prior("student_t(3, 0, 1)", class = "b")
  )
  prior_table <- as.data.frame(priors)
  missing_columns <- setdiff(names(template), names(prior_table))

  for (column in missing_columns) {
    prior_table[[column]] <- if (column %in% c("lb", "ub")) {
      NA_character_
    } else if (identical(column, "source")) {
      "user"
    } else {
      ""
    }
  }

  prior_table <- prior_table[, names(template), drop = FALSE]
  structure(
    prior_table,
    class = c("brmsprior", "data.frame")
  )
}

add_primary_piecewise_slope_priors <- function(
    priors,
    coefficient_scale = 0.5,
    random_slope_scale = 0.5
) {
  priors <- normalize_primary_brms_priors(priors)
  coefficient_prior <- paste0(
    "student_t(3, 0, ", coefficient_scale, ")"
  )
  random_slope_prior <- paste0(
    "student_t(3, 0, ", random_slope_scale, ")"
  )
  c(
    priors,
    brms::set_prior(
      coefficient_prior,
      class = "b",
      coef = "lsplineprotein_intake0.51"
    ),
    brms::set_prior(
      coefficient_prior,
      class = "b",
      coef = "lsplineprotein_intake0.52"
    ),
    brms::set_prior(
      random_slope_prior,
      class = "sd",
      group = "study_no",
      coef = "lsplineprotein_intake0.51"
    ),
    brms::set_prior(
      random_slope_prior,
      class = "sd",
      group = "study_no",
      coef = "lsplineprotein_intake0.52"
    )
  )
}

add_primary_function_tail_priors <- function(
    priors,
    include_interactions = FALSE,
    coefficient_scale = 0.5,
    first_random_slope_scale = 0.5,
    second_random_slope_scale = 0.25
) {
  priors <- normalize_primary_brms_priors(priors)
  fixed_coefficients <- c(
    "lsplineprotein_intake0.51",
    "lsplineprotein_intake0.52"
  )
  interaction_coefficients <- c(
    "resistance_exercise_code:lsplineprotein_intake0.51",
    "resistance_exercise_code:lsplineprotein_intake0.52"
  )
  coefficients_to_replace <- c(
    fixed_coefficients,
    if (isTRUE(include_interactions)) interaction_coefficients
  )
  replace <-
    (priors$class == "b" & priors$coef %in% coefficients_to_replace) |
    (
      priors$class == "sd" &
        priors$group == "study_no" &
        priors$coef %in% fixed_coefficients
    )
  priors <- priors[!replace, , drop = FALSE]
  priors <- structure(priors, class = c("brmsprior", "data.frame"))

  coefficient_prior <- paste0(
    "student_t(3, 0, ", coefficient_scale, ")"
  )
  priors <- c(
    priors,
    brms::set_prior(
      coefficient_prior,
      class = "b",
      coef = fixed_coefficients[[1L]]
    ),
    brms::set_prior(
      coefficient_prior,
      class = "b",
      coef = fixed_coefficients[[2L]]
    ),
    brms::set_prior(
      paste0("student_t(3, 0, ", first_random_slope_scale, ")"),
      class = "sd",
      group = "study_no",
      coef = fixed_coefficients[[1L]]
    ),
    brms::set_prior(
      paste0("normal(0, ", second_random_slope_scale, ")"),
      class = "sd",
      group = "study_no",
      coef = fixed_coefficients[[2L]]
    )
  )

  if (isTRUE(include_interactions)) {
    priors <- c(
      priors,
      brms::set_prior(
        coefficient_prior,
        class = "b",
        coef = interaction_coefficients[[1L]]
      ),
      brms::set_prior(
        coefficient_prior,
        class = "b",
        coef = interaction_coefficients[[2L]]
      )
    )
  }

  priors
}

# Primary candidate registry ----

# Turn the model names and formulas into indexed specifications with the
# settings needed by fitting and diagnostic targets.
get_primary_model_fit_specifications <- function(
    model_names,
    model_formulas
) {
  if (length(model_names) != 26L || length(model_formulas) != 26L) {
    stop("Expected exactly 26 candidate models.", call. = FALSE)
  }

  numeric_rt_smooth <- stringr::str_detect(
    model_formulas,
    "by\\s*=\\s*resistance_exercise_code"
  )
  factor_rt_smooth <- stringr::str_detect(
    model_formulas,
    "by\\s*=\\s*resistance_exercise_factor"
  )
  numeric_study_smooth <- stringr::str_detect(
    model_formulas,
    "s\\(\\s*protein_intake\\s*,\\s*study_no\\s*,\\s*bs\\s*=\\s*['\"]fs['\"]"
  )
  factor_study_smooth <- stringr::str_detect(
    model_formulas,
    "s\\(\\s*protein_intake\\s*,\\s*study_factor\\s*,\\s*bs\\s*=\\s*['\"]fs['\"]"
  )

  if (
    any(numeric_rt_smooth) ||
      sum(factor_rt_smooth) != 4L ||
      any(numeric_study_smooth) ||
      sum(factor_study_smooth) != 4L
  ) {
    stop(
      "The primary formulas must contain four factor-by RT smooths and ",
      "four factor study smooths, with no numeric-by equivalents.",
      call. = FALSE
    )
  }

  control_groups <- tibble::tibble(
    outcome = c(
      rep("whole_body", 4),
      rep("hypertrophy", 3),
      rep("strength", 4),
      rep("function", 4)
    ),
    adapt_delta = c(
      0.99, 0.995, 0.995, 0.999,
      0.99, 0.995, 0.999,
      0.99, 0.995, 0.995, 0.999,
      0.99, 0.995, 0.995, 0.999
    ),
    max_treedepth = c(
      10L, 12L, 15L, 15L,
      10L, 12L, 15L,
      10L, 12L, 15L, 15L,
      10L, 12L, 15L, 15L
    ),
    model_index = list(
      c(1, 2, 3, 9, 15, 16, 21, 22),
      c(4, 5, 6, 7, 8, 10, 12, 13, 14, 18, 19, 20, 24, 25),
      11,
      c(17, 23, 26),
      c(1, 2, 3, 9, 10, 15, 16, 21, 22),
      c(4, 5, 6, 7, 8, 11, 12, 13, 14, 18, 19, 20, 24, 25, 26),
      c(17, 23),
      c(1, 2),
      c(6, 7, 8, 9, 10, 13, 14, 15, 16, 19, 22, 25),
      c(3, 5, 12, 18, 21, 24),
      c(4, 11, 17, 20, 23, 26),
      c(1, 16),
      c(2, 3, 4, 9, 11, 13, 14, 15, 21, 22, 25, 26),
      c(5, 6, 7, 8, 12, 18, 19, 20, 24),
      c(10, 17, 23)
    )
  ) |>
    tidyr::unnest_longer(model_index) |>
    dplyr::mutate(model_index = as.integer(model_index))

  model_table <- tibble::tibble(
    model_index = seq_along(model_names),
    model_name = model_names,
    preregistered_formula = model_formulas
  )

  specifications <- purrr::map_dfr(
    c("whole_body", "hypertrophy", "strength", "function"),
    ~ dplyr::mutate(model_table, outcome = .x, .before = 1)
  ) |>
    dplyr::left_join(
      control_groups,
      by = c("outcome", "model_index")
    ) |>
    dplyr::mutate(
      pair = paste(outcome, model_index, sep = "_"),
      uncorrelated_study_effects = pair %in% c(
        "whole_body_10",
        "strength_12", "strength_13", "strength_21", "strength_24",
        "function_10", "function_12", "function_24"
      ),
      regularised_variance_priors = pair %in% c(
        "whole_body_11", "whole_body_17", "whole_body_23", "whole_body_26",
        "hypertrophy_11", "hypertrophy_23",
        "strength_4", "strength_11", "strength_17", "strength_23", "strength_26",
        "function_11", "function_17", "function_23"
      ),
      full_schedule_trial = pair %in% c(
        "whole_body_3", "whole_body_16",
        "hypertrophy_17", "hypertrophy_20"
      ),
      formula = purrr::map2_chr(
        preregistered_formula,
        uncorrelated_study_effects,
        ~ if (.y) make_uncorrelated_study_formula(.x) else .x
      )
    ) |>
    dplyr::select(-pair)

  if (
    nrow(specifications) != 104L ||
      anyNA(specifications$adapt_delta) ||
      anyNA(specifications$max_treedepth) ||
      sum(specifications$uncorrelated_study_effects) != 8L ||
      sum(specifications$regularised_variance_priors) != 14L ||
      sum(specifications$full_schedule_trial) != 4L ||
      any(
        specifications$uncorrelated_study_effects &
          specifications$regularised_variance_priors
      )
  ) {
    stop("The primary model fitting registry failed validation.", call. = FALSE)
  }

  specification_list <- purrr::pmap(specifications, list)
  names(specification_list) <- paste0(
    specifications$outcome,
    "_model_",
    stringr::str_pad(specifications$model_index, 2, pad = "0")
  )

  specification_list
}

# Diagnostic remedies ----

# Build the first refit specifications from the primary registry. Subsequent
# builders carry forward model indices so refitted candidates can be overlaid
# in their original slots.
get_primary_model_remedy_specifications <- function(specifications) {
  specification_table <- dplyr::bind_rows(specifications) |>
    dplyr::mutate(
      pair = paste(outcome, model_index, sep = "_"),
      standardised_moderators = stringr::str_detect(
        formula,
        "\\bage_centre\\b|\\bduration_centre\\b"
      ),
      tps_model = stringr::str_detect(
        formula,
        "s\\(\\s*protein_intake\\s*,"
      ),
      natural_spline = outcome == "strength" & model_index == 4L,
      half_unit_protein = tps_model & !natural_spline,
      targeted_smooth_prior = tps_model & !natural_spline,
      newly_uncorrelated_study_effects = pair %in% c(
        "function_13", "function_22", "function_25"
      ),
      uncorrelated_study_effects =
        uncorrelated_study_effects |
          newly_uncorrelated_study_effects,
      requires_refit =
        standardised_moderators |
          tps_model |
          newly_uncorrelated_study_effects,
      formula = purrr::map2_chr(
        formula,
        standardised_moderators,
        ~ if (.y) {
          replace_primary_standardised_moderators(.x)
        } else {
          .x
        }
      ),
      formula = purrr::map2_chr(
        formula,
        half_unit_protein,
        ~ if (.y) {
          replace_primary_half_unit_protein(.x)
        } else {
          .x
        }
      ),
      formula = purrr::map2_chr(
        formula,
        natural_spline,
        ~ if (.y) {
          replace_primary_tps_with_natural_spline(.x)
        } else {
          .x
        }
      ),
      formula = purrr::map2_chr(
        formula,
        newly_uncorrelated_study_effects,
        ~ if (.y) make_uncorrelated_study_formula(.x) else .x
      ),
      adapt_delta = dplyr::case_when(
        natural_spline ~ 0.995,
        pair %in% c("function_13", "function_22") ~ 0.999,
        targeted_smooth_prior ~ 0.999,
        TRUE ~ adapt_delta
      ),
      max_treedepth = dplyr::case_when(
        natural_spline ~ 12L,
        pair %in% c("function_13", "function_22") ~ 15L,
        targeted_smooth_prior ~ 15L,
        TRUE ~ max_treedepth
      ),
      smooth_prior_scale = dplyr::if_else(
        targeted_smooth_prior,
        0.5,
        NA_real_
      )
    ) |>
    dplyr::filter(requires_refit) |>
    dplyr::select(-pair)

  outcome_counts <- specification_table |>
    dplyr::count(outcome) |>
    tibble::deframe()

  if (
    nrow(specification_table) != 69L ||
      !identical(
        outcome_counts[c(
          "whole_body", "hypertrophy", "strength", "function"
        )],
        c(
          whole_body = 17L,
          hypertrophy = 17L,
          strength = 17L,
          "function" = 18L
        )
      ) ||
      sum(specification_table$standardised_moderators) != 52L ||
      sum(specification_table$half_unit_protein) != 31L ||
      sum(specification_table$natural_spline) != 1L ||
      sum(specification_table$targeted_smooth_prior) != 31L ||
      sum(
        specification_table$newly_uncorrelated_study_effects
      ) != 3L ||
      any(
        specification_table$standardised_moderators &
          stringr::str_detect(
            specification_table$formula,
            "\\bage_centre\\b|\\bduration_centre\\b"
          )
      ) ||
      any(
        specification_table$half_unit_protein &
          stringr::str_detect(
            specification_table$formula,
            "\\bprotein_intake\\b"
          )
      )
  ) {
    stop(
      "The primary post-pilot remedy registry failed validation.",
      call. = FALSE
    )
  }

  specification_list <- purrr::pmap(specification_table, list)
  names(specification_list) <- paste0(
    specification_table$outcome,
    "_model_",
    stringr::str_pad(
      specification_table$model_index,
      2,
      pad = "0"
    )
  )

  specification_list
}

select_primary_model_remedy_specifications <- function(
    specifications,
    outcome
) {
  selected <- purrr::keep(
    specifications,
    ~ identical(.x$outcome, outcome)
  )

  expected <- if (identical(outcome, "function")) 18L else 17L

  if (length(selected) != expected) {
    stop(
      "Expected ", expected,
      " remedied specifications for outcome: ",
      outcome,
      call. = FALSE
    )
  }

  selected
}

get_primary_model_second_remedy_specifications <- function(
    primary_specifications,
    first_remedy_specifications
) {
  primary_table <- dplyr::bind_rows(primary_specifications) |>
    dplyr::mutate(remedy_stage = 0L)
  first_remedy_table <- dplyr::bind_rows(first_remedy_specifications) |>
    dplyr::mutate(remedy_stage = 1L)

  current_table <- dplyr::bind_rows(
    primary_table,
    first_remedy_table
  ) |>
    dplyr::arrange(outcome, model_index, remedy_stage) |>
    dplyr::group_by(outcome, model_index) |>
    dplyr::slice_tail(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      pair = paste(outcome, model_index, sep = "_"),
      standardised_moderators = dplyr::coalesce(
        standardised_moderators,
        FALSE
      ),
      half_unit_protein = dplyr::coalesce(
        half_unit_protein,
        FALSE
      ),
      natural_spline = dplyr::coalesce(natural_spline, FALSE),
      targeted_smooth_prior = dplyr::coalesce(
        targeted_smooth_prior,
        FALSE
      ),
      newly_uncorrelated_study_effects = dplyr::coalesce(
        newly_uncorrelated_study_effects,
        FALSE
      )
    )

  remedy_pairs <- c(
    "whole_body_13",
    "strength_1", "strength_2", "strength_9", "strength_16",
    "strength_19",
    "function_1", "function_13", "function_16", "function_21",
    "function_25"
  )

  specification_table <- current_table |>
    dplyr::filter(pair %in% remedy_pairs) |>
    dplyr::mutate(
      adapt_delta = dplyr::case_when(
        pair %in% c("whole_body_13", "function_25") ~ 0.999,
        TRUE ~ adapt_delta
      ),
      max_treedepth = dplyr::case_when(
        pair %in% c(
          "function_1", "function_16",
          "strength_1", "strength_2"
        ) ~ 12L,
        TRUE ~ 15L
      ),
      regularised_variance_priors = dplyr::if_else(
        pair == "function_13",
        TRUE,
        regularised_variance_priors
      ),
      targeted_piecewise_slope_priors = pair == "function_13",
      second_remedy = TRUE,
      requires_refit = TRUE,
      second_remedy_reason = dplyr::case_when(
        pair == "whole_body_13" ~ paste(
          "one divergence; first sampler-only escalation before changing",
          "the random-effect covariance structure"
        ),
        pair == "function_13" ~ paste(
          "persistent divergences after stringent sampler controls;",
          "regularise piecewise slopes and their study-level SDs"
        ),
        pair == "function_25" ~ paste(
          "divergences; first escalate sampler controls before targeted",
          "interaction/random-slope regularisation"
        ),
        TRUE ~ "maximum-treedepth saturation; depth-only escalation"
      )
    ) |>
    dplyr::select(-pair, -remedy_stage)

  specification_pairs <- paste(
    specification_table$outcome,
    specification_table$model_index,
    sep = "_"
  )
  expected_adapt_delta <- c(
    whole_body_13 = 0.999,
    strength_1 = 0.99,
    strength_2 = 0.99,
    strength_9 = 0.995,
    strength_16 = 0.995,
    strength_19 = 0.995,
    function_1 = 0.99,
    function_13 = 0.999,
    function_16 = 0.99,
    function_21 = 0.995,
    function_25 = 0.999
  )

  expected_max_treedepth <- c(
    whole_body_13 = 15L,
    strength_1 = 12L,
    strength_2 = 12L,
    strength_9 = 15L,
    strength_16 = 15L,
    strength_19 = 15L,
    function_1 = 12L,
    function_13 = 15L,
    function_16 = 12L,
    function_21 = 15L,
    function_25 = 15L
  )

  function_13 <- specification_table |>
    dplyr::filter(outcome == "function", model_index == 13L)

  if (
    nrow(specification_table) != 11L ||
      !setequal(specification_pairs, remedy_pairs) ||
      anyDuplicated(specification_pairs) ||
      !identical(
        unname(specification_table$adapt_delta),
        unname(expected_adapt_delta[specification_pairs])
      ) ||
      !identical(
        unname(specification_table$max_treedepth),
        unname(expected_max_treedepth[specification_pairs])
      ) ||
      nrow(function_13) != 1L ||
      !isTRUE(function_13$uncorrelated_study_effects) ||
      !isTRUE(function_13$standardised_moderators) ||
      !isTRUE(function_13$regularised_variance_priors) ||
      !isTRUE(function_13$targeted_piecewise_slope_priors) ||
      sum(specification_table$targeted_piecewise_slope_priors) != 1L ||
      !all(specification_table$second_remedy) ||
      !all(specification_table$requires_refit)
  ) {
    stop(
      "The second-stage primary remedy registry failed validation.",
      call. = FALSE
    )
  }

  specification_list <- purrr::pmap(specification_table, list)
  names(specification_list) <- paste0(
    specification_table$outcome,
    "_model_",
    stringr::str_pad(
      specification_table$model_index,
      2,
      pad = "0"
    )
  )

  specification_list
}

select_primary_model_second_remedy_specifications <- function(
    specifications,
    outcome
) {
  selected <- purrr::keep(
    specifications,
    ~ identical(.x$outcome, outcome)
  )
  expected <- switch(
    outcome,
    whole_body = 1L,
    strength = 5L,
    "function" = 5L,
    0L
  )

  if (length(selected) != expected) {
    stop(
      "Expected ", expected,
      " second-stage remedies for outcome: ", outcome,
      call. = FALSE
    )
  }

  selected
}

get_primary_model_third_remedy_specifications <- function(
    second_remedy_specifications
) {
  specification_table <- dplyr::bind_rows(
    second_remedy_specifications
  ) |>
    dplyr::filter(
      outcome == "function",
      model_index %in% c(13L, 25L)
    ) |>
    dplyr::mutate(
      regularised_variance_priors = TRUE,
      targeted_piecewise_slope_priors = TRUE,
      targeted_function_tail_priors = TRUE,
      targeted_piecewise_interaction_priors = model_index == 25L,
      third_remedy = TRUE,
      requires_refit = TRUE,
      third_remedy_reason = dplyr::case_when(
        model_index == 13L ~ paste(
          "one residual divergence localised to the weakly identified",
          "second-segment study-slope SD; strengthen that SD prior"
        ),
        model_index == 25L ~ paste(
          "one residual divergence localised chiefly to the second",
          "piecewise coefficient and study-slope SD; apply targeted",
          "tail-piece and interaction regularisation"
        )
      )
    )

  expected_indices <- c(13L, 25L)
  if (
    nrow(specification_table) != 2L ||
      !identical(sort(specification_table$model_index), expected_indices) ||
      any(specification_table$outcome != "function") ||
      any(specification_table$adapt_delta != 0.999) ||
      any(specification_table$max_treedepth != 15L) ||
      !all(specification_table$uncorrelated_study_effects) ||
      !all(specification_table$standardised_moderators) ||
      !all(specification_table$regularised_variance_priors) ||
      !all(specification_table$targeted_piecewise_slope_priors) ||
      !all(specification_table$targeted_function_tail_priors) ||
      !identical(
        specification_table$targeted_piecewise_interaction_priors[
          order(specification_table$model_index)
        ],
        c(FALSE, TRUE)
      ) ||
      !all(specification_table$third_remedy) ||
      !all(specification_table$requires_refit)
  ) {
    stop(
      "The third-stage primary remedy registry failed validation.",
      call. = FALSE
    )
  }

  specification_list <- purrr::pmap(specification_table, list)
  names(specification_list) <- paste0(
    specification_table$outcome,
    "_model_",
    stringr::str_pad(specification_table$model_index, 2, pad = "0")
  )
  specification_list
}

select_primary_model_third_remedy_specifications <- function(
    specifications,
    outcome
) {
  selected <- purrr::keep(
    specifications,
    ~ identical(.x$outcome, outcome)
  )
  expected <- if (identical(outcome, "function")) 2L else 0L

  if (length(selected) != expected) {
    stop(
      "Expected ", expected,
      " third-stage remedies for outcome: ", outcome,
      call. = FALSE
    )
  }

  selected
}

select_primary_model_fit_specifications <- function(
    specifications,
    outcome
) {
  selected <- purrr::keep(
    specifications,
    ~ identical(.x$outcome, outcome)
  )

  if (length(selected) != 26L) {
    stop("Expected 26 specifications for outcome: ", outcome, call. = FALSE)
  }

  selected
}

# Fit and assemble candidates ----

# Fit a single specification against the prepared data and priors. Keep the
# model identity and fit settings alongside the fitted object for downstream
# diagnostics.
fit_primary_candidate_model <- function(
    specification,
    priors,
    data,
    chains = 4,
    cores = 4,
    warmup = 2000,
    iter = 8000,
    seed = 1988
) {
  if (iter - warmup != 6000L) {
    stop("Primary fits must retain 6000 iterations per chain.", call. = FALSE)
  }

  required_variables <- c(
    "yi", "vi", "protein_intake", "resistance_exercise_code",
    "resistance_exercise_factor", "study_factor", "study_no",
    "arm_no", "outcome_no", "age_centre", "duration_centre"
  )

  missing_variables <- setdiff(required_variables, names(data))
  if (length(missing_variables) > 0L) {
    stop(
      "The model data are missing: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  if (
    !is.factor(data$resistance_exercise_factor) ||
      !is.factor(data$study_factor)
  ) {
    stop(
      "resistance_exercise_factor and study_factor must both be factors.",
      call. = FALSE
    )
  }

  if (isTRUE(specification$regularised_variance_priors)) {
    priors <- add_regularised_variance_priors(priors)
  }

  model <- brms::brm(
    formula = stats::as.formula(specification$formula),
    data = data,
    prior = priors,
    chains = chains,
    cores = cores,
    seed = seed,
    warmup = warmup,
    iter = iter,
    control = list(
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )

  list(
    outcome = specification$outcome,
    model_index = specification$model_index,
    model_name = specification$model_name,
    specification = specification,
    model = model
  )
}

prepare_primary_remedied_candidate <- function(
    specification,
    priors,
    data
) {
  if (!isTRUE(specification$requires_refit)) {
    stop(
      "Only post-pilot replacement models can use this function.",
      call. = FALSE
    )
  }

  if (isTRUE(specification$regularised_variance_priors)) {
    priors <- add_regularised_variance_priors(priors)
  }

  if (isTRUE(specification$standardised_moderators)) {
    data <- add_primary_standardised_moderators(data)
    priors <- add_primary_standardised_moderator_priors(priors)
  }

  if (isTRUE(specification$half_unit_protein)) {
    data <- add_primary_half_unit_protein(data)
  }

  if (isTRUE(specification$natural_spline)) {
    priors <- remove_global_primary_prior(priors, "sds")
  }

  if (isTRUE(specification$targeted_smooth_prior)) {
    priors <- add_primary_targeted_smooth_prior(
      priors,
      scale = specification$smooth_prior_scale
    )
  }

  list(
    formula = specification$formula,
    priors = priors,
    data = data
  )
}

fit_primary_remedied_candidate_model <- function(
    specification,
    priors,
    data,
    chains = 4,
    cores = 4,
    warmup = 2000,
    iter = 8000,
    seed = 1988
) {
  if (iter - warmup != 6000L) {
    stop("Primary fits must retain 6000 iterations per chain.", call. = FALSE)
  }

  required_variables <- c(
    "yi", "vi", "protein_intake", "resistance_exercise_code",
    "resistance_exercise_factor", "study_factor", "study_no",
    "arm_no", "outcome_no", "age_centre", "duration_centre"
  )

  missing_variables <- setdiff(required_variables, names(data))
  if (length(missing_variables) > 0L) {
    stop(
      "The model data are missing: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  if (
    !is.factor(data$resistance_exercise_factor) ||
      !is.factor(data$study_factor)
  ) {
    stop(
      "resistance_exercise_factor and study_factor must both be factors.",
      call. = FALSE
    )
  }

  prepared <- prepare_primary_remedied_candidate(
    specification = specification,
    priors = priors,
    data = data
  )

  model <- brms::brm(
    formula = stats::as.formula(prepared$formula),
    data = prepared$data,
    prior = prepared$priors,
    chains = chains,
    cores = cores,
    seed = seed,
    warmup = warmup,
    iter = iter,
    control = list(
      adapt_delta = specification$adapt_delta,
      max_treedepth = specification$max_treedepth
    ),
    save_pars = brms::save_pars(all = TRUE)
  )

  list(
    outcome = specification$outcome,
    model_index = specification$model_index,
    model_name = specification$model_name,
    specification = specification,
    model = model
  )
}

fit_primary_second_remedied_candidate_model <- function(
    specification,
    priors,
    data,
    chains = 4,
    cores = 4,
    warmup = 2000,
    iter = 8000,
    seed = 1988
) {
  if (!isTRUE(specification$second_remedy)) {
    stop(
      "Only second-stage primary remedies can use this function.",
      call. = FALSE
    )
  }

  priors <- normalize_primary_brms_priors(priors)

  if (isTRUE(specification$targeted_piecewise_slope_priors)) {
    priors <- add_primary_piecewise_slope_priors(priors)
  }

  # Function model 13 has random-effect SDs but no smooth-term SD. Apply the
  # regularising global SD prior directly so the generic first-stage helper
  # does not add an inapplicable class = "sds" prior.
  fitting_specification <- specification
  if (
    isTRUE(specification$targeted_piecewise_slope_priors) &&
      isTRUE(specification$regularised_variance_priors)
  ) {
    priors <- c(
      priors,
      brms::set_prior("student_t(3, 0, 1)", class = "sd")
    )
    fitting_specification$regularised_variance_priors <- FALSE
  }

  fitted <- fit_primary_remedied_candidate_model(
    specification = fitting_specification,
    priors = priors,
    data = data,
    chains = chains,
    cores = cores,
    warmup = warmup,
    iter = iter,
    seed = seed
  )

  fitted$specification <- specification
  fitted
}

fit_primary_third_remedied_candidate_model <- function(
    specification,
    priors,
    data,
    chains = 4,
    cores = 4,
    warmup = 2000,
    iter = 8000,
    seed = 1988
) {
  if (
    !isTRUE(specification$third_remedy) ||
      !isTRUE(specification$targeted_function_tail_priors)
  ) {
    stop(
      "Only third-stage function remedies can use this function.",
      call. = FALSE
    )
  }

  priors <- add_primary_function_tail_priors(
    priors = priors,
    include_interactions =
      specification$targeted_piecewise_interaction_priors
  )
  priors <- c(
    priors,
    brms::set_prior("student_t(3, 0, 1)", class = "sd")
  )

  fitting_specification <- specification
  fitting_specification$regularised_variance_priors <- FALSE
  fitted <- fit_primary_remedied_candidate_model(
    specification = fitting_specification,
    priors = priors,
    data = data,
    chains = chains,
    cores = cores,
    warmup = warmup,
    iter = iter,
    seed = seed
  )

  fitted$specification <- specification
  fitted
}

merge_primary_candidate_fits <- function(
    existing_candidates,
    replacement_candidates
) {
  if (length(existing_candidates) != 26L) {
    stop("Expected 26 existing candidate models.", call. = FALSE)
  }

  existing_indices <- purrr::map_int(
    existing_candidates,
    "model_index"
  )
  replacement_indices <- purrr::map_int(
    replacement_candidates,
    "model_index"
  )
  outcomes <- unique(c(
    purrr::map_chr(existing_candidates, "outcome"),
    purrr::map_chr(replacement_candidates, "outcome")
  ))
  expected_replacements <- if (
    identical(outcomes, "function")
  ) {
    18L
  } else {
    17L
  }

  if (
    length(outcomes) != 1L ||
      !setequal(existing_indices, seq_len(26L)) ||
      anyDuplicated(replacement_indices) ||
      length(replacement_indices) != expected_replacements ||
      !all(replacement_indices %in% existing_indices) ||
      !all(purrr::map_lgl(
        replacement_candidates,
        ~ isTRUE(.x$specification$requires_refit)
      ))
  ) {
    stop(
      "The primary replacement-model set failed validation.",
      call. = FALSE
    )
  }

  candidates <- stats::setNames(
    existing_candidates,
    existing_indices
  )
  candidates[as.character(replacement_indices)] <-
    replacement_candidates
  candidates <- candidates[order(as.integer(names(candidates)))]

  if (!identical(
    unname(purrr::map_int(candidates, "model_index")),
    1:26
  )) {
    stop(
      "The merged primary model set failed validation.",
      call. = FALSE
    )
  }

  unname(candidates)
}

overlay_primary_candidate_fits <- function(
    existing_candidates,
    replacement_candidates,
    expected_replacements
) {
  if (length(existing_candidates) != 26L) {
    stop("Expected 26 existing candidate models.", call. = FALSE)
  }

  existing_indices <- purrr::map_int(existing_candidates, "model_index")
  replacement_indices <- purrr::map_int(
    replacement_candidates,
    "model_index"
  )
  outcomes <- unique(c(
    purrr::map_chr(existing_candidates, "outcome"),
    purrr::map_chr(replacement_candidates, "outcome")
  ))

  if (
    length(outcomes) != 1L ||
      !setequal(existing_indices, seq_len(26L)) ||
      anyDuplicated(replacement_indices) ||
      length(replacement_indices) != expected_replacements ||
      !all(replacement_indices %in% existing_indices) ||
      !all(purrr::map_lgl(
        replacement_candidates,
        ~ isTRUE(.x$specification$second_remedy)
      ))
  ) {
    stop(
      "The second-stage replacement-model set failed validation.",
      call. = FALSE
    )
  }

  candidates <- stats::setNames(existing_candidates, existing_indices)
  candidates[as.character(replacement_indices)] <- replacement_candidates
  candidates <- candidates[order(as.integer(names(candidates)))]

  if (!identical(
    unname(purrr::map_int(candidates, "model_index")),
    1:26
  )) {
    stop(
      "The second-stage merged model set failed validation.",
      call. = FALSE
    )
  }

  unname(candidates)
}

overlay_primary_third_remedied_candidate_fits <- function(
    existing_candidates,
    replacement_candidates
) {
  if (length(existing_candidates) != 26L) {
    stop("Expected 26 existing candidate models.", call. = FALSE)
  }

  # Dynamic `targets` branches arrive as named lists. Remove those storage
  # names before strict value comparisons so valid model indices are not
  # rejected solely because they carry branch-name attributes.
  existing_indices <- unname(
    purrr::map_int(existing_candidates, "model_index")
  )
  replacement_indices <- unname(purrr::map_int(
    replacement_candidates,
    "model_index"
  ))
  outcomes <- unname(unique(c(
    purrr::map_chr(existing_candidates, "outcome"),
    purrr::map_chr(replacement_candidates, "outcome")
  )))

  if (
    !identical(outcomes, "function") ||
      !setequal(existing_indices, seq_len(26L)) ||
      !identical(sort(replacement_indices), c(13L, 25L)) ||
      anyDuplicated(replacement_indices) ||
      !all(purrr::map_lgl(
        replacement_candidates,
        ~ isTRUE(.x$specification$third_remedy)
      ))
  ) {
    stop(
      "The third-stage replacement-model set failed validation.",
      call. = FALSE
    )
  }

  candidates <- stats::setNames(existing_candidates, existing_indices)
  candidates[as.character(replacement_indices)] <- replacement_candidates
  candidates <- candidates[order(as.integer(names(candidates)))]

  if (!identical(
    unname(purrr::map_int(candidates, "model_index")),
    1:26
  )) {
    stop(
      "The third-stage merged model set failed validation.",
      call. = FALSE
    )
  }

  unname(candidates)
}

assemble_primary_candidate_models <- function(fitted_candidates) {
  if (length(fitted_candidates) != 26L) {
    stop("Expected 26 completed candidate models.", call. = FALSE)
  }

  model_indices <- purrr::map_int(fitted_candidates, "model_index")
  fitted_candidates <- fitted_candidates[order(model_indices)]

  models <- purrr::map(fitted_candidates, "model")
  names(models) <- purrr::map_chr(fitted_candidates, "model_name")

  if (anyDuplicated(names(models)) || !all(purrr::map_lgl(models, inherits, "brmsfit"))) {
    stop("The assembled candidate model list failed validation.", call. = FALSE)
  }

  models
}

safe_diagnostic_minimum <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }

  min(x, na.rm = TRUE)
}

safe_diagnostic_maximum <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }

  max(x, na.rm = TRUE)
}

# Combine sampler diagnostics and posterior predictive checks for one
# candidate, retaining its identity so the results can be compared across
# refit stages.
get_primary_candidate_model_diagnostics <- function(
    fitted_candidate,
    ppc_draws = 400
) {
  model <- fitted_candidate$model
  specification <- fitted_candidate$specification

  nuts <- brms::nuts_params(model) |>
    tibble::as_tibble()

  energy_diagnostics <- nuts |>
    dplyr::filter(Parameter == "energy__") |>
    dplyr::group_by(Chain) |>
    dplyr::summarise(
      ebfmi = mean(diff(Value)^2) / stats::var(Value),
      .groups = "drop"
    )

  divergent_transitions <- nuts |>
    dplyr::filter(Parameter == "divergent__") |>
    dplyr::summarise(value = sum(Value)) |>
    dplyr::pull(value)

  maximum_treedepth_transitions <- nuts |>
    dplyr::filter(Parameter == "treedepth__") |>
    dplyr::summarise(
      value = sum(Value >= specification$max_treedepth)
    ) |>
    dplyr::pull(value)

  model_variables <- posterior::variables(model)
  monitored_variables <- model_variables[
    !stringr::str_detect(
      model_variables,
      "^(lp__|lprior$|log_lik\\[)"
    )
  ]

  parameter_summary <- posterior::as_draws_array(
    model,
    variable = monitored_variables
  ) |>
    posterior::summarise_draws(
      "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"
    ) |>
    tibble::as_tibble() |>
    dplyr::mutate(mcse_sd_ratio = mcse_mean / sd)

  replicated <- brms::posterior_predict(model, ndraws = ppc_draws)
  expected <- brms::posterior_epred(model, ndraws = ppc_draws)
  observed <- model$data$yi
  predictive_sd <- apply(replicated, 2, stats::sd)
  standardised_residual <- (
    observed - colMeans(expected)
  ) / predictive_sd
  pit <- colMeans(sweep(replicated, 2, observed, FUN = "<="))
  predictive_interval <- apply(
    replicated,
    2,
    stats::quantile,
    probs = c(0.05, 0.95)
  )

  summary <- tibble::tibble(
    outcome = fitted_candidate$outcome,
    model_index = fitted_candidate$model_index,
    model = fitted_candidate$model_name,
    source = infer_primary_candidate_source(specification),
    uncorrelated_study_effects =
      specification$uncorrelated_study_effects,
    regularised_variance_priors =
      specification$regularised_variance_priors,
    full_schedule_trial = specification$full_schedule_trial,
    standardised_moderators =
      isTRUE(specification$standardised_moderators),
    half_unit_protein =
      isTRUE(specification$half_unit_protein),
    natural_spline =
      isTRUE(specification$natural_spline),
    targeted_smooth_prior =
      isTRUE(specification$targeted_smooth_prior),
    smooth_prior_scale = if (
      is.null(specification$smooth_prior_scale)
    ) {
      NA_real_
    } else {
      specification$smooth_prior_scale
    },
    adapt_delta = specification$adapt_delta,
    max_treedepth = specification$max_treedepth,
    divergent_transitions = divergent_transitions,
    maximum_treedepth_transitions = maximum_treedepth_transitions,
    minimum_ebfmi = safe_diagnostic_minimum(energy_diagnostics$ebfmi),
    maximum_rhat = safe_diagnostic_maximum(parameter_summary$rhat),
    minimum_bulk_ess = safe_diagnostic_minimum(parameter_summary$ess_bulk),
    minimum_tail_ess = safe_diagnostic_minimum(parameter_summary$ess_tail),
    maximum_mcse_sd_ratio =
      safe_diagnostic_maximum(parameter_summary$mcse_sd_ratio),
    maximum_absolute_standardised_residual =
      safe_diagnostic_maximum(abs(standardised_residual)),
    pit_below_0.025 = sum(pit < 0.025),
    pit_above_0.975 = sum(pit > 0.975),
    predictive_90_interval_coverage = mean(
      observed >= predictive_interval[1, ] &
        observed <= predictive_interval[2, ]
    ),
    geometry_pass = divergent_transitions == 0L &&
      maximum_treedepth_transitions == 0L &&
      minimum_ebfmi >= 0.30,
    convergence_pass = maximum_rhat <= 1.01 &&
      minimum_bulk_ess >= 400 &&
      minimum_tail_ess >= 400 &&
      maximum_mcse_sd_ratio <= 0.10,
    diagnostic_pass = geometry_pass && convergence_pass
  )

  parameter_failures <- parameter_summary |>
    dplyr::filter(
      rhat > 1.01 |
        ess_bulk < 400 |
        ess_tail < 400 |
        mcse_sd_ratio > 0.10
    )

  list(
    summary = summary,
    parameter_summary = parameter_summary,
    parameter_failures = parameter_failures,
    energy_diagnostics = energy_diagnostics,
    standardised_residuals = dplyr::bind_cols(
      tibble::tibble(
        outcome = fitted_candidate$outcome,
        model_index = fitted_candidate$model_index,
        observation = seq_along(observed)
      ),
      tibble::as_tibble(model$data[, intersect(
        c("study_no", "arm_no", "outcome_no"),
        names(model$data)
      ), drop = FALSE]),
      tibble::tibble(
        observed = observed,
        standardised_residual = standardised_residual,
        pit = pit
      )
    )
  )
}

# Historical model averaging and plots ----

# These helpers compute weights, predictions and slopes for the earlier
# workflow. The final primary posterior mixture is built in
# primary_model_averaging.R.

# Functions to generate weights and obtain model averaged preds/slopes
get_model_weights <- function(models) {

  models_nunes <- models

  # pass each to object as loo doesn't like the list for some reason
  m_intercept <- models_nunes$`Intercept only model`
  m_linear <- models_nunes$`Linear model`
  m_lspline <- models_nunes$`Linear spline (1.6g/kg/day) model`
  m_tps <- models_nunes$`Thin plate spline model`
  m_intercept_mods <- models_nunes$`Intercept only model (+ moderators)`
  m_linear_mods <- models_nunes$`Linear model (+ moderators)`
  m_lspline_mods <- models_nunes$`Linear spline (1.6g/kg/day) model (+ moderators)`
  m_tps_mods <- models_nunes$ `Thin plate spline model (+ moderators)`
  m_linear_ranef <- models_nunes$`Linear model (+ random slopes)`
  m_lspline_ranef <- models_nunes$`Linear spline (1.6g/kg/day) model (+ random slopes)`
  m_tps_ranef <- models_nunes$ `Thin plate spline model (+ random smooths)`
  m_linear_mods_ranef <- models_nunes$`Linear model (+ moderators & random slopes)`
  m_lspline_mods_ranef <- models_nunes$`Linear spline (1.6g/kg/day) model (+ moderators & random slopes)`
  m_tps_mods_ranef <- models_nunes$`Thin plate spline model (+ moderators & random smooths)`
  m_linear_interact <- models_nunes$`Linear model (+ RT interaction)`
  m_lspline_interact <- models_nunes$`Linear spline (1.6g/kg/day) model (+ RT interaction)`
  m_tps_interact <- models_nunes$`Thin plate spline model (+ RT interaction smooth)`
  m_linear_mods_interact <- models_nunes$`Linear model (+ moderators & RT interaction)`
  m_lspline_mods_interact <- models_nunes$`Linear spline (1.6g/kg/day) model (+ moderators & RT interaction)`
  m_tps_mods_interact <- models_nunes$`Thin plate spline model (+ moderators & RT interaction)`
  m_linear_ranef_interact <- models_nunes$`Linear model (+ random slopes & RT interaction)`
  m_lspline_ranef_interact <- models_nunes$`Linear spline (1.6g/kg/day) model (+ random slopes & RT interaction)`
  m_tps_ranef_interact <- models_nunes$`Thin plate spline model (+ random smooths & RT interaction smooth)`
  m_linear_mods_ranef_interact <- models_nunes$`Linear model (+ moderators & random slopes & RT interaction)`
  m_lspline_mods_ranef_interact <- models_nunes$`Linear spline (1.6g/kg/day) model (+ moderators & random slopes & RT interaction)`
  m_tps_mods_ranef_interact <- models_nunes$`Thin plate spline model (+ moderators & random smooths & RT interaction smooth)`



  weights_stacking <- loo::loo_model_weights(m_intercept,
                                             m_linear,
                                             m_lspline,
                                             m_tps,
                                             m_intercept_mods,
                                             m_linear_mods,
                                             m_lspline_mods,
                                             m_tps_mods,
                                             m_linear_ranef,
                                             m_lspline_ranef,
                                             m_tps_ranef,
                                             m_linear_mods_ranef,
                                             m_lspline_mods_ranef,
                                             m_tps_mods_ranef,
                                             m_linear_interact,
                                             m_lspline_interact,
                                             m_tps_interact,
                                             m_linear_mods_interact,
                                             m_lspline_mods_interact,
                                             m_tps_mods_interact,
                                             m_linear_ranef_interact,
                                             m_lspline_ranef_interact,
                                             m_tps_ranef_interact,
                                             m_linear_mods_ranef_interact,
                                             m_lspline_mods_ranef_interact,
                                             m_tps_mods_ranef_interact,
                                             method = "stacking")

  return(weights_stacking)
}

get_model_avg_preds <- function(models, weights) {

  models_nunes <- models
  weights_stacking <- weights

  preds <- map(models_nunes, ~ predictions(.x,
                                           re_formula= NA,
                                           newdata   = datagrid(
                                             resistance_exercise_code = c(0,1),
                                             protein_intake = seq(-0.32, 3.28, by = 0.01),
                                             age_centre = 0,
                                             duration_centre = 0,
                                             vi = 0
                                           )))

  weighted_preds <- map2(preds, weights_stacking, function(pred_df, w) {
    pred_df %>%
      mutate(
        .pred_weighted = estimate * w,
        .lower_weighted = conf.low * w,
        .upper_weighted = conf.high * w
      )
  })

  combined_preds <- reduce(weighted_preds, function(df1, df2) {
    df1 %>%
      mutate(
        .pred_weighted = .pred_weighted + df2$.pred_weighted,
        .lower_weighted = .lower_weighted + df2$.lower_weighted,
        .upper_weighted = .upper_weighted + df2$.upper_weighted
      )
  })

  final_predictions <- combined_preds %>%
    mutate(
      estimate = .pred_weighted,
      conf.low = .lower_weighted,
      conf.high = .upper_weighted
    ) |>
    mutate(
      resistance_exercise_label = "Resistance Exercise?",
      resistance_exercise = case_when(
        resistance_exercise_code == 0 ~ "YES",
        resistance_exercise_code == 1 ~ "NO"
      ),
      protein_intake = protein_intake + 1.12
    )

  return(final_predictions)
}

get_model_avg_slopes <- function(models, weights) {

  models_nunes <- models
  weights_stacking <- weights

  # Because can't run slopes(model, variables = "protein_intake") on intercept model with no protein_intake
  safe_slopes <- possibly(slopes, otherwise = NULL)

  slopes <- map(models_nunes, ~ safe_slopes(.x,
                                            re_formula= NA,
                                            variables = "protein_intake",
                                            newdata   = datagrid(
                                              resistance_exercise_code = c(0,1),
                                              protein_intake = seq(-0.32, 3.28, by = 0.42),
                                              age_centre = 0,
                                              duration_centre = 0,
                                              vi = 0
                                            )) |> get_draws())


  # remove nulls from slopes and weights
  slopes_clean <- compact(slopes)
  weights_stacking_clean <- weights_stacking[!map_lgl(slopes, is.null)]

  # Convert each to a matrix (rowid-draw order must match)
  draw_matrices <- map(slopes_clean, ~ matrix(.x$draw, ncol = 1))

  # Stack into matrix: rows = rowid × draw, cols = models
  draw_mat <- do.call(cbind, draw_matrices)  # fast column binding

  # Multiply each column by stacking weight
  weighted_draw_mat <- sweep(draw_mat, 2, weights_stacking_clean, FUN = "*")

  # Sum across models
  stacked_draws <- rowSums(weighted_draw_mat)

  # Reattach draw + rowid info
  meta <- slopes_clean[[1]] |> select(rowid, draw, resistance_exercise_code, protein_intake)

  stacked_df <- bind_cols(meta, estimate = stacked_draws) |>
    mutate(
      resistance_exercise_label = "Resistance Exercise?",
      resistance_exercise = case_when(
        resistance_exercise_code == 0 ~ "YES",
        resistance_exercise_code == 1 ~ "NO"
      ),
      protein_intake = protein_intake + 1.12
    )

  stacked_df
}

get_rope_pd <- function(slopes) {

  rope <- slopes |>
    group_by(protein_intake, resistance_exercise) |>
    nest() |>
    mutate(
      rope_percent = map(data, ~ rope(.x$estimate*0.42, range = c(-0.1, 0.1), ci = 1))
    ) |>
    unnest(rope_percent) |>
    select(protein_intake, resistance_exercise, ROPE_low, ROPE_high, ROPE_Percentage) |>
    filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
    mutate(type = "rope")


  pd <- slopes |>
    group_by(protein_intake, resistance_exercise) |>
    nest() |>
    mutate(
      rope_percent = map(data, ~ rope(.x$estimate*0.42, range = c(0.1, Inf), ci = 1))
    ) |>
    unnest(rope_percent) |>
    select(protein_intake, resistance_exercise, ROPE_low, ROPE_high, ROPE_Percentage) |>
    filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
    mutate(type = "pd")

  rope_pd <- bind_rows(rope, pd)

  return(rope_pd)


}

plot_preds_slopes <- function(data, preds, slopes) {

  preds_plot <- preds |>
    filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
    ggplot(aes(x = protein_intake, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
    geom_point(data = data |>
                 mutate(resistance_exercise_label = "Resistance Exercise?"),
               aes(x = protein_intake + 1.12, y = yi, size = size, color = condition), alpha = 0.5) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.5, color = "black", size = 0.25) +
    geom_line(aes(y = estimate), size = 0.75, color = "black") +
    scale_color_manual(values = c("darkorange", "darkgreen")) +
    ggh4x::facet_nested(.~ resistance_exercise_label + resistance_exercise,
                        scales = "free") +
    labs(
      x = "Protein Intake (g/kg/day)",
      y = "Lean Mass/Muscle Change (SMD)",
      color = "Condition",
      title = "Predicted lean mass/muscle change",
      # subtitle = "Note, red band indicates prior interval for resistance exercise at habitual protein intake, and blue band indicates non-training controls"
    ) +
    guides(
      size = "none"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 10)
    )


  slopes_plot <- slopes |>
    filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
    ggplot(aes(x = protein_intake, y = estimate*0.42)) +
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
    geom_hline(yintercept = c(-0.1,0.1), linetype = "dashed", alpha = 0.75) +
    stat_slabinterval(.width = .95) +
    ggh4x::facet_nested(.~ resistance_exercise_label + resistance_exercise,
                        scales = "free") +
    labs(
    x = "Protein Intake (g/kg/day)",
    y = "Marginal Slope for Lean Mass Change (SMD)",
    title = "Incremental marginal slope of lean mass/muscle change for 0.42 g/kg/day increase",
    subtitle = "Note, 0.42 g/kg/day increase reflects the difference in median protein intake between intervention and control conditions\nDashed lines reflect a smallest effect size of interest of 0.1 SMD"
    ) +
    guides(
      size = "none"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 8),
      axis.title = element_text(size = 10)
    )

  combined_plots <- (preds_plot / slopes_plot) +
    plot_annotation(
      tag_levels = "i",
      title = "Dose-response association of protein intake with lean mass change",
      subtitle = "Bayesian model averaging with LOO-CV stacking weights",
      caption = "Note, reanalysis of Nunes et al. (DOI: 10.1002/jcsm.12922) for purpose of pre-registration of analysis code pipelines"
    ) &
    theme(
      legend.position = "bottom"
    )

  return(combined_plots)
}

# Updated analysis with all re-extracted data ----

## Read and prepare new data ----

# Prepare and audit the extracted workbooks ----

# Use shared exclusion labels for preprocessing and the systematic-review flow
# tables so eligibility decisions can be traced across outputs.
primary_exclusion_labels <- function() {
  c(
    exclusion_missing_mean_change = "missing mean change",
    exclusion_missing_baseline_sd =
      "missing baseline standardisation SD",
    exclusion_missing_sample_size = "missing outcome sample size",
    exclusion_missing_pre_post_correlation =
      "missing pre-post correlation",
    exclusion_effect_size_not_estimable =
      "effect size not estimable",
    exclusion_missing_protein_intake = "missing protein intake",
    exclusion_missing_resistance_exercise_code =
      "missing resistance exercise code",
    exclusion_missing_mean_age = "missing mean age",
    exclusion_missing_duration = "missing duration"
  )
}

read_protein_workbook <- function(path){

  sheets <- readxl::excel_sheets(path)

  data <- purrr::map(
    sheets,
    ~ readxl::read_excel(path, sheet = .x,
                         guess_max = 100000) |>
      janitor::clean_names()
  )

  names(data) <- janitor::make_clean_names(sheets)

  data

}

# Convert the extracted workbook sheets into analysis contrasts. Recover
# missing uncertainty where possible, join baseline and follow-up records, and
# retain the information needed to audit eligibility.
prepare_protein_data <- function(protein_data){

  # rename any variables needed
  studies      <- protein_data$study_characteristics |>
    rename(authors = "authors_e_g_smith_or_smith_and_jones_or_smith_et_al")

  arms         <- protein_data$arm_characteristics |>
    rename(condition = "intervention_i_e_protein_supplement_or_comparator_i_e_placebo_or_control",
           protein_source = "protein_source_placebo_compound_note_copy_paste_directly_from_study_methods",
           absolute_dose = "absolute_protein_dose_g_day_note_if_placebo_condition_then_check_whether_it_includes_any_protein",
           prop_male = "proportion_male_0_1_i_e_all_male_group_1_all_female_group_0")

  observations <- protein_data$arm_based_data_by_timepoint |>
    rename(outcome_group = "outcome_group_whole_body_lean_mass_direct_muscle_size_mass_strength_function",
           week = "timepoint_week_note_baseline_0")

  contrasts    <- protein_data$contrasts_between_timepoints |>
    rename(contrast_weeks = "contrast_e_g_baseline_to_12_weeks_0_12_baseline_to_24_weeks_0_24_12_weeks_to_24_weeks_12_24",
           outcome_group = "outcome_group_whole_body_lean_mass_direct_muscle_size_mass_strength_function") |>

    # separate contrast timepoints
    separate(contrast_weeks, into = c("time1", "time2"), sep = "_") |>
    mutate(time1 = as.numeric(time1),
           time2 = as.numeric(time2))

  # calculate missing sds from ses, p and t stats

  # for arm characteristics
  # variables that have both _sd/_se counterparts
  bases <- names(arms) |>
    str_subset("_(sd|se)$") |>
    str_remove("_(sd|se)$") |>
    unique()

  for (v in bases) {

    sd_col <- paste0(v, "_sd")
    se_col <- paste0(v, "_se")

    # Calculate missing SDs
    arms[[sd_col]] <-
      if_else(
        is.na(arms[[sd_col]]) &
          !is.na(arms[[se_col]]) &
          !is.na(arms$sample_size_at_baseline),
        arms[[se_col]] * sqrt(arms$sample_size_at_baseline),
        arms[[sd_col]]
      )

    # Calculate missing SEs
    arms[[se_col]] <-
      if_else(
        is.na(arms[[se_col]]) &
          !is.na(arms[[sd_col]]) &
          !is.na(arms$sample_size_at_baseline),
        arms[[sd_col]] / sqrt(arms$sample_size_at_baseline),
        arms[[se_col]]
      )

  }

  # for observations
  observations <- observations |>
    mutate(
      # Recover standard_error from 95% CI
      standard_error = coalesce(
        standard_error,
        (upper_95_percent_confidence_interval - lower_95_percent_confidence_interval) /
          (2 * qt(0.975, df = sample_size_at_timepoint - 1))
      ),

      # Recover standard_deviation from standard_error
      standard_deviation = coalesce(
        standard_deviation,
        standard_error * sqrt(sample_size_at_timepoint)
      ),

      # Recover standard_error from standard_deviation (for completeness)
      standard_error = coalesce(
        standard_error,
        standard_deviation / sqrt(sample_size_at_timepoint)
      )
    ) |>

    select(-upper_95_percent_confidence_interval, -lower_95_percent_confidence_interval)

  # Note we do contrasts after joining with observations for the corresponding sample sizes

  # First we generate all possible contrasts with baseline from the observations and join with contrasts extracted

  ## Generate all baseline -> follow-up contrasts
  contrast_grid <-
    observations |>
    filter(!is.na(week)) |>
    group_by(study_no, arm_no, outcome_group, measurement_method, unit_of_measurement) |>
    summarise(
      baseline = list(week[week == 0]),
      followup = list(sort(unique(week[week > 0]))),
      .groups = "drop"
    ) |>
    unnest(baseline) |>
    unnest(followup) |>
    transmute(
      study_no,
      arm_no,
      outcome_group,
      measurement_method,
      unit_of_measurement,
      time1 = baseline,
      time2 = followup
    )

  contrast_keys <- c(
    "study_no",
    "arm_no",
    "outcome_group",
    "measurement_method",
    "unit_of_measurement",
    "time1",
    "time2"
  )

  duplicate_contrasts <-
    contrasts |>
    count(across(all_of(contrast_keys))) |>
    filter(n > 1)

  duplicate_grid_rows <-
    contrast_grid |>
    count(across(all_of(contrast_keys))) |>
    filter(n > 1)

  if (nrow(duplicate_contrasts) > 0) {
    stop(
      "Duplicate reported contrast keys found in ",
      nrow(duplicate_contrasts),
      " rows.",
      call. = FALSE
    )
  }

  if (nrow(duplicate_grid_rows) > 0) {
    stop(
      "Duplicate observation-derived contrast keys found in ",
      nrow(duplicate_grid_rows),
      " rows.",
      call. = FALSE
    )
  }

  contrasts_join <-
    contrast_grid |>
    mutate(arm_timepoint_rows_present = TRUE) |>
    full_join(
      contrasts |>
        mutate(reported_contrast_present = TRUE),
      by = contrast_keys
    ) |>
    mutate(
      arm_timepoint_rows_present = coalesce(
        arm_timepoint_rows_present,
        FALSE
      ),
      reported_contrast_present = coalesce(
        reported_contrast_present,
        FALSE
      ),
      contrast_source = case_when(
        arm_timepoint_rows_present & reported_contrast_present ~
          "reported contrast and arm timepoint rows",
        reported_contrast_present ~ "reported contrast only",
        arm_timepoint_rows_present ~ "arm timepoint rows only"
      )
    )

  lost_reported_contrasts <-
    contrasts |>
    anti_join(contrasts_join, by = contrast_keys)

  if (nrow(lost_reported_contrasts) > 0) {
    stop(
      nrow(lost_reported_contrasts),
      " reported contrasts were lost during data preparation.",
      call. = FALSE
    )
  }

  # Link each side of the contrast to its timepoint so missing uncertainty can
  # be recovered using the corresponding sample sizes.
  baseline_obs <-
    observations |>
    rename(
      time1 = week,
      baseline_mean = mean,
      baseline_sd = standard_deviation,
      baseline_se = standard_error,
      baseline_n = sample_size_at_timepoint
    )


  contrasts_join <-
    contrasts_join |>
    left_join(
      baseline_obs,
      by = c(
        "study_no",
        "arm_no",
        "outcome_group",
        "measurement_method",
        "unit_of_measurement",
        "time1"
      )
    )


  followup_obs <-
    observations |>
    rename(
      time2 = week,
      followup_mean = mean,
      followup_sd = standard_deviation,
      followup_se = standard_error,
      followup_n = sample_size_at_timepoint
    )

  contrasts_join <-
    contrasts_join |>
    left_join(
      followup_obs,
      by = c(
        "study_no",
        "arm_no",
        "outcome_group",
        "measurement_method",
        "unit_of_measurement",
        "time2"
      )
    )

  # for contrasts
  contrasts_analysis <- contrasts_join |>

    # Explicitly nest arm identifiers within studies while retaining the
    # extracted arm number for joins and audit trails.
    group_by(study_no, arm_no) |>
    mutate(arm_no_nested = cur_group_id()) |>
    ungroup() |>

    # Generate one stable outcome identifier and retain this definition
    # throughout preparation.
    group_by(
      study_no,
      arm_no,
      outcome_group,
      measurement_method,
      unit_of_measurement
    ) |>
    mutate(outcome_no = cur_group_id()) |>
    ungroup() |>

    # Use the smaller outcome-specific sample size when both timepoints are
    # available. If one is absent, retain the available timepoint sample size.
    mutate(
      n = case_when(
        !is.na(baseline_n) & !is.na(followup_n) ~
          pmin(baseline_n, followup_n),
        !is.na(followup_n) ~ followup_n,
        !is.na(baseline_n) ~ baseline_n,
        TRUE ~ NA_real_
      )
    ) |>

    # Add missing contrast means
    mutate(
      mean = coalesce(
        mean,
        followup_mean - baseline_mean
      )
    ) |>
    rename(
      mean_change = "mean"
    ) |>

    mutate(
      # Recover standard_error from 95% CI
      standard_error = coalesce(
        standard_error,
        (upper_95_percent_confidence_interval - lower_95_percent_confidence_interval) /
          (2 * qt(0.975, df = n - 1))
      ),

      # # Recover SE from reported t statistic
      # standard_error = coalesce(
      #   standard_error,
      #   abs(mean / t_statistic)
      # ),
      #
      # # Recover SE from reported two-sided p value
      # standard_error = coalesce(
      #   standard_error,
      #   abs(mean) /
      #     qt(
      #       p_value / 2,
      #       df = n - 1,
      #       lower.tail = FALSE
      #     )
      # ),

      # Recover standard_deviation from standard_error
      standard_deviation = coalesce(
        standard_deviation,
        standard_error * sqrt(n)
      ),

      # Recover standard_error from standard_deviation (for completeness)
      standard_error = coalesce(
        standard_error,
        standard_deviation / sqrt(n)
      )
    ) |>

    # Note whether change SD is present
    mutate(
      has_change_sd = !is.na(standard_deviation)
    ) |>

    # Calculate pre-post correlations
    mutate(
      r = case_when(
        !is.na(baseline_sd) & !is.na(followup_sd) & !is.na(standard_deviation) ~
        (
          baseline_sd^2 +
            followup_sd^2 -
            standard_deviation^2
        ) /
        (
          2 * baseline_sd * followup_sd
        )
    )
    ) |>
    mutate(
      r = pmin(0.999, pmax(-0.999, r)) # Constrain it in case of errors
    )

  # Meta-analyse pre-post correlations to impute missing ones
  # convert using Fishers r to z
  contrasts_analysis <- escalc(measure = "ZCOR", ri = r, ni = n, data = contrasts_analysis)

  meta_r <- rma.mv(yi, V=vi, data=contrasts_analysis,
                         random = list(~ 1 | study_no, ~ 1 | arm_no_nested, ~ 1 | outcome_no), method="REML", test="t",
                         control=list(optimizer="optim", optmethod="Nelder-Mead"))

  robu_meta_r <- robust(meta_r, contrasts_analysis$study_no)

  z2r_RT <- psych::fisherz2r(robu_meta_r$b)
  z2r_RT_lower <- psych::fisherz2r(robu_meta_r$ci.lb)
  z2r_RT_upper <- psych::fisherz2r(robu_meta_r$ci.ub)

  contrasts_analysis <- contrasts_analysis |>
    mutate(
      r = case_when(
        is.na(r) ~ z2r_RT,
        .default = r
      )
    )

  # Calculate pooled sd for each outcome and then becker's d for each contrast

  baseline_pool <-
    contrasts_analysis |>
    filter(
      !is.na(baseline_n),
      baseline_n > 1,
      !is.na(baseline_sd),
      baseline_sd > 0
    ) |>
    distinct(
      study_no,
      outcome_no,
      arm_no,
      outcome_group,
      measurement_method,
      unit_of_measurement,
      baseline_n,
      baseline_sd
    ) |>
    group_by(study_no,
             outcome_group,
             measurement_method,
             unit_of_measurement) |>
    summarise(
      baseline_sd_pooled = sqrt(
        sum((baseline_n - 1) * baseline_sd^2) /
          sum(baseline_n - 1)
      ),
      .groups = "drop"
    )

  contrasts_analysis <-
    contrasts_analysis |>
    left_join(
      baseline_pool,
      by = c("study_no",
             "outcome_group",
             "measurement_method",
             "unit_of_measurement")
    ) |>

    # add zero to just use reported or imputed change for numerator
    mutate(
      zero = 0
    )

  contrasts_analysis <- escalc(
    measure="SMCR",
    m1i=mean_change,
    m2i=zero,
    sd1i=baseline_sd_pooled,
    ni=n,
    ri=r,
    data = contrasts_analysis
  )

  # Now rejoin with study and arm based data for dose and other characteristics
  contrasts_analysis <-
    contrasts_analysis |>
    left_join(
      studies,
      by = "study_no"
    )

  contrasts_analysis <-
    contrasts_analysis |>
    left_join(
      arms |>
        select(
          study_no,
          arm_no,
          condition,
          trained_or_untrained,
          training_experience_years,
          prop_male,
          mean_age_years,
          sd_age_years,
          mean_body_mass_kg,
          sd_body_mass_kg,
          resistance_exercise_yes_no,
          intervention_dietary_protein_intake_g_kg_day,
          intervention_dietary_protein_intake_g_kg_day_se
        ),
      by = c("study_no", "arm_no")
    ) |>
    rename(protein_intake_raw = "intervention_dietary_protein_intake_g_kg_day",
           protein_intake_se = "intervention_dietary_protein_intake_g_kg_day_se")

  # centre  predictors

  # first get protein median from controls/placebos
  protein_median_controls_placebos <- contrasts_analysis |>
    filter(condition == "control" | condition == "placebo") |>
    group_by(study_no, arm_no) |>
    slice_head(n = 1) |>
    ungroup() |>
    summarise(
      protein_median = median(protein_intake_raw, na.rm = TRUE)
    )

  # Retain a reason for each excluded row so the same eligibility decisions
  # feed the model data and review flow counts.
  exclusion_labels <- primary_exclusion_labels()


  contrasts_analysis <- contrasts_analysis |>
  mutate(
    protein_median = protein_median_controls_placebos$protein_median[1],
    protein_intake = protein_intake_raw - protein_median,
    duration_centre = time2 - 12,
    age_centre = mean_age_years - 25,
    resistance_exercise_code = case_when(
      resistance_exercise_yes_no == "yes" ~ 0,
      resistance_exercise_yes_no == "no" ~ 1,
    ),
    # Add explicit factor coding for factor-by smooths while retaining the
    # numeric coding used for the other interaction models
    resistance_exercise_factor = factor(
      resistance_exercise_yes_no,
      levels = c("yes", "no")
    )
  ) |>

    # Make the nested arm identifier the analysis identifier while preserving
    # the workbook arm number explicitly for traceability.
    rename(
      arm_no_extracted = arm_no,
      arm_no = arm_no_nested
    ) |>

    # Record every reason retained rows are not eligible for the primary
    # complete-case analysis. The individual flags support totals by reason;
    # the combined field supports row-level auditing.
    mutate(
      exclusion_missing_mean_change = is.na(mean_change),
      exclusion_missing_baseline_sd = is.na(baseline_sd_pooled),
      exclusion_missing_sample_size = is.na(n),
      exclusion_missing_pre_post_correlation = is.na(r),
      exclusion_effect_size_not_estimable =
        is.na(yi) | is.na(vi) | coalesce(vi <= 0, FALSE),
      exclusion_missing_protein_intake = is.na(protein_intake),
      exclusion_missing_resistance_exercise_code =
        is.na(resistance_exercise_code),
      exclusion_missing_mean_age = is.na(age_centre),
      exclusion_missing_duration = is.na(duration_centre)
    ) |>
    rowwise() |>
    mutate(
      analysis_exclusion_count = sum(
        c_across(all_of(names(exclusion_labels)))
      ),
      analysis_exclusion_reasons = {
        exclusion_flags <- c_across(all_of(names(exclusion_labels)))
        reasons <- unname(exclusion_labels[exclusion_flags])

        if (length(reasons) == 0) {
          NA_character_
        } else {
          paste(reasons, collapse = "; ")
        }
      },
      analysis_eligible_primary = analysis_exclusion_count == 0
    ) |>
    ungroup()

  nested_arm_count <- n_distinct(contrasts_analysis$arm_no)
  extracted_study_arm_count <-
    contrasts_analysis |>
    distinct(study_no, arm_no_extracted) |>
    nrow()

  if (nested_arm_count != extracted_study_arm_count) {
    stop(
      "Nested arm identifiers do not map one-to-one to extracted study-arms.",
      call. = FALSE
    )
  }

  return(contrasts_analysis)

}

# Aggregate the prepared observations and exclusions into review counts and
# outcome-level audit tables.
summarise_systematic_review_data <- function(
  data,
  additional_exclusion_labels = NULL
) {
  required_columns <- c(
    "study_no", "arm_no_extracted", "outcome_group", "condition", "n",
    "time2", "analysis_eligible_primary", "analysis_exclusion_reasons"
  )
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns)) {
    stop(
      "Systematic-review summary is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  outcome_labels <- c(
    "whole body lean mass" = "Whole-body lean/fat-free mass",
    "direct muscle size/mass" = "Measures of muscle size or mass",
    "strength" = "Strength",
    "function" = "Physical function"
  )
  outcome_levels <- names(outcome_labels)
  condition_levels <- c("control", "placebo", "protein")

  if (
    !is.null(additional_exclusion_labels) &&
      (is.null(names(additional_exclusion_labels)) ||
        any(names(additional_exclusion_labels) == "") ||
        anyDuplicated(names(additional_exclusion_labels)))
  ) {
    stop("Additional exclusion labels must be a uniquely named vector.",
         call. = FALSE)
  }

  exclusion_labels <- c(
    primary_exclusion_labels(),
    additional_exclusion_labels
  )
  exclusion_flags <- names(exclusion_labels)

  missing_flags <- setdiff(exclusion_flags, names(data))
  if (length(missing_flags)) {
    stop(
      "Systematic-review summary is missing exclusion flags: ",
      paste(missing_flags, collapse = ", "),
      call. = FALSE
    )
  }

  eligible <- data |>
    dplyr::filter(analysis_eligible_primary)

  if (!nrow(eligible)) {
    stop("No contrasts are eligible for the primary analysis.", call. = FALSE)
  }

  # A participant is counted once per contributing arm and outcome using the
  # largest analysed sample size for that arm. This avoids multiplying the
  # same participants by measurement, unit, or follow-up while retaining an
  # outcome-specific count when outcome availability differs within an arm.
  arm_outcome <- eligible |>
    dplyr::group_by(
      outcome_group, study_no, arm_no_extracted, condition
    ) |>
    dplyr::summarise(
      participants = max(n, na.rm = TRUE),
      .groups = "drop"
    )

  arm_overall <- arm_outcome |>
    dplyr::group_by(study_no, arm_no_extracted, condition) |>
    dplyr::summarise(
      participants = max(participants),
      .groups = "drop"
    )

  outcome_core <- eligible |>
    dplyr::group_by(outcome_group) |>
    dplyr::summarise(
      studies = dplyr::n_distinct(study_no),
      contrasts = dplyr::n(),
      .groups = "drop"
    )

  outcome_arm_counts <- arm_outcome |>
    dplyr::count(outcome_group, condition, name = "arms") |>
    tidyr::complete(
      outcome_group = outcome_levels,
      condition = condition_levels,
      fill = list(arms = 0L)
    )

  outcome_participant_counts <- arm_outcome |>
    dplyr::group_by(outcome_group, condition) |>
    dplyr::summarise(
      participants = sum(participants),
      .groups = "drop"
    ) |>
    tidyr::complete(
      outcome_group = outcome_levels,
      condition = condition_levels,
      fill = list(participants = 0)
    )

  outcome_duration <- eligible |>
    dplyr::group_by(outcome_group, study_no) |>
    dplyr::summarise(
      duration_weeks = max(time2),
      .groups = "drop"
    ) |>
    dplyr::group_by(outcome_group) |>
    dplyr::summarise(
      duration_median_weeks = stats::median(duration_weeks),
      duration_q1_weeks = as.numeric(
        stats::quantile(duration_weeks, probs = 0.25, names = FALSE)
      ),
      duration_q3_weeks = as.numeric(
        stats::quantile(duration_weeks, probs = 0.75, names = FALSE)
      ),
      .groups = "drop"
    )

  by_outcome <- outcome_core |>
    dplyr::left_join(
      outcome_arm_counts |>
        tidyr::pivot_wider(
          names_from = condition,
          values_from = arms,
          names_prefix = "arms_"
        ),
      by = "outcome_group"
    ) |>
    dplyr::left_join(
      outcome_participant_counts |>
        tidyr::pivot_wider(
          names_from = condition,
          values_from = participants,
          names_prefix = "participants_"
        ),
      by = "outcome_group"
    ) |>
    dplyr::left_join(outcome_duration, by = "outcome_group") |>
    dplyr::mutate(
      outcome = unname(outcome_labels[outcome_group]),
      arms = arms_control + arms_placebo + arms_protein,
      participants =
        participants_control + participants_placebo + participants_protein
    ) |>
    dplyr::arrange(match(outcome_group, outcome_levels)) |>
    dplyr::select(
      outcome_group, outcome, studies, arms,
      arms_control, arms_placebo, arms_protein,
      participants,
      participants_control, participants_placebo, participants_protein,
      contrasts,
      duration_median_weeks, duration_q1_weeks, duration_q3_weeks
    )

  overall_duration <- eligible |>
    dplyr::group_by(study_no) |>
    dplyr::summarise(
      duration_weeks = max(time2),
      .groups = "drop"
    ) |>
    dplyr::summarise(
      duration_median_weeks = stats::median(duration_weeks),
      duration_q1_weeks = as.numeric(
        stats::quantile(duration_weeks, probs = 0.25, names = FALSE)
      ),
      duration_q3_weeks = as.numeric(
        stats::quantile(duration_weeks, probs = 0.75, names = FALSE)
      )
    )

  overall_condition <- arm_overall |>
    dplyr::group_by(condition) |>
    dplyr::summarise(
      arms = dplyr::n(),
      participants = sum(participants),
      .groups = "drop"
    ) |>
    tidyr::complete(
      condition = condition_levels,
      fill = list(arms = 0L, participants = 0)
    )

  condition_value <- function(variable, condition) {
    overall_condition[[variable]][
      match(condition, overall_condition$condition)
    ]
  }

  overall <- tibble::tibble(
    outcome_group = "overall",
    outcome = "Overall",
    studies = dplyr::n_distinct(eligible$study_no),
    arms = nrow(arm_overall),
    arms_control = condition_value("arms", "control"),
    arms_placebo = condition_value("arms", "placebo"),
    arms_protein = condition_value("arms", "protein"),
    participants = sum(arm_overall$participants),
    participants_control = condition_value("participants", "control"),
    participants_placebo = condition_value("participants", "placebo"),
    participants_protein = condition_value("participants", "protein"),
    contrasts = nrow(eligible),
    duration_median_weeks = overall_duration$duration_median_weeks,
    duration_q1_weeks = overall_duration$duration_q1_weeks,
    duration_q3_weeks = overall_duration$duration_q3_weeks
  )

  by_outcome_display <- by_outcome |>
    dplyr::transmute(
      Outcome = outcome,
      Studies = studies,
      Arms = sprintf(
        "%s (%s/%s/%s)",
        arms, arms_control, arms_placebo, arms_protein
      ),
      Participants = sprintf(
        "%s (%s/%s/%s)",
        participants,
        participants_control,
        participants_placebo,
        participants_protein
      ),
      Contrasts = contrasts,
      `Duration, weeks` = sprintf(
        "%s [%s, %s]",
        duration_median_weeks, duration_q1_weeks, duration_q3_weeks
      )
    )

  candidate_rows <- data |>
    dplyr::count(outcome_group, name = "candidate_contrasts")

  excluded_rows <- data |>
    dplyr::arrange(
      match(outcome_group, outcome_levels),
      study_no, arm_no_extracted, measurement_method,
      unit_of_measurement, time1, time2
    ) |>
    dplyr::mutate(
      contrast_candidate_id = sprintf("C%04d", dplyr::row_number())
    ) |>
    dplyr::filter(!analysis_eligible_primary) |>
    dplyr::select(
      contrast_candidate_id,
      dplyr::any_of(c(
        "study_no", "authors", "year_of_publication", "doi",
        "arm_no_extracted", "condition", "outcome_group",
        "measurement_method", "unit_of_measurement", "time1", "time2",
        "contrast_source", "n", "mean_change", "baseline_sd_pooled", "r",
        "yi", "vi", "protein_intake_raw", "resistance_exercise_yes_no",
        "mean_age_years", "analysis_exclusion_count",
        "analysis_exclusion_reasons", exclusion_flags
      ))
    ) |>
    # metafor may retain one-column matrix attributes on calculated columns.
    # Flatten them only in the publication-facing audit so CSV output remains
    # rectangular without changing the analysis object or model inputs.
    dplyr::mutate(
      r = as.numeric(r),
      exclusion_missing_pre_post_correlation =
        as.logical(exclusion_missing_pre_post_correlation)
    )

  exclusion_reasons <- excluded_rows |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(exclusion_flags),
      names_to = "reason_flag",
      values_to = "applies"
    ) |>
    dplyr::filter(applies) |>
    dplyr::mutate(
      reason = unname(exclusion_labels[reason_flag])
    )

  exclusion_reason_overall <- exclusion_reasons |>
    dplyr::count(reason_flag, name = "excluded_contrasts") |>
    tidyr::complete(
      reason_flag = exclusion_flags,
      fill = list(excluded_contrasts = 0L)
    ) |>
    dplyr::mutate(
      reason = unname(exclusion_labels[reason_flag]),
      percent_of_excluded_contrasts =
        100 * excluded_contrasts / nrow(excluded_rows)
    ) |>
    dplyr::select(
      reason_flag, reason, excluded_contrasts,
      percent_of_excluded_contrasts
    ) |>
    dplyr::arrange(dplyr::desc(excluded_contrasts), reason)

  exclusion_reason_by_outcome <- exclusion_reasons |>
    dplyr::count(
      outcome_group, reason_flag,
      name = "excluded_contrasts"
    ) |>
    tidyr::complete(
      outcome_group = outcome_levels,
      reason_flag = exclusion_flags,
      fill = list(excluded_contrasts = 0L)
    ) |>
    dplyr::left_join(
      excluded_rows |>
        dplyr::count(outcome_group, name = "all_excluded_contrasts"),
      by = "outcome_group"
    ) |>
    dplyr::mutate(
      outcome = unname(outcome_labels[outcome_group]),
      reason = unname(exclusion_labels[reason_flag]),
      percent_of_excluded_contrasts =
        100 * excluded_contrasts / all_excluded_contrasts
    ) |>
    dplyr::select(
      outcome_group, reason_flag, reason, excluded_contrasts,
      all_excluded_contrasts, outcome, percent_of_excluded_contrasts
    ) |>
    dplyr::arrange(
      match(outcome_group, outcome_levels),
      dplyr::desc(excluded_contrasts), reason
    )

  flow_by_outcome <- candidate_rows |>
    dplyr::left_join(
      eligible |>
        dplyr::count(outcome_group, name = "included_contrasts"),
      by = "outcome_group"
    ) |>
    dplyr::left_join(
      excluded_rows |>
        dplyr::count(outcome_group, name = "excluded_contrasts"),
      by = "outcome_group"
    ) |>
    dplyr::mutate(
      outcome = unname(outcome_labels[outcome_group]),
      included_contrasts = tidyr::replace_na(included_contrasts, 0L),
      excluded_contrasts = tidyr::replace_na(excluded_contrasts, 0L)
    ) |>
    dplyr::arrange(match(outcome_group, outcome_levels))

  list(
    overall = overall,
    by_outcome = by_outcome,
    by_outcome_display = by_outcome_display,
    flow_by_outcome = flow_by_outcome,
    exclusion_reason_overall = exclusion_reason_overall,
    exclusion_reason_by_outcome = exclusion_reason_by_outcome,
    excluded_rows = excluded_rows,
    definitions = c(
      "Studies are trials with at least one contrast eligible for the primary analysis.",
      paste(
        "Arms are distinct study-arm combinations with at least one eligible",
        "contrast and are classified as control, placebo, or protein."
      ),
      paste(
        "Participants are counted once per contributing arm using the largest",
        "analysed sample for that arm; outcome-specific participant counts may",
        "therefore overlap and should not be summed across outcomes."
      ),
      paste(
        "Contrasts are eligible baseline-to-follow-up outcome rows; multiple",
        "measurements or follow-ups from one arm are counted separately."
      ),
      paste(
        "Duration is the longest eligible follow-up per study, summarised as",
        "median [first quartile, third quartile] in weeks."
      ),
      paste(
        "Exclusion reasons are non-mutually-exclusive, so reason-specific",
        "counts can sum to more than the number of excluded contrasts."
      )
    )
  )
}

write_systematic_review_audit_files <- function(summary) {
  output_dir <- here::here("supplement", "data")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_files <- c(
    summary = file.path(output_dir, "systematic-review-summary.csv"),
    flow = file.path(output_dir, "contrast-flow-by-outcome.csv"),
    exclusion_reasons =
      file.path(output_dir, "contrast-exclusion-reasons.csv"),
    exclusion_reasons_by_outcome =
      file.path(output_dir, "contrast-exclusion-reasons-by-outcome.csv"),
    excluded_contrasts =
      file.path(output_dir, "excluded-contrast-audit.csv")
  )

  readr::write_csv(
    dplyr::bind_rows(summary$overall, summary$by_outcome),
    output_files[["summary"]]
  )
  readr::write_csv(summary$flow_by_outcome, output_files[["flow"]])
  readr::write_csv(
    summary$exclusion_reason_overall,
    output_files[["exclusion_reasons"]]
  )
  readr::write_csv(
    summary$exclusion_reason_by_outcome,
    output_files[["exclusion_reasons_by_outcome"]]
  )
  readr::write_csv(
    summary$excluded_rows,
    output_files[["excluded_contrasts"]]
  )

  unname(output_files)
}

# Use the historical strength and function results to construct prior sets on
# the corresponding outcome scales.

## Fitting old pre-reg brms models ----
get_steele_strength_priors <- function(data) {

  # uses data from Steele et al. (https://doi.org/10.1080/02640414.2023.2286748) for strength

  steele_data <- data |>
    mutate(
      # calculate pre-post SDs from SEs
      RT_pre_sd = if_else(is.na(RT_pre_se), RT_pre_sd, RT_pre_se * sqrt(RT_n)),
      CON_pre_sd = if_else(is.na(CON_pre_se), CON_pre_sd, CON_pre_se * sqrt(CON_n)),
      RT_post_sd = if_else(is.na(RT_post_se), RT_post_sd, RT_post_se * sqrt(RT_n)),
      CON_post_sd = if_else(is.na(CON_post_se), CON_post_sd, CON_post_se * sqrt(CON_n)),

      # convert p to t (Change scores)
      RT_delta_t_value = if_else(is.na(RT_delta_t_value), qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE), RT_delta_t_value),
      CON_delta_t_value = if_else(is.na(CON_delta_t_value), qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE), CON_delta_t_value),

      # convert t to SE (Change scores)
      RT_delta_se = if_else(is.na(RT_delta_se), if_else(is.na(RT_delta_m), (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value), RT_delta_se),
      CON_delta_se = if_else(is.na(CON_delta_se), if_else(is.na(CON_delta_m), (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value), CON_delta_se),

      # make positive
      RT_delta_se = if_else(RT_delta_se < 0, RT_delta_se * -1, RT_delta_se),
      CON_delta_se = if_else(CON_delta_se < 0, CON_delta_se * -1, CON_delta_se),

      # convert CI to SE (Change scores)
      RT_delta_se = if_else(is.na(RT_delta_se), (RT_delta_CI_upper - RT_delta_CI_lower)/3.92, RT_delta_se),
      CON_delta_se = if_else(is.na(CON_delta_se), (CON_delta_CI_upper - CON_delta_CI_lower)/3.92, CON_delta_se),

      # convert SE to SD (Change scores)
      RT_delta_sd = if_else(is.na(RT_delta_sd), RT_delta_se * sqrt(RT_n), RT_delta_sd),
      CON_delta_sd = if_else(is.na(CON_delta_sd), CON_delta_se * sqrt(CON_n), CON_delta_sd),

      # calculate pre-post correlation coefficient for those with pre, post, and delta SDs
      RT_ri = (RT_pre_sd^2 + RT_post_sd^2 - RT_delta_sd^2)/(2 * RT_pre_sd * RT_post_sd),
      CON_ri = (CON_pre_sd^2 + CON_post_sd^2 - CON_delta_sd^2)/(2 * CON_pre_sd * CON_post_sd)
    ) |>
    mutate(
      # Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
      RT_ri = if_else(between(RT_ri,-1,1) == FALSE, NA, RT_ri),
      CON_ri = if_else(between(CON_ri,-1,1) == FALSE, NA, CON_ri)
    ) |>
    mutate(pre_sd_pool = sqrt(((RT_n - 1) * RT_pre_sd ^ 2 +
                                 (CON_n - 1) * CON_pre_sd ^ 2) /  (RT_n + CON_n - 2)))

  # Calculate pre-post correlations for SMCR calculations

  steele_data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = steele_data)

  meta_RT_ri <- rma.mv(yi, V=vi, data=steele_data,
                       random = list(~ 1 | study, ~1 | arm, ~1 | es), method="REML", test="t",
                       control=list(optimizer="optim", optmethod="Nelder-Mead"))

  rob_meta_RT_ri <- robust(meta_RT_ri, steele_data$study)

  z2r_RT <- psych::fisherz2r(rob_meta_RT_ri$b[1])

  steele_data$RT_ri <- if_else(is.na(steele_data$RT_ri), z2r_RT, steele_data$RT_ri)

  steele_data_CON <- steele_data |>
    select(study, arm, es, CON_ri, CON_n) |>
    distinct()

  steele_data_CON <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = steele_data_CON)

  meta_CON_ri <- rma.mv(yi, V=vi, data=steele_data_CON,
                        random = list(~ 1 | study, ~1 | arm, ~1 | es), method="REML", test="t",
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))

  rob_meta_CON_ri <- robust(meta_CON_ri, steele_data_CON$study)

  z2r_CON <- psych::fisherz2r(rob_meta_CON_ri$b[1])

  steele_data$CON_ri <- if_else(is.na(steele_data$CON_ri), z2r_CON, steele_data$CON_ri)


  # calculate for pre-post and delta separately

  # Pivot longer so in arm based format

  steele_data <- steele_data |>
    select(-yi, -vi, -RT_only) |>
    pivot_longer(
      cols = contains(c("RT_", "CON_")),
      names_to = "what",
      values_to = "value"
    ) |>
    separate("what", into = c("condition", "what"), sep = "_", extra = "merge") |>
    pivot_wider(names_from = "what",
                values_from = "value")

  # calculate for pre-post and delta separately
  steele_data_pre_post <- steele_data |>
    filter(!is.na(post_m))

  steele_data_pre_post <- escalc(
    measure = "SMCR",
    m1i = post_m,
    m2i = pre_m,
    sd1i = pre_sd_pool,
    # sd2i = pre_sd,
    ni = n,
    ri = ri,
    data = steele_data_pre_post
  )

  steele_data_joined_delta <- steele_data |>
    filter(is.na(post_m)) |>
    mutate(ref = 0)

  data_joined_delta <- escalc(
    measure = "SMCR",
    m1i = delta_m,
    m2i = ref,
    sd1i = pre_sd_pool,
    # sd2i = pre_sd,
    ni = n,
    ri = ri,
    data = steele_data_joined_delta
  )

  steele_data_joined <- bind_rows(steele_data_pre_post, steele_data_joined_delta) |>
    filter(outcome == "strength") |>
    unite("arm", c(outcome,study,arm,condition), sep = "_", remove = FALSE) |>
    mutate(arm = dense_rank(arm)) |>
    rowid_to_column("effect") |>
    mutate(condition_code = case_when(
      condition == "RT" ~ 0,
      condition == "CON" ~ 1,
    ))

  # frequentist meta-analysis estimates
  steele_meta <- rma.mv(yi, vi,
                        random = list(~ 1 | study, ~1 | arm, ~1 | effect),
                        mods = ~ condition_code,
                        data = steele_data_joined,
                        method="REML", test="t",
                        # control=list(optimizer="optim", optmethod="Nelder-Mead")
  )

  rob_steele_meta <- robust(steele_meta, cluster = steele_data_joined$study, clubSandwich = TRUE)

  # set priors based on estimates (assumes habitual protein intake to get effects of RT and CON at that level)

  priors <-
    c(
      # # We set priors based on the overall estimates for both control and RT conditions in an arm based model using Steele et al. data
      set_prior(paste("student_t(3,", rob_steele_meta$b[1],",", rob_steele_meta$se[1],")"), class = "b", coef = "Intercept"),
      set_prior(paste("student_t(3,", rob_steele_meta$b[2],",", rob_steele_meta$se[2],")"), class = "b", coef = "resistance_exercise_code"),
      # All other b parameters are kept as uninformative but weakly regularising
      set_prior("student_t(3,0,1)", class = "b")
      # All other priors for variance parameters are kept as default weakly regularising
    )

  priors

}

get_kirk_function_priors <- function() {

  # Taken directly from posterior estimates reported in Kirk et al. (https://doi.org/10.3390/jfmk9040239) for functional outcomes
  # as arm based model already fit
  # set priors based on estimates (assumes habitual protein intake to get effects of RT and CON at that level)

  priors <-
    c(
      # # We set priors based on the overall estimates for both control and RT conditions in an arm based model using Steele et al. data
        set_prior("student_t(3, 0.72, 0.18)", class = "b", coef = "Intercept"),
      set_prior("student_t(3, -0.63, 0.205)", class = "b", coef = "resistance_exercise_code"),
      # All other b parameters are kept as uninformative but weakly regularising
      set_prior("student_t(3,0,1)", class = "b")
      # All other priors for variance parameters are kept as default weakly regularising
    )

  priors

}



## Fitting error in variables models using stan ----

# Stan inputs for measurement-error models ----

# Translate the prior and analysis representations into the numeric fields
# expected by the supplementary Stan programs.

# Helper functions to take old brms prior and extract to save rerunning earlier parts of pipeline
extract_student_t_prior <- function(prior_text) {
  matches <- stringr::str_match(
    stringr::str_remove_all(prior_text, "\\s+"),
    "^student_t\\(([-+0-9.eE]+),([-+0-9.eE]+),([-+0-9.eE]+)\\)$"
  )

  if (anyNA(matches)) {
    stop("Expected a student_t(df, location, scale) prior.")
  }

  c(
    df = as.numeric(matches[2]),
    location = as.numeric(matches[3]),
    scale = as.numeric(matches[4])
  )
}

make_stan_priors <- function(steele_priors) {
  prior_table <- as.data.frame(steele_priors)

  get_prior <- function(coefficient) {
    prior_text <- prior_table |>
      dplyr::filter(
        class == "b",
        coef == coefficient
      ) |>
      dplyr::pull(prior)

    if (length(prior_text) != 1L) {
      stop(
        "Expected exactly one prior for coefficient: ",
        coefficient
      )
    }

    extract_student_t_prior(prior_text)
  }

  intercept_prior <- get_prior("Intercept")
  rt_prior <- get_prior("resistance_exercise_code")

  list(
    prior_intercept_location = unname(intercept_prior["location"]),
    prior_intercept_scale = unname(intercept_prior["scale"]),
    prior_rt_location = unname(rt_prior["location"]),
    prior_rt_scale = unname(rt_prior["scale"])
  )
}

# Prepare data for stan

make_stan_data <- function(
    data,
    study_id_col = "study_no"
) {
  required_cols <- c(
    study_id_col,
    "arm_no",
    "outcome_no",
    "yi",
    "vi",
    "protein_intake",
    "protein_intake_se",
    "resistance_exercise_code",
    "duration_centre",
    "age_centre"
  )

  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "data is missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  data <- data |>
    dplyr::mutate(study_id = .data[[study_id_col]])

  # Each study–arm supplies one observed intake and one intake SE.
  # Do not silently continue if duplicated rows disagree.
  arm_checks <- data |>
    dplyr::group_by(study_id, arm_no) |>
    dplyr::summarise(
      n_protein_intake = dplyr::n_distinct(protein_intake),
      n_protein_intake_se = dplyr::n_distinct(protein_intake_se),
      .groups = "drop"
    )

  bad_arms <- arm_checks |>
    dplyr::filter(
      n_protein_intake != 1L |
        n_protein_intake_se != 1L
    )

  if (nrow(bad_arms) > 0) {
    stop(
      "Each study–arm must have one protein_intake and one ",
      "protein_intake_se. Found inconsistent values in ",
      nrow(bad_arms),
      " study–arms.",
      call. = FALSE
    )
  }

  if (anyNA(data$protein_intake_se) || any(data$protein_intake_se <= 0)) {
    stop(
      "protein_intake_se must be present and strictly positive.",
      call. = FALSE
    )
  }

  if (anyNA(data$vi) || any(data$vi <= 0)) {
    stop("vi must be present and strictly positive.", call. = FALSE)
  }

  data_indexed <- data |>
    dplyr::arrange(study_id, arm_no, outcome_no) |>
    dplyr::group_by(study_id) |>
    dplyr::mutate(study_index = dplyr::cur_group_id()) |>
    dplyr::ungroup() |>
    dplyr::group_by(study_id, arm_no) |>
    dplyr::mutate(arm_index = dplyr::cur_group_id()) |>
    dplyr::ungroup() |>
    dplyr::group_by(study_id, arm_no, outcome_no) |>
    dplyr::mutate(outcome_index = dplyr::cur_group_id()) |>
    dplyr::ungroup()

  arm_data <- data_indexed |>
    dplyr::distinct(
      study_id,
      arm_no,
      arm_index,
      protein_intake,
      protein_intake_se
    ) |>
    dplyr::arrange(arm_index)

  list(
    # Outcome-effect-size level
    N = nrow(data_indexed),
    y = data_indexed$yi,
    se_y = sqrt(data_indexed$vi),

    # Latent protein intake: one value per nested study–arm
    N_arm = nrow(arm_data),
    arm = data_indexed$arm_index,
    x_obs = arm_data$protein_intake,
    x_se = arm_data$protein_intake_se,

    # Nested meta-analytic random-intercept structure
    N_study = dplyr::n_distinct(data_indexed$study_index),
    study = data_indexed$study_index,

    N_outcome = dplyr::n_distinct(data_indexed$outcome_index),
    outcome = data_indexed$outcome_index,

    # Fixed-effect moderators
    resistance_exercise_code = data_indexed$resistance_exercise_code,
    duration_centre = data_indexed$duration_centre,
    age_centre = data_indexed$age_centre
  )
}

