# Functions for targets pipeline

# Misc
set_rstan_options <- function() {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)
}

# Preliminary re-analysis with Nunes et al. (https://doi.org/10.1002/jcsm.12922) lean soft tissue/muscle data ----

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
      set_prior("student_t(3,0,10)", class = "b")
      # All other priors for variance parameters are kept as default weakly regularising
    )
  
  priors
  
}

# Determining baseline standard deviations to provide reference values for smallest effect size of interest of 0.1

  ### ADD FUNCTIONS TO CALCULATE META-ANALYTIC ESTIMATES OF STANDARD DEVIATIONS

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
    filter(!is.na(duration))
  
  benito_data <- escalc(
    measure = "SDLN",
    sdi = pre_sd,
    ni = pre_n,
    data = benito_data
  )
  
  benito_meta_sd <- rma.mv(yi, vi,
                           random = list(~ 1 | study, ~ 1 | arm),
                           mods = ~ 0 + outcome,
                           data = benito_data,
                           method="REML", test="t")
  
  benito_sds <- broom::tidy(benito_meta_sd) |>
    mutate(across(c(estimate, std.error), exp))
  
  return(benito_sds)
}

estimate_strength_sd <- function(data) {
  
  open_powerlifting_sd <- data |>
    filter(Equipment == "Raw" & Tested == "Yes" & Division == "Open") |>
    mutate(Best3SquatKg = if_else(Best3SquatKg > 0, Best3SquatKg, NA), # We also remove any failed lifts as it impacts the totals
           Best3BenchKg = if_else(Best3BenchKg > 0, Best3BenchKg, NA),
           Best3DeadliftKg = if_else(Best3DeadliftKg > 0, Best3DeadliftKg, NA)) |>
    group_by(Name) |>
    arrange(Date) |>
    slice_head(n=1) |>
    ungroup() |>
    summarise(
      squat_sd = sd(Best3SquatKg, na.rm=TRUE),
      bench_sd = sd(Best3BenchKg, na.rm=TRUE),
      deadlift_sd = sd(Best3DeadliftKg, na.rm=TRUE)
    )
  
  return(open_powerlifting_sd)
}

# Function to fit all candidate models
fit_candidate_models <- function(formula, priors, data) {
  formula <- as.formula(paste0(formula))
  model <- brm(
    formula,
    data = data,
    prior = priors,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    control = list(adapt_delta = 0.99),
    save_pars = save_pars(all = TRUE) 
  )
  
  return(model)
}