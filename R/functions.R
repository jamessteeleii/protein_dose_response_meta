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
  
  
  buckner_sd <- as_tibble(predict(rob_buckner_meta_sd, 
          newmods = log(median(buckner_data$pre_mean)))) |>
    mutate(across(where(is.double), exp))
  
  
  return(buckner_sd)
  
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
    iter = 6000, # just for pre-reg so it doesn't take so long to run
    # iter = 42000, # so as to have 40k iterations for BF calculations
    control = list(adapt_delta = 0.99),
    save_pars = save_pars(all = TRUE) 
  )
  
  return(model)
}

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
      title = "Dose-response effects of protein intake on lean mass change",
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
  
  contrasts_join <-
    contrast_grid |>
    left_join(
      contrasts,
      by = c(
        "study_no",
        "arm_no",
        "outcome_group",
        "measurement_method",
        "unit_of_measurement",
        "time1",
        "time2"
      )
    )
  
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
    
    mutate(
      outcome_no = consecutive_id(
        study_no,
        arm_no,
        outcome_group,
        measurement_method,
        unit_of_measurement
      )
    ) |>
    
    # Sample size as smallest of baseline and followup
    mutate(n = pmin(baseline_n, followup_n)) |>
    
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
                         random = list(~ 1 | study_no, ~ 1 | arm_no, ~ 1 | outcome_no), method="REML", test="t",
                         # mods = ~ 0 + outcome_grou,
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
    distinct(
      study_no,
      outcome_no,
      arm_no,
      baseline_n,
      baseline_sd
    ) |>
    group_by(study_no, outcome_no) |>
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
      by = c("study_no", "outcome_no")
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
    rename(protein_intake = "intervention_dietary_protein_intake_g_kg_day",
           protein_intake_se = "intervention_dietary_protein_intake_g_kg_day_se")
  
  # centre  predictors
  contrasts_analysis <- contrasts_analysis |>
  mutate(
    protein_intake = protein_intake - median(subset(contrasts_analysis$protein_intake, condition == "control" | condition == "placebo"),
                                             na.rm=TRUE),
    duration_centre = time2 - 12,
    age_centre = mean_age_years - 25,
    resistance_exercise_code = case_when(
      resistance_exercise_yes_no == "yes" ~ 0,
      resistance_exercise_yes_no == "no" ~ 1,
    )
  ) |>
    
    group_by(study_no, arm_no, outcome_group, measurement_method) |>
    mutate(
      outcome_number = cur_group_id()
    )
  
  return(contrasts_analysis)

}

## Fitting models using stan ----

make_stan_data <- function(data) {
  arm_data <- data |>
    dplyr::distinct(
      arm_no,
      protein_intake,
      protein_intake_se
    ) |>
    dplyr::arrange(arm_no) |>
    dplyr::mutate(arm_index = dplyr::row_number())
  
  data <- data |>
    dplyr::left_join(
      dplyr::select(arm_data, arm_no, arm_index),
      by = "arm_no"
    )
  
  list(
    N = nrow(data),
    y = data$yi,
    se_y = sqrt(data$vi),
    
    A = nrow(arm_data),
    arm = data$arm_index,
    x_obs = arm_data$protein_intake,
    x_se = arm_data$protein_intake_se,
    
    J_study = dplyr::n_distinct(data$study_no),
    study = match(data$id, unique(data$study_no)),
    
    J_arm = dplyr::n_distinct(data$arm_no),
    meta_arm = data$arm_index,
    
    J_outcome = dplyr::n_distinct(data$outcome_no),
    outcome = match(data$outcome_no, unique(data$outcome_no)),
    
    resistance_exercise_code = data$resistance_exercise_code,
    age_centre = data$age_centre,
    duration_centre = data$duration_centre
  )
}

      # # Functions to generate weights and obtain model averaged preds/slopes
      # get_model_weights_new <- function(models) {
      #   
      #   models_new <- models
      #   
      #   # pass each to object as loo doesn't like the list for some reason
      #   m_intercept <- models_new$`Intercept only model`
      #   m_linear <- models_new$`Linear model`
      #   m_lspline <- models_new$`Linear spline (1.6g/kg/day) model`
      #   m_tps <- models_new$`Thin plate spline model`
      #   m_intercept_mods <- models_new$`Intercept only model (+ moderators)`
      #   m_linear_mods <- models_new$`Linear model (+ moderators)`
      #   m_lspline_mods <- models_new$`Linear spline (1.6g/kg/day) model (+ moderators)`
      #   m_tps_mods <- models_new$ `Thin plate spline model (+ moderators)`
      #   m_linear_ranef <- models_new$`Linear model (+ random slopes)`
      #   m_lspline_ranef <- models_new$`Linear spline (1.6g/kg/day) model (+ random slopes)`
      #   m_tps_ranef <- models_new$ `Thin plate spline model (+ random smooths)`
      #   m_linear_mods_ranef <- models_new$`Linear model (+ moderators & random slopes)`
      #   m_lspline_mods_ranef <- models_new$`Linear spline (1.6g/kg/day) model (+ moderators & random slopes)`
      #   m_tps_mods_ranef <- models_new$`Thin plate spline model (+ moderators & random smooths)`
      #   m_linear_interact <- models_new$`Linear model (+ RT interaction)`
      #   m_lspline_interact <- models_new$`Linear spline (1.6g/kg/day) model (+ RT interaction)`
      #   m_tps_interact <- models_new$`Thin plate spline model (+ RT interaction smooth)`
      #   m_linear_mods_interact <- models_new$`Linear model (+ moderators & RT interaction)`
      #   m_lspline_mods_interact <- models_new$`Linear spline (1.6g/kg/day) model (+ moderators & RT interaction)` 
      #   m_tps_mods_interact <- models_new$`Thin plate spline model (+ moderators & RT interaction)`
      #   m_linear_ranef_interact <- models_new$`Linear model (+ random slopes & RT interaction)`
      #   m_lspline_ranef_interact <- models_new$`Linear spline (1.6g/kg/day) model (+ random slopes & RT interaction)`
      #   m_tps_ranef_interact <- models_new$`Thin plate spline model (+ random smooths & RT interaction smooth)`
      #   m_linear_mods_ranef_interact <- models_new$`Linear model (+ moderators & random slopes & RT interaction)`
      #   m_lspline_mods_ranef_interact <- models_new$`Linear spline (1.6g/kg/day) model (+ moderators & random slopes & RT interaction)`
      #   m_tps_mods_ranef_interact <- models_new$`Thin plate spline model (+ moderators & random smooths & RT interaction smooth)`
      #   
      #   
      #   
      #   weights_stacking <- loo::loo_model_weights(m_intercept, 
      #                                              m_linear, 
      #                                              m_lspline, 
      #                                              m_tps,
      #                                              m_intercept_mods, 
      #                                              m_linear_mods, 
      #                                              m_lspline_mods, 
      #                                              m_tps_mods,
      #                                              m_linear_ranef, 
      #                                              m_lspline_ranef,
      #                                              m_tps_ranef,
      #                                              m_linear_mods_ranef, 
      #                                              m_lspline_mods_ranef, 
      #                                              m_tps_mods_ranef,
      #                                              m_linear_interact, 
      #                                              m_lspline_interact, 
      #                                              m_tps_interact,
      #                                              m_linear_mods_interact, 
      #                                              m_lspline_mods_interact, 
      #                                              m_tps_mods_interact,
      #                                              m_linear_ranef_interact, 
      #                                              m_lspline_ranef_interact,
      #                                              m_tps_ranef_interact,
      #                                              m_linear_mods_ranef_interact, 
      #                                              m_lspline_mods_ranef_interact, 
      #                                              m_tps_mods_ranef_interact,
      #                                              method = "stacking")
      #   
      #   return(weights_stacking)
      # }