library(tidyverse)
library(metafor)
library(emmeans)
library(mgcv)
library(patchwork)
library(splines)
library(lspline)
library(rstan)
library(brms)


study_characteristics <- read_csv("nunes_study_characteristics.csv") |>
  janitor::clean_names() |>
  select(id, ref, sex, age, duration, resistance_exercise, int_group, int_protein_intake, con_protein_intake) |>
  pivot_longer(
    8:9,
    names_to = "condition",
    values_to = "protein_intake"
  ) |>
  mutate(
    protein_intake = if_else(protein_intake == 10000, NA, protein_intake),
    condition = case_when(
      condition == "int_protein_intake" ~ "int",
      condition == "con_protein_intake" ~ "con"
    )
  ) |>
  rename(intervention = "int_group")

data <- read_csv("nunes_lean_mass_data.csv") |>
  janitor::clean_names() |>
  mutate(
    # Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
    int_ri = (int_pre_sd^2 + int_post_sd^2 - int_delta_sd^2)/(2 * int_pre_sd * int_post_sd),
    con_ri = (con_pre_sd^2 + con_post_sd^2 - con_delta_sd^2)/(2 * con_pre_sd * con_post_sd)
  ) |>
  mutate(
    # Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
    int_ri = ifelse(between(int_ri,-1,1) == FALSE, NA, int_ri),
    con_ri = ifelse(between(con_ri,-1,1) == FALSE, NA, con_ri)
  ) |>
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

data$int_ri <- ifelse(is.na(data$int_ri), z2r_int, data$int_ri)

data_con <- data |>
  group_by(id) |>
  slice_head()

data_con <- escalc(measure = "ZCOR", ri = con_ri, ni = con_n, data = data_con)

meta_con_ri <- rma.mv(yi, V=vi, data=data_con,
                      random = list(~ 1 | id, ~1 | arm, ~1 | effect), method="REML", test="t",
                      control=list(optimizer="optim", optmethod="Nelder-Mead"))

rob_meta_con_ri <- robust(meta_con_ri, data_con$id)

z2r_con <- psych::fisherz2r(rob_meta_con_ri$b[1])

data$con_ri <- ifelse(is.na(data$con_ri), z2r_con, data$con_ri)

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
  )


data_joined |>
  ggplot(aes(x = protein_intake + 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_smooth(aes(y = yi, group = id), se=FALSE, method = "lm") +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition"
  ) +
  guides(
    size = "none"
  ) +
  facet_grid(. ~ resistance_exercise) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )
  
##### brms models ----

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

# get prior meta-analysis data for set priors

steele_data <- read_csv(url("https://github.com/jamessteeleii/Meta-Analysis-of-Variation-in-Resistance-Training/raw/refs/heads/main/data/Polito%20et%20al.%20RT%20Extracted%20Data.csv"))

# calculate pre-post SDs from SEs
steele_data$RT_pre_sd <- ifelse(is.na(steele_data$RT_pre_se), steele_data$RT_pre_sd, steele_data$RT_pre_se * sqrt(steele_data$RT_n))
steele_data$CON_pre_sd <- ifelse(is.na(steele_data$CON_pre_se), steele_data$CON_pre_sd, steele_data$CON_pre_se * sqrt(steele_data$CON_n))
steele_data$RT_post_sd <- ifelse(is.na(steele_data$RT_post_se), steele_data$RT_post_sd, steele_data$RT_post_se * sqrt(steele_data$RT_n))
steele_data$CON_post_sd <- ifelse(is.na(steele_data$CON_post_se), steele_data$CON_post_sd, steele_data$CON_post_se * sqrt(steele_data$CON_n))

# convert p to t (Change scores)
steele_data$RT_delta_t_value <- replmiss(steele_data$RT_delta_t_value, with(steele_data, qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE)))
steele_data$CON_delta_t_value <- replmiss(steele_data$CON_delta_t_value, with(steele_data, qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE)))

# convert t to SE (Change scores)
steele_data$RT_delta_se <- replmiss(steele_data$RT_delta_se, with(steele_data, ifelse(is.na(RT_delta_m), 
                                                                 (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value)))
steele_data$CON_delta_se <- replmiss(steele_data$CON_delta_se, with(steele_data, ifelse(is.na(CON_delta_m), 
                                                                   (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value)))

# make positive
steele_data$RT_delta_se <- ifelse(steele_data$RT_delta_se < 0, steele_data$RT_delta_se * -1, steele_data$RT_delta_se)
steele_data$CON_delta_se <- ifelse(steele_data$CON_delta_se < 0, steele_data$CON_delta_se * -1, steele_data$CON_delta_se)

# convert CI to SE (Change scores)
steele_data$RT_delta_se <- replmiss(steele_data$RT_delta_se, with(steele_data, (RT_delta_CI_upper - RT_delta_CI_lower)/3.92))
steele_data$CON_delta_se <- replmiss(steele_data$CON_delta_se, with(steele_data, (CON_delta_CI_upper - CON_delta_CI_lower)/3.92))

# convert SE to SD (Change scores)
steele_data$RT_delta_sd <- replmiss(steele_data$RT_delta_sd, with(steele_data, RT_delta_se * sqrt(RT_n)))
steele_data$CON_delta_sd <- replmiss(steele_data$CON_delta_sd, with(steele_data, CON_delta_se * sqrt(CON_n)))

# calculate pre-post correlation coefficient for those with pre, post, and delta SDs
steele_data$RT_ri <- (steele_data$RT_pre_sd^2 + steele_data$RT_post_sd^2 - steele_data$RT_delta_sd^2)/(2 * steele_data$RT_pre_sd * steele_data$RT_post_sd)
steele_data$CON_ri <- (steele_data$CON_pre_sd^2 + steele_data$CON_post_sd^2 - steele_data$CON_delta_sd^2)/(2 * steele_data$CON_pre_sd * steele_data$CON_post_sd)

# impute median within group pre-post correlations assumptions to studies
# we'll filter to > -1 < 1 as there are some odd values likely due to misreporting, miscalcuations in original studies, or dropouts
ri <- c(steele_data$RT_ri, steele_data$CON_ri)
ri <- subset(ri, ri >-1 & ri <1)
steele_data$ri_avg <- as.numeric(strrep(median(ri, na.rm = TRUE), 1))

steele_data_RT <- escalc(measure="SMCR", m1i=RT_post_m, 
                         m2i=RT_pre_m, sd1i=RT_pre_sd, ni=RT_n, ri=ri_avg, data = steele_data) |>
  mutate(condition = 0) |>
  select(outcome, study, arm, condition, yi, vi)

steele_data_CON <- escalc(measure="SMCR", m1i=CON_post_m, 
                         m2i=CON_pre_m, sd1i=CON_pre_sd, ni=CON_n, ri=ri_avg, data = steele_data) |>
  mutate(condition = 1) |>
  select(outcome, study, condition, yi, vi) |>
  distinct() |>
  rowid_to_column("arm") |>
  relocate(outcome, study, arm, condition, yi, vi)


steele_data <- bind_rows(steele_data_RT, steele_data_CON) |>
  filter(outcome == "hypertrophy") |>
  unite("arm", c(outcome,study,arm,condition), sep = "_", remove = FALSE) |>
  mutate(arm = dense_rank(arm)) |>
  rowid_to_column("effect")

# frequentist meta-analysis estimates
steele_meta <- rma.mv(yi, vi,
                           random = list(~ 1 | study, ~1 | arm, ~1 | effect),
                           mods = ~ condition,
                           data = steele_data,
                           method="REML", test="t",
                           # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_steele_meta <- robust(steele_meta, cluster = steele_data$study, clubSandwich = TRUE)

priors <- 
  c(
    # # We set priors based on the overall estimates for both control and RT conditions in an arm based model using Steele et al. data
    set_prior(paste("student_t(3,", rob_steele_meta$b[1],",", rob_steele_meta$se[1],")"), class = "b", coef = "Intercept"),
    set_prior(paste("student_t(3,", rob_steele_meta$b[2],",", rob_steele_meta$se[2],")"), class = "b", coef = "resistance_exercise_code"),
    # All other b parameters are kept as uninformative but weakly regularising
    set_prior("student_t(3,0,10)", class = "b")
    # All other priors for variance parameters are kept as default weakly regularising
  )

  ##### Fit all dose function models ----

##### Intercept only model model

meta_arm_intercept <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                            (1 | id) + 
                            (1 | arm) + 
                            (1 | effect),
                          data = data_joined,
                          prior = priors,
                          chains = 4,
                          cores = 4,
                          seed = 1988,
                          warmup = 2000,
                          iter = 8000,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_intercept)
# 
# plot(meta_arm_intercept)

##### Linear model
meta_arm_linear <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect),
                       data = data_joined,
                       prior = priors,
                       chains = 4,
                       cores = 4,
                       seed = 1988,
                       warmup = 2000,
                       iter = 8000,
                       control = list(adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_linear)
# 
# plot(meta_arm_linear)

##### Linear model
meta_arm_lspline <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect),
                       data = data_joined,
                       prior = priors,
                       chains = 4,
                       cores = 4,
                       seed = 1988,
                       warmup = 2000,
                       iter = 8000,
                       control = list(adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_lspline)
# 
# plot(meta_arm_lspline)

##### Linear and log terms model
meta_arm_log <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                      protein_intake + log1p(protein_intake) +
                          (1 | id) + 
                          (1 | arm) + 
                          (1 | effect),
                        data = data_joined,
                        prior = priors,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        warmup = 2000,
                        iter = 8000,
                        control = list(adapt_delta = 0.99),
                    save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_log)
# 
# plot(meta_arm_log)

##### TPS model
meta_arm_tps <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                      s(protein_intake, bs="tp", k=4) +
                      (1 | id) + 
                      (1 | arm) + 
                      (1 | effect),
                    data = data_joined,
                    prior = priors,
                    chains = 4,
                    cores = 4,
                    seed = 1988,
                    warmup = 2000,
                    iter = 8000,
                    control = list(adapt_delta = 0.99),
                    save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_tps)
# 
# plot(meta_arm_tps)


##### Fit all dose function models with additional centred moderators ----

##### Intercept only model model

meta_arm_mods_intercept <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                 duration_centre + age_centre +
                            (1 | id) + 
                            (1 | arm) + 
                            (1 | effect),
                          data = data_joined,
                          prior = priors,
                          chains = 4,
                          cores = 4,
                          seed = 1988,
                          warmup = 2000,
                          iter = 8000,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_mods_intercept)
# 
# plot(meta_arm_mods_intercept)

##### Linear model
meta_arm_mods_linear <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake +
                           duration_centre + age_centre +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect),
                       data = data_joined,
                       prior = priors,
                       chains = 4,
                       cores = 4,
                       seed = 1988,
                       warmup = 2000,
                       iter = 8000,
                       control = list(adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_mods_linear)
# 
# plot(meta_arm_mods_linear)

##### Linear model
meta_arm_mods_lspline <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) +
                            duration_centre + age_centre +
                            (1 | id) + 
                          (1 | arm) + 
                          (1 | effect),
                        data = data_joined,
                        prior = priors,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        warmup = 2000,
                        iter = 8000,
                        control = list(adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_mods_lspline)
# 
# plot(meta_arm_mods_lspline)

##### Linear and log terms model
meta_arm_mods_log <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                      protein_intake + log1p(protein_intake) +
                        duration_centre + age_centre +
                        (1 | id) + 
                      (1 | arm) + 
                      (1 | effect),
                    data = data_joined,
                    prior = priors,
                    chains = 4,
                    cores = 4,
                    seed = 1988,
                    warmup = 2000,
                    iter = 8000,
                    control = list(adapt_delta = 0.99),
                    save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_mods_log)
# 
# plot(meta_arm_mods_log)

##### TPS model
meta_arm_mods_tps <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                           s(protein_intake, bs="tp", k=4) +
                           duration_centre + age_centre +
                        (1 | id) + 
                      (1 | arm) + 
                      (1 | effect),
                    data = data_joined,
                    prior = priors,
                    chains = 4,
                    cores = 4,
                    seed = 1988,
                    warmup = 2000,
                    iter = 8000,
                    control = list(adapt_delta = 0.99),
                    save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_mods_tps)
# 
# plot(meta_arm_mods_tps)



##### Fit all dose function models with random slopes ----

##### Linear model
meta_arm_slopes_linear <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake +
                         (protein_intake | id) + 
                         (1 | arm) + 
                         (1 | effect),
                       data = data_joined,
                       prior = priors,
                       chains = 4,
                       cores = 4,
                       seed = 1988,
                       warmup = 2000,
                       iter = 8000,
                       control = list(adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_linear)
# 
# plot(meta_arm_slopes_linear)

##### Linear spline model
meta_arm_slopes_lspline <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) +
                          (lspline(protein_intake, 0.48) | id) + 
                          (1 | arm) + 
                          (1 | effect),
                        data = data_joined,
                        prior = priors,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        warmup = 2000,
                        iter = 8000,
                        control = list(adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_lspline)
# 
# plot(meta_arm_slopes_lspline)

##### Linear and log terms model
meta_arm_slopes_log <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                      protein_intake + log1p(protein_intake) +
                      (protein_intake + log1p(protein_intake) | id) + 
                      (1 | arm) + 
                      (1 | effect),
                    data = data_joined,
                    prior = priors,
                    chains = 4,
                    cores = 4,
                    seed = 1988,
                    warmup = 2000,
                    iter = 8000,
                    control = list(adapt_delta = 0.99),
                    save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_log)
# 
# plot(meta_arm_slopes_log)


##### TPS model
meta_arm_slopes_tps <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                             s(protein_intake, bs="tp", k=4) +
                             s(protein_intake, id, bs="fs") +
                      (1 | arm) + 
                      (1 | effect),
                    data = data_joined,
                    prior = priors,
                    chains = 4,
                    cores = 4,
                    seed = 1988,
                    warmup = 2000,
                    iter = 8000,
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 12),
                    save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_tps)
# 
# plot(meta_arm_slopes_tps)

##### Fit all dose function models with additional centred moderators with random slopes ----

##### Linear model
meta_arm_slopes_mods_linear <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                protein_intake +
                                  duration_centre + age_centre +
                                  (protein_intake | id) + 
                                (1 | arm) + 
                                (1 | effect),
                              data = data_joined,
                              prior = priors,
                              chains = 4,
                              cores = 4,
                              seed = 1988,
                              warmup = 2000,
                              iter = 8000,
                              control = list(adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_mods_linear)
# 
# plot(meta_arm_slopes_mods_linear)

##### Linear spline model
meta_arm_slopes_mods_lspline <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                 lspline(protein_intake, 0.48) +
                                   duration_centre + age_centre +
                                   (lspline(protein_intake, 0.48) | id) + 
                                 (1 | arm) + 
                                 (1 | effect),
                               data = data_joined,
                               prior = priors,
                               chains = 4,
                               cores = 4,
                               seed = 1988,
                               warmup = 2000,
                               iter = 8000,
                               control = list(adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_mods_lspline)
# 
# plot(meta_arm_slopes_mods_lspline)

##### Linear and log terms model
meta_arm_slopes_mods_log <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                             protein_intake + log1p(protein_intake) +
                               duration_centre + age_centre +
                               (protein_intake + log1p(protein_intake) | id) + 
                             (1 | arm) + 
                             (1 | effect),
                           data = data_joined,
                           prior = priors,
                           chains = 4,
                           cores = 4,
                           seed = 1988,
                           warmup = 2000,
                           iter = 8000,
                           control = list(adapt_delta = 0.99),
                           save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_mods_log)
# 
# plot(meta_arm_slopes_mods_log)


##### TPS model
meta_arm_slopes_mods_tps <- brm(yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                  s(protein_intake, bs="tp", k=4) +
                                  s(protein_intake, id, bs="fs") +
                               duration_centre + age_centre +
                             (1 | arm) + 
                             (1 | effect),
                           data = data_joined,
                           prior = priors,
                           chains = 4,
                           cores = 4,
                           seed = 1988,
                           warmup = 2000,
                           iter = 8000,
                           control = list(adapt_delta = 0.99,
                                          max_treedepth = 12),
                           save_pars = save_pars(all = TRUE)
)

# pp_check(meta_arm_slopes_mods_tps)
# 
# plot(meta_arm_slopes_mods_tps)
















##### Compare models ----

library(bayestestR)

BF_mean_models <- bayesfactor_models(
  meta_arm_intercept,
  meta_arm_linear,
  meta_arm_lspline,
  meta_arm_log,
  meta_arm_tps,
  meta_arm_mods_intercept,
  meta_arm_mods_linear,
  meta_arm_mods_lspline,
  meta_arm_mods_log,
  meta_arm_mods_tps,
  meta_arm_slopes_linear,
  meta_arm_slopes_lspline,
  meta_arm_slopes_log,
  meta_arm_slopes_tps,
  meta_arm_slopes_mods_linear,
  meta_arm_slopes_mods_lspline,
  meta_arm_slopes_mods_log,
  meta_arm_slopes_mods_tps
)

BF_2log <- function(x) (2*x)

BF_mean_models <- as_tibble(as.matrix(BF_mean_models))  |>
  mutate_at(1:18, BF_2log) |>
  rowid_to_column("Denominator") |>
  mutate(Denominator = case_when(
    Denominator == 1 ~ "Intercept only model",
    Denominator == 2 ~ "Linear model",
    Denominator == 3 ~ "Linear spline (1.6g/kg/day) model",
    Denominator == 4 ~ "Linear and log term model",
    Denominator == 5 ~  "Thin plate spline model",
    Denominator == 6 ~ "Intercept only model (+ moderators)",
    Denominator == 7 ~ "Linear model (+ moderators)",
    Denominator == 8 ~ "Linear spline (1.6g/kg/day) model (+ moderators)",
    Denominator == 9 ~ "Linear and log term model (+ moderators)",
    Denominator == 10 ~  "Thin plate spline model (+ moderators)",
    Denominator == 11 ~ "Linear model (+ random slopes)",
    Denominator == 12 ~ "Linear spline (1.6g/kg/day) model (+ random slopes)",
    Denominator == 13 ~ "Linear and log term model (+ random slopes)",
    Denominator == 14 ~  "Thin plate spline model (+ random smooths)",
    Denominator == 15 ~ "Linear model (+ moderators & random slopes)",
    Denominator == 16 ~ "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
    Denominator == 17 ~ "Linear and log term model (+ moderators & random slopes)",
    Denominator == 18 ~  "Thin plate spline model (+ moderators & random smooths)"
    
  )) |>
  rename("Intercept only model" = 2,
         "Linear model" = 3,
         "Linear spline (1.6g/kg/day) model" = 4,
         "Linear and log term model" = 5,
         "Thin plate spline model" = 6,
         "Intercept only model (+ moderators)" = 7,
         "Linear model (+ moderators)" = 8,
         "Linear spline (1.6g/kg/day) model (+ moderators)" = 9,
         "Linear and log term model (+ moderators)" = 10,
         "Thin plate spline model (+ moderators)" = 11,
         "Linear model (+ random slopes)" = 12,
         "Linear spline (1.6g/kg/day) model (+ random slopes)" = 13,
         "Linear and log term model (+ random slopes)" = 14,
         "Thin plate spline model (+ random smooths)" = 15,
         "Linear model (+ moderators & random slopes)" = 16,
         "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)" = 17,
         "Linear and log term model (+ moderators & random slopes)" = 18,
         "Thin plate spline model (+ moderators & random smooths)" = 19) |>
  pivot_longer(2:19, names_to = "Numerator", values_to = "logBF")

model_comparisons <- BF_mean_models |> 
  mutate(Denominator = factor(Denominator, levels= c(
    "Intercept only model",
    "Linear model",
    "Linear spline (1.6g/kg/day) model",
    "Linear and log term model",
    "Thin plate spline model",
    "Intercept only model (+ moderators)",
    "Linear model (+ moderators)",
    "Linear spline (1.6g/kg/day) model (+ moderators)" ,
    "Linear and log term model (+ moderators)",
    "Thin plate spline model (+ moderators)",
    "Linear model (+ random slopes)",
    "Linear spline (1.6g/kg/day) model (+ random slopes)",
    "Linear and log term model (+ random slopes)",
    "Thin plate spline model (+ random smooths)",
    "Linear model (+ moderators & random slopes)",
    "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
    "Linear and log term model (+ moderators & random slopes)",
    "Thin plate spline model (+ moderators & random smooths)")),
    Numerator = factor(Numerator, levels= c( 
      "Intercept only model",
      "Linear model",
      "Linear spline (1.6g/kg/day) model",
      "Linear and log term model",
      "Thin plate spline model",
      "Intercept only model (+ moderators)",
      "Linear model (+ moderators)",
      "Linear spline (1.6g/kg/day) model (+ moderators)" ,
      "Linear and log term model (+ moderators)",
      "Thin plate spline model (+ moderators)",
      "Linear model (+ random slopes)",
      "Linear spline (1.6g/kg/day) model (+ random slopes)",
      "Linear and log term model (+ random slopes)",
      "Thin plate spline model (+ random smooths)",
      "Linear model (+ moderators & random slopes)",
      "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
      "Linear and log term model (+ moderators & random slopes)",
      "Thin plate spline model (+ moderators & random smooths)")),
    logBF = as.numeric(logBF)) |>
  ggplot(aes(x=Numerator, y=Denominator, fill=logBF)) +
  geom_tile() +
  geom_raster() +
  geom_text(aes(label = round(logBF,2))) +
  scale_fill_gradient2(low = "#E69F00", mid="white", high = "#56B4E9") +
  scale_y_discrete(limits=rev, labels = function(x) str_wrap(x, width = 25)) +
  scale_x_discrete(position = "top", labels = function(x) str_wrap(x, width = 25)) +
  labs(title = "Comparing models using 2×log(BF)",
       fill = "2×log(BF)",
       caption = "Kass and Raferty (1995) scale:
       -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2.5),
        axis.text = element_text(size = 6))

ggsave(plot = model_comparisons, filename = "plots/model_comparisons.tiff", device = "tiff", dpi = 300, w = 25, h = 12.5)



library(tidybayes)
library(marginaleffects)

preds <- crossing(
  resistance_exercise_code = c(0,1),
  protein_intake = seq(-0.32, 3.28, by = 0.01),
  vi = 0
) |>
  add_epred_draws(meta_arm_tps, re_formula = NA) |>
  mutate(
    protein_intake = protein_intake + 1.12
  )


summary <- preds |>
  group_by(resistance_exercise_code, protein_intake) |>
  mean_qi() |>
  mutate(
    resistance_exercise_label = "Resistance Exercise?",
    resistance_exercise = case_when(
      resistance_exercise_code == 0 ~ "YES",
      resistance_exercise_code == 1 ~ "NO"
    )
  )

preds_plot <- summary |>
  ggplot(aes(x = protein_intake, y = .epred)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
  geom_ribbon(aes(ymin = steele_meta$ci.lb[1], ymax = steele_meta$ci.ub[1]), 
             color = NA, fill = "darkred", alpha = 0.25) +
  geom_ribbon(aes(ymin = steele_meta$ci.lb[1] + steele_meta$ci.lb[2], ymax = steele_meta$ci.ub[1] + steele_meta$ci.ub[2]), 
              color = NA, fill = "#56B4E9", alpha = 0.25) +
  geom_point(data = data_joined |>
               mutate(resistance_exercise_label = "Resistance Exercise?"), 
             aes(x = protein_intake + 1.12, y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = .epred.lower, ymax = .epred.upper),
              alpha = 0.5, color = "black", size = 0.25) +
  geom_line(aes(y = .epred), size = 0.75, color = "black") +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  ggh4x::facet_nested(.~ resistance_exercise_label + resistance_exercise) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass/Muscle Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass/muscle change",
    subtitle = "Note, red band indicates prior interval for resistance exercise at habitual protein intake, and blue band indicates non-training controls"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slopes <- slopes(
  meta_arm_tps,
  variables = "protein_intake",
  newdata = datagrid(
    protein_intake = seq(-0.32, 3.28, by = 0.01),
    vi = 0
  ),
  re_formula = NA
) |>
  mutate(protein_intake = protein_intake + 1.12)


slopes_plot <- slopes |>
  ggplot(aes(x = protein_intake, y = estimate*0.42)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
  geom_hline(yintercept = c(-0.1,0.1), linetype = "dashed", alpha = 0.75) +
  geom_ribbon(aes(ymin = conf.low*0.42, ymax = conf.high*0.42),
              alpha = 0.5, color = "black", size = 0.25) +
  geom_line(size = 1, color = "black") +
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
    axis.title = element_text(size = 10)
  )

tps_plots <- (preds_plot / slopes_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based thin plate spline meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study, arm, and effect"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = tps_plots, filename = "plots/tps_plots.tiff", device = "tiff", dpi = 300, w = 10, h = 10)






##### Best fit model ----
# preds <- predict(rob_meta_arm_tps, newmods=PredictMat(sm, data.frame(protein_intake=seq(-0.7,2.9,0.1))), addx = TRUE)

# #### g marginal slope
# rob_meta_arm_tps_slope_prep <- emmprep(rob_meta_arm_tps,
#                                        at=PredictMat(sm, data.frame(protein_intake=seq(-0.7,2.9,0.1))))
# 
# rob_meta_arm_tps_slope_prep <- confint(emmeans(rob_meta_arm_tps_slope_prep,consec~protein_intake,
#                                                weights = "prop")$contrasts)
# 
# rob_meta_arm_tps_slope_prep2<-qt(.975,rob_meta_arm_tps_slope_prep$df)
# rob_meta_arm_tps_slope_prep3<-sqrt((rob_meta_arm_tps_slope_prep$SE^2)+sum(rob_meta_arm_tps$sigma2))
# 
# rob_meta_arm_tps_slope_prep$lower.PI<-rob_meta_arm_tps_slope_prep$estimate-(rob_meta_arm_tps_slope_prep2*rob_meta_arm_tps_slope_prep3)
# rob_meta_arm_tps_slope_prep$upper.PI<-rob_meta_arm_tps_slope_prep$estimate+(rob_meta_arm_tps_slope_prep2*rob_meta_arm_tps_slope_prep3)
# 


preds <- predict(rob_meta_arm_slopes_linear, newmods=cbind(seq(-0.7,2.9,0.1)), addx = TRUE)


#### g marginal slope
rob_meta_arm_slopes_linear_slope_prep <- emmprep(rob_meta_arm_slopes_linear,
                                                 at=list(protein_intake=seq(-0.7,2.9,0.1)))

rob_meta_arm_slopes_linear_slope_prep <- confint(emmeans(rob_meta_arm_slopes_linear_slope_prep,consec~protein_intake,
                                                         weights = "prop")$contrasts)

rob_meta_arm_slopes_linear_slope_prep2<-qt(.975,rob_meta_arm_slopes_linear_slope_prep$df)
rob_meta_arm_slopes_linear_slope_prep3<-sqrt((rob_meta_arm_slopes_linear_slope_prep$SE^2)+sum(rob_meta_arm_slopes_linear$sigma2))

rob_meta_arm_slopes_linear_slope_prep$lower.PI<-rob_meta_arm_slopes_linear_slope_prep$estimate-(rob_meta_arm_slopes_linear_slope_prep2*rob_meta_arm_slopes_linear_slope_prep3)
rob_meta_arm_slopes_linear_slope_prep$upper.PI<-rob_meta_arm_slopes_linear_slope_prep$estimate+(rob_meta_arm_slopes_linear_slope_prep2*rob_meta_arm_slopes_linear_slope_prep3)





pred_plot <- preds |>
  as.data.frame() |>
  mutate(protein_intake = seq(-0.7,2.9,0.1)) |>
  ggplot(aes(x = protein_intake + 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data_joined, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass/Muscle Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass/muscle change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )

slope_plot <- rob_meta_arm_slopes_linear_slope_prep |>
  mutate(protein_intake = seq(-0.6,2.9,0.1)) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
  geom_line(aes(y=estimate)) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Marginal Slope for Lean Mass Change (SMD)",
    # y = "Marginal Slope for Lean Mass Change (SMD)",
    title = "Incremental marginal slope of lean mass change for 0.1 g/kg/day increase"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )

plots <- (pred_plot / slope_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based linear meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect"
  ) &
  theme(
    legend.position = "bottom"
  )

# slope_plot <- preds |>
#   as.data.frame()  |>
#   mutate(protein_intake = seq(0,3.6,0.1),
#          slope = c(NA, diff(pred) / diff(protein_intake)),
#          lower.CL = c(NA, diff(ci.lb) / diff(protein_intake)),
#          upper.CL = c(NA, diff(ci.ub) / diff(protein_intake))) |>
#   ggplot(aes(x = protein_intake + 0.8)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
#   geom_line(aes(y=slope)) +
#   labs(
#     x = "Protein Intake (g/kg/day)",
#     y = "Marginal Slope for Lean Mass Change (SMD)",
#     # y = "Marginal Slope for Lean Mass Change (SMD)",
#     title = "Incremental marginal slope of lean mass change for 0.1 g/kg/day increase"
#   ) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(size = 12),
#     axis.title = element_text(size = 10)
#   )
# 
# 
# tps_plots <- (pred_tps_plot / slope_tps_plot) +
#   plot_annotation(
#     tag_levels = "i",
#     title = "Arm based thin plate spline meta-regression of protein intake on lean mass change",
#     subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
#   ) &
#   theme(
#     legend.position = "bottom"
#   )

ggsave(plot = pred_plot, filename = "plots/linear_slopes_plot.tiff", device = "tiff", dpi = 300, w = 7.5, h = 3.75)





##### Old metafor models ----

##### Fit all dose function models ----

##### Intercept only model model
meta_arm_intercept <- rma.mv(yi, vi,
                             random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                             # struct = "GEN",
                             data = data_joined,
                             method="REML", test="t",
                             # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_intercept <- robust(meta_arm_intercept, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_linear <- rma.mv(yi, vi,
                          random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                          mods = ~ protein_intake,
                          # struct = "GEN",
                          data = data_joined,
                          method="REML", test="t",
                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_linear <- robust(meta_arm_linear, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_lspline <- rma.mv(yi, vi,
                           random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                           mods = ~ lspline(protein_intake, 0.215),
                           # struct = "GEN",,
                           # struct = "GEN",
                           data = data_joined,
                           method="REML", test="t",
                           # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_lspline <- robust(meta_arm_lspline, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear and log terms model
meta_arm_log <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                       mods = ~ protein_intake + log1p(protein_intake),
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_log <- robust(meta_arm_log, cluster = data_joined$id, clubSandwich = TRUE)

##### NCS model
meta_arm_ncs <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                       mods = ~ ns(protein_intake, knots=3),
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_ncs <- robust(meta_arm_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=4), data=data_joined, absorb.cons=TRUE)[[1]]

meta_arm_tps <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                       mods = ~ sm$X,
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_tps <- robust(meta_arm_tps, cluster = data_joined$id, clubSandwich = TRUE)

##### Fit all dose function models with additional centred moderators ----

##### Intercept only model model
meta_arm_mods_intercept <- rma.mv(yi, vi,
                             random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                             mods = ~ duration_centre + resistance_exercise_code + age_centre,
                             # struct = "GEN",
                             data = data_joined,
                             method="REML", test="t",
                             # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_mods_intercept <- robust(meta_arm_mods_intercept, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_mods_linear <- rma.mv(yi, vi,
                             random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                             mods = ~ protein_intake +
                              duration_centre + resistance_exercise_code + age_centre,
                             # struct = "GEN",
                             data = data_joined,
                             method="REML", test="t",
                             # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_mods_linear <- robust(meta_arm_mods_linear, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_mods_lspline <- rma.mv(yi, vi,
                           random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                           mods = ~ lspline(protein_intake, 0.215) + 
                             duration_centre + resistance_exercise_code + age_centre,
                           # struct = "GEN",,
                           # struct = "GEN",
                           data = data_joined,
                           method="REML", test="t",
                           # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_mods_lspline <- robust(meta_arm_mods_lspline, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear and log terms model
meta_arm_mods_log <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                       mods = ~ protein_intake + log1p(protein_intake) +
                         duration_centre + resistance_exercise_code + age_centre,
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_mods_log <- robust(meta_arm_mods_log, cluster = data_joined$id, clubSandwich = TRUE)

##### NCS model
meta_arm_mods_ncs <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                       mods = ~ ns(protein_intake, knots=3) +
                         duration_centre + resistance_exercise_code + age_centre,
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_mods_ncs <- robust(meta_arm_mods_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=4), data=data_joined, absorb.cons=TRUE)[[1]]

meta_arm_mods_tps <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                       mods = ~ sm$X +
                         duration_centre + resistance_exercise_code + age_centre,
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_mods_tps <- robust(meta_arm_mods_tps, cluster = data_joined$id, clubSandwich = TRUE)




##### Fit all dose function models with random slopes ----

##### Linear model
meta_arm_slopes_linear <- rma.mv(yi, vi,
                          random = list(~ ~ protein_intake | id, ~1 | effect),
                          mods = ~ protein_intake,
                          struct = "GEN",
                          data = data_joined,
                          method="REML", test="t",
                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_linear <- robust(meta_arm_slopes_linear, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_slopes_lspline <- rma.mv(yi, vi,
                           random = list(~ lspline(protein_intake, 0.215) | id, ~1 | effect),
                           mods = ~ lspline(protein_intake, 0.215),
                           struct = "GEN",
                           data = data_joined,
                           method="REML", test="t",
                           # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_lspline <- robust(meta_arm_slopes_lspline, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear and log terms model
meta_arm_slopes_log <- rma.mv(yi, vi,
                       random = list(~ protein_intake + log1p(protein_intake) | id, ~1 | effect),
                       mods = ~ protein_intake + log1p(protein_intake),
                       struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_log <- robust(meta_arm_slopes_log, cluster = data_joined$id, clubSandwich = TRUE)

##### NCS model
meta_arm_slopes_ncs <- rma.mv(yi, vi,
                       random = list(~ ns(protein_intake, knots=3) | id, ~1 | effect),
                       mods = ~ ns(protein_intake, knots=3),
                       struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_slopes_ncs <- robust(meta_arm_slopes_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=4), data=data_joined, absorb.cons=TRUE)[[1]]

data_joined <- bind_cols(data_joined, sm$X)

meta_arm_slopes_tps <- rma.mv(yi, vi,
                       random = list(~ `...34` + `...35` | id, ~1 | effect),
                       mods = ~ sm$X,
                       struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_tps <- robust(meta_arm_slopes_tps, cluster = data_joined$id, clubSandwich = TRUE)

##### Fit all dose function models with additional centred moderators with random slopes ----

##### Linear model
meta_arm_slopes_mods_linear <- rma.mv(yi, vi,
                               random = list(~ protein_intake | id, ~1 | effect),
                               mods = ~ protein_intake +
                                 duration_centre + resistance_exercise_code + age_centre,
                               struct = "GEN",
                               data = data_joined,
                               method="REML", test="t",
                               # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_mods_linear <- robust(meta_arm_slopes_mods_linear, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_slopes_mods_lspline <- rma.mv(yi, vi,
                                random = list(~ lspline(protein_intake, 0.215) | id, ~1 | effect),
                                mods = ~ lspline(protein_intake, 0.215) + 
                                  duration_centre + resistance_exercise_code + age_centre,
                                struct = "GEN",
                                data = data_joined,
                                method="REML", test="t",
                                # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_mods_lspline <- robust(meta_arm_slopes_mods_lspline, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear and log terms model
meta_arm_slopes_mods_log <- rma.mv(yi, vi,
                            random = list(~ protein_intake + log1p(protein_intake) | id, ~1 | effect),
                            mods = ~ protein_intake + log1p(protein_intake) +
                              duration_centre + resistance_exercise_code + age_centre,
                            struct = "GEN",
                            data = data_joined,
                            method="REML", test="t",
                            # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_mods_log <- robust(meta_arm_slopes_mods_log, cluster = data_joined$id, clubSandwich = TRUE)

##### NCS model
meta_arm_slopes_mods_ncs <- rma.mv(yi, vi,
                            random = list(~ ns(protein_intake, knots=3) | id, ~1 | effect),
                            mods = ~ ns(protein_intake, knots=3) +
                              duration_centre + resistance_exercise_code + age_centre,
                            struct = "GEN",
                            data = data_joined,
                            method="REML", test="t",
                            # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_slopes_mods_ncs <- robust(meta_arm_slopes_mods_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=4), data=data_joined, absorb.cons=TRUE)[[1]]

meta_arm_slopes_mods_tps <- rma.mv(yi, vi,
                            random = list(~ `...34` + `...35`  | id, ~1 | effect),
                            mods = ~ sm$X +
                              duration_centre + resistance_exercise_code + age_centre,
                            struct = "GEN",
                            data = data_joined,
                            method="REML", test="t",
                            # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_slopes_mods_tps <- robust(meta_arm_slopes_mods_tps, cluster = data_joined$id, clubSandwich = TRUE)


















##### Compare models ----

library(bayestestR)

BF_mean_models <- bayesfactor_models(rob_meta_arm_intercept,
                                     rob_meta_arm_linear,
                                     rob_meta_arm_lspline,
                                     rob_meta_arm_log,
                                     rob_meta_arm_ncs,
                                     rob_meta_arm_tps,
                                     rob_meta_arm_mods_intercept,
                                     rob_meta_arm_mods_linear,
                                     rob_meta_arm_mods_lspline,
                                     rob_meta_arm_mods_log,
                                     rob_meta_arm_mods_ncs,
                                     rob_meta_arm_mods_tps,
                                     rob_meta_arm_slopes_linear,
                                     rob_meta_arm_slopes_lspline,
                                     rob_meta_arm_slopes_log,
                                     rob_meta_arm_slopes_ncs,
                                     rob_meta_arm_slopes_tps,
                                     rob_meta_arm_slopes_mods_linear,
                                     rob_meta_arm_slopes_mods_lspline,
                                     rob_meta_arm_slopes_mods_log,
                                     rob_meta_arm_slopes_mods_ncs,
                                     rob_meta_arm_slopes_mods_tps)
BF_2log <- function(x) (2*x)

BF_mean_models <- as_tibble(as.matrix(BF_mean_models))  |>
  mutate_at(1:22, BF_2log) |>
  rowid_to_column("Denominator") |>
  mutate(Denominator = case_when(
    Denominator == 1 ~ "Intercept only model",
    Denominator == 2 ~ "Linear model",
    Denominator == 3 ~ "Linear spline (1.6g/kg/day) model",
    Denominator == 4 ~ "Linear and log term model",
    Denominator == 5 ~  "Natural cubic spline model",
    Denominator == 6 ~  "Thin plate spline model",
    Denominator == 7 ~ "Intercept only model (+ moderators)",
    Denominator == 8 ~ "Linear model (+ moderators)",
    Denominator == 9 ~ "Linear spline (1.6g/kg/day) model (+ moderators)",
    Denominator == 10 ~ "Linear and log term model (+ moderators)",
    Denominator == 11 ~  "Natural cubic spline model (+ moderators)",
    Denominator == 12 ~  "Thin plate spline model (+ moderators)",
    Denominator == 13 ~ "Linear model (+ random slopes)",
    Denominator == 14 ~ "Linear spline (1.6g/kg/day) model (+ random slopes)",
    Denominator == 15 ~ "Linear and log term model (+ random slopes)",
    Denominator == 16 ~  "Natural cubic spline model (+ random slopes)",
    Denominator == 17 ~  "Thin plate spline model (+ random slopes)",
    Denominator == 18 ~ "Linear model (+ moderators & random slopes)",
    Denominator == 19 ~ "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
    Denominator == 20 ~ "Linear and log term model (+ moderators & random slopes)",
    Denominator == 21 ~  "Natural cubic spline model (+ moderators & random slopes)",
    Denominator == 22 ~  "Thin plate spline model (+ moderators & random slopes)"
    
  )) |>
  rename("Intercept only model" = 2,
         "Linear model" = 3,
         "Linear spline (1.6g/kg/day) model" = 4,
         "Linear and log term model" = 5,
         "Natural cubic spline model" = 6,
         "Thin plate spline model" = 7,
         "Intercept only model (+ moderators)" = 8,
         "Linear model (+ moderators)" = 9,
         "Linear spline (1.6g/kg/day) model (+ moderators)" = 10,
         "Linear and log term model (+ moderators)" = 11,
         "Natural cubic spline model (+ moderators)" = 12,
         "Thin plate spline model (+ moderators)" = 13,
         "Linear model (+ random slopes)" = 14,
         "Linear spline (1.6g/kg/day) model (+ random slopes)" = 15,
         "Linear and log term model (+ random slopes)" = 16,
         "Natural cubic spline model (+ random slopes)" = 17,
         "Thin plate spline model (+ random slopes)" = 18,
         "Linear model (+ moderators & random slopes)" = 19,
         "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)" = 20,
         "Linear and log term model (+ moderators & random slopes)" = 21,
         "Natural cubic spline model (+ moderators & random slopes)" = 22,
         "Thin plate spline model (+ moderators & random slopes)" = 23) |>
  pivot_longer(2:23, names_to = "Numerator", values_to = "logBF")

model_comparisons <- BF_mean_models |> 
  mutate(Denominator = factor(Denominator, levels= c(
    "Intercept only model",
    "Linear model",
    "Linear spline (1.6g/kg/day) model",
    "Linear and log term model",
    "Natural cubic spline model",
    "Thin plate spline model",
    "Intercept only model (+ moderators)",
    "Linear model (+ moderators)",
    "Linear spline (1.6g/kg/day) model (+ moderators)" ,
    "Linear and log term model (+ moderators)",
    "Natural cubic spline model (+ moderators)",
    "Thin plate spline model (+ moderators)",
    "Linear model (+ random slopes)",
    "Linear spline (1.6g/kg/day) model (+ random slopes)",
    "Linear and log term model (+ random slopes)",
    "Natural cubic spline model (+ random slopes)",
    "Thin plate spline model (+ random slopes)",
    "Linear model (+ moderators & random slopes)",
    "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
    "Linear and log term model (+ moderators & random slopes)",
    "Natural cubic spline model (+ moderators & random slopes)",
    "Thin plate spline model (+ moderators & random slopes)")),
    Numerator = factor(Numerator, levels= c( 
      "Intercept only model",
      "Linear model",
      "Linear spline (1.6g/kg/day) model",
      "Linear and log term model",
      "Natural cubic spline model",
      "Thin plate spline model",
      "Intercept only model (+ moderators)",
      "Linear model (+ moderators)",
      "Linear spline (1.6g/kg/day) model (+ moderators)" ,
      "Linear and log term model (+ moderators)",
      "Natural cubic spline model (+ moderators)",
      "Thin plate spline model (+ moderators)",
      "Linear model (+ random slopes)",
      "Linear spline (1.6g/kg/day) model (+ random slopes)",
      "Linear and log term model (+ random slopes)",
      "Natural cubic spline model (+ random slopes)",
      "Thin plate spline model (+ random slopes)",
      "Linear model (+ moderators & random slopes)",
      "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
      "Linear and log term model (+ moderators & random slopes)",
      "Natural cubic spline model (+ moderators & random slopes)",
      "Thin plate spline model (+ moderators & random slopes)")),
    logBF = as.numeric(logBF)) |>
  ggplot(aes(x=Numerator, y=Denominator, fill=logBF)) +
  geom_tile() +
  geom_raster() +
  geom_text(aes(label = round(logBF,2))) +
  scale_fill_gradient2(low = "#E69F00", mid="white", high = "#56B4E9") +
  scale_y_discrete(limits=rev, labels = function(x) str_wrap(x, width = 25)) +
  scale_x_discrete(position = "top", labels = function(x) str_wrap(x, width = 25)) +
  labs(title = "Comparing models using 2×log(BF)",
       fill = "2×log(BF)",
       caption = "Kass and Raferty (1995) scale:
       -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2.5),
        axis.text = element_text(size = 6))



ggsave(plot = model_comparisons, filename = "plots/model_comparisons.tiff", device = "tiff", dpi = 300, w = 25, h = 12.5)





##### Best fit model ----
# preds <- predict(rob_meta_arm_tps, newmods=PredictMat(sm, data.frame(protein_intake=seq(-0.7,2.9,0.1))), addx = TRUE)

# #### g marginal slope
# rob_meta_arm_tps_slope_prep <- emmprep(rob_meta_arm_tps,
#                                        at=PredictMat(sm, data.frame(protein_intake=seq(-0.7,2.9,0.1))))
# 
# rob_meta_arm_tps_slope_prep <- confint(emmeans(rob_meta_arm_tps_slope_prep,consec~protein_intake,
#                                                weights = "prop")$contrasts)
# 
# rob_meta_arm_tps_slope_prep2<-qt(.975,rob_meta_arm_tps_slope_prep$df)
# rob_meta_arm_tps_slope_prep3<-sqrt((rob_meta_arm_tps_slope_prep$SE^2)+sum(rob_meta_arm_tps$sigma2))
# 
# rob_meta_arm_tps_slope_prep$lower.PI<-rob_meta_arm_tps_slope_prep$estimate-(rob_meta_arm_tps_slope_prep2*rob_meta_arm_tps_slope_prep3)
# rob_meta_arm_tps_slope_prep$upper.PI<-rob_meta_arm_tps_slope_prep$estimate+(rob_meta_arm_tps_slope_prep2*rob_meta_arm_tps_slope_prep3)
# 


preds <- predict(rob_meta_arm_slopes_linear, newmods=cbind(seq(-0.7,2.9,0.1)), addx = TRUE)


#### g marginal slope
rob_meta_arm_slopes_linear_slope_prep <- emmprep(rob_meta_arm_slopes_linear,
                                          at=list(protein_intake=seq(-0.7,2.9,0.1)))

rob_meta_arm_slopes_linear_slope_prep <- confint(emmeans(rob_meta_arm_slopes_linear_slope_prep,consec~protein_intake,
                                                  weights = "prop")$contrasts)

rob_meta_arm_slopes_linear_slope_prep2<-qt(.975,rob_meta_arm_slopes_linear_slope_prep$df)
rob_meta_arm_slopes_linear_slope_prep3<-sqrt((rob_meta_arm_slopes_linear_slope_prep$SE^2)+sum(rob_meta_arm_slopes_linear$sigma2))

rob_meta_arm_slopes_linear_slope_prep$lower.PI<-rob_meta_arm_slopes_linear_slope_prep$estimate-(rob_meta_arm_slopes_linear_slope_prep2*rob_meta_arm_slopes_linear_slope_prep3)
rob_meta_arm_slopes_linear_slope_prep$upper.PI<-rob_meta_arm_slopes_linear_slope_prep$estimate+(rob_meta_arm_slopes_linear_slope_prep2*rob_meta_arm_slopes_linear_slope_prep3)





pred_plot <- preds |>
  as.data.frame() |>
  mutate(protein_intake = seq(-0.7,2.9,0.1)) |>
  ggplot(aes(x = protein_intake + 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data_joined, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass/Muscle Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass/muscle change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )

slope_plot <- rob_meta_arm_slopes_linear_slope_prep |>
  mutate(protein_intake = seq(-0.6,2.9,0.1)) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
  geom_line(aes(y=estimate)) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Marginal Slope for Lean Mass Change (SMD)",
    # y = "Marginal Slope for Lean Mass Change (SMD)",
    title = "Incremental marginal slope of lean mass change for 0.1 g/kg/day increase"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )

plots <- (pred_plot / slope_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based linear meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect"
  ) &
  theme(
    legend.position = "bottom"
  )

# slope_plot <- preds |>
#   as.data.frame()  |>
#   mutate(protein_intake = seq(0,3.6,0.1),
#          slope = c(NA, diff(pred) / diff(protein_intake)),
#          lower.CL = c(NA, diff(ci.lb) / diff(protein_intake)),
#          upper.CL = c(NA, diff(ci.ub) / diff(protein_intake))) |>
#   ggplot(aes(x = protein_intake + 0.8)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
#   geom_line(aes(y=slope)) +
#   labs(
#     x = "Protein Intake (g/kg/day)",
#     y = "Marginal Slope for Lean Mass Change (SMD)",
#     # y = "Marginal Slope for Lean Mass Change (SMD)",
#     title = "Incremental marginal slope of lean mass change for 0.1 g/kg/day increase"
#   ) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(size = 12),
#     axis.title = element_text(size = 10)
#   )
# 
# 
# tps_plots <- (pred_tps_plot / slope_tps_plot) +
#   plot_annotation(
#     tag_levels = "i",
#     title = "Arm based thin plate spline meta-regression of protein intake on lean mass change",
#     subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
#   ) &
#   theme(
#     legend.position = "bottom"
#   )

ggsave(plot = pred_plot, filename = "plots/linear_slopes_plot.tiff", device = "tiff", dpi = 300, w = 7.5, h = 3.75)


##### Linear model brms
meta_arm_linear <- rma.mv(yi, vi,
                          random = list(~ 1 | id, ~1 | arm, ~1 | effect),
                          mods = ~ protein_intake,
                          # struct = "GEN",
                          data = data_joined,
                          method="REML", test="t",
                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)



meta_arm_linear_brms <- brm(yi|se(sqrt(vi)) ~ 1 + protein_intake + (1 | id) + (1|effect))