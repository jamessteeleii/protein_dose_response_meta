library(tidyverse)
library(metafor)
library(emmeans)
library(mgcv)
library(patchwork)
library(splines)
library(lspline)

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
  rowid_to_column("effect")

# Calculate pre-post correlations for SMCR calculations

data <- escalc(measure = "ZCOR", ri = int_ri, ni = int_n, data = data)

meta_int_ri <- rma.mv(yi, V=vi, data=data,
                      random = list(~ 1 | id, ~1 | effect), method="REML", test="t",
                      control=list(optimizer="optim", optmethod="Nelder-Mead"))

rob_meta_int_ri <- robust(meta_int_ri, data$id)

z2r_int <- psych::fisherz2r(rob_meta_int_ri$b[1])

data$int_ri <- ifelse(is.na(data$int_ri), z2r_int, data$int_ri)

data_con <- data |>
  group_by(id) |>
  slice_head()

data_con <- escalc(measure = "ZCOR", ri = con_ri, ni = con_n, data = data_con)

meta_con_ri <- rma.mv(yi, V=vi, data=data_con,
                      random = list(~ 1 | id, ~1 | effect), method="REML", test="t",
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
  select(-effect) |>
  rowid_to_column("effect") |>
  # add study weights/sizes
  mutate(
    wi = 1/sqrt(vi),
    size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE))) |>
  
  # centre  predictors
  mutate(
    protein_intake = protein_intake - 1.385,
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
  # geom_smooth(aes(y = yi, group = id), se=FALSE, method = "lm") +
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
  # facet_grid(. ~ resistance_exercise) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )
  

##### Fit all dose function models ----

##### Intercept only model model
meta_arm_intercept <- rma.mv(yi, vi,
                             random = list(~ 1 | id, ~1 | effect),
                             # struct = "GEN",
                             data = data_joined,
                             method="REML", test="t",
                             # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_intercept <- robust(meta_arm_intercept, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_linear <- rma.mv(yi, vi,
                          random = list(~ 1 | id, ~1 | effect),
                          mods = ~ protein_intake,
                          # struct = "GEN",
                          data = data_joined,
                          method="REML", test="t",
                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_linear <- robust(meta_arm_linear, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_lspline <- rma.mv(yi, vi,
                           random = list(~ 1 | id, ~1 | effect),
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
                       random = list(~ 1 | id, ~1 | effect),
                       mods = ~ protein_intake + log1p(protein_intake),
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_log <- robust(meta_arm_log, cluster = data_joined$id, clubSandwich = TRUE)

##### NCS model
meta_arm_ncs <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | effect),
                       mods = ~ ns(protein_intake, df=3),
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_ncs <- robust(meta_arm_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=3), data=data_joined, absorb.cons=TRUE)[[1]]

meta_arm_tps <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | effect),
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
                             random = list(~ 1 | id, ~1 | effect),
                             mods = ~ duration_centre + resistance_exercise_code + age_centre,
                             # struct = "GEN",
                             data = data_joined,
                             method="REML", test="t",
                             # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_mods_intercept <- robust(meta_arm_mods_intercept, cluster = data_joined$id, clubSandwich = TRUE)

##### Linear model
meta_arm_mods_linear <- rma.mv(yi, vi,
                             random = list(~ 1 | id, ~1 | effect),
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
                           random = list(~ 1 | id, ~1 | effect),
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
                       random = list(~ 1 | id, ~1 | effect),
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
                       random = list(~ 1 | id, ~1 | effect),
                       mods = ~ ns(protein_intake, df=3) +
                         duration_centre + resistance_exercise_code + age_centre,
                       # struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_mods_ncs <- robust(meta_arm_mods_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=3), data=data_joined, absorb.cons=TRUE)[[1]]

meta_arm_mods_tps <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | effect),
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
                       random = list(~ ns(protein_intake, df=3) | id, ~1 | effect),
                       mods = ~ ns(protein_intake, df=3),
                       struct = "GEN",
                       data = data_joined,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_slopes_ncs <- robust(meta_arm_slopes_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=3), data=data_joined, absorb.cons=TRUE)[[1]]

data_joined <- bind_cols(data_joined, sm$X)

meta_arm_slopes_tps <- rma.mv(yi, vi,
                       random = list(~ `...33` + `...34` | id, ~1 | effect),
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
                            random = list(~ ns(protein_intake, df=3) | id, ~1 | effect),
                            mods = ~ ns(protein_intake, df=3) +
                              duration_centre + resistance_exercise_code + age_centre,
                            struct = "GEN",
                            data = data_joined,
                            method="REML", test="t",
                            # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_slopes_mods_ncs <- robust(meta_arm_slopes_mods_ncs, cluster = data_joined$id, clubSandwich = TRUE)

##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=3), data=data_joined, absorb.cons=TRUE)[[1]]

meta_arm_slopes_mods_tps <- rma.mv(yi, vi,
                            random = list(~ `...33` + `...34`  | id, ~1 | effect),
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
preds <- predict(rob_meta_arm_tps, newmods=PredictMat(sm, data.frame(protein_intake=seq(-0.7,2.9,0.1))), addx = TRUE)

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




pred_tps_plot <- preds |>
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


# slope_tps_plot <- preds |>
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

ggsave(plot = pred_tps_plot, filename = "plots/tps_plot.tiff", device = "tiff", dpi = 300, w = 7.5, h = 3.75)
