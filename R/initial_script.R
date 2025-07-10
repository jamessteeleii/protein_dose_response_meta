library(tidyverse)
library(metafor)
library(emmeans)
library(mgcv)
library(patchwork)
library(splines)
library(lspline)

study_characteristics <- read_csv("nunes_study_characteristics.csv") |>
  janitor::clean_names() |>
  select(id, ref, sex, age, duration, resistance_exercise, int_protein_intake, con_protein_intake) |>
  pivot_longer(
    7:8,
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
  group_by(id, condition) |>
  slice_head()

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

data <- escalc(measure = "ZCOR", ri = int_ri, ni = int_n, data = data)

meta_int_ri <- rma.mv(yi, V=vi, data=data,
                      random = ~ 1 | id / intervention / effect, method="REML", test="t",
                      control=list(optimizer="optim", optmethod="Nelder-Mead"))

rob_meta_int_ri <- robust(meta_int_ri, data$id)

z2r_int <- psych::fisherz2r(rob_meta_int_ri$b[1])

data$int_ri <- ifelse(is.na(data$int_ri), z2r_int, data$int_ri)

data_con <- data |>
  group_by(id) |>
  slice_head()

data_con <- escalc(measure = "ZCOR", ri = con_ri, ni = con_n, data = data_con)

### Note, data is coded with study and arm as having explicit nesting so all random effects are (~ 1 | study, ~ 1 | arm)
meta_con_ri <- rma.mv(yi, V=vi, data=data_con,
                      random = ~ 1 | id / effect, method="REML", test="t",
                      control=list(optimizer="optim", optmethod="Nelder-Mead"))

rob_meta_con_ri <- robust(meta_con_ri, data_con$id)

z2r_con <- psych::fisherz2r(rob_meta_con_ri$b[1])

data$con_ri <- ifelse(is.na(data$con_ri), z2r_con, data$con_ri)

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



data <- left_join(data, study_characteristics, by = c("id", "condition")) |>
  filter(!is.na(protein_intake))

# data <- escalc(
#   measure = "MC",
#   m1i = post_m,
#   m2i = pre_m,
#   sd1i = post_sd,
#   sd2i = pre_sd,
#   ni = n,
#   ri = ri,
#   data = data
# )

data <- escalc(
  measure = "SMCR",
  m1i = post_m,
  m2i = pre_m,
  sd1i = pre_sd_pool,
  # sd2i = pre_sd,
  ni = n,
  ri = ri,
  data = data
)

data <- data |>
  select(-effect) |>
  rowid_to_column("effect") |>
  # add study weights/sizes
  mutate(
    wi = 1/sqrt(vi),
    size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE))) |>
  
  # centre  predictors
  mutate(
    protein_intake = protein_intake - 0.8,
    duration_centre = duration - 12,
    age_centre = age - 40,
    resistance_exercise_code = case_when(
      resistance_exercise == "YES" ~ 0,
      resistance_exercise == "NO" ~ 1,
    )
  )

# data|>
#   ggplot(aes(x = protein_intake)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
#   geom_smooth(aes(y = yi, group = id), se=FALSE) +
#   scale_color_manual(values = c("darkorange", "darkgreen")) +
#   labs(
#     x = "Protein Intake (g/kg/day)",
#     y = "Lean Mass Change (SMD)",
#     # y = "Lean Mass Change (SMD)",
#     color = "Condition"
#   ) +
#   guides(
#     size = "none"
#   ) +
#   facet_grid(. ~ resistance_exercise) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(size = 12),
#     axis.title = element_text(size = 10)
#   )

##### Intercept only model model
meta_arm_intercept <- rma.mv(yi, vi,
                          random = list(~ 1 | id, ~1 | effect),
                          struct = "GEN",
                          data = data,
                          method="REML", test="t",
                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_intercept <- robust(meta_arm_intercept, cluster = data$id, clubSandwich = TRUE)

##### Linear model
meta_arm_linear <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | effect),
                       mods = ~ protein_intake,
                       struct = "GEN",
                       data = data,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_linear <- robust(meta_arm_linear, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_linear, newmods=cbind(seq(0.8,4.4,0.1)), addx = TRUE)


#### g marginal slope
rob_meta_arm_linear_slope_prep <- emmprep(rob_meta_arm_linear,
                                       at=list(protein_intake=seq(0.8,4.4,0.1)))

rob_meta_arm_linear_slope_prep <- confint(emmeans(rob_meta_arm_linear_slope_prep,consec~protein_intake,
                                               weights = "prop")$contrasts)

rob_meta_arm_linear_slope_prep2<-qt(.975,rob_meta_arm_linear_slope_prep$df)
rob_meta_arm_linear_slope_prep3<-sqrt((rob_meta_arm_linear_slope_prep$SE^2)+sum(rob_meta_arm_linear$sigma2))

rob_meta_arm_linear_slope_prep$lower.PI<-rob_meta_arm_linear_slope_prep$estimate-(rob_meta_arm_linear_slope_prep2*rob_meta_arm_linear_slope_prep3)
rob_meta_arm_linear_slope_prep$upper.PI<-rob_meta_arm_linear_slope_prep$estimate+(rob_meta_arm_linear_slope_prep2*rob_meta_arm_linear_slope_prep3)




pred_linear_plot <- preds |>
  as.data.frame() |>
  rename(protein_intake = "X.protein_intake") |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_linear_plot <- rob_meta_arm_linear_slope_prep |>
  mutate(protein_intake = seq(0.8,4.3,0.1)) |>
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


linear_plots <- (pred_linear_plot / slope_linear_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based linear meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = pred_linear_plot, filename = "plots/linear_plot.tiff", device = "tiff", dpi = 300, w = 7.5, h = 3.75)



stdres_plot_linear <- data.frame(fit = fitted(meta_arm_linear), rstand = rstandard(meta_arm_linear)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_linear <- data.frame(fit = fitted(meta_arm_linear), rstand = rstandard(meta_arm_linear)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_linear / res_plot_linear) 


##### Linear model
meta_arm_lspline <- rma.mv(yi, vi,
                          random = list(~ 1 | id, ~1 | effect),
                          mods = ~ lspline(protein_intake, 0.8),
                          struct = "GEN",
                          data = data,
                          method="REML", test="t",
                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_lspline <- robust(meta_arm_lspline, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_lspline, newmods=cbind(seq(0,3.6,0.1),seq(0,3.6,0.1), 0, 0, 0), addx = TRUE)

preds <- predict(rob_meta_arm_lspline, addx = TRUE)
# 
# 
# #### g marginal slope
# rob_meta_arm_lspline_slope_prep <- emmprep(rob_meta_arm_lspline,
#                                           at=list(protein_intake=seq(0.8,4.4,0.1)))
# 
# rob_meta_arm_lspline_slope_prep <- confint(emmeans(rob_meta_arm_lspline_slope_prep,consec~protein_intake,
#                                                   weights = "prop")$contrasts)
# 
# rob_meta_arm_lspline_slope_prep2<-qt(.975,rob_meta_arm_lspline_slope_prep$df)
# rob_meta_arm_lspline_slope_prep3<-sqrt((rob_meta_arm_lspline_slope_prep$SE^2)+sum(rob_meta_arm_lspline$sigma2))
# 
# rob_meta_arm_lspline_slope_prep$lower.PI<-rob_meta_arm_lspline_slope_prep$estimate-(rob_meta_arm_lspline_slope_prep2*rob_meta_arm_lspline_slope_prep3)
# rob_meta_arm_lspline_slope_prep$upper.PI<-rob_meta_arm_lspline_slope_prep$estimate+(rob_meta_arm_lspline_slope_prep2*rob_meta_arm_lspline_slope_prep3)
# 
# 
# 
# 
pred_lspline_plot <- preds |>
  as.data.frame() |>
  rename(protein_intake = "X.lspline.protein_intake..0.8.1") |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )
# 
# 
# slope_lspline_plot <- rob_meta_arm_lspline_slope_prep |>
#   mutate(protein_intake = seq(0.8,4.3,0.1)) |>
#   ggplot(aes(x = protein_intake)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
#   geom_line(aes(y=estimate)) +
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
# linear_plots <- (pred_lspline_plot / slope_lspline_plot) +
#   plot_annotation(
#     tag_levels = "i",
#     title = "Arm based linear meta-regression of protein intake on lean mass change",
#     subtitle = "Model includes random intercepts for study and effect"
#   ) &
#   theme(
#     legend.position = "bottom"
#   )
# 
# ggsave(plot = pred_lspline_plot, filename = "plots/linear_plot.tiff", device = "tiff", dpi = 300, w = 7.5, h = 3.75)
# 
# 
# 
# stdres_plot_lspline <- data.frame(fit = fitted(meta_arm_lspline), rstand = rstandard(meta_arm_lspline)$z) %>%
#   ggplot(aes(x=fit, y=rstand)) +
#   geom_hline(yintercept=0, linetype = "dashed") +
#   geom_point(alpha=0.5) +
#   geom_smooth() +
#   labs(x = "Fitted values",
#        y = "Standardised Residuals",
#        title = "Linearity",
#        subtitle = "Reference line should be flat and horizontal") +
#   theme_classic()
# 
# res_plot_lspline <- data.frame(fit = fitted(meta_arm_lspline), rstand = rstandard(meta_arm_lspline)$resid) %>%
#   ggplot(aes(x=fit, y=rstand)) +
#   geom_hline(yintercept=0, linetype = "dashed") +
#   geom_point(alpha=0.5) +
#   geom_smooth() +
#   labs(x = "Fitted values",
#        y = "Observed Residuals",
#        title = "Linearity",
#        subtitle = "Reference line should be flat and horizontal") +
#   theme_classic()
# 
# (stdres_plot_lspline / res_plot_lspline) 

##### Linear and log terms model
meta_arm_log <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | effect),
                       mods = ~ protein_intake + log(protein_intake),
                       struct = "GEN",
                       data = data,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_log <- robust(meta_arm_log, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_log, newmods=cbind(seq(0.8,4.4,0.1), log(seq(0.8,4.4,0.1))), addx = TRUE)


#### g marginal slope
rob_meta_arm_log_slope_prep <- emmprep(rob_meta_arm_log,
                                       at=list(protein_intake=seq(0.8,4.4,0.1)))

rob_meta_arm_log_slope_prep <- confint(emmeans(rob_meta_arm_log_slope_prep,consec~protein_intake,
                                               weights = "prop")$contrasts)

rob_meta_arm_log_slope_prep2<-qt(.975,rob_meta_arm_log_slope_prep$df)
rob_meta_arm_log_slope_prep3<-sqrt((rob_meta_arm_log_slope_prep$SE^2)+sum(rob_meta_arm_log$sigma2))

rob_meta_arm_log_slope_prep$lower.PI<-rob_meta_arm_log_slope_prep$estimate-(rob_meta_arm_log_slope_prep2*rob_meta_arm_log_slope_prep3)
rob_meta_arm_log_slope_prep$upper.PI<-rob_meta_arm_log_slope_prep$estimate+(rob_meta_arm_log_slope_prep2*rob_meta_arm_log_slope_prep3)




pred_log_plot <- preds |>
  as.data.frame() |>
  rename(protein_intake = "X.protein_intake") |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_log_plot <- rob_meta_arm_log_slope_prep |>
  mutate(protein_intake = seq(0.8,4.3,0.1)) |>
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


log_plots <- (pred_log_plot / slope_log_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based linear and log term meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = log_plots, filename = "plots/log_plots.tiff", device = "tiff", dpi = 300, w = 7.5, h = 7.5)



stdres_plot_log <- data.frame(fit = fitted(meta_arm_log), rstand = rstandard(meta_arm_log)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_log <- data.frame(fit = fitted(meta_arm_log), rstand = rstandard(meta_arm_log)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_log / res_plot_log) 


##### NCS model
meta_arm_ncs <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | effect),
                       mods = ~ ns(protein_intake, df=3),
                       struct = "GEN",
                       data = data,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


rob_meta_arm_ncs <- robust(meta_arm_ncs, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_ncs, newmods=unname(ns(seq(0.8,4.4,0.1), knots=c(1.23,1.56), Boundary.knots=c(0.8,4.4))), addx = TRUE)


#### g marginal slope
rob_meta_arm_ncs_slope_prep <- emmprep(rob_meta_arm_ncs,
                                       at=list(protein_intake=seq(0.8,4.4,0.1)))

rob_meta_arm_ncs_slope_prep <- confint(emmeans(rob_meta_arm_ncs_slope_prep,consec~protein_intake,
                                               weights = "prop")$contrasts)

rob_meta_arm_ncs_slope_prep2<-qt(.975,rob_meta_arm_ncs_slope_prep$df)
rob_meta_arm_ncs_slope_prep3<-sqrt((rob_meta_arm_ncs_slope_prep$SE^2)+sum(rob_meta_arm_ncs$sigma2))

rob_meta_arm_ncs_slope_prep$lower.PI<-rob_meta_arm_ncs_slope_prep$estimate-(rob_meta_arm_ncs_slope_prep2*rob_meta_arm_ncs_slope_prep3)
rob_meta_arm_ncs_slope_prep$upper.PI<-rob_meta_arm_ncs_slope_prep$estimate+(rob_meta_arm_ncs_slope_prep2*rob_meta_arm_ncs_slope_prep3)




pred_ncs_plot <- preds |>
  as.data.frame() |>
  mutate(protein_intake = seq(0.8,4.4,0.1)) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_ncs_plot <- rob_meta_arm_ncs_slope_prep |>
  mutate(protein_intake = seq(0.8,4.3,0.1)) |>
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


ncs_plots <- (pred_ncs_plot / slope_ncs_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based natural cubic spline meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = ncs_plots, filename = "plots/ncs_plots.tiff", device = "tiff", dpi = 300, w = 7.5, h = 7.5)


stdres_plot_ncs <- data.frame(fit = fitted(meta_arm_ncs), rstand = rstandard(meta_arm_ncs)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_ncs <- data.frame(fit = fitted(meta_arm_ncs), rstand = rstandard(meta_arm_ncs)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_ncs / res_plot_ncs) 


##### TPS model
sm <- smoothCon(s(protein_intake, bs="tp", k=3), data=data, absorb.cons=TRUE)[[1]]

meta_arm_tps <- rma.mv(yi, vi,
                       random = list(~ 1 | id, ~1 | effect),
                       mods = ~ sm$X,
                       struct = "GEN",
                       data = data,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_tps <- robust(meta_arm_tps, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_tps, newmods=PredictMat(sm, data.frame(protein_intake=seq(0.8,4.4,0.1))), addx = TRUE)


# #### g marginal slope
# rob_meta_arm_tps_slope_prep <- emmprep(rob_meta_arm_tps,
#                                        at=list(protein_intake=seq(0.8,4.4,0.1)))
# 
# rob_meta_arm_tps_slope_prep <- confint(emmeans(rob_meta_arm_tps_slope_prep,consec~protein_intake,
#                                                weights = "prop")$contrasts)
# 
# rob_meta_arm_tps_slope_prep2<-qt(.975,rob_meta_arm_tps_slope_prep$df)
# rob_meta_arm_tps_slope_prep3<-sqrt((rob_meta_arm_tps_slope_prep$SE^2)+sum(rob_meta_arm_tps$sigma2))
# 
# rob_meta_arm_tps_slope_prep$lower.PI<-rob_meta_arm_tps_slope_prep$estimate-(rob_meta_arm_tps_slope_prep2*rob_meta_arm_tps_slope_prep3)
# rob_meta_arm_tps_slope_prep$upper.PI<-rob_meta_arm_tps_slope_prep$estimate+(rob_meta_arm_tps_slope_prep2*rob_meta_arm_tps_slope_prep3)




pred_tps_plot <- preds |>
  as.data.frame() |>
  mutate(protein_intake = seq(0.8,4.4,0.1)) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_tps_plot <- preds |>
  as.data.frame()  |>
  mutate(protein_intake = seq(0.8,4.4,0.1),
         slope = c(NA, diff(pred) / diff(protein_intake)),
         lower.CL = c(NA, diff(ci.lb) / diff(protein_intake)),
         upper.CL = c(NA, diff(ci.ub) / diff(protein_intake))) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
  geom_line(aes(y=slope)) +
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


tps_plots <- (pred_tps_plot / slope_tps_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based thin plate spline meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = tps_plots, filename = "plots/tps_plots.tiff", device = "tiff", dpi = 300, w = 7.5, h = 7.5)


stdres_plot_tps <- data.frame(fit = fitted(meta_arm_tps), rstand = rstandard(meta_arm_tps)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_tps <- data.frame(fit = fitted(meta_arm_tps), rstand = rstandard(meta_arm_tps)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_tps / res_plot_tps) 



##### Linear model
meta_arm_linear_slopes <- rma.mv(yi, vi,
                          random = list(~ protein_intake | id, ~1 | effect),
                          mods = ~ protein_intake,
                          struct = "GEN",
                          data = data,
                          method="REML", test="t",
                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_linear_slopes <- robust(meta_arm_linear_slopes, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_linear_slopes, newmods=cbind(seq(0.8,4.4,0.1)), addx = TRUE)


#### g marginal slope
rob_meta_arm_linear_slopes_slope_prep <- emmprep(rob_meta_arm_linear_slopes,
                                          at=list(protein_intake=seq(0.8,4.4,0.1)))

rob_meta_arm_linear_slopes_slope_prep <- confint(emmeans(rob_meta_arm_linear_slopes_slope_prep,consec~protein_intake,
                                                  weights = "prop")$contrasts)

rob_meta_arm_linear_slopes_slope_prep2<-qt(.975,rob_meta_arm_linear_slopes_slope_prep$df)
rob_meta_arm_linear_slopes_slope_prep3<-sqrt((rob_meta_arm_linear_slopes_slope_prep$SE^2)+sum(rob_meta_arm_linear_slopes$sigma2))

rob_meta_arm_linear_slopes_slope_prep$lower.PI<-rob_meta_arm_linear_slopes_slope_prep$estimate-(rob_meta_arm_linear_slopes_slope_prep2*rob_meta_arm_linear_slopes_slope_prep3)
rob_meta_arm_linear_slopes_slope_prep$upper.PI<-rob_meta_arm_linear_slopes_slope_prep$estimate+(rob_meta_arm_linear_slopes_slope_prep2*rob_meta_arm_linear_slopes_slope_prep3)




pred_linear_slopes_plot <- preds |>
  as.data.frame() |>
  rename(protein_intake = "X.protein_intake") |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_linear_slopes_plot <- rob_meta_arm_linear_slopes_slope_prep |>
  mutate(protein_intake = seq(0.8,4.3,0.1)) |>
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


linear_plots <- (pred_linear_slopes_plot / slope_linear_slopes_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based linear meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = linear_plots, filename = "plots/linear_plots.tiff", device = "tiff", dpi = 300, w = 7.5, h = 7.5)



stdres_plot_linear_slopes <- data.frame(fit = fitted(meta_arm_linear_slopes), rstand = rstandard(meta_arm_linear_slopes)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_linear_slopes <- data.frame(fit = fitted(meta_arm_linear_slopes), rstand = rstandard(meta_arm_linear_slopes)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_linear_slopes / res_plot_linear_slopes) 

##### Linear model
meta_arm_lspline_slopes <- rma.mv(yi, vi,
                           random = list(~ lspline(protein_intake, 0.8) | id, ~1 | effect),
                           mods = ~ lspline(protein_intake, 0.8),
                           struct = "GEN",
                           data = data,
                           method="REML", test="t",
                           # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_lspline_slopes <- robust(meta_arm_lspline_slopes, cluster = data$id, clubSandwich = TRUE)

# preds <- predict(rob_meta_arm_lspline_slopes, newmods=cbind(seq(0.8,4.4,0.1)), addx = TRUE)
# 
# 
# #### g marginal slope
# rob_meta_arm_lspline_slopes_slope_prep <- emmprep(rob_meta_arm_lspline_slopes,
#                                           at=list(protein_intake=seq(0.8,4.4,0.1)))
# 
# rob_meta_arm_lspline_slopes_slope_prep <- confint(emmeans(rob_meta_arm_lspline_slopes_slope_prep,consec~protein_intake,
#                                                   weights = "prop")$contrasts)
# 
# rob_meta_arm_lspline_slopes_slope_prep2<-qt(.975,rob_meta_arm_lspline_slopes_slope_prep$df)
# rob_meta_arm_lspline_slopes_slope_prep3<-sqrt((rob_meta_arm_lspline_slopes_slope_prep$SE^2)+sum(rob_meta_arm_lspline_slopes$sigma2))
# 
# rob_meta_arm_lspline_slopes_slope_prep$lower.PI<-rob_meta_arm_lspline_slopes_slope_prep$estimate-(rob_meta_arm_lspline_slopes_slope_prep2*rob_meta_arm_lspline_slopes_slope_prep3)
# rob_meta_arm_lspline_slopes_slope_prep$upper.PI<-rob_meta_arm_lspline_slopes_slope_prep$estimate+(rob_meta_arm_lspline_slopes_slope_prep2*rob_meta_arm_lspline_slopes_slope_prep3)
# 
# 
# 
# 
# pred_lspline_slopes_plot <- preds |>
#   as.data.frame() |>
#   rename(protein_intake = "X.protein_intake") |>
#   ggplot(aes(x = protein_intake)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
#   geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
#   geom_line(aes(y=pred)) +
#   scale_color_manual(values = c("darkorange", "darkgreen")) +
#   labs(
#     x = "Protein Intake (g/kg/day)",
#     y = "Lean Mass Change (SMD)",
#     # y = "Lean Mass Change (SMD)",
#     color = "Condition",
#     title = "Predicted lean mass change"
#   ) +
#   guides(
#     size = "none"
#   ) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(size = 12),
#     axis.title = element_text(size = 10)
#   )
# 
# 
# slope_lspline_slopes_plot <- rob_meta_arm_lspline_slopes_slope_prep |>
#   mutate(protein_intake = seq(0.8,4.3,0.1)) |>
#   ggplot(aes(x = protein_intake)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
#   geom_line(aes(y=estimate)) +
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
# linear_plots <- (pred_lspline_slopes_plot / slope_lspline_slopes_plot) +
#   plot_annotation(
#     tag_levels = "i",
#     title = "Arm based linear meta-regression of protein intake on lean mass change",
#     subtitle = "Model includes random intercepts for study and effect"
#   ) &
#   theme(
#     legend.position = "bottom"
#   )
# 
# ggsave(plot = pred_lspline_slopes_plot, filename = "plots/linear_plot.tiff", device = "tiff", dpi = 300, w = 7.5, h = 3.75)
# 
# 
# 
# stdres_plot_lspline_slopes <- data.frame(fit = fitted(meta_arm_lspline_slopes), rstand = rstandard(meta_arm_lspline_slopes)$z) %>%
#   ggplot(aes(x=fit, y=rstand)) +
#   geom_hline(yintercept=0, linetype = "dashed") +
#   geom_point(alpha=0.5) +
#   geom_smooth() +
#   labs(x = "Fitted values",
#        y = "Standardised Residuals",
#        title = "Linearity",
#        subtitle = "Reference line should be flat and horizontal") +
#   theme_classic()
# 
# res_plot_lspline_slopes <- data.frame(fit = fitted(meta_arm_lspline_slopes), rstand = rstandard(meta_arm_lspline_slopes)$resid) %>%
#   ggplot(aes(x=fit, y=rstand)) +
#   geom_hline(yintercept=0, linetype = "dashed") +
#   geom_point(alpha=0.5) +
#   geom_smooth() +
#   labs(x = "Fitted values",
#        y = "Observed Residuals",
#        title = "Linearity",
#        subtitle = "Reference line should be flat and horizontal") +
#   theme_classic()
# 
# (stdres_plot_lspline_slopes / res_plot_lspline_slopes) 

##### Linear and log terms model
meta_arm_log_slopes <- rma.mv(yi, vi,
                       random = list(~ protein_intake + log(protein_intake) | id, ~1 | effect),
                       mods = ~ protein_intake + log(protein_intake),
                       struct = "GEN",
                       data = data,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_log_slopes <- robust(meta_arm_log_slopes, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_log_slopes, newmods=cbind(seq(0.8,4.4,0.1), log(seq(0.8,4.4,0.1))), addx = TRUE)


#### g marginal slope
rob_meta_arm_log_slopes_slope_prep <- emmprep(rob_meta_arm_log_slopes,
                                       at=list(protein_intake=seq(0.8,4.4,0.1)))

rob_meta_arm_log_slopes_slope_prep <- confint(emmeans(rob_meta_arm_log_slopes_slope_prep,consec~protein_intake,
                                               weights = "prop")$contrasts)

rob_meta_arm_log_slopes_slope_prep2<-qt(.975,rob_meta_arm_log_slopes_slope_prep$df)
rob_meta_arm_log_slopes_slope_prep3<-sqrt((rob_meta_arm_log_slopes_slope_prep$SE^2)+sum(rob_meta_arm_log_slopes$sigma2))

rob_meta_arm_log_slopes_slope_prep$lower.PI<-rob_meta_arm_log_slopes_slope_prep$estimate-(rob_meta_arm_log_slopes_slope_prep2*rob_meta_arm_log_slopes_slope_prep3)
rob_meta_arm_log_slopes_slope_prep$upper.PI<-rob_meta_arm_log_slopes_slope_prep$estimate+(rob_meta_arm_log_slopes_slope_prep2*rob_meta_arm_log_slopes_slope_prep3)




pred_log_slopes_plot <- preds |>
  as.data.frame() |>
  rename(protein_intake = "X.protein_intake") |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_log_slopes_plot <- rob_meta_arm_log_slopes_slope_prep |>
  mutate(protein_intake = seq(0.8,4.3,0.1)) |>
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


log_plots <- (pred_log_slopes_plot / slope_log_slopes_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based linear and log term meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = log_plots, filename = "plots/log_plots.tiff", device = "tiff", dpi = 300, w = 7.5, h = 7.5)



stdres_plot_log_slopes <- data.frame(fit = fitted(meta_arm_log_slopes), rstand = rstandard(meta_arm_log_slopes)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_log_slopes <- data.frame(fit = fitted(meta_arm_log_slopes), rstand = rstandard(meta_arm_log_slopes)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_log_slopes / res_plot_log_slopes) 





##### ncs_slopes model
meta_arm_ncs_slopes <- rma.mv(yi, vi,
                       random = list(~ ns(protein_intake, df=3) | id, ~1 | effect),
                       mods = ~ ns(protein_intake, df=3),
                       struct = "GEN",
                       data = data,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_ncs_slopes <- robust(meta_arm_ncs_slopes, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_ncs_slopes, newmods=unname(ns(seq(0.8,4.4,0.1), knots=c(1.23,1.56), Boundary.knots=c(0.8,4.4))), addx = TRUE)


#### g marginal slope
rob_meta_arm_ncs_slopes_slope_prep <- emmprep(rob_meta_arm_ncs_slopes,
                                       at=list(protein_intake=seq(0.8,4.4,0.1)))

rob_meta_arm_ncs_slopes_slope_prep <- confint(emmeans(rob_meta_arm_ncs_slopes_slope_prep,consec~protein_intake,
                                               weights = "prop")$contrasts)

rob_meta_arm_ncs_slopes_slope_prep2<-qt(.975,rob_meta_arm_ncs_slopes_slope_prep$df)
rob_meta_arm_ncs_slopes_slope_prep3<-sqrt((rob_meta_arm_ncs_slopes_slope_prep$SE^2)+sum(rob_meta_arm_ncs_slopes$sigma2))

rob_meta_arm_ncs_slopes_slope_prep$lower.PI<-rob_meta_arm_ncs_slopes_slope_prep$estimate-(rob_meta_arm_ncs_slopes_slope_prep2*rob_meta_arm_ncs_slopes_slope_prep3)
rob_meta_arm_ncs_slopes_slope_prep$upper.PI<-rob_meta_arm_ncs_slopes_slope_prep$estimate+(rob_meta_arm_ncs_slopes_slope_prep2*rob_meta_arm_ncs_slopes_slope_prep3)




pred_ncs_slopes_plot <- preds |>
  as.data.frame() |>
  mutate(protein_intake = seq(0.8,4.4,0.1)) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_ncs_slopes_plot <- rob_meta_arm_ncs_slopes_slope_prep |>
  mutate(protein_intake = seq(0.8,4.3,0.1)) |>
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


ncs_slopes_plots <- (pred_ncs_slopes_plot / slope_ncs_slopes_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based natural cubic spline meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = ncs_slopes_plots, filename = "plots/ncs_slopes_plots.tiff", device = "tiff", dpi = 300, w = 7.5, h = 7.5)


stdres_plot_ncs_slopes <- data.frame(fit = fitted(meta_arm_ncs_slopes), rstand = rstandard(meta_arm_ncs_slopes)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_ncs_slopes <- data.frame(fit = fitted(meta_arm_ncs_slopes), rstand = rstandard(meta_arm_ncs_slopes)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_ncs_slopes / res_plot_ncs_slopes) 


##### tps_slopes model
sm <- smoothCon(s(protein_intake, bs="tp", k=3), data=data, absorb.cons=TRUE)[[1]]

data <- bind_cols(data, sm$X)

meta_arm_tps_slopes <- rma.mv(yi, vi,
                       random = list(~ `...32` + `...33` | id, ~1 | effect),
                       mods = ~ sm$X,
                       struct = "GEN",
                       data = data,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

rob_meta_arm_tps_slopes <- robust(meta_arm_tps_slopes, cluster = data$id, clubSandwich = TRUE)

preds <- predict(rob_meta_arm_tps_slopes, newmods=PredictMat(sm, data.frame(protein_intake=seq(0.8,4.4,0.1))), addx = TRUE)


# #### g marginal slope
# rob_meta_arm_tps_slopes_slope_prep <- emmprep(rob_meta_arm_tps_slopes,
#                                        at=list(protein_intake=seq(0.8,4.4,0.1)))
# 
# rob_meta_arm_tps_slopes_slope_prep <- confint(emmeans(rob_meta_arm_tps_slopes_slope_prep,consec~protein_intake,
#                                                weights = "prop")$contrasts)
# 
# rob_meta_arm_tps_slopes_slope_prep2<-qt(.975,rob_meta_arm_tps_slopes_slope_prep$df)
# rob_meta_arm_tps_slopes_slope_prep3<-sqrt((rob_meta_arm_tps_slopes_slope_prep$SE^2)+sum(rob_meta_arm_tps_slopes$sigma2))
# 
# rob_meta_arm_tps_slopes_slope_prep$lower.PI<-rob_meta_arm_tps_slopes_slope_prep$estimate-(rob_meta_arm_tps_slopes_slope_prep2*rob_meta_arm_tps_slopes_slope_prep3)
# rob_meta_arm_tps_slopes_slope_prep$upper.PI<-rob_meta_arm_tps_slopes_slope_prep$estimate+(rob_meta_arm_tps_slopes_slope_prep2*rob_meta_arm_tps_slopes_slope_prep3)




pred_tps_slopes_plot <- preds |>
  as.data.frame() |>
  mutate(protein_intake = seq(0.8,4.4,0.1)) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = data, aes(y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.5) +
  geom_line(aes(y=pred)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  labs(
    x = "Protein Intake (g/kg/day)",
    y = "Lean Mass Change (SMD)",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Predicted lean mass change"
  ) +
  guides(
    size = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10)
  )


slope_tps_slopes_plot <- preds |>
  as.data.frame()  |>
  mutate(protein_intake = seq(0.8,4.4,0.1),
         slope = c(NA, diff(pred) / diff(protein_intake)),
         lower.CL = c(NA, diff(ci.lb) / diff(protein_intake)),
         upper.CL = c(NA, diff(ci.ub) / diff(protein_intake))) |>
  ggplot(aes(x = protein_intake)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.5) +
  geom_line(aes(y=slope)) +
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


tps_slopes_plots <- (pred_tps_slopes_plot / slope_tps_slopes_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Arm based thin plate spline meta-regression of protein intake on lean mass change",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for protein intake at study level"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = tps_slopes_plots, filename = "plots/tps_slopes_plots.tiff", device = "tiff", dpi = 300, w = 7.5, h = 7.5)


stdres_plot_tps_slopes <- data.frame(fit = fitted(meta_arm_tps_slopes), rstand = rstandard(meta_arm_tps_slopes)$z) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Standardised Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

res_plot_tps_slopes <- data.frame(fit = fitted(meta_arm_tps_slopes), rstand = rstandard(meta_arm_tps_slopes)$resid) %>%
  ggplot(aes(x=fit, y=rstand)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(alpha=0.5) +
  geom_smooth() +
  labs(x = "Fitted values",
       y = "Observed Residuals",
       title = "Linearity",
       subtitle = "Reference line should be flat and horizontal") +
  theme_classic()

(stdres_plot_tps_slopes / res_plot_tps_slopes) 



# Compare models

library(bayestestR)

BF_mean_models <- bayesfactor_models(rob_meta_arm_intercept,
                                     rob_meta_arm_linear,
                                     rob_meta_arm_linear_slopes,
                                     rob_meta_arm_lspline,
                                     rob_meta_arm_lspline_slopes,
                                     rob_meta_arm_log,
                                     rob_meta_arm_log_slopes,
                                     rob_meta_arm_ncs,
                                     rob_meta_arm_ncs_slopes,
                                     rob_meta_arm_tps,
                                     rob_meta_arm_tps_slopes)
BF_2log <- function(x) (2*x)

BF_mean_models <- as_tibble(as.matrix(BF_mean_models))  |>
  mutate_at(1:11, BF_2log) |>
  rowid_to_column("Denominator") |>
  mutate(Denominator = case_when(
    Denominator == 1 ~ "Intercept only model",
    Denominator == 2 ~ "Linear model",
    Denominator == 3 ~ "Linear model (random slopes)",
    Denominator == 4 ~ "Linear spline (1.6g/kg/day) model",
    Denominator == 5 ~ "Linear spline (1.6g/kg/day) model (random slopes)",
    Denominator == 6 ~ "Linear and log term model",
    Denominator == 7 ~ "Linear and log term model (random slopes)",
    Denominator == 8 ~  "Natural cubic spline model",
    Denominator == 9 ~  "Natural cubic spline model (random slopes)",
    Denominator == 10 ~  "Thin plate spline model",
    Denominator == 11 ~  "Thin plate spline model (random slopes)"
    
  )) |>
  rename("Intercept only model" = 2,
         "Linear model" = 3,
         "Linear model (random slopes)" = 4,
         "Linear spline (1.6g/kg/day) model" = 5,
         "Linear spline (1.6g/kg/day) model (random slopes)" = 6,
         "Linear and log term model" = 7,
         "Linear and log term model (random slopes)" = 8,
         "Natural cubic spline model" = 9,
         "Natural cubic spline model (random slopes)" = 10,
         "Thin plate spline model" = 11,
         "Thin plate spline model (random slopes)" = 12) |>
  pivot_longer(2:12, names_to = "Numerator", values_to = "logBF")

model_comparisons <- BF_mean_models |> 
  mutate(Denominator = factor(Denominator, levels= c(
    "Intercept only model",
    "Linear model",
    "Linear model (random slopes)",
    "Linear spline (1.6g/kg/day) model",
    "Linear spline (1.6g/kg/day) model (random slopes)",
    "Linear and log term model",
    "Linear and log term model (random slopes)",
    "Natural cubic spline model",
    "Natural cubic spline model (random slopes)",
    "Thin plate spline model",
    "Thin plate spline model (random slopes)")),
    Numerator = factor(Numerator, levels= c( 
      "Intercept only model",
      "Linear model",
      "Linear model (random slopes)",
      "Linear spline (1.6g/kg/day) model",
      "Linear spline (1.6g/kg/day) model (random slopes)",
      "Linear and log term model",
      "Linear and log term model (random slopes)",
      "Natural cubic spline model",
      "Natural cubic spline model (random slopes)",
      "Thin plate spline model",
      "Thin plate spline model (random slopes)")),
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



ggsave(plot = model_comparisons, filename = "plots/model_comparisons.tiff", device = "tiff", dpi = 300, w = 12.5, h = 10)


exp(3.09/2)



##### Variance model
data <- escalc(
  measure = "SDLN",
  mi = delta_m,
  sdi = delta_sd,
  ni = n,
  data = data
)

data_var <- data |>
  filter(delta_m != 0) |>
  select(-effect) |>
  rowid_to_column("effect") |>
  # add study weights/sizes
  mutate(
    wi = 1/sqrt(vi),
    size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))

meta_arm_var <- rma.mv(yi, vi,
                       random =  list(~ log(abs(delta_m)) | id, ~1|effect),
                       mods = ~ log(abs(delta_m)) + condition,
                       # struct = "GEN",
                       data = data_var,
                       method="REML", test="t",
                       # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

meta_arm_var <- robust(meta_arm_var, cluster = data_var$id, clubSandwich = TRUE)

preds <- predict(meta_arm_var, addx = TRUE)

# 
# #### g marginal slope
# meta_arm_var_slope_prep <- emmprep(meta_arm_var,
#                                    at=list(protein_intake=seq(-4.7,2,0.1)))
# 
# meta_arm_var_slope_prep <- confint(emmeans(meta_arm_var_slope_prep,consec~protein_intake,
#                                            weights = "prop")$contrasts)
# 
# meta_arm_var_slope_prep2<-qt(.975,meta_arm_var_slope_prep$df)
# meta_arm_var_slope_prep3<-sqrt((meta_arm_var_slope_prep$SE^2)+sum(meta_arm_var$sigma2))
# 
# meta_arm_var_slope_prep$lower.PI<-meta_arm_var_slope_prep$estimate-(meta_arm_var_slope_prep2*meta_arm_var_slope_prep3)
# meta_arm_var_slope_prep$upper.PI<-meta_arm_var_slope_prep$estimate+(meta_arm_var_slope_prep2*meta_arm_var_slope_prep3)
# 
# 


pred_var_plot <- preds |>
  as.data.frame() |>
  rename(log_delta_m = "X.log.abs.delta_m..",
         condition = "X.conditionint") |>
  mutate(
    condition = case_when(
      condition == 1 ~ "int",
      condition == 0 ~ "con"
    )
  ) |>
  ggplot(aes(x = log_delta_m)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = data_var, aes(x = log(delta_m), y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub, fill = condition), alpha = 0.5) +
  geom_line(aes(y=pred, colour = condition)) +
  scale_color_manual(values = c("gold", "black")) +
  scale_fill_manual(values = c("gold", "black")) +
  labs(
    x = "Log Standard Deviation of Change",
    y = "Log Mean of Change",
    # y = "Lean Mass Change (SMD)",
    color = "Condition",
    title = "Arm based meta-regression of log mean change on log standard deviation change in lean mass",
    subtitle = "Model includes random intercepts for study and effect, and random slopes for log mean at study level"
  ) +
  guides(
    size = "none",
    color = "none"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    legend.position = "bottom"
  )


pred_var_plot

ggsave(plot = pred_var_plot, filename = "plots/var_plot.tiff", device = "tiff", dpi = 300, w = 7.5, h = 5)



data |>
  ggplot(aes(x = protein_intake, y=n)) +
  geom_point()