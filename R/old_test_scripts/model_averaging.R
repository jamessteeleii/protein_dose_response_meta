targets::tar_load(c(nunes_data_prepared,models_nunes))

library(tidyverse)
library(bayestestR)
library(brms)
library(lspline)
library(marginaleffects)
library(patchwork)
library(tidybayes)


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

# Extract posterior draws

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

rope <- stacked_df |>
  group_by(protein_intake, resistance_exercise) |>
  nest() |>
  mutate(
    rope_percent = map(data, ~ rope(.x$estimate*0.42, range = c(-0.1, 0.1), ci = 1))
  ) |>
  unnest(rope_percent) |>
  select(protein_intake, resistance_exercise, ROPE_low, ROPE_high, ROPE_Percentage) |>
  filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
  mutate(type = "rope")
  

pd <- stacked_df |>
  group_by(protein_intake, resistance_exercise) |>
  nest() |>
  mutate(
    rope_percent = map(data, ~ rope(.x$estimate*0.42, range = c(0.1, Inf), ci = 1))
  ) |>
  unnest(rope_percent) |>
  select(protein_intake, resistance_exercise, ROPE_low, ROPE_high, ROPE_Percentage) |>
  filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
  mutate(type = "pd")
  



slopes_plot <- stacked_df |>
  filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
  ggplot(aes(x = protein_intake, y = estimate*0.42)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
  geom_hline(yintercept = c(-0.1,0.1), linetype = "dashed", alpha = 0.75) +
  stat_slabinterval() +
  ggh4x::facet_nested(.~ resistance_exercise_label + resistance_exercise,
                      scales = "free") +
  # add mean and qi
  annotate("text", (summary$site_centred*100)+45, summary$estimate*0.218, label = round(summary$estimate*0.218,2), size = 3) +
  annotate("text", (summary$site_centred*100)+45, summary$conf.low*0.218, label = round(summary$conf.low*0.218,2), size = 3) +
  annotate("text", (summary$site_centred*100)+45, summary$conf.high*0.218, label = round(summary$conf.high*0.218,2), size = 3) +
  # add percent in rope
  annotate("rect", xmin = 7.5, xmax = 92.5, ymin = -0.525, ymax = -0.3583332, color = "black", fill = "white") +
  annotate("text", (rope_percents$site_centred*100)+50, -0.4833332, label = glue::glue("{round(rope_percents$rope_percent*100,2)}%"), size = 3) +
  annotate("text", x=50, y=-0.4, label="Percentage of Posterior Distibution Within ROPE [-0.1,0.1]", size = 3) +
  # add probability of positive effect
  annotate("rect", xmin = 7.5, xmax = 92.5, ymin = 0.3583332, ymax = 0.525, color = "black", fill = "white") +
  annotate("text", (rope_percents$site_centred*100)+50, 0.4833332, label = glue::glue("{round(rope_percents$pd*100,2)}%"), size = 3) +
  annotate("text", x=50, y=0.4, label="Probability of Meaningful Positive Effect (i.e., >0.1)", size = 3)
  
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



# weighted_preds <- map2(preds, weights_stacking, function(pred_df, w) {
#   pred_df %>%
#     mutate(
#       .pred_weighted = estimate * w,
#       .lower_weighted = conf.low * w,
#       .upper_weighted = conf.high * w
#     )
# })
# 
# combined_preds <- reduce(weighted_preds, function(df1, df2) {
#   df1 %>%
#     mutate(
#       .pred_weighted = .pred_weighted + df2$.pred_weighted,
#       .lower_weighted = .lower_weighted + df2$.lower_weighted,
#       .upper_weighted = .upper_weighted + df2$.upper_weighted
#     )
# }) 
# 
# final_predictions <- combined_preds %>%
#   mutate(
#     estimate = .pred_weighted,
#     conf.low = .lower_weighted,
#     conf.high = .upper_weighted
#   ) |>
#   mutate(
#     resistance_exercise_label = "Resistance Exercise?",
#     resistance_exercise = case_when(
#       resistance_exercise_code == 0 ~ "YES",
#       resistance_exercise_code == 1 ~ "NO"
#     ),
#     protein_intake = protein_intake + 1.12
#   )




preds_plot <- final_predictions |> 
  filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
  ggplot(aes(x = protein_intake, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
  geom_point(data = nunes_data_prepared |>
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





safe_slopes <- possibly(slopes, otherwise = NULL) # because can't run on intercept model with no protein_intake

slopes <- map(models_nunes, ~ safe_slopes(.x, 
                                         re_formula= NA,
                                     variables = "protein_intake",
                                     newdata   = datagrid(
                                           resistance_exercise_code = c(0,1),
                                           protein_intake = seq(-0.32, 3.28, by = 0.01),
                                           age_centre = 0,
                                           duration_centre = 0,
                                           vi = 0
                                         )))


# Step 1: Identify a template (non-null slopes data frame) to use as a dummy
template <- compact(slopes)[[1]]  # first non-null slopes result

# Step 2: Replace NULLs with zero-weighted dummy versions
slopes_filled <- map2(slopes, weights_stacking, function(s, w) {
  if (is.null(s)) {
    template %>%
      mutate(
        estimate = 0,
        conf.low = 0,
        conf.high = 0
      )
  } else {
    s
  }
})

weighted_slopes <- map2(slopes_filled, weights_stacking, function(slope_df, w) {
  slope_df %>%
    mutate(
      .slope_weighted = estimate * w,
      .lower_weighted = conf.low * w,
      .upper_weighted = conf.high * w
    )
})

combined_slopes <- reduce(weighted_slopes, function(df1, df2) {
  df1 %>%
    mutate(
      .slope_weighted = .slope_weighted + df2$.slope_weighted,
      .lower_weighted = .lower_weighted + df2$.lower_weighted,
      .upper_weighted = .upper_weighted + df2$.upper_weighted
    )
}) 

final_slopes <- combined_slopes %>%
  mutate(
    estimate = .slope_weighted,
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


slopes_plot <- final_slopes |>
  filter(resistance_exercise == "YES" | protein_intake < 1.6) |>
  ggplot(aes(x = protein_intake, y = estimate*0.42)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
  geom_hline(yintercept = c(-0.1,0.1), linetype = "dashed", alpha = 0.75) +
  geom_ribbon(aes(ymin = conf.low*0.42, ymax = conf.high*0.42),
              alpha = 0.5, color = "black", size = 0.25) +
  geom_line(size = 1, color = "black") +
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
    axis.title = element_text(size = 10)
  )     


nunes_plots <- (preds_plot / slopes_plot) +
  plot_annotation(
    tag_levels = "i",
    title = "Dose-response effects of protein intake on lean mass change",
    subtitle = "Bayesian model averaging with LOO-CV stacking weights"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(plot = nunes_plots, filename = "plots/nunes_plots.tiff", device = "tiff", dpi = 300, w = 10, h = 10)











# example with just three of the models
m_intercept_preds <- predictions(m_intercept,
                                 # variables = "protein_intake",
                                 re_formula= NA,
                                 newdata   = datagrid(
                                   resistance_exercise_code = c(0,1),
                                   protein_intake = seq(-0.32, 3.28, by = 0.01),
                                   age_centre = 0,
                                   duration_centre = 0,
                                   vi = 0
                                 )) 
  


m_lspline_mods_preds <- predictions(m_lspline_mods,
                                 # variables = "protein_intake",
                                 re_formula= NA,
                                 newdata   = datagrid(
                                   resistance_exercise_code = c(0,1),
                                   protein_intake = seq(-0.32, 3.28, by = 0.01),
                                   age_centre = 0,
                                   duration_centre = 0,
                                   vi = 0
                                 ))

m_tps_ranef_preds <- predictions(m_tps_ranef,
                                    # variables = "protein_intake",
                                    re_formula= NA,
                                    newdata   = datagrid(
                                      resistance_exercise_code = c(0,1),
                                      protein_intake = seq(-0.32, 3.28, by = 0.01),
                                      age_centre = 0,
                                      duration_centre = 0,
                                      vi = 0
                                      # id = NA, arm = NA, effect = NA
                                    ))


weights_keep <- weights_stacking[c(1,7,11)] / sum(weights_stacking[c(1,7,11)])


avg_preds <- tibble(
  estimate = (weights_keep[1] * m_intercept_preds$estimate +
     weights_keep[2] * m_lspline_mods_preds$estimate +
     weights_keep[3] * m_tps_ranef_preds$estimate),
  conf.low = (weights_keep[1] * m_intercept_preds$conf.low +
                weights_keep[2] * m_lspline_mods_preds$conf.low +
                weights_keep[3] * m_tps_ranef_preds$conf.low),
  conf.high = (weights_keep[1] * m_intercept_preds$conf.high +
                weights_keep[2] * m_lspline_mods_preds$conf.high +
                weights_keep[3] * m_tps_ranef_preds$conf.high)
) |>
  bind_cols(crossing(
    resistance_exercise_code = c(0,1),
    protein_intake = seq(-0.32, 3.28, by = 0.01),
    age_centre = 0,
    duration_centre = 0,
    vi = 0,
    id = NA, arm = NA, effect = NA
  )) |>
  mutate(
    resistance_exercise_label = "Resistance Exercise?",
    resistance_exercise = case_when(
      resistance_exercise_code == 0 ~ "YES",
      resistance_exercise_code == 1 ~ "NO"
    ),
    protein_intake = protein_intake + 1.12
  )

preds_plot <- avg_preds |>
  ggplot(aes(x = protein_intake, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.75) +
  geom_point(data = nunes_data_prepared |>
               mutate(resistance_exercise_label = "Resistance Exercise?"), 
             aes(x = protein_intake + 1.12, y = yi, size = size, color = condition), alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.5, color = "black", size = 0.25) +
  geom_line(aes(y = estimate), size = 0.75, color = "black") +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  ggh4x::facet_nested(.~ resistance_exercise_label + resistance_exercise) +
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

m_intercept_slopes <- slopes(m_intercept,
                                 variables = "protein_intake",
                                 re_formula= NA,
                                 newdata   = datagrid(
                                   resistance_exercise_code = c(0,1),
                                   protein_intake = seq(-0.32, 3.28, by = 0.01),
                                   age_centre = 0,
                                   duration_centre = 0,
                                   vi = 0
                                 )) 



m_lspline_mods_slopes <- slopes(m_lspline_mods,
                                    variables = "protein_intake",
                                    re_formula= NA,
                                    newdata   = datagrid(
                                      resistance_exercise_code = c(0,1),
                                      protein_intake = seq(-0.32, 3.28, by = 0.01),
                                      age_centre = 0,
                                      duration_centre = 0,
                                      vi = 0
                                    ))

m_tps_ranef_slopes <- slopes(m_tps_ranef,
                                 variables = "protein_intake",
                                 re_formula= NA,
                                 newdata   = datagrid(
                                   resistance_exercise_code = c(0,1),
                                   protein_intake = seq(-0.32, 3.28, by = 0.01),
                                   age_centre = 0,
                                   duration_centre = 0,
                                   vi = 0
                                 ))


weights_keep <- weights_stacking[c(1,7,11)] / sum(weights_stacking[c(1,7,11)])


avg_slopes <- tibble(
  estimate = (weights_keep[1] * 0 +
                weights_keep[2] * m_lspline_mods_slopes$estimate +
                weights_keep[3] * m_tps_ranef_slopes$estimate),
  conf.low = (weights_keep[1] * 0 +
                weights_keep[2] * m_lspline_mods_slopes$conf.low +
                weights_keep[3] * m_tps_ranef_slopes$conf.low),
  conf.high = (weights_keep[1] * 0 +
                 weights_keep[2] * m_lspline_mods_slopes$conf.high +
                 weights_keep[3] * m_tps_ranef_slopes$conf.high)
) |>
  bind_cols(crossing(
    resistance_exercise_code = c(0,1),
    protein_intake = seq(-0.32, 3.28, by = 0.01),
    age_centre = 0,
    duration_centre = 0,
    vi = 0,
    id = NA, arm = NA, effect = NA
  )) |>
  mutate(
    resistance_exercise_label = "Resistance Exercise?",
    resistance_exercise = case_when(
      resistance_exercise_code == 0 ~ "YES",
      resistance_exercise_code == 1 ~ "NO"
    ),
    protein_intake = protein_intake + 1.12
  )


slopes_plot <- avg_slopes |>
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


