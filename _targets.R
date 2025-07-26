# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(crew)


# Set target options:
tar_option_set(packages = c("here",
                            "tidyverse",
                            "metafor",
                            "mgcv",
                            "patchwork",
                            "splines",
                            "lspline",
                            "rstan",
                            "brms",
                            "loo",
                            "bayestestR",
                            "tidybayes",
                            "marginaleffects",
                            "grateful",
                            "progressr"),
               seed = 1988,  # <-- GLOBAL reproducible seed
               memory = "transient",
               format = "qs",
               garbage_collection = TRUE,
               storage = "worker",
               retrieval = "worker",
               controller = crew_controller_local(workers = 10)) 

# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions.R")

# List of targets
list(
  tar_target(
    rstan_options,
    set_rstan_options()
  ),
  
  ##### Analysis pipeline for pre-registration to demonstrate models ----
  # i.e., Preliminary re-analysis with Nunes et al.
  
  # Load and prepare data
  tar_target(
    nunes_study_characteristics_file,
    here("data", "nunes_study_characteristics.csv"),
    format = "file"
  ),
  
  tar_target(
    nunes_lean_mass_data_file,
    here("data", "nunes_lean_mass_data.csv"),
    format = "file"
  ),
  
  tar_target(
    nunes_data_prepared,
    read_prepare_nunes_data(nunes_study_characteristics_file, 
                            nunes_lean_mass_data_file)
  ),
  
  tar_target(
    steele_data,
    read_csv(url("https://github.com/jamessteeleii/Meta-Analysis-of-Variation-in-Resistance-Training/raw/refs/heads/main/data/Polito%20et%20al.%20RT%20Extracted%20Data.csv"))
  ),
  
  # Get priors from arm based meta-analysis of Steele et al. data 
  tar_target(
    steele_priors,
    get_steele_priors(steele_data)
  ),
  
  # Fit all candidate models to Nunes et al. data
  # Fit all CGM models, generate predictions, slopes, and test slopes ----
  tar_target(
    model_names,
    c(
      "Intercept only model",
      "Linear model",
      "Linear spline (1.6g/kg/day) model",
      "Thin plate spline model",
      "Intercept only model (+ moderators)",
      "Linear model (+ moderators)",
      "Linear spline (1.6g/kg/day) model (+ moderators)" ,
      "Thin plate spline model (+ moderators)",
      "Linear model (+ random slopes)",
      "Linear spline (1.6g/kg/day) model (+ random slopes)",
      "Thin plate spline model (+ random smooths)",
      "Linear model (+ moderators & random slopes)",
      "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
      "Thin plate spline model (+ moderators & random smooths)",
      
      # all models with the inclusion of resistance training and protein interaction
      "Linear model (+ RT interaction)",
      "Linear spline (1.6g/kg/day) model (+ RT interaction)",
      "Thin plate spline model (+ RT interaction smooth)",
      "Linear model (+ moderators & RT interaction)",
      "Linear spline (1.6g/kg/day) model (+ moderators & RT interaction)" ,
      "Thin plate spline model (+ moderators & RT interaction)",
      "Linear model (+ random slopes & RT interaction)",
      "Linear spline (1.6g/kg/day) model (+ random slopes & RT interaction)",
      "Thin plate spline model (+ random smooths & RT interaction smooth)",
      "Linear model (+ moderators & random slopes & RT interaction)",
      "Linear spline (1.6g/kg/day) model (+ moderators & random slopes & RT interaction)",
      "Thin plate spline model (+ moderators & random smooths & RT interaction smooth)"
    )
  ),
  
  tar_target(
    model_formulas,
    c(
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
        (1 | id) + 
        (1 | arm) + 
        (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                      s(protein_intake, bs='tp', k=4) +
                      (1 | id) + 
                      (1 | arm) + 
                      (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                 duration_centre + age_centre +
                            (1 | id) + 
                            (1 | arm) + 
                            (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake +
                           duration_centre + age_centre +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) +
                            duration_centre + age_centre +
                            (1 | id) + 
                          (1 | arm) + 
                          (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                           s(protein_intake, bs='tp', k=4) +
                           duration_centre + age_centre +
                        (1 | id) + 
                      (1 | arm) + 
                      (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake +
                         (protein_intake | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) +
                          (lspline(protein_intake, 0.48) | id) + 
                          (1 | arm) + 
                          (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                             s(protein_intake, bs='tp', k=4) +
                             s(protein_intake, id, bs='fs') +
                      (1 | arm) + 
                      (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                protein_intake +
                                  duration_centre + age_centre +
                                  (protein_intake | id) + 
                                (1 | arm) + 
                                (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                 lspline(protein_intake, 0.48) +
                                   duration_centre + age_centre +
                                   (lspline(protein_intake, 0.48) | id) + 
                                 (1 | arm) + 
                                 (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                  s(protein_intake, bs='tp', k=4) +
                                  s(protein_intake, id, bs='fs') +
                               duration_centre + age_centre +
                             (1 | arm) + 
                             (1 | effect)",
      
      # all models with the inclusion of resistance training and protein interaction
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake + protein_intake:resistance_exercise_code +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) + lspline(protein_intake, 0.48):resistance_exercise_code +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                      s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                      (1 | id) + 
                      (1 | arm) + 
                      (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake + protein_intake:resistance_exercise_code +
                           duration_centre + age_centre +
                         (1 | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) + lspline(protein_intake, 0.48):resistance_exercise_code +
                            duration_centre + age_centre +
                            (1 | id) + 
                          (1 | arm) + 
                          (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                           s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                           duration_centre + age_centre +
                        (1 | id) + 
                      (1 | arm) + 
                      (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                         protein_intake + protein_intake:resistance_exercise_code +
                         (protein_intake | id) + 
                         (1 | arm) + 
                         (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                          lspline(protein_intake, 0.48) + lspline(protein_intake, 0.48):resistance_exercise_code +
                          (lspline(protein_intake, 0.48) | id) + 
                          (1 | arm) + 
                          (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                             s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                             s(protein_intake, id, bs='fs') +
                      (1 | arm) + 
                      (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                protein_intake + protein_intake:resistance_exercise_code +
                                  duration_centre + age_centre +
                                  (protein_intake | id) + 
                                (1 | arm) + 
                                (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                 lspline(protein_intake, 0.48) + lspline(protein_intake, 0.48):resistance_exercise_code +
                                   duration_centre + age_centre +
                                   (lspline(protein_intake, 0.48) | id) + 
                                 (1 | arm) + 
                                 (1 | effect)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code + 
                                  s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                                  s(protein_intake, id, bs='fs') +
                               duration_centre + age_centre +
                             (1 | arm) + 
                             (1 | effect)"
    )
  ),
  
  tar_target(
    models_nunes,
     purrr::map2(
       model_formulas,
       model_names,
       ~ fit_candidate_models(
         formula = .x,
         priors = steele_priors,
         data = nunes_data_prepared
       )
     ) |> set_names(model_names)
  ),

  tar_target(
    models_nunes_weights,
    get_model_weights(models_nunes)
  ),

  tar_target(
    preds_nunes,
    get_model_avg_preds(models_nunes, models_nunes_weights)
  ),
  
  tar_target(
    slopes_nunes,
    get_model_avg_slopes(models_nunes, models_nunes_weights)
  ),
  
  tar_target(
    ropes_pd_nunes,
    get_rope_pd(slopes_nunes)
  ),
  
  tar_target(
    plot_nunes,
    plot_preds_slopes(nunes_data_prepared, preds_nunes, slopes_nunes)
  ),
  
  tar_target(
    plot_nunes_tiff,
    ggsave(
      plot = plot_nunes,
      filename = "plots/plot_nunes.tiff",
      device = "tiff",
      dpi = 300,
      w = 8,
      h = 8
    )
  ),
  
  # Determining baseline standard deviations to provide reference values for smallest effect size of interest of 0.1
  # Lean soft tissue/fat free mass from Benito et al. (https://doi.org/10.3390/ijerph17041285)
  tar_target(
    benito_data,
    read_csv(url("https://github.com/jamessteeleii/glp1_rt_cits/raw/refs/heads/master/data/benito_study_data.csv"))
  ),
  
  tar_target(
    lean_mass_sd,
    estimate_lean_mass_sd(benito_data)
  ),
  
  tar_target(
    buckner_data_file,
    here("data", 
         "buckner_data.csv"),
    format = "file"
  ),
  
  tar_target(
    buckner_data,
    read_csv(buckner_data_file)
  ),
  
  tar_target(
    muscle_thickness_sd,
    estimate_muscle_thickness_sd(buckner_data)
  ),
  
  tar_target(
    open_powerlifting_data_file,
    here("data", 
         "openpowerlifting-latest",
         "openpowerlifting-2025-07-12",
         "openpowerlifting-2025-07-12-6c3c0797.csv"),
    format = "file"
  ),
  
  tar_target(
    open_powerlifting_data,
    read_csv(open_powerlifting_data_file)
  ),
  
  tar_target(
    strength_sd,
    estimate_strength_sd(open_powerlifting_data)
  ),
  
  # Grateful report ----
  tar_target(grateful_report, cite_packages(out.dir = "grateful/", cite.tidyverse = TRUE, out.format = "pdf"))
  
)
