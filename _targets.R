# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(crew)


# Set target options:
tar_option_set(packages = c("here",
                            "qs",
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
               retrieval = "worker"
               # controller = crew_controller_local(workers = 4)
               ) 

# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions.R")

# List of targets
list(
  tar_target(
    rstan_options,
    set_rstan_options()
  ),
  
  # Analysis pipeline for pre-registration to demonstrate models ----
  # i.e., Preliminary re-analysis with Nunes et al.
  
  ## Load and prepare data ----
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
    steele_data_file,
    here(
      "data",
      "external",
      "steele_rt_extracted_data_289a48f.csv"
    ),
    format = "file"
  ),

  tar_target(
    steele_data,
    read_csv(steele_data_file)
  ),
  
  ## Get priors from arm based meta-analysis of Steele et al. data ----
  # Note these are for hypertrophy (lean mass/muscle mass/size) only as Nunes reanalysis was only with these outcomes
  tar_target(
    steele_priors,
    get_steele_priors(steele_data)
  ),
  
  ## Fit all candidate models to Nunes et al. data ----
  # Fit all CGM models, generate predictions, slopes, and test slopes 
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
    {
      output_file <- here("plots", "plot_nunes.tiff")
      ggsave(
        plot = plot_nunes,
        filename = output_file,
        device = "tiff",
        dpi = 300,
        w = 8,
        h = 8
      )
      output_file
    },
    format = "file"
  ),
  
  # Original brms pre-registered updated analysis with all re-extracted data ----
  
  ## Read and prepare all new data ----
  
  tar_target(
    protein_data_file,
    here("data", "protein_meta_data_all.xlsx"),
    format = "file"
  ),
  
  tar_target(
    protein_data,
    read_protein_workbook(protein_data_file)
  ),
  
  # Note this is the data set where centering of protein intake is performed
  tar_target(
    protein_data_prepared,
    prepare_protein_data(protein_data)
  ),
  
  # Note this is the complete row wise data set for analysis
  tar_target(
    protein_data_complete,
    protein_data_prepared |>
      dplyr::filter(
        dplyr::if_all(
          c(
            yi, vi, protein_intake,
            resistance_exercise_code, age_centre, duration_centre
          ),
          ~ !is.na(.x)
        )
      )
  ),
  
  ## Get additional priors for strength and function outcomes ----
  # Note these are for hypertrophy (lean mass/muscle mass/size) only as Nunes reanalysis was only with these outcomes
  tar_target(
    steele_strength_priors,
    get_steele_strength_priors(steele_data)
  ),
  
  tar_target(
    kirk_function_priors,
    get_kirk_function_priors()
  ),
  
  
  ## Fit all candidate models to new data ----
  # Fit all CGM models, generate predictions, slopes, and test slopes
  
  # Note, as control median intake slightly different (1.1 g/kg/day) in the new data
  # The knot centred protein intake for lspline has been shifted to 0.5 to reflect the 1.6 g/kg/day
  tar_target(
    model_formulas_new,
    c(
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
        (1 | study_no) +
        (1 | arm_no) +
        (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                         protein_intake +
                         (1 | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                          lspline(protein_intake, 0.5) +
                         (1 | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                      s(protein_intake, bs='tp', k=4) +
                      (1 | study_no) +
                      (1 | arm_no) +
                      (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                                 duration_centre + age_centre +
                            (1 | study_no) +
                            (1 | arm_no) +
                            (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                         protein_intake +
                           duration_centre + age_centre +
                         (1 | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                          lspline(protein_intake, 0.5) +
                            duration_centre + age_centre +
                            (1 | study_no) +
                          (1 | arm_no) +
                          (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                           s(protein_intake, bs='tp', k=4) +
                           duration_centre + age_centre +
                        (1 | study_no) +
                      (1 | arm_no) +
                      (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                         protein_intake +
                         (protein_intake | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                          lspline(protein_intake, 0.5) +
                          (lspline(protein_intake, 0.5) | study_no) +
                          (1 | arm_no) +
                          (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                             s(protein_intake, bs='tp', k=4) +
                             s(protein_intake, study_no, bs='fs') +
                      (1 | arm_no) +
                      (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                                protein_intake +
                                  duration_centre + age_centre +
                                  (protein_intake | study_no) +
                                (1 | arm_no) +
                                (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                                 lspline(protein_intake, 0.5) +
                                   duration_centre + age_centre +
                                   (lspline(protein_intake, 0.5) | study_no) +
                                 (1 | arm_no) +
                                 (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                                  s(protein_intake, bs='tp', k=4) +
                                  s(protein_intake, study_no, bs='fs') +
                               duration_centre + age_centre +
                             (1 | arm_no) +
                             (1 | outcome_no)",

      # all models with the inclusion of resistance training and protein interaction
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                         protein_intake + protein_intake:resistance_exercise_code +
                         (1 | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                          lspline(protein_intake, 0.5) + lspline(protein_intake, 0.5):resistance_exercise_code +
                         (1 | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                      s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                      (1 | study_no) +
                      (1 | arm_no) +
                      (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                         protein_intake + protein_intake:resistance_exercise_code +
                           duration_centre + age_centre +
                         (1 | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                          lspline(protein_intake, 0.5) + lspline(protein_intake, 0.5):resistance_exercise_code +
                            duration_centre + age_centre +
                            (1 | study_no) +
                          (1 | arm_no) +
                          (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                           s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                           duration_centre + age_centre +
                        (1 | study_no) +
                      (1 | arm_no) +
                      (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                         protein_intake + protein_intake:resistance_exercise_code +
                         (protein_intake | study_no) +
                         (1 | arm_no) +
                         (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                          lspline(protein_intake, 0.5) + lspline(protein_intake, 0.5):resistance_exercise_code +
                          (lspline(protein_intake, 0.5) | study_no) +
                          (1 | arm_no) +
                          (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                             s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                             s(protein_intake, study_no, bs='fs') +
                      (1 | arm_no) +
                      (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                                protein_intake + protein_intake:resistance_exercise_code +
                                  duration_centre + age_centre +
                                  (protein_intake | study_no) +
                                (1 | arm_no) +
                                (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                                 lspline(protein_intake, 0.5) + lspline(protein_intake, 0.5):resistance_exercise_code +
                                   duration_centre + age_centre +
                                   (lspline(protein_intake, 0.5) | study_no) +
                                 (1 | arm_no) +
                                 (1 | outcome_no)",
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +
                                  s(protein_intake, by=resistance_exercise_code, m=1, bs='tp', k=4) +
                                  s(protein_intake, study_no, bs='fs') +
                               duration_centre + age_centre +
                             (1 | arm_no) +
                             (1 | outcome_no)"
    )
  ),

  ## Whole body lean mass ----
  tar_target(
    protein_data_complete_whole_body,
    {
      protein_data_complete |> filter(outcome_group == "whole body lean mass") |>
        # add study weights/sizes
        mutate(
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),

  tar_target(
    models_new_whole_body,
    purrr::map2(
      model_formulas_new,
      model_names,
      ~ fit_candidate_models(
        formula = .x,
        priors = steele_priors,
        data = protein_data_complete_whole_body
      )
    ) |> set_names(model_names)
  ),

  tar_target(
    models_new_weights_whole_body,
    get_model_weights(models_new_whole_body)
  ),

  tar_target(
    preds_new_whole_body,
    get_model_avg_preds_new(models_new_whole_body, models_new_weights_whole_body, protein_data_complete_whole_body)
  ),

  tar_target(
    contrasts_new_whole_body_0.8,
    get_model_avg_contrasts_new(models_new_whole_body, models_new_weights_whole_body, protein_data_complete_whole_body)
  ),
  
  tar_target(
    contrasts_new_whole_body_1.1,
    get_model_avg_contrasts_new(models_new_whole_body, models_new_weights_whole_body, protein_data_complete_whole_body,
                                reference = 1.1)
  ),

  tar_target(
    ropes_pd_new_whole_body_0.8,
    get_rope_pd_new(contrasts_new_whole_body_0.8)
  ),
  
  tar_target(
    ropes_pd_new_whole_body_1.1,
    get_rope_pd_new(contrasts_new_whole_body_1.1)
  ),

  tar_target(
    plot_new_whole_body,
    plot_preds_contrasts_new(protein_data_complete_whole_body,
                      preds_new_whole_body, contrasts_new_whole_body_0.8, contrasts_new_whole_body_1.1,
                      "whole body lean soft tissue change")
  ),

  tar_target(
    plot_new_whole_body_tiff,
    {
      output_file <- here("plots", "plot_new_whole_body.tiff")
      ggsave(
        plot = plot_new_whole_body,
        filename = output_file,
        device = "tiff",
        dpi = 300,
        w = 8,
        h = 12
      )
      output_file
    },
    format = "file"
  ),

  ## Hypertrophy ----
  tar_target(
    protein_data_complete_hypertrophy,
    {
      protein_data_complete |> filter(outcome_group == "direct muscle size/mass") |>
        # add study weights/sizes
        mutate(
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),

  tar_target(
    models_new_hypertrophy,
    purrr::map2(
      model_formulas_new,
      model_names,
      ~ fit_candidate_models(
        formula = .x,
        priors = steele_priors,
        data = protein_data_complete_hypertrophy
      )
    ) |> set_names(model_names)
  ),

  tar_target(
    models_new_weights_hypertrophy,
    get_model_weights(models_new_hypertrophy)
  ),

  tar_target(
    preds_new_hypertrophy,
    get_model_avg_preds_new(models_new_hypertrophy, models_new_weights_hypertrophy, protein_data_complete_hypertrophy)
  ),

  tar_target(
    contrasts_new_hypertrophy_0.8,
    get_model_avg_contrasts_new(models_new_hypertrophy, models_new_weights_hypertrophy, protein_data_complete_hypertrophy)
  ),
  
  tar_target(
    contrasts_new_hypertrophy_1.1,
    get_model_avg_contrasts_new(models_new_hypertrophy, models_new_weights_hypertrophy, protein_data_complete_hypertrophy,
                                reference = 1.1)
  ),

  tar_target(
    ropes_pd_new_hypertrophy_0.8,
    get_rope_pd_new(contrasts_new_hypertrophy_0.8)
  ),
  
  tar_target(
    ropes_pd_new_hypertrophy_1.1,
    get_rope_pd_new(contrasts_new_hypertrophy_1.1)
  ),

  tar_target(
    plot_new_hypertrophy,
    plot_preds_contrasts_new(protein_data_complete_hypertrophy,
                      preds_new_hypertrophy, contrasts_new_hypertrophy_0.8, contrasts_new_hypertrophy_1.1,
                      "direct muscle mass/size change")
  ),

  tar_target(
    plot_new_hypertrophy_tiff,
    {
      output_file <- here("plots", "plot_new_hypertrophy.tiff")
      ggsave(
        plot = plot_new_hypertrophy,
        filename = output_file,
        device = "tiff",
        dpi = 300,
        w = 8,
        h = 12
      )
      output_file
    },
    format = "file"
  ),

  ## Strength ----
  tar_target(
    protein_data_complete_strength,
    {
      protein_data_complete |> filter(outcome_group == "strength") |>
        # add study weights/sizes
        mutate(
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),

  tar_target(
    models_new_strength,
    purrr::map2(
      model_formulas_new,
      model_names,
      ~ fit_candidate_models(
        formula = .x,
        priors = steele_strength_priors,
        data = protein_data_complete_strength
      )
    ) |> set_names(model_names)
  ),

  tar_target(
    models_new_weights_strength,
    get_model_weights(models_new_strength)
  ),

  tar_target(
    preds_new_strength,
    get_model_avg_preds_new(models_new_strength, models_new_weights_strength, protein_data_complete_strength)
  ),

  tar_target(
    contrasts_new_strength_0.8,
    get_model_avg_contrasts_new(models_new_strength, models_new_weights_strength, protein_data_complete_strength)
  ),
  
  tar_target(
    contrasts_new_strength_1.1,
    get_model_avg_contrasts_new(models_new_strength, models_new_weights_strength, protein_data_complete_strength,
                                reference = 1.1)
  ),

  tar_target(
    ropes_pd_new_strength_0.8,
    get_rope_pd_new(contrasts_new_strength_0.8)
  ),
  
  tar_target(
    ropes_pd_new_strength_1.1,
    get_rope_pd_new(contrasts_new_strength_1.1)
  ),

  tar_target(
    plot_new_strength,
    plot_preds_contrasts_new(protein_data_complete_strength,
                          preds_new_strength, contrasts_new_strength_0.8, contrasts_new_strength_1.1,
                          "strength change")
  ),

  tar_target(
    plot_new_strength_tiff,
    {
      output_file <- here("plots", "plot_new_strength.tiff")
      ggsave(
        plot = plot_new_strength,
        filename = output_file,
        device = "tiff",
        dpi = 300,
        w = 8,
        h = 12
      )
      output_file
    },
    format = "file"
  ),
  
  ## Function ----
  tar_target(
    protein_data_complete_function,
    {
      protein_data_complete |> filter(outcome_group == "function") |>
        
        # add sign for function outcomes and convert effect sizes to correct sign
        mutate(
          sign = case_when(
           unit_of_measurement == "score" ~ "increase",
           unit_of_measurement == "seconds" ~ "decrease",
           unit_of_measurement == "m/s" ~ "increase",
           unit_of_measurement == "watts" ~ "increase",
           unit_of_measurement == "number" ~ "increase",
           unit_of_measurement == "metres" ~ "increase",
           unit_of_measurement == "cm" ~ "increase",
           unit_of_measurement == "reps" ~ "increase",
           unit_of_measurement == "duration" ~ "increase",
           unit_of_measurement == "minutes" ~ "decrease"
          ),
          yi = yi * if_else(sign == "increase", 1, -1)
        ) |>
        
        # add study weights/sizes
        mutate(
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),
  
  tar_target(
    models_new_function,
    purrr::map2(
      model_formulas_new,
      model_names,
      ~ fit_candidate_models(
        formula = .x,
        priors = kirk_function_priors,
        data = protein_data_complete_function
      )
    ) |> set_names(model_names)
  ),
  
  tar_target(
    models_new_weights_function,
    get_model_weights(models_new_function)
  ),
  
  tar_target(
    preds_new_function,
    get_model_avg_preds_new(models_new_function, models_new_weights_function, protein_data_complete_function)
  ),
  
  tar_target(
    contrasts_new_function_0.8,
    get_model_avg_contrasts_new(models_new_function, models_new_weights_function, protein_data_complete_function)
  ),
  
  tar_target(
    contrasts_new_function_1.1,
    get_model_avg_contrasts_new(models_new_function, models_new_weights_function, protein_data_complete_function,
                                reference = 1.1)
  ),
  
  tar_target(
    ropes_pd_new_function_0.8,
    get_rope_pd_new(contrasts_new_function_0.8)
  ),
  
  tar_target(
    ropes_pd_new_function_1.1,
    get_rope_pd_new(contrasts_new_function_1.1)
  ),
  
  tar_target(
    plot_new_function,
    plot_preds_contrasts_new(protein_data_complete_function,
                             preds_new_function, contrasts_new_function_0.8, contrasts_new_function_1.1,
                             "function change")
  ),
  
  tar_target(
    plot_new_function_tiff,
    {
      output_file <- here("plots", "plot_new_function.tiff")
      ggsave(
        plot = plot_new_function,
        filename = output_file,
        device = "tiff",
        dpi = 300,
        w = 8,
        h = 12
      )
      output_file
    },
    format = "file"
  ),
  
  # Amended pre-registration updated reanalysis converting brms models to stan to enable error on predictors ----
  
  # Note this is the complete row wise data set for analysis including only arms with error on protein intake
  tar_target(
    protein_data_complete_me,
    protein_data_prepared |>
      dplyr::filter(
        dplyr::if_all(
          c(
            yi, vi, protein_intake, protein_intake_se,
            resistance_exercise_code, age_centre, duration_centre
          ),
          ~ !is.na(.x)
        )
      )
  ),
  
  tar_target(
    stan_data_whole_body,
    make_stan_data(protein_data_complete_me |> filter(outcome_group == "whole body lean mass"))
  ),
  
  tar_target(
    stan_data_hypertrophy,
    make_stan_data(protein_data_complete_me |> filter(outcome_group == "direct muscle size/mass"))
  ),
  
  tar_target(
    stan_data_strength,
    make_stan_data(protein_data_complete_me |> filter(outcome_group == "strength"))
  ),
  
  tar_target(
    stan_data_function,
    make_stan_data(protein_data_complete_me |> filter(outcome_group == "function"))
  ),
  
  ## Get priors from old brms prior target ----
  tar_target(
    stan_priors,
    make_stan_priors(steele_priors)
  ),

  
  
  # Determining baseline standard deviations to provide reference values for smallest effect size of interest of 0.1 ----
  # Lean soft tissue/fat free mass from Benito et al. (https://doi.org/10.3390/ijerph17041285)
  tar_target(
    benito_data_file,
    here(
      "data",
      "external",
      "benito_study_data_18f3fc7.csv"
    ),
    format = "file"
  ),

  tar_target(
    benito_data,
    read_csv(benito_data_file)
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
    open_powerlifting_url,
    paste0(
      "https://openpowerlifting.gitlab.io/opl-csv/files/",
      "openpowerlifting-latest.zip"
    ),
    format = "url"
  ),

  tar_target(
    open_powerlifting_snapshot,
    download_open_powerlifting_strength_sd(open_powerlifting_url)
  ),

  tar_target(
    strength_sd,
    open_powerlifting_snapshot |>
      select(squat_sd, bench_sd, deadlift_sd)
  ),
  
  # Grateful report ----
  tar_target(
    grateful_report,
    cite_packages(
      out.dir = "grateful/",
      cite.tidyverse = TRUE,
      out.format = "pdf"
    ),
    format = "file"
  )
  
)
