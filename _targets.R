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
tar_source("R/primary_model_diagnostics.R")
tar_source("R/primary_model_loo.R")
tar_source("R/primary_model_averaging.R")

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
    model_names_nunes,
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
    model_formulas_nunes,
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
      model_formulas_nunes,
      model_names_nunes,
      ~ fit_candidate_models(
        formula = .x,
        priors = steele_priors,
        data = nunes_data_prepared
      )
    ) |> set_names(model_names_nunes)
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
  tar_target(
    steele_strength_priors,
    get_steele_strength_priors(steele_data)
  ),
  
  tar_target(
    kirk_function_priors,
    get_kirk_function_priors()
  ),
  
  
  ## Fit all candidate models to new data ----

  # Note, as control median intake slightly different (1.1 g/kg/day) in the new data
  # The knot centred protein intake for lspline has been shifted to 0.5 to reflect the 1.6 g/kg/day
  tar_target(
    model_names_primary,
    c(
      "Intercept only model",
      "Linear model",
      "Linear spline (1.6g/kg/day) model",
      "Thin plate spline model",
      "Intercept only model (+ moderators)",
      "Linear model (+ moderators)",
      "Linear spline (1.6g/kg/day) model (+ moderators)",
      "Thin plate spline model (+ moderators)",
      "Linear model (+ random slopes)",
      "Linear spline (1.6g/kg/day) model (+ random slopes)",
      "Thin plate spline model (+ random smooths)",
      "Linear model (+ moderators & random slopes)",
      "Linear spline (1.6g/kg/day) model (+ moderators & random slopes)",
      "Thin plate spline model (+ moderators & random smooths)",
      "Linear model (+ RT interaction)",
      "Linear spline (1.6g/kg/day) model (+ RT interaction)",
      "Thin plate spline model (+ RT interaction smooth)",
      "Linear model (+ moderators & RT interaction)",
      "Linear spline (1.6g/kg/day) model (+ moderators & RT interaction)",
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
    model_formulas_primary,
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
                             s(protein_intake, study_factor, bs = 'fs', k = 4) +
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
                                  s(protein_intake, study_factor, bs = 'fs', k = 4) +
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
                      s(protein_intake, by=resistance_exercise_factor, m=1, bs='tp', k=4) +
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
                           s(protein_intake, by=resistance_exercise_factor, m=1, bs='tp', k=4) +
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
                             s(protein_intake, by=resistance_exercise_factor, m=1, bs='tp', k=4) +
                             s(protein_intake, study_factor, bs = 'fs', k = 4) +
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
                                  s(protein_intake, by=resistance_exercise_factor, m=1, bs='tp', k=4) +
                                  s(protein_intake, study_factor, bs = 'fs', k = 4) +
                               duration_centre + age_centre +
                             (1 | arm_no) +
                             (1 | outcome_no)"
    )
  ),

  # Store the earlier pilot-supported fitting decisions in one validated
  # registry. The completed model_fit_new_* branches are retained so that
  # unchanged models can be reused. Post-pilot remedies are fitted separately
  # by model_fit_remedied_* and merged by model index below.
  # Corrected LOO, stacking, posterior-mixture predictions, and contrasts are
  # connected downstream of the final audited LOO target.
  tar_target(
    primary_model_fit_specifications,
    get_primary_model_fit_specifications(
      model_names = model_names_primary,
      model_formulas = model_formulas_primary
    ),
    iteration = "list"
  ),

  tar_target(
    primary_model_fit_specification_table,
    dplyr::bind_rows(primary_model_fit_specifications)
  ),

  # Post-pilot replacements are kept separate from the completed primary
  # branches. This allows unchanged fits to be reused and limits refitting to
  # the 69 outcome-model combinations whose formula, priors, or controls have
  # changed.
  tar_target(
    primary_model_remedy_specifications,
    get_primary_model_remedy_specifications(
      primary_model_fit_specifications
    ),
    iteration = "list"
  ),

  tar_target(
    primary_model_remedy_specification_table,
    dplyr::bind_rows(primary_model_remedy_specifications)
  ),

  # Second-stage replacements are restricted to the 11 full models that
  # retained sampling problems after the first post-pilot remedies. Run the
  # model_fit_second_remedied_* targets explicitly so that the other 93
  # diagnostically acceptable full fits remain untouched.
  tar_target(
    primary_model_second_remedy_specifications,
    get_primary_model_second_remedy_specifications(
      primary_model_fit_specifications,
      primary_model_remedy_specifications
    ),
    iteration = "list"
  ),

  tar_target(
    primary_model_second_remedy_specification_table,
    dplyr::bind_rows(primary_model_second_remedy_specifications)
  ),

  # Final targeted function-model remedies are isolated from every completed
  # fit. Only function models 13 and 25 enter this registry.
  tar_target(
    primary_model_third_remedy_specifications,
    get_primary_model_third_remedy_specifications(
      primary_model_second_remedy_specifications
    ),
    iteration = "list"
  ),

  tar_target(
    primary_model_third_remedy_specification_table,
    dplyr::bind_rows(primary_model_third_remedy_specifications)
  ),

  ### Whole body lean mass ----
  tar_target(
    protein_data_complete_whole_body,
    {
      protein_data_complete |> filter(outcome_group == "whole body lean mass") |>
        # add study weights/sizes
        mutate(
          study_factor = factor(study_no),
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),

  tar_target(
    model_fit_specifications_whole_body,
    select_primary_model_fit_specifications(
      primary_model_fit_specifications,
      outcome = "whole_body"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_remedied_whole_body,
    select_primary_model_remedy_specifications(
      primary_model_remedy_specifications,
      outcome = "whole_body"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_second_remedied_whole_body,
    select_primary_model_second_remedy_specifications(
      primary_model_second_remedy_specifications,
      outcome = "whole_body"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_new_whole_body,
    fit_primary_candidate_model(
      specification = model_fit_specifications_whole_body,
      priors = steele_priors,
      data = protein_data_complete_whole_body
    ),
    pattern = map(model_fit_specifications_whole_body),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_remedied_whole_body,
    fit_primary_remedied_candidate_model(
      specification =
        model_fit_specifications_remedied_whole_body,
      priors = steele_priors,
      data = protein_data_complete_whole_body
    ),
    pattern = map(model_fit_specifications_remedied_whole_body),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_second_remedied_whole_body,
    fit_primary_second_remedied_candidate_model(
      specification =
        model_fit_specifications_second_remedied_whole_body,
      priors = steele_priors,
      data = protein_data_complete_whole_body
    ),
    pattern = map(
      model_fit_specifications_second_remedied_whole_body
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_final_whole_body,
    overlay_primary_candidate_fits(
      existing_candidates = merge_primary_candidate_fits(
        existing_candidates = model_fit_new_whole_body,
        replacement_candidates = model_fit_remedied_whole_body
      ),
      replacement_candidates = model_fit_second_remedied_whole_body,
      expected_replacements = 1L
    ),
    iteration = "list"
  ),

  tar_target(
    models_new_whole_body,
    assemble_primary_candidate_models(model_fit_final_whole_body)
  ),

  tar_target(
    model_diagnostics_new_whole_body,
    diagnose_primary_candidate_set(model_fit_final_whole_body),
    iteration = "list"
  ),

  tar_target(
    model_diagnostics_summary_new_whole_body,
    dplyr::bind_rows(
      purrr::map(model_diagnostics_new_whole_body, "summary")
    )
  ),

  tar_target(
    primary_psis_loo_initial_whole_body,
    compute_primary_candidate_psis_loo(model_fit_final_whole_body),
    pattern = map(model_fit_final_whole_body),
    iteration = "list",
    error = "continue"
  ),

  ### Hypertrophy ----
  tar_target(
    protein_data_complete_hypertrophy,
    {
      protein_data_complete |> filter(outcome_group == "direct muscle size/mass") |>
        # add study weights/sizes
        mutate(
          study_factor = factor(study_no),
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),

  tar_target(
    model_fit_specifications_hypertrophy,
    select_primary_model_fit_specifications(
      primary_model_fit_specifications,
      outcome = "hypertrophy"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_remedied_hypertrophy,
    select_primary_model_remedy_specifications(
      primary_model_remedy_specifications,
      outcome = "hypertrophy"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_new_hypertrophy,
    fit_primary_candidate_model(
      specification = model_fit_specifications_hypertrophy,
      priors = steele_priors,
      data = protein_data_complete_hypertrophy
    ),
    pattern = map(model_fit_specifications_hypertrophy),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_remedied_hypertrophy,
    fit_primary_remedied_candidate_model(
      specification =
        model_fit_specifications_remedied_hypertrophy,
      priors = steele_priors,
      data = protein_data_complete_hypertrophy
    ),
    pattern = map(model_fit_specifications_remedied_hypertrophy),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_final_hypertrophy,
    merge_primary_candidate_fits(
      existing_candidates = model_fit_new_hypertrophy,
      replacement_candidates = model_fit_remedied_hypertrophy
    ),
    iteration = "list"
  ),

  tar_target(
    models_new_hypertrophy,
    assemble_primary_candidate_models(model_fit_final_hypertrophy)
  ),

  tar_target(
    model_diagnostics_new_hypertrophy,
    diagnose_primary_candidate_set(model_fit_final_hypertrophy),
    iteration = "list"
  ),

  tar_target(
    model_diagnostics_summary_new_hypertrophy,
    dplyr::bind_rows(
      purrr::map(model_diagnostics_new_hypertrophy, "summary")
    )
  ),

  tar_target(
    primary_psis_loo_initial_hypertrophy,
    compute_primary_candidate_psis_loo(model_fit_final_hypertrophy),
    pattern = map(model_fit_final_hypertrophy),
    iteration = "list",
    error = "continue"
  ),

  ### Strength ----
  tar_target(
    protein_data_complete_strength,
    {
      protein_data_complete |> filter(outcome_group == "strength") |>
        # add study weights/sizes
        mutate(
          study_factor = factor(study_no),
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),

  tar_target(
    model_fit_specifications_strength,
    select_primary_model_fit_specifications(
      primary_model_fit_specifications,
      outcome = "strength"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_remedied_strength,
    select_primary_model_remedy_specifications(
      primary_model_remedy_specifications,
      outcome = "strength"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_second_remedied_strength,
    select_primary_model_second_remedy_specifications(
      primary_model_second_remedy_specifications,
      outcome = "strength"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_new_strength,
    fit_primary_candidate_model(
      specification = model_fit_specifications_strength,
      priors = steele_strength_priors,
      data = protein_data_complete_strength
    ),
    pattern = map(model_fit_specifications_strength),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_remedied_strength,
    fit_primary_remedied_candidate_model(
      specification = model_fit_specifications_remedied_strength,
      priors = steele_strength_priors,
      data = protein_data_complete_strength
    ),
    pattern = map(model_fit_specifications_remedied_strength),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_second_remedied_strength,
    fit_primary_second_remedied_candidate_model(
      specification =
        model_fit_specifications_second_remedied_strength,
      priors = steele_strength_priors,
      data = protein_data_complete_strength
    ),
    pattern = map(
      model_fit_specifications_second_remedied_strength
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_final_strength,
    overlay_primary_candidate_fits(
      existing_candidates = merge_primary_candidate_fits(
        existing_candidates = model_fit_new_strength,
        replacement_candidates = model_fit_remedied_strength
      ),
      replacement_candidates = model_fit_second_remedied_strength,
      expected_replacements = 5L
    ),
    iteration = "list"
  ),

  tar_target(
    models_new_strength,
    assemble_primary_candidate_models(model_fit_final_strength)
  ),

  tar_target(
    model_diagnostics_new_strength,
    diagnose_primary_candidate_set(model_fit_final_strength),
    iteration = "list"
  ),

  tar_target(
    model_diagnostics_summary_new_strength,
    dplyr::bind_rows(
      purrr::map(model_diagnostics_new_strength, "summary")
    )
  ),

  tar_target(
    primary_psis_loo_initial_strength,
    compute_primary_candidate_psis_loo(model_fit_final_strength),
    pattern = map(model_fit_final_strength),
    iteration = "list",
    error = "continue"
  ),

  ### Function ----
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
          study_factor = factor(study_no),
          wi = 1/sqrt(vi),
          size = 0.5 + 3.0 * (wi - min(wi, na.rm=TRUE))/(max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)))
    }
  ),
  
  tar_target(
    model_fit_specifications_function,
    select_primary_model_fit_specifications(
      primary_model_fit_specifications,
      outcome = "function"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_remedied_function,
    select_primary_model_remedy_specifications(
      primary_model_remedy_specifications,
      outcome = "function"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_second_remedied_function,
    select_primary_model_second_remedy_specifications(
      primary_model_second_remedy_specifications,
      outcome = "function"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_specifications_third_remedied_function,
    select_primary_model_third_remedy_specifications(
      primary_model_third_remedy_specifications,
      outcome = "function"
    ),
    iteration = "list"
  ),

  tar_target(
    model_fit_new_function,
    fit_primary_candidate_model(
      specification = model_fit_specifications_function,
      priors = kirk_function_priors,
      data = protein_data_complete_function
    ),
    pattern = map(model_fit_specifications_function),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_remedied_function,
    fit_primary_remedied_candidate_model(
      specification = model_fit_specifications_remedied_function,
      priors = kirk_function_priors,
      data = protein_data_complete_function
    ),
    pattern = map(model_fit_specifications_remedied_function),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_second_remedied_function,
    fit_primary_second_remedied_candidate_model(
      specification =
        model_fit_specifications_second_remedied_function,
      priors = kirk_function_priors,
      data = protein_data_complete_function
    ),
    pattern = map(
      model_fit_specifications_second_remedied_function
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_third_remedied_function,
    fit_primary_third_remedied_candidate_model(
      specification =
        model_fit_specifications_third_remedied_function,
      priors = kirk_function_priors,
      data = protein_data_complete_function
    ),
    pattern = map(
      model_fit_specifications_third_remedied_function
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    model_fit_final_function,
    overlay_primary_third_remedied_candidate_fits(
      existing_candidates = overlay_primary_candidate_fits(
        existing_candidates = merge_primary_candidate_fits(
          existing_candidates = model_fit_new_function,
          replacement_candidates = model_fit_remedied_function
        ),
        replacement_candidates = model_fit_second_remedied_function,
        expected_replacements = 5L
      ),
      replacement_candidates = model_fit_third_remedied_function
    ),
    iteration = "list"
  ),

  tar_target(
    models_new_function,
    assemble_primary_candidate_models(model_fit_final_function)
  ),

  tar_target(
    model_diagnostics_new_function,
    diagnose_primary_candidate_set(model_fit_final_function),
    iteration = "list"
  ),

  tar_target(
    model_diagnostics_summary_new_function,
    dplyr::bind_rows(
      purrr::map(model_diagnostics_new_function, "summary")
    )
  ),

  tar_target(
    primary_model_diagnostic_audit,
    compile_primary_model_diagnostic_audit(
      model_diagnostics_new_whole_body,
      model_diagnostics_new_hypertrophy,
      model_diagnostics_new_strength,
      model_diagnostics_new_function
    )
  ),

  tar_target(
    primary_model_diagnostic_audit_files,
    write_primary_model_diagnostic_audit(
      primary_model_diagnostic_audit
    ),
    format = "file"
  ),

  tar_target(
    primary_psis_loo_initial_function,
    compute_primary_candidate_psis_loo(model_fit_final_function),
    pattern = map(model_fit_final_function),
    iteration = "list",
    error = "continue"
  ),

  ## LOO audit ----
  ### Integrated observation-level PSIS-LOO ----
  # Each candidate remains a
  # separate branch so successful work is cached if another model requires
  # attention. The final audit is the decision gate for stacking.
  tar_target(
    primary_psis_loo_integrated_whole_body,
    compute_primary_candidate_integrated_psis_loo(
      fitted_candidate = model_fit_final_whole_body,
      initial_loo_result = primary_psis_loo_initial_whole_body
    ),
    pattern = map(
      model_fit_final_whole_body,
      primary_psis_loo_initial_whole_body
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    primary_psis_loo_integrated_hypertrophy,
    compute_primary_candidate_integrated_psis_loo(
      fitted_candidate = model_fit_final_hypertrophy,
      initial_loo_result = primary_psis_loo_initial_hypertrophy
    ),
    pattern = map(
      model_fit_final_hypertrophy,
      primary_psis_loo_initial_hypertrophy
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    primary_psis_loo_integrated_strength,
    compute_primary_candidate_integrated_psis_loo(
      fitted_candidate = model_fit_final_strength,
      initial_loo_result = primary_psis_loo_initial_strength
    ),
    pattern = map(
      model_fit_final_strength,
      primary_psis_loo_initial_strength
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    primary_psis_loo_integrated_function,
    compute_primary_candidate_integrated_psis_loo(
      fitted_candidate = model_fit_final_function,
      initial_loo_result = primary_psis_loo_initial_function
    ),
    pattern = map(
      model_fit_final_function,
      primary_psis_loo_initial_function
    ),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    primary_psis_loo_integrated_audit,
    combine_primary_integrated_psis_loo_results(
      primary_psis_loo_integrated_whole_body,
      primary_psis_loo_integrated_hypertrophy,
      primary_psis_loo_integrated_strength,
      primary_psis_loo_integrated_function
    )
  ),

  tar_target(
    primary_psis_loo_integrated_audit_files,
    write_primary_integrated_psis_loo_audit(
      primary_psis_loo_integrated_audit
    ),
    format = "file"
  ),

  ### Exact observation-level LOO  ----
  # restricted to the three remaining
  # integrated-PSIS warnings. Refit and scoring are separate dynamic stages
  # so each expensive refit is cached before the inexpensive exact score and
  # final audit are assembled.
  tar_target(
    primary_psis_loo_exact_jobs,
    select_primary_integrated_exact_loo_jobs(
      integrated_audit = primary_psis_loo_integrated_audit,
      fitted_whole_body = model_fit_final_whole_body,
      fitted_hypertrophy = model_fit_final_hypertrophy
    ),
    iteration = "list"
  ),

  tar_target(
    primary_psis_loo_exact_refits,
    fit_primary_integrated_exact_loo_case(
      primary_psis_loo_exact_jobs
    ),
    pattern = map(primary_psis_loo_exact_jobs),
    iteration = "list",
    error = "continue"
  ),

  tar_target(
    primary_psis_loo_exact_corrections,
    compute_primary_integrated_exact_loo_case(
      primary_psis_loo_exact_refits
    ),
    pattern = map(primary_psis_loo_exact_refits),
    iteration = "list",
    error = "continue"
  ),

  # Only the hypertrophy model-20 exact refit had a residual geometry issue
  # (one divergence in 24,000 retained draws). Escalate its sampler controls
  # without invalidating or rerunning either successful whole-body refit.
  tar_target(
    primary_psis_loo_exact_hypertrophy_model20_remedy_job,
    select_primary_exact_loo_geometry_remedy_job(
      primary_psis_loo_exact_jobs
    )
  ),

  tar_target(
    primary_psis_loo_exact_hypertrophy_model20_remedy_refit,
    fit_primary_integrated_exact_loo_geometry_remedy(
      primary_psis_loo_exact_hypertrophy_model20_remedy_job
    )
  ),

  tar_target(
    primary_psis_loo_exact_hypertrophy_model20_remedy_correction,
    compute_primary_integrated_exact_loo_case(
      primary_psis_loo_exact_hypertrophy_model20_remedy_refit
    )
  ),

  tar_target(
    primary_psis_loo_exact_corrections_final,
    replace_primary_exact_loo_geometry_remedy(
      exact_results = primary_psis_loo_exact_corrections,
      remedied_result =
        primary_psis_loo_exact_hypertrophy_model20_remedy_correction
    )
  ),

  tar_target(
    primary_psis_loo_final_audit,
    apply_primary_integrated_exact_loo_corrections(
      integrated_audit = primary_psis_loo_integrated_audit,
      exact_results = primary_psis_loo_exact_corrections_final
    )
  ),

  tar_target(
    primary_psis_loo_final_audit_files,
    write_primary_final_loo_audit(
      primary_psis_loo_final_audit
    ),
    format = "file"
  ),

  ## Corrected stacking and posterior-mixture analysis ----
  # The only LOO inputs below are the corrected objects in the final audit.
  # No downstream target recomputes ordinary or conditional LOO.
  tar_target(
    primary_stacking_weights,
    compute_primary_stacking_weights(
      primary_psis_loo_final_audit$results
    )
  ),

  tar_target(
    primary_stacking_weights_summary,
    summarise_primary_stacking_weights(primary_stacking_weights)
  ),

  tar_target(
    primary_stacking_weights_file,
    {
      output_dir <- here("_targets", "user", "primary_model_averaging")
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      output_file <- file.path(output_dir, "stacking_weights.csv")
      readr::write_csv(primary_stacking_weights_summary, output_file)
      output_file
    },
    format = "file"
  ),

  # Retain the established outcome-specific target names for the manuscript.
  tar_target(
    models_new_weights_whole_body,
    primary_stacking_weights$whole_body
  ),
  tar_target(
    models_new_weights_hypertrophy,
    primary_stacking_weights$hypertrophy
  ),
  tar_target(
    models_new_weights_strength,
    primary_stacking_weights$strength
  ),
  tar_target(
    models_new_weights_function,
    primary_stacking_weights[["function"]]
  ),

  # A model-averaged draw is a complete dose-response curve sampled from one
  # candidate model according to its corrected stacking weight.
  tar_target(
    primary_model_average_whole_body,
    build_primary_model_average(
      models = models_new_whole_body,
      weights = models_new_weights_whole_body,
      data = protein_data_complete_whole_body,
      seed = 1989L
    )
  ),
  tar_target(
    primary_model_average_hypertrophy,
    build_primary_model_average(
      models = models_new_hypertrophy,
      weights = models_new_weights_hypertrophy,
      data = protein_data_complete_hypertrophy,
      seed = 1990L
    )
  ),
  tar_target(
    primary_model_average_strength,
    build_primary_model_average(
      models = models_new_strength,
      weights = models_new_weights_strength,
      data = protein_data_complete_strength,
      seed = 1991L
    )
  ),
  tar_target(
    primary_model_average_function,
    build_primary_model_average(
      models = models_new_function,
      weights = models_new_weights_function,
      data = protein_data_complete_function,
      seed = 1992L
    )
  ),

  #### Predictions ----
  # Dose-response point estimates are posterior-mixture means and interval
  # endpoints are mixture quantiles, not weighted model-specific endpoints.
  tar_target(
    preds_new_whole_body,
    summarise_primary_dose_response(primary_model_average_whole_body)
  ),
  tar_target(
    preds_new_hypertrophy,
    summarise_primary_dose_response(primary_model_average_hypertrophy)
  ),
  tar_target(
    preds_new_strength,
    summarise_primary_dose_response(primary_model_average_strength)
  ),
  tar_target(
    preds_new_function,
    summarise_primary_dose_response(primary_model_average_function)
  ),

  #### Contrasts ----
  # Every current-to-proposed heatmap contrast uses paired values from the
  # same complete model-averaged curve.
  tar_target(
    primary_heatmap_contrasts_whole_body,
    compute_primary_heatmap_contrasts(
      primary_model_average_whole_body,
      protein_data_complete_whole_body
    )
  ),
  tar_target(
    primary_heatmap_contrasts_hypertrophy,
    compute_primary_heatmap_contrasts(
      primary_model_average_hypertrophy,
      protein_data_complete_hypertrophy
    )
  ),
  tar_target(
    primary_heatmap_contrasts_strength,
    compute_primary_heatmap_contrasts(
      primary_model_average_strength,
      protein_data_complete_strength
    )
  ),
  tar_target(
    primary_heatmap_contrasts_function,
    compute_primary_heatmap_contrasts(
      primary_model_average_function,
      protein_data_complete_function
    )
  ),

  # Publication-facing tidy heatmap data. Unsupported cells are retained and
  # explicitly flagged so the plot can mask them without discarding the three
  # posterior probabilities relative to the ROPE boundaries.
  tar_target(
    primary_heatmap_data_whole_body,
    primary_heatmap_contrasts_whole_body$contrasts
  ),
  tar_target(
    primary_heatmap_data_hypertrophy,
    primary_heatmap_contrasts_hypertrophy$contrasts
  ),
  tar_target(
    primary_heatmap_data_strength,
    primary_heatmap_contrasts_strength$contrasts
  ),
  tar_target(
    primary_heatmap_data_function,
    primary_heatmap_contrasts_function$contrasts
  ),

  tar_target(
    plot_new_whole_body,
    plot_primary_predictions_heatmap(
      protein_data_complete_whole_body,
      preds_new_whole_body,
      primary_heatmap_contrasts_whole_body,
      "whole body lean soft tissue change"
    )
  ),
  tar_target(
    plot_new_hypertrophy,
    plot_primary_predictions_heatmap(
      protein_data_complete_hypertrophy,
      preds_new_hypertrophy,
      primary_heatmap_contrasts_hypertrophy,
      "direct muscle mass/size change"
    )
  ),
  tar_target(
    plot_new_strength,
    plot_primary_predictions_heatmap(
      protein_data_complete_strength,
      preds_new_strength,
      primary_heatmap_contrasts_strength,
      "strength change"
    )
  ),
  tar_target(
    plot_new_function,
    plot_primary_predictions_heatmap(
      protein_data_complete_function,
      preds_new_function,
      primary_heatmap_contrasts_function,
      "function change"
    )
  ),

  tar_target(
    plot_new_whole_body_tiff,
    {
      output_file <- here("plots", "plot_new_whole_body.tiff")
      ggplot2::ggsave(output_file, plot_new_whole_body, device = "tiff",
                      dpi = 300, width = 8, height = 12)
      output_file
    },
    format = "file"
  ),
  tar_target(
    plot_new_hypertrophy_tiff,
    {
      output_file <- here("plots", "plot_new_hypertrophy.tiff")
      ggplot2::ggsave(output_file, plot_new_hypertrophy, device = "tiff",
                      dpi = 300, width = 8, height = 12)
      output_file
    },
    format = "file"
  ),
  tar_target(
    plot_new_strength_tiff,
    {
      output_file <- here("plots", "plot_new_strength.tiff")
      ggplot2::ggsave(output_file, plot_new_strength, device = "tiff",
                      dpi = 300, width = 8, height = 12)
      output_file
    },
    format = "file"
  ),
  tar_target(
    plot_new_function_tiff,
    {
      output_file <- here("plots", "plot_new_function.tiff")
      ggplot2::ggsave(output_file, plot_new_function, device = "tiff",
                      dpi = 300, width = 8, height = 12)
      output_file
    },
    format = "file"
  ),

  # Publication-ready deviation record for the preregistered brms primary
  # analysis. The rendered page reads the definitive sampling and LOO audit
  # artifacts directly, so changes to either final audit invalidate the HTML.
  tar_target(
    preregistration_deviations_inputs,
    list(
      diagnostics = primary_model_diagnostic_audit_files,
      loo = primary_psis_loo_final_audit_files,
      stacking = primary_stacking_weights_file
    )
  ),
  tar_target(
    preregistration_deviations_html,
    render_quarto_document(
      input = "supplement/preregistration-deviations.qmd",
      output = "supplement/preregistration-deviations.html",
      dependencies = unlist(
        preregistration_deviations_inputs,
        use.names = FALSE
      )
    ),
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
    open_powerlifting_snapshot_file,
    here(
      "data",
      "external",
      "openpowerlifting_strength_reference_2026-08-29.csv"
    ),
    format = "file"
  ),

  tar_target(
    open_powerlifting_snapshot,
    read_open_powerlifting_strength_snapshot(
      open_powerlifting_snapshot_file
    )
  ),

  tar_target(
    strength_sd,
    open_powerlifting_snapshot |>
      select(squat_sd, bench_sd, deadlift_sd)
  ),

  tar_target(
    primary_reference_values,
    make_primary_reference_values(
      protein_data = protein_data_prepared,
      lean_mass_sd = lean_mass_sd,
      muscle_thickness_sd = muscle_thickness_sd,
      open_powerlifting_snapshot = open_powerlifting_snapshot
    )
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
