# Rebuild the reported analyses

Start with the [clean-clone instructions](../REPRODUCIBILITY.md#start-from-a-clean-clone).
The main `_targets.R` graph has three exposure groups. The preregistered primary
analysis is the relative-intake model average. In the absolute groups,
“Primary” labels the base model family for that exploratory exposure.
The tables give the switch and target for each phase. Work through one family
and one row at a time. Run each base family before its sensitivity families.
Keep all accepted fits, final LOO
reductions, and target stores through measurement, GRADE, age/duration, and
training-status reporting.

For each row, substitute its switch, value, and target in this R pattern.
Use the same worker count in the environment and the command.

```r
switch <- "PROTEIN_UPDATE_2026_PRIMARY"
target <- "update_2026_primary_fit_file"
workers <- 2L
Sys.setenv(PROTEIN_MODEL_WORKERS = workers)
do.call(Sys.setenv, setNames(list("fit"), switch))
targets::tar_make_future(names = tidyselect::all_of(target), workers = workers)
Sys.unsetenv(switch)
```

For a row with several targets, pass them as a character vector to
`tidyselect::all_of()`. Complete the family-specific preflight, then inspect
its diagnostic and LOO audits after each phase. The file targets that assemble
audits and manifests run as dependencies of the next phase. An audit failure
calls for reviewing the new run and its diagnostics before proceeding. The
recorded remedies are the ones used for the published results.

The [2026 runbook](../data/update_2026/README.md#local-primary-refit-runbook)
explains the relative primary review gates and the total-intake common-sample
refit. Use the stage tables below for families whose historical runbook does
not list all phases. The existing guides describe the statistical checks and
reporting outputs. The initial preparation for absolute exposures is
`update_2026_absolute_exposure_audit_files`.

## Relative total protein (g/kg/day)

### Primary

Start with `update_2026_primary_input_preflight_file`. After the phases below, build `update_2026_primary_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_PRIMARY` | `fit` | `update_2026_primary_fit_file` |
| `PROTEIN_UPDATE_2026_PRIMARY_REMEDY` | `fit` | `update_2026_primary_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_PRIMARY_LOO` | `run` | `update_2026_primary_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_PRIMARY_EXACT` | `fit` | `update_2026_primary_exact_loo_file` |
| `PROTEIN_UPDATE_2026_PRIMARY_REPORTING` | `run` | `update_2026_primary_model_average_file` |

### Mundlak

Start with `update_2026_mundlak_input_preflight_file`. After the phases below, build `update_2026_mundlak_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_MUNDLAK` | `fit` | `update_2026_mundlak_fit_file` |
| `PROTEIN_UPDATE_2026_MUNDLAK_REMEDY` | `fit` | `update_2026_mundlak_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_MUNDLAK_SECOND_REMEDY` | `fit` | `update_2026_mundlak_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_MUNDLAK_LOO` | `run` | `update_2026_mundlak_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_MUNDLAK_EXACT` | `fit` | `update_2026_mundlak_exact_file` |
| `PROTEIN_UPDATE_2026_MUNDLAK_REPORTING` | `run` | `update_2026_mundlak_prediction_draw_file` |

### Strict within-study

Start with `update_2026_within_input_preflight_file`. After the phases below, build `update_2026_within_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_WITHIN` | `fit` | `update_2026_within_fit_file` |
| `PROTEIN_UPDATE_2026_WITHIN_REMEDY` | `fit` | `update_2026_within_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_WITHIN_LOO` | `run` | `update_2026_within_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_WITHIN_EXACT` | `fit` | `update_2026_within_exact_file` |
| `PROTEIN_UPDATE_2026_WITHIN_EXACT_REMEDY` | `fit` | `update_2026_within_exact_remedy_file` |
| `PROTEIN_UPDATE_2026_WITHIN_REPORTING` | `run` | `update_2026_within_prediction_draw_file` |

### Achieved intake

Start with `update_2026_achieved_input_preflight_file`. After the phases below, build `update_2026_achieved_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_ACHIEVED` | `fit` | `update_2026_achieved_fit_file` |
| `PROTEIN_UPDATE_2026_ACHIEVED_REMEDY` | `fit` | `update_2026_achieved_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_ACHIEVED_LOO` | `run` | `update_2026_achieved_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_ACHIEVED_EXACT` | `fit` | `update_2026_achieved_exact_file` |
| `PROTEIN_UPDATE_2026_ACHIEVED_EXACT_REMEDY` | `fit` | `update_2026_achieved_exact_remedy_file` |
| `PROTEIN_UPDATE_2026_ACHIEVED_REPORTING` | `run` | `update_2026_achieved_prediction_draw_file` |

### Assignment-derived intake

Start with `update_2026_assignment_input_preflight_file`. After the phases below, build `update_2026_assignment_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_ASSIGNMENT` | `fit` | `update_2026_assignment_fit_file` |
| `PROTEIN_UPDATE_2026_ASSIGNMENT_REMEDY` | `fit` | `update_2026_assignment_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_ASSIGNMENT_LOO` | `run` | `update_2026_assignment_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_ASSIGNMENT_EXACT` | `fit` | `update_2026_assignment_exact_file` |
| `PROTEIN_UPDATE_2026_ASSIGNMENT_EXACT_REMEDY` | `fit` | `update_2026_assignment_exact_remedy_file` |
| `PROTEIN_UPDATE_2026_ASSIGNMENT_REPORTING` | `run` | `update_2026_assignment_prediction_draw_file` |

### Errors in variables

Start with `update_2026_eiv_input_preflight_file`. After the phases below, build `update_2026_eiv_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_EIV` | `fit` | `update_2026_eiv_fit_file` |
| `PROTEIN_UPDATE_2026_EIV_REMEDY` | `fit` | `update_2026_eiv_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_EIV_LOO` | `run` | `update_2026_eiv_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_EIV_EXACT` | `fit` | `update_2026_eiv_exact_file` |
| `PROTEIN_UPDATE_2026_EIV_REPORTING` | `run` | `update_2026_eiv_prediction_draw_file` |

### Residual correlation

Start with `update_2026_rescor_input_preflight_file`. After the phases below, build `update_2026_rescor_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_RESCOR` | `fit` | `update_2026_rescor_fit_file` |
| `PROTEIN_UPDATE_2026_RESCOR_REMEDY` | `fit` | `update_2026_rescor_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_RESCOR_SECOND_REMEDY` | `fit` | `update_2026_rescor_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_RESCOR_LOO` | `run` | `update_2026_rescor_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_RESCOR_OBSERVATION_LOO` | `run` | `update_2026_rescor_observation_loo_file` |
| `PROTEIN_UPDATE_2026_RESCOR_OBSERVATION_EXACT` | `fit` | `update_2026_rescor_observation_exact_result_file` |
| `PROTEIN_UPDATE_2026_RESCOR_OBSERVATION_EXACT_REMEDY` | `fit` | `update_2026_rescor_observation_exact_remedy_result_file` |
| `PROTEIN_UPDATE_2026_RESCOR_REPORTING` | `run` | `update_2026_rescor_prediction_draw_file` |

### Protein manipulation strategy

Start with `update_2026_strategy_input_preflight_file`. After the phases below, build `update_2026_strategy_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_STRATEGY` | `fit` | `update_2026_strategy_fit_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_REMEDY` | `fit` | `update_2026_strategy_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_LOO` | `run` | `update_2026_strategy_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_EXACT` | `fit` | `update_2026_strategy_exact_fit_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_REPORTING` | `run` | `update_2026_strategy_prediction_draw_file` |

### Strategy within study

Start with `update_2026_strategy_within_input_preflight_file`. After the phases below, build `update_2026_strategy_within_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_STRATEGY_WITHIN` | `fit` | `update_2026_strategy_within_fit_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REMEDY` | `fit` | `update_2026_strategy_within_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_WITHIN_LOO` | `run` | `update_2026_strategy_within_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT` | `fit` | `update_2026_strategy_within_exact_fit_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT_REMEDY` | `fit` | `update_2026_strategy_within_exact_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REPORTING` | `run` | `update_2026_strategy_within_prediction_draw_file` |


## Absolute total protein (g/day)

### Primary

Start with `update_2026_total_primary_input_preflight_file`. After the phases below, build `update_2026_total_primary_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_PRIMARY` | `fit` | `update_2026_total_primary_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY` | `fit` | `update_2026_total_primary_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_PRIMARY_LOO` | `run` | `update_2026_total_primary_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE` | `fit` | `update_2026_total_primary_common_sample_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_REMEDY` | `fit` | `update_2026_total_primary_common_sample_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_LOO` | `run` | `update_2026_total_primary_common_sample_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REPORTING` | `run` | `update_2026_total_primary_prediction_draw_file` |

### Mundlak

Start with `update_2026_total_mundlak_input_preflight_file`. After the phases below, build `update_2026_total_mundlak_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_MUNDLAK` | `fit` | `update_2026_total_mundlak_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_REMEDY` | `fit` | `update_2026_total_mundlak_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_LOO` | `run` | `update_2026_total_mundlak_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_EXACT` | `fit` | `update_2026_total_mundlak_exact_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_CORRECTION` | `run` | `update_2026_total_mundlak_exact_correction_file` |
| `PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_REPORTING` | `run` | `update_2026_total_mundlak_prediction_draw_file` |

### Strict within-study

Start with `update_2026_total_within_input_preflight_file`. After the phases below, build `update_2026_total_within_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_WITHIN` | `fit` | `update_2026_total_within_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_WITHIN_REMEDY` | `fit` | `update_2026_total_within_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_WITHIN_LOO` | `run` | `update_2026_total_within_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_WITHIN_EXACT` | `fit` | `update_2026_total_within_exact_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_WITHIN_CORRECTION` | `run` | `update_2026_total_within_exact_correction_file` |
| `PROTEIN_UPDATE_2026_TOTAL_WITHIN_REPORTING` | `run` | `update_2026_total_within_prediction_draw_file` |

### Achieved intake

Start with `update_2026_total_achieved_input_preflight_file`. After the phases below, build `update_2026_total_achieved_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_ACHIEVED` | `fit` | `update_2026_total_achieved_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ACHIEVED_REMEDY` | `fit` | `update_2026_total_achieved_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ACHIEVED_LOO` | `run` | `update_2026_total_achieved_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ACHIEVED_EXACT` | `fit` | `update_2026_total_achieved_exact_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ACHIEVED_CORRECTION` | `run` | `update_2026_total_achieved_exact_correction_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ACHIEVED_REPORTING` | `run` | `update_2026_total_achieved_prediction_draw_file` |

### Assignment-derived intake

Start with `update_2026_total_assignment_input_preflight_file`. After the phases below, build `update_2026_total_assignment_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT` | `fit` | `update_2026_total_assignment_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT_REMEDY` | `fit` | `update_2026_total_assignment_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT_SECOND_REMEDY` | `fit` | `update_2026_total_assignment_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT_THIRD_REMEDY` | `fit` | `update_2026_total_assignment_third_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT_LOO` | `run` | `update_2026_total_assignment_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT_EXACT` | `fit` | `update_2026_total_assignment_exact_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT_CORRECTION` | `run` | `update_2026_total_assignment_exact_correction_file` |
| `PROTEIN_UPDATE_2026_TOTAL_ASSIGNMENT_REPORTING` | `run` | `update_2026_total_assignment_prediction_draw_file` |

### Errors in variables

Start with `update_2026_total_eiv_input_preflight_file`. After the phases below, build `update_2026_total_eiv_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_EIV` | `fit` | `update_2026_total_eiv_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_EIV_REMEDY` | `fit` | `update_2026_total_eiv_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_EIV_LOO` | `run` | `update_2026_total_eiv_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_EIV_EXACT` | `fit` | `update_2026_total_eiv_exact_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_EIV_REPORTING` | `run` | `update_2026_total_eiv_prediction_draw_file` |

### Residual correlation

Start with `update_2026_total_rescor_input_preflight_file`. After the phases below, build `update_2026_total_rescor_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR` | `fit` | `update_2026_total_rescor_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_REMEDY` | `fit` | `update_2026_total_rescor_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_SECOND_REMEDY` | `fit` | `update_2026_total_rescor_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_THIRD_REMEDY` | `fit` | `update_2026_total_rescor_third_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_LOO` | `run` | `update_2026_total_rescor_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_EXACT` | `fit` | `update_2026_total_rescor_exact_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_EXACT_REMEDY` | `fit` | `update_2026_total_rescor_exact_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_EXACT_SECOND_REMEDY` | `fit` | `update_2026_total_rescor_exact_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_EXACT_THIRD_REMEDY` | `fit` | `update_2026_total_rescor_exact_third_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_CORRECTION` | `run` | `update_2026_total_rescor_exact_correction_file` |
| `PROTEIN_UPDATE_2026_TOTAL_RESCOR_REPORTING` | `run` | `update_2026_total_rescor_prediction_draw_file` |

### Protein manipulation strategy

Start with `update_2026_total_strategy_input_preflight_file`. After the phases below, build `update_2026_total_strategy_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY` | `fit` | `update_2026_total_strategy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_REMEDY` | `fit` | `update_2026_total_strategy_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_LOO` | `run` | `update_2026_total_strategy_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_COMMON_SAMPLE` | `fit` | `update_2026_total_strategy_common_sample_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_COMMON_SAMPLE_LOO` | `run` | `update_2026_total_strategy_common_sample_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_REPORTING` | `run` | `update_2026_total_strategy_prediction_draw_file` |

### Strategy within study

Start with `update_2026_total_strategy_within_input_preflight_file`. After the phases below, build `update_2026_total_strategy_within_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_WITHIN` | `fit` | `update_2026_total_strategy_within_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_WITHIN_REMEDY` | `fit` | `update_2026_total_strategy_within_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_WITHIN_LOO` | `run` | `update_2026_total_strategy_within_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_WITHIN_EXACT` | `fit` | `update_2026_total_strategy_within_exact_fit_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_WITHIN_CORRECTION` | `run` | `update_2026_total_strategy_within_correction_file` |
| `PROTEIN_UPDATE_2026_TOTAL_STRATEGY_WITHIN_REPORTING` | `run` | `update_2026_total_strategy_within_prediction_draw_file` |


## Prescribed supplementation (g/day)

### Primary

Start with `update_2026_supplement_primary_input_preflight_file`. After the phases below, build `update_2026_supplement_primary_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_PRIMARY` | `fit` | `update_2026_supplement_primary_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_PRIMARY_REMEDY` | `fit` | `update_2026_supplement_primary_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_PRIMARY_LOO` | `run` | `update_2026_supplement_primary_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_PRIMARY_COMMON_SAMPLE` | `fit` | `update_2026_supplement_primary_common_sample_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_PRIMARY_COMMON_SAMPLE_LOO` | `run` | `update_2026_supplement_primary_common_sample_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_PRIMARY_EXACT` | `fit` | `update_2026_supplement_primary_exact_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_PRIMARY_REPORTING` | `run` | `update_2026_supplement_primary_prediction_draw_file` |

### Mundlak

Start with `update_2026_supplement_mundlak_input_preflight_file`. After the phases below, build `update_2026_supplement_mundlak_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_MUNDLAK` | `fit` | `update_2026_supplement_mundlak_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_MUNDLAK_REMEDY` | `fit` | `update_2026_supplement_mundlak_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_MUNDLAK_LOO` | `run` | `update_2026_supplement_mundlak_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_MUNDLAK_EXACT` | `fit` | `update_2026_supplement_mundlak_exact_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_MUNDLAK_CORRECTION` | `run` | `update_2026_supplement_mundlak_exact_correction_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_MUNDLAK_REPORTING` | `run` | `update_2026_supplement_mundlak_prediction_draw_file` |

### Strict within-study

Start with `update_2026_supplement_within_input_preflight_file`. After the phases below, build `update_2026_supplement_within_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN` | `fit` | `update_2026_supplement_within_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN_REMEDY` | `fit` | `update_2026_supplement_within_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN_LOO` | `run` | `update_2026_supplement_within_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN_EXACT` | `fit` | `update_2026_supplement_within_exact_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN_EXACT_REMEDY` | `fit` | `update_2026_supplement_within_exact_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN_EXACT_SECOND_REMEDY` | `fit` | `update_2026_supplement_within_exact_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN_CORRECTION` | `run` | `update_2026_supplement_within_exact_correction_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_WITHIN_REPORTING` | `run` | `update_2026_supplement_within_prediction_draw_file` |

### Residual correlation

Start with `update_2026_supplement_rescor_input_preflight_file`. After the phases below, build `update_2026_supplement_rescor_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR` | `fit` | `update_2026_supplement_rescor_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_REMEDY` | `fit` | `update_2026_supplement_rescor_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_SECOND_REMEDY` | `fit` | `update_2026_supplement_rescor_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_THIRD_REMEDY` | `fit` | `update_2026_supplement_rescor_third_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_LOO` | `run` | `update_2026_supplement_rescor_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_EXACT` | `fit` | `update_2026_supplement_rescor_exact_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_EXACT_REMEDY` | `fit` | `update_2026_supplement_rescor_exact_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_EXACT_SECOND_REMEDY` | `fit` | `update_2026_supplement_rescor_exact_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_EXACT_THIRD_REMEDY` | `fit` | `update_2026_supplement_rescor_exact_third_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_CORRECTION` | `run` | `update_2026_supplement_rescor_exact_correction_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_RESCOR_REPORTING` | `run` | `update_2026_supplement_rescor_prediction_draw_file` |

### Protein manipulation strategy

Start with `update_2026_supplement_strategy_input_preflight_file`. After the phases below, build `update_2026_supplement_strategy_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY` | `fit` | `update_2026_supplement_strategy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_REMEDY` | `fit` | `update_2026_supplement_strategy_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_LOO` | `run` | `update_2026_supplement_strategy_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_COMMON_SAMPLE` | `fit` | `update_2026_supplement_strategy_common_sample_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_COMMON_SAMPLE_LOO` | `run` | `update_2026_supplement_strategy_common_sample_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_EXACT` | `fit` | `update_2026_supplement_strategy_exact_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_REPORTING` | `run` | `update_2026_supplement_strategy_prediction_draw_file` |

### Strategy within study

Start with `update_2026_supplement_strategy_within_input_preflight_file`. After the phases below, build `update_2026_supplement_strategy_within_verification_files`.

| Switch | Value | Target |
| --- | --- | --- |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN` | `fit` | `update_2026_supplement_strategy_within_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_REMEDY` | `fit` | `update_2026_supplement_strategy_within_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_REMEDY_2` | `fit` | `update_2026_supplement_strategy_within_second_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_STRUCTURAL` | `fit` | `update_2026_supplement_strategy_within_structural_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_LOO` | `run` | `update_2026_supplement_strategy_within_integrated_loo_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_EXACT` | `fit` | `update_2026_supplement_strategy_within_exact_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_EXACT_REMEDY` | `fit` | `update_2026_supplement_strategy_within_exact_remedy_fit_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_CORRECTION` | `run` | `update_2026_supplement_strategy_within_exact_correction_file` |
| `PROTEIN_UPDATE_2026_SUPPLEMENT_STRATEGY_WITHIN_REPORTING` | `run` | `update_2026_supplement_strategy_within_prediction_draw_file` |

## Publication bias and GRADE

1. Build `update_2026_publication_bias_robma_preflight_file` to check JAGS and
   the pinned RoBMA package. Set `PROTEIN_PUBLICATION_BIAS_ROBMA=fit` and build
   `update_2026_publication_bias_fit_file`. Clear the switch and build
   `update_2026_publication_bias_reporting_bundle_files`,
   `update_2026_publication_bias_extended_data_files`,
   `update_2026_publication_bias_funnel_plot_files`,
   `update_2026_publication_bias_all_contrast_model_plot_files` and
   `update_2026_publication_bias_outcome_model_estimate_files`.
2. With the accepted relative primary fits, corrected LOO and reporting
   complete, set `PROTEIN_UPDATE_2026_GRADE=run`. Build
   `update_2026_grade_preflight_file`, then `update_2026_grade_analysis_files`.
   Clear the switch. See [GRADE inputs](../data/GRADE/README.md).

## Measurement, age and duration, and training status

| Analysis | Prerequisites from the main graph | Stages |
| --- | --- | --- |
| Both muscle-measurement analyses | `update_2026_primary_verification_files` and the reference-prior inputs | [Measurement stages](measurement-analyses.md#run-in-reviewed-stages). Run through `report` for both variants. |
| Age and intervention duration | Relative, total-intake, and supplementation primary verification bundles and retained fits/LOO, plus both measurement bundles | [Age/duration stages](age-duration-interactions.md). Start with `Rscript sensitivity/run_age_duration_batch.R prepare` and retain the freshly built total-intake parents. |
| Reported training status | Same three parent verification bundles, both measurement bundles, and the final age/duration reports | [Training-status stages](training-status.md), including the [marginal exact-refit pilot](training-status.md#11-november-marginal-exact-refit-pilot). |

## Rebuild the reporting interface

In the main graph, build `systematic_review_audit_files`,
`update_2026_primary_publication_figure_files`,
`update_2026_primary_mundlak_within_figure_files`, and
`update_2026_exposure_comparison_figure_files` after their parent reports.
The measurement, age/duration, and training-status report stages write their
own outputs. From the project root, refresh the portable tables, regenerate
the PRISMA includes, and render the publication package:

```sh
Rscript scripts/materialize_preprint_reporting.R
Rscript supplement/build-prisma-includes.R
Rscript R/rob2_final_figures.R
python -m pip install -r scripts/requirements-release.txt
python scripts/release_bundle.py outputs/preprint-rebuilt --render
```

Compare the rebuilt analytical tables and diagnostic dispositions with the
published interfaces. Sampler and numerical-library differences may affect
Monte Carlo results. Preserve the historical release files for exact
comparison.
