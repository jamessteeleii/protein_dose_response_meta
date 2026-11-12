# Sensitivity-analysis pilots

The separate post-registration [age/duration interaction workflow](age-duration-interactions.md)
and [RT training-status workflow](training-status.md) have their own target
scripts, stores, phase switches and stage runners. They are not dispatched by
the main `_targets.R` pipeline or the pilot commands below. See the
[repository workflow overview](../README.md#separate-ageduration-and-training-status-workflows)
for the rationale, prerequisite order, storage locations and standalone
supplement renders. Keep their stores separate and run each stage in its
documented clean process.

The two muscle-measurement analyses are documented separately in
[measurement-analyses.md](measurement-analyses.md). They exclude trunk and
skinfold estimates for the registered-method rerun, then restrict to
MRI/CT/ultrasound/biopsy for the direct-method sensitivity. Both remain
isolated from the other extended model families and are now complete.

This directory documents the opt-in pilot targets added for the supplementary
analyses. The pilot code is intentionally isolated from the primary model
fits. With the default environment (`PROTEIN_SENSITIVITY_PILOTS=none`), no
sensitivity target compiles Stan or samples a model.

All documented final analysis families are now complete. The full suite used
2,662 distinct candidate models but 4,266 completed fitting jobs after
diagnostic remedies, common-sample replacements, recovery runs, and exact
case-deletion refits. Recorded fitting alone exceeded 1,659 model-job hours.
Read the root [reproducibility and fitting-time audit](../REPRODUCIBILITY.md)
before enabling any historical phase. The commands below remain a staged
provenance record; they must not be executed as one batch.

## Storage design for remaining full analyses

The completed primary analysis and the first full Mundlak expansion predate
the storage design below. They collect fitted objects into list targets, which
causes the same posterior draws to be serialised in branch, collected-fit, and
ensemble targets. After their reporting bundles are frozen, those redundant
objects can be removed with `prune_targets.R`.

Every remaining full sensitivity family must instead follow these invariants:

1. **One file per fitted model.** A dynamic fit branch writes one persistent,
   ignored fit file and returns only its path as a `format = "file"` target.
   The full fit object is never returned as an ordinary `qs` target.
2. **Path-only ensemble manifests.** Remedy selection produces a small table
   containing outcome, model index, stage, accepted fit path, file size, and
   checksum. It must not collect fitted objects into `*_fits`, `*_new_fits`,
   or `*_ensemble_fits` lists.
3. **Branch-local reduction.** Diagnostics, integrated LOO contributions,
   posterior predictions, and exact-refit preparation load one accepted fit,
   calculate the required reduced object, and release the fit before the next
   branch. Aggregation occurs only after those reductions.
4. **No fit objects in job records.** Integrated- and exact-LOO job targets
   carry a fit path plus row/model metadata. Embedding a fit in each job would
   serialise it again.
5. **Average reduced predictions, not models.** Model averaging combines
   outcome-specific prediction-draw matrices and weights. No downstream
   target needs the complete ensemble resident or serialised at once.
6. **Freeze, verify, then prune.** The final family target writes a
   checksummed bundle under `results/verification/<family>/`. Only after that
   bundle passes its checksum audit may the single-copy fit directory and any
   temporary prediction draws be deleted.
7. **Primary-analysis parity.** Sampling diagnostics use the same unchanged
   gates, remedies make the smallest candidate-specific change, and model
   comparison follows the same conditional audit, integrated
   observation-level PSIS, exact case deletion above the Pareto-k threshold,
   outcome-specific stacking, and posterior-mixture reporting sequence used
   by the final primary analysis. A sensitivity-specific departure must be
   explicit and documented.
8. **Performance gate before fan-out.** Run and time one simple and one complex
   branch before dispatching a full family. LOO should use the same efficient
   matrix calculation as the primary workflow when the harmonised draw-by-row
   matrix is comfortably bounded. Do not use the slower observation-at-a-time
   `brms::loo(pointwise = TRUE)` route merely to obtain pointwise diagnostics;
   the ordinary LOO result already contains those diagnostics. Prevent hidden
   nested parallelism inside branch-local reductions, and parallelise at only
   one level.

This makes peak storage approximately the sum of the unique accepted fits,
rather than that sum multiplied by several collection and job targets. It
also keeps peak memory close to one fitted model plus its temporary reduction.
The fit-file directories are computational cache, not release artifacts; they
must remain ignored by Git. If complete fits are retained after publication,
archive them externally with the ensemble manifest and checksums.

## Safe structural check

Run this first:

```r
targets::tar_make(sensitivity_pilot_structure_audit)
targets::tar_read(sensitivity_pilot_structure_audit)
```

This builds the final-specification registry and the pilot data/formulas only.
It does not compile or sample a model. The default pilot uses whole-body lean
mass and protein models 2, 3, 4, 13, and 26. These span linear,
segmented-linear, nonlinear, moderator/random-slope, and maximally structured
nonlinear forms. Models 1 and 5 are also fitted unchanged on the same rows so
the pilot ensemble retains both null families. Change
`sensitivity_pilot_config` in `_targets.R` to alter the outcome, protein-model
indices, chains, or iterations.

## Explicitly enable one family

Enable only the family being tested, then make its fit or summary target:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "mundlak")
targets::tar_make(sensitivity_mundlak_pilot_fit)
targets::tar_make(sensitivity_mundlak_pilot_diagnostics)

Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "within")
targets::tar_make(sensitivity_within_pilot_fit)
targets::tar_make(sensitivity_within_pilot_diagnostics)
targets::tar_read(sensitivity_within_pilot_diagnostics)

Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(sensitivity_exposure_pilot_diagnostics)
targets::tar_read(sensitivity_exposure_pilot_diagnostics) |>
  dplyr::select(
    exposure,
    model_index,
    divergent_transitions,
    maximum_treedepth_transitions,
    minimum_ebfmi,
    maximum_rhat,
    minimum_bulk_ess,
    minimum_tail_ess,
    maximum_mcse_sd_ratio,
    geometry_pass,
    convergence_pass,
    diagnostic_pass
  )

Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "rescor")
targets::tar_make(sensitivity_rescor_pilot_diagnostics)
targets::tar_make(sensitivity_rescor_pilot_parameter_summary)
targets::tar_read(sensitivity_rescor_pilot_diagnostics)
targets::tar_read(sensitivity_rescor_pilot_parameter_summary)

Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "eiv")
targets::tar_make(sensitivity_eiv_smoke_diagnostics)
targets::tar_make(sensitivity_eiv_smoke_parameter_summary)
targets::tar_read(sensitivity_eiv_smoke_diagnostics)
targets::tar_read(sensitivity_eiv_smoke_parameter_summary)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

## Protein addition versus substitution strategy

Two linked pilots treat protein addition and macronutrient substitution as
distinct study-level randomized intervention strategies. The first retains the
arm-based outcome model; the second uses strict matched protein-arm minus
comparator-arm contrasts. They are joint meta-regressions, not unrelated
subgroup fits. Every one of the 26 candidate forms can receive a strategy main
effect. Each protein-containing form also receives the matching strategy-by-
dose basis: linear, two segmented-linear columns, natural-spline columns, or a
strategy-difference thin-plate smooth in the arm model; and the corresponding
strategy-by-protein-contrast columns in the strict within-study model. The
original age, duration, resistance-training, random-effect, remedy, and prior
structures are retained. Neither pilot introduces a strategy-by-dose-by-
resistance-training term.

The human audit remains in
`outputs/01a07ac4-523d-74d2-927d-bf4fe7aa0475/protein_manipulation_strategy_first_pass_audit.xlsx`.
Its versionable machine input is
`data/protein_manipulation_strategy_audit.csv`. This CSV is joined only to
`sensitivity_primary_data_by_outcome`, after the primary analysis datasets
have already been prepared. Changing a classification therefore invalidates
only the strategy family, not the completed primary or earlier sensitivity
fits. The modelling data retain only the strategy code and factor; evidence
excerpts and other audit prose are not copied into every fitted model.
The final CSV contains the study identifiers, classification, confidence, and
supporting evidence needed for analysis and sharing. Completed review fields
and derived eligibility flags are intentionally omitted. The reader derives
the principal and strict eligibility flags each time the file is loaded.

The principal sample includes resolved Addition and Substitution studies. The
strict sample is prepared in parallel and requires a High-confidence
classification from the completed audit. Mixed trials are excluded from the
arm-level joint analysis:
their shared comparator cannot be assigned uniquely to both strategies without
duplicating observations. Protein-only comparisons and unclear studies are
also excluded. They remain visible in the audit CSV.

Before sampling, run both structure and readiness audits. These do not fit a
model:

```r
targets::tar_make(
  sensitivity_strategy_pilot_structure_audit,
  reporter = "balanced"
)
targets::tar_make(
  sensitivity_strategy_within_pilot_structure_audit,
  reporter = "balanced"
)
targets::tar_make(
  sensitivity_strategy_full_formula_audit,
  reporter = "balanced"
)

strategy_audit <- targets::tar_read(
  sensitivity_strategy_pilot_structure_audit
)

strategy_audit$ready_to_sample
strategy_audit$strategy_counts
strategy_audit$outcome_summary
strategy_audit$support
strategy_audit$blocking_studies
strategy_audit$formula_summary

strategy_within_audit <- targets::tar_read(
  sensitivity_strategy_within_pilot_structure_audit
)

strategy_within_audit$ready_to_sample
strategy_within_audit$contrast_rows
strategy_within_audit$study_counts
strategy_within_audit$support
strategy_within_audit$formula_summary

targets::tar_read(sensitivity_strategy_full_formula_audit)[c(
  "models", "outcomes", "all_formulas_pass"
)]
```

The completed audit contains 99 resolved Addition or Substitution studies, of
which 70 are High confidence. All 109 audit rows have been reviewed. The lean
final CSV therefore represents a confirmed audit; schema, study-linkage, and
classification-value failures still stop the sampling gate before any branch
starts.

The formula-only unit test is independent of the targets store and never
samples:

```r
source("sensitivity/test_strategy_formulas.R")
```

Once the shared audit gate is open, each whole-body pilot uses models 1, 2, 3,
4, 5, 13, and 26, mirroring the other pilot families. Run the arm-based pilot
first:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "strategy")

targets::tar_make(
  sensitivity_strategy_pilot_diagnostics,
  reporter = "balanced"
)
targets::tar_read(sensitivity_strategy_pilot_diagnostics)
```

If models 1--4 reproduce the initial pilot pattern of clean geometry but low
ESS or high relative MCSE, run the focused longer-sampling ensemble:

```r
targets::tar_make(
  sensitivity_strategy_remedy_diagnostics,
  reporter = "balanced"
)

targets::tar_read(sensitivity_strategy_remedy_diagnostics) |>
  purrr::map_dfr("summary") |>
  dplyr::select(
    model_index,
    divergent_transitions,
    maximum_treedepth_transitions,
    minimum_ebfmi,
    maximum_rhat,
    minimum_bulk_ess,
    minimum_tail_ess,
    maximum_mcse_sd_ratio,
    geometry_pass,
    convergence_pass,
    diagnostic_pass
  )
```

This refits only models 1--4 with 1,500 warmup and 3,000 sampling iterations.
Models 5, 13, and 26 are retained from the initial run. All subsequent arm-
strategy LOO and model averaging targets use this combined ensemble.

Then run the strict within-study strategy pilot separately:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "strategy_within")

targets::tar_make(
  sensitivity_strategy_within_pilot_diagnostics,
  reporter = "balanced"
)
targets::tar_read(sensitivity_strategy_within_pilot_diagnostics)
```

Only after the sampling diagnostics pass should model comparison and reporting
be requested. Conditional LOO is an audit comparator; the final weights use
integrated observation-level PSIS and exact case deletion for Pareto-k values
above 0.7, as in the primary analysis:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "strategy")

targets::tar_make(
  sensitivity_strategy_pilot_integrated_loo_diagnostics,
  reporter = "balanced"
)
targets::tar_read(
  sensitivity_strategy_pilot_integrated_loo_diagnostics
)

targets::tar_make(
  sensitivity_strategy_pilot_final_loo_audit,
  reporter = "balanced"
)
targets::tar_make(
  sensitivity_strategy_pilot_reporting,
  reporter = "balanced"
)

strategy_reporting <- targets::tar_read(
  sensitivity_strategy_pilot_reporting
)
strategy_reporting$dose_response
strategy_reporting$addition_minus_substitution
strategy_reporting$common_support

Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "strategy_within")

targets::tar_make(
  sensitivity_strategy_within_pilot_integrated_loo_diagnostics,
  reporter = "balanced"
)
targets::tar_make(
  sensitivity_strategy_within_pilot_final_loo_audit,
  reporter = "balanced"
)
targets::tar_make(
  sensitivity_strategy_within_pilot_reporting,
  reporter = "balanced"
)

strategy_within_reporting <- targets::tar_read(
  sensitivity_strategy_within_pilot_reporting
)
strategy_within_reporting$dose_response_contrasts
strategy_within_reporting$addition_minus_substitution
strategy_within_reporting$common_support

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Both reporting objects contain paired Addition and Substitution curves, their
draw-wise Addition-minus-Substitution curve, and observed/common-support flags.
The arm object reports predicted arm changes; the within object reports the
protein-arm minus comparator-arm effect with comparator intake held at 0.8
g/kg/day.

## Full protein-strategy expansions

Both 104-model sampling stages use the same one-fit-per-file and path-only
manifest design as the completed EIV and residual-correlation families. The
seven accepted whole-body pilot fits are archived once and reused, leaving 97
new fits per analysis. Run the arm-based strategy family first:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "strategy")

targets::tar_make(
  sensitivity_strategy_full_structure_audit,
  reporter = "balanced"
)
targets::tar_read(sensitivity_strategy_full_structure_audit)

targets::tar_make(
  sensitivity_strategy_full_diagnostic_audit,
  reporter = "balanced"
)
```

Each successful branch is retained if the run is interrupted. Inspect the
104-model audit before constructing any focused remedies or downstream LOO:

```r
strategy_full_audit <- targets::tar_read(
  sensitivity_strategy_full_diagnostic_audit
)
strategy_full_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]
strategy_full_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

The observed initial audit identified 12 candidates requiring remediation.
Refit only those candidates with 2,000 warmup and 4,000 retained iterations,
an adapt-delta floor of 0.999, and a maximum-treedepth floor of 15:

```r
targets::tar_make(
  sensitivity_strategy_full_final_diagnostic_audit,
  reporter = "balanced"
)

strategy_final_audit <- targets::tar_read(
  sensitivity_strategy_full_final_diagnostic_audit
)
strategy_final_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]
strategy_final_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

The 92 initially passing fit files are left untouched. Remedy files replace
their matching initial files only in
`sensitivity_strategy_full_accepted_fit_manifest`, and only after every final
sampling diagnostic passes.

After all 104 models pass, run the branch-local integrated LOO stage. This
uses up to 4,000 posterior draws per model, keeps computation single-core
inside each branch to avoid repeated Windows clusters, and writes only reduced
LOO objects:

```r
targets::tar_make(
  sensitivity_strategy_full_integrated_loo_audit,
  reporter = "balanced"
)

strategy_loo <- targets::tar_read(
  sensitivity_strategy_full_integrated_loo_audit
)
strategy_loo$outcome_summary
strategy_loo$exact_refits_required
strategy_loo$high_k_cases
```

Completed LOO branches are restart-safe and do not duplicate the accepted fit
files. The observed audit selected eight exact deletions: three whole-body and
five hypertrophy cases. Fit and diagnose only those cases next:

```r
targets::tar_make(
  sensitivity_strategy_full_exact_loo_diagnostic_audit,
  reporter = "balanced"
)

strategy_exact <- targets::tar_read(
  sensitivity_strategy_full_exact_loo_diagnostic_audit
)
strategy_exact[c(
  "exact_refits",
  "exact_refits_passing",
  "exact_refits_requiring_remediation",
  "all_exact_refit_diagnostics_pass"
)]
strategy_exact$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

After every exact refit passes, request
`sensitivity_strategy_full_final_loo_audit` and
`sensitivity_strategy_full_stacking_weights`. Exact corrections use the
arm-model integration path rather than the Mundlak basis reconstruction.

Once the final LOO audit is stacking-ready, freeze the completed arm analysis
before constructing reduced predictions:

```r
source("sensitivity/freeze_strategy_full_checkpoint.R")

targets::tar_make(
  sensitivity_strategy_verification_files,
  reporter = "balanced"
)

strategy_reporting <- targets::tar_read(
  sensitivity_strategy_full_reporting_audit
)
strategy_reporting$summary
strategy_reporting$all_reporting_checks_pass
```

The reporting stage allocates 4,000 posterior-mixture draws per outcome and
retains strategy-specific dose-response curves, increments from 0.8 g/kg/day,
dose-pair contrasts, support flags, and the direct draw-wise Addition-minus-
Substitution curve. The verification bundle is written to
`results/verification/strategy/`; it contains the reduced scientific results,
diagnostics, LOO summaries, weights, file checksums, and the direct strategy
contrast needed for inspection after the large fit cache is pruned.

After `all_reporting_checks_pass` is `TRUE`, preview the storage cleanup:

```sh
Rscript sensitivity/prune_targets.R --family=strategy
```

Review the printed targets, directories, and total size. If they are correct,
execute the same verified cleanup with:

```sh
Rscript sensitivity/prune_targets.R --family=strategy --execute
```

The strategy-within pilot and full-stage targets are deliberately excluded
from this family pattern, so pruning the arm-strategy analysis cannot remove
the inputs needed for the next expansion.

Only after the arm family is complete should the strategy-modified strict
within-study family be started:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "strategy_within")

targets::tar_make(
  sensitivity_strategy_within_full_structure_audit,
  reporter = "balanced"
)
targets::tar_read(sensitivity_strategy_within_full_structure_audit)

targets::tar_make(
  sensitivity_strategy_within_full_diagnostic_audit,
  reporter = "balanced"
)
```

The initial 104-model audit selected six focused remedies: hypertrophy models
10 and 22 for isolated divergences, and strength models 9, 10, 12, and 13 for
convergence or efficiency thresholds. Refit only those candidates with 2,000
warmup and 4,000 retained iterations, an `adapt_delta` floor of 0.999, and a
maximum-treedepth floor of 15:

```r
targets::tar_make(
  sensitivity_strategy_within_full_final_diagnostic_audit,
  reporter = "balanced"
)
```

The remaining 98 initial fits are retained unchanged. Integrated LOO begins
only after the final 104-model diagnostic audit passes. Run the restart-safe,
one-file-per-model integrated LOO stage with:

```r
targets::tar_make(
  sensitivity_strategy_within_full_integrated_loo_audit,
  reporter = "balanced"
)

strategy_within_loo <- targets::tar_read(
  sensitivity_strategy_within_full_integrated_loo_audit
)
strategy_within_loo$outcome_summary
strategy_within_loo$exact_refits_required
strategy_within_loo$high_k_cases
```

No exact deletion refits are launched until this audit has identified the
specific model-observation cases with integrated Pareto-k above 0.7.

The full integrated audit selected 67 exact observation-deletion refits: 31
for function, 25 for hypertrophy, 11 for whole-body lean mass, and none for
strength. Run the restart-safe exact stage to its diagnostic gate with:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "strategy_within")

targets::tar_make(
  sensitivity_strategy_within_full_exact_loo_diagnostic_audit,
  reporter = "balanced"
)

strategy_within_exact <- targets::tar_read(
  sensitivity_strategy_within_full_exact_loo_diagnostic_audit
)
strategy_within_exact[c(
  "exact_refits",
  "exact_refits_passing",
  "exact_refits_requiring_remediation",
  "all_exact_refit_diagnostics_pass"
)]
strategy_within_exact$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

Each refit uses 2,000 warmup and 6,000 retained iterations and is written to
its own file. Completed branches therefore survive an interruption. Corrected
LOO and stacking targets remain downstream of the diagnostic gate and must not
be used until every exact refit passes.

All 67 exact refits passed without remediation. After the corrected final LOO
audit is stacking-ready, freeze the compact checkpoint before constructing
reduced posterior predictions:

```r
source("sensitivity/freeze_strategy_within_full_checkpoint.R")

targets::tar_make(
  sensitivity_strategy_within_verification_files,
  reporter = "balanced"
)

strategy_within_reporting <- targets::tar_read(
  sensitivity_strategy_within_full_reporting_audit
)
strategy_within_reporting$summary
strategy_within_reporting$all_reporting_checks_pass
```

The resulting files preserve strategy-specific within-study contrast curves,
addition-minus-substitution contrasts, supported increments, pairwise dose
contrasts, model weights, sampling diagnostics, and exact-LOO provenance. Only
after every reporting check passes, preview and then execute the verified prune:

```sh
Rscript sensitivity/prune_targets.R --family=strategy_within
Rscript sensitivity/prune_targets.R --family=strategy_within --execute
```

## Full Mundlak expansion

After the whole-body pilot passed, the full Mundlak stage reuses its seven
accepted fits and independently checkpoints the remaining 97 candidates. The
result is the complete 26-model candidate set for each of whole-body lean mass,
hypertrophy, strength, and function. Pending fits use 1,500 warmup and 2,000
retained iterations per chain. Following the pilot remedy, model 4 receives the
regularised variance prior prospectively in the three new outcomes; all other
candidate-specific forms, covariates, priors, and primary-analysis remedies are
retained.

First run the inexpensive expansion audit:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "mundlak")
targets::tar_make(sensitivity_mundlak_full_structure_audit)
targets::tar_read(sensitivity_mundlak_full_structure_audit)
```

Then this single command can be left running unattended. Each model branch is
saved as soon as it completes, and an error in one branch does not prevent the
remaining branches from being attempted:

```r
targets::tar_make(
  sensitivity_mundlak_full_diagnostic_audit,
  reporter = "balanced"
)
```

After it finishes, inspect the sampling gate with:

```r
mundlak_audit <- targets::tar_read(
  sensitivity_mundlak_full_diagnostic_audit
)

mundlak_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]

mundlak_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass) |>
  dplyr::select(
    outcome,
    model_index,
    divergent_transitions,
    maximum_treedepth_transitions,
    minimum_ebfmi,
    maximum_rhat,
    minimum_bulk_ess,
    minimum_tail_ess,
    maximum_mcse_sd_ratio,
    status
  )

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

The completed expansion produced 96 sampling-diagnostic passes from 104 fits.
The remaining eight cases had no divergences and otherwise healthy geometry:
whole-body model 15; hypertrophy model 1; strength models 1, 5, and 6; and
function models 2, 5, and 15. Seven were marginal R-hat or variance-component
ESS misses. Strength model 6 additionally had two maximum-treedepth events.

The focused remedy retains every model, functional form, covariate, and prior.
All eight cases use 2,000 warmup and 4,000 retained iterations per chain;
strength model 6 also raises `max_treedepth` to 15. Remedy fits are written
directly as one persistent file per model, and the targets store only their
paths. This avoids creating another collected-fit or ensemble copy of the
posterior draws.

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "mundlak")
targets::tar_make(
  sensitivity_mundlak_full_final_diagnostic_audit,
  reporter = "balanced"
)

mundlak_final_audit <- targets::tar_read(
  sensitivity_mundlak_full_final_diagnostic_audit
)

mundlak_final_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]

mundlak_final_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass) |>
  dplyr::select(
    outcome,
    model_index,
    divergent_transitions,
    maximum_treedepth_transitions,
    minimum_ebfmi,
    maximum_rhat,
    minimum_bulk_ess,
    minimum_tail_ess,
    maximum_mcse_sd_ratio,
    diagnostic_pass
  )

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Do not prune the original full-fit targets at this point. Once the eight
remedies pass, materialise `sensitivity_mundlak_full_accepted_fit_manifest`.
That target makes one persistent, checksummed fit file for each accepted model
and replaces failed originals with their remedies in a path-only manifest.
This conversion is post-processing only, but temporarily requires roughly one
additional unique-fit-sized copy on disk. Run it as its own checkpoint:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "mundlak")
targets::tar_make(
  sensitivity_mundlak_full_accepted_fit_manifest,
  reporter = "balanced"
)

mundlak_fit_manifest <- targets::tar_read(
  sensitivity_mundlak_full_accepted_fit_manifest
)

mundlak_fit_manifest |>
  dplyr::count(outcome, stage)

sum(mundlak_fit_manifest$bytes) / 1024^3
```

The full integrated-LOO stage then reads one accepted fit file per dynamic
branch. Conditional and integrated PSIS are calculated while that single fit
is resident, and only their reduced LOO objects and pointwise diagnostics are
saved. All candidates use the same deterministic 4,000-draw LOO budget.

```r
targets::tar_make(
  sensitivity_mundlak_full_integrated_loo_audit,
  reporter = "balanced"
)

mundlak_loo_audit <- targets::tar_read(
  sensitivity_mundlak_full_integrated_loo_audit
)

mundlak_loo_audit$outcome_summary
mundlak_loo_audit$exact_refits_required
mundlak_loo_audit$high_k_cases

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

The completed integrated audit found 13 Pareto-k cases above 0.7 across 11
models. Eleven cases belong to nine whole-body models and two belong to two
hypertrophy models; strength and function require no exact deletion. No case
has Pareto k above 1. The exact stage follows the primary-analysis schedule of
2,000 warmup and 6,000 retained iterations per chain. Each deletion refit is
written directly to one persistent file, while its target record contains only
the accepted-fit path, integrated-LOO path, and scalar case metadata.

Run only the diagnostic gate first:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "mundlak")
targets::tar_make(
  sensitivity_mundlak_full_exact_loo_diagnostic_audit,
  reporter = "balanced"
)

mundlak_exact_audit <- targets::tar_read(
  sensitivity_mundlak_full_exact_loo_diagnostic_audit
)

mundlak_exact_audit[c(
  "exact_refits",
  "exact_refits_passing",
  "exact_refits_requiring_remediation",
  "all_exact_refit_diagnostics_pass"
)]

mundlak_exact_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Exact scoring, correction, and outcome-specific stacking remain downstream of
this gate. Only after those results, curve summaries, and a verified family
bundle are complete should `prune_targets.R` be used, first in dry-run mode.

After the corrected LOO audit is stacking-ready, freeze its lightweight hand-
off checkpoint once. This reads cached results only and does not run the
pipeline:

```r
source("sensitivity/freeze_mundlak_checkpoint.R")
```

The checkpoint contains the analysis data, specifications, accepted-fit paths,
reduced diagnostics, corrected LOO summaries, and stacking weights, but no
fitted objects, LOO objects, or posterior draws. Downstream reporting treats it
as an independent input. Consequently, changes to the opt-in environment flag
cannot accidentally send reporting back through the completed fitting or LOO
graph.

Then generate the full Mundlak reporting layer with sensitivity fitting
disabled. This stage does not sample or refit models. It allocates 4,000
complete posterior dose-response curves per outcome according to the exact-
corrected stacking weights, evaluates each selected fit in its own branch on
the same 0.01 g/kg/day grid used by the primary analysis, and writes only
reduced prediction files. It then produces absolute curves, paired increments
from 0.8 g/kg/day, observed-support flags, pairwise SESOI contrasts, and local-
equivalence summaries.

```r
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")

targets::tar_make(
  sensitivity_mundlak_full_reporting_audit,
  reporter = "balanced"
)

mundlak_reporting_audit <- targets::tar_read(
  sensitivity_mundlak_full_reporting_audit
)

mundlak_reporting_audit$summary
mundlak_reporting_audit$all_reporting_checks_pass
```

The implementation never constructs a 104-fit ensemble target. Each selected
fit is read once by a dynamic prediction branch, its allocated draws are
checkpointed to a small file, and only four outcome-level 4,000-curve mixtures
are assembled. Once the reporting audit passes, freeze the durable,
checksummed family bundle:

```r
targets::tar_make(
  sensitivity_mundlak_verification_files,
  reporter = "balanced"
)

targets::tar_read(sensitivity_mundlak_verification_files)
read.csv("results/verification/mundlak/manifest.csv")
```

The verification bundle deliberately excludes fitted model objects and draw
matrices. It retains the analysis data, candidate registry, sampling and exact
LOO diagnostics, final corrected LOO summaries, stacking weights, prediction
allocation, absolute and increment curves, exposure support, pairwise
contrasts, and local-equivalence tables. Do not prune until both the reporting
audit and this checksum manifest have been inspected.

The EIV smoke stage fits only model 2, the unmoderated linear candidate, using
the whole-body pilot data. It contains 233 outcome rows from 181 arms in 69
studies. If it passes, its saved fit is reused as model 2 when the other six
pilot candidates are fitted, so the smoke test is not duplicated.

The initial smoke fit was geometrically clean and met all strict thresholds
except the tail ESS for `sd_arm` (370 rather than 400). The focused remedy
keeps the model and priors unchanged and uses 1,500 warmup and 2,000 retained
iterations per chain. This remedied fit is the model-2 fit reused downstream:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "eiv")
targets::tar_make(sensitivity_eiv_smoke_remedy_diagnostics)
targets::tar_make(sensitivity_eiv_smoke_remedy_parameter_summary)
targets::tar_read(sensitivity_eiv_smoke_remedy_diagnostics)
targets::tar_read(sensitivity_eiv_smoke_remedy_parameter_summary)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

The first remedy produced an `sd_arm` tail ESS of 398.716, just below the
unchanged threshold of 400. To avoid a post-hoc threshold waiver, the second
remedy uses 2,000 warmup and 4,000 retained iterations per chain. This is the
model-2 fit reused in the complete EIV ensemble:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "eiv")
targets::tar_make(sensitivity_eiv_smoke_second_remedy_diagnostics)
targets::tar_make(sensitivity_eiv_smoke_second_remedy_parameter_summary)
targets::tar_read(sensitivity_eiv_smoke_second_remedy_diagnostics)
targets::tar_read(sensitivity_eiv_smoke_second_remedy_parameter_summary)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

### Focused EIV ensemble remedies

The initial remaining-candidate pilot had clean geometry throughout. Models
1, 3, 4, 5, and 13 first use a 1,500-warmup/2,000-retained remedy for marginal
R-hat or variance-component ESS misses. Models 1, 3, 4, and 13 pass at that
stage. Model 5 retains variance-component mixing limitations and joins the
more complex model 26 at the 2,000-warmup/4,000-retained stage. Models 1, 2,
3, 4, and 13 are not refitted again.

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "eiv")
targets::tar_make(sensitivity_eiv_remedy_diagnostics)
targets::tar_read(sensitivity_eiv_remedy_diagnostics)

targets::tar_make(sensitivity_eiv_second_remedy_diagnostics)
targets::tar_read(sensitivity_eiv_second_remedy_diagnostics)

targets::tar_make(sensitivity_eiv_final_diagnostics)
targets::tar_read(sensitivity_eiv_final_diagnostics)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

All seven final EIV pilot candidates now pass the unchanged geometry and
convergence gates. The final ensemble retains models 1, 2, 3, 4, 5, 13, and
26; the focused remedy fits are reused rather than sampled again.

Before calculating stacking weights, run the outcome-random-effects integrated
PSIS audit. This is post-processing of the saved fits and does not sample any
new model:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "eiv")
targets::tar_make(sensitivity_eiv_integrated_loo_diagnostics)
targets::tar_read(sensitivity_eiv_integrated_loo_diagnostics) |>
  dplyr::select(
    model_index,
    model_name,
    method,
    elpd_loo,
    se_elpd_loo,
    p_loo,
    maximum_pareto_k,
    n_pareto_k_gt_threshold,
    n_pareto_k_gt_1
  )
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

The integrated calculation uses 4,000 harmonised posterior draws per model,
retains each draw's latent true protein intake, and analytically integrates the
Gaussian study, arm, and outcome effects. The conditional `log_lik_y` result is
reported alongside it only as a diagnostic. Stacking remains blocked if the
integrated result contains any Pareto-k value above 0.7; such cases will be
given exact leave-one-out refits before weights are calculated.

In the completed whole-body pilot, integration leaves one marginal warning:
model 3 has one observation with Pareto k = 0.742 and no integrated case has
k above 1. Inspect and run that single exact refit with:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "eiv")
targets::tar_make(sensitivity_eiv_exact_loo_cases)
targets::tar_read(sensitivity_eiv_exact_loo_cases)

targets::tar_make(sensitivity_eiv_exact_loo_refit_diagnostics)
targets::tar_read(sensitivity_eiv_exact_loo_refit_diagnostics)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Unlike the integrated audit, the second command samples one leave-one-out
model. It retains the held-out arm's reported intake and intake SE, while
excluding only the selected outcome from the likelihood. The exact correction,
final stacking audit, weights, and averaged curve remain downstream of the
refit's diagnostic gate.

The hand-written Stan pilots save their CmdStan CSV files under
`_targets/user/sensitivity_cmdstan/`. This persistent location is required
because fitted `CmdStanMCMC` objects retain paths to those files; the default R
temporary directory can disappear before `{targets}` serialises or reloads the
fit.

Both hand-written Stan likelihoods retain the preregistered meta-analytic
sampling model: the reported outcome standard errors are treated as known and
no additional observation-level residual standard deviation is estimated.
Heterogeneity remains represented by the candidate's study, arm, and outcome
random effects. The residual-correlation pilot changes only the correlation of
the known sampling errors within its selected cluster. Its logit-normal prior
is centred at zero correlation after accounting for the cluster-size-dependent
lower bound required by an exchangeable correlation matrix.

### Full residual-correlation expansion

The residual-correlation sensitivity is expanded to the complete final
candidate registry: 26 models for each of whole-body lean mass, hypertrophy,
strength, and physical function. It retains the primary functional forms,
moderators, resistance-training interactions, random study slopes or smooths,
candidate-specific remedies, and prior scales. The sole substantive change is
that known sampling errors within the same study have a common estimated
exchangeable correlation. The two null candidates are retained so all 26
models enter outcome-specific stacking.

Two predictive audits are retained because they answer different questions.
Study-level LOO assesses transport to a completely new study and is kept as a
cluster-influence diagnostic. It proved intentionally severe in this data:
large recurring study clusters produced 166 high-k cases, including 100 above
1, so it is not used automatically for stacking.

The like-for-like comparison with the primary analysis instead uses a
correlated conditional observation-level likelihood. For every posterior draw
and outcome row, a standalone generated-quantities pass reconstructs the
study's marginal multivariate-normal covariance (including residual
correlation and integrated study, arm, and outcome heterogeneity) and evaluates
`p(y_i | y_-i, theta)`. No posterior model is refitted at this stage. Only the
compact LOO reduction is archived; temporary CmdStan output is deleted.

The no-sampling structure gate is:

```r
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
targets::tar_make(
  sensitivity_rescor_full_structure_audit,
  reporter = "balanced"
)
targets::tar_read(sensitivity_rescor_full_structure_audit)
```

After it passes, the initial 104-model sampling stage can be left running and
resumed safely after interruption:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "rescor")
targets::tar_make(
  sensitivity_rescor_full_diagnostic_audit,
  reporter = "balanced"
)

rescor_audit <- targets::tar_read(
  sensitivity_rescor_full_diagnostic_audit
)
rescor_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]
rescor_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Stop at that audit so any failures can receive the smallest necessary remedy.
After the sampling remedies pass, run the correlated conditional
observation-level LOO audit before starting any exact refits.

If isolated divergences remain after the first remedy, the second remedy
retains the model unchanged and uses 3,000 warmup, 8,000 retained iterations,
`adapt_delta = 0.9999`, and `max_treedepth = 16`. It is requested through
`sensitivity_rescor_full_complete_diagnostic_audit` and selects only candidates
that still fail the first-remedy audit.

If an isolated divergence persists, the final sampling-only remedy again
retains the entire candidate unchanged and uses a fresh seed, 4,000 warmup,
12,000 retained iterations, `adapt_delta = 0.99999`, and
`max_treedepth = 17`. It is requested through
`sensitivity_rescor_full_resolved_diagnostic_audit` and refits only candidates
that still fail the second remedy.

For function model 23, the single divergent draw persisted and repeatedly
occurred at the boundary of its weakly identified study-level random-smooth
correlation matrix. The targeted structural remedy retains the nonlinear
population curve, RT interaction smooth, study-specific random smooths, and
all coefficient and variance priors, but estimates the random-smooth
coefficients independently. This is the same class of minimal computational
remedy represented by `uncorrelated_study_effects` in the final candidate
registry. The final audit is
`sensitivity_rescor_full_fully_resolved_diagnostic_audit`.

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "rescor")
targets::tar_make(
  sensitivity_rescor_full_observation_loo_audit,
  reporter = "balanced"
)

rescor_observation_loo <- targets::tar_read(
  sensitivity_rescor_full_observation_loo_audit
)
rescor_observation_loo$loo_unit
rescor_observation_loo$outcome_summary
rescor_observation_loo$exact_refits_required
rescor_observation_loo$high_k_cases

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Do not run the existing study-level exact-LOO target. Its 166 requested
leave-study-out fits belong to the retained new-study influence analysis, not
the like-for-like stacking route. The observation audit requires 48 exact
deletions: 45 are confined to five hypertrophy rows from study 124, one is a
strength row, and two are the same whole-body row under two models. Run the
storage-minimal diagnostic gate with:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "rescor")
targets::tar_make(
  sensitivity_rescor_full_observation_exact_diagnostic_audit,
  reporter = "balanced"
)

rescor_observation_exact <- targets::tar_read(
  sensitivity_rescor_full_observation_exact_diagnostic_audit
)
rescor_observation_exact[c(
  "exact_refits",
  "exact_refits_passing",
  "exact_refits_requiring_remediation",
  "all_exact_refit_diagnostics_pass"
)]
rescor_observation_exact$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Each branch immediately diagnoses and scores its refit, archives only the
compact result, and deletes its temporary CmdStan sampler files. If the run is
interrupted, completed compact branches are reused. Do not proceed to final
stacking until this gate passes; any failures receive a targeted second
schedule rather than rerunning all 48 cases.

The initial exact gate retained 40 of 48 cases. The eight failures were all
hypertrophy deletions with exactly one divergence, no treedepth events, and
otherwise passing convergence and efficiency diagnostics. The targeted second
schedule uses a new seed, 3,000 warmup iterations, 8,000 retained iterations,
`adapt_delta = 0.9999`, and `max_treedepth = 16`. It replaces only those eight
compact records:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "rescor")
targets::tar_make(
  sensitivity_rescor_full_observation_exact_complete_diagnostic_audit,
  reporter = "balanced"
)

rescor_observation_exact_complete <- targets::tar_read(
  sensitivity_rescor_full_observation_exact_complete_diagnostic_audit
)
rescor_observation_exact_complete[c(
  "exact_refits",
  "exact_refits_passing",
  "exact_refits_requiring_remediation",
  "all_exact_refit_diagnostics_pass"
)]
rescor_observation_exact_complete$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Six of the eight first-remedy cases passed. Cases 3 and 19 (hypertrophy
observation 252 under models 1 and 9) each retained one divergence, with no
treedepth event and otherwise passing convergence and efficiency diagnostics.
The final sampling-only remedy retains both models unchanged, uses a fresh
seed, 4,000 warmup iterations, 12,000 retained iterations,
`adapt_delta = 0.99999`, and `max_treedepth = 17`, and replaces only those two
compact records:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "rescor")
targets::tar_make(
  sensitivity_rescor_full_observation_exact_fully_resolved_diagnostic_audit,
  reporter = "balanced"
)

rescor_observation_exact_resolved <- targets::tar_read(
  sensitivity_rescor_full_observation_exact_fully_resolved_diagnostic_audit
)
rescor_observation_exact_resolved[c(
  "exact_refits",
  "exact_refits_passing",
  "exact_refits_requiring_remediation",
  "all_exact_refit_diagnostics_pass"
)]
rescor_observation_exact_resolved$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

After the observation-level exact-refit and stacking gates pass, freeze the
compact hand-off and build the reporting bundle with fitting disabled:

```r
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
source("sensitivity/freeze_rescor_full_checkpoint.R")
targets::tar_make(
  sensitivity_rescor_verification_files,
  reporter = "balanced"
)
```

The final bundle must pass before previewing and executing
`sensitivity/prune_targets.R --family=rescor`.

### Focused within-contrast convergence follow-up

The initial whole-body within-study contrast pilot passed all geometry checks.
Models 2, 3, 4, 5, 13, and 26 also passed the strict convergence thresholds.
Only model 1 showed marginal R-hat failures for three nuisance hierarchical
parameters, with adequate ESS and no divergences or treedepth events. Its
focused remedy retains the same model and priors and extends the schedule to
1,500 warmup and 2,000 retained iterations per chain:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "within")
targets::tar_make(sensitivity_within_remedy_fit)
targets::tar_make(sensitivity_within_remedy_diagnostics)
targets::tar_read(sensitivity_within_remedy_diagnostics)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Once all seven remedied candidates pass, compare conditional and integrated
observation-level PSIS. This stage does not refit the models. Keep the family
enabled so the saved fitted branches remain active dependencies:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "within")
targets::tar_make(sensitivity_within_integrated_loo_diagnostics)
targets::tar_read(sensitivity_within_integrated_loo_diagnostics)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Do not calculate within-contrast stacking weights from the conditional LOO
objects. Any integrated Pareto-k values above 0.7 must receive exact
case-deletion correction first.

If all integrated Pareto-k values are at or below 0.7, materialise the final
stacking weights and support-aware model-averaged reporting tables:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "within")
targets::tar_make(sensitivity_within_reporting_outputs)
targets::tar_read(sensitivity_within_integrated_weights)$weights
targets::tar_read(sensitivity_within_support_by_group)
targets::tar_read(sensitivity_within_common_support_curve)
targets::tar_read(sensitivity_within_common_support_increments)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

The absolute curve is the predicted protein-arm minus comparator-arm outcome
contrast while comparator intake is held at 0.8 g/kg/day. The increment table
subtracts the 0.8 reference draw by draw, separating dose-response shape from
any dose-independent arm contrast. Direct resistance/no-resistance comparison
is restricted to common observed treatment-intake support.

### Full within-contrast expansion

The passing seven-model whole-body pilot is frozen in
`results/checkpoints/within_pilot_completed.qs`. The full expansion reuses
those seven fits and fits the remaining 97 candidates: 19 whole-body models
and all 26 models for hypertrophy, strength, and physical function. Every new
fit uses the pilot-informed schedule of 1,500 warmup and 2,000 retained
iterations per chain.

Unlike the original Mundlak expansion, every new fitted object is written
directly to one persistent file. No target collects the 97 new fits or the
104-model ensemble, so the target store does not create redundant copies of
the posterior objects.

First run the no-sampling structure gate:

```r
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")

targets::tar_make(
  sensitivity_within_full_structure_audit,
  reporter = "balanced"
)

within_structure <- targets::tar_read(
  sensitivity_within_full_structure_audit
)
within_structure$checks
within_structure$contrast_summary
within_structure$structure_pass
```

After that gate passes, enable only the within family and request the final
diagnostic audit. This is the long-running sampling command. Branch errors use
`error = "continue"`, so one failed candidate will not prevent later models
from being checkpointed.

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "within")

targets::tar_make(
  sensitivity_within_full_diagnostic_audit,
  reporter = "balanced"
)

within_diagnostic_audit <- targets::tar_read(
  sensitivity_within_full_diagnostic_audit
)

within_diagnostic_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]

within_diagnostic_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Stop at this diagnostic gate. Any failed candidates will receive focused
remedies based on their observed failure mode before integrated LOO,
case-deletion correction, stacking, and reporting are added.

After all full within-study fits and exact case-deletion refits pass, the
final corrected LOO audit must be stacking-ready before reporting. Freeze the
completed analysis state into a lightweight checkpoint:

```r
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
source("sensitivity/freeze_within_checkpoint.R")
```

Then materialise the storage-aware 4,000-curve mixtures, support-aware
reporting summaries, and durable verification bundle:

```r
targets::tar_make(
  sensitivity_within_verification_files,
  reporter = "balanced"
)

within_reporting <- targets::tar_read(
  sensitivity_within_full_reporting_audit
)
within_reporting$summary
within_reporting$all_reporting_checks_pass
```

This stage does not sample or refit models. It opens each selected accepted fit
separately, writes only its allocated prediction curves, and never constructs
a combined 104-fit object. Every one of the 26 candidates remains in each
outcome-specific stacking calculation; finite 4,000-draw allocation naturally
gives prediction draws only to models whose weights are large enough. The
within-study absolute curve is the predicted protein-arm minus comparator-arm
contrast with comparator intake fixed at 0.8 g/kg/day. Its increment curve
subtracts the treatment prediction at 0.8 draw by draw.

Do not prune the accepted or exact-refit files until the reporting audit and
the checksums in `results/verification/within/manifest.csv` have been checked.

### Full achieved-versus-assignment expansion

The passing whole-body pilots are frozen in
`results/checkpoints/exposure_pilot_completed.qs`. The full expansion retains
the final 26-model candidate registry and outcome-specific priors for both
exposure definitions. It reuses 7 achieved and 7 assignment-derived pilot
fits, then writes the remaining 97 fits per exposure directly to independent
files. Achieved and assignment models use identical rows within each outcome.

The no-sampling structure gate is:

```r
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
targets::tar_make(
  sensitivity_exposure_full_structure_audit,
  reporter = "balanced"
)
targets::tar_read(sensitivity_exposure_full_structure_audit)
```

For a resumable expansion, fit and diagnose achieved intake across all four
outcomes first:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(
  sensitivity_achieved_full_diagnostic_audit,
  reporter = "balanced"
)

achieved_audit <- targets::tar_read(
  sensitivity_achieved_full_diagnostic_audit
)
achieved_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]
achieved_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

In the initial achieved-intake expansion, 98/104 candidates passed. The six
focused remedies are hypertrophy models 11, 14, 18, and 23, strength model 6,
and physical-function model 14. Five failures were marginal convergence/MCSE
shortfalls in hierarchical parameters; strength model 6 had one maximum-
treedepth transition. There were no divergences. The remedy keeps all formulas
and priors unchanged and uses 2,000 warmup, 4,000 retained iterations,
`adapt_delta = 0.999`, and maximum treedepth 15:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(
  sensitivity_achieved_full_final_diagnostic_audit,
  reporter = "balanced"
)

achieved_final_audit <- targets::tar_read(
  sensitivity_achieved_full_final_diagnostic_audit
)
achieved_final_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]
achieved_final_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

After reviewing that gate, fit and diagnose the assignment-derived ensemble:

```r
targets::tar_make(
  sensitivity_assignment_full_diagnostic_audit,
  reporter = "balanced"
)

assignment_audit <- targets::tar_read(
  sensitivity_assignment_full_diagnostic_audit
)
assignment_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]
assignment_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

In the initial assignment-derived expansion, 101/104 candidates passed. The
three focused remedies are hypertrophy model 25 and strength models 2 and 22.
The first had one divergence, the second narrowly missed the convergence
threshold, and the third had one maximum-treedepth transition. Formulas and
priors remain unchanged; the remedy uses 2,000 warmup, 4,000 retained
iterations, `adapt_delta = 0.999`, and maximum treedepth 15:

```r
targets::tar_make(
  sensitivity_assignment_full_final_diagnostic_audit,
  reporter = "balanced"
)

assignment_final_audit <- targets::tar_read(
  sensitivity_assignment_full_final_diagnostic_audit
)
assignment_final_audit[c(
  "models",
  "models_passing_sampling",
  "models_requiring_remediation",
  "all_sampling_diagnostics_pass"
)]
assignment_final_audit$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

After all 104 assignment-derived fits pass, materialise the accepted manifest
without refitting any model:

```r
targets::tar_make(
  sensitivity_assignment_full_accepted_fit_manifest,
  reporter = "balanced"
)

Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

### Full achieved/assignment integrated LOO

After both 104-model sampling gates pass, compute the integrated observation-
level PSIS-LOO reductions one exposure definition at a time. Each branch opens
one accepted fit and writes one small reduction directly to disk; the fitted
models are never collected into another large target. Begin with achieved
intake:

```r
targets::tar_make(
  sensitivity_achieved_full_integrated_loo_audit,
  reporter = "balanced"
)

achieved_loo <- targets::tar_read(
  sensitivity_achieved_full_integrated_loo_audit
)
achieved_loo$outcome_summary
achieved_loo$exact_refits_required
```

This stage uses 4,000 posterior draws, integrates out study, arm, and outcome
effects, and only identifies the observation deletions requiring exact
correction. Exact refits and assignment-derived LOO are deliberately gated on
review of this audit.

The achieved-intake audit selected 17 exact deletions: 14 whole-body and 3
hypertrophy cases, across 13 candidate models. Fit and diagnose only those
cases with the conservative exact-refit schedule:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")

targets::tar_make(
  sensitivity_achieved_full_exact_loo_diagnostic_audit,
  reporter = "balanced"
)

achieved_exact <- targets::tar_read(
  sensitivity_achieved_full_exact_loo_diagnostic_audit
)
achieved_exact[c(
  "exact_refits",
  "exact_refits_passing",
  "exact_refits_requiring_remediation",
  "all_exact_refit_diagnostics_pass"
)]
achieved_exact$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass)
```

Each exact refit uses 4 chains, 2,000 warmup and 6,000 retained iterations,
with `adapt_delta` at least 0.999 and maximum treedepth at least 15. No model
without a selected high-k observation is refitted.

Alternatively, requesting `sensitivity_exposure_full_diagnostic_audit` fits
both 97-model expansions in one pipeline invocation. The split commands are
preferred because they provide an earlier diagnostic checkpoint and reduce
the amount of unfinished work present at one time. The pilot-informed initial
schedule is 1,500 warmup plus 2,000 retained iterations per chain. Any failures
are remedied selectively after inspection; they do not cause the full 208-fit
set to be repeated.

### Focused achieved/assignment convergence follow-up

The initial common-sample exposure pilot had clean geometry for all 14 fits.
Longer chains are used only for achieved models 1, 2, 3, and 26 and assignment
models 2 and 26. Model 1 contains no exposure term and is identical across the
two exposure definitions, so its single remedied fit is reused in both
ensembles. Functional forms and priors are unchanged.

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(sensitivity_exposure_remedy_diagnostics)
targets::tar_read(sensitivity_exposure_remedy_diagnostics) |>
  dplyr::select(
    exposure,
    model_index,
    divergent_transitions,
    maximum_treedepth_transitions,
    minimum_ebfmi,
    maximum_rhat,
    minimum_bulk_ess,
    minimum_tail_ess,
    maximum_mcse_sd_ratio,
    geometry_pass,
    convergence_pass,
    diagnostic_pass
  )
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

If the first exposure remedy leaves only achieved-intake model 3 outside the
strict Monte Carlo thresholds, its second-stage remedy retains the same model
and priors and uses 2,000 warmup and 4,000 retained iterations per chain:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(sensitivity_exposure_final_diagnostics)
targets::tar_read(sensitivity_exposure_final_diagnostics) |>
  dplyr::select(
    exposure,
    model_index,
    divergent_transitions,
    maximum_treedepth_transitions,
    minimum_ebfmi,
    maximum_rhat,
    minimum_bulk_ess,
    minimum_tail_ess,
    maximum_mcse_sd_ratio,
    geometry_pass,
    convergence_pass,
    diagnostic_pass
  )
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Once all 14 models pass, run integrated observation-level PSIS for both final
ensembles. This reuses the saved fits and does not sample new posteriors:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(sensitivity_exposure_integrated_loo_diagnostics)
targets::tar_read(sensitivity_exposure_integrated_loo_diagnostics)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

Model 26 is integrated using its study factor-smooth design rather than being
treated as an ordinary random-slope model. Integrated Pareto-k values above
0.7 require exact case-deletion correction before either ensemble is stacked.

The whole-body pilot leaves four such cases: achieved models 3 and 4 at
observation 87, and assignment model 4 at observations 87 and 131. Run only
those exact refits and inspect their diagnostics before applying corrections:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(sensitivity_exposure_exact_loo_refit_diagnostics)
targets::tar_read(sensitivity_exposure_exact_loo_refit_diagnostics)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

After all exact refits pass, apply their cumulative corrections (including
both assignment model-4 observations), calculate final stacking weights, and
materialise the two model-averaged curves:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "exposure")
targets::tar_make(sensitivity_exposure_final_outputs)
targets::tar_read(sensitivity_achieved_integrated_weights)$weights
targets::tar_read(sensitivity_assignment_integrated_weights)$weights
targets::tar_read(sensitivity_achieved_final_loo_audit)$model_summary
targets::tar_read(sensitivity_assignment_final_loo_audit)$model_summary
targets::tar_read(sensitivity_exposure_final_comparison)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

### Focused Mundlak convergence follow-up

After reviewing the first whole-body pilot, the follow-up stage refits models
1, 2, 3, 4, 13, and 26 with 1,500 warmup and 2,000 retained iterations per
chain. Model 4 also receives the same `student_t(3, 0, 1)` regularisation of
hierarchical standard deviations already used by remedied primary candidates.
Its thin-plate basis, covariates, and Mundlak decomposition are unchanged.
The passing model-5 fit is reused, and the follow-up diagnostics and averaging
targets reconstruct the complete seven-model ensemble.

For LOO comparison, the extended refits are deterministically thinned to the
same 4,000 total posterior draws available for model 5. This is required by
the `loo` model-weight interface and affects only the LOO calculation; curve
generation can still sample from every saved posterior draw.

The first LOO target is intentionally conditional and diagnostic. Multilevel
effects can make its observation-level Pareto values unstable. The integrated
follow-up analytically marginalises the Gaussian study, arm, and outcome
effects—including correlated study-level Mundlak basis slopes—using the same
strategy as the final primary-model LOO workflow:

```r
targets::tar_make(sensitivity_mundlak_integrated_loo_diagnostics)
targets::tar_read(sensitivity_mundlak_integrated_loo_diagnostics)
```

`sensitivity_mundlak_integrated_weights` refuses to calculate stacking weights
if the integrated pass retains any Pareto-k value above 0.7. Such cases need
exact case-deletion correction before model averaging.

In the initial whole-body pilot, integration left one warning for each of
models 3, 4, and 26, all for observation 3 (study 8, arm 7, outcome 43).
The exact stage therefore contains only three cached refits:

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "mundlak")
targets::tar_make(sensitivity_mundlak_exact_loo_refit)
targets::tar_make(sensitivity_mundlak_exact_loo_refit_diagnostics)
targets::tar_read(sensitivity_mundlak_exact_loo_refit_diagnostics)
```

After those diagnostics pass, making
`sensitivity_mundlak_integrated_weights` applies the exact corrections first
and then calculates final stacking weights. The downstream averaged curve is
still blocked if an exact refit fails its sampler diagnostics.

The reporting pilot retains the full prediction grid but marks every point as
common central support, common observed-tail support, group-only central
support, group-only observed-tail support, or extrapolation. It also reports
posterior changes from 0.8 g/kg/day, which separate dose-response magnitude
from the fitted outcome intercept. Direct resistance/no-resistance comparison
uses only their common observed range.

```r
targets::tar_make(sensitivity_mundlak_reporting_outputs)
targets::tar_read(sensitivity_mundlak_support_by_group)
targets::tar_read(sensitivity_mundlak_common_support_increments)
targets::tar_read(sensitivity_mundlak_group_supported_increments)
```

`sensitivity_mundlak_reporting_pilot` is the composite object used as an
upstream dependency. Making it alone does not materialise the downstream
convenience targets listed above.

```r
Sys.setenv(PROTEIN_SENSITIVITY_PILOTS = "mundlak")
targets::tar_make(sensitivity_mundlak_remedy_fit)
targets::tar_make(sensitivity_mundlak_remedy_diagnostics)
targets::tar_read(sensitivity_mundlak_remedy_diagnostics)
Sys.unsetenv("PROTEIN_SENSITIVITY_PILOTS")
```

The equivalent convenience commands are:

```sh
Rscript sensitivity/run_pilot.R mundlak-remedy fit
Rscript sensitivity/run_pilot.R mundlak-remedy diagnostics
```

The convenience runner provides the same guard from a fresh R process:

```sh
Rscript sensitivity/run_pilot.R structure
Rscript sensitivity/run_pilot.R mundlak fit
Rscript sensitivity/run_pilot.R mundlak summary
Rscript sensitivity/run_pilot.R eiv diagnostics
```

Never set the flag to `all` before an unqualified `tar_make()` unless all five
pilot families are genuinely intended to run.

## What each pilot implements

### Mundlak

Every column of the protein basis is decomposed into its study mean and its
arm-specific deviation from that mean. Study means are calculated from unique
arms so arms with more reported outcomes do not receive extra weight. Linear
and segmented-linear spans are exact. Thin-plate and natural-spline candidates
use their corresponding basis spans. The thin-plate basis is rotated into its
unpenalized and penalty-scaled components so the targeted smooth prior can be
applied to the penalized columns; study-specific smooths are represented as
random coefficients on the within-study basis. Predictions vary the
within-study basis while holding the between-study basis at the habitual
protein reference.

### Strict within-study contrasts

Protein arms are matched to comparator arms within study, outcome,
measurement, time contrast, and resistance-training condition. The response is
the difference in standardised change and each protein basis is the treatment
minus comparator basis. The initial variance is the sum of the two arm-level
variances, which assumes independent arms. Multi-arm trials that reuse a
comparator therefore require special attention in the pilot diagnostics before
this becomes a final sensitivity analysis. The primary intercept and RT prior
scales are retained but their locations are recentered at zero because the
response is now an arm contrast rather than an arm-level change.

### Achieved versus assignment-derived intake

Assignment-derived total intake is baseline dietary protein plus prescribed
protein grams per day divided by baseline body mass. A missing prescribed dose
is set to zero only for control/placebo arms; it remains missing for protein
arms. Achieved and assignment-derived candidate sets use the same complete-case
rows, but each exposure is centred on its own comparator-arm median. The final
primary formulas, outcome-specific priors, remedy decisions, and pilot stacking
and posterior-mixture summaries are retained.

### Residual correlation

The representative default is candidate 6: a linear protein model with age and
duration moderators. The hand-written Stan model estimates a common
within-study residual correlation in the known sampling-error component while
retaining study, arm, outcome, and residual heterogeneity terms. The correlation
is constrained to the valid compound-symmetry range for the largest study.

### Errors in variables

The generic hand-written Stan program supports the protein-containing candidate
architecture: linear, segmented-linear, and nonlinear cubic bases; age and
duration moderators; RT interactions; study random basis slopes; and arm and
outcome random intercepts. It models one latent true intake per study-arm from
the reported intake and its SE, adds residual heterogeneity to the known outcome
SE, and emits outcome log-likelihood and population prediction draws for pilot
stacking and posterior-mixture curves.

Linear and segmented-linear conversions preserve their original functional
forms. For latent nonlinear intake, the current Stan pilot uses an orthogonal
three-degree-of-freedom restricted/natural cubic representation in place of the
`mgcv` thin-plate parameterisation, because an `mgcv` design matrix cannot be
precomputed at an unknown latent exposure. The strength model-4 natural-spline
remedy uses its boundary and one-third/two-third quantile knot convention; the
other nonlinear models use four interior-range quantiles. This preserves
nonlinear flexibility and the natural-spline function space, while the
thin-plate conversion is intentionally marked as an approximation to audit
during piloting, not silently treated as an exact conversion.

## Pilot outputs

Each family has separate fit, diagnostic, integrated/exact PSIS-LOO,
stacking-weight, posterior-mixture, and curve-summary targets. The original
residual-correlation pilot remains available as a one-model smoke test; the
full residual-correlation expansion now carries the complete 104-candidate
registry through study-level LOO and the storage-aware reporting workflow.

## Cross-sensitivity supplementary comparisons

After a full sensitivity family has written a verified
`results/verification/<family>/` bundle, rebuild the supplementary comparison
figures and agreement tables with:

```r
source("sensitivity/plot_full_sensitivity_comparisons.R")
```

The script discovers completed bundles that contain both
`dose_response_predictions.csv.gz` and `pairwise_contrasts.csv.gz`. It reads
cached primary summaries and verification files only and never fits a model.
The protein-intake strategy families will therefore be added automatically
once they use the same verified reporting schema.

The overview figures use the primary-analysis `theme_classic()` styling,
Okabe-Ito colours, posterior means, and marginal 95% posterior quantile
intervals. Absolute curves omit the strict within-study family because its
level is an arm contrast rather than an absolute outcome prediction. The
0.8-g/kg/day reference-increment figures include every family; proposed doses
are restricted to their observed support, and captions identify that the 0.8
reference itself may be extrapolated for some outcome/RT panels.

Full heatmap comparisons are limited to dose pairs supported by both the
primary and sensitivity analyses and exclude diagonal/duplicated decreases
from numeric summaries. Reported agreement measures comprise differences in
posterior mean contrasts, overlap coefficients for marginal 95% intervals,
four-way decision-classification agreement (benefit, equivalence, detriment,
or uncertain), and total-variation distance between the three posterior
classification-probability vectors. Marginal interval overlap is descriptive;
it is not a posterior probability that separately fitted analyses agree.
