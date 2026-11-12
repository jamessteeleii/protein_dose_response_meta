# Updated-search extraction workbook

`protein_meta_data_update_2026.xlsx` is the isolated extraction workbook for
studies found by the 2026 search update. Do not enter updated-search data in
`data/protein_meta_data_all.xlsx`; that workbook remains the frozen canonical
input for the completed primary and sensitivity analyses.

## Provenance

The workbook was created on 21 September 2026 from
`data/protein_meta_data_all.xlsx` (SHA-256
`79b9ca399b09af20a31e50bdda0ea3e49b33b2682bfbfcdc30148185976371b4`).
It retains the four original worksheet names, column headings, ordering,
formatting, and workbook annotations. It was subsequently populated with the
updated-search records and passed its final extraction audit on 28 September
2026.

## Final extraction status

- 10 study records, 29 arm records, 324 arm-by-timepoint rows, and 149 reported
  contrast rows are present.
- All 149 reported contrast keys match their corresponding timepoint keys.
  Thirteen additional contrasts will be generated from complete baseline and
  follow-up rows during preparation.
- No duplicate or orphan keys, invalid sample sizes, negative SDs,
  out-of-range correlations, or spreadsheet formula errors were found.
- Ceglia et al. is retained at study and intended-arm level but has no
  extractable randomized-arm outcome or contrast data.
- Seyhan et al. is retained for provenance but is designated for exclusion
  from quantitative synthesis by author-group decision.
- Huijing et al. uses the reported longitudinal change estimates and the
  week-24 longitudinal protein-intake table. Differences between reported
  change and arithmetic differences of rounded marginal means are intentional.
- The resistance and concurrent-training arms from Bagheri et al. and
  Vaithanomsat et al. are all retained for consistency with the existing
  eligibility approach.

## Completed-suite computational warning

The updated-search analyses are complete. This file remains the detailed,
checkpoint-by-checkpoint historical runbook; it is not a script to execute
from top to bottom. The main graph alone contributed 2,078 distinct final
candidate models and 3,073 completed model-fitting/refitting jobs. The full
three-graph project total is 2,662 candidates, 4,266 fitting jobs, at least
1,659 model-job hours, and about 747 GiB of retained local computational state
at the 12 November 2026 audit.

Before rerunning any stage, read the root
[`REPRODUCIBILITY.md`](../../REPRODUCIBILITY.md). Enable only the switch named
for the reviewed checkpoint, preserve successful branches and fit archives,
and stop for the documented diagnostic review. Do not run all commands in
this file as a batch, do not point two processes at the same target store, and
do not delete `_targets/` or `model_archive/` to obtain a clean restart.

## Extraction rules

- Do not rename worksheets or columns, and do not change column order.
- Confirm that a report is not a duplicate or companion report of an existing
  study before assigning a new `Study No`.
- The largest study number in the frozen canonical workbook is 162. If all ten
  update records are genuinely new studies, assign 163 onward; companion
  reports should share the appropriate study number instead.
- Use the same `Study No` consistently on all four worksheets. Number arms
  within each study and use the same `Arm No` everywhere that arm appears.
- Preserve the existing outcome-group labels and record the reported unit of
  measurement exactly. Do not convert or standardise values in this workbook.
- Leave unavailable quantities blank rather than entering zero. Put any
  extraction decisions or ambiguities in cell comments so they can be reviewed
  before integration.
- Add rows as needed while copying the formatting from the immediately
  preceding template row.

This file is the sole new data input to the updated-analysis pipeline.
Integration with the canonical dataset occurs only through the explicit,
audited in-memory combination described below.

## Combined-data refit interface

The update branch is prefixed `update_2026_` in `_targets.R`. It appends the
canonical and updated-search workbooks in memory, then performs one preparation
pass across all studies. Neither source workbook is modified. This common
prepared object supplies the updated primary analysis and will subsequently
supply each sensitivity family.

Study-level quantitative decisions are recorded in `review_decisions.csv`.
Ceglia et al. remains review-only and Seyhan et al. remains excluded from
quantitative synthesis. Preparation checks and comparisons with the frozen
analysis are written to `results/update_2026/primary_refit/`.

The primary refit starts from the final accepted specification and sampler
settings for each of the 104 frozen candidate models. Each fitted object is
written immediately beneath
`model_archive/update_2026/primary/candidates/`, so a
completed model survives an interrupted run without duplicating the full fit
inside the targets object store. All 104 fits receive a fresh diagnostic audit
before any LOO calculation or model averaging proceeds.

## Local primary-refit runbook

Large fitting stages are intentionally never launched during pipeline
development. The environment switch below is required before the 104-model
fit target can run; an ordinary `tar_make()` leaves those fits frozen.

Run the following commands from an R session whose working directory is the
project root. The examples use two simultaneous models while each model runs
four parallel chains. Reduce `workers` if memory pressure is high.

1. Run the data audit and the lightweight all-model preflight. This constructs
   Stan data for every final specification but does not compile or sample a
   model.

   ```r
   targets::tar_make(
     update_2026_combined_primary_preparation_audit_files,
     reporter = "balanced"
   )
   targets::tar_make(
     update_2026_primary_input_preflight_file,
     reporter = "balanced"
   )
   ```

   Confirm that `results/update_2026/primary_refit/input_preflight.csv`
   contains 104 rows and that every `valid` value is `TRUE` before moving on.

2. Explicitly enable and run the expensive primary fits.

   Before starting, check free space on the project drive. The frozen primary
   archive is the best empirical size guide for this candidate set;
   retain comfortable additional headroom for temporary compilation files,
   any diagnostic remedies, and later exact-LOO refits. Do not start if the
   projected archive would leave less than roughly one additional archive's
   worth of free space.

   ```r
   Sys.setenv(PROTEIN_UPDATE_2026_PRIMARY = "fit")
   targets::tar_make_future(
     update_2026_primary_fit_file,
     workers = 2,
     reporter = "balanced"
   )
   ```

   The command is resumable: rerun the same command after an interruption and
   `targets` will retain completed branches. A complete first pass contains 104
   non-empty `model_*.qs` files beneath
   `model_archive/update_2026/primary/candidates/`. These archive files are the
   single authoritative posterior copies; the target store contains paths,
   not duplicate fitted objects.

3. Once all fits are complete, create the fresh diagnostic audit.

   ```r
   targets::tar_make_future(
     update_2026_primary_diagnostic_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   The durable audit is written beneath
   `results/update_2026/primary_refit/diagnostics/`. Stop at this review
   gate: any newly failing model should receive the smallest appropriate
   remedy and be refitted before updated LOO or stacking is prepared.

4. The first-pass audit found one sampling failure: function model 10 had one
   divergent transition at `adapt_delta = 0.999`, with all convergence and
   efficiency criteria passing. Its first remedy retains the model and draw
   schedule and changes only `adapt_delta` to `0.9995`. Run this isolated
   branch with the initial-fit switch disabled:

   ```r
   Sys.unsetenv("PROTEIN_UPDATE_2026_PRIMARY")
   Sys.setenv(PROTEIN_UPDATE_2026_PRIMARY_REMEDY = "fit")
   targets::tar_make_future(
     update_2026_primary_remedy_diagnostic_audit_files,
     workers = 1,
     reporter = "balanced"
   )
   ```

   Inspect `results/update_2026/primary_refit/remedy_1/diagnostic_summary.csv`
   before proceeding. If `diagnostic_pass` is `TRUE`, create the authoritative
   104-fit path registry:

   ```r
   targets::tar_make(
     update_2026_primary_accepted_fit_manifest_file,
     reporter = "balanced"
   )
   Sys.unsetenv("PROTEIN_UPDATE_2026_PRIMARY_REMEDY")
   ```

   This does not duplicate posterior objects. The accepted manifest points to
   103 initial files and the one accepted remedy file.

5. After the accepted manifest is complete, run the integrated PSIS-LOO stage.
   This computes the conventional conditional PSIS result and the final
   study-effect-integrated observation-level PSIS result one accepted model at
   a time. Only the compact integrated reductions are retained.

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_PRIMARY",
     "PROTEIN_UPDATE_2026_PRIMARY_REMEDY"
   ))
   Sys.setenv(PROTEIN_UPDATE_2026_PRIMARY_LOO = "run")
   targets::tar_make_future(
     update_2026_primary_integrated_loo_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   The run is resumable. Each completed reduction is written beneath
   `model_archive/update_2026/primary/loo/integrated/`; the target store keeps
   paths rather than copies. The audit is written beneath
   `results/update_2026/primary_refit/loo_psis/`.

   Review `outcome_summary.csv`, `model_summary.csv`, and `high_k_cases.csv`,
   then stop. Exact case-deletion refits are generated only for the newly
   observed Pareto-k warnings and therefore require a separate, data-driven
   target prepared after this audit.

6. The updated integrated audit identified three exact deletions: whole-body
   models 3 and 19 at observation 3, and hypertrophy model 13 at observation
   253. First verify the held-out rows and random-slope design without fitting:

   ```r
   targets::tar_make(
     update_2026_primary_exact_loo_input_preflight_file,
     reporter = "balanced"
   )
   ```

   Confirm the three rows in
   `results/update_2026/primary_refit/loo_exact/input_preflight.csv` are valid.
   Then run the storage-minimal exact stage with all earlier switches disabled:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_PRIMARY",
     "PROTEIN_UPDATE_2026_PRIMARY_REMEDY",
     "PROTEIN_UPDATE_2026_PRIMARY_LOO"
   ))
   Sys.setenv(PROTEIN_UPDATE_2026_PRIMARY_EXACT = "fit")
   targets::tar_make_future(
     update_2026_primary_exact_loo_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   Each branch retains only its compact corrected LOO object beneath
   `model_archive/update_2026/primary/loo/exact/`; the temporary case-deletion
   `brmsfit` is released rather than archived. Review the diagnostic audit at
   `results/update_2026/primary_refit/loo_exact/` and stop again before final
   correction assembly or stacking.

7. Once all exact-refit diagnostics pass, assemble the final 104-path LOO
   manifest and calculate stacking weights. This stage does not fit models and
   does not duplicate the compact LOO objects.

   ```r
   Sys.unsetenv("PROTEIN_UPDATE_2026_PRIMARY_EXACT")
   targets::tar_make(
     update_2026_primary_stacking_checkpoint_files,
     reporter = "balanced"
   )
   ```

   Final LOO audits are written beneath
   `results/update_2026/primary_refit/loo_final/`; stacking weights are written
   to `results/update_2026/primary_refit/stacking_weights.csv`. Review this
   checkpoint before model-averaged predictions and reporting are prepared.

8. Once the final LOO and stacking checkpoint passes, generate the updated
   model-averaged reporting data. This stage does not refit models. Each of the
   four outcome branches opens its 26 accepted fits one at a time, draws only
   its allocated posterior predictions, and retains one compact 4,000-draw
   mixture. Use no more than two workers to keep peak memory bounded:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_PRIMARY",
     "PROTEIN_UPDATE_2026_PRIMARY_REMEDY",
     "PROTEIN_UPDATE_2026_PRIMARY_LOO",
     "PROTEIN_UPDATE_2026_PRIMARY_EXACT"
   ))
   Sys.setenv(
     PROTEIN_UPDATE_2026_PRIMARY_REPORTING = "run",
     PROTEIN_MODEL_WORKERS = "2"
   )
   targets::tar_make_future(
     update_2026_primary_reporting_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   The run is resumable by outcome. Compact posterior mixtures and reporting
   objects are written beneath
   `model_archive/update_2026/primary/reporting/`; auditable flat tables are
   written beneath `results/update_2026/primary_refit/reporting/`. Review
   `summary.csv`, `model_weights.csv`, and `draw_allocation.csv`, then stop
   before preparing figures, the final verification bundle, or pruning.

9. Once every reporting check passes, create the final updated-primary
   verification bundle. This is a lightweight assembly step: it substitutes
   the accepted function-model remedy into the diagnostic and specification
   registries, checks the complete 104-model lineage, hashes the three source
   data files, and writes compact reporting tables. It does not open or refit
   any model.

   ```r
   Sys.unsetenv("PROTEIN_UPDATE_2026_PRIMARY_REPORTING")
   targets::tar_make(
     update_2026_primary_verification_files,
     reporter = "balanced"
   )
   ```

   The durable bundle is written beneath
   `results/verification/primary_update_2026/`. Review the checkpoint target
   and bundle manifest before preparing figures or a separate pruning plan.

10. Render the four updated primary figures only after the verification
    checkpoint passes. These plots read the verified reporting tables and
    combined analysis rows; they do not load model fits or posterior mixtures.

    ```r
    targets::tar_make(
      update_2026_primary_figure_files,
      reporter = "balanced"
    )
    ```

    TIFF and manuscript PDF versions, together with a checksum manifest, are
    written beneath `plots/update_2026/primary/`. Inspect these before any
    updated-primary model archive is pruned.

11. Generate a read-only pruning plan after the verification bundle and all
    eight figures have passed inspection:

    ```r
    targets::tar_make(
      update_2026_primary_pruning_plan_files,
      reporter = "balanced"
    )
    ```

    The plan compares two explicit retention policies. `conservative` removes
    only the superseded initial function-model fit and keeps all 104 accepted
    fits locally reusable. `lean` mirrors the completed sensitivity-family
    policy: it removes every full fit after reporting is frozen, while keeping
    compact LOO/reporting objects, all audits, the verification bundle, and
    figures. This target is dry-run only and never deletes a file.

12. The author-selected policy for the primary analysis is `conservative`.
    After reviewing the plan, execute it explicitly from R:

    ```r
    pruning_plan <- targets::tar_read(
      update_2026_primary_pruning_plan
    )
    source("R/incremental_update.R")
    execute_update_2026_primary_pruning(
      plan = pruning_plan,
      mode = "conservative",
      confirmation = "DELETE SUPERSEDED UPDATED PRIMARY FIT"
    )
    ```

    The helper revalidates the verification and figure checksums, deletes only
    the superseded initial function model 10, confirms that all 104 accepted
    fits remain, and writes an execution receipt beneath
    `results/update_2026/primary_refit/pruning/`. Updated sensitivity families
    will instead use their `lean` policy only after their own verification
    bundles pass.

Closing the R session clears its opt-in environment switches. They can also be
cleared explicitly with `Sys.unsetenv()` as shown above.

The same gated, file-backed, review-between-stages pattern will be used for
each sensitivity family and later for the separate absolute-total-intake and
supplement-dose model tracks.

## Shared absolute-exposure audit

The absolute-total-intake and prescribed-supplement programmes begin with one
sampling-free audit. It uses reported total protein g/day where available and
otherwise derives total intake as g/kg/day multiplied by mean body mass. The
audit retains the source of every selected value and propagates uncertainty
for derived values only when both intake and body-mass uncertainty are
available, using a first-order delta method that treats the two arm means as
independent. Missing supplement doses are set to zero only for control or
placebo arms; reported nonzero comparator protein is retained, and missing
protein-arm doses remain missing.

The audit also freezes the candidate-model decision. It retains the original
model indices but omits indices 3, 7, 10, 13, 16, 19, 22, and 25 because their
fixed 1.6 g/kg/day linear-spline knot has no unique g/day equivalent. This
leaves 18 candidates per outcome and 72 initial fits per applicable family.

Build the review files without enabling any fitting switch:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_PRIMARY",
  "PROTEIN_UPDATE_2026_ACHIEVED",
  "PROTEIN_UPDATE_2026_ASSIGNMENT",
  "PROTEIN_UPDATE_2026_EIV",
  "PROTEIN_UPDATE_2026_RESCOR",
  "PROTEIN_UPDATE_2026_STRATEGY",
  "PROTEIN_UPDATE_2026_STRATEGY_WITHIN"
))

targets::tar_make(
  update_2026_absolute_exposure_audit_files,
  reporter = "balanced"
)
```

Review `results/update_2026/absolute_exposure/audit/` before constructing the
first total-g/day primary model. In particular, confirm `checks.csv`,
`outcome_coverage.csv`, `reported_derived_agreement.csv`, `model_registry.csv`,
`exposure_provenance.csv`, `support_by_outcome.csv`, `scaling_plan.csv`, and
`family_plan.csv`.

After that review is confirmed, build the sampling-free total-intake primary
preflight:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_PRIMARY",
  "PROTEIN_UPDATE_2026_PRIMARY_REMEDY",
  "PROTEIN_UPDATE_2026_PRIMARY_LOO",
  "PROTEIN_UPDATE_2026_PRIMARY_EXACT",
  "PROTEIN_UPDATE_2026_PRIMARY_REPORTING",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_input_preflight_file,
  reporter = "balanced"
)
```

This constructs Stan data but does not sample. It must contain 72 valid model
inputs: 18 candidates for each outcome. Total intake is centred at the audited
78.8 g/day control/placebo median and represented per 50 g/day. The nine
moderator candidates per outcome also include outcome-standardised body mass;
the other nine retain every otherwise eligible contrast. The preflight records
the resulting rows, studies, arms, exposure range, sampler controls, and draw
schedule in
`results/update_2026/absolute_exposure/total_g_day/primary/input_preflight.csv`.
The registry starts from the final accepted g/kg/day sampling lineage rather
than the obsolete original controls. All 32 retained thin-plate candidates use
the accepted smooth prior, `adapt_delta = 0.999`, and maximum tree depth 15;
their exposure is already scaled per 50 g/day, so the earlier half-unit
g/kg/day transformation is not carried over. The strength model-4 thin-plate
form is restored because the exposure and scale have changed.
The primary fit targets remain deliberately disabled until this checkpoint is
reviewed. Once it passes, enable only this family and run four resumable model
branches at a time:

```r
Sys.setenv(
  PROTEIN_MODEL_WORKERS = "4",
  PROTEIN_UPDATE_2026_TOTAL_PRIMARY = "fit"
)

targets::tar_make_future(
  update_2026_total_primary_diagnostic_audit_files,
  workers = 4,
  reporter = "balanced"
)
```

Each branch fits four chains and writes its accepted candidate immediately to
`model_archive/update_2026/absolute_exposure/total_g_day/primary/initial/`.
Completed branches are retained after an interruption. The stage stops at the
72-model sampling-diagnostic review; it does not start LOO, remediation, or
reporting. After completion, clear the fitting switches and inspect the audit:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY",
  "PROTEIN_MODEL_WORKERS"
))

total_primary_files <- unlist(
  targets::tar_read(update_2026_total_primary_fit_file),
  use.names = FALSE
)
total_primary_manifest <- targets::tar_read(
  update_2026_total_primary_fit_manifest
)
total_primary_diagnostics <- targets::tar_read(
  update_2026_total_primary_diagnostic_audit
)

c(
  files = length(total_primary_files),
  manifest_rows = nrow(total_primary_manifest),
  all_files_exist = all(file.exists(total_primary_files)),
  all_files_nonempty = all(file.info(total_primary_files)$size > 0),
  all_hashes_present = all(nchar(total_primary_manifest$sha256) == 64L),
  models = total_primary_diagnostics$models,
  models_passing = total_primary_diagnostics$models_passing_sampling,
  models_requiring_remediation =
    total_primary_diagnostics$models_requiring_remediation
)

total_primary_diagnostics$outcome_summary |>
  print(n = Inf, width = Inf)

total_primary_diagnostics$diagnostic_summary |>
  dplyr::filter(!diagnostic_pass) |>
  dplyr::select(
    outcome,
    model_index,
    model,
    divergent_transitions,
    maximum_treedepth_transitions,
    minimum_ebfmi,
    maximum_rhat,
    minimum_bulk_ess,
    minimum_tail_ess,
    maximum_mcse_sd_ratio,
    diagnostic_pass,
    status
  ) |>
  print(n = Inf, width = Inf)

stopifnot(
  length(total_primary_files) == 72L,
  nrow(total_primary_manifest) == 72L,
  all(file.exists(total_primary_files)),
  all(file.info(total_primary_files)$size > 0),
  all(total_primary_manifest$bytes > 0),
  all(nchar(total_primary_manifest$sha256) == 64L),
  total_primary_diagnostics$models == 72L
)
```

The initial run identified one isolated sampling issue: hypertrophy model 20
had one divergent transition, with no tree-depth, E-BFMI, R-hat, ESS, or MCSE
failure. Build and inspect its sampling-free remedy preflight:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_remedy_preflight_file,
  reporter = "balanced"
)

total_primary_remedy_preflight <- targets::tar_read(
  update_2026_total_primary_remedy_preflight
)
total_primary_remedy_preflight |>
  print(n = Inf, width = Inf)

stopifnot(
  nrow(total_primary_remedy_preflight) == 1L,
  total_primary_remedy_preflight$outcome == "hypertrophy",
  total_primary_remedy_preflight$model_index == 20L,
  total_primary_remedy_preflight$adapt_delta == 0.9995,
  total_primary_remedy_preflight$max_treedepth == 15L,
  total_primary_remedy_preflight$iter_warmup == 2000L,
  total_primary_remedy_preflight$iter_sampling == 6000L,
  total_primary_remedy_preflight$chains == 4L,
  total_primary_remedy_preflight$parallel_chains == 4L,
  total_primary_remedy_preflight$initial_fit_exists,
  total_primary_remedy_preflight$valid
)
```

The remedy retains the formula, data, priors, draw schedule, and tree-depth
limit and changes only `adapt_delta` from 0.999 to 0.9995. Sampling remains
disabled until this one-row checkpoint is reviewed.

After review, fit only the one remedy branch:

```r
Sys.setenv(
  PROTEIN_MODEL_WORKERS = "1",
  PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY = "fit"
)

targets::tar_make(
  update_2026_total_primary_final_diagnostic_audit_files,
  reporter = "balanced"
)
```

Once it finishes, clear all total-intake switches and build the sampling-free
accepted-fit and LOO checkpoint:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_LOO",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_loo_preflight_files,
  reporter = "balanced"
)

total_primary_accepted <- targets::tar_read(
  update_2026_total_primary_accepted_fit_manifest
)
total_primary_loo_preflight <- targets::tar_read(
  update_2026_total_primary_loo_preflight
)

total_primary_loo_preflight$checks |>
  print(n = Inf, width = Inf)

total_primary_accepted |>
  dplyr::count(outcome, stage) |>
  print(n = Inf)

c(
  models = nrow(total_primary_accepted),
  initial_fits = sum(total_primary_accepted$stage == "initial"),
  remedy_fits = sum(total_primary_accepted$stage == "remedy_1"),
  all_files_exist = all(file.exists(total_primary_accepted$fit_path)),
  all_files_nonempty =
    all(file.info(total_primary_accepted$fit_path)$size > 0),
  all_hashes_present =
    all(nchar(total_primary_accepted$sha256) == 64L),
  ready = total_primary_loo_preflight$ready
)

stopifnot(
  nrow(total_primary_accepted) == 72L,
  sum(total_primary_accepted$stage == "initial") == 71L,
  sum(total_primary_accepted$stage == "remedy_1") == 1L,
  all(file.exists(total_primary_accepted$fit_path)),
  all(file.info(total_primary_accepted$fit_path)$size > 0),
  all(nchar(total_primary_accepted$sha256) == 64L),
  all(total_primary_loo_preflight$checks$passed),
  isTRUE(total_primary_loo_preflight$ready)
)
```

This checkpoint must preserve 18 accepted models per outcome, with 71 initial
fits and the accepted hypertrophy-model-20 remedy. It does not compute LOO.
After it is reviewed, enable only integrated observation-level PSIS-LOO. Each
branch uses one core internally; use two concurrent branches initially because
the accepted 6,000-draw fit objects are large:

```r
Sys.setenv(
  PROTEIN_MODEL_WORKERS = "2",
  PROTEIN_UPDATE_2026_TOTAL_PRIMARY_LOO = "run"
)

targets::tar_make_future(
  update_2026_total_primary_integrated_loo_audit_files,
  workers = 2,
  reporter = "balanced"
)
```

Completed LOO branches are retained after interruption. This stage writes
compact reductions under
`model_archive/update_2026/absolute_exposure/total_g_day/primary/loo/` and
stops at the integrated PSIS audit. Any observations above Pareto-k 0.7 are
handled by a later exact-refit checkpoint; this command cannot start those
refits.

The completed LOO audit had no Pareto-k values above 0.7. It also exposed one
cross-model comparability issue that must be resolved before stacking:
hypertrophy moderator candidates used the 421-row body-mass-complete sample,
whereas the nine non-moderator candidates used 429 rows. Stacking requires
the same pointwise observations for every candidate within an outcome. The
other three outcomes already have identical observation sets. Build the
sampling-free preflight for the nine necessary hypertrophy common-sample
refits:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_LOO",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_common_sample_preflight_file,
  reporter = "balanced"
)

total_primary_common_preflight <- targets::tar_read(
  update_2026_total_primary_common_sample_preflight
)

total_primary_common_preflight |>
  dplyr::select(
    outcome,
    model_index,
    model_name,
    contrasts,
    studies,
    common_observation_set,
    adapt_delta,
    max_treedepth,
    iter_warmup,
    iter_sampling,
    chains,
    parallel_chains,
    superseded_fit_exists,
    valid
  ) |>
  print(n = Inf, width = Inf)

stopifnot(
  nrow(total_primary_common_preflight) == 9L,
  identical(
    total_primary_common_preflight$model_index,
    c(1L, 2L, 4L, 9L, 11L, 15L, 17L, 21L, 23L)
  ),
  all(total_primary_common_preflight$contrasts == 421L),
  all(total_primary_common_preflight$studies == 53L),
  all(total_primary_common_preflight$common_observation_set),
  all(total_primary_common_preflight$superseded_fit_exists),
  all(total_primary_common_preflight$valid),
  all(total_primary_common_preflight$chains == 4L),
  all(total_primary_common_preflight$parallel_chains == 4L)
)
```

After review, enable only this stage. These are replacement fits for stacking,
not sampling remedies; they preserve each model's existing priors and sampler
controls while changing only the analysis rows:

```r
Sys.setenv(
  PROTEIN_MODEL_WORKERS = "3",
  PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE = "fit"
)

targets::tar_make_future(
  update_2026_total_primary_common_sample_diagnostic_audit_files,
  workers = 3,
  reporter = "balanced"
)
```

The nine original 429-row hypertrophy fits and their LOO reductions remain
preserved as superseded audit artifacts. No existing accepted fit or compact
LOO result is overwritten.

The common-sample diagnostic audit retained eight fits and identified one
convergence-only remedy: hypertrophy model 11 had clean geometry but two
random-smooth parameters below the tail-ESS threshold. Build its sampling-free
preflight:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_REMEDY",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_common_sample_remedy_preflight_file,
  reporter = "balanced"
)
```

The remedy retains the formula, 421-row common sample, priors,
`adapt_delta = 0.999`, tree depth 15, and 2,000 warmup iterations. It doubles
sampling from 6,000 to 12,000 iterations per chain. After review, enable only
the one remedy fit:

```r
Sys.setenv(
  PROTEIN_MODEL_WORKERS = "1",
  PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_REMEDY = "fit"
)

targets::tar_make(
  update_2026_total_primary_common_sample_final_diagnostic_audit_files,
  reporter = "balanced"
)
```

The remedy audit must accept all nine common-sample models. Then clear every
total-intake fitting switch and build the sampling-free replacement-LOO
checkpoint:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_LOO",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_LOO",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_common_sample_loo_preflight_files,
  reporter = "balanced"
)

total_common_final <- targets::tar_read(
  update_2026_total_primary_common_sample_final_diagnostic_audit
)
total_common_accepted <- targets::tar_read(
  update_2026_total_primary_common_sample_accepted_fit_manifest
)
total_common_loo_preflight <- targets::tar_read(
  update_2026_total_primary_common_sample_loo_preflight
)

c(
  models = total_common_final$models,
  models_passing = total_common_final$models_passing_sampling,
  models_requiring_remediation =
    total_common_final$models_requiring_remediation,
  accepted_models = nrow(total_common_accepted),
  preflight_ready = total_common_loo_preflight$ready
)

total_common_accepted |>
  dplyr::count(stage) |>
  print(n = Inf)

total_common_loo_preflight$checks |>
  print(n = Inf, width = Inf)

stopifnot(
  total_common_final$models == 9L,
  total_common_final$models_passing_sampling == 9L,
  total_common_final$models_requiring_remediation == 0L,
  total_common_final$all_sampling_diagnostics_pass,
  nrow(total_common_accepted) == 9L,
  sum(total_common_accepted$stage == "common_sample_1") == 8L,
  sum(
    total_common_accepted$stage == "common_sample_remedy_1"
  ) == 1L,
  all(file.exists(total_common_accepted$fit_path)),
  all(file.info(total_common_accepted$fit_path)$size > 0),
  all(nchar(total_common_accepted$sha256) == 64L),
  all(total_common_loo_preflight$checks$passed),
  isTRUE(total_common_loo_preflight$ready)
)
```

This checkpoint keeps the eight accepted first-pass common-sample fits and
the model-11 remedy, verifies that all nine use the same 421 hypertrophy
observations, and does not sample. Once reviewed, calculate only their nine
replacement integrated observation-level PSIS-LOO reductions. Each branch
uses one core internally; two concurrent branches limit memory pressure from
the large fit archives:

```r
Sys.setenv(
  PROTEIN_MODEL_WORKERS = "2",
  PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_LOO = "run"
)

targets::tar_make_future(
  update_2026_total_primary_common_sample_integrated_loo_audit_files,
  workers = 2,
  reporter = "balanced"
)

Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_LOO",
  "PROTEIN_MODEL_WORKERS"
))
```

The replacement LOO files are written separately under
`model_archive/update_2026/absolute_exposure/total_g_day/primary/loo/common_sample_final/`.
The audit verifies that all nine pointwise vectors contain the identical 421
observations and stops before stacking. Any Pareto-k values above 0.7 must be
handled through a later exact-refit checkpoint; this command cannot start an
exact refit or overwrite the original 429-row LOO lineage.

The completed replacement audit contains nine aligned 421-observation LOO
objects, a maximum Pareto-k of 0.405, and no cases above 0.7. No exact refit is
therefore required. With every fitting switch disabled, assemble the final
72-model fit and LOO registries and calculate the four stacking vectors:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_LOO",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_LOO",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_stacking_checkpoint_files,
  reporter = "balanced"
)

total_final_loo <- targets::tar_read(
  update_2026_total_primary_final_loo_audit
)
total_weights <- targets::tar_read(
  update_2026_total_primary_stacking_weights
)

total_final_loo$outcome_summary |>
  print(n = Inf, width = Inf)

total_final_loo$final_loo_manifest |>
  dplyr::count(outcome, loo_source) |>
  print(n = Inf)

total_weights |>
  dplyr::group_by(outcome) |>
  dplyr::summarise(
    models = dplyr::n(),
    weight_sum = sum(weight),
    nonzero_weights = sum(weight > 0),
    maximum_weight = max(weight),
    .groups = "drop"
  ) |>
  print(n = Inf, width = Inf)

stopifnot(
  nrow(total_final_loo$model_summary) == 72L,
  nrow(total_final_loo$final_loo_manifest) == 72L,
  total_final_loo$exact_refits == 0L,
  isTRUE(total_final_loo$stacking_ready),
  all(total_final_loo$outcome_summary$models == 18L),
  all(total_final_loo$outcome_summary$high_k_cases == 0L),
  all(total_final_loo$outcome_summary$observation_sets_aligned),
  sum(
    total_final_loo$final_loo_manifest$loo_source ==
      "common_sample_integrated_psis"
  ) == 9L,
  nrow(total_weights) == 72L,
  all(
    abs(tapply(total_weights$weight, total_weights$outcome, sum) - 1) <
      1e-8
  )
)
```

This stage performs no sampling. It replaces only the nine superseded
hypertrophy fits and LOO reductions, retains the other 63 accepted entries,
and refuses to stack unless all 18 models within every outcome use identical
pointwise observations. The model-11 convergence remedy retained 48,000
posterior draws versus 24,000 for the other candidates. Because stacking uses
the already-computed pointwise ELPD values, the implementation standardises
only the draw-count metadata on temporary in-memory LOO copies before calling
`loo_model_weights()`; it does not truncate draws or modify any archived LOO
object.

After the stacking checkpoint passes, build the sampling-free reporting
preflight. It allocates 4,000 posterior-mixture draws per outcome, uses a
50 g/day lower grid limit, adds the audited 78.785 g/day reference explicitly,
and advances in 1 g/day steps through each outcome's observed maximum:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REPORTING",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_reporting_preflight_file,
  reporter = "balanced"
)

total_reporting_preflight <- targets::tar_read(
  update_2026_total_primary_reporting_preflight
)

total_reporting_preflight |>
  dplyr::summarise(
    selected_models = dplyr::n(),
    minimum_draws = min(allocated_draws),
    maximum_draws = max(allocated_draws),
    allocated_draws = sum(allocated_draws),
    retained_weight = sum(weight),
    reference_g_day = dplyr::first(reference_g_day),
    exposure_minimum_g_day = dplyr::first(exposure_minimum_g_day),
    exposure_maximum_g_day = dplyr::first(exposure_maximum_g_day),
    all_fits_exist = all(fit_exists),
    all_valid = all(valid),
    .by = outcome
  ) |>
  print(n = Inf, width = Inf)

stopifnot(
  all(total_reporting_preflight$valid),
  all(
    tapply(
      total_reporting_preflight$allocated_draws,
      total_reporting_preflight$outcome,
      sum
    ) == 4000L
  ),
  all(total_reporting_preflight$reference_g_day == 78.785),
  all(total_reporting_preflight$exposure_minimum_g_day == 50),
  all(total_reporting_preflight$grid_step_g_day == 1)
)
```

Once reviewed, enable only posterior prediction and reporting. Two concurrent
branches limit memory pressure while opening the large accepted fit archives;
each branch uses one core and performs no Stan sampling:

```r
Sys.setenv(
  PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REPORTING = "run",
  PROTEIN_MODEL_WORKERS = "2"
)

targets::tar_make_future(
  update_2026_total_primary_reporting_audit_files,
  workers = 2,
  reporter = "balanced"
)

Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REPORTING",
  "PROTEIN_MODEL_WORKERS"
))
```

Prediction draws are written per selected model so completed branches remain
resumable. Each prediction grid inherits grouping-factor levels and moderator
coding from that accepted fit's own training data. This is required for the
hypertrophy common-sample replacements, whose study-factor levels deliberately
exclude observations removed by the common-sample alignment; the exposure
range still comes from the full outcome dataset. The four model-average files
retain 4,000 draws each and report absolute curves, changes from 78.785 g/day,
pairwise contrasts, support, and local-equivalence summaries in grams per day.

After all four reporting audits pass, freeze the total-intake primary family in
an immutable verification bundle. This is sampling-free and reads only the
accepted fit manifests and compact diagnostic, LOO, prediction, and reporting
artifacts:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_LOO",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_REMEDY",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_COMMON_SAMPLE_LOO",
  "PROTEIN_UPDATE_2026_TOTAL_PRIMARY_REPORTING",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_primary_verification_files,
  reporter = "balanced"
)

total_primary_verification <- targets::tar_read(
  update_2026_total_primary_verification_checkpoint
)
total_primary_verification_files <- targets::tar_read(
  update_2026_total_primary_verification_files
)

total_primary_verification$checks |>
  print(n = Inf, width = Inf)

c(
  files = length(total_primary_verification_files),
  all_files_exist = all(file.exists(total_primary_verification_files)),
  all_files_nonempty =
    all(file.info(total_primary_verification_files)$size > 0)
)

stopifnot(
  nrow(total_primary_verification$checks) == 14L,
  all(total_primary_verification$checks$passed),
  length(total_primary_verification_files) == 39L,
  all(file.exists(total_primary_verification_files)),
  all(file.info(total_primary_verification_files)$size > 0)
)
```

The bundle freezes the 72-model aligned LOO lineage, four 18-model stacking
vectors, 16,000 allocated prediction draws, all four reporting files, and the
shared absolute-exposure provenance. It contains 37 audit artifacts plus its
manifest and README. Only after that checkpoint is reviewed should the
sampling-free lean-pruning plan be built:

```r
targets::tar_make(
  update_2026_total_primary_pruning_plan_files,
  reporter = "balanced"
)

total_primary_pruning <- targets::tar_read(
  update_2026_total_primary_pruning_plan
)

total_primary_pruning$summary |>
  print(width = Inf)

total_primary_pruning$fit_file_manifest |>
  dplyr::count(archive_stage, accepted_fit) |>
  print(n = Inf)

total_primary_pruning$compact_archive_manifest |>
  dplyr::count(archive_component) |>
  print(n = Inf)

stopifnot(
  nrow(total_primary_pruning$fit_file_manifest) == 83L,
  sum(total_primary_pruning$fit_file_manifest$accepted_fit) == 72L,
  sum(!total_primary_pruning$fit_file_manifest$accepted_fit) == 11L,
  nrow(total_primary_pruning$deletion_file_manifest) == 83L,
  all(file.exists(total_primary_pruning$deletion_file_manifest$path)),
  nrow(total_primary_pruning$compact_archive_manifest) == 113L,
  all(file.exists(total_primary_pruning$compact_archive_manifest$path)),
  all(nchar(total_primary_pruning$compact_archive_manifest$sha256) == 64L),
  file.exists(total_primary_pruning$verification_manifest)
)
```

The plan deletes nothing. Its audited lean policy proposes deleting the 72
accepted and 11 superseded fit archives while retaining 81 compact LOO files,
32 prediction/reporting files, and the complete verification bundle. Do not
execute it until every row has been reviewed and the irreversible loss of the
fit objects has been explicitly confirmed. After that separate confirmation,
the guarded executor requires the exact phrase below and writes a checksum
receipt before deleting only the audited fit files:

```r
total_primary_pruning_receipt <-
  execute_update_2026_total_primary_pruning(
    plan = total_primary_pruning,
    confirmation = "DELETE ALL TOTAL-PRIMARY FIT ARCHIVES"
  )
```

## Absolute total-intake Mundlak family

The next total-intake sensitivity family is a fresh 72-model `brms` Mundlak
analysis: 18 retained model forms for each of four outcomes. The eight fixed
1.6 g/kg/day spline models remain omitted because they have no unique g/day
knot. Model 4 retains the prospectively accepted variance regularisation. Five
models inherit longer sampling from the accepted relative-intake Mundlak
lineage: whole-body models 4 and 6, strength models 5 and 6, and function model
5. The other 67 models use the base Mundlak schedule.

All 18 candidates within an outcome use the same observation set from the
outset. In particular, all hypertrophy models use the 421-row, 53-study
body-mass-complete sample. This preserves pointwise comparability for later
stacking and avoids repeating the nine common-sample replacement fits needed
by the total-intake primary family.

First build and inspect the sampling-free preflight. It constructs all 72 Stan
data lists through `brms` but does not compile or sample:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_mundlak_input_preflight_file,
  reporter = "balanced"
)

total_mundlak_preflight <- targets::tar_read(
  update_2026_total_mundlak_input_preflight
)

total_mundlak_preflight |>
  dplyr::count(
    outcome,
    historical_schedule,
    contrasts,
    studies
  ) |>
  print(n = Inf, width = Inf)

total_mundlak_preflight |>
  dplyr::filter(historical_schedule != "base") |>
  dplyr::select(
    outcome,
    model_index,
    model_name,
    historical_schedule,
    adapt_delta,
    max_treedepth,
    iter_warmup,
    iter_sampling,
    chains,
    valid
  ) |>
  print(n = Inf, width = Inf)

schedule_counts <- table(
  total_mundlak_preflight$historical_schedule
)

stopifnot(
  nrow(total_mundlak_preflight) == 72L,
  all(total_mundlak_preflight$valid),
  all(total_mundlak_preflight$common_sample),
  !anyDuplicated(paste(
    total_mundlak_preflight$outcome,
    total_mundlak_preflight$model_index
  )),
  all(table(total_mundlak_preflight$outcome) == 18L),
  identical(unname(schedule_counts["base"]), 67L),
  identical(unname(schedule_counts["longer_sampling"]), 5L),
  all(
    total_mundlak_preflight$contrasts[
      total_mundlak_preflight$outcome == "hypertrophy"
    ] == 421L
  ),
  sum(
    total_mundlak_preflight$model_index == 4L &
      total_mundlak_preflight$regularised_variance_priors
  ) == 4L,
  all(total_mundlak_preflight$chains == 4L),
  all(total_mundlak_preflight$parallel_chains == 4L)
)
```

After that checkpoint is reviewed, enable only this family and fit up to four
candidates concurrently. Each model runs four chains internally and writes its
fit to an independent resumable file:

```r
Sys.setenv(
  PROTEIN_UPDATE_2026_TOTAL_MUNDLAK = "fit",
  PROTEIN_MODEL_WORKERS = "4"
)

targets::tar_make_future(
  update_2026_total_mundlak_diagnostic_audit_files,
  workers = 4,
  reporter = "balanced"
)

Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK",
  "PROTEIN_MODEL_WORKERS"
))
```

This first run stops at the 72-model sampling-diagnostic audit. Remedies, LOO,
stacking, reporting, verification, and pruning remain separate later gates.

The initial audit identified four convergence-only failures: whole-body model
2, strength models 2 and 17, and function model 6. All four passed the geometry
checks with no divergences or maximum-treedepth transitions. Their first remedy
therefore preserves the `brms` model and sampler controls while increasing each
fit to 2,000 warmup and 4,000 retained iterations per chain.

Build and inspect the remedy preflight without sampling:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK",
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_REMEDY",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_mundlak_remedy_preflight_file,
  reporter = "balanced"
)

total_mundlak_remedy_preflight <- targets::tar_read(
  update_2026_total_mundlak_remedy_preflight
)

total_mundlak_remedy_preflight |>
  print(n = Inf, width = Inf)

remedy_keys <- paste(
  total_mundlak_remedy_preflight$outcome,
  total_mundlak_remedy_preflight$model_index,
  sep = "::"
)

stopifnot(
  nrow(total_mundlak_remedy_preflight) == 4L,
  setequal(
    remedy_keys,
    c("whole_body::2", "strength::2", "strength::17", "function::6")
  ),
  all(total_mundlak_remedy_preflight$initial_fit_exists),
  all(total_mundlak_remedy_preflight$valid),
  all(total_mundlak_remedy_preflight$iter_warmup == 2000L),
  all(total_mundlak_remedy_preflight$iter_sampling == 4000L),
  all(total_mundlak_remedy_preflight$chains == 4L),
  all(total_mundlak_remedy_preflight$parallel_chains == 4L)
)
```

After that checkpoint passes, enable only the four remedy branches:

```r
Sys.setenv(
  PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_REMEDY = "fit",
  PROTEIN_MODEL_WORKERS = "4"
)

targets::tar_make_future(
  update_2026_total_mundlak_final_diagnostic_audit_files,
  workers = 4,
  reporter = "balanced"
)

Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_REMEDY",
  "PROTEIN_MODEL_WORKERS"
))
```

When all four remedies pass, materialise the accepted-fit registry and the
sampling-free LOO checkpoint:

```r
targets::tar_make(
  update_2026_total_mundlak_loo_preflight_files,
  reporter = "balanced"
)

total_mundlak_final <- targets::tar_read(
  update_2026_total_mundlak_final_diagnostic_audit
)
total_mundlak_accepted <- targets::tar_read(
  update_2026_total_mundlak_accepted_fit_manifest
)
total_mundlak_loo_preflight <- targets::tar_read(
  update_2026_total_mundlak_loo_preflight
)

total_mundlak_accepted |>
  dplyr::count(outcome, stage) |>
  print(n = Inf)

total_mundlak_loo_preflight$checks |>
  print(n = Inf, width = Inf)

stopifnot(
  total_mundlak_final$models == 72L,
  total_mundlak_final$models_passing_sampling == 72L,
  total_mundlak_final$models_requiring_remediation == 0L,
  total_mundlak_final$all_sampling_diagnostics_pass,
  nrow(total_mundlak_accepted) == 72L,
  sum(total_mundlak_accepted$stage == "initial") == 68L,
  sum(total_mundlak_accepted$stage == "remedy_1") == 4L,
  all(file.exists(total_mundlak_accepted$fit_path)),
  all(nchar(total_mundlak_accepted$sha256) == 64L),
  all(total_mundlak_loo_preflight$checks$passed),
  isTRUE(total_mundlak_loo_preflight$ready)
)
```

After the LOO checkpoint passes, enable only the 72 integrated PSIS-LOO
branches. Each branch uses 4,000 evenly spaced posterior draws and one core;
the outer pipeline may run four independent models concurrently:

```r
Sys.setenv(
  PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_LOO = "run",
  PROTEIN_MODEL_WORKERS = "4"
)

targets::tar_make_future(
  update_2026_total_mundlak_integrated_loo_audit_files,
  workers = 4,
  reporter = "balanced"
)

Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_LOO",
  "PROTEIN_MODEL_WORKERS"
))
```

This gate ends at the integrated LOO audit. Any observations with Pareto
`k > 0.7` are handled by a later exact-deletion stage; stacking is not enabled
until all such cases are resolved and the within-outcome observation sets are
confirmed identical.

The integrated audit identified two exact deletions, both for hypertrophy model
20: observations 248 and 249 from study 124, arm 189, outcome 1023. Build and
inspect their sampling-free preflight first:

```r
Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_LOO",
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_EXACT",
  "PROTEIN_MODEL_WORKERS"
))

targets::tar_make(
  update_2026_total_mundlak_exact_preflight_file,
  reporter = "balanced"
)

total_mundlak_exact_preflight <- targets::tar_read(
  update_2026_total_mundlak_exact_preflight
)

total_mundlak_exact_preflight |>
  print(n = Inf, width = Inf)

stopifnot(
  nrow(total_mundlak_exact_preflight) == 2L,
  identical(total_mundlak_exact_preflight$case_id, 1:2),
  all(total_mundlak_exact_preflight$outcome == "hypertrophy"),
  all(total_mundlak_exact_preflight$model_index == 20L),
  setequal(total_mundlak_exact_preflight$observation, c(248L, 249L)),
  all(total_mundlak_exact_preflight$pareto_k > 0.7),
  all(total_mundlak_exact_preflight$iter_warmup == 2000L),
  all(total_mundlak_exact_preflight$iter_sampling == 6000L),
  all(total_mundlak_exact_preflight$original_fit_exists),
  all(total_mundlak_exact_preflight$integrated_loo_exists),
  all(total_mundlak_exact_preflight$valid)
)
```

After review, enable only the two exact-refit branches. Run the outer branches
sequentially because each refit already opens four internal chain workers. A
valid existing exact archive is reused on restart. This gate archives the
refits and checks their sampling diagnostics; it does not yet alter the final
LOO vectors:

```r
Sys.setenv(
  PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_EXACT = "fit",
  PROTEIN_MODEL_WORKERS = "1"
)

targets::tar_make(
  update_2026_total_mundlak_exact_diagnostic_audit_file,
  reporter = "balanced"
)

Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_EXACT",
  "PROTEIN_MODEL_WORKERS"
))
```

Once both exact refits pass, enable the post-fit correction gate. This opens
the completed refits sequentially, writes two compact exact corrections,
patches the final 72 LOO vectors, and calculates the four stacking vectors; it
does not sample any model:

```r
Sys.unsetenv("PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_EXACT")
Sys.setenv(
  PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_CORRECTION = "run",
  PROTEIN_MODEL_WORKERS = "1"
)

targets::tar_make(
  update_2026_total_mundlak_stacking_weights_file,
  reporter = "balanced"
)

Sys.unsetenv(c(
  "PROTEIN_UPDATE_2026_TOTAL_MUNDLAK_CORRECTION",
  "PROTEIN_MODEL_WORKERS"
))

total_mundlak_final_loo <- targets::tar_read(
  update_2026_total_mundlak_final_loo_audit
)
total_mundlak_weights <- targets::tar_read(
  update_2026_total_mundlak_stacking_weights
)

total_mundlak_final_loo$outcome_summary |>
  print(n = Inf, width = Inf)

total_mundlak_weights |>
  dplyr::summarise(
    models = dplyr::n(),
    weight_sum = sum(weight),
    nonzero_weights = sum(weight > 0),
    maximum_weight = max(weight),
    .by = outcome
  ) |>
  print(n = Inf, width = Inf)

stopifnot(
  total_mundlak_final_loo$exact_refits == 2L,
  total_mundlak_final_loo$all_exact_refit_diagnostics_pass,
  isTRUE(total_mundlak_final_loo$stacking_ready),
  nrow(total_mundlak_final_loo$model_summary) == 72L,
  all(
    total_mundlak_final_loo$outcome_summary$
      unresolved_high_k_cases == 0L
  ),
  nrow(total_mundlak_weights) == 72L,
  all(abs(
    tapply(
      total_mundlak_weights$weight,
      total_mundlak_weights$outcome,
      sum
    ) - 1
  ) < 1e-8)
)
```

## Local combined-data Mundlak runbook

The Mundlak update is independent of the frozen sensitivity targets and never
reuses their seven pilot fits. It fits all 104 outcome-model combinations on
the combined data. It starts from the final accepted primary specifications,
including the accepted `adapt_delta = 0.9995` setting for function model 10,
and retains the prospectively selected model-4 variance regularisation from
the original Mundlak workflow.

This first hand-off deliberately ends at the diagnostic review gate. Run it
from an R session at the project root, with no more than two model workers.

1. Clear any old fitting switches and run the non-sampling preflight:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_PRIMARY",
     "PROTEIN_UPDATE_2026_PRIMARY_REMEDY",
     "PROTEIN_UPDATE_2026_PRIMARY_LOO",
     "PROTEIN_UPDATE_2026_PRIMARY_EXACT",
     "PROTEIN_UPDATE_2026_PRIMARY_REPORTING",
     "PROTEIN_SENSITIVITY_PILOTS",
     "PROTEIN_UPDATE_2026_MUNDLAK"
   ))
   Sys.setenv(PROTEIN_MODEL_WORKERS = "2")
   targets::tar_make(
     update_2026_mundlak_input_preflight_file,
     reporter = "balanced"
   )
   ```

   Confirm that
   `results/update_2026/sensitivity/mundlak/input_preflight.csv` has 104 rows,
   all `valid` values are `TRUE`, four model-4 rows have regularised variance
   priors, and function model 10 has `adapt_delta` 0.9995.

2. Enable only the combined-data Mundlak fits and run at most two models at a
   time:

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_MUNDLAK = "fit"
   )
   targets::tar_make_future(
     update_2026_mundlak_fit_file,
     workers = 2,
     reporter = "balanced"
   )
   ```

   This is resumable. Successful branches are stored beneath
   `model_archive/update_2026/sensitivity/mundlak/candidates/`, and the target
   store retains paths rather than duplicate posterior objects.

3. When all 104 fits finish, build the diagnostic audit with the same
   two-worker ceiling:

   ```r
   targets::tar_make_future(
     update_2026_mundlak_diagnostic_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   Review the files beneath
   `results/update_2026/sensitivity/mundlak/diagnostics/` and stop. The next
   remedy stage will be prepared only from these combined-data results. Do not
   carry forward the eight historical Mundlak remedies automatically.

4. The combined-data audit identified eight sampling failures: whole-body
   models 3, 4, 6, and 13; strength models 5, 6, and 25; and function model 5.
   None had divergences. Seven receive only the longer historical remedy
   schedule of 2,000 warmup and 4,000 retained iterations per chain. Strength
   model 25 additionally receives `max_treedepth = 15` because it recorded 11
   saturated transitions. First run the non-sampling remedy preflight:

   ```r
   Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK")
   targets::tar_make(
     update_2026_mundlak_remedy_preflight_file,
     reporter = "balanced"
   )
   ```

   Confirm that the preflight contains exactly those eight models, all rows
   are valid, and only strength model 25 has the tree-depth remedy.

5. Enable only the eight remedy branches, retaining the two-worker ceiling:

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_MUNDLAK_REMEDY = "fit"
   )
   targets::tar_make_future(
     update_2026_mundlak_remedy_diagnostic_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   The replacement fits are stored beneath
   `model_archive/update_2026/sensitivity/mundlak/remedies/remedy_1/`.
   Review the final 104-model audit beneath
   `results/update_2026/sensitivity/mundlak/remedy_1/`, then disable the switch
   with `Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_REMEDY")`. Do not proceed to
   LOO unless all eight replacements pass.

6. Seven first remedies passed. Strength model 6 retained three saturated
   transitions at `max_treedepth = 12`, while its convergence, ESS, MCSE, and
   divergence checks passed. Its second remedy therefore keeps the 2,000/
   4,000 schedule and changes only `max_treedepth` to 15. Run its non-sampling
   preflight with the earlier switches disabled:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_MUNDLAK",
     "PROTEIN_UPDATE_2026_MUNDLAK_REMEDY"
   ))
   targets::tar_make(
     update_2026_mundlak_second_remedy_preflight_file,
     reporter = "balanced"
   )
   ```

7. After confirming the single valid strength-model-6 row, enable only its
   second remedy and retain the two-worker ceiling (only one branch will run):

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_MUNDLAK_SECOND_REMEDY = "fit"
   )
   targets::tar_make_future(
     update_2026_mundlak_second_remedy_diagnostic_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   Review the final replacement-aware audit beneath
   `results/update_2026/sensitivity/mundlak/remedy_2/`, then disable the switch
   with `Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_SECOND_REMEDY")`. LOO remains
   blocked until the final 104-model diagnostic audit passes.

8. After all 104 final diagnostics pass, create the authoritative accepted-fit
   manifest. This hashes but does not refit the 96 accepted initial files,
   seven accepted first remedies, and one accepted second remedy:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_MUNDLAK",
     "PROTEIN_UPDATE_2026_MUNDLAK_REMEDY",
     "PROTEIN_UPDATE_2026_MUNDLAK_SECOND_REMEDY"
   ))
   targets::tar_make(
     update_2026_mundlak_accepted_fit_manifest_file,
     reporter = "balanced"
   )
   ```

   Confirm that the manifest contains 104 unique outcome-model keys and stage
   counts of 96 initial, seven `remedy_1`, and one `remedy_2` fit.

9. Enable the integrated observation-level PSIS-LOO stage only after the
   accepted manifest passes. Each branch opens one accepted fit, uses up to
   4,000 draws, and writes one compact reduction; no posterior fit is copied.

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_MUNDLAK_LOO = "run"
   )
   targets::tar_make_future(
     update_2026_mundlak_integrated_loo_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   Review `outcome_summary.csv`, `model_summary.csv`, and `high_k_cases.csv`
   beneath `results/update_2026/sensitivity/mundlak/loo_psis/`, then disable
   the switch with `Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_LOO")`. Stop at
   this gate: exact case-deletion refits must be generated from the new Pareto-k
   audit rather than inherited from the frozen Mundlak analysis.

10. The combined-data integrated audit identified 15 exact deletions: nine
    whole-body cases across eight models and six hypertrophy cases across four
    models. Function and strength require none. First verify all held-out rows,
    accepted-fit paths, integrated reductions, and the 2,000-warmup/6,000-
    sampling schedule without fitting:

    ```r
    Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_LOO")
    targets::tar_make(
      update_2026_mundlak_exact_preflight_file,
      reporter = "balanced"
    )
    ```

11. Once all 15 preflight rows are valid, enable only the exact stage and keep
    the two-worker ceiling:

    ```r
    Sys.setenv(
      PROTEIN_MODEL_WORKERS = "2",
      PROTEIN_UPDATE_2026_MUNDLAK_EXACT = "fit"
    )
    targets::tar_make_future(
      update_2026_mundlak_exact_audit_files,
      workers = 2,
      reporter = "balanced"
    )
    ```

    Each branch refits one deletion, calculates its exact integrated-LOO
    correction, and discards the temporary `brmsfit`. Only compact corrections
    are retained beneath
    `model_archive/update_2026/sensitivity/mundlak/loo/exact/`. Review the
    diagnostic audit beneath
    `results/update_2026/sensitivity/mundlak/loo_exact/`, then disable the
    switch with `Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_EXACT")`. Exact
    corrections must not be applied or used for stacking until every refit
    diagnostic passes.

12. Once all 15 exact diagnostics pass, assemble the corrected 104-model LOO
    registry and calculate four 26-model stacking vectors. This is a
    non-sampling step that reads only compact integrated reductions and exact
    corrections:

    ```r
    Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_EXACT")
    targets::tar_make(
      update_2026_mundlak_stacking_checkpoint_files,
      reporter = "balanced"
    )
    ```

    Final LOO tables are written beneath
    `results/update_2026/sensitivity/mundlak/loo_final/`, and weights are
    written to
    `results/update_2026/sensitivity/mundlak/stacking_weights.csv`. Confirm
    that no high-k cases remain, every outcome contains 26 models, and each
    outcome's weights sum to one before preparing model-averaged reporting.

13. Build the reporting allocation preflight. This allocates 4,000 mixture
    draws per outcome according to the corrected stacking weights and checks
    every selected accepted-fit path without opening a model:

    ```r
    targets::tar_make(
      update_2026_mundlak_reporting_preflight_file,
      reporter = "balanced"
    )
    ```

    Confirm four outcomes, 4,000 allocated draws per outcome, at least one
    selected model per outcome, and no invalid paths.

14. Enable the storage-aware reporting stage with two workers. Each selected
    model is opened in its own branch and writes only its allocated prediction
    draws. Four outcome-level branches then create compact 4,000-draw mixtures
    and reporting tables:

    ```r
    Sys.setenv(
      PROTEIN_MODEL_WORKERS = "2",
      PROTEIN_UPDATE_2026_MUNDLAK_REPORTING = "run"
    )
    targets::tar_make_future(
      update_2026_mundlak_reporting_audit_files,
      workers = 2,
      reporter = "balanced"
    )
    ```

    Compact reporting artifacts are stored beneath
    `model_archive/update_2026/sensitivity/mundlak/reporting/`; flat audit
    tables are written beneath
    `results/update_2026/sensitivity/mundlak/reporting/`. Review the summary,
    prediction allocation, weights, dose-response tables, contrasts, support,
    and local-equivalence output, then disable the switch with
    `Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_REPORTING")`. Do not prune fits
    until the final family verification bundle passes.

15. Write and inspect the final updated Mundlak verification bundle. This is a
    non-sampling step and does not require a runtime switch:

    ```r
    Sys.unsetenv("PROTEIN_UPDATE_2026_MUNDLAK_REPORTING")
    targets::tar_make(
      update_2026_mundlak_verification_files,
      reporter = "balanced"
    )
    ```

    The bundle is written beneath
    `results/verification/mundlak_update_2026/`. Confirm that every checkpoint
    passes and every bundle file exists before preparing the lean-pruning
    plan. The established Mundlak overlay and increment figures compare the
    primary, Mundlak, and within-study families, so they are deliberately
    deferred until the combined-data within-study refit is complete.

16. Generate the dry-run lean-pruning plan. This validates the verification
    bundle and inventories the complete fit lineage, but deletes nothing:

    ```r
    targets::tar_make(
      update_2026_mundlak_pruning_plan_files,
      reporter = "balanced"
    )
    ```

    The plan must contain exactly 113 full fits: 104 initial candidates, eight
    first remedies, and one second remedy. Of those, 104 are accepted fits and
    nine are superseded fits. Review the summary and exact path manifest beneath
    `results/update_2026/sensitivity/mundlak/pruning/` before proceeding.

17. After reviewing the plan, execute the selected lean policy explicitly from
    R. This is the only step that removes files:

    ```r
    mundlak_pruning_plan <- targets::tar_read(
      update_2026_mundlak_pruning_plan
    )
    source("R/incremental_update.R")
    source("R/update_2026_sensitivity.R")
    execute_update_2026_mundlak_pruning(
      plan = mundlak_pruning_plan,
      confirmation = "DELETE ALL UPDATED MUNDLAK FIT ARCHIVES"
    )
    ```

    The helper revalidates the frozen verification checksums, restricts every
    deletion to the updated Mundlak `candidates/` or `remedies/` subtree,
    records the pre-deletion SHA-256 hash of each fit, preserves the compact
    LOO and reporting archives, and writes `executed_lean.csv` as its receipt.

## Lean retention and pruning rule

Do not prune any updated candidate fit at the first-pass diagnostic stage.
After diagnostic remedies, LOO corrections, stacking, predictions, and the
updated verification bundle have all passed, retain only:

- the compact diagnostics, exact-corrected LOO reductions, stacking weights,
  predictions, reporting objects, specifications, provenance manifests, and
  final verification bundle needed to inspect and reproduce reported
  summaries from retained reductions; and
- the source data and code required for a future full refit.

For updated sensitivity families, all full candidate and remedy fit archives
are intentionally disposable once that gate passes; this includes the 104
accepted fits as well as superseded fits. Transient draw matrices and temporary
exact-refit objects are also disposable. Pruning is prepared as a separate,
explicit local step with an exact-path manifest and receipt; it is never
performed automatically by fitting or verification commands. The updated
primary analysis remains under its separate conservative-retention policy.

## Local combined-data within-study runbook

The strict within-study family is the next updated sensitivity analysis. It
reconstructs protein-arm minus matched comparator-arm contrasts from the
combined arm data and refits all 104 candidates. No frozen within-study fit is
reused. The historically remedied hypertrophy and strength model 9 candidates
start with their accepted 2,000-warmup/4,000-sampling schedule,
`adapt_delta = 0.999`, and `max_treedepth = 15`; all other candidates start at
1,500/2,000 with their final primary specification settings.

1. Rebuild every transformed within-study input and Stan data set without
   compiling or fitting a model:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_MUNDLAK",
     "PROTEIN_UPDATE_2026_MUNDLAK_REMEDY",
     "PROTEIN_UPDATE_2026_MUNDLAK_SECOND_REMEDY",
     "PROTEIN_UPDATE_2026_MUNDLAK_LOO",
     "PROTEIN_UPDATE_2026_MUNDLAK_EXACT",
     "PROTEIN_UPDATE_2026_MUNDLAK_REPORTING",
     "PROTEIN_UPDATE_2026_WITHIN"
   ))
   targets::tar_make(
     update_2026_within_input_preflight_file,
     reporter = "balanced"
   )
   ```

   Confirm 26 valid candidates per outcome, four non-empty contrast data sets,
   and exactly the two historical remedy schedules described above.

2. After the preflight passes, explicitly enable the initial fit stage with no
   more than two simultaneous models:

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_WITHIN = "fit"
   )
   targets::tar_make_future(
     update_2026_within_fit_file,
     workers = 2,
     reporter = "balanced"
   )
   ```

   The 104 full fits are written directly beneath
   `model_archive/update_2026/sensitivity/within/candidates/`. Stop after the
   files complete, disable the switch, and run the diagnostic audit as a
   separate review gate. Remedy, LOO, reporting, verification, and lean-pruning
   targets are added only after their preceding audit determines the required
   work.

3. The initial combined-data diagnostic audit selected four remedies:
   hypertrophy model 10 and strength models 12, 23, and 25. Generate their
   non-sampling preflight first:

   ```r
   targets::tar_make(
     update_2026_within_remedy_preflight_file,
     reporter = "balanced"
   )
   ```

   All four use 2,000 warmup and 4,000 retained iterations. Hypertrophy model
   10 alone raises `adapt_delta` from 0.99 to 0.999 because it had one
   divergence. Existing tree-depth ceilings are retained because none of the
   four initial fits had a saturated transition.

4. After reviewing the remedy preflight, enable only the four remedy branches:

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_WITHIN_REMEDY = "fit"
   )
   targets::tar_make_future(
     update_2026_within_remedy_fit_file,
     workers = 2,
     reporter = "balanced"
   )
   ```

   Remedy fits are written beneath
   `model_archive/update_2026/sensitivity/within/remedies/remedy_1/`.
   Disable the remedy switch after all four files exist, then run the combined
   104-model remedy diagnostic audit as the next review gate.

5. Once the combined diagnostic audit contains 104 passing models, create and
   inspect the accepted-fit registry. This hashes the 100 accepted initial fits
   and four accepted remedies but performs no sampling:

   ```r
   targets::tar_make(
     update_2026_within_accepted_fit_manifest_file,
     reporter = "balanced"
   )
   ```

6. After the accepted registry passes, run storage-aware integrated PSIS-LOO
   with no more than two model workers:

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_WITHIN_LOO = "run"
   )
   targets::tar_make_future(
     update_2026_within_integrated_loo_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   Each branch opens one accepted fit and writes only its compact integrated
   reduction beneath
   `model_archive/update_2026/sensitivity/within/loo/integrated/`. Disable the
   switch after completion and review all Pareto-k cases before exact refits
   are prepared.

7. The integrated audit selected 59 exact observation deletions: six
   whole-body, 26 hypertrophy, one strength, and 26 function cases. Build the
   exact-refit preflight without sampling:

   ```r
   targets::tar_make(
     update_2026_within_exact_preflight_file,
     reporter = "balanced"
   )
   ```

8. After confirming all 59 records, run the exact refits with two workers:

   ```r
   Sys.setenv(
     PROTEIN_MODEL_WORKERS = "2",
     PROTEIN_UPDATE_2026_WITHIN_EXACT = "fit"
   )
   targets::tar_make_future(
     update_2026_within_exact_audit_files,
     workers = 2,
     reporter = "balanced"
   )
   ```

   Each deletion uses 2,000 warmup and 6,000 retained iterations. The
   temporary refitted model is discarded immediately after its exact score and
   diagnostics are calculated; only compact corrections are retained beneath
   `model_archive/update_2026/sensitivity/within/loo/exact/`. Disable the exact
   switch after completion. Do not apply corrections or compute stacking
   weights until every exact-refit diagnostic passes.

9. The first exact-refit audit produced one divergence in whole-body model 21,
   observation 2 (case 6); the other 58 exact refits passed. Build the isolated
   remedy preflight:

   ```r
   targets::tar_make(
     update_2026_within_exact_remedy_preflight_file,
     reporter = "balanced"
   )
   ```

   The remedy retains the 2,000/6,000 schedule and tree-depth ceiling of 12,
   raises only `adapt_delta` from 0.99 to 0.999, and uses a fresh seed.

10. Enable and run only that exact-remedy branch:

    ```r
    Sys.setenv(
      PROTEIN_MODEL_WORKERS = "2",
      PROTEIN_UPDATE_2026_WITHIN_EXACT_REMEDY = "fit"
    )
    targets::tar_make_future(
      update_2026_within_exact_remedy_audit_files,
      workers = 2,
      reporter = "balanced"
    )
    ```

    Disable the switch after completion. The remedy again retains only its
    compact correction; final LOO correction and stacking remain blocked until
    its sampling diagnostics pass.

11. Once the exact remedy passes, replace only case 6 in the compact
    correction registry, assemble the final exact-corrected LOO results, and
    calculate four 26-model stacking vectors:

    ```r
    Sys.unsetenv("PROTEIN_UPDATE_2026_WITHIN_EXACT_REMEDY")
    targets::tar_make(
      update_2026_within_stacking_checkpoint_files,
      reporter = "balanced"
    )
    ```

    This step performs no sampling and does not open the full fit archives.
    Confirm that all 59 final exact diagnostics pass, no unresolved high-k
    cases remain, and each outcome's 26 stacking weights sum to one before
    preparing model-averaged reporting.

12. Build the within-study reporting allocation without opening any fitted
    model:

    ```r
    targets::tar_make(
      update_2026_within_reporting_preflight_file,
      reporter = "balanced"
    )
    ```

    The preflight deterministically allocates 4,000 posterior-mixture draws
    per outcome from the four 26-model stacking vectors. Only models allocated
    at least one draw proceed to prediction, and every accepted fit path must
    exist before the reporting stage can be enabled.

13. After confirming the preflight, generate compact prediction reductions,
    four posterior mixtures, and the audited reporting files with two workers:

    ```r
    Sys.setenv(
      PROTEIN_MODEL_WORKERS = "2",
      PROTEIN_UPDATE_2026_WITHIN_REPORTING = "run"
    )
    targets::tar_make_future(
      update_2026_within_reporting_audit_files,
      workers = 2,
      reporter = "balanced"
    )
    ```

    Each selected full fit is opened only for its allocated predictions. The
    retained artifacts are compact model-specific draws, four 4,000-draw
    mixtures, reporting objects, and flat audit tables under
    `model_archive/update_2026/sensitivity/within/reporting/` and
    `results/update_2026/sensitivity/within/reporting/`. Disable the reporting
    switch after completion.

14. Freeze the complete within-study lineage in its verification bundle:

    ```r
    Sys.unsetenv("PROTEIN_UPDATE_2026_WITHIN_REPORTING")
    targets::tar_make(
      update_2026_within_verification_files,
      reporter = "balanced"
    )
    ```

    This lightweight checkpoint verifies the immutable inputs, contrast data,
    104-model registry, accepted-fit lineage, final diagnostics, all 59 exact
    corrections, corrected LOO objects, stacking weights, prediction
    allocation, and reporting outputs. It must pass before cross-family
    figures are built or any full within-study fit archive is pruned.

15. Build the updated cross-family figures and their exact plotted-data
    archive:

    ```r
    targets::tar_make(
      update_2026_primary_mundlak_within_figure_files,
      reporter = "balanced"
    )
    ```

    This opens only the four compact primary posterior mixtures and the
    compact Mundlak/within reporting objects. It creates an absolute
    primary-versus-Mundlak overlay and a 0.8 g/kg/day-referenced
    primary/Mundlak/strict-within overlay restricted to the observed
    within-study treatment-intake support. PNG and PDF figures, plotted CSV
    data, validation checks, sizes, and SHA-256 hashes are written beneath
    `results/update_2026/`; the earlier frozen figures are not overwritten.

16. Generate the within-study dry-run lean-pruning plan. This validates the
    verification bundle and cross-family figure archive, and inventories the
    complete fit lineage, but deletes nothing:

    ```r
    targets::tar_make(
      update_2026_within_pruning_plan_files,
      reporter = "balanced"
    )
    ```

    The plan must contain exactly 108 full fits: 104 initial candidates and
    four first remedies. Of those, 104 are accepted fits and four are
    superseded initial fits. Review the summary and exact path manifest beneath
    `results/update_2026/sensitivity/within/pruning/` before proceeding.

17. After reviewing the plan, execute the selected lean policy explicitly
    from R. This is the only step that removes files:

    ```r
    within_pruning_plan <- targets::tar_read(
      update_2026_within_pruning_plan
    )
    source("R/incremental_update.R")
    source("R/update_2026_sensitivity.R")
    execute_update_2026_within_pruning(
      plan = within_pruning_plan,
      confirmation = "DELETE ALL UPDATED WITHIN FIT ARCHIVES"
    )
    ```

    The helper revalidates the frozen verification bundle and figure hashes,
    restricts every deletion to the updated within-study `candidates/` or
    `remedies/` subtree, records each pre-deletion SHA-256 hash, retains the
    compact LOO and reporting archives, and writes `executed_lean.csv` as its
    receipt.

## Updated protein-manipulation strategy audit

The completed 109-study strategy audit remains unchanged at
`data/protein_manipulation_strategy_audit.csv`. The six studies newly included
in the combined 2026 analysis are recorded separately in
`protein_manipulation_strategy_audit_additions.csv`. Each addition identifies
its Zotero storage key, PDF filename, and SHA-256 hash and records a short
source passage and page. All six classifications have now been author-reviewed
and have `review_status = Confirmed`.

Build the sampling-free review checkpoint with:

```r
targets::tar_make(
  update_2026_strategy_audit_review_files,
  reporter = "balanced"
)
```

The checkpoint combines the historical and updated records, confirms that all
95 studies in the updated analysis are covered, and writes the six-row review
table beneath `results/update_2026/sensitivity/strategy/audit/`. It requires
`audit_ready = TRUE` before any updated strategy model is fitted.

## Updated protein-manipulation strategy family

The updated family refits all 104 candidates on the combined arm data using
the standard `brms::brm()` fitting path. There is no handwritten Stan model in
this family. The sampling-free preflight calls `brms::make_standata()` only to
validate every formula, prior set, transformed dataset, and grouping structure
before expensive fitting is enabled.

The accepted historical strategy checkpoint supplies sampling schedules only;
no historical fit is reused. Twelve candidates that previously needed longer
sampling start at 2,000 warmup and 4,000 retained iterations with
`adapt_delta >= 0.999` and `max_treedepth >= 15`. The other 92 candidates use
1,500 warmup and 2,000 retained iterations. Every candidate runs four chains.

1. Rebuild the confirmed audit and run the sampling-free 104-model preflight:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY",
     "PROTEIN_MODEL_WORKERS"
   ))

   targets::tar_make(
     update_2026_strategy_audit_review_files,
     reporter = "balanced"
   )
   targets::tar_make(
     update_2026_strategy_input_preflight_file,
     reporter = "balanced"
   )
   ```

   Inspect the checkpoint before sampling:

   ```r
   strategy_review <- targets::tar_read(
     update_2026_strategy_audit_review
   )
   strategy_preflight <- targets::tar_read(
     update_2026_strategy_input_preflight
   )

   strategy_review$structural_checks |>
     print(n = Inf, width = Inf)

   strategy_preflight |>
     dplyr::count(outcome, historical_schedule) |>
     print(n = Inf)

   stopifnot(
     isTRUE(strategy_review$audit_ready),
     nrow(strategy_preflight) == 104L,
     all(strategy_preflight$valid),
     all(table(strategy_preflight$outcome) == 26L),
     sum(strategy_preflight$historical_schedule == "base") == 92L,
     sum(
       strategy_preflight$historical_schedule == "longer_sampling"
     ) == 12L,
     all(strategy_preflight$chains == 4L),
     all(strategy_preflight$parallel_chains == 4L)
   )
   ```

2. After the preflight passes, fit four `brms` models concurrently. Each model
   uses its own four parallel chains, so this permits up to 16 sampling chains
   at once; reduce `workers` if memory pressure is high:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY = "fit",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_diagnostic_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

   This target retains each completed `brmsfit` as an independent file, so a
   later restart keeps completed branches. The completed initial audit has 92
   passing models and 12 models requiring a focused remedy.

3. Fit only those 12 remedy models with four workers:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_REMEDY = "fit",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_final_diagnostic_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_REMEDY",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

   Eleven models use 2,000 warmup and 4,000 retained iterations with
   `adapt_delta >= 0.999` and `max_treedepth >= 15`. Function model 11 already
   used that schedule initially, so its remedy advances to 3,000 warmup and
   6,000 retained iterations. Review the final diagnostics before accepting
   remedy fits or beginning LOO.

4. Once all remedy diagnostics pass, calculate integrated observation-level
   PSIS-LOO for the 104 accepted fits. Each branch uses one core and writes an
   update-specific compact reduction, leaving the historical strategy archive
   untouched:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY",
     "PROTEIN_UPDATE_2026_STRATEGY_REMEDY"
   ))
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_LOO = "run",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_integrated_loo_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_LOO",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

   The integrated audit identifies every observation with Pareto k above 0.7.
   Exact deletion refits are defined only after that audited case set is known.

5. The completed audit selected 16 observations for exact deletion refits:
   seven whole-body, seven hypertrophy, two function, and none for strength.
   Run the sampling-free exact-refit preflight first:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY",
     "PROTEIN_UPDATE_2026_STRATEGY_REMEDY",
     "PROTEIN_UPDATE_2026_STRATEGY_LOO",
     "PROTEIN_UPDATE_2026_STRATEGY_EXACT",
     "PROTEIN_MODEL_WORKERS"
   ))

   targets::tar_make(
     update_2026_strategy_exact_preflight_file,
     reporter = "balanced"
   )
   ```

   After reviewing the 16-row preflight, fit four exact refits concurrently.
   Each is an ordinary `brms::update()` fit with four chains, 2,000 warmup
   iterations, 6,000 retained iterations, `adapt_delta >= 0.999`, and
   `max_treedepth >= 15`:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_EXACT = "fit",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_exact_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_EXACT",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

   Completed branches are retained independently in the update-specific model
   archive, so restarting the command does not repeat finished refits.

6. After all 16 exact-refit diagnostics pass, convert the saved fits to compact
   pointwise corrections, patch the final LOO objects, and compute the four
   stacking vectors. This performs no additional sampling:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY",
     "PROTEIN_UPDATE_2026_STRATEGY_REMEDY",
     "PROTEIN_UPDATE_2026_STRATEGY_LOO",
     "PROTEIN_UPDATE_2026_STRATEGY_EXACT"
   ))
   Sys.setenv(PROTEIN_MODEL_WORKERS = "4")

   targets::tar_make_future(
     update_2026_strategy_stacking_checkpoint_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv("PROTEIN_MODEL_WORKERS")
   ```

7. Build and inspect the sampling-free reporting preflight, then generate the
   allocated posterior predictions with four workers:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY",
     "PROTEIN_UPDATE_2026_STRATEGY_REMEDY",
     "PROTEIN_UPDATE_2026_STRATEGY_LOO",
     "PROTEIN_UPDATE_2026_STRATEGY_EXACT",
     "PROTEIN_UPDATE_2026_STRATEGY_REPORTING",
     "PROTEIN_MODEL_WORKERS"
   ))

   targets::tar_make(
     update_2026_strategy_reporting_preflight_file,
     reporter = "balanced"
   )

   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_REPORTING = "run",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_reporting_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_REPORTING",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

8. Freeze the completed family as a verification bundle before planning any
   storage pruning:

   ```r
   targets::tar_make(
     update_2026_strategy_verification_files,
     reporter = "balanced"
   )
   ```

9. Generate the read-only lean-pruning plan. This inventories the 104 accepted
   model fits, 12 superseded initial fits, and 16 full exact-refit objects, but
   does not delete anything. Compact PSIS reductions, exact corrections,
   reporting objects, and the immutable verification bundle are retained:

   ```r
   targets::tar_make(
     update_2026_strategy_pruning_plan_files,
     reporter = "balanced"
   )
   ```

   Review `update_2026_strategy_pruning_plan` before separately authorising
   destructive lean pruning.

10. After explicit review and confirmation, execute pruning outside the target
    graph. The exact confirmation string is deliberately required, and a
    checksum receipt is written before any archive is deleted:

    ```r
    strategy_pruning <- targets::tar_read(
      update_2026_strategy_pruning_plan
    )

    strategy_pruning_receipt <-
      execute_update_2026_strategy_pruning(
        plan = strategy_pruning,
        confirmation = "DELETE ALL UPDATED STRATEGY FIT ARCHIVES"
      )
    ```

## Updated strategy-modified strict within-study family

This family is the next combined-data sensitivity analysis after the
arm-based strategy family. It fits the same 104 primary candidate structures
to strict matched protein-arm minus comparator-arm contrasts, with Addition
versus Substitution as a main effect and as an interaction with every protein
contrast basis. These are ordinary `brms` fits; handwritten Stan is not used.

The historical checkpoint supplies sampling schedules only. Six previously
remedied candidates begin with 2,000 warmup and 4,000 retained iterations,
`adapt_delta >= 0.999`, and `max_treedepth >= 15`. The other 98 candidates use
1,500 warmup and 2,000 retained iterations. Every model is refitted because
the combined dataset has changed.

1. Clear fitting switches and run the sampling-free 104-model input preflight:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN",
     "PROTEIN_MODEL_WORKERS"
   ))

   targets::tar_make(
     update_2026_strategy_within_input_preflight_file,
     reporter = "balanced"
   )
   ```

   Review all four outcome samples, the six inherited longer schedules, and
   confirm that every `valid` value is `TRUE` before sampling.

2. Fit four models concurrently. Each model runs four chains internally and
   each completed fit is retained as an independent update-specific file:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_WITHIN = "fit",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_within_diagnostic_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

   Stop at the diagnostic audit. Any updated failures receive the smallest
   data-driven remedy before LOO or model averaging is added.

3. The completed initial audit retains 101 passing models and selects only
   strength models 14, 21, and 22 for a first remedy. All three have zero
   divergences and zero maximum-treedepth transitions but narrowly miss the
   convergence or precision thresholds. Their remedy therefore keeps the
   original `brms` structures and increases sampling to 2,000 warmup and 4,000
   retained iterations per chain, with `adapt_delta >= 0.999` and
   `max_treedepth >= 15`.

   First run and inspect the sampling-free remedy preflight:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REMEDY",
     "PROTEIN_MODEL_WORKERS"
   ))

   targets::tar_make(
     update_2026_strategy_within_remedy_preflight_file,
     reporter = "balanced"
   )
   ```

   After confirming all three rows, run the three independent remedies in
   parallel and stop again at the final sampling audit:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REMEDY = "fit",
     PROTEIN_MODEL_WORKERS = "3"
   )

   targets::tar_make_future(
     update_2026_strategy_within_accepted_fit_manifest_file,
     workers = 3,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REMEDY",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

4. After the final audit accepts all 104 models, run the sampling-free LOO
   preflight. The accepted registry contains 101 initial fits and three remedy
   fits, and each LOO branch is limited to one core:

   ```r
   targets::tar_make(
     update_2026_strategy_within_loo_preflight_files,
     reporter = "balanced"
   )
   ```

   Once reviewed, calculate four integrated observation-level PSIS-LOO
   reductions concurrently. This opens the accepted fits but does not refit
   them:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_WITHIN_LOO = "run",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_within_integrated_loo_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_LOO",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

5. The completed integrated audit selects 65 high-Pareto-k observation cases:
   10 whole-body, 26 hypertrophy, none for strength, and 29 function. This is
   close to the 67 exact cases in the frozen historical strategy-within family
   and reflects the comparatively small strict within-study samples. Run the
   sampling-free exact-refit preflight before starting any deletions:

   ```r
   targets::tar_make(
     update_2026_strategy_within_exact_preflight_file,
     reporter = "balanced"
   )
   ```

   After confirming all 65 records, run four exact deletion refits at a time.
   Each is an ordinary `brms::update()` with 2,000 warmup and 6,000 retained
   iterations per chain:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT = "fit",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_within_exact_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

6. The first exact-refit audit accepts 64 cases. Hypertrophy model 22,
   observation 158 (case 30), has one divergence but otherwise passes all
   convergence and efficiency checks. Preserve the 64 accepted fits and run a
   sampling-free preflight for the single geometry remedy. The remedy changes
   only `adapt_delta` to 0.9995 and uses a fresh seed:

   ```r
   targets::tar_make(
     update_2026_strategy_within_exact_remedy_preflight_file,
     reporter = "balanced"
   )
   ```

   After review, run the isolated exact-refit remedy:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT_REMEDY = "fit",
     PROTEIN_MODEL_WORKERS = "1"
   )

   targets::tar_make(
     update_2026_strategy_within_exact_remedy_audit_file,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT_REMEDY",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

7. Once case 30 passes, combine it with the other 64 accepted exact fits and
   generate the compact exact corrections. This stage performs no sampling;
   four workers may read and reduce four saved fits concurrently. The final
   targets apply all 65 corrections and calculate the four stacking vectors:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT",
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_EXACT_REMEDY"
   ))
   Sys.setenv(PROTEIN_MODEL_WORKERS = "4")

   targets::tar_make_future(
     update_2026_strategy_within_stacking_checkpoint_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv("PROTEIN_MODEL_WORKERS")
   ```

8. After the corrected LOO audit is stacking-ready, run the sampling-free
   reporting preflight. It allocates exactly 4,000 posterior-mixture draws per
   outcome and retains only models allocated at least one draw:

   ```r
   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REPORTING",
     "PROTEIN_MODEL_WORKERS"
   ))

   targets::tar_make(
     update_2026_strategy_within_reporting_preflight_file,
     reporter = "balanced"
   )
   ```

   Once reviewed, generate posterior predictions for the selected models and
   write the four model-averaged reporting objects. Four workers open four
   accepted fits concurrently; this stage performs no sampling:

   ```r
   Sys.setenv(
     PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REPORTING = "run",
     PROTEIN_MODEL_WORKERS = "4"
   )

   targets::tar_make_future(
     update_2026_strategy_within_reporting_audit_files,
     workers = 4,
     reporter = "balanced"
   )

   Sys.unsetenv(c(
     "PROTEIN_UPDATE_2026_STRATEGY_WITHIN_REPORTING",
     "PROTEIN_MODEL_WORKERS"
   ))
   ```

9. Freeze the completed family in an immutable verification bundle. This
   rechecks the source hashes, all five 104-model registries, the 101+3
   accepted-fit lineage, 65 exact corrections, four stacking vectors, and all
   reporting files:

   ```r
   targets::tar_make(
     update_2026_strategy_within_verification_files,
     reporter = "balanced"
   )
   ```

10. Build and inspect the lean-pruning plan. This target is non-destructive:
    it only inventories full fit archives proposed for deletion and compact
    artifacts that must remain:

    ```r
    targets::tar_make(
      update_2026_strategy_within_pruning_plan_files,
      reporter = "balanced"
    )
    ```

    Do not delete anything until the verification bundle and every row of the
    pruning plan have been reviewed and explicitly confirmed.

11. After explicit confirmation, execute the guarded lean-pruning operation.
    The executor revalidates the immutable verification manifest, the complete
    173-file fit inventory, and all 212 retained compact artifacts before it
    deletes anything. It first writes a pending receipt containing a fresh
    SHA-256 checksum for every fit archive, then deletes only the audited files
    and promotes the pending receipt to `executed_lean.csv` after the retained
    artifacts pass a second integrity check:

    ```r
    strategy_within_pruning <- targets::tar_read(
      update_2026_strategy_within_pruning_plan
    )

    strategy_within_pruning_receipt <-
      execute_update_2026_strategy_within_pruning(
        plan = strategy_within_pruning,
        confirmation =
          "DELETE ALL UPDATED STRATEGY-WITHIN FIT ARCHIVES"
      )
    ```

    This removes the 104 accepted model-fit objects, three superseded initial
    model-fit objects, 65 accepted exact-refit objects, and one superseded
    exact-refit object. The frozen verification bundle, corrected LOO objects,
    exact corrections, stacking weights, predictions, and reporting outputs
    remain; recovering the full fit objects later requires refitting this
    family.
