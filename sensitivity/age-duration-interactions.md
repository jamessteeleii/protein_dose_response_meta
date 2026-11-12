# Primary-family protein interactions with age and duration

Prepared 2 November and completed 9 November 2026. This post-registration
extension produced 336 new interaction candidates, required recovery of 72
previously pruned absolute-primary fits, and completed 761 interaction fitting
jobs plus 73 recovery jobs after remedies and exact deletions. The recorded
model-job time was 300.55 hours. See the root
[reproducibility and fitting-time audit](../REPRODUCIBILITY.md) before any
rerun. Both muscle-measurement analyses remain prerequisites.

## Scope and design

Only the frozen 2026 **primary** families are expanded. Neither Mundlak,
within-study, achieved/assignment, EIV, residual-correlation nor intervention
strategy families enter these new averages.

Training-status interactions were initially deferred on 2 November 2026
because the reported classifications may poorly measure prior effective
training exposure; see [the author's methodological discussion](https://steele-research.com/essays/how-trained-are-trained-participants/).
They were subsequently completed as a separate exploratory workflow and are
not part of the age/duration model averages. This extension remains
protein-by-age and protein-by-duration only.
Training-status labels do not determine eligibility, so no additional
common-sample restrictions or primary-model refits are introduced for that
covariate. The standalone exploratory audit in
`_targets_training_status.R` is not part of this pipeline.

| Exposure | Original candidates per outcome | New candidates per outcome | Expanded candidates per outcome | New fits across four outcomes |
| --- | ---: | ---: | ---: | ---: |
| Relative intake | 26 | 36 | 62 | 144 |
| Absolute total intake | 18 | 24 | 42 | 96 |
| Absolute supplementation | 18 | 24 | 42 | 96 |

For every covariate-adjusted dose candidate, add age-only, duration-only and
joint protein interactions. Both additive moderators remain in every new
candidate. Intercept-only candidates remain in the original set but do not
generate interaction duplicates. The absolute families retain their original
exclusion of the fixed 1.6-g/kg/day piecewise candidates. No absolute-dose
breakpoint is invented.

Protein-by-moderator terms follow the inherited dose basis. Thin-plate smooths
use numeric-by varying coefficients, constrained to zero at the frozen centred
protein reference to avoid duplicating additive moderator effects. Existing
natural-spline remedies retain their fixed basis. RT-specific parents allow
the moderator interactions to differ by RT and include moderator-by-RT
lower-order terms. Joint candidates contain both protein interactions, without
an age-by-duration or protein-by-age-by-duration term. Random study slopes and
study factor smooths retain their original structure.

The original calibrated coefficient priors are frozen in
`sensitivity/age-duration-base-priors.csv`; their four source fit hashes are in
`sensitivity/age-duration-prior-provenance.csv`. This portable recipe avoids
depending on prior targets in the main cache or recalibrating priors with
different package availability. `freeze_age_duration_priors.R` documents its
extraction from accepted relative-primary intercept-only fits. Inherited
candidate-specific priors are reconstructed by the existing preparation
helpers. New fixed coefficients use Student-t(3, 0, 0.5) priors; added smooth
deviations inherit the parent's smooth-SD prior, or half-normal(0, 0.5) if none
was specified. Age and duration are standardised on the exact outcome/common
sample, using the existing transformations. New fits use four chains, 2,000
warm-up and 6,000 retained iterations per chain, adapt_delta at least 0.999,
and maximum tree depth at least 15. Prior predictive/empirical adequacy remains
a matter for pilot review; a valid design matrix does not establish it.

## Isolation and gates

- Script: `_targets_age_duration.R`; store: `_targets_age_duration/`.
- New archives: `model_archive/update_2026/age_duration/`.
- Preparation, diagnostics, LOO and reporting: `results/update_2026/age_duration/`.
- Separate public report: `supplement/age-duration-interactions.qmd`.
- Portable report interface: `manuscript/generated/age-duration-interactions/`.
- Default phase: `PROTEIN_INTERACTION_PHASE=none`. Every sampling, LOO,
  correction, combination and prediction helper also has an execution guard.
- Every expensive phase requires both completed, checksummed
  `measurement_prereg_methods_update_2026` and
  `measurement_direct_only_update_2026` verification bundles. It never runs
  those analyses, edits their targets, or requests their main pipeline.
- Existing manuscript figures and original primary averages are not replaced.
- The frozen likelihood rows, their order and their numerical values must
  agree across candidates. In particular, total-intake muscle models use the
  accepted **421-row common sample**, rather than mixing 425-row and 421-row
  fits. Supplementation likewise retains its final body-mass-complete sample.
- Older saved model frames contain only formula variables, and older LOO
  manifests omit time points. The source common sample must match the final
  LOO order exactly. Available moderator fits additionally verify their actual
  predictor values, including standardised duration. No partial-key join can
  interchange identical effects measured at different follow-ups.
- New fit archives preserve source identity columns in the model data.
  Predictor fingerprints use 13 significant digits to tolerate historical
  CSV round-trip precision; raw data values are never rounded or modified.
- All candidates must pass the existing geometry, convergence, ESS and MCSE
  gates before LOO. PPC warnings and observation-level diagnostics are saved.
  Models are not dropped to bypass failures. Remedy decisions go in the
  `age-duration-remedies.csv`. The remedy hook supports sampler control
  escalation and an optional longer retained sampling schedule. Other changes
  require targeted review.
- Conditional PSIS is audited, then the existing Gaussian observation-level
  integration is used with one internal core per branch. Its adapter uses the
  original parent index solely to reconstruct random-effect structure; saved
  new candidates keep distinct IDs. All new Pareto-k > 0.7 cases require exact
  refits. The dose-family batch scores each passing refit in the same worker
  and saves its diagnostics and correction together under `loo/exact_checkpoints/`.
  Failed diagnostics are retained for remedy review. The subsequent `correction`
  checkpoint validates and merges the saved corrections. An interrupted or
  errored case before its checkpoint is written needs to be refitted. The
  individual-stage workflow retains full exact fits for separate scoring.
- `combine` requires complete, corrected LOO on matching rows, checks actual
  saved fit identities and predictors, and estimates fresh stacking weights
  within each exposure/outcome. Original weights are not appended or reused
  as weights for the expanded ensemble. Weight sums are not posterior
  probabilities that a moderator hypothesis is true.

## Current storage and recovery requirement

On 2 November, the relative-primary 104 accepted fits and supplement-primary
72 accepted fits were present. The **72 accepted total-primary fits were
absent**, with an existing lean-pruning receipt dated 13 October. Compact LOO,
diagnostics, specifications, data and reported results remain intact.

Preparation can audit all three datasets without those absent posteriors.
Interaction fitting/combination is blocked until the required primary fits
are available. Restore a verified external backup if one exists. Otherwise,
the isolated recovery workflow refits the 72 total-primary candidates on their
frozen final common samples, without touching the old outputs. It inherits
the accepted diagnostic controls and preserves the recorded 12,000 retained
iterations for muscle candidate 11. The one sampling-length exception is
documented in `age-duration-total-recovery.csv`.

Recovery receives fresh diagnostics, integrated LOO and any newly needed exact
corrections. It does **not** silently associate new posterior draws with old
LOO reductions. Only `recover_accept` writes the replacement mapping consumed
by this extension. Historical checkpoints and their now-missing archive paths
remain unchanged. New recovery failures require reviewed remedies before
advancing; the original numerical findings are not guaranteed to be identical
to a new Monte Carlo run.

There are 336 new interaction fits, plus 72 recovery fits if a backup cannot
be restored. Preparation estimates archive sizes from corresponding parent
fits and records a 1.5-times planning allowance. This is not a space guarantee:
remedies, exact refits, model compilation and reporting need extra room.

### Safe completed-family pruning

The read-only audit verified all retained bundle and compact-file SHA-256
checksums for **13 non-primary absolute families**. Their existing planned
deletion files occupied 566.9 GiB. The largest choices were:

| Family | Reclaimable GiB |
| --- | ---: |
| Supplement residual correlation | 119.3 |
| Total-intake residual correlation | 116.4 |
| Supplement strategy | 73.4 |
| Total-intake strategy | 68.6 |

The first two provide 235.7 GiB; all four provide 377.7 GiB. Counts describe
the audit time, not a future guarantee. The complete table is generated at
`results/update_2026/age_duration/storage/pruning_audit.csv`.

In an R 4.4.3 console at the project root, this rechecks the audit and **does
not delete anything**:

```r
source("sensitivity/age-duration-storage.R")
iad_storage_audit()
```

If the author chooses to reclaim the two largest eligible families, these
are the explicit execution commands (not run by the preparing agent):

```r
iad_prune_completed_family("supplement_rescor", execute = TRUE)
iad_prune_completed_family("total_rescor", execute = TRUE)
```

For the next two eligible families, if more space is needed:

```r
iad_prune_completed_family("supplement_strategy", execute = TRUE)
iad_prune_completed_family("total_strategy", execute = TRUE)
```

Other eligible keys are `total_assignment`, `total_eiv`, `total_mundlak`,
`total_strategy_within`, `total_within`, `supplement_mundlak`,
`supplement_within`, `supplement_strategy_within`, and `total_achieved`.
Use `iad_prune_completed_family("<key>")` first for an individual dry run.

The helper accepts only those 13 named families, checks current sizes and
every retained hash, rejects overlap with protected primary fits/LOO and
compact outputs, checks resolved paths stay inside the named family archive,
and deletes only listed leaf files. It records hashes and a deletion receipt.
It never deletes directories, target-store contents, raw data or manuscript
interfaces. Pruning preserves existing reporting but removes the ability to
make new predictions from those sensitivity posteriors without refitting.
The non-primary absolute-family helper excludes the primary and recovery
archives, the interaction fits, and the muscle-measurement analyses.

### Lean plan for the completed age and duration fits

The completed interaction models have a separate dry-run plan. It lists only
the full fits recorded in the initial and accepted diagnostic manifests for
the three interaction families. The 11 November audit found 337 files using
286.04 GiB, comprising 336 accepted fits and the superseded initial fit for
relative-intake physical-function model 32. The original primary fits and the
recovered absolute total-intake primary fits are outside the deletion list.

From the project root in R 4.4.3, regenerate the review files with:

```r
source("sensitivity/age-duration-pruning.R")
iad_plan_interaction_pruning()
```

The script writes `summary.csv`, `deletion_file_manifest.csv`, and
`retained_file_manifest.csv` under
`results/update_2026/age_duration/storage/interaction_pruning/`. It checks
the fit paths and sizes, records the fit hashes from the diagnostics, and
checks the twelve reporting manifests and their artifacts against their saved
hashes. It verifies the final LOO reductions against the accepted fits and
exact-case records, then records checksums for the compact LOO files, the
exact-refit checkpoints, the report files, and the static supplement. The
training-status workflow uses the retained reports and supplementary interface
for its age/duration completion gate.

This is a plan for review. The script only writes the three plan files. Before
executing any deletion, compare each candidate file with its recorded SHA-256
and recheck the retained manifest. The current predictions, contrasts, plots,
weights, and diagnostics can be used from the retained reports. Generating
new conditional profiles or posterior summaries after pruning requires
restoring the full fits or refitting. The standalone supplement can still be
rendered from the static outputs.

## Reviewed-stage R handoff

### Unattended initial-fit batch (author request, 3 November 2026)

The author requested one launch for all three exposure definitions, replacing
manual pilot-to-fit handoffs with independent automated pilot gates. Use
`sensitivity/run_age_duration_batch.R` for this route; do not run a legacy
interaction target or another copy of the batch against this store concurrently.
The original scientific specifications, priors, seeds, common samples and
four-chain schedules are unchanged.

| Dose family | Pilot fits | Remaining fits | Total initial fits |
| --- | ---: | ---: | ---: |
| Relative total intake | 8 | 136 | 144 |
| Absolute total intake | 8 | 88 | 96 |
| Absolute supplementation | 8 | 88 | 96 |
| Total | 24 | 312 | 336 |

There is one four-model scheduler (at most 16 chains), not three independent
four-worker jobs. Each family's eight pilots span all four outcomes and must
pass the existing geometry, convergence, ESS and MCSE gates on the exact
expected data before that family's remaining fits are released. A failure
blocks that family's descendants while independent families may finish;
errors are reported at the end, with successful archives retained. PPC warnings
remain recorded for later review and are not silently turned into exclusions.
No remedies, new-model LOO, stacking or final reporting are automatically run
after initial fitting. Full diagnostic review follows separately by dose family.

The launcher first retries/completes the recovered-primary reference reporting
checkpoint, then builds all 336 non-sampling Stan inputs and profile/provenance
tables. Recovery reporting had failed because its model-average object omitted
the settings list needed by the contrast helper; this is corrected and covered
by a synthetic regression test. The recovered posteriors and LOO are unchanged.

At preparation it checks free space on the project drive against **1.5 times
the parent-sized estimate for still-missing fits plus a 20 GiB reserve**.
Already saved fits are credited when resuming. This is a planning gate, not a
guarantee about archive size. The 3 November audit estimated 318.2 GiB for all
336 fits, plus reserve (about 338 GiB), versus 251.3 GiB free. The reverified
completed `supplement_strategy` and `total_strategy` leaf-file plans reclaim
73.4 and 68.6 GiB (142.0 GiB together); their retained public bundles and compact
outputs passed SHA-256 checks. No primary, recovery, measurement, training-status
or new interaction file is eligible. Pruning is never automatic. If choosing
those two plans, they can be executed by the approved initial-batch launcher
below, or with the existing individual commands:

```r
source("sensitivity/age-duration-storage.R")
iad_prune_completed_family("supplement_strategy", execute = TRUE)
iad_prune_completed_family("total_strategy", execute = TRUE)
```

The author approved preparing the two strategy-family pruning plans to run
before the interaction fits. In R 4.4.3 at the project root, one explicit call
rechecks both retained-artifact plans, executes only their named leaf-file
deletions, retries/completes recovery reporting, checks preparation/storage,
and queues all 336 initial fits with separate dose-family pilot gates:

```r
source("sensitivity/run_age_duration_batch.R")
run_age_duration_initial_batch(workers = 4L, prune_completed_strategies = TRUE)
```

This includes all age-only, duration-only and joint age/duration candidates;
it does not include training-status fits. On success inspect the three files
`results/update_2026/age_duration/{relative,total,supplement}/checkpoints/initial_fits.csv`.
Pilot audits and pass records are under each family's `diagnostics/pilot/`.
If interrupted or errored, review the failure, then use the same launcher;
valid saved fits are reused and no cache is erased.

Pruning defaults to FALSE in the wrapper; the explicit TRUE above is the
author-selected deletion action. Both strategies are verified before deletion,
and each execution rechecks its plan and records checksums/deletion receipts.
Already completed pruning is a no-op on restart; a pending partial receipt
requires review rather than silent continuation. If space was reclaimed
separately, use `run_age_duration_initial_batch(workers = 4L)` or
`run_age_duration_checkpoint("initial_fit", workers = 4L)`.

For subsequent checkpoints choose **one dose family at a time**, for example:

```r
run_age_duration_checkpoint("diagnostics", families = "relative", workers = 4L)
```

The same entry point supports `remedy`, `loo`, `exact`, `exact_remedy`, `correction`, `combine`
and `report`, using `relative`, `total` or `supplement`. Follow checkpoint review
before each new stage. These map to `iad_<family>_initial_audit_files`,
`accepted_audit_file`, `integrated_audit_files`, `exact_audit_file`,
`final_loo_files`, `combination_audit_files`, and `report_checkpoint_file`.
All checks, none-needed markers and audit paths are isolated by family; no
weights are pooled across exposure definitions. The final `supplement` stage
requires all twelve exposure/outcome reports before updating the shared static
supplement interface. Do not request it until all three reporting checkpoints
have been reviewed.

The exact-refit run completed on 6 November. All 422 cases finished, with
420 passing sampling diagnostics. Relative-intake physical-function models
44 (protein × duration) and 56 (protein × age and duration), both omitting
observation 193, had maximum Rhat 1.010826 and 1.015271 respectively.
Both passed geometry and ESS gates. The reviewed ledger
`sensitivity/age-duration-exact-remedies.csv` records 12,000 retained
iterations per chain for these two cases, retaining four chains, 2,000 warmup
iterations, seed 1988, formula, data, priors and sampler controls. The shared
exact-fit helper permits this longer schedule only through the explicit
`iter_sampling = 12000L` argument; its default remains 6,000.

Run the retry checkpoint across the three dose families with:

```r
source("sensitivity/run_age_duration_batch.R")
run_age_duration_checkpoint("exact_remedy",
  families = c("relative", "total", "supplement"), workers = 4L)
```

Only the two failed cases are sampled. The other 420 retain their existing
checkpoint paths. Initial checkpoints remain untouched, and compact retries
are written under `loo/exact_checkpoints/remedy_1/`. Each retry records the
initial checkpoint checksum, reviewed settings and diagnostics. A completed
matching retry is reused on restart. As with the initial compact run, an
interruption before a case saves its checkpoint requires rerunning that case.
Review the three `loo/exact_review_diagnostics.csv` files before advancing to
`correction`. Failed retries save their diagnostic evidence and remain blocked
from scoring. No retry or other pipeline stage was executed by the agent.

The 9 November combination attempt exposed a legacy-frame compatibility
issue. Some saved primary spline fits retain only formula variables, including
`study_factor`, and omit the redundant `study_no` column. The observation
identity helper now recognises the study factor's character labels when the
study ID column is absent. It continues to check study, arm and outcome IDs,
row order, effects and variances; the existing predictor and formula checks
remain in place. Synthetic checks reject altered IDs, reordered rows and
changed likelihood inputs, including factors whose level codes differ from
their study labels. No fit or LOO archive is modified by this adapter. The
combination error now names the family, outcome, model and failed check.
Resume with the same `combine` checkpoint command after reviewing this fix.

The subsequent reporting attempt exposed a tibble data-mask collision in
`iad_profile_support()`. Its output column named `profile` shadowed the input
record before later columns accessed the record's age and duration. The helper
now uses the distinct argument name `profile_row`. The reporting regression
script `sensitivity/test_age_duration_reporting.R` checks support counts and
empty neighbourhoods, then exercises the complete reporting formatter with
synthetic predictions in a temporary project. It checks paired age/duration
contrasts against known values and verifies output checksums. Its optional
`--review-current` mode checks the grids and support tables for all 72 saved
profiles without posterior prediction or pipeline execution. Rerun only the
`report` checkpoint after this correction; fitted models and stacking weights
remain saved.

Synthetic batch tests cover missing/duplicated/failed pilots, data fingerprints,
independence of other dose families, phase boundaries, storage/resume accounting
and the recovered reporting settings. No batch was run by the preparing agent.

### Original individual-stage handoff

The 4 November initial diagnostic review passed 335 of 336 candidates.
Relative-intake physical function model 32 (protein × age, parent 14) had
R-hat 1.01126 for the arm intercept SD against the 1.01 gate. Geometry, ESS
and MCSE passed. The ledger records a `remedy_1` retry with 12,000 retained
iterations per chain, up from 6,000, using the original warmup, seed, model,
data, priors and sampler controls. Both absolute families passed every
sampling gate and need only the no-remedy marker at this checkpoint.

The author confirmed the execution order on 2 November: **complete recovery
of the missing absolute total-intake primary family first**. This includes
its diagnostics, any required remedies, corrected LOO, primary-only stacking
and reference-profile reporting checks. Only after reviewing that completed
family should interaction preparation/fitting begin. Expanded averaging is
a subsequent, separate supplementary step within each dose type; exposure
scales are not pooled into one ensemble.

Use the **R 4.4.3 console and project renv environment** from
the repository root. Run only the stage currently supplied
after checkpoint review. No `callr_function=NULL` or in-session pipeline runs.
Do not run all stages together. Default concurrency is four models, each with
four chains; lower `PROTEIN_MODEL_WORKERS` if other runs are still active.

The next command prepares **only the 72 total-primary recovery models**,
without sampling or preparing the 336 interaction candidates:

```r
stopifnot(getRversion() == package_version("4.4.3"))
Sys.setenv(PROTEIN_MODEL_WORKERS = "4", PROTEIN_INTERACTION_PHASE = "none")
targets::tar_make_future(
  names = tidyselect::all_of("iad_recovery_preparation_files"),
  script = "_targets_age_duration.R", store = "_targets_age_duration",
  workers = 4L, reporter = "balanced"
)
```

Return the printed contents of:

```r
readr::read_csv("results/update_2026/age_duration/recovery/preflight.csv", show_col_types = FALSE) |>
  dplyr::count(outcome, valid)
```

Or reply **“recovery preparation complete”** so the agent can inspect the
preflight and specifications files directly.
Do not begin a recovery pilot until both measurement bundles are complete,
storage is sufficient, and the preparation checkpoint has been reviewed.

For later stages, use this same command form, changing only the phase and
target from the table below:

```r
Sys.setenv(PROTEIN_INTERACTION_PHASE = "recover_pilot")
targets::tar_make_future(
  names = tidyselect::all_of("iad_recovery_pilot_audit_files"),
  script = "_targets_age_duration.R", store = "_targets_age_duration",
  workers = 4L, reporter = "balanced"
)
```

| Stage | Phase | Target | Checkpoint relative to `results/update_2026/age_duration/` |
| --- | --- | --- | --- |
| Recovery pilot (8 fits, reused) | `recover_pilot` | `iad_recovery_pilot_audit_files` | `diagnostics/recovery_pilot/summary.csv` |
| Remaining recovery fits (64) | `recover_fit` | `iad_recovery_fit_file` | Then diagnostics stage |
| Recovery diagnostics | `none` | `iad_recovery_audit_files` | `diagnostics/recovery_initial/summary.csv` |
| Reviewed recovery remedies | `recover_remedy` | `iad_recovery_remedy_audit_files` | `diagnostics/recovery_remedy/summary.csv` |
| Recovery integrated LOO | `recover_loo` | `iad_recovery_loo_audit_files` | `recovery/loo_summary.csv`, `recovery/exact_cases.csv` |
| Recovery exact fits | `recover_exact` | `iad_recovery_exact_audit_file` | `recovery/exact_diagnostics.csv` |
| Recovery scoring/final LOO | `recover_correction` | `iad_recovery_final_loo_files` | Compact corrected files; inspect exact diagnostics |
| Accept recovered primary family | `recover_accept` | `iad_recovery_accepted_file` | `recovery/accepted_manifest.csv` |
| Recovered primary-only weights and reference-profile reporting | `recover_report` | `iad_recovery_reporting_audit_files` | `recovery/reporting/audit.csv`, four outcome directories |
| Interaction pilot (24 fits, reused) | `pilot` | `iad_pilot_audit_files` | `diagnostics/pilot/summary.csv` |
| Remaining interaction fits (312) | `fit` | `iad_fit_file` | Then diagnostics stage |
| Initial interaction diagnostics | `none` | `iad_initial_audit_files` | `diagnostics/initial/summary.csv` |
| Reviewed remedies or no-remedy marker | `remedy` | `iad_accepted_audit_file` | `diagnostics/accepted.csv` |
| Interaction integrated LOO | `loo` | `iad_integrated_audit_files` | `loo/summary.csv`, `loo/exact_cases.csv` |
| Interaction exact fits | `exact` | `iad_exact_audit_file` | `loo/exact_diagnostics.csv` |
| Exact scoring/final LOO | `correction` | `iad_final_loo_files` | Compact corrected files; inspect exact diagnostics |
| Gated expanded averages | `combine` | `iad_combination_audit_files` | `combination/stacking_weights.csv`, `combination/summary.csv` |
| Supplementary reporting | `report` | `iad_supplement_file` | Twelve reporting manifests, static interface and bundles |

For shell users, `sensitivity/run_age_duration.R <stage>` provides the same
isolated clean-process workflow. Use `recover_prepare` for the current
recovery-only handoff. Its `prepare` stage belongs to the later interaction
work. A failure is returned for review; do not manually edit `_targets*`
stores or erase successful branches.

Recovery diagnostic review, 2 November 2026: all 72 fits were saved and 71
passed. Strength candidate 24 had one divergent transition, with convergence,
ESS, E-BFMI and tree-depth checks passing. The reviewed ledger
`age-duration-recovery-remedies.csv` schedules only that candidate at
adapt_delta 0.999 (previously 0.995), retaining tree depth 15, four chains,
2,000 warmup and 6,000 retained iterations, the same seed, formula, priors,
and data. The new archive is under `recovery/remedy_1/strength/`; the original
recovery file and its initial diagnostics are preserved. Remedy diagnostics
must pass before LOO; accepted selection retains the other 71 fits and records
the selected stage. A failed remedy requires further checkpoint review.

Recovery acceptance review, 3 November: all 72 accepted fits and final LOO
reductions passed identity, observation-order and checksum checks. Strength
candidate 24 uses `remedy_1`; the other 71 retain `recovery_1`. No integrated
Pareto-k exceeded 0.7 (maximum 0.69444), so the exact stage recorded no cases.

The `recover_report` checkpoint now reconstructs four primary-only ensembles
from those new LOO reductions, with 4,000 paired mixture draws on the original
1 g/day grid, age 25, duration 12 weeks, and mean body mass. It saves fresh
weights, conditional means, contrasts from the frozen reference dose, full
pairwise contrasts/support, and local-equivalence summaries under
`recovery/reporting/`. Historical weight and curve comparisons are descriptive
checks of the new Monte Carlo run, not automatic equivalence tests; review
them before starting interaction preparation. The original reports and
manuscript interface are unchanged. Interaction preparation, pilots and
combination require the completed recovery reporting checkpoint, tied to the
accepted manifest. No fitting, LOO or rendering is repeated in this stage.

Temporary stacking inputs standardise draw-count metadata, following the
original absolute-primary helper, because accepted muscle candidate 11 has
twice the retained draws. Only pointwise ELPD enters stacking; saved LOO
objects and posterior draws are not changed. The expanded ensembles use the
same adapter to handle mixed original and interaction sampling lengths.

## Planned reporting

The six profiles are fixed from data coverage before examining new posterior
results: original age 25/week 12; outcome-specific median age/duration;
age quartiles at median duration; duration quartiles at median age. Age
quartiles count each distinct arm once. Duration quartiles count each distinct
study/follow-up once. They are not participant-weighted population summaries.
Repeated quartiles give coincident profiles and zero between-profile contrasts
by construction; they provide no separate test of moderation at those values.

Reporting preserves the relative 0.01-g/kg/day and absolute 1-g/day dose grids,
the frozen reference dose, the primary 0.1-SD SESOI and 95% probability rule.
All profiles use coherent component/posterior draws. Paired differences in
protein contrasts compare age Q3 versus Q1 and duration Q3 versus Q1 for the
same current/proposed intake and RT condition. The same-threshold probabilities
for these differences are descriptive; the primary SESOI was not separately
validated for effect modification.

Support tables show the number of studies/arms within 10 years and four weeks
of each profile, and that neighbourhood's observed dose range. Inferential
flags require at least three such studies and both doses within that range.
This transparent screen is not a convex-hull check or proof of causal joint
support. Plots show the entire prediction grid, including poorly supported
regions, which must be read alongside the support tables. Predictions concern
study/arm mean age, not individual age effects.

For absolute-exposure moderator models, body mass stays at the outcome mean.
Each scope includes a primary-only reconstruction on the **same frozen common
sample and same profile**, so changes from the expanded average are not
confounded with changing the profile. Original 25/12 figures remain intact.

The separate Quarto supplement reads only committed static derivatives and
the `_results.qmd` interface. It never reads a target store or fits. The initial
interface clearly says results are pending. Once the report gate has passed,
review and commit the generated interface; render locally with:

```r
quarto::quarto_render("supplement/age-duration-interactions.qmd")
```

If the Quarto R package is unavailable, use the documented installed CLI:
`quarto render supplement/age-duration-interactions.qmd`. No new R dependency
is required for that standalone render.

## Verification record

- All 336 interaction Stan inputs built without compilation or sampling.
- Default-off and missing-stage guards, common-sample recovery, predictor
  transformations, covariance-preserving contrasts, and corrected-LOO
  observation/diagnostic checks have non-sampling coverage.
- The separate target graph validates. No main `_targets.R`, existing model
  function, dependency lockfile, manuscript result or main target cache was
  changed by this work.
- The pruning audit verified all 13 eligible retained bundles/compact archives;
  it deleted no files.
- Posterior diagnostics, actual integrated/exact LOO compatibility of the new
  fitted smooths, runtime, disk use, new predictions and the final rendered
  results remain to be checked after the author's staged runs.
