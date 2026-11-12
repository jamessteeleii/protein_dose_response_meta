# Reported training status in RT studies: staged supplementary workflow

Completed 11 November 2026. The final workflow contains 248 candidate models
and required 359 fitting/refitting jobs after two candidate remedies, 102
initial exact deletions, conventional exact retries, and four marginal exact
fits. Recorded model-job time, including the standalone marginal pilot, was
101.42 hours. The complete fitting and storage context is in the root
[reproducibility and fitting-time audit](../REPRODUCIBILITY.md). The stage
commands below are retained for provenance and controlled reproduction; do
not execute them as a batch.

## Historical preparation and pruning record

On 9 November 2026, all twelve age/duration reporting manifests and both
completed muscle-measurement bundles were checksum-verified. Training
preparation was then the next user-run checkpoint. The training adapter
supplied its own family pilot gate to match the shared fitter updated during
the age/duration runs. The fresh storage estimate showed 288.2 GiB free after
the pruning below.

The 9 November lean-pruning audit reverified every retained compact artifact and
public bundle for the 13 eligible non-primary absolute families. At that audit,
nine had planned files present, totalling 189.2188 GiB: `total_assignment`, `total_eiv`,
`total_mundlak`, `total_strategy_within`, `total_within`, `supplement_mundlak`,
`supplement_within`, `supplement_strategy_within`, and `total_achieved`. The two
residual-correlation and two strategy families were already pruned. All four
older relative-sensitivity deletion manifests also have zero files remaining.
No deletion was performed by the audit. Existing `iad_prune_completed_family()`
commands in `sensitivity/age-duration-storage.R` recheck each plan immediately
before an explicit user-requested execution. These nine plans preserve all
training source inputs, primary/recovery/age-duration archives and measurement
bundles. They preserve current sensitivity reporting but remove the listed full
posteriors; new predictions from those sensitivities would require refitting or
restoring an external backup. Review the fresh free space after pruning before
moving to training preparation.

Later on 9 November, the author authorised execution of those nine prepared
plans after confirming that the existing sensitivity prediction profiles would
be retained. All nine completed: 1,185 planned leaf files were removed,
reclaiming 189.2188 GiB and leaving 288.21 GiB free on the project drive.
Each family's `pruning/executed_for_age_duration.csv` records the removed
files, their SHA-256 checksums and deletion times; no pending receipts remain.
The refreshed storage audit has zero planned files remaining across all 13
eligible families, and all retained bundle and compact-archive checksums pass.
Post-pruning checks also verified the three primary reporting bundles, all
twelve age/duration reports, both measurement bundles and all 21 frozen
training source inputs. The 104 relative-primary fits, 72 supplement-primary
fits and 72 recovered total-primary fits remained present. At that checkpoint,
the training pipeline had not yet been run; it was completed subsequently as
recorded at the top of this document.

Prepared 2 November 2026. Execute only after the age/duration analyses finish.
This is an exploratory post-registration extension of the three frozen primary
families. It does not expand the other sensitivity families or cross training
status with the new age/duration interactions.

## Frozen design

Use each exposure/outcome's frozen final primary common sample, then retain RT
arms (`resistance_exercise_code == 0`, checked against the yes/no label) with
reported trained or untrained status. Active and unreported records are excluded;
no source labels are reclassified. All candidates within a scope use exactly the
same rows, order, outcomes, variances and exposure definition. Exposure-specific
samples are retained; the three exposure scales are not averaged together.

Replace the original RT main effect with `training_status_code` (trained=0,
untrained=1). Replace protein-by-RT terms and RT-specific factor smooths with
protein-by-status terms and status-specific factor smooths. This gives the exact
analogue of the original family: common versus status-specific dose functions,
with/without the original additive covariates, common versus study-varying slopes
or smooths, plus no-dose models. Every candidate includes the status main effect.
The actual fitted formulas, including inherited diagnostic remedies, are audited.

Relative total intake has 26 candidates/outcome (linear, piecewise at 1.6 g/kg/day,
smooth and no-dose forms); absolute total and supplementation each have 18
(linear, smooth and no-dose forms). Total: **248 fits**, including **24 reusable
pilots** (indices 6 and 26 in each of twelve scopes), then **224 remaining fits**.
Original study, arm and outcome dependence structures are retained. Physical
function remains included for completeness, with descriptive interpretation.

| Outcome | Relative trained/untrained studies | Total | Supplement |
| --- | ---: | ---: | ---: |
| Whole body | 22/45 | 21/45 | 18/36 |
| Muscle size/mass | 7/33 | 7/33 | 8/26 |
| Strength | 21/49 | 20/49 | 17/38 |
| Physical function | 2/17 | 2/16 | 1/13 |

Status has no within-study variation in these frozen data. Comparisons are
between-study moderation associations; they cannot identify individual training
history effects. Range coverage and repeated outcome rows do not increase the
number of independent trained studies. Small interactions cannot establish that
actual accumulated training exposure is irrelevant.

## Priors, fitting and safeguards

Use the existing frozen coefficient-prior recipe in
`age-duration-base-priors.csv`, derived by `freeze_age_duration_priors.R`, plus
applicable parent variance, protein and covariate priors. Remove priors on old RT
coefficients rather than interpreting them as priors on training status. New
status main effects and status-by-dose fixed coefficients use Student-t(3,0,0.5).
Smooth SDs inherit applicable parent priors or use half-normal(0,0.5) when none
was specified. Fixed effects are checked for full rank by `make_standata`.
Age/duration and body-mass covariates are standardised within each restricted
outcome/exposure dataset. Dose centring and the original exposure reference are
retained. These choices are recorded before posterior results are available.

Four chains, 2,000 warmup and 6,000 retained draws per chain; adapt_delta at least
0.999 and maximum tree depth at least 15. Reuse accepted formula remedies as the
starting design, then evaluate fresh diagnostics. No candidate is silently dropped.
Failed diagnostics require a reviewed remedy. `training-status-remedies.csv`
supports documented sampler-control escalation and reviewed retained-sampling
length increases through its `iter_sampling` column. Structural or prior
remedies require a targeted implementation after inspection. Keep successful
fits and exact refits. Do not erase target stores to address failures.

### 10 November initial diagnostic review

All 248 initial fit archives are present (150.6 GiB); diagnostics pass for 246.
There are no tree-depth hits. Two candidates require the recorded `remedy_1`:

- Total intake / whole body / 17: one divergent transition, with convergence
  passing. Raise adapt_delta from 0.999 to 0.9995; retain 6,000 draws per chain.
- Supplementation / strength / 14: geometry passes, but minimum tail ESS is
  226.5 and maximum MCSE/SD is 0.1307; 58 parameters fail the parameter screen.
  Retain 12,000 draws per chain with the existing adapt_delta 0.999.

Both retain four chains, 2,000 warmup, maximum depth 15, identical formulas,
priors and restricted samples. The shared fitter already honours the reviewed
sampling length. New fits use separate `remedy_1` paths; preserve all initial
fits and the 246 passing candidates. Fresh remedy diagnostics must pass before
the 248-row accepted manifest releases LOO. These are new full runs, not chains
appended to the previous posterior. No fitting was executed by the preparing
agent. Free space at review was 135.4 GiB.

The broad predictive-check screen flags 244 candidates. Across candidates,
90% predictive interval coverage ranges from 0.897 to 1.000; relative-intake
physical function has the largest standardised residual (7.54). These checks
remain reported and reviewed alongside observation-level LOO; they are not
sampling failures and do not justify dropping candidates or observations.

Every expensive stage checks all twelve completed age/duration reporting manifests
and their file checksums, the final age/duration supplement interface, and both
completed muscle-measurement verification bundles. This gate reads completion
artifacts; it never dispatches another analysis. Preparation is non-sampling and
can run independently, but the agreed schedule is to wait for age/duration.

The graph `_targets_training_status.R` uses store `_targets_training_status/`.
Archives live under `model_archive/update_2026/training_status/`; checkpoints
under `results/update_2026/training_status/`. It reuses the diagnostic and exact-LOO
helpers through adapters sourced only in this separate clean process. Do not
source `R/training_status_workflow.R` into the other analysis sessions. No original
primary posterior is used as an RT-only combination candidate; all 248 are refit.
Missing historical total-intake posteriors therefore do not block this workflow.
Frozen compact LOO/source checksums are still required for source-sample auditing.

Storage estimates use corresponding parent archive sizes and a 1.5 multiplier;
they are planning estimates, not bounds. Remedies, exact refits, compilation and
reporting require additional space. Protect all new archives until reporting and
the verification review are complete. This workflow performs no pruning.

## One reviewed checkpoint at a time

Use R 4.4.3 with the restored project environment. After age/duration completion,
the first handoff is preparation only:

```r
stopifnot(getRversion() == package_version("4.4.3"))
Sys.setenv(PROTEIN_MODEL_WORKERS = "4", PROTEIN_TRAINING_PHASE = "none")
targets::tar_make_future(
  names = tidyselect::all_of("tst_preparation_files"),
  script = "_targets_training_status.R",
  store = "_targets_training_status",
  workers = 4L, reporter = "balanced"
)
```

Reply **training preparation complete** for inspection of `preparation/preflight.csv`,
`counts.csv`, `status_counts.csv`, `specifications.csv`, and `storage_estimate.csv`.
All checkpoint paths below are relative to `results/update_2026/training_status/`.
Do not run later stages until the preceding checkpoint has been reviewed.

| Stage | Phase | Target | Checkpoint |
| --- | --- | --- | --- |
| Prepare | none | tst_preparation_files | preparation/preflight.csv |
| Pilot, 24 fits | pilot | tst_pilot_audit_files | diagnostics/pilot/summary.csv |
| Remaining 224 fits | fit | tst_fit_file | Then diagnostics |
| Diagnostics | none | tst_initial_audit_files | diagnostics/initial/summary.csv |
| Reviewed remedies/no-remedy marker | remedy | tst_accepted_audit_file | diagnostics/accepted.csv |
| Integrated LOO | loo | tst_loo_audit_files | loo/summary.csv, loo/exact_cases.csv |
| Exact case-deletion fits | exact | tst_exact_audit_file | loo/exact_diagnostics.csv |
| Reviewed exact remedies | exact_remedy | tst_exact_review_audit_file | loo/exact_review_diagnostics.csv |
| Marginal exact remedies | exact_marginal | tst_marginal_audit_file | loo/marginal_review_diagnostics.csv |
| Exact scoring/final LOO | correction | tst_final_loo_files | Corrected archive files |
| Fresh stacking | combine | tst_combination_files | combination/stacking_weights.csv |
| Supplement | report | tst_supplement_file | Twelve reporting/file_manifest.csv files |

For subsequent console runs use the same call, replacing only the phase and
target. The shell runner `sensitivity/run_training_status.R <stage>` provides
the same clean-process `tar_make_future()` handoff. Four workers means up to
16 chains; lower the count if other computations are running.

LOO uses the existing Gaussian observation-level integration of multilevel terms.
All integrated Pareto-k >0.7 cases receive exact refits. The 10 November review
selected compact exact checkpoints, as used in the age/duration workflow: the
refit exists in memory, its diagnostics are saved, and passing cases are scored
before the compact checkpoint is written. Full exact-refit posteriors are not
archived. The correction stage validates/copies the passing scores and assembles
final LOO separately. Missing, duplicate, failed or mismatched corrections
block stacking. Every corrected LOO must use the same likelihood rows/order.
Fresh weights are estimated over the complete 26/18/18 candidate family within
each outcome/exposure. Weight on interaction candidates is not a probability
that moderation exists. Uncertainty uses a coherent 4,000-draw posterior mixture.

### 10 November integrated LOO review

All 248 integrated LOO files passed SHA-256, identity, finite-score and exact-case
coverage checks. There are 102 cases above k=0.7 across 45 models (maximum k=1.469):
relative strength 48 / whole body 10; supplementation strength 19 / whole body 6;
total intake strength 17 / whole body 2. No muscle or physical-function refits are
required. Each case is an observation deletion within its candidate model, not a
separate excluded study.

Parent-sized full exact archives estimate 98.51 GiB; a 1.5 allowance is 147.76 GiB
against 133.12 GiB free. Consequently `tst_exact_file` uses the existing shared
compact exact implementation. No initial exact archives existed when this choice
was made. Four concurrent in-memory refits, temporary files, compilation and
possible retries still need space. Passing compact checkpoints resume without
refitting. A scoring error before checkpoint creation requires that case to be
rerun; no full posterior can be restored from a compact checkpoint.

The training-only exact adapter now forwards the accepted candidate's retained
sampling length to the shared exact fitter. This gives 6000 draws per chain for
99 cases and 12000 for the three supplementation-strength model-14 cases. Both
use four chains, 2000 warmup and inherited sampler controls, priors and formula.
The exact diagnostic table records both warmup and retained length. Failed exact
diagnostics require a reviewed remedy before final LOO or stacking.

### 10 November exact diagnostic review

All 102 initial compact checkpoints were written. Sampling diagnostics pass for
98; four fail, all supplementation strength, model 14 or 26, observation 80 or
81. Model 14 has Rhat 1.100/1.068 and bulk ESS 25/50; observation 81 also has ten
divergences. Model 26 has Rhat 1.019/1.049 and bulk ESS 93/65, with no divergences.
No candidate or observation is removed.

`training-status-exact-remedies.csv` binds each reviewed retry to the SHA-256 of
its failed initial checkpoint. All four use 4000 warmup, 12000 retained draws,
adapt_delta 0.9995, maximum tree depth 15, four chains and seed 1988. Longer
adaptation and tighter integrator control address the mixing/geometry failures;
model 26 also doubles retained draws. Formula and prior changes for an isolated
case would change its predictive target and are not part of these remedies.

`tst_exact_review_audit_file` returns all 102 diagnostic records: 98 passing initial
paths unchanged, plus four compact retries in `exact_checkpoints/remedy_1/`.
Initial artifacts remain immutable. Passing retries resume by matching their
entire reviewed plan. Final correction now consumes these reviewed paths and
still refuses any failed or missing score. Failed retries require further review;
no exact fits were executed by the preparing agent.

### 11 November marginal exact-refit pilot

The first retry stopped with one unfinished case and three persistent
convergence failures. The source audit matched the original timepoint means
and SDs. The author confirmed that Figure 3 is the sole source of the two
digitised change SDs. The frozen observations and sampling variances are the
inputs to the pilot.

Run one case first, supplementation-strength model 26 at observation 80.
The new likelihood integrates the Gaussian arm, outcome and study-spline
coefficients over their existing priors within each study. The saved Stan
prior statements for the remaining parameters are retained. The spline basis
comes from the accepted full-fit design, as in integrated PSIS, and is fixed
when the held-out response is removed. Rebuilding the basis on the reduced
data can change the spline matrices. Prediction uses the conditional Gaussian
density for the held-out response given the remaining responses in its study.

`sensitivity/run_training_status_marginal_pilot.R --check` verifies the saved
case and fit checksums, design, held-out log-density agreement with the
production integrated scorer, and Stan syntax. `--fit` runs four chains with
2000 warmup, 6000 retained draws, adapt_delta 0.999 and maximum tree depth 15.
The pilot saves diagnostics and held-out log-density draws to
`results/update_2026/training_status/loo/marginal_pilot/`. Review the pilot
before adopting the marginal refit for the other cases or assembling final LOO.

The pilot passed with maximum Rhat 1.00393, minimum bulk ESS 2035, minimum tail
ESS 4327 and zero divergences or treedepth hits. The exact held-out score was
-9.4976 (Monte Carlo SE 0.1194 after restoring the saved draws' four-chain
layout). The `exact_marginal` stage validates and reuses
that checkpoint, then fits model 26 observation 81 and model 14 observations 80
and 81. Model 14 retains 12000 draws per chain, model 26 retains 6000, and all
use the pilot's 2000 warmup and sampler controls. The stage combines these four
marginal checkpoints with the 98 passing initial exact results for review.
Original archives, including the interrupted sampler-only retry, are retained.
Final LOO assembly consumes the reviewed selection and still requires a finite
passing exact score for every high-k case.

## Reporting and reproducibility

Report two shared profiles: age 25/week 12 and the restricted sample's median
arm age and median study/follow-up duration. Hold standardised body mass at zero
for absolute covariate-adjusted candidates. Report expected outcomes, protein
contrasts from the frozen reference dose, and trained minus untrained differences
in those contrasts. Population predictions exclude ordinary group effects and
subtract the study-specific factor smooth using the same posterior draw IDs.
Checks on relative model 14 and supplementation model 26 verified that changing
the placeholder study changes the raw prediction but leaves the population
prediction invariant to numerical precision. All status/profile/dose comparisons preserve posterior draw
and candidate pairing. The reference-based contrast curve does not represent all
possible current/proposed intake pairs. Retain full empirical group dose ranges
and study counts. A descriptive shared-range flag requires at least three studies
per group and both doses within both ranges; it does not guarantee joint support
at either age/duration profile. Physical function is always explicitly labelled
sparse and potentially prior-sensitive, even if a model happens to converge.

Generate the static public interface under `manuscript/generated/training-status/`:
predictions, protein contrasts, paired interaction contrasts, group support,
profiles, weights and figures. `supplement/training-status.qmd` includes those
files without reading fitted objects or the target store. Once numerical and
visual checks pass, review the generated diff and render locally:

```r
quarto::quarto_render("supplement/training-status.qmd")
```

Equivalent CLI: `quarto render supplement/training-status.qmd`. Review the render
before manuscript inclusion. No analytical pipeline or render was run during
preparation. `sensitivity/test_training_status.R` checks all actual Stan inputs
without compilation/sampling; `targets::tar_validate()` checks the separate graph.

Preparation verification: all 248 actual Stan inputs passed the finite/full-rank
fixed-design checks in R 4.4.3. Restricted-sample grids, absolute exposure scaling,
status-sensitive fingerprints, posterior contrast pairing and disabled-stage
guards passed. The separate target graph validated and all new R files parsed.

## Lean pruning plan for publication review

11 November 2026. Training-status fitting, corrected LOO, stacking and reporting
are complete. The dry-run manifests in `results/update_2026/training_status/pruning/`
cover the training-status model archive.

| Plan | Files | Size | Manifest |
| --- | ---: | ---: | --- |
| Propose removal of superseded initial fits | 2 | 1.53 GiB | `deletion_file_manifest.csv` |
| Retain accepted fits | 248 | 151.33 GiB | `retained_fit_manifest.csv` |
| Retain compact LOO archives and exact checkpoints | 707 | 18.77 MiB | `compact_archive_manifest.csv` |

The two candidates are the initial fits for absolute total intake, whole-body
model 17, and absolute supplementation, strength model 14. Both failed initial
diagnostics and have accepted replacements in `remedy_1/`. The plan rechecked
the four files' SHA-256 checksums against the diagnostic records. It also
checked all 248 accepted-fit paths against the combination and integrated LOO
manifests, verified all 12 reporting and verification bundles, and checked the
36 links in the generated supplement. The compact manifest records fresh
checksums for the retained archives. The accepted-fit manifest carries forward
the checksums verified at combination.

Keep the accepted fits through publication review so the prediction profiles,
contrasts or model mixtures can be revised. The saved tables and figures can
render from the static supplement interface. If a later plan removes accepted
fits, generating new posterior predictions would require restoring those
archives or refitting. Removing the two superseded files would also mean
restoring or refitting them to repeat their original fit-level diagnostics.
The plan preserves the diagnostic tables, checkpoints, target store, and
manuscript outputs. Its scope excludes the primary, age/duration and
muscle-measurement analyses.

The manifests are a snapshot for review. Recheck the proposed files and their
accepted replacements against the recorded sizes and checksums before any
deletion. The plan removed no files. Standalone rendering remains part of the
publication review.
