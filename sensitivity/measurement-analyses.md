# Muscle measurement analyses

Prepared 27 October and completed 2 November 2026, after the project author
confirmed that the other extended model families had finished. The two final
26-candidate ensembles required 64 fitting/refitting jobs, including two
diagnostic remedies and ten exact deletions, and 7.93 recorded model-job hours.
See the root [reproducibility and fitting-time audit](../REPRODUCIBILITY.md)
before rerunning the staged commands below.

## Completed state

- `prepare`: complete for both variants. The audited datasets contain 413
  contrasts from 52 studies and 143 arms (`prereg_methods`) and 255 contrasts
  from 29 studies and 84 arms (`direct_only`). All 52 Stan inputs passed.
- `pilot`: complete; all four full-schedule fits passed sampling and convergence
  checks and were reused in the full ensembles.
- `fit` and initial diagnostics: complete for all 52 candidates. All 26
  `prereg_methods` candidates passed. In `direct_only`, candidates 10 and 25 had
  two and one divergent transitions respectively; all convergence checks passed.
- `remedy`: complete for `direct_only` candidates 10 and 25 only, with
  `adapt_delta = 0.999` and `max_treedepth = 15`. The formula, priors, and full
  sampling schedule are retained. The decision ledger records the initial
  controls and observed reason for each refit. All 52 accepted candidates passed
  sampling/convergence checks with zero divergences or treedepth hits.
- PPC warnings remain in the audit. Initial 90% predictive interval coverage
  ranged from 93.7% to 95.2% for `prereg_methods` and 94.5% to 96.5% for
  `direct_only`; 20 and 12 observations respectively were flagged in at least
  one candidate. These flags are not sampler-control failures.
- LOO and exact corrections: complete. Five model-observation cases per variant
  required exact refits; all ten passed. Corrections cover precisely the flagged
  cases. Maximum remaining PSIS k is 0.6833 (`prereg_methods`) and 0.6651
  (`direct_only`).
- Reporting: complete. Both new 24-artifact bundles and the original 28-artifact
  primary bundle passed file-size and SHA-256 checks. Both new bundles passed all
  five verification checks, including complete candidates, accepted sampling
  diagnostics, corrected LOO, stacking weights, and support flags.
- Figure review: the registered-method figure passed visual inspection. The
  direct-only title was clipped; its source title has been shortened to
  "changes in muscle size", with MRI/CT/ultrasound/biopsy details retained in
  the caption. Regenerate the direct-only figure target and review it before
  manuscript use. This change does not require posterior reporting or refits.

The first accepted-manifest attempt stopped because the audit expected the full
extraction table inside each saved brms fit. brms retains only model-frame
columns. The audit now reconstructs predictors from the frozen selection and
checks every retained column, including transformed moderators and piecewise
basis matrices, row order, responses, variances, and grouping identities. The
complete selection fingerprint remains in the manifest. Saved fit files are
reused; this audit correction does not alter model specifications or sampling.

## Two prespecified-for-this-rerun datasets

Both analyses reuse the completed 2026 primary **26-candidate muscle-outcome
framework**, including its final accepted formulas, priors, sampler controls,
four chains, 2,000 warm-up iterations, and 6,000 retained iterations per chain.
They estimate new stacking weights on their own retained observations. The
additional measurement analyses were decided after registration; that timing
should remain explicit even where eligibility restores the registered method
list.

| Key | Retained measurements | Current contrasts |
| --- | --- | ---: |
| `prereg_methods` | Original muscle-group imaging/biopsy and limb DXA/BIA estimates, excluding every trunk and skinfold row | 413 |
| `direct_only` | Only muscle-group MRI, CT, ultrasound, and muscle biopsy measurements | 255 |

The original 425-row group remains intact for comparison. The first selection
excludes six trunk rows (studies 167 and 172) and six skinfold rows (studies
123 and 124). The second excludes all DXA/BIA and skinfold proxies. No
whole-body, strength, function, or extended-family fits are altered. Eligible
direct measurements are retained from studies that also report excluded
proxies.

The first key refers to the requested registered-method restriction. The two
generically labelled study-11 DXA rows (outcomes 66 and 72; arms 14 and 15)
were confirmed appendicular by the project author on 27 October 2026.
[measurement-annotations.csv](measurement-annotations.csv) records that
confirmation and its provenance. It does not alter either the workbook or
the frozen checkpoint. Both rows remain included in `prereg_methods` and
excluded from `direct_only`.

## Inputs and isolation

- File targets declare the manifest, analysis rows, and final specifications
  in `results/verification/primary_update_2026/`. Both input files must match
  their manifest SHA-256 checksums. Parsing considers all rows and fails on
  data loss. The checkpoint is the frozen parent analysis, not a mutable raw
  workbook or a lookup of the latest result.
- A separate file target declares `sensitivity/measurement-annotations.csv`.
  Annotations match study, arm, outcome, and both time points, require the
  expected original method label, and fail on missing or ambiguous matches.
  The new selection audit and analysis copy carry both the original label
  and `measurement_method_reviewed`, plus confirmation provenance. The
  original `measurement_method` column is preserved. The annotation table is
  included in preparation reports and final verification bundles. Existing
  target definitions and original input files do not depend on this new file.
- Retained `yi`, `vi`, study/arm/outcome identifiers, time contrasts, protein
  intake values, and original intake centring are unchanged. Factors are
  reconstructed after CSV import. Moderator transformations are built by the
  existing primary preparation functions on the retained dataset, and the
  age-25/week-12 prediction profile is retained.
- Fit branches return file paths only. Fit objects are stored once under
  ignored `model_archive/update_2026/measurement/<key>/`. Pilot and remedy
  stages have distinct subdirectories. Accepted manifests contain paths,
  byte sizes, and checksums; no target serialises a collected model ensemble.
- Preparatory audits and timings go under
  `results/update_2026/measurement/<key>/`. Final reporting bundles go under
  `results/verification/measurement_<key>_update_2026/`; figures go under
  `plots/update_2026/measurement/<key>/`. Existing outputs are not overwritten.
- `PROTEIN_MEASUREMENT_PHASES=none` is the default. Each expensive phase has
  its own cue and an execution guard inside the helper. An absent phase
  cannot silently launch a fit, exact refit, LOO reduction, or posterior
  prediction extraction.

## Run in reviewed stages

Use **R 4.4.3 and the project renv library**, as for the other families. On
the preparation machine the unqualified shell `Rscript` points to R 4.6.1;
select the 4.4.3 executable explicitly rather than changing the lockfile or
loading a different package environment.

From the project root, with the appropriate executable selected:

```text
Rscript sensitivity/run_measurement.R prepare
Rscript sensitivity/run_measurement.R pilot
```

`prepare` builds audited selections and all 52 Stan input matrices without
sampling. Inspect the selection and warning ledgers. The input checks retain
the primary model-specific priors. Some candidates report that a redundant
global coefficient prior is unused because all coefficients already have
individual priors; this is recorded in the preflight CSV rather than silently
changing those priors.

`pilot` runs **candidates 2 and 26 for each dataset** on the full primary
sampling schedule, records wall time, and applies the unchanged primary
diagnostic gates. Inspect the performance and diagnostic audit before full
fan-out. These four fits become part of the final ensembles and are not fit
again during the full phase. A failing pilot blocks full fan-out until its
specification is appropriately remedied and the decision documented.

For the current handoff, `prepare` is already complete. Run only `pilot`, then
return these two checkpoint files for review before running `fit`:

```text
results/update_2026/measurement/prereg_methods/performance_pilot/audit.csv
results/update_2026/measurement/direct_only/performance_pilot/audit.csv
```

```text
Rscript sensitivity/run_measurement.R fit
Rscript sensitivity/run_measurement.R diagnostics
```

`fit` runs the other 24 candidates per set. The runner uses
`tar_make_future()` with its normal clean child process and
`PROTEIN_MODEL_WORKERS` (default four), matching the existing four-model,
four-chain schedule. Thus there are 52 initial model fits in total, including
the four pilot fits. Successful branches are reusable after interruptions.
Read the initial diagnostic/PPC audit before continuing; do not discard
candidates or reuse the broad-group stacking weights to avoid a problem.

```text
Rscript sensitivity/run_measurement.R remedy
Rscript sensitivity/run_measurement.R loo
```

`remedy` reads the initially empty tracked
`sensitivity/measurement-remedies.csv`. With no new remedies required it
writes a lightweight no-remedy marker. Otherwise, add a separate row for
each documented candidate-specific increase in sampler controls:

```text
variant,model_index,stage,adapt_delta,max_treedepth,rationale
prereg_methods,13,remedy_1,0.9995,15,<actual diagnostic reason>
```

The example is a format illustration, **not a prescribed remedy**. Do not
add it in advance. New stages preserve earlier fits. The current hook handles
sampler-control remedies; geometry, prior, formula, or sampling-length changes
must be implemented and documented after inspecting the actual diagnostics,
using the same smallest-change process as the other families. No such
unobserved diagnostic decision has been guessed during preparation.

Accepted fits receive a fresh complete diagnostic audit. `loo` is blocked
unless all 26 candidates pass sampling/convergence checks. PPC warnings
remain visible rather than automatically excluding observations or models.
The LOO stage uses the primary conditional audit and Gaussian integrated
observation-level PSIS implementation with one internal core per branch.
The retained observation manifests must agree across all candidates.

```text
Rscript sensitivity/run_measurement.R exact
Rscript sensitivity/run_measurement.R report
```

`exact` selects the newly observed Pareto-k > 0.7 cases, performs exact
case-deletion refits, and assembles final LOO reductions. It supports zero,
one, or multiple cases per model rather than inheriting the old analysis's
case IDs or counts. Exact-refit diagnostics must pass before stacking. A
failed exact refit requires a documented targeted remedy; the system stops
rather than presenting an incomplete ensemble as final.

`report` estimates new 26-model stacking weights, streams accepted fits one
at a time into a 4,000-draw posterior mixture, and writes curves, pairwise
SESOI probabilities, empirical support, local equivalence, and figures. It
uses the primary 0.1-SD SESOI, 95% criterion, 0.01-g/kg/day grid, and reference
profile. It then writes and checksum-verifies each reporting checkpoint.
New figure titles distinguish the registered-method group from imaging and
biopsy measures. No current manuscript result is automatically replaced.

The existing extended muscle-model families, publication-bias analysis, and
GRADE appraisal describe the original broader outcome set. Keep the original
425-row primary report as the matched baseline for those extended-family
comparisons. A comparison with the corrected 413-row primary would otherwise
change both measurement eligibility and the modelling assumption. Identify
the outcome set in every table/caption, and qualify robustness statements
accordingly. Before replacing primary manuscript results, update the muscle
dataset counts and assess the implications for its GRADE and publication-bias
reporting. These downstream reporting/appraisal updates do not schedule
refits of the other extended families here.

After completion, compare the original 425-row, 413-row, and 255-row reports
within common RT-specific exposure support. Report study/arm/contrast counts,
curve differences and uncertainty, SESOI probabilities, and changed decision
regions. The existing broad cross-family plotting script has a fixed family
registry; it will not silently insert these single-outcome bundles. A final
manuscript comparison should be added after the two reports are validated.

Run the existing checksum audit before any later pruning:

```text
Rscript sensitivity/audit_verification_bundles.R
```

No pruning or cache deletion is scheduled by these new targets.

## Preparation checks

```text
Rscript sensitivity/test_measurement_sensitivity.R
Rscript sensitivity/test_measurement_sensitivity.R --preflight
```

The tests check both eligibility boundaries, unchanged numerical inputs and
identifiers, unknown-method failures, all candidate identities, disabled
execution guards, and synthetic zero/multiple-case exact-LOO correction and
stacking. `--preflight` reads the small cached `steele_priors` object without
building its target and constructs the 52 model inputs without compilation
or posterior sampling. For a clone without that prior object, use the future
`prepare` stage to build the declared upstream prior dependency.

During preparation, these sampling-free tests and all 52 input constructions
passed under R 4.4.3. The full targets graph also passed `tar_validate()`.
Posterior fitting, actual fit/LOO diagnostics, posterior reports, and rendered
figures remain to be checked after the deferred runs.

## Implementation files

- [R/measurement_sensitivity.R](../R/measurement_sensitivity.R): audited
  selection, isolated fits, diagnostic gates, LOO/exact corrections, streaming
  reporting, figures, and verified checkpoints.
- [R/measurement_targets.R](../R/measurement_targets.R) and
  [_targets.R](../_targets.R): identical staged target definitions for both sets.
- [R/incremental_update.R](../R/incremental_update.R): optional output paths
  for reused primary reporting helpers and an optional general design input
  for exact scoring. Existing defaults and estimands are unchanged.
- [run_measurement.R](run_measurement.R),
  [measurement-remedies.csv](measurement-remedies.csv),
  [measurement-annotations.csv](measurement-annotations.csv), and
  [test_measurement_sensitivity.R](test_measurement_sensitivity.R): deferred
  runner, initially empty remedy ledger, author-confirmation ledger, and sampling-free checks.
- [README.md](README.md) links this guide; [.gitignore](../.gitignore) excludes
  new local working reports and plots while preserving final verification bundles.
- [Manuscript source](../manuscript/manuscript.qmd) and
  [terminology review](../manuscript/outcome-terminology-review.md): source-only
  writing suggestions, plot-label checklist, and the agreed deferred analysis plan.
