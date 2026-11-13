# Reproducibility and fitting-time audit

This document records the completed 2026 updated-search analysis suite and the
computational resources observed while producing it. It is the starting point
for anyone considering a complete rerun. The stage-specific scientific and
diagnostic decisions remain in the linked runbooks; this document does not
replace them.

## Do not start a full rerun casually

The completed suite is much larger than a single `targets::tar_make()` call.
It uses three independent target graphs, opt-in phase switches, reviewed
diagnostic remedies, exact case-deletion refits, and deliberate checkpoints
between stages.

The completed run involved:

- **2,662 distinct candidate models** contributing to the final analysis
  suite;
- **4,266 completed model-fitting jobs** after common-sample replacements,
  recovery fits, diagnostic remedies, and exact case-deletion refits are
  included;
- **1,659.04 model-job hours** (69.1 sequential days) of recorded fitting;
- **1,765.57 target-job hours** (73.6 sequential days) after diagnostics,
  LOO, stacking, predictions, reporting, and other recorded processing are
  included; and
- approximately **602.33 hours (25.1 days) of non-overlapping active target
  time** after concurrent branches are collapsed.

These are lower-bound observations, not a promise that another machine will
finish in the same time. Do not erase a target store, delete a fit archive, or
enable every phase switch in an attempt to obtain a one-command rebuild.

## Audit scope and definitions

The audit was completed on 12 November 2026 from the final local metadata in:

- `_targets/meta/meta` for the updated main, absolute-exposure,
  publication-bias, and muscle-measurement analyses;
- `_targets_age_duration/meta/meta` for primary recovery and age/duration
  interactions; and
- `_targets_training_status/meta/meta` for training-status interactions.

The standalone training-status marginal exact-refit pilot was added from
`results/update_2026/training_status/loo/marginal_pilot/diagnostics.csv`
because it ran outside the target graph and records 3,475.01 seconds.

The audited source identities were:

| Timing source | SHA-256 |
| --- | --- |
| Main metadata | `4268abb81e0b5f7de9b53f97e7ff71037c0af832804d0e171c68a299f6400377` |
| Age/duration metadata | `2a841bf2f4168f5923bc2ff1af6b5ccf88ef24678d50a342a165ceced34a6803` |
| Training-status metadata | `6fbb5f351e105aee436887e8ddae7dccff73d47db6ec6d211adca7f3028b84e9` |
| Marginal-pilot diagnostics | `142cda6cd1d41b08307239a8e3e77ead73e90966cf366200c7bd06da21409ad3` |

Dynamic-pattern summary rows were excluded because their `seconds` values are
aggregates of branch rows and would double count work. Function and object
definition rows were excluded. Completed branch and ordinary target timings
were summed. Intervals derived as `completion time - seconds` through
`completion time` were merged to estimate non-overlapping active wall time.

The audit covers the final updated-search suite from 28 September through
11 November 2026. It excludes superseded pre-update analyses completed before
28 September, package installation, source-data extraction, document editing,
and human review time. The `targets` metadata retains the latest successful
execution of a target, so failed, aborted, or overwritten attempts are not
fully represented. One known interrupted training-status retry is therefore
not timed. Treat every runtime below as a minimum.

### Distinct candidate models

"Candidate models" counts each model specification in the final scientific
ensembles once. It does not count a diagnostic remedy or an exact deletion as
a new scientific candidate.

| Analysis block | Candidate models |
| --- | ---: |
| Relative-intake primary plus eight sensitivity families | 936 |
| Absolute total-intake primary plus eight sensitivity families | 648 |
| Absolute supplementation primary plus five sensitivity families | 432 |
| Two muscle-measurement sensitivity ensembles | 52 |
| Publication-bias models | 10 |
| New age/duration interaction candidates | 336 |
| Training-status interaction candidates | 248 |
| **Total** | **2,662** |

The age/duration ensembles also reuse 248 primary candidates. Seventy-two
absolute-primary candidates had to be recovered because their full fit
archives had previously been pruned. Those recovery runs are included in the
job and time totals below but not counted as new candidate specifications.

### Completed fitting jobs and elapsed branch time

"Fitting jobs" counts every completed model execution, including initial
candidates, replacement fits on aligned samples, diagnostic remedies, recovery
runs, exact case deletions, and reviewed exact-refit remedies.

| Workstream | Completed fitting jobs | Model-job hours |
| --- | ---: | ---: |
| Relative-intake primary and sensitivities | 1,262 | 302.69 |
| Absolute total intake | 924 | 525.76 |
| Absolute supplementation | 813 | 420.62 |
| Muscle-measurement sensitivities | 64 | 7.93 |
| Publication-bias models | 10 | 0.07 |
| Absolute-primary recovery for age/duration work | 73 | 19.05 |
| Age/duration interaction workflow | 761 | 281.50 |
| Training-status workflow, including the standalone marginal pilot | 359 | 101.42 |
| **Total** | **4,266** | **1,659.04** |

The large difference between candidate models and fitting jobs is expected.
Pareto-k diagnostics required many observation-deletion fits, especially for
the residual-correlation and interaction analyses. Remedy fits preserve the
candidate identity while changing only reviewed sampler controls or retained
sampling length unless the relevant runbook explicitly records another
decision.

### All recorded target processing

| Target graph | Recorded target-job hours |
| --- | ---: |
| Main graph (`_targets.R`) | 1,325.69 |
| Age/duration graph (`_targets_age_duration.R`) | 328.11 |
| Training-status graph (`_targets_training_status.R`) plus standalone pilot | 111.78 |
| **Total** | **1,765.57** |

Model fitting accounts for 1,659.04 hours. The remaining 106.53 hours cover
recorded diagnostics, integrated LOO, exact scoring, stacking, posterior
prediction, reporting, verification, and other pipeline processing.

Several model branches ran concurrently, and each Bayesian fit ordinarily ran
four chains in parallel. Summed branch time is therefore not calendar time and
is not CPU-core time. Collapsing overlapping target intervals gives 602.33
hours of active target time. The observed calendar window was 44.2 days,
including review pauses and periods with no computation.

## Reference hardware and storage

The timed work ran under R 4.4.3 on Windows on an AMD Ryzen 9 3900XT
(12 physical cores, 24 logical processors). The usual schedule allowed four
model branches, each with four parallel chains, for up to 16 simultaneous
chains. Runtime will vary materially with processor speed, memory, storage,
operating system, compiler, and the selected worker count.

At the audit, the retained local computational footprint was:

| Location | Files | GiB |
| --- | ---: | ---: |
| `model_archive/` | 7,231 | 702.53 |
| Main target store | 10,820 | 43.42 |
| Age/duration target store | 1,323 | 0.04 |
| Training-status target store | 473 | 0.02 |
| Updated results, verification bundles, and generated interfaces | 3,226 | 1.44 |
| **Observed retained total** | **23,073** | **747.44** |

This is not peak usage. Several completed family archives had already been
pruned, while a clean rebuild may temporarily hold initial fits, remedies,
exact refits, compiler output, and replacement fits together. A complete rerun
should use a fast local disk with **at least 1.2 TiB free** at the start. Do not
place the repository, `model_archive/`, any `_targets*/` store, `renv/library/`,
the global `renv` cache, or temporary Stan files in Dropbox, OneDrive, iCloud,
Google Drive, or another filesystem-sync directory.

At least 64 GiB RAM is recommended for the documented four-worker schedule.
On a smaller machine, reduce `PROTEIN_MODEL_WORKERS` to one or two. This lowers
peak memory and simultaneous chains but increases wall time. Never run two
processes against the same target store concurrently.

## Software and input requirements

Use the versions and inputs recorded by the repository, not current latest
versions:

1. R 4.4.3, matching `renv.lock`. Confirm `R.version.string` before restoring;
   a different `Rscript` earlier on `PATH` is not an acceptable substitute.
2. Restore R packages with `renv::restore()` from the project root. Keep the
   project library and global cache outside synced storage.
3. Use a C/C++ toolchain compatible with R 4.4.3 (Rtools44 on Windows).
4. Use CmdStan 2.39.0 for the CmdStanR analyses and verify the toolchain with
   `cmdstanr::check_cmdstan_toolchain()`.
5. Use Quarto 1.7 or later and a TeX installation providing XeLaTeX for PDF
   rendering. RStudio's bundled Quarto is acceptable when detected.
6. Use the frozen workbooks and reviewed CSVs bundled under `data/`.
   The updated workbook is appended in memory. See the
   [data guide](data/README.md) for file identities and provenance.
7. Use a normal local clone with long-path support and sufficient temporary
   space. Machine-local cache and temporary roots may be set using the safe
   placeholders in `.Renviron.example`; never commit the real `.Renviron`.

At the 12 November documentation audit, the generic `Rscript` on this machine's
`PATH` resolved to R 4.6.1 even though the completed target processes recorded
R 4.4.3. XeLaTeX and `g++` were also not visible on that shell's `PATH`. These
observations do not invalidate the completed stored analyses, but they mean
the current shell is not evidence of a clean reproducibility environment.
Select and verify the required R 4.4.3 toolchain explicitly.

## Start from a clean clone

Clone this repository into a local project directory. Use the R 4.4.3
executable for every command. The project startup file bootstraps the pinned
`renv` version. Install the system prerequisites above and the JAGS
installation required by RoBMA before restoring the package library.

```sh
git clone https://github.com/jamessteeleii/protein_dose_response_meta.git
cd protein_dose_response_meta
Rscript -e 'renv::restore(); renv::status()'
Rscript scripts/verify_release_inputs.R
Rscript scripts/check_analysis_inputs.R
```

The checks read the public workbooks, prepare the four primary outcome
datasets, compare the effect sizes, variances, and key model inputs with
the published data, and validate the strategy and RoB 2 inputs without sampling.
The expected eligible contrast counts are 283 whole-body, 425 direct muscle,
794 strength, and 362 physical function.

For the first pipeline checkpoint, start R in the project root and run:

```r
Sys.setenv(PROTEIN_MODEL_WORKERS = "1")
for (script in c("_targets.R", "_targets_age_duration.R", "_targets_training_status.R")) {
  targets::tar_validate(script = script)
}
targets::tar_make(c(
  update_2026_workbook_checks,
  update_2026_combined_primary_preparation_audit_files,
  update_2026_strategy_audit_review_files
))
```

The first results are written under `results/update_2026/primary_refit/`
and `results/update_2026/sensitivity/strategy/audit/`. Continue with the
stages below. Keep successful fits and reductions until all dependent
analyses and reports have finished.

## Archived analysis inputs

The [input manifest](config/data-manifest.yml) identifies the bundled files
and records the checksums of the historical sources. The full original
extraction and RoB workbooks, source-text strategy audit, and exact historical
fit and verification bundles are held in the research archive. These are
useful for provenance and resuming the historical analyses. The public
workflow uses the bundled extraction data, frozen pre-update reference
sample, reviewed strategy classifications, and final RoB 2 judgments.
`config/verification-manifest-index.csv` identifies the original verification
bundles. Each new run writes its own manifests and checksums.

## Required execution order

There is deliberately no supported "run everything" command. Build and review
one gated stage at a time in this order:

1. Build the [updated relative-intake primary analysis](data/update_2026/README.md#local-primary-refit-runbook),
   including corrected LOO, reporting, and
   `update_2026_primary_verification_files`. That target creates the parent
   checkpoint for [both muscle-measurement analyses](sensitivity/measurement-analyses.md).
   Build the relative, total-intake, and supplementation sensitivity families
   using the [analysis stage index](#analysis-stage-index) and the documented
   per-family review gates. Run publication bias before
   [GRADE](data/GRADE/README.md#run-after-the-current-model-pipeline-has-finished).
2. Retain the relative, total-intake, and supplementation primary fit files
   and final LOO reductions. Their three `*_primary_verification_files`
   targets generate the manifests read by the
   [age/duration workflow](sensitivity/age-duration-interactions.md).
   On a clean rebuild with all accepted total-intake fits present, run
   interaction preparation directly. The recovery stages reproduce the
   historical refit of the pruned total-intake family when those fits are absent.
3. Training-status fitting, exact corrections, combination, and reporting
   through `_targets_training_status.R`. Follow
   `sensitivity/training-status.md`.
4. Verify the new bundles in `results/verification/`. Refresh the reporting
   tables with `Rscript scripts/materialize_preprint_reporting.R` and the
   PRISMA includes with `Rscript supplement/build-prisma-includes.R`. Compare
   the rebuilt values and diagnostic dispositions with the committed reports.
5. Render manuscripts and supplements from the committed static interfaces in
   a clean clone without fitted models or target stores.

Every expensive phase is default-off. Set only the switch for the reviewed
stage, preserve successful branches, inspect its diagnostic checkpoint, and
unset the switch before moving on. Do not skip a preflight, diagnostic audit,
exact-refit review, verification bundle, or checksum gate.

## Analysis stage index

Use the [stage index](sensitivity/analysis-stages.md) to find the targets and
environment switches for each reported analysis. The index includes the
relative-intake, absolute total-intake, and prescribed supplementation
families, and the reporting and verification checkpoints for the separate
measurement, age/duration, and training-status workflows.

## Fresh-clone release test

Before release, record a clean-clone test that demonstrates:

- exact package restoration under R 4.4.3;
- discovery of the compiler, CmdStan 2.39.0, Quarto, and XeLaTeX;
- validation of all three target graphs without executing expensive stages;
- identity and checksum verification of required source data and frozen
  prerequisite bundles;
- successful representative lightweight targets;
- successful standalone rendering from committed static interfaces; and
- either a complete staged rebuild or an explicit account of which expensive
  stages were verified from frozen bundles instead.

### Evidence from the candidate checkout

An isolated checkout of public master plus this data and reproducibility
update started with empty project libraries, target stores, model archives,
and analysis output directories. Under R 4.4.3, the pinned library restored
and `renv::status()` reported a consistent project. Rtools44, CmdStan 2.39.0,
JAGS 4.3.2, Quarto and XeLaTeX were detected. All three target graphs validated.
The input check reproduced the four eligible contrast counts and matched
the published primary effect sizes, variances, intake, RT coding, and sample
sizes for all 1,864 rows. It covered all 95 studies in the updated strategy
audit and verified RoB 2 coverage for GRADE. The recomputed coefficient priors
and 312 primary and strategy specifications matched the accepted recipes.
All 248 initial primary Stan inputs across the three exposure scales passed
the input checks without compilation or sampling.
The PRISMA includes rebuilt to the same contents. The manuscript and all 15
HTML supplements rendered, and the package checker verified the linked
downloads, including the two workbooks.

### Verify the frozen results or rerun the analysis

The repository provides the extraction data, analysis code, and frozen
reporting files used by the manuscript and supplements. Readers can check
the reported values against those files and render the documents from them.
The research archive holds the complete historical verification bundles and
retained fitted model objects. The candidate checkout was checked through
data preparation, initial primary model inputs, and standalone rendering.
A full clean-clone refit and post-sampling rebuild were outside that check.

For an independent rerun, follow the
[analysis stage index](sensitivity/analysis-stages.md) from an empty analysis
checkout. Run `Rscript sensitivity/audit_verification_bundles.R` to verify the
new bundles, then compare the analysis samples, diagnostic decisions,
exact-deletion coverage, stacking weights, and reported tables with the
published outputs.

The workflow sets seeds for sampling and reporting. Exact agreement also
depends on the software, compiler, hardware, and sampling configuration.
See the [Stan reproducibility guide](https://mc-stan.org/docs/reference-manual/reproducibility.html).

The recorded remedy ledgers describe specific historical diagnostic failures.
If a fresh run produces a different failure set, retain the run outputs and
review the discrepancy before updating any ledger. Exact historical posterior
draws and render bytes remain research-archive artifacts.
