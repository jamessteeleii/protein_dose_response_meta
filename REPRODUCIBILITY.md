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
6. Restore the archived canonical and updated-search workbooks unchanged. The
   updated workbook is appended in memory; neither raw workbook is modified.
   Follow `data/update_2026/README.md` and `data/external/README.md` for input
   provenance and refresh rules.
7. Use a normal local clone with long-path support and sufficient temporary
   space. Machine-local cache and temporary roots may be set using the safe
   placeholders in `.Renviron.example`; never commit the real `.Renviron`.

At the 12 November documentation audit, the generic `Rscript` on this machine's
`PATH` resolved to R 4.6.1 even though the completed target processes recorded
R 4.4.3. XeLaTeX and `g++` were also not visible on that shell's `PATH`. These
observations do not invalidate the completed stored analyses, but they mean
the current shell is not evidence of a clean reproducibility environment.
Select and verify the required R 4.4.3 toolchain explicitly.

## Archived analysis inputs

The public release tree contains portable reporting files and curated inputs.
The original extraction and RoB workbooks, strategy evidence audit and full
frozen verification bundles remain in the research archive. Contact the
corresponding author for access and restore files at their recorded
project-relative paths before analytical
work. The filenames and checksums of key inputs are recorded in
`config/data-manifest.yml`. Verify the restored copies before running a stage:

```sh
Rscript scripts/verify_release_inputs.R --archived
```

`config/verification-manifest-index.csv` identifies the frozen bundle manifests
by SHA-256. Each full bundle's own manifest records its member checksums.
The three analysis workflows also need their archived accepted fits or the
documented recovery steps for any phase that consumes posterior objects.
The published reports can be rendered without those materials. Their portable
study-level effects exclude correspondence fields, and their public RoB
domain audit retains decisions without the original assessor notes.

## Required execution order

There is deliberately no supported "run everything" command. Build and review
one gated stage at a time in this order:

1. Main updated relative-intake primary and sensitivity families, absolute
   total-intake families, absolute supplementation families, publication-bias
   models, and both measurement sensitivities. Follow
   `data/update_2026/README.md`, `sensitivity/README.md`, and
   `sensitivity/measurement-analyses.md`.
2. Age/duration primary recovery, interaction fitting, LOO, combination, and
   reporting through `_targets_age_duration.R`. Follow
   `sensitivity/age-duration-interactions.md`.
3. Training-status fitting, exact corrections, combination, and reporting
   through `_targets_training_status.R`. Follow
   `sensitivity/training-status.md`.
4. Verify every frozen bundle in `results/verification/`, then review and
   commit the static public interfaces under `manuscript/generated/`.
5. Render manuscripts and supplements from the committed static interfaces in
   a clean clone without fitted models or target stores.

Every expensive phase is default-off. Set only the switch for the reviewed
stage, preserve successful branches, inspect its diagnostic checkpoint, and
unset the switch before moving on. Do not skip a preflight, diagnostic audit,
exact-refit review, verification bundle, or checksum gate.

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

Do not describe the project as independently reproducible until that
fresh-clone test has been completed and recorded. The local audit above proves
what this development machine executed; it does not itself prove that a new
machine can reproduce the complete suite without undeclared state.
