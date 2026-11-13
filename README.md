# Protein dose-response meta-analysis

This repository contains the data, analysis code, figures, and manuscript
sources for *Protein intake and changes in whole-body lean/fat-free mass,
measures of muscle size or mass, strength, and physical function in healthy
adults: an updated systematic review and Bayesian model averaged dose-response
meta-analysis*.

The project was preregistered on the
[Open Science Framework](https://osf.io/w2qpx/overview). The originally
preregistered models remain the primary analyses. Measurement-error Stan
models were completed as supplementary analyses and remain clearly separated
in the code and manuscript.

## Citation

Please cite as:

> Steele, J., Nunes, E.A., Moore, M., Lees, M., Roxburgh, L., Viloria Medina, A., and Phillips, S. (2026). Protein intake and changes in whole-body lean/fat-free mass, measures of muscle size or mass, strength, and physical function in healthy adults: an updated systematic review and Bayesian model averaged dose-response meta-analysis. SportRxiv. <https://doi.org/10.51224/SportRxiv.1166>

## Supplementary materials

The [supplementary materials](https://jamessteeleii.github.io/protein_dose_response_meta/) include the analyses, plots, and supporting data.

## Repository structure

- `_targets.R` defines the main analysis pipeline.
- `_targets_age_duration.R` and `_targets_training_status.R` define separate,
  opt-in supplementary pipelines; see the workflow overview below.
- `R/functions.R` contains the functions used by the pipeline.
- `config/data-manifest.yml` records bundled and archived analysis inputs.
- `data/graphical_data/` contains digitised numerical extracts.
- `data/external/` contains small pinned external inputs and their provenance.
- `stan/` contains supplementary Stan models and supporting scripts.
- `sensitivity/` documents the opt-in sensitivity pilot targets and runner.
- `results/freeze/` records the earlier analysis freeze. Full verification
  bundles are part of the research archive.
- `manuscript/` contains the Quarto manuscript and portable reporting files.
- `supplement/` contains the supplementary Quarto sources and public audit tables.
- `renv.lock` records the R package environment.

Git contains the source and the compact public reporting interface. The
release package contains the rendered PDF, supplements, and linked downloads.
See [the archive instructions](REPRODUCIBILITY.md#archived-analysis-inputs)
for the original workbooks and frozen analytical checkpoints.

> **Full-rerun warning:** the completed suite contains **2,662 distinct
> candidate models** and required **4,266 completed fitting/refitting jobs**.
> The audit records at least **1,659 model-job hours**, **1,766 total target-job
> hours**, about **25 days of active wall-clock computation**, and a retained
> local footprint of about **747 GiB**. Do not attempt an indiscriminate full
> rebuild. Read [the reproducibility and fitting-time audit](REPRODUCIBILITY.md)
> and follow the reviewed stage order below.

## Requirements

- R 4.4.3, matching `renv.lock`.
- A C/C++ toolchain compatible with R 4.4.3. On Windows, use Rtools44.
- Quarto 1.7 or later, as required by the bundled manuscript extensions.
- A TeX distribution with XeLaTeX to render the manuscript PDF; TinyTeX is a
  suitable option.
- Internet access for package restoration and external data downloads.
- CmdStan 2.39.0 to run the supplementary CmdStanR scripts.
- JAGS 4.3.1 or later for the recorded RoBMA environment.
- A fast local filesystem with at least 1.2 TiB free for a clean full rerun.
  The observed retained footprint was about 747 GiB after some completed fit
  archives had already been pruned; peak usage is higher.
- At least 64 GiB RAM is recommended for the documented four-worker schedule.
  Use one or two workers on smaller machines.

Keep the clone, target stores, model archive, temporary directory, project
library, and global `renv` cache outside all filesystem-sync directories. See
[`.Renviron.example`](.Renviron.example) for optional machine-local cache and
temporary-path placeholders.

RStudio is optional, but opening `protein_dose_response_meta.Rproj` is a
convenient way to start in the correct project directory. Recent RStudio
versions include a compatible Quarto installation.

## Restore the R environment

Clone or download the repository, start R 4.4.3 in the project root, and run:

```r
install.packages("renv")
renv::restore()
```

Direct R package dependencies are declared in `DESCRIPTION`, and `renv` uses
explicit snapshots so it does not scan data and generated outputs for package
references. When adding or removing a direct dependency, update `DESCRIPTION`,
test the change, and then run `renv::snapshot()` to refresh `renv.lock`.

The restore may compile packages and can take some time. The supplementary
Stan scripts additionally require a local CmdStan installation. To reproduce
the development environment used here:

```r
cmdstanr::install_cmdstan(version = "2.39.0")
cmdstanr::check_cmdstan_toolchain()
```

## Run the analysis pipeline

The project uses [`targets`](https://books.ropensci.org/targets/) to track
dependencies and avoid unnecessary recomputation.

```r
targets::tar_manifest()
```

Manifest inspection is non-sampling. An ordinary `targets::tar_make()` is not
a complete reproduction command: expensive families are deliberately
default-off, and two completed workflows use separate scripts and stores. A
clean full reproduction is a staged operation measured in weeks and hundreds
of gigabytes. Follow the execution order in
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md), enable only one reviewed phase,
and use the `names` argument to request its documented checkpoint. Preserve
successful branches and never erase a store to resume an interrupted run.

Model-fitting targets support two levels of parallelism. `_targets.R` declares
a multisession plan with four model workers by default, while each Bayesian
model retains its four parallel chains. Keep the `tar_make_future()` worker
limit aligned with the plan. The following is the historical pattern for the
updated-primary stage; do not run it unless that specific refit has been
deliberately selected from the reproducibility runbook:

```r
model_workers <- 4L
Sys.setenv(
  PROTEIN_MODEL_WORKERS = model_workers,
  PROTEIN_UPDATE_2026_PRIMARY = "fit"
)
targets::tar_make_future(
  update_2026_primary_fit_file,
  workers = model_workers,
  reporter = "balanced"
)
```

This runs at most four model branches simultaneously (16 chains in total).
Completed and errored models remain separate dynamic branches, so later runs
can skip successful fits and revisit failures independently.

## Separate age/duration and training-status workflows

These post-registration supplementary analyses intentionally retain separate
`targets` scripts, local stores, and clean R processes. This arrangement was
retained on 3 November 2026 to preserve the existing workflows and completed
computational state. Their execution controls keep supplementary work separate
from the frozen primary analyses. Training-status adapters also replace some
age/duration helper functions within their own process, so do not source
`R/training_status_workflow.R` into the main or age/duration pipeline session.

| Workflow | Script | Local store | Stage instructions |
| --- | --- | --- | --- |
| Main analyses and measurement sensitivities | `_targets.R` | `_targets/` | This README and [measurement analyses](sensitivity/measurement-analyses.md) |
| Age and duration interactions, including total-primary recovery | `_targets_age_duration.R` | `_targets_age_duration/` | [Age/duration workflow](sensitivity/age-duration-interactions.md) |
| Reported training status in RT studies | `_targets_training_status.R` | `_targets_training_status/` | [Training-status workflow](sensitivity/training-status.md) |

An ordinary `targets::tar_make()` uses the main script and does **not** run
either separate workflow. Full reproduction of these supplementary analyses
requires their documented stage commands as well as the main pipeline. This
is a documented exception to the project's default single-entry-point rebuild.
Dependencies between workflows are checked through frozen files and SHA-256
manifests; they are not automatic scheduling edges in the main target graph.
A missing prerequisite stops execution rather than launching another pipeline.

The execution order is:

1. Produce the frozen relative, absolute total-intake and supplementation
   primary checkpoints, and complete both muscle-measurement verification
   bundles using the main pipeline.
2. Restore missing total-primary fit archives from a verified backup, or
   complete the reviewed recovery stages in the age/duration workflow before
   interaction fitting and combination. See its recovery instructions for the
   required diagnostic and LOO checks.
3. Complete the age/duration stages through reporting.
4. Complete the training-status stages; expensive stages verify the completed
   age/duration reports as well as both measurement bundles.

Use R 4.4.3 with the restored project environment, from the repository root.
The shell runners select the correct script, store and phase. For example,
these commands perform preparation only, at the appropriate point above:

```sh
Rscript sensitivity/run_age_duration.R recover_prepare
Rscript sensitivity/run_age_duration.R prepare
Rscript sensitivity/run_training_status.R prepare
```

Run one reviewed stage at a time using the linked stage tables, rather than
running this example block as a batch. Both phase switches default to `none`;
later fitting, LOO, combination and reporting stages require explicit selection.
Resume through the same runner and store, preserving successful branches and
archives; do not merge, rename, manually edit or erase the stores to resume.

Each workflow writes its single-copy fit archives under
`model_archive/update_2026/{age_duration,training_status}/` and its diagnostic
and reporting checkpoints under `results/update_2026/{age_duration,training_status}/`.
Those local products and all three stores are ignored by Git. Public report
inputs are generated under `manuscript/generated/age-duration-interactions/`
and `manuscript/generated/training-status/`; review and commit validated
changes to these static interfaces. With those interfaces present, render
the supplements without fitted models or target stores:

```sh
quarto render supplement/age-duration-interactions.qmd
quarto render supplement/training-status.qmd
```

## Inspect results without refitting

The public reporting interface supports inspection without fitting. Full
checkpoints in the research archive restore to `results/verification/` and
can be read without a local target store. Each family has a manifest with
dimensions, descriptions, byte sizes and SHA-256 checksums. Tabular artifacts are gzip-compressed CSV files; non-tabular R
objects use `qs`.

After completing or changing the primary analysis, refresh its checkpoint:

```r
targets::tar_make(primary_verification_files)
```

The primary bundle includes the complete-case analysis rows, candidate-model
registry, sampler and posterior-predictive diagnostics, corrected LOO audit,
stacking weights, dose-response summaries, every plotted pairwise contrast
and SESOI probability, local-equivalence summaries, and raw-unit reference
values. It is sufficient to audit manuscript numbers and reconstruct the
reporting layer. It is not a substitute for an independent posterior refit.

An equivalent bundle should be added at the end of each completed sensitivity
family before that family's fit cache is removed. Large fitted-model archives,
if retained for exact computational reproduction, should be deposited in
versioned external storage with checksums rather than committed to Git.
Remaining full sensitivity families use file-backed fit branches and
path-only ensemble manifests so that a fit is not serialised repeatedly in
collected-fit, ensemble, and LOO-job targets. The detailed storage invariants
are documented in [`sensitivity/README.md`](sensitivity/README.md).

## Prune the local cache

`sensitivity/prune_targets.R` is dry-run-first. It refuses family-specific
pruning unless the corresponding reporting manifest exists and every bundle
checksum is valid. Preview obsolete objects or completed primary-analysis fit
objects with:

```sh
Rscript sensitivity/prune_targets.R --family=obsolete
Rscript sensitivity/prune_targets.R --family=primary
```

Add `--execute` only after reviewing the printed targets and sizes. The
`obsolete` mode wraps `targets::tar_prune()`; family modes use
`targets::tar_delete()` to remove fitted models, posterior-mixture draw
objects, and exact-refit job objects while leaving the versioned bundle
untouched. Deleted targets will be recomputed if the full pipeline is later
requested.

Sensitivity-analysis pilots are disabled by default and therefore do not
compile or sample during an ordinary run. Start with the non-sampling
`sensitivity_pilot_structure_audit` target and see
[`sensitivity/README.md`](sensitivity/README.md) before enabling a pilot family.

The OpenPowerlifting reference values are read from a compact, versioned
summary of the snapshot dated 29 August 2026 and accessed on 2 September
2026. The summary records the source filename, archive SHA-256 hash, row
counts, filtering yield, and calculated reference standard deviations. A
later nightly release therefore cannot silently change the manuscript. See
[`data/external/README.md`](data/external/README.md) for full provenance and
the explicit refresh procedure.

## Render the documents

The main manuscript and all 15 HTML supplements render from the checked-in
reporting files. The rendering checkout can omit target stores, model archives,
and the original workbooks. From a clean checkout with R 4.4.3, the restored R
rendering packages, Quarto, XeLaTeX, and Python 3.12:

```sh
python -m pip install -r scripts/requirements-release.txt
python scripts/release_bundle.py outputs/preprint --render
```

The builder writes an index, PDF, HTML supplements, linked data, and a manifest
in `outputs/preprint/`, with a ZIP beside it. It checks the PDF link annotations
and the supplement filenames, fragments, and download paths. For an offline
review copy, run the builder against the rendered files with `--review` and a
new destination. Open the review bundle using `Open review.cmd` on Windows or
`python open-review.py` from the extracted bundle. The launcher serves only
that bundle over loopback and opens the index in the browser.

When approved analytical results change, restore the corresponding frozen
outputs and run `Rscript scripts/materialize_preprint_reporting.R` to refresh
the portable inputs. The main graph tracks these files as a file target. The
age/duration, training-status, and primary figure targets write their own
interfaces. Review changes in reported values and figures before committing.
The PRISMA include builder uses the registers, archived workbooks and portable
primary study effects. Check the checklist locations against the final PDF.

## Publish the supplements on GitHub Pages

In the repository, select **Settings → Pages → Build and deployment →
Source → GitHub Actions**. The **Preprint package** workflow deploys to
GitHub Pages after a push to `master` passes the package checks. Pull
requests build a reviewable package. To redeploy manually, open **Actions →
Preprint package → Run workflow**, select `master`, and enable **Publish the
checked manuscript and supplements to GitHub Pages**. Approve the
`github-pages` environment if GitHub requests approval.

The workflow publishes the checked package with its `supplement/` and
`manuscript/generated/` layout. Check the deployed site and sample manuscript
links before circulating the submission PDF.

## Final reproducibility cleanup

Before the repository is released, perform and document a fresh-clone test on
a machine or clean environment with no project-specific state. The release
check must include:

- restoring the exact R environment from `renv.lock` without relying on an
  existing user library;
- detecting required external executables such as Quarto, Rtools, XeLaTeX,
  and CmdStan by verifying that the configured files actually exist;
- avoiding hard-coded paths and stale environment variables, including a
  missing standalone Quarto installation when RStudio's bundled Quarto is
  available;
- failing early with clear installation or configuration instructions when a
  required dependency is genuinely unavailable;
- checking that targeted pipeline runs do not require unrelated document or
  analysis dependencies merely to construct the pipeline; and
- running the documented restore, manifest, representative lightweight
  targets, document rendering, and—where practical—the full pipeline from the
  clean clone.

Include the separate age/duration and training-status workflows in that
reproducibility record: record stages completed, prerequisite manifest
identities, and standalone supplement-render results. A successful main
pipeline run alone does not verify these supplementary analyses.

Record the model-count, runtime, storage, hardware, and known-limitations audit
described in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). Confirm explicitly
that the test began with no model archive or pre-existing target store.

This checklist is a release requirement rather than evidence that the current
development environment has already passed a clean-clone reproducibility
test.

## License

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
License](https://creativecommons.org/licenses/by-nc-sa/4.0/).

