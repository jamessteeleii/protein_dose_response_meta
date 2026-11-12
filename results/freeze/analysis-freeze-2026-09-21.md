# Analysis freeze: 21 September 2026

## Scope

This freeze captures the completed primary analyses and the completed
sensitivity analyses before data from the updated systematic search are added.
The canonical extraction workbook has not been changed for the update.

## Verification outcome

- The primary analysis contains 104 accepted models: 26 candidate models for
  each of whole-body lean mass, direct measures of muscle size, strength, and
  physical function.
- All 104 accepted primary models pass the recorded sampling-diagnostic rules.
- Final primary LOO diagnostics contain no remaining Pareto-k warnings, and all
  four outcomes are marked ready for stacking.
- The 26 stacking weights within each primary outcome sum to one.
- Versioned verification bundles exist for the primary analysis and nine
  sensitivity/reporting families. Every file listed in their manifests was
  present and matched its recorded SHA-256 checksum at the freeze.
- The strategy-formula regression test passes for all 26 arm-based model
  formulas and their within-contrast counterparts.
- All 30 R pipeline and support files parse successfully, and the pipeline
  manifest resolves to 830 targets.
- `renv::status()` reports a consistent project library and lockfile.

The machine-readable audit is in `verification-bundle-audit.csv`. Individual
artifact descriptions, sizes, and checksums are in
`results/verification/<family>/manifest.csv`.

## Preserved fit objects and pruning

Before pruning, the four accepted primary model ensembles were copied from the
targets store to `model_archive/primary/` and verified byte-for-byte. These
files total 58,233,219,439 bytes (about 54.23 GiB). The directory is excluded
from Git because of its size; `primary-fit-archive-manifest.csv` is versioned
and records the filename, size, source-target timestamp, and SHA-256 checksum
for each archive. This local archive should be included in the project's
external research-data backup.

After that verification, the guarded primary prune removed 223 redundant
cached target records. The `_targets` directory decreased from 269,423,360,522
bytes (250.92 GiB) to 45,903,521,872 bytes (42.75 GiB), releasing
223,519,838,650 bytes (208.17 GiB). Versioned results, verification bundles,
raw inputs, sensitivity checkpoints, and the separate primary fit archive were
not removed.

The sensitivity fit-retention audit records an important limitation. Compact
reporting checkpoints survive for every completed family, but full posterior
fit archives are currently present only for the residual-correlation family
(104 of 104 fits). The other completed sensitivity families have their 104
accepted fits documented but not retained locally, so they require refitting
if posterior draws are needed again.

## Pipeline state and interpretation

The reporting outputs are frozen and verified; the live target graph is not a
claim of a completely up-to-date cache. It reports 627 outdated targets because
the analysis code and reporting pipeline were developed after the cached fits,
and because the deliberately pruned primary fit targets are now absent. Running
the full pipeline would therefore initiate substantial model fitting. Future
work should consume the versioned verification bundles or the preserved primary
fit archives unless an intentional refit is planned.

This distinction is central to the update workflow: the freeze preserves the
submitted-analysis state, while new studies should be entered in the separate
incremental workbook and incorporated through a dedicated update pipeline.
