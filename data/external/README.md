# External analysis inputs

This directory contains small external inputs that are required by the
analysis pipeline. The copies are pinned to exact source revisions so a fresh
clone does not depend on mutable branch contents.

## Resistance-training meta-analysis data

`steele_rt_extracted_data_289a48f.csv` was retrieved from commit
`289a48f79da71db0c374bb64fb1df33e32b7185d` of
[Meta-Analysis-of-Variation-in-Resistance-Training](https://github.com/jamessteeleii/Meta-Analysis-of-Variation-in-Resistance-Training).

Source file:
`data/Polito et al. RT Extracted Data.csv`

SHA-256:
`2AE9EF09C7343FCB32272A2FECA942DE61C23F08D6E7B9C38E8179F4A474FCD8`

## Benito study data

`benito_study_data_18f3fc7.csv` was retrieved from commit
`18f3fc73a65d101597c3d77c8a06f7d9483fa10b` of
[glp1_rt_cits](https://github.com/jamessteeleii/glp1_rt_cits).

Source file: `data/benito_study_data.csv`

SHA-256:
`474B23E015E510FF87692E751DE28B87209609CB2BC297AB34E769358E7E4C02`

## OpenPowerlifting data

The manuscript uses the versioned summary
`openpowerlifting_strength_reference_2026-08-29.csv`. It was calculated from
the OpenPowerlifting snapshot dated 29 August 2026 and accessed from the
official nightly bulk download on 2 September 2026. The source archive
contained `openpowerlifting-2026-08-29-04805e21.csv` and had SHA-256:

`9B2AC7AEE514AB6A8454456435F4DAF49F91EFF6A2911B39B2F97CFDD4939996`

The full archive is intentionally not stored in Git because it contains more
than four million rows. The committed summary records the source identity,
row counts, filtering yield, exercise-specific sample sizes, and calculated
standard deviations needed by the manuscript. The pipeline reads this pinned
summary, so a later nightly release cannot silently change the reported
values. `download_open_powerlifting_strength_sd()` remains available as an
explicit refresh utility; updating the reference values requires creating and
reviewing a newly dated summary file.

Official download page:
[OpenPowerlifting Bulk CSV Downloads](https://openpowerlifting.gitlab.io/opl-csv/bulk-csv.html)

OpenPowerlifting releases its CSV data to the public domain. Attribution:
This project uses data from the OpenPowerlifting project,
<https://www.openpowerlifting.org/>.
