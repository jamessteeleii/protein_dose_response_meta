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

The full OpenPowerlifting CSV is intentionally not stored in this repository.
The `open_powerlifting_url` target tracks the official nightly ZIP using its
HTTP ETag/Last-Modified metadata. When the snapshot changes, the downstream
target downloads it to a temporary directory, reads only the columns needed
to estimate the reference strength standard deviations, records the archive
filename and SHA-256 hash, and removes the temporary download.

The values are used as indicative reference values rather than as estimates
whose exact replication requires a permanently fixed OpenPowerlifting
snapshot. Accordingly, they may change slightly when the nightly dataset is
updated.

Official download page:
[OpenPowerlifting Bulk CSV Downloads](https://openpowerlifting.gitlab.io/opl-csv/bulk-csv.html)

OpenPowerlifting releases its CSV data to the public domain. Attribution:
This project uses data from the OpenPowerlifting project,
<https://www.openpowerlifting.org/>.
