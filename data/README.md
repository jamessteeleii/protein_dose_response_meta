# Data and extraction workbooks

| Excel workbook | Scope | Download |
| --- | --- | --- |
| [`protein_meta_data_all.xlsx`](protein_meta_data_all.xlsx) | Canonical extraction, 109 study records | [Download canonical workbook](https://raw.githubusercontent.com/jamessteeleii/protein_dose_response_meta/master/data/protein_meta_data_all.xlsx) |
| [`protein_meta_data_update_2026.xlsx`](update_2026/protein_meta_data_update_2026.xlsx) | 2026 search update, 10 study records | [Download 2026 workbook](https://raw.githubusercontent.com/jamessteeleii/protein_dose_response_meta/master/data/update_2026/protein_meta_data_update_2026.xlsx) |

Open either file in Excel or LibreOffice. Each file has four sheets.

| Sheet | Records |
| --- | --- |
| Study characteristics | Study identifiers, bibliographic details and study-level information |
| Arm characteristics | Intervention, comparator, participant characteristics and protein intake by arm |
| Arm based data by timepoint | Outcome measurements by study, arm, and timepoint |
| Contrasts between timepoints | Reported changes and their uncertainty |

Join the sheets by `Study No` and, where available, `Arm No`. Arm numbers are
within-study identifiers. Use both workbooks and the
[2026 review decisions](update_2026/review_decisions.csv) for the complete
updated review. The pipeline combines the workbooks in memory and applies
eligibility and analysis rules. Blank extraction cells indicate unavailable
or inapplicable entries. Check the column headings for outcome and dose units.
The [2026 data guide](update_2026/README.md) describes update-specific decisions.

## Analysis-ready data and supporting files

- [Updated primary study-level effects](../manuscript/generated/reporting/primary-study-effects.csv.gz)
  provides the effect sizes, sampling variances, and selected input fields
  used in the updated primary models. Extract the gzip file to open the CSV.
- [Systematic-review audit downloads](https://jamessteeleii.github.io/protein_dose_response_meta/supplement/systematic-review-data-audit.html#downloadable-audit-files)
  includes the summary, contrast flow, and reasons for exclusion.
- [Digitised figure data](graphical_data/) and
  [external reference inputs](external/README.md) document the supporting data.
- [Strategy classifications](protein_manipulation_strategy_audit.csv)
  records the historical study-level coding and rationale.
- [Risk-of-bias and GRADE inputs](GRADE/README.md) document the reviewed judgments.

The [frozen reporting files](../manuscript/generated/reporting/) contain
the results used in the manuscript and supplements. They can be inspected
without running the analysis. See the [results and rendering guide](../README.md#inspect-results-without-refitting).

The public canonical workbook preserves all cell values and formulas,
including study author emails. Its working comments and local Excel metadata
were removed. The 2026 workbook matches the frozen original byte for byte.
Use the [input manifest](../config/data-manifest.yml) to identify the downloaded
versions and check their SHA-256 hashes. For the staged computational workflow,
see the [clean-clone instructions](../REPRODUCIBILITY.md#start-from-a-clean-clone).
