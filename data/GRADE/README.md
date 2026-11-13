# GRADE decision sheet

[`rob2_overall_judgments.csv`](rob2_overall_judgments.csv) provides the final
overall risk-of-bias judgment by study number used by the GRADE analysis.
See the [risk-of-bias supplement](https://jamessteeleii.github.io/protein_dose_response_meta/supplement/risk-of-bias.html)
for the domain-level judgments.

`grade_domain_judgments.csv` is the reviewer-controlled decision layer for the
GRADE evidence profile. Automated analyses provide evidence and suggested
flags, but do not write final domain judgments.

There is one row for every retained primary outcome, resistance-training
condition, and planned contrast. For RT, the fixed 0.8-g/kg/day family contains
1.2, 1.6, 2.0, and 2.4 versus 0.8 g/kg/day. The fixed 1.6-g/kg/day RT family
begins at 2.0 and increases by 0.4 g/kg/day to the highest increment within
each outcome's pooled observed range: 4.4 for whole-body lean soft tissue mass,
2.4 for direct muscle mass/size, and 3.2 for strength and physical function.

For No RT, proposed intakes are restricted to the observed upper intake within
each outcome: 1.80 g/kg/day for whole-body lean soft tissue mass and about 1.82
g/kg/day for the other outcomes. On the planned 0.4-g/kg/day grid, this retains
only 1.2 and 1.6 versus 0.8 g/kg/day. No fixed-1.6-reference No-RT contrast is
therefore reported.

The report's strict spanning count is a count of unique studies, not arms. A
study is counted once when its within-study achieved-intake range brackets both
contrast endpoints (study minimum no greater than the lower endpoint and study
maximum no less than the upper endpoint). No tolerance is used for this count.

For a completed row:

- `starting_certainty` is an integer from 1 (very low) to 4 (high).
- Each `*_downgrade` is 0, 1, or 2.
- `other_adjustment` is normally 0; positive values downgrade and negative
  values upgrade. Any upgrade requires an explicit rationale.
- Domain labels should use `not serious`, `serious`, `very serious`, or a
  clearly explained alternative such as `unable to assess`.
- Publication bias can use `not suspected`, `strongly suspected`, or another
  explicitly justified formulation.
- `importance` should be `critical`, `important`, or `not important`.

Blank numeric fields intentionally keep overall certainty as `Pending`. This
prevents the pipeline from turning diagnostic thresholds into GRADE judgments
without reviewer sign-off.

## Decision rules used for the completed appraisal

- Start at high certainty because the evidence derives from randomised trials.
- Rate risk of bias as serious for every outcome/RT stratum because most
  contributing studies had some concerns or high RoB 2 judgments, RoB 2 was
  assessed for an index result rather than every model-derived contrast, and
  modelling achieved intake may retain residual confounding.
- Use ASI as the primary inconsistency prompt. Rate inconsistency serious at
  ASI >=25%; do not downgrade for DI alone because decision-category
  dispersion is addressed under imprecision. For whole-body mass with RT,
  2.0 versus 1.6 g/kg/day is also rated serious: overall ASI was borderline
  (24.9%) but direct-support ASI was 42.4%.
- Rate indirectness serious when either target dose lies outside the
  outcome/RT-specific observed range; otherwise do not downgrade. This avoids
  counting achieved-intake residual confounding in two domains.
- Rate imprecision against the preregistered +/-0.10-SMD SESOI: no downgrade
  when the 95% CrI remains in one decision region, one level when it crosses
  one boundary, and two levels when it crosses both boundaries.
- Rate publication bias as strongly suspected (one level) for all outcomes,
  based on the outcome-level RoBMA analyses and attenuation after adjustment.
- Apply no dose-response upgrade. The evidence already starts high, the
  updated guidance focuses the upgrade on non-randomised evidence, and GRADE
  discourages upgrading after downgrading.
- Treat all four outcomes as critical. The completed judgments are reported
  in the [GRADE supplement](https://jamessteeleii.github.io/protein_dose_response_meta/supplement/grade-appraisal.html).

## Run after the current model pipeline has finished

The GRADE branch is opt-in and does not refit any model. It reuses the accepted
primary fits, corrected LOO objects, primary population-effect mixtures, RoB 2
results, and publication-bias results. Run the final report target to build all
of its missing dependencies in sequence:

```r
Sys.setenv(
  PROTEIN_MODEL_WORKERS = "1",
  PROTEIN_UPDATE_2026_GRADE = "run"
)
targets::tar_make(update_2026_grade_supplement_html)
```

One worker is the conservative default because the inconsistency stage opens
large accepted primary fits one at a time. Do not apply the later lean-pruning
step to those fits until this extraction has completed. The finished report is
written to `supplement/grade-appraisal.html`, with its underlying tables and
checksums in `results/update_2026/grade/`.
