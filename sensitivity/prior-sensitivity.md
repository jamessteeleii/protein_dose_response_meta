# Prior sensitivity of the primary relative-intake candidates

Run these commands from the project root with the restored R 4.4.3 environment
and Quarto on the path:

```sh
Rscript sensitivity/run_prior_sensitivity.R
quarto render supplement/prior-sensitivity.qmd
```

The runner uses four independent workers by default and writes its cache to
`_targets_prior_sensitivity/`. Set `PROTEIN_MODEL_WORKERS` to override the
worker count. Each candidate branch reads one accepted fit and returns compact
summaries. The runner uses the archived posterior draws without refitting
models or invoking the main analysis graph.

Restore the frozen `results/verification/primary_update_2026/` bundle and
accepted relative-primary fits under `model_archive/update_2026/primary/`
from the research archive. Contact the corresponding author for access.
The workflow verifies the frozen bundle manifest and inputs against the
manifest index, then reads all 104 accepted primary candidates, 26 per
outcome. Missing accepted fits produce an error naming the exact path
expected in this clone.

The workflow applies power-scaling to the saved Stan `lprior` terms and,
separately, to the empirical intercept and resistance-training priors. It
compares powers 0.8 and 1.2 with the original posterior on the same outcome
population prediction profile, for seven contrasts within each observed
outcome and training-condition range. Failed importance-sampling checks stay
in the diagnostics and have missing reweighted estimates. The report compares
the priors within each candidate and lists the original stacking weights.

The Quarto report can render in a fresh clone once the five small files in
`manuscript/generated/prior-sensitivity/` are available. Review the generated
diagnostics, contrasts, support table, candidate table and figure before
including them in a handoff. The report uses only those portable results.
