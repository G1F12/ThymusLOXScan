# v5.4 Baran-Gale and falsification update summary

## Scope

This update adds two computational analyses without modifying the manuscript, README, one-page summary, releases, or tags.

## mTEC1 Loxl2 matched-gene falsification

Status: completed.

Classification: Loxl2-like collapse uncommon among matched genes.

GSE240016 annotated mTEC1 cells were summarized at the biological-sample level using raw counts from `adata.raw.X`. Samples with fewer than 10 annotated mTEC1 cells were excluded to match the prior dropout/depth audit. Loxl2 showed young mean detection of 0.090, aged mean detection of 0.006, and aged-minus-young detection delta of -0.084.

Matched-gene background rates were descriptive:

| window | matched genes | fraction detection delta <= Loxl2 | fraction detection delta <= Loxl2 plus ordered young>aged | fraction log2CPM delta <= Loxl2 |
|---|---:|---:|---:|---:|
| strict | 2453 | 0.009 | 0.009 | 0.007 |
| medium | 4937 | 0.024 | 0.024 | 0.005 |
| broad | 7228 | 0.058 | 0.058 | 0.006 |

Interpretation: Matched-gene falsification suggested that the Loxl2 detection pattern was uncommon among similarly expressed mTEC1 genes, supporting prioritization as a transcript-level candidate while not excluding dropout/depth effects.

## Baran-Gale / MouseThymusAgeing reanalysis

Status: completed through cached MouseThymusAgeing ExperimentHub SMART-seq2 RDS files. The requested R exporter was created, but `Rscript` was not available on PATH in this environment, so the Python script used the cached public RDS fallback.

Classification: completed: Loxl2 aged-lower in TEC/mTEC-like group.

Key Loxl2 direction in TEC/mTEC-like groups: aged-lower. In the mTEC-like grouping, oldest-minus-youngest log2CPM delta was -2.168 with Spearman rho -0.565 over sample-proxy summaries. Broad TEC, mTEClo, mTEChi, and cTEC also showed aged-lower oldest-minus-youngest direction; post-Aire mTEC was flat by oldest-minus-youngest log2CPM.

Interpretation: Baran-Gale / MouseThymusAgeing provides independent TEC/mTEC-like directional support for aged-lower Loxl2 at the transcript level. This is not exact GSE240016 mTEC1 validation.

## Manuscript and release recommendation

Manuscript update: recommended only after human review of whether these analyses materially improve the current cautious v5.2 interpretation. No manuscript file was changed here.

New release: not recommended automatically. A release should be considered only after reviewing whether the new analyses warrant a versioned scientific update.
