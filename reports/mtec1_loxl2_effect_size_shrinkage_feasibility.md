# mTEC1 Loxl2 conservative effect-size audit

## Summary

DESeq2 LFC shrinkage was not run in this environment. `Rscript` was not available, and no apeglm/ashr-based shrinkage result was generated. The audit therefore reports conservative sample-level and downsampled effect-size summaries alongside the raw DESeq2 LFC.

## Effect sizes

- Raw pseudobulk DESeq2 log2FC, 18mo vs 02mo: -3.291
- Sample-level mean log2(CPM+1) delta, aged minus young: -2.521
- Median young log2(CPM+1): 3.464
- Median aged log2(CPM+1): 0.943
- Detection-rate delta, aged minus young: -0.0836
- Prior depth-matched downsampled log2(CPM+1) delta, aged minus young: -2.766

## Interpretation

The direction remains aged-lower across raw pseudobulk, direct sample-level log2(CPM+1), detection rate, and prior downsampling. However, shrinkage was unavailable and the biological design remains n=2 vs n=2, so the raw DESeq2 LFC should be treated as a prioritization effect size rather than a stable quantitative estimate.

Recommended reporting: show the raw DESeq2 LFC together with the sample-level log2(CPM+1) delta, detection delta, and downsampled delta. State that shrinkage was not run in this environment.

## Output

- `results/tables/mtec1_loxl2_effect_size_conservative_summary.tsv`
