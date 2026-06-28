# Fibroblast subtype-composition adjustment for LOX-family age effects

## Scope

This is a small-sample sensitivity analysis using sample-level broad fibroblast pseudobulk expression. The adjusted model regresses log2(CPM + 1) on age and the capsFB, intFB, medFB, and Fat fractions within broad FB. With only six broad-FB samples, adjusted coefficients are descriptive rather than definitive inference.

## Overall conclusion

The directional broad FB conclusions can be retained, but their inferential strength should be downgraded to a small-n composition-sensitivity result rather than treated as definitive broad-FB-intrinsic effects.

Because the adjusted model has only one residual degree of freedom after including age plus four subtype fractions, p-values are not interpreted for the adjusted fits. The useful signal is whether the age coefficient keeps its sign compared with the unadjusted sample-level broad FB effect and the subtype-stratified DESeq2 effects.

## Sample subtype composition

| sample | stage | FB cells | capsFB | intFB | medFB | Fat |
|---|---:|---:|---:|---:|---:|---:|
| mo02_CD45neg1_d0 | 02mo | 1709 | 0.401 | 0.406 | 0.188 | 0.005 |
| mo02_CD45neg2_d0 | 02mo | 1473 | 0.392 | 0.416 | 0.188 | 0.004 |
| mo02_FB_d0 | 02mo | 3977 | 0.336 | 0.348 | 0.312 | 0.004 |
| mo18_CD45neg1_d0 | 18mo | 431 | 0.399 | 0.181 | 0.401 | 0.019 |
| mo18_CD45neg2_d0 | 18mo | 775 | 0.354 | 0.227 | 0.404 | 0.015 |
| mo18_FB_d0 | 18mo | 4875 | 0.169 | 0.334 | 0.491 | 0.007 |

## Gene-level comparison

### Lox

- Broad FB DESeq2 effect: log2FC=-1.29, padj=4.36e-05, down_in_aged.
- Sample-level unadjusted age coefficient: -1.19 log2(CPM+1).
- Composition-adjusted age coefficient: -3.88 log2(CPM+1); extremely small-n OLS sensitivity model; df_resid=1.
- Subtype-stratified effects: 3:capsFB: log2FC=-1.46, padj=7.52e-05, down_in_aged; 4:intFB: log2FC=-0.626, padj=0.0716, down_in_aged; 5:medFB: log2FC=1.1, padj=0.0747, up_in_aged; 9:Fat: log2FC=-0.384, padj=1, down_in_aged.
- Conclusion: retain: adjusted coefficient keeps the broad FB direction.

### Loxl1

- Broad FB DESeq2 effect: log2FC=-0.616, padj=0.0144, down_in_aged.
- Sample-level unadjusted age coefficient: -0.413 log2(CPM+1).
- Composition-adjusted age coefficient: -3.34 log2(CPM+1); extremely small-n OLS sensitivity model; df_resid=1.
- Subtype-stratified effects: 3:capsFB: log2FC=-0.254, padj=0.482, down_in_aged; 4:intFB: log2FC=-0.512, padj=0.0541, down_in_aged; 5:medFB: log2FC=0.751, padj=0.0226, up_in_aged; 9:Fat: log2FC=-0.586, padj=1, down_in_aged.
- Conclusion: retain: adjusted coefficient keeps the broad FB direction.

### Loxl2

- Broad FB DESeq2 effect: log2FC=-0.517, padj=0.0208, down_in_aged.
- Sample-level unadjusted age coefficient: -0.317 log2(CPM+1).
- Composition-adjusted age coefficient: -2.83 log2(CPM+1); extremely small-n OLS sensitivity model; df_resid=1.
- Subtype-stratified effects: 3:capsFB: log2FC=0.147, padj=0.866, up_in_aged; 4:intFB: log2FC=-0.578, padj=0.19, down_in_aged; 5:medFB: log2FC=-0.995, padj=0.00672, down_in_aged; 9:Fat: log2FC=-0.0326, padj=1, down_in_aged.
- Conclusion: retain: adjusted coefficient keeps the broad FB direction.

## Interpretation

This analysis should be treated as a sensitivity check. The broad FB DESeq2 model remains the primary broad-cell-type result, but the adjusted coefficients test whether those signals plausibly survive subtype-mixture imbalance. Any downgraded gene should be described as composition-sensitive rather than a clean broad-FB-intrinsic age effect.