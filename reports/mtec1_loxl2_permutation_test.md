# mTEC1 Loxl2 exact permutation test

## Scope

This audit uses the four GSE240016 biological samples with at least 10 annotated `13:mTEC1` cells. Loxl2 values were recomputed from `adata.raw.X` in `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad`.

With n=2 versus n=2 biological samples, the exact permutation test has only six possible labelings, so p-value resolution is limited. This test is used to calibrate the parametric DESeq2 result rather than replace it.

## Sample values

| sample | age | cells | detecting cells | detection rate | log2(CPM+1) |
|---|---:|---:|---:|---:|---:|
| mo02_CD45neg1_d0 | 02mo | 300 | 26 | 0.086667 | 3.443872 |
| mo02_CD45neg2_d0 | 02mo | 279 | 26 | 0.093190 | 3.483679 |
| mo18_CD45neg1_d0 | 18mo | 289 | 1 | 0.003460 | 0.596223 |
| mo18_CD45neg2_d0 | 18mo | 764 | 7 | 0.009162 | 1.289430 |

## Exact permutation results

- detection_rate: observed aged-minus-young statistic = -0.083617.
- detection_rate: two-sided exact p-value = 0.333333.
- detection_rate: one-sided exact p-value for aged-lower direction = 0.166667.
- log2CPM1: observed aged-minus-young statistic = -2.520949.
- log2CPM1: two-sided exact p-value = 0.333333.
- log2CPM1: one-sided exact p-value for aged-lower direction = 0.166667.

The minimum possible one-sided p-value is 1/6 = 0.166667; the minimum possible two-sided p-value is 2/6 = 0.333333 when both directions are counted.

## Interpretation

Both young samples remain higher than both aged samples for detection rate and log2(CPM+1). The exact test supports the direction as sample-level ordered, but the p-values cannot be small with only four samples.

## Output

- `results/tables/mtec1_loxl2_permutation_test.tsv`