# Hypothesis Falsification Report

Audit date: 2026-06-26

Hypothesis under attack: LOX-family transcript changes reflect genuine age-associated biology in thymic stromal cells.

Adversarial alternative: the observed LOX-family changes are artifacts of cell composition, pseudobulk aggregation, normalization, sequencing depth, QC filtering, dropout, batch/sample imbalance, annotation, outliers, or gene detection.

Bottom line: the main highlighted directions do **not** disappear under the controls that can be run from this repository. However, the controls do **not** prove the biological hypothesis. They show that the signals are heavily entangled with detection rate, sample depth, sample type, and cell composition. The broad fibroblast LOX-family result is especially vulnerable to composition artifacts. The mTEC1 `Loxl2` result survives simple controls but remains non-confirmatory because it is only n=2 vs n=2 and mostly a near-absence/detection-rate result in aged cells.

## Controls Run

Inputs inspected:

- `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad`
- `data/processed/thymus_annotated.h5ad`
- `results/tables/lox_pseudobulk_complete_results.csv`
- `results/figures/per_sample/lox_per_sample_pseudobulk_values.csv`
- existing reports and scripts under `reports/` and `scripts/`

Controls performed:

- sample-by-stage count audit;
- cell-type and subtype composition audit;
- per-sample raw-count CPM summaries for highlighted LOX findings;
- per-sample gene detection-rate summaries;
- comparison of raw-file and processed-file normalized `adata.X` summaries;
- mitochondrial and gene-count QC sensitivity checks;
- stricter QC sensitivity checks;
- exclusion of sorted/enriched `FB` and `EC` sample types, retaining CD45-negative samples only;
- leave-one-sample-out direction checks;
- sample-level sequencing depth and QC summary;
- subtype-composition decomposition of broad fibroblast signals;
- low-expression and dropout review.

## Key Findings Under Attack

Highlighted findings from the manuscript/result tables:

| Finding | DESeq2 result | Biological samples |
|---|---:|---:|
| pooled FB `Lox` | log2FC=-1.29, padj=4.36e-5 | 3 young vs 3 aged |
| pooled FB `Loxl1` | log2FC=-0.62, padj=0.014 | 3 young vs 3 aged |
| pooled FB `Loxl2` | log2FC=-0.52, padj=0.021 | 3 young vs 3 aged |
| capsFB `Lox` | log2FC=-1.46, padj=7.52e-5 | 3 young vs 3 aged |
| medFB `Loxl1` | log2FC=+0.75, padj=0.023 | 3 young vs 3 aged |
| medFB `Loxl2` | log2FC=-1.00, padj=0.0067 | 3 young vs 3 aged |
| mTEC1 `Loxl2` | log2FC=-3.29, padj=9.56e-4 | 2 young vs 2 aged |

## Sample Imbalance And Batch-Like Structure

### Result

Strong artifact risk remains. Age is confounded with sample identity and likely with sample preparation context. The dataset contains distinct sample types such as CD45-negative, FB-enriched, and EC-enriched samples, and these are not modeled.

### Evidence

Sample counts:

| sample | stage | cells |
|---|---:|---:|
| `mo02_CD45neg1_d0` | 02mo | 3,381 |
| `mo02_CD45neg2_d0` | 02mo | 3,015 |
| `mo02_EC_d0` | 02mo | 572 |
| `mo02_FB_d0` | 02mo | 3,994 |
| `mo18_CD45neg1_d0` | 18mo | 1,997 |
| `mo18_CD45neg2_d0` | 18mo | 4,289 |
| `mo18_EC_d0` | 18mo | 710 |
| `mo18_FB_d0` | 18mo | 4,974 |

Median total counts differ strongly by sample and age context:

| sample | stage | median total counts | median genes |
|---|---:|---:|---:|
| `mo02_CD45neg1_d0` | 02mo | 5,474 | 2,182 |
| `mo02_CD45neg2_d0` | 02mo | 5,945 | 2,336 |
| `mo02_FB_d0` | 02mo | 7,572 | 2,659 |
| `mo18_CD45neg1_d0` | 18mo | 4,857 | 2,113 |
| `mo18_CD45neg2_d0` | 18mo | 4,303 | 1,907 |
| `mo18_FB_d0` | 18mo | 3,069 | 1,527 |

The aged FB-enriched sample has much lower median depth and gene detection than the young FB-enriched sample.

### Interpretation

This is a plausible artifact mechanism. Lower depth and lower detected genes in aged samples can reduce low/moderate-expression gene detection and create apparent downregulation.

### Did the LOX signal survive?

Partly. The highlighted directions survived exclusion of FB/EC-enriched samples, but this leaves only 2 vs 2 CD45-negative samples for many comparisons. Surviving with n=2 vs n=2 is not robust validation.

### Fix

Model batch/sample-preparation covariates if available. If age and batch/preparation are inseparable, state that causal age inference is not possible. Replicate in an independently balanced dataset.

## Sequencing Depth And Normalization Artifacts

### Result

The main directions are not solely artifacts of DESeq2 or pseudobulk CPM normalization, because raw-count CPM, normalized `adata.X`, and processed-object summaries agree in direction. But the signal is still depth-sensitive because LOX detection is much lower in aged samples.

### Evidence

Per-sample pseudobulk `log2(CPM+1)` deltas, aged minus young:

| Finding | delta log2(CPM+1) | detection-rate delta | Direction survives CPM? |
|---|---:|---:|---|
| pooled FB `Lox` | -1.190 | -0.192 | yes |
| pooled FB `Loxl1` | -0.413 | -0.197 | yes |
| pooled FB `Loxl2` | -0.317 | -0.130 | yes |
| capsFB `Lox` | -1.494 | -0.223 | yes |
| medFB `Loxl1` | +0.821 | +0.048 | yes |
| medFB `Loxl2` | -0.907 | -0.129 | yes |
| mTEC1 `Loxl2` | -2.521 | -0.084 | yes |

Normalized `adata.X` deltas in the raw file and processed file were nearly identical:

| Finding | normalized mean delta, raw file | normalized mean delta, processed file |
|---|---:|---:|
| pooled FB `Lox` | -0.2650 | -0.2648 |
| pooled FB `Loxl1` | -0.3921 | -0.3923 |
| pooled FB `Loxl2` | -0.1283 | -0.1281 |
| capsFB `Lox` | -0.3067 | -0.3067 |
| medFB `Loxl1` | +0.1174 | +0.1174 |
| medFB `Loxl2` | -0.1464 | -0.1464 |
| mTEC1 `Loxl2` | -0.0678 | -0.0678 |

### Interpretation

The result is not created by one normalization representation. However, lower detection in aged cells remains a serious artifact pathway. For several genes, especially `Lox`, `Loxl2`, and mTEC1 `Loxl2`, expression and detection collapse together.

### Did the LOX signal survive?

Yes for direction, across CPM and normalized expression summaries.

### Fix

Report sample-level detection rates beside expression. Use sample-level detection models as sensitivity analyses. Validate with methods insensitive to scRNA-seq dropout, such as RNA FISH, spatial transcriptomics, or targeted qPCR on sorted populations.

## Mitochondrial Filtering And QC Thresholds

### Result

The highlighted directions survive the repository QC threshold and a stricter QC threshold.

### Evidence

Repository-style QC filter: `n_genes_by_counts >= 200`, `n_genes_by_counts <= 6000`, `mito_frac <= 0.20`. This removed 55 cells.

Strict QC filter tested: `n_genes_by_counts >= 500`, `n_genes_by_counts <= 5000`, `mito_frac <= 0.10`. This removed 451 cells.

Under both filters, all highlighted directions remained the same:

| Finding | original delta log2(CPM+1) | repository QC delta | strict QC delta |
|---|---:|---:|---:|
| pooled FB `Lox` | -1.190 | -1.189 | -1.187 |
| pooled FB `Loxl1` | -0.413 | -0.413 | -0.413 |
| pooled FB `Loxl2` | -0.317 | -0.316 | -0.315 |
| capsFB `Lox` | -1.494 | -1.494 | -1.493 |
| medFB `Loxl1` | +0.821 | +0.821 | +0.833 |
| medFB `Loxl2` | -0.907 | -0.907 | -0.904 |
| mTEC1 `Loxl2` | -2.521 | -2.514 | -2.537 |

### Interpretation

Mitochondrial or simple gene-count filtering does not explain the highlighted directions.

### Did the LOX signal survive?

Yes.

### Fix

Still include this sensitivity table in the supplement, because the primary pseudobulk script currently reads the raw annotated object rather than the filtered processed object.

## Cell Composition Effects

### Result

This is one of the strongest artifact explanations for the broad fibroblast findings.

### Evidence

Fibroblast subtype fractions changed by age:

| stage | capsFB fraction | intFB fraction | medFB fraction | Fat fraction |
|---|---:|---:|---:|---:|
| 02mo | 0.376 | 0.390 | 0.230 | 0.004 |
| 18mo | 0.307 | 0.247 | 0.432 | 0.014 |

The aged fibroblast compartment is more medFB-heavy and less intFB-heavy. This matters because LOX-family expression differs strongly by subtype.

Approximate decomposition of broad fibroblast shifts:

| Gene | within-subtype expression shift at young composition | composition-only shift using young expression | combined approximate shift |
|---|---:|---:|---:|
| `Lox` | -0.517 | -0.803 | -0.936 |
| `Loxl1` | +0.001 | -0.640 | -0.421 |
| `Loxl2` | -0.317 | -0.134 | -0.578 |

### Interpretation

The broad fibroblast `Lox` and `Loxl1` decreases can be substantially explained by subtype composition. For `Loxl1`, the decomposition is especially damaging: the within-subtype shift is essentially zero at young composition, while composition alone predicts a broad decrease. Yet medFB `Loxl1` increases with age. A broad pooled fibroblast statement therefore hides opposing subtype effects.

For `Loxl2`, both composition and within-subtype shifts contribute, so composition alone does not fully explain the signal.

### Did the LOX signal survive?

Subtype-specific highlighted signals survived. Broad fibroblast signals did **not** survive as clean biological interpretations because composition is a major alternative explanation.

### Fix

Do not use broad fibroblast pseudobulk as primary biological evidence. Make subtype-level analysis primary and explicitly model/sample-report composition. For broad fibroblast claims, use composition-adjusted models or stratified summaries.

## Pseudobulk Aggregation Artifacts

### Result

The key directions survive alternate aggregation summaries, but pseudobulk aggregation still hides detection-rate and composition effects.

### Evidence

The same directions appear in:

- DESeq2 pseudobulk log2FC;
- sample-level `log2(CPM+1)` from raw counts;
- sample-level normalized `adata.X` mean expression;
- per-cell detection fractions.

Examples:

- capsFB `Lox`: all aged samples have lower `log2(CPM+1)` than all young samples.
- medFB `Loxl2`: all aged samples have lower `log2(CPM+1)` than all young samples.
- medFB `Loxl1`: aged samples are higher on average, though the ranges overlap.
- mTEC1 `Loxl2`: both aged samples are lower than both young samples.

### Interpretation

The directions are not obvious artifacts of one aggregation formula. The more serious issue is that pseudobulk sums cannot distinguish fewer expressing cells from lower expression among expressing cells, and they cannot separate age from sample preparation with this small n.

### Did the LOX signal survive?

Yes for direction. No for strong causal interpretation.

### Fix

Report sample-level CPM, detection rate, and nonzero expression separately. Add mixed/hurdle sensitivity models only if enough samples exist; here, the sample count is too low for many groups.

## Dropout And Gene Detection Rates

### Result

Dropout/detection is a major artifact threat. Several highlighted findings are primarily detection-rate differences.

### Evidence

Average detection-rate changes, aged minus young:

| Finding | detection-rate change |
|---|---:|
| pooled FB `Lox` | -0.192 |
| pooled FB `Loxl1` | -0.197 |
| pooled FB `Loxl2` | -0.130 |
| capsFB `Lox` | -0.223 |
| medFB `Loxl1` | +0.048 |
| medFB `Loxl2` | -0.129 |
| mTEC1 `Loxl2` | -0.084 |

mTEC1 `Loxl2` is particularly sparse:

| sample | stage | mTEC1 cells | `Loxl2` raw counts | detection |
|---|---:|---:|---:|---:|
| `mo02_CD45neg1_d0` | 02mo | 300 | 27 | 0.0867 |
| `mo02_CD45neg2_d0` | 02mo | 279 | 28 | 0.0932 |
| `mo18_CD45neg1_d0` | 18mo | 289 | 1 | 0.0035 |
| `mo18_CD45neg2_d0` | 18mo | 764 | 7 | 0.0092 |

### Interpretation

The mTEC1 `Loxl2` signal is effectively a detection/near-absence result in aged samples, not a precise continuous expression estimate. That may be real biology, but it is also exactly where scRNA-seq dropout, low RNA content, annotation shifts, and library-depth differences are most dangerous.

### Did the LOX signal survive?

Yes for direction. No for confidence that the signal is not dropout-sensitive.

### Fix

Validate `Loxl2` detection in mTEC1 by orthogonal methods. At minimum, report sample-level detection rates and raw counts as the main evidence, not only DESeq2 fold changes.

## Low-Expression Genes

### Result

The most dramatic epithelial result is low-count and unstable by design.

### Evidence

`13:mTEC1` `Loxl2` has DESeq2 `baseMean=16.42`, and the raw sample counts are only 27, 28, 1, and 7. The result has log2FC=-3.29, but the absolute counts are tiny.

Several other nominal LOX results have very low baseMean and do not survive FDR correction, for example `Loxl2` in `7:vSMC/PC` with baseMean 2.45 and padj approximately 1.0.

### Interpretation

Large fold changes for low-count genes are easy to overstate. The mTEC1 result should be described as low detection in aged mTEC1, not as a robust quantitative expression decrease.

### Did the LOX signal survive?

Direction survives, quantitative effect size is fragile.

### Fix

Use minimum count/detection thresholds for claims. Report raw counts and detection rates. Use shrunken LFCs or avoid ranking by unshrunken LFC for low-count genes.

## Outlier Samples

### Result

The highlighted directions are not driven by a single obvious sample in simple leave-one-sample-out direction checks.

### Evidence

Leave-one-sample-out aged-minus-young `log2(CPM+1)` deltas:

| Finding | Leave-one-sample-out direction |
|---|---|
| pooled FB `Lox` | always negative; range -1.49 to -1.00 |
| pooled FB `Loxl1` | always negative; range -0.62 to -0.29 |
| pooled FB `Loxl2` | always negative; range -0.43 to -0.25 |
| capsFB `Lox` | always negative; range -1.81 to -1.18 |
| medFB `Loxl1` | always positive; range +0.55 to +0.97 |
| medFB `Loxl2` | always negative; range -1.03 to -0.77 |
| mTEC1 `Loxl2` | always negative; range -2.87 to -2.17 |

### Interpretation

No single obvious sample explains the highlighted directions. However, leave-one-out on n=2 or n=3 per group is a weak control; it cannot estimate variance reliably.

### Did the LOX signal survive?

Yes for direction.

### Fix

Obtain more biological samples. Do not use leave-one-out as proof of robustness in n=2/n=3 designs.

## Annotation Mistakes

### Result

Annotation remains a serious unresolved artifact risk. This repository does not independently rederive cell labels or test label uncertainty.

### Evidence

The workflow retains published `cell_type` and `cell_type_subset` labels. The analyses depend heavily on fine labels such as `3:capsFB`, `5:medFB`, and `13:mTEC1`. Age changes subtype proportions, and some aged samples have very different subtype distributions.

Specific vulnerability:

- mTEC1 `Loxl2` is based on fine epithelial annotation and only two samples per age.
- Broad fibroblast findings change when subtype composition changes.
- Aged mTEC1 has very low `Loxl2` detection; mislabeling of a small epithelial state could create an apparent subtype-specific absence.

### Interpretation

The controls cannot rule out annotation artifacts. This is a major remaining failure mode.

### Did the LOX signal survive?

Survives within the existing labels. It has not been tested against alternative labels.

### Fix

Repeat the analysis under alternative annotation resolutions, marker-gated subsets, and original-publication annotation confidence scores if available. Show marker distributions by sample and age for mTEC1 and fibroblast subtypes.

## Batch Effects

### Result

Batch effects cannot be ruled out from this repository.

### Evidence

The DESeq2 design is `~ stage`. No batch term is modeled. Sample names indicate different sorted/enriched sample types (`CD45neg`, `EC`, `FB`) and age-specific sample identifiers. Sequencing depth and gene detection differ by age/sample.

### Interpretation

Batch/sample-preparation confounding remains a plausible artifact explanation. The controls show direction stability, but they do not separate age from technical structure.

### Did the LOX signal survive?

Survives simple stratification excluding FB/EC enriched samples, but then many comparisons are only 2 vs 2.

### Fix

Use balanced samples and model batch/preparation. If batch is not estimable, state that batch-confounded age effects cannot be excluded.

## Summary By Artifact Mechanism

| Artifact mechanism | Did it explain away the highlighted directions? | Remaining risk |
|---|---|---|
| cell composition | partially, especially broad FB `Lox` and `Loxl1` | high |
| pseudobulk aggregation | no, directions survive CPM and normalized summaries | medium |
| normalization artifact | no, raw-file and processed-file normalized summaries agree | low-medium |
| sequencing depth | not fully, but strongly confounded with age/sample | high |
| mitochondrial filtering | no, directions survive tested QC filters | low |
| low expression | does not erase direction, but weakens mTEC1 `Loxl2` strongly | high |
| dropout/detection | does not erase direction, but may be the dominant observed phenomenon | high |
| batch effects | cannot be ruled out | high |
| sample imbalance | cannot be fixed with current data | high |
| annotation mistakes | cannot be ruled out | high |
| outlier samples | no single obvious outlier | medium because n is tiny |
| gene detection rates | direction survives, but detection shifts are large | high |

## Final Falsification Verdict

The observed LOX-family changes are **not proven to be artifacts** by the controls available in this repository. The highlighted directions survive:

- sample-level CPM checks;
- normalized-expression checks;
- repository and strict QC filters;
- exclusion of FB/EC-enriched samples;
- leave-one-sample-out direction checks;
- raw-file versus processed-file comparisons.

But the hypothesis is **not proven biologically true**. The strongest surviving alternative explanations are:

1. broad fibroblast effects are partly or mostly cell-composition artifacts;
2. aged samples have lower sequencing depth/gene detection in relevant contexts;
3. many LOX results are detection/dropout phenomena rather than stable quantitative expression shifts;
4. batch/sample-preparation confounding is not modeled;
5. mTEC1 `Loxl2` is a low-count n=2 vs n=2 result and cannot support strong inference;
6. annotation uncertainty remains untested.

The honest conclusion is narrow: within the existing annotations of this processed public dataset, the sample-level direction of several LOX-family changes is stable to simple technical controls. The result does not yet survive the controls that matter most for causal biological inference: independent validation, balanced replication, explicit batch modeling, and annotation-uncertainty sensitivity.
