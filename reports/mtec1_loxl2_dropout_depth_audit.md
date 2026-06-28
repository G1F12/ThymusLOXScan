# mTEC1 Loxl2 dropout/depth audit

## Executive summary

Dropout/depth artifact risk classification: **Partial / unresolved**.

Annotated mTEC1 cells were identified as `13:mTEC1` in `cell_type_subset`. Raw counts came from `adata.raw.X`. This audit is descriptive because the mTEC1 comparison has only 2 young and 2 aged biological samples, and cell-level models do not provide independent biological replication.

The Loxl2 signal remains directionally young-higher/aged-lower in direct detection, pseudobulk CPM, pairwise sample checks, and the implemented depth-matched downsampling. However, dropout/depth artifact is not fully ruled out because the biological replicate count is small and Loxl2 detection is sparse in aged mTEC1.

## Direct detection-rate result

Young mean mTEC1 Loxl2-positive percentage: 9.0%.
Aged mean mTEC1 Loxl2-positive percentage: 0.6%.

| sample | age | mTEC1 cells | Loxl2-positive cells | detection % | Loxl2 log2(CPM+1) |
|---|---:|---:|---:|---:|---:|
| mo02_CD45neg1_d0 | 02mo | 300 | 26 | 8.7 | 3.444 |
| mo02_CD45neg2_d0 | 02mo | 279 | 26 | 9.3 | 3.484 |
| mo18_CD45neg1_d0 | 18mo | 288 | 1 | 0.3 | 0.604 |
| mo18_CD45neg2_d0 | 18mo | 763 | 7 | 0.9 | 1.296 |

## Sequencing-depth/QC result

Aged mTEC1 samples have lower mean median total_counts or detected-gene depth than young samples.
Young mean median total_counts was 8535.2; aged mean median total_counts was 5547.0. Young mean median detected genes was 2912.8; aged mean median detected genes was 2244.8.
Lower depth alone does not obviously explain the Loxl2 drop if aged depth is comparable or higher; if aged depth is lower, it remains a plausible partial contributor.

## Downsampling result

Downsampling used a fixed seed and target per-cell raw library depth of 5457, chosen as the minimum median raw library size across the four mTEC1 samples. Cells already below the target were retained at their observed depth and counted in the output tables.
All young samples remained higher than all aged samples in 100.0% of iterations for detection rate and 100.0% for log2(CPM+1).
The mean downsampled aged-minus-young difference was -0.0465 for detection rate and -2.766 for log2(CPM+1).

## Pairwise sample-level result

Both young samples are higher than both aged samples for Loxl2 detection rate, log2(CPM+1), and CPM.

## Descriptive QC-adjusted model

The descriptive cell-level logistic model was fit with QC covariates; cells are not independent biological replicates.

## External context

GSE223049 provides broad sorted thymic epithelial context in which Loxl2 is directionally aged-lower, but that bulk epithelial comparison cannot test mTEC1 specificity or single-cell dropout.
E-MTAB-8560 provides TEC and mTEC-like single-cell age-series context with aged-lower Loxl2 tendency in broad TEC, mTEC-like, mTEClo, and mTEChi groupings. Those labels are not one-to-one with the current mTEC1 annotation.
Together, these external results reduce the chance that the entire observation is a one-off internal artifact, but they do not rule out dropout or depth sensitivity in GSE240016 annotated mTEC1 cells.

## Final interpretation

Dropout/depth artifact is not fully ruled out. The signal remains directionally supported after several checks, including sample-level pairwise comparisons and depth-matched downsampling, but the evidence remains limited by sparse Loxl2 detection and n=2 young versus n=2 aged biological samples.

## Manuscript recommendation

The v5.2 manuscript should retain cautious wording: mTEC1 `Loxl2` remains a transcript-level candidate with partial/unresolved dropout/depth risk, not a protein-level or functional finding.

## Output files

- `results\tables\mtec1_loxl2_dropout_depth_audit.tsv`
- `results\tables\mtec1_loxl2_pairwise_sample_direction.tsv`
- `results\tables\mtec1_loxl2_downsampling_iterations.tsv`
- `results\tables\mtec1_loxl2_downsampling_summary.tsv`
- `results\tables\mtec1_loxl2_depth_adjusted_logistic.tsv`
