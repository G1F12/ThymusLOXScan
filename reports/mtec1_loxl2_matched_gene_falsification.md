# mTEC1 Loxl2 matched-gene falsification

## Scope

This is a matched-gene descriptive falsification analysis for GSE240016 annotated mTEC1 cells. It asks how often genes with broadly similar young mTEC1 detection and expression show a Loxl2-like young-high / aged-low collapse. It is not proof that the Loxl2 pattern is biological, and it does not exclude dropout or depth effects.

## Inputs and metadata

- mTEC1 label: `13:mTEC1` from `cell_type_subset`
- age/stage column: `age_group`
- biological sample column: `sample`
- count matrix: `adata.raw.X`
- LOX-family genes, mitochondrial genes, ribosomal genes, and genes with zero young detection were excluded from the matched background.
- Ribosomal genes were excluded from the background.
- Samples with fewer than 10 annotated mTEC1 cells were excluded to match the prior dropout/depth audit design.

| sample | age | mTEC1 cells |
|---|---:|---:|
| mo02_CD45neg1_d0 | 02mo | 300 |
| mo02_CD45neg2_d0 | 02mo | 279 |
| mo18_CD45neg1_d0 | 18mo | 288 |
| mo18_CD45neg2_d0 | 18mo | 763 |

Excluded low-cell samples:

| sample | age | mTEC1 cells |
|---|---:|---:|
| mo18_EC_d0 | 18mo | 3 |
| mo18_FB_d0 | 18mo | 2 |

## Loxl2 observed statistics

| statistic | value |
|---|---:|
| young mean detection | 0.090 |
| aged mean detection | 0.006 |
| detection delta, aged minus young | -0.084 |
| detection ratio, aged / young | 0.070 |
| young mean log2(CPM + 1) | 3.464 |
| aged mean log2(CPM + 1) | 0.950 |
| log2(CPM + 1) delta, aged minus young | -2.514 |
| all young samples above all aged samples, detection | yes |
| all young samples above all aged samples, log2CPM | yes |

## Matched-gene background

Unique genes appearing in at least one matched set: 7228.

| window | matched genes | frac detection delta <= Loxl2 | empirical p | frac delta + ordered | frac log2CPM delta <= Loxl2 |
|---|---:|---:|---:|---:|---:|
| strict | 2453 | 0.009 | 0.009 | 0.009 | 0.007 |
| medium | 4937 | 0.024 | 0.024 | 0.024 | 0.005 |
| broad | 7228 | 0.058 | 0.059 | 0.058 | 0.006 |

## Classification

Classification: **Loxl2-like collapse uncommon among matched genes**.

Matched-gene falsification suggested that the Loxl2 detection pattern was uncommon among similarly expressed mTEC1 genes, supporting prioritization as a transcript-level candidate while not excluding dropout/depth effects.

## Output files

- `results/tables/mtec1_loxl2_matched_gene_falsification_all_genes.tsv`
- `results/tables/mtec1_loxl2_matched_gene_falsification_matched_genes.tsv`
- `results/tables/mtec1_loxl2_matched_gene_falsification_summary.tsv`
- `results/figures/falsification/mtec1_loxl2_matched_gene_detection_delta.png`
- `results/figures/falsification/mtec1_loxl2_matched_gene_detection_delta.pdf`
- `results/figures/falsification/mtec1_loxl2_matched_gene_log2cpm_delta.png`
- `results/figures/falsification/mtec1_loxl2_matched_gene_log2cpm_delta.pdf`