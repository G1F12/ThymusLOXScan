# Highest-impact next analysis after falsification

## Decision

Chosen analysis: external validation in an independent public thymic aging dataset.

Dataset used: GSE223049, a mouse sorted bulk RNA-seq dataset with young 2-month and aged 22-24-month samples across thymic and immune cell types.

Script run:

`python scripts/external_validation_gse223049_lox.py`

Primary outputs:

- `results/tables/external_gse223049_lox_validation.tsv`
- `results/tables/external_gse223049_lox_validation_summary.tsv`

## Why this was the highest-impact choice

The falsification report concluded that the internal LOX-family directions survived many simple controls, but that the study still did not survive the controls most important for biological confidence: independent validation, balanced replication, explicit batch modeling, and annotation-uncertainty sensitivity.

Several strong internal sensitivity analyses have now been run:

- fibroblast subtype-composition adjustment;
- detection-rate versus nonzero-expression decomposition;
- annotation marker sanity checks and stricter marker-positive subsets.

Those analyses increase confidence that the current dataset was handled cautiously, but they remain constrained by the same original samples, same sample-preparation structure, same published annotations, and same scRNA-seq detection limitations. An independent dataset is therefore the single analysis that most increases confidence because it asks whether a related LOX-family aging pattern appears outside the original GSE240016 single-cell dataset.

## Why the other candidates were not chosen

Fibroblast subtype-composition adjustment was important for claim safety, but it still used only six broad fibroblast samples from the original dataset and could not definitively separate age from subtype mixture.

Detection-rate decomposition clarified whether signals reflect fewer expressing cells or lower nonzero expression, but it did not address independent replication.

Annotation uncertainty checks helped test whether labels were marker-consistent, but they still depended on the same original annotated object and could not prove annotation correctness.

External validation directly addresses the falsification report's strongest remaining weakness: whether the signal generalizes beyond the current dataset.

## Analysis Performed

The script downloaded the processed count matrix:

`GSE223049_RNA_seq_counts_23_cell_types.txt.gz`

from GEO and summarized LOX-family expression for the sorted thymic populations most relevant to the current study:

- `Thymic_fibroblasts`
- `Thymic_epithelial`

For each LOX-family gene, the script computed:

- raw counts;
- library size;
- CPM;
- `log2(CPM + 1)`;
- mean young expression;
- mean aged expression;
- aged-minus-young delta.

This is a descriptive external validation analysis, not a full differential-expression model.

## Results

| external cell type | gene | young n | aged n | delta log2(CPM+1), aged minus young | interpretation |
|---|---|---:|---:|---:|---|
| Thymic fibroblasts | `Lox` | 4 | 4 | -0.354 | Same broad aged-lower direction as current fibroblast/capsFB `Lox` findings |
| Thymic fibroblasts | `Loxl1` | 4 | 4 | -0.112 | Weak broad aged-lower direction; does not validate medFB `Loxl1` increase because fibroblast subtypes are pooled |
| Thymic fibroblasts | `Loxl2` | 4 | 4 | -1.206 | Strong broad aged-lower direction, externally consistent with current fibroblast `Loxl2` decrease |
| Thymic epithelial | `Loxl2` | 5 | 5 | -0.588 | Broad epithelial aged-lower direction, externally consistent with current mTEC1 `Loxl2` direction but not mTEC1-specific |

## Impact On Confidence

This analysis increases confidence that at least part of the LOX-family aging pattern is not unique to the original single-cell dataset:

- broad thymic fibroblast `Lox` is lower in aged samples in GSE223049;
- broad thymic fibroblast `Loxl2` is substantially lower in aged samples in GSE223049;
- broad thymic epithelial `Loxl2` is lower in aged samples in GSE223049.

The strongest external support is for fibroblast `Loxl2` and epithelial `Loxl2` directionality. The support for fibroblast `Lox` is directionally consistent but smaller. The external dataset does not support or refute the medFB `Loxl1` increase because it lacks fibroblast subtype resolution.

## What This Does Not Prove

This does not validate the precise subtype-specific claims:

- capsFB `Lox` decrease remains only partially externally validated through broad thymic fibroblast `Lox`;
- medFB `Loxl1` increase is not validated;
- medFB `Loxl2` decrease is only broadly supported through thymic fibroblast `Loxl2`;
- mTEC1 `Loxl2` decrease is only broadly supported through thymic epithelial `Loxl2`.

The dataset is sorted bulk RNA-seq rather than single-cell RNA-seq, so it cannot test detection-rate mechanisms or marker-defined mTEC1/capsFB/medFB subtypes. It also uses older mice, 22-24 months, compared with the current study's 18-month aged group.

## Bottom Line

The highest-impact next analysis was independent external validation. GSE223049 provides meaningful external support for aged-lower thymic fibroblast `Lox`, aged-lower thymic fibroblast `Loxl2`, and aged-lower thymic epithelial `Loxl2`. This strengthens confidence in the broad direction of selected LOX-family aging signals, especially `Loxl2`, while preserving the safer manuscript framing: subtype-specific claims still require subtype-resolved independent validation.
