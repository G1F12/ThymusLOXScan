# External dataset search for LOX-family thymic aging validation

## Search Scope

Searched GEO, ArrayExpress/BioStudies, CELLxGENE-facing resources, and PubMed-linked supplementary data for public thymus datasets with aging or age-series structure, mouse or human samples, RNA-seq or single-cell RNA-seq, and stromal, epithelial, or fibroblast-like compartments where possible.

Search date: 2026-06-27.

## Candidate Datasets

| accession | source | species | age groups | sample count | data type | cell types available | LOX-family validation feasible? | limitations |
|---|---|---|---|---:|---|---|---|---|
| GSE223049 | GEO | Mouse | Young 2 months vs old 22-24 months | 170 processed bulk RNA-seq samples across 23 immune/stromal cell types; target thymus subsets include 4 young + 4 aged thymic fibroblast samples and 5 young + 5 aged thymic epithelial samples | Bulk RNA-seq of sorted cell types | Includes `Thymic_fibroblasts` and `Thymic_epithelial` among many sorted cell types | Yes, immediately feasible for broad thymic fibroblast and broad thymic epithelial LOX-family validation | No fibroblast subtype resolution; no mTEC1-specific analysis; sorted bulk cell-type definitions may differ from GSE240016 annotations |
| E-MTAB-8560 | ArrayExpress/BioStudies | Mouse | 1, 4, 16, 32, and 52 weeks; metadata also contains NA week | 11,682 single-cell samples | Smart-seq2 single-cell RNA-seq | TEC-enriched: cTEC, mTEClo, mTEChi, Dsg3+TEC, minibulk samples | Feasible for epithelial LOX-family age-series, especially broad cTEC/mTEC validation | No fibroblast/stromal fibroblast compartment; age range is first year of life, not very old mouse aging; processing requires ArrayExpress SDRF/count retrieval and cell-level parsing |
| E-MTAB-8737 | ArrayExpress/BioStudies | Mouse | Age-altered TEC differentiation context; metadata exposed as replicate factor | 6 samples | Single-cell RNA-seq / lineage tracing-related TEC study | TEC lineage-tracing samples from beta5t-expressing progenitors | Possibly useful as supportive TEC aging context | Small sample count; not a straightforward LOX validation dataset from metadata alone; no fibroblast compartment |
| GSE231906 | GEO, PubMed-linked supplement | Human | Broad human age range in metadata, including infancy, children, adults, and older adults up to 70 years; PBMC metadata includes 105 years | GEO lists 59 samples; downloaded cell metadata includes 213,620 thymocyte/stromal-associated cells across 37 thymus donors plus 92,218 integrated stromal cells | 10x single-cell RNA-seq | Metadata includes hematopoietic, mesenchymal, epithelial, endothelial, fibroblast-like (`Fb_1`, `Fb_2`), TEC, mTEC, cTEC, endothelial, VSMC-like labels | Feasible after full matrix download for human age-associated LOX-family analysis in fibroblast-like, epithelial, and endothelial compartments | Large expression matrix tarball (~3.7 GB); stromal metadata merges sources and not all integrated stromal records have donor ages in the same table; age, sampling, and source effects will need careful modeling |
| GSE147520 / E-MTAB-8581 | GEO / ArrayExpress | Human and mouse | Human fetal/postnatal/adult thymus; ages include fetal weeks, infants, adolescents, adults up to 35-40 years in ArrayExpress metadata | ArrayExpress lists 74 samples; GEO GSE147520 is linked to the human thymus cell atlas | Single-cell RNA-seq | Rich thymic epithelial, mesenchymal/fibroblast, endothelial, immune compartments in human thymus atlas | Partially feasible for developmental/adult baseline LOX expression across stromal compartments | Not a true aging dataset; mostly developmental and young/adult atlas; limited older-age validation |
| E-MTAB-11341 | ArrayExpress/BioStudies | Human | Fetal 11, 12, 15, 16, 18 weeks | 15 samples | 10x spatial transcriptomics | Human fetal liver, spleen, thymus spatial sections | Not suitable for aging validation; possible spatial developmental context only | Fetal only; spatial spots, not adult/aged thymus; no aging contrast |
| GSE240017 | GEO | Mouse | 2-month steady state and post-total-body-irradiation context from Kousa et al. | Bulk RNA-seq study linked to current source paper | Bulk RNA-seq of CD45-negative thymic cells | Broad CD45-negative thymic stromal cells | Not independent for aging validation; may help technical cross-checks only | Same project family as current dataset; not a clean young-vs-aged independent validation dataset |

## Already Processed Candidate

`GSE223049` was feasible and was processed immediately with `scripts/external_validation_gse223049_lox.py`.

Generated outputs:

- `data/external/GSE223049/GSE223049_RNA_seq_counts_23_cell_types.txt.gz`
- `data/external/GSE223049/GSE223049_RNA_seq_counts_23_cell_types.txt`
- `results/tables/external_gse223049_lox_validation.tsv`
- `results/tables/external_gse223049_lox_validation_summary.tsv`

Key first-pass external signals from GSE223049:

| cell type | gene | age contrast | delta log2(CPM+1), aged minus young | relevance |
|---|---|---|---:|---|
| Thymic fibroblasts | `Lox` | 22-24mo vs 2mo | -0.354 | Same broad direction as current pooled/capsFB `Lox`, but not subtype-specific |
| Thymic fibroblasts | `Loxl1` | 22-24mo vs 2mo | -0.112 | Weak broad decrease; does not validate medFB `Loxl1` increase because subtype resolution is absent |
| Thymic fibroblasts | `Loxl2` | 22-24mo vs 2mo | -1.206 | Strong broad fibroblast support for decreased `Loxl2` |
| Thymic epithelial | `Loxl2` | 22-24mo vs 2mo | -0.588 | Supports a broad epithelial decrease, but not mTEC1-specific |

## Best Validation Uses

1. Use GSE223049 for independent mouse bulk validation of broad thymic fibroblast and thymic epithelial LOX-family direction.
2. Use E-MTAB-8560 for mouse TEC age-series validation of epithelial LOX-family changes across cTEC/mTEC compartments.
3. Use GSE231906 for human single-cell validation, prioritizing fibroblast-like, mesenchymal, epithelial, cTEC, and mTEC labels after full matrix download.
4. Use GSE147520/E-MTAB-8581 as a human thymus atlas reference, not as primary aging validation.

## Bottom Line

Independent validation is feasible, but it will be uneven by claim:

- Broad fibroblast `Lox`/`Loxl2` direction can be checked immediately in GSE223049.
- medFB-specific `Loxl1` increase and `Loxl2` decrease require single-cell or subtype-resolved data; GSE231906 is the most promising human candidate, while E-MTAB-8560 is epithelial-only.
- mTEC1 `Loxl2` can be approximated in TEC datasets, but matching the exact mTEC1 annotation will require reannotation or mapping to dataset-specific TEC labels.
