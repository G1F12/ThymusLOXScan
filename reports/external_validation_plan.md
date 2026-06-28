# External validation plan for LOX-family thymic aging findings

## Priority Order

### 1. GSE223049: sorted mouse cell-type bulk RNA-seq

Status: processed in this repository.

Script:

- `scripts/external_validation_gse223049_lox.py`

Outputs:

- `results/tables/external_gse223049_lox_validation.tsv`
- `results/tables/external_gse223049_lox_validation_summary.tsv`

Validation target:

- Broad thymic fibroblast `Lox`, `Loxl1`, `Loxl2`
- Broad thymic epithelial `Loxl2`

Interpretation:

- Supports lower aged thymic fibroblast `Lox` and especially `Loxl2`.
- Supports lower aged thymic epithelial `Loxl2`.
- Does not test capsFB, medFB, or mTEC1 specificity.
- Does not validate medFB `Loxl1` increase because fibroblast subtypes are pooled.

Next step:

- Add a simple external validation figure showing per-sample log2(CPM+1) for `Lox`, `Loxl1`, and `Loxl2` in thymic fibroblasts and `Loxl2` in thymic epithelial cells.

### 2. E-MTAB-8560: mouse TEC single-cell age-series

Status: metadata inspected through BioStudies/ArrayExpress API.

Known metadata:

- Species: mouse.
- Sample count: 11,682 single-cell samples.
- Ages: 1, 4, 16, 32, and 52 weeks, plus NA week.
- Cell types include cTEC, mTEClo, mTEChi, Dsg3+TEC, and mini-bulk categories.

Validation target:

- TEC and mTEC-like LOX-family age trends, especially `Loxl2`.

Planned processing:

1. Download SDRF and processed expression files from ArrayExpress/BioStudies.
2. Parse cell-level metadata for age, individual, and TEC subtype.
3. Create pseudobulk profiles by individual x age x TEC subtype.
4. Test `Lox`, `Loxl1`, `Loxl2`, `Loxl3`, and `Loxl4` trends across age, using age as ordered or continuous where metadata supports it.
5. Compare broad mTEC-like `Loxl2` trend with current mTEC1 `Loxl2` decrease, while noting annotation mismatch.

Limitations:

- TEC-only; no fibroblast validation.
- First-year mouse age-series is not the same as 2mo vs 18mo.
- Exact mTEC1 mapping is not guaranteed.

### 3. GSE231906: human thymus single-cell aging dataset

Status: cell-level metadata downloaded and inspected; full expression matrix not downloaded.

Local metadata:

- `tmp/external_metadata/GSE231906_cell-level_metadata.xlsx`

Useful metadata:

- Thymocyte/stromal-associated sheet: 213,620 cells across 37 thymus donors.
- Age range includes infancy through older adults up to 70 years.
- Major labels include `Mes`, `Epi`, `Endo`, and hematopoietic compartments.
- Fine labels include fibroblast-like `Fb_1`, `Fb_2`, cTEC, mTEC, post-AIRE mTEC, endothelial, VSMCs, and other TEC states.
- Integrated stromal sheet has 92,218 cells with identities including `Fb_1`, `Fb_2`, `Endo`, `VSMCs`, `cTEC_lo`, `cTEC_hi`, `mTEC_lo`, `mTEC_hi`, and `post_AIRE_mTEC`.

Validation target:

- Human fibroblast-like LOX-family age associations.
- Human epithelial/mTEC-like `LOXL2` age association.

Planned processing:

1. Download `GSE231906_RAW.tar` only when sufficient disk/time is available.
2. Extract expression matrices and match barcodes to `thymocyte_metadata` and `stromal_cell_metadata`.
3. Normalize gene symbols to human LOX family: `LOX`, `LOXL1`, `LOXL2`, `LOXL3`, `LOXL4`.
4. Build donor-level pseudobulk summaries for fibroblast-like (`Fb_1`, `Fb_2`, `Mes`) and epithelial (`cTEC`, `mTEC`, `post_AIRE_mTEC`) compartments.
5. Model age as continuous years and, secondarily, binned age groups.
6. Include source/dataset covariates where possible because integrated stromal metadata pools Deng, Bautista, and Park sources.

Limitations:

- Human thymus samples vary by donor source and age distribution.
- Integrated stromal metadata does not uniformly expose age in the stromal sheet itself; donor/source matching is required.
- LOX-family expression in human thymus may be development-, disease-, or source-sensitive.
- This can validate conservation or directionality only cautiously.

### 4. GSE147520 / E-MTAB-8581: human thymus atlas

Status: metadata-level candidate.

Validation target:

- Human thymic stromal baseline expression and developmental/adult comparison.

Use:

- Reference atlas for which compartments express LOX-family genes.
- Not a primary aging validation.

Limitations:

- Mostly developmental, fetal, pediatric, and young/adult samples.
- Does not directly test old-age thymic involution.

## External Claim Mapping

| Current finding | Best external validation source | Expected validation level |
|---|---|---|
| capsFB `Lox` decrease | GSE223049 thymic fibroblasts for broad direction; GSE231906 fibroblast-like subsets for subtype-like human check | Partial |
| medFB `Loxl1` increase | GSE231906 fibroblast-like subsets if subtype mapping is possible | Uncertain |
| medFB `Loxl2` decrease | GSE223049 broad thymic fibroblasts; GSE231906 fibroblast-like subsets | Moderate for broad fibroblast, uncertain for medFB |
| mTEC1 `Loxl2` decrease/low detection | E-MTAB-8560 mTEC-like age-series; GSE223049 broad thymic epithelial; GSE231906 human mTEC-like subsets | Moderate for epithelial direction, uncertain for exact mTEC1 |

## Immediate Result From Processed External Dataset

GSE223049 supports broad external consistency for several directions:

- Thymic fibroblast `Lox`: aged lower than young, delta log2(CPM+1) = -0.354.
- Thymic fibroblast `Loxl2`: aged lower than young, delta log2(CPM+1) = -1.206.
- Thymic epithelial `Loxl2`: aged lower than young, delta log2(CPM+1) = -0.588.

It does not support or refute medFB `Loxl1` increase because the dataset lacks fibroblast subtype resolution.

## Recommended Next Implementation

1. Add plotting for the processed GSE223049 table.
2. Add an E-MTAB-8560 downloader/parser for SDRF plus processed count files.
3. Add a guarded GSE231906 downloader with a disk-space warning before downloading the ~3.7 GB expression tarball.
4. Reuse the existing project output pattern: sample-level LOX tables, age-summary TSVs, and a short Markdown interpretation report per external dataset.
