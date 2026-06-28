# GSE231906 human external validation

## Guard status

- GEO raw archive: `GSE231906_RAW.tar`, expected size 3.69 GB.
- Available disk at data directory: 231.15 GB.
- Metadata workbook expected size: 0.017 GB.
- Raw archive present locally: False.
- Extracted matrix candidates found: 0.
- Expression parser status: no local expression matrix candidate found.

## Metadata-derived candidate groups

| candidate group | sheet | cells | donors/samples | min age years | max age years |
|---|---|---:|---:|---:|---:|
| cTEC_like | stromal_cell_metadata | 10683 | 45 | NA | NA |
| cTEC_like | thymocyte_metadata | 6000 | 28 | 1.75 | 70 |
| endothelial | stromal_cell_metadata | 17665 | 55 | NA | NA |
| endothelial | thymocyte_metadata | 1949 | 28 | 4 | 70 |
| epithelial | stromal_cell_metadata | 20025 | 55 | NA | NA |
| epithelial | thymocyte_metadata | 10637 | 30 | 1.75 | 70 |
| fibroblast_like_mesenchymal | stromal_cell_metadata | 53880 | 59 | NA | NA |
| fibroblast_like_mesenchymal | thymocyte_metadata | 36659 | 31 | 1.75 | 70 |
| mTEC_like | stromal_cell_metadata | 5268 | 48 | NA | NA |
| mTEC_like | thymocyte_metadata | 2053 | 27 | 1.75 | 70 |
| other | ETP_TP_metadata | 36368 | 49 | NA | NA |
| other | PBMC_metadata | 174142 | 22 | 0.833 | 105 |
| other | stromal_cell_metadata | 648 | 40 | NA | NA |
| other | thymocyte_metadata | 164375 | 37 | 0.0329 | 70 |

## Claim Classification

| group | gene | classification | reason |
|---|---|---|---|
| fibroblast_like_mesenchymal | LOX | not testable | Only metadata was parsed; LOX-family expression was not available for donor-level pseudobulk analysis. |
| fibroblast_like_mesenchymal | LOXL1 | not testable | Only metadata was parsed; LOX-family expression was not available for donor-level pseudobulk analysis. |
| fibroblast_like_mesenchymal | LOXL2 | not testable | Only metadata was parsed; LOX-family expression was not available for donor-level pseudobulk analysis. |
| mTEC_like | LOXL2 | not testable | Only metadata was parsed; LOX-family expression was not available for donor-level pseudobulk analysis. |

## Interpretation

No human LOX-family expression analysis was performed. This run is metadata-only and must not be used to claim human conservation of LOX-family aging patterns. GSE231906 is only a candidate future validation dataset. The metadata identifies promising fibroblast-like/mesenchymal, epithelial, mTEC-like, cTEC-like, and endothelial compartments, but expression is required for validation.

If expression is later parsed, the analysis should use donor-level pseudobulk summaries, model age as a continuous variable, include source/study covariates where possible, and avoid treating cells as independent biological replicates.