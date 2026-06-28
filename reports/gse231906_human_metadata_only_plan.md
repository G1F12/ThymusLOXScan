# GSE231906 human metadata-only plan

## Why this stopped before full validation

The full GEO raw archive is 3.69 GB. The guarded pipeline does not download it automatically unless `--allow-large-download` is used and disk space is adequate. This avoids silently pulling a large TAR of CSV/MTX/TSV files.

## Manual download instructions

1. Download `GSE231906_RAW.tar` from `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE231906&format=file`.
2. Save it to `D:/ThymusLOXScan/data/external/GSE231906/GSE231906_RAW.tar`.
3. Inspect the TAR members before extraction, for example with Python `tarfile` or `tar -tf GSE231906_RAW.tar`.
4. Extract only the expression matrix, feature/gene, and barcode files needed for thymic stromal/epithelial cells.
5. Match barcodes to the metadata workbook by donor/sample/source columns.
6. Build donor-level pseudobulk summaries for LOX, LOXL1, LOXL2, LOXL3, and LOXL4 in fibroblast-like/mesenchymal, epithelial, mTEC-like, cTEC-like, and endothelial groups.
7. Analyze LOX-family expression versus age as donor-level continuous-age trends, with source/study covariates when possible.

## Parsed metadata outputs

- `results/tables/external_gse231906_human_lox_summary.tsv`
- `results/tables/external_gse231906_human_lox_by_donor.tsv`

## Current claim status

No human LOX-family expression analysis was performed. No human conservation claim is supported. GSE231906 is only a candidate future validation dataset.