# Supplementary tables

These tables are generated from existing repository data and analysis outputs. They do not rerun statistical models or alter scientific results.

Regenerate all supplementary tables from the repository root with:

```bash
python supplementary_tables/make_supplementary_tables.py
```

## Supplementary_Table_1_cell_counts.tsv

Generated from `data/processed/thymus_annotated.h5ad` using the annotated cell metadata columns `sample`, `age_group`, `cell_type`, and `cell_type_subset`. Rows are sample-by-age-by-cell-type-by-subtype counts. The `subtype` column corresponds to `cell_type_subset` in the AnnData object.

## Supplementary_Table_2_pseudobulk_LOX_results.tsv

Generated from `results/tables/lox_pseudobulk_complete_results.csv`, which contains the LOX-family DESeq2 pseudobulk results. The table preserves the reported DESeq2 statistics and replicate counts. Rows analyzed at the broad cell-type level have an empty `subtype`; rows analyzed at the subtype level have `cell_type` inferred from the majority parent annotation in `data/processed/thymus_annotated.h5ad`.

## Supplementary_Table_3_single_cell_tests.tsv

Generated from `results/sc_mannwhitney_FB_combined.csv` and `results/sc_mannwhitney_mTEC1.csv`. These are descriptive single-cell-level Mann-Whitney U tests for pooled fibroblasts and mTEC1, respectively. The source files do not include multiple-testing adjusted p-values, so `adjusted_pvalue` is left blank.

## Supplementary_Table_4_correlations.tsv

Generated from `results/sc_spearman_correlations.csv`. These are descriptive Spearman correlations reported for pooled fibroblast analyses, with the original `stage` column encoded in `compartment` as `pooled fibroblasts; stage=<stage>`. The source file does not include multiple-testing adjusted p-values, so `adjusted_pvalue` is left blank.
