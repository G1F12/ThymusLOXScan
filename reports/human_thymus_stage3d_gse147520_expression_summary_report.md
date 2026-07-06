# GSE147520 Human Thymic Epithelial LOX-Family Expression Summary

Search date: 2026-07-06

## Dataset and Local File Status

- Local H5AD path: `data\external\human_thymus\GSE147520\GSE147520_epithelial_cells.h5ad`
- AnnData shape: 14217 cells x 804 features
- Local H5AD/H5AD.GZ files are under ignored external data path.

## Matrix Semantics Decision

- matrix_semantics_clear: partial
- detection matrix used: raw.X
- mean-value matrix used: not_selected
- mean-value summaries generated: no
- reason: Target genes are absent from compact X but present in raw.X; detection summaries use raw.X > 0, while mean-value summaries are not generated because value semantics remain partial.

## Fields Used

- donor field: not available in observed metadata; output uses `not_available` placeholder
- sample field: `samples`
- age/development field: `age_development`, derived from `samples`
- sex field: not available in observed metadata; output uses `not_available` placeholder
- fine annotation field: `cell_types_epith`
- broad annotation field: `cell_type_broad_derived`

## Dataset Breadth

- donors: 1 placeholder group
- samples/age groups: 5
- development-or-age groups: 5
- fine epithelial labels: 9
- broad derived epithelial labels: 6

## Gene Coverage

- LOX-family genes present: LOX, LOXL1, LOXL2, LOXL3, LOXL4 in raw.var/raw.X

## Cell-Count Balance

Cell counts are uneven across sample and epithelial-label groups. Outputs are grouped at sample/label level and do not use cells as independent biological replicates.

## Detection-Summary Overview

Target-gene detection summaries were generated from `raw.X > 0` because the target genes are present in raw.var/raw.X but absent from compact X.

## Mean-Value-Summary Status

Mean-value summaries were not generated.

## Park/Yayon Workflow Context

The same cautious workflow was applied at the file and metadata level, but GSE147520 differs from Park and Yayon because compact X omits the target genes while raw.X contains them.

## Stage 3E Recommendation

Stage 3E cross-human context matrix can proceed, but GSE147520 should be represented with clear layer-specific notes: compact X omits the target genes, while detection context can be summarized from raw.X.

## Limitations

- Donor and sex fields are not available in the observed H5AD metadata.
- Age/development information is embedded in the `samples` labels.
- Compact X lacks the target LOX-family genes; target-gene detection uses raw.X only.
- No statistical tests were run.
- No mouse-human comparison was performed.

This Stage 3D summary provides exploratory donor-aware transcript-level context for LOX-family genes in the GSE147520 human thymic epithelial dataset. It does not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, does not validate the mouse result, and does not provide mechanism, protein-level evidence, functional evidence, or therapeutic relevance.

Figures created: results\figures\human_thymus_stage3d_gse147520\gse147520_cell_counts_by_development.png, results\figures\human_thymus_stage3d_gse147520\gse147520_lox_detection_by_development_fine.png, results\figures\human_thymus_stage3d_gse147520\gse147520_lox_detection_by_celltype_fine.png