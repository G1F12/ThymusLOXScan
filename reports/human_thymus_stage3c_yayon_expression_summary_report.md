# Yayon Human TEC LOX-Family Expression Summary

Search date: 2026-07-06

## Dataset and Local File Status

- Local H5AD path: `data\external\human_thymus\cellxgene_yayon_tec\yayon_tec.h5ad`
- AnnData shape: 25726 cells x 35477 features
- Local H5AD is under ignored external data path.

## Matrix Semantics Decision

- matrix_semantics_clear: partial
- detection matrix used: X
- mean-value matrix used: not_selected
- mean-value summaries generated: no
- reason: The object has no raw layer and no local matrix annotation sufficient to label X values for mean-value summaries.

## Fields Used

- donor/sample fields: `donor_id`, `sample`
- age/development field: `development_stage`
- sex field: `sex`
- fine annotation field: `cell_type_level_4_explore`
- broad annotation field: `cell_type_level_2`

## Dataset Breadth

- donors: 18
- samples: 30
- development stages: 13
- fine TEC labels: 13
- broad TEC labels: 5
- fine donor/sample/subtype groups: 261
- broad donor/sample/subtype groups: 114

## Gene Coverage

- LOX-family genes present: LOX, LOXL1, LOXL2, LOXL3, LOXL4

## Cell-Count Balance

Cell counts are uneven across donor/sample/development/subtype groups. Stage 3C outputs therefore summarize donor/sample groups and do not use cells as independent biological replicates.

## Detection-Summary Overview

Detection summaries were generated for fine and broad TEC annotation levels using target-gene values greater than zero in `X`.

## Mean-Value-Summary Status

Mean-value summaries were not generated because local matrix semantics remain partial.

## Park Workflow Context

The same cautious workflow used for Park was applied here: schema inspection, matrix-semantics audit, donor/sample-aware detection summaries, no statistical tests, and no cross-dataset biological conclusion.

## Stage 3D Recommendation

Stage 3D should proceed to the GSE147520 epithelial H5AD if local access and schema checks remain feasible.

## Limitations

- The Yayon TEC dataset is useful as human TEC context, with fetal and early postnatal labels rather than an aged-adult-focused design.
- Matrix semantics remain partial for mean-value summaries.
- No statistical tests were run.
- No mouse-human comparison was performed.

This Stage 3C summary provides exploratory donor-aware transcript-level context for LOX-family genes in the Yayon human TEC dataset. It does not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, does not validate the mouse result, and does not provide mechanism, protein-level evidence, functional evidence, or therapeutic relevance.

Figures created: results\figures\human_thymus_stage3c_yayon\yayon_lox_detection_by_development_fine.png, results\figures\human_thymus_stage3c_yayon\yayon_lox_detection_by_celltype_fine.png, results\figures\human_thymus_stage3c_yayon\yayon_cell_counts_by_development.png