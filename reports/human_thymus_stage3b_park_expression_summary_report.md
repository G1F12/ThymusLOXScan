# Park Human TEC LOX-Family Expression Summary

Search date: 2026-07-06

## Dataset and Local File Status

- Local H5AD path: `data\external\human_thymus\cellxgene_park_tec\park_tec.h5ad`
- AnnData shape: 18524 cells x 32088 features
- Local H5AD is under ignored external data path.

## Matrix Semantics Decision

- matrix_semantics_clear: partial
- detection matrix used: X
- mean-value matrix used: not_selected
- mean-value summaries generated: no
- reason: X and raw.X are present and target genes are available, but the H5AD does not provide enough local matrix annotation to label values as raw counts or normalized expression for mean-value summaries.

## Fields Used

- donor/sample fields: `donor_id`, `sample_id`
- age/development field: `development_stage`
- sex field: `sex`
- fine annotation field: `celltypes`
- broad annotation field: `cell_type`

## Dataset Breadth

- donors: 16
- samples: 33
- development stages: 12
- fine TEC labels: 9
- broad TEC labels: 3
- fine donor/sample/subtype groups: 156
- broad donor/sample/subtype groups: 58

## Gene Coverage

- LOX-family genes present: LOX, LOXL1, LOXL2, LOXL3, LOXL4

## Cell-Count Balance

Cell counts are uneven across donor/sample/development/subtype groups. Stage 3B outputs therefore summarize donor/sample groups and do not use cells as independent biological replicates.

## Detection-Summary Overview

Detection summaries were generated for fine and broad TEC annotation levels using target-gene values greater than zero in `X`.

## Mean-Value-Summary Status

Mean-value summaries were not generated because local matrix semantics remain partial.

## Stage 3C Recommendation

Stage 3C should proceed to the Yayon TEC dataset to test whether the same donor-aware parsing workflow can be applied in an independent human TEC resource.

## Limitations

- The Park dataset is primarily developmental/TEC context, with fetal and selected postnatal stages rather than an aged-adult-focused design.
- Matrix semantics remain partial for mean-value summaries.
- No statistical tests were run.
- No mouse-human comparison was performed.

This Stage 3B summary provides exploratory donor-aware transcript-level context for LOX-family genes in the Park human TEC dataset. It does not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, does not validate the mouse result, and does not provide mechanism, protein-level evidence, functional evidence, or therapeutic relevance.

Figures created: results\figures\human_thymus_stage3b_park\park_lox_detection_by_development_fine.png, results\figures\human_thymus_stage3b_park\park_lox_detection_by_celltype_fine.png, results\figures\human_thymus_stage3b_park\park_cell_counts_by_development.png