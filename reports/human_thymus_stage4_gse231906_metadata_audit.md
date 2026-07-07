# GSE231906 Metadata Audit

## Purpose

Audit the GSE231906 cell-level metadata for donor/sample-aware human thymus aging context before expression parsing.

## Local Data

- Metadata workbook: `data/external/human_thymus/GSE231906/GSE231906_cell-level_metadata.xlsx`
- Workbook parsed: yes
- Cells represented in metadata: 516352
- Donor/sample labels available: yes
- Age field usable: yes
- Sex field usable: yes
- Target compartment rows identified: 140851

## Target Compartments

- ctec_like: 16683 cells, 45 donors/samples by donor field, 45 samples by sample field
- endothelial: 19614 cells, 55 donors/samples by donor field, 55 samples by sample field
- epithelial_or_epi: 30662 cells, 55 donors/samples by donor field, 55 samples by sample field
- immature_TEC: 4113 cells, 35 donors/samples by donor field, 35 samples by sample field
- mesenchymal_or_fibroblast_like: 90575 cells, 59 donors/samples by donor field, 59 samples by sample field
- mtec_like: 7321 cells, 48 donors/samples by donor field, 48 samples by sample field
- other: 375501 cells, 85 donors/samples by donor field, 85 samples by sample field
- post_AIRE_mTEC: 1609 cells, 33 donors/samples by donor field, 33 samples by sample field

## Interpretation

The metadata are technically usable for donor/sample-aware inventory work if the expression barcodes can be linked back to the workbook. Expression feasibility is not concluded by this metadata-only step.

This GSE231906 analysis provides aged-human thymus transcript-level context only. It does not establish human conservation of the mouse GSE240016 mTEC1 Loxl2 candidate signal, does not validate the mouse result, and does not provide mechanism, protein-level evidence, functional evidence, LOX activity evidence, or therapeutic relevance.
