# GSE147520 Human Thymic Epithelial H5AD Schema Report

Search date: 2026-07-06

The local GSE147520 epithelial H5AD was parsed for schema and gene-presence readiness only.

- AnnData shape: 14217 cells x 804 features
- Raw object exists: yes
- Layers found: none
- Donor field: not available in observed metadata
- Sample field: samples
- Age/development field: samples
- Sex field: not available in observed metadata
- Fine cell-type field: cell_types_epith
- Broad cell-type field: derived from cell_types_epith for summary grouping
- LOX-family genes detected: LOX, LOXL1, LOXL2, LOXL3, LOXL4

The target genes are absent from compact X but present in raw.var/raw.X, so downstream Stage 3D detection summaries should use raw.X while avoiding mean-value summaries unless value semantics are clear.

This report does not make a human expression conclusion.