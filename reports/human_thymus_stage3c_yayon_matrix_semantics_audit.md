# Yayon Human TEC Matrix Semantics Audit

Search date: 2026-07-06

The Yayon TEC H5AD was inspected for target-gene matrix behavior without writing large derived matrices.

- X matrix: {"class": "anndata._core.sparse_dataset._CSRDataset", "dtype": "float32", "shape": "(25726, 35477)", "storage": "sparse"}
- raw.X matrix: {"class": "not_available", "dtype": "not_available", "shape": "not_available", "storage": "not_available"}
- Layers: none
- Schema reference: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md
- Schema version: 7.1.0

## Conclusion

- matrix_semantics_clear: partial
- recommended_detection_matrix: X
- recommended_mean_value_matrix: not_selected
- mean_value_allowed: no
- reason: The object has no raw layer and no local matrix annotation sufficient to label X values for mean-value summaries.

Detection summaries can proceed using `X > 0`. Mean-value summaries are not generated in Stage 3C because local value semantics remain insufficient for that output.