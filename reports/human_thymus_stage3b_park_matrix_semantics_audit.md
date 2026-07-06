# Park Human TEC Matrix Semantics Audit

Search date: 2026-07-06

The Park TEC H5AD was inspected for target-gene matrix behavior without writing large derived matrices.

- X matrix: {"class": "anndata._core.sparse_dataset._CSRDataset", "dtype": "float32", "shape": "(18524, 32088)", "storage": "sparse"}
- raw.X matrix: {"class": "anndata._core.sparse_dataset._CSRDataset", "dtype": "float32", "shape": "(18524, 32088)", "storage": "sparse"}
- Layers: none
- Schema reference: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md
- Schema version: 7.1.0
- Target genes compared: LOX, LOXL1, LOXL2, LOXL3, LOXL4

## Conclusion

- matrix_semantics_clear: partial
- recommended_detection_matrix: X
- recommended_mean_value_matrix: not_selected
- mean_value_allowed: no
- reason: X and raw.X are present and target genes are available, but the H5AD does not provide enough local matrix annotation to label values as raw counts or normalized expression for mean-value summaries.

Detection summaries can proceed using `X > 0`. Mean-value summaries are not generated in Stage 3B because value semantics remain insufficiently documented locally.