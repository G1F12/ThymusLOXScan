# GSE147520 Human Thymic Epithelial Matrix Semantics Audit

Search date: 2026-07-06

The GSE147520 epithelial H5AD was inspected for target-gene matrix behavior without writing large derived matrices.

- X matrix: {'shape': '(14217, 804)', 'dtype': 'float32', 'storage': 'dense', 'class': 'h5py._hl.dataset.Dataset'}
- raw.X matrix: {'shape': '(14217, 26587)', 'dtype': 'float32', 'storage': 'sparse', 'class': 'anndata._core.sparse_dataset._CSRDataset'}
- Layers: none

## Conclusion

- matrix_semantics_clear: partial
- recommended_detection_matrix: raw.X
- recommended_mean_value_matrix: not_selected
- mean_value_allowed: no
- reason: Target genes are absent from compact X but present in raw.X; detection summaries use raw.X > 0, while mean-value summaries are not generated because value semantics remain partial.

Target-gene detection can proceed from `raw.X > 0` when target genes are present only in raw.var. Mean-value summaries are not generated in Stage 3D.