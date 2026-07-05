# E-MTAB-8560 mTEC LOX Reanalysis

Final classification: mixed/inconclusive.

This is an independent mTEC-focused transcript-level test. It is not exact GSE240016 mTEC1 replication.

## Biological Unit

One dot/table row at the inferential level is one biological mouse ID defined as `age_week + Characteristics[individual]` after joining official SDRF metadata to MouseThymusAgeing colData. Cells are not treated as independent biological replicates.

## Loxl2 Summary

| group | 4w vs 52w delta | 16w vs 52w delta | 4+16w vs 52w delta |
|---|---:|---:|---:|
| mTEClo | -2.53 | -1.269 | -1.9 |
| mTEChi | -0.9707 | -0.1622 | -0.5665 |
| combined_mTEClo_mTEChi | -1.304 | -0.4735 | -0.8887 |
| cTEC | -2.847 | -0.7336 | -1.79 |

## Batch Status

Batch-adjusted models were not estimated because SortDay, PlateID, and run fields are complex cell/acquisition-level variables rather than simple mouse-level covariates. This is reported as a limitation.

## Output Files

- `results/tables/external_emtab8560_mtec_lox_pseudobulk.tsv`
- `results/tables/external_emtab8560_mtec_lox_by_age.tsv`
- `results/tables/external_emtab8560_mtec_lox_models.tsv`
- `results/tables/external_emtab8560_mtec_lox_permutation.tsv`
- `results/tables/external_emtab8560_mtec_lox_global_fdr.tsv`
