# E-MTAB-8560 Biological Unit Audit

Classification: partially verified: mouse IDs available but batch/age partly confounded.

## Biological Unit

The true biological unit used for downstream analysis is `age_week + Characteristics[individual]`, joined from official ArrayExpress/BioStudies SDRF metadata onto R-exported MouseThymusAgeing colData by normalized cell/source/assay identifiers.

SortDay and PlateID were not treated as mouse IDs.

Cell-to-mouse mapped fraction: 1.000.

## Mice Per Group

| group | age_week | n mice | n cells |
|---|---:|---:|---:|
| cTEC | 1 | 5 | 141 |
| cTEC | 4 | 5 | 175 |
| cTEC | 16 | 5 | 81 |
| cTEC | 32 | 5 | 144 |
| cTEC | 52 | 5 | 124 |
| combined_mTEClo_mTEChi | 1 | 5 | 221 |
| combined_mTEClo_mTEChi | 4 | 5 | 325 |
| combined_mTEClo_mTEChi | 16 | 5 | 192 |
| combined_mTEClo_mTEChi | 32 | 5 | 272 |
| combined_mTEClo_mTEChi | 52 | 5 | 287 |
| mTEChi | 1 | 5 | 129 |
| mTEChi | 4 | 5 | 193 |
| mTEChi | 16 | 5 | 99 |
| mTEChi | 32 | 5 | 143 |
| mTEChi | 52 | 5 | 151 |
| mTEClo | 1 | 5 | 92 |
| mTEClo | 4 | 5 | 132 |
| mTEClo | 16 | 5 | 93 |
| mTEClo | 32 | 5 | 129 |
| mTEClo | 52 | 4 | 136 |

## Batch Fields

- SortDay: {1: 5, 4: 5, 16: 5, 32: 5, 52: 5}
- PlateID: {1: 13, 4: 16, 16: 14, 32: 12, 52: 11}
- Comment[ENA_RUN]: {1: 433, 4: 599, 16: 318, 32: 496, 52: 481}

Batch-aware modeling is limited because SortDay, PlateID, and ENA run are not simple mouse-level covariates. The downstream analysis uses per-mouse aggregation and reports batch limitations explicitly.

Cells are not treated as independent biological replicates.
