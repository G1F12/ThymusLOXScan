# Cross-dataset LOX-family validation analysis

## Scope

This report combines internal GSE240016 pseudobulk results with available external outputs from GSE223049, E-MTAB-8560, and GSE231906. Broad sorted bulk evidence is kept broad and is not counted as subtype-specific validation.

## Supported Across Datasets

- `broad_FB_Lox_aged_lower`: 1 external exact or appropriately resolved row(s) match the expected direction; net score 1.
- `broad_FB_Loxl2_aged_lower`: 1 external exact or appropriately resolved row(s) match the expected direction; net score 1.
- `epithelial_Loxl2_aged_lower`: 2 external exact or appropriately resolved row(s) match the expected direction; net score 2.

## Only Internally Observed Or Externally Not Testable

- `capsFB_Lox_aged_lower`: internal direction is `aged-lower`, but external rows are not exact tests or are unavailable.
- `medFB_Loxl1_aged_higher`: internal direction is `aged-higher`, but external rows are not exact tests or are unavailable.
- `medFB_Loxl2_aged_lower`: internal direction is `aged-lower`, but external rows are not exact tests or are unavailable.

## Contradicted Findings

No claim is contradicted by the currently available external rows.

## Approximate Directional Support

- `mTEC1_Loxl2_aged_lower`: 3 mTEC-like/TEC-subset row(s) match the expected direction, but these are not exact mTEC1 replication.

## Why Exact Subtype Validation Remains Limited

- GSE223049 is sorted bulk RNA-seq with broad thymic fibroblast and thymic epithelial populations, so it cannot validate capsFB, medFB, or mTEC1 specificity.
- E-MTAB-8560 is TEC-only and useful for mTEC-like `Loxl2`, but it has no fibroblast compartment and does not share the exact mTEC1 annotation.
- GSE231906 currently contributes metadata-only candidate groups in this repository; expression was not parsed because the raw matrix archive is large and requires guarded donor/source matching.
- Therefore, broad fibroblast and epithelial directions can be compared across datasets, while medFB and mTEC1 specificity remain internally observed or approximately tested rather than independently subtype-validated.

## Matrix Outputs

- `results/tables/cross_dataset_lox_validation_matrix.tsv`
- `results/tables/cross_dataset_lox_consistency_summary.tsv`