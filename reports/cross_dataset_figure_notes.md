# Cross-dataset LOX figure notes

## Scope

These figures summarize direction consistency only. They do not label any comparison as validated unless subtype-resolved independent replication is present; none of the current external rows meet that stronger standard for capsFB, medFB, or exact mTEC1 specificity.

## Generated figures

- `results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.png`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.pdf`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.png`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.pdf`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.png`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.pdf`

## Label definitions

- internally observed: the direction is observed in GSE240016 but exact external testing is unavailable.
- directionally consistent: an external row matches the expected direction at an appropriate broad resolution.
- approximate external context: an external row has a matching direction but not exact subtype resolution.
- opposite direction: an external row has the opposite direction at an appropriate resolution.
- not testable: the dataset lacks the needed compartment, expression matrix, or resolution.
- metadata only: candidate groups are present, but no expression analysis was performed.

## Claim status counts

- internally observed: 3
- directionally consistent: 3
- approximate external context: 1
- not testable: 0
- opposite direction: 0
- metadata only: 0

## Cautions

GSE223049 broad fibroblast and epithelial bulk rows are not subtype-specific validation. E-MTAB-8560 provides TEC and mTEC-like directional support, not exact mTEC1 replication. GSE231906 remains metadata-only in this repository, so its candidate human groups are not expression validation.