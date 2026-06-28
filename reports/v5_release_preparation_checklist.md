# v5 release preparation checklist

Final v5 manuscript candidate: `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md`

Suggested tag: `v5.0-external-validation`

Suggested title: External-validation update v5.0

Suggested release description:

"External-validation update adding GSE223049 descriptive statistics, E-MTAB-8560 TEC age-series context, guarded GSE231906 metadata-only planning, and cross-dataset LOX-family direction summaries. The manuscript remains hypothesis-generating and does not claim causal, protein-level, functional, therapeutic, human-conservation, or exact subtype-resolved external-validation evidence."

## README Public Manuscript Check

- README points to v5.1 as the current external-validation manuscript candidate.
- README describes `v4.2-final-safe` as the previous stable cautious release.

## Include In v5 Release Candidate

- FOUND: `scripts/external_validation_gse223049_lox.py`
- FOUND: `scripts/external_validation_emtab8560_tec_lox.py`
- FOUND: `scripts/external_validation_gse231906_human_lox.py`
- FOUND: `scripts/build_cross_dataset_lox_matrix.py`
- FOUND: `scripts/figures/plot_cross_dataset_lox_validation.py`
- FOUND: `reports/gse223049_external_validation_v2.md`
- FOUND: `reports/emtab8560_tec_external_validation.md`
- FOUND: `reports/gse231906_human_external_validation.md`
- FOUND: `reports/gse231906_human_metadata_only_plan.md`
- FOUND: `reports/external_dataset_search_v2.md`
- FOUND: `reports/external_validation_priority_matrix.md`
- FOUND: `reports/cross_dataset_lox_analysis.md`
- FOUND: `reports/cross_dataset_figure_notes.md`
- FOUND: `reports/v5_file_and_reproducibility_check.md`
- FOUND: `reports/v5_hostile_overclaim_audit.md`
- FOUND: `reports/v5_1_overclaim_fixes_applied.md`
- FOUND: `reports/v5_1_safety_check.md`
- FOUND: `reports/v5_1_file_and_safety_final_check.md`
- FOUND: `reports/v5_release_preparation_checklist.md`
- FOUND: `results/tables/external_gse223049_lox_validation_stats.tsv`
- FOUND: `results/tables/external_gse223049_lox_validation_summary_v2.tsv`
- FOUND: `results/tables/external_emtab8560_tec_lox_by_age.tsv`
- FOUND: `results/tables/external_emtab8560_tec_lox_summary.tsv`
- FOUND: `results/tables/external_gse231906_human_lox_summary.tsv`
- FOUND: `results/tables/external_gse231906_human_lox_by_donor.tsv`
- FOUND: `results/tables/cross_dataset_lox_validation_matrix.tsv`
- FOUND: `results/tables/cross_dataset_lox_consistency_summary.tsv`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.png`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.pdf`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.png`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.pdf`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.png`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.pdf`
- FOUND: `manuscript/LOX_thymus_aging_public_preprint_v5_external_validation.md`
- FOUND: `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md`

## Exclude From Release Staging

- Large downloaded raw data under `data/external/`.
- Raw or processed `.h5ad` files.
- `tmp/`.
- Old drafts not intended for release.

## Files Over 50 MB

- `data/processed/thymus_annotated.h5ad` (763.5 MB)
- `data/processed/thymus_preprocessed.h5ad` (763.5 MB)
- `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad` (765.7 MB)

## Suggested Git Commands Only

```powershell
git status --short
git add manuscript/LOX_thymus_aging_public_preprint_v5_external_validation.md manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md reports/gse223049_external_validation_v2.md reports/emtab8560_tec_external_validation.md reports/gse231906_human_external_validation.md reports/gse231906_human_metadata_only_plan.md reports/external_dataset_search_v2.md reports/external_validation_priority_matrix.md reports/cross_dataset_lox_analysis.md reports/cross_dataset_figure_notes.md reports/v5_file_and_reproducibility_check.md reports/v5_hostile_overclaim_audit.md reports/v5_1_overclaim_fixes_applied.md reports/v5_1_safety_check.md reports/v5_1_file_and_safety_final_check.md reports/v5_release_preparation_checklist.md scripts/external_validation_gse223049_lox.py scripts/external_validation_emtab8560_tec_lox.py scripts/external_validation_gse231906_human_lox.py scripts/build_cross_dataset_lox_matrix.py scripts/figures/plot_cross_dataset_lox_validation.py results/tables/external_gse223049_lox_validation_stats.tsv results/tables/external_gse223049_lox_validation_summary_v2.tsv results/tables/external_emtab8560_tec_lox_by_age.tsv results/tables/external_emtab8560_tec_lox_summary.tsv results/tables/external_gse231906_human_lox_summary.tsv results/tables/external_gse231906_human_lox_by_donor.tsv results/tables/cross_dataset_lox_validation_matrix.tsv results/tables/cross_dataset_lox_consistency_summary.tsv results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.png results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.pdf results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.png results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.pdf results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.png results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.pdf
git commit -m "Prepare v5 external-validation update"
git tag v5.0-external-validation
# push only after final review:
git push origin HEAD
git push origin v5.0-external-validation
```

No commit, push, deletion, or v4.2 overwrite was performed by this checklist task.
