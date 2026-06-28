# v5 file and reproducibility check

Input manuscript: `manuscript/LOX_thymus_aging_public_preprint_v5_external_validation.md`

## Referenced Local Files Found

- `data/external/GSE223049/`
- `data/processed/thymus_annotated.h5ad`
- `data/processed/thymus_preprocessed.h5ad`
- `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad`
- `notebooks/03_preprocessing_thymus.ipynb`
- `notebooks/04_cell_type_annotation.ipynb`
- `notebooks/06_GSEA_TGFb_analysis.ipynb`
- `reports/annotation_uncertainty_check.md`
- `reports/cross_dataset_figure_notes.md`
- `reports/cross_dataset_lox_analysis.md`
- `reports/emtab8560_tec_external_validation.md`
- `reports/external_dataset_search.md`
- `reports/external_dataset_search_v2.md`
- `reports/external_validation_plan.md`
- `reports/fb_composition_adjustment.md`
- `reports/gse223049_external_validation_v2.md`
- `reports/gse231906_human_external_validation.md`
- `reports/gse231906_human_metadata_only_plan.md`
- `reports/highest_impact_next_analysis.md`
- `reports/hypothesis_falsification.md`
- `reports/lox_detection_rate_analysis.md`
- `reports/mtec1_loxl2_robustness.md`
- `results/figures/annotation_sanity/`
- `results/figures/annotation_sanity/marker_expression_long.tsv`
- `results/figures/annotation_sanity/marker_summary_by_sample.tsv`
- `results/figures/annotation_sanity/strict_marker_positive_lox_age_summary.tsv`
- `results/figures/annotation_sanity/strict_marker_positive_lox_by_sample.tsv`
- `results/figures/annotation_sanity/strict_marker_positive_lox_summary.png`
- `results/figures/external_validation/cross_dataset/`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.pdf`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.pdf`
- `results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.pdf`
- `results/figures/external_validation/emtab8560/`
- `results/figures/final/`
- `results/figures/final/Fig1_pseudobulk_volcano_final.pdf`
- `results/figures/final/Fig1_pseudobulk_volcano_final.png`
- `results/figures/final/Fig2_summary_4panel_final.pdf`
- `results/figures/final/Fig2_summary_4panel_final.png`
- `results/figures/per_sample/lox_per_sample_pseudobulk_expression.pdf`
- `results/figures/per_sample/mtec1_loxl2_per_sample.pdf`
- `results/magic_imputation_LOX.csv`
- `results/sc_mannwhitney_FB_combined.csv`
- `results/sc_mannwhitney_mTEC1.csv`
- `results/sc_spearman_correlations.csv`
- `results/tables/cross_dataset_lox_consistency_summary.tsv`
- `results/tables/cross_dataset_lox_validation_matrix.tsv`
- `results/tables/external_emtab8560_tec_lox_by_age.tsv`
- `results/tables/external_emtab8560_tec_lox_summary.tsv`
- `results/tables/external_gse223049_lox_validation.tsv`
- `results/tables/external_gse223049_lox_validation_stats.tsv`
- `results/tables/external_gse223049_lox_validation_summary.tsv`
- `results/tables/external_gse223049_lox_validation_summary_v2.tsv`
- `results/tables/external_gse231906_human_lox_by_donor.tsv`
- `results/tables/external_gse231906_human_lox_summary.tsv`
- `results/tables/fb_composition_adjusted_lox.tsv`
- `results/tables/lox_detection_rates_by_sample.tsv`
- `results/tables/lox_pseudobulk_complete_results.csv`
- `results/tables/lox_pseudobulk_complete_results.tsv`
- `results/tables/mtec1_loxl2_per_sample_expression.tsv`
- `scripts/annotation_uncertainty_check.py`
- `scripts/build_cross_dataset_lox_matrix.py`
- `scripts/external_validation_emtab8560_tec_lox.py`
- `scripts/external_validation_gse223049_lox.py`
- `scripts/external_validation_gse231906_human_lox.py`
- `scripts/fb_composition_adjusted_lox.py`
- `scripts/figures/plot_cross_dataset_lox_validation.py`
- `scripts/figures/plot_final_summary.py`
- `scripts/figures/plot_final_volcano.py`
- `scripts/lox_detection_rates_by_sample.py`
- `scripts/magic_imputation_lox.py`
- `scripts/make_pseudobulk_results_table.py`
- `scripts/plot_per_sample_lox_expression.py`
- `scripts/pseudobulk_deseq2_lox.py`
- `supplementary_tables/make_supplementary_tables.py`

## Missing Referenced Local Files

None.

## Script Compile Checks Passed

- `python -m py_compile scripts/external_validation_gse223049_lox.py`
- `python -m py_compile scripts/external_validation_emtab8560_tec_lox.py`
- `python -m py_compile scripts/external_validation_gse231906_human_lox.py`
- `python -m py_compile scripts/build_cross_dataset_lox_matrix.py`
- `python -m py_compile scripts/figures/plot_cross_dataset_lox_validation.py`

## Script Compile Checks Failed

None.

## Outputs Not Reproducible From Current Local Data Without Large Downloads

- `GSE231906` human LOX-family expression analysis is not reproducible from current local data because no expression matrix/raw archive is present locally; the guarded script intentionally produced metadata-only outputs.

## Files Over 50 MB

- `data/processed/thymus_annotated.h5ad` (763.5 MB)
- `data/processed/thymus_preprocessed.h5ad` (763.5 MB)
- `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad` (765.7 MB)

## Release-Readiness Assessment

v5 file/reproducibility checks are structurally passable, with the explicit caveat that GSE231906 remains metadata-only and not human expression validation.