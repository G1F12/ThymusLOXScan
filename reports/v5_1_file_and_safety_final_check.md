# v5.1 file and safety final check

Final manuscript candidate: `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md`

## Required Safety Checks

- PASS: v5.1 manuscript exists. D:/ThymusLOXScan/manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md
- PASS: README points to v5.1 as current external-validation manuscript. README current manuscript section names v5.1.
- PASS: v4.2 described only as previous stable release. README keeps v4.2 as previous stable cautious release.
- PASS: no TODO. No TODO/TBD tokens in v5.1 manuscript.
- PASS: no placeholder DOI. No placeholder DOI pattern.
- PASS: no fake bioRxiv ID. No fake bioRxiv ID pattern.
- PASS: no missing referenced local files. All referenced local files exist.
- PASS: no causal claim. Causal language appears only as negation or unresolved limitation.
- PASS: no protein-level claim. Protein-level mentions are negations or future-work requirements.
- PASS: no enzymatic activity claim. Enzymatic activity mentions are negations or future-work requirements.
- PASS: no ECM-crosslinking current-data claim. Crosslinking mentions are background, negations, or future-work requirements.
- PASS: no thymic functional current-data claim. Functional mentions are negations or future-work requirements.
- PASS: no therapeutic claim. No therapeutic language in v5.1 manuscript.
- PASS: no human-conservation claim. Human-conservation mentions explicitly deny support.
- PASS: no exact subtype-resolved external-validation overclaim. No exact subtype-resolved external-validation claim.
- PASS: no E-MTAB-8560 exact mTEC1 replication claim. E-MTAB-8560 is framed as approximate TEC/mTEC-like context.
- PASS: no GSE223049 medFB/capsFB replication claim. GSE223049 is framed as broad sorted bulk only.
- PASS: no GSE231906 expression-level human validation claim. GSE231906 is metadata-only.
- PASS: no public-facing workflow traces. No configured trace terms found in release files.
- PASS: all release scripts compile. Compiled: scripts/external_validation_gse223049_lox.py, scripts/external_validation_emtab8560_tec_lox.py, scripts/external_validation_gse231906_human_lox.py, scripts/build_cross_dataset_lox_matrix.py, scripts/figures/plot_cross_dataset_lox_validation.py

## Referenced Results Files

- FOUND: `results/figures/annotation_sanity/`
- FOUND: `results/figures/annotation_sanity/marker_expression_long.tsv`
- FOUND: `results/figures/annotation_sanity/marker_summary_by_sample.tsv`
- FOUND: `results/figures/annotation_sanity/strict_marker_positive_lox_age_summary.tsv`
- FOUND: `results/figures/annotation_sanity/strict_marker_positive_lox_by_sample.tsv`
- FOUND: `results/figures/annotation_sanity/strict_marker_positive_lox_summary.png`
- FOUND: `results/figures/external_validation/cross_dataset/`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.pdf`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_claim_consistency.png`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.pdf`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_direction_heatmap.png`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.pdf`
- FOUND: `results/figures/external_validation/cross_dataset/cross_dataset_lox_external_effect_points.png`
- FOUND: `results/figures/external_validation/emtab8560/`
- FOUND: `results/figures/final/`
- FOUND: `results/figures/final/Fig1_pseudobulk_volcano_final.pdf`
- FOUND: `results/figures/final/Fig1_pseudobulk_volcano_final.png`
- FOUND: `results/figures/final/Fig2_summary_4panel_final.pdf`
- FOUND: `results/figures/final/Fig2_summary_4panel_final.png`
- FOUND: `results/figures/per_sample/lox_per_sample_pseudobulk_expression.pdf`
- FOUND: `results/figures/per_sample/mtec1_loxl2_per_sample.pdf`
- FOUND: `results/magic_imputation_LOX.csv`
- FOUND: `results/sc_mannwhitney_FB_combined.csv`
- FOUND: `results/sc_mannwhitney_mTEC1.csv`
- FOUND: `results/sc_spearman_correlations.csv`
- FOUND: `results/tables/cross_dataset_lox_consistency_summary.tsv`
- FOUND: `results/tables/cross_dataset_lox_validation_matrix.tsv`
- FOUND: `results/tables/external_emtab8560_tec_lox_by_age.tsv`
- FOUND: `results/tables/external_emtab8560_tec_lox_summary.tsv`
- FOUND: `results/tables/external_gse223049_lox_validation.tsv`
- FOUND: `results/tables/external_gse223049_lox_validation_stats.tsv`
- FOUND: `results/tables/external_gse223049_lox_validation_summary.tsv`
- FOUND: `results/tables/external_gse223049_lox_validation_summary_v2.tsv`
- FOUND: `results/tables/external_gse231906_human_lox_by_donor.tsv`
- FOUND: `results/tables/external_gse231906_human_lox_summary.tsv`
- FOUND: `results/tables/fb_composition_adjusted_lox.tsv`
- FOUND: `results/tables/lox_detection_rates_by_sample.tsv`
- FOUND: `results/tables/lox_pseudobulk_complete_results.csv`
- FOUND: `results/tables/lox_pseudobulk_complete_results.tsv`
- FOUND: `results/tables/mtec1_loxl2_per_sample_expression.tsv`

## Script Compile Checks

- PASS: `python -m py_compile scripts/external_validation_gse223049_lox.py`
- PASS: `python -m py_compile scripts/external_validation_emtab8560_tec_lox.py`
- PASS: `python -m py_compile scripts/external_validation_gse231906_human_lox.py`
- PASS: `python -m py_compile scripts/build_cross_dataset_lox_matrix.py`
- PASS: `python -m py_compile scripts/figures/plot_cross_dataset_lox_validation.py`

## Final Determination

v5.1 passes the final file, safety, public-cleanup, and script compile checks. It remains a computational, hypothesis-generating transcript-level reanalysis.
