# v5.2 patch release file check

Date: 2026-06-28

## Required files
- Present: `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`
- Present: `reports/v5_2_dropout_protein_update_summary.md`
- Present: `reports/v5_2_safety_check.md`
- Present: `reports/v5_2_file_check.md`
- Present: `one_page_summary.md`
- Present: `wet_lab_validation_plan.md`
- Present: `reports/v5_2_outreach_docs_safety_check.md`
- Present: `scripts/mtec1_loxl2_dropout_depth_audit.py`
- Present: `reports/mtec1_loxl2_dropout_depth_audit.md`
- Present: `reports/loxl2_protein_spatial_validation_feasibility.md`
- Present: `results/tables/mtec1_loxl2_dropout_depth_audit.tsv`
- Present: `results/tables/mtec1_loxl2_pairwise_sample_direction.tsv`
- Present: `results/tables/mtec1_loxl2_downsampling_summary.tsv`
- Present: `results/tables/mtec1_loxl2_depth_adjusted_logistic.tsv`
- Present: `results/tables/loxl2_protein_spatial_evidence_audit.tsv`
- Present: `results/tables/loxl2_ihc_antibody_feasibility.tsv`

## Optional staged figure assets requested for the v5.2 patch
- Present: `results/figures/dropout/mtec1_loxl2_detection_by_sample.png`
- Present: `results/figures/dropout/mtec1_loxl2_detection_by_sample.pdf`
- Present: `results/figures/dropout/mtec1_loxl2_depth_qc_by_sample.png`
- Present: `results/figures/dropout/mtec1_loxl2_depth_qc_by_sample.pdf`
- Present: `results/figures/dropout/mtec1_loxl2_detection_vs_depth.png`
- Present: `results/figures/dropout/mtec1_loxl2_detection_vs_depth.pdf`
- Present: `results/figures/dropout/mtec1_loxl2_downsampling_sensitivity.png`
- Present: `results/figures/dropout/mtec1_loxl2_downsampling_sensitivity.pdf`

## Exclusion checks
- No `.h5ad` files are required for this release.
- No `data/` or `data/external/` files are required for this release.
- No large raw or processed data files are required for this release.
- Existing v5.0 and v5.1 release files are not overwritten by this patch file set.

## Conclusion
All required v5.2 patch-release files are present. The intended release payload is manuscript, report, script, small table, small figure, and outreach-document focused.
