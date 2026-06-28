# Final public release checklist

Checked final manuscript:

- `manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md`

Supporting files checked:

- `reports/v4_2_final_safety_check.md`
- `README.md`
- `environment.yml`
- `requirements.txt`
- `results/tables/`
- `results/figures/`
- `supplementary_tables/`

## Verification summary

Status: ready for a cautious public code/manuscript release.

The final manuscript is the safest public version currently prepared. It preserves the subtype-focused LOX-family transcript findings while avoiding causal, protein-level, enzymatic, ECM-crosslinking, thymic-functional, therapeutic, mediator, promotional-ranking, and subtype-specific external-validation overclaims.

## Manuscript text scan

| check | result |
|---|---|
| TODO/TBD markers | Pass: no matches in final manuscript. |
| Placeholder text | Pass: no placeholder/fake-template strings detected. |
| Fake DOI / placeholder DOI | Pass: no fake DOI or placeholder DOI strings detected. |
| Unresolved public claims | Pass: current-data claims are framed as transcript-level, computational, descriptive, or hypothesis-generating where appropriate. |
| Risky final wording from hostile review | Pass: no `most prominent`, `most recurrent`, `functional validation`, `possible mediator`, `strongest`, `validated`, `proves`, or `therapeutic` matches. |

Remaining cautious language such as "not subtype-resolved validation", "future experimental validation", and "does not establish" appears only as caveat/limitation wording or in existing file names.

## Referenced file check

All local repository paths referenced by the final manuscript were checked.

- Local paths checked: 52
- Missing local paths: 0

Result: pass.

## Figures and supplementary tables

All manuscript-referenced figures exist:

- `results/figures/final/Fig1_pseudobulk_volcano_final.pdf`
- `results/figures/final/Fig1_pseudobulk_volcano_final.png`
- `results/figures/final/Fig2_summary_4panel_final.pdf`
- `results/figures/final/Fig2_summary_4panel_final.png`
- `results/figures/per_sample/lox_per_sample_pseudobulk_expression.pdf`
- `results/figures/per_sample/mtec1_loxl2_per_sample.pdf`
- `results/figures/annotation_sanity/strict_marker_positive_lox_summary.png`

All manuscript-referenced annotation-sanity figure tables exist:

- `results/figures/annotation_sanity/marker_summary_by_sample.tsv`
- `results/figures/annotation_sanity/strict_marker_positive_lox_age_summary.tsv`

All supplementary tables referenced by table number exist:

- `supplementary_tables/Supplementary_Table_1_cell_counts.tsv`
- `supplementary_tables/Supplementary_Table_2_pseudobulk_LOX_results.tsv`
- `supplementary_tables/Supplementary_Table_3_single_cell_tests.tsv`
- `supplementary_tables/Supplementary_Table_4_correlations.tsv`

Result: pass.

## README check

`README.md` now points to the final public manuscript:

- `manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md`

Result: pass.

## Files to include in the release

Recommended core release contents:

- `README.md`
- `environment.yml`
- `requirements.txt`
- `manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md`
- `reports/v4_2_final_safety_check.md`
- `reports/final_public_release_checklist.md`
- `results/tables/`
- `results/figures/final/`
- `results/figures/per_sample/`
- `results/figures/annotation_sanity/`
- `supplementary_tables/`
- `scripts/`
- `notebooks/`, if publishing the exploratory record is desired

Recommended analysis/audit reports to include for transparency:

- `reports/fb_composition_adjustment.md`
- `reports/lox_detection_rate_analysis.md`
- `reports/annotation_uncertainty_check.md`
- `reports/external_dataset_search.md`
- `reports/external_validation_plan.md`
- `reports/highest_impact_next_analysis.md`
- `reports/hypothesis_falsification.md`
- `reports/v4_1_hostile_reaudit.md`
- `reports/v4_1_overclaim_fixes_applied.md`
- `reports/v4_1_final_overclaim_scan.md`
- `reports/v4_overclaim_audit.md`

## Files to exclude from a lightweight public release

Recommended exclusions:

- `.git/`
- local execution-cache directories
- `tmp/`
- `output/`, unless specific outputs are intentionally curated
- `__pycache__/` and `*.pyc`
- large raw or processed data files under `data/raw/` and `data/processed/`, unless redistribution is explicitly allowed and desired
- older manuscript drafts not intended as public-facing text, especially if they contain TODOs or outdated claims
- ad hoc intermediate notebook outputs not referenced by the final manuscript

If the repository release includes `data/external/GSE223049/`, note that it is derived from public GEO data and should be distributed only if compatible with the source terms and repository size constraints.

## Suggested GitHub release metadata

Suggested tag:

- `v4.2-final-safe`

Suggested release title:

- `Final safe public preprint v4.2: subtype-dependent LOX-family transcript changes in aging murine thymic stroma`

Short release description:

> Public release of the final cautious computational preprint and reproducibility materials for a GSE240016 thymic stromal scRNA-seq reanalysis focused on subtype-dependent LOX-family transcript changes, with broad sorted-bulk directional comparison in GSE223049. The manuscript is hypothesis-generating and does not claim protein-level, enzymatic, ECM-crosslinking, functional, therapeutic, causal, or subtype-specific external-validation evidence.

## Final release checklist

- [x] Final manuscript exists.
- [x] README points to the final manuscript version.
- [x] Environment files exist.
- [x] Referenced local files exist.
- [x] Referenced final figures exist.
- [x] Referenced supplementary tables exist.
- [x] Final safety check exists.
- [x] No TODOs/placeholders/fake DOI strings detected in the final manuscript.
- [x] No unresolved overclaim wording detected in the final manuscript scan.
- [ ] Confirm whether raw/processed data files should be excluded from the GitHub release because of size or redistribution constraints.
- [ ] Optionally create a clean release branch/tag after final repository review.
