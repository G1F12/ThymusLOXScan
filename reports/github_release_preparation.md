# GitHub release preparation

Prepared branch:

- `release/v4.2-final-safe`

Suggested commit message:

- `Prepare v4.2 final safe public release`

Suggested tag:

- `v4.2-final-safe`

Suggested GitHub release title:

- `Final safe public preprint v4.2`

Suggested GitHub release description:

> Public release of the final cautious computational preprint and reproducibility materials for a GSE240016 thymic stromal scRNA-seq reanalysis focused on subtype-dependent LOX-family transcript changes, with broad sorted-bulk directional comparison in GSE223049. The manuscript is hypothesis-generating and does not claim protein-level, enzymatic, ECM-crosslinking, functional, therapeutic, causal, or subtype-specific external-validation evidence.

## Current branch

Current branch after preparation:

- `release/v4.2-final-safe`

## Repository status summary

Release files were staged. Older drafts and extra reports remain untracked and were not deleted.

Modified tracked files staged for release:

- `.gitignore`
- `README.md`
- `requirements.txt`

New files staged for release:

- `environment.yml`
- `manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md`
- selected release/audit reports
- public result tables
- final, per-sample, and annotation-sanity figures
- release scripts
- supplementary tables

Untracked files intentionally not staged include older manuscript drafts and additional internal reports, such as:

- `manuscript/LOX_thymus_aging_public_preprint.md`
- `manuscript/LOX_thymus_aging_public_preprint_v3_subtype_focused.md`
- `manuscript/LOX_thymus_aging_public_preprint_v4.md`
- `manuscript/LOX_thymus_aging_public_preprint_v4_1_safe.md`
- `manuscript/LOX_thymus_aging_revised.md`
- `manuscript/abstract_revised.md`
- `manuscript/discussion_revised.md`
- `manuscript/methods_detailed.md`
- `results_final_v3.md`
- extra internal reports not required for the release

## Files staged for release

Exact staged release set at preparation time:

```text
.gitignore
README.md
environment.yml
manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md
reports/annotation_uncertainty_check.md
reports/external_dataset_search.md
reports/external_validation_plan.md
reports/fb_composition_adjustment.md
reports/final_public_release_checklist.md
reports/github_release_preparation.md
reports/highest_impact_next_analysis.md
reports/hypothesis_falsification.md
reports/lox_detection_rate_analysis.md
reports/v4_1_hostile_reaudit.md
reports/v4_2_final_safety_check.md
reports/v4_overclaim_audit.md
requirements.txt
results/figures/annotation_sanity/capsFB_marker_detection_heatmap.png
results/figures/annotation_sanity/capsFB_marker_expression_by_sample.png
results/figures/annotation_sanity/capsFB_marker_score_by_sample.png
results/figures/annotation_sanity/intFB_marker_detection_heatmap.png
results/figures/annotation_sanity/intFB_marker_expression_by_sample.png
results/figures/annotation_sanity/intFB_marker_score_by_sample.png
results/figures/annotation_sanity/mTEC1_marker_detection_heatmap.png
results/figures/annotation_sanity/mTEC1_marker_expression_by_sample.png
results/figures/annotation_sanity/mTEC1_marker_score_by_sample.png
results/figures/annotation_sanity/marker_expression_long.tsv
results/figures/annotation_sanity/marker_summary_by_sample.tsv
results/figures/annotation_sanity/medFB_marker_detection_heatmap.png
results/figures/annotation_sanity/medFB_marker_expression_by_sample.png
results/figures/annotation_sanity/medFB_marker_score_by_sample.png
results/figures/annotation_sanity/strict_marker_positive_lox_age_summary.tsv
results/figures/annotation_sanity/strict_marker_positive_lox_by_sample.tsv
results/figures/annotation_sanity/strict_marker_positive_lox_summary.png
results/figures/final/Fig1_pseudobulk_volcano_final.pdf
results/figures/final/Fig1_pseudobulk_volcano_final.png
results/figures/final/Fig2_summary_4panel_final.pdf
results/figures/final/Fig2_summary_4panel_final.png
results/figures/per_sample/lox_per_sample_pseudobulk_expression.pdf
results/figures/per_sample/lox_per_sample_pseudobulk_expression.png
results/figures/per_sample/lox_per_sample_pseudobulk_values.csv
results/figures/per_sample/mtec1_loxl2_per_sample.pdf
results/tables/external_gse223049_lox_validation.tsv
results/tables/external_gse223049_lox_validation_summary.tsv
results/tables/fb_composition_adjusted_lox.tsv
results/tables/lox_detection_rates_by_sample.tsv
results/tables/lox_pseudobulk_complete_results.csv
results/tables/lox_pseudobulk_complete_results.tsv
results/tables/mtec1_loxl2_per_sample_expression.tsv
scripts/annotation_uncertainty_check.py
scripts/external_validation_gse223049_lox.py
scripts/fb_composition_adjusted_lox.py
scripts/figures/plot_final_summary.py
scripts/figures/plot_final_volcano.py
scripts/generate_single_cell_tests.py
scripts/generate_spearman_correlations.py
scripts/lox_detection_rates_by_sample.py
scripts/make_pseudobulk_results_table.py
scripts/plot_per_sample_lox_expression.py
supplementary_tables/README.md
supplementary_tables/Supplementary_Table_1_cell_counts.tsv
supplementary_tables/Supplementary_Table_2_pseudobulk_LOX_results.tsv
supplementary_tables/Supplementary_Table_3_single_cell_tests.tsv
supplementary_tables/Supplementary_Table_4_correlations.tsv
supplementary_tables/make_supplementary_tables.py
```

## Files excluded

Excluded by `.gitignore`:

- `.codex/`
- `.agents/`
- `tmp/`
- `output/`
- `__pycache__/`
- `*.pyc`, `*.pyo`, and related Python bytecode patterns
- `.ipynb_checkpoints/`
- large raw files under `data/raw/`
- large processed files under `data/processed/`
- local downloaded external data under `data/external/`
- local environment folders such as `.venv/`, `venv/`, `env/`, `ENV/`, `.conda/`
- OS/editor junk including `.DS_Store`, `Thumbs.db`, `Desktop.ini`, editor swap files, and backup files

Excluded by not staging:

- older manuscript drafts
- extra internal reports not requested for release
- local raw/processed data files
- local external GEO count downloads
- scratch/output directories

No files were deleted.

## Large files detected

Tracked data/large-file check:

- `data/raw/.gitkeep`
- `data/processed/.gitkeep`
- No tracked files over 50 MB were detected.
- No staged files over 50 MB were detected.

Large local untracked/ignored data files detected:

```text
802886492 data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad
800565374 data/processed/thymus_preprocessed.h5ad
800565374 data/processed/thymus_annotated.h5ad
22088141  data/external/GSE223049/GSE223049_RNA_seq_counts_23_cell_types.txt
5312111   data/external/GSE223049/GSE223049_RNA_seq_counts_23_cell_types.txt.gz
```

Recommendation:

- Do not add raw or processed `.h5ad` files to Git.
- Keep download instructions in `README.md`.
- Keep small derived tables in `results/tables/`.
- Git LFS is not needed for the staged release set because no staged file exceeds 50 MB.
- If future releases require raw or processed AnnData objects, use external archival storage or Git LFS rather than ordinary Git history.

## .gitignore update

`.gitignore` was updated to exclude:

- Codex/local agent state
- scratch/output directories
- Python bytecode/cache files
- Jupyter checkpoints
- raw, processed, and external local data downloads
- local environment folders
- OS/editor junk

## README check

`README.md` now points near the top to:

- `manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md`

It also states:

- "Older manuscript files are working drafts."

## Final safety checks

Final checks passed:

- no TODO in final manuscript;
- no placeholder DOI;
- no fake bioRxiv ID;
- no missing referenced local files;
- README points to final manuscript;
- `reports/v4_2_final_safety_check.md` exists;
- `reports/final_public_release_checklist.md` exists.

The final manuscript safety check confirms:

- no causal claim;
- no protein-level claim;
- no enzymatic activity claim;
- no ECM crosslinking claim;
- no thymic functional claim;
- no therapeutic claim;
- no subtype-specific external validation claim;
- no `possible mediator` language;
- no promotional ranking language.

## Commands to commit, tag, and push

Do not run these until release scope is approved.

Review staged changes:

```bash
git status --short
git diff --cached --stat
git diff --cached --name-status
```

Commit:

```bash
git commit -m "Prepare v4.2 final safe public release"
```

Create annotated tag:

```bash
git tag -a v4.2-final-safe -m "Final safe public preprint v4.2"
```

Push branch and tag:

```bash
git push -u origin release/v4.2-final-safe
git push origin v4.2-final-safe
```

Suggested GitHub release title:

```text
Final safe public preprint v4.2
```

Suggested GitHub release description:

```text
Public release of the final cautious computational preprint and reproducibility materials for a GSE240016 thymic stromal scRNA-seq reanalysis focused on subtype-dependent LOX-family transcript changes, with broad sorted-bulk directional comparison in GSE223049. The manuscript is hypothesis-generating and does not claim protein-level, enzymatic, ECM-crosslinking, functional, therapeutic, causal, or subtype-specific external-validation evidence.
```
