# Lab-facing file hygiene audit

Date: 2026-06-28

## Working tree context

Pre-existing local changes were present before this cleanup, including modified `environment.yml` and `requirements.txt`, untracked older manuscript drafts, untracked older reports, untracked external-validation figure folders, and an untracked downsampling-iterations table. These were not created by this cleanup and are not part of the intended public cleanup staging set.

## Files intended for staging

- `.gitignore`
- `README.md`
- `one_page_summary.md`
- `wet_lab_validation_plan.md`
- `scripts/download_data.py`
- `scripts/dropout_analysis_lox.py`
- `scripts/magic_imputation_lox.py`
- `supplementary_tables/make_supplementary_tables.py`
- `reports/final_public_release_checklist.md`
- `reports/github_release_preparation.md`
- `reports/mtec1_loxl2_dropout_depth_audit.md`
- `reports/v5_2_dropout_protein_update_summary.md`
- `reports/v5_2_file_check.md`
- `reports/lab_facing_ai_trace_audit.md`
- `reports/lab_facing_style_cleanup_audit.md`
- `reports/lab_facing_version_consistency_audit.md`
- `reports/lab_facing_overclaim_safety_audit.md`
- `reports/lab_facing_file_hygiene_audit.md`

## Exclusion checks

- No `.h5ad` files are intended for staging.
- No `data/`, `data/raw/`, `data/processed/`, or `data/external/` files are intended for staging.
- No large raw or processed files are intended for staging.
- Pre-existing modified `environment.yml` and `requirements.txt` are not intended for staging.
- Unrelated old manuscript drafts are not intended for staging.
- v4.2, v5.0, and v5.1 manuscripts are not intended for staging.
- v5.2 remains the current public manuscript/release.

## Final status

Pass, provided staging is restricted to the files listed above and excludes pre-existing unrelated worktree changes.
