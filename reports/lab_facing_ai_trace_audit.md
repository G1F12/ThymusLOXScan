# Lab-facing AI trace audit

Date: 2026-06-28

## Files scanned

Primary public-facing files:

- `README.md`
- `one_page_summary.md`
- `wet_lab_validation_plan.md`
- `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`
- `reports/v5_2_dropout_protein_update_summary.md`
- `reports/v5_2_safety_check.md`
- `reports/v5_2_file_check.md`
- `reports/v5_2_outreach_docs_safety_check.md`
- `reports/v5_2_patch_release_file_check.md`
- `reports/mtec1_loxl2_dropout_depth_audit.md`
- `reports/loxl2_protein_spatial_validation_feasibility.md`

Repository-wide text/code scan:

- 156 non-binary paths were scanned after excluding `.git/`, virtual environments, cache folders, `data/`, `.h5ad`, images, PDFs, archives, Word files, and notebooks with embedded image payloads.

## Search terms

Case-insensitive search covered `Codex`, `ChatGPT`, `OpenAI`, `GPT`, `LLM`, AI-generated wording, assistant/user/prompt/conversation wording, usage-limit wording, `.codex`, and `.agents`.

## Matches found

- No matches were found in the primary public-facing file set after cleanup.
- Tracked script/report matches were found before cleanup in:
  - `scripts/dropout_analysis_lox.py`
  - `scripts/magic_imputation_lox.py`
  - `supplementary_tables/make_supplementary_tables.py`
  - `scripts/download_data.py`
  - `reports/final_public_release_checklist.md`
  - `reports/github_release_preparation.md`
- Remaining matches are limited to untracked local `reports/repo_audit.md`, which is not part of the staged public cleanup set.

## Replacements made

- Replaced `.codex/pydeps` runtime fallbacks with neutral `local_pydeps` paths.
- Replaced `.codex/matplotlib` with `local_cache/matplotlib`.
- Replaced "project prompt" wording in `scripts/download_data.py` with neutral dataset-configuration language.
- Replaced release-report mentions of `.codex/`, `.agents/`, and local-agent state with neutral local execution-cache wording.
- Added `local_pydeps/` and `local_cache/` to `.gitignore`.

## Matches intentionally kept

- `.codex/` and `.agents/` remain only in `.gitignore`, where they define ignored local execution-cache directories.
- Untracked `reports/repo_audit.md` was not edited or staged because it is a pre-existing local draft outside the public cleanup commit.

## Final status

Pass for staged/public-facing repository content.
