# v5.1 public trace cleanup

## Files scanned

- `README.md`
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md`
- `reports/gse223049_external_validation_v2.md`
- `reports/emtab8560_tec_external_validation.md`
- `reports/gse231906_human_external_validation.md`
- `reports/gse231906_human_metadata_only_plan.md`
- `reports/external_dataset_search_v2.md`
- `reports/external_validation_priority_matrix.md`
- `reports/cross_dataset_lox_analysis.md`
- `reports/cross_dataset_figure_notes.md`
- `reports/v5_file_and_reproducibility_check.md`
- `reports/v5_hostile_overclaim_audit.md`
- `reports/v5_1_overclaim_fixes_applied.md`
- `reports/v5_1_safety_check.md`
- `reports/v5_1_file_and_safety_final_check.md`
- `reports/v5_release_preparation_checklist.md`
- `scripts/external_validation_gse223049_lox.py`
- `scripts/external_validation_emtab8560_tec_lox.py`
- `scripts/external_validation_gse231906_human_lox.py`
- `scripts/build_cross_dataset_lox_matrix.py`
- `scripts/figures/plot_cross_dataset_lox_validation.py`

## Risky strings found

None found in release-candidate files after cleanup scan.

## Replacements made

- Updated `README.md` so v5.1 is the current external-validation manuscript and v4.2 is the previous stable cautious release.
- Updated `reports/v5_release_preparation_checklist.md` to reflect the README v5.1 pointer.
- No public workflow trace replacements were required because none were present in the scanned release files.

## Strings intentionally kept and why

- `README.md:76`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:9`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:15`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:31`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:49`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:103`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:111`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:135`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:149`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:153`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:155`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:190`: `model` kept as scientific/statistical use.
- `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md:195`: `model` kept as scientific/statistical use.
- `reports/emtab8560_tec_external_validation.md:24`: `model` kept as scientific/statistical use.
- `reports/gse231906_human_external_validation.md:44`: `model` kept as scientific/statistical use.
- `scripts/external_validation_emtab8560_tec_lox.py:455`: `model` kept as scientific/statistical use.
- `scripts/external_validation_gse231906_human_lox.py:494`: `model` kept as scientific/statistical use.
- `scripts/external_validation_gse231906_human_lox.py:561`: `model` kept as scientific/statistical use.
- `scripts/external_validation_gse231906_human_lox.py:566`: `model` kept as scientific/statistical use.
- `scripts/external_validation_gse231906_human_lox.py:596`: `model` kept as scientific/statistical use.

## Assessment

The release-candidate files do not contain public-facing workflow traces from the configured scan list. Legitimate scientific and statistical uses of `model` were retained.