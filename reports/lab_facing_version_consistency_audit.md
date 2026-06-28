# Lab-facing version consistency audit

Date: 2026-06-28

## Required public state

- Current manuscript: `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`
- Current release: `v5.2-dropout-protein-feasibility`
- Previous external-validation update: `v5.0-external-validation`
- Previous stable cautious release: `v4.2-final-safe`
- v5.1 is historical/intermediate, not current.

## Files checked

- `README.md`
- `one_page_summary.md`
- `wet_lab_validation_plan.md`
- `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`
- v5.2 release and safety reports under `reports/`

## Findings

- `README.md` previously pointed to `manuscript/LOX_thymus_aging_public_preprint_v5_1_external_validation_safe.md` as current and emphasized `v4.2-final-safe`.
- `one_page_summary.md` already pointed to the v5.2 manuscript.
- `wet_lab_validation_plan.md` already pointed to the v5.2 manuscript.
- Historical v5.1, v5.0, and v4.2 references in v5.2 reports were acceptable when clearly describing manuscript lineage or previous releases.

## Replacements made

- Updated `README.md` to identify the v5.2 manuscript as current.
- Updated `README.md` to identify `v5.2-dropout-protein-feasibility` as the current release.
- Added previous-release context for `v5.0-external-validation` and `v4.2-final-safe`.
- Updated the mTEC1 dropout/depth report recommendation from v5.1 wording to v5.2 wording.

## Final checks

- README points to v5.2: yes.
- One-page summary points to v5.2: yes.
- Wet-lab validation plan points to v5.2: yes.
- README does not imply v4.2, v5.0, or v5.1 is current: yes.

## Final status

Pass.
