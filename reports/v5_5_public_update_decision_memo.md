# v5.5 Public Update Decision Memo

## Scope

This memo decides whether and how to incorporate the E-MTAB-8560 mTEC LOX reanalysis into public project materials. It is a planning memo only. It does not update the manuscript, README, one-page summary, release notes, or repository tags.

## Sources Reviewed

- `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`
- `reports/v5_4_methodological_audit_update_summary.md`
- `reports/mtec1_loxl2_matched_gene_falsification.md`
- `reports/mtec1_loxl2_permutation_test.md`
- `reports/mtec1_loxl2_effect_size_shrinkage_feasibility.md`
- `reports/lox_family_global_fdr_matrix.md`
- `reports/gse240016_sample_metadata_batch_audit.md`
- `reports/gse240016_geo_sra_batch_audit.md`
- `reports/gse240016_geo_sra_batch_audit_update_summary.md`
- `reports/external_emtab8560_biological_unit_audit.md`
- `reports/external_emtab8560_mtec_lox_reanalysis.md`
- `reports/external_emtab8560_vs_gse240016_interpretation.md`
- `reports/v5_5_emtab8560_mtec_update_summary.md`
- `reports/v5_5_emtab8560_mtec_safety_check.md`
- `reports/v5_5_emtab8560_mtec_file_check.md`
- Requested E-MTAB-8560, GSE240016, and LOX-family result tables under `results/tables/`.

## Decision Summary

1. Should the manuscript be updated now?

Recommended choice: **minor cautious update only**.

Rationale: the manuscript can eventually include a short external-context update if reviewed by a human author, but the E-MTAB-8560 result is mixed/inconclusive and should not change the title, central claim, or primary Results framing. No manuscript edit should be made until the patch plan is reviewed.

2. Should README be updated?

Recommended choice: **small "latest analyses" note**.

Rationale: README is the appropriate public place for a brief pointer to the new reports without implying a new release or stronger result.

3. Should one-page summary be updated?

Recommended choice: **one cautious sentence only**.

Rationale: the one-page summary already includes E-MTAB-8560 as TEC/mTEC-like context. If changed, it should add only the mixed/inconclusive mTEC-focused classification.

4. Should a new release/tag be created now?

Recommended choice: **no**.

Rationale: the update is interpretive and planning-oriented. No release or tag should be created before human review and before any agreed public-text patch is made.

5. Should E-MTAB-8560 be included in public project materials?

Recommended placement: **Supplementary/External context**.

Not recommended for main Results. Discussion limitations may mention it briefly if the manuscript is later patched. It should not be held back entirely because the workflow and outputs are already useful context, but the placement must remain cautious.

6. Does E-MTAB-8560 increase confidence in mTEC/TEC Loxl2?

Recommended choice: **mixed/inconclusive**.

Rationale: E-MTAB-8560 provides independent mTEC-focused transcript-level context and shows aged-lower Loxl2 patterns most clearly in mTEClo and cTEC. However, mTEChi support is weaker, combined mTEClo+mTEChi is directional but not strong, global FDR support for Loxl2 rows is not strong, batch-aware modeling was limited/not estimated, and the dataset does not match the GSE240016 mTEC1 label exactly.

## Required Safe Conclusion

E-MTAB-8560 provides independent mTEC-focused transcript-level context, with the clearest aged-lower Loxl2 pattern in mTEClo and cTEC, weaker support in mTEChi, and overall mixed/inconclusive classification. This does not constitute exact GSE240016 mTEC1 validation.

## Safe Wording

- "independent mTEC-focused transcript-level context"
- "aged-lower Loxl2 pattern in mTEClo and cTEC"
- "weaker support in mTEChi"
- "directional but not strong in combined mTEClo+mTEChi"
- "mixed/inconclusive classification"
- "not exact GSE240016 mTEC1 validation"
- "batch-aware modeling was limited/not estimated"
- "cells were not treated as biological replicates"
- "GSE240016 mTEC1 Loxl2 remains a candidate-level transcript signal"
- "dropout/depth and sample or batch confounding remain unresolved concerns"

## Forbidden Wording

Do not use the following as project claims:

- validated
- confirmed
- proved
- robustly replicated
- exact mTEC1 replication
- batch ruled out
- dropout ruled out
- protein-level validation
- mechanism
- driver
- therapeutic
- human conservation
- thymic rejuvenation
- functional restoration

## Recommended Public-Materials Scope

- README: add a small latest-analysis note pointing to the E-MTAB-8560 mTEC reports and stating mixed/inconclusive classification.
- Manuscript: defer editing until review of `reports/v5_5_manuscript_patch_plan.md`; if edited, use one short external-context paragraph and one limitation sentence only.
- One-page summary: add at most one cautious sentence.
- Release/tag: do not create now.

## Decision

Proceed with planning documents only. Do not create a release or tag. Do not strengthen the main manuscript claim. Do not edit public-facing core text until the exact scope is reviewed.
