# v5.5 Public Update Proposed Wording

This file proposes text blocks only. It does not edit the README, manuscript, one-page summary, release notes, or tags.

## 1. README Latest Analysis Paragraph

Status: **recommended**

Suggested placement: after "Current Manuscript And Release" or within "Key Outputs".

Proposed text:

> Latest analysis note: E-MTAB-8560 MouseThymusAgeing was reanalyzed as independent mTEC-focused transcript-level context for LOX-family patterns. The E-MTAB-8560 Loxl2 pattern is clearest as aged-lower in mTEClo and cTEC, weaker in mTEChi, and mixed/inconclusive overall; it is not exact GSE240016 mTEC1 validation. See `reports/external_emtab8560_biological_unit_audit.md`, `reports/external_emtab8560_mtec_lox_reanalysis.md`, and `reports/external_emtab8560_vs_gse240016_interpretation.md`.

## 2. Manuscript Results Paragraph

Status: **not recommended**

Rationale: the main Results should not be expanded now because E-MTAB-8560 is mixed/inconclusive, not exact mTEC1 validation, and should remain external context unless a reviewed manuscript patch is approved.

If a short external-context paragraph is later required, use only a cautious version:

> In an independent E-MTAB-8560 MouseThymusAgeing reanalysis, per-mouse summaries provided mTEC-focused transcript-level context for LOX-family genes. Loxl2 showed the clearest aged-lower pattern in mTEClo and cTEC, weaker support in mTEChi, and directional but not strong support in combined mTEClo+mTEChi. The overall E-MTAB-8560 classification was mixed/inconclusive, with limited batch-aware modeling and no exact match to the GSE240016 mTEC1 label.

## 3. Manuscript Discussion Paragraph

Status: **optional**

Suggested placement: external comparisons paragraph in Discussion, not central Results.

Proposed text:

> The updated E-MTAB-8560 mTEC-focused reanalysis is best interpreted as transcript-level external context rather than a direct test of the GSE240016 mTEC1 result. E-MTAB-8560 showed aged-lower Loxl2 patterns most clearly in mTEClo and cTEC, weaker support in mTEChi, and mixed/inconclusive support overall. These data slightly broaden the epithelial context in which Loxl2 can be prioritized for follow-up, but they do not change the candidate-level status of the GSE240016 mTEC1 observation.

## 4. Manuscript Limitations Paragraph

Status: **recommended if the manuscript is patched**

Suggested placement: Limitations.

Proposed text:

> The E-MTAB-8560 reanalysis does not resolve the main limitations of the GSE240016 mTEC1 Loxl2 finding. GSE240016 remains limited by n=2 versus n=2 biological samples, sparse Loxl2 detection, lower aged mTEC1 depth, limited exact-permutation resolution, and unresolved sample or batch confounding. E-MTAB-8560 uses a different age series and annotation framework, lacks an exact GSE240016 mTEC1 counterpart, and had limited/not estimated batch-aware modeling; cells were not treated as biological replicates.

## 5. One-Page Summary Optional Sentence

Status: **recommended**

Suggested placement: under "External context".

Proposed text:

> Updated E-MTAB-8560 mTEC-focused reanalysis remains mixed/inconclusive, with the clearest aged-lower Loxl2 pattern in mTEClo and cTEC and weaker support in mTEChi.

## 6. Future Release Notes Draft

Status: **optional for a later reviewed release; not recommended now**

Use only if a future human-reviewed public-text update is approved.

Proposed text:

> Added a cautious v5.5 public update incorporating E-MTAB-8560 mTEC-focused transcript-level context. The update keeps the GSE240016 mTEC1 Loxl2 result at candidate level, reports E-MTAB-8560 as mixed/inconclusive, and does not change the main manuscript claim.
