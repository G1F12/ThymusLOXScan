# v5.5 Manuscript Patch Plan

This is a patch plan only. It does not modify the manuscript.

## Recommendation

E-MTAB-8560 should be handled as **Supplementary/External context**, not as a main Results finding. The manuscript should receive at most a minor cautious update after human review. No release or tag should be created before that review.

## 1. Exact Sections That Would Need Updating

- Abstract: optional one-clause adjustment only if the manuscript is otherwise patched. Do not strengthen the claim.
- Results, "Cross-dataset external comparisons": optional short external-context sentence or paragraph.
- Discussion, external comparisons paragraph: preferred location for one cautious E-MTAB-8560 update.
- Limitations: required if any E-MTAB-8560 update is added.
- Data and Code Availability: optional addition of the new report/table paths if the manuscript is patched.
- Figure legends: optional only if a supplementary/external figure is referenced.

## 2. Exact Figures/Tables That Could Be Referenced

Tables:

- `results/tables/external_emtab8560_mtec_lox_models.tsv`
- `results/tables/external_emtab8560_mtec_lox_permutation.tsv`
- `results/tables/external_emtab8560_mtec_lox_global_fdr.tsv`
- `results/tables/external_emtab8560_biological_unit_audit.tsv`

Reports:

- `reports/external_emtab8560_biological_unit_audit.md`
- `reports/external_emtab8560_mtec_lox_reanalysis.md`
- `reports/external_emtab8560_vs_gse240016_interpretation.md`

Figures, if present and reviewed:

- `results/figures/external_validation/emtab8560_mtec/`

## 3. Main Text or Supplement/Context

E-MTAB-8560 belongs in supplement/context. It should not be framed as a main Results anchor because the final classification is mixed/inconclusive, global FDR support for Loxl2 rows is not strong, mTEChi support is weaker, and the dataset does not provide an exact GSE240016 mTEC1 counterpart.

## 4. Previous Statements That Should Remain Unchanged

- GSE240016 mTEC1 Loxl2 remains a candidate-level transcript signal.
- The mTEC1 comparison remains n=2 versus n=2 biological samples.
- Exact permutation p-value resolution remains limited.
- Aged mTEC1 samples have lower depth/detected genes.
- Dropout/depth artifact remains partial/unresolved.
- Sample or batch confounding remains unresolved.
- Matched-gene falsification supports prioritization only and does not resolve technical confounding.
- Cells are not biological replicates.
- Public protein/spatial resources provide weak indirect context only.
- No protein, functional, causal, therapeutic, human-conservation, or thymic restoration conclusion is supported.

## 5. Statements That Should Be Softened

- Any prior phrase implying E-MTAB-8560 provides stronger support should be softened to "external context" or "mTEC-focused transcript-level context".
- Any "mTEC-like" language should specify that it is not exact GSE240016 mTEC1.
- Any broad "epithelial Loxl2 support" phrase should note that support is strongest in mTEClo and cTEC, weaker in mTEChi, and mixed/inconclusive overall.
- Any statement that sounds like technical concerns were resolved should be softened to say that concerns remain.

## 6. New Caveats That Must Be Added

- E-MTAB-8560 is not exact GSE240016 mTEC1 validation.
- The E-MTAB-8560 biological unit required SDRF joining; mouse IDs were not directly available in R colData alone.
- Batch and age were partly confounded; batch-aware modeling was limited/not estimated.
- Global FDR support for E-MTAB-8560 Loxl2 rows is not strong.
- mTEChi support is weaker or absent compared with mTEClo and cTEC.
- No p-values should be combined across datasets.

## 7. Title/Abstract

Title should not change.

Abstract should change only if the manuscript receives a reviewed minor update. If changed, it should use one cautious phrase such as "E-MTAB-8560 provided mixed/inconclusive mTEC-focused transcript-level context" and should not alter the central conclusion.

## 8. Main Claim

The main claim should not change and should not be strengthened. The safest main claim remains that subtype-dependent LOX-family transcript differences in murine thymic stroma are hypothesis-generating and require orthogonal follow-up.

## Proposed Patch Scope

Recommended patch scope after review:

- Add a small README note.
- Add one optional Discussion paragraph.
- Add one required Limitations sentence or paragraph if the Discussion is updated.
- Optionally add report/table paths to Data and Code Availability.
- Do not add a new main figure.
- Do not create a release/tag until human review approves the final wording.
