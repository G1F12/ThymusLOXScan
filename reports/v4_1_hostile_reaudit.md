# v4.1 hostile re-audit

Compared:

- `manuscript/LOX_thymus_aging_public_preprint_v4_1_safe.md`
- `reports/v4_overclaim_audit.md`

## Overall verdict

The serious overclaims from the previous audit were fixed. The v4.1 manuscript no longer presents LOX-family findings as biological ECM remodeling, does not describe `Loxl2` as a possible mediator, does not claim GSE223049 subtype-specific validation, and replaces most confidence-inflating words such as "supported", "survived", "retained", "strongly reduced", and "strongest" with descriptive or sign-based language.

Remaining risky-word hits are mostly acceptable because they occur in prior-literature background, explicit limitations, explicit caveats, file names, or future-validation language. I found two wording issues that are still somewhat unsafe or stylistically risky:

1. "requiring protein, spatial, enzymatic, and functional validation" in the Abstract could still imply that the current candidate is close to functional interpretation.
2. "the most prominent epithelial LOX-family candidate" and "most recurrent LOX-family candidate transcript" are acceptable but slightly promotional; "candidate" language is safer than v4, but could be made even flatter.

## Previous audit completion check

| previous audit area | status in v4.1 | hostile-reviewer comment |
|---|---|---|
| Title used "remodeling" | Fixed | Title now says "transcript changes". |
| Abstract used "remodeling" for current findings | Fixed | Current results now use "transcript differences" and "transcript marker". |
| mTEC1 "strongly reduced" | Fixed | Now says aged-lower estimate and low detection. |
| Composition model "retained" | Fixed | Now says coefficients remained aged-lower, with small-model caveat. |
| Marker-positive "retained" | Fixed | Now says directions had the same sign. |
| GSE223049 "validation" overclaim | Fixed | Now sorted bulk, broad, directional-only, and not subtype-resolved. |
| "possible mediator" | Fixed | Removed. |
| Functional/regeneration priming from prior public dataset | Fixed | Prior dataset sentence no longer says may influence function/regeneration. |
| "reproducible candidate" | Fixed | Now candidate patterns motivating validation/replication. |
| Broad FB "supported" | Fixed | Now "estimated". |
| "supportive of subtype-resolved findings" | Fixed | Now compartment-level context. |
| "primary fibroblast evidence" | Fixed | Now main fibroblast transcript summaries. |
| "subtype-dependent remodeling" in Results | Fixed | Now transcript pattern. |
| mTEC1 heading | Fixed | Heading now includes aged-lower/low detection/n=2 vs n=2. |
| "strongest" ranking | Fixed | Replaced with "most prominent" or "most recurrent". Still slightly promotional but safer. |
| "robustness inspection" | Fixed | Now per-sample inspection. |
| "robustness controls" heading | Fixed | Now sensitivity analyses and unresolved concerns. |
| "driven by fewer cells" | Fixed | Now coincided with detectable-count differences. |
| "survived" | Fixed | Now had the same sign. |
| marker label "supported" | Fixed | Now more consistent with expected labels. |
| "argue against artifact" | Fixed | Now less obvious under filtering strategy. |
| GSE223049 "reproduced" | Fixed in Results | Future-work sentence still uses "reproduce" for future validation, which is acceptable. |
| "substantially lower" | Fixed | Now mean lower. |
| "provides broad external support" | Fixed | Now directionally consistent. |
| correlation biological feature tracking | Fixed | Now weak co-variation with transcript markers. |
| Methods "robust to" | Fixed | Now limited directional sensitivity check. |
| Methods "external validation" | Fixed in prose | File names still contain "validation", acceptable as file paths. |
| Discussion "supports/remodeling" | Fixed | Now suggests transcript differences. |
| fibrosis/tissue-remodeling overinterpretation | Fixed | Now explicitly transcript-level. |
| stromal-state remodeling | Fixed | Now altered proportions or expression levels. |
| "possible mediator" in Discussion | Fixed | Removed. |
| external comparison "strengthens claims" | Fixed | Now adds broad directional context. |
| Limitations "validates" | Fixed | Now "is consistent only with". |
| Data availability "validation analyses" | Fixed | Now sensitivity, annotation-sanity, external-comparison analyses. |
| Figure legends "FDR-significant" without caveat | Fixed | Now small-sample pseudobulk threshold. |
| Figure 2 "subtype-specific" | Fixed | Now selected annotated subtypes. |
| Supplementary "robustness visualizations" | Fixed | Now per-sample descriptive visualizations. |

## Remaining risky-word phrases

### 1. "extracellular matrix remodeling" in Abstract background

**Quote**

> The lysyl oxidase (LOX) family has been implicated in extracellular matrix remodeling and aging in multiple tissues, but its role in thymic stromal aging remains unclear.

**Assessment**

Acceptable. This is prior-literature background, not a claim about the current data. It does not say the current analysis demonstrates ECM remodeling.

**Replacement if extra conservative**

> The lysyl oxidase (LOX) family has been implicated in extracellular matrix biology and aging in multiple tissues, but its role in thymic stromal aging remains unclear.

### 2. "not subtype-resolved validation" in Abstract

**Quote**

> As a sorted bulk RNA-seq, broad-cell-type external comparison, GSE223049 was directionally consistent only with broad aged-lower thymic fibroblast `Lox` and `Loxl2` and broad thymic epithelial `Loxl2`; it was not subtype-resolved validation and could not test medullary fibroblast-specific `Loxl1`.

**Assessment**

Acceptable. The risky word "validation" is explicitly negated and caveated. This sentence is actually protective.

**Replacement**

No replacement needed.

### 3. "functional validation" in Abstract

**Quote**

> Together, these analyses nominate subtype-dependent LOX-family transcript differences as a hypothesis-generating feature of murine thymic stromal aging, with `Loxl2` emerging as the most recurrent candidate transcript marker requiring protein, spatial, enzymatic, and functional validation.

**Assessment**

Mostly acceptable, but still mildly unsafe. "Functional validation" is framed as required future work, not current evidence. However, in the Abstract it may still make the candidate sound closer to a functional hypothesis than the data justify.

**Safer replacement**

> Together, these analyses nominate subtype-dependent LOX-family transcript differences as a hypothesis-generating feature of murine thymic stromal aging, with `Loxl2` emerging as the most recurrent candidate transcript marker requiring orthogonal experimental follow-up.

### 4. "support T cell development" in Introduction

**Quote**

> Age-associated thymic involution is accompanied by reduced thymic output and changes in the thymic stromal microenvironment, including epithelial and mesenchymal compartments that support T cell development [1-4].

**Assessment**

Acceptable. This is general thymus biology background, not a claim that LOX-family transcript changes alter T cell development.

**Replacement**

No replacement needed.

### 5. Prior-literature LOX remodeling and crosslinking background

**Quote**

> The lysyl oxidase family, including LOX and LOX-like enzymes LOXL1-4, encodes secreted copper-dependent amine oxidases involved in extracellular matrix biology, particularly collagen and elastin crosslinking [7,8].

**Assessment**

Acceptable. This is prior biological background about the LOX family, not a current-data claim.

**Replacement**

No replacement needed.

### 6. Prior-literature "tissue remodeling"

**Quote**

> LOX-family genes have also been linked to tissue remodeling, fibrosis, cancer-associated matrix organization, and age-related extracellular matrix changes in multiple contexts [9,10].

**Assessment**

Acceptable, but it is close to the manuscript's risk boundary. It is clearly prior-literature context, but a hostile reviewer may prefer less loaded language.

**Safer replacement**

> LOX-family genes have also been linked to extracellular-matrix-related gene programs, fibrosis, cancer-associated matrix organization, and age-related extracellular matrix changes in multiple contexts [9,10].

### 7. Prior-literature "epithelial or endothelial tissue remodeling"

**Quote**

> Individual family members can have distinct biological associations: LOXL1 is required for elastic fiber homeostasis [11], while LOXL2 has been implicated in collagen IV crosslinking, basement membrane organization, and epithelial or endothelial tissue remodeling [12,13].

**Assessment**

Acceptable as prior-literature background. It does not claim current protein or ECM effects. Still, the sentence is biologically loaded and could prime readers.

**Safer replacement**

> Individual family members can have distinct biological associations: LOXL1 is required for elastic fiber homeostasis [11], while LOXL2 has been implicated in collagen IV crosslinking, basement membrane organization, and epithelial or endothelial biology [12,13].

### 8. "supporting context" in Introduction

**Quote**

> We used pseudobulk differential expression with biological samples as replicates as the primary inferential framework, with descriptive single-cell summaries, per-sample checks, and correlation analyses as supporting context.

**Assessment**

Acceptable. "Supporting context" does not overclaim; it correctly subordinates descriptive analyses to pseudobulk.

**Replacement**

No replacement needed.

### 9. "future experimental validation"

**Quote**

> The goal was not to prove a mechanism, but to identify candidate stromal expression patterns that may motivate future experimental validation and independent subtype-resolved replication.

**Assessment**

Acceptable. "Validation" refers to future work and is paired with "not to prove a mechanism."

**Replacement**

No replacement needed.

### 10. "functional extracellular matrix changes" in Results

**Quote**

> These results show a subtype-dependent LOX-family transcript pattern rather than a single pooled-fibroblast decrease, but they do not establish functional extracellular matrix changes.

**Assessment**

Acceptable. The risky phrase is explicitly negated. This is a useful safeguard.

**Replacement**

No replacement needed.

### 11. "requiring further validation" in mTEC1 Results

**Quote**

> Because n=2 per age group limits inference, this result should be treated as a candidate mTEC1-associated transcript pattern requiring further validation rather than as functional evidence.

**Assessment**

Acceptable. It explicitly limits inference and negates functional evidence.

**Replacement**

No replacement needed.

### 12. "marker support was partial"

**Quote**

> Marker expression was more consistent with the expected medFB and mTEC1 labels than with capsFB and intFB labels, where marker support was partial, especially in aged samples.

**Assessment**

Acceptable. "Support" here refers to marker evidence for labels and is explicitly partial.

**Replacement**

No replacement needed.

### 13. "not subtype-resolved validation" in GSE223049 Results

**Quote**

> Thus, GSE223049 sorted bulk RNA-seq is directionally consistent with aged-lower broad fibroblast `Lox`/`Loxl2` and broad epithelial `Loxl2`, especially `Loxl2`, but it does not test exact fibroblast-subtype or mTEC1-specific claims and is not subtype-resolved validation.

**Assessment**

Acceptable. This is the correct caveat. No subtype-specific external validation is claimed.

**Replacement**

No replacement needed.

### 14. "descriptive supporting summaries" in Methods

**Quote**

> Because these files are derived outputs without a standalone generator script identified in the current repository, they are used only as descriptive supporting summaries.

**Assessment**

Acceptable. This is cautious.

**Replacement**

No replacement needed.

### 15. "not as subtype-resolved validation" in Methods

**Quote**

> This was treated as a descriptive external comparison of broad cell-type directionality, not as subtype-resolved validation.

**Assessment**

Acceptable. The risky word is negated and correctly caveated.

**Replacement**

No replacement needed.

### 16. File names containing "validation"

**Quote**

> `scripts/external_validation_gse223049_lox.py`

**Quote**

> `results/tables/external_gse223049_lox_validation.tsv`

**Quote**

> `results/tables/external_gse223049_lox_validation_summary.tsv`

**Quote**

> `reports/external_validation_plan.md`

**Assessment**

Acceptable. These are existing file names and local paths, not scientific claims. Do not rename just to satisfy prose style unless the repository naming convention is being revised.

**Replacement**

No replacement needed.

### 17. "subtype-resolved validation" in Discussion

**Quote**

> It appears in multiple contexts: broad fibroblast summaries, medullary fibroblast pseudobulk results, mTEC1 low-detection results, detection-rate decomposition, marker-positive annotation sanity checks, and a sorted bulk RNA-seq GSE223049 comparison that is directionally consistent only for broad fibroblast and broad epithelial patterns and is not subtype-resolved validation.

**Assessment**

Acceptable. This is heavily caveated and explicitly denies subtype-resolved validation.

**Replacement**

No replacement needed.

### 18. "requires subtype-resolved validation" in Discussion

**Quote**

> The safest interpretation is that `Loxl2` shows a candidate age-associated transcript difference in annotated mTEC1 cells that requires subtype-resolved validation.

**Assessment**

Acceptable. This calls for future validation and does not claim current validation.

**Replacement**

No replacement needed.

### 19. "functional role" in future work

**Quote**

> Perturbation experiments would be needed before asking whether `Loxl2` or other LOX-family members have any functional role in thymic aging rather than merely marking transcript differences among stromal states.

**Assessment**

Acceptable but mildly risky. It is framed as a future question requiring perturbation experiments, not a current claim. If the goal is maximal conservatism, remove "functional role".

**Safer replacement**

> Perturbation experiments would be needed before asking whether `Loxl2` or other LOX-family members do anything beyond marking transcript differences among stromal states.

### 20. Limitations: functional consequences

**Quote**

> Fourth, transcript abundance does not establish protein abundance, protein localization, secretion, LOX enzymatic activity, extracellular matrix composition, collagen or elastin crosslinking, tissue mechanics, thymic rejuvenation, or functional consequences.

**Assessment**

Acceptable. This is a direct limitation and does not overclaim.

**Replacement**

No replacement needed.

### 21. Limitations: marker support

**Quote**

> The analysis relies on published annotations; marker-positive subset checks reduce but do not eliminate annotation uncertainty, with partial marker support for capsFB and intFB, especially in aged samples.

**Assessment**

Acceptable. It is cautious.

**Replacement**

No replacement needed.

### 22. Reference titles with mechanisms/remodeling

**Quote**

> Cox TR, Erler JT. Remodeling and homeostasis of the extracellular matrix: implications for fibrotic diseases and cancer.

**Quote**

> Hamazaki Y. Adult thymic epithelial cell (TEC) progenitors and TEC stem cells: models and mechanisms for TEC development and maintenance.

**Assessment**

Acceptable. These are reference titles, not manuscript claims.

**Replacement**

No replacement needed.

## Remaining unsafe overclaims

### Unsafe or near-unsafe item A: Abstract "functional validation"

**Quote**

> ...with `Loxl2` emerging as the most recurrent candidate transcript marker requiring protein, spatial, enzymatic, and functional validation.

**Why still risky**

It is a future-validation statement, but in the Abstract it links the candidate marker to functional validation and may make the reader infer a functional hypothesis. This is not a severe overclaim, but it is the one remaining phrase I would soften before public posting.

**Recommended replacement**

> ...with `Loxl2` emerging as the most recurrent candidate transcript marker requiring orthogonal experimental follow-up.

### Unsafe or near-unsafe item B: "most prominent" / "most recurrent"

**Quote**

> The most prominent epithelial LOX-family candidate in this analysis was aged-lower `Loxl2` in mTEC1.

**Quote**

> `Loxl2` is the most recurrent LOX-family candidate transcript in this analysis.

**Why still risky**

These are much safer than "strongest", but they still rank the signal. The ranking is defensible descriptively, yet a hostile reviewer could ask for fully neutral wording.

**Recommended replacement**

> An epithelial LOX-family candidate in this analysis was aged-lower `Loxl2` in mTEC1.

> `Loxl2` appeared repeatedly across the transcript-level summaries in this analysis.

## Final hostile-reviewer conclusion

Every serious overclaim from `reports/v4_overclaim_audit.md` was addressed. v4.1 no longer claims or implies causality, protein-level validation, LOX enzymatic activity change, ECM crosslinking change, thymic functional change, therapeutic relevance, or subtype-specific external validation from GSE223049.

The manuscript is now appropriately cautious for a computational preprint. The only remaining edits I would recommend are stylistic softening of Abstract "functional validation" and optional removal of ranking language such as "most prominent" and "most recurrent".
