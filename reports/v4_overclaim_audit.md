# v4 overclaim audit

Hostile-reviewer audit of `manuscript/LOX_thymus_aging_public_preprint_v4.md`.

## Overall assessment

The manuscript is generally cautious, but several sentences still use language that can sound stronger than the data support. The main recurring issues are:

- "remodeling" can imply biological ECM or tissue remodeling rather than transcript differences.
- "strongly reduced", "strongest", "supported", "survived", "retained", and "validation" can imply more confidence than a small-n computational reanalysis warrants.
- "possible mediator" implies mechanistic involvement that is not tested.
- External GSE223049 evidence is broad sorted bulk RNA-seq and should not be framed as validation of subtype-specific results.

## Flagged sentences

### 1. Title implies biological remodeling

**Quoted sentence**

> Subtype-dependent LOX-family remodeling in aging murine thymic stroma

**Problem**

"Remodeling" can be read as biological stromal or ECM remodeling, not just transcript-level differences. The study does not measure protein, enzymatic activity, ECM structure, tissue mechanics, or function.

**Safer replacement**

> Subtype-dependent LOX-family transcript changes in aging murine thymic stroma

### 2. Abstract overuses remodeling language

**Quoted sentence**

> Pseudobulk and sample-level summaries indicated subtype-dependent LOX-family remodeling rather than a uniform fibroblast-wide response.

**Problem**

"Indicated" plus "remodeling" may imply a demonstrated biological process. The analyses show transcript abundance differences in annotated cell groups.

**Safer replacement**

> Pseudobulk and sample-level summaries showed subtype-dependent LOX-family transcript differences rather than a uniform fibroblast-wide pattern.

### 3. mTEC1 effect size language is too strong

**Quoted sentence**

> In the mTEC1 epithelial compartment, `Loxl2` was strongly reduced and detected in very few aged cells, but this comparison included only two young and two aged biological samples.

**Problem**

"Strongly reduced" overemphasizes the estimate despite n=2 vs n=2 and low detection. The caveat helps, but the first clause still sounds too confident.

**Safer replacement**

> In the mTEC1 epithelial compartment, `Loxl2` showed an aged-lower pseudobulk estimate and very low aged-cell detection, but this comparison included only two young and two aged biological samples.

### 4. "Retained" may imply robustness beyond what the model supports

**Quoted sentence**

> Broad fibroblast sensitivity analysis retained aged-lower `Lox`, `Loxl1`, and `Loxl2` directions after subtype-composition adjustment, although the adjusted model had only six samples and one residual degree of freedom.

**Problem**

"Retained" sounds like a robustness pass. With one residual degree of freedom, the result is only a sign check.

**Safer replacement**

> In a very small fibroblast composition-adjusted sensitivity model, the age coefficients for `Lox`, `Loxl1`, and `Loxl2` remained aged-lower, but the model had only six samples and one residual degree of freedom.

### 5. Marker-positive subsets framed too strongly

**Quoted sentence**

> Annotation sanity checks showed that capsFB, intFB, medFB, and mTEC1 labels expressed expected markers to varying degrees, and all four focal directions were retained in stricter marker-positive subsets.

**Problem**

"Showed" and "retained" imply annotation robustness. Marker-positive subset checks reduce one concern but do not validate labels, especially for partial capsFB/intFB marker support.

**Safer replacement**

> Annotation sanity checks found expected marker expression to varying degrees in capsFB, intFB, medFB, and mTEC1 labels, and the four focal transcript directions had the same sign in stricter marker-positive subsets.

### 6. External validation language can be read as too definitive

**Quoted sentence**

> Independent sorted bulk RNA-seq from GSE223049 supported broad aged-lower directionality for thymic fibroblast `Lox` and `Loxl2` and thymic epithelial `Loxl2`, while it did not validate the subtype-specific medullary fibroblast `Loxl1` increase.

**Problem**

"Supported" is acceptable but "validate" in the same sentence risks suggesting true validation is possible from this broad sorted bulk dataset. The dataset is descriptive, bulk, and not subtype-resolved.

**Safer replacement**

> Independent sorted bulk RNA-seq from GSE223049 showed the same broad aged-lower direction for thymic fibroblast `Lox` and `Loxl2` and thymic epithelial `Loxl2`, but it could not test medullary fibroblast-specific `Loxl1`.

### 7. Candidate mediator language overreaches

**Quoted sentence**

> Together, these analyses identify subtype-dependent LOX-family transcript remodeling as a cautious candidate feature of murine thymic stromal aging, with `Loxl2` emerging as the strongest candidate marker and possible mediator requiring protein, spatial, enzymatic, and functional validation.

**Problem**

"Possible mediator" implies a mechanistic role. No perturbation, protein, enzyme, ECM, or functional data are shown. "Strongest candidate" also sounds like ranking confidence beyond the evidence.

**Safer replacement**

> Together, these analyses nominate subtype-dependent LOX-family transcript differences as a hypothesis-generating feature of murine thymic stromal aging, with `Loxl2` emerging as the most recurrent candidate transcript marker requiring protein, spatial, enzymatic, and functional validation.

### 8. Introduction invokes thymic functional change near the study rationale

**Quoted sentence**

> Recent single-cell and spatial transcriptomic analyses of murine thymic aging identified age-associated epithelial states and stromal changes that may influence thymic function and regeneration [5,6].

**Problem**

This sentence refers to prior work, but in the current manuscript it may prime readers to infer functional implications from the LOX-family reanalysis. The present study does not test function or regeneration.

**Safer replacement**

> Recent single-cell and spatial transcriptomic analyses of murine thymic aging identified age-associated epithelial states and stromal changes in public datasets [5,6].

### 9. "Reproducible candidate" is too strong for one main dataset plus broad bulk support

**Quoted sentence**

> The goal was not to prove a mechanism, but to identify reproducible candidate stromal expression patterns that may motivate future experimental validation.

**Problem**

"Reproducible" is too strong before independent subtype-resolved replication. GSE223049 supports only broad directions.

**Safer replacement**

> The goal was not to prove a mechanism, but to identify candidate stromal expression patterns that may motivate future experimental validation and independent subtype-resolved replication.

### 10. "Supported lower expression" overstates broad fibroblast inference

**Quoted sentence**

> At the broad fibroblast level, pseudobulk DESeq2 supported lower expression of several LOX-family genes in aged samples.

**Problem**

"Supported" can sound inferentially strong despite n=3 vs n=3 and known subtype-composition sensitivity.

**Safer replacement**

> At the broad fibroblast level, pseudobulk DESeq2 estimated lower expression of several LOX-family genes in aged samples.

### 11. "Supportive of subtype-resolved findings" can overstate pooled fibroblast value

**Quoted sentence**

> Thus, reduced broad fibroblast expression is best presented as partly influenced by subtype composition and supportive of, rather than primary to, the subtype-resolved findings.

**Problem**

"Supportive of" still implies that the broad fibroblast result strengthens the subtype findings. Because broad effects can be composition artifacts, it should be presented as context.

**Safer replacement**

> Thus, reduced broad fibroblast expression is best presented as partly influenced by subtype composition and as compartment-level context rather than primary evidence for the subtype-resolved findings.

### 12. "Primary fibroblast evidence" is too confident for n=3 vs n=3

**Quoted sentence**

> Subtype-resolved pseudobulk analysis provided the primary fibroblast evidence and indicated that LOX-family aging patterns were not uniform across fibroblast populations.

**Problem**

"Evidence" and "indicated" sound strong. The analysis is still small-n and annotation-dependent.

**Safer replacement**

> Subtype-resolved pseudobulk analysis provided the main fibroblast transcript summaries and suggested that LOX-family age-associated patterns were not uniform across annotated fibroblast populations.

### 13. "Support subtype-dependent remodeling" overstates biology

**Quoted sentence**

> These results support subtype-dependent LOX-family remodeling rather than a single pooled-fibroblast decrease, but they do not establish functional extracellular matrix remodeling.

**Problem**

The second clause is cautious, but "support subtype-dependent LOX-family remodeling" still implies biology beyond transcript differences.

**Safer replacement**

> These results support a subtype-dependent LOX-family transcript pattern rather than a single pooled-fibroblast decrease, but they do not establish functional extracellular matrix remodeling.

### 14. Section heading overstates mTEC1 result

**Quoted sentence**

> Loxl2 is strongly reduced in mTEC1, with limited replicate support

**Problem**

The heading foregrounds "strongly reduced" and softens the caveat. For n=2 vs n=2 and near-absence/detection-rate behavior, the heading should be more descriptive.

**Safer replacement**

> mTEC1 Loxl2 shows aged-lower expression and low aged detection in a n=2 vs n=2 comparison

### 15. "Strongest subtype-level epithelial change" overstates ranking

**Quoted sentence**

> The strongest subtype-level epithelial LOX-family change was observed for `Loxl2` in mTEC1.

**Problem**

"Strongest" ranks an unstable low-count result and may invite overinterpretation.

**Safer replacement**

> The most prominent epithelial LOX-family candidate in this analysis was aged-lower `Loxl2` in mTEC1.

### 16. Candidate mTEC1 aging marker language is still a bit strong

**Quoted sentence**

> Because n=2 per age group limits inference, this result should be treated as a candidate mTEC1 aging marker requiring validation rather than as evidence of a proven functional mechanism.

**Problem**

"Candidate mTEC1 aging marker" may imply subtype specificity and biological marker status. The finding is a low-detection transcript pattern in a small comparison.

**Safer replacement**

> Because n=2 per age group limits inference, this result should be treated as a candidate mTEC1-associated transcript pattern requiring validation rather than as evidence of a functional mechanism.

### 17. "Robustness inspection" sounds too reassuring

**Quoted sentence**

> Per-sample robustness inspection showed that both young mTEC1 pseudobulk samples had higher normalized `Loxl2` expression than both aged mTEC1 samples.

**Problem**

"Robustness inspection" may imply robustness despite only four total samples.

**Safer replacement**

> Per-sample inspection showed that both young mTEC1 pseudobulk samples had higher normalized `Loxl2` expression than both aged mTEC1 samples.

### 18. "Robustness analyses" heading may overstate sensitivity checks

**Quoted sentence**

> Robustness controls and sensitivity analyses

**Problem**

"Robustness controls" can suggest the concerns are controlled. Many concerns remain unresolved.

**Safer replacement**

> Sensitivity analyses and unresolved robustness concerns

### 19. "Model showed" too strong for one residual degree of freedom

**Quoted sentence**

> A fibroblast subtype-composition sensitivity model showed that broad fibroblast aged-lower directions for `Lox`, `Loxl1`, and `Loxl2` were retained after adjustment for capsFB, intFB, medFB, and Fat fractions.

**Problem**

With one residual degree of freedom, the model cannot meaningfully establish adjusted effects. "Showed" and "retained" are too strong.

**Safer replacement**

> In a fibroblast subtype-composition sensitivity model, broad fibroblast age coefficients for `Lox`, `Loxl1`, and `Loxl2` remained aged-lower after adjustment for capsFB, intFB, medFB, and Fat fractions.

### 20. "Driven" implies mechanism of the observed signal

**Quoted sentence**

> Detection-rate decomposition indicated that several signals were driven partly or largely by fewer cells expressing the relevant gene.

**Problem**

"Driven by" can imply a resolved causal decomposition. Detection-rate summaries are descriptive and cannot separate true biological absence from dropout or depth.

**Safer replacement**

> Detection-rate decomposition suggested that several transcript differences coincided partly or largely with fewer cells having detectable counts for the relevant gene.

### 21. "Survived" and "reduce likelihood" overstate internal checks

**Quoted sentence**

> The highlighted directions survived quality-control threshold checks, raw-count CPM summaries, normalized-expression summaries, exclusion of sorted/enriched sample preparations, and leave-one-sample-out direction checks where feasible.

**Problem**

"Survived" sounds like rigorous falsification. With n=2 or n=3 per group, these are limited direction checks.

**Safer replacement**

> The highlighted directions had the same sign in quality-control threshold checks, raw-count CPM summaries, normalized-expression summaries, exclusion of sorted/enriched sample preparations, and leave-one-sample-out direction checks where feasible.

### 22. "Reduce likelihood" is too reassuring

**Quoted sentence**

> These controls reduce the likelihood that the results are artifacts of one normalization representation, one QC filter, or one obvious outlier sample.

**Problem**

This may overstate the strength of small-n checks. They show consistency across summaries but cannot estimate artifact probability.

**Safer replacement**

> These checks show that the reported directions are not unique to one normalization representation, one QC filter, or one obvious outlier sample.

### 23. Marker expression "supported label consistency" overstates annotation check

**Quoted sentence**

> Marker expression broadly supported medFB and mTEC1 label consistency, while capsFB and intFB marker support was partial, especially in aged samples.

**Problem**

"Supported label consistency" can sound like annotation validation. Marker expression is only a sanity check.

**Safer replacement**

> Marker expression was more consistent with the expected medFB and mTEC1 labels than with capsFB and intFB labels, where marker support was partial, especially in aged samples.

### 24. "Argue against artifact explanation" is too strong

**Quoted sentence**

> These findings argue against the simplest marker-threshold artifact explanation, but they do not prove that the published annotations are correct or that annotation uncertainty is fully resolved.

**Problem**

"Argue against" is acceptable but still a bit adversarially strong; the check only shows same-sign summaries after one marker-positive filtering strategy.

**Safer replacement**

> These findings make a simple marker-threshold artifact less obvious under this filtering strategy, but they do not prove that the published annotations are correct or that annotation uncertainty is fully resolved.

### 25. External validation "reproduced" overstates what bulk RNA-seq can do

**Quoted sentence**

> It can nevertheless ask whether broad thymic fibroblast and epithelial LOX-family directions are reproduced outside the original single-cell dataset.

**Problem**

"Reproduced" implies replication. The external dataset differs in age, sorting, technology, and resolution and is summarized descriptively.

**Safer replacement**

> It can nevertheless provide a descriptive check of whether broad thymic fibroblast and epithelial LOX-family directions are similar outside the original single-cell dataset.

### 26. "Substantially lower" overstates descriptive external result

**Quoted sentence**

> In GSE223049 thymic fibroblasts, `Lox` was lower in aged samples, with delta log2(CPM+1) = -0.354, and `Loxl2` was substantially lower, with delta log2(CPM+1) = -1.206.

**Problem**

"Substantially lower" is interpretive and may imply statistical testing or biological magnitude. This external analysis is descriptive.

**Safer replacement**

> In GSE223049 thymic fibroblasts, mean `Lox` was lower in aged samples, with delta log2(CPM+1) = -0.354, and mean `Loxl2` was lower, with delta log2(CPM+1) = -1.206.

### 27. "Provides broad external support" is too strong

**Quoted sentence**

> Thus, GSE223049 provides broad external support for fibroblast `Lox`/`Loxl2` and epithelial `Loxl2` aged-lower directionality, especially `Loxl2`, but not for exact fibroblast-subtype or mTEC1-specific claims.

**Problem**

"Provides broad external support" may be acceptable, but a hostile reviewer could read it as validation. Because this is descriptive sorted bulk RNA-seq, "is directionally consistent" is safer.

**Safer replacement**

> Thus, GSE223049 is directionally consistent with aged-lower broad fibroblast `Lox`/`Loxl2` and broad epithelial `Loxl2`, especially `Loxl2`, but it does not test exact fibroblast-subtype or mTEC1-specific claims.

### 28. Correlation interpretation could imply biological feature tracking

**Quoted sentence**

> These descriptive correlations are consistent with `Loxl2` varying with structural fibroblast features more than with a `Snai1`-associated EMT-like signal, but they do not establish pathway direction or protein-level co-regulation.

**Problem**

"Structural fibroblast features" is biologically loaded and not directly measured; only transcript correlations with marker genes were measured.

**Safer replacement**

> These descriptive correlations are consistent with `Loxl2` co-varying weakly with selected structural or mesenchymal transcript markers more than with `Snai1`, but they do not establish pathway direction or protein-level co-regulation.

### 29. Methods uses "robust" for an underpowered adjustment

**Quoted sentence**

> To evaluate whether broad fibroblast LOX-family age effects were robust to fibroblast subtype mixture, sample-level broad fibroblast pseudobulk expression was modeled with age while adjusting for capsFB, intFB, medFB, and Fat fractions within each biological sample.

**Problem**

"Robust to" is too strong because the model has six samples and one residual degree of freedom.

**Safer replacement**

> To perform a limited directional sensitivity check for fibroblast subtype mixture, sample-level broad fibroblast pseudobulk expression was modeled with age while adjusting for capsFB, intFB, medFB, and Fat fractions within each biological sample.

### 30. Methods labels GSE223049 as validation

**Quoted sentence**

> Independent external validation used GSE223049, a mouse sorted bulk RNA-seq dataset with 2-month and 22-24-month samples across sorted cell types.

**Problem**

"Validation" overstates the descriptive broad-cell-type check.

**Safer replacement**

> Independent external comparison used GSE223049, a mouse sorted bulk RNA-seq dataset with 2-month and 22-24-month samples across sorted cell types.

### 31. Methods repeats validation terminology

**Quoted sentence**

> This was treated as descriptive external validation of broad cell-type directionality, not as subtype-resolved validation.

**Problem**

Even with "descriptive", "validation" may be too strong. The data are an external directionality comparison.

**Safer replacement**

> This was treated as a descriptive external comparison of broad cell-type directionality, not as subtype-resolved validation.

### 32. Discussion conclusion uses "supports"

**Quoted sentence**

> This computational reanalysis supports a narrow central conclusion: age-associated LOX-family remodeling in thymic stroma is subtype-dependent and partly externally supported, especially for `Loxl2` directionality.

**Problem**

"Supports" and "remodeling" imply more than transcript-level associations. "Partly externally supported" may also overstate GSE223049.

**Safer replacement**

> This computational reanalysis suggests a narrow central conclusion: age-associated LOX-family transcript differences in thymic stroma are subtype-dependent, with broad external directional consistency for selected `Loxl2` comparisons.

### 33. Canonical fibrosis contrast risks biological overinterpretation

**Quoted sentence**

> The overall pattern differs from a simple canonical fibrosis-like model in which LOX-family expression uniformly increases with age or tissue remodeling.

**Problem**

This is a useful discussion point, but it may imply the current dataset tests a fibrosis model or tissue remodeling. It only assesses transcripts in thymic stromal cells.

**Safer replacement**

> At the transcript level, the observed pattern does not resemble a simple model in which all LOX-family genes increase uniformly across stromal compartments with age.

### 34. Stromal-state remodeling language implies biology

**Quoted sentence**

> This suggests that thymic aging may involve loss, gain, or remodeling of specific LOX-positive stromal states rather than a single monotonic LOX-family response across all fibroblasts.

**Problem**

"May involve loss, gain, or remodeling" implies cellular-state changes that are not directly resolved. Composition and annotation artifacts remain possible.

**Safer replacement**

> This suggests that the annotated aged thymic stroma may contain altered proportions or expression levels of LOX-positive stromal states rather than a single monotonic LOX-family response across all fibroblasts.

### 35. "Strongest cautious candidate" and "appears in multiple contexts" overstate convergence

**Quoted sentence**

> `Loxl2` is the strongest cautious candidate signal in this analysis.

**Problem**

"Strongest" is subjective and may sound like a ranked biological conclusion.

**Safer replacement**

> `Loxl2` is the most recurrent LOX-family candidate transcript in this analysis.

### 36. "Independent validation" in context list overstates GSE223049

**Quoted sentence**

> It appears in multiple contexts: broad fibroblast summaries, medullary fibroblast pseudobulk results, mTEC1 low-detection results, detection-rate decomposition, marker-positive annotation sanity checks, and independent GSE223049 broad fibroblast and epithelial validation.

**Problem**

"Independent GSE223049 ... validation" overstates the broad sorted bulk comparison and may imply validation of medFB/mTEC1.

**Safer replacement**

> It appears in multiple contexts: broad fibroblast summaries, medullary fibroblast pseudobulk results, mTEC1 low-detection results, detection-rate decomposition, marker-positive annotation sanity checks, and a broad external GSE223049 fibroblast and epithelial directionality comparison.

### 37. "Possible mediator" is mechanistic overclaim

**Quoted sentence**

> These convergent observations make `Loxl2` a candidate marker and possible mediator of age-associated thymic stromal remodeling.

**Problem**

"Possible mediator" implies causal mechanism. "Convergent" and "remodeling" further strengthen the biological implication.

**Safer replacement**

> These observations make `Loxl2` a recurrent candidate transcript marker of age-associated thymic stromal differences.

### 38. "Strong" mTEC1 direction is too much

**Quoted sentence**

> The direction was strong, appeared in per-sample summaries, and survived stricter marker-positive filtering, while GSE223049 supported a broad epithelial aged-lower `Loxl2` direction.

**Problem**

"Strong", "survived", and "supported" compound confidence in a fragile n=2 vs n=2 low-detection result, and GSE223049 is broad epithelial only.

**Safer replacement**

> The direction was aged-lower in per-sample summaries and had the same sign after stricter marker-positive filtering, while GSE223049 showed a broad epithelial aged-lower `Loxl2` direction.

### 39. "Marks a candidate epithelial state change" may overstate state biology

**Quoted sentence**

> The safest interpretation is that `Loxl2` marks a candidate age-associated epithelial state change that requires subtype-resolved validation.

**Problem**

"Marks" and "state change" imply a biological epithelial state not directly validated. The result may reflect dropout, annotation, or sampling.

**Safer replacement**

> The safest interpretation is that `Loxl2` shows a candidate age-associated transcript difference in annotated mTEC1 cells that requires subtype-resolved validation.

### 40. "External validation strengthens claims" overstates comparison

**Quoted sentence**

> The external validation strengthens only the broadest claims.

**Problem**

"Validation" and "strengthens claims" can sound too assertive for descriptive sorted bulk RNA-seq.

**Safer replacement**

> The external comparison only adds broad directional context.

### 41. "Supports" repeated in discussion

**Quoted sentence**

> GSE223049 supports aged-lower thymic fibroblast `Lox` and `Loxl2` and aged-lower thymic epithelial `Loxl2`, but it is sorted bulk RNA-seq and cannot test capsFB, medFB, or mTEC1 specificity.

**Problem**

"Supports" is stronger than necessary and may imply validation.

**Safer replacement**

> GSE223049 shows aged-lower mean thymic fibroblast `Lox` and `Loxl2` and aged-lower mean thymic epithelial `Loxl2`, but it is sorted bulk RNA-seq and cannot test capsFB, medFB, or mTEC1 specificity.

### 42. "Independent support is strongest" is too confident

**Quoted sentence**

> Thus, independent support is strongest for `Loxl2` directionality and weaker for exact subtype assignments.

**Problem**

"Independent support is strongest" may over-rank external evidence. There is only one broad sorted bulk external dataset.

**Safer replacement**

> Thus, the external comparison is most directionally consistent for broad `Loxl2` patterns and does not establish exact subtype assignments.

### 43. Perturbation sentence implies plausible contribution

**Quoted sentence**

> Perturbation experiments would be needed to determine whether `Loxl2` or other LOX-family members contribute to thymic aging rather than merely marking altered stromal states.

**Problem**

This is mostly cautious, but "contribute to thymic aging" can sound like the manuscript has raised a therapeutic or mechanistic hypothesis. It is acceptable in future-work framing, but safer wording is possible.

**Safer replacement**

> Perturbation experiments would be needed before asking whether `Loxl2` or other LOX-family members have any functional role in thymic aging rather than merely marking transcript differences among stromal states.

### 44. Limitation says GSE223049 "validates"

**Quoted sentence**

> Second, the independent support from GSE223049 comes from sorted bulk RNA-seq rather than subtype-resolved single-cell data, so it validates only broad fibroblast and epithelial directionality.

**Problem**

"Validates" is too strong even with "only broad". Use "is consistent with".

**Safer replacement**

> Second, the independent comparison from GSE223049 comes from sorted bulk RNA-seq rather than subtype-resolved single-cell data, so it is consistent only with broad fibroblast and epithelial directionality.

### 45. Data availability labels scripts as validation analyses

**Quoted sentence**

> New robustness and validation analyses for this version include `scripts/fb_composition_adjusted_lox.py`, `scripts/lox_detection_rates_by_sample.py`, `scripts/annotation_uncertainty_check.py`, and `scripts/external_validation_gse223049_lox.py`; reports in `reports/fb_composition_adjustment.md`, `reports/lox_detection_rate_analysis.md`, `reports/annotation_uncertainty_check.md`, `reports/external_dataset_search.md`, `reports/external_validation_plan.md`, `reports/highest_impact_next_analysis.md`, and `reports/hypothesis_falsification.md`; and tables in `results/tables/fb_composition_adjusted_lox.tsv`, `results/tables/lox_detection_rates_by_sample.tsv`, `results/tables/external_gse223049_lox_validation.tsv`, `results/tables/external_gse223049_lox_validation_summary.tsv`, `results/figures/annotation_sanity/strict_marker_positive_lox_by_sample.tsv`, and `results/figures/annotation_sanity/strict_marker_positive_lox_age_summary.tsv`.

**Problem**

"Validation analyses" again overstates external descriptive comparison and annotation checks.

**Safer replacement**

> New sensitivity, annotation-sanity, and external-comparison analyses for this version include `scripts/fb_composition_adjusted_lox.py`, `scripts/lox_detection_rates_by_sample.py`, `scripts/annotation_uncertainty_check.py`, and `scripts/external_validation_gse223049_lox.py`; reports in `reports/fb_composition_adjustment.md`, `reports/lox_detection_rate_analysis.md`, `reports/annotation_uncertainty_check.md`, `reports/external_dataset_search.md`, `reports/external_validation_plan.md`, `reports/highest_impact_next_analysis.md`, and `reports/hypothesis_falsification.md`; and tables in `results/tables/fb_composition_adjusted_lox.tsv`, `results/tables/lox_detection_rates_by_sample.tsv`, `results/tables/external_gse223049_lox_validation.tsv`, `results/tables/external_gse223049_lox_validation_summary.tsv`, `results/figures/annotation_sanity/strict_marker_positive_lox_by_sample.tsv`, and `results/figures/annotation_sanity/strict_marker_positive_lox_age_summary.tsv`.

### 46. Figure legend says "FDR-significant" without small-n caveat

**Quoted sentence**

> Labeled points indicate key FDR-significant LOX-family results.

**Problem**

This is statistically true if the table reports padj, but in figure-legend context it may imply robust inference despite small biological sample counts and unmodeled preparation effects.

**Safer replacement**

> Labeled points indicate LOX-family results meeting the table's FDR threshold in this small-sample pseudobulk analysis.

### 47. Figure legend says "Subtype-specific" without annotation caveat

**Quoted sentence**

> (A) Subtype-specific pseudobulk log2 fold changes for LOX-family genes in selected fibroblast and mTEC1 subtypes.

**Problem**

"Subtype-specific" can sound like independently validated subtype biology. The subtypes are published annotations with remaining uncertainty.

**Safer replacement**

> (A) Pseudobulk log2 fold changes for LOX-family genes in selected annotated fibroblast and mTEC1 subtypes.

### 48. Supplementary figure legend uses robustness language

**Quoted sentence**

> These plots are robustness visualizations and do not replace pseudobulk DESeq2 inference.

**Problem**

"Robustness visualizations" is too reassuring for simple per-sample plots.

**Safer replacement**

> These plots are per-sample descriptive visualizations and do not replace pseudobulk DESeq2 inference.

## Highest-priority edits before public posting

1. Replace "remodeling" with "transcript changes" or "transcript differences" in the title, Abstract, Results, and Discussion unless explicitly referring to prior literature.
2. Replace "validation" with "external comparison" or "directional consistency" for GSE223049, except when explicitly saying it is not subtype-resolved validation.
3. Remove "possible mediator" from the Abstract and Discussion unless a direct caveat is placed in the same sentence that no functional role is tested.
4. Replace "strongly reduced", "strongest", "survived", "retained", and "supported" with sign-based or descriptive wording in low-n and sensitivity-analysis contexts.
5. Make every GSE223049 sentence specify "broad" fibroblast or epithelial directionality and avoid language implying capsFB, medFB, or mTEC1 replication.
