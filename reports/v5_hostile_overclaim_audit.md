# v5 hostile overclaim audit

Input manuscript: `manuscript/LOX_thymus_aging_public_preprint_v5_external_validation.md`

Safety standard compared against: `reports/v4_2_final_safety_check.md` and `reports/v5_safety_check.md`.

## Summary

The audit found no direct claims of causality, protein-level validation, LOX enzymatic activity changes, ECM crosslinking changes, thymic functional change, therapeutic relevance, human conservation, or exact subtype-resolved external validation as established results. The main risk was repeated use of validation/support/replication language around descriptive external comparisons.

## Flagged Sentences and Safer Replacements

### 1

Quoted sentence: "The goal was not to prove a mechanism, but to identify candidate stromal expression patterns that may motivate future experimental validation and independent subtype-resolved replication."

Risk: Uses validation/replication language that could imply stronger confirmatory evidence than the current computational reanalysis provides.

Safer replacement: "The goal was not to prove a mechanism, but to identify candidate stromal expression patterns that may motivate future experimental testing and independent subtype-resolved reanalysis."

### 2

Quoted sentence: "Because n=2 per age group limits inference, this result should be treated as a candidate mTEC1-associated transcript pattern requiring further validation rather than as functional evidence."

Risk: Validation language is stronger than warranted for a two-versus-two internal comparison.

Safer replacement: "Because n=2 per age group limits inference, this result should be treated as a candidate mTEC1-associated transcript pattern requiring further independent testing rather than as functional evidence."

### 3

Quoted sentence: "Thus, GSE223049 sorted bulk RNA-seq is directionally consistent with aged-lower broad fibroblast `Lox`/`Loxl2` and broad epithelial `Loxl2`, especially `Loxl2`, but it does not test exact fibroblast-subtype or mTEC1-specific claims and is not subtype-resolved validation."

Risk: Mostly cautious, but ending with validation language invites ambiguity.

Safer replacement: "Thus, GSE223049 sorted bulk RNA-seq is directionally similar to aged-lower broad fibroblast `Lox`/`Loxl2` and broad epithelial `Loxl2`, especially `Loxl2`, but it does not test exact fibroblast-subtype or mTEC1-specific claims and is not subtype-resolved evidence."

### 4

Quoted sentence: "We next assembled a cross-dataset validation matrix across the internal GSE240016 pseudobulk results, GSE223049 sorted bulk RNA-seq, E-MTAB-8560 TEC Smart-seq2 age-series summaries, and GSE231906 human thymus metadata."

Risk: The matrix is descriptive and includes metadata-only rows, so calling it a validation matrix overstates the evidentiary level.

Safer replacement: "We next assembled a cross-dataset comparison matrix across the internal GSE240016 pseudobulk results, GSE223049 sorted bulk RNA-seq, E-MTAB-8560 TEC Smart-seq2 age-series summaries, and GSE231906 human thymus metadata."

### 5

Quoted sentence: "The cross-dataset matrix supported three broad or appropriately resolved directional conclusions."

Risk: Supported can be read as confirmatory; the matrix provides directional context only.

Safer replacement: "The cross-dataset matrix provided descriptive directional context for three broad or appropriately resolved transcript-level comparisons."

### 6

Quoted sentence: "Third, the mTEC1 `Loxl2` direction had approximate mTEC-like support from E-MTAB-8560 mTEC-like, mTEClo, and mTEChi summaries, but this was not exact mTEC1 replication."

Risk: Support/replication language could imply exact mTEC1 validation despite the caveat.

Safer replacement: "Third, the mTEC1 `Loxl2` direction had approximate mTEC-like directional context from E-MTAB-8560 mTEC-like, mTEClo, and mTEChi summaries, but this was not exact mTEC1 recapitulation."

### 7

Quoted sentence: "GSE231906 contributed metadata-defined candidate human fibroblast-like and epithelial groups but no expression-based validation in this repository."

Risk: Validation language is inappropriate for metadata-only output.

Safer replacement: "GSE231906 contributed metadata-defined candidate human fibroblast-like and epithelial groups but no expression-based evidence in this repository."

### 8

Quoted sentence: "Because expression was not parsed, GSE231906 was not used to claim human conservation or expression-level validation."

Risk: Good caveat, but expression-level validation language is best avoided entirely.

Safer replacement: "Because expression was not parsed, GSE231906 was not used to claim human conservation or expression-level evidence."

### 9

Quoted sentence: "This computational reanalysis suggests a narrow central conclusion: age-associated LOX-family transcript differences in thymic stroma are subtype-dependent, with broad external directional consistency for selected fibroblast and epithelial directions."

Risk: Directional consistency is acceptable but can read too strong in the Discussion opening.

Safer replacement: "This computational reanalysis suggests a narrow central conclusion: age-associated LOX-family transcript differences in thymic stroma are subtype-dependent, with broad external directional context for selected fibroblast and epithelial directions."

### 10

Quoted sentence: "Fourth, GSE231906 was handled as metadata-only in this repository because the expression archive is large and requires guarded donor/source matching; it therefore does not provide expression-level human validation here."

Risk: Metadata-only data cannot provide validation; safer to say evidence.

Safer replacement: "Fourth, GSE231906 was handled as metadata-only in this repository because the expression archive is large and requires guarded donor/source matching; it therefore does not provide expression-level human evidence here."

### 11

Quoted sentence: "Fifth, the mTEC1 `Loxl2` comparison in the original single-cell dataset includes only two biological samples per age group, and medFB `Loxl1` increase lacks subtype-resolved external replication."

Risk: Replication implies a precise external test; none exists for medFB.

Safer replacement: "Fifth, the mTEC1 `Loxl2` comparison in the original single-cell dataset includes only two biological samples per age group, and medFB `Loxl1` increase lacks subtype-resolved external reanalysis."

### 12

Quoted sentence: "External datasets also use different age ranges, species, technologies, sorting strategies, annotation systems, and compartment definitions, so cross-dataset agreement should be read as directional context rather than validation of identical cell states."

Risk: The sentence is already cautious, but removes validation framing from the final phrase.

Safer replacement: "External datasets also use different age ranges, species, technologies, sorting strategies, annotation systems, and compartment definitions, so cross-dataset agreement should be read as directional context rather than evidence for identical cell states."

## Keyword Scan Terms

`validated`, `validation`, `support`, `supported`, `replication`, `replicated`, `conserved`, `conservation`, `mechanism`, `functional`, `therapeutic`, `rejuvenation`, `crosslinking`, `causality`, `causal`, `strongest`, `robust`, `proves`.

Remaining occurrences in v5.1 are limited to cautious negations, background citations, method/report filenames, or statements that required future assays rather than present evidence.