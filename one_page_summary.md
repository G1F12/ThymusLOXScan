# One-page summary: LOX-family transcript changes in aging murine thymic stroma

## Project
Computational reanalysis of public murine thymic stromal scRNA-seq data to examine LOX-family transcript patterns during aging.

## Dataset
- Primary dataset: GSE240016.
- 22,932 CD45-negative thymic stromal cells.
- Young 2-month vs aged 18-month female C57BL/6 mice.
- Biological-sample-aware pseudobulk was used as the main inferential framework.

## Main finding
LOX-family transcript differences in aging thymic stroma appear subtype-dependent rather than uniformly fibroblast-wide.

## Key internal observations
- capsFB: aged-lower Lox.
- medFB: aged-higher Loxl1 and aged-lower Loxl2.
- mTEC1: aged-lower Loxl2 candidate signal, but n=2 vs n=2 and sparse detection.

## External context
- GSE223049: broad sorted mouse thymic fibroblast and epithelial directional context.
- E-MTAB-8560: TEC/mTEC-like Loxl2 age-series context.
- E-MTAB-8560 adds independent mTEC-focused transcript-level context, with an aged-lower Loxl2 pattern clearest in mTEClo/cTEC but overall mixed/inconclusive and not an exact GSE240016 mTEC1 validation.
- GSE231906: metadata-only human candidate dataset; no expression-level human conservation claim.

## Technical caveats
- mTEC1 Loxl2 dropout/depth risk is Partial / unresolved.
- Aged mTEC1 samples have lower sequencing depth.
- The mTEC1 Loxl2 direction persists after pairwise sample checks and depth-matched downsampling.
- Public protein/spatial resources provide weak indirect support only.
- No direct mouse mTEC1 LOXL2 protein localization or age-aware thymic epithelial protein evidence was found.

## What this does not claim
- No causal claim.
- No protein-level validation.
- No LOX enzymatic activity claim.
- No ECM-crosslinking claim.
- No thymic functional or rejuvenation claim.
- No therapeutic claim.
- No human conservation claim.
- No exact subtype-resolved external validation.

## Most useful next experiment
LOXL2 IHC or IF on young versus aged mouse thymus sections, ideally co-stained with epithelial/medullary markers such as EPCAM, KRT5, KRT8, AIRE, CCL21A, or related TEC markers.

## Public repository
https://github.com/G1F12/ThymusLOXScan

## Current manuscript
manuscript/LOX_thymus_aging_public_preprint_v5_5_external_context.md

## Current release
v5.5-external-mtec-context
