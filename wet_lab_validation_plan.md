# Wet-lab validation plan: LOXL2 and LOX-family candidates in aging thymic stroma

## Goal
Evaluate whether transcript-level LOX-family candidates observed in public scRNA-seq data correspond to reproducible RNA or protein differences in young versus aged thymic tissue.

## Priority 1: LOXL2 protein localization in thymic epithelium
Question:
Is LOXL2 protein detectable in young mouse thymic medullary/epithelial regions, and is it reduced or redistributed in aged thymus?

Recommended assay:
- Immunohistochemistry or immunofluorescence for LOXL2.
- Young vs aged mouse thymus sections.
- Co-staining with epithelial/medullary markers if possible.

Suggested markers:
- EPCAM for epithelial cells.
- KRT8 and KRT5 for TEC compartments.
- AIRE / CCL21A / H2-Aa depending on available antibodies and lab focus.
- DAPI/nuclear counterstain.

Minimum design:
- At least 3 young and 3 aged mice if feasible.
- If only pilot testing is possible: 1 young and 1 aged mouse can test feasibility but cannot support strong inference.
- Multiple sections per animal.
- Matched staining and imaging conditions.

Controls:
- no-primary control;
- secondary-only control;
- positive-control tissue recommended by antibody vendor;
- antibody specificity check using an orthogonal antibody if feasible;
- same exposure/imaging settings across young and aged groups.

Readouts:
- LOXL2 signal intensity in epithelial/medullary regions;
- fraction of epithelial/medullary cells or regions with detectable LOXL2;
- co-localization with epithelial markers;
- comparison of young vs aged signal under identical imaging conditions.

Interpretation:
- Young-higher epithelial LOXL2 protein would support the mTEC/epithelial transcript candidate.
- Similar LOXL2 protein in young and aged would suggest transcript-only, dropout-sensitive, or post-transcriptionally buffered signal.
- No reliable LOXL2 staining would make the antibody/protein-level question unresolved rather than negative, unless controls are strong.

## Priority 2: RNA-level spatial validation
Question:
Does Loxl2 RNA localize to medullary/epithelial thymic regions and decrease with age?

Possible assays:
- RNAscope or RNA in situ hybridization for Loxl2.
- Spatial transcriptomics if available.
- Co-detection with epithelial/medullary markers.

Why useful:
This tests whether the single-cell RNA signal has tissue-spatial support while avoiding some antibody specificity limitations.

## Priority 3: Sorted-cell qPCR or bulk validation
Question:
Do sorted TEC or mTEC-enriched populations show aged-lower Loxl2?

Possible assays:
- FACS-sort TEC or mTEC-enriched populations.
- qPCR for Loxl2 and selected LOX-family genes.
- Include housekeeping genes and biological replicates.

Limitations:
Sorted-cell qPCR supports compartment-level RNA direction but does not prove protein localization or function.

## Priority 4: Fibroblast candidates
Questions:
- Is capsFB-associated Lox aged-lower?
- Is medFB-associated Loxl1 aged-higher and Loxl2 aged-lower?

Possible assays:
- Spatial RNA/protein staining with fibroblast markers.
- Sorted fibroblast subtype qPCR if markers and sorting strategy are available.
- Multiplex RNA in situ hybridization if feasible.

## What would count as strong support
- Reproducible young-vs-aged difference across biological animals.
- Same direction as computational result.
- Signal located in the relevant epithelial or fibroblast compartment.
- Appropriate technical controls.
- Ideally, agreement across RNA and protein-level assays.

## What would refute or weaken the candidate
- No difference across biological replicates despite reliable detection.
- Opposite direction in well-controlled tissue assays.
- LOXL2 protein signal not localized to the expected epithelial/medullary regions.
- Antibody staining fails specificity controls.
- RNA signal disappears in spatial/RNAscope validation.

## What this plan does not test
- Causality.
- LOX enzymatic activity.
- ECM crosslinking.
- Thymic functional output.
- Therapeutic potential.

## Recommended first experiment
A small feasibility IHC/IF pilot:
- young and aged mouse thymus sections;
- LOXL2 antibody;
- epithelial marker co-staining;
- identical imaging settings;
- qualitative plus simple quantitative comparison.

## Link to computational work
Repository:
https://github.com/G1F12/ThymusLOXScan

Current manuscript:
manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md
