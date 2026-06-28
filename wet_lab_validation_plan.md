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

## Candidate LOXL2 antibodies to evaluate

These candidates are starting points for discussion with a histology or imaging lab. None should be treated as guaranteed for the exact mouse thymus aging question, and mouse tissue compatibility needs to be checked before interpreting a young/aged comparison.

| Antibody / catalogue | Vendor/resource | Host/clonality | Listed applications | Listed species reactivity | IHC/IF relevance | Caveats | Suggested use in pilot |
|---|---|---|---|---|---|---|---|
| CAB025848 / MAB2639 | Human Protein Atlas / R&D Systems | Mouse monoclonal | HPA lists tissue IHC with IHC reliability approved; R&D product page lists IHC among applications | HPA context is human protein atlas use; R&D page lists human as the main species reactivity and mouse under cited use rather than as the primary listed species | Most relevant public IHC-associated LOXL2 starting point from HPA | Mouse-on-mouse staining may be difficult; mouse thymus and mTEC specificity are not established by these public pages | Discuss first as an FFPE IHC feasibility candidate, with positive-control tissue and species-specific controls |
| HPA036257 | Human Protein Atlas / Atlas Antibodies | Rabbit polyclonal | HPA lists ICC/IF enhanced and protein array supported; HPA IHC is not available | Human antigen context on HPA | Potential IF-oriented reagent to discuss if tissue IF optimization is planned | HPA does not provide tissue IHC support for this antibody; polyclonal specificity and mouse thymus performance need local checks | Consider only as an exploratory IF/orthogonal candidate after checking mouse cross-reactivity and control tissue |
| HPA056542 | Human Protein Atlas / Atlas Antibodies | Rabbit polyclonal | HPA lists ICC/IF enhanced, WB uncertain, and protein array supported; HPA IHC is not available | Human antigen context on HPA | Potential IF-oriented reagent to discuss if a second rabbit polyclonal is useful | HPA does not provide tissue IHC support for this antibody; WB uncertainty and mouse thymus performance are caveats | Consider as a secondary/orthogonal feasibility reagent rather than the primary FFPE IHC choice |
| ab96233 | Abcam | Rabbit polyclonal | Vendor page lists WB, IHC-P, ICC/IF, and ELISA-style use | Vendor page is human-focused on the accessed product information | Relevant to human IHC-P/IF examples, not direct mouse thymus support | Not a direct mouse thymus solution unless mouse reactivity and tissue performance are separately supported | Discuss only if the lab has prior mouse experience with this reagent or can run appropriate controls |

Source pages checked for this planning table: HPA LOXL2 antibody summary (`https://www.proteinatlas.org/ENSG00000134013-LOXL2/summary/antibody`), R&D Systems MAB2639 product page (`https://www.rndsystems.com/products/human-lysyl-oxidase-homolog-2-loxl2-antibody-262418_mab2639`), and Abcam ab96233 product page (`https://www.abcam.com/en-us/products/primary-antibodies/loxl2-antibody-ab96233`).

## Suggested pilot staining conditions to discuss with the lab

- FFPE versus frozen sections should be chosen based on the antibody page, lab experience, and control tissue availability.
- Start from vendor or HPA antigen retrieval and dilution guidance rather than inventing new conditions.
- Include no-primary and secondary-only controls.
- Include a positive-control tissue recommended by the vendor or HPA when feasible.
- Use identical imaging settings for young and aged sections.
- Prefer co-staining, or adjacent sections, with epithelial and medullary markers such as EPCAM, KRT8, KRT5, AIRE, CCL21A, or H2-Aa depending on antibody compatibility.
- Consider an orthogonal antibody if the initial signal is ambiguous.

## Expected staining interpretation

- Young-higher epithelial or medullary LOXL2 staining would support the transcript candidate.
- Similar young and aged staining would weaken the protein-level version of the candidate.
- Failed or nonspecific staining is inconclusive, not negative.
- Antibody specificity controls determine whether the result is interpretable.

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
