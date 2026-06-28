# LOXL2 protein/spatial validation feasibility audit

## A. Executive summary

Public resources provide weak indirect support for the mTEC1/epithelial Loxl2 transcript finding, but they do not provide protein-level validation of mouse mTEC1 LOXL2. The best public support is spatial RNA-level context from a human thymus Visium collection and human thymic epithelial single-cell RNA data, where LOXL2 is detectable in medullary/epithelial contexts. This is not a substitute for IHC and does not validate mouse mTEC1 LOXL2 protein.

The strongest protein-level public context is general LOXL2 protein evidence in HPA and ProteomicsDB, plus HPA antibody feasibility. However, HPA thymus specifically reports no protein expression annotation, ProteomicsDB does not return a thymus LOXL2 row, and checked mouse thymus/TEC proteomics result files did not contain Loxl2.

Final classification: weak indirect support.

## B. What can be checked on PC

- HPA LOXL2 tissue and antibody pages can be checked for thymus protein annotation, antibody IDs, and assay reliability.
- ProteomicsDB can be queried for human LOXL2 protein expression across tissues.
- PRIDE/ProteomeXchange metadata and processed search-result files can be checked where available.
- Public H5AD files from CELLxGENE can be downloaded and queried for LOXL2 spatial RNA and thymic epithelial RNA context.

## C. What cannot be checked without wet-lab IHC

- Mouse mTEC1 LOXL2 protein localization.
- Whether LOXL2 protein is present in the same cells that carry the Loxl2 transcript signal.
- Whether antibody staining is specific in aged mouse thymus.
- Whether LOXL2 protein is extracellular, stromal-adjacent, epithelial, vascular, or background in the relevant mouse tissue sections.

## D. HPA / public IHC evidence

HPA lists LOXL2 as having protein-level evidence globally and lists antibodies CAB025848, HPA036257, and HPA056542. The global HPA row reports IH reliability as Approved and IF reliability as Enhanced. On the thymus tissue page, however, the expression summary reports protein expression as No data. The thymus page provides RNA context, including consensus RNA around 1.7 nTPM and low HPA single-cell RNA signal in human medullary thymic epithelial cells, but this is RNA-level context.

Interpretation: HPA is useful for antibody feasibility and general protein-level public context. It does not provide thymus-specific protein annotation that can be used as mouse mTEC1 protein support.

## E. Public proteomics evidence

ProteomicsDB returned human LOXL2 protein-expression rows for several tissues, including lymph node, spleen, and tonsil, but no thymus row was returned in the API result queried during this audit. This supports only broad protein-level public context.

PXD067907 is a mouse bulk thymus TMT dataset. Its processed protein table was downloaded and searched for Loxl2, LOXL2, and Q9R0B9; no Loxl2 hit was found. PXD042241 is a mouse cTEC/mTEC MHC-I immunopeptidomics dataset. Four downloaded mzIdentML result files were searched; no Loxl2 hit was found. PXD040722 is a mouse aging proteomic atlas but does not include thymus in the project description. PXD053154 is a promising 41-organ mouse aging/protein-restriction proteomics project, but this audit did not identify a processed thymus matrix suitable for confirming Loxl2 detection.

Interpretation: no true public protein dataset found here supports mouse mTEC1 LOXL2 protein localization or age-associated protein change.

## F. Spatial transcriptomics evidence

The CELLxGENE human thymus collection contains multiple human thymus Visium datasets and thymic epithelial single-cell subsets. In the merged pediatric human thymus Visium H5AD, LOXL2 was detected in 849 of 25,593 spots (3.32%); medulla-annotated spots showed 371 of 6,770 spots with signal (5.48%). In the human thymic epithelial subset H5AD, LOXL2 was detected in 597 of 25,726 TECs (2.32%), including mTEC and mTEC-mimetic categories.

Interpretation: this is useful spatial RNA-level context and epithelial RNA-level context in human thymus. It is not protein-level evidence and does not validate mouse mTEC1 LOXL2 protein.

## G. Antibody feasibility

HPA provides the most relevant public antibody starting points. CAB025848 / R&D Systems MAB2639 has HPA IHC support, but mouse-on-mouse IHC may be technically difficult and mouse thymus specificity would need controls. HPA036257 and HPA056542 have HPA ICC/IF support but no HPA IHC support. Abcam ab96233 is a human-reactive rabbit polyclonal with vendor-listed IHC-P and IF support, but the accessed product page lists human reactivity, so it is not a direct mouse thymus solution.

Feasible wet-lab follow-up would require positive-control tissue, no-primary/secondary controls, species-appropriate controls, and ideally an orthogonal antibody or peptide/blocking/knockdown-style specificity evidence.

## H. Final verdict

Computational/public evidence classification: weak indirect support.

Public resources provide protein-level public context and spatial RNA-level context, but they are not a substitute for IHC. No public dataset found in this audit provides direct mouse mTEC1 LOXL2 protein localization or an age-aware mouse thymic epithelial LOXL2 protein measurement. Wet-lab IHC or another direct protein assay is still required for protein-level support of the mouse mTEC1/epithelial Loxl2 transcript finding.

## Key sources

- HPA LOXL2 thymus page: https://www.proteinatlas.org/ENSG00000134013-LOXL2/tissue/thymus
- HPA LOXL2 antibody page: https://www.proteinatlas.org/ENSG00000134013-LOXL2/summary/antibody
- HPA downloadable proteinatlas.tsv: https://www.proteinatlas.org/download/proteinatlas.tsv.zip
- ProteomicsDB expression API: https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinexpression.xsodata/
- PRIDE PXD067907: https://www.ebi.ac.uk/pride/archive/projects/PXD067907
- PRIDE PXD042241: https://www.ebi.ac.uk/pride/archive/projects/PXD042241
- PRIDE PXD040722: https://www.ebi.ac.uk/pride/archive/projects/PXD040722
- PRIDE PXD053154: https://www.ebi.ac.uk/pride/archive/projects/PXD053154
- CELLxGENE human thymus collection: https://cellxgene.cziscience.com/collections/fc19ae6c-d7c1-4dce-b703-62c5d52061b4
- UniProt LOXL2 Q9Y4K0: https://www.uniprot.org/uniprotkb/Q9Y4K0/entry
- Abcam ab96233: https://www.abcam.com/en-us/products/primary-antibodies/loxl2-antibody-ab96233
