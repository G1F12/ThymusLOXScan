# Mouse Thymus LOXL2 Wet-Lab Validation Proposal

## Background

ThymusLOXScan is an independent computational reanalysis of public thymus aging datasets focused on LOX-family transcript patterns in stromal and epithelial compartments.

## Current Computational Finding

The primary mouse observation is a candidate-level `Loxl2` transcript signal in annotated GSE240016 mTEC1 cells, with lower detection in aged samples than young samples.

## Key Limitations

- n=2 young and n=2 aged biological samples.
- Sparse aged mTEC1 detection.
- Aged mTEC1 samples have lower depth/detected genes.
- Dropout/depth contribution remains partial and unresolved.
- Batch/sample confounding is not ruled out.
- No protein, spatial, activity, function, or causal assay has been performed.

## Experimental Question

Is LOXL2 protein or `Loxl2` RNA detectable in mouse thymic epithelial/medullary regions, and does the signal differ between young and aged animals?

## Preferred Experiment: LOXL2 IHC/IF

Perform LOXL2 IHC or IF on young and aged mouse thymus sections with epithelial and medullary co-markers where feasible.

## Alternative Experiment: Loxl2 RNA-ISH / RNAscope

Use RNA-ISH/RNAscope to localize `Loxl2` transcript in thymic epithelial/medullary regions if antibody performance is uncertain or nonspecific.

## Alternative Experiment: Sorted TEC/mTEC qPCR/ddPCR

Sort TEC/mTEC-enriched populations and test `Loxl2` transcript abundance by qPCR or ddPCR as an independent transcript-level check.

## Sample Design

- Preferred: n >= 3 young + n >= 3 aged mouse thymus sections from independent biological animals.
- Pilot: n=1 young + n=1 aged only as technical feasibility.
- Technical sections alone should not be treated as biological replication.

## Marker Panel

- LOXL2
- EPCAM
- KRT5
- KRT8
- AIRE
- CCL21A if feasible
- DAPI

## Controls

- No-primary control.
- Secondary-only control.
- Positive-control tissue if antibody datasheet supports it.
- Antibody specificity check or orthogonal antibody if feasible.
- Identical imaging settings across groups.
- Blinded quantification if possible.

## Imaging/Quantification Plan

- Acquire matched cortical/medullary regions per animal.
- Segment or annotate epithelial-rich and medullary marker-positive regions.
- Quantify LOXL2 signal intensity or positive area by animal-level summaries.
- Report animal-level summaries, not cell-level tests as biological evidence.

## Expected Outcomes

The experiment can show whether LOXL2 protein or `Loxl2` RNA is detectable in relevant regions and whether the candidate direction is supported at the animal level.

## Interpretation Table

| Outcome | Interpretation |
|---|---|
| Young higher than aged in epithelial/mTEC-like regions | Supportive orthogonal mouse protein/spatial context, not causality. |
| No difference | Weakens protein-level candidate interpretation. |
| Nonspecific staining | Inconclusive; move to RNA-ISH/qPCR. |
| RNA-ISH supports transcript localization | Supports transcript-level spatial context, not protein/function. |
| qPCR supports sorted TEC decrease | Independent transcript-level support, not protein/function. |

## What the Experiment Cannot Prove

It cannot prove mechanism, causality, thymic function, LOX activity, ECM crosslinking, therapeutic relevance, rejuvenation, or human conservation.

## Estimated Timeline

- Week 1: assay feasibility and antibody/probe review.
- Weeks 2-4: pilot staining or RNA-ISH.
- Weeks 5-6: animal-level quantification and interpretation.

## What Collaborator Needs To Provide

- Mouse thymus tissue or access to sections.
- Histology/IHC/IF/RNA-ISH workflow advice or execution.
- Marker-panel feasibility input.
- Imaging and quantification guidance.

## What the Computational Side Can Provide

- Candidate-gene rationale.
- Public dataset summaries and caveat tables.
- Suggested analysis plan for animal-level quantification.
- Figure-ready context summaries after experimental data are generated.
