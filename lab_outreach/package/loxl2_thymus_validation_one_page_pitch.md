# LOXL2 in Mouse Thymic Epithelial Aging: Collaboration Pitch

## One-Sentence Summary

LOXL2 is a recurrent candidate transcript-level signal in mouse thymic epithelial/stromal aging analyses, but it remains unvalidated at the protein/spatial level and requires orthogonal experimental testing.

## Current Computational Evidence

- Primary signal: public mouse thymus single-cell reanalysis identified lower `Loxl2` transcript detection in aged annotated mTEC1 cells in GSE240016.
- Biological unit: n=2 young and n=2 aged samples.
- Directional checks: sample-level comparisons and depth-matched downsampling remain directionally young-higher/aged-lower.
- Current status: candidate-level transcript signal.

## Why This Remains Candidate-Level

- `Loxl2` detection is sparse, especially in aged mTEC1.
- Aged mTEC1 samples have lower sequencing depth/detected genes.
- Dropout/depth contribution remains partial and unresolved.
- Batch/sample confounding is not ruled out.
- No protein, spatial, functional, activity, or causal experiment has been performed.

## External Context

E-MTAB-8560 provides independent mTEC-focused transcript-level context. The clearest aged-lower `Loxl2` patterns are in mTEClo and cTEC, with weaker support in mTEChi. Overall classification remains mixed/inconclusive, and this is not exact GSE240016 mTEC1 validation.

## Human Exploratory Context

Exploratory public human TEC/epithelial datasets provide transcript-level context for LOX-family genes, but they do not establish human conservation or validate the mouse mTEC1 candidate signal.

## Critical Missing Experiment

Mouse thymus LOXL2 IHC/IF or RNA-ISH across young and aged biological animals, ideally with epithelial/medullary markers.

## Proposed Collaboration

We are seeking a wet-lab collaborator or core facility to advise on and potentially run a pilot mouse thymus LOXL2 protein/spatial or RNA-localization assay.

## What Result Would Strengthen the Project

Detectable LOXL2 protein or `Loxl2` RNA signal in epithelial/mTEC-like regions, with a reproducible young-vs-aged difference across independent animals, would provide supportive orthogonal mouse context.

## What Result Would Weaken or Falsify the Candidate

No detectable LOXL2 signal, no young-vs-aged difference in epithelial/mTEC-like regions, or nonspecific staining that cannot be resolved would weaken the protein/spatial interpretation and motivate RNA-ISH or sorted TEC qPCR/ddPCR instead.

## What This Does Not Claim

This project does not claim human conservation, validation of the mouse result, mechanism, causality, thymic function, LOX activity, ECM crosslinking, therapeutic relevance, rejuvenation, or protein-level validation.

## Links

- GitHub repository: https://github.com/G1F12/ThymusLOXScan
- GitHub release tag: `v5.5-external-mtec-context`
- Public manuscript source: `manuscript/LOX_thymus_aging_public_preprint_v5_5_external_context.md`
- Public manuscript PDF: `manuscript/pdf/LOX_thymus_aging_public_preprint_v5_5_external_context.pdf`
