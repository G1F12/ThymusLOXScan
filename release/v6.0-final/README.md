# ThymusLOXScan v6.0 final

This public package contains a frozen computational preprint on subtype-dependent LOX-family transcript changes in aging murine thymic stroma. Its primary dataset is GSE240016. The principal candidate observation is aged-lower `Loxl2` transcript detection/expression in annotated mTEC1, based on 2+2 sample labels; it is hypothesis-generating and not a protein, functional, causal, or therapeutic result.

External datasets provide broad, related-subtype, mixed, or detection-only context. GSE231906 is `DETECTION_CONTEXT_ONLY`: it supports within-unit detection summaries among uniquely joined barcodes, not donor-level aging inference or cross-unit mean-expression comparisons.

## Contents

- `manuscript/`: final Markdown, PDF, HTML, and editable DOCX.
- `tables/`: the three compact authoritative tables cited in the manuscript and the GSE231906 strict join audit.
- `governance/`: compact datasets, claims, output-authority, and GSE231906 authority records.
- `reproducibility/`: environment specifications, requirements, input acquisition guidance, and frozen-output checks.

## Inspection and reuse

Read the PDF first, then compare its tables with the compact source tables and governance records. Obtain public inputs from GEO or BioStudies using the accessions in `governance/datasets_v6.tsv`. The release provides environment specifications, scripts, input provenance and frozen output checks. A fully independent clean-room reproduction was not completed before this release.

## Citation and contact

See `CITATION.cff`. The repository is https://github.com/G1F12/ThymusLOXScan. License: MIT.
