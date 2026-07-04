# GSE240016 GEO/SRA Batch Audit Update Summary

## Previous Local-Metadata Conclusion

The prior local-metadata conclusion was that the AnnData object exposed sample, stage, day, QC, and cell annotation fields, but not mouse ID, sequencing run, accession, library prep date, processing batch, sex, litter, lane, or flowcell. Residual sample or batch confounding therefore remained possible.

## New Official Metadata Findings

Official GEO/SRA metadata maps the four local steady-state CD45-negative samples to GSM, SRX, SRR, and BioSample accessions with high confidence. The available fields show the same library strategy/source/selection, instrument model, center, BioProject, enrichment condition, SRA load calendar date, and SRA release calendar date across the four samples.

However, accession blocks are sample-specific and separate by age, and the metadata still do not expose mouse ID, litter, processing date, library prep date, lane, flowcell, or explicit batch.

## Conclusion Change

Conclusion changed: no. The official metadata reduces concern about a simple center/instrument/library-type confound but does not resolve sample-level or unreported batch variables.

## Remaining Unresolved

- Mouse ID and litter.
- Processing date and library prep date.
- Lane and flowcell.
- Explicit batch labels.
- Whether local lower aged mTEC1 depth reflects biology, cell state, technical quality, or a mixture.

## Effect on mTEC1 Loxl2 Interpretation

The mTEC1 Loxl2 result should remain candidate-level and potentially confounded. Official metadata does not validate the biology and does not resolve the dropout/depth concern.

Batch classification: batch concern remains.
