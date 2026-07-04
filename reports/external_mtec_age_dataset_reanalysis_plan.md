# External mTEC Age Dataset Reanalysis Plan

## P1 Candidate: E-MTAB-8560

Accession: E-MTAB-8560.

Files needed:
- MouseThymusAgeing SMART-seq2 processed SingleCellExperiment objects.
- MouseThymusAgeing SMARTseq metadata.
- ArrayExpress/BioStudies SDRF and IDF metadata for run, plate, sort, age, and individual fields if needed.

File size/type:
- Use processed count matrices and metadata through Bioconductor.
- Do not download FASTQ, BAM, SRA, or raw alignment files.
- SDRF/IDF metadata can be inspected or cached locally if needed, but raw sequencing files should not be downloaded.

Processed-data access:
- Bioconductor package: `MouseThymusAgeing`.
- Documentation: https://bioconductor.org/packages/release/data/experiment/html/MouseThymusAgeing.html
- Study metadata: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8560

Sample metadata join strategy:
- Use cell-level metadata fields for Age, SortType, SubType, SortDay, PlateID, CellID, and sizeFactor.
- Use ArrayExpress SDRF fields for individual/mouse, ENA run, library, and FACS marker metadata when needed.
- Define the biological unit as mouse/individual within age and sort group, not cell.
- Exclude cells lacking confident age, individual, or mTEC-relevant group labels.

Cell groups to test:
- Primary: mTEClo and mTEChi.
- Secondary: inferred Mature.mTEC, Post-Aire.mTEC, Tuft-like.mTEC, and combined mTEC-like groups if enough mice contribute per age.
- Broad TEC/cTEC groups should be reported separately and not used as mTEC-specific evidence.

LOX-family genes to test:
- Lox
- Loxl1
- Loxl2
- Loxl3
- Loxl4

Model:
- Pseudobulk per mouse x age x group.
- Primary model: expression ~ age.
- If estimable: expression ~ SortDay + age or expression ~ PlateID + age.
- If batch is partially confounded, report unadjusted and batch-aware sensitivity results without overinterpreting.
- Use exact/permutation tests where replicate counts and group definitions support it.
- Apply global FDR across the LOX-family gene by TEC-group matrix.

Expected output files:
- `results/tables/external_emtab8560_mtec_lox_pseudobulk.tsv`
- `results/tables/external_emtab8560_mtec_lox_models.tsv`
- `results/tables/external_emtab8560_mtec_lox_batch_audit.tsv`
- `reports/external_emtab8560_mtec_lox_reanalysis.md`

Claim wording if result supports Loxl2:
- "E-MTAB-8560 provides an independent mTEC-focused transcript-level test with directional support for aged-lower Loxl2 in [specific group]."
- "This is not exact GSE240016 mTEC1 replication unless labels match."

Claim wording if result does not support Loxl2:
- "E-MTAB-8560 does not support aged-lower Loxl2 under this annotation and age contrast."
- "The GSE240016 mTEC1 Loxl2 signal should remain candidate-level and dataset-specific unless additional evidence emerges."

## P2 Candidate: E-MTAB-8737

Accession: E-MTAB-8737.

Files needed:
- MouseThymusAgeing droplet processed SingleCellExperiment objects.
- MouseThymusAgeing droplet metadata.
- ArrayExpress/BioStudies SDRF and IDF metadata for source/run/replicate/HTO fields.

File size/type:
- Use processed count matrices and metadata only.
- Do not download FASTQ, BAM, SRA, or raw alignment files.

Processed-data access:
- Bioconductor package: `MouseThymusAgeing`.
- Documentation: https://bioconductor.org/packages/release/data/experiment/html/MouseThymusAgeing.html
- Study metadata: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8737

Sample metadata join strategy:
- Use droplet metadata to identify multiplexed samples, HTO tags, age/induction group, ZsG group, run, and replicate.
- Use cell-level labels for mTEC/cTEC/mTEC-like subsets after demultiplexing.
- Define biological unit conservatively from sample/replicate metadata; do not treat cells as independent biological replicates.

Cell groups to test:
- mTEC and mTEC-like droplet annotations.
- cTEC only as a separate context group.

LOX-family genes to test:
- Lox
- Loxl1
- Loxl2
- Loxl3
- Loxl4

Model:
- Pseudobulk per biological replicate/sample x age/induction group x cell group.
- Primary model: expression ~ age.
- If run or HTO batch is estimable: expression ~ run + age.
- Use exact/permutation tests where n allows.
- Apply global FDR across the LOX-family gene by group matrix.

Expected output files:
- `results/tables/external_emtab8737_mtec_lox_pseudobulk.tsv`
- `results/tables/external_emtab8737_mtec_lox_models.tsv`
- `results/tables/external_emtab8737_mtec_lox_batch_audit.tsv`
- `reports/external_emtab8737_mtec_lox_reanalysis.md`

Claim wording if result supports Loxl2:
- "E-MTAB-8737 provides secondary directional support for aged-lower Loxl2 in a mTEC-like lineage-tracing context."
- "This is not exact GSE240016 mTEC1 replication."

Claim wording if result does not support Loxl2:
- "E-MTAB-8737 does not support aged-lower Loxl2 under this lineage-tracing annotation."
- "The result should be interpreted as a sensitivity check, not as definitive evidence against the GSE240016 candidate signal."
