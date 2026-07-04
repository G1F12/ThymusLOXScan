# E-MTAB-8560 Metadata Source Notes

## Official URLs Used

- BioStudies/ArrayExpress study page: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8560
- SDRF metadata: https://www.ebi.ac.uk/biostudies/files/E-MTAB-8560/E-MTAB-8560.sdrf.txt
- IDF metadata: https://www.ebi.ac.uk/biostudies/files/E-MTAB-8560/E-MTAB-8560.idf.txt
- ENA portal API queried from `Comment[ENA_RUN]` accessions.

## Metadata Files Fetched

- `data/external/metadata/emtab8560/emtab8560_sdrf.tsv`
- `data/external/metadata/emtab8560/emtab8560_idf.tsv`
- ENA run metadata was queried but no table was written.

## Key Fields Found

- Individual/mouse candidate field: `Characteristics[individual]` with values 1, 2, 3, 4, 5.
- Age field: `Characteristics[age]` / `Factor Value[age]` with values , 1, 16, 32, 4, 52.
- Cell-type/sort fields include `Characteristics[cell type]`, `Characteristics[FACS marker]`, and `Characteristics[inferred cell type]`; observed cell types include Dsg3+TEC, Dsg3+TEC.MiniBulk, EMPTY, cTEC, cTECMiniBulk, mTEChi, mTEChi.MiniBulk, mTEClo, mTEClo.MiniBulk.
- Run/batch-related fields include `Comment[ENA_RUN]`, `Assay Name`, source names, and library metadata.
- PlateID and SortDay are expected in MouseThymusAgeing colData rather than the raw SDRF field names.

## Missing or Unresolved

- The metadata fetch alone does not verify that exported R colData can map every analyzed cell to a biological individual.
- Batch-aware modeling requires successful R export and a downstream biological-unit audit.

## SDRF Fields

- `Source Name`
- `Comment[ENA_SAMPLE]`
- `Comment[BioSD_SAMPLE]`
- `Characteristics[organism]`
- `Characteristics[strain]`
- `Characteristics[genotype]`
- `Characteristics[individual]`
- `Characteristics[age]`
- `Unit[time unit]`
- `Term Source REF`
- `Term Accession Number`
- `Characteristics[developmental stage]`
- `Characteristics[sex]`
- `Characteristics[organism part]`
- `Characteristics[cell type]`
- `Characteristics[FACS marker]`
- `Characteristics[post analysis well quality]`
- `Characteristics[inferred cell type]`
- `Characteristics[single cell well quality]`
- `Material Type`
- `Protocol REF`
- `Performer`
- `Protocol REF.1`
- `Performer.1`
- `Protocol REF.2`
- `Performer.2`
- `Extract Name`
- `Comment[LIBRARY_LAYOUT]`
- `Comment[LIBRARY_SELECTION]`
- `Comment[LIBRARY_SOURCE]`
- `Comment[LIBRARY_STRAND]`
- `Comment[LIBRARY_STRATEGY]`
- `Comment[NOMINAL_LENGTH]`
- `Comment[NOMINAL_SDEV]`
- `Comment[ORIENTATION]`
- `Comment[end bias]`
- `Comment[input molecule]`
- `Comment[library construction]`
- `Comment[primer]`
- `Comment[single cell isolation]`
- `Comment[spike in]`
- `Comment[spike in dilution]`
- `Protocol REF.3`
- `Performer.3`
- `Assay Name`
- `Technology Type`
- `Comment[ENA_EXPERIMENT]`
- `Scan Name`
- `Comment[SUBMITTED_FILE_NAME]`
- `Comment[ENA_RUN]`
- `Comment[FASTQ_URI]`
- `Factor Value[single cell identifier]`
- `Factor Value[age]`
- `Unit[time unit].1`
- `Term Source REF.1`
- `Term Accession Number.1`
- `Factor Value[cell type]`
