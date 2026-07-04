# GSE240016 Metadata Source Notes

## Official URLs Inspected

- GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240016
- GEO SOFT metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240016/soft/GSE240016_family.soft.gz
- GEO MINiML metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240016/miniml/GSE240016_family.xml.tgz
- NCBI SRA E-utilities RunInfo queried from SRX accessions linked in GEO sample records.
- ENA browser API queried for SRR run accessions when available.
- NCBI BioSample E-utilities queried for BioSample accessions from SRA RunInfo.

## Field Sources

- GEO: GSM accession, sample title, source name, sample characteristics, library strategy/source/selection, protocol text, and SRX relation.
- SRA RunInfo: SRR run accessions, SRX experiment, BioSample, BioProject, SRA study, library layout, platform/model, center, load date, release date, spots, and bases.
- BioSample: BioSample attributes where available.
- ENA: run-level cross-check fields where the API returned records.

## Not Found

- Mouse ID.
- Litter.
- Exact processing date.
- Library prep date.
- Sequencing lane.
- Flowcell.
- Explicit batch field.

## Manual Mapping Assumptions

- Local `mo02_CD45neg1_d0` maps to GEO steady-state 02-month CD45-negative thymic stromal replicate 1.
- Local `mo02_CD45neg2_d0` maps to GEO steady-state 02-month CD45-negative thymic stromal replicate 2.
- Local `mo18_CD45neg1_d0` maps to GEO steady-state 18-month CD45-negative thymic stromal replicate 1.
- Local `mo18_CD45neg2_d0` maps to GEO steady-state 18-month CD45-negative thymic stromal replicate 2.

These assumptions are high confidence because the age, enrichment, steady-state/untreated timepoint, and replicate tokens match directly.
