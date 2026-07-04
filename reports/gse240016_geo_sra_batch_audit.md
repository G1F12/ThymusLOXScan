# GSE240016 GEO/SRA Batch Metadata Audit

## Scope

This audit joins local AnnData sample labels to official GEO/SRA metadata for GSE240016. It evaluates whether the four steady-state CD45-negative mTEC1-contributing samples are visibly confounded with available accession, run, library, center, instrument, or timing metadata.

## Mapping Result

| local sample | stage | mapped GSM | SRR runs | BioSample | confidence | basis |
|---|---|---|---|---|---|---|
| mo02_CD45neg1_d0 | 02mo | GSM7679869 | SRR25509581;SRR25509582;SRR25509583;SRR25509584 | SAMN36828570 | high | 02-mo CD45- thymic stroma steady state rep1 |
| mo02_CD45neg2_d0 | 02mo | GSM7679870 | SRR25509585;SRR25509586;SRR25509587;SRR25509588 | SAMN36828569 | high | 02-mo CD45- thymic stroma steady state rep2 |
| mo18_CD45neg1_d0 | 18mo | GSM7679885 | SRR25509527;SRR25509528 | SAMN36828554 | high | 18-mo CD45- thymic stroma steady state rep1 |
| mo18_CD45neg2_d0 | 18mo | GSM7679886 | SRR25509525;SRR25509526 | SAMN36828553 | high | 18-mo CD45- thymic stroma steady state rep2 |

All four local mTEC1 sample labels were mapped to official GSM/SRX/SRR/BioSample records with high confidence using age, CD45-negative stromal sample title, steady-state/untreated timepoint, and replicate tokens.

## Local mTEC1 QC

| local sample | stage | total cells | mTEC1 cells | median total counts mTEC1 | median detected genes mTEC1 | mean total counts mTEC1 | mean detected genes mTEC1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| mo02_CD45neg1_d0 | 02mo | 3381 | 300 | 8038.5 | 2819.5 | 9107.5 | 2847.5 |
| mo02_CD45neg2_d0 | 02mo | 3015 | 279 | 9032.0 | 3006.0 | 9852.2 | 3061.1 |
| mo18_CD45neg1_d0 | 18mo | 1997 | 289 | 5641.0 | 2317.0 | 6761.5 | 2462.6 |
| mo18_CD45neg2_d0 | 18mo | 4289 | 764 | 5457.0 | 2175.5 | 6343.7 | 2322.5 |

Aged mTEC1 cells remain lower for local median total counts: yes.
Aged mTEC1 cells remain lower for local median detected genes: yes.

## Confound Matrix Summary

| technical field | age coincides? | assessment |
|---|---|---|
| GSM accession | yes_or_sample_specific | values differ between age groups or are sample-specific |
| SRR run accessions | yes_or_sample_specific | values differ between age groups or are sample-specific |
| BioSample | yes_or_sample_specific | values differ between age groups or are sample-specific |
| sequencing run date | unavailable | field not present in GEO/SRA/BioSample metadata |
| SRA load timestamp | yes_or_sample_specific | values differ between age groups or are sample-specific |
| SRA load calendar date | no | available values do not separate young and aged samples |
| SRA release calendar date | no | available values do not separate young and aged samples |
| library prep date | unavailable | field not present in GEO/SRA/BioSample metadata |
| library strategy | no | available values do not separate young and aged samples |
| library source | no | available values do not separate young and aged samples |
| library selection | no | available values do not separate young and aged samples |
| instrument | no | available values do not separate young and aged samples |
| center | no | available values do not separate young and aged samples |
| lane/flowcell | unavailable | field not present in GEO/SRA/BioSample metadata |
| sort/enrichment condition | no | available values do not separate young and aged samples |
| day/timepoint | no | available values do not separate young and aged samples |

## Explicit Answers

1. Each of the four mTEC1 local sample labels can be mapped to a GSM/SRX/SRR/BioSample record: yes, with high confidence.
2. The young and aged samples have distinct GSM/SRX/BioSample/SRR accessions. SRA lists multiple runs per sample. The available metadata do not provide lane or flowcell fields, so same-lane or different-lane status cannot be determined.
3. Age coincides with GSM, SRX, BioSample, and SRR accession blocks because each biological sample has its own accessions and the accession ranges separate by age. Age does not obviously coincide with instrument, center, library strategy/source/selection, enrichment condition, SRA load calendar date, or SRA release calendar date. SRA load timestamps differ by sample at minute-scale resolution, but these are administrative load timestamps, not sequencing run dates. Sequencing run date, library prep date, processing date, lane, and flowcell were not found.
4. CD45-negative stromal enrichment status is consistent across the four samples: yes (CD45- thymic stroma cells).
5. Aged samples are systematically lower for local mTEC1 depth/detected-gene medians: yes; this remains a local QC concern.
6. Official metadata reduces concern about a simple center, instrument, library-type, or public-load-date confound, but it does not resolve sample-level, processing-date, lane/flowcell, or library-prep confounding.
7. Unresolved fields include mouse ID, litter, exact processing date, library prep date, sequencing lane, flowcell, and whether all four steady-state CD45-negative libraries were processed in the same technical batch.

Classification: batch concern remains.

This classification is cautious. It does not state that batch confounding is ruled out.
