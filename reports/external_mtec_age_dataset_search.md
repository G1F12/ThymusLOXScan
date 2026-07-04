# External mTEC Age Dataset Search

Search date: 2026-07-05.

## Search Sources Checked

- GEO and linked SRA metadata pages.
- ArrayExpress/BioStudies study pages and MAGE-TAB metadata.
- ENA-linked run metadata where exposed through BioStudies.
- PubMed/journal pages.
- Official lab/data GitHub repositories.
- Bioconductor MouseThymusAgeing package documentation.
- Existing local repository reports for previously used external datasets.

## Search Queries Used

- mouse thymus aging mTEC RNA-seq GEO
- mTEC aging bulk RNA-seq mouse GEO
- Aire positive mTEC aged mouse RNA-seq
- thymic epithelial cell aging single cell RNA-seq mouse mTEC GEO
- sorted mTEC old young mouse RNA-seq
- thymus aging mTEC ArrayExpress
- mouse thymic epithelial cells age RNA-seq mTEChi mTEClo
- site:ncbi.nlm.nih.gov/geo mouse mTEC aging RNA-seq
- site:ebi.ac.uk/biostudies thymus mTEC aging mouse RNA-seq
- site:github.com mTEC aging RNA-seq thymus mouse

## Top Candidate Datasets

### P1: E-MTAB-8560 / MouseThymusAgeing

E-MTAB-8560 is the strongest candidate found. It is an independent mouse thymic epithelial ageing dataset from Baran-Gale et al. with SMART-seq2 data from FACS-purified TEC populations across the first year of mouse life. The Bioconductor MouseThymusAgeing documentation states that cells were sorted from 5 mice at each age, with ages 1, 4, 16, 32, and 52 weeks, and sorted groups including mTEClo, mTEChi, cTEC, and Dsg3+ TEC. It also states that data are provided as processed count matrices with sample-level and feature-level metadata.

This meets the target better than GSE240016 for an independent mTEC-focused transcript-level test: mTEC-relevant sort labels are present, n=5 per age is available, mouse/individual fields are present in ArrayExpress metadata, processed counts are available, and run metadata exists. It is still not exact GSE240016 mTEC1 replication because the labels are mTEClo/mTEChi and inferred TEC subtypes rather than the Kousa et al. mTEC1 annotation.

### P2: E-MTAB-8737 / MouseThymusAgeing

E-MTAB-8737 is a secondary candidate from the same Baran-Gale study ecosystem. It profiles TECs from a beta5t-lineage tracing model and includes mTEC/cTEC structure with cell hashing. The package documentation describes 3 replicates per age for the droplet experiment. It is useful as a sensitivity dataset, but it is less direct than E-MTAB-8560 because the design is transgenic/lineage-tracing, replicate count is lower than the preferred n>=4/group, and source/run metadata require a careful join.

## Rejected or Lower-Priority Datasets

- GSE56928: sorted cTEC, mTEClo, and mTEChi microarray data across 1, 3, and 6 months. It is mTEC-specific but has only two replicates per age/subset, so it does not improve biological power over the current n=2 vs n=2 issue.
- GSE223049: has 5 young and 5 aged thymic epithelial bulk RNA-seq samples with processed counts, but the thymic epithelial population is broad and not mTEC-specific.
- GSE240016: current primary dataset; not independent and retains the n=2 vs n=2 mTEC1 limitation.
- GSE240017: same Kousa et al. study family, broad CD45-negative bulk, and no aged steady-state samples.
- GSE137699: mTEC single-cell dataset with processed matrices and metadata, but not an old-vs-young ageing dataset.
- GSE114651: sorted mTEC subset RNA-seq, but not an ageing comparison.
- GSE103967: thymic stromal reference with TEC/mTEC states, but not an ageing comparison.
- GSE236070: young mimetic TEC subset bulk RNA-seq, but no ageing comparison.

## Target Criteria Assessment

- Sorted mTEC / mTEC subtype: met by E-MTAB-8560; partially met by E-MTAB-8737.
- n>=4 young and aged: met by E-MTAB-8560 using selected younger age versus 52-week groups; not met by E-MTAB-8737.
- Processed counts available: met by E-MTAB-8560 and E-MTAB-8737 through MouseThymusAgeing.
- Biological sample IDs available: met by E-MTAB-8560 via individual/mouse fields; partial for E-MTAB-8737 because multiplexed sample/run structure needs careful handling.
- Batch/run metadata available: available in ArrayExpress/BioStudies SDRF for both, with E-MTAB-8560 also documenting block sorting intended to reduce age/batch confounding.

## Recommendation

Proceed with reanalysis of E-MTAB-8560 first. It is the only P1 dataset found and is the best available independent mTEC-focused transcript-level test for Loxl2 with better biological power than the current GSE240016 n=2 vs n=2 mTEC1 signal.

Use E-MTAB-8737 only as a secondary P2 sensitivity dataset. Treat GSE56928 and GSE223049 as directional context only.

## Concrete Next-Step Plan for E-MTAB-8560

- Download processed MouseThymusAgeing SMART-seq2 count matrices and cell metadata through Bioconductor, not raw FASTQ/BAM/SRA files.
- Parse `Age`, `SortType`, `SubType`, `SortDay`, `PlateID`, `CellID`, `sizeFactor`, and mouse/individual metadata.
- Aggregate counts per biological mouse within mTEClo, mTEChi, and mTEC-like inferred groups.
- Test LOX-family genes: Lox, Loxl1, Loxl2, Loxl3, Loxl4.
- Primary comparisons: 4 or 16 weeks versus 52 weeks, with 1 week handled cautiously because it is juvenile.
- Fit expression ~ age for each mTEC group and LOX-family gene; use expression ~ batch + age if SortDay/PlateID is estimable without perfect confounding.
- Apply global FDR across tested gene-by-group comparisons.
- Include exact/permutation tests where the per-mouse replicate structure supports it.
- Report whether E-MTAB-8560 provides independent mTEC-focused transcript-level directional support, no support, or mixed evidence for aged-lower Loxl2.

Do not describe any result as exact GSE240016 mTEC1 replication unless the label and biological unit can be matched more closely than currently indicated.
