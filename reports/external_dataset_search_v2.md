# External Dataset Search v2

Search date: 2026-06-28.

Scope: second-pass search of GEO, ArrayExpress/BioStudies, CELLxGENE, and PubMed-linked supplementary resources for public thymus datasets that could test LOX-family thymic aging directions. No large expression archive was downloaded during this search.

## Scoring Rubric

| score | meaning |
|---:|---|
| 5 | Immediately usable processed expression, relevant age contrast, relevant thymic stromal/TEC groups. |
| 4 | Usable with moderate parsing; relevant age series/expression, but compartment or replicate caveats. |
| 3 | Possible but heavy; useful labels or age range, but large files or donor/source matching required. |
| 2 | Context only; useful thymus expression labels but no clean aging contrast. |
| 1 | Limited or indirect; not a practical validation dataset for current claims. |
| 0 | Not suitable for LOX-family aging validation. |

## Immediately Usable

| accession/resource | source | species | data type | age groups | n biological samples/donors if known | stromal/TEC/fibroblast coverage | expression matrix available | approximate file size | processed metadata available | feasibility | reason | claims it can test | claims it cannot test |
|---|---|---|---|---|---:|---|---|---|---|---:|---|---|---|
| GSE223049 | GEO; sample record confirms FACS purified thymic fibroblast RNA-seq with pooled animals per replicate: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8377364 | Mouse | Sorted bulk RNA-seq | 2 months vs 22-24 months | Thymic fibroblast: 4 young + 4 aged; thymic epithelial: 5 young + 5 aged in local parsed table | Broad `Thymic_fibroblasts` and broad `Thymic_epithelial`; no subtype resolution | Yes, processed counts are local | `GSE223049_RNA_seq_counts_23_cell_types.txt.gz`, 5.1 MB local copy | Yes | 5 | Best independent, small, directly parsed mouse comparison for broad fibroblast and epithelial directions. | Broad fibroblast `Lox` aged-lower; broad fibroblast `Loxl2` aged-lower; broad epithelial `Loxl2` aged-lower. | capsFB, medFB, exact mTEC1, human conservation, protein/activity/functional claims. |

## Usable With Moderate Work

| accession/resource | source | species | data type | age groups | n biological samples/donors if known | stromal/TEC/fibroblast coverage | expression matrix available | approximate file size | processed metadata available | feasibility | reason | claims it can test | claims it cannot test |
|---|---|---|---|---|---:|---|---|---|---|---:|---|---|---|
| E-MTAB-8560 / MouseThymusAgeing SMART-seq2 | ArrayExpress/BioStudies and Bioconductor MouseThymusAgeing vignette: https://bioconductor.org/packages/release/data/experiment/vignettes/MouseThymusAgeing/inst/doc/MouseThymusAgeing.html | Mouse | Smart-seq2 single-cell RNA-seq | 1, 4, 16, 32, 52 weeks | Local pipeline summarizes 5 acquisition/sort-day files; exact donor independence remains cautious | TEC-enriched: cTEC, mTEClo, mTEChi, mTEC-like groups; no fibroblasts | Yes, local processed RDS files | SDRF 18.2 MB; needed local processed RDS files about 38 MB total | Yes | 4 | Strong feasible TEC age-series context, but first-year age range and TEC labels are not exact mTEC1 replication. | Broad TEC/epithelial `Loxl2` aged-lower directional context; mTEC-like `Loxl2` approximate context. | Fibroblast claims; exact mTEC1 validation; human conservation; protein/activity/functional claims. |

## Possible But Heavy

| accession/resource | source | species | data type | age groups | n biological samples/donors if known | stromal/TEC/fibroblast coverage | expression matrix available | approximate file size | processed metadata available | feasibility | reason | claims it can test | claims it cannot test |
|---|---|---|---|---|---:|---|---|---|---|---:|---|---|---|
| GSE231906 | GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231906 | Human | 10x scRNA-seq plus TCR-seq | Human thymus/PBMC metadata spans infancy to older adults; local thymocyte metadata max 70 years | GEO lists 59 samples; local metadata includes thymus-associated donor/sample labels but ages are not uniform across integrated sheets | Fibroblast-like/mesenchymal, epithelial, mTEC-like, cTEC-like, endothelial candidate labels | Raw archive available, not downloaded automatically | `GSE231906_RAW.tar`, expected 3.69 GB; metadata workbook 17.2 MB | Yes, local metadata workbook | 3 | Most promising human candidate, but expression requires guarded download/extraction and donor/source matching. | Candidate future human donor-level fibroblast-like/epithelial LOX-family age analysis after expression parsing. | Current repository cannot claim human conservation; cannot test exact mouse medFB/capsFB/mTEC1 replication; no protein/activity/functional claims. |

## Context Only

| accession/resource | source | species | data type | age groups | n biological samples/donors if known | stromal/TEC/fibroblast coverage | expression matrix available | approximate file size | processed metadata available | feasibility | reason | claims it can test | claims it cannot test |
|---|---|---|---|---|---:|---|---|---|---|---:|---|---|---|
| CELLxGENE collection `fc19ae6c-d7c1-4dce-b703-62c5d52061b4`, "A spatial human thymus cell atlas mapped to a continuous tissue axis" | CELLxGENE: https://cellxgene.cziscience.com/collections/fc19ae6c-d7c1-4dce-b703-62c5d52061b4; Nature/PubMed-linked article reports fetal through pediatric thymus atlas: https://www.nature.com/articles/s41586-024-07944-6 | Human | scRNA-seq atlas subsets plus spatial data | Fetal to pediatric/early childhood, not old-age aging | Multiple samples; not a clean young-vs-aged adult contrast | Fibroblast, TEC, vascular and whole-thymus atlas assets | Yes through CELLxGENE assets | Subset H5ADs tens to hundreds of MB; whole atlas about 2 GB | Yes | 2 | Useful expression/annotation context, but mainly developmental/pediatric. | Compartment expression context only. | Aged-lower directions; medFB/capsFB/mTEC1 aging validation; human conservation. |
| GSE147520 / E-MTAB-8581 | GEO/ArrayExpress; BioStudies record: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8581/ | Human | 10x scRNA-seq thymus atlas | Fetal, postnatal/pediatric and adult samples, not a clean old-age contrast | GEO primary sample set about 5 donors; ArrayExpress atlas has broader developmental/adult metadata | Rich stromal and TEC labels | Yes | Representative epithelial H5AD about 99 MB; combined H5AD hundreds of MB | Yes | 2 | Good human thymus atlas context, but not a direct aging validation. | Baseline compartment expression context. | Mouse aging directions, exact subtype-resolved external validation, human conservation. |
| SCP3197 / Zenodo-linked "Tissue architecture dynamics underlying immune development and decline in the thymus" | PubMed/supplementary/portal-linked resource; manual portal inspection required | Mouse | Spatial transcriptomics/TCR-focused lifespan thymus resource | Reported 4, 8, 16, 32, 78 weeks | 21 thymus samples reported in resource descriptions | Tissue-level architecture; cell-group extraction depends on spatial annotations | Public resource described; LOX-ready matrix inventory not confirmed locally | Not confirmed; no auto-download | Unclear | 2 | Potentially relevant lifespan context, but not a ready subtype-resolved scRNA validation dataset. | Possible future spatial directional context after manual extraction. | Immediate fibroblast/TEC validation; exact medFB/capsFB/mTEC1 validation. |

## Not Suitable

| accession/resource | source | species | data type | age groups | n biological samples/donors if known | stromal/TEC/fibroblast coverage | expression matrix available | approximate file size | processed metadata available | feasibility | reason | claims it can test | claims it cannot test |
|---|---|---|---|---|---:|---|---|---|---|---:|---|---|---|
| E-MTAB-8737 | ArrayExpress/BioStudies; also referenced by MouseThymusAgeing vignette | Mouse | TEC lineage tracing/single-cell RNA-seq | Age-altered TEC differentiation context, not a clean validation series | 6 sample assays visible from public metadata | TEC-only; no fibroblasts | Top-level MAGE-TAB metadata available; expression retrieval not direct from local search | IDF/SDRF small; expression not confirmed | Yes | 1 | TEC context only, not a straightforward young-vs-aged LOX validation dataset. | Limited TEC marker/context checks. | Fibroblast claims; exact mTEC1 replication; broad aging validation. |
| GSE240017 | GEO | Mouse | Bulk RNA-seq CD45-negative thymic cells | 2-month steady state and post-irradiation time course | Multiple samples but not independent aged contrast | Broad CD45-negative thymic stroma | Yes | Local processed table 4.2 MB if retrieved | Yes | 1 | Same project family as internal data and lacks aged comparison. | None for independent aging validation. | All v5 external validation claims. |
| E-MTAB-11341 | ArrayExpress/BioStudies | Human | 10x spatial transcriptomics | Fetal 11-18 weeks | 15 fetal spatial samples across tissues | Fetal thymus spatial sections | Yes but image/spatial-heavy | TIFFs tens to hundreds of MB each; SDRF small | Yes | 0 | Fetal-only spatial data; no aging direction. | None for current aging claims. | All v5 external validation claims. |
| GSE114651 | GEO | Mouse | Bulk RNA-seq sorted mTEC subsets | Adult pulse-chase lineage tracing, not age series | 4 biological replicates per mTEC subset | mTEC subset context only | Yes | Raw tar about 2.5 MB | Limited | 0 | Useful marker context, but no young-vs-aged contrast. | None for current aging claims. | All v5 external validation claims. |

## Search Conclusion

The most immediately usable external context remains broad and descriptive: GSE223049 for broad fibroblast and broad epithelial directions, and E-MTAB-8560 for TEC/mTEC-like `Loxl2` age-series context. No public dataset found in this second pass provides independent, exact subtype-resolved mouse fibroblast evidence for capsFB `Lox` or medFB `Loxl1`/`Loxl2`.

GSE231906 remains a candidate future human validation dataset only. In the current repository state, only metadata was parsed; no human LOX-family expression analysis was performed and no human conservation claim is supported.
