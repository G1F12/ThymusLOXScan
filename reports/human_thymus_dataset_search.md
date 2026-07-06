# Human thymus dataset search

Search date: 2026-07-06

This report summarizes dataset discovery and triage only. It does not provide evidence for human conservation of any mouse finding, and it does not validate LOX-family expression changes in human thymus.

## Sources searched

- GEO and SRA via NCBI accession pages.
- ArrayExpress/BioStudies and ENA links surfaced through HCA/CELLxGENE and article data-availability sections.
- CELLxGENE Discover collections and CELLxGENE public curation API metadata.
- Human Cell Atlas Data Explorer.
- PubMed and primary journal pages.
- Zenodo records linked from primary papers.
- Lab GitHub repositories linked from primary papers.
- Human Protein Atlas/GTEx-style bulk tissue pages as broad context only.

## Search terms used

- human thymus single-cell RNA-seq age
- human thymus epithelial mTEC cTEC scRNA-seq
- human thymus spatial transcriptomics age
- human thymus aging scRNA-seq epithelial
- human thymus LOXL2
- GSE231906 thymus expression matrix
- CELLxGENE human thymus age epithelial
- human thymus stromal fibroblast single cell RNA sequencing
- human thymus developmental aging single cell
- GSE147520 human thymic stroma epithelial cells H5AD
- human thymus CELLxGENE Park epithelial mTEC cTEC
- spatial human thymus atlas CELLxGENE TEC subset

## Top P1/P2 candidates

Detailed triage is in `results/tables/human_thymus_dataset_search_candidates.tsv`.

P1 candidates:

- CELLxGENE/HCA Park thymus atlas: `de13e3e2-23b6-40ed-a413-e9e12d7d3910`, HCA `c1810dbc-16d2-45c3-b45e-3e675f88d87b`, `E-MTAB-8581`. The CELLxGENE TEC subset is about 252 MB and has donor IDs, development-stage labels, sex labels, and cTEC/mTEC/TEC annotations. This is the cleanest first processed matrix for a donor-aware Stage 2 audit.
- `GSE147520`: human thymic stromal scRNA-seq with processed H5AD files, including a 99 MB epithelial-cell file. It includes fetal, postnatal, and one adult specimen and has epithelial, mesenchymal, endothelial, cTEC, and mTEC-relevant annotations in the published analysis.
- CELLxGENE Yayon spatial thymus atlas: `fc19ae6c-d7c1-4dce-b703-62c5d52061b4`. The thymic epithelial subset is about 128.8 MB, with donor IDs, development-stage labels, sex labels, and fine TEC annotations including mTEC type 1/2/3, cTEC, and corticomedullary TEC labels. Fibroblast and vascular subsets are also feasible downloads.

P2 candidates:

- `GSE231906`: the most directly age-relevant human thymus candidate. GEO provides a 3.7 GB processed archive and a 17.2 MB cell-level metadata workbook. The metadata workbook contains age, sex, donor/sample labels, and Epi/Mes/Endo annotations, including cTEC and mTEC-like labels. The expression archive still needs parsing and matrix-to-metadata checks.
- Li et al. 2024 spatial/multi-omics thymus resource: Zenodo `13207776`, with processed RDS objects totaling about 5.1 GB. This is potentially useful spatial context but is large and thymocyte/spatial focused.
- Kamaraj et al. 2026 spatial cartography resource: analyzed data are stated to be available at Zenodo `12595241`, with raw data in BioProject-linked repositories. This is relevant for fetal/pediatric spatial TEC context, but file inventory and metadata fields need a follow-up audit.

## GSE231906 status

`GSE231906` should be treated as P2, not metadata-only. GEO now lists processed data as a 3.7 GB `GSE231906_RAW.tar` archive containing CSV/MTX/TSV files and a separate 17.2 MB `GSE231906_cell-level_metadata.xlsx` file. The metadata workbook was inspected in a temp location only.

Useful metadata fields observed:

- thymus-associated donor/sample IDs such as `donor1`, `donor7-2`, `donor40-1`, and `donor47`;
- sex labels;
- age labels ranging from days/months/years to older adult ages;
- broad labels `Epi`, `Mes`, `Endo`, and `Hema`;
- epithelial labels including `cTEC_lo`, `cTEC_hi`, `mTEC_lo`, `mTEC_hi`, `post_AIRE_mTEC`, `Immature_TEC`, `MKI67+mTEC`, `MKI67+cTEC`, `TEC(neuron)`, `TEC(myo)`, `Ionocyte`, `Ciliated`, and `Tuft`;
- mesenchymal/fibroblast labels including `Fb_1` and `Fb_2`;
- endothelial labels.

Immediate blocker: expression parsing was not completed in Stage 1 because the expression archive is large. Stage 2 should inspect the archive manifest and test one thymus matrix-to-metadata join before any broader analysis.

## CELLxGENE status

Two CELLxGENE collections are strong candidates:

- Park thymus atlas collection `de13e3e2-23b6-40ed-a413-e9e12d7d3910`: immediately usable for a Stage 2 audit, especially its 252 MB human TEC subset.
- Yayon spatial thymus atlas collection `fc19ae6c-d7c1-4dce-b703-62c5d52061b4`: immediately usable for a Stage 2 audit of curated TEC, fibroblast, and vascular subsets, with the TEC subset being the most relevant first file.

CELLxGENE should be audited before large GEO archives because its curated H5AD assets expose donor, development-stage, sex, cell-type, and file-size metadata through a stable public API.

## Rejected and lower-priority datasets

- `GSE195812`: human thymocyte-focused sorted populations; large RDS files and no target epithelial/stromal compartment focus. P3.
- `GSE144870`: postnatal thymocyte-seeding progenitor dataset; useful thymocyte context but lacks target epithelial/stromal compartments. P3.
- `GSE206710`: pediatric thymocyte dataset; no target epithelial/stromal compartments. P3.
- `GSE271304`: postnatal thymocyte CITE-seq; valuable for thymocyte mapping but not a TEC/fibroblast reanalysis target. P3.
- `GSE201718`: mTEC-specific 5'Cap/RNA-seq abundance context with paired immature/mature mTEC donor labels, but age/sex metadata were not evident from the series summary and the assay is not directly comparable to scRNA-seq matrices. P3 until metadata are confirmed.
- Human Protein Atlas/GTEx-style bulk thymus expression: broad tissue context only, not compartment-resolved. P3.
- `GSE215418`: mouse TEC CITE-seq/scRNA-seq, rejected for this human dataset search.

## Stage 2 audit recommendation

Recommended first audit:

1. CELLxGENE Park human thymic epithelial cell H5AD.
2. CELLxGENE Yayon thymic epithelial cell H5AD.
3. `GSE147520_epithelial_cells.h5ad.gz`.
4. `GSE231906` archive manifest and one pilot matrix-to-metadata join.

If the Stage 2 goal is specifically donor-age coverage across pediatric-to-adult thymus, prioritize `GSE231906` after the smaller CELLxGENE/H5AD audit establishes the analysis pattern.

## Reanalysis readiness

At least three public processed-matrix candidates are ready for Stage 2 audit: CELLxGENE Park, CELLxGENE Yayon, and `GSE147520`. `GSE231906` is likely useful but needs manual large-archive handling before it can be called ready for expression reanalysis.

No expression reanalysis was performed in Stage 1. No human conservation, validation, mechanism, causality, protein-level, functional, rejuvenation, or intervention-related conclusion is supported by this report.
