# Human Thymus Expression Feasibility Source Notes

Search date: 2026-07-06

This file records source-level observations from metadata and file inventories only. No expression matrix was parsed.

## CELLxGENE Park TEC

- Accession/collection inspected: CELLxGENE `de13e3e2-23b6-40ed-a413-e9e12d7d3910`; HCA `c1810dbc-16d2-45c3-b45e-3e675f88d87b`; E-MTAB-8581.
- URL or endpoint inspected: `https://api.cellxgene.cziscience.com/curation/v1/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910`; `https://cellxgene.cziscience.com/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910`; `https://explore.data.humancellatlas.org/projects/c1810dbc-16d2-45c3-b45e-3e675f88d87b`; `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8581`.
- What was found: preferred TEC H5AD asset and curated donor, development-stage, sex, assay, and cell-type metadata.
- Processed matrix exists: yes, H5AD asset about 252.1 MB.
- Metadata exists: yes, through CELLxGENE collection API.
- Too large for Stage 2 download: yes, above 100 MB.
- Suitability: suitable for Stage 3 parsing after manual download outside tracked files.

## CELLxGENE Yayon TEC

- Accession/collection inspected: CELLxGENE `fc19ae6c-d7c1-4dce-b703-62c5d52061b4`.
- URL or endpoint inspected: `https://api.cellxgene.cziscience.com/curation/v1/collections/fc19ae6c-d7c1-4dce-b703-62c5d52061b4`; `https://cellxgene.cziscience.com/collections/fc19ae6c-d7c1-4dce-b703-62c5d52061b4`.
- What was found: TEC, fibroblast, and vascular subset H5AD assets with donor, development-stage, sex, assay, and cell-type metadata.
- Processed matrix exists: yes. TEC H5AD about 128.8 MB; fibroblast H5AD about 71.0 MB; vascular H5AD about 89.9 MB.
- Metadata exists: yes, through CELLxGENE collection API.
- Too large for Stage 2 download: TEC subset is above 100 MB; companion subset files are smaller but still not downloaded in this feasibility stage.
- Suitability: suitable for Stage 3 parsing after manual download outside tracked files.

## GSE147520

- Accession/collection inspected: GSE147520; SRP255854; PRJNA624123.
- URL or endpoint inspected: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147520`; `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147520/suppl/`.
- What was found: processed H5AD files and a raw archive in the GEO supplementary directory.
- Processed matrix exists: yes. Epithelial H5AD is about 99.0 MB by binary units but 103,768,949 bytes; all-cell and combined H5AD files are larger.
- Metadata exists: GEO series metadata and likely H5AD obs metadata; age/donor fields were not confirmed from the local Stage 2 inventory.
- Too large for Stage 2 download: yes under the byte threshold rule for expression candidates, and expression files were intentionally not downloaded.
- Suitability: suitable for Stage 3 obs/schema inspection after manual download of the epithelial H5AD.

## GSE231906

- Accession/collection inspected: GSE231906; PRJNA970218.
- URL or endpoint inspected: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231906`; `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/`; `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/GSE231906_cell-level_metadata.xlsx`.
- What was found: a 3.7 GB expression-data archive and a 17.2 MB cell-level metadata workbook.
- Processed matrix exists: likely yes in the large supplementary archive, but archive structure was not inspected locally.
- Metadata exists: yes. Workbook fields include donor/sample, sex, age, broad cell class, and fine epithelial annotations including mTEC/cTEC examples.
- Too large for Stage 2 download: yes, the expression archive is about 3.7 GB.
- Suitability: biologically relevant and donor-aware metadata are present, but expression parsing feasibility remains constrained by archive size and untested barcode joins.
