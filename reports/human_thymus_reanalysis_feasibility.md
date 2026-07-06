# Human thymus reanalysis feasibility

Overall readiness classification: `ready_for_stage2_audit`

This classification means public datasets were found that are suitable for a focused feasibility audit. It does not mean that any human LOX-family result has been established.

## P1 candidates

### CELLxGENE/HCA Park thymus atlas

Useful because it provides curated H5AD files with donor IDs, development-stage labels, sex labels, and TEC/mTEC/cTEC annotations. The human TEC subset is about 252 MB, which is feasible for a targeted Stage 2 audit.

Blocks immediate analysis:

- confirm whether usable counts are in `raw.X` or `X`;
- confirm exact observation fields for donor, sample/library, age/development stage, sex, and assay;
- check per-donor cell counts for LOX/LOXL1/LOXL2/LOXL3/LOXL4 in each TEC label.

Files to check next:

- CELLxGENE H5AD for `Human Thymic Epithelial Cells - A cell atlas of human thymic development defines T cell repertoire formation`;
- Zenodo `sample_metadata_fix.xlsx` if donor/sample mapping needs cross-checking.

Large download likely: no for the TEC subset; yes if using the full atlas or Zenodo ZIP.

Donor-aware analysis seems possible: yes, pending field and replicate-balance audit.

### GSE147520 human thymic stroma

Useful because GEO provides a small processed epithelial H5AD and all-cell H5ADs. The published dataset targets thymic stromal cells and includes fetal, postnatal, and one adult sample.

Blocks immediate analysis:

- confirm donor/sample fields and sex fields inside H5AD;
- decide whether to use epithelial-only or all-cell object for fibroblast/endothelial context;
- evaluate whether the small donor count is adequate for the intended donor-aware summary.

Files to check next:

- `GSE147520_epithelial_cells.h5ad.gz` (99.0 MB);
- `GSE147520_all_cells.h5ad.gz` (465.0 MB) if non-epithelial stroma are needed.

Large download likely: no for epithelial-only; moderate for all-cell.

Donor-aware analysis seems possible: yes, but with a small donor set and limited older adult coverage.

### CELLxGENE Yayon spatial thymus atlas

Useful because the CELLxGENE TEC subset is about 128.8 MB and includes donor IDs, development-stage labels, sex labels, and fine TEC labels including mTEC type 1/2/3, cTEC, and corticomedullary TEC. Separate fibroblast and vascular subsets are also small enough for focused audit.

Blocks immediate analysis:

- confirm expression layer and normalization status;
- confirm donor and sample/library fields across subset files;
- avoid mixing spatial and single-cell assets without a defined audit plan.

Files to check next:

- CELLxGENE `thymus scRNA-seq atlas - thymic epithelial cell subset`;
- CELLxGENE `thymus scRNA-seq atlas - fibroblast subset`;
- CELLxGENE `thymus scRNA-seq atlas - vascular cell subset`.

Large download likely: no for subset files; yes for the full atlas and multiple spatial files.

Donor-aware analysis seems possible: yes, for fetal/pediatric developmental stages.

## P2 candidates

### GSE231906 human thymus aging dataset

Useful because it has donor/sample IDs, age labels, sex labels, broad Epi/Mes/Endo labels, TEC sublabels, and many thymus samples spanning a broad age range. It is the most age-relevant public candidate found in Stage 1.

Blocks immediate analysis:

- processed expression archive is 3.7 GB;
- expression parsing and matrix-to-metadata joining are unverified;
- archive contents need inspection before deciding whether a full download/extraction is justified.

Files to check next:

- `GSE231906_RAW.tar` manifest or partial listing if possible;
- one pilot MTX/CSV/TSV expression unit from a thymus sample;
- `GSE231906_cell-level_metadata.xlsx` for barcode and annotation mapping.

Large download likely: yes.

Donor-aware analysis seems possible: yes, pending successful expression parsing and barcode/sample mapping.

### Li et al. 2024 spatial/multi-omics thymus resource

Useful because it provides public processed RDS objects for human thymus scRNA/spatial analysis and is linked to primary article data availability.

Blocks immediate analysis:

- processed objects total about 5.1 GB;
- exact donor/sample and TEC/fibroblast annotation fields need inspection;
- primary focus is spatial/thymocyte organization rather than a lightweight stromal expression matrix.

Files to check next:

- `thymus.sc.RDS` (4.4 GB);
- `thymus.st.rds` (667.5 MB);
- GitHub metadata/code to identify object fields before download.

Large download likely: yes.

Donor-aware analysis seems possible: possible, but unverified.

### Kamaraj et al. 2026 spatial cartography resource

Useful because the article reports human fetal and pediatric thymus spatial transcriptomics with TEC, cTEC, mTEC, fibroblast, and endothelial spatial markers and states analyzed data are available through Zenodo.

Blocks immediate analysis:

- Zenodo file inventory and file sizes were not resolved during Stage 1;
- exact donor/sample and matrix formats need audit;
- spatial format may require a different analysis plan than single-cell H5ADs.

Files to check next:

- Zenodo `12595241` analyzed-data files;
- BioProject-linked raw-data records only if processed data are insufficient;
- GitHub notebooks for object schema.

Large download likely: unknown, but likely moderate to large.

Donor-aware analysis seems possible: possible for developmental/fetal-pediatric context, pending metadata audit.

## Readiness summary

- `ready_for_stage2_audit`: yes.
- `possible_but_needs_manual_download`: yes, especially `GSE231906`, Li 2024, and possibly Kamaraj 2026.
- `metadata_only_for_now`: not for the whole search; only applies to datasets whose expression objects are too large or unparsed at Stage 1.
- `no_suitable_dataset_found`: no.

Recommended Stage 2 order: CELLxGENE Park TEC subset, CELLxGENE Yayon TEC subset, `GSE147520` epithelial H5AD, then `GSE231906` large-archive pilot.
