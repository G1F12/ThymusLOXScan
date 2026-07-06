# Human Thymus Expression Feasibility Audit

Search date: 2026-07-06

This Stage 2 audit inspected metadata endpoints, file inventories, and small metadata files only. It did not run LOX-family expression analysis, did not download expression matrices into the repository, and does not support any human expression or conservation conclusion.

## Datasets Audited

| Candidate | Accession or portal | Stage 2 classification |
|---|---|---|
| CELLxGENE Park TEC subset | CELLxGENE de13e3e2-23b6-40ed-a413-e9e12d7d3910; HCA c1810dbc-16d2-45c3-b45e-3e675f88d87b; E-MTAB-8581 | possible_but_needs_manual_download |
| CELLxGENE Yayon TEC subset | CELLxGENE fc19ae6c-d7c1-4dce-b703-62c5d52061b4 | possible_but_needs_manual_download |
| GSE147520 thymic stromal data | GSE147520; SRP255854; PRJNA624123 | possible_but_needs_manual_download |
| GSE231906 aging thymus data | GSE231906; PRJNA970218 | possible_but_large_data |

## CELLxGENE Park Status

CELLxGENE collection metadata were reachable through the public collection API. The preferred dataset is "Human Thymic Epithelial Cells - A cell atlas of human thymic development defines T cell repertoire formation".

- Expression matrix availability: processed H5AD asset recorded, about 252.1 MB.
- Donor/sample ID availability: donor identifiers are exposed in CELLxGENE metadata.
- Age or development-stage metadata availability: development-stage labels are exposed in CELLxGENE metadata.
- Epithelial/mTEC/cTEC annotation availability: curated cell-type labels include thymic epithelial compartments and mTEC/cTEC-relevant labels.
- Batch/run/library metadata availability: assay and dataset identifiers are present; run-level fields still need local H5AD obs inspection.
- Download burden: moderate. Manual Stage 3 download is required outside tracked files.
- Stage 2 blocker: expression layer, raw/count availability, and exact obs schema were not inspected locally.

## CELLxGENE Yayon Status

CELLxGENE collection metadata were reachable through the public collection API. The relevant subset files include thymic epithelial, fibroblast, and vascular subsets.

- Expression matrix availability: processed H5AD assets recorded. TEC subset is about 128.8 MB; fibroblast subset is about 71.0 MB; vascular subset is about 89.9 MB.
- Donor/sample ID availability: donor identifiers are exposed in CELLxGENE metadata.
- Age or development-stage metadata availability: development-stage labels are exposed in CELLxGENE metadata.
- Epithelial/mTEC/cTEC annotation availability: curated cell-type labels include thymic epithelial, mTEC, and cTEC-relevant compartments.
- Batch/run/library metadata availability: assay and dataset identifiers are present; run-level fields still need local H5AD obs inspection.
- Download burden: moderate. Manual Stage 3 download is required outside tracked files.
- Stage 2 blocker: expression layer and harmonized metadata across subset files were not inspected locally.

## GSE147520 Status

GEO supplementary-file inventory was reachable. The inventory includes processed H5AD files relevant to thymic stromal and epithelial analysis.

- Expression matrix availability: processed H5AD files recorded, including `GSE147520_epithelial_cells.h5ad.gz` at about 99.0 MB, `GSE147520_all_cells.h5ad.gz` at about 465.0 MB, and `GSE147520_Bautista_Park_combined.h5ad.gz` at about 398.4 MB.
- Donor/sample ID availability: likely available through sample-level metadata or H5AD obs, but not confirmed in a local expression object.
- Age or development-stage metadata availability: not confirmed from the Stage 2 metadata fetch; exact H5AD obs or sample-table mapping still needs inspection.
- Epithelial/mTEC/cTEC annotation availability: expected from processed epithelial H5AD, but exact annotation fields were not inspected locally.
- Batch/run/library metadata availability: GEO and H5AD-derived sample metadata are expected to be sufficient for an audit, but local obs inspection is needed.
- Download burden: low to moderate for the epithelial H5AD, larger for all-cell or combined files. Manual Stage 3 download is required outside tracked files.
- Stage 2 blocker: exact donor/sample, age, and annotation fields inside the H5AD were not inspected locally.

## GSE231906 Status

GEO supplementary-file inventory and the small metadata workbook were reachable. The metadata workbook was read from the operating-system temporary directory only and was not written into the repository.

- Expression matrix availability: the supplementary directory includes `GSE231906_RAW.tar`, about 3.7 GB. Stage 1 and Stage 2 treat this as the large expression-data archive to inspect later; it was not downloaded.
- Donor/sample ID availability: workbook fields include donor/sample identifiers such as `geo_sample_id`.
- Age or development-stage metadata availability: workbook fields include age labels with examples spanning pediatric and adult donors.
- Epithelial/mTEC/cTEC annotation availability: workbook fields include broad Epi labels and fine labels including mTEC and cTEC examples.
- Batch/run/library metadata availability: source and sample fields are available or inferable from workbook and archive metadata.
- Download burden: large. A controlled Stage 3 archive inspection is required outside tracked files.
- Stage 2 blocker: matrix layout, barcode linkage, and matrix-to-metadata joins remain untested.

## Matrix and Metadata Summary

| Candidate | Matrix status | Metadata status | Donor/age status | TEC annotation status |
|---|---|---|---|---|
| CELLxGENE Park TEC | curated_h5ad_available_metadata_only_not_downloaded | usable_metadata_inventory | available | available |
| CELLxGENE Yayon TEC | curated_h5ad_available_metadata_only_not_downloaded | usable_metadata_inventory | available | available |
| GSE147520 | processed_h5ad_available_metadata_only_not_downloaded | usable_metadata_inventory | likely/partial until local obs inspection | expected but not locally confirmed |
| GSE231906 | large_processed_archive_available_not_parsed | usable_metadata_inventory | available | available in metadata |

## Final Recommendation for Stage 3

- proceed_to_stage3_expression_parsing: yes
- recommended_dataset_for_stage3: CELLxGENE Park TEC subset first; CELLxGENE Yayon TEC subset second; GSE147520 epithelial H5AD third; GSE231906 only as a controlled large-data parsing pilot.
- recommended_stage3_order: 1. CELLxGENE Park TEC H5AD; 2. CELLxGENE Yayon TEC H5AD; 3. GSE147520 epithelial H5AD; 4. GSE231906 archive manifest and one join test.
- manual_download_required: yes
- expected_main_blocker: manual expression-file download outside tracked files, followed by local obs/layer/schema inspection; GSE231906 additionally has a 3.7 GB archive burden and untested barcode joins.
- confidence: moderate-high for CELLxGENE metadata feasibility; moderate for GSE147520; moderate for GSE231906 metadata feasibility but lower for expression parsing until archive structure is inspected.

No human expression conclusion is made in this stage.
