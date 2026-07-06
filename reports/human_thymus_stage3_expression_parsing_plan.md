# Human Thymus Stage 3 Expression Parsing Plan

This is a plan only. No LOX-family human expression analysis was performed in Stage 2.

## Selected Datasets

1. CELLxGENE Park TEC subset: `de13e3e2-23b6-40ed-a413-e9e12d7d3910`.
2. CELLxGENE Yayon TEC subset: `fc19ae6c-d7c1-4dce-b703-62c5d52061b4`.
3. GSE147520 epithelial H5AD: `GSE147520_epithelial_cells.h5ad.gz`.
4. GSE231906: large archive only after smaller candidates are parsed and only outside tracked files.

## Manual Download Targets

- Park TEC H5AD: `https://datasets.cellxgene.cziscience.com/0af6c752-f1f4-44b3-8186-32f590fdd531.h5ad`.
- Yayon TEC H5AD: `https://datasets.cellxgene.cziscience.com/ad4813c6-33ce-48e6-8505-8bbaca5408c5.h5ad`.
- Optional Yayon fibroblast H5AD: `https://datasets.cellxgene.cziscience.com/6e0c3112-bbaa-4702-a2bd-1fe921ab1cc6.h5ad`.
- Optional Yayon vascular H5AD: `https://datasets.cellxgene.cziscience.com/f483be56-16e7-42ee-b476-132854b14813.h5ad`.
- GSE147520 epithelial H5AD: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147520/suppl/GSE147520_epithelial_cells.h5ad.gz`.
- GSE231906 large archive, deferred: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/GSE231906_RAW.tar`.

## Local Path Convention

Use:

```text
data/external/human_thymus/<accession>/
```

Examples:

```text
data/external/human_thymus/cellxgene_park_tec/
data/external/human_thymus/cellxgene_yayon_tec/
data/external/human_thymus/GSE147520/
data/external/human_thymus/GSE231906/
```

Files under `data/external` should remain untracked unless a tiny metadata file is intentionally allowed. Large expression files, H5AD files, archives, and derived processed data should not be committed.

## Target Genes for Later Analysis

- LOX
- LOXL1
- LOXL2
- LOXL3
- LOXL4

## Target Compartments for Later Analysis

- epithelial / TEC-like
- mTEC-like
- cTEC-like
- fibroblast-like / mesenchymal
- endothelial as optional context

## Stage 3 Parsing Steps

1. Download the Park TEC H5AD manually to the untracked external-data path.
2. Inspect `obs`, `var`, available layers, raw matrix availability, donor/development-stage fields, and cell-type fields.
3. Confirm that target genes exist in `var` using stable gene symbols or mapped feature names.
4. Produce a donor-aware pseudobulk feasibility table before any expression interpretation.
5. Repeat for the Yayon TEC H5AD and optional fibroblast/vascular subset files.
6. Inspect GSE147520 epithelial H5AD obs schema and compare donor/sample and age fields.
7. Defer GSE231906 until a controlled archive manifest or small pilot extraction can be performed outside tracked files.

## Donor-Aware Pseudobulk Plan

- Group cells by dataset, donor/sample ID, age or development-stage label, sex if available, and compartment/subtype.
- Summarize target-gene counts or normalized values only after confirming matrix layer semantics.
- Keep per-donor sample counts and cell counts in output tables.
- Flag groups with too few cells or unclear donor labels.
- Avoid cross-dataset harmonization until each dataset has a local schema audit.

## No-Overclaim Wording

Use language such as "human dataset parsing feasibility", "candidate expression context", and "transcript-level exploratory reanalysis". Do not describe Stage 3 outputs as conservation evidence or as functional, protein-level, activity, causal, or intervention evidence.

## Expected Stage 3 Outputs

- `results/tables/human_thymus_stage3_dataset_schema.tsv`
- `results/tables/human_thymus_stage3_gene_presence.tsv`
- `results/tables/human_thymus_stage3_pseudobulk_input_summary.tsv`
- `results/tables/human_thymus_stage3_expression_summary.tsv`
- `reports/human_thymus_stage3_expression_parsing_report.md`
- `reports/human_thymus_stage3_safety_check.md`
