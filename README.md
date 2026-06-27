# ThymusLOXScan

Computational reanalysis of LOX-family gene expression in aging murine thymic stromal cells.

**Final public manuscript:** `manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md`

Older manuscript files are working drafts.

## Scientific Summary

This repository analyzes LOX-family gene expression in CD45-negative thymic stromal cells from young and aged mice. The analysis focuses on pseudobulk differential expression across biological samples, descriptive single-cell tests, per-sample robustness checks, and manuscript-ready figures/tables for LOX-family genes (`Lox`, `Loxl1`, `Loxl2`, `Loxl3`, `Loxl4`). The current manuscript interpretation is hypothesis-generating and does not claim protein-level validation or direct functional mechanism.

## Dataset

- Accession: `GSE240016`
- Dataset: Kousa et al. thymic aging single-cell RNA-seq dataset
- Species: mouse
- Comparison used here: young `02mo` versus aged `18mo`
- Main raw input expected by this repository:
  - `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad`

The full raw `.h5ad` file may be too large to redistribute in a lightweight code repository. External reviewers should download the public data from GEO accession `GSE240016` or the linked public data source for the study, then place the annotated CD45-negative thymic stroma AnnData file at the path above. If the file name differs after download, rename it or update script input paths accordingly.

## Repository Structure

```text
ThymusLOXScan/
|-- data/
|   |-- raw/                         # downloaded public input data
|   |-- processed/                   # processed AnnData objects and intermediate CSVs
|-- notebooks/                       # exploratory preprocessing, annotation, and analysis notebooks
|-- scripts/                         # reproducible analysis and validation scripts
|   |-- figures/                     # final figure-generation scripts
|-- results/
|   |-- tables/                      # reviewer-facing pseudobulk result tables
|   |-- figures/final/               # final manuscript figure PNG/PDF files
|   |-- figures/per_sample/          # per-sample robustness figures
|   |-- pseudobulk_partial/          # partial pseudobulk DESeq2 outputs
|-- supplementary_tables/            # supplementary TSV tables and generator
|-- reports/                         # audit, methods, and robustness reports
|-- manuscript/                      # manuscript text and revised sections
|-- figures/                         # older notebook-generated figures
|-- requirements.txt                 # base Python dependencies
```

## Required Software

Recommended:

- Python 3.11
- R is not required for the current Python-based DESeq2 workflow, which uses `pydeseq2`
- JupyterLab or Jupyter Notebook for running notebooks

Core packages are specified in both `environment.yml` and `requirements.txt`:

- `scanpy`
- `anndata`
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `leidenalg`
- `python-igraph`
- `pydeseq2`
- `magic-impute`
- `statsmodels`
- `gseapy`
- `jupyter`

The repository does not use R, `rpy2`, or an R/Bioconductor DESeq2 installation in the current scripted workflow.

## Environment Setup

Preferred conda/mamba setup:

```bash
conda env create -f environment.yml
conda activate thymus-loxscan
```

Equivalent mamba setup:

```bash
mamba env create -f environment.yml
mamba activate thymus-loxscan
```

Pip/venv setup:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Optional validation:

```bash
python scripts/validate_pipeline.py
```

Note: `scripts/validate_pipeline.py` contains older checks that expect some legacy processed outputs and may warn about human-style LOX symbols. Treat validation output as a sanity check, not as the definitive manuscript pipeline.

## Data Preparation

Expected raw input:

```text
data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad
```

To prepare the repository:

```bash
mkdir -p data/raw data/processed results results/tables results/figures/final
```

Then download the public `GSE240016` thymic stromal AnnData file from GEO or the linked public source and place it in `data/raw/` with the expected filename.

Important: `scripts/download_data.py` currently references an older or incorrect accession in its internal text and should not be used as the authoritative download route for this manuscript without manual confirmation. The authoritative accession for this project is `GSE240016`.

## Full Pipeline

The notebooks record the original exploratory workflow, while the scripts regenerate the main reviewer-facing outputs. Run commands from the repository root.

1. Inspect raw data:

```bash
jupyter notebook notebooks/02_data_inspection.ipynb
```

2. Preprocess and annotate cells:

```bash
jupyter notebook notebooks/03_preprocessing_thymus.ipynb
jupyter notebook notebooks/04_cell_type_annotation.ipynb
```

Expected processed outputs include:

```text
data/processed/thymus_preprocessed.h5ad
data/processed/thymus_annotated.h5ad
```

3. Run original LOX exploratory analyses:

```bash
jupyter notebook notebooks/05_LOX_expression_analysis.ipynb
jupyter notebook notebooks/05b_followup_analysis.ipynb
```

4. Run pseudobulk DESeq2 analysis:

```bash
python scripts/pseudobulk_deseq2_lox.py
python scripts/make_pseudobulk_results_table.py
```

The complete pseudobulk table script sums raw counts within biological sample and annotation group before fitting DESeq2 with design `~ stage` and contrast `18mo` versus `02mo`.

5. Run descriptive single-cell and auxiliary analyses:

```bash
python scripts/dropout_analysis_lox.py
python scripts/magic_imputation_lox.py
python scripts/summarize_results.py
```

MAGIC-related outputs should be interpreted cautiously. In the current manuscript framing, MAGIC is visualization/supporting analysis only and should not be treated as the primary inferential test.

6. Generate per-sample robustness plots:

```bash
python scripts/plot_per_sample_lox_expression.py
```

7. Generate final manuscript figures:

```bash
python scripts/figures/plot_final_volcano.py
python scripts/figures/plot_final_summary.py
```

8. Generate supplementary tables:

```bash
python supplementary_tables/make_supplementary_tables.py
```

## Regenerating Tables

Main pseudobulk LOX result tables:

```bash
python scripts/make_pseudobulk_results_table.py
```

Outputs:

```text
results/tables/lox_pseudobulk_complete_results.csv
results/tables/lox_pseudobulk_complete_results.tsv
```

Supplementary tables:

```bash
python supplementary_tables/make_supplementary_tables.py
```

Outputs:

```text
supplementary_tables/Supplementary_Table_1_cell_counts.tsv
supplementary_tables/Supplementary_Table_2_pseudobulk_LOX_results.tsv
supplementary_tables/Supplementary_Table_3_single_cell_tests.tsv
supplementary_tables/Supplementary_Table_4_correlations.tsv
```

Per-sample robustness values:

```text
results/figures/per_sample/lox_per_sample_pseudobulk_values.csv
results/tables/mtec1_loxl2_per_sample_expression.tsv
```

## Regenerating Figures

Final manuscript figures:

```bash
python scripts/figures/plot_final_volcano.py
python scripts/figures/plot_final_summary.py
```

Outputs:

```text
results/figures/final/Fig1_pseudobulk_volcano_final.png
results/figures/final/Fig1_pseudobulk_volcano_final.pdf
results/figures/final/Fig2_summary_4panel_final.png
results/figures/final/Fig2_summary_4panel_final.pdf
```

Per-sample robustness figures:

```bash
python scripts/plot_per_sample_lox_expression.py
```

Outputs:

```text
results/figures/per_sample/lox_per_sample_pseudobulk_expression.png
results/figures/per_sample/lox_per_sample_pseudobulk_expression.pdf
results/figures/per_sample/mtec1_loxl2_per_sample.pdf
```

Older exploratory figures are generated from notebooks and stored under `figures/`.

## Manuscript Files

The final public manuscript for this release is:

```text
manuscript/LOX_thymus_aging_public_preprint_v4_2_final_safe.md
```

Manuscript-related working files are stored in `manuscript/`, including:

```text
manuscript/manuscript.md
manuscript/manuscript.docx
manuscript/methods.md
manuscript/methods_detailed.md
manuscript/results_draft.md
manuscript/discussion.md
manuscript/discussion_revised.md
manuscript/abstract_revised.md
manuscript/data_availability.md
```

References and remaining citation gaps:

```text
manuscript/references.bib        # create or include when final references are added
reports/reference_gaps.md        # create or include when citation-gap tracking is available
```

There is currently no automated manuscript-build script. Manuscript outputs should be edited or exported from the files in `manuscript/`.

## Known Limitations

- The analysis is a computational reanalysis of a public dataset and does not include new experimental validation.
- RNA expression changes do not establish protein abundance, enzymatic activity, ECM crosslinking changes, or causal mechanisms of thymic involution.
- Some subtype-level pseudobulk comparisons may have limited biological replicate counts; these should be interpreted cautiously.
- Single-cell-level Mann-Whitney tests are descriptive and do not replace sample-level pseudobulk inference.
- MAGIC/imputation outputs are not the primary inferential evidence.
- Some older notebooks and helper scripts are exploratory and may contain legacy assumptions. The reviewer-facing scripted outputs are the preferred reproducibility path.
- The raw public `.h5ad` file may need to be downloaded manually because redistribution and file hosting depend on the original data provider.

## Citation

Preprint in preparation.
