# ThymusLOXScan

Computational reanalysis of LOX-family transcript patterns in aging murine thymic stromal cells.

This repository analyzes public single-cell RNA-seq data from CD45-negative murine thymic stroma to evaluate whether LOX-family transcript changes vary by stromal subtype during aging. The analysis uses biological-sample-aware pseudobulk summaries as the main inferential framework, with descriptive single-cell, detection-rate, per-sample, dropout/depth, annotation-sanity, and external directional-context analyses as supporting checks. The interpretation is hypothesis-generating and does not claim protein-level validation, causality, enzymatic activity, ECM crosslinking changes, thymic rejuvenation, therapeutic relevance, human conservation, or exact subtype-resolved external validation.

## Current Manuscript And Release

Current manuscript:

`manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`

Current release:

`v5.2-dropout-protein-feasibility`

Previous releases:

- `v5.0-external-validation`: previous external-validation update.
- `v4.2-final-safe`: previous stable cautious release.

Older manuscript files are retained as working drafts or historical versions.

## Key Outputs

- Current public manuscript: `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`
- One-page collaborator summary: `one_page_summary.md`
- Wet-lab validation plan: `wet_lab_validation_plan.md`
- Main pseudobulk table: `results/tables/lox_pseudobulk_complete_results.tsv`
- Cross-dataset direction matrix: `results/tables/cross_dataset_lox_validation_matrix.tsv`
- Dropout/depth audit: `reports/mtec1_loxl2_dropout_depth_audit.md`
- LOXL2 protein/spatial feasibility audit: `reports/loxl2_protein_spatial_validation_feasibility.md`
- Final figures: `results/figures/final/`
- Supplementary tables: `supplementary_tables/`

## Dataset

- Accession: `GSE240016`
- Dataset: Kousa et al. thymic aging single-cell RNA-seq dataset
- Species: mouse
- Comparison used here: young `02mo` versus aged `18mo`
- Main raw input expected by this repository:
  - `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad`

The full raw `.h5ad` file may be too large to redistribute in a lightweight code repository. Reviewers and collaborators should download the public data from GEO accession `GSE240016` or the linked public data source for the study, then place the annotated CD45-negative thymic stroma AnnData file at the path above. If the file name differs after download, rename it or update script input paths accordingly.

## What This Project Does Not Claim

- It does not establish protein abundance, LOXL2 protein localization, secretion, enzymatic activity, or ECM crosslinking changes.
- It does not establish causality, mechanism, thymic functional effects, rejuvenation, or therapeutic relevance.
- It does not claim that dropout/depth artifacts are ruled out for sparse mTEC1 `Loxl2`.
- It does not claim human conservation or exact subtype-resolved external validation.
- It does not replace direct wet-lab validation in independent cohorts.

## Reproducibility

The notebooks record the exploratory workflow, while the scripts regenerate the main reviewer-facing outputs. Run commands from the repository root after preparing the environment and placing the raw input file in `data/raw/`.

Recommended setup:

```bash
conda env create -f environment.yml
conda activate thymus-loxscan
```

Equivalent pip/venv setup:

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

`scripts/validate_pipeline.py` includes older checks for some legacy processed outputs and should be treated as a repository sanity check rather than the definitive manuscript pipeline.

## Main Analysis Commands

Prepare the expected directories:

```bash
mkdir -p data/raw data/processed results results/tables results/figures/final
```

Run the main scripted analyses:

```bash
python scripts/pseudobulk_deseq2_lox.py
python scripts/make_pseudobulk_results_table.py
python scripts/dropout_analysis_lox.py
python scripts/magic_imputation_lox.py
python scripts/summarize_results.py
python scripts/plot_per_sample_lox_expression.py
python scripts/figures/plot_final_volcano.py
python scripts/figures/plot_final_summary.py
python supplementary_tables/make_supplementary_tables.py
```

Additional v5.x analyses include external directional-context scripts, annotation-sanity checks, fibroblast composition sensitivity analysis, mTEC1 `Loxl2` dropout/depth audit, and LOXL2 protein/spatial feasibility summaries. The relevant scripts are in `scripts/`, `scripts/figures/`, and `supplementary_tables/`; associated outputs are in `results/`, `supplementary_tables/`, and `reports/`.

MAGIC/imputation outputs are visualization or qualitative sensitivity outputs only and are not the primary inferential evidence.

## Repository Structure

```text
ThymusLOXScan/
|-- data/                         # local raw and processed data; large files are not tracked
|-- notebooks/                    # exploratory preprocessing, annotation, and analysis notebooks
|-- scripts/                      # reproducible analysis and validation scripts
|-- results/
|   |-- tables/                   # pseudobulk, external-context, and audit tables
|   |-- figures/final/            # final manuscript figures
|   |-- figures/per_sample/       # per-sample checks
|-- supplementary_tables/         # supplementary TSV tables and generator
|-- reports/                      # audit, methods, sensitivity, and safety reports
|-- manuscript/                   # manuscript text and historical drafts
|-- requirements.txt              # base Python dependencies
|-- environment.yml               # conda environment
```

## Links

- One-page summary: `one_page_summary.md`
- Wet-lab validation plan: `wet_lab_validation_plan.md`
- Current v5.2 manuscript: `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`
- Releases: https://github.com/G1F12/ThymusLOXScan/releases
- Current release: https://github.com/G1F12/ThymusLOXScan/releases/tag/v5.2-dropout-protein-feasibility
- Previous external-validation update: https://github.com/G1F12/ThymusLOXScan/releases/tag/v5.0-external-validation
- Previous stable cautious release: https://github.com/G1F12/ThymusLOXScan/releases/tag/v4.2-final-safe

## Citation

Preprint in preparation.
