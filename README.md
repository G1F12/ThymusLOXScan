# ThymusLOXScan

## Abstract
The lysyl oxidase (LOX) family of extracellular matrix crosslinking enzymes 
has been implicated in tissue fibrosis and aging across multiple organs, yet 
its role in thymic involution remains unexplored. Here we analyze LOX-family 
gene expression dynamics during murine thymic aging using a publicly available 
single-cell RNA-seq dataset of 22,932 CD45-negative thymic stromal cells from 
young (2-month) and aged (18-month) C57BL/6 mice. Pseudobulk differential 
expression analysis and single-cell level testing reveal that LOX-family members 
are broadly downregulated across thymic fibroblasts with age, with the strongest 
effects observed for Lox and Loxl1. Fibroblast subtype-resolved analysis uncovers 
spatially divergent isoform dynamics: capsular fibroblasts show a robust reduction 
in Lox (log2FC = -1.46, padj = 7.5e-5), while medullary fibroblasts display 
opposing trajectories for Loxl1 (increased) and Loxl2 (decreased), suggesting 
niche-specific ECM remodeling programs. In medullary thymic epithelial cells, 
Loxl2 is selectively and substantially reduced in aged mice (log2FC = -3.29, 
padj = 9.6e-4). Loxl2 expression correlates with structural ECM components 
Col1a1 and Vim but not with the EMT driver Snai1, pointing to a structural 
rather than fibrotic role. These findings identify LOX-family downregulation 
as a feature of thymic stromal aging and nominate Loxl2 as a candidate mediator 
of age-associated ECM remodeling in both fibroblast and epithelial compartments.

## Data Source
Dataset: GSE240016 (Kousa et al., Nature Immunology 2024)
Species: Mus musculus (female C57BL/6)
Age groups: 2 months (young) vs 18 months (aged)
Cells: 22,932 CD45-negative thymic stromal cells

## Repository Structure
```text
ThymusLOXScan/
+-- README.md
+-- requirements.txt
+-- data/
    +-- processed/
        +-- LOXL2_EMT_spearman_correlations.csv
        +-- LOX_differential_results.csv
        +-- fibroblast_subtype_LOX_breakdown.csv
        +-- mTEC_Lox_Loxl2_stage_correlations.csv
        +-- thymus_annotated.h5ad
        +-- thymus_preprocessed.h5ad
    +-- raw/
        +-- GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad
+-- figures/
    +-- Fig2_effect_size_bubbleplot.png
    +-- Fig_pseudobulk_volcano.png
    +-- Fig_summary_4panel.png
    +-- LOX_summary_figure.png
    +-- existing_annotation_marker_dotplot_by_subset.png
    +-- existing_umap_annotations_stage.png
    +-- fibroblast_subtype_heatmap.png
    +-- mTEC_Lox_Loxl2_scatter.png
    +-- mTEC_Lox_vs_Loxl2_divergence.png
    +-- qc_after_filtering_violin.png
    +-- qc_before_filtering_violin.png
    +-- umap_cell_type_stage.png
    +-- umap_cell_type_subset.png
    +-- LOX_EMT_correlation/
        +-- Loxl2_vs_Col1a1_fibroblasts.png
        +-- Loxl2_vs_Snai1_fibroblasts.png
        +-- Loxl2_vs_Vim_fibroblasts.png
    +-- LOX_basic/
        +-- LOX_dotplot_cTEC_by_stage.png
        +-- LOX_dotplot_endothelial_by_stage.png
        +-- LOX_dotplot_fibroblasts_by_stage.png
        +-- LOX_dotplot_mTEC_by_stage.png
        +-- LOX_violin_all_cells_by_stage.png
        +-- LOX_violin_by_cell_type.png
+-- manuscript/
    +-- discussion.md
    +-- results_draft.md
+-- notebooks/
    +-- 01_scanpy_tutorial_PBMC.ipynb
    +-- 02_data_inspection.ipynb
    +-- 03_preprocessing_thymus.ipynb
    +-- 04_cell_type_annotation.ipynb
    +-- 05_LOX_expression_analysis.ipynb
    +-- 05b_followup_analysis.ipynb
    +-- 06_GSEA_TGFb_analysis.ipynb
+-- results/
    +-- LOX_master_summary.csv
    +-- dropout_analysis_LOX.csv
    +-- magic_imputation_LOX.csv
    +-- pseudobulk_deseq2_LOX.csv
    +-- sc_mannwhitney_FB_combined.csv
    +-- sc_mannwhitney_mTEC1.csv
    +-- sc_spearman_correlations.csv
    +-- pseudobulk_partial/
        +-- 0_arEC.csv
        +-- 10_aaTEC1.csv
        +-- 11_aaTEC2.csv
        +-- 12_cTEC.csv
        +-- 12_early_Pr.csv
        +-- 13_mTEC1.csv
        +-- 14_mTEC-prol.csv
        +-- 15_mTEC2.csv
        +-- 16_mimetic_basal_.csv
        +-- 17_mimetic_tuft_.csv
        +-- 18_mimetic_neuroendo_.csv
        +-- 19_mimetic_goblet_.csv
        +-- 1_capEC.csv
        +-- 20_mimetic_microfold_.csv
        +-- 2_venEC.csv
        +-- 3_capsFB.csv
        +-- 4_intFB.csv
        +-- 5_medFB.csv
        +-- 6_MEC.csv
        +-- 7_vSMC_PC.csv
        +-- 8_nmSC.csv
        +-- 9_Fat.csv
+-- scripts/
    +-- download_data.py
    +-- dropout_analysis_lox.py
    +-- magic_imputation_lox.py
    +-- pseudobulk_deseq2_lox.py
    +-- summarize_results.py
    +-- validate_pipeline.py
```

## Key Results
- Lox and Loxl1 show the strongest age-associated reduction across fibroblasts
  (rank-biserial effect sizes: -0.226 and -0.154 respectively)
- Capsular fibroblasts: Lox log2FC = -1.46 (padj = 7.5e-5)
- Medullary fibroblasts: isoform divergence — Loxl1 up (+0.75), Loxl2 down (-1.00)
- mTEC1: Loxl2 log2FC = -3.29 (padj = 9.6e-4)
- Loxl2 co-expresses with Col1a1 and Vim but not Snai1

## Reproducing the Analysis
1. git clone https://github.com/G1F12/ThymusLOXScan
2. pip install -r requirements.txt
3. Download GSE240016 h5ad from GEO and place in data/raw/
4. Run scripts in order:
   python scripts/pseudobulk_deseq2_lox.py
   python scripts/dropout_analysis_lox.py
   python scripts/magic_imputation_lox.py
   python scripts/summarize_results.py
5. Figures generated by running notebooks 03-05b in order

## Software Versions
anndata 0.12.18 | pydeseq2 0.5.4 | scipy 1.16.1
pandas 2.3.1 | numpy 2.2.6 | matplotlib 3.10.5

## Citation
[PLACEHOLDER — fill after bioRxiv upload]

## Author
Aliaksandr Karatseyeu
