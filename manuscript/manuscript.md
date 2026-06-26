# LOX-family gene expression is reduced in thymic fibroblasts during aging and shows isoform-divergent dynamics in medullary niches

**Author:** Aliaksandr Karatseyeu

**Correspondence:** https://github.com/G1F12/ThymusLOXScan

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

## Results

## Dataset and cell type composition

To investigate age-dependent changes in LOX-family expression within the thymic stroma, we analyzed the publicly available GSE240016 single-cell RNA-seq dataset, comprising 22,932 CD45-negative thymic stromal cells isolated from young (2-month-old) and aged (18-month-old) female C57BL/6 mice (Kousa et al., Nat Immunol, 2024). Cell type annotations were retained from the original study and spanned fibroblasts, thymic epithelial cells (TECs), endothelial cells, vascular smooth muscle cells/pericytes, mesenchymal epithelial cells, and non-myelinating Schwann cells. Fine-grained annotations further resolved fibroblast subtypes — including capsular (capsFB), interstitial (intFB), and medullary (medFB) populations — as well as cortical and medullary TEC states, enabling cell-type-resolved analysis of LOX-family dynamics.

## LOX-family expression is broadly reduced in thymic fibroblasts during aging

Examination of LOX-family expression across the fibroblast compartment revealed a consistent downward shift in aged relative to young cells. Single-cell level analysis (Mann–Whitney U test across all fibroblasts pooled; n = 7,159 young, 6,079 aged cells) showed significant reductions for all five family members: Loxl1 (p = 2.6×10⁻¹²¹, rank-biserial effect size = −0.226), Lox (p = 4.9×10⁻⁹⁸, effect size = −0.154), Loxl2 (p = 3.2×10⁻⁷², effect size = −0.131), Loxl3 (p = 1.4×10⁻²⁶, effect size = −0.070), and Loxl4 (p = 7.7×10⁻⁴, effect size = −0.004). Notably, the Loxl4 effect size was negligible despite nominal significance, reflecting the large cell numbers rather than a biologically meaningful reduction. Pseudobulk DESeq2 analysis, which aggregates counts per biological sample and controls for inter-animal variability, confirmed statistically significant age-dependent changes for individual isoforms within specific fibroblast subtypes (described below), consistent with the direction observed at the single-cell level. The coordinated downregulation of the LOX family in fibroblasts argues against a simple pro-fibrotic remodeling response during thymic involution, and instead suggests a reduction in constitutive ECM crosslinking capacity in the aging thymic stroma.

## Fibroblast subtype analysis reveals spatially divergent LOX dynamics

Disaggregation by fibroblast subtype revealed that the apparent global LOX decline was heterogeneous across anatomically distinct populations. The strongest pseudobulk-confirmed effect was observed in capsular fibroblasts: capsFB exhibited a robust age-associated decrease in Lox expression (log2FC = −1.46, padj = 7.5×10⁻⁵). Interstitial fibroblasts showed consistent directional reductions in both Lox (log2FC = −0.63, padj = 0.072) and Loxl1 (log2FC = −0.51, padj = 0.054), which approached but did not reach FDR significance under the conservative pseudobulk framework.

Medullary fibroblasts displayed a markedly different pattern. Loxl1 was significantly elevated in aged medFB (log2FC = +0.75, padj = 0.023), while Loxl2 was concomitantly reduced (log2FC = −1.00, padj = 0.007). Lox showed a directional trend toward increased expression in aged medFB (log2FC = +1.10, padj = 0.075) that did not survive correction for multiple comparisons. These opposing isoform trajectories within the medullary niche indicate that fibroblast subtype identity, rather than a uniform aging program, determines the direction and magnitude of LOX-family regulation during thymic involution.

## Loxl2 is selectively reduced in medullary TECs during aging

Among thymic epithelial populations, the most pronounced LOX-family change was observed for Loxl2 in medullary TECs (mTEC1 subtype): pseudobulk DESeq2 identified a large and statistically significant reduction in aged relative to young mTECs (log2FC = −3.29, padj = 9.6×10⁻⁴). No other LOX-family member reached significance in this population at either the pseudobulk or single-cell level. Lox showed a directional trend toward upregulation (log2FC = +0.79, padj = 0.130; single-cell Mann–Whitney p = 0.123), consistent in direction with the medFB pattern but not statistically supported. This reduction corresponded to an approximately 9.8-fold decrease in linear expression (2^3.29), representing the largest magnitude change observed across all LOX-family members and cell type combinations tested. The selective reduction of Loxl2 in mTEC1 is the primary finding in the epithelial compartment.

## Loxl2 co-expression with structural ECM components in fibroblasts

Within the fibroblast compartment, Loxl2 expression showed weak but statistically robust co-expression with structural extracellular matrix genes. Significant positive correlations were observed with Col1a1 (Spearman rho = 0.19, p = 3.2×10⁻¹¹¹) and Vim (rho = 0.14, p = 8.0×10⁻⁶³), consistent with a role for Loxl2 in fibroblast-associated matrix remodeling programs. Loxl2 did not significantly correlate with Snai1 (rho = 0.005, p = 0.578), suggesting that its co-expression with Col1a1 and Vim reflects structural ECM organization rather than Snai1-dependent epithelial-to-mesenchymal transition.

## Methods

Analyses were run with Python 3.11.0. MAGIC imputation used magic-impute v3.0.0.

Analysis code is available at https://github.com/G1F12/ThymusLOXScan.

## Discussion
The present reanalysis identified a coordinated reduction of LOX-family gene expression in thymic fibroblasts during aging. In pooled fibroblasts, single-cell Mann-Whitney testing showed significant decreases for Loxl1, Lox, Loxl2, and Loxl3, with the strongest rank-biserial effect observed for Loxl1, whereas Loxl4 showed statistical significance with a negligible effect size. Pseudobulk analysis supported this pattern in specific isoform-subtype combinations, including reduced Lox in capsular fibroblasts and reduced Loxl2 in medullary fibroblasts. This direction was notable because aging in several peripheral organs, including lung, liver, and kidney, has often been associated with fibrotic extracellular matrix stiffening in which LOX-family enzymes are typically upregulated. Thus, thymic fibroblast aging appeared to involve a distinct remodeling state characterized by reduced expression of multiple LOX-family members rather than a canonical fibrotic LOX increase.

The pseudobulk results further indicated subtype-specific and isoform-divergent remodeling within the medullary fibroblast compartment. In medullary fibroblasts, Loxl2 was reduced, whereas Loxl1 was increased, suggesting that aging did not uniformly suppress all LOX-family activity within this niche. Previous studies have shown that LOX-family members can differ in substrate preference, with Loxl1 linked to elastin crosslinking and Loxl2 linked to collagen IV organization. The opposing behavior of Loxl1 and Loxl2 in medullary fibroblasts therefore may have reflected distinct extracellular matrix remodeling programs rather than a single global loss of matrix-modifying capacity. Kousa et al. (2024) reported that thymic fibroblasts upregulated inflammaging programs with age in the same dataset, and the present findings placed LOX-family remodeling alongside that inflammatory shift. However, these data did not establish whether LOX changes preceded, followed, or occurred in parallel with fibroblast inflammaging.

A second notable finding was the strong reduction of Loxl2 in the mTEC1 compartment by pseudobulk analysis. Loxl2 was reduced with a log2 fold change of -3.29, and this change reached adjusted significance. Because Loxl2 has been implicated in basement membrane organization in epithelial tissues, its loss in aging mTECs could compromise the structural integrity of the epithelial-mesenchymal interface and influence maintenance of the thymic epithelial cell niche. This interpretation remained provisional because the mTEC1 pseudobulk comparison included only two biological replicates per age group, making independent validation essential.

Co-expression patterns in fibroblasts supported a structural extracellular matrix interpretation. Loxl2 expression correlated positively with Col1a1 and Vim, whereas it showed no significant association with Snai1. This pattern argued that Loxl2 variation was more closely aligned with matrix and mesenchymal structural maintenance than with an active fibrotic or epithelial-mesenchymal transition-like program. In this respect, thymic aging appeared to differ from fibrotic pathologies in which LOX-family induction can accompany broader collagen remodeling and profibrotic transcriptional programs.

Several limitations constrained the interpretation of these findings. The study was a computational reanalysis of an existing single-cell dataset and lacked direct wet-lab validation of transcript-level changes. The comparison involved a single young and aged timepoint, 2 months versus 18 months, and therefore did not resolve the temporal onset or progression of LOX-family remodeling. The original dataset from Kousa et al. (2024) used female C57BL/6 mice, so the results may not generalize across sex, particularly given established sex differences in thymic aging. Future work should validate Loxl2 protein localization by immunofluorescence at the mTEC basement membrane, test functional consequences of LOX inhibition with BAPN in aged mice, and evaluate whether bulk RNA-seq data from human thymus in GTEx support conserved age-associated changes in LOX-family expression.

## Data Availability
Raw single-cell RNA-seq data analyzed in this study are publicly available 
at the NCBI Gene Expression Omnibus under accession number GSE240016 
(Kousa et al., 2024). All analysis code and intermediate result files 
are available at https://github.com/G1F12/ThymusLOXScan.

## Citation
Karatseyeu A. LOX-family gene expression is reduced in thymic fibroblasts 
during aging and shows isoform-divergent dynamics in medullary niches. 
bioRxiv (2025). https://github.com/G1F12/ThymusLOXScan 
[FILL DOI AFTER BIORXIV UPLOAD]

