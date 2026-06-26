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
# Results

## Dataset and cell type composition
The analysis used the annotated GSE240016 single-cell RNA-seq dataset, comprising 22,932 CD45-negative thymic stromal cells from young (02mo) and aged (18mo) C57BL/6 mice [REF]. The dataset contained existing annotations spanning fibroblasts (FB), thymic epithelial cells (TEC), endothelial cells (EC), vascular smooth muscle/pericyte-like cells, mesenchymal epithelial cells, and non-myelinating Schwann cells. Fine annotations resolved fibroblast, endothelial, and TEC subsets, including capsular, interstitial, and medullary fibroblasts as well as cortical and medullary TEC states.

## LOX-family expression is broadly reduced in thymic fibroblasts during aging
LOX-family expression was broadly reduced in thymic fibroblasts during aging. All five assayed genes showed lower mean expression in 18mo fibroblasts relative to 02mo fibroblasts, with the strongest statistical evidence observed for Loxl1 (p = 2.6e-121, rank-biserial effect size = -0.226) and Lox (p = 4.9e-98, effect size = -0.154). Loxl2 and Loxl3 were also significantly reduced with age (p = 3.2e-72 and p = 1.4e-26, respectively), whereas Loxl4 showed a smaller but significant decrease (p = 7.7e-4). These results indicated that aging was associated with a coordinated reduction of lysyl oxidase family expression in the fibroblast compartment rather than a generalized fibroblast-associated induction of matrix-remodeling enzymes.

## Fibroblast subtype analysis reveals spatially divergent LOX dynamics
Fibroblast subtype analysis revealed that the age-associated LOX-family pattern was spatially heterogeneous across annotated fibroblast states. Capsular fibroblasts and interstitial fibroblasts showed reduced Lox and Loxl1 expression in 18mo mice, with capsular fibroblasts showing a pronounced reduction in Lox (log2FC = -1.21) and interstitial fibroblasts showing reductions in both Lox (log2FC = -0.53) and Loxl1 (log2FC = -0.37). In contrast, medullary fibroblasts showed increased Lox and Loxl1 expression with age (log2FC = 0.86 and 0.78, respectively), despite reduced Loxl2 and Loxl4. These data suggested that the global fibroblast LOX decline reflected subtype-specific remodeling rather than a uniform response across fibroblast niches (Table: fibroblast_subtype_LOX_breakdown.csv).

## Lox and Loxl2 show opposing age-dependent dynamics in medullary TECs
Medullary TECs displayed isoform-specific LOX-family remodeling with age. Lox expression increased in aged mTECs (p = 2.2e-7, effect size = 0.047), whereas Loxl2 decreased substantially (p = 1.8e-13, effect size = -0.071). This opposing behavior indicated that mTEC aging was not characterized by uniform LOX-family activation or repression, but instead by divergent regulation of individual isoforms. Consistent with this divergence, Lox and Loxl2 were not significantly correlated in 02mo mTECs (rho = 0.028, p = 0.248), whereas a weak but significant negative correlation emerged in 18mo mTECs (rho = -0.059, p = 0.017).

## Loxl2 co-expression with structural ECM genes in fibroblasts
In fibroblasts, Loxl2 expression showed weak but highly significant co-expression with structural extracellular matrix genes. Loxl2 correlated with Col1a1 (rho = 0.19, p = 3.2e-111) and Vim (rho = 0.14, p = 8.0e-63), supporting a relationship between Loxl2 and fibroblast matrix/mesenchymal programs. In contrast, Loxl2 did not significantly correlate with Snai1 (rho = 0.005, p = 0.578), suggesting that the association was more consistent with structural ECM remodeling than with a Snai1-linked EMT program.

## Methods
# Methods

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

