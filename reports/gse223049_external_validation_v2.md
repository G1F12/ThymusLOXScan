# GSE223049 external validation v2

## Scope

GSE223049 is used here as an independent, broad sorted bulk RNA-seq comparison of young 2-month and aged 22-24-month mouse thymic samples. The relevant sorted populations are thymic fibroblasts and thymic epithelial cells.

All statistics in this report are descriptive. This dataset is not subtype-resolved and therefore cannot validate capsFB, medFB, or mTEC1 specificity. It should be interpreted as broad directional external context, not as a replacement for subtype-resolved validation.

## Directional consistency

| external comparison | young n | aged n | aged-minus-young delta log2(CPM+1) | descriptive p | interpretation |
|---|---:|---:|---:|---:|---|
| Thymic fibroblast Lox | 4 | 4 | -0.354 | 0.114 | Directionally consistent with broad/capsFB aged-lower Lox, but only broad fibroblast support. |
| Thymic fibroblast Loxl1 | 4 | 4 | -0.112 | 1.000 | Not consistent with a medFB Loxl1 increase; medFB Loxl1 increase is not validated. |
| Thymic fibroblast Loxl2 | 4 | 4 | -1.206 | 0.029 | Directionally consistent with aged-lower fibroblast/medFB Loxl2, but only broad fibroblast support. |
| Thymic epithelial Loxl2 | 5 | 5 | -0.588 | 0.151 | Directionally consistent with aged-lower epithelial/mTEC1 Loxl2, but not mTEC1-specific. |

## Findings not testable in GSE223049

- capsFB specificity is not testable because thymic fibroblasts are pooled in sorted bulk RNA-seq.
- medFB specificity is not testable because thymic fibroblast subtypes are not separated.
- mTEC1 specificity is not testable because thymic epithelial cells are pooled.
- Detection-rate or single-cell subtype mechanisms are not testable from this bulk count matrix.
- The medFB Loxl1 increase is not validated; broad thymic fibroblast Loxl1 is slightly lower in aged samples in this dataset.

## Figures

- `results/figures/external_validation/gse223049/gse223049_thymic_fibroblast_lox_dot_box.png`
- `results/figures/external_validation/gse223049/gse223049_thymic_fibroblast_lox_dot_box.pdf`
- `results/figures/external_validation/gse223049/gse223049_thymic_epithelial_loxl2_dot_box.png`
- `results/figures/external_validation/gse223049/gse223049_thymic_epithelial_loxl2_dot_box.pdf`
- `results/figures/external_validation/gse223049/gse223049_lox_family_delta_summary.png`
- `results/figures/external_validation/gse223049/gse223049_lox_family_delta_summary.pdf`

## Bottom line

GSE223049 supports broad external directionality for aged-lower thymic fibroblast Lox, aged-lower thymic fibroblast Loxl2, and aged-lower thymic epithelial Loxl2. It does not validate subtype-specific capsFB, medFB, or mTEC1 claims, and it does not validate the medFB Loxl1 increase.