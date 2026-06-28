# E-MTAB-8560 TEC external validation

## Scope

E-MTAB-8560 is a mouse TEC Smart-seq2 single-cell age-series spanning 1, 4, 16, 32, and 52 weeks. This script uses public ArrayExpress/BioStudies metadata plus processed Bioconductor MouseThymusAgeing SMART-seq2 resources.

Cells are summarized descriptively. The parsed metadata provides cell IDs, age, sort type, sort day, and subtype labels, but not a clean donor-level design suitable for treating cells as independent biological replicates. Spearman and age-bin comparisons therefore use sort_day x age summaries and should be read as trend checks only.

## Loxl2 interpretation

| group | cells | sort-day samples | oldest-minus-youngest mean delta | Spearman rho | direction |
|---|---:|---:|---:|---:|---|
| broad TEC | 2327 | 5 | -1.041 | -0.820 | aged-lower |
| mTEC-like | 1735 | 5 | -1.287 | -0.800 | aged-lower |
| mTEClo | 582 | 5 | -1.529 | -0.810 | aged-lower |
| mTEChi | 715 | 5 | -0.859 | -0.479 | aged-lower |

Overall, E-MTAB-8560 supports a broad aged-lower TEC/mTEC-like tendency, but cannot specifically validate the current mTEC1 Loxl2 finding.

## Relationship to current mTEC1 Loxl2

- This dataset can test broad TEC, cTEC, mTEClo, mTEChi, and mTEC-like age trends.
- It cannot directly test the current mTEC1 label because the annotation systems do not match one-to-one.
- The age range ends at 52 weeks, whereas the current aged comparison is older; this may capture maturation and midlife remodeling more than late thymic involution.
- Any agreement should be framed as directional external context, not subtype-specific validation.

## Figures

- `results/figures/external_validation/emtab8560/emtab8560_lox_family_expression_by_age.png`
- `results/figures/external_validation/emtab8560/emtab8560_lox_family_expression_by_age.pdf`
- `results/figures/external_validation/emtab8560/emtab8560_loxl2_broad_tec_mtec_like_by_age.png`
- `results/figures/external_validation/emtab8560/emtab8560_loxl2_broad_tec_mtec_like_by_age.pdf`
- `results/figures/external_validation/emtab8560/emtab8560_lox_family_detection_by_age.png`
- `results/figures/external_validation/emtab8560/emtab8560_lox_family_detection_by_age.pdf`