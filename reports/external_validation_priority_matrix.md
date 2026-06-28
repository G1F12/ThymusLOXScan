# External Validation Priority Matrix

Purpose: rank public thymus datasets by practical value for cautious LOX-family aging validation. Scores reflect feasibility for transcript-level directional checks only.

## Priority By Validation Target

| priority | target claim | best dataset | score | validation level | next action | reason |
|---:|---|---|---:|---|---|---|
| 1 | Broad fibroblast `Lox` aged-lower | GSE223049 | 5 | Independent mouse broad sorted bulk directional context | Already processed; include as broad context only | Small processed count matrix, young-vs-old mice, sorted broad thymic fibroblast samples. |
| 2 | Broad fibroblast `Loxl2` aged-lower | GSE223049 | 5 | Independent mouse broad sorted bulk directional context | Already processed; include in matrix and figures | Strong broad fibroblast aged-lower direction, but no fibroblast subtype resolution. |
| 3 | Broad epithelial `Loxl2` aged-lower | GSE223049 and E-MTAB-8560 | 5 / 4 | Broad epithelial and TEC age-series context | Already processed | GSE223049 tests broad epithelial bulk; E-MTAB-8560 adds TEC age-series context. |
| 4 | mTEC-like `Loxl2` age trend | E-MTAB-8560 | 4 | Approximate mTEC-like directional context | Already processed; never label exact mTEC1 replication | cTEC/mTEClo/mTEChi and inferred groups are useful but not exact mTEC1. |
| 5 | Human fibroblast-like/epithelial LOX-family age association | GSE231906 | 3 | Candidate future human donor-level expression analysis | Keep metadata-only until expression parsing is deliberately scheduled | Broad age range and useful labels, but raw archive is large and donor/source matching matters. |
| 6 | Human stromal/TEC baseline context | CELLxGENE human thymus atlas / GSE147520-E-MTAB-8581 | 2 | Developmental/pediatric context only | Do not use for aging claims | Useful labels and expression assets, but not old-age validation. |
| 7 | Mouse lifespan spatial context | SCP3197 / Zenodo-linked resource | 2 | Possible spatial context | Manual portal inspection only if spatial validation is prioritized | File inventory and cell-group extraction are not ready for immediate LOX analysis. |
| 8 | TEC marker/context only | E-MTAB-8737 or GSE114651 | 1 / 0 | Not suitable for current validation | Do not prioritize | No clean aged-vs-young validation structure. |

## Claim-By-Dataset Feasibility

| claim | immediately usable | moderate work | possible but heavy | not suitable / cannot test | current conclusion |
|---|---|---|---|---|---|
| `broad_FB_Lox_aged_lower` | GSE223049 | None | GSE231906 after guarded human matrix parsing | E-MTAB-8560 lacks fibroblasts | Broad mouse direction can be described from GSE223049; human remains future-only. |
| `broad_FB_Loxl2_aged_lower` | GSE223049 | None | GSE231906 after guarded human matrix parsing | E-MTAB-8560 lacks fibroblasts | Broad mouse direction can be described from GSE223049; not subtype-resolved. |
| `capsFB_Lox_aged_lower` | None | None | Possibly human fibroblast-like approximation after GSE231906 parsing, not exact | GSE223049 broad bulk; E-MTAB-8560 TEC-only | No exact external validation found. |
| `medFB_Loxl1_aged_higher` | None | None | Possibly human fibroblast-like approximation after GSE231906 parsing, not exact | GSE223049 broad fibroblast is not subtype-specific and trends weakly lower | Internally observed and externally unresolved. |
| `medFB_Loxl2_aged_lower` | None for exact medFB | None | GSE231906 possible approximation after parsing | GSE223049 broad fibroblast cannot validate medFB | Broad direction exists, exact medFB validation absent. |
| `epithelial_Loxl2_aged_lower` | GSE223049 | E-MTAB-8560 | GSE231906 after expression parsing | Developmental atlases only context | Broad external directional context is available, still transcript-level only. |
| `mTEC1_Loxl2_aged_lower` | None exact | E-MTAB-8560 approximate mTEC-like context | GSE231906 possible mTEC-like human approximation after parsing | GSE223049 broad epithelial cannot validate mTEC1 | Directionally contextualized, not exact mTEC1 replication. |

## Bottom Line

The release should present GSE223049 and E-MTAB-8560 as descriptive external direction checks, not definitive validation. GSE231906 should remain metadata-only in v5.1 unless a future deliberate expression-analysis task is performed. Exact subtype-resolved external validation of capsFB, medFB, and mTEC1 remains unavailable.
