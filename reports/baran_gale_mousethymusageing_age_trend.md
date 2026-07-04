# Baran-Gale Loxl2 age trend audit

## Scope

This trend audit reuses the current Baran-Gale fallback summaries but downgrades inference because true biological mouse/sample IDs were not verified. P-values are therefore not reported.

## Classification

Baran-Gale Loxl2 trend classification: **mostly decreasing descriptive trend in broad TEC/mTEC-like sample-proxy summaries, with sample identity uncertain**.

The aged-lower direction appears in broad TEC, mTEC-like, mTEClo, mTEChi, and cTEC summaries. Post-Aire mTEC is inconsistent and sparse.

## Trend summary

| group | ages | sample proxies per age | oldest-minus-youngest log2CPM+1 | descriptive rho | trend |
|---|---|---|---:|---:|---|
| broad_TEC | 1;4;16;32;52 | 1:5;4:5;16:5;32:5;52:5 | -2.027 | -0.900 | mostly decreasing |
| mTEC_like | 1;4;16;32;52 | 1:5;4:5;16:5;32:5;52:5 | -2.168 | -0.900 | mostly decreasing |
| mTEClo | 1;4;16;32;52 | 1:5;4:5;16:5;32:5;52:4 | -3.536 | -0.900 | mostly decreasing |
| mTEChi | 1;4;16;32;52 | 1:5;4:5;16:5;32:5;52:5 | -1.259 | -0.900 | mostly decreasing |
| post_Aire_mTEC | 1;4;16;32;52 | 1:4;4:5;16:5;32:5;52:4 | 0.000 | -0.335 | inconsistent |
| cTEC | 1;4;16;32;52 | 1:5;4:5;16:5;32:5;52:5 | -1.582 | -0.900 | mostly decreasing |

Descriptive rho and linear coefficients are calculated across age-level means only. They are not biological-sample inference.

## Interpretation

The Baran-Gale age series is useful as provisional context that Loxl2 can be lower at older ages in several TEC/mTEC-like summaries. It is not a substitute for per-animal validation of the GSE240016 mTEC1 result.

## Output

- `results/tables/external_baran_gale_loxl2_age_trend.tsv`
