# Baran-Gale / MouseThymusAgeing LOX-family reanalysis

## Scope

This reanalysis inspects LOX-family transcript abundance in public Baran-Gale / MouseThymusAgeing mouse thymic epithelial ageing data. It is framed as independent TEC/mTEC-like directional context, not exact subtype replication of GSE240016 mTEC1.

Input source used by this run: cached MouseThymusAgeing ExperimentHub SMART-seq2 RDS files.

## Metadata inspection

- dataset / source: MouseThymusAgeing ExperimentHub SMART-seq2 cached RDS
- dataset / technology: SMART-seq2
- dataset / n_cells: 2327
- metadata / age_weeks: 1;4;16;32;52
- metadata / sample_proxy: SortDay; no clean donor-level replicate field was available in parsed colData
- metadata / sort_type_labels: cTEC;gmTEC;mTEChi;mTEClo
- metadata / subtype_labels: Intertypical TEC;Mature cTEC;Mature mTEC;Perinatal cTEC;Post-Aire mTEC;Proliferating TEC;Tuft-like mTEC;nTEC;sTEC
- loading / r_export_status: not used; FileNotFoundError: R export not found: D:\ThymusLOXScan\data\external\BaranGale\baran_gale_lox_cell_values_export.tsv

Candidate groups tested: broad TEC, mTEC-like, mTEClo, mTEChi, post-Aire mTEC when labeled, and cTEC.
Sample-level summaries use SortDay as a sample proxy because parsed metadata did not expose a clean donor-level replicate field.

## Loxl2 direction

| group | sample proxies | age points | oldest-minus-youngest log2CPM delta | Spearman rho | direction |
|---|---:|---:|---:|---:|---|
| broad_TEC | 5 | 5 | -2.027 | -0.612 | aged-lower |
| mTEC_like | 5 | 5 | -2.168 | -0.565 | aged-lower |
| mTEClo | 5 | 5 | -3.536 | -0.762 | aged-lower |
| mTEChi | 5 | 5 | -1.259 | -0.184 | aged-lower |
| post_Aire_mTEC | 5 | 5 | 0.000 | -0.335 | flat |
| cTEC | 5 | 5 | -1.582 | -0.582 | aged-lower |

## Classification

Classification: **completed: Loxl2 aged-lower in TEC/mTEC-like group**.
Key Loxl2 direction in TEC/mTEC-like groups: **aged-lower**.

Baran-Gale / MouseThymusAgeing provides independent TEC/mTEC-like directional support for aged-lower Loxl2 at the transcript level. This is not exact GSE240016 mTEC1 validation.

## Output files

- `results/tables/external_baran_gale_metadata_summary.tsv`
- `results/tables/external_baran_gale_lox_by_sample.tsv`
- `results/tables/external_baran_gale_lox_age_summary.tsv`
- `results/tables/external_baran_gale_lox_models.tsv`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_detection_by_age.png`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_detection_by_age.pdf`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_log2cpm_by_age.png`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_log2cpm_by_age.pdf`