# v5.4 Baran-Gale and falsification file check

## Created outputs

Scripts:

- `scripts/mtec1_loxl2_matched_gene_falsification.py`
- `scripts/external_validation_baran_gale_mousethymusageing_export.R`
- `scripts/external_validation_baran_gale_mousethymusageing_lox.py`

Reports:

- `reports/mtec1_loxl2_matched_gene_falsification.md`
- `reports/baran_gale_mousethymusageing_lox_reanalysis.md`
- `reports/v5_4_baran_gale_falsification_update_summary.md`
- `reports/v5_4_baran_gale_falsification_safety_check.md`
- `reports/v5_4_baran_gale_falsification_file_check.md`

Tables:

- `results/tables/mtec1_loxl2_matched_gene_falsification_all_genes.tsv`
- `results/tables/mtec1_loxl2_matched_gene_falsification_matched_genes.tsv`
- `results/tables/mtec1_loxl2_matched_gene_falsification_summary.tsv`
- `results/tables/external_baran_gale_metadata_summary.tsv`
- `results/tables/external_baran_gale_lox_by_sample.tsv`
- `results/tables/external_baran_gale_lox_age_summary.tsv`
- `results/tables/external_baran_gale_lox_models.tsv`

Figures:

- `results/figures/falsification/mtec1_loxl2_matched_gene_detection_delta.png`
- `results/figures/falsification/mtec1_loxl2_matched_gene_detection_delta.pdf`
- `results/figures/falsification/mtec1_loxl2_matched_gene_log2cpm_delta.png`
- `results/figures/falsification/mtec1_loxl2_matched_gene_log2cpm_delta.pdf`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_detection_by_age.png`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_detection_by_age.pdf`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_log2cpm_by_age.png`
- `results/figures/external_validation/baran_gale/baran_gale_loxl2_log2cpm_by_age.pdf`

## Staged-file hygiene

Checked staged paths with `git diff --cached --name-only`.

- No `.h5ad` files staged: yes.
- No `data/raw` files staged: yes.
- No `data/processed` files staged: yes.
- No `data/external` files staged: yes.
- No environment files staged: yes.
- No private outreach files staged: yes.
- No unrelated manuscript drafts staged: yes.
- Text-only staged trace scan passed: yes.

Large data note: the staged files include derived TSV summaries and figure files only. No raw matrices or large external data files are staged.
