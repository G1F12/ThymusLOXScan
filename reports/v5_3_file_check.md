# v5.3 file check

## Required files and outputs

| Check | Status |
|---|---|
| `LICENSE` exists | yes |
| GTEx script exists | yes: `scripts/external_validation_gtex_thymus_lox.py` |
| GTEx output exists | yes: feasibility output created |
| GTEx feasibility report exists | yes: `reports/gtex_thymus_lox_age_feasibility.md` |
| mTEC1 per-sample figure exists | yes: PNG and PDF outputs created |
| mTEC1 per-sample report exists | yes: `reports/mtec1_loxl2_per_sample_raw_figure.md` |
| `wet_lab_validation_plan.md` includes candidate antibody table | yes |
| README includes MIT License section | yes |
| Current version remains v5.2 | yes |

## Data hygiene

No large raw data files are required by this update. The GTEx script does not download large GTEx archives automatically.

The update intentionally does not stage:

- `.h5ad` files
- `data/`
- `data/raw/`
- `data/processed/`
- `data/external/`
- large GTEx raw archives

At final staging, the staged file list contains no `.h5ad` paths and no paths under `data/`, `data/raw/`, `data/processed/`, or `data/external/`.

## Expected staged files for this update

- `LICENSE`
- `README.md`
- `wet_lab_validation_plan.md`
- `scripts/external_validation_gtex_thymus_lox.py`
- `scripts/figures/plot_mtec1_loxl2_per_sample_raw.py`
- `reports/gtex_thymus_lox_age_feasibility.md`
- `results/tables/external_gtex_thymus_lox_feasibility.tsv`
- `results/figures/per_sample/mtec1_loxl2_detection_per_sample_raw.png`
- `results/figures/per_sample/mtec1_loxl2_detection_per_sample_raw.pdf`
- `results/figures/per_sample/mtec1_loxl2_log2cpm_per_sample_raw.png`
- `results/figures/per_sample/mtec1_loxl2_log2cpm_per_sample_raw.pdf`
- `results/figures/per_sample/mtec1_loxl2_per_sample_raw_2panel.png`
- `results/figures/per_sample/mtec1_loxl2_per_sample_raw_2panel.pdf`
- `reports/mtec1_loxl2_per_sample_raw_figure.md`
- `reports/v5_3_license_gtex_antibody_update_summary.md`
- `reports/v5_3_safety_check.md`
- `reports/v5_3_file_check.md`

## Final staged-file confirmation

`git diff --cached --name-only` was checked after staging. The staged files are limited to license/readme/wet-lab planning text, scripts, GTEx feasibility output, mTEC1 per-sample figure outputs, and v5.3 reports.

Additional staged-path checks:

- `.h5ad` staged: no
- `data/` staged: no
- `data/raw/` staged: no
- `data/processed/` staged: no
- `data/external/` staged: no
- large raw GTEx archive staged: no
