# v5.3 license, GTEx, antibody, and per-sample update summary

## Scope

This update adds repository licensing, a guarded GTEx whole-thymus LOX-family age-context script, a transparent four-sample mTEC1 `Loxl2` figure, and more concrete LOXL2 antibody planning details.

No manuscript file was modified. The current public manuscript remains `manuscript/LOX_thymus_aging_public_preprint_v5_2_dropout_protein_feasibility.md`, and the current public release remains `v5.2-dropout-protein-feasibility`.

## Files added or updated

- Added `LICENSE` with the MIT License for code and scripts.
- Updated `README.md` with a license section that separates code, non-raw-data project text/report reuse, and third-party dataset terms.
- Added `scripts/external_validation_gtex_thymus_lox.py`.
- Added GTEx feasibility outputs because no local donor-level GTEx thymus expression and age metadata were found:
  - `reports/gtex_thymus_lox_age_feasibility.md`
  - `results/tables/external_gtex_thymus_lox_feasibility.tsv`
- Added `scripts/figures/plot_mtec1_loxl2_per_sample_raw.py`.
- Added mTEC1 per-sample raw figures:
  - `results/figures/per_sample/mtec1_loxl2_detection_per_sample_raw.png`
  - `results/figures/per_sample/mtec1_loxl2_detection_per_sample_raw.pdf`
  - `results/figures/per_sample/mtec1_loxl2_log2cpm_per_sample_raw.png`
  - `results/figures/per_sample/mtec1_loxl2_log2cpm_per_sample_raw.pdf`
  - `results/figures/per_sample/mtec1_loxl2_per_sample_raw_2panel.png`
  - `results/figures/per_sample/mtec1_loxl2_per_sample_raw_2panel.pdf`
- Added `reports/mtec1_loxl2_per_sample_raw_figure.md`.
- Updated `wet_lab_validation_plan.md` with candidate LOXL2 antibody catalogue details and high-level pilot planning text.

## GTEx status

Classification: age analysis not possible from available public/local files.

The script checked `data/external/GTEx/` and `data/external/gtex/`. No compatible local GTEx gene TPM matrix, sample attributes, or subject phenotype age metadata were present. The feasibility report states that GTEx remains a planned human whole-tissue analysis requiring donor-level expression and age metadata.

## mTEC1 per-sample figure status

The per-sample figure uses `results/tables/mtec1_loxl2_dropout_depth_audit.tsv`. It plots exactly four biological samples: two young `02mo` and two aged `18mo`. Both young samples are higher than both aged samples for detection rate and log2(CPM+1). The report describes this as descriptive visualization, not strong statistical inference.

## Antibody planning status

The wet-lab plan now lists CAB025848/MAB2639, HPA036257, HPA056542, and ab96233 as candidates to evaluate. The table includes vendor/resource, host/clonality, listed applications, listed species reactivity, IHC/IF relevance, caveats, and suggested pilot use. It does not state that any antibody is already suitable for the exact mouse thymus aging question.

## Version and claim posture

This is a repository improvement commit, not a new scientific release. No new release or tag is recommended from this update alone because GTEx did not produce an interpretable donor-level age analysis and the antibody work is planning-only.
