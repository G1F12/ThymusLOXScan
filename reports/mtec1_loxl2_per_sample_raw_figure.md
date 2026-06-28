# mTEC1 Loxl2 per-sample raw figure

Each point is one biological sample/animal, not one cell.

The figure shows n=2 young `02mo` and n=2 aged `18mo` samples. Both young samples are higher than both aged samples for Loxl2 detection rate and for Loxl2 log2(CPM+1). This is a descriptive visualization, not strong statistical inference.

## Sample values

| sample | age | mTEC1 cells | detection rate (%) | log2(CPM+1) |
|---|---|---:|---:|---:|
| mo02_CD45neg1_d0 | 02mo | 300 | 8.667 | 3.444 |
| mo02_CD45neg2_d0 | 02mo | 279 | 9.319 | 3.484 |
| mo18_CD45neg1_d0 | 18mo | 288 | 0.347 | 0.604 |
| mo18_CD45neg2_d0 | 18mo | 763 | 0.917 | 1.296 |

## Direction check

- Lowest young detection rate: 8.667%; highest aged detection rate: 0.917%.
- Lowest young log2(CPM+1): 3.444; highest aged log2(CPM+1): 1.296.
- No statistical test is added because the comparison has two biological samples per age group.

## Output files

- `results/figures/per_sample/mtec1_loxl2_detection_per_sample_raw.png`
- `results/figures/per_sample/mtec1_loxl2_detection_per_sample_raw.pdf`
- `results/figures/per_sample/mtec1_loxl2_log2cpm_per_sample_raw.png`
- `results/figures/per_sample/mtec1_loxl2_log2cpm_per_sample_raw.pdf`
- `results/figures/per_sample/mtec1_loxl2_per_sample_raw_2panel.png`
- `results/figures/per_sample/mtec1_loxl2_per_sample_raw_2panel.pdf`