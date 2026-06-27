# Annotation uncertainty check for key LOX-family results

## Scope

This analysis inspects whether cells labeled as capsFB, intFB, medFB, and mTEC1 show expected marker expression and whether key LOX-family summaries retain their direction in stricter marker-positive subsets. It is a sanity check, not proof that the original annotations are correct.

## Marker sources

The original Kousa et al. Nature Immunology paper states that major stromal lineages were annotated using canonical markers and maps fibroblast and TEC subsets to published signatures. The paper lists capsFB markers including `Dpp4`, `Smpd3`, and `Pi16`; medFB markers `Ptn` and `Postn`; intFB markers `Inmt` and `Gpx3`; and mTEC1 markers `Ccl21a`, `Itgb4`, and `Ly6a`. It also lists broad lineage markers including FB `Pdgfra`, TEC `Epcam` and `H2-Aa`, EC `Pecam1` and `Cdh5`, MEC `Upk3b` and `Nkain4`, vSMC/PC `Acta2` and `Myl9`, and nmSC `Gfap`, `Ngfr`, and `S100b`.

## Marker sanity summary

| group | markers used | mean marker score 02mo | mean marker score 18mo | strict marker-positive fraction 02mo | strict marker-positive fraction 18mo | interpretation |
|---|---|---:|---:|---:|---:|---|
| capsFB | Dpp4, Smpd3, Pi16, Pdgfra | 2.449 | 2.597 | 0.609 | 0.407 | marker support is partial; annotation should be treated cautiously |
| intFB | Inmt, Gpx3, Pdgfra | 2.096 | 1.787 | 0.796 | 0.455 | marker support is partial; annotation should be treated cautiously |
| mTEC1 | Ccl21a, Itgb4, Ly6a, Epcam, H2-Aa, Krt8, Krt5 | 2.505 | 2.871 | 0.876 | 0.976 | marker expression broadly supports label consistency |
| medFB | Ptn, Postn, Pdgfra | 2.068 | 2.131 | 0.728 | 0.517 | marker expression broadly supports label consistency |

## Stricter marker-positive LOX summaries

| group | gene | subset | n 02mo | n 18mo | delta log2(CPM+1) | delta detection | direction check |
|---|---|---|---:|---:|---:|---:|---|
| capsFB | Lox | annotated_all | 3 | 3 | -1.494 | -0.223 | expected_direction_retained |
| capsFB | Lox | strict_marker_positive | 3 | 3 | -1.974 | -0.257 | expected_direction_retained |
| mTEC1 | Loxl2 | annotated_all | 2 | 2 | -2.521 | -0.084 | expected_direction_retained |
| mTEC1 | Loxl2 | strict_marker_positive | 2 | 2 | -2.657 | -0.091 | expected_direction_retained |
| medFB | Loxl1 | annotated_all | 3 | 3 | 0.821 | 0.048 | expected_direction_retained |
| medFB | Loxl1 | strict_marker_positive | 3 | 3 | 1.036 | 0.070 | expected_direction_retained |
| medFB | Loxl2 | annotated_all | 3 | 3 | -0.907 | -0.129 | expected_direction_retained |
| medFB | Loxl2 | strict_marker_positive | 3 | 3 | -0.958 | -0.141 | expected_direction_retained |

## Interpretation

- capsFB `Lox`: annotated delta=-1.494, strict delta=-1.974; direction is retained in the stricter marker-positive subset. Mean strict-subset cells per sample were 569.7 at 02mo and 199.7 at 18mo, so this remains descriptive.
- medFB `Loxl1`: annotated delta=0.821, strict delta=1.036; direction is retained in the stricter marker-positive subset. Mean strict-subset cells per sample were 491.3 at 02mo and 581.7 at 18mo, so this remains descriptive.
- medFB `Loxl2`: annotated delta=-0.907, strict delta=-0.958; direction is retained in the stricter marker-positive subset. Mean strict-subset cells per sample were 491.3 at 02mo and 581.7 at 18mo, so this remains descriptive.
- mTEC1 `Loxl2`: annotated delta=-2.521, strict delta=-2.657; direction is retained in the stricter marker-positive subset. Mean strict-subset cells per sample were 253.5 at 02mo and 501.5 at 18mo, so this remains descriptive.

## Figures

- `D:/ThymusLOXScan/results/figures/annotation_sanity/capsFB_marker_score_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/capsFB_marker_detection_heatmap.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/capsFB_marker_expression_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/intFB_marker_score_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/intFB_marker_detection_heatmap.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/intFB_marker_expression_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/medFB_marker_score_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/medFB_marker_detection_heatmap.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/medFB_marker_expression_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/mTEC1_marker_score_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/mTEC1_marker_detection_heatmap.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/mTEC1_marker_expression_by_sample.png`
- `D:/ThymusLOXScan/results/figures/annotation_sanity/strict_marker_positive_lox_summary.png`

## Bottom line

The marker checks provide evidence about consistency with expected marker expression, but they do not prove annotation correctness. Key LOX-family directions that are retained in stricter marker-positive subsets are more robust to simple marker-threshold uncertainty. Any result that loses direction or relies on small strict-subset cell counts should be described as annotation-sensitive.