#!/usr/bin/env python
"""Marker sanity checks and stricter marker-positive LOX summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "annotation_sanity"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "annotation_uncertainty_check.md"

YOUNG_LABEL = "02mo"
AGED_LABEL = "18mo"
CPM_SCALE = 1_000_000
MIN_CELLS_FOR_LOX_SUMMARY = 10

GROUPS = {
    "capsFB": {
        "annotation": "3:capsFB",
        "lineage_annotation": "FB",
        "markers": ["Dpp4", "Smpd3", "Pi16"],
        "lineage_markers": ["Pdgfra"],
        "key_gene": "Lox",
        "expected_direction": "down",
    },
    "intFB": {
        "annotation": "4:intFB",
        "lineage_annotation": "FB",
        "markers": ["Inmt", "Gpx3"],
        "lineage_markers": ["Pdgfra"],
        "key_gene": None,
        "expected_direction": None,
    },
    "medFB": {
        "annotation": "5:medFB",
        "lineage_annotation": "FB",
        "markers": ["Ptn", "Postn"],
        "lineage_markers": ["Pdgfra"],
        "key_gene": None,
        "expected_direction": None,
    },
    "mTEC1": {
        "annotation": "13:mTEC1",
        "lineage_annotation": "TEC",
        "markers": ["Ccl21a", "Itgb4", "Ly6a"],
        "lineage_markers": ["Epcam", "H2-Aa", "Krt8", "Krt5"],
        "key_gene": "Loxl2",
        "expected_direction": "down",
    },
}

KEY_TESTS = [
    {"group": "capsFB", "gene": "Lox", "expected_direction": "down"},
    {"group": "medFB", "gene": "Loxl1", "expected_direction": "up"},
    {"group": "medFB", "gene": "Loxl2", "expected_direction": "down"},
    {"group": "mTEC1", "gene": "Loxl2", "expected_direction": "down"},
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def as_array(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def matrix_gene_vector(matrix, var_names: pd.Index, gene: str) -> np.ndarray:
    idx = var_names.get_loc(gene)
    return as_array(matrix[:, idx])


def available_genes(var_names: pd.Index, genes: list[str]) -> list[str]:
    return [gene for gene in genes if gene in var_names]


def expression_frame(adata, genes: list[str]) -> pd.DataFrame:
    var_names = pd.Index(adata.var_names.astype(str))
    obs = adata.obs[["sample", "stage", "cell_type", "cell_type_subset"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)
    for gene in available_genes(var_names, genes):
        obs[gene] = matrix_gene_vector(adata.X, var_names, gene)
    return obs


def raw_expression_frame(adata, genes: list[str]) -> pd.DataFrame:
    if adata.raw is None:
        raise ValueError("Input AnnData must contain adata.raw.X raw counts.")
    raw_names = pd.Index(adata.raw.var_names.astype(str))
    obs = adata.obs[["sample", "stage", "cell_type", "cell_type_subset"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)
    obs["_row_pos"] = np.arange(adata.n_obs)
    for gene in available_genes(raw_names, genes):
        obs[gene] = matrix_gene_vector(adata.raw.X, raw_names, gene)
    obs["_total_counts"] = as_array(adata.raw.X.sum(axis=1))
    return obs


def add_marker_metrics(df: pd.DataFrame, group: str) -> pd.DataFrame:
    markers = available_genes(pd.Index(df.columns), GROUPS[group]["markers"])
    lineage = available_genes(pd.Index(df.columns), GROUPS[group]["lineage_markers"])
    out = df.copy()
    all_markers = markers + lineage
    for gene in all_markers:
        out[f"{gene}_detected"] = out[gene] > 0
    out["marker_score"] = out[markers].mean(axis=1) if markers else np.nan
    out["n_subtype_markers_detected"] = out[[f"{gene}_detected" for gene in markers]].sum(axis=1) if markers else 0
    out["n_lineage_markers_detected"] = out[[f"{gene}_detected" for gene in lineage]].sum(axis=1) if lineage else 0
    min_subtype = 2 if len(markers) >= 3 else 1
    min_lineage = 1 if lineage else 0
    out["strict_marker_positive"] = (
        (out["n_subtype_markers_detected"] >= min_subtype)
        & (out["n_lineage_markers_detected"] >= min_lineage)
    )
    return out


def marker_summary_for_group(expr_df: pd.DataFrame, group: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = GROUPS[group]
    focus = expr_df.loc[expr_df["cell_type_subset"].eq(cfg["annotation"])].copy()
    focus = add_marker_metrics(focus, group)
    markers = available_genes(pd.Index(focus.columns), cfg["markers"] + cfg["lineage_markers"])

    rows = []
    for sample, sample_df in focus.groupby("sample", sort=True, observed=True):
        stage = sample_df["stage"].iloc[0]
        base = {
            "group": group,
            "sample_id": sample,
            "age_group": stage,
            "n_cells": int(len(sample_df)),
            "marker_score_mean": float(sample_df["marker_score"].mean()),
            "strict_marker_positive_fraction": float(sample_df["strict_marker_positive"].mean()),
        }
        for marker in markers:
            base[f"{marker}_mean"] = float(sample_df[marker].mean())
            base[f"{marker}_detection_rate"] = float((sample_df[marker] > 0).mean())
        rows.append(base)
    sample_summary = pd.DataFrame(rows)

    cell_rows = []
    for marker in markers:
        cell_rows.append(
            focus[["sample", "stage", "cell_type_subset", marker]]
            .rename(columns={"sample": "sample_id", "stage": "age_group", marker: "expression"})
            .assign(group=group, marker=marker)
        )
    long_cells = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    return sample_summary, long_cells


def plot_marker_score(sample_summary: pd.DataFrame, group: str, out_dir: Path) -> Path:
    sns.set_theme(style="whitegrid", context="paper")
    df = sample_summary.loc[sample_summary["group"].eq(group)].copy()
    df["age_group"] = pd.Categorical(df["age_group"], [YOUNG_LABEL, AGED_LABEL], ordered=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    sns.stripplot(data=df, x="age_group", y="marker_score_mean", hue="age_group", size=7, ax=axes[0], legend=False)
    sns.lineplot(data=df, x="age_group", y="marker_score_mean", units="sample_id", estimator=None, color="0.55", ax=axes[0])
    axes[0].set_title(f"{group}: marker score by sample")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Mean normalized marker expression")

    sns.stripplot(data=df, x="age_group", y="strict_marker_positive_fraction", hue="age_group", size=7, ax=axes[1], legend=False)
    sns.lineplot(data=df, x="age_group", y="strict_marker_positive_fraction", units="sample_id", estimator=None, color="0.55", ax=axes[1])
    axes[1].set_title(f"{group}: strict marker-positive fraction")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Fraction")
    axes[1].set_ylim(-0.03, 1.03)
    fig.tight_layout()
    path = out_dir / f"{group}_marker_score_by_sample.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_marker_heatmap(sample_summary: pd.DataFrame, group: str, out_dir: Path) -> Path:
    df = sample_summary.loc[sample_summary["group"].eq(group)].copy()
    markers = [
        col.removesuffix("_detection_rate")
        for col in df.columns
        if col.endswith("_detection_rate")
    ]
    heat = df.set_index("sample_id")[[f"{marker}_detection_rate" for marker in markers]]
    heat.columns = markers
    fig, ax = plt.subplots(figsize=(max(5, len(markers) * 0.75), 3.8))
    sns.heatmap(heat, vmin=0, vmax=1, cmap="viridis", annot=True, fmt=".2f", ax=ax)
    ax.set_title(f"{group}: marker detection rates by sample")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    path = out_dir / f"{group}_marker_detection_heatmap.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_marker_violin(long_cells: pd.DataFrame, group: str, out_dir: Path) -> Path:
    sns.set_theme(style="whitegrid", context="paper")
    df = long_cells.loc[long_cells["group"].eq(group)].copy()
    df["age_group"] = pd.Categorical(df["age_group"], [YOUNG_LABEL, AGED_LABEL], ordered=True)
    n_markers = df["marker"].nunique()
    fig, axes = plt.subplots(n_markers, 1, figsize=(8.4, max(2.3 * n_markers, 3.0)), squeeze=False)
    for ax, marker in zip(axes.ravel(), sorted(df["marker"].unique())):
        mdf = df.loc[df["marker"].eq(marker)]
        sns.violinplot(data=mdf, x="sample_id", y="expression", hue="age_group", cut=0, inner=None, ax=ax)
        ax.set_title(f"{group}: {marker}")
        ax.set_xlabel("")
        ax.set_ylabel("Normalized expression")
        ax.tick_params(axis="x", labelrotation=35)
        ax.legend(title="", loc="upper right", frameon=True)
    fig.tight_layout()
    path = out_dir / f"{group}_marker_expression_by_sample.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return path


def summarize_lox(raw_df: pd.DataFrame, expr_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for test in KEY_TESTS:
        group = test["group"]
        gene = test["gene"]
        cfg = GROUPS[group]
        if gene not in raw_df.columns:
            continue
        annotated_raw = raw_df.loc[raw_df["cell_type_subset"].eq(cfg["annotation"])].copy()
        annotated_expr = expr_df.loc[expr_df["cell_type_subset"].eq(cfg["annotation"])].copy()
        annotated_expr = add_marker_metrics(annotated_expr, group)
        annotated_raw["strict_marker_positive"] = annotated_expr["strict_marker_positive"].to_numpy()

        for subset_name, subset_df in [
            ("annotated_all", annotated_raw),
            ("strict_marker_positive", annotated_raw.loc[annotated_raw["strict_marker_positive"]].copy()),
        ]:
            for sample, sample_df in subset_df.groupby("sample", sort=True, observed=True):
                if len(sample_df) < MIN_CELLS_FOR_LOX_SUMMARY:
                    continue
                values = sample_df[gene].to_numpy()
                total = float(sample_df["_total_counts"].sum())
                gene_sum = float(values.sum())
                n_cells = int(len(values))
                n_detecting = int((values > 0).sum())
                rows.append(
                    {
                        "group": group,
                        "gene": gene,
                        "subset": subset_name,
                        "sample_id": sample,
                        "age_group": sample_df["stage"].iloc[0],
                        "n_cells": n_cells,
                        "raw_counts_sum": gene_sum,
                        "total_counts_sum": total,
                        "cpm": gene_sum / total * CPM_SCALE if total > 0 else np.nan,
                        "log2_cpm_plus1": np.log2(gene_sum / total * CPM_SCALE + 1) if total > 0 else np.nan,
                        "n_detecting_cells": n_detecting,
                        "detection_rate": n_detecting / n_cells if n_cells else np.nan,
                        "mean_nonzero_expression": float(values[values > 0].mean()) if n_detecting else np.nan,
                        "expected_direction": test["expected_direction"],
                    }
                )
    return pd.DataFrame(rows)


def compare_age(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, df in summary.groupby(["group", "gene", "subset"], observed=True):
        group, gene, subset = keys
        means = df.groupby("age_group", observed=True)[["log2_cpm_plus1", "detection_rate", "n_cells"]].mean()
        young = means.loc[YOUNG_LABEL] if YOUNG_LABEL in means.index else pd.Series(dtype=float)
        aged = means.loc[AGED_LABEL] if AGED_LABEL in means.index else pd.Series(dtype=float)
        rows.append(
            {
                "group": group,
                "gene": gene,
                "subset": subset,
                "n_samples_02mo": int(df["age_group"].eq(YOUNG_LABEL).sum()),
                "n_samples_18mo": int(df["age_group"].eq(AGED_LABEL).sum()),
                "mean_log2_cpm_02mo": float(young.get("log2_cpm_plus1", np.nan)),
                "mean_log2_cpm_18mo": float(aged.get("log2_cpm_plus1", np.nan)),
                "delta_log2_cpm_18mo_minus_02mo": float(
                    aged.get("log2_cpm_plus1", np.nan) - young.get("log2_cpm_plus1", np.nan)
                ),
                "mean_detection_02mo": float(young.get("detection_rate", np.nan)),
                "mean_detection_18mo": float(aged.get("detection_rate", np.nan)),
                "delta_detection_18mo_minus_02mo": float(
                    aged.get("detection_rate", np.nan) - young.get("detection_rate", np.nan)
                ),
                "mean_n_cells_02mo": float(young.get("n_cells", np.nan)),
                "mean_n_cells_18mo": float(aged.get("n_cells", np.nan)),
            }
        )
    return pd.DataFrame(rows)


def plot_lox_summaries(lox_summary: pd.DataFrame, out_dir: Path) -> Path:
    sns.set_theme(style="whitegrid", context="paper")
    tests = lox_summary[["group", "gene"]].drop_duplicates().to_dict("records")
    fig, axes = plt.subplots(len(tests), 2, figsize=(8.6, 2.7 * len(tests)), squeeze=False)
    for row_idx, test in enumerate(tests):
        df = lox_summary.loc[
            lox_summary["group"].eq(test["group"]) & lox_summary["gene"].eq(test["gene"])
        ].copy()
        df["age_group"] = pd.Categorical(df["age_group"], [YOUNG_LABEL, AGED_LABEL], ordered=True)
        sns.stripplot(data=df, x="age_group", y="log2_cpm_plus1", hue="subset", dodge=True, ax=axes[row_idx, 0])
        axes[row_idx, 0].set_title(f"{test['group']} {test['gene']}: expression")
        axes[row_idx, 0].set_xlabel("")
        axes[row_idx, 0].set_ylabel("log2(CPM + 1)")
        axes[row_idx, 0].legend(title="", fontsize=7)
        sns.stripplot(data=df, x="age_group", y="detection_rate", hue="subset", dodge=True, ax=axes[row_idx, 1])
        axes[row_idx, 1].set_title(f"{test['group']} {test['gene']}: detection")
        axes[row_idx, 1].set_xlabel("")
        axes[row_idx, 1].set_ylabel("Detection rate")
        axes[row_idx, 1].set_ylim(-0.03, 1.03)
        axes[row_idx, 1].legend(title="", fontsize=7)
    fig.tight_layout()
    path = out_dir / "strict_marker_positive_lox_summary.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return path


def direction_call(delta: float, expected: str | None) -> str:
    if not np.isfinite(delta) or expected is None:
        return "not_evaluable"
    if expected == "down":
        return "expected_direction_retained" if delta < 0 else "expected_direction_not_retained"
    if expected == "up":
        return "expected_direction_retained" if delta > 0 else "expected_direction_not_retained"
    return "not_evaluable"


def write_report(
    marker_summary: pd.DataFrame,
    lox_summary: pd.DataFrame,
    lox_compare: pd.DataFrame,
    figure_paths: list[Path],
    report_path: Path,
) -> None:
    marker_age = []
    for group, df in marker_summary.groupby("group", observed=True):
        means = df.groupby("age_group", observed=True)[["marker_score_mean", "strict_marker_positive_fraction"]].mean()
        marker_age.append(
            {
                "group": group,
                "marker_score_02mo": means.loc[YOUNG_LABEL, "marker_score_mean"] if YOUNG_LABEL in means.index else np.nan,
                "marker_score_18mo": means.loc[AGED_LABEL, "marker_score_mean"] if AGED_LABEL in means.index else np.nan,
                "strict_fraction_02mo": means.loc[YOUNG_LABEL, "strict_marker_positive_fraction"] if YOUNG_LABEL in means.index else np.nan,
                "strict_fraction_18mo": means.loc[AGED_LABEL, "strict_marker_positive_fraction"] if AGED_LABEL in means.index else np.nan,
            }
        )
    marker_age_df = pd.DataFrame(marker_age)

    lines = [
        "# Annotation uncertainty check for key LOX-family results",
        "",
        "## Scope",
        "",
        "This analysis inspects whether cells labeled as capsFB, intFB, medFB, and mTEC1 show expected marker expression and whether key LOX-family summaries retain their direction in stricter marker-positive subsets. It is a sanity check, not proof that the original annotations are correct.",
        "",
        "## Marker sources",
        "",
        "The original Kousa et al. Nature Immunology paper states that major stromal lineages were annotated using canonical markers and maps fibroblast and TEC subsets to published signatures. The paper lists capsFB markers including `Dpp4`, `Smpd3`, and `Pi16`; medFB markers `Ptn` and `Postn`; intFB markers `Inmt` and `Gpx3`; and mTEC1 markers `Ccl21a`, `Itgb4`, and `Ly6a`. It also lists broad lineage markers including FB `Pdgfra`, TEC `Epcam` and `H2-Aa`, EC `Pecam1` and `Cdh5`, MEC `Upk3b` and `Nkain4`, vSMC/PC `Acta2` and `Myl9`, and nmSC `Gfap`, `Ngfr`, and `S100b`.",
        "",
        "## Marker sanity summary",
        "",
        "| group | markers used | mean marker score 02mo | mean marker score 18mo | strict marker-positive fraction 02mo | strict marker-positive fraction 18mo | interpretation |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for _, row in marker_age_df.sort_values("group").iterrows():
        group = row["group"]
        markers = ", ".join(GROUPS[group]["markers"] + GROUPS[group]["lineage_markers"])
        strict_min = min(row["strict_fraction_02mo"], row["strict_fraction_18mo"])
        interp = "marker expression broadly supports label consistency" if strict_min >= 0.50 else "marker support is partial; annotation should be treated cautiously"
        lines.append(
            f"| {group} | {markers} | {row['marker_score_02mo']:.3f} | {row['marker_score_18mo']:.3f} | "
            f"{row['strict_fraction_02mo']:.3f} | {row['strict_fraction_18mo']:.3f} | {interp} |"
        )

    lines.extend(
        [
            "",
            "## Stricter marker-positive LOX summaries",
            "",
            "| group | gene | subset | n 02mo | n 18mo | delta log2(CPM+1) | delta detection | direction check |",
            "|---|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for _, row in lox_compare.sort_values(["group", "gene", "subset"]).iterrows():
        expected = next((test["expected_direction"] for test in KEY_TESTS if test["group"] == row["group"] and test["gene"] == row["gene"]), None)
        call = direction_call(row["delta_log2_cpm_18mo_minus_02mo"], expected)
        lines.append(
            f"| {row['group']} | {row['gene']} | {row['subset']} | {int(row['n_samples_02mo'])} | "
            f"{int(row['n_samples_18mo'])} | {row['delta_log2_cpm_18mo_minus_02mo']:.3f} | "
            f"{row['delta_detection_18mo_minus_02mo']:.3f} | {call} |"
        )

    lines.extend(["", "## Interpretation", ""])
    for test in KEY_TESTS:
        df = lox_compare.loc[lox_compare["group"].eq(test["group"]) & lox_compare["gene"].eq(test["gene"])]
        annotated = df.loc[df["subset"].eq("annotated_all")]
        strict = df.loc[df["subset"].eq("strict_marker_positive")]
        if annotated.empty or strict.empty:
            lines.append(f"- {test['group']} `{test['gene']}`: not evaluable in both annotated and strict subsets.")
            continue
        a_delta = float(annotated.iloc[0]["delta_log2_cpm_18mo_minus_02mo"])
        s_delta = float(strict.iloc[0]["delta_log2_cpm_18mo_minus_02mo"])
        s_n_y = float(strict.iloc[0]["mean_n_cells_02mo"])
        s_n_a = float(strict.iloc[0]["mean_n_cells_18mo"])
        call = direction_call(s_delta, test["expected_direction"])
        if call == "expected_direction_retained":
            statement = "direction is retained in the stricter marker-positive subset"
        else:
            statement = "direction is not retained in the stricter marker-positive subset"
        lines.append(
            f"- {test['group']} `{test['gene']}`: annotated delta={a_delta:.3f}, strict delta={s_delta:.3f}; "
            f"{statement}. Mean strict-subset cells per sample were {s_n_y:.1f} at 02mo and {s_n_a:.1f} at 18mo, so this remains descriptive."
        )

    lines.extend(
        [
            "",
            "## Figures",
            "",
        ]
    )
    for path in figure_paths:
        lines.append(f"- `{path.as_posix()}`")

    lines.extend(
        [
            "",
            "## Bottom line",
            "",
            "The marker checks provide evidence about consistency with expected marker expression, but they do not prove annotation correctness. Key LOX-family directions that are retained in stricter marker-positive subsets are more robust to simple marker-threshold uncertainty. Any result that loses direction or relies on small strict-subset cell counts should be described as annotation-sensitive.",
        ]
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input_h5ad}", flush=True)
    adata = ad.read_h5ad(args.input_h5ad)
    all_marker_genes = sorted({gene for cfg in GROUPS.values() for gene in cfg["markers"] + cfg["lineage_markers"]})
    key_genes = sorted({test["gene"] for test in KEY_TESTS})
    expr_df = expression_frame(adata, all_marker_genes)
    raw_df = raw_expression_frame(adata, all_marker_genes + key_genes)

    marker_summaries = []
    marker_long = []
    figure_paths = []
    for group in GROUPS:
        sample_summary, long_cells = marker_summary_for_group(expr_df, group)
        marker_summaries.append(sample_summary)
        marker_long.append(long_cells)
        figure_paths.append(plot_marker_score(sample_summary, group, args.output_dir))
        figure_paths.append(plot_marker_heatmap(sample_summary, group, args.output_dir))
        if not long_cells.empty:
            figure_paths.append(plot_marker_violin(long_cells, group, args.output_dir))

    marker_summary = pd.concat(marker_summaries, ignore_index=True)
    marker_long_df = pd.concat(marker_long, ignore_index=True)
    lox_summary = summarize_lox(raw_df, expr_df)
    lox_compare = compare_age(lox_summary)
    figure_paths.append(plot_lox_summaries(lox_summary, args.output_dir))

    marker_summary.to_csv(args.output_dir / "marker_summary_by_sample.tsv", sep="\t", index=False)
    marker_long_df.to_csv(args.output_dir / "marker_expression_long.tsv", sep="\t", index=False)
    lox_summary.to_csv(args.output_dir / "strict_marker_positive_lox_by_sample.tsv", sep="\t", index=False)
    lox_compare.to_csv(args.output_dir / "strict_marker_positive_lox_age_summary.tsv", sep="\t", index=False)
    write_report(marker_summary, lox_summary, lox_compare, figure_paths, args.report)

    print(f"Saved figures and tables under: {args.output_dir}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
