#!/usr/bin/env python
"""Plot cross-dataset LOX-family direction consistency summaries."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MATRIX = PROJECT_ROOT / "results" / "tables" / "cross_dataset_lox_validation_matrix.tsv"
DEFAULT_SUMMARY = PROJECT_ROOT / "results" / "tables" / "cross_dataset_lox_consistency_summary.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "external_validation" / "cross_dataset"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "cross_dataset_figure_notes.md"

STATUS_ORDER = [
    "internally observed",
    "directionally consistent",
    "approximate external context",
    "not testable",
    "opposite direction",
    "metadata only",
]
STATUS_TO_VALUE = {
    "internally observed": 0.5,
    "directionally consistent": 1,
    "approximate external context": 0.75,
    "metadata only": 0.1,
    "opposite direction": -1,
    "not testable": 0,
}
STATUS_COLORS = {
    "internally observed": "#8E8E8E",
    "directionally consistent": "#2C7FB8",
    "approximate external context": "#7BCCC4",
    "metadata only": "#F0E442",
    "opposite direction": "#C43C39",
    "not testable": "#D9D9D9",
}
DATASET_ORDER = ["GSE240016_internal", "GSE223049", "E-MTAB-8560", "GSE231906"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def status_for_row(row: pd.Series) -> str:
    if row["dataset"] == "GSE240016_internal" and row["confidence_level"] == "internal_primary":
        return "internally observed"
    if row["dataset"] == "GSE231906":
        return "metadata only"
    if row["confidence_level"] in {"external_broad_context_not_subtype", "external_TEC_age_series_descriptive"}:
        if row["consistency_score"] == -1:
            return "opposite direction"
        if row["consistency_score"] == 1:
            return "approximate external context"
        return "not testable"
    if row["consistency_score"] == 1:
        return "directionally consistent"
    if row["consistency_score"] == -1:
        return "opposite direction"
    return "not testable"


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_heatmap(matrix: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    df = matrix.copy()
    df["status"] = df.apply(status_for_row, axis=1)
    df["row_label"] = df["dataset"] + " | " + df["cell_group"].astype(str)
    df["column_label"] = df["claim_mapping"].map(
        lambda value: textwrap.fill(str(value).replace("_", " "), width=22)
    )
    df["status_value"] = df["status"].map(STATUS_TO_VALUE)

    row_order = []
    for dataset in DATASET_ORDER:
        row_order.extend(df.loc[df["dataset"].eq(dataset), "row_label"].drop_duplicates().tolist())
    col_order = df["column_label"].drop_duplicates().tolist()

    heat = (
        df.pivot_table(
            index="row_label",
            columns="column_label",
            values="status_value",
            aggfunc="max",
            fill_value=0,
        )
        .reindex(index=row_order, columns=col_order)
        .fillna(0)
    )

    cmap = ListedColormap([STATUS_COLORS[label] for label in STATUS_ORDER])
    norm = BoundaryNorm(np.arange(-0.5, len(STATUS_ORDER) + 0.5, 1), cmap.N)
    value_to_code = {STATUS_TO_VALUE[label]: idx for idx, label in enumerate(STATUS_ORDER)}
    codes = heat.map(lambda value: value_to_code.get(value, 1))

    fig, ax = plt.subplots(figsize=(11.5, max(5.0, 0.34 * len(heat) + 1.4)))
    sns.heatmap(
        codes,
        cmap=cmap,
        norm=norm,
        cbar=False,
        linewidths=0.6,
        linecolor="white",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Cross-dataset LOX-family direction consistency", fontsize=13)
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)

    handles = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=STATUS_COLORS[label], markersize=9)
        for label in STATUS_ORDER
    ]
    ax.legend(handles, STATUS_ORDER, title="", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return save_figure(fig, output_dir, "cross_dataset_lox_direction_heatmap")


def plot_effects(matrix: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    df = matrix.loc[
        matrix["effect_metric"].isin(
            [
                "delta_log2_CPM_plus1_aged_minus_young",
                "mean_delta_oldest_minus_youngest_log2_norm_plus1",
            ]
        )
        & matrix["effect_value"].notna()
    ].copy()
    df["status"] = df.apply(status_for_row, axis=1)
    df["label"] = df["dataset"] + " | " + df["cell_group"].astype(str)
    df["claim_gene"] = df["claim_mapping"] + "\n" + df["gene"].astype(str)

    fig, ax = plt.subplots(figsize=(9.6, max(4.2, 0.34 * len(df) + 1.2)))
    sns.scatterplot(
        data=df,
        x="effect_value",
        y="label",
        hue="status",
        style="gene",
        palette=STATUS_COLORS,
        hue_order=STATUS_ORDER,
        s=95,
        edgecolor="black",
        linewidth=0.4,
        ax=ax,
    )
    ax.axvline(0, color="0.25", linestyle="--", linewidth=0.9)
    ax.set_xlabel("Comparable descriptive delta, high age minus low age")
    ax.set_ylabel("")
    ax.set_title("External descriptive effect sizes where comparable", fontsize=13)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return save_figure(fig, output_dir, "cross_dataset_lox_external_effect_points")


def summary_status(row: pd.Series) -> str:
    if row["n_external_opposite"] > 0:
        return "opposite direction"
    if str(row["conclusion"]).startswith("externally_directionally_supported"):
        if row["claim_mapping"] == "mTEC1_Loxl2_aged_lower":
            return "approximate external context"
        return "directionally consistent"
    if row["internal_direction"] in {"aged-lower", "aged-higher"}:
        return "internally observed"
    return "not testable"


def plot_claim_summary(summary: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    df = summary.copy()
    df["status"] = df.apply(summary_status, axis=1)
    df["claim_label"] = df["claim_mapping"].str.replace("_", " ", regex=False)
    df["claim_label"] = [textwrap.fill(label, width=28) for label in df["claim_label"]]
    df["plot_value"] = df["status"].map(STATUS_TO_VALUE)

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    sns.barplot(
        data=df,
        x="plot_value",
        y="claim_label",
        hue="status",
        dodge=False,
        palette=STATUS_COLORS,
        hue_order=STATUS_ORDER,
        ax=ax,
    )
    for container in ax.containers:
        ax.bar_label(container, labels=[""] * len(container), padding=2)
    ax.set_xlim(-1.1, 1.35)
    ax.set_xticks([-1, 0, 0.1, 0.5, 0.75, 1])
    ax.set_xticklabels(["opposite", "not testable", "metadata only", "internally observed", "approximate", "consistent"], rotation=0)
    ax.axvline(0, color="0.4", linewidth=0.8)
    ax.set_xlabel("Cross-dataset direction category")
    ax.set_ylabel("")
    ax.set_title("Claim-level cross-dataset consistency", fontsize=13)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return save_figure(fig, output_dir, "cross_dataset_lox_claim_consistency")


def write_report(report: Path, figure_paths: list[Path], matrix: pd.DataFrame, summary: pd.DataFrame) -> None:
    report.parent.mkdir(parents=True, exist_ok=True)
    counts = summary.apply(summary_status, axis=1).value_counts().to_dict()
    lines = [
        "# Cross-dataset LOX figure notes",
        "",
        "## Scope",
        "",
        "These figures summarize direction consistency only. They do not label any comparison as validated unless subtype-resolved independent replication is present; none of the current external rows meet that stronger standard for capsFB, medFB, or exact mTEC1 specificity.",
        "",
        "## Generated figures",
        "",
    ]
    lines.extend(f"- `{path.relative_to(PROJECT_ROOT).as_posix()}`" for path in figure_paths)
    lines.extend(
        [
            "",
            "## Label definitions",
            "",
            "- internally observed: the direction is observed in GSE240016 but exact external testing is unavailable.",
            "- directionally consistent: an external row matches the expected direction at an appropriate broad resolution.",
            "- approximate external context: an external row has a matching direction but not exact subtype resolution.",
            "- opposite direction: an external row has the opposite direction at an appropriate resolution.",
            "- not testable: the dataset lacks the needed compartment, expression matrix, or resolution.",
            "- metadata only: candidate groups are present, but no expression analysis was performed.",
            "",
            "## Claim status counts",
            "",
        ]
    )
    for label in STATUS_ORDER:
        lines.append(f"- {label}: {int(counts.get(label, 0))}")
    lines.extend(
        [
            "",
            "## Cautions",
            "",
            "GSE223049 broad fibroblast and epithelial bulk rows are not subtype-specific validation. E-MTAB-8560 provides TEC and mTEC-like directional support, not exact mTEC1 replication. GSE231906 remains metadata-only in this repository, so its candidate human groups are not expression validation.",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    matrix = pd.read_csv(args.matrix, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")
    sns.set_theme(style="whitegrid", context="paper")

    figure_paths = []
    for paths in [
        plot_heatmap(matrix, args.output_dir),
        plot_effects(matrix, args.output_dir),
        plot_claim_summary(summary, args.output_dir),
    ]:
        figure_paths.extend(paths)
    write_report(args.report, figure_paths, matrix, summary)
    print(f"Saved figures: {args.output_dir}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
