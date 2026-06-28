#!/usr/bin/env python
"""Plot raw per-biological-sample mTEC1 Loxl2 summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "results" / "tables" / "mtec1_loxl2_dropout_depth_audit.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "per_sample"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "mtec1_loxl2_per_sample_raw_figure.md"

AGE_ORDER = ["02mo", "18mo"]
AGE_LABELS = {"02mo": "young 02mo", "18mo": "aged 18mo"}
AGE_PALETTE = {"02mo": "#0072B2", "18mo": "#D55E00"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def load_values(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {
        "sample_id",
        "age_group",
        "n_mTEC1_cells",
        "Loxl2_detection_rate",
        "Loxl2_log2CPM1",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input table is missing required columns: {sorted(missing)}")
    df = df.loc[df["age_group"].isin(AGE_ORDER)].copy()
    if len(df) != 4:
        raise ValueError(f"Expected exactly four mTEC1 biological samples, found {len(df)}.")
    if df.groupby("age_group")["sample_id"].nunique().to_dict() != {"02mo": 2, "18mo": 2}:
        raise ValueError("Expected two young and two aged biological samples.")
    df["age_group"] = pd.Categorical(df["age_group"], categories=AGE_ORDER, ordered=True)
    df["age_label"] = df["age_group"].map(AGE_LABELS)
    df["detection_percent"] = df["Loxl2_detection_rate"] * 100
    df["log2_cpm_plus1"] = df["Loxl2_log2CPM1"]
    return df.sort_values(["age_group", "sample_id"])


def plot_metric(df: pd.DataFrame, metric: str, ylabel: str, title: str, output_dir: Path, stem: str) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(4.8, 3.9))
    sns.pointplot(
        data=df,
        x="age_label",
        y=metric,
        order=[AGE_LABELS[a] for a in AGE_ORDER],
        color="0.55",
        errorbar=None,
        markers="_",
        linestyles="",
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="age_label",
        y=metric,
        order=[AGE_LABELS[a] for a in AGE_ORDER],
        hue="age_group",
        palette=AGE_PALETTE,
        jitter=0.06,
        size=7,
        edgecolor="black",
        linewidth=0.55,
        ax=ax,
        legend=False,
    )
    for x_pos, age in enumerate(AGE_ORDER):
        sub = df.loc[df["age_group"].astype(str).eq(age)].sort_values("sample_id")
        offsets = [-0.055, 0.055]
        for offset, (_, row) in zip(offsets, sub.iterrows()):
            ax.annotate(
                str(row["sample_id"]),
                xy=(x_pos + offset, row[metric]),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                color="0.25",
            )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_combined(df: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
    specs = [
        ("detection_percent", "Detection rate (% Loxl2-positive mTEC1 cells)", "Detection"),
        ("log2_cpm_plus1", "Loxl2 log2(CPM + 1)", "Expression"),
    ]
    for ax, (metric, ylabel, title) in zip(axes, specs):
        sns.pointplot(
            data=df,
            x="age_label",
            y=metric,
            order=[AGE_LABELS[a] for a in AGE_ORDER],
            color="0.55",
            errorbar=None,
            markers="_",
            linestyles="",
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x="age_label",
            y=metric,
            order=[AGE_LABELS[a] for a in AGE_ORDER],
            hue="age_group",
            palette=AGE_PALETTE,
            jitter=0.06,
            size=7,
            edgecolor="black",
            linewidth=0.55,
            ax=ax,
            legend=False,
        )
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / "mtec1_loxl2_per_sample_raw_2panel.png"
    pdf = output_dir / "mtec1_loxl2_per_sample_raw_2panel.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def write_report(df: pd.DataFrame, paths: list[Path], report: Path) -> None:
    young_min_detection = df.loc[df["age_group"].astype(str).eq("02mo"), "detection_percent"].min()
    aged_max_detection = df.loc[df["age_group"].astype(str).eq("18mo"), "detection_percent"].max()
    young_min_expr = df.loc[df["age_group"].astype(str).eq("02mo"), "log2_cpm_plus1"].min()
    aged_max_expr = df.loc[df["age_group"].astype(str).eq("18mo"), "log2_cpm_plus1"].max()
    lines = [
        "# mTEC1 Loxl2 per-sample raw figure",
        "",
        "Each point is one biological sample/animal, not one cell.",
        "",
        "The figure shows n=2 young `02mo` and n=2 aged `18mo` samples. Both young samples are higher than both aged samples for Loxl2 detection rate and for Loxl2 log2(CPM+1). This is a descriptive visualization, not strong statistical inference.",
        "",
        "## Sample values",
        "",
        "| sample | age | mTEC1 cells | detection rate (%) | log2(CPM+1) |",
        "|---|---|---:|---:|---:|",
    ]
    for _, row in df.iterrows():
        lines.append(
            f"| {row['sample_id']} | {row['age_group']} | {int(row['n_mTEC1_cells'])} | "
            f"{row['detection_percent']:.3f} | {row['log2_cpm_plus1']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Direction check",
            "",
            f"- Lowest young detection rate: {young_min_detection:.3f}%; highest aged detection rate: {aged_max_detection:.3f}%.",
            f"- Lowest young log2(CPM+1): {young_min_expr:.3f}; highest aged log2(CPM+1): {aged_max_expr:.3f}.",
            "- No statistical test is added because the comparison has two biological samples per age group.",
            "",
            "## Output files",
            "",
        ]
    )
    lines.extend(f"- `{path.relative_to(PROJECT_ROOT).as_posix()}`" for path in paths)
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    df = load_values(args.input)
    paths = []
    paths.extend(
        plot_metric(
            df,
            "detection_percent",
            "Detection rate (% Loxl2-positive mTEC1 cells)",
            "mTEC1 Loxl2 detection by biological sample",
            args.output_dir,
            "mtec1_loxl2_detection_per_sample_raw",
        )
    )
    paths.extend(
        plot_metric(
            df,
            "log2_cpm_plus1",
            "Loxl2 log2(CPM + 1)",
            "mTEC1 Loxl2 expression by biological sample",
            args.output_dir,
            "mtec1_loxl2_log2cpm_per_sample_raw",
        )
    )
    paths.extend(plot_combined(df, args.output_dir))
    write_report(df, paths, args.report)
    print(f"Saved figures to: {args.output_dir}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
