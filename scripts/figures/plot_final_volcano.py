#!/usr/bin/env python
"""Generate final publication-quality LOX pseudobulk volcano plot."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_TABLE = PROJECT_ROOT / "results" / "tables" / "lox_pseudobulk_complete_results.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "final"

BLUE = "#0072B2"
ORANGE = "#D55E00"
GRAY = "#7A7A7A"


def label_for(row: pd.Series) -> str:
    group = row["cell_type_or_subtype"]
    return f"{row['gene']} {group}"


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(RESULTS_TABLE)
    plot_df = df.loc[
        df["status"].eq("ok") & df["padj"].notna(),
        [
            "gene",
            "annotation_level",
            "cell_type_or_subtype",
            "log2FoldChange",
            "padj",
            "significance_label",
            "direction",
        ],
    ].copy()
    plot_df["minus_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
    plot_df["is_fdr"] = plot_df["padj"] < 0.05
    plot_df["point_class"] = np.where(
        plot_df["is_fdr"] & (plot_df["log2FoldChange"] > 0),
        "FDR < 0.05, higher in aged",
        np.where(
            plot_df["is_fdr"] & (plot_df["log2FoldChange"] < 0),
            "FDR < 0.05, lower in aged",
            "Not FDR-significant",
        ),
    )

    label_keys = {
        ("subtype", "3:capsFB", "Lox"),
        ("subtype", "5:medFB", "Loxl1"),
        ("subtype", "5:medFB", "Loxl2"),
        ("subtype", "13:mTEC1", "Loxl2"),
        ("cell_type", "FB", "Lox"),
        ("cell_type", "FB", "Loxl1"),
        ("cell_type", "FB", "Loxl2"),
    }
    label_df = plot_df.loc[
        plot_df.apply(
            lambda row: (row["annotation_level"], row["cell_type_or_subtype"], row["gene"])
            in label_keys,
            axis=1,
        )
        & plot_df["is_fdr"]
    ].copy()

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    palette = {
        "Not FDR-significant": GRAY,
        "FDR < 0.05, lower in aged": BLUE,
        "FDR < 0.05, higher in aged": ORANGE,
    }
    for point_class, sub in plot_df.groupby("point_class", sort=False):
        ax.scatter(
            sub["log2FoldChange"],
            sub["minus_log10_padj"],
            s=62 if point_class != "Not FDR-significant" else 34,
            c=palette[point_class],
            edgecolor="black" if point_class != "Not FDR-significant" else "none",
            linewidth=0.4,
            alpha=0.9 if point_class != "Not FDR-significant" else 0.45,
            label=point_class,
            zorder=3 if point_class != "Not FDR-significant" else 2,
        )

    ax.axvline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.8)
    ax.axhline(-np.log10(0.05), color="black", linestyle=":", linewidth=1.0, alpha=0.85)
    x_min = min(-6.1, plot_df["log2FoldChange"].min() - 0.3)
    x_max = max(3.1, plot_df["log2FoldChange"].max() + 0.4)
    ax.set_xlim(x_min, x_max)
    ax.text(
        x_min + 0.08,
        -np.log10(0.05) + 0.04,
        "padj = 0.05",
        ha="left",
        va="bottom",
        fontsize=9,
        color="0.25",
    )

    offsets = {
        ("subtype", "3:capsFB", "Lox"): (-52, 12),
        ("subtype", "5:medFB", "Loxl1"): (14, 16),
        ("subtype", "5:medFB", "Loxl2"): (12, 10),
        ("subtype", "13:mTEC1", "Loxl2"): (14, -18),
        ("cell_type", "FB", "Lox"): (-46, 15),
        ("cell_type", "FB", "Loxl1"): (-50, -18),
        ("cell_type", "FB", "Loxl2"): (-50, -20),
    }
    for _, row in label_df.iterrows():
        key = (row["annotation_level"], row["cell_type_or_subtype"], row["gene"])
        xytext = offsets.get(key, (10, 10))
        ax.annotate(
            label_for(row),
            xy=(row["log2FoldChange"], row["minus_log10_padj"]),
            xytext=xytext,
            textcoords="offset points",
            ha="left" if xytext[0] >= 0 else "right",
            va="center",
            fontsize=9,
            arrowprops={"arrowstyle": "-", "color": "0.35", "linewidth": 0.8},
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "0.85", "alpha": 0.92},
        )

    ax.set_xlabel("DESeq2 log2 fold change (aged 18mo vs young 02mo)", fontsize=11)
    ax.set_ylabel("-log10 adjusted p-value", fontsize=11)
    ax.set_title("LOX-family pseudobulk differential expression", fontsize=13, pad=12)
    ax.legend(frameon=True, fontsize=9, loc="upper right", title=None)
    ax.tick_params(axis="both", labelsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()

    png = OUTPUT_DIR / "Fig1_pseudobulk_volcano_final.png"
    pdf = OUTPUT_DIR / "Fig1_pseudobulk_volcano_final.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png}")
    print(f"Saved {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
