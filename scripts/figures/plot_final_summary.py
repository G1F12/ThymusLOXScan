#!/usr/bin/env python
"""Generate final four-panel LOX-family manuscript summary figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PSEUDOBULK = PROJECT_ROOT / "results" / "tables" / "lox_pseudobulk_complete_results.csv"
CORRELATIONS = PROJECT_ROOT / "results" / "sc_spearman_correlations.csv"
FB_EFFECTS = PROJECT_ROOT / "results" / "sc_mannwhitney_FB_combined.csv"
MTEC1_EFFECTS = PROJECT_ROOT / "results" / "sc_mannwhitney_mTEC1.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "final"

LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
DISPLAY_GROUPS = ["3:capsFB", "4:intFB", "5:medFB", "13:mTEC1"]
GROUP_LABELS = {
    "3:capsFB": "capsFB",
    "4:intFB": "intFB",
    "5:medFB": "medFB",
    "13:mTEC1": "mTEC1",
}
GENE_DISPLAY = {
    "Lox": "Lox",
    "Loxl1": "Loxl1",
    "Loxl2": "Loxl2",
    "Loxl3": "Loxl3",
    "Loxl4": "Loxl4",
}
BLUE = "#0072B2"
ORANGE = "#D55E00"
GRAY = "#777777"


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="left",
    )


def p_text(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def plot_heatmap(ax, pb: pd.DataFrame) -> None:
    sub = pb.loc[
        pb["annotation_level"].eq("subtype")
        & pb["cell_type_or_subtype"].isin(DISPLAY_GROUPS)
        & pb["gene"].isin(LOX_GENES)
        & pb["status"].eq("ok")
    ].copy()
    heat = (
        sub.pivot(index="cell_type_or_subtype", columns="gene", values="log2FoldChange")
        .reindex(index=DISPLAY_GROUPS, columns=LOX_GENES)
        .rename(index=GROUP_LABELS, columns=GENE_DISPLAY)
    )
    annot = heat.copy().astype(object)
    raw_group_by_label = {v: k for k, v in GROUP_LABELS.items()}
    for idx in annot.index:
        for col in annot.columns:
            raw_group = raw_group_by_label[idx]
            raw_gene = col
            row = sub.loc[(sub["cell_type_or_subtype"] == raw_group) & (sub["gene"] == raw_gene)]
            if row.empty or pd.isna(row.iloc[0]["log2FoldChange"]):
                annot.loc[idx, col] = ""
            else:
                star = "*" if row.iloc[0]["padj"] < 0.05 else ""
                annot.loc[idx, col] = f"{row.iloc[0]['log2FoldChange']:.2f}{star}"

    sns.heatmap(
        heat,
        cmap="vlag",
        center=0,
        vmin=-3.5,
        vmax=3.5,
        annot=annot,
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "DESeq2 log2FC\naged vs young"},
        ax=ax,
    )
    ax.set_xlabel("LOX-family gene", fontsize=10)
    ax.set_ylabel("Subtype", fontsize=10)
    ax.set_title("Subtype-specific pseudobulk effects", fontsize=12)
    ax.tick_params(axis="x", labelrotation=0, labelsize=9)
    ax.tick_params(axis="y", labelrotation=0, labelsize=9)


def plot_significant_bubbles(ax, pb: pd.DataFrame) -> None:
    sig = pb.loc[pb["significance_label"].eq("FDR<0.05")].copy()
    sig["label"] = sig["cell_type_or_subtype"].replace(GROUP_LABELS)
    sig["minus_log10_padj"] = -np.log10(sig["padj"].clip(lower=1e-300))
    sig = sig.sort_values(["annotation_level", "cell_type_or_subtype", "gene"])
    y_order = sig["label"].drop_duplicates().tolist()
    y_map = {label: i for i, label in enumerate(y_order)}
    x_map = {gene: i for i, gene in enumerate(LOX_GENES)}

    sizes = 80 + 90 * sig["minus_log10_padj"]
    scatter = ax.scatter(
        sig["gene"].map(x_map),
        sig["label"].map(y_map),
        s=sizes,
        c=sig["log2FoldChange"],
        cmap="vlag",
        vmin=-3.5,
        vmax=3.5,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_xticks(range(len(LOX_GENES)))
    ax.set_xticklabels(LOX_GENES, fontsize=9)
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels(y_order, fontsize=9)
    ax.set_xlabel("Gene", fontsize=10)
    ax.set_ylabel("Cell type/subtype", fontsize=10)
    ax.set_title("FDR-significant pseudobulk results", fontsize=12)
    ax.grid(True, axis="both", color="0.9", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_xlim(-0.35, len(LOX_GENES) - 0.65)
    ax.set_ylim(-0.35, len(y_order) - 0.65)
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("log2FC aged vs young", fontsize=9)
    for _, row in sig.iterrows():
        ax.text(
            x_map[row["gene"]],
            y_map[row["label"]],
            f"{row['log2FoldChange']:.1f}",
            ha="center",
            va="center",
            fontsize=8,
            color="black",
        )


def plot_correlation_table(ax) -> None:
    corr = pd.read_csv(CORRELATIONS)
    sub = corr.loc[
        corr["gene1"].eq("Loxl2") & corr["gene2"].isin(["Col1a1", "Vim", "Snai1"])
    ].copy()
    sub["gene2"] = pd.Categorical(sub["gene2"], categories=["Col1a1", "Vim", "Snai1"], ordered=True)
    sub = sub.sort_values("gene2")
    ax.axis("off")
    rows = [
        [row["gene2"], f"{row['spearman_rho']:.3f}", p_text(row["pvalue"]), f"{int(row['n_cells']):,}"]
        for _, row in sub.iterrows()
    ]
    table = ax.table(
        cellText=rows,
        colLabels=["Marker", "Spearman rho", "p-value", "cells"],
        loc="center",
        cellLoc="center",
        colLoc="center",
        colWidths=[0.22, 0.28, 0.27, 0.18],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.05, 1.55)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("0.75")
        if row == 0:
            cell.set_facecolor("#EAEAEA")
            cell.set_text_props(weight="bold")
        else:
            cell.set_facecolor("white")
    ax.set_title("Fibroblast Loxl2 co-expression", fontsize=12, pad=10)
    ax.text(
        0.5,
        0.08,
        "Per-cell Spearman correlations; descriptive only",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=9,
        color="0.3",
    )


def plot_effect_sizes(ax) -> None:
    fb = pd.read_csv(FB_EFFECTS)
    fb["analysis"] = "Pooled fibroblasts"
    mtec = pd.read_csv(MTEC1_EFFECTS)
    mtec["analysis"] = "mTEC1"
    effects = pd.concat([fb, mtec], ignore_index=True)
    effects["gene"] = pd.Categorical(effects["gene"], categories=LOX_GENES, ordered=True)
    effects = effects.sort_values(["analysis", "gene"])

    y_positions = np.arange(len(LOX_GENES))
    offsets = {"Pooled fibroblasts": -0.16, "mTEC1": 0.16}
    colors = {"Pooled fibroblasts": BLUE, "mTEC1": ORANGE}
    for analysis, sub in effects.groupby("analysis", sort=False):
        ax.scatter(
            sub["rank_biserial"],
            y_positions + offsets[analysis],
            s=64,
            color=colors[analysis],
            edgecolor="black",
            linewidth=0.45,
            label=analysis,
            zorder=3,
        )
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(LOX_GENES, fontsize=9)
    ax.set_xlabel("Rank-biserial effect size\n(aged vs young cells)", fontsize=10)
    ax.set_title("Descriptive single-cell effect sizes", fontsize=12)
    ax.legend(frameon=True, fontsize=9, loc="lower right")
    ax.grid(True, axis="x", color="0.9")
    ax.set_xlim(-0.26, 0.08)
    sns.despine(ax=ax)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pb = pd.read_csv(PSEUDOBULK)

    sns.set_theme(style="white", context="paper")
    fig = plt.figure(figsize=(11.2, 8.4))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.36, hspace=0.42)

    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]
    plot_heatmap(axes[0], pb)
    plot_significant_bubbles(axes[1], pb)
    plot_correlation_table(axes[2])
    plot_effect_sizes(axes[3])
    for ax, label in zip(axes, ["A", "B", "C", "D"]):
        panel_label(ax, label)

    fig.suptitle("LOX-family expression changes in aging thymic stroma", fontsize=15, y=0.99)
    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.08, top=0.91, wspace=0.40, hspace=0.44)

    png = OUTPUT_DIR / "Fig2_summary_4panel_final.png"
    pdf = OUTPUT_DIR / "Fig2_summary_4panel_final.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png}")
    print(f"Saved {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
