#!/usr/bin/env python
"""Plot per-biological-sample pseudobulk LOX expression for key findings."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
DEFAULT_RESULTS_TABLE = PROJECT_ROOT / "results" / "tables" / "lox_pseudobulk_complete_results.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "per_sample"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "per_sample_outlier_check.md"

YOUNG_LABEL = "02mo"
AGED_LABEL = "18mo"
AGE_ORDER = [YOUNG_LABEL, AGED_LABEL]
AGE_PALETTE = {YOUNG_LABEL: "#0072B2", AGED_LABEL: "#D55E00"}
MIN_CELLS_PER_SAMPLE_GROUP = 10
CPM_SCALE = 1_000_000

REQUESTED_FINDINGS = [
    ("subtype", "3:capsFB", "Lox"),
    ("subtype", "5:medFB", "Loxl1"),
    ("subtype", "5:medFB", "Loxl2"),
    ("subtype", "13:mTEC1", "Loxl2"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create per-sample pseudobulk expression plots for LOX findings."
    )
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--results-table", type=Path, default=DEFAULT_RESULTS_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_SAMPLE_GROUP)
    return parser.parse_args()


def selected_findings(results: pd.DataFrame) -> pd.DataFrame:
    requested = pd.DataFrame(
        [
            {
                "annotation_level": level,
                "cell_type_or_subtype": annotation,
                "gene": gene,
                "selection_reason": "requested_key_finding",
            }
            for level, annotation, gene in REQUESTED_FINDINGS
        ]
    )
    significant = results.loc[
        results["significance_label"].eq("FDR<0.05"),
        ["annotation_level", "cell_type_or_subtype", "gene"],
    ].copy()
    significant["selection_reason"] = "other_FDR_significant"

    combined = pd.concat([requested, significant], ignore_index=True)
    combined = combined.drop_duplicates(
        ["annotation_level", "cell_type_or_subtype", "gene"], keep="first"
    )
    return combined.merge(
        results,
        on=["annotation_level", "cell_type_or_subtype", "gene"],
        how="left",
        suffixes=("", "_result"),
    )


def group_column(annotation_level: str) -> str:
    if annotation_level == "cell_type":
        return "cell_type"
    if annotation_level == "subtype":
        return "cell_type_subset"
    raise ValueError(f"Unknown annotation_level: {annotation_level}")


def vector_to_array(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def pseudobulk_values_for_finding(adata, finding: pd.Series, min_cells: int) -> pd.DataFrame:
    if adata.raw is None:
        raise ValueError("Input AnnData must contain adata.raw.X raw counts.")

    annotation_level = str(finding["annotation_level"])
    annotation = str(finding["cell_type_or_subtype"])
    gene = str(finding["gene"])
    group_col = group_column(annotation_level)

    raw_genes = pd.Index(adata.raw.var_names.astype(str))
    if gene not in raw_genes:
        raise KeyError(f"Gene {gene} is absent from adata.raw.var_names.")

    obs = adata.obs[["sample", "stage", group_col]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs[group_col] = obs[group_col].astype(str)
    obs["_row_pos"] = np.arange(adata.n_obs)
    obs = obs.loc[obs[group_col] == annotation].copy()

    raw_x = adata.raw.X
    if not sparse.issparse(raw_x):
        raw_x = np.asarray(raw_x)
    gene_idx = raw_genes.get_loc(gene)

    rows = []
    for sample, sample_obs in obs.groupby("sample", sort=True, observed=True):
        n_cells = len(sample_obs)
        if n_cells < min_cells:
            continue

        stages = sample_obs["stage"].unique()
        if len(stages) != 1:
            raise ValueError(f"Sample {sample} has multiple stages for {annotation}: {stages}")

        row_idx = sample_obs["_row_pos"].to_numpy()
        sample_matrix = raw_x[row_idx, :]
        total_counts = float(np.asarray(sample_matrix.sum()).ravel()[0])
        gene_counts = float(vector_to_array(raw_x[row_idx, gene_idx]).sum())
        cpm = (gene_counts / total_counts * CPM_SCALE) if total_counts > 0 else np.nan
        rows.append(
            {
                "annotation_level": annotation_level,
                "cell_type_or_subtype": annotation,
                "gene": gene,
                "sample": sample,
                "stage": stages[0],
                "n_cells": n_cells,
                "gene_counts": gene_counts,
                "total_counts": total_counts,
                "cpm": cpm,
                "log2_cpm_plus1": np.log2(cpm + 1) if np.isfinite(cpm) else np.nan,
                "log2FoldChange": finding.get("log2FoldChange", np.nan),
                "padj": finding.get("padj", np.nan),
                "direction": finding.get("direction", "not_available"),
                "selection_reason": finding.get("selection_reason", "not_available"),
            }
        )
    return pd.DataFrame(rows)


def panel_title(row: pd.Series) -> str:
    annotation = row["cell_type_or_subtype"]
    gene = row["gene"]
    padj = row["padj"]
    log2fc = row["log2FoldChange"]
    padj_text = "NA" if pd.isna(padj) else f"{padj:.2g}"
    lfc_text = "NA" if pd.isna(log2fc) else f"{log2fc:.2f}"
    return f"{gene} in {annotation}\nDESeq2 log2FC={lfc_text}, padj={padj_text}"


def plot_values(values: pd.DataFrame, findings: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    n_panels = len(findings)
    n_cols = 2
    n_rows = math.ceil(n_panels / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8.4, 3.35 * n_rows), squeeze=False)
    axes_flat = axes.ravel()

    for ax, (_, finding) in zip(axes_flat, findings.iterrows()):
        mask = (
            values["annotation_level"].eq(finding["annotation_level"])
            & values["cell_type_or_subtype"].eq(finding["cell_type_or_subtype"])
            & values["gene"].eq(finding["gene"])
        )
        df = values.loc[mask].copy()
        df["stage"] = pd.Categorical(df["stage"], categories=AGE_ORDER, ordered=True)

        sns.boxplot(
            data=df,
            x="stage",
            y="log2_cpm_plus1",
            order=AGE_ORDER,
            width=0.48,
            showcaps=True,
            showfliers=False,
            boxprops={"facecolor": "white", "edgecolor": "0.35", "linewidth": 1.0},
            whiskerprops={"color": "0.35", "linewidth": 1.0},
            medianprops={"color": "0.15", "linewidth": 1.2},
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x="stage",
            y="log2_cpm_plus1",
            order=AGE_ORDER,
            hue="stage",
            palette=AGE_PALETTE,
            size=7,
            jitter=0.08,
            edgecolor="black",
            linewidth=0.45,
            ax=ax,
            legend=False,
        )

        y_min = float(df["log2_cpm_plus1"].min())
        y_max = float(df["log2_cpm_plus1"].max())
        y_range = max(y_max - y_min, 0.25)
        ax.set_ylim(y_min - 0.12 * y_range, y_max + 0.24 * y_range)

        for x_pos, stage in enumerate(AGE_ORDER):
            stage_df = df.loc[df["stage"].astype(str).eq(stage)]
            offsets = np.linspace(-0.075, 0.075, max(len(stage_df), 1))
            for offset, (_, sample_row) in zip(offsets, stage_df.sort_values("sample").iterrows()):
                label = f"n={int(sample_row['n_cells'])}"
                ax.annotate(
                    label,
                    xy=(x_pos + offset, sample_row["log2_cpm_plus1"]),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7.0,
                    color="0.25",
                )

        ax.set_title(panel_title(finding), fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("Pseudobulk log2(CPM + 1)")
        ax.tick_params(axis="x", labelrotation=0)

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    fig.suptitle("Per-sample pseudobulk LOX-family expression", fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "lox_per_sample_pseudobulk_expression.png"
    pdf_path = output_dir / "lox_per_sample_pseudobulk_expression.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def monotonic_support(df: pd.DataFrame, direction: str) -> tuple[str, str]:
    young = df.loc[df["stage"].eq(YOUNG_LABEL), "log2_cpm_plus1"].dropna().to_numpy()
    aged = df.loc[df["stage"].eq(AGED_LABEL), "log2_cpm_plus1"].dropna().to_numpy()
    if len(young) == 0 or len(aged) == 0:
        return "insufficient samples", "manual review needed"

    young_mean = float(np.mean(young))
    aged_mean = float(np.mean(aged))
    observed = "up_in_aged" if aged_mean > young_mean else "down_in_aged" if aged_mean < young_mean else "no_change"

    if direction == "down_in_aged":
        complete = bool(np.max(aged) < np.min(young))
        supporting_aged = int(np.sum(aged < np.median(young)))
        total_aged = len(aged)
        support = f"{supporting_aged}/{total_aged} aged samples below the young median"
    elif direction == "up_in_aged":
        complete = bool(np.min(aged) > np.max(young))
        supporting_aged = int(np.sum(aged > np.median(young)))
        total_aged = len(aged)
        support = f"{supporting_aged}/{total_aged} aged samples above the young median"
    else:
        complete = False
        support = f"observed mean direction: {observed}"

    if observed != direction:
        statement = "not consistent in normalized pseudobulk means"
    elif complete:
        statement = "consistent across samples with no overlap between age groups"
    elif supporting_aged == total_aged:
        statement = "consistent across aged samples, with some overlap in the full sample ranges"
    elif supporting_aged >= max(1, total_aged - 1):
        statement = "mostly consistent, but one aged sample is close to or overlaps young values"
    else:
        statement = "possible single-sample sensitivity; do not overinterpret"
    return statement, support


def leave_one_out_direction(df: pd.DataFrame) -> str:
    rows = []
    for sample in df["sample"].unique():
        kept = df.loc[~df["sample"].eq(sample)]
        young = kept.loc[kept["stage"].eq(YOUNG_LABEL), "log2_cpm_plus1"].dropna()
        aged = kept.loc[kept["stage"].eq(AGED_LABEL), "log2_cpm_plus1"].dropna()
        if young.empty or aged.empty:
            continue
        diff = float(aged.mean() - young.mean())
        direction = "up_in_aged" if diff > 0 else "down_in_aged" if diff < 0 else "no_change"
        rows.append(f"drop {sample}: {direction} (delta={diff:.2f})")
    return "; ".join(rows) if rows else "not available"


def write_report(values: pd.DataFrame, findings: pd.DataFrame, report_path: Path) -> None:
    lines = [
        "# Per-sample LOX pseudobulk outlier check",
        "",
        "This report summarizes per-biological-sample pseudobulk normalized expression for key LOX-family findings. Values are log2(CPM + 1) from raw counts summed within each biological sample and annotation group. The plots are robustness visualizations and do not replace the DESeq2 pseudobulk inference.",
        "",
        "## Overall interpretation",
        "",
        "The requested key findings are shown with individual biological samples separated by age. Several comparisons have only two or three samples per age group, so the plots should be read as an outlier check rather than independent confirmation. No conclusion below should be interpreted as stronger than the underlying pseudobulk DESeq2 model.",
        "",
        "## Finding-level notes",
        "",
    ]

    for _, finding in findings.iterrows():
        mask = (
            values["annotation_level"].eq(finding["annotation_level"])
            & values["cell_type_or_subtype"].eq(finding["cell_type_or_subtype"])
            & values["gene"].eq(finding["gene"])
        )
        df = values.loc[mask].copy()
        direction = str(finding.get("direction", "not_available"))
        statement, support = monotonic_support(df, direction)
        loo = leave_one_out_direction(df)
        young_n = int(df["stage"].eq(YOUNG_LABEL).sum())
        aged_n = int(df["stage"].eq(AGED_LABEL).sum())

        lines.extend(
            [
                f"### {finding['gene']} in {finding['cell_type_or_subtype']} ({finding['annotation_level']})",
                "",
                f"- DESeq2 direction: {direction}; log2FoldChange={finding.get('log2FoldChange', np.nan):.3g}; padj={finding.get('padj', np.nan):.3g}.",
                f"- Samples plotted: {young_n} young and {aged_n} aged biological samples.",
                f"- Outlier check: {statement}.",
                f"- Sample support: {support}.",
                f"- Leave-one-sample-out visual direction check: {loo}.",
                "- Interpretation: This result does not appear to be driven by one sample." if "consistent" in statement else "- Interpretation: This result may be sensitive to individual samples; interpret cautiously.",
                "",
            ]
        )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    results = pd.read_csv(args.results_table)
    findings = selected_findings(results)

    print(f"Loading {args.input_h5ad}", flush=True)
    adata = sc.read_h5ad(args.input_h5ad)

    value_frames = []
    for _, finding in findings.iterrows():
        print(
            f"Collecting {finding['gene']} in {finding['cell_type_or_subtype']} "
            f"({finding['annotation_level']})",
            flush=True,
        )
        value_frames.append(pseudobulk_values_for_finding(adata, finding, args.min_cells))

    values = pd.concat(value_frames, ignore_index=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    values_path = args.output_dir / "lox_per_sample_pseudobulk_values.csv"
    values.to_csv(values_path, index=False)

    png_path, pdf_path = plot_values(values, findings, args.output_dir)
    write_report(values, findings, args.report)

    print(f"Saved values: {values_path}", flush=True)
    print(f"Saved PNG: {png_path}", flush=True)
    print(f"Saved PDF: {pdf_path}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
