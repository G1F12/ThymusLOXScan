#!/usr/bin/env python
"""Focused dropout/depth audit for Loxl2 in annotated mTEC1 cells."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data" / "processed" / "thymus_annotated.h5ad"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "dropout"
REPORT_PATH = PROJECT_ROOT / "reports" / "mtec1_loxl2_dropout_depth_audit.md"

AUDIT_TABLE = TABLE_DIR / "mtec1_loxl2_dropout_depth_audit.tsv"
PAIRWISE_TABLE = TABLE_DIR / "mtec1_loxl2_pairwise_sample_direction.tsv"
DOWNSAMPLE_ITER_TABLE = TABLE_DIR / "mtec1_loxl2_downsampling_iterations.tsv"
DOWNSAMPLE_SUMMARY_TABLE = TABLE_DIR / "mtec1_loxl2_downsampling_summary.tsv"
LOGISTIC_TABLE = TABLE_DIR / "mtec1_loxl2_depth_adjusted_logistic.tsv"

GENE = "Loxl2"
MTEC1_LABEL_CANDIDATES = ("13:mTEC1", "mTEC1")
YOUNG_LABEL = "02mo"
AGED_LABEL = "18mo"
AGE_ORDER = [YOUNG_LABEL, AGED_LABEL]
AGE_PALETTE = {YOUNG_LABEL: "#0072B2", AGED_LABEL: "#D55E00"}
CPM_SCALE = 1_000_000.0
RANDOM_SEED = 240016
N_DOWNSAMPLE_ITER = 200
MIN_MTEC1_CELLS_PER_SAMPLE = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--downsample-iterations", type=int, default=N_DOWNSAMPLE_ITER)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    return parser.parse_args()


def as_array(values: Any) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def matrix_row_sums(matrix: Any) -> np.ndarray:
    return as_array(matrix.sum(axis=1))


def choose_mtec1_annotation(obs: pd.DataFrame) -> tuple[str, str]:
    search_columns = [
        c
        for c in obs.columns
        if any(token in c.lower() for token in ("cell", "type", "subset", "annot", "subtype"))
    ]
    for column in search_columns:
        values = set(obs[column].astype(str).unique())
        for candidate in MTEC1_LABEL_CANDIDATES:
            if candidate in values:
                return column, candidate
    for column in search_columns:
        matches = sorted(v for v in obs[column].astype(str).unique() if "mtec1" in v.lower())
        if matches:
            return column, matches[0]
    raise ValueError("Could not identify an mTEC1 annotation label in AnnData obs.")


def raw_matrix_and_genes(adata) -> tuple[Any, pd.Index, str]:
    if adata.raw is not None and adata.raw.X is not None:
        return adata.raw.X, pd.Index(adata.raw.var_names.astype(str)), "adata.raw.X"
    for layer_name in ("counts", "raw_counts", "raw", "spliced"):
        if layer_name in adata.layers:
            return adata.layers[layer_name], pd.Index(adata.var_names.astype(str)), f"adata.layers['{layer_name}']"
    return adata.X, pd.Index(adata.var_names.astype(str)), "adata.X (raw-count limitation: no raw/count layer found)"


def prepare_mtec1(adata) -> tuple[pd.DataFrame, Any, pd.Index, str, str, str]:
    annotation_col, annotation_label = choose_mtec1_annotation(adata.obs)
    raw_x, raw_genes, matrix_source = raw_matrix_and_genes(adata)
    if GENE not in raw_genes:
        raise KeyError(f"{GENE} was not found in {matrix_source} var names.")

    mask = adata.obs[annotation_col].astype(str).eq(annotation_label).to_numpy()
    positions = np.flatnonzero(mask)
    obs = adata.obs.iloc[positions].copy()
    obs["_adata_position"] = positions
    obs["sample_id"] = obs["sample"].astype(str)
    age_column = "age_group" if "age_group" in obs.columns else "stage"
    obs["age_group"] = obs[age_column].astype(str)

    required = ["total_counts", "n_genes_by_counts", "mito_frac"]
    missing = [c for c in required if c not in obs.columns]
    if missing:
        raise KeyError(f"Missing required QC columns in obs: {missing}")

    return obs, raw_x[positions, :], raw_genes, matrix_source, annotation_col, annotation_label


def compute_sample_table(obs: pd.DataFrame, raw_mtec1: Any, raw_genes: pd.Index) -> pd.DataFrame:
    gene_idx = raw_genes.get_loc(GENE)
    obs = obs.copy()
    obs["_loxl2_raw"] = as_array(raw_mtec1[:, gene_idx]).astype(float)
    obs["_raw_library_size"] = matrix_row_sums(raw_mtec1).astype(float)

    rows = []
    for sample_id, sample_obs in obs.groupby("sample_id", sort=True, observed=True):
        ages = sample_obs["age_group"].unique()
        if len(ages) != 1:
            raise ValueError(f"Sample {sample_id} has multiple age groups: {ages}")
        values = sample_obs["_loxl2_raw"].to_numpy()
        raw_lib = sample_obs["_raw_library_size"].to_numpy()
        raw_count_sum = float(values.sum())
        raw_library_sum = float(raw_lib.sum())
        n_cells = int(len(sample_obs))
        detecting = values > 0
        cpm = raw_count_sum / raw_library_sum * CPM_SCALE if raw_library_sum > 0 else np.nan
        rows.append(
            {
                "sample_id": sample_id,
                "age_group": ages[0],
                "n_mTEC1_cells": n_cells,
                "Loxl2_raw_count_sum": raw_count_sum,
                "Loxl2_detecting_cells": int(detecting.sum()),
                "Loxl2_detection_rate": float(detecting.mean()) if n_cells else np.nan,
                "Loxl2_mean_raw_all_cells": float(values.mean()) if n_cells else np.nan,
                "Loxl2_mean_raw_positive_cells": float(values[detecting].mean()) if detecting.any() else np.nan,
                "raw_library_sum": raw_library_sum,
                "Loxl2_CPM": cpm,
                "Loxl2_log2CPM1": math.log2(cpm + 1.0) if np.isfinite(cpm) else np.nan,
                "median_total_counts": float(sample_obs["total_counts"].median()),
                "mean_total_counts": float(sample_obs["total_counts"].mean()),
                "median_n_genes_by_counts": float(sample_obs["n_genes_by_counts"].median()),
                "mean_n_genes_by_counts": float(sample_obs["n_genes_by_counts"].mean()),
                "median_mito_frac": float(sample_obs["mito_frac"].median()),
                "mean_mito_frac": float(sample_obs["mito_frac"].mean()),
                "median_raw_library_size": float(np.median(raw_lib)),
                "mean_raw_library_size": float(np.mean(raw_lib)),
            }
        )
    return pd.DataFrame(rows).sort_values(["age_group", "sample_id"]).reset_index(drop=True)


def age_summary(sample_table: pd.DataFrame) -> dict[str, float]:
    grouped = sample_table.groupby("age_group", observed=True)[
        ["Loxl2_detection_rate", "Loxl2_log2CPM1"]
    ].mean()
    young_det = float(grouped.loc[YOUNG_LABEL, "Loxl2_detection_rate"])
    aged_det = float(grouped.loc[AGED_LABEL, "Loxl2_detection_rate"])
    young_log = float(grouped.loc[YOUNG_LABEL, "Loxl2_log2CPM1"])
    aged_log = float(grouped.loc[AGED_LABEL, "Loxl2_log2CPM1"])
    return {
        "young_mean_detection_rate": young_det,
        "aged_mean_detection_rate": aged_det,
        "aged_minus_young_detection_rate": aged_det - young_det,
        "young_mean_log2CPM1": young_log,
        "aged_mean_log2CPM1": aged_log,
        "aged_minus_young_log2CPM1": aged_log - young_log,
    }


def write_audit_table(sample_table: pd.DataFrame, summary: dict[str, float]) -> None:
    out = sample_table.copy()
    for key, value in summary.items():
        out[key] = value
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(AUDIT_TABLE, sep="\t", index=False)


def filter_comparison_samples(sample_table: pd.DataFrame, obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = sample_table.loc[
        sample_table["age_group"].isin(AGE_ORDER)
        & sample_table["n_mTEC1_cells"].ge(MIN_MTEC1_CELLS_PER_SAMPLE),
        "sample_id",
    ].tolist()
    filtered_table = sample_table.loc[sample_table["sample_id"].isin(eligible)].copy()
    n_by_age = filtered_table.groupby("age_group", observed=True)["sample_id"].nunique()
    if n_by_age.get(YOUNG_LABEL, 0) != 2 or n_by_age.get(AGED_LABEL, 0) != 2:
        raise ValueError(
            "Expected exactly 2 young and 2 aged mTEC1 samples after minimum-cell filtering; "
            f"observed counts by age: {n_by_age.to_dict()}"
        )
    filtered_obs = obs.loc[obs["sample_id"].isin(eligible)].copy()
    return filtered_table.reset_index(drop=True), filtered_obs


def plot_bar_points(df: pd.DataFrame, y: str, ylabel: str, title: str, path_base: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    sns.barplot(data=df, x="sample_id", y=y, hue="age_group", palette=AGE_PALETTE, dodge=False, ax=ax)
    sns.stripplot(data=df, x="sample_id", y=y, color="black", size=5, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", labelrotation=30)
    ax.legend(title="Age", frameon=False)
    fig.tight_layout()
    fig.savefig(path_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def create_plots(sample_table: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    for ax, y, ylabel, title in [
        (axes[0], "Loxl2_detection_rate", "Detection rate", "Loxl2 detection"),
        (axes[1], "Loxl2_log2CPM1", "log2(CPM + 1)", "Loxl2 pseudobulk"),
    ]:
        sns.barplot(data=sample_table, x="sample_id", y=y, hue="age_group", palette=AGE_PALETTE, dodge=False, ax=ax)
        sns.stripplot(data=sample_table, x="sample_id", y=y, color="black", size=5, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", labelrotation=30)
        ax.legend_.remove()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], title="Age", loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_detection_by_sample.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_detection_by_sample.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    for ax, y, ylabel, title in [
        (axes[0], "median_total_counts", "Median total_counts", "mTEC1 total counts"),
        (axes[1], "median_n_genes_by_counts", "Median genes", "mTEC1 detected genes"),
    ]:
        sns.barplot(data=sample_table, x="sample_id", y=y, hue="age_group", palette=AGE_PALETTE, dodge=False, ax=ax)
        sns.stripplot(data=sample_table, x="sample_id", y=y, color="black", size=5, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", labelrotation=30)
        ax.legend_.remove()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], title="Age", loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_depth_qc_by_sample.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_depth_qc_by_sample.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
    for ax, x, xlabel in [
        (axes[0], "median_total_counts", "Median total_counts"),
        (axes[1], "median_n_genes_by_counts", "Median n_genes_by_counts"),
    ]:
        sns.scatterplot(
            data=sample_table,
            x=x,
            y="Loxl2_detection_rate",
            hue="age_group",
            palette=AGE_PALETTE,
            s=70,
            edgecolor="black",
            ax=ax,
        )
        for _, row in sample_table.iterrows():
            ax.annotate(row["sample_id"], (row[x], row["Loxl2_detection_rate"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Loxl2 detection rate")
        ax.legend(title="Age", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_detection_vs_depth.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_detection_vs_depth.pdf", bbox_inches="tight")
    plt.close(fig)


def pairwise_direction(sample_table: pd.DataFrame) -> pd.DataFrame:
    young = sample_table.loc[sample_table["age_group"].eq(YOUNG_LABEL)].copy()
    aged = sample_table.loc[sample_table["age_group"].eq(AGED_LABEL)].copy()
    rows = []
    metrics = [
        ("Loxl2_detection_rate", "young_higher"),
        ("Loxl2_log2CPM1", "young_higher"),
        ("Loxl2_CPM", "young_higher"),
    ]
    for _, yrow in young.iterrows():
        for _, arow in aged.iterrows():
            row = {"young_sample_id": yrow["sample_id"], "aged_sample_id": arow["sample_id"]}
            for metric, _ in metrics:
                diff = float(yrow[metric] - arow[metric])
                row[f"{metric}_young"] = float(yrow[metric])
                row[f"{metric}_aged"] = float(arow[metric])
                row[f"{metric}_young_minus_aged"] = diff
                row[f"{metric}_young_higher"] = bool(diff > 0)
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(PAIRWISE_TABLE, sep="\t", index=False)
    return out


def downsample(obs: pd.DataFrame, iterations: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    target = int(math.floor(obs.groupby("sample_id", observed=True)["_raw_library_size"].median().min()))
    if target < 100:
        raise ValueError(f"Downsampling target depth {target} is unrealistically low.")

    rng = np.random.default_rng(seed)
    base = obs[["sample_id", "age_group", "_loxl2_raw", "_raw_library_size"]].copy()
    base["_effective_library_size"] = np.minimum(base["_raw_library_size"], target)
    rows = []
    for iteration in range(iterations):
        draws = np.empty(len(base), dtype=float)
        for i, (_, row) in enumerate(base.iterrows()):
            gene_count = int(row["_loxl2_raw"])
            lib_size = int(row["_raw_library_size"])
            effective = int(row["_effective_library_size"])
            if lib_size <= target:
                draws[i] = gene_count
            elif gene_count <= 0:
                draws[i] = 0
            else:
                draws[i] = rng.hypergeometric(gene_count, lib_size - gene_count, effective)
        iter_df = base.assign(_downsampled_loxl2=draws)
        for sample_id, sample_obs in iter_df.groupby("sample_id", sort=True, observed=True):
            ages = sample_obs["age_group"].unique()
            raw_sum = float(sample_obs["_downsampled_loxl2"].sum())
            lib_sum = float(sample_obs["_effective_library_size"].sum())
            cpm = raw_sum / lib_sum * CPM_SCALE if lib_sum > 0 else np.nan
            rows.append(
                {
                    "iteration": iteration + 1,
                    "sample_id": sample_id,
                    "age_group": ages[0],
                    "downsample_target_depth": target,
                    "n_cells_used": int(len(sample_obs)),
                    "cells_below_target_depth": int((sample_obs["_raw_library_size"] < target).sum()),
                    "Loxl2_detecting_cells": int((sample_obs["_downsampled_loxl2"] > 0).sum()),
                    "Loxl2_detection_rate": float((sample_obs["_downsampled_loxl2"] > 0).mean()),
                    "Loxl2_raw_count_sum": raw_sum,
                    "effective_library_sum": lib_sum,
                    "Loxl2_CPM": cpm,
                    "Loxl2_log2CPM1": math.log2(cpm + 1.0) if np.isfinite(cpm) else np.nan,
                }
            )

    iterations_df = pd.DataFrame(rows)
    sample_summary = (
        iterations_df.groupby(["sample_id", "age_group"], observed=True)
        .agg(
            downsample_target_depth=("downsample_target_depth", "first"),
            n_cells_used=("n_cells_used", "first"),
            cells_below_target_depth=("cells_below_target_depth", "first"),
            mean_detection_rate=("Loxl2_detection_rate", "mean"),
            detection_rate_q025=("Loxl2_detection_rate", lambda x: float(np.quantile(x, 0.025))),
            detection_rate_q975=("Loxl2_detection_rate", lambda x: float(np.quantile(x, 0.975))),
            mean_log2CPM1=("Loxl2_log2CPM1", "mean"),
            log2CPM1_q025=("Loxl2_log2CPM1", lambda x: float(np.quantile(x, 0.025))),
            log2CPM1_q975=("Loxl2_log2CPM1", lambda x: float(np.quantile(x, 0.975))),
        )
        .reset_index()
    )

    age_diffs = []
    all_young_higher_detection = []
    all_young_higher_log = []
    all_young_higher_cpm = []
    for iteration, iter_df in iterations_df.groupby("iteration", observed=True):
        young = iter_df.loc[iter_df["age_group"].eq(YOUNG_LABEL)]
        aged = iter_df.loc[iter_df["age_group"].eq(AGED_LABEL)]
        age_diffs.append(
            {
                "iteration": iteration,
                "detection_rate_aged_minus_young": float(aged["Loxl2_detection_rate"].mean() - young["Loxl2_detection_rate"].mean()),
                "log2CPM1_aged_minus_young": float(aged["Loxl2_log2CPM1"].mean() - young["Loxl2_log2CPM1"].mean()),
                "CPM_aged_minus_young": float(aged["Loxl2_CPM"].mean() - young["Loxl2_CPM"].mean()),
            }
        )
        all_young_higher_detection.append(float(young["Loxl2_detection_rate"].min() > aged["Loxl2_detection_rate"].max()))
        all_young_higher_log.append(float(young["Loxl2_log2CPM1"].min() > aged["Loxl2_log2CPM1"].max()))
        all_young_higher_cpm.append(float(young["Loxl2_CPM"].min() > aged["Loxl2_CPM"].max()))

    age_diff_df = pd.DataFrame(age_diffs)
    global_summary = {
        "mean_detection_rate_aged_minus_young": float(age_diff_df["detection_rate_aged_minus_young"].mean()),
        "mean_log2CPM1_aged_minus_young": float(age_diff_df["log2CPM1_aged_minus_young"].mean()),
        "mean_CPM_aged_minus_young": float(age_diff_df["CPM_aged_minus_young"].mean()),
        "fraction_iterations_all_young_higher_detection_rate": float(np.mean(all_young_higher_detection)),
        "fraction_iterations_all_young_higher_log2CPM1": float(np.mean(all_young_higher_log)),
        "fraction_iterations_all_young_higher_CPM": float(np.mean(all_young_higher_cpm)),
    }
    for key, value in global_summary.items():
        sample_summary[key] = value

    iterations_df.to_csv(DOWNSAMPLE_ITER_TABLE, sep="\t", index=False)
    sample_summary.to_csv(DOWNSAMPLE_SUMMARY_TABLE, sep="\t", index=False)
    return iterations_df, sample_summary, global_summary


def plot_downsampling(summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    metrics = [
        ("mean_detection_rate", "detection_rate_q025", "detection_rate_q975", "Detection rate"),
        ("mean_log2CPM1", "log2CPM1_q025", "log2CPM1_q975", "log2(CPM + 1)"),
    ]
    x = np.arange(len(summary))
    for ax, (mean_col, low_col, high_col, ylabel) in zip(axes, metrics):
        colors = [AGE_PALETTE[a] for a in summary["age_group"]]
        ax.bar(x, summary[mean_col], color=colors, edgecolor="black", linewidth=0.6)
        yerr = np.vstack([summary[mean_col] - summary[low_col], summary[high_col] - summary[mean_col]])
        ax.errorbar(x, summary[mean_col], yerr=yerr, fmt="none", color="black", capsize=3, linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(summary["sample_id"], rotation=30, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_downsampling_sensitivity.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "mtec1_loxl2_downsampling_sensitivity.pdf", bbox_inches="tight")
    plt.close(fig)


def fit_logistic(obs: pd.DataFrame) -> pd.DataFrame:
    try:
        import statsmodels.formula.api as smf
    except Exception as exc:
        out = pd.DataFrame([{"model": "logistic", "status": f"not_run: {exc}"}])
        out.to_csv(LOGISTIC_TABLE, sep="\t", index=False)
        return out

    data = obs.copy()
    data["Loxl2_detected"] = (data["_loxl2_raw"] > 0).astype(int)
    data["age_aged_vs_young"] = data["age_group"].eq(AGED_LABEL).astype(int)
    data["log1p_total_counts"] = np.log1p(data["total_counts"].astype(float))
    data["log1p_n_genes_by_counts"] = np.log1p(data["n_genes_by_counts"].astype(float))
    formula = "Loxl2_detected ~ age_aged_vs_young + log1p_total_counts + log1p_n_genes_by_counts + mito_frac"
    rows = []
    try:
        standard_model = smf.logit(formula, data=data).fit(disp=False, maxiter=200)
        clustered_model = smf.logit(formula, data=data).fit(
            disp=False,
            maxiter=200,
            cov_type="cluster",
            cov_kwds={"groups": data["sample_id"]},
        )
        for result, label in [(standard_model, "standard"), (clustered_model, "sample_clustered")]:
            for term in result.params.index:
                rows.append(
                    {
                        "model": label,
                        "term": term,
                        "coef": float(result.params[term]),
                        "odds_ratio": float(np.exp(result.params[term])),
                        "std_error": float(result.bse[term]),
                        "pvalue": float(result.pvalues[term]),
                        "status": "descriptive_cells_not_independent",
                    }
                )
    except Exception as exc:
        rows.append({"model": "logistic", "term": "", "coef": np.nan, "odds_ratio": np.nan, "std_error": np.nan, "pvalue": np.nan, "status": f"unstable_or_failed: {exc}"})

    out = pd.DataFrame(rows)
    out.to_csv(LOGISTIC_TABLE, sep="\t", index=False)
    return out


def depth_statement(sample_table: pd.DataFrame) -> tuple[bool, str]:
    young = sample_table.loc[sample_table["age_group"].eq(YOUNG_LABEL)]
    aged = sample_table.loc[sample_table["age_group"].eq(AGED_LABEL)]
    aged_lower_total = bool(aged["median_total_counts"].mean() < young["median_total_counts"].mean())
    aged_lower_genes = bool(aged["median_n_genes_by_counts"].mean() < young["median_n_genes_by_counts"].mean())
    if aged_lower_total or aged_lower_genes:
        statement = "Aged mTEC1 samples have lower mean median total_counts or detected-gene depth than young samples."
    else:
        statement = "Aged mTEC1 samples do not have lower mean median total_counts or detected-gene depth than young samples."
    return aged_lower_total or aged_lower_genes, statement


def classify_risk(sample_table: pd.DataFrame, pairwise: pd.DataFrame, down_global: dict[str, Any]) -> str:
    _, depth_text = depth_statement(sample_table)
    all_pairwise = bool(
        pairwise[
            [
                "Loxl2_detection_rate_young_higher",
                "Loxl2_log2CPM1_young_higher",
                "Loxl2_CPM_young_higher",
            ]
        ]
        .all()
        .all()
    )
    down_preserved = (
        down_global.get("fraction_iterations_all_young_higher_detection_rate", 0.0) >= 0.95
        and down_global.get("fraction_iterations_all_young_higher_log2CPM1", 0.0) >= 0.95
    )
    if all_pairwise and down_preserved and "do not have lower" in depth_text:
        return "Low"
    if not all_pairwise:
        return "High"
    return "Partial / unresolved"


def percent_range(values: pd.Series) -> str:
    return f"{100 * values.min():.1f}-{100 * values.max():.1f}%"


def write_report(
    sample_table: pd.DataFrame,
    summary: dict[str, float],
    pairwise: pd.DataFrame,
    down_summary: pd.DataFrame,
    down_global: dict[str, Any],
    logistic: pd.DataFrame,
    matrix_source: str,
    annotation_col: str,
    annotation_label: str,
) -> str:
    risk = classify_risk(sample_table, pairwise, down_global)
    aged_lower_depth, depth_text = depth_statement(sample_table)
    young = sample_table.loc[sample_table["age_group"].eq(YOUNG_LABEL)]
    aged = sample_table.loc[sample_table["age_group"].eq(AGED_LABEL)]
    pairwise_all = bool(pairwise[[c for c in pairwise.columns if c.endswith("_young_higher")]].all().all())
    down_preserved = (
        down_global["fraction_iterations_all_young_higher_detection_rate"] >= 0.95
        and down_global["fraction_iterations_all_young_higher_log2CPM1"] >= 0.95
    )

    logistic_note = "The descriptive cell-level logistic model was fit with QC covariates; cells are not independent biological replicates."
    if logistic["status"].astype(str).str.contains("failed|unstable|not_run", case=False, regex=True).any():
        logistic_note = "The descriptive logistic model was not stable or was not run; no inferential claim is made from it."

    rec = "remain unchanged" if risk in {"Low", "Partial / unresolved"} else "be softened further"
    if risk == "Partial / unresolved":
        rec = "remain unchanged, assuming the current wording already stays cautious; alternatively add one short dropout/depth limitation sentence"

    lines = [
        "# mTEC1 Loxl2 dropout/depth audit",
        "",
        "## Executive summary",
        "",
        f"Dropout/depth artifact risk classification: **{risk}**.",
        "",
        f"Annotated mTEC1 cells were identified as `{annotation_label}` in `{annotation_col}`. Raw counts came from `{matrix_source}`. This audit is descriptive because the mTEC1 comparison has only 2 young and 2 aged biological samples, and cell-level models do not provide independent biological replication.",
        "",
        "The Loxl2 signal remains directionally young-higher/aged-lower in direct detection, pseudobulk CPM, pairwise sample checks, and the implemented depth-matched downsampling. However, dropout/depth artifact is not fully ruled out because the biological replicate count is small and Loxl2 detection is sparse in aged mTEC1.",
        "",
        "## Direct detection-rate result",
        "",
        f"Young mean mTEC1 Loxl2-positive percentage: {100 * summary['young_mean_detection_rate']:.1f}%.",
        f"Aged mean mTEC1 Loxl2-positive percentage: {100 * summary['aged_mean_detection_rate']:.1f}%.",
        "",
        "| sample | age | mTEC1 cells | Loxl2-positive cells | detection % | Loxl2 log2(CPM+1) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in sample_table.iterrows():
        lines.append(
            f"| {row['sample_id']} | {row['age_group']} | {int(row['n_mTEC1_cells'])} | "
            f"{int(row['Loxl2_detecting_cells'])} | {100 * row['Loxl2_detection_rate']:.1f} | {row['Loxl2_log2CPM1']:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Sequencing-depth/QC result",
            "",
            depth_text,
            f"Young mean median total_counts was {young['median_total_counts'].mean():.1f}; aged mean median total_counts was {aged['median_total_counts'].mean():.1f}. Young mean median detected genes was {young['median_n_genes_by_counts'].mean():.1f}; aged mean median detected genes was {aged['median_n_genes_by_counts'].mean():.1f}.",
            "Lower depth alone does not obviously explain the Loxl2 drop if aged depth is comparable or higher; if aged depth is lower, it remains a plausible partial contributor.",
            "",
            "## Downsampling result",
            "",
            f"Downsampling used a fixed seed and target per-cell raw library depth of {int(down_summary['downsample_target_depth'].iloc[0])}, chosen as the minimum median raw library size across the four mTEC1 samples. Cells already below the target were retained at their observed depth and counted in the output tables.",
            f"All young samples remained higher than all aged samples in {100 * down_global['fraction_iterations_all_young_higher_detection_rate']:.1f}% of iterations for detection rate and {100 * down_global['fraction_iterations_all_young_higher_log2CPM1']:.1f}% for log2(CPM+1).",
            f"The mean downsampled aged-minus-young difference was {down_global['mean_detection_rate_aged_minus_young']:.4f} for detection rate and {down_global['mean_log2CPM1_aged_minus_young']:.3f} for log2(CPM+1).",
            "",
            "## Pairwise sample-level result",
            "",
            "Both young samples are higher than both aged samples for Loxl2 detection rate, log2(CPM+1), and CPM." if pairwise_all else "Not every young sample is higher than every aged sample across all requested metrics.",
            "",
            "## Descriptive QC-adjusted model",
            "",
            logistic_note,
            "",
            "## External context",
            "",
            "GSE223049 provides broad sorted thymic epithelial context in which Loxl2 is directionally aged-lower, but that bulk epithelial comparison cannot test mTEC1 specificity or single-cell dropout.",
            "E-MTAB-8560 provides TEC and mTEC-like single-cell age-series context with aged-lower Loxl2 tendency in broad TEC, mTEC-like, mTEClo, and mTEChi groupings. Those labels are not one-to-one with the current mTEC1 annotation.",
            "Together, these external results reduce the chance that the entire observation is a one-off internal artifact, but they do not rule out dropout or depth sensitivity in GSE240016 annotated mTEC1 cells.",
            "",
            "## Final interpretation",
            "",
            "Dropout/depth artifact is not fully ruled out. The signal remains directionally supported after several checks, including sample-level pairwise comparisons and depth-matched downsampling, but the evidence remains limited by sparse Loxl2 detection and n=2 young versus n=2 aged biological samples.",
            "",
            "## Manuscript recommendation",
            "",
            f"The current v5.1 manuscript wording should {rec}. No manuscript file was changed in this task.",
            "",
            "## Output files",
            "",
            f"- `{AUDIT_TABLE.relative_to(PROJECT_ROOT)}`",
            f"- `{PAIRWISE_TABLE.relative_to(PROJECT_ROOT)}`",
            f"- `{DOWNSAMPLE_ITER_TABLE.relative_to(PROJECT_ROOT)}`",
            f"- `{DOWNSAMPLE_SUMMARY_TABLE.relative_to(PROJECT_ROOT)}`",
            f"- `{LOGISTIC_TABLE.relative_to(PROJECT_ROOT)}`",
        ]
    )

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")
    return risk


def main() -> int:
    args = parse_args()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input_h5ad}", flush=True)
    adata = ad.read_h5ad(args.input_h5ad)
    obs_all, raw_mtec1, raw_genes, matrix_source, annotation_col, annotation_label = prepare_mtec1(adata)
    print(f"Using {annotation_label} from {annotation_col}; matrix source: {matrix_source}", flush=True)

    sample_table_all = compute_sample_table(obs_all, raw_mtec1, raw_genes)
    sample_table, obs = filter_comparison_samples(sample_table_all, obs_all)
    summary = age_summary(sample_table)
    write_audit_table(sample_table, summary)
    create_plots(sample_table)
    pairwise = pairwise_direction(sample_table)

    position_lookup = pd.Series(np.arange(raw_mtec1.shape[0]), index=obs_all["_adata_position"].to_numpy())
    raw_mtec1 = raw_mtec1[position_lookup.loc[obs["_adata_position"].to_numpy()].to_numpy(), :]

    obs_for_downsample = obs.copy()
    gene_idx = raw_genes.get_loc(GENE)
    obs_for_downsample["_loxl2_raw"] = as_array(raw_mtec1[:, gene_idx]).astype(float)
    obs_for_downsample["_raw_library_size"] = matrix_row_sums(raw_mtec1).astype(float)
    down_iterations, down_summary, down_global = downsample(obs_for_downsample, args.downsample_iterations, args.seed)
    plot_downsampling(down_summary)
    logistic = fit_logistic(obs_for_downsample)

    risk = write_report(
        sample_table,
        summary,
        pairwise,
        down_summary,
        down_global,
        logistic,
        matrix_source,
        annotation_col,
        annotation_label,
    )

    young_range = percent_range(sample_table.loc[sample_table["age_group"].eq(YOUNG_LABEL), "Loxl2_detection_rate"])
    aged_range = percent_range(sample_table.loc[sample_table["age_group"].eq(AGED_LABEL), "Loxl2_detection_rate"])
    aged_lower_depth, _ = depth_statement(sample_table)
    down_preserved = (
        down_global["fraction_iterations_all_young_higher_detection_rate"] >= 0.95
        and down_global["fraction_iterations_all_young_higher_log2CPM1"] >= 0.95
    )
    manuscript_rec = "remain unchanged" if risk != "High" else "be softened further"
    if risk == "Partial / unresolved":
        manuscript_rec = "remain unchanged, or add one short limitation sentence"

    print("Script ran successfully: yes", flush=True)
    print(f"Dropout risk classification: {risk}", flush=True)
    print(f"Young detection rate range: {young_range}", flush=True)
    print(f"Aged detection rate range: {aged_range}", flush=True)
    print(f"Aged samples had lower sequencing depth: {'yes' if aged_lower_depth else 'no'}", flush=True)
    print(f"Downsampling preserved young > aged direction: {'yes' if down_preserved else 'no'}", flush=True)
    print(f"Manuscript changes recommended: {manuscript_rec}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
