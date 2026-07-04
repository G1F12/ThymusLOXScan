#!/usr/bin/env python
"""Matched-gene falsification analysis for sparse mTEC1 Loxl2 detection."""

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
RAW_FALLBACK_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"

TABLE_DIR = PROJECT_ROOT / "results" / "tables"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "falsification"
REPORT_PATH = PROJECT_ROOT / "reports" / "mtec1_loxl2_matched_gene_falsification.md"

ALL_GENES_TABLE = TABLE_DIR / "mtec1_loxl2_matched_gene_falsification_all_genes.tsv"
MATCHED_GENES_TABLE = TABLE_DIR / "mtec1_loxl2_matched_gene_falsification_matched_genes.tsv"
SUMMARY_TABLE = TABLE_DIR / "mtec1_loxl2_matched_gene_falsification_summary.tsv"

GENE = "Loxl2"
LOX_FAMILY = {"Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"}
MTEC1_LABEL_CANDIDATES = ("13:mTEC1", "mTEC1")
YOUNG_LABEL = "02mo"
AGED_LABEL = "18mo"
CPM_SCALE = 1_000_000.0
MIN_MTEC1_CELLS_PER_SAMPLE = 10

WINDOWS = {
    "strict": {"det_min": 0.05, "det_max": 0.15, "log2cpm_window": 0.5},
    "medium": {"det_min": 0.03, "det_max": 0.20, "log2cpm_window": 1.0},
    "broad": {"det_min": 0.01, "det_max": 0.30, "log2cpm_window": 1.5},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--min-cells-per-sample", type=int, default=MIN_MTEC1_CELLS_PER_SAMPLE)
    return parser.parse_args()


def as_array(values: Any) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def choose_input(path: Path) -> Path:
    if path.exists():
        return path
    if RAW_FALLBACK_H5AD.exists():
        return RAW_FALLBACK_H5AD
    raise FileNotFoundError(f"Neither {path} nor {RAW_FALLBACK_H5AD} exists.")


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


def choose_age_column(obs: pd.DataFrame) -> str:
    for column in ("age_group", "stage", "age", "Age"):
        if column in obs.columns:
            values = set(obs[column].astype(str).unique())
            if YOUNG_LABEL in values and AGED_LABEL in values:
                return column
    for column in obs.columns:
        values = set(obs[column].astype(str).unique())
        if YOUNG_LABEL in values and AGED_LABEL in values:
            return column
    raise ValueError(f"Could not identify an age/stage column containing {YOUNG_LABEL} and {AGED_LABEL}.")


def choose_sample_column(obs: pd.DataFrame) -> str:
    preferred = ("sample", "sample_id", "mouse_id", "donor", "replicate", "library")
    for column in preferred:
        if column in obs.columns and obs[column].nunique() >= 2:
            return column
    candidates = [c for c in obs.columns if any(token in c.lower() for token in ("sample", "mouse", "donor"))]
    for column in candidates:
        if obs[column].nunique() >= 2:
            return column
    raise ValueError("Could not identify a biological sample column in AnnData obs.")


def raw_matrix_and_genes(adata: ad.AnnData) -> tuple[Any, pd.Index, str]:
    if adata.raw is not None and adata.raw.X is not None:
        return adata.raw.X, pd.Index(adata.raw.var_names.astype(str)), "adata.raw.X"
    for layer_name in ("counts", "raw_counts", "raw", "spliced"):
        if layer_name in adata.layers:
            return adata.layers[layer_name], pd.Index(adata.var_names.astype(str)), f"adata.layers['{layer_name}']"
    return adata.X, pd.Index(adata.var_names.astype(str)), "adata.X (raw-count limitation: no raw/count layer found)"


def prepare_mtec1(adata: ad.AnnData) -> tuple[pd.DataFrame, Any, pd.Index, dict[str, str]]:
    annotation_col, annotation_label = choose_mtec1_annotation(adata.obs)
    age_col = choose_age_column(adata.obs)
    sample_col = choose_sample_column(adata.obs)
    raw_x, raw_genes, matrix_source = raw_matrix_and_genes(adata)
    if GENE not in raw_genes:
        raise KeyError(f"{GENE} was not found in {matrix_source} var names.")

    mask = adata.obs[annotation_col].astype(str).eq(annotation_label).to_numpy()
    positions = np.flatnonzero(mask)
    obs = adata.obs.iloc[positions].copy()
    obs["_adata_position"] = positions
    obs["sample_id"] = obs[sample_col].astype(str)
    obs["age_group"] = obs[age_col].astype(str)
    obs = obs.loc[obs["age_group"].isin([YOUNG_LABEL, AGED_LABEL])].copy()
    positions = obs["_adata_position"].to_numpy()
    obs["_matrix_position"] = np.arange(len(obs))

    metadata = {
        "annotation_col": annotation_col,
        "annotation_label": annotation_label,
        "age_col": age_col,
        "sample_col": sample_col,
        "matrix_source": matrix_source,
    }
    return obs, raw_x[positions, :], raw_genes, metadata


def summarize_per_sample(obs: pd.DataFrame, matrix: Any, genes: pd.Index) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_rows = []
    gene_frames = []
    for sample_id, sample_obs in obs.groupby("sample_id", sort=True, observed=True):
        ages = sample_obs["age_group"].unique()
        if len(ages) != 1:
            raise ValueError(f"Sample {sample_id} has multiple age labels: {ages}")
        sample_positions = sample_obs["_matrix_position"].to_numpy()
        sub = matrix[sample_positions, :]
        n_cells = sub.shape[0]
        counts = as_array(sub.sum(axis=0)).astype(float)
        detecting = as_array((sub > 0).sum(axis=0)).astype(float)
        library_sum = float(counts.sum())
        cpm = counts / library_sum * CPM_SCALE if library_sum > 0 else np.full_like(counts, np.nan)
        df = pd.DataFrame(
            {
                "gene": genes.to_numpy(),
                "sample_id": sample_id,
                "age_group": ages[0],
                "n_cells": int(n_cells),
                "n_detecting_cells": detecting.astype(int),
                "detection_rate": detecting / n_cells if n_cells else np.nan,
                "summed_raw_counts": counts,
                "raw_library_sum": library_sum,
                "CPM": cpm,
                "log2CPM1": np.log2(cpm + 1.0),
            }
        )
        gene_frames.append(df)
        sample_rows.append({"sample_id": sample_id, "age_group": ages[0], "n_cells": int(n_cells)})
    return pd.concat(gene_frames, ignore_index=True), pd.DataFrame(sample_rows)


def gene_statistics(per_sample: pd.DataFrame) -> pd.DataFrame:
    metric_frames = []
    sample_vectors = []
    for gene, df in per_sample.groupby("gene", sort=False, observed=True):
        young = df.loc[df["age_group"].eq(YOUNG_LABEL)].sort_values("sample_id")
        aged = df.loc[df["age_group"].eq(AGED_LABEL)].sort_values("sample_id")
        if young.empty or aged.empty:
            continue
        young_det = float(young["detection_rate"].mean())
        aged_det = float(aged["detection_rate"].mean())
        young_log = float(young["log2CPM1"].mean())
        aged_log = float(aged["log2CPM1"].mean())
        detection_ratio = aged_det / young_det if young_det > 0 else np.nan
        metric_frames.append(
            {
                "gene": gene,
                "young_mean_detection": young_det,
                "aged_mean_detection": aged_det,
                "detection_delta": aged_det - young_det,
                "detection_ratio": detection_ratio,
                "young_mean_log2CPM": young_log,
                "aged_mean_log2CPM": aged_log,
                "log2CPM_delta": aged_log - young_log,
                "all_young_above_all_aged_detection": bool(young["detection_rate"].min() > aged["detection_rate"].max()),
                "all_young_above_all_aged_log2CPM": bool(young["log2CPM1"].min() > aged["log2CPM1"].max()),
                "n_young_samples": int(young["sample_id"].nunique()),
                "n_aged_samples": int(aged["sample_id"].nunique()),
            }
        )
        sample_vectors.append(
            {
                "gene": gene,
                "young_detection_values": ";".join(f"{x:.6g}" for x in young["detection_rate"]),
                "aged_detection_values": ";".join(f"{x:.6g}" for x in aged["detection_rate"]),
                "young_log2CPM_values": ";".join(f"{x:.6g}" for x in young["log2CPM1"]),
                "aged_log2CPM_values": ";".join(f"{x:.6g}" for x in aged["log2CPM1"]),
            }
        )
    return pd.DataFrame(metric_frames).merge(pd.DataFrame(sample_vectors), on="gene", how="left")


def is_mito_gene(gene: str) -> bool:
    g = gene.upper()
    return g.startswith("MT-") or g.startswith("MT.")


def is_ribo_gene(gene: str) -> bool:
    g = gene.upper()
    return g.startswith("RPS") or g.startswith("RPL") or g.startswith("MRPS") or g.startswith("MRPL")


def build_background(all_genes: pd.DataFrame) -> pd.DataFrame:
    out = all_genes.copy()
    out["excluded_lox_family"] = out["gene"].isin(LOX_FAMILY)
    out["excluded_mitochondrial"] = out["gene"].map(is_mito_gene)
    out["excluded_ribosomal"] = out["gene"].map(is_ribo_gene)
    out["background_included"] = (
        ~out["excluded_lox_family"]
        & ~out["excluded_mitochondrial"]
        & ~out["excluded_ribosomal"]
        & out["young_mean_detection"].gt(0)
    )
    return out


def empirical_p(mask: pd.Series) -> float:
    return float((int(mask.sum()) + 1) / (len(mask) + 1)) if len(mask) else np.nan


def summarize_windows(all_genes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    if GENE not in set(all_genes["gene"]):
        raise KeyError(f"{GENE} was not found in per-gene statistics.")
    loxl2 = all_genes.loc[all_genes["gene"].eq(GENE)].iloc[0]
    background = all_genes.loc[all_genes["background_included"]].copy()

    matched_frames = []
    summary_rows = []
    for name, spec in WINDOWS.items():
        match = background.loc[
            background["young_mean_detection"].between(spec["det_min"], spec["det_max"], inclusive="both")
            & background["young_mean_log2CPM"].between(
                loxl2["young_mean_log2CPM"] - spec["log2cpm_window"],
                loxl2["young_mean_log2CPM"] + spec["log2cpm_window"],
                inclusive="both",
            )
        ].copy()
        match.insert(0, "match_window", name)
        matched_frames.append(match)

        n = len(match)
        det_delta_extreme = match["detection_delta"].le(loxl2["detection_delta"])
        det_ratio_extreme = match["detection_ratio"].le(loxl2["detection_ratio"])
        log_delta_extreme = match["log2CPM_delta"].le(loxl2["log2CPM_delta"])
        det_order = match["all_young_above_all_aged_detection"].astype(bool)
        log_order = match["all_young_above_all_aged_log2CPM"].astype(bool)
        det_both = det_delta_extreme & det_order
        log_both = log_delta_extreme & log_order

        summary_rows.append(
            {
                "match_window": name,
                "n_matched": n,
                "young_detection_min": spec["det_min"],
                "young_detection_max": spec["det_max"],
                "young_log2CPM_window": spec["log2cpm_window"],
                "loxl2_young_mean_detection": loxl2["young_mean_detection"],
                "loxl2_aged_mean_detection": loxl2["aged_mean_detection"],
                "loxl2_detection_delta": loxl2["detection_delta"],
                "loxl2_detection_ratio": loxl2["detection_ratio"],
                "loxl2_young_mean_log2CPM": loxl2["young_mean_log2CPM"],
                "loxl2_aged_mean_log2CPM": loxl2["aged_mean_log2CPM"],
                "loxl2_log2CPM_delta": loxl2["log2CPM_delta"],
                "loxl2_all_young_above_all_aged_detection": loxl2["all_young_above_all_aged_detection"],
                "loxl2_all_young_above_all_aged_log2CPM": loxl2["all_young_above_all_aged_log2CPM"],
                "fraction_detection_delta_le_loxl2": float(det_delta_extreme.mean()) if n else np.nan,
                "fraction_detection_ratio_le_loxl2": float(det_ratio_extreme.mean()) if n else np.nan,
                "fraction_log2CPM_delta_le_loxl2": float(log_delta_extreme.mean()) if n else np.nan,
                "fraction_all_young_above_all_aged_detection": float(det_order.mean()) if n else np.nan,
                "fraction_all_young_above_all_aged_log2CPM": float(log_order.mean()) if n else np.nan,
                "fraction_detection_delta_le_loxl2_and_ordered": float(det_both.mean()) if n else np.nan,
                "fraction_log2CPM_delta_le_loxl2_and_ordered": float(log_both.mean()) if n else np.nan,
                "empirical_p_detection_delta_le_loxl2": empirical_p(det_delta_extreme),
                "empirical_p_detection_ratio_le_loxl2": empirical_p(det_ratio_extreme),
                "empirical_p_log2CPM_delta_le_loxl2": empirical_p(log_delta_extreme),
                "empirical_p_detection_delta_le_loxl2_and_ordered": empirical_p(det_both),
                "empirical_p_log2CPM_delta_le_loxl2_and_ordered": empirical_p(log_both),
            }
        )
    matched = pd.concat(matched_frames, ignore_index=True) if matched_frames else pd.DataFrame()
    return pd.DataFrame(summary_rows), matched, loxl2


def classify(summary: pd.DataFrame) -> tuple[str, str]:
    usable = summary.loc[summary["n_matched"].ge(20)].copy()
    if usable.empty:
        return "matched-gene set too small / inconclusive", "The matched-gene background was too small for stable falsification inference; results are reported descriptively."
    row = usable.loc[usable["match_window"].eq("medium")]
    if row.empty:
        row = usable.iloc[[0]]
    ref = row.iloc[0]
    primary = ref["fraction_detection_delta_le_loxl2_and_ordered"]
    if pd.isna(primary):
        return "matched-gene set too small / inconclusive", "The matched-gene background was too small for stable falsification inference; results are reported descriptively."
    if primary <= 0.10:
        return (
            "Loxl2-like collapse uncommon among matched genes",
            "Matched-gene falsification suggested that the Loxl2 detection pattern was uncommon among similarly expressed mTEC1 genes, supporting prioritization as a transcript-level candidate while not excluding dropout/depth effects.",
        )
    if primary >= 0.25:
        return (
            "Loxl2-like collapse common among matched genes",
            "Matched-gene falsification showed that Loxl2-like detection collapse occurs frequently among similarly expressed mTEC1 genes, supporting the possibility that the mTEC1 Loxl2 pattern is partly technical.",
        )
    return "matched-gene set too small / inconclusive", "The matched-gene background was too small for stable falsification inference; results are reported descriptively."


def plot_histograms(matched: pd.DataFrame, loxl2: pd.Series) -> list[Path]:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    figure_paths: list[Path] = []
    for metric, xlabel, stem in [
        ("detection_delta", "Aged minus young detection rate", "mtec1_loxl2_matched_gene_detection_delta"),
        ("log2CPM_delta", "Aged minus young log2(CPM + 1)", "mtec1_loxl2_matched_gene_log2cpm_delta"),
    ]:
        fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4), sharey=False)
        for ax, window in zip(axes, WINDOWS):
            df = matched.loc[matched["match_window"].eq(window)]
            if df.empty:
                ax.text(0.5, 0.5, "No matched genes", ha="center", va="center", transform=ax.transAxes)
            else:
                sns.histplot(df[metric], bins=min(30, max(5, int(math.sqrt(len(df))))), color="#4C78A8", ax=ax)
            ax.axvline(float(loxl2[metric]), color="#D55E00", linestyle="--", linewidth=1.5, label="Loxl2")
            ax.set_title(f"{window} matched genes")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Gene count")
        axes[0].legend(frameon=False)
        fig.suptitle("Matched-gene descriptive falsification, not proof", y=1.04)
        fig.tight_layout()
        png = FIGURE_DIR / f"{stem}.png"
        pdf = FIGURE_DIR / f"{stem}.pdf"
        fig.savefig(png, dpi=300, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)
        figure_paths.extend([png, pdf])
    return figure_paths


def fmt(value: Any, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return "yes" if value else "no"
    return f"{float(value):.{digits}f}"


def write_report(
    summary: pd.DataFrame,
    matched: pd.DataFrame,
    loxl2: pd.Series,
    metadata: dict[str, str],
    sample_table: pd.DataFrame,
    excluded_samples: pd.DataFrame,
    classification: str,
    interpretation: str,
    figure_paths: list[Path],
) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    n_background = int(matched[["match_window", "gene"]].drop_duplicates()["gene"].nunique()) if not matched.empty else 0
    lines = [
        "# mTEC1 Loxl2 matched-gene falsification",
        "",
        "## Scope",
        "",
        "This is a matched-gene descriptive falsification analysis for GSE240016 annotated mTEC1 cells. It asks how often genes with broadly similar young mTEC1 detection and expression show a Loxl2-like young-high / aged-low collapse. It is not proof that the Loxl2 pattern is biological, and it does not exclude dropout or depth effects.",
        "",
        "## Inputs and metadata",
        "",
        f"- mTEC1 label: `{metadata['annotation_label']}` from `{metadata['annotation_col']}`",
        f"- age/stage column: `{metadata['age_col']}`",
        f"- biological sample column: `{metadata['sample_col']}`",
        f"- count matrix: `{metadata['matrix_source']}`",
        "- LOX-family genes, mitochondrial genes, ribosomal genes, and genes with zero young detection were excluded from the matched background.",
        "- Ribosomal genes were excluded from the background.",
        f"- Samples with fewer than {MIN_MTEC1_CELLS_PER_SAMPLE} annotated mTEC1 cells were excluded to match the prior dropout/depth audit design.",
        "",
        "| sample | age | mTEC1 cells |",
        "|---|---:|---:|",
    ]
    for _, row in sample_table.sort_values(["age_group", "sample_id"]).iterrows():
        lines.append(f"| {row['sample_id']} | {row['age_group']} | {int(row['n_cells'])} |")
    if not excluded_samples.empty:
        lines.extend(["", "Excluded low-cell samples:", "", "| sample | age | mTEC1 cells |", "|---|---:|---:|"])
        for _, row in excluded_samples.sort_values(["age_group", "sample_id"]).iterrows():
            lines.append(f"| {row['sample_id']} | {row['age_group']} | {int(row['n_cells'])} |")

    lines.extend(
        [
            "",
            "## Loxl2 observed statistics",
            "",
            "| statistic | value |",
            "|---|---:|",
            f"| young mean detection | {fmt(loxl2['young_mean_detection'])} |",
            f"| aged mean detection | {fmt(loxl2['aged_mean_detection'])} |",
            f"| detection delta, aged minus young | {fmt(loxl2['detection_delta'])} |",
            f"| detection ratio, aged / young | {fmt(loxl2['detection_ratio'])} |",
            f"| young mean log2(CPM + 1) | {fmt(loxl2['young_mean_log2CPM'])} |",
            f"| aged mean log2(CPM + 1) | {fmt(loxl2['aged_mean_log2CPM'])} |",
            f"| log2(CPM + 1) delta, aged minus young | {fmt(loxl2['log2CPM_delta'])} |",
            f"| all young samples above all aged samples, detection | {fmt(loxl2['all_young_above_all_aged_detection'])} |",
            f"| all young samples above all aged samples, log2CPM | {fmt(loxl2['all_young_above_all_aged_log2CPM'])} |",
            "",
            "## Matched-gene background",
            "",
            f"Unique genes appearing in at least one matched set: {n_background}.",
            "",
            "| window | matched genes | frac detection delta <= Loxl2 | empirical p | frac delta + ordered | frac log2CPM delta <= Loxl2 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['match_window']} | {int(row['n_matched'])} | "
            f"{fmt(row['fraction_detection_delta_le_loxl2'])} | "
            f"{fmt(row['empirical_p_detection_delta_le_loxl2'])} | "
            f"{fmt(row['fraction_detection_delta_le_loxl2_and_ordered'])} | "
            f"{fmt(row['fraction_log2CPM_delta_le_loxl2'])} |"
        )

    lines.extend(
        [
            "",
            "## Classification",
            "",
            f"Classification: **{classification}**.",
            "",
            interpretation,
            "",
            "## Output files",
            "",
            f"- `{ALL_GENES_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
            f"- `{MATCHED_GENES_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
            f"- `{SUMMARY_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
        ]
    )
    lines.extend(f"- `{path.relative_to(PROJECT_ROOT).as_posix()}`" for path in figure_paths)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    input_h5ad = choose_input(args.input_h5ad)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading {input_h5ad}", flush=True)
    adata = ad.read_h5ad(input_h5ad)
    obs, raw_mtec1, raw_genes, metadata = prepare_mtec1(adata)
    print(
        f"Using {metadata['annotation_label']} from {metadata['annotation_col']}; matrix source: {metadata['matrix_source']}",
        flush=True,
    )

    sample_counts = (
        obs.groupby(["sample_id", "age_group"], sort=True, observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    excluded_samples = sample_counts.loc[sample_counts["n_cells"].lt(args.min_cells_per_sample)].copy()
    eligible_samples = sample_counts.loc[sample_counts["n_cells"].ge(args.min_cells_per_sample), "sample_id"]
    obs = obs.loc[obs["sample_id"].isin(eligible_samples)].copy()
    if obs.empty:
        raise ValueError("No mTEC1 samples remained after minimum-cell filtering.")

    per_sample, sample_table = summarize_per_sample(obs, raw_mtec1, raw_genes)
    all_genes = gene_statistics(per_sample)
    all_genes = build_background(all_genes)
    summary, matched, loxl2 = summarize_windows(all_genes)
    classification, interpretation = classify(summary)
    figure_paths = plot_histograms(matched, loxl2)

    all_genes.to_csv(ALL_GENES_TABLE, sep="\t", index=False)
    matched.to_csv(MATCHED_GENES_TABLE, sep="\t", index=False)
    summary.to_csv(SUMMARY_TABLE, sep="\t", index=False)
    write_report(summary, matched, loxl2, metadata, sample_table, excluded_samples, classification, interpretation, figure_paths)

    print("Falsification test completed: yes", flush=True)
    print(f"Matched-gene classification: {classification}", flush=True)
    print(f"Loxl2 detection delta: {loxl2['detection_delta']:.6f}", flush=True)
    print(f"Loxl2 log2CPM delta: {loxl2['log2CPM_delta']:.6f}", flush=True)
    print(f"Saved report: {REPORT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
