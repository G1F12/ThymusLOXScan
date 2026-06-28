#!/usr/bin/env python
"""Summarize LOX-family expression in external GSE223049 thymus samples.

The external dataset is sorted bulk RNA-seq, so all statistics are descriptive
and are used only for broad directional validation.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data" / "external" / "GSE223049"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "tables" / "external_gse223049_lox_validation.tsv"
DEFAULT_STATS_OUTPUT = PROJECT_ROOT / "results" / "tables" / "external_gse223049_lox_validation_stats.tsv"
DEFAULT_SUMMARY_V2_OUTPUT = (
    PROJECT_ROOT / "results" / "tables" / "external_gse223049_lox_validation_summary_v2.tsv"
)
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "external_validation" / "gse223049"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "gse223049_external_validation_v2.md"

COUNTS_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223049/suppl/GSE223049_RNA_seq_counts_23_cell_types.txt.gz"
LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
TARGET_PREFIXES = ["Thymic_epithelial", "Thymic_fibroblasts"]
YOUNG_LABEL = "young_2mo"
AGED_LABEL = "aged_22_24mo"
AGE_ORDER = [YOUNG_LABEL, AGED_LABEL]
AGE_DISPLAY = {YOUNG_LABEL: "Young 2 mo", AGED_LABEL: "Aged 22-24 mo"}
AGE_PALETTE = {YOUNG_LABEL: "#0072B2", AGED_LABEL: "#D55E00"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stats-output", type=Path, default=DEFAULT_STATS_OUTPUT)
    parser.add_argument("--summary-v2-output", type=Path, default=DEFAULT_SUMMARY_V2_OUTPUT)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def download_counts(data_dir: Path) -> Path:
    data_dir.mkdir(parents=True, exist_ok=True)
    gz_path = data_dir / "GSE223049_RNA_seq_counts_23_cell_types.txt.gz"
    txt_path = data_dir / "GSE223049_RNA_seq_counts_23_cell_types.txt"
    if not gz_path.exists() and not txt_path.exists():
        print(f"Downloading {COUNTS_URL}", flush=True)
        urllib.request.urlretrieve(COUNTS_URL, gz_path)
    if not txt_path.exists():
        print(f"Decompressing {gz_path}", flush=True)
        with gzip.open(gz_path, "rb") as src, txt_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    return txt_path


def parse_sample(sample: str) -> tuple[str | None, str | None, str | None]:
    for prefix in TARGET_PREFIXES:
        if sample.startswith(prefix + "_"):
            suffix = sample.removeprefix(prefix + "_")
            parts = suffix.split("_")
            if len(parts) >= 2 and parts[0] in {"Y", "O"}:
                age = YOUNG_LABEL if parts[0] == "Y" else AGED_LABEL
                replicate = parts[1]
                return prefix, age, replicate
    return None, None, None


def standard_error(values: pd.Series) -> float:
    values = values.dropna().astype(float)
    if len(values) < 2:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def confidence_interval_delta(young: np.ndarray, aged: np.ndarray) -> tuple[float, float]:
    if len(young) < 2 or len(aged) < 2:
        return np.nan, np.nan
    young_var = float(np.var(young, ddof=1))
    aged_var = float(np.var(aged, ddof=1))
    se_delta = np.sqrt(young_var / len(young) + aged_var / len(aged))
    delta = float(np.mean(aged) - np.mean(young))
    if not np.isfinite(se_delta):
        return np.nan, np.nan
    if se_delta == 0:
        return delta, delta
    numerator = (young_var / len(young) + aged_var / len(aged)) ** 2
    denominator = (
        (young_var / len(young)) ** 2 / (len(young) - 1)
        + (aged_var / len(aged)) ** 2 / (len(aged) - 1)
    )
    if denominator <= 0:
        return np.nan, np.nan
    df = numerator / denominator
    critical = stats.t.ppf(0.975, df)
    return float(delta - critical * se_delta), float(delta + critical * se_delta)


def hedges_g(young: np.ndarray, aged: np.ndarray) -> float:
    if len(young) < 2 or len(aged) < 2:
        return np.nan
    young_var = float(np.var(young, ddof=1))
    aged_var = float(np.var(aged, ddof=1))
    pooled_var = ((len(young) - 1) * young_var + (len(aged) - 1) * aged_var) / (
        len(young) + len(aged) - 2
    )
    if pooled_var <= 0 or not np.isfinite(pooled_var):
        return np.nan
    cohen_d = (float(np.mean(aged)) - float(np.mean(young))) / np.sqrt(pooled_var)
    correction = 1 - (3 / (4 * (len(young) + len(aged)) - 9))
    return float(cohen_d * correction)


def mann_whitney_pvalue(young: np.ndarray, aged: np.ndarray) -> float:
    if len(young) < 2 or len(aged) < 2:
        return np.nan
    try:
        return float(stats.mannwhitneyu(aged, young, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def build_stats(sample_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cell_type, gene), df in sample_df.groupby(["external_cell_type", "gene"], sort=True):
        young = df.loc[df["age_group"].eq(YOUNG_LABEL), "log2_cpm_plus1"].dropna().to_numpy(float)
        aged = df.loc[df["age_group"].eq(AGED_LABEL), "log2_cpm_plus1"].dropna().to_numpy(float)
        ci_low, ci_high = confidence_interval_delta(young, aged)
        young_sd = float(np.std(young, ddof=1)) if len(young) >= 2 else np.nan
        aged_sd = float(np.std(aged, ddof=1)) if len(aged) >= 2 else np.nan
        young_mean = float(np.mean(young)) if len(young) else np.nan
        aged_mean = float(np.mean(aged)) if len(aged) else np.nan
        delta = aged_mean - young_mean if np.isfinite(young_mean) and np.isfinite(aged_mean) else np.nan
        rows.append(
            {
                "external_cell_type": cell_type,
                "gene": gene,
                "young_n": int(len(young)),
                "aged_n": int(len(aged)),
                "mean_young_log2_cpm_plus1": young_mean,
                "mean_aged_log2_cpm_plus1": aged_mean,
                "delta_log2_cpm_aged_minus_young": delta,
                "sd_young_log2_cpm_plus1": young_sd,
                "sd_aged_log2_cpm_plus1": aged_sd,
                "se_young_log2_cpm_plus1": young_sd / np.sqrt(len(young)) if len(young) >= 2 else np.nan,
                "se_aged_log2_cpm_plus1": aged_sd / np.sqrt(len(aged)) if len(aged) >= 2 else np.nan,
                "delta_95ci_low": ci_low,
                "delta_95ci_high": ci_high,
                "nonparametric_pvalue_mann_whitney": mann_whitney_pvalue(young, aged),
                "effect_size_hedges_g_aged_minus_young": hedges_g(young, aged),
                "statistics_scope": (
                    "descriptive only; sorted bulk RNA-seq and not subtype-resolved"
                ),
            }
        )
    return pd.DataFrame(rows)


def build_summary_v2(sample_df: pd.DataFrame, stats_df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        sample_df.groupby(["external_cell_type", "gene", "age_group"], observed=True)
        .agg(
            n_samples=("sample_id", "nunique"),
            mean_log2_cpm_plus1=("log2_cpm_plus1", "mean"),
            sd_log2_cpm_plus1=("log2_cpm_plus1", lambda x: x.dropna().std(ddof=1)),
            se_log2_cpm_plus1=("log2_cpm_plus1", standard_error),
            mean_cpm=("cpm", "mean"),
            mean_raw_counts=("raw_counts", "mean"),
        )
        .reset_index()
    )
    return summary.merge(stats_df, on=["external_cell_type", "gene"], how="left")


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> tuple[Path, Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    png = figure_dir / f"{stem}.png"
    pdf = figure_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_dot_box(
    sample_df: pd.DataFrame,
    figure_dir: Path,
    cell_type: str,
    genes: list[str],
    stem: str,
    title: str,
) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    df = sample_df.loc[
        sample_df["external_cell_type"].eq(cell_type) & sample_df["gene"].isin(genes)
    ].copy()
    df["gene"] = pd.Categorical(df["gene"], categories=genes, ordered=True)
    df["age_group"] = pd.Categorical(df["age_group"], categories=AGE_ORDER, ordered=True)

    fig_width = max(4.6, 2.55 * len(genes) + 1.2)
    fig, ax = plt.subplots(figsize=(fig_width, 4.2))
    sns.boxplot(
        data=df,
        x="gene",
        y="log2_cpm_plus1",
        hue="age_group",
        hue_order=AGE_ORDER,
        palette={key: "white" for key in AGE_ORDER},
        width=0.62,
        fliersize=0,
        linewidth=1.0,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="gene",
        y="log2_cpm_plus1",
        hue="age_group",
        hue_order=AGE_ORDER,
        dodge=True,
        jitter=0.08,
        palette=AGE_PALETTE,
        size=6,
        edgecolor="black",
        linewidth=0.45,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    legend_kwargs = {"frameon": False}
    if len(genes) == 1:
        legend_kwargs |= {"loc": "upper left", "bbox_to_anchor": (1.02, 1.0), "borderaxespad": 0}
    else:
        legend_kwargs |= {"loc": "best"}
    ax.legend(
        handles[-2:],
        [AGE_DISPLAY[label] for label in labels[-2:]],
        title="Age group",
        **legend_kwargs,
    )
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("log2(CPM + 1)")
    fig.tight_layout()
    return save_figure(fig, figure_dir, stem)


def plot_delta_summary(stats_df: pd.DataFrame, figure_dir: Path) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    df = stats_df.copy()
    df["external_cell_type_label"] = df["external_cell_type"].str.replace("_", " ", regex=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    sns.pointplot(
        data=df,
        x="gene",
        y="delta_log2_cpm_aged_minus_young",
        hue="external_cell_type_label",
        dodge=0.35,
        linestyle="none",
        errorbar=None,
        markers=["o", "s"],
        palette=["#009E73", "#CC79A7"],
        ax=ax,
    )
    for _, row in df.iterrows():
        x_base = LOX_GENES.index(row["gene"])
        offset = -0.12 if row["external_cell_type"] == "Thymic_epithelial" else 0.12
        if np.isfinite(row["delta_95ci_low"]) and np.isfinite(row["delta_95ci_high"]):
            ax.vlines(
                x_base + offset,
                row["delta_95ci_low"],
                row["delta_95ci_high"],
                color="0.35",
                linewidth=1.0,
                zorder=0,
            )
    ax.axhline(0, color="0.25", linewidth=0.9, linestyle="--")
    ax.set_title("GSE223049 LOX-family aged-minus-young deltas", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("Delta log2(CPM + 1), aged minus young")
    ax.legend(title="", frameon=False, loc="best")
    fig.tight_layout()
    return save_figure(fig, figure_dir, "gse223049_lox_family_delta_summary")


def fmt(value: float, digits: int = 3) -> str:
    return "NA" if pd.isna(value) else f"{value:.{digits}f}"


def write_report(stats_df: pd.DataFrame, figure_paths: list[Path], report_path: Path) -> None:
    def row(cell_type: str, gene: str) -> pd.Series:
        hit = stats_df.loc[stats_df["external_cell_type"].eq(cell_type) & stats_df["gene"].eq(gene)]
        if hit.empty:
            raise KeyError(f"Missing stats for {cell_type} {gene}")
        return hit.iloc[0]

    fb_lox = row("Thymic_fibroblasts", "Lox")
    fb_loxl1 = row("Thymic_fibroblasts", "Loxl1")
    fb_loxl2 = row("Thymic_fibroblasts", "Loxl2")
    epi_loxl2 = row("Thymic_epithelial", "Loxl2")

    lines = [
        "# GSE223049 external validation v2",
        "",
        "## Scope",
        "",
        "GSE223049 is used here as an independent, broad sorted bulk RNA-seq comparison of young 2-month and aged 22-24-month mouse thymic samples. The relevant sorted populations are thymic fibroblasts and thymic epithelial cells.",
        "",
        "All statistics in this report are descriptive. This dataset is not subtype-resolved and therefore cannot validate capsFB, medFB, or mTEC1 specificity. It should be interpreted as broad directional external context, not as a replacement for subtype-resolved validation.",
        "",
        "## Directional consistency",
        "",
        "| external comparison | young n | aged n | aged-minus-young delta log2(CPM+1) | descriptive p | interpretation |",
        "|---|---:|---:|---:|---:|---|",
        f"| Thymic fibroblast Lox | {int(fb_lox['young_n'])} | {int(fb_lox['aged_n'])} | {fmt(fb_lox['delta_log2_cpm_aged_minus_young'])} | {fmt(fb_lox['nonparametric_pvalue_mann_whitney'])} | Directionally consistent with broad/capsFB aged-lower Lox, but only broad fibroblast support. |",
        f"| Thymic fibroblast Loxl1 | {int(fb_loxl1['young_n'])} | {int(fb_loxl1['aged_n'])} | {fmt(fb_loxl1['delta_log2_cpm_aged_minus_young'])} | {fmt(fb_loxl1['nonparametric_pvalue_mann_whitney'])} | Not consistent with a medFB Loxl1 increase; medFB Loxl1 increase is not validated. |",
        f"| Thymic fibroblast Loxl2 | {int(fb_loxl2['young_n'])} | {int(fb_loxl2['aged_n'])} | {fmt(fb_loxl2['delta_log2_cpm_aged_minus_young'])} | {fmt(fb_loxl2['nonparametric_pvalue_mann_whitney'])} | Directionally consistent with aged-lower fibroblast/medFB Loxl2, but only broad fibroblast support. |",
        f"| Thymic epithelial Loxl2 | {int(epi_loxl2['young_n'])} | {int(epi_loxl2['aged_n'])} | {fmt(epi_loxl2['delta_log2_cpm_aged_minus_young'])} | {fmt(epi_loxl2['nonparametric_pvalue_mann_whitney'])} | Directionally consistent with aged-lower epithelial/mTEC1 Loxl2, but not mTEC1-specific. |",
        "",
        "## Findings not testable in GSE223049",
        "",
        "- capsFB specificity is not testable because thymic fibroblasts are pooled in sorted bulk RNA-seq.",
        "- medFB specificity is not testable because thymic fibroblast subtypes are not separated.",
        "- mTEC1 specificity is not testable because thymic epithelial cells are pooled.",
        "- Detection-rate or single-cell subtype mechanisms are not testable from this bulk count matrix.",
        "- The medFB Loxl1 increase is not validated; broad thymic fibroblast Loxl1 is slightly lower in aged samples in this dataset.",
        "",
        "## Figures",
        "",
    ]
    for path in figure_paths:
        lines.append(f"- `{path.relative_to(PROJECT_ROOT).as_posix()}`")

    lines.extend(
        [
            "",
            "## Bottom line",
            "",
            "GSE223049 supports broad external directionality for aged-lower thymic fibroblast Lox, aged-lower thymic fibroblast Loxl2, and aged-lower thymic epithelial Loxl2. It does not validate subtype-specific capsFB, medFB, or mTEC1 claims, and it does not validate the medFB Loxl1 increase.",
        ]
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    counts_path = download_counts(args.data_dir)
    counts = pd.read_csv(counts_path, sep="\t")
    gene_col = counts.columns[0]
    counts[gene_col] = counts[gene_col].astype(str)
    counts = counts.set_index(gene_col)

    target_samples = []
    sample_meta = {}
    for sample in counts.columns:
        cell_type, age, replicate = parse_sample(sample)
        if cell_type:
            target_samples.append(sample)
            sample_meta[sample] = {
                "external_cell_type": cell_type,
                "age_group": age,
                "replicate": replicate,
            }

    rows = []
    lib_sizes = counts[target_samples].sum(axis=0)
    for gene in LOX_GENES:
        if gene not in counts.index:
            for sample in target_samples:
                rows.append({"gene": gene, "sample_id": sample, "present": False, **sample_meta[sample]})
            continue
        for sample in target_samples:
            raw = float(counts.loc[gene, sample])
            cpm = raw / float(lib_sizes[sample]) * 1_000_000 if lib_sizes[sample] else np.nan
            rows.append(
                {
                    "gene": gene,
                    "sample_id": sample,
                    "present": True,
                    **sample_meta[sample],
                    "raw_counts": raw,
                    "library_size": float(lib_sizes[sample]),
                    "cpm": cpm,
                    "log2_cpm_plus1": np.log2(cpm + 1) if np.isfinite(cpm) else np.nan,
                }
            )

    sample_df = pd.DataFrame(rows)
    legacy_summary = (
        sample_df.groupby(["external_cell_type", "gene", "age_group"], observed=True)
        .agg(
            n_samples=("sample_id", "nunique"),
            mean_log2_cpm_plus1=("log2_cpm_plus1", "mean"),
            mean_cpm=("cpm", "mean"),
            mean_raw_counts=("raw_counts", "mean"),
        )
        .reset_index()
    )
    wide = legacy_summary.pivot_table(
        index=["external_cell_type", "gene"],
        columns="age_group",
        values="mean_log2_cpm_plus1",
        aggfunc="first",
    ).reset_index()
    wide["delta_log2_cpm_aged_minus_young"] = wide.get("aged_22_24mo", np.nan) - wide.get("young_2mo", np.nan)

    out = sample_df.merge(
        wide[["external_cell_type", "gene", "delta_log2_cpm_aged_minus_young"]],
        on=["external_cell_type", "gene"],
        how="left",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary_path = args.output.with_name(args.output.stem + "_summary.tsv")
    legacy_summary.merge(wide, on=["external_cell_type", "gene"], how="left").to_csv(
        summary_path, sep="\t", index=False
    )

    stats_df = build_stats(sample_df)
    summary_v2 = build_summary_v2(sample_df, stats_df)
    args.stats_output.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(args.stats_output, sep="\t", index=False)
    summary_v2.to_csv(args.summary_v2_output, sep="\t", index=False)

    figure_paths = []
    for paths in [
        plot_dot_box(
            sample_df,
            args.figure_dir,
            "Thymic_fibroblasts",
            ["Lox", "Loxl1", "Loxl2"],
            "gse223049_thymic_fibroblast_lox_dot_box",
            "Thymic fibroblast LOX-family expression",
        ),
        plot_dot_box(
            sample_df,
            args.figure_dir,
            "Thymic_epithelial",
            ["Loxl2"],
            "gse223049_thymic_epithelial_loxl2_dot_box",
            "Thymic epithelial Loxl2 expression",
        ),
        plot_delta_summary(stats_df, args.figure_dir),
    ]:
        figure_paths.extend(paths)
    write_report(stats_df, figure_paths, args.report)

    print(f"Saved sample-level table: {args.output}", flush=True)
    print(f"Saved summary table: {summary_path}", flush=True)
    print(f"Saved stats table: {args.stats_output}", flush=True)
    print(f"Saved v2 summary table: {args.summary_v2_output}", flush=True)
    print(f"Saved figures: {args.figure_dir}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
