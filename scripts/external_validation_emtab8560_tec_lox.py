#!/usr/bin/env python
"""External TEC age-series validation for LOX-family genes in E-MTAB-8560.

This script uses the public ArrayExpress/BioStudies SDRF plus the Bioconductor
MouseThymusAgeing ExperimentHub SMART-seq2 resources. Cell-level values are
summarized descriptively because cells are not independent biological
replicates.
"""

from __future__ import annotations

import argparse
import sqlite3
import textwrap
import urllib.request
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

try:
    import pyreadr
    import rdata
except ImportError:  # pragma: no cover - handled by failure report at runtime.
    pyreadr = None
    rdata = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data" / "external" / "E-MTAB-8560"
DEFAULT_BY_AGE = PROJECT_ROOT / "results" / "tables" / "external_emtab8560_tec_lox_by_age.tsv"
DEFAULT_SUMMARY = PROJECT_ROOT / "results" / "tables" / "external_emtab8560_tec_lox_summary.tsv"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "external_validation" / "emtab8560"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "emtab8560_tec_external_validation.md"
DEFAULT_FAILURE_REPORT = PROJECT_ROOT / "reports" / "emtab8560_download_or_parse_failure.md"

ARRAYEXPRESS_BASE = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8560"
EXPERIMENTHUB_SQLITE = "https://experimenthub.bioconductor.org/metadata/experimenthub.sqlite3"
SMARTSEQ_PATH_PREFIX = "MouseThymusAgeing/SMARTseq/1.0.0"

LOX_GENE_IDS = {
    "Lox": "ENSMUSG00000030084",
    "Loxl1": "ENSMUSG00000032383",
    "Loxl2": "ENSMUSG00000034205",
    "Loxl3": "ENSMUSG00000026922",
    "Loxl4": "ENSMUSG00000029723",
}
AGE_ORDER = [1, 4, 16, 32, 52]
AGE_LABELS = {age: f"{age} wk" for age in AGE_ORDER}
GROUP_ORDER = ["broad_TEC", "cTEC", "mTEClo", "mTEChi", "mTEC_like"]
PLOT_PALETTE = {
    "broad_TEC": "#0072B2",
    "cTEC": "#009E73",
    "mTEClo": "#CC79A7",
    "mTEChi": "#D55E00",
    "mTEC_like": "#6A3D9A",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--by-age-output", type=Path, default=DEFAULT_BY_AGE)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--failure-report", type=Path, default=DEFAULT_FAILURE_REPORT)
    return parser.parse_args()


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        print(f"Downloading {url}", flush=True)
        urllib.request.urlretrieve(url, dest)
    return dest


def write_failure_report(report_path: Path, reason: str, steps: list[str]) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# E-MTAB-8560 download or parse failure",
        "",
        "## Exact reason",
        "",
        reason,
        "",
        "## Next manual steps",
        "",
    ]
    lines.extend(f"{i}. {step}" for i, step in enumerate(steps, start=1))
    report_path.write_text("\n".join(lines), encoding="utf-8")


def dframe_to_df(obj) -> pd.DataFrame:
    data = {str(key): value for key, value in obj.listData.items()}
    df = pd.DataFrame(data)
    df.index = [str(x) for x in obj.rownames]
    return df


def read_r_dframe(path: Path) -> pd.DataFrame:
    if rdata is None:
        raise RuntimeError("Python package 'rdata' is required to parse Bioconductor DFrame RDS files.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parsed = rdata.parser.parse_file(path)
        obj = rdata.conversion.convert(parsed)
    return dframe_to_df(obj)


def read_r_dataframe(path: Path) -> pd.DataFrame:
    if pyreadr is None:
        raise RuntimeError("Python package 'pyreadr' is required to parse count matrix RDS files.")
    result = pyreadr.read_r(str(path))
    if not result:
        raise RuntimeError(f"No object found in RDS file: {path}")
    return next(iter(result.values()))


def experimenthub_resources(sqlite_path: Path) -> pd.DataFrame:
    con = sqlite3.connect(sqlite_path)
    try:
        query = """
        select resources.ah_id, resources.title, location_prefixes.location_prefix,
               rdatapaths.rdatapath
        from resources
        join rdatapaths on resources.id = rdatapaths.resource_id
        join location_prefixes on resources.location_prefix_id = location_prefixes.id
        where rdatapaths.rdatapath like ?
        order by resources.ah_id
        """
        return pd.read_sql_query(query, con, params=(f"%{SMARTSEQ_PATH_PREFIX}%",))
    finally:
        con.close()


def download_public_inputs(data_dir: Path) -> tuple[Path, pd.DataFrame]:
    for name in ["E-MTAB-8560.idf.txt", "E-MTAB-8560.sdrf.txt"]:
        download(f"{ARRAYEXPRESS_BASE}/{name}", data_dir / name)

    sqlite_path = download(EXPERIMENTHUB_SQLITE, data_dir / "experimenthub.sqlite3")
    resources = experimenthub_resources(sqlite_path)
    if resources.empty:
        raise RuntimeError("ExperimentHub metadata did not contain MouseThymusAgeing SMARTseq resources.")

    rds_dir = data_dir / "experimenthub_smartseq_rds"
    wanted = resources.loc[
        resources["rdatapath"].str.contains("rowdata|coldata-day|counts-processed-day", regex=True)
    ].copy()
    for _, row in wanted.iterrows():
        url = row["location_prefix"] + row["rdatapath"]
        download(url, rds_dir / Path(row["rdatapath"]).name)
    return rds_dir, wanted


def load_smartseq(data_dir: Path) -> pd.DataFrame:
    rds_dir, _ = download_public_inputs(data_dir)
    rowdata = read_r_dframe(rds_dir / "rowdata.rds")
    gene_to_row = {
        gene: int(rowdata.index[rowdata["Geneid"].eq(ensembl_id)][0]) - 1
        for gene, ensembl_id in LOX_GENE_IDS.items()
        if rowdata["Geneid"].eq(ensembl_id).any()
    }
    missing = sorted(set(LOX_GENE_IDS) - set(gene_to_row))
    if missing:
        raise RuntimeError(f"LOX-family Ensembl IDs absent from rowData: {missing}")

    frames = []
    for day in range(1, 6):
        counts = read_r_dataframe(rds_dir / f"counts-processed-day{day}.rds")
        coldata = read_r_dframe(rds_dir / f"coldata-day{day}.rds")
        coldata = coldata.copy()
        coldata["cell_id"] = coldata["CellID"].astype(str)
        coldata["age_week"] = coldata["Age"].astype(str).str.extract(r"(\d+)").astype(int)
        coldata["sort_day"] = coldata["SortDay"].astype(int)
        coldata["sort_type"] = coldata["SortType"].astype(str)
        coldata["subtype"] = coldata["SubType"].astype(str)
        coldata["size_factor"] = coldata["sizeFactor"].astype(float)

        for gene, row_idx in gene_to_row.items():
            raw = counts.iloc[row_idx].astype(float)
            values = pd.DataFrame(
                {
                    "cell_id": raw.index.astype(str),
                    "gene": gene,
                    "raw_counts": raw.to_numpy(float),
                }
            ).merge(
                coldata[
                    [
                        "cell_id",
                        "age_week",
                        "sort_day",
                        "sort_type",
                        "subtype",
                        "size_factor",
                    ]
                ],
                on="cell_id",
                how="left",
                validate="one_to_one",
            )
            values["normalized_expression"] = values["raw_counts"] / values["size_factor"]
            values["log2_norm_plus1"] = np.log2(values["normalized_expression"] + 1)
            values["detected"] = values["raw_counts"] > 0
            frames.append(values)
    return pd.concat(frames, ignore_index=True)


def add_group_rows(values: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_masks = {
        "broad_TEC": pd.Series(True, index=values.index),
        "cTEC": values["sort_type"].eq("cTEC") | values["subtype"].str.contains("cTEC", case=False, na=False),
        "mTEClo": values["sort_type"].eq("mTEClo"),
        "mTEChi": values["sort_type"].eq("mTEChi"),
        "mTEC_like": values["sort_type"].isin(["mTEClo", "mTEChi", "gmTEC"])
        | values["subtype"].str.contains("mTEC", case=False, na=False),
    }
    for group, mask in group_masks.items():
        df = values.loc[mask].copy()
        df["analysis_group"] = group
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def summarize_by_age(grouped_values: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_level = (
        grouped_values.groupby(["analysis_group", "gene", "age_week", "sort_day"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_log2_norm_plus1=("log2_norm_plus1", "mean"),
            detection_rate=("detected", "mean"),
        )
        .reset_index()
    )

    by_age = (
        grouped_values.groupby(["analysis_group", "gene", "age_week"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            n_sort_day_samples=("sort_day", "nunique"),
            mean_log2_norm_plus1=("log2_norm_plus1", "mean"),
            median_log2_norm_plus1=("log2_norm_plus1", "median"),
            detection_rate=("detected", "mean"),
        )
        .reset_index()
    )
    by_age["n_biological_samples_if_available"] = pd.NA
    by_age["age_label"] = by_age["age_week"].map(AGE_LABELS)
    by_age["replicate_definition"] = (
        "sort_day/acquisition-day summaries; exact donor-level metadata is not available in parsed colData"
    )

    summary_rows = []
    for (group, gene), df in sample_level.groupby(["analysis_group", "gene"], observed=True):
        df = df.sort_values("age_week")
        age_means = df.groupby("age_week", observed=True)["mean_log2_norm_plus1"].mean()
        first_age = int(age_means.index.min())
        last_age = int(age_means.index.max())
        delta_old_minus_young = float(age_means.loc[last_age] - age_means.loc[first_age])
        direction = (
            "aged-lower"
            if delta_old_minus_young < 0
            else "aged-higher"
            if delta_old_minus_young > 0
            else "flat"
        )
        if df["age_week"].nunique() >= 3:
            rho, spearman_p = stats.spearmanr(df["age_week"], df["mean_log2_norm_plus1"])
        else:
            rho, spearman_p = np.nan, np.nan

        early = df.loc[df["age_week"].le(4), "mean_log2_norm_plus1"].dropna()
        old = df.loc[df["age_week"].ge(32), "mean_log2_norm_plus1"].dropna()
        if len(early) >= 2 and len(old) >= 2:
            bin_delta = float(old.mean() - early.mean())
            bin_p = float(stats.mannwhitneyu(old, early, alternative="two-sided").pvalue)
        else:
            bin_delta, bin_p = np.nan, np.nan

        summary_rows.append(
            {
                "analysis_group": group,
                "gene": gene,
                "n_cells_total": int(
                    grouped_values.loc[
                        grouped_values["analysis_group"].eq(group)
                        & grouped_values["gene"].eq(gene),
                        "cell_id",
                    ].nunique()
                ),
                "n_biological_samples_if_available": pd.NA,
                "n_sort_day_samples_available": int(df["sort_day"].nunique()),
                "youngest_age_week": first_age,
                "oldest_age_week": last_age,
                "mean_delta_oldest_minus_youngest": delta_old_minus_young,
                "age_trend_direction": direction,
                "spearman_rho_age_vs_sample_mean": float(rho),
                "spearman_pvalue": float(spearman_p),
                "old_32_52wk_minus_young_1_4wk_delta": bin_delta,
                "old_vs_young_bin_mann_whitney_pvalue": bin_p,
                "replicate_caution": (
                    "Cells are not independent biological replicates; trend tests use sort_day x age summaries."
                ),
            }
        )
    return by_age, pd.DataFrame(summary_rows)


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_expression(by_age: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    df = by_age.loc[by_age["analysis_group"].isin(["broad_TEC", "cTEC", "mTEC_like"])].copy()
    grid = sns.relplot(
        data=df,
        x="age_week",
        y="mean_log2_norm_plus1",
        hue="analysis_group",
        col="gene",
        col_wrap=3,
        kind="line",
        marker="o",
        palette=PLOT_PALETTE,
        height=3.0,
        aspect=1.05,
        facet_kws={"sharey": False},
    )
    grid.set_axis_labels("Age (weeks)", "Mean log2(size-factor normalized counts + 1)")
    grid.set_titles("{col_name}")
    for ax in grid.axes.flat:
        ax.set_xticks(AGE_ORDER)
    grid.fig.suptitle("E-MTAB-8560 LOX-family expression by TEC age series", y=1.03)
    return save_figure(grid.fig, output_dir, "emtab8560_lox_family_expression_by_age")


def plot_loxl2_focus(by_age: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    df = by_age.loc[
        by_age["gene"].eq("Loxl2") & by_age["analysis_group"].isin(["broad_TEC", "mTEC_like", "mTEClo", "mTEChi"])
    ].copy()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    sns.lineplot(
        data=df,
        x="age_week",
        y="mean_log2_norm_plus1",
        hue="analysis_group",
        marker="o",
        palette=PLOT_PALETTE,
        ax=ax,
    )
    ax.set_xticks(AGE_ORDER)
    ax.set_xlabel("Age (weeks)")
    ax.set_ylabel("Mean log2(size-factor normalized counts + 1)")
    ax.set_title("E-MTAB-8560 Loxl2 in broad TEC and mTEC-like groups")
    ax.legend(title="", frameon=False)
    fig.tight_layout()
    return save_figure(fig, output_dir, "emtab8560_loxl2_broad_tec_mtec_like_by_age")


def plot_detection(by_age: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    sns.set_theme(style="whitegrid", context="paper")
    df = by_age.loc[by_age["analysis_group"].isin(["broad_TEC", "cTEC", "mTEC_like"])].copy()
    grid = sns.relplot(
        data=df,
        x="age_week",
        y="detection_rate",
        hue="analysis_group",
        col="gene",
        col_wrap=3,
        kind="line",
        marker="o",
        palette=PLOT_PALETTE,
        height=3.0,
        aspect=1.05,
    )
    grid.set_axis_labels("Age (weeks)", "Detection rate")
    grid.set_titles("{col_name}")
    for ax in grid.axes.flat:
        ax.set_xticks(AGE_ORDER)
        ax.set_ylim(-0.02, 1.02)
    grid.fig.suptitle("E-MTAB-8560 LOX-family detection rate by age", y=1.03)
    return save_figure(grid.fig, output_dir, "emtab8560_lox_family_detection_by_age")


def fmt(value: float, digits: int = 3) -> str:
    return "NA" if pd.isna(value) else f"{value:.{digits}f}"


def write_report(summary: pd.DataFrame, figure_paths: list[Path], report_path: Path) -> None:
    def get(group: str, gene: str) -> pd.Series:
        hit = summary.loc[summary["analysis_group"].eq(group) & summary["gene"].eq(gene)]
        if hit.empty:
            raise KeyError(f"Missing summary row for {group} {gene}")
        return hit.iloc[0]

    broad = get("broad_TEC", "Loxl2")
    mtec = get("mTEC_like", "Loxl2")
    mteclo = get("mTEClo", "Loxl2")
    mtechi = get("mTEChi", "Loxl2")

    broad_direction = broad["age_trend_direction"]
    mtec_direction = mtec["age_trend_direction"]
    if broad_direction == "aged-lower" or mtec_direction == "aged-lower":
        support = (
            "supports a broad aged-lower TEC/mTEC-like tendency, but cannot specifically validate "
            "the current mTEC1 Loxl2 finding."
        )
    elif broad_direction == "aged-higher" and mtec_direction == "aged-higher":
        support = (
            "contradicts a broad aged-lower tendency in this external TEC age series, while still "
            "not directly testing mTEC1."
        )
    else:
        support = "cannot establish a clear broad TEC or mTEC-like aged-lower tendency."

    lines = [
        "# E-MTAB-8560 TEC external validation",
        "",
        "## Scope",
        "",
        "E-MTAB-8560 is a mouse TEC Smart-seq2 single-cell age-series spanning 1, 4, 16, 32, and 52 weeks. This script uses public ArrayExpress/BioStudies metadata plus processed Bioconductor MouseThymusAgeing SMART-seq2 resources.",
        "",
        "Cells are summarized descriptively. The parsed metadata provides cell IDs, age, sort type, sort day, and subtype labels, but not a clean donor-level design suitable for treating cells as independent biological replicates. Spearman and age-bin comparisons therefore use sort_day x age summaries and should be read as trend checks only.",
        "",
        "## Loxl2 interpretation",
        "",
        "| group | cells | sort-day samples | oldest-minus-youngest mean delta | Spearman rho | direction |",
        "|---|---:|---:|---:|---:|---|",
        f"| broad TEC | {int(broad['n_cells_total'])} | {int(broad['n_sort_day_samples_available'])} | {fmt(broad['mean_delta_oldest_minus_youngest'])} | {fmt(broad['spearman_rho_age_vs_sample_mean'])} | {broad['age_trend_direction']} |",
        f"| mTEC-like | {int(mtec['n_cells_total'])} | {int(mtec['n_sort_day_samples_available'])} | {fmt(mtec['mean_delta_oldest_minus_youngest'])} | {fmt(mtec['spearman_rho_age_vs_sample_mean'])} | {mtec['age_trend_direction']} |",
        f"| mTEClo | {int(mteclo['n_cells_total'])} | {int(mteclo['n_sort_day_samples_available'])} | {fmt(mteclo['mean_delta_oldest_minus_youngest'])} | {fmt(mteclo['spearman_rho_age_vs_sample_mean'])} | {mteclo['age_trend_direction']} |",
        f"| mTEChi | {int(mtechi['n_cells_total'])} | {int(mtechi['n_sort_day_samples_available'])} | {fmt(mtechi['mean_delta_oldest_minus_youngest'])} | {fmt(mtechi['spearman_rho_age_vs_sample_mean'])} | {mtechi['age_trend_direction']} |",
        "",
        f"Overall, E-MTAB-8560 {support}",
        "",
        "## Relationship to current mTEC1 Loxl2",
        "",
        "- This dataset can test broad TEC, cTEC, mTEClo, mTEChi, and mTEC-like age trends.",
        "- It cannot directly test the current mTEC1 label because the annotation systems do not match one-to-one.",
        "- The age range ends at 52 weeks, whereas the current aged comparison is older; this may capture maturation and midlife remodeling more than late thymic involution.",
        "- Any agreement should be framed as directional external context, not subtype-specific validation.",
        "",
        "## Figures",
        "",
    ]
    lines.extend(f"- `{path.relative_to(PROJECT_ROOT).as_posix()}`" for path in figure_paths)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    try:
        values = load_smartseq(args.data_dir)
        grouped_values = add_group_rows(values)
        by_age, summary = summarize_by_age(grouped_values)

        args.by_age_output.parent.mkdir(parents=True, exist_ok=True)
        by_age.to_csv(args.by_age_output, sep="\t", index=False)
        summary.to_csv(args.summary_output, sep="\t", index=False)

        figure_paths = []
        for paths in [
            plot_expression(by_age, args.figure_dir),
            plot_loxl2_focus(by_age, args.figure_dir),
            plot_detection(by_age, args.figure_dir),
        ]:
            figure_paths.extend(paths)
        write_report(summary, figure_paths, args.report)

        print(f"Saved by-age table: {args.by_age_output}", flush=True)
        print(f"Saved summary table: {args.summary_output}", flush=True)
        print(f"Saved figures: {args.figure_dir}", flush=True)
        print(f"Saved report: {args.report}", flush=True)
        return 0
    except Exception as exc:  # noqa: BLE001 - failure report should catch exact parse/download issue.
        reason = f"{type(exc).__name__}: {exc}"
        steps = [
            "Confirm that ArrayExpress SDRF/IDF files are reachable from https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8560/.",
            "Confirm that ExperimentHub metadata contains MouseThymusAgeing/SMARTseq/1.0.0 resources.",
            "Install Python packages pyreadr and rdata, or use R with MouseThymusAgeing::MouseSMARTseqData() to export counts, rowData, and colData.",
            "If Bioconductor object formats changed, manually export the SMART-seq2 SingleCellExperiment to TSV or H5AD and rerun this analysis from those exported files.",
        ]
        write_failure_report(args.failure_report, textwrap.dedent(reason), steps)
        print(f"Failed: {reason}", flush=True)
        print(f"Saved failure report: {args.failure_report}", flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
