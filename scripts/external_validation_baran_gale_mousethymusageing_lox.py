#!/usr/bin/env python
"""Baran-Gale / MouseThymusAgeing LOX-family TEC reanalysis."""

from __future__ import annotations

import argparse
import math
import textwrap
import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

try:
    import pyreadr
    import rdata
except ImportError:  # pragma: no cover - handled as feasibility at runtime.
    pyreadr = None
    rdata = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPORT_DIR = PROJECT_ROOT / "data" / "external" / "BaranGale"
FALLBACK_RDS_DIR = PROJECT_ROOT / "data" / "external" / "E-MTAB-8560" / "experimenthub_smartseq_rds"

TABLE_DIR = PROJECT_ROOT / "results" / "tables"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "external_validation" / "baran_gale"
REPORT_PATH = PROJECT_ROOT / "reports" / "baran_gale_mousethymusageing_lox_reanalysis.md"
FEASIBILITY_REPORT = PROJECT_ROOT / "reports" / "baran_gale_mousethymusageing_feasibility.md"

METADATA_SUMMARY_TABLE = TABLE_DIR / "external_baran_gale_metadata_summary.tsv"
BY_SAMPLE_TABLE = TABLE_DIR / "external_baran_gale_lox_by_sample.tsv"
AGE_SUMMARY_TABLE = TABLE_DIR / "external_baran_gale_lox_age_summary.tsv"
MODELS_TABLE = TABLE_DIR / "external_baran_gale_lox_models.tsv"
FEASIBILITY_TABLE = TABLE_DIR / "external_baran_gale_feasibility.tsv"

LOX_GENE_IDS = {
    "Lox": "ENSMUSG00000030084",
    "Loxl1": "ENSMUSG00000032383",
    "Loxl2": "ENSMUSG00000034205",
    "Loxl3": "ENSMUSG00000026922",
    "Loxl4": "ENSMUSG00000029723",
}
LOX_GENES = list(LOX_GENE_IDS)
AGE_ORDER = [1, 4, 16, 32, 52]
GROUP_ORDER = ["broad_TEC", "mTEC_like", "mTEClo", "mTEChi", "post_Aire_mTEC", "cTEC"]
GROUP_PALETTE = {
    "broad_TEC": "#4C78A8",
    "mTEC_like": "#D55E00",
    "mTEClo": "#CC79A7",
    "mTEChi": "#009E73",
    "post_Aire_mTEC": "#7F3C8D",
    "cTEC": "#F2B701",
}
CPM_SCALE = 1_000_000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", type=Path, default=EXPORT_DIR)
    parser.add_argument("--fallback-rds-dir", type=Path, default=FALLBACK_RDS_DIR)
    return parser.parse_args()


def write_feasibility(reason: str, detail: str) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FEASIBILITY_REPORT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "status": "metadata/data loading not possible",
                "reason": reason,
                "detail": detail,
            }
        ]
    ).to_csv(FEASIBILITY_TABLE, sep="\t", index=False)
    FEASIBILITY_REPORT.write_text(
        "\n".join(
            [
                "# Baran-Gale / MouseThymusAgeing feasibility",
                "",
                "Baran-Gale / MouseThymusAgeing remains a planned external TEC-focused analysis because the required expression/metadata could not be loaded in the current environment.",
                "",
                "## Reason",
                "",
                reason,
                "",
                "## Detail",
                "",
                detail,
            ]
        ),
        encoding="utf-8",
    )


def dframe_to_df(obj: Any) -> pd.DataFrame:
    data = {str(key): value for key, value in obj.listData.items()}
    df = pd.DataFrame(data)
    df.index = [str(x) for x in obj.rownames]
    return df


def read_r_dframe(path: Path) -> pd.DataFrame:
    if rdata is None:
        raise RuntimeError("Python package rdata is required to parse MouseThymusAgeing metadata RDS files.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parsed = rdata.parser.parse_file(path)
        obj = rdata.conversion.convert(parsed)
    return dframe_to_df(obj)


def read_r_dataframe(path: Path) -> pd.DataFrame:
    if pyreadr is None:
        raise RuntimeError("Python package pyreadr is required to parse MouseThymusAgeing count RDS files.")
    result = pyreadr.read_r(str(path))
    if not result:
        raise RuntimeError(f"No object found in RDS file: {path}")
    return next(iter(result.values()))


def parse_age_week(values: pd.Series) -> pd.Series:
    return values.astype(str).str.extract(r"(\d+)")[0].astype(int)


def load_from_r_export(export_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    values_path = export_dir / "baran_gale_lox_cell_values_export.tsv"
    metadata_path = export_dir / "baran_gale_metadata_summary_export.tsv"
    if not values_path.exists():
        raise FileNotFoundError(f"R export not found: {values_path}")
    values = pd.read_csv(values_path, sep="\t")
    metadata = pd.read_csv(metadata_path, sep="\t") if metadata_path.exists() else pd.DataFrame()
    values["age_week"] = parse_age_week(values["age"])
    values["sample_id"] = values["sample_id"].astype(str)
    values["sort_day"] = values["sort_day"].astype(str)
    values["sort_type"] = values["sort_type"].astype(str)
    values["subtype"] = values["subtype"].astype(str)
    values["cell_id"] = values["cell_id"].astype(str)
    values["raw_counts"] = values["raw_counts"].astype(float)
    values["raw_library_size"] = values["raw_library_size"].astype(float)
    values["detected"] = values["raw_counts"] > 0
    return values, metadata, "R export TSVs from MouseThymusAgeing::MouseSMARTseqData"


def load_from_cached_experimenthub(rds_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    if not rds_dir.exists():
        raise FileNotFoundError(f"Cached MouseThymusAgeing RDS directory not found: {rds_dir}")
    rowdata_path = rds_dir / "rowdata.rds"
    if not rowdata_path.exists():
        raise FileNotFoundError(f"Missing rowData RDS: {rowdata_path}")

    rowdata = read_r_dframe(rowdata_path)
    gene_to_row: dict[str, int] = {}
    for gene, ensembl_id in LOX_GENE_IDS.items():
        if "Geneid" in rowdata.columns and rowdata["Geneid"].astype(str).eq(ensembl_id).any():
            gene_to_row[gene] = int(rowdata.index[rowdata["Geneid"].astype(str).eq(ensembl_id)][0]) - 1
        elif gene in rowdata.index:
            gene_to_row[gene] = int(rowdata.index.get_loc(gene))
    missing = sorted(set(LOX_GENES) - set(gene_to_row))
    if missing:
        raise RuntimeError(f"LOX-family genes absent from rowData: {missing}")

    frames = []
    metadata_frames = []
    for day in range(1, 6):
        counts_path = rds_dir / f"counts-processed-day{day}.rds"
        coldata_path = rds_dir / f"coldata-day{day}.rds"
        if not counts_path.exists() or not coldata_path.exists():
            raise FileNotFoundError(f"Missing MouseThymusAgeing day {day} count or metadata RDS.")

        counts = read_r_dataframe(counts_path)
        coldata = read_r_dframe(coldata_path).copy()
        coldata["cell_id"] = coldata["CellID"].astype(str)
        coldata["age_week"] = parse_age_week(coldata["Age"])
        coldata["sample_id"] = "sort_day_" + coldata["SortDay"].astype(str)
        coldata["sort_day"] = coldata["SortDay"].astype(str)
        coldata["sort_type"] = coldata["SortType"].astype(str)
        coldata["subtype"] = coldata["SubType"].astype(str)
        library_size = counts.sum(axis=0).astype(float)
        library_size.index = library_size.index.astype(str)
        coldata["raw_library_size"] = coldata["cell_id"].map(library_size)
        metadata_frames.append(
            coldata[
                [
                    "cell_id",
                    "age_week",
                    "sample_id",
                    "sort_day",
                    "sort_type",
                    "subtype",
                    "raw_library_size",
                ]
            ]
        )

        for gene, row_idx in gene_to_row.items():
            raw = counts.iloc[row_idx].astype(float)
            df = pd.DataFrame({"cell_id": raw.index.astype(str), "gene": gene, "raw_counts": raw.to_numpy(float)})
            df = df.merge(
                coldata[
                    [
                        "cell_id",
                        "age_week",
                        "sample_id",
                        "sort_day",
                        "sort_type",
                        "subtype",
                        "raw_library_size",
                    ]
                ],
                on="cell_id",
                how="left",
                validate="one_to_one",
            )
            df["detected"] = df["raw_counts"] > 0
            frames.append(df)

    values = pd.concat(frames, ignore_index=True)
    metadata = pd.concat(metadata_frames, ignore_index=True).drop_duplicates("cell_id")
    metadata_summary = pd.DataFrame(
        [
            {"section": "dataset", "field": "source", "value": "MouseThymusAgeing ExperimentHub SMART-seq2 cached RDS"},
            {"section": "dataset", "field": "technology", "value": "SMART-seq2"},
            {"section": "dataset", "field": "n_cells", "value": str(metadata["cell_id"].nunique())},
            {"section": "metadata", "field": "age_weeks", "value": ";".join(map(str, sorted(metadata["age_week"].unique())))},
            {"section": "metadata", "field": "sample_proxy", "value": "SortDay; no clean donor-level replicate field was available in parsed colData"},
            {"section": "metadata", "field": "sort_type_labels", "value": ";".join(sorted(metadata["sort_type"].unique()))},
            {"section": "metadata", "field": "subtype_labels", "value": ";".join(sorted(metadata["subtype"].unique()))},
        ]
    )
    return values, metadata_summary, "cached MouseThymusAgeing ExperimentHub SMART-seq2 RDS files"


def load_values(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    try:
        return load_from_r_export(args.export_dir)
    except Exception as export_exc:
        try:
            values, metadata, source = load_from_cached_experimenthub(args.fallback_rds_dir)
            metadata = pd.concat(
                [
                    metadata,
                    pd.DataFrame(
                        [
                            {
                                "section": "loading",
                                "field": "r_export_status",
                                "value": f"not used; {type(export_exc).__name__}: {export_exc}",
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            return values, metadata, source
        except Exception as fallback_exc:
            raise RuntimeError(
                "Could not load R export TSVs or cached MouseThymusAgeing RDS files. "
                f"R export issue: {type(export_exc).__name__}: {export_exc}; "
                f"fallback issue: {type(fallback_exc).__name__}: {fallback_exc}"
            ) from fallback_exc


def add_analysis_groups(values: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subtype = values["subtype"].fillna("").astype(str)
    sort_type = values["sort_type"].fillna("").astype(str)
    group_masks = {
        "broad_TEC": pd.Series(True, index=values.index),
        "mTEC_like": sort_type.isin(["mTEClo", "mTEChi", "gmTEC"]) | subtype.str.contains("mTEC", case=False, na=False),
        "mTEClo": sort_type.eq("mTEClo"),
        "mTEChi": sort_type.eq("mTEChi"),
        "post_Aire_mTEC": subtype.str.contains("Aire", case=False, na=False) & subtype.str.contains("mTEC", case=False, na=False),
        "cTEC": sort_type.eq("cTEC") | subtype.str.contains("cTEC", case=False, na=False),
    }
    for group, mask in group_masks.items():
        df = values.loc[mask].copy()
        if df.empty:
            continue
        df["analysis_group"] = group
        rows.append(df)
    if not rows:
        raise RuntimeError("No TEC/mTEC-like analysis groups could be formed from metadata labels.")
    return pd.concat(rows, ignore_index=True)


def summarize_by_sample(grouped: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["analysis_group", "gene", "age_week", "sample_id"]
    for keys, df in grouped.groupby(group_cols, sort=True, observed=True):
        analysis_group, gene, age_week, sample_id = keys
        first = df.iloc[0]
        n_cells = int(df["cell_id"].nunique())
        raw_count_sum = float(df["raw_counts"].sum())
        library_sum = float(df.drop_duplicates("cell_id")["raw_library_size"].sum())
        cpm = raw_count_sum / library_sum * CPM_SCALE if library_sum > 0 else np.nan
        rows.append(
            {
                "analysis_group": analysis_group,
                "gene": gene,
                "age_week": int(age_week),
                "sample_id": sample_id,
                "sort_day": first["sort_day"],
                "sort_type_labels_in_group": ";".join(sorted(df["sort_type"].dropna().astype(str).unique())),
                "subtype_labels_in_group": ";".join(sorted(df["subtype"].dropna().astype(str).unique())),
                "n_cells": n_cells,
                "n_detecting_cells": int(df.loc[df["raw_counts"].gt(0), "cell_id"].nunique()),
                "detection_rate": float(df["raw_counts"].gt(0).mean()) if len(df) else np.nan,
                "summed_raw_counts": raw_count_sum,
                "raw_library_sum": library_sum,
                "CPM": cpm,
                "log2CPM1": math.log2(cpm + 1.0) if np.isfinite(cpm) else np.nan,
                "replicate_caution": "SortDay is used as a sample-level proxy; parsed metadata did not expose a clean donor-level replicate field.",
            }
        )
    return pd.DataFrame(rows)


def summarize_by_age(by_sample: pd.DataFrame) -> pd.DataFrame:
    return (
        by_sample.groupby(["analysis_group", "gene", "age_week"], sort=True, observed=True)
        .agg(
            n_sample_proxy=("sample_id", "nunique"),
            n_cells=("n_cells", "sum"),
            mean_detection_rate=("detection_rate", "mean"),
            median_detection_rate=("detection_rate", "median"),
            mean_log2CPM1=("log2CPM1", "mean"),
            median_log2CPM1=("log2CPM1", "median"),
        )
        .reset_index()
    )


def model_summaries(by_sample: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (group, gene), df in by_sample.groupby(["analysis_group", "gene"], sort=True, observed=True):
        df = df.sort_values(["age_week", "sample_id"])
        age_means = df.groupby("age_week", observed=True)[["detection_rate", "log2CPM1"]].mean()
        if age_means.empty:
            continue
        first_age = int(age_means.index.min())
        last_age = int(age_means.index.max())
        det_delta = float(age_means.loc[last_age, "detection_rate"] - age_means.loc[first_age, "detection_rate"])
        log_delta = float(age_means.loc[last_age, "log2CPM1"] - age_means.loc[first_age, "log2CPM1"])
        if df["age_week"].nunique() >= 3:
            det_rho, det_p = stats.spearmanr(df["age_week"], df["detection_rate"])
            log_rho, log_p = stats.spearmanr(df["age_week"], df["log2CPM1"])
        else:
            det_rho, det_p, log_rho, log_p = np.nan, np.nan, np.nan, np.nan
        early = df.loc[df["age_week"].le(4)]
        old = df.loc[df["age_week"].ge(32)]
        if len(early) >= 2 and len(old) >= 2:
            old_minus_young_log = float(old["log2CPM1"].mean() - early["log2CPM1"].mean())
            old_minus_young_det = float(old["detection_rate"].mean() - early["detection_rate"].mean())
            age_bin_p = float(stats.mannwhitneyu(old["log2CPM1"], early["log2CPM1"], alternative="two-sided").pvalue)
        else:
            old_minus_young_log, old_minus_young_det, age_bin_p = np.nan, np.nan, np.nan
        direction = "aged-lower" if log_delta < 0 else "aged-higher" if log_delta > 0 else "flat"
        rows.append(
            {
                "analysis_group": group,
                "gene": gene,
                "n_sample_proxy": int(df["sample_id"].nunique()),
                "n_age_points": int(df["age_week"].nunique()),
                "youngest_age_week": first_age,
                "oldest_age_week": last_age,
                "detection_delta_oldest_minus_youngest": det_delta,
                "log2CPM_delta_oldest_minus_youngest": log_delta,
                "detection_spearman_rho": float(det_rho),
                "detection_spearman_pvalue": float(det_p),
                "log2CPM_spearman_rho": float(log_rho),
                "log2CPM_spearman_pvalue": float(log_p),
                "old_32_52wk_minus_young_1_4wk_detection_delta": old_minus_young_det,
                "old_32_52wk_minus_young_1_4wk_log2CPM_delta": old_minus_young_log,
                "age_bin_log2CPM_mann_whitney_pvalue": age_bin_p,
                "direction_by_oldest_minus_youngest_log2CPM": direction,
                "model_caution": "Descriptive trend over sample-proxy summaries; donor-level replicate metadata was not available in parsed colData.",
            }
        )
    return pd.DataFrame(rows)


def classify(models: pd.DataFrame) -> tuple[str, str]:
    loxl2 = models.loc[models["gene"].eq("Loxl2") & models["analysis_group"].isin(["mTEC_like", "broad_TEC"])]
    if loxl2.empty:
        return "labels not compatible with TEC/mTEC-like comparison", "not available"
    mtec = loxl2.loc[loxl2["analysis_group"].eq("mTEC_like")]
    ref = mtec.iloc[0] if not mtec.empty else loxl2.iloc[0]
    delta = ref["log2CPM_delta_oldest_minus_youngest"]
    rho = ref["log2CPM_spearman_rho"]
    if pd.isna(delta):
        return "inconclusive due to insufficient biological replicate metadata", "not available"
    if delta < 0 and (pd.isna(rho) or rho <= 0):
        return "completed: Loxl2 aged-lower in TEC/mTEC-like group", "aged-lower"
    if delta > 0 and (pd.isna(rho) or rho >= 0):
        return "completed: Loxl2 aged-higher", "aged-higher"
    return "completed: no clear Loxl2 age trend", "mixed/no clear trend"


def save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    png = FIGURE_DIR / f"{stem}.png"
    pdf = FIGURE_DIR / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return [png, pdf]


def plot_loxl2(age_summary: pd.DataFrame) -> list[Path]:
    sns.set_theme(style="whitegrid", context="paper")
    paths: list[Path] = []
    df = age_summary.loc[age_summary["gene"].eq("Loxl2") & age_summary["analysis_group"].isin(GROUP_ORDER)].copy()
    df["analysis_group"] = pd.Categorical(df["analysis_group"], GROUP_ORDER, ordered=True)
    for metric, ylabel, stem in [
        ("mean_detection_rate", "Mean detection rate", "baran_gale_loxl2_detection_by_age"),
        ("mean_log2CPM1", "Mean log2(CPM + 1)", "baran_gale_loxl2_log2cpm_by_age"),
    ]:
        fig, ax = plt.subplots(figsize=(6.6, 4.2))
        sns.lineplot(
            data=df.sort_values(["analysis_group", "age_week"]),
            x="age_week",
            y=metric,
            hue="analysis_group",
            marker="o",
            palette=GROUP_PALETTE,
            ax=ax,
        )
        ax.set_xticks(AGE_ORDER)
        ax.set_xlabel("Age (weeks)")
        ax.set_ylabel(ylabel)
        ax.set_title("Baran-Gale / MouseThymusAgeing Loxl2 TEC/mTEC-like context")
        ax.legend(title="", frameon=False, fontsize=8)
        fig.tight_layout()
        paths.extend(save_figure(fig, stem))
    return paths


def fmt(value: Any, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def write_report(
    metadata: pd.DataFrame,
    by_sample: pd.DataFrame,
    age_summary: pd.DataFrame,
    models: pd.DataFrame,
    classification: str,
    key_direction: str,
    source: str,
    figure_paths: list[Path],
) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    loxl2 = models.loc[models["gene"].eq("Loxl2") & models["analysis_group"].isin(GROUP_ORDER)].copy()
    loxl2["analysis_group"] = pd.Categorical(loxl2["analysis_group"], GROUP_ORDER, ordered=True)
    loxl2 = loxl2.sort_values("analysis_group")

    if classification == "completed: Loxl2 aged-lower in TEC/mTEC-like group":
        interpretation = "Baran-Gale / MouseThymusAgeing provides independent TEC/mTEC-like directional support for aged-lower Loxl2 at the transcript level. This is not exact GSE240016 mTEC1 validation."
    elif classification == "completed: no clear Loxl2 age trend":
        interpretation = "No supportive aged-lower Loxl2 trend was observed in the Baran-Gale TEC-focused dataset under the tested annotation and sample-level summaries."
    elif classification == "completed: Loxl2 aged-higher":
        interpretation = "The tested Baran-Gale TEC/mTEC-like summaries showed aged-higher Loxl2, so this dataset does not provide aged-lower directional context under the current grouping."
    else:
        interpretation = "The Baran-Gale / MouseThymusAgeing reanalysis was inconclusive under the available metadata and grouping."

    metadata_lines = []
    if not metadata.empty:
        for _, row in metadata.iterrows():
            section = row.get("section", "metadata")
            field = row.get("field", row.iloc[0])
            value = row.get("value", row.iloc[-1])
            metadata_lines.append(f"- {section} / {field}: {value}")

    lines = [
        "# Baran-Gale / MouseThymusAgeing LOX-family reanalysis",
        "",
        "## Scope",
        "",
        "This reanalysis inspects LOX-family transcript abundance in public Baran-Gale / MouseThymusAgeing mouse thymic epithelial ageing data. It is framed as independent TEC/mTEC-like directional context, not exact subtype replication of GSE240016 mTEC1.",
        "",
        f"Input source used by this run: {source}.",
        "",
        "## Metadata inspection",
        "",
    ]
    lines.extend(metadata_lines or ["- Metadata summary was not available."])
    lines.extend(
        [
            "",
            "Candidate groups tested: broad TEC, mTEC-like, mTEClo, mTEChi, post-Aire mTEC when labeled, and cTEC.",
            "Sample-level summaries use SortDay as a sample proxy because parsed metadata did not expose a clean donor-level replicate field.",
            "",
            "## Loxl2 direction",
            "",
            "| group | sample proxies | age points | oldest-minus-youngest log2CPM delta | Spearman rho | direction |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    for _, row in loxl2.iterrows():
        lines.append(
            f"| {row['analysis_group']} | {int(row['n_sample_proxy'])} | {int(row['n_age_points'])} | "
            f"{fmt(row['log2CPM_delta_oldest_minus_youngest'])} | {fmt(row['log2CPM_spearman_rho'])} | "
            f"{row['direction_by_oldest_minus_youngest_log2CPM']} |"
        )
    lines.extend(
        [
            "",
            "## Classification",
            "",
            f"Classification: **{classification}**.",
            f"Key Loxl2 direction in TEC/mTEC-like groups: **{key_direction}**.",
            "",
            interpretation,
            "",
            "## Output files",
            "",
            f"- `{METADATA_SUMMARY_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
            f"- `{BY_SAMPLE_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
            f"- `{AGE_SUMMARY_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
            f"- `{MODELS_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
        ]
    )
    lines.extend(f"- `{path.relative_to(PROJECT_ROOT).as_posix()}`" for path in figure_paths)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    try:
        values, metadata, source = load_values(args)
        grouped = add_analysis_groups(values)
        by_sample = summarize_by_sample(grouped)
        age_summary = summarize_by_age(by_sample)
        models = model_summaries(by_sample)
        classification, key_direction = classify(models)
        figure_paths = plot_loxl2(age_summary)

        metadata.to_csv(METADATA_SUMMARY_TABLE, sep="\t", index=False)
        by_sample.to_csv(BY_SAMPLE_TABLE, sep="\t", index=False)
        age_summary.to_csv(AGE_SUMMARY_TABLE, sep="\t", index=False)
        models.to_csv(MODELS_TABLE, sep="\t", index=False)
        write_report(metadata, by_sample, age_summary, models, classification, key_direction, source, figure_paths)

        print("Baran-Gale analysis completed: yes", flush=True)
        print(f"Baran-Gale classification: {classification}", flush=True)
        print(f"Key Loxl2 direction: {key_direction}", flush=True)
        print(f"Saved report: {REPORT_PATH}", flush=True)
        return 0
    except Exception as exc:  # noqa: BLE001 - feasibility output should capture exact reason.
        reason = f"{type(exc).__name__}: {exc}"
        write_feasibility(
            "metadata/data loading not possible",
            textwrap.dedent(reason),
        )
        print("Baran-Gale analysis completed: no", flush=True)
        print(f"Baran-Gale classification: metadata/data loading not possible ({reason})", flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
