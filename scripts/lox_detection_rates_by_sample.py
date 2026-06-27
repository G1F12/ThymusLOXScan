#!/usr/bin/env python
"""Compute sample-level LOX-family detection rates in selected fibroblast/TEC groups."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
DEFAULT_OUTPUT_TABLE = PROJECT_ROOT / "results" / "tables" / "lox_detection_rates_by_sample.tsv"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "lox_detection_rate_analysis.md"

GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
YOUNG_LABEL = "02mo"
AGED_LABEL = "18mo"
GROUPS = [
    {"label": "pooled FB", "level": "cell_type", "column": "cell_type", "value": "FB"},
    {"label": "capsFB", "level": "subtype", "column": "cell_type_subset", "value": "3:capsFB"},
    {"label": "intFB", "level": "subtype", "column": "cell_type_subset", "value": "4:intFB"},
    {"label": "medFB", "level": "subtype", "column": "cell_type_subset", "value": "5:medFB"},
    {"label": "mTEC1", "level": "subtype", "column": "cell_type_subset", "value": "13:mTEC1"},
]
KEY_FINDINGS = [
    ("pooled FB", "Lox", "down"),
    ("pooled FB", "Loxl1", "down"),
    ("pooled FB", "Loxl2", "down"),
    ("capsFB", "Lox", "down"),
    ("medFB", "Loxl1", "up"),
    ("medFB", "Loxl2", "down"),
    ("mTEC1", "Loxl2", "down"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--output-table", type=Path, default=DEFAULT_OUTPUT_TABLE)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def vector_to_array(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def gene_vector(raw_x, gene_idx: int) -> np.ndarray:
    values = raw_x[:, gene_idx]
    return vector_to_array(values)


def compute_detection_table(adata) -> pd.DataFrame:
    if adata.raw is None:
        raise ValueError("Input AnnData must contain adata.raw.X raw counts.")

    obs = adata.obs[["sample", "stage", "cell_type", "cell_type_subset"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)

    raw_genes = pd.Index(adata.raw.var_names.astype(str))
    raw_x = adata.raw.X
    rows = []

    for group in GROUPS:
        group_mask = obs[group["column"]].eq(group["value"]).to_numpy()
        group_obs = obs.loc[group_mask].copy()
        group_positions = np.flatnonzero(group_mask)

        for gene in GENES:
            if gene not in raw_genes:
                continue
            expr = gene_vector(raw_x[group_positions, :], raw_genes.get_loc(gene))
            group_obs = group_obs.assign(_expr=expr)

            for sample_id, sample_obs in group_obs.groupby("sample", sort=True, observed=True):
                stages = sample_obs["stage"].unique()
                if len(stages) != 1:
                    raise ValueError(f"Sample {sample_id} has multiple age groups: {stages}")
                values = sample_obs["_expr"].to_numpy()
                n_cells = int(len(values))
                detecting = values > 0
                n_detecting = int(detecting.sum())
                raw_counts_sum = float(values.sum())
                rows.append(
                    {
                        "sample_id": sample_id,
                        "age_group": stages[0],
                        "cell_type_or_subtype": group["label"],
                        "annotation_level": group["level"],
                        "source_annotation": group["value"],
                        "gene": gene,
                        "n_cells": n_cells,
                        "raw_counts_sum": raw_counts_sum,
                        "n_detecting_cells": n_detecting,
                        "detection_rate": n_detecting / n_cells if n_cells else np.nan,
                        "mean_nonzero_expression": float(values[detecting].mean())
                        if n_detecting
                        else np.nan,
                    }
                )

    return pd.DataFrame(rows)


def summarize_age_change(table: pd.DataFrame, group: str, gene: str) -> dict[str, float | str]:
    subset = table.loc[table["cell_type_or_subtype"].eq(group) & table["gene"].eq(gene)].copy()
    means = subset.groupby("age_group", observed=True)[
        ["detection_rate", "mean_nonzero_expression", "n_cells"]
    ].mean()
    out: dict[str, float | str] = {"cell_type_or_subtype": group, "gene": gene}
    for metric in ["detection_rate", "mean_nonzero_expression", "n_cells"]:
        out[f"{metric}_02mo_mean"] = float(means.loc[YOUNG_LABEL, metric]) if YOUNG_LABEL in means.index else np.nan
        out[f"{metric}_18mo_mean"] = float(means.loc[AGED_LABEL, metric]) if AGED_LABEL in means.index else np.nan
        out[f"{metric}_delta_18mo_minus_02mo"] = (
            out[f"{metric}_18mo_mean"] - out[f"{metric}_02mo_mean"]
            if np.isfinite(out[f"{metric}_18mo_mean"]) and np.isfinite(out[f"{metric}_02mo_mean"])
            else np.nan
        )
    return out


def classify_mechanism(summary: dict[str, float | str], expected_direction: str) -> str:
    det_delta = float(summary["detection_rate_delta_18mo_minus_02mo"])
    nonzero_delta = float(summary["mean_nonzero_expression_delta_18mo_minus_02mo"])
    det_threshold = 0.03
    expr_threshold = 0.10

    if not np.isfinite(det_delta) or not np.isfinite(nonzero_delta):
        return "unclear"

    if expected_direction == "up":
        more_cells = det_delta >= det_threshold
        higher_nonzero = nonzero_delta >= expr_threshold
        if more_cells and higher_nonzero:
            return "unclear for fewer/lower framework; age effect is upward with both detection and nonzero expression higher"
        if more_cells:
            return "unclear for fewer/lower framework; detection increases while nonzero expression does not"
        if higher_nonzero:
            return "unclear for fewer/lower framework; nonzero expression increases without higher detection"
        return "unclear"

    fewer_cells = det_delta <= -det_threshold
    lower_nonzero = nonzero_delta <= -expr_threshold

    if fewer_cells and lower_nonzero:
        return "both fewer expressing cells and lower expression among expressing cells"
    if fewer_cells:
        return "fewer expressing cells"
    if lower_nonzero:
        return "lower expression among expressing cells"
    return "unclear"


def write_report(table: pd.DataFrame, report_path: Path) -> None:
    summaries = []
    for group, gene, expected_direction in KEY_FINDINGS:
        summary = summarize_age_change(table, group, gene)
        summary["expected_direction"] = expected_direction
        summaries.append(summary)

    lines = [
        "# LOX-family detection-rate analysis by sample",
        "",
        "## Scope",
        "",
        "Detection was computed per biological sample using raw counts in `adata.raw.X`. A cell is counted as detecting a gene when its raw count for that gene is greater than zero. Mean nonzero expression is the mean raw count among detecting cells only.",
        "",
        "This is descriptive and sample-level. It is intended to separate changes in the fraction of expressing cells from changes in expression intensity among cells that still express the gene.",
        "",
        "## Key finding interpretation",
        "",
        "| group | gene | detection 02mo | detection 18mo | nonzero mean 02mo | nonzero mean 18mo | interpretation |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]

    for summary in summaries:
        mechanism = classify_mechanism(summary, str(summary["expected_direction"]))
        lines.append(
            f"| {summary['cell_type_or_subtype']} | {summary['gene']} | "
            f"{summary['detection_rate_02mo_mean']:.3f} | {summary['detection_rate_18mo_mean']:.3f} | "
            f"{summary['mean_nonzero_expression_02mo_mean']:.3f} | "
            f"{summary['mean_nonzero_expression_18mo_mean']:.3f} | {mechanism} |"
        )

    lines.extend(
        [
            "",
            "## Notes by finding",
            "",
        ]
    )
    for summary in summaries:
        mechanism = classify_mechanism(summary, str(summary["expected_direction"]))
        det_delta = summary["detection_rate_delta_18mo_minus_02mo"]
        expr_delta = summary["mean_nonzero_expression_delta_18mo_minus_02mo"]
        lines.extend(
            [
                f"### {summary['gene']} in {summary['cell_type_or_subtype']}",
                "",
                f"- Detection-rate change, 18mo minus 02mo: {det_delta:.3f}.",
                f"- Mean-nonzero-expression change, 18mo minus 02mo: {expr_delta:.3f} raw counts among detecting cells.",
                f"- Expected pseudobulk direction: {summary['expected_direction']}_in_aged.",
                f"- Classification: {mechanism}.",
                "",
            ]
        )

    lines.extend(
        [
            "## Caveats",
            "",
            "Detection rates are sensitive to sequencing depth, library complexity, and the small number of biological samples. These summaries should be read alongside the pseudobulk age models rather than as standalone inference.",
        ]
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    print(f"Loading {args.input_h5ad}", flush=True)
    adata = ad.read_h5ad(args.input_h5ad)
    table = compute_detection_table(adata)

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_table, sep="\t", index=False)
    write_report(table, args.report)

    print(f"Saved table: {args.output_table}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
