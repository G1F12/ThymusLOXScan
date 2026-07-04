#!/usr/bin/env python
"""Exact label-permutation audit for GSE240016 mTEC1 Loxl2."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
TABLE_PATH = PROJECT_ROOT / "results" / "tables" / "mtec1_loxl2_permutation_test.tsv"
REPORT_PATH = PROJECT_ROOT / "reports" / "mtec1_loxl2_permutation_test.md"

GENE = "Loxl2"
SUBTYPE = "13:mTEC1"
MIN_CELLS = 10
CPM_SCALE = 1_000_000.0


def row_sums(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    return np.asarray(matrix).sum(axis=1)


def per_sample_values() -> pd.DataFrame:
    adata = sc.read_h5ad(INPUT_H5AD)
    if adata.raw is None:
        raise RuntimeError("adata.raw is required for raw-count pseudobulk summaries.")
    if GENE not in adata.raw.var_names:
        raise RuntimeError(f"{GENE} is absent from adata.raw.var_names.")

    obs = adata.obs[["sample", "stage", "cell_type_subset"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)
    obs["_row_pos"] = np.arange(adata.n_obs)
    obs = obs.loc[obs["cell_type_subset"].eq(SUBTYPE)].copy()

    gene_idx = int(pd.Index(adata.raw.var_names.astype(str)).get_loc(GENE))
    rows: list[dict[str, object]] = []
    for sample, sample_obs in obs.groupby("sample", sort=True, observed=True):
        if len(sample_obs) < MIN_CELLS:
            continue
        stages = sorted(sample_obs["stage"].unique())
        if len(stages) != 1:
            raise RuntimeError(f"Expected one stage for {sample}; observed {stages}.")
        row_idx = sample_obs["_row_pos"].to_numpy()
        raw_block = adata.raw.X[row_idx, :]
        gene_values = raw_block[:, gene_idx]
        gene_values = np.asarray(gene_values.toarray()).ravel() if sparse.issparse(gene_values) else np.asarray(gene_values).ravel()
        library_sum = float(row_sums(raw_block).sum())
        raw_count_sum = float(gene_values.sum())
        cpm = raw_count_sum / library_sum * CPM_SCALE if library_sum > 0 else math.nan
        rows.append(
            {
                "sample_id": sample,
                "age_group": stages[0],
                "n_mTEC1_cells": int(len(sample_obs)),
                "Loxl2_raw_count_sum": raw_count_sum,
                "Loxl2_detecting_cells": int((gene_values > 0).sum()),
                "Loxl2_detection_rate": float((gene_values > 0).mean()),
                "raw_library_sum": library_sum,
                "Loxl2_CPM": cpm,
                "Loxl2_log2CPM1": float(math.log2(cpm + 1.0)) if math.isfinite(cpm) else math.nan,
            }
        )
    values = pd.DataFrame(rows).sort_values(["age_group", "sample_id"]).reset_index(drop=True)
    stage_counts = values["age_group"].value_counts().to_dict()
    if stage_counts.get("02mo", 0) != 2 or stage_counts.get("18mo", 0) != 2:
        raise RuntimeError(f"Expected 2 young and 2 aged samples after filtering; observed {stage_counts}.")
    return values


def stat(values: pd.DataFrame, aged_samples: set[str], metric: str) -> float:
    aged = values.loc[values["sample_id"].isin(aged_samples), metric].mean()
    young = values.loc[~values["sample_id"].isin(aged_samples), metric].mean()
    return float(aged - young)


def permutation_table(values: pd.DataFrame) -> pd.DataFrame:
    observed_aged = set(values.loc[values["age_group"].eq("18mo"), "sample_id"])
    samples = list(values["sample_id"])
    rows: list[dict[str, object]] = []
    for metric, label in [
        ("Loxl2_detection_rate", "detection_rate"),
        ("Loxl2_log2CPM1", "log2CPM1"),
    ]:
        observed = stat(values, observed_aged, metric)
        perm_stats = []
        for i, aged_tuple in enumerate(itertools.combinations(samples, 2), start=1):
            aged = set(aged_tuple)
            current = stat(values, aged, metric)
            perm_stats.append(current)
            rows.append(
                {
                    "metric": label,
                    "statistic_direction": "aged_minus_young",
                    "permutation_id": i,
                    "assigned_aged_samples": ";".join(sorted(aged)),
                    "assigned_young_samples": ";".join(sorted(set(samples) - aged)),
                    "permuted_statistic": current,
                    "observed_statistic": observed,
                    "is_observed_labeling": aged == observed_aged,
                }
            )
        perm = np.asarray(perm_stats)
        two_sided = float(np.mean(np.abs(perm) >= abs(observed) - 1e-12))
        one_sided = float(np.mean(perm <= observed + 1e-12))
        for row in rows:
            if row["metric"] == label:
                row["two_sided_exact_pvalue"] = two_sided
                row["one_sided_exact_pvalue_aged_lower"] = one_sided
                row["minimum_one_sided_pvalue"] = 1.0 / len(perm)
                row["minimum_two_sided_pvalue"] = 2.0 / len(perm)
    return pd.DataFrame(rows)


def write_report(values: pd.DataFrame, tests: pd.DataFrame) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# mTEC1 Loxl2 exact permutation test",
        "",
        "## Scope",
        "",
        "This audit uses the four GSE240016 biological samples with at least 10 annotated `13:mTEC1` cells. Loxl2 values were recomputed from `adata.raw.X` in `data/raw/GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad`.",
        "",
        "With n=2 versus n=2 biological samples, the exact permutation test has only six possible labelings, so p-value resolution is limited. This test is used to calibrate the parametric DESeq2 result rather than replace it.",
        "",
        "## Sample values",
        "",
        "| sample | age | cells | detecting cells | detection rate | log2(CPM+1) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in values.iterrows():
        lines.append(
            f"| {row['sample_id']} | {row['age_group']} | {int(row['n_mTEC1_cells'])} | "
            f"{int(row['Loxl2_detecting_cells'])} | {row['Loxl2_detection_rate']:.6f} | "
            f"{row['Loxl2_log2CPM1']:.6f} |"
        )
    lines.extend(["", "## Exact permutation results", ""])
    for metric in ["detection_rate", "log2CPM1"]:
        one = tests.loc[tests["metric"].eq(metric)].iloc[0]
        lines.extend(
            [
                f"- {metric}: observed aged-minus-young statistic = {one['observed_statistic']:.6f}.",
                f"- {metric}: two-sided exact p-value = {one['two_sided_exact_pvalue']:.6f}.",
                f"- {metric}: one-sided exact p-value for aged-lower direction = {one['one_sided_exact_pvalue_aged_lower']:.6f}.",
            ]
        )
    lines.extend(
        [
            "",
            "The minimum possible one-sided p-value is 1/6 = 0.166667; the minimum possible two-sided p-value is 2/6 = 0.333333 when both directions are counted.",
            "",
            "## Interpretation",
            "",
            "Both young samples remain higher than both aged samples for detection rate and log2(CPM+1). The exact test supports the direction as sample-level ordered, but the p-values cannot be small with only four samples.",
            "",
            "## Output",
            "",
            f"- `{TABLE_PATH.relative_to(PROJECT_ROOT).as_posix()}`",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    values = per_sample_values()
    tests = permutation_table(values)
    tests.to_csv(TABLE_PATH, sep="\t", index=False)
    write_report(values, tests)
    print(f"Saved {TABLE_PATH}")
    print(f"Saved {REPORT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
