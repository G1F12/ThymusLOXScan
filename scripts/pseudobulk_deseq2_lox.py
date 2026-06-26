#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
OUTPUT_CSV = PROJECT_ROOT / "results" / "pseudobulk_deseq2_LOX.csv"
PARTIAL_DIR = PROJECT_ROOT / "results" / "pseudobulk_partial"
MIN_CELLS = 10
LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
PRIORITY_SUBTYPES = ["3:capsFB", "4:intFB", "5:medFB", "13:mTEC1"]


def sparse_row_sum(matrix, row_indices: np.ndarray) -> np.ndarray:
    """Sum selected rows from a sparse/dense expression matrix."""
    sub = matrix[row_indices, :]
    summed = np.asarray(sub.sum(axis=0)).ravel()
    return np.rint(summed).astype(np.int64)


def safe_filename(value: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in value)


def partial_csv_path(subtype: str) -> Path:
    return PARTIAL_DIR / f"{safe_filename(subtype)}.csv"


def get_subtype_order(adata: AnnData) -> list[str]:
    subtypes = sorted(adata.obs["cell_type_subset"].astype(str).unique())
    priority = [subtype for subtype in PRIORITY_SUBTYPES if subtype in subtypes]
    remaining = [subtype for subtype in subtypes if subtype not in priority]
    return priority + remaining


def build_pseudobulk_for_subtype(adata: AnnData, subtype: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if adata.raw is None:
        raise ValueError("adata.raw is missing; this workflow requires adata.raw.X counts.")

    obs = adata.obs[["sample", "cell_type_subset", "stage"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs["_row_pos"] = np.arange(adata.n_obs)
    obs = obs.loc[obs["cell_type_subset"] == subtype].copy()

    raw_x = adata.raw.X
    if not sparse.issparse(raw_x):
        raw_x = np.asarray(raw_x)
    genes = pd.Index(adata.raw.var_names.astype(str))

    count_rows: list[np.ndarray] = []
    meta_rows: list[dict[str, object]] = []

    grouped = obs.groupby("sample", sort=True, observed=True)
    for sample, group in grouped:
        n_cells = len(group)
        if n_cells < MIN_CELLS:
            continue

        stages = group["stage"].unique()
        if len(stages) != 1:
            raise ValueError(f"Expected one stage per sample; got {stages} for sample={sample}")

        sample_id = f"{sample}__{subtype}".replace(":", "_").replace("/", "_")
        row_sum = sparse_row_sum(raw_x, group["_row_pos"].to_numpy())
        count_rows.append(row_sum)
        meta_rows.append(
            {
                "sample_id": sample_id,
                "sample": sample,
                "cell_type_subset": subtype,
                "stage": stages[0],
                "n_cells": n_cells,
            }
        )

    counts = pd.DataFrame(count_rows, index=[r["sample_id"] for r in meta_rows], columns=genes)
    metadata = pd.DataFrame(meta_rows)
    if metadata.empty:
        metadata = pd.DataFrame(columns=["sample", "cell_type_subset", "stage", "n_cells"])
        metadata.index.name = "sample_id"
        return counts, metadata

    metadata = metadata.set_index("sample_id")
    metadata["stage"] = pd.Categorical(metadata["stage"], categories=["02mo", "18mo"], ordered=False)
    return counts, metadata


def empty_or_skipped_results(subtype: str, lox_genes_present: list[str], status: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "cell_type_subset": subtype,
                "gene": gene,
                "present": gene in lox_genes_present,
                "status": status,
                "n_pseudobulk_samples": 0,
                "n_02mo": 0,
                "n_18mo": 0,
                "cells_02mo": 0,
                "cells_18mo": 0,
                "baseMean": np.nan,
                "log2FoldChange": np.nan,
                "lfcSE": np.nan,
                "pvalue": np.nan,
                "padj": np.nan,
            }
            for gene in LOX_GENES
        ]
    )


def run_deseq_for_subtype(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    subtype: str,
    lox_genes_present: list[str],
) -> pd.DataFrame:
    if metadata.empty:
        return empty_or_skipped_results(subtype, lox_genes_present, "skipped_no_pseudobulk_samples")

    subtype_meta = metadata.loc[metadata["cell_type_subset"] == subtype].copy()
    subtype_counts = counts.loc[subtype_meta.index].copy()

    stage_counts = subtype_meta["stage"].value_counts()
    result_rows: list[dict[str, object]] = []
    if set(stage_counts.index[stage_counts > 0]) != {"02mo", "18mo"}:
        for gene in LOX_GENES:
            result_rows.append(
                {
                    "cell_type_subset": subtype,
                    "gene": gene,
                    "present": gene in lox_genes_present,
                    "status": "skipped_missing_stage",
                    "n_pseudobulk_samples": len(subtype_meta),
                    "n_02mo": int(stage_counts.get("02mo", 0)),
                    "n_18mo": int(stage_counts.get("18mo", 0)),
                    "cells_02mo": int(subtype_meta.loc[subtype_meta["stage"] == "02mo", "n_cells"].sum()),
                    "cells_18mo": int(subtype_meta.loc[subtype_meta["stage"] == "18mo", "n_cells"].sum()),
                    "baseMean": np.nan,
                    "log2FoldChange": np.nan,
                    "lfcSE": np.nan,
                    "pvalue": np.nan,
                    "padj": np.nan,
                }
            )
        return pd.DataFrame(result_rows)

    nonzero_genes = subtype_counts.sum(axis=0) > 0
    subtype_counts = subtype_counts.loc[:, nonzero_genes]

    dds = DeseqDataSet(
        counts=subtype_counts,
        metadata=subtype_meta[["stage"]],
        design="~ stage",
        ref_level=["stage", "02mo"],
        refit_cooks=False,
        quiet=True,
    )
    dds.deseq2()

    stats = DeseqStats(dds, contrast=["stage", "18mo", "02mo"], quiet=True)
    stats.summary()
    results = stats.results_df.copy()

    for gene in LOX_GENES:
        present = gene in lox_genes_present
        if gene in results.index:
            row = results.loc[gene]
            status = "ok"
            values = {
                "baseMean": row.get("baseMean", np.nan),
                "log2FoldChange": row.get("log2FoldChange", np.nan),
                "lfcSE": row.get("lfcSE", np.nan),
                "pvalue": row.get("pvalue", np.nan),
                "padj": row.get("padj", np.nan),
            }
        else:
            status = "not_tested_zero_counts" if present else "gene_absent"
            values = {
                "baseMean": np.nan,
                "log2FoldChange": np.nan,
                "lfcSE": np.nan,
                "pvalue": np.nan,
                "padj": np.nan,
            }

        result_rows.append(
            {
                "cell_type_subset": subtype,
                "gene": gene,
                "present": present,
                "status": status,
                "n_pseudobulk_samples": len(subtype_meta),
                "n_02mo": int(stage_counts.get("02mo", 0)),
                "n_18mo": int(stage_counts.get("18mo", 0)),
                "cells_02mo": int(subtype_meta.loc[subtype_meta["stage"] == "02mo", "n_cells"].sum()),
                "cells_18mo": int(subtype_meta.loc[subtype_meta["stage"] == "18mo", "n_cells"].sum()),
                **values,
            }
        )

    return pd.DataFrame(result_rows)


def failure_result(subtype: str, metadata: pd.DataFrame, lox_genes_present: list[str], exc: Exception) -> pd.DataFrame:
    subtype_meta = metadata.loc[metadata["cell_type_subset"] == subtype] if not metadata.empty else metadata
    stage_counts = subtype_meta["stage"].value_counts() if not subtype_meta.empty else pd.Series(dtype=int)
    return pd.DataFrame(
        [
            {
                "cell_type_subset": subtype,
                "gene": gene,
                "present": gene in lox_genes_present,
                "status": f"failed: {exc}",
                "n_pseudobulk_samples": len(subtype_meta),
                "n_02mo": int(stage_counts.get("02mo", 0)),
                "n_18mo": int(stage_counts.get("18mo", 0)),
                "cells_02mo": int(subtype_meta.loc[subtype_meta["stage"] == "02mo", "n_cells"].sum())
                if not subtype_meta.empty
                else 0,
                "cells_18mo": int(subtype_meta.loc[subtype_meta["stage"] == "18mo", "n_cells"].sum())
                if not subtype_meta.empty
                else 0,
                "baseMean": np.nan,
                "log2FoldChange": np.nan,
                "lfcSE": np.nan,
                "pvalue": np.nan,
                "padj": np.nan,
            }
            for gene in LOX_GENES
        ]
    )


def main() -> int:
    print(f"Loading: {INPUT_H5AD}", flush=True)
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.raw.n_vars:,} raw genes", flush=True)

    raw_genes = pd.Index(adata.raw.var_names.astype(str))
    lox_genes_present = [gene for gene in LOX_GENES if gene in raw_genes]
    print("\nLOX genes present in adata.raw.var_names:", flush=True)
    for gene in LOX_GENES:
        print(f"  {gene}: {gene in raw_genes}", flush=True)

    subtype_order = get_subtype_order(adata)
    print(f"\nFound {len(subtype_order)} cell_type_subset values.", flush=True)
    print("Subtype run order:", flush=True)
    for i, subtype in enumerate(subtype_order, start=1):
        print(f"  {i}: {subtype}", flush=True)

    PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
    partial_paths: list[Path] = []
    for i, subtype in enumerate(subtype_order, start=1):
        path = partial_csv_path(subtype)
        print(f"\nProcessing subtype {i} of {len(subtype_order)}: {subtype}", flush=True)
        if path.exists():
            print(f"  Existing partial found, reusing: {path}", flush=True)
            partial_paths.append(path)
            continue

        metadata = pd.DataFrame()
        try:
            counts, metadata = build_pseudobulk_for_subtype(adata, subtype)
            print(f"  Pseudobulk profiles with >= {MIN_CELLS} cells: {len(metadata)}", flush=True)
            if not metadata.empty:
                print(metadata.reset_index().to_string(index=False), flush=True)
            subtype_results = run_deseq_for_subtype(counts, metadata, subtype, lox_genes_present)
        except Exception as exc:
            print(f"  FAILED: {exc}", flush=True)
            subtype_results = failure_result(subtype, metadata, lox_genes_present, exc)

        subtype_results.to_csv(path, index=False)
        partial_paths.append(path)
        print(f"  Saved partial: {path}", flush=True)

    result_df = pd.concat([pd.read_csv(path) for path in partial_paths], ignore_index=True)
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(OUTPUT_CSV, index=False)

    print("\nSummary table: DESeq2 LOX-family results, contrast 18mo vs 02mo", flush=True)
    print(result_df.to_string(index=False), flush=True)
    print(f"\nSaved CSV: {OUTPUT_CSV}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
