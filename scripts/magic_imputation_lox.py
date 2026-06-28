#!/usr/bin/env python
from __future__ import annotations

import os
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYDEPS = PROJECT_ROOT / "local_pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / "local_cache" / "matplotlib"))

import anndata as ad
import magic
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu


INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
OUTPUT_CSV = PROJECT_ROOT / "results" / "magic_imputation_LOX.csv"
LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]


def matrix_to_gene_frame(matrix, var_names: pd.Index, genes: list[str]) -> pd.DataFrame:
    present = [gene for gene in genes if gene in var_names]
    indices = [var_names.get_loc(gene) for gene in present]
    values = matrix[:, indices]
    if sparse.issparse(values):
        values = values.toarray()
    else:
        values = np.asarray(values)
    return pd.DataFrame(values, columns=present)


def direction(delta: float) -> str:
    if not np.isfinite(delta) or np.isclose(delta, 0):
        return "flat"
    return "up" if delta > 0 else "down"


def main() -> int:
    print(f"Loading: {INPUT_H5AD}")
    adata = ad.read_h5ad(INPUT_H5AD)
    print(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print("Using normalized matrix: adata.X")

    obs = adata.obs[["stage", "cell_type", "cell_type_subset"]].copy()
    obs["stage"] = obs["stage"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)

    fb_mask = obs["cell_type"] == "FB"
    fb_obs = obs.loc[fb_mask, ["stage", "cell_type_subset"]].reset_index(drop=True)
    fb_x = adata.X[fb_mask.to_numpy(), :]
    var_names = pd.Index(adata.var_names.astype(str))

    print(f"Fibroblast subset: {fb_x.shape[0]:,} cells x {fb_x.shape[1]:,} genes")
    print("\nFB cell_type_subset x stage counts:")
    print(pd.crosstab(fb_obs["cell_type_subset"], fb_obs["stage"], margins=True, margins_name="Total").to_string())

    print("\nGene presence in adata.var_names:")
    for gene in LOX_GENES:
        print(f"  {gene}: {gene in var_names}")

    before_df = matrix_to_gene_frame(fb_x, var_names, LOX_GENES)
    missing_genes = [gene for gene in LOX_GENES if gene not in before_df.columns]
    if missing_genes:
        print(f"\nMissing genes skipped: {missing_genes}")

    # MAGIC accepts sparse inputs with gene names via a DataFrame; converting the full
    # normalized FB matrix to float32 dense is the most compatible path for this version.
    print("\nPreparing FB expression matrix for MAGIC...")
    if sparse.issparse(fb_x):
        fb_dense = fb_x.toarray().astype(np.float32, copy=False)
    else:
        fb_dense = np.asarray(fb_x, dtype=np.float32)
    fb_df = pd.DataFrame(fb_dense, columns=var_names)

    print("Running MAGIC: t=3, k=10, default knn_dist...")
    magic_operator = magic.MAGIC(knn=10, t=3, verbose=1)
    magic_operator.fit(fb_df)
    after_df = magic_operator.transform(genes=list(before_df.columns))
    if not isinstance(after_df, pd.DataFrame):
        after_df = pd.DataFrame(after_df, columns=before_df.columns)
    after_df = after_df.loc[:, list(before_df.columns)].reset_index(drop=True)

    rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []

    for gene in before_df.columns:
        for subtype in sorted(fb_obs["cell_type_subset"].unique()):
            subtype_mask = fb_obs["cell_type_subset"] == subtype
            stage_02_mask = subtype_mask & (fb_obs["stage"] == "02mo")
            stage_18_mask = subtype_mask & (fb_obs["stage"] == "18mo")

            before_02 = before_df.loc[stage_02_mask, gene].to_numpy()
            before_18 = before_df.loc[stage_18_mask, gene].to_numpy()
            after_02 = after_df.loc[stage_02_mask, gene].to_numpy()
            after_18 = after_df.loc[stage_18_mask, gene].to_numpy()

            mean_before_02 = float(np.mean(before_02)) if len(before_02) else np.nan
            mean_before_18 = float(np.mean(before_18)) if len(before_18) else np.nan
            mean_after_02 = float(np.mean(after_02)) if len(after_02) else np.nan
            mean_after_18 = float(np.mean(after_18)) if len(after_18) else np.nan

            before_delta = mean_before_18 - mean_before_02
            after_delta = mean_after_18 - mean_after_02
            before_direction = direction(before_delta)
            after_direction = direction(after_delta)
            direction_flip = (
                before_direction in {"up", "down"}
                and after_direction in {"up", "down"}
                and before_direction != after_direction
            )

            if len(after_02) and len(after_18):
                test = mannwhitneyu(after_18, after_02, alternative="two-sided")
                statistic = float(test.statistic)
                pvalue = float(test.pvalue)
                rank_biserial = float((2 * statistic / (len(after_18) * len(after_02))) - 1)
            else:
                statistic = np.nan
                pvalue = np.nan
                rank_biserial = np.nan

            comparison_rows.append(
                {
                    "gene": gene,
                    "cell_type_subset": subtype,
                    "n_02mo": int(len(after_02)),
                    "n_18mo": int(len(after_18)),
                    "mean_before_02mo": mean_before_02,
                    "mean_before_18mo": mean_before_18,
                    "mean_after_02mo": mean_after_02,
                    "mean_after_18mo": mean_after_18,
                    "before_direction": before_direction,
                    "after_direction": after_direction,
                    "direction_flip": direction_flip,
                    "mannwhitneyu_statistic": statistic,
                    "mannwhitneyu_pvalue": pvalue,
                    "rank_biserial_18mo_vs_02mo": rank_biserial,
                }
            )

            for stage, n_cells, mean_before, mean_after in [
                ("02mo", len(before_02), mean_before_02, mean_after_02),
                ("18mo", len(before_18), mean_before_18, mean_after_18),
            ]:
                rows.append(
                    {
                        "gene": gene,
                        "cell_type_subset": subtype,
                        "stage": stage,
                        "n_cells": int(n_cells),
                        "mean_before": mean_before,
                        "mean_after": mean_after,
                        "before_direction": before_direction,
                        "after_direction": after_direction,
                        "direction_flip": direction_flip,
                        "mannwhitneyu_statistic": statistic,
                        "mannwhitneyu_pvalue": pvalue,
                        "rank_biserial_18mo_vs_02mo": rank_biserial,
                    }
                )

    result_df = pd.DataFrame(rows)
    comparison_df = pd.DataFrame(comparison_rows)

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(OUTPUT_CSV, index=False)

    print("\nMAGIC LOX comparison table, one row per gene x FB subtype:")
    print(comparison_df.to_string(index=False))

    print("\nDirection flips:")
    flips = comparison_df.loc[comparison_df["direction_flip"]]
    if flips.empty:
        print("No direction flips detected.")
    else:
        print(flips.to_string(index=False))

    print(f"\nSaved CSV: {OUTPUT_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
