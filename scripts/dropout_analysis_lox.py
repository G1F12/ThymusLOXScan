#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYDEPS = PROJECT_ROOT / "local_pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact


INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
OUTPUT_CSV = PROJECT_ROOT / "results" / "dropout_analysis_LOX.csv"
LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    pvals = np.asarray(pvalues, dtype=float)
    padj = np.full(pvals.shape, np.nan, dtype=float)
    valid = np.isfinite(pvals)
    if not valid.any():
        return padj.tolist()

    valid_p = pvals[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    valid_indices = np.flatnonzero(valid)
    padj[valid_indices[order]] = adjusted
    return padj.tolist()


def get_gene_vector(raw_x, gene_idx: int) -> np.ndarray:
    values = raw_x[:, gene_idx]
    if sparse.issparse(values):
        values = values.toarray().ravel()
    else:
        values = np.asarray(values).ravel()
    return values


def main() -> int:
    print(f"Loading: {INPUT_H5AD}")
    adata = ad.read_h5ad(INPUT_H5AD)
    if adata.raw is None:
        raise ValueError("adata.raw is missing; this analysis requires adata.raw.X.")

    obs = adata.obs[["cell_type_subset", "stage"]].copy()
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)
    obs["stage"] = obs["stage"].astype(str)

    raw_genes = pd.Index(adata.raw.var_names.astype(str))
    raw_x = adata.raw.X

    print(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.raw.n_vars:,} raw genes")
    print("Using raw count matrix: adata.raw.X")
    print("\nGene presence:")
    for gene in LOX_GENES:
        print(f"  {gene}: {gene in raw_genes}")

    summary_rows: list[dict[str, object]] = []
    test_rows: list[dict[str, object]] = []

    group_order = (
        obs[["cell_type_subset", "stage"]]
        .drop_duplicates()
        .sort_values(["cell_type_subset", "stage"])
        .itertuples(index=False, name=None)
    )
    groups = {(subtype, stage): np.flatnonzero((obs["cell_type_subset"] == subtype) & (obs["stage"] == stage)) for subtype, stage in group_order}

    for gene in LOX_GENES:
        if gene not in raw_genes:
            for (subtype, stage), idx in groups.items():
                summary_rows.append(
                    {
                        "gene": gene,
                        "cell_type_subset": subtype,
                        "stage": stage,
                        "n_cells": len(idx),
                        "dropout_rate": np.nan,
                        "mean_nonzero": np.nan,
                    }
                )
            continue

        expr = get_gene_vector(raw_x, raw_genes.get_loc(gene))
        per_group_counts: dict[str, dict[str, int]] = {}

        for (subtype, stage), idx in groups.items():
            group_expr = expr[idx]
            n_cells = len(group_expr)
            zero_count = int(np.sum(group_expr == 0))
            nonzero_count = int(n_cells - zero_count)
            mean_nonzero = float(group_expr[group_expr > 0].mean()) if nonzero_count else np.nan
            dropout_rate = float(zero_count / n_cells) if n_cells else np.nan

            summary_rows.append(
                {
                    "gene": gene,
                    "cell_type_subset": subtype,
                    "stage": stage,
                    "n_cells": n_cells,
                    "dropout_rate": dropout_rate,
                    "mean_nonzero": mean_nonzero,
                }
            )
            per_group_counts[f"{subtype}|{stage}"] = {
                "zero_count": zero_count,
                "nonzero_count": nonzero_count,
                "n_cells": n_cells,
            }

        for subtype in sorted(obs["cell_type_subset"].unique()):
            g02 = per_group_counts.get(f"{subtype}|02mo")
            g18 = per_group_counts.get(f"{subtype}|18mo")
            if g02 is None or g18 is None:
                odds_ratio = np.nan
                pvalue = np.nan
                zero_02mo = nonzero_02mo = zero_18mo = nonzero_18mo = np.nan
            else:
                zero_02mo = g02["zero_count"]
                nonzero_02mo = g02["nonzero_count"]
                zero_18mo = g18["zero_count"]
                nonzero_18mo = g18["nonzero_count"]
                odds_ratio, pvalue = fisher_exact(
                    [[zero_18mo, nonzero_18mo], [zero_02mo, nonzero_02mo]],
                    alternative="two-sided",
                )

            test_rows.append(
                {
                    "gene": gene,
                    "cell_type_subset": subtype,
                    "odds_ratio_18mo_vs_02mo": odds_ratio,
                    "pvalue": pvalue,
                    "zero_02mo": zero_02mo,
                    "nonzero_02mo": nonzero_02mo,
                    "zero_18mo": zero_18mo,
                    "nonzero_18mo": nonzero_18mo,
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    tests_df = pd.DataFrame(test_rows)
    tests_df["padj"] = benjamini_hochberg(tests_df["pvalue"].tolist())

    full_df = summary_df.merge(
        tests_df[["gene", "cell_type_subset", "odds_ratio_18mo_vs_02mo", "pvalue", "padj"]],
        on=["gene", "cell_type_subset"],
        how="left",
    )

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    full_df.to_csv(OUTPUT_CSV, index=False)

    print("\nDropout summary by gene x cell_type_subset x stage:")
    print(summary_df.to_string(index=False))

    print("\nFisher exact tests: dropout odds in 18mo vs 02mo")
    print(tests_df.to_string(index=False))

    print(f"\nSaved CSV: {OUTPUT_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
