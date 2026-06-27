#!/usr/bin/env python
"""Regenerate descriptive single-cell Mann-Whitney LOX-family test tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "data" / "processed" / "thymus_annotated.h5ad"
DEFAULT_FB_OUTPUT = PROJECT_ROOT / "results" / "sc_mannwhitney_FB_combined.csv"
DEFAULT_MTEC1_OUTPUT = PROJECT_ROOT / "results" / "sc_mannwhitney_mTEC1.csv"

LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
YOUNG = "02mo"
AGED = "18mo"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate descriptive single-cell Mann-Whitney LOX-family tables."
    )
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--fb-output", type=Path, default=DEFAULT_FB_OUTPUT)
    parser.add_argument("--mtec1-output", type=Path, default=DEFAULT_MTEC1_OUTPUT)
    return parser.parse_args()


def expression_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene not found in adata.var_names: {gene}")
    values = adata[:, gene].X
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def make_table(adata: ad.AnnData, mask: pd.Series) -> pd.DataFrame:
    obs = adata.obs.loc[mask, ["stage"]].copy()
    young_mask = obs["stage"].astype(str).eq(YOUNG).to_numpy()
    aged_mask = obs["stage"].astype(str).eq(AGED).to_numpy()

    rows: list[dict[str, object]] = []
    subset = adata[mask.to_numpy(), :].copy()
    for gene in LOX_GENES:
        values = expression_vector(subset, gene)
        young = values[young_mask]
        aged = values[aged_mask]
        # Existing repository outputs store the U statistic from young-vs-aged
        # ordering, while the rank-biserial effect is oriented aged-vs-young.
        test = mannwhitneyu(young, aged, alternative="two-sided")
        statistic = float(test.statistic)
        rank_biserial = float(1 - (2 * statistic / (len(aged) * len(young))))
        rows.append(
            {
                "gene": gene,
                "n_02mo": int(len(young)),
                "n_18mo": int(len(aged)),
                "U_statistic": statistic,
                "pvalue": float(test.pvalue),
                "rank_biserial": rank_biserial,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    adata = ad.read_h5ad(args.input_h5ad)
    required = {"stage", "cell_type", "cell_type_subset"}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise KeyError(f"Input AnnData missing required obs columns: {missing}")

    obs = adata.obs
    fb_mask = obs["cell_type"].astype(str).eq("FB")
    mtec1_mask = obs["cell_type_subset"].astype(str).eq("13:mTEC1")

    fb_table = make_table(adata, fb_mask)
    mtec1_table = make_table(adata, mtec1_mask)

    args.fb_output.parent.mkdir(parents=True, exist_ok=True)
    args.mtec1_output.parent.mkdir(parents=True, exist_ok=True)
    fb_table.to_csv(args.fb_output, index=False)
    mtec1_table.to_csv(args.mtec1_output, index=False)
    print(f"Saved {args.fb_output}")
    print(f"Saved {args.mtec1_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
