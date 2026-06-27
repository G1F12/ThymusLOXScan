#!/usr/bin/env python
"""Regenerate descriptive fibroblast Spearman correlation table."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "data" / "processed" / "thymus_annotated.h5ad"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "sc_spearman_correlations.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate descriptive fibroblast Spearman correlation table."
    )
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def expression_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene not found in adata.var_names: {gene}")
    values = adata[:, gene].X
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def correlation_row(adata: ad.AnnData, gene1: str, gene2: str, stage: str) -> dict[str, object]:
    rho, pvalue = spearmanr(expression_vector(adata, gene1), expression_vector(adata, gene2))
    return {
        "gene1": gene1,
        "gene2": gene2,
        "stage": stage,
        "n_cells": int(adata.n_obs),
        "spearman_rho": float(rho),
        "pvalue": float(pvalue),
    }


def main() -> int:
    args = parse_args()
    adata = ad.read_h5ad(args.input_h5ad)
    required = {"stage", "cell_type"}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise KeyError(f"Input AnnData missing required obs columns: {missing}")

    fb = adata[adata.obs["cell_type"].astype(str).eq("FB").to_numpy(), :].copy()

    rows = [
        correlation_row(fb, "Loxl2", "Col1a1", "all"),
        correlation_row(fb, "Loxl2", "Vim", "all"),
        correlation_row(fb, "Loxl2", "Snai1", "all"),
    ]
    for stage in ["02mo", "18mo"]:
        stage_fb = fb[fb.obs["stage"].astype(str).eq(stage).to_numpy(), :].copy()
        rows.append(correlation_row(stage_fb, "Lox", "Loxl2", stage))

    table = pd.DataFrame(
        rows,
        columns=["gene1", "gene2", "stage", "n_cells", "spearman_rho", "pvalue"],
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    print(f"Saved {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
