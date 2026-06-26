#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PSEUDOBULK_CSV = PROJECT_ROOT / "results" / "pseudobulk_deseq2_LOX.csv"
DROPOUT_CSV = PROJECT_ROOT / "results" / "dropout_analysis_LOX.csv"
MAGIC_CSV = PROJECT_ROOT / "results" / "magic_imputation_LOX.csv"
OUTPUT_CSV = PROJECT_ROOT / "results" / "LOX_master_summary.csv"


def direction_consistent(log2fc: float, dropout_or: float) -> object:
    if pd.isna(log2fc) or pd.isna(dropout_or) or log2fc == 0 or dropout_or == 1:
        return pd.NA
    return (log2fc < 0 and dropout_or > 1) or (log2fc > 0 and dropout_or < 1)


def main() -> int:
    deseq = pd.read_csv(PSEUDOBULK_CSV)
    dropout = pd.read_csv(DROPOUT_CSV)
    magic = pd.read_csv(MAGIC_CSV)

    deseq_summary = deseq.rename(
        columns={
            "log2FoldChange": "log2FC_deseq2",
            "padj": "padj_deseq2",
        }
    )[["gene", "cell_type_subset", "log2FC_deseq2", "padj_deseq2"]]

    dropout_summary = (
        dropout.sort_values(["gene", "cell_type_subset", "stage"])
        .drop_duplicates(["gene", "cell_type_subset"])
        .rename(
            columns={
                "odds_ratio_18mo_vs_02mo": "dropout_OR",
                "padj": "dropout_padj",
            }
        )[["gene", "cell_type_subset", "dropout_OR", "dropout_padj"]]
    )

    magic_summary = (
        magic.sort_values(["gene", "cell_type_subset", "stage"])
        .drop_duplicates(["gene", "cell_type_subset"])
        .rename(columns={"direction_flip": "imputation_flip"})[
            ["gene", "cell_type_subset", "imputation_flip"]
        ]
    )

    summary = (
        deseq_summary.merge(dropout_summary, on=["gene", "cell_type_subset"], how="outer")
        .merge(magic_summary, on=["gene", "cell_type_subset"], how="outer")
        .sort_values(["cell_type_subset", "gene"])
        .reset_index(drop=True)
    )

    summary["direction_consistent"] = [
        direction_consistent(log2fc, dropout_or)
        for log2fc, dropout_or in zip(summary["log2FC_deseq2"], summary["dropout_OR"])
    ]

    significant = (summary["padj_deseq2"] < 0.05) | (summary["dropout_padj"] < 0.05)
    summary = summary.loc[
        significant.fillna(False),
        [
            "gene",
            "cell_type_subset",
            "log2FC_deseq2",
            "padj_deseq2",
            "dropout_OR",
            "dropout_padj",
            "imputation_flip",
            "direction_consistent",
        ],
    ]

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT_CSV, index=False)

    print(summary.to_string(index=False))
    print(f"\nSaved CSV: {OUTPUT_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
