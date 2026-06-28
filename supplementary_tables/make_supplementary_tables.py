#!/usr/bin/env python
"""Generate manuscript supplementary tables from existing repository outputs."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYDEPS = PROJECT_ROOT / "local_pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

import anndata as ad  # noqa: E402


INPUT_H5AD = PROJECT_ROOT / "data" / "processed" / "thymus_annotated.h5ad"
PSEUDOBULK_RESULTS = PROJECT_ROOT / "results" / "tables" / "lox_pseudobulk_complete_results.csv"
FB_SINGLE_CELL = PROJECT_ROOT / "results" / "sc_mannwhitney_FB_combined.csv"
MTEC1_SINGLE_CELL = PROJECT_ROOT / "results" / "sc_mannwhitney_mTEC1.csv"
CORRELATIONS = PROJECT_ROOT / "results" / "sc_spearman_correlations.csv"
OUTPUT_DIR = PROJECT_ROOT / "supplementary_tables"


def empty_string_for_missing(df: pd.DataFrame) -> pd.DataFrame:
    return df.where(pd.notna(df), "")


def load_obs() -> pd.DataFrame:
    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    obs = adata.obs[["sample", "age_group", "cell_type", "cell_type_subset"]].copy()
    adata.file.close()
    return obs.astype(str)


def write_cell_counts(obs: pd.DataFrame) -> None:
    table = (
        obs.groupby(["sample", "age_group", "cell_type", "cell_type_subset"], observed=True)
        .size()
        .reset_index(name="cell_count")
        .rename(columns={"sample": "sample_id", "cell_type_subset": "subtype"})
        .sort_values(["sample_id", "cell_type", "subtype"])
    )
    table = table[["sample_id", "age_group", "cell_type", "subtype", "cell_count"]]
    table.to_csv(OUTPUT_DIR / "Supplementary_Table_1_cell_counts.tsv", sep="\t", index=False)


def subtype_to_cell_type(obs: pd.DataFrame) -> dict[str, str]:
    mapping: dict[str, str] = {}
    counts = (
        obs.groupby(["cell_type_subset", "cell_type"], observed=True)
        .size()
        .reset_index(name="n")
        .sort_values(["cell_type_subset", "n"], ascending=[True, False])
    )
    for subtype, group in counts.groupby("cell_type_subset", sort=False):
        mapping[subtype] = group.iloc[0]["cell_type"]
    return mapping


def write_pseudobulk_results(obs: pd.DataFrame) -> None:
    pb = pd.read_csv(PSEUDOBULK_RESULTS)
    parent_cell_type = subtype_to_cell_type(obs)

    pb["cell_type"] = ""
    pb["subtype"] = ""
    is_cell_type = pb["annotation_level"].eq("cell_type")
    pb.loc[is_cell_type, "cell_type"] = pb.loc[is_cell_type, "cell_type_or_subtype"]
    pb.loc[~is_cell_type, "subtype"] = pb.loc[~is_cell_type, "cell_type_or_subtype"]
    pb.loc[~is_cell_type, "cell_type"] = pb.loc[~is_cell_type, "cell_type_or_subtype"].map(parent_cell_type)

    table = pb.rename(columns={"age_comparison": "comparison"})[
        [
            "gene",
            "cell_type",
            "subtype",
            "comparison",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "pvalue",
            "padj",
            "n_young_samples",
            "n_aged_samples",
        ]
    ].sort_values(["cell_type", "subtype", "gene"])
    empty_string_for_missing(table).to_csv(
        OUTPUT_DIR / "Supplementary_Table_2_pseudobulk_LOX_results.tsv",
        sep="\t",
        index=False,
    )


def single_cell_table(path: Path, compartment: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    return pd.DataFrame(
        {
            "gene": df["gene"],
            "compartment": compartment,
            "n_young_cells": df["n_02mo"],
            "n_aged_cells": df["n_18mo"],
            "test": "Mann-Whitney U test",
            "pvalue": df["pvalue"],
            "adjusted_pvalue": "",
            "rank_biserial_effect_size": df["rank_biserial"],
        }
    )


def write_single_cell_tests() -> None:
    table = pd.concat(
        [
            single_cell_table(FB_SINGLE_CELL, "pooled fibroblasts"),
            single_cell_table(MTEC1_SINGLE_CELL, "mTEC1"),
        ],
        ignore_index=True,
    ).sort_values(["compartment", "gene"])
    table.to_csv(OUTPUT_DIR / "Supplementary_Table_3_single_cell_tests.tsv", sep="\t", index=False)


def write_correlations() -> None:
    corr = pd.read_csv(CORRELATIONS)
    table = pd.DataFrame(
        {
            "gene1": corr["gene1"],
            "gene2": corr["gene2"],
            "compartment": "pooled fibroblasts; stage=" + corr["stage"].astype(str),
            "method": "Spearman correlation",
            "rho": corr["spearman_rho"],
            "pvalue": corr["pvalue"],
            "adjusted_pvalue": "",
        }
    ).sort_values(["compartment", "gene1", "gene2"])
    table.to_csv(OUTPUT_DIR / "Supplementary_Table_4_correlations.tsv", sep="\t", index=False)


def write_readme() -> None:
    readme = """# Supplementary tables

These tables are generated from existing repository data and analysis outputs. They do not rerun statistical models or alter scientific results.

Regenerate all supplementary tables from the repository root with:

```bash
python supplementary_tables/make_supplementary_tables.py
```

## Supplementary_Table_1_cell_counts.tsv

Generated from `data/processed/thymus_annotated.h5ad` using the annotated cell metadata columns `sample`, `age_group`, `cell_type`, and `cell_type_subset`. Rows are sample-by-age-by-cell-type-by-subtype counts. The `subtype` column corresponds to `cell_type_subset` in the AnnData object.

## Supplementary_Table_2_pseudobulk_LOX_results.tsv

Generated from `results/tables/lox_pseudobulk_complete_results.csv`, which contains the LOX-family DESeq2 pseudobulk results. The table preserves the reported DESeq2 statistics and replicate counts. Rows analyzed at the broad cell-type level have an empty `subtype`; rows analyzed at the subtype level have `cell_type` inferred from the majority parent annotation in `data/processed/thymus_annotated.h5ad`.

## Supplementary_Table_3_single_cell_tests.tsv

Generated from `results/sc_mannwhitney_FB_combined.csv` and `results/sc_mannwhitney_mTEC1.csv`. These are descriptive single-cell-level Mann-Whitney U tests for pooled fibroblasts and mTEC1, respectively. The source files do not include multiple-testing adjusted p-values, so `adjusted_pvalue` is left blank.

## Supplementary_Table_4_correlations.tsv

Generated from `results/sc_spearman_correlations.csv`. These are descriptive Spearman correlations reported for pooled fibroblast analyses, with the original `stage` column encoded in `compartment` as `pooled fibroblasts; stage=<stage>`. The source file does not include multiple-testing adjusted p-values, so `adjusted_pvalue` is left blank.
"""
    (OUTPUT_DIR / "README.md").write_text(readme, encoding="utf-8")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    obs = load_obs()
    write_cell_counts(obs)
    write_pseudobulk_results(obs)
    write_single_cell_tests()
    write_correlations()
    write_readme()
    print(f"Supplementary tables written to {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
