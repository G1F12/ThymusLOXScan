#!/usr/bin/env python
"""Global Benjamini-Hochberg correction across GSE240016 LOX-family tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_TABLE = PROJECT_ROOT / "results" / "tables" / "lox_pseudobulk_complete_results.tsv"
FALLBACK_INPUT = PROJECT_ROOT / "results" / "pseudobulk_deseq2_LOX.csv"
OUTPUT_TABLE = PROJECT_ROOT / "results" / "tables" / "lox_family_global_fdr_matrix.tsv"
REPORT_PATH = PROJECT_ROOT / "reports" / "lox_family_global_fdr_matrix.md"

LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pvalues.astype(float)
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return out
    ordered = valid.sort_values()
    n = len(ordered)
    adjusted = ordered * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted.iloc[::-1]).iloc[::-1].clip(upper=1.0)
    out.loc[adjusted.index] = adjusted
    return out


def load_results() -> pd.DataFrame:
    if INPUT_TABLE.exists():
        df = pd.read_csv(INPUT_TABLE, sep="\t")
        group_col = "cell_type_or_subtype"
    elif FALLBACK_INPUT.exists():
        df = pd.read_csv(FALLBACK_INPUT)
        df["annotation_level"] = "subtype"
        group_col = "cell_type_subset"
    else:
        raise FileNotFoundError(f"Missing {INPUT_TABLE} and {FALLBACK_INPUT}")
    df = df.loc[df["gene"].isin(LOX_GENES)].copy()
    df = df.rename(
        columns={
            group_col: "group_or_subtype",
            "pvalue": "raw_pvalue",
            "padj": "original_padj",
            "log2FoldChange": "log2FC_18mo_vs_02mo",
        }
    )
    for col in ["raw_pvalue", "original_padj", "log2FC_18mo_vs_02mo"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "annotation_level" not in df.columns:
        df["annotation_level"] = "subtype"
    return df


def make_matrix(df: pd.DataFrame) -> pd.DataFrame:
    out_cols = [
        "gene",
        "annotation_level",
        "group_or_subtype",
        "status",
        "n_young_samples",
        "n_aged_samples",
        "young_cell_count",
        "aged_cell_count",
        "raw_pvalue",
        "original_padj",
        "log2FC_18mo_vs_02mo",
    ]
    for col in out_cols:
        if col not in df.columns:
            df[col] = np.nan
    matrix = df[out_cols].copy()
    matrix["global_bh_fdr_lox_x_subtype_matrix"] = bh_fdr(matrix["raw_pvalue"])
    matrix["global_fdr_classification"] = np.select(
        [
            matrix["global_bh_fdr_lox_x_subtype_matrix"].lt(0.05),
            matrix["global_bh_fdr_lox_x_subtype_matrix"].lt(0.10),
            matrix["raw_pvalue"].notna(),
        ],
        ["global_FDR_lt_0.05", "global_FDR_lt_0.10", "not_significant_after_global_matrix_correction"],
        default="not_tested",
    )
    return matrix.sort_values(["annotation_level", "group_or_subtype", "gene"]).reset_index(drop=True)


def write_report(matrix: pd.DataFrame) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    tested = matrix.loc[matrix["raw_pvalue"].notna()].copy()
    sig05 = tested.loc[tested["global_bh_fdr_lox_x_subtype_matrix"].lt(0.05)]
    sig10 = tested.loc[
        tested["global_bh_fdr_lox_x_subtype_matrix"].ge(0.05)
        & tested["global_bh_fdr_lox_x_subtype_matrix"].lt(0.10)
    ]
    lines = [
        "# LOX-family global FDR matrix",
        "",
        "## Scope",
        "",
        "This report applies a single Benjamini-Hochberg correction across the internal GSE240016 LOX-family pseudobulk matrix only. External datasets are not mixed into this FDR calculation.",
        "",
        f"Tested rows with raw p-values: {len(tested)}.",
        "",
        "## Global FDR < 0.05",
        "",
    ]
    if sig05.empty:
        lines.append("- None.")
    else:
        for _, row in sig05.sort_values("global_bh_fdr_lox_x_subtype_matrix").iterrows():
            lines.append(
                f"- {row['group_or_subtype']} / {row['gene']}: log2FC={row['log2FC_18mo_vs_02mo']:.3f}, "
                f"raw p={row['raw_pvalue']:.3g}, global FDR={row['global_bh_fdr_lox_x_subtype_matrix']:.3g}."
            )
    lines.extend(["", "## Global FDR < 0.10 only", ""])
    if sig10.empty:
        lines.append("- None.")
    else:
        for _, row in sig10.sort_values("global_bh_fdr_lox_x_subtype_matrix").iterrows():
            lines.append(
                f"- {row['group_or_subtype']} / {row['gene']}: log2FC={row['log2FC_18mo_vs_02mo']:.3f}, "
                f"raw p={row['raw_pvalue']:.3g}, global FDR={row['global_bh_fdr_lox_x_subtype_matrix']:.3g}."
            )
    nonsig = tested.loc[tested["global_bh_fdr_lox_x_subtype_matrix"].ge(0.10)]
    lines.extend(
        [
            "",
            "## Not significant after global matrix correction",
            "",
            f"{len(nonsig)} tested rows were not significant at global FDR < 0.10.",
            "",
            "## Output",
            "",
            f"- `{OUTPUT_TABLE.relative_to(PROJECT_ROOT).as_posix()}`",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUTPUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    matrix = make_matrix(load_results())
    matrix.to_csv(OUTPUT_TABLE, sep="\t", index=False)
    write_report(matrix)
    print(f"Saved {OUTPUT_TABLE}")
    print(f"Saved {REPORT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
