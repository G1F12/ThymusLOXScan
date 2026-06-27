#!/usr/bin/env python
"""Sensitivity analysis for broad fibroblast LOX age effects adjusted for subtype mix."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
DEFAULT_RESULTS_TABLE = PROJECT_ROOT / "results" / "tables" / "lox_pseudobulk_complete_results.csv"
DEFAULT_OUTPUT_TABLE = PROJECT_ROOT / "results" / "tables" / "fb_composition_adjusted_lox.tsv"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "fb_composition_adjustment.md"

GENES = ["Lox", "Loxl1", "Loxl2"]
YOUNG_LABEL = "02mo"
AGED_LABEL = "18mo"
FB_CELL_TYPE = "FB"
SUBTYPE_LABELS = {
    "capsFB": "3:capsFB",
    "intFB": "4:intFB",
    "medFB": "5:medFB",
    "Fat": "9:Fat",
}
CPM_SCALE = 1_000_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--results-table", type=Path, default=DEFAULT_RESULTS_TABLE)
    parser.add_argument("--output-table", type=Path, default=DEFAULT_OUTPUT_TABLE)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def vector_to_array(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def sample_summaries(adata) -> pd.DataFrame:
    if adata.raw is None:
        raise ValueError("Input AnnData must contain adata.raw.X raw counts.")

    raw_genes = pd.Index(adata.raw.var_names.astype(str))
    missing = [gene for gene in GENES if gene not in raw_genes]
    if missing:
        raise KeyError(f"Genes absent from adata.raw.var_names: {missing}")

    obs = adata.obs[["sample", "stage", "cell_type", "cell_type_subset"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["cell_type_subset"] = obs["cell_type_subset"].astype(str)
    obs["_row_pos"] = np.arange(adata.n_obs)
    fb_obs = obs.loc[obs["cell_type"].eq(FB_CELL_TYPE)].copy()

    raw_x = adata.raw.X
    if not sparse.issparse(raw_x):
        raw_x = np.asarray(raw_x)

    rows = []
    for sample, sample_obs in fb_obs.groupby("sample", sort=True, observed=True):
        stages = sample_obs["stage"].unique()
        if len(stages) != 1:
            raise ValueError(f"Sample {sample} has multiple stage labels: {stages}")

        row_idx = sample_obs["_row_pos"].to_numpy()
        sample_matrix = raw_x[row_idx, :]
        total_counts = float(np.asarray(sample_matrix.sum()).ravel()[0])
        n_fb_cells = int(len(sample_obs))
        subtype_counts = sample_obs["cell_type_subset"].value_counts()

        base = {
            "sample": sample,
            "stage": stages[0],
            "age_18mo": 1 if stages[0] == AGED_LABEL else 0,
            "n_fb_cells": n_fb_cells,
            "fb_total_counts": total_counts,
        }
        for short_name, full_label in SUBTYPE_LABELS.items():
            count = int(subtype_counts.get(full_label, 0))
            base[f"{short_name}_cells"] = count
            base[f"{short_name}_fraction"] = count / n_fb_cells if n_fb_cells else np.nan

        for gene in GENES:
            gene_counts = float(vector_to_array(raw_x[row_idx, raw_genes.get_loc(gene)]).sum())
            cpm = (gene_counts / total_counts * CPM_SCALE) if total_counts > 0 else np.nan
            row = base | {
                "gene": gene,
                "gene_counts": gene_counts,
                "cpm": cpm,
                "log2_cpm_plus1": np.log2(cpm + 1) if np.isfinite(cpm) else np.nan,
            }
            rows.append(row)

    return pd.DataFrame(rows)


def fit_age_model(df: pd.DataFrame, covariates: list[str]) -> dict[str, float | str]:
    fit_df = df[["log2_cpm_plus1", "age_18mo", *covariates]].dropna().copy()
    y = fit_df["log2_cpm_plus1"].astype(float)
    x = sm.add_constant(fit_df[["age_18mo", *covariates]].astype(float), has_constant="add")
    model = sm.OLS(y, x).fit()
    coef = float(model.params.get("age_18mo", np.nan))
    pvalue = float(model.pvalues.get("age_18mo", np.nan)) if model.df_resid > 0 else np.nan
    se = float(model.bse.get("age_18mo", np.nan)) if model.df_resid > 0 else np.nan
    return {
        "age_coef_log2_cpm_plus1": coef,
        "age_se": se,
        "age_pvalue": pvalue,
        "n_samples": int(model.nobs),
        "model_rank": int(model.df_model + 1),
        "df_resid": float(model.df_resid),
        "r_squared": float(model.rsquared),
        "covariates": ",".join(covariates) if covariates else "none",
        "model_note": (
            f"extremely small-n OLS sensitivity model; df_resid={model.df_resid:g}"
            if covariates
            else f"sample-level unadjusted OLS comparison; df_resid={model.df_resid:g}"
        ),
    }


def broad_deseq2_effects(results_table: Path) -> pd.DataFrame:
    results = pd.read_csv(results_table)
    mask = (
        results["annotation_level"].eq("cell_type")
        & results["cell_type_or_subtype"].eq(FB_CELL_TYPE)
        & results["gene"].isin(GENES)
    )
    cols = [
        "gene",
        "log2FoldChange",
        "lfcSE",
        "pvalue",
        "padj",
        "significance_label",
        "direction",
    ]
    return results.loc[mask, cols].rename(
        columns={
            "log2FoldChange": "broad_fb_deseq2_log2fc",
            "lfcSE": "broad_fb_deseq2_lfcSE",
            "pvalue": "broad_fb_deseq2_pvalue",
            "padj": "broad_fb_deseq2_padj",
            "significance_label": "broad_fb_deseq2_significance",
            "direction": "broad_fb_deseq2_direction",
        }
    )


def subtype_effects(results_table: Path) -> pd.DataFrame:
    results = pd.read_csv(results_table)
    rows = []
    wanted = set(SUBTYPE_LABELS.values())
    for gene in GENES:
        parts = []
        for subtype in SUBTYPE_LABELS.values():
            hit = results.loc[
                results["annotation_level"].eq("subtype")
                & results["cell_type_or_subtype"].eq(subtype)
                & results["gene"].eq(gene)
            ]
            if hit.empty:
                parts.append(f"{subtype}: not_available")
                continue
            row = hit.iloc[0]
            parts.append(
                f"{subtype}: log2FC={row['log2FoldChange']:.3g}, "
                f"padj={row['padj']:.3g}, {row['direction']}"
            )
        rows.append({"gene": gene, "subtype_stratified_effects": "; ".join(parts)})
    return pd.DataFrame(rows)


def classify(row: pd.Series) -> str:
    broad = row["broad_fb_deseq2_log2fc"]
    adjusted = row["adjusted_age_coef_log2_cpm_plus1"]
    if pd.isna(adjusted):
        return "downgrade: adjusted model unavailable"
    if np.sign(broad) != np.sign(adjusted):
        return "downgrade: adjusted sample-level age coefficient changes sign"
    ratio = abs(adjusted) / abs(row["unadjusted_age_coef_log2_cpm_plus1"])
    if ratio < 0.5:
        return "downgrade: adjusted age coefficient is strongly attenuated"
    return "retain: adjusted coefficient keeps the broad FB direction"


def build_results(summaries: pd.DataFrame, results_table: Path) -> pd.DataFrame:
    fraction_covariates = [f"{name}_fraction" for name in SUBTYPE_LABELS]
    rows = []
    for gene, gene_df in summaries.groupby("gene", sort=False):
        unadjusted = fit_age_model(gene_df, [])
        adjusted = fit_age_model(gene_df, fraction_covariates)
        rows.append(
            {
                "gene": gene,
                "unadjusted_age_coef_log2_cpm_plus1": unadjusted["age_coef_log2_cpm_plus1"],
                "unadjusted_age_se": unadjusted["age_se"],
                "unadjusted_age_pvalue": unadjusted["age_pvalue"],
                "unadjusted_r_squared": unadjusted["r_squared"],
                "adjusted_age_coef_log2_cpm_plus1": adjusted["age_coef_log2_cpm_plus1"],
                "adjusted_age_se": adjusted["age_se"],
                "adjusted_age_pvalue": adjusted["age_pvalue"],
                "adjusted_r_squared": adjusted["r_squared"],
                "adjusted_n_samples": adjusted["n_samples"],
                "adjusted_model_rank": adjusted["model_rank"],
                "adjusted_df_resid": adjusted["df_resid"],
                "adjusted_covariates": adjusted["covariates"],
                "adjusted_model_note": adjusted["model_note"],
            }
        )

    out = pd.DataFrame(rows)
    out = out.merge(broad_deseq2_effects(results_table), on="gene", how="left")
    out = out.merge(subtype_effects(results_table), on="gene", how="left")
    out["composition_sensitivity_conclusion"] = out.apply(classify, axis=1)

    fraction_summary = (
        summaries.drop_duplicates("sample")
        .sort_values("sample")
        [["sample", "stage", "n_fb_cells", *fraction_covariates]]
    )
    for _, sample_row in fraction_summary.iterrows():
        suffix = sample_row["sample"]
        for col in ["stage", "n_fb_cells", *fraction_covariates]:
            out[f"{suffix}_{col}"] = sample_row[col]
    return out


def write_report(table: pd.DataFrame, summaries: pd.DataFrame, report_path: Path) -> None:
    sample_rows = summaries.drop_duplicates("sample").sort_values("sample")
    fraction_cols = [f"{name}_fraction" for name in SUBTYPE_LABELS]

    retained = table["composition_sensitivity_conclusion"].str.startswith("retain").sum()
    downgraded = len(table) - retained
    overall = (
        "Broad FB conclusions should be downgraded where the adjusted age coefficient changes sign "
        "or is materially attenuated; retained where direction is stable."
    )
    if downgraded == 0:
        overall = (
            "The directional broad FB conclusions can be retained, but their inferential strength should be "
            "downgraded to a small-n composition-sensitivity result rather than treated as definitive "
            "broad-FB-intrinsic effects."
        )
    elif retained == 0:
        overall = "Broad FB conclusions should be downgraded because adjustment materially changes all focal genes."

    lines = [
        "# Fibroblast subtype-composition adjustment for LOX-family age effects",
        "",
        "## Scope",
        "",
        "This is a small-sample sensitivity analysis using sample-level broad fibroblast pseudobulk expression. The adjusted model regresses log2(CPM + 1) on age and the capsFB, intFB, medFB, and Fat fractions within broad FB. With only six broad-FB samples, adjusted coefficients are descriptive rather than definitive inference.",
        "",
        "## Overall conclusion",
        "",
        overall,
        "",
        "Because the adjusted model has only one residual degree of freedom after including age plus four subtype fractions, p-values are not interpreted for the adjusted fits. The useful signal is whether the age coefficient keeps its sign compared with the unadjusted sample-level broad FB effect and the subtype-stratified DESeq2 effects.",
        "",
        "## Sample subtype composition",
        "",
        "| sample | stage | FB cells | capsFB | intFB | medFB | Fat |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in sample_rows.iterrows():
        lines.append(
            f"| {row['sample']} | {row['stage']} | {int(row['n_fb_cells'])} | "
            f"{row['capsFB_fraction']:.3f} | {row['intFB_fraction']:.3f} | "
            f"{row['medFB_fraction']:.3f} | {row['Fat_fraction']:.3f} |"
        )

    lines.extend(["", "## Gene-level comparison", ""])
    for _, row in table.iterrows():
        lines.extend(
            [
                f"### {row['gene']}",
                "",
                f"- Broad FB DESeq2 effect: log2FC={row['broad_fb_deseq2_log2fc']:.3g}, padj={row['broad_fb_deseq2_padj']:.3g}, {row['broad_fb_deseq2_direction']}.",
                f"- Sample-level unadjusted age coefficient: {row['unadjusted_age_coef_log2_cpm_plus1']:.3g} log2(CPM+1).",
                f"- Composition-adjusted age coefficient: {row['adjusted_age_coef_log2_cpm_plus1']:.3g} log2(CPM+1); {row['adjusted_model_note']}.",
                f"- Subtype-stratified effects: {row['subtype_stratified_effects']}.",
                f"- Conclusion: {row['composition_sensitivity_conclusion']}.",
                "",
            ]
        )

    lines.extend(
        [
            "## Interpretation",
            "",
            "This analysis should be treated as a sensitivity check. The broad FB DESeq2 model remains the primary broad-cell-type result, but the adjusted coefficients test whether those signals plausibly survive subtype-mixture imbalance. Any downgraded gene should be described as composition-sensitive rather than a clean broad-FB-intrinsic age effect.",
        ]
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    print(f"Loading {args.input_h5ad}", flush=True)
    adata = sc.read_h5ad(args.input_h5ad)
    summaries = sample_summaries(adata)
    table = build_results(summaries, args.results_table)

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_table, sep="\t", index=False)
    write_report(table, summaries, args.report)

    print(f"Saved table: {args.output_table}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
