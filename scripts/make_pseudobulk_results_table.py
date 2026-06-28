#!/usr/bin/env python
"""Create complete LOX-family pseudobulk DESeq2 result tables.

This script intentionally uses biological samples as replicates. For each
annotation group, raw counts are summed within sample before DESeq2 is run.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "tables"

LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
YOUNG_LABEL = "02mo"
AGED_LABEL = "18mo"
AGE_COMPARISON = f"{AGED_LABEL}_vs_{YOUNG_LABEL}"
MIN_CELLS_PER_SAMPLE_GROUP = 10
MIN_REPLICATES_PER_GROUP = 3
LOW_EXPRESSION_BASEMEAN = 10.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate complete LOX-family pseudobulk DESeq2 tables."
    )
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=DEFAULT_INPUT_H5AD,
        help=f"Input AnnData file. Default: {DEFAULT_INPUT_H5AD}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=MIN_CELLS_PER_SAMPLE_GROUP,
        help="Minimum cells required for one sample-by-annotation pseudobulk replicate.",
    )
    parser.add_argument(
        "--low-expression-basemean",
        type=float,
        default=LOW_EXPRESSION_BASEMEAN,
        help="Warn when a tested LOX gene has baseMean below this threshold.",
    )
    return parser.parse_args()


def sparse_row_sum(matrix, row_indices: np.ndarray) -> np.ndarray:
    selected = matrix[row_indices, :]
    summed = np.asarray(selected.sum(axis=0)).ravel()
    return np.rint(summed).astype(np.int64)


def clean_label(value: str) -> str:
    return value.replace(":", "_").replace("/", "_").replace(" ", "_")


def empty_result_row(
    *,
    annotation_level: str,
    annotation: str,
    gene: str,
    present: bool,
    status: str,
    n_young_samples: int,
    n_aged_samples: int,
    young_cell_count: int,
    aged_cell_count: int,
) -> dict[str, object]:
    return {
        "annotation_level": annotation_level,
        "cell_type_or_subtype": annotation,
        "gene": gene,
        "age_comparison": AGE_COMPARISON,
        "n_young_samples": n_young_samples,
        "n_aged_samples": n_aged_samples,
        "young_cell_count": young_cell_count,
        "aged_cell_count": aged_cell_count,
        "baseMean": np.nan,
        "log2FoldChange": np.nan,
        "lfcSE": np.nan,
        "stat": np.nan,
        "pvalue": np.nan,
        "padj": np.nan,
        "significance_label": status,
        "direction": "not_tested",
        "status": status,
        "gene_present": present,
    }


def build_pseudobulk(
    adata: AnnData,
    group_col: str,
    group_value: str,
    min_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if adata.raw is None:
        raise ValueError("adata.raw is missing; pseudobulk requires raw counts in adata.raw.X.")

    obs = adata.obs[["sample", "stage", group_col]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["stage"] = obs["stage"].astype(str)
    obs[group_col] = obs[group_col].astype(str)
    obs["_row_pos"] = np.arange(adata.n_obs)
    obs = obs.loc[obs[group_col] == group_value].copy()

    raw_x = adata.raw.X
    if not sparse.issparse(raw_x):
        raw_x = np.asarray(raw_x)
    genes = pd.Index(adata.raw.var_names.astype(str))

    count_rows: list[np.ndarray] = []
    meta_rows: list[dict[str, object]] = []
    for sample, sample_obs in obs.groupby("sample", sort=True, observed=True):
        n_cells = len(sample_obs)
        if n_cells < min_cells:
            continue

        stages = sample_obs["stage"].unique()
        if len(stages) != 1:
            raise ValueError(
                f"Expected one stage per sample for {group_col}={group_value}; "
                f"got {stages} for sample={sample}."
            )

        sample_id = f"{sample}__{clean_label(group_value)}"
        count_rows.append(sparse_row_sum(raw_x, sample_obs["_row_pos"].to_numpy()))
        meta_rows.append(
            {
                "sample_id": sample_id,
                "sample": sample,
                "stage": stages[0],
                "annotation": group_value,
                "n_cells": n_cells,
            }
        )

    counts = pd.DataFrame(count_rows, index=[row["sample_id"] for row in meta_rows], columns=genes)
    metadata = pd.DataFrame(meta_rows)
    if metadata.empty:
        metadata = pd.DataFrame(columns=["sample", "stage", "annotation", "n_cells"])
        metadata.index.name = "sample_id"
        return counts, metadata

    metadata = metadata.set_index("sample_id")
    metadata["stage"] = pd.Categorical(
        metadata["stage"], categories=[YOUNG_LABEL, AGED_LABEL], ordered=False
    )
    return counts, metadata


def significance_label(padj: object, pvalue: object, status: str) -> str:
    if status != "ok":
        return status
    if pd.isna(padj):
        return "padj_NA"
    if padj < 0.05:
        return "FDR<0.05"
    if pd.notna(pvalue) and pvalue < 0.05:
        return "nominal_p<0.05"
    return "not_significant"


def direction(log2fc: object, status: str) -> str:
    if status != "ok" or pd.isna(log2fc):
        return "not_tested"
    if np.isclose(float(log2fc), 0.0):
        return "no_change"
    return "up_in_aged" if float(log2fc) > 0 else "down_in_aged"


def warn_replicates(annotation_level: str, annotation: str, n_young: int, n_aged: int) -> None:
    if n_young < MIN_REPLICATES_PER_GROUP or n_aged < MIN_REPLICATES_PER_GROUP:
        warnings.warn(
            f"{annotation_level}={annotation}: fewer than {MIN_REPLICATES_PER_GROUP} "
            f"biological replicates in at least one group "
            f"({YOUNG_LABEL}={n_young}, {AGED_LABEL}={n_aged}).",
            RuntimeWarning,
            stacklevel=2,
        )


def run_deseq_for_annotation(
    *,
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    annotation_level: str,
    annotation: str,
    raw_genes: pd.Index,
    low_expression_basemean: float,
) -> pd.DataFrame:
    stage_counts = metadata["stage"].value_counts() if not metadata.empty else pd.Series(dtype=int)
    n_young = int(stage_counts.get(YOUNG_LABEL, 0))
    n_aged = int(stage_counts.get(AGED_LABEL, 0))
    young_cells = (
        int(metadata.loc[metadata["stage"] == YOUNG_LABEL, "n_cells"].sum())
        if not metadata.empty
        else 0
    )
    aged_cells = (
        int(metadata.loc[metadata["stage"] == AGED_LABEL, "n_cells"].sum())
        if not metadata.empty
        else 0
    )

    warn_replicates(annotation_level, annotation, n_young, n_aged)

    if metadata.empty:
        return pd.DataFrame(
            [
                empty_result_row(
                    annotation_level=annotation_level,
                    annotation=annotation,
                    gene=gene,
                    present=gene in raw_genes,
                    status="skipped_no_pseudobulk_samples",
                    n_young_samples=n_young,
                    n_aged_samples=n_aged,
                    young_cell_count=young_cells,
                    aged_cell_count=aged_cells,
                )
                for gene in LOX_GENES
            ]
        )

    observed_stages = set(stage_counts.index[stage_counts > 0].astype(str))
    if observed_stages != {YOUNG_LABEL, AGED_LABEL}:
        return pd.DataFrame(
            [
                empty_result_row(
                    annotation_level=annotation_level,
                    annotation=annotation,
                    gene=gene,
                    present=gene in raw_genes,
                    status="skipped_missing_stage",
                    n_young_samples=n_young,
                    n_aged_samples=n_aged,
                    young_cell_count=young_cells,
                    aged_cell_count=aged_cells,
                )
                for gene in LOX_GENES
            ]
        )

    nonzero_genes = counts.sum(axis=0) > 0
    deseq_counts = counts.loc[:, nonzero_genes].copy()

    dds = DeseqDataSet(
        counts=deseq_counts,
        metadata=metadata[["stage"]],
        design="~ stage",
        ref_level=["stage", YOUNG_LABEL],
        refit_cooks=False,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["stage", AGED_LABEL, YOUNG_LABEL], quiet=True)
    stats.summary()
    deseq_results = stats.results_df.copy()

    rows: list[dict[str, object]] = []
    for gene in LOX_GENES:
        present = gene in raw_genes
        if gene in deseq_results.index:
            result = deseq_results.loc[gene]
            status = "ok"
            base_mean = result.get("baseMean", np.nan)
            padj = result.get("padj", np.nan)
            pvalue = result.get("pvalue", np.nan)
            if pd.isna(padj):
                warnings.warn(
                    f"{annotation_level}={annotation}, gene={gene}: padj is NA.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            if pd.notna(base_mean) and base_mean < low_expression_basemean:
                warnings.warn(
                    f"{annotation_level}={annotation}, gene={gene}: low expression "
                    f"baseMean={base_mean:.3g} < {low_expression_basemean}.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            row = {
                "annotation_level": annotation_level,
                "cell_type_or_subtype": annotation,
                "gene": gene,
                "age_comparison": AGE_COMPARISON,
                "n_young_samples": n_young,
                "n_aged_samples": n_aged,
                "young_cell_count": young_cells,
                "aged_cell_count": aged_cells,
                "baseMean": base_mean,
                "log2FoldChange": result.get("log2FoldChange", np.nan),
                "lfcSE": result.get("lfcSE", np.nan),
                "stat": result.get("stat", np.nan),
                "pvalue": pvalue,
                "padj": padj,
                "status": status,
                "gene_present": present,
            }
            row["significance_label"] = significance_label(row["padj"], row["pvalue"], status)
            row["direction"] = direction(row["log2FoldChange"], status)
        else:
            status = "not_tested_zero_counts" if present else "gene_absent"
            row = empty_result_row(
                annotation_level=annotation_level,
                annotation=annotation,
                gene=gene,
                present=present,
                status=status,
                n_young_samples=n_young,
                n_aged_samples=n_aged,
                young_cell_count=young_cells,
                aged_cell_count=aged_cells,
            )
        rows.append(row)

    return pd.DataFrame(rows)


def annotation_values(adata: AnnData, group_col: str) -> list[str]:
    values = sorted(adata.obs[group_col].astype(str).dropna().unique())
    return [value for value in values if value and value.lower() != "nan"]


def main() -> int:
    args = parse_args()
    print(f"Loading input AnnData: {args.input_h5ad}", flush=True)
    adata = sc.read_h5ad(args.input_h5ad)
    if adata.raw is None:
        raise ValueError("Input AnnData must contain adata.raw with raw counts.")

    required_obs = {"sample", "stage", "cell_type", "cell_type_subset"}
    missing = sorted(required_obs - set(adata.obs.columns))
    if missing:
        raise KeyError(f"Input AnnData is missing required obs columns: {missing}")

    raw_genes = pd.Index(adata.raw.var_names.astype(str))
    print(f"Loaded {adata.n_obs:,} cells x {adata.raw.n_vars:,} raw genes.", flush=True)
    print("Running pseudobulk DESeq2 with biological samples as replicates.", flush=True)

    tables: list[pd.DataFrame] = []
    for annotation_level, group_col in [
        ("cell_type", "cell_type"),
        ("subtype", "cell_type_subset"),
    ]:
        values = annotation_values(adata, group_col)
        print(f"\nProcessing {annotation_level}: {len(values)} groups", flush=True)
        for i, value in enumerate(values, start=1):
            print(f"  [{i}/{len(values)}] {value}", flush=True)
            try:
                counts, metadata = build_pseudobulk(adata, group_col, value, args.min_cells)
                result = run_deseq_for_annotation(
                    counts=counts,
                    metadata=metadata,
                    annotation_level=annotation_level,
                    annotation=value,
                    raw_genes=raw_genes,
                    low_expression_basemean=args.low_expression_basemean,
                )
            except Exception as exc:
                warnings.warn(
                    f"{annotation_level}={value}: DESeq2 failed: {exc}",
                    RuntimeWarning,
                    stacklevel=2,
                )
                result = pd.DataFrame(
                    [
                        empty_result_row(
                            annotation_level=annotation_level,
                            annotation=value,
                            gene=gene,
                            present=gene in raw_genes,
                            status=f"failed: {exc}",
                            n_young_samples=0,
                            n_aged_samples=0,
                            young_cell_count=0,
                            aged_cell_count=0,
                        )
                        for gene in LOX_GENES
                    ]
                )
            tables.append(result)

    complete = pd.concat(tables, ignore_index=True)
    ordered_columns = [
        "gene",
        "annotation_level",
        "cell_type_or_subtype",
        "age_comparison",
        "n_young_samples",
        "n_aged_samples",
        "young_cell_count",
        "aged_cell_count",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "significance_label",
        "direction",
        "status",
        "gene_present",
    ]
    complete = complete.loc[:, ordered_columns].sort_values(
        ["annotation_level", "cell_type_or_subtype", "gene"]
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "lox_pseudobulk_complete_results.csv"
    tsv_path = args.output_dir / "lox_pseudobulk_complete_results.tsv"
    complete.to_csv(csv_path, index=False)
    complete.to_csv(tsv_path, index=False, sep="\t")

    print(f"\nSaved CSV: {csv_path}", flush=True)
    print(f"Saved TSV: {tsv_path}", flush=True)
    print(f"Rows: {len(complete):,}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
