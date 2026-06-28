#!/usr/bin/env python
"""Build a cross-dataset LOX-family validation matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "tables" / "cross_dataset_lox_validation_matrix.tsv"
DEFAULT_SUMMARY = PROJECT_ROOT / "results" / "tables" / "cross_dataset_lox_consistency_summary.tsv"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "cross_dataset_lox_analysis.md"

CLAIM_EXPECTED = {
    "broad_FB_Lox_aged_lower": "aged-lower",
    "broad_FB_Loxl2_aged_lower": "aged-lower",
    "capsFB_Lox_aged_lower": "aged-lower",
    "medFB_Loxl1_aged_higher": "aged-higher",
    "medFB_Loxl2_aged_lower": "aged-lower",
    "epithelial_Loxl2_aged_lower": "aged-lower",
    "mTEC1_Loxl2_aged_lower": "aged-lower",
}

CLAIM_LABELS = {
    "broad_FB_Lox_aged_lower": ("FB", "Lox", ["FB", "pooled FB"]),
    "broad_FB_Loxl2_aged_lower": ("FB", "Loxl2", ["FB", "pooled FB"]),
    "capsFB_Lox_aged_lower": ("3:capsFB", "Lox", ["3:capsFB"]),
    "medFB_Loxl1_aged_higher": ("5:medFB", "Loxl1", ["5:medFB"]),
    "medFB_Loxl2_aged_lower": ("5:medFB", "Loxl2", ["5:medFB"]),
    "epithelial_Loxl2_aged_lower": ("TEC/Epi", "Loxl2", ["TEC", "Epi"]),
    "mTEC1_Loxl2_aged_lower": ("13:mTEC1", "Loxl2", ["13:mTEC1"]),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def direction_from_value(value: float) -> str:
    if pd.isna(value):
        return "not_testable"
    if value < 0:
        return "aged-lower"
    if value > 0:
        return "aged-higher"
    return "flat"


def consistency(direction: str, expected: str, confidence: str) -> int:
    if (
        confidence in {"not_testable", "inconclusive", "external_broad_context_not_subtype"}
        or direction in {"not_testable", "inconclusive"}
    ):
        return 0
    if direction == expected:
        return 1
    if direction in {"aged-lower", "aged-higher"} and expected in {"aged-lower", "aged-higher"}:
        return -1
    return 0


def base_row(**kwargs) -> dict[str, object]:
    row = {
        "dataset": "",
        "species": "",
        "data_type": "",
        "age_contrast": "",
        "biological_resolution": "",
        "cell_group": "",
        "gene": "",
        "young_or_low_age_n": np.nan,
        "aged_or_high_age_n": np.nan,
        "effect_metric": "",
        "effect_value": np.nan,
        "direction": "not_testable",
        "pvalue": np.nan,
        "confidence_level": "inconclusive",
        "claim_mapping": "",
        "notes": "",
    }
    row.update(kwargs)
    return row


def load_table(path: Path, sep: str = "\t") -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path, sep=sep)


def add_internal_rows(rows: list[dict[str, object]]) -> None:
    path = PROJECT_ROOT / "results" / "tables" / "lox_pseudobulk_complete_results.tsv"
    df = load_table(path)
    if df is None:
        return
    wanted = {
        "broad_FB_Lox_aged_lower": ("cell_type", "FB", "Lox"),
        "broad_FB_Loxl2_aged_lower": ("cell_type", "FB", "Loxl2"),
        "capsFB_Lox_aged_lower": ("subtype", "3:capsFB", "Lox"),
        "medFB_Loxl1_aged_higher": ("subtype", "5:medFB", "Loxl1"),
        "medFB_Loxl2_aged_lower": ("subtype", "5:medFB", "Loxl2"),
        "epithelial_Loxl2_aged_lower": ("cell_type", "TEC", "Loxl2"),
        "mTEC1_Loxl2_aged_lower": ("subtype", "13:mTEC1", "Loxl2"),
    }
    for claim, (level, group, gene) in wanted.items():
        hit = df.loc[
            df["annotation_level"].eq(level)
            & df["cell_type_or_subtype"].eq(group)
            & df["gene"].eq(gene)
        ]
        if hit.empty and group == "FB":
            hit = df.loc[
                df["annotation_level"].eq(level)
                & df["cell_type_or_subtype"].isin(["FB", "pooled FB"])
                & df["gene"].eq(gene)
            ]
        if hit.empty:
            rows.append(
                base_row(
                    dataset="GSE240016_internal",
                    species="Mouse",
                    data_type="single-cell pseudobulk DESeq2",
                    age_contrast="18mo vs 02mo",
                    biological_resolution=level,
                    cell_group=group,
                    gene=gene,
                    claim_mapping=claim,
                    confidence_level="not_testable",
                    notes="Requested internal row not found.",
                )
            )
            continue
        r = hit.iloc[0]
        direction = "aged-higher" if r["log2FoldChange"] > 0 else "aged-lower" if r["log2FoldChange"] < 0 else "flat"
        rows.append(
            base_row(
                dataset="GSE240016_internal",
                species="Mouse",
                data_type="single-cell pseudobulk DESeq2",
                age_contrast=str(r["age_comparison"]),
                biological_resolution=str(r["annotation_level"]),
                cell_group=str(r["cell_type_or_subtype"]),
                gene=str(r["gene"]),
                young_or_low_age_n=r["n_young_samples"],
                aged_or_high_age_n=r["n_aged_samples"],
                effect_metric="DESeq2_log2FoldChange_aged_vs_young",
                effect_value=r["log2FoldChange"],
                direction=direction,
                pvalue=r["pvalue"],
                confidence_level="internal_primary",
                claim_mapping=claim,
                notes=f"Internal observation; padj={r['padj']}; significance={r['significance_label']}.",
            )
        )


def add_gse223049_rows(rows: list[dict[str, object]]) -> None:
    path = PROJECT_ROOT / "results" / "tables" / "external_gse223049_lox_validation_summary_v2.tsv"
    if not path.exists():
        path = PROJECT_ROOT / "results" / "tables" / "external_gse223049_lox_validation_summary.tsv"
    df = load_table(path)
    if df is None:
        return
    df = df.drop_duplicates(["external_cell_type", "gene"])
    mapping = [
        ("Thymic_fibroblasts", "Lox", "broad_FB_Lox_aged_lower", "Broad sorted bulk fibroblast; not capsFB-specific."),
        ("Thymic_fibroblasts", "Loxl2", "broad_FB_Loxl2_aged_lower", "Broad sorted bulk fibroblast; not medFB-specific."),
        ("Thymic_fibroblasts", "Loxl1", "medFB_Loxl1_aged_higher", "Broad fibroblast only; cannot validate or contradict medFB specificity."),
        ("Thymic_fibroblasts", "Loxl2", "medFB_Loxl2_aged_lower", "Broad fibroblast context only; not subtype-specific validation."),
        ("Thymic_epithelial", "Loxl2", "epithelial_Loxl2_aged_lower", "Broad sorted thymic epithelial validation."),
        ("Thymic_epithelial", "Loxl2", "mTEC1_Loxl2_aged_lower", "Broad epithelial context only; cannot validate mTEC1 specificity."),
    ]
    for cell, gene, claim, note in mapping:
        hit = df.loc[df["external_cell_type"].eq(cell) & df["gene"].eq(gene)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        effect = r["delta_log2_cpm_aged_minus_young"]
        pval = r.get("nonparametric_pvalue_mann_whitney", np.nan)
        rows.append(
            base_row(
                dataset="GSE223049",
                species="Mouse",
                data_type="sorted bulk RNA-seq",
                age_contrast="22-24mo vs 2mo",
                biological_resolution="broad sorted cell type",
                cell_group=cell,
                gene=gene,
                young_or_low_age_n=r.get("young_n", np.nan),
                aged_or_high_age_n=r.get("aged_n", np.nan),
                effect_metric="delta_log2_CPM_plus1_aged_minus_young",
                effect_value=effect,
                direction=direction_from_value(effect),
                pvalue=pval,
                confidence_level=(
                    "external_broad_context_not_subtype"
                    if claim in {"medFB_Loxl1_aged_higher", "medFB_Loxl2_aged_lower", "mTEC1_Loxl2_aged_lower"}
                    else "external_broad_descriptive"
                ),
                claim_mapping=claim,
                notes=note,
            )
        )


def add_emtab8560_rows(rows: list[dict[str, object]]) -> None:
    path = PROJECT_ROOT / "results" / "tables" / "external_emtab8560_tec_lox_summary.tsv"
    df = load_table(path)
    if df is None:
        return
    mapping = [
        ("broad_TEC", "Loxl2", "epithelial_Loxl2_aged_lower", "TEC-only age series; sort-day summaries, not donor-level biological replication."),
        ("mTEC_like", "Loxl2", "mTEC1_Loxl2_aged_lower", "mTEC-like groups only; not exact mTEC1 annotation."),
        ("mTEClo", "Loxl2", "mTEC1_Loxl2_aged_lower", "mTEClo approximation; not exact mTEC1 annotation."),
        ("mTEChi", "Loxl2", "mTEC1_Loxl2_aged_lower", "mTEChi approximation; not exact mTEC1 annotation."),
    ]
    for group, gene, claim, note in mapping:
        hit = df.loc[df["analysis_group"].eq(group) & df["gene"].eq(gene)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        effect = r["mean_delta_oldest_minus_youngest"]
        rows.append(
            base_row(
                dataset="E-MTAB-8560",
                species="Mouse",
                data_type="Smart-seq2 single-cell RNA-seq",
                age_contrast=f"{int(r['oldest_age_week'])}wk vs {int(r['youngest_age_week'])}wk",
                biological_resolution="TEC sorted/inferred group",
                cell_group=group,
                gene=gene,
                young_or_low_age_n=np.nan,
                aged_or_high_age_n=np.nan,
                effect_metric="mean_delta_oldest_minus_youngest_log2_norm_plus1",
                effect_value=effect,
                direction=direction_from_value(effect),
                pvalue=r["spearman_pvalue"],
                confidence_level="external_TEC_age_series_descriptive",
                claim_mapping=claim,
                notes=note,
            )
        )


def add_gse231906_rows(rows: list[dict[str, object]]) -> None:
    path = PROJECT_ROOT / "results" / "tables" / "external_gse231906_human_lox_summary.tsv"
    df = load_table(path)
    if df is None:
        return
    mapping = [
        ("fibroblast_like_mesenchymal", "LOX", "broad_FB_Lox_aged_lower"),
        ("fibroblast_like_mesenchymal", "LOXL1", "medFB_Loxl1_aged_higher"),
        ("fibroblast_like_mesenchymal", "LOXL2", "broad_FB_Loxl2_aged_lower"),
        ("fibroblast_like_mesenchymal", "LOXL2", "medFB_Loxl2_aged_lower"),
        ("epithelial", "LOXL2", "epithelial_Loxl2_aged_lower"),
        ("mTEC_like", "LOXL2", "mTEC1_Loxl2_aged_lower"),
    ]
    for group, gene, claim in mapping:
        hit = df.loc[df["candidate_group"].eq(group)]
        n_cells = hit["n_cells"].sum() if not hit.empty else np.nan
        n_donors = hit["n_donors_or_samples"].max() if not hit.empty else np.nan
        rows.append(
            base_row(
                dataset="GSE231906",
                species="Human",
                data_type="10x scRNA-seq",
                age_contrast="continuous human age metadata; expression not parsed",
                biological_resolution="metadata candidate group",
                cell_group=group,
                gene=gene,
                young_or_low_age_n=np.nan,
                aged_or_high_age_n=np.nan,
                effect_metric="not_available_metadata_only",
                effect_value=np.nan,
                direction="not_testable",
                pvalue=np.nan,
                confidence_level="not_testable",
                claim_mapping=claim,
                notes=f"Metadata-only candidate group; expression matrix not parsed. Candidate cells={n_cells}, donors/samples={n_donors}.",
            )
        )


def add_missing_claim_rows(rows: list[dict[str, object]]) -> None:
    existing = {(r["dataset"], r["claim_mapping"]) for r in rows}
    datasets = [
        ("GSE223049", "Mouse", "sorted bulk RNA-seq", "external broad bulk"),
        ("E-MTAB-8560", "Mouse", "Smart-seq2 single-cell RNA-seq", "TEC-only"),
        ("GSE231906", "Human", "10x scRNA-seq", "metadata-only"),
    ]
    for dataset, species, dtype, resolution in datasets:
        for claim, (cell_group, gene, _) in CLAIM_LABELS.items():
            if (dataset, claim) in existing:
                continue
            reason = "not applicable to dataset compartment or only metadata available"
            rows.append(
                base_row(
                    dataset=dataset,
                    species=species,
                    data_type=dtype,
                    age_contrast="not tested",
                    biological_resolution=resolution,
                    cell_group=cell_group,
                    gene=gene,
                    confidence_level="not_testable",
                    claim_mapping=claim,
                    notes=reason,
                )
            )


def build_summary(matrix: pd.DataFrame) -> pd.DataFrame:
    df = matrix.copy()
    df["expected_direction"] = df["claim_mapping"].map(CLAIM_EXPECTED)
    df["consistency_score"] = [
        consistency(direction, expected, confidence)
        for direction, expected, confidence in zip(df["direction"], df["expected_direction"], df["confidence_level"])
    ]
    external = df.loc[~df["dataset"].eq("GSE240016_internal")]
    rows = []
    for claim, expected in CLAIM_EXPECTED.items():
        all_rows = df.loc[df["claim_mapping"].eq(claim)]
        ext_rows = external.loc[external["claim_mapping"].eq(claim)]
        supported = int((ext_rows["consistency_score"] == 1).sum())
        contradicted = int((ext_rows["consistency_score"] == -1).sum())
        not_testable = int((ext_rows["consistency_score"] == 0).sum())
        internal = all_rows.loc[all_rows["dataset"].eq("GSE240016_internal")]
        internal_direction = internal["direction"].iloc[0] if not internal.empty else "not_available"
        if contradicted:
            conclusion = "contradicted_by_external"
        elif supported and claim == "mTEC1_Loxl2_aged_lower":
            conclusion = "externally_directionally_supported_approximate_mTEC_like"
        elif supported:
            conclusion = "externally_directionally_supported"
        elif internal_direction in {"aged-lower", "aged-higher"}:
            conclusion = "only_internally_observed_or_not_testable_external"
        else:
            conclusion = "not_testable"
        rows.append(
            {
                "claim_mapping": claim,
                "expected_direction": expected,
                "internal_direction": internal_direction,
                "n_external_supporting": supported,
                "n_external_opposite": contradicted,
                "n_external_not_testable_or_inconclusive": not_testable,
                "net_external_consistency_score": int(ext_rows["consistency_score"].sum()),
                "conclusion": conclusion,
            }
        )
    return pd.DataFrame(rows), df


def write_report(matrix: pd.DataFrame, summary: pd.DataFrame, report_path: Path) -> None:
    supported = summary.loc[summary["conclusion"].eq("externally_directionally_supported")]
    approximate = summary.loc[
        summary["conclusion"].eq("externally_directionally_supported_approximate_mTEC_like")
    ]
    internal_only = summary.loc[summary["conclusion"].eq("only_internally_observed_or_not_testable_external")]
    contradicted = summary.loc[summary["conclusion"].eq("contradicted_by_external")]

    lines = [
        "# Cross-dataset LOX-family validation analysis",
        "",
        "## Scope",
        "",
        "This report combines internal GSE240016 pseudobulk results with available external outputs from GSE223049, E-MTAB-8560, and GSE231906. Broad sorted bulk evidence is kept broad and is not counted as subtype-specific validation.",
        "",
        "## Supported Across Datasets",
        "",
    ]
    if supported.empty:
        lines.append("No claim has external directional support in the available matrix.")
    else:
        for _, row in supported.iterrows():
            lines.append(
                f"- `{row['claim_mapping']}`: {row['n_external_supporting']} external exact or appropriately resolved row(s) match the expected direction; net score {row['net_external_consistency_score']}."
            )
    lines.extend(["", "## Only Internally Observed Or Externally Not Testable", ""])
    if internal_only.empty:
        lines.append("None.")
    else:
        for _, row in internal_only.iterrows():
            lines.append(
                f"- `{row['claim_mapping']}`: internal direction is `{row['internal_direction']}`, but external rows are not exact tests or are unavailable."
            )
    lines.extend(["", "## Contradicted Findings", ""])
    if contradicted.empty:
        lines.append("No claim is contradicted by the currently available external rows.")
    else:
        for _, row in contradicted.iterrows():
            lines.append(
                f"- `{row['claim_mapping']}`: {row['n_external_opposite']} external row(s) have the opposite direction."
            )
    lines.extend(
        [
            "",
            "## Approximate Directional Support",
            "",
        ]
    )
    if approximate.empty:
        lines.append("No claim is currently limited to approximate external directional support.")
    else:
        for _, row in approximate.iterrows():
            lines.append(
                f"- `{row['claim_mapping']}`: {row['n_external_supporting']} mTEC-like/TEC-subset row(s) match the expected direction, but these are not exact mTEC1 replication."
            )
    lines.extend(
        [
            "",
            "## Why Exact Subtype Validation Remains Limited",
            "",
            "- GSE223049 is sorted bulk RNA-seq with broad thymic fibroblast and thymic epithelial populations, so it cannot validate capsFB, medFB, or mTEC1 specificity.",
            "- E-MTAB-8560 is TEC-only and useful for mTEC-like `Loxl2`, but it has no fibroblast compartment and does not share the exact mTEC1 annotation.",
            "- GSE231906 currently contributes metadata-only candidate groups in this repository; expression was not parsed because the raw matrix archive is large and requires guarded donor/source matching.",
            "- Therefore, broad fibroblast and epithelial directions can be compared across datasets, while medFB and mTEC1 specificity remain internally observed or approximately tested rather than independently subtype-validated.",
            "",
            "## Matrix Outputs",
            "",
            "- `results/tables/cross_dataset_lox_validation_matrix.tsv`",
            "- `results/tables/cross_dataset_lox_consistency_summary.tsv`",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows: list[dict[str, object]] = []
    add_internal_rows(rows)
    add_gse223049_rows(rows)
    add_emtab8560_rows(rows)
    add_gse231906_rows(rows)
    add_missing_claim_rows(rows)
    matrix = pd.DataFrame(rows)
    summary, matrix = build_summary(matrix)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t", index=False)
    summary.to_csv(args.summary_output, sep="\t", index=False)
    write_report(matrix, summary, args.report)
    print(f"Saved matrix: {args.output}", flush=True)
    print(f"Saved summary: {args.summary_output}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
