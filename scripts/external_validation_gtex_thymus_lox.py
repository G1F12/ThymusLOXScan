#!/usr/bin/env python
"""Guarded GTEx thymus LOX-family age-context analysis.

The script uses local GTEx files only. It does not download full GTEx archives.
If donor-level thymus expression and age metadata are not available locally, it
writes a feasibility report instead of attempting an age analysis.
"""

from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

try:
    import statsmodels.formula.api as smf
except ImportError:  # pragma: no cover - optional fallback
    smf = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEARCH_DIRS = [
    PROJECT_ROOT / "data" / "external" / "GTEx",
    PROJECT_ROOT / "data" / "external" / "gtex",
]
DEFAULT_BY_SAMPLE = PROJECT_ROOT / "results" / "tables" / "external_gtex_thymus_lox_by_sample.tsv"
DEFAULT_AGE_SUMMARY = PROJECT_ROOT / "results" / "tables" / "external_gtex_thymus_lox_age_summary.tsv"
DEFAULT_MODELS = PROJECT_ROOT / "results" / "tables" / "external_gtex_thymus_lox_age_models.tsv"
DEFAULT_FEASIBILITY = PROJECT_ROOT / "results" / "tables" / "external_gtex_thymus_lox_feasibility.tsv"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "gtex_thymus_lox_age_analysis.md"
DEFAULT_FEASIBILITY_REPORT = PROJECT_ROOT / "reports" / "gtex_thymus_lox_age_feasibility.md"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "external_validation" / "gtex"

LOX_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]
CLASS_NOT_POSSIBLE = "age analysis not possible from available public/local files"
CLASS_FAILED = "analysis failed due to missing/invalid metadata"


@dataclass
class LocalFiles:
    expression: Path | None = None
    sample_attributes: Path | None = None
    subject_phenotypes: Path | None = None
    tissue_annotation: Path | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", action="append", type=Path, default=[])
    parser.add_argument("--by-sample-output", type=Path, default=DEFAULT_BY_SAMPLE)
    parser.add_argument("--age-summary-output", type=Path, default=DEFAULT_AGE_SUMMARY)
    parser.add_argument("--models-output", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--feasibility-output", type=Path, default=DEFAULT_FEASIBILITY)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--feasibility-report", type=Path, default=DEFAULT_FEASIBILITY_REPORT)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz" or path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def delimiter(path: Path) -> str:
    with open_text(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                return "\t" if line.count("\t") >= line.count(",") else ","
    return "\t"


def candidate_files(search_dirs: list[Path]) -> list[Path]:
    files: list[Path] = []
    for folder in search_dirs:
        if folder.exists():
            for path in folder.rglob("*"):
                suffixes = [suffix.lower() for suffix in path.suffixes]
                effective_suffix = suffixes[-2] if suffixes[-1:] == [".gz"] and len(suffixes) >= 2 else path.suffix.lower()
                if path.is_file() and effective_suffix in {".txt", ".tsv", ".csv", ".gct"}:
                    files.append(path)
    return sorted(files)


def pick_files(files: list[Path]) -> LocalFiles:
    picked = LocalFiles()
    for path in files:
        name = path.name.lower()
        if picked.sample_attributes is None and (
            "sampleattributes" in name or "sample_attributes" in name or "sample-attributes" in name
        ):
            picked.sample_attributes = path
        elif picked.subject_phenotypes is None and (
            "subjectphenotypes" in name or "subject_phenotypes" in name or "phenotype" in name
        ):
            picked.subject_phenotypes = path
        elif picked.tissue_annotation is None and "tissue" in name and "annotation" in name:
            picked.tissue_annotation = path
        elif picked.expression is None and (
            "gene_tpm" in name
            or "gene_tpm" in str(path).lower()
            or ("tpm" in name and "gene" in name)
            or ("expression" in name and "gene" in name)
            or name.endswith(".gct")
        ):
            picked.expression = path
    return picked


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    sep = delimiter(path)
    skiprows = 2 if path.name.endswith(".gct") or path.name.endswith(".gct.gz") else 0
    return pd.read_csv(path, sep=sep, skiprows=skiprows, low_memory=False, **kwargs)


def subject_from_sample(sample_id: str) -> str:
    parts = str(sample_id).split("-")
    return "-".join(parts[:2]) if len(parts) >= 2 else str(sample_id)


def midpoint_age(value: object) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value)
    nums = re.findall(r"\d+(?:\.\d+)?", text)
    if len(nums) >= 2:
        return (float(nums[0]) + float(nums[1])) / 2
    if len(nums) == 1:
        return float(nums[0])
    return np.nan


def ordered_age_label(value: object) -> str:
    if pd.isna(value):
        return "unknown"
    return str(value).strip()


def find_column(columns: list[str], candidates: list[str]) -> str | None:
    lower = {col.lower(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def thymus_samples(sample_attributes: pd.DataFrame) -> pd.DataFrame:
    cols = list(sample_attributes.columns)
    sample_col = find_column(cols, ["SAMPID", "sample_id", "sample", "run_s"])
    tissue_cols = [col for col in cols if col.lower() in {"smtsd", "smts", "tissue", "tissue_site_detail"}]
    if sample_col is None or not tissue_cols:
        raise ValueError("Sample attributes must include a sample ID and tissue column.")
    tissue_text = sample_attributes[tissue_cols].astype(str).agg(" ".join, axis=1).str.lower()
    out = sample_attributes.loc[tissue_text.str.contains("thymus", na=False)].copy()
    out = out.rename(columns={sample_col: "sample_id"})
    out["subject_id"] = out["sample_id"].map(subject_from_sample)
    return out


def annotate_samples(thymus: pd.DataFrame, phenotypes: pd.DataFrame) -> pd.DataFrame:
    cols = list(phenotypes.columns)
    subject_col = find_column(cols, ["SUBJID", "subject_id", "donor_id"])
    age_col = find_column(cols, ["AGE", "age"])
    sex_col = find_column(cols, ["SEX", "sex"])
    if subject_col is None or age_col is None:
        raise ValueError("Subject phenotype table must include subject ID and age columns.")
    pheno = phenotypes.rename(columns={subject_col: "subject_id", age_col: "age"})
    if sex_col:
        pheno = pheno.rename(columns={sex_col: "sex"})
    else:
        pheno["sex"] = pd.NA
    keep = ["subject_id", "age", "sex"]
    merged = thymus.merge(pheno[keep], on="subject_id", how="left")
    merged["age_bin"] = merged["age"].map(ordered_age_label)
    merged["age_midpoint"] = merged["age"].map(midpoint_age)

    rin_col = find_column(list(merged.columns), ["SMRIN", "RIN", "rna_integrity_number"])
    ischemic_col = find_column(list(merged.columns), ["SMTSISCH", "TRISCHD", "ischemic_time"])
    batch_col = find_column(list(merged.columns), ["SMGEBTCH", "SMCENTER", "batch", "platform"])
    merged["rin"] = pd.to_numeric(merged[rin_col], errors="coerce") if rin_col else np.nan
    merged["ischemic_time"] = pd.to_numeric(merged[ischemic_col], errors="coerce") if ischemic_col else np.nan
    merged["batch_or_platform"] = merged[batch_col].astype(str) if batch_col else pd.NA
    return merged


def expression_for_genes(expression_path: Path, sample_ids: list[str]) -> pd.DataFrame:
    header = read_table(expression_path, nrows=0)
    cols = list(header.columns)
    gene_id_col = find_column(cols, ["Name", "gene_id", "Geneid", "gene"])
    symbol_col = find_column(cols, ["Description", "gene_name", "gene_symbol", "symbol"])
    sample_cols = [sample for sample in sample_ids if sample in cols]
    if not sample_cols:
        raise ValueError("No thymus sample IDs from metadata were found in the expression matrix columns.")

    usecols = []
    if gene_id_col:
        usecols.append(gene_id_col)
    if symbol_col and symbol_col not in usecols:
        usecols.append(symbol_col)
    usecols.extend(sample_cols)
    expr = read_table(expression_path, usecols=usecols)

    if symbol_col:
        expr["_symbol"] = expr[symbol_col].astype(str).str.upper()
    elif gene_id_col:
        expr["_symbol"] = expr[gene_id_col].astype(str).str.split(".").str[0].str.upper()
    else:
        first = expr.columns[0]
        expr["_symbol"] = expr[first].astype(str).str.upper()

    subset = expr.loc[expr["_symbol"].isin(LOX_GENES), ["_symbol", *sample_cols]].copy()
    if subset.empty:
        raise ValueError("LOX-family gene symbols were not found in the expression matrix.")
    subset = subset.groupby("_symbol", as_index=False)[sample_cols].mean(numeric_only=True)
    long = subset.melt(id_vars="_symbol", var_name="sample_id", value_name="tpm")
    long = long.rename(columns={"_symbol": "gene"})
    long["tpm"] = pd.to_numeric(long["tpm"], errors="coerce")
    long["log2_tpm_plus1"] = np.log2(long["tpm"] + 1)
    return long


def build_by_sample(local: LocalFiles) -> pd.DataFrame:
    assert local.expression and local.sample_attributes and local.subject_phenotypes
    samples = thymus_samples(read_table(local.sample_attributes))
    annotated = annotate_samples(samples, read_table(local.subject_phenotypes))
    if annotated.empty:
        raise ValueError("No thymus samples found in local GTEx sample attributes.")
    if annotated["age_midpoint"].notna().sum() == 0:
        raise ValueError("Thymus samples were found, but no usable donor age metadata was available.")
    expr = expression_for_genes(local.expression, annotated["sample_id"].astype(str).tolist())
    cols = [
        "sample_id",
        "subject_id",
        "age",
        "age_bin",
        "age_midpoint",
        "sex",
        "rin",
        "ischemic_time",
        "batch_or_platform",
    ]
    return expr.merge(annotated[cols], on="sample_id", how="left")


def summarize_age(by_sample: pd.DataFrame) -> pd.DataFrame:
    return (
        by_sample.groupby(["gene", "age_bin"], dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            n_subjects=("subject_id", "nunique"),
            mean_age_midpoint=("age_midpoint", "mean"),
            mean_log2_tpm_plus1=("log2_tpm_plus1", "mean"),
            median_log2_tpm_plus1=("log2_tpm_plus1", "median"),
            sd_log2_tpm_plus1=("log2_tpm_plus1", "std"),
            mean_tpm=("tpm", "mean"),
            detected_fraction=("tpm", lambda x: float((pd.to_numeric(x, errors="coerce") > 0).mean())),
        )
        .reset_index()
        .sort_values(["gene", "mean_age_midpoint"])
    )


def fit_models(by_sample: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, df in by_sample.groupby("gene"):
        df = df.dropna(subset=["age_midpoint", "log2_tpm_plus1"]).copy()
        if len(df) < 4 or df["age_midpoint"].nunique() < 2:
            rows.append(
                {
                    "gene": gene,
                    "model": "not_fit",
                    "age_term": "age_midpoint",
                    "age_estimate": np.nan,
                    "age_pvalue": np.nan,
                    "spearman_rho": np.nan,
                    "spearman_pvalue": np.nan,
                    "n_samples": len(df),
                    "note": "Insufficient samples or age variation.",
                }
            )
            continue
        rho, rho_p = stats.spearmanr(df["age_midpoint"], df["log2_tpm_plus1"])
        note = "Spearman trend plus linear model where available."
        estimate = pvalue = np.nan
        model_name = "spearman_only"
        if smf is not None:
            terms = ["age_midpoint"]
            if df["sex"].notna().nunique() > 1:
                terms.append("C(sex)")
            if df["rin"].notna().sum() >= max(4, int(len(df) * 0.5)):
                terms.append("rin")
            if df["ischemic_time"].notna().sum() >= max(4, int(len(df) * 0.5)):
                terms.append("ischemic_time")
            if df["batch_or_platform"].notna().nunique() > 1 and len(df) >= 10:
                terms.append("C(batch_or_platform)")
            formula = "log2_tpm_plus1 ~ " + " + ".join(terms)
            try:
                result = smf.ols(formula, data=df).fit()
                estimate = float(result.params.get("age_midpoint", np.nan))
                pvalue = float(result.pvalues.get("age_midpoint", np.nan))
                model_name = formula
            except Exception as exc:  # noqa: BLE001 - report exact issue without failing all genes
                note = f"OLS failed; Spearman retained. Reason: {type(exc).__name__}: {exc}"
        rows.append(
            {
                "gene": gene,
                "model": model_name,
                "age_term": "age_midpoint",
                "age_estimate": estimate,
                "age_pvalue": pvalue,
                "spearman_rho": float(rho),
                "spearman_pvalue": float(rho_p),
                "n_samples": len(df),
                "note": note,
            }
        )
    return pd.DataFrame(rows)


def classify(models: pd.DataFrame) -> str:
    hit = models.loc[models["gene"].eq("LOXL2")]
    if hit.empty:
        return CLASS_FAILED
    row = hit.iloc[0]
    estimate = row["age_estimate"]
    pvalue = row["age_pvalue"]
    if pd.isna(estimate):
        estimate = row["spearman_rho"]
        pvalue = row["spearman_pvalue"]
    if pd.isna(estimate) or pd.isna(pvalue) or pvalue >= 0.05:
        return "age analysis completed: LOXL2 no clear age trend"
    if estimate < 0:
        return "age analysis completed: LOXL2 aged-lower"
    return "age analysis completed: LOXL2 aged-higher"


def save_figures(by_sample: pd.DataFrame, figure_dir: Path) -> list[Path]:
    sns.set_theme(style="whitegrid", context="paper")
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    lo = by_sample.loc[by_sample["gene"].eq("LOXL2")].copy()
    lo = lo.sort_values("age_midpoint")
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    sns.scatterplot(data=lo, x="age_midpoint", y="log2_tpm_plus1", color="#0072B2", s=42, ax=ax)
    if lo["age_midpoint"].nunique() >= 2:
        sns.regplot(
            data=lo,
            x="age_midpoint",
            y="log2_tpm_plus1",
            scatter=False,
            color="0.25",
            ci=None,
            ax=ax,
        )
    ax.set_xlabel("Donor age midpoint from GTEx age bin (years)")
    ax.set_ylabel("LOXL2 log2(TPM + 1)")
    ax.set_title("GTEx whole-thymus LOXL2 by donor age bin")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        path = figure_dir / f"gtex_thymus_loxl2_by_age.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    sns.lineplot(
        data=by_sample.sort_values("age_midpoint"),
        x="age_midpoint",
        y="log2_tpm_plus1",
        hue="gene",
        estimator="mean",
        errorbar=None,
        marker="o",
        ax=ax,
    )
    ax.set_xlabel("Donor age midpoint from GTEx age bin (years)")
    ax.set_ylabel("Mean log2(TPM + 1)")
    ax.set_title("GTEx whole-thymus LOX-family age trends")
    ax.legend(title="", frameon=False)
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        path = figure_dir / f"gtex_thymus_lox_family_age_trends.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def feasibility_rows(search_dirs: list[Path], files: list[Path], local: LocalFiles, reason: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "requirement": "local_search_paths",
                "status": "checked",
                "detail": "; ".join(str(path) for path in search_dirs),
            },
            {
                "requirement": "candidate_files_found",
                "status": str(len(files)),
                "detail": "; ".join(str(path.relative_to(PROJECT_ROOT)) for path in files[:25]),
            },
            {
                "requirement": "gene_tpm_expression_matrix",
                "status": "present" if local.expression else "missing",
                "detail": str(local.expression or ""),
            },
            {
                "requirement": "sample_attributes",
                "status": "present" if local.sample_attributes else "missing",
                "detail": str(local.sample_attributes or ""),
            },
            {
                "requirement": "subject_phenotype_age_metadata",
                "status": "present" if local.subject_phenotypes else "missing",
                "detail": str(local.subject_phenotypes or ""),
            },
            {
                "requirement": "classification",
                "status": CLASS_NOT_POSSIBLE if reason == "missing_files" else CLASS_FAILED,
                "detail": reason,
            },
        ]
    )


def write_feasibility_report(path: Path, table: pd.DataFrame, reason: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# GTEx thymus LOX-family age feasibility",
        "",
        f"Classification: {CLASS_NOT_POSSIBLE if reason == 'missing_files' else CLASS_FAILED}",
        "",
        "GTEx remains a planned human whole-tissue analysis requiring donor-level expression and age metadata.",
        "",
        "## What was checked",
        "",
        "| requirement | status | detail |",
        "|---|---|---|",
    ]
    for _, row in table.iterrows():
        detail = str(row["detail"]).replace("|", "\\|")
        lines.append(f"| {row['requirement']} | {row['status']} | {detail} |")
    lines.extend(
        [
            "",
            "## Manual download instructions",
            "",
            "Place the GTEx v8 gene TPM matrix and annotation files under `data/external/GTEx/` or `data/external/gtex/`.",
            "",
            "Expected local files include:",
            "",
            "- `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` or another gene-level TPM/expression matrix.",
            "- `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`.",
            "- `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt`.",
            "",
            "Do not commit the full GTEx matrix or other large third-party data files to this repository.",
            "",
            "## Methods-ready note",
            "",
            "If donor-level files are added locally, this analysis will subset GTEx whole-thymus samples, summarize LOX, LOXL1, LOXL2, LOXL3, and LOXL4 as log2(TPM+1), and fit descriptive age associations with available covariates. GTEx thymus is whole tissue, not TEC/mTEC subtype-resolved, and cannot validate mouse mTEC1.",
            "",
            "## Discussion-ready note",
            "",
            "At this stage, GTEx provides no completed age-context result for this repository. Any future GTEx result should be described as broad human thymus whole-tissue directional context only, and only if donor-level expression and age metadata support that description.",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_analysis_report(path: Path, classification: str, models: pd.DataFrame, figures: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lo = models.loc[models["gene"].eq("LOXL2")].iloc[0]
    if classification.endswith("aged-lower"):
        interpretation = "LOXL2 shows broad human thymus whole-tissue directional context, not validation."
    elif classification.endswith("aged-higher") or classification.endswith("no clear age trend"):
        interpretation = "No supportive broad GTEx whole-tissue age trend was observed."
    else:
        interpretation = "The GTEx result could not be classified from the available files."
    lines = [
        "# GTEx thymus LOX-family age analysis",
        "",
        f"Classification: {classification}",
        "",
        "## Scope",
        "",
        "GTEx thymus is whole tissue and is not TEC, mTEC, fibroblast, or mTEC1 subtype-resolved. Donor age, tissue composition, RNA quality, ischemic time, sex, and technical covariates may confound descriptive age trends.",
        "",
        "## LOXL2 result",
        "",
        f"- Age estimate: {lo['age_estimate']}",
        f"- Age p-value: {lo['age_pvalue']}",
        f"- Spearman rho: {lo['spearman_rho']}",
        f"- Spearman p-value: {lo['spearman_pvalue']}",
        f"- Interpretation: {interpretation}",
        "",
        "## Methods-ready note",
        "",
        "Local GTEx v8-style files were parsed for whole-thymus samples. LOX-family gene expression was analyzed as log2(TPM+1). Age was modeled from GTEx donor age bins using bin midpoints, with available sex, RIN, ischemic time, and batch/platform covariates included when present.",
        "",
        "## Discussion-ready note",
        "",
        "This analysis can provide broad human thymus whole-tissue age context, but it does not validate mouse mTEC1, does not prove human conservation, and does not address protein-level LOXL2.",
        "",
        "## Figures",
        "",
    ]
    lines.extend(f"- `{path.relative_to(PROJECT_ROOT).as_posix()}`" for path in figures)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    search_dirs = [*DEFAULT_SEARCH_DIRS, *args.data_dir]
    files = candidate_files(search_dirs)
    local = pick_files(files)
    required = [local.expression, local.sample_attributes, local.subject_phenotypes]
    if any(path is None for path in required):
        table = feasibility_rows(search_dirs, files, local, "missing_files")
        args.feasibility_output.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.feasibility_output, sep="\t", index=False)
        write_feasibility_report(args.feasibility_report, table, "missing_files")
        print(f"Saved feasibility table: {args.feasibility_output}", flush=True)
        print(f"Saved feasibility report: {args.feasibility_report}", flush=True)
        return 0

    try:
        by_sample = build_by_sample(local)
        age_summary = summarize_age(by_sample)
        models = fit_models(by_sample)
        classification = classify(models)
        figures = save_figures(by_sample, args.figure_dir)
        args.by_sample_output.parent.mkdir(parents=True, exist_ok=True)
        by_sample.to_csv(args.by_sample_output, sep="\t", index=False)
        age_summary.to_csv(args.age_summary_output, sep="\t", index=False)
        models.to_csv(args.models_output, sep="\t", index=False)
        write_analysis_report(args.report, classification, models, figures)
        print(f"Saved by-sample table: {args.by_sample_output}", flush=True)
        print(f"Saved age summary: {args.age_summary_output}", flush=True)
        print(f"Saved models: {args.models_output}", flush=True)
        print(f"Saved report: {args.report}", flush=True)
        return 0
    except Exception as exc:  # noqa: BLE001 - failure report should capture exact issue
        reason = f"{type(exc).__name__}: {exc}"
        table = feasibility_rows(search_dirs, files, local, reason)
        args.feasibility_output.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.feasibility_output, sep="\t", index=False)
        write_feasibility_report(args.feasibility_report, table, reason)
        print(f"Analysis failed; saved feasibility report: {reason}", flush=True)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
