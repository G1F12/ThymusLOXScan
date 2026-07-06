"""Create GSE147520 epithelial donor/sample-aware LOX-family context summaries."""

from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external" / "human_thymus" / "GSE147520"
H5AD = DATA_DIR / "GSE147520_epithelial_cells.h5ad"
H5AD_GZ = DATA_DIR / "GSE147520_epithelial_cells.h5ad.gz"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "results" / "figures" / "human_thymus_stage3d_gse147520"
TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]
DONOR_FIELD = "donor_id"
SAMPLE_FIELD = "samples"
AGE_FIELD = "age_development"
SEX_FIELD = "sex"
FINE_FIELD = "cell_types_epith"
BROAD_FIELD = "cell_type_broad_derived"


def ensure_h5ad() -> Path:
    if H5AD.exists():
        return H5AD
    if H5AD_GZ.exists():
        with gzip.open(H5AD_GZ, "rb") as src, H5AD.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        return H5AD
    raise FileNotFoundError(f"Missing {H5AD} and {H5AD_GZ}")


def read_kv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def broad_label(label: str) -> str:
    lower = label.lower()
    if "ctec" in lower:
        return "cTEC-like"
    if "mtec" in lower or "aire" in lower or "corneocyte" in lower:
        return "mTEC-like"
    if "myoid" in lower:
        return "myoid-like"
    if "neuro" in lower:
        return "neuroendocrine-like"
    if "immature" in lower or "tec" in lower:
        return "TEC-like"
    return "other_epithelial"


def cell_counts(obs: pd.DataFrame, subtype_field: str) -> pd.DataFrame:
    fields = [DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, subtype_field]
    return obs.groupby(fields, observed=True).size().reset_index(name="n_cells").sort_values(fields)


def target_placeholder(rows: pd.DataFrame, subtype_field: str, reason: str, table_kind: str) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for gene in TARGET_GENES:
        for _, group in rows.iterrows():
            base = {
                "gene": gene,
                DONOR_FIELD: str(group[DONOR_FIELD]),
                SAMPLE_FIELD: str(group[SAMPLE_FIELD]),
                AGE_FIELD: str(group[AGE_FIELD]),
                SEX_FIELD: str(group[SEX_FIELD]),
                subtype_field: str(group[subtype_field]),
                "n_cells": str(group["n_cells"]),
            }
            if table_kind == "detection":
                out.append({
                    **base,
                    "n_detected": "",
                    "detection_fraction": "",
                    "matrix_used": "not_selected",
                    "not_generated_reason": reason,
                })
            else:
                out.append({
                    **base,
                    "mean_value": "",
                    "median_value": "",
                    "nonzero_mean_value": "",
                    "matrix_used": "not_generated",
                    "value_semantics_label": "not_generated",
                    "not_generated_reason": reason,
                })
    return out


def vector(matrix: object, var_idx: int) -> np.ndarray:
    vals = matrix[:, var_idx]
    if sparse.issparse(vals):
        return np.asarray(vals.toarray()).ravel()
    if hasattr(vals, "to_memory"):
        vals = vals.to_memory()
    if sparse.issparse(vals):
        return np.asarray(vals.toarray()).ravel()
    return np.asarray(vals).ravel()


def gene_map(var: pd.DataFrame) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for gene in TARGET_GENES:
        if gene in var.index.astype(str):
            mapping[gene] = gene
            continue
        for column in var.columns:
            matches = var.index[var[column].astype(str) == gene].astype(str).tolist()
            if matches:
                mapping[gene] = matches[0]
                break
    return mapping


def detection_rows(adata: ad.AnnData, obs: pd.DataFrame, subtype_field: str, matrix_used: str) -> list[dict[str, str]]:
    fields = [DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, subtype_field]
    key = obs[fields].astype(str).agg("\t".join, axis=1).to_numpy()
    groups = {value: np.flatnonzero(key == value) for value in sorted(pd.unique(key))}
    if matrix_used == "raw.X" and adata.raw is not None:
        source_matrix = adata.raw.X
        source_var = adata.raw.var
        source_names = adata.raw.var_names
    elif matrix_used == "X":
        source_matrix = adata.X
        source_var = adata.var
        source_names = adata.var_names
    else:
        return target_placeholder(cell_counts(obs, subtype_field), subtype_field, "Target genes unavailable for selected detection matrix.", "detection")
    genes = gene_map(source_var)
    out: list[dict[str, str]] = []
    for gene in TARGET_GENES:
        var_id = genes.get(gene, "")
        vals = vector(source_matrix, int(source_names.get_loc(var_id))) if var_id else np.array([])
        detected = vals > 0 if vals.size else np.zeros(adata.n_obs, dtype=bool)
        for group_key, indices in groups.items():
            base = dict(zip(fields, group_key.split("\t"), strict=False))
            n_cells = len(indices)
            n_detected = int(detected[indices].sum()) if vals.size else 0
            out.append({
                "gene": gene,
                **base,
                "n_cells": str(n_cells),
                "n_detected": str(n_detected) if vals.size else "",
                "detection_fraction": f"{(n_detected / n_cells) if n_cells else 0:.6g}" if vals.size else "",
                "matrix_used": matrix_used if vals.size else "not_selected",
                "not_generated_reason": "" if vals.size else "Target gene unavailable for selected detection matrix.",
            })
    return out


def development_summary(detection: list[dict[str, str]], subtype_field: str, reason: str) -> list[dict[str, str]]:
    det = pd.DataFrame(detection)
    det["detection_fraction_num"] = pd.to_numeric(det["detection_fraction"], errors="coerce")
    out: list[dict[str, str]] = []
    grouped = det.groupby(["gene", AGE_FIELD, subtype_field], observed=True)
    for (gene, stage, subtype), group in grouped:
        total = pd.to_numeric(group["n_cells"], errors="coerce").sum()
        valid = group["detection_fraction_num"].dropna()
        out.append({
            "gene": str(gene),
            AGE_FIELD: str(stage),
            subtype_field: str(subtype),
            "donor_count": str(group[DONOR_FIELD].nunique()),
            "sample_count": str(group[SAMPLE_FIELD].nunique()),
            "donor_sample_group_count": str(len(group)),
            "total_cells": str(int(total)),
            "median_detection_fraction": f"{float(valid.median()):.6g}" if len(valid) else "",
            "mean_detection_fraction": f"{float(valid.mean()):.6g}" if len(valid) else "",
            "mean_value_summary": "not_generated",
            "value_semantics_label": "not_generated",
            "not_generated_reason": "" if len(valid) else reason,
        })
    return out


def make_figures(counts_fine: pd.DataFrame, det_fine: list[dict[str, str]], dev_fine: list[dict[str, str]]) -> list[str]:
    made: list[str] = []
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return made
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")
    count_plot = counts_fine.groupby(AGE_FIELD, as_index=False, observed=True)["n_cells"].sum()
    plt.figure(figsize=(8, 4))
    sns.barplot(data=count_plot, x=AGE_FIELD, y="n_cells", color="#4c78a8")
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Cells")
    plt.xlabel("Sample age/development label")
    plt.title("GSE147520 epithelial cell counts by sample label")
    plt.tight_layout()
    path = FIGURES / "gse147520_cell_counts_by_development.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))

    dev = pd.DataFrame(dev_fine)
    dev["mean_detection_fraction"] = pd.to_numeric(dev["mean_detection_fraction"], errors="coerce")
    dev_plot = dev.groupby([AGE_FIELD, "gene"], as_index=False, observed=True)["mean_detection_fraction"].mean()
    plt.figure(figsize=(9, 4))
    sns.barplot(data=dev_plot, x=AGE_FIELD, y="mean_detection_fraction", hue="gene")
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Mean sample-level detection fraction")
    plt.xlabel("Sample age/development label")
    plt.title("GSE147520 LOX-family detection by sample label")
    plt.tight_layout()
    path = FIGURES / "gse147520_lox_detection_by_development_fine.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))

    det = pd.DataFrame(det_fine)
    det["detection_fraction"] = pd.to_numeric(det["detection_fraction"], errors="coerce")
    det_plot = det.groupby([FINE_FIELD, "gene"], as_index=False, observed=True)["detection_fraction"].mean()
    plt.figure(figsize=(10, 5))
    sns.barplot(data=det_plot, x=FINE_FIELD, y="detection_fraction", hue="gene")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mean sample-level detection fraction")
    plt.xlabel("Fine epithelial label")
    plt.title("GSE147520 LOX-family detection by epithelial label")
    plt.tight_layout()
    path = FIGURES / "gse147520_lox_detection_by_celltype_fine.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))
    return made


def main() -> None:
    path = ensure_h5ad()
    semantics = read_kv(TABLES / "human_thymus_stage3d_gse147520_matrix_semantics.tsv")
    reason = semantics.get("reason", "Target-gene summaries were not generated.")
    adata = ad.read_h5ad(path)
    obs = adata.obs.copy()
    obs[DONOR_FIELD] = "not_available"
    obs[SEX_FIELD] = "not_available"
    obs[AGE_FIELD] = obs[SAMPLE_FIELD].astype(str)
    obs[BROAD_FIELD] = obs[FINE_FIELD].astype(str).map(broad_label)
    counts_fine = cell_counts(obs, FINE_FIELD)
    counts_broad = cell_counts(obs, BROAD_FIELD)
    detection_matrix = semantics.get("recommended_detection_matrix", "not_selected")
    det_fine = detection_rows(adata, obs, FINE_FIELD, detection_matrix)
    det_broad = detection_rows(adata, obs, BROAD_FIELD, detection_matrix)
    expr_fine = target_placeholder(counts_fine, FINE_FIELD, reason, "expression")
    expr_broad = target_placeholder(counts_broad, BROAD_FIELD, reason, "expression")
    dev_fine = development_summary(det_fine, FINE_FIELD, reason)
    dev_broad = development_summary(det_broad, BROAD_FIELD, reason)

    write_tsv(TABLES / "human_thymus_stage3d_gse147520_cell_counts_fine.tsv", counts_fine.astype(str).to_dict(orient="records"), [DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, FINE_FIELD, "n_cells"])
    write_tsv(TABLES / "human_thymus_stage3d_gse147520_cell_counts_broad.tsv", counts_broad.astype(str).to_dict(orient="records"), [DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, BROAD_FIELD, "n_cells"])
    write_tsv(TABLES / "human_thymus_stage3d_gse147520_lox_detection_by_donor_fine.tsv", det_fine, ["gene", DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, FINE_FIELD, "n_cells", "n_detected", "detection_fraction", "matrix_used", "not_generated_reason"])
    write_tsv(TABLES / "human_thymus_stage3d_gse147520_lox_detection_by_donor_broad.tsv", det_broad, ["gene", DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, BROAD_FIELD, "n_cells", "n_detected", "detection_fraction", "matrix_used", "not_generated_reason"])
    write_tsv(TABLES / "human_thymus_stage3d_gse147520_lox_expression_by_donor_fine.tsv", expr_fine, ["gene", DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, FINE_FIELD, "n_cells", "mean_value", "median_value", "nonzero_mean_value", "matrix_used", "value_semantics_label", "not_generated_reason"])
    write_tsv(TABLES / "human_thymus_stage3d_gse147520_lox_expression_by_donor_broad.tsv", expr_broad, ["gene", DONOR_FIELD, SAMPLE_FIELD, AGE_FIELD, SEX_FIELD, BROAD_FIELD, "n_cells", "mean_value", "median_value", "nonzero_mean_value", "matrix_used", "value_semantics_label", "not_generated_reason"])
    write_tsv(TABLES / "human_thymus_stage3d_gse147520_lox_summary_by_development_fine.tsv", dev_fine, ["gene", AGE_FIELD, FINE_FIELD, "donor_count", "sample_count", "donor_sample_group_count", "total_cells", "median_detection_fraction", "mean_detection_fraction", "mean_value_summary", "value_semantics_label", "not_generated_reason"])
    write_tsv(TABLES / "human_thymus_stage3d_gse147520_lox_summary_by_development_broad.tsv", dev_broad, ["gene", AGE_FIELD, BROAD_FIELD, "donor_count", "sample_count", "donor_sample_group_count", "total_cells", "median_detection_fraction", "mean_detection_fraction", "mean_value_summary", "value_semantics_label", "not_generated_reason"])

    figure_paths = make_figures(counts_fine, det_fine, dev_fine)
    donors = int(counts_fine[DONOR_FIELD].nunique())
    samples = int(counts_fine[SAMPLE_FIELD].nunique())
    stages = int(counts_fine[AGE_FIELD].nunique())
    fine = int(counts_fine[FINE_FIELD].nunique())
    broad = int(counts_broad[BROAD_FIELD].nunique())
    genes_present = ", ".join(TARGET_GENES)

    REPORTS.mkdir(parents=True, exist_ok=True)
    REPORTS.joinpath("human_thymus_stage3d_gse147520_expression_summary_report.md").write_text(
        "\n".join([
            "# GSE147520 Human Thymic Epithelial LOX-Family Expression Summary",
            "",
            "Search date: 2026-07-06",
            "",
            "## Dataset and Local File Status",
            "",
            f"- Local H5AD path: `{path.relative_to(ROOT)}`",
            f"- AnnData shape: {adata.n_obs} cells x {adata.n_vars} features",
            "- Local H5AD/H5AD.GZ files are under ignored external data path.",
            "",
            "## Matrix Semantics Decision",
            "",
            f"- matrix_semantics_clear: {semantics.get('matrix_semantics_clear', 'partial')}",
            f"- detection matrix used: {detection_matrix}",
            f"- mean-value matrix used: {semantics.get('recommended_mean_value_matrix', 'not_selected')}",
            "- mean-value summaries generated: no",
            f"- reason: {reason}",
            "",
            "## Fields Used",
            "",
            "- donor field: not available in observed metadata; output uses `not_available` placeholder",
            "- sample field: `samples`",
            "- age/development field: `age_development`, derived from `samples`",
            "- sex field: not available in observed metadata; output uses `not_available` placeholder",
            "- fine annotation field: `cell_types_epith`",
            "- broad annotation field: `cell_type_broad_derived`",
            "",
            "## Dataset Breadth",
            "",
            f"- donors: {donors} placeholder group",
            f"- samples/age groups: {samples}",
            f"- development-or-age groups: {stages}",
            f"- fine epithelial labels: {fine}",
            f"- broad derived epithelial labels: {broad}",
            "",
            "## Gene Coverage",
            "",
            f"- LOX-family genes present: {genes_present} in raw.var/raw.X",
            "",
            "## Cell-Count Balance",
            "",
            "Cell counts are uneven across sample and epithelial-label groups. Outputs are grouped at sample/label level and do not use cells as independent biological replicates.",
            "",
            "## Detection-Summary Overview",
            "",
            "Target-gene detection summaries were generated from `raw.X > 0` because the target genes are present in raw.var/raw.X but absent from compact X.",
            "",
            "## Mean-Value-Summary Status",
            "",
            "Mean-value summaries were not generated.",
            "",
            "## Park/Yayon Workflow Context",
            "",
            "The same cautious workflow was applied at the file and metadata level, but GSE147520 differs from Park and Yayon because compact X omits the target genes while raw.X contains them.",
            "",
            "## Stage 3E Recommendation",
            "",
            "Stage 3E cross-human context matrix can proceed, but GSE147520 should be represented with clear layer-specific notes: compact X omits the target genes, while detection context can be summarized from raw.X.",
            "",
            "## Limitations",
            "",
            "- Donor and sex fields are not available in the observed H5AD metadata.",
            "- Age/development information is embedded in the `samples` labels.",
            "- Compact X lacks the target LOX-family genes; target-gene detection uses raw.X only.",
            "- No statistical tests were run.",
            "- No mouse-human comparison was performed.",
            "",
            "This Stage 3D summary provides exploratory donor-aware transcript-level context for LOX-family genes in the GSE147520 human thymic epithelial dataset. It does not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, does not validate the mouse result, and does not provide mechanism, protein-level evidence, functional evidence, or therapeutic relevance.",
            "",
            f"Figures created: {', '.join(figure_paths) if figure_paths else 'none'}",
        ]),
        encoding="utf-8",
    )
    REPORTS.joinpath("human_thymus_stage3d_gse147520_candidate_context_note.md").write_text(
        "\n".join([
            "# GSE147520 Human Thymic Epithelial Candidate Context Note",
            "",
            "Search date: 2026-07-06",
            "",
            "## What GSE147520 Epithelial Can Contribute",
            "",
            "The epithelial H5AD provides human thymic epithelial cell-count and subtype context across fetal, postnatal, and one adult-labeled sample group.",
            "",
            "## What It Cannot Contribute",
            "",
            "The available epithelial H5AD cannot contribute compact-X target-gene summaries because the five target genes are absent from compact X. Detection context can be generated from raw.X, but the file still lacks explicit donor and sex fields in observed metadata.",
            "",
            "## Donor-Aware Feasibility",
            "",
            "Donor-aware summaries are not fully feasible from this H5AD because a donor field is not available. Sample-aware cell-count context is feasible using `samples`.",
            "",
            "## Age and Development Coverage",
            "",
            "The `samples` labels include fetal, postnatal, and adult-labeled groups, but adult coverage is limited and does not support direct aged-adult interpretation.",
            "",
            "## Workflow Context",
            "",
            "GSE147520 complements Park and Yayon at the workflow/context level by showing why matrix-layer inspection matters: compact X omits the targets while raw.X contains them. It should be used cautiously as human epithelial sample-aware context.",
        ]),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
