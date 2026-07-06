"""Create donor-aware Yayon human TEC LOX-family expression-context summaries."""

from __future__ import annotations

import csv
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data" / "external" / "human_thymus" / "cellxgene_yayon_tec" / "yayon_tec.h5ad"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "results" / "figures" / "human_thymus_stage3c_yayon"
TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]
BASE_FIELDS = ["donor_id", "sample", "development_stage", "sex"]
FINE_FIELD = "cell_type_level_4_explore"
BROAD_FIELD = "cell_type_level_2"


def read_kv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


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
    if "feature_name" in var.columns:
        for gene in TARGET_GENES:
            matches = var.index[var["feature_name"].astype(str) == gene].astype(str).tolist()
            if matches:
                mapping[gene] = matches[0]
    for gene in TARGET_GENES:
        if gene in var.index.astype(str):
            mapping.setdefault(gene, gene)
    return mapping


def cell_counts(obs: pd.DataFrame, subtype_field: str) -> pd.DataFrame:
    fields = [*BASE_FIELDS, subtype_field]
    return obs.groupby(fields, observed=True).size().reset_index(name="n_cells").sort_values(fields)


def detection_rows(adata: ad.AnnData, subtype_field: str, gene_to_var: dict[str, str], matrix_used: str) -> list[dict[str, str]]:
    obs = adata.obs[[*BASE_FIELDS, subtype_field]].astype(str)
    group_fields = [*BASE_FIELDS, subtype_field]
    key = obs[group_fields].agg("\t".join, axis=1).to_numpy()
    groups = {value: np.flatnonzero(key == value) for value in sorted(pd.unique(key))}
    rows: list[dict[str, str]] = []
    for gene in TARGET_GENES:
        var_id = gene_to_var.get(gene, "")
        vals = vector(adata.X, int(adata.var_names.get_loc(var_id))) if var_id else np.array([])
        detected = vals > 0 if vals.size else np.zeros(adata.n_obs, dtype=bool)
        for group_key, indices in groups.items():
            base = dict(zip(group_fields, group_key.split("\t"), strict=False))
            n_cells = len(indices)
            n_detected = int(detected[indices].sum())
            rows.append({
                "gene": gene,
                **base,
                "n_cells": str(n_cells),
                "n_detected": str(n_detected),
                "detection_fraction": f"{(n_detected / n_cells) if n_cells else 0:.6g}",
                "matrix_used": matrix_used,
            })
    return rows


def expression_placeholder(subtype_field: str, reason: str) -> list[dict[str, str]]:
    return [{
        "gene": gene,
        "donor_id": "",
        "sample": "",
        "development_stage": "",
        "sex": "",
        subtype_field: "",
        "n_cells": "",
        "mean_value": "",
        "median_value": "",
        "nonzero_mean_value": "",
        "matrix_used": "not_generated",
        "value_semantics_label": "not_generated",
        "not_generated_reason": reason,
    } for gene in TARGET_GENES]


def dev_summary(rows: list[dict[str, str]], subtype_field: str) -> list[dict[str, str]]:
    df = pd.DataFrame(rows)
    df["detection_fraction"] = pd.to_numeric(df["detection_fraction"], errors="coerce")
    out: list[dict[str, str]] = []
    for (gene, stage, subtype), group in df.groupby(["gene", "development_stage", subtype_field], observed=True):
        out.append({
            "gene": str(gene),
            "development_stage": str(stage),
            subtype_field: str(subtype),
            "donor_count": str(group["donor_id"].nunique()),
            "sample_count": str(group["sample"].nunique()),
            "donor_sample_group_count": str(len(group)),
            "total_cells": str(int(pd.to_numeric(group["n_cells"], errors="coerce").sum())),
            "median_detection_fraction": f"{float(group['detection_fraction'].median()):.6g}",
            "mean_detection_fraction": f"{float(group['detection_fraction'].mean()):.6g}",
            "mean_value_summary": "not_generated",
            "value_semantics_label": "not_generated",
        })
    return out


def figures(dev_rows: list[dict[str, str]], det_rows: list[dict[str, str]], counts: pd.DataFrame) -> list[str]:
    made: list[str] = []
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return made
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")
    dev = pd.DataFrame(dev_rows)
    dev["mean_detection_fraction"] = pd.to_numeric(dev["mean_detection_fraction"], errors="coerce")
    plot = dev.groupby(["development_stage", "gene"], as_index=False, observed=True)["mean_detection_fraction"].mean()
    plt.figure(figsize=(11, 5))
    sns.barplot(data=plot, x="development_stage", y="mean_detection_fraction", hue="gene")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mean donor/sample detection fraction")
    plt.xlabel("Development stage")
    plt.title("Yayon TEC LOX-family detection by development stage")
    plt.tight_layout()
    path = FIGURES / "yayon_lox_detection_by_development_fine.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))

    det = pd.DataFrame(det_rows)
    det["detection_fraction"] = pd.to_numeric(det["detection_fraction"], errors="coerce")
    plot = det.groupby([FINE_FIELD, "gene"], as_index=False, observed=True)["detection_fraction"].mean()
    plt.figure(figsize=(10, 5))
    sns.barplot(data=plot, x=FINE_FIELD, y="detection_fraction", hue="gene")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mean donor/sample detection fraction")
    plt.xlabel("Fine TEC label")
    plt.title("Yayon TEC LOX-family detection by fine TEC label")
    plt.tight_layout()
    path = FIGURES / "yayon_lox_detection_by_celltype_fine.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))

    count_plot = counts.groupby("development_stage", as_index=False, observed=True)["n_cells"].sum()
    plt.figure(figsize=(9, 4))
    sns.barplot(data=count_plot, x="development_stage", y="n_cells", color="#4c78a8")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Cells")
    plt.xlabel("Development stage")
    plt.title("Yayon TEC cell counts by development stage")
    plt.tight_layout()
    path = FIGURES / "yayon_cell_counts_by_development.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))
    return made


def main() -> None:
    if not H5AD.exists():
        raise FileNotFoundError(f"Missing Yayon TEC H5AD: {H5AD}")
    semantics = read_kv(TABLES / "human_thymus_stage3c_yayon_matrix_semantics.tsv")
    detection_matrix = semantics.get("recommended_detection_matrix", "X")
    mean_allowed = semantics.get("mean_value_allowed", "no") == "yes"
    mean_matrix = semantics.get("recommended_mean_value_matrix", "not_selected")
    reason = semantics.get("reason", "Mean-value summaries were not allowed by matrix semantics audit.")
    adata = ad.read_h5ad(H5AD)
    gene_to_var = gene_map(adata.var)
    obs = adata.obs[[*BASE_FIELDS, FINE_FIELD, BROAD_FIELD]].copy()
    counts_fine = cell_counts(obs, FINE_FIELD)
    counts_broad = cell_counts(obs, BROAD_FIELD)
    det_fine = detection_rows(adata, FINE_FIELD, gene_to_var, detection_matrix)
    det_broad = detection_rows(adata, BROAD_FIELD, gene_to_var, detection_matrix)
    expr_fine = [] if mean_allowed else expression_placeholder(FINE_FIELD, reason)
    expr_broad = [] if mean_allowed else expression_placeholder(BROAD_FIELD, reason)
    dev_fine = dev_summary(det_fine, FINE_FIELD)
    dev_broad = dev_summary(det_broad, BROAD_FIELD)

    write_tsv(TABLES / "human_thymus_stage3c_yayon_cell_counts_fine.tsv", counts_fine.astype(str).to_dict(orient="records"), [*BASE_FIELDS, FINE_FIELD, "n_cells"])
    write_tsv(TABLES / "human_thymus_stage3c_yayon_cell_counts_broad.tsv", counts_broad.astype(str).to_dict(orient="records"), [*BASE_FIELDS, BROAD_FIELD, "n_cells"])
    write_tsv(TABLES / "human_thymus_stage3c_yayon_lox_detection_by_donor_fine.tsv", det_fine, ["gene", *BASE_FIELDS, FINE_FIELD, "n_cells", "n_detected", "detection_fraction", "matrix_used"])
    write_tsv(TABLES / "human_thymus_stage3c_yayon_lox_detection_by_donor_broad.tsv", det_broad, ["gene", *BASE_FIELDS, BROAD_FIELD, "n_cells", "n_detected", "detection_fraction", "matrix_used"])
    expr_cols_f = ["gene", *BASE_FIELDS, FINE_FIELD, "n_cells", "mean_value", "median_value", "nonzero_mean_value", "matrix_used", "value_semantics_label", "not_generated_reason"]
    expr_cols_b = ["gene", *BASE_FIELDS, BROAD_FIELD, "n_cells", "mean_value", "median_value", "nonzero_mean_value", "matrix_used", "value_semantics_label", "not_generated_reason"]
    write_tsv(TABLES / "human_thymus_stage3c_yayon_lox_expression_by_donor_fine.tsv", expr_fine, expr_cols_f)
    write_tsv(TABLES / "human_thymus_stage3c_yayon_lox_expression_by_donor_broad.tsv", expr_broad, expr_cols_b)
    write_tsv(TABLES / "human_thymus_stage3c_yayon_lox_summary_by_development_fine.tsv", dev_fine, ["gene", "development_stage", FINE_FIELD, "donor_count", "sample_count", "donor_sample_group_count", "total_cells", "median_detection_fraction", "mean_detection_fraction", "mean_value_summary", "value_semantics_label"])
    write_tsv(TABLES / "human_thymus_stage3c_yayon_lox_summary_by_development_broad.tsv", dev_broad, ["gene", "development_stage", BROAD_FIELD, "donor_count", "sample_count", "donor_sample_group_count", "total_cells", "median_detection_fraction", "mean_detection_fraction", "mean_value_summary", "value_semantics_label"])

    figure_paths = figures(dev_fine, det_fine, counts_fine)
    donors = int(obs["donor_id"].nunique())
    samples = int(obs["sample"].nunique())
    stages = int(obs["development_stage"].nunique())
    fine = int(obs[FINE_FIELD].nunique())
    broad = int(obs[BROAD_FIELD].nunique())
    fetal = [x for x in sorted(obs["development_stage"].astype(str).unique()) if "week post-fertilization" in x]
    postnatal = [x for x in sorted(obs["development_stage"].astype(str).unique()) if "month" in x or "year" in x or "newborn" in x]
    genes_present = [g for g in TARGET_GENES if g in gene_to_var]

    REPORTS.mkdir(parents=True, exist_ok=True)
    REPORTS.joinpath("human_thymus_stage3c_yayon_expression_summary_report.md").write_text(
        "\n".join([
            "# Yayon Human TEC LOX-Family Expression Summary",
            "",
            "Search date: 2026-07-06",
            "",
            "## Dataset and Local File Status",
            "",
            f"- Local H5AD path: `{H5AD.relative_to(ROOT)}`",
            f"- AnnData shape: {adata.n_obs} cells x {adata.n_vars} features",
            "- Local H5AD is under ignored external data path.",
            "",
            "## Matrix Semantics Decision",
            "",
            f"- matrix_semantics_clear: {semantics.get('matrix_semantics_clear', 'partial')}",
            f"- detection matrix used: {detection_matrix}",
            f"- mean-value matrix used: {mean_matrix}",
            f"- mean-value summaries generated: {'yes' if mean_allowed else 'no'}",
            f"- reason: {reason}",
            "",
            "## Fields Used",
            "",
            "- donor/sample fields: `donor_id`, `sample`",
            "- age/development field: `development_stage`",
            "- sex field: `sex`",
            f"- fine annotation field: `{FINE_FIELD}`",
            f"- broad annotation field: `{BROAD_FIELD}`",
            "",
            "## Dataset Breadth",
            "",
            f"- donors: {donors}",
            f"- samples: {samples}",
            f"- development stages: {stages}",
            f"- fine TEC labels: {fine}",
            f"- broad TEC labels: {broad}",
            f"- fine donor/sample/subtype groups: {len(counts_fine)}",
            f"- broad donor/sample/subtype groups: {len(counts_broad)}",
            "",
            "## Gene Coverage",
            "",
            f"- LOX-family genes present: {', '.join(genes_present)}",
            "",
            "## Cell-Count Balance",
            "",
            "Cell counts are uneven across donor/sample/development/subtype groups. Stage 3C outputs therefore summarize donor/sample groups and do not use cells as independent biological replicates.",
            "",
            "## Detection-Summary Overview",
            "",
            "Detection summaries were generated for fine and broad TEC annotation levels using target-gene values greater than zero in `X`.",
            "",
            "## Mean-Value-Summary Status",
            "",
            "Mean-value summaries were not generated because local matrix semantics remain partial.",
            "",
            "## Park Workflow Context",
            "",
            "The same cautious workflow used for Park was applied here: schema inspection, matrix-semantics audit, donor/sample-aware detection summaries, no statistical tests, and no cross-dataset biological conclusion.",
            "",
            "## Stage 3D Recommendation",
            "",
            "Stage 3D should proceed to the GSE147520 epithelial H5AD if local access and schema checks remain feasible.",
            "",
            "## Limitations",
            "",
            "- The Yayon TEC dataset is useful as human TEC context, with fetal and early postnatal labels rather than an aged-adult-focused design.",
            "- Matrix semantics remain partial for mean-value summaries.",
            "- No statistical tests were run.",
            "- No mouse-human comparison was performed.",
            "",
            "This Stage 3C summary provides exploratory donor-aware transcript-level context for LOX-family genes in the Yayon human TEC dataset. It does not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, does not validate the mouse result, and does not provide mechanism, protein-level evidence, functional evidence, or therapeutic relevance.",
            "",
            f"Figures created: {', '.join(figure_paths) if figure_paths else 'none'}",
        ]),
        encoding="utf-8",
    )
    REPORTS.joinpath("human_thymus_stage3c_yayon_candidate_context_note.md").write_text(
        "\n".join([
            "# Yayon Human TEC Candidate Context Note",
            "",
            "Search date: 2026-07-06",
            "",
            "## What Yayon TEC Can Contribute",
            "",
            "The Yayon TEC subset provides a locally parseable human TEC resource with donor, sample, development-stage, sex, fine TEC label, and broader TEC label fields. It supports donor/sample-aware detection summaries for LOX-family transcripts.",
            "",
            "## What Yayon TEC Cannot Contribute",
            "",
            "This dataset cannot by itself establish aged-adult thymus behavior, protein abundance, activity, functional outcome, or cross-species interpretation. It should not be used as a direct aging validation dataset.",
            "",
            "## Donor-Aware Feasibility",
            "",
            f"Donor-aware summaries are feasible across {donors} donors, {samples} samples, {stages} development stages, {fine} fine TEC labels, and {broad} broader TEC labels.",
            "",
            "## Age and Development Coverage",
            "",
            f"Development-stage labels include {len(fetal)} fetal-stage labels and {len(postnatal)} postnatal or early-life labels. This supports developmental and TEC context more than aging-oriented interpretation.",
            "",
            "## Workflow Context",
            "",
            "Yayon complements Park at the workflow and context level because the same donor-aware parsing approach can be applied to another human TEC resource. Any broader cross-dataset synthesis should remain exploratory unless later stages add suitable aged-adult datasets.",
        ]),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
