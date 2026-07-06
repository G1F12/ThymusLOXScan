"""Create donor-aware Park human TEC LOX-family expression-context summaries."""

from __future__ import annotations

import csv
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data" / "external" / "human_thymus" / "cellxgene_park_tec" / "park_tec.h5ad"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "results" / "figures" / "human_thymus_stage3b_park"

SEMANTICS_TSV = TABLES / "human_thymus_stage3b_park_matrix_semantics.tsv"

CELL_COUNTS_FINE = TABLES / "human_thymus_stage3b_park_cell_counts_fine.tsv"
CELL_COUNTS_BROAD = TABLES / "human_thymus_stage3b_park_cell_counts_broad.tsv"
DETECTION_FINE = TABLES / "human_thymus_stage3b_park_lox_detection_by_donor_fine.tsv"
DETECTION_BROAD = TABLES / "human_thymus_stage3b_park_lox_detection_by_donor_broad.tsv"
EXPRESSION_FINE = TABLES / "human_thymus_stage3b_park_lox_expression_by_donor_fine.tsv"
EXPRESSION_BROAD = TABLES / "human_thymus_stage3b_park_lox_expression_by_donor_broad.tsv"
DEV_FINE = TABLES / "human_thymus_stage3b_park_lox_summary_by_development_fine.tsv"
DEV_BROAD = TABLES / "human_thymus_stage3b_park_lox_summary_by_development_broad.tsv"

REPORT = REPORTS / "human_thymus_stage3b_park_expression_summary_report.md"
CONTEXT_NOTE = REPORTS / "human_thymus_stage3b_park_candidate_context_note.md"

TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]
BASE_FIELDS = ["donor_id", "sample_id", "development_stage", "sex"]


def read_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            values[row["metric"]] = row["value"]
    return values


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def as_vector(matrix: object, var_idx: int) -> np.ndarray:
    values = matrix[:, var_idx]
    if sparse.issparse(values):
        return np.asarray(values.toarray()).ravel()
    if hasattr(values, "to_memory"):
        values = values.to_memory()
    if sparse.issparse(values):
        return np.asarray(values.toarray()).ravel()
    return np.asarray(values).ravel()


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
    return (
        obs.groupby(fields, observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(fields)
    )


def detection_rows(adata: ad.AnnData, subtype_field: str, gene_to_var: dict[str, str]) -> list[dict[str, str]]:
    obs = adata.obs[[*BASE_FIELDS, subtype_field]].astype(str).copy()
    group_fields = [*BASE_FIELDS, subtype_field]
    group_key = obs[group_fields].agg("\t".join, axis=1).to_numpy()
    unique_keys = sorted(pd.unique(group_key))
    grouped_indices = {key: np.flatnonzero(group_key == key) for key in unique_keys}
    rows: list[dict[str, str]] = []
    for gene in TARGET_GENES:
        var_id = gene_to_var.get(gene, "")
        vector = as_vector(adata.X, int(adata.var_names.get_loc(var_id))) if var_id else np.array([])
        detected = vector > 0 if vector.size else np.zeros(adata.n_obs, dtype=bool)
        for key, indices in grouped_indices.items():
            values = key.split("\t")
            row = dict(zip(group_fields, values, strict=False))
            n_cells = len(indices)
            n_detected = int(detected[indices].sum())
            rows.append({
                "gene": gene,
                **row,
                "n_cells": str(n_cells),
                "n_detected": str(n_detected),
                "detection_fraction": f"{(n_detected / n_cells) if n_cells else 0:.6g}",
            })
    return rows


def expression_not_generated_rows(subtype_field: str, reason: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for gene in TARGET_GENES:
        rows.append({
            "gene": gene,
            "donor_id": "",
            "sample_id": "",
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
        })
    return rows


def development_summary(detection: list[dict[str, str]], subtype_field: str) -> list[dict[str, str]]:
    df = pd.DataFrame(detection)
    df["detection_fraction"] = pd.to_numeric(df["detection_fraction"], errors="coerce")
    grouped = df.groupby(["gene", "development_stage", subtype_field], observed=True)
    rows: list[dict[str, str]] = []
    for (gene, stage, subtype), group in grouped:
        rows.append({
            "gene": str(gene),
            "development_stage": str(stage),
            subtype_field: str(subtype),
            "donor_count": str(group["donor_id"].nunique()),
            "sample_count": str(group["sample_id"].nunique()),
            "donor_sample_group_count": str(len(group)),
            "total_cells": str(int(pd.to_numeric(group["n_cells"], errors="coerce").sum())),
            "median_detection_fraction": f"{float(group['detection_fraction'].median()):.6g}",
            "mean_detection_fraction": f"{float(group['detection_fraction'].mean()):.6g}",
            "mean_value_summary": "not_generated",
            "value_semantics_label": "not_generated",
        })
    return rows


def make_figures(dev_fine_rows: list[dict[str, str]], det_fine_rows: list[dict[str, str]], counts_fine: pd.DataFrame) -> list[str]:
    created: list[str] = []
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return created
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    dev_df = pd.DataFrame(dev_fine_rows)
    dev_df["mean_detection_fraction"] = pd.to_numeric(dev_df["mean_detection_fraction"], errors="coerce")
    for_plot = dev_df.groupby(["development_stage", "gene"], as_index=False)["mean_detection_fraction"].mean()
    plt.figure(figsize=(11, 5))
    sns.barplot(data=for_plot, x="development_stage", y="mean_detection_fraction", hue="gene")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mean donor/sample detection fraction")
    plt.xlabel("Development stage")
    plt.title("Park TEC LOX-family detection by development stage")
    plt.tight_layout()
    path = FIGURES / "park_lox_detection_by_development_fine.png"
    plt.savefig(path, dpi=200)
    plt.close()
    created.append(str(path.relative_to(ROOT)))

    det_df = pd.DataFrame(det_fine_rows)
    det_df["detection_fraction"] = pd.to_numeric(det_df["detection_fraction"], errors="coerce")
    celltype_plot = det_df.groupby(["celltypes", "gene"], as_index=False)["detection_fraction"].mean()
    plt.figure(figsize=(10, 5))
    sns.barplot(data=celltype_plot, x="celltypes", y="detection_fraction", hue="gene")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mean donor/sample detection fraction")
    plt.xlabel("Fine TEC label")
    plt.title("Park TEC LOX-family detection by fine TEC label")
    plt.tight_layout()
    path = FIGURES / "park_lox_detection_by_celltype_fine.png"
    plt.savefig(path, dpi=200)
    plt.close()
    created.append(str(path.relative_to(ROOT)))

    count_plot = counts_fine.groupby("development_stage", as_index=False, observed=True)["n_cells"].sum()
    plt.figure(figsize=(9, 4))
    sns.barplot(data=count_plot, x="development_stage", y="n_cells", color="#4c78a8")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Cells")
    plt.xlabel("Development stage")
    plt.title("Park TEC cell counts by development stage")
    plt.tight_layout()
    path = FIGURES / "park_cell_counts_by_development.png"
    plt.savefig(path, dpi=200)
    plt.close()
    created.append(str(path.relative_to(ROOT)))
    return created


def main() -> None:
    if not H5AD.exists():
        raise FileNotFoundError(f"Missing Park TEC H5AD: {H5AD}")
    if not SEMANTICS_TSV.exists():
        raise FileNotFoundError("Run matrix semantics audit before expression summary.")

    semantics = read_key_values(SEMANTICS_TSV)
    detection_matrix = semantics.get("recommended_detection_matrix", "X")
    mean_allowed = semantics.get("mean_value_allowed", "no") == "yes"
    mean_matrix = semantics.get("recommended_mean_value_matrix", "not_selected")
    reason = semantics.get("reason", "Mean-value summaries were not allowed by matrix semantics audit.")

    adata = ad.read_h5ad(H5AD)
    gene_to_var = gene_map(adata.var)
    obs = adata.obs[[*BASE_FIELDS, "celltypes", "cell_type"]].copy()
    counts_fine = cell_counts(obs, "celltypes")
    counts_broad = cell_counts(obs, "cell_type")
    det_fine = detection_rows(adata, "celltypes", gene_to_var)
    det_broad = detection_rows(adata, "cell_type", gene_to_var)
    expr_fine = expression_not_generated_rows("celltypes", reason) if not mean_allowed else []
    expr_broad = expression_not_generated_rows("cell_type", reason) if not mean_allowed else []
    dev_fine = development_summary(det_fine, "celltypes")
    dev_broad = development_summary(det_broad, "cell_type")

    write_tsv(CELL_COUNTS_FINE, counts_fine.astype(str).to_dict(orient="records"), [*BASE_FIELDS, "celltypes", "n_cells"])
    write_tsv(CELL_COUNTS_BROAD, counts_broad.astype(str).to_dict(orient="records"), [*BASE_FIELDS, "cell_type", "n_cells"])
    write_tsv(DETECTION_FINE, det_fine, ["gene", *BASE_FIELDS, "celltypes", "n_cells", "n_detected", "detection_fraction"])
    write_tsv(DETECTION_BROAD, det_broad, ["gene", *BASE_FIELDS, "cell_type", "n_cells", "n_detected", "detection_fraction"])
    expression_fine_cols = ["gene", *BASE_FIELDS, "celltypes", "n_cells", "mean_value", "median_value", "nonzero_mean_value", "matrix_used", "value_semantics_label", "not_generated_reason"]
    expression_broad_cols = ["gene", *BASE_FIELDS, "cell_type", "n_cells", "mean_value", "median_value", "nonzero_mean_value", "matrix_used", "value_semantics_label", "not_generated_reason"]
    write_tsv(EXPRESSION_FINE, expr_fine, expression_fine_cols)
    write_tsv(EXPRESSION_BROAD, expr_broad, expression_broad_cols)
    write_tsv(DEV_FINE, dev_fine, ["gene", "development_stage", "celltypes", "donor_count", "sample_count", "donor_sample_group_count", "total_cells", "median_detection_fraction", "mean_detection_fraction", "mean_value_summary", "value_semantics_label"])
    write_tsv(DEV_BROAD, dev_broad, ["gene", "development_stage", "cell_type", "donor_count", "sample_count", "donor_sample_group_count", "total_cells", "median_detection_fraction", "mean_detection_fraction", "mean_value_summary", "value_semantics_label"])

    figure_paths = make_figures(dev_fine, det_fine, counts_fine)

    donors = int(obs["donor_id"].nunique())
    samples = int(obs["sample_id"].nunique())
    stages = int(obs["development_stage"].nunique())
    fine_labels = int(obs["celltypes"].nunique())
    broad_labels = int(obs["cell_type"].nunique())
    total_groups_fine = len(counts_fine)
    total_groups_broad = len(counts_broad)
    present_genes = [gene for gene in TARGET_GENES if gene in gene_to_var]
    fetal_like = [stage for stage in sorted(obs["development_stage"].astype(str).unique()) if "week post-fertilization" in stage]
    postnatal_like = [stage for stage in sorted(obs["development_stage"].astype(str).unique()) if "year" in stage or "month" in stage]

    REPORTS.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "\n".join([
            "# Park Human TEC LOX-Family Expression Summary",
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
            "- donor/sample fields: `donor_id`, `sample_id`",
            "- age/development field: `development_stage`",
            "- sex field: `sex`",
            "- fine annotation field: `celltypes`",
            "- broad annotation field: `cell_type`",
            "",
            "## Dataset Breadth",
            "",
            f"- donors: {donors}",
            f"- samples: {samples}",
            f"- development stages: {stages}",
            f"- fine TEC labels: {fine_labels}",
            f"- broad TEC labels: {broad_labels}",
            f"- fine donor/sample/subtype groups: {total_groups_fine}",
            f"- broad donor/sample/subtype groups: {total_groups_broad}",
            "",
            "## Gene Coverage",
            "",
            f"- LOX-family genes present: {', '.join(present_genes)}",
            "",
            "## Cell-Count Balance",
            "",
            "Cell counts are uneven across donor/sample/development/subtype groups. Stage 3B outputs therefore summarize donor/sample groups and do not use cells as independent biological replicates.",
            "",
            "## Detection-Summary Overview",
            "",
            "Detection summaries were generated for fine and broad TEC annotation levels using target-gene values greater than zero in `X`.",
            "",
            "## Mean-Value-Summary Status",
            "",
            "Mean-value summaries were not generated because local matrix semantics remain partial.",
            "",
            "## Stage 3C Recommendation",
            "",
            "Stage 3C should proceed to the Yayon TEC dataset to test whether the same donor-aware parsing workflow can be applied in an independent human TEC resource.",
            "",
            "## Limitations",
            "",
            "- The Park dataset is primarily developmental/TEC context, with fetal and selected postnatal stages rather than an aged-adult-focused design.",
            "- Matrix semantics remain partial for mean-value summaries.",
            "- No statistical tests were run.",
            "- No mouse-human comparison was performed.",
            "",
            "This Stage 3B summary provides exploratory donor-aware transcript-level context for LOX-family genes in the Park human TEC dataset. It does not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, does not validate the mouse result, and does not provide mechanism, protein-level evidence, functional evidence, or therapeutic relevance.",
            "",
            f"Figures created: {', '.join(figure_paths) if figure_paths else 'none'}",
        ]),
        encoding="utf-8",
    )

    CONTEXT_NOTE.write_text(
        "\n".join([
            "# Park Human TEC Candidate Context Note",
            "",
            "Search date: 2026-07-06",
            "",
            "## What Park TEC Can Contribute",
            "",
            "The Park TEC dataset provides a locally parseable human TEC-focused resource with donor, sample, development-stage, sex, fine TEC label, and broad TEC label fields. It supports donor/sample-aware detection summaries for LOX-family transcripts across human TEC labels.",
            "",
            "## What Park TEC Cannot Contribute",
            "",
            "This dataset cannot by itself establish aged-adult thymus behavior, protein abundance, activity, functional outcome, or cross-species interpretation. It is not a direct aged adult thymus validation dataset.",
            "",
            "## Donor-Aware Feasibility",
            "",
            f"Donor-aware summaries are feasible across {donors} donors, {samples} samples, {stages} development stages, {fine_labels} fine TEC labels, and {broad_labels} broad TEC labels.",
            "",
            "## Age and Development Coverage",
            "",
            f"Development-stage labels include {len(fetal_like)} fetal-stage labels and {len(postnatal_like)} postnatal-stage labels. Because the dataset mainly spans developmental/fetal/pediatric stages, it should be used as human developmental and TEC context rather than direct aging validation.",
            "",
            "## Recommended Use",
            "",
            "Use Park TEC as the first human matrix-parsing reference and as a donor-aware developmental TEC context dataset. Cross-dataset work should next test the same workflow on Yayon TEC before any broader synthesis.",
        ]),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
