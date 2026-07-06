"""Create donor-aware input summaries for the Park human TEC H5AD.

This script creates cell-count and target-gene detection previews. It does not
run statistical tests and does not treat cells as biological replicates.
"""

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

OBS_TSV = TABLES / "human_thymus_stage3a_park_obs_fields.tsv"
GENE_TSV = TABLES / "human_thymus_stage3a_park_gene_presence.tsv"
SUMMARY_TSV = TABLES / "human_thymus_stage3a_park_pseudobulk_input_summary.tsv"
COUNTS_TSV = TABLES / "human_thymus_stage3a_park_cell_counts_by_donor_subtype.tsv"
DETECTION_TSV = TABLES / "human_thymus_stage3a_park_lox_detection_preview.tsv"
REPORT = REPORTS / "human_thymus_stage3a_park_pseudobulk_input_summary.md"

TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def choose_field(obs_rows: list[dict[str, str]], role: str, preferred: list[str]) -> tuple[str, str]:
    available = {row["field_name"] for row in obs_rows}
    candidates = [row["field_name"] for row in obs_rows if role in row["candidate_role"].split(";")]
    preferred_hits = [field for field in preferred if field in available]
    if len(preferred_hits) == 1:
        return preferred_hits[0], "selected_preferred"
    if len(preferred_hits) > 1:
        return preferred_hits[0], "selected_first_preferred_with_alternatives:" + ",".join(preferred_hits[1:])
    if len(candidates) == 1:
        return candidates[0], "selected_single_candidate"
    if len(candidates) > 1:
        return candidates[0], "ambiguous_selected_first:" + ",".join(candidates[1:])
    return "", "missing"


def vector_to_1d(matrix: object) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def main() -> None:
    if not H5AD.exists():
        raise FileNotFoundError(f"Missing Park TEC H5AD: {H5AD}")
    if not OBS_TSV.exists() or not GENE_TSV.exists():
        raise FileNotFoundError("Run schema inspection before pseudobulk input summary.")

    obs_rows = read_tsv(OBS_TSV)
    donor_field, donor_status = choose_field(obs_rows, "donor", ["donor_id"])
    sample_field, sample_status = choose_field(obs_rows, "sample", ["sample_id"])
    age_field, age_status = choose_field(obs_rows, "age_or_development", ["development_stage", "development_stage_ontology_term_id"])
    sex_field, sex_status = choose_field(obs_rows, "sex", ["sex"])
    cell_field, cell_status = choose_field(obs_rows, "cell_type", ["celltypes", "cell_type"])

    adata = ad.read_h5ad(H5AD)

    obs = adata.obs.copy()
    group_columns = [donor_field, sample_field, age_field, cell_field]
    if sex_field:
        group_columns.insert(3, sex_field)
    group_columns = [column for index, column in enumerate(group_columns) if column and column not in group_columns[:index]]

    counts = (
        obs.groupby(group_columns, observed=True)
        .size()
        .reset_index(name="cell_count")
        .sort_values(group_columns)
    )
    counts.insert(0, "dataset", "CELLxGENE_Park_TEC")
    counts_rows = counts.astype(str).to_dict(orient="records")

    gene_rows = read_tsv(GENE_TSV)
    gene_to_var_id = {
        row["target_gene"]: row["matching_var_ids"].split(";")[0]
        for row in gene_rows
        if row.get("present") == "yes" and row.get("matching_var_ids")
    }

    detection_rows: list[dict[str, str]] = []
    matrix_semantics = "partial"
    matrix_note = "X is present and raw exists, but layer/count semantics were not asserted by this Stage 3A audit; detection fractions use X > 0 only."
    grouping_index = obs[group_columns].astype(str).agg("\t".join, axis=1)
    grouped_indices = {key: np.flatnonzero(grouping_index.to_numpy() == key) for key in sorted(grouping_index.unique())}
    group_values = {row_key: row_key.split("\t") for row_key in grouped_indices}

    for gene in TARGET_GENES:
        var_id = gene_to_var_id.get(gene, "")
        if not var_id:
            for key, indices in grouped_indices.items():
                values = group_values[key]
                base = dict(zip(group_columns, values, strict=False))
                detection_rows.append({
                    "dataset": "CELLxGENE_Park_TEC",
                    **base,
                    "target_gene": gene,
                    "gene_present": "no",
                    "n_cells": str(len(indices)),
                    "n_detected": "",
                    "detection_fraction": "",
                    "mean_value": "",
                    "matrix_used": "X",
                    "matrix_semantics": matrix_semantics,
                    "notes": "Target gene was not detected in schema table.",
                })
            continue
        var_idx = int(adata.var_names.get_loc(var_id))
        vector = vector_to_1d(adata.X[:, var_idx])
        detected = vector > 0
        for key, indices in grouped_indices.items():
            values = group_values[key]
            base = dict(zip(group_columns, values, strict=False))
            n_cells = len(indices)
            n_detected = int(detected[indices].sum())
            detection_rows.append({
                "dataset": "CELLxGENE_Park_TEC",
                **base,
                "target_gene": gene,
                "gene_present": "yes",
                "n_cells": str(n_cells),
                "n_detected": str(n_detected),
                "detection_fraction": f"{(n_detected / n_cells) if n_cells else 0:.6g}",
                "mean_value": "",
                "matrix_used": "X",
                "matrix_semantics": matrix_semantics,
                "notes": "Detection preview only; no statistical test was run.",
            })

    summary_rows = [
        {"metric": "dataset", "value": "CELLxGENE_Park_TEC"},
        {"metric": "h5ad_path", "value": str(H5AD.relative_to(ROOT))},
        {"metric": "n_cells", "value": str(adata.n_obs)},
        {"metric": "n_features", "value": str(adata.n_vars)},
        {"metric": "selected_donor_field", "value": donor_field},
        {"metric": "donor_field_status", "value": donor_status},
        {"metric": "selected_sample_field", "value": sample_field},
        {"metric": "sample_field_status", "value": sample_status},
        {"metric": "selected_age_development_field", "value": age_field},
        {"metric": "age_development_field_status", "value": age_status},
        {"metric": "selected_sex_field", "value": sex_field},
        {"metric": "sex_field_status", "value": sex_status},
        {"metric": "selected_cell_type_field", "value": cell_field},
        {"metric": "cell_type_field_status", "value": cell_status},
        {"metric": "target_genes_present", "value": "; ".join(gene_to_var_id.keys())},
        {"metric": "cell_count_groups", "value": str(len(counts_rows))},
        {"metric": "detection_preview_rows", "value": str(len(detection_rows))},
        {"metric": "matrix_semantics", "value": matrix_semantics},
        {"metric": "matrix_note", "value": matrix_note},
        {"metric": "statistical_tests_run", "value": "no"},
    ]

    write_tsv(SUMMARY_TSV, summary_rows, ["metric", "value"])
    write_tsv(COUNTS_TSV, counts_rows, list(counts.columns.astype(str)))
    detection_columns = [
        "dataset",
        *group_columns,
        "target_gene",
        "gene_present",
        "n_cells",
        "n_detected",
        "detection_fraction",
        "mean_value",
        "matrix_used",
        "matrix_semantics",
        "notes",
    ]
    write_tsv(DETECTION_TSV, detection_rows, detection_columns)

    REPORTS.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "\n".join([
            "# Park Human TEC Pseudobulk Input Summary",
            "",
            "Search date: 2026-07-06",
            "",
            "The Park TEC H5AD was summarized into donor-aware input tables for a later expression-summary stage.",
            "",
            f"- Selected donor field: {donor_field} ({donor_status})",
            f"- Selected sample field: {sample_field} ({sample_status})",
            f"- Selected age/development field: {age_field} ({age_status})",
            f"- Selected sex field: {sex_field} ({sex_status})",
            f"- Selected cell-type field: {cell_field} ({cell_status})",
            f"- Cell-count groups: {len(counts_rows)}",
            f"- Target genes present: {', '.join(gene_to_var_id.keys())}",
            "- Expression preview generated: yes, detection fraction only.",
            f"- Matrix/layer semantics: {matrix_semantics}; {matrix_note}",
            "",
            "No statistical tests were run, and cells were not treated as biological replicates.",
        ]),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
