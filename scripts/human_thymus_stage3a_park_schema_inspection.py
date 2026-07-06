"""Inspect the local CELLxGENE Park TEC H5AD schema.

The script writes small schema tables and a report only. It does not write
large matrix-derived artifacts.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import anndata as ad


ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data" / "external" / "human_thymus" / "cellxgene_park_tec" / "park_tec.h5ad"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"

TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]

SCHEMA_TSV = TABLES / "human_thymus_stage3a_park_schema.tsv"
OBS_TSV = TABLES / "human_thymus_stage3a_park_obs_fields.tsv"
VAR_TSV = TABLES / "human_thymus_stage3a_park_var_fields.tsv"
LAYERS_TSV = TABLES / "human_thymus_stage3a_park_layers.tsv"
GENE_TSV = TABLES / "human_thymus_stage3a_park_gene_presence.tsv"
REPORT = REPORTS / "human_thymus_stage3a_park_schema_report.md"


def matrix_type(value: Any) -> str:
    dtype = getattr(value, "dtype", None)
    shape = getattr(value, "shape", None)
    return f"{type(value).__module__}.{type(value).__name__}; shape={shape}; dtype={dtype}"


def unique_examples(series: Any, limit: int = 12) -> tuple[int, str]:
    values = series.dropna().astype(str)
    uniques = list(values.unique())
    return len(uniques), "; ".join(uniques[:limit])


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def candidate_role(column: str) -> str:
    lower = column.lower()
    roles: list[str] = []
    if lower in {"donor_id", "donor", "participant_id"} or "donor" in lower:
        roles.append("donor")
    if lower in {"sample_id", "library", "file"} or "sample" in lower:
        roles.append("sample")
    if lower in {"development_stage", "age", "age_group"} or "development" in lower or "age" in lower:
        roles.append("age_or_development")
    if lower == "sex" or "sex" in lower:
        roles.append("sex")
    if lower in {"celltypes", "cell_type", "cell_type_ontology_term_id"} or "cell" in lower:
        roles.append("cell_type")
    if lower in {"assay", "assay_ontology_term_id", "method"} or "assay" in lower or "method" in lower:
        roles.append("assay")
    if lower in {"tissue", "tissue_ontology_term_id", "tissue_type"} or "tissue" in lower:
        roles.append("tissue")
    if lower in {"disease", "disease_ontology_term_id"} or "disease" in lower:
        roles.append("disease")
    if lower in {"suspension_type"}:
        roles.append("suspension_type")
    if lower in {"file", "method", "sort"}:
        roles.append("batch_like")
    return ";".join(roles) if roles else "other"


def main() -> None:
    if not H5AD.exists():
        raise FileNotFoundError(f"Missing Park TEC H5AD: {H5AD}")

    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        schema_rows = [
            {"metric": "h5ad_path", "value": str(H5AD.relative_to(ROOT))},
            {"metric": "file_size_bytes", "value": str(H5AD.stat().st_size)},
            {"metric": "n_obs", "value": str(adata.n_obs)},
            {"metric": "n_vars", "value": str(adata.n_vars)},
            {"metric": "raw_exists", "value": "yes" if adata.raw is not None else "no"},
            {"metric": "x_matrix_type", "value": matrix_type(adata.X)},
            {"metric": "raw_x_matrix_type", "value": matrix_type(adata.raw.X) if adata.raw is not None else "not_available"},
            {"metric": "n_obs_columns", "value": str(len(adata.obs.columns))},
            {"metric": "n_var_columns", "value": str(len(adata.var.columns))},
            {"metric": "n_layers", "value": str(len(adata.layers.keys()))},
            {"metric": "uns_keys", "value": "; ".join(map(str, adata.uns.keys()))},
        ]

        obs_rows: list[dict[str, str]] = []
        for column in adata.obs.columns:
            n_unique, examples = unique_examples(adata.obs[column])
            obs_rows.append({
                "field_name": column,
                "dtype": str(adata.obs[column].dtype),
                "non_null_count": str(int(adata.obs[column].notna().sum())),
                "unique_count": str(n_unique),
                "examples": examples,
                "candidate_role": candidate_role(column),
            })

        var_rows: list[dict[str, str]] = []
        for column in adata.var.columns:
            n_unique, examples = unique_examples(adata.var[column])
            var_rows.append({
                "field_name": column,
                "dtype": str(adata.var[column].dtype),
                "non_null_count": str(int(adata.var[column].notna().sum())),
                "unique_count": str(n_unique),
                "examples": examples,
                "candidate_role": "gene_symbol" if column == "feature_name" else "feature_metadata",
            })

        layer_rows: list[dict[str, str]] = [{
            "matrix_name": "X",
            "available": "yes",
            "matrix_type": matrix_type(adata.X),
            "notes": "Primary matrix present; semantic meaning must be interpreted from dataset documentation before mean-value summaries.",
        }]
        if adata.raw is not None:
            layer_rows.append({
                "matrix_name": "raw.X",
                "available": "yes",
                "matrix_type": matrix_type(adata.raw.X),
                "notes": "Raw object exists; exact count/normalization semantics are not asserted by this audit.",
            })
        else:
            layer_rows.append({
                "matrix_name": "raw.X",
                "available": "no",
                "matrix_type": "not_available",
                "notes": "No raw object present.",
            })
        for layer in adata.layers.keys():
            layer_rows.append({
                "matrix_name": f"layers/{layer}",
                "available": "yes",
                "matrix_type": matrix_type(adata.layers[layer]),
                "notes": "Layer present; semantics require local documentation check.",
            })

        gene_rows: list[dict[str, str]] = []
        var_names = adata.var_names.astype(str)
        symbol_columns = [column for column in adata.var.columns if column.lower() in {"feature_name", "gene_symbol", "symbol"}]
        for gene in TARGET_GENES:
            present_in_var_names = gene in set(var_names)
            matching_columns: list[str] = []
            matching_var_ids: list[str] = []
            for column in symbol_columns:
                matches = adata.var.index[adata.var[column].astype(str) == gene].astype(str).tolist()
                if matches:
                    matching_columns.append(column)
                    matching_var_ids.extend(matches)
            gene_rows.append({
                "target_gene": gene,
                "present": "yes" if present_in_var_names or matching_var_ids else "no",
                "match_location": "var_names" if present_in_var_names else ";".join(matching_columns),
                "matching_var_ids": ";".join(matching_var_ids[:20]),
                "notes": "Detected by exact feature-name match." if matching_var_ids else "No exact target-gene match detected.",
            })

        write_tsv(SCHEMA_TSV, schema_rows, ["metric", "value"])
        write_tsv(OBS_TSV, obs_rows, ["field_name", "dtype", "non_null_count", "unique_count", "examples", "candidate_role"])
        write_tsv(VAR_TSV, var_rows, ["field_name", "dtype", "non_null_count", "unique_count", "examples", "candidate_role"])
        write_tsv(LAYERS_TSV, layer_rows, ["matrix_name", "available", "matrix_type", "notes"])
        write_tsv(GENE_TSV, gene_rows, ["target_gene", "present", "match_location", "matching_var_ids", "notes"])

        donor_candidates = [row["field_name"] for row in obs_rows if "donor" in row["candidate_role"]]
        age_candidates = [row["field_name"] for row in obs_rows if "age_or_development" in row["candidate_role"]]
        sex_candidates = [row["field_name"] for row in obs_rows if "sex" in row["candidate_role"]]
        cell_candidates = [row["field_name"] for row in obs_rows if "cell_type" in row["candidate_role"]]

        REPORTS.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(
            "\n".join([
                "# Park Human TEC H5AD Schema Report",
                "",
                "Search date: 2026-07-06",
                "",
                "The local CELLxGENE Park TEC H5AD was parsed for schema and gene-presence readiness only.",
                "",
                f"- AnnData shape: {adata.n_obs} cells x {adata.n_vars} features",
                f"- Raw object exists: {'yes' if adata.raw is not None else 'no'}",
                f"- Layers found: {', '.join(map(str, adata.layers.keys())) if adata.layers.keys() else 'none'}",
                f"- Donor candidate fields: {', '.join(donor_candidates) if donor_candidates else 'none'}",
                f"- Age/development candidate fields: {', '.join(age_candidates) if age_candidates else 'none'}",
                f"- Sex candidate fields: {', '.join(sex_candidates) if sex_candidates else 'none'}",
                f"- Cell-type candidate fields: {', '.join(cell_candidates) if cell_candidates else 'none'}",
                f"- LOX-family genes detected: {', '.join(row['target_gene'] for row in gene_rows if row['present'] == 'yes')}",
                "",
                "Matrix/layer semantics remain partial because no explicit layer annotation was available in the object.",
                "",
                "This report does not make a human expression conclusion.",
            ]),
            encoding="utf-8",
        )
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
