"""Inspect the local CELLxGENE Yayon TEC H5AD schema."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import anndata as ad


ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data" / "external" / "human_thymus" / "cellxgene_yayon_tec" / "yayon_tec.h5ad"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def matrix_type(value: Any) -> str:
    return f"{type(value).__module__}.{type(value).__name__}; shape={getattr(value, 'shape', '')}; dtype={getattr(value, 'dtype', '')}"


def examples(series: Any, limit: int = 12) -> tuple[int, str]:
    vals = series.dropna().astype(str)
    uniq = list(vals.unique())
    return len(uniq), "; ".join(uniq[:limit])


def role(column: str) -> str:
    lower = column.lower()
    roles: list[str] = []
    if lower == "donor_id" or "donor" in lower:
        roles.append("donor")
    if lower in {"sample", "sample_id", "sample_uuid", "library_uuid"} or "sample" in lower or "library" in lower:
        roles.append("sample")
    if lower in {"donor_age", "development_stage", "age_group"} or "age" in lower or "development" in lower:
        roles.append("age_or_development")
    if lower == "sex" or "sex" in lower:
        roles.append("sex")
    if "cell_type" in lower or lower in {"author_cell_type", "cell_state"}:
        roles.append("cell_type")
    if lower in {"assay", "assay_ontology_term_id", "sequencing_platform"}:
        roles.append("assay")
    if "tissue" in lower:
        roles.append("tissue")
    if "disease" in lower:
        roles.append("disease")
    if lower == "suspension_type":
        roles.append("suspension_type")
    return ";".join(roles) if roles else "other"


def gene_mapping(var: Any) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {gene: [] for gene in TARGET_GENES}
    for gene in TARGET_GENES:
        if gene in var.index.astype(str):
            out[gene].append(gene)
        for column in var.columns:
            matches = var.index[var[column].astype(str) == gene].astype(str).tolist()
            out[gene].extend(matches)
    return out


def main() -> None:
    if not H5AD.exists():
        raise FileNotFoundError(f"Missing Yayon TEC H5AD: {H5AD}")
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
        obs_rows = []
        for column in adata.obs.columns:
            n_unique, ex = examples(adata.obs[column])
            obs_rows.append({
                "field_name": column,
                "dtype": str(adata.obs[column].dtype),
                "non_null_count": str(int(adata.obs[column].notna().sum())),
                "unique_count": str(n_unique),
                "examples": ex,
                "candidate_role": role(column),
            })
        var_rows = []
        for column in adata.var.columns:
            n_unique, ex = examples(adata.var[column])
            var_rows.append({
                "field_name": column,
                "dtype": str(adata.var[column].dtype),
                "non_null_count": str(int(adata.var[column].notna().sum())),
                "unique_count": str(n_unique),
                "examples": ex,
                "candidate_role": "gene_symbol" if column == "feature_name" else "feature_metadata",
            })
        layer_rows = [{
            "matrix_name": "X",
            "available": "yes",
            "matrix_type": matrix_type(adata.X),
            "notes": "Primary matrix present; value semantics require audit before mean summaries.",
        }, {
            "matrix_name": "raw.X",
            "available": "yes" if adata.raw is not None else "no",
            "matrix_type": matrix_type(adata.raw.X) if adata.raw is not None else "not_available",
            "notes": "Raw object status recorded from AnnData.",
        }]
        for layer in adata.layers.keys():
            layer_rows.append({
                "matrix_name": f"layers/{layer}",
                "available": "yes",
                "matrix_type": matrix_type(adata.layers[layer]),
                "notes": "Layer present; semantics require audit.",
            })
        mapping = gene_mapping(adata.var)
        gene_rows = []
        for gene, matches in mapping.items():
            gene_rows.append({
                "target_gene": gene,
                "present": "yes" if matches else "no",
                "match_location": "feature_name" if matches else "",
                "matching_var_ids": ";".join(dict.fromkeys(matches)),
                "notes": "Detected by exact feature-name match." if matches else "No exact target-gene match detected.",
            })

        write_tsv(TABLES / "human_thymus_stage3c_yayon_schema.tsv", schema_rows, ["metric", "value"])
        write_tsv(TABLES / "human_thymus_stage3c_yayon_obs_fields.tsv", obs_rows, ["field_name", "dtype", "non_null_count", "unique_count", "examples", "candidate_role"])
        write_tsv(TABLES / "human_thymus_stage3c_yayon_var_fields.tsv", var_rows, ["field_name", "dtype", "non_null_count", "unique_count", "examples", "candidate_role"])
        write_tsv(TABLES / "human_thymus_stage3c_yayon_layers.tsv", layer_rows, ["matrix_name", "available", "matrix_type", "notes"])
        write_tsv(TABLES / "human_thymus_stage3c_yayon_gene_presence.tsv", gene_rows, ["target_gene", "present", "match_location", "matching_var_ids", "notes"])

        REPORTS.mkdir(parents=True, exist_ok=True)
        REPORTS.joinpath("human_thymus_stage3c_yayon_schema_report.md").write_text(
            "\n".join([
                "# Yayon Human TEC H5AD Schema Report",
                "",
                "Search date: 2026-07-06",
                "",
                "The local CELLxGENE Yayon TEC H5AD was parsed for schema and gene-presence readiness only.",
                "",
                f"- AnnData shape: {adata.n_obs} cells x {adata.n_vars} features",
                f"- Raw object exists: {'yes' if adata.raw is not None else 'no'}",
                f"- Layers found: {', '.join(map(str, adata.layers.keys())) if adata.layers.keys() else 'none'}",
                "- Selected donor candidate: donor_id",
                "- Selected sample candidate: sample",
                "- Selected age/development candidate: development_stage",
                "- Selected sex candidate: sex",
                "- Selected fine cell-type candidate: cell_type_level_4_explore",
                "- Selected broad cell-type candidate: cell_type_level_2",
                f"- LOX-family genes detected: {', '.join(row['target_gene'] for row in gene_rows if row['present'] == 'yes')}",
                "",
                "This report does not make a human expression conclusion.",
            ]),
            encoding="utf-8",
        )
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
