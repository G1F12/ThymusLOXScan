"""Inspect the local GSE147520 epithelial H5AD schema."""

from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path
from typing import Any

import anndata as ad


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external" / "human_thymus" / "GSE147520"
H5AD = DATA_DIR / "GSE147520_epithelial_cells.h5ad"
H5AD_GZ = DATA_DIR / "GSE147520_epithelial_cells.h5ad.gz"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]


def ensure_h5ad() -> Path:
    if H5AD.exists():
        return H5AD
    if H5AD_GZ.exists():
        with gzip.open(H5AD_GZ, "rb") as src, H5AD.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        return H5AD
    raise FileNotFoundError(f"Missing {H5AD} and {H5AD_GZ}")


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
    if "donor" in lower:
        roles.append("donor")
    if lower in {"samples", "sample", "sample_id"} or "sample" in lower:
        roles.append("sample")
    if "age" in lower or "development" in lower or lower == "samples":
        roles.append("age_or_development")
    if lower == "sex" or "sex" in lower:
        roles.append("sex")
    if "cell_type" in lower or "cell_types" in lower or "leiden" in lower:
        roles.append("cell_type")
    return ";".join(roles) if roles else "other"


def find_gene_matches(var: Any, gene: str) -> list[str]:
    matches: list[str] = []
    if gene in var.index.astype(str):
        matches.append(gene)
    for column in var.columns:
        column_matches = var.index[var[column].astype(str) == gene].astype(str).tolist()
        matches.extend(column_matches)
    return list(dict.fromkeys(matches))


def main() -> None:
    path = ensure_h5ad()
    adata = ad.read_h5ad(path, backed="r")
    try:
        schema_rows = [
            {"metric": "h5ad_path", "value": str(path.relative_to(ROOT))},
            {"metric": "h5ad_gz_path", "value": str(H5AD_GZ.relative_to(ROOT)) if H5AD_GZ.exists() else "not_available"},
            {"metric": "file_size_bytes", "value": str(path.stat().st_size)},
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
                "candidate_role": "feature_metadata",
            })
        layer_rows = [{
            "matrix_name": "X",
            "available": "yes",
            "matrix_type": matrix_type(adata.X),
            "notes": "Primary matrix present; target-gene availability and value semantics require audit.",
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
        gene_rows = []
        raw_var = adata.raw.var if adata.raw is not None else None
        for gene in TARGET_GENES:
            var_matches = find_gene_matches(adata.var, gene)
            raw_matches = find_gene_matches(raw_var, gene) if raw_var is not None else []
            gene_rows.append({
                "target_gene": gene,
                "present": "yes" if var_matches or raw_matches else "no",
                "match_location": ";".join(x for x, matches in [("var", var_matches), ("raw.var", raw_matches)] if matches),
                "matching_var_ids": ";".join(var_matches),
                "matching_raw_var_ids": ";".join(raw_matches),
                "notes": "Target gene detected." if var_matches or raw_matches else "Target gene not present in this epithelial H5AD feature set.",
            })

        write_tsv(TABLES / "human_thymus_stage3d_gse147520_schema.tsv", schema_rows, ["metric", "value"])
        write_tsv(TABLES / "human_thymus_stage3d_gse147520_obs_fields.tsv", obs_rows, ["field_name", "dtype", "non_null_count", "unique_count", "examples", "candidate_role"])
        write_tsv(TABLES / "human_thymus_stage3d_gse147520_var_fields.tsv", var_rows, ["field_name", "dtype", "non_null_count", "unique_count", "examples", "candidate_role"])
        write_tsv(TABLES / "human_thymus_stage3d_gse147520_layers.tsv", layer_rows, ["matrix_name", "available", "matrix_type", "notes"])
        write_tsv(TABLES / "human_thymus_stage3d_gse147520_gene_presence.tsv", gene_rows, ["target_gene", "present", "match_location", "matching_var_ids", "matching_raw_var_ids", "notes"])

        REPORTS.mkdir(parents=True, exist_ok=True)
        REPORTS.joinpath("human_thymus_stage3d_gse147520_schema_report.md").write_text(
            "\n".join([
                "# GSE147520 Human Thymic Epithelial H5AD Schema Report",
                "",
                "Search date: 2026-07-06",
                "",
                "The local GSE147520 epithelial H5AD was parsed for schema and gene-presence readiness only.",
                "",
                f"- AnnData shape: {adata.n_obs} cells x {adata.n_vars} features",
                f"- Raw object exists: {'yes' if adata.raw is not None else 'no'}",
                f"- Layers found: {', '.join(map(str, adata.layers.keys())) if adata.layers.keys() else 'none'}",
                "- Donor field: not available in observed metadata",
                "- Sample field: samples",
                "- Age/development field: samples",
                "- Sex field: not available in observed metadata",
                "- Fine cell-type field: cell_types_epith",
                "- Broad cell-type field: derived from cell_types_epith for summary grouping",
                f"- LOX-family genes detected: {', '.join(row['target_gene'] for row in gene_rows if row['present'] == 'yes') or 'none'}",
                "",
                "The target genes are absent from compact X but present in raw.var/raw.X, so downstream Stage 3D detection summaries should use raw.X while avoiding mean-value summaries unless value semantics are clear.",
                "",
                "This report does not make a human expression conclusion.",
            ]),
            encoding="utf-8",
        )
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
