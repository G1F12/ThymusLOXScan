"""Audit X/raw.X semantics for Yayon human TEC target-gene summaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse


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


def desc(matrix: Any) -> dict[str, str]:
    return {
        "shape": str(getattr(matrix, "shape", "")),
        "dtype": str(getattr(matrix, "dtype", "")),
        "storage": "sparse" if sparse.issparse(matrix) or "sparse" in type(matrix).__module__.lower() else "dense",
        "class": f"{type(matrix).__module__}.{type(matrix).__name__}",
    }


def vector(matrix: Any, var_idx: int) -> np.ndarray:
    vals = matrix[:, var_idx]
    if sparse.issparse(vals):
        return np.asarray(vals.toarray()).ravel()
    if hasattr(vals, "to_memory"):
        vals = vals.to_memory()
    if sparse.issparse(vals):
        return np.asarray(vals.toarray()).ravel()
    return np.asarray(vals).ravel()


def mostly_integer(vals: np.ndarray) -> str:
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return "unknown"
    return "yes" if float(np.mean(np.isclose(finite, np.round(finite), atol=1e-6))) >= 0.995 else "no"


def unique_capped(vals: np.ndarray, cap: int = 5000) -> str:
    if vals.size > cap:
        idx = np.linspace(0, vals.size - 1, cap, dtype=int)
        return f"{len(np.unique(vals[idx]))}_in_{cap}_sample"
    return str(len(np.unique(vals)))


def gene_map(var: Any) -> dict[str, str]:
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


def main() -> None:
    if not H5AD.exists():
        raise FileNotFoundError(f"Missing Yayon TEC H5AD: {H5AD}")
    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        raw_exists = adata.raw is not None
        x_desc = desc(adata.X)
        raw_desc = desc(adata.raw.X) if raw_exists else {"shape": "not_available", "dtype": "not_available", "storage": "not_available", "class": "not_available"}
        x_map = gene_map(adata.var)
        raw_map = gene_map(adata.raw.var) if raw_exists else {}
        rows: list[dict[str, str]] = []
        for gene in TARGET_GENES:
            x_var = x_map.get(gene, "")
            raw_var = raw_map.get(gene, "")
            x_vals = vector(adata.X, int(adata.var_names.get_loc(x_var))) if x_var else np.array([])
            raw_vals = vector(adata.raw.X, int(adata.raw.var_names.get_loc(raw_var))) if raw_exists and raw_var else np.array([])
            identical = bool(x_vals.size and raw_vals.size and np.array_equal(x_vals, raw_vals))
            rows.append({
                "gene": gene,
                "x_var_id": x_var,
                "raw_var_id": raw_var,
                "x_min": f"{float(np.min(x_vals)):.6g}" if x_vals.size else "",
                "x_max": f"{float(np.max(x_vals)):.6g}" if x_vals.size else "",
                "raw_min": f"{float(np.min(raw_vals)):.6g}" if raw_vals.size else "",
                "raw_max": f"{float(np.max(raw_vals)):.6g}" if raw_vals.size else "",
                "x_nonzero_fraction": f"{float(np.mean(x_vals > 0)):.6g}" if x_vals.size else "",
                "raw_nonzero_fraction": f"{float(np.mean(raw_vals > 0)):.6g}" if raw_vals.size else "",
                "x_unique_values": unique_capped(x_vals) if x_vals.size else "",
                "raw_unique_values": unique_capped(raw_vals) if raw_vals.size else "",
                "x_mostly_integer_like": mostly_integer(x_vals) if x_vals.size else "unknown",
                "raw_mostly_integer_like": mostly_integer(raw_vals) if raw_vals.size else "unknown",
                "x_raw_identical": "yes" if identical else "no",
                "notes": "Target-gene vector compared where available.",
            })

        matrix_semantics_clear = "partial"
        detection_matrix = "X"
        mean_matrix = "not_selected"
        mean_allowed = "no"
        reason = "The object has no raw layer and no local matrix annotation sufficient to label X values for mean-value summaries."
        if raw_exists:
            reason = "X and raw.X were inspected, but local matrix annotation is not sufficient for mean-value summaries."

        semantic_rows = [
            {"metric": "h5ad_path", "value": str(H5AD.relative_to(ROOT))},
            {"metric": "n_obs", "value": str(adata.n_obs)},
            {"metric": "n_vars", "value": str(adata.n_vars)},
            {"metric": "x_shape", "value": x_desc["shape"]},
            {"metric": "x_dtype", "value": x_desc["dtype"]},
            {"metric": "x_storage", "value": x_desc["storage"]},
            {"metric": "x_class", "value": x_desc["class"]},
            {"metric": "raw_exists", "value": "yes" if raw_exists else "no"},
            {"metric": "raw_x_shape", "value": raw_desc["shape"]},
            {"metric": "raw_x_dtype", "value": raw_desc["dtype"]},
            {"metric": "raw_x_storage", "value": raw_desc["storage"]},
            {"metric": "raw_x_class", "value": raw_desc["class"]},
            {"metric": "layers", "value": "; ".join(map(str, adata.layers.keys())) if adata.layers.keys() else "none"},
            {"metric": "uns_keys", "value": "; ".join(map(str, adata.uns.keys()))},
            {"metric": "schema_reference", "value": str(adata.uns.get("schema_reference", ""))},
            {"metric": "schema_version", "value": str(adata.uns.get("schema_version", ""))},
            {"metric": "title", "value": str(adata.uns.get("title", ""))},
            {"metric": "var_columns", "value": "; ".join(map(str, adata.var.columns))},
            {"metric": "raw_var_columns", "value": "; ".join(map(str, adata.raw.var.columns)) if raw_exists else "not_available"},
            {"metric": "matrix_semantics_clear", "value": matrix_semantics_clear},
            {"metric": "recommended_detection_matrix", "value": detection_matrix},
            {"metric": "recommended_mean_value_matrix", "value": mean_matrix},
            {"metric": "mean_value_allowed", "value": mean_allowed},
            {"metric": "reason", "value": reason},
        ]
        write_tsv(TABLES / "human_thymus_stage3c_yayon_matrix_semantics.tsv", semantic_rows, ["metric", "value"])
        write_tsv(TABLES / "human_thymus_stage3c_yayon_x_raw_comparison.tsv", rows, [
            "gene", "x_var_id", "raw_var_id", "x_min", "x_max", "raw_min", "raw_max",
            "x_nonzero_fraction", "raw_nonzero_fraction", "x_unique_values", "raw_unique_values",
            "x_mostly_integer_like", "raw_mostly_integer_like", "x_raw_identical", "notes",
        ])
        REPORTS.mkdir(parents=True, exist_ok=True)
        REPORTS.joinpath("human_thymus_stage3c_yayon_matrix_semantics_audit.md").write_text(
            "\n".join([
                "# Yayon Human TEC Matrix Semantics Audit",
                "",
                "Search date: 2026-07-06",
                "",
                "The Yayon TEC H5AD was inspected for target-gene matrix behavior without writing large derived matrices.",
                "",
                f"- X matrix: {json.dumps(x_desc, sort_keys=True)}",
                f"- raw.X matrix: {json.dumps(raw_desc, sort_keys=True)}",
                f"- Layers: {'none' if not adata.layers.keys() else ', '.join(map(str, adata.layers.keys()))}",
                f"- Schema reference: {str(adata.uns.get('schema_reference', ''))}",
                f"- Schema version: {str(adata.uns.get('schema_version', ''))}",
                "",
                "## Conclusion",
                "",
                f"- matrix_semantics_clear: {matrix_semantics_clear}",
                f"- recommended_detection_matrix: {detection_matrix}",
                f"- recommended_mean_value_matrix: {mean_matrix}",
                f"- mean_value_allowed: {mean_allowed}",
                f"- reason: {reason}",
                "",
                "Detection summaries can proceed using `X > 0`. Mean-value summaries are not generated in Stage 3C because local value semantics remain insufficient for that output.",
            ]),
            encoding="utf-8",
        )
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
