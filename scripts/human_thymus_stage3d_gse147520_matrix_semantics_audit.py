"""Audit matrix semantics for GSE147520 epithelial target-gene summaries."""

from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse


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


def find_match(var: Any, gene: str) -> str:
    if gene in var.index.astype(str):
        return gene
    for column in var.columns:
        matches = var.index[var[column].astype(str) == gene].astype(str).tolist()
        if matches:
            return matches[0]
    return ""


def main() -> None:
    path = ensure_h5ad()
    adata = ad.read_h5ad(path, backed="r")
    try:
        raw_exists = adata.raw is not None
        x_desc = desc(adata.X)
        raw_desc = desc(adata.raw.X) if raw_exists else {"shape": "not_available", "dtype": "not_available", "storage": "not_available", "class": "not_available"}
        rows: list[dict[str, str]] = []
        present_any = False
        for gene in TARGET_GENES:
            x_var = find_match(adata.var, gene)
            raw_var = find_match(adata.raw.var, gene) if raw_exists else ""
            present_any = present_any or bool(x_var or raw_var)
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
                "x_mostly_integer_like": mostly_integer(x_vals) if x_vals.size else "not_available",
                "raw_mostly_integer_like": mostly_integer(raw_vals) if raw_vals.size else "not_available",
                "x_raw_identical": "yes" if identical else "no",
                "notes": "Target gene not present in available feature set." if not x_var and not raw_var else "Target-gene vector inspected where available.",
            })
        semantics_clear = "partial"
        raw_has_targets = any(row["raw_var_id"] for row in rows)
        x_has_targets = any(row["x_var_id"] for row in rows)
        detection_matrix = "raw.X" if raw_has_targets else ("X" if x_has_targets else "not_selected")
        mean_matrix = "not_selected"
        mean_allowed = "no"
        reason = "The target LOX-family genes are not present in the available epithelial H5AD feature set; detection and mean-value summaries are therefore not generated for this dataset."
        if raw_has_targets and not x_has_targets:
            reason = "Target genes are absent from compact X but present in raw.X; detection summaries use raw.X > 0, while mean-value summaries are not generated because value semantics remain partial."
        elif present_any:
            reason = "At least one target gene is present, but local matrix semantics are insufficient for mean-value summaries."

        semantic_rows = [
            {"metric": "h5ad_path", "value": str(path.relative_to(ROOT))},
            {"metric": "h5ad_gz_path", "value": str(H5AD_GZ.relative_to(ROOT)) if H5AD_GZ.exists() else "not_available"},
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
            {"metric": "matrix_semantics_clear", "value": semantics_clear},
            {"metric": "recommended_detection_matrix", "value": detection_matrix},
            {"metric": "recommended_mean_value_matrix", "value": mean_matrix},
            {"metric": "mean_value_allowed", "value": mean_allowed},
            {"metric": "reason", "value": reason},
        ]
        write_tsv(TABLES / "human_thymus_stage3d_gse147520_matrix_semantics.tsv", semantic_rows, ["metric", "value"])
        write_tsv(TABLES / "human_thymus_stage3d_gse147520_x_raw_comparison.tsv", rows, [
            "gene", "x_var_id", "raw_var_id", "x_min", "x_max", "raw_min", "raw_max",
            "x_nonzero_fraction", "raw_nonzero_fraction", "x_unique_values", "raw_unique_values",
            "x_mostly_integer_like", "raw_mostly_integer_like", "x_raw_identical", "notes",
        ])
        REPORTS.mkdir(parents=True, exist_ok=True)
        REPORTS.joinpath("human_thymus_stage3d_gse147520_matrix_semantics_audit.md").write_text(
            "\n".join([
                "# GSE147520 Human Thymic Epithelial Matrix Semantics Audit",
                "",
                "Search date: 2026-07-06",
                "",
                "The GSE147520 epithelial H5AD was inspected for target-gene matrix behavior without writing large derived matrices.",
                "",
                f"- X matrix: {x_desc}",
                f"- raw.X matrix: {raw_desc}",
                f"- Layers: {'none' if not adata.layers.keys() else ', '.join(map(str, adata.layers.keys()))}",
                "",
                "## Conclusion",
                "",
                f"- matrix_semantics_clear: {semantics_clear}",
                f"- recommended_detection_matrix: {detection_matrix}",
                f"- recommended_mean_value_matrix: {mean_matrix}",
                f"- mean_value_allowed: {mean_allowed}",
                f"- reason: {reason}",
                "",
                "Target-gene detection can proceed from `raw.X > 0` when target genes are present only in raw.var. Mean-value summaries are not generated in Stage 3D.",
            ]),
            encoding="utf-8",
        )
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
