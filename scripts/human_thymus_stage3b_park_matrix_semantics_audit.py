"""Audit X/raw.X semantics for Park human TEC target-gene summaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data" / "external" / "human_thymus" / "cellxgene_park_tec" / "park_tec.h5ad"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"

TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]

SEMANTICS_TSV = TABLES / "human_thymus_stage3b_park_matrix_semantics.tsv"
COMPARISON_TSV = TABLES / "human_thymus_stage3b_park_x_raw_comparison.tsv"
REPORT = REPORTS / "human_thymus_stage3b_park_matrix_semantics_audit.md"


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def matrix_descriptor(matrix: Any) -> dict[str, str]:
    return {
        "shape": str(getattr(matrix, "shape", "")),
        "dtype": str(getattr(matrix, "dtype", "")),
        "storage": "sparse" if sparse.issparse(matrix) or "sparse" in type(matrix).__module__.lower() else "dense",
        "class": f"{type(matrix).__module__}.{type(matrix).__name__}",
    }


def as_vector(matrix: Any, var_idx: int) -> np.ndarray:
    values = matrix[:, var_idx]
    if sparse.issparse(values):
        return np.asarray(values.toarray()).ravel()
    if hasattr(values, "to_memory"):
        values = values.to_memory()
    if sparse.issparse(values):
        return np.asarray(values.toarray()).ravel()
    return np.asarray(values).ravel()


def mostly_integer(values: np.ndarray) -> str:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return "unknown"
    return "yes" if float(np.mean(np.isclose(finite, np.round(finite), atol=1e-6))) >= 0.995 else "no"


def unique_count_capped(values: np.ndarray, cap: int = 5000) -> str:
    if values.size > cap:
        idx = np.linspace(0, values.size - 1, cap, dtype=int)
        return f"{len(np.unique(values[idx]))}_in_{cap}_sample"
    return str(len(np.unique(values)))


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
        raise FileNotFoundError(f"Missing Park TEC H5AD: {H5AD}")

    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        raw_exists = adata.raw is not None
        x_desc = matrix_descriptor(adata.X)
        raw_desc = matrix_descriptor(adata.raw.X) if raw_exists else {
            "shape": "not_available",
            "dtype": "not_available",
            "storage": "not_available",
            "class": "not_available",
        }
        uns_keys = "; ".join(map(str, adata.uns.keys()))
        schema_reference = str(adata.uns.get("schema_reference", ""))
        schema_version = str(adata.uns.get("schema_version", ""))
        title = str(adata.uns.get("title", ""))

        x_gene_map = gene_map(adata.var)
        raw_gene_map = gene_map(adata.raw.var) if raw_exists else {}
        comparison_rows: list[dict[str, str]] = []
        any_difference = False
        all_x_integer = True
        all_raw_integer = True
        for gene in TARGET_GENES:
            x_var = x_gene_map.get(gene, "")
            raw_var = raw_gene_map.get(gene, "")
            x_values = as_vector(adata.X, int(adata.var_names.get_loc(x_var))) if x_var else np.array([])
            raw_values = (
                as_vector(adata.raw.X, int(adata.raw.var_names.get_loc(raw_var)))
                if raw_exists and raw_var else np.array([])
            )
            identical = bool(x_values.size and raw_values.size and np.array_equal(x_values, raw_values))
            any_difference = any_difference or (x_values.size > 0 and raw_values.size > 0 and not identical)
            x_integer = mostly_integer(x_values) if x_values.size else "unknown"
            raw_integer = mostly_integer(raw_values) if raw_values.size else "unknown"
            all_x_integer = all_x_integer and x_integer == "yes"
            all_raw_integer = all_raw_integer and raw_integer == "yes"
            comparison_rows.append({
                "gene": gene,
                "x_var_id": x_var,
                "raw_var_id": raw_var,
                "x_min": f"{float(np.min(x_values)):.6g}" if x_values.size else "",
                "x_max": f"{float(np.max(x_values)):.6g}" if x_values.size else "",
                "raw_min": f"{float(np.min(raw_values)):.6g}" if raw_values.size else "",
                "raw_max": f"{float(np.max(raw_values)):.6g}" if raw_values.size else "",
                "x_nonzero_fraction": f"{float(np.mean(x_values > 0)):.6g}" if x_values.size else "",
                "raw_nonzero_fraction": f"{float(np.mean(raw_values > 0)):.6g}" if raw_values.size else "",
                "x_unique_values": unique_count_capped(x_values) if x_values.size else "",
                "raw_unique_values": unique_count_capped(raw_values) if raw_values.size else "",
                "x_mostly_integer_like": x_integer,
                "raw_mostly_integer_like": raw_integer,
                "x_raw_identical": "yes" if identical else "no",
                "notes": "Target-gene vectors are identical in X and raw.X." if identical else "Target-gene vectors differ or raw vector is unavailable.",
            })

        matrix_semantics_clear = "partial"
        recommended_detection_matrix = "X"
        recommended_mean_value_matrix = "not_selected"
        mean_value_allowed = "no"
        reason = (
            "X and raw.X are present and target genes are available, but the H5AD does not provide enough "
            "local matrix annotation to label values as raw counts or normalized expression for mean-value summaries."
        )
        if all_x_integer and all_raw_integer and raw_exists and not any_difference:
            reason = (
                "X and raw.X target-gene vectors are identical and mostly integer-like, but the object still lacks "
                "explicit local matrix semantics. Detection summaries are the safest Stage 3B output."
            )

        semantics_rows = [
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
            {"metric": "uns_keys", "value": uns_keys},
            {"metric": "schema_reference", "value": schema_reference},
            {"metric": "schema_version", "value": schema_version},
            {"metric": "title", "value": title},
            {"metric": "var_columns", "value": "; ".join(map(str, adata.var.columns))},
            {"metric": "raw_var_columns", "value": "; ".join(map(str, adata.raw.var.columns)) if raw_exists else "not_available"},
            {"metric": "matrix_semantics_clear", "value": matrix_semantics_clear},
            {"metric": "recommended_detection_matrix", "value": recommended_detection_matrix},
            {"metric": "recommended_mean_value_matrix", "value": recommended_mean_value_matrix},
            {"metric": "mean_value_allowed", "value": mean_value_allowed},
            {"metric": "reason", "value": reason},
        ]

        write_tsv(SEMANTICS_TSV, semantics_rows, ["metric", "value"])
        write_tsv(COMPARISON_TSV, comparison_rows, [
            "gene", "x_var_id", "raw_var_id", "x_min", "x_max", "raw_min", "raw_max",
            "x_nonzero_fraction", "raw_nonzero_fraction", "x_unique_values", "raw_unique_values",
            "x_mostly_integer_like", "raw_mostly_integer_like", "x_raw_identical", "notes",
        ])

        REPORTS.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(
            "\n".join([
                "# Park Human TEC Matrix Semantics Audit",
                "",
                "Search date: 2026-07-06",
                "",
                "The Park TEC H5AD was inspected for target-gene matrix behavior without writing large derived matrices.",
                "",
                f"- X matrix: {json.dumps(x_desc, sort_keys=True)}",
                f"- raw.X matrix: {json.dumps(raw_desc, sort_keys=True)}",
                f"- Layers: {'none' if not adata.layers.keys() else ', '.join(map(str, adata.layers.keys()))}",
                f"- Schema reference: {schema_reference}",
                f"- Schema version: {schema_version}",
                f"- Target genes compared: {', '.join(TARGET_GENES)}",
                "",
                "## Conclusion",
                "",
                f"- matrix_semantics_clear: {matrix_semantics_clear}",
                f"- recommended_detection_matrix: {recommended_detection_matrix}",
                f"- recommended_mean_value_matrix: {recommended_mean_value_matrix}",
                f"- mean_value_allowed: {mean_value_allowed}",
                f"- reason: {reason}",
                "",
                "Detection summaries can proceed using `X > 0`. Mean-value summaries are not generated in Stage 3B because value semantics remain insufficiently documented locally.",
            ]),
            encoding="utf-8",
        )
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
