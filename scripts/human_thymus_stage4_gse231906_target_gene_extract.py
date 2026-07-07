"""Stream LOX-family target-gene detection summaries from GSE231906."""

from __future__ import annotations

import gzip
import io
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from human_thymus_stage4_gse231906_join_pilot import core_barcode, strip_suffix
from human_thymus_stage4_gse231906_metadata_audit import (
    ARCHIVE_TAR,
    REPORTS,
    ROOT,
    TABLES,
    TARGET_GENES,
    load_metadata,
    rel,
    write_tsv,
)


CANDIDATES = TABLES / "human_thymus_stage4_gse231906_matrix_file_candidates.tsv"
JOIN_SUMMARY = TABLES / "human_thymus_stage4_gse231906_join_pilot_summary.tsv"


def text_member(tar: tarfile.TarFile, member_name: str):
    raw = tar.extractfile(member_name)
    if raw is None:
        raise RuntimeError(f"Could not open archive member {member_name}")
    if member_name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8", errors="replace")
    return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")


def read_tsv_first_columns(tar: tarfile.TarFile, member_name: str) -> list[list[str]]:
    rows: list[list[str]] = []
    with text_member(tar, member_name) as handle:
        for line in handle:
            if line.strip():
                rows.append(line.rstrip("\n").split("\t"))
    return rows


def gene_names(feature_rows: list[list[str]]) -> list[str]:
    names = []
    for row in feature_rows:
        names.append(row[1] if len(row) > 1 and row[1] else row[0])
    return names


def barcodes(barcode_rows: list[list[str]]) -> list[str]:
    return [row[0] for row in barcode_rows if row]


def filter_metadata_for_unit(metadata: pd.DataFrame, unit_key: str, files: list[str]) -> pd.DataFrame:
    text = " ".join([unit_key, *files]).lower()
    sample_hits = []
    for sample in metadata["sample_id"].dropna().astype(str).unique():
        if sample and sample.lower() in text:
            sample_hits.append(sample)
    if sample_hits:
        return metadata.loc[metadata["sample_id"].astype(str).isin(sample_hits)].copy()
    donor_hits = []
    for donor in metadata["donor_id"].dropna().astype(str).unique():
        if donor and donor.lower() in text:
            donor_hits.append(donor)
    if donor_hits:
        return metadata.loc[metadata["donor_id"].astype(str).isin(donor_hits)].copy()
    return metadata.copy()


def match_cells(unit_barcodes: list[str], metadata: pd.DataFrame) -> tuple[pd.DataFrame, str, float]:
    expr = pd.DataFrame({"expression_barcode": unit_barcodes})
    expr["exact"] = expr["expression_barcode"].map(core_barcode)
    expr["stripped"] = expr["expression_barcode"].map(strip_suffix)
    meta = metadata.copy()
    meta["exact"] = meta["cell_id"].map(core_barcode)
    meta["stripped"] = meta["cell_id"].map(strip_suffix)
    best = pd.DataFrame()
    best_name = "none"
    best_rate = -1.0
    for name in ["exact", "stripped"]:
        joined = expr.merge(meta, on=name, how="left", suffixes=("", "_metadata"))
        matched = joined.loc[joined["cell_id"].notna(), "expression_barcode"].nunique()
        rate = matched / len(expr) if len(expr) else 0.0
        if rate > best_rate:
            best = joined
            best_name = name
            best_rate = rate
    best = best.loc[best["cell_id"].notna()].drop_duplicates("expression_barcode").copy()
    return best, best_name, best_rate


def stream_target_values(
    tar: tarfile.TarFile,
    matrix_file: str,
    genes: list[str],
    unit_barcodes: list[str],
) -> tuple[dict[str, dict[str, float]], bool, str]:
    target_rows = {idx + 1: gene.upper() for idx, gene in enumerate(genes) if gene.upper() in TARGET_GENES}
    values: dict[str, dict[str, float]] = {gene: defaultdict(float) for gene in TARGET_GENES}
    integer_like = True
    orientation = "unknown"
    with text_member(tar, matrix_file) as handle:
        dims = None
        for line in handle:
            if line.startswith("%") or not line.strip():
                continue
            parts = line.strip().split()
            if dims is None:
                dims = tuple(map(int, parts[:3]))
                if dims[0] == len(genes) and dims[1] == len(unit_barcodes):
                    orientation = "genes_by_cells"
                elif dims[1] == len(genes) and dims[0] == len(unit_barcodes):
                    orientation = "cells_by_genes"
                else:
                    orientation = f"unmatched_dims_{dims[0]}x{dims[1]}"
                continue
            if len(parts) < 3:
                continue
            i, j = int(parts[0]), int(parts[1])
            value = float(parts[2])
            if abs(value - round(value)) > 1e-8:
                integer_like = False
            if orientation == "genes_by_cells" and i in target_rows:
                if 1 <= j <= len(unit_barcodes):
                    values[target_rows[i]][unit_barcodes[j - 1]] += value
            elif orientation == "cells_by_genes" and j in target_rows:
                if 1 <= i <= len(unit_barcodes):
                    values[target_rows[j]][unit_barcodes[i - 1]] += value
    return values, integer_like, orientation


def aggregate_unit(unit: pd.Series, metadata: pd.DataFrame, tar: tarfile.TarFile) -> tuple[list[dict[str, object]], dict[str, object], bool]:
    feature_rows = read_tsv_first_columns(tar, unit["features_file"])
    barcode_rows = read_tsv_first_columns(tar, unit["barcodes_file"])
    genes = gene_names(feature_rows)
    unit_barcodes = barcodes(barcode_rows)
    unit_meta = filter_metadata_for_unit(metadata, unit["unit_key"], [unit["matrix_file"], unit["features_file"], unit["barcodes_file"]])
    joined, strategy, rate = match_cells(unit_barcodes, unit_meta)
    values, integer_like, orientation = stream_target_values(tar, unit["matrix_file"], genes, unit_barcodes)
    rows: list[dict[str, object]] = []
    if not joined.empty:
        group_cols = ["donor_id", "sample_id", "age_years", "age_group", "sex", "broad_compartment", "fine_cell_type", "target_compartments"]
        for keys, group in joined.groupby(group_cols, dropna=False):
            base = dict(zip(group_cols, keys))
            cell_ids = group["expression_barcode"].astype(str).tolist()
            for gene in TARGET_GENES:
                vals = np.array([values[gene].get(cell, 0.0) for cell in cell_ids], dtype=float)
                rows.append({
                    **base,
                    "unit_key": unit["unit_key"],
                    "gene": gene,
                    "n_cells": int(len(vals)),
                    "n_detected_cells": int((vals > 0).sum()),
                    "detection_fraction": float((vals > 0).mean()) if len(vals) else np.nan,
                    "sum_value": float(vals.sum()),
                    "mean_value": float(vals.mean()) if len(vals) else np.nan,
                })
    summary = {
        "unit_key": unit["unit_key"],
        "matrix_file": unit["matrix_file"],
        "features_file": unit["features_file"],
        "barcodes_file": unit["barcodes_file"],
        "n_expression_cells": len(unit_barcodes),
        "n_matched_expression_cells": int(joined["expression_barcode"].nunique()) if not joined.empty else 0,
        "match_rate": f"{rate:.6f}",
        "join_strategy": strategy,
        "n_features": len(genes),
        "target_genes_present": ";".join([gene for gene in TARGET_GENES if gene in set(map(str.upper, genes))]),
        "matrix_orientation": orientation,
        "matrix_integer_like_for_targets": "yes" if integer_like else "no",
    }
    count_like = integer_like and any(token in str(unit["matrix_file"]).lower() for token in ["raw", "count", "matrix.mtx"])
    return rows, summary, count_like


def weighted_summary(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for keys, group in df.groupby(group_cols + ["gene"], dropna=False):
        detected = group["n_detected_cells"].sum()
        cells = group["n_cells"].sum()
        rows.append({
            **dict(zip(group_cols + ["gene"], keys)),
            "n_cells": int(cells),
            "n_detected_cells": int(detected),
            "detection_fraction": float(detected / cells) if cells else np.nan,
            "n_donor_sample_units": int(group[["donor_id", "sample_id"]].drop_duplicates().shape[0]) if {"donor_id", "sample_id"}.issubset(group.columns) else "",
        })
    return pd.DataFrame(rows)


def main() -> int:
    if not JOIN_SUMMARY.exists():
        raise RuntimeError("Join pilot summary is missing.")
    join_summary = pd.read_csv(JOIN_SUMMARY, sep="\t")
    if join_summary.empty or str(join_summary.loc[0, "pilot_pass"]).lower() != "yes":
        print("join_pilot_pass=no; target extraction skipped")
        return 0
    candidates = pd.read_csv(CANDIDATES, sep="\t").fillna("")
    units = candidates.loc[candidates["complete_mtx_triplet"].eq("yes")].copy()
    if units.empty:
        raise RuntimeError("No complete MTX units available for target extraction.")
    metadata, _ = load_metadata()
    all_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    count_like_all = True
    with tarfile.open(ARCHIVE_TAR, "r") as tar:
        for _, unit in units.iterrows():
            rows, summary, count_like = aggregate_unit(unit, metadata, tar)
            all_rows.extend(rows)
            summaries.append(summary)
            count_like_all = count_like_all and count_like
    if not all_rows:
        raise RuntimeError("No matched target-gene rows were generated.")
    df = pd.DataFrame(all_rows)
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_target_gene_cell_metadata_join_summary.tsv",
        summaries,
        ["unit_key", "matrix_file", "features_file", "barcodes_file", "n_expression_cells", "n_matched_expression_cells", "match_rate", "join_strategy", "n_features", "target_genes_present", "matrix_orientation", "matrix_integer_like_for_targets"],
    )
    donor = weighted_summary(df, ["donor_id", "sample_id", "age_years", "age_group", "sex", "broad_compartment", "fine_cell_type", "target_compartments"])
    fine = weighted_summary(df, ["fine_cell_type", "target_compartments"])
    broad = weighted_summary(df, ["broad_compartment"])
    age = weighted_summary(df, ["age_group", "broad_compartment", "fine_cell_type", "target_compartments"])
    write_tsv(TABLES / "human_thymus_stage4_gse231906_lox_detection_by_donor_subtype.tsv", donor.to_dict("records"), donor.columns.tolist())
    write_tsv(TABLES / "human_thymus_stage4_gse231906_lox_detection_by_age_subtype.tsv", age.to_dict("records"), age.columns.tolist())
    write_tsv(TABLES / "human_thymus_stage4_gse231906_lox_detection_by_fine_celltype.tsv", fine.to_dict("records"), fine.columns.tolist())
    write_tsv(TABLES / "human_thymus_stage4_gse231906_lox_detection_by_broad_compartment.tsv", broad.to_dict("records"), broad.columns.tolist())
    expression_note = "Mean-value summaries withheld because matrix semantics are not clear enough for expression-value interpretation."
    if count_like_all:
        expr = donor.copy()
        expr["mean_value_note"] = "Matrix values were integer-like in target rows; use as count-like descriptive context only."
        write_tsv(TABLES / "human_thymus_stage4_gse231906_lox_expression_by_donor_subtype.tsv", expr.to_dict("records"), expr.columns.tolist())
        expression_note = "Count-like descriptive mean-value summaries were generated, but detection fractions remain the primary context."
    print(f"target_units_processed={len(units)}")
    print(f"donor_subtype_rows={len(donor)}")
    print(f"mean_values_generated={'yes' if count_like_all else 'no'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
