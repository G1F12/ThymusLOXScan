"""Build a non-extracting manifest for the GSE231906 raw archive."""

from __future__ import annotations

import re
import tarfile
from pathlib import Path

import pandas as pd

from human_thymus_stage4_gse231906_metadata_audit import (
    ARCHIVE_TAR,
    ARCHIVE_URL,
    DATA_DIR,
    REPORTS,
    TABLES,
    download_if_missing,
    free_gb,
    rel,
    write_disk_report,
    write_tsv,
)


def classify_member(name: str) -> str:
    lower = name.lower()
    if lower.endswith((".mtx", ".mtx.gz")):
        return "matrix_mtx"
    if "feature" in lower and lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "features_tsv"
    if "gene" in lower and lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "features_tsv"
    if "barcode" in lower and lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "barcodes_tsv"
    if "metadata" in lower and lower.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz")):
        return "metadata_csv"
    if any(token in lower for token in ["count", "counts"]) and lower.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz")):
        return "count_csv"
    if any(token in lower for token in ["expression", "expr", "matrix"]) and lower.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz")):
        return "expression_csv"
    if lower.endswith((".gz", ".zip", ".bz2", ".xz")):
        return "compressed_file"
    return "unknown"


def unit_key(name: str, classification: str) -> str:
    path = Path(name)
    base = path.name
    lower = base.lower()
    for suffix in [
        ".gz", ".mtx", ".tsv", ".txt", ".csv",
    ]:
        if lower.endswith(suffix):
            base = base[: -len(suffix)]
            lower = base.lower()
    base = re.sub(r"(^|[_\-.])(matrix|features|genes|barcodes|counts?|expression|expr)$", "", base, flags=re.I)
    base = re.sub(r"[_\-.]+$", "", base)
    if not base or base.lower() in {"matrix", "features", "genes", "barcodes"}:
        return path.parent.as_posix()
    return (path.parent / base).as_posix()


def manifest_rows() -> list[dict[str, object]]:
    download_if_missing(ARCHIVE_URL, ARCHIVE_TAR)
    rows: list[dict[str, object]] = []
    with tarfile.open(ARCHIVE_TAR, "r") as tar:
        for member in tar:
            classification = classify_member(member.name)
            rows.append({
                "member_name": member.name,
                "size_bytes": int(member.size),
                "is_file": "yes" if member.isfile() else "no",
                "classification": classification,
                "unit_key": unit_key(member.name, classification) if member.isfile() else "",
            })
    return rows


def candidate_rows(manifest: list[dict[str, object]]) -> list[dict[str, object]]:
    df = pd.DataFrame(manifest)
    df = df.loc[df["is_file"].eq("yes")].copy()
    rows: list[dict[str, object]] = []
    for key, group in df.groupby("unit_key", dropna=False):
        files = {kind: group.loc[group["classification"].eq(kind), "member_name"].tolist() for kind in group["classification"].unique()}
        matrices = files.get("matrix_mtx", [])
        features = files.get("features_tsv", [])
        barcodes = files.get("barcodes_tsv", [])
        count_csv = files.get("count_csv", [])
        expr_csv = files.get("expression_csv", [])
        complete = bool(matrices and features and barcodes)
        if complete or count_csv or expr_csv:
            rows.append({
                "unit_key": key,
                "unit_type": "mtx_triplet" if complete else "dense_table_candidate",
                "complete_mtx_triplet": "yes" if complete else "no",
                "matrix_file": matrices[0] if matrices else "",
                "features_file": features[0] if features else "",
                "barcodes_file": barcodes[0] if barcodes else "",
                "dense_candidate_file": (count_csv + expr_csv)[0] if count_csv or expr_csv else "",
                "total_unit_size_bytes": int(group["size_bytes"].sum()),
                "n_files_in_unit": int(len(group)),
            })
    return sorted(rows, key=lambda row: (row["complete_mtx_triplet"] != "yes", row["total_unit_size_bytes"]))


def write_report(manifest: list[dict[str, object]], candidates: list[dict[str, object]]) -> None:
    class_counts = pd.DataFrame(manifest).groupby("classification").size().reset_index(name="n_files")
    lines = [
        "# GSE231906 Archive Manifest Report",
        "",
        "## Archive",
        "",
        f"- Local archive path: `{rel(ARCHIVE_TAR)}`",
        f"- Archive exists: {'yes' if ARCHIVE_TAR.exists() else 'no'}",
        f"- Archive size: {ARCHIVE_TAR.stat().st_size if ARCHIVE_TAR.exists() else 0} bytes",
        f"- Free space after manifest step: {free_gb(DATA_DIR):.2f} GB",
        "",
        "## Manifest Summary",
        "",
        f"- Files/directories listed: {len(manifest)}",
        f"- Candidate expression units: {len(candidates)}",
        "- Full extraction performed: no",
        "",
        "## File Classes",
        "",
    ]
    for row in class_counts.to_dict("records"):
        lines.append(f"- {row['classification']}: {row['n_files']}")
    if candidates:
        first = candidates[0]
        lines.extend([
            "",
            "## Smallest Candidate Unit",
            "",
            f"- Unit key: `{first['unit_key']}`",
            f"- Unit type: {first['unit_type']}",
            f"- Complete MTX triplet: {first['complete_mtx_triplet']}",
            f"- Unit size: {first['total_unit_size_bytes']} bytes",
        ])
    else:
        lines.extend(["", "No complete expression unit was identified from member names alone."])
    (REPORTS / "human_thymus_stage4_gse231906_archive_manifest_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if not write_disk_report():
        return 0
    rows = manifest_rows()
    candidates = candidate_rows(rows)
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_archive_manifest.tsv",
        rows,
        ["member_name", "size_bytes", "is_file", "classification", "unit_key"],
    )
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_matrix_file_candidates.tsv",
        candidates,
        ["unit_key", "unit_type", "complete_mtx_triplet", "matrix_file", "features_file", "barcodes_file", "dense_candidate_file", "total_unit_size_bytes", "n_files_in_unit"],
    )
    write_report(rows, candidates)
    print(f"archive_manifest_rows={len(rows)}")
    print(f"candidate_units={len(candidates)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
