"""Extract one minimal GSE231906 pilot expression unit."""

from __future__ import annotations

import shutil
import tarfile
from pathlib import Path
import re
import gzip
import io

import pandas as pd

from human_thymus_stage4_gse231906_metadata_audit import (
    ARCHIVE_TAR,
    PILOT_DIR,
    REPORTS,
    TABLES,
    ARCHIVE_TAR,
    load_metadata,
    rel,
    write_tsv,
)


CANDIDATES = TABLES / "human_thymus_stage4_gse231906_matrix_file_candidates.tsv"


def sample_token(unit_key: str) -> str:
    match = re.search(r"_(donor[0-9]+(?:-[0-9]+)?)$", str(unit_key))
    return match.group(1) if match else str(unit_key)


def core_barcode(value: object) -> str:
    text = str(value).strip()
    if "_donor" in text:
        text = text.split("_donor", 1)[0]
    return text


def estimate_match_rates(candidates: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    if not ARCHIVE_TAR.exists() or candidates.empty:
        candidates["estimated_match_rate"] = 0.0
        return candidates
    meta = metadata.copy()
    meta["core_barcode"] = meta["cell_id"].map(core_barcode)
    rates: list[float] = []
    with tarfile.open(ARCHIVE_TAR, "r") as tar:
        for _, row in candidates.iterrows():
            sample = sample_token(row["unit_key"])
            sample_meta = meta.loc[meta["sample_id"].astype(str).eq(sample)]
            barcode_set = set(sample_meta["core_barcode"].dropna().astype(str))
            if not barcode_set:
                rates.append(0.0)
                continue
            raw = tar.extractfile(row["barcodes_file"])
            if raw is None:
                rates.append(0.0)
                continue
            handle = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8", errors="replace") if str(row["barcodes_file"]).endswith(".gz") else io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            barcodes = [line.strip().split("\t")[0] for line in handle if line.strip()]
            matched = sum(1 for barcode in barcodes if core_barcode(barcode) in barcode_set)
            rates.append(matched / len(barcodes) if barcodes else 0.0)
    candidates = candidates.copy()
    candidates["estimated_match_rate"] = rates
    return candidates


def safe_extract_member(tar: tarfile.TarFile, member_name: str, dest_dir: Path) -> Path:
    member = tar.getmember(member_name)
    if not member.isfile():
        raise RuntimeError(f"Refusing to extract non-file archive member: {member_name}")
    target = (dest_dir / member_name).resolve()
    root = dest_dir.resolve()
    if root not in target.parents and target != root:
        raise RuntimeError(f"Refusing path outside pilot directory: {member_name}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with tar.extractfile(member) as src, target.open("wb") as dst:
        if src is None:
            raise RuntimeError(f"Could not open archive member: {member_name}")
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    return target


def choose_unit(candidates: pd.DataFrame) -> pd.Series | None:
    complete = candidates.loc[candidates["complete_mtx_triplet"].eq("yes")].copy()
    if complete.empty:
        dense = candidates.loc[candidates["dense_candidate_file"].astype(str).ne("")].copy()
        if dense.empty:
            return None
        dense = dense.sort_values(["total_unit_size_bytes", "unit_key"])
        return dense.iloc[0]
    metadata, _ = load_metadata()
    target_counts = (
        metadata.loc[
            metadata["target_compartments"].str.contains("epithelial|mtec|ctec|post_AIRE|immature_TEC", case=False, na=False)
        ]
        .groupby("sample_id", dropna=False)
        .size()
    )
    target_sample_set = set(target_counts.loc[target_counts >= 100].index.astype(str))
    if not target_sample_set:
        target_sample_set = set(target_counts.index.astype(str))
    tec_mask = complete["unit_key"].map(sample_token).isin(target_sample_set)
    tec = complete.loc[tec_mask].copy()
    selected = estimate_match_rates(tec if not tec.empty else complete, metadata)
    passing = selected.loc[selected["estimated_match_rate"] >= 0.80].copy()
    selected = passing if not passing.empty else selected
    selected = selected.sort_values(["total_unit_size_bytes", "estimated_match_rate", "unit_key"], ascending=[True, False, True])
    return selected.iloc[0]


def failure_report(reason: str) -> None:
    (REPORTS / "human_thymus_stage4_gse231906_pilot_extract_failed.md").write_text(
        "\n".join([
            "# GSE231906 Pilot Extraction Failed",
            "",
            reason,
            "",
            "No expression analysis was run from this archive because no safe pilot unit could be identified.",
        ]) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    if not ARCHIVE_TAR.exists():
        failure_report(f"Archive missing at `{rel(ARCHIVE_TAR)}`. Run the archive manifest step first.")
        return 0
    if not CANDIDATES.exists():
        failure_report(f"Candidate table missing at `{rel(CANDIDATES)}`. Run the archive manifest step first.")
        return 0
    candidates = pd.read_csv(CANDIDATES, sep="\t").fillna("")
    unit = choose_unit(candidates)
    if unit is None:
        failure_report("No complete MTX triplet or dense count/expression table candidate was identified.")
        return 0

    PILOT_DIR.mkdir(parents=True, exist_ok=True)
    members = []
    if unit["complete_mtx_triplet"] == "yes":
        members = [unit["matrix_file"], unit["features_file"], unit["barcodes_file"]]
    elif unit["dense_candidate_file"]:
        members = [unit["dense_candidate_file"]]
    members = [member for member in members if isinstance(member, str) and member]

    rows = []
    with tarfile.open(ARCHIVE_TAR, "r") as tar:
        for member in members:
            path = safe_extract_member(tar, member, PILOT_DIR)
            rows.append({
                "unit_key": unit["unit_key"],
                "member_name": member,
                "local_path": rel(path),
                "size_bytes": path.stat().st_size,
                "extraction_role": (
                    "matrix" if member == unit.get("matrix_file", "") else
                    "features" if member == unit.get("features_file", "") else
                    "barcodes" if member == unit.get("barcodes_file", "") else
                    "dense_candidate"
                ),
            })

    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_pilot_extracted_files.tsv",
        rows,
        ["unit_key", "member_name", "local_path", "size_bytes", "extraction_role"],
    )
    (REPORTS / "human_thymus_stage4_gse231906_pilot_extract_report.md").write_text(
        "\n".join([
            "# GSE231906 Pilot Extraction Report",
            "",
            f"- Selected unit key: `{unit['unit_key']}`",
            f"- Unit type: {unit['unit_type']}",
            f"- Selection rule: smallest complete TEC/thymus-labeled unit if present; otherwise smallest complete expression unit.",
            f"- Extracted file count: {len(rows)}",
            f"- Extraction directory: `{rel(PILOT_DIR)}`",
            "- Full archive extraction performed: no",
            "",
            "The extracted expression files remain under the untracked external-data path and are not intended for staging.",
        ]) + "\n",
        encoding="utf-8",
    )
    print(f"pilot_extracted_files={len(rows)}")
    print(f"pilot_unit={unit['unit_key']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
