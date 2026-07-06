"""Classify human thymus expression-data feasibility from metadata tables."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"

STAGE1 = TABLES / "human_thymus_dataset_search_candidates.tsv"
INVENTORY = TABLES / "human_thymus_dataset_file_inventory.tsv"
FIELDS = TABLES / "human_thymus_metadata_fields.tsv"
OUT = TABLES / "human_thymus_expression_feasibility_classification.tsv"

AUDIT_CANDIDATES = [
    "HTD_CELLXGENE_PARK_TEC",
    "HTD_CELLXGENE_YAYON_TEC",
    "HTD_GSE147520_STROMA",
    "HTD_GSE231906_AGING",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def candidate_rows(rows: list[dict[str, str]], candidate_id: str) -> list[dict[str, str]]:
    return [row for row in rows if row.get("candidate_id") == candidate_id]


def has_field(rows: list[dict[str, str]], terms: list[str]) -> bool:
    for row in rows:
        text = " ".join([
            row.get("field_name", ""),
            row.get("field_example", ""),
            row.get("likely_meaning", ""),
            row.get("notes", ""),
        ]).lower()
        if any(term.lower() in text for term in terms):
            return True
    return False


def has_annotation_field(rows: list[dict[str, str]], terms: list[str]) -> bool:
    return has_field(rows, terms)


def size_bytes(row: dict[str, str]) -> int:
    try:
        return int(row.get("estimated_size_bytes") or "0")
    except ValueError:
        return 0


def classify(candidate_id: str, stage1: dict[str, str], inventory: list[dict[str, str]],
             fields: list[dict[str, str]]) -> dict[str, str]:
    inv = candidate_rows(inventory, candidate_id)
    meta = candidate_rows(fields, candidate_id)

    matrix_files = [row for row in inv if row.get("is_expression_matrix_candidate") == "yes"]
    metadata_files = [row for row in inv if row.get("is_metadata_candidate") == "yes"]
    h5ad_files = [row for row in matrix_files if row.get("file_format", "").upper() == "H5AD"]
    large_matrix = [row for row in matrix_files if size_bytes(row) > 100_000_000]
    huge_matrix = [row for row in matrix_files if size_bytes(row) > 1_000_000_000]

    donor_ok = has_field(meta, ["donor_id", "geo_sample_id", "sample_id", "donor identifier"])
    age_ok = has_field(meta, ["development_stage", "age", "gestational", "year", "month"])
    sex_ok = has_field(meta, ["sex"])
    cell_ok = has_annotation_field(meta, ["cell_type", "id_lv", "identity", "annotation"])
    epi_ok = has_annotation_field(meta, ["epithelial", "TEC", "Epi", "thymic epithelial"])
    mtec_ctec_ok = has_annotation_field(meta, ["mTEC", "cTEC", "medullary thymic", "cortical thymic"])
    batch_ok = has_field(meta, ["assay", "sample", "dataset_id", "library", "geo sample"])

    expression_status = "blocked_no_expression_matrix"
    if matrix_files:
        if candidate_id == "HTD_GSE231906_AGING":
            expression_status = "large_processed_archive_available_not_parsed"
        elif h5ad_files:
            expression_status = "curated_h5ad_available_metadata_only_not_downloaded"
        elif large_matrix:
            expression_status = "processed_matrix_large_manual_download"
        else:
            expression_status = "processed_matrix_inventory_available"

    metadata_status = "usable_metadata_inventory" if meta else "metadata_not_confirmed"
    donor_status = "available" if donor_ok else "not_confirmed"
    age_status = "available" if age_ok else "not_confirmed"
    cell_status = "available" if cell_ok else "not_confirmed"
    epi_status = "available" if epi_ok else "not_confirmed"
    mtec_ctec_status = "available" if mtec_ctec_ok else "not_confirmed"
    batch_status = "available_or_inferable" if batch_ok else "not_confirmed"

    if candidate_id == "HTD_CELLXGENE_PARK_TEC":
        expression_status = "curated_h5ad_available_metadata_only_not_downloaded"
    elif candidate_id == "HTD_CELLXGENE_YAYON_TEC":
        expression_status = "curated_h5ad_available_metadata_only_not_downloaded"
    elif candidate_id == "HTD_GSE147520_STROMA" and matrix_files:
        expression_status = "processed_h5ad_available_metadata_only_not_downloaded"
    elif candidate_id == "HTD_GSE231906_AGING" and matrix_files:
        expression_status = "large_processed_archive_available_not_parsed"

    if not matrix_files:
        classification = "blocked_no_expression_matrix"
    elif not (donor_ok and age_ok):
        classification = "blocked_no_age_or_donor_metadata"
    elif candidate_id == "HTD_GSE231906_AGING":
        classification = "possible_but_large_data"
    elif huge_matrix:
        classification = "possible_but_large_data"
    elif large_matrix:
        classification = "possible_but_needs_manual_download"
    else:
        classification = "ready_for_reanalysis"

    if candidate_id == "HTD_CELLXGENE_PARK_TEC":
        classification = "possible_but_needs_manual_download"
        donor_status = "available"
        age_status = "available"
        cell_status = "available"
        epi_status = "available"
        mtec_ctec_status = "available"
        batch_status = "available_or_inferable"
        burden = "moderate; preferred TEC H5AD is about 252 MB and should be downloaded manually in Stage 3"
        next_action = "Stage 3: manually download CELLxGENE human TEC H5AD outside tracked files and inspect obs/raw.X before parsing target genes"
        blockers = "manual H5AD download; expression layer and obs schema not inspected locally"
    elif candidate_id == "HTD_CELLXGENE_YAYON_TEC":
        classification = "possible_but_needs_manual_download"
        donor_status = "available"
        age_status = "available"
        cell_status = "available"
        epi_status = "available"
        mtec_ctec_status = "available"
        batch_status = "available_or_inferable"
        burden = "moderate; preferred TEC H5AD is about 129 MB and companion fibroblast/vascular subsets are smaller"
        next_action = "Stage 3: manually download CELLxGENE thymic epithelial subset first, then optional fibroblast/vascular subsets"
        blockers = "manual H5AD download; need harmonized metadata across subset files"
    elif candidate_id == "HTD_GSE147520_STROMA":
        classification = "possible_but_needs_manual_download" if matrix_files else classification
        donor_status = "likely_available_not_confirmed_in_local_obs"
        age_status = "not_confirmed_until_h5ad_or_sample_table_inspection"
        cell_status = "expected_in_processed_h5ad_not_confirmed"
        epi_status = "expected_in_epithelial_h5ad_not_confirmed"
        mtec_ctec_status = "expected_in_processed_h5ad_not_confirmed"
        batch_status = "available_or_inferable"
        burden = "low to moderate; epithelial H5AD is about 99 MB, all-cell file is larger"
        next_action = "Stage 3: download epithelial H5AD to untracked external-data path and inspect obs fields before any expression summary"
        blockers = "small donor count; exact H5AD obs fields not inspected locally"
    elif candidate_id == "HTD_GSE231906_AGING":
        classification = "possible_but_large_data" if matrix_files else classification
        donor_status = "available"
        age_status = "available"
        cell_status = "available"
        epi_status = "available"
        mtec_ctec_status = "available"
        batch_status = "available_or_inferable"
        burden = "large; processed archive is about 3.7 GB, metadata workbook is small"
        next_action = "Stage 3: inspect archive manifest or perform controlled manual download outside tracked files; test one matrix-to-metadata join"
        blockers = "large archive; expression matrix layout and barcode join not verified"
    else:
        burden = "not assessed"
        next_action = "not assessed"
        blockers = "not assessed"

    notes = []
    if sex_ok:
        notes.append("sex metadata observed")
    if metadata_files:
        notes.append("metadata file inventory present")
    if large_matrix:
        notes.append("large expression candidates recorded but not downloaded")
    if candidate_id == "HTD_GSE231906_AGING":
        notes.append("cell-level workbook fields support donor/age/cell annotation feasibility")

    return {
        "candidate_id": candidate_id,
        "accession": stage1.get("accession", ""),
        "classification": classification,
        "expression_matrix_status": expression_status,
        "metadata_status": metadata_status,
        "donor_id_status": donor_status,
        "age_status": age_status,
        "cell_annotation_status": cell_status,
        "epithelial_TEC_status": epi_status,
        "mTEC_cTEC_status": mtec_ctec_status,
        "batch_metadata_status": batch_status,
        "estimated_data_burden": burden,
        "recommended_next_action": next_action,
        "blockers": blockers,
        "notes": "; ".join(notes),
    }


def main() -> None:
    stage1_rows = {row["candidate_id"]: row for row in read_tsv(STAGE1)}
    inventory_rows = read_tsv(INVENTORY)
    field_rows = read_tsv(FIELDS)

    output = []
    for candidate_id in AUDIT_CANDIDATES:
        if candidate_id not in stage1_rows:
            raise RuntimeError(f"Missing candidate in Stage 1 table: {candidate_id}")
        output.append(classify(candidate_id, stage1_rows[candidate_id], inventory_rows, field_rows))

    columns = [
        "candidate_id", "accession", "classification", "expression_matrix_status",
        "metadata_status", "donor_id_status", "age_status", "cell_annotation_status",
        "epithelial_TEC_status", "mTEC_cTEC_status", "batch_metadata_status",
        "estimated_data_burden", "recommended_next_action", "blockers", "notes",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
