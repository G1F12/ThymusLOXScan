"""Stage 4 GSE231906 metadata audit for human thymus LOX-family context."""

from __future__ import annotations

import csv
import re
import shutil
import urllib.request
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external" / "human_thymus" / "GSE231906"
PILOT_DIR = DATA_DIR / "pilot_extract"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "results" / "figures" / "human_thymus_stage4_gse231906"

METADATA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/"
    "GSE231906_cell-level_metadata.xlsx"
)
ARCHIVE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/GSE231906_RAW.tar"
METADATA_XLSX = DATA_DIR / "GSE231906_cell-level_metadata.xlsx"
ARCHIVE_TAR = DATA_DIR / "GSE231906_RAW.tar"

TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]
MIN_ARCHIVE_FREE_GB = 12.0
REQUIRED_NO_OVERCLAIM = (
    "This GSE231906 analysis provides aged-human thymus transcript-level context only. "
    "It does not establish human conservation of the mouse GSE240016 mTEC1 Loxl2 candidate signal, "
    "does not validate the mouse result, and does not provide mechanism, protein-level evidence, "
    "functional evidence, LOX activity evidence, or therapeutic relevance."
)


def ensure_dirs() -> None:
    for path in [DATA_DIR, TABLES, REPORTS, FIGURES]:
        path.mkdir(parents=True, exist_ok=True)


def rel(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def write_tsv(path: Path, rows: Iterable[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: "" if pd.isna(row.get(column, "")) else row.get(column, "") for column in columns})


def download_if_missing(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    temp = dest.with_suffix(dest.suffix + ".part")
    if temp.exists():
        temp.unlink()
    urllib.request.urlretrieve(url, temp)
    temp.replace(dest)
    return dest


def free_gb(path: Path) -> float:
    path.mkdir(parents=True, exist_ok=True)
    return shutil.disk_usage(path).free / 1024**3


def write_disk_report() -> bool:
    ensure_dirs()
    gb = free_gb(DATA_DIR)
    ok = gb >= MIN_ARCHIVE_FREE_GB
    lines = [
        "# GSE231906 Disk-Space Check",
        "",
        f"- Data directory: `{rel(DATA_DIR)}`",
        f"- Free disk space: {gb:.2f} GB",
        f"- Minimum requested before archive download: {MIN_ARCHIVE_FREE_GB:.2f} GB",
        f"- Disk space sufficient for archive step: {'yes' if ok else 'no'}",
    ]
    (REPORTS / "human_thymus_stage4_gse231906_disk_space_check.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if not ok:
        (REPORTS / "human_thymus_stage4_gse231906_blocked_insufficient_disk.md").write_text(
            "\n".join([
                "# GSE231906 Blocked: Insufficient Disk",
                "",
                f"Free disk space was {gb:.2f} GB, below the requested 12 GB guard. The raw archive was not downloaded.",
                "",
                "Manual plan for a larger disk:",
                "",
                "1. Create `data/external/human_thymus/GSE231906/` on a filesystem with more free space.",
                "2. Download `GSE231906_cell-level_metadata.xlsx` first and rerun the metadata audit.",
                "3. Download `GSE231906_RAW.tar` only under that GSE231906 external-data directory.",
                "4. Rerun the archive manifest script before any extraction.",
                "5. Extract only a minimal pilot expression unit for matrix-to-metadata linkage testing.",
            ]) + "\n",
            encoding="utf-8",
        )
    return ok


def norm_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def parse_age_years(value: object) -> float:
    text = norm_text(value).lower()
    if not text:
        return np.nan
    match = re.search(r"(\d+(?:\.\d+)?)", text)
    if not match:
        return np.nan
    age = float(match.group(1))
    if "week" in text:
        return age / 52.1775
    if "month" in text or re.search(r"\bmo\b", text):
        return age / 12.0
    if "day" in text:
        return age / 365.25
    if "fetal" in text or "gest" in text:
        return age / 52.1775
    return age


def age_group(age: float) -> str:
    if pd.isna(age):
        return "age_unknown"
    if age < 0:
        return "age_unknown"
    if age < 1:
        return "fetal_or_infant_lt1y"
    if age < 18:
        return "pediatric_1_17y"
    if age < 40:
        return "young_adult_18_39y"
    if age < 65:
        return "adult_40_64y"
    return "older_adult_65plus"


def choose_column(columns: list[str], preferred: list[str], contains: list[str], avoid: list[str] | None = None) -> str:
    avoid = avoid or []
    lower_map = {column.lower(): column for column in columns}
    for name in preferred:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    for column in columns:
        lower = column.lower()
        if any(token in lower for token in contains) and not any(token in lower for token in avoid):
            return column
    return ""


def infer_columns(df: pd.DataFrame) -> dict[str, str]:
    columns = [str(column) for column in df.columns]
    return {
        "cell_id": choose_column(columns, ["barcode", "cell_id", "cell", "CellID"], ["barcode", "cell"], ["type", "identity"]),
        "donor_id": choose_column(columns, ["donor_id", "donor", "Donor"], ["donor", "individual"]),
        "sample_id": choose_column(columns, ["sample_id", "geo_sample_id", "sample", "GSM"], ["sample", "geo"]),
        "age": choose_column(columns, ["Age", "age", "age_years"], ["age"]),
        "sex": choose_column(columns, ["Sex", "sex", "gender"], ["sex", "gender"]),
        "broad": choose_column(columns, ["id_lv1", "broad_cell_type", "major_cell_type", "cell_type"], ["lv1", "broad", "major"]),
        "fine": choose_column(
            columns,
            ["id_lv4", "id_lv3", "id_lv2", "phenotypical_id", "fine_cell_type", "identity"],
            ["lv4", "lv3", "lv2", "phenotyp", "fine", "identity", "annotation"],
        ),
        "batch": choose_column(columns, ["batch", "library", "source", "orig.ident"], ["batch", "library", "source", "study"]),
    }


def compartment_flags(broad: object, fine: object) -> list[str]:
    text = f"{norm_text(broad)} {norm_text(fine)}".lower()
    flags: list[str] = []
    if any(token in text for token in ["epi", "tec", "ctec", "mtec", "aire"]):
        flags.append("epithelial_or_epi")
    if "mtec" in text or "post_aire" in text or "post-aire" in text or "aire" in text:
        flags.append("mtec_like")
    if "ctec" in text:
        flags.append("ctec_like")
    if "post_aire" in text or "post-aire" in text:
        flags.append("post_AIRE_mTEC")
    if "immature" in text and ("tec" in text or "epi" in text):
        flags.append("immature_TEC")
    if any(token in text for token in ["mes", "fibro", "fb", "vsmc", "smooth"]):
        flags.append("mesenchymal_or_fibroblast_like")
    if "endo" in text or "endothelial" in text:
        flags.append("endothelial")
    return flags or ["other"]


def load_metadata() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    workbook = download_if_missing(METADATA_URL, METADATA_XLSX)
    xls = pd.ExcelFile(workbook)
    normalized: list[pd.DataFrame] = []
    field_rows: list[dict[str, object]] = []
    for sheet in xls.sheet_names:
        df = pd.read_excel(workbook, sheet_name=sheet)
        df.columns = [str(column) for column in df.columns]
        if df.empty:
            continue
        inferred = infer_columns(df)
        for column in df.columns:
            lower = column.lower()
            roles = [role for role, selected in inferred.items() if selected == column]
            values = df[column].dropna().astype(str).unique()[:10]
            field_rows.append({
                "sheet": sheet,
                "field_name": column,
                "candidate_role": ";".join(roles) if roles else "other",
                "dtype": str(df[column].dtype),
                "non_null_count": int(df[column].notna().sum()),
                "unique_count": int(df[column].nunique(dropna=True)),
                "example_values": "; ".join(map(str, values)),
                "field_hint": lower,
            })
        if inferred["cell_id"]:
            cell = df[inferred["cell_id"]].map(norm_text)
        else:
            cell = pd.Series([f"{sheet}:{idx}" for idx in range(len(df))], index=df.index)
        sample_col = inferred["sample_id"] or inferred["donor_id"]
        donor_col = inferred["donor_id"] or inferred["sample_id"]
        norm = pd.DataFrame({
            "sheet": sheet,
            "row_index": np.arange(len(df)),
            "cell_id": cell,
            "donor_id": df[donor_col].map(norm_text) if donor_col else "",
            "sample_id": df[sample_col].map(norm_text) if sample_col else "",
            "age_raw": df[inferred["age"]].map(norm_text) if inferred["age"] else "",
            "sex": df[inferred["sex"]].map(norm_text) if inferred["sex"] else "",
            "broad_compartment": df[inferred["broad"]].map(norm_text) if inferred["broad"] else "",
            "fine_cell_type": df[inferred["fine"]].map(norm_text) if inferred["fine"] else "",
            "batch_source_library": df[inferred["batch"]].map(norm_text) if inferred["batch"] else "",
        })
        norm["age_years"] = norm["age_raw"].map(parse_age_years)
        norm["age_group"] = norm["age_years"].map(age_group)
        norm["target_compartments"] = [
            ";".join(compartment_flags(row.broad_compartment, row.fine_cell_type))
            for row in norm.itertuples(index=False)
        ]
        normalized.append(norm)
    if not normalized:
        raise RuntimeError(f"No usable metadata sheets found in {workbook}")
    return pd.concat(normalized, ignore_index=True), field_rows


def unique_join(values: Iterable[object], limit: int = 20) -> str:
    items = [norm_text(value) for value in values if norm_text(value)]
    return "; ".join(sorted(dict.fromkeys(items))[:limit])


def write_metadata_outputs(metadata: pd.DataFrame, field_rows: list[dict[str, object]]) -> None:
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_metadata_fields.tsv",
        field_rows,
        ["sheet", "field_name", "candidate_role", "dtype", "non_null_count", "unique_count", "example_values", "field_hint"],
    )
    counts = (
        metadata.groupby(["sheet", "broad_compartment", "fine_cell_type"], dropna=False)
        .agg(
            n_cells=("cell_id", "count"),
            n_unique_cell_ids=("cell_id", "nunique"),
            n_donors=("donor_id", "nunique"),
            n_samples=("sample_id", "nunique"),
            age_groups=("age_group", unique_join),
            compartments=("target_compartments", unique_join),
        )
        .reset_index()
    )
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_celltype_counts.tsv",
        counts.to_dict("records"),
        ["sheet", "broad_compartment", "fine_cell_type", "n_cells", "n_unique_cell_ids", "n_donors", "n_samples", "age_groups", "compartments"],
    )
    donor = (
        metadata.groupby(["donor_id", "sample_id"], dropna=False)
        .agg(
            n_cells=("cell_id", "count"),
            age_raw=("age_raw", unique_join),
            age_years_min=("age_years", "min"),
            age_years_max=("age_years", "max"),
            age_group=("age_group", unique_join),
            sex=("sex", unique_join),
            broad_compartments=("broad_compartment", unique_join),
            fine_cell_type_count=("fine_cell_type", "nunique"),
            target_compartments=("target_compartments", unique_join),
        )
        .reset_index()
    )
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_donor_age_summary.tsv",
        donor.to_dict("records"),
        ["donor_id", "sample_id", "n_cells", "age_raw", "age_years_min", "age_years_max", "age_group", "sex", "broad_compartments", "fine_cell_type_count", "target_compartments"],
    )
    exploded = metadata.assign(target_compartment=metadata["target_compartments"].str.split(";")).explode("target_compartment")
    inventory = (
        exploded.groupby(["target_compartment", "broad_compartment", "fine_cell_type"], dropna=False)
        .agg(
            n_cells=("cell_id", "count"),
            n_donors=("donor_id", "nunique"),
            n_samples=("sample_id", "nunique"),
            age_groups=("age_group", unique_join),
        )
        .reset_index()
        .sort_values(["target_compartment", "n_cells"], ascending=[True, False])
    )
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_target_compartment_inventory.tsv",
        inventory.to_dict("records"),
        ["target_compartment", "broad_compartment", "fine_cell_type", "n_cells", "n_donors", "n_samples", "age_groups"],
    )


def metadata_report(metadata: pd.DataFrame) -> None:
    target_labels = metadata.loc[metadata["target_compartments"].str.contains("epithelial|mtec|ctec|mesenchymal|endothelial", case=False, na=False)]
    lines = [
        "# GSE231906 Metadata Audit",
        "",
        "## Purpose",
        "",
        "Audit the GSE231906 cell-level metadata for donor/sample-aware human thymus aging context before expression parsing.",
        "",
        "## Local Data",
        "",
        f"- Metadata workbook: `{rel(METADATA_XLSX)}`",
        f"- Workbook parsed: {'yes' if METADATA_XLSX.exists() else 'no'}",
        f"- Cells represented in metadata: {len(metadata)}",
        f"- Donor/sample labels available: {'yes' if metadata['donor_id'].astype(bool).any() or metadata['sample_id'].astype(bool).any() else 'no'}",
        f"- Age field usable: {'yes' if metadata['age_years'].notna().any() else 'no'}",
        f"- Sex field usable: {'yes' if metadata['sex'].astype(bool).any() else 'no'}",
        f"- Target compartment rows identified: {len(target_labels)}",
        "",
        "## Target Compartments",
        "",
    ]
    for compartment, group in metadata.assign(target_compartment=metadata["target_compartments"].str.split(";")).explode("target_compartment").groupby("target_compartment"):
        lines.append(f"- {compartment}: {len(group)} cells, {group['donor_id'].nunique()} donors/samples by donor field, {group['sample_id'].nunique()} samples by sample field")
    lines.extend([
        "",
        "## Interpretation",
        "",
        "The metadata are technically usable for donor/sample-aware inventory work if the expression barcodes can be linked back to the workbook. Expression feasibility is not concluded by this metadata-only step.",
        "",
        REQUIRED_NO_OVERCLAIM,
    ])
    (REPORTS / "human_thymus_stage4_gse231906_metadata_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if not write_disk_report():
        return 0
    metadata, field_rows = load_metadata()
    write_metadata_outputs(metadata, field_rows)
    metadata_report(metadata)
    print(f"metadata_rows={len(metadata)}")
    print(f"metadata_workbook={METADATA_XLSX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
