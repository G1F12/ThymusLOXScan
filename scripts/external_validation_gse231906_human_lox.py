#!/usr/bin/env python
"""Guarded human thymus LOX-family validation scaffold for GSE231906.

The GSE231906 raw expression archive is large, so this script downloads/parses
metadata by default and refuses to download the full matrix automatically unless
explicitly allowed. If a compatible local matrix is already available, the
script can be extended from the discovered metadata paths without treating cells
as independent biological replicates.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

try:
    import anndata as ad
except ImportError:  # pragma: no cover - optional dependency in minimal environments
    ad = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data" / "external" / "GSE231906"
DEFAULT_EXISTING_METADATA = PROJECT_ROOT / "tmp" / "external_metadata" / "GSE231906_cell-level_metadata.xlsx"
DEFAULT_SUMMARY = PROJECT_ROOT / "results" / "tables" / "external_gse231906_human_lox_summary.tsv"
DEFAULT_BY_DONOR = PROJECT_ROOT / "results" / "tables" / "external_gse231906_human_lox_by_donor.tsv"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "external_validation" / "gse231906"
DEFAULT_REPORT = PROJECT_ROOT / "reports" / "gse231906_human_external_validation.md"
DEFAULT_METADATA_PLAN = PROJECT_ROOT / "reports" / "gse231906_human_metadata_only_plan.md"

GEO_RAW_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE231906&format=file"
GEO_METADATA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/"
    "GSE231906_cell-level_metadata.xlsx"
)
RAW_ARCHIVE_NAME = "GSE231906_RAW.tar"
RAW_ARCHIVE_EXPECTED_BYTES = 3_961_333_760
METADATA_EXPECTED_BYTES = 18_065_742
DEFAULT_MAX_AUTO_DOWNLOAD_BYTES = 2_000_000_000

LOX_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]
CLAIMS = [
    ("fibroblast_like_mesenchymal", "LOX", "current fibroblast/capsFB LOX aged-lower direction"),
    ("fibroblast_like_mesenchymal", "LOXL1", "current medFB LOXL1/LOXl1 aged-higher direction"),
    ("fibroblast_like_mesenchymal", "LOXL2", "current medFB LOXL2/Loxl2 aged-lower direction"),
    ("mTEC_like", "LOXL2", "current mTEC1 LOXL2/Loxl2 aged-lower direction"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--existing-metadata", type=Path, default=DEFAULT_EXISTING_METADATA)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--by-donor-output", type=Path, default=DEFAULT_BY_DONOR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--metadata-plan", type=Path, default=DEFAULT_METADATA_PLAN)
    parser.add_argument("--max-auto-download-bytes", type=int, default=DEFAULT_MAX_AUTO_DOWNLOAD_BYTES)
    parser.add_argument("--allow-large-download", action="store_true")
    return parser.parse_args()


def bytes_to_gb(value: int) -> float:
    return value / 1024**3


def disk_free_bytes(path: Path) -> int:
    path.mkdir(parents=True, exist_ok=True)
    return shutil.disk_usage(path).free


def ensure_metadata(data_dir: Path, existing_metadata: Path) -> Path:
    data_dir.mkdir(parents=True, exist_ok=True)
    dest = data_dir / "GSE231906_cell-level_metadata.xlsx"
    if dest.exists():
        return dest
    if existing_metadata.exists():
        shutil.copy2(existing_metadata, dest)
        return dest
    print(f"Downloading metadata workbook: {GEO_METADATA_URL}", flush=True)
    urllib.request.urlretrieve(GEO_METADATA_URL, dest)
    return dest


def parse_age_years(value) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    match = re.search(r"(\d+(?:\.\d+)?)", text)
    if not match:
        return np.nan
    age = float(match.group(1))
    if "month" in text:
        return age / 12
    if "week" in text:
        return age / 52
    if "day" in text:
        return age / 365.25
    return age


def age_bin(age: float) -> str:
    if pd.isna(age):
        return "age_unknown"
    if age < 1:
        return "infant_lt1y"
    if age < 12:
        return "child_1_11y"
    if age < 18:
        return "adolescent_12_17y"
    if age < 40:
        return "adult_18_39y"
    if age < 65:
        return "adult_40_64y"
    return "older_65plus"


def classify_group(row: pd.Series) -> list[str]:
    labels = " ".join(str(row.get(col, "")) for col in row.index).lower()
    groups = []
    if any(token in labels for token in ["fb_1", "fb_2", "fibro", "mes", "vsmc"]):
        groups.append("fibroblast_like_mesenchymal")
    if any(token in labels for token in ["epi", "tec", "mtec", "ctec", "post_aire"]):
        groups.append("epithelial")
    if any(token in labels for token in ["mtec", "post_aire"]):
        groups.append("mTEC_like")
    if "ctec" in labels:
        groups.append("cTEC_like")
    if any(token in labels for token in ["endo", "endothelial"]):
        groups.append("endothelial")
    return groups or ["other"]


def load_metadata(workbook: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    xl = pd.ExcelFile(workbook)
    all_rows = []
    group_rows = []
    for sheet in xl.sheet_names:
        if sheet == "Table Caption":
            continue
        df = pd.read_excel(workbook, sheet_name=sheet)
        df.columns = [str(col) for col in df.columns]
        donor_col = "geo_sample_id" if "geo_sample_id" in df.columns else "sample_id" if "sample_id" in df.columns else None
        age_col = "Age" if "Age" in df.columns else None
        sex_col = "Sex" if "Sex" in df.columns else None
        source_col = "source" if "source" in df.columns else None
        broad_col = "id_lv1" if "id_lv1" in df.columns else "identity" if "identity" in df.columns else None
        fine_candidates = [col for col in ["id_lv4", "id_lv3", "id_lv2", "identity", "phenotypical_id"] if col in df.columns]
        fine_col = fine_candidates[0] if fine_candidates else None

        normalized = pd.DataFrame(
            {
                "sheet": sheet,
                "barcode": df["barcode"].astype(str) if "barcode" in df.columns else pd.Series(pd.NA, index=df.index),
                "donor_id": df[donor_col].astype(str) if donor_col else pd.Series(pd.NA, index=df.index),
                "age": df[age_col] if age_col else pd.Series(pd.NA, index=df.index),
                "sex": df[sex_col] if sex_col else pd.Series(pd.NA, index=df.index),
                "source_study": df[source_col] if source_col else pd.Series(pd.NA, index=df.index),
                "broad_cell_type": df[broad_col] if broad_col else pd.Series(pd.NA, index=df.index),
                "fine_cell_type": df[fine_col] if fine_col else pd.Series(pd.NA, index=df.index),
            }
        )
        normalized["age_years"] = normalized["age"].map(parse_age_years)
        normalized["age_bin"] = normalized["age_years"].map(age_bin)
        normalized["candidate_groups"] = [
            ";".join(classify_group(row)) for _, row in normalized.iterrows()
        ]
        all_rows.append(normalized)

        exploded = normalized.assign(candidate_group=normalized["candidate_groups"].str.split(";")).explode(
            "candidate_group"
        )
        group_rows.append(exploded)

    return pd.concat(all_rows, ignore_index=True), pd.concat(group_rows, ignore_index=True)


def build_metadata_tables(metadata: pd.DataFrame, grouped: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_donor = (
        grouped.groupby(["candidate_group", "sheet", "donor_id"], dropna=False)
        .agg(
            n_cells=("barcode", "nunique"),
            age_years=("age_years", "first"),
            age=("age", "first"),
            age_bin=("age_bin", "first"),
            sex=("sex", "first"),
            source_study=("source_study", "first"),
            broad_cell_types=("broad_cell_type", lambda x: ";".join(sorted({str(v) for v in x.dropna()}))),
            fine_cell_types=("fine_cell_type", lambda x: ";".join(sorted({str(v) for v in x.dropna()}))[:1000]),
        )
        .reset_index()
    )
    for gene in LOX_GENES:
        by_donor[f"{gene}_pseudobulk_expression"] = pd.NA
        by_donor[f"{gene}_detection_rate"] = pd.NA
    by_donor["expression_status"] = "not_available_metadata_only"

    summary = (
        grouped.groupby(["candidate_group", "sheet"], dropna=False)
        .agg(
            n_cells=("barcode", "nunique"),
            n_donors_or_samples=("donor_id", "nunique"),
            n_donors_with_age=("donor_id", lambda x: int(grouped.loc[x.index].dropna(subset=["age_years"])["donor_id"].nunique())),
            min_age_years=("age_years", "min"),
            max_age_years=("age_years", "max"),
            sex_values=("sex", lambda x: ";".join(sorted({str(v) for v in x.dropna()}))),
            source_values=("source_study", lambda x: ";".join(sorted({str(v) for v in x.dropna()}))),
            broad_cell_types=("broad_cell_type", lambda x: ";".join(sorted({str(v) for v in x.dropna()}))[:1000]),
            fine_cell_type_examples=("fine_cell_type", lambda x: ";".join(sorted({str(v) for v in x.dropna()}))[:1000]),
        )
        .reset_index()
    )
    summary["human_lox_genes_targeted"] = ",".join(LOX_GENES)
    summary["expression_status"] = "not_available_metadata_only"
    summary["analysis_note"] = (
        "Metadata parsed; expression matrix not analyzed. Do not treat cells as biological replicates."
    )
    return summary, by_donor


def is_gzip(path: Path) -> bool:
    return path.suffix == ".gz" or path.name.endswith(".gz")


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if is_gzip(path) else path.open("rt", encoding="utf-8", errors="replace")


def detect_delimiter(path: Path) -> str:
    with open_text(path) as handle:
        line = handle.readline()
    return "\t" if line.count("\t") >= line.count(",") else ","


def local_expression_status(data_dir: Path) -> dict[str, object]:
    raw_path = data_dir / RAW_ARCHIVE_NAME
    extracted_candidates = (
        list(data_dir.glob("**/*.h5ad"))
        + list(data_dir.glob("**/matrix.mtx*"))
        + list(data_dir.glob("**/*matrix*.csv*"))
        + list(data_dir.glob("**/*matrix*.tsv*"))
        + list(data_dir.glob("**/*counts*.csv*"))
        + list(data_dir.glob("**/*counts*.tsv*"))
        + list(data_dir.glob("**/*expression*.csv*"))
        + list(data_dir.glob("**/*expression*.tsv*"))
    )
    status: dict[str, object] = {
        "raw_archive_present": raw_path.exists(),
        "raw_archive_path": str(raw_path) if raw_path.exists() else "",
        "raw_archive_size_bytes": raw_path.stat().st_size if raw_path.exists() else 0,
        "extracted_matrix_candidates": [str(path) for path in extracted_candidates[:20]],
    }
    if raw_path.exists():
        try:
            with tarfile.open(raw_path) as tar:
                status["raw_archive_member_examples"] = tar.getnames()[:25]
        except tarfile.TarError as exc:
            status["raw_archive_member_examples"] = [f"Could not inspect tar: {exc}"]
    else:
        status["raw_archive_member_examples"] = []
    return status


def find_feature_barcode_files(matrix_path: Path) -> tuple[Path | None, Path | None]:
    folder = matrix_path.parent
    feature_patterns = ["features.tsv*", "genes.tsv*", "*features*.tsv*", "*genes*.tsv*"]
    barcode_patterns = ["barcodes.tsv*", "*barcodes*.tsv*"]
    features = next((p for pat in feature_patterns for p in folder.glob(pat)), None)
    barcodes = next((p for pat in barcode_patterns for p in folder.glob(pat)), None)
    return features, barcodes


def read_single_column(path: Path) -> list[str]:
    values = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            values.append(line.rstrip("\n").split("\t")[0])
    return values


def expression_candidates(data_dir: Path) -> list[Path]:
    status = local_expression_status(data_dir)
    return [Path(path) for path in status["extracted_matrix_candidates"]]


def pseudobulk_from_dense_table(matrix_path: Path, donor_map: pd.DataFrame) -> pd.DataFrame | None:
    sep = detect_delimiter(matrix_path)
    header = pd.read_csv(matrix_path, sep=sep, nrows=0)
    columns = [str(col) for col in header.columns]
    lox_cols = [col for col in columns if col.upper() in LOX_GENES]

    if lox_cols:
        usecols = [columns[0], *lox_cols]
        df = pd.read_csv(matrix_path, sep=sep, usecols=usecols)
        barcode_col = columns[0]
        df = df.rename(columns={barcode_col: "barcode"})
    else:
        first_col = columns[0]
        gene_hits = []
        for chunk in pd.read_csv(matrix_path, sep=sep, chunksize=20000):
            gene_col = first_col
            subset = chunk.loc[chunk[gene_col].astype(str).str.upper().isin(LOX_GENES)]
            if not subset.empty:
                gene_hits.append(subset)
        if not gene_hits:
            return None
        gene_df = pd.concat(gene_hits, ignore_index=True)
        gene_col = first_col
        value_cols = [col for col in gene_df.columns if col != gene_col]
        df = gene_df.set_index(gene_col)[value_cols].T.reset_index().rename(columns={"index": "barcode"})

    df["barcode"] = df["barcode"].astype(str)
    merged = donor_map.merge(df, on="barcode", how="inner")
    if merged.empty:
        return None
    rows = []
    for keys, group in merged.groupby(["candidate_group", "sheet", "donor_id"], dropna=False):
        row = {
            "candidate_group": keys[0],
            "sheet": keys[1],
            "donor_id": keys[2],
            "n_cells_with_expression": int(group["barcode"].nunique()),
        }
        for gene in LOX_GENES:
            col = next((c for c in group.columns if c.upper() == gene), None)
            if col is None:
                row[f"{gene}_pseudobulk_expression"] = pd.NA
                row[f"{gene}_detection_rate"] = pd.NA
            else:
                values = pd.to_numeric(group[col], errors="coerce").fillna(0)
                row[f"{gene}_pseudobulk_expression"] = float(values.sum())
                row[f"{gene}_detection_rate"] = float((values > 0).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def pseudobulk_from_mtx(matrix_path: Path, donor_map: pd.DataFrame) -> pd.DataFrame | None:
    features_path, barcodes_path = find_feature_barcode_files(matrix_path)
    if not features_path or not barcodes_path:
        return None
    features = read_single_column(features_path)
    barcodes = read_single_column(barcodes_path)
    gene_indices = [i for i, gene in enumerate(features) if str(gene).upper() in LOX_GENES]
    if not gene_indices:
        return None
    mat = mmread(str(matrix_path)).tocsr()
    if mat.shape[0] == len(features) and mat.shape[1] == len(barcodes):
        sub = mat[gene_indices, :].T
    elif mat.shape[1] == len(features) and mat.shape[0] == len(barcodes):
        sub = mat[:, gene_indices]
    else:
        return None
    expr = pd.DataFrame.sparse.from_spmatrix(sub, index=barcodes, columns=[features[i].upper() for i in gene_indices])
    expr = expr.reset_index().rename(columns={"index": "barcode"})
    return pseudobulk_from_expression_frame(expr, donor_map)


def pseudobulk_from_expression_frame(expr: pd.DataFrame, donor_map: pd.DataFrame) -> pd.DataFrame | None:
    expr["barcode"] = expr["barcode"].astype(str)
    merged = donor_map.merge(expr, on="barcode", how="inner")
    if merged.empty:
        return None
    rows = []
    for keys, group in merged.groupby(["candidate_group", "sheet", "donor_id"], dropna=False):
        row = {
            "candidate_group": keys[0],
            "sheet": keys[1],
            "donor_id": keys[2],
            "n_cells_with_expression": int(group["barcode"].nunique()),
        }
        for gene in LOX_GENES:
            if gene not in group.columns:
                row[f"{gene}_pseudobulk_expression"] = pd.NA
                row[f"{gene}_detection_rate"] = pd.NA
            else:
                values = pd.to_numeric(group[gene], errors="coerce").fillna(0)
                row[f"{gene}_pseudobulk_expression"] = float(values.sum())
                row[f"{gene}_detection_rate"] = float((values > 0).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def pseudobulk_from_h5ad(matrix_path: Path, donor_map: pd.DataFrame) -> pd.DataFrame | None:
    if ad is None:
        return None
    adata = ad.read_h5ad(matrix_path, backed="r")
    var_names = pd.Index([str(v).upper() for v in adata.var_names])
    gene_indices = [i for i, gene in enumerate(var_names) if gene in LOX_GENES]
    if not gene_indices:
        adata.file.close()
        return None
    barcodes = pd.Index([str(v) for v in adata.obs_names])
    matched = donor_map.loc[donor_map["barcode"].astype(str).isin(barcodes)].copy()
    if matched.empty:
        adata.file.close()
        return None
    barcode_to_pos = pd.Series(np.arange(len(barcodes)), index=barcodes)
    positions = matched["barcode"].astype(str).map(barcode_to_pos).dropna().astype(int).to_numpy()
    sub = adata.X[positions, :][:, gene_indices]
    if sparse.issparse(sub):
        sub = sub.toarray()
    expr = pd.DataFrame(sub, columns=[var_names[i] for i in gene_indices])
    expr.insert(0, "barcode", matched["barcode"].astype(str).to_numpy())
    adata.file.close()
    return pseudobulk_from_expression_frame(expr, donor_map)


def parse_local_expression(data_dir: Path, grouped: pd.DataFrame) -> tuple[pd.DataFrame | None, str]:
    donor_map = grouped[
        ["candidate_group", "sheet", "donor_id", "barcode", "age_years", "age", "age_bin", "sex", "source_study"]
    ].dropna(subset=["barcode", "donor_id"])
    donor_map["barcode"] = donor_map["barcode"].astype(str)
    for candidate in expression_candidates(data_dir):
        try:
            if candidate.suffix == ".h5ad":
                parsed = pseudobulk_from_h5ad(candidate, donor_map)
            elif candidate.name.startswith("matrix.mtx"):
                parsed = pseudobulk_from_mtx(candidate, donor_map)
            else:
                parsed = pseudobulk_from_dense_table(candidate, donor_map)
        except Exception as exc:  # keep the guard from failing on unknown huge layouts
            parsed = None
            last_error = f"{candidate}: {exc}"
        else:
            last_error = f"{candidate}: incompatible layout or no matching barcodes/LOX genes"
        if parsed is not None and not parsed.empty:
            return parsed, str(candidate)
    return None, locals().get("last_error", "no local expression matrix candidate found")


def merge_expression_into_by_donor(by_donor: pd.DataFrame, expression: pd.DataFrame | None) -> pd.DataFrame:
    if expression is None:
        return by_donor
    base = by_donor.drop(columns=[col for col in by_donor.columns if col.endswith("_pseudobulk_expression") or col.endswith("_detection_rate")], errors="ignore")
    merged = base.merge(expression, on=["candidate_group", "sheet", "donor_id"], how="left")
    for gene in LOX_GENES:
        for suffix in ["pseudobulk_expression", "detection_rate"]:
            col = f"{gene}_{suffix}"
            if col not in merged.columns:
                merged[col] = pd.NA
    merged["expression_status"] = np.where(
        merged["n_cells_with_expression"].fillna(0).astype(float) > 0,
        "donor_level_pseudobulk_available",
        "not_available_for_donor_group",
    )
    return merged


def update_summary_expression_status(summary: pd.DataFrame, by_donor: pd.DataFrame, expression_available: bool) -> pd.DataFrame:
    summary = summary.copy()
    if not expression_available:
        return summary
    expr_counts = (
        by_donor.loc[by_donor["expression_status"].eq("donor_level_pseudobulk_available")]
        .groupby(["candidate_group", "sheet"], dropna=False)
        .agg(n_donors_with_expression=("donor_id", "nunique"), n_cells_with_expression=("n_cells_with_expression", "sum"))
        .reset_index()
    )
    summary = summary.merge(expr_counts, on=["candidate_group", "sheet"], how="left")
    summary["expression_status"] = np.where(
        summary["n_donors_with_expression"].fillna(0).astype(float) > 0,
        "donor_level_pseudobulk_available",
        "not_available_for_group",
    )
    summary["analysis_note"] = (
        "Local expression matrix parsed only to donor/sample-level LOX-family pseudobulk summaries; cells are not treated as biological replicates."
    )
    return summary


def claim_classification(expression_available: bool) -> pd.DataFrame:
    rows = []
    for group, gene, claim in CLAIMS:
        rows.append(
            {
                "claim_group": group,
                "gene": gene,
                "current_claim": claim,
                "classification": "inconclusive" if expression_available else "not testable",
                "reason": (
                    "A local expression matrix was parsed to donor-level summaries, but no formal human age model is claimed by this guarded script."
                    if expression_available
                    else "Only metadata was parsed; LOX-family expression was not available for donor-level pseudobulk analysis."
                ),
            }
        )
    return pd.DataFrame(rows)


def fmt_number(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def write_reports(
    report: Path,
    metadata_plan: Path,
    summary: pd.DataFrame,
    by_donor: pd.DataFrame,
    expression_status: dict[str, object],
    free_bytes: int,
    expression_available: bool,
    expression_source: str,
) -> None:
    report.parent.mkdir(parents=True, exist_ok=True)
    metadata_plan.parent.mkdir(parents=True, exist_ok=True)
    claims = claim_classification(expression_available)
    group_lines = []
    for _, row in summary.sort_values(["candidate_group", "sheet"]).iterrows():
        group_lines.append(
            f"| {row['candidate_group']} | {row['sheet']} | {int(row['n_cells'])} | "
            f"{int(row['n_donors_or_samples'])} | {fmt_number(row['min_age_years'])} | "
            f"{fmt_number(row['max_age_years'])} |"
        )

    lines = [
        "# GSE231906 human external validation",
        "",
        "## Guard status",
        "",
        f"- GEO raw archive: `{RAW_ARCHIVE_NAME}`, expected size {bytes_to_gb(RAW_ARCHIVE_EXPECTED_BYTES):.2f} GB.",
        f"- Available disk at data directory: {bytes_to_gb(free_bytes):.2f} GB.",
        f"- Metadata workbook expected size: {bytes_to_gb(METADATA_EXPECTED_BYTES):.3f} GB.",
        f"- Raw archive present locally: {expression_status['raw_archive_present']}.",
        f"- Extracted matrix candidates found: {len(expression_status['extracted_matrix_candidates'])}.",
        f"- Expression parser status: {expression_source}.",
        "",
        "## Metadata-derived candidate groups",
        "",
        "| candidate group | sheet | cells | donors/samples | min age years | max age years |",
        "|---|---|---:|---:|---:|---:|",
        *group_lines,
        "",
        "## Claim Classification",
        "",
        "| group | gene | classification | reason |",
        "|---|---|---|---|",
    ]
    for _, row in claims.iterrows():
        lines.append(f"| {row['claim_group']} | {row['gene']} | {row['classification']} | {row['reason']} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                "A compatible local expression matrix was parsed into donor/sample-level LOX-family pseudobulk summaries. This guarded script still does not fit a formal human age model or support a human conservation claim."
                if expression_available
                else "No human LOX-family expression analysis was performed. This run is metadata-only and must not be used to claim human conservation of LOX-family aging patterns. GSE231906 is only a candidate future validation dataset. The metadata identifies promising fibroblast-like/mesenchymal, epithelial, mTEC-like, cTEC-like, and endothelial compartments, but expression is required for validation."
            ),
            "",
            "If expression is later parsed, the analysis should use donor-level pseudobulk summaries, model age as a continuous variable, include source/study covariates where possible, and avoid treating cells as independent biological replicates.",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")

    manual = [
        "# GSE231906 human metadata-only plan",
        "",
        "## Why this stopped before full validation",
        "",
        f"The full GEO raw archive is {bytes_to_gb(RAW_ARCHIVE_EXPECTED_BYTES):.2f} GB. The guarded pipeline does not download it automatically unless `--allow-large-download` is used and disk space is adequate. This avoids silently pulling a large TAR of CSV/MTX/TSV files.",
        "",
        "## Manual download instructions",
        "",
        f"1. Download `{RAW_ARCHIVE_NAME}` from `{GEO_RAW_URL}`.",
        f"2. Save it to `{(DEFAULT_DATA_DIR / RAW_ARCHIVE_NAME).as_posix()}`.",
        "3. Inspect the TAR members before extraction, for example with Python `tarfile` or `tar -tf GSE231906_RAW.tar`.",
        "4. Extract only the expression matrix, feature/gene, and barcode files needed for thymic stromal/epithelial cells.",
        "5. Match barcodes to the metadata workbook by donor/sample/source columns.",
        "6. Build donor-level pseudobulk summaries for LOX, LOXL1, LOXL2, LOXL3, and LOXL4 in fibroblast-like/mesenchymal, epithelial, mTEC-like, cTEC-like, and endothelial groups.",
        "7. Analyze LOX-family expression versus age as donor-level continuous-age trends, with source/study covariates when possible.",
        "",
        "## Parsed metadata outputs",
        "",
        f"- `{DEFAULT_SUMMARY.relative_to(PROJECT_ROOT).as_posix()}`",
        f"- `{DEFAULT_BY_DONOR.relative_to(PROJECT_ROOT).as_posix()}`",
        "",
        "## Current claim status",
        "",
        (
            "A local expression matrix was detected and summarized, but no human conservation claim is supported without a prespecified donor-level age model."
            if expression_available
            else "No human LOX-family expression analysis was performed. No human conservation claim is supported. GSE231906 is only a candidate future validation dataset."
        ),
    ]
    metadata_plan.write_text("\n".join(manual), encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    free_bytes = disk_free_bytes(args.data_dir)
    metadata_path = ensure_metadata(args.data_dir, args.existing_metadata)
    metadata, grouped = load_metadata(metadata_path)
    summary, by_donor = build_metadata_tables(metadata, grouped)
    expression_status = local_expression_status(args.data_dir)

    raw_path = args.data_dir / RAW_ARCHIVE_NAME
    if (
        not raw_path.exists()
        and RAW_ARCHIVE_EXPECTED_BYTES <= args.max_auto_download_bytes
        and free_bytes > RAW_ARCHIVE_EXPECTED_BYTES * 2
        and args.allow_large_download
    ):
        urllib.request.urlretrieve(GEO_RAW_URL, raw_path)
        expression_status = local_expression_status(args.data_dir)

    expression_table, expression_source = parse_local_expression(args.data_dir, grouped)
    expression_available = expression_table is not None
    by_donor = merge_expression_into_by_donor(by_donor, expression_table)
    summary = update_summary_expression_status(summary, by_donor, expression_available)

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)
    by_donor.to_csv(args.by_donor_output, sep="\t", index=False)
    write_reports(
        args.report,
        args.metadata_plan,
        summary,
        by_donor,
        expression_status,
        free_bytes,
        expression_available,
        expression_source,
    )

    print(f"Saved summary table: {args.summary_output}", flush=True)
    print(f"Saved donor metadata table: {args.by_donor_output}", flush=True)
    print(f"Saved report: {args.report}", flush=True)
    print(f"Saved metadata-only plan: {args.metadata_plan}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
