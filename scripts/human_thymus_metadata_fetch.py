"""Fetch public metadata inventories for human thymus feasibility audit.

This script intentionally records metadata and file inventories only. It does
not download expression matrices into the repository and does not fetch large
archives. The only optional file download is the small GSE231906 metadata
workbook to a temporary directory for column inspection.
"""

from __future__ import annotations

import csv
import html.parser
import json
import os
import re
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "tables"
STAGE1_TABLE = RESULTS / "human_thymus_dataset_search_candidates.tsv"

CELLXGENE_API = "https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}"
GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
GEO_SUPPL = {
    "GSE147520": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147520/suppl/",
    "GSE231906": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/",
}
HCA_PARK = "https://explore.data.humancellatlas.org/projects/c1810dbc-16d2-45c3-b45e-3e675f88d87b"
ARRAYEXPRESS_PARK = "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8581"

TOP_CANDIDATES = {
    "HTD_CELLXGENE_PARK_TEC",
    "HTD_CELLXGENE_YAYON_TEC",
    "HTD_GSE147520_STROMA",
    "HTD_GSE231906_AGING",
}

CELLXGENE_COLLECTIONS = {
    "HTD_CELLXGENE_PARK_TEC": {
        "accession": "CELLxGENE de13e3e2-23b6-40ed-a413-e9e12d7d3910",
        "collection_id": "de13e3e2-23b6-40ed-a413-e9e12d7d3910",
        "preferred_title_pattern": "Human Thymic Epithelial Cells",
    },
    "HTD_CELLXGENE_YAYON_TEC": {
        "accession": "CELLxGENE fc19ae6c-d7c1-4dce-b703-62c5d52061b4",
        "collection_id": "fc19ae6c-d7c1-4dce-b703-62c5d52061b4",
        "preferred_title_pattern": "thymic epithelial cell subset",
    },
}


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[tuple[str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        attr = dict(attrs)
        href = attr.get("href")
        if href:
            self.links.append((href, ""))


def fetch_bytes(url: str, limit: int = 20_000_000) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "ThymusLOXScan-metadata-audit"})
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                content_length = response.headers.get("Content-Length")
                if content_length and int(content_length) > limit:
                    raise RuntimeError(f"Refusing to fetch >{limit} bytes from {url}")
                data = response.read(limit + 1)
            break
        except Exception as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(2)
    else:
        raise RuntimeError(f"Failed to fetch {url}: {last_error}")
    if len(data) > limit:
        raise RuntimeError(f"Refusing to keep >{limit} bytes from {url}")
    return data


def fetch_text(url: str, limit: int = 20_000_000) -> str:
    data = fetch_bytes(url, limit=limit)
    return data.decode("utf-8", errors="replace")


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "ThymusLOXScan-metadata-audit"})
    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            length = response.headers.get("Content-Length")
            return int(length) if length else None
    except Exception:
        return None


def human_size(size: int | None) -> str:
    if size is None:
        return "unknown"
    value = float(size)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024 or unit == "TB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return f"{size} B"


def classify_file(name: str, filetype: str = "") -> dict[str, str]:
    fmt = filetype or Path(urllib.parse.urlparse(name).path).suffix.lstrip(".").upper() or "unknown"
    lower = f"{name} {fmt}".lower()
    is_raw = any(token in lower for token in ["raw", "fastq", ".bam", ".sra"])
    is_metadata = any(token in lower for token in ["metadata", "sample", "feature_reference", "soft", "miniml"])
    is_matrix = fmt.upper() in {"H5AD", "H5", "RDS", "MTX"} or any(
        token in lower for token in [".h5ad", ".mtx", ".rds", ".h5", "matrix", "counts", "abundance"]
    )
    is_processed = is_matrix or is_metadata
    return {
        "file_format": fmt,
        "is_expression_matrix_candidate": "yes" if is_matrix and not is_metadata else "no",
        "is_metadata_candidate": "yes" if is_metadata else "no",
        "is_raw_archive": "yes" if is_raw else "no",
        "is_processed": "yes" if is_processed and not is_raw else "no",
    }


def download_recommendation(size: int | None, is_raw: str, is_metadata: str, is_matrix: str) -> tuple[str, str]:
    if is_raw == "yes":
        return "no", "raw archive or raw-like file is outside Stage 2 scope"
    if is_metadata == "yes":
        return "yes", "metadata file or endpoint is small enough for audit"
    if size is not None and size > 100_000_000:
        return "no", "expression candidate is above 100 MB; record only, do not download in Stage 2"
    if is_matrix == "yes":
        return "no", "expression candidate should not be downloaded into repository during feasibility audit"
    return "yes", "small non-expression metadata or landing page"


def read_stage1_candidates() -> dict[str, dict[str, str]]:
    if not STAGE1_TABLE.exists():
        raise FileNotFoundError(f"Missing Stage 1 table: {STAGE1_TABLE}")
    with STAGE1_TABLE.open(newline="", encoding="utf-8") as handle:
        return {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def add_inventory(rows: list[dict[str, str]], candidate_id: str, accession: str, source: str,
                  file_name: str, file_url: str, file_type: str, size: int | None,
                  reason_override: str | None = None) -> None:
    flags = classify_file(file_name, file_type)
    if candidate_id == "HTD_GSE231906_AGING" and file_name == "GSE231906_RAW.tar":
        flags["is_expression_matrix_candidate"] = "yes"
        flags["is_metadata_candidate"] = "no"
        flags["is_raw_archive"] = "no"
        flags["is_processed"] = "yes"
    recommend, reason = download_recommendation(
        size,
        flags["is_raw_archive"],
        flags["is_metadata_candidate"],
        flags["is_expression_matrix_candidate"],
    )
    if reason_override:
        reason = reason_override
    rows.append({
        "candidate_id": candidate_id,
        "accession": accession,
        "source": source,
        "file_name": file_name,
        "file_url": file_url,
        "file_type": file_type or "unknown",
        "file_format": flags["file_format"],
        "estimated_size_bytes": "" if size is None else str(size),
        "estimated_size_human": human_size(size),
        "is_expression_matrix_candidate": flags["is_expression_matrix_candidate"],
        "is_metadata_candidate": flags["is_metadata_candidate"],
        "is_raw_archive": flags["is_raw_archive"],
        "is_processed": flags["is_processed"],
        "download_recommended_now": recommend,
        "reason": reason,
    })


def add_access(rows: list[dict[str, str]], candidate_id: str, accession: str, source: str,
               url_type: str, url: str, manual: str, notes: str) -> None:
    rows.append({
        "candidate_id": candidate_id,
        "accession": accession,
        "source": source,
        "url_type": url_type,
        "url": url,
        "requires_manual_download": manual,
        "notes": notes,
    })


def add_field(rows: list[dict[str, str]], candidate_id: str, accession: str, source: str,
              endpoint: str, field: str, example: str, meaning: str,
              needed: str, notes: str) -> None:
    rows.append({
        "candidate_id": candidate_id,
        "accession": accession,
        "source": source,
        "metadata_file_or_endpoint": endpoint,
        "field_name": field,
        "field_example": example,
        "likely_meaning": meaning,
        "needed_for_reanalysis": needed,
        "notes": notes,
    })


def cellxgene_metadata(inventory: list[dict[str, str]], fields: list[dict[str, str]],
                       access: list[dict[str, str]]) -> None:
    for candidate_id, info in CELLXGENE_COLLECTIONS.items():
        accession = info["accession"]
        api_url = CELLXGENE_API.format(collection_id=info["collection_id"])
        collection = json.loads(fetch_text(api_url, limit=50_000_000))
        add_access(access, candidate_id, accession, "CELLxGENE", "collection_api", api_url, "no",
                   "Public metadata endpoint used for file and field inventory")
        add_access(access, candidate_id, accession, "CELLxGENE", "collection_page", collection["collection_url"], "no",
                   "Interactive collection page")

        pattern = info["preferred_title_pattern"].lower()
        preferred_seen = False
        for dataset in collection.get("datasets", []):
            title = dataset.get("title", "")
            is_preferred = pattern in title.lower()
            if is_preferred:
                preferred_seen = True
            if is_preferred or candidate_id == "HTD_CELLXGENE_YAYON_TEC" and any(
                token in title.lower() for token in ["fibroblast subset", "vascular cell subset"]
            ):
                for asset in dataset.get("assets", []):
                    add_inventory(
                        inventory,
                        candidate_id,
                        accession,
                        "CELLxGENE",
                        title,
                        asset.get("url", ""),
                        asset.get("filetype", "H5AD"),
                        int(asset["filesize"]) if asset.get("filesize") is not None else None,
                    )
                add_access(access, candidate_id, accession, "CELLxGENE", "dataset_explorer",
                           dataset.get("explorer_url", ""), "no", f"Explorer for {title}")

                field_specs = [
                    ("dataset_id", dataset.get("dataset_id", ""), "CELLxGENE dataset identifier", "yes"),
                    ("dataset_version_id", dataset.get("dataset_version_id", ""), "versioned dataset asset identifier", "yes"),
                    ("title", title, "dataset title", "yes"),
                    ("donor_id", "; ".join(dataset.get("donor_id", [])[:5]), "donor identifier list", "yes"),
                    ("development_stage", "; ".join(x["label"] for x in dataset.get("development_stage", [])[:5]), "age/development stage", "yes"),
                    ("sex", "; ".join(x["label"] for x in dataset.get("sex", [])[:5]), "donor sex", "yes"),
                    ("cell_type", "; ".join(x["label"] for x in dataset.get("cell_type", [])[:8]), "curated cell-type labels", "yes"),
                    ("assay", "; ".join(x["label"] for x in dataset.get("assay", [])[:5]), "library/assay type", "yes"),
                    ("raw_data_location", str(dataset.get("raw_data_location", "")), "expression layer pointer", "yes"),
                    ("suspension_type", "; ".join(dataset.get("suspension_type", [])), "cell/nucleus/spatial suspension type", "no"),
                ]
                for field_name, example, meaning, needed in field_specs:
                    add_field(fields, candidate_id, accession, "CELLxGENE", api_url, field_name,
                              example, meaning, needed, f"Metadata for {title}")

        if not preferred_seen:
            add_field(fields, candidate_id, accession, "CELLxGENE", api_url, "preferred_dataset",
                      "not found", "preferred Stage 2 dataset title match", "yes",
                      "Collection API was reachable but expected title pattern was not found")

    add_access(access, "HTD_CELLXGENE_PARK_TEC", "HCA c1810dbc-16d2-45c3-b45e-3e675f88d87b",
               "Human Cell Atlas", "project_page", HCA_PARK, "no",
               "HCA project page for Park atlas metadata and accessions")
    add_access(access, "HTD_CELLXGENE_PARK_TEC", "E-MTAB-8581", "ArrayExpress/BioStudies",
               "study_page", ARRAYEXPRESS_PARK, "no",
               "ArrayExpress/BioStudies landing page; CELLxGENE metadata is more directly useful for Stage 2")


def parse_geo_supplement_listing(url: str) -> list[tuple[str, str, int | None]]:
    text = fetch_text(url, limit=5_000_000)
    parser = LinkParser()
    parser.feed(text)
    files: list[tuple[str, str, int | None]] = []
    for href, _ in parser.links:
        name = urllib.parse.unquote(href.split("/")[-1])
        if not name or name in {"../", "/"} or name.startswith("?"):
            continue
        if name.endswith("/"):
            continue
        file_url = urllib.parse.urljoin(url, href)
        if not file_url.startswith(url):
            continue
        size = head_size(file_url)
        files.append((name, file_url, size))
    return files


def extract_geo_fields(accession: str, text: str) -> list[tuple[str, str, str, str]]:
    fields: list[tuple[str, str, str, str]] = []
    patterns = {
        "title": r"Title\s+([^\n]+)",
        "organism": r"Organism\s+([^\n]+)",
        "experiment_type": r"Experiment type\s+([^\n]+)",
        "overall_design": r"Overall design\s+([^\n]+)",
        "samples": r"Samples \((\d+)\)",
        "bioproject": r"BioProject\s+([A-Z0-9]+)",
    }
    for field, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            fields.append((field, match.group(1).strip(), "GEO series metadata", "yes" if field in {"title", "samples", "overall_design"} else "no"))
    for sample in re.findall(r"(GSM\d+)\s+([^\n]+)", text)[:8]:
        fields.append(("sample_title", f"{sample[0]} {sample[1].strip()}", "sample identifier and title", "yes"))
    return fields


def gse231906_metadata_workbook(fields: list[dict[str, str]], inventory: list[dict[str, str]]) -> None:
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231906/suppl/GSE231906%5Fcell%2Dlevel%5Fmetadata%2Exlsx"
    size = head_size(url)
    add_inventory(
        inventory,
        "HTD_GSE231906_AGING",
        "GSE231906",
        "GEO",
        "GSE231906_cell-level_metadata.xlsx",
        url,
        "XLSX",
        size,
        "metadata workbook is below 100 MB and may be fetched to a temp directory for schema inspection",
    )
    if size is not None and size > 100_000_000:
        return
    try:
        import pandas as pd  # type: ignore
    except Exception as exc:
        add_field(fields, "HTD_GSE231906_AGING", "GSE231906", "GEO", url, "workbook_inspection",
                  "skipped", "metadata workbook field inspection", "yes", f"pandas unavailable: {exc}")
        return
    try:
        temp_path = Path(tempfile.gettempdir()) / "GSE231906_cell-level_metadata_stage2.xlsx"
        if not temp_path.exists() or (size and temp_path.stat().st_size != size):
            temp_path.write_bytes(fetch_bytes(url, limit=100_000_000))
        xls = pd.ExcelFile(temp_path)
        add_field(fields, "HTD_GSE231906_AGING", "GSE231906", "GEO", url, "sheets",
                  "; ".join(xls.sheet_names), "workbook sheet inventory", "yes",
                  "Downloaded to OS temp directory only, not repository")
        for sheet in xls.sheet_names:
            if sheet == "Table Caption":
                continue
            df = pd.read_excel(temp_path, sheet_name=sheet, nrows=20)
            for col in df.columns:
                examples = [str(x) for x in df[col].dropna().astype(str).unique()[:3]]
                needed = "yes" if col.lower() in {
                    "barcode", "geo_sample_id", "sample_id", "sex", "age", "id_lv1",
                    "id_lv2", "id_lv3", "id_lv4", "identity", "source"
                } else "no"
                meaning = {
                    "barcode": "cell barcode",
                    "geo_sample_id": "donor/sample identifier",
                    "sample_id": "sample identifier",
                    "Sex": "donor sex",
                    "Age": "donor age",
                    "id_lv1": "broad cell class",
                    "id_lv2": "intermediate cell class",
                    "id_lv3": "fine cell annotation",
                    "id_lv4": "fine cell annotation",
                    "identity": "cell identity annotation",
                    "source": "source dataset label",
                }.get(col, "metadata field")
                add_field(fields, "HTD_GSE231906_AGING", "GSE231906", "GEO", f"{url}#{sheet}",
                          col, "; ".join(examples), meaning, needed, "Field observed in metadata workbook")
        summary_cols = ["geo_sample_id", "Sex", "Age", "id_lv1", "id_lv2", "id_lv3", "id_lv4"]
        thymocyte = pd.read_excel(temp_path, sheet_name="thymocyte_metadata", usecols=summary_cols)
        for col in summary_cols:
            examples = [str(x) for x in thymocyte[col].dropna().astype(str).unique()[:20]]
            add_field(fields, "HTD_GSE231906_AGING", "GSE231906", "GEO",
                      f"{url}#thymocyte_metadata", f"{col}_unique_examples",
                      "; ".join(examples), f"unique examples for {col}", "yes",
                      "Whole-sheet column scan of small metadata workbook in OS temp")
        epi = thymocyte[thymocyte["id_lv1"].astype(str).str.contains("Epi", case=False, na=False)]
        for col in ["id_lv2", "id_lv3", "id_lv4"]:
            examples = [str(x) for x in epi[col].dropna().astype(str).unique()[:20]]
            add_field(fields, "HTD_GSE231906_AGING", "GSE231906", "GEO",
                      f"{url}#thymocyte_metadata", f"epithelial_{col}_examples",
                      "; ".join(examples), "epithelial compartment annotation examples", "yes",
                      "Epithelial rows detected from id_lv1 in metadata workbook")
        stromal = pd.read_excel(temp_path, sheet_name="stromal_cell_metadata")
        for col in ["identity", "source", "sample_id"]:
            if col in stromal.columns:
                examples = [str(x) for x in stromal[col].dropna().astype(str).unique()[:20]]
                add_field(fields, "HTD_GSE231906_AGING", "GSE231906", "GEO",
                          f"{url}#stromal_cell_metadata", f"stromal_{col}_examples",
                          "; ".join(examples), f"stromal metadata examples for {col}", "yes",
                          "Whole-sheet column scan of small metadata workbook in OS temp")
    except Exception as exc:
        add_field(fields, "HTD_GSE231906_AGING", "GSE231906", "GEO", url, "workbook_inspection",
                  "failed", "metadata workbook field inspection", "yes", str(exc))


def geo_metadata(inventory: list[dict[str, str]], fields: list[dict[str, str]],
                 access: list[dict[str, str]]) -> None:
    geo_candidates = {
        "HTD_GSE147520_STROMA": "GSE147520",
        "HTD_GSE231906_AGING": "GSE231906",
    }
    for candidate_id, accession in geo_candidates.items():
        geo_url = GEO_BASE.format(accession=accession)
        add_access(access, candidate_id, accession, "GEO", "series_page", geo_url, "no",
                   "Official GEO series metadata page")
        try:
            text = fetch_text(geo_url, limit=5_000_000)
            for field, example, meaning, needed in extract_geo_fields(accession, text):
                add_field(fields, candidate_id, accession, "GEO", geo_url, field, example, meaning, needed,
                          "Parsed from GEO series page")
        except Exception as exc:
            add_field(fields, candidate_id, accession, "GEO", geo_url, "series_page_fetch",
                      "failed", "GEO series page reachability", "yes", str(exc))
        suppl_url = GEO_SUPPL[accession]
        add_access(access, candidate_id, accession, "GEO", "supplementary_directory", suppl_url, "no",
                   "FTP/HTTPS supplementary-file directory; inventory only")
        try:
            for name, file_url, size in parse_geo_supplement_listing(suppl_url):
                add_inventory(inventory, candidate_id, accession, "GEO", name, file_url,
                              Path(name).suffix.lstrip(".").upper() or "unknown", size)
        except Exception as exc:
            add_field(fields, candidate_id, accession, "GEO", suppl_url, "supplementary_directory_fetch",
                      "failed", "GEO supplementary directory reachability", "yes", str(exc))
    gse231906_metadata_workbook(fields, inventory)


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def main() -> None:
    candidates = read_stage1_candidates()
    missing = sorted(TOP_CANDIDATES - set(candidates))
    if missing:
        raise RuntimeError(f"Stage 1 candidate table missing expected candidates: {', '.join(missing)}")

    inventory: list[dict[str, str]] = []
    fields: list[dict[str, str]] = []
    access: list[dict[str, str]] = []

    cellxgene_metadata(inventory, fields, access)
    geo_metadata(inventory, fields, access)

    write_tsv(RESULTS / "human_thymus_dataset_file_inventory.tsv", inventory, [
        "candidate_id", "accession", "source", "file_name", "file_url", "file_type",
        "file_format", "estimated_size_bytes", "estimated_size_human",
        "is_expression_matrix_candidate", "is_metadata_candidate", "is_raw_archive",
        "is_processed", "download_recommended_now", "reason",
    ])
    write_tsv(RESULTS / "human_thymus_metadata_fields.tsv", fields, [
        "candidate_id", "accession", "source", "metadata_file_or_endpoint", "field_name",
        "field_example", "likely_meaning", "needed_for_reanalysis", "notes",
    ])
    write_tsv(RESULTS / "human_thymus_candidate_access_urls.tsv", access, [
        "candidate_id", "accession", "source", "url_type", "url",
        "requires_manual_download", "notes",
    ])


if __name__ == "__main__":
    main()
