"""Fetch small official GEO/SRA metadata tables for GSE240016.

This script intentionally downloads metadata only. It does not download FASTQ,
BAM, SRA, or count matrix files.
"""

from __future__ import annotations

import csv
import gzip
import io
import re
import sys
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "data" / "external" / "metadata"

GEO_SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/"
    "GSE240016/soft/GSE240016_family.soft.gz"
)
GEO_MINIML_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/"
    "GSE240016/miniml/GSE240016_family.xml.tgz"
)
EUTILS_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EUTILS_FETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def fetch_url(url: str) -> bytes:
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.content


def clean_text(value: str | None) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def parse_characteristics(lines: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in lines:
        if ":" in item:
            key, value = item.split(":", 1)
            values[key.strip().lower().replace(" ", "_")] = clean_text(value)
    return values


def parse_geo_soft(text: str) -> pd.DataFrame:
    rows = []
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        accession = clean_text(lines[0])
        multi: dict[str, list[str]] = {
            "characteristics": [],
            "relations": [],
            "supplementary_files": [],
        }
        row = {"gsm_accession": accession}

        for line in lines[1:]:
            if not line.startswith("!Sample_"):
                continue
            key, value = line.split("=", 1)
            key = key.strip().replace("!Sample_", "").replace("_ch1", "")
            value = clean_text(value)
            if key == "characteristics":
                multi["characteristics"].append(value)
            elif key == "relation":
                multi["relations"].append(value)
            elif key == "supplementary_file":
                multi["supplementary_files"].append(value)
            else:
                row[key] = value

        characteristics = parse_characteristics(multi["characteristics"])
        for key, value in characteristics.items():
            row[f"characteristics_{key}"] = value
        row["all_characteristics"] = "; ".join(multi["characteristics"])
        row["relations"] = "; ".join(multi["relations"])
        row["supplementary_files"] = "; ".join(multi["supplementary_files"])
        srx = sorted(set(re.findall(r"SRX\d+", row["relations"])))
        row["srx_accession"] = ";".join(srx)
        rows.append(row)

    return pd.DataFrame(rows).sort_values("gsm_accession")


def fetch_sra_runinfo(srx_accessions: list[str]) -> pd.DataFrame:
    if not srx_accessions:
        return pd.DataFrame()

    search = requests.get(
        EUTILS_SEARCH,
        params={
            "db": "sra",
            "term": " OR ".join(srx_accessions),
            "retmode": "json",
            "retmax": 1000,
        },
        timeout=60,
    )
    search.raise_for_status()
    ids = search.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return pd.DataFrame()

    fetched = requests.get(
        EUTILS_FETCH,
        params={
            "db": "sra",
            "id": ",".join(ids),
            "rettype": "runinfo",
            "retmode": "text",
        },
        timeout=120,
    )
    fetched.raise_for_status()
    df = pd.read_csv(io.StringIO(fetched.text))

    # Keep metadata only. Exclude raw-data URL and hashes.
    drop_cols = [c for c in ["download_path", "RunHash", "ReadHash"] if c in df.columns]
    df = df.drop(columns=drop_cols)
    return df.sort_values(["LibraryName", "Run"], na_position="last")


def fetch_ena_run_metadata(run_accessions: list[str]) -> pd.DataFrame:
    if not run_accessions:
        return pd.DataFrame()
    fields = [
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "secondary_sample_accession",
        "study_accession",
        "secondary_study_accession",
        "sample_title",
        "experiment_title",
        "first_public",
        "last_updated",
        "instrument_platform",
        "instrument_model",
        "library_name",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "center_name",
        "fastq_bytes",
    ]
    response = requests.get(
        ENA_PORTAL,
        params={
            "accession": ",".join(run_accessions),
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "tsv",
            "download": "false",
        },
        timeout=120,
    )
    if response.status_code != 200 or not response.text.strip():
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO(response.text), sep="\t")
    # Keep size metadata but not file URLs or raw data paths.
    return df.sort_values("run_accession")


def fetch_biosample_metadata(biosample_accessions: list[str]) -> pd.DataFrame:
    rows = []
    for biosample in sorted(set(x for x in biosample_accessions if x)):
        response = requests.get(
            EUTILS_FETCH,
            params={"db": "biosample", "id": biosample, "retmode": "xml"},
            timeout=60,
        )
        if response.status_code != 200:
            rows.append({"BioSample": biosample, "fetch_status": f"http_{response.status_code}"})
            continue
        try:
            root = ET.fromstring(response.content)
        except ET.ParseError:
            rows.append({"BioSample": biosample, "fetch_status": "parse_error"})
            continue
        sample = root.find(".//BioSample")
        row = {"BioSample": biosample, "fetch_status": "ok"}
        if sample is not None:
            row["accession"] = sample.attrib.get("accession", "")
            row["publication_date"] = sample.attrib.get("publication_date", "")
            row["last_update"] = sample.attrib.get("last_update", "")
            for attr in sample.findall(".//Attribute"):
                key = attr.attrib.get("attribute_name") or attr.attrib.get("harmonized_name") or "attribute"
                key = key.lower().replace(" ", "_")
                row[f"biosample_{key}"] = clean_text(attr.text)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("BioSample")


def inspect_miniml_available() -> dict[str, str]:
    """Record whether MINiML is accessible without saving the XML payload."""
    try:
        payload = fetch_url(GEO_MINIML_URL)
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
            members = [m.name for m in archive.getmembers()]
        return {"miniml_url": GEO_MINIML_URL, "status": "ok", "members": ";".join(members)}
    except Exception as exc:  # noqa: BLE001 - diagnostic metadata script
        return {"miniml_url": GEO_MINIML_URL, "status": f"failed: {exc}", "members": ""}


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    soft_payload = fetch_url(GEO_SOFT_URL)
    soft_text = gzip.decompress(soft_payload).decode("utf-8", errors="replace")
    geo_samples = parse_geo_soft(soft_text)
    geo_samples.to_csv(OUT_DIR / "gse240016_geo_samples.tsv", sep="\t", index=False)

    srx_accessions = sorted(
        set(
            item
            for values in geo_samples["srx_accession"].dropna().astype(str)
            for item in values.split(";")
            if item
        )
    )
    runinfo = fetch_sra_runinfo(srx_accessions)
    runinfo.to_csv(OUT_DIR / "gse240016_sra_runinfo.tsv", sep="\t", index=False)

    run_accessions = runinfo["Run"].dropna().astype(str).tolist() if "Run" in runinfo else []
    ena = fetch_ena_run_metadata(run_accessions)
    if not ena.empty:
        ena.to_csv(OUT_DIR / "gse240016_ena_run_metadata.tsv", sep="\t", index=False)

    biosamples = runinfo["BioSample"].dropna().astype(str).tolist() if "BioSample" in runinfo else []
    biosample_df = fetch_biosample_metadata(biosamples)
    biosample_df.to_csv(OUT_DIR / "gse240016_biosample_metadata.tsv", sep="\t", index=False)

    source_notes = pd.DataFrame([inspect_miniml_available()])
    source_notes.to_csv(OUT_DIR / "gse240016_metadata_fetch_notes.tsv", sep="\t", index=False)

    print(f"Wrote {len(geo_samples)} GEO sample rows")
    print(f"Wrote {len(runinfo)} SRA run rows")
    print(f"Wrote {len(biosample_df)} BioSample rows")
    if not ena.empty:
        print(f"Wrote {len(ena)} ENA run rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
