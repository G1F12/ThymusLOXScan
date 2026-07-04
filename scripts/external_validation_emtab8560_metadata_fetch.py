"""Fetch official small E-MTAB-8560 metadata from BioStudies/ArrayExpress.

No raw sequencing files are downloaded.
"""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import requests


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "data" / "external" / "metadata" / "emtab8560"
REPORT = PROJECT_ROOT / "reports" / "external_emtab8560_metadata_source_notes.md"

BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8560"
SDRF_URL = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-8560/E-MTAB-8560.sdrf.txt"
IDF_URL = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-8560/E-MTAB-8560.idf.txt"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def fetch_text(url: str) -> str:
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    return response.text


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)

    sdrf_text = fetch_text(SDRF_URL)
    idf_text = fetch_text(IDF_URL)
    (OUT_DIR / "emtab8560_sdrf.tsv").write_text(sdrf_text, encoding="utf-8")
    (OUT_DIR / "emtab8560_idf.tsv").write_text(idf_text, encoding="utf-8")

    sdrf = pd.read_csv(io.StringIO(sdrf_text), sep="\t", dtype=str).fillna("")
    runs = sorted(x for x in sdrf.get("Comment[ENA_RUN]", pd.Series(dtype=str)).unique() if x)
    ena_written = False
    if runs:
        response = requests.get(
            ENA_URL,
            params={
                "accession": ",".join(runs),
                "result": "read_run",
                "fields": ",".join(
                    [
                        "run_accession",
                        "experiment_accession",
                        "sample_accession",
                        "secondary_sample_accession",
                        "study_accession",
                        "instrument_platform",
                        "instrument_model",
                        "library_strategy",
                        "library_source",
                        "library_selection",
                        "library_layout",
                        "center_name",
                        "first_public",
                        "last_updated",
                    ]
                ),
                "format": "tsv",
                "download": "false",
            },
            timeout=180,
        )
        if response.status_code == 200 and response.text.strip():
            (OUT_DIR / "emtab8560_ena_runs.tsv").write_text(response.text, encoding="utf-8")
            ena_written = True

    fields = list(sdrf.columns)
    field_lines = "\n".join(f"- `{field}`" for field in fields)
    unique_age = sorted(sdrf.get("Characteristics[age]", pd.Series(dtype=str)).unique())
    unique_individual = sorted(sdrf.get("Characteristics[individual]", pd.Series(dtype=str)).unique())
    unique_cell_type = sorted(sdrf.get("Characteristics[cell type]", pd.Series(dtype=str)).unique())

    notes = [
        "# E-MTAB-8560 Metadata Source Notes",
        "",
        "## Official URLs Used",
        "",
        f"- BioStudies/ArrayExpress study page: {BIOSTUDIES}",
        f"- SDRF metadata: {SDRF_URL}",
        f"- IDF metadata: {IDF_URL}",
        "- ENA portal API queried from `Comment[ENA_RUN]` accessions.",
        "",
        "## Metadata Files Fetched",
        "",
        "- `data/external/metadata/emtab8560/emtab8560_sdrf.tsv`",
        "- `data/external/metadata/emtab8560/emtab8560_idf.tsv`",
        "- `data/external/metadata/emtab8560/emtab8560_ena_runs.tsv`" if ena_written else "- ENA run metadata was queried but no table was written.",
        "",
        "## Key Fields Found",
        "",
        f"- Individual/mouse candidate field: `Characteristics[individual]` with values {', '.join(unique_individual[:20])}.",
        f"- Age field: `Characteristics[age]` / `Factor Value[age]` with values {', '.join(unique_age)}.",
        f"- Cell-type/sort fields include `Characteristics[cell type]`, `Characteristics[FACS marker]`, and `Characteristics[inferred cell type]`; observed cell types include {', '.join(unique_cell_type[:20])}.",
        "- Run/batch-related fields include `Comment[ENA_RUN]`, `Assay Name`, source names, and library metadata.",
        "- PlateID and SortDay are expected in MouseThymusAgeing colData rather than the raw SDRF field names.",
        "",
        "## Missing or Unresolved",
        "",
        "- The metadata fetch alone does not verify that exported R colData can map every analyzed cell to a biological individual.",
        "- Batch-aware modeling requires successful R export and a downstream biological-unit audit.",
        "",
        "## SDRF Fields",
        "",
        field_lines,
    ]
    REPORT.write_text("\n".join(notes) + "\n", encoding="utf-8")

    summary = pd.DataFrame(
        [
            {
                "field": "Characteristics[individual]",
                "n_unique": len(unique_individual),
                "values": ";".join(unique_individual),
            },
            {
                "field": "Characteristics[age]",
                "n_unique": len(unique_age),
                "values": ";".join(unique_age),
            },
            {
                "field": "Characteristics[cell type]",
                "n_unique": len(unique_cell_type),
                "values": ";".join(unique_cell_type),
            },
            {
                "field": "Comment[ENA_RUN]",
                "n_unique": len(runs),
                "values": ";".join(runs[:20]),
            },
        ]
    )
    summary.to_csv(OUT_DIR / "emtab8560_metadata_source_notes.tsv", sep="\t", index=False)
    print(f"Wrote SDRF rows: {len(sdrf)}")
    print(f"Wrote ENA metadata: {'yes' if ena_written else 'no'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
