#!/usr/bin/env python
"""
Download raw input datasets for ThymusLOXScan.

Dataset notes
-------------
Primary requested GEO accession: GSE198256
    Current GEO metadata: Homo sapiens monocyte RNA-seq from acute and
    convalescent severe COVID-19 patients.
    Important: this accession does not appear to be the mouse aging thymus
    scRNA-seq dataset described in the project prompt. The script still
    downloads the available GEO supplementary files for auditability, but users
    should verify the intended thymus accession before biological analysis.
    Species: Homo sapiens according to GEO.
    Age groups: not a young-vs-old mouse thymus design according to GEO.
    Cell numbers: not provided as a thymus scRNA-seq AnnData file by GEO.
    Sequencing method: RNA-seq, Illumina NovaSeq 6000 according to GEO.

Backup dataset URL: Human Cell Atlas / CELLxGENE thymus URL
    URL: https://datasets.cellxgene.cziscience.com/thymus
    Expected contents: human thymus single-cell data if the endpoint resolves
    to a downloadable AnnData object.
    Species: Homo sapiens.
    Age groups: human developmental, pediatric, and/or adult thymus depending
    on the specific CELLxGENE-hosted object served by the endpoint.
    Cell numbers: determined after download by reading the .h5ad file.
    Sequencing method: single-cell RNA-seq.
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"

GEO_ACCESSION = "GSE198256"
GEO_RECORD_URL = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GEO_ACCESSION}"
GEO_SUPPL_DIR_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/{GEO_ACCESSION}/suppl/"
HCA_THYMUS_URL = "https://datasets.cellxgene.cziscience.com/thymus"


def human_size(num_bytes: int) -> str:
    """Return a compact human-readable file size."""
    size = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{num_bytes} B"


def download_url(url: str, output_path: Path, timeout: int = 60) -> Path | None:
    """Download one URL with clear error messages."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading: {url}")
    print(f"Destination: {output_path}")

    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "ThymusLOXScan/0.1 (+https://www.ncbi.nlm.nih.gov/geo/)"},
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            with output_path.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
    except urllib.error.HTTPError as exc:
        print(f"ERROR: HTTP {exc.code} while downloading {url}: {exc.reason}", file=sys.stderr)
        return None
    except urllib.error.URLError as exc:
        print(f"ERROR: Network failure while downloading {url}: {exc.reason}", file=sys.stderr)
        return None
    except TimeoutError:
        print(f"ERROR: Timed out while downloading {url}", file=sys.stderr)
        return None
    except OSError as exc:
        print(f"ERROR: Could not write {output_path}: {exc}", file=sys.stderr)
        return None

    if output_path.exists():
        print(f"Downloaded {output_path.name} ({human_size(output_path.stat().st_size)})")
        return output_path

    print(f"ERROR: Download finished but {output_path} was not created.", file=sys.stderr)
    return None


def list_geo_supplementary_files() -> list[str]:
    """Return supplementary filenames available in the GEO FTP directory."""
    print(f"Checking GEO supplementary directory: {GEO_SUPPL_DIR_URL}")
    try:
        with urllib.request.urlopen(GEO_SUPPL_DIR_URL, timeout=60) as response:
            html = response.read().decode("utf-8", errors="replace")
    except Exception as exc:
        print(f"ERROR: Could not list GEO supplementary files: {exc}", file=sys.stderr)
        return []

    hrefs = re.findall(r'href="([^"]+)"', html)
    files = []
    for href in hrefs:
        name = urllib.parse.unquote(href)
        if name in {"../", "./"} or name.endswith("/"):
            continue
        files.append(name)
    return sorted(set(files))


def download_geo_accession() -> list[Path]:
    """Download all supplementary files for the requested GEO accession."""
    print("\n=== GEO primary requested dataset ===")
    print(f"Accession: {GEO_ACCESSION}")
    print(f"GEO record: {GEO_RECORD_URL}")
    print(
        "WARNING: GEO currently describes GSE198256 as human monocyte COVID-19 RNA-seq, "
        "not mouse aging thymus scRNA-seq."
    )

    downloaded = []
    for filename in list_geo_supplementary_files():
        url = GEO_SUPPL_DIR_URL + urllib.parse.quote(filename)
        path = download_url(url, RAW_DIR / filename)
        if path is not None:
            downloaded.append(path)

    if not downloaded:
        print("No GEO supplementary files were downloaded.", file=sys.stderr)

    return downloaded


def download_hca_backup() -> Path | None:
    """Attempt to download the backup CELLxGENE/HCA thymus AnnData file."""
    print("\n=== Backup Human Cell Atlas / CELLxGENE thymus dataset ===")
    print("Attempting backup download from the URL supplied in the project prompt.")
    return download_url(HCA_THYMUS_URL, RAW_DIR / "hca_thymus.h5ad")


def print_h5ad_info(path: Path) -> None:
    """Print basic AnnData information if a downloaded file is readable as .h5ad."""
    try:
        import scanpy as sc
    except ImportError:
        print("scanpy is not installed; skipping .h5ad inspection.")
        return

    try:
        adata = sc.read_h5ad(path)
    except Exception as exc:
        print(f"{path.name}: not readable as .h5ad ({exc})")
        return

    print(f"{path.name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"obs columns: {list(adata.obs.columns)}")


def print_download_summary(paths: list[Path]) -> None:
    """Print file sizes and lightweight metadata after downloading."""
    print("\n=== Download summary ===")
    if not paths:
        print("No files were downloaded.")
        return

    for path in paths:
        size = human_size(path.stat().st_size)
        print(f"{path.relative_to(PROJECT_ROOT)} - {size}")
        if path.suffix.lower() == ".h5ad":
            print_h5ad_info(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download raw data for ThymusLOXScan.")
    parser.add_argument(
        "--include-backup",
        action="store_true",
        help="Always attempt the Human Cell Atlas / CELLxGENE thymus backup download.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    downloaded = download_geo_accession()
    has_h5ad = any(path.suffix.lower() == ".h5ad" for path in downloaded)

    if args.include_backup or not has_h5ad:
        backup = download_hca_backup()
        if backup is not None:
            downloaded.append(backup)

    print_download_summary(downloaded)

    if not any(path.suffix.lower() == ".h5ad" for path in downloaded):
        print(
            "\nWARNING: No .h5ad file was downloaded. The inspection notebook expects an .h5ad "
            "file in data/raw/. Verify the correct aging thymus accession or provide a direct "
            "CELLxGENE .h5ad download URL.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
