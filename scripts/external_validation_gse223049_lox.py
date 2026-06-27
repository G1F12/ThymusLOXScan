#!/usr/bin/env python
"""Download and summarize LOX-family expression in external GSE223049 thymus samples."""

from __future__ import annotations

import argparse
import gzip
import shutil
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data" / "external" / "GSE223049"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "tables" / "external_gse223049_lox_validation.tsv"

COUNTS_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223049/suppl/GSE223049_RNA_seq_counts_23_cell_types.txt.gz"
LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
TARGET_PREFIXES = ["Thymic_epithelial", "Thymic_fibroblasts"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def download_counts(data_dir: Path) -> Path:
    data_dir.mkdir(parents=True, exist_ok=True)
    gz_path = data_dir / "GSE223049_RNA_seq_counts_23_cell_types.txt.gz"
    txt_path = data_dir / "GSE223049_RNA_seq_counts_23_cell_types.txt"
    if not gz_path.exists() and not txt_path.exists():
        print(f"Downloading {COUNTS_URL}", flush=True)
        urllib.request.urlretrieve(COUNTS_URL, gz_path)
    if not txt_path.exists():
        print(f"Decompressing {gz_path}", flush=True)
        with gzip.open(gz_path, "rb") as src, txt_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    return txt_path


def parse_sample(sample: str) -> tuple[str | None, str | None, str | None]:
    for prefix in TARGET_PREFIXES:
        if sample.startswith(prefix + "_"):
            suffix = sample.removeprefix(prefix + "_")
            parts = suffix.split("_")
            if len(parts) >= 2 and parts[0] in {"Y", "O"}:
                age = "young_2mo" if parts[0] == "Y" else "aged_22_24mo"
                replicate = parts[1]
                return prefix, age, replicate
    return None, None, None


def main() -> int:
    args = parse_args()
    counts_path = download_counts(args.data_dir)
    counts = pd.read_csv(counts_path, sep="\t")
    gene_col = counts.columns[0]
    counts[gene_col] = counts[gene_col].astype(str)
    counts = counts.set_index(gene_col)

    target_samples = []
    sample_meta = {}
    for sample in counts.columns:
        cell_type, age, replicate = parse_sample(sample)
        if cell_type:
            target_samples.append(sample)
            sample_meta[sample] = {
                "external_cell_type": cell_type,
                "age_group": age,
                "replicate": replicate,
            }

    rows = []
    lib_sizes = counts[target_samples].sum(axis=0)
    for gene in LOX_GENES:
        if gene not in counts.index:
            for sample in target_samples:
                rows.append({"gene": gene, "sample_id": sample, "present": False, **sample_meta[sample]})
            continue
        for sample in target_samples:
            raw = float(counts.loc[gene, sample])
            cpm = raw / float(lib_sizes[sample]) * 1_000_000 if lib_sizes[sample] else np.nan
            rows.append(
                {
                    "gene": gene,
                    "sample_id": sample,
                    "present": True,
                    **sample_meta[sample],
                    "raw_counts": raw,
                    "library_size": float(lib_sizes[sample]),
                    "cpm": cpm,
                    "log2_cpm_plus1": np.log2(cpm + 1) if np.isfinite(cpm) else np.nan,
                }
            )

    sample_df = pd.DataFrame(rows)
    summary = (
        sample_df.groupby(["external_cell_type", "gene", "age_group"], observed=True)
        .agg(
            n_samples=("sample_id", "nunique"),
            mean_log2_cpm_plus1=("log2_cpm_plus1", "mean"),
            mean_cpm=("cpm", "mean"),
            mean_raw_counts=("raw_counts", "mean"),
        )
        .reset_index()
    )
    wide = summary.pivot_table(
        index=["external_cell_type", "gene"],
        columns="age_group",
        values="mean_log2_cpm_plus1",
        aggfunc="first",
    ).reset_index()
    wide["delta_log2_cpm_aged_minus_young"] = wide.get("aged_22_24mo", np.nan) - wide.get("young_2mo", np.nan)

    out = sample_df.merge(
        wide[["external_cell_type", "gene", "delta_log2_cpm_aged_minus_young"]],
        on=["external_cell_type", "gene"],
        how="left",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary_path = args.output.with_name(args.output.stem + "_summary.tsv")
    summary.merge(wide, on=["external_cell_type", "gene"], how="left").to_csv(summary_path, sep="\t", index=False)
    print(f"Saved sample-level table: {args.output}", flush=True)
    print(f"Saved summary table: {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
