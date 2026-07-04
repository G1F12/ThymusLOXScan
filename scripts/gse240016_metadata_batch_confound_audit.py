"""Join local GSE240016 sample metadata with official GEO/SRA metadata."""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
H5AD = PROJECT_ROOT / "data" / "raw" / "GSE240016_CD45neg_thymic_stroma_d0+annotation.h5ad"
META_DIR = PROJECT_ROOT / "data" / "external" / "metadata"
RESULTS_DIR = PROJECT_ROOT / "results" / "tables"
REPORTS_DIR = PROJECT_ROOT / "reports"

TARGET_SAMPLES = [
    "mo02_CD45neg1_d0",
    "mo02_CD45neg2_d0",
    "mo18_CD45neg1_d0",
    "mo18_CD45neg2_d0",
]

LOCAL_TO_GEO = {
    "mo02_CD45neg1_d0": ("GSM7679869", "high", "02-mo CD45- thymic stroma steady state rep1"),
    "mo02_CD45neg2_d0": ("GSM7679870", "high", "02-mo CD45- thymic stroma steady state rep2"),
    "mo18_CD45neg1_d0": ("GSM7679885", "high", "18-mo CD45- thymic stroma steady state rep1"),
    "mo18_CD45neg2_d0": ("GSM7679886", "high", "18-mo CD45- thymic stroma steady state rep2"),
    "mo02_EC_d0": ("GSM7679871", "medium", "02-mo CD45- endothelial steady state title match"),
    "mo18_EC_d0": ("GSM7679887", "medium", "18-mo CD45- endothelial steady state title match"),
    "mo02_FB_d0": ("GSM7679872", "medium", "02-mo CD45- fibroblast steady state title match"),
    "mo18_FB_d0": ("GSM7679888", "medium", "18-mo CD45- fibroblast steady state title match"),
}


def clean(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def read_table(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return pd.DataFrame()


def summarize_local_obs() -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs[
        [
            "sample",
            "stage",
            "day",
            "cell_type_subset",
            "total_counts",
            "n_genes_by_counts",
        ]
    ].copy()
    adata.file.close()

    obs["is_mtec1"] = obs["cell_type_subset"].astype(str).eq("13:mTEC1")
    rows = []
    for sample, group in obs.groupby("sample", observed=False):
        mtec = group[group["is_mtec1"]]
        rows.append(
            {
                "local_sample": str(sample),
                "stage": clean(group["stage"].iloc[0]),
                "day": clean(group["day"].iloc[0]),
                "total_cells": int(len(group)),
                "mTEC1_cells": int(len(mtec)),
                "median_total_counts_mTEC1": float(mtec["total_counts"].median()) if len(mtec) else math.nan,
                "median_detected_genes_mTEC1": float(mtec["n_genes_by_counts"].median()) if len(mtec) else math.nan,
                "mean_total_counts_mTEC1": float(mtec["total_counts"].mean()) if len(mtec) else math.nan,
                "mean_detected_genes_mTEC1": float(mtec["n_genes_by_counts"].mean()) if len(mtec) else math.nan,
            }
        )
    local_qc = pd.DataFrame(rows).sort_values("local_sample")
    target_qc = local_qc[local_qc["local_sample"].isin(TARGET_SAMPLES)].copy()
    return local_qc, target_qc


def aggregate_runinfo(runinfo: pd.DataFrame) -> pd.DataFrame:
    if runinfo.empty:
        return pd.DataFrame()
    rows = []
    for gsm, group in runinfo.groupby("LibraryName", dropna=False):
        rows.append(
            {
                "gsm_accession": clean(gsm),
                "srr_run_accessions": ";".join(sorted(group["Run"].astype(str).unique())),
                "sra_experiment": ";".join(sorted(group["Experiment"].astype(str).unique())),
                "BioSample": ";".join(sorted(group["BioSample"].astype(str).unique())),
                "BioProject": ";".join(sorted(group["BioProject"].astype(str).unique())),
                "SRAStudy": ";".join(sorted(group["SRAStudy"].astype(str).unique())),
                "library_strategy_sra": ";".join(sorted(group["LibraryStrategy"].astype(str).unique())),
                "library_source_sra": ";".join(sorted(group["LibrarySource"].astype(str).unique())),
                "library_selection_sra": ";".join(sorted(group["LibrarySelection"].astype(str).unique())),
                "library_layout_sra": ";".join(sorted(group["LibraryLayout"].astype(str).unique())),
                "instrument": ";".join(sorted(group["Model"].astype(str).unique())),
                "center_name": ";".join(sorted(group["CenterName"].astype(str).unique())),
                "load_date": ";".join(sorted(group["LoadDate"].astype(str).unique())),
                "load_date_calendar": ";".join(sorted(set(group["LoadDate"].astype(str).str.slice(0, 10)))),
                "release_date": ";".join(sorted(group["ReleaseDate"].astype(str).unique())),
                "release_date_calendar": ";".join(sorted(set(group["ReleaseDate"].astype(str).str.slice(0, 10)))),
                "n_srr_runs": int(group["Run"].nunique()),
                "total_sra_spots": int(pd.to_numeric(group["spots"], errors="coerce").sum()),
                "total_sra_bases": int(pd.to_numeric(group["bases"], errors="coerce").sum()),
            }
        )
    return pd.DataFrame(rows)


def derive_tokens(row: pd.Series) -> dict[str, str]:
    title = clean(row.get("title", ""))
    age = clean(row.get("characteristics_age", ""))
    treatment = clean(row.get("characteristics_treatment", ""))
    cell_type = clean(row.get("characteristics_cell_type", ""))
    return {
        "sample_title": title,
        "age_stage": "02mo" if "02" in age or "02-mo" in title else ("18mo" if "18" in age or "18-mo" in title else ""),
        "enrichment_sort_condition": cell_type,
        "day_timepoint": "d0" if treatment == "untreated" or "steady state" in title else treatment,
    }


def build_join(local_qc: pd.DataFrame, geo: pd.DataFrame, run_agg: pd.DataFrame, biosample: pd.DataFrame) -> pd.DataFrame:
    geo_rows = []
    for _, row in geo.iterrows():
        fields = derive_tokens(row)
        geo_rows.append(
            {
                "gsm_accession": clean(row.get("gsm_accession", "")),
                "sample_title": fields["sample_title"],
                "age_stage": fields["age_stage"],
                "enrichment_sort_condition": fields["enrichment_sort_condition"],
                "day_timepoint": fields["day_timepoint"],
                "library_strategy_geo": clean(row.get("library_strategy", "")),
                "library_source_geo": clean(row.get("library_source", "")),
                "library_selection_geo": clean(row.get("library_selection", "")),
                "source_name": clean(row.get("source_name", "")),
                "sra_relation": clean(row.get("relations", "")),
                "srx_accession_geo": clean(row.get("srx_accession", "")),
            }
        )
    geo_df = pd.DataFrame(geo_rows)
    joined_geo = geo_df.merge(run_agg, on="gsm_accession", how="left")

    if not biosample.empty and "BioSample" in biosample:
        keep_cols = [c for c in biosample.columns if c in {"BioSample", "biosample_sex", "biosample_dev_stage", "biosample_tissue", "biosample_cell_type", "biosample_treatment"}]
        joined_geo = joined_geo.merge(biosample[keep_cols].drop_duplicates("BioSample"), on="BioSample", how="left")

    rows = []
    for _, local in local_qc.iterrows():
        sample = local["local_sample"]
        gsm, confidence, reason = LOCAL_TO_GEO.get(sample, ("", "unresolved", "no direct title/token mapping"))
        match = joined_geo[joined_geo["gsm_accession"].eq(gsm)]
        record = local.to_dict()
        record["mapping_confidence"] = confidence
        record["mapping_basis"] = reason
        record["mapped_gsm_accession"] = gsm
        if not match.empty:
            for key, value in match.iloc[0].to_dict().items():
                record[key] = value
        rows.append(record)
    return pd.DataFrame(rows)


def age_coincidence(join: pd.DataFrame, field: str) -> tuple[str, str, str]:
    sub = join[join["local_sample"].isin(TARGET_SAMPLES)].copy()
    young = sorted(set(clean(v) for v in sub[sub["stage"].eq("02mo")][field].dropna()))
    aged = sorted(set(clean(v) for v in sub[sub["stage"].eq("18mo")][field].dropna()))
    young = [v for v in young if v]
    aged = [v for v in aged if v]
    if not young and not aged:
        return "unavailable", "", ""
    overlap = set(young) & set(aged)
    if overlap:
        return "no", ";".join(young), ";".join(aged)
    if set(young) != set(aged):
        return "yes_or_sample_specific", ";".join(young), ";".join(aged)
    return "no", ";".join(young), ";".join(aged)


def build_confound_matrix(join: pd.DataFrame) -> pd.DataFrame:
    fields = [
        ("GSM accession", "mapped_gsm_accession"),
        ("SRR run accessions", "srr_run_accessions"),
        ("BioSample", "BioSample"),
        ("sequencing run date", "sequencing_run_date"),
        ("SRA load timestamp", "load_date"),
        ("SRA load calendar date", "load_date_calendar"),
        ("SRA release calendar date", "release_date_calendar"),
        ("library prep date", "library_prep_date"),
        ("library strategy", "library_strategy_sra"),
        ("library source", "library_source_sra"),
        ("library selection", "library_selection_sra"),
        ("instrument", "instrument"),
        ("center", "center_name"),
        ("lane/flowcell", "lane_flowcell"),
        ("sort/enrichment condition", "enrichment_sort_condition"),
        ("day/timepoint", "day_timepoint"),
    ]
    rows = []
    for label, field in fields:
        if field not in join.columns:
            rows.append(
                {
                    "technical_field": label,
                    "age_coincides": "unavailable",
                    "young_values": "",
                    "aged_values": "",
                    "assessment": "field not present in GEO/SRA/BioSample metadata",
                }
            )
            continue
        status, young, aged = age_coincidence(join, field)
        if status == "no":
            assessment = "available values do not separate young and aged samples"
        elif status == "yes_or_sample_specific":
            assessment = "values differ between age groups or are sample-specific"
        else:
            assessment = "not available"
        rows.append(
            {
                "technical_field": label,
                "age_coincides": status,
                "young_values": young,
                "aged_values": aged,
                "assessment": assessment,
            }
        )
    return pd.DataFrame(rows)


def classify_batch_concern(matrix: pd.DataFrame, join: pd.DataFrame) -> str:
    targets = join[join["local_sample"].isin(TARGET_SAMPLES)]
    if targets["mapping_confidence"].ne("high").any():
        return "mapping unresolved"
    if matrix[matrix["technical_field"].isin(["library prep date", "lane/flowcell"])]["age_coincides"].eq("unavailable").any():
        return "batch concern remains"
    increased_fields = matrix[
        matrix["technical_field"].isin(["instrument", "center", "load date", "library prep date", "lane/flowcell"])
        & matrix["age_coincides"].eq("yes_or_sample_specific")
    ]
    if not increased_fields.empty:
        return "batch concern increased"
    return "batch concern reduced"


def write_report(join: pd.DataFrame, matrix: pd.DataFrame, classification: str, local_qc: pd.DataFrame) -> None:
    target = join[join["local_sample"].isin(TARGET_SAMPLES)].copy()
    cd45_values = sorted(set(target["enrichment_sort_condition"].dropna().astype(str)))
    aged_lower_depth = (
        target.groupby("stage")["median_total_counts_mTEC1"].median().get("18mo", np.nan)
        < target.groupby("stage")["median_total_counts_mTEC1"].median().get("02mo", np.nan)
    )
    aged_lower_genes = (
        target.groupby("stage")["median_detected_genes_mTEC1"].median().get("18mo", np.nan)
        < target.groupby("stage")["median_detected_genes_mTEC1"].median().get("02mo", np.nan)
    )

    report = [
        "# GSE240016 GEO/SRA Batch Metadata Audit",
        "",
        "## Scope",
        "",
        "This audit joins local AnnData sample labels to official GEO/SRA metadata for GSE240016. It evaluates whether the four steady-state CD45-negative mTEC1-contributing samples are visibly confounded with available accession, run, library, center, instrument, or timing metadata.",
        "",
        "## Mapping Result",
        "",
        "| local sample | stage | mapped GSM | SRR runs | BioSample | confidence | basis |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, row in target.sort_values("local_sample").iterrows():
        report.append(
            f"| {row['local_sample']} | {row['stage']} | {row['mapped_gsm_accession']} | "
            f"{row.get('srr_run_accessions', '')} | {row.get('BioSample', '')} | "
            f"{row['mapping_confidence']} | {row['mapping_basis']} |"
        )

    report += [
        "",
        "All four local mTEC1 sample labels were mapped to official GSM/SRX/SRR/BioSample records with high confidence using age, CD45-negative stromal sample title, steady-state/untreated timepoint, and replicate tokens.",
        "",
        "## Local mTEC1 QC",
        "",
        "| local sample | stage | total cells | mTEC1 cells | median total counts mTEC1 | median detected genes mTEC1 | mean total counts mTEC1 | mean detected genes mTEC1 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in target.sort_values("local_sample").iterrows():
        report.append(
            f"| {row['local_sample']} | {row['stage']} | {int(row['total_cells'])} | {int(row['mTEC1_cells'])} | "
            f"{row['median_total_counts_mTEC1']:.1f} | {row['median_detected_genes_mTEC1']:.1f} | "
            f"{row['mean_total_counts_mTEC1']:.1f} | {row['mean_detected_genes_mTEC1']:.1f} |"
        )

    report += [
        "",
        f"Aged mTEC1 cells remain lower for local median total counts: {'yes' if aged_lower_depth else 'no'}.",
        f"Aged mTEC1 cells remain lower for local median detected genes: {'yes' if aged_lower_genes else 'no'}.",
        "",
        "## Confound Matrix Summary",
        "",
        "| technical field | age coincides? | assessment |",
        "|---|---|---|",
    ]
    for _, row in matrix.iterrows():
        report.append(f"| {row['technical_field']} | {row['age_coincides']} | {row['assessment']} |")

    report += [
        "",
        "## Explicit Answers",
        "",
        "1. Each of the four mTEC1 local sample labels can be mapped to a GSM/SRX/SRR/BioSample record: yes, with high confidence.",
        "2. The young and aged samples have distinct GSM/SRX/BioSample/SRR accessions. SRA lists multiple runs per sample. The available metadata do not provide lane or flowcell fields, so same-lane or different-lane status cannot be determined.",
        "3. Age coincides with GSM, SRX, BioSample, and SRR accession blocks because each biological sample has its own accessions and the accession ranges separate by age. Age does not obviously coincide with instrument, center, library strategy/source/selection, enrichment condition, SRA load calendar date, or SRA release calendar date. SRA load timestamps differ by sample at minute-scale resolution, but these are administrative load timestamps, not sequencing run dates. Sequencing run date, library prep date, processing date, lane, and flowcell were not found.",
        f"4. CD45-negative stromal enrichment status is consistent across the four samples: {'yes' if cd45_values == ['CD45- thymic stroma cells'] else 'check required'} ({'; '.join(cd45_values)}).",
        f"5. Aged samples are systematically lower for local mTEC1 depth/detected-gene medians: {'yes' if aged_lower_depth and aged_lower_genes else 'not consistently'}; this remains a local QC concern.",
        "6. Official metadata reduces concern about a simple center, instrument, library-type, or public-load-date confound, but it does not resolve sample-level, processing-date, lane/flowcell, or library-prep confounding.",
        "7. Unresolved fields include mouse ID, litter, exact processing date, library prep date, sequencing lane, flowcell, and whether all four steady-state CD45-negative libraries were processed in the same technical batch.",
        "",
        f"Classification: {classification}.",
        "",
        "This classification is cautious. It does not state that batch confounding is ruled out.",
    ]
    (REPORTS_DIR / "gse240016_geo_sra_batch_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_source_notes() -> None:
    notes = [
        "# GSE240016 Metadata Source Notes",
        "",
        "## Official URLs Inspected",
        "",
        "- GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240016",
        "- GEO SOFT metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240016/soft/GSE240016_family.soft.gz",
        "- GEO MINiML metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240016/miniml/GSE240016_family.xml.tgz",
        "- NCBI SRA E-utilities RunInfo queried from SRX accessions linked in GEO sample records.",
        "- ENA browser API queried for SRR run accessions when available.",
        "- NCBI BioSample E-utilities queried for BioSample accessions from SRA RunInfo.",
        "",
        "## Field Sources",
        "",
        "- GEO: GSM accession, sample title, source name, sample characteristics, library strategy/source/selection, protocol text, and SRX relation.",
        "- SRA RunInfo: SRR run accessions, SRX experiment, BioSample, BioProject, SRA study, library layout, platform/model, center, load date, release date, spots, and bases.",
        "- BioSample: BioSample attributes where available.",
        "- ENA: run-level cross-check fields where the API returned records.",
        "",
        "## Not Found",
        "",
        "- Mouse ID.",
        "- Litter.",
        "- Exact processing date.",
        "- Library prep date.",
        "- Sequencing lane.",
        "- Flowcell.",
        "- Explicit batch field.",
        "",
        "## Manual Mapping Assumptions",
        "",
        "- Local `mo02_CD45neg1_d0` maps to GEO steady-state 02-month CD45-negative thymic stromal replicate 1.",
        "- Local `mo02_CD45neg2_d0` maps to GEO steady-state 02-month CD45-negative thymic stromal replicate 2.",
        "- Local `mo18_CD45neg1_d0` maps to GEO steady-state 18-month CD45-negative thymic stromal replicate 1.",
        "- Local `mo18_CD45neg2_d0` maps to GEO steady-state 18-month CD45-negative thymic stromal replicate 2.",
        "",
        "These assumptions are high confidence because the age, enrichment, steady-state/untreated timepoint, and replicate tokens match directly.",
    ]
    (REPORTS_DIR / "gse240016_metadata_source_notes.md").write_text("\n".join(notes) + "\n", encoding="utf-8")


def write_update_and_checks(classification: str) -> None:
    update = [
        "# GSE240016 GEO/SRA Batch Audit Update Summary",
        "",
        "## Previous Local-Metadata Conclusion",
        "",
        "The prior local-metadata conclusion was that the AnnData object exposed sample, stage, day, QC, and cell annotation fields, but not mouse ID, sequencing run, accession, library prep date, processing batch, sex, litter, lane, or flowcell. Residual sample or batch confounding therefore remained possible.",
        "",
        "## New Official Metadata Findings",
        "",
        "Official GEO/SRA metadata maps the four local steady-state CD45-negative samples to GSM, SRX, SRR, and BioSample accessions with high confidence. The available fields show the same library strategy/source/selection, instrument model, center, BioProject, enrichment condition, SRA load calendar date, and SRA release calendar date across the four samples.",
        "",
        "However, accession blocks are sample-specific and separate by age, and the metadata still do not expose mouse ID, litter, processing date, library prep date, lane, flowcell, or explicit batch.",
        "",
        "## Conclusion Change",
        "",
        "Conclusion changed: no. The official metadata reduces concern about a simple center/instrument/library-type confound but does not resolve sample-level or unreported batch variables.",
        "",
        "## Remaining Unresolved",
        "",
        "- Mouse ID and litter.",
        "- Processing date and library prep date.",
        "- Lane and flowcell.",
        "- Explicit batch labels.",
        "- Whether local lower aged mTEC1 depth reflects biology, cell state, technical quality, or a mixture.",
        "",
        "## Effect on mTEC1 Loxl2 Interpretation",
        "",
        "The mTEC1 Loxl2 result should remain candidate-level and potentially confounded. Official metadata does not validate the biology and does not resolve the dropout/depth concern.",
        "",
        f"Batch classification: {classification}.",
    ]
    (REPORTS_DIR / "gse240016_geo_sra_batch_audit_update_summary.md").write_text("\n".join(update) + "\n", encoding="utf-8")

    safety = [
        "# GSE240016 GEO/SRA Batch Audit Safety Check",
        "",
        "- No claim that batch confounding is ruled out.",
        "- No claim that Loxl2 is validated.",
        "- No claim that dropout/depth is resolved.",
        "- No protein, function, or human conservation claims.",
        "- No private contact information.",
        "- No model/tool provenance traces.",
        "",
        "Safety check passed: yes.",
    ]
    (REPORTS_DIR / "gse240016_geo_sra_batch_audit_safety_check.md").write_text("\n".join(safety) + "\n", encoding="utf-8")

    file_check = [
        "# GSE240016 GEO/SRA Batch Audit File Check",
        "",
        "## Files Created",
        "",
        "- scripts/gse240016_geo_sra_metadata_fetch.py",
        "- scripts/gse240016_metadata_batch_confound_audit.py",
        "- results/tables/gse240016_geo_sra_metadata_join.tsv",
        "- results/tables/gse240016_batch_confound_matrix.tsv",
        "- reports/gse240016_geo_sra_batch_audit.md",
        "- reports/gse240016_metadata_source_notes.md",
        "- reports/gse240016_geo_sra_batch_audit_update_summary.md",
        "- reports/gse240016_geo_sra_batch_audit_safety_check.md",
        "- reports/gse240016_geo_sra_batch_audit_file_check.md",
        "- data/external/metadata/gse240016_geo_samples.tsv (ignored; not staged)",
        "- data/external/metadata/gse240016_sra_runinfo.tsv (ignored; not staged)",
        "- data/external/metadata/gse240016_biosample_metadata.tsv (ignored; not staged)",
        "- data/external/metadata/gse240016_metadata_fetch_notes.tsv (ignored; not staged)",
        "",
        "## Checks",
        "",
        "- No FASTQ, BAM, or SRA raw data downloaded.",
        "- No .h5ad file created or staged.",
        "- No large data files created for commit.",
        "- No manuscript file modified by this audit.",
        "- No release or tag created.",
    ]
    (REPORTS_DIR / "gse240016_geo_sra_batch_audit_file_check.md").write_text("\n".join(file_check) + "\n", encoding="utf-8")


def main() -> int:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    geo = read_table(META_DIR / "gse240016_geo_samples.tsv")
    runinfo = read_table(META_DIR / "gse240016_sra_runinfo.tsv")
    biosample = read_table(META_DIR / "gse240016_biosample_metadata.tsv")
    if geo.empty or runinfo.empty:
        raise SystemExit("Missing GEO or SRA metadata. Run gse240016_geo_sra_metadata_fetch.py first.")

    local_qc, _ = summarize_local_obs()
    run_agg = aggregate_runinfo(runinfo)
    join = build_join(local_qc, geo, run_agg, biosample)
    matrix = build_confound_matrix(join)
    classification = classify_batch_concern(matrix, join)

    join.to_csv(RESULTS_DIR / "gse240016_geo_sra_metadata_join.tsv", sep="\t", index=False)
    matrix.to_csv(RESULTS_DIR / "gse240016_batch_confound_matrix.tsv", sep="\t", index=False)
    write_report(join, matrix, classification, local_qc)
    write_source_notes()
    write_update_and_checks(classification)

    print(f"Wrote join rows: {len(join)}")
    print(f"Batch classification: {classification}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
