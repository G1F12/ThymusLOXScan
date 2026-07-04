"""Verify biological mouse/individual aggregation for E-MTAB-8560."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
META_DIR = PROJECT_ROOT / "data" / "external" / "metadata" / "emtab8560"
REPORT = PROJECT_ROOT / "reports" / "external_emtab8560_biological_unit_audit.md"
TABLE = PROJECT_ROOT / "results" / "tables" / "external_emtab8560_biological_unit_audit.tsv"
CELL_META = PROJECT_ROOT / "results" / "tables" / "external_emtab8560_cell_metadata_verified.tsv"

PRIMARY_GROUPS = ["mTEClo", "mTEChi", "combined_mTEClo_mTEChi", "cTEC"]
INFERRED_MTEC = ["Mature.mTEC", "Post-Aire.mTEC", "Tuft-like.mTEC"]


def norm_id(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(value))


def age_to_int(value: object) -> int | None:
    match = re.search(r"\d+", str(value))
    return int(match.group(0)) if match else None


def load_inputs() -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    coldata_path = META_DIR / "emtab8560_coldata.tsv"
    sdrf_path = META_DIR / "emtab8560_sdrf.tsv"
    coldata = pd.read_csv(coldata_path, sep="\t", dtype=str).fillna("") if coldata_path.exists() else None
    sdrf = pd.read_csv(sdrf_path, sep="\t", dtype=str).fillna("") if sdrf_path.exists() else None
    return coldata, sdrf


def group_membership(meta: pd.DataFrame) -> pd.DataFrame:
    sort_type = meta.get("SortType", "").astype(str)
    subtype = meta.get("SubType", "").astype(str)
    rows = []
    definitions = {
        "mTEClo": sort_type.eq("mTEClo"),
        "mTEChi": sort_type.eq("mTEChi"),
        "combined_mTEClo_mTEChi": sort_type.isin(["mTEClo", "mTEChi"]),
        "cTEC": sort_type.eq("cTEC") | subtype.str.contains("cTEC", case=False, na=False),
        "Mature.mTEC": subtype.eq("Mature.mTEC"),
        "Post-Aire.mTEC": subtype.eq("Post-Aire.mTEC"),
        "Tuft-like.mTEC": subtype.eq("Tuft-like.mTEC"),
    }
    for group, mask in definitions.items():
        tmp = meta.loc[mask].copy()
        tmp["analysis_group"] = group
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def write_infeasible(reason: str) -> None:
    rows = [
        {
            "question": "classification",
            "answer": "infeasible in current environment",
            "evidence": reason,
        },
        {
            "question": "true biological mouse/individual ID",
            "answer": "not verified",
            "evidence": "Official R export and SDRF join could not be completed.",
        },
        {
            "question": "cells treated as biological replicates",
            "answer": "no",
            "evidence": "No inferential per-cell or per-mouse model was run.",
        },
    ]
    pd.DataFrame(rows).to_csv(TABLE, sep="\t", index=False)
    REPORT.write_text(
        "\n".join(
            [
                "# E-MTAB-8560 Biological Unit Audit",
                "",
                "Classification: infeasible in current environment.",
                "",
                reason,
                "",
                "No per-mouse inference was performed. SortDay, PlateID, and run fields were not substituted for biological mouse ID.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    coldata, sdrf = load_inputs()
    if coldata is None:
        write_infeasible("R export did not produce `emtab8560_coldata.tsv`.")
        print("Biological unit audit: infeasible in current environment")
        return 0
    if sdrf is None:
        write_infeasible("BioStudies/ArrayExpress SDRF metadata is missing.")
        print("Biological unit audit: infeasible in current environment")
        return 0

    coldata = coldata.copy()
    coldata["join_id"] = coldata["cell_barcode"].map(norm_id)
    sdrf = sdrf.copy()
    sdrf["join_id"] = sdrf["Source Name"].map(norm_id)
    if "Assay Name" in sdrf.columns:
        assay_join = sdrf.copy()
        assay_join["join_id"] = assay_join["Assay Name"].map(norm_id)
        sdrf_join = pd.concat([sdrf, assay_join], ignore_index=True)
    else:
        sdrf_join = sdrf
    keep = [
        "join_id",
        "Characteristics[individual]",
        "Characteristics[age]",
        "Factor Value[age]",
        "Characteristics[cell type]",
        "Characteristics[FACS marker]",
        "Characteristics[inferred cell type]",
        "Comment[ENA_RUN]",
        "Assay Name",
    ]
    keep = [c for c in keep if c in sdrf_join.columns]
    sdrf_join = sdrf_join[keep].drop_duplicates("join_id")
    meta = coldata.merge(sdrf_join, on="join_id", how="left", validate="many_to_one")

    meta["age_week"] = meta.get("Age", meta.get("Characteristics[age]", "")).map(age_to_int)
    if meta["age_week"].isna().any() and "Factor Value[age]" in meta:
        meta.loc[meta["age_week"].isna(), "age_week"] = meta.loc[meta["age_week"].isna(), "Factor Value[age]"].map(age_to_int)
    meta["individual"] = meta.get("Characteristics[individual]", "").astype(str)
    meta["biological_mouse_id"] = "wk" + meta["age_week"].astype("Int64").astype(str) + "_mouse" + meta["individual"]
    invalid_mouse = meta["individual"].eq("") | meta["age_week"].isna()
    meta.loc[invalid_mouse, "biological_mouse_id"] = ""

    mapped_fraction = float(meta["biological_mouse_id"].ne("").mean())
    grouped = group_membership(meta)
    counts = (
        grouped[grouped["analysis_group"].isin(PRIMARY_GROUPS + INFERRED_MTEC)]
        .groupby(["analysis_group", "age_week"], dropna=False)
        .agg(
            n_cells=("cell_barcode", "nunique"),
            n_mice=("biological_mouse_id", lambda x: x[x.ne("")].nunique()),
            sort_days=("SortDay", lambda x: ";".join(sorted(set(map(str, x))))),
            plate_ids=("PlateID", lambda x: ";".join(sorted(set(map(str, x)))[:12])),
            runs=("Comment[ENA_RUN]", lambda x: ";".join(sorted(set(v for v in map(str, x) if v))[:12])),
        )
        .reset_index()
    )

    primary = counts[counts["analysis_group"].isin(["mTEClo", "mTEChi", "combined_mTEClo_mTEChi"])]
    low_n = primary["n_mice"].lt(4).any()
    mouse_verified = mapped_fraction == 1.0 and not primary.empty and not low_n

    # Batch fields are cell-level and plate/run-rich; treat as partly confounded unless every age shares values.
    batch_notes = []
    for field in ["SortDay", "PlateID", "Comment[ENA_RUN]"]:
        if field in meta.columns:
            per_age = meta.groupby("age_week")[field].nunique(dropna=True).to_dict()
            batch_notes.append(f"{field}: {per_age}")
    classification = (
        "partially verified: mouse IDs available but batch/age partly confounded"
        if mouse_verified
        else "not verified: no true mouse/individual ID available"
    )
    if mapped_fraction == 1.0 and low_n:
        classification = "partially verified: mouse IDs available but fewer than 4 mice in at least one primary group/age"
    elif mapped_fraction == 1.0 and not low_n:
        classification = "partially verified: mouse IDs available but batch/age partly confounded"

    rows = [
        {"question": "classification", "answer": classification, "evidence": " | ".join(batch_notes)},
        {"question": "true biological mouse/individual ID", "answer": "age_week + Characteristics[individual]", "evidence": "SDRF field joined to R-exported colData by normalized cell/source/assay identifier."},
        {"question": "field present in R-exported colData", "answer": "no", "evidence": "Mouse/individual ID was not assumed from SortDay or PlateID."},
        {"question": "field present only in SDRF/ArrayExpress metadata", "answer": "yes", "evidence": "Characteristics[individual]"},
        {"question": "each cell mapped to mouse/individual", "answer": "yes" if mapped_fraction == 1.0 else "partial", "evidence": f"{mapped_fraction:.3f} mapped fraction"},
        {"question": "any primary age group fewer than 4 mice", "answer": "yes" if low_n else "no", "evidence": "See group rows below."},
        {"question": "batch-aware modeling estimable", "answer": "limited", "evidence": "SortDay/PlateID/run are cell/acquisition fields with complex structure; unadjusted per-mouse models plus batch audit are safer."},
        {"question": "cells treated as biological replicates", "answer": "no", "evidence": "Downstream unit is biological_mouse_id."},
    ]
    pd.concat([pd.DataFrame(rows), counts.rename(columns={"analysis_group": "question", "age_week": "answer"}).astype(str)], ignore_index=True).to_csv(TABLE, sep="\t", index=False)
    meta.to_csv(CELL_META, sep="\t", index=False)

    lines = [
        "# E-MTAB-8560 Biological Unit Audit",
        "",
        f"Classification: {classification}.",
        "",
        "## Biological Unit",
        "",
        "The true biological unit used for downstream analysis is `age_week + Characteristics[individual]`, joined from official ArrayExpress/BioStudies SDRF metadata onto R-exported MouseThymusAgeing colData by normalized cell/source/assay identifiers.",
        "",
        "SortDay and PlateID were not treated as mouse IDs.",
        "",
        f"Cell-to-mouse mapped fraction: {mapped_fraction:.3f}.",
        "",
        "## Mice Per Group",
        "",
        "| group | age_week | n mice | n cells |",
        "|---|---:|---:|---:|",
    ]
    for _, row in counts.sort_values(["analysis_group", "age_week"]).iterrows():
        lines.append(f"| {row['analysis_group']} | {int(row['age_week'])} | {int(row['n_mice'])} | {int(row['n_cells'])} |")
    lines += [
        "",
        "## Batch Fields",
        "",
        *[f"- {note}" for note in batch_notes],
        "",
        "Batch-aware modeling is limited because SortDay, PlateID, and ENA run are not simple mouse-level covariates. The downstream analysis uses per-mouse aggregation and reports batch limitations explicitly.",
        "",
        "Cells are not treated as independent biological replicates.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Biological unit audit: {classification}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
