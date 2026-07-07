"""Pilot matrix-to-metadata join for GSE231906."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from human_thymus_stage4_gse231906_metadata_audit import (
    PILOT_DIR,
    REPORTS,
    REQUIRED_NO_OVERCLAIM,
    ROOT,
    TABLES,
    TARGET_GENES,
    load_metadata,
    rel,
    write_tsv,
)


PILOT_FILES = TABLES / "human_thymus_stage4_gse231906_pilot_extracted_files.tsv"


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def read_barcodes(path: Path) -> list[str]:
    values: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.strip():
                values.append(line.rstrip("\n").split("\t")[0].split(",")[0])
    return values


def read_features(path: Path) -> list[str]:
    genes: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            gene = parts[1] if len(parts) > 1 and parts[1] else parts[0]
            genes.append(gene)
    return genes


def norm_barcode(value: object) -> str:
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def core_barcode(value: object) -> str:
    text = norm_barcode(value)
    if "_donor" in text:
        text = text.split("_donor", 1)[0]
    return text


def strip_suffix(value: object) -> str:
    text = core_barcode(value)
    if "-" in text:
        head, tail = text.rsplit("-", 1)
        if tail.isdigit():
            return head
    return text


def load_pilot_paths() -> dict[str, Path]:
    if not PILOT_FILES.exists():
        return {}
    files = pd.read_csv(PILOT_FILES, sep="\t").fillna("")
    paths: dict[str, Path] = {}
    for _, row in files.iterrows():
        role = str(row["extraction_role"])
        path = Path(row["local_path"])
        if not path.is_absolute():
            path = ROOT / path
        paths[role] = path
    return paths


def match_metadata(barcodes: list[str], metadata: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    expr = pd.DataFrame({"expression_barcode": barcodes})
    expr["barcode_exact"] = expr["expression_barcode"].map(core_barcode)
    expr["barcode_stripped"] = expr["expression_barcode"].map(strip_suffix)
    meta = metadata.copy()
    meta["metadata_barcode_exact"] = meta["cell_id"].map(core_barcode)
    meta["metadata_barcode_stripped"] = meta["cell_id"].map(strip_suffix)

    strategies = [
        ("exact", "barcode_exact", "metadata_barcode_exact"),
        ("strip_numeric_suffix", "barcode_stripped", "metadata_barcode_stripped"),
    ]
    best_join = pd.DataFrame()
    best_rate = -1.0
    best_summary: dict[str, object] = {
        "join_strategy": "none",
        "total_expression_cells": len(expr),
        "matched_metadata_cells": 0,
        "unmatched_expression_cells": len(expr),
        "unmatched_metadata_cells": len(meta),
        "match_rate": 0.0,
    }
    for name, expr_col, meta_col in strategies:
        dedup_meta = meta.loc[meta[meta_col].ne("")].drop_duplicates([meta_col, "sample_id", "donor_id", "fine_cell_type"])
        joined = expr.merge(dedup_meta, left_on=expr_col, right_on=meta_col, how="left")
        matched_expr = joined.loc[joined["cell_id"].notna(), "expression_barcode"].nunique()
        rate = matched_expr / len(expr) if len(expr) else 0.0
        if rate > best_rate:
            best_join = joined
            best_rate = rate
            best_summary = {
                "join_strategy": name,
                "total_expression_cells": len(expr),
                "matched_metadata_cells": int(matched_expr),
                "unmatched_expression_cells": int(len(expr) - matched_expr),
                "unmatched_metadata_cells": int(max(len(meta) - joined["cell_id"].nunique(), 0)),
                "match_rate": rate,
            }
    return best_join, best_summary


def write_preliminary_reports(summary: dict[str, object], genes_present: list[str], target_compartments: list[str], pilot_pass: bool) -> None:
    grade = "weak_or_limited_context" if pilot_pass else "not_evaluable"
    recommendation = (
        "Continue to target-gene extraction before considering any future supplement wording."
        if pilot_pass
        else "Archive/schema follow-up is recommended before expression interpretation."
    )
    lines = [
        "# GSE231906 Balanced Report",
        "",
        "## Purpose",
        "",
        "Assess whether GSE231906 can provide donor/sample-aware aged-human thymus transcript-level context for LOX-family genes.",
        "",
        "## Metadata Audit",
        "",
        "Metadata workbook parsing was attempted before expression interpretation.",
        "",
        "## Matrix-to-Metadata Join Result",
        "",
        f"- Pilot join successful: {'yes' if pilot_pass else 'no'}",
        f"- Join match rate: {float(summary.get('match_rate', 0.0)):.3f}",
        f"- Join strategy: {summary.get('join_strategy', 'none')}",
        "",
        "## Gene Presence",
        "",
        f"- LOX-family genes present in pilot feature list: {', '.join(genes_present) if genes_present else 'none'}",
        "",
        "## Target Compartments Found",
        "",
        f"- {', '.join(target_compartments) if target_compartments else 'none after pilot join'}",
        "",
        "## Strengths",
        "",
        "- The metadata and archive workflow is reproducible and keeps raw files in the external-data path.",
        "- The pilot join directly tests barcode linkage before any donor/sample-aware summaries.",
        "",
        "## Limitations",
        "",
        "- Full LOX-family donor/sample-aware detection summaries have not yet been generated at this checkpoint.",
        "- Mean-value interpretation remains withheld until matrix semantics are clearer.",
        "",
        "## Evidence Grade",
        "",
        grade,
        "",
        "## Recommendation",
        "",
        recommendation,
        "",
        "## No-Overclaim Statement",
        "",
        REQUIRED_NO_OVERCLAIM,
    ]
    (REPORTS / "human_thymus_stage4_gse231906_balanced_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (REPORTS / "human_thymus_stage4_gse231906_decision_memo.md").write_text(
        "\n".join([
            "# GSE231906 Decision Memo",
            "",
            f"- Was GSE231906 technically parsed? {'partial' if pilot_pass else 'no'}",
            f"- Was metadata join successful? {'yes' if pilot_pass else 'no'}",
            f"- Were LOX-family genes detected? {'yes: ' + ', '.join(genes_present) if genes_present else 'no'}",
            f"- Were mTEC-like/cTEC-like/Epi compartments available? {'yes: ' + ', '.join(target_compartments) if target_compartments else 'no'}",
            f"- Was donor/sample-aware aging context feasible? {'pending target extraction' if pilot_pass else 'no'}",
            f"- Evidence grade: {grade}",
            "- Should this be mentioned in future manuscript supplement? conditional",
            "- Should this trigger immediate manuscript update? no",
            "- Should this trigger release now? no",
            "- Should wet-lab proposal continue? yes",
            "- Best balanced wording for outreach: GSE231906 is being evaluated as aged-human thymus transcript-level context, with no validation claim.",
            "- Best balanced wording for a future manuscript supplement: pending donor/sample-aware LOX-family summaries.",
        ]) + "\n",
        encoding="utf-8",
    )
    (REPORTS / "human_thymus_stage4_gse231906_one_page_summary.md").write_text(
        "\n".join([
            "# GSE231906 One-Page Summary",
            "",
            "## What was tested",
            "Metadata, archive structure, pilot extraction, and matrix-to-metadata barcode linkage.",
            "",
            "## What was newly possible",
            "A controlled pilot join checkpoint was created before any biological interpretation.",
            "",
            "## Main result",
            f"Pilot join pass status: {'yes' if pilot_pass else 'no'}, match rate {float(summary.get('match_rate', 0.0)):.3f}.",
            "",
            "## Evidence grade",
            grade,
            "",
            "## How it affects human relevance",
            "This checkpoint affects technical feasibility only.",
            "",
            "## What it does not prove",
            REQUIRED_NO_OVERCLAIM,
            "",
            "## Next step",
            recommendation,
        ]) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    paths = load_pilot_paths()
    metadata, _ = load_metadata()
    barcodes_path = paths.get("barcodes")
    features_path = paths.get("features")
    if not barcodes_path or not features_path or not barcodes_path.exists() or not features_path.exists():
        summary = {
            "join_strategy": "missing_pilot_barcodes_or_features",
            "total_expression_cells": 0,
            "matched_metadata_cells": 0,
            "unmatched_expression_cells": 0,
            "unmatched_metadata_cells": len(metadata),
            "match_rate": 0.0,
        }
        genes_present: list[str] = []
        target_compartments: list[str] = []
        pilot_pass = False
    else:
        barcodes = read_barcodes(barcodes_path)
        genes = read_features(features_path)
        genes_present = [gene for gene in TARGET_GENES if gene in set(map(str.upper, genes))]
        joined, summary = match_metadata(barcodes, metadata)
        target_compartments = sorted(
            joined.loc[joined["cell_id"].notna(), "target_compartments"]
            .dropna()
            .astype(str)
            .str.split(";")
            .explode()
            .loc[lambda s: s.ne("other")]
            .unique()
            .tolist()
        )
        required_fields = ["donor_id", "sample_id", "age_years", "fine_cell_type"]
        usable_fields = all(field in joined.columns and joined[field].notna().any() for field in required_fields)
        pilot_pass = bool(summary["match_rate"] >= 0.80 and genes_present and usable_fields and target_compartments)

    summary_row = {
        **summary,
        "match_rate": f"{float(summary.get('match_rate', 0.0)):.6f}",
        "fields_available_after_join": "donor_id;sample_id;age_years;sex;broad_compartment;fine_cell_type;target_compartments",
        "target_compartment_availability_after_join": ";".join(target_compartments),
        "pilot_pass": "yes" if pilot_pass else "no",
    }
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_join_pilot_summary.tsv",
        [summary_row],
        ["join_strategy", "total_expression_cells", "matched_metadata_cells", "unmatched_expression_cells", "unmatched_metadata_cells", "match_rate", "fields_available_after_join", "target_compartment_availability_after_join", "pilot_pass"],
    )
    gene_rows = [{
        "target_gene": gene,
        "present_in_pilot_expression_unit": "yes" if gene in genes_present else "no",
        "feature_match": gene if gene in genes_present else "",
    } for gene in TARGET_GENES]
    write_tsv(
        TABLES / "human_thymus_stage4_gse231906_gene_presence_pilot.tsv",
        gene_rows,
        ["target_gene", "present_in_pilot_expression_unit", "feature_match"],
    )
    (REPORTS / "human_thymus_stage4_gse231906_join_pilot_report.md").write_text(
        "\n".join([
            "# GSE231906 Join Pilot Report",
            "",
            f"- Total expression cells: {summary_row['total_expression_cells']}",
            f"- Matched metadata cells: {summary_row['matched_metadata_cells']}",
            f"- Unmatched expression cells: {summary_row['unmatched_expression_cells']}",
            f"- Unmatched metadata cells: {summary_row['unmatched_metadata_cells']}",
            f"- Match rate: {summary_row['match_rate']}",
            f"- Join strategy: {summary_row['join_strategy']}",
            f"- LOX-family genes present: {', '.join(genes_present) if genes_present else 'none'}",
            f"- Target compartments after join: {', '.join(target_compartments) if target_compartments else 'none'}",
            f"- Proceed to target-gene extraction: {'yes' if pilot_pass else 'no'}",
            "",
            "No biological inference or cell-level statistical testing was performed in this join pilot.",
        ]) + "\n",
        encoding="utf-8",
    )
    write_preliminary_reports(summary, genes_present, target_compartments, pilot_pass)
    print(f"join_match_rate={summary_row['match_rate']}")
    print(f"pilot_pass={summary_row['pilot_pass']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
