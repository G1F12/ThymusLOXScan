"""Audit whether E-MTAB-8560 per-mouse LOX reanalysis is feasible."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
META_DIR = PROJECT_ROOT / "data" / "external" / "metadata" / "emtab8560"
REPORT = PROJECT_ROOT / "reports" / "external_emtab8560_biological_unit_audit.md"
TABLE = PROJECT_ROOT / "results" / "tables" / "external_emtab8560_biological_unit_audit.tsv"


def main() -> int:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    TABLE.parent.mkdir(parents=True, exist_ok=True)

    coldata_path = META_DIR / "emtab8560_coldata.tsv"
    sdrf_path = META_DIR / "emtab8560_sdrf.tsv"

    if not coldata_path.exists():
        rows = [
            {
                "question": "classification",
                "answer": "infeasible in current environment",
                "evidence": "Rscript was unavailable or R export did not produce emtab8560_coldata.tsv.",
            },
            {
                "question": "true biological mouse/individual ID",
                "answer": "not verified",
                "evidence": "Official R object colData could not be exported through MouseThymusAgeing API.",
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
                    "Rscript was not available, so the official MouseThymusAgeing R API could not be used to export colData and assay layers. The SDRF metadata can identify candidate fields such as `Characteristics[individual]`, but this does not verify that every analyzed cell in the official processed object maps to a true biological mouse/individual.",
                    "",
                    "No per-mouse inference was performed. SortDay, PlateID, and run fields were not substituted for biological mouse ID.",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        print("Biological unit audit: infeasible in current environment")
        return 0

    coldata = pd.read_csv(coldata_path, sep="\t", dtype=str).fillna("")
    sdrf = pd.read_csv(sdrf_path, sep="\t", dtype=str).fillna("") if sdrf_path.exists() else pd.DataFrame()
    candidate_fields = [c for c in coldata.columns if any(k in c.lower() for k in ["individual", "mouse", "animal"])]
    classification = "not verified: no true mouse/individual ID available"
    if candidate_fields:
        classification = "partially verified: mouse IDs available but batch/age partly confounded"

    rows = [
        {"question": "classification", "answer": classification, "evidence": ";".join(candidate_fields)},
        {"question": "field present in R-exported colData", "answer": "yes" if candidate_fields else "no", "evidence": ";".join(candidate_fields)},
        {"question": "field present in SDRF", "answer": "yes" if "Characteristics[individual]" in sdrf.columns else "no", "evidence": "Characteristics[individual]"},
    ]
    pd.DataFrame(rows).to_csv(TABLE, sep="\t", index=False)
    REPORT.write_text(
        "# E-MTAB-8560 Biological Unit Audit\n\n"
        f"Classification: {classification}.\n\n"
        "A downstream model should only proceed after confirming the candidate field maps every analyzed cell to a biological mouse/individual.\n",
        encoding="utf-8",
    )
    print(f"Biological unit audit: {classification}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
