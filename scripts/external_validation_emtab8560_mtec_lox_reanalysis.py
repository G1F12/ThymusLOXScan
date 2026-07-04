"""Run E-MTAB-8560 mTEC LOX analysis only after verified R export."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
META_DIR = PROJECT_ROOT / "data" / "external" / "metadata" / "emtab8560"
REPORT = PROJECT_ROOT / "reports" / "external_emtab8560_mtec_lox_reanalysis.md"


def main() -> int:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    required = [
        META_DIR / "emtab8560_coldata.tsv",
        META_DIR / "emtab8560_lox_counts.tsv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        REPORT.write_text(
            "\n".join(
                [
                    "# E-MTAB-8560 mTEC LOX Reanalysis",
                    "",
                    "Final classification: not reliable/infeasible.",
                    "",
                    "The official MouseThymusAgeing R export was not available in this environment, so per-mouse pseudobulk, permutation tests, trend models, FDR, and figures were not generated.",
                    "",
                    "Missing required files:",
                    *[f"- `{item}`" for item in missing],
                    "",
                    "No cells were treated as independent biological replicates.",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        print("E-MTAB-8560 mTEC LOX reanalysis infeasible: missing R export")
        return 0
    raise SystemExit("R export exists, but modeling implementation is intentionally not run without biological-unit verification.")


if __name__ == "__main__":
    raise SystemExit(main())
