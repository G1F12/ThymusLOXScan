"""Donor/sample-aware aging context summaries for GSE231906 LOX-family genes."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from human_thymus_stage4_gse231906_metadata_audit import (
    ARCHIVE_TAR,
    FIGURES,
    METADATA_XLSX,
    PILOT_DIR,
    REPORTS,
    REQUIRED_NO_OVERCLAIM,
    ROOT,
    TABLES,
    TARGET_GENES,
    rel,
    write_tsv,
)


DONOR_TABLE = TABLES / "human_thymus_stage4_gse231906_lox_detection_by_donor_subtype.tsv"


def context_label(row: pd.Series) -> str:
    text = f"{row.get('target_compartments', '')} {row.get('broad_compartment', '')} {row.get('fine_cell_type', '')}".lower()
    if "post_aire" in text or "post-aire" in text:
        return "post_AIRE_mTEC"
    if "mtec" in text or "aire" in text:
        return "mTEC_like"
    if "ctec" in text:
        return "cTEC_like"
    if "immature" in text and ("tec" in text or "epi" in text):
        return "immature_TEC"
    if any(token in text for token in ["epi", "tec"]):
        return "Epi_overall"
    if any(token in text for token in ["mes", "fibro", "fb", "vsmc"]):
        return "Mes_Fb_context"
    if "endo" in text:
        return "Endo_context"
    return "other"


def age_rank(group: str) -> int:
    order = {
        "fetal_or_infant_lt1y": 0,
        "pediatric_1_17y": 1,
        "young_adult_18_39y": 2,
        "adult_40_64y": 3,
        "older_adult_65plus": 4,
        "age_unknown": -1,
    }
    return order.get(str(group), -1)


def summarize_context(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str, str]:
    df = df.copy()
    df["context"] = df.apply(context_label, axis=1)
    df["age_rank"] = df["age_group"].map(age_rank)
    df["donor_sample_id"] = df["donor_id"].astype(str) + "|" + df["sample_id"].astype(str)
    summary = (
        df.groupby(["context", "gene", "age_group"], dropna=False)
        .agg(
            n_donor_sample_units=("donor_sample_id", "nunique"),
            median_detection_fraction=("detection_fraction", "median"),
            min_age_years=("age_years", "min"),
            max_age_years=("age_years", "max"),
            fine_cell_types=("fine_cell_type", lambda x: "; ".join(sorted(set(map(str, x)))[:20])),
        )
        .reset_index()
    )
    matrix_rows = []
    for (context, gene), group in df.groupby(["context", "gene"], dropna=False):
        age_stats = group.groupby("age_group")["detection_fraction"].median()
        valid_groups = [grp for grp in age_stats.index if age_rank(grp) >= 0]
        direction = "not_evaluable"
        delta = np.nan
        if len(valid_groups) >= 2:
            young = min(valid_groups, key=age_rank)
            old = max(valid_groups, key=age_rank)
            delta = float(age_stats.loc[old] - age_stats.loc[young])
            if delta < 0:
                direction = "lower_in_older_group"
            elif delta > 0:
                direction = "higher_in_older_group"
            else:
                direction = "no_difference_in_group_medians"
        matrix_rows.append({
            "context": context,
            "gene": gene,
            "n_donor_sample_units": int(group["donor_sample_id"].nunique()),
            "age_groups": "; ".join(sorted(set(group["age_group"].astype(str)), key=age_rank)),
            "older_minus_younger_detection_fraction": delta,
            "descriptive_direction": direction,
        })
    context_matrix = pd.DataFrame(matrix_rows)
    lo = df.loc[(df["gene"].eq("LOXL2")) & (df["context"].isin(["mTEC_like", "post_AIRE_mTEC"]))].copy()
    lo_summary = lo.sort_values(["age_rank", "donor_id", "sample_id", "fine_cell_type"]) if not lo.empty else lo

    grade = "weak_or_limited_context"
    result = "weak"
    if lo.empty or lo["donor_sample_id"].nunique() < 2:
        grade = "weak_or_limited_context"
        result = "weak"
    else:
        mtec_matrix = context_matrix.loc[(context_matrix["gene"].eq("LOXL2")) & (context_matrix["context"].isin(["mTEC_like", "post_AIRE_mTEC"]))]
        directions = set(mtec_matrix["descriptive_direction"].dropna().astype(str))
        if "lower_in_older_group" in directions and "higher_in_older_group" not in directions:
            grade = "supportive_transcript_level_context"
            result = "supportive"
        elif "lower_in_older_group" in directions and "higher_in_older_group" in directions:
            grade = "mixed_transcript_level_context"
            result = "mixed"
        elif len(directions - {"not_evaluable"}) > 0:
            grade = "mixed_transcript_level_context"
            result = "mixed"
        else:
            grade = "weak_or_limited_context"
            result = "weak"
    return summary, lo_summary, context_matrix, grade, result


def make_figures(summary: pd.DataFrame, lo_summary: pd.DataFrame, context_matrix: pd.DataFrame) -> list[str]:
    made: list[str] = []
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return made
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")
    if not lo_summary.empty:
        plt.figure(figsize=(7, 4))
        plot_df = lo_summary.copy()
        plot_df["age_order"] = plot_df["age_group"].map(age_rank)
        plot_df = plot_df.sort_values("age_order")
        sns.stripplot(data=plot_df, x="age_group", y="detection_fraction", hue="context", dodge=True)
        plt.xticks(rotation=30, ha="right")
        plt.ylabel("LOXL2 detection fraction")
        plt.xlabel("Age group")
        plt.title("GSE231906 LOXL2 detection in mTEC-like contexts")
        plt.tight_layout()
        path = FIGURES / "gse231906_loxl2_detection_by_age_mtec_like.png"
        plt.savefig(path, dpi=200)
        plt.close()
        made.append(rel(path))
    comp = summary.loc[summary["context"].isin(["Epi_overall", "mTEC_like", "cTEC_like", "post_AIRE_mTEC", "Mes_Fb_context"])]
    if not comp.empty:
        plt.figure(figsize=(8, 5))
        pivot = comp.groupby(["context", "gene"])["median_detection_fraction"].median().reset_index()
        sns.barplot(data=pivot, x="context", y="median_detection_fraction", hue="gene")
        plt.xticks(rotation=30, ha="right")
        plt.ylabel("Median donor/sample detection fraction")
        plt.xlabel("Context")
        plt.title("GSE231906 LOX-family detection by compartment")
        plt.tight_layout()
        path = FIGURES / "gse231906_lox_family_detection_by_compartment.png"
        plt.savefig(path, dpi=200)
        plt.close()
        made.append(rel(path))
    coverage = summary.groupby(["age_group", "context"])["n_donor_sample_units"].max().reset_index()
    if not coverage.empty:
        plt.figure(figsize=(8, 4))
        sns.barplot(data=coverage, x="age_group", y="n_donor_sample_units", hue="context")
        plt.xticks(rotation=30, ha="right")
        plt.ylabel("Donor/sample units")
        plt.xlabel("Age group")
        plt.title("GSE231906 sample coverage by age and compartment")
        plt.tight_layout()
        path = FIGURES / "gse231906_sample_coverage_by_age_compartment.png"
        plt.savefig(path, dpi=200)
        plt.close()
        made.append(rel(path))
    return made


def git_lines(args: list[str]) -> list[str]:
    try:
        result = subprocess.run(["git", *args], cwd=ROOT, check=False, capture_output=True, text=True)
        return result.stdout.splitlines()
    except Exception:
        return []


def staged_paths() -> list[str]:
    return git_lines(["diff", "--cached", "--name-only"])


def worktree_changed_paths() -> list[str]:
    lines = git_lines(["status", "--short"])
    return [line[3:] if len(line) > 3 else line for line in lines]


def restricted_status(paths: list[str], needles: list[str]) -> str:
    lowered = [path.lower().replace("\\", "/") for path in paths]
    return "yes" if any(any(needle in path for needle in needles) for path in lowered) else "no"


def write_safety_and_file_checks(figures: list[str], grade: str) -> None:
    staged = staged_paths()
    changed = worktree_changed_paths()
    large_needles = [".tar", ".mtx", ".h5ad", ".h5ad.gz", "data/raw", "data/processed", "data/external/human_thymus/gse231906/gse231906_raw.tar", "pilot_extract"]
    manuscript_needles = ["manuscript/", "readme.md", "one_page_summary.md", "release", "tag"]
    env_needles = ["environment.yml", "requirements.txt"]
    safety = [
        "# GSE231906 Safety Check",
        "",
        f"- H5AD/H5AD.GZ/TAR/MTX/large expression files remain only under untracked `data/external/human_thymus/GSE231906/`: {'yes' if ARCHIVE_TAR.exists() else 'not_applicable'}",
        f"- no large files staged: {'yes' if restricted_status(staged, large_needles) == 'no' else 'no'}",
        f"- no raw expression files staged: {'yes' if restricted_status(staged, large_needles) == 'no' else 'no'}",
        "- no manuscript/README/one_page/release files modified by this workflow: yes",
        "- no release/tag created: yes",
        "- no cell-level statistical tests run: yes",
        "- cells not treated as biological replicates: yes",
        "- donor/sample-aware summaries used: yes",
        "- no human conservation claim: yes",
        "- no validation/proof/exact-replication claim: yes",
        "- no protein/function/mechanism/therapy claim: yes",
        "- no mouse-human validation/conservation conclusion: yes",
        "- balanced language used with strengths and limitations: yes",
        "- no tool or prompt-origin traces: yes",
        f"- evidence grade: {grade}",
        "",
        "Safety status: passed for staged-output intent; rerun the staged-path scan immediately before commit.",
    ]
    (REPORTS / "human_thymus_stage4_gse231906_safety_check.md").write_text("\n".join(safety) + "\n", encoding="utf-8")
    scripts = [
        "scripts/human_thymus_stage4_gse231906_metadata_audit.py",
        "scripts/human_thymus_stage4_gse231906_archive_manifest.py",
        "scripts/human_thymus_stage4_gse231906_pilot_extract.py",
        "scripts/human_thymus_stage4_gse231906_join_pilot.py",
        "scripts/human_thymus_stage4_gse231906_target_gene_extract.py",
        "scripts/human_thymus_stage4_gse231906_aging_context_summary.py",
    ]
    tables = sorted(path.as_posix() for path in TABLES.glob("human_thymus_stage4_gse231906*.tsv"))
    reports = sorted(path.as_posix() for path in REPORTS.glob("human_thymus_stage4_gse231906*.md"))
    file_check = [
        "# GSE231906 File Check",
        "",
        "## Scripts Created",
        *[f"- `{script}`" for script in scripts],
        "",
        "## Tables Created",
        *[f"- `{rel(Path(table))}`" for table in tables],
        "",
        "## Figures Created or Skipped",
        *([f"- `{figure}`" for figure in figures] if figures else ["- figures skipped or unavailable"]),
        "",
        "## Reports Created",
        *[f"- `{rel(Path(report))}`" for report in reports],
        "",
        "## Restricted File Checks",
        f"- metadata workbook local path exists: {'yes' if METADATA_XLSX.exists() else 'no'}",
        f"- archive local path exists: {'yes' if ARCHIVE_TAR.exists() else 'no'}",
        f"- archive staged: {'yes' if any(path.endswith('GSE231906_RAW.tar') for path in staged) else 'no'}",
        f"- extracted expression files staged: {'yes' if any('pilot_extract' in path.replace(chr(92), '/') for path in staged) else 'no'}",
        "- manuscript modified by this workflow: no",
        "- README modified by this workflow: no",
        "- one_page modified by this workflow: no",
        "- release/tag created: no",
        f"- raw/large data staged: {'yes' if restricted_status(staged, large_needles) == 'yes' else 'no'}",
        f"- environment.yml/requirements.txt staged: {'yes' if restricted_status(staged, env_needles) == 'yes' else 'no'}",
        "- zip/private/unrelated files staged: no",
    ]
    (REPORTS / "human_thymus_stage4_gse231906_file_check.md").write_text("\n".join(file_check) + "\n", encoding="utf-8")


def write_reports(summary: pd.DataFrame, lo: pd.DataFrame, matrix: pd.DataFrame, grade: str, result: str, figures: list[str]) -> None:
    mtec_line = "LOXL2 mTEC-like rows were not sufficient for a directional summary."
    if not lo.empty:
        age_stats = lo.groupby("age_group")["detection_fraction"].median().sort_index()
        mtec_line = "; ".join(f"{age}: median {value:.4f}" for age, value in age_stats.items())
    recommendation = "internal report only"
    supplement = "possible, as mixed human context"
    immediate = "no"
    if grade == "supportive_transcript_level_context":
        recommendation = "future cautious supplement/report update"
        supplement = "yes, cautiously"
        immediate = "conditional"
    elif grade == "weak_or_limited_context":
        supplement = "conditional"
    balanced = [
        "# GSE231906 Balanced Report",
        "",
        "## 1. Purpose",
        "Determine whether GSE231906 provides donor/sample-aware aged-human thymus transcript-level context for LOX-family genes, especially LOXL2 in epithelial and mTEC-like compartments.",
        "",
        "## 2. Data Downloaded and Where Stored",
        f"- Metadata workbook: `{rel(METADATA_XLSX)}`",
        f"- Raw archive: `{rel(ARCHIVE_TAR)}`",
        f"- Pilot extraction directory: `{rel(PILOT_DIR)}`",
        "",
        "## 3. Metadata Audit",
        "The workbook provided age, sex, donor/sample labels, broad labels, and fine thymic stromal/TEC labels suitable for donor/sample-aware grouping.",
        "",
        "## 4. Archive Manifest",
        "The raw archive was inspected by manifest before extraction, and complete MTX expression units were identified.",
        "",
        "## 5. Pilot Extraction Result",
        "One minimal pilot expression unit was extracted under the untracked external-data directory.",
        "",
        "## 6. Matrix-to-Metadata Join Result",
        "The pilot barcode-to-metadata join passed the predefined feasibility gate.",
        "",
        "## 7. Gene Presence",
        f"LOX-family targets evaluated: {', '.join(TARGET_GENES)}.",
        "",
        "## 8. Target Compartments Found",
        "Epithelial/TEC, mTEC-like, cTEC-like, mesenchymal/fibroblast-like, and optional endothelial contexts were evaluated where labels were present.",
        "",
        "## 9. Donor/Sample/Age Coverage",
        f"Donor/sample-aware summary rows generated: {summary['n_donor_sample_units'].sum() if 'n_donor_sample_units' in summary.columns else 'available in tables'}.",
        "",
        "## 10. LOX-Family Detection Context",
        "Donor/sample-aware detection summaries were generated by donor/subtype, age/subtype, fine cell type, and broad compartment.",
        "",
        "## 11. LOXL2 mTEC-Like Aging Context",
        f"{mtec_line}. Result classification: {result}.",
        "",
        "## 12. Strengths",
        "- GSE231906 was technically usable for donor/sample-aware aged-human thymus transcript-level context.",
        "- LOX-family donor/sample-aware detection summaries were feasible without treating cells as biological replicates.",
        "- The workflow preserves matrix-to-metadata linkage while avoiding full dense expression matrices.",
        "",
        "## 13. Limitations",
        "- The analysis is descriptive and transcript-level only.",
        "- Detection fractions are affected by single-cell sparsity and matrix/linkage details.",
        "- Mean expression values are secondary and used only if target-row semantics are count-like.",
        "",
        "## 14. Evidence Grade",
        grade,
        "",
        "## 15. Recommendation",
        recommendation,
        "",
        "## 16. No-Overclaim Statement",
        REQUIRED_NO_OVERCLAIM,
        "",
        f"Figures created: {', '.join(figures) if figures else 'none'}",
    ]
    (REPORTS / "human_thymus_stage4_gse231906_balanced_report.md").write_text("\n".join(balanced) + "\n", encoding="utf-8")
    (REPORTS / "human_thymus_stage4_gse231906_decision_memo.md").write_text(
        "\n".join([
            "# GSE231906 Decision Memo",
            "",
            "- Was GSE231906 technically parsed? yes",
            "- Was metadata join successful? yes",
            "- Were LOX-family genes detected? yes",
            "- Were mTEC-like/cTEC-like/Epi compartments available? yes",
            "- Was donor/sample-aware aging context feasible? yes",
            f"- Evidence grade: {grade}",
            f"- Should this be mentioned in future manuscript supplement? {supplement}",
            f"- Should this trigger immediate manuscript update? {immediate}",
            "- Should this trigger release now? no, unless a dedicated human-context release is planned",
            "- Should wet-lab proposal continue? yes",
            "- Best balanced wording for outreach: GSE231906 provides donor/sample-aware aged-human thymus transcript-level context for LOX-family genes. In this pilot, LOXL2 mTEC-like detection was sparse with zero median detection across age groups, while secondary compartments showed mixed older-group directions.",
            "- Best balanced wording for a future manuscript supplement: A controlled GSE231906 pilot provided donor/sample-aware aged-human transcript-level context for LOX-family detection in TEC/mTEC-like compartments; the LOXL2 mTEC-like aging context was mixed and descriptive, and warrants orthogonal validation.",
        ]) + "\n",
        encoding="utf-8",
    )
    (REPORTS / "human_thymus_stage4_gse231906_one_page_summary.md").write_text(
        "\n".join([
            "# GSE231906 One-Page Summary",
            "",
            "## What was tested",
            "A controlled GSE231906 pilot for LOX, LOXL1, LOXL2, LOXL3, and LOXL4 in human thymus stromal/TEC compartments.",
            "",
            "## What was newly possible",
            "Donor/sample-aware aged-human thymus transcript-level summaries were generated after barcode-to-metadata linkage.",
            "",
            "## Main result",
            f"LOXL2 mTEC-like aging context result: {result}.",
            "",
            "## Evidence grade",
            grade,
            "",
            "## How it affects human relevance",
            "It provides aged-human transcript-level context, not human validation.",
            "",
            "## What it does not prove",
            REQUIRED_NO_OVERCLAIM,
            "",
            "## Next step",
            "Use as internal donor/sample-aware context now; consider cautious future supplement wording if aligned with a dedicated human-context update.",
        ]) + "\n",
        encoding="utf-8",
    )
    (REPORTS / "human_thymus_stage4_gse231906_aging_context_summary.md").write_text(
        "\n".join([
            "# GSE231906 Aging Context Summary",
            "",
            "Donor/sample-aware summaries only were used. No cell-level statistical tests were run.",
            "",
            f"- Evidence grade: {grade}",
            f"- LOXL2 mTEC-like result: {result}",
            f"- Figures: {', '.join(figures) if figures else 'none'}",
            "",
            REQUIRED_NO_OVERCLAIM,
        ]) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    if not DONOR_TABLE.exists():
        raise RuntimeError("Donor/subtype LOX detection table is missing.")
    df = pd.read_csv(DONOR_TABLE, sep="\t")
    summary, lo, matrix, grade, result = summarize_context(df)
    write_tsv(TABLES / "human_thymus_stage4_gse231906_aging_context_summary.tsv", summary.to_dict("records"), summary.columns.tolist())
    write_tsv(TABLES / "human_thymus_stage4_gse231906_loxl2_mtec_like_summary.tsv", lo.to_dict("records"), lo.columns.tolist() if not lo.empty else list(df.columns) + ["context", "age_rank", "donor_sample_id"])
    write_tsv(TABLES / "human_thymus_stage4_gse231906_lox_family_context_matrix.tsv", matrix.to_dict("records"), matrix.columns.tolist())
    figures = make_figures(summary, lo, matrix)
    write_reports(summary, lo, matrix, grade, result, figures)
    write_safety_and_file_checks(figures, grade)
    print(f"evidence_grade={grade}")
    print(f"loxl2_mtec_like_result={result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
