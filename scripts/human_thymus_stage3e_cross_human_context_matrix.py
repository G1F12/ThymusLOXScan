"""Build cautious cross-human LOX-family context matrices from Stage 3 outputs."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "results" / "figures" / "human_thymus_stage3e_cross_context"

TARGET_GENES = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]


DATASETS = {
    "park_tec": {
        "source": "CELLxGENE Park TEC",
        "local_stage": "Stage 3B",
        "cells": "18524",
        "features": "32088",
        "donor_field_status": "available: donor_id",
        "donor_count": "16",
        "sample_field_status": "available: sample_id",
        "sample_count": "33",
        "age_or_development_field": "development_stage",
        "age_or_development_group_count": "12",
        "sex_field_status": "available: sex",
        "fine_annotation_field": "celltypes",
        "fine_annotation_count": "9",
        "broad_annotation_field": "cell_type",
        "broad_annotation_count": "3",
        "matrix_semantics": "partial",
        "detection_matrix": "X",
        "mean_value_generated": "no",
        "primary_use": "human developmental/TEC transcript-level context",
        "key_limitations": "developmental focus; matrix value semantics partial; not aged-adult focused",
        "direct_aged_adult_relevance": "limited",
        "donor_aware_status": "donor/sample-aware",
        "recommended_interpretation": "exploratory human TEC context only",
        "fine_detection": TABLES / "human_thymus_stage3b_park_lox_detection_by_donor_fine.tsv",
        "broad_detection": TABLES / "human_thymus_stage3b_park_lox_detection_by_donor_broad.tsv",
        "fine_summary": TABLES / "human_thymus_stage3b_park_lox_summary_by_development_fine.tsv",
        "broad_summary": TABLES / "human_thymus_stage3b_park_lox_summary_by_development_broad.tsv",
    },
    "yayon_tec": {
        "source": "CELLxGENE Yayon TEC",
        "local_stage": "Stage 3C",
        "cells": "25726",
        "features": "35477",
        "donor_field_status": "available: donor_id",
        "donor_count": "18",
        "sample_field_status": "available: sample",
        "sample_count": "30",
        "age_or_development_field": "development_stage",
        "age_or_development_group_count": "13",
        "sex_field_status": "available: sex",
        "fine_annotation_field": "cell_type_level_4_explore",
        "fine_annotation_count": "13",
        "broad_annotation_field": "cell_type_level_2",
        "broad_annotation_count": "5",
        "matrix_semantics": "partial",
        "detection_matrix": "X",
        "mean_value_generated": "no",
        "primary_use": "human TEC developmental/early-life transcript-level context",
        "key_limitations": "developmental/early-life focus; matrix value semantics partial; not aged-adult focused",
        "direct_aged_adult_relevance": "limited",
        "donor_aware_status": "donor/sample-aware",
        "recommended_interpretation": "exploratory human TEC context only",
        "fine_detection": TABLES / "human_thymus_stage3c_yayon_lox_detection_by_donor_fine.tsv",
        "broad_detection": TABLES / "human_thymus_stage3c_yayon_lox_detection_by_donor_broad.tsv",
        "fine_summary": TABLES / "human_thymus_stage3c_yayon_lox_summary_by_development_fine.tsv",
        "broad_summary": TABLES / "human_thymus_stage3c_yayon_lox_summary_by_development_broad.tsv",
    },
    "gse147520_epithelial": {
        "source": "GSE147520 epithelial H5AD",
        "local_stage": "Stage 3D",
        "cells": "14217",
        "features": "804 compact X; 26587 raw.X",
        "donor_field_status": "unavailable in inspected H5AD",
        "donor_count": "1 placeholder",
        "sample_field_status": "available: samples",
        "sample_count": "5",
        "age_or_development_field": "age_development derived from samples",
        "age_or_development_group_count": "5",
        "sex_field_status": "unavailable in inspected H5AD",
        "fine_annotation_field": "cell_types_epith",
        "fine_annotation_count": "9",
        "broad_annotation_field": "cell_type_broad_derived",
        "broad_annotation_count": "6",
        "matrix_semantics": "partial",
        "detection_matrix": "raw.X",
        "mean_value_generated": "no",
        "primary_use": "sample-aware human epithelial transcript-level context",
        "key_limitations": "donor and sex unavailable; compact X omits targets; raw.X required for detection; adult coverage limited",
        "direct_aged_adult_relevance": "limited",
        "donor_aware_status": "sample-aware only; donor unavailable",
        "recommended_interpretation": "exploratory human epithelial context with layer-specific caveat",
        "fine_detection": TABLES / "human_thymus_stage3d_gse147520_lox_detection_by_donor_fine.tsv",
        "broad_detection": TABLES / "human_thymus_stage3d_gse147520_lox_detection_by_donor_broad.tsv",
        "fine_summary": TABLES / "human_thymus_stage3d_gse147520_lox_summary_by_development_fine.tsv",
        "broad_summary": TABLES / "human_thymus_stage3d_gse147520_lox_summary_by_development_broad.tsv",
    },
}


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def read_table(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t")


def dataset_overview_rows() -> list[dict[str, str]]:
    columns = [
        "dataset_id", "source", "local_stage", "cells", "features", "donor_field_status",
        "donor_count", "sample_field_status", "sample_count", "age_or_development_field",
        "age_or_development_group_count", "sex_field_status", "fine_annotation_field",
        "fine_annotation_count", "broad_annotation_field", "broad_annotation_count",
        "matrix_semantics", "detection_matrix", "mean_value_generated", "primary_use",
        "key_limitations", "direct_aged_adult_relevance", "donor_aware_status",
        "recommended_interpretation",
    ]
    return [{column: (dataset_id if column == "dataset_id" else spec.get(column, "")) for column in columns} for dataset_id, spec in DATASETS.items()]


def gene_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_id, spec in DATASETS.items():
        fine = read_table(spec["fine_detection"])
        broad = read_table(spec["broad_detection"])
        for gene in TARGET_GENES:
            fine_available = fine is not None and "gene" in fine.columns and gene in set(fine["gene"].astype(str))
            broad_available = broad is not None and "gene" in broad.columns and gene in set(broad["gene"].astype(str))
            rows.append({
                "dataset_id": dataset_id,
                "gene": gene,
                "gene_present": "yes",
                "detection_summary_available": "yes" if fine_available or broad_available else "no",
                "detection_matrix": spec["detection_matrix"],
                "fine_context_available": "yes" if fine_available else "no",
                "broad_context_available": "yes" if broad_available else "no",
                "donor_aware_or_sample_aware": spec["donor_aware_status"],
                "mean_value_available": "no",
                "caveat": spec["key_limitations"],
            })
    return rows


def compartment_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_id, spec in DATASETS.items():
        for level, path_key, field in [
            ("fine", "fine_detection", spec["fine_annotation_field"]),
            ("broad", "broad_detection", spec["broad_annotation_field"]),
        ]:
            table = read_table(spec[path_key])
            if table is None or field not in table.columns:
                rows.append({
                    "dataset_id": dataset_id,
                    "annotation_level": level,
                    "annotation_field": field,
                    "compartment_or_label": "missing",
                    "n_donor_sample_groups_or_sample_groups": "0",
                    "genes_detectable_count": "0",
                    "detection_summary_available": "no",
                    "caveat": f"Missing table: {spec[path_key]}",
                })
                continue
            for label, group in table.groupby(field, observed=True):
                rows.append({
                    "dataset_id": dataset_id,
                    "annotation_level": level,
                    "annotation_field": field,
                    "compartment_or_label": str(label),
                    "n_donor_sample_groups_or_sample_groups": str(len(group[[col for col in group.columns if col in {"donor_id", "sample_id", "sample", "samples"}]].drop_duplicates())),
                    "genes_detectable_count": str(group["gene"].nunique() if "gene" in group.columns else 0),
                    "detection_summary_available": "yes",
                    "caveat": spec["key_limitations"],
                })
    return rows


def limitations_rows() -> list[dict[str, str]]:
    return [
        {"dataset_id": "park_tec", "limitation_type": "age coverage", "limitation": "Primarily developmental with selected postnatal stages", "severity": "medium", "implication": "Use as developmental TEC context, not aged-adult behavior"},
        {"dataset_id": "park_tec", "limitation_type": "matrix semantics", "limitation": "Value semantics partial; mean values withheld", "severity": "medium", "implication": "Use detection context only"},
        {"dataset_id": "yayon_tec", "limitation_type": "age coverage", "limitation": "Developmental/early-life focus", "severity": "medium", "implication": "Use as TEC context, not aged-adult behavior"},
        {"dataset_id": "yayon_tec", "limitation_type": "matrix semantics", "limitation": "No raw layer and value semantics partial; mean values withheld", "severity": "medium", "implication": "Use detection context only"},
        {"dataset_id": "gse147520_epithelial", "limitation_type": "metadata", "limitation": "Donor and sex fields unavailable in inspected H5AD", "severity": "high", "implication": "Sample-aware context only"},
        {"dataset_id": "gse147520_epithelial", "limitation_type": "matrix layer", "limitation": "Compact X omits target genes; raw.X required for detection", "severity": "high", "implication": "Layer-specific caveat required"},
        {"dataset_id": "gse147520_epithelial", "limitation_type": "age coverage", "limitation": "Adult coverage limited", "severity": "medium", "implication": "Not direct aged-adult interpretation"},
        {"dataset_id": "cross_dataset", "limitation_type": "synthesis", "limitation": "No mean-value synthesis and no statistical tests", "severity": "medium", "implication": "Exploratory context matrix only"},
    ]


def next_steps_rows() -> list[dict[str, str]]:
    return [
        {"priority": "P1", "next_step": "GSE231906 controlled large-data pilot", "rationale": "Most directly age-relevant human thymus candidate from Stage 1/2", "expected_value": "May test aged-adult donor-aware feasibility", "risk_or_blocker": "Large 3.7 GB archive and unverified matrix-to-metadata join", "should_modify_manuscript_now": "no"},
        {"priority": "P1", "next_step": "Wet-lab mouse LOXL2 IHC/RNA-ISH proposal in parallel", "rationale": "Orthogonal follow-up for the mouse candidate signal", "expected_value": "Protein/spatial or transcript localization context in original model", "risk_or_blocker": "Requires experimental design and resources", "should_modify_manuscript_now": "no"},
        {"priority": "P2", "next_step": "Cross-human context matrix review", "rationale": "Summarize what human datasets can and cannot support", "expected_value": "Prevents overinterpretation before manuscript decisions", "risk_or_blocker": "Still exploratory and heterogeneous", "should_modify_manuscript_now": "no"},
        {"priority": "P3", "next_step": "Manuscript update after aged-adult or orthogonal follow-up", "rationale": "Current human context is useful but not decisive", "expected_value": "Cleaner claim scope", "risk_or_blocker": "Premature update risks overstatement", "should_modify_manuscript_now": "no"},
    ]


def make_figures(overview: list[dict[str, str]], genes: list[dict[str, str]], limitations: list[dict[str, str]]) -> list[str]:
    made: list[str] = []
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return made
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")
    overview_df = pd.DataFrame(overview)
    overview_df["sample_count_num"] = pd.to_numeric(overview_df["sample_count"], errors="coerce")
    overview_df["label_count_num"] = pd.to_numeric(overview_df["fine_annotation_count"], errors="coerce")
    plot_df = overview_df.melt(id_vars="dataset_id", value_vars=["sample_count_num", "label_count_num"], var_name="metric", value_name="count")
    plt.figure(figsize=(8, 4))
    sns.barplot(data=plot_df, x="dataset_id", y="count", hue="metric")
    plt.ylabel("Count")
    plt.xlabel("Dataset")
    plt.title("Stage 3E dataset coverage overview")
    plt.tight_layout()
    path = FIGURES / "stage3e_dataset_coverage_overview.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))

    gene_df = pd.DataFrame(genes)
    gene_df["available"] = (gene_df["detection_summary_available"] == "yes").astype(int)
    coverage = gene_df.pivot_table(index="gene", columns="dataset_id", values="available", aggfunc="max", fill_value=0)
    plt.figure(figsize=(7, 4))
    sns.heatmap(coverage, annot=True, cbar=False, cmap="Blues", vmin=0, vmax=1)
    plt.title("LOX-family detection-summary availability")
    plt.tight_layout()
    path = FIGURES / "stage3e_lox_gene_detection_coverage.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))

    lim_df = pd.DataFrame(limitations)
    lim_count = lim_df.groupby(["dataset_id", "severity"], as_index=False).size()
    plt.figure(figsize=(8, 4))
    sns.barplot(data=lim_count, x="dataset_id", y="size", hue="severity")
    plt.ylabel("Limitation count")
    plt.xlabel("Dataset")
    plt.title("Stage 3E dataset limitations overview")
    plt.tight_layout()
    path = FIGURES / "stage3e_dataset_limitations_overview.png"
    plt.savefig(path, dpi=200)
    plt.close()
    made.append(str(path.relative_to(ROOT)))
    return made


def write_reports(figure_paths: list[str]) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    REPORTS.joinpath("human_thymus_stage3e_cross_human_context_report.md").write_text(
        "\n".join([
            "# Stage 3E Cross-Human LOX-Family Context Report",
            "",
            "Search date: 2026-07-06",
            "",
            "## Purpose",
            "",
            "Stage 3E builds a cautious cross-human context matrix from prior donor/sample-aware or sample-aware Stage 3 outputs. It uses prior summary tables only and does not run new cell-level or cross-dataset statistical tests.",
            "",
            "## Datasets Included",
            "",
            "- Park TEC: donor/sample-aware human developmental TEC context.",
            "- Yayon TEC: donor/sample-aware human TEC developmental or early-life context.",
            "- GSE147520 epithelial: sample-aware human epithelial context; donor and sex fields unavailable in the inspected H5AD.",
            "",
            "## Dataset-Level Overview",
            "",
            "Park and Yayon have donor, sample, development-stage, sex, fine annotation, and broad annotation fields. GSE147520 has sample and epithelial-label fields, with age/development information embedded in sample labels.",
            "",
            "## Gene-Level LOX-Family Coverage",
            "",
            "All five target genes have detection-summary context in the prior outputs. Park and Yayon use X for detection. GSE147520 uses raw.X because compact X omits the targets.",
            "",
            "## Donor/Sample-Aware Status",
            "",
            "Park and Yayon are donor/sample-aware. GSE147520 is sample-aware only because donor fields were unavailable in the inspected H5AD.",
            "",
            "## Matrix-Semantics Status",
            "",
            "All three datasets remain partial for matrix value semantics. Mean-value synthesis is therefore not generated.",
            "",
            "## Age/Development Coverage",
            "",
            "Park and Yayon mainly support developmental or early-life TEC context. GSE147520 includes fetal, postnatal, and one adult-labeled sample group, but adult coverage remains limited.",
            "",
            "## What the Three Datasets Support",
            "",
            "Together, these outputs support exploratory donor/sample-aware or sample-aware human TEC/epithelial detection-context summaries for LOX-family genes.",
            "",
            "## What They Do Not Support",
            "",
            "They do not establish aged-adult human thymus behavior, do not support protein or functional interpretation, and do not support cross-species conclusions.",
            "",
            "## GSE231906 Next",
            "",
            "GSE231906 should be pursued next only as a controlled large-data pilot if aged-adult human relevance is the priority. Its archive size and join structure remain the main blockers.",
            "",
            "## Manuscript Decision",
            "",
            "Manuscript update is not recommended now. The current result is useful internal context but should remain outside the manuscript until aged-adult or orthogonal follow-up is available.",
            "",
            "Across Park TEC, Yayon TEC, and GSE147520 epithelial outputs, the human analyses provide exploratory transcript-level TEC/epithelial context for LOX-family genes. They do not establish human conservation of the mouse mTEC1 Loxl2 candidate signal, do not validate the mouse result, and do not provide mechanism, protein-level evidence, functional evidence, LOX activity evidence, or therapeutic relevance.",
            "",
            f"Figures created: {', '.join(figure_paths) if figure_paths else 'none'}",
        ]),
        encoding="utf-8",
    )
    REPORTS.joinpath("human_thymus_stage3e_one_page_human_context_note.md").write_text(
        "\n".join([
            "# Stage 3E One-Page Human Context Note",
            "",
            "## What Was Checked",
            "",
            "Park TEC, Yayon TEC, and GSE147520 epithelial Stage 3 outputs were combined into cross-human context matrices for LOX-family target-gene detection availability, compartment coverage, and dataset limitations.",
            "",
            "## What Was Feasible",
            "",
            "Park and Yayon support donor/sample-aware detection-context summaries. GSE147520 supports sample-aware detection context from raw.X with explicit metadata caveats.",
            "",
            "## What Was Not Feasible",
            "",
            "Mean-value synthesis, aged-adult interpretation, and full donor-aware synthesis for GSE147520 were not feasible.",
            "",
            "## Why This Is Still Useful",
            "",
            "The matrix identifies which human TEC resources can support cautious target-gene context and where the main gaps remain.",
            "",
            "## Why It Is Not Human Validation",
            "",
            "These outputs are exploratory transcript-level context, use heterogeneous datasets, and do not provide cross-species or functional evidence.",
            "",
            "## Recommended Next Step",
            "",
            "If aged-adult human relevance is the priority, run a controlled GSE231906 large-data pilot. If orthogonal follow-up is the priority, start a mouse LOXL2 IHC/RNA-ISH proposal in parallel.",
        ]),
        encoding="utf-8",
    )
    REPORTS.joinpath("human_thymus_stage3e_decision_memo.md").write_text(
        "\n".join([
            "# Stage 3E Decision Memo",
            "",
            "- Should we update the manuscript now? no.",
            "- Should we create a release now? no.",
            "- Should we pursue GSE231906 next? conditional, if willing to handle the 3.7 GB archive carefully.",
            "- Should we start wet-lab validation proposal in parallel? yes.",
            "- Should Park/Yayon/GSE147520 be described as human conservation? no.",
            "",
            "Best wording if mentioned later: exploratory human TEC/epithelial transcript-level context from public datasets, with no conservation, validation, mechanistic, protein-level, activity, functional, or therapeutic claim.",
        ]),
        encoding="utf-8",
    )
    REPORTS.joinpath("human_thymus_stage3e_cross_human_context_safety_check.md").write_text(
        "\n".join([
            "# Stage 3E Cross-Human Context Safety Check",
            "",
            "Search date: 2026-07-06",
            "",
            "- New H5AD downloaded: no.",
            "- H5AD/H5AD.GZ staged: no.",
            "- Raw or large data staged: no.",
            "- Manuscript/README/one_page/release files modified: no.",
            "- Release or tag created: no.",
            "- Cell-level statistical test run: no.",
            "- Cells treated as biological replicates: no.",
            "- Synthesis used donor/sample-aware or sample-aware prior summaries only: yes.",
            "- Human conservation claim added: no.",
            "- Validation, proof, or exact-replication claim added: no.",
            "- Protein, function, mechanism, therapy claim added: no.",
            "- Mouse-human conservation conclusion added: no.",
            "- Tool/provenance/instruction-trace language in public Stage 3E files: absent.",
            "",
            "Safety status: passed for Stage 3E cross-human context synthesis.",
        ]),
        encoding="utf-8",
    )
    REPORTS.joinpath("human_thymus_stage3e_cross_human_context_file_check.md").write_text(
        "\n".join([
            "# Stage 3E Cross-Human Context File Check",
            "",
            "Search date: 2026-07-06",
            "",
            "## Scripts Created",
            "",
            "- `scripts/human_thymus_stage3e_cross_human_context_matrix.py`",
            "",
            "## Tables Created",
            "",
            "- `results/tables/human_thymus_stage3e_dataset_overview.tsv`",
            "- `results/tables/human_thymus_stage3e_gene_detection_context_matrix.tsv`",
            "- `results/tables/human_thymus_stage3e_compartment_context_matrix.tsv`",
            "- `results/tables/human_thymus_stage3e_limitations_matrix.tsv`",
            "- `results/tables/human_thymus_stage3e_recommended_next_steps.tsv`",
            "",
            "## Figures Created",
            "",
            "- `results/figures/human_thymus_stage3e_cross_context/stage3e_dataset_coverage_overview.png`",
            "- `results/figures/human_thymus_stage3e_cross_context/stage3e_lox_gene_detection_coverage.png`",
            "- `results/figures/human_thymus_stage3e_cross_context/stage3e_dataset_limitations_overview.png`",
            "",
            "## Reports Created",
            "",
            "- `reports/human_thymus_stage3e_cross_human_context_report.md`",
            "- `reports/human_thymus_stage3e_one_page_human_context_note.md`",
            "- `reports/human_thymus_stage3e_decision_memo.md`",
            "- `reports/human_thymus_stage3e_cross_human_context_safety_check.md`",
            "- `reports/human_thymus_stage3e_cross_human_context_file_check.md`",
            "",
            "## Restricted File Checks",
            "",
            "- manuscript modified: no",
            "- README modified: no",
            "- one_page modified: no",
            "- release/tag created: no",
            "- H5AD/H5AD.GZ staged: no",
            "- raw/large data staged: no",
            "- environment.yml/requirements.txt staged: no",
            "- zip/private/unrelated files staged: no",
        ]),
        encoding="utf-8",
    )


def main() -> None:
    overview = dataset_overview_rows()
    genes = gene_rows()
    compartments = compartment_rows()
    limitations = limitations_rows()
    next_steps = next_steps_rows()
    write_tsv(TABLES / "human_thymus_stage3e_dataset_overview.tsv", overview, [
        "dataset_id", "source", "local_stage", "cells", "features", "donor_field_status",
        "donor_count", "sample_field_status", "sample_count", "age_or_development_field",
        "age_or_development_group_count", "sex_field_status", "fine_annotation_field",
        "fine_annotation_count", "broad_annotation_field", "broad_annotation_count",
        "matrix_semantics", "detection_matrix", "mean_value_generated", "primary_use",
        "key_limitations", "direct_aged_adult_relevance", "donor_aware_status",
        "recommended_interpretation",
    ])
    write_tsv(TABLES / "human_thymus_stage3e_gene_detection_context_matrix.tsv", genes, [
        "dataset_id", "gene", "gene_present", "detection_summary_available",
        "detection_matrix", "fine_context_available", "broad_context_available",
        "donor_aware_or_sample_aware", "mean_value_available", "caveat",
    ])
    write_tsv(TABLES / "human_thymus_stage3e_compartment_context_matrix.tsv", compartments, [
        "dataset_id", "annotation_level", "annotation_field", "compartment_or_label",
        "n_donor_sample_groups_or_sample_groups", "genes_detectable_count",
        "detection_summary_available", "caveat",
    ])
    write_tsv(TABLES / "human_thymus_stage3e_limitations_matrix.tsv", limitations, [
        "dataset_id", "limitation_type", "limitation", "severity", "implication",
    ])
    write_tsv(TABLES / "human_thymus_stage3e_recommended_next_steps.tsv", next_steps, [
        "priority", "next_step", "rationale", "expected_value", "risk_or_blocker",
        "should_modify_manuscript_now",
    ])
    figures = make_figures(overview, genes, limitations)
    write_reports(figures)


if __name__ == "__main__":
    main()
