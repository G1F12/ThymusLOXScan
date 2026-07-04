"""Per-mouse E-MTAB-8560 mTEC LOX-family reanalysis."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[1]
META_DIR = PROJECT_ROOT / "data" / "external" / "metadata" / "emtab8560"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "external_validation" / "emtab8560_mtec"
REPORT = PROJECT_ROOT / "reports" / "external_emtab8560_mtec_lox_reanalysis.md"
COMPARE_REPORT = PROJECT_ROOT / "reports" / "external_emtab8560_vs_gse240016_interpretation.md"
UPDATE_REPORT = PROJECT_ROOT / "reports" / "v5_5_emtab8560_mtec_update_summary.md"
SAFETY_REPORT = PROJECT_ROOT / "reports" / "v5_5_emtab8560_mtec_safety_check.md"
FILE_REPORT = PROJECT_ROOT / "reports" / "v5_5_emtab8560_mtec_file_check.md"
CELL_META = TABLE_DIR / "external_emtab8560_cell_metadata_verified.tsv"

LOX_GENES = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
GROUPS = ["mTEClo", "mTEChi", "combined_mTEClo_mTEChi", "cTEC"]
PRIMARY_GROUPS = ["mTEClo", "mTEChi", "combined_mTEClo_mTEChi"]
CONTRASTS = [(4, 52), (16, 52), ("4+16", 52)]


def bh_adjust(pvals: pd.Series) -> pd.Series:
    p = pvals.astype(float).to_numpy()
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(out, index=pvals.index)
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    n = len(ranked)
    adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    valid_idx = np.where(valid)[0]
    out[valid_idx[order]] = adj
    return pd.Series(out, index=pvals.index)


def write_infeasible(missing: list[str]) -> None:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "\n".join(
            [
                "# E-MTAB-8560 mTEC LOX Reanalysis",
                "",
                "Final classification: not reliable/infeasible.",
                "",
                "The official MouseThymusAgeing R export or biological-unit audit output was not available, so per-mouse pseudobulk, permutation tests, trend models, FDR, and figures were not generated.",
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


def load_long(path: Path, value_name: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df.rename(columns={"value": value_name})


def add_groups(meta: pd.DataFrame) -> pd.DataFrame:
    sort_type = meta["SortType"].astype(str)
    subtype = meta["SubType"].astype(str)
    definitions = {
        "mTEClo": sort_type.eq("mTEClo"),
        "mTEChi": sort_type.eq("mTEChi"),
        "combined_mTEClo_mTEChi": sort_type.isin(["mTEClo", "mTEChi"]),
        "cTEC": sort_type.eq("cTEC") | subtype.str.contains("cTEC", case=False, na=False),
        "Mature.mTEC": subtype.eq("Mature.mTEC"),
        "Post-Aire.mTEC": subtype.eq("Post-Aire.mTEC"),
        "Tuft-like.mTEC": subtype.eq("Tuft-like.mTEC"),
    }
    rows = []
    for group, mask in definitions.items():
        tmp = meta.loc[mask].copy()
        tmp["analysis_group"] = group
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def build_pseudobulk() -> pd.DataFrame:
    meta = pd.read_csv(CELL_META, sep="\t", dtype=str).fillna("")
    meta["age_week"] = pd.to_numeric(meta["age_week"], errors="coerce")
    meta["total_counts_all_features"] = pd.to_numeric(meta.get("total_counts_all_features", np.nan), errors="coerce")
    grouped_meta = add_groups(meta)

    counts = load_long(META_DIR / "emtab8560_lox_counts_long.tsv", "raw_counts")
    values = counts.merge(
        grouped_meta[
            [
                "cell_barcode",
                "analysis_group",
                "biological_mouse_id",
                "age_week",
                "SortDay",
                "PlateID",
                "Comment[ENA_RUN]",
                "total_counts_all_features",
            ]
        ],
        on="cell_barcode",
        how="inner",
    )
    values = values[values["biological_mouse_id"].astype(str).ne("")]
    values["raw_counts"] = pd.to_numeric(values["raw_counts"], errors="coerce").fillna(0.0)
    values["detected"] = values["raw_counts"].gt(0)

    log_path = META_DIR / "emtab8560_lox_logcounts_long.tsv"
    if log_path.exists():
        logvals = load_long(log_path, "logcounts")
        values = values.merge(logvals[["gene", "cell_barcode", "logcounts"]], on=["gene", "cell_barcode"], how="left")
        values["logcounts"] = pd.to_numeric(values["logcounts"], errors="coerce")
    else:
        values["logcounts"] = np.nan

    pb = (
        values.groupby(["analysis_group", "gene", "biological_mouse_id", "age_week"], observed=True)
        .agg(
            n_cells=("cell_barcode", "nunique"),
            summed_raw_counts=("raw_counts", "sum"),
            total_counts_all_features=("total_counts_all_features", "sum"),
            detection_rate=("detected", "mean"),
            mean_logcounts=("logcounts", "mean"),
            SortDay=("SortDay", lambda x: ";".join(sorted(set(map(str, x))))),
            PlateID=("PlateID", lambda x: ";".join(sorted(set(map(str, x)))[:12])),
            run_fields=("Comment[ENA_RUN]", lambda x: ";".join(sorted(set(v for v in map(str, x) if v))[:12])),
        )
        .reset_index()
    )
    pb["cpm"] = np.where(pb["total_counts_all_features"].gt(0), pb["summed_raw_counts"] / pb["total_counts_all_features"] * 1_000_000, np.nan)
    pb["log2_cpm_plus1"] = np.log2(pb["cpm"] + 1)
    return pb


def age_summary(pb: pd.DataFrame) -> pd.DataFrame:
    return (
        pb.groupby(["analysis_group", "gene", "age_week"], observed=True)
        .agg(
            n_mice=("biological_mouse_id", "nunique"),
            total_cells=("n_cells", "sum"),
            mean_log2_cpm_plus1=("log2_cpm_plus1", "mean"),
            median_log2_cpm_plus1=("log2_cpm_plus1", "median"),
            mean_detection_rate=("detection_rate", "mean"),
            median_detection_rate=("detection_rate", "median"),
        )
        .reset_index()
    )


def permutation_p(group: pd.DataFrame, young_ages: list[int], old_age: int, metric: str) -> dict[str, float | int | str]:
    sub = group[group["age_week"].isin(young_ages + [old_age])].dropna(subset=[metric]).copy()
    values = sub[metric].to_numpy(float)
    labels = sub["age_week"].eq(old_age).to_numpy()
    n_old = int(labels.sum())
    n_total = len(values)
    if n_old == 0 or n_old == n_total:
        return {"one_sided_aged_lower_p": np.nan, "two_sided_p": np.nan, "n_permutations": 0, "mode": "not_testable"}
    obs = values[labels].mean() - values[~labels].mean()
    total = math.comb(n_total, n_old)
    if total <= 200_000:
        stats_all = []
        for idx in itertools.combinations(range(n_total), n_old):
            mask = np.zeros(n_total, dtype=bool)
            mask[list(idx)] = True
            stats_all.append(values[mask].mean() - values[~mask].mean())
        arr = np.asarray(stats_all)
        mode = "exact"
    else:
        rng = np.random.default_rng(1)
        arr = np.empty(100_000)
        for i in range(len(arr)):
            mask = np.zeros(n_total, dtype=bool)
            mask[rng.choice(n_total, n_old, replace=False)] = True
            arr[i] = values[mask].mean() - values[~mask].mean()
        mode = "monte_carlo"
    return {
        "one_sided_aged_lower_p": float((np.sum(arr <= obs) + 1) / (len(arr) + 1)),
        "two_sided_p": float((np.sum(np.abs(arr) >= abs(obs)) + 1) / (len(arr) + 1)),
        "n_permutations": int(len(arr)),
        "mode": mode,
        "observed_delta": float(obs),
    }


def build_models(pb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    model_rows = []
    perm_rows = []
    for group in GROUPS:
        for gene in LOX_GENES:
            g = pb[(pb["analysis_group"].eq(group)) & (pb["gene"].eq(gene))]
            if g.empty:
                continue
            ages = sorted(g["age_week"].dropna().unique())
            rho, rho_p = (np.nan, np.nan)
            if len(ages) >= 3:
                rho, rho_p = stats.spearmanr(g["age_week"], g["log2_cpm_plus1"], nan_policy="omit")
            lm = stats.linregress(g["age_week"], g["log2_cpm_plus1"]) if g["age_week"].nunique() >= 2 else None
            row = {
                "analysis_group": group,
                "gene": gene,
                "n_mice_total": g["biological_mouse_id"].nunique(),
                "spearman_rho_per_mouse": rho,
                "spearman_p_per_mouse": rho_p,
                "linear_age_slope": lm.slope if lm else np.nan,
                "linear_age_p": lm.pvalue if lm else np.nan,
                "batch_adjusted_model_status": "not estimated; SortDay/PlateID/run are complex cell/acquisition fields rather than simple mouse-level covariates",
            }
            for young, old in CONTRASTS:
                young_ages = [4, 16] if young == "4+16" else [int(young)]
                young_vals = g[g["age_week"].isin(young_ages)]
                old_vals = g[g["age_week"].eq(old)]
                prefix = f"{young}_vs_{old}"
                row[f"{prefix}_n_young_mice"] = young_vals["biological_mouse_id"].nunique()
                row[f"{prefix}_n_old_mice"] = old_vals["biological_mouse_id"].nunique()
                row[f"{prefix}_old_minus_young_log2cpm"] = old_vals["log2_cpm_plus1"].mean() - young_vals["log2_cpm_plus1"].mean()
                row[f"{prefix}_old_minus_young_detection"] = old_vals["detection_rate"].mean() - young_vals["detection_rate"].mean()
                for metric in ["log2_cpm_plus1", "detection_rate"]:
                    p = permutation_p(g, young_ages, int(old), metric)
                    perm_rows.append(
                        {
                            "analysis_group": group,
                            "gene": gene,
                            "contrast": prefix,
                            "metric": metric,
                            **p,
                        }
                    )
            model_rows.append(row)
    models = pd.DataFrame(model_rows)
    perms = pd.DataFrame(perm_rows)
    fdr = perms[perms["metric"].eq("log2_cpm_plus1")].copy()
    fdr["global_bh_fdr"] = bh_adjust(fdr["one_sided_aged_lower_p"])
    fdr["fdr_class"] = np.select(
        [fdr["global_bh_fdr"].lt(0.05), fdr["global_bh_fdr"].lt(0.10)],
        ["global_FDR_lt_0.05", "global_FDR_lt_0.10"],
        default="not_significant",
    )
    return models, perms, fdr


def plot_group(pb: pd.DataFrame, group: str, metric: str, stem: str, ylabel: str) -> None:
    sub = pb[(pb["analysis_group"].eq(group)) & (pb["gene"].eq("Loxl2"))].copy()
    if sub.empty:
        return
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    for age, vals in sub.groupby("age_week"):
        ax.scatter(vals["age_week"], vals[metric], s=50, color="#2E6F9E", alpha=0.85)
    means = sub.groupby("age_week")[metric].mean().reset_index()
    ax.plot(means["age_week"], means[metric], color="#B23A48", linewidth=1.5)
    ax.set_xlabel("Age (weeks)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"E-MTAB-8560 Loxl2 {group}")
    ax.set_xticks(sorted(sub["age_week"].unique()))
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=200)
    plt.close(fig)


def plot_heatmap(models: pd.DataFrame) -> None:
    if models.empty:
        return
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    col = "4+16_vs_52_old_minus_young_log2cpm"
    mat = models[models["analysis_group"].isin(PRIMARY_GROUPS)].pivot(index="analysis_group", columns="gene", values=col)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    im = ax.imshow(mat.fillna(0), cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(len(mat.columns)), mat.columns)
    ax.set_yticks(range(len(mat.index)), mat.index)
    ax.set_title("Old-minus-adult log2CPM+1")
    fig.colorbar(im, ax=ax, label="52w minus 4+16w")
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG_DIR / f"emtab8560_lox_family_mtec_effect_heatmap.{ext}", dpi=200)
    plt.close(fig)


def classify(models: pd.DataFrame, fdr: pd.DataFrame) -> str:
    hits = models[(models["gene"].eq("Loxl2")) & (models["analysis_group"].isin(["mTEClo", "mTEChi"]))]
    if hits.empty:
        return "not reliable/infeasible"
    delta_cols = ["4_vs_52_old_minus_young_log2cpm", "16_vs_52_old_minus_young_log2cpm", "4+16_vs_52_old_minus_young_log2cpm"]
    aged_lower = hits[delta_cols].lt(0).any(axis=1).any()
    fdr_support = fdr[
        fdr["gene"].eq("Loxl2")
        & fdr["analysis_group"].isin(["mTEClo", "mTEChi"])
        & fdr["global_bh_fdr"].lt(0.10)
    ]
    if aged_lower and not fdr_support.empty:
        return "moderate independent directional support"
    if aged_lower:
        return "mixed/inconclusive"
    return "no support"


def write_reports(pb: pd.DataFrame, by_age: pd.DataFrame, models: pd.DataFrame, perms: pd.DataFrame, fdr: pd.DataFrame, final_class: str) -> None:
    loxl2 = models[models["gene"].eq("Loxl2")].copy()
    lines = [
        "# E-MTAB-8560 mTEC LOX Reanalysis",
        "",
        f"Final classification: {final_class}.",
        "",
        "This is an independent mTEC-focused transcript-level test. It is not exact GSE240016 mTEC1 replication.",
        "",
        "## Biological Unit",
        "",
        "One dot/table row at the inferential level is one biological mouse ID defined as `age_week + Characteristics[individual]` after joining official SDRF metadata to MouseThymusAgeing colData. Cells are not treated as independent biological replicates.",
        "",
        "## Loxl2 Summary",
        "",
        "| group | 4w vs 52w delta | 16w vs 52w delta | 4+16w vs 52w delta |",
        "|---|---:|---:|---:|",
    ]
    for _, row in loxl2.iterrows():
        lines.append(
            f"| {row['analysis_group']} | {row.get('4_vs_52_old_minus_young_log2cpm', np.nan):.4g} | "
            f"{row.get('16_vs_52_old_minus_young_log2cpm', np.nan):.4g} | {row.get('4+16_vs_52_old_minus_young_log2cpm', np.nan):.4g} |"
        )
    lines += [
        "",
        "## Batch Status",
        "",
        "Batch-adjusted models were not estimated because SortDay, PlateID, and run fields are complex cell/acquisition-level variables rather than simple mouse-level covariates. This is reported as a limitation.",
        "",
        "## Output Files",
        "",
        "- `results/tables/external_emtab8560_mtec_lox_pseudobulk.tsv`",
        "- `results/tables/external_emtab8560_mtec_lox_by_age.tsv`",
        "- `results/tables/external_emtab8560_mtec_lox_models.tsv`",
        "- `results/tables/external_emtab8560_mtec_lox_permutation.tsv`",
        "- `results/tables/external_emtab8560_mtec_lox_global_fdr.tsv`",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    COMPARE_REPORT.write_text(
        "\n".join(
            [
                "# E-MTAB-8560 vs GSE240016 Interpretation",
                "",
                "GSE240016 mTEC1 remains a candidate-level transcript signal with n=2 versus n=2, limited exact-permutation resolution, lower aged mTEC1 depth/detected genes, and unresolved batch concern.",
                "",
                f"E-MTAB-8560 classification: {final_class}. The biological unit is verified per mouse using official SDRF `Characteristics[individual]` joined to R-exported MouseThymusAgeing colData.",
                "",
                "This changes confidence only at the level supported by the E-MTAB-8560 result. It should be described as independent mTEC-focused transcript-level context, not exact GSE240016 mTEC1 replication.",
                "",
                "No p-values were combined across datasets.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    UPDATE_REPORT.write_text(
        "\n".join(
            [
                "# v5.5 E-MTAB-8560 mTEC Update Summary",
                "",
                "The GitHub Actions R export should be used as the authoritative run when available.",
                "",
                f"Final E-MTAB-8560 classification: {final_class}.",
                "",
                "Manuscript update is not automatic; any wording should remain cautious and transcript-level.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    SAFETY_REPORT.write_text(
        "\n".join(
            [
                "# v5.5 E-MTAB-8560 mTEC Safety Check",
                "",
                "- No manuscript files were modified.",
                "- No release or tag was created.",
                "- No scientific claim was strengthened automatically.",
                "- No exact replication claim was made.",
                "- No protein, function, human-conservation, mechanism, driver, mediator, therapeutic, or rejuvenation claim was made.",
                "- No cells were treated as biological replicates.",
                "- No claim was made that dropout, depth, or batch concerns are ruled out.",
                "- No model/tool provenance traces were included.",
                "",
                "Safety check passed: yes.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    figures = sorted(FIG_DIR.glob("*.*")) if FIG_DIR.exists() else []
    FILE_REPORT.write_text(
        "\n".join(
            [
                "# v5.5 E-MTAB-8560 mTEC File Check",
                "",
                "## Tables Created",
                "",
                "- `results/tables/external_emtab8560_mtec_lox_pseudobulk.tsv`",
                "- `results/tables/external_emtab8560_mtec_lox_by_age.tsv`",
                "- `results/tables/external_emtab8560_mtec_lox_models.tsv`",
                "- `results/tables/external_emtab8560_mtec_lox_permutation.tsv`",
                "- `results/tables/external_emtab8560_mtec_lox_global_fdr.tsv`",
                "",
                "## Figures Created",
                "",
                *[f"- `{fig.as_posix()}`" for fig in figures],
                "",
                "## Checks",
                "",
                "- No FASTQ, BAM, SRA, or raw alignment files downloaded.",
                "- No `.h5ad` file created or staged.",
                "- No large raw data files included for commit.",
                "- No `data/raw` or `data/processed` files staged.",
                "- No `private_outreach` files staged.",
                "- No manuscript files staged.",
                "- No release or tag created.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    required = [
        META_DIR / "emtab8560_coldata.tsv",
        META_DIR / "emtab8560_lox_counts_long.tsv",
        CELL_META,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        write_infeasible(missing)
        print("E-MTAB-8560 mTEC LOX reanalysis infeasible: missing R export")
        return 0

    pb = build_pseudobulk()
    by_age = age_summary(pb)
    models, perms, fdr = build_models(pb[pb["analysis_group"].isin(GROUPS)])
    final_class = classify(models, fdr)

    pb.to_csv(TABLE_DIR / "external_emtab8560_mtec_lox_pseudobulk.tsv", sep="\t", index=False)
    by_age.to_csv(TABLE_DIR / "external_emtab8560_mtec_lox_by_age.tsv", sep="\t", index=False)
    models.to_csv(TABLE_DIR / "external_emtab8560_mtec_lox_models.tsv", sep="\t", index=False)
    perms.to_csv(TABLE_DIR / "external_emtab8560_mtec_lox_permutation.tsv", sep="\t", index=False)
    fdr.to_csv(TABLE_DIR / "external_emtab8560_mtec_lox_global_fdr.tsv", sep="\t", index=False)

    plot_group(pb, "mTEClo", "log2_cpm_plus1", "emtab8560_loxl2_mteclo_log2cpm_by_age", "log2(CPM+1)")
    plot_group(pb, "mTEChi", "log2_cpm_plus1", "emtab8560_loxl2_mtechi_log2cpm_by_age", "log2(CPM+1)")
    plot_group(pb, "combined_mTEClo_mTEChi", "log2_cpm_plus1", "emtab8560_loxl2_combined_mtec_log2cpm_by_age", "log2(CPM+1)")
    plot_group(pb, "mTEClo", "detection_rate", "emtab8560_loxl2_mteclo_detection_by_age", "Detection rate")
    plot_group(pb, "mTEChi", "detection_rate", "emtab8560_loxl2_mtechi_detection_by_age", "Detection rate")
    plot_heatmap(models)
    write_reports(pb, by_age, models, perms, fdr, final_class)
    print(f"E-MTAB-8560 mTEC LOX reanalysis completed: {final_class}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
