#!/usr/bin/env python
"""Validate the ThymusLOXScan analysis pipeline end to end."""

from __future__ import annotations

import importlib
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REQUIREMENTS_PATH = PROJECT_ROOT / "requirements.txt"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURES_DIR = PROJECT_ROOT / "figures"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

ANNOTATED_H5AD = PROCESSED_DIR / "thymus_annotated.h5ad"
PREPROCESSED_H5AD = PROCESSED_DIR / "thymus_preprocessed.h5ad"
DIFF_RESULTS_CSV = PROCESSED_DIR / "LOX_differential_results.csv"

LOX_HUMAN = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"]
LOX_MOUSE = ["Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4"]
EXPECTED_CELL_TYPES = ["cTEC", "mTEC", "fibroblasts"]

IMPORT_NAME_OVERRIDES = {
    "python-igraph": "igraph",
}


@dataclass
class CheckResult:
    name: str
    passed: bool
    messages: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    fixes: list[str] = field(default_factory=list)

    def status(self) -> str:
        if self.passed and self.warnings:
            return "PASS with WARNINGS"
        return "PASS" if self.passed else "FAIL"


def print_result(result: CheckResult) -> None:
    print(f"\n[{result.status()}] {result.name}")
    for message in result.messages:
        print(f"  - {message}")
    for warning in result.warnings:
        print(f"  WARNING: {warning}")
    for fix in result.fixes:
        print(f"  FIX: {fix}")


def read_requirements() -> list[str]:
    requirements = []
    for line in REQUIREMENTS_PATH.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        package = stripped.split("==")[0].split(">=")[0].split("<=")[0].split("~=")[0].strip()
        requirements.append(package)
    return requirements


def import_name(package: str) -> str:
    return IMPORT_NAME_OVERRIDES.get(package, package.replace("-", "_"))


def check_environment() -> CheckResult:
    result = CheckResult("Check 1 - Environment", passed=True)
    missing = []

    for package in read_requirements():
        module_name = import_name(package)
        try:
            importlib.import_module(module_name)
            result.messages.append(f"Import OK: {package} as {module_name}")
        except Exception as exc:
            missing.append((package, module_name, exc))

    for package, module_name, exc in missing:
        result.passed = False
        result.messages.append(f"Import failed: {package} as {module_name} ({exc})")

    for module_name in ["scanpy", "anndata", "pandas", "numpy"]:
        try:
            module = importlib.import_module(module_name)
            result.messages.append(f"{module_name} version: {getattr(module, '__version__', 'unknown')}")
        except Exception:
            pass

    if missing:
        result.fixes.append("Install dependencies with: pip install -r requirements.txt")

    return result


def load_annotated(result: CheckResult):
    if not ANNOTATED_H5AD.exists():
        result.passed = False
        result.messages.append(f"Missing file: {ANNOTATED_H5AD}")
        result.fixes.append("Run notebooks/03_preprocessing_thymus.ipynb and notebooks/04_cell_type_annotation.ipynb.")
        return None

    try:
        import scanpy as sc

        return sc.read_h5ad(ANNOTATED_H5AD)
    except Exception as exc:
        result.passed = False
        result.messages.append(f"Could not load {ANNOTATED_H5AD}: {exc}")
        result.fixes.append("Regenerate thymus_annotated.h5ad from notebook 04 and confirm the file is a valid AnnData object.")
        return None


def check_data_integrity() -> tuple[CheckResult, object | None]:
    result = CheckResult("Check 2 - Data integrity", passed=True)
    adata = load_annotated(result)
    if adata is None:
        return result, None

    result.messages.append(f"Loaded annotated object: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    conditions = [
        (adata.n_obs > 1000, "n_obs > 1000", "Use the full thymus dataset or revisit QC thresholds if too many cells were filtered."),
        (adata.n_vars > 5000, "n_vars > 5000", "Confirm that the saved annotated object retains enough genes, ideally via adata.raw or by saving before HVG-only subsetting."),
        ("age_group" in adata.obs.columns, "'age_group' column exists in adata.obs", "Create or rename an age metadata column to adata.obs['age_group']."),
        ("cell_type" in adata.obs.columns, "'cell_type' column exists in adata.obs", "Run notebook 04 and create adata.obs['cell_type']."),
        (len(adata.obs["cell_type"].dropna().unique()) >= 4 if "cell_type" in adata.obs else False, "at least 4 distinct cell types annotated", "Review CellTypist and marker-based annotation in notebook 04."),
        (len(adata.obs["age_group"].dropna().unique()) >= 2 if "age_group" in adata.obs else False, "at least 2 distinct age groups present", "Map the dataset age metadata into young/old labels in adata.obs['age_group']."),
    ]

    for passed, label, fix in conditions:
        result.messages.append(f"{'OK' if passed else 'FAILED'}: {label}")
        if not passed:
            result.passed = False
            result.fixes.append(fix)

    present_human = [gene for gene in LOX_HUMAN if gene in adata.var_names]
    missing_human = [gene for gene in LOX_HUMAN if gene not in adata.var_names]
    present_mouse = [gene for gene in LOX_MOUSE if gene in adata.var_names]
    result.messages.append(f"Human-style LOX symbols present in adata.var_names: {present_human}")
    result.messages.append(f"Mouse-style LOX symbols present in adata.var_names: {present_mouse}")

    if missing_human:
        result.passed = False
        result.messages.append(f"FAILED: missing exact requested LOX symbols: {missing_human}")
        result.fixes.append(
            "If this is mouse data, either convert gene symbols to uppercase human-style aliases for this check or update downstream notebooks to use mouse symbols consistently."
        )
    else:
        result.messages.append("OK: all exact requested LOX symbols are present in adata.var_names")

    return result, adata


def check_cell_type_proportions(adata) -> CheckResult:
    result = CheckResult("Check 3 - Cell type proportions sanity", passed=True)
    if adata is None:
        result.passed = False
        result.fixes.append("Fix Check 2 first so the annotated AnnData object can be loaded.")
        return result

    if "cell_type" not in adata.obs or "age_group" not in adata.obs:
        result.passed = False
        result.messages.append("Missing cell_type or age_group column.")
        result.fixes.append("Run annotation notebook and ensure adata.obs has both 'cell_type' and 'age_group'.")
        return result

    counts = adata.obs["cell_type"].astype(str).value_counts()
    percentages = counts / counts.sum() * 100
    result.messages.append("Cell type counts and percentages:")
    for cell_type, count in counts.items():
        result.messages.append(f"{cell_type}: {count:,} cells ({percentages[cell_type]:.2f}%)")

    lower_counts = {cell_type.lower(): count for cell_type, count in counts.items()}
    for expected in EXPECTED_CELL_TYPES:
        count = lower_counts.get(expected.lower(), 0)
        if count < 50:
            result.passed = False
            result.warnings.append(f"{expected} has fewer than 50 cells ({count}); differential expression may be unreliable.")
            result.fixes.append(f"Review annotation and filtering for {expected}; consider broader labels or obtain more cells.")

    age_counts = adata.obs["age_group"].astype(str).value_counts()
    result.messages.append("Age group counts:")
    for age_group, count in age_counts.items():
        result.messages.append(f"{age_group}: {count:,} cells")

    if len(age_counts) >= 2:
        imbalance_ratio = age_counts.max() / max(age_counts.min(), 1)
        result.messages.append(f"Age imbalance ratio: {imbalance_ratio:.2f}x")
        if imbalance_ratio > 10:
            result.passed = False
            result.messages.append("FAILED: age groups are severely imbalanced (>10x difference).")
            result.fixes.append("Check age_group mapping, sample inclusion, and QC filtering by age group.")
    else:
        result.passed = False
        result.fixes.append("Ensure at least two age groups exist in adata.obs['age_group'].")

    return result


def expression_vector(adata, gene: str) -> np.ndarray:
    matrix = adata[:, gene].X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix).ravel()


def check_lox_expression(adata) -> CheckResult:
    result = CheckResult("Check 4 - LOX expression sanity", passed=True)
    if adata is None:
        result.passed = False
        result.fixes.append("Fix Check 2 first so LOX expression can be evaluated.")
        return result

    expr_source = adata.raw.to_adata() if adata.raw is not None else adata
    gene_lookup = {gene.upper(): gene for gene in expr_source.var_names.astype(str)}
    human_present = [gene for gene in LOX_HUMAN if gene in expr_source.var_names]
    mouse_present = [gene for gene in LOX_MOUSE if gene in expr_source.var_names]
    auto_genes = [gene_lookup.get(gene.upper()) for gene in LOX_HUMAN]
    auto_genes = [gene for gene in auto_genes if gene is not None]

    result.messages.append(f"Human capitalization present: {human_present}")
    result.messages.append(f"Mouse capitalization present: {mouse_present}")
    result.messages.append(f"Auto-selected LOX genes for expression sanity: {auto_genes}")

    if not auto_genes:
        result.passed = False
        result.messages.append("FAILED: no LOX family symbols found under human or mouse capitalization.")
        result.fixes.append("Inspect adata.var_names and map Ensembl IDs or alternate gene symbols to LOX family symbols.")
        return result

    expressed_percentages = []
    for gene in auto_genes:
        values = expression_vector(expr_source, gene)
        pct_expressing = float(np.mean(values > 0) * 100)
        mean_expression = float(np.mean(values))
        max_expression = float(np.max(values))
        expressed_percentages.append(pct_expressing)
        result.messages.append(
            f"{gene}: {pct_expressing:.2f}% cells expressing, mean={mean_expression:.4f}, max={max_expression:.4f}"
        )

    if len(auto_genes) == 5 and all(pct == 0 for pct in expressed_percentages):
        result.warnings.append(
            "All 5 LOX family genes show 0% expression; this is likely a gene-name mismatch or expression matrix issue."
        )
        result.fixes.append("Check whether the dataset uses mouse symbols, Ensembl IDs, or whether adata.raw contains the full gene matrix.")

    return result


def check_output_files() -> CheckResult:
    result = CheckResult("Check 5 - Figures and processed outputs", passed=True)
    png_files = sorted(FIGURES_DIR.rglob("*.png"))
    result.messages.append(f"Found {len(png_files)} PNG files under {FIGURES_DIR}")

    if len(png_files) < 5:
        result.passed = False
        result.messages.append("FAILED: fewer than 5 .png files found.")
        result.fixes.append("Run notebooks 03, 04, and 05 to generate QC, annotation, and LOX analysis figures.")

    required_files = [PREPROCESSED_H5AD, ANNOTATED_H5AD, DIFF_RESULTS_CSV]
    for path in required_files:
        if path.exists():
            result.messages.append(f"OK: {path.relative_to(PROJECT_ROOT)} exists")
        else:
            result.passed = False
            result.messages.append(f"FAILED: missing {path.relative_to(PROJECT_ROOT)}")

    if not PREPROCESSED_H5AD.exists():
        result.fixes.append("Run notebooks/03_preprocessing_thymus.ipynb to create thymus_preprocessed.h5ad.")
    if not ANNOTATED_H5AD.exists():
        result.fixes.append("Run notebooks/04_cell_type_annotation.ipynb to create thymus_annotated.h5ad.")
    if not DIFF_RESULTS_CSV.exists():
        result.fixes.append("Run notebooks/05_LOX_expression_analysis.ipynb to create LOX_differential_results.csv.")

    return result


def check_reproducibility() -> CheckResult:
    result = CheckResult("Check 6 - Reproducibility", passed=True)
    notebooks = sorted(NOTEBOOKS_DIR.glob("*.ipynb"))
    missing_seed = []

    for notebook in notebooks:
        text = notebook.read_text(encoding="utf-8")
        has_random_state = "random_state" in text
        has_np_seed = "np.random.seed" in text
        result.messages.append(
            f"{notebook.name}: random_state={'yes' if has_random_state else 'no'}, np.random.seed={'yes' if has_np_seed else 'no'}"
        )
        if not (has_random_state or has_np_seed):
            missing_seed.append(notebook.name)

    if missing_seed:
        result.warnings.append(
            "Random seeds are missing in notebooks: " + ", ".join(missing_seed) + ". Results may not be fully reproducible."
        )
        result.fixes.append("Add np.random.seed(0) near notebook setup and pass random_state=0 to UMAP, Leiden, and other stochastic steps.")

    return result


def main() -> int:
    print("ThymusLOXScan pipeline validation")
    print(f"Project root: {PROJECT_ROOT}")

    results = []
    env_result = check_environment()
    results.append(env_result)

    data_result, adata = check_data_integrity()
    results.append(data_result)
    results.append(check_cell_type_proportions(adata))
    results.append(check_lox_expression(adata))
    results.append(check_output_files())
    results.append(check_reproducibility())

    for result in results:
        print_result(result)

    passed_count = sum(result.passed for result in results)
    total_count = len(results)
    print(f"\nPIPELINE STATUS: {passed_count}/{total_count} checks passed")

    failed_results = [result for result in results if not result.passed]
    if failed_results:
        print("\nFailures to fix:")
        for result in failed_results:
            print(f"- {result.name}")
            for fix in dict.fromkeys(result.fixes):
                print(f"  FIX: {fix}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
