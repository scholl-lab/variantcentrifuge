# File: variantcentrifuge/association/diagnostics.py
# Location: variantcentrifuge/variantcentrifuge/association/diagnostics.py
"""
Diagnostics module for association analysis QC.

Provides:
- compute_lambda_gc(): Genomic inflation factor from p-value distribution
- compute_qq_data(): Observed vs expected -log10(p) for QQ plots
- emit_sample_size_warnings(): Cohort-level power/balance warnings
- compute_per_gene_warnings(): Per-gene carrier count warnings
- write_diagnostics(): Write all diagnostics files to output directory

Usage
-----
After engine.run_all() completes, call write_diagnostics() with the results
DataFrame and diagnostics directory path. Lazy-imported by AssociationAnalysisStage
to avoid unnecessary scipy/numpy overhead when diagnostics are not requested.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from variantcentrifuge.association.base import AssociationConfig

logger = logging.getLogger("variantcentrifuge")
from variantcentrifuge.association.diagnostics_qq import (  # noqa: F401
    _EXPECTED_CHI2_MEDIAN,
    compute_lambda_gc,
    compute_qq_data,
    write_qq_plot,
)


def emit_sample_size_warnings(
    n_cases: int,
    n_controls: int,
    config: AssociationConfig,
) -> list[str]:
    """
    Check cohort-level sample size and balance; emit warnings and return flag strings.

    Parameters
    ----------
    n_cases : int
        Total number of case samples.
    n_controls : int
        Total number of control samples.
    config : AssociationConfig
        Configuration with min_cases and max_case_control_ratio thresholds.

    Returns
    -------
    list[str]
        Warning flag strings for any triggered conditions:
        - "LOW_CASE_COUNT" if n_cases < config.min_cases
        - "IMBALANCED_COHORT" if max(n_cases, n_controls)/min(n_cases, n_controls) > config.max_case_control_ratio
    """
    warnings: list[str] = []

    if n_cases < config.min_cases:
        msg = (
            f"Association diagnostics: n_cases={n_cases} is below the recommended "
            f"minimum of {config.min_cases}. Results may be underpowered."
        )
        logger.warning(msg)
        warnings.append("LOW_CASE_COUNT")

    if n_cases > 0 and n_controls > 0:
        ratio = max(n_cases, n_controls) / min(n_cases, n_controls)
        if ratio > config.max_case_control_ratio:
            larger = "cases" if n_cases > n_controls else "controls"
            smaller = "controls" if larger == "cases" else "cases"
            msg = (
                f"Association diagnostics: {larger} outnumber {smaller} {ratio:.1f}:1, "
                f"exceeding the {config.max_case_control_ratio:.0f}:1 threshold. "
                "Type I error inflation risk without SPA/Firth correction."
            )
            logger.warning(msg)
            warnings.append("IMBALANCED_COHORT")

    return warnings


def compute_per_gene_warnings(
    gene: str,
    case_carriers: int,
    config: AssociationConfig,
) -> list[str]:
    """
    Compute per-gene warnings based on carrier counts.

    Parameters
    ----------
    gene : str
        Gene symbol (used for logging).
    case_carriers : int
        Number of case samples carrying at least one qualifying variant.
    config : AssociationConfig
        Configuration with min_case_carriers threshold.

    Returns
    -------
    list[str]
        Warning flag strings for any triggered conditions:
        - "LOW_CARRIER_COUNT" if case_carriers < config.min_case_carriers
    """
    warnings: list[str] = []

    if case_carriers < config.min_case_carriers:
        logger.warning(
            f"Association diagnostics: gene {gene} has only {case_carriers} "
            f"case carrier(s) (threshold: {config.min_case_carriers}). "
            "Gene-level result may be unreliable."
        )
        warnings.append("LOW_CARRIER_COUNT")

    return warnings


def write_diagnostics(
    results_df: pd.DataFrame,
    diagnostics_dir: str | Path,
    test_names: list[str],
    n_cases: int,
    n_controls: int,
    cohort_warnings: list[str],
) -> None:
    """
    Write diagnostics files to the specified directory.

    Creates three files:
    - lambda_gc.tsv: lambda_GC per test with n_tests column
    - qq_data.tsv: concatenated QQ data for all tests
    - summary.txt: human-readable summary with sample sizes and lambda_GC values

    Parameters
    ----------
    results_df : pd.DataFrame
        Association results DataFrame from engine.run_all(). Must have a "gene"
        column plus p-value columns named ``{test_name}_pvalue`` and optionally
        ``acat_o_pvalue``.
    diagnostics_dir : str | Path
        Directory to write output files. Created if it does not exist.
    test_names : list[str]
        List of active test names (e.g. ["fisher", "skat_python"]).
    n_cases : int
        Total case sample count.
    n_controls : int
        Total control sample count.
    cohort_warnings : list[str]
        Cohort-level warning flags from emit_sample_size_warnings().
    """
    diag_dir = Path(diagnostics_dir)
    diag_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Determine which p-value columns are present in results_df
    # ------------------------------------------------------------------
    # Primary tests: {test_name}_pvalue
    # ACAT-O: acat_o_pvalue (always added if ACAT-O was run)
    all_test_ids: list[str] = list(test_names)
    if "acat_o_pvalue" in results_df.columns and "acat_o" not in all_test_ids:
        all_test_ids.append("acat_o")

    # Build mapping: test_id -> column name
    p_col_map: dict[str, str] = {}
    for tid in all_test_ids:
        col = f"{tid}_pvalue"
        if col in results_df.columns:
            p_col_map[tid] = col

    # ------------------------------------------------------------------
    # Compute lambda_GC and QQ data per test
    # Note: lambda_GC for Fisher is included here as a diagnostic metric
    # only. Do NOT use it to inflate-correct Fisher p-values — Fisher's
    # exact test is not asymptotic and GC correction is statistically
    # invalid for exact tests. See compute_lambda_gc() Notes for details.
    # ------------------------------------------------------------------
    lambda_rows: list[dict] = []
    qq_frames: list[pd.DataFrame] = []

    for tid, col in p_col_map.items():
        p_values: list[float | None] = results_df[col].tolist()
        lam = compute_lambda_gc(p_values)
        n_valid = sum(1 for p in p_values if p is not None and not np.isnan(p))

        if lam is not None:
            lambda_rows.append({"test_name": tid, "lambda_gc": lam, "n_tests": n_valid})

        qq_df = compute_qq_data(p_values, test_name=tid)
        if not qq_df.empty:
            qq_frames.append(qq_df)

    # ------------------------------------------------------------------
    # Write lambda_gc.tsv
    # ------------------------------------------------------------------
    lambda_df = pd.DataFrame(lambda_rows, columns=["test_name", "lambda_gc", "n_tests"])
    lambda_df.to_csv(diag_dir / "lambda_gc.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------
    # Write qq_data.tsv
    # ------------------------------------------------------------------
    if qq_frames:
        qq_combined = pd.concat(qq_frames, ignore_index=True)
    else:
        qq_combined = pd.DataFrame(columns=["test", "expected_neg_log10_p", "observed_neg_log10_p"])
    qq_combined.to_csv(diag_dir / "qq_data.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------
    # Write summary.txt
    # ------------------------------------------------------------------
    lines: list[str] = [
        "Association Analysis Diagnostics Summary",
        "=" * 42,
        "",
        "Sample sizes:",
        f"  n_cases    = {n_cases}",
        f"  n_controls = {n_controls}",
        "",
        "Active tests:",
    ]
    for tid in all_test_ids:
        lines.append(f"  - {tid}")

    lines.append("")
    lines.append("Lambda_GC (genomic inflation factor):")

    if lambda_rows:
        for row in lambda_rows:
            lam_str = f"{row['lambda_gc']:.4f}"
            n_str = f"(n={row['n_tests']})"
            flag = ""
            if row["n_tests"] < 100:
                flag = " [UNRELIABLE: n < 100]"
            lines.append(f"  {row['test_name']:<20s} {lam_str}  {n_str}{flag}")
    else:
        lines.append("  (no valid p-values available)")

    if cohort_warnings:
        lines.append("")
        lines.append("Cohort-level warnings:")
        for w in cohort_warnings:
            lines.append(f"  WARNING: {w}")

    # Count genes with per-gene warnings
    if "warnings" in results_df.columns:
        n_flagged = int((results_df["warnings"] != "").sum())
        lines.append("")
        lines.append(f"Genes with per-gene warnings: {n_flagged}")

    lines.append("")

    summary_path = diag_dir / "summary.txt"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ------------------------------------------------------------------
    # Optional QQ plot PNG (DIAG-04): requires matplotlib; skipped if absent
    # ------------------------------------------------------------------
    qq_plot_path = diag_dir / "qq_plot.png"
    write_qq_plot(qq_combined, qq_plot_path)

    logger.info(f"Diagnostics written to {diag_dir}")


def write_fdr_weight_diagnostics(
    genes: list[str],
    raw_weights: dict[str, float],
    normalized_weights: np.ndarray,
    unweighted_corrected: np.ndarray,
    weighted_corrected: np.ndarray,
    fdr_threshold: float,
    diagnostics_dir: str | Path,
) -> None:
    """
    Write per-gene FDR weight diagnostics file.

    Creates ``fdr_weight_diagnostics.tsv`` in the diagnostics directory with
    one row per gene, detailing raw weight, normalized weight, unweighted q-value,
    weighted q-value, and whether significance changed at the given threshold.

    Parameters
    ----------
    genes : list of str
        Gene names aligned 1:1 with ``normalized_weights``, ``unweighted_corrected``,
        and ``weighted_corrected``.
    raw_weights : dict[str, float]
        Mapping from gene name to raw (pre-normalization) weight from the weight file.
        Genes not in the map received the default weight=1.0.
    normalized_weights : np.ndarray
        Per-gene normalized weights (mean=1.0) aligned with ``genes``.
    unweighted_corrected : np.ndarray
        Unweighted BH/Bonferroni corrected q-values aligned with ``genes``.
    weighted_corrected : np.ndarray
        Weighted BH/Bonferroni corrected q-values aligned with ``genes``.
    fdr_threshold : float
        Significance threshold (e.g. 0.05) used to classify significance changes.
    diagnostics_dir : str or Path
        Directory to write the output file. Created if it does not exist.
    """
    diag_dir = Path(diagnostics_dir)
    diag_dir.mkdir(parents=True, exist_ok=True)

    m = len(genes)
    rows: list[dict] = []
    n_gained = 0
    n_lost = 0
    n_with_file_weight = 0

    for i, gene in enumerate(genes):
        raw_w = raw_weights.get(gene, 1.0)
        if gene in raw_weights:
            n_with_file_weight += 1
        norm_w = float(normalized_weights[i])
        unw_q = float(unweighted_corrected[i])
        wt_q = float(weighted_corrected[i])

        was_sig = unw_q < fdr_threshold
        is_sig = wt_q < fdr_threshold

        if not was_sig and is_sig:
            change = "gained"
            n_gained += 1
        elif was_sig and not is_sig:
            change = "lost"
            n_lost += 1
        else:
            change = ""

        rows.append(
            {
                "gene": gene,
                "raw_weight": raw_w,
                "normalized_weight": norm_w,
                "unweighted_q": unw_q,
                "weighted_q": wt_q,
                "significance_change": change,
            }
        )

    diag_df = pd.DataFrame(
        rows,
        columns=[
            "gene",
            "raw_weight",
            "normalized_weight",
            "unweighted_q",
            "weighted_q",
            "significance_change",
        ],
    )
    out_path = diag_dir / "fdr_weight_diagnostics.tsv"
    diag_df.to_csv(out_path, sep="\t", index=False)

    # Summary log
    logger.info(
        f"FDR weighting impact: {n_gained} genes gained significance, "
        f"{n_lost} genes lost significance (threshold={fdr_threshold})"
    )

    # Weight distribution log
    norm_arr = np.asarray(normalized_weights, dtype=float)
    if len(norm_arr) > 0:
        logger.info(
            f"Normalized weight distribution (n={m}): "
            f"min={float(norm_arr.min()):.3f}, "
            f"max={float(norm_arr.max()):.3f}, "
            f"median={float(np.median(norm_arr)):.3f}, "
            f"mean={float(norm_arr.mean()):.3f}"
        )

    # Coverage log
    n_default = m - n_with_file_weight
    logger.info(
        f"FDR weight coverage: {n_with_file_weight} genes with weights from file, "
        f"{n_default} genes at default weight=1.0"
    )

    # Effective number of tests
    if len(norm_arr) > 0 and (norm_arr**2).sum() > 0:
        n_eff = float(norm_arr.sum() ** 2 / (norm_arr**2).sum())
        logger.info(f"Effective number of tests (sum(w)^2/sum(w^2)): {n_eff:.1f}")

    logger.info(f"FDR weight diagnostics written to {out_path}")
