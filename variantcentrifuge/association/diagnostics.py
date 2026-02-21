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
from scipy.stats import chi2 as chi2_dist

if TYPE_CHECKING:
    from variantcentrifuge.association.base import AssociationConfig

logger = logging.getLogger("variantcentrifuge")

# Hardcoded median of chi2(df=1) distribution.
# Equivalent to: from scipy.stats import chi2; chi2.ppf(0.5, df=1)
# Hardcoded to avoid calling chi2.ppf at import time (cold import cost).
_EXPECTED_CHI2_MEDIAN: float = 0.45493642311957174


def compute_lambda_gc(p_values: list[float | None]) -> float | None:
    """
    Compute genomic inflation factor (lambda_GC) from a list of p-values.

    Lambda_GC is the ratio of the observed median chi2(1) statistic to the
    expected median under the null. A well-calibrated analysis has lambda_GC
    near 1.0. Values > 1.05 suggest inflation (population stratification,
    cryptic relatedness) or true signal.

    Parameters
    ----------
    p_values : list[float | None]
        Raw (uncorrected) p-values. None and NaN values are filtered out.

    Returns
    -------
    float | None
        lambda_GC, or None if fewer than 2 valid p-values are available.
    """
    # Filter out None and NaN
    valid = np.array(
        [p for p in p_values if p is not None and not np.isnan(p)],
        dtype=float,
    )

    if len(valid) < 2:
        return None

    # Clip to avoid -inf in log / chi2 conversion
    valid = np.clip(valid, 1e-300, 1.0 - 1e-15)

    n_valid = len(valid)
    if n_valid < 100:
        logger.warning(f"lambda_GC computed on {n_valid} tests — unreliable for n < 100")

    # Convert p-values to chi2(df=1) statistics via survival function inverse
    chi2_obs = chi2_dist.isf(valid, df=1)

    return float(np.median(chi2_obs) / _EXPECTED_CHI2_MEDIAN)


def compute_qq_data(p_values: list[float | None], test_name: str) -> pd.DataFrame:
    """
    Compute observed vs expected -log10(p) data for QQ plots.

    Uses the Hazen quantile formula for expected quantiles: i/(n+1) for
    rank i in 1..n. Returns data sorted ascending by expected_neg_log10_p
    (non-significant end first, significant end last) — correct order for
    sequential QQ plot rendering.

    Parameters
    ----------
    p_values : list[float | None]
        Raw p-values. None and NaN are excluded.
    test_name : str
        Test identifier written to the "test" column.

    Returns
    -------
    pd.DataFrame
        Columns: "test", "expected_neg_log10_p", "observed_neg_log10_p".
        Empty DataFrame with same columns if no valid p-values.
    """
    _empty = pd.DataFrame(columns=["test", "expected_neg_log10_p", "observed_neg_log10_p"])

    valid = np.array(
        [p for p in p_values if p is not None and not np.isnan(p)],
        dtype=float,
    )

    if len(valid) == 0:
        return _empty

    n = len(valid)
    # Sort ascending so smallest observed p is at the end of QQ plot
    observed_sorted = np.sort(valid)

    # Hazen quantile formula: expected quantile for rank i (1-indexed) is i/(n+1)
    ranks = np.arange(1, n + 1, dtype=float)
    expected_p = ranks / (n + 1)

    # Convert to -log10, clip observed to avoid +inf
    expected_neg_log10 = -np.log10(expected_p)
    observed_clipped = np.clip(observed_sorted, 1e-300, 1.0)
    observed_neg_log10 = -np.log10(observed_clipped)

    df = pd.DataFrame(
        {
            "test": test_name,
            "expected_neg_log10_p": expected_neg_log10,
            "observed_neg_log10_p": observed_neg_log10,
        }
    )

    # Sort by expected ascending (smallest first = non-significant end)
    df = df.sort_values("expected_neg_log10_p", ascending=True).reset_index(drop=True)

    return df


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
        - "IMBALANCED_COHORT" if n_controls/n_cases > config.max_case_control_ratio
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
        ratio = n_controls / n_cases
        if ratio > config.max_case_control_ratio:
            msg = (
                f"Association diagnostics: case:control ratio is 1:{ratio:.1f}, "
                f"exceeding the 1:{config.max_case_control_ratio:.0f} threshold. "
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


def write_qq_plot(
    qq_data: pd.DataFrame,
    output_path: str | Path,
) -> bool:
    """Write QQ plot as PNG. Returns True if successful, False if matplotlib absent.

    Uses lazy import with matplotlib.use("Agg") for headless HPC environments.
    The Agg backend MUST be set before importing matplotlib.pyplot to avoid
    display errors on systems without a graphical environment.

    Parameters
    ----------
    qq_data : pd.DataFrame
        QQ data with columns: test, expected_neg_log10_p, observed_neg_log10_p.
        Typically produced by concatenating compute_qq_data() results across tests.
    output_path : str | Path
        Output file path (e.g. /path/to/diagnostics/qq_plot.png).

    Returns
    -------
    bool
        True if plot was written successfully, False if matplotlib is not
        installed or qq_data is empty.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")  # headless HPC — MUST be before pyplot import
        import matplotlib.pyplot as plt
    except ImportError:
        logger.info("matplotlib not installed — QQ plot skipped")
        return False

    if qq_data.empty:
        logger.info("No QQ data available — QQ plot skipped")
        return False

    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot each test as a separate series
    for test_name, group in qq_data.groupby("test"):
        ax.scatter(
            group["expected_neg_log10_p"],
            group["observed_neg_log10_p"],
            s=4,
            alpha=0.6,
            label=str(test_name),
        )

    # Identity line: expected == observed under the null hypothesis
    max_val = (
        max(
            qq_data["expected_neg_log10_p"].max(),
            qq_data["observed_neg_log10_p"].max(),
        )
        + 0.5
    )
    ax.plot([0, max_val], [0, max_val], "k--", linewidth=0.8, label="Expected")

    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title("QQ Plot")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()

    plt.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info(f"QQ plot written to {output_path}")
    return True


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
        column plus p-value columns named ``{test_name}_p_value`` and optionally
        ``acat_o_p_value``.
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
    # Primary tests: {test_name}_p_value
    # ACAT-O: acat_o_p_value (always added if ACAT-O was run)
    all_test_ids: list[str] = list(test_names)
    if "acat_o_p_value" in results_df.columns and "acat_o" not in all_test_ids:
        all_test_ids.append("acat_o")

    # Build mapping: test_id -> column name
    p_col_map: dict[str, str] = {}
    for tid in all_test_ids:
        col = f"{tid}_p_value"
        if col in results_df.columns:
            p_col_map[tid] = col

    # ------------------------------------------------------------------
    # Compute lambda_GC and QQ data per test
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
