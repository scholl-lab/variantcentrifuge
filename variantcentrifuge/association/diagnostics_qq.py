# File: variantcentrifuge/association/diagnostics_qq.py
"""lambda_GC and QQ statistics and plot (split from diagnostics.py; re-exported there)."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist

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

    Notes
    -----
    **Fisher's exact test exemption:** Fisher p-values must NOT be corrected
    via lambda_GC. Fisher's exact test is a permutation-based exact test, not
    an asymptotic chi2 approximation — the GC correction formula assumes the
    chi2(1) asymptotic distribution and is statistically invalid for exact
    tests. This function computes lambda_GC for Fisher as a *diagnostic only*;
    callers must not apply lambda_GC correction to Fisher p-values.

    **Applicable tests:** GC correction is valid only for score-based
    asymptotic tests (e.g., SKAT, burden) whose statistics follow chi2(1)
    under the null.

    **Inflation thresholds (approximate guidance):**
    - lambda_GC ~ 1.0: well-calibrated
    - lambda_GC > 1.2: notable inflation — investigate population stratification
      or cryptic relatedness before interpreting results

    **Over-correction risk:** In studies with genuine rare-variant signal or
    small sample sizes, lambda_GC deflation from true associations can reduce
    power (Devlin & Roeder 1999, Gorroochurn et al. 2025). Use GC correction
    only when inflation is clearly present and attributable to confounding.
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
