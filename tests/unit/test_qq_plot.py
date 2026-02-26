"""
Unit tests for write_qq_plot() in variantcentrifuge.association.diagnostics.

Tests:
- QQ plot written when matplotlib is available
- Returns False when matplotlib is absent (mocked ImportError)
- Returns False for empty DataFrame
- write_diagnostics() calls write_qq_plot() and produces qq_plot.png
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.diagnostics import write_diagnostics, write_qq_plot

# Check if matplotlib is available (optional dependency)
try:
    import matplotlib  # noqa: F401

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_qq_data() -> pd.DataFrame:
    """Sample QQ data with two test groups (fisher and skat_python)."""
    rng = np.random.default_rng(seed=42)
    n = 50

    # Group 1: fisher — null p-values
    null_p = rng.uniform(0, 1, size=n)
    ranks1 = np.arange(1, n + 1, dtype=float)
    expected1 = ranks1 / (n + 1)
    df_fisher = pd.DataFrame(
        {
            "test": "fisher",
            "expected_neg_log10_p": -np.log10(expected1),
            "observed_neg_log10_p": -np.log10(np.sort(null_p)),
        }
    )

    # Group 2: skat_python — slightly inflated p-values
    inflated_p = rng.beta(0.5, 1.0, size=n)
    df_skat = pd.DataFrame(
        {
            "test": "skat_python",
            "expected_neg_log10_p": -np.log10(ranks1 / (n + 1)),
            "observed_neg_log10_p": -np.log10(np.clip(np.sort(inflated_p), 1e-300, 1.0)),
        }
    )

    return pd.concat([df_fisher, df_skat], ignore_index=True)


@pytest.fixture
def minimal_results_df() -> pd.DataFrame:
    """Minimal association results DataFrame for write_diagnostics() tests."""
    rng = np.random.default_rng(seed=7)
    p_vals = rng.uniform(0, 1, size=20).tolist()
    return pd.DataFrame(
        {
            "gene": [f"GENE{i}" for i in range(20)],
            "fisher_pvalue": p_vals,
        }
    )


# ---------------------------------------------------------------------------
# Tests: write_qq_plot with matplotlib available
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not installed")
class TestWriteQQPlotWithMatplotlib:
    """write_qq_plot() when matplotlib is installed."""

    def test_produces_png_file(self, sample_qq_data, tmp_path):
        """Output PNG file exists and has non-zero size."""
        output_path = tmp_path / "qq_plot.png"
        result = write_qq_plot(sample_qq_data, output_path)

        assert result is True, "write_qq_plot() should return True when matplotlib available"
        assert output_path.exists(), "qq_plot.png should exist after successful write"
        assert output_path.stat().st_size > 0, "qq_plot.png should have non-zero size"

    def test_returns_true(self, sample_qq_data, tmp_path):
        """Function returns True when plot is written successfully."""
        output_path = tmp_path / "qq.png"
        assert write_qq_plot(sample_qq_data, output_path) is True

    def test_accepts_string_path(self, sample_qq_data, tmp_path):
        """output_path as string (not Path) is accepted."""
        output_path = str(tmp_path / "qq_str.png")
        result = write_qq_plot(sample_qq_data, output_path)
        assert result is True
        assert Path(output_path).exists()

    def test_single_test_group(self, tmp_path):
        """Single test group (no groupby edge cases)."""
        n = 20
        ranks = np.arange(1, n + 1, dtype=float)
        df = pd.DataFrame(
            {
                "test": "fisher",
                "expected_neg_log10_p": -np.log10(ranks / (n + 1)),
                "observed_neg_log10_p": -np.log10(ranks / (n + 1)),  # perfect calibration
            }
        )
        output_path = tmp_path / "single_test.png"
        result = write_qq_plot(df, output_path)
        assert result is True
        assert output_path.exists()

    def test_logs_written_path(self, sample_qq_data, tmp_path, caplog):
        """Log message includes info about the written file."""
        import logging

        output_path = tmp_path / "qq_log_test.png"
        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            write_qq_plot(sample_qq_data, output_path)
        assert any("QQ plot" in record.message for record in caplog.records), (
            f"Expected QQ plot log message, got: {[r.message for r in caplog.records]}"
        )


# ---------------------------------------------------------------------------
# Tests: write_qq_plot with matplotlib absent
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestWriteQQPlotWithoutMatplotlib:
    """write_qq_plot() when matplotlib is not installed."""

    def test_returns_false_when_matplotlib_missing(self, sample_qq_data, tmp_path):
        """Returns False and does not raise when matplotlib is absent."""
        output_path = tmp_path / "qq_absent.png"

        # Remove matplotlib from sys.modules and patch the diagnostics module's
        # __import__ by making matplotlib raise ImportError on import.

        original_modules = sys.modules.copy()
        # Temporarily suppress matplotlib from sys.modules to force ImportError
        # on the lazy import inside write_qq_plot()
        try:
            # Remove matplotlib and its submodules from sys.modules
            mpl_keys = [k for k in sys.modules if k == "matplotlib" or k.startswith("matplotlib.")]
            for k in mpl_keys:
                sys.modules.pop(k, None)

            # Reload to clear cached module state in diagnostics
            # (matplotlib was already imported above in the test process)
            # Instead, just patch sys.modules["matplotlib"] to trigger ImportError
            sys.modules["matplotlib"] = None  # type: ignore[assignment]

            result = write_qq_plot(sample_qq_data, output_path)
        finally:
            # Restore all matplotlib modules
            sys.modules.update(original_modules)

        assert result is False
        assert not output_path.exists(), "No file should be written when matplotlib is absent"

    def test_logs_info_when_matplotlib_missing(self, sample_qq_data, tmp_path, caplog):
        """Logs INFO message about skipped QQ plot when matplotlib is absent."""
        import logging

        output_path = tmp_path / "qq_log_absent.png"

        original_modules = sys.modules.copy()
        try:
            mpl_keys = [k for k in sys.modules if k == "matplotlib" or k.startswith("matplotlib.")]
            for k in mpl_keys:
                sys.modules.pop(k, None)
            sys.modules["matplotlib"] = None  # type: ignore[assignment]

            with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
                write_qq_plot(sample_qq_data, output_path)
        finally:
            sys.modules.update(original_modules)

        messages = [r.message for r in caplog.records]
        assert any("matplotlib" in msg and "skipped" in msg for msg in messages), (
            f"Expected matplotlib skip message, got: {messages}"
        )


# ---------------------------------------------------------------------------
# Tests: write_qq_plot with empty DataFrame
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestWriteQQPlotEmptyData:
    """write_qq_plot() when qq_data is empty."""

    def test_returns_false_for_empty_dataframe(self, tmp_path):
        """Empty DataFrame: returns False without writing file."""
        empty_df = pd.DataFrame(columns=["test", "expected_neg_log10_p", "observed_neg_log10_p"])
        output_path = tmp_path / "qq_empty.png"
        result = write_qq_plot(empty_df, output_path)

        assert result is False
        assert not output_path.exists(), "No file should be written for empty DataFrame"

    def test_logs_info_for_empty_dataframe(self, tmp_path, caplog):
        """Logs INFO about no QQ data when DataFrame is empty."""
        import logging

        empty_df = pd.DataFrame(columns=["test", "expected_neg_log10_p", "observed_neg_log10_p"])
        output_path = tmp_path / "qq_empty_log.png"

        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            write_qq_plot(empty_df, output_path)

        messages = [r.message for r in caplog.records]
        assert any("QQ" in msg and "skipped" in msg for msg in messages), (
            f"Expected QQ skip message, got: {messages}"
        )


# ---------------------------------------------------------------------------
# Tests: write_diagnostics produces qq_plot.png
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestWriteDiagnosticsQQPlot:
    """write_diagnostics() produces qq_plot.png alongside existing TSV files."""

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not installed")
    def test_produces_qq_plot_alongside_tsv(self, minimal_results_df, tmp_path):
        """write_diagnostics() writes qq_plot.png, lambda_gc.tsv, qq_data.tsv, summary.txt."""
        diag_dir = tmp_path / "diagnostics"

        write_diagnostics(
            results_df=minimal_results_df,
            diagnostics_dir=diag_dir,
            test_names=["fisher"],
            n_cases=10,
            n_controls=10,
            cohort_warnings=[],
        )

        # Check all expected files exist
        assert (diag_dir / "lambda_gc.tsv").exists()
        assert (diag_dir / "qq_data.tsv").exists()
        assert (diag_dir / "summary.txt").exists()
        # QQ plot should be present when matplotlib is installed
        assert (diag_dir / "qq_plot.png").exists(), (
            "qq_plot.png should be produced by write_diagnostics() when matplotlib is available"
        )
        assert (diag_dir / "qq_plot.png").stat().st_size > 0

    def test_diagnostics_still_succeeds_without_matplotlib(self, minimal_results_df, tmp_path):
        """write_diagnostics() completes successfully even when matplotlib is absent.

        The TSV and TXT files are written regardless; only qq_plot.png is skipped.
        Uses mock.patch to replace write_qq_plot with a version that returns False,
        simulating the matplotlib-absent behavior without patching builtins.__import__.
        """
        diag_dir = tmp_path / "diagnostics_no_mpl"

        # Patch write_qq_plot to simulate matplotlib-absent behavior
        with patch(
            "variantcentrifuge.association.diagnostics.write_qq_plot",
            return_value=False,
        ):
            write_diagnostics(
                results_df=minimal_results_df,
                diagnostics_dir=diag_dir,
                test_names=["fisher"],
                n_cases=10,
                n_controls=10,
                cohort_warnings=[],
            )

        # TSV and TXT files should still be present (written before write_qq_plot call)
        assert (diag_dir / "lambda_gc.tsv").exists()
        assert (diag_dir / "qq_data.tsv").exists()
        assert (diag_dir / "summary.txt").exists()
        # PNG should NOT exist (write_qq_plot returned False and wrote nothing)
        assert not (diag_dir / "qq_plot.png").exists()
