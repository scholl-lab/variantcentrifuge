"""
Unit tests for the association diagnostics module.

Tests lambda_GC calibration, QQ data computation, sample size warnings,
per-gene warnings, write_diagnostics file output, and stage-level integration
(config -> write_diagnostics -> files on disk).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.diagnostics import (
    compute_lambda_gc,
    compute_per_gene_warnings,
    compute_qq_data,
    emit_sample_size_warnings,
    write_diagnostics,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def default_config():
    """AssociationConfig with default thresholds."""
    return AssociationConfig(
        min_cases=200,
        max_case_control_ratio=20.0,
        min_case_carriers=10,
    )


# ---------------------------------------------------------------------------
# Lambda_GC tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestComputeLambdaGC:
    """Tests for compute_lambda_gc() genomic inflation factor computation."""

    def test_lambda_gc_uniform(self):
        """Lambda_GC on 1000 uniform(0,1) p-values should be within [0.85, 1.15].

        Under the null, the genomic inflation factor should be near 1.0.
        A generous tolerance [0.85, 1.15] accounts for random sampling variance.
        """
        rng = np.random.default_rng(seed=42)
        p_values = list(rng.uniform(0, 1, size=1000))
        lam = compute_lambda_gc(p_values)
        assert lam is not None
        assert 0.85 <= lam <= 1.15, f"Lambda_GC on null data was {lam:.4f}, expected near 1.0"

    def test_lambda_gc_inflated(self):
        """Lambda_GC on inflated chi2 statistics should return ~1.5.

        Simulate inflation by multiplying chi2(1) random variables by 1.5,
        then converting back to p-values. Lambda_GC should detect the inflation.
        """
        from scipy.stats import chi2 as chi2_dist

        rng = np.random.default_rng(seed=123)
        n = 1000
        inflation_factor = 1.5
        # Generate chi2(1) statistics with inflation
        chi2_vals = rng.chisquare(df=1, size=n) * inflation_factor
        # Convert to p-values
        p_values = list(chi2_dist.sf(chi2_vals, df=1))
        lam = compute_lambda_gc(p_values)
        assert lam is not None
        assert abs(lam - inflation_factor) < 0.15, (
            f"Lambda_GC={lam:.4f} should be near inflation factor {inflation_factor}"
        )

    def test_lambda_gc_empty_returns_none(self):
        """compute_lambda_gc([]) returns None (no valid p-values)."""
        result = compute_lambda_gc([])
        assert result is None

    def test_lambda_gc_single_returns_none(self):
        """compute_lambda_gc([0.5]) returns None (fewer than 2 valid p-values)."""
        result = compute_lambda_gc([0.5])
        assert result is None

    def test_lambda_gc_two_nones_returns_none(self):
        """compute_lambda_gc([None, None]) returns None."""
        result = compute_lambda_gc([None, None])
        assert result is None

    def test_lambda_gc_with_nones_filters_correctly(self):
        """Lambda_GC filters None values before computation.

        Result with None values interspersed should equal result on the
        non-None values only.
        """
        rng = np.random.default_rng(seed=99)
        clean_p = list(rng.uniform(0, 1, size=500))
        # Intersperse None values
        mixed_p: list[float | None] = []
        for i, p in enumerate(clean_p):
            mixed_p.append(p)
            if i % 5 == 0:
                mixed_p.append(None)

        lam_clean = compute_lambda_gc(clean_p)
        lam_mixed = compute_lambda_gc(mixed_p)

        assert lam_clean is not None
        assert lam_mixed is not None
        assert lam_clean == pytest.approx(lam_mixed, rel=1e-9)

    def test_lambda_gc_returns_float(self):
        """compute_lambda_gc returns a Python float (not numpy scalar)."""
        rng = np.random.default_rng(seed=7)
        p_values = list(rng.uniform(0, 1, size=200))
        lam = compute_lambda_gc(p_values)
        assert isinstance(lam, float)

    def test_lambda_gc_well_calibrated_range(self):
        """Lambda_GC on large null sample is within tighter range [0.9, 1.1]."""
        rng = np.random.default_rng(seed=2024)
        p_values = list(rng.uniform(0, 1, size=5000))
        lam = compute_lambda_gc(p_values)
        assert lam is not None
        assert 0.9 <= lam <= 1.1, f"Lambda_GC={lam:.4f} outside [0.9, 1.1] for n=5000 null"


# ---------------------------------------------------------------------------
# QQ data tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestComputeQQData:
    """Tests for compute_qq_data() observed vs expected -log10(p) data."""

    def test_qq_data_shape(self):
        """100 p-values produces DataFrame with exactly 100 rows."""
        rng = np.random.default_rng(seed=10)
        p_values = list(rng.uniform(0, 1, size=100))
        qq = compute_qq_data(p_values, test_name="fisher")
        assert qq.shape[0] == 100

    def test_qq_data_columns(self):
        """QQ DataFrame has exactly the required three columns."""
        qq = compute_qq_data([0.1, 0.05, 0.01], test_name="fisher")
        assert list(qq.columns) == ["test", "expected_neg_log10_p", "observed_neg_log10_p"]

    def test_qq_data_sort_order_ascending(self):
        """expected_neg_log10_p is sorted ascending (non-significant end first)."""
        rng = np.random.default_rng(seed=20)
        p_values = list(rng.uniform(0, 1, size=100))
        qq = compute_qq_data(p_values, test_name="test")
        assert qq["expected_neg_log10_p"].is_monotonic_increasing, (
            "expected_neg_log10_p should be sorted ascending (smallest value first)"
        )

    def test_qq_data_empty_input(self):
        """Empty input produces empty DataFrame with correct columns."""
        qq = compute_qq_data([], test_name="fisher")
        assert qq.empty
        assert list(qq.columns) == ["test", "expected_neg_log10_p", "observed_neg_log10_p"]

    def test_qq_data_hazen_quantile_formula(self):
        """Expected quantile values match the Hazen formula i/(n+1).

        For n=100 p-values, sorted ascending by expected_neg_log10_p:
        - First row: expected_p = n/(n+1) -> expected_neg_log10_p = -log10(n/(n+1))
        - Last row: expected_p = 1/(n+1) -> expected_neg_log10_p = -log10(1/(n+1))
        """
        n = 100
        rng = np.random.default_rng(seed=30)
        p_values = list(rng.uniform(0, 1, size=n))
        qq = compute_qq_data(p_values, test_name="fisher")

        # Sorted ascending: first = smallest expected_neg_log10_p = -log10(n/(n+1))
        expected_first = -np.log10(n / (n + 1))
        expected_last = -np.log10(1 / (n + 1))

        assert qq.iloc[0]["expected_neg_log10_p"] == pytest.approx(expected_first, rel=1e-9)
        assert qq.iloc[-1]["expected_neg_log10_p"] == pytest.approx(expected_last, rel=1e-9)

    def test_qq_data_test_name_column(self):
        """test column contains the provided test_name for all rows."""
        qq = compute_qq_data([0.1, 0.05, 0.01], test_name="skat_python")
        assert (qq["test"] == "skat_python").all()

    def test_qq_data_nones_excluded(self):
        """None values in input are excluded from the QQ computation."""
        p_clean = [0.1, 0.05, 0.01]
        p_with_none = [0.1, None, 0.05, None, 0.01]

        qq_clean = compute_qq_data(p_clean, test_name="fisher")
        qq_mixed = compute_qq_data(p_with_none, test_name="fisher")

        # Both should have same number of rows (3)
        assert qq_clean.shape == qq_mixed.shape

    def test_qq_data_observed_values_finite(self):
        """All observed_neg_log10_p values should be finite (no +inf from p=0)."""
        import math

        # Include very small p-value to test clipping
        p_values = [1e-200, 0.05, 0.1]
        qq = compute_qq_data(p_values, test_name="fisher")
        for val in qq["observed_neg_log10_p"]:
            assert math.isfinite(val), f"Non-finite observed_neg_log10_p: {val}"


# ---------------------------------------------------------------------------
# Sample size warning tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestEmitSampleSizeWarnings:
    """Tests for emit_sample_size_warnings() cohort-level warning detection."""

    def test_warnings_low_cases(self, default_config):
        """n_cases=100 below min_cases=200 triggers LOW_CASE_COUNT warning."""
        warnings = emit_sample_size_warnings(100, 200, default_config)
        assert "LOW_CASE_COUNT" in warnings

    def test_warnings_imbalanced_cohort(self, default_config):
        """n_controls/n_cases > max_case_control_ratio triggers IMBALANCED_COHORT."""
        # n_controls/n_cases = 5000/200 = 25 > threshold of 20
        warnings = emit_sample_size_warnings(200, 5000, default_config)
        assert "IMBALANCED_COHORT" in warnings

    def test_warnings_none_when_adequate(self, default_config):
        """Adequate sample sizes with balanced ratio produce no warnings."""
        # n_cases=300 >= 200, ratio = 600/300 = 2.0 < 20.0
        warnings = emit_sample_size_warnings(300, 600, default_config)
        assert len(warnings) == 0

    def test_warnings_both_triggered(self, default_config):
        """Both LOW_CASE_COUNT and IMBALANCED_COHORT can trigger simultaneously."""
        # n_cases=50 < 200 (LOW_CASE_COUNT), 50 * 25 = 1250 controls -> ratio=25 (IMBALANCED)
        warnings = emit_sample_size_warnings(50, 1250, default_config)
        assert "LOW_CASE_COUNT" in warnings
        assert "IMBALANCED_COHORT" in warnings

    def test_warnings_returns_list(self, default_config):
        """emit_sample_size_warnings always returns a list."""
        result = emit_sample_size_warnings(300, 300, default_config)
        assert isinstance(result, list)

    def test_warnings_exactly_at_threshold_no_warning(self, default_config):
        """n_cases exactly at min_cases=200 should not trigger LOW_CASE_COUNT."""
        warnings = emit_sample_size_warnings(200, 200, default_config)
        assert "LOW_CASE_COUNT" not in warnings

    def test_warnings_custom_min_cases(self):
        """Custom min_cases threshold is respected."""
        cfg = AssociationConfig(min_cases=50, max_case_control_ratio=20.0, min_case_carriers=10)
        warnings = emit_sample_size_warnings(40, 100, cfg)
        assert "LOW_CASE_COUNT" in warnings

        warnings2 = emit_sample_size_warnings(60, 100, cfg)
        assert "LOW_CASE_COUNT" not in warnings2


# ---------------------------------------------------------------------------
# Per-gene warning tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestComputePerGeneWarnings:
    """Tests for compute_per_gene_warnings() carrier count warning detection."""

    def test_per_gene_warnings_low_carriers(self, default_config):
        """case_carriers=5 below min_case_carriers=10 triggers LOW_CARRIER_COUNT."""
        warnings = compute_per_gene_warnings("BRCA1", 5, default_config)
        assert "LOW_CARRIER_COUNT" in warnings

    def test_per_gene_warnings_none_when_adequate(self, default_config):
        """case_carriers >= min_case_carriers produces no warnings."""
        warnings = compute_per_gene_warnings("BRCA1", 15, default_config)
        assert len(warnings) == 0

    def test_per_gene_warnings_exactly_at_threshold_no_warning(self, default_config):
        """case_carriers exactly at min_case_carriers should not trigger warning."""
        warnings = compute_per_gene_warnings("BRCA1", 10, default_config)
        assert "LOW_CARRIER_COUNT" not in warnings

    def test_per_gene_warnings_returns_list(self, default_config):
        """compute_per_gene_warnings always returns a list."""
        result = compute_per_gene_warnings("GENE1", 20, default_config)
        assert isinstance(result, list)

    def test_per_gene_warnings_zero_carriers(self, default_config):
        """Zero case carriers triggers LOW_CARRIER_COUNT."""
        warnings = compute_per_gene_warnings("BRCA1", 0, default_config)
        assert "LOW_CARRIER_COUNT" in warnings

    def test_per_gene_warnings_custom_threshold(self):
        """Custom min_case_carriers threshold is respected."""
        cfg = AssociationConfig(min_cases=200, max_case_control_ratio=20.0, min_case_carriers=5)
        warnings = compute_per_gene_warnings("GENE1", 3, cfg)
        assert "LOW_CARRIER_COUNT" in warnings

        warnings2 = compute_per_gene_warnings("GENE1", 6, cfg)
        assert "LOW_CARRIER_COUNT" not in warnings2


# ---------------------------------------------------------------------------
# write_diagnostics tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestWriteDiagnostics:
    """Tests for write_diagnostics() file output."""

    def _make_results_df(self, n_genes: int = 5) -> pd.DataFrame:
        """Build a minimal results DataFrame for write_diagnostics() testing."""
        rng = np.random.default_rng(seed=42)
        genes = [f"GENE{i}" for i in range(n_genes)]
        fisher_p = list(rng.uniform(0.001, 0.9, size=n_genes))
        acat_p = list(rng.uniform(0.001, 0.9, size=n_genes))
        return pd.DataFrame({"gene": genes, "fisher_p_value": fisher_p, "acat_o_p_value": acat_p})

    def test_write_diagnostics_creates_files(self, tmp_path):
        """write_diagnostics creates all three output files."""
        results_df = self._make_results_df(n_genes=10)
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path / "diag",
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )
        diag_dir = tmp_path / "diag"
        assert (diag_dir / "lambda_gc.tsv").exists()
        assert (diag_dir / "qq_data.tsv").exists()
        assert (diag_dir / "summary.txt").exists()

    def test_write_diagnostics_creates_directory(self, tmp_path):
        """write_diagnostics creates the directory if it doesn't exist."""
        results_df = self._make_results_df(n_genes=5)
        nested_dir = tmp_path / "a" / "b" / "c"
        assert not nested_dir.exists()
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=nested_dir,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )
        assert nested_dir.exists()

    def test_lambda_gc_tsv_columns(self, tmp_path):
        """lambda_gc.tsv has correct columns: test_name, lambda_gc, n_tests."""
        results_df = self._make_results_df(n_genes=10)
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )
        lambda_df = pd.read_csv(tmp_path / "lambda_gc.tsv", sep="\t")
        assert set(lambda_df.columns) == {"test_name", "lambda_gc", "n_tests"}

    def test_lambda_gc_tsv_has_rows_for_each_test(self, tmp_path):
        """lambda_gc.tsv has rows for each primary test and ACAT-O."""
        results_df = self._make_results_df(n_genes=10)
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )
        lambda_df = pd.read_csv(tmp_path / "lambda_gc.tsv", sep="\t")
        test_names_in_file = set(lambda_df["test_name"])
        assert "fisher" in test_names_in_file
        assert "acat_o" in test_names_in_file

    def test_qq_data_tsv_content(self, tmp_path):
        """qq_data.tsv has rows for each p-value from each test."""
        n_genes = 10
        results_df = self._make_results_df(n_genes=n_genes)
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )
        qq_df = pd.read_csv(tmp_path / "qq_data.tsv", sep="\t")
        assert not qq_df.empty
        # Should have rows for fisher and acat_o (n_genes each = 2*n_genes)
        assert len(qq_df) == 2 * n_genes

    def test_summary_txt_content(self, tmp_path):
        """summary.txt contains sample sizes and lambda_GC values."""
        results_df = self._make_results_df(n_genes=10)
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path,
            test_names=["fisher"],
            n_cases=250,
            n_controls=400,
            cohort_warnings=[],
        )
        summary_text = (tmp_path / "summary.txt").read_text()
        assert "n_cases" in summary_text
        assert "n_controls" in summary_text
        assert "250" in summary_text
        assert "400" in summary_text

    def test_summary_txt_contains_lambda_gc(self, tmp_path):
        """summary.txt contains lambda_GC section header."""
        results_df = self._make_results_df(n_genes=10)
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )
        summary_text = (tmp_path / "summary.txt").read_text()
        assert "Lambda_GC" in summary_text or "lambda_GC" in summary_text.lower()

    def test_summary_txt_shows_cohort_warnings(self, tmp_path):
        """summary.txt lists cohort-level warnings when present."""
        results_df = self._make_results_df(n_genes=10)
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=["LOW_CASE_COUNT", "IMBALANCED_COHORT"],
        )
        summary_text = (tmp_path / "summary.txt").read_text()
        assert "LOW_CASE_COUNT" in summary_text
        assert "IMBALANCED_COHORT" in summary_text

    def test_write_diagnostics_with_warnings_column(self, tmp_path):
        """summary.txt reports gene warning count when 'warnings' column present."""
        n_genes = 5
        results_df = self._make_results_df(n_genes=n_genes)
        results_df["warnings"] = ["LOW_CARRIER_COUNT", "", "LOW_CARRIER_COUNT", "", ""]
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=tmp_path,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )
        summary_text = (tmp_path / "summary.txt").read_text()
        # 2 genes with warnings should be reported
        assert "Genes with per-gene warnings: 2" in summary_text


# ---------------------------------------------------------------------------
# Stage-level diagnostics integration test
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestStageDiagnosticsIntegration:
    """Tests for the full wiring from config -> write_diagnostics -> files on disk.

    These tests simulate what AssociationAnalysisStage does when diagnostics_output
    is configured: they call write_diagnostics() with a realistic results DataFrame
    and verify all expected files are created with correct content.

    This validates the wiring path that the stage will use without needing to
    instantiate the full stage or PipelineContext.
    """

    def test_stage_diagnostics_integration(self, tmp_path):
        """Full path: config with diagnostics_output -> write_diagnostics -> files on disk.

        Simulates the exact sequence AssociationAnalysisStage uses:
        1. Read config.diagnostics_output path
        2. Build results_df from engine output (with fisher_p_value + acat_o_p_value)
        3. Call write_diagnostics() with all required parameters
        4. Assert all three files exist in the specified directory
        5. Assert file contents are meaningful
        """
        # Step 1: Create config as the stage would read it
        config = AssociationConfig(
            min_cases=200,
            max_case_control_ratio=20.0,
            min_case_carriers=10,
            diagnostics_output=str(tmp_path / "diag"),
        )

        # Step 2: Build a mock results_df from engine output (at least 2 genes)
        rng = np.random.default_rng(seed=2024)
        n_genes = 5
        genes = ["BRCA1", "TP53", "MYH7", "LMNA", "RYR2"]
        fisher_p = list(rng.uniform(0.001, 0.9, size=n_genes))
        acat_p = list(rng.uniform(0.001, 0.9, size=n_genes))
        results_df = pd.DataFrame(
            {
                "gene": genes,
                "n_cases": [300] * n_genes,
                "n_controls": [300] * n_genes,
                "n_variants": [3, 2, 5, 1, 4],
                "fisher_p_value": fisher_p,
                "fisher_or": [1.5, 2.0, 0.8, 1.2, 3.1],
                "acat_o_p_value": acat_p,
                "acat_o_corrected_p_value": list(rng.uniform(0.01, 0.9, size=n_genes)),
            }
        )

        # Step 3: Call write_diagnostics as the stage would
        assert config.diagnostics_output is not None  # guard for type checker
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=config.diagnostics_output,
            test_names=["fisher"],
            n_cases=300,
            n_controls=300,
            cohort_warnings=[],
        )

        # Step 4: Assert directory was created
        diag_dir = tmp_path / "diag"
        assert diag_dir.exists(), "Diagnostics directory was not created"

        # Step 5: Assert all three files exist
        assert (diag_dir / "lambda_gc.tsv").exists(), "lambda_gc.tsv not created"
        assert (diag_dir / "qq_data.tsv").exists(), "qq_data.tsv not created"
        assert (diag_dir / "summary.txt").exists(), "summary.txt not created"

        # Step 6: Assert lambda_gc.tsv has rows for "fisher" and "acat_o"
        lambda_df = pd.read_csv(diag_dir / "lambda_gc.tsv", sep="\t")
        test_names_in_file = set(lambda_df["test_name"])
        assert "fisher" in test_names_in_file, (
            f"'fisher' not in lambda_gc.tsv test_name column: {test_names_in_file}"
        )
        assert "acat_o" in test_names_in_file, (
            f"'acat_o' not in lambda_gc.tsv test_name column: {test_names_in_file}"
        )

        # Step 7: Assert qq_data.tsv has rows (at least n_genes rows for "fisher" alone)
        qq_df = pd.read_csv(diag_dir / "qq_data.tsv", sep="\t")
        assert not qq_df.empty, "qq_data.tsv is empty"
        assert len(qq_df) >= n_genes, (
            f"qq_data.tsv has {len(qq_df)} rows, expected at least {n_genes}"
        )

        # Step 8: Assert summary.txt contains "n_cases" (key content marker)
        summary_text = (diag_dir / "summary.txt").read_text()
        assert "n_cases" in summary_text, "summary.txt does not contain 'n_cases'"
        assert "300" in summary_text, "summary.txt does not contain n_cases value '300'"

    def test_stage_diagnostics_integration_with_per_gene_warnings(self, tmp_path):
        """Integration test with warnings column in results_df.

        Verifies that per-gene warnings from the engine output propagate correctly
        to the summary.txt file when the 'warnings' column is present.
        """
        config = AssociationConfig(
            diagnostics_output=str(tmp_path / "diag2"),
        )

        rng = np.random.default_rng(seed=99)
        n_genes = 4
        genes = ["BRCA1", "TP53", "MYH7", "LMNA"]
        results_df = pd.DataFrame(
            {
                "gene": genes,
                "fisher_p_value": list(rng.uniform(0.001, 0.9, size=n_genes)),
                "acat_o_p_value": list(rng.uniform(0.001, 0.9, size=n_genes)),
                # Two genes with warnings, two without
                "warnings": ["LOW_CARRIER_COUNT", "", "LOW_CARRIER_COUNT", ""],
            }
        )

        assert config.diagnostics_output is not None
        write_diagnostics(
            results_df=results_df,
            diagnostics_dir=config.diagnostics_output,
            test_names=["fisher"],
            n_cases=250,
            n_controls=500,
            cohort_warnings=[],
        )

        diag_dir = tmp_path / "diag2"
        summary_text = (diag_dir / "summary.txt").read_text()
        assert "Genes with per-gene warnings: 2" in summary_text
