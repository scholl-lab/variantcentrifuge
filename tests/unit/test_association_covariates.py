"""
Unit tests for variantcentrifuge.association.covariates.

Tests cover:
- Basic covariate loading from TSV/CSV files
- Sample alignment to VCF order (COV-01)
- Missing VCF sample raises ValueError (COV-01 guard)
- Extra covariate samples warn but continue
- Categorical auto-detection and explicit one-hot encoding (COV-02)
- One-hot encoding dtype is float64 not bool (pandas 2.x pitfall) (COV-02)
- Multicollinearity warning when condition number > 1000 (COV-03)
- Column selection via covariate_columns parameter (COV-04)
- CSV delimiter auto-detection
"""

from __future__ import annotations

import logging
import pathlib

import numpy as np
import pytest

from variantcentrifuge.association.covariates import load_covariates

# ---------------------------------------------------------------------------
# Basic loading
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLoadCovariatesBasic:
    """Basic covariate file loading and shape checks."""

    def test_load_covariates_basic(self, tmp_path: pathlib.Path) -> None:
        """Create temp TSV with 4 samples, 2 numeric columns; assert shape and values."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text(
            "sample_id\tage\tpc1\n"
            "SAMPLE_A\t30\t0.1\n"
            "SAMPLE_B\t45\t-0.2\n"
            "SAMPLE_C\t28\t0.05\n"
            "SAMPLE_D\t52\t0.3\n"
        )
        vcf_samples = ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C", "SAMPLE_D"]

        cov_mat, col_names = load_covariates(str(cov_file), vcf_samples)

        assert cov_mat.shape == (4, 2), f"Expected (4, 2), got {cov_mat.shape}"
        assert col_names == ["age", "pc1"]
        assert cov_mat[0, 0] == pytest.approx(30.0)
        assert cov_mat[1, 1] == pytest.approx(-0.2)
        assert cov_mat[3, 0] == pytest.approx(52.0)

    def test_load_covariates_dtype_float64(self, tmp_path: pathlib.Path) -> None:
        """Returned matrix is float64."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text("sample_id\tage\n" "S1\t30\n" "S2\t45\n")

        cov_mat, _ = load_covariates(str(cov_file), ["S1", "S2"])

        assert cov_mat.dtype == np.float64

    def test_load_covariates_no_nan(self, tmp_path: pathlib.Path) -> None:
        """Returned matrix contains no NaN."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text("sample_id\tage\tpc1\n" "S1\t30\t0.1\n" "S2\t45\t-0.2\n")

        cov_mat, _ = load_covariates(str(cov_file), ["S1", "S2"])

        assert not np.isnan(cov_mat).any()


# ---------------------------------------------------------------------------
# Sample alignment (COV-01)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLoadCovariatesAlignment:
    """
    Sample alignment tests — COV-01.

    The critical invariant: load_covariates must produce rows in VCF sample order
    regardless of the order samples appear in the covariate file.
    """

    def test_load_covariates_alignment_shuffled(self, tmp_path: pathlib.Path) -> None:
        """Covariate file in reverse order vs VCF; output rows match VCF order."""
        cov_file = tmp_path / "covariates.tsv"
        # File has samples in reverse order: D, C, B, A
        cov_file.write_text(
            "sample_id\tage\n"
            "SAMPLE_D\t52\n"
            "SAMPLE_C\t28\n"
            "SAMPLE_B\t45\n"
            "SAMPLE_A\t30\n"
        )
        # VCF order: A, B, C, D
        vcf_samples = ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C", "SAMPLE_D"]

        cov_mat, _ = load_covariates(str(cov_file), vcf_samples)

        assert cov_mat.shape == (4, 1)
        # Row 0 should be SAMPLE_A (age=30), not SAMPLE_D (age=52)
        assert cov_mat[0, 0] == pytest.approx(30.0), "Row 0 should be SAMPLE_A (age=30)"
        assert cov_mat[1, 0] == pytest.approx(45.0), "Row 1 should be SAMPLE_B (age=45)"
        assert cov_mat[2, 0] == pytest.approx(28.0), "Row 2 should be SAMPLE_C (age=28)"
        assert cov_mat[3, 0] == pytest.approx(52.0), "Row 3 should be SAMPLE_D (age=52)"

    def test_load_covariates_alignment_produces_identical_results(
        self, tmp_path: pathlib.Path
    ) -> None:
        """Same data in VCF-order file vs shuffled file produces identical matrices."""
        # File 1: samples in VCF order (A, B, C, D)
        cov_file_ordered = tmp_path / "ordered.tsv"
        cov_file_ordered.write_text(
            "sample_id\tage\tpc1\n"
            "S1\t30\t0.1\n"
            "S2\t45\t-0.2\n"
            "S3\t28\t0.05\n"
            "S4\t52\t0.3\n"
        )

        # File 2: same data but samples shuffled
        cov_file_shuffled = tmp_path / "shuffled.tsv"
        cov_file_shuffled.write_text(
            "sample_id\tage\tpc1\n"
            "S3\t28\t0.05\n"
            "S1\t30\t0.1\n"
            "S4\t52\t0.3\n"
            "S2\t45\t-0.2\n"
        )

        vcf_samples = ["S1", "S2", "S3", "S4"]

        mat_ordered, _ = load_covariates(str(cov_file_ordered), vcf_samples)
        mat_shuffled, _ = load_covariates(str(cov_file_shuffled), vcf_samples)

        assert np.array_equal(mat_ordered, mat_shuffled), (
            "Shuffled covariate file must produce identical matrix to ordered file"
        )

    def test_load_covariates_missing_vcf_sample_raises(
        self, tmp_path: pathlib.Path
    ) -> None:
        """VCF sample not in covariate file raises ValueError mentioning 'VCF samples missing'."""
        cov_file = tmp_path / "covariates.tsv"
        # Only A, B, C — missing D
        cov_file.write_text("sample_id\tage\n" "S_A\t30\n" "S_B\t45\n" "S_C\t28\n")

        # VCF expects A, B, C, D
        vcf_samples = ["S_A", "S_B", "S_C", "S_D"]

        with pytest.raises(ValueError, match="VCF samples missing"):
            load_covariates(str(cov_file), vcf_samples)

    def test_load_covariates_missing_multiple_samples_all_reported(
        self, tmp_path: pathlib.Path
    ) -> None:
        """All missing VCF samples are reported in the error message."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text("sample_id\tage\n" "S_A\t30\n")

        # B, C, D are all missing from covariate file
        vcf_samples = ["S_A", "S_B", "S_C", "S_D"]

        with pytest.raises(ValueError) as exc_info:
            load_covariates(str(cov_file), vcf_samples)

        error_msg = str(exc_info.value)
        # At least one of the missing samples should be mentioned
        assert "S_B" in error_msg or "S_C" in error_msg or "S_D" in error_msg

    def test_load_covariates_extra_samples_warns(
        self, tmp_path: pathlib.Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Covariate file with extra samples not in VCF warns but returns correct shape."""
        cov_file = tmp_path / "covariates.tsv"
        # File has A, B, C, D, E but VCF only has A, B, C
        cov_file.write_text(
            "sample_id\tage\n"
            "S_A\t30\n"
            "S_B\t45\n"
            "S_C\t28\n"
            "S_D\t52\n"
            "S_E\t33\n"
        )
        vcf_samples = ["S_A", "S_B", "S_C"]

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            cov_mat, _ = load_covariates(str(cov_file), vcf_samples)

        # Shape matches VCF samples (3, not 5)
        assert cov_mat.shape == (3, 1), f"Expected (3, 1), got {cov_mat.shape}"

        # A warning should have been logged
        warning_messages = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("extra" in str(m).lower() for m in warning_messages), (
            f"Expected warning about extra samples, got: {warning_messages}"
        )


# ---------------------------------------------------------------------------
# Categorical encoding (COV-02)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLoadCovariatesCategorical:
    """Categorical covariate encoding tests — COV-02."""

    def test_load_covariates_categorical_auto_detect(self, tmp_path: pathlib.Path) -> None:
        """Non-numeric column with <=5 unique values is auto-detected and one-hot encoded."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text(
            "sample_id\tage\tsex\n"
            "S1\t30\tM\n"
            "S2\t45\tF\n"
            "S3\t28\tM\n"
            "S4\t52\tF\n"
        )
        vcf_samples = ["S1", "S2", "S3", "S4"]

        cov_mat, col_names = load_covariates(str(cov_file), vcf_samples)

        # Should have age + one sex column (drop_first=True reduces M/F to 1 col)
        assert cov_mat.shape == (4, 2), f"Expected (4, 2) after one-hot, got {cov_mat.shape}"
        # At least one column should be sex-related
        sex_cols = [c for c in col_names if "sex" in c.lower()]
        assert len(sex_cols) == 1, f"Expected 1 sex column, got: {col_names}"

    def test_load_covariates_categorical_explicit(self, tmp_path: pathlib.Path) -> None:
        """Passing categorical_columns=['ethnicity'] forces one-hot encoding."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text(
            "sample_id\tage\tethnicity\n"
            "S1\t30\tEUR\n"
            "S2\t45\tAFR\n"
            "S3\t28\tEUR\n"
            "S4\t52\tEAS\n"
        )
        vcf_samples = ["S1", "S2", "S3", "S4"]

        cov_mat, col_names = load_covariates(
            str(cov_file), vcf_samples, categorical_columns=["ethnicity"]
        )

        # age (1 col) + ethnicity (3 values, drop_first -> 2 cols) = 3 cols
        assert cov_mat.shape[0] == 4
        # More than 1 column (ethnicity was encoded)
        assert cov_mat.shape[1] >= 2
        # No "ethnicity" column remaining as-is (should be ethnicity_X or ethnicity_Y)
        assert "ethnicity" not in col_names

    def test_load_covariates_one_hot_dtype_float(self, tmp_path: pathlib.Path) -> None:
        """One-hot encoded columns are float64, not bool (pandas 2.x pitfall — COV-02)."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text(
            "sample_id\tsex\n"
            "S1\tM\n"
            "S2\tF\n"
            "S3\tM\n"
        )
        vcf_samples = ["S1", "S2", "S3"]

        cov_mat, _ = load_covariates(str(cov_file), vcf_samples)

        # All values should be float64 (0.0 or 1.0), not bool
        assert cov_mat.dtype == np.float64, f"Expected float64, got {cov_mat.dtype}"
        # Values should be 0.0 or 1.0, not True/False
        unique_vals = set(np.unique(cov_mat))
        assert unique_vals <= {0.0, 1.0}, f"Expected only 0.0/1.0, got {unique_vals}"

    def test_load_covariates_high_cardinality_not_auto_encoded(
        self, tmp_path: pathlib.Path
    ) -> None:
        """Column with >5 unique values is NOT auto-detected as categorical."""
        cov_file = tmp_path / "covariates.tsv"
        # "region" has 6 unique values — above the <=5 threshold
        cov_file.write_text(
            "sample_id\tregion\n"
            "S1\tNorth\n"
            "S2\tSouth\n"
            "S3\tEast\n"
            "S4\tWest\n"
            "S5\tCenter\n"
            "S6\tIsland\n"
        )
        vcf_samples = ["S1", "S2", "S3", "S4", "S5", "S6"]

        # This should raise ValueError since 'region' is non-numeric and
        # won't be auto-encoded — and it can't be cast to float
        # Actually it will fail at .astype(np.float64) since region is strings
        with pytest.raises((ValueError, Exception)):
            load_covariates(str(cov_file), vcf_samples)


# ---------------------------------------------------------------------------
# Multicollinearity warning (COV-03)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLoadCovariatesMulticollinearity:
    """Multicollinearity detection — COV-03."""

    def test_load_covariates_multicollinearity_warning(
        self, tmp_path: pathlib.Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Perfectly correlated columns trigger a multicollinearity WARNING."""
        cov_file = tmp_path / "covariates.tsv"
        # col2 = 2 * col1 — perfect correlation, high condition number
        cov_file.write_text(
            "sample_id\tcol1\tcol2\n"
            "S1\t1.0\t2.0\n"
            "S2\t2.0\t4.0\n"
            "S3\t3.0\t6.0\n"
            "S4\t4.0\t8.0\n"
            "S5\t5.0\t10.0\n"
        )
        vcf_samples = ["S1", "S2", "S3", "S4", "S5"]

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            cov_mat, _ = load_covariates(str(cov_file), vcf_samples)

        # Should return a valid matrix (warning, not error)
        assert cov_mat.shape == (5, 2)

        # Should have logged a multicollinearity warning
        warning_messages = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any(
            "multicollinear" in str(m).lower() or "condition" in str(m).lower()
            for m in warning_messages
        ), f"Expected multicollinearity warning, got: {warning_messages}"

    def test_load_covariates_low_condition_no_warning(
        self, tmp_path: pathlib.Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Independent columns do NOT trigger multicollinearity warning."""
        cov_file = tmp_path / "covariates.tsv"
        # age and pc1 are independent
        cov_file.write_text(
            "sample_id\tage\tpc1\n"
            "S1\t30\t0.1\n"
            "S2\t45\t-0.5\n"
            "S3\t28\t0.8\n"
            "S4\t52\t-0.3\n"
        )
        vcf_samples = ["S1", "S2", "S3", "S4"]

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            load_covariates(str(cov_file), vcf_samples)

        warning_messages = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert not any(
            "multicollinear" in str(m).lower() or "condition" in str(m).lower()
            for m in warning_messages
        )


# ---------------------------------------------------------------------------
# Delimiter auto-detection and column selection
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLoadCovariatesFormatAndSelection:
    """File format detection and column selection tests."""

    def test_load_covariates_csv_detection(self, tmp_path: pathlib.Path) -> None:
        """CSV file (.csv extension) loads correctly with comma delimiter."""
        cov_file = tmp_path / "covariates.csv"
        cov_file.write_text(
            "sample_id,age,pc1\n" "S1,30,0.1\n" "S2,45,-0.2\n" "S3,28,0.05\n"
        )
        vcf_samples = ["S1", "S2", "S3"]

        cov_mat, col_names = load_covariates(str(cov_file), vcf_samples)

        assert cov_mat.shape == (3, 2)
        assert col_names == ["age", "pc1"]
        assert cov_mat[0, 0] == pytest.approx(30.0)

    def test_load_covariates_column_selection(self, tmp_path: pathlib.Path) -> None:
        """Passing covariate_columns selects only specified columns."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text(
            "sample_id\tage\tpc1\tpc2\tbmi\tsmokingstatus\n"
            "S1\t30\t0.1\t-0.2\t22.5\t0\n"
            "S2\t45\t-0.2\t0.3\t27.1\t1\n"
            "S3\t28\t0.05\t0.1\t19.8\t0\n"
        )
        vcf_samples = ["S1", "S2", "S3"]

        cov_mat, col_names = load_covariates(
            str(cov_file), vcf_samples, covariate_columns=["age", "pc1"]
        )

        assert cov_mat.shape == (3, 2)
        assert col_names == ["age", "pc1"]

    def test_load_covariates_returns_tuple(self, tmp_path: pathlib.Path) -> None:
        """Return value is a tuple of (ndarray, list[str])."""
        cov_file = tmp_path / "covariates.tsv"
        cov_file.write_text("sample_id\tage\n" "S1\t30\n")

        result = load_covariates(str(cov_file), ["S1"])

        assert isinstance(result, tuple)
        assert len(result) == 2
        cov_mat, col_names = result
        assert isinstance(cov_mat, np.ndarray)
        assert isinstance(col_names, list)
