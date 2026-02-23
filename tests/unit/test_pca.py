# File: tests/unit/test_pca.py
"""
Unit tests for variantcentrifuge/association/pca.py.

Covers:
- _detect_pca_format() on all three format families
- load_pca_file() with PLINK .eigenvec (with / without header)
- load_pca_file() with AKT stdout format (headerless, single ID column)
- load_pca_file() with generic TSV (header row)
- n_components slicing and "fewer than requested" warning
- >20 components warning
- Sample alignment in different order
- Missing sample raises ValueError
- merge_pca_covariates() with None covariate matrix
- merge_pca_covariates() with existing covariate matrix
"""

from __future__ import annotations

import textwrap

import numpy as np
import pytest

from variantcentrifuge.association.pca import (
    _detect_pca_format,
    load_pca_file,
    merge_pca_covariates,
)

# ---------------------------------------------------------------------------
# Helper: write a temp file from multi-line string
# ---------------------------------------------------------------------------


def _write(tmp_path, filename: str, content: str) -> str:
    p = tmp_path / filename
    p.write_text(textwrap.dedent(content))
    return str(p)


# ---------------------------------------------------------------------------
# _detect_pca_format tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestDetectPcaFormat:
    def test_plink_with_hash_fid_header(self):
        lines = ["#FID\tIID\tPC1\tPC2", "FAM1\tSAMPLE1\t0.1\t0.2"]
        assert _detect_pca_format(lines) == "plink_header"

    def test_plink_with_fid_header_no_hash(self):
        lines = ["FID\tIID\tPC1\tPC2", "FAM1\tSAMPLE1\t0.1\t0.2"]
        assert _detect_pca_format(lines) == "plink_header"

    def test_plink_fid_case_insensitive(self):
        lines = ["fid\tIID\tPC1\tPC2", "FAM1\tSAMPLE1\t0.1\t0.2"]
        assert _detect_pca_format(lines) == "plink_header"

    def test_plink_without_header(self):
        # Two non-numeric columns then numeric PCs
        lines = ["FAM1 SAMPLE1 0.1 0.2 0.3", "FAM2 SAMPLE2 0.4 0.5 0.6"]
        assert _detect_pca_format(lines) == "plink_nohdr"

    def test_akt_format(self):
        # Single ID column then numeric PCs (no header)
        lines = ["SAMPLE1 0.1 0.2 0.3", "SAMPLE2 0.4 0.5 0.6"]
        assert _detect_pca_format(lines) == "akt_or_generic"

    def test_empty_lines(self):
        assert _detect_pca_format([]) == "generic"

    def test_generic_fallback(self):
        # All numeric columns — fall through to generic
        lines = ["0.1 0.2 0.3", "0.4 0.5 0.6"]
        assert _detect_pca_format(lines) == "generic"


# ---------------------------------------------------------------------------
# load_pca_file: PLINK .eigenvec with header
# ---------------------------------------------------------------------------

PLINK_HDR_CONTENT = """\
    FID\tIID\tPC1\tPC2\tPC3
    FAM1\tSAMPLE_A\t0.11\t0.21\t0.31
    FAM1\tSAMPLE_B\t0.12\t0.22\t0.32
    FAM2\tSAMPLE_C\t0.13\t0.23\t0.33
"""

PLINK_HASH_HDR_CONTENT = """\
    #FID\tIID\tPC1\tPC2\tPC3
    FAM1\tSAMPLE_A\t0.11\t0.21\t0.31
    FAM1\tSAMPLE_B\t0.12\t0.22\t0.32
    FAM2\tSAMPLE_C\t0.13\t0.23\t0.33
"""

VCF_SAMPLES = ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]


@pytest.mark.unit
class TestLoadPcaFilePlinkHeader:
    def test_plink_header_basic(self, tmp_path):
        p = _write(tmp_path, "test.eigenvec", PLINK_HDR_CONTENT)
        matrix, col_names = load_pca_file(p, VCF_SAMPLES, n_components=3)
        assert matrix.shape == (3, 3)
        assert col_names == ["PC1", "PC2", "PC3"]

    def test_plink_hash_header_basic(self, tmp_path):
        p = _write(tmp_path, "test.eigenvec", PLINK_HASH_HDR_CONTENT)
        matrix, col_names = load_pca_file(p, VCF_SAMPLES, n_components=3)
        assert matrix.shape == (3, 3)
        assert col_names == ["PC1", "PC2", "PC3"]

    def test_plink_header_uses_iid_not_fid(self, tmp_path):
        """IID (column 2) should be the sample index, not FID (column 1)."""
        content = "FID\tIID\tPC1\nFAM1\tIID_SAMPLE_A\t0.11\nFAM2\tIID_SAMPLE_B\t0.22\n"
        p = _write(tmp_path, "test.eigenvec", content)
        matrix, _ = load_pca_file(p, ["IID_SAMPLE_A", "IID_SAMPLE_B"], n_components=1)
        assert matrix.shape == (2, 1)
        assert float(matrix[0, 0]) == pytest.approx(0.11)
        assert float(matrix[1, 0]) == pytest.approx(0.22)

    def test_plink_header_alignment(self, tmp_path):
        """Matrix rows must match vcf_samples order even if file order differs."""
        p = _write(tmp_path, "test.eigenvec", PLINK_HDR_CONTENT)
        vcf = ["SAMPLE_C", "SAMPLE_A", "SAMPLE_B"]
        matrix, _ = load_pca_file(p, vcf, n_components=1)
        # SAMPLE_C should come first in output
        assert float(matrix[0, 0]) == pytest.approx(0.13)
        assert float(matrix[1, 0]) == pytest.approx(0.11)
        assert float(matrix[2, 0]) == pytest.approx(0.12)


# ---------------------------------------------------------------------------
# load_pca_file: PLINK .eigenvec without header
# ---------------------------------------------------------------------------

PLINK_NOHDR_CONTENT = """\
    FAM1 SAMPLE_A 0.11 0.21 0.31
    FAM1 SAMPLE_B 0.12 0.22 0.32
    FAM2 SAMPLE_C 0.13 0.23 0.33
"""


@pytest.mark.unit
class TestLoadPcaFilePlinkNoHeader:
    def test_plink_nohdr_basic(self, tmp_path):
        p = _write(tmp_path, "test.eigenvec", PLINK_NOHDR_CONTENT)
        matrix, col_names = load_pca_file(p, VCF_SAMPLES, n_components=3)
        assert matrix.shape == (3, 3)
        assert col_names == ["PC1", "PC2", "PC3"]

    def test_plink_nohdr_uses_second_column(self, tmp_path):
        """Ensure IID (2nd column) is used, not FID (1st column)."""
        content = "FAM1 SAMPLE_A 0.55\nFAM2 SAMPLE_B 0.66\n"
        p = _write(tmp_path, "test.eigenvec", content)
        matrix, _ = load_pca_file(p, ["SAMPLE_A", "SAMPLE_B"], n_components=1)
        assert float(matrix[0, 0]) == pytest.approx(0.55)

    def test_plink_nohdr_alignment(self, tmp_path):
        p = _write(tmp_path, "test.eigenvec", PLINK_NOHDR_CONTENT)
        vcf = ["SAMPLE_C", "SAMPLE_B", "SAMPLE_A"]
        matrix, _ = load_pca_file(p, vcf, n_components=1)
        assert float(matrix[0, 0]) == pytest.approx(0.13)  # SAMPLE_C
        assert float(matrix[2, 0]) == pytest.approx(0.11)  # SAMPLE_A


# ---------------------------------------------------------------------------
# load_pca_file: AKT format (single ID column, no header)
# ---------------------------------------------------------------------------

AKT_CONTENT = """\
    SAMPLE_A 0.11 0.21 0.31
    SAMPLE_B 0.12 0.22 0.32
    SAMPLE_C 0.13 0.23 0.33
"""


@pytest.mark.unit
class TestLoadPcaFileAktFormat:
    def test_akt_basic(self, tmp_path):
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        matrix, col_names = load_pca_file(p, VCF_SAMPLES, n_components=3)
        assert matrix.shape == (3, 3)
        assert col_names == ["PC1", "PC2", "PC3"]

    def test_akt_correct_values(self, tmp_path):
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        matrix, _ = load_pca_file(p, VCF_SAMPLES, n_components=2)
        assert float(matrix[0, 0]) == pytest.approx(0.11)
        assert float(matrix[0, 1]) == pytest.approx(0.21)

    def test_akt_alignment_reorder(self, tmp_path):
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        vcf = ["SAMPLE_C", "SAMPLE_A", "SAMPLE_B"]
        matrix, _ = load_pca_file(p, vcf, n_components=1)
        assert float(matrix[0, 0]) == pytest.approx(0.13)  # SAMPLE_C first


# ---------------------------------------------------------------------------
# load_pca_file: generic TSV (header row)
# ---------------------------------------------------------------------------

GENERIC_TSV_CONTENT = """\
    sample_id\tPC1\tPC2\tPC3
    SAMPLE_A\t0.11\t0.21\t0.31
    SAMPLE_B\t0.12\t0.22\t0.32
    SAMPLE_C\t0.13\t0.23\t0.33
"""


@pytest.mark.unit
class TestLoadPcaFileGenericTsv:
    def test_generic_tsv_basic(self, tmp_path):
        p = _write(tmp_path, "pca.tsv", GENERIC_TSV_CONTENT)
        matrix, col_names = load_pca_file(p, VCF_SAMPLES, n_components=3)
        assert matrix.shape == (3, 3)
        assert col_names == ["PC1", "PC2", "PC3"]

    def test_generic_tsv_correct_values(self, tmp_path):
        p = _write(tmp_path, "pca.tsv", GENERIC_TSV_CONTENT)
        matrix, _ = load_pca_file(p, VCF_SAMPLES, n_components=1)
        assert float(matrix[2, 0]) == pytest.approx(0.13)  # SAMPLE_C


# ---------------------------------------------------------------------------
# n_components slicing
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestNComponentsSlicing:
    def test_slice_to_fewer_components(self, tmp_path):
        """Request 2 of 3 available -> shape is (n_samples, 2)."""
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        matrix, col_names = load_pca_file(p, VCF_SAMPLES, n_components=2)
        assert matrix.shape == (3, 2)
        assert col_names == ["PC1", "PC2"]

    def test_fewer_available_than_requested_warns(self, tmp_path, caplog):
        """Request 10, file has 3 -> use 3, emit warning."""
        import logging

        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            matrix, col_names = load_pca_file(p, VCF_SAMPLES, n_components=10)
        assert matrix.shape[1] == 3
        assert col_names == ["PC1", "PC2", "PC3"]
        assert any("only has" in r.message or "only" in r.message for r in caplog.records)

    def test_more_than_20_components_warns(self, tmp_path, caplog):
        """Requesting >20 PCA components logs a warning."""
        import logging

        # Build a file with 25 numeric columns
        header = "sample_id\t" + "\t".join(f"PC{i}" for i in range(1, 26))
        samples = ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
        rows = [
            f"{s}\t" + "\t".join(f"0.{j}{i:02d}" for i in range(1, 26))
            for j, s in enumerate(samples, 1)
        ]
        content = "\n".join([header, *rows]) + "\n"
        p = _write(tmp_path, "pca.tsv", content)
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            load_pca_file(p, samples, n_components=25)
        assert any(
            "25" in r.message or ">20" in r.message or "20" in r.message for r in caplog.records
        )


# ---------------------------------------------------------------------------
# Sample alignment and missing sample errors
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSampleAlignment:
    def test_missing_sample_raises_value_error(self, tmp_path):
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        with pytest.raises(ValueError, match="missing from PCA file"):
            load_pca_file(p, ["SAMPLE_A", "SAMPLE_MISSING"], n_components=1)

    def test_extra_samples_in_pca_file_are_ignored(self, tmp_path):
        """Extra samples in PCA file (not in VCF) are silently dropped."""
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        # Only request SAMPLE_A — SAMPLE_B and SAMPLE_C are extras
        matrix, _ = load_pca_file(p, ["SAMPLE_A"], n_components=1)
        assert matrix.shape == (1, 1)
        assert float(matrix[0, 0]) == pytest.approx(0.11)

    def test_alignment_matches_vcf_order(self, tmp_path):
        """When VCF order is reversed, matrix rows must follow VCF order."""
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        vcf = ["SAMPLE_C", "SAMPLE_B", "SAMPLE_A"]
        matrix, _ = load_pca_file(p, vcf, n_components=1)
        # File order: A=0.11, B=0.12, C=0.13
        # Expected output order: C, B, A
        assert float(matrix[0, 0]) == pytest.approx(0.13)
        assert float(matrix[1, 0]) == pytest.approx(0.12)
        assert float(matrix[2, 0]) == pytest.approx(0.11)


# ---------------------------------------------------------------------------
# merge_pca_covariates
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestMergePcaCovariates:
    def test_merge_with_none_covariate_matrix(self):
        """When covariate_matrix is None, return pca_matrix unchanged."""
        pca = np.array([[0.1, 0.2], [0.3, 0.4]])
        pca_names = ["PC1", "PC2"]
        merged, names = merge_pca_covariates(pca, pca_names, None, [])
        assert merged is pca
        assert names == ["PC1", "PC2"]

    def test_merge_with_existing_covariate_matrix(self):
        """PCA columns are appended after existing covariate columns."""
        cov = np.array([[1.0, 2.0], [3.0, 4.0]])
        cov_names = ["AGE", "SEX"]
        pca = np.array([[0.1, 0.2], [0.3, 0.4]])
        pca_names = ["PC1", "PC2"]
        merged, names = merge_pca_covariates(pca, pca_names, cov, cov_names)
        assert merged.shape == (2, 4)
        assert names == ["AGE", "SEX", "PC1", "PC2"]
        np.testing.assert_array_almost_equal(merged[:, :2], cov)
        np.testing.assert_array_almost_equal(merged[:, 2:], pca)

    def test_merge_shape_mismatch_raises(self):
        """Shape mismatch on axis 0 raises AssertionError."""
        cov = np.zeros((3, 2))
        pca = np.zeros((2, 1))
        with pytest.raises(AssertionError):
            merge_pca_covariates(pca, ["PC1"], cov, ["A", "B"])

    def test_merge_returns_float64(self):
        pca = np.array([[1, 2]], dtype=np.float32)
        merged, _ = merge_pca_covariates(pca, ["PC1", "PC2"], None, [])
        # pca itself is returned when covariate_matrix is None
        assert merged.dtype == np.float32  # returned as-is without copy

    def test_merge_col_names_concatenated(self):
        cov = np.zeros((2, 1))
        pca = np.zeros((2, 3))
        _, names = merge_pca_covariates(pca, ["PC1", "PC2", "PC3"], cov, ["COVAR1"])
        assert names == ["COVAR1", "PC1", "PC2", "PC3"]
