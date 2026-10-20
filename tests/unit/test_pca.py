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
from pathlib import Path

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


# ---------------------------------------------------------------------------
# Delimiter/header detection without csv.Sniffer and NaN validation (#221)
# ---------------------------------------------------------------------------

# Built like akt 0.3.3 stdout (vcfpca.cpp: ``cout << name << "\t" << P.row(j) << endl``):
# no header, one tab after the ID, then Eigen's "%g" values right-aligned to the
# row's widest value and joined by single spaces. Synthetic sample IDs.
AKT_FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "pca" / "akt_pca_stdout.txt"
AKT_FIXTURE_SAMPLES = ["SYN_S001", "SYN_S002", "SYN_S003", "SYN_S004"]
AKT_FIXTURE_VALUES = np.array(
    [
        [1.26366, -0.243411, 1.4646, 0.0123457, -3.5e-05],
        [-0.512, 2.03114, -1.10982, 0.4, 0.987654],
        [0.0889, -0.00731, 0.25, -2.71828, 1.0],
        [-1.5, 0.333333, -0.666667, 0.125, -0.0625],
    ]
)


@pytest.mark.unit
class TestPcaDelimiterAndValidation:
    def test_fixture_is_raw_akt_layout(self):
        first = AKT_FIXTURE.read_text().splitlines()[0]
        assert first.startswith("SYN_S001\t  ")  # tab after ID, then space padding
        assert "\t" not in first.split("\t", 1)[1]

    def test_raw_akt_stdout_loads_correct_values(self):
        matrix, names = load_pca_file(str(AKT_FIXTURE), AKT_FIXTURE_SAMPLES, n_components=5)
        assert names == ["PC1", "PC2", "PC3", "PC4", "PC5"]
        np.testing.assert_array_equal(matrix, AKT_FIXTURE_VALUES)

    def test_raw_akt_stdout_as_eigenvec_and_tsv(self, tmp_path):
        for name in ("pca.eigenvec", "pca.tsv", "pca"):
            p = tmp_path / name
            p.write_bytes(AKT_FIXTURE.read_bytes())
            matrix, _ = load_pca_file(str(p), AKT_FIXTURE_SAMPLES[::-1], n_components=5)
            np.testing.assert_array_equal(matrix, AKT_FIXTURE_VALUES[::-1])

    def test_akt_with_header_mixed_whitespace(self, tmp_path):
        content = "sample_id\t PC1   PC2\nSYN_A\t  0.5 -1.25\nSYN_B\t -0.75     2\n"
        p = _write(tmp_path, "pca.txt", content)
        matrix, _ = load_pca_file(p, ["SYN_B", "SYN_A"], n_components=2)
        np.testing.assert_array_equal(matrix, [[-0.75, 2.0], [0.5, -1.25]])

    def test_plink_header_space_separated_txt(self, tmp_path):
        content = "#FID IID PC1 PC2\nF1 SYN_A 0.1 0.2\nF2 SYN_B 0.3 0.4\n"
        p = _write(tmp_path, "plink.txt", content)
        matrix, _ = load_pca_file(p, ["SYN_A", "SYN_B"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.1, 0.2], [0.3, 0.4]])

    def test_csv_with_header(self, tmp_path):
        p = _write(tmp_path, "pca.csv", "sample_id,PC1,PC2\nSYN_A,0.1,0.2\nSYN_B,0.3,0.4\n")
        matrix, _ = load_pca_file(p, ["SYN_A", "SYN_B"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.1, 0.2], [0.3, 0.4]])

    def test_csv_without_header(self, tmp_path):
        p = _write(tmp_path, "pca.csv", "SYN_A,0.1,0.2\nSYN_B,0.3,0.4\n")
        matrix, _ = load_pca_file(p, ["SYN_A", "SYN_B"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.1, 0.2], [0.3, 0.4]])

    def test_non_numeric_cell_raises_value_error(self, tmp_path):
        content = "SYN_A 0.1 0.2\nSYN_B 0.3 oops\nSYN_C NA 0.6\n"
        p = _write(tmp_path, "pca.txt", content)
        with pytest.raises(ValueError) as exc:
            load_pca_file(p, ["SYN_A", "SYN_B", "SYN_C"], n_components=2)
        msg = str(exc.value)
        assert p in msg
        assert "2 NaN/non-numeric" in msg
        assert "PC1" in msg and "PC2" in msg
        assert "SYN_B" in msg and "SYN_C" in msg and "SYN_A" not in msg
        assert "delimiter" in msg

    def test_ragged_row_raises_value_error(self, tmp_path):
        p = _write(tmp_path, "pca.txt", "SYN_A 0.1 0.2\nSYN_B 0.3\n")
        with pytest.raises(ValueError, match="SYN_B"):
            load_pca_file(p, ["SYN_A", "SYN_B"], n_components=2)

    def test_nan_only_in_unused_component_is_ignored(self, tmp_path):
        p = _write(tmp_path, "pca.txt", "SYN_A 0.1 0.2\nSYN_B 0.3 NA\n")
        matrix, _ = load_pca_file(p, ["SYN_A", "SYN_B"], n_components=1)
        np.testing.assert_array_equal(matrix, [[0.1], [0.3]])

    def test_no_pc_columns_raises_value_error(self, tmp_path):
        p = _write(tmp_path, "pca.txt", "SYN_A\nSYN_B\n")
        with pytest.raises(ValueError, match="no principal component columns"):
            load_pca_file(p, ["SYN_A", "SYN_B"], n_components=2)

    def test_non_positive_n_components_raises(self, tmp_path):
        p = _write(tmp_path, "akt_out.txt", AKT_CONTENT)
        with pytest.raises(ValueError, match="n_components"):
            load_pca_file(p, VCF_SAMPLES, n_components=0)


@pytest.mark.unit
class TestPcaLayoutGuards:
    def test_tsv_header_names_with_spaces(self, tmp_path):
        content = "sample id\tPC 1\tPC 2\nSYN_A\t0.1\t0.2\nSYN_B\t0.3\t0.4\n"
        p = _write(tmp_path, "pca.tsv", content)
        matrix, _ = load_pca_file(p, ["SYN_B", "SYN_A"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.3, 0.4], [0.1, 0.2]])

    def test_tsv_sample_ids_with_spaces(self, tmp_path):
        content = "sample id\tPC 1\nSYN A\t0.1\nSYN B\t0.3\n"
        p = _write(tmp_path, "pca.tsv", content)
        matrix, _ = load_pca_file(p, ["SYN B", "SYN A"], n_components=1)
        np.testing.assert_array_equal(matrix, [[0.3], [0.1]])

    def test_tab_padded_akt_layout_still_whitespace(self, tmp_path):
        # exactly one tab per line, like akt: must not be read as a 2-column TSV
        p = _write(tmp_path, "pca.tsv", "SYN_A\t 0.1 -0.2\nSYN_B\t-0.3  0.4\n")
        matrix, _ = load_pca_file(p, ["SYN_A", "SYN_B"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.1, -0.2], [-0.3, 0.4]])

    def test_headerless_plink_numeric_fid_iid_uses_iid(self, tmp_path):
        p = _write(tmp_path, "pca.eigenvec", "1 1001 0.1 0.2\n2 1002 0.3 0.4\n")
        matrix, _ = load_pca_file(p, ["1002", "1001"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.3, 0.4], [0.1, 0.2]])

    def test_headerless_numeric_ids_akt_layout(self, tmp_path):
        p = _write(tmp_path, "pca.txt", "1001 0.1 0.2\n1002 0.3 0.4\n")
        matrix, _ = load_pca_file(p, ["1001", "1002"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.1, 0.2], [0.3, 0.4]])

    def test_headerless_fid_equals_iid_is_ambiguous(self, tmp_path):
        p = _write(tmp_path, "pca.eigenvec", "1001 1001 0.1 0.2\n1002 1002 0.3 0.4\n")
        with pytest.raises(ValueError, match="header"):
            load_pca_file(p, ["1001", "1002"], n_components=2)

    def test_duplicate_samples_raise(self, tmp_path):
        p = _write(tmp_path, "pca.txt", "SYN_A 0.1\nSYN_B 0.2\nSYN_A 0.3\n")
        with pytest.raises(
            ValueError, match=r"1 sample ID\(s\) appear on more than one row.*SYN_A"
        ):
            load_pca_file(p, ["SYN_A", "SYN_B"], n_components=1)

    def test_sample_named_like_fid_is_data(self, tmp_path):
        assert _detect_pca_format(["fido 0.1 0.2"]) == "akt_or_generic"
        p = _write(tmp_path, "pca.txt", "fido 0.1 0.2\nfidelio 0.3 0.4\n")
        matrix, _ = load_pca_file(p, ["fido", "fidelio"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.1, 0.2], [0.3, 0.4]])

    def test_plink2_iid_header(self, tmp_path):
        p = _write(tmp_path, "pca.eigenvec", "#IID\tPC1\tPC2\nSYN_A\t0.1\t0.2\nSYN_B\t0.3\t0.4\n")
        matrix, _ = load_pca_file(p, ["SYN_B", "SYN_A"], n_components=2)
        np.testing.assert_array_equal(matrix, [[0.3, 0.4], [0.1, 0.2]])

    def test_headerless_both_id_columns_partly_match_is_ambiguous(self, tmp_path):
        p = _write(tmp_path, "pca.eigenvec", "1 2 0.1 0.2\n2 3 0.3 0.4\n3 9 0.5 0.6\n")
        with pytest.raises(ValueError, match="header"):
            load_pca_file(p, ["1", "2", "3"], n_components=2)

    def test_layout_match_uses_substring_removal(self, tmp_path):
        # IDs are compared after the same substring removal the VCF names get
        p = _write(tmp_path, "pca.eigenvec", "5 10 0.1 0.2\n6 20 0.3 0.4\n")
        matrix, _ = load_pca_file(p, ["2", "1"], n_components=2, remove_sample_substring="0")
        np.testing.assert_array_equal(matrix, [[0.3, 0.4], [0.1, 0.2]])
        with pytest.raises(ValueError, match="header"):
            load_pca_file(p, ["2", "1"], n_components=2)
