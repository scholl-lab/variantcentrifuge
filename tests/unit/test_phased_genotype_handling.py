"""Tests for phased genotype (pipe-separated) handling across the codebase.

VCF files can contain both unphased (0/1) and phased (0|1) genotypes.
bcftools query preserves the original format, so all genotype-handling
functions must support both separators.

Regression test for: phased genotypes being silently ignored, causing
incorrect allele counts, het/hom classification, and gene burden results.
"""

import pandas as pd
import pytest

from variantcentrifuge.genotype_utils import (
    get_allele_count,
    is_het,
    is_hom_alt,
    is_ref,
    is_variant,
    parse_genotype,
)
from variantcentrifuge.helpers import genotype_to_allele_count
from variantcentrifuge.stats import compute_basic_stats


# ---------------------------------------------------------------------------
# helpers.genotype_to_allele_count
# ---------------------------------------------------------------------------
class TestGenotypeToAlleleCount:
    """Test that genotype_to_allele_count handles both phased and unphased."""

    @pytest.mark.parametrize(
        "gt,expected",
        [
            # Unphased
            ("0/0", 0),
            ("0/1", 1),
            ("1/0", 1),
            ("1/1", 2),
            ("./.", 0),
            # Phased
            ("0|0", 0),
            ("0|1", 1),
            ("1|0", 1),
            ("1|1", 2),
            # Edge cases
            ("", 0),
            (".", 0),
        ],
    )
    def test_allele_count(self, gt, expected):
        assert genotype_to_allele_count(gt) == expected


# ---------------------------------------------------------------------------
# genotype_utils (parse_genotype, is_het, is_hom_alt, etc.)
# ---------------------------------------------------------------------------
class TestGenotypeUtils:
    """Test that genotype_utils handles both phased and unphased."""

    @pytest.mark.parametrize(
        "gt,expected",
        [
            ("0/1", (0, 1)),
            ("1/0", (1, 0)),
            ("0|1", (0, 1)),
            ("1|0", (1, 0)),
            ("1/1", (1, 1)),
            ("1|1", (1, 1)),
            ("0/0", (0, 0)),
            ("0|0", (0, 0)),
            ("./.", (None, None)),
        ],
    )
    def test_parse_genotype(self, gt, expected):
        assert parse_genotype(gt) == expected

    @pytest.mark.parametrize("gt", ["0/1", "1/0", "0|1", "1|0"])
    def test_is_het(self, gt):
        assert is_het(gt) is True

    @pytest.mark.parametrize("gt", ["0/0", "1/1", "0|0", "1|1", "./."])
    def test_is_not_het(self, gt):
        assert is_het(gt) is False

    @pytest.mark.parametrize("gt", ["1/1", "1|1"])
    def test_is_hom_alt(self, gt):
        assert is_hom_alt(gt) is True

    @pytest.mark.parametrize("gt", ["0/0", "0/1", "0|0", "0|1", "./."])
    def test_is_not_hom_alt(self, gt):
        assert is_hom_alt(gt) is False

    @pytest.mark.parametrize("gt", ["0/0", "0|0"])
    def test_is_ref(self, gt):
        assert is_ref(gt) is True

    @pytest.mark.parametrize("gt", ["0/1", "1/1", "0|1", "1|1"])
    def test_is_variant(self, gt):
        assert is_variant(gt) is True

    @pytest.mark.parametrize("gt", ["0/0", "0|0", "./."])
    def test_is_not_variant(self, gt):
        assert is_variant(gt) is False

    @pytest.mark.parametrize(
        "gt,expected",
        [
            ("0/1", 1),
            ("0|1", 1),
            ("1/0", 1),
            ("1|0", 1),
            ("1/1", 2),
            ("1|1", 2),
            ("0/0", 0),
            ("0|0", 0),
            ("./.", 0),
        ],
    )
    def test_get_allele_count(self, gt, expected):
        assert get_allele_count(gt) == expected


# ---------------------------------------------------------------------------
# stats.compute_summary_stats — het/hom counting
# ---------------------------------------------------------------------------
class TestStatsPhasedGenotypes:
    """Test that compute_basic_stats counts phased genotypes correctly."""

    def test_het_hom_counts_with_phased_genotypes(self):
        """Phased genotypes 0|1 and 1|1 must be counted as het and hom."""
        df = pd.DataFrame(
            {
                "GENE": ["A", "A", "A", "A"],
                "GT": [
                    "S1(0/1)",  # unphased het
                    "S2(0|1)",  # phased het
                    "S3(1/1)",  # unphased hom
                    "S4(1|1)",  # phased hom
                ],
            }
        )
        result = compute_basic_stats(df, {"S1", "S2", "S3", "S4"})
        metrics = dict(zip(result["metric"], result["value"], strict=False))
        assert metrics["Het counts"] == "2"
        assert metrics["Hom counts"] == "2"

    def test_only_unphased_still_works(self):
        """Ensure unphased-only data still works correctly."""
        df = pd.DataFrame(
            {
                "GENE": ["A", "A"],
                "GT": [
                    "S1(0/1);S2(1/0)",
                    "S3(1/1)",
                ],
            }
        )
        result = compute_basic_stats(df, {"S1", "S2", "S3"})
        metrics = dict(zip(result["metric"], result["value"], strict=False))
        assert metrics["Het counts"] == "2"
        assert metrics["Hom counts"] == "1"


# ---------------------------------------------------------------------------
# Comp het vectorized — GENOTYPE_ENCODING
# ---------------------------------------------------------------------------
class TestCompHetVectorizedEncoding:
    """Test that GENOTYPE_ENCODING includes phased genotypes."""

    def test_encoding_has_phased_entries(self):
        from variantcentrifuge.inheritance.comp_het_vectorized import GENOTYPE_ENCODING

        assert GENOTYPE_ENCODING["0|1"] == 1
        assert GENOTYPE_ENCODING["1|0"] == 1
        assert GENOTYPE_ENCODING["1|1"] == 2
        assert GENOTYPE_ENCODING["0|0"] == 0

    def test_encoding_matches_unphased(self):
        from variantcentrifuge.inheritance.comp_het_vectorized import GENOTYPE_ENCODING

        assert GENOTYPE_ENCODING["0/1"] == GENOTYPE_ENCODING["0|1"]
        assert GENOTYPE_ENCODING["1/0"] == GENOTYPE_ENCODING["1|0"]
        assert GENOTYPE_ENCODING["1/1"] == GENOTYPE_ENCODING["1|1"]
        assert GENOTYPE_ENCODING["0/0"] == GENOTYPE_ENCODING["0|0"]


# ---------------------------------------------------------------------------
# gene_burden._gt_to_dosage
# ---------------------------------------------------------------------------
class TestGeneBurdenDosage:
    """Test that _gt_to_dosage handles phased genotypes."""

    def test_phased_het(self):
        from variantcentrifuge.gene_burden import _gt_to_dosage

        assert _gt_to_dosage("0|1") == 1
        assert _gt_to_dosage("1|0") == 1

    def test_phased_hom(self):
        from variantcentrifuge.gene_burden import _gt_to_dosage

        assert _gt_to_dosage("1|1") == 2

    def test_phased_ref(self):
        from variantcentrifuge.gene_burden import _gt_to_dosage

        assert _gt_to_dosage("0|0") == 0
