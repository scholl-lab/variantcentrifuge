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

    def test_encode_genotypes_categorical_without_missing(self):
        """Categorical series without ./. in categories encodes correctly (issue #81)."""
        import numpy as np
        import pandas as pd

        from variantcentrifuge.inheritance.comp_het_vectorized import encode_genotypes

        cat_series = pd.Categorical(["0/0", "0/1", "1/1", "0/1"])
        series = pd.Series(cat_series)
        assert isinstance(series.dtype, pd.CategoricalDtype)

        result = encode_genotypes(series)
        np.testing.assert_array_equal(result, np.array([0, 1, 2, 1], dtype=np.int8))

    def test_encode_genotypes_categorical_with_nan(self):
        """Categorical series with NaN fills ./. and encodes as -1 (issue #81)."""
        import numpy as np
        import pandas as pd

        from variantcentrifuge.inheritance.comp_het_vectorized import encode_genotypes

        cat_series = pd.Categorical(["0/0", None, "1/1", None])
        series = pd.Series(cat_series)
        assert isinstance(series.dtype, pd.CategoricalDtype)

        result = encode_genotypes(series)
        np.testing.assert_array_equal(result, np.array([0, -1, 2, -1], dtype=np.int8))

    def test_encode_genotypes_regular_series(self):
        """Non-categorical series still works (regression guard)."""
        import numpy as np
        import pandas as pd

        from variantcentrifuge.inheritance.comp_het_vectorized import encode_genotypes

        series = pd.Series(["0/0", "0/1", None, "1/1"])
        assert not isinstance(series.dtype, pd.CategoricalDtype)

        result = encode_genotypes(series)
        np.testing.assert_array_equal(result, np.array([0, 1, -1, 2], dtype=np.int8))

    def test_encode_genotypes_categorical_full_pipeline_path(self):
        """Categorical GT columns work through vectorized_deduce_patterns (issue #81)."""
        import pandas as pd

        from variantcentrifuge.inheritance.vectorized_deducer import vectorized_deduce_patterns

        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 200],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["BRCA1", "BRCA1"],
                "SAMPLE1": pd.Categorical(["0/1", "0/0"]),
                "SAMPLE2": pd.Categorical(["0/0", "0/1"]),
            }
        )
        pedigree = {
            "SAMPLE1": {
                "family_id": "FAM1",
                "individual_id": "SAMPLE1",
                "paternal_id": "0",
                "maternal_id": "0",
                "sex": "1",
                "phenotype": "2",
            },
            "SAMPLE2": {
                "family_id": "FAM1",
                "individual_id": "SAMPLE2",
                "paternal_id": "0",
                "maternal_id": "0",
                "sex": "2",
                "phenotype": "1",
            },
        }
        # Should not raise TypeError
        result = vectorized_deduce_patterns(df, pedigree, ["SAMPLE1", "SAMPLE2"])
        assert len(result) == 2
        # Patterns should be real, not empty
        for patterns in result:
            assert isinstance(patterns, list)
            assert len(patterns) > 0

    def test_encode_genotypes_categorical_analyze_inheritance(self):
        """Full analyze_inheritance with Categorical GT produces no 'error' patterns (issue #81)."""
        import pandas as pd

        from variantcentrifuge.inheritance.analyzer import analyze_inheritance

        df = pd.DataFrame(
            {
                "CHROM": ["chr17", "chr17"],
                "POS": [100, 200],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["BRCA1", "BRCA1"],
                "PROBAND": pd.Categorical(["0/1", "1/1"]),
                "FATHER": pd.Categorical(["0/1", "0/1"]),
                "MOTHER": pd.Categorical(["0/0", "0/1"]),
            }
        )
        pedigree = {
            "PROBAND": {
                "family_id": "FAM1",
                "individual_id": "PROBAND",
                "paternal_id": "FATHER",
                "maternal_id": "MOTHER",
                "sex": "1",
                "phenotype": "2",
            },
            "FATHER": {
                "family_id": "FAM1",
                "individual_id": "FATHER",
                "paternal_id": "0",
                "maternal_id": "0",
                "sex": "1",
                "phenotype": "1",
            },
            "MOTHER": {
                "family_id": "FAM1",
                "individual_id": "MOTHER",
                "paternal_id": "0",
                "maternal_id": "0",
                "sex": "2",
                "phenotype": "1",
            },
        }
        result_df = analyze_inheritance(df, pedigree, ["PROBAND", "FATHER", "MOTHER"])
        assert "Inheritance_Pattern" in result_df.columns
        assert (result_df["Inheritance_Pattern"] != "error").all()


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
