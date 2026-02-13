"""Regression tests for compound het false positives with split-snpeff-lines.

When --split-snpeff-lines creates multiple rows per variant (one per transcript),
the compound het analysis must count unique genomic variants, not rows.
See: https://github.com/scholl-lab/variantcentrifuge/issues/33
"""

import pandas as pd
import pytest

from variantcentrifuge.inheritance.comp_het import (
    analyze_gene_for_compound_het as analyze_original,
)
from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
)


@pytest.fixture
def single_variant_three_transcripts():
    """One genomic variant with 3 transcript annotation rows (split-snpeff-lines).

    All rows share the same CHROM/POS/REF/ALT but differ in annotation fields.
    This must NOT produce compound het calls.
    """
    return pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "BRCA1",
                "ANN": "transcript1|missense_variant",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "BRCA1",
                "ANN": "transcript2|synonymous_variant",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "BRCA1",
                "ANN": "transcript3|intron_variant",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
        ]
    )


@pytest.fixture
def two_variants_with_transcript_duplicates():
    """Two unique genomic variants, each with 2 transcript rows.

    Variant 1 (chr1:1000 A>T): from father
    Variant 2 (chr1:2000 C>G): from mother
    This IS a true compound het and must be detected.
    """
    return pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "BRCA1",
                "ANN": "transcript1|missense_variant",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "BRCA1",
                "ANN": "transcript2|synonymous_variant",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "1",
                "POS": "2000",
                "REF": "C",
                "ALT": "G",
                "GENE": "BRCA1",
                "ANN": "transcript1|frameshift_variant",
                "child": "0/1",
                "father": "0/0",
                "mother": "0/1",
            },
            {
                "CHROM": "1",
                "POS": "2000",
                "REF": "C",
                "ALT": "G",
                "GENE": "BRCA1",
                "ANN": "transcript2|splice_variant",
                "child": "0/1",
                "father": "0/0",
                "mother": "0/1",
            },
        ]
    )


@pytest.fixture
def trio_pedigree():
    """Standard trio pedigree with affected child."""
    return {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
    }


SAMPLE_LIST = ["child", "father", "mother"]


class TestCompHetSplitLinesFalsePositive:
    """Test that transcript-duplicated rows do not cause false compound het calls."""

    @pytest.mark.comp_het
    def test_original_no_false_positive(self, single_variant_three_transcripts, trio_pedigree):
        """Original implementation must not call comp het for a single variant."""
        results = analyze_original(
            single_variant_three_transcripts, trio_pedigree.copy(), SAMPLE_LIST
        )
        assert results == {}, (
            f"Expected no compound het calls for a single variant with "
            f"3 transcript rows, got: {results}"
        )

    @pytest.mark.comp_het
    def test_vectorized_no_false_positive(self, single_variant_three_transcripts, trio_pedigree):
        """Vectorized implementation must not call comp het for a single variant."""
        results = analyze_gene_for_compound_het_vectorized(
            single_variant_three_transcripts, trio_pedigree.copy(), SAMPLE_LIST
        )
        assert results == {}, (
            f"Expected no compound het calls for a single variant with "
            f"3 transcript rows, got: {results}"
        )


class TestCompHetSplitLinesTruePositive:
    """Test that true compound hets are still detected with transcript duplicates."""

    @pytest.mark.comp_het
    def test_original_true_positive(self, two_variants_with_transcript_duplicates, trio_pedigree):
        """Original implementation must still detect true comp het with duplicated rows."""
        results = analyze_original(
            two_variants_with_transcript_duplicates, trio_pedigree.copy(), SAMPLE_LIST
        )
        assert len(results) > 0, "Expected compound het detection for 2 unique variants"
        # Both variant keys should be present
        assert "1:1000:A>T" in results
        assert "1:2000:C>G" in results

    @pytest.mark.comp_het
    def test_vectorized_true_positive(self, two_variants_with_transcript_duplicates, trio_pedigree):
        """Vectorized implementation must still detect true comp het with duplicated rows."""
        results = analyze_gene_for_compound_het_vectorized(
            two_variants_with_transcript_duplicates, trio_pedigree.copy(), SAMPLE_LIST
        )
        assert len(results) > 0, "Expected compound het detection for 2 unique variants"
        # Both variant keys should be present
        assert "1:1000:A>T" in results
        assert "1:2000:C>G" in results
