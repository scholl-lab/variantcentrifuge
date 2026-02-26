"""Tests for inheritance analysis performance optimizations.

These tests verify that the performance optimizations produce identical
results to the original code paths, ensuring correctness is preserved.
"""

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.inheritance.analyzer import _create_inheritance_details
from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
    build_variant_keys_array,
    encode_genotypes,
)
from variantcentrifuge.inheritance.vectorized_deducer import vectorized_deduce_patterns


@pytest.fixture
def trio_df():
    """Create a trio DataFrame with compound het variants."""
    return pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "GENE1",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "1",
                "POS": "2000",
                "REF": "C",
                "ALT": "G",
                "GENE": "GENE1",
                "child": "0/1",
                "father": "0/0",
                "mother": "0/1",
            },
            {
                "CHROM": "1",
                "POS": "3000",
                "REF": "G",
                "ALT": "A",
                "GENE": "GENE1",
                "child": "0/0",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "2",
                "POS": "4000",
                "REF": "T",
                "ALT": "C",
                "GENE": "GENE2",
                "child": "0/1",
                "father": "0/0",
                "mother": "0/0",
            },
        ]
    )


@pytest.fixture
def trio_pedigree():
    """Standard trio pedigree data."""
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


@pytest.fixture
def sample_list():
    return ["child", "father", "mother"]


class TestGenotypeMatrixReturn:
    """Test that vectorized_deduce_patterns can return its genotype matrix."""

    def test_returns_patterns_only_by_default(self, trio_df, trio_pedigree, sample_list):
        """Default behavior returns only patterns list."""
        result = vectorized_deduce_patterns(trio_df, trio_pedigree, sample_list)
        assert isinstance(result, list)
        assert len(result) == len(trio_df)

    def test_returns_tuple_when_requested(self, trio_df, trio_pedigree, sample_list):
        """return_genotype_matrix=True returns (patterns, matrix, mapping)."""
        result = vectorized_deduce_patterns(
            trio_df, trio_pedigree, sample_list, return_genotype_matrix=True
        )
        assert isinstance(result, tuple)
        assert len(result) == 3
        patterns, gt_matrix, sample_to_idx = result
        assert isinstance(patterns, list)
        assert isinstance(gt_matrix, np.ndarray)
        assert isinstance(sample_to_idx, dict)

    def test_genotype_matrix_shape(self, trio_df, trio_pedigree, sample_list):
        """Matrix shape is (n_variants x n_samples)."""
        _, gt_matrix, _sample_to_idx = vectorized_deduce_patterns(
            trio_df, trio_pedigree, sample_list, return_genotype_matrix=True
        )
        assert gt_matrix.shape == (len(trio_df), len(sample_list))
        assert gt_matrix.dtype == np.int8

    def test_sample_to_idx_correct(self, trio_df, trio_pedigree, sample_list):
        """sample_to_idx maps all samples correctly."""
        _, gt_matrix, sample_to_idx = vectorized_deduce_patterns(
            trio_df, trio_pedigree, sample_list, return_genotype_matrix=True
        )
        for sample_id in sample_list:
            assert sample_id in sample_to_idx
            idx = sample_to_idx[sample_id]
            # Verify genotype values match encode_genotypes output
            expected = encode_genotypes(trio_df[sample_id])
            np.testing.assert_array_equal(gt_matrix[:, idx], expected)

    def test_patterns_identical_both_modes(self, trio_df, trio_pedigree, sample_list):
        """Patterns are identical whether or not matrix is returned."""
        patterns_only = vectorized_deduce_patterns(
            trio_df, trio_pedigree, sample_list, return_genotype_matrix=False
        )
        patterns_with_matrix, _, _ = vectorized_deduce_patterns(
            trio_df, trio_pedigree, sample_list, return_genotype_matrix=True
        )
        assert patterns_only == patterns_with_matrix

    def test_single_sample_returns_matrix(self, trio_df):
        """Single-sample / no-pedigree path also returns matrix when requested."""
        # Empty pedigree triggers single-sample path
        result = vectorized_deduce_patterns(trio_df, {}, ["child"], return_genotype_matrix=True)
        assert isinstance(result, tuple)
        _patterns, gt_matrix, sample_to_idx = result
        assert gt_matrix.shape[0] == len(trio_df)
        assert "child" in sample_to_idx


class TestCompHetWithPrebuiltMatrix:
    """Test compound het analysis with pre-built genotype matrix."""

    def test_identical_results_with_prebuilt(self, trio_df, trio_pedigree, sample_list):
        """Results are identical with and without pre-built matrix."""
        # Get matrix from Pass 1
        _, gt_matrix, sample_to_idx = vectorized_deduce_patterns(
            trio_df, trio_pedigree, sample_list, return_genotype_matrix=True
        )

        gene_df = trio_df[trio_df["GENE"] == "GENE1"]
        gene_iloc = trio_df.index.get_indexer(gene_df.index)
        gene_row_indices = np.arange(len(trio_df))[gene_iloc]

        # Without pre-built matrix
        result_without = analyze_gene_for_compound_het_vectorized(
            gene_df, trio_pedigree.copy(), sample_list
        )

        # With pre-built matrix
        result_with = analyze_gene_for_compound_het_vectorized(
            gene_df,
            trio_pedigree.copy(),
            sample_list,
            gt_matrix=gt_matrix,
            sample_to_idx=sample_to_idx,
            gene_row_indices=gene_row_indices,
        )

        # Same variant keys
        assert set(result_without.keys()) == set(result_with.keys())

        # Same comp het info per variant
        for vk in result_without:
            assert set(result_without[vk].keys()) == set(result_with[vk].keys())
            for sample_id in result_without[vk]:
                assert (
                    result_without[vk][sample_id]["comp_het_type"]
                    == result_with[vk][sample_id]["comp_het_type"]
                )
                assert (
                    result_without[vk][sample_id]["partner_variants"]
                    == result_with[vk][sample_id]["partner_variants"]
                )

    def test_prebuilt_with_duplicates(self, trio_pedigree, sample_list):
        """Pre-built matrix handles duplicate rows (split snpEff lines)."""
        df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "G1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/0",
                },
                # Duplicate of above (different transcript)
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "G1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/0",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "G1",
                    "child": "0/1",
                    "father": "0/0",
                    "mother": "0/1",
                },
            ]
        )
        _, gt_matrix, sample_to_idx = vectorized_deduce_patterns(
            df, trio_pedigree, sample_list, return_genotype_matrix=True
        )
        gene_row_indices = np.arange(len(df))

        # Should not crash despite duplicates
        result = analyze_gene_for_compound_het_vectorized(
            df,
            trio_pedigree.copy(),
            sample_list,
            gt_matrix=gt_matrix,
            sample_to_idx=sample_to_idx,
            gene_row_indices=gene_row_indices,
        )
        assert isinstance(result, dict)


class TestCreateInheritanceDetailsDict:
    """Test that _create_inheritance_details works with dict input."""

    def test_dict_input_produces_same_output(self, trio_pedigree, sample_list):
        """Dict produces identical output to pd.Series."""
        row_data = {"child": "0/1", "father": "0/0", "mother": "0/1"}

        series_input = pd.Series(row_data)
        dict_input = row_data

        details_series = _create_inheritance_details(
            series_input,
            "autosomal_dominant_possible",
            ["autosomal_dominant_possible"],
            0.5,
            None,
            trio_pedigree,
            sample_list,
        )
        details_dict = _create_inheritance_details(
            dict_input,
            "autosomal_dominant_possible",
            ["autosomal_dominant_possible"],
            0.5,
            None,
            trio_pedigree,
            sample_list,
        )

        assert details_series == details_dict

    def test_dict_with_none_values(self, trio_pedigree, sample_list):
        """Dict with None values (missing genotypes) doesn't crash."""
        row_data = {"child": "0/1", "father": None, "mother": "0/1"}

        details = _create_inheritance_details(
            row_data,
            "unknown",
            ["unknown"],
            0.3,
            None,
            trio_pedigree,
            sample_list,
        )
        # father should be skipped (None)
        sample_ids = [s["sample_id"] for s in details["samples_with_pattern"]]
        assert "father" not in sample_ids
        assert "child" in sample_ids


class TestVariantKeysPrecomputed:
    """Test build_variant_keys_array helper."""

    def test_keys_match_per_row(self, trio_df):
        """Pre-computed keys match per-row create_variant_key_fast output."""
        from variantcentrifuge.inheritance.comp_het_vectorized import create_variant_key_fast

        keys = build_variant_keys_array(trio_df)
        assert len(keys) == len(trio_df)

        for i in range(len(trio_df)):
            expected = create_variant_key_fast(trio_df, i)
            assert keys[i] == expected

    def test_keys_format(self, trio_df):
        """Keys have expected format chr:pos:ref>alt."""
        keys = build_variant_keys_array(trio_df)
        assert keys[0] == "1:1000:A>T"
        assert keys[1] == "1:2000:C>G"
        assert keys[2] == "1:3000:G>A"
        assert keys[3] == "2:4000:T>C"

    def test_comp_het_with_precomputed_keys(self, trio_df, trio_pedigree, sample_list):
        """Comp het results are identical with pre-computed keys."""
        _, gt_matrix, sample_to_idx = vectorized_deduce_patterns(
            trio_df, trio_pedigree, sample_list, return_genotype_matrix=True
        )
        gene_df = trio_df[trio_df["GENE"] == "GENE1"]
        gene_iloc = trio_df.index.get_indexer(gene_df.index)
        gene_row_indices = np.arange(len(trio_df))[gene_iloc]
        gene_vkeys = build_variant_keys_array(trio_df)[gene_iloc]

        # Without pre-computed keys
        result_without = analyze_gene_for_compound_het_vectorized(
            gene_df,
            trio_pedigree.copy(),
            sample_list,
            gt_matrix=gt_matrix,
            sample_to_idx=sample_to_idx,
            gene_row_indices=gene_row_indices,
        )

        # With pre-computed keys
        result_with = analyze_gene_for_compound_het_vectorized(
            gene_df,
            trio_pedigree.copy(),
            sample_list,
            gt_matrix=gt_matrix,
            sample_to_idx=sample_to_idx,
            gene_row_indices=gene_row_indices,
            variant_keys=gene_vkeys,
        )

        assert result_without == result_with
