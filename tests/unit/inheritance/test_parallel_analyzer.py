"""Test parallel inheritance analyzer."""

import pytest
import pandas as pd
from unittest.mock import patch, Mock
import json

from variantcentrifuge.inheritance.parallel_analyzer import (
    analyze_inheritance_parallel,
    _process_gene_group,
)


class TestParallelInheritanceAnalyzer:
    """Test parallel inheritance analysis functionality."""

    @pytest.fixture
    def sample_df(self):
        """Create a sample DataFrame with multiple genes."""
        data = {
            "CHROM": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3"],
            "POS": [100, 200, 300, 400, 500, 600],
            "REF": ["A", "T", "G", "C", "A", "T"],
            "ALT": ["T", "C", "A", "T", "G", "A"],
            "GENE": ["GENE1", "GENE1", "GENE2", "GENE2", "GENE3", "GENE3"],
            "proband": ["0/1", "0/1", "0/1", "1/1", "0/1", "0/0"],
            "mother": ["0/1", "0/0", "0/1", "0/1", "0/0", "0/1"],
            "father": ["0/0", "0/1", "0/0", "0/1", "0/1", "0/0"],
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def pedigree_data(self):
        """Create sample pedigree data."""
        return {
            "proband": {
                "sample_id": "proband",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "1",
                "affected_status": "2",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "father": {
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            },
        }

    def test_process_gene_group_basic(self, sample_df, pedigree_data):
        """Test processing a single gene group."""
        gene_df = sample_df[sample_df["GENE"] == "GENE1"]
        sample_list = ["proband", "mother", "father"]

        gene, results = _process_gene_group(
            "GENE1", gene_df, pedigree_data, sample_list, use_vectorized=False
        )

        assert gene == "GENE1"
        assert isinstance(results, dict)

    def test_process_gene_group_empty(self, pedigree_data):
        """Test processing empty gene group."""
        empty_df = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "GENE"])
        sample_list = ["proband", "mother", "father"]

        gene, results = _process_gene_group("GENE1", empty_df, pedigree_data, sample_list)

        assert gene == "GENE1"
        assert results == {}

    def test_process_gene_group_single_variant(self, sample_df, pedigree_data):
        """Test processing gene with single variant (no compound het possible)."""
        single_df = sample_df.iloc[0:1]
        sample_list = ["proband", "mother", "father"]

        gene, results = _process_gene_group("GENE1", single_df, pedigree_data, sample_list)

        assert gene == "GENE1"
        assert results == {}

    def test_analyze_inheritance_parallel_basic(self, sample_df, pedigree_data):
        """Test basic parallel inheritance analysis."""
        sample_list = ["proband", "mother", "father"]

        result = analyze_inheritance_parallel(
            sample_df,
            pedigree_data,
            sample_list,
            use_vectorized_comp_het=False,
            n_workers=2,
            min_variants_for_parallel=1,
        )

        # Check required columns exist
        assert "Inheritance_Pattern" in result.columns
        assert "Inheritance_Details" in result.columns

        # Check all rows have patterns
        assert result["Inheritance_Pattern"].notna().all()

        # Check details are valid JSON
        for details in result["Inheritance_Details"]:
            parsed = json.loads(details)
            assert "primary_pattern" in parsed
            assert "all_patterns" in parsed
            assert "confidence" in parsed

    def test_analyze_inheritance_parallel_empty_df(self, pedigree_data):
        """Test parallel analysis with empty DataFrame."""
        empty_df = pd.DataFrame()
        sample_list = ["proband", "mother", "father"]

        result = analyze_inheritance_parallel(empty_df, pedigree_data, sample_list)

        assert len(result) == 0
        assert "Inheritance_Pattern" in result.columns
        assert "Inheritance_Details" in result.columns

    def test_analyze_inheritance_parallel_no_pedigree(self, sample_df):
        """Test parallel analysis without pedigree data."""
        sample_list = ["proband", "mother", "father"]

        result = analyze_inheritance_parallel(sample_df, {}, sample_list)

        # Should still work but treat as unrelated individuals
        assert "Inheritance_Pattern" in result.columns
        assert len(result) == len(sample_df)

    @patch("variantcentrifuge.inheritance.parallel_analyzer.ProcessPoolExecutor")
    def test_parallel_execution(self, mock_executor_class, sample_df, pedigree_data):
        """Test that parallel execution is used when appropriate."""
        sample_list = ["proband", "mother", "father"]

        # Mock the executor
        mock_executor = Mock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor
        mock_future = Mock()
        mock_future.result.return_value = ("GENE1", {})
        mock_executor.submit.return_value = mock_future

        # Mock as_completed to return our futures
        with patch(
            "variantcentrifuge.inheritance.parallel_analyzer.as_completed"
        ) as mock_completed:
            mock_completed.return_value = [mock_future]

            analyze_inheritance_parallel(
                sample_df,
                pedigree_data,
                sample_list,
                n_workers=2,
                min_variants_for_parallel=1,
            )

        # Check executor was created with correct workers
        mock_executor_class.assert_called_once_with(max_workers=2)

        # Check genes were submitted for processing
        assert mock_executor.submit.called

    def test_sequential_fallback_small_dataset(self, sample_df, pedigree_data):
        """Test that sequential processing is used for small datasets."""
        sample_list = ["proband", "mother", "father"]

        # Set high threshold so sequential processing is used
        result = analyze_inheritance_parallel(
            sample_df,
            pedigree_data,
            sample_list,
            n_workers=2,
            min_variants_for_parallel=1000,  # High threshold
        )

        # Should still produce valid results
        assert "Inheritance_Pattern" in result.columns
        assert len(result) == len(sample_df)

    def test_timing_information_logged(self, sample_df, pedigree_data, caplog):
        """Test that timing information is logged."""
        import logging

        sample_list = ["proband", "mother", "father"]

        with caplog.at_level(logging.INFO):
            analyze_inheritance_parallel(
                sample_df,
                pedigree_data,
                sample_list,
                min_variants_for_parallel=1,
            )

        # Check for timing logs
        log_messages = [record.message for record in caplog.records]
        timing_logs = [msg for msg in log_messages if "complete in" in msg]
        assert len(timing_logs) > 0

        # Check for pattern distribution log
        pattern_logs = [msg for msg in log_messages if "Pattern distribution" in msg]
        assert len(pattern_logs) > 0

    def test_compound_het_detection_parallel(self, pedigree_data):
        """Test compound heterozygous detection works correctly in parallel."""
        # Create data with potential compound het
        data = {
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "REF": ["A", "T"],
            "ALT": ["T", "C"],
            "GENE": ["GENE1", "GENE1"],
            "proband": ["0/1", "0/1"],  # Both heterozygous in proband
            "mother": ["0/1", "0/0"],  # First from mother
            "father": ["0/0", "0/1"],  # Second from father
        }
        df = pd.DataFrame(data)
        sample_list = ["proband", "mother", "father"]

        result = analyze_inheritance_parallel(
            df,
            pedigree_data,
            sample_list,
            use_vectorized_comp_het=False,
            n_workers=2,
            min_variants_for_parallel=1,
        )

        # Check that compound het was detected
        patterns = result["Inheritance_Pattern"].tolist()
        assert any("compound" in pattern for pattern in patterns)

        # Check details contain compound het info
        for idx, row in result.iterrows():
            details = json.loads(row["Inheritance_Details"])
            if "compound" in row["Inheritance_Pattern"]:
                samples = details.get("samples_with_pattern", [])
                proband_info = next((s for s in samples if s["sample_id"] == "proband"), None)
                if proband_info:
                    assert (
                        "compound_het_partner" in proband_info
                        or "compound_het_partners" in proband_info
                    )

    def test_error_handling_in_worker(self, sample_df, pedigree_data):
        """Test error handling when worker process fails."""
        sample_list = ["proband", "mother", "father"]

        with patch(
            "variantcentrifuge.inheritance.parallel_analyzer" ".analyze_gene_for_compound_het"
        ) as mock_analyze:
            # Make the function raise an error
            mock_analyze.side_effect = Exception("Test error")

            # Should still complete but log error
            result = analyze_inheritance_parallel(
                sample_df,
                pedigree_data,
                sample_list,
                use_vectorized_comp_het=False,
                n_workers=2,
                min_variants_for_parallel=1,
            )

            # Should still return valid DataFrame
            assert isinstance(result, pd.DataFrame)
            assert "Inheritance_Pattern" in result.columns

    def test_vectorized_vs_original_consistency(self, sample_df, pedigree_data):
        """Test that vectorized and original implementations produce consistent results."""
        sample_list = ["proband", "mother", "father"]

        # Run with original implementation
        result_original = analyze_inheritance_parallel(
            sample_df,
            pedigree_data,
            sample_list,
            use_vectorized_comp_het=False,
            n_workers=1,  # Single worker to ensure deterministic
        )

        # Run with vectorized implementation (if available)
        try:
            result_vectorized = analyze_inheritance_parallel(
                sample_df,
                pedigree_data,
                sample_list,
                use_vectorized_comp_het=True,
                n_workers=1,  # Single worker to ensure deterministic
            )

            # Compare patterns (exact match)
            assert result_original["Inheritance_Pattern"].equals(
                result_vectorized["Inheritance_Pattern"]
            )
        except ImportError:
            # Vectorized implementation not available, skip comparison
            pytest.skip("Vectorized implementation not available")
