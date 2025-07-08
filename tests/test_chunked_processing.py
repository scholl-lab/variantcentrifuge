"""Test the chunked processing functionality for large files."""

import os
import tempfile
import pandas as pd
import pytest
from unittest.mock import Mock, patch

from variantcentrifuge.pipeline import (
    read_tsv_in_gene_chunks,
    sort_tsv_by_gene,
    process_chunked_pipeline,
)


class TestGeneAwareChunking:
    """Test the gene-aware chunking reader."""

    def test_read_tsv_in_gene_chunks_basic(self):
        """Test basic chunking with complete genes."""
        # Create test data
        data = pd.DataFrame(
            {
                "CHROM": ["chr1"] * 10 + ["chr2"] * 10,
                "POS": list(range(100, 110)) + list(range(200, 210)),
                "GENE": ["GENE1"] * 5 + ["GENE2"] * 5 + ["GENE3"] * 5 + ["GENE4"] * 5,
                "REF": ["A"] * 20,
                "ALT": ["T"] * 20,
            }
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            data.to_csv(f, sep="\t", index=False)
            temp_file = f.name

        try:
            # Read in chunks of size 6
            chunks = list(read_tsv_in_gene_chunks(temp_file, gene_column="GENE", chunksize=6))

            # The chunking algorithm keeps genes together, but the exact chunking
            # depends on the algorithm's logic. With chunksize=6:
            # - It will try to keep at least 6 rows per chunk
            # - It won't split genes across chunks
            # - GENE1 (5 rows) + GENE2 (5 rows) = 10 rows in first chunk
            # - GENE3 (5 rows) + GENE4 (5 rows) = 10 rows in second chunk
            assert len(chunks) == 2

            # First chunk should contain GENE1 and GENE2
            assert set(chunks[0]["GENE"].unique()) == {"GENE1", "GENE2"}
            assert len(chunks[0]) == 10

            # Second chunk should contain GENE3 and GENE4
            assert set(chunks[1]["GENE"].unique()) == {"GENE3", "GENE4"}
            assert len(chunks[1]) == 10

        finally:
            os.unlink(temp_file)

    def test_read_tsv_in_gene_chunks_large_gene(self):
        """Test chunking when a single gene exceeds chunk size."""
        # Create test data with one very large gene
        data = pd.DataFrame(
            {
                "CHROM": ["chr1"] * 20,
                "POS": list(range(100, 120)),
                "GENE": ["LARGE_GENE"] * 15 + ["SMALL_GENE"] * 5,
                "REF": ["A"] * 20,
                "ALT": ["T"] * 20,
            }
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            data.to_csv(f, sep="\t", index=False)
            temp_file = f.name

        try:
            # Read in chunks of size 5 (smaller than LARGE_GENE)
            chunks = list(read_tsv_in_gene_chunks(temp_file, gene_column="GENE", chunksize=5))

            # Should have 2 chunks despite small chunksize
            assert len(chunks) == 2

            # First chunk should contain only LARGE_GENE (15 rows)
            assert set(chunks[0]["GENE"].unique()) == {"LARGE_GENE"}
            assert len(chunks[0]) == 15

            # Second chunk should contain only SMALL_GENE (5 rows)
            assert set(chunks[1]["GENE"].unique()) == {"SMALL_GENE"}
            assert len(chunks[1]) == 5

        finally:
            os.unlink(temp_file)

    def test_read_tsv_in_gene_chunks_gzipped(self):
        """Test chunking with gzipped input file."""
        # Create test data
        data = pd.DataFrame(
            {
                "CHROM": ["chr1"] * 10,
                "POS": list(range(100, 110)),
                "GENE": ["GENE1"] * 5 + ["GENE2"] * 5,
                "REF": ["A"] * 10,
                "ALT": ["T"] * 10,
            }
        )

        with tempfile.NamedTemporaryFile(mode="wb", suffix=".tsv.gz", delete=False) as f:
            data.to_csv(f, sep="\t", index=False, compression="gzip")
            temp_file = f.name

        try:
            # Read in chunks
            chunks = list(read_tsv_in_gene_chunks(temp_file, gene_column="GENE", chunksize=6))

            # With chunksize=6 and GENE1 (5 rows) + GENE2 (5 rows) = 10 rows total,
            # the algorithm will keep both genes together in one chunk since:
            # - Gene boundary at index 5 is less than chunksize 6
            # - No other boundaries exist after the chunksize threshold
            assert len(chunks) == 1

            # Verify both genes are in the single chunk
            assert set(chunks[0]["GENE"].unique()) == {"GENE1", "GENE2"}
            assert len(chunks[0]) == 10

        finally:
            os.unlink(temp_file)


class TestSortByGene:
    """Test the memory-efficient gene sorting."""

    @patch("subprocess.run")
    def test_sort_tsv_by_gene_basic(self, mock_run):
        """Test basic sorting functionality."""
        # Create test data
        data = pd.DataFrame(
            {
                "CHROM": ["chr1"] * 5,
                "POS": [100, 101, 102, 103, 104],
                "GENE": ["GENE3", "GENE1", "GENE2", "GENE1", "GENE2"],
                "REF": ["A"] * 5,
                "ALT": ["T"] * 5,
            }
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            data.to_csv(f, sep="\t", index=False)
            input_file = f.name

        output_file = input_file + ".sorted"

        # Mock successful execution and create the output file
        def mock_run_side_effect(cmd, **kwargs):
            # Create the output file with sorted data
            sorted_data = data.sort_values("GENE")
            sorted_data.to_csv(output_file, sep="\t", index=False)
            return Mock(returncode=0, stderr="")

        mock_run.side_effect = mock_run_side_effect

        try:
            # Call sort function
            sort_tsv_by_gene(
                input_file, output_file, gene_column="GENE", memory_limit="100M", parallel=2
            )

            # Verify subprocess was called
            assert mock_run.called

            # Check the command includes our parameters
            cmd = mock_run.call_args[0][0]
            assert "--buffer-size=100M" in cmd
            assert "--parallel=2" in cmd
            assert "-k3,3" in cmd  # GENE is column 3

            # Verify output file was created
            assert os.path.exists(output_file)

        finally:
            os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    @patch("subprocess.run")
    def test_sort_tsv_by_gene_gzipped(self, mock_run):
        """Test sorting with gzipped input/output."""
        # Create test data
        data = pd.DataFrame(
            {
                "CHROM": ["chr1"] * 3,
                "POS": [100, 101, 102],
                "GENE": ["GENE2", "GENE1", "GENE3"],
            }
        )

        with tempfile.NamedTemporaryFile(mode="wb", suffix=".tsv.gz", delete=False) as f:
            data.to_csv(f, sep="\t", index=False, compression="gzip")
            input_file = f.name

        output_file = input_file.replace(".tsv.gz", ".sorted.tsv.gz")

        # Mock successful execution and create the output file
        def mock_run_side_effect(cmd, **kwargs):
            # Create the output file with sorted data
            sorted_data = data.sort_values("GENE")
            sorted_data.to_csv(output_file, sep="\t", index=False, compression="gzip")
            return Mock(returncode=0, stderr="")

        mock_run.side_effect = mock_run_side_effect

        try:
            # Call sort function
            sort_tsv_by_gene(input_file, output_file)

            # Verify subprocess was called
            assert mock_run.called

            # Check the command handles gzip properly
            cmd = mock_run.call_args[0][0]
            assert "gzip -cd" in cmd  # Decompress input
            assert "gzip -c >" in cmd  # Compress output

            # Verify output file was created
            assert os.path.exists(output_file)

        finally:
            os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)


class TestChunkedPipeline:
    """Test the full chunked pipeline processing."""

    @patch("variantcentrifuge.pipeline_core.analyze_variants")
    @patch("variantcentrifuge.pipeline_core.annotate_dataframe_with_features")
    def test_process_chunked_pipeline_basic(self, mock_annotate, mock_analyze):
        """Test basic chunked pipeline processing."""
        # Create test data
        data = pd.DataFrame(
            {
                "CHROM": ["chr1"] * 10,
                "POS": list(range(100, 110)),
                "GENE": ["GENE1"] * 5 + ["GENE2"] * 5,
                "REF": ["A"] * 10,
                "ALT": ["T"] * 10,
                "IMPACT": ["HIGH"] * 10,
            }
        )

        # Mock analyze_variants to return processed lines
        def mock_analyze_func(file_handle, cfg):
            lines = file_handle.readlines()
            if lines:
                yield lines[0].strip()  # Header
                for line in lines[1:]:
                    yield line.strip()

        mock_analyze.side_effect = mock_analyze_func
        mock_annotate.side_effect = lambda df, features: df

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            data.to_csv(f, sep="\t", index=False)
            input_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_file = f.name

        try:
            # Create minimal config and args
            cfg = {
                "keep_intermediates": False,
                "calculate_inheritance": False,
                "late_filtering": False,
                "no_links": True,
            }

            args = Mock()
            args.add_column = None

            custom_features = {}

            # Process pipeline
            success = process_chunked_pipeline(
                final_tsv=input_file,
                final_output=output_file,
                cfg=cfg,
                custom_features=custom_features,
                scoring_config=None,
                pedigree_data=None,
                args=args,
                base_name="test",
                intermediate_dir=tempfile.gettempdir(),
                chunksize=6,
            )

            assert success is True

            # Verify output was created
            assert os.path.exists(output_file)

            # Read output and verify VAR_ID was added
            with open(output_file, "r") as f:
                lines = f.readlines()

            # Should have header + 10 data rows
            assert len(lines) == 11

            # Check header has VAR_ID
            assert lines[0].startswith("VAR_ID\t")

            # Check data rows have VAR_ID
            for i, line in enumerate(lines[1:], 1):
                assert line.startswith(f"var_{i:04d}_")

        finally:
            os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_process_chunked_pipeline_no_gene_column(self):
        """Test chunked processing fails gracefully without GENE column."""
        # Create test data without GENE column
        data = pd.DataFrame(
            {
                "CHROM": ["chr1"] * 5,
                "POS": list(range(100, 105)),
                "REF": ["A"] * 5,
                "ALT": ["T"] * 5,
            }
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            data.to_csv(f, sep="\t", index=False)
            input_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_file = f.name

        try:
            cfg = {"keep_intermediates": False}
            args = Mock()

            # Process pipeline
            success = process_chunked_pipeline(
                final_tsv=input_file,
                final_output=output_file,
                cfg=cfg,
                custom_features={},
                scoring_config=None,
                pedigree_data=None,
                args=args,
                base_name="test",
                intermediate_dir=tempfile.gettempdir(),
            )

            # Should return False (fall back to regular processing)
            assert success is False

        finally:
            os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
