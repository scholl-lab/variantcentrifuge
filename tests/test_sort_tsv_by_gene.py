#!/usr/bin/env python3
"""
Tests for the sort_tsv_by_gene function.

These tests verify the functionality of the memory-efficient
TSV sorting implementation.
"""

import os
import tempfile
import gzip
import pytest
from unittest.mock import patch, Mock
from variantcentrifuge.pipeline import sort_tsv_by_gene


class TestSortTsvByGene:
    """Test the sort_tsv_by_gene function with various scenarios."""

    def test_sort_basic_tsv(self):
        """Test sorting a basic TSV file."""
        # Create test data
        test_data = """CHROM\tPOS\tGENE\tREF\tALT
chr1\t100\tGENE3\tA\tT
chr1\t200\tGENE1\tC\tG
chr1\t300\tGENE2\tG\tA
chr1\t400\tGENE1\tT\tC
chr1\t500\tGENE3\tA\tG"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(test_data)
            input_file = f.name

        output_file = input_file + ".sorted"

        try:
            # Sort the file
            _ = sort_tsv_by_gene(input_file, output_file, gene_column="GENE")

            # Verify output file exists
            assert os.path.exists(output_file)
            assert _ == output_file

            # Read and verify content
            with open(output_file, "r") as f:
                lines = f.readlines()

            # Check header is preserved
            assert lines[0].strip() == "CHROM\tPOS\tGENE\tREF\tALT"

            # Check data is sorted by gene
            genes = [line.split("\t")[2] for line in lines[1:]]
            assert genes == ["GENE1", "GENE1", "GENE2", "GENE3", "GENE3"]

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_sort_gzipped_files(self):
        """Test sorting with gzipped input and output."""
        test_data = """CHROM\tPOS\tGENE\tREF\tALT
chr1\t100\tZZZ\tA\tT
chr1\t200\tAAA\tC\tG
chr1\t300\tMMM\tG\tA"""

        with tempfile.NamedTemporaryFile(mode="wb", suffix=".tsv.gz", delete=False) as f:
            f.write(gzip.compress(test_data.encode("utf-8")))
            input_file = f.name

        output_file = input_file.replace(".tsv.gz", ".sorted.tsv.gz")

        try:
            # Sort the file
            _ = sort_tsv_by_gene(input_file, output_file, gene_column="GENE")

            # Verify output file exists
            assert os.path.exists(output_file)
            assert output_file.endswith(".gz")

            # Read and verify content
            with gzip.open(output_file, "rt") as f:
                lines = f.readlines()

            # Check data is sorted
            genes = [line.split("\t")[2] for line in lines[1:]]
            assert genes == ["AAA", "MMM", "ZZZ"]

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_sort_with_custom_parameters(self):
        """Test sorting with custom memory limit and parallel threads."""
        test_data = """CHROM\tPOS\tGENE
chr1\t100\tGENE1
chr1\t200\tGENE2"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(test_data)
            input_file = f.name

        output_file = input_file + ".sorted"
        temp_dir = tempfile.gettempdir()

        try:
            # Sort with custom parameters
            _ = sort_tsv_by_gene(
                input_file,
                output_file,
                gene_column="GENE",
                temp_dir=temp_dir,
                memory_limit="500M",
                parallel=2,
            )

            # Verify output file exists
            assert os.path.exists(output_file)

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_empty_file_error(self):
        """Test error handling for empty files."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            # Write empty file
            f.write("")
            input_file = f.name

        output_file = input_file + ".sorted"

        try:
            # Should raise ValueError for empty file
            with pytest.raises(ValueError, match="empty or has no header"):
                sort_tsv_by_gene(input_file, output_file)

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_missing_gene_column_error(self):
        """Test error handling when gene column is missing."""
        test_data = """CHROM\tPOS\tREF\tALT
chr1\t100\tA\tT
chr1\t200\tC\tG"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(test_data)
            input_file = f.name

        output_file = input_file + ".sorted"

        try:
            # Should raise ValueError for missing GENE column
            with pytest.raises(ValueError, match="Gene column 'GENE' not found"):
                sort_tsv_by_gene(input_file, output_file, gene_column="GENE")

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_sort_command_failure(self):
        """Test error handling when sort command fails."""
        test_data = """CHROM\tPOS\tGENE
chr1\t100\tGENE1"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(test_data)
            input_file = f.name

        output_file = "/invalid/path/output.tsv"  # Invalid output path

        try:
            # Should raise RuntimeError when sort fails
            with pytest.raises(RuntimeError, match="Failed to sort TSV file"):
                sort_tsv_by_gene(input_file, output_file)

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)

    def test_special_characters_in_filename(self):
        """Test handling of filenames with special characters."""
        test_data = """CHROM\tPOS\tGENE
chr1\t100\tGENE1
chr1\t200\tGENE2"""

        # Create file with spaces and special characters
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=" test file (v1).tsv", delete=False, prefix="data "
        ) as f:
            f.write(test_data)
            input_file = f.name

        output_file = input_file + ".sorted"

        try:
            # Should handle special characters properly
            _ = sort_tsv_by_gene(input_file, output_file)
            assert os.path.exists(output_file)

            # Verify content
            with open(output_file, "r") as f:
                lines = f.readlines()
            assert len(lines) == 3  # Header + 2 data lines

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_alternate_gene_column_name(self):
        """Test sorting with non-standard gene column name."""
        test_data = """CHROM\tPOS\tANN[0].GENE\tREF\tALT
chr1\t100\tGENE3\tA\tT
chr1\t200\tGENE1\tC\tG
chr1\t300\tGENE2\tG\tA"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(test_data)
            input_file = f.name

        output_file = input_file + ".sorted"

        try:
            # Sort using alternate column name
            _ = sort_tsv_by_gene(input_file, output_file, gene_column="ANN[0].GENE")

            # Verify sorting worked
            with open(output_file, "r") as f:
                lines = f.readlines()

            genes = [line.split("\t")[2] for line in lines[1:]]
            assert genes == ["GENE1", "GENE2", "GENE3"]

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    @patch("subprocess.run")
    def test_bash_not_found_fallback(self, mock_run):
        """Test fallback when bash is not found."""
        # Mock subprocess to simulate successful sort
        mock_run.return_value = Mock(returncode=0, stderr="")

        test_data = """CHROM\tPOS\tGENE
chr1\t100\tGENE1"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(test_data)
            input_file = f.name

        output_file = input_file + ".sorted"

        # Create output file to simulate successful sort
        with open(output_file, "w") as f:
            f.write(test_data)

        try:
            with patch("shutil.which", return_value=None):
                with patch(
                    "os.path.exists", side_effect=lambda x: x == output_file or x == input_file
                ):
                    sort_tsv_by_gene(input_file, output_file)

            # Should complete successfully with None executable
            assert mock_run.called
            assert mock_run.call_args[1]["executable"] is None

        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
