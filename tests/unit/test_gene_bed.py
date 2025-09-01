"""Unit tests for gene BED creation with merging functionality."""

import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from variantcentrifuge.gene_bed import get_gene_bed


class TestGeneBedMerging:
    """Test BED file creation with overlapping interval merging."""

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmp:
            yield tmp

    def create_mock_genes2bed_output(self, bed_path: str):
        """Create mock snpEff genes2bed output with overlapping intervals."""
        # Simulate LIPE and LIPE-AS1 overlapping intervals (with 1000bp padding)
        overlapping_bed_content = """chr19	42905179	42907179	LIPE
chr19	42906179	42908179	LIPE-AS1
chr20	12345000	12347000	OTHER_GENE
"""
        with open(bed_path, "w") as f:
            f.write(overlapping_bed_content)

    @patch("subprocess.run")
    def test_bed_merging_removes_overlaps(self, mock_subprocess, temp_dir):
        """Test that bedtools merge is called and removes overlapping intervals."""
        # Setup test data
        reference = "GRCh37.p13"
        gene_name = "LIPE LIPE-AS1 OTHER_GENE"
        output_dir = temp_dir

        # Mock subprocess calls
        def mock_run(cmd, **kwargs):
            """Mock subprocess.run calls."""
            if "genes2bed" in cmd:
                # Mock snpEff genes2bed - create overlapping intervals
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    self.create_mock_genes2bed_output(output_file.name)
            elif "sortBed" in cmd:
                # Mock sortBed - copy input to output (already sorted in mock)
                input_file = cmd[cmd.index('-i') + 1]
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    with open(input_file, 'r') as inf, open(output_file.name, 'w') as outf:
                        outf.write(inf.read())
            elif "bedtools" in cmd and "merge" in cmd:
                # Mock bedtools merge - create merged output
                input_file = cmd[cmd.index('-i') + 1]
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    # Simulate merging: chr19 intervals should be merged
                    merged_content = """chr19	42905179	42908179
chr20	12345000	12347000
"""
                    with open(output_file.name, 'w') as f:
                        f.write(merged_content)

        mock_subprocess.side_effect = mock_run

        # Call function
        result_bed = get_gene_bed(
            reference=reference,
            gene_name=gene_name,
            interval_expand=1000,
            add_chr=False,
            output_dir=output_dir
        )

        # Verify bedtools merge was called
        merge_calls = [call for call in mock_subprocess.call_args_list 
                      if len(call[0]) > 0 and "bedtools" in call[0][0] and "merge" in call[0][0]]
        assert len(merge_calls) == 1, "bedtools merge should be called once"
        
        # Verify merge command structure
        merge_cmd = merge_calls[0][0][0]
        assert merge_cmd[0] == "bedtools", "First argument should be bedtools"
        assert merge_cmd[1] == "merge", "Second argument should be merge"
        assert merge_cmd[2] == "-i", "Third argument should be -i"
        assert len(merge_cmd) == 4, "Should have 4 arguments total"

        # Verify result file exists
        assert os.path.exists(result_bed)

    @patch("subprocess.run")
    def test_bed_merging_with_chr_prefix(self, mock_subprocess, temp_dir):
        """Test that chr prefix addition works correctly with merged intervals."""
        # Setup
        reference = "GRCh37.p13"
        gene_name = "LIPE LIPE-AS1"
        
        def mock_run(cmd, **kwargs):
            if "genes2bed" in cmd:
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    # Create intervals without chr prefix
                    with open(output_file.name, 'w') as f:
                        f.write("19\t42905179\t42907179\tLIPE\n19\t42906179\t42908179\tLIPE-AS1\n")
            elif "sortBed" in cmd:
                input_file = cmd[cmd.index('-i') + 1]
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    with open(input_file, 'r') as inf, open(output_file.name, 'w') as outf:
                        outf.write(inf.read())
            elif "bedtools" in cmd and "merge" in cmd:
                input_file = cmd[cmd.index('-i') + 1] 
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    # Merged interval without chr prefix
                    with open(output_file.name, 'w') as f:
                        f.write("19\t42905179\t42908179\n")

        mock_subprocess.side_effect = mock_run

        # Call with add_chr=True
        result_bed = get_gene_bed(
            reference=reference,
            gene_name=gene_name,
            interval_expand=1000,
            add_chr=True,
            output_dir=temp_dir
        )

        # Read result file and verify chr prefix was added to merged intervals
        with open(result_bed, 'r') as f:
            content = f.read()
            assert content.startswith("chr19"), "Chr prefix should be added to merged intervals"

    @patch("subprocess.run")
    def test_bed_merging_preserves_sorting(self, mock_subprocess, temp_dir):
        """Test that merging happens after sorting."""
        reference = "GRCh37.p13"
        gene_name = "TEST_GENE"
        
        # Track call order
        call_order = []
        
        def mock_run(cmd, **kwargs):
            if "genes2bed" in cmd:
                call_order.append("genes2bed")
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    with open(output_file.name, 'w') as f:
                        f.write("chr1\t1000\t2000\tGENE1\n")
            elif "sortBed" in cmd:
                call_order.append("sortBed")
                input_file = cmd[cmd.index('-i') + 1]
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    with open(input_file, 'r') as inf, open(output_file.name, 'w') as outf:
                        outf.write(inf.read())
            elif "bedtools" in cmd and "merge" in cmd:
                call_order.append("bedtools_merge")
                input_file = cmd[cmd.index('-i') + 1]
                output_file = kwargs.get('stdout') 
                if hasattr(output_file, 'name'):
                    with open(input_file, 'r') as inf, open(output_file.name, 'w') as outf:
                        outf.write(inf.read())

        mock_subprocess.side_effect = mock_run

        # Call function
        get_gene_bed(
            reference=reference,
            gene_name=gene_name,
            interval_expand=0,
            add_chr=False,
            output_dir=temp_dir
        )

        # Verify correct order: genes2bed -> sortBed -> bedtools merge
        expected_order = ["genes2bed", "sortBed", "bedtools_merge"]
        assert call_order == expected_order, f"Expected {expected_order}, got {call_order}"

    @patch("subprocess.run")
    def test_bed_merging_caching(self, mock_subprocess, temp_dir):
        """Test that merged BED files are properly cached."""
        reference = "GRCh37.p13"
        gene_name = "TEST_GENE"
        
        def mock_run(cmd, **kwargs):
            if "genes2bed" in cmd:
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    with open(output_file.name, 'w') as f:
                        f.write("chr1\t1000\t2000\tGENE\n")
            elif any(x in cmd for x in ["sortBed", "merge"]):
                input_file = cmd[cmd.index('-i') + 1]
                output_file = kwargs.get('stdout')
                if hasattr(output_file, 'name'):
                    with open(input_file, 'r') as inf, open(output_file.name, 'w') as outf:
                        outf.write(inf.read())

        mock_subprocess.side_effect = mock_run

        # First call - should generate BED file
        result1 = get_gene_bed(
            reference=reference,
            gene_name=gene_name,
            interval_expand=0,
            add_chr=False,
            output_dir=temp_dir
        )

        # Reset mock call count
        mock_subprocess.reset_mock()

        # Second call - should use cached file
        result2 = get_gene_bed(
            reference=reference,
            gene_name=gene_name,
            interval_expand=0,
            add_chr=False,
            output_dir=temp_dir
        )

        # Verify same result and no subprocess calls for cached version
        assert result1 == result2, "Cached result should be identical"
        assert not mock_subprocess.called, "No subprocess calls should be made for cached file"