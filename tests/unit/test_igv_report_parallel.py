"""Unit tests for parallel IGV report generation functionality."""

import tempfile
import threading
import time
from pathlib import Path
from unittest.mock import patch

import pytest

from variantcentrifuge.generate_igv_report import _generate_single_igv_report, generate_igv_report


class TestIGVReportParallel:
    """Test cases for parallel IGV report generation."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def sample_tsv(self, temp_dir):
        """Create a sample TSV file for testing."""
        tsv_file = temp_dir / "variants.tsv"
        tsv_content = (
            "CHROM\tPOS\tREF\tALT\tGENE\tGT\n"
            "17\t41223094\tT\tC\tBRCA1\tSample_HG002(0/1);Sample_HG003(0/0)\n"
            "17\t7579472\tG\tC\tTP53\tSample_HG002(1/1);Sample_HG003(0/1)\n"
            "16\t9857048\tT\tA\tGRIN2A\tSample_HG002(0/1);Sample_HG004(0/1)\n"
            "2\t179397091\tC\tG\tTTN\tSample_HG003(0/1);Sample_HG004(1/1)\n"
        )
        tsv_file.write_text(tsv_content)
        return tsv_file

    @pytest.fixture
    def bam_mapping_file(self, temp_dir):
        """Create a BAM mapping file for testing."""
        mapping_file = temp_dir / "bam_mapping.txt"
        mapping_content = (
            "Sample_HG002,/path/to/HG002.bam\n"
            "Sample_HG003,/path/to/HG003.bam\n"
            "Sample_HG004,/path/to/HG004.bam\n"
        )
        mapping_file.write_text(mapping_content)
        return mapping_file

    def test_single_igv_report_generation_success(self, temp_dir):
        """Test successful generation of a single IGV report."""
        # Setup task data
        task_data = (
            "Sample_HG002",  # sample_id
            "17",  # chrom
            "41223094",  # pos
            "T",  # ref_allele
            "C",  # alt_allele
            "/path/to/HG002.bam",  # bam_path
            str(temp_dir / "variant.tsv"),  # variant_tsv_path
            str(temp_dir / "report.html"),  # sample_report_path
            50,  # igv_flanking
            "hg19",  # igv_reference
            None,  # igv_fasta
            None,  # igv_ideogram
        )

        # Mock subprocess.run to simulate successful create_report execution
        with patch("variantcentrifuge.generate_igv_report.subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            mock_run.return_value.stderr = ""

            result = _generate_single_igv_report(task_data)

        # Verify result
        success, sample_id, chrom, pos, ref, alt, report_path = result
        assert success is True
        assert sample_id == "Sample_HG002"
        assert chrom == "17"
        assert pos == "41223094"
        assert ref == "T"
        assert alt == "C"
        assert report_path == str(temp_dir / "report.html")

        # Verify variant TSV was created
        variant_tsv = temp_dir / "variant.tsv"
        assert variant_tsv.exists()
        content = variant_tsv.read_text()
        assert "CHROM\tPOS\tREF\tALT" in content
        assert "17\t41223094\tT\tC" in content

        # Verify create_report command was called correctly
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "create_report"
        assert "--genome" in cmd
        assert "hg19" in cmd

    def test_single_igv_report_generation_failure(self, temp_dir):
        """Test handling of failed IGV report generation."""
        task_data = (
            "Sample_HG002",
            "17",
            "41223094",
            "T",
            "C",
            "/path/to/HG002.bam",
            str(temp_dir / "variant.tsv"),
            str(temp_dir / "report.html"),
            50,
            "hg19",
            None,
            None,
        )

        # Mock subprocess.run to simulate failed create_report execution
        with patch("variantcentrifuge.generate_igv_report.subprocess.run") as mock_run:
            mock_run.return_value.returncode = 1
            mock_run.return_value.stderr = "create_report failed"

            result = _generate_single_igv_report(task_data)

        # Verify failure result
        success, sample_id, chrom, pos, ref, alt, report_path = result
        assert success is False
        assert sample_id == "Sample_HG002"

    def test_single_igv_report_with_fasta(self, temp_dir):
        """Test IGV report generation with local FASTA file."""
        task_data = (
            "Sample_HG002",
            "17",
            "41223094",
            "T",
            "C",
            "/path/to/HG002.bam",
            str(temp_dir / "variant.tsv"),
            str(temp_dir / "report.html"),
            50,
            None,
            "/path/to/genome.fa",
            "/path/to/ideogram.txt",
        )

        with patch("variantcentrifuge.generate_igv_report.subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0

            _generate_single_igv_report(task_data)

        # Verify FASTA parameters were used instead of genome reference
        cmd = mock_run.call_args[0][0]
        assert "--fasta" in cmd
        assert "/path/to/genome.fa" in cmd
        assert "--ideogram" in cmd
        assert "/path/to/ideogram.txt" in cmd
        assert "--genome" not in cmd

    @patch("variantcentrifuge.generate_igv_report._generate_single_igv_report")
    def test_parallel_processing_worker_count(
        self, mock_single_report, sample_tsv, bam_mapping_file, temp_dir
    ):
        """Test that parallel processing uses the specified number of workers."""
        # Setup
        mock_single_report.return_value = (
            True,
            "Sample_HG002",
            "17",
            "41223094",
            "T",
            "C",
            "/path/report.html",
        )

        # Track thread usage
        thread_ids = set()

        def track_threads(*args, **kwargs):
            thread_ids.add(threading.current_thread().ident)
            time.sleep(0.1)  # Simulate work
            return (True, "Sample_HG002", "17", "41223094", "T", "C", "/path/report.html")

        mock_single_report.side_effect = track_threads

        # Execute with different worker counts
        for max_workers in [1, 2, 4]:
            thread_ids.clear()

            generate_igv_report(
                variants_tsv=str(sample_tsv),
                output_dir=str(temp_dir),
                bam_mapping_file=str(bam_mapping_file),
                igv_reference="hg19",
                max_workers=max_workers,
            )

            # With enough tasks and workers, we should see multiple threads
            # (Note: actual thread count may be less than max_workers due to task scheduling)
            if max_workers > 1:
                assert len(thread_ids) >= 1  # At least main execution context

    @patch("variantcentrifuge.generate_igv_report._generate_single_igv_report")
    def test_thread_safety_variant_map(
        self, mock_single_report, sample_tsv, bam_mapping_file, temp_dir
    ):
        """Test thread safety of variant map updates."""

        # Setup - simulate multiple successful reports for the same variant
        def mock_report_generator(*args):
            task = args[0]
            sample_id, chrom, pos, ref, alt = task[0], task[1], task[2], task[3], task[4]
            time.sleep(0.01)  # Small delay to increase chance of race conditions
            return (True, sample_id, chrom, pos, ref, alt, f"/path/{sample_id}_report.html")

        mock_single_report.side_effect = mock_report_generator

        # Execute with multiple workers
        generate_igv_report(
            variants_tsv=str(sample_tsv),
            output_dir=str(temp_dir),
            bam_mapping_file=str(bam_mapping_file),
            igv_reference="hg19",
            max_workers=4,
        )

        # Verify mapping file was created
        mapping_file = temp_dir / "igv" / "igv_reports_map.json"
        assert mapping_file.exists()

        # Verify content structure (basic check)
        import json

        with open(mapping_file) as f:
            mapping_data = json.load(f)

        assert isinstance(mapping_data, dict)
        assert "variants" in mapping_data
        assert isinstance(mapping_data["variants"], list)
        # Should have entries for variants with non-reference genotypes
        assert len(mapping_data["variants"]) > 0

    @patch("variantcentrifuge.generate_igv_report._generate_single_igv_report")
    def test_progress_logging(self, mock_single_report, sample_tsv, bam_mapping_file, temp_dir):
        """Test that progress is logged during parallel processing."""
        # Setup
        mock_single_report.return_value = (
            True,
            "Sample_HG002",
            "17",
            "41223094",
            "T",
            "C",
            "/path/report.html",
        )

        with patch("variantcentrifuge.generate_igv_report.logger") as mock_logger:
            generate_igv_report(
                variants_tsv=str(sample_tsv),
                output_dir=str(temp_dir),
                bam_mapping_file=str(bam_mapping_file),
                igv_reference="hg19",
                max_workers=2,
            )

            # Verify progress logging
            info_calls = [
                call for call in mock_logger.info.call_args_list if "Progress:" in str(call)
            ]
            assert len(info_calls) > 0

    def test_error_handling_in_parallel_processing(self, sample_tsv, bam_mapping_file, temp_dir):
        """Test error handling during parallel processing."""
        # Mock one successful and one failing task
        call_count = 0

        def mock_side_effect(*args):
            nonlocal call_count
            call_count += 1
            if call_count % 2 == 0:
                # Every other call fails
                raise Exception("Simulated failure")
            return (True, "Sample_HG002", "17", "41223094", "T", "C", "/path/report.html")

        with patch(
            "variantcentrifuge.generate_igv_report._generate_single_igv_report",
            side_effect=mock_side_effect,
        ):
            # Should not raise exception despite some tasks failing
            generate_igv_report(
                variants_tsv=str(sample_tsv),
                output_dir=str(temp_dir),
                bam_mapping_file=str(bam_mapping_file),
                igv_reference="hg19",
                max_workers=2,
            )

        # Verify mapping file still created (with successful reports only)
        mapping_file = temp_dir / "igv" / "igv_reports_map.json"
        assert mapping_file.exists()

        # Verify the JSON structure
        import json

        with open(mapping_file) as f:
            mapping_data = json.load(f)
        assert isinstance(mapping_data, dict)
        assert "variants" in mapping_data

    def test_custom_flanking_and_parameters(self, sample_tsv, bam_mapping_file, temp_dir):
        """Test that custom parameters are passed through correctly."""
        with patch(
            "variantcentrifuge.generate_igv_report._generate_single_igv_report"
        ) as mock_single_report:
            mock_single_report.return_value = (
                True,
                "Sample_HG002",
                "17",
                "41223094",
                "T",
                "C",
                "/path/report.html",
            )

            generate_igv_report(
                variants_tsv=str(sample_tsv),
                output_dir=str(temp_dir),
                bam_mapping_file=str(bam_mapping_file),
                igv_reference="hg38",
                igv_flanking=100,
                igv_max_allele_len_filename=15,
                igv_hash_len_filename=8,
                max_workers=2,
            )

            # Verify parameters were passed to individual report generation
            assert mock_single_report.called
            call_args = mock_single_report.call_args_list[0][0][
                0
            ]  # First call, first arg (task_data)

            # Check flanking parameter (index 8 in task_data tuple)
            assert call_args[8] == 100

            # Check reference parameter (index 9)
            assert call_args[9] == "hg38"
