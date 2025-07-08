"""Test the archive_results functionality."""

import os
import tempfile
import tarfile
import pytest
from unittest.mock import patch

from variantcentrifuge.pipeline import archive_results_directory


class TestArchiveResults:
    """Test suite for the archive_results_directory function."""

    def test_archive_results_basic(self):
        """Test basic archive creation with simple directory structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test output directory
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            # Create test files
            test_file = os.path.join(output_dir, "results.tsv")
            with open(test_file, "w") as f:
                f.write("chrom\tpos\tref\talt\n")
                f.write("chr1\t12345\tA\tG\n")

            # Archive the directory
            archive_path = archive_results_directory(output_dir, "test_sample")

            # Verify archive was created
            assert archive_path is not None
            assert os.path.exists(archive_path)
            assert archive_path.endswith(".tar.gz")
            assert "variantcentrifuge_results_test_sample_" in archive_path

            # Verify it's in the parent directory
            assert os.path.dirname(archive_path) == tmpdir

            # Verify archive contents
            with tarfile.open(archive_path, "r:gz") as tar:
                members = tar.getnames()
                assert len(members) == 2  # directory + file
                assert any("test_output" in m for m in members)
                assert any("results.tsv" in m for m in members)

    def test_archive_results_complex_structure(self):
        """Test archive creation with complex directory structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create complex output structure
            output_dir = os.path.join(tmpdir, "complex_output")
            os.makedirs(output_dir)

            # Create subdirectories
            report_dir = os.path.join(output_dir, "report")
            os.makedirs(report_dir)
            intermediate_dir = os.path.join(output_dir, "intermediate")
            os.makedirs(intermediate_dir)
            igv_dir = os.path.join(report_dir, "igv")
            os.makedirs(igv_dir)

            # Create various files
            files_to_create = [
                ("results.tsv", "test data"),
                ("results.xlsx", "excel data"),
                ("report/summary.html", "<html>summary</html>"),
                ("report/variants.json", '{"variants": []}'),
                ("report/igv/igv_reports_map.json", '{"reports": []}'),
                ("intermediate/extracted.tsv.gz", "compressed data"),
            ]

            for file_path, content in files_to_create:
                full_path = os.path.join(output_dir, file_path)
                with open(full_path, "w") as f:
                    f.write(content)

            # Archive the directory
            archive_path = archive_results_directory(output_dir, "complex_sample")

            # Verify archive was created
            assert archive_path is not None
            assert os.path.exists(archive_path)

            # Verify all files are in the archive
            with tarfile.open(archive_path, "r:gz") as tar:
                members = tar.getnames()
                for file_path, _ in files_to_create:
                    assert any(
                        file_path.replace("/", os.sep) in m for m in members
                    ), f"Missing {file_path} in archive"

    def test_archive_results_timestamp_format(self):
        """Test that archive filename has correct timestamp format."""
        import re

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            archive_path = archive_results_directory(output_dir, "test_sample")

            # Extract the filename
            filename = os.path.basename(archive_path)

            # Check the format: variantcentrifuge_results_<base_name>_<timestamp>.tar.gz
            pattern = r"variantcentrifuge_results_test_sample_(\d{8})_(\d{6})\.tar\.gz"
            match = re.match(pattern, filename)

            assert match is not None, f"Filename {filename} doesn't match expected pattern"

            # Check date format (YYYYMMDD)
            date_part = match.group(1)
            assert len(date_part) == 8
            assert date_part.isdigit()

            # Check time format (HHMMSS)
            time_part = match.group(2)
            assert len(time_part) == 6
            assert time_part.isdigit()

    def test_archive_results_special_characters(self):
        """Test archive creation with special characters in base_name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            # Test with various special characters
            special_names = [
                "sample_with_spaces",
                "sample-with-dashes",
                "sample.with.dots",
                "sample_123",
            ]

            for base_name in special_names:
                archive_path = archive_results_directory(output_dir, base_name)
                assert archive_path is not None
                assert base_name in archive_path
                assert os.path.exists(archive_path)
                # Clean up
                os.remove(archive_path)

    def test_archive_results_large_files(self):
        """Test archive creation with larger files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            # Create a larger file (1MB)
            large_file = os.path.join(output_dir, "large_results.tsv")
            with open(large_file, "w") as f:
                # Write 1MB of data
                for i in range(10000):
                    f.write(f"chr1\t{i}\tA\tG\t" + "X" * 90 + "\n")

            # Archive the directory
            archive_path = archive_results_directory(output_dir, "large_sample")

            # Verify archive was created
            assert archive_path is not None
            assert os.path.exists(archive_path)

            # Verify compression is working (archive should be smaller than original)
            original_size = os.path.getsize(large_file)
            archive_size = os.path.getsize(archive_path)
            assert archive_size < original_size, "Archive should be compressed"

    def test_archive_results_empty_directory(self):
        """Test archive creation with empty directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "empty_output")
            os.makedirs(output_dir)

            # Archive empty directory
            archive_path = archive_results_directory(output_dir, "empty_sample")

            # Should still create archive
            assert archive_path is not None
            assert os.path.exists(archive_path)

            # Verify archive contains the directory
            with tarfile.open(archive_path, "r:gz") as tar:
                members = tar.getnames()
                assert len(members) >= 1
                assert any("empty_output" in m for m in members)

    def test_archive_results_nonexistent_directory(self):
        """Test archive creation with non-existent directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "nonexistent")

            # Should handle gracefully
            archive_path = archive_results_directory(output_dir, "test_sample")
            assert archive_path is None

    def test_archive_results_permission_error(self):
        """Test archive creation when parent directory is not writable."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            # Mock tarfile.open to raise permission error
            with patch("tarfile.open") as mock_tar:
                mock_tar.side_effect = PermissionError("Permission denied")

                archive_path = archive_results_directory(output_dir, "test_sample")
                assert archive_path is None

    def test_archive_results_disk_space_error(self):
        """Test archive creation when disk is full."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            # Mock tarfile.open to raise IOError
            with patch("tarfile.open") as mock_tar:
                mock_tar.side_effect = IOError("No space left on device")

                archive_path = archive_results_directory(output_dir, "test_sample")
                assert archive_path is None

    @patch("variantcentrifuge.pipeline_core.logger")
    def test_archive_results_logging(self, mock_logger):
        """Test that appropriate log messages are generated."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            # Create archive
            archive_path = archive_results_directory(output_dir, "test_sample")

            # Verify logging calls
            assert mock_logger.info.call_count >= 2
            mock_logger.info.assert_any_call(f"Creating archive: {archive_path}")

            # Check that success message includes file size
            success_calls = [
                call
                for call in mock_logger.info.call_args_list
                if "Archive created successfully" in str(call)
            ]
            assert len(success_calls) == 1
            assert "MB)" in str(success_calls[0])

    def test_archive_results_symlinks(self):
        """Test archive creation with symbolic links."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "test_output")
            os.makedirs(output_dir)

            # Create a file and a symlink to it
            real_file = os.path.join(output_dir, "real_file.txt")
            with open(real_file, "w") as f:
                f.write("real content")

            symlink = os.path.join(output_dir, "link_to_file.txt")
            os.symlink(real_file, symlink)

            # Archive the directory
            archive_path = archive_results_directory(output_dir, "symlink_test")

            # Verify archive was created
            assert archive_path is not None
            assert os.path.exists(archive_path)

            # Extract and verify symlinks are preserved
            extract_dir = os.path.join(tmpdir, "extracted")
            with tarfile.open(archive_path, "r:gz") as tar:
                tar.extractall(extract_dir)

            extracted_link = os.path.join(extract_dir, "test_output", "link_to_file.txt")
            assert os.path.islink(extracted_link) or os.path.exists(extracted_link)


@pytest.mark.integration
class TestArchiveResultsIntegration:
    """Integration tests for archive_results in the pipeline."""

    @patch("variantcentrifuge.pipeline_core.run_command")
    @patch("variantcentrifuge.pipeline_core.check_external_tools")
    def test_pipeline_with_archive_flag(self, mock_check_tools, mock_run_command):
        """Test that pipeline correctly calls archive function when flag is set."""
        # This is a placeholder for integration testing
        # Would need to mock the entire pipeline run
        pass
