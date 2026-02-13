"""Tests for CLI debug logging functionality."""

import contextlib
import sys
from unittest.mock import patch

from variantcentrifuge.cli import main


class TestCLIDebugLogging:
    """Test CLI debug logging functionality."""

    @patch("variantcentrifuge.cli.run_refactored_pipeline")
    @patch("variantcentrifuge.cli.validate_vcf_file")
    @patch("variantcentrifuge.cli.validate_mandatory_parameters")
    def test_command_line_logging_in_debug_mode(
        self, mock_validate_params, mock_validate_vcf, mock_run_pipeline, caplog
    ):
        """Test that command line invocation is logged in debug mode."""
        import logging

        caplog.set_level(logging.DEBUG)

        # Mock sys.argv to simulate command line arguments
        test_args = [
            "variantcentrifuge",
            "--log-level",
            "DEBUG",
            "-v",
            "/tmp/test.vcf",
            "-r",
            "GRCh37",
            "--gene-name",
            "BRCA1",
            "--case-samples",
            "Sample 1,Sample2",
        ]

        # Mock the pipeline to avoid actually running it
        mock_run_pipeline.return_value = 0
        mock_validate_vcf.return_value = True
        mock_validate_params.return_value = True

        with patch.object(sys, "argv", test_args), contextlib.suppress(SystemExit):
            main()

        # Check that debug messages were logged
        debug_messages = [
            record.message for record in caplog.records if record.levelname == "DEBUG"
        ]

        # Find the command line invocation message
        command_line_messages = [msg for msg in debug_messages if "Command line invocation:" in msg]
        assert len(command_line_messages) == 1

        command_line_msg = command_line_messages[0]

        # Verify the command line was reconstructed correctly
        assert "variantcentrifuge" in command_line_msg
        assert "--log-level DEBUG" in command_line_msg
        assert "-v /tmp/test.vcf" in command_line_msg
        assert "-r GRCh37" in command_line_msg
        assert "--gene-name BRCA1" in command_line_msg
        # The space in "Sample 1,Sample2" should be properly quoted
        assert "'Sample 1,Sample2'" in command_line_msg or '"Sample 1,Sample2"' in command_line_msg

        # Check that Python executable was logged
        python_messages = [msg for msg in debug_messages if "Python executable:" in msg]
        assert len(python_messages) == 1
        assert sys.executable in python_messages[0]

        # Check that working directory was logged
        working_dir_messages = [msg for msg in debug_messages if "Working directory:" in msg]
        assert len(working_dir_messages) == 1

    @patch("variantcentrifuge.cli.run_refactored_pipeline")
    @patch("variantcentrifuge.cli.validate_vcf_file")
    @patch("variantcentrifuge.cli.validate_mandatory_parameters")
    def test_no_command_line_logging_in_info_mode(
        self, mock_validate_params, mock_validate_vcf, mock_run_pipeline, caplog
    ):
        """Test that command line invocation is NOT logged in INFO mode."""
        import logging

        caplog.set_level(logging.INFO)

        # Mock sys.argv to simulate command line arguments
        test_args = [
            "variantcentrifuge",
            "--log-level",
            "INFO",  # INFO level, not DEBUG
            "-v",
            "/tmp/test.vcf",
            "-r",
            "GRCh37",
        ]

        # Mock the pipeline to avoid actually running it
        mock_run_pipeline.return_value = 0
        mock_validate_vcf.return_value = True
        mock_validate_params.return_value = True

        with patch.object(sys, "argv", test_args), contextlib.suppress(SystemExit):
            main()

        # Check that debug messages were NOT logged
        debug_messages = [
            record.message for record in caplog.records if record.levelname == "DEBUG"
        ]

        # No debug messages should be present in INFO mode
        command_line_messages = [msg for msg in debug_messages if "Command line invocation:" in msg]
        assert len(command_line_messages) == 0

        python_messages = [msg for msg in debug_messages if "Python executable:" in msg]
        assert len(python_messages) == 0

        working_dir_messages = [msg for msg in debug_messages if "Working directory:" in msg]
        assert len(working_dir_messages) == 0
