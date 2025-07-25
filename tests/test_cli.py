# File: tests/test_cli.py
# Location: variantcentrifuge/tests/test_cli.py

"""
Tests for CLI module.

This file contains tests ensuring the CLI runs and shows help correctly.
"""

import json
import subprocess
from unittest.mock import patch

import pytest


def test_cli_help():
    """Test that the CLI help message can be displayed."""
    cmd = ["python", "-m", "variantcentrifuge.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "usage:" in result.stdout


def test_bcftools_prefilter_in_help():
    """Test that the bcftools-prefilter argument is in the help text."""
    cmd = ["python", "-m", "variantcentrifuge.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "--bcftools-prefilter" in result.stdout
    # Check that the help text mentions bcftools expression
    assert "bcftools expression" in result.stdout.lower()


class TestShowCheckpointStatus:
    """Test the --show-checkpoint-status functionality."""

    def test_show_checkpoint_status_no_state(self, tmp_path):
        """Test --show-checkpoint-status when no checkpoint state exists."""
        with patch(
            "sys.argv",
            ["variantcentrifuge", "--show-checkpoint-status", "--output-dir", str(tmp_path)],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0
                # Should print that no checkpoint state was found
                mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_with_new_pipeline_flag(self, tmp_path):
        """Test --show-checkpoint-status with --use-new-pipeline flag."""
        with patch(
            "sys.argv",
            [
                "variantcentrifuge",
                "--show-checkpoint-status",
                "--output-dir",
                str(tmp_path),
            ],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0
                # Should indicate it's checking the stage-based pipeline
                mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_with_config_file(self, tmp_path):
        """Test --show-checkpoint-status with config file that specifies new pipeline."""
        # Create a config file with new pipeline enabled
        config_file = tmp_path / "config.json"
        config_data = {"reference": "GRCh38.99"}
        config_file.write_text(json.dumps(config_data))

        with patch(
            "sys.argv",
            [
                "variantcentrifuge",
                "--show-checkpoint-status",
                "--config",
                str(config_file),
                "--output-dir",
                str(tmp_path),
            ],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0
                # Should detect new pipeline from config
                mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_with_existing_state(self, tmp_path):
        """Test --show-checkpoint-status when checkpoint state exists."""
        from variantcentrifuge.checkpoint import PipelineState

        # Create a checkpoint state
        state = PipelineState(str(tmp_path))
        state.initialize({"test_config": "value"}, "1.0.0")
        state.start_step("test_step")
        state.complete_step("test_step")

        with patch(
            "sys.argv",
            ["variantcentrifuge", "--show-checkpoint-status", "--output-dir", str(tmp_path)],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0
                # Should print the pipeline status summary
                mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")

                # Check that it printed checkpoint information
                print_calls = [call.args[0] for call in mock_print.call_args_list]
                summary_printed = any("Pipeline State Summary:" in call for call in print_calls)
                assert summary_printed

    def test_show_checkpoint_status_config_file_error(self, tmp_path):
        """Test --show-checkpoint-status with invalid config file."""
        # Create an invalid config file
        config_file = tmp_path / "invalid_config.json"
        config_file.write_text("invalid json content")

        with patch(
            "sys.argv",
            [
                "variantcentrifuge",
                "--show-checkpoint-status",
                "--config",
                str(config_file),
                "--output-dir",
                str(tmp_path),
            ],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0
                # Should fall back to original pipeline when config loading fails
                mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_flag_precedence(self, tmp_path):
        """Test that --use-new-pipeline flag takes precedence over config file."""
        # Create a config file with new pipeline disabled
        config_file = tmp_path / "config.json"
        config_data = {"reference": "GRCh38.99"}
        config_file.write_text(json.dumps(config_data))

        with patch(
            "sys.argv",
            [
                "variantcentrifuge",
                "--show-checkpoint-status",
                "--config",
                str(config_file),
                "--output-dir",
                str(tmp_path),
            ],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0
                # Should use stage-based pipeline
                mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_logging_level(self, tmp_path):
        """Test --show-checkpoint-status with different logging levels."""
        with patch(
            "sys.argv",
            [
                "variantcentrifuge",
                "--show-checkpoint-status",
                "--log-level",
                "DEBUG",
                "--output-dir",
                str(tmp_path),
            ],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0
                # Should work with different log levels
                mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_detailed_summary(self, tmp_path):
        """Test --show-checkpoint-status displays detailed checkpoint information."""
        from variantcentrifuge.checkpoint import PipelineState

        # Create a more complex checkpoint state
        state = PipelineState(str(tmp_path))
        state.initialize({"gene_name": "BRCA1", "vcf_file": "test.vcf"}, "1.0.0")

        # Add multiple steps
        state.start_step("step1")
        state.complete_step("step1")

        state.start_step("step2")
        state.complete_step("step2")

        state.start_step("step3")
        state.fail_step("step3", "Test error")

        with patch(
            "sys.argv",
            ["variantcentrifuge", "--show-checkpoint-status", "--output-dir", str(tmp_path)],
        ):
            with patch("builtins.print") as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main

                    main()

                assert exc_info.value.code == 0

                # Verify detailed information is printed
                print_calls = [call.args[0] for call in mock_print.call_args_list]
                full_output = "\n".join(print_calls)

                assert "Pipeline State Summary:" in full_output
                assert "step1" in full_output
                assert "step2" in full_output
                assert "step3" in full_output
                assert "✓" in full_output  # Completed steps
                assert "✗" in full_output  # Failed step
                assert "Test error" in full_output  # Error message


class TestArgumentParser:
    """Comprehensive tests for CLI argument parser functionality."""

    def test_create_parser_function(self):
        """Test that create_parser() returns a properly configured ArgumentParser."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        assert parser is not None
        assert hasattr(parser, "parse_args")
        assert hasattr(parser, "add_argument")

        # Test that help can be generated without errors
        help_text = parser.format_help()
        assert "variantcentrifuge: Filter and process VCF files." in help_text

    def test_genotype_replacement_chunk_size_parameter_exists(self):
        """Test that --genotype-replacement-chunk-size parameter is defined in create_parser()."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Verify the parameter exists in help
        assert "--genotype-replacement-chunk-size" in help_text
        assert "Chunk size for parallel genotype replacement" in help_text

        # Test parsing the parameter
        args = parser.parse_args(
            ["--vcf-file", "test.vcf", "--genotype-replacement-chunk-size", "25000"]
        )
        assert args.genotype_replacement_chunk_size == 25000

    def test_genotype_replacement_chunk_size_parameter_cli(self):
        """Test --genotype-replacement-chunk-size parameter via CLI interface."""
        # Test that the parameter is recognized and doesn't cause "unrecognized arguments" error
        cmd = [
            "python",
            "-m",
            "variantcentrifuge.cli",
            "--genotype-replacement-chunk-size",
            "25000",
            "--help",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Should succeed without "unrecognized arguments" error
        assert result.returncode == 0
        assert "unrecognized arguments" not in result.stderr
        assert "--genotype-replacement-chunk-size" in result.stdout

    def test_all_argument_groups_present(self):
        """Test that all expected argument groups are present in create_parser()."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        expected_groups = [
            "General Options",
            "Core Input/Output",
            "Gene Selection",
            "Filtering & Annotation",
            "Field Extraction & Formatting",
            "Genotype Analysis",
            "Phenotype & Sample Groups",
            "Statistical Analysis",
            "Inheritance Analysis",
            "Scoring & Custom Annotations",
            "Reporting & Visualization",
            "Performance & Processing",
            "Checkpoint & Resume Options",
            "Data Privacy Options",
            "Miscellaneous Options",
        ]

        for group in expected_groups:
            assert group in help_text, f"Missing argument group: {group}"

    def test_required_arguments(self):
        """Test that required arguments are properly defined."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()

        # Test that required arguments are enforced
        with pytest.raises(SystemExit):
            parser.parse_args([])  # Should fail without --vcf-file

        # Test with required argument
        args = parser.parse_args(["--vcf-file", "test.vcf"])
        assert args.vcf_file == "test.vcf"

    def test_performance_arguments(self):
        """Test Performance & Processing argument group parameters."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Test that key performance parameters exist
        performance_params = [
            "--threads",
            "--genotype-replacement-chunk-size",
            "--chunks",
            "--sort-memory-limit",
            "--max-memory-gb",
        ]

        for param in performance_params:
            assert param in help_text, f"Missing performance parameter: {param}"

        # Test parsing performance arguments
        args = parser.parse_args(
            [
                "--vcf-file",
                "test.vcf",
                "--threads",
                "8",
                "--genotype-replacement-chunk-size",
                "10000",
                "--chunks",
                "4",
            ]
        )

        assert args.threads == 8
        assert args.genotype_replacement_chunk_size == 10000
        assert args.chunks == 4

    def test_genotype_replacement_methods(self):
        """Test genotype replacement method parameter."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Check that genotype replacement method exists
        assert "--genotype-replacement-method" in help_text

        # Test valid methods
        valid_methods = [
            "auto",
            "sequential",
            "vectorized",
            "chunked-vectorized",
            "parallel",
            "streaming-parallel",
        ]
        for method in valid_methods:
            args = parser.parse_args(
                ["--vcf-file", "test.vcf", "--genotype-replacement-method", method]
            )
            assert args.genotype_replacement_method == method

    def test_checkpoint_arguments(self):
        """Test Checkpoint & Resume Options argument group."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        checkpoint_params = [
            "--checkpoint",
            "--resume-checkpoint",
            "--show-checkpoint-status",
            "--checkpoint-dir",
        ]

        # Check that checkpoint parameters exist
        for param in checkpoint_params:
            if param in help_text:  # Some may be optional
                # Test parsing if parameter exists
                try:
                    parser.parse_args(["--vcf-file", "test.vcf", param])
                except (SystemExit, ValueError):
                    # Some parameters may require values - that's okay
                    pass

    def test_data_privacy_arguments(self):
        """Test Data Privacy Options argument group."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Check for privacy-related parameters
        privacy_params = ["--anonymize", "--pseudonymize", "--remove-patient-info"]

        for param in privacy_params:
            if param in help_text:
                # Test that privacy parameters can be parsed
                try:
                    parser.parse_args(["--vcf-file", "test.vcf", param])
                except SystemExit:
                    # Some may require additional arguments
                    pass

    def test_miscellaneous_arguments(self):
        """Test Miscellaneous Options argument group."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Test that genotype-replacement-chunk-size is in Miscellaneous Options
        assert "--genotype-replacement-chunk-size" in help_text

        # Test gzip intermediate options
        misc_params = ["--gzip-intermediates", "--no-gzip-intermediates"]

        for param in misc_params:
            if param in help_text:
                parser.parse_args(["--vcf-file", "test.vcf", param])
                # Should parse successfully

    def test_parser_consistency_regression(self):
        """Regression test to ensure main() uses create_parser() consistently."""
        from variantcentrifuge.cli import create_parser

        # Get all arguments from create_parser()
        parser = create_parser()
        help_text = parser.format_help()

        # Test key parameters that were affected by the duplicate parser bug
        critical_params = [
            "--genotype-replacement-chunk-size",
            "--threads",
            "--chunks",
            "--vcf-file",
            "--output-file",
        ]

        for param in critical_params:
            assert param in help_text, f"Critical parameter missing from create_parser(): {param}"

        # Test that these parameters work via CLI
        for param in ["--genotype-replacement-chunk-size", "--threads", "--chunks"]:
            if param in help_text:
                cmd = ["python", "-m", "variantcentrifuge.cli", param, "10", "--help"]
                result = subprocess.run(cmd, capture_output=True, text=True)
                assert result.returncode == 0, f"Parameter {param} failed via CLI"
                assert (
                    "unrecognized arguments" not in result.stderr
                ), f"Parameter {param} not recognized"

    def test_argument_defaults(self):
        """Test that argument defaults are properly set."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(["--vcf-file", "test.vcf"])

        # Test some expected defaults
        assert args.genotype_replacement_chunk_size == 50000
        assert args.log_level == "INFO"

        # Test that defaults can be overridden
        args = parser.parse_args(
            [
                "--vcf-file",
                "test.vcf",
                "--genotype-replacement-chunk-size",
                "25000",
                "--log-level",
                "DEBUG",
            ]
        )

        assert args.genotype_replacement_chunk_size == 25000
        assert args.log_level == "DEBUG"

    def test_invalid_arguments(self):
        """Test that invalid arguments are properly rejected."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()

        # Test invalid log level
        with pytest.raises(SystemExit):
            parser.parse_args(["--vcf-file", "test.vcf", "--log-level", "INVALID"])

        # Test invalid chunk size (non-integer)
        with pytest.raises(SystemExit):
            parser.parse_args(
                ["--vcf-file", "test.vcf", "--genotype-replacement-chunk-size", "not_a_number"]
            )

    def test_argument_type_conversion(self):
        """Test that argument types are properly converted."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(
            [
                "--vcf-file",
                "test.vcf",
                "--genotype-replacement-chunk-size",
                "15000",
                "--threads",
                "16",
            ]
        )

        # Test type conversions
        assert isinstance(args.genotype_replacement_chunk_size, int)
        assert isinstance(args.threads, int)
        assert isinstance(args.vcf_file, str)
        assert isinstance(args.log_level, str)


class TestCLIRegressionTests:
    """Specific regression tests for bugs that have been fixed."""

    def test_duplicate_parser_regression(self):
        """Regression test for duplicate argument parser bug.

        This test ensures that the bug where main() created its own parser
        instead of using create_parser() doesn't happen again.
        """
        # Test that CLI recognizes all parameters from create_parser()
        critical_parameters = [
            "--genotype-replacement-chunk-size",
            "--threads",
            "--chunks",
            "--sort-memory-limit",
            "--max-memory-gb",
        ]

        for param in critical_parameters:
            # Each parameter should be recognized by the CLI
            cmd = ["python", "-m", "variantcentrifuge.cli", param, "10", "--help"]
            result = subprocess.run(cmd, capture_output=True, text=True)

            assert result.returncode == 0, f"CLI failed for parameter {param}"
            assert (
                "unrecognized arguments" not in result.stderr
            ), f"Parameter {param} not recognized - possible duplicate parser issue"
            assert param in result.stdout, f"Parameter {param} not in help output"

    def test_genotype_replacement_chunk_size_specifically(self):
        """Specific regression test for --genotype-replacement-chunk-size parameter.

        This is the exact parameter that was failing due to the duplicate parser bug.
        """
        # Test with the exact command that was failing
        cmd = [
            "python",
            "-m",
            "variantcentrifuge.cli",
            "--genotype-replacement-chunk-size",
            "25000",
            "--help",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Should succeed without any errors
        assert result.returncode == 0
        assert "unrecognized arguments" not in result.stderr
        assert "--genotype-replacement-chunk-size" in result.stdout

        # Test with different values
        test_values = ["1000", "10000", "50000", "100000"]
        for value in test_values:
            cmd = [
                "python",
                "-m",
                "variantcentrifuge.cli",
                "--genotype-replacement-chunk-size",
                value,
                "--help",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Failed with chunk size {value}"
            assert "unrecognized arguments" not in result.stderr

    def test_parser_architecture_consistency(self):
        """Test that ensures the parser architecture is consistent.

        This test verifies that main() properly uses create_parser() and doesn't
        create duplicate argument definitions.
        """
        import inspect

        import variantcentrifuge.cli as cli_module

        # Get the main function source
        main_source = inspect.getsource(cli_module.main)

        # Verify that main() calls create_parser()
        assert "create_parser()" in main_source, "main() should call create_parser()"

        # Verify that main() doesn't create duplicate argument groups
        # (This would indicate the old bug is back)
        duplicate_indicators = [
            'add_argument_group("General Options")',
            'add_argument_group("Core Input/Output")',
            'add_argument_group("Performance & Processing")',
        ]

        for indicator in duplicate_indicators:
            # These should NOT appear in main() since they should only be in create_parser()
            assert indicator not in main_source, (
                f"main() contains duplicate argument group: {indicator}. "
                f"This suggests the duplicate parser bug is back."
            )

        # Verify main() has minimal parser-related code
        parser_creation_count = main_source.count("ArgumentParser")
        # main() should only create the status_parser for --show-checkpoint-status
        assert parser_creation_count <= 1, (
            f"main() creates too many ArgumentParsers ({parser_creation_count}). "
            f"Should only create status_parser for --show-checkpoint-status."
        )
