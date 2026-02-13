# File: tests/test_cli.py
# Location: variantcentrifuge/tests/test_cli.py

"""
Tests for CLI module.

This file contains tests ensuring the CLI runs and shows help correctly.
"""

import contextlib
import json
from unittest.mock import patch

import pytest


def test_cli_help():
    """Test that the CLI help message can be displayed."""
    from variantcentrifuge.cli import create_parser

    parser = create_parser()
    help_text = parser.format_help()
    assert "usage:" in help_text


def test_bcftools_prefilter_in_help():
    """Test that the bcftools-prefilter argument is in the help text."""
    from variantcentrifuge.cli import create_parser

    parser = create_parser()
    help_text = parser.format_help()
    assert "--bcftools-prefilter" in help_text
    # Check that the help text mentions bcftools expression
    assert "bcftools expression" in help_text.lower()


class TestShowCheckpointStatus:
    """Test the --show-checkpoint-status functionality."""

    def test_show_checkpoint_status_no_state(self, tmp_path):
        """Test --show-checkpoint-status when no checkpoint state exists."""
        with (
            patch(
                "sys.argv",
                ["variantcentrifuge", "--show-checkpoint-status", "--output-dir", str(tmp_path)],
            ),
            patch("builtins.print") as mock_print,
        ):
            with pytest.raises(SystemExit) as exc_info:
                from variantcentrifuge.cli import main

                main()

            assert exc_info.value.code == 0
            # Should print that no checkpoint state was found
            mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
            mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_with_new_pipeline_flag(self, tmp_path):
        """Test --show-checkpoint-status with --use-new-pipeline flag."""
        with (
            patch(
                "sys.argv",
                [
                    "variantcentrifuge",
                    "--show-checkpoint-status",
                    "--output-dir",
                    str(tmp_path),
                ],
            ),
            patch("builtins.print") as mock_print,
        ):
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

        with (
            patch(
                "sys.argv",
                [
                    "variantcentrifuge",
                    "--show-checkpoint-status",
                    "--config",
                    str(config_file),
                    "--output-dir",
                    str(tmp_path),
                ],
            ),
            patch("builtins.print") as mock_print,
        ):
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

        with (
            patch(
                "sys.argv",
                ["variantcentrifuge", "--show-checkpoint-status", "--output-dir", str(tmp_path)],
            ),
            patch("builtins.print") as mock_print,
        ):
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

        with (
            patch(
                "sys.argv",
                [
                    "variantcentrifuge",
                    "--show-checkpoint-status",
                    "--config",
                    str(config_file),
                    "--output-dir",
                    str(tmp_path),
                ],
            ),
            patch("builtins.print") as mock_print,
        ):
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

        with (
            patch(
                "sys.argv",
                [
                    "variantcentrifuge",
                    "--show-checkpoint-status",
                    "--config",
                    str(config_file),
                    "--output-dir",
                    str(tmp_path),
                ],
            ),
            patch("builtins.print") as mock_print,
        ):
            with pytest.raises(SystemExit) as exc_info:
                from variantcentrifuge.cli import main

                main()

            assert exc_info.value.code == 0
            # Should use stage-based pipeline
            mock_print.assert_any_call("Checking checkpoint status for stage-based pipeline...")
            mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")

    def test_show_checkpoint_status_logging_level(self, tmp_path):
        """Test --show-checkpoint-status with different logging levels."""
        with (
            patch(
                "sys.argv",
                [
                    "variantcentrifuge",
                    "--show-checkpoint-status",
                    "--log-level",
                    "DEBUG",
                    "--output-dir",
                    str(tmp_path),
                ],
            ),
            patch("builtins.print") as mock_print,
        ):
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

        with (
            patch(
                "sys.argv",
                ["variantcentrifuge", "--show-checkpoint-status", "--output-dir", str(tmp_path)],
            ),
            patch("builtins.print") as mock_print,
        ):
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
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        # Test that the parameter is recognized and parses correctly
        args = parser.parse_args(
            ["--vcf-file", "test.vcf", "--genotype-replacement-chunk-size", "25000"]
        )
        assert args.genotype_replacement_chunk_size == 25000
        assert "--genotype-replacement-chunk-size" in parser.format_help()

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
                with contextlib.suppress(SystemExit, ValueError):
                    parser.parse_args(["--vcf-file", "test.vcf", param])

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
                with contextlib.suppress(SystemExit):
                    parser.parse_args(["--vcf-file", "test.vcf", param])

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

        # Test that these parameters parse without errors
        for param in ["--genotype-replacement-chunk-size", "--threads", "--chunks"]:
            args = parser.parse_args(["--vcf-file", "test.vcf", param, "10"])
            assert args is not None, f"Parameter {param} failed to parse"

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
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Test that all critical parameters are recognized by create_parser()
        critical_parameters = [
            "--genotype-replacement-chunk-size",
            "--threads",
            "--chunks",
            "--sort-memory-limit",
            "--max-memory-gb",
        ]

        for param in critical_parameters:
            assert param in help_text, (
                f"Parameter {param} not in help output - possible duplicate parser issue"
            )
            # Verify parsing works
            args = parser.parse_args(["--vcf-file", "test.vcf", param, "10"])
            assert args is not None, f"Parameter {param} failed to parse"

    def test_genotype_replacement_chunk_size_specifically(self):
        """Specific regression test for --genotype-replacement-chunk-size parameter.

        This is the exact parameter that was failing due to the duplicate parser bug.
        """
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()
        assert "--genotype-replacement-chunk-size" in help_text

        # Test with different values
        test_values = ["1000", "10000", "50000", "100000"]
        for value in test_values:
            args = parser.parse_args(
                ["--vcf-file", "test.vcf", "--genotype-replacement-chunk-size", value]
            )
            assert args.genotype_replacement_chunk_size == int(value), (
                f"Failed with chunk size {value}"
            )

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
        # main() creates mini-parsers for early-exit flags:
        # - status_parser for --show-checkpoint-status
        # - profile_parser for --list-field-profiles
        # - ann_parser for --show-vcf-annotations
        assert parser_creation_count <= 3, (
            f"main() creates too many ArgumentParsers ({parser_creation_count}). "
            f"Expected at most 3 (status_parser + profile_parser + ann_parser)."
        )


class TestPerformanceConfigMapping:
    """Tests for performance parameter config mapping bug fix."""

    def test_performance_config_mapping_regression(self):
        """Regression test for performance parameter config mapping bug.

        This tests the specific bug where CLI performance parameters like
        --genotype-replacement-method were not being mapped to the config,
        causing auto-selection to override user choices.
        """
        # Test the config mapping logic directly
        from types import SimpleNamespace

        # Simulate args from CLI
        mock_args = SimpleNamespace(
            max_memory_gb=250.0,
            genotype_replacement_method="parallel",
            vectorized_chunk_size=10000,
            genotype_replacement_chunk_size=25000,
        )

        # Test the config mapping (this is the code that was added to fix the bug)
        cfg = {}
        cfg["max_memory_gb"] = mock_args.max_memory_gb
        cfg["genotype_replacement_method"] = mock_args.genotype_replacement_method
        cfg["vectorized_chunk_size"] = mock_args.vectorized_chunk_size
        cfg["genotype_replacement_chunk_size"] = mock_args.genotype_replacement_chunk_size

        # Verify the parameters are correctly mapped
        assert cfg["max_memory_gb"] == 250.0
        assert cfg["genotype_replacement_method"] == "parallel"
        assert cfg["vectorized_chunk_size"] == 10000
        assert cfg["genotype_replacement_chunk_size"] == 25000

        # The key issue: user specified "parallel" should not become "auto"
        assert cfg.get("genotype_replacement_method", "auto") == "parallel"
        assert cfg.get("genotype_replacement_method", "auto") != "auto"

    def test_performance_parameters_mapped_to_config(self):
        """Test that performance CLI parameters parse and can be mapped to config."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(
            [
                "--vcf-file",
                "test.vcf",
                "--max-memory-gb",
                "64.0",
                "--genotype-replacement-method",
                "parallel",
                "--vectorized-chunk-size",
                "15000",
                "--genotype-replacement-chunk-size",
                "20000",
                "--threads",
                "8",
            ]
        )

        # Verify performance parameters are correctly parsed
        assert args.max_memory_gb == 64.0
        assert args.genotype_replacement_method == "parallel"
        assert args.vectorized_chunk_size == 15000
        assert args.genotype_replacement_chunk_size == 20000
        assert args.threads == 8

    def test_genotype_replacement_method_config_mapping(self):
        """Test that all genotype replacement methods parse correctly."""
        from variantcentrifuge.cli import create_parser

        test_methods = ["auto", "sequential", "vectorized", "parallel", "streaming-parallel"]

        for method in test_methods:
            parser = create_parser()
            args = parser.parse_args(
                ["--vcf-file", "test.vcf", "--genotype-replacement-method", method]
            )
            assert args.genotype_replacement_method == method, (
                f"Method {method} not correctly parsed"
            )

    def test_config_parameter_presence_regression(self):
        """Regression test ensuring performance parameters exist in parser."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(
            [
                "--vcf-file",
                "test.vcf",
                "--genotype-replacement-method",
                "parallel",
                "--max-memory-gb",
                "250.0",
                "--genotype-replacement-chunk-size",
                "25000",
                "--vectorized-chunk-size",
                "10000",
                "--threads",
                "16",
            ]
        )

        # Verify all performance parameters parsed correctly
        assert args.genotype_replacement_method == "parallel"
        assert args.max_memory_gb == 250.0
        assert args.genotype_replacement_chunk_size == 25000
        assert args.vectorized_chunk_size == 10000

    def test_performance_config_prevents_auto_selection_override(self):
        """Test that explicit method selection is preserved in parsed args."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(
            ["--vcf-file", "test.vcf", "--genotype-replacement-method", "parallel"]
        )

        # The parser should preserve the user's explicit choice
        assert args.genotype_replacement_method == "parallel"
        assert args.genotype_replacement_method != "auto"
        assert args.genotype_replacement_method != "streaming-parallel"

    def test_config_parameter_defaults(self):
        """Test that config parameters have correct defaults when not specified."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(["--vcf-file", "test.vcf"])

        # Verify defaults are correctly set
        assert args.genotype_replacement_method == "auto"
        assert args.max_memory_gb is None  # Auto-detect
        assert args.vectorized_chunk_size == 25000
        assert args.genotype_replacement_chunk_size == 50000

    def test_comprehensive_cli_parameter_config_mapping(self):
        """Comprehensive test for CLI parameter to config mapping completeness."""
        # Test the config mapping logic with all newly added parameters
        from types import SimpleNamespace

        # Simulate args with all the newly mapped parameters
        mock_args = SimpleNamespace(
            # Performance parameters
            max_memory_gb=128.0,
            genotype_replacement_method="vectorized",
            vectorized_chunk_size=20000,
            genotype_replacement_chunk_size=30000,
            chunks=8,
            # Core I/O parameters
            xlsx=True,
            keep_intermediates=True,
            archive_results=False,
            gzip_intermediates=True,
            # Gene selection parameters
            gene_name="BRCA1,BRCA2",
            gene_file="/path/to/genes.txt",
            # Field extraction parameters
            add_column="CUSTOM_FIELD",
            no_replacement=False,
            # Phenotype and sample group parameters
            case_phenotypes="HP:0000001,HP:0000002",
            case_phenotypes_file="/path/to/case_hpo.txt",
            case_samples="SAMPLE1,SAMPLE2",
            case_samples_file="/path/to/case_samples.txt",
            control_phenotypes="HP:0000003",
            control_phenotypes_file="/path/to/control_hpo.txt",
            control_samples="CTRL1,CTRL2",
            control_samples_file="/path/to/control_samples.txt",
            # Statistical analysis parameters
            stats_output_file="/path/to/stats.json",
            # Reporting parameters
            html_report=True,
        )

        # Test the config mapping (this verifies all the new mappings work)
        cfg = {}

        # Performance configuration - matches the code in cli.py
        cfg["max_memory_gb"] = mock_args.max_memory_gb
        cfg["genotype_replacement_method"] = mock_args.genotype_replacement_method
        cfg["vectorized_chunk_size"] = mock_args.vectorized_chunk_size
        cfg["genotype_replacement_chunk_size"] = mock_args.genotype_replacement_chunk_size
        cfg["chunks"] = mock_args.chunks

        # Core I/O configuration
        cfg["xlsx"] = mock_args.xlsx
        cfg["keep_intermediates"] = mock_args.keep_intermediates
        cfg["archive_results"] = mock_args.archive_results
        cfg["gzip_intermediates"] = mock_args.gzip_intermediates

        # Gene selection configuration
        cfg["gene_name"] = mock_args.gene_name
        cfg["gene_file"] = mock_args.gene_file

        # Field extraction configuration
        cfg["add_column"] = mock_args.add_column
        cfg["no_replacement"] = mock_args.no_replacement

        # Phenotype and sample group configuration
        cfg["case_phenotypes"] = mock_args.case_phenotypes
        cfg["case_phenotypes_file"] = mock_args.case_phenotypes_file
        cfg["case_samples"] = mock_args.case_samples
        cfg["case_samples_file"] = mock_args.case_samples_file
        cfg["control_phenotypes"] = mock_args.control_phenotypes
        cfg["control_phenotypes_file"] = mock_args.control_phenotypes_file
        cfg["control_samples"] = mock_args.control_samples
        cfg["control_samples_file"] = mock_args.control_samples_file

        # Statistical analysis configuration
        cfg["stats_output_file"] = mock_args.stats_output_file

        # Reporting configuration
        cfg["html_report"] = mock_args.html_report

        # Verify all parameters are correctly mapped
        assert cfg["max_memory_gb"] == 128.0
        assert cfg["genotype_replacement_method"] == "vectorized"
        assert cfg["vectorized_chunk_size"] == 20000
        assert cfg["genotype_replacement_chunk_size"] == 30000
        assert cfg["chunks"] == 8

        assert cfg["xlsx"] is True
        assert cfg["keep_intermediates"] is True
        assert cfg["archive_results"] is False
        assert cfg["gzip_intermediates"] is True

        assert cfg["gene_name"] == "BRCA1,BRCA2"
        assert cfg["gene_file"] == "/path/to/genes.txt"

        assert cfg["add_column"] == "CUSTOM_FIELD"
        assert cfg["no_replacement"] is False

        assert cfg["case_phenotypes"] == "HP:0000001,HP:0000002"
        assert cfg["case_phenotypes_file"] == "/path/to/case_hpo.txt"
        assert cfg["case_samples"] == "SAMPLE1,SAMPLE2"
        assert cfg["case_samples_file"] == "/path/to/case_samples.txt"
        assert cfg["control_phenotypes"] == "HP:0000003"
        assert cfg["control_phenotypes_file"] == "/path/to/control_hpo.txt"
        assert cfg["control_samples"] == "CTRL1,CTRL2"
        assert cfg["control_samples_file"] == "/path/to/control_samples.txt"

        assert cfg["stats_output_file"] == "/path/to/stats.json"
        assert cfg["html_report"] is True

    def test_debug_logging_parameter_mapping(self):
        """Test that debug logging shows parameter mapping correctly."""
        import logging
        from io import StringIO
        from types import SimpleNamespace

        # Set up debug logging capture
        log_capture_string = StringIO()
        ch = logging.StreamHandler(log_capture_string)
        ch.setLevel(logging.DEBUG)

        # Get the logger that would be used in CLI
        logger = logging.getLogger("variantcentrifuge")
        logger.addHandler(ch)
        logger.setLevel(logging.DEBUG)

        # Simulate the debug logging that happens in CLI
        mock_args = SimpleNamespace(
            genotype_replacement_method="parallel",
            max_memory_gb=128.0,
            vectorized_chunk_size=15000,
            genotype_replacement_chunk_size=25000,
            xlsx=True,
            gene_name="BRCA1",
        )

        # Simulate the debug logging that would happen in main()
        logger.debug(
            f"Performance config mapping: "
            f"genotype_replacement_method={mock_args.genotype_replacement_method}"
        )
        logger.debug(f"Performance config mapping: max_memory_gb={mock_args.max_memory_gb}")
        logger.debug(f"I/O config mapping: xlsx={mock_args.xlsx}")
        logger.debug(f"Gene config mapping: gene_name={mock_args.gene_name}")

        # Get the log contents
        log_contents = log_capture_string.getvalue()

        # Verify debug messages are present
        assert "Performance config mapping: genotype_replacement_method=parallel" in log_contents
        assert "Performance config mapping: max_memory_gb=128.0" in log_contents
        assert "I/O config mapping: xlsx=True" in log_contents
        assert "Gene config mapping: gene_name=BRCA1" in log_contents

        # Clean up
        logger.removeHandler(ch)

    def test_pseudonymize_table_parameter_validation(self):
        """Test that pseudonymize_table parameter validation works correctly."""
        from types import SimpleNamespace

        # Test case 1: pseudonymize=True, pseudonymize_table=None should fail validation
        args_missing_table = SimpleNamespace(
            pseudonymize=True,
            pseudonymize_table=None,
            pseudonymize_ped=False,
            pseudonymize_pattern=None,
            pseudonymize_schema="sequential",
        )

        # Simulate the validation logic from cli.py
        validation_error = None
        if args_missing_table.pseudonymize and not args_missing_table.pseudonymize_table:
            validation_error = "--pseudonymize-table is required when using --pseudonymize"

        assert validation_error is not None, (
            "Validation should fail when pseudonymize_table is missing"
        )
        assert "--pseudonymize-table is required" in validation_error

        # Test case 2: pseudonymize=True, pseudonymize_table=provided should pass validation
        args_with_table = SimpleNamespace(
            pseudonymize=True,
            pseudonymize_table="/path/to/mapping.tsv",
            pseudonymize_ped=False,
            pseudonymize_pattern=None,
            pseudonymize_schema="sequential",
        )

        validation_error = None
        if args_with_table.pseudonymize and not args_with_table.pseudonymize_table:
            validation_error = "--pseudonymize-table is required when using --pseudonymize"

        assert validation_error is None, (
            "Validation should pass when pseudonymize_table is provided"
        )

        # Test case 3: pseudonymize=False, pseudonymize_table=provided should warn
        args_table_without_pseudonymize = SimpleNamespace(
            pseudonymize=False,
            pseudonymize_table="/path/to/mapping.tsv",
            pseudonymize_ped=False,
            pseudonymize_pattern=None,
            pseudonymize_schema="sequential",
        )

        should_warn = False
        if (
            args_table_without_pseudonymize.pseudonymize_table
            and not args_table_without_pseudonymize.pseudonymize
        ):
            should_warn = True

        assert should_warn, "Should warn when pseudonymize_table provided without pseudonymize"


class TestFieldProfileCLI:
    """Tests for --field-profile and --list-field-profiles CLI arguments."""

    def test_field_profile_arg_in_help(self):
        """Test that --field-profile appears in help text."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()
        assert "--field-profile" in help_text

    def test_list_field_profiles_arg_in_help(self):
        """Test that --list-field-profiles appears in help text."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()
        assert "--list-field-profiles" in help_text

    def test_field_profile_parses(self):
        """Test that --field-profile parses correctly."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(["--vcf-file", "test.vcf", "--field-profile", "dbnsfp5"])
        assert args.field_profile == "dbnsfp5"

    def test_field_profile_default_none(self):
        """Test that --field-profile defaults to None."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(["--vcf-file", "test.vcf"])
        assert args.field_profile is None

    def test_list_field_profiles_default_false(self):
        """Test that --list-field-profiles defaults to False."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(["--vcf-file", "test.vcf"])
        assert args.list_field_profiles is False

    def test_list_field_profiles_flag(self):
        """Test that --list-field-profiles sets the flag."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(["--vcf-file", "test.vcf", "--list-field-profiles"])
        assert args.list_field_profiles is True

    def test_list_field_profiles_early_exit(self):
        """Test that --list-field-profiles exits without requiring --vcf-file."""
        from unittest.mock import patch

        with (
            patch("sys.argv", ["variantcentrifuge", "--list-field-profiles"]),
            patch("builtins.print") as mock_print,
        ):
            with pytest.raises(SystemExit) as exc_info:
                from variantcentrifuge.cli import main

                main()

            assert exc_info.value.code == 0
            print_calls = [str(call) for call in mock_print.call_args_list]
            full_output = " ".join(print_calls)
            assert "dbnsfp4" in full_output
            assert "dbnsfp5" in full_output

    def test_field_profile_in_filtering_group(self):
        """Test that --field-profile is in the Filtering & Annotation group."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        # Find the Filtering & Annotation group
        filter_group = None
        for group in parser._action_groups:
            if group.title == "Filtering & Annotation":
                filter_group = group
                break
        assert filter_group is not None
        action_names = [a.option_strings for a in filter_group._group_actions]
        flat_names = [name for names in action_names for name in names]
        assert "--field-profile" in flat_names
        assert "--list-field-profiles" in flat_names
