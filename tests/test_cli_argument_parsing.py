# File: tests/test_cli_argument_parsing.py
# Location: variantcentrifuge/tests/test_cli_argument_parsing.py

"""
Focused tests for CLI argument parsing functionality.

This module contains targeted tests for argument parsing to prevent regression
of the duplicate argument parser bug and ensure all critical parameters work.
"""

import subprocess


class TestCriticalArgumentParsing:
    """Tests for critical CLI parameters that must always work."""

    def test_genotype_replacement_chunk_size_regression(self):
        """Regression test for --genotype-replacement-chunk-size parameter.

        This parameter was broken due to duplicate argument parser bug.
        This test ensures it never breaks again.
        """
        test_values = ["1000", "25000", "50000", "100000"]

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

            # Must succeed without "unrecognized arguments" error
            assert result.returncode == 0, f"Failed with chunk size {value}"
            assert (
                "unrecognized arguments" not in result.stderr
            ), f"Parameter not recognized with value {value}"
            assert "--genotype-replacement-chunk-size" in result.stdout

    def test_performance_parameters_work(self):
        """Test that all critical performance parameters are recognized."""
        performance_params = {
            "--threads": "8",
            "--chunks": "4",
            "--genotype-replacement-chunk-size": "25000",
            "--sort-memory-limit": "4G",
            "--max-memory-gb": "16",
        }

        for param, value in performance_params.items():
            cmd = ["python", "-m", "variantcentrifuge.cli", param, value, "--help"]
            result = subprocess.run(cmd, capture_output=True, text=True)

            assert result.returncode == 0, f"Parameter {param} failed"
            assert (
                "unrecognized arguments" not in result.stderr
            ), f"Parameter {param} not recognized"

    def test_genotype_replacement_methods_work(self):
        """Test that all genotype replacement methods are recognized."""
        methods = [
            "auto",
            "sequential",
            "vectorized",
            "chunked-vectorized",
            "parallel",
            "streaming-parallel",
        ]

        for method in methods:
            cmd = [
                "python",
                "-m",
                "variantcentrifuge.cli",
                "--genotype-replacement-method",
                method,
                "--help",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)

            assert result.returncode == 0, f"Method {method} failed"
            assert "unrecognized arguments" not in result.stderr, f"Method {method} not recognized"

    def test_no_duplicate_parser_indicators(self):
        """Test that ensures no signs of duplicate argument parser architecture."""
        import inspect
        import variantcentrifuge.cli as cli_module

        # Get main function source to check for architectural issues
        main_source = inspect.getsource(cli_module.main)

        # These patterns indicate duplicate parser bug is back
        forbidden_patterns = [
            'add_argument_group("General Options")',
            'add_argument_group("Core Input/Output")',
            'add_argument_group("Performance & Processing")',
            'add_argument_group("Miscellaneous Options")',
        ]

        for pattern in forbidden_patterns:
            assert (
                pattern not in main_source
            ), f"main() contains duplicate parser pattern: {pattern}"

        # main() should call create_parser()
        assert "create_parser()" in main_source, "main() should call create_parser()"

    def test_parameter_consistency_between_help_and_parsing(self):
        """Test that parameters shown in help actually work when used."""
        # Get help output
        result = subprocess.run(
            ["python", "-m", "variantcentrifuge.cli", "--help"], capture_output=True, text=True
        )

        assert result.returncode == 0
        help_text = result.stdout

        # Test that key parameters from help actually work
        critical_params = ["--genotype-replacement-chunk-size", "--threads", "--chunks"]

        for param in critical_params:
            # Parameter should be in help
            assert param in help_text, f"Parameter {param} missing from help"

            # Parameter should work when used
            cmd = ["python", "-m", "variantcentrifuge.cli", param, "10", "--help"]
            test_result = subprocess.run(cmd, capture_output=True, text=True)
            assert test_result.returncode == 0, f"Parameter {param} doesn't work"
            assert "unrecognized arguments" not in test_result.stderr


class TestArgumentParserArchitecture:
    """Tests for the overall argument parser architecture."""

    def test_create_parser_contains_all_critical_parameters(self):
        """Test that create_parser() contains all critical parameters."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Critical parameters that caused the original bug
        critical_params = [
            "--genotype-replacement-chunk-size",
            "--threads",
            "--vcf-file",
            "--output-file",
            "--chunks",
        ]

        for param in critical_params:
            assert param in help_text, f"Critical parameter {param} missing from create_parser()"

    def test_main_function_architecture(self):
        """Test that main() function has correct architecture."""
        import inspect
        import variantcentrifuge.cli as cli_module

        # Get main function source
        main_source = inspect.getsource(cli_module.main)

        # main() should call create_parser() exactly once
        create_parser_calls = main_source.count("create_parser()")
        assert (
            create_parser_calls == 1
        ), f"main() should call create_parser() exactly once, found {create_parser_calls}"

        # main() should create minimal additional parsers (only status_parser)
        parser_creations = main_source.count("ArgumentParser")
        assert parser_creations <= 1, (
            f"main() should create at most 1 ArgumentParser (status_parser), "
            f"found {parser_creations}"
        )

    def test_argument_groups_only_in_create_parser(self):
        """Test that argument groups are only defined in create_parser()."""
        import inspect
        import variantcentrifuge.cli as cli_module

        # Get both function sources
        create_parser_source = inspect.getsource(cli_module.create_parser)
        main_source = inspect.getsource(cli_module.main)

        # Argument groups should exist in create_parser
        expected_groups = ["General Options", "Performance & Processing", "Miscellaneous Options"]

        for group in expected_groups:
            assert (
                f'"{group}"' in create_parser_source
            ), f"Argument group '{group}' should be in create_parser()"

            # But NOT in main() (except for status_parser which is different)
            if group != "General Options":  # status_parser might have general options
                assert (
                    f'add_argument_group("{group}")' not in main_source
                ), f"Argument group '{group}' should NOT be duplicated in main()"
