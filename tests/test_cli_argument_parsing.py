# File: tests/test_cli_argument_parsing.py
# Location: variantcentrifuge/tests/test_cli_argument_parsing.py

"""
Focused tests for CLI argument parsing functionality.

This module contains targeted tests for argument parsing to prevent regression
of the duplicate argument parser bug and ensure all critical parameters work.
"""

import subprocess

import pytest


def test_cli_annotate_bed_populates_canonical_and_legacy_config_keys(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    captured_configs = []

    def fake_run_pipeline(args):
        captured_configs.append(args.config.copy())

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t99\t101\tregion\n")

    monkeypatch.setattr(cli_module, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--annotate-bed",
            str(bed),
            "--output-dir",
            str(tmp_path),
        ],
    )

    exit_code = cli_module.main()

    assert exit_code == 0
    assert len(captured_configs) == 1
    cfg = captured_configs[0]
    assert cfg["annotate_bed_files"] == [str(bed)]
    assert cfg["annotate_bed"] == [str(bed)]


def test_transcript_filtering_parameters_are_documented():
    from variantcentrifuge.cli import create_parser

    parser = create_parser()
    help_text = parser.format_help()

    assert "--transcript-list" in help_text
    assert "--transcript-file" in help_text
    assert "MANE" in help_text
    assert "GRCh37" in help_text
    assert "RefSeq Select" in help_text


def test_variant_weight_column_cli_populates_config(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    captured_configs = []

    def fake_run_pipeline(args):
        captured_configs.append(args.config.copy())

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    monkeypatch.setattr(cli_module, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--output-dir",
            str(tmp_path),
            "--perform-association",
            "--association-tests",
            "logistic_burden",
            "--variant-weights",
            "score_column",
            "--variant-weight-column",
            "nephro_candidate_score",
            "--variant-weight-params",
            '{"score_min":0,"score_max":10,"floor":0.1}',
        ],
    )

    exit_code = cli_module.main()

    assert exit_code == 0
    assert captured_configs[0]["variant_weights"] == "score_column"
    assert captured_configs[0]["variant_weight_column"] == "nephro_candidate_score"
    assert captured_configs[0]["variant_weight_params"] == {
        "score_min": 0,
        "score_max": 10,
        "floor": 0.1,
    }


def test_variant_weight_params_cli_requires_json_object(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--output-dir",
            str(tmp_path),
            "--perform-association",
            "--association-tests",
            "skat",
            "--variant-weight-params",
            '["score_min",0]',
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli_module.main()

    assert exc_info.value.code == 2


def test_score_column_cli_requires_variant_weight_column(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--output-dir",
            str(tmp_path),
            "--perform-association",
            "--association-tests",
            "skat",
            "--variant-weights",
            "score_column",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli_module.main()

    assert exc_info.value.code == 2


class TestCriticalArgumentParsing:
    """Tests for critical CLI parameters that must always work."""

    def test_genotype_replacement_chunk_size_regression(self):
        """Regression test for --max-memory-gb parameter (replaced removed flags).

        Ensures critical memory configuration parameter always works.
        """
        test_values = ["8.0", "16.0", "64.0", "250.0"]

        for value in test_values:
            cmd = [
                "python",
                "-m",
                "variantcentrifuge.cli",
                "--max-memory-gb",
                value,
                "--help",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)

            # Must succeed without "unrecognized arguments" error
            assert result.returncode == 0, f"Failed with memory {value}"
            assert "unrecognized arguments" not in result.stderr, (
                f"Parameter not recognized with value {value}"
            )
            assert "--max-memory-gb" in result.stdout

    def test_performance_parameters_work(self):
        """Test that all critical performance parameters are recognized."""
        performance_params = {
            "--threads": "8",
            "--sort-memory-limit": "4G",
            "--max-memory-gb": "16",
        }

        for param, value in performance_params.items():
            cmd = ["python", "-m", "variantcentrifuge.cli", param, value, "--help"]
            result = subprocess.run(cmd, capture_output=True, text=True)

            assert result.returncode == 0, f"Parameter {param} failed"
            assert "unrecognized arguments" not in result.stderr, (
                f"Parameter {param} not recognized"
            )

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
            assert pattern not in main_source, (
                f"main() contains duplicate parser pattern: {pattern}"
            )

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
        critical_params = ["--max-memory-gb", "--threads", "--genotype-replacement-method"]

        for param in critical_params:
            # Parameter should be in help
            assert param in help_text, f"Parameter {param} missing from help"

            # Parameter should work when used (use appropriate test value per param)
            if param == "--genotype-replacement-method":
                test_value = "auto"
            elif param == "--max-memory-gb":
                test_value = "16.0"
            else:
                test_value = "10"
            cmd = ["python", "-m", "variantcentrifuge.cli", param, test_value, "--help"]
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

        # Critical parameters that must always exist
        critical_params = [
            "--max-memory-gb",
            "--threads",
            "--vcf-file",
            "--output-file",
            "--genotype-replacement-method",
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
        assert create_parser_calls == 1, (
            f"main() should call create_parser() exactly once, found {create_parser_calls}"
        )

        # main() should create minimal additional parsers for early-exit flags:
        # - status_parser for --show-checkpoint-status
        # - profile_parser for --list-field-profiles
        # - ann_parser for --show-vcf-annotations
        parser_creations = main_source.count("ArgumentParser")
        assert parser_creations <= 3, (
            f"main() should create at most 3 ArgumentParsers "
            f"(status_parser + profile_parser + ann_parser), found {parser_creations}"
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
            assert f'"{group}"' in create_parser_source, (
                f"Argument group '{group}' should be in create_parser()"
            )

            # But NOT in main() (except for status_parser which is different)
            if group != "General Options":  # status_parser might have general options
                assert f'add_argument_group("{group}")' not in main_source, (
                    f"Argument group '{group}' should NOT be duplicated in main()"
                )
