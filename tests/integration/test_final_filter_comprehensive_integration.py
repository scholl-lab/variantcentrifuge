"""Comprehensive integration tests for --final-filter functionality.

This test file consolidates testing for final filtering across both legacy and
stage-based pipelines, including scoring integration and various filter scenarios.
"""

import json
from argparse import Namespace
from unittest.mock import patch

import pandas as pd
import pytest

from variantcentrifuge.cli import main
from variantcentrifuge.pipeline import build_pipeline_stages, run_refactored_pipeline


@pytest.fixture
def test_vcf_content():
    """Standard VCF content for final filter testing."""
    return """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Description="CADD Phred Score">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr17\t43044295\t.\tA\tG\t100\tPASS\tCADD_PHRED=25;AF=0.01\tGT:DP\t0/1:30
chr17\t43044296\t.\tC\tT\t90\tPASS\tCADD_PHRED=15;AF=0.05\tGT:DP\t0/1:25
chr17\t43044297\t.\tG\tA\t80\tPASS\tCADD_PHRED=35;AF=0.001\tGT:DP\t1/1:40
chr1\t1000\t.\tA\tT\t100\tPASS\tAF=0.01\tGT\t0/1
chr1\t2000\t.\tG\tC\t100\tPASS\tAF=0.5\tGT\t1/1
"""


@pytest.fixture
def scoring_config(tmp_path):
    """Create scoring configuration for testing."""
    scoring_dir = tmp_path / "scoring"
    scoring_dir.mkdir()

    # Variable assignment config
    var_config = scoring_dir / "variable_assignment_config.json"
    var_config.write_text(
        json.dumps(
            {
                "variables": {
                    "CADD_PHRED": "CADD_PHRED_value|default:0",
                    "AF": "AF_value|default:1",
                }
            }
        )
    )

    # Formula config
    formula_config = scoring_dir / "formula_config.json"
    formula_config.write_text(
        json.dumps(
            {
                "scores": {
                    "combined_score": "CADD_PHRED / 40.0 + (1 - AF) * 0.5",
                    "cadd_high": "CADD_PHRED > 20",
                }
            }
        )
    )

    return str(scoring_dir)


@pytest.fixture
def mock_external_tools():
    """Mock all external tools for pipeline testing."""

    def mock_check_tools(*args, **kwargs):
        return True

    def mock_run_command(cmd, output_file=None, **kwargs):
        """Mock external tool commands."""
        cmd_str = " ".join(cmd)

        if "genes2bed" in cmd_str:
            # Create gene BED file
            bed_file = cmd[-1]
            with open(bed_file, "w") as f:
                f.write("chr17\t43044000\t43045000\tBRCA1\n")
                f.write("chr1\t500\t3000\tTEST_GENE\n")

        elif "bcftools" in cmd_str and "view" in cmd_str:
            # Simulate variant extraction - copy input to output
            # In real scenario, would filter VCF, here we just pass through
            pass

        elif "SnpSift" in cmd_str and "filter" in cmd_str:
            # Pass through filtering (no actual filtering in mock)
            pass

        elif "SnpSift" in cmd_str and "extractFields" in cmd_str:
            # Create TSV output with test data
            output_idx = cmd.index("-o") if "-o" in cmd else cmd.index(">") + 1
            output_file = cmd[output_idx]

            with open(output_file, "w") as f:
                f.write("CHROM\tPOS\tREF\tALT\tGENE\tCADD_PHRED\tAF\tGT\n")
                f.write("chr17\t43044295\tA\tG\tBRCA1\t25\t0.01\t0/1\n")
                f.write("chr17\t43044296\tC\tT\tBRCA1\t15\t0.05\t0/1\n")
                f.write("chr17\t43044297\tG\tA\tBRCA1\t35\t0.001\t1/1\n")
                f.write("chr1\t1000\tA\tT\tTEST_GENE\t10\t0.01\t0/1\n")
                f.write("chr1\t2000\tG\tC\tTEST_GENE\t8\t0.5\t1/1\n")

    return mock_run_command, mock_check_tools


class TestFinalFilterLegacyPipeline:
    """Test final filtering in the legacy pipeline."""

    def test_final_filter_with_score(
        self, monkeypatch, tmp_path, test_vcf_content, mock_external_tools
    ):
        """Test final filter with computed score columns (legacy pipeline)."""
        mock_run_command, mock_check_tools = mock_external_tools

        # Create test files
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(test_vcf_content)

        output_dir = tmp_path / "output"
        output_file = output_dir / "test.final.tsv"

        # Apply mocks
        monkeypatch.setattr("variantcentrifuge.utils.check_external_tools", mock_check_tools)
        monkeypatch.setattr("variantcentrifuge.utils.run_command", mock_run_command)

        # Mock CLI arguments
        args = [
            "--vcf-file",
            str(vcf_file),
            "--gene-file",
            "dummy.txt",  # Will be mocked
            "--output-dir",
            str(output_dir),
            "--output-file",
            "test.tsv",
            "--preset",
            "moderate",
            "--final-filter",
            "combined_score > 0.5",
            "--add-variant-id",
        ]

        # Patch sys.argv and run
        with patch("sys.argv", ["variantcentrifuge"] + args):
            # Mock the analyze_variants function to add scores
            def mock_analyze_variants(file_handle, config):
                # Simulate scored output
                yield "CHROM\tPOS\tREF\tALT\tGENE\tcombined_score\tInheritance_Pattern"
                yield "chr17\t43044295\tA\tG\tBRCA1\t0.8\tde_novo"
                yield "chr17\t43044296\tC\tT\tBRCA1\t0.3\tinherited"
                yield "chr17\t43044297\tG\tA\tBRCA1\t0.9\tcompound_heterozygous"

            monkeypatch.setattr(
                "variantcentrifuge.analyze_variants.analyze_variants", mock_analyze_variants
            )

            # Run the main CLI
            result = main()

        # Verify results
        assert result == 0
        assert output_file.exists()

        # Check filtered results
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2  # Only scores > 0.5 should remain
        assert all(float(score) > 0.5 for score in df["combined_score"])

    def test_final_filter_with_string_columns(
        self, monkeypatch, tmp_path, test_vcf_content, mock_external_tools
    ):
        """Test final filter with string column filtering."""
        mock_run_command, mock_check_tools = mock_external_tools

        # Create test files
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(test_vcf_content)

        output_dir = tmp_path / "output"

        # Apply mocks
        monkeypatch.setattr("variantcentrifuge.utils.check_external_tools", mock_check_tools)
        monkeypatch.setattr("variantcentrifuge.utils.run_command", mock_run_command)

        # Mock CLI arguments
        args = [
            "--vcf-file",
            str(vcf_file),
            "--gene-file",
            "dummy.txt",
            "--output-dir",
            str(output_dir),
            "--output-file",
            "test.tsv",
            "--preset",
            "moderate",
            "--final-filter",
            'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]',
        ]

        with patch("sys.argv", ["variantcentrifuge"] + args):

            def mock_analyze_variants(file_handle, config):
                yield "CHROM\tPOS\tREF\tALT\tGENE\tInheritance_Pattern"
                yield "chr17\t43044295\tA\tG\tBRCA1\tde_novo"
                yield "chr17\t43044296\tC\tT\tBRCA1\tinherited"
                yield "chr17\t43044297\tG\tA\tBRCA1\tcompound_heterozygous"

            monkeypatch.setattr(
                "variantcentrifuge.analyze_variants.analyze_variants", mock_analyze_variants
            )
            result = main()

        assert result == 0
        output_file = output_dir / "test.final.tsv"
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2  # Only de_novo and compound_heterozygous
        assert set(df["Inheritance_Pattern"]) == {"de_novo", "compound_heterozygous"}

    def test_final_filter_complex_expression(
        self, monkeypatch, tmp_path, test_vcf_content, mock_external_tools
    ):
        """Test final filter with complex boolean expressions."""
        mock_run_command, mock_check_tools = mock_external_tools

        # Create test files
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(test_vcf_content)

        output_dir = tmp_path / "output"

        # Apply mocks
        monkeypatch.setattr("variantcentrifuge.utils.check_external_tools", mock_check_tools)
        monkeypatch.setattr("variantcentrifuge.utils.run_command", mock_run_command)

        # Complex filter with AND/OR logic
        args = [
            "--vcf-file",
            str(vcf_file),
            "--gene-file",
            "dummy.txt",
            "--output-dir",
            str(output_dir),
            "--output-file",
            "test.tsv",
            "--preset",
            "moderate",
            "--final-filter",
            '(score > 0.7 and GENE == "BRCA1") or '
            '(score > 0.5 and Inheritance_Pattern == "de_novo")',
        ]

        with patch("sys.argv", ["variantcentrifuge"] + args):

            def mock_analyze_variants(file_handle, config):
                yield "CHROM\tPOS\tREF\tALT\tGENE\tscore\tInheritance_Pattern"
                yield "chr17\t43044295\tA\tG\tBRCA1\t0.8\tde_novo"
                yield "chr17\t43044296\tC\tT\tBRCA1\t0.3\tinherited"
                yield "chr1\t1000\tA\tT\tTEST_GENE\t0.6\tde_novo"
                yield "chr1\t2000\tG\tC\tTEST_GENE\t0.4\tinherited"

            monkeypatch.setattr(
                "variantcentrifuge.analyze_variants.analyze_variants", mock_analyze_variants
            )
            result = main()

        assert result == 0
        output_file = output_dir / "test.final.tsv"
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2  # Should match the complex logic


class TestFinalFilterStagePipeline:
    """Test final filtering in the stage-based pipeline."""

    @patch("variantcentrifuge.utils.check_external_tools")
    @patch("variantcentrifuge.utils.run_command")
    def test_final_filter_on_computed_score(
        self, mock_run_command, mock_check_tools, tmp_path, test_vcf_content, scoring_config
    ):
        """Test that final filter can filter on computed scores in stage-based pipeline."""
        mock_check_tools.return_value = True

        # Setup mock for external tools
        def mock_command_side_effect(cmd, **kwargs):
            cmd_str = " ".join(cmd)
            if "SnpSift" in cmd_str and "extractFields" in cmd_str:
                output_file = cmd[cmd.index("-o") + 1]
                with open(output_file, "w") as f:
                    f.write("CHROM\tPOS\tREF\tALT\tGENE\tCADD_PHRED\tAF\n")
                    f.write("chr17\t43044295\tA\tG\tBRCA1\t25\t0.01\n")
                    f.write("chr17\t43044296\tC\tT\tBRCA1\t15\t0.05\n")
                    f.write("chr17\t43044297\tG\tA\tBRCA1\t35\t0.001\n")

        mock_run_command.side_effect = mock_command_side_effect

        # Create test files
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(test_vcf_content)

        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTEST_GENE\n")

        output_dir = tmp_path / "output"
        output_dir.mkdir()

        # Create args for stage-based pipeline
        args = Namespace(
            vcf_file=str(vcf_file),
            gene_file=str(gene_file),
            output_dir=str(output_dir),
            output_file="test.tsv",
            scoring_config=scoring_config,
            final_filter="combined_score > 0.6",
            preset="rare",
            
            log_level="WARNING",
            add_variant_id=True,
        )

        # Run pipeline
        result = run_refactored_pipeline(args)
        assert result == 0

        # Verify final output exists and contains filtered results
        final_output = output_dir / "test.final.tsv"
        assert final_output.exists()

    def test_stage_execution_order(self, tmp_path, test_vcf_content, scoring_config):
        """Test that FinalFilteringStage runs after VariantScoringStage."""
        # Create minimal test setup
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(test_vcf_content)

        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\n")

        args = Namespace(
            vcf_file=str(vcf_file),
            gene_file=str(gene_file),
            output_dir=str(tmp_path / "output"),
            scoring_config=scoring_config,
            final_filter="combined_score > 0.5",
            
        )

        # Build stages and check dependencies
        stages = build_pipeline_stages(args)

        # Find final filtering and scoring stages
        final_filter_stage = None
        scoring_stage = None

        for stage in stages:
            if stage.name == "final_filtering":
                final_filter_stage = stage
            elif stage.name == "variant_scoring":
                scoring_stage = stage

        assert final_filter_stage is not None
        assert scoring_stage is not None

        # Verify that final filtering depends on or runs after scoring
        # This could be through dependencies or soft_dependencies
        scoring_deps = getattr(final_filter_stage, "soft_dependencies", set())
        hard_deps = getattr(final_filter_stage, "dependencies", set())

        # Final filtering should run after scoring through soft dependencies
        assert "variant_scoring" in scoring_deps or "variant_scoring" in hard_deps
