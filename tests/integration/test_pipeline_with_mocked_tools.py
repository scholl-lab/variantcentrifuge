"""Integration tests for the refactored pipeline with mocked external tools.

These tests verify the complete pipeline flow without requiring actual bioinformatics tools.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd
import pytest
from argparse import Namespace

from variantcentrifuge.pipeline_refactored import run_refactored_pipeline, build_pipeline_stages
from variantcentrifuge.pipeline_core import PipelineContext, PipelineRunner
from variantcentrifuge.pipeline_core.workspace import Workspace


class MockedToolsTestCase:
    """Base class for tests with mocked external tools."""

    @pytest.fixture(autouse=True)
    def setup_mocks(self):
        """Set up mocks for external tools."""
        self.bcftools_patcher = patch("variantcentrifuge.stages.processing_stages.run_command")
        self.snpeff_patcher = patch("variantcentrifuge.gene_bed.run_command")
        self.snpsift_patcher = patch("variantcentrifuge.filters.run_command")

        self.mock_bcftools = self.bcftools_patcher.start()
        self.mock_snpeff = self.snpeff_patcher.start()
        self.mock_snpsift = self.snpsift_patcher.start()

        # Configure default behaviors
        self._configure_tool_mocks()

        yield

        self.bcftools_patcher.stop()
        self.snpeff_patcher.stop()
        self.snpsift_patcher.stop()

    def _configure_tool_mocks(self):
        """Configure default behaviors for mocked tools."""

        # Mock bcftools responses
        def bcftools_side_effect(cmd, *args, **kwargs):
            result = Mock()
            result.returncode = 0

            if "view" in cmd:
                # Mock variant extraction
                result.stdout = ""
            elif "index" in cmd:
                # Mock index creation
                result.stdout = ""

            return result

        # Mock snpEff responses
        def snpeff_side_effect(cmd, *args, **kwargs):
            result = Mock()
            result.returncode = 0

            if "genes2bed" in " ".join(cmd):
                # Mock BED file generation
                result.stdout = "chr17\t43044294\t43125483\tBRCA1\n"

            return result

        # Mock SnpSift responses
        def snpsift_side_effect(cmd, *args, **kwargs):
            result = Mock()
            result.returncode = 0

            if "filter" in cmd:
                # Mock filtering
                result.stdout = ""
            elif "extractFields" in cmd:
                # Mock field extraction
                result.stdout = "CHROM\tPOS\tREF\tALT\nchr17\t43044295\tA\tG\n"

            return result

        self.mock_bcftools.side_effect = bcftools_side_effect
        self.mock_snpeff.side_effect = snpeff_side_effect
        self.mock_snpsift.side_effect = snpsift_side_effect


class TestBasicPipelineFlow(MockedToolsTestCase):
    """Test basic pipeline flow with single gene."""

    def test_single_gene_extraction(self, tmp_path):
        """Test complete pipeline for single gene extraction."""
        # Create test VCF file
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr17\t43044295\t.\tA\tG\t100\tPASS\t.\n"
        )

        # Create arguments
        args = Namespace(
            vcf_file=str(vcf_file),
            gene_name="BRCA1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=str(tmp_path),
            config=None,
            log_level="INFO",
            reference="GRCh37",
            preset=None,
            filter=None,
            late_filtering=False,
            final_filter=None,
            fields_to_extract=["CHROM", "POS", "REF", "ALT"],
            threads=1,
            no_stats=True,
            xlsx=False,
            html_report=False,
            phenotype_file=None,
            scoring_config_path=None,
            ped_file=None,
            calculate_inheritance=False,
            annotate_bed=None,
            annotate_gene_list=None,
            annotate_json_genes=None,
            case_samples=None,
            control_samples=None,
            case_samples_file=None,
            control_samples_file=None,
            case_phenotypes=None,
            control_phenotypes=None,
            case_phenotypes_file=None,
            control_phenotypes_file=None,
            pseudonymize=False,
            use_new_pipeline=True,
            start_time=None,
        )

        # Mock file operations
        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.stages.output_stages.pd.DataFrame.to_csv"
        ):

            # Run pipeline
            run_refactored_pipeline(args)

        # Verify stages were built correctly
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        assert "configuration_loading" in stage_names
        assert "gene_bed_creation" in stage_names
        assert "variant_extraction" in stage_names
        assert "field_extraction" in stage_names
        assert "tsv_output" in stage_names


class TestComplexPipelineFlow(MockedToolsTestCase):
    """Test complex pipeline flows with multiple features."""

    def test_pipeline_with_inheritance_and_scoring(self, tmp_path):
        """Test pipeline with inheritance analysis and scoring."""
        # Create test files
        vcf_file = tmp_path / "cohort.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tChild\tFather\tMother\n"
            "chr17\t43044295\t.\tA\tG\t100\tPASS\t.\tGT:DP:GQ\t0/1:30:99\t0/0:25:99\t0/0:28:99\n"
        )

        ped_file = tmp_path / "family.ped"
        ped_file.write_text(
            "FAM001\tFather\t0\t0\t1\t1\n"
            "FAM001\tMother\t0\t0\t2\t1\n"
            "FAM001\tChild\tFather\tMother\t1\t2\n"
        )

        # Create scoring config
        scoring_dir = tmp_path / "scoring"
        scoring_dir.mkdir()

        var_config = scoring_dir / "variable_assignment_config.json"
        var_config.write_text(json.dumps({"CADD_PHRED": {"field": "CADD_PHRED", "default": 0}}))

        formula_config = scoring_dir / "formula_config.json"
        formula_config.write_text(
            json.dumps({"formulas": [{"name": "TestScore", "formula": "CADD_PHRED * 2"}]})
        )

        # Create arguments
        args = Namespace(
            vcf_file=str(vcf_file),
            gene_name="BRCA1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=str(tmp_path),
            config=None,
            log_level="INFO",
            reference="GRCh37",
            preset=None,
            filter=None,
            late_filtering=False,
            final_filter=None,
            fields_to_extract=["CHROM", "POS", "REF", "ALT", "GT"],
            threads=1,
            no_stats=False,
            xlsx=True,
            html_report=True,
            phenotype_file=None,
            scoring_config_path=str(scoring_dir),
            ped_file=str(ped_file),
            calculate_inheritance=True,
            inheritance_mode="simple",
            annotate_bed=None,
            annotate_gene_list=None,
            annotate_json_genes=None,
            case_samples=None,
            control_samples=None,
            case_samples_file=None,
            control_samples_file=None,
            case_phenotypes=None,
            control_phenotypes=None,
            case_phenotypes_file=None,
            control_phenotypes_file=None,
            pseudonymize=False,
            use_new_pipeline=True,
            start_time=None,
            no_replacement=False,
            perform_gene_burden=False,
        )

        # Mock DataFrame operations
        mock_df = pd.DataFrame(
            {
                "CHROM": ["chr17"],
                "POS": [43044295],
                "REF": ["A"],
                "ALT": ["G"],
                "GT": ["Child(0/1);Father(0/0);Mother(0/0)"],
                "Gene_Name": ["BRCA1"],
            }
        )

        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.stages.analysis_stages.pd.read_csv", return_value=mock_df
        ), patch(
            "variantcentrifuge.stages.output_stages.pd.DataFrame.to_csv"
        ), patch(
            "variantcentrifuge.converter.convert_to_excel"
        ), patch(
            "variantcentrifuge.converter.produce_report_json"
        ), patch(
            "variantcentrifuge.generate_html_report.generate_html_report"
        ):

            # Run pipeline
            run_refactored_pipeline(args)

        # Verify analysis stages were included
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        assert "pedigree_loading" in stage_names
        assert "scoring_config_loading" in stage_names
        assert "inheritance_analysis" in stage_names
        assert "variant_scoring" in stage_names
        assert "excel_report" in stage_names
        assert "html_report" in stage_names


class TestErrorHandling(MockedToolsTestCase):
    """Test error handling in the pipeline."""

    def test_missing_gene_handling(self, tmp_path):
        """Test handling of missing gene in snpEff."""

        # Configure snpEff to return error for missing gene
        def snpeff_error_side_effect(cmd, *args, **kwargs):
            result = Mock()
            if "genes2bed" in " ".join(cmd):
                result.returncode = 1
                result.stderr = "ERROR: Gene INVALID_GENE not found in database"
            else:
                result.returncode = 0
                result.stdout = ""
            return result

        self.mock_snpeff.side_effect = snpeff_error_side_effect

        # Create test VCF
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n" "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )

        args = Namespace(
            vcf_file=str(vcf_file),
            gene_name="INVALID_GENE",
            gene_file=None,
            output_file="output.tsv",
            output_dir=str(tmp_path),
            config=None,
            log_level="INFO",
            reference="GRCh37",
            use_new_pipeline=True,
            start_time=None,
        )

        # Should raise an error
        with pytest.raises(Exception) as exc_info:
            run_refactored_pipeline(args)

        assert "Gene INVALID_GENE not found" in str(exc_info.value)

    def test_vcf_processing_error(self, tmp_path):
        """Test handling of VCF processing errors."""

        # Configure bcftools to return error
        def bcftools_error_side_effect(cmd, *args, **kwargs):
            result = Mock()
            result.returncode = 1
            result.stderr = "ERROR: Invalid VCF format"
            return result

        self.mock_bcftools.side_effect = bcftools_error_side_effect

        # Create malformed VCF
        vcf_file = tmp_path / "bad.vcf"
        vcf_file.write_text("INVALID VCF CONTENT")

        args = Namespace(
            vcf_file=str(vcf_file),
            gene_name="BRCA1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=str(tmp_path),
            config=None,
            log_level="INFO",
            reference="GRCh37",
            use_new_pipeline=True,
            start_time=None,
        )

        # Should handle error gracefully
        with pytest.raises(Exception) as exc_info:
            run_refactored_pipeline(args)

        assert "Invalid VCF format" in str(exc_info.value)


class TestParallelProcessing(MockedToolsTestCase):
    """Test parallel processing capabilities."""

    def test_parallel_variant_extraction(self, tmp_path):
        """Test parallel extraction with multiple genes."""
        # Create test VCF
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n" "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )

        # Multiple genes for parallel processing
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTP53\nEGFR\nKRAS\n")

        args = Namespace(
            vcf_file=str(vcf_file),
            gene_name=None,
            gene_file=str(gene_file),
            output_file="output.tsv",
            output_dir=str(tmp_path),
            config=None,
            log_level="INFO",
            reference="GRCh37",
            threads=4,  # Enable parallel processing
            use_new_pipeline=True,
            start_time=None,
            preset=None,
            filter=None,
            late_filtering=False,
            final_filter=None,
            fields_to_extract=["CHROM", "POS", "REF", "ALT"],
            no_stats=True,
            xlsx=False,
            html_report=False,
        )

        # Mock BED generation for multiple genes
        bed_counter = 0

        def multi_gene_snpeff(cmd, *args, **kwargs):
            nonlocal bed_counter
            result = Mock()
            result.returncode = 0

            if "genes2bed" in " ".join(cmd):
                # Return different BED regions for different genes
                genes = ["BRCA1", "TP53", "EGFR", "KRAS"]
                chroms = ["chr17", "chr17", "chr7", "chr12"]
                starts = ["43044294", "7571719", "55086724", "25357722"]
                ends = ["43125483", "7590863", "55279321", "25403870"]

                if bed_counter < len(genes):
                    result.stdout = f"{chroms[bed_counter]}\t{starts[bed_counter]}\t{ends[bed_counter]}\t{genes[bed_counter]}\n"
                    bed_counter += 1

            return result

        self.mock_snpeff.side_effect = multi_gene_snpeff

        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.stages.output_stages.pd.DataFrame.to_csv"
        ):

            # Run pipeline
            run_refactored_pipeline(args)

        # Verify parallel extraction stage was used
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        # With multiple threads, should use parallel extraction
        assert "parallel_variant_extraction" in stage_names
