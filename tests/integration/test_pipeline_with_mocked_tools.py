"""Integration tests for the refactored pipeline with mocked external tools.

These tests verify the complete pipeline flow without requiring actual bioinformatics tools.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
import pandas as pd
import pytest
from argparse import Namespace

from variantcentrifuge.pipeline_refactored import run_refactored_pipeline, build_pipeline_stages


class MockedToolsTestCase:
    """Base class for tests with mocked external tools."""

    @pytest.fixture(autouse=True)
    def setup_mocks(self):
        """Set up mocks for external tools."""
        self.bcftools_patcher = patch("variantcentrifuge.stages.processing_stages.run_command")
        self.snpeff_patcher = patch("variantcentrifuge.gene_bed.subprocess.run")
        self.snpsift_patcher = patch("variantcentrifuge.filters.run_command")
        # Also patch run_command in extractor module
        self.extractor_patcher = patch("variantcentrifuge.extractor.run_command")

        self.mock_bcftools = self.bcftools_patcher.start()
        self.mock_snpeff = self.snpeff_patcher.start()
        self.mock_snpsift = self.snpsift_patcher.start()
        self.mock_extractor = self.extractor_patcher.start()

        # Configure default behaviors
        self._configure_tool_mocks()

        yield

        self.bcftools_patcher.stop()
        self.snpeff_patcher.stop()
        self.snpsift_patcher.stop()
        self.extractor_patcher.stop()

    def _configure_tool_mocks(self):
        """Configure default behaviors for mocked tools."""
        # Mock utils run_command for SnpSift extractFields
        def utils_side_effect(cmd, output_file=None, *args, **kwargs):
            if "SnpSift" in cmd and "extractFields" in cmd:
                # Extract field names from the command (last N args after the filename)
                field_idx = None
                for i, arg in enumerate(cmd):
                    if arg.endswith('.vcf') or arg.endswith('.vcf.gz'):
                        field_idx = i + 1
                        break
                
                if field_idx and field_idx < len(cmd):
                    fields = cmd[field_idx:]
                    header = "\t".join(fields) + "\n"
                    # Generate dummy data based on field count
                    data_values = []
                    for field in fields:
                        if field == "CHROM":
                            data_values.append("chr17")
                        elif field == "POS":
                            data_values.append("43044295")
                        elif field == "REF":
                            data_values.append("A")
                        elif field == "ALT":
                            data_values.append("G")
                        elif field == "QUAL":
                            data_values.append("30")
                        else:
                            data_values.append("NA")
                    data = "\t".join(data_values) + "\n"
                    content = header + data
                else:
                    content = "CHROM\tPOS\tREF\tALT\nchr17\t43044295\tA\tG\n"
                
                if output_file:
                    with open(output_file, 'w') as f:
                        f.write(content)
                    return output_file
                else:
                    return content
            # Default behavior
            result = Mock()
            result.returncode = 0
            result.stdout = ""
            return result

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
            # subprocess.run returns a CompletedProcess object
            from subprocess import CompletedProcess

            # Check if this is a snpEff genes2bed command
            if isinstance(cmd, list) and "genes2bed" in cmd:
                # Write BED content to the output file if stdout is redirected
                if "stdout" in kwargs and hasattr(kwargs["stdout"], "write"):
                    kwargs["stdout"].write("chr17\t43044294\t43125483\tBRCA1\n")
                return CompletedProcess(cmd, 0)

            # For sortBed commands
            if isinstance(cmd, list) and "sortBed" in cmd:
                # Just return success
                if "stdout" in kwargs and hasattr(kwargs["stdout"], "write"):
                    kwargs["stdout"].write("chr17\t43044294\t43125483\tBRCA1\n")
                return CompletedProcess(cmd, 0)

            return CompletedProcess(cmd if isinstance(cmd, list) else [cmd], 0)

        # Mock SnpSift responses
        def snpsift_side_effect(cmd, *args, **kwargs):
            result = Mock()
            result.returncode = 0

            if "filter" in cmd:
                # Mock filtering
                result.stdout = ""
            elif "extractFields" in cmd:
                # Mock field extraction - write to output file if specified
                if "-o" in cmd:
                    output_idx = cmd.index("-o") + 1
                    if output_idx < len(cmd):
                        output_file = cmd[output_idx]
                        # Write TSV output
                        with open(output_file, 'w') as f:
                            f.write("CHROM\tPOS\tREF\tALT\tQUAL\nchr17\t43044295\tA\tG\t30\n")
                result.stdout = "CHROM\tPOS\tREF\tALT\tQUAL\nchr17\t43044295\tA\tG\t30\n"

            return result

        self.mock_bcftools.side_effect = bcftools_side_effect
        self.mock_snpeff.side_effect = snpeff_side_effect
        self.mock_snpsift.side_effect = snpsift_side_effect
        self.mock_extractor.side_effect = utils_side_effect


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
            fields_to_extract="CHROM POS REF ALT",  # Should be a string, not a list
            extract=["CHROM", "POS", "REF", "ALT"],  # Also need extract field for config
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
            keep_intermediates=False,
            enable_checkpoint=False,
            no_replacement=False,
        )

        # Mock analyze_variants to just return input data
        def mock_analyze_variants(inp, config):
            # Read lines from input file
            lines = list(inp)
            # Add Gene_Name column if not present
            if lines:
                header = lines[0].strip()
                if "Gene_Name" not in header:
                    lines[0] = header + "\tGene_Name\n"
                    # Add BRCA1 to all data lines
                    for i in range(1, len(lines)):
                        lines[i] = lines[i].strip() + "\tBRCA1\n"
            # Return as generator
            for line in lines:
                yield line.strip()
        
        # Mock file operations
        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.stages.analysis_stages.analyze_variants", side_effect=mock_analyze_variants
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
            fields_to_extract="CHROM POS REF ALT GT",  # Should be a string, not a list
            extract=["CHROM", "POS", "REF", "ALT", "GT"],  # Also need extract field
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
            keep_intermediates=False,
            enable_checkpoint=False,
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

        # Mock read_pedigree to return expected data
        mock_pedigree_data = {
            "Father": {
                "family_id": "FAM001",
                "sample_id": "Father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1"
            },
            "Mother": {
                "family_id": "FAM001",
                "sample_id": "Mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1"
            },
            "Child": {
                "family_id": "FAM001",
                "sample_id": "Child",
                "father_id": "Father",
                "mother_id": "Mother",
                "sex": "1",
                "affected_status": "2"
            }
        }
        
        # Mock analyze_variants to just return input data
        def mock_analyze_variants(inp, config):
            # Read lines from input file
            lines = list(inp)
            # Add Gene_Name column if not present
            if lines:
                header = lines[0].strip()
                if "Gene_Name" not in header:
                    lines[0] = header + "\tGene_Name\n"
                    # Add BRCA1 to all data lines
                    for i in range(1, len(lines)):
                        lines[i] = lines[i].strip() + "\tBRCA1\n"
            # Return as generator
            for line in lines:
                yield line.strip()
        
        # Mock functions need to accept **kwargs for extra arguments
        def mock_convert_to_excel(*args, **kwargs):
            pass
        
        def mock_produce_report_json(*args, **kwargs):
            return {}
        
        def mock_generate_html_report(*args, **kwargs):
            pass
        
        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.stages.analysis_stages.pd.read_csv", return_value=mock_df
        ), patch(
            "variantcentrifuge.converter.convert_to_excel", side_effect=mock_convert_to_excel
        ), patch(
            "variantcentrifuge.converter.produce_report_json", side_effect=mock_produce_report_json
        ), patch(
            "variantcentrifuge.generate_html_report.generate_html_report", side_effect=mock_generate_html_report
        ), patch(
            "variantcentrifuge.stages.setup_stages.read_pedigree", return_value=mock_pedigree_data
        ), patch(
            "variantcentrifuge.stages.analysis_stages.analyze_variants", side_effect=mock_analyze_variants
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
        # When multiple reports are requested, a combined stage is used
        assert "parallel_report_generation" in stage_names


class TestErrorHandling(MockedToolsTestCase):
    """Test error handling in the pipeline."""

    def test_missing_gene_handling(self, tmp_path):
        """Test handling of missing gene in snpEff."""
        # Configure snpEff to return error for missing gene
        def snpeff_error_side_effect(cmd, *args, **kwargs):
            from subprocess import CompletedProcess, CalledProcessError
            
            if isinstance(cmd, list) and "genes2bed" in cmd and "INVALID_GENE" in cmd:
                # Raise CalledProcessError when check=True and returncode != 0
                raise CalledProcessError(
                    1, cmd, stderr="ERROR: Gene INVALID_GENE not found in database"
                )
            # For other commands, return success
            return CompletedProcess(cmd if isinstance(cmd, list) else [cmd], 0)

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
            phenotype_file=None,
            fields_to_extract="CHROM POS REF ALT",
            threads=1,
            preset=None,
            filter=None,
            late_filtering=False,
            final_filter=None,
            no_stats=True,
            xlsx=False,
            html_report=False,
            keep_intermediates=False,
            enable_checkpoint=False,
            no_replacement=False,
        )

        # Should raise an error
        with pytest.raises(Exception) as exc_info:
            run_refactored_pipeline(args)

        # Check for either the original error message or the stage execution error
        error_msg = str(exc_info.value)
        assert ("INVALID_GENE" in error_msg and "genes2bed" in error_msg) or "Gene INVALID_GENE not found" in error_msg

    def test_vcf_processing_error(self, tmp_path):
        """Test handling of VCF processing errors."""
        # Configure bcftools to return error
        def bcftools_error_side_effect(cmd, *args, **kwargs):
            from variantcentrifuge.utils import CommandError
            # Raise CommandError which is what run_command raises on failure
            raise CommandError("ERROR: Invalid VCF format")

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
            phenotype_file=None,
            fields_to_extract="CHROM POS REF ALT",
            threads=1,
            preset=None,
            filter=None,
            late_filtering=False,
            final_filter=None,
            no_stats=True,
            xlsx=False,
            html_report=False,
            keep_intermediates=False,
            enable_checkpoint=False,
            no_replacement=False,
        )

        # Should handle error gracefully
        with pytest.raises(Exception) as exc_info:
            run_refactored_pipeline(args)

        # Check for error related to VCF processing
        error_msg = str(exc_info.value)
        assert "Invalid VCF format" in error_msg or "Variant extraction failed" in error_msg


class TestParallelProcessing(MockedToolsTestCase):
    """Test parallel processing capabilities."""

    @pytest.mark.xfail(reason="Parallel processing with temporary directories requires complex mocking")
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
            fields_to_extract="CHROM POS REF ALT",  # Should be a string, not a list
            extract=["CHROM", "POS", "REF", "ALT"],  # Also need extract field
            no_stats=True,
            xlsx=False,
            html_report=False,
            phenotype_file=None,
            keep_intermediates=False,
            enable_checkpoint=False,
            no_replacement=False,
        )

        # Mock BED generation for multiple genes
        bed_counter = 0

        def multi_gene_snpeff(cmd, *args, **kwargs):
            nonlocal bed_counter
            from subprocess import CompletedProcess
            
            if isinstance(cmd, list) and "genes2bed" in cmd:
                # Return different BED regions for different genes
                genes = ["BRCA1", "TP53", "EGFR", "KRAS"]
                chroms = ["chr17", "chr17", "chr7", "chr12"]
                starts = ["43044294", "7571719", "55086724", "25357722"]
                ends = ["43125483", "7590863", "55279321", "25403870"]
                
                # Find which gene is being requested
                gene_idx = -1
                for i, gene in enumerate(genes):
                    if gene in cmd:
                        gene_idx = i
                        break
                
                if gene_idx >= 0 and "stdout" in kwargs and hasattr(kwargs["stdout"], "write"):
                    kwargs["stdout"].write(
                        f"{chroms[gene_idx]}\t{starts[gene_idx]}\t"
                        f"{ends[gene_idx]}\t{genes[gene_idx]}\n"
                    )
                
                return CompletedProcess(cmd, 0)
            
            # For sortBed
            if isinstance(cmd, list) and "sortBed" in cmd:
                if "stdout" in kwargs and hasattr(kwargs["stdout"], "write"):
                    # Just echo back some BED content
                    kwargs["stdout"].write("chr17\t43044294\t43125483\tBRCA1\n")
                return CompletedProcess(cmd, 0)
                
            return CompletedProcess(cmd if isinstance(cmd, list) else [cmd], 0)

        self.mock_snpeff.side_effect = multi_gene_snpeff

        # Mock analyze_variants to just return input data
        def mock_analyze_variants(inp, config):
            # Read lines from input file
            lines = list(inp)
            # Add Gene_Name column if not present
            if lines:
                header = lines[0].strip()
                if "Gene_Name" not in header:
                    lines[0] = header + "\tGene_Name\n"
                    # Add gene names to all data lines
                    for i in range(1, len(lines)):
                        lines[i] = lines[i].strip() + "\tBRCA1\n"
            # Return as generator
            for line in lines:
                yield line.strip()

        # Mock chunk file creation
        original_bcftools = self.mock_bcftools.side_effect
        
        def bcftools_chunk_side_effect(cmd, *args, **kwargs):
            # If this is creating a chunk file, create it
            if isinstance(cmd, list) and "view" in cmd and "-o" in cmd:
                output_idx = cmd.index("-o") + 1
                if output_idx < len(cmd):
                    output_file = cmd[output_idx]
                    # Write valid VCF content - need to handle .gz files
                    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
                    if output_file.endswith('.gz'):
                        import gzip
                        content = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        with gzip.open(output_file, 'wb') as f:
                            f.write(content)
                    else:
                        Path(output_file).write_text(
                            "##fileformat=VCFv4.2\n"
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        )
            
            # Call original if it exists
            if original_bcftools:
                return original_bcftools(cmd, *args, **kwargs)
            else:
                result = Mock()
                result.returncode = 0
                result.stdout = ""
                return result
        
        self.mock_bcftools.side_effect = bcftools_chunk_side_effect
        
        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.stages.analysis_stages.analyze_variants", side_effect=mock_analyze_variants
        ):

            # Run pipeline
            run_refactored_pipeline(args)

        # Verify parallel extraction stage was used
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        # With multiple threads, should use parallel extraction
        assert "parallel_variant_extraction" in stage_names


class TestBCFToolsPrefilter(MockedToolsTestCase):
    """Test bcftools prefilter functionality in the pipeline."""

    @pytest.fixture
    def temp_files(self):
        """Create temporary files for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            vcf_file = tmp_path / "test.vcf"
            vcf_file.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr17\t43044295\t.\tA\tG\t30\tPASS\tAC=2\n"
            )
            yield {"vcf": str(vcf_file), "output": str(tmp_path / "output.tsv"), "tmpdir": tmp_path}

    def test_bcftools_prefilter_applied(self, temp_files):
        """Test that bcftools prefilter is correctly applied during extraction."""
        args = Namespace(
            gene_name="BRCA1",
            gene_file=None,
            vcf_file=temp_files["vcf"],
            output_file=temp_files["output"],
            output_dir=str(temp_files["tmpdir"]),
            bcftools_prefilter='FILTER="PASS" & INFO/AC[0] < 10',
            preset=None,
            filter=None,
            filters=None,
            late_filtering=False,
            final_filter=None,
            fields_to_extract="CHROM POS REF ALT QUAL",
            extract=["CHROM", "POS", "REF", "ALT", "QUAL"],
            keep_intermediates=False,
            enable_checkpoint=False,
            no_replacement=False,
            threads=1,
            quiet=True,
            log_level="INFO",
            reference="GRCh37",
            config=None,
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

        # Track bcftools calls
        bcftools_calls = []

        # Store the original side effect from the base class
        original_bcftools_side_effect = self.mock_bcftools.side_effect
        
        def track_bcftools_calls(cmd, *args, **kwargs):
            # Track the call
            if isinstance(cmd, list) and len(cmd) > 0 and cmd[0] == "bcftools":
                bcftools_calls.append(cmd)
            
            # Create expected output files
            if isinstance(cmd, list) and "view" in cmd and "-o" in cmd:
                output_idx = cmd.index("-o") + 1
                output_file = cmd[output_idx]
                # Write valid VCF content
                Path(output_file).parent.mkdir(parents=True, exist_ok=True)
                Path(output_file).write_text(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr17\t43044295\t.\tA\tG\t30\tPASS\tAC=2\n"
                )
            
            # Call the original side effect for default behavior
            if original_bcftools_side_effect:
                return original_bcftools_side_effect(cmd, *args, **kwargs)
            else:
                result = Mock()
                result.returncode = 0
                result.stdout = ""
                return result

        # bcftools commands in filters.py are caught by snpsift_patcher
        self.mock_snpsift.side_effect = track_bcftools_calls

        # Mock analyze_variants to just return input data
        def mock_analyze_variants(inp, config):
            # Read lines from input file
            lines = list(inp)
            # Add Gene_Name column if not present
            if lines:
                header = lines[0].strip()
                if "Gene_Name" not in header:
                    lines[0] = header + "\tGene_Name\n"
                    # Add BRCA1 to all data lines
                    for i in range(1, len(lines)):
                        lines[i] = lines[i].strip() + "\tBRCA1\n"
            # Return as generator
            for line in lines:
                yield line.strip()

        # Set up other mocks
        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.helpers.get_vcf_samples", return_value=set()
        ), patch(
            "variantcentrifuge.stages.analysis_stages.analyze_variants", side_effect=mock_analyze_variants
        ):
            # Run pipeline
            run_refactored_pipeline(args)

        # Verify bcftools was called with the prefilter
        view_calls = [call for call in bcftools_calls if "view" in call]
        assert len(view_calls) > 0, "bcftools view should have been called"

        # Check that at least one view call has the prefilter
        prefilter_found = False
        for call in view_calls:
            if "-i" in call:
                idx = call.index("-i")
                if idx + 1 < len(call) and 'FILTER="PASS"' in call[idx + 1]:
                    prefilter_found = True
                    break

        assert prefilter_found, "bcftools prefilter expression should have been applied"

    def test_bcftools_prefilter_with_complex_expression(self, temp_files):
        """Test bcftools prefilter with complex expression."""
        args = Namespace(
            gene_name="BRCA1",
            gene_file=None,
            vcf_file=temp_files["vcf"],
            output_file=temp_files["output"],
            output_dir=str(temp_files["tmpdir"]),
            bcftools_prefilter="(QUAL >= 30) & (FORMAT/DP[*] >= 10) & (INFO/AC[0] < 5)",
            preset=None,
            filter=None,
            filters=None,
            late_filtering=False,
            final_filter=None,
            fields_to_extract="CHROM POS REF ALT QUAL",
            extract=["CHROM", "POS", "REF", "ALT", "QUAL"],
            keep_intermediates=False,
            enable_checkpoint=False,
            no_replacement=False,
            threads=1,
            quiet=True,
            log_level="INFO",
            reference="GRCh37",
            config=None,
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

        # Track bcftools calls
        bcftools_calls = []

        # Store the original side effect from the base class
        original_bcftools_side_effect = self.mock_bcftools.side_effect
        
        def track_bcftools_calls(cmd, *args, **kwargs):
            # Track the call
            if isinstance(cmd, list) and len(cmd) > 0 and cmd[0] == "bcftools":
                bcftools_calls.append(cmd)
            
            # Create expected output files
            if isinstance(cmd, list) and "view" in cmd and "-o" in cmd:
                output_idx = cmd.index("-o") + 1
                output_file = cmd[output_idx]
                # Write valid VCF content
                Path(output_file).parent.mkdir(parents=True, exist_ok=True)
                Path(output_file).write_text(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr17\t43044295\t.\tA\tG\t30\tPASS\tAC=2\n"
                )
            
            # Call the original side effect for default behavior
            if original_bcftools_side_effect:
                return original_bcftools_side_effect(cmd, *args, **kwargs)
            else:
                result = Mock()
                result.returncode = 0
                result.stdout = ""
                return result

        # bcftools commands in filters.py are caught by snpsift_patcher
        self.mock_snpsift.side_effect = track_bcftools_calls

        # Mock analyze_variants to just return input data
        def mock_analyze_variants(inp, config):
            # Read lines from input file
            lines = list(inp)
            # Add Gene_Name column if not present
            if lines:
                header = lines[0].strip()
                if "Gene_Name" not in header:
                    lines[0] = header + "\tGene_Name\n"
                    # Add BRCA1 to all data lines
                    for i in range(1, len(lines)):
                        lines[i] = lines[i].strip() + "\tBRCA1\n"
            # Return as generator
            for line in lines:
                yield line.strip()

        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.helpers.get_vcf_samples", return_value=set()
        ), patch(
            "variantcentrifuge.stages.analysis_stages.analyze_variants", side_effect=mock_analyze_variants
        ):
            # Run pipeline
            run_refactored_pipeline(args)

        # Verify the complex expression was passed correctly
        view_calls = [call for call in bcftools_calls if "view" in call]
        complex_filter_found = False

        for call in view_calls:
            if "-i" in call:
                idx = call.index("-i")
                if idx + 1 < len(call) and "(QUAL >= 30)" in call[idx + 1]:
                    complex_filter_found = True
                    assert "FORMAT/DP[*] >= 10" in call[idx + 1]
                    assert "INFO/AC[0] < 5" in call[idx + 1]
                    break

        assert complex_filter_found, "Complex bcftools filter expression should have been applied"
