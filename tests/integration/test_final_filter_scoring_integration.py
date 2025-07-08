"""Integration test for final filtering with scoring.

This test verifies that the FinalFilteringStage correctly runs after
scoring stages and can filter on computed score columns.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
import pandas as pd
import pytest
from argparse import Namespace

from variantcentrifuge.pipeline_refactored import run_refactored_pipeline, build_pipeline_stages


class TestFinalFilterWithScoring:
    """Test final filtering with scoring in the refactored pipeline."""

    @pytest.fixture
    def setup_files(self):
        """Create test files for the pipeline."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            
            # Create test VCF file
            vcf_file = tmp_path / "test.vcf"
            vcf_file.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                "chr17\t43044295\t.\tA\tG\t100\tPASS\tCADD_PHRED=25;AF=0.01\tGT:DP\t0/1:30\n"
                "chr17\t43044296\t.\tC\tT\t90\tPASS\tCADD_PHRED=15;AF=0.05\tGT:DP\t0/1:25\n"
                "chr17\t43044297\t.\tG\tA\t80\tPASS\tCADD_PHRED=35;AF=0.001\tGT:DP\t1/1:40\n"
            )

            # Create scoring config
            scoring_dir = tmp_path / "scoring"
            scoring_dir.mkdir()

            var_config = scoring_dir / "variable_assignment_config.json"
            var_config.write_text(json.dumps({
                "CADD_PHRED": {"field": "CADD_PHRED", "default": 0},
                "AF": {"field": "AF", "default": 1}
            }))

            formula_config = scoring_dir / "formula_config.json"
            formula_config.write_text(json.dumps({
                "formulas": [
                    {
                        "name": "base_score",
                        "formula": "CADD_PHRED"
                    },
                    {
                        "name": "nephro_candidate_score",
                        "formula": "CADD_PHRED * (1 if AF < 0.01 else 0.5)"
                    }
                ]
            }))

            # Create output directory and intermediate directory
            output_dir = tmp_path / "output"
            output_dir.mkdir(exist_ok=True)
            (output_dir / "intermediate").mkdir(exist_ok=True)
            
            yield {
                "vcf": str(vcf_file),
                "scoring_dir": str(scoring_dir),
                "output_dir": str(output_dir)
            }

    @patch("variantcentrifuge.utils.run_command")
    @patch("variantcentrifuge.gene_bed.subprocess.run")
    @patch("variantcentrifuge.helpers.get_vcf_samples")
    def test_final_filter_on_computed_score(self, mock_get_samples, mock_snpeff, mock_run_command, setup_files):
        """Test that final filter can filter on computed scores."""
        # Mock VCF samples
        mock_get_samples.return_value = ["Sample1"]

        # Mock snpEff for gene BED creation
        def snpeff_side_effect(cmd, *args, **kwargs):
            from subprocess import CompletedProcess
            if isinstance(cmd, list) and "genes2bed" in cmd:
                if "stdout" in kwargs and hasattr(kwargs["stdout"], "write"):
                    kwargs["stdout"].write("chr17\t43044294\t43125483\tBRCA1\n")
            return CompletedProcess(cmd, 0)
        mock_snpeff.side_effect = snpeff_side_effect

        # Mock run_command for bcftools and SnpSift
        def run_command_side_effect(cmd, output_file=None, *args, **kwargs):
            result = Mock()
            result.returncode = 0
            result.stdout = ""

            if "bcftools" in cmd and "view" in cmd:
                # Mock variant extraction
                if output_file:
                    Path(output_file).write_text(
                        "##fileformat=VCFv4.2\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                        "chr17\t43044295\t.\tA\tG\t100\tPASS\tCADD_PHRED=25;AF=0.01\tGT:DP\t0/1:30\n"
                        "chr17\t43044296\t.\tC\tT\t90\tPASS\tCADD_PHRED=15;AF=0.05\tGT:DP\t0/1:25\n"
                        "chr17\t43044297\t.\tG\tA\t80\tPASS\tCADD_PHRED=35;AF=0.001\tGT:DP\t1/1:40\n"
                    )
            elif "SnpSift" in cmd and "extractFields" in cmd:
                # Mock field extraction
                if output_file:
                    with open(output_file, 'w') as f:
                        f.write("CHROM\tPOS\tREF\tALT\tQUAL\tCADD_PHRED\tAF\tGT\n")
                        f.write("chr17\t43044295\tA\tG\t100\t25\t0.01\tSample1(0/1)\n")
                        f.write("chr17\t43044296\tC\tT\t90\t15\t0.05\tSample1(0/1)\n")
                        f.write("chr17\t43044297\tG\tA\t80\t35\t0.001\tSample1(1/1)\n")

            return result

        mock_run_command.side_effect = run_command_side_effect

        # Create arguments with final filter on computed score
        args = Namespace(
            vcf_file=setup_files["vcf"],
            gene_name="BRCA1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=setup_files["output_dir"],
            config=None,
            log_level="DEBUG",  # Enable debug logging
            reference="GRCh37",
            preset=None,
            filter=None,
            late_filtering=False,
            final_filter="nephro_candidate_score > 20",  # Filter on computed score
            fields_to_extract="CHROM POS REF ALT QUAL CADD_PHRED AF GT",
            extract=["CHROM", "POS", "REF", "ALT", "QUAL", "CADD_PHRED", "AF", "GT"],
            threads=1,
            no_stats=True,
            xlsx=False,
            html_report=False,
            phenotype_file=None,
            scoring_config_path=setup_files["scoring_dir"],
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
            perform_gene_burden=False,
            skip_variant_analysis=False,
        )

        # Mock analyze_variants to add Gene_Name
        def mock_analyze_variants(inp, config):
            lines = list(inp)
            if lines:
                header = lines[0].strip()
                if "Gene_Name" not in header:
                    lines[0] = header + "\tGene_Name\n"
                    for i in range(1, len(lines)):
                        lines[i] = lines[i].strip() + "\tBRCA1\n"
            for line in lines:
                yield line.strip()

        # Track DataFrame writes to verify filtering
        written_dfs = []
        def track_df_write(df, *args, **kwargs):
            written_dfs.append(df.copy())

        with patch("variantcentrifuge.stages.processing_stages.Path.exists", return_value=True), \
             patch("variantcentrifuge.stages.processing_stages.Path.touch"), \
             patch("variantcentrifuge.stages.output_stages.pd.DataFrame.to_csv", side_effect=track_df_write), \
             patch("variantcentrifuge.stages.analysis_stages.analyze_variants", side_effect=mock_analyze_variants), \
             patch("variantcentrifuge.pipeline_core.workspace.Path.mkdir"):

            # Run pipeline
            run_refactored_pipeline(args)

        # Verify stages were built correctly
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        # Verify all necessary stages are present
        assert "scoring_config_loading" in stage_names
        assert "variant_scoring" in stage_names
        assert "final_filtering" in stage_names
        assert "tsv_output" in stage_names

        # Verify filtering happened
        assert len(written_dfs) > 0, "No DataFrames were written"
        final_df = written_dfs[-1]  # Last written DataFrame should be the final output

        # Check that scoring columns exist
        assert "base_score" in final_df.columns, f"Columns: {list(final_df.columns)}"
        assert "nephro_candidate_score" in final_df.columns, f"Columns: {list(final_df.columns)}"

        # Verify filtering worked:
        # - chr17:43044295: CADD=25, AF=0.01 -> score = 25 * 1 = 25 (PASS)
        # - chr17:43044296: CADD=15, AF=0.05 -> score = 15 * 0.5 = 7.5 (FAIL)
        # - chr17:43044297: CADD=35, AF=0.001 -> score = 35 * 1 = 35 (PASS)
        assert len(final_df) == 2, f"Expected 2 variants after filtering, got {len(final_df)}"
        
        # Check the remaining variants have high scores
        assert all(final_df["nephro_candidate_score"] > 20), \
            f"Some variants have low scores: {final_df['nephro_candidate_score'].tolist()}"

    def test_stage_execution_order(self, setup_files):
        """Test that FinalFilteringStage runs after VariantScoringStage."""
        args = Namespace(
            vcf_file=setup_files["vcf"],
            gene_name="BRCA1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=setup_files["output_dir"],
            config=None,
            log_level="INFO",
            reference="GRCh37",
            scoring_config_path=setup_files["scoring_dir"],
            final_filter="nephro_candidate_score > 20",
            fields_to_extract="CHROM POS REF ALT",
            extract=["CHROM", "POS", "REF", "ALT"],
            no_stats=True,
            xlsx=False,
            html_report=False,
            use_new_pipeline=True,
            phenotype_file=None,
            ped_file=None,
            annotate_bed=None,
            annotate_gene_list=None,
            annotate_json_genes=None,
            threads=1,
        )

        # Build stages
        stages = build_pipeline_stages(args)
        
        # Find the stages we care about
        scoring_stage = None
        filtering_stage = None
        for stage in stages:
            if stage.name == "variant_scoring":
                scoring_stage = stage
            elif stage.name == "final_filtering":
                filtering_stage = stage

        assert scoring_stage is not None, "VariantScoringStage not found"
        assert filtering_stage is not None, "FinalFilteringStage not found"

        # Check dependencies
        assert "variant_scoring" in filtering_stage.soft_dependencies, \
            "FinalFilteringStage should have variant_scoring as soft dependency"

        # Create a runner to check execution order
        from variantcentrifuge.pipeline_core import PipelineRunner
        runner = PipelineRunner()
        
        # Get execution plan
        plan = runner.dry_run(stages)
        
        # Find which level each stage is at
        scoring_level = None
        filtering_level = None
        for level_idx, level in enumerate(plan):
            if "variant_scoring" in level:
                scoring_level = level_idx
            if "final_filtering" in level:
                filtering_level = level_idx

        assert scoring_level is not None, "variant_scoring not in execution plan"
        assert filtering_level is not None, "final_filtering not in execution plan"
        assert filtering_level > scoring_level, \
            f"final_filtering (level {filtering_level}) should run after variant_scoring (level {scoring_level})"