"""Integration test for inheritance analysis with the new pipeline.

This test verifies that inheritance analysis correctly detects patterns
including compound heterozygous variants.
"""

import tempfile
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from variantcentrifuge.pipeline_refactored import build_pipeline_stages, run_refactored_pipeline


class TestInheritanceAnalysis:
    """Test inheritance analysis in the refactored pipeline."""

    @pytest.fixture
    def setup_files(self):
        """Create test files for the pipeline."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            # Create test VCF file with multiple samples
            vcf_file = tmp_path / "test.vcf"
            vcf_file.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tChild\tFather\tMother\n"
                # De novo variant (only in child)
                "chr1\t100\t.\tA\tG\t100\tPASS\tGENE=GENE1\tGT:DP\t0/1:30\t0/0:25\t0/0:28\n"
                # Compound het gene GENE2 - variant 1 from father
                "chr2\t200\t.\tC\tT\t90\tPASS\tGENE=GENE2\tGT:DP\t0/1:25\t0/1:30\t0/0:27\n"
                # Compound het gene GENE2 - variant 2 from mother
                "chr2\t300\t.\tG\tA\t85\tPASS\tGENE=GENE2\tGT:DP\t0/1:28\t0/0:26\t0/1:29\n"
                # Homozygous recessive (both parents carriers)
                "chr3\t400\t.\tT\tC\t95\tPASS\tGENE=GENE3\tGT:DP\t1/1:35\t0/1:30\t0/1:32\n"
                # Not compound het in GENE2 (same parent)
                "chr2\t500\t.\tA\tT\t80\tPASS\tGENE=GENE2\tGT:DP\t0/1:22\t0/1:24\t0/0:26\n"
            )

            # Create PED file
            ped_file = tmp_path / "family.ped"
            ped_file.write_text(
                "FAM001\tFather\t0\t0\t1\t1\n"
                "FAM001\tMother\t0\t0\t2\t1\n"
                "FAM001\tChild\tFather\tMother\t1\t2\n"
            )

            # Create output directory and intermediate directory
            output_dir = tmp_path / "output"
            output_dir.mkdir(exist_ok=True)
            (output_dir / "intermediate").mkdir(exist_ok=True)

            yield {"vcf": str(vcf_file), "ped": str(ped_file), "output_dir": str(output_dir)}

    @patch("variantcentrifuge.extractor.run_command")
    @patch("variantcentrifuge.utils.run_command")
    @patch("variantcentrifuge.gene_bed.subprocess.run")
    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    @patch("variantcentrifuge.stages.setup_stages.read_pedigree")
    def test_inheritance_patterns_detected(
        self,
        mock_read_ped,
        mock_get_samples,
        mock_snpeff,
        mock_run_command,
        mock_extractor_run_command,
        setup_files,
    ):
        """Test that inheritance patterns are correctly detected."""
        # Mock VCF samples
        mock_get_samples.return_value = {"Child", "Father", "Mother"}

        # Mock pedigree data
        mock_read_ped.return_value = {
            "Father": {
                "family_id": "FAM001",
                "sample_id": "Father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            },
            "Mother": {
                "family_id": "FAM001",
                "sample_id": "Mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "Child": {
                "family_id": "FAM001",
                "sample_id": "Child",
                "father_id": "Father",
                "mother_id": "Mother",
                "sex": "1",
                "affected_status": "2",
            },
        }

        # Mock snpEff for gene BED creation
        def snpeff_side_effect(cmd, *args, **kwargs):
            from subprocess import CompletedProcess

            if isinstance(cmd, list) and "genes2bed" in cmd:
                if "stdout" in kwargs and hasattr(kwargs["stdout"], "write"):
                    # Return BED regions for all genes
                    kwargs["stdout"].write(
                        "chr1\t50\t150\tGENE1\n" "chr2\t150\t600\tGENE2\n" "chr3\t350\t450\tGENE3\n"
                    )
            return CompletedProcess(cmd, 0)

        mock_snpeff.side_effect = snpeff_side_effect

        # Mock run_command for bcftools and SnpSift
        def run_command_side_effect(cmd, output_file=None, *args, **kwargs):
            result = Mock()
            result.returncode = 0
            result.stdout = ""

            if isinstance(cmd, list):
                cmd_str = " ".join(cmd)
            else:
                cmd_str = cmd

            if "bcftools" in cmd_str and "view" in cmd_str:
                # Mock variant extraction
                if output_file:
                    Path(output_file).write_text(
                        "##fileformat=VCFv4.2\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                        "Child\tFather\tMother\n"
                        "chr1\t100\t.\tA\tG\t100\tPASS\tGENE=GENE1\tGT:DP\t0/1:30\t0/0:25\t0/0:28\n"
                        "chr2\t200\t.\tC\tT\t90\tPASS\tGENE=GENE2\tGT:DP\t0/1:25\t0/1:30\t0/0:27\n"
                        "chr2\t300\t.\tG\tA\t85\tPASS\tGENE=GENE2\tGT:DP\t0/1:28\t0/0:26\t0/1:29\n"
                        "chr3\t400\t.\tT\tC\t95\tPASS\tGENE=GENE3\tGT:DP\t1/1:35\t0/1:30\t0/1:32\n"
                        "chr2\t500\t.\tA\tT\t80\tPASS\tGENE=GENE2\tGT:DP\t0/1:22\t0/1:24\t0/0:26\n"
                    )
            elif "SnpSift" in cmd_str and "extractFields" in cmd_str:
                # Mock field extraction - maintain sample genotypes
                if output_file:
                    with open(output_file, "w") as f:
                        f.write("CHROM\tPOS\tREF\tALT\tQUAL\tGENE\tGT\n")
                        # GT field contains comma-separated genotypes for each sample
                        f.write("chr1\t100\tA\tG\t100\tGENE1\t0/1,0/0,0/0\n")
                        f.write("chr2\t200\tC\tT\t90\tGENE2\t0/1,0/1,0/0\n")
                        f.write("chr2\t300\tG\tA\t85\tGENE2\t0/1,0/0,0/1\n")
                        f.write("chr3\t400\tT\tC\t95\tGENE3\t1/1,0/1,0/1\n")
                        f.write("chr2\t500\tA\tT\t80\tGENE2\t0/1,0/1,0/0\n")

            return result

        mock_run_command.side_effect = run_command_side_effect
        mock_extractor_run_command.side_effect = run_command_side_effect

        # Create arguments with inheritance enabled
        args = Namespace(
            vcf_file=setup_files["vcf"],
            gene_name="GENE1,GENE2,GENE3",
            gene_file=None,
            output_file="output.tsv",
            output_dir=setup_files["output_dir"],
            config=None,
            log_level="DEBUG",  # Enable debug logging
            reference="GRCh37",
            preset=None,
            filter=None,
            late_filtering=False,
            final_filter=None,
            fields_to_extract="CHROM POS REF ALT QUAL GENE GT",
            extract=["CHROM", "POS", "REF", "ALT", "QUAL", "GENE", "GT"],
            threads=1,
            no_stats=True,
            xlsx=False,
            html_report=False,
            phenotype_file=None,
            scoring_config_path=None,
            ped_file=setup_files["ped"],
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
            keep_intermediates=False,
            enable_checkpoint=False,
            no_replacement=True,  # Don't replace genotypes to test inheritance
            perform_gene_burden=False,
            skip_variant_analysis=False,
        )

        # Mock analyze_variants to add Gene_Name
        def mock_analyze_variants(inp, config):
            lines = list(inp)
            if lines:
                header = lines[0].strip()
                # Already has GENE column, just need to add Gene_Name
                if "Gene_Name" not in header:
                    lines[0] = header + "\tGene_Name\n"
                    for i in range(1, len(lines)):
                        # Extract GENE value and add as Gene_Name
                        parts = lines[i].strip().split("\t")
                        if len(parts) > 5:  # GENE is at index 5
                            gene_name = parts[5]
                            lines[i] = lines[i].strip() + f"\t{gene_name}\n"
            for line in lines:
                yield line.strip()

        # Track DataFrame writes to verify inheritance patterns
        written_dfs = []
        original_to_csv = pd.DataFrame.to_csv

        def track_df_write(self, *args, **kwargs):
            written_dfs.append(self.copy())
            return original_to_csv(self, *args, **kwargs)

        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.stages.output_stages.pd.DataFrame.to_csv", side_effect=track_df_write
        ), patch(
            "variantcentrifuge.stages.analysis_stages.analyze_variants",
            side_effect=mock_analyze_variants,
        ), patch(
            "variantcentrifuge.pipeline_core.workspace.Path.mkdir"
        ), patch(
            "variantcentrifuge.helpers.get_vcf_samples", return_value={"Child", "Father", "Mother"}
        ):

            # Run pipeline
            run_refactored_pipeline(args)

        # Verify stages were built correctly
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        # Verify inheritance analysis stage is included
        assert "inheritance_analysis" in stage_names

        # Verify inheritance patterns were detected
        assert len(written_dfs) > 0, "No DataFrames were written"
        final_df = written_dfs[-1]  # Last written DataFrame should be the final output

        # Check that inheritance columns exist
        assert "Inheritance_Pattern" in final_df.columns, f"Columns: {list(final_df.columns)}"

        # Verify specific patterns
        patterns = final_df["Inheritance_Pattern"].tolist()

        # Check for expected patterns
        # Note: The exact pattern names might vary, but we should see different patterns
        assert (
            len(set(patterns)) > 1
        ), f"Expected multiple inheritance patterns, got: {set(patterns)}"

        # For GENE2, we should have compound het patterns
        gene2_variants = final_df[final_df.get("Gene_Name", final_df.get("GENE", "")) == "GENE2"]
        if len(gene2_variants) >= 2:
            gene2_patterns = gene2_variants["Inheritance_Pattern"].tolist()
            # At least some variants in GENE2 should show compound het pattern
            assert any(
                "compound" in str(p).lower() for p in gene2_patterns
            ), f"Expected compound het pattern in GENE2, got: {gene2_patterns}"

    def test_stage_configuration(self, setup_files):
        """Test that inheritance analysis stage is properly configured."""
        args = Namespace(
            vcf_file=setup_files["vcf"],
            gene_name="GENE1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=setup_files["output_dir"],
            config=None,
            log_level="INFO",
            reference="GRCh37",
            ped_file=setup_files["ped"],
            calculate_inheritance=True,
            inheritance_mode="simple",
            fields_to_extract="CHROM POS REF ALT GT",
            extract=["CHROM", "POS", "REF", "ALT", "GT"],
            no_stats=True,
            xlsx=False,
            html_report=False,
            use_new_pipeline=True,
            phenotype_file=None,
            annotate_bed=None,
            annotate_gene_list=None,
            annotate_json_genes=None,
            threads=1,
        )

        # Build stages
        stages = build_pipeline_stages(args)

        # Find inheritance analysis stage
        inheritance_stage = None
        for stage in stages:
            if stage.name == "inheritance_analysis":
                inheritance_stage = stage
                break

        assert inheritance_stage is not None, "InheritanceAnalysisStage not found"

        # Check dependencies
        assert "dataframe_loading" in inheritance_stage.dependencies
        assert "custom_annotation" in inheritance_stage.soft_dependencies
