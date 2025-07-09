"""Integration test for parallel processing with the new pipeline.

This test verifies that the ParallelCompleteProcessingStage correctly
processes data in parallel and produces the same results as sequential processing.
"""

import tempfile
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from variantcentrifuge.pipeline_refactored import build_pipeline_stages, run_refactored_pipeline


@pytest.mark.slow
class TestParallelProcessing:
    """Test parallel processing in the refactored pipeline."""

    @pytest.fixture
    def setup_files(self):
        """Create test files for the pipeline."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            # Create test VCF file with multiple genes
            vcf_file = tmp_path / "test.vcf"
            vcf_file.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                # GENE1 variants
                "chr1\t100\t.\tA\tG\t100\tPASS\tGENE=GENE1;AC=2\tGT:DP\t0/1:30\n"
                "chr1\t200\t.\tC\tT\t90\tPASS\tGENE=GENE1;AC=1\tGT:DP\t0/1:25\n"
                "chr1\t300\t.\tG\tA\t80\tPASS\tGENE=GENE1;AC=3\tGT:DP\t1/1:40\n"
                # GENE2 variants
                "chr2\t1000\t.\tT\tC\t95\tPASS\tGENE=GENE2;AC=2\tGT:DP\t0/1:35\n"
                "chr2\t2000\t.\tA\tG\t85\tPASS\tGENE=GENE2;AC=1\tGT:DP\t0/1:28\n"
                # GENE3 variants
                "chr3\t5000\t.\tC\tG\t110\tPASS\tGENE=GENE3;AC=2\tGT:DP\t0/1:45\n"
                "chr3\t6000\t.\tT\tA\t75\tPASS\tGENE=GENE3;AC=4\tGT:DP\t1/1:30\n"
            )

            # Create output directory and intermediate directory
            output_dir = tmp_path / "output"
            output_dir.mkdir(exist_ok=True)
            (output_dir / "intermediate").mkdir(exist_ok=True)

            yield {"vcf": str(vcf_file), "output_dir": str(output_dir), "tmp_path": tmp_path}

    @patch("variantcentrifuge.extractor.run_command")
    @patch("variantcentrifuge.utils.run_command")
    @patch("variantcentrifuge.gene_bed.subprocess.run")
    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_parallel_vs_sequential_results(
        self,
        mock_get_samples,
        mock_snpeff,
        mock_run_command,
        mock_extractor_run_command,
        setup_files,
    ):
        """Test that parallel processing produces same results as sequential."""
        # Mock VCF samples
        mock_get_samples.return_value = ["Sample1"]

        # Mock snpEff for gene BED creation - create separate regions for each gene
        def snpeff_side_effect(cmd, *args, **kwargs):
            from subprocess import CompletedProcess

            if isinstance(cmd, list) and "genes2bed" in cmd:
                if "stdout" in kwargs and hasattr(kwargs["stdout"], "write"):
                    # Create BED regions for three genes
                    kwargs["stdout"].write(
                        "chr1\t50\t350\tGENE1\n"
                        "chr2\t950\t2050\tGENE2\n"
                        "chr3\t4950\t6050\tGENE3\n"
                    )
            return CompletedProcess(cmd, 0)

        mock_snpeff.side_effect = snpeff_side_effect

        # Track which files were processed
        processed_files = {"sequential": [], "parallel": []}

        # Mock run_command for bcftools and SnpSift
        def run_command_side_effect(cmd, output_file=None, *args, **kwargs):
            result = Mock()
            result.returncode = 0
            result.stdout = ""

            if isinstance(cmd, list):
                cmd_str = " ".join(cmd)
            else:
                cmd_str = cmd

            # Track which mode we're in based on thread count
            mode = "parallel" if "chunk_" in cmd_str else "sequential"

            if "bcftools" in cmd_str and "view" in cmd_str:
                # Mock variant extraction
                if output_file:
                    # Extract which region is being processed
                    bed_file = None
                    for i, arg in enumerate(cmd):
                        if arg == "-R" and i + 1 < len(cmd):
                            bed_file = cmd[i + 1]
                            break

                    # Write appropriate variants based on BED file
                    if bed_file and "chunk_" in bed_file:
                        # Parallel mode - write subset based on chunk
                        chunk_num = int(Path(bed_file).stem.split("_")[1])
                        if chunk_num == 0:
                            # First chunk - GENE1
                            Path(output_file).write_text(
                                "##fileformat=VCFv4.2\n"
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                                "chr1\t100\t.\tA\tG\t100\tPASS\tGENE=GENE1;AC=2\tGT:DP\t0/1:30\n"
                                "chr1\t200\t.\tC\tT\t90\tPASS\tGENE=GENE1;AC=1\tGT:DP\t0/1:25\n"
                                "chr1\t300\t.\tG\tA\t80\tPASS\tGENE=GENE1;AC=3\tGT:DP\t1/1:40\n"
                            )
                        elif chunk_num == 1:
                            # Second chunk - GENE2
                            Path(output_file).write_text(
                                "##fileformat=VCFv4.2\n"
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                                "chr2\t1000\t.\tT\tC\t95\tPASS\tGENE=GENE2;AC=2\tGT:DP\t0/1:35\n"
                                "chr2\t2000\t.\tA\tG\t85\tPASS\tGENE=GENE2;AC=1\tGT:DP\t0/1:28\n"
                            )
                        else:
                            # Third chunk - GENE3
                            Path(output_file).write_text(
                                "##fileformat=VCFv4.2\n"
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                                "chr3\t5000\t.\tC\tG\t110\tPASS\tGENE=GENE3;AC=2\tGT:DP\t0/1:45\n"
                                "chr3\t6000\t.\tT\tA\t75\tPASS\tGENE=GENE3;AC=4\tGT:DP\t1/1:30\n"
                            )
                    else:
                        # Sequential mode - write all variants
                        Path(output_file).write_text(
                            "##fileformat=VCFv4.2\n"
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                            "chr1\t100\t.\tA\tG\t100\tPASS\tGENE=GENE1;AC=2\tGT:DP\t0/1:30\n"
                            "chr1\t200\t.\tC\tT\t90\tPASS\tGENE=GENE1;AC=1\tGT:DP\t0/1:25\n"
                            "chr1\t300\t.\tG\tA\t80\tPASS\tGENE=GENE1;AC=3\tGT:DP\t1/1:40\n"
                            "chr2\t1000\t.\tT\tC\t95\tPASS\tGENE=GENE2;AC=2\tGT:DP\t0/1:35\n"
                            "chr2\t2000\t.\tA\tG\t85\tPASS\tGENE=GENE2;AC=1\tGT:DP\t0/1:28\n"
                            "chr3\t5000\t.\tC\tG\t110\tPASS\tGENE=GENE3;AC=2\tGT:DP\t0/1:45\n"
                            "chr3\t6000\t.\tT\tA\t75\tPASS\tGENE=GENE3;AC=4\tGT:DP\t1/1:30\n"
                        )
                    processed_files[mode].append(("extract", output_file))

            elif "SnpSift" in cmd_str:
                if "filter" in cmd_str and output_file:
                    # Mock filtering - filter by QUAL >= 80 and AC < 3
                    input_file = None
                    for arg in cmd:
                        if arg.endswith(".vcf.gz") and arg != output_file:
                            input_file = arg
                            break

                    if input_file and Path(input_file).exists():
                        # Read input and filter
                        lines = Path(input_file).read_text().strip().split("\n")
                        output_lines = []
                        for line in lines:
                            if line.startswith("#"):
                                output_lines.append(line)
                            else:
                                parts = line.split("\t")
                                if len(parts) > 7:
                                    qual = float(parts[5])
                                    info = parts[7]
                                    ac = int(info.split("AC=")[1].split(";")[0])
                                    if qual >= 80 and ac < 3:
                                        output_lines.append(line)
                        Path(output_file).write_text("\n".join(output_lines) + "\n")
                        processed_files[mode].append(("filter", output_file))

                elif "extractFields" in cmd_str and output_file:
                    # Mock field extraction
                    with open(output_file, "w") as f:
                        f.write("CHROM\tPOS\tREF\tALT\tQUAL\tGENE\tAC\tGT\n")
                        # Parse the input VCF to extract fields
                        input_file = None
                        for arg in cmd:
                            if arg.endswith(".vcf.gz") and arg != output_file:
                                input_file = arg
                                break

                        if input_file and Path(input_file).exists():
                            lines = Path(input_file).read_text().strip().split("\n")
                            for line in lines:
                                if not line.startswith("#"):
                                    parts = line.split("\t")
                                    if len(parts) >= 10:
                                        chrom = parts[0]
                                        pos = parts[1]
                                        ref = parts[3]
                                        alt = parts[4]
                                        qual = parts[5]
                                        info = parts[7]
                                        gt = parts[9].split(":")[0]

                                        # Extract GENE and AC from INFO
                                        gene = (
                                            info.split("GENE=")[1].split(";")[0]
                                            if "GENE=" in info
                                            else ""
                                        )
                                        ac = (
                                            info.split("AC=")[1].split(";")[0]
                                            if "AC=" in info
                                            else ""
                                        )

                                        f.write(
                                            f"{chrom}\t{pos}\t{ref}\t{alt}\t{qual}\t{gene}\t"
                                            f"{ac}\tSample1:{gt}\n"
                                        )
                        processed_files[mode].append(("extract", output_file))

            return result

        mock_run_command.side_effect = run_command_side_effect
        mock_extractor_run_command.side_effect = run_command_side_effect

        # Common arguments
        base_args = {
            "vcf_file": setup_files["vcf"],
            "gene_name": "GENE1,GENE2,GENE3",
            "gene_file": None,
            "output_file": "output.tsv",
            "output_dir": setup_files["output_dir"],
            "config": None,
            "log_level": "DEBUG",
            "reference": "GRCh37",
            "preset": None,
            "filter": "(QUAL >= 80) & (AC < 3)",  # Filter expression
            "late_filtering": False,
            "final_filter": None,
            "fields_to_extract": "CHROM POS REF ALT QUAL GENE AC GT",
            "extract": ["CHROM", "POS", "REF", "ALT", "QUAL", "GENE", "AC", "GT"],
            "no_stats": True,
            "xlsx": False,
            "html_report": False,
            "phenotype_file": None,
            "scoring_config_path": None,
            "ped_file": None,
            "calculate_inheritance": False,
            "no_replacement": True,
            "use_new_pipeline": True,
            "keep_intermediates": True,
        }

        # Test sequential processing (threads=1)
        sequential_args = Namespace(**{**base_args, "threads": 1})

        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.pipeline_core.workspace.Path.mkdir"
        ), patch(
            "variantcentrifuge.helpers.get_vcf_samples", return_value=["Sample1"]
        ):

            # Reset processed files
            processed_files["sequential"] = []
            run_refactored_pipeline(sequential_args)

        # Read sequential output
        sequential_output = Path(setup_files["output_dir"]) / "output.tsv"
        assert sequential_output.exists()
        sequential_df = pd.read_csv(sequential_output, sep="\t")

        # Test parallel processing (threads=3)
        parallel_args = Namespace(**{**base_args, "threads": 3})

        # Create new output file
        parallel_args.output_file = "output_parallel.tsv"

        with patch(
            "variantcentrifuge.stages.processing_stages.Path.exists", return_value=True
        ), patch("variantcentrifuge.stages.processing_stages.Path.touch"), patch(
            "variantcentrifuge.pipeline_core.workspace.Path.mkdir"
        ), patch(
            "variantcentrifuge.helpers.get_vcf_samples", return_value=["Sample1"]
        ):

            # Reset processed files
            processed_files["parallel"] = []
            run_refactored_pipeline(parallel_args)

        # Read parallel output
        parallel_output = Path(setup_files["output_dir"]) / "output_parallel.tsv"
        assert parallel_output.exists()
        parallel_df = pd.read_csv(parallel_output, sep="\t")

        # Compare results
        # Sort both DataFrames by CHROM and POS to ensure consistent ordering
        sequential_df = sequential_df.sort_values(["CHROM", "POS"]).reset_index(drop=True)
        parallel_df = parallel_df.sort_values(["CHROM", "POS"]).reset_index(drop=True)

        # Should have same number of rows
        assert len(sequential_df) == len(
            parallel_df
        ), f"Different row counts: sequential={len(sequential_df)}, parallel={len(parallel_df)}"

        # Should have same columns
        assert list(sequential_df.columns) == list(parallel_df.columns)

        # Should have same data
        pd.testing.assert_frame_equal(sequential_df, parallel_df)

        # Verify filtering worked (QUAL >= 80 and AC < 3)
        assert all(sequential_df["QUAL"].astype(float) >= 80)
        assert all(sequential_df["AC"].astype(int) < 3)

        # Verify we have variants from all genes
        genes = set(sequential_df["GENE"])
        assert "GENE1" in genes
        assert "GENE2" in genes
        assert "GENE3" in genes

    def test_stage_configuration_parallel(self, setup_files):
        """Test that parallel processing uses correct stages."""
        args = Namespace(
            vcf_file=setup_files["vcf"],
            gene_name="GENE1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=setup_files["output_dir"],
            config=None,
            log_level="INFO",
            reference="GRCh37",
            threads=4,  # Parallel processing
            filter="QUAL >= 30",
            extract=["CHROM", "POS", "REF", "ALT"],
            no_stats=True,
            xlsx=False,
            html_report=False,
            use_new_pipeline=True,
        )

        # Build stages
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        # Should use ParallelCompleteProcessingStage instead of individual stages
        assert "parallel_complete_processing" in stage_names
        assert "variant_extraction" not in stage_names
        assert "parallel_variant_extraction" not in stage_names
        assert "snpsift_filtering" not in stage_names
        assert "field_extraction" not in stage_names

    def test_stage_configuration_sequential(self, setup_files):
        """Test that sequential processing uses correct stages."""
        args = Namespace(
            vcf_file=setup_files["vcf"],
            gene_name="GENE1",
            gene_file=None,
            output_file="output.tsv",
            output_dir=setup_files["output_dir"],
            config=None,
            log_level="INFO",
            reference="GRCh37",
            threads=1,  # Sequential processing
            filter="QUAL >= 30",
            extract=["CHROM", "POS", "REF", "ALT"],
            no_stats=True,
            xlsx=False,
            html_report=False,
            use_new_pipeline=True,
            late_filtering=False,
        )

        # Build stages
        stages = build_pipeline_stages(args)
        stage_names = [s.name for s in stages]

        # Should use individual stages, not parallel processing
        assert "parallel_complete_processing" not in stage_names
        assert "variant_extraction" in stage_names
        assert "snpsift_filtering" in stage_names
        assert "field_extraction" in stage_names
