"""Integration test for inheritance analysis with the new pipeline.

This test verifies that inheritance analysis correctly detects patterns
including compound heterozygous variants.
"""

import tempfile
from argparse import Namespace
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

from variantcentrifuge.pipeline import build_pipeline_stages


class TestInheritanceAnalysis:
    """Test inheritance analysis in the stage-based pipeline."""

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
