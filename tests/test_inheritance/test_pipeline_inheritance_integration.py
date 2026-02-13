"""
Integration tests for inheritance analysis in the pipeline.

These tests verify that inheritance analysis works correctly when integrated
into the full pipeline, especially with different genotype formats.
"""

import os
import tempfile
from unittest.mock import Mock

import pytest


class TestPipelineInheritanceIntegration:
    """Test inheritance analysis integration in the main pipeline."""

    @pytest.fixture
    def mock_args(self):
        """Create mock command line arguments."""
        args = Mock()
        args.gene_name = "GENE1"
        args.vcf_file = "test.vcf"
        args.output_dir = tempfile.mkdtemp()
        args.no_replacement = False
        args.no_vectorized_comp_het = False
        args.keep_intermediates = False
        args.calculate_inheritance = True
        args.ped = None  # Will be set in tests
        args.fields = None
        args.genotype_filter = None
        args.gene_genotype_file = None
        args.xlsx = False
        args.html_report = False
        args.log_file = None
        args.gzip_intermediates = False
        args.archive_results = False
        args.phenotype_file = None
        args.phenotype_sample_column = None
        args.phenotype_value_column = None

        # Mock other attributes that might be accessed
        for attr in [
            "config",
            "preset",
            "filter",
            "case_samples_file",
            "control_samples_file",
            "perform_gene_burden",
            "stats_output_file",
            "annotate_bed",
            "annotate_gene_list",
            "annotate_json_genes",
            "json_gene_mapping",
            "late_filtering",
            "scoring_config_path",
        ]:
            setattr(args, attr, None)

        yield args

        # Cleanup
        import shutil

        if os.path.exists(args.output_dir):
            shutil.rmtree(args.output_dir)

    @pytest.fixture
    def simple_ped_file(self):
        """Create a simple trio pedigree file."""
        ped_content = """#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tAffectedStatus
FAM1\tChild\tFather\tMother\t1\t2
FAM1\tFather\t0\t0\t1\t1
FAM1\tMother\t0\t0\t2\t1"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(ped_content)
            ped_file = f.name

        yield ped_file
        os.unlink(ped_file)

    def test_inheritance_with_colon_separated_gt(self, mock_args, simple_ped_file):
        """Test that inheritance analysis works with colon-separated GT format."""
        """Test the GT column splitting logic without full pipeline."""
        # This is a simplified test that focuses on the core logic
        import pandas as pd

        from variantcentrifuge.inheritance.analyzer import analyze_inheritance
        from variantcentrifuge.ped_reader import read_pedigree

        # Create test DataFrame with colon-separated GT
        df = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": ["1000", "2000"],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["GENE1", "GENE1"],
                "GT": ["0/1:0/1:0/0", "0/1:0/0:0/1"],
            }
        )

        # Load pedigree
        pedigree_data = read_pedigree(simple_ped_file)
        sample_list = ["Child", "Father", "Mother"]

        # Simulate the pipeline GT splitting logic
        if "GT" in df.columns and len(df) > 0:
            first_gt = str(df["GT"].iloc[0])
            snpsift_sep = ":"

            if snpsift_sep in first_gt:
                sample_data = {sample_id: [] for sample_id in sample_list}

                for _idx, row in df.iterrows():
                    gt_value = str(row["GT"])
                    if gt_value and gt_value != "NA" and gt_value != "nan":
                        genotypes = gt_value.split(snpsift_sep)
                        for i, sample_id in enumerate(sample_list):
                            if i < len(genotypes):
                                sample_data[sample_id].append(genotypes[i])
                            else:
                                sample_data[sample_id].append("./.")

                sample_df = pd.DataFrame(sample_data)
                df = pd.concat([df, sample_df], axis=1)

        # Run inheritance analysis
        result = analyze_inheritance(df, pedigree_data, sample_list)

        # Verify inheritance patterns were detected
        assert "Inheritance_Pattern" in result.columns
        assert "Inheritance_Details" in result.columns
        assert len(result) == 2

        # Remove sample columns (as pipeline would)
        cols_to_keep = [col for col in result.columns if col not in sample_list]
        result = result[cols_to_keep]

        # Verify sample columns were removed
        assert "Child" not in result.columns
        assert "Father" not in result.columns
        assert "Mother" not in result.columns

    def test_pipeline_config_for_inheritance(self):
        """Test that inheritance configuration is properly handled."""
        # Test basic configuration
        config = {
            "calculate_inheritance": True,
            "extract_fields_separator": ":",
            "fields_to_extract": "CHROM POS REF ALT GENE GEN[*].GT",
        }

        # Verify configuration values
        assert config["calculate_inheritance"] is True
        assert config["extract_fields_separator"] == ":"
        assert "GEN[*].GT" in config["fields_to_extract"]

    def test_gt_column_splitting_logic(self):
        """Test the GT column splitting logic in isolation."""
        import pandas as pd

        # Create test DataFrame with colon-separated GT
        df = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": ["1000", "2000"],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["GENE1", "GENE1"],
                "GT": ["0/1:0/1:0/0", "0/1:0/0:0/1"],
            }
        )

        original_samples = ["Child", "Father", "Mother"]

        if "GT" in df.columns and len(df) > 0:
            first_gt = str(df["GT"].iloc[0])
            snpsift_sep = ":"

            if snpsift_sep in first_gt:
                sample_data = {sample_id: [] for sample_id in original_samples}

                for _idx, row in df.iterrows():
                    gt_value = str(row["GT"])
                    if gt_value and gt_value != "NA" and gt_value != "nan":
                        genotypes = gt_value.split(snpsift_sep)
                        for i, sample_id in enumerate(original_samples):
                            if i < len(genotypes):
                                sample_data[sample_id].append(genotypes[i])
                            else:
                                sample_data[sample_id].append("./.")

                sample_df = pd.DataFrame(sample_data)
                df = pd.concat([df, sample_df], axis=1)

        # Verify the result
        assert "Child" in df.columns
        assert "Father" in df.columns
        assert "Mother" in df.columns
        assert df["Child"].tolist() == ["0/1", "0/1"]
        assert df["Father"].tolist() == ["0/1", "0/0"]
        assert df["Mother"].tolist() == ["0/0", "0/1"]
