"""
Test inheritance analysis with colon-separated genotype format.

This tests the scenario where GEN[*].GT is used in field extraction,
resulting in a single GT column with colon-separated values for all samples.
"""

import os
import tempfile

import pandas as pd
import pytest

from variantcentrifuge.inheritance.analyzer import analyze_inheritance
from variantcentrifuge.ped_reader import read_pedigree


class TestColonSeparatedGenotypes:
    """Test inheritance analysis with colon-separated GT column format."""

    @pytest.fixture
    def simple_pedigree_file(self):
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

    @pytest.fixture
    def large_pedigree_file(self):
        """Create a pedigree with many samples to test performance."""
        n_samples = 100
        ped_lines = ["#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tAffectedStatus"]

        for i in range(n_samples):
            sample_id = f"Sample{i}"
            ped_lines.append(f"FAM1\t{sample_id}\t0\t0\t1\t2")

        ped_content = "\n".join(ped_lines)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(ped_content)
            ped_file = f.name

        yield ped_file
        os.unlink(ped_file)

    def test_compound_het_with_colon_separated_gt(self, simple_pedigree_file):
        """Test compound heterozygous detection with colon-separated GT column."""
        # Create test data as it comes from GEN[*].GT extraction
        data = {
            "CHROM": ["1", "1"],
            "POS": ["1000", "2000"],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "GENE": ["GENE1", "GENE1"],
            "GT": ["0/1:0/1:0/0", "0/1:0/0:0/1"],  # Child:Father:Mother
        }
        df = pd.DataFrame(data)

        # Load pedigree
        pedigree_data = read_pedigree(simple_pedigree_file)
        sample_list = ["Child", "Father", "Mother"]

        # Split GT column (simulating pipeline behavior)
        sample_data = {sample_id: [] for sample_id in sample_list}
        for idx, row in df.iterrows():
            gt_value = str(row["GT"])
            if gt_value and gt_value != "NA":
                genotypes = gt_value.split(":")
                for i, sample_id in enumerate(sample_list):
                    if i < len(genotypes):
                        sample_data[sample_id].append(genotypes[i])
                    else:
                        sample_data[sample_id].append("./.")

        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        # Run inheritance analysis
        result = analyze_inheritance(df, pedigree_data, sample_list)

        # Check results
        assert len(result) == 2
        assert all(result["Inheritance_Pattern"] == "compound_heterozygous")

    def test_de_novo_with_colon_separated_gt(self, simple_pedigree_file):
        """Test de novo detection with colon-separated GT column."""
        data = {
            "CHROM": ["1"],
            "POS": ["1000"],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["GENE1"],
            "GT": ["0/1:0/0:0/0"],  # Child:Father:Mother - de novo in child
        }
        df = pd.DataFrame(data)

        # Load pedigree
        pedigree_data = read_pedigree(simple_pedigree_file)
        sample_list = ["Child", "Father", "Mother"]

        # Split GT column
        sample_data = {sample_id: [] for sample_id in sample_list}
        for idx, row in df.iterrows():
            genotypes = str(row["GT"]).split(":")
            for i, sample_id in enumerate(sample_list):
                sample_data[sample_id].append(genotypes[i] if i < len(genotypes) else "./.")

        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        # Run inheritance analysis
        result = analyze_inheritance(df, pedigree_data, sample_list)

        assert len(result) == 1
        assert result.iloc[0]["Inheritance_Pattern"] == "de_novo"

    def test_recessive_with_colon_separated_gt(self, simple_pedigree_file):
        """Test autosomal recessive detection with colon-separated GT column."""
        data = {
            "CHROM": ["1"],
            "POS": ["1000"],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["GENE1"],
            "GT": ["1/1:0/1:0/1"],  # Child:Father:Mother - recessive
        }
        df = pd.DataFrame(data)

        pedigree_data = read_pedigree(simple_pedigree_file)
        sample_list = ["Child", "Father", "Mother"]

        # Split GT column
        sample_data = {sample_id: [] for sample_id in sample_list}
        for idx, row in df.iterrows():
            genotypes = str(row["GT"]).split(":")
            for i, sample_id in enumerate(sample_list):
                sample_data[sample_id].append(genotypes[i] if i < len(genotypes) else "./.")

        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        result = analyze_inheritance(df, pedigree_data, sample_list)

        assert len(result) == 1
        assert result.iloc[0]["Inheritance_Pattern"] == "autosomal_recessive"

    def test_missing_genotypes_in_colon_separated(self, simple_pedigree_file):
        """Test handling of missing genotypes in colon-separated format."""
        data = {
            "CHROM": ["1", "2"],
            "POS": ["1000", "2000"],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "GENE": ["GENE1", "GENE2"],
            "GT": ["0/1:0/0", "0/1:0/1:./."],  # Missing mother in second variant
        }
        df = pd.DataFrame(data)

        pedigree_data = read_pedigree(simple_pedigree_file)
        sample_list = ["Child", "Father", "Mother"]

        # Split GT column
        sample_data = {sample_id: [] for sample_id in sample_list}
        for idx, row in df.iterrows():
            genotypes = str(row["GT"]).split(":")
            for i, sample_id in enumerate(sample_list):
                sample_data[sample_id].append(genotypes[i] if i < len(genotypes) else "./.")

        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        result = analyze_inheritance(df, pedigree_data, sample_list)

        assert len(result) == 2
        # First variant should have a pattern, second might be unknown due to missing data
        assert result.iloc[0]["Inheritance_Pattern"] != "error"

    def test_large_sample_count_performance(self, large_pedigree_file):
        """Test that analysis works efficiently with many samples."""
        import warnings

        n_samples = 100
        sample_list = [f"Sample{i}" for i in range(n_samples)]

        # Create data with many samples
        data = {
            "CHROM": ["1", "1"],
            "POS": ["1000", "2000"],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "GENE": ["GENE1", "GENE1"],
            "GT": [
                ":".join(["0/1" if i < 2 else "0/0" for i in range(n_samples)]),
                ":".join(["0/1" if i == 0 or i == 50 else "0/0" for i in range(n_samples)]),
            ],
        }
        df = pd.DataFrame(data)

        pedigree_data = read_pedigree(large_pedigree_file)

        # Split GT column
        sample_data = {sample_id: [] for sample_id in sample_list}
        for idx, row in df.iterrows():
            genotypes = str(row["GT"]).split(":")
            for i, sample_id in enumerate(sample_list):
                sample_data[sample_id].append(genotypes[i] if i < len(genotypes) else "./.")

        # This should not produce fragmentation warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sample_df = pd.DataFrame(sample_data)
            df = pd.concat([df, sample_df], axis=1)

            # Check for performance warnings
            perf_warnings = [
                warning for warning in w if "PerformanceWarning" in str(warning.category)
            ]
            assert len(perf_warnings) == 0, f"Got {len(perf_warnings)} performance warnings"

        # Run analysis
        result = analyze_inheritance(df, pedigree_data, sample_list)

        assert len(result) == 2
        assert "Inheritance_Pattern" in result.columns
        assert "Inheritance_Details" in result.columns

    def test_different_separators(self, simple_pedigree_file):
        """Test handling of different field separators."""
        # Test with comma separator
        data = {
            "CHROM": ["1"],
            "POS": ["1000"],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["GENE1"],
            "GT": ["0/1,0/0,0/0"],  # Using comma instead of colon
        }
        df = pd.DataFrame(data)

        pedigree_data = read_pedigree(simple_pedigree_file)
        sample_list = ["Child", "Father", "Mother"]

        # Split with comma
        sample_data = {sample_id: [] for sample_id in sample_list}
        for idx, row in df.iterrows():
            genotypes = str(row["GT"]).split(",")
            for i, sample_id in enumerate(sample_list):
                sample_data[sample_id].append(genotypes[i] if i < len(genotypes) else "./.")

        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        result = analyze_inheritance(df, pedigree_data, sample_list)

        assert len(result) == 1
        assert result.iloc[0]["Inheritance_Pattern"] == "de_novo"

    def test_empty_gt_column(self, simple_pedigree_file):
        """Test handling of empty or NA GT values."""
        data = {
            "CHROM": ["1", "2", "3"],
            "POS": ["1000", "2000", "3000"],
            "REF": ["A", "G", "C"],
            "ALT": ["T", "C", "T"],
            "GENE": ["GENE1", "GENE2", "GENE3"],
            "GT": ["0/1:0/0:0/0", "NA", ""],  # Various empty values
        }
        df = pd.DataFrame(data)

        pedigree_data = read_pedigree(simple_pedigree_file)
        sample_list = ["Child", "Father", "Mother"]

        # Split GT column with proper handling of NA/empty
        sample_data = {sample_id: [] for sample_id in sample_list}
        for idx, row in df.iterrows():
            gt_value = str(row["GT"])
            if gt_value and gt_value != "NA" and gt_value != "nan" and gt_value != "":
                genotypes = gt_value.split(":")
                for i, sample_id in enumerate(sample_list):
                    sample_data[sample_id].append(genotypes[i] if i < len(genotypes) else "./.")
            else:
                for sample_id in sample_list:
                    sample_data[sample_id].append("./.")

        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        result = analyze_inheritance(df, pedigree_data, sample_list)

        assert len(result) == 3
        # First variant should have a pattern, others should be unknown
        assert result.iloc[0]["Inheritance_Pattern"] == "de_novo"
        assert result.iloc[1]["Inheritance_Pattern"] == "unknown"
        assert result.iloc[2]["Inheritance_Pattern"] == "unknown"
