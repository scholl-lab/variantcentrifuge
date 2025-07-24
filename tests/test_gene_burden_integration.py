"""
Integration tests for gene burden analysis with sample assignment validation.

Tests the complete pipeline flow including:
- Sample loading from case/control files
- Sample assignment in helpers.assign_case_control_counts
- Gene burden analysis calculations
- Deterministic results across multiple runs
"""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def test_data_paths():
    """Return paths to test data files."""
    fixtures_dir = Path(__file__).parent / "fixtures" / "geneburden" / "output"
    return {
        "vc": fixtures_dir / "enhanced_test_data.vcf.gz",
        "case_samples": fixtures_dir / "case_samples.txt",
        "control_samples": fixtures_dir / "control_samples.txt",
        "genes": fixtures_dir / "test_genes.txt",
    }


@pytest.mark.integration
@pytest.mark.gene_burden
def test_gene_burden_sample_assignment(test_data_paths):
    """Test that case/control samples are correctly assigned and counted."""
    # Skip if test data not available
    for name, path in test_data_paths.items():
        if not path.exists():
            pytest.skip(f"Test data not found: {name} at {path}")

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        # Run gene burden analysis
        cmd = [
            "variantcentrifuge",
            "--vcf-file",
            str(test_data_paths["vc"]),
            "--gene-file",
            str(test_data_paths["genes"]),
            "--case-samples-file",
            str(test_data_paths["case_samples"]),
            "--control-samples-file",
            str(test_data_paths["control_samples"]),
            "--perform-gene-burden",
            "--preset",
            "high_or_moderate",
            "--output-dir",
            str(output_dir),
            "--output-file",
            "test_results.tsv",
            "--use-new-pipeline",
            "--log-level",
            "INFO",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Check gene burden results file exists
        gene_burden_file = output_dir / "test_results.gene_burden.tsv"
        assert gene_burden_file.exists(), "Gene burden results file not created"

        # Load and validate results
        df = pd.read_csv(gene_burden_file, sep="\t")

        # Validate basic structure
        assert len(df) > 0, "No gene burden results generated"
        required_columns = [
            "GENE",
            "proband_count",
            "control_count",
            "proband_allele_count",
            "control_allele_count",
            "raw_p_value",
            "corrected_p_value",
            "odds_ratio",
        ]
        for col in required_columns:
            assert col in df.columns, f"Missing required column: {col}"

        # Validate sample counts are non-zero
        assert df["proband_count"].iloc[0] > 0, "Proband count is zero - case samples not assigned"
        assert (
            df["control_count"].iloc[0] > 0
        ), "Control count is zero - control samples not assigned"

        # Validate reasonable sample counts (based on test data)
        expected_case_count = 40  # From test data
        expected_control_count = 60  # From test data

        assert (
            df["proband_count"].iloc[0] == expected_case_count
        ), f"Expected {expected_case_count} cases, got {df['proband_count'].iloc[0]}"
        assert (
            df["control_count"].iloc[0] == expected_control_count
        ), f"Expected {expected_control_count} controls, got {df['control_count'].iloc[0]}"

        # Validate no negative values in contingency tables
        for _, row in df.iterrows():
            gene = row["GENE"]
            p_count = row["proband_count"]
            c_count = row["control_count"]
            p_alleles = row["proband_allele_count"]
            c_alleles = row["control_allele_count"]

            # Check allele counts don't exceed theoretical maximum
            assert (
                p_alleles <= p_count * 2
            ), f"Gene {gene}: proband alleles ({p_alleles}) > 2 * proband count ({p_count})"
            assert (
                c_alleles <= c_count * 2
            ), f"Gene {gene}: control alleles ({c_alleles}) > 2 * control count ({c_count})"

            # Check reference allele counts would be non-negative
            p_ref_alleles = p_count * 2 - p_alleles
            c_ref_alleles = c_count * 2 - c_alleles
            assert (
                p_ref_alleles >= 0
            ), f"Gene {gene}: negative proband reference alleles ({p_ref_alleles})"
            assert (
                c_ref_alleles >= 0
            ), f"Gene {gene}: negative control reference alleles ({c_ref_alleles})"


@pytest.mark.integration
@pytest.mark.gene_burden
def test_gene_burden_determinism(test_data_paths):
    """Test that gene burden analysis produces identical results across multiple runs."""
    # Skip if test data not available
    for name, path in test_data_paths.items():
        if not path.exists():
            pytest.skip(f"Test data not found: {name} at {path}")

    results = []

    # Run analysis 3 times
    for run_num in range(3):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(test_data_paths["vc"]),
                "--gene-file",
                str(test_data_paths["genes"]),
                "--case-samples-file",
                str(test_data_paths["case_samples"]),
                "--control-samples-file",
                str(test_data_paths["control_samples"]),
                "--perform-gene-burden",
                "--preset",
                "high_or_moderate",
                "--output-dir",
                str(output_dir),
                "--output-file",
                f"run_{run_num}_results.tsv",
                "--use-new-pipeline",
                "--log-level",
                "WARNING",  # Reduce output
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Run {run_num} failed: {result.stderr}"

            # Load results
            gene_burden_file = output_dir / f"run_{run_num}_results.gene_burden.tsv"
            assert gene_burden_file.exists(), f"Run {run_num} gene burden file not created"

            df = pd.read_csv(gene_burden_file, sep="\t")
            df = df.sort_values("GENE").reset_index(drop=True)  # Ensure consistent order
            results.append(df)

    # Compare all results for determinism
    base_df = results[0]

    for i, df in enumerate(results[1:], 1):
        # Check same number of genes
        assert len(base_df) == len(
            df
        ), f"Run 1 and run {i + 1} have different number of genes: {len(base_df)} vs {len(df)}"

        # Check identical gene order
        pd.testing.assert_series_equal(
            base_df["GENE"],
            df["GENE"],
            check_names=False,
            obj=f"Gene order differs between run 1 and run {i + 1}",
        )

        # Check identical statistical results
        for col in [
            "proband_count",
            "control_count",
            "proband_allele_count",
            "control_allele_count",
            "raw_p_value",
            "corrected_p_value",
        ]:
            pd.testing.assert_series_equal(
                base_df[col],
                df[col],
                check_names=False,
                obj=f"Column {col} differs between run 1 and run {i + 1}",
            )

        # Check odds ratios (may have floating point differences)
        pd.testing.assert_series_equal(
            base_df["odds_ratio"],
            df["odds_ratio"],
            check_names=False,
            check_exact=False,
            rtol=1e-10,  # Very tight tolerance
            obj=f"Odds ratios differ between run 1 and run {i + 1}",
        )


@pytest.mark.unit
def test_assign_case_control_counts_unit():
    """Unit test for the assign_case_control_counts helper function."""
    import pandas as pd

    from variantcentrifuge.helpers import assign_case_control_counts

    # Create test DataFrame with GT column
    test_data = {
        "GENE": ["GENE1", "GENE1", "GENE2"],
        "GT": [
            "CASE_001(1/1);CASE_002(0/1);CTRL_001(0/0);CTRL_002(0/1)",
            "CASE_001(0/1);CASE_003(0/0);CTRL_001(1/1);CTRL_003(0/1)",
            "CASE_002(1/1);CTRL_001(0/1);CTRL_002(0/0)",
        ],
    }
    df = pd.DataFrame(test_data)

    # Define sample groups
    case_samples = {"CASE_001", "CASE_002", "CASE_003"}
    control_samples = {"CTRL_001", "CTRL_002", "CTRL_003"}
    all_samples = case_samples | control_samples

    # Run function
    result_df = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    # Validate structure
    expected_columns = [
        "proband_count",
        "control_count",
        "proband_variant_count",
        "control_variant_count",
        "proband_allele_count",
        "control_allele_count",
        "proband_homozygous_count",
        "control_homozygous_count",
    ]
    for col in expected_columns:
        assert col in result_df.columns, f"Missing column: {col}"

    # Validate counts
    assert result_df["proband_count"].iloc[0] == 3, "Wrong proband count"
    assert result_df["control_count"].iloc[0] == 3, "Wrong control count"

    # Validate first row: CASE_001(1/1), CASE_002(0/1), CTRL_002(0/1)
    # Cases: 1 homozygous (2 alleles) + 1 heterozygous (1 allele) = 3 alleles,
    # 2 samples with variants
    # Controls: 1 heterozygous (1 allele) = 1 allele, 1 sample with variant
    assert result_df["proband_variant_count"].iloc[0] == 2, (
        f"Row 0: Expected 2 case samples with variants, "
        f"got {result_df['proband_variant_count'].iloc[0]}"
    )
    assert result_df["control_variant_count"].iloc[0] == 1, (
        f"Row 0: Expected 1 control sample with variant, "
        f"got {result_df['control_variant_count'].iloc[0]}"
    )
    assert (
        result_df["proband_allele_count"].iloc[0] == 3
    ), f"Row 0: Expected 3 case alleles, got {result_df['proband_allele_count'].iloc[0]}"
    assert (
        result_df["control_allele_count"].iloc[0] == 1
    ), f"Row 0: Expected 1 control allele, got {result_df['control_allele_count'].iloc[0]}"


@pytest.mark.unit
def test_gene_burden_analysis_unit():
    """Unit test for the perform_gene_burden_analysis function."""
    import pandas as pd

    from variantcentrifuge.gene_burden import perform_gene_burden_analysis

    # Create test data with known counts
    test_data = {
        "GENE": ["GENE1", "GENE1", "GENE2", "GENE2"],
        "proband_count": [10, 10, 10, 10],
        "control_count": [10, 10, 10, 10],
        "proband_variant_count": [2, 1, 3, 2],
        "control_variant_count": [1, 0, 1, 1],
        "proband_allele_count": [3, 1, 4, 3],
        "control_allele_count": [1, 0, 2, 1],
    }
    df = pd.DataFrame(test_data)

    # Test configuration
    config = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    # Run analysis
    result = perform_gene_burden_analysis(df, config)

    # Validate structure
    assert len(result) == 2, "Should have 2 genes"
    assert "GENE" in result.columns
    assert "odds_ratio" in result.columns
    assert "raw_p_value" in result.columns
    assert "corrected_p_value" in result.columns

    # Validate gene aggregation
    gene1_row = result[result["GENE"] == "GENE1"].iloc[0]
    gene2_row = result[result["GENE"] == "GENE2"].iloc[0]

    # GENE1: sum of variants = 3+1=4 case alleles, 1+0=1 control allele
    assert (
        gene1_row["proband_allele_count"] == 4
    ), f"GENE1: Expected 4 case alleles, got {gene1_row['proband_allele_count']}"
    assert (
        gene1_row["control_allele_count"] == 1
    ), f"GENE1: Expected 1 control allele, got {gene1_row['control_allele_count']}"

    # GENE2: sum of variants = 4+3=7 case alleles, 2+1=3 control alleles
    assert (
        gene2_row["proband_allele_count"] == 7
    ), f"GENE2: Expected 7 case alleles, got {gene2_row['proband_allele_count']}"
    assert (
        gene2_row["control_allele_count"] == 3
    ), f"GENE2: Expected 3 control alleles, got {gene2_row['control_allele_count']}"

    # Validate no negative reference allele counts
    for _, row in result.iterrows():
        p_count = row["proband_count"]
        c_count = row["control_count"]
        p_alleles = row["proband_allele_count"]
        c_alleles = row["control_allele_count"]

        p_ref = p_count * 2 - p_alleles
        c_ref = c_count * 2 - c_alleles

        assert p_ref >= 0, f"Negative case reference alleles for {row['GENE']}: {p_ref}"
        assert c_ref >= 0, f"Negative control reference alleles for {row['GENE']}: {c_ref}"
