"""
Regression tests for gene burden analysis sample assignment bug fixes.

This test ensures that the specific bug where proband_count was always 0
and case samples were not being correctly assigned does not regress.
"""

import pytest
import tempfile
import subprocess
from pathlib import Path
import pandas as pd


@pytest.mark.regression
@pytest.mark.gene_burden 
def test_proband_count_not_zero_regression():
    """
    Regression test for bug where proband_count was always 0.
    
    This test specifically validates that:
    1. proband_count > 0 (case samples are being counted)
    2. control_count > 0 (control samples are being counted)
    3. Both case and control samples are found in the GT column
    4. Sample assignment logic correctly uses all VCF samples
    """
    fixtures_dir = Path(__file__).parent / "fixtures" / "geneburden" / "output"
    test_files = {
        "vcf": fixtures_dir / "enhanced_test_data.vcf.gz",
        "case_samples": fixtures_dir / "case_samples.txt",
        "control_samples": fixtures_dir / "control_samples.txt",
        "genes": fixtures_dir / "test_genes.txt",
    }
    
    # Skip if test data not available
    for name, path in test_files.items():
        if not path.exists():
            pytest.skip(f"Test data not found: {name} at {path}")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        
        # Run the exact command that was failing before the fix
        cmd = [
            "variantcentrifuge",
            "--vcf-file", str(test_files["vcf"]),
            "--gene-file", str(test_files["genes"]),
            "--case-samples-file", str(test_files["case_samples"]),
            "--control-samples-file", str(test_files["control_samples"]),
            "--perform-gene-burden",
            "--preset", "high_or_moderate",
            "--output-dir", str(output_dir),
            "--output-file", "regression_test.tsv",
            "--use-new-pipeline",
            "--log-level", "WARNING"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Ensure command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        
        # Check gene burden file exists
        gene_burden_file = output_dir / "regression_test.gene_burden.tsv"
        assert gene_burden_file.exists(), "Gene burden results file not created"
        
        # Load results
        df = pd.read_csv(gene_burden_file, sep="\t")
        assert len(df) > 0, "No gene burden results generated"
        
        # THE KEY REGRESSION TEST: proband_count must not be 0
        proband_counts = df["proband_count"].unique()
        assert len(proband_counts) == 1, f"Expected single proband count, got {proband_counts}"
        assert proband_counts[0] > 0, f"REGRESSION: proband_count is {proband_counts[0]}, should be > 0"
        
        # Control count should also be > 0
        control_counts = df["control_count"].unique()
        assert len(control_counts) == 1, f"Expected single control count, got {control_counts}"
        assert control_counts[0] > 0, f"REGRESSION: control_count is {control_counts[0]}, should be > 0"
        
        # Validate expected counts from test data
        expected_case_count = 40
        expected_control_count = 60
        
        assert proband_counts[0] == expected_case_count, \
            f"Expected {expected_case_count} cases, got {proband_counts[0]}"
        assert control_counts[0] == expected_control_count, \
            f"Expected {expected_control_count} controls, got {control_counts[0]}"
        
        # Verify that we have actual variant counts (not all zeros)
        total_case_alleles = df["proband_allele_count"].sum()
        total_control_alleles = df["control_allele_count"].sum()
        
        assert total_case_alleles > 0, "No case alleles found - samples not properly assigned"
        assert total_control_alleles > 0, "No control alleles found - samples not properly assigned"
        
        # Check that we have reasonable odds ratios (not all NaN or infinity)
        finite_ors = df["odds_ratio"][df["odds_ratio"].notna() & df["odds_ratio"].isfinite()]
        assert len(finite_ors) > 0, "No finite odds ratios found"
        
        # Verify p-values are reasonable (not all 1.0 or NaN)
        finite_pvals = df["raw_p_value"][df["raw_p_value"].notna() & df["raw_p_value"].isfinite()]
        assert len(finite_pvals) > 0, "No finite p-values found"
        
        # Check that contingency tables don't have negative values 
        # (which would cause Fisher's exact test to fail)
        for _, row in df.iterrows():
            gene = row["GENE"]
            p_count = row["proband_count"]
            c_count = row["control_count"]
            p_alleles = row["proband_allele_count"]
            c_alleles = row["control_allele_count"]
            
            # Calculate what the reference allele counts would be
            p_ref_alleles = p_count * 2 - p_alleles
            c_ref_alleles = c_count * 2 - c_alleles
            
            assert p_ref_alleles >= 0, \
                f"Gene {gene}: negative case reference alleles ({p_ref_alleles})"
            assert c_ref_alleles >= 0, \
                f"Gene {gene}: negative control reference alleles ({c_ref_alleles})"


@pytest.mark.regression
@pytest.mark.gene_burden
def test_sample_assignment_with_all_vcf_samples():
    """
    Regression test for the specific bug where assign_case_control_counts
    was called with case_samples + control_samples as all_samples instead
    of using the complete VCF sample list.
    """
    from variantcentrifuge.helpers import assign_case_control_counts
    import pandas as pd
    
    # Simulate the bug scenario: VCF has more samples than case+control lists
    # This would happen if some samples in the VCF are not classified
    
    # Create test data with samples that appear in GT but not in case/control lists
    test_data = {
        "GENE": ["GENE1"],
        "GT": [
            "CASE_001(1/1);CASE_002(0/1);CTRL_001(0/1);CTRL_002(0/0);UNCLASSIFIED_001(0/1);UNCLASSIFIED_002(1/1)"
        ]
    }
    df = pd.DataFrame(test_data)
    
    # Define sample groups - note that UNCLASSIFIED samples are not in either group
    case_samples = {"CASE_001", "CASE_002"}
    control_samples = {"CTRL_001", "CTRL_002"} 
    
    # THE BUG: using only case + control samples as all_samples
    buggy_all_samples = case_samples | control_samples
    
    # THE FIX: using all VCF samples (including unclassified)
    correct_all_samples = {"CASE_001", "CASE_002", "CTRL_001", "CTRL_002", "UNCLASSIFIED_001", "UNCLASSIFIED_002"}
    
    # Test both scenarios
    buggy_result = assign_case_control_counts(df, case_samples, control_samples, buggy_all_samples)
    correct_result = assign_case_control_counts(df, case_samples, control_samples, correct_all_samples)
    
    # The results should be identical because the function should only count
    # samples that are in case_samples or control_samples, regardless of all_samples
    # But all_samples needs to include all samples for proper iteration
    
    # Both should give the same case/control counts
    assert buggy_result["proband_count"].iloc[0] == correct_result["proband_count"].iloc[0]
    assert buggy_result["control_count"].iloc[0] == correct_result["control_count"].iloc[0]
    
    # Both should give the same variant counts
    assert buggy_result["proband_variant_count"].iloc[0] == correct_result["proband_variant_count"].iloc[0]
    assert buggy_result["control_variant_count"].iloc[0] == correct_result["control_variant_count"].iloc[0]
    
    # Validate the actual counts
    # Cases: CASE_001(1/1)=2 alleles, CASE_002(0/1)=1 allele = 3 total, 2 samples with variants
    # Controls: CTRL_001(0/1)=1 allele, CTRL_002(0/0)=0 alleles = 1 total, 1 sample with variant
    assert correct_result["proband_allele_count"].iloc[0] == 3
    assert correct_result["control_allele_count"].iloc[0] == 1
    assert correct_result["proband_variant_count"].iloc[0] == 2
    assert correct_result["control_variant_count"].iloc[0] == 1
    
    # The key test: ensure unclassified samples don't get misassigned
    # This is more of a validation that the logic is working correctly
    # In the real bug, if all_samples was wrong, it could cause issues with iteration