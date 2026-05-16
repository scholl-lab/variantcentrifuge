"""
Simple integration test for --final-filter functionality.

This test focuses on verifying the filter_dataframe_with_query function
works correctly when integrated into the pipeline flow.
"""

import pandas as pd
import pytest

from variantcentrifuge.filters import filter_dataframe_with_query, filter_tsv_with_expression


def test_final_filter_integration():
    """Test that final filter correctly filters a DataFrame like the pipeline would produce."""
    # Create a sample DataFrame that mimics what the pipeline would produce
    df = pd.DataFrame(
        {
            "VAR_ID": ["var_0001_abc1", "var_0002_def2", "var_0003_ghi3", "var_0004_jkl4"],
            "CHROM": ["chr1", "chr1", "chr2", "chr3"],
            "POS": ["1000", "2000", "3000", "4000"],
            "REF": ["A", "G", "C", "T"],
            "ALT": ["T", "C", "G", "A"],
            "GENE": ["BRCA1", "BRCA1", "TP53", "EGFR"],
            "IMPACT": ["HIGH", "MODERATE", "HIGH", "LOW"],
            "inheritance_score": ["0.95", "0.4", "0.8", "0.2"],
            "CADD_phred": ["25.5", "15.0", "30.0", "10.0"],
            "Inheritance_Pattern": ["de_novo", "inherited", "compound_heterozygous", "unknown"],
        }
    )

    # Test 1: Filter by numeric score
    result1 = filter_dataframe_with_query(df, "inheritance_score > 0.5")
    assert len(result1) == 2
    assert set(result1["GENE"]) == {"BRCA1", "TP53"}

    # Test 2: Filter by string field
    result2 = filter_dataframe_with_query(df, 'IMPACT == "HIGH"')
    assert len(result2) == 2
    assert set(result2["CHROM"]) == {"chr1", "chr2"}

    # Test 3: Complex filter with AND
    result3 = filter_dataframe_with_query(df, 'IMPACT == "HIGH" and inheritance_score > 0.9')
    assert len(result3) == 1
    assert result3.iloc[0]["GENE"] == "BRCA1"

    # Test 4: Filter using IN operator
    result4 = filter_dataframe_with_query(
        df, 'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]'
    )
    assert len(result4) == 2
    assert set(result4["GENE"]) == {"BRCA1", "TP53"}

    # Test 5: Filter with OR condition
    result5 = filter_dataframe_with_query(df, "inheritance_score > 0.9 or CADD_phred > 28")
    assert len(result5) == 2
    # POS might be converted to numeric, so convert to string for comparison
    assert set(result5["POS"].astype(str)) == {"1000", "3000"}

    # Test 6: Verify that VAR_ID column is preserved
    result6 = filter_dataframe_with_query(df, 'GENE == "BRCA1"')
    assert "VAR_ID" in result6.columns
    assert len(result6) == 2
    assert all(vid.startswith("var_") for vid in result6["VAR_ID"])


def test_final_filter_with_empty_result():
    """Test that filter handles cases where no rows match."""
    df = pd.DataFrame({"CHROM": ["chr1", "chr2"], "POS": ["1000", "2000"], "score": ["0.1", "0.2"]})

    # Filter that matches nothing
    result = filter_dataframe_with_query(df, "score > 0.5")
    assert len(result) == 0
    assert list(result.columns) == list(df.columns)


def test_final_filter_preserves_all_columns():
    """Test that filtering preserves all columns from the original DataFrame."""
    # Create a DataFrame with many columns like the real pipeline output
    df = pd.DataFrame(
        {
            "VAR_ID": ["var_0001_xyz"],
            "CHROM": ["chr1"],
            "POS": ["1000"],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["BRCA1"],
            "IMPACT": ["HIGH"],
            "FILTER": ["PASS"],
            "QUAL": ["100"],
            "Custom_Annotation": ["InGeneList=cancer_genes"],
            "inheritance_score": ["0.95"],
            "Inheritance_Pattern": ["de_novo"],
            "Inheritance_Confidence": ["0.9"],
            "Link_gnomAD": ["https://gnomad.broadinstitute.org/..."],
            "Link_ClinVar": ["https://www.ncbi.nlm.nih.gov/clinvar/..."],
        }
    )

    # Apply a filter
    result = filter_dataframe_with_query(df, "inheritance_score > 0.5")

    # Check all columns are preserved
    assert list(result.columns) == list(df.columns)
    assert len(result) == 1


def test_final_filter_invalid_expression_raises():
    df = pd.DataFrame({"CHROM": ["chr1"], "score": ["0.1"]})

    with pytest.raises(ValueError, match="Invalid final filter expression"):
        filter_dataframe_with_query(df, 'CHROM === "chr1"')


def test_tsv_filter_invalid_expression_raises(tmp_path):
    input_tsv = tmp_path / "input.tsv"
    output_tsv = tmp_path / "output.tsv"
    input_tsv.write_text("CHROM\tscore\nchr1\t0.1\n")

    with pytest.raises(ValueError, match="Invalid TSV filter expression"):
        filter_tsv_with_expression(str(input_tsv), str(output_tsv), 'CHROM === "chr1"')

    assert not output_tsv.exists()


if __name__ == "__main__":
    test_final_filter_integration()
    test_final_filter_with_empty_result()
    test_final_filter_preserves_all_columns()
    print("All tests passed!")
