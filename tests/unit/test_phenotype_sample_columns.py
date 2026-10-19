"""
Tests for phenotype extraction from per-sample columns (Phase 11).

These tests prove that extract_phenotypes_from_sample_columns() produces
equivalent output to the old extract_phenotypes_for_gt_row() function.
"""

from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.phenotype import (
    extract_phenotypes_for_gt_row,
    extract_phenotypes_from_sample_columns,
)


@pytest.mark.unit
def test_basic_extraction():
    """Test basic phenotype extraction from per-sample columns."""
    row = pd.Series({"Sample1": "0/1", "Sample2": "0/0", "Sample3": "1/1"})
    vcf_samples = ["Sample1", "Sample2", "Sample3"]
    phenotypes = {"Sample1": {"HPO:001"}, "Sample3": {"HPO:002"}}

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)

    # Sample1 has variant and phenotype HPO:001
    # Sample2 has no variant (0/0) -> excluded
    # Sample3 has variant and phenotype HPO:002
    assert result == "Sample1(HPO:001);Sample3(HPO:002)"


@pytest.mark.unit
def test_multiple_phenotypes_per_sample():
    """Test extraction with multiple phenotypes per sample."""
    row = pd.Series({"Sample1": "0/1", "Sample2": "1/1"})
    vcf_samples = ["Sample1", "Sample2"]
    phenotypes = {
        "Sample1": {"HPO:001", "HPO:002", "HPO:003"},
        "Sample2": {"HPO:010"},
    }

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)

    # Multiple phenotypes should be sorted and comma-separated
    assert result == "Sample1(HPO:001,HPO:002,HPO:003);Sample2(HPO:010)"


@pytest.mark.unit
def test_multiple_samples_with_variants():
    """Test extraction with multiple samples having variants."""
    row = pd.Series(
        {
            "Sample1": "0/1",
            "Sample2": "0/0",
            "Sample3": "1/1",
            "Sample4": "0/1",
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3", "Sample4"]
    phenotypes = {
        "Sample1": {"HPO:A"},
        "Sample3": {"HPO:B", "HPO:C"},
        "Sample4": {"HPO:D"},
    }

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)

    # Sample2 excluded (0/0), others included
    assert result == "Sample1(HPO:A);Sample3(HPO:B,HPO:C);Sample4(HPO:D)"


@pytest.mark.unit
def test_sample_with_variant_but_no_phenotypes():
    """Test extraction when sample has variant but no phenotypes."""
    row = pd.Series({"Sample1": "0/1", "Sample2": "1/1"})
    vcf_samples = ["Sample1", "Sample2"]
    phenotypes = {"Sample1": {"HPO:001"}}  # Sample2 not in phenotypes dict

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)

    # Sample2 has variant but no phenotypes -> empty parentheses
    assert result == "Sample1(HPO:001);Sample2()"


@pytest.mark.unit
def test_all_reference():
    """Test extraction when all samples are reference."""
    row = pd.Series({"Sample1": "0/0", "Sample2": "./.", "Sample3": "0/0"})
    vcf_samples = ["Sample1", "Sample2", "Sample3"]
    phenotypes = {"Sample1": {"HPO:001"}, "Sample2": {"HPO:002"}}

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)

    # All reference -> empty string
    assert result == ""


@pytest.mark.unit
def test_missing_na_values():
    """Test extraction with missing/NA values."""
    row = pd.Series(
        {
            "Sample1": "0/1",
            "Sample2": pd.NA,
            "Sample3": "",
            "Sample4": ".",
            "Sample5": "1/1",
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]
    phenotypes = {
        "Sample1": {"HPO:A"},
        "Sample2": {"HPO:B"},
        "Sample3": {"HPO:C"},
        "Sample4": {"HPO:D"},
        "Sample5": {"HPO:E"},
    }

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)

    # Only Sample1 and Sample5 have valid variants
    assert result == "Sample1(HPO:A);Sample5(HPO:E)"


@pytest.mark.unit
def test_equivalence_with_old_function():
    """
    Critical test: Verify equivalence between old and new phenotype extraction.

    This proves that extract_phenotypes_from_sample_columns() produces the same
    output as extract_phenotypes_for_gt_row() for the same variant data.
    """
    # Define test data
    vcf_samples = ["Sample1", "Sample2", "Sample3"]
    phenotypes = {
        "Sample1": {"HPO:001", "HPO:002"},
        "Sample2": {"HPO:010"},
        # Sample3 has no phenotypes
    }

    # Test case 1: Mixed genotypes
    # Old format: packed GT string
    gt_old = "Sample1(0/1);Sample3(1/1)"
    result_old = extract_phenotypes_for_gt_row(gt_old, phenotypes)

    # New format: per-sample columns
    row_new = pd.Series({"Sample1": "0/1", "Sample2": "0/0", "Sample3": "1/1"})
    result_new = extract_phenotypes_from_sample_columns(row_new, vcf_samples, phenotypes)

    # Should produce identical output
    assert result_old == result_new == "Sample1(HPO:001,HPO:002);Sample3()"

    # Test case 2: All reference
    gt_old2 = ""
    result_old2 = extract_phenotypes_for_gt_row(gt_old2, phenotypes)

    row_new2 = pd.Series({"Sample1": "0/0", "Sample2": "0/0", "Sample3": "./."})
    result_new2 = extract_phenotypes_from_sample_columns(row_new2, vcf_samples, phenotypes)

    assert result_old2 == result_new2 == ""

    # Test case 3: Single sample with phenotypes
    gt_old3 = "Sample2(0/1)"
    result_old3 = extract_phenotypes_for_gt_row(gt_old3, phenotypes)

    row_new3 = pd.Series({"Sample1": "0/0", "Sample2": "0/1", "Sample3": "0/0"})
    result_new3 = extract_phenotypes_from_sample_columns(row_new3, vcf_samples, phenotypes)

    assert result_old3 == result_new3 == "Sample2(HPO:010)"


@pytest.mark.unit
def test_namedtuple_access():
    """Test that function works with namedtuple (from itertuples)."""
    # Create a DataFrame and convert row to namedtuple
    df = pd.DataFrame({"Sample1": ["0/1"], "Sample2": ["1/1"]})
    row = next(df.itertuples(index=False))  # namedtuple

    vcf_samples = ["Sample1", "Sample2"]
    phenotypes = {"Sample1": {"HPO:A"}, "Sample2": {"HPO:B"}}

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)

    assert result == "Sample1(HPO:A);Sample2(HPO:B)"


@pytest.mark.unit
def test_phased_genotypes():
    """Test extraction with phased genotypes (|)."""
    row = pd.Series({"Sample1": "0|1", "Sample2": "1|0", "Sample3": "0|0"})
    vcf_samples = ["Sample1", "Sample2", "Sample3"]
    phenotypes = {"Sample1": {"HPO:X"}, "Sample2": {"HPO:Y"}}

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    assert result == "Sample1(HPO:X);Sample2(HPO:Y)"


@pytest.mark.unit
def test_indexed_genotype_columns():
    """Test extraction from standard bcftools GEN[i].GT indexed columns."""
    row = pd.Series({"GEN[0].GT": "0/1", "GEN[1].GT": "0/0", "GEN[2].GT": "1/1"})
    vcf_samples = ["SampleA", "SampleB", "SampleC"]
    phenotypes = {"SampleA": {"HPO:001"}, "SampleC": {"HPO:002"}}

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    assert result == "SampleA(HPO:001);SampleC(HPO:002)"


@pytest.mark.unit
def test_reversed_and_double_digit_indexed_columns():
    """Test extraction with reversed column order and indices above 9."""
    samples = [f"Sample{i}" for i in range(12)]
    # Reversed dictionary order
    row_data = {f"GEN[{i}].GT": "0/0" for i in range(11, -1, -1)}
    row_data["GEN[10].GT"] = "0/1"
    row_data["GEN[0].GT"] = "1/1"
    row = pd.Series(row_data)

    phenotypes = {"Sample0": {"HPO:ZERO"}, "Sample10": {"HPO:TEN"}}

    result = extract_phenotypes_from_sample_columns(row, samples, phenotypes)
    assert result == "Sample0(HPO:ZERO);Sample10(HPO:TEN)"


@pytest.mark.unit
def test_reversed_indexed_columns_order():
    """Test extraction when GEN[i].GT columns appear in reversed order (Series and namedtuple)."""
    vcf_samples = ["SampleA", "SampleB", "SampleC"]
    phenotypes = {"SampleA": {"HPO:001"}, "SampleC": {"HPO:002"}}

    # Series with reversed raw bracket column order (as passed by df.apply)
    row_raw = pd.Series({"GEN[2].GT": "1/1", "GEN[1].GT": "0/0", "GEN[0].GT": "0/1"})
    result_raw = extract_phenotypes_from_sample_columns(row_raw, vcf_samples, phenotypes)
    assert result_raw == "SampleA(HPO:001);SampleC(HPO:002)"

    # Series with reversed sanitized column order
    row_san = pd.Series({"GEN_2__GT": "1/1", "GEN_1__GT": "0/0", "GEN_0__GT": "0/1"})
    result_san = extract_phenotypes_from_sample_columns(row_san, vcf_samples, phenotypes)
    assert result_san == "SampleA(HPO:001);SampleC(HPO:002)"

    # Namedtuple with reversed sanitized column order (from itertuples)
    df = pd.DataFrame([{"GEN_2__GT": "1/1", "GEN_1__GT": "0/0", "GEN_0__GT": "0/1"}])
    nt_row = next(df.itertuples(index=False))
    nt_result = extract_phenotypes_from_sample_columns(nt_row, vcf_samples, phenotypes)
    assert nt_result == "SampleA(HPO:001);SampleC(HPO:002)"


@pytest.mark.unit
def test_sanitized_indexed_columns():
    """Test extraction from sanitized GEN_i__GT columns."""
    row = pd.Series({"GEN_0__GT": "0/1", "GEN_1__GT": "0/0"})
    vcf_samples = ["SampleX", "SampleY"]
    phenotypes = {"SampleX": {"HPO:X"}}

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    assert result == "SampleX(HPO:X)"


@pytest.mark.unit
def test_phenotype_extraction_colliding_sample_names_and_tuple_methods():
    """Verify sample names matching columns (CHROM) or tuple methods (count) are handled safely."""
    vcf_samples = ["CHROM", "count"]
    phenotypes = {"CHROM": {"HPO:CHR"}, "count": {"HPO:CNT"}}

    # Test 1: Series with CHROM column and indexed GEN columns
    row_series = pd.Series({"CHROM": "chr1", "GEN[0].GT": "0/0", "GEN[1].GT": "0/1"})
    res_series = extract_phenotypes_from_sample_columns(row_series, vcf_samples, phenotypes)
    assert res_series == "count(HPO:CNT)"

    # Test 2: namedtuple from DataFrame with CHROM and sanitized GEN columns
    df = pd.DataFrame([{"CHROM": "chr1", "GEN_0__GT": "0/0", "GEN_1__GT": "0/1"}])
    nt_row = next(df.itertuples(index=False))
    res_nt = extract_phenotypes_from_sample_columns(nt_row, vcf_samples, phenotypes)
    assert res_nt == "count(HPO:CNT)"

    # Test 3: namedtuple where sample "count" is NOT in DataFrame and has no GEN column
    df_no_count = pd.DataFrame([{"CHROM": "chr1", "GEN_0__GT": "0/0"}])
    nt_no_count = next(df_no_count.itertuples(index=False))
    res_no_count = extract_phenotypes_from_sample_columns(nt_no_count, vcf_samples, phenotypes)
    assert res_no_count == ""


@pytest.mark.unit
def test_series_with_fields_column_name_not_treated_as_namedtuple():
    """Verify pd.Series with a sample column named '_fields' is not mistaken for a namedtuple."""
    vcf_samples = ["_fields", "S2"]
    phenotypes = {"_fields": {"H1"}, "S2": {"H2"}}
    row = pd.Series({"_fields": "0/1", "S2": "1/1"})

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    assert result == "_fields(H1);S2(H2)"


@pytest.mark.unit
@pytest.mark.parametrize(
    "gt,is_carrier",
    [
        ("0|0", False),
        ("0/0", False),
        (".|.", False),
        ("./.", False),
        ("0/.", False),
        ("./1", False),
        ("0", False),
        (".", False),
        ("NA", False),
        ("", False),
        ("0/0:30", False),
        ("0|1", True),
        ("1|0", True),
        ("1/1", True),
        ("0/1", True),
        ("1", True),
        ("0/1:20", True),
        ("1/2", True),
    ],
)
def test_genotype_allele_calling_and_exclusion(gt: str, is_carrier: bool) -> None:
    """Verify carrier detection and exclusion for phased, missing, and haploid GTs."""
    vcf_samples = ["S1"]
    phenotypes = {"S1": {"HPO:1"}}
    row = pd.Series({"S1": gt})

    result_cols = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    result_old = extract_phenotypes_for_gt_row(f"S1({gt})", phenotypes)

    expected = "S1(HPO:1)" if is_carrier else ""
    assert result_cols == expected
    assert result_old == expected


@pytest.mark.unit
def test_sample_named_column_colliding_with_indexed_gen_bracket(tmp_path: Path) -> None:
    """Verify sample named GEN[1].GT does not cause collision for subsequent sample S2."""
    from argparse import Namespace
    from unittest.mock import Mock

    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace
    from variantcentrifuge.stages.post_processing_stages import PhenotypeIntegrationStage

    vcf_samples = ["GEN[1].GT", "S2"]
    phenotypes = {"GEN[1].GT": {"H1"}, "S2": {"H2"}}
    row = pd.Series({"GEN[1].GT": "0/1", "S2": "0/0"})

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    assert result == "GEN[1].GT(H1)"

    tsv_file = tmp_path / "test_bracket.tsv"
    tsv_file.write_text("CHROM\tPOS\tGEN[1].GT\tS2\nchr1\t100\t0/1\t0/0\n")

    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = tmp_path
    ws.get_intermediate_path = lambda name: tmp_path / name

    ctx = PipelineContext(args=Namespace(no_replacement=True), config={}, workspace=ws)
    ctx.extracted_tsv = tsv_file
    ctx.data = tsv_file
    ctx.vcf_samples = vcf_samples
    ctx.phenotype_data = phenotypes

    stage = PhenotypeIntegrationStage()
    out_ctx = stage._process(ctx)
    df_out = pd.read_csv(out_ctx.phenotypes_added_tsv, sep="\t")
    assert df_out["Phenotypes"].iloc[0] == "GEN[1].GT(H1)"


@pytest.mark.unit
def test_sample_named_column_colliding_with_indexed_gen_underscore(tmp_path: Path) -> None:
    """Verify sample named GEN_1__GT does not cause collision for subsequent sample S2."""
    from argparse import Namespace
    from unittest.mock import Mock

    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace
    from variantcentrifuge.stages.post_processing_stages import PhenotypeIntegrationStage

    vcf_samples = ["GEN_1__GT", "S2"]
    phenotypes = {"GEN_1__GT": {"H1"}, "S2": {"H2"}}
    row = pd.Series({"GEN_1__GT": "0/1", "S2": "0/0"})

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    assert result == "GEN_1__GT(H1)"

    tsv_file = tmp_path / "test_underscore.tsv"
    tsv_file.write_text("CHROM\tPOS\tGEN_1__GT\tS2\nchr1\t100\t0/1\t0/0\n")

    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = tmp_path
    ws.get_intermediate_path = lambda name: tmp_path / name

    ctx = PipelineContext(args=Namespace(no_replacement=True), config={}, workspace=ws)
    ctx.extracted_tsv = tsv_file
    ctx.data = tsv_file
    ctx.vcf_samples = vcf_samples
    ctx.phenotype_data = phenotypes

    stage = PhenotypeIntegrationStage()
    out_ctx = stage._process(ctx)
    df_out = pd.read_csv(out_ctx.phenotypes_added_tsv, sep="\t")
    assert df_out["Phenotypes"].iloc[0] == "GEN_1__GT(H1)"


@pytest.mark.unit
def test_sample_named_columns_matching_indexed_pattern_precedence():
    """Verify that when vcf_samples = ['GEN[1].GT', 'GEN[0].GT'] and columns are indexed
    GEN[0].GT, GEN[1].GT, indexed matching takes precedence so sample 0 ('GEN[1].GT')
    maps to GEN[0].GT ('0/1') and gets annotated.
    """
    vcf_samples = ["GEN[1].GT", "GEN[0].GT"]
    phenotypes = {"GEN[1].GT": {"H1"}, "GEN[0].GT": {"H0"}}
    row = pd.Series({"GEN[0].GT": "0/1", "GEN[1].GT": "0/0"})

    result = extract_phenotypes_from_sample_columns(row, vcf_samples, phenotypes)
    assert result == "GEN[1].GT(H1)"
