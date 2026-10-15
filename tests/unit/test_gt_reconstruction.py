"""
Tests for GT column reconstruction from per-sample columns (Phase 11).

These tests prove that reconstruct_gt_column() produces identical output format
to the old genotype replacement pipeline. Per-sample columns use the naming
convention from bcftools query: GEN[0].GT, GEN[1].GT, etc. (or sanitized
GEN_0__GT, GEN_1__GT after DataFrame loading).
"""

import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import (
    _find_per_sample_gt_columns,
    reconstruct_gt_column,
)


def _make_gt_col(idx: int, sanitized: bool = True) -> str:
    """Helper to create GT column name."""
    return f"GEN_{idx}__GT" if sanitized else f"GEN[{idx}].GT"


@pytest.mark.unit
def test_find_per_sample_gt_columns_sanitized():
    """Test detection of sanitized GT column names."""
    df = pd.DataFrame({"CHROM": [], "GEN_0__GT": [], "GEN_1__GT": [], "GEN_2__GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT"]


@pytest.mark.unit
def test_find_per_sample_gt_columns_original():
    """Test detection of original (unsanitized) GT column names."""
    df = pd.DataFrame({"CHROM": [], "GEN[0].GT": [], "GEN[1].GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == ["GEN[0].GT", "GEN[1].GT"]


@pytest.mark.unit
def test_find_per_sample_gt_columns_sorted():
    """Test that columns are sorted by index, not lexicographically."""
    df = pd.DataFrame({"GEN_10__GT": [], "GEN_2__GT": [], "GEN_0__GT": [], "GEN_1__GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_10__GT"]


@pytest.mark.unit
def test_find_per_sample_gt_columns_no_match():
    """Test that non-GT columns are not matched."""
    df = pd.DataFrame({"CHROM": [], "POS": [], "GENE": [], "GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == []


@pytest.mark.unit
def test_basic_reconstruction():
    """Test basic GT reconstruction with mixed genotypes."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["0/1", "1/1"],
            "GEN_1__GT": ["0/0", "0/1"],
            "GEN_2__GT": ["1/1", "./."],
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Row 0: Sample1(0/1) and Sample3(1/1) - skip Sample2 (0/0)
    # Row 1: Sample1(1/1) and Sample2(0/1) - skip Sample3 (./.)
    assert "GT" in result.columns
    assert result["GT"].iloc[0] == "Sample1(0/1);Sample3(1/1)"
    assert result["GT"].iloc[1] == "Sample1(1/1);Sample2(0/1)"

    # Per-sample columns should be dropped
    assert "GEN_0__GT" not in result.columns
    assert "GEN_1__GT" not in result.columns
    assert "GEN_2__GT" not in result.columns


@pytest.mark.unit
def test_all_reference_genotypes():
    """Test GT reconstruction when all samples are reference."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["0/0", "./."],
            "GEN_1__GT": ["0/0", "0/0"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "GT" in result.columns
    assert result["GT"].iloc[0] == ""
    assert result["GT"].iloc[1] == ""


@pytest.mark.unit
def test_missing_values():
    """Test GT reconstruction with NA/NaN/empty values."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [111, 222, 333],
            "GEN_0__GT": ["0/1", pd.NA, ""],
            "GEN_1__GT": [None, "1/1", "."],
            "GEN_2__GT": ["1/0", "0/0", "0/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "Sample1(0/1);Sample3(1/0)"
    assert result["GT"].iloc[1] == "Sample2(1/1)"
    assert result["GT"].iloc[2] == "Sample3(0/1)"


@pytest.mark.unit
def test_single_sample():
    """Test GT reconstruction with only one sample."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["0/1", "1/1"],
        }
    )
    vcf_samples = ["SampleA"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "SampleA(0/1)"
    assert result["GT"].iloc[1] == "SampleA(1/1)"


@pytest.mark.unit
def test_many_samples():
    """Test GT reconstruction with 12 samples."""
    samples = [f"Sample{i:02d}" for i in range(1, 13)]
    data = {
        "CHROM": ["chr1"],
        "POS": [12345],
    }
    # Samples 1, 3, 5, 7, 9, 11 have variants; others are 0/0 or ./.
    gts = ["0/1", "0/0", "1/1", "0/0", "0/1", "./.", "1/0", "0/0", "1/1", "0/0", "0/1", "0/0"]
    for i, gt in enumerate(gts):
        data[f"GEN_{i}__GT"] = [gt]

    df = pd.DataFrame(data)
    result = reconstruct_gt_column(df, samples)

    expected = "Sample01(0/1);Sample03(1/1);Sample05(0/1);Sample07(1/0);Sample09(1/1);Sample11(0/1)"
    assert result["GT"].iloc[0] == expected


@pytest.mark.unit
def test_phased_genotypes():
    """Test GT reconstruction with phased genotypes (|)."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [12345, 67890, 11111],
            "GEN_0__GT": ["0|1", "1|0", "0|0"],
            "GEN_1__GT": ["0|0", "1|1", ".|."],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "Sample1(0|1)"
    assert result["GT"].iloc[1] == "Sample1(1|0);Sample2(1|1)"
    assert result["GT"].iloc[2] == ""


@pytest.mark.unit
def test_hemizygous_genotypes():
    """Test GT reconstruction with hemizygous genotypes (single allele)."""
    df = pd.DataFrame(
        {
            "CHROM": ["chrX", "chrX"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["1", "0"],
            "GEN_1__GT": ["0/1", "0/0"],
        }
    )
    vcf_samples = ["MaleSample", "FemaleSample"]

    result = reconstruct_gt_column(df, vcf_samples)

    # "0" is not in skip list, so hemizygous ref is included
    assert result["GT"].iloc[0] == "MaleSample(1);FemaleSample(0/1)"
    assert result["GT"].iloc[1] == "MaleSample(0)"


@pytest.mark.unit
def test_per_sample_columns_dropped():
    """Test that per-sample columns are removed after reconstruction."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "REF": ["A"],
            "ALT": ["G"],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": ["1/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "CHROM" in result.columns
    assert "POS" in result.columns
    assert "REF" in result.columns
    assert "ALT" in result.columns
    assert "GT" in result.columns
    assert "GEN_0__GT" not in result.columns
    assert "GEN_1__GT" not in result.columns


@pytest.mark.unit
def test_column_ordering():
    """Test that GT column is added at the end."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "REF": ["A"],
            "ALT": ["G"],
            "GEN_0__GT": ["0/1"],
        }
    )
    vcf_samples = ["Sample1"]

    result = reconstruct_gt_column(df, vcf_samples)

    columns = result.columns.tolist()
    assert columns == ["CHROM", "POS", "REF", "ALT", "GT"]


@pytest.mark.unit
def test_empty_dataframe():
    """Test GT reconstruction with zero rows."""
    df = pd.DataFrame(
        {
            "CHROM": pd.Series(dtype=str),
            "POS": pd.Series(dtype=int),
            "GEN_0__GT": pd.Series(dtype=str),
            "GEN_1__GT": pd.Series(dtype=str),
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "GT" in result.columns
    assert len(result) == 0
    assert "GEN_0__GT" not in result.columns
    assert "GEN_1__GT" not in result.columns


@pytest.mark.unit
def test_original_column_names():
    """Test GT reconstruction with unsanitized GEN[N].GT column names."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "GEN[0].GT": ["0/1"],
            "GEN[1].GT": ["1/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "Sample1(0/1);Sample2(1/1)"
    assert "GEN[0].GT" not in result.columns
    assert "GEN[1].GT" not in result.columns


@pytest.mark.unit
def test_sample_count_mismatch_warning(caplog):
    """Test warning when GT column count differs from sample count."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": ["1/1"],
        }
    )
    # More samples than GT columns
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "GT column count (2) != sample count (3)" in caplog.text
    # Should still work with available columns
    assert result["GT"].iloc[0] == "Sample1(0/1);Sample2(1/1)"


@pytest.mark.unit
def test_reconstruct_gt_column_with_extra_sample_fields_retains_normal():
    """reconstruct_gt_column appends extra_sample_fields (GT:DP:AD:AF) and retains 0/0 normal."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "GEN_0__GT": ["0/0"],
            "GEN_0__DP": ["42"],
            "GEN_0__AD": ["42,0"],
            "GEN_0__AF": ["0.00"],
            "GEN_1__GT": ["0/1"],
            "GEN_1__DP": ["85"],
            "GEN_1__AD": ["60,25"],
            "GEN_1__AF": ["0.29"],
        }
    )
    res = reconstruct_gt_column(
        df,
        ["NORMAL", "TUMOR"],
        extra_sample_fields=["GEN[*].DP", "GEN[*].AD", "GEN[*].AF"],
    )
    assert res["GT"].iloc[0] == "NORMAL(0/0:42:42,0:0.00);TUMOR(0/1:85:60,25:0.29)"
    assert "GEN_0__DP" not in res.columns
    assert "GEN_1__AF" not in res.columns


@pytest.mark.unit
def test_cli_preserves_config_json_options_and_final_filtering_reads_filters(tmp_path, monkeypatch):
    """CLI -c config.json preserves options and TSVOutputStage packs fields."""
    import json
    import logging
    from pathlib import Path
    from unittest.mock import Mock

    from variantcentrifuge import cli
    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace
    from variantcentrifuge.stages.output_stages import FinalFilteringStage, TSVOutputStage

    vcf = tmp_path / "test.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n"
    )
    tx_file = tmp_path / "tx.txt"
    tx_file.write_text("ENST00000269305.4\n")
    cfg_file = tmp_path / "custom_config.json"
    base_cfg = json.loads((Path(cli.__file__).parent / "config.json").read_text(encoding="utf-8"))
    base_cfg.update(
        {
            "filters": "(FILTER = 'PASS')",
            "transcript_file": str(tx_file),
            "final_filter": "POS > 50",
            "late_filtering": True,
            "extra_sample_fields": ["GEN[*].DP", "GEN[*].AD", "GEN[*].AF"],
            "extra_sample_field_delimiter": "|",
        }
    )
    cfg_file.write_text(json.dumps(base_cfg), encoding="utf-8")

    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda a: captured.update(a.config))
    monkeypatch.setattr(
        "sys.argv",
        ["variantcentrifuge", "-v", str(vcf), "-g", "TP53", "-c", str(cfg_file)],
    )
    prev_level = logging.getLogger("variantcentrifuge").level
    try:
        cli.main()
    finally:
        logging.getLogger("variantcentrifuge").setLevel(prev_level)

    assert captured["transcript_file"] == str(tx_file)
    assert captured["final_filter"] == "POS > 50"
    assert captured["late_filtering"] is True
    assert captured["append_extra_sample_fields"] is True
    assert captured["extra_sample_fields"] == ["GEN[*].DP", "GEN[*].AD", "GEN[*].AF"]
    assert captured["extra_sample_field_delimiter"] == "|"

    stage = FinalFilteringStage()
    ctx = Mock(spec=PipelineContext)
    ctx.config = {"late_filtering": True, "filters": "POS > 100"}
    ctx.current_dataframe = pd.DataFrame({"CHROM": ["chr1", "chr1"], "POS": [50, 200]})
    out_ctx = stage._process(ctx)
    assert len(out_ctx.current_dataframe) == 1
    assert out_ctx.current_dataframe["POS"].iloc[0] == 200

    # End-to-end test through TSVOutputStage:
    # 1. With append_extra_sample_fields=False (default), 0/0 is excluded and GEN[i].X dropped
    raw_df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [200],
            "REF": ["C"],
            "ALT": ["T"],
            "GEN_0__GT": ["0/0"],
            "GEN_0__DP": ["42"],
            "GEN_0__AD": ["42,0"],
            "GEN_0__AF": ["0.00"],
            "GEN_1__GT": ["0/1"],
            "GEN_1__DP": ["85"],
            "GEN_1__AD": ["60,25"],
            "GEN_1__AF": ["0.29"],
        }
    )
    ws = Workspace(output_dir=tmp_path, base_name="test")
    tsv_ctx_default = PipelineContext(
        args=Mock(output_file=str(tmp_path / "default.tsv"), output_dir=str(tmp_path)),
        config={"output_file": str(tmp_path / "default.tsv"), "no_links": True},
        workspace=ws,
    )
    tsv_ctx_default.current_dataframe = raw_df.copy()
    tsv_ctx_default.vcf_samples = ["NORMAL", "TUMOR"]
    TSVOutputStage()._process(tsv_ctx_default)
    default_out = pd.read_csv(tmp_path / "default.tsv", sep="\t")
    assert default_out["GT"].iloc[0] == "TUMOR(0/1)"
    assert not any(c.startswith("GEN") for c in default_out.columns)

    # 2. With append_extra_sample_fields=True through VariantAnalysisStage -> TSVOutputStage,
    # configured fields and 0/0 reference normal are preserved through analysis and packed at output
    from variantcentrifuge.stages.analysis_stages import VariantAnalysisStage

    analysis_df = raw_df.copy()
    analysis_df["POS"] = analysis_df["POS"].astype(str)
    analysis_df["GENE"] = ["KCNJ5"]
    tsv_ctx_extra = PipelineContext(
        args=Mock(output_file=str(tmp_path / "extra.tsv"), output_dir=str(tmp_path)),
        config={
            "output_file": str(tmp_path / "extra.tsv"),
            "no_links": True,
            "no_stats": True,
            "gzip_intermediates": False,
            "append_extra_sample_fields": True,
            "extra_sample_fields": ["GEN[*].DP", "GEN[*].AD", "GEN[*].AF"],
            "extra_sample_field_delimiter": "|",
        },
        workspace=ws,
    )
    tsv_ctx_extra.current_dataframe = analysis_df
    tsv_ctx_extra.vcf_samples = ["NORMAL", "TUMOR"]
    tsv_ctx_extra.column_rename_map = {
        "GEN[0].GT": "GEN_0__GT",
        "GEN[0].DP": "GEN_0__DP",
        "GEN[0].AD": "GEN_0__AD",
        "GEN[0].AF": "GEN_0__AF",
        "GEN[1].GT": "GEN_1__GT",
        "GEN[1].DP": "GEN_1__DP",
        "GEN[1].AD": "GEN_1__AD",
        "GEN[1].AF": "GEN_1__AF",
    }
    tsv_ctx_extra = VariantAnalysisStage()._process(tsv_ctx_extra)
    TSVOutputStage()._process(tsv_ctx_extra)
    extra_out = pd.read_csv(tmp_path / "extra.tsv", sep="\t")
    assert extra_out["GT"].iloc[0] == "NORMAL(0/0|42|42,0|0.00);TUMOR(0/1|85|60,25|0.29)"
    assert not any(c.startswith(("GEN_", "GEN[")) for c in extra_out.columns)
