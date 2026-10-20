"""Unit tests for preserving per-sample genotypes without multiplying variant rows (#152)."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.analysis_stages import VariantAnalysisStage


def _create_context(tmp_path: Path, df: pd.DataFrame, samples: list[str]) -> PipelineContext:
    ws = Workspace(tmp_path, "test_ws")
    context = PipelineContext(
        args=Namespace(
            output_file=str(tmp_path / "out.tsv"),
            output_dir=str(tmp_path),
        ),
        config={
            "output_file": str(tmp_path / "out.tsv"),
            "no_links": True,
            "no_stats": True,
            "gzip_intermediates": False,
            "case_samples": [samples[0]],
            "control_samples": [samples[1]] if len(samples) > 1 else [],
        },
        workspace=ws,
    )
    context.current_dataframe = df
    context.vcf_samples = samples
    return context


@pytest.mark.parametrize("n_rows", [2, 3, 5])
def test_duplicate_keys_row_count_preservation(tmp_path: Path, n_rows: int):
    """Test that 2, 3, and 5 duplicate-key input rows produce exactly 2, 3, and 5 output rows."""
    samples = ["SAMPLE1", "SAMPLE2"]
    # Identical genomic keys across all rows
    data = {
        "CHROM": ["chr1"] * n_rows,
        "POS": ["1000"] * n_rows,
        "REF": ["A"] * n_rows,
        "ALT": ["T"] * n_rows,
        "GENE": ["BRCA1"] * n_rows,
        "annotation_tag": [f"tag_{i}" for i in range(n_rows)],
        "GEN[0].GT": ["0/1" if i % 2 == 0 else "1/1" for i in range(n_rows)],
        "GEN[1].GT": ["0/0" for _ in range(n_rows)],
    }
    input_df = pd.DataFrame(data)
    input_copy = input_df.copy(deep=True)

    ctx = _create_context(tmp_path, input_df, samples)
    processed = VariantAnalysisStage()._process(ctx)
    result_df = processed.current_dataframe

    # Verify exact row count preservation (no k^2 Cartesian explosion)
    assert len(result_df) == n_rows
    # Verify input DataFrame was not mutated in place
    pd.testing.assert_frame_equal(input_df, input_copy)

    # Verify per-sample genotypes and tags remain correctly paired
    for i in range(n_rows):
        assert result_df["annotation_tag"].iloc[i] == f"tag_{i}"
        expected_gt = "0/1" if i % 2 == 0 else "1/1"
        assert result_df["GEN[0].GT"].iloc[i] == expected_gt
        assert result_df["GEN[1].GT"].iloc[i] == "0/0"


def test_genotype_and_allele_count_pairing(tmp_path: Path):
    """Verify distinct genotypes (0/1 vs 1/1) and allele counts stay correctly paired per row."""
    samples = ["CASE1", "CTRL1"]
    # Row 0 has 0/1 (1 allele), Row 1 has 1/1 (2 alleles)
    input_df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1"],
            "POS": ["2000", "2000"],
            "REF": ["C", "C"],
            "ALT": ["G", "G"],
            "GENE": ["TP53", "TP53"],
            "custom_note": ["heterozygous_row", "homozygous_row"],
            "GEN[0].GT": ["0/1", "1/1"],
            "GEN[1].GT": ["0/0", "0/0"],
        }
    )

    ctx = _create_context(tmp_path, input_df, samples)
    processed = VariantAnalysisStage()._process(ctx)
    result = processed.current_dataframe

    assert len(result) == 2
    row0 = result.iloc[0]
    row1 = result.iloc[1]

    # Row 0: 0/1
    assert row0["custom_note"] == "heterozygous_row"
    assert row0["GEN[0].GT"] == "0/1"
    assert row0["GT"] == "CASE1(0/1)"
    if "proband_allele_count" in result.columns:
        assert int(row0["proband_allele_count"]) == 1

    # Row 1: 1/1
    assert row1["custom_note"] == "homozygous_row"
    assert row1["GEN[0].GT"] == "1/1"
    assert row1["GT"] == "CASE1(1/1)"
    if "proband_allele_count" in result.columns:
        assert int(row1["proband_allele_count"]) == 2


@pytest.mark.parametrize("featureid_scenario", ["absent", "repeated", "distinct"])
def test_featureid_handling_with_duplicates(tmp_path: Path, featureid_scenario: str):
    """Cover absent, repeated, and distinct FEATUREID values with duplicate genomic keys."""
    samples = ["SAMPLE1", "SAMPLE2"]
    base_data = {
        "CHROM": ["chr7", "chr7"],
        "POS": ["500", "500"],
        "REF": ["G", "G"],
        "ALT": ["A", "A"],
        "GENE": ["CFTR", "CFTR"],
        "GEN[0].GT": ["0/1", "1/1"],
        "GEN[1].GT": ["0/0", "0/1"],
    }
    if featureid_scenario == "repeated":
        base_data["FEATUREID"] = ["ENST00000001", "ENST00000001"]
    elif featureid_scenario == "distinct":
        base_data["FEATUREID"] = ["ENST00000001", "ENST00000002"]

    input_df = pd.DataFrame(base_data)
    ctx = _create_context(tmp_path, input_df, samples)
    processed = VariantAnalysisStage()._process(ctx)
    result = processed.current_dataframe

    assert len(result) == 2
    assert result["GEN[0].GT"].tolist() == ["0/1", "1/1"]
    assert result["GEN[1].GT"].tolist() == ["0/0", "0/1"]


def test_sanitized_column_names_with_rename_map(tmp_path: Path):
    """Test preservation with sanitized column names (GEN_0__GT) and column_rename_map."""
    samples = ["CASE", "CTRL"]
    input_df = pd.DataFrame(
        {
            "CHROM": ["chr2", "chr2"],
            "POS": ["3000", "3000"],
            "REF": ["T", "T"],
            "ALT": ["C", "C"],
            "GENE": ["MUC1", "MUC1"],
            "ANN_0__FEATUREID": ["TRANS1", "TRANS1"],
            "GEN_0__GT": ["0/1", "1/1"],
            "GEN_1__GT": ["0/0", "0/1"],
        }
    )
    ctx = _create_context(tmp_path, input_df, samples)
    ctx.column_rename_map = {
        "ANN[0].FEATUREID": "ANN_0__FEATUREID",
        "GEN[0].GT": "GEN_0__GT",
        "GEN[1].GT": "GEN_1__GT",
    }

    processed = VariantAnalysisStage()._process(ctx)
    result = processed.current_dataframe

    assert len(result) == 2
    # Verify sanitized names are preserved
    assert "GEN_0__GT" in result.columns
    assert "GEN_1__GT" in result.columns
    assert result["GEN_0__GT"].tolist() == ["0/1", "1/1"]
    assert result["GEN_1__GT"].tolist() == ["0/0", "0/1"]


def test_existing_packed_gt_and_missing_genotypes(tmp_path: Path):
    """Test behavior when packed GT exists and missing genotypes (./.) are present."""
    samples = ["S1", "S2"]
    input_df = pd.DataFrame(
        {
            "CHROM": ["chr3", "chr3"],
            "POS": ["4000", "4000"],
            "REF": ["A", "A"],
            "ALT": ["C", "C"],
            "GENE": ["GENE1", "GENE1"],
            "GT": ["S1(0/1);S2(./.)", "S1(./.);S2(1/1)"],
            "GEN[0].GT": ["0/1", "./."],
            "GEN[1].GT": ["./.", "1/1"],
        }
    )
    ctx = _create_context(tmp_path, input_df, samples)
    processed = VariantAnalysisStage()._process(ctx)
    result = processed.current_dataframe

    assert len(result) == 2
    assert result["GT"].tolist() == ["S1(0/1);S2(./.)", "S1(./.);S2(1/1)"]
    assert result["GEN[0].GT"].tolist() == ["0/1", "./."]
    assert result["GEN[1].GT"].tolist() == ["./.", "1/1"]
