"""Tests for Issue #150: Missing cohort rejection and TSV fallback fail-on-empty handling."""

from __future__ import annotations

import gzip
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.pipeline import _validate_requested_outputs
from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage
from variantcentrifuge.stages.gene_burden_stages import GeneBurdenAnalysisStage
from variantcentrifuge.stages.tsv_output_stages import TSVOutputStage, _tsv_has_data_rows


def test_tsv_has_data_rows_helper(tmp_path: Path):
    """Test _tsv_has_data_rows on plain and gzip files with and without data rows."""
    # Plain empty
    p_empty = tmp_path / "empty.tsv"
    p_empty.write_text("")
    assert not _tsv_has_data_rows(p_empty)

    # Plain header-only
    p_header = tmp_path / "header_only.tsv"
    p_header.write_text("CHROM\tPOS\tREF\tALT\n")
    assert not _tsv_has_data_rows(p_header)

    # Plain with data
    p_data = tmp_path / "with_data.tsv"
    p_data.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\n")
    assert _tsv_has_data_rows(p_data)

    # Gzip header-only
    gz_header = tmp_path / "header_only.tsv.gz"
    with gzip.open(gz_header, "wt") as f:
        f.write("CHROM\tPOS\tREF\tALT\n")
    assert not _tsv_has_data_rows(gz_header)

    # Gzip with data
    gz_data = tmp_path / "with_data.tsv.gz"
    with gzip.open(gz_data, "wt") as f:
        f.write("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\n")
    assert _tsv_has_data_rows(gz_data)


def test_tsv_fallback_fail_on_empty_exits_2(tmp_path: Path):
    """Chunked analysis TSV fallback exits 2 when fail_on_empty is True and output has no rows."""
    workspace = Workspace(tmp_path / "work", "test")
    empty_chunked = tmp_path / "chunked.tsv"
    empty_chunked.write_text("CHROM\tPOS\tREF\tALT\n")

    config = {
        "use_chunked_processing": True,
        "fail_on_empty": True,
        "output_file": str(tmp_path / "output.tsv"),
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.chunked_analysis_tsv = empty_chunked
    context.current_dataframe = None

    stage = TSVOutputStage()
    with pytest.raises(SystemExit) as excinfo:
        stage._process(context)
    assert excinfo.value.code == 2


def test_tsv_fallback_allows_empty_when_fail_on_empty_disabled(tmp_path: Path):
    """Chunked analysis TSV fallback succeeds when fail_on_empty is False."""
    workspace = Workspace(tmp_path / "work", "test")
    empty_chunked = tmp_path / "chunked.tsv"
    empty_chunked.write_text("CHROM\tPOS\tREF\tALT\n")

    out_file = tmp_path / "output.tsv"
    config = {
        "use_chunked_processing": True,
        "fail_on_empty": False,
        "output_file": str(out_file),
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.chunked_analysis_tsv = empty_chunked
    context.current_dataframe = None

    stage = TSVOutputStage()
    stage._process(context)
    assert out_file.exists()
    assert out_file.read_text() == "CHROM\tPOS\tREF\tALT\n"


def test_missing_cohorts_rejected_for_gene_burden(tmp_path: Path):
    """Missing cohorts produce explicit ValueError when gene burden is requested."""
    workspace = Workspace(tmp_path / "work", "test")
    stage = GeneBurdenAnalysisStage()

    # Both missing
    ctx_both = PipelineContext(
        Namespace(),
        {"perform_gene_burden": True, "case_samples": [], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="both case and control cohorts are missing"):
        stage._process(ctx_both)

    # Missing cases
    ctx_no_cases = PipelineContext(
        Namespace(),
        {"perform_gene_burden": True, "case_samples": [], "control_samples": ["CTRL1"]},
        workspace,
    )
    with pytest.raises(ValueError, match="case cohort is missing"):
        stage._process(ctx_no_cases)

    # Missing controls
    ctx_no_controls = PipelineContext(
        Namespace(),
        {"perform_gene_burden": True, "case_samples": ["CASE1"], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="control cohort is missing"):
        stage._process(ctx_no_controls)

    # Disabled burden skips without error
    ctx_disabled = PipelineContext(Namespace(), {"perform_gene_burden": False}, workspace)
    stage._process(ctx_disabled)


def test_missing_cohorts_rejected_for_association(tmp_path: Path):
    """Missing cohorts produce explicit ValueError when association is requested."""
    workspace = Workspace(tmp_path / "work", "test")
    stage = AssociationAnalysisStage()

    # Both missing
    ctx_both = PipelineContext(
        Namespace(),
        {"perform_association": True, "case_samples": [], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="both case and control cohorts are missing"):
        stage._process(ctx_both)

    # Missing cases
    ctx_no_cases = PipelineContext(
        Namespace(),
        {"perform_association": True, "case_samples": [], "control_samples": ["CTRL1"]},
        workspace,
    )
    with pytest.raises(ValueError, match="case cohort is missing"):
        stage._process(ctx_no_cases)

    # Missing controls
    ctx_no_controls = PipelineContext(
        Namespace(),
        {"perform_association": True, "case_samples": ["CASE1"], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="control cohort is missing"):
        stage._process(ctx_no_controls)

    # Disabled association skips without error
    ctx_disabled = PipelineContext(Namespace(), {"perform_association": False}, workspace)
    stage._process(ctx_disabled)


@pytest.mark.parametrize(
    "src_gz,dst_gz,has_data",
    [
        (False, False, True),
        (False, True, True),
        (True, False, True),
        (True, True, True),
        (True, True, False),
    ],
)
def test_chunked_analysis_tsv_roundtrips(
    tmp_path: Path, src_gz: bool, dst_gz: bool, has_data: bool
):
    """Test chunked_analysis_tsv fallback across plain/gz and data/header-only combinations."""
    workspace = Workspace(tmp_path / "work", "test")
    src_file = tmp_path / ("chunked.tsv.gz" if src_gz else "chunked.tsv")
    content = "CHROM\tPOS\tREF\tALT\n" + ("chr1\t100\tA\tT\n" if has_data else "")

    if src_gz:
        with gzip.open(src_file, "wt") as f:
            f.write(content)
    else:
        src_file.write_text(content)

    out_file = tmp_path / ("output.tsv.gz" if dst_gz else "output.tsv")
    config = {
        "use_chunked_processing": True,
        "fail_on_empty": False,
        "output_file": str(out_file),
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.chunked_analysis_tsv = src_file
    context.current_dataframe = None

    stage = TSVOutputStage()
    stage._process(context)
    assert out_file.exists()

    if dst_gz:
        with gzip.open(out_file, "rt") as f:
            read_back = f.read()
    else:
        read_back = out_file.read_text()
    assert read_back == content


def test_chunked_output_file_fallback_branches(tmp_path: Path):
    """Test chunked_output_file fallback branch when chunked_analysis_tsv is absent."""
    workspace = Workspace(tmp_path / "work", "test")
    workspace.output_dir.mkdir(parents=True, exist_ok=True)
    stage = TSVOutputStage()

    # 1. Non-empty file
    chunked_file = workspace.output_dir / "results.tsv"
    chunked_file.write_text("CHROM\tPOS\n1\t100\n")
    ctx_ok = PipelineContext(
        Namespace(),
        {
            "use_chunked_processing": True,
            "chunked_output_file": "results.tsv",
            "fail_on_empty": True,
        },
        workspace,
    )
    ctx_ok.current_dataframe = None
    stage._process(ctx_ok)
    assert ctx_ok.final_output_path == chunked_file

    # 2. Header-only with fail_on_empty=True exits 2
    chunked_empty = workspace.output_dir / "empty_results.tsv"
    chunked_empty.write_text("CHROM\tPOS\n")
    ctx_fail = PipelineContext(
        Namespace(),
        {
            "use_chunked_processing": True,
            "chunked_output_file": "empty_results.tsv",
            "fail_on_empty": True,
        },
        workspace,
    )
    ctx_fail.current_dataframe = None
    with pytest.raises(SystemExit) as excinfo:
        stage._process(ctx_fail)
    assert excinfo.value.code == 2

    # 3. Header-only with fail_on_empty=False succeeds
    ctx_pass_empty = PipelineContext(
        Namespace(),
        {
            "use_chunked_processing": True,
            "chunked_output_file": "empty_results.tsv",
            "fail_on_empty": False,
        },
        workspace,
    )
    ctx_pass_empty.current_dataframe = None
    stage._process(ctx_pass_empty)
    assert ctx_pass_empty.final_output_path == chunked_empty


def test_zero_variant_cohort_produces_sidecars_validating_outputs(tmp_path: Path):
    """Valid zero-variant cohorts produce header-only sidecars that pass output validation."""
    workspace = Workspace(tmp_path / "work", "test")
    workspace.output_dir.mkdir(parents=True, exist_ok=True)

    burden_out = str(workspace.output_dir / "burden.tsv")
    assoc_out = str(workspace.output_dir / "assoc.tsv")
    final_out = workspace.output_dir / "final.tsv"
    final_out.write_text("CHROM\tPOS\tREF\tALT\n")

    config = {
        "output_file": str(final_out),
        "case_samples": ["CASE1"],
        "control_samples": ["CTRL1"],
        "perform_gene_burden": True,
        "gene_burden_output": burden_out,
        "perform_association": True,
        "association_output": assoc_out,
        "gzip_intermediates": False,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.current_dataframe = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "GENE", "GT"])
    context.final_output_path = final_out

    # Process gene burden on empty df -> writes header-only sidecar
    GeneBurdenAnalysisStage()._process(context)
    assert Path(burden_out).exists()
    assert _tsv_has_data_rows(Path(burden_out)) is False

    # Process association on empty df -> writes header-only sidecar
    AssociationAnalysisStage()._process(context)
    assert Path(assoc_out).exists()
    assert _tsv_has_data_rows(Path(assoc_out)) is False

    # Output validation must pass without raising RuntimeError
    _validate_requested_outputs(context)


def test_restored_header_only_checkpoint_fail_on_empty_exits_2(tmp_path: Path):
    """Restored header-only checkpoint DataFrame exits 2 under fail_on_empty=True."""
    workspace = Workspace(tmp_path / "work", "test")
    config = {
        "output_file": str(tmp_path / "final.tsv"),
        "fail_on_empty": True,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.current_dataframe = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT"])

    stage = TSVOutputStage()
    with pytest.raises(SystemExit) as excinfo:
        stage._process(context)
    assert excinfo.value.code == 2
