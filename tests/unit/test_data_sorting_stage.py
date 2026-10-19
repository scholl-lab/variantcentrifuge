"""Unit tests for DataSortingStage and _sort_tsv_by_gene."""

import contextlib
import gzip
from pathlib import Path

import pytest

from variantcentrifuge.stages.post_processing_stages import DataSortingStage


def test_sort_tsv_by_gene_prevents_command_injection(tmp_path: Path) -> None:
    """Verify command substitutions in memory_limit and temp_dir are not executed."""
    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "output.tsv"
    marker_file = tmp_path / "injected.marker"

    # Minimal TSV with header and two rows
    input_file.write_text("GENE\tVARIANT\nBRCA2\tvar1\nBRCA1\tvar2\n")

    injected_cmd = f"$(touch {marker_file})"
    malicious_memory = f"1G{injected_cmd}"
    malicious_temp = str(tmp_path / f"tmp{injected_cmd}")

    stage = DataSortingStage()

    # The function should run sort safely or fail in sort, but NOT execute the subshell
    with contextlib.suppress(Exception):
        stage._sort_tsv_by_gene(
            input_file=str(input_file),
            output_file=str(output_file),
            gene_column="GENE",
            memory_limit=malicious_memory,
            parallel=1,
            temp_dir=malicious_temp,
        )

    assert not marker_file.exists(), (
        "Command substitution was executed via unquoted shell argument!"
    )


@pytest.mark.parametrize(
    "input_gz,output_gz",
    [
        (False, False),
        (False, True),
        (True, False),
        (True, True),
    ],
)
def test_sort_tsv_by_gene_all_compression_combinations(
    tmp_path: Path, input_gz: bool, output_gz: bool
) -> None:
    """Verify that all four plain/gzip combinations preserve header, rows, and sort order."""
    header = "GENE\tCHROM\tPOS\n"
    rows = ["TP53\tchr17\t7673800\n", "BRCA1\tchr17\t43044295\n", "ATM\tchr11\t108223000\n"]
    raw_content = header + "".join(rows)

    in_ext = ".tsv.gz" if input_gz else ".tsv"
    out_ext = ".tsv.gz" if output_gz else ".tsv"
    input_file = tmp_path / f"input{in_ext}"
    output_file = tmp_path / f"output{out_ext}"

    if input_gz:
        with gzip.open(input_file, "wt") as f:
            f.write(raw_content)
    else:
        input_file.write_text(raw_content)

    stage = DataSortingStage()
    stage._sort_tsv_by_gene(
        input_file=str(input_file),
        output_file=str(output_file),
        gene_column="GENE",
        memory_limit="1G",
        parallel=1,
    )

    if output_gz:
        with gzip.open(output_file, "rt") as f:
            out_content = f.read()
    else:
        out_content = output_file.read_text()

    out_lines = out_content.splitlines()
    assert len(out_lines) == 4
    assert out_lines[0] == "GENE\tCHROM\tPOS"
    # Alphabetical order by GENE: ATM, BRCA1, TP53
    assert out_lines[1].startswith("ATM\t")
    assert out_lines[2].startswith("BRCA1\t")
    assert out_lines[3].startswith("TP53\t")


def test_data_sorting_stage_process_resolves_ann_zero_gene(tmp_path: Path) -> None:
    """Verify DataSortingStage._process recognizes ANN[0].GENE and sorts correctly."""
    from argparse import Namespace
    from unittest.mock import Mock

    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace

    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir()
    ws = Mock(spec=Workspace)
    ws.base_name = "test_sample"
    ws.intermediate_dir = inter_dir
    ws.get_intermediate_path = lambda name: inter_dir / name

    input_file = inter_dir / "test_sample.extracted.tsv"
    input_file.write_text("CHROM\tPOS\tANN[0].GENE\nchr17\t100\tTP53\nchr17\t200\tBRCA1\n")

    context = PipelineContext(
        args=Namespace(),
        config={"sort_parallel": 1, "sort_memory_limit": "1G"},
        workspace=ws,
    )
    context.extracted_tsv = input_file

    stage = DataSortingStage()
    updated_ctx = stage._process(context)

    assert updated_ctx.extracted_tsv != input_file
    assert updated_ctx.extracted_tsv.exists()

    content = updated_ctx.extracted_tsv.read_text().splitlines()
    assert len(content) == 3
    assert content[0] == "CHROM\tPOS\tANN[0].GENE"
    assert content[1] == "chr17\t200\tBRCA1"
    assert content[2] == "chr17\t100\tTP53"


def test_data_sorting_stage_process_custom_gene_column(tmp_path: Path) -> None:
    """Verify DataSortingStage._process supports custom gene_column_name."""
    from argparse import Namespace
    from unittest.mock import Mock

    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace

    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir()
    ws = Mock(spec=Workspace)
    ws.base_name = "test_sample"
    ws.intermediate_dir = inter_dir
    ws.get_intermediate_path = lambda name: inter_dir / name

    input_file = inter_dir / "test_sample.extracted.tsv"
    input_file.write_text("CHROM\tPOS\tCUSTOM_GENE\nchr17\t100\tTP53\nchr17\t200\tBRCA1\n")

    context = PipelineContext(
        args=Namespace(),
        config={"gene_column_name": "CUSTOM_GENE", "sort_parallel": 1, "sort_memory_limit": "1G"},
        workspace=ws,
    )
    context.extracted_tsv = input_file

    stage = DataSortingStage()
    updated_ctx = stage._process(context)

    assert updated_ctx.extracted_tsv.exists()
    content = updated_ctx.extracted_tsv.read_text().splitlines()
    assert content[1].endswith("BRCA1")
    assert content[2].endswith("TP53")


def test_sort_tsv_by_gene_stable_order_with_ties(tmp_path: Path) -> None:
    """Verify that GNU sort --stable preserves relative order for tied gene keys."""
    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "output.tsv"

    input_file.write_text(
        "GENE\tPOS\tORDER\nTP53\t100\t1\nBRCA1\t200\ttie_A\nBRCA1\t300\ttie_B\nATM\t400\t4\n"
    )

    stage = DataSortingStage()
    stage._sort_tsv_by_gene(
        input_file=str(input_file),
        output_file=str(output_file),
        gene_column="GENE",
        memory_limit="1G",
        parallel=1,
    )

    lines = output_file.read_text().splitlines()
    assert len(lines) == 5
    assert lines[1].startswith("ATM\t")
    # Tied BRCA1 rows must preserve input order: tie_A before tie_B
    assert lines[2] == "BRCA1\t200\ttie_A"
    assert lines[3] == "BRCA1\t300\ttie_B"
    assert lines[4].startswith("TP53\t")


def test_sort_tsv_by_gene_paths_with_spaces_quotes_semicolons(tmp_path: Path) -> None:
    """Verify sorting safely handles paths with spaces, quotes, and semicolons."""
    special_dir = tmp_path / "dir with spaces 'quoted' and;semicolon"
    special_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = special_dir / "tmp dir;nested"
    temp_dir.mkdir(parents=True, exist_ok=True)

    input_file = special_dir / "input with spaces.tsv"
    output_file = special_dir / "output with spaces.tsv"

    input_file.write_text("GENE\tPOS\nTP53\t100\nBRCA1\t200\n")

    stage = DataSortingStage()
    stage._sort_tsv_by_gene(
        input_file=str(input_file),
        output_file=str(output_file),
        gene_column="GENE",
        memory_limit="1G",
        parallel=1,
        temp_dir=str(temp_dir),
    )

    assert output_file.exists()
    lines = output_file.read_text().splitlines()
    assert lines[0] == "GENE\tPOS"
    assert lines[1] == "BRCA1\t200"
    assert lines[2] == "TP53\t100"
