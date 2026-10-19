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
