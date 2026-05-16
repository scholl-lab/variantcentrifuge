import os
import subprocess
import sys
from unittest.mock import patch

import pytest

from variantcentrifuge.filters import (
    _raise_for_snpsift_filter_stderr,
    apply_snpsift_filter,
)
from variantcentrifuge.utils import run_command


def test_run_command_can_return_completed_process():
    result = run_command([sys.executable, "-c", "print('ok')"], return_result=True)

    assert isinstance(result, subprocess.CompletedProcess)
    assert result.returncode == 0
    assert result.stdout.strip() == "ok"


@pytest.mark.parametrize(
    "stderr",
    [
        "line 1:126 mismatched input ')' expecting {'*'}",
        "line 1:293 missing ')' at '<EOF>'",
        "line 1:10 token recognition error at: '@'",
        "Error parsing expression",
        "LexerNoViableAltException('@')",
        "Unknown parameter 'BAD'",
        "INFO field 'ClinVar_CLNSIG' not found",
        "Exception in thread \"main\" java.lang.RuntimeException: Cannot parse EffectType 'HIGH'",
    ],
)
def test_snpsift_filter_stderr_parser_diagnostics_raise(stderr):
    with pytest.raises(RuntimeError, match="SnpSift filter reported parser diagnostics"):
        _raise_for_snpsift_filter_stderr(stderr, "BAD_EXPR")


def test_snpsift_filter_stderr_allows_empty_or_nonfatal_warnings():
    _raise_for_snpsift_filter_stderr("", "EXPR")
    _raise_for_snpsift_filter_stderr("Some non-fatal note", "EXPR")


@patch("variantcentrifuge.filters.run_command")
@patch("variantcentrifuge.filters.tempfile.mkstemp")
@patch("variantcentrifuge.filters._snpsift_memory_flag", return_value="-Xmx1g")
def test_apply_snpsift_filter_raises_before_compressing_when_snpsift_reports_parser_error(
    mock_memory,
    mock_mkstemp,
    mock_run_command,
    tmp_path,
):
    input_vcf = tmp_path / "input.vcf.gz"
    output_vcf = tmp_path / "output.vcf.gz"
    tmp_vcf = tmp_path / "tmp-filter.vcf"
    input_vcf.write_text("placeholder")
    tmp_vcf.write_text("")

    tmp_fd = os.open(tmp_vcf, os.O_RDWR)
    mock_mkstemp.return_value = (tmp_fd, str(tmp_vcf))

    mock_run_command.return_value = subprocess.CompletedProcess(
        args=["SnpSift"],
        returncode=0,
        stdout="",
        stderr="line 1:126 mismatched input ')' expecting {'*'}",
    )

    with pytest.raises(RuntimeError, match="SnpSift filter reported parser diagnostics"):
        apply_snpsift_filter(str(input_vcf), "BAD_EXPR", {"threads": 1}, str(output_vcf))

    assert mock_run_command.call_count == 1
    assert not tmp_vcf.exists()
