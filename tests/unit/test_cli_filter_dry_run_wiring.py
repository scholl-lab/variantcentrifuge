"""cli.main must call the SnpSift filter dry run before the pipeline starts (#211)."""

import shutil
from pathlib import Path

from variantcentrifuge import cli as cli_module

SRC = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")


def _rec(calls):
    return lambda *a: calls.append(a)


def _run(monkeypatch, tmp_path, extra, dry_run):
    vcf = tmp_path / "input.vcf.gz"
    shutil.copy(SRC, vcf)
    shutil.copy(str(SRC) + ".tbi", str(vcf) + ".tbi")
    order = []
    # overrides the autouse conftest stub
    monkeypatch.setattr("variantcentrifuge.filter_dry_run.validate_snpsift_filter", dry_run)
    monkeypatch.setattr(cli_module, "run_pipeline", lambda args: order.append("pipeline"))
    monkeypatch.setattr(cli_module, "validate_required_tools", lambda *a, **k: None)
    argv = ["vc", "--vcf-file", str(vcf), "--gene-name", "all", "--output-dir", str(tmp_path)]
    monkeypatch.setattr("sys.argv", argv + extra)
    return cli_module.main(), order, vcf


def test_dry_run_called_before_pipeline(monkeypatch, tmp_path):
    calls = []
    rc, order, _ = _run(monkeypatch, tmp_path, ["--filters", "(AF > 0.1)"], _rec(calls))
    assert rc == 0 and order == ["pipeline"]
    assert len(calls) == 1 and calls[0][0].endswith("input.vcf.gz") and calls[0][1] == "(AF > 0.1)"


def test_bad_filter_exits_before_pipeline(monkeypatch, tmp_path):
    def dry(*a):
        raise ValueError("SnpSift rejected the filter expression")

    rc, order, _ = _run(monkeypatch, tmp_path, ["--filters", "(AF > 0.1)"], dry)
    assert rc != 0 and order == []


def test_late_filtering_skips_dry_run(monkeypatch, tmp_path):
    calls = []
    _run(
        monkeypatch,
        tmp_path,
        ["--filters", "(AF > 0.1)", "--late-filtering"],
        _rec(calls),
    )
    assert not calls
