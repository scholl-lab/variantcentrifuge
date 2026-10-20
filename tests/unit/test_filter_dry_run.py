"""SnpSift filter dry run on header + synthetic record (#211)."""

import shutil
import subprocess
from pathlib import Path

import pytest

from variantcentrifuge import filter_dry_run
from variantcentrifuge.filter_dry_run import _synthetic_vcf, validate_snpsift_filter

HDR = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr7,length=10>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)
FIXTURE = Path(__file__).parents[1] / "fixtures/1000g_cohort/1000g_nephro.GRCh38.annotated.vcf.gz"


def test_synthetic_record_with_samples():
    last = _synthetic_vcf(HDR).splitlines()[-1]
    assert last == "chr7\t1\t.\tA\tC\t.\t.\t.\tGT\t./.\t./."


def test_synthetic_record_no_contig_no_samples():
    h = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    assert _synthetic_vcf(h).splitlines()[-1] == "1\t1\t.\tA\tC\t.\t.\t."


def test_extra_info_declared_in_header():
    assert "##INFO=<ID=RAW_FILTER," in _synthetic_vcf(HDR, ["RAW_FILTER"])
    assert "RAW_FILTER" not in _synthetic_vcf(HDR)


def _mock(monkeypatch, rc, stderr=""):
    calls = []

    def fake(cmd, **kw):
        calls.append(cmd)
        if cmd[0] == "bcftools":
            return subprocess.CompletedProcess(cmd, 0, HDR, "")
        return subprocess.CompletedProcess(cmd, rc, "", stderr)

    monkeypatch.setattr(filter_dry_run.subprocess, "run", fake)
    return calls


def test_accepts(monkeypatch):
    calls = _mock(monkeypatch, 0)
    validate_snpsift_filter("x.vcf.gz", "(AF > 1)", ["-Xmx1g"])
    assert calls[1][:4] == ["SnpSift", "-Xmx1g", "filter", "(AF > 1)"]


def test_rejects_with_snpsift_message_and_expression(monkeypatch):
    _mock(monkeypatch, 1, "RuntimeException: INFO field 'HIGH' not found in VCF header")
    with pytest.raises(ValueError) as e:
        validate_snpsift_filter("x.vcf.gz", "(ANN[0].IMPACT = HIGH)", [])
    assert "INFO field 'HIGH' not found" in str(e.value)
    assert "(ANN[0].IMPACT = HIGH)" in str(e.value)


@pytest.mark.skipif(
    not (shutil.which("SnpSift") and shutil.which("bcftools")), reason="SnpSift/bcftools missing"
)
class TestRealSnpSift:
    @pytest.mark.parametrize(
        "expr",
        ["(ANN[0].IMPACT = HIGH)", "(FOO > 1)", "(AF >"],
        ids=["unquoted", "unknown", "syntax"],
    )
    def test_bad(self, expr):
        with pytest.raises(ValueError, match="SnpSift rejected"):
            validate_snpsift_filter(str(FIXTURE), expr, ["-Xmx1g"])

    @pytest.mark.parametrize(
        "expr",
        [
            "((ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE'))",
            "(ANN[ANY].IMPACT has 'HIGH') & (AF < 0.5)",
            "(GEN[*].GT = '0/1') | (GEN[0].GT = '0/1')",
            "(countHom() >= 1) | (countHet() > 0)",
            "(QUAL >= 30) | (na AC)",
        ],
    )
    def test_good(self, expr):
        validate_snpsift_filter(str(FIXTURE), expr, ["-Xmx1g"])

    def test_raw_filter_declared_when_requested(self):
        validate_snpsift_filter(str(FIXTURE), "(exists RAW_FILTER)", ["-Xmx1g"], ["RAW_FILTER"])
        with pytest.raises(ValueError):
            validate_snpsift_filter(str(FIXTURE), "(exists RAW_FILTER)", ["-Xmx1g"])
