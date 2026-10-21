"""#233: not_artefact_af5 preset and the missing-INFO/AF warning."""

import argparse
import logging
import shutil

import pytest

from variantcentrifuge.cli_config import resolve_presets_and_tumor_normal
from variantcentrifuge.config import load_config

AF_HEADER = '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'


def _vcf(tmp_path, with_af: bool):
    path = tmp_path / ("af.vcf" if with_af else "no_af.vcf")
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + (AF_HEADER if with_af else "")
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )
    return str(path)


def _resolve(vcf: str, preset: str):
    args = argparse.Namespace(preset=[preset], final_filter=None, vcf_file=vcf)
    return resolve_presets_and_tumor_normal(args, load_config(), None)


@pytest.mark.unit
def test_presets_keep_not_artefact_and_add_af5():
    presets = load_config()["presets"]
    assert presets["not_artefact"] == "((AC[ANY] <= 513) & (QUAL[ANY] >= 100))"
    assert presets["not_artefact_af5"] == "((AF[ANY] <= 0.05) & (QUAL[ANY] >= 100))"


@pytest.mark.unit
@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools required")
@pytest.mark.parametrize("with_af", [False, True])
def test_af_preset_warns_only_without_info_af(tmp_path, caplog, with_af):
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        filters, _ = _resolve(_vcf(tmp_path, with_af), "not_artefact_af5")
    assert "AF[ANY] <= 0.05" in filters
    warned = any("INFO/AF" in r.getMessage() for r in caplog.records)
    assert warned is not with_af


@pytest.mark.unit
@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools required")
def test_non_af_preset_does_not_warn(tmp_path, caplog):
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        _resolve(_vcf(tmp_path, False), "not_artefact")
    assert not any("INFO/AF" in r.getMessage() for r in caplog.records)


@pytest.mark.unit
def test_af_preset_with_missing_vcf_defers_to_not_found_error(tmp_path):
    from variantcentrifuge.validators import validate_vcf_file

    missing = str(tmp_path / "missing.vcf.gz")
    _resolve(missing, "not_artefact_af5")  # must not crash reading the header
    with pytest.raises(FileNotFoundError, match="VCF file not found"):
        validate_vcf_file(missing)


@pytest.mark.unit
@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools required")
@pytest.mark.parametrize("kind", ["filter", "final_filter"])
def test_af_in_filter_or_final_filter_warns_without_info_af(tmp_path, caplog, kind):
    expr = "(AF[ANY] <= 0.01)"
    args = argparse.Namespace(
        preset=None,
        final_filter=expr if kind == "final_filter" else None,
        vcf_file=_vcf(tmp_path, False),
    )
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        resolve_presets_and_tumor_normal(args, load_config(), expr if kind == "filter" else None)
    assert any("INFO/AF" in r.getMessage() for r in caplog.records)
