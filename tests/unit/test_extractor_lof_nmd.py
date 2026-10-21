"""LOF[0].* and NMD[0].* subfield extraction (#234)."""

import shutil
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.extractor import (
    build_bcftools_format_string,
    extract_fields_bcftools,
    extract_fields_snpsift,
    get_expected_extraction_columns,
)

SUBFIELDS = ["GENE", "GENEID", "NUMTR", "PERC"]
FIELDS = ["CHROM", "POS"] + [f"{tag}[0].{s}" for tag in ("LOF", "NMD") for s in SUBFIELDS]
EXPECTED = {
    "LOF": ["GENE_A", "ENSGA", "1", "0.75"],
    "NMD": ["GENE_N", "ENSGN", "4", "0.25"],
}


@pytest.mark.unit
def test_lof_subfields_use_raw_lof_column_once():
    format_str, columns = build_bcftools_format_string(FIELDS)
    assert format_str.count("%INFO/LOF") == 1 and "LOF[0]" not in format_str
    assert columns == ["CHROM", "POS", "LOF", "NMD"]
    assert get_expected_extraction_columns(FIELDS) == FIELDS


@pytest.mark.unit
def test_bare_lof_and_nmd_fields_keep_raw_column():
    assert get_expected_extraction_columns(["CHROM", "LOF", "NMD"]) == ["CHROM", "LOF", "NMD"]


def _write_vcf(path: Path) -> None:
    info = "LOF=(GENE_A|ENSGA|1|0.75),(GENE_B|ENSGB|2|0.50);NMD=(GENE_N|ENSGN|4|0.25)"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects">\n'
        '##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t100\t.\tA\tT\t50\tPASS\t{info}\n"
        "chr1\t200\t.\tG\tC\t50\tPASS\t.\n",
        encoding="utf-8",
    )


@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools required")
def test_bcftools_extracts_every_lof_and_nmd_subfield(tmp_path):
    vcf = tmp_path / "in.vcf"
    _write_vcf(vcf)
    out = tmp_path / "out.tsv"
    extract_fields_bcftools(str(vcf), " ".join(FIELDS), {"temp_dir": str(tmp_path)}, str(out))
    df = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
    assert list(df.columns) == FIELDS
    for tag, values in EXPECTED.items():
        assert df.loc[0, [f"{tag}[0].{s}" for s in SUBFIELDS]].tolist() == values
    assert df.loc[1, FIELDS[2:]].tolist() == ["NA"] * 8


@pytest.mark.skipif(not shutil.which("SnpSift"), reason="SnpSift required")
def test_snpsift_extracts_every_lof_and_nmd_subfield(tmp_path):
    # SnpSift parses LOF/NMD natively; values chosen so Java float formatting is verbatim.
    vcf = tmp_path / "in.vcf"
    _write_vcf(vcf)
    out = tmp_path / "out.tsv"
    extract_fields_snpsift(str(vcf), " ".join(FIELDS), {}, str(out))
    df = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
    for tag, values in EXPECTED.items():
        assert df.loc[0, [f"{tag}_{s}" for s in SUBFIELDS]].tolist() == values


@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools required")
@pytest.mark.parametrize("tag", ["LOF", "NMD"])
def test_bare_tag_plus_subfield_adds_raw_column_once(tmp_path, tag):
    fields = ["CHROM", tag, f"{tag}[0].PERC"]
    assert build_bcftools_format_string(fields)[1] == ["CHROM", tag]
    assert get_expected_extraction_columns(fields) == fields
    vcf = tmp_path / "in.vcf"
    _write_vcf(vcf)
    out = tmp_path / "out.tsv"
    extract_fields_bcftools(str(vcf), " ".join(fields), {"temp_dir": str(tmp_path)}, str(out))
    df = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
    assert list(df.columns) == fields
    raw = {"LOF": "(GENE_A|ENSGA|1|0.75),(GENE_B|ENSGB|2|0.50)", "NMD": "(GENE_N|ENSGN|4|0.25)"}
    assert df.loc[0].tolist() == ["chr1", raw[tag], EXPECTED[tag][3]]
