"""Mask label + annotation provenance recorded in artifacts and compared by meta (#197)."""

from __future__ import annotations

import gzip
import hashlib
import json
import logging

import pytest

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.meta.artifact import load_artifact, save_artifact
from variantcentrifuge.association.meta.manifest import (
    MetaCompatibilityError,
    validate_comparability,
)
from variantcentrifuge.association.meta.provenance import (
    PROVENANCE_FIELDS,
    collect_provenance,
)

KEYS = {"G": [f"1:{100 + j}:A:G" for j in range(5)]}

HEADER = (
    "##fileformat=VCFv4.2\n"
    '##SnpEffVersion="5.2 (build 2023-09-29 06:17), by Pablo Cingolani"\n'
    '##SnpEffCmd="SnpEff  -noStats -lof GRCh38.p14 cohort.vcf "\n'
    '##SnpSiftCmd="SnpSift dbnsfp -f CADD_phred -db /db/dbNSFP4.4a.txt.gz -"\n'
    '##SnpSiftCmd="SnpSift annotate -info AF /db/gnomad.exomes.v4.1.sites.vcf.gz -"\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _vcf(tmp_path, text=HEADER, gz=True):
    path = tmp_path / ("in.vcf.gz" if gz else "in.vcf")
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return str(path)


def test_collect_provenance_from_vcf_header_and_files(tmp_path):
    bed = tmp_path / "r.bed"
    bed.write_text("1\t1\t100\n")
    prov = collect_provenance({"vcf_file": _vcf(tmp_path), "regions_bed": str(bed)})
    assert prov == {
        "snpeff_version": "5.2 (build 2023-09-29 06:17), by Pablo Cingolani",
        "snpeff_genome_db": "GRCh38.p14",
        "dbnsfp_version": "dbNSFP4.4a.txt.gz",
        "gnomad_source": "gnomad.exomes.v4.1.sites.vcf.gz",
        "regions_bed_sha256": hashlib.sha256(b"1\t1\t100\n").hexdigest(),
        "gene_file_sha256": "none",
    }
    assert set(prov) == set(PROVENANCE_FIELDS)


def test_collect_provenance_absent_header_lines_are_unknown(tmp_path):
    vcf = _vcf(tmp_path, "##fileformat=VCFv4.2\n#CHROM\tPOS\n", gz=False)
    prov = collect_provenance({"vcf_file": vcf})
    assert prov["snpeff_version"] is None and prov["dbnsfp_version"] is None
    assert prov["gnomad_source"] is None and prov["snpeff_genome_db"] is None
    assert collect_provenance({})["snpeff_version"] is None  # no VCF configured


def test_snpeff_genome_from_stdin_command(tmp_path):
    vcf = _vcf(tmp_path, '##SnpEffCmd="SnpEff ann -noStats GRCh37.75 "\n#CHROM\n', gz=False)
    assert collect_provenance({"vcf_file": vcf})["snpeff_genome_db"] == "GRCh37.75"


def test_label_and_provenance_round_trip(tmp_path):
    art = make_stratum("A", seed=1, gene_keys=KEYS)
    art.mask_label = "pLoF"
    art.provenance = {"dbnsfp_version": "dbNSFP4.4a.txt.gz"}
    save_artifact(art, str(tmp_path / "a"))
    back = load_artifact(str(tmp_path / "a"))
    assert back.mask_label == "pLoF"
    assert back.provenance == {"dbnsfp_version": "dbNSFP4.4a.txt.gz"}


def test_artifact_without_provenance_loads_as_unknown(tmp_path):
    save_artifact(make_stratum("A", seed=1, gene_keys=KEYS), str(tmp_path / "a"))
    meta_path = tmp_path / "a.meta.json"
    meta = json.loads(meta_path.read_text())
    meta.pop("provenance"), meta.pop("mask_label")
    meta_path.write_text(json.dumps(meta))
    back = load_artifact(str(tmp_path / "a"))
    assert back.provenance == {} and back.mask_label is None


def _pair(prov_a, prov_b):
    a = make_stratum("A", seed=1, gene_keys=KEYS)
    b = make_stratum("B", seed=2, gene_keys=KEYS)
    a.provenance, b.provenance = prov_a, prov_b
    return [a, b]


def test_provenance_mismatch_refused_naming_field_and_strata():
    arts = _pair({"dbnsfp_version": "dbNSFP4.4a"}, {"dbnsfp_version": "dbNSFP4.3a"})
    with pytest.raises(MetaCompatibilityError, match=r"dbnsfp_version.*A=dbNSFP4.4a.*B=dbNSFP4.3a"):
        validate_comparability(arts)


def test_allow_annotation_mismatch_downgrades_to_warning(caplog):
    arts = _pair({"regions_bed_sha256": "aa"}, {"regions_bed_sha256": "bb"})
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        validate_comparability(arts, allow_annotation_mismatch=True)
    assert "regions_bed_sha256" in caplog.text


def test_unknown_provenance_warns_but_pools(caplog):
    arts = _pair({"snpeff_version": "5.2"}, {})
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        validate_comparability(arts)
    assert "snpeff_version" in caplog.text and "unknown" in caplog.text


def test_cli_allow_annotation_mismatch(tmp_path):
    from variantcentrifuge.association.meta.cli import main

    paths = []
    for art in _pair({"gnomad_source": "v4.1"}, {"gnomad_source": "v2.1"}):
        save_artifact(art, str(tmp_path / art.stratum_id))
        paths.append(art.stratum_id)
    manifest = tmp_path / "m.tsv"
    manifest.write_text("stratum_id\tartifact_path\n" + "".join(f"{p}\t{p}\n" for p in paths))
    out = str(tmp_path / "o.tsv")
    assert main(["--manifest", str(manifest), "--output", out]) == 2
    assert main(["--manifest", str(manifest), "--output", out, "--allow-annotation-mismatch"]) == 0
