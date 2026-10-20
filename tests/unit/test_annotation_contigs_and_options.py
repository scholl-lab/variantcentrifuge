"""Unit tests for contig prefix normalization and workflow option forwarding (#151)."""

import gzip
import shlex
from pathlib import Path
from unittest.mock import Mock

import pytest

from variantcentrifuge.annotator import (
    _find_region_overlaps,
    _parse_bed_line,
)
from variantcentrifuge.cli import create_parser
from variantcentrifuge.filters import validate_bed_vcf_contigs
from workflow.rules.helpers import build_vc_flags


# ---------------------------------------------------------------------------
# 1. Contig prefix normalization tests
# ---------------------------------------------------------------------------
def test_parse_bed_line_removeprefix() -> None:
    """Verify that _parse_bed_line only removes leading 'chr' prefix."""
    # Standard prefixes
    assert _parse_bed_line("chr1\t100\t200\treg1") == ("1", 100, 200, "reg1")
    assert _parse_bed_line("1\t100\t200\treg1") == ("1", 100, 200, "reg1")
    assert _parse_bed_line("chrX\t100\t200\tregX") == ("X", 100, 200, "regX")
    assert _parse_bed_line("chr18\t100\t200\treg18") == ("18", 100, 200, "reg18")

    # Double prefix or non-prefix 'chr' inside contig name
    assert _parse_bed_line("chrchr1\t100\t200\tdouble") == ("chr1", 100, 200, "double")


def test_find_region_overlaps_removeprefix() -> None:
    """Verify _find_region_overlaps handles contig prefixes symmetrically."""
    tree_mock = Mock()
    overlap_mock = Mock()
    overlap_mock.data = "custom_reg"
    tree_mock.overlap = Mock(return_value=[overlap_mock])
    regions_by_chrom = {"1": tree_mock}

    import pandas as pd

    row_chr1 = pd.Series({"CHROM": "chr1", "POS": 150})
    row_1 = pd.Series({"CHROM": "1", "POS": 150})

    assert _find_region_overlaps(row_chr1, regions_by_chrom) == ["Region=custom_reg"]
    assert _find_region_overlaps(row_1, regions_by_chrom) == ["Region=custom_reg"]


# ---------------------------------------------------------------------------
# 2. BED vs VCF contig prefix compatibility tests
# ---------------------------------------------------------------------------
def _write_vcf(path: Path, contigs: list[str], data_lines: list[str] | None = None) -> None:
    lines = ["##fileformat=VCFv4.2"]
    for c in contigs:
        lines.append(f"##contig=<ID={c},length=10000>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    if data_lines:
        lines.extend(data_lines)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
    else:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_bed(path: Path, lines: list[str]) -> None:
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
    else:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_contig_prefix_mismatch_rejected(tmp_path: Path) -> None:
    """Clear prefix mismatch (e.g. BED has '1' and VCF has 'chr1') must be rejected."""
    vcf = tmp_path / "test.vcf.gz"
    bed = tmp_path / "test.bed"

    _write_vcf(vcf, ["chr1", "chr2"])
    _write_bed(bed, ["1\t100\t200\tregion1"])

    with pytest.raises(ValueError, match="Contig prefix mismatch"):
        validate_bed_vcf_contigs(str(vcf), str(bed))


def test_contig_prefix_mismatch_bed_chr_vcf_nochr_rejected(tmp_path: Path) -> None:
    """Clear prefix mismatch (BED has 'chr1' and VCF has '1') must be rejected."""
    vcf = tmp_path / "test.vcf.gz"
    bed = tmp_path / "test.bed"

    _write_vcf(vcf, ["1", "2"])
    _write_bed(bed, ["chr1\t100\t200\tregion1"])

    with pytest.raises(ValueError, match="Contig prefix mismatch"):
        validate_bed_vcf_contigs(str(vcf), str(bed))


def test_contig_matching_legitimate_empty_allowed(tmp_path: Path) -> None:
    """Matching contigs with no variants in region must not raise validation error."""
    vcf = tmp_path / "test.vcf.gz"
    bed = tmp_path / "test.bed"

    _write_vcf(vcf, ["chr1", "chr2"])
    _write_bed(bed, ["chr1\t500\t600\tempty_region"])

    # Must return without error
    validate_bed_vcf_contigs(str(vcf), str(bed))


def test_contig_unmatched_warning_logged(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Unmatched contigs (some match, some do not) must log a warning."""
    vcf = tmp_path / "test.vcf.gz"
    bed = tmp_path / "test.bed"

    _write_vcf(vcf, ["chr1", "chr2"])
    _write_bed(bed, ["chr1\t100\t200\tregion1", "chr99\t100\t200\tunmatched"])

    with caplog.at_level("WARNING"):
        validate_bed_vcf_contigs(str(vcf), str(bed))

    assert any("not found in VCF" in record.message for record in caplog.records)


# ---------------------------------------------------------------------------
# 3. Workflow helpers option forwarding tests
# ---------------------------------------------------------------------------
def test_workflow_interval_expand_rejected() -> None:
    """Setting interval_expand at workflow level must raise actionable ValueError."""
    vc_config = {"interval_expand": 1000}
    with pytest.raises(ValueError, match="interval_expand"):
        build_vc_flags(vc_config, {})


def test_workflow_sample_indices_index_zero_preserved() -> None:
    """Sample indices 0 must be forwarded (--tumor-sample-index 0)."""
    vc_config = {
        "tumor_sample_index": 0,
        "normal_sample_index": 1,
    }
    flags = build_vc_flags(vc_config, {})
    assert "--tumor-sample-index 0" in flags
    assert "--normal-sample-index 1" in flags


def test_workflow_normal_sample_index_zero_preserved() -> None:
    """Normal sample index 0 must be forwarded (--normal-sample-index 0)."""
    vc_config = {
        "tumor_sample_index": 1,
        "normal_sample_index": 0,
    }
    flags = build_vc_flags(vc_config, {})
    assert "--tumor-sample-index 1" in flags
    assert "--normal-sample-index 0" in flags


def test_workflow_reference_and_regions_bed_forwarded() -> None:
    """--reference and --regions-bed must be forwarded and quoted."""
    vc_config = {
        "reference": "GRCh38",
        "regions_bed": "/path with spaces/my regions.bed",
    }
    flags = build_vc_flags(vc_config, {})
    assert "--reference GRCh38" in flags
    assert "--regions-bed '/path with spaces/my regions.bed'" in flags


def test_partial_matching_contig_does_not_mask_autosome_mismatch(tmp_path: Path) -> None:
    """A matching contig like MT must NOT mask chr prefix mismatches on autosomes."""
    vcf = tmp_path / "test.vcf.gz"
    bed = tmp_path / "test.bed"

    _write_vcf(vcf, ["chr1", "chr2", "MT"])
    _write_bed(bed, ["1\t100\t200\treg1", "2\t100\t200\treg2", "MT\t10\t50\tregMT"])

    with pytest.raises(ValueError, match="Contig prefix mismatch"):
        validate_bed_vcf_contigs(str(vcf), str(bed))


def test_scaffold_names_distinct() -> None:
    """Scaffold contigs like scaffold_chr1 and scaffold_1 must remain distinct."""
    assert _parse_bed_line("scaffold_chr1\t100\t200\treg")[0] == "scaffold_chr1"
    assert _parse_bed_line("scaffold_1\t100\t200\treg")[0] == "scaffold_1"


def test_extract_variants_raises_on_contig_mismatch(tmp_path: Path) -> None:
    """extract_variants must validate contigs and fail before running bcftools."""
    vcf = tmp_path / "test.vcf.gz"
    bed = tmp_path / "test.bed"
    out_vcf = tmp_path / "out.vcf.gz"

    _write_vcf(vcf, ["chr1", "chr2"])
    _write_bed(bed, ["1\t100\t200\treg1"])

    from variantcentrifuge.filters import extract_variants

    with pytest.raises(ValueError, match="Contig prefix mismatch"):
        extract_variants(str(vcf), str(bed), {"threads": 1}, str(out_vcf))


def test_roundtrip_workflow_flags_through_create_parser() -> None:
    """Flags built by build_vc_flags must round-trip cleanly through create_parser."""
    vc_config = {
        "tumor_sample_index": 0,
        "normal_sample_index": 1,
        "reference": "GRCh38",
        "regions_bed": "/path with spaces/my regions.bed",
    }
    flags = build_vc_flags(vc_config, {})
    cmd_str = " ".join(flags)
    split_args = shlex.split(cmd_str)

    parser = create_parser()
    # Provide required arguments for CLI parsing
    args = parser.parse_args(["--vcf-file", "in.vcf", "--gene-name", "all", *split_args])

    assert args.tumor_sample_index == 0
    assert args.normal_sample_index == 1
    assert args.reference == "GRCh38"
    assert args.regions_bed == "/path with spaces/my regions.bed"
