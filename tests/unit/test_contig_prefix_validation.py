"""Regression tests for VCF contig extraction and prefix validation (#151)."""

from pathlib import Path

import pytest

from variantcentrifuge.filters import validate_bed_vcf_contigs


def test_contig_fallback_scan_detects_mismatch_beyond_1000_records(tmp_path: Path) -> None:
    """Verify fallback body scan does not truncate at 1,000 records on headerless-contig VCF."""
    vcf = tmp_path / "headerless_contigs.vcf"
    bed = tmp_path / "test.bed"

    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for i in range(1001):
        lines.append(f"chr1\t{i + 1}\t.\tA\tT\t100\tPASS\t.")
    lines.append("chr2\t100\t.\tA\tT\t100\tPASS\t.")
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")

    bed.write_text("2\t90\t110\tregion_2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Contig prefix mismatch"):
        validate_bed_vcf_contigs(str(vcf), str(bed))


def test_partial_contig_header_scan_detects_mismatch(tmp_path: Path) -> None:
    """Verify that a VCF declaring only chrM in header but having data records on chr1
    triggers contig prefix mismatch when BED uses '1'."""
    vcf = tmp_path / "partial_header_contigs.vcf"
    bed = tmp_path / "test.bed"

    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chrM,length=16569>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t100\t.\tA\tT\t100\tPASS\t.",
    ]
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bed.write_text("1\t90\t110\tregion_1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Contig prefix mismatch"):
        validate_bed_vcf_contigs(str(vcf), str(bed))


def test_partial_contig_header_fallback_scan_beyond_100_records(tmp_path: Path) -> None:
    """Verify that a VCF declaring only chrM in header and having >100 chrM records before chr1
    triggers contig prefix mismatch via full records scan when BED uses '1'."""
    vcf = tmp_path / "partial_header_many_records.vcf"
    bed = tmp_path / "test.bed"

    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chrM,length=16569>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for i in range(120):
        lines.append(f"chrM\t{i + 1}\t.\tA\tT\t100\tPASS\t.")
    lines.append("chr1\t100\t.\tA\tT\t100\tPASS\t.")
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bed.write_text("1\t90\t110\tregion_1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Contig prefix mismatch"):
        validate_bed_vcf_contigs(str(vcf), str(bed))


def test_mixed_order_does_not_falsely_mismatch_chr_to_plain(tmp_path: Path) -> None:
    """Verify 101 records on chr1 followed by record 102 on 1 with BED interval on 1
    does not raise ValueError.
    """
    vcf = tmp_path / "chr1_then_1.vcf"
    bed = tmp_path / "test.bed"

    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for i in range(101):
        lines.append(f"chr1\t{i + 1}\t.\tA\tT\t100\tPASS\t.")
    lines.append("1\t150\t.\tA\tT\t100\tPASS\t.")
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bed.write_text("1\t140\t160\tregion_1\n", encoding="utf-8")

    # Should not raise ValueError because 1 is present in VCF at record 102
    validate_bed_vcf_contigs(str(vcf), str(bed))


def test_mixed_order_does_not_falsely_mismatch_plain_to_chr(tmp_path: Path) -> None:
    """Verify 101 records on 1 followed by record 102 on chr1 with BED interval on chr1
    does not raise ValueError.
    """
    vcf = tmp_path / "1_then_chr1.vcf"
    bed = tmp_path / "test.bed"

    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for i in range(101):
        lines.append(f"1\t{i + 1}\t.\tA\tT\t100\tPASS\t.")
    lines.append("chr1\t150\t.\tA\tT\t100\tPASS\t.")
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bed.write_text("chr1\t140\t160\tregion_1\n", encoding="utf-8")

    # Should not raise ValueError because chr1 is present in VCF at record 102
    validate_bed_vcf_contigs(str(vcf), str(bed))


def test_validate_bed_vcf_contigs_reads_at_most_100_records_with_matched_and_unmatched_contig(
    tmp_path: Path,
) -> None:
    """validate_bed_vcf_contigs must read at most 100 data records even when BED contains
    an unmatched contig (e.g. chrMT) alongside a matched contig (chr1).
    """
    from unittest.mock import patch

    from variantcentrifuge.filters import _VCF_CONTIG_CACHE, _VCF_CONTIG_FULL_SCANNED

    vcf = tmp_path / "large_cohort.vcf"
    bed = tmp_path / "genes_with_mt.bed"

    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=248956422>",
        "##contig=<ID=chrM,length=16569>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for i in range(500):
        lines.append(f"chr1\t{i + 1}\t.\tA\tT\t100\tPASS\t.")
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bed.write_text("chr1\t10\t200\tGENE1\nchrMT\t1\t100\tMT_GENE\n", encoding="utf-8")

    _VCF_CONTIG_CACHE.pop(str(vcf), None)
    _VCF_CONTIG_FULL_SCANNED.discard(str(vcf))

    real_open = open
    data_lines_read = 0

    class CountingReader:
        def __init__(self, fh):
            self._fh = fh

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return self._fh.__exit__(*args)

        def __iter__(self):
            return self

        def __next__(self):
            nonlocal data_lines_read
            line = next(self._fh)
            if not line.startswith("#"):
                data_lines_read += 1
            return line

    def wrapped_open(file, *args, **kwargs):
        fh = real_open(file, *args, **kwargs)
        if str(file) == str(vcf):
            return CountingReader(fh)
        return fh

    with patch("builtins.open", side_effect=wrapped_open):
        validate_bed_vcf_contigs(str(vcf), str(bed))

    assert data_lines_read <= 100
    assert str(vcf) not in _VCF_CONTIG_FULL_SCANNED
