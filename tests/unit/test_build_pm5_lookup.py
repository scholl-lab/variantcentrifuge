"""Unit tests for variantcentrifuge.build_pm5_lookup."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from variantcentrifuge.build_pm5_lookup import (
    _review_stars,
    _simplify_clnsig,
    build_pm5_lookup,
    create_parser,
    write_lookup,
)

# ---------------------------------------------------------------------------
# Fixtures: minimal variant_summary.txt.gz
# ---------------------------------------------------------------------------
HEADER = (
    "#AlleleID\tType\tName\tGeneID\tGeneSymbol\tHGNC_ID\t"
    "ClinicalSignificance\tClinSigSimple\tLastEvaluated\t"
    "RS# (dbSNP)\tnsv/esv (dbVar)\tRCVaccession\tPhenotypeIDS\t"
    "PhenotypeList\tOrigin\tOriginSimple\tAssembly\t"
    "ChromosomeAccession\tChromosome\tStart\tStop\tReferenceAllele\t"
    "AlternateAllele\tCytogenetic\tReviewStatus\tNumberSubmitters\t"
    "Guidelines\tTestedInGTR\tOtherIDs\tSubmitterCategories\t"
    "VariationID\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF\t"
    "RS# (dbSNP2)\tGeneSymbol2\tHGNC_ID2\tAlleleID2\t"
    "ChromosomeAccession2\tVariationType"
)


def _row(*parts: str) -> str:
    """Join tab-delimited ClinVar row parts into a single line."""
    return "\t".join(parts)


ROWS = [
    # Row 1: Valid pathogenic missense, GRCh38, germline
    _row(
        "1",
        "single nucleotide variant",
        "NM_000059.4(BRCA2):c.8167G>C (p.Asp2723His)",
        "675",
        "BRCA2",
        "HGNC:1101",
        "Pathogenic",
        "1",
        "2020-01-01",
        "rs123",
        "-",
        "RCV001",
        "MedGen:C01",
        "Breast cancer",
        "germline",
        "germline",
        "GRCh38",
        "NC_000013.11",
        "13",
        "32914438",
        "32914438",
        "G",
        "C",
        "13q13.1",
        "criteria provided, multiple submitters, no conflicts",
        "3",
        "-",
        "N",
        "-",
        "1",
        "100",
        "32914438",
        "G",
        "C",
        "rs123",
        "BRCA2",
        "HGNC:1101",
        "1",
        "NC_000013.11",
        "SNV",
    ),
    # Row 2: Another pathogenic missense at same position, different AA
    _row(
        "2",
        "single nucleotide variant",
        "NM_000059.4(BRCA2):c.8168A>G (p.Asp2723Gly)",
        "675",
        "BRCA2",
        "HGNC:1101",
        "Likely_pathogenic",
        "1",
        "2020-01-01",
        "rs456",
        "-",
        "RCV002",
        "MedGen:C01",
        "Breast cancer",
        "germline",
        "germline",
        "GRCh38",
        "NC_000013.11",
        "13",
        "32914439",
        "32914439",
        "A",
        "G",
        "13q13.1",
        "criteria provided, single submitter",
        "1",
        "-",
        "N",
        "-",
        "1",
        "101",
        "32914439",
        "A",
        "G",
        "rs456",
        "BRCA2",
        "HGNC:1101",
        "2",
        "NC_000013.11",
        "SNV",
    ),
    # Row 3: TP53 pathogenic
    _row(
        "3",
        "single nucleotide variant",
        "NM_000546.6(TP53):c.742C>T (p.Arg248Trp)",
        "7157",
        "TP53",
        "HGNC:11998",
        "Pathogenic",
        "1",
        "2020-01-01",
        "rs789",
        "-",
        "RCV003",
        "MedGen:C02",
        "Li-Fraumeni",
        "germline",
        "germline",
        "GRCh38",
        "NC_000017.11",
        "17",
        "7674220",
        "7674220",
        "C",
        "T",
        "17p13.1",
        "reviewed by expert panel",
        "5",
        "-",
        "N",
        "-",
        "1",
        "102",
        "7674220",
        "C",
        "T",
        "rs789",
        "TP53",
        "HGNC:11998",
        "3",
        "NC_000017.11",
        "SNV",
    ),
    # Row 4: Should be excluded - Conflicting
    _row(
        "4",
        "single nucleotide variant",
        "NM_000546.6(TP53):c.743G>A (p.Arg248Gln)",
        "7157",
        "TP53",
        "HGNC:11998",
        "Conflicting classifications of pathogenicity",
        "0",
        "2020-01-01",
        "rs101",
        "-",
        "RCV004",
        "MedGen:C02",
        "Li-Fraumeni",
        "germline",
        "germline",
        "GRCh38",
        "NC_000017.11",
        "17",
        "7674221",
        "7674221",
        "G",
        "A",
        "17p13.1",
        "criteria provided, conflicting classifications",
        "2",
        "-",
        "N",
        "-",
        "1",
        "103",
        "7674221",
        "G",
        "A",
        "rs101",
        "TP53",
        "HGNC:11998",
        "4",
        "NC_000017.11",
        "SNV",
    ),
    # Row 5: Should be excluded - not germline (somatic only)
    _row(
        "5",
        "single nucleotide variant",
        "NM_000059.4(BRCA2):c.8100T>A (p.Ser2700Arg)",
        "675",
        "BRCA2",
        "HGNC:1101",
        "Pathogenic",
        "1",
        "2020-01-01",
        "rs102",
        "-",
        "RCV005",
        "MedGen:C01",
        "Cancer",
        "somatic",
        "somatic",
        "GRCh38",
        "NC_000013.11",
        "13",
        "32914370",
        "32914370",
        "T",
        "A",
        "13q13.1",
        "criteria provided, single submitter",
        "1",
        "-",
        "N",
        "-",
        "1",
        "104",
        "32914370",
        "T",
        "A",
        "rs102",
        "BRCA2",
        "HGNC:1101",
        "5",
        "NC_000013.11",
        "SNV",
    ),
    # Row 6: Should be excluded - GRCh37
    _row(
        "6",
        "single nucleotide variant",
        "NM_000059.4(BRCA2):c.8167G>C (p.Asp2723His)",
        "675",
        "BRCA2",
        "HGNC:1101",
        "Pathogenic",
        "1",
        "2020-01-01",
        "rs123",
        "-",
        "RCV001",
        "MedGen:C01",
        "Breast cancer",
        "germline",
        "germline",
        "GRCh37",
        "NC_000013.10",
        "13",
        "32930639",
        "32930639",
        "G",
        "C",
        "13q13.1",
        "criteria provided, multiple submitters, no conflicts",
        "3",
        "-",
        "N",
        "-",
        "1",
        "100",
        "32930639",
        "G",
        "C",
        "rs123",
        "BRCA2",
        "HGNC:1101",
        "1",
        "NC_000013.10",
        "SNV",
    ),
    # Row 7: No review stars (should be excluded at min_stars=1)
    _row(
        "7",
        "single nucleotide variant",
        "NM_000000.1(TESTGENE):c.100A>G (p.Met1Val)",
        "999",
        "TESTGENE",
        "HGNC:999",
        "Pathogenic",
        "1",
        "2020-01-01",
        "rs999",
        "-",
        "RCV006",
        "MedGen:C03",
        "Test",
        "germline",
        "germline",
        "GRCh38",
        "NC_000001.11",
        "1",
        "100",
        "100",
        "A",
        "G",
        "1p36.33",
        "no assertion criteria provided",
        "1",
        "-",
        "N",
        "-",
        "1",
        "105",
        "100",
        "A",
        "G",
        "rs999",
        "TESTGENE",
        "HGNC:999",
        "7",
        "NC_000001.11",
        "SNV",
    ),
]


@pytest.fixture()
def fixture_variant_summary(tmp_path: Path) -> str:
    """Create a minimal variant_summary.txt.gz fixture."""
    p = tmp_path / "variant_summary.txt.gz"
    content = HEADER + "\n" + "\n".join(ROWS) + "\n"
    with gzip.open(str(p), "wt", encoding="utf-8") as f:
        f.write(content)
    return str(p)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
class TestReviewStars:
    def test_known_values(self) -> None:
        assert _review_stars("criteria provided, single submitter") == 1
        assert _review_stars("reviewed by expert panel") == 3
        assert _review_stars("practice guideline") == 4
        assert _review_stars("no assertion criteria provided") == 0

    def test_unknown_defaults_to_zero(self) -> None:
        assert _review_stars("something_unknown") == 0


class TestSimplifyClnsig:
    def test_pathogenic(self) -> None:
        assert _simplify_clnsig("Pathogenic") == "P"

    def test_likely_pathogenic(self) -> None:
        assert _simplify_clnsig("Likely_pathogenic") == "LP"

    def test_pathogenic_likely_pathogenic(self) -> None:
        assert _simplify_clnsig("Pathogenic/Likely_pathogenic") == "P"


# ---------------------------------------------------------------------------
# build_pm5_lookup
# ---------------------------------------------------------------------------
class TestBuildPM5Lookup:
    def test_basic_build(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38")
        assert "_meta" in data
        assert "lookup" in data
        assert data["_meta"]["assembly"] == "GRCh38"

    def test_correct_filtering(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38")
        meta = data["_meta"]
        # Should have rows 1, 2, 3 (row 4=conflicting, 5=somatic, 6=GRCh37, 7=0-star excluded)
        assert meta["rows_after_filter"] == 3

    def test_lookup_contents(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38")
        lookup = data["lookup"]
        assert "BRCA2:2723" in lookup
        assert "TP53:248" in lookup

    def test_brca2_entries(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38")
        entries = data["lookup"]["BRCA2:2723"]
        assert len(entries) == 2
        alt_aas = {e["alt_aa"] for e in entries}
        assert "His" in alt_aas
        assert "Gly" in alt_aas

    def test_deduplication(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38")
        # Row 1 and row 6 are the same variant but different assemblies
        # Only row 1 (GRCh38) should be included, so no duplicates
        entries = data["lookup"]["BRCA2:2723"]
        his_entries = [e for e in entries if e["alt_aa"] == "His"]
        assert len(his_entries) == 1

    def test_grch37_filter(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh37")
        lookup = data["lookup"]
        # Only row 6 is GRCh37
        assert "BRCA2:2723" in lookup
        assert "TP53:248" not in lookup

    def test_exclude_lp(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38", include_lp=False)
        meta = data["_meta"]
        # Row 2 (LP) should be excluded
        assert meta["rows_after_filter"] == 2

    def test_min_stars_3(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38", min_stars=3)
        # Only row 3 (reviewed by expert panel, 3 stars) passes
        assert data["_meta"]["rows_after_filter"] == 1

    def test_meta_counts(self, fixture_variant_summary: str) -> None:
        data = build_pm5_lookup(fixture_variant_summary, assembly="GRCh38")
        meta = data["_meta"]
        assert meta["position_count"] == 2  # BRCA2:2723 + TP53:248
        assert meta["entry_count"] == 3  # 2 at BRCA2 + 1 at TP53


# ---------------------------------------------------------------------------
# write_lookup
# ---------------------------------------------------------------------------
class TestWriteLookup:
    def _sample_data(self) -> dict:
        entry = {"ref_aa": "R", "alt_aa": "H", "clnsig": "P", "review": "ep"}
        return {"_meta": {}, "lookup": {"A:1": [entry]}}

    def test_write_json(self, tmp_path: Path) -> None:
        data = self._sample_data()
        p = tmp_path / "out.json"
        write_lookup(data, str(p))
        loaded = json.loads(p.read_text())
        assert "A:1" in loaded["lookup"]

    def test_write_gzipped_json(self, tmp_path: Path) -> None:
        data = self._sample_data()
        p = tmp_path / "out.json.gz"
        write_lookup(data, str(p))
        with gzip.open(str(p), "rt", encoding="utf-8") as f:
            loaded = json.load(f)
        assert "A:1" in loaded["lookup"]


# ---------------------------------------------------------------------------
# CLI parser
# ---------------------------------------------------------------------------
class TestBuildPM5Parser:
    def test_requires_output(self) -> None:
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--download"])

    def test_download_flag(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["--download", "--output", "pm5.json"])
        assert args.download is True

    def test_input_and_assembly(self) -> None:
        parser = create_parser()
        args = parser.parse_args(
            [
                "--input",
                "variant_summary.txt.gz",
                "--assembly",
                "GRCh37",
                "--output",
                "pm5.json",
            ]
        )
        assert args.input == "variant_summary.txt.gz"
        assert args.assembly == "GRCh37"

    def test_include_lp_default_true(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["--download", "--output", "pm5.json"])
        assert args.include_lp is True

    def test_no_include_lp(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["--download", "--output", "pm5.json", "--no-include-lp"])
        assert args.include_lp is False
