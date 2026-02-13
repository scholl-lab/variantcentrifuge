"""Unit tests for variantcentrifuge.clinvar_pm5."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from variantcentrifuge.clinvar_pm5 import (
    AA_3TO1,
    annotate_pm5,
    check_pm5,
    load_pm5_lookup,
    normalize_aa,
    parse_hgvs_p,
)


# ---------------------------------------------------------------------------
# parse_hgvs_p
# ---------------------------------------------------------------------------
class TestParseHgvsP:
    def test_3letter_standard(self) -> None:
        assert parse_hgvs_p("p.Arg123His") == ("Arg", 123, "His")

    def test_1letter_standard(self) -> None:
        assert parse_hgvs_p("p.R123H") == ("R", 123, "H")

    def test_3letter_parenthesized(self) -> None:
        assert parse_hgvs_p("p.(Gly480Cys)") == ("Gly", 480, "Cys")

    def test_1letter_parenthesized(self) -> None:
        assert parse_hgvs_p("p.(G480C)") == ("G", 480, "C")

    def test_nonsense_3letter(self) -> None:
        assert parse_hgvs_p("p.Arg123*") == ("Arg", 123, "*")

    def test_nonsense_1letter(self) -> None:
        assert parse_hgvs_p("p.R123*") == ("R", 123, "*")

    def test_unknown_effect(self) -> None:
        assert parse_hgvs_p("p.?") is None

    def test_synonymous(self) -> None:
        assert parse_hgvs_p("p.=") is None

    def test_na(self) -> None:
        assert parse_hgvs_p("NA") is None

    def test_dot(self) -> None:
        assert parse_hgvs_p(".") is None

    def test_empty_string(self) -> None:
        assert parse_hgvs_p("") is None

    def test_embedded_in_longer_string(self) -> None:
        # Like ClinVar Name field
        result = parse_hgvs_p("NM_000059.4(BRCA2):c.8167G>C (p.Asp2723His)")
        assert result == ("Asp", 2723, "His")

    def test_frameshift_returns_none(self) -> None:
        assert parse_hgvs_p("p.Arg123Profs*45") is None


# ---------------------------------------------------------------------------
# normalize_aa
# ---------------------------------------------------------------------------
class TestNormalizeAA:
    def test_3letter_to_1letter(self) -> None:
        assert normalize_aa("Arg") == "R"
        assert normalize_aa("His") == "H"
        assert normalize_aa("Ter") == "*"

    def test_1letter_passthrough(self) -> None:
        assert normalize_aa("R") == "R"
        assert normalize_aa("H") == "H"
        assert normalize_aa("*") == "*"

    def test_all_standard_amino_acids(self) -> None:
        for three, one in AA_3TO1.items():
            assert normalize_aa(three) == one


# ---------------------------------------------------------------------------
# check_pm5
# ---------------------------------------------------------------------------
SAMPLE_LOOKUP: dict[str, list[dict[str, str]]] = {
    "BRCA2:2723": [
        {"ref_aa": "Asp", "alt_aa": "His", "clnsig": "P", "review": "expert_panel"},
        {"ref_aa": "Asp", "alt_aa": "Gly", "clnsig": "LP", "review": "single_submitter"},
    ],
    "TP53:248": [
        {"ref_aa": "Arg", "alt_aa": "Trp", "clnsig": "P", "review": "expert_panel"},
        {"ref_aa": "Arg", "alt_aa": "Gln", "clnsig": "P", "review": "expert_panel"},
        {"ref_aa": "Arg", "alt_aa": "Leu", "clnsig": "P", "review": "multiple_submitters"},
    ],
    "SCN5A:1795": [
        {"ref_aa": "Arg", "alt_aa": "His", "clnsig": "LP", "review": "single_submitter"},
    ],
}


class TestCheckPM5:
    def test_pm5_triggered_one_pathogenic(self) -> None:
        # Query: BRCA2 Asp2723Ala (different from His and Gly)
        result = check_pm5("BRCA2", "Asp", 2723, "Ala", SAMPLE_LOOKUP)
        assert result["triggered"] is True
        assert result["strength"] == "PM5"
        assert result["evidence_count"] == 2  # His(P) + Gly(LP)

    def test_pm5_strong_multiple_pathogenic(self) -> None:
        # Query: TP53 Arg248Pro (different from Trp, Gln, Leu)
        result = check_pm5("TP53", "Arg", 248, "Pro", SAMPLE_LOOKUP)
        assert result["triggered"] is True
        assert result["strength"] == "PM5_Strong"
        assert result["evidence_count"] == 3

    def test_pm5_supporting_only_lp(self) -> None:
        # Query: SCN5A Arg1795Cys (different from His which is LP)
        result = check_pm5("SCN5A", "Arg", 1795, "Cys", SAMPLE_LOOKUP)
        assert result["triggered"] is True
        assert result["strength"] == "PM5_Supporting"
        assert result["evidence_count"] == 1

    def test_ps1_not_pm5_same_aa_change(self) -> None:
        # Query: BRCA2 Asp2723His — same change as known pathogenic → PS1, not PM5
        result = check_pm5("BRCA2", "Asp", 2723, "His", SAMPLE_LOOKUP)
        # Should still trigger PM5 from the Gly(LP) entry
        assert result["triggered"] is True
        assert result["evidence_count"] == 1  # only Gly(LP), His excluded
        assert result["strength"] == "PM5_Supporting"

    def test_not_triggered_no_match(self) -> None:
        result = check_pm5("NONEXISTENT", "Arg", 999, "His", SAMPLE_LOOKUP)
        assert result["triggered"] is False
        assert result["strength"] is None

    def test_not_triggered_all_same_change(self) -> None:
        # All known entries are the query change
        lookup = {
            "GENE1:100": [
                {"ref_aa": "Arg", "alt_aa": "His", "clnsig": "P", "review": "ep"},
            ]
        }
        result = check_pm5("GENE1", "Arg", 100, "His", lookup)
        assert result["triggered"] is False

    def test_ref_aa_mismatch_skipped(self) -> None:
        # Lookup has Asp at 2723, but query has Glu → ref mismatch
        result = check_pm5("BRCA2", "Glu", 2723, "Ala", SAMPLE_LOOKUP)
        assert result["triggered"] is False

    def test_mixed_codes_normalized(self) -> None:
        # Query uses 1-letter, lookup uses 3-letter
        result = check_pm5("BRCA2", "D", 2723, "A", SAMPLE_LOOKUP)
        assert result["triggered"] is True
        assert result["evidence_count"] == 2

    def test_known_variants_string(self) -> None:
        result = check_pm5("TP53", "Arg", 248, "Pro", SAMPLE_LOOKUP)
        known = result["known_variants"]
        assert "Arg248Trp(P)" in known
        assert "Arg248Gln(P)" in known
        assert "Arg248Leu(P)" in known


# ---------------------------------------------------------------------------
# load_pm5_lookup
# ---------------------------------------------------------------------------
class TestLoadPM5Lookup:
    def test_load_json(self, tmp_path: Path) -> None:
        meta = {
            "position_count": 1,
            "entry_count": 1,
            "assembly": "GRCh38",
            "build_date": "2026-01-01",
        }
        entry = {"ref_aa": "Arg", "alt_aa": "His", "clnsig": "P", "review": "ep"}
        data = {"_meta": meta, "lookup": {"GENE1:100": [entry]}}
        p = tmp_path / "pm5.json"
        p.write_text(json.dumps(data))
        result = load_pm5_lookup(str(p))
        assert "GENE1:100" in result

    def test_load_gzipped_json(self, tmp_path: Path) -> None:
        import gzip

        meta = {
            "position_count": 1,
            "entry_count": 1,
            "assembly": "GRCh38",
            "build_date": "2026-01-01",
        }
        entry = {"ref_aa": "Arg", "alt_aa": "His", "clnsig": "P", "review": "ep"}
        data = {"_meta": meta, "lookup": {"GENE1:100": [entry]}}
        p = tmp_path / "pm5.json.gz"
        with gzip.open(str(p), "wt", encoding="utf-8") as f:
            json.dump(data, f)
        result = load_pm5_lookup(str(p))
        assert "GENE1:100" in result


# ---------------------------------------------------------------------------
# annotate_pm5
# ---------------------------------------------------------------------------
class TestAnnotatePM5:
    def test_adds_three_columns(self) -> None:
        df = pd.DataFrame(
            {
                "GENE": ["BRCA2", "TP53", "UNKNOWN"],
                "HGVS_P": ["p.Asp2723Ala", "p.Arg248Pro", "p.?"],
            }
        )
        result = annotate_pm5(df, SAMPLE_LOOKUP)
        assert "PM5" in result.columns
        assert "PM5_evidence_count" in result.columns
        assert "PM5_known_variants" in result.columns

    def test_correct_values(self) -> None:
        df = pd.DataFrame(
            {
                "GENE": ["BRCA2", "TP53", "UNKNOWN"],
                "HGVS_P": ["p.Asp2723Ala", "p.Arg248Pro", "p.?"],
            }
        )
        result = annotate_pm5(df, SAMPLE_LOOKUP)
        assert result.iloc[0]["PM5"] == "PM5"
        assert result.iloc[1]["PM5"] == "PM5_Strong"
        assert result.iloc[2]["PM5"] == ""

    def test_empty_dataframe(self) -> None:
        df = pd.DataFrame({"GENE": [], "HGVS_P": []})
        result = annotate_pm5(df, SAMPLE_LOOKUP)
        assert len(result) == 0
        assert "PM5" in result.columns

    def test_missing_hgvs_column(self) -> None:
        df = pd.DataFrame(
            {
                "GENE": ["BRCA2"],
                "HGVS_P": ["NA"],
            }
        )
        result = annotate_pm5(df, SAMPLE_LOOKUP)
        assert result.iloc[0]["PM5"] == ""

    def test_preserves_existing_columns(self) -> None:
        df = pd.DataFrame(
            {
                "GENE": ["BRCA2"],
                "HGVS_P": ["p.Asp2723Ala"],
                "CHROM": ["13"],
                "POS": [32914438],
            }
        )
        result = annotate_pm5(df, SAMPLE_LOOKUP)
        assert "CHROM" in result.columns
        assert "POS" in result.columns
        assert result.iloc[0]["CHROM"] == "13"
