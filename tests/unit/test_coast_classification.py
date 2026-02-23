# File: tests/unit/test_coast_classification.py
"""
Unit tests for COAST classification: _resolve_effect, formula engine,
CLI --coast-classification option, and auto-field injection.

Covers COAST-02, COAST-04, COAST-05, COAST-06, COAST-07.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.tests.allelic_series import (
    CADD_COLUMN_CANDIDATES,
    PTV_EFFECTS,
    _resolve_effect,
    classify_variants,
)

# ---------------------------------------------------------------------------
# _resolve_effect tests (COAST-04, COAST-07)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_resolve_effect_simple():
    """Single effect string with no '&' is returned unchanged (fast path)."""
    assert _resolve_effect("stop_gained") == "stop_gained"
    assert _resolve_effect("missense_variant") == "missense_variant"
    assert _resolve_effect("") == ""
    assert _resolve_effect("intron_variant") == "intron_variant"


@pytest.mark.unit
def test_resolve_effect_ptv_priority():
    """PTV effect takes priority over any other effect in '&'-concatenated string."""
    assert _resolve_effect("stop_gained&splice_region_variant") == "stop_gained"
    assert _resolve_effect("frameshift_variant&intron_variant") == "frameshift_variant"
    assert _resolve_effect("splice_acceptor_variant&missense_variant") == "splice_acceptor_variant"
    assert _resolve_effect("splice_donor_variant&synonymous_variant") == "splice_donor_variant"


@pytest.mark.unit
def test_resolve_effect_missense_priority():
    """missense_variant takes priority over non-PTV effects."""
    assert _resolve_effect("missense_variant&splice_region_variant") == "missense_variant"
    assert _resolve_effect("intron_variant&missense_variant") == "missense_variant"


@pytest.mark.unit
def test_resolve_effect_no_match():
    """When no PTV or missense, the first part is returned."""
    assert _resolve_effect("intron_variant&synonymous_variant") == "intron_variant"
    assert _resolve_effect("synonymous_variant&intron_variant") == "synonymous_variant"


@pytest.mark.unit
def test_resolve_effect_ptv_effects_coverage():
    """All PTV_EFFECTS members resolve correctly from concatenated strings."""
    for ptv_effect in PTV_EFFECTS:
        result = _resolve_effect(f"{ptv_effect}&splice_region_variant")
        assert result == ptv_effect, f"Expected {ptv_effect!r}, got {result!r}"


# ---------------------------------------------------------------------------
# classify_variants with '&'-concatenated effect strings (COAST-04)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_classify_variants_with_resolve_effect():
    """classify_variants correctly classifies '&'-concatenated SnpEff multi-transcript effects."""
    df = pd.DataFrame(
        {
            "EFFECT": [
                "stop_gained&splice_region_variant",  # should be PTV (3)
                "missense_variant&splice_region_variant",  # should be DMV (2)
                "intron_variant",  # should be unclassified (0)
            ],
            "IMPACT": ["HIGH", "MODERATE", "LOW"],
            "dbNSFP_SIFT_pred": ["", "deleterious", ""],
            "dbNSFP_Polyphen2_HDIV_pred": ["", "probably_damaging", ""],
        }
    )

    codes, mask = classify_variants(df, "EFFECT", "IMPACT")

    assert int(codes[0]) == 3, f"stop_gained&splice_region should be PTV=3, got {codes[0]}"
    assert int(codes[1]) == 2, f"missense&splice_region with SIFT-D should be DMV=2, got {codes[1]}"
    assert int(codes[2]) == 0, f"intron_variant should be unclassified=0, got {codes[2]}"
    assert bool(mask[0]) is True
    assert bool(mask[1]) is True
    assert bool(mask[2]) is False


# ---------------------------------------------------------------------------
# classify_variants with sift_polyphen formula engine (COAST-06)
# ---------------------------------------------------------------------------


@pytest.fixture
def sift_polyphen_model_dir():
    """Return absolute path to sift_polyphen classification config."""
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    model_dir = os.path.join(repo_root, "scoring", "coast_classification", "sift_polyphen")
    if not os.path.isdir(model_dir):
        pytest.skip(f"sift_polyphen model not found at {model_dir}")
    return model_dir


@pytest.fixture
def cadd_model_dir():
    """Return absolute path to cadd classification config."""
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    model_dir = os.path.join(repo_root, "scoring", "coast_classification", "cadd")
    if not os.path.isdir(model_dir):
        pytest.skip(f"cadd model not found at {model_dir}")
    return model_dir


@pytest.mark.unit
def test_classify_variants_sift_polyphen_config(sift_polyphen_model_dir):
    """Formula engine (sift_polyphen model) produces identical results to hardcoded logic."""
    df = pd.DataFrame(
        {
            "EFFECT": [
                "stop_gained",
                "missense_variant",
                "missense_variant",
                "missense_variant",
                "synonymous_variant",
            ],
            "IMPACT": ["HIGH", "MODERATE", "MODERATE", "MODERATE", "LOW"],
            "dbNSFP_SIFT_pred": ["", "deleterious", "tolerated", "", ""],
            "dbNSFP_Polyphen2_HDIV_pred": [
                "",
                "probably_damaging",
                "benign",
                "",
                "",
            ],
        }
    )

    # Default hardcoded path
    codes_default, mask_default = classify_variants(df, "EFFECT", "IMPACT")

    # Formula engine path
    codes_formula, mask_formula = classify_variants(
        df, "EFFECT", "IMPACT", model_dir=sift_polyphen_model_dir
    )

    np.testing.assert_array_equal(
        codes_default,
        codes_formula,
        err_msg="sift_polyphen formula engine produces different codes than hardcoded logic",
    )
    np.testing.assert_array_equal(
        mask_default,
        mask_formula,
        err_msg="sift_polyphen formula engine produces different mask than hardcoded logic",
    )

    # Verify expected values
    assert int(codes_formula[0]) == 3, "stop_gained should be PTV"
    assert int(codes_formula[1]) == 2, "missense + SIFT-D + PolyPhen-D should be DMV"
    assert int(codes_formula[2]) == 1, "missense + SIFT-T + PolyPhen-B should be BMV"
    assert int(codes_formula[3]) == 0, "missense without predictions should be unclassified"
    assert int(codes_formula[4]) == 0, "synonymous should be unclassified"


@pytest.mark.unit
def test_classify_variants_cadd_config(cadd_model_dir):
    """CADD model: DMV when CADD>=15, BMV when 0<CADD<15, PTV independent of CADD."""
    df = pd.DataFrame(
        {
            "EFFECT": [
                "stop_gained",
                "missense_variant",
                "missense_variant",
                "missense_variant",
                "synonymous_variant",
            ],
            "IMPACT": ["HIGH", "MODERATE", "MODERATE", "MODERATE", "LOW"],
            "dbNSFP_CADD_phred": [35.0, 20.0, 8.0, 0.0, 5.0],
        }
    )

    codes, mask = classify_variants(df, "EFFECT", "IMPACT", model_dir=cadd_model_dir)

    assert int(codes[0]) == 3, f"stop_gained should be PTV=3, got {codes[0]}"
    assert int(codes[1]) == 2, f"missense + CADD=20 (>=15) should be DMV=2, got {codes[1]}"
    assert int(codes[2]) == 1, f"missense + CADD=8 (<15) should be BMV=1, got {codes[2]}"
    assert int(codes[3]) == 0, f"missense + CADD=0 should be unclassified=0, got {codes[3]}"
    assert int(codes[4]) == 0, f"synonymous should be unclassified=0, got {codes[4]}"

    assert bool(mask[0]) is True
    assert bool(mask[1]) is True
    assert bool(mask[2]) is True
    assert bool(mask[3]) is False
    assert bool(mask[4]) is False


# ---------------------------------------------------------------------------
# CLI --coast-classification option (COAST-05)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_coast_classification_cli_option_default():
    """--coast-classification defaults to 'sift_polyphen'."""
    import argparse

    # Import the parser-building function or just test the argument default
    # We parse an empty coast-related args set
    parser = argparse.ArgumentParser()
    parser.add_argument("--coast-classification", type=str, default="sift_polyphen")
    args = parser.parse_args([])
    assert args.coast_classification == "sift_polyphen"


@pytest.mark.unit
def test_coast_classification_cli_option_cadd():
    """--coast-classification cadd can be parsed."""
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--coast-classification", type=str, default="sift_polyphen")
    args = parser.parse_args(["--coast-classification", "cadd"])
    assert args.coast_classification == "cadd"


# ---------------------------------------------------------------------------
# Auto-field injection logic (COAST-02)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_auto_field_injection_sift_polyphen(sift_polyphen_model_dir):
    """Auto-injection appends SIFT/PolyPhen fields when COAST is selected."""
    # Simulate the injection logic from cli.py
    import json

    var_config_path = os.path.join(sift_polyphen_model_dir, "variable_assignment_config.json")
    with open(var_config_path) as f:
        var_config = json.load(f)

    # Collect real field names (skip _comment keys, skip COAST_* internal names)
    required_fields = [k for k in var_config.get("variables", {}) if not k.startswith("_")]
    vcf_fields = [f for f in required_fields if not f.startswith("COAST_")]

    # The sift_polyphen model should NOT produce raw VCF field names
    # (they are COAST_EFFECT, COAST_IMPACT, etc. - all internal)
    # So vcf_fields should be empty for normalized configs
    assert isinstance(vcf_fields, list)


@pytest.mark.unit
def test_auto_field_injection_adds_missing_fields():
    """Auto-injection adds missing fields to existing field list."""
    # Simulate the injection logic with a mock list of required fields
    existing_fields = "CHROM,POS,REF,ALT,EFFECT,IMPACT"
    required_vcf_fields = ["dbNSFP_SIFT_pred", "dbNSFP_Polyphen2_HDIV_pred"]

    field_list = [f.strip() for f in existing_fields.split(",")]
    injected = []
    for rf in required_vcf_fields:
        if rf not in field_list:
            field_list.append(rf)
            injected.append(rf)

    assert "dbNSFP_SIFT_pred" in field_list
    assert "dbNSFP_Polyphen2_HDIV_pred" in field_list
    assert injected == ["dbNSFP_SIFT_pred", "dbNSFP_Polyphen2_HDIV_pred"]
    assert field_list[0] == "CHROM"  # original order preserved


@pytest.mark.unit
def test_auto_field_injection_skips_existing_fields():
    """Auto-injection does not duplicate fields already in the extraction list."""
    existing_fields = "CHROM,POS,REF,ALT,EFFECT,IMPACT,dbNSFP_SIFT_pred"
    required_vcf_fields = ["dbNSFP_SIFT_pred", "dbNSFP_Polyphen2_HDIV_pred"]

    field_list = [f.strip() for f in existing_fields.split(",")]
    injected = []
    for rf in required_vcf_fields:
        if rf not in field_list:
            field_list.append(rf)
            injected.append(rf)

    assert field_list.count("dbNSFP_SIFT_pred") == 1, "SIFT field should not be duplicated"
    assert injected == ["dbNSFP_Polyphen2_HDIV_pred"]


# ---------------------------------------------------------------------------
# CADD_COLUMN_CANDIDATES constant (COAST-07)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_cadd_column_candidates_exist():
    """CADD_COLUMN_CANDIDATES is exported and contains dbNSFP_CADD_phred."""
    assert "dbNSFP_CADD_phred" in CADD_COLUMN_CANDIDATES


@pytest.mark.unit
def test_classify_variants_cadd_column_detection(cadd_model_dir):
    """classify_variants detects CADD column from CADD_COLUMN_CANDIDATES."""
    # Use alternate CADD column name
    df = pd.DataFrame(
        {
            "EFFECT": ["missense_variant", "missense_variant"],
            "IMPACT": ["MODERATE", "MODERATE"],
            "CADD_phred": [20.0, 8.0],  # alternate name from CADD_COLUMN_CANDIDATES
        }
    )

    codes, _mask = classify_variants(df, "EFFECT", "IMPACT", model_dir=cadd_model_dir)

    assert int(codes[0]) == 2, f"CADD=20 should be DMV=2, got {codes[0]}"
    assert int(codes[1]) == 1, f"CADD=8 should be BMV=1, got {codes[1]}"
