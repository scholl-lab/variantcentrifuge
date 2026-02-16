"""
Pytest-based golden file validation tests.

These tests ensure inheritance analysis produces clinically equivalent output
after vectorization by comparing against golden reference files.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

# Add project root and scripts to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "scripts"))

from scripts.validate_inheritance import get_all_scenarios
from variantcentrifuge.inheritance.analyzer import analyze_inheritance

GOLDEN_DIR = project_root / "tests" / "fixtures" / "golden"

# Extract scenario names for parametrization
ALL_SCENARIOS = get_all_scenarios()
SCENARIO_NAMES = [scenario[3] for scenario in ALL_SCENARIOS]


@pytest.fixture(scope="module")
def check_golden_files_exist():
    """Check if golden files exist, skip tests if not."""
    if not GOLDEN_DIR.exists():
        pytest.skip(
            f"Golden directory does not exist: {GOLDEN_DIR}. "
            "Run 'python scripts/validate_inheritance.py generate' first."
        )

    missing_files = []
    for scenario_name in SCENARIO_NAMES:
        golden_path = GOLDEN_DIR / f"{scenario_name}.parquet"
        if not golden_path.exists():
            missing_files.append(str(golden_path))

    if missing_files:
        pytest.skip(
            f"Golden files missing: {missing_files}. "
            "Run 'python scripts/validate_inheritance.py generate' first."
        )


@pytest.mark.unit
@pytest.mark.inheritance
@pytest.mark.parametrize("scenario_name", SCENARIO_NAMES)
def test_golden_file_match(scenario_name, check_golden_files_exist):
    """Test that inheritance analysis matches golden reference for each scenario."""
    # Get the scenario data
    scenario = next(s for s in ALL_SCENARIOS if s[3] == scenario_name)
    df, pedigree_data, sample_list, _ = scenario

    # Load golden file
    golden_path = GOLDEN_DIR / f"{scenario_name}.parquet"
    golden_df = pd.read_parquet(golden_path)

    # Run analyze_inheritance with current code
    result_df = analyze_inheritance(df, pedigree_data, sample_list)

    # Sort both by stable key (CHROM, POS, REF, ALT)
    sort_cols = ["CHROM", "POS", "REF", "ALT"]
    golden_sorted = golden_df.sort_values(sort_cols).reset_index(drop=True)
    result_sorted = result_df.sort_values(sort_cols).reset_index(drop=True)

    # Assert row count matches
    assert len(golden_sorted) == len(result_sorted), (
        f"Row count mismatch: golden={len(golden_sorted)}, result={len(result_sorted)}"
    )

    # Assert Inheritance_Pattern matches exactly
    assert "Inheritance_Pattern" in golden_sorted.columns, "Golden file missing Inheritance_Pattern"
    assert "Inheritance_Pattern" in result_sorted.columns, "Result missing Inheritance_Pattern"

    for idx in range(len(golden_sorted)):
        golden_pattern = golden_sorted.loc[idx, "Inheritance_Pattern"]
        result_pattern = result_sorted.loc[idx, "Inheritance_Pattern"]

        variant_key = (
            f"{golden_sorted.loc[idx, 'CHROM']}:{golden_sorted.loc[idx, 'POS']}"
            f":{golden_sorted.loc[idx, 'REF']}>{golden_sorted.loc[idx, 'ALT']}"
        )

        assert golden_pattern == result_pattern, (
            f"Pattern mismatch at {variant_key}: golden={golden_pattern}, result={result_pattern}"
        )

    # Assert Inheritance_Details matches (parse JSON and compare key fields)
    assert "Inheritance_Details" in golden_sorted.columns, "Golden file missing Inheritance_Details"
    assert "Inheritance_Details" in result_sorted.columns, "Result missing Inheritance_Details"

    for idx in range(len(golden_sorted)):
        golden_details = json.loads(golden_sorted.loc[idx, "Inheritance_Details"])
        result_details = json.loads(result_sorted.loc[idx, "Inheritance_Details"])

        variant_key = (
            f"{golden_sorted.loc[idx, 'CHROM']}:{golden_sorted.loc[idx, 'POS']}"
            f":{golden_sorted.loc[idx, 'REF']}>{golden_sorted.loc[idx, 'ALT']}"
        )

        # Compare primary_pattern
        assert golden_details.get("primary_pattern") == result_details.get("primary_pattern"), (
            f"primary_pattern mismatch at {variant_key}"
        )

        # Compare all_patterns (as sets for order-independence)
        golden_patterns = set(golden_details.get("all_patterns", []))
        result_patterns = set(result_details.get("all_patterns", []))
        assert golden_patterns == result_patterns, (
            f"all_patterns mismatch at {variant_key}: "
            f"golden={golden_patterns}, result={result_patterns}"
        )

        # Compare confidence (within epsilon=0.001)
        golden_conf = golden_details.get("confidence", 0.0)
        result_conf = result_details.get("confidence", 0.0)
        assert abs(golden_conf - result_conf) <= 0.001, (
            f"confidence mismatch at {variant_key}: "
            f"golden={golden_conf:.3f}, result={result_conf:.3f}"
        )

        # Compare samples_with_pattern (sample IDs should match)
        golden_samples = {s["sample_id"] for s in golden_details.get("samples_with_pattern", [])}
        result_samples = {s["sample_id"] for s in result_details.get("samples_with_pattern", [])}
        assert golden_samples == result_samples, (
            f"samples_with_pattern mismatch at {variant_key}: "
            f"golden={golden_samples}, result={result_samples}"
        )


@pytest.mark.unit
@pytest.mark.inheritance
def test_scenario_determinism():
    """Test that all scenarios produce identical output when run twice."""
    for df, pedigree_data, sample_list, scenario_name in ALL_SCENARIOS:
        # Run twice
        result1 = analyze_inheritance(df.copy(), pedigree_data, sample_list)
        result2 = analyze_inheritance(df.copy(), pedigree_data, sample_list)

        # Sort both by stable key
        sort_cols = ["CHROM", "POS", "REF", "ALT"]
        result1_sorted = result1.sort_values(sort_cols).reset_index(drop=True)
        result2_sorted = result2.sort_values(sort_cols).reset_index(drop=True)

        # Compare Inheritance_Pattern
        assert result1_sorted["Inheritance_Pattern"].equals(
            result2_sorted["Inheritance_Pattern"]
        ), f"Inheritance_Pattern not deterministic for {scenario_name}"

        # Compare Inheritance_Details (JSON parse to handle order differences)
        for idx in range(len(result1_sorted)):
            details1 = json.loads(result1_sorted.loc[idx, "Inheritance_Details"])
            details2 = json.loads(result2_sorted.loc[idx, "Inheritance_Details"])

            # Compare key fields (not full JSON string due to possible key ordering)
            assert details1.get("primary_pattern") == details2.get("primary_pattern"), (
                f"primary_pattern not deterministic for {scenario_name} row {idx}"
            )

            assert set(details1.get("all_patterns", [])) == set(details2.get("all_patterns", [])), (
                f"all_patterns not deterministic for {scenario_name} row {idx}"
            )

            confidence_diff = abs(details1.get("confidence", 0.0) - details2.get("confidence", 0.0))
            assert confidence_diff <= 0.001, (
                f"confidence not deterministic for {scenario_name} row {idx}"
            )


@pytest.mark.unit
@pytest.mark.inheritance
def test_scenario_data_types():
    """Test that all scenarios have expected data types in output."""
    for df, pedigree_data, sample_list, scenario_name in ALL_SCENARIOS:
        result = analyze_inheritance(df, pedigree_data, sample_list)

        # Check required columns exist
        assert "Inheritance_Pattern" in result.columns, (
            f"Missing Inheritance_Pattern in {scenario_name}"
        )
        assert "Inheritance_Details" in result.columns, (
            f"Missing Inheritance_Details in {scenario_name}"
        )

        # Check Inheritance_Pattern is string type
        assert result["Inheritance_Pattern"].dtype == object, (
            f"Inheritance_Pattern should be string type in {scenario_name}"
        )

        # Check Inheritance_Details is valid JSON
        for idx, details_str in enumerate(result["Inheritance_Details"]):
            try:
                details = json.loads(details_str)
                assert isinstance(details, dict), (
                    f"Inheritance_Details should be dict in {scenario_name}"
                )
                assert "primary_pattern" in details, (
                    f"Missing primary_pattern in Inheritance_Details for {scenario_name} row {idx}"
                )
                assert "all_patterns" in details, (
                    f"Missing all_patterns in Inheritance_Details for {scenario_name} row {idx}"
                )
                assert "confidence" in details, (
                    f"Missing confidence in Inheritance_Details for {scenario_name} row {idx}"
                )
            except json.JSONDecodeError:
                pytest.fail(f"Invalid JSON in Inheritance_Details for {scenario_name} row {idx}")


@pytest.mark.unit
@pytest.mark.inheritance
def test_all_scenarios_covered():
    """Test that all expected scenarios are present."""
    expected_scenarios = {
        "trio_denovo",
        "trio_dominant",
        "trio_recessive",
        "trio_denovo_candidate",
        "single_sample",
        "extended_family",
        "x_linked",
        "mitochondrial",
        "compound_het",
        "edge_cases",
    }

    actual_scenarios = set(SCENARIO_NAMES)

    assert expected_scenarios == actual_scenarios, (
        f"Scenario mismatch: expected={expected_scenarios}, actual={actual_scenarios}"
    )


@pytest.mark.unit
@pytest.mark.inheritance
def test_golden_files_readable():
    """Test that all golden files can be read without errors."""
    if not GOLDEN_DIR.exists():
        pytest.skip("Golden directory does not exist")

    for scenario_name in SCENARIO_NAMES:
        golden_path = GOLDEN_DIR / f"{scenario_name}.parquet"
        if not golden_path.exists():
            pytest.skip(f"Golden file missing: {golden_path}")

        # Should be able to read without error
        df = pd.read_parquet(golden_path)

        # Should have expected columns
        assert "CHROM" in df.columns, f"Missing CHROM in {scenario_name}"
        assert "POS" in df.columns, f"Missing POS in {scenario_name}"
        assert "REF" in df.columns, f"Missing REF in {scenario_name}"
        assert "ALT" in df.columns, f"Missing ALT in {scenario_name}"
        assert "Inheritance_Pattern" in df.columns, (
            f"Missing Inheritance_Pattern in {scenario_name}"
        )
        assert "Inheritance_Details" in df.columns, (
            f"Missing Inheritance_Details in {scenario_name}"
        )

        # Should have at least one row
        assert len(df) > 0, f"Golden file is empty: {scenario_name}"
