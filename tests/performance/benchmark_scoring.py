"""
Scoring benchmarks.

Tests variant scoring formula application using pd.eval at multiple scales.
Scoring performance is primarily variant-count sensitive.
"""

import pytest

from variantcentrifuge.scoring import apply_scoring


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000, 10000])
def test_scoring_apply_scaling(benchmark, synthetic_variants, synthetic_scoring_config, n_variants):
    """
    Benchmark scoring formula application at multiple scales.

    Scoring uses pd.eval to apply formulas over DataFrame columns. Performance
    scales with variant count (row operations) rather than sample count.
    """
    # Use moderate sample count (scoring is not sample-count sensitive)
    n_samples = 10

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Get scoring config
    scoring_config = synthetic_scoring_config

    # Count formulas for metadata
    n_formulas = len(scoring_config.get("formulas", {}))

    # Benchmark scoring application
    def run_scoring():
        return apply_scoring(df, scoring_config)

    if n_variants >= 10000:
        result = benchmark.pedantic(run_scoring, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_scoring)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_formulas"] = n_formulas
    benchmark.extra_info["component"] = "scoring"

    # Verify result has new score columns
    assert result is not None
    # formulas is a list of dicts, extract all score names
    for formula_dict in scoring_config["formulas"]:
        for score_name in formula_dict.keys():
            assert score_name in result.columns


@pytest.mark.performance
@pytest.mark.parametrize("n_formulas", [1, 2, 5])
def test_scoring_formula_scaling(benchmark, synthetic_variants, n_formulas):
    """
    Benchmark scoring with varying formula counts.

    Tests whether performance scales linearly with number of formulas or
    has overhead from formula compilation/evaluation setup.
    """
    n_variants = 1000
    n_samples = 10

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Build scoring config with variable formula count
    scoring_config = {
        "variables": {
            "impact_var": {
                "IMPACT": {
                    "HIGH": 3,
                    "MODERATE": 2,
                    "LOW": 1,
                    "MODIFIER": 0,
                }
            },
            "effect_var": {
                "EFFECT": {
                    "stop_gained": 5,
                    "frameshift_variant": 4,
                    "missense_variant": 2,
                    "splice_region_variant": 1,
                    "synonymous_variant": 0,
                }
            },
        },
        "formulas": [],
    }

    # Add requested number of formulas (formulas is a list of dicts)
    for i in range(n_formulas):
        formula_name = f"score_{i}"
        if i == 0:
            scoring_config["formulas"].append({formula_name: "impact_var"})
        elif i == 1:
            scoring_config["formulas"].append({formula_name: "effect_var"})
        else:
            # Combine for additional formulas
            scoring_config["formulas"].append({formula_name: f"impact_var + effect_var * {i}"})

    # Benchmark
    result = benchmark(apply_scoring, df, scoring_config)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_formulas"] = n_formulas
    benchmark.extra_info["component"] = "scoring_formula_scaling"

    # Verify
    assert result is not None
    assert len([col for col in result.columns if col.startswith("score_")]) == n_formulas
