# Pure Python COAST Implementation Plan

**Author:** Senior Data Scientist / Software Engineer / Mathematician assessment
**Date:** 2026-02-22
**Status:** Draft — ready for GSD phase integration

## Problem Statement

The COAST (COding-variant Allelic Series Test) is currently implemented as an R wrapper via rpy2, calling `AllelicSeries::COAST()`. This creates:

1. **R dependency** — inconsistent with SKAT Python backend philosophy (Phase 21)
2. **Thread-safety limitation** — `parallel_safe=False` due to rpy2 GIL
3. **Deployment friction** — requires R + AllelicSeries CRAN package installed
4. **Testing complexity** — requires elaborate rpy2 mocking infrastructure

## Mathematical Analysis

COAST (McCaw et al., AJHG 2023) decomposes into **three computational layers**, all of which we already have infrastructure for.

### Layer 1: Burden Component (6 tests)

COAST tests **three models** × **two encodings** = 6 burden tests:

| Model | Encoding | Degrees of Freedom | Formula |
|-------|----------|-------------------|---------|
| Baseline | Count | 3 | `E(Y) = α + N₁β₁ + N₂β₂ + N₃β₃ + X'γ` |
| Baseline | Indicator | 3 | `E(Y) = α + I(N₁>0)β₁ + I(N₂>0)β₂ + I(N₃>0)β₃ + X'γ` |
| Sum | Count | 1 | `E(Y) = α + (w₁N₁ + w₂N₂ + w₃N₃)β + X'γ` |
| Sum | Indicator | 1 | `E(Y) = α + (w₁I₁ + w₂I₂ + w₃I₃)β + X'γ` |
| Max | Count | 1 | `E(Y) = α + max(w₁N₁, w₂N₂, w₃N₃)β + X'γ` |
| Max | Indicator | 1 | `E(Y) = α + max(w₁I₁, w₂I₂, w₃I₃)β + X'γ` |

Where:
- `Nₗ = Σⱼ Gⱼ` for variants in category l (allele count)
- `Iₗ = I(Nₗ > 0)` (carrier indicator)
- `w = (1, 2, 3)` default weights for BMV, DMV, PTV

**Test statistics:**
- Continuous traits: Wald chi-squared from OLS regression
- Binary traits: Likelihood ratio test from logistic regression (deviance difference)

**What we already have:** `statsmodels.OLS` and `statsmodels.Logit` are used in `linear_burden.py` and `logistic_burden.py`. The regression infrastructure is proven.

### Layer 2: Allelic Series SKAT (1 test)

Standard SKAT score test with **annotation-aware weights**:

```
w_SKAT,j = sqrt(w_{A_j} / var(G_j))
```

Where `w_{A_j}` is the allelic series weight for variant j's annotation category. This encodes: rarer + more deleterious variants → larger expected effects.

The test statistic `Q = Σⱼ w²ⱼ Sⱼ²` follows a mixture of chi-squareds under H₀. P-values via Davies → saddlepoint → Liu fallback chain.

**What we already have:** `PythonSKATBackend._test_skat()` implements exactly this. The only change is substituting annotation-aware weights for Beta(MAF) weights.

### Layer 3: Cauchy Omnibus Combination

Combine 7 component p-values (6 burden + 1 SKAT) via ACAT:

```
T = Σₘ wₘ tan(π(0.5 - pₘ))
p_omni = 0.5 - arctan(T / Σwₘ) / π
```

Default weights: each burden test gets weight 1, SKAT gets weight 6 (achieving 50/50 burden vs SKAT evidence split).

**What we already have:** `cauchy_combination()` in `association/acat.py` — proven, handles edge cases, tiny-p approximation.

## Existing Infrastructure Reuse

| Component | Source | Status |
|-----------|--------|--------|
| Davies p-value (qfc.c CFFI) | `backends/davies.py` | Production, validated against R |
| SKAT score test | `backends/python_backend.py::_test_skat()` | Production, validated against R |
| Eigenvalue decomposition | `backends/python_backend.py::_get_lambda()` | Production |
| Null model fitting (GLM) | `backends/python_backend.py::fit_null_model()` | Production |
| Cauchy combination | `tests/acat.py::cauchy_combination()` | Production |
| Variant classification | `tests/allelic_series.py::classify_variants()` | Production, already pure Python |
| OLS regression | `tests/linear_burden.py` | Production |
| Logistic regression | `tests/logistic_burden.py` | Production |
| Beta(MAF) weights | `weights.py::beta_maf_weights()` | Production |

## Implementation Plan

### Task 1: PythonCOASTBackend class (~150 lines)

New file: `variantcentrifuge/association/backends/coast_python.py`

```python
class PythonCOASTBackend:
    """Pure Python COAST implementation matching R AllelicSeries::COAST()."""

    def fit_null_model(self, phenotype, covariates, trait_type):
        """Reuse PythonSKATBackend.fit_null_model() — identical null model."""

    def _aggregate_by_category(self, geno, anno, weights):
        """Compute N_l (count) and I_l (indicator) per category.
        Returns: count_matrix (n_samples × 3), indicator_matrix (n_samples × 3)
        """

    def _burden_baseline(self, agg_matrix, null_model, trait_type):
        """3-df chi-squared test: all three category coefficients jointly zero."""

    def _burden_sum(self, agg_matrix, weights, null_model, trait_type):
        """1-df test: weighted sum predictor."""

    def _burden_max(self, agg_matrix, weights, null_model, trait_type):
        """1-df test: max weighted predictor."""

    def _allelic_skat(self, geno, anno, weights, null_model):
        """SKAT with annotation-aware weights: w_j = sqrt(w_{A_j} / var(G_j))."""

    def test_gene(self, gene, geno, anno, phenotype, covariates, weights, trait_type):
        """Run all 7 components, combine via Cauchy, return omnibus p-value."""
```

### Task 2: PurePythonCOASTTest wrapper (~100 lines)

New file or extend: `variantcentrifuge/association/tests/allelic_series_python.py`

- Subclass `AssociationTest` (like `PurePythonSKATTest` wraps `PythonSKATBackend`)
- `parallel_safe = True` (no rpy2!)
- Reuses existing `classify_variants()` from `allelic_series.py`
- Delegates to `PythonCOASTBackend.test_gene()`
- Returns `TestResult` with omnibus p-value + component p-values in `extra`

### Task 3: Engine/Registry integration (~30 lines)

- Register `PurePythonCOASTTest` in `_build_registry()` as `"coast"` when R unavailable
- Add backend selection logic: `--coast-backend python|r|auto` (mirroring SKAT pattern)
- Default: `auto` (try R first, fall back to Python)

### Task 4: Validation test suite (~200 lines)

- **Unit tests:** Each burden model variant against manual statsmodels reference
- **SKAT component:** Against existing `PythonSKATBackend` with same weights
- **Cauchy combination:** Already tested in `test_acat.py`
- **Integration:** Compare Python vs R COAST p-values on synthetic data (tiered tolerance, same as SKAT validation)
- **Edge cases:** Missing categories, single variant per category, all-zero genotypes

### Task 5: Cleanup and documentation (~20 lines)

- Update `COASTTest` to attempt Python backend before R
- Update CLI `--help` text
- Log backend choice at INFO level

## Estimated Complexity

| Task | New Lines | Difficulty |
|------|-----------|-----------|
| PythonCOASTBackend | ~150 | Moderate (regression + SKAT reuse) |
| PurePythonCOASTTest | ~100 | Easy (wrapper pattern proven) |
| Registry integration | ~30 | Easy |
| Validation tests | ~200 | Moderate (reference values) |
| Cleanup | ~20 | Easy |
| **Total** | **~500** | **Moderate overall** |

## Risk Assessment

| Risk | Mitigation |
|------|-----------|
| Numerical divergence from R | Tiered tolerance (same as SKAT: <1e-4 relative for p>1e-4) |
| Max model non-differentiability | Use indicator-based Wald test (not score); max() is piecewise linear |
| 3-df baseline test p-value method | Use chi2.sf on Wald statistic (identical to R's anova/LRT) |
| Annotation-aware SKAT weight edge cases | Guard div-by-zero when var(G_j)=0 (monomorphic variant → skip) |

## Success Criteria

1. `--coast-backend python` produces p-values within tiered tolerance of R `AllelicSeries::COAST()`
2. `--coast-backend python` works without R installed
3. `parallel_safe=True` on PurePythonCOASTTest (thread-safe)
4. All existing COAST tests pass with Python backend
5. No regression in other test results (Fisher, burden, SKAT, ACAT-O)

## References

- McCaw et al. (2023) "An allelic-series rare-variant association test for candidate-gene discovery" AJHG
- Liu & Xie (2020) Cauchy combination test
- Lee et al. (2012) SKAT-O optimal unified test
- Wu et al. (2011) SKAT rare variant association test
