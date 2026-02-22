# Phase 24: Pure Python COAST Backend - Research

**Researched:** 2026-02-22
**Domain:** Pure Python implementation of R AllelicSeries::COAST() allelic series test
**Confidence:** HIGH (all key findings verified directly from source code)

---

## Summary

Phase 24 implements a pure Python COAST backend that mirrors the R `AllelicSeries::COAST()`
function without requiring rpy2 or R. All three computational layers of COAST are already
implemented in the codebase — the implementation is primarily assembly work, not new
algorithms. The SKAT backend pattern (Phase 21) provides the exact structural template to
follow. The mathematical plan in `.planning/PURE_PYTHON_COAST_PLAN.md` is accurate and ready
to execute.

**Primary recommendation:** Follow the SKAT backend pattern exactly. Create
`coast_python.py` (backend), `allelic_series_python.py` (test wrapper), update `engine.py`
and `base.py`, and add `--coast-backend` to `cli.py`. Reuse `PythonSKATBackend`,
`cauchy_combination()`, OLS/logistic regression infrastructure — all already proven against R.

---

## Standard Stack

### Core (all already installed — no new dependencies)

| Library | Version | Purpose | Why Used |
|---------|---------|---------|----------|
| numpy | >=1.24 | Array ops, matrix algebra | Already a core dependency |
| scipy | >=1.10 | `chi2.sf()`, `norm.sf()` for p-values | Already used by SKAT backend |
| statsmodels | >=0.14 | GLM null model, OLS, Logit | Already used by burden tests |

### Infrastructure Already in Codebase

| Module | Location | What It Provides |
|--------|----------|-----------------|
| `PythonSKATBackend` | `backends/python_backend.py` | `fit_null_model()`, `_test_skat()`, `_compute_eigenvalues_filtered()` |
| `cauchy_combination()` | `tests/acat.py` | Cauchy/ACAT combination with edge-case handling |
| `classify_variants()` | `tests/allelic_series.py` | BMV/DMV/PTV annotation from EFFECT/IMPACT columns |
| `beta_maf_weights()` | `weights.py` | Beta(MAF) weight computation |
| `AssociationTest` | `base.py` | Abstract base class for all tests |
| `AssociationConfig` | `base.py` | Config dataclass (add `coast_backend` field here) |
| Davies p-value chain | `backends/davies.py` | `compute_pvalue()`: Davies → saddlepoint → Liu |

**No new libraries required.** All dependencies are already in `pyproject.toml`.

---

## Architecture Patterns

### Pattern 1: SKAT Backend Split (proven in Phase 20-21)

COAST follows the identical split used by SKAT:

```
tests/allelic_series.py        ← COASTTest (R backend, existing)
tests/allelic_series_python.py ← PurePythonCOASTTest (new, wraps backend)
backends/coast_python.py        ← PythonCOASTBackend (new, pure math)
engine.py                       ← backend-aware swap at from_names() time
```

### Pattern 2: Backend-Aware Registry Swap (existing in engine.py)

```python
# engine.py from_names() — EXISTING SKAT PATTERN
skat_backend = getattr(config, "skat_backend", "r")
if skat_backend == "python":
    registry["skat"] = PurePythonSKATTest
elif skat_backend == "auto":
    try:
        get_skat_backend("r")
    except Exception:
        registry["skat"] = PurePythonSKATTest

# NEW: identical pattern for COAST
coast_backend = getattr(config, "coast_backend", "auto")
if coast_backend == "python":
    registry["coast"] = PurePythonCOASTTest
elif coast_backend == "auto":
    try:
        import rpy2.robjects  # noqa: F401
        importr("AllelicSeries")
    except Exception:
        registry["coast"] = PurePythonCOASTTest
```

### Pattern 3: PurePythonCOASTTest class structure

Directly mirrors `PurePythonSKATTest` from `tests/skat_python.py`:

```python
class PurePythonCOASTTest(AssociationTest):
    parallel_safe: bool = True  # Thread-safe — no rpy2

    @property
    def name(self) -> str:
        return "coast"  # Same name as COASTTest — engine swaps implementation

    def check_dependencies(self) -> None:
        # Just verify numpy/scipy/statsmodels (always present)
        pass

    def effect_column_names(self) -> dict[str, str | None]:
        return {"effect": None, "se": None, "ci_lower": None, "ci_upper": None}

    def prepare(self, gene_count: int) -> None:
        # Identical to COASTTest.prepare() — just no R gc()

    def run(self, gene, contingency_data, config) -> TestResult:
        # 1. Classify variants (reuse classify_variants())
        # 2. Filter to COAST-eligible variants
        # 3. Fit null model (lazy, via PythonCOASTBackend)
        # 4. Run PythonCOASTBackend.test_gene()
        # 5. Return TestResult with omnibus p-value + component p-values in extra
```

### Pattern 4: COAST Mathematical Decomposition (verified from AllelicSeries source)

COAST produces **7 component p-values** combined via ACAT with specific weights:

```
6 burden p-values:
  base_count   — Baseline count encoding (3-df Wald/LRT)
  base_ind     — Baseline indicator encoding (3-df Wald/LRT)
  max_count    — Max allele count across categories (1-df)
  max_ind      — Max indicator across categories (1-df)
  sum_count    — Weighted sum of allele counts (1-df)
  sum_ind      — Weighted sum of indicators (1-df)

1 SKAT p-value:
  allelic_skat — Standard SKAT with annotation-aware weights

Cauchy weights:
  Each burden test: weight = 1
  SKAT test: weight = 6  (= n_burden / n_skat = 6/1, achieving 50/50 split)
```

**Verified from AllelicSeries R source:**
```r
omni_weights <- c(rep(1, n_burden), rep(n_burden / n_skat, n_skat))
# n_burden = 6, n_skat = 1 → SKAT weight = 6
```

### Pattern 5: Burden Aggregation

Each burden test aggregates genotype matrix `G` (n_samples × n_variants) by category:

```python
def _aggregate_by_category(geno, anno_codes, weights, indicator=False, method="none"):
    """
    Aggregate variants by COAST category (BMV=1, DMV=2, PTV=3).

    Parameters
    ----------
    geno : np.ndarray, shape (n_samples, n_variants_filtered)
        Genotype dosage for COAST-eligible variants only.
    anno_codes : np.ndarray, shape (n_variants_filtered,)
        Annotation codes 1/2/3 for BMV/DMV/PTV.
    weights : list[float]
        Category weights (default [1.0, 2.0, 3.0]).
    indicator : bool
        If True, convert allele counts to 0/1 presence indicators.
    method : str
        "none" → return category-level aggregates (3-column matrix)
        "sum"  → return weighted sum across categories (1-column vector)
        "max"  → return max weighted value across categories (1-column vector)

    Returns
    -------
    np.ndarray, shape (n_samples, 1) or (n_samples, 3)
        Aggregated burden predictor(s).
    """
    n_samples = geno.shape[0]
    category_agg = np.zeros((n_samples, 3), dtype=np.float64)

    for i, cat in enumerate([1, 2, 3]):  # BMV, DMV, PTV
        cat_mask = (anno_codes == cat)
        if cat_mask.any():
            cat_sum = geno[:, cat_mask].sum(axis=1)  # allele count per sample
            if indicator:
                cat_sum = (cat_sum > 0).astype(np.float64)
            category_agg[:, i] = cat_sum * weights[i]

    if method == "none":
        return category_agg  # (n_samples, 3) — for baseline 3-df test
    elif method == "sum":
        return category_agg.sum(axis=1, keepdims=True)  # (n_samples, 1)
    elif method == "max":
        return category_agg.max(axis=1, keepdims=True)  # (n_samples, 1)
```

**Key distinction for baseline model:**
- `method="none"` returns a 3-column matrix → baseline uses all 3 columns simultaneously
- The baseline burden test is a multi-df joint test (H₀: β₁ = β₂ = β₃ = 0)
- P-value: chi-squared Wald test on 3 parameters from OLS/Logit

### Pattern 6: Allelic Series SKAT Weights (verified from AllelicSeries R source)

```python
def _compute_allelic_skat_weights(geno_filtered, anno_codes_filtered, coast_weights):
    """
    Annotation-aware SKAT weights for allelic series test.

    From AllelicSeries source:
        aaf = apply(geno, 2, mean) / 2   # alternate allele frequency
        w[anno == i] = weights[i]         # assign category weight
        v = aaf * (1 - aaf)               # allele freq variance
        skat_weights = sqrt(w / v)        # annotation-aware weight

    Parameters
    ----------
    geno_filtered : np.ndarray, shape (n_samples, n_variants_filtered)
    anno_codes_filtered : np.ndarray, shape (n_variants_filtered,) — values 1/2/3
    coast_weights : list[float] — category weights [w_bmv, w_dmv, w_ptv]

    Returns
    -------
    np.ndarray, shape (n_variants_filtered,)
        SKAT weight per variant.
    """
    # Alternate allele frequency per variant
    aaf = geno_filtered.mean(axis=0) / 2.0  # (n_variants,)
    aaf = np.clip(aaf, 1e-8, 1 - 1e-8)

    # Category weight per variant
    w = np.zeros(len(anno_codes_filtered), dtype=np.float64)
    for i, cat in enumerate([1, 2, 3]):  # BMV, DMV, PTV
        w[anno_codes_filtered == cat] = coast_weights[i]

    # Variance of genotype (for Hardy-Weinberg: var(G) = 2 * p * (1-p))
    v = aaf * (1.0 - aaf)  # Note: R uses p*(1-p) not 2p(1-p) — see below

    # Annotation-aware SKAT weight
    skat_weights = np.sqrt(w / v)  # (n_variants,)
    return skat_weights
```

**Important:** The R source uses `aaf * (1 - aaf)` for variance (not `2 * aaf * (1-aaf)`). This is the convention used by the AllelicSeries package. Use the same convention.

### Pattern 7: Burden Test P-value Computation

```python
def _run_burden_test(predictor, phenotype, covariates, null_model, trait_type):
    """
    Run OLS (continuous) or Logit LRT (binary) burden test.

    Returns p-value as float or None on failure.
    """
    import statsmodels.api as sm

    design = sm.add_constant(predictor)  # predictor is (n_samples, 1) or (n_samples, 3)
    if covariates is not None:
        design = np.column_stack([design, covariates])

    if trait_type == "quantitative":
        # OLS Wald chi-squared on burden coefficient(s)
        fit = sm.OLS(phenotype, design).fit()
        # For 1-df: p = fit.pvalues[1]
        # For 3-df: Wald chi-sq test H0: β₁=β₂=β₃=0
        if predictor.shape[1] == 1:
            return float(fit.pvalues[1])
        else:
            # 3-df joint Wald test: fit.f_test(R) where R selects burden coefficients
            r_matrix = np.eye(3, design.shape[1], k=1)  # rows selecting β1, β2, β3
            f_test = fit.f_test(r_matrix)
            return float(f_test.pvalue)
    else:
        # Binary: LRT = deviance difference between null and full model
        # Null model: intercept + covariates (no burden)
        null_design = design[:, [0]]  # intercept only
        if covariates is not None:
            null_design = np.column_stack([null_design, covariates])

        null_fit = sm.Logit(phenotype, null_design).fit(disp=False)
        full_fit = sm.Logit(phenotype, design).fit(disp=False)

        lrt = 2 * (full_fit.llf - null_fit.llf)
        df = predictor.shape[1]  # 1 or 3
        return float(scipy.stats.chi2.sf(lrt, df=df))
```

### Pattern 8: AssociationConfig Extension

Add `coast_backend` field to `AssociationConfig` in `base.py`:

```python
# Phase 24: COAST backend field
coast_backend: str = "auto"
"""COAST computation backend: "r" (R via rpy2), "python", or "auto" (try r first)."""
```

### Recommended Project Structure (new files)

```
variantcentrifuge/association/
├── backends/
│   ├── coast_python.py          ← NEW: PythonCOASTBackend (pure math layer)
│   └── __init__.py              ← No change needed (coast backend not factory-ized)
├── tests/
│   ├── allelic_series.py        ← EXISTING: COASTTest (R backend)
│   └── allelic_series_python.py ← NEW: PurePythonCOASTTest (wraps backend)
├── base.py                      ← ADD: coast_backend field to AssociationConfig
└── engine.py                    ← ADD: coast backend-aware swap logic

variantcentrifuge/
└── cli.py                       ← ADD: --coast-backend argument

tests/unit/
├── test_coast.py                ← EXISTING: classify_variants + COASTTest tests
└── test_coast_python.py         ← NEW: PurePythonCOASTTest + validation tests
```

### Anti-Patterns to Avoid

- **Don't rebuild null model per gene.** Fit once per cohort (same as `PurePythonSKATTest`), store as `self._null_model`.
- **Don't import the COAST backend lazily inside `run()`.** Initialize in `check_dependencies()`, store as `self._backend`.
- **Don't use `parallel_safe = False` on the Python test.** The whole point is `parallel_safe = True`.
- **Don't call `_test_skat()` directly from `PythonSKATBackend`.** The allelic SKAT needs custom annotation weights, not Beta(MAF) weights. Inline the SKAT score test with the new weight vector.
- **Don't use the 3-df joint Wald test via `f_test()` for logistic regression.** Use LRT (deviance difference) for binary traits — that's what AllelicSeries does.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Cauchy combination | Custom ACAT | `cauchy_combination()` in `tests/acat.py` | Already handles tiny-p approximation, None filtering, weight normalization |
| GLM null model | Custom GLM | `PythonSKATBackend.fit_null_model()` | Already validated against R SKAT null model |
| Davies p-value chain | Custom eigenvalue p-value | `compute_pvalue()` in `backends/davies.py` | Three-tier Davies→saddlepoint→Liu, already validated |
| OLS regression | Custom OLS | `statsmodels.OLS` | Already used in `linear_burden.py` |
| Logistic regression | Custom Logit | `statsmodels.Logit` | Already used in `logistic_burden.py` |
| Variant classification | Custom EFFECT/IMPACT parser | `classify_variants()` in `allelic_series.py` | Already handles multi-value fields, auto-detect columns |
| Beta(MAF) weights | Custom MAF weight | `beta_maf_weights()` in `weights.py` | For validation; COAST uses annotation weights, not MAF weights for SKAT |
| Eigenvalue computation | Custom eigh | `scipy.linalg.eigh(driver='evr')` | Matches R's DSYEVR, already used by SKAT backend |

**Key insight:** The pure Python COAST backend is an assembly task. All primitives exist and are validated. The new code connects them using the AllelicSeries::COAST() algorithm specification.

---

## Common Pitfalls

### Pitfall 1: Allelic SKAT Weight Convention
**What goes wrong:** Using `2 * aaf * (1-aaf)` for variance (diploid convention) instead of `aaf * (1-aaf)` (R AllelicSeries convention).
**Why it happens:** Hardy-Weinberg variance of diploid genotype G ∈ {0,1,2} is `2p(1-p)`, but the R source uses `p*(1-p)` directly.
**How to avoid:** Use `v = aaf * (1 - aaf)` exactly as in AllelicSeries source. This is a constant factor difference that shifts all SKAT p-values — matching R output requires matching this convention.
**Warning signs:** SKAT component p-values diverge systematically (not just numerical noise) from R.

### Pitfall 2: Null Model for Burden Test vs SKAT
**What goes wrong:** Using the wrong null model for burden tests. The burden tests use OLS/Logit fitted *with the burden predictor* as full model and *without* as null model (for LRT). The SKAT null model is fitted without any genotype predictors.
**Why it happens:** The SKAT null model residuals are not the right residuals for burden LRT.
**How to avoid:** Fit a fresh null model for each burden test configuration (or use the shared null model from `PythonSKATBackend.fit_null_model()` only for the SKAT component, not for burden p-values).
**Correction:** Actually — for burden tests, use standard OLS/Logit directly (no pre-fitted null model needed). For SKAT, reuse `PythonSKATBackend.fit_null_model()`.

### Pitfall 3: Baseline Test Degrees of Freedom
**What goes wrong:** Running the baseline burden test as 1-df instead of 3-df.
**Why it happens:** The baseline uses all three category columns simultaneously, so it's a joint test.
**How to avoid:** For the baseline burden test (method="none"), pass all 3 category columns as predictors and test their joint significance (chi2 with df=3). This requires the Wald F-test (quantitative) or LRT with df=3 (binary).

### Pitfall 4: Max Predictor Gradient
**What goes wrong:** Trying to compute a gradient-based test statistic for the max predictor.
**Why it happens:** `max()` is not differentiable, but regression still works because you're regressing on the max value, not differentiating through it.
**How to avoid:** The max predictor `max(w₁N₁, w₂N₂, w₃N₃)` per sample is just a scalar computed before fitting. Fit standard OLS/Logit on this scalar predictor — 1 burden coefficient, 1-df test.

### Pitfall 5: Category Weights in Baseline Test
**What goes wrong:** Applying category weights to the baseline model (they should use UNWEIGHTED category aggregates for the baseline).
**Why it happens:** Misreading the AllelicSeries algorithm.
**How to avoid:** Baseline uses `weights = unif_weights` (all 1s) — the three columns are raw per-category allele counts/indicators without the 1/2/3 scaling. Only sum and max models apply the user-specified weights.

**Verified from AllelicSeries R source:**
```r
burden_results$base <- BurdenWrap(indicator = FALSE, method = "none", weights = unif_weights)
burden_results$ind  <- BurdenWrap(indicator = TRUE,  method = "none", weights = unif_weights)
# Note: unif_weights = c(1, 1, 1), NOT the user's coast_weights
```

### Pitfall 6: Cauchy Combination Weight Ordering
**What goes wrong:** Passing SKAT weight before burden weights, or using wrong proportionality.
**Why it happens:** The ordering must match the p-value array order passed to `cauchy_combination()`.
**How to avoid:** Follow the pattern: collect 6 burden p-values first, then SKAT p-value last. Pass `weights=[1,1,1,1,1,1,6]` (6 burden with weight 1, 1 SKAT with weight 6).

### Pitfall 7: Engine Backend Swap Order
**What goes wrong:** The `coast_backend` check in `engine.py` must come BEFORE the unknown-name check.
**Why it happens:** Same pitfall as SKAT (already documented in existing code comments).
**How to avoid:** Add coast backend swap logic immediately after SKAT backend swap logic in `engine.py::from_names()`.

---

## Code Examples

### Example 1: Full PythonCOASTBackend.test_gene() skeleton

```python
# Source: .planning/PURE_PYTHON_COAST_PLAN.md + AllelicSeries R source verification

def test_gene(self, gene, geno_filtered, anno_codes_filtered, phenotype, covariates,
              coast_weights, trait_type):
    """
    Run all 7 COAST components and return omnibus p-value.

    Returns dict with keys: p_value, burden_p_values (list of 6), skat_p_value
    """
    from variantcentrifuge.association.tests.acat import cauchy_combination

    # Ensure null model is fitted
    if self._null_model is None:
        self._null_model = self._skat_backend.fit_null_model(phenotype, covariates, trait_type)

    # Compute per-sample category aggregates (for burden)
    # Baseline: unweighted (weights=[1,1,1]), count and indicator
    cat_count_uniform = _aggregate_by_category(geno_filtered, anno_codes_filtered,
                                               [1.0, 1.0, 1.0], indicator=False, method="none")
    cat_ind_uniform = _aggregate_by_category(geno_filtered, anno_codes_filtered,
                                             [1.0, 1.0, 1.0], indicator=True, method="none")

    # Sum/Max: with coast_weights (e.g. [1,2,3])
    sum_count = _aggregate_by_category(geno_filtered, anno_codes_filtered,
                                       coast_weights, indicator=False, method="sum")
    sum_ind   = _aggregate_by_category(geno_filtered, anno_codes_filtered,
                                       coast_weights, indicator=True,  method="sum")
    max_count = _aggregate_by_category(geno_filtered, anno_codes_filtered,
                                       coast_weights, indicator=False, method="max")
    max_ind   = _aggregate_by_category(geno_filtered, anno_codes_filtered,
                                       coast_weights, indicator=True,  method="max")

    # Run 6 burden tests
    p_base_count = _run_burden_test(cat_count_uniform, phenotype, covariates, trait_type, df=3)
    p_base_ind   = _run_burden_test(cat_ind_uniform,   phenotype, covariates, trait_type, df=3)
    p_sum_count  = _run_burden_test(sum_count,         phenotype, covariates, trait_type, df=1)
    p_sum_ind    = _run_burden_test(sum_ind,            phenotype, covariates, trait_type, df=1)
    p_max_count  = _run_burden_test(max_count,         phenotype, covariates, trait_type, df=1)
    p_max_ind    = _run_burden_test(max_ind,            phenotype, covariates, trait_type, df=1)

    burden_pvals = [p_base_count, p_base_ind, p_sum_count, p_sum_ind, p_max_count, p_max_ind]

    # Run allelic SKAT with annotation-aware weights
    skat_weights = _compute_allelic_skat_weights(geno_filtered, anno_codes_filtered,
                                                  coast_weights)
    p_skat = _run_allelic_skat(geno_filtered, skat_weights, self._null_model)

    # Cauchy combination: 6 burden (weight=1) + 1 SKAT (weight=6)
    all_pvals = burden_pvals + [p_skat]
    cauchy_weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 6.0]
    p_omni = cauchy_combination(all_pvals, weights=cauchy_weights)

    return {
        "p_value": p_omni,
        "burden_p_values": burden_pvals,
        "skat_p_value": p_skat,
    }
```

### Example 2: Allelic SKAT with custom weights (inline, not via PythonSKATBackend.test_gene)

```python
# Source: PythonSKATBackend._test_skat() adapted with annotation weights

def _run_allelic_skat(geno_filtered, skat_weights, null_model):
    """
    SKAT score test with annotation-aware weights (not Beta(MAF) weights).

    skat_weights : np.ndarray shape (n_variants,) — pre-computed annotation weights
    """
    residuals = null_model.extra["residuals"]

    # Apply custom weights (replaces beta_maf_weights call in _test_skat)
    geno_weighted = geno_filtered * skat_weights[np.newaxis, :]

    # Score statistic Q = (score_vec @ score_vec) / 2
    score_vec = geno_weighted.T @ residuals
    q_stat = float(score_vec @ score_vec) / 2.0

    # Eigenvalues (reuse projection-adjusted approach from PythonSKATBackend)
    lambdas = skat_backend._compute_eigenvalues_filtered(geno_weighted, null_model)

    if len(lambdas) == 0:
        return 1.0

    # Davies → saddlepoint → Liu p-value
    p_value, _, _ = compute_pvalue(q_stat, lambdas)
    return float(p_value)
```

### Example 3: Engine registry addition (in engine.py)

```python
# After SKAT backend swap, add COAST backend swap:
coast_backend = getattr(config, "coast_backend", "auto")
if coast_backend == "python":
    from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest
    registry["coast"] = PurePythonCOASTTest
elif coast_backend == "auto":
    # Try R; fall back to Python
    try:
        import rpy2.robjects  # noqa: F401
        from rpy2.robjects.packages import importr
        importr("AllelicSeries")
    except Exception:
        from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest
        registry["coast"] = PurePythonCOASTTest
        logger.info("R/AllelicSeries unavailable; using Python COAST backend (auto mode)")
```

### Example 4: CLI argument (in cli.py, after --skat-backend)

```python
stats_group.add_argument(
    "--coast-backend",
    choices=["auto", "r", "python"],
    default="auto",
    help="COAST computation backend: auto (prefer R, fall back to Python), r, or python",
)
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| R-only COAST via rpy2 | Python-first with R fallback | Phase 24 | `parallel_safe=True`, no R dependency |
| `parallel_safe=False` on COAST | `parallel_safe=True` on Python backend | Phase 24 | COAST can run in parallel gene analysis |

**No deprecated approaches** — this is purely additive. The R backend (`COASTTest`) remains in the registry as the fallback and for `--coast-backend r`.

---

## Open Questions

1. **Null model sharing between burden and SKAT components**
   - What we know: `PurePythonSKATTest` fits the null model once per cohort, cached as `self._null_model`
   - What's unclear: The SKAT component uses the GLM null model residuals. Burden tests do NOT use this null model — they fit fresh OLS/Logit per test. This is slightly wasteful but correct.
   - Recommendation: Do NOT share the null model across burden tests. The SKAT null model extracts residuals for the score test; burden tests need full regression, not just residuals.

2. **Indicator aggregation for baseline model when weights=[1,1,1]**
   - What we know: `BurdenWrap(method="none", weights=unif_weights)` uses uniform weights
   - What's unclear: Does `weights=[1,1,1]` in the baseline change anything vs not applying weights?
   - Recommendation: Since unif_weights = [1,1,1], `category_agg[:, i] *= 1.0` is a no-op. Baseline is just raw per-category sums/indicators.

3. **Tolerance validation setup**
   - What we know: Success criterion requires < 1e-4 relative for p > 1e-4; < 0.5 on log10 for p ≤ 1e-4
   - What's unclear: Need synthetic data where R is available to generate ground-truth p-values
   - Recommendation: Create a tolerance validation test that is `pytest.mark.slow` and skipped when R is unavailable. Generate synthetic data in the test, run both R and Python COAST, compare.

---

## Sources

### Primary (HIGH confidence)

- **Codebase: `tests/skat_python.py`** — `PurePythonSKATTest` class (exact structural template)
- **Codebase: `tests/allelic_series.py`** — `COASTTest`, `classify_variants()` (full implementation read)
- **Codebase: `backends/python_backend.py`** — `PythonSKATBackend` including `fit_null_model()`, `_test_skat()`, `_compute_eigenvalues_filtered()` (full implementation read)
- **Codebase: `tests/acat.py`** — `cauchy_combination()` with weights support (full implementation read)
- **Codebase: `tests/linear_burden.py`** — OLS regression pattern (full implementation read)
- **Codebase: `tests/logistic_burden.py`** — Logit LRT pattern (full implementation read)
- **Codebase: `weights.py`** — `beta_maf_weights()` (full implementation read)
- **Codebase: `base.py`** — `AssociationConfig`, `AssociationTest` ABC (full implementation read)
- **Codebase: `engine.py`** — `AssociationEngine.from_names()` backend swap pattern (full implementation read)
- **Codebase: `backends/__init__.py`** — `get_skat_backend()` factory (full implementation read)
- **AllelicSeries R source via WebFetch** — Exact burden model specification, SKAT weight formula, Cauchy weights (verified from raw GitHub source)

### Secondary (MEDIUM confidence)

- [CRAN AllelicSeries vignette](https://cran.r-project.org/web/packages/AllelicSeries/vignettes/coast.html) — Algorithm overview, p-value component names: base, ind, max_count, max_ind, sum_count, sum_ind, allelic_skat, omni
- [McCaw et al. AJHG 2023 PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10432147/) — Original paper describing COAST
- [GitHub insitro/AllelicSeries](https://github.com/insitro/AllelicSeries) — Repository structure confirmation

---

## Metadata

**Confidence breakdown:**

| Area | Level | Reason |
|------|-------|--------|
| Standard stack | HIGH | All dependencies already in codebase, no new deps needed |
| Architecture | HIGH | SKAT backend pattern read directly from source; fully verified |
| COAST algorithm (burden models) | HIGH | AllelicSeries R source fetched and parsed directly |
| COAST algorithm (SKAT weights) | HIGH | Exact formula from AllelicSeries R source: `sqrt(w / v)` where `v = aaf*(1-aaf)` |
| COAST algorithm (Cauchy weights) | HIGH | Exact formula from AllelicSeries R source: `c(rep(1,6), 6)` |
| Pitfalls | HIGH | Derived from direct code analysis and R source comparison |

**Research date:** 2026-02-22
**Valid until:** 2026-04-22 (stable domain; AllelicSeries algorithm unlikely to change)
