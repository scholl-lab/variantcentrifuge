# Phase 22: ACAT-O and Diagnostics - Research

**Researched:** 2026-02-21
**Domain:** Statistical omnibus testing (ACAT), genomic diagnostics (lambda_GC, QQ plots), Python/scipy
**Confidence:** HIGH

## Summary

Phase 22 adds three capabilities to the association framework: (1) ACAT-V and ACAT-O omnibus p-value computation per gene, (2) a single FDR correction applied to ACAT-O p-values only (not per-test), and (3) a diagnostics module producing lambda_GC per test and QQ plot data TSVs.

ACAT (Aggregated Cauchy Association Test) is a well-defined algorithm from Liu & Xie (2020, AJHG). The formula is exact and simple: transform each p-value via `tan((0.5 - p) * pi)`, sum weighted transforms, read off a Cauchy survival probability. The only edge case requiring special handling is p < 1e-16 (use `1/(p*pi)` instead of `tan` to avoid floating-point overflow). ACAT-O in this codebase means: combine burden_p + skat_p + (optionally acat_v_p) with equal weights using the Cauchy combination formula, then apply a single FDR pass across all genes on the ACAT-O p-values.

Lambda_GC and QQ data are standard genomics diagnostics computed entirely from the p-value vectors already produced by the engine — no new statistical computation is needed, just numpy/scipy operations on existing results.

**Primary recommendation:** Implement `association/tests/acat.py` as an `AssociationTest` subclass that reads pre-computed p-values from `contingency_data` (not from the engine's per-gene loop), plus a standalone `association/diagnostics.py` module with functions `compute_lambda_gc()`, `compute_qq_data()`, and `emit_sample_size_warnings()`.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| numpy | >=1.24 (already in deps) | Array operations for ACAT, lambda_GC, QQ | Already used throughout codebase |
| scipy.stats | >=1.10 (already in deps) | `cauchy.sf()` for ACAT p-value; `chi2.ppf()` for lambda_GC | Already used in genotype_matrix.py and backends |
| statsmodels | >=0.14 (already in deps) | `multipletests` for FDR in correction.py | Already used in correction.py |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pathlib.Path | stdlib | Diagnostics output directory | Always — matches existing codebase style |
| logging | stdlib | Sample size warnings, debug output | Always — matches existing logging pattern |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| scipy.stats.cauchy.sf | 0.5 - arctan(T/w)/pi | Both are identical (verified numerically); use cauchy.sf — clearer intent |
| chi2.ppf(1-p, df=1) | chi2.isf(p, df=1) | Both equivalent; use chi2.isf — more direct |

**Installation:** No new dependencies. All required libraries already in pyproject.toml.

## Architecture Patterns

### How ACAT Fits in the Engine

ACAT-O is NOT a primary test that computes its own statistic from raw genotype data. It is a **meta-test** that combines p-values already computed by other tests (fisher, logistic_burden, linear_burden, skat). This means:

1. The engine runs all primary tests in the per-gene loop (as today).
2. After the per-gene loop (but before FDR), the engine (or a post-processing step) computes ACAT-O per gene by reading existing p-values from `results_by_test`.
3. A single FDR pass is applied to ACAT-O p-values only (ARCH-03).

**Two design options:**

**Option A (recommended): ACAT as AssociationTest that reads from contingency_data**
The engine passes pre-computed p-values from prior tests into `contingency_data` for each gene before the ACAT test runs. Since the engine processes tests in the order they are registered, if `acat_o` is registered after `burden` and `skat`, those p-values can be injected into contingency_data by the engine or stage before calling `acat_o.run()`.

**Option B: Post-loop computation in engine.run_all()**
Compute ACAT-O after the main gene loop in a dedicated post-processing block in `AssociationEngine.run_all()`, directly reading from `results_by_test`.

Option B is simpler and cleaner architecturally — ACAT-O is fundamentally different from primary tests (it has no contingency_data input). The engine already has post-loop logic (FDR correction). Adding ACAT-O as a post-loop step in `run_all()` avoids shoehorning a meta-test into the AssociationTest ABC.

**Recommendation: Use Option B.** Add `_compute_acat_o()` as a private method on `AssociationEngine`, called from `run_all()` after the gene loop, before existing FDR correction. ACAT-O results are added as a synthetic TestResult per gene in `results_by_test["acat_o"]`. FDR is then applied to `acat_o` results only (replacing the current per-test FDR).

### Recommended Project Structure
```
variantcentrifuge/association/
├── tests/
│   ├── acat.py           # ACAT-V + ACAT-O Cauchy combination + AcatTest (if needed as ABC)
├── diagnostics.py        # compute_lambda_gc(), compute_qq_data(), emit_sample_size_warnings()
├── engine.py             # Add _compute_acat_o() + modify FDR strategy (ARCH-03)
└── base.py               # No changes needed
```

### Pattern 1: ACAT Cauchy Combination Formula

**What:** Transform p-values to Cauchy variables, sum with weights, recover p-value via Cauchy CDF.
**When to use:** Any time multiple independent or correlated p-values need to be combined.

```python
# Source: Liu & Xie (2020, AJHG), doi:10.1016/j.ajhg.2019.01.002
# Verified against R ACAT package (github.com/yaowuliu/ACAT)
import numpy as np
from scipy.stats import cauchy

def acat_combine(p_values: np.ndarray, weights: np.ndarray | None = None) -> float:
    """
    Aggregate Cauchy Association Test (ACAT) p-value combination.

    Parameters
    ----------
    p_values : np.ndarray
        Array of p-values in [0, 1]. None values are dropped before computation.
    weights : np.ndarray or None
        Non-negative weights. None means equal weights. Normalized internally.

    Returns
    -------
    float
        Combined p-value in [0, 1]. Returns 1.0 if all inputs are None or 1.0.
    """
    p_arr = np.asarray([p for p in p_values if p is not None], dtype=float)
    if len(p_arr) == 0:
        return None  # No valid inputs

    if weights is None:
        w = np.ones(len(p_arr))
    else:
        w = np.asarray(weights, dtype=float)
    w = w / w.sum()  # normalize

    # Edge case: any p >= 1.0 → skip (treat as non-informative)
    keep = p_arr < 1.0
    if not keep.any():
        return 1.0
    p_arr = p_arr[keep]
    w = w[keep]
    w = w / w.sum()  # re-normalize after filtering

    # Edge case: p == 0 or p < 1e-16 → use 1/(p*pi) approximation
    # (tan((0.5 - p)*pi) overflows for p → 0)
    is_tiny = p_arr < 1e-16
    transformed = np.where(
        is_tiny,
        1.0 / (p_arr * np.pi),                 # approximation for p << 1
        np.tan((0.5 - p_arr) * np.pi)           # standard formula
    )

    T_stat = float(np.sum(w * transformed))
    p_combined = float(cauchy.sf(T_stat))
    return max(0.0, min(1.0, p_combined))
```

**Numerical verification:**
- `acat_combine([0.05, 0.05, 0.05])` = 0.05 (equal inputs → same output)
- `acat_combine([0.5, 0.5])` = 0.5 (null inputs → null output)
- `acat_combine([0.001, 0.01, 0.05])` ≈ 0.00268 (drives toward smallest)

### Pattern 2: ACAT-V (Variant-Level Cauchy Combination)

**What:** Per-gene combination of per-variant marginal score test p-values. Weights are `w_i = w_skat_i * MAF_i * (1 - MAF_i)`.
**When to use:** When per-variant marginal p-values are available (from null model residuals).

```python
# Source: Liu & Xie (2020, AJHG) — equation for w_{i,ACAT-V}
# w_i_acat_v = w_i_skat * MAF_i * (1 - MAF_i)
# where w_i_skat = dbeta(MAF_i; a1, a2) (Beta distribution density)

from scipy.stats import beta as beta_dist

def compute_acat_v_weights(mafs: np.ndarray, beta_a: float = 1.0, beta_b: float = 25.0) -> np.ndarray:
    """ACAT-V variant weights from MAF using SKAT Beta weight convention."""
    skat_weights = beta_dist.pdf(mafs, a=beta_a, b=beta_b)
    acat_v_weights = skat_weights * mafs * (1 - mafs)
    return acat_v_weights
```

**Note on marginal score test p-values:** Computing per-variant p-values requires null model residuals `r = Y - mu_hat`. For binary traits, the score statistic for variant j is `S_j = G_j' r` and `Var(S_j) = G_j' diag(mu*(1-mu)) G_j`. This requires the fitted null model — same null model as SKAT. This computation is NOT trivial for binary traits (saddlepoint or moment-match needed for accuracy). See OMNI-01 note below.

### Pattern 3: ACAT-O (Omnibus Test)

**What:** Combine burden_p + skat_p (and optionally acat_v_p) per gene with equal weights.
**When to use:** After all primary tests have run for a gene.

```python
# ACAT-O for this codebase: combine available test p-values with equal weights
# Requirements say: burden + SKAT + ACAT-V
# But ACAT-V requires per-variant p-values; if not available, combine what exists

def compute_acat_o(
    burden_p: float | None,
    skat_p: float | None,
    acat_v_p: float | None = None,
) -> float | None:
    """Omnibus ACAT-O combining burden + SKAT + ACAT-V p-values."""
    inputs = [p for p in [burden_p, skat_p, acat_v_p] if p is not None]
    if len(inputs) < 2:
        return None  # Need at least 2 tests to combine
    return acat_combine(np.array(inputs))
```

### Pattern 4: Lambda_GC Computation

**What:** Genomic inflation factor — ratio of observed median chi2 to expected median chi2(1).
**When to use:** Per test, after all genes tested.

```python
# Source: Devlin & Roeder (1999); standard GWAS diagnostic
# Formula: lambda_GC = median(chi2_observed) / qchisq(0.5, df=1)
# where chi2_observed = chi2.isf(p_values, df=1)
# and qchisq(0.5, df=1) = 0.4549364... (verified: scipy.stats.chi2.ppf(0.5, df=1))

from scipy.stats import chi2 as chi2_dist

def compute_lambda_gc(p_values: np.ndarray) -> float | None:
    """
    Compute genomic inflation factor lambda_GC from a vector of p-values.

    Parameters
    ----------
    p_values : np.ndarray
        Raw (uncorrected) p-values. None/NaN values are excluded.

    Returns
    -------
    float or None
        Lambda_GC. None if fewer than 2 valid p-values.
    """
    _EXPECTED_MEDIAN_CHI2_1 = 0.45493642311957174  # chi2.ppf(0.5, df=1)
    valid = np.asarray([p for p in p_values if p is not None and not np.isnan(p)], dtype=float)
    valid = np.clip(valid, 1e-300, 1.0 - 1e-15)  # guard against 0 and 1
    if len(valid) < 2:
        return None
    chi2_obs = chi2_dist.isf(valid, df=1)
    return float(np.median(chi2_obs) / _EXPECTED_MEDIAN_CHI2_1)
```

**Key constant:** `qchisq(0.5, 1) = chi2.ppf(0.5, df=1) = 0.45493642311957174`
(not 0.456 as sometimes cited — use the exact scipy value).

### Pattern 5: QQ Plot Data Generation

**What:** Observed vs expected -log10(p) for each test; output as TSV.
**When to use:** After lambda_GC, same p-value vector.

```python
def compute_qq_data(p_values: np.ndarray, test_name: str) -> pd.DataFrame:
    """
    Compute QQ plot data (observed vs expected -log10(p)).

    Uses Hazen quantile formula for expected: i/(n+1) for rank i in 1..n.
    Returns DataFrame with columns: test, expected_log10p, observed_log10p.
    """
    valid = np.sort([p for p in p_values if p is not None and not np.isnan(p)])
    n = len(valid)
    if n == 0:
        return pd.DataFrame(columns=["test", "expected_log10p", "observed_log10p"])

    expected = np.arange(1, n + 1) / (n + 1)  # Hazen quantile formula
    # Both sorted ascending → smallest p at index 0 → largest -log10(p) at index 0
    obs_log = -np.log10(np.clip(valid, 1e-300, 1.0))
    exp_log = -np.log10(expected)

    return pd.DataFrame({
        "test": test_name,
        "expected_log10p": exp_log[::-1],   # ascending expected (small first)
        "observed_log10p": obs_log[::-1],   # ascending observed
    })
```

### Pattern 6: Sample Size Warnings

**What:** Warn when analysis power is likely inadequate.
**When to use:** Before tests run (at stage initialization), and per-gene for case_carriers.

```python
def emit_sample_size_warnings(
    n_cases: int,
    n_controls: int,
    gene: str | None = None,
    case_carriers: int | None = None,
) -> list[str]:
    """
    Emit standard sample size warnings per DIAG-06.

    Thresholds:
    - n_cases < 200: low power warning
    - case:control ratio > 1:20: imbalanced cohort warning
    - case_carriers < 10 (per gene): low carrier count flag
    """
    warnings_out = []
    if n_cases < 200:
        msg = f"Low case count (n_cases={n_cases} < 200): power may be insufficient"
        logger.warning(msg)
        warnings_out.append("LOW_CASE_COUNT")
    if n_controls > 0 and n_cases > 0 and (n_controls / n_cases) > 20:
        msg = f"Imbalanced cohort (1:{n_controls // n_cases} case:control ratio > 1:20)"
        logger.warning(msg)
        warnings_out.append("IMBALANCED_COHORT")
    if gene is not None and case_carriers is not None and case_carriers < 10:
        msg = f"Gene {gene}: case_carriers={case_carriers} < 10 (low power for this gene)"
        logger.warning(msg)
        warnings_out.append("LOW_CARRIER_COUNT")
    return warnings_out
```

### Pattern 7: FDR Strategy Change (ARCH-03)

The current engine applies FDR per-test (lines 182-198 of engine.py). Phase 22 changes this:
- Primary tests (fisher, burden, skat): corrected_p_value set to None (no per-test FDR).
- ACAT-O only: single FDR pass across all genes.

```python
# In engine.run_all() — after ACAT-O computed, before building DataFrame:
acat_o_genes = [g for g in all_genes if results_by_test["acat_o"][g].p_value is not None]
if acat_o_genes:
    raw_acat_o = [results_by_test["acat_o"][g].p_value for g in acat_o_genes]
    corrected_acat_o = apply_correction(raw_acat_o, self._config.correction_method)
    for gene, corr_p in zip(acat_o_genes, corrected_acat_o):
        results_by_test["acat_o"][gene].corrected_p_value = float(corr_p)
```

**IMPORTANT:** The current engine's per-test FDR loop (lines 182-198) must be removed or scoped to skip ACAT-O (which computes its own FDR separately). Primary tests will have `corrected_p_value = None` in the output going forward (per ARCH-03). This is a breaking change to output schema — plan accordingly.

### Anti-Patterns to Avoid

- **Don't register ACAT-O as an AssociationTest ABC subclass that runs in the gene loop.** It has no contingency_data — it reads from other tests' results. The ABC contract requires contingency_data; violating it leads to a hack.
- **Don't apply FDR separately per test.** ARCH-03 is explicit: single FDR on ACAT-O. The current per-test FDR in engine.py must be removed for primary tests.
- **Don't use `math.tan` for numpy arrays.** Use `np.tan` — `math.tan` cannot vectorize.
- **Don't clip p=0 to a small epsilon.** Use the `1/(p*pi)` branch for p < 1e-16 as the ACAT paper specifies.
- **Don't use 0.456 as the lambda_GC denominator.** The exact scipy value is 0.45493642... — hardcode from `chi2.ppf(0.5, df=1)` once at module level.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Cauchy CDF/SF | Custom CDF | `scipy.stats.cauchy.sf()` | Already verified equivalent to formula; handles numerical edge cases |
| Chi2 inverse CDF | Newton's method or lookup | `scipy.stats.chi2.isf(p, df=1)` | Exact, fast, handles boundary cases |
| FDR correction | Benjamini-Hochberg from scratch | `correction.apply_correction()` (already exists) | Already tested, guarded for empty input |
| Directory creation | Manual path checks | `pathlib.Path.mkdir(parents=True, exist_ok=True)` | Idiomatic Python |

**Key insight:** The entire ACAT, lambda_GC, and QQ computation is pure numpy/scipy — no new dependencies, no new algorithms, just the formulas applied to existing p-value arrays.

## Common Pitfalls

### Pitfall 1: ACAT-V Per-Variant P-Value Complexity

**What goes wrong:** The requirements state OMNI-01 = "ACAT-V per-variant Cauchy combination of marginal score test p-values per gene." Computing per-variant marginal score test p-values for binary traits requires saddlepoint approximation or moment-matching (normal approximation is anti-conservative for rare variants with MAC < 10). This is the same complexity as one step of SKAT's Davies computation.

**Why it happens:** Implementers treat `z = S/sqrt(Var(S))` → `p = 2*norm.sf(|z|)` as sufficient. It is not for rare binary variants.

**How to avoid:** Check whether the R backend exposes per-variant p-values from `SKAT_Null_Model` already. The R SKAT package's `Get.marginal.pval()` function does this correctly. If using R backend: call `Get.marginal.pval()` via rpy2. If pure Python: use normal approximation but flag as approximate for MAC < 10.

**Warning signs:** Lambda_GC on ACAT-V p-values >> 1.0 on a permuted null phenotype (should be ~1.0).

### Pitfall 2: FDR Behavior Change is Breaking

**What goes wrong:** Removing per-test FDR and applying it only to ACAT-O changes all existing output columns. Tests that previously checked `fisher_corrected_p_value` will find `None`.

**Why it happens:** ARCH-03 is a policy decision that touches engine.py column generation.

**How to avoid:** Update tests that check `corrected_p_value` on primary tests. The output schema change is intentional — document it in the plan.

**Warning signs:** Test failures on `test_association_engine.py` asserting `corrected_p_value is not None` for Fisher/burden.

### Pitfall 3: ACAT-O Requires at Least 2 Non-None P-Values

**What goes wrong:** Gene has only SKAT result (burden skipped due to convergence failure). `acat_combine([None, 0.03])` with filtering returns `acat_combine([0.03])` — a single p-value. ACAT of 1 value = that value unchanged, which is mathematically correct but misleading.

**Why it happens:** Not all tests produce valid p-values for all genes.

**How to avoid:** Require at least 2 non-None inputs to compute ACAT-O; return `None` otherwise. Document this skip condition.

### Pitfall 4: Lambda_GC on Small Gene Panels is Unreliable

**What goes wrong:** With < 100 genes, lambda_GC is extremely noisy. A value of 0.95–1.05 does NOT confirm calibration.

**Why it happens:** The test requirement DIAG-02 says "lambda_GC within [0.95, 1.05] on permuted null phenotype." This is only reliable for genome-wide studies (10,000+ tests).

**How to avoid:** Warn in logs when n_tests < 100 that lambda_GC is unreliable. Test lambda_GC with 1000+ permuted genes, not 10.

### Pitfall 5: QQ Data Sort Order

**What goes wrong:** QQ data written in wrong order makes the QQ plot look like a scatter. The convention is x-axis = expected (ascending), y-axis = observed (ascending), so the diagonal runs bottom-left to top-right.

**Why it happens:** Sorting by -log10(p) ascending means the largest values (most significant) come last — correct for QQ convention.

**How to avoid:** Sort both expected and observed ascending (by expected), so the output TSV rows go from small-expected (non-significant) to large-expected (most significant). Verify first few rows have small values.

### Pitfall 6: Missing ACAT-O When SKAT Not Run

**What goes wrong:** User runs only `--association-tests fisher,logistic_burden` (no SKAT). ACAT-O cannot be computed because SKAT p-value is missing.

**Why it happens:** ACAT-O is defined as combining burden + SKAT. If SKAT isn't run, ACAT-O is undefined.

**How to avoid:** Guard ACAT-O computation: skip if fewer than 2 of the required p-values (burden_p, skat_p) are non-None. Log an INFO message explaining why ACAT-O was skipped.

## Code Examples

### Full ACAT Module Pattern
```python
# Source: Pattern from existing tests/acat.py plan + Liu & Xie (2020)
# File: variantcentrifuge/association/tests/acat.py

import numpy as np
from scipy.stats import cauchy

_TINY_P_THRESHOLD = 1e-16

def cauchy_combination(p_values: np.ndarray, weights: np.ndarray | None = None) -> float | None:
    """
    Cauchy combination test. Core of both ACAT-V and ACAT-O.

    Returns None if fewer than 2 valid (non-None, < 1.0) p-values.
    Returns 1.0 if all p-values are >= 1.0.
    """
    valid_mask = np.array([(p is not None and not np.isnan(p) and 0 <= p < 1.0)
                           for p in p_values])
    p_arr = np.asarray(p_values, dtype=float)[valid_mask]

    if len(p_arr) < 2:
        return None if len(p_arr) < 1 else float(p_arr[0])

    if weights is not None:
        w = np.asarray(weights, dtype=float)[valid_mask]
    else:
        w = np.ones(len(p_arr))
    w = w / w.sum()

    transformed = np.where(
        p_arr < _TINY_P_THRESHOLD,
        1.0 / (p_arr * np.pi),
        np.tan((0.5 - p_arr) * np.pi)
    )
    T = float(np.dot(w, transformed))
    return float(np.clip(cauchy.sf(T), 0.0, 1.0))
```

### Lambda_GC and QQ in diagnostics.py
```python
# Source: Devlin & Roeder (1999); standard genomics diagnostic
# File: variantcentrifuge/association/diagnostics.py

from scipy.stats import chi2 as chi2_dist
import numpy as np
import pandas as pd

_EXPECTED_CHI2_MEDIAN = chi2_dist.ppf(0.5, df=1)  # 0.45493642311957174

def compute_lambda_gc(p_values: list[float | None]) -> float | None:
    valid = np.array([p for p in p_values if p is not None and 0 < p <= 1], dtype=float)
    if len(valid) < 2:
        return None
    return float(np.median(chi2_dist.isf(valid, df=1)) / _EXPECTED_CHI2_MEDIAN)

def compute_qq_data(p_values: list[float | None], test_name: str) -> pd.DataFrame:
    valid = np.sort([p for p in p_values if p is not None and 0 < p <= 1])
    n = len(valid)
    if n == 0:
        return pd.DataFrame(columns=["test", "expected_neg_log10_p", "observed_neg_log10_p"])
    expected = np.arange(1, n + 1) / (n + 1)
    return pd.DataFrame({
        "test": test_name,
        "expected_neg_log10_p": -np.log10(expected),
        "observed_neg_log10_p": -np.log10(np.clip(valid, 1e-300, 1.0)),
    }).sort_values("expected_neg_log10_p").reset_index(drop=True)
```

### Diagnostics Stage / Post-Processing Pattern
```python
# After engine.run_all() returns results_df, in AssociationAnalysisStage._process():
diagnostics_output = context.config.get("diagnostics_output")
if diagnostics_output:
    from variantcentrifuge.association.diagnostics import (
        compute_lambda_gc, compute_qq_data, emit_sample_size_warnings
    )
    diag_dir = Path(diagnostics_output)
    diag_dir.mkdir(parents=True, exist_ok=True)

    # Lambda_GC per test
    lambda_records = []
    p_col_tests = [c.replace("_p_value", "") for c in results_df.columns
                   if c.endswith("_p_value")]
    for test in p_col_tests:
        p_col = f"{test}_p_value"
        lam = compute_lambda_gc(results_df[p_col].tolist())
        lambda_records.append({"test": test, "lambda_gc": lam})
    pd.DataFrame(lambda_records).to_csv(diag_dir / "lambda_gc.tsv", sep="\t", index=False)

    # QQ data per test
    for test in p_col_tests:
        p_col = f"{test}_p_value"
        qq_df = compute_qq_data(results_df[p_col].tolist(), test)
        qq_df.to_csv(diag_dir / f"qq_data_{test}.tsv", sep="\t", index=False)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Per-test FDR in engine.py | Single FDR on ACAT-O only | Phase 22 (ARCH-03) | Output schema change: primary tests lose `corrected_p_value` |
| No omnibus test | ACAT-O omnibus combining burden + SKAT | Phase 22 | Adds ACAT-O p-value per gene as the primary result |
| No diagnostics | lambda_GC + QQ TSV per test | Phase 22 | Enables QC of association results |

**Note on DIAG-01 "per-gene TSV with all standard output columns":** This is already produced by the engine's `run_all()` → `results_df.to_csv()`. Phase 22 does not add a new output file here — it ensures the existing output file includes `acat_o_p`, `acat_o_corrected_p`, and sample-size warning flags as new columns.

## Open Questions

1. **ACAT-V per-variant p-value source (OMNI-01)**
   - What we know: ACAT-V requires marginal score test p-values per variant per gene. The R ACAT package's `Get.marginal.pval()` computes these from the SKAT null model.
   - What's unclear: Does the R SKAT backend (rpy2) expose `Get.marginal.pval()` results, or does the planner need to add a separate rpy2 call? Is ACAT-V computable with only the pure Python burden test (no SKAT null model)?
   - Recommendation: For Phase 22, compute ACAT-V using normal approximation (z = score/sqrt(var_score)) for a first implementation, document that it is approximate for MAC < 10. Exact saddlepoint implementation deferred to Phase 23 or later.

2. **FDR on ACAT-O only vs maintaining per-test corrected_p_value**
   - What we know: ARCH-03 says single FDR on ACAT-O. Current engine applies FDR per-test.
   - What's unclear: Do downstream stages (Excel report, HTML report) consume `corrected_p_value` columns from primary tests? If so, removing them is a breaking UI change.
   - Recommendation: Keep per-test `corrected_p_value` columns populated (for display) but make ACAT-O the canonical significance column. Clarify with orchestrator before plan.

3. **ACAT-O when only 1 test produces a valid p-value**
   - What we know: With n=1 non-None p-value, ACAT of 1 value = that p-value unchanged.
   - What's unclear: Is this correct behavior or should it return None?
   - Recommendation: Return None if fewer than 2 non-None inputs (require at least burden + SKAT).

## Sources

### Primary (HIGH confidence)
- Liu & Xie (2020, AJHG) — "ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies" — PMC6407498 — ACAT formula, ACAT-V weights, ACAT-O definition
- github.com/yaowuliu/ACAT — R reference implementation — edge case handling (p=0: `1/p/pi`, p=1: skip, mixed 0+1: error)
- scipy.stats.cauchy documentation — `cauchy.sf()` equivalent to Cauchy CDF formula
- scipy.stats.chi2 documentation — `chi2.ppf(0.5, df=1)` = 0.45493642... for lambda_GC

### Secondary (MEDIUM confidence)
- lawlessgenomics.com/topic/acat — Summary of ACAT algorithm with formula verification
- felixfan.github.io/Genomic-Inflation-Factor — Lambda_GC formula (verified against scipy)
- Devlin & Roeder (1999) — Original genomic control paper (lambda_GC formula)

### Tertiary (LOW confidence)
- gist.github.com/ryananeff/c66cdf086979b13e855f2c3d0f3e54e1 — Python ACAT port (not official; formula correct but not authoritative)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all libraries already in deps, formulas verified numerically
- ACAT formula + edge cases: HIGH — verified against R reference implementation
- Lambda_GC formula: HIGH — verified with scipy (exact constant: 0.45493642...)
- QQ data generation: HIGH — standard GWAS practice, verified in code
- ACAT-V per-variant p-values: MEDIUM — computation path with R backend not verified; normal approximation is approximate
- Architecture (Option B post-loop ACAT-O): MEDIUM — fits existing engine structure; alternative ABC approach also viable
- FDR strategy change: HIGH — ARCH-03 is explicit in STATE.md

**Research date:** 2026-02-21
**Valid until:** 2026-03-21 (stable algorithms, no expiration concern)
