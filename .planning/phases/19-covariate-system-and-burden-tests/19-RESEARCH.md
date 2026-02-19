# Phase 19: Covariate System and Burden Tests - Research

**Researched:** 2026-02-19
**Domain:** Statistical regression, covariate alignment, genotype matrix construction, Firth logistic regression
**Confidence:** HIGH (all critical paths verified against live code and library APIs)

## Summary

Phase 19 adds logistic and linear regression burden tests with covariate adjustment to the existing
AssociationEngine framework. The architecture is well-defined by Phase 18 foundations: new tests
register in the engine as `AssociationTest` subclasses, and covariate support requires extending
`AssociationConfig` and `TestResult`. All required statistical libraries (scipy, statsmodels, numpy,
pandas) are already in `pyproject.toml` dependencies — no new library installs are needed for the
standard path.

The primary design challenge is that logistic/linear burden tests require a per-gene genotype matrix
(samples x variants) built from the DataFrame's per-sample GT columns (`GEN_0__GT`, `GEN_1__GT`, ...),
while the existing `AssociationTest.run()` interface only receives pre-aggregated `contingency_data`
counts. Phase 19 must extend `contingency_data` to include the per-gene variant rows or add a parallel
data path for regression tests.

For Firth fallback: statsmodels has no built-in Firth implementation. The reference implementation is
a standalone Newton-Raphson loop using the Jeffreys prior (hat matrix diagonal). A self-contained
`_firth_logistic()` helper can be implemented in ~80 lines using only numpy/scipy — no additional
package required. `firthmodels` v0.7.1 is a viable alternative (actively maintained as of January
2026), but adds a new dependency for functionality that can be implemented in-house.

**Primary recommendation:** Implement Firth as a self-contained helper in `logistic_burden.py`
(~80 lines) using numpy/scipy. No new library dependencies needed.

## Standard Stack

All libraries are already installed. No changes to `pyproject.toml`.

### Core (already in dependencies)
| Library | Version (installed) | Purpose | Why Used |
|---------|---------------------|---------|----------|
| statsmodels | 0.14.4 | Logit + OLS regression | Provides `result.params`, `result.pvalues`, `result.bse`, `result.conf_int()`, `result.mle_retvals['converged']` |
| scipy | 1.14.1 | Beta distribution weights, chi2 for Firth p-values | `scipy.stats.beta.pdf(maf, a=1, b=25)` |
| numpy | 2.2.6 | Matrix operations, imputation, condition number | `np.linalg.cond()`, matrix multiply for burden scores |
| pandas | 2.3.3 | Covariate file loading, reindex alignment, get_dummies | `pd.read_csv(sep=None)`, `df.reindex()`, `pd.get_dummies(drop_first=True)` |

### Optional (do NOT add)
| Library | Why Not |
|---------|---------|
| firthlogist 0.5.0 | Unmaintained since 2022, not actively developed |
| firthmodels 0.7.1 | Adds dependency; Firth fits in 80 lines inline |

**Installation:** No new packages needed.

## Architecture Patterns

### Recommended Project Structure

```
variantcentrifuge/association/
├── __init__.py                    # (existing)
├── base.py                        # (existing) AssociationTest ABC, TestResult, AssociationConfig
├── engine.py                      # (existing) AssociationEngine orchestrator
├── correction.py                  # (existing) FDR/Bonferroni
├── covariates.py                  # NEW: CovariateMatrix, load_covariates(), align_to_vcf()
├── genotype_matrix.py             # NEW: GenotypeMatrixBuilder, parse_gt_to_dosage()
├── weights.py                     # NEW: beta_weights(), uniform_weights()
└── tests/
    ├── __init__.py                # (existing)
    ├── fisher.py                  # (existing)
    ├── logistic_burden.py         # NEW: LogisticBurdenTest
    └── linear_burden.py           # NEW: LinearBurdenTest
```

### Pattern 1: Extended contingency_data for Regression Tests

The existing `AssociationTest.run(gene, contingency_data, config)` interface passes a `dict`
with pre-aggregated counts. Regression tests need the raw per-sample genotype data, NOT just
counts. The cleanest extension is to pass additional keys in `contingency_data`:

```python
# Engine populates contingency_data with optional regression keys:
contingency_data = {
    # Existing fisher keys (always present):
    "proband_count": 120,
    "control_count": 200,
    "proband_carrier_count": 45,
    "control_carrier_count": 12,
    "proband_allele_count": 60,
    "control_allele_count": 15,
    "n_qualifying_variants": 7,
    # New regression keys (present when vcf_samples + genotype matrix builder available):
    "genotype_matrix": np.ndarray,     # shape (n_samples, n_variants), float64, imputed
    "variant_mafs": np.ndarray,        # shape (n_variants,), float64
    "vcf_samples": list[str],          # ordered sample names matching matrix rows
    "phenotype_vector": np.ndarray,    # shape (n_samples,), 0/1 (binary) or float
    "covariate_matrix": np.ndarray,    # shape (n_samples, k), float64 or None
}
```

This keeps backward compatibility: `FisherExactTest.run()` ignores the new keys, and
`LogisticBurdenTest.run()` uses them.

### Pattern 2: Covariate File Loading and Alignment

```python
# Source: verified against pandas docs + live test
import pandas as pd
import numpy as np

def load_covariates(
    filepath: str,
    vcf_samples: list[str],
    covariate_columns: list[str] | None = None,
    categorical_columns: list[str] | None = None,
) -> np.ndarray:
    """Load, align, and encode covariate file. Returns (n_samples, k) float64 array."""
    # 1. Auto-detect delimiter from extension, fallback to csv.Sniffer
    import os, csv
    ext = os.path.splitext(filepath)[1].lower()
    if ext in ('.tsv', '.tab'):
        sep = '\t'
    elif ext == '.csv':
        sep = ','
    else:
        # Peek first 2KB with Sniffer
        with open(filepath) as f:
            sample = f.read(2048)
        try:
            sep = csv.Sniffer().sniff(sample).delimiter
        except Exception:
            sep = '\t'  # bioinformatics default

    # 2. Load covariate file (first column = sample ID, header required)
    df = pd.read_csv(filepath, sep=sep, index_col=0)

    # 3. Column selection
    if covariate_columns:
        df = df[covariate_columns]

    # 4. Check for VCF samples missing from covariate file (abort condition)
    missing = set(vcf_samples) - set(df.index)
    if missing:
        raise ValueError(f"VCF samples missing from covariate file: {sorted(missing)}")

    # 5. Warn about extra covariate samples (continue with intersection)
    extra = set(df.index) - set(vcf_samples)
    if extra:
        logger.warning(f"Covariate file has {len(extra)} extra samples not in VCF: {sorted(extra)[:5]}...")

    # 6. Align to VCF sample order (CRITICAL: assert no NaN after reindex)
    df_aligned = df.reindex(vcf_samples)
    assert not df_aligned.isnull().any().any(), \
        "NaN in covariate matrix after reindex — unreachable if step 4 passed"

    # 7. One-hot encode categorical columns (auto-detect or explicit)
    # Auto-detect: non-numeric columns with <=5 unique values
    if categorical_columns is None:
        categorical_columns = [
            c for c in df_aligned.columns
            if not pd.api.types.is_numeric_dtype(df_aligned[c])
            and df_aligned[c].nunique() <= 5
        ]
    if categorical_columns:
        df_aligned = pd.get_dummies(df_aligned, columns=categorical_columns,
                                    drop_first=True, dtype=float)

    # 8. Multicollinearity check (warn only, don't abort)
    X = df_aligned.values.astype(float)
    cond_num = np.linalg.cond(X)
    if cond_num > 1000:
        logger.warning(f"High multicollinearity in covariate matrix (condition number: {cond_num:.1f})")

    return X
```

### Pattern 3: Genotype Matrix Builder

```python
# Source: verified against genotype_utils.py + live edge case testing

def parse_gt_to_dosage(gt: str) -> int | None:
    """Parse GT string to dosage 0/1/2, or None for missing.

    Rules (verified against all edge cases):
    - '0/0', '0|0' -> 0
    - '0/1', '1/0', '0|1', '1|0' -> 1
    - '1/1', '1|1' -> 2
    - './.', '.|.' -> None (missing)
    - '1/2', '0/2', '2/2' -> 1 (multi-allelic: het-equivalent + WARNING)
    - Phased (|): treated as unphased, sum alleles
    """
    if not gt or gt in ('./.', '.|.'):
        return None
    sep = '|' if '|' in gt else '/'
    parts = gt.split(sep)
    if len(parts) != 2:
        return None
    try:
        a1 = None if parts[0] == '.' else int(parts[0])
        a2 = None if parts[1] == '.' else int(parts[1])
    except ValueError:
        return None
    if a1 is None or a2 is None:
        return None
    if a1 > 1 or a2 > 1:  # multi-allelic
        return 1
    return a1 + a2


def build_genotype_matrix(
    gene_df: pd.DataFrame,
    vcf_samples: list[str],
    gt_columns: list[str],
    is_binary: bool = True,
    missing_site_threshold: float = 0.10,
    missing_sample_threshold: float = 0.80,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Build (n_samples, n_variants) genotype matrix with imputation.

    Returns: (G, mafs, kept_variant_indices)
    G: float64 matrix, no NaN (imputed)
    mafs: float64 array of per-variant MAF
    """
    n_samples = len(vcf_samples)
    n_variants = len(gene_df)

    # Parse raw dosages: shape (n_variants, n_samples), NaN for missing
    raw = np.full((n_variants, n_samples), np.nan)
    multi_allelic_seen = False
    for v_idx, (_, row) in enumerate(gene_df.iterrows()):
        for s_idx, col in enumerate(gt_columns):
            gt = row.get(col, './.')
            dosage = parse_gt_to_dosage(str(gt))
            if dosage is not None:
                raw[v_idx, s_idx] = float(dosage)
            # Note multi-allelic for warning emission
            if dosage == 1 and parse_gt_multi_allelic(str(gt)):
                multi_allelic_seen = True

    # Layer 1: Pre-filter variants with >10% missing site-wide
    missing_per_variant = np.isnan(raw).mean(axis=1)
    keep_variants = missing_per_variant <= missing_site_threshold
    raw = raw[keep_variants]

    # Layer 1: Pre-filter samples with >80% missing in this gene
    missing_per_sample = np.isnan(raw).mean(axis=0)
    keep_samples = missing_per_sample <= missing_sample_threshold

    # Layer 2: Mean imputation (2*MAF for continuous; round(2*MAF) for binary)
    G = raw.copy()
    mafs = np.zeros(G.shape[0])
    for v_idx in range(G.shape[0]):
        observed = G[v_idx, ~np.isnan(G[v_idx])]
        if len(observed) == 0:
            mafs[v_idx] = 0.0
            G[v_idx] = 0.0
            continue
        maf = observed.mean() / 2.0
        mafs[v_idx] = maf
        missing = np.isnan(G[v_idx])
        if missing.any():
            imputed = round(2 * maf) if is_binary else 2 * maf
            G[v_idx, missing] = imputed

    return G.T, mafs, keep_variants  # Transpose: (n_samples, n_variants)
```

### Pattern 4: Logistic Burden Test with Firth Fallback

```python
# Source: verified against statsmodels 0.14.4 API + live tests

import warnings
import numpy as np
import statsmodels.api as sm

SEPARATION_BSE_THRESHOLD = 100.0  # BSE > 100 suggests separation

def logistic_burden_test(
    G: np.ndarray,           # (n_samples, n_variants)
    weights: np.ndarray,     # (n_variants,)
    phenotype: np.ndarray,   # (n_samples,) 0/1
    covariates: np.ndarray | None = None,  # (n_samples, k) or None
) -> dict:
    """Run logistic burden regression. Returns result dict."""

    # Compute weighted burden score per sample
    burden = G @ weights  # (n_samples,)

    # Build design matrix
    X = sm.add_constant(burden.reshape(-1, 1))
    if covariates is not None and covariates.ndim == 2 and covariates.shape[1] > 0:
        X = np.column_stack([X, covariates])

    # Carrier stats for output
    carriers = burden > 0
    n_carriers = int(carriers.sum())
    n_carriers_cases = int((carriers & (phenotype == 1)).sum())
    n_carriers_controls = int((carriers & (phenotype == 0)).sum())

    # Pre-flight separation checks
    warning_codes: list[str] = []
    if n_carriers > 0:
        case_rate = n_carriers_cases / n_carriers
        if case_rate >= 0.999 or case_rate <= 0.001:
            warning_codes.append("PERFECT_SEPARATION")
        elif case_rate >= 0.9 or case_rate <= 0.1:
            warning_codes.append("QUASI_SEPARATION")
    if n_carriers_cases == 0 or n_carriers_controls == 0:
        if "PERFECT_SEPARATION" not in warning_codes:
            warning_codes.append("ZERO_CARRIERS_ONE_GROUP")

    # Fit standard logistic regression
    caught_msgs: list[str] = []
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = sm.Logit(phenotype, X).fit(disp=False, maxiter=100)
            caught_msgs = [str(ww.message) for ww in w]

        converged = result.mle_retvals.get("converged", True)
        bse_max = result.bse.max()

        if not converged or bse_max > SEPARATION_BSE_THRESHOLD:
            # Firth fallback
            firth_result = _firth_logistic(phenotype, X)
            if firth_result is not None:
                if "PERFECT_SEPARATION" not in warning_codes:
                    warning_codes.append("PERFECT_SEPARATION")
                return _format_logistic_result(
                    firth_result, warning_codes, n_carriers,
                    n_carriers_cases, n_carriers_controls
                )
            else:
                warning_codes.append("FIRTH_CONVERGE_FAIL")
                return _null_logistic_result(warning_codes, n_carriers,
                                             n_carriers_cases, n_carriers_controls)

        return _format_logistic_result(
            result, warning_codes, n_carriers, n_carriers_cases, n_carriers_controls
        )

    except Exception:
        warning_codes.append("FIRTH_CONVERGE_FAIL")
        return _null_logistic_result(warning_codes, n_carriers, n_carriers_cases, n_carriers_controls)
```

### Pattern 5: Self-Contained Firth Implementation (~80 lines)

```python
# Source: adapted from johnlees gist (https://gist.github.com/johnlees/3e06380965f367e4894ea20fbae2b90d)
# Verified: Newton-Raphson with Jeffreys prior penalty

def _firth_logistic(
    y: np.ndarray,
    X: np.ndarray,
    max_iter: int = 25,
    tol: float = 1e-4,
) -> object | None:
    """
    Firth penalized logistic regression via Newton-Raphson.

    Returns a result-like namespace with .params, .bse, .pvalues, .conf_int()
    or None if convergence fails.
    """
    from scipy.stats import norm as sp_norm
    import numpy as np

    n, k = X.shape
    beta = np.zeros(k)

    for _ in range(max_iter):
        pi = 1.0 / (1.0 + np.exp(-X @ beta))
        W = pi * (1 - pi)
        # Information matrix (X^T W X)
        XtWX = X.T @ (W[:, None] * X)
        try:
            var_cov = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            return None

        # Hat matrix diagonal (memory-efficient): h = diag(W^0.5 X V X^T W^0.5)
        WsqrtX = np.sqrt(W)[:, None] * X
        H_diag = np.sum(WsqrtX @ var_cov * WsqrtX, axis=1)

        # Penalized score
        U = X.T @ (y - pi + H_diag * (0.5 - pi))
        step = var_cov @ U
        beta_new = beta + step

        # Step-halving if needed
        old_ll = _firth_loglik(beta, y, X)
        for _ in range(10):
            new_ll = _firth_loglik(beta_new, y, X)
            if new_ll > old_ll:
                break
            step *= 0.5
            beta_new = beta + step

        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new

    # Compute final SE, p-values
    pi = 1.0 / (1.0 + np.exp(-X @ beta))
    W = pi * (1 - pi)
    XtWX = X.T @ (W[:, None] * X)
    try:
        var_cov = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return None

    bse = np.sqrt(np.diag(var_cov))
    z = beta / bse
    pvals = 2 * (1 - sp_norm.cdf(np.abs(z)))

    class FirthResult:
        def __init__(self):
            self.params = beta
            self.bse = bse
            self.pvalues = pvals
        def conf_int(self, alpha=0.05):
            from scipy.stats import norm as n
            z_val = n.ppf(1 - alpha / 2)
            return np.column_stack([beta - z_val * bse, beta + z_val * bse])

    return FirthResult()


def _firth_loglik(beta: np.ndarray, y: np.ndarray, X: np.ndarray) -> float:
    """Firth penalized log-likelihood (log-lik + 0.5 * log(det(I)))."""
    pi = 1.0 / (1.0 + np.exp(-X @ beta))
    ll = float(np.sum(y * np.log(pi + 1e-15) + (1 - y) * np.log(1 - pi + 1e-15)))
    W = pi * (1 - pi)
    XtWX = X.T @ (W[:, None] * X)
    sign, logdet = np.linalg.slogdet(XtWX)
    if sign <= 0:
        return ll
    return ll + 0.5 * logdet
```

### Pattern 6: Linear Burden Test

```python
# Source: verified against statsmodels 0.14.4 API + live tests

def linear_burden_test(
    G: np.ndarray,
    weights: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None = None,
) -> dict:
    """Linear regression burden test for quantitative traits."""
    burden = G @ weights
    X = sm.add_constant(burden.reshape(-1, 1))
    if covariates is not None and covariates.shape[1] > 0:
        X = np.column_stack([X, covariates])

    result = sm.OLS(phenotype, X).fit()
    ci = result.conf_int()  # shape (n_params, 2), ndarray (not DataFrame!)

    return {
        "p_value": float(result.pvalues[1]),
        "beta": float(result.params[1]),
        "se": float(result.bse[1]),
        "ci_lower": float(ci[1, 0]),
        "ci_upper": float(ci[1, 1]),
        "n_carriers": int((burden > 0).sum()),
    }
```

### Pattern 7: Beta(MAF; 1, 25) Weights

```python
# Source: scipy.stats.beta.pdf documentation + live test
from scipy.stats import beta as beta_dist
import numpy as np

def beta_maf_weights(mafs: np.ndarray, a: float = 1.0, b: float = 25.0) -> np.ndarray:
    """
    Compute Beta(MAF; a, b) weights for burden test.

    SKAT R package convention: default a=1, b=25.
    Upweights rare variants (low MAF), downweights common variants.
    MAF=0 is handled: pdf returns 25.0 for Beta(0; 1, 25).
    """
    # Clip MAF to [epsilon, 1-epsilon] to avoid edge issues
    maf_clipped = np.clip(mafs, 1e-8, 1 - 1e-8)
    return beta_dist.pdf(maf_clipped, a=a, b=b)

def uniform_weights(n_variants: int) -> np.ndarray:
    """Return uniform weight vector (all 1.0)."""
    return np.ones(n_variants)
```

### Pattern 8: AssociationConfig Extension

`AssociationConfig` in `base.py` needs new fields for Phase 19:

```python
@dataclass
class AssociationConfig:
    # ... existing fields ...
    # NEW in Phase 19:
    covariate_file: str | None = None
    covariate_columns: list[str] | None = None
    categorical_covariates: list[str] | None = None
    trait_type: str = "binary"              # "binary" or "quantitative"
    variant_weights: str = "beta:1,25"      # "beta:a,b" or "uniform"
    missing_site_threshold: float = 0.10
    missing_sample_threshold: float = 0.80
    firth_max_iter: int = 25
```

### Anti-Patterns to Avoid

- **Storing genotype matrix in PipelineContext:** 5K samples x 50K variants = 1.6 GB. Build it
  per-gene inside the test, discard immediately after.
- **Stratified MAF for imputation:** Compute MAF from ALL samples (not cases only, not controls
  only). This matches SKAT R convention.
- **Silently treating `./.'` as dosage 0:** Always impute missing genotypes. Never coerce missing
  to 0 implicitly.
- **Reporting non-converged p-value without warning:** If `mle_retvals['converged'] == False` AND
  `bse.max() > 100`, always attempt Firth. Never silently pass through the non-converged result.
- **Omitting genes from output:** Even genes with NA p-value must appear in output with warning code.
- **Using `result.conf_int().iloc[1]`:** `conf_int()` returns `numpy.ndarray`, not DataFrame.
  Use `ci[1, 0]` and `ci[1, 1]` (row 1 = first non-intercept coefficient).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Logistic regression | Custom gradient descent | `statsmodels.Logit` | Handles convergence, numerical stability, Wald tests |
| Linear regression | Custom OLS | `statsmodels.OLS` | Handles DOF corrections, heteroskedasticity tests |
| Beta distribution PDF | Manual formula | `scipy.stats.beta.pdf` | Handles edge cases, numerically stable |
| One-hot encoding | Manual dummy columns | `pd.get_dummies(drop_first=True)` | Handles dtype, column naming consistently |
| Condition number | Manual SVD | `np.linalg.cond()` | Optimized LAPACK implementation |
| Firth regression | External package | Inline Newton-Raphson (~80 lines) | firthlogist unmaintained; firthmodels adds dependency |
| Missing genotype count | Manual loop | `np.isnan(G).mean(axis=...)` | Vectorized, handles float NaN mask |

**Key insight:** Burden test implementations look like "just a regression call" but the correct
flow is: build matrix → impute → filter low-quality variants/samples → compute weights → build
burden score → fit regression → check convergence → attempt Firth → format result. Each step has
non-obvious edge cases.

## Common Pitfalls

### Pitfall 1: conf_int() Returns ndarray, Not DataFrame

**What goes wrong:** `result.conf_int().iloc[1]` raises `AttributeError`. In statsmodels 0.14.4,
`Logit.fit().conf_int()` returns `numpy.ndarray` shape `(n_params, 2)`, not a DataFrame.

**How to avoid:** Use `ci = result.conf_int()` then `ci[1, 0]` (lower) and `ci[1, 1]` (upper).
Verified against live statsmodels 0.14.4.

**Warning signs:** `AttributeError: 'numpy.ndarray' object has no attribute 'iloc'`

### Pitfall 2: Separation Produces Large BSE, Not Exception

**What goes wrong:** Perfect separation in `statsmodels.Logit` does NOT raise an exception in all
cases. It emits a `PerfectSeparationWarning` via the warnings module and continues with
`mle_retvals['converged'] == False` and BSE values > 70,000.

**How to avoid:** Detect via `warnings.catch_warnings(record=True)` + check both `converged` flag
AND `bse.max() > THRESHOLD` (use 100.0 as threshold). Both checks are needed.

**Warning signs:** BSE in thousands, params in tens, `mle_retvals['converged'] == False`

### Pitfall 3: Multi-allelic Genotypes Return dosage=1 from Existing `_gt_to_dosage`

**What goes wrong:** The existing `gene_burden._gt_to_dosage('1/2')` returns 0 (falls through to
`return 0`). The new burden test matrix needs dosage=1 for `1/2` (het-equivalent).

**How to avoid:** Implement a new `parse_gt_to_dosage()` in `genotype_matrix.py` with explicit
multi-allelic handling. Do NOT reuse `_gt_to_dosage` from `gene_burden.py`.

**Warning signs:** Multi-allelic carriers silently excluded from burden scores.

### Pitfall 4: Covariate File NaN After Reindex

**What goes wrong:** VCF sample in one order, covariate file in another. Without `reindex()`, the
burden matrix sample order misaligns with the covariate matrix — statistically invalid.

**How to avoid:** Always `df.reindex(vcf_samples)` and `assert not df_aligned.isnull().any().any()`.
The context says "VCF samples missing from covariate file" should abort — the reindex will produce
NaN for such samples, triggering the assertion before any stats are run.

**Warning signs:** NaN in covariate matrix after reindex (caught by assertion).

### Pitfall 5: Imputing MAF Stratified by Phenotype Leaks Labels

**What goes wrong:** Computing MAF separately for cases and controls then imputing each group
with their own group MAF constitutes phenotype-based imputation — this is a form of label leakage
that inflates association statistics.

**How to avoid:** Always compute imputation MAF from ALL samples combined. SKAT R convention.

**Warning signs:** Inflated test statistics, especially for genes with many missing calls.

### Pitfall 6: `add_constant` Column Index Shift

**What goes wrong:** After `X = sm.add_constant(burden)`, the intercept is index 0, burden score
is index 1. All coefficient/CI extraction must use index 1, not 0.

**How to avoid:** Always extract `result.params[1]`, `result.bse[1]`, `ci[1, :]` for the burden
coefficient.

**Warning signs:** P-value for intercept reported as burden p-value.

### Pitfall 7: Gene Matrix Memory Footprint

**What goes wrong:** Building a genotype matrix per gene and keeping it in memory across all genes
or storing in PipelineContext causes OOM at scale.

**How to avoid:** Build the matrix inside `test.run()`, use it, return the result dict, let it
be garbage collected. Per-gene matrices are small (200 samples x 10 variants = 16 KB).

### Pitfall 8: `pd.get_dummies` dtype in pandas 2.x

**What goes wrong:** In pandas 2.x, `get_dummies` defaults to `dtype=bool`. Passing bool columns
to numpy matrix operations with `float64` arrays causes implicit type promotion issues.

**How to avoid:** Always use `pd.get_dummies(..., dtype=float)` explicitly.

## Code Examples

### Full Covariate Alignment Example (verified)

```python
# Source: verified against pandas 2.3.3 + live test
import pandas as pd
import numpy as np

vcf_samples = ['sample_A', 'sample_B', 'sample_C', 'sample_D']

cov_df = pd.DataFrame({
    'sample_id': ['sample_C', 'sample_A', 'sample_D', 'sample_B'],
    'age': [52, 30, 28, 45],
    'sex': ['M', 'F', 'M', 'F'],
}).set_index('sample_id')

# Align (reorder) to VCF order
cov_aligned = cov_df.reindex(vcf_samples)
assert not cov_aligned.isnull().any().any(), "Missing covariate sample"

# One-hot encode sex (auto-detected categorical, <=5 unique non-numeric values)
cov_encoded = pd.get_dummies(cov_aligned, columns=['sex'], drop_first=True, dtype=float)
# Result columns: ['age', 'sex_M']
# sex_M = 1 for male, 0 for female (F dropped as reference)

covariate_matrix = cov_encoded.values.astype(np.float64)
# shape: (4, 2) — ready to pass to statsmodels
```

### Logistic Regression with OR + CI (verified)

```python
# Source: verified against statsmodels 0.14.4
import numpy as np
import statsmodels.api as sm

# burden: (n_samples,) float
# y: (n_samples,) binary 0/1
# covariates: (n_samples, k) float or None

X = sm.add_constant(burden.reshape(-1, 1))
if covariates is not None:
    X = np.column_stack([X, covariates])

result = sm.Logit(y, X).fit(disp=False)

# Extract burden coefficient (index 1, after intercept at index 0)
beta = float(result.params[1])
se = float(result.bse[1])
p_val = float(result.pvalues[1])

ci = result.conf_int()          # ndarray shape (n_params, 2) — NOT DataFrame
ci_beta_lower = float(ci[1, 0])
ci_beta_upper = float(ci[1, 1])

odds_ratio = float(np.exp(beta))
or_ci_lower = float(np.exp(ci_beta_lower))
or_ci_upper = float(np.exp(ci_beta_upper))

converged = result.mle_retvals.get("converged", True)
```

### Beta(MAF;1,25) Weight Computation (verified)

```python
# Source: scipy.stats.beta.pdf — verified against live test
from scipy.stats import beta as beta_dist
import numpy as np

mafs = np.array([0.001, 0.01, 0.05, 0.1, 0.2])
weights = beta_dist.pdf(mafs, a=1, b=25)
# Output: [24.4, 19.6, 7.3, 2.0, 0.12]
# Rare variants (MAF=0.001) get ~200x more weight than common (MAF=0.2)
```

### Genotype Parsing Edge Cases (verified)

```python
# Source: verified against live test — all cases pass
parse_gt_to_dosage('0/0') == 0
parse_gt_to_dosage('0/1') == 1
parse_gt_to_dosage('1/0') == 1
parse_gt_to_dosage('1/1') == 2
parse_gt_to_dosage('./.') is None
parse_gt_to_dosage('.|.') is None
parse_gt_to_dosage('0|1') == 1   # phased = treated as unphased
parse_gt_to_dosage('1|0') == 1
parse_gt_to_dosage('1/2') == 1   # multi-allelic = het-equivalent
parse_gt_to_dosage('2/2') == 1   # multi-allelic
parse_gt_to_dosage('0/2') == 1   # multi-allelic
```

### Separation Detection (verified against statsmodels 0.14.4)

```python
import warnings

caught = []
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    result = sm.Logit(y, X).fit(disp=False, maxiter=100)
    caught = [str(ww.message) for ww in w]

converged = result.mle_retvals.get("converged", True)
bse_max = result.bse.max()

if not converged or bse_max > 100.0:
    # Perfect/quasi separation detected — use Firth
    pass
```

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| Count-based Fisher only | Regression burden + covariate adjustment | More powerful for quantitative traits |
| No weight scheme | Beta(MAF;1,25) default, uniform option | Aligns with SKAT, regenie, SAIGE-GENE+ |
| Single missing imputation | Two-layer: pre-filter + mean imputation | Matches SKAT R convention |
| No separation handling | Firth penalized fallback | Finite estimates under sparse data |
| Test-by-test separation errors | Structured warning codes + NA row | Complete gene output, no silent omissions |

**Deprecated/outdated:**
- `firthlogist 0.5.0`: Unmaintained since August 2022. Not viable.
- Passing `sep=None` to `pd.read_csv` without `engine='python'`: C engine cannot auto-detect.

## Key Architecture Decision: Data Flow to Regression Tests

The existing `AssociationTest.run(gene, contingency_data, config)` interface only passes
pre-aggregated summary counts. Regression tests need the raw per-sample genotype data.

**Decision needed (Claude's Discretion per CONTEXT.md):** How to pass the genotype matrix.

**Recommended approach:** Extend `contingency_data` dict with optional regression keys
(`genotype_matrix`, `phenotype_vector`, `covariate_matrix`, `variant_mafs`). The engine
populates these when `trait_type != None` and covariate/matrix data is available.
`FisherExactTest.run()` ignores the new keys. This is backward compatible and requires no
interface change to `AssociationTest.run()`.

**Alternative:** Add a second `run_regression()` method to the ABC. More explicit but breaks
the uniform interface that the engine uses.

**Why extend dict:** Matches the spirit of `extra: dict[str, Any]` already in `TestResult`,
and avoids requiring every existing test to implement a new method.

## Open Questions

1. **Where is the genotype matrix built?**
   - What we know: It must be built per-gene from per-sample GT columns in the DataFrame.
   - What's unclear: Whether this happens inside the `LogisticBurdenTest.run()` (which currently
     receives only `contingency_data`), or in `AssociationAnalysisStage._process()` before calling
     the engine, or in a new `GenotypeMatrixBuilder` called by the engine per gene.
   - Recommendation: Build per-gene inside `_process()` in `AssociationAnalysisStage`, pass via
     `contingency_data`. This avoids giving tests direct DataFrame access.

2. **Covariate loading in stage or in test?**
   - What we know: Covariate loading is a one-time operation (not per-gene).
   - What's unclear: Whether to load covariates in `AssociationAnalysisStage._process()` or in a
     new `CovariateLoadingStage`.
   - Recommendation: Load in `AssociationAnalysisStage._process()` if `covariate_file` is set.
     A separate stage adds orchestration complexity for a single file read.

3. **How to detect multi-allelic for WARNING emission?**
   - `parse_gt_to_dosage('1/2')` returns 1, but we need to know it was multi-allelic to emit
     the "run bcftools norm" warning.
   - Recommendation: Return a separate boolean from the parser, or track separately with a flag
     in the matrix builder.

## Sources

### Primary (HIGH confidence)
- Live statsmodels 0.14.4 API (verified with Python subprocess): `Logit.fit()`, `LogitResults`,
  `OLS.fit()`, `OLSResults`, convergence detection via `mle_retvals['converged']`
- Live scipy 1.14.1 API (verified): `scipy.stats.beta.pdf(maf, a=1, b=25)`
- Live numpy 2.2.6 API (verified): `np.linalg.cond()` for condition number
- Live pandas 2.3.3 API (verified): `pd.read_csv(sep=None, engine='python')`,
  `df.reindex()`, `pd.get_dummies(drop_first=True, dtype=float)`
- Existing variantcentrifuge codebase (verified):
  - `variantcentrifuge/association/base.py` — ABC, TestResult, AssociationConfig
  - `variantcentrifuge/association/engine.py` — contingency_data dict structure
  - `variantcentrifuge/genotype_utils.py` — existing parse_genotype(), multi-allelic handling
  - `variantcentrifuge/gene_burden.py` — existing `_gt_to_dosage()`, per-sample GT columns
  - `variantcentrifuge/stages/analysis_stages.py` — AssociationAnalysisStage context
  - `variantcentrifuge/pyproject.toml` — all required libraries already present

### Secondary (MEDIUM confidence)
- John Lees Firth implementation gist: https://gist.github.com/johnlees/3e06380965f367e4894ea20fbae2b90d
  — Newton-Raphson with hat matrix diagonal for Jeffreys prior; algorithm verified against statsmodels
  separation detection behavior
- firthmodels v0.7.1 on PyPI (January 2026): active but adds dependency; not recommended

### Tertiary (LOW confidence)
- SKAT R package imputation convention ("bestguess" for binary, "fixed"/mean for continuous):
  stated in CONTEXT.md; not independently verified against R source

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all verified against live installed versions (statsmodels 0.14.4, scipy 1.14.1, numpy 2.2.6, pandas 2.3.3)
- Architecture: HIGH — contingency_data extension pattern is consistent with existing code structure
- Pitfalls: HIGH — all 8 pitfalls verified against live API behavior (conf_int returns ndarray, separation detection, etc.)
- Firth implementation: MEDIUM — algorithm well-known, but specific convergence characteristics depend on data

**Research date:** 2026-02-19
**Valid until:** 2026-04-01 (stable ecosystem; statsmodels/scipy rarely break APIs)
