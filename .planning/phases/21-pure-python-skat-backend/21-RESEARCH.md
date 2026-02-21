# Phase 21: Pure Python SKAT Backend - Research

**Researched:** 2026-02-21
**Domain:** Statistical genetics / quadratic form p-value computation / Python C extension packaging
**Confidence:** HIGH (primary decisions from CONTEXT.md verified against authoritative sources)

## Summary

Phase 21 implements a pure Python SKAT backend that mirrors the R SKAT package's mathematical
formulation closely enough to produce p-values within tiered tolerance thresholds. The three
primary technical challenges are: (1) correctly implementing the score test statistic and null
model in Python using statsmodels/scipy to match R's GLM residuals and LAPACK eigendecomposition;
(2) bundling the Davies qfc.cpp C++ source (from CompQuadForm CRAN) via CFFI API mode with
graceful fallback; and (3) implementing the Kuonen saddlepoint (Lugannani-Rice) and Liu
moment-matching fallbacks that match the GENESIS/GMMAT hybrid triggering strategy.

The codebase already has the full SKATBackend ABC, NullModelResult dataclass, RSKATBackend
as oracle, and a test infrastructure with mocked backends. Phase 21 slots in as a new file
`backends/python_backend.py` plus `backends/davies.py` (p-value layer), with
`tests/skat_python.py` as the engine-layer wrapper. The build system must migrate from
hatchling to setuptools to support cffi_modules C extension compilation.

**Primary recommendation:** Implement in sub-plan order: Liu fallback first (pure scipy, always
available), then Davies CFFI layer, then full PythonSKATBackend, then SKAT-O rho grid search.
This order lets each sub-plan be independently testable with no external dependencies.

---

## Standard Stack

All packages below are already installed in the project environment or required as new dependencies.

### Core (already in pyproject.toml)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| numpy | 2.2.6 | Matrix ops, eigenvalue filtering | Required throughout |
| scipy | 1.14.1 | `linalg.eigh(driver='evr')`, `integrate.quad`, `stats.chi2` | Matches R's LAPACK/QUADPACK |
| statsmodels | 0.14.4 | GLM null model, `resid_response` | R-matching GLM residuals |

### New Dependencies
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| cffi | 2.0.0 | CFFI API mode compilation of qfc.cpp | argon2-cffi / cryptography pattern; already installed |

### Build System Change
| From | To | Why |
|------|----|-----|
| hatchling | setuptools | hatchling has no C extension support; setuptools is the standard for cffi_modules |

**Installation (new runtime dependency):**
```bash
# cffi is already installed (2.0.0 confirmed in environment)
# No new pip installs needed for runtime beyond cffi
# Build-time: setuptools, cffi>=1.0.0
```

---

## Architecture Patterns

### Recommended Project Structure
```
variantcentrifuge/association/
├── backends/
│   ├── __init__.py         # get_skat_backend() factory — update "python" case
│   ├── base.py             # SKATBackend ABC, NullModelResult (unchanged)
│   ├── r_backend.py        # RSKATBackend (unchanged, Phase 20)
│   ├── python_backend.py   # PythonSKATBackend (Phase 21-02)
│   └── davies.py           # p-value layer: Davies + saddlepoint + Liu (Phase 21-01)
├── tests/
│   ├── skat_r.py           # RSKATTest (unchanged)
│   └── skat_python.py      # PurePythonSKATTest (Phase 21-02)
└── data/
    └── qfc.cpp             # Bundled CompQuadForm C++ source
variantcentrifuge/_davies_build.py  # CFFI ffibuilder for qfc.cpp
tests/unit/
├── test_davies_pvalue.py   # Unit tests for p-value layer
├── test_skat_python_backend.py  # Unit tests for PythonSKATBackend
└── test_skat_python_comparison.py  # Validation: Python vs R (50+ genes)
```

### Pattern 1: Liu Moment-Matching (last resort, always available)

**What:** Match first 4 cumulants of the weighted chi-squared sum to a scaled non-central
chi-squared distribution. Liu et al. 2009.

**When to use:** Both Davies and saddlepoint have failed, or as development scaffold.

**Algorithm (directly from R source `Get_Liu_Params_Mod_Lambda`):**
```python
# Source: https://rdrr.io/cran/SKAT/src/R/Function.R
def _liu_pvalue(q: float, lambdas: np.ndarray) -> float:
    """Liu moment-matching p-value for Q ~ sum(lambda_j * chi2_1)."""
    from scipy.stats import ncx2

    # Cumulants from eigenvalues (central chi2, df=1 each)
    c = np.array([
        np.sum(lambdas),           # c1
        np.sum(lambdas**2),        # c2
        np.sum(lambdas**3),        # c3
        np.sum(lambdas**4),        # c4
    ])
    mu_q = c[0]
    sigma_q = np.sqrt(2.0 * c[1])
    s1 = c[2] / c[1] ** 1.5
    s2 = c[3] / c[1] ** 2.0

    if s1 ** 2 > s2:
        a = 1.0 / (s1 - np.sqrt(s1 ** 2 - s2))
        d = s1 * a ** 3 - a ** 2
        l = a ** 2 - 2.0 * d
    else:
        l = 1.0 / s2
        a = np.sqrt(l)
        d = 0.0

    mu_x = l + d
    sigma_x = np.sqrt(2.0) * a

    # Standardise Q to match scaled noncentralchi2(l, ncp=d)
    q_norm = (q - mu_q) / sigma_q * sigma_x + mu_x
    p = ncx2.sf(q_norm, df=l, nc=d)
    return float(p)
```

### Pattern 2: Kuonen Saddlepoint (secondary, ~30 lines of math)

**What:** Lugannani-Rice saddlepoint approximation using CGF of weighted chi-squared sum.
K(t) = -1/2 * sum log(1 - 2t*lambda_j). Newton-Raphson (or brentq) solves K'(t_hat) = q.

**When to use:** Davies ifault > 0, p >= 1.0, p < 1000*eps, or p <= 1e-5.

**Example (Kuonen 1999 Lugannani-Rice):**
```python
# Source: GENESIS variantSetTests.R + Kuonen 1999 derivation
from scipy.stats import norm
from scipy.optimize import brentq

def _kuonen_pvalue(q: float, lambdas: np.ndarray) -> float | None:
    """Lugannani-Rice saddlepoint p-value."""
    lam_pos = lambdas[lambdas > 0]
    if len(lam_pos) == 0:
        return None
    t_max = 0.4999 / lam_pos.max()  # slightly inside singularity

    def K(t):
        return -0.5 * np.sum(np.log(1.0 - 2.0 * t * lambdas))

    def Kprime(t):
        return np.sum(lambdas / (1.0 - 2.0 * t * lambdas))

    def Kdprime(t):
        return 2.0 * np.sum(lambdas ** 2 / (1.0 - 2.0 * t * lambdas) ** 2)

    # Check q is in reachable range
    mu = Kprime(0.0)  # = sum(lambdas)
    if q <= mu:
        return None  # below mean; saddlepoint not valid for lower tail here

    try:
        t_hat = brentq(lambda t: Kprime(t) - q, 1e-10, t_max, xtol=1e-12)
    except ValueError:
        return None

    w = np.sign(t_hat) * np.sqrt(2.0 * (t_hat * q - K(t_hat)))
    u = t_hat * np.sqrt(Kdprime(t_hat))

    # Near-singularity guard (Taylor expansion when u ≈ w)
    if abs(u) < 1e-8:
        return float(norm.sf(w))

    p = float(norm.sf(w + np.log(u / w) / w))
    return p if 0.0 < p < 1.0 else None
```

### Pattern 3: Davies via CFFI (primary method)

**What:** Compile qfc.cpp (from CompQuadForm) at wheel-build time via CFFI API mode.
At runtime: try/except ImportError falls back to Liu.

**qfc() C++ signature (from CompQuadForm src/qfc.cpp):**
```c
// Parameters: lb1=eigenvalues, nc1=noncentrality(zeros), n1=df(ones),
//             r1=count, sigma=0.0, c1=test_stat, lim1=100000, acc=1e-9,
//             trace[7]=diagnostics, ifault=error_code, res=p_value
void qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma,
         double *c1, int *lim1, double *acc, double* trace,
         int* ifault, double *res);
```

**davies.py public API:**
```python
# Source: CONTEXT.md decisions + verified qfc() signature
def davies_pvalue(q: float, lambdas: np.ndarray, acc: float = 1e-9,
                  lim: int = 100_000) -> tuple[float | None, int]:
    """
    Compute P[Q > q] using Davies qfc C extension.

    Returns (p_value, ifault). ifault=0 means success.
    p_value is None if Davies unavailable or ifault > 0.
    """
    ...

def compute_pvalue(q: float, lambdas: np.ndarray) -> tuple[float, str, bool]:
    """
    Full fallback chain: Davies -> saddlepoint -> Liu.

    Returns (p_value, p_method, p_converged).
    p_method: 'davies' | 'saddlepoint' | 'liu'
    p_converged: True only if Davies succeeded (ifault=0)
    """
    ...
```

**CFFI build script (`_davies_build.py`):**
```python
# Pattern from argon2-cffi-bindings, cffi 2.0 compatible
from cffi import FFI
import os

ffibuilder = FFI()
ffibuilder.cdef("""
    void qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma,
             double *c1, int *lim1, double *acc, double* trace,
             int* ifault, double *res);
""")
_src = os.path.join(os.path.dirname(__file__), "variantcentrifuge", "association", "data", "qfc.cpp")
ffibuilder.set_source(
    "variantcentrifuge._qfc",
    '#include "qfc.h"',
    sources=[_src],
    source_extension='.cpp',
    libraries=[],
)
if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
```

**OptionalBuildExt pattern (setup.py):**
```python
from setuptools import setup
from setuptools.command.build_ext import build_ext

class OptionalBuildExt(build_ext):
    def build_extension(self, ext):
        try:
            super().build_extension(ext)
        except Exception as e:
            import warnings
            warnings.warn(f"Davies extension build failed: {e}. Will use Liu fallback.")

setup(
    setup_requires=["cffi>=1.0.0"],
    cffi_modules=["_davies_build.py:ffibuilder"],
    cmdclass={"build_ext": OptionalBuildExt},
)
```

### Pattern 4: PythonSKATBackend Score Test

**What:** SKAT score statistic Q = r^T K r where K = Z W^2 Z^T (linear weighted kernel),
W = diag(sqrt(beta_maf_weights)), Z is n x m genotype matrix, r is null model residuals.

**Critical R-matching details (from CONTEXT.md + SKAT source):**
```python
# Source: CONTEXT.md decisions verified against SKAT R source (leelabsg/SKAT)

def _compute_skat_stat(
    Z: np.ndarray,      # (n_samples, n_variants) genotype matrix
    weights: np.ndarray,  # (n_variants,) from beta_maf_weights()
    residuals: np.ndarray,  # (n_samples,) null model response residuals
    trait_type: str,
) -> tuple[float, np.ndarray]:
    """Compute Q and eigenvalues of the score matrix."""
    # Step 1: weighted genotype matrix
    W = np.diag(weights)
    ZW = Z @ W  # (n_samples, n_variants)

    # Step 2: Score statistic Q = r^T (ZW)(ZW)^T r / 2
    # (division by 2 matches R convention for linear kernel)
    Zr = ZW.T @ residuals  # (n_variants,)
    Q = float(Zr @ Zr) / 2.0

    # Step 3: Residual variance for linear trait (binary: sigma2=1)
    sigma2 = np.var(residuals, ddof=1) if trait_type == "quantitative" else 1.0

    # Step 4: Kernel matrix eigenvalues for p-value computation
    # R uses K = ZW P ZW^T where P = I - X(X^T X)^-1 X^T
    # With intercept-only null: P ~ I - 1*1^T/n
    # Use scipy.linalg.eigh(driver='evr') to match R's DSYEVR
    K = ZW @ ZW.T  # (n_samples, n_samples), only eigenvalues needed
    lambdas_all = scipy.linalg.eigh(K / (2.0 * sigma2), eigvals_only=True, driver='evr')

    # Step 5: Eigenvalue filtering — EXACT R threshold
    # Source: SKAT R/Function.R Get_Lambda()
    pos_mask = lambdas_all >= 0
    if not pos_mask.any():
        return Q, np.array([])
    threshold = lambdas_all[pos_mask].mean() / 100_000.0
    lambdas = lambdas_all[lambdas_all > threshold]

    return Q, lambdas
```

**Null model fitting (statsmodels GLM):**
```python
# Source: statsmodels GLMResults docs + CONTEXT.md decision
import statsmodels.api as sm

def _fit_null_model(phenotype, covariates, trait_type):
    X = sm.add_constant(covariates) if covariates is not None else np.ones((len(phenotype), 1))
    if trait_type == "binary":
        family = sm.families.Binomial()
    else:
        family = sm.families.Gaussian()
    model = sm.GLM(phenotype, X, family=family)
    result = model.fit()
    # resid_response = endog - fittedvalues (y - mu_hat)
    # This matches R's SKAT_Null_Model response residuals
    residuals = result.resid_response
    return result, residuals
```

### Pattern 5: SKAT-O Rho Grid Search

**What:** Unified test Q_rho = (1-rho)*Q_SKAT + rho*Q_Burden over rho grid.
Optimal rho minimizes p-value (or equivalently maximizes test statistic).

**Fixed rho grid (matches R exactly):**
```python
# Source: CONTEXT.md decisions + leelabsg/SKAT SKAT-O documentation
SKATO_RHO_GRID = [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0]

# SKAT-O p-value integration:
from scipy.integrate import quad
p_skato, _ = quad(integrand_func, 0, np.inf, epsabs=1.49e-8)
# epsabs=1.49e-8 matches R's QUADPACK default tolerance
```

### Anti-Patterns to Avoid
- **Using resid_deviance or resid_pearson:** R uses response residuals (y - mu). Using other
  residual types breaks R agreement.
- **scipy.linalg.eigh with driver='ev':** Use `driver='evr'` to match R's DSYEVR algorithm.
- **Negative eigenvalue threshold:** Use `max(eigenvalues, 0)` only for the filtering step,
  not as a clamp. The R threshold is `lambda > mean(positive_lambdas) / 100_000`.
- **acc=1e-6 for Davies:** R SKAT default is 1e-6 but GENESIS/SKATh use 1e-9. Use 1e-9 per
  CONTEXT.md decision for tighter accuracy.
- **cffi ABI mode:** Use API mode (set_source + cdef). ABI mode is slower and harder to
  distribute.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| GLM null model | Custom logistic regression | `statsmodels.api.GLM` | R agreement on residuals |
| Eigendecomposition | numpy.linalg.eigh or custom | `scipy.linalg.eigh(driver='evr')` | Matches R's DSYEVR exactly |
| Davies integration | Pure Python characteristic function inversion | `qfc.cpp` via CFFI | Same binary as R CompQuadForm |
| Numerical integration | Trapezoid/Riemann for SKAT-O | `scipy.integrate.quad(epsabs=1.49e-8)` | Matches R's QUADPACK |
| Beta MAF weights | Custom MAF weight formula | `variantcentrifuge.association.weights.beta_maf_weights` | Already implemented in Phase 19 |

**Key insight:** SKAT's mathematical correctness depends on using the exact same LAPACK driver
(DSYEVR) and the same C code (qfc.cpp) as R does. Any pure-Python reimplementation of these
numerical kernels will diverge in the extreme tails, which is precisely where GWAS associations
of interest lie.

---

## Common Pitfalls

### Pitfall 1: Response vs. Deviance Residuals
**What goes wrong:** Using `resid_deviance` or `resid_pearson` instead of `resid_response`.
**Why it happens:** Deviance residuals are the statsmodels default for diagnostics.
**How to avoid:** Always use `result.resid_response` which equals `endog - fittedvalues` (y - mu_hat).
**Warning signs:** p-values that are plausible but consistently off by 2-5x from R oracle.

### Pitfall 2: Davies False Convergence (ifault=0 but Wrong p-value)
**What goes wrong:** Davies returns ifault=0 but p-value near 0 or above 1.
**Why it happens:** Numerical cancellation near the singularity at t -> 1/(2*lambda_max).
**How to avoid:** Implement all four GENESIS+GMMAT triggering conditions:
  1. `ifault > 0`
  2. `p >= 1.0`
  3. `p < 1000 * np.finfo(float).eps` (~ 2.22e-13)
  4. `p <= 1e-5` (proactive saddlepoint — catches false convergence zone)
**Warning signs:** Genes with very extreme p-values (< 1e-10) where saddlepoint would give different result.

### Pitfall 3: cffi_modules Not Supported by hatchling
**What goes wrong:** `pyproject.toml` with `build-backend = "hatchling.build"` ignores cffi_modules.
**Why it happens:** hatchling is a pure-Python build backend with no C extension hook.
**How to avoid:** Migrate build-backend to setuptools. Verified: the project has no hatchling
  plugins that would break migration.
**Warning signs:** Wheel builds succeed but `_qfc` module is never compiled.

### Pitfall 4: cffi 2.0 API Changes
**What goes wrong:** Code written for cffi 1.x may use deprecated `ffi.verify()`.
**Why it happens:** `ffi.verify()` (API level + in-line mode) was deprecated in 1.x and removed in 2.0.
**How to avoid:** Use out-of-line API mode: `ffibuilder.set_source()` + `ffibuilder.cdef()`.
  The project already has cffi 2.0.0 installed — no 1.x patterns will work.
**Warning signs:** ImportError with "verify" not found.

### Pitfall 5: Matrix Rank Computation Timing
**What goes wrong:** `numpy.linalg.matrix_rank` returns 0 after eigenvalue filtering, skipping
  valid genes.
**Why it happens:** Rank is checked after filtering, but the filter itself can reduce rank.
**How to avoid:** Check `matrix_rank >= 2` on the FULL eigenvalue set (before filtering),
  not on the filtered subset.
**Warning signs:** Genes with 3+ variants being skipped with `rank_deficient`.

### Pitfall 6: SKAT-O Before SKAT Validates
**What goes wrong:** Implementing SKAT-O (rho grid search) before basic SKAT p-values match R.
**Why it happens:** SKAT-O depends on SKAT and Burden statistics being correct.
**How to avoid:** Sub-plan 21-02 explicitly gates SKAT-O behind SKAT validation per CONTEXT.md.
**Warning signs:** SKAT-O rho=0 result differs from pure SKAT result.

### Pitfall 7: Liu Anti-Conservative Tail Behavior
**What goes wrong:** Using Liu as primary method produces inflated type I error at genome-wide
  significance (alpha ~ 1e-8).
**Why it happens:** Liu moment-matching fails in extreme tails (known since Wu et al. 2011).
  Can inflate type I error ~25x at alpha=1e-8 per CONTEXT.md.
**How to avoid:** Liu is last resort only. Saddlepoint must catch p <= 1e-5 before Liu is used.
**Warning signs:** p_method = "liu" for genes with small p-values (< 1e-4).

### Pitfall 8: pyproject.toml cffi_modules Not a Standard Key
**What goes wrong:** `cffi_modules` in `[tool.setuptools]` section is not recognized.
**Why it happens:** cffi_modules is a setup() parameter, not a pyproject.toml table key.
  It requires setup.py (with cffi_modules= argument) or setup.cfg.
**How to avoid:** Keep a minimal `setup.py` alongside `pyproject.toml`:
  ```python
  from setuptools import setup
  setup(cffi_modules=["_davies_build.py:ffibuilder"])
  ```
**Warning signs:** Extension compiles locally but not in CI or pip install.

---

## Code Examples

### Complete davies.py Fallback Chain
```python
# Source: CONTEXT.md + GENESIS variantSetTests.R + verified Python math

import logging
import sys
import numpy as np

logger = logging.getLogger("variantcentrifuge")

# Lazy import: try compiled extension first, fall back to Liu-only mode
_qfc = None
_DAVIES_AVAILABLE = False

def _try_load_davies() -> bool:
    global _qfc, _DAVIES_AVAILABLE
    if _DAVIES_AVAILABLE:
        return True
    if os.environ.get("VARIANTCENTRIFUGE_NO_C_EXT"):
        return False
    try:
        from variantcentrifuge import _qfc as _qfc_mod
        _qfc = _qfc_mod
        _DAVIES_AVAILABLE = True
        return True
    except ImportError:
        return False

_EPS1000 = 1000.0 * np.finfo(float).eps  # ~2.22e-13
_PROACTIVE_THRESHOLD = 1e-5  # GMMAT-style proactive saddlepoint

def compute_pvalue(
    q: float,
    lambdas: np.ndarray,
    acc: float = 1e-9,
    lim: int = 100_000,
) -> tuple[float, str, bool]:
    """
    Three-tier fallback: Davies -> saddlepoint -> Liu.

    Returns
    -------
    (p_value, p_method, p_converged)
    """
    if len(lambdas) == 0:
        return float("nan"), "liu", False

    # Tier 1: Davies (requires compiled _qfc)
    if _try_load_davies():
        p_davies, ifault = _davies_call(q, lambdas, acc, lim)
        needs_fallback = (
            ifault > 0
            or p_davies is None
            or p_davies >= 1.0
            or p_davies < _EPS1000
            or p_davies <= _PROACTIVE_THRESHOLD
        )
        if not needs_fallback:
            return p_davies, "davies", True
        logger.info(
            f"Davies fallback triggered: ifault={ifault}, p={p_davies} "
            f"(thresholds: >=1.0 | <{_EPS1000:.2e} | <={_PROACTIVE_THRESHOLD})"
        )

    # Tier 2: Kuonen saddlepoint
    p_sp = _kuonen_pvalue(q, lambdas)
    if p_sp is not None and 0.0 < p_sp < 1.0:
        return p_sp, "saddlepoint", False

    # Tier 3: Liu moment-matching (last resort)
    logger.info("Saddlepoint failed, using Liu moment-matching (last resort)")
    return _liu_pvalue(q, lambdas), "liu", False
```

### Eigenvalue Computation Matching R
```python
# Source: SKAT R/Function.R Get_Lambda() + CONTEXT.md decision

import scipy.linalg

def compute_eigenvalues(
    Z: np.ndarray,           # (n_samples, n_variants) weighted genotype matrix
    residuals: np.ndarray,   # (n_samples,) response residuals
    sigma2: float,           # residual variance (1.0 for binary)
) -> np.ndarray:
    """
    Compute eigenvalues of score kernel matrix for p-value computation.
    Uses scipy.linalg.eigh(driver='evr') to match R's DSYEVR.
    """
    n = len(residuals)

    # P matrix = I - X(X^TX)^{-1}X^T (simplified for intercept-only null)
    # K = Z^T P Z weighted by W
    # Efficient: compute ZZ^T eigenvalues (same nonzero spectrum as Z^TZ)
    K = Z @ Z.T  # (n_samples, n_samples)

    # Scale by residual variance
    K_scaled = K / (2.0 * sigma2) if sigma2 > 0 else K / 2.0

    # Eigendecomposition: driver='evr' matches DSYEVR used by R's eigen()
    lambdas_all = scipy.linalg.eigh(K_scaled, eigvals_only=True, driver='evr')

    # Threshold: keep eigenvalues > mean(positive_eigenvalues) / 100_000
    # EXACT match to R: Get_Lambda() in Function.R
    pos = lambdas_all[lambdas_all >= 0]
    if len(pos) == 0:
        return np.array([])
    threshold = pos.mean() / 100_000.0
    return lambdas_all[lambdas_all > threshold]
```

### PythonSKATBackend skeleton
```python
# Source: SKATBackend ABC in backends/base.py + CONTEXT.md

from variantcentrifuge.association.backends.base import NullModelResult, SKATBackend

class PythonSKATBackend(SKATBackend):
    """Pure Python SKAT backend (Phase 21)."""

    def detect_environment(self) -> None:
        """Verify numpy, scipy, statsmodels are importable. Always succeeds."""
        import numpy, scipy, statsmodels  # noqa: F401
        from variantcentrifuge.association.backends.davies import _try_load_davies
        if not _try_load_davies():
            logger.info("Davies C extension unavailable; Liu fallback active")

    def fit_null_model(self, phenotype, covariates, trait_type) -> NullModelResult:
        import statsmodels.api as sm
        X = sm.add_constant(covariates) if covariates is not None else np.ones((len(phenotype), 1))
        family = sm.families.Binomial() if trait_type == "binary" else sm.families.Gaussian()
        glm_result = sm.GLM(phenotype, X, family=family).fit()
        residuals = glm_result.resid_response  # y - mu_hat (response residuals)
        return NullModelResult(
            model=glm_result,
            trait_type=trait_type,
            n_samples=len(phenotype),
            adjustment=False,
            extra={"residuals": residuals},
        )

    def test_gene(self, gene, genotype_matrix, null_model, method, weights_beta) -> dict:
        from variantcentrifuge.association.weights import beta_maf_weights
        from variantcentrifuge.association.backends.davies import compute_pvalue

        Z = genotype_matrix  # (n_samples, n_variants)
        n_variants = Z.shape[1]

        # Rank check BEFORE eigenvalue filtering
        if np.linalg.matrix_rank(Z) < 2:
            return {"p_value": None, "rho": None, "n_variants": n_variants,
                    "n_marker_test": 0, "warnings": [],
                    "skip_reason": "rank_deficient", "p_method": None,
                    "p_converged": False}

        residuals = null_model.extra["residuals"]
        mafs = Z.mean(axis=0) / 2.0
        a1, a2 = weights_beta
        weights = beta_maf_weights(mafs, a=a1, b=a2)

        ZW = Z * weights[np.newaxis, :]  # element-wise broadcast
        Q, lambdas = _compute_skat_stat(ZW, residuals)

        p_value, p_method, p_converged = compute_pvalue(Q, lambdas)
        return {
            "p_value": p_value,
            "rho": None,
            "n_variants": n_variants,
            "n_marker_test": int((weights > 0).sum()),
            "warnings": [],
            "p_method": p_method,
            "p_converged": p_converged,
        }
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Davies only | Davies + saddlepoint + Liu fallback chain | GENESIS/GMMAT ~2016-2018 | False convergence detection |
| acc=1e-6 | acc=1e-9 (SKATh recommendation) | Phase 21 | Tighter convergence guard |
| verify() CFFI mode | set_source() out-of-line API mode | cffi 2.0 (2024) | verify() removed in 2.0 |
| hatchling build | setuptools + cffi_modules | Phase 21 | C extension support required |
| qfc.c (C) | qfc.cpp (C++) | CompQuadForm current CRAN | Same algorithm, C++ linkage |

**Deprecated/outdated:**
- `ffi.verify()`: Removed in cffi 2.0. Must use out-of-line API mode.
- `SKAT_Null_Model_MomentAdjust()`: Deprecated R function; SKAT 2.x uses `SKAT_Null_Model(Adjustment=TRUE)`.
- Davies with lim=10000: Current SKATh recommendation is lim=100_000 for better convergence at small p-values.

---

## Open Questions

1. **qfc.cpp vs qfc.c: C++ linkage requirement**
   - What we know: CompQuadForm src/qfc.cpp is C++ (confirmed via GitHub). The file uses C++ linkage (`extern "C"` may or may not be declared internally).
   - What's unclear: Whether extern "C" is declared in the header, or if the CFFI set_source needs `extern "C" { void qfc(...); }` in the cdef.
   - Recommendation: When downloading qfc.cpp, inspect for extern "C". If absent, wrap in `extern "C" {}` block in CFFI cdef or add `extern "C"` to the set_source C source snippet.

2. **cffi_modules in setup.py vs pyproject.toml**
   - What we know: `cffi_modules` is a setup() parameter; setuptools pyproject.toml `[tool.setuptools]` does not have a cffi_modules key. Must use setup.py.
   - What's unclear: Whether setuptools 70+ supports cffi_modules purely via pyproject.toml (it does not as of Feb 2026).
   - Recommendation: Keep a minimal `setup.py` alongside `pyproject.toml` that only contains `setup(cffi_modules=...)`. Everything else stays in pyproject.toml.

3. **Binary trait P matrix: projection vs. none**
   - What we know: For SKAT with binary traits, R uses modified residuals that account for the projected P matrix. For the intercept-only null, this simplifies but is not exactly I.
   - What's unclear: Whether the Python implementation needs to explicitly construct P = I - H for covariate case or if statsmodels GLM residuals already encode this.
   - Recommendation: statsmodels resid_response = y - mu_hat already accounts for the fitted model (covariates included in X). This should match R's residuals. Validate against R oracle in sub-plan 21-03.

---

## Sources

### Primary (HIGH confidence)
- SKAT R source `R/Function.R` on rdrr.io — Liu algorithm with exact cumulant formulas and eigenvalue threshold `mean(lambda)/100000`
- statsmodels GLMResults documentation — `resid_response = endog - fittedvalues` confirmed
- scipy.linalg.eigh documentation — `driver='evr'` maps to DSYEVR; confirmed in environment (scipy 1.14.1)
- CompQuadForm GitHub `src/qfc.cpp` — function signature `void qfc(double* lb1, ...)` confirmed
- cffi 2.0 documentation — `set_source()` out-of-line API mode; `verify()` removed
- GENESIS `R/variantSetTests.R` — Davies fallback conditions: ifault>0, p<1000*eps, p>1, p<=1e-5
- Phase 21 CONTEXT.md — All locked decisions confirmed against sources above

### Secondary (MEDIUM confidence)
- GENESIS variantSetTests.R via GitHub (WebFetch verified) — Davies acc=1e-9 and CompQuadForm package call pattern
- setuptools documentation — `optional=True` on Extension marks as non-fatal; `cffi_modules` requires setup.py
- cffi 2.0 whatsnew (WebFetch verified) — No breaking changes to set_source/compile API; only Python 3.8 dropped

### Tertiary (LOW confidence)
- momentchi2py PyPI — Alternative Liu/HBE implementation; not recommended over inline math (adds dependency for ~10 lines of math)

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all packages verified in environment (numpy 2.2.6, scipy 1.14.1, statsmodels 0.14.4, cffi 2.0.0)
- Architecture: HIGH — SKATBackend ABC, NullModelResult, get_skat_backend() already exist; patterns confirmed
- Liu algorithm: HIGH — exact R source code obtained from rdrr.io, verified Python implementation computes correctly
- Kuonen saddlepoint: HIGH — Lugannani-Rice formula verified in Python, math confirmed numerically
- Davies CFFI: MEDIUM — qfc.cpp signature confirmed, but C++ linkage details (extern "C") need verification at extraction time
- Pitfalls: HIGH — sourced from GENESIS R code, cffi 2.0 docs, and SKAT R source
- Build system: MEDIUM — cffi_modules + setuptools pattern is standard but exact pyproject.toml co-existence details need testing

**Research date:** 2026-02-21
**Valid until:** 2026-03-21 (cffi, setuptools, scipy stable; 30-day estimate)
