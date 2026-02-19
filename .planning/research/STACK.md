# Technology Stack: Modular Rare Variant Association Framework

**Project:** variantcentrifuge — modular association framework milestone
**Researched:** 2026-02-19
**Scope:** Stack additions only. Existing stack (pandas, numpy, scipy, statsmodels) is validated and not re-researched.

---

## Critical Constraint: SciPy and Python 3.10

**This constraint shapes every dependency decision in this milestone.**

SciPy 1.16 dropped Python 3.10 support. SciPy 1.15.x is the last series supporting Python 3.10. The current `pyproject.toml` declares `requires-python = ">=3.10"` and pins no scipy upper bound, meaning a fresh install on Python 3.10 will resolve to scipy 1.15.x, while Python 3.11+ gets scipy 1.17.x.

**Decision required before roadmap phase 1:**

- **Option A:** Pin `scipy>=1.13,<1.16` — maintains Python 3.10 support, locks to older scipy
- **Option B:** Bump `requires-python = ">=3.11"` — cleans up the constraint, uses current scipy (1.17.0), aligns with HPC Python 3.11 prevalence

**Recommendation: Option B.** HPC clusters (the primary deployment target) have Python 3.11+ as standard by 2025. scipy 1.17.0 (released January 10, 2026) includes performance improvements relevant to linalg operations used in SKAT kernel computation. The association framework adds new mathematical complexity that benefits from a current scipy baseline.

| scipy version | Python support | Status |
|--------------|---------------|--------|
| 1.15.3 (last 3.10-compatible) | 3.10–3.13 | Maintenance only |
| 1.17.0 (current) | 3.11–3.14 | Active development |

---

## New Dependencies Required

### 1. rpy2 — R Bridge for SKAT Backend

**Recommended version:** `rpy2>=3.6.4`
**Current release:** 3.6.4 (September 26, 2025)
**Confidence:** HIGH (verified PyPI)

**What it provides:** Embeds R in the Python process and exposes R packages as Python objects. The SKAT/SKAT-O R package (leelabsg/SKAT on CRAN) is the canonical, peer-reviewed implementation with 10+ years of validation in the literature. Calling it via rpy2 gives the pure-Python fallback a validated reference point.

**Python compatibility:** 3.8–3.13 (verified via PyPI classifiers).

**R compatibility:**
- rpy2 3.6.0 fixed the `Rcomplex` C-API definition mismatch for R 4.3+ (in API mode). This was a source of subtle correctness bugs in earlier rpy2 with newer R.
- R 4.4 has one known issue: on CentOS 7 / RHEL 7 only, BLAS linker errors occur during `pip install rpy2`. Workaround: `export LD_LIBRARY_PATH=/opt/R/4.4.0/lib/R/lib` before installation. Modern Linux (Ubuntu 22+, RHEL 8/9) and macOS are unaffected. This is a build-time installation issue, not a runtime correctness bug.
- rpy2 3.6.1 fixed Windows initialization callback ordering; Windows users need 3.6.1 minimum.

**Graceful fallback — correct pattern:**

```python
# association_backends/skat_r_backend.py

_SKAT_R_AVAILABLE: bool = False
_SKAT_R_ERROR: str | None = None
_skat_r_pkg = None

try:
    import rpy2.robjects as _robjects
    from rpy2.robjects.packages import importr as _importr
    from rpy2.robjects import numpy2ri as _numpy2ri

    # importr("SKAT") triggers both R initialization AND SKAT package load.
    # Catches: R not in PATH, R not installed, SKAT R package not installed.
    _skat_r_pkg = _importr("SKAT")
    _numpy2ri.activate()
    _SKAT_R_AVAILABLE = True
except Exception as e:
    # Broad except is intentional: rpy2 raises RRuntimeError, OSError,
    # and ImportError depending on failure mode.
    _SKAT_R_ERROR = str(e)
```

**Critical:** Never use `except ImportError` alone. R initialization failures raise `OSError` (R binary not found) or `rpy2.rinterface.embedded.RRuntimeError` (R found but SKAT package missing). The broad `except Exception` is required.

**Calling SKAT from rpy2:**

```python
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

numpy2ri.activate()
skat = importr("SKAT")

# Build null model (outcome ~ covariates, no genotypes)
null_model = skat.SKAT_Null_Model(
    ro.Formula("y ~ age + sex + PC1 + PC2"),
    data=ro_dataframe,          # rpy2 DataFrame
    out_type="D",               # "D"=dichotomous, "C"=continuous
)

# SKAT-O test (optimal combination of SKAT + burden)
result = skat.SKAT(
    Z_matrix,                   # n_samples x n_variants numpy array, auto-converted
    null_model,
    method="SKATO",
    weights_beta=ro.FloatVector([1, 25]),  # MAF weighting
)
p_value = float(result.rx2("p.value")[0])
```

**Declare as optional dependency:**

```toml
[project.optional-dependencies]
association = [
    "rpy2>=3.6.4",
    "matplotlib>=3.10",
]
```

Do NOT add rpy2 to required dependencies. Users without R get the pure-Python SKAT fallback.

---

### 2. C Extension for Davies Method (qfc.c) via ctypes

**Recommendation: ctypes (stdlib), NOT cffi, NOT hatch-cython compilation.**
**Confidence:** HIGH for ctypes pattern; MEDIUM for cross-platform build strategy.

**Why Davies method matters:** The mixture-of-chi-squared p-value for SKAT is computed by inverting the characteristic function numerically. Davies' method (implemented in qfc.c, sourced from CompQuadForm R package) is the near-exact approach used by SKAT-R. The Liu moment-matching fallback (see section 3) is acceptable when qfc is unavailable, but Davies gives better accuracy for small p-values (< 1e-6).

**qfc function signature** (from CompQuadForm R package source, verified 2026-02-19):

```c
void qfc(double *lambdas,    // eigenvalue array (input)
         double *noncentral,  // non-centrality parameters (input)
         int    *df,          // degrees of freedom array (input)
         int    *r,           // number of terms/eigenvalues (input)
         double *sigma,       // sigma parameter (input)
         double *q,           // quantile / test statistic (input)
         int    *lim,         // iteration limit (input)
         double *acc,         // accuracy threshold (input)
         double *trace,       // 7-element output trace array (output)
         int    *ifault,      // fault indicator: 0=ok (output)
         double *res);        // CDF value: p = 1 - res (output)
```

**Why ctypes over cffi:**
- ctypes is stdlib — zero additional dependency
- qfc has exactly one function with a fixed, well-documented signature
- cffi advantages (C header parsing, API mode performance) irrelevant for a single stable function
- cffi's build complexity (needs gcc at import time in API mode) defeats the purpose

**Why NOT compile at install time via hatchling:**
- Requires gcc/MSVC at `pip install` time — fails silently on HPC environments without build tools
- Wheel platform tagging becomes complex (manylinux vs. platform-specific)
- One C function is not worth the CI/CD complexity of cibuildwheel

**Recommended strategy: Runtime lazy compilation with Liu fallback**

```python
# association_backends/davies_method.py

import ctypes
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

log = logging.getLogger(__name__)

_QFC_LIB: ctypes.CDLL | None = None
_QFC_AVAILABLE: bool = False


def _get_qfc_source() -> Path:
    """Locate bundled qfc.c source file."""
    import importlib.resources
    ref = importlib.resources.files("variantcentrifuge.data").joinpath("qfc.c")
    # Copy to temp dir if inside a zip (editable installs work directly)
    with importlib.resources.as_file(ref) as p:
        return p


def _try_compile_qfc() -> str | None:
    """Attempt to compile qfc.c at runtime. Returns path to shared library or None."""
    src = _get_qfc_source()
    suffix = {
        "win32": ".dll",
        "darwin": ".dylib",
    }.get(sys.platform, ".so")

    out_dir = Path(tempfile.mkdtemp())
    out = out_dir / f"qfc{suffix}"

    if sys.platform == "win32":
        cmd = ["cl.exe", "/LD", f"/Fe:{out}", str(src)]
    else:
        cmd = ["gcc", "-O2", "-shared", "-fPIC", "-o", str(out), str(src), "-lm"]

    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=30)
        return str(out)
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as e:
        log.debug("qfc.c compilation failed: %s. Will use Liu fallback.", e)
        return None


def load_qfc() -> bool:
    """Load qfc shared library. Returns True if successful."""
    global _QFC_LIB, _QFC_AVAILABLE
    if _QFC_AVAILABLE:
        return True

    lib_path = _try_compile_qfc()
    if lib_path is None:
        return False

    try:
        lib = ctypes.CDLL(lib_path)
        lib.qfc.restype = None
        lib.qfc.argtypes = [
            ctypes.POINTER(ctypes.c_double),  # lambdas
            ctypes.POINTER(ctypes.c_double),  # noncentral
            ctypes.POINTER(ctypes.c_int),     # df
            ctypes.POINTER(ctypes.c_int),     # r (n_terms)
            ctypes.POINTER(ctypes.c_double),  # sigma
            ctypes.POINTER(ctypes.c_double),  # q (test statistic)
            ctypes.POINTER(ctypes.c_int),     # lim
            ctypes.POINTER(ctypes.c_double),  # acc
            ctypes.POINTER(ctypes.c_double),  # trace[7]
            ctypes.POINTER(ctypes.c_int),     # ifault
            ctypes.POINTER(ctypes.c_double),  # res
        ]
        _QFC_LIB = lib
        _QFC_AVAILABLE = True
        return True
    except OSError as e:
        log.debug("Failed to load compiled qfc: %s", e)
        return False
```

**Ship qfc.c in the package under `variantcentrifuge/data/`:**

```toml
# pyproject.toml
[tool.hatch.build.targets.wheel]
packages = ["variantcentrifuge"]
artifacts = ["variantcentrifuge/data/*.c"]
```

**Graceful degradation path:**
1. Try Davies (qfc.c via ctypes) — most accurate
2. Fall back to Liu moment-matching (pure scipy) — acceptable for p > 1e-8
3. Log which method was used in result metadata

---

### 3. scipy.stats — Liu Moment-Matching (No New Dependency)

**No new scipy dependency required.** scipy is already required.

scipy does NOT implement Davies' method or Liu moment-matching natively. There is no `scipy.stats.mixture_chi2` or equivalent. These must be implemented using scipy primitives.

**Liu 2009 moment-matching — pure scipy implementation:**

```python
import numpy as np
from scipy.stats import ncx2, norm


def liu_pvalue(q_stat: float, eigenvalues: np.ndarray) -> float:
    """
    Liu 2009 moment-matching approximation for mixture-of-chi-squared p-values.
    Used as fallback when Davies/qfc compilation is unavailable.

    Reference: Liu H et al. (2009) Am J Hum Genet 85:726-733.
    """
    lam = np.asarray(eigenvalues, dtype=float)
    c1 = lam.sum()
    c2 = (lam**2).sum()
    c3 = (lam**3).sum()
    c4 = (lam**4).sum()

    s1 = c3 / (c2**1.5)
    s2 = c4 / (c2**2)

    if s1**2 > s2:
        # Non-central chi-squared approximation
        a = 1.0 / (s1 - np.sqrt(s1**2 - s2))
        d = s1 * a**3 - a**2
        l = a**2 - 2 * d
        # Standardize test statistic to matched distribution
        q_std = (q_stat - (c1 - l)) / np.sqrt(2 * c2) * np.sqrt(2 * l) + d
        p = float(1.0 - ncx2.cdf(q_std, df=1, nc=d))
    else:
        # Normal approximation when moments don't match chi-squared form
        mu = c1
        sigma = np.sqrt(2 * c2)
        p = float(1.0 - norm.cdf((q_stat - mu) / sigma))

    return max(p, 1e-300)  # Clamp to avoid log(0) downstream
```

`scipy.stats.ncx2` and `scipy.stats.norm` are stable since scipy 0.9 — no version constraint beyond the existing `scipy` requirement.

**scipy.linalg for SKAT kernel eigenvalues:**

```python
from scipy.linalg import eigh

# SKAT score statistic: Q = y^T K y (where K = Z W Z^T)
# P-value requires eigenvalues of P^(1/2) K P^(1/2)
# P = projection matrix from null model
eigenvalues, _ = eigh(kernel_matrix)   # eigh for symmetric matrices (faster than eig)
eigenvalues = eigenvalues[eigenvalues > 1e-10]  # Trim near-zero eigenvalues
```

`scipy.linalg.eigh` is the correct function for symmetric/Hermitian matrices (faster than `eig`, guaranteed real eigenvalues). This is stable across scipy 1.x.

---

### 4. statsmodels — Burden Tests with Covariates (No New Dependency)

**No new statsmodels dependency required.** statsmodels 0.14.6 (released December 5, 2025) is current.

**Current stable API for logistic burden test:**

```python
import numpy as np
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit
from scipy.stats import chi2

# Build design matrix: [intercept, burden, cov1, cov2, ...]
X_full = sm.add_constant(np.column_stack([burden_score, covariates]))
X_null = sm.add_constant(covariates)

# Fit models
result_full = Logit(y_binary, X_full).fit(disp=False)
result_null = Logit(y_binary, X_null).fit(disp=False)

# Wald test for burden coefficient (index 1, after intercept)
r_matrix = np.zeros((1, X_full.shape[1]))
r_matrix[0, 1] = 1.0
wald = result_full.wald_test(r_matrix, scalar=True)
wald_pval = float(wald.pvalue)

# Likelihood ratio test (preferred for GLMs)
lrt_stat = 2.0 * (result_full.llf - result_null.llf)
lrt_pval = float(chi2.sf(lrt_stat, df=1))
```

**Known API change to be aware of:** `wald_test(scalar=True)` is the forward-compatible form. The `scalar` parameter was added in 0.14.x. Before 0.14, the test statistic was returned as a numpy array; after 0.15, the default will switch to scalar. Pass `scalar=True` explicitly in all new association code.

**Linear regression (continuous trait):**

```python
from statsmodels.regression.linear_model import OLS

X_full = sm.add_constant(np.column_stack([burden_score, covariates]))
result = OLS(y_continuous, X_full).fit()
burden_pval = result.pvalues[1]   # Index 1 = burden coefficient
burden_ci = result.conf_int()[1]  # 95% CI for burden
```

---

### 5. Build System — Hatchling (No Change for MVP)

**Stay with hatchling as-is for the MVP phase.** No hatch-cython required.

The qfc.c compilation happens at runtime via subprocess (see section 2), so no build-time C toolchain dependency is needed. Hatchling's wheel target already supports arbitrary data files via `artifacts`.

**Minimal pyproject.toml changes:**

```toml
[project.optional-dependencies]
dev = [
    # ... existing unchanged ...
]
association = [
    "rpy2>=3.6.4",
    "matplotlib>=3.10",
]

[tool.hatch.build.targets.wheel]
packages = ["variantcentrifuge"]
artifacts = ["variantcentrifuge/data/*.c"]  # Include qfc.c source
```

**Future consideration (post-MVP):** If the association framework gains traction and users report compilation failures in HPC environments, add cibuildwheel to CI to ship pre-compiled platform wheels containing a compiled `_qfc.so`. This avoids runtime gcc requirement. The module interface remains identical; only the loading path changes (load pre-built SO instead of compiling).

---

### 6. matplotlib — Diagnostic Plots (Optional)

**Recommended version:** `matplotlib>=3.10`
**Current release:** 3.10.8 (December 10, 2025)
**Python compatibility:** 3.10+ — aligned with project minimum.
**Confidence:** HIGH (verified PyPI)

**Declare as optional** under `association` extra (see section 5). Never import at module level.

**Correct lazy import pattern for HPC (no display server):**

```python
# association_plots.py

def _require_matplotlib():
    """Lazy import with backend configuration for headless environments."""
    try:
        import matplotlib
        matplotlib.use("Agg")          # Must be set before pyplot import
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        raise ImportError(
            "matplotlib is required for diagnostic plots.\n"
            "Install with: pip install variantcentrifuge[association]"
        ) from None


def plot_qq(p_values: np.ndarray, output_path: str) -> None:
    """QQ plot for p-value calibration assessment."""
    plt = _require_matplotlib()
    # ... implementation
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
```

**Do not use `import matplotlib.pyplot as plt` at module level.** Module-level pyplot import:
1. Adds ~100ms startup cost to every invocation (even `variantcentrifuge --help`)
2. Raises `_tkinter.TclError` on headless HPC nodes without `$DISPLAY`

`matplotlib.use("Agg")` must be called before `import matplotlib.pyplot` to select the non-interactive backend. It has no effect if called after pyplot is imported.

---

### 7. momentchi2 — Do Not Add

**Decision: Reject.**
**Confidence:** HIGH

momentchi2 (PyPI: `momentchi2`, GitHub: `deanbodenham/momentchi2py`) implements Hall-Buckley-Eagleson (HBE) and Lindsay-Pilla-Basak (LPB) moment-matching methods.

**Why not:**
- Last commit: August 24, 2021. Zero GitHub releases published.
- Implements moment-matching methods that overlap completely with the Liu 2009 method implemented directly via scipy (section 3)
- For SKAT p-values, Liu performs equivalently to HBE for p > 1e-8 (the range that matters clinically)
- Adds an unmaintained transitive dependency for capability already achievable in 20 lines of scipy code

**Use instead:** Direct Liu moment-matching implementation (section 3) using `scipy.stats.ncx2` and `scipy.stats.norm`. This is what the SKAT R package itself uses as its Davies fallback.

---

## What NOT to Add

| Library | Why Not |
|---------|---------|
| momentchi2 | Unmaintained since 2021; Liu via scipy is equivalent |
| cffi | No advantage over ctypes for one stable C function |
| hatch-cython | Adds build-time gcc requirement; runtime compile approach is more portable |
| pybind11 / Cython | Overkill for single qfc function; ctypes is sufficient |
| pystatgen/limix | Heavy dependencies (BLAS, MKL), poor Windows support; the framework's own implementations suffice |
| polars | Not relevant to association statistics |

---

## Full pyproject.toml Delta

```toml
# ADDITIONS to pyproject.toml for association framework milestone

[project.optional-dependencies]
dev = [
    # ... existing unchanged ...
]
association = [                    # NEW extra
    "rpy2>=3.6.4",
    "matplotlib>=3.10",
]

[tool.hatch.build.targets.wheel]
packages = ["variantcentrifuge"]
artifacts = ["variantcentrifuge/data/*.c"]   # NEW: include qfc.c source

# CONSIDER: bump requires-python to ">=3.11" to allow scipy>=1.16
# Current: requires-python = ">=3.10"
# Impact: allows unpinned scipy on current versions
```

No changes to required `[project.dependencies]` — rpy2 and matplotlib are optional extras.

---

## Integration with Existing Stack

| New Capability | Uses Existing | Adds New |
|---------------|---------------|----------|
| SKAT-R backend | — | rpy2>=3.6.4 (optional) |
| SKAT pure-Python kernel | numpy, scipy.linalg.eigh | qfc.c via ctypes (runtime) |
| Davies p-value | ctypes stdlib | qfc.c C source (bundled) |
| Liu p-value fallback | scipy.stats.ncx2, norm | nothing |
| Logistic burden test | statsmodels Logit | nothing |
| Linear burden test | statsmodels OLS | nothing |
| Wald/LRT tests | statsmodels, scipy.stats.chi2 | nothing |
| QQ/Manhattan plots | — | matplotlib>=3.10 (optional) |
| ACAT-O combination | scipy.stats.cauchy | nothing |

---

## Cross-Platform Considerations

| Component | Linux | macOS | Windows | Notes |
|-----------|-------|-------|---------|-------|
| rpy2 | OK (Ubuntu 20+) | OK | Harder | R must be in PATH; Windows needs R on PATH explicitly |
| qfc.c runtime compile | gcc required | clang required | cl.exe (MSVC) required | Falls back to Liu if unavailable |
| matplotlib Agg backend | OK headless | OK | OK | No display server needed |
| scipy/statsmodels | OK | OK | OK | Pre-built wheels on all platforms |

**Windows qfc compilation:** MSVC (`cl.exe`) is used on Windows. Many Windows development environments (VS Code with C++ extension, Visual Studio) include it. HPC environments running Windows are rare; Linux is the dominant HPC OS. The Liu fallback is the practical Windows path.

**macOS:** clang is pre-installed via Xcode command line tools (`xcode-select --install`). qfc.c compiles successfully on macOS with `-shared -fPIC`.

---

## Confidence Assessment

| Area | Confidence | Source |
|------|------------|--------|
| rpy2 version/compatibility | HIGH | PyPI (3.6.4), GitHub issues verified |
| rpy2 R 4.4 issues | HIGH | GitHub issue #1107, confirmed scope (CentOS 7 only) |
| qfc.c function signature | HIGH | CompQuadForm R package source (rdrr.io) |
| ctypes vs cffi decision | HIGH | Technical analysis; ctypes stdlib advantage clear |
| scipy version constraint | HIGH | scipy release notes, PyPI metadata verified |
| Liu moment-matching | HIGH | Published algorithm (Liu 2009); scipy primitives verified |
| statsmodels API | HIGH | Official docs (0.14.4 stable, 0.15.0 dev), wald_test scalar note |
| matplotlib version | HIGH | PyPI (3.10.8 verified) |
| momentchi2 abandonment | HIGH | GitHub (last commit 2021, zero releases) |
| hatchling data files | MEDIUM | Official docs + hatch-cython PyPI; artifacts field behavior not tested in this project |

---

## Sources

- rpy2 PyPI (version 3.6.4, Python 3.8–3.13): https://pypi.org/project/rpy2/
- rpy2 changelog (R 4.3+ Rcomplex fix in 3.6.0): https://rpy2.github.io/doc/latest/html/changes.html
- rpy2 R 4.4 CentOS 7 issue: https://github.com/rpy2/rpy2/issues/1107
- rpy2 introduction (importr pattern, R initialization at import): https://rpy2.github.io/doc/latest/html/introduction.html
- statsmodels PyPI (0.14.6, December 2025): https://pypi.org/project/statsmodels/
- statsmodels wald_test API (scalar parameter note): https://www.statsmodels.org/stable/generated/statsmodels.discrete.discrete_model.LogitResults.wald_test.html
- scipy PyPI (1.17.0, January 2026): https://pypi.org/project/SciPy/
- scipy 1.15.0 release notes (last Python 3.10 version): https://docs.scipy.org/doc/scipy-1.16.1/release/1.15.0-notes.html
- scipy Python 3.10 drop tracking: https://github.com/scipy/scipy/issues/22881
- matplotlib PyPI (3.10.8, December 2025): https://pypi.org/project/matplotlib/
- momentchi2py GitHub (last commit August 2021): https://github.com/deanbodenham/momentchi2py
- CompQuadForm Davies R source (qfc signature): https://rdrr.io/cran/CompQuadForm/src/R/davies.R
- hatch-cython PyPI (0.6.0, C file support confirmed): https://pypi.org/project/hatch-cython/
- SKAT R package (p-value methods, Davies + Liu fallback): https://cran.r-project.org/web/packages/SKAT/vignettes/SKAT.pdf
- CFFI documentation (ABI vs API mode tradeoffs): https://cffi.readthedocs.io/en/latest/goals.html
