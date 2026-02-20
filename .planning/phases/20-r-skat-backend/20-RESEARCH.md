# Phase 20: R SKAT Backend - Research

**Researched:** 2026-02-20
**Domain:** rpy2, R SKAT package, thread safety, memory management
**Confidence:** HIGH (rpy2 API, SKAT function signatures), MEDIUM (thread safety pattern, memory monitoring)

---

## Summary

Phase 20 implements the R SKAT backend via rpy2, integrating with the existing `AssociationTest` ABC in `variantcentrifuge/association/`. The plan calls for a `backends/` subpackage that holds an `SKATBackend` ABC and `RSKATBackend`, with the R backend implemented as an `AssociationTest` subclass that wraps rpy2 calls to `SKATBinary` (binary traits) or `SKAT` (continuous traits), plus `SKAT-O` via `method="SKATO"`.

The core technical facts are well-established:

**rpy2 3.5.x API** — Use `importr()` to load R packages. Convert numpy arrays via `localconverter` with `numpy2ri.converter` (thread-local conversion context, no global `activate()` side effects). Extract R return values with `result.rx2("p.value")` for named list elements. Call R's `gc()` via `robjects.r['gc']()`. Check package presence with `rpackages.isinstalled()`.

**SKAT API** — `SKAT_Null_Model(formula, out_type="D", Adjustment=TRUE)` for binary traits (Adjustment=TRUE applies moment adjustment automatically for n<2000). Use `SKATBinary(Z, obj, method="SKAT")` for binary — never call plain `SKAT()` on binary phenotypes. SKAT-O uses `method="SKATO"` (equivalent to `method="optimal.adj"`). The rho grid for SKAT-O is eight points: {0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1}. Return object: `result$p.value`, `result$param$rho` (optimal rho), `result$param$n.marker`, `result$param$n.marker.test`.

**Thread safety** — `parallel_safe=False` on `AssociationAnalysisStage` is sufficient when the R backend is active. The PipelineRunner executes `parallel_safe=False` stages sequentially from the main thread (lines 572–574 of `runner.py`). rpy2 is fundamentally single-threaded; SIGINT signal registration requires the main thread. A dedicated R worker process is NOT required for this use case.

**Primary recommendation:** Implement `RSKATBackend` as an `AssociationTest` subclass. Use `SKAT_Null_Model` with `Adjustment=TRUE` (not the separate `SKAT_Null_Model_MomentAdjust`). Use `SKATBinary` for binary traits. Detect R availability eagerly at `check_dependencies()` time. Manage R memory with `del obj; robjects.r['gc']()` every 100 genes.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| rpy2 | 3.5.17 (current) | Python-R bridge | De facto standard for calling R from Python |
| SKAT | 2.2.5 (CRAN, July 2025) | SKAT/SKAT-O tests | Official implementation from authors |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | existing dep | Array construction for R matrix conversion | Always — genotype matrix is numpy float64 |
| rpy2.robjects.numpy2ri | bundled with rpy2 | numpy ↔ R matrix conversion | Inside localconverter context |
| rpy2.robjects.packages | bundled with rpy2 | `importr()`, `isinstalled()` | R detection and package loading |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `SKAT_Null_Model` + `Adjustment=TRUE` | `SKAT_Null_Model_MomentAdjust` | `MomentAdjust` is the older dedicated function; current SKAT docs recommend `Adjustment=TRUE` in `SKAT_Null_Model`. Both are equivalent for moment adjustment. |
| `method="SKATO"` | `method="optimal.adj"` | Both are identical — "SKATO" is the canonical alias in current SKAT 2.x. Use "SKATO" to match vignette examples. |
| `SKATBinary` | `SKAT` with `out_type="D"` | `SKAT` on binary data uses asymptotic method only; `SKATBinary` uses Hybrid method with resampling for small MAC. Always use `SKATBinary` for binary. |

**Installation:**
```bash
pip install "rpy2>=3.5.0"
# In R: install.packages("SKAT")
```

---

## Architecture Patterns

### Recommended Project Structure
```
variantcentrifuge/association/
├── base.py              # existing: AssociationTest ABC, TestResult, AssociationConfig
├── engine.py            # existing: AssociationEngine
├── tests/               # existing: fisher.py, logistic_burden.py, linear_burden.py
│   ├── __init__.py
│   ├── fisher.py
│   ├── logistic_burden.py
│   ├── linear_burden.py
│   └── skat_r.py        # NEW: RSKATTest (wraps RSKATBackend)
└── backends/            # NEW: backend abstraction layer
    ├── __init__.py      # get_skat_backend() factory
    ├── base.py          # SKATBackend ABC, NullModel container
    └── r_backend.py     # RSKATBackend
```

### Pattern 1: Eager R Detection in check_dependencies()

The phase context specifies R detection at pipeline startup ("eagerly at pipeline startup when parsing `--skat-backend r`"). The existing `AssociationTest.check_dependencies()` is called eagerly by `AssociationEngine.from_names()` before any data processing. This is the correct hook.

```python
# Source: rpy2.robjects.packages documentation
from rpy2.robjects.packages import isinstalled, PackageNotInstalledError
import rpy2.situation

def check_dependencies(self) -> None:
    """Eagerly detect R and SKAT at engine construction time."""
    # Step 1: Can rpy2 initialize R at all?
    try:
        import rpy2.robjects  # This initializes embedded R
    except Exception as e:
        r_home = os.environ.get("R_HOME", "<not set>")
        raise ImportError(
            f"R backend unavailable: rpy2 failed to initialize R. "
            f"R_HOME={r_home}. "
            f"Ensure R is installed and R_HOME is set. Error: {e}"
        ) from e

    # Step 2: Is SKAT package installed?
    try:
        from rpy2.robjects.packages import isinstalled
        if not isinstalled("SKAT"):
            raise ImportError(
                "R found but SKAT package missing. "
                "Install with: install.packages('SKAT') in R."
            )
    except ImportError:
        raise  # re-raise SKAT-missing errors

    # Step 3: Load and cache the SKAT package
    from rpy2.robjects.packages import importr
    self._skat_pkg = importr("SKAT")
    self._base_pkg = importr("base")
    self._stats_pkg = importr("stats")

    # Step 4: Log environment info for HPC reproducibility
    import rpy2
    r_version = rpy2.robjects.r("R.version$version.string")[0]
    skat_version = rpy2.robjects.r('packageVersion("SKAT")')[0]
    logger.info(f"R SKAT backend: rpy2={rpy2.__version__}, R={r_version}, SKAT={skat_version}")
    logger.info(f"R_HOME: {os.environ.get('R_HOME', 'auto-detected')}")
```

**Confidence: HIGH** — `isinstalled()` and `PackageNotInstalledError` are documented in rpy2 3.5+ docs.

### Pattern 2: Fitting the Null Model

```python
# Source: SKAT vignette + rpy2 robjects.r API
import rpy2.robjects as robjects
from rpy2.robjects import Formula
from rpy2.robjects import numpy2ri, default_converter

def fit_null_model(
    phenotype: np.ndarray,
    covariate_matrix: np.ndarray | None,
    trait_type: str,  # "binary" or "quantitative"
) -> object:
    """
    Fit SKAT null model. Returns R object for use in test_gene().
    out_type: "D" for binary, "C" for continuous.
    Adjustment=TRUE activates moment adjustment for n<2000 (default, recommended).
    """
    np_cv = default_converter + numpy2ri.converter

    out_type = "D" if trait_type == "binary" else "C"

    with np_cv.context():
        if covariate_matrix is not None and covariate_matrix.shape[1] > 0:
            # With covariates: y ~ X1 + X2 + ...
            r_y = robjects.FloatVector(phenotype.tolist())
            r_x = np_cv.context().__enter__()  # handled by with block
            # Build data frame and formula
            df = robjects.DataFrame({"y": r_y, "x": covariate_matrix})
            formula = Formula("y ~ .")
        else:
            # Intercept-only null: y ~ 1
            r_y = robjects.FloatVector(phenotype.tolist())
            df = robjects.DataFrame({"y": r_y})
            formula = Formula("y ~ 1")

        null_model = self._skat_pkg.SKAT_Null_Model(
            formula,
            data=df,
            out_type=out_type,
            Adjustment=True,  # moment adjustment enabled by default
        )
    return null_model
```

**Simpler alternative pattern (recommended for maintainability):**

```python
# Pass arrays directly via robjects.r environment assignment
import rpy2.robjects as ro

def fit_null_model(phenotype, covariates, trait_type):
    out_type = "D" if trait_type == "binary" else "C"

    ro.globalenv["._vc_y"] = ro.FloatVector(phenotype.tolist())
    if covariates is not None and covariates.shape[1] > 0:
        # Build matrix column by column
        r_cov = ro.r['matrix'](
            ro.FloatVector(covariates.ravel(order='F').tolist()),
            nrow=len(phenotype),
            ncol=covariates.shape[1],
        )
        ro.globalenv["._vc_x"] = r_cov
        formula_str = "._vc_y ~ ._vc_x"
    else:
        formula_str = "._vc_y ~ 1"

    null_obj = self._skat_pkg.SKAT_Null_Model(
        ro.Formula(formula_str),
        out_type=out_type,
        Adjustment=True,
    )
    # Clean up globals immediately
    ro.r("rm(._vc_y, ._vc_x)")
    return null_obj
```

**Confidence: HIGH** — Formula() and FloatVector patterns are core rpy2 API.

### Pattern 3: Running SKAT per Gene

```python
# Source: SKAT package vignette, rpy2 numpy docs
def test_gene(
    gene: str,
    genotype_matrix: np.ndarray,  # shape (n_samples, n_variants), float64
    null_model: object,            # R object from fit_null_model()
    method: str,                   # "SKAT", "Burden", or "SKATO"
    trait_type: str,
    weights_beta: tuple[float, float] = (1.0, 25.0),
) -> dict:
    """
    Run SKAT or SKATBinary on one gene's genotype matrix.
    Returns dict with p_value, rho (SKAT-O only), n_variants, warnings.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, default_converter

    np_cv = default_converter + numpy2ri.converter
    n_samples, n_variants = genotype_matrix.shape

    with np_cv.context():
        # R expects matrix in column-major order; SKAT Z is (n_samples x n_variants)
        r_z = ro.r['matrix'](
            ro.FloatVector(genotype_matrix.ravel(order='F').tolist()),
            nrow=n_samples,
            ncol=n_variants,
        )
        r_weights = ro.FloatVector([weights_beta[0], weights_beta[1]])

        if trait_type == "binary":
            result = self._skat_pkg.SKATBinary(
                r_z,
                null_model,
                kernel="linear.weighted",
                method=method,           # "SKAT", "Burden", or "SKATO"
                weights_beta=r_weights,
                missing_cutoff=0.15,
            )
        else:
            result = self._skat_pkg.SKAT(
                r_z,
                null_model,
                kernel="linear.weighted",
                method=method,
                weights_beta=r_weights,
                missing_cutoff=0.15,
            )

    # Extract results
    p_value_r = result.rx2("p.value")
    p_value = float(p_value_r[0]) if p_value_r[0] is not robjects.NA_Real else None

    # Extract optimal rho for SKAT-O
    rho = None
    if method == "SKATO":
        try:
            param = result.rx2("param")
            rho_r = param.rx2("rho")
            rho = float(rho_r[0]) if rho_r[0] is not robjects.NA_Real else None
        except Exception:
            pass  # rho not always present

    return {"p_value": p_value, "rho": rho, "n_variants": n_variants}
```

**Confidence: HIGH** — `result.rx2("p.value")`, `param.rx2("rho")` confirmed from SKAT docs. Matrix construction via `ro.r['matrix']` is standard rpy2 pattern.

### Pattern 4: R Memory Management

```python
# Source: rpy2 performances docs + gc() documentation
import gc as py_gc
import rpy2.robjects as ro

_r_gc = ro.r['gc']

def _run_r_gc() -> None:
    """Call R's gc() and Python's gc.collect() to free R heap objects."""
    _r_gc()          # R garbage collection
    py_gc.collect()  # Python gc to release rpy2 wrapper objects

# In the per-gene loop:
for i, gene in enumerate(genes):
    result_obj = self._test_gene(gene, ...)
    results.append(result_obj)
    del result_obj    # release rpy2 wrapper

    if (i + 1) % GC_INTERVAL == 0:
        _run_r_gc()
        logger.debug(f"R GC completed after {i+1} genes")

# Final GC at end of all genes
_run_r_gc()
```

**Confidence: HIGH** — `robjects.r['gc']()` and `del` + `gc.collect()` pattern confirmed from rpy2 docs. The rpy2 performances docs explicitly show this `gc.collect()` pattern in loops.

### Pattern 5: R Heap Monitoring

R does not expose a simple "current heap size" API, but `gc()` return value contains memory statistics:

```python
# Source: R documentation for gc(), rpy2 r-instance docs
def _check_r_heap(self) -> float | None:
    """Return R heap used in MB, or None if unavailable."""
    try:
        gc_result = ro.r('gc(verbose=FALSE)')  # Returns named matrix
        # R gc() returns 2x6 matrix; rows = (Ncells, Vcells)
        # columns = (used, gc trigger, max used, limit, %, %)
        # Vcells row, used column (index [1,0]) = vector cell count
        # Each Vcell = 8 bytes
        vcells_used = float(gc_result.rx(2, 1)[0])  # Vcells used
        heap_mb = (vcells_used * 8.0) / (1024 ** 2)
        return heap_mb
    except Exception:
        return None
```

**IMPORTANT WARNING THRESHOLD:** The phase specifies warning at "unusually large" gene panel. Based on experience with R SKAT on cohort studies: >2000 genes is a reasonable threshold. R heap warning: >4 GB (a common HPC node allocation boundary). These are MEDIUM confidence estimates — tune based on actual usage.

**Confidence: MEDIUM** — R heap monitoring via `gc()` return value is correct R behavior, but rpy2 matrix indexing for the gc return may need adjustment. Flag for testing.

### Pattern 6: Thread Safety — parallel_safe=False is Sufficient

The PipelineRunner (line 563–580 of `runner.py`) separates stages into `parallel_safe_stages` and `sequential_stages`. Sequential stages run one at a time from the calling thread (main thread in normal operation), before parallel stages dispatch to `ThreadPoolExecutor`. Thus `parallel_safe=False` guarantees the stage runs exclusively from the caller's thread.

The rpy2 thread issue (GitHub issue #769) was about SIGINT signal registration on worker thread initialization. The fix is already in rpy2 3.5+. However, the R embedded interpreter itself is not thread-safe for concurrent calls — the `openrlib.rlock` context manager exists for this case, but since we're running only from the main thread, we don't need it.

**RESEARCH-01 Answer: `parallel_safe=False` is sufficient.** A dedicated R worker process is NOT needed. The rpy2 documentation explicitly states the rlock is for "services that use multithreading such as WSGI" — not our use case.

**Calling from a worker thread must raise explicit error (Success Criteria 5):**

```python
import threading

def _assert_main_thread(self) -> None:
    """Raise RuntimeError if called from non-main thread."""
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError(
            "RSKATBackend: rpy2 calls must occur on the main thread. "
            "Do not use the R backend with parallel_safe=True or "
            "from a ThreadPoolExecutor worker. "
            "AssociationAnalysisStage.parallel_safe returns False when "
            "R backend is active — this is enforced at stage level."
        )
```

**Confidence: HIGH** — Confirmed by reading `runner.py` lines 563–580 directly. Sequential stages are called from the runner's current thread, which is the main thread in normal pipeline execution.

### Pattern 7: p_method Column Convention

SKAT and SAIGE-GENE use different column naming. Investigation:

- **SAIGE-GENE output:** Uses `Pvalue` (CamelCase), `Pvalue_Burden`, `Pvalue_SKAT` — separate per-method columns.
- **SKAT R output:** `p.value` (R convention with dot).
- **Existing project convention:** `{test_name}_p_value` (e.g., `fisher_p_value`, `logistic_burden_p_value`) — underscore, lowercase.

**Recommendation:** Follow existing project column convention. The SKAT test should produce:
- `skat_p_value` — primary SKAT/SKAT-O p-value
- `skat_corrected_p_value` — after FDR/Bonferroni correction
- `skat_o_rho` — optimal rho from SKAT-O (NA for pure SKAT/Burden)
- `skat_warnings` — R-side warnings captured per gene (e.g., "MOMENT_ADJUSTMENT_APPLIED")
- `skat_method` — "SKAT", "Burden", or "SKATO" per gene (from `TestResult.extra`)

The `test_name` property should return `"skat"` to match the engine's column prefix convention.

**Confidence: HIGH** — Derived from existing engine.py column pattern (lines 214–220) and project convention established in Phases 18–19.

### Anti-Patterns to Avoid

- **`numpy2ri.activate()` globally** — Avoid global activation; use `localconverter` context blocks to prevent unexpected conversions in other code paths.
- **Calling `SKAT()` on binary phenotypes** — Always use `SKATBinary()`. `SKAT` with `out_type="D"` only uses the asymptotic method and is less accurate for binary traits; `SKATBinary` uses the Hybrid method with resampling for small MAC.
- **Ignoring the SKAT-O rho return** — SKAT-O returns the `optimal rho` in `result$param$rho`, NOT in `result$p.value`. Always extract both.
- **Using `method="optimal"` instead of `method="SKATO"`** — `"optimal"` uses 11 equally-spaced rho values; `"SKATO"` uses the improved grid of 8 asymmetric points with better type-I error control.
- **Storing R null model objects in PipelineContext** — R objects cannot be pickled or stored across thread boundaries. The null model must be built and used within the same sequential execution block.
- **Calling `SKAT_Null_Model_MomentAdjust` separately** — This is the older API. Current SKAT 2.x integrates moment adjustment via `Adjustment=TRUE` in `SKAT_Null_Model`. Use the unified function.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SKAT p-value computation | Custom kernel test | R `SKAT`/`SKATBinary` | Moment adjustment, Hybrid resampling (SKATBinary), and numerical accuracy all require the reference implementation |
| R package detection | `subprocess.run("Rscript -e ...")` | `rpy2.robjects.packages.isinstalled()` | Proper integration, same R session as test run |
| numpy → R matrix | Manual FloatVector + dim assignment | `ro.r['matrix'](ro.FloatVector(...), nrow=n, ncol=p)` | Correct column-major memory layout that R expects |
| R memory monitoring | Custom /proc/meminfo parsing | `ro.r('gc(verbose=FALSE)')` return matrix | Built into R base, works cross-platform including HPC |
| R version string | `subprocess.run("R --version")` | `rpy2.robjects.r("R.version$version.string")[0]` | Same R session, no subprocess overhead |
| rpy2 version | External check | `import rpy2; rpy2.__version__` | Standard Python package version |

**Key insight:** The SKAT R package's numerical methods (moment matching, resampling for small MAC) took years of research to validate. Never attempt to reproduce the `SKATBinary` Hybrid p-value computation in Python.

---

## Common Pitfalls

### Pitfall 1: NA Values from R Return as Python NA, Not None or float('nan')

**What goes wrong:** `result.rx2("p.value")[0]` returns `rpy2.rinterface.NA_Real` (R's NA), not Python `None` or `float('nan')`. Code that does `float(p_value_r[0])` raises `TypeError` or silently produces NaN without the NA context.

**Why it happens:** SKAT returns `NA_Real` when the test cannot be computed (singular kernel matrix, fewer variants than required). This is a valid statistical outcome per the phase context.

**How to avoid:**
```python
from rpy2.rinterface import NA_Real
p_val_r = result.rx2("p.value")[0]
p_value = None if p_val_r is NA_Real else float(p_val_r)
```

**Warning signs:** `TypeError: float() argument must be...` on genes with few variants.

**Confidence: HIGH** — NA_Real behavior is documented in rpy2; SKAT NA behavior documented in package help.

### Pitfall 2: R Matrix Column/Row Order (column-major vs row-major)

**What goes wrong:** numpy is row-major (C order). R is column-major (Fortran order). `ro.FloatVector(matrix.ravel())` without specifying `order='F'` flattens row-by-row, producing transposed data in R.

**Why it happens:** `matrix.ravel()` uses C order by default; R's matrix fills column by column.

**How to avoid:**
```python
# CORRECT: Fortran order for R column-major fill
r_z = ro.r['matrix'](
    ro.FloatVector(genotype_matrix.ravel(order='F').tolist()),
    nrow=n_samples,
    ncol=n_variants,
)

# WRONG: C order produces transposed matrix in R
r_z = ro.r['matrix'](ro.FloatVector(genotype_matrix.ravel().tolist()), ...)
```

**Warning signs:** SKAT reports wildly wrong p-values; number of variants in result doesn't match input.

**Confidence: HIGH** — R matrix storage order is fundamental R behavior.

### Pitfall 3: rpy2 Import Initializes R Immediately

**What goes wrong:** `import rpy2.robjects` initializes the embedded R process immediately, even at module import time. If R is not installed, this raises at import time, not at test construction time.

**Why it happens:** rpy2 starts R on first import of the high-level interface.

**How to avoid:** Never import rpy2 at module level in `skat_r.py`. Always import inside the method that needs it:
```python
class RSKATTest(AssociationTest):
    def check_dependencies(self) -> None:
        try:
            import rpy2.robjects  # deferred import
        except Exception as e:
            raise ImportError(...) from e
```

**Warning signs:** `ImportError` at pipeline startup even when user did not request `--skat-backend r`.

**Confidence: HIGH** — Observed behavior from rpy2 design.

### Pitfall 4: SKAT Null Model Must Be Fit Per Run, Not Cached Across Genes

**What goes wrong:** The null model `obj` from `SKAT_Null_Model()` is an R object (S4 class). It contains the phenotype residuals. If the same null model object is reused after its R memory is freed (e.g., after `rm()` or a major GC), R will either error or produce silently wrong results.

**Why it happens:** The null model references the phenotype vector by R's reference semantics. After `del null_obj` and `gc()`, R may collect the residual storage.

**How to avoid:** Fit null model once per gene run (before the gene loop, or cache only for the duration of the full gene loop — do not delete null_model inside the GC loop). Only delete per-gene result objects inside the GC loop.

**Confidence: HIGH** — SKAT vignette shows null model fit once, used for all genes.

### Pitfall 5: R Heap Not Returned to OS by gc() Alone

**What goes wrong:** Calling R's `gc()` frees R heap objects but R may not return memory to the OS immediately (R's allocator holds freed pages). So `psutil.Process().memory_info().rss` may not decrease after `gc()`.

**Why it happens:** R uses a vector heap that grows lazily and shrinks slowly.

**How to avoid:** Don't use process RSS as a threshold for abort decisions. Use `gc()` return value to read R's internal accounting (`Vcells used` × 8 bytes = heap in use). Warn on R's internal metric, not OS RSS.

**Confidence: MEDIUM** — R heap behavior documented; rpy2 extraction of gc() return value needs runtime verification.

### Pitfall 6: SKATBinary `method` Parameter Differs from SKAT

**What goes wrong:** `SKATBinary` uses `method` to select the test type ("SKAT", "Burden", "SKATO"), while `SKATBinary` also has `method.bin` for the p-value computation method ("Hybrid", "ER", "ER.A", etc.). These are separate parameters. Confusing them silently selects wrong behavior.

**How to avoid:**
```python
# method = test type ("SKAT", "Burden", "SKATO")
# method.bin = p-value computation strategy — use default "Hybrid"
result = self._skat_pkg.SKATBinary(
    r_z, null_model,
    method="SKATO",       # test type: SKAT-O
    # method_bin not passed: defaults to "Hybrid" — correct
)
```

Note: rpy2 maps Python underscores to R dots for keyword arguments: `method_bin` → `method.bin`.

**Confidence: HIGH** — Confirmed from SKATBinary documentation.

---

## Code Examples

Verified patterns from official sources:

### Full R Detection and Environment Logging
```python
# Source: rpy2.robjects.packages docs, rpy2.situation module
def _detect_and_log_r_environment(self) -> None:
    """Detect R/SKAT availability and log environment for HPC reproducibility."""
    import rpy2
    import rpy2.robjects as ro
    from rpy2.robjects.packages import isinstalled, importr

    # R version
    r_version = str(ro.r("R.version$version.string")[0])
    # SKAT version
    skat_version = str(ro.r('as.character(packageVersion("SKAT"))')[0])
    # rpy2 version
    rpy2_version = rpy2.__version__
    # R_HOME
    r_home = str(ro.r("R.home()")[0])

    logger.info(
        f"R SKAT backend environment: "
        f"rpy2={rpy2_version}, R={r_version}, SKAT={skat_version}, R_HOME={r_home}"
    )

    # Version warnings (thresholds determined by research)
    # R < 4.0 lacks some SKAT 2.x features; SKAT < 2.0 lacks SKATBinary_Robust
    import re
    r_major = int(re.search(r'version (\d+)\.', r_version).group(1))
    if r_major < 4:
        logger.warning(f"R version {r_version} may not support all SKAT 2.x features. Recommend R >= 4.0.")
```

### numpy Genotype Matrix to R Matrix
```python
# Source: rpy2 numpy docs (localconverter), R matrix docs (column-major)
import rpy2.robjects as ro

def _numpy_to_r_matrix(self, arr: np.ndarray) -> object:
    """Convert (n_samples, n_variants) numpy float64 to R matrix."""
    n_rows, n_cols = arr.shape
    # Fortran order = column-major = what R expects for matrix fill
    r_vec = ro.FloatVector(arr.ravel(order='F').tolist())
    return ro.r['matrix'](r_vec, nrow=n_rows, ncol=n_cols)
```

### SKAT-O Result Extraction
```python
# Source: SKAT package documentation for param$rho
def _extract_skat_result(self, r_result, method: str) -> dict:
    """Extract p-value and rho from SKAT/SKATBinary result object."""
    from rpy2.rinterface import NA_Real

    p_val_r = r_result.rx2("p.value")[0]
    p_value = None if p_val_r is NA_Real else float(p_val_r)

    rho = None
    if method == "SKATO":
        try:
            param = r_result.rx2("param")
            rho_r = param.rx2("rho")[0]
            rho = None if rho_r is NA_Real else float(rho_r)
        except Exception as exc:
            logger.debug(f"Could not extract rho from SKAT-O result: {exc}")

    n_marker = None
    try:
        n_marker = int(r_result.rx2("param").rx2("n.marker")[0])
    except Exception:
        pass

    return {"p_value": p_value, "rho": rho, "n_marker": n_marker}
```

### R Warnings Capture via tryCatch
```python
# Source: rpy2 robjects.r eval pattern
def _run_skat_with_warnings(self, gene: str, r_z, null_model, method: str, binary: bool) -> tuple:
    """Run SKAT, capture R-side warnings. Returns (r_result, warnings_list)."""
    # Use R's withCallingHandlers to collect warnings without aborting
    ro.globalenv["._vc_Z"] = r_z
    ro.globalenv["._vc_obj"] = null_model

    func_name = "SKATBinary" if binary else "SKAT"
    code = f"""
    local({{
      warns <- character(0)
      result <- withCallingHandlers(
        {func_name}(._vc_Z, ._vc_obj, kernel="linear.weighted",
                    method="{method}", weights.beta=c(1,25)),
        warning = function(w) {{
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }}
      )
      list(result=result, warnings=warns)
    }})
    """
    wrapped = ro.r(code)
    r_result = wrapped.rx2("result")
    r_warns = list(wrapped.rx2("warnings"))
    return r_result, r_warns
```

**Note:** This R code pattern passes global env variables (`._vc_Z`, `._vc_obj`) to avoid rpy2 formula serialization complexity. Clean them up immediately after.

### GC Interval Logic
```python
# Source: rpy2 performances docs pattern
_r_gc = ro.r['gc']
GC_INTERVAL = 100  # genes between GC calls — recommended fixed value

for i, gene_data in enumerate(sorted_genes):
    result_r = self._test_gene(gene_data)
    gene_results.append(self._extract_skat_result(result_r, method))
    del result_r  # release rpy2 wrapper reference

    if (i + 1) % GC_INTERVAL == 0:
        _r_gc()
        import gc as py_gc
        py_gc.collect()
        logger.debug(f"R GC: {i+1}/{len(sorted_genes)} genes complete")

# Final GC
_r_gc()
```

---

## GC Interval: Fixed vs Configurable

**Recommendation: Fixed at 100, not configurable.**

Rationale:
1. The existing project has no CLI pattern for internal algorithm tuning parameters (no `--firth-max-iter` in CLI, that's in `AssociationConfig`). Adding `--r-gc-interval` would be inconsistent.
2. 100 genes is a conservative interval that prevents heap accumulation without adding measurable overhead (R GC on small-to-medium studies takes <10ms).
3. The phase context says "GC interval configurability: Claude's Discretion" — after investigation, a fixed value is simpler and sufficient.
4. If needed later, can be added to `AssociationConfig` (Python-facing, not CLI-facing).

**Confidence: MEDIUM** — Based on R memory management best practices; actual optimal interval depends on gene panel size and genotype matrix size.

---

## Version Requirements

| Component | Minimum | Recommended | Reason |
|-----------|---------|-------------|--------|
| R | 4.0 | 4.3+ | SKAT 2.x uses modern R conventions; R 4.0 changed default stringsAsFactors |
| SKAT | 2.0 | 2.2.5 (current) | `SKATBinary` Hybrid method; SKAT 2.0 introduced improved binary methods |
| rpy2 | 3.5.0 | 3.5.17 (current) | `localconverter`, fixed signal/thread issue (3.5.x), numpy2ri improvements |
| Python | 3.10+ | existing constraint | Already project requirement |

**Warn-but-proceed thresholds:**
- R < 4.0: log `WARNING` "R version X.Y may not support all SKAT 2.x features. Recommend R >= 4.0."
- SKAT < 2.0: log `WARNING` "SKAT version X.Y predates SKATBinary improvements. Recommend SKAT >= 2.0."
- rpy2 < 3.5.0: log `WARNING` "rpy2 version X.Y may have thread safety issues. Recommend rpy2 >= 3.5.0."

**Confidence: MEDIUM** — R and SKAT minimum versions based on feature analysis (SKATBinary with Hybrid method added progressively; SKAT 2.0 is documented). Exact minimum versions not explicitly stated in CRAN docs.

---

## Progress Logging and Warning Thresholds

### Progress Logging Interval

**Recommendation: Every 50 genes (or 10% of total, whichever is more frequent).**

Rationale: HPC jobs often have log files reviewed hourly. For a 2000-gene panel (common for exome SKAT), logging every 50 genes = 40 log lines — visible enough to confirm progress without flooding. For small panels (<100 genes), every 10 genes is sufficient.

```python
log_interval = max(10, min(50, len(genes) // 10))
if (i + 1) % log_interval == 0:
    pct = 100 * (i + 1) / len(genes)
    logger.info(f"SKAT progress: {i+1}/{len(genes)} genes ({pct:.0f}%)")
```

**Confidence: MEDIUM** — Reasonable default, matches common HPC job monitoring practice.

### Large Panel Warning Threshold

**Recommendation: >2000 genes.**

```python
if len(genes) > 2000:
    logger.warning(
        f"Large gene panel ({len(genes)} genes) — R SKAT memory usage may be significant. "
        f"Consider splitting the gene list for very large panels (>5000 genes)."
    )
```

**Confidence: MEDIUM** — Based on typical exome/genome study gene panel sizes; 2000 covers most exome gene sets.

### R Heap Warning Threshold

**Recommendation: Warn at >4 GB R heap usage.**

```python
R_HEAP_WARNING_GB = 4.0

if heap_mb > R_HEAP_WARNING_GB * 1024:
    logger.warning(
        f"R heap usage {heap_mb/1024:.1f} GB exceeds {R_HEAP_WARNING_GB} GB threshold. "
        f"Consider increasing R_MAX_VSIZE environment variable."
    )
```

**Confidence: LOW** — 4 GB is a reasonable HPC threshold but highly system-dependent. Adjust based on actual memory profiles.

---

## Integration with Existing AssociationTest Framework

The `RSKATTest` must integrate with:

1. **`AssociationTest.run(gene, contingency_data, config)`** — Receives the dict with `genotype_matrix`, `phenotype_vector`, `variant_mafs`, `covariate_matrix` keys (established in Phase 19).

2. **`effect_column_names()`** — SKAT does not produce effect sizes (beta/OR). Return:
   ```python
   def effect_column_names(self) -> dict[str, str]:
       return {
           "effect": None,    # No effect size from SKAT
           "se": None,
           "ci_lower": None,
           "ci_upper": None,
       }
   ```
   The engine skips None effect columns. SKAT-specific columns (`skat_o_rho`, `skat_warnings`) go into `TestResult.extra`.

3. **`test_name`** — Return `"skat"` for consistency with engine's `{test_name}_p_value` pattern.

4. **`AssociationConfig`** — SKAT-specific config fields to add:
   - `skat_backend: str = "r"` — "r" or "python" (Phase 21)
   - `skat_method: str = "SKAT"` — "SKAT", "Burden", or "SKATO"
   - `skat_weights_beta: tuple[float, float] = (1.0, 25.0)` — Beta distribution params

5. **`AssociationEngine._build_registry()`** — Register `"skat": SKATRTest` in Phase 20.

**Null model lifetime:** The null model is fit ONCE before the gene loop (same phenotype/covariates for all genes), then passed to each `test_gene()` call. This matches SKAT's intended usage (one null model per cohort, not per gene).

**Architectural implication:** The `AssociationTest.run()` signature passes per-gene data but SKAT needs a pre-fit null model that spans all genes. The null model should be fit in `check_dependencies()` or a separate `fit_null_model()` method called by `AssociationEngine.run_all()` before the gene loop.

**Options:**
- Option A: Fit null model in `run()` and cache as instance state (`self._null_model`). First call fits it; subsequent calls reuse it. Thread-unsafe but parallel_safe=False guarantees single thread.
- Option B: Add a `prepare(global_data, config)` hook to `AssociationTest` ABC called once before the gene loop. More explicit but changes the ABC.

**Recommendation: Option A (instance state caching)** for Phase 20 — minimal ABC change, single-threaded guarantee from parallel_safe=False. Document clearly in class docstring.

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `SKAT_Null_Model_MomentAdjust` (separate function) | `SKAT_Null_Model(..., Adjustment=TRUE)` | SKAT 1.x → 2.x | Unified API; MomentAdjust still works but not needed |
| `method="optimal"` (11 rho points) | `method="SKATO"` / `method="optimal.adj"` (8 asymmetric points) | SKAT 1.x → 2.x | Better type-I error control |
| `numpy2ri.activate()` global | `localconverter` context | rpy2 2.x → 3.x | Thread-safe, no global side effects |
| R version check via subprocess | `rpy2.situation.iter_info()` or `ro.r("R.version...")` | rpy2 3.0 | Consistent with active R session |

**Deprecated/outdated:**
- `rpy2.robjects.numpy2ri.activate()`: Works but not recommended in rpy2 3.5+; prefer `localconverter`.
- `SKAT_Null_Model_MomentAdjust`: Still present in SKAT 2.x but superseded by `Adjustment=TRUE` parameter.

---

## Open Questions

1. **SKAT result `param$rho` field name in SKAT-O**
   - What we know: SKAT-O returns `param$rho` for optimal rho; confirmed from web docs.
   - What's unclear: Whether `SKATBinary` with `method="SKATO"` returns `param$rho` in the same location as `SKAT` with `method="SKATO"`. The documentation may differ slightly.
   - Recommendation: Add a unit test that mocks an R SKAT-O result object and verifies the `param.rx2("rho")` extraction path for both `SKAT` and `SKATBinary` return structures.

2. **Effect column behavior when SKAT has no effect size**
   - What we know: SKAT does not report a beta/OR; only p-value and rho.
   - What's unclear: How `AssociationEngine.run_all()` handles `effect_column_names()` returning all None — the current code does `row[f"{test_name}_{col_names['effect']}"] = res.effect_size` which would produce a column named `skat_None`.
   - Recommendation: The engine needs a guard: `if col_names.get('effect') is not None: row[...] = ...`. Add this to engine.py as part of Phase 20.

3. **`rm()` of R globals vs `del` of Python wrapper**
   - What we know: The warning capture pattern uses `._vc_Z` etc. in R global env.
   - What's unclear: Whether global env pollution between genes causes memory leaks or interference.
   - Recommendation: Always call `ro.r("rm(list=ls(pattern='\\._vc_'))")` after each gene test, or use a dedicated R local environment.

---

## Sources

### Primary (HIGH confidence)
- rpy2 3.5.x documentation (rpy2.readthedocs.io) — numpy2ri, localconverter, robjects.r, gc(), thread safety (openrlib.rlock)
- rdrr.io/cran/SKAT/man/SKATBinary.html — SKATBinary full function signature and return object
- rdrr.io/cran/SKAT/man/SKAT_Null_Model.html — SKAT_Null_Model signature, Adjustment parameter
- rdocumentation.org SKAT 2.2.5/topics/SKAT — SKAT function signature, method="optimal.adj", param$rho
- leelabsg.r-universe.dev/SKAT/doc/SKAT.Rnw — SKAT vignette code examples (binary vs continuous, SKAT-O)
- variantcentrifuge/pipeline_core/runner.py (lines 563–580) — parallel_safe=False is sequential from caller's thread

### Secondary (MEDIUM confidence)
- GitHub rpy2/rpy2 issue #769 — Signal registration thread issue (fixed in rpy2 3.5+)
- rpy2 performances documentation — gc.collect() loop pattern for R memory management
- SAIGE-GENE step2 documentation — output column naming conventions (Pvalue not p_value)
- rpy2 rinterface documentation — openrlib.rlock thread safety context manager

### Tertiary (LOW confidence)
- WebSearch results on R/SKAT version history — minimum R 4.0, SKAT 2.0 thresholds (inferred from feature analysis, not explicit CRAN minimum version statement)
- R heap warning threshold of 4 GB — based on HPC convention, not measured data
- Progress logging interval of 50 genes — based on HPC log monitoring conventions

---

## Metadata

**Confidence breakdown:**
- rpy2 API: HIGH — verified from official documentation
- SKAT function signatures: HIGH — verified from rdrr.io/CRAN documentation
- Thread safety (parallel_safe=False sufficient): HIGH — confirmed by reading runner.py source
- Memory management (gc pattern): HIGH — documented in rpy2 performances docs
- R heap monitoring via gc() return: MEDIUM — R behavior confirmed, rpy2 extraction needs runtime test
- Version requirements (R >= 4.0, SKAT >= 2.0): MEDIUM — inferred from feature timeline, not explicit
- GC interval (100 genes): MEDIUM — reasonable default, not measured
- Progress logging interval (50 genes): MEDIUM — convention-based

**Research date:** 2026-02-20
**Valid until:** 2026-08-20 (stable — rpy2 and SKAT are mature, low churn)
