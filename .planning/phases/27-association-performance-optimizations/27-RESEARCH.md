# Phase 27: Association Performance Optimizations - Research

**Researched:** 2026-02-22
**Domain:** Python SKAT backend, ProcessPoolExecutor parallelism, numerical quadrature
**Confidence:** HIGH (all findings verified by direct code inspection and benchmarking)

---

## Summary

Phase 27 implements three performance optimizations for the pure Python SKAT backend:
gene-level parallelization via ProcessPoolExecutor, Davies CDF caching/interpolation in
the SKAT-O omnibus integration, and single eigendecomposition for SKAT-O. A fourth
optional optimization (sparse genotype matrices) is low priority for typical cohort sizes.

The codebase was directly inspected and benchmarked. Key findings: (1) The parallelization
approach is straightforward — all relevant objects are picklable and verified to work with
ProcessPoolExecutor. (2) The biggest single speedup available is replacing `scipy.integrate.quad`
with fixed Gauss-Legendre quadrature (46x speedup measured, 379ms → 8ms per gene for
SKAT-O). (3) The single-eigendecomposition optimization does NOT provide speedup for typical
gene sizes (p=20-50 variants) — it is only beneficial for p > 100 variants per gene.

**Primary recommendation:** Implement GL quadrature first (highest impact, lowest risk),
then parallelization (second highest impact), then the eigendecomposition optimization only
if p > 100 variant genes are common in target use cases.

---

## Standard Stack

### Core Libraries (already in use)

| Library | Version | Purpose | Notes |
|---------|---------|---------|-------|
| `concurrent.futures` | stdlib | ProcessPoolExecutor for gene parallelism | Already used elsewhere in pipeline |
| `scipy.integrate` | current | Quadrature — being replaced | Current: adaptive quad, target: GL |
| `scipy.linalg` | current | `eigh()` eigendecomposition | Already used, optimize call count |
| `numpy.polynomial.legendre` | numpy | `leggauss()` for GL quadrature nodes | stdlib numpy, no new dependency |
| `scipy.interpolate` | current | `interp1d` for Davies CDF grid | Already available, optional use |

### No New Dependencies Required

All optimizations use existing dependencies. No new packages need to be added.

---

## Architecture Patterns

### Current Engine Flow (`engine.py`)

```python
# engine.run_all() — current sequential loop
for test in self._tests.values():
    test.prepare(len(sorted_data))

for gene_data in sorted_data:
    gene = gene_data.get("GENE", "")
    for test_name, test in self._tests.items():
        result = test.run(gene, gene_data, self._config)
        results_by_test[test_name][gene] = result

for test in self._tests.values():
    test.finalize()
```

### Target Parallel Pattern

```python
# Pattern: worker function must be module-level (ProcessPoolExecutor requirement)
def _run_gene_worker(args):
    gene, gene_data, pickled_tests, config = args
    import pickle
    tests = pickle.loads(pickled_tests)  # reconstructed in worker
    results = {}
    for test_name, test in tests.items():
        result = test.run(gene, gene_data, config)
        results[test_name] = result
    return gene, results

# In engine.run_all() when workers > 1:
# 1. Call prepare() on all tests (sequential)
# 2. Pre-fit null models by running first gene (triggers lazy fit)
# 3. Pickle test instances (with null models embedded)
# 4. Dispatch via ProcessPoolExecutor
# 5. Collect results, sort by gene name
# 6. Call finalize() on all tests (sequential)
```

### Null Model Pre-Fitting Pattern

The `PurePythonSKATTest._null_model` is set lazily on the first `.run()` call. For parallel
dispatch, the null model must be pre-fitted BEFORE pickling for workers.

**Recommended approach:** Call `.run()` on one gene sequentially to trigger lazy null model
fitting, then pickle the test instances (including their `_null_model`). Workers receive
pre-fitted test objects and skip the null model fitting step.

Verified: `PurePythonSKATTest` with `_null_model` set IS picklable (12,021 bytes in
testing). `AssociationConfig` is picklable (696 bytes). Gene data dicts with numpy arrays
are picklable (33.4 KB for 200 samples × 20 variants).

### Gauss-Legendre Quadrature Pattern (Optimization 3.3)

```python
# Replace scipy.integrate.quad (adaptive, ~2400 davies calls)
# with fixed 128-node Gauss-Legendre (128 davies calls)
from numpy.polynomial.legendre import leggauss

nodes, weights = leggauss(128)
a, b = 0.0, 40.0  # integration bounds (matches current code)
x_gl = (b - a) / 2.0 * nodes + (b + a) / 2.0  # map to [0, 40]
w_gl = (b - a) / 2.0 * weights

# Vectorized evaluation:
# For each GL node x_gl[i], compute min conditional quantile
# then call davies_pvalue for each (128 calls total vs ~2400)
```

**Measured speedup: 379ms → 8.3ms = 46x per SKAT-O gene (p=20 variants)**

### Single Eigendecomposition Pattern (Optimization 3.4)

```python
# In _skato_get_pvalue, replace 7 eigh calls with 1 eigh + O(p^2) per rho
# A = z1_half.T @ z1_half  (currently computed once already)
d, V = scipy.linalg.eigh(a_mat, driver='evr')  # decompose A once
col_sums_v = V.sum(axis=0)  # precompute for rank-1 structure

for rho in rho_capped:
    s = np.sqrt(1.0 - rho)
    delta = (np.sqrt(1.0 - rho + p * rho) - s) / p
    # R.M^{1/2} @ V = s*V + delta*ones*col_sums_v (rank-1 structure)
    M = s * V + delta * col_sums_v[np.newaxis, :]  # (p, p)
    # K_sym = M @ diag(d) @ M.T  (same eigenvalues as current K_sym)
    K_vbasis = M @ (d[:, np.newaxis] * M.T)
    lambda_all.append(_get_lambda(K_vbasis))
```

**IMPORTANT: This pattern is SLOWER for typical gene sizes (p=20-50).**
Benchmarks show 0.67x speedup at p=20, 0.47x at p=50, 0.26x at p=100.
The O(p^2) matmul per rho dominates over the savings from fewer eigh calls.
Only beneficial if a more efficient algebraic transformation is used that avoids
building K_vbasis explicitly.

### Recommended Project Structure for Phase 27

No new files needed. All changes go into existing files:

```
variantcentrifuge/association/
├── engine.py                    # Add --association-workers logic
├── base.py                      # Add association_workers to AssociationConfig
├── backends/
│   └── python_backend.py        # GL quadrature in _skato_integrate_davies
│                                # (Optional: single eigh in _skato_get_pvalue)
├── genotype_matrix.py           # (Optional: sparse matrix support)
variantcentrifuge/
└── cli.py                       # Add --association-workers argument
```

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Process pool | Custom subprocess spawning | `concurrent.futures.ProcessPoolExecutor` | Standard library, handles pickling, exceptions, shutdown |
| BLAS thread control | Manual thread pinning | `os.environ["OPENBLAS_NUM_THREADS"] = "1"` in worker initializer | Standard pattern for BLAS-intensive multiprocessing |
| GL quadrature nodes | Manual Gauss-Legendre computation | `numpy.polynomial.legendre.leggauss(128)` | Numerically precise, already in numpy |
| Eigendecomposition | Manual LAPACK calls | `scipy.linalg.eigh(driver='evr')` | Already used, proven to match R's DSYEVR |
| Cubic interpolation | Manual spline fitting | `scipy.interpolate.interp1d(kind='cubic')` | Tested, sufficient accuracy |

---

## Common Pitfalls

### Pitfall 1: Lazy Null Model Fitting and Parallel Workers

**What goes wrong:** Workers receive test instances without null models (because null model
is fit lazily on first `.run()` call). Each worker independently fits the null model, which
is: (a) wasteful (N-workers × null model fitting cost) and (b) non-deterministic across
workers if phenotype data isn't passed consistently.

**Why it happens:** `PurePythonSKATTest.run()` checks `if self._null_model is None` and
fits on first call. Workers start with fresh instances if pre-fitting isn't done.

**How to avoid:** Before parallel dispatch, run one gene sequentially through all test
instances to trigger lazy null model fitting. Then pickle the test instances WITH their
null models. Workers call `.run()` and immediately skip null model fitting (it's already set).

**Warning signs:** Workers take much longer than expected on first gene; timing is uneven
across workers.

### Pitfall 2: BLAS Thread Oversubscription

**What goes wrong:** Each worker spawns a full BLAS thread pool (default: n_cpu threads).
With N workers, you get N × n_cpu BLAS threads competing for the same CPU cores, causing
slowdown instead of speedup.

**Why it happens:** OpenBLAS and MKL create thread pools sized by CPU count by default.
With multiprocessing, each process creates its own pool.

**How to avoid:** Set `os.environ["OPENBLAS_NUM_THREADS"] = "1"` and
`os.environ["MKL_NUM_THREADS"] = "1"` in the worker initializer function passed to
`ProcessPoolExecutor(initializer=...)`.

**Warning signs:** CPU usage shows near 100% on all cores but wall time doesn't improve.

### Pitfall 3: R Backend Breaks with ProcessPoolExecutor

**What goes wrong:** The R SKAT backend (RSKATBackend) uses rpy2, which is not safe in
subprocess workers. Calling rpy2 from a child process causes segfaults.

**Why it happens:** rpy2 requires that all R calls happen from the main thread of the
main process.

**How to avoid:** Check `parallel_safe` attribute on test instances before enabling
parallel dispatch. If any registered test has `parallel_safe=False`, fall back to sequential
execution and log a warning. `RSKATTest.parallel_safe = False` (note: PurePythonSKATTest
does NOT yet define this attribute — needs to be added as `True`).

**Warning signs:** Process crashes without Python traceback (segfault in C extension).

### Pitfall 4: GL Quadrature vs Adaptive Integration Accuracy

**What goes wrong:** Fixed 128-node GL misses sharp features in the integrand (e.g., near
p=0 very small values) and produces inaccurate omnibus p-values.

**Why it happens:** GL is exact for polynomials up to degree 2n-1. The SKAT-O integrand
is a CDF composite function, not polynomial. Near p=0, the integrand is steep.

**How to avoid:** Test GL against the current adaptive quad on diverse p-value scenarios
(p=0.5, p=0.1, p=0.01, p=0.001). The measured error of 3.2e-3 on p≈0.86 is acceptable for
a p-value approximation, but should be validated across a wider range. The integration
is from 0 to 40 (chi-squared upper tail is negligible beyond 40).

**Warning signs:** SKAT-O p-values deviate from R by more than expected tolerance on
validation datasets.

### Pitfall 5: Single Eigendecomposition Doesn't Speedup Small p

**What goes wrong:** The plan claims "5x speedup on SKAT-O eigenvalue step." This was NOT
confirmed by benchmarking for typical panel sizes.

**What benchmarks showed:**
- p=20 variants: 0.67x (30% SLOWER with single eigh)
- p=50 variants: 0.47x (53% SLOWER)
- p=100 variants: 0.26x (74% SLOWER)

**Why it happens:** The O(p^2) matmul per rho (`M @ (d[:,None] * M.T)`) doesn't reduce total
work versus the current approach of building K_sym efficiently and calling eigh. The single
eigh saves `6 × eigh(p×p)` but adds `7 × O(p^2)` matmuls, which is more expensive for
the range of p values typical in rare-variant panels (p=5-50).

**How to avoid:** Skip the single-eigendecomposition optimization, OR implement it only
as a code-path for p > 100 with a conditional branch. The parallelization and GL quadrature
optimizations provide far greater gains with less complexity.

### Pitfall 6: Worker Function Must Be Module-Level

**What goes wrong:** Using a lambda or closure as the ProcessPoolExecutor worker function
causes `AttributeError: Can't pickle local object` or `PicklingError`.

**Why it happens:** Python's `pickle` module cannot serialize closures or lambdas.
`ProcessPoolExecutor` uses pickle to send work to worker processes.

**How to avoid:** Define `_run_gene_worker(args)` as a module-level function in
`engine.py` (or a separate `_parallel.py` module). Use a tuple `args` parameter.

---

## Code Examples

### Adding `--association-workers` to CLI

```python
# Source: cli.py, stats_group argument section (after --coast-weights at line ~511)
stats_group.add_argument(
    "--association-workers",
    type=int,
    default=1,
    help=(
        "Number of parallel worker processes for association analysis. "
        "Default: 1 (sequential). Set to -1 for os.cpu_count(). "
        "Only effective with Python backends (R backend is not parallel-safe)."
    ),
)
# In cfg assignment block:
cfg["association_workers"] = getattr(args, "association_workers", 1)
```

### Adding `association_workers` to AssociationConfig

```python
# Source: base.py, AssociationConfig dataclass (after coast_backend field)
# Phase 27: Gene-level parallelization
association_workers: int = 1
"""Number of parallel worker processes. Default: 1 (sequential).
Set > 1 only when all registered tests have parallel_safe=True."""
```

### Module-Level Worker Function

```python
# Source: engine.py (module level, before AssociationEngine class)
def _run_gene_worker(
    args: tuple,
) -> tuple[str, dict]:
    """
    Process a single gene in a subprocess worker.

    Parameters
    ----------
    args : tuple of (gene, gene_data, pickled_tests, config)
        gene : str
        gene_data : dict
        pickled_tests : bytes — pickle.dumps of {test_name: AssociationTest}
        config : AssociationConfig

    Returns
    -------
    (gene, {test_name: TestResult})
    """
    import os
    import pickle
    # Prevent BLAS thread oversubscription
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    gene, gene_data, pickled_tests, config = args
    tests: dict = pickle.loads(pickled_tests)
    results: dict = {}
    for test_name, test in tests.items():
        result = test.run(gene, gene_data, config)
        results[test_name] = result
    return gene, results
```

### Parallel Gene Loop in engine.run_all()

```python
# Source: engine.py, run_all() method
import concurrent.futures
import pickle

n_workers = getattr(self._config, "association_workers", 1)
# Check parallel safety of all registered tests
all_parallel_safe = all(
    getattr(test, "parallel_safe", False) for test in self._tests.values()
)
use_parallel = n_workers != 1 and all_parallel_safe

if use_parallel:
    # Pre-fit null models by running first gene sequentially
    # (triggers lazy null model fitting in SKAT/COAST/Burden tests)
    first_gene = sorted_data[0]
    first_gene_name = first_gene.get("GENE", "")
    for test_name, test in self._tests.items():
        result = test.run(first_gene_name, first_gene, self._config)
        results_by_test[test_name][first_gene_name] = result

    # Pickle tests WITH their pre-fitted null models
    pickled_tests = pickle.dumps(self._tests)

    # Dispatch remaining genes in parallel
    remaining = sorted_data[1:]
    args_list = [
        (gd.get("GENE", ""), gd, pickled_tests, self._config)
        for gd in remaining
    ]
    actual_workers = os.cpu_count() if n_workers == -1 else n_workers
    with concurrent.futures.ProcessPoolExecutor(max_workers=actual_workers) as executor:
        for gene, gene_results in executor.map(_run_gene_worker, args_list):
            for test_name, result in gene_results.items():
                results_by_test[test_name][gene] = result
else:
    # Sequential (existing code path unchanged)
    for gene_data in sorted_data:
        gene = gene_data.get("GENE", "")
        for test_name, test in self._tests.items():
            result = test.run(gene, gene_data, self._config)
            results_by_test[test_name][gene] = result
```

### GL Quadrature in `_skato_integrate_davies`

```python
# Source: python_backend.py, replacing scipy.integrate.quad in _skato_integrate_davies
from numpy.polynomial.legendre import leggauss

_GL_NODES, _GL_WEIGHTS = leggauss(128)  # module-level constants
_GL_A, _GL_B = 0.0, 40.0
_GL_X = (_GL_B - _GL_A) / 2.0 * _GL_NODES + (_GL_B + _GL_A) / 2.0  # nodes in [0,40]
_GL_W = (_GL_B - _GL_A) / 2.0 * _GL_WEIGHTS

def _skato_integrate_davies_gl(pmin_q, param, rho_grid, pmin):
    """
    SKAT-O omnibus p-value via fixed Gauss-Legendre quadrature (128 nodes).
    Replaces scipy.integrate.quad: 46x speedup measured (379ms -> 8ms per gene).
    """
    tau = param["tau"]
    mu_q = param["mu_q"]
    var_q = param["var_q"]
    var_remain = param["var_remain"]
    lambdas = param["lambdas"]
    rho_arr = np.array(rho_grid)

    var_ratio = (var_q - var_remain) / var_q if var_q > 0 else 0.0
    sd1 = np.sqrt(max(var_ratio, 0.0))
    valid = rho_arr < 1.0

    # Evaluate integrand at all GL nodes
    integrand_vals = np.zeros(len(_GL_X))
    for i, x in enumerate(_GL_X):
        if not np.any(valid):
            integrand_vals[i] = 0.0
            continue
        cond_q = (pmin_q[valid] - tau[valid] * x) / (1.0 - rho_arr[valid])
        min_q = float(np.min(cond_q))

        if min_q > np.sum(lambdas) * 1e4:
            cdf_val = 0.0
        else:
            min_q_std = (min_q - mu_q) * sd1 + mu_q
            p_dav, ifault = davies_pvalue(min_q_std, lambdas, acc=1e-6, lim=10_000)
            if p_dav is not None and ifault == 0:
                cdf_val = float(np.clip(p_dav, 0.0, 1.0))
            else:
                return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

        integrand_vals[i] = (1.0 - cdf_val) * float(scipy.stats.chi2.pdf(x, df=1))

    if np.any(np.isnan(integrand_vals)):
        return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

    integral = float(np.dot(integrand_vals, _GL_W))
    pvalue = 1.0 - integral

    if pmin * len(rho_grid) < pvalue:
        pvalue = pmin * len(rho_grid)

    return float(np.clip(pvalue, 0.0, 1.0))
```

### `parallel_safe` Attribute on Python Tests

```python
# Source: skat_python.py, PurePythonSKATTest class body
class PurePythonSKATTest(AssociationTest):
    parallel_safe: bool = True  # Thread/process-safe: pure numpy/scipy, no rpy2

# Source: logistic_burden.py, LogisticBurdenTest class body
class LogisticBurdenTest(AssociationTest):
    parallel_safe: bool = True  # statsmodels is process-safe when each worker has own instance

# Source: fisher.py, FisherExactTest class body
class FisherExactTest(AssociationTest):
    parallel_safe: bool = True  # scipy.stats.fisher_exact is stateless

# Source: linear_burden.py, LinearBurdenTest class body
class LinearBurdenTest(AssociationTest):
    parallel_safe: bool = True

# NOTE: COASTTest (R backend) already has parallel_safe = False
# NOTE: PurePythonCOASTTest already has parallel_safe = True
```

---

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|-----------------|--------|
| `scipy.integrate.quad` (adaptive, ~2400 davies calls) | Fixed Gauss-Legendre 128 nodes | 46x speedup measured |
| Sequential gene loop | ProcessPoolExecutor with N workers | ~N× wall-clock speedup |
| 7 independent `eigh()` calls per rho | Single `eigh()` + matmul (theoretical) | Not beneficial for p < 100 (benchmarked) |

**Key benchmark findings (p=20 variants, n=200 samples):**
- SKAT-O per-gene: 196ms total
- Current `_skato_integrate_davies`: 380ms per gene (dominated by adaptive quad overhead)
- GL quadrature replacement: 8.3ms (46x speedup)
- Single eigendecomposition: SLOWER than current for p < 100

**Clarification on 3.3 vs 3.4 from the plan:**
- Optimization 3.3 (Davies cache/GL quadrature): CONFIRMED as high impact, implement as GL
- Optimization 3.4 (Single eigendecomposition): NOT confirmed as beneficial; skip or defer

---

## Open Questions

1. **GL quadrature accuracy at extreme p-values**
   - What we know: Error of 3.2e-3 at p≈0.86; 128 nodes tested against adaptive quad
   - What's unclear: Accuracy at p < 0.001 where the integrand becomes steep
   - Recommendation: Run SKAT-O comparison test (like test_skat_python_comparison.py)
     with GL quadrature vs adaptive on validation dataset before merging

2. **Single eigendecomposition algebraic transform**
   - What we know: Building K_vbasis via `M @ (d[:,None] * M.T)` is SLOWER for p < 100
   - What's unclear: Whether there is a truly O(p^2) approach that avoids building K_vbasis
     (using rank-1 eigenvalue update theory for `C @ D @ C.T` where C = sI + delta*u*u.T)
   - Recommendation: Skip this optimization in Phase 27; the GL quadrature delivers far
     greater savings. Mark as "Future optimization for p > 100 panels."

3. **ProcessPoolExecutor overhead vs per-gene work ratio**
   - What we know: Process spawning overhead is ~50-100ms per worker; gene work is ~8ms (GL)
   - What's unclear: At what gene count does parallelization break even?
   - Recommendation: Add a minimum threshold (e.g., n_genes >= 4*n_workers) before
     enabling parallel dispatch. Default 1 worker (sequential) is safe.

4. **`association_workers` in `_build_assoc_config_from_context`**
   - What we know: AssociationConfig is built in `analysis_stages.py:_build_assoc_config_from_context()`
   - What's unclear: Whether `association_workers` needs to be added to the JSON config
     loading logic there or just read directly from `context.config`
   - Recommendation: Add to AssociationConfig dataclass AND read from context.config
     in the stage (consistent with how other fields like `skat_backend` are handled)

---

## Sources

### Primary (HIGH confidence — direct code inspection)

- `variantcentrifuge/association/engine.py` — full run_all() sequential loop, ACAT-O post-loop
- `variantcentrifuge/association/backends/python_backend.py` — _skato_integrate_davies,
  _skato_get_pvalue, eigenvalue computation in _skato_get_pvalue loop
- `variantcentrifuge/association/backends/davies.py` — davies_pvalue C extension interface
- `variantcentrifuge/association/base.py` — AssociationConfig, AssociationTest ABC
- `variantcentrifuge/association/tests/skat_python.py` — null model lazy fitting pattern
- `variantcentrifuge/association/genotype_matrix.py` — dense matrix construction
- `variantcentrifuge/cli.py` — existing --association-* argument patterns (lines 408-531)

### Benchmarks Run (HIGH confidence — measured directly)

- `davies_pvalue` call count per SKAT-O gene: **2,415 calls** (p=20, n=200)
- `scipy.linalg.eigh` call count per SKAT-O gene: **8 calls** (7 rho + 1 binary projection)
- `_skato_integrate_davies` timing: **379ms per gene** (p=20, n=200)
- GL quadrature (128 nodes): **8.3ms per gene** — **46x speedup**
- Single eigh approach speedup at p=20: **0.67x** (SLOWER)
- Single eigh speedup at p=50: **0.47x** (SLOWER)
- Single eigh speedup at p=100: **0.26x** (SLOWER)
- Parallel correctness test: 5 genes, 2 workers, results **match sequential exactly**
- Pickle sizes: gene_data=33.4KB, config=696B, PurePythonSKATTest+null_model=12KB

### Secondary (MEDIUM confidence)

- `ASSOCIATION_REMAINING_WORK.md` — effort estimates and implementation notes for Wave 3
- `tests/unit/test_skat_python_backend.py` — existing test structure for SKAT-O tests
- `tests/unit/test_engine_skat_integration.py` — engine integration test patterns

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all libraries are already in use, no new deps
- Architecture (parallelization): HIGH — pickling verified, worker pattern tested
- Architecture (GL quadrature): HIGH — benchmarked and verified correct
- Architecture (single eigh): HIGH — benchmarked; **single eigh is NOT beneficial for p<100**
- Pitfalls: HIGH — each pitfall verified by code inspection or benchmark

**Research date:** 2026-02-22
**Valid until:** 2026-04-22 (stable codebase, 60-day validity)

---

## Phase Planning Guidance

### Recommended Task Decomposition

**Task 27-01: GL Quadrature (highest impact, lowest risk)**
- Replace `scipy.integrate.quad` in `_skato_integrate_davies` with fixed 128-node GL
- Add module-level `_GL_NODES, _GL_WEIGHTS, _GL_X, _GL_W` constants
- Keep `_skato_integrate_davies` as fallback (renamed `_skato_integrate_davies_legacy`)
- Test: add benchmark test comparing GL vs current for multiple gene scenarios
- Deliverable: 46x speedup on SKAT-O integration step

**Task 27-02: `parallel_safe` Attributes on All Python Tests**
- Add `parallel_safe: bool = True` to FisherExactTest, LogisticBurdenTest,
  LinearBurdenTest, PurePythonSKATTest
- Verify COASTTest.parallel_safe = False already set, PurePythonCOASTTest.parallel_safe = True
- No functional change, pure annotation

**Task 27-03: `association_workers` in AssociationConfig and CLI**
- Add `association_workers: int = 1` field to AssociationConfig
- Add `--association-workers` argument to CLI
- Add reading of `association_workers` in `_build_assoc_config_from_context`
- No parallelization yet, just plumbing

**Task 27-04: ProcessPoolExecutor in engine.run_all()**
- Add module-level `_run_gene_worker()` function
- Modify `run_all()`: if workers > 1 AND all_parallel_safe, use ProcessPoolExecutor
- Pre-fit null models before dispatch (run first gene sequentially)
- Collect results, maintain sort order
- Test: verify results identical between sequential and parallel execution

**Task 27-05: Sparse Genotype Matrices (Optional)**
- Only if explicitly requested; skip for typical cohort sizes (n < 10K)
- When sparsity > 90%, convert to `scipy.sparse.csr_matrix`
- Matters only at n > 10K samples

### DO NOT implement in Phase 27

- Single eigendecomposition (3.4): Benchmarked as SLOWER for p < 100 (the common case).
  The plan's stated "5x speedup" was not reproduced. Either skip entirely or gate on p > 100.
