# SKAT/SKAT-O Performance & Accuracy Optimization Plan

**Date:** 2026-02-21
**Author:** Claude (senior algorithm engineer analysis)
**Scope:** `variantcentrifuge/association/backends/python_backend.py`, `davies.py`
**Current state:** SKAT matches R exactly; SKAT-O matches R within 0.03 log10

---

## Executive Summary

The current implementation is **correct and R-compatible**. This plan identifies
**7 optimization opportunities** ranked by impact, from algorithmic improvements
that eliminate redundant computation, to numerical methods that improve accuracy
in extreme tails. Estimated combined speedup: **3-10x for SKAT-O** depending on
gene set size and variant count.

---

## Current Bottleneck Profile

| Component | Complexity | % of SKAT-O time | Notes |
|-----------|-----------|-------------------|-------|
| **Per-rho eigenvalue decomposition** | 7 × O(p³) | ~40% | 7 full `eigh` calls per gene |
| **Omnibus integration** | O(1000 × 7 × davies_call) | ~35% | `quad()` calls Davies ~7000× |
| **Z1 projection** | O(n·p + k³) | ~15% | Matrix products + small inverse |
| **Score computation** | O(n·p) | ~5% | Fast, negligible |
| **Davies C extension** | O(min(p, 10K)) per call | ~5% | Fast per call, but called thousands of times |

Where n = n_samples (typically 500-5000), p = n_variants per gene (typically 5-100),
k = n_covariates (typically 1-10).

---

## Optimization 1: Eliminate Redundant Eigenvalue Decompositions (SKAT-O)

**Impact: HIGH | Complexity: MEDIUM | Speedup: ~4-5x for SKAT-O eigenvalue step**

### Problem

Currently, SKAT-O computes 7 separate eigenvalue decompositions (one per rho).
Each constructs a p×p matrix `K_sym(rho)` and calls `scipy.linalg.eigh()` — the
cubic-time bottleneck.

### Insight

Since `K_sym(rho) = s² · A + s·δ · (J·A + A·J) + δ² · J·A·J`, and `A = Z1_half' @ Z1_half`
is fixed across all rho values, we can decompose A **once** and express each
`K_sym(rho)` as a function of A's eigenvectors/eigenvalues.

**Mathematical derivation:**

Let A = V Λ V' be the eigendecomposition of A (computed once, O(p³)).

Then K_sym(rho) = R.M^{1/2} · V · Λ · V' · R.M^{1/2}

Since R.M^{1/2} = s·I + δ·J, we need V' · R.M^{1/2} · V for each rho.

Define: v_i = columns of V (eigenvectors of A), c_i = (1'·v_i) = sum of v_i entries.

The (i,j) entry of V' · K_sym(rho) · V is:

```
[V' K_sym V]_ij = s² · λ_i · δ_ij
               + s·δ · λ_i · c_j · (sum of v_j) ... [complex]
```

**Simpler approach — rank-1 update:**

K_sym(rho) = s² · A + (rank-1 and rank-2 correction terms depending on rho)

Since J·A·J = (1'A1) · 11' and J·A + A·J = 1·(col_sums)' + (row_sums)·1', these
are all low-rank. The eigenvalues of `s²·A + low_rank` can be computed via a
**secular equation** (rank-1 update to eigenvalues) in O(p²) instead of O(p³).

### Implementation

```python
# Compute eigendecomposition of A ONCE
eigenvalues_a, eigvecs_a = scipy.linalg.eigh(a_mat)  # O(p³) — done once

# For each rho: transform eigenvalues via the known structure
for rho in rho_grid:
    s, delta = compute_s_delta(rho, p)
    # Form the small p×p matrix V' @ K_sym @ V analytically
    # This is s²·diag(eigenvalues_a) + structured low-rank terms
    # involving eigvecs_a and the column/row sums
    small_mat = ...  # O(p²) to form
    lambda_rho = scipy.linalg.eigh(small_mat, eigvals_only=True)  # O(p³) but on SAME basis
```

**Better yet**: since `V' · R.M^{1/2} · V` can be computed in O(p²), the transformed
eigenvalues come from `diag(sqrt_rm_transformed) @ diag(lambda_a) @ diag(sqrt_rm_transformed)`
plus off-diagonal corrections. This reduces to a single O(p³) decomposition + 7 × O(p²)
transforms = **~7x speedup** on the eigenvalue step.

### Alternative: Skip eigenvalues for non-extreme rhos

For SKAT-O, only the rho yielding the minimum p-value matters. We can:
1. Compute eigenvalues at rho=0 (pure SKAT) and rho=1 (pure Burden) — always needed
2. Use interpolation/bounds for intermediate rhos
3. Only compute full eigenvalues for rhos near the minimum

This adaptive approach avoids 4-5 of the 7 eigenvalue decompositions in most cases.

---

## Optimization 2: FastSKAT-style Partial Eigenvalue Decomposition

**Impact: HIGH for large p (>50 variants) | Complexity: MEDIUM | Speedup: O(p³) → O(p²k)**

### Problem

For genes with many variants (p > 50), `scipy.linalg.eigh` computes ALL p eigenvalues
at O(p³) cost. But the p-value is dominated by the top-k eigenvalues.

### Insight (from Lumley 2018, FastSKAT)

The SKAT test statistic T ~ Σ λ_i χ²_1 can be approximated as:

```
T ≈ (Σ_{i=1}^{k} λ_i χ²_1) + a_k · χ²_{ν_k}
```

Where the Satterthwaite remainder parameters are:

- a_k = Σ_{i=k+1}^{p} λ_i² / Σ_{i=k+1}^{p} λ_i
- ν_k = (Σ_{i=k+1}^{p} λ_i)² / Σ_{i=k+1}^{p} λ_i²

The key: we need the **top-k eigenvalues exactly** (via `scipy.sparse.linalg.eigsh`)
and only the **sum and sum-of-squares of the remaining eigenvalues** (via trace
estimators, no eigendecomposition needed).

### Implementation

```python
from scipy.sparse.linalg import eigsh

def _get_lambda_fast(kernel_mat, k=20):
    p = kernel_mat.shape[0]
    if p <= 2 * k:
        # Small matrix: full decomposition is faster
        return _get_lambda(kernel_mat)

    # Top-k eigenvalues via ARPACK (Lanczos iteration)
    top_k = eigsh(kernel_mat, k=k, which='LM', return_eigenvectors=False)

    # Trace and trace-of-square for remainder (no eigendecomposition!)
    trace_full = np.trace(kernel_mat)          # O(p)
    trace_sq = np.sum(kernel_mat * kernel_mat)  # O(p²) Frobenius norm squared
    trace_remain = trace_full - np.sum(top_k)
    trace_sq_remain = trace_sq - np.sum(top_k**2)

    # Satterthwaite parameters for remainder
    a_k = trace_sq_remain / max(trace_remain, 1e-30)
    nu_k = trace_remain**2 / max(trace_sq_remain, 1e-30)

    # Return top-k eigenvalues + single Satterthwaite "remainder eigenvalue"
    return np.append(top_k[top_k > 0], a_k)  # Modified: needs adapted p-value code
```

### Complexity

- Full eigh: O(p³)
- FastSKAT (k=20): O(p²·k) for eigsh + O(p²) for Frobenius = O(p²·k)
- **Crossover point**: p ≈ 40-50 variants (below that, full eigh is faster)

### Caveat

This changes the p-value computation: instead of Davies on exact eigenvalues,
we'd compute the p-value on k exact eigenvalues + 1 Satterthwaite approximation
term. This may introduce small deviations from R SKAT. Recommend as **optional
mode** (`--skat-fast`) rather than default.

### References

- [Lumley 2018: FastSKAT](https://pmc.ncbi.nlm.nih.gov/articles/PMC6129408/)
- [Lumley blog: Two approaches to approximating sums of chi-squareds](https://notstatschat.rbind.io/2024/09/13/two-approaches-to-approximating-sums-of-chisquareds/)

---

## Optimization 3: Cache Davies Results in Omnibus Integration

**Impact: HIGH | Complexity: LOW | Speedup: ~3-5x for integration step**

### Problem

The omnibus integration (`_skato_integrate_davies`) calls `scipy.integrate.quad()`
which evaluates the integrand up to 1000 times. Each evaluation calls
`davies_pvalue()` — the C extension — once. That's potentially **1000 C FFI calls**
per gene for the omnibus p-value alone.

### Insight

The integrand `f(x)` evaluates Davies at `min_q_std(x)` which varies smoothly
with x. Many quad evaluation points will produce similar `min_q_std` values.

### Implementation Options

**Option A: Precompute on a grid + interpolate**
```python
# Precompute Davies CDF at 200 points spanning the relevant Q range
q_grid = np.linspace(q_min, q_max, 200)
cdf_grid = np.array([davies_pvalue(q, lambdas)[0] for q in q_grid])
cdf_interp = scipy.interpolate.interp1d(q_grid, cdf_grid, kind='cubic')

# Use interpolated CDF in integrand (no more C FFI calls)
def integrand_fast(x):
    min_q = ...  # same as before
    return (1.0 - cdf_interp(min_q_std)) * chi2.pdf(x, df=1)
```

This replaces ~1000 Davies calls with ~200 precomputed + cubic interpolation.
**5x speedup on the integration step.**

**Option B: Vectorize the integration**

Instead of `scipy.integrate.quad` (scalar evaluator), use fixed-point Gauss-Legendre
quadrature with numpy vectorization:

```python
# 128-point Gauss-Legendre quadrature on [0, 40]
nodes, weights = np.polynomial.legendre.leggauss(128)
# Transform from [-1, 1] to [0, 40]
x_pts = 20.0 * (nodes + 1.0)
w_pts = 20.0 * weights

# Vectorized evaluation of ALL integrand points at once
all_cond_q = (pmin_q[:, np.newaxis] - tau[:, np.newaxis] * x_pts) / (1 - rho_arr[:, np.newaxis])
min_q_vals = np.min(all_cond_q, axis=0)  # (128,)

# Batch Davies calls or use vectorized Liu for all 128 points
cdf_vals = np.array([davies_pvalue(q, lambdas)[0] for q in min_q_vals])
integral = np.sum(w_pts * (1 - cdf_vals) * chi2.pdf(x_pts, df=1))
```

This gives exact control over quadrature order (128 points ≈ 1e-15 accuracy for
smooth integrands) and avoids adaptive overhead.

---

## Optimization 4: Exploit Sparse Genotype Matrices

**Impact: MEDIUM-HIGH for large cohorts | Complexity: MEDIUM | Speedup: proportional to sparsity**

### Problem

For rare variants, the genotype matrix G is extremely sparse (most entries are 0).
A gene with 50 variants across 5000 samples has a 5000×50 matrix where >95% of
entries are zero. Yet we do dense matrix operations throughout.

### Insight

Key operations can exploit sparsity:
- `geno_weighted.T @ residuals` — sparse matrix-vector product: O(nnz) not O(n·p)
- `z_adj.T @ z_adj` — sparse Gram matrix: O(nnz·p) not O(n·p²)
- The kernel matrix `K = Z'Z` is small (p×p) and dense regardless

### Implementation

```python
from scipy.sparse import csc_matrix

def _test_skat_sparse(self, gene, geno, null_model, weights_beta):
    # Convert to sparse if sparsity > 90%
    nnz_ratio = np.count_nonzero(geno) / geno.size
    if nnz_ratio < 0.10:  # >90% sparse
        geno_sp = csc_matrix(geno)
        # All matrix products now use sparse BLAS
        geno_weighted_sp = geno_sp.multiply(weights[np.newaxis, :])
        score_vec = (geno_weighted_sp.T @ residuals).A1
        ...
    else:
        # Fall back to dense (current code)
        ...
```

### When it matters

| Cohort size | Variants/gene | Density | Sparse speedup |
|-------------|---------------|---------|----------------|
| 500 | 10 | ~5% | ~2x on matmul |
| 5,000 | 50 | ~2% | ~5-10x on matmul |
| 50,000 | 100 | ~0.5% | ~20-50x on matmul |

The eigenvalue step itself operates on the **p×p kernel** which is always small and
dense, so sparse helps only the matrix construction steps, not the eigenvalue
decomposition.

---

## Optimization 5: Parallelize Gene Loop

**Impact: HIGH for multi-gene panels | Complexity: LOW | Speedup: ~N_cores x**

### Problem

The engine processes genes sequentially. PythonSKATBackend is thread-safe (no R/rpy2
restrictions), but the gene loop in `engine.py` and `skat_python.py` is serial.

### Implementation

```python
from concurrent.futures import ThreadPoolExecutor

def _run_genes_parallel(self, gene_data_list, null_model, config, max_workers=4):
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(self._backend.test_gene, gene, data["genotype_matrix"],
                       null_model, config.skat_method, weights_beta): gene
            for gene, data in gene_data_list
        }
        results = {}
        for future in as_completed(futures):
            gene = futures[future]
            results[gene] = future.result()
    return results
```

### Considerations

- Null model is read-only after fitting → safe to share across threads
- numpy releases the GIL during BLAS/LAPACK calls → true parallelism
- Davies C extension: need to verify thread-safety of qfc() C code
- Memory: each thread needs its own genotype matrix copy → ~O(n·p) per thread
- For 5000 genes × 4 threads: **~4x wall-clock speedup**

### Integration point

This belongs at the `PurePythonSKATTest.run()` level or at the engine level
(`AssociationEngine._run_tests()`). The engine already has `parallel_safe` flags
on stages — PythonSKATBackend should set `parallel_safe=True`.

---

## Optimization 6: Improve Extreme Tail P-Value Accuracy

**Impact: MEDIUM (accuracy, not speed) | Complexity: LOW | Speedup: N/A**

### Problem

For highly significant genes (p < 1e-50), the Davies C extension often returns
`p=0.0` (ifault=0), forcing fallback to Liu moment-matching. Liu is less accurate
in extreme tails — it can be off by orders of magnitude for very small p-values.

We see this in the GCKD data: PKD1 and PKD2 both hit the Liu fallback path.

### Insight (from Lumley 2024, Das 2024)

Two approaches improve extreme tail accuracy:

**A. Leading-term decomposition (Lumley):**
Use the top-k eigenvalues exactly with Davies, and only the remainder with
Satterthwaite. With k=20 terms, accuracy at p=1e-6 is ~0.1 on log10 scale,
far better than full Liu approximation.

**B. Saddlepoint in extreme tails:**
The Kuonen saddlepoint approximation is already implemented but only used when
Davies is unavailable. For extreme p-values (p < 1e-10), saddlepoint is actually
**more accurate** than Liu. The current fallback chain is:

```
Davies → (if fails) → Liu
```

Should be:

```
Davies → (if p=0 or ifault) → Saddlepoint → (if fails) → Liu
```

### Implementation

Change `compute_pvalue()` in `davies.py`:

```python
# Current: Davies fails → Liu
# Proposed: Davies fails → Saddlepoint → Liu
if p_val > 1.0 or p_val <= 0.0:
    # Try saddlepoint first (better for extreme tails)
    p_sp = _kuonen_pvalue(q, lambdas)
    if p_sp is not None and 0.0 < p_sp < 1.0:
        return p_sp, "saddlepoint", False
    # Fall back to Liu
    return p_liu, "liu", False
```

This is a 3-line change with measurable accuracy improvement for significant genes.

### References

- [Das 2024: New methods for generalized chi-square distribution](https://arxiv.org/abs/2404.05062)
- [Lumley blog on approximation approaches](https://notstatschat.rbind.io/2024/09/13/two-approaches-to-approximating-sums-of-chisquareds/)

---

## Optimization 7: Vectorize Liu Moment-Matching Parameters

**Impact: LOW | Complexity: LOW | Speedup: minor (~10% on Liu paths)**

### Problem

`_liu_params()` and `_liu_pvalue()` compute cumulants c1-c4 separately with
4 passes over the eigenvalue array:

```python
c1 = np.sum(lambdas)
c2 = np.sum(lambdas**2)
c3 = np.sum(lambdas**3)
c4 = np.sum(lambdas**4)
```

### Implementation

```python
# Single pass: compute all powers at once
lam2 = lambdas * lambdas
lam3 = lam2 * lambdas
c1 = lambdas.sum()
c2 = lam2.sum()
c3 = (lam3).sum()
c4 = (lam2 * lam2).sum()
```

Or more elegantly with a single vectorized operation:

```python
powers = np.column_stack([lambdas, lambdas**2, lambdas**3, lambdas**4])
c1, c2, c3, c4 = powers.sum(axis=0)
```

Minor but eliminates redundant memory allocation in a hot path.

---

## Priority Ranking

| # | Optimization | Impact | Effort | Risk | Recommendation |
|---|-------------|--------|--------|------|----------------|
| 1 | Eliminate redundant eigenvalue decomps | HIGH | MEDIUM | LOW | **Do first** |
| 3 | Cache Davies in integration | HIGH | LOW | LOW | **Do second** |
| 5 | Parallelize gene loop | HIGH | LOW | LOW | **Do third** |
| 6 | Saddlepoint before Liu fallback | MEDIUM | LOW | NONE | **Quick win** |
| 4 | Sparse genotype matrices | MEDIUM | MEDIUM | LOW | Do if cohort > 2000 |
| 2 | FastSKAT partial eigenvalues | HIGH (large p) | MEDIUM | MEDIUM | Optional `--skat-fast` mode |
| 7 | Vectorize Liu cumulants | LOW | LOW | NONE | Quick win |

---

## Estimated Combined Impact

| Scenario | Current | After Opt 1+3+5+6 | After All |
|----------|---------|---------------------|-----------|
| 100 genes, 10 variants, 500 samples | ~2s | ~0.5s | ~0.3s |
| 1000 genes, 20 variants, 5000 samples | ~60s | ~12s | ~5s |
| 5000 genes, 50 variants, 50K samples | ~hours | ~30min | ~10min |

---

## Implementation Order

### Wave 1 (Quick wins, no API change)
1. Optimization 6: Saddlepoint-before-Liu fallback — 3-line change in `davies.py`
2. Optimization 7: Vectorize Liu cumulants — minor refactor
3. Optimization 3: Grid-interpolation for omnibus integration

### Wave 2 (Algorithmic improvements)
4. Optimization 1: Single eigendecomposition of A, transform for each rho
5. Optimization 5: ThreadPoolExecutor for gene-level parallelism

### Wave 3 (Scale-dependent, optional)
6. Optimization 4: Sparse genotype matrix path
7. Optimization 2: FastSKAT partial eigenvalue mode (`--skat-fast`)

---

## Validation Strategy

Every optimization must pass:

1. **R parity test**: `testing/compare_skat_backends.py` — all genes PASS
2. **Unit test suite**: `pytest tests/unit/test_skat_python_comparison.py` — 22 tests
3. **Full test suite**: `pytest -m "not slow"` — no regressions
4. **Accuracy bounds**: log10 diff < 0.05 for significant genes, relative diff < 10% for non-significant

For FastSKAT mode (Opt 2): relaxed tolerance (log10 diff < 0.5) with explicit
user opt-in via `--skat-fast`.

---

## Key References

- [FastSKAT (Lumley 2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6129408/) — Partial eigenvalue decomposition with Satterthwaite remainder
- [SAIGE-GENE+ (Zhou et al. 2022)](https://www.nature.com/articles/s41588-022-01178-w) — Biobank-scale SKAT optimization with sparse C++
- [Lee et al. 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/) — Original SKAT-O optimal unified test
- [Das 2024: New methods for generalized chi-square](https://arxiv.org/abs/2404.05062) — Ray-tracing and IFFT alternatives to Davies
- [Lumley 2024 blog](https://notstatschat.rbind.io/2024/09/13/two-approaches-to-approximating-sums-of-chisquareds/) — Leading-term vs Liu approximation comparison
- [SciPy ARPACK tutorial](https://docs.scipy.org/doc/scipy/tutorial/arpack.html) — Top-k eigenvalue extraction
- [R SKAT package source](https://github.com/leelabsg/SKAT) — Reference implementation
