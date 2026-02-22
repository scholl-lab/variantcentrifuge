# Association Framework — Remaining Work & Polish Plan

**Date:** 2026-02-22
**Scope:** Post-Phase 24 work to finalize v0.15.0 milestone
**Consolidates:** ASSOCIATION-FRAMEWORK-DESIGN.md (open items), SKAT_OPTIMIZATION_PLAN.md (all 7 opts)

---

## Current State

Phases 18-24 are complete. The association framework is **functionally complete**:
- 6 tests: Fisher, logistic burden, linear burden, SKAT, SKAT-O, COAST
- ACAT-O omnibus combination
- Dual backends: R (rpy2) and pure Python for both SKAT and COAST
- Covariates, PCA file loading, Beta/CADD/REVEL weights
- Diagnostics: lambda_GC, QQ data, QQ plot (matplotlib)
- FDR/Bonferroni correction
- 1932+ tests passing, CI clean

**What's missing:** R deprecation path, user documentation, performance optimizations,
and a few feature gaps (ACAT-V, JSON config, AKT/PLINK wrappers).

---

## Part 1: Make Python Default, R Optional

**Rationale:** Both SKAT and COAST now have fully validated pure Python backends.
Python backends are thread-safe (`parallel_safe=True`), have no external dependency,
and match R within 0.003 log10 on real data. R backends add deployment friction
(R + rpy2 + SKAT/AllelicSeries packages) with no accuracy advantage.

### 1.1 Swap Default Backend to Python

**Files:** `backends/__init__.py`, `engine.py`, `base.py`

Current `--skat-backend auto` and `--coast-backend auto` both try R first, fall back
to Python. Reverse this:

```
auto (current):  try R → fall back to Python
auto (new):      use Python (log "using Python backend")
                 if --skat-backend r: try R → error if unavailable
```

- Change `get_skat_backend("auto")` to return `PythonSKATBackend()` directly
- Change COAST auto-detection in `engine.py` similarly
- Keep `--skat-backend r` and `--coast-backend r` as explicit opt-in
- Add deprecation warning when R backend is used:
  `"R SKAT backend is deprecated and will be removed in v0.17.0. Use --skat-backend python (default)."`
- Update `AssociationConfig` default: `skat_backend="python"`, `coast_backend="python"`

### 1.2 Mark R Backend as Deprecated

**Files:** `backends/r_backend.py`, `tests/skat_r.py`, `tests/allelic_series.py`

- Add `warnings.warn("RSKATBackend is deprecated...", DeprecationWarning)` in
  `RSKATBackend.__init__()` and `COASTTest.__init__()`
- Add `# DEPRECATED` docstring header to both classes
- Keep all R backend code functional — do not remove yet

### 1.3 Future: Remove R Backends (v0.17.0)

Not for this milestone. Plan for v0.17.0:
- Delete `r_backend.py`, `skat_r.py`, `allelic_series.py` (R wrapper)
- Remove rpy2 from optional dependencies
- Remove `--skat-backend r` and `--coast-backend r` CLI args

---

## Part 2: Documentation

**Current state:** Association testing is entirely undocumented outside of CLI `--help`
text and source docstrings. No guide, no README mention, no API docs, no examples.

### 2.1 Association Testing Guide

**New file:** `docs/source/guides/association_testing.md`

Sections:
1. **Overview** — what the framework does, when to use `--perform-association` vs
   `--perform-gene-burden`, available tests
2. **Quick Start** — minimal CLI examples:
   - Fisher only (backward compatible)
   - Logistic burden with covariates
   - SKAT-O with PCA
   - COAST allelic series
   - Full pipeline: all tests + ACAT-O omnibus
3. **Available Tests** — table: test name, what it tests, supports covariates?,
   supports weights?, output columns
4. **Covariates** — file format (TSV, first col = sample ID), `--covariates` column
   selection, categorical auto-detection, `--categorical-covariates`
5. **PCA Integration** — `--pca-file` formats (PLINK .eigenvec, AKT, generic TSV),
   `--pca-components`, how PCs merge with covariates
6. **Variant Weights** — `--variant-weights beta:1,25` / `uniform` / `cadd` / `revel` /
   `combined`, when to use each
7. **COAST Allelic Series** — what BMV/DMV/PTV classification means, required annotation
   columns (`dbNSFP_SIFT_pred`, `dbNSFP_Polyphen2_HDIV_pred`), `--coast-weights`
8. **Diagnostics** — `--diagnostics-output`, lambda_GC interpretation, QQ plots
9. **Output Format** — column naming convention (`{test}_p_value`, `{test}_effect_size`),
   ACAT-O columns, warnings column
10. **Backend Selection** — Python (default), R (deprecated opt-in), `--skat-backend`,
    `--coast-backend`

### 2.2 Update Existing Docs

**Files to update:**

| File | Change |
|------|--------|
| `docs/source/usage.md` | Add `--perform-association` to Statistical Analysis table, link to guide |
| `docs/source/guides/cohort_analysis.md` | Add section on modular association tests after Gene Burden Testing |
| `docs/source/faq.md` | Add FAQ entries: "How do I run SKAT-O?", "Do I need R?", "How do I add covariates?" |
| `docs/source/api/index.md` | Add `association` module to toctree |
| `docs/source/index.md` | Update feature list: "Modular association testing: SKAT-O, COAST, burden, ACAT-O" |
| `README.md` | Update feature list, add association testing bullet |

### 2.3 API Reference Stubs

**New files:**
- `docs/source/api/association.md` — autodoc for `variantcentrifuge.association`
- Cover: `engine.py`, `base.py`, `correction.py`, `covariates.py`, `pca.py`,
  `weights.py`, `diagnostics.py`, `backends/`, `tests/`

### 2.4 Changelog

**File:** `docs/source/changelog.md`

Add v0.15.0 section covering all phases 18-24:
- Modular association engine with 6 test types
- Pure Python SKAT/SKAT-O backend (no R required)
- Pure Python COAST backend (no R required)
- ACAT-O omnibus p-value combination
- Covariate adjustment with auto-categorical detection
- PCA file loading (PLINK, AKT, generic)
- CADD/REVEL/combined variant weights
- Diagnostics: lambda_GC, QQ plots
- R backends deprecated in favor of Python

---

## Part 3: Performance Optimizations

Based on SKAT_OPTIMIZATION_PLAN.md and research into current best practices
(regenie, SAIGE-GENE+, FastSKAT, Das 2024, MCMC-CE 2025).

### 3.1 Saddlepoint Before Liu Fallback (Quick Win)

**Impact:** Accuracy improvement for extreme tails (p < 1e-8)
**Effort:** ~10 lines in `davies.py`
**Risk:** None

Current fallback when Davies C ext is available: `Davies → Liu`.
`_kuonen_pvalue()` already exists but is only used when C ext is absent.

Change `compute_pvalue()` so when Davies returns out-of-range (p > 1 or p <= 0):
1. Try saddlepoint first
2. Fall back to Liu only if saddlepoint also fails

This matches SAIGE-GENE+ best practice: SPA is preferred over Liu for extreme tails
because Liu systematically overestimates p-values when eigenvalues are non-uniform.

### 3.2 Gene-Level Parallelization

**Impact:** ~Nx wall-clock speedup for multi-gene panels (N = worker count)
**Effort:** ~50 lines in `engine.py`
**Risk:** Low — Python SKAT/COAST are already `parallel_safe=True`

Current `engine.run_all()` loops genes sequentially.

Best practice from regenie/STAARpipeline: **process-level parallelism over genes**,
with `OPENBLAS_NUM_THREADS=1` per worker to avoid BLAS thread oversubscription.

Implementation:
- Add `--association-workers N` CLI arg (default: 1, no parallelism)
- Use `concurrent.futures.ProcessPoolExecutor(max_workers=N)` in `engine.run_all()`
- Set `os.environ["OPENBLAS_NUM_THREADS"] = "1"` in worker initializer
- Null model is fit once (serializable via pickle), shared across workers
- Each worker gets: gene name + genotype matrix + null model + config
- Collect results, sort by gene name, continue to ACAT-O/FDR

ProcessPoolExecutor (not ThreadPoolExecutor) because:
- Avoids GIL contention on Python-level code between BLAS calls
- Avoids BLAS thread pool oversubscription
- Workers are fully independent (no shared mutable state)
- Matches how regenie and STAARpipeline parallelize gene tests

### 3.3 Cache/Interpolate Davies in Omnibus Integration

**Impact:** ~3-5x speedup on SKAT-O omnibus integration step
**Effort:** ~40 lines in `python_backend.py`
**Risk:** Low — interpolation error < 1e-10 for smooth CDF

Current `_skato_integrate_davies()` calls `davies_pvalue()` up to ~1000 times
per gene (once per `scipy.integrate.quad` function evaluation).

Replace with:
1. Precompute Davies CDF at 200 points spanning the relevant Q range
2. Build cubic interpolant via `scipy.interpolate.interp1d`
3. Use interpolated CDF in integrand (zero FFI calls during integration)

Alternative (simpler): fixed Gauss-Legendre quadrature with 128 nodes,
vectorized evaluation of all nodes at once. Eliminates adaptive overhead
and gives exact control over evaluation count.

### 3.4 Single Eigendecomposition for SKAT-O

**Impact:** ~5x speedup on SKAT-O eigenvalue step (7 eigh calls → 1)
**Effort:** ~80 lines in `python_backend.py`
**Risk:** Medium — requires careful linear algebra

Current code computes 7 separate `scipy.linalg.eigh()` calls (one per rho).
Since `K_sym(rho) = R.M(rho)^{1/2} · A · R.M(rho)^{1/2}` and A = Z1'Z1 is fixed,
decompose A once and transform eigenvalues algebraically for each rho.

The analytical R.M^{1/2} formula (already implemented) means the transformed
matrix `V' · R.M^{1/2} · V` can be computed in O(p²) per rho, where V are
eigenvectors of A.

### 3.5 Sparse Genotype Matrices (Optional, Large Cohorts)

**Impact:** Memory and speed for n > 10K samples with rare variants (MAF < 1%)
**Effort:** ~60 lines in `genotype_matrix.py`, `python_backend.py`
**Risk:** Low — opt-in path, dense path unchanged

When sparsity > 90%, convert to `scipy.sparse.csr_matrix`. Benefits:
- Score vector `Z'r` is sparse mat-vec: O(nnz) not O(np)
- Kernel `Z'Z` construction: O(nnz·p) not O(n·p²)
- Memory: O(nnz) not O(np)

Per research: matters at n > 10K samples. Below that, dense is faster due
to sparse overhead. For biobank scale (n > 100K), essential.

### 3.6 ACAT-V Per-Variant Score Test

**Impact:** Completes the ACAT-O omnibus (adds sparse-signal detection)
**Effort:** ~40 lines in `acat.py`
**Risk:** None

ACAT-V combines per-variant marginal score test p-values via Cauchy combination.
Powerful when only a small fraction of variants are causal (where SKAT and burden
both lose power). Computationally trivial — O(M) where M = variants.

Per research (Liu 2019, STAARpipeline): ACAT-V is recommended as a component
of ACAT-O, not standalone. The omnibus becomes:
`ACAT-O = cauchy_combine(p_fisher, p_burden, p_skat, p_acat_v)`

Implementation:
- For each variant j: score test p-value `p_j = 2 * Phi(-|S_j|/sqrt(V_j))`
  where S_j = G_j' @ residuals and V_j = G_j' @ V @ G_j
- `p_acat_v = cauchy_combination(p_1, ..., p_M, weights=beta(MAF))`
- Add as a component in `compute_acat_o()`

---

## Part 4: Minor Feature Gaps

### 4.1 JSON Config for Association (Low Priority)

Phase 23-04 planned JSON config mode but the implementation scope was reduced.
Currently `--association-config` is not a CLI arg. The `config.json` has no
`association` section.

Defer to post-v0.15.0. CLI args are sufficient for all current use cases.
Users who want reproducible configs can use shell scripts.

### 4.2 AKT/PLINK Pipeline Stage Wrappers (Low Priority)

`pca.py` loads pre-computed PCA files from any tool. The `--pca-tool akt` flag
exists in CLI but the actual AKT invocation stage is minimal. Users run AKT/PLINK
externally and provide the eigenvec file.

Defer to post-v0.15.0. The file loader covers the primary use case.

### 4.3 Kinship Matrix / Mixed Models (Future)

Not in scope for v0.15.0. Would require SAIGE-style sparse GRM support.
Only relevant for biobank-scale with related individuals.

---

## Execution Order

### Wave 1: Defaults + Quick Wins (no API change, low risk)

| # | Task | Effort | Files |
|---|------|--------|-------|
| 1.1 | Swap default backend to Python | 30 min | `backends/__init__.py`, `engine.py`, `base.py` |
| 1.2 | Mark R backends deprecated | 15 min | `r_backend.py`, `skat_r.py`, `allelic_series.py` |
| 3.1 | Saddlepoint-before-Liu fallback | 15 min | `davies.py` |
| 3.6 | ACAT-V per-variant score test | 1 hr | `acat.py`, `engine.py` |

### Wave 2: Documentation

| # | Task | Effort | Files |
|---|------|--------|-------|
| 2.1 | Association testing guide | 2 hr | `docs/source/guides/association_testing.md` |
| 2.2 | Update existing docs | 1 hr | `usage.md`, `cohort_analysis.md`, `faq.md`, `index.md`, `README.md` |
| 2.3 | API reference stubs | 30 min | `docs/source/api/association.md`, `api/index.md` |
| 2.4 | Changelog v0.15.0 | 30 min | `docs/source/changelog.md` |

### Wave 3: Performance (independent, can be deferred)

| # | Task | Effort | Files |
|---|------|--------|-------|
| 3.2 | Gene-level parallelization | 2 hr | `engine.py`, `cli.py` |
| 3.3 | Davies cache/interpolation | 1 hr | `python_backend.py` |
| 3.4 | Single eigendecomposition | 3 hr | `python_backend.py` |
| 3.5 | Sparse genotype matrices | 2 hr | `genotype_matrix.py`, `python_backend.py` |

### Deferred (post-v0.15.0)

- 4.1 JSON config mode
- 4.2 AKT/PLINK wrappers
- 4.3 Kinship / mixed models
- 1.3 Remove R backends (v0.17.0)

---

## Total Estimated Effort

| Wave | Effort | Priority |
|------|--------|----------|
| Wave 1: Defaults + Quick Wins | ~2 hr | **Must-have for v0.15.0** |
| Wave 2: Documentation | ~4 hr | **Must-have for v0.15.0** |
| Wave 3: Performance | ~8 hr | Nice-to-have (defer if time-constrained) |

---

## References

- [regenie (Mbatchou 2021)](https://www.nature.com/articles/s41588-021-00870-7) — gene-level parallelization pattern
- [SAIGE-GENE+ (Zhou 2022)](https://www.nature.com/articles/s41588-022-01178-w) — SPA for extreme tails, sparse genotypes
- [FastSKAT (Lumley 2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6129408/) — partial eigenvalue decomposition
- [Das 2024](https://arxiv.org/abs/2404.05062) — new methods for generalized chi-square distribution
- [MCMC-CE (2025)](https://www.biorxiv.org/content/10.1101/2025.03.16.643492v1) — extreme tail p-values
- [Lumley 2024 blog](https://notstatschat.rbind.io/2024/09/13/two-approaches-to-approximating-sums-of-chisquareds/) — SPA vs Liu guidance
- [Liu & Xie 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/) — ACAT/ACAT-V
- [STAARpipeline](https://www.nature.com/articles/s41592-022-01641-w) — ACAT-O as recommended omnibus
