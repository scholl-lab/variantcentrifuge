# Stack Research: v0.16.0 Association Hardening & Multi-Cohort Features

**Project:** VariantCentrifuge
**Researched:** 2026-02-23
**Scope:** Stack additions/changes ONLY for NEW v0.16.0 features. Existing stack
(scipy 1.14.1, statsmodels 0.14.4, numpy 2.2.6, pandas, rpy2>=3.6.4) is validated
and not re-evaluated here.

---

## Current Pinned Versions (Installed as of v0.15.0)

| Library | Installed | Role |
|---------|-----------|------|
| scipy | 1.14.1 | linalg, stats, sparse |
| numpy | 2.2.6 | array operations |
| statsmodels | 0.14.4 | GLM null model, multipletests |

---

## Feature 1: Gene-Level Prior Weighting for FDR

### Weighted Benjamini-Hochberg

**Verdict: No new dependency. Implement directly in numpy (~12 lines).**

`statsmodels.stats.multitest.multipletests` does NOT have a `weights` parameter.
Confirmed against statsmodels 0.15.0 dev documentation. Parameters accepted: `pvals`,
`alpha`, `method`, `maxiter`, `is_sorted`, `returnsorted` — no weights.

`scipy.stats.false_discovery_control` (added in scipy 1.11) also has no weights parameter.

Weighted BH is a short, well-understood algorithm (Genovese et al. 2006):

1. Normalize weights to sum to `m` (number of tests): `w_norm = w * m / sum(w)`
2. Scale p-values: `p_adj_input = p * m / w_norm` (equivalent: `p / (w / sum(w))`)
3. Apply standard BH to scaled p-values
4. Cap at 1.0

This is ~12 lines of numpy. No new dependency is warranted.

**Integration point:** Extend `variantcentrifuge/association/correction.py`.
Add `apply_weighted_bh(pvals, weights)` alongside existing `apply_correction()`.
The `weights` array should represent gene-level priors (e.g., derived from pathway
membership, prior association evidence, or gene size). Weights are normalized
internally — callers pass un-normalized priors.

**Confidence:** HIGH. Verified against statsmodels and scipy official documentation.

### IHW (Independent Hypothesis Weighting)

**Verdict: No Python implementation exists. Do NOT add for v0.16.0.**

IHW exists only as a Bioconductor R package (current: Bioconductor 3.20). As of
February 2026, there is no maintained Python port on PyPI or conda-forge.

The IHW algorithm requires iterative isotonic regression optimization and cross-validation
folds — approximately 200 lines of non-trivial statistical code with its own test burden.
Implementing it from scratch risks divergence from the canonical Bioconductor version.

| Option | Assessment |
|--------|-----------|
| Skip for v0.16.0 | Recommended. Weighted BH covers 90% of the use case with zero risk. |
| rpy2 bridge | Possible via existing optional `rpy2>=3.6.4`. Adds R + Bioconductor install requirement. Viable if demanded. |
| Custom Python port | Not recommended. Maintenance burden, correctness risk. |

**Recommendation:** Weighted BH with configurable gene-level priors is sufficient for
v0.16.0. Document IHW as a future rpy2-bridge option if users request it.

---

## Feature 2: Sample-Level Case-Confidence Weights

### Weighted GLM Null Model

**Verdict: No new dependency. `statsmodels.GLM` already supports sample weights.**

`statsmodels.GLM` accepts two weight parameters, confirmed in statsmodels 0.15.0 docs:

| Parameter | Type | Use Case |
|-----------|------|----------|
| `freq_weights` | 1D integer array | Replicate observations by count |
| `var_weights` | 1D float array | Analytic/inverse-variance weights |

For case-confidence weights (continuous values in (0, 1] representing certainty of
case/control status), `var_weights` is the correct parameter. An observation with
higher confidence gets lower variance and therefore more weight in the fit.

**Usage (no backend changes needed for the GLM itself):**

```python
glm = sm.GLM(
    phenotype,
    design_matrix,
    family=sm.families.Binomial(),
    var_weights=case_confidence_weights,  # shape (n_samples,), values in (0, 1]
)
fit_result = glm.fit()
```

**Official caveat:** "Using weights is not verified yet for all possible options and
results." Robust covariance types (HC0, HC3) are unverified. Standard covariance
(the default used by the null model) is verified.

**Integration point:** `PythonSKATBackend.fit_null_model()` in
`variantcentrifuge/association/backends/python_backend.py`. Add optional
`sample_weights: np.ndarray | None = None` parameter to `fit_null_model()`. When
provided, pass as `var_weights=sample_weights` to `sm.GLM`. The fitted `mu_hat`
and residuals propagate through to `test_gene()` unchanged — the eigenvalue
computation uses fitted values which will reflect the weighting automatically.

### Weighted SKAT Score Statistic

For SKAT with per-sample weights, the score statistic is modified per the SKAT
literature (Wu et al. 2011):

```
Q_w = r_w^T G W_v G^T r_w / 2
```

Where `r_w = sqrt(sample_weights) * residuals` and the kernel eigenvalues are
computed on `Z_w = diag(sqrt(sample_weights)) @ G_weighted`. This is pure numpy
broadcasting — no new dependencies.

**Confidence:** HIGH for `statsmodels.GLM` `var_weights` (verified in docs).
MEDIUM for weighted SKAT score modification (literature-based, requires validation
against R SKAT reference implementation).

---

## Feature 3: COAST Classification Scoring (Configurable SIFT/PolyPhen/CADD/REVEL)

**Verdict: No new dependency. Extend existing `weights.py` and `coast_python.py`.**

COAST classification maps per-variant annotations to severity tiers. The existing
`weights.py` already handles CADD/REVEL functional scoring with fallback logic.
SIFT and PolyPhen-2 annotation fields are available as dbNSFP columns already
extracted by SnpSift into the gene DataFrame.

The classification logic is pure Python/numpy: threshold lookups against annotation
columns with configurable cutoffs. No new library is needed.

**Integration point:** Add `coast_classification_weights()` to
`variantcentrifuge/association/weights.py`. Classification thresholds (e.g.,
SIFT < 0.05 = damaging, PolyPhen-2 > 0.908 = probably_damaging) should be
configurable via a new `coast_classification_config.json` under
`variantcentrifuge/scoring/coast_classification/`, parallel to the existing
`nephro_candidate_score/` pattern.

---

## Feature 4: Sparse Genotype Matrices for Large Cohorts

### Recommended Format: `scipy.sparse.csr_array`

**Verdict: Use `scipy.sparse.csr_array`. Available in currently installed scipy 1.14.1.
No new dependency.**

`csr_array` is available as of scipy 1.7 (new array API, confirmed present in scipy 1.14.1
from `import scipy.sparse; scipy.sparse.csr_array` — works without error).

**Format comparison for SKAT genotype matrix (n_samples x n_variants):**

| Format | G @ r (score vec) | G * w (variant weights) | Build from dosage data | Notes |
|--------|-------------------|------------------------|------------------------|-------|
| `csr_array` | Optimal (`@` op) | Efficient via `diags(w)` | Via `coo_array.tocsr()` | **Recommended** |
| `csc_array` | Slow | Fast column ops | Via `coo_array.tocsc()` | Use only if column-slice dominant |
| `coo_array` | Slow | Slow | Direct construction | Construction only, convert before use |
| dense `np.ndarray` | Fast (BLAS) | Fast (broadcasting) | Direct | Better for small windows (<200 samples x 2000 variants) |

**Key SKAT operation in sparse form:**

```python
from scipy.sparse import csr_array, diags

# G: csr_array, shape (n_samples, n_variants)
# weights: np.ndarray, shape (n_variants,)
# residuals: np.ndarray, shape (n_samples,)

# Score vector — central SKAT Q computation (G_weighted.T @ residuals)
W_diag = diags(weights)          # (n_variants, n_variants) diagonal sparse
G_weighted = G @ W_diag          # csr @ dia -> csr, efficient
score_vec = G_weighted.T @ residuals  # csr.T @ dense -> dense, optimal for CSR

# Q statistic
q_stat = float(score_vec @ score_vec) / 2.0
```

**API note — new array vs legacy matrix interface:**

scipy 1.14+ recommends `csr_array` (new array API) over `csr_matrix` (legacy).
Key difference: `*` is element-wise in the new API (NumPy semantics). Use `@` for
matrix multiplication. From scipy 1.15 release notes: "sparse.linalg and sparse.csgraph
now work with sparse arrays" — confirming that linear algebra operations are compatible
with the new API in the version immediately after currently installed scipy.

**Sparsity and breakeven analysis:**

For rare variants (MAF < 1%), expected carrier fraction ≈ 2*MAF ≈ 2% of entries are
non-zero. Memory comparison at float64:

| Cohort Size | Variants/Gene | Dense Memory | Sparse Memory (2% fill) | Ratio |
|-------------|---------------|-------------|------------------------|-------|
| 200 samples | 100 | 160 KB | ~7 KB | 23x |
| 1,000 samples | 1,000 | 8 MB | ~330 KB | 24x |
| 5,000 samples | 5,000 | 200 MB | ~8 MB | 25x |

Dense is faster for small matrices (BLAS-optimized contiguous memory). Sparse wins
on memory and is competitive on speed only when fill fraction < ~10%.

**Recommendation:** Implement sparse as an opt-in path controlled by a CLI flag
(`--sparse-genotype-matrix`) and/or an automatic threshold (enable sparse when
`n_samples * n_variants > 500_000`). Do NOT make sparse the default for v0.16.0 —
profiling data is needed first.

**Integration point:** `build_genotype_matrix()` in
`variantcentrifuge/association/genotype_matrix.py`. Add `sparse: bool = False`
parameter. When `True`, return `scipy.sparse.csr_array` instead of `np.ndarray`.
The SKAT backends need a type branch to handle both.

**Confidence:** HIGH for format choice and API availability. MEDIUM for breakeven
threshold (based on typical sparsity estimate, not benchmarked on this codebase).

---

## Feature 5: Single Eigendecomposition Optimization for SKAT-O

**Verdict: Algorithmic change only. No new dependency. Defer to post-profiling.**

### Current Behavior (v0.15.0)

`_skato_get_pvalue()` in `python_backend.py` calls `_get_lambda(k_sym)` once per
rho in the fixed 7-point grid, totaling 7 `scipy.linalg.eigh(driver='evr')` calls
per gene. The analytical `R.M^{1/2}` approach already avoids Cholesky instability.

### Potential Optimization

The base kernel `A = Z1_half.T @ Z1_half` is rho-independent. Computing `eigh(A)` once
and using rank-1 eigenvalue perturbation theory for each rho mixture could reduce
eigendecompositions from 7 to 1 per gene.

### Why This Should Be Deferred

From the v0.15.0 performance analysis (issue #76):

| Stage | Wall Time | Fraction |
|-------|-----------|----------|
| field_extraction (SnpSift) | 9,664s | 26% |
| genotype_replacement | 25,293s | 69% |
| All other stages | ~1,400s | 5% |

SKAT-O eigendecomposition on typical gene windows (p = 5–50 variants) involves
`eigh()` on a `(50, 50)` matrix — this takes ~10 microseconds. Even at 20,000 genes,
7 calls per gene = 140,000 calls * 10µs = 1.4 seconds total. This is not a measurable
bottleneck.

**Recommendation:** Note the optimization opportunity in a code comment. Implement
only after profiling shows SKAT-O eigendecomposition as a real bottleneck (i.e., after
the genotype_replacement bottleneck is addressed in a future milestone).

---

## Feature 6: PCAComputationStage Wiring (AKT/PLINK2 Subprocess)

**Verdict: No new Python dependency. Subprocess pattern already established in pipeline.**

### Current State

`variantcentrifuge/association/pca.py` already parses PCA output from three formats:
- PLINK `.eigenvec` with/without header
- AKT stdout / generic TSV

PCA file loading (`load_pca_file()`) and covariate merging (`merge_pca_covariates()`)
are fully implemented. What is missing: a pipeline stage that computes the PCA
file from the VCF before loading.

### Tool Choice: AKT as Primary, PLINK2 as Alternative

| Tool | Input | Notes |
|------|-------|-------|
| AKT | VCF directly | BCFtools plugin. Likely available if bcftools is installed. PCA output already parsed by `pca.py`. |
| PLINK2 | VCF or PLINK binary | Widely used in GWAS. Requires conversion step for VCF input. |

AKT is the right first implementation: it takes a VCF directly (no intermediate
conversion), its output format is already handled by the parser, and it follows the
same external-tool pattern as bcftools in the pipeline.

PLINK2 is a reasonable second backend for cohorts already in PLINK format.

**No new Python dependencies.** AKT and PLINK2 are external binaries, same pattern
as bcftools/SnpSift. Document as optional external tools, not Python packages.

**Integration point:** New `PCAComputationStage` in `variantcentrifuge/stages/analysis_stages.py`
(or a dedicated `stages/pca_stages.py`). The stage:

1. Declares dependency on `ConfigurationLoadingStage`
2. Runs `akt pca <vcf> ...` or `plink2 --pca ...` via `subprocess.run()`
3. Writes PCA output to workspace
4. Sets `context.pca_file_path` for downstream consumption by `GeneBurdenAnalysisStage`

Follow the subprocess exception handling pattern from `BCFToolsPrefilterStage` in
`processing_stages.py`.

**Confidence:** HIGH for subprocess approach and AKT format compatibility. MEDIUM
for AKT command-line API stability (AKT is a BCFtools plugin with less formal
versioning than PLINK2).

---

## Summary: Net Stack Changes for v0.16.0

**Zero new Python package dependencies required.**

| Feature | Change | Confidence |
|---------|--------|-----------|
| Weighted BH | Add `apply_weighted_bh()` to `correction.py` (numpy only) | HIGH |
| IHW | Defer; use rpy2 bridge later if needed | HIGH |
| Weighted GLM null model | Use existing `statsmodels.GLM(var_weights=...)` | HIGH |
| Weighted SKAT score | Numpy broadcasting extension to `_test_skat()` | MEDIUM |
| COAST classification | Extend `weights.py` + new config JSON | HIGH |
| Sparse genotype matrix | `scipy.sparse.csr_array` (in scipy 1.14.1) opt-in | HIGH |
| SKAT-O single eigh | Defer; document as future optimization | HIGH |
| PCA computation stage | New stage using subprocess (AKT/PLINK2 as external tools) | MEDIUM |

### Version Pinning

No version bumps required for any new feature. The installed stack (scipy 1.14.1,
statsmodels 0.14.4) supports all needed APIs:

- `scipy.sparse.csr_array`: available since scipy 1.7, confirmed in 1.14.1
- `statsmodels.GLM(var_weights=...)`: available since statsmodels 0.9.0
- `statsmodels.stats.multitest.multipletests`: no `weights` param, confirmed; weighted BH implemented in numpy instead

**If scipy is bumped in the future:**
- scipy 1.15.0 adds full `sparse.linalg` compatibility with the new array API
- scipy 1.16.0 requires Python 3.11+ — do not bump to 1.16 while `requires-python = ">=3.10"`

---

## What NOT to Add

| Considered | Reason to Reject |
|------------|-----------------|
| Any IHW Python package | None exist on PyPI as of February 2026 |
| `scikit-allel` | ~50MB dep; `scipy.sparse.csr_array` already handles the sparse genotype matrix need |
| `plink2` Python bindings | No stable Python bindings exist; subprocess is the correct pattern |
| `polars` | Pandas deeply embedded; mixed-frame code would be a maintenance burden |
| New rpy2 version bump | Current `rpy2>=3.6.4` constraint is sufficient |
| `pyarrow` version bump | Already pinned `pyarrow>=14.0`; no new association feature needs a newer version |

---

## Sources

- statsmodels GLM `var_weights` documentation (verified 0.15.0):
  https://www.statsmodels.org/dev/generated/statsmodels.genmod.generalized_linear_model.GLM.html
- statsmodels `multipletests` — no weights parameter (verified 0.15.0):
  https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
- scipy `false_discovery_control` — no weights parameter (verified 1.17.0):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html
- scipy sparse docs — CSR recommendation, new array API (1.17.0):
  https://docs.scipy.org/doc/scipy/reference/sparse.html
- scipy 1.15.0 release notes — `sparse.linalg` array API compatibility:
  https://docs.scipy.org/doc/scipy-1.15.0/release/1.15.0-notes.html
- scipy 1.16.0 release notes — Python 3.11+ requirement:
  https://docs.scipy.org/doc/scipy/release/1.16.0-notes.html
- IHW Bioconductor (R only, no Python port, current release 3.20):
  https://bioconductor.org/packages/3.20/bioc/html/IHW.html
- IHW GitHub (nignatiadis/IHW — R only):
  https://github.com/nignatiadis/IHW
- statsmodels weighted GLM examples:
  https://www.statsmodels.org/dev/examples/notebooks/generated/glm_weights.html
