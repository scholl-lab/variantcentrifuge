# Domain Pitfalls: Rare Variant Association Testing (SKAT/SKAT-O/ACAT-O)

**Domain:** Statistical rare variant association testing in clinical genomics
**Focus:** Adding SKAT/SKAT-O, burden tests, ACAT-O, covariate adjustment, PCA integration to VariantCentrifuge
**Researched:** 2026-02-19
**Context:** GCKD cohort (5,125 samples, 22 GB VCF). Existing Fisher's exact pipeline. Adding dual-backend (R gold standard + Python fallback), compiled Davies method, covariates, PCA.

---

## Executive Summary

Implementing rare variant association tests correctly is deceptively hard. The statistical methods look simple in papers but have layers of edge cases: numerical precision failures that silently return wrong p-values, binary trait corrections that are absolutely required but often omitted, rpy2 integration that breaks in multi-threaded environments, and genotype matrix construction mistakes that invalidate entire analyses. Most pitfalls do not raise errors — they return plausible-looking but wrong results.

The three most dangerous failure modes for this milestone:

1. **Silent numerical failure:** Davies method returns `0.0` or `1.0` for extreme p-values without any error. Results look valid but are wrong.
2. **Binary trait type I error inflation:** Using continuous-trait SKAT on binary phenotypes with small/imbalanced samples produces inflated false positive rates — up to 5-fold above nominal level.
3. **Population stratification overcorrection:** Using too many PCs (>10) inflates type I error rates in burden and SKAT tests, counteracting the purpose of stratification adjustment.

---

## Critical Pitfalls

### Pitfall 1: Davies Method Silent Failure for Extreme P-Values

**What goes wrong:**

The Davies method for computing p-values from mixtures of chi-squared distributions can "falsely converge" and return incorrect results in the range 10^-6 to 10^-8 when using default accuracy settings (10^-6 accuracy, 10^4 max integration terms). Crucially, the function does not raise an error — it returns a plausible-looking p-value that is simply wrong.

Additionally, when the true p-value is extremely small, the method may return exactly `0.0` due to floating-point underflow. This looks like "very significant" but is uninformative — you cannot distinguish p=10^-20 from p=10^-300.

```python
# DANGEROUS: Default Davies accuracy — silently wrong for p < 1e-6
p_value = davies_method(Q_stat, eigenvalues)  # May return wrong value

# CORRECT: Higher accuracy settings with saddlepoint fallback
p_value = davies_method(Q_stat, eigenvalues,
                        accuracy=1e-9,
                        lim=1e6)            # More terms
if p_value == 0.0:
    # Fallback: saddlepoint approximation gives bounded estimate
    p_value = saddlepoint_approximation(Q_stat, eigenvalues)
```

**Why it happens:**

Davies algorithm inverts a characteristic function integral. For very small tail probabilities, the integral requires many terms and tight accuracy. The SKAT R package default (10^-6 accuracy, 10^4 terms) was calibrated for genome-wide thresholds, not for the 10^-8 to 10^-10 range where exome studies operate.

**Consequences:**

- Type I error inflation at stringent significance thresholds
- Gene-level p-values that appear significant but are numerical artifacts
- Analysis at genome-wide significance (5×10^-8) is unreliable with default settings

**Prevention:**

1. **Use stricter Davies settings:** accuracy=1e-9, lim=1e6 (from Liu et al. 2016 PMC4761292)
2. **Implement hybrid fallback:** Davies → saddlepoint approximation → Satterthwaite (last resort)
3. **Never use p=0.0 directly:** Log "p < X" where X is the precision floor; do not rank genes by exact 0.0 values
4. **Validate with R SKAT:** The R SKAT package v2.0+ uses the corrected hybrid approach — compare Python implementation against R results on the same data

**Warning signs:**

- p-values of exactly `0.0` appearing frequently
- P-value distributions showing unexpected spikes at 10^-6 or 10^-8
- Results disagree with R SKAT on the same input

**Implementation step:** Davies wrapper in Python backend — implement with correct defaults from day one.

**Severity:** CRITICAL — silent wrong p-values

**Sources:**
- [Liu et al. 2016 — efficient p-value calculation for SKAT (PMC4761292)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761292/) — PRIMARY source for hybrid Davies/saddlepoint approach
- [SKAT CRAN package documentation](https://cran.r-project.org/web/packages/SKAT/SKAT.pdf)

---

### Pitfall 2: Using Continuous-Trait SKAT for Binary Phenotypes Without Small-Sample Correction

**What goes wrong:**

SKAT was developed for quantitative traits. When applied to binary (case/control) traits with small samples, the asymptotic null distribution is inaccurate because the sparse genotype matrix makes the "large-sample" assumption invalid. The result is conservative type I error — empirical type I error is 5-fold below nominal. This means you lose substantial power.

More dangerously: **case-control imbalance causes the opposite effect** — inflated type I error. With very few cases (e.g., 50 cases, 5000 controls), standard SKAT produces false positives at 5-fold above nominal.

```python
# WRONG for binary traits, even if "it runs"
null_model = fit_null_model(y_binary, covariates, family="binomial")
p_value = skat(G, null_model, method="davies")  # Anti-conservative with imbalanced data

# CORRECT: Use SKATBinary equivalent with small-sample adjustment
null_model = fit_null_model(y_binary, covariates, family="binomial")
p_value = skat_binary(G, null_model,
                      method="Hybrid",          # Selects MA/QA/ER based on MAC
                      is_accurate_small_sample=True)
```

**Why it happens:**

The SKAT statistic's variance under the null is estimated from the asymptotic formula. For binary traits with sparse genotype data, the actual small-sample variance is much smaller than the asymptotic estimate, making the test either over- or under-conservative depending on the direction of imbalance.

**Consequences:**

- Conservative SKAT: miss real associations (power loss up to 80%)
- Inflated SKAT with case-control imbalance: false discoveries in top results
- Passing p-values to researchers who report associations that do not replicate

**Prevention:**

1. **Always use `SKATBinary` for binary traits, never `SKAT`:** The R package enforces this distinction — mirror it in Python
2. **Implement moment-matching (MA) adjustment:** Uses eigenvalue decomposition to estimate exact small-sample variance and kurtosis
3. **Check effective sample size before running:** Warn if n_cases < 200 or case:control ratio > 1:20
4. **Default to "Hybrid" method:** Automatically selects MA/QA/ER based on total minor allele count (MAC), number of carriers, and case-control imbalance

**Warning signs:**

- P-value histogram is sharply L-shaped (too conservative) or U-shaped (inflated)
- QQ-plot lambda > 1.1 under null simulation
- R SKAT result disagrees by >2 orders of magnitude from Python result on binary traits

**Implementation step:** Phase 1 of association framework — trait type detection must gate which null model and test function is called.

**Severity:** CRITICAL — produces wrong p-values for the primary use case (case/control)

**Sources:**
- [Lee et al. 2012 — SKAT-O with small-sample adjustment (PMC3415556)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
- [SKAT package — SKATBinary documentation](https://www.rdocumentation.org/packages/SKAT/versions/2.0.0/topics/SKATBinary)

---

### Pitfall 3: Population Stratification Overcorrection with Excessive PCs

**What goes wrong:**

Adding PCs as covariates corrects for population stratification — but adding too many PCs causes its own type I error inflation. Using 50 PCs inflates type I error regardless of the level of actual stratification, because the logistic regression model becomes overfit. The optimal number of PCs is usually 2-10, not 20-50.

An additional complication: rare variants are more geographically clustered than common variants (due to recency of mutations), so PCA computed on common variants may not adequately capture stratification for rare variant tests. SKAT tests are particularly sensitive to this because population-specific rare variants create spurious associations.

```python
# DANGEROUS: Too many PCs
covariates = pd.concat([age, sex, pcs[:50]], axis=1)  # Inflated type I error

# CORRECT: Conservative number of PCs
covariates = pd.concat([age, sex, pcs[:10]], axis=1)  # Tested in literature
# AND: Validate PC count using null simulations or genomic inflation lambda
```

**Why it happens:**

PCA on genotype matrices identifies population structure axes. The first few PCs capture continental ancestry; later PCs capture LD blocks and local structure. Including PCs that capture LD structure as covariates removes legitimate signal from the test statistic.

**Consequences:**

- Type I error inflation: false positives in fine-scale stratified cohorts
- Reduced power: PCs eat degrees of freedom
- Inconsistent results across cohorts if different PC counts used

**Prevention:**

1. **Compute genomic inflation factor (lambda GC) under null:** Using permuted phenotypes, compute SKAT p-values and measure lambda. Should be ~1.0 for correct number of PCs
2. **Start with 5 PCs, validate with simulation:** Increase only if lambda > 1.05 under null
3. **Use PCA on common variants (MAF > 0.01) for stratification correction**, not rare variants
4. **Document PC computation method:** LD pruning, MAF threshold, exclusion of high-LD regions (MHC, inversions) — differences here invalidate cross-cohort comparison
5. **Warn if requested PCs > 20** with a message explaining the inflation risk

**Warning signs:**

- QQ-plot inflation (lambda > 1.1) after adding more PCs
- Results change substantially when varying PC count from 5 to 10 to 20
- Genomic inflation factor increases (gets worse) as PC count increases

**Implementation step:** PCA module and covariate preparation stage — enforce PC count guidance and provide lambda diagnostic.

**Severity:** CRITICAL — systematic false positive inflation in stratified cohorts

**Sources:**
- [Genome-wide stratification impact on rare variant tests (PMC6283567)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283567/) — empirical evidence for overcorrection with 50 PCs
- [Controlling stratification in rare variant studies (PMC8463695)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463695/)

---

### Pitfall 4: rpy2 Thread Safety — R is Single-Threaded

**What goes wrong:**

R itself is not thread-safe. The rpy2 bridge inherits this limitation. Importing rpy2 in any non-main thread causes `ValueError: signal only works in main thread`. Calling R functions from multiple threads simultaneously corrupts the R interpreter state, causing segfaults or silent wrong results.

VariantCentrifuge uses `ThreadPoolExecutor` and `ProcessPoolExecutor` in `pipeline_core/runner.py`. If the SKAT R backend is called within a thread pool stage, the pipeline will crash or produce wrong results.

```python
# WRONG: R backend called from ThreadPoolExecutor worker thread
with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(run_skat_r, gene_data) for gene in genes]
    # Crashes with ValueError or segfault

# CORRECT: R backend only in main thread; use ProcessPoolExecutor for Python backend
# Option A: Serial R calls from main thread
for gene in genes:
    run_skat_r(gene_data)  # Main thread only

# Option B: ProcessPoolExecutor for Python backend (separate interpreter per process)
with ProcessPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(run_skat_python, gene_data) for gene in genes]
```

**Why it happens:**

R's signal handling is registered to the main thread at startup. rpy2 wraps R's C interface, inheriting all thread restrictions. Python's GIL does not help here because the unsafe code is inside the R C library, not Python.

**Consequences:**

- Crash on multi-gene parallel analysis (which is the common case)
- Segfaults that corrupt output files without error messages
- Silent wrong results if R state is corrupted but interpreter continues

**Prevention:**

1. **Establish architecture invariant:** R backend calls are ONLY from main thread or a dedicated R worker process
2. **Route R calls through a serializing queue:** Main thread owns R; parallel worker threads submit jobs to a queue and wait for results
3. **Use ProcessPoolExecutor for Python backend only** (separate interpreter = no shared R state)
4. **Test in multi-gene mode explicitly:** Run with 100+ genes and verify no crashes
5. **In PipelineStage declaration:** Mark any stage using R backend as `parallel_safe = False`

**Warning signs:**

- Intermittent crashes only when running with multiple genes
- `ValueError: signal only works in main thread` during initialization
- Segfaults with no Python traceback (crash inside R C library)

**Implementation step:** Architecture decision in Phase 1 — must be decided before writing any SKAT code.

**Severity:** CRITICAL — crashes parallel pipeline

**Sources:**
- [Run R in multi-threaded environment — BioStars discussion](https://www.biostars.org/p/174459/)
- [rpy2 thread safety — Streamlit issue](https://github.com/streamlit/streamlit/issues/2618)

---

## High Severity Pitfalls

### Pitfall 5: rpy2 Memory Leaks with Large R Objects in Loops

**What goes wrong:**

When calling R functions from Python in a gene-by-gene loop, R objects created during each call are not freed because Python does not know the size of R objects (they live in R's memory space). Python's garbage collector may not trigger often enough, causing transient memory explosions.

For VariantCentrifuge with 5,125 samples and thousands of genes, each SKAT call allocates an n×p kernel matrix. Without explicit cleanup, memory grows unboundedly across genes.

```python
# DANGEROUS: R objects accumulate without cleanup
for gene in genes:
    result = r_skat_function(G_matrix, null_model, weights)
    p_values[gene] = result  # R objects held alive via references

# CORRECT: Explicit R garbage collection per batch
import gc
import rpy2.robjects as ro

r_gc = ro.r('gc')  # R's garbage collector
for i, gene in enumerate(genes):
    result = r_skat_function(G_matrix, null_model, weights)
    p_values[gene] = float(result[0])  # Convert to Python primitive immediately
    del result                          # Release rpy2 reference
    if i % 100 == 0:
        r_gc()                          # Trigger R GC every 100 genes
        gc.collect()                    # Trigger Python GC too
```

**Why it happens:**

R uses NAMED-based reference tracking rather than Python's reference counting. rpy2 bridges the two, but the size of R objects is not visible to Python's GC. When large matrices (5125×50 = 256K floats per gene) are created in a loop without cleanup, peak memory can exceed available RAM.

**Consequences:**

- OOM killer terminates the pipeline after hours of computation
- Swap thrashing on HPC nodes slows analysis to unusable speeds
- Memory profiling shows monotonically increasing usage

**Prevention:**

1. **Convert R results to Python primitives immediately:** `float(r_result[0])` not `r_result[0]`
2. **Explicit `del` of R object references after each gene**
3. **Call `rpy2.robjects.r('gc()')` every N genes** (N=50-100 depending on cohort size)
4. **Process genes in batches:** Complete and clear each batch before starting next
5. **Profile memory on a 100-gene run** before declaring production-ready

**Warning signs:**

- Memory usage grows linearly with number of genes processed
- `htop` shows increasing RSS without plateau
- Jobs killed with signal 9 (OOM) on HPC

**Implementation step:** R backend wrapper — build cleanup into the gene iteration loop.

**Severity:** HIGH — OOM crash on large cohorts

**Sources:**
- [rpy2 memory management documentation](https://rpy2.github.io/doc/v3.2.x/html/rinterface-memorymanagement.html)
- [rpy2 performance and memory guidelines](https://rpy2.github.io/doc/latest/html/performances.html)

---

### Pitfall 6: Eigenvalue Computation Instability from Near-Singular Kernel Matrices

**What goes wrong:**

The SKAT test statistic Q has a null distribution that is a mixture of chi-squared distributions, where mixing weights are the eigenvalues of a matrix derived from the genotype covariance structure. Computing eigenvalues of a near-singular (low-rank) matrix produces numerical instability: negative eigenvalues that should be zero, eigenvalues differing by orders of magnitude, and non-convergence in iterative solvers.

This happens routinely in rare variant analysis because:
- Many genes have only 2-5 qualifying variants (very low-rank genotype matrix)
- High LD between variants makes the covariance matrix near-singular
- With all-rare-allele variants, most eigenvalues are effectively zero

```python
# PROBLEM: Eigenvalues with numerical noise
eigenvalues = np.linalg.eigvalsh(K_matrix)
# eigenvalues = [-1e-14, 2.3e-12, 0.45, 1.2, 3.7]  # Negative eigenvalues!

# CORRECT: Threshold small/negative eigenvalues
eigenvalues = np.linalg.eigvalsh(K_matrix)
eigenvalues = eigenvalues[eigenvalues > 1e-10 * eigenvalues.max()]  # Trim near-zero
# Use only positive eigenvalues for Davies method
```

**Why it happens:**

Floating-point arithmetic in eigendecomposition produces small negative values for theoretically zero eigenvalues. The Davies method requires nonnegative eigenvalues (kernel matrices are positive semidefinite by construction, but numerically may not be). Passing negative eigenvalues to Davies produces undefined behavior.

**Consequences:**

- Davies method returns NaN or garbage p-values
- Crash in eigenvalue computation for degenerate cases
- Different results across platforms (LAPACK implementations differ)

**Prevention:**

1. **Use `np.linalg.eigvalsh` not `np.linalg.eig`:** `eigvalsh` is specialized for symmetric matrices and is more numerically stable
2. **Threshold negative eigenvalues to zero:** `eigenvalues = np.maximum(eigenvalues, 0)` before Davies
3. **Use SVD as alternative for near-singular matrices:** More stable than eigendecomposition for rank-deficient cases
4. **Pre-check matrix rank:** If `np.linalg.matrix_rank(G, tol=1e-10) < 2`, skip SKAT and return NA (too few effective variants)
5. **Use `scipy.linalg.eigh` with `driver='evd'`** for full stability on small matrices

**Warning signs:**

- NaN p-values for specific genes
- Results that differ between platforms (Windows vs Linux) due to LAPACK differences
- Genes with 1-2 variants consistently failing

**Implementation step:** Eigenvalue computation in Python SKAT backend — validate inputs and handle degenerate cases.

**Severity:** HIGH — silent NaN results or platform-dependent wrong answers

**Sources:**
- [fastSKAT — eigenvalue approximation for large n (PMC4375394)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4375394/)
- [NumPy eigvalsh documentation](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eig.html)

---

### Pitfall 7: Genotype Matrix Construction — Multi-Allelic Sites and Missing Genotypes

**What goes wrong:**

Building the n_samples × n_variants dosage matrix (0, 1, 2) from VCF genotypes has multiple failure modes specific to rare variant analysis:

**Multi-allelic sites:** A site with REF=A, ALT=T,G has genotypes coded as 0 (ref), 1 (first alt T), 2 (second alt G). A genotype of `1/2` means heterozygous for two different alternate alleles. If you decode `1/2` as dosage=2, you overstate the alternate allele burden. The correct approach for SKAT is to treat multi-allelic sites as dosage=1 (carrier) or split into separate rows per ALT.

**Missing genotypes:** `./. ` means unknown, not reference. Mean-imputation (impute to 2×MAF) is standard for SKAT but changes the effective sample size calculation. Setting missing to 0 (reference) systematically biases the test toward the null.

**Phase (| vs /):** For SKAT, phase information is irrelevant — `0|1` and `0/1` are both dosage=1. Code that doesn't handle `|` as a separator silently treats phased genotypes as unreadable.

```python
# WRONG: Does not handle multi-allelic or missing
def gt_to_dosage(gt):
    return gt.count("1")  # "1/2" returns 2 (WRONG), "./." returns 0 (WRONG)

# CORRECT:
def gt_to_dosage(gt):
    if gt in (".", "./.", ".|."):
        return np.nan  # Missing — impute separately
    # Normalize separator
    alleles = gt.replace("|", "/").split("/")
    # Count non-reference alleles (any value != "0" or ".")
    return sum(1 for a in alleles if a not in ("0", "."))
    # For multi-allelic 1/2: returns 2 which is correct (2 non-ref alleles total)
    # But for SKAT burden collapse, may want: min(count, 1) for carrier-only
```

**Why it happens:**

VariantCentrifuge already handles this for Fisher's exact test, but SKAT needs a true dosage matrix rather than binary carrier calls. The existing `_gt_to_dosage` function in `gene_burden.py` handles `1/1`, `0/1`, `1|1`, `0|1` correctly but may not handle `1/2` (multi-allelic het) or missing values correctly for SKAT input.

**Consequences:**

- Incorrect dosage values for multi-allelic sites inflate burden statistics
- Missing genotypes coded as reference (dosage=0) bias test toward null
- Phase separator bugs produce zero-dosage for all phased genotypes

**Prevention:**

1. **Audit existing `_gt_to_dosage`** — verify it handles: `1/2`, `2/2`, `.|.`, `.`, `0|1`
2. **For SKAT, use mean imputation** for missing genotypes: `G[G.isna()] = 2 * MAF_of_variant`
3. **For multi-allelic sites in SKAT**: either (a) split into separate columns per ALT allele or (b) count all non-REF alleles (dosage=2 for `1/2`)
4. **Log the count of missing genotypes** per gene — if >10% missing, warn
5. **Unit test with synthetic VCF** covering all genotype formats

**Warning signs:**

- All-zero rows in dosage matrix (indicates phase-separator bug)
- Different carrier counts between burden analysis and SKAT
- Genes near multi-allelic sites showing unusual results

**Implementation step:** Genotype matrix builder — new function that extends existing `_gt_to_dosage`.

**Severity:** HIGH — incorrect input matrix invalidates all test results

**Sources:**
- [rvtests missing genotype mean imputation](https://github.com/zhanxw/rvtests)
- [Multi-allelic variant representation (PMC7288273)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7288273/)

---

### Pitfall 8: SKAT-O Rho Grid and Type I Error Control

**What goes wrong:**

SKAT-O optimizes over a mixing parameter rho (ρ) that interpolates between SKAT (ρ=0) and burden test (ρ=1). The original implementation uses 11 equally spaced ρ values from 0 to 1. The p-value of SKAT-O is computed by integrating over this grid — but the integration itself requires careful implementation to maintain proper type I error control.

Two specific bugs appear in re-implementations:

1. **Using `SKAT_O` instead of `optimal.adj`:** The "optimal.adj" method applies a small-sample adjustment to the combined statistic. The raw "SKAT_O" without adjustment has conservative type I error for small samples.

2. **Wrong variance estimation in the combination step:** When combining SKAT and burden p-values for the optimal ρ, the variance-covariance matrix of the component statistics must be estimated correctly. Errors here cause the rho-specific p-value transformation to be wrong, breaking the uniformity guarantee.

```python
# WRONG: Naive rho grid without small-sample adjustment
rho_values = np.linspace(0, 1, 11)
pvals_per_rho = [skat_with_rho(G, residuals, rho) for rho in rho_values]
p_skatO = min_pval_combination(pvals_per_rho)  # WRONG — no correlation structure

# CORRECT: Account for correlation between rho-specific statistics
# The key step is computing the covariance of (T_rho1, T_rho2, ...) statistics
# This requires the kernel matrix eigenvalues at each rho value
```

**Why it happens:**

SKAT-O is computationally intensive and the math is non-trivial. Re-implementations often simplify the combination step incorrectly. The R SKAT package source code is the ground truth — any Python reimplementation must produce identical results on test cases.

**Consequences:**

- SKAT-O p-values that are uniformly too small or too large
- Loss of power advantage over SKAT (SKAT-O should be more powerful)
- Inconsistent comparison with published results using R SKAT

**Prevention:**

1. **Use R backend as ground truth for SKAT-O** — do not reimplement from scratch for the Python fallback
2. **For Python fallback, implement SKAT only** (not SKAT-O) initially; add SKAT-O only after extensive validation
3. **Validate on simulated null data:** 10,000 permutations should give uniform p-value distribution
4. **Cross-validate with R SKAT** on 50+ real genes: Python result should match R within 10% for p > 0.001

**Warning signs:**

- Python SKAT-O gives systematically different results from R SKAT-O
- P-value histogram shows non-uniform distribution under null
- Power advantage of SKAT-O over SKAT disappears

**Implementation step:** Python SKAT-O backend — validate before deploying.

**Severity:** HIGH — wrong test statistic for the omnibus test

**Sources:**
- [SKAT-O optimal unified approach (PMC3415556)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
- [SKAT R package GitHub — leelabsg/SKAT](https://github.com/leelabsg/SKAT)

---

### Pitfall 9: ctypes C Compilation Failures at Install Time

**What goes wrong:**

Compiling `qfc.c` (Davies' method C implementation) at install time fails on:

- **Windows:** No C compiler in PATH by default; MSVC requires separate installation; MinGW has inconsistent ABI
- **macOS ARM64 (Apple Silicon):** ctypes on ARM64 has known bugs with variadic function calls (CPython issue #42880); needs explicit architecture flag `-arch arm64`
- **HPC nodes:** Often lack write permission to `/tmp` or block `gcc` during job execution; compilation must happen at install time, not at runtime
- **conda environments:** Package-provided GCC may conflict with system GCC

```python
# DANGEROUS: Compile at first use (runtime compilation on HPC)
try:
    lib = ctypes.CDLL("/tmp/qfc.so")
except OSError:
    os.system("gcc -shared -o /tmp/qfc.so qfc.c")  # Fails on HPC!
    lib = ctypes.CDLL("/tmp/qfc.so")

# CORRECT: Check at module import time with informative error
import importlib.resources
_QFC_LIB = None

def _load_qfc():
    global _QFC_LIB
    try:
        lib_path = _find_compiled_qfc()  # Look in package data first
        _QFC_LIB = ctypes.CDLL(lib_path)
    except OSError:
        _QFC_LIB = None
        logger.warning("qfc.c not compiled — Davies method unavailable; using saddlepoint")
```

**Why it happens:**

C compilation is platform-specific and requires a toolchain. Users may not have the required tools, or compilation may fail silently. The qfc.c file from the CompQuadForm R package is the standard Davies implementation used by SKAT.

**Consequences:**

- Silent fallback to less accurate approximation (user doesn't know Davies failed)
- Installation fails on user machines, blocking all SKAT functionality
- Different results between users who successfully compiled and those who didn't

**Prevention:**

1. **Make Davies method optional, never required:** Saddlepoint fallback must always work
2. **Compile as part of package setup.py/build** with clear error on failure
3. **Provide pre-compiled wheels for common platforms** (Linux x86_64, macOS arm64/x86_64, Windows x64)
4. **Add `check_dependencies` function** that tests compilation at startup and logs the result
5. **Never compile at runtime** — always at install time or at package import in development mode
6. **Test CI on all three platforms:** GitHub Actions matrix with ubuntu, macos, windows

**Warning signs:**

- `ImportError` or `OSError` when loading qfc.so
- Users report different p-values across machines (some have Davies, some don't)
- CI passes but user machine fails to compile

**Implementation step:** Package build system — must be designed before first release.

**Severity:** HIGH — breaks installation on common platforms

**Sources:**
- [CPython ARM64 ctypes variadic bug (issue #42880)](https://bugs.python.org/issue42880)
- [Python ctypes cross-platform package sample](https://github.com/joerick/python-ctypes-package-sample)

---

### Pitfall 10: Multiple Testing Correction Applied to the Wrong Set

**What goes wrong:**

When reporting SKAT + burden + ACAT-O p-values for the same genes, applying multiple testing correction to each test separately (N genes per test) is incorrect — you are running 3 tests per gene and must account for the additional multiplicity. However, combining all three tests into one correction (3N hypotheses) is also wrong because the tests are correlated (all test the same gene).

Additionally, the Benjamini-Hochberg FDR procedure assumes continuous p-values with uniform null distribution. For discrete tests (Fisher's exact, burden with few carriers), this assumption fails and BH can be anti-conservative.

```python
# WRONG: Independent correction per test
skat_fdr = multipletests(skat_pvals, method='fdr_bh')[1]
burden_fdr = multipletests(burden_pvals, method='fdr_bh')[1]
# This pretends SKAT and burden are testing different genes

# WRONG: Naive combination of all tests
all_pvals = np.concatenate([skat_pvals, burden_pvals, acato_pvals])
all_fdr = multipletests(all_pvals, method='fdr_bh')[1]
# This penalizes each test 3x by inflating denominator

# CORRECT for ACAT-O approach: combine p-values first, then correct once
acato_pvals = [acat_combine([skat_p, burden_p, acatv_p]) for gene in genes]
gene_level_fdr = multipletests(acato_pvals, method='fdr_bh')[1]
```

**Why it happens:**

Multiple testing correction was designed for independent tests. Gene-based rare variant tests running SKAT + burden on the same gene are neither independent nor identically distributed. The ACAT-O framework exists precisely to address this — combine within-gene tests using Cauchy combination, then correct across genes.

**Consequences:**

- Overly conservative results (too few discoveries)
- Inconsistent results reported for SKAT vs burden from same analysis
- Different FDR cutoffs needed for different test types confuses users

**Prevention:**

1. **Use ACAT-O as the primary test** — it combines SKAT, burden, and variant-level results in one framework, then apply single FDR correction across genes
2. **Report raw p-values for all tests** alongside the combined ACAT-O p-value — researchers can evaluate separately
3. **Document which p-value was corrected** explicitly in output headers
4. **For exploratory analysis, use FDR (BH)** for hypothesis generation; for confirmatory analysis, use Bonferroni
5. **Warn when both SKAT and burden p-values are available** without ACAT-O combination

**Warning signs:**

- SKAT FDR and burden FDR give different significant genes (expected sometimes, but should be investigated)
- User reports "my gene was significant for SKAT but not after correction"
- Corrected p-values that are lower than raw p-values (mathematical error)

**Implementation step:** Results reporting stage — multiple testing correction strategy must be decided in design.

**Severity:** HIGH — incorrect statistical inference

**Sources:**
- [ACAT p-value combination for rare variants (PMC6407498)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/)
- [Leveraging gene-level prediction for rare variant testing (PMC8872452)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8872452/)

---

## Moderate Pitfalls

### Pitfall 11: R_HOME Detection Failures in Conda and HPC Environments

**What goes wrong:**

rpy2 requires `R_HOME` to be set to the R installation directory. In conda environments, R may be installed in a non-standard location. On HPC, module systems may change `PATH` but not `R_HOME`. The failure mode is an unhelpful `ImportError` or `RRuntimeError` at startup rather than at the point of first R call.

Common failure: user has both system R (/usr/bin/R) and conda R (/opt/conda/bin/R). rpy2 finds system R via R_HOME but conda SKAT package is installed in conda R. The mismatch causes `library(SKAT)` to fail with "package not found."

**Prevention:**

1. **Set R_HOME explicitly** in the VariantCentrifuge launcher based on which R is in PATH
2. **Validate SKAT is loadable before any analysis:** `r('library(SKAT)')` in a startup check function
3. **Provide clear error message** if SKAT is not found: "Install SKAT R package in the R found at [path]"
4. **Document conda setup steps** explicitly: `conda install -c bioconda r-skat` vs `R -e 'install.packages("SKAT")'`
5. **Fall back to Python backend automatically** if R/SKAT is unavailable (with warning)

**Warning signs:**

- Works in development environment but fails on CI or user machine
- `RRuntimeError: Error in library(SKAT)` at startup
- rpy2 imports successfully but R is wrong version

**Implementation step:** R backend initialization — check at import time with informative failure.

**Severity:** MEDIUM — blocks R backend but Python fallback handles it

**Sources:**
- [rpy2 + conda/mamba setup guide (2025)](https://sychen9584.github.io/posts/2025/06/rpy2/)
- [rpy2 RRuntimeError known issues](https://github.com/rpy2/rpy2/issues/911)

---

### Pitfall 12: Covariate Multicollinearity Causing Model Fitting Failures

**What goes wrong:**

PCs 1-10 may be highly correlated with categorical covariates like ancestry group or geographic region. Including both PCs and ancestry categories in the covariate matrix causes multicollinearity, making the null model matrix singular and the logistic regression unfit.

More subtly: if all cases have missing values for one covariate and all controls have values, the covariate perfectly predicts case status — the model is unidentifiable.

```python
# PROBLEM: PC1 correlates r=0.99 with "ancestry" categorical variable
covariates = [age, sex, pc1, pc2, ..., pc10, ancestry]  # Singular design matrix

# DETECT:
cond_number = np.linalg.cond(covariate_matrix)
if cond_number > 1e8:
    logger.warning(f"Covariate matrix near-singular (condition number: {cond_number:.2e})")
```

**Prevention:**

1. **Check condition number** of the covariate matrix before fitting null model
2. **Never include both PCs and ancestry categories** — use one or the other
3. **Detect perfect separation** (categorical variable perfectly predicts outcome) before fitting
4. **Log covariate missingness by case/control status** — unequal missingness is a red flag
5. **Handle factor encoding explicitly:** Dummy-encode categoricals with drop_first=True to avoid perfect collinearity

**Warning signs:**

- Null model fitting fails to converge
- Extreme coefficient estimates in logistic regression (>100)
- `np.linalg.LinAlgError: Singular matrix`

**Implementation step:** Covariate preparation module.

**Severity:** MEDIUM — causes analysis failure with cryptic error

---

### Pitfall 13: Sample Mismatch Between Genotype Matrix and Covariate Data

**What goes wrong:**

Genotype data comes from VCF (sample order determined by VCF header). Covariate data comes from a separate file with arbitrary row order. If sample IDs are matched incorrectly, each sample's covariates are assigned to the wrong genotype row, corrupting the entire analysis.

This is particularly dangerous because: (1) the analysis runs without error, (2) the p-values look plausible, and (3) the error is invisible without an explicit validation check.

**Prevention:**

1. **Always reindex covariate DataFrame** to VCF sample order explicitly
2. **Assert no NaN after reindex** — if a sample in VCF has no covariate row, that is an error
3. **Log sample counts:** "Matched X of Y VCF samples to covariate file"
4. **Include sample ID as an explicit check column** in the covariate file validation
5. **Test with shuffled covariate order** — results must be identical to correctly-ordered version

**Warning signs:**

- Different results when covariate file rows are reordered
- Number of analyzed samples differs from expected

**Implementation step:** Data loading stage — validation before any analysis.

**Severity:** MEDIUM — silent analysis corruption

---

### Pitfall 14: Small Gene Edge Cases — Zero Variants, Monomorphic Variants, MAC=1

**What goes wrong:**

Rare variant analysis has three common degenerate cases that must be handled explicitly:

**Zero qualifying variants in gene:** SKAT cannot run. Return NA, not 0 or 1.

**All variants monomorphic in analyzed samples:** After applying MAF filter and extracting the gene cohort, all variants may be fixed (everyone has the same genotype). The genotype matrix is all-zero — SKAT returns a test statistic of 0 and p=1.0 trivially, which is technically correct but may be misleading.

**Singletons (MAC=1):** A variant carried by exactly one individual. For binary traits, there are only 2 possible configurations (singleton in case or singleton in control). The SKAT contribution from a singleton is essentially a single Bernoulli trial, which has very limited statistical power. Including many singletons can dilute SKAT power.

```python
# Check before running SKAT
if G_matrix.shape[1] == 0:
    return {"gene": gene, "p_skat": np.nan, "n_variants": 0, "note": "no_qualifying_variants"}

if G_matrix.sum() == 0:
    return {"gene": gene, "p_skat": 1.0, "n_variants": G_matrix.shape[1], "note": "all_monomorphic"}

mac_per_variant = G_matrix.sum(axis=0)
n_singletons = (mac_per_variant == 1).sum()
if n_singletons > G_matrix.shape[1] * 0.8:
    logger.warning(f"Gene {gene}: {n_singletons}/{G_matrix.shape[1]} variants are singletons")
```

**Prevention:**

1. **Return NA (not 1.0 or error) for genes with no qualifying variants** — downstream correction methods need to know vs. skip
2. **Log monomorphic and singleton counts** per gene for QC purposes
3. **Provide option to exclude singletons** from SKAT analysis (they dilute power)
4. **Test all three edge cases** with explicit unit tests

**Warning signs:**

- Many genes with p=1.0 exactly (monomorphic)
- Crash on genes with 0 variants (index error in matrix operations)

**Implementation step:** Gene-level SKAT dispatcher — pre-flight checks before matrix construction.

**Severity:** MEDIUM — crashes or misleading results for common gene classes

**Sources:**
- [How singletons are treated in SKAT (SKAT Google Group)](https://groups.google.com/g/skat_slee/c/amQZMDTk7CE)
- [SKAT — handling monomorphic and high-missingness SNPs](https://cran.r-project.org/web/packages/SKAT/SKAT.pdf)

---

### Pitfall 15: Liu Moment-Matching Approximation Anti-Conservation

**What goes wrong:**

The Liu moment-matching approximation (non-central chi-squared) is used as a fallback when Davies method fails. However, it is systematically anti-conservative for small significance levels (p < 10^-4). This means using Liu as the primary method for extreme p-values inflates type I error — the exact opposite of what you want for rare variant association at stringent thresholds.

The Liu approximation is appropriate as a sanity check or for p-values near 0.05, but should never be the final answer for results below 10^-3.

**Prevention:**

1. **Fallback order:** Davies → saddlepoint approximation → Liu (last resort only)
2. **When Liu is used, flag in output:** Add a column `p_method` = "davies"/"saddlepoint"/"liu" so downstream analysis knows confidence level
3. **Saddlepoint approximation is strictly better than Liu** for small tails — implement it first

**Warning signs:**

- P-value histogram shows too many extreme values under null
- Results with Liu method systematically smaller than Davies for same data

**Implementation step:** P-value computation backend — approximation hierarchy.

**Severity:** MEDIUM — inflated type I error for exact p-values

**Sources:**
- [Satterthwaite vs Davies type I error comparison (PMC4761292)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761292/)

---

### Pitfall 16: Cryptic Relatedness in SKAT Null Model

**What goes wrong:**

SKAT's standard null model assumes independent samples. The GCKD cohort (kidney disease) may include related individuals (family members, duplicate samples). Including related samples in standard SKAT violates the independence assumption, inflating type I error.

For the GCKD cohort specifically: check if kinship matrix or identity-by-descent (IBD) estimates are available. If any pairs have kinship > 0.05, they are considered related and standard SKAT is invalid.

**Prevention:**

1. **Compute pairwise IBD/kinship** using PLINK or similar before running SKAT
2. **If relatedness exists, use SKAT with mixed model null** (emmax-style) or GENESIS package's SKAT implementation
3. **Alternatively, prune related samples** (keep one per related pair) — simpler but loses data
4. **Document relatedness assumption** in analysis report

**Warning signs:**

- Known to have family structure in the cohort
- Inflated genomic inflation factor (lambda > 1.0) even after PC correction

**Implementation step:** Data preparation — must be assessed before analysis, not during.

**Severity:** MEDIUM — systematic type I error inflation for related cohorts

**Sources:**
- [GENESIS package for SKAT with relatedness](https://rdrr.io/bioc/GENESIS/man/assocTestAggregate.html)

---

## Low Severity Pitfalls

### Pitfall 17: Variant Weighting Default Beta(1,25) Not Always Optimal

**What goes wrong:**

SKAT by default uses a beta(1,25) weighting function that up-weights very rare variants (MAF near 0) and down-weights common variants. This is appropriate when you expect rarer variants to have larger effects (the common assumption in rare disease genetics). However, for complex traits where effect size does not strongly correlate with rarity, beta(1,1) (equal weights) may be more powerful.

**Prevention:**

1. **Report results for both beta(1,25) and beta(1,1)** weighting schemes
2. **Allow weight scheme to be configured** in the analysis parameters
3. **Document which weights were used** in the output

**Severity:** LOW — affects power but not correctness

---

### Pitfall 18: ACAT-O Requires Uniform Null Distribution of Input P-Values

**What goes wrong:**

ACAT-O combines p-values using the Cauchy distribution, which requires that input p-values are uniformly distributed under the null. If the input SKAT p-values have been discretized (e.g., from a permutation test with 1000 permutations), the distribution is not uniform and ACAT-O is invalid.

**Prevention:**

1. **Use analytical p-values** (from Davies/saddlepoint) as ACAT-O inputs, not permutation p-values
2. **Validate input p-values are not all 0 or 1** before ACAT-O combination
3. **If any input test returns NA, use only available tests** in ACAT-O combination

**Severity:** LOW — only relevant if using permutation-based p-values

**Sources:**
- [ACAT methodology and assumptions (PMC6407498)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/)

---

### Pitfall 19: Validation Against Wrong Reference

**What goes wrong:**

When validating the Python SKAT implementation against R SKAT, small differences in:
- MAF computation (from VCF allele count vs from genotype matrix)
- Missing genotype handling (exclude vs mean-impute)
- Weight precision

...can produce p-value differences of up to 2x for the same gene. Treating any difference as a bug wastes time; treating real bugs as "acceptable tolerance" misses errors.

**Prevention:**

1. **Establish explicit tolerance:** p-value difference < 10% for p > 0.001 is acceptable; absolute difference in log10(p) < 0.05 is the standard
2. **Use identical input processing** — same MAF computation, same missing handling
3. **Test on the SKAT package's example dataset** (`SKAT.haplotypes`) which has known correct answers
4. **For p < 0.001, expect larger differences** between methods — this is normal numerical behavior

**Severity:** LOW — validation process pitfall, not a correctness pitfall in production

**Sources:**
- [SKAT.haplotypes reference dataset](https://cran.r-project.org/web/packages/SKAT/SKAT.pdf)

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| Dual backend architecture | R thread safety (#4) | Enforce main-thread-only for R; ProcessPool for Python |
| Null model fitting for binary traits | Binary trait correction (#2) | Use SKATBinary variant from day one |
| Genotype matrix extraction | Multi-allelic and missing (#7) | Extend existing `_gt_to_dosage` with explicit tests |
| PCA integration | PC count overcorrection (#3) | Default to 10 PCs max; lambda GC diagnostic |
| Davies method compilation | C compilation failures (#9) | Optional compile; saddlepoint always available |
| Eigenvalue computation | Near-singular matrix (#6) | Threshold and use `eigvalsh`, not `eig` |
| Multiple testing | Wrong correction set (#10) | ACAT-O combination then single FDR |
| R memory management | Memory leaks in loops (#5) | Explicit R GC every 100 genes |
| SKAT-O implementation | Rho grid correctness (#8) | Use R backend; validate Python with extensive null sim |
| Sample/covariate joining | Sample mismatch (#13) | Explicit reindex with assertion |
| Gene edge cases | Zero/singleton variants (#14) | Pre-flight checks before matrix construction |
| Validation against R | Wrong tolerance (#19) | Document acceptable tolerance before testing |

---

## Summary: Risk Matrix

| Pitfall | Severity | Detection Difficulty | Implementation Phase |
|---------|----------|---------------------|---------------------|
| Davies silent failure (#1) | CRITICAL | Hard — no error raised | P-value computation |
| Binary trait type I error (#2) | CRITICAL | Medium — requires simulation | Null model design |
| PC overcorrection (#3) | CRITICAL | Medium — QQ-plot needed | PCA module |
| R thread safety (#4) | CRITICAL | Easy — crashes immediately | Architecture Phase 1 |
| R memory leaks (#5) | HIGH | Medium — slow accumulation | R backend iteration |
| Eigenvalue instability (#6) | HIGH | Medium — NaN for some genes | SKAT computation |
| Genotype matrix construction (#7) | HIGH | Hard — plausible wrong values | Input preparation |
| SKAT-O rho correctness (#8) | HIGH | Medium — requires validation | SKAT-O backend |
| C compilation failures (#9) | HIGH | Easy — install-time failure | Package build |
| Multiple testing strategy (#10) | HIGH | Medium — review needed | Results reporting |
| R_HOME detection (#11) | MEDIUM | Easy — clear error message | R initialization |
| Covariate multicollinearity (#12) | MEDIUM | Easy — condition number check | Covariate prep |
| Sample mismatch (#13) | MEDIUM | Hard — silent corruption | Data loading |
| Edge case genes (#14) | MEDIUM | Easy — explicit pre-checks | Gene dispatcher |
| Liu anti-conservation (#15) | MEDIUM | Medium — requires simulation | P-value backend |
| Cryptic relatedness (#16) | MEDIUM | Medium — requires IBD analysis | Study design |

---

## Key References

**Numerical Precision:**
- [Liu et al. 2016 — hybrid Davies/saddlepoint for SKAT (PMC4761292)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761292/)
- [SKAT CRAN package — CompQuadForm Davies implementation](https://cran.r-project.org/web/packages/CompQuadForm/index.html)

**Binary Traits and Small Samples:**
- [Lee et al. 2012 — SKAT-O with small-sample adjustment (PMC3415556)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
- [SKATBinary documentation](https://www.rdocumentation.org/packages/SKAT/versions/2.0.0/topics/SKATBinary)

**Population Stratification:**
- [Fine-scale stratification impact on rare variant tests (PMC6283567)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283567/)
- [Controlling stratification in rare variant studies (PMC8463695)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463695/)

**rpy2 Integration:**
- [rpy2 performance and memory documentation](https://rpy2.github.io/doc/latest/html/performances.html)
- [rpy2 memory management](https://rpy2.github.io/doc/v3.2.x/html/rinterface-memorymanagement.html)
- [rpy2 thread safety — Streamlit issue](https://github.com/streamlit/streamlit/issues/2618)

**ACAT-O:**
- [ACAT: Cauchy combination for rare variant p-values (PMC6407498)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/)

**Genotype Matrix:**
- [Multi-allelic variant association analysis (PMC7288273)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7288273/)
- [rvtests — mean imputation for missing genotypes](https://github.com/zhanxw/rvtests)

**SKAT Implementation:**
- [SKAT R package GitHub — leelabsg/SKAT](https://github.com/leelabsg/SKAT)
- [SKAT original paper (PMC3135811)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3135811/)

---

## Research Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Davies numerical precision | HIGH | Primary literature (Liu et al. 2016) verified with official SKAT docs |
| Binary trait correction | HIGH | Original SKAT-O paper (Lee et al. 2012) + SKAT package docs |
| PC overcorrection | HIGH | Empirical study (PMC6283567) with specific PC counts tested |
| rpy2 thread safety | HIGH | Multiple independent sources confirm; documented in rpy2 code |
| rpy2 memory leaks | HIGH | Official rpy2 performance documentation |
| Genotype matrix construction | MEDIUM | Cross-referenced rvtests and VCF specifications |
| Eigenvalue stability | MEDIUM | NumPy documentation + general numerical linear algebra |
| ctypes compilation | MEDIUM | CPython issue tracker + cross-platform package examples |
| SKAT-O rho correctness | MEDIUM | Based on SKAT paper + general knowledge; specific bugs LOW confidence |
| Multiple testing strategy | HIGH | ACAT paper + FDR literature well-established |

**Overall confidence:** HIGH for critical pitfalls (numerical, thread safety, binary traits); MEDIUM for implementation-specific details (SKAT-O rho bugs, eigenvalue thresholds).
