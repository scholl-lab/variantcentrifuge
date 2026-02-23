# Variantcentrifuge Inheritance Analysis: Performance Deep-Dive & Optimization Report

## Executive Summary

The inheritance analysis in variantcentrifuge v0.14.3 takes **~10-12 hours** for the full AGDE cohort (294,776 variants × 1,452 samples × 18,629 genes) on 16 cores. The v0.14.3 `ThreadPoolExecutor` fix solved the OOM problem but provides **zero actual parallelism** because the hot path is 100% GIL-bound Python code. With the optimizations proposed here, the runtime could drop to **5-30 minutes**.

---

## 1. Measured Timings

| Run | Samples | Variants | Genes | Pass 1 | Pass 2 | Pass 3 | Total |
|-----|---------|----------|-------|--------|--------|--------|-------|
| Research broad | 603 | 140,570 | 16,580 | 62s | 2,206s (37 min) | 207s | 2,484s (41 min) |
| Full broad | 1,452 | 294,776 | 18,629 | 228s | ~36,000s est (10h) | ~800s est | ~37,000s est (10.3h) |

**Scaling**: 2.4× more samples → **16× longer** Pass 2. This super-linear scaling confirms an algorithmic problem, not just more data.

---

## 2. Architecture & Call Graph

```
analyze_inheritance_parallel()                          [parallel_analyzer.py]
  │
  ├── Pass 1: vectorized_deduce_patterns()              [vectorized_deducer.py]
  │     ├── _build_genotype_matrix()                    ← Builds int8 (294k × 1452) matrix
  │     │     └── encode_genotypes() × 1,452 samples    ← 228s, done ONCE
  │     └── per-sample pattern checks (numpy)           ← 42s, vectorized, fast
  │
  ├── Pass 2: ThreadPoolExecutor(16 workers)            [parallel_analyzer.py]  ★ BOTTLENECK
  │     └── per gene (×18,629): _process_gene_group()
  │           └── analyze_gene_for_compound_het_vectorized()  [comp_het_vectorized.py]
  │                 ├── encode_genotypes() × 1,452 samples    ← REDUNDANT (already done in Pass 1)
  │                 ├── per-sample loop (×1,452):             ← GIL-bound, not parallelized
  │                 │     ├── het detection (tiny numpy)
  │                 │     ├── get_parents() (dict lookup)
  │                 │     └── find_potential_partners_vectorized() (Python loop)
  │                 └── create_variant_key_fast() × M         ← df.iloc[] per variant, expensive
  │
  └── Pass 3: _finalize_inheritance_patterns()          [analyzer.py]
        └── per variant (×294,776):                     ← Pure Python loop
              ├── calculate_segregation_score()          ← loops over 1,452 samples per pattern
              ├── prioritize_patterns()                 ← dict lookups
              ├── pd.Series() construction              ← ~50KB allocation per variant
              ├── create_inheritance_details()           ← loops over samples again
              └── json.dumps()                          ← per variant
```

---

## 3. Where Time Is Spent

| Component | Est. Time | % Total | GIL Status | Nature |
|-----------|-----------|---------|------------|--------|
| **Pass 2: encode_genotypes() per gene** | 4-5h | ~45% | **Held** | Redundant re-encoding of all 27M sample×gene pairs |
| **Pass 2: per-sample inner loop** | 3-4h | ~35% | **Held** | Python loop over 1,452 samples per gene |
| **Pass 3: segregation + finalize** | ~500-800s | ~7% | **Held** | Python loop, pd.Series construction per row |
| **Pass 1: build genotype matrix** | ~185s | ~2% | **Held** | encode_genotypes() ×1,452 (done once, OK) |
| **Pass 1: pattern checks** | ~42s | ~0.4% | Released | NumPy vectorized, fast |
| **Pass 2: apply comp het results** | ~200-400s | ~4% | **Held** | Nested Python loops |

**Key finding**: >95% of runtime is GIL-bound Python code. ThreadPoolExecutor's 16 threads achieve **effectively single-threaded** execution.

---

## 4. Root Causes

### 4.1 Redundant Genotype Encoding (Pass 2 repeats Pass 1)

Pass 1 builds a complete genotype matrix via `_build_genotype_matrix()`:
```python
gt_matrix = np.full((294776, 1452), -1, dtype=np.int8)  # 409 MB, built once
for i, sample_id in enumerate(sample_list):
    gt_matrix[:, i] = encode_genotypes(df[sample_id])   # 1,452 calls total
```

But this matrix is **discarded**. Pass 2 rebuilds it per-gene:
```python
# Inside analyze_gene_for_compound_het_vectorized(), called 18,629 times:
genotype_matrix = {}
for sample_id in sample_list:                            # 1,452 iterations
    genotype_matrix[sample_id] = encode_genotypes(gene_df_unique[sample_id])
```

Total `encode_genotypes()` calls: **1,452 × 18,629 = 27,009,108** (vs 1,452 in Pass 1).

Each call creates a pandas Series, converts Categorical to object, calls `fillna`, does string comparisons — for typically **2-50 values**. The overhead of Series creation dominates actual computation by 100×+.

### 4.2 ThreadPoolExecutor Cannot Parallelize GIL-Bound Code

The v0.14.3 fix replaced `ProcessPoolExecutor` with `ThreadPoolExecutor` to solve OOM. This fixed memory but eliminated any chance of parallelism. The compound het inner loop is:
- `encode_genotypes()`: pandas Series operations → GIL held
- `for sample_id in sample_list:` → Python loop → GIL held
- `get_parents()`, `is_affected()`: dict lookups → GIL held
- `find_potential_partners_vectorized()`: Python loop over `het_indices` → GIL held
- `create_variant_key_fast()`: `df.iloc[idx]` + f-string → GIL held

The only NumPy operations that release the GIL are `sample_genotypes == 1` and `np.where(het_mask)` on arrays of 2-50 elements — nanoseconds of work.

**Result**: 16 threads contend for the GIL, adding overhead with zero parallel gain. Likely **slower** than single-threaded.

### 4.3 Per-Variant Python Loop in Pass 3

Pass 3 iterates 294,776 times in pure Python. Per variant:
- `pd.Series()` construction from 1,452 sample columns: ~50 KB allocation
- `calculate_segregation_score()`: loops over samples per pattern
- `create_inheritance_details()`: loops over samples again
- `json.dumps()`: serialization per variant

Total transient allocations: ~14 GB churning through Python's garbage collector.

### 4.4 Super-Linear Scaling Mechanism

Scaling from 603→1,452 samples (2.4×) causes 16× slowdown because:
1. `encode_genotypes()` overhead per gene grows linearly in samples (2.4×)
2. More samples → more heterozygous variants per gene → more partner pairs → O(het²) pair enumeration
3. GIL contention worsens with more work per gene
4. Pass 3 segregation checking: 1,452 vs 603 samples per variant

At **10,000 samples**, estimated Pass 2 time: **days to weeks**.

---

## 5. Proposed Optimizations (Ranked by Impact)

### Optimization 1: Reuse Pass 1 Genotype Matrix in Pass 2

**Impact: ~10× speedup for Pass 2 (eliminates 27M redundant encode_genotypes calls)**

The genotype matrix built in Pass 1 should be passed to Pass 2. For each gene, slice the pre-computed matrix by the gene's variant indices instead of re-encoding from DataFrame columns.

```python
# Pass 1: Build matrix and KEEP it
gt_matrix, sample_to_idx = _build_genotype_matrix(df, sample_list)

# Pass 2: Slice per gene (zero-copy numpy view if contiguous)
for gene, gene_df in genes_with_multiple_variants:
    gene_row_indices = gene_df.index.values
    gene_gt = gt_matrix[gene_row_indices, :]  # int8 subarray, ~16×1452 bytes
    # Pass gene_gt and sample_to_idx to worker instead of gene_df
```

This eliminates:
- 27 million `encode_genotypes()` calls
- All pandas Series creation in Pass 2
- All Categorical→object dtype conversions
- All `fillna` and string comparison overhead

### Optimization 2: Vectorize Het Detection Across All Samples

**Impact: ~50× speedup for het identification step**

Replace the per-sample Python loop with a single matrix operation:

```python
# Instead of:
for sample_id in sample_list:                    # 1,452 Python iterations
    sample_gt = genotype_matrix[sample_id]
    het_mask = sample_gt == 1
    het_indices = np.where(het_mask)[0]
    if len(het_indices) < 2: continue

# Use:
het_mask = (gene_gt == 1)                        # (n_vars, 1452) bool, one numpy op
het_counts = het_mask.sum(axis=0)                # (1452,) int, one numpy op
candidate_samples = np.where(het_counts >= 2)[0] # typically small subset
# Only loop over candidate_samples (e.g., 50-200 instead of 1,452)
```

This reduces the Python loop from 1,452 iterations to only the few samples that actually have ≥2 het variants in the gene.

### Optimization 3: Sparse Matrix Multiply for Gene-Level Candidate Finding

**Impact: Replaces the entire per-gene loop for candidate identification**

Instead of iterating 18,629 genes and checking each sample, use one sparse matrix multiply:

```python
from scipy.sparse import csr_matrix

# Build gene→variant mapping matrix once (n_genes × n_variants, sparse)
gene_variant_map = csr_matrix(
    (np.ones(n_variants), (gene_ids, np.arange(n_variants))),
    shape=(n_genes, n_variants)
)

# het_mask: (n_variants, n_samples) bool from gt_matrix == 1
# One matrix multiply gives het counts per gene per sample:
het_per_gene_sample = gene_variant_map @ het_mask.astype(np.int32)
# Shape: (n_genes, n_samples) — takes seconds

# All compound het candidates at once:
candidate_genes, candidate_samples = np.where(het_per_gene_sample >= 2)
```

This replaces **18,629 gene iterations × 1,452 sample iterations** with a single sparse matrix multiply that completes in seconds.

### Optimization 4: Numba JIT for Pair Enumeration

**Impact: 30-100× speedup for pair detection + parent checking**

After identifying candidate (gene, sample) pairs, the remaining work is enumerating variant pairs and checking parent genotypes. This is inherently iterative but can be compiled with Numba:

```python
from numba import njit, prange

@njit(parallel=True)
def find_compound_het_pairs(gt_matrix, het_mask, candidate_genes, candidate_samples,
                            gene_variant_ranges, parent1_col, parent2_col):
    """
    gt_matrix: (n_variants, n_samples) int8
    het_mask: (n_variants, n_samples) bool
    candidate_genes/samples: 1D arrays from np.where()
    gene_variant_ranges: (n_genes, 2) int — start/end indices
    parent1_col, parent2_col: (n_samples,) int — column indices, -1 if no parent
    """
    n_candidates = len(candidate_genes)
    # Process candidates in parallel across cores (GIL-free)
    for ci in prange(n_candidates):
        gene = candidate_genes[ci]
        sample = candidate_samples[ci]
        start, end = gene_variant_ranges[gene]

        het_vars = np.where(het_mask[start:end, sample])[0] + start
        p1 = parent1_col[sample]
        p2 = parent2_col[sample]

        for i in range(len(het_vars)):
            for j in range(i + 1, len(het_vars)):
                vi, vj = het_vars[i], het_vars[j]
                if p1 >= 0 and p2 >= 0:
                    # Check trans configuration
                    trans = ((gt_matrix[vi, p1] == 1 and gt_matrix[vj, p2] == 1) or
                             (gt_matrix[vi, p2] == 1 and gt_matrix[vj, p1] == 1))
                    # ... store result
```

Key benefits:
- `@njit(parallel=True)` + `prange` bypasses the GIL entirely
- Compiles to LLVM machine code — no Python interpreter overhead
- Automatic thread-level parallelism across candidate (gene, sample) pairs
- int8 array access runs at memory bandwidth speed

### Optimization 5: ProcessPoolExecutor + Shared Memory (Alternative to Numba)

**Impact: True 16× parallelism for GIL-bound code**

If Numba is not desired as a dependency, use `ProcessPoolExecutor` with `multiprocessing.shared_memory` to get real parallelism without pickle overhead:

```python
from multiprocessing.shared_memory import SharedMemory

# Create shared memory genotype matrix (built once, read by all workers)
shm = SharedMemory(create=True, size=gt_matrix.nbytes)
shared_arr = np.ndarray(gt_matrix.shape, dtype=gt_matrix.dtype, buffer=shm.buf)
np.copyto(shared_arr, gt_matrix)

# Workers reconstruct numpy view from shared memory name (zero copy)
def worker(shm_name, shape, dtype, gene_start, gene_end, sample_to_idx, ...):
    existing_shm = SharedMemory(name=shm_name)
    gt = np.ndarray(shape, dtype=dtype, buffer=existing_shm.buf)
    gene_gt = gt[gene_start:gene_end, :]  # zero-copy slice
    # ... analyze compound hets using numpy operations
    existing_shm.close()
    return results  # only return small result dict
```

Each worker reads from the same 409 MB shared memory buffer — no pickling, no copying. Only the small result dicts are serialized back.

### Optimization 6: Vectorize Pass 3

**Impact: ~10× speedup for Pass 3 (from ~800s to ~80s)**

Replace `pd.Series` construction and per-row Python loops:

```python
# Instead of constructing pd.Series per row for create_inheritance_details():
# Pre-extract sample genotypes as numpy columns (already available from gt_matrix)

# Instead of per-row json.dumps():
# Build result columns as numpy string arrays, assign in bulk

# Instead of per-row calculate_segregation_score() with sample loop:
# Vectorize across all variants at once using the genotype matrix
```

At minimum, replace `pd.Series(row_dict)` with a plain dict — `create_inheritance_details()` only does `row[sample_id]` lookups which a dict handles fine.

---

## 6. Projected Performance

| Optimization | Current | Optimized | Speedup |
|-------------|---------|-----------|---------|
| Pass 2 genotype encoding | 27M encode_genotypes() calls | Reuse Pass 1 matrix (0 calls) | **~10×** |
| Het identification | Python loop ×1,452 per gene | Vectorized matrix op | **~50×** |
| Candidate finding | 18,629 gene iterations | One sparse matrix multiply | **~100-500×** |
| Pair enumeration | Python nested loops | Numba @njit(parallel=True) | **~30-100×** |
| Multi-core parallelism | ThreadPoolExecutor (GIL = 1 core) | ProcessPool+SharedMem or Numba prange | **~8-16×** |
| Pass 3 finalization | pd.Series per row, Python loops | Vectorized + dict | **~10×** |

**Conservative overall estimate**: 10-12 hours → **5-30 minutes**

Even just Optimization 1 alone (reuse genotype matrix) would reduce runtime to **~1-2 hours**.

---

## 7. Implementation Priority

| Priority | What | Effort | Impact | Dependencies |
|----------|------|--------|--------|-------------|
| **P0** | Reuse Pass 1 genotype matrix in Pass 2 | Small | 10× | None — add parameter passing |
| **P1** | Vectorize het detection (matrix == 1) | Small | 50× | P0 (needs matrix) |
| **P2** | Sparse matrix candidate finding | Medium | 100×+ | P1 (needs het_mask), scipy dep |
| **P3** | Numba pair enumeration OR ProcessPool+SharedMem | Medium | 30-100× | P0 (needs matrix) |
| **P4** | Vectorize Pass 3 (eliminate pd.Series per row) | Medium | 10× | None |
| **P5** | Pre-compute variant keys once | Small | 2-3× | None |

P0 and P1 together would bring the full AGDE run from 10+ hours to under 1 hour with minimal code changes.

---

## 8. References

- [NumPy boolean array operations](https://jakevdp.github.io/PythonDataScienceHandbook/02.06-boolean-arrays-and-masks.html)
- [scikit-allel GenotypeArray](https://scikit-allel.readthedocs.io/en/stable/model/ndarray.html)
- [Numba parallel prange](https://numba.readthedocs.io/en/stable/user/parallel.html)
- [Python multiprocessing.shared_memory](https://docs.python.org/3/library/multiprocessing.shared_memory.html)
- [Shared memory for NumPy arrays](https://mingze-gao.com/posts/python-shared-memory-in-multiprocessing/)
- [SciPy sparse matrix operations](https://docs.scipy.org/doc/scipy/reference/sparse.html)
- [slivar compound het detection](https://github.com/brentp/slivar) — reference implementation in Nim
- [cyvcf2 fast VCF parsing](https://academic.oup.com/bioinformatics/article/33/12/1867/2971439)
- [Polars vs pandas benchmarks](https://pola.rs/posts/benchmarks/)
- [ThreadPoolExecutor vs GIL](https://superfastpython.com/threadpoolexecutor-vs-gil/)
- [Faster multiprocessing without pickle](https://pythonspeed.com/articles/faster-multiprocessing-pickle/)
- [Free-threaded Python 3.13](https://realpython.com/python313-free-threading-jit/)

---

## 9. Relation to Existing Issues

- **Issue #81**: Categorical dtype on GT columns broke inheritance (fixed v0.14.2) — unmasked this performance issue
- **Issue #82**: ProcessPoolExecutor OOM (fixed v0.14.3 with ThreadPoolExecutor) — solved memory but not performance
- **This analysis**: ThreadPoolExecutor provides zero parallelism; architectural changes needed for scalability

---

*Analysis performed 2026-02-21 on variantcentrifuge v0.14.3 (commit 1d3952a), AGDE exomes full broad cohort.*
