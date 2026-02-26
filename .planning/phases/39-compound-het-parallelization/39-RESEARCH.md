# Phase 39: Compound Het Parallelization - Research

**Researched:** 2026-02-26
**Domain:** Python GIL, ThreadPoolExecutor, NumPy threading, compound het analysis optimization
**Confidence:** HIGH

## Summary

The compound het Pass 2 already uses ThreadPoolExecutor in `parallel_analyzer.py`. The current implementation
passes `gene_df` (a pandas DataFrame slice) to each worker, which triggers GIL-holding operations inside
the workers: `drop_duplicates()`, `duplicated()`, and `iloc[]`. These pandas operations prevent threads
from running truly in parallel despite NumPy being present.

The optimization strategy — locked by CONTEXT.md decisions — keeps ThreadPoolExecutor but eliminates the
GIL-bound Python work from the hot path. Two complementary changes achieve this:

1. **Pre-dispatch transformation**: Move DataFrame operations (`drop_duplicates`, `duplicated`, `iloc`)
   from workers to the main thread's pre-dispatch loop. Workers receive only NumPy arrays.
2. **Pedigree pre-computation**: Replace per-worker pedigree dict lookups (`get_parents`, `is_affected`)
   with three NumPy integer arrays (`father_idx_arr`, `mother_idx_arr`, `affected_arr`) computed once
   before dispatching any work.

After these changes, the worker function's hot path consists almost entirely of NumPy array operations
(GIL-releasing), enabling true thread-level parallelism.

**Primary recommendation:** Implement pre-dispatch deduplication + pedigree arrays, add
`VARIANTCENTRIFUGE_BATCH_SIZE` env var, create a standalone synthetic benchmark, and profile on the GCKD
dataset using py-spy to validate and document the speedup.

---

## Standard Stack

The established tools for this optimization domain:

### Core (already in project)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| numpy | 2.2.6 | GIL-releasing array ops | Already used in gt_matrix; releases GIL for ufuncs, indexing |
| pandas | 2.3.3 | DataFrame pre-processing (main thread only) | Pre-dispatch ops stay in main thread |
| concurrent.futures.ThreadPoolExecutor | stdlib | Thread pool management | Already in place, handles Windows/single-core gracefully |
| psutil | already installed | ResourceManager CPU/memory detection | Already used in ResourceManager |

### Profiling (install before benchmarking)
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| py-spy | 0.4.1 | Sampling profiler for live Python processes | Profile BEFORE optimization to find actual hotspots |

**Installation:**
```bash
pip install py-spy
```

### Supporting (existing project infrastructure)
| Component | Location | Purpose |
|-----------|----------|---------|
| ResourceManager | variantcentrifuge/memory/resource_manager.py | auto_workers(), should_parallelize() |
| synthetic data fixtures | tests/performance/conftest.py | benchmark data generation |
| benchmark_parallelism.py | tests/performance/benchmark_parallelism.py | existing parallel benchmarks |

---

## Architecture Patterns

### Recommended Code Organization

The changes are surgical — no new files needed beyond the standalone benchmark script.

```
variantcentrifuge/inheritance/
├── parallel_analyzer.py         # MODIFY: pre-dispatch, pedigree arrays, env var
├── comp_het_vectorized.py       # MODIFY: new worker signature (numpy-only hot path)
├── analyzer.py                  # NO CHANGE: sequential path untouched
tests/performance/
├── standalone_bench_comp_het_parallel.py   # NEW: standalone synthetic benchmark
```

### Pattern 1: Pre-Dispatch Deduplication

**What:** Compute the dedup mask in the main thread's pre-dispatch loop, pass deduplicated row indices
and variant keys directly to workers. Workers never call `drop_duplicates()` or `duplicated()`.

**When to use:** Any DataFrame operation that is per-gene but identical for all genes can be moved here.

**Current code (parallel_analyzer.py lines 206-211):**
```python
# CURRENT - worker receives gene_df (DataFrame)
for gene, group_df in gene_groups:
    if len(group_df) > 1 and not pd.isna(gene) and gene != "":
        gene_iloc = df.index.get_indexer(group_df.index)
        gene_row_idx = row_positions[gene_iloc]
        gene_vkeys = all_variant_keys[gene_iloc]
        genes_with_multiple_variants.append((gene, group_df, gene_row_idx, gene_vkeys))
```

**Proposed code (pre-dispatch with dedup):**
```python
# PROPOSED - worker receives only numpy arrays
for gene, group_df in gene_groups:
    if len(group_df) <= 1 or pd.isna(gene) or gene == "":
        continue
    gene_iloc = df.index.get_indexer(group_df.index)

    # Dedup in main thread (GIL-holding pandas, but sequential - not in workers)
    dedup_mask = ~group_df.duplicated(subset=["CHROM", "POS", "REF", "ALT"])
    n_unique = int(dedup_mask.sum())
    if n_unique < 2:
        continue

    # Pass only deduplicated numpy arrays to worker (no DataFrame reference)
    dedup_iloc = gene_iloc[dedup_mask.values]
    gene_row_idx_deduped = row_positions[dedup_iloc]
    gene_vkeys_deduped = all_variant_keys[dedup_iloc]
    gene_name_str = str(gene)

    genes_to_dispatch.append(
        (gene_name_str, n_unique, gene_row_idx_deduped, gene_vkeys_deduped)
    )
```

### Pattern 2: Pedigree Array Pre-Computation

**What:** Build three NumPy integer arrays once in the main thread, replacing per-worker dict lookups.
Arrays are indexed by position in `sample_list`.

**When to use:** Any per-worker access to `pedigree_data` dict can be replaced with array indexing.

**Implementation (new function in parallel_analyzer.py):**
```python
def build_pedigree_arrays(
    sample_list: list[str],
    pedigree_data: dict[str, dict[str, Any]],
    sample_to_idx: dict[str, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Pre-compute pedigree relationships as integer arrays for GIL-free worker access.

    Returns arrays indexed by position in sample_list:
      father_idx_arr[i] = index of sample_list[i]'s father in sample_list, or -1
      mother_idx_arr[i] = index of sample_list[i]'s mother in sample_list, or -1
      affected_arr[i]   = 1 if sample_list[i] is affected, else 0
    """
    n = len(sample_list)
    father_idx_arr = np.full(n, -1, dtype=np.int32)
    mother_idx_arr = np.full(n, -1, dtype=np.int32)
    affected_arr = np.zeros(n, dtype=np.int8)

    for i, sample_id in enumerate(sample_list):
        if sample_id not in pedigree_data:
            continue
        info = pedigree_data[sample_id]
        if info.get("affected_status") == "2":
            affected_arr[i] = 1
        father_id = info.get("father_id", "0")
        if father_id and father_id != "0" and father_id in sample_to_idx:
            father_idx_arr[i] = sample_to_idx[father_id]
        mother_id = info.get("mother_id", "0")
        if mother_id and mother_id != "0" and mother_id in sample_to_idx:
            mother_idx_arr[i] = sample_to_idx[mother_id]

    return father_idx_arr, mother_idx_arr, affected_arr
```

**Dtype rationale:**
- `np.int32` for indices: supports cohorts up to 2.1B samples (int16 would overflow at 32767)
- `np.int8` for affected flags: binary value, minimal memory, no overflow risk

### Pattern 3: NumPy-Only Worker Function

**What:** The worker function (`_process_gene_group`) receives no DataFrame. All pedigree lookups
replaced with array indexing. Worker produces the same output dict format.

**Proposed new signature:**
```python
def _process_gene_group_arrays(
    gene: str,
    n_unique: int,
    gene_row_indices: np.ndarray,   # deduplicated row indices into gt_matrix
    variant_keys: np.ndarray,       # deduplicated variant key strings
    gt_matrix: np.ndarray,          # shared read-only (n_total_variants x n_samples)
    sample_to_idx: dict[str, int],  # dict lookup for column position (fast, small)
    father_idx_arr: np.ndarray,     # pre-computed father indices
    mother_idx_arr: np.ndarray,     # pre-computed mother indices
    affected_arr: np.ndarray,       # pre-computed affected flags
    sample_list: list[str],         # needed for result dict keys (sample_id strings)
) -> tuple[str, dict[str, Any]]:
    ...
```

**Worker hot path after optimization:**
1. `gt_matrix[gene_row_indices, :]` — numpy fancy indexing, GIL-releasing
2. `gene_gt_mat == 1` — numpy comparison, GIL-releasing
3. `het_mat.sum(axis=0)` — numpy reduction, GIL-releasing
4. Per-sample: `father_idx_arr[pos]` — numpy scalar indexing, GIL-releasing
5. `find_potential_partners_vectorized()` — already all-numpy, GIL-releasing
6. Result dict building — Python object creation, GIL-holding but O(n_results), not O(n_genes*n_samples)

### Pattern 4: VARIANTCENTRIFUGE_BATCH_SIZE Environment Variable

**What:** Allow HPC users to tune the batch size via environment variable, consistent with the project's
`VARIANTCENTRIFUGE_NO_C_EXT` env var pattern.

**Implementation:**
```python
def _get_batch_size(n_workers: int) -> int:
    """Get batch size from env var or compute default."""
    env_val = os.environ.get("VARIANTCENTRIFUGE_BATCH_SIZE")
    if env_val is not None:
        try:
            size = int(env_val)
            if size >= 1:
                return size
            logger.warning(
                f"VARIANTCENTRIFUGE_BATCH_SIZE={env_val!r} must be >= 1, using auto"
            )
        except ValueError:
            logger.warning(
                f"Invalid VARIANTCENTRIFUGE_BATCH_SIZE={env_val!r}, using auto"
            )
    return min(4 * n_workers, 500)
```

**Validation choice (Claude's discretion):** Warn and fall back on invalid values (don't raise).
Users on HPC may set env vars in cluster config files and a silent fallback is more robust than a crash.

### Pattern 5: Standalone Synthetic Benchmark Script

**What:** A standalone Python script (not a pytest test) that generates synthetic data, runs Pass 2
multiple times with different parameters, and prints timing comparison.

**Location:** `tests/performance/standalone_bench_comp_het_parallel.py`

**Structure:**
```python
#!/usr/bin/env python3
"""
Standalone synthetic benchmark for compound het Pass 2 parallelization.

Run: python tests/performance/standalone_bench_comp_het_parallel.py

Generates 1000 genes / 10000 variants, times Pass 2 with different
min_variants_for_parallel and worker count settings.
"""
import time
import numpy as np
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from variantcentrifuge.inheritance.parallel_analyzer import analyze_inheritance_parallel
from tests.performance.helpers.synthetic_data import (
    generate_synthetic_variants,
    generate_synthetic_pedigree,
)

# ... benchmark loop over [25, 50, 100, 200, 500] thresholds and worker counts
# ... print results table
# No assertions - log timing for human review
```

**Key design decision (Claude's discretion):** Use realistic skewed gene distribution (a few large genes,
many small genes) rather than uniform distribution. Real genomic data has PKD1 with ~4000 transcripts,
BRCA1 with ~200 variants, and many genes with 2-5 variants. Skew reveals load balancing behavior.

### Anti-Patterns to Avoid

- **Passing DataFrames to workers**: `gene_df` in worker = GIL-holding pandas ops in parallel threads
- **Per-worker dict.get() calls**: `pedigree_data[sample_id]` in each worker = O(n_genes * n_samples) Python
- **ProcessPoolExecutor**: Adds pickle overhead for shared `gt_matrix`, which is large. ThreadPoolExecutor
  with numpy ops is both simpler and faster for this use case
- **Dataclass wrapper for worker args**: CONTEXT.md says individual parameters (current style). A dataclass
  would require an extra object allocation and attribute access per worker invocation
- **Removing `gene_df` parameter before profiling confirms benefit**: Profile first, then decide whether
  the sequential path (`analyze_gene_for_compound_het_vectorized` in `analyzer.py`) also needs updating

---

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CPU count detection | Custom platform code | ResourceManager.cpu_cores | Already handles SLURM, PBS, cgroups, psutil |
| Worker count tuning | Custom formula | ResourceManager.auto_workers() | Already accounts for memory + CPU |
| Parallel threshold | Hardcoded constant | ResourceManager.should_parallelize() | Already in use, already parameterized |
| Memory detection | Custom / psutil direct | ResourceManager._detect_memory() | HPC-aware: SLURM, PBS, cgroups, psutil |
| Synthetic test data | New fixture | tests/performance/helpers/synthetic_data.py | generate_synthetic_variants + generate_synthetic_pedigree already exist |
| Profiling infrastructure | Custom timing code | py-spy (external) + existing pass_times dict | pass_times already logs per-pass timing; py-spy for GIL visualization |

**Key insight:** `ResourceManager` already encapsulates all parallelism decisions. The optimization is
surgical — eliminate GIL contention in the already-correct ThreadPoolExecutor architecture.

---

## Common Pitfalls

### Pitfall 1: Breaking the Sequential Path (analyzer.py)

**What goes wrong:** `analyze_gene_for_compound_het_vectorized` is called from BOTH `analyzer.py`
(sequential path) and `parallel_analyzer.py` (parallel path). Changing the signature to remove
`gene_df` would break `analyzer.py`.

**Why it happens:** The two paths share the vectorized worker function.

**How to avoid:** Create a NEW worker function `_process_gene_group_arrays()` for the parallel path.
Keep `analyze_gene_for_compound_het_vectorized()` unchanged for the sequential path. Shared logic
(e.g., `find_potential_partners_vectorized`) remains in `comp_het_vectorized.py`.

**Warning signs:** If `tests/test_inheritance/` fails after the change, the sequential path was modified.

### Pitfall 2: Forgetting the Dedup Mask Alignment

**What goes wrong:** When computing `dedup_mask = ~group_df.duplicated(...)` in the main thread,
the resulting indices must be aligned with `gene_row_idx` (which is already based on `gene_iloc`).
Using `dedup_mask.values` on a non-contiguous index produces correct results only if the index
positions align with `gene_iloc`.

**Why it happens:** `group_df` has the original DataFrame index, not positional. `gene_iloc` is
already positional. Double-indexing `gene_iloc[dedup_mask.values]` is correct because both are
positional within the gene group.

**How to avoid:** Always slice `gene_row_idx` and `gene_vkeys` using `dedup_mask.values` (numpy bool
array), not the pandas index.

**Warning signs:** Test `test_comp_het_split_lines.py` catches false positives from duplicated
transcript rows — run this test specifically after the change.

### Pitfall 3: int16 Overflow for Large Cohorts

**What goes wrong:** `np.int16` max is 32767. A cohort with >32767 samples would silently overflow
father/mother indices, producing wrong pedigree lookups.

**Why it happens:** int16 seems natural for indices but the limit is too low for biobank-scale data.

**How to avoid:** Use `np.int32` for father/mother index arrays. (See Pattern 2 above.)

**Warning signs:** Not caught by existing tests (which use 3-sample trios). Would only manifest with
extremely large cohorts.

### Pitfall 4: Profiling Without a Baseline

**What goes wrong:** Implementing the optimization without capturing a baseline py-spy profile means
there's no evidence of what changed or whether it helped.

**Why it happens:** Temptation to go straight to implementation.

**How to avoid:** The CONTEXT.md explicitly says "Profile first, then optimize." Capture py-spy
flamegraph BEFORE any code changes. Commands:
```bash
# Profile BEFORE optimization
py-spy record -o before_opt.svg -- python standalone_bench_comp_het_parallel.py

# Profile AFTER optimization
py-spy record -o after_opt.svg -- python standalone_bench_comp_het_parallel.py
```

**Warning signs:** If the plan doesn't have a "capture baseline" task before any implementation task.

### Pitfall 5: min_variants_for_parallel Threshold Selection

**What goes wrong:** Choosing 100 (current default) without evidence may mean parallelism activates
too early (ThreadPoolExecutor overhead dominates for small datasets) or too late (loses speedup for
medium datasets).

**Why it happens:** Current default is arbitrary.

**How to avoid:** Benchmark at 25, 50, 100, 200, 500 on GCKD dataset. Find empirical crossover where
parallel is faster than sequential. Set default to that value.

**Warning signs:** If benchmark shows sequential is faster at the default threshold.

### Pitfall 6: Worker Still Holds gene_df Reference

**What goes wrong:** After optimization, if `gene_df` is still captured in a closure or passed as an
argument, threads still compete for the GIL when any pandas operation accesses it.

**Why it happens:** Easy to miss one access path.

**How to avoid:** After removing `gene_df` from the worker signature, verify with a grep:
```bash
grep "gene_df" variantcentrifuge/inheritance/parallel_analyzer.py
```
The pre-dispatch loop may reference `group_df` for dedup, but this happens in the main thread.

---

## Code Examples

Verified patterns from codebase inspection:

### Existing: Pre-built Matrix Slicing (already GIL-releasing)
```python
# Source: variantcentrifuge/inheritance/comp_het_vectorized.py lines 145-152
# This path is already in use when has_prebuilt=True
if has_prebuilt:
    kept_indices = gene_row_indices[dedup_mask.values]  # numpy fancy index
    for sample_id in sample_list:
        if sample_id in sample_to_idx:
            genotype_matrix[sample_id] = gt_matrix[kept_indices, sample_to_idx[sample_id]]
        else:
            genotype_matrix[sample_id] = np.full(len(gene_df_unique), -1, dtype=np.int8)
```

### Existing: Candidate Detection (already GIL-releasing)
```python
# Source: variantcentrifuge/inheritance/comp_het_vectorized.py lines 162-170
# np.column_stack, == operator, .sum(axis=0) — all GIL-releasing
if sample_ids_present:
    gene_gt_mat = np.column_stack([genotype_matrix[s] for s in sample_ids_present])
    het_mat = gene_gt_mat == 1          # bool mask, GIL-releasing
    het_counts = het_mat.sum(axis=0)    # reduction, GIL-releasing
    candidate_mask = het_counts >= 2
    candidate_samples = [sample_ids_present[i] for i in np.where(candidate_mask)[0]]
```

### Existing: VARIANTCENTRIFUGE_NO_C_EXT Pattern (precedent for env vars)
```python
# Source: variantcentrifuge/association/backends/davies.py lines 73-75
if os.environ.get("VARIANTCENTRIFUGE_NO_C_EXT"):
    logger.warning(
        "VARIANTCENTRIFUGE_NO_C_EXT set: Davies C extension disabled. "
    )
```

### Existing: Deduplication in Sequential Path (move this to pre-dispatch in parallel path)
```python
# Source: variantcentrifuge/inheritance/comp_het_vectorized.py lines 123-128
# CURRENTLY in worker - needs to move to parallel_analyzer pre-dispatch loop
gene_df_unique = gene_df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
if len(gene_df_unique) < 2:
    return comp_het_results
gene_name = gene_df_unique.iloc[0].get("GENE", "Unknown")
dedup_mask = ~gene_df.duplicated(subset=["CHROM", "POS", "REF", "ALT"])
```

### Existing: ResourceManager.auto_workers() (already correct, use as-is)
```python
# Source: variantcentrifuge/memory/resource_manager.py lines 281-326
# Already handles: memory constraints, cpu constraints, task count constraints
n_workers = rm.auto_workers(
    task_count=num_genes, memory_per_task_gb=memory_per_gene_gb
)
```

### Existing: Pass Timing Infrastructure (already present)
```python
# Source: parallel_analyzer.py — pass_times dict already tracks per-pass timing
pass_times = {}
pass_times["pattern_deduction"] = time.time() - pass1_start
pass_times["compound_het_analysis"] = time.time() - pass2_start
# Add more granular sub-timing inside Pass 2 if needed for profiling evidence
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| df.apply(per-row Python) | vectorized_deduce_patterns (Pass 1) | Phase 09 | 4-8x Pass 1 speedup |
| Sequential gene loop | ThreadPoolExecutor (Pass 2) | Phase 12 | Parallelism infrastructure in place |
| Per-worker dict lookups | (target) numpy pedigree arrays | Phase 39 | Eliminate GIL contention |
| gene_df passed to workers | (target) pre-dispatch dedup + numpy only | Phase 39 | Workers hold GIL less |

**Deprecated/outdated:**
- `use_vectorized_comp_het` parameter: Deprecated in both `analyzer.py` and `parallel_analyzer.py`,
  kept for backward compatibility. Do not use it for new code.
- `comp_het.py` (original non-vectorized): Superseded by `comp_het_vectorized.py`. Not in scope for
  this phase.

---

## Open Questions

1. **Actual GIL contention fraction in current code**
   - What we know: `drop_duplicates`, `duplicated`, `iloc` hold GIL; numpy ops release it
   - What's unclear: What percentage of worker wall-time is GIL-holding in practice on real data?
   - Recommendation: py-spy flamegraph BEFORE optimization will answer this definitively

2. **Optimal min_variants_for_parallel threshold**
   - What we know: Current default is 100; ResourceManager.min_items_for_parallel=100
   - What's unclear: The actual crossover point on GCKD data (100? 50? 200?)
   - Recommendation: Sweep 25, 50, 100, 200, 500 in the benchmark; set default to empirical crossover

3. **Whether the apply phase benefits from optimization**
   - What we know: CONTEXT.md explicitly defers apply phase optimization to profiling evidence
   - What's unclear: How much time the apply phase takes relative to Pass 2 computation
   - Recommendation: Profile it; if it's <5% of Pass 2 time, leave it alone

4. **py-spy on WSL2 environment**
   - What we know: py-spy 0.4.1 is available for Linux x86_64 (confirmed via pip dry-run)
   - What's unclear: Whether sudo is required for py-spy on WSL2 to access process memory
   - Recommendation: Try `py-spy record --subprocesses` first; use `sudo` if permission denied

---

## Sources

### Primary (HIGH confidence)
- Codebase: `variantcentrifuge/inheritance/parallel_analyzer.py` — full inspection of current Pass 2 implementation
- Codebase: `variantcentrifuge/inheritance/comp_het_vectorized.py` — full inspection of worker function
- Codebase: `variantcentrifuge/memory/resource_manager.py` — ResourceManager API
- Codebase: `variantcentrifuge/ped_reader.py` — `get_parents()`, `is_affected()` implementation
- Codebase: `.planning/phases/39-compound-het-parallelization/39-CONTEXT.md` — locked decisions
- Runtime: `numpy.__version__ == 2.2.6`, `pandas.__version__ == 2.3.3` — confirmed
- Runtime: `pip show py-spy` — version 0.4.1 available for install

### Secondary (MEDIUM confidence)
- Live micro-test: DataFrame ops (100x sequential) = 0.134s; NumPy equivalent = 0.001s
  — confirms orders-of-magnitude difference confirming GIL-holding vs GIL-releasing
- Live thread test: NumPy parallel (4 threads) achieves ~2.2x vs sequential; Pandas ~1.5x
  — confirms NumPy has better GIL behavior in threads (though test scale was small)
- Codebase: `tests/performance/benchmark_parallelism.py` — existing benchmark patterns and fixtures
- Codebase: `.planning/REQUIREMENTS.md` — PERF-01 wording updated by CONTEXT.md decisions

### Tertiary (LOW confidence)
- CONTEXT.md statement: "research confirmed: pandas DataFrame operations hold the GIL even for
  simple indexing; NumPy array ops release it" — stated as confirmed but original source not cited

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — verified by inspection of actual installed versions and existing code
- Architecture patterns: HIGH — directly derived from reading current implementation
- Don't hand-roll: HIGH — ResourceManager, synthetic data fixtures confirmed to exist and work
- Pitfalls: HIGH (code correctness), MEDIUM (performance estimates) — correctness pitfalls from
  code inspection; performance estimates from micro-benchmarks
- Profiling approach: HIGH — py-spy confirmed available; commands from official py-spy docs pattern

**Research date:** 2026-02-26
**Valid until:** 2026-09-26 (stable Python stdlib + NumPy, unlikely to change)
