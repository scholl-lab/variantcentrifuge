# Phase 12: Parallelization & Chunking - Research

**Researched:** 2026-02-16
**Domain:** Python multiprocessing, pandas DataFrame chunking, memory management
**Confidence:** HIGH

## Summary

Phase 12 requires careful investigation before implementation. Phases 7-11 transformed the pipeline dramatically:

- **Phase 7:** 48-98% gene burden speedup, removed dead GT parsing loop
- **Phase 8:** 82-84% memory reduction via dtype optimization
- **Phase 11:** Eliminated 7-hour genotype replacement stage and 2.7-hour SnpSift (replaced with bcftools 19.4x faster)

**Current state analysis:**
- `inheritance_memory_manager.py` (385 lines) provides comprehensive HPC environment detection (SLURM, PBS, cgroups), chunk size calculation, and parallelism tuning
- `parallel_analyzer.py` (333 lines) implements gene-level parallel processing for compound het detection with ProcessPoolExecutor
- `--chunks` CLI flag exists but is NOT used by inheritance analysis (only for legacy file processing)
- `--threads` CLI flag controls ProcessPoolExecutor max_workers (default: 1)
- No memory pools, no async I/O, no work stealing currently implemented

**Primary recommendation:** Investigate real bottlenecks FIRST via profiling, then implement only what's proven necessary. Don't optimize based on assumptions from pre-Phase-11 pipeline.

## Standard Stack

### Core Libraries (Already in Use)

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `multiprocessing` | stdlib | Process-based parallelism | Python standard for CPU-bound parallel tasks, used by ProcessPoolExecutor |
| `concurrent.futures` | stdlib | High-level threading/process interface | Industry-standard abstraction over multiprocessing, thread-safe |
| `psutil` | 5.9+ | System/process resource monitoring | De facto standard for cross-platform memory/CPU detection |
| `pandas` | 2.0+ | DataFrame processing | Core data structure, chunksize support built-in |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `numpy` | 1.24+ | Numerical operations | Already used for vectorized genotype operations |
| `tracemalloc` | stdlib | Python memory profiling | Used in benchmark_memory_budgets.py for peak memory tracking |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| ProcessPoolExecutor | multiprocessing.Pool | Pool requires manual chunk management, ProcessPoolExecutor higher-level |
| psutil | /proc parsing | psutil is cross-platform (Linux/Mac/Windows), /proc is Linux-only |
| Manual chunking | Dask/Polars | Adds heavy dependencies for marginal benefit on already-optimized pipeline |

**Installation:** No new dependencies required (all stdlib or already installed).

## Architecture Patterns

### Current Parallelism Implementation

**Gene-level parallel processing (parallel_analyzer.py):**
```
analyze_inheritance_parallel()
  ├── Pass 1: vectorized_deduce_patterns() [sequential]
  ├── Pass 2: Compound Het [PARALLEL]
  │   ├── Group by GENE
  │   ├── Submit to ProcessPoolExecutor
  │   └── Collect results via as_completed()
  └── Pass 3: Prioritization [sequential]
```

**Parallelism control:**
- `n_workers` parameter (defaults to None → uses CPU count)
- `min_variants_for_parallel` threshold (default: 100)
- Auto-disables for small datasets

**Current worker detection:**
```python
# InheritanceMemoryManager.calculate_optimal_parallelism()
cpu_cores = psutil.cpu_count(logical=False) or 4
max_workers = min(max_parallel_chunks, cpu_cores, len(chunk_sizes))
```

### Current Memory Management Implementation

**Memory detection hierarchy (inheritance_memory_manager.py):**
```
1. CLI parameter (--max-memory-gb)
2. SLURM_MEM_PER_NODE environment variable
3. PBS_RESC_MEM environment variable
4. cgroup limits (/sys/fs/cgroup/memory/memory.limit_in_bytes)
5. psutil.virtual_memory().available (fallback)
```

**Chunk size calculation:**
```python
# Factors:
bytes_per_sample_column = 8  # float64
pandas_overhead_factor = 2.5
inheritance_analysis_factor = 3.0

# Formula:
memory_per_variant = (num_samples * 8 * 2.5 * 3.0) / (1024**3)
max_variants = int(target_memory_gb / memory_per_variant)
max_variants = max(100, min(max_variants, 1_000_000))  # Bounds
```

**Memory strategy selection:**
```python
if full_dataset_fits_in_memory:
    return {"approach": "full_dataset", "chunks": 1}
else:
    return {"approach": "chunked", "chunks": N, "max_workers": W}
```

### Pattern 1: ProcessPoolExecutor with Dynamic Worker Count

**What:** Calculate optimal worker count based on CPU cores and memory constraints
**When to use:** CPU-bound tasks with memory considerations (inheritance analysis)
**Example:**
```python
# From parallel_analyzer.py lines 163-174
with ProcessPoolExecutor(max_workers=n_workers) as executor:
    future_to_gene = {
        executor.submit(
            _process_gene_group,
            gene,
            gene_df,
            pedigree_data,
            sample_list,
        ): gene
        for gene, gene_df in genes_with_multiple_variants
    }

    for future in as_completed(future_to_gene):
        gene_name, comp_het_results = future.result()
```

**Source:** Standard concurrent.futures pattern, see [Python ProcessPoolExecutor docs](https://docs.python.org/3/library/concurrent.futures.html)

### Pattern 2: Environment-Aware Memory Detection

**What:** Detect memory limits from HPC schedulers before falling back to system total
**When to use:** Code running in SLURM/PBS/containerized environments
**Example:**
```python
# From inheritance_memory_manager.py lines 78-87
slurm_mem = os.getenv("SLURM_MEM_PER_NODE")
if slurm_mem:
    try:
        slurm_gb = float(slurm_mem) / 1024  # SLURM uses MB
        logger.info(f"Using SLURM allocated memory: {slurm_gb:.1f}GB")
        return slurm_gb
    except (ValueError, TypeError):
        logger.warning(f"Invalid SLURM_MEM_PER_NODE value: {slurm_mem}")
```

**Source:** HPC best practices for memory-aware applications

### Pattern 3: Pandas Chunked Processing

**What:** Process large DataFrames in chunks to stay within memory limits
**When to use:** When full dataset exceeds available memory
**Example:**
```python
# Pattern from pandas official docs (not currently used for inheritance)
for chunk in pd.read_csv(large_file, chunksize=10000):
    process_chunk(chunk)
    results.append(chunk_result)
combined = pd.concat(results)
```

**Source:** [Pandas scaling docs](https://pandas.pydata.org/docs/user_guide/scale.html)

### Anti-Patterns to Avoid

- **Premature parallelization:** Adding multiprocessing to tasks that aren't CPU-bound (I/O waits dominate)
- **Too-small chunks:** Overhead of inter-process communication dominates actual work
- **Too-large chunks:** Memory pressure causes swapping, negating performance gains
- **Ignoring GIL:** Threading won't help CPU-bound Python code, ProcessPoolExecutor required

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Worker pool management | Custom process spawning | ProcessPoolExecutor | Handles worker lifecycle, exception propagation, resource cleanup |
| Memory detection | /proc file parsing | psutil.virtual_memory() | Cross-platform, handles edge cases (containers, cgroups) |
| Chunk size tuning | Fixed hardcoded values | Dynamic calculation based on memory/data | Data size varies 1000x between small and large cohorts |
| Process communication | Pipes/queues | submit() + as_completed() | Automatic result collection, progress tracking |

**Key insight:** Python stdlib concurrent.futures provides robust abstractions that handle edge cases (process crashes, pickling errors, resource leaks) that naive implementations miss.

## Common Pitfalls

### Pitfall 1: Default Chunk Size Assumes Equal Work Distribution

**What goes wrong:** ProcessPoolExecutor's default chunksize calculation (`len(iterable) / (4 * num_workers)`) assumes uniform task duration. For inheritance analysis, genes have vastly different variant counts (1 variant vs 5000 variants), causing load imbalance.

**Why it happens:** Default formula optimized for homogeneous tasks. Gene processing time is O(n²) for compound het detection, so one large gene blocks all workers.

**How to avoid:** Pre-sort genes by variant count descending, submit largest first. Alternatively, split large genes into sub-chunks.

**Warning signs:**
- `htop` shows most workers idle while 1-2 busy
- Total time dominated by longest-running gene
- Log messages show uneven "completed N/M genes" progress

**Code location:** `parallel_analyzer.py:148-174` groups genes but doesn't sort by size.

### Pitfall 2: Memory Detection Fails in Containers Without cgroup v2

**What goes wrong:** `psutil.virtual_memory()` returns host system memory (256 GB) when running in a Docker container with 8 GB limit. Memory manager calculates chunk sizes assuming 256 GB, causing OOM kills.

**Why it happens:** cgroup v1 limits not exposed via psutil on all platforms. Code checks `/sys/fs/cgroup/memory/memory.limit_in_bytes` but not all container runtimes expose this.

**How to avoid:**
1. Check BOTH cgroup v1 and v2 paths
2. Validate detected memory against `ulimit -m` if available
3. Add `--max-memory-gb` CLI override documentation for container users

**Warning signs:**
- Logs show allocated memory > container limit
- OOM killer messages in `dmesg`
- Process RSS exceeds expected values

**Code location:** `inheritance_memory_manager.py:124-142` cgroup detection needs v2 path validation.

### Pitfall 3: Pickling Overhead Dominates Small Tasks

**What goes wrong:** For small datasets (<100 variants), the time to pickle/unpickle DataFrames and send to worker processes exceeds the actual analysis time. Parallel version slower than sequential.

**Why it happens:** ProcessPoolExecutor serializes all arguments via pickle. Large DataFrames (with nested dicts in columns like pedigree_data) are expensive to serialize.

**How to avoid:**
1. Enforce minimum dataset size threshold (already implemented: `min_variants_for_parallel=100`)
2. For very large pedigrees, consider shared memory (but adds complexity)
3. Benchmark sequential vs parallel and choose dynamically

**Warning signs:**
- Parallel version slower than sequential for small datasets
- High system CPU time (vs user CPU time) in profiling
- `strace` shows lots of pipe/socket I/O

**Code location:** `parallel_analyzer.py:140-145` checks threshold but value may need tuning.

### Pitfall 4: Phase 11 Eliminated the Chunking Need

**What goes wrong:** Implementing complex chunking infrastructure for genotype processing when Phase 11 already eliminated the 7-hour bottleneck.

**Why it happens:** Original Phase 12 requirements written before Phase 11's bcftools transformation. Old assumptions about bottlenecks no longer valid.

**How to avoid:** Profile current pipeline FIRST. If inheritance analysis is now <5% of total time, chunking won't help.

**Warning signs:**
- Profiling shows inheritance <10% of wall time
- I/O (bcftools) dominates execution time
- Chunking adds complexity for <5% speedup

**Code location:** Phase 12 should START with profiling before implementing anything.

## Code Examples

### Auto-Detecting Optimal Worker Count

```python
# Source: inheritance_memory_manager.py:221-264
def calculate_optimal_parallelism(
    self, chunk_sizes: list, num_samples: int
) -> tuple[int, dict[str, Any]]:
    """Calculate optimal number of parallel workers based on memory constraints."""
    safe_memory_gb = self._allocated_memory_gb * self.memory_safety_factor
    inheritance_memory_gb = safe_memory_gb * self.inheritance_memory_fraction

    # Memory constraint: how many chunks fit in parallel?
    max_chunk_size = max(chunk_sizes) if chunk_sizes else 1000
    memory_per_chunk = self.estimate_chunk_memory_requirement(max_chunk_size, num_samples)
    max_parallel_chunks = max(1, int(inheritance_memory_gb / memory_per_chunk))

    # CPU constraint: physical cores
    cpu_cores = psutil.cpu_count(logical=False) or 4

    # Final decision: minimum of all constraints
    max_workers = min(max_parallel_chunks, cpu_cores, len(chunk_sizes))

    logger.info(
        f"Optimal parallelism: {max_workers} workers "
        f"(memory limit: {max_parallel_chunks}, CPU cores: {cpu_cores})"
    )

    return max_workers, memory_info
```

**Key principles:**
1. Memory constraint calculated from chunk size × parallelism
2. CPU constraint from physical cores (not logical — avoids hyperthreading)
3. Task constraint from number of chunks (no point in 8 workers for 3 chunks)
4. Take minimum of all constraints

### Adaptive Small-Dataset Detection

```python
# Source: parallel_analyzer.py:140-145
use_parallel = (
    (n_workers is None or n_workers > 1)
    and len(df) >= min_variants_for_parallel
    and "GENE" in df.columns
)

if use_parallel and "GENE" in df.columns:
    # Use ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        ...
else:
    # Fall back to sequential
    for gene, gene_df in df.groupby("GENE", observed=True):
        comp_het_results = analyze_gene_for_compound_het_vectorized(...)
```

**Why this matters:** Avoids pickling overhead for datasets where parallel isn't beneficial.

### HPC Environment Detection

```python
# Source: inheritance_memory_manager.py:58-123
def _get_allocated_memory(self) -> float:
    """Get allocated memory limit in GB.

    Priority:
    1. CLI parameter (--max-memory-gb)
    2. SLURM_MEM_PER_NODE
    3. PBS_RESC_MEM
    4. cgroup limits
    5. Available system memory
    """
    # 1. CLI override
    cli_memory = self.config.get("max_memory_gb")
    if cli_memory:
        logger.info(f"Using CLI-specified memory limit: {cli_memory:.1f}GB")
        return float(cli_memory)

    # 2. SLURM
    slurm_mem = os.getenv("SLURM_MEM_PER_NODE")
    if slurm_mem:
        slurm_gb = float(slurm_mem) / 1024  # SLURM uses MB
        logger.info(f"Using SLURM allocated memory: {slurm_gb:.1f}GB")
        return slurm_gb

    # 3. PBS (handles "16gb", "16000mb" formats)
    pbs_mem = os.getenv("PBS_RESC_MEM")
    if pbs_mem:
        pbs_mem_lower = pbs_mem.lower()
        if pbs_mem_lower.endswith("gb"):
            pbs_gb = float(pbs_mem_lower.replace("gb", ""))
        elif pbs_mem_lower.endswith("mb"):
            pbs_gb = float(pbs_mem_lower.replace("mb", "")) / 1024
        logger.info(f"Using PBS allocated memory: {pbs_gb:.1f}GB")
        return pbs_gb

    # 4. cgroup (containers)
    cgroup_limit = self._get_cgroup_memory_limit()
    if cgroup_limit:
        return cgroup_limit

    # 5. Fallback
    available_gb = psutil.virtual_memory().available / (1024**3)
    return float(available_gb)
```

**Robustness:** Falls back gracefully through 5 detection methods.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Fixed chunk sizes | Dynamic calculation based on memory + data size | v0.12.0 | Adaptive to environment (HPC vs laptop) |
| multiprocessing.Pool | concurrent.futures.ProcessPoolExecutor | v0.11.0+ | Better exception handling, cleaner API |
| Manual worker management | Auto-detect CPU cores + memory limit | v0.12.0 | No tuning required for different systems |
| Hardcoded parallelism | Adaptive thresholds (min_variants_for_parallel) | v0.12.0 | Small datasets use sequential automatically |

**Deprecated/outdated:**
- **--chunks CLI flag:** Only used for legacy file processing, NOT used by inheritance analysis (lines: cli.py:753, 1303)
- **Fixed SLURM detection:** Old versions only checked SLURM, new code checks PBS/cgroups too
- **Threading for analysis:** GIL makes threading useless for CPU-bound Python; ProcessPoolExecutor standard now

## Open Questions

### 1. What are the ACTUAL bottlenecks post-Phase-11?

**What we know:**
- Phase 11 eliminated 7-hour genotype replacement
- Phase 11 eliminated 2.7-hour SnpSift (replaced with 8-min bcftools)
- Phase 7-9 optimized inheritance analysis significantly
- Old profiling data shows inheritance was 5-30 min of 10-hour pipeline

**What's unclear:**
- Where does time go NOW in the transformed pipeline?
- Is inheritance analysis still a bottleneck or is it bcftools I/O?
- Do large cohorts (5000+ samples) hit memory limits with current code?

**Recommendation:** START Phase 12 with profiling runs on realistic data (GCKD 5125-sample cohort) to identify what actually needs optimization.

### 2. Is memory pooling actually needed given Phase 8 gains?

**What we know:**
- Phase 8 achieved 82-84% memory reduction via dtype optimization
- Current memory manager uses conservative safety factors (92% of allocated, 85% for inheritance)
- No memory pools currently implemented

**What's unclear:**
- Does Phase 8 dtype optimization eliminate need for reusable buffers?
- Are allocation/deallocation overhead significant in profiling?
- Would memory pools add complexity for marginal benefit?

**Recommendation:** Use `tracemalloc` to measure allocation patterns. If <5% of time is spent in malloc/free, skip memory pools.

### 3. Should --chunks CLI flag be removed or repurposed?

**What we know:**
- `--chunks` exists (cli.py:753) but is NOT used by InheritanceMemoryManager
- Chunk size calculated dynamically based on memory/data size
- CONTEXT.md says "remove --chunks CLI flag"

**What's unclear:**
- Does legacy code still use `config["chunks"]`?
- Safe to remove or does it break backward compatibility?

**Recommendation:** Grep for all `config["chunks"]` usage, verify none are critical, then remove CLI flag and document auto-chunking in help text.

### 4. Is async I/O beneficial when bcftools is the I/O bottleneck?

**What we know:**
- Phase 11 replaced SnpSift with bcftools query (19.4x faster)
- bcftools is a subprocess call, not Python I/O
- Async I/O helps when Python reads/writes are interleaved with processing

**What's unclear:**
- Can we overlap bcftools subprocess with downstream processing?
- Would async file reading help for TSV/Excel output generation?
- Is I/O bound or CPU bound in current pipeline?

**Recommendation:** Profile I/O wait time. If processes spend >30% in iowait, investigate async. Otherwise, skip.

## Sources

### Primary (HIGH confidence)

- **variantcentrifuge codebase analysis:**
  - `variantcentrifuge/memory/inheritance_memory_manager.py` (385 lines) — comprehensive memory detection + chunking
  - `variantcentrifuge/inheritance/parallel_analyzer.py` (333 lines) — ProcessPoolExecutor implementation
  - `variantcentrifuge/stages/analysis_stages.py:840-1090` — InheritanceAnalysisStage integration
  - `.planning/phases/07-quick-wins-tier-1/07-VERIFICATION.md` — Phase 7 performance gains (48-98% speedup)
  - `.planning/phases/11-pipeline-io-elimination/11-03-SUMMARY.md` — Phase 11 transformation details
  - `.planning/performance-analysis-report.md` — Pre-Phase-7 baseline (GCKD cohort analysis)

- **Python official documentation:**
  - [concurrent.futures — Launching parallel tasks](https://docs.python.org/3/library/concurrent.futures.html) — ProcessPoolExecutor API
  - [multiprocessing — Process-based parallelism](https://docs.python.org/3/library/multiprocessing.html) — Process pool fundamentals

- **Pandas official documentation:**
  - [Scaling to large datasets](https://pandas.pydata.org/docs/user_guide/scale.html) — Chunked processing patterns

### Secondary (MEDIUM confidence)

- [How to Configure Multiprocessing Pool.map() Chunksize - Super Fast Python](https://superfastpython.com/multiprocessing-pool-map-chunksize/) — Chunksize calculation best practices
- [Configure Max Workers For The ProcessPoolExecutor - Super Fast Python](https://superfastpython.com/processpoolexecutor-number-of-workers/) — Worker count tuning guidance
- [Data and chunk sizes matter when using multiprocessing.Pool.map() in Python - Medium](https://rvprasad.medium.com/data-and-chunk-sizes-matter-when-using-multiprocessing-pool-map-in-python-5023c96875ef) — Empirical chunking analysis
- [Scalable Data Processing with Pandas: Handling Large CSV Files in Chunks - Medium](https://medium.com/itversity/scalable-data-processing-with-pandas-handling-large-csv-files-in-chunks-b15c5a79a3e3) — Pandas chunking patterns
- [Efficient Pandas: Using Chunksize for Large Datasets - Towards AI](https://towardsai.net/p/data-science/efficient-pandas-using-chunksize-for-large-data-sets-c66bf3037f93) — Chunksize parameter usage

### Tertiary (LOW confidence)

- WebSearch results from February 2026 queries about multiprocessing best practices — general guidance, not project-specific

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — All libraries already in use, well-documented
- Architecture: HIGH — Codebase analysis complete, patterns verified in production code
- Pitfalls: HIGH — Based on actual code inspection and Phase 11 context
- Memory pools/async I/O necessity: LOW — Requires profiling to determine

**Research date:** 2026-02-16
**Valid until:** 90 days (stable Python stdlib APIs, pandas patterns)

**Critical insight:** Phase 12 scope MUST be driven by profiling, not assumptions. Original requirements (PARLZ-02, PARLZ-03, PARLZ-04) were written before Phase 11's 19.4x bcftools speedup and genotype replacement elimination. **Investigate first, implement only what profiling proves necessary.**

**Test coverage findings:**
- InheritanceMemoryManager: NO dedicated unit tests found
- parallel_analyzer.py: 1 test file (`tests/unit/inheritance/test_parallel_analyzer.py`, 11,619 bytes)
- Memory budgets: `tests/performance/benchmark_memory_budgets.py` tracks peak memory with tracemalloc
- Integration tests: Multiple files test parallel processing end-to-end

**Gap:** InheritanceMemoryManager lacks unit tests for environment detection (SLURM, PBS, cgroups). Should add tests before refactoring to pipeline-wide module.
