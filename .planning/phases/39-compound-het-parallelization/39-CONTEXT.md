# Phase 39: Compound Het Parallelization - Context

**Gathered:** 2026-02-26
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace the GIL-bound compound het Pass 2 loop with optimized parallelism for measurable speedup on multi-core systems. The existing ThreadPoolExecutor architecture is correct — optimization focuses on eliminating GIL contention in worker code, not switching executor type. All existing compound het tests must pass unchanged.

</domain>

<decisions>
## Implementation Decisions

### Parallelism strategy
- Keep ThreadPoolExecutor (NOT ProcessPoolExecutor or Numba)
- ThreadPoolExecutor is correct because: gt_matrix is shared read-only NumPy (GIL-releasing), ProcessPoolExecutor adds serialization overhead for small per-gene DataFrames, Numba can't handle pedigree dicts/pandas and has silent-failure risk
- Focus optimization on eliminating Python-object GIL contention inside workers
- Pre-compute pedigree lookups as integer arrays (father_indices, mother_indices, affected_flags as int8/int16 NumPy arrays) — replaces O(n_genes x n_samples) dict accesses in workers
- Move DataFrame operations (drop_duplicates, iloc, duplicated) from workers to main thread pre-dispatch
- Update roadmap success criteria to reflect ThreadPoolExecutor optimization (not replacement with ProcessPoolExecutor)

### Data transfer to workers
- Workers receive only NumPy arrays — no DataFrame references in the hot path
- Individual array parameters (current function signature style), not a dataclass wrapper
- Pre-dispatch loop in main thread handles: deduplication, variant key slicing, row index computation, pedigree array pre-computation
- Result format stays as dict (gene_name -> {variant_key: info_dict}) — apply phase is sequential, not a bottleneck
- Apply phase optimization deferred to profiling evidence

### Fallback thresholds
- No platform-specific fallback needed — ThreadPoolExecutor handles single-core (1 worker) and Windows gracefully; ResourceManager already adapts
- min_variants_for_parallel threshold: benchmark at multiple values (25, 50, 100, 200, 500) on GCKD dataset to find empirical crossover point
- ResourceManager auto_workers(): tune only if profiling shows suboptimal worker count
- Batch size (currently min(4*n_workers, 500)): benchmark at different values + make configurable via VARIANTCENTRIFUGE_BATCH_SIZE environment variable for HPC users

### Benchmark approach
- Profile with py-spy before AND after optimization on GCKD real dataset (testing/gckd_all.GRCh37.annotated.vcf.gz with testing/morbidgenes_500plusPKD1u2.txt gene list)
- Synthetic benchmark: standalone script in tests/performance/ (NOT pytest test), 1000 genes / 10000 variants
- Real benchmark: GCKD dataset wall-time comparison of Pass 2 before/after
- No assertion thresholds — log timing comparison for human review, avoid flaky test issues from hardware variance

### Claude's Discretion
- Exact pedigree array layout and dtype choices
- Whether to split pre-computation into a separate function or inline in analyze_inheritance_parallel
- py-spy vs scalene for profiling tool choice
- Synthetic benchmark gene/variant distribution (uniform vs realistic)
- Logging granularity for timing breakdown
- Whether batch_size env var needs validation/bounds

</decisions>

<specifics>
## Specific Ideas

- "Profile first, then optimize" — user wants empirical evidence driving each optimization decision
- GCKD test dataset available at: testing/gckd_all.GRCh37.annotated.vcf.gz with gene list testing/morbidgenes_500plusPKD1u2.txt, samples in testing/PKD_samples.txt and testing/nonPKD_samples.txt
- Environment variable pattern for HPC configurability (like scikit-learn's LOKY_MAX_CPU_COUNT, NumPy's OMP_NUM_THREADS)
- Research confirmed: pandas DataFrame operations hold the GIL even for simple indexing; NumPy array ops release it

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 39-compound-het-parallelization*
*Context gathered: 2026-02-26*
