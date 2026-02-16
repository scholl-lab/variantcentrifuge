# Requirements: VariantCentrifuge v0.13.0 Performance Optimization

**Defined:** 2026-02-14
**Core Value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs — now 3-4x faster on large cohorts

## v1 Requirements

Requirements for this milestone. Each maps to roadmap phases.

### Benchmarking Framework (BENCH)

- [x] **BENCH-01**: Performance benchmark suite exists with pytest-benchmark covering inheritance analysis, genotype replacement, gene burden, scoring, and DataFrame I/O
- [x] **BENCH-02**: Synthetic data generators produce reproducible variant DataFrames at standard sizes (100, 1K, 10K, 50K variants) — benchmarks must use synthetic or public datasets only, no private cohort data in committed code
- [x] **BENCH-03**: Performance canary assertions detect 20%+ regressions without flaking on CI noise
- [x] **BENCH-04**: Ratio assertions compare vectorized vs sequential implementations within the same run (zero flakiness)
- [x] **BENCH-05**: Memory budget assertions via tracemalloc enforce per-function memory limits
- [x] **BENCH-06**: Benchmark results include extra metadata (n_variants, n_genes, peak_memory_mb)

### Quick Wins — Tier 1 (QWIN)

- [x] **QWIN-01**: Dead GT parsing loop removed from gene_burden.py (lines 220-249) with regression test confirming identical output
- [x] **QWIN-02**: `observed=True` added to all groupby calls across the codebase (12+ call sites)
- [x] **QWIN-03**: Temp file leak in gene_bed.py fixed (bed_path removed after merge completes)
- [x] **QWIN-04**: `gc.collect()` added between pipeline stages in runner.py to free intermediate memory

### DataFrame Optimization (DFOPT)

- [x] **DFOPT-01**: PyArrow engine used for hot-path `pd.read_csv()` calls with `engine="pyarrow"`
- [x] **DFOPT-02**: Categorical dtypes applied to low-cardinality columns (CHROM, IMPACT, FILTER, EFFECT, GENE) via optimization utility
- [x] **DFOPT-03**: iterrows() replaced with itertuples() in inheritance analysis Pass 2-3 (analyzer.py lines 128, 138)
- [x] **DFOPT-04**: DataFrame passed directly from pipeline context to Excel stage, eliminating redundant TSV re-read from disk

### Inheritance Analysis Optimization (INHER)

- [x] **INHER-01**: `df.apply(axis=1)` in inheritance Pass 1 (analyzer.py line 86) replaced with vectorized NumPy boolean mask operations, achieving 3.9-7.1x Pass 1 speedup (40-47% full analysis improvement)
- [x] **INHER-02**: Compound het result application (Pass 2) vectorized using column operations instead of iterrows loop
- [x] **INHER-03**: Full vectorization of `deduce_patterns_for_variant` using NumPy matrix operations on genotype columns, achieving 3.9-7.1x speedup on Pass 1
- [x] **INHER-04**: Vectorized compound het detection made the default implementation (comp_het_vectorized.py)

### Output Optimization (OUTPT)

- [x] **OUTPT-01**: xlsxwriter used for initial Excel write (2-5x faster), with openpyxl used only for finalization (hyperlinks, freeze panes, auto-filters)
- [x] **OUTPT-02**: GT column pre-parsed once at DataFrame load time into structured data, eliminating per-gene re-parsing

### Pipeline I/O Elimination (PIPEIO) — Supersedes PIPLN

- [x] **PIPEIO-01**: GenotypeReplacementStage skipped during pipeline — raw genotype columns flow directly to analysis stages
- [x] **PIPEIO-02**: `create_sample_columns_from_gt_intelligent()` handles raw SnpSift format without re-parsing replaced format
- [x] **PIPEIO-03**: Genotype replacement deferred to output time — TSV/Excel produces identical `"Sample(0/1)"` format
- [x] **PIPEIO-04**: SnpSift extractFields (2.7 hrs) replaced with `bcftools query` for VCF-to-TSV field extraction (C-based, 19.4x faster measured). Note: SnpSift filter stays — bcftools cannot replicate `na`, `has`, `ANN[ANY]` operators

### Parallelization & Chunking — Issue #62 (PARLZ)

**Note:** Scope revised based on Phase 7-11 transformations. Research showed memory pools and async I/O unnecessary after Phase 8 (82-84% memory reduction) and Phase 11 (I/O bottleneck elimination).

- [x] **PARLZ-01**: Dynamic content-aware chunking calculates optimal chunk sizes based on variant count, sample count, and available memory via pipeline-wide ResourceManager
- [x] **PARLZ-02**: Gene sorting (largest-first) for inheritance analysis load balancing, preventing straggler effects in parallel compound het detection (replaces work stealing — simpler, equally effective)
- [x] **PARLZ-03**: Pipeline-wide ResourceManager replaces inheritance-specific InheritanceMemoryManager with zero dead code (replaces memory pools — unnecessary after Phase 8 memory reduction)
- [x] **PARLZ-04**: Per-stage memory reporting at INFO level shows peak RSS and memory breakdown after pipeline completes (replaces async I/O — unnecessary after Phase 11 I/O elimination)

## Future Requirements

Deferred to later milestones. Tracked but not in current roadmap.

### CI Benchmark Integration

- **CIBENCH-01**: GitHub Actions benchmark.yml workflow with github-action-benchmark for historical tracking
- **CIBENCH-02**: PR regression alerts at 150% threshold with comment-on-alert
- **CIBENCH-03**: Performance charts at GitHub Pages URL

### Advanced Profiling

- **PROF-01**: pytest-memray integration for Linux CI (C-extension-aware memory tracking)
- **PROF-02**: Scalene line-level profiling integration for development-time analysis

### Architecture Cleanup

- **DEPR-01**: Deprecate classic pipeline mode (`pipeline.py`) — stage-based pipeline (`pipeline_core/`) should become the sole architecture. Classic mode adds maintenance burden, duplicates code paths, and prevents unified optimization. Migration: add deprecation warning when `--use-new-pipeline` is absent, then remove classic mode in a future milestone.

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Polars migration | 3-6 week rewrite, breaks scikit-learn/matplotlib integration, unacceptable risk for clinical tool |
| Free-threaded Python 3.13+ | numpy/pandas lack support for nogil builds |
| SIMD/cache-aware layouts | Diminishing returns, high complexity, pipeline is I/O and logic-bound not compute-bound |
| GPU acceleration | Insufficient numeric computation; pipeline is string-processing and I/O bound |
| SharedMemory for parallel stages | Complexity not justified for small gene-level DataFrames (1-20 variants each) |
| Modin | Drop-in parallelized pandas but doesn't fix root cause (inefficient apply/iterrows) |
| DuckDB/SQLite intermediates | Query layer overhead; pandas already in-memory and optimized |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| BENCH-01 | Phase 6 | Complete |
| BENCH-02 | Phase 6 | Complete |
| BENCH-03 | Phase 6 | Complete |
| BENCH-04 | Phase 6 | Complete |
| BENCH-05 | Phase 6 | Complete |
| BENCH-06 | Phase 6 | Complete |
| QWIN-01 | Phase 7 | Complete |
| QWIN-02 | Phase 7 | Complete |
| QWIN-03 | Phase 7 | Complete |
| QWIN-04 | Phase 7 | Complete |
| DFOPT-01 | Phase 8 | Complete |
| DFOPT-02 | Phase 8 | Complete |
| DFOPT-03 | Phase 8 | Complete |
| DFOPT-04 | Phase 8 | Complete |
| INHER-01 | Phase 9 | Complete |
| INHER-02 | Phase 9 | Complete |
| INHER-03 | Phase 9 | Complete |
| INHER-04 | Phase 9 | Complete |
| OUTPT-01 | Phase 10 | Complete |
| OUTPT-02 | Phase 10 | Complete |
| PIPEIO-01 | Phase 11 | Complete |
| PIPEIO-02 | Phase 11 | Complete |
| PIPEIO-03 | Phase 11 | Complete |
| PIPEIO-04 | Phase 11 | Complete |
| PARLZ-01 | Phase 12 | Complete |
| PARLZ-02 | Phase 12 | Complete |
| PARLZ-03 | Phase 12 | Complete |
| PARLZ-04 | Phase 12 | Complete |

**Coverage:**
- v1 requirements: 28 total (PIPLN-01/02 superseded by PIPEIO-01/02/03/04)
- Mapped to phases: 28
- Unmapped: 0 (100% coverage)

**Coverage Validation:**
- Phase 6: 6 requirements (BENCH-01 through BENCH-06)
- Phase 7: 4 requirements (QWIN-01 through QWIN-04)
- Phase 8: 4 requirements (DFOPT-01 through DFOPT-04)
- Phase 9: 4 requirements (INHER-01 through INHER-04)
- Phase 10: 2 requirements (OUTPT-01, OUTPT-02)
- Phase 11: 4 requirements (PIPEIO-01 through PIPEIO-04)
- Phase 12: 4 requirements (PARLZ-01 through PARLZ-04)
- **Total: 28/28 requirements mapped (100%)**

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-16 — All requirements complete (PARLZ scope revised, PIPEIO marked complete)*
