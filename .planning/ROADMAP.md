# Roadmap: VariantCentrifuge

## Milestones

- ✅ **v0.12.1 Baseline** - Phases 1-5 (shipped 2026-02-14, pre-GSD)
- 🚧 **v0.13.0 Performance Optimization** - Phases 6-12 (in progress)

## Phases

<details>
<summary>✅ v0.12.1 Baseline (Phases 1-5) - SHIPPED 2026-02-14</summary>

**Milestone Goal:** Full-featured variant analysis pipeline with two architectures, 40+ stages, inheritance analysis, gene burden, scoring, multi-format output. All 30 historical issues resolved. 1035 tests passing. CI/CD with Docker, docs, and multi-platform testing.

**Note:** Phases 1-5 represent completed work before GSD tracking. Details available in git history.

</details>

### 🚧 v0.13.0 Performance Optimization (In Progress)

**Milestone Goal:** 3-4x full pipeline speedup on large cohorts through benchmarking infrastructure and systematic optimization of DataFrame operations, inheritance analysis, and output generation.

#### ✅ Phase 6: Benchmark Framework (Complete)
**Goal**: Performance benchmarking infrastructure exists with synthetic data and regression detection
**Depends on**: Nothing (first phase of milestone)
**Requirements**: BENCH-01, BENCH-02, BENCH-03, BENCH-04, BENCH-05, BENCH-06
**Success Criteria** (what must be TRUE):
  1. ✅ Benchmark suite runs via `pytest tests/performance/` covering inheritance, genotype replacement, gene burden, scoring, and DataFrame I/O
  2. ✅ Synthetic data generators produce reproducible DataFrames at 100, 1K, 10K, 50K variants without referencing private cohort names
  3. ✅ Ratio assertions compare vectorized vs sequential implementations within same run with zero CI flakiness
  4. ✅ Performance canary assertions detect 20%+ regressions and fail tests when exceeded
  5. ✅ Memory budget assertions enforce per-function limits via tracemalloc
**Plans:** 4 plans (4/4 complete)
**Verified:** 5/5 must-haves passed

Plans:
- [x] 06-01-PLAN.md -- Foundation: pytest-benchmark dep, synthetic data generators, memory helpers, conftest fixtures
- [x] 06-02-PLAN.md -- Component benchmarks: inheritance, comp_het, genotype replacement, gene burden, scoring, DataFrame I/O
- [x] 06-03-PLAN.md -- Ratio assertions (vectorized vs sequential) and memory budget enforcement via tracemalloc
- [x] 06-04-PLAN.md -- Macro pipeline benchmarks, result diff helper, end-to-end suite verification

#### Phase 7: Quick Wins - Tier 1
**Goal**: Immediate 30-40% speedup through zero-risk standard optimizations
**Depends on**: Phase 6 (benchmarks prove improvement)
**Requirements**: QWIN-01, QWIN-02, QWIN-03, QWIN-04
**Success Criteria** (what must be TRUE):
  1. Dead GT parsing loop (gene_burden.py:220-249) removed with regression test confirming identical output
  2. All 12+ groupby call sites use `observed=True` to prevent categorical dtype slowdown
  3. Temp file leak in gene_bed.py fixed (bed_path removed after merge)
  4. Memory freed between pipeline stages via `gc.collect()` in runner.py
  5. Benchmarks show 30-40% speedup on gene burden analysis with identical output
**Plans**: TBD

Plans:
- [ ] 07-01: [TBD during planning]

#### Phase 8: DataFrame Optimization
**Goal**: 50-70% memory reduction and 2-3x I/O speedup through optimal DataFrame loading
**Depends on**: Phase 7 (observed=True must exist before categorical dtypes)
**Requirements**: DFOPT-01, DFOPT-02, DFOPT-03, DFOPT-04
**Success Criteria** (what must be TRUE):
  1. PyArrow engine used for hot-path CSV reads with `engine="pyarrow"`
  2. Categorical dtypes applied to low-cardinality columns (CHROM, IMPACT, FILTER, EFFECT, GENE) without breaking comparisons
  3. iterrows replaced with itertuples in inheritance Pass 2-3 achieving 10-13x speedup
  4. DataFrame passed directly from pipeline context to Excel stage without redundant disk read
  5. Memory profiling shows 50-70% reduction, benchmarks show 2-3x I/O speedup
**Plans**: TBD

Plans:
- [ ] 08-01: [TBD during planning]

#### Phase 9: Inheritance Analysis Optimization
**Goal**: 10-100x speedup on inheritance analysis (40-60% of total pipeline time)
**Depends on**: Phase 8 (categorical dtypes enable vectorization)
**Requirements**: INHER-01, INHER-02, INHER-03, INHER-04
**Success Criteria** (what must be TRUE):
  1. `df.apply(axis=1)` in Pass 1 replaced with itertuples or vectorized operations achieving 10-100x speedup
  2. Compound het result application (Pass 2) vectorized using column operations instead of iterrows
  3. Full vectorization of deduce_patterns_for_variant implemented using NumPy matrix operations on genotype columns
  4. Vectorized compound het detection is default implementation (comp_het_vectorized.py)
  5. Benchmarks show 10-100x Pass 1 speedup with byte-identical output preserved
**Plans**: TBD

Plans:
- [ ] 09-01: [TBD during planning]

#### Phase 10: Output Optimization
**Goal**: 2-5x faster Excel generation and eliminate redundant GT parsing
**Depends on**: Phase 8 (DataFrame in-memory pass-through must exist)
**Requirements**: OUTPT-01, OUTPT-02
**Success Criteria** (what must be TRUE):
  1. xlsxwriter used for initial Excel write with openpyxl for finalization (hyperlinks, freeze panes, auto-filters)
  2. GT column pre-parsed once at DataFrame load time into structured data
  3. Per-gene GT re-parsing eliminated across gene burden and statistics stages
  4. Benchmarks show 2-5x Excel generation speedup
  5. Output Excel file has hyperlinks, freeze panes, and auto-filters verified by tests
**Plans**: TBD

Plans:
- [ ] 10-01: [TBD during planning]

#### Phase 11: Pipeline & Cython Optimization
**Goal**: Eliminate intermediate temp files and accelerate genotype replacement kernel
**Depends on**: Phase 9 (inheritance must be optimized first, higher impact)
**Requirements**: PIPLN-01, PIPLN-02
**Success Criteria** (what must be TRUE):
  1. Pipe fusion implemented for SnpSift filter → bgzip chain without intermediate temp VCF
  2. Both process return codes checked with proper error propagation
  3. Cython extension created for genotype replacement kernel (hot path string processing)
  4. Benchmarks show 10-30s saved on filtering, 2-3x speedup on genotype replacement
  5. Cross-platform tests pass on Windows and Linux
**Plans**: TBD

Plans:
- [ ] 11-01: [TBD during planning]

#### Phase 12: Parallelization & Chunking
**Goal**: Dynamic work distribution and memory-efficient processing at scale
**Depends on**: Phase 11 (all other optimizations complete first)
**Requirements**: PARLZ-01, PARLZ-02, PARLZ-03, PARLZ-04
**Success Criteria** (what must be TRUE):
  1. Dynamic chunking calculates optimal sizes based on variant count, sample count, and genotype density
  2. Adaptive work stealing redistributes work from large gene groups to idle workers during inheritance analysis
  3. Memory pool management reduces allocation overhead with reusable pre-allocated buffers
  4. Asynchronous I/O with memory mapping overlaps reads with processing for I/O-bound workloads
  5. Benchmarks show improved scaling on large cohorts (5000+ samples) with no memory regressions
**Plans**: TBD

Plans:
- [ ] 12-01: [TBD during planning]

## Progress

**Execution Order:** Phases execute in numeric order: 6 → 7 → 8 → 9 → 10 → 11 → 12

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1-5. Baseline | v0.12.1 | N/A | Complete | 2026-02-14 |
| 6. Benchmark Framework | v0.13.0 | 4/4 | Complete | 2026-02-14 |
| 7. Quick Wins - Tier 1 | v0.13.0 | 0/TBD | Not started | - |
| 8. DataFrame Optimization | v0.13.0 | 0/TBD | Not started | - |
| 9. Inheritance Analysis Optimization | v0.13.0 | 0/TBD | Not started | - |
| 10. Output Optimization | v0.13.0 | 0/TBD | Not started | - |
| 11. Pipeline & Cython Optimization | v0.13.0 | 0/TBD | Not started | - |
| 12. Parallelization & Chunking | v0.13.0 | 0/TBD | Not started | - |
