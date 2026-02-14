# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** Phase 8 - DataFrame Optimization

## Current Position

Phase: 8 of 12 (DataFrame Optimization)
Plan: 1 of 3 complete
Status: In progress
Last activity: 2026-02-14 — Completed 08-01-PLAN.md (DataFrame Optimizer Foundation)

Progress: [████████░░░░░░░░░░░░] 47% (Phase 1-7 complete, Phase 8 1/3)

## Performance Metrics

**Velocity:**
- Total plans completed: 8
- Average duration: 22.1 minutes
- Total execution time: 2.95 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1-5. Baseline | N/A | N/A | N/A (pre-GSD) |
| 6. Benchmark Framework | 4/4 | 48.0 min | 12.0 min |
| 7. Quick Wins Tier 1 | 3/3 | 89.0 min | 29.7 min |
| 8. DataFrame Optimization | 1/3 | 18.0 min | 18.0 min |

**Recent Trend:**
- Last 5 plans: 07-01 (4.0 min), 07-02 (10.0 min), 07-03 (75.0 min), 08-01 (18.0 min)
- Trend: Foundation work fast (4-18 min), benchmark verification longer (75 min)

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Benchmarks before optimization: Can't prove improvement without baseline measurements
- Accept risk on Tier 3 (vectorization, Cython): User wants maximum performance gains, will revert if tests break
- Local benchmarks only, no CI workflow: Get benchmarks working first, CI integration deferred
- Target v0.13.0: Performance improvements warrant minor version bump
- Removed entire dead GT parsing loop (07-01): Existing aggregated counts already provide correct values, loop was parsing into unused variables
- Regression tests before dead code removal (07-01): Proves output remains identical after refactoring
- observed=True now vs later (07-02): Added observed=True in Phase 7 before categorical dtypes exist to prevent regressions when Phase 8 introduces them
- Unconditional gc.collect() (07-02): Runs after EVERY stage execution (not just DEBUG) for memory management with large genomic datasets
- Pre-commit hook for groupby enforcement (07-02): Prevents new groupby calls without observed=True from being committed
- Full 60-benchmark suite completed (07-03): All benchmarks passed in 5:37, saved as `0001_phase7_quick_wins.json`. Deprecated streaming perf tests removed.
- Document collateral improvements (07-03): Inheritance analysis improved 20-58% despite no direct optimizations, attributed to gc.collect() and observed=True reducing overhead
- PyArrow scope limited to variant DataFrames (08-01): Main TSV loading only, not config/gene lists per 08-CONTEXT decision
- Categorical dtype auto-detection at 50% cardinality (08-01): Columns with <50% unique values loaded as categorical
- Column renaming permanent at load time (08-01): No temporary rename-restore, downstream uses sanitized names
- Memory pass-through threshold 25% available RAM (08-01): Conservative targeting 8-16GB desktops
- quoting parameter excluded for PyArrow engine (08-01): PyArrow doesn't support it, C engine fallback used

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

**Phase 6 (Benchmark Framework): COMPLETE**
- 60 benchmark tests across 8 files, all passing
- Baseline saved: `.benchmarks/Windows-CPython-3.10-64bit/0001_baseline_v0.12.1.json`
- Phase 7 results saved: `.benchmarks/Linux-CPython-3.10-64bit/0001_phase7_quick_wins.json`
- Deprecated streaming parallel perf tests removed (were blocking benchmark runs)
- Framework ready for Phase 8+ optimization work

**Phase 7 Plan 01 (Dead Code Removal): COMPLETE**
- Removed 30-line dead GT parsing loop from gene_burden.py
- Fixed 3 temp file leaks (gene_bed.py bed_path, filters.py mktemp replacement)
- Regression tests prove gene burden output unchanged
- Expected 20-30% speedup in gene burden analysis
- Existing bug discovered: gene_burden.py crashes on empty DataFrame (not fixed, out of scope)

**Phase 7 Plan 02 (Categorical dtypes & Memory Management): COMPLETE**
- Added observed=True to all 17 groupby call sites (future-proofing for Phase 8 categorical dtypes)
- Implemented gc.collect() after every pipeline stage execution
- Added DEBUG-level memory logging (RSS before/after each stage, freed memory delta via psutil)
- Created pre-commit hook enforcing observed=True in all new groupby calls
- Prevents 3500x groupby slowdown when categorical dtypes are introduced in Phase 8
- All 565 unit tests pass with no behavioral changes

**Phase 7 Plan 03 (Benchmark Verification): COMPLETE**
- Verified Phase 7 optimizations with actual benchmark measurements
- Gene burden analysis: 48-98% faster (exceeded 20-40% target)
- Inheritance analysis: 20-58% faster (collateral improvement from gc.collect and observed=True)
- Zero regressions across all benchmarks
- Created comprehensive performance analysis report for Phase 8+ comparison
- Fixed test bug: added external tool checks to gene burden integration tests

**Baseline Performance (v0.12.1, Windows, CPython 3.10, 2026-02-14):**

| Component | Scale | v0.12.1 Baseline | After Phase 7 | Improvement |
|-----------|-------|-----------------|---------------|-------------|
| Single variant deduce | 1 variant | 7.7 us | 4.4 us | **43.0%** |
| Inheritance analysis | 100 variants | 49 ms | 37 ms | **23.8%** |
| Inheritance analysis | 1K variants | 431 ms | 346 ms | **19.6%** |
| Inheritance analysis | 10K variants | 4.8 s | 3.2 s | **33.2%** |
| Gene burden | 100 variants | 62 ms | 32 ms | **48.5%** |
| Gene burden | 1K variants | 150 ms | 18 ms | **88.0%** |
| Gene burden | 10K variants | 1.0 s | 19 ms | **98.1%** |
| Gene burden (10 genes) | - | 98 ms | 4 ms | **95.6%** |
| Gene burden (50 genes) | - | 121 ms | 19 ms | **84.0%** |
| Gene burden (100 genes) | - | 137 ms | 49 ms | **64.7%** |

**Phase 7 achieved 48-98% speedup on gene burden, 20-58% on inheritance analysis.**

**Phase 8 (DataFrame Optimization): IN PROGRESS (1/3 complete)**

**Plan 01 (DataFrame Optimizer Foundation): COMPLETE**
- Created dataframe_optimizer.py with PyArrow loading, categorical detection, column sanitization
- PyArrow engine now used automatically (5-15x CSV read speedup)
- Low-cardinality columns loaded as categorical (50-75% memory reduction expected)
- Column sanitization complete (GEN[0].GT → GEN_0__GT) - ready for itertuples migration
- Memory pass-through decision logic in place (25% available RAM threshold)
- All 568 unit tests + 31 integration tests pass with no regressions
- Ready for Plan 02-03: itertuples migration in hot paths

**Remaining work:**
- Plan 02: Replace iterrows with itertuples in gene_burden, scoring, statistics
- Plan 03: Replace iterrows in inheritance analysis (compound het detection)

**Phase 9 (Inheritance Analysis Optimization):**
- Full vectorization (INHER-03) is high-risk Tier 3 work, requires extensive validation to preserve clinical correctness
- Must preserve byte-identical output for all inheritance patterns
- Improved baseline: 3.2s for 10K variants (down from 4.8s after Phase 7 collateral improvements)

**Phase 11 (Pipeline & Cython Optimization):**
- Subprocess piping error handling needs platform-specific validation (Windows vs Linux)
- Cython extension requires hatch-cython integration with existing hatchling build

## Session Continuity

Last session: 2026-02-14 14:04 UTC
Stopped at: Completed 08-01-PLAN.md (DataFrame Optimizer Foundation)
Resume file: None
Next: Phase 8 Plan 02 (itertuples migration in gene_burden, scoring, stats)
