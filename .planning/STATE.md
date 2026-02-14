# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** Phase 7 - Quick Wins Tier 1

## Current Position

Phase: 7 of 12 (Quick Wins - Tier 1)
Plan: 2 of 3 complete
Status: In progress
Last activity: 2026-02-14 — Completed 07-02-PLAN.md (Categorical dtypes & Memory Management)

Progress: [██████░░░░░░░░░░░░░░] 38% (Phase 1-6 complete, 07-01 and 07-02 complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 6
- Average duration: 10.3 minutes
- Total execution time: 1.03 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1-5. Baseline | N/A | N/A | N/A (pre-GSD) |
| 6. Benchmark Framework | 4/4 | 48.0 min | 12.0 min |
| 7. Quick Wins Tier 1 | 2/3 | 14.0 min | 7.0 min |

**Recent Trend:**
- Last 5 plans: 06-02 (12.0 min), 06-04 (26.0 min), 07-01 (4.0 min), 07-02 (10.0 min)
- Trend: Quick wins are fast (4-10 min), benchmark tests are slower (12-26 min)

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

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

**Phase 6 (Benchmark Framework): COMPLETE**
- 60 benchmark tests across 8 files, all passing
- Baseline saved: `.benchmarks/Windows-CPython-3.10-64bit/0001_baseline_v0.12.1.json`
- Framework ready for Phase 7+ optimization work

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

**Baseline Performance (v0.12.1, Windows, CPython 3.10, 2026-02-14):**

| Component | Scale | Mean | Notes |
|-----------|-------|------|-------|
| Single variant deduce | 1 variant | 7.7 us | Inner loop function |
| Inheritance analysis | 100 variants | 49 ms | Full orchestrator |
| Inheritance analysis | 1K variants | 431 ms | |
| Inheritance analysis | 10K variants | 4.8 s | Main optimization target |
| Full inheritance cohort | 5K×100 samples | 2.3 s | Macro benchmark |
| Comp het vectorized | 1K variants | 2.1 ms | Per-gene |
| Comp het original | 1K variants | 1.6 ms | Vectorized NOT faster at small scale |
| Gene burden | 100 variants | 62 ms | |
| Gene burden | 1K variants | 150 ms | |
| Gene burden | 10K variants | 1.0 s | |
| Gene burden full | 5K×200 samples | 745 ms | Macro benchmark |
| Genotype replacement vec | 100 variants | 297 ms | File-based, includes I/O |
| Genotype replacement vec | 1K variants | 2.7 s | |
| Genotype replacement vec | 10K variants | 21.9 s | Very slow — optimization target |
| Genotype replacement seq | 1K variants | 11.5 ms | Line-based iterator |
| CSV read | 50K variants | 57 ms | pandas default engine |
| PyArrow read | 50K variants | 35 ms | 1.6x faster than default |
| CSV write | 50K variants | 156 ms | |

**Phase 8 (DataFrame Optimization): READY**
- PyArrow engine + pandas 3.0 string dtype changes require careful testing to avoid `pd.NA` comparison breakage
- All groupby calls now have observed=True (Phase 7 Plan 02 complete) - no risk of 3500x slowdown when categorical dtypes are introduced
- Pre-commit hook enforces observed=True in new groupby calls (prevents regressions)

**Phase 9 (Inheritance Analysis Optimization):**
- Full vectorization (INHER-03) is high-risk Tier 3 work, requires extensive validation to preserve clinical correctness
- Must preserve byte-identical output for all inheritance patterns

**Phase 11 (Pipeline & Cython Optimization):**
- Subprocess piping error handling needs platform-specific validation (Windows vs Linux)
- Cython extension requires hatch-cython integration with existing hatchling build

## Session Continuity

Last session: 2026-02-14 10:37 UTC
Stopped at: Completed 07-02-PLAN.md
Resume file: None
Next: 07-03-PLAN.md (Cython CI/CD Integration)
