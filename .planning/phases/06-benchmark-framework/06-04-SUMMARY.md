---
phase: 06-benchmark-framework
plan: 04
subsystem: testing
tags: [pytest, pytest-benchmark, performance, benchmarking, synthetic-data]

# Dependency graph
requires:
  - phase: 06-02
    provides: Component-level benchmarks (inheritance, comp_het, genotype replacement, gene burden, scoring, DataFrame I/O)
  - phase: 06-03
    provides: Ratio assertions and memory budgets
provides:
  - Result diff helper for comparing pytest-benchmark JSON files
  - Macro-level pipeline benchmarks (full inheritance, gene burden at cohort scale)
  - Verified end-to-end benchmark framework workflow
  - Confirmed all BENCH-01 through BENCH-06 requirements met
affects: [07-optimize-inheritance, 08-optimize-genotype-replacement, 09-optimize-gene-burden, 10-optimize-dataframe-io, 11-cython-critical-paths]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Result diff utility: color-coded comparison tables with ANSI codes"
    - "Macro benchmarks use pedantic mode (rounds=3, iterations=1) for expensive operations"
    - "GT transformation helpers convert comma-separated to post-replacement format"

key-files:
  created:
    - tests/performance/helpers/result_diff.py
    - tests/performance/benchmark_pipeline_macro.py
  modified: []

key-decisions:
  - "Result diff helper uses simple color-coded table output (green=faster, red=slower) for sprint tooling"
  - "Macro benchmarks test Python code paths only (no external tools) per CONTEXT.md requirements"
  - "Large cohort test (10K x 500) marked @pytest.mark.slow for skipability during quick iterations"
  - "Benchmark comparison uses within-run ratio assertions to eliminate flakiness from system variance"

patterns-established:
  - "Macro benchmarks use _transform_to_postprocessed_gt helper to convert synthetic GT data"
  - "All benchmarks include extra_info metadata: level (micro/meso/macro), component, scale parameters"
  - "Result diff helper loads pytest-benchmark JSON, matches by name, shows ±% change with status"

# Metrics
duration: 26min
completed: 2026-02-14
---

# Phase 6 Plan 4: Result Diff Helper and Macro Pipeline Benchmarks Summary

**Result diff utility with color-coded comparison tables and macro pipeline benchmarks exercising full inheritance/gene burden at 100-500 sample cohort scale**

## Performance

- **Duration:** 26 min
- **Started:** 2026-02-14T07:24:26Z
- **Completed:** 2026-02-14T07:50:26Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Result diff helper compares pytest-benchmark JSON files with color-coded output (green=faster, red=slower, yellow=new/removed)
- Macro benchmarks test full inheritance analysis (5K-10K variants, 100-500 samples) and gene burden (5K variants, 200 samples, 100 genes)
- End-to-end benchmark suite verified: all 60 tests pass, save/compare workflow works, all BENCH requirements met
- Framework ready for Phase 7+ optimization work with baseline measurements

## Task Commits

Each task was committed atomically:

1. **Task 1: Create result diff helper and macro pipeline benchmarks** - `24bcce4` (feat)

Task 2 was verification-only (no new files).

**Plan metadata:** (to be added in final commit)

## Files Created/Modified
- `tests/performance/helpers/result_diff.py` - Compares pytest-benchmark JSON files, prints color-coded table showing faster/slower/same benchmarks with percentage changes
- `tests/performance/benchmark_pipeline_macro.py` - Macro-level benchmarks for full inheritance analysis (cohort scale: 100-500 samples) and gene burden (200 samples, 100 genes)

## Decisions Made

**Result diff simplicity:** Kept the helper dead simple per CONTEXT.md ("Simple diff/comparison helper: load two JSON result files, show what got faster/slower"). No histograms, no fancy formatting, just clear table output. This is sprint tooling, not permanent infrastructure.

**Macro benchmark scope:** Focused on Python code paths only (inheritance analyzer, gene burden) without external tools. The inheritance analysis orchestrator accounts for 40-60% of total pipeline time, making it the primary optimization target.

**Large cohort test marked slow:** The 10K x 500 sample benchmark is marked `@pytest.mark.slow` so it can be skipped with `pytest -m "performance and not slow"` during quick iterations. Full suite still runs it when needed for comprehensive baseline measurements.

**GT format transformation:** Created `_transform_to_postprocessed_gt` helper to convert synthetic GT data (comma-separated) to post-replacement format (SAMPLE(genotype);...) that the inheritance analyzer expects. This matches the format after genotype replacement stage in the real pipeline.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**pytest-benchmark compare requires exact filename:** The `--benchmark-compare` flag needs the exact filename without path, not just the save name. However, the result diff helper provides better functionality anyway (color coding, clear percentage changes, new/removed tracking), so this becomes the primary comparison tool.

**Background task output files empty:** When running long benchmarks in background, output files weren't being written. Switched to foreground execution with extended timeout for reliability. Background execution works for status tracking but not for immediate output capture.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Benchmark framework complete and verified:**
- BENCH-01: Coverage ✓ (inheritance, comp_het, genotype replacement, gene burden, scoring, DataFrame I/O, macro)
- BENCH-02: Synthetic data generators ✓ (reproducible, seeded)
- BENCH-03: Regression detection ✓ (--benchmark-compare-fail=mean:20%)
- BENCH-04: Ratio assertions ✓ (zero-flakiness within-run comparisons)
- BENCH-05: Memory budgets ✓ (tracemalloc with warning-only enforcement)
- BENCH-06: Metadata ✓ (all 60 benchmarks include extra_info)

**Ready for Phase 7 (Optimize Inheritance):**
- Baseline measurements can be captured with: `pytest tests/performance/benchmark_inheritance.py tests/performance/benchmark_comp_het.py --benchmark-save=baseline_inheritance`
- Result diff helper ready to compare before/after optimization runs
- Ratio assertions will catch regressions in vectorized vs sequential performance
- Memory budgets will detect memory issues from optimization changes

**Framework characteristics:**
- 60 total benchmark tests across 9 files
- Tests run in 60-120 seconds with `--benchmark-disable` (correctness verification)
- Full benchmark suite takes longer (with timing overhead) but can be filtered by scale
- All tests pass, framework is stable and reproducible

**No blockers.** Optimization phases can begin immediately.

---
*Phase: 06-benchmark-framework*
*Completed: 2026-02-14*
