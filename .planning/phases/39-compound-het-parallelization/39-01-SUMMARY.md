---
phase: 39-compound-het-parallelization
plan: 01
subsystem: testing
tags: [benchmark, performance, compound-het, parallelization, synthetic-data, inheritance]

# Dependency graph
requires:
  - phase: 38-codebase-cleanup
    provides: clean codebase starting point for optimization work
provides:
  - Standalone synthetic benchmark script for compound het parallelization
  - Pre-optimization baseline timing (sequential: ~2.3s for 1000 genes / 10000 variants)
  - Confirmed GIL contention problem: all parallel configurations slower than sequential
affects:
  - 39-02 (optimization plan uses this baseline to validate speedup)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Standalone benchmark scripts (not pytest) for manual timing review in tests/performance/
    - sys.path.insert + noqa: I001 pattern for runnable scripts outside package install

key-files:
  created:
    - tests/performance/standalone_bench_comp_het_parallel.py
    - .planning/phases/39-compound-het-parallelization/39-BASELINE.txt
  modified: []

key-decisions:
  - "Benchmark measures total wall-time (all 3 passes), not isolated Pass 2, since Pass 1+3 are fast relative to Pass 2 for large datasets"
  - "Baseline confirms GIL contention hypothesis: parallel always slower than sequential pre-optimization"
  - "Synthetic data uses per-sample columns (CHILD_001, FATHER_001, etc.) not GEN_N__GT sanitized format - analyze_inheritance_parallel expects sample IDs directly as column names"

patterns-established:
  - "Standalone benchmark pattern: sys.path manipulation + noqa: I001 for import ordering, no assertions, exit 0 always"

# Metrics
duration: 17min
completed: 2026-02-26
---

# Phase 39 Plan 01: Compound Het Parallelization Benchmark Summary

**Standalone synthetic benchmark measuring GIL-bound compound het Pass 2 performance; baseline confirms parallel path is 0.21x-0.92x slower than sequential due to GIL contention**

## Performance

- **Duration:** 17 min
- **Started:** 2026-02-26T21:44:29Z
- **Completed:** 2026-02-26T22:02:23Z
- **Tasks:** 2
- **Files modified:** 2 created

## Accomplishments

- Created `tests/performance/standalone_bench_comp_het_parallel.py` — self-contained benchmark generating 1000 genes / 10000 variants with realistic skewed distribution (few large genes, many small), 2 trios (6 samples), sweeping 5 thresholds x 4 worker counts
- Captured pre-optimization baseline to `.planning/phases/39-compound-het-parallelization/39-BASELINE.txt` — sequential median 2.337s, all parallel configs slower (worst: 0.21x at workers=16, threshold=50 and workers=2, threshold=100)
- Confirmed codebase is clean before optimization: all 2228 tests pass, `make ci-check` passes

## Task Commits

Each task was committed atomically:

1. **Task 1: Create standalone synthetic benchmark script** - `d0dc17a` (feat)
2. **Task 2: Capture baseline timing and verify existing tests pass** - `f3e0199` (chore)

**Plan metadata:** (created after this summary)

## Files Created/Modified

- `tests/performance/standalone_bench_comp_het_parallel.py` — Standalone benchmark script; generates synthetic data, sweeps configurations, prints timing table to stdout; no assertions; exit 0
- `.planning/phases/39-compound-het-parallelization/39-BASELINE.txt` — Pre-optimization timing baseline; sequential ~2.3s median for 1000 genes/10000 variants on WSL2 with 16 CPUs

## Decisions Made

- Benchmark measures total `analyze_inheritance_parallel()` wall-time rather than isolated Pass 2, because Pass 1 (vectorized deduction) and Pass 3 (prioritization) are fast relative to Pass 2 for large datasets, making total time a reasonable proxy
- Synthetic data uses sample IDs directly as DataFrame column names (e.g., `CHILD_001`, `FATHER_001`), matching what `analyze_inheritance_parallel()` expects — not the `GEN_N__GT` sanitized format used in pipeline output
- Used `# noqa: I001` on the first third-party import after sys.path manipulation to satisfy ruff's import sorting rules while keeping the script self-contained and runnable without installation

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Ruff I001 import sorting violation in standalone script**

- **Found during:** Task 2 (CI check run)
- **Issue:** sys.path.insert before third-party imports triggers ruff I001 (unsorted import block). E402 noqa doesn't apply (that rule is disabled). Need I001 noqa.
- **Fix:** Added `# noqa: I001` on first post-path import (`import numpy as np`), then ran `ruff format` to fix whitespace
- **Files modified:** `tests/performance/standalone_bench_comp_het_parallel.py`
- **Verification:** `make ci-check` passes with all lint, format, typecheck, and test stages green
- **Committed in:** `d0dc17a` (Task 1 commit, file was amended before final commit)

---

**Total deviations:** 1 auto-fixed (1 bug - lint compliance)
**Impact on plan:** Necessary for CI compliance. No scope creep.

## Issues Encountered

None - plan executed as specified. Benchmark ran successfully on first attempt after lint fix.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Baseline timing captured and committed: Plan 02 can reference `39-BASELINE.txt` to validate any speedup
- Key finding for Plan 02: All parallel configurations are slower than sequential. The GIL contention hypothesis from CONTEXT.md is confirmed empirically. Optimization target is to make parallel >= sequential, then faster.
- Sequential baseline (median 2.337s) establishes the floor; any optimization yielding >1.0x ratio compared to sequential will represent genuine improvement
- All 2228 existing tests pass — clean starting point for optimization work

---
*Phase: 39-compound-het-parallelization*
*Completed: 2026-02-26*
