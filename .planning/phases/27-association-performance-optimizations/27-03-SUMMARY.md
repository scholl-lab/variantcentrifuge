---
phase: 27-association-performance-optimizations
plan: "03"
subsystem: association
tags: [concurrent.futures, ProcessPoolExecutor, parallelization, BLAS, numpy, scipy]

# Dependency graph
requires:
  - phase: 27-01
    provides: parallel_safe=True on all Python-backend test classes (Fisher, LogisticBurden, LinearBurden, PurePythonSKAT)
  - phase: 27-02
    provides: AssociationConfig.association_workers field and --association-workers CLI arg

provides:
  - ProcessPoolExecutor-based gene-level parallelization in AssociationEngine.run_all()
  - Module-level _run_gene_worker() worker function for subprocess dispatch
  - Module-level _worker_initializer() to prevent BLAS thread oversubscription
  - Null-model pre-fitting via sequential first-gene before parallel dispatch
  - Sequential fallback with warning when any test has parallel_safe=False
  - Unit tests verifying parallel/sequential result identity and fallback behavior

affects:
  - future performance benchmarking phases
  - association analysis runtime for large multi-gene panels

# Tech tracking
tech-stack:
  added: [concurrent.futures.ProcessPoolExecutor, pickle]
  patterns:
    - "Module-level worker functions required for ProcessPoolExecutor pickle compatibility"
    - "First-gene sequential execution triggers lazy null model fitting before parallel dispatch"
    - "BLAS thread oversubscription prevention via initializer in worker processes"

key-files:
  created:
    - tests/unit/test_parallel_association.py
  modified:
    - variantcentrifuge/association/engine.py

key-decisions:
  - "PERF-03: First gene runs sequentially to trigger lazy null model fitting in SKAT/COAST tests before pickle-based dispatch to workers"
  - "PERF-04: Worker processes set OPENBLAS/MKL/OMP NUM_THREADS=1 to prevent BLAS oversubscription in ProcessPoolExecutor"
  - "PERF-05: parallel path only taken when n_workers != 1 AND all tests parallel_safe AND len(sorted_data) > 1"
  - "PERF-06: Worker under-provisioning guard: actual_workers = max(1, len(remaining)//2) when remaining < actual_workers * 2"

patterns-established:
  - "Parallel dispatch via module-level function + pickle of test instances with pre-fitted null models"
  - "Sequential fallback log warning pattern: association_workers=%d requested but not all tests are parallel_safe"

# Metrics
duration: 10min
completed: 2026-02-22
---

# Phase 27 Plan 03: Gene-Level Parallelization Summary

**ProcessPoolExecutor gene-level parallelization in AssociationEngine.run_all() with BLAS oversubscription prevention, lazy null model pre-fitting, and parallel/sequential result identity verified by unit tests**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-22T20:27:09Z
- **Completed:** 2026-02-22T20:37:50Z
- **Tasks:** 2/2
- **Files modified:** 2 (1 created, 1 modified)

## Accomplishments

- Added `_worker_initializer()` and `_run_gene_worker()` as module-level functions (required for ProcessPoolExecutor pickle compatibility)
- Modified `run_all()` to dispatch genes to worker processes when `association_workers > 1` and all tests are `parallel_safe=True`
- First gene always runs sequentially before parallel dispatch to trigger lazy null model fitting in SKAT/COAST tests
- Non-parallel-safe tests (R backend) trigger sequential fallback with explicit warning message
- Unit tests confirm parallel/sequential result identity, fallback behavior, single-gene guard, and BLAS env var setting

## Task Commits

1. **Task 1: ProcessPoolExecutor parallel gene loop in engine.py** - `3acec44` (feat)
2. **Task 2: Parallel execution tests** - `210cd9a` (test)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `variantcentrifuge/association/engine.py` - Added `_worker_initializer()`, `_run_gene_worker()`, and parallel path in `run_all()` with sequential fallback logic
- `tests/unit/test_parallel_association.py` - 4 unit tests: parallel/sequential identity, fallback warning, single-gene guard, worker initializer env vars

## Decisions Made

- **PERF-03:** First gene runs sequentially before parallel dispatch. This is necessary because SKAT and COAST tests fit their null models lazily on the first `run()` call. If the null model is not fitted before pickling, each worker would refit it — negating the performance benefit and potentially causing inconsistency.
- **PERF-04:** Worker initializer sets OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, and OMP_NUM_THREADS to "1". This prevents N_workers * BLAS_threads oversubscription which would cause CPU thrashing on multi-core machines.
- **PERF-05:** Three conditions must all be true for parallel path: `n_workers != 1`, `all_parallel_safe`, and `len(sorted_data) > 1`. Single-gene panels skip parallel overhead entirely.
- **PERF-06:** Worker count under-provisioning guard prevents spawning more workers than useful: `actual_workers = max(1, len(remaining) // 2)` when remaining genes < `actual_workers * 2`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Type annotation quotes removed for UP037 lint compliance**

- **Found during:** Task 1 (engine.py modification)
- **Issue:** String-quoted type annotations `"tuple[str, dict, bytes, AssociationConfig]"` trigger ruff UP037 (remove-quotes) even with `from __future__ import annotations` present
- **Fix:** Removed quotes from `_run_gene_worker()` parameter and return type annotations
- **Files modified:** `variantcentrifuge/association/engine.py`
- **Verification:** `make lint` passes
- **Committed in:** `3acec44` (Task 1 commit, formatting pass)

**2. [Rule 1 - Bug] Ruff format reformatted log f-string and list comprehension**

- **Found during:** Task 2 verification (`make ci-check`)
- **Issue:** f-string continuation `f"Gene {first_gene} | ... "` and multi-line list comprehension `[(gd.get("GENE",...) for gd in remaining]` were reformatted by ruff
- **Fix:** Applied `python -m ruff format variantcentrifuge/association/engine.py`
- **Files modified:** `variantcentrifuge/association/engine.py`
- **Verification:** `make format-check` passes
- **Committed in:** `210cd9a` (Task 2 commit staging)

---

**Total deviations:** 2 auto-fixed (both Rule 1 - lint/format compliance)
**Impact on plan:** Both required for CI compliance. No functional changes, no scope creep.

## Issues Encountered

None — plan executed exactly as specified. The parallel path correctness was confirmed by `pd.testing.assert_frame_equal` in the unit tests.

## Next Phase Readiness

- Phase 27 is complete (all 3 plans done: GL quadrature, workers CLI plumbing, parallelization)
- Milestone v0.15.0 is code-complete — all phases 18-27 implemented
- The association engine now supports ~Nx wall-clock speedup for multi-gene panels when using Python-backend tests
- R-backend tests retain sequential execution path (parallel_safe=False unchanged)
- No blockers for release preparation

---
*Phase: 27-association-performance-optimizations*
*Completed: 2026-02-22*
