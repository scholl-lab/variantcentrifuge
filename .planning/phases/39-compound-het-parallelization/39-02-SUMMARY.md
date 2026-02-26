---
phase: 39-compound-het-parallelization
plan: 02
subsystem: performance
tags: [compound-het, parallelization, numpy, gil, inheritance, pedigree, threading]

# Dependency graph
requires:
  - phase: 39-01-compound-het-parallelization
    provides: baseline benchmark confirming GIL contention (parallel 0.21x-0.92x vs sequential)
  - phase: 38-codebase-cleanup
    provides: clean codebase starting point for optimization work
provides:
  - numpy-only parallel worker (_process_gene_group_arrays) with no DataFrame in hot path
  - pre-computed pedigree integer arrays (father_idx_arr, mother_idx_arr, affected_arr)
  - VARIANTCENTRIFUGE_BATCH_SIZE env var for batch size control
  - pre-dispatch deduplication in main thread before worker submission
  - post-optimization benchmark showing 1.06x-1.66x speedup over sequential
affects:
  - future parallelization work (pattern established for numpy-only worker dispatch)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Pre-dispatch deduplication: main thread deduplicates before submitting to workers
    - Pedigree arrays pattern: pre-compute int32/int8 arrays indexed by sample position for O(1) worker lookups
    - Numpy-only worker: workers receive arrays only, no DataFrame references in hot path
    - Environment variable override: VARIANTCENTRIFUGE_BATCH_SIZE with validation and fallback

key-files:
  created: []
  modified:
    - variantcentrifuge/inheritance/parallel_analyzer.py
    - variantcentrifuge/inheritance/comp_het_vectorized.py
    - .planning/phases/39-compound-het-parallelization/39-AFTER.txt
    - .planning/ROADMAP.md

key-decisions:
  - "Pre-dispatch dedup in main thread eliminates per-worker drop_duplicates() DataFrame call"
  - "Pedigree arrays use sample_to_idx for parent column lookups (int32 -1 sentinel for absent)"
  - "eff_sample_to_idx defensive fallback handles None gt_matrix edge case without assert"
  - "Tasks 1+2 committed together since import in Task 1 requires usage from Task 2 (lint would fail)"

patterns-established:
  - "Numpy-only worker pattern: receive (gene, n_unique, row_indices, variant_keys, gt_matrix, arrays...) not DataFrames"
  - "Pedigree array builder: build_pedigree_arrays() called once in main thread, result shared read-only to workers"

# Metrics
duration: 15min
completed: 2026-02-26
---

# Phase 39 Plan 02: Compound Het Parallelization Optimization Summary

**GIL contention eliminated in compound het Pass 2 by pre-computing pedigree integer arrays and dispatching numpy-only workers; parallel path achieves 1.06x-1.66x speedup (vs baseline 0.21x-0.92x slower than sequential)**

## Performance

- **Duration:** 15 min
- **Started:** 2026-02-26T22:05:24Z
- **Completed:** 2026-02-26T22:20:46Z
- **Tasks:** 3
- **Files modified:** 4

## Accomplishments

- Added `build_pedigree_arrays()` to `parallel_analyzer.py`: pre-computes `father_idx_arr` (int32), `mother_idx_arr` (int32), `affected_arr` (int8) indexed by sample position, eliminating per-worker `dict.get()` pedigree lookups
- Added `_get_batch_size()` to `parallel_analyzer.py`: respects `VARIANTCENTRIFUGE_BATCH_SIZE` env var with integer validation and auto-fallback to `min(4*n_workers, 500)`
- Added `_process_gene_group_arrays()` to `comp_het_vectorized.py`: numpy-only worker receiving pre-deduplicated row indices and pre-computed pedigree arrays; no DataFrame references, no `drop_duplicates()`, no `dict.get()` calls in hot path
- Rewired parallel dispatch in `analyze_inheritance_parallel()`: deduplication now done in main thread (pre-dispatch), pedigree arrays built once before `ThreadPoolExecutor` block, workers dispatched with numpy arrays only
- Benchmark improvement: all parallel configurations now faster than sequential (1.06x-1.66x) vs baseline (0.21x-0.92x); workers=2-16 at threshold=25-200 consistently outperform sequential

## Benchmark Comparison

| Config | Baseline (s) | After (s) | Speedup change |
|--------|-------------|-----------|----------------|
| sequential | 2.337 | 1.762 | baseline (-25% absolute improvement) |
| workers=2, threshold=25 | 6.024 (0.39x) | 1.246 (1.41x) | 0.39x → 1.41x |
| workers=4, threshold=25 | 6.681 (0.35x) | 1.096 (1.61x) | 0.35x → 1.61x |
| workers=16, threshold=25 | 5.681 (0.41x) | 1.069 (1.65x) | 0.41x → 1.65x |
| workers=2, threshold=50 | 4.319 (0.54x) | 1.086 (1.62x) | 0.54x → 1.62x |
| workers=4, threshold=200 | 5.285 (0.44x) | 1.065 (1.65x) | 0.44x → 1.65x |

Best result: workers=4, threshold=50 at 1.66x speedup over sequential.

## Task Commits

Each task was committed atomically:

1. **Task 1+2: Add pedigree array builder, numpy-only worker, and wire parallel dispatch** - `a987215` (feat)
2. **Task 3: Capture post-optimization benchmark and update roadmap** - `017554f` (feat)

## Files Created/Modified

- `variantcentrifuge/inheritance/parallel_analyzer.py` - Added `build_pedigree_arrays()`, `_get_batch_size()`, imported `_process_gene_group_arrays`; rewired parallel dispatch to use numpy-only workers with pre-computed pedigree arrays and pre-dispatch dedup
- `variantcentrifuge/inheritance/comp_het_vectorized.py` - Added `_process_gene_group_arrays()` numpy-only worker (165 lines); no DataFrame ops in hot path
- `.planning/phases/39-compound-het-parallelization/39-AFTER.txt` - Post-optimization benchmark results
- `.planning/ROADMAP.md` - Marked Phase 39 plans complete, v0.17.0 milestone shipped

## Decisions Made

- Pre-dispatch deduplication in main thread eliminates per-worker `drop_duplicates()` DataFrame call — each worker receives only the de-duplicated row indices as a numpy array
- Pedigree arrays use `sample_to_idx` for parent column lookups; int32 `-1` sentinel signals absent parent (avoids conditional dict lookups in worker hot path)
- `eff_sample_to_idx` defensive fallback handles `None gt_matrix` edge case gracefully without raising `AssertionError`
- Tasks 1 and 2 committed together since the import added in Task 1 requires usage wiring from Task 2 (lint would fail on unused import between separate commits)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Added ruff formatting fix for _process_gene_group_arrays**
- **Found during:** Task 3 (make ci-check)
- **Issue:** One line in `_process_gene_group_arrays` exceeded 100 chars (ruff format complaint)
- **Fix:** Applied `ruff format` to reformat the chained index + tolist() call
- **Files modified:** `variantcentrifuge/inheritance/comp_het_vectorized.py`
- **Verification:** `make ci-check` passes
- **Committed in:** `017554f` (Task 3 commit)

---

**Total deviations:** 1 auto-fixed (formatting)
**Impact on plan:** Minor formatting only. No scope creep.

## Issues Encountered

- Tasks 1 and 2 could not be committed separately because the import of `_process_gene_group_arrays` in Task 1 triggers lint failure (F401 unused import) until Task 2 wires usage. Combined into a single commit.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 39 complete; v0.17.0 milestone (Tech Debt Cleanup & Compound Het Parallelization) shipped
- Sequential branch and `_process_gene_group()` preserved for backward compatibility
- All 128 inheritance tests pass, 2066 non-slow tests pass, `make ci-check` passes
- No blockers for future work

---
*Phase: 39-compound-het-parallelization*
*Completed: 2026-02-26*
