---
phase: 31-coast-fix
plan: 01
subsystem: association
tags: [coast, allelic-series, genotype-matrix, partial-category, regression]

# Dependency graph
requires:
  - phase: 24-pure-python-coast-backend
    provides: PurePythonCOASTTest and PythonCOASTBackend used by both fixes
  - phase: 19-covariate-system-and-burden-tests
    provides: AssociationAnalysisStage._process() where GT fallback was added
provides:
  - COAST-01: GT matrix fallback in branch 2 of GT reconstruction (variants_df is None)
  - COAST-03: partial-category fallback in both R and Python COAST backends
  - coast_status and coast_missing_categories in TestResult.extra
  - finalize() INFO summary: complete/partial/skipped counts per backend
  - Unit tests for both fixes (18+11 tests)
affects: [31-02, 35-weighted-coast, 36-sparse-matrices]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Partial-category tolerance: test proceeds with N-of-M categories present, returns None only when all empty"
    - "coast_status field in TestResult.extra: 'complete'/'partial'/'skipped' for downstream diagnostics"
    - "GT matrix fallback chain: variants_df -> current df -> None (log at INFO)"

key-files:
  created:
    - tests/unit/test_coast_partial_category.py
    - tests/unit/test_genotype_matrix_fallback.py
  modified:
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/association/tests/allelic_series_python.py
    - variantcentrifuge/association/tests/allelic_series.py
    - tests/unit/test_coast_python.py

key-decisions:
  - "COAST-03: skip only when ALL categories empty (len(present)==0), not when any single category missing"
  - "coast_status added to successful TestResult.extra: 'complete' when all 3 present, 'partial' when 1-2 present"
  - "Both R and Python COAST backends receive identical partial-category fix for behavioral consistency"
  - "Existing test_missing_categories_returns_none renamed and inverted to test_partial_categories_proceeds"

patterns-established:
  - "GT fallback: always check variants_df first, then current df, guard with if fallback_df is not None else []"
  - "COAST counters: _n_complete/_n_partial/_n_skipped initialized in prepare(), logged in finalize()"

# Metrics
duration: 9min
completed: 2026-02-23
---

# Phase 31 Plan 01: COAST-01 and COAST-03 Bug Fixes Summary

**COAST now produces p-values for real-world genes with 1-2 variant categories, and recovers the genotype matrix when variants_df is None by falling back to the current DataFrame**

## Performance

- **Duration:** ~9 min
- **Started:** 2026-02-23T17:53:26Z
- **Completed:** 2026-02-23T18:03:02Z
- **Tasks:** 2
- **Files modified:** 6 (3 source + 3 test)

## Accomplishments

- COAST-01 fixed: genotype matrix is now available when GT is already reconstructed and `context.variants_df` is None — the `elif` branch in `AssociationAnalysisStage._process()` now checks `df` as a second fallback
- COAST-03 fixed: both Python and R COAST backends now run COAST with 1 or 2 categories present; only skip when ALL three are empty (changed from skip-if-any-missing to skip-if-all-missing)
- `coast_status` ("complete"/"partial"/"skipped") and `coast_missing_categories` added to `TestResult.extra` for downstream diagnostics
- `finalize()` logs INFO summary: "COAST: N genes tested (X complete, Y partial, Z skipped)"
- All 1985 tests pass (including 18 new tests and 1 updated COAST test)

## Task Commits

1. **Task 1: COAST-01 GT matrix fallback and COAST-03 partial-category logic** - `68247a3` (fix)
2. **Task 2: Unit tests for partial-category and GT matrix fallback** - `6686b7b` (test)

**Plan metadata:** (see docs commit below)

## Files Created/Modified

- `variantcentrifuge/stages/analysis_stages.py` - Added `df` fallback in `elif` branch when `variants_df` has no GT columns; updated regression comment to include "coast"
- `variantcentrifuge/association/tests/allelic_series_python.py` - Replaced `if n_bmv==0 or n_dmv==0 or n_ptv==0` guard with `len(present)==0`; added `_n_complete/_n_partial/_n_skipped` counters; updated `finalize()` summary; added `coast_status` and `coast_missing_categories` to returned `extra`
- `variantcentrifuge/association/tests/allelic_series.py` - Same partial-category fix applied to R-backend `COASTTest.run()`; same counters and `finalize()` update
- `tests/unit/test_coast_partial_category.py` - 7 tests: all categories present (complete), PTV-only (partial), BMV+DMV only (partial), BMV only (partial), all empty (skips), finalize counters, pre-category guard
- `tests/unit/test_genotype_matrix_fallback.py` - 11 tests: `_find_per_sample_gt_columns` recognition + 7 fallback branch logic tests
- `tests/unit/test_coast_python.py` - Renamed `test_missing_categories_returns_none` -> `test_partial_categories_proceeds`; inverted assertion to expect numeric p_value with `coast_status=="partial"`

## Decisions Made

- COAST-03 skip threshold changed from "any category missing" to "ALL categories missing" — matches the R reference implementation (insitro/AllelicSeries `drop_empty=TRUE` behavior)
- `coast_status` field added to both success and skip outcomes for diagnostic transparency
- Both R and Python backends receive identical fix so behavior is consistent regardless of backend selection
- `variants_df` is checked before `df` in the fallback chain (preferred, as it preserves the pre-reconstruction state)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Updated existing test that asserted the old skip-on-missing behavior**

- **Found during:** Task 2 (writing unit tests)
- **Issue:** `test_missing_categories_returns_none` in `test_coast_python.py` asserted `p_value is None` for a PTV-only gene — this was now wrong after COAST-03 fix
- **Fix:** Renamed test to `test_partial_categories_proceeds`, inverted assertion to expect numeric p_value and `coast_status=="partial"`
- **Files modified:** `tests/unit/test_coast_python.py`
- **Verification:** All COAST tests pass (127 passed, 0 failed)
- **Committed in:** `6686b7b` (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug: outdated test asserting old skip behavior)
**Impact on plan:** Required for test suite correctness; the fix was minimal and the new assertion is more meaningful.

## Issues Encountered

- Ruff format check failed on `allelic_series_python.py` (multi-line logger.info that formatter collapsed to single line) — fixed with `ruff format`
- Ruff lint errors in test files (unused imports, unsorted imports, line too long) — fixed with `ruff check --fix` + `ruff format`

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- COAST now produces valid p-values for real-world cohorts where most genes lack some variant categories
- `coast_status` and `coast_missing_categories` available in `TestResult.extra` for diagnostics
- Ready for Phase 31 Plan 02 (remaining COAST bug fixes: multi-transcript SnpEff effect strings, classification scoring config)
- Phase 35 (weighted COAST) dependency on Phase 31 now unblocked for COAST-01/COAST-03 portion

---
*Phase: 31-coast-fix*
*Completed: 2026-02-23*
