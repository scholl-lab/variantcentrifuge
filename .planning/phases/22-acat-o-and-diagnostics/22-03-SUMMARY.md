---
phase: 22-acat-o-and-diagnostics
plan: "03"
subsystem: testing
tags: [acat, cauchy-combination, lambda-gc, qq-plot, diagnostics, association, pytest]

# Dependency graph
requires:
  - phase: 22-01
    provides: cauchy_combination(), compute_acat_o() in acat.py; ARCH-03 FDR on ACAT-O only
  - phase: 22-02
    provides: diagnostics.py with compute_lambda_gc(), compute_qq_data(), write_diagnostics()

provides:
  - Comprehensive unit tests for ACAT Cauchy formula (numerical accuracy, edge cases)
  - Lambda_GC calibration tests on uniform null and inflated distributions
  - QQ data shape/sort/Hazen quantile formula tests
  - write_diagnostics file output tests (all 3 files: lambda_gc.tsv, qq_data.tsv, summary.txt)
  - Stage-level diagnostics integration test (config.diagnostics_output -> files on disk)
  - ACAT-O engine integration tests (column presence, k=1 pass-through, FDR, ARCH-03)

affects:
  - Phase 23 (can build on established test patterns)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Test classes named TestXxx group related assertions under pytest.mark.unit"
    - "np.random.default_rng(seed) used for reproducible random test data"
    - "Stage integration tests: call write_diagnostics() directly with mock results_df to verify wiring"
    - "Published-value tests pin implementation to specific reference numbers from papers"

key-files:
  created:
    - tests/unit/test_association_acat.py
    - tests/unit/test_association_diagnostics.py
  modified:
    - tests/unit/test_association_engine.py

key-decisions:
  - "cauchy_combination([0.01, 0.05]) result < 0.05 (not < 0.01) — Cauchy combination is not guaranteed to be smaller than the best individual p-value"
  - "Stage integration test calls write_diagnostics() directly (not via stage) — sufficient to validate config-to-disk wiring without full PipelineContext setup"

patterns-established:
  - "ACAT equal-input symmetry: cauchy_combination([p,p,p]) == p (proven numerically)"
  - "Lambda_GC bounds: [0.85, 1.15] for n=1000, [0.9, 1.1] for n=5000 null distributions"

# Metrics
duration: 10min
completed: 2026-02-21
---

# Phase 22 Plan 03: ACAT-O + Diagnostics Test Suite Summary

**91 unit tests validating ACAT Cauchy formula against Liu & Xie 2020 published values, lambda_GC null calibration, diagnostics file output, and ARCH-03 single-FDR engine integration**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-21T18:07:49Z
- **Completed:** 2026-02-21T18:18:48Z
- **Tasks:** 2
- **Files modified:** 3 (2 created, 1 updated)

## Accomplishments

- 27 ACAT tests: equal-input symmetry for p in {0.01, 0.05, 0.1, 0.5}, published reference value `0.00268` from Liu & Xie 2020 Table 1, tiny-p branch validation (`1e-20` → finite non-NaN result), weighted combination, and all edge cases (None, empty, single value, all-1.0)
- 40 diagnostics tests: lambda_GC on 1000 uniform null p-values within [0.85, 1.15] and inflated chi2 near 1.5, QQ sort order ascending, Hazen quantile formula verification, all three output files content, stage integration test confirming `config.diagnostics_output → write_diagnostics() → files on disk`
- 7 new ACAT-O engine integration tests: column presence, k=1 pass-through equality, FDR applied to corrected p-values, ARCH-03 enforcement (no `fisher_corrected_p_value` column), zero-variant gene exclusion, range validation in [0,1]

## Task Commits

1. **Task 1: ACAT formula, diagnostics tests, and stage integration** - `10ae6fa` (test)
2. **Task 2: Engine integration tests update for ARCH-03** - `1850cb2` (test)

## Files Created/Modified

- `tests/unit/test_association_acat.py` - 27 tests for `cauchy_combination()` and `compute_acat_o()`: numerical accuracy, published value, edge cases, pass-through logic
- `tests/unit/test_association_diagnostics.py` - 40 tests for diagnostics module: lambda_GC calibration, QQ data, warnings, file output, stage integration
- `tests/unit/test_association_engine.py` - Added `TestAssociationEngineAcatO` class with 7 tests for ACAT-O column presence, FDR strategy, and ARCH-03 enforcement; existing 17 tests unchanged

## Decisions Made

- `cauchy_combination([0.01, 0.05])` yields `0.01668`, not `< 0.01` — the Cauchy combination with equal weights is bounded between the two inputs, not guaranteed to be smaller than the best. Corrected test assertion from `< 0.01` to `< 0.05` (less significant input).
- Stage integration tests call `write_diagnostics()` directly rather than instantiating `AssociationAnalysisStage` — sufficient to validate config-to-disk wiring without requiring a full `PipelineContext` + workspace setup.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Incorrect test assertion for combined p-value magnitude**
- **Found during:** Task 1 (ACAT formula tests)
- **Issue:** `test_cauchy_combined_smaller_than_inputs` asserted `result < 0.01` for `cauchy_combination([0.01, 0.05])`. Actual result is `0.01668` (between the two inputs, not smaller than the best).
- **Fix:** Corrected assertion to `result < 0.05` (less significant input) with improved documentation explaining Cauchy combination behavior.
- **Files modified:** `tests/unit/test_association_acat.py`
- **Verification:** Test passes with corrected assertion.
- **Committed in:** `10ae6fa` (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug in test assertion)
**Impact on plan:** Corrected a faulty assumption about Cauchy combination magnitude. The formula correctly combines p-values but the combined result for equal weights on [0.01, 0.05] is 0.01668 (between inputs), not below the smallest input. No scope creep.

## Issues Encountered

None — plan executed cleanly after correcting the test assertion about Cauchy combination magnitude.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- Phase 22 fully complete: ACAT-O + diagnostics implementation (Plans 01-02) and validation (Plan 03) done.
- 1249 unit tests passing with zero regressions.
- CI passes: lint, format-check, typecheck, and test-fast all green.
- Phase 23 (PCA + Functional Weights + Allelic Series + JSON Config) can begin independently.

---
*Phase: 22-acat-o-and-diagnostics*
*Completed: 2026-02-21*
