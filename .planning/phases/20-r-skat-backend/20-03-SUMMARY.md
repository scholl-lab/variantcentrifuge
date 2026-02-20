---
phase: 20-r-skat-backend
plan: 03
subsystem: testing
tags: [pytest, rpy2, mocking, skat, unit-tests, association]

# Dependency graph
requires:
  - phase: 20-r-skat-backend
    provides: "RSKATBackend (r_backend.py), RSKATTest (skat_r.py), engine lifecycle hooks"

provides:
  - "72 unit tests covering RSKATBackend, RSKATTest, and engine SKAT integration"
  - "sys.modules mocking pattern for testing deferred rpy2 imports"
  - "Verified: binary/continuous dispatch, NA_Real handling, GC interval, thread safety"
  - "Verified: engine None-effect guard, extra column naming, prepare/finalize lifecycle"

affects:
  - "Phase 21 (Pure Python SKAT): can use same test patterns for Davies/saddlepoint backend"
  - "Phase 22 (ACAT-O): engine column structure tests serve as regression guard"

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "sys.modules injection for deferred rpy2 imports: mock parent module attributes (mock_rpy2.robjects = mock_ro) for correct submodule resolution"
    - "NA_Real sentinel pattern: use unique object() as sentinel, inject via rpy2.rinterface.NA_Real in sys.modules"
    - "Backend mock pattern: inject _backend directly on RSKATTest to bypass check_dependencies"

key-files:
  created:
    - tests/unit/test_skat_r_backend.py
    - tests/unit/test_skat_r_test.py
    - tests/unit/test_engine_skat_integration.py
  modified:
    - variantcentrifuge/association/tests/skat_r.py
    - variantcentrifuge/association/backends/r_backend.py

key-decisions:
  - "sys.modules patch requires linking parent mock attributes (mock_rpy2.robjects = mock_ro); bare sys.modules injection does not work for nested submodule imports"
  - "Unused unpacked variables prefixed with _ per RUF059 (ruff); tests using mock_backend kept explicitly named"
  - "Engine lifecycle hook tests use MagicMock.wraps to preserve real Fisher behavior while spying on prepare/finalize calls"

patterns-established:
  - "rpy2 mock pattern: _patch_rpy2() returns (mods, mock_rpy2, mock_ro, mock_rpacks, na_sentinel) with parent attributes linked"
  - "Outer R result mock: _make_r_outer_result() returns outer.rx2('result')/outer.rx2('warnings') structure matching real SKAT output"
  - "Mock SKAT test factory: _make_mock_skat_test() creates AssociationTest with all-None effect columns and configurable extra dict"

# Metrics
duration: 16min
completed: 2026-02-20
---

# Phase 20 Plan 03: R SKAT Unit Tests Summary

**72 mocked unit tests verifying RSKATBackend (rpy2 dispatch, NA_Real, GC, thread safety), RSKATTest (null model caching, lifecycle hooks), and engine SKAT column output — all without requiring R installed**

## Performance

- **Duration:** 16 min
- **Started:** 2026-02-20T15:36:39Z
- **Completed:** 2026-02-20T15:53:05Z
- **Tasks:** 2/2
- **Files modified:** 5 (3 created, 2 reformatted)

## Accomplishments

- Created 30 RSKATBackend tests: detect_environment (3 scenarios), log_environment, fit_null_model (binary/continuous/thread guard), test_gene (SKATBinary vs SKAT dispatch, SKAT-O rho, NA_Real, warnings, column-major matrix, R exception propagation, globals cleanup), memory management (GC at 100-gene intervals, heap check), thread safety (4 tests)
- Created 27 RSKATTest tests: basic properties (name, effect_column_names all None), null model caching (fitted once across 3 genes), skip conditions (no genotype matrix, zero variants, no phenotype), extra columns (skat_o_rho, skat_method, skat_n_marker_test, skat_warnings), lifecycle hooks (prepare/finalize), GC triggering, progress tracking
- Created 15 engine integration tests: SKAT column structure (no skat_None, no double-prefix), extra columns written, Fisher+SKAT multi-test, correction behavior (None excluded from FDR), engine lifecycle hooks (prepare/finalize called with correct args)

## Task Commits

Each task was committed atomically:

1. **Task 1: RSKATBackend unit tests with mocked rpy2** - `3b796fc` (test)
2. **Task 2: RSKATTest + engine integration tests** - `e7eacaa` (test)

**Plan metadata:** (this commit)

## Files Created/Modified

- `tests/unit/test_skat_r_backend.py` — 30 tests for RSKATBackend with sys.modules mocking
- `tests/unit/test_skat_r_test.py` — 27 tests for RSKATTest with mocked backend
- `tests/unit/test_engine_skat_integration.py` — 15 engine integration tests with mock SKAT test
- `variantcentrifuge/association/tests/skat_r.py` — ruff format (no logic change)
- `variantcentrifuge/association/backends/r_backend.py` — ruff format (no logic change)

## Decisions Made

- **sys.modules mock hierarchy**: Discovered that `import rpy2.robjects.packages as rpacks` inside method bodies resolves via the parent `rpy2` module object, not directly from `sys.modules['rpy2.robjects.packages']`. Required `mock_rpy2.robjects = mock_ro` and `mock_ro.packages = mock_rpacks` to make submodule imports work correctly. This pattern is documented in test helpers.
- **NA_Real sentinel testing**: Each test creates a unique `object()` as the NA_Real sentinel and passes it to both the outer result builder and sys.modules patch, ensuring the `is NA_Real` identity check in `r_backend.py` works correctly.
- **Engine mock SKAT test**: Used `_make_mock_skat_test()` factory with configurable `extra_default` and `skip_genes` parameters to cover all engine column scenarios without any rpy2 dependency.

## Deviations from Plan

None — plan executed exactly as written. The mocking strategy evolved during implementation (Rule 1 style — working toward correct solution), but this is normal for mock-heavy test writing, not a deviation from scope.

## Issues Encountered

One non-trivial implementation challenge: the deferred `import rpy2.robjects as ro` pattern in `r_backend.py` requires that the parent `rpy2` mock object's `.robjects` attribute point to the `mock_ro` object. Simple `sys.modules` injection without parent-attribute linking causes `import rpy2.robjects.packages as rpacks` to resolve to `mock_rpy2.robjects.packages` (MagicMock auto-attribute) rather than our explicit `mock_rpacks`. Fixed by building a properly-linked hierarchy in `_patch_rpy2()`.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- Phase 20 complete: RSKATBackend skeleton (20-01), full implementation (20-02), unit tests (20-03)
- Phase 21 (Pure Python SKAT) can begin: same `_patch_rpy2()` helper can be adapted for Davies/saddlepoint backend tests
- Engine column structure tests in `test_engine_skat_integration.py` serve as regression guard for Phase 22 ACAT-O output

---
*Phase: 20-r-skat-backend*
*Completed: 2026-02-20*
