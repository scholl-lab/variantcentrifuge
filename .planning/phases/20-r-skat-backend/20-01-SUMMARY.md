---
phase: 20
plan: 01
subsystem: association-backends
tags: [skat, rpy2, r-backend, abc, factory, engine]

dependency-graph:
  requires:
    - "19-01: genotype_matrix.py, covariates.py (matrix keys reused by RSKATTest)"
    - "18-01: AssociationTest ABC, TestResult, AssociationEngine (extended here)"
  provides:
    - "backends/ subpackage with SKATBackend ABC and NullModelResult"
    - "RSKATBackend skeleton with deferred rpy2 imports and eager R detection"
    - "get_skat_backend() factory with r/python/auto selector"
    - "RSKATTest AssociationTest wrapper registered as 'skat' in engine"
    - "Engine None-effect column guard (no more skat_None column names)"
    - "Engine extra dict loop for SKAT-specific output columns"
  affects:
    - "20-02: fills RSKATBackend.fit_null_model() and test_gene() stubs"
    - "21: PurePythonSKATBackend will implement the same SKATBackend ABC"
    - "22: ACAT-O engine integration reuses the extra column loop"

tech-stack:
  added: []
  patterns:
    - "Deferred rpy2 imports: all rpy2 usage gated behind method bodies"
    - "Factory function with deferred class imports (get_skat_backend)"
    - "Null model lifecycle: fit-once-reuse-per-gene caching in RSKATTest"
    - "Thread safety guard: _assert_main_thread() on all rpy2 entry points"

key-files:
  created:
    - variantcentrifuge/association/backends/__init__.py
    - variantcentrifuge/association/backends/base.py
    - variantcentrifuge/association/backends/r_backend.py
    - variantcentrifuge/association/tests/skat_r.py
  modified:
    - variantcentrifuge/association/base.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/association/tests/__init__.py
    - tests/unit/test_association_engine.py

decisions:
  - id: IMPL-17
    decision: "RSKATTest.check_dependencies() hardcodes backend='r'; no auto-detection at test level"
    rationale: "RSKATTest IS the R SKAT test; auto-selection between R/Python backends is at the factory level, not the test class level. A separate PurePythonSKATTest will exist in Phase 21."
  - id: IMPL-18
    decision: "Extra columns written with bare key names (e.g. skat_o_rho, not skat_skat_o_rho)"
    rationale: "Keys in TestResult.extra are already namespaced by the test (skat_o_rho). Adding test_name prefix would double-namespace producing skat_skat_o_rho, which is unreadable. Bare keys match the research column naming recommendation."
  - id: IMPL-19
    decision: "effect_column_names() return type is dict[str, str | None] not dict[str, str]"
    rationale: "SKAT has no effect size; all four effect slots are None. The base class default returns dict[str, str] but SKAT overrides with None values. Type broadened to accommodate this without breaking existing tests."

metrics:
  duration: "~26 minutes"
  completed: "2026-02-20"
  tests-passing: 1006
  files-created: 4
  files-modified: 4
---

# Phase 20 Plan 01: SKAT Backend Abstraction Layer Summary

**One-liner:** SKAT backend ABC + RSKATBackend skeleton with deferred rpy2 and eager R detection; RSKATTest registered in engine with None-effect guard.

## What Was Built

### backends/ Subpackage

Three new files establishing the SKAT computation layer:

**`backends/base.py`** — `SKATBackend` ABC with five abstract methods (`detect_environment`, `log_environment`, `fit_null_model`, `test_gene`, `cleanup`) and `NullModelResult` dataclass holding the fitted R null model object.

**`backends/r_backend.py`** — `RSKATBackend(SKATBackend)` with:
- ALL rpy2 imports deferred to method bodies (never at module level)
- `detect_environment()`: probes rpy2 availability, checks SKAT is installed in R, loads packages, extracts version strings — produces actionable ImportError messages with R_HOME diagnostics
- `log_environment()`: INFO line + version warnings for R < 4.0, SKAT < 2.0, rpy2 < 3.5.0
- `_assert_main_thread()`: RuntimeError guard on all rpy2 entry points
- Stub `fit_null_model()` and `test_gene()` raising `NotImplementedError("Implemented in Plan 20-02")`
- `cleanup()`: calls `ro.r["gc"]()` + Python `gc.collect()`

**`backends/__init__.py`** — `get_skat_backend(backend_name)` factory with "r"/"python"/"auto" modes; all class imports deferred inside function body.

### AssociationConfig Extension (base.py)

Added two fields to `AssociationConfig`:
- `skat_backend: str = "auto"` — selects R, Python, or auto-detect
- `skat_method: str = "SKAT"` — "SKAT", "Burden", or "SKATO"

### RSKATTest (tests/skat_r.py)

`RSKATTest(AssociationTest)` registered as `"skat"` in the engine:
- `effect_column_names()` returns all-None (SKAT has no effect size)
- `check_dependencies()` calls `detect_environment()` + `log_environment()` eagerly
- `run()` handles no-genotype-matrix and empty-matrix early returns; lazy null model caching
- `_parse_weights_beta()` converts "beta:1,25" format to (1.0, 25.0) tuple

### Engine Fixes (engine.py)

Two fixes applied to `run_all()`:

**None-effect column guard:** Replaced unconditional `row[f"{test_name}_{col_names['effect']}"]` with four guarded `if col_names.get("x") is not None:` checks. Without this fix, SKAT's all-None `effect_column_names()` would produce `skat_None` as a DataFrame column name.

**Extra columns loop:** Added `for extra_key, extra_val in res.extra.items(): row[extra_key] = extra_val` after the effect block. Writes SKAT-specific output columns (`skat_o_rho`, `skat_warnings`, `skat_method`, `skat_n_marker_test`) without double-prefixing. No-op for Fisher and burden tests (empty extra dicts).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Pre-existing test used "skat" as unknown-test sentinel**

- **Found during:** Task 2 verification (`pytest tests/unit/`)
- **Issue:** `test_from_names_unknown_test_raises_value_error` and `test_from_names_unknown_test_error_includes_test_name` both called `AssociationEngine.from_names(["skat"], ...)` expecting `ValueError`. Now that "skat" is registered, the engine instantiates `RSKATTest` and calls `check_dependencies()`, which raises `ImportError` (rpy2 not installed in this environment). The tests were placeholders written when "skat" was the canonical unknown test name.
- **Fix:** Changed both tests to use `"nonexistent_test_xyz"` as the unknown test sentinel. This correctly exercises the ValueError path without depending on a specific real test name being absent from the registry.
- **Files modified:** `tests/unit/test_association_engine.py`
- **Commit:** 67fe194

## Verification Results

All plan verification checks pass:

1. `from variantcentrifuge.association.backends import get_skat_backend` — OK
2. `sorted(_build_registry().keys())` — `['fisher', 'linear_burden', 'logistic_burden', 'skat']`
3. `RSKATTest().effect_column_names()` — `{'effect': None, 'se': None, 'ci_lower': None, 'ci_upper': None}`
4. `make lint && make format` — clean
5. `pytest tests/unit/` — 1006 passed (excluding pre-existing flaky timing test in `test_inheritance_parallel_integration.py`)

## Next Phase Readiness

Plan 20-02 can proceed immediately. It will:
1. Fill `RSKATBackend.fit_null_model()` — call `SKAT.SKAT_Null_Model()` or `SKAT.SKAT_Null_Model_Robust()` depending on trait_type
2. Fill `RSKATBackend.test_gene()` — call `SKAT.SKAT()` or `SKAT.SKATBinary()` with genotype matrix passed as R matrix
3. Handle numpy <-> R array conversion (rpy2 FloatVector/IntVector)
4. Extract p_value, rho, n_marker_test from R result list

No blockers. All invariants from ARCH-01 (parallel_safe=False) are correctly anticipated.
