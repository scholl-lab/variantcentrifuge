---
phase: 23-pca-functional-weights-allelic-series-json-config
plan: 04
subsystem: association
tags: [json-config, association, matplotlib, qq-plot, validation, cli-override]

# Dependency graph
requires:
  - phase: 23-01
    provides: PCA integration, AssociationConfig PCA fields
  - phase: 23-02
    provides: Functional weights (CADD/REVEL), AssociationConfig weight fields
  - phase: 23-03
    provides: COAST allelic series test, AssociationConfig coast_weights field
  - phase: 22
    provides: ACAT-O, diagnostics framework, write_diagnostics(), QQ TSV output

provides:
  - JSON config mode: "association" section in config.json validated at startup
  - _validate_association_config_dict(): fail-fast, collects all errors before raising
  - _build_assoc_config_from_context(): CLI > JSON > default precedence builder
  - VALID_ASSOCIATION_KEYS: frozenset of all valid keys for "association" section
  - write_qq_plot(): lazy matplotlib import, Agg backend, headless HPC safe
  - write_diagnostics() now produces qq_plot.png alongside lambda_gc.tsv/qq_data.tsv
  - 35 new unit tests (24 json_config + 11 qq_plot); 1854 total passing

affects:
  - future association features: any new AssociationConfig fields need VALID_ASSOCIATION_KEYS update
  - HPC batch workflows: can now configure all association parameters via config.json

# Tech tracking
tech-stack:
  added: [matplotlib (optional, lazy import)]
  patterns:
    - JSON config validation with fail-fast + all-errors-collected approach
    - CLI > JSON > default precedence via _get() helper closure in builder function
    - Lazy matplotlib import with matplotlib.use("Agg") before pyplot for headless HPC
    - Refactor inline config construction into testable module-level builder function

key-files:
  created:
    - tests/unit/test_json_config.py
    - tests/unit/test_qq_plot.py
  modified:
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/association/diagnostics.py

key-decisions:
  - "VALID_ASSOCIATION_KEYS uses unprefixed field names matching AssociationConfig fields; CLI prefixed keys (association_min_cases) mapped via _get(cli_key, json_key) pairs"
  - "nullable=True/False in _get() closure: None in context.config means 'not set' for nullable fields, allowing JSON fallback; for non-nullable fields, key presence determines CLI override"
  - "matplotlib mock in absent-matplotlib tests uses sys.modules[matplotlib]=None trick (simpler than builtins.__import__ patch which causes recursion via pandas imports)"
  - "Integration test for matplotlib-absent in write_diagnostics() patches write_qq_plot directly (cleaner than sys.modules manipulation at that level)"
  - "association_tests (test names list) not in AssociationConfig; reads with JSON fallback separately from test_names resolution line in _process()"

patterns-established:
  - "Pattern: _validate_*_dict() collects all errors then raises with count in header"
  - "Pattern: _build_*_from_context() isolates config construction for testability"
  - "Pattern: lazy matplotlib import with use('Agg') in write_*() functions for optional HPC-safe visualization"

# Metrics
duration: 15min
completed: 2026-02-22
---

# Phase 23 Plan 04: JSON Config + QQ Plot Summary

**JSON config mode for association analysis via "association" section in config.json, with fail-fast validation and CLI override precedence; matplotlib QQ plot written alongside diagnostic TSVs**

## Performance

- **Duration:** 15 min
- **Started:** 2026-02-22T02:38:55Z
- **Completed:** 2026-02-22T02:54:08Z
- **Tasks:** 2
- **Files modified:** 4 (analysis_stages.py, diagnostics.py, test_json_config.py, test_qq_plot.py)

## Accomplishments

- JSON "association" section in config.json validated at startup with VALID_ASSOCIATION_KEYS frozenset, type checking, and enum validation — all errors collected before raising
- _build_assoc_config_from_context() replaces ~25-line inline AssociationConfig block in _process() with CLI > JSON > default precedence via _get() closure
- write_qq_plot() with lazy matplotlib import and Agg backend; returns bool; called automatically from write_diagnostics()
- 35 new unit tests (1854 total passing); CI fully clean

## Task Commits

Each task was committed atomically:

1. **Task 1: JSON config loading, validation, CLI override** - `701fdc7` (feat)
2. **Task 2: Matplotlib QQ plot in diagnostics.py and unit tests** - `2a33386` (feat)

**Plan metadata:** (see final commit below)

## Files Created/Modified

- `variantcentrifuge/stages/analysis_stages.py` - Added VALID_ASSOCIATION_KEYS, _validate_association_config_dict(), _build_assoc_config_from_context(); refactored AssociationAnalysisStage._process() AssociationConfig construction
- `variantcentrifuge/association/diagnostics.py` - Added write_qq_plot() with lazy matplotlib import; write_diagnostics() now calls write_qq_plot() after qq_data.tsv
- `tests/unit/test_json_config.py` - 24 unit tests: validation (15), config building (9)
- `tests/unit/test_qq_plot.py` - 11 unit tests: matplotlib-present, matplotlib-absent, empty-data, write_diagnostics integration

## Decisions Made

- **IMPL-51:** VALID_ASSOCIATION_KEYS uses unprefixed field names matching AssociationConfig fields; CLI prefixed keys (e.g. `association_min_cases`) are mapped via `_get(cli_key="association_min_cases", json_key="min_cases")` pairs — separating CLI and JSON key naming conventions cleanly
- **IMPL-52:** `nullable=True/False` parameter in `_get()` closure handles the distinction between "CLI wrote None (not set)" vs "CLI wrote non-None (set)"; for non-nullable fields, key presence in context.config determines CLI win, enabling JSON fallback for fields without a CLI flag
- **IMPL-53:** matplotlib mock for absent-matplotlib tests uses `sys.modules["matplotlib"] = None` approach; patching `builtins.__import__` globally causes RecursionError via pandas internal imports
- **IMPL-54:** `association_tests` (the test names list) is NOT in AssociationConfig; read with JSON fallback separately in the test_names resolution line in `_process()` using `context.config.get("association_tests") or _json_assoc.get("association_tests") or ["fisher"]`

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed matplotlib absent test using builtins.__import__ patch**

- **Found during:** Task 2 (test_qq_plot.py implementation)
- **Issue:** `patch("builtins.__import__", side_effect=...)` approach for simulating absent matplotlib caused `RecursionError: maximum recursion depth exceeded in __instancecheck__` because pandas internal imports also triggered the side effect during `write_diagnostics()` execution
- **Fix:** Used `sys.modules["matplotlib"] = None` to force ImportError on lazy import; for the write_diagnostics integration test, patched `write_qq_plot` directly with `return_value=False` to simulate absent matplotlib cleanly
- **Files modified:** tests/unit/test_qq_plot.py
- **Verification:** All 11 QQ plot tests pass
- **Committed in:** 2a33386

---

**Total deviations:** 1 auto-fixed (1 bug in test implementation)
**Impact on plan:** Test implementation corrected to use reliable import mocking strategy. No scope creep.

## Issues Encountered

- The `association_tests` field is in VALID_ASSOCIATION_KEYS (valid JSON config key) but is NOT a field on AssociationConfig — it feeds the test_names list for AssociationEngine. Required a separate JSON fallback resolution line in _process() beyond what _build_assoc_config_from_context() handles. Documented and implemented correctly.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

Phase 23 Plan 04 is the final plan in Phase 23. The full v0.15.0 Modular Rare Variant Association Framework is now complete:

- Phase 18: Foundation + Fisher refactor ✓
- Phase 19: Covariate system + burden tests ✓
- Phase 20: R SKAT backend ✓
- Phase 21: Pure Python SKAT backend ✓
- Phase 22: ACAT-O + diagnostics (lambda_GC, QQ TSV, sample size warnings) ✓
- Phase 23: PCA + Functional Weights + Allelic Series + JSON Config ✓

All 47 requirements across 6 phases are complete. The milestone is ready for v0.15.0 release tagging.

---
*Phase: 23-pca-functional-weights-allelic-series-json-config*
*Completed: 2026-02-22*
