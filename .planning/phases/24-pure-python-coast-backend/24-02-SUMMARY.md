---
phase: 24-pure-python-coast-backend
plan: 02
subsystem: association
tags: [coast, allelic-series, backend-selection, cli, config]

# Dependency graph
requires:
  - phase: 24-01
    provides: PurePythonCOASTTest class in allelic_series_python.py
  - phase: 23-04
    provides: JSON config system, VALID_ASSOCIATION_KEYS, _build_assoc_config_from_context
  - phase: 21-02
    provides: skat_backend swap pattern (IMPL-29) in engine.from_names()
provides:
  - coast_backend field on AssociationConfig (default "auto")
  - COAST backend-aware swap in AssociationEngine.from_names()
  - --coast-backend CLI argument (choices: auto/r/python)
  - coast_backend wired through cfg dict into AssociationConfig
  - coast_backend in VALID_ASSOCIATION_KEYS for JSON config support
affects:
  - 24-03 (validation tests will verify backend selection behavior)
  - future users: --coast-backend python enables R-free COAST analysis

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Backend-aware swap in from_names() mirrors SKAT swap (IMPL-29): runs BEFORE unknown-name check"
    - "coast_backend follows exact skat_backend pattern in all three files (base.py, engine.py, cli.py)"

key-files:
  created: []
  modified:
    - variantcentrifuge/association/base.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/cli.py
    - variantcentrifuge/stages/analysis_stages.py

key-decisions:
  - "IMPL-59: coast_backend swap follows SKAT swap pattern exactly (post-skat, pre-unknown-name-check)"
  - "IMPL-60: coast auto mode probes rpy2 AND AllelicSeries importr() — more specific than SKAT auto"
  - "IMPL-61: coast_backend wired through analysis_stages.py _build_assoc_config_from_context() using _get() with nullable=False"

patterns-established:
  - "Backend selection pattern: CLI arg -> cfg dict -> AssociationConfig field -> engine swap (same for SKAT and COAST)"
  - "VALID_ASSOCIATION_KEYS + str_keys + enum validation must all be updated together for new string config fields"

# Metrics
duration: 3min
completed: 2026-02-22
---

# Phase 24 Plan 02: COAST Backend Selection Summary

**`--coast-backend python|r|auto` CLI flag wiring COAST backend swap into engine registry, matching the SKAT backend pattern exactly across base.py, engine.py, cli.py, and analysis_stages.py**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-22T03:03:52Z
- **Completed:** 2026-02-22T03:07:00Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Added `coast_backend` field to `AssociationConfig` (default `"auto"`) mirroring `skat_backend`
- Added backend-aware COAST swap in `AssociationEngine.from_names()`: python routes to PurePythonCOASTTest, auto tries R/AllelicSeries first then falls back, r keeps COASTTest — swap runs BEFORE unknown-name validation (IMPL-29 ordering invariant)
- Added `--coast-backend` CLI argument to stats group with choices [auto, r, python]
- Wired `coast_backend` through cfg dict into `AssociationConfig` via `_build_assoc_config_from_context()`
- Added `coast_backend` to `VALID_ASSOCIATION_KEYS`, `str_keys`, and enum validation in `_validate_association_config_dict()`

## Task Commits

Each task was committed atomically:

1. **Task 1: AssociationConfig coast_backend field + engine backend swap** - `4f06bee` (feat)
2. **Task 2: CLI --coast-backend argument + config wiring** - `6d6d0f7` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `variantcentrifuge/association/base.py` - Added `coast_backend: str = "auto"` field to AssociationConfig
- `variantcentrifuge/association/engine.py` - Added COAST backend swap block after SKAT swap in from_names()
- `variantcentrifuge/cli.py` - Added `--coast-backend` argument; wires `cfg["coast_backend"]`
- `variantcentrifuge/stages/analysis_stages.py` - Added coast_backend to VALID_ASSOCIATION_KEYS, str_keys, enum validation, and _build_assoc_config_from_context()

## Decisions Made

- **IMPL-59:** COAST backend swap placement: immediately after SKAT swap block, before `available = sorted(registry.keys())`. Follows IMPL-29 (backend-aware swap BEFORE unknown-name check).
- **IMPL-60:** COAST auto mode is more specific than SKAT auto: probes both `rpy2.robjects` AND `importr("AllelicSeries")` since rpy2 presence alone doesn't guarantee AllelicSeries is installed. SKAT auto only probes `get_skat_backend("r")`.
- **IMPL-61:** `coast_backend` uses `nullable=False` in `_get()` (same as `skat_backend`). The CLI always writes the key with its default, so non-nullable is correct for typed string fields with explicit choices.

## Deviations from Plan

None — plan executed exactly as written. The four-file scope (base.py, engine.py, cli.py, analysis_stages.py) matched the plan specification. The analysis_stages.py changes (VALID_ASSOCIATION_KEYS + str_keys + enum validation + AssociationConfig construction) were implied by the plan's "ensure the AssociationConfig construction picks up coast_backend" instruction.

## Issues Encountered

None.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- COAST backend selection fully wired; ready for Phase 24 plan 03 (validation tests)
- Plan 03 should verify: from_names() with coast_backend="python" instantiates PurePythonCOASTTest; from_names() with coast_backend="r" keeps COASTTest; auto mode falls back correctly when R/AllelicSeries unavailable (mock test)
- No blockers.

---
*Phase: 24-pure-python-coast-backend*
*Completed: 2026-02-22*
