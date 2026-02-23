---
phase: 28-tech-debt-cleanup
plan: 02
subsystem: cli
tags: [argparse, association, skat, diagnostics, config]

# Dependency graph
requires:
  - phase: 23-pca-functional-weights
    provides: _build_assoc_config_from_context with _get() calls for skat_method and diagnostic thresholds
  - phase: 27-performance-optimizations
    provides: association_workers CLI arg pattern (last written cfg assignment)
provides:
  - Four new CLI arguments exposing previously JSON-only AssociationConfig fields
  - --skat-method (choices: SKAT/Burden/SKATO, default: SKAT)
  - --min-cases (int, default: 200)
  - --max-case-control-ratio (float, default: 20.0)
  - --min-case-carriers (int, default: 10)
  - cfg dict wiring for all four args matching _get() cli_keys in analysis_stages.py
  - 12 unit tests validating propagation from cfg to AssociationConfig
affects: [29-classic-pipeline-deprecation, future CLI documentation, user guides]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "CLI arg -> cfg[key] -> _build_assoc_config_from_context() -> AssociationConfig"
    - "Always-write pattern: write cfg key with argparse default; _get() nullable=False handles precedence"
    - "Diagnostic threshold args use association_ prefix in cfg; skat_method uses plain key"

key-files:
  created:
    - tests/unit/test_cli_association_args.py
  modified:
    - variantcentrifuge/cli.py

key-decisions:
  - "Always-write cfg pattern for all four new args (matches skat_backend, coast_backend, association_workers)"
  - "skat_method uses plain cfg key (no association_ prefix); diagnostic thresholds use association_ prefix"
  - "Unused AssociationConfig import removed from test file after lint failure (auto-fixed)"

patterns-established:
  - "Phase 28 CLI args follow identical pattern to Phase 27 --association-workers"
  - "Test pattern: _make_context(config_dict) -> _build_assoc_config_from_context(ctx) -> assert fields"

# Metrics
duration: 4min
completed: 2026-02-23
---

# Phase 28 Plan 02: CLI Association Arguments Summary

**Four JSON-only association config fields exposed as CLI args (--skat-method, --min-cases, --max-case-control-ratio, --min-case-carriers) with 12 unit tests validating cfg propagation to AssociationConfig**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-23T06:34:38Z
- **Completed:** 2026-02-23T06:38:51Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Added four `add_argument()` calls in `stats_group` of `cli.py` (after `--association-workers` block)
- Wired four `cfg["key"] = getattr(args, ...)` assignments after `cfg["association_workers"]`
- Created `tests/unit/test_cli_association_args.py` with 12 unit tests covering all four new args
- CI check passes: 2001 tests passed, lint clean, format clean

## Task Commits

Each task was committed atomically:

1. **Task 1: Add CLI arguments and config wiring in cli.py** - `18769c6` (feat)
2. **Task 2: Write tests for CLI arg propagation** - `6b754f3` (test)

## Files Created/Modified

- `variantcentrifuge/cli.py` - Added four `add_argument()` calls in stats_group, four `cfg[key]` assignments in config assembly block
- `tests/unit/test_cli_association_args.py` - 12 unit tests covering CLI->cfg->AssociationConfig propagation

## Decisions Made

- **Always-write pattern for cfg assignments:** All four new args write their cfg key unconditionally (using argparse defaults), matching the existing pattern for `skat_backend`, `coast_backend`, and `association_workers`. The `_get()` helper in `_build_assoc_config_from_context` handles precedence correctly since non-nullable fields check `cli_key in cfg`.
- **Key naming:** `skat_method` uses the plain cfg key (no `association_` prefix) to match the `_get("skat_method", ...)` call in analysis_stages.py line 2278. The three diagnostic threshold args use the `association_` prefix to match their respective `_get("association_min_cases", json_key="min_cases", ...)` calls.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Removed unused AssociationConfig import from test file**

- **Found during:** Task 2 verification (`make ci-check`)
- **Issue:** Plan specified `from variantcentrifuge.association.base import AssociationConfig` in imports, but tests assert field values directly without `isinstance()` checks; ruff F401 unused import error
- **Fix:** Removed the unused import; tests work correctly without it
- **Files modified:** `tests/unit/test_cli_association_args.py`
- **Verification:** `make ci-check` passed after fix (2001 tests passed, lint clean)
- **Committed in:** `6b754f3` (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (unused import lint error)
**Impact on plan:** Minor — import not needed; all 12 tests pass and cover the same assertions.

## Issues Encountered

None - plan executed cleanly with one minor lint fix.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 28 Plan 02 complete. Phase 28 is now fully complete (both 28-01 and 28-02 done).
- Phase 29 (Classic Pipeline Deprecation) can proceed — no blockers from this phase.
- All four new CLI args are fully wired and tested.

---
*Phase: 28-tech-debt-cleanup*
*Completed: 2026-02-23*
