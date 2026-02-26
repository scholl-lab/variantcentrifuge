---
phase: 38-codebase-cleanup
plan: 03
subsystem: pipeline
tags: [stages, exports, pipeline-core, cleanup]

# Dependency graph
requires: []
provides:
  - Complete __all__ exports for stages/__init__.py (6 previously missing stage classes)
  - Resolved TODO intelligent batching comment in pipeline_core/runner.py
affects: [any code doing "from variantcentrifuge.stages import *" or introspecting __all__]

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - variantcentrifuge/stages/__init__.py
    - variantcentrifuge/pipeline_core/runner.py

key-decisions:
  - "Sorted __all__ alphabetically (no section comments) to satisfy RUF022 isort-style lint rule"
  - "TODO intelligent batching replaced with deferred-work NOTE explaining why it is not yet needed"

patterns-established: []

# Metrics
duration: 5min
completed: 2026-02-26
---

# Phase 38 Plan 03: Missing Stage Exports and TODO Cleanup Summary

**Added 6 missing stage classes to stages/__init__.py __all__ and replaced TODO intelligent batching in runner.py with a deferred-work note**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-26T17:15:55Z
- **Completed:** 2026-02-26T17:21:22Z
- **Tasks:** 1
- **Files modified:** 2

## Accomplishments

- Added imports and __all__ entries for AssociationAnalysisStage, ClinVarPM5Stage, VariantAnalysisStage, DataSortingStage, ParallelCompleteProcessingStage, PCAComputationStage
- Sorted __all__ alphabetically to comply with ruff RUF022 rule (section-comments removed to allow full sort)
- Replaced `TODO: Implement intelligent batching` with `NOTE: Intelligent batching deferred — current fixed-size batching is sufficient for observed workloads. Revisit if stage count exceeds ~100.`

## Task Commits

Each task was committed atomically:

1. **Task 1: Add missing __all__ exports and resolve TODO comments** - `24a80a5` (chore)

**Plan metadata:** (see final docs commit below)

## Files Created/Modified

- `variantcentrifuge/stages/__init__.py` - Added 6 missing stage imports and __all__ entries; alphabetical sort
- `variantcentrifuge/pipeline_core/runner.py` - TODO replaced with deferred-work NOTE

## Decisions Made

- Removed section comments (# Setup stages, # Processing stages, etc.) from `__all__` because ruff RUF022 requires isort-style alphabetical sort and comments break it. The imports at the top of the file already group by module so the information is not lost.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed pre-existing formatting violation in analysis_stages.py**
- **Found during:** Task 1 verification (make ci-check)
- **Issue:** `analysis_stages.py` line 2197 had a multi-line f-string that ruff format wanted to collapse to a single line
- **Fix:** Ran `ruff format variantcentrifuge/stages/analysis_stages.py`
- **Files modified:** variantcentrifuge/stages/analysis_stages.py (the formatting fix was applied in-place by ruff and absorbed into an earlier stash cycle; it was not part of my commit since git diff showed it was already applied by the prior committed plan)
- **Verification:** `ruff format --check` passes with 255 files already formatted
- **Committed in:** Pre-existing in working tree; not a new introduction

---

**Total deviations:** 1 auto-fixed (pre-existing formatting issue discovered during ci-check)
**Impact on plan:** Necessary for CI to pass. No scope creep.

## Issues Encountered

- ruff RUF022 rejects `__all__` lists with section comments when they break alphabetical order. The fix was to remove section comments and sort purely alphabetically. The `ruff check --unsafe-fixes --fix` command was used to apply the canonical sort order directly.
- Several pre-existing test failures exist in the working tree (from prior plans 38-01/38-02): `test_inheritance` tests importing `create_inheritance_details` and `adjust_pattern_score`, and `test_chunked_loading` unit test. These are not caused by this plan's changes and are tracked as part of the broader 38 cleanup effort.

## Next Phase Readiness

- stages/__init__.py is now complete and discoverable
- runner.py TODO backlog is clear for this item
- Pre-existing test failures from 38-01/38-02 work need resolution before CI passes fully

---
*Phase: 38-codebase-cleanup*
*Completed: 2026-02-26*
