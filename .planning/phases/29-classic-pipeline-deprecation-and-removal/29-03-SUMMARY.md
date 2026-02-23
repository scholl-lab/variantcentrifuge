---
phase: 29-classic-pipeline-deprecation-and-removal
plan: "03"
subsystem: docs
tags: [documentation, architecture, pipeline, deprecation]

# Dependency graph
requires:
  - phase: 29-01
    provides: pipeline.py refactored as stage-based entry point (classic code removed)
  - phase: 29-02
    provides: --use-new-pipeline CLI flag removed
provides:
  - CLAUDE.md with single stage-based pipeline architecture description
  - README.md without --use-new-pipeline feature flag reference
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Single-pipeline documentation: docs describe one architecture, not selectable modes"

key-files:
  created: []
  modified:
    - CLAUDE.md
    - README.md

key-decisions:
  - "Replace 'Two Pipeline Modes' section with 'Pipeline Architecture' describing pipeline.py as the stage-based entry point"
  - "Rename 'Stage-Based Pipeline (pipeline_core/)' subsection to 'Pipeline Core (pipeline_core/)' for clarity"

patterns-established:
  - "Documentation accuracy: CLAUDE.md reflects only what exists in code"

# Metrics
duration: 1min
completed: 2026-02-23
---

# Phase 29 Plan 03: Documentation Cleanup Summary

**CLAUDE.md and README.md updated to reflect single stage-based pipeline: removed "Two Pipeline Modes" section, --use-new-pipeline references, and classic pipeline descriptions**

## Performance

- **Duration:** 1 min
- **Started:** 2026-02-23T06:39:20Z
- **Completed:** 2026-02-23T06:40:27Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Replaced CLAUDE.md "Two Pipeline Modes" section with accurate "Pipeline Architecture" section describing pipeline.py as the stage-based entry point
- Removed --use-new-pipeline from README.md features list; stage-based pipeline now described as the architecture rather than an option
- Documentation now accurately reflects post-29-01/29-02 codebase where only one pipeline exists

## Task Commits

Each task was committed atomically:

1. **Task 1: Update CLAUDE.md architecture section** - `96c26b2` (docs)
2. **Task 2: Update README.md** - `6039a0f` (docs)

## Files Created/Modified
- `CLAUDE.md` - "Two Pipeline Modes" replaced with "Pipeline Architecture"; "Stage-Based Pipeline" subsection renamed to "Pipeline Core"
- `README.md` - Removed `(--use-new-pipeline)` parenthetical from feature list item

## Decisions Made
- Renamed the pipeline_core/ subsection heading from "Stage-Based Pipeline (pipeline_core/)" to "Pipeline Core (pipeline_core/)" to avoid implying there is a counterpart "non-stage-based pipeline"

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 29 is complete: pipeline.py is the stage-based pipeline, --use-new-pipeline flag removed, documentation updated
- All three plans in Phase 29 are done; the classic pipeline deprecation milestone is fully delivered
- Ready to update STATE.md and audit milestone v0.15.0 completion

---
*Phase: 29-classic-pipeline-deprecation-and-removal*
*Completed: 2026-02-23*
