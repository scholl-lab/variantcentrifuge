---
phase: 37-association-resource-management-memory-streaming
plan: 01
subsystem: pipeline
tags: [resource-manager, memory, pipeline-context, performance]

# Dependency graph
requires:
  - phase: 34-tech-debt
    provides: analysis_stages.py with ChunkedInheritanceAnalysisStage and InheritanceAnalysisStage
provides:
  - PipelineContext.resource_manager field (Optional, defaults to None)
  - One-time ResourceManager initialization in pipeline.py
  - Fallback pattern for all 4 local instantiations in analysis_stages.py
  - 8 unit tests for shared ResourceManager pattern
affects:
  - 37-02: genotype matrix streaming (Fix 5)
  - 37-03: memory reporting (Fix 6)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Shared singleton pattern: ResourceManager initialized once in pipeline.py, accessed via context.resource_manager"
    - "Fallback pattern: rm = context.resource_manager; if rm is None: rm = ResourceManager(config=context.config)"
    - "TYPE_CHECKING import for forward references in dataclass fields"

key-files:
  created:
    - tests/unit/test_pipeline_context_resource_manager.py
  modified:
    - variantcentrifuge/pipeline_core/context.py
    - variantcentrifuge/pipeline.py
    - variantcentrifuge/stages/analysis_stages.py

key-decisions:
  - "resource_manager absent from merge_from() — parent context's instance must not be overwritten by parallel contexts"
  - "Fallback pattern (if rm is None) ensures backward compatibility for tests that create PipelineContext without pipeline.py"
  - "cli.py and filters.py ResourceManager instances deliberately NOT touched — they run before context exists"
  - "parallel_analyzer.py deliberately NOT touched — receives DataFrame not PipelineContext"

patterns-established:
  - "Resource access: rm = context.resource_manager; if rm is None: from ..memory import ResourceManager; rm = ResourceManager(config=context.config)"
  - "TYPE_CHECKING block used for forward-reference imports in context.py to avoid circular imports"

# Metrics
duration: 17min
completed: 2026-02-25
---

# Phase 37 Plan 01: Shared ResourceManager on PipelineContext Summary

**Single ResourceManager initialized once in pipeline.py and shared across all stages via context.resource_manager, eliminating 4 redundant psutil.virtual_memory() calls in analysis_stages.py**

## Performance

- **Duration:** ~17 min
- **Started:** 2026-02-25T09:50:38Z
- **Completed:** 2026-02-25T10:07:10Z
- **Tasks:** 2
- **Files modified:** 3 modified, 1 created

## Accomplishments

- Added `resource_manager: "ResourceManager | None" = None` field to PipelineContext dataclass
- Initialized ResourceManager exactly once in pipeline.py (after context creation, before checkpoint init)
- Replaced all 4 local `ResourceManager(config=context.config)` instantiations in analysis_stages.py with context access + fallback
- Field absent from `merge_from()` — parallel context merges cannot overwrite parent's instance
- 8 unit tests cover: default None, accepts instance, set after creation, merge_from preservation (2 tests), fallback pattern, mock injection

## Task Commits

Each task was committed atomically:

1. **Task 1: Add resource_manager field to PipelineContext and initialize in pipeline.py** - `724be89` (feat)
2. **Task 2: Replace local ResourceManager instantiations with context access** - `ae8f5cc` (feat)

## Files Created/Modified

- `variantcentrifuge/pipeline_core/context.py` - Added TYPE_CHECKING import for ResourceManager, added resource_manager field
- `variantcentrifuge/pipeline.py` - One-time ResourceManager initialization after context creation
- `variantcentrifuge/stages/analysis_stages.py` - 4 local instantiations replaced with fallback pattern; 2 pre-existing SIM102 lint fixes
- `tests/unit/test_pipeline_context_resource_manager.py` - 8 unit tests for shared ResourceManager pattern

## Decisions Made

- `resource_manager` is intentionally absent from `merge_from()`. Parallel contexts share the parent's ResourceManager by design; allowing a child to overwrite it would break consistency and could lose the initialized instance.
- The fallback `if rm is None: rm = ResourceManager(config=context.config)` is kept in all 4 locations. This ensures stages work correctly in unit tests that create `PipelineContext` directly without going through `pipeline.py`.
- `cli.py` and `filters.py` ResourceManager instances were deliberately left unchanged — they run during CLI argument resolution before `PipelineContext` even exists.
- `parallel_analyzer.py` was deliberately left unchanged — it receives a DataFrame, not a `PipelineContext`, so it must create its own instance.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed pre-existing SIM102 lint errors in analysis_stages.py**

- **Found during:** Task 2 verification (make lint)
- **Issue:** Two nested `if` statements flagged by SIM102 at lines 1479 and 2478 in analysis_stages.py, pre-existing from working tree changes. Also one E501 line too long at line 1473.
- **Fix:** Flattened nested ifs into combined `if A and B and C` conditions; wrapped long list comprehension to fit 100-char limit
- **Files modified:** `variantcentrifuge/stages/analysis_stages.py`
- **Verification:** `make lint` passes
- **Committed in:** `ae8f5cc` (Task 2 commit)

**2. [Rule 1 - Bug] Applied ruff format to pre-existing formatting issues in 4 other files**

- **Found during:** Task 2 verification (make ci-check)
- **Issue:** 5 files had format violations (pre-existing in working tree): `test_association_engine.py`, `test_coast_python_comparison.py`, `test_engine_skat_integration.py`, `association/engine.py`, `stages/analysis_stages.py`
- **Fix:** Ran `python -m ruff format variantcentrifuge/ tests/` to apply all formatting fixes
- **Files modified:** The 5 files listed above (staged with analysis_stages.py in Task 2 commit)
- **Verification:** `make format-check` passes
- **Committed in:** `ae8f5cc` (Task 2 commit — formatting applied to analysis_stages.py; other 4 files staged)

---

**Total deviations:** 2 auto-fixed (2 Rule 1 - Bug)
**Impact on plan:** Both fixes were necessary for CI to pass. No scope creep.

## Issues Encountered

None — plan executed cleanly. The import ordering in the new test file and context.py TYPE_CHECKING block required a ruff fix (I001) but was caught immediately by lint check.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- PERF-04 (shared ResourceManager) complete
- context.resource_manager is the single source of truth for memory/CPU allocation
- Ready for 37-02: genotype matrix streaming (Fix 5) — context.resource_manager available for streaming decisions
- Ready for 37-03: memory reporting (Fix 6) — can access rm.memory_gb, rm.cpu_count from context

---
*Phase: 37-association-resource-management-memory-streaming*
*Completed: 2026-02-25*
