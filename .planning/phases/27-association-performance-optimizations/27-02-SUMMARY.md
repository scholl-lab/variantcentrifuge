---
phase: 27-association-performance-optimizations
plan: 02
subsystem: association
tags: [association, config, cli, parallelization, plumbing, dataclass]

# Dependency graph
requires:
  - phase: 24-pure-python-coast-backend
    provides: AssociationConfig.coast_backend field (pattern followed)
  - plan: 27-01
    provides: parallel_safe=True on Python-backend tests (consumed by 27-03)
provides:
  - AssociationConfig.association_workers field (int, default=1)
  - --association-workers CLI argument (type=int, default=1, -1=auto)
  - _build_assoc_config_from_context reads association_workers from context.config
  - Unit test suite confirming the config plumbing end-to-end
affects:
  - 27-03: ProcessPoolExecutor parallelization (reads association_workers from AssociationConfig)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Config plumbing pattern: CLI arg -> cfg dict -> _build_assoc_config_from_context -> AssociationConfig field

key-files:
  created:
    - tests/unit/test_association_workers_config.py
  modified:
    - variantcentrifuge/association/base.py
    - variantcentrifuge/cli.py
    - variantcentrifuge/stages/analysis_stages.py

decisions:
  - id: IMPL-70
    decision: association_workers uses nullable=False in _get() (same as coast_backend, skat_backend)
    rationale: CLI always writes the key with its default (1); non-nullable is correct for typed int fields

metrics:
  duration: "~5 minutes"
  completed: "2026-02-22"
---

# Phase 27 Plan 02: Association Workers Config Plumbing Summary

**One-liner:** Pure config plumbing for --association-workers CLI arg flowing into AssociationConfig.association_workers int field (default=1, -1=auto).

## What Was Built

Added the configuration interface for gene-level parallelization that Plan 27-03's
ProcessPoolExecutor logic will read from. This is intentionally pure plumbing with
no parallelization logic — separating interface from implementation keeps each plan focused.

### Changes Made

**variantcentrifuge/association/base.py**
- Added `association_workers: int = 1` field to `AssociationConfig` dataclass after `coast_backend`
- Docstring explains: default=1 (sequential), >1 uses ProcessPoolExecutor, -1=os.cpu_count()

**variantcentrifuge/cli.py**
- Added `--association-workers` argument to `stats_group` (type=int, default=1)
- Added `cfg["association_workers"] = getattr(args, "association_workers", 1)` after coast_weights block

**variantcentrifuge/stages/analysis_stages.py**
- Added `"association_workers"` to `VALID_ASSOCIATION_KEYS` frozenset
- Added `"association_workers"` to `int_keys` validation set in `_validate_association_config_dict`
- Added `association_workers=_get("association_workers", default=1, nullable=False)` to `_build_assoc_config_from_context` return

**tests/unit/test_association_workers_config.py** (new)
- `test_association_config_default_workers`: AssociationConfig() has association_workers == 1
- `test_association_config_custom_workers`: AssociationConfig(association_workers=4) stores 4
- `test_association_config_auto_workers`: -1 sentinel stored as-is (cpu_count resolution in engine)
- `test_build_assoc_config_reads_workers`: stage builder reads from context.config correctly

## Decisions Made

| ID | Decision | Rationale |
|----|----------|-----------|
| IMPL-70 | association_workers uses nullable=False in _get() | CLI always writes the key with its default (1); non-nullable is correct for typed int fields. Same pattern as coast_backend and skat_backend. |

## Verification Results

- `pytest tests/unit/test_association_workers_config.py -v`: 4/4 passed
- `pytest -m "unit and not slow"`: 1565 passed, no regressions
- `make lint`: All checks passed

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Line too long in analysis_stages.py int_keys**

- **Found during:** Task 1 lint check
- **Issue:** Inline set `{"pca_components", "min_cases", "min_case_carriers", "firth_max_iter", "association_workers"}` exceeded 100-char line limit (108 chars)
- **Fix:** Expanded to multi-line set literal matching ruff E501 requirement
- **Files modified:** variantcentrifuge/stages/analysis_stages.py

## Next Phase Readiness

Plan 27-03 can now proceed. The ProcessPoolExecutor logic reads `config.association_workers`
from the `AssociationConfig` dataclass that `_build_assoc_config_from_context` populates.
The -1 sentinel convention (means os.cpu_count()) is established and documented.

## Task Commits

| Task | Name | Commit | Files |
|------|------|--------|-------|
| 1 | Add association_workers to config, CLI, and stage builder | 5210f84 | base.py, cli.py, analysis_stages.py |
| 2 | Unit test for association_workers config flow | 75e5f6b | test_association_workers_config.py |
