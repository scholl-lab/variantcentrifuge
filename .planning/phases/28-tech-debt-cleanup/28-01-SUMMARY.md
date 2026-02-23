---
phase: 28-tech-debt-cleanup
plan: 01
subsystem: association-tests
tags: [tech-debt, parallel-safe, skat-r, verification, documentation]

dependency-graph:
  requires:
    - 20-r-skat-backend  # RSKATTest was originally written here
    - 22-acat-o-and-diagnostics  # VERIFICATION.md being updated
    - 27-performance-optimizations  # ARCH-05 added parallel_safe to Python backends
  provides:
    - RSKATTest explicit parallel_safe class attribute
    - Phase 22 lambda_gc verification gap closed
  affects:
    - 29-classic-pipeline-deprecation  # cleaner codebase going forward

tech-stack:
  added: []
  patterns:
    - parallel_safe class-level typed attribute (consistent across all AssociationTest subclasses)

file-tracking:
  created: []
  modified:
    - variantcentrifuge/association/tests/skat_r.py
    - .planning/phases/22-acat-o-and-diagnostics/22-VERIFICATION.md

decisions:
  - id: TECH-01
    decision: "parallel_safe: bool = False added as class attribute (not instance attribute) in RSKATTest"
    rationale: "Matches pattern of all other AssociationTest subclasses; engine uses getattr(test, 'parallel_safe', False) so class-level is correct"

metrics:
  duration: "5 minutes"
  completed: "2026-02-23"
  tasks: 2/2
---

# Phase 28 Plan 01: RSKATTest parallel_safe + lambda_gc gap closure Summary

**One-liner:** Added `parallel_safe: bool = False` class attribute to RSKATTest (last missing subclass) and closed Phase 22 lambda_gc.tsv/.txt documentation gap.

## What Was Built

Two targeted tech-debt items completed:

1. **RSKATTest parallel_safe attribute** — RSKATTest was the only `AssociationTest` subclass without an explicit `parallel_safe` class-level typed attribute. All other subclasses (FisherExactTest, LogisticBurdenTest, LinearBurdenTest, PurePythonSKATTest, COASTTest, PurePythonCOASTTest) declare it explicitly. This inconsistency is now resolved. The engine already uses `getattr(test, 'parallel_safe', False)` so the attribute was effectively False before — this change makes it explicit and type-annotated.

2. **Phase 22 lambda_gc verification gap** — The Phase 22 VERIFICATION.md reported Gap 2 as a filename mismatch: criterion said `lambda_gc.txt` but implementation created `lambda_gc.tsv`. On re-inspection, ROADMAP criterion #3 (line 168) already says `lambda_gc.tsv`. The gap was a documentation artifact from the original verification. The gap is now marked resolved, criterion 3 updated from PARTIAL to VERIFIED, and the score updated from 2/5 to 3/5.

## Tasks Completed

| Task | Name | Commit | Files |
|------|------|--------|-------|
| 1 | Add parallel_safe class attribute to RSKATTest | 5a2416f | variantcentrifuge/association/tests/skat_r.py |
| 2 | Close lambda_gc verification gap | 8c40157 | .planning/phases/22-acat-o-and-diagnostics/22-VERIFICATION.md |

## Decisions Made

| ID | Decision | Rationale |
|----|----------|-----------|
| TECH-01 | `parallel_safe: bool = False` as class attribute (not instance) | Matches pattern of all other subclasses; engine uses getattr() so class-level is semantically correct |

## Verification Results

1. `grep -n "parallel_safe" variantcentrifuge/association/tests/skat_r.py` — shows attribute at line 82 (class body, before `__init__`)
2. `grep "lambda_gc.tsv" .planning/ROADMAP.md` — confirmed at line 168 (criterion already correct)
3. `make ci-check` — 1989 passed, 3 skipped, 47 warnings — ALL CI CHECKS PASSED

## Deviations from Plan

None — plan executed exactly as written.

## Next Phase Readiness

Plan 02 of Phase 28 (CLI args for skat_method and diagnostic thresholds) can proceed immediately. No blockers introduced.

The RSKATTest subclass consistency gap is closed. The Phase 22 verification record is accurate.
