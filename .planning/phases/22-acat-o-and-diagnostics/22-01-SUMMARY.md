---
phase: 22-acat-o-and-diagnostics
plan: "01"
subsystem: association
tags: [acat, cauchy, omnibus, fdr, scipy, rare-variant, p-value-combination]

# Dependency graph
requires:
  - phase: 21-pure-python-skat
    provides: PurePythonSKATTest backend and TestResult infrastructure used by engine
  - phase: 18-foundation
    provides: AssociationEngine, AssociationConfig, AssociationTest ABC, correction.py
provides:
  - cauchy_combination() implementing Liu & Xie (2020) formula with numerical stability guard
  - compute_acat_o() per-gene omnibus p-value combinator
  - AssociationConfig warning threshold fields (min_cases, max_case_control_ratio, min_case_carriers, diagnostics_output)
  - _compute_acat_o() post-loop meta-test in engine
  - Single FDR pass on ACAT-O p-values only (ARCH-03)
  - acat_o_p_value and acat_o_corrected_p_value columns in engine output
affects:
  - 22-02: Diagnostics plan uses AssociationConfig.diagnostics_output and acat_o columns
  - 23-pca-weights: ACAT-V (variant-level) will use cauchy_combination() from this phase

# Tech tracking
tech-stack:
  added:
    - scipy.stats.cauchy (already in deps; first use for SF computation)
  patterns:
    - Post-loop meta-test pattern: ACAT-O runs after per-gene loop, not in the gene loop
    - Lazy import pattern: compute_acat_o() imported inside _compute_acat_o() method to avoid circular imports
    - ARCH-03 single-FDR pattern: one correction pass on omnibus p-values, primary tests uncorrected

key-files:
  created:
    - variantcentrifuge/association/tests/acat.py
  modified:
    - variantcentrifuge/association/base.py
    - variantcentrifuge/association/tests/__init__.py
    - variantcentrifuge/association/engine.py
    - tests/unit/test_association_engine.py
    - tests/unit/test_association_fisher.py
    - tests/unit/test_engine_skat_integration.py

key-decisions:
  - "cauchy_combination() uses 1/(p*pi) approximation for p < 1e-16 to avoid tan() overflow (Liu & Xie 2020)"
  - "Single valid p-value returns as pass-through (k=1 case per CONTEXT.md decision)"
  - "ACAT-O is NOT in _TEST_REGISTRY — it is a post-loop meta-test, not a primary gene-loop test"
  - "FDR applied only to ACAT-O across genes (ARCH-03); primary test corrected_p_value is always None"
  - "acat_o variable named t_stat in code (T is conventional in literature but fails ruff N806)"

patterns-established:
  - "Post-loop meta-test: call after finalize(), create separate results dict, apply FDR, add ACAT-O columns to rows"
  - "ARCH-03 FDR: single apply_correction() call on acat_o p-values; no per-test correction loops"
  - "Test updates for ARCH-03: replace {test}_corrected_p_value assertions with acat_o_corrected_p_value"

# Metrics
duration: 9min
completed: 2026-02-21
---

# Phase 22 Plan 01: ACAT-O + ARCH-03 FDR Summary

**Liu & Xie (2020) Cauchy combination formula implemented with overflow guard; ACAT-O post-loop meta-test added to engine; single FDR pass on ACAT-O omnibus p-values replaces per-test FDR (ARCH-03)**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-21T17:41:42Z
- **Completed:** 2026-02-21T17:51:00Z
- **Tasks:** 2
- **Files modified:** 6 (1 created)

## Accomplishments
- Created `variantcentrifuge/association/tests/acat.py` with `cauchy_combination()` (numerical stability guard for p < 1e-16) and `compute_acat_o()` (gene-level combinator)
- Integrated ACAT-O as a post-loop meta-test in `AssociationEngine._compute_acat_o()`: runs after all primary tests finalize, combines per-test p-values via Cauchy formula, applies single FDR correction
- Implemented ARCH-03 breaking change: removed per-test FDR loop, FDR applied only to `acat_o_p_value` across all genes; primary tests have no `corrected_p_value` column
- Added Phase 22 warning threshold fields to `AssociationConfig`: `min_cases`, `max_case_control_ratio`, `min_case_carriers`, `diagnostics_output`
- Updated 3 test files (17 individual test assertions) for ARCH-03 schema change; all 1175 unit tests pass

## Task Commits

1. **Task 1: ACAT Cauchy combination module + AssociationConfig warning fields** - `ddddaa6` (feat)
2. **Task 2: Engine ACAT-O post-loop computation + FDR strategy change** - `e6b7dcb` (feat)

**Plan metadata:** (included in task 2 commit)

## Files Created/Modified
- `variantcentrifuge/association/tests/acat.py` - cauchy_combination(), compute_acat_o(), _TINY_P_THRESHOLD constant
- `variantcentrifuge/association/base.py` - Added min_cases, max_case_control_ratio, min_case_carriers, diagnostics_output fields to AssociationConfig
- `variantcentrifuge/association/tests/__init__.py` - Lazy import block for cauchy_combination and compute_acat_o
- `variantcentrifuge/association/engine.py` - _compute_acat_o() method, ARCH-03 FDR change, acat_o columns in output
- `tests/unit/test_association_engine.py` - Updated 2 tests: expect acat_o columns, no fisher_corrected_p_value
- `tests/unit/test_association_fisher.py` - Updated bit-identity test: ACAT-O FDR instead of fisher FDR
- `tests/unit/test_engine_skat_integration.py` - Updated 4 tests for ARCH-03 column schema

## Decisions Made
- ACAT-O is NOT registered in `_TEST_REGISTRY` and NOT callable via `from_names(['acat_o'])`. It is computed post-loop as a meta-test. Adding it to the registry would allow users to request it as a primary gene-loop test, which makes no sense (it has no genotype input, only combines results).
- `cauchy_combination()` uses the `1/(p*pi)` approximation (Liu & Xie 2020) for `p < 1e-16` to avoid `tan((0.5 - p) * pi)` overflowing to ±inf for extremely small p-values.
- The Cauchy statistic variable is named `t_stat` in code (convention is `T` in literature, but ruff N806 enforces lowercase variables in functions).
- `p >= 1.0` values are filtered before combining (they contribute `tan(-pi/2) = -inf` which is numerically problematic and non-informative).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Three ruff lint errors in acat.py fixed before commit**
- **Found during:** Task 1 (lint check)
- **Issue:** `from typing import Sequence` → must use `collections.abc`; if-else weight normalization → ternary; variable `T` → lowercase `t_stat`
- **Fix:** Applied all three fixes: collections.abc import, ternary form, renamed variable
- **Files modified:** variantcentrifuge/association/tests/acat.py
- **Verification:** `ruff check variantcentrifuge/` returns "All checks passed!"
- **Committed in:** ddddaa6 (Task 1 commit includes the clean version)

**2. [Rule 1 - Bug] ARCH-03 breaking change affected 3 test files (not just the 1 specified in plan)**
- **Found during:** Task 2 verification
- **Issue:** Plan specified updating `tests/unit/test_association_engine.py` only; but `test_association_fisher.py` and `test_engine_skat_integration.py` also had `{test}_corrected_p_value` assertions that now fail
- **Fix:** Updated all affected tests — 4 additional test methods across 2 files
- **Files modified:** tests/unit/test_association_fisher.py, tests/unit/test_engine_skat_integration.py
- **Verification:** All 1175 unit tests pass
- **Committed in:** e6b7dcb (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 lint/style, 1 test scope undercount)
**Impact on plan:** Both auto-fixes required for correct operation. No scope creep.

## Issues Encountered
- None beyond deviations documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- `cauchy_combination()` is ready for reuse in Phase 23 ACAT-V (variant-level Cauchy)
- `AssociationConfig.diagnostics_output` field is ready for Phase 22 plan 02 (diagnostics output)
- `AssociationConfig.min_cases`, `max_case_control_ratio`, `min_case_carriers` ready for sample size warning logic
- Engine output now has `acat_o_p_value` and `acat_o_corrected_p_value` columns — downstream stages (Excel/HTML report) will need these added to their column handling in Phase 22 plan 02

---
*Phase: 22-acat-o-and-diagnostics*
*Completed: 2026-02-21*
