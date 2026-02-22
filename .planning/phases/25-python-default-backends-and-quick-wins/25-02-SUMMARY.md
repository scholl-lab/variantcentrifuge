---
phase: 25-python-default-backends-and-quick-wins
plan: 02
subsystem: association
tags: [acat-v, skat, cauchy-combination, score-test, omnibus, rare-variant]

# Dependency graph
requires:
  - phase: 22-acat-o-diagnostics
    provides: cauchy_combination() and compute_acat_o() in acat.py
  - phase: 21-pure-python-skat-backend
    provides: PurePythonSKATTest, PythonSKATBackend with null model extras (residuals, sigma2, mu_hat)
provides:
  - compute_acat_v() per-variant score test in acat.py (Liu & Xie 2020, STAARpipeline)
  - ACAT-V wired into PurePythonSKATTest.run() stored in extra["acat_v_p"]
  - ACAT-V included in ACAT-O omnibus combination in engine._compute_acat_o()
  - 12 unit tests covering all edge cases and end-to-end wiring
affects:
  - phase: 26-documentation
  - phase: 27-performance-optimizations

# Tech tracking
tech-stack:
  added: [scipy.stats.norm (for normal distribution SF in score test p-values)]
  patterns:
    - "ACAT-V: per-variant score test with Cauchy combination; O(M) per gene"
    - "Single-pass loop: score vector + variance vector computed vectorized, then per-variant loop for valid p/weight pairs"
    - "Binary: diag(mu*(1-mu)) variance; quantitative: sigma2*I variance"
    - "mu_hat: np.ndarray | None — guarded by 'if trait_type == binary and mu_hat is not None'"
    - "ACAT-V stored in TestResult.extra dict (key: acat_v_p) — written to output TSV automatically"
    - "ACAT-V block in engine inserted AFTER test_pvals collection loop, BEFORE compute_acat_o() call"

key-files:
  created:
    - tests/unit/test_acat_v.py
  modified:
    - variantcentrifuge/association/tests/acat.py
    - variantcentrifuge/association/tests/skat_python.py
    - variantcentrifuge/association/engine.py

key-decisions:
  - "IMPL-62: compute_acat_v uses single-pass loop (vectorized score/variance, then per-variant p collection)"
  - "IMPL-63: mu_hat declared np.ndarray | None; all uses guarded by binary+not-None check"
  - "IMPL-64: acat_v_p=None set in ALL early-return extra dicts for schema consistency"
  - "IMPL-65: ACAT-V block in engine after collection loop, before compute_acat_o() to ensure inclusion"

patterns-established:
  - "ACAT-V uses same Beta(MAF; 1, 25) weights as SKAT for consistency"
  - "Extra keys in TestResult propagate to output TSV automatically (no engine changes needed for column)"

# Metrics
duration: 30min
completed: 2026-02-22
---

# Phase 25 Plan 02: ACAT-V Implementation Summary

**ACAT-V per-variant Cauchy score test added to acat.py and wired into PurePythonSKATTest + ACAT-O omnibus for sparse-signal detection power**

## Performance

- **Duration:** ~30 min
- **Started:** 2026-02-22T18:21:17Z
- **Completed:** 2026-02-22T18:50:00Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Implemented compute_acat_v() in acat.py: vectorized score/variance computation, single-pass loop for per-variant z-statistics and p-values, Cauchy combination with Beta(MAF) weights
- Wired ACAT-V into PurePythonSKATTest.run() using same beta weights as SKAT; stored in extra["acat_v_p"]
- Added ACAT-V block in _compute_acat_o() (after test_pvals collection, before compute_acat_o call) so it feeds into omnibus
- All early-return code paths set acat_v_p=None for schema consistency in output TSV
- 12 comprehensive unit tests: binary/quantitative traits, mu_hat=None, edge cases, signal detection, weight effects, ACAT-O integration, end-to-end wiring

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement compute_acat_v() and wire into SKAT + ACAT-O** - `f721ea8` (feat)
2. **Task 2: ACAT-V unit tests** - `3b8b0f2` (test)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `variantcentrifuge/association/tests/acat.py` - Added compute_acat_v() after cauchy_combination(); added `from scipy.stats import norm`
- `variantcentrifuge/association/tests/skat_python.py` - ACAT-V computation after test_gene(); acat_v_p in extra dict; early-return extra dicts get acat_v_p=None
- `variantcentrifuge/association/engine.py` - ACAT-V block in _compute_acat_o() after test_pvals loop, before compute_acat_o()
- `tests/unit/test_acat_v.py` - 12 unit tests for compute_acat_v and end-to-end PurePythonSKATTest wiring

## Decisions Made

- **IMPL-62**: compute_acat_v uses single-pass loop — vectorized score_vec and var_vec computed via numpy, then per-variant loop collects valid (p, weight) pairs. This avoids intermediate array allocation and is clear about the filtering logic.
- **IMPL-63**: mu_hat typed as `np.ndarray | None`; all uses inside function guarded by `if trait_type == "binary" and mu_hat is not None:`. Quantitative path always uses sigma2*I variance regardless of whether mu_hat was passed.
- **IMPL-64**: acat_v_p=None set in early-return extra dicts (NO_GENOTYPE_MATRIX and NO_PHENOTYPE_OR_EMPTY_MATRIX) so the key is always present regardless of code path, preventing KeyError when engine reads res.extra.
- **IMPL-65**: ACAT-V block in engine inserted AFTER the test_pvals collection loop and BEFORE compute_acat_o() call — critical ordering to ensure ACAT-V is included in the omnibus combination.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed lint error: unused n_samples variable**
- **Found during:** Task 1 verification (make lint)
- **Issue:** `n_samples, n_variants = geno.shape` — ruff RUF059: unpacked variable `n_samples` never used
- **Fix:** Changed to `_n_samples, n_variants = geno.shape` (prefix underscore convention)
- **Files modified:** variantcentrifuge/association/tests/acat.py
- **Verification:** `make lint` passes with no errors
- **Committed in:** f721ea8 (Task 1 commit)

**2. [Rule 1 - Bug] Fixed test data singular matrix in end-to-end test**
- **Found during:** Task 2 (test_pure_python_skat_stores_acat_v_in_extra)
- **Issue:** Test used `covariates = np.ones((n_samples, 1))` — backend calls `sm.add_constant(covariates)` which adds intercept column, making design matrix with two identical intercept columns = singular LinAlgError
- **Fix:** Removed explicit covariate_matrix from test (intercept-only null model via None); used balanced phenotype array (n_cases=[1]*25 + [0]*25) and p=[0.6,0.3,0.1] genotype probs matching existing SKAT test conventions
- **Files modified:** tests/unit/test_acat_v.py
- **Verification:** All 12 tests pass; no LinAlgError
- **Committed in:** 3b8b0f2 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 lint bug, 1 test data bug)
**Impact on plan:** Both auto-fixes necessary for correctness. No scope creep.

## Issues Encountered

None beyond the two auto-fixed items above.

## Next Phase Readiness
- ACAT-V fully implemented and tested; acat_v_p column appears in association output TSV for any gene tested with PurePythonSKATTest
- Phase 26 (Documentation) can reference ACAT-V as part of the omnibus pipeline description
- Phase 27 (Performance) unaffected — ACAT-V is O(M) per gene, computationally trivial

---
*Phase: 25-python-default-backends-and-quick-wins*
*Completed: 2026-02-22*
