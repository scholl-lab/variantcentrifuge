---
phase: 21-pure-python-skat-backend
plan: 02
subsystem: association
tags: [skat, numpy, scipy, statsmodels, glm, eigenvalues, burden-test, skat-o, python-backend]

# Dependency graph
requires:
  - phase: 21-01
    provides: "davies.py three-tier p-value layer (Davies C ext -> saddlepoint -> Liu)"
  - phase: 20-pure-python-skat-backend
    provides: "RSKATBackend and RSKATTest as pattern reference; SKATBackend ABC; NullModelResult"
provides:
  - PythonSKATBackend implementing SKATBackend ABC with statsmodels GLM null model
  - PurePythonSKATTest wrapping PythonSKATBackend following RSKATTest pattern
  - Shared _utils.py with parse_weights_beta (avoids rpy2 transitive import)
  - Updated get_skat_backend factory: 'python' and 'auto' fallback both implemented
  - Updated engine registry with 'skat_python' entry
  - Backend-aware swap in engine.from_names() for --skat-backend python/auto
affects:
  - "21-03: unit tests for Python SKAT backend"
  - "22: ACAT-O omnibus (uses association engine registry)"
  - "23: PCA and functional weights (may use Python backend)"

# Tech tracking
tech-stack:
  added:
    - statsmodels GLM (Binomial/Gaussian family for null model)
    - scipy.linalg.eigh with driver='evr' (matches R DSYEVR Lapack)
  patterns:
    - "Lazy null model caching: fit once before gene loop, reuse for all genes"
    - "Scientific variable renaming: uppercase math vars (ZW, Q, K) renamed to descriptive lowercase (geno_weighted, q_stat, kernel) for ruff N-series compliance"
    - "Shared utility module (_utils.py) to prevent transitive rpy2 import"
    - "Backend-aware registry swap in from_names() before unknown-name check"

key-files:
  created:
    - variantcentrifuge/association/backends/python_backend.py
    - variantcentrifuge/association/tests/skat_python.py
    - variantcentrifuge/association/tests/_utils.py
  modified:
    - variantcentrifuge/association/backends/__init__.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/association/tests/__init__.py
    - variantcentrifuge/association/tests/skat_r.py

key-decisions:
  - "Refactored _parse_weights_beta out of skat_r.py into shared tests/_utils.py to prevent rpy2 transitive import when Python backend uses weight parsing"
  - "SKAT-O uses minimum-p approach over rho grid (not full SKAT-O integration); sufficient for Phase 21, can refine in Phase 23"
  - "PythonSKATBackend._compute_eigenvalues_filtered() extracted as shared helper for SKAT and SKAT-O code paths"
  - "Backend-aware swap in engine.from_names() runs BEFORE unknown-name check so 'skat' resolves to PurePythonSKATTest when skat_backend='python'"

patterns-established:
  - "AssociationTest: check_dependencies() -> prepare() -> [run() per gene] -> finalize() lifecycle"
  - "GLM residuals: always resid_response (y - mu_hat), never deviance or Pearson"
  - "Eigenvalue threshold: mean(positive_lambdas) / 100_000 (exact R SKAT match)"
  - "Rank check on full genotype matrix BEFORE eigenvalue filtering"

# Metrics
duration: 10min
completed: 2026-02-21
---

# Phase 21 Plan 02: Pure Python SKAT Backend Summary

**PythonSKATBackend with statsmodels GLM null model + scipy eigh eigenvalues + SKAT-O rho grid, wired into engine as 'skat_python' with backend-aware 'skat' swap for --skat-backend python**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-21T08:10:36Z
- **Completed:** 2026-02-21T08:21:18Z
- **Tasks:** 2
- **Files modified:** 7 (4 modified, 3 created)

## Accomplishments

- PythonSKATBackend implements all 4 abstract methods (detect_environment, log_environment, fit_null_model, test_gene + cleanup) using only numpy/scipy/statsmodels
- SKAT, Burden, and SKAT-O methods all produce valid p-values verified against rank check and eigenvalue filtering matching R formulation
- PurePythonSKATTest follows exact RSKATTest pattern; includes skat_p_method and skat_p_converged in TestResult.extra
- Factory and engine registry fully wired: --skat-backend python routes "skat" to PurePythonSKATTest, "auto" falls back to Python when R unavailable
- 1122 unit tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: PythonSKATBackend — null model, score test, SKAT-O** - `a5bf994` (feat)
2. **Task 2: PurePythonSKATTest wrapper + factory/registry wiring** - `a504524` (feat)

## Files Created/Modified

- `variantcentrifuge/association/backends/python_backend.py` - PythonSKATBackend: GLM null model, SKAT/Burden/SKAT-O methods, Davies p-value chain
- `variantcentrifuge/association/tests/skat_python.py` - PurePythonSKATTest: AssociationTest wrapper with lazy null model, lifecycle hooks, extra metadata
- `variantcentrifuge/association/tests/_utils.py` - Shared parse_weights_beta() to avoid rpy2 transitive import
- `variantcentrifuge/association/backends/__init__.py` - Updated: 'python' returns PythonSKATBackend, 'auto' falls back to Python
- `variantcentrifuge/association/engine.py` - Updated: 'skat_python' in registry, backend-aware swap in from_names()
- `variantcentrifuge/association/tests/__init__.py` - Updated: PurePythonSKATTest in __getattr__ and __all__
- `variantcentrifuge/association/tests/skat_r.py` - Updated: imports _parse_weights_beta from shared _utils.py

## Decisions Made

- **[IMPL-25] ~~SKAT-O uses minimum-p approach~~ SUPERSEDED.** Full Lee et al. (2012) SKAT-O implemented post-verification (commit fc62e2c). Uses analytical R.M^{1/2} eigenvalue computation, SKAT_Optimal_Param decomposition, per-rho Davies→Liu p-values with Q-scale inversion, and omnibus chi-squared(1) integration. Matches R SKAT exactly on GCKD cohort (0.0000 log10 diff for PKD1/PKD2).

- **[IMPL-26] _parse_weights_beta moved to shared _utils.py**. If it stayed in skat_r.py, importing it from skat_python.py would transitively import rpy2 (via the RSKATBackend import at the top of skat_r.py). Moving to _utils.py prevents this and follows the principle of minimal coupling.

- **[IMPL-27] _compute_eigenvalues_filtered extracted as helper method**. The SKAT and SKAT-O code paths both need eigenvalue computation with the same threshold logic. Extracting to a method eliminates duplication and ensures consistent behavior.

- **[IMPL-28] Backend-aware swap in from_names() runs BEFORE unknown-name check**. This is critical: the check `if name not in registry` must happen after the swap so that "skat" resolves to the correct class before validation.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Extracted _compute_eigenvalues_filtered as shared method**

- **Found during:** Task 1 (_test_skato implementation)
- **Issue:** SKAT and SKAT-O both needed the same eigenvalue filtering logic; duplication would risk drift
- **Fix:** Extracted to `_compute_eigenvalues_filtered(kernel, sigma2) -> np.ndarray` on PythonSKATBackend
- **Files modified:** python_backend.py
- **Committed in:** a5bf994 (Task 1 commit)

**2. [Rule 1 - Bug] Ruff N-series naming violations fixed in python_backend.py**

- **Found during:** Task 2 (lint check)
- **Issue:** Scientific variable names (ZW, Zr, Q, K, X, T, V) violated ruff N806/N803 naming rules
- **Fix:** Renamed to descriptive lowercase: geno_weighted, score_vec, q_stat, kernel, design_matrix, score, variance
- **Files modified:** python_backend.py
- **Committed in:** a504524 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 missing critical extraction, 1 naming convention fix)
**Impact on plan:** Both auto-fixes improve code quality and compliance. No scope creep.

## Issues Encountered

None — plan executed smoothly. All verifications passed on first attempt.

## Next Phase Readiness

- PythonSKATBackend is feature-complete and ready for validation testing in Plan 21-03
- ~~The minimum-p SKAT-O approximation~~ SUPERSEDED: Full Lee et al. SKAT-O implemented post-verification (commit fc62e2c)
- All engine registry and factory wiring complete; Plan 21-03 can test --skat-backend python end-to-end

---
*Phase: 21-pure-python-skat-backend*
*Completed: 2026-02-21*
