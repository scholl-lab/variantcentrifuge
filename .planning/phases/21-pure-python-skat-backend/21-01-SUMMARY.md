---
phase: 21-pure-python-skat-backend
plan: 01
subsystem: statistics
tags: [skat, pvalue, davies, cffi, scipy, setuptools, c-extension, statistical-genetics]

# Dependency graph
requires:
  - phase: 20-r-skat-backend
    provides: SKATBackend ABC, NullModelResult, get_skat_backend() factory — all unchanged

provides:
  - davies.py: three-tier p-value fallback chain (Davies->Kuonen saddlepoint->Liu)
  - qfc.cpp: standalone CompQuadForm C++ source (Davies 1980) without R headers
  - _davies_build.py: CFFI out-of-line API build script for qfc C++ extension
  - setup.py: OptionalBuildExt + cffi_modules for C extension compilation
  - pyproject.toml: build backend migrated from hatchling to setuptools

affects:
  - 21-02 (PythonSKATBackend): calls compute_pvalue() from davies.py
  - 21-03 (validation tests): tests both Davies and fallback paths

# Tech tracking
tech-stack:
  added:
    - cffi 2.0.0 (runtime dep, CFFI out-of-line API mode for C++ extension)
    - setuptools>=64 (build dep, replaces hatchling for C extension support)
  patterns:
    - Three-tier p-value fallback chain with GENESIS+GMMAT hybrid triggering
    - CFFI out-of-line API mode with extern "C" forward declaration to bridge C++/C linkage
    - OptionalBuildExt pattern for non-fatal C extension builds

key-files:
  created:
    - variantcentrifuge/association/backends/davies.py
    - variantcentrifuge/association/data/qfc.cpp
    - variantcentrifuge/association/data/__init__.py
    - variantcentrifuge/_davies_build.py
    - setup.py
  modified:
    - pyproject.toml

key-decisions:
  - "CFFI set_source header must use extern C brackets to match C-linked qfc() symbol from C++ compiled qfc.cpp"
  - "qfc.cpp R headers replaced with standard <cmath>/<cstdlib>/<csetjmp> — math is identical"
  - "davies_pvalue returns (None, -1) when extension unavailable; compute_pvalue skips tier 1 cleanly"
  - "Davies triggers fallback at p<=1e-5 (proactive saddlepoint per GMMAT) even when ifault=0"

patterns-established:
  - "P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly in SKAT backend"
  - "Davies C extension: lazy-loaded once, cached in module globals, respects VARIANTCENTRIFUGE_NO_C_EXT"
  - "Fallback logging: INFO level at each tier transition (Davies->saddlepoint, saddlepoint->Liu)"

# Metrics
duration: 11min
completed: 2026-02-21
---

# Phase 21 Plan 01: P-Value Layer Summary

**Three-tier Davies->Kuonen saddlepoint->Liu p-value chain with compiled qfc C++ extension and setuptools build migration**

## Performance

- **Duration:** 11 min
- **Started:** 2026-02-21T07:53:52Z
- **Completed:** 2026-02-21T08:05:10Z
- **Tasks:** 2/2
- **Files modified:** 6 (4 created, 2 modified)

## Accomplishments

- Implemented `compute_pvalue()` API with GENESIS+GMMAT hybrid fallback triggering (Davies->saddlepoint->Liu)
- Compiled real CompQuadForm qfc.cpp via CFFI with correct C/C++ linkage bridging
- Migrated build system from hatchling to setuptools with non-fatal OptionalBuildExt
- All 1120 existing unit tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: davies.py — Liu, Kuonen saddlepoint, and Davies with fallback chain** - `2922f69` (feat)
2. **Task 2: Build system migration — hatchling to setuptools + cffi_modules** - `36d8450` (feat)

## Files Created/Modified

- `variantcentrifuge/association/backends/davies.py` - Three-tier p-value fallback chain with compute_pvalue() public API
- `variantcentrifuge/association/data/qfc.cpp` - Standalone CompQuadForm C++ source (R headers replaced with standard headers)
- `variantcentrifuge/association/data/__init__.py` - Empty package marker
- `variantcentrifuge/_davies_build.py` - CFFI out-of-line API build script; extern "C" forward declaration fixes C++/C linkage
- `setup.py` - OptionalBuildExt (non-fatal C extension build) + cffi_modules hook
- `pyproject.toml` - Build backend hatchling -> setuptools>=64; cffi added to runtime deps

## Decisions Made

- **CFFI C++/C linkage fix:** The CFFI-generated wrapper is compiled as C++ (due to source_extension='.cpp'). Without `extern "C"` in the forward declaration inside set_source, the wrapper produced C++ name-mangled symbol references (`_Z3qfcPdS_...`) that couldn't find the C-linked `qfc` from qfc.cpp. Fix: wrapped the forward declaration in `#ifdef __cplusplus extern "C" {}` guards.

- **qfc.cpp standalone version:** Original CompQuadForm qfc.cpp includes `<R.h>` and `"Rmath.h"`. Inspection showed no R-specific functions are called — only standard C math. Replaced with `<cmath>`, `<cstdlib>`, `<csetjmp>`. Math is byte-for-byte identical.

- **Davies proactive fallback threshold:** At p<=1e-5, Davies is suspected to have false convergence (GMMAT pattern). Even when ifault=0, compute_pvalue routes to saddlepoint for p<=1e-5 to avoid potential numerical cancellation errors near the integration singularity.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] CFFI C extension missing extern "C" forward declaration**

- **Found during:** Task 2 (build system migration, install verification)
- **Issue:** `_qfc.abi3.so` compiled but failed to load with `undefined symbol: _Z3qfcPdS_...` (C++ mangled). CFFI generated C++ wrapper code that called qfc() without C linkage, while qfc.cpp exported with C linkage.
- **Fix:** Added `#ifdef __cplusplus extern "C" {}` guards around the void qfc() forward declaration in `_davies_build.py`'s set_source argument.
- **Files modified:** `variantcentrifuge/_davies_build.py`
- **Verification:** `python -c "import variantcentrifuge._qfc; print('loaded')"` succeeds; `davies_pvalue(50.0, array)` returns ifault=0
- **Committed in:** `36d8450` (Task 2 commit)

**2. [Rule 1 - Bug] Ruff N802 lint errors in _kuonen_pvalue nested functions**

- **Found during:** Task 2 lint verification
- **Issue:** Nested functions `K()`, `Kprime()`, `Kdprime()` violated N802 (uppercase function names). Also SIM108 for if-else vs ternary.
- **Fix:** Renamed to `cgf()`, `cgf_prime()`, `cgf_dprime()`; converted if-else to ternary.
- **Files modified:** `variantcentrifuge/association/backends/davies.py`
- **Verification:** `make lint` passes cleanly
- **Committed in:** `36d8450` (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 Rule 1 bugs)
**Impact on plan:** Both fixes necessary for functional correctness and code quality. No scope creep.

## Issues Encountered

- **qfc.cpp originally had R headers:** `<R.h>` and `"Rmath.h"` were included but unused (no R-specific API calls in the file). Replaced with standard headers — trivial fix that enabled standalone compilation.
- **Davies returns p=-1.0 for small q:** For q=15 with lambdas=[5,3,1], Davies returns ifault=1 and p=-1.0 (a known behavior when the distribution is out of the integration range). The fallback to saddlepoint (p=0.1726) correctly handles this case.

## Next Phase Readiness

- `compute_pvalue()` is fully functional with all three tiers
- Davies C extension compiles and loads on Linux (tested with Python 3.10)
- Build system migration complete; `pip install -e ".[dev]"` succeeds
- Plan 21-02 (PythonSKATBackend) can import `from variantcentrifuge.association.backends.davies import compute_pvalue` immediately
- No blockers for Phase 21-02

---
*Phase: 21-pure-python-skat-backend*
*Completed: 2026-02-21*
