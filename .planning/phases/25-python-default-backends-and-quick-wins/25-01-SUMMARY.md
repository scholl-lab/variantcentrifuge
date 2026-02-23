---
phase: 25
plan: 01
subsystem: association-testing
tags: [python-backend, skat, coast, deprecation, saddlepoint, davies]

dependency-graph:
  requires:
    - "24-pure-python-coast-backend"
    - "21-pure-python-skat-backend"
  provides:
    - "Python as default backend for SKAT and COAST"
    - "Deprecation warnings on R backend wrappers"
    - "Saddlepoint-before-Liu fallback in Davies compute_pvalue"
  affects:
    - "25-02-acat-v"
    - "26-documentation"

tech-stack:
  added: []
  patterns:
    - "Unified in('python','auto') pattern for backend swap"
    - "Three-tier p-value chain: Davies -> saddlepoint -> Liu"
    - "DeprecationWarning with stacklevel=2 on R test wrappers"

file-tracking:
  created:
    - tests/unit/test_backend_defaults.py
    - tests/unit/test_davies_fallback.py
  modified:
    - variantcentrifuge/association/backends/__init__.py
    - variantcentrifuge/association/base.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/association/tests/allelic_series.py
    - variantcentrifuge/association/tests/skat_r.py
    - variantcentrifuge/cli.py
    - variantcentrifuge/stages/analysis_stages.py
    - tests/unit/test_json_config.py

decisions:
  - id: IMPL-62
    decision: "auto backend selector maps directly to PythonSKATBackend (no R probe)"
    rationale: "R-probe in auto was expensive and error-prone; Python is always available"
  - id: IMPL-63
    decision: "Unified in('python','auto') pattern in engine.py for backend swap"
    rationale: "Single code path for auto and python avoids divergence if default changes"
  - id: IMPL-64
    decision: "DeprecationWarning only on RSKATTest/COASTTest __init__, not RSKATBackend"
    rationale: "RSKATTest.check_dependencies() constructs RSKATBackend; double-warning avoided"
  - id: IMPL-65
    decision: "Saddlepoint fallback inserted between Davies and Liu (not proactive threshold)"
    rationale: "Phase 25 scope: out-of-range only. Proactive threshold (IMPL-26) deferred to Phase 26+"

metrics:
  duration: "~20 minutes"
  completed: "2026-02-22"
  tasks-completed: 3
  tests-added: 24
  tests-total: 1968
---

# Phase 25 Plan 01: Python Default Backends and Quick Wins Summary

**One-liner:** Python backends as default for SKAT/COAST with R deprecation warnings and saddlepoint-before-Liu Davies fallback.

## What Was Built

Three changes that improve the out-of-box experience for users without R installed:

1. **Python-first backend defaults**: Both SKAT and COAST now default to the pure Python backend (validated in Phases 21+24). The `auto` selector and explicit `python` selector use identical code paths via the unified `in ("python", "auto")` pattern. Users with R can still use `--skat-backend r` / `--coast-backend r`.

2. **R backend deprecation warnings**: `RSKATTest()` and `COASTTest()` now emit `DeprecationWarning` on instantiation with `stacklevel=2` (warning points to caller code, not internal implementation). The warnings reference v0.17.0 as the removal target. R backends remain fully functional — deprecated, not removed.

3. **Saddlepoint-before-Liu Davies fallback**: When Davies returns an out-of-range p-value (p > 1.0 or p <= 0.0), the code now tries Kuonen saddlepoint approximation before falling through to Liu moment-matching. This matches SAIGE-GENE+ best practice for extreme-tail accuracy. The three-tier chain when Davies is available: Davies -> saddlepoint -> Liu.

## Tasks Completed

| Task | Name | Commit | Files |
|------|------|--------|-------|
| 1 | Swap default backends to Python + R deprecation warnings | 8eeffbe | backends/__init__.py, engine.py, base.py, cli.py, skat_r.py, allelic_series.py, analysis_stages.py, test_json_config.py |
| 2 | Saddlepoint-before-Liu fallback in compute_pvalue | 061664a | davies.py |
| 3 | New unit tests for defaults and fallback behavior | 479846f | test_backend_defaults.py, test_davies_fallback.py |

## Decisions Made

| Decision | Rationale |
|----------|-----------|
| `auto` maps directly to Python (no R probe) | R probe was expensive and error-prone; Python always available |
| Unified `in ("python", "auto")` pattern in engine.py | Single code path avoids divergence if default changes later |
| DeprecationWarning only on test wrappers, not RSKATBackend | RSKATTest.check_dependencies() constructs RSKATBackend — would double-warn |
| Out-of-range saddlepoint fallback only (no proactive threshold) | IMPL-26 proactive threshold deferred; Phase 25 scope = out-of-range only |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] `_build_assoc_config_from_context` had hardcoded "auto" defaults**

- **Found during:** Task 1 verification (test_json_config.py assertions)
- **Issue:** `analysis_stages.py` `_build_assoc_config_from_context()` had its own `default="auto"` for skat_backend/coast_backend independent of `AssociationConfig` defaults. Changing `AssociationConfig` defaults alone was insufficient.
- **Fix:** Updated `_get("skat_backend", default="python")` and `_get("coast_backend", default="python")` in `analysis_stages.py`
- **Files modified:** `variantcentrifuge/stages/analysis_stages.py`
- **Commit:** 8eeffbe

**2. [Rule 1 - Bug] Unused import in test_davies_fallback.py**

- **Found during:** Task 3 CI check
- **Issue:** `MagicMock` imported but not used (ruff F401)
- **Fix:** Removed `MagicMock` from import; `patch` alone was sufficient
- **Commit:** 479846f (applied via `make format`)

## Architecture Invariants Preserved

- R backends remain fully functional (not removed, just deprecated)
- `get_skat_backend("r")` still returns `RSKATBackend`
- `--skat-backend r` / `--coast-backend r` still work end-to-end
- Davies in-range behavior completely unchanged
- No proactive saddlepoint threshold implemented (Phase 26+ scope)

## Next Phase Readiness

Phase 25-02 (ACAT-V): Engine ACAT-V integration code was pre-seeded in engine.py during checker phase. Ready to implement `compute_acat_v()` in acat.py and add to `PurePythonSKATTest.run()`.
