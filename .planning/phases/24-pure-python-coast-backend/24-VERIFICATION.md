---
phase: 24-pure-python-coast-backend
verified: 2026-02-22T03:34:48Z
status: gaps_found
score: 4/5 must-haves verified
gaps:
  - truth: "Python COAST produces p-values within tiered tolerance of R AllelicSeries::COAST()"
    status: partial
    reason: >
      The tiered tolerance test infrastructure exists (tier1: rel < 1e-4, tier2: log10 < 0.5)
      but both tolerance tests compare Python output against itself (two identical runs),
      not against actual R AllelicSeries::COAST() reference values. The R golden script
      (scripts/generate_coast_golden.R) is implemented and would generate the reference values,
      but the output has NOT been run and the values have NOT been hardcoded into
      test_coast_python_comparison.py. The SUMMARY explicitly labels this as a
      "Placeholder pattern: comparison tests use statistical ranges (no hardcoded R values)
      until user runs R script offline." The implementation framework is correct; the
      R-vs-Python numeric comparison is absent.
    artifacts:
      - path: "tests/unit/test_coast_python_comparison.py"
        issue: >
          test_tier1_tolerance_relative_check and test_tier2_tolerance_log10_check both
          compute r1 and r2 from the same Python call (not from R). The comparison verifies
          Python self-consistency, not agreement with R AllelicSeries::COAST(). The comment
          at line 256 says 'R golden values (PLACEHOLDER - update after running
          scripts/generate_coast_golden.R)'.
      - path: "tests/fixtures/coast_golden/README.md"
        issue: >
          README says values are 'embedded as constants' in the comparison test file,
          but no such constants exist (no _R_GOLDEN_* dicts). The directory contains
          only the README.
    missing:
      - "Run scripts/generate_coast_golden.R with R + AllelicSeries installed"
      - "Copy _R_GOLDEN_* p-value dicts from R script output into test_coast_python_comparison.py"
      - "Replace self-consistency tolerance tests with actual Python-vs-R comparison assertions"
human_verification:
  - test: "Run R golden script and compare Python output"
    expected: >
      After running `Rscript scripts/generate_coast_golden.R`, Python COAST p-values for
      the 5 scenarios match within tiered tolerance: rel < 1e-4 for p > 1e-4;
      log10 diff < 0.5 (i.e. within half an order of magnitude) for p <= 1e-4.
    why_human: >
      Requires R + AllelicSeries package to be installed. Cannot be verified in CI
      on this machine (AllelicSeries not installed). The R script exists and is ready;
      the numeric comparison requires a human with the R environment to execute it and
      populate the test constants.
---

# Phase 24: Pure Python COAST Backend Verification Report

**Phase Goal:** Users without R can run the COAST allelic series test via a pure Python implementation that matches R AllelicSeries::COAST() output within tiered tolerance, reusing existing SKAT, burden regression, and Cauchy combination infrastructure — enabling `parallel_safe=True` and eliminating the R/rpy2 dependency for COAST.

**Verified:** 2026-02-22T03:34:48Z
**Status:** gaps_found (4/5 success criteria verified; SC1 partial)
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Python COAST p-values match R AllelicSeries::COAST() within tiered tolerance | PARTIAL | Tolerance test framework exists but compares Python-to-Python; no hardcoded R golden values |
| 2 | --coast-backend python works without R installed | VERIFIED | PurePythonCOASTTest uses only numpy/scipy/statsmodels; no rpy2 import in source; check_dependencies() succeeds with mocked rpy2 |
| 3 | PurePythonCOASTTest.parallel_safe is True | VERIFIED | `PurePythonCOASTTest.__dict__['parallel_safe'] = True`; COASTTest has False |
| 4 | All existing COAST tests pass; no regressions | VERIFIED | 1916 tests pass, 43 existing test_coast.py tests pass, 62 new COAST Python tests pass |
| 5 | --coast-backend python/r/auto backend selection mirrors SKAT pattern | VERIFIED | CLI, AssociationConfig, engine.from_names(), analysis_stages.py all wired; auto probes rpy2+AllelicSeries then falls back |

**Score:** 4/5 truths verified (1 partial)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/association/backends/coast_python.py` | PythonCOASTBackend with 6 burden + 1 SKAT + Cauchy | VERIFIED | 472 lines; real implementation; 3 layers (_aggregate_by_category, _run_burden_test, _compute_allelic_skat_weights, test_gene); no rpy2 imports |
| `variantcentrifuge/association/tests/allelic_series_python.py` | PurePythonCOASTTest wrapper with parallel_safe=True | VERIFIED | 480 lines; parallel_safe class attribute = True; name = "coast"; all 7 skip guards; extra keys match COASTTest |
| `variantcentrifuge/association/base.py` | coast_backend field on AssociationConfig | VERIFIED | Line 173: `coast_backend: str = "auto"` |
| `variantcentrifuge/association/engine.py` | Backend swap in from_names() | VERIFIED | Lines 150-173: python/auto swap before unknown-name check; auto probes importr("AllelicSeries") |
| `variantcentrifuge/cli.py` | --coast-backend CLI argument | VERIFIED | Lines 429-434: `--coast-backend` with choices=[auto, r, python]; line 1146 wires to cfg |
| `variantcentrifuge/stages/analysis_stages.py` | coast_backend in VALID_ASSOCIATION_KEYS, str_keys, validation, config build | VERIFIED | Lines 2053, 2110, 2157-2160, 2272: fully wired in all 4 required locations |
| `tests/unit/test_coast_python.py` | 41 unit tests for COAST backend layers | VERIFIED | 41 tests pass across 5 classes (aggregation, burden, SKAT weights, omnibus, lifecycle) |
| `tests/unit/test_coast_python_comparison.py` | R reference comparison tests | PARTIAL | 21 tests exist and pass; tolerance tests verify self-consistency, not Python-vs-R comparison |
| `scripts/generate_coast_golden.R` | R script for offline golden value generation | VERIFIED | 238 lines; 5 scenarios; structured for copy-paste into Python test constants |
| `tests/fixtures/coast_golden/README.md` | Documentation for golden value workflow | VERIFIED | Exists; correctly documents the approach |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| PurePythonCOASTTest.run() | PythonCOASTBackend.test_gene() | self._backend.test_gene() | WIRED | allelic_series_python.py line 398 |
| PythonCOASTBackend.test_gene() | cauchy_combination() | from acat import cauchy_combination | WIRED | coast_python.py line 392 |
| PythonCOASTBackend._run_allelic_skat() | PythonSKATBackend._compute_eigenvalues_filtered() | self._skat_backend._compute_eigenvalues_filtered() | WIRED | coast_python.py line 335 |
| PythonCOASTBackend._run_allelic_skat() | compute_pvalue() | from davies import compute_pvalue | WIRED | coast_python.py line 320 |
| AssociationEngine.from_names() | PurePythonCOASTTest | registry["coast"] = PurePythonCOASTTest | WIRED | engine.py lines 154-158 (python), 168-172 (auto) |
| CLI --coast-backend | AssociationConfig.coast_backend | cfg["coast_backend"] -> _build_assoc_config_from_context() | WIRED | cli.py line 1146 + analysis_stages.py line 2272 |
| test_coast_python_comparison.py | R golden values | hardcoded _R_GOLDEN_* constants | NOT WIRED | No R golden values exist; PLACEHOLDER comment at line 256 |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| COAST-PY-01: P-values match R within tiered tolerance | PARTIAL | Tolerance tests compare Python-vs-Python; no R reference values hardcoded |
| COAST-PY-02: R-free COAST pipeline | SATISFIED | No rpy2 in source; check_dependencies() works without R |
| COAST-PY-03: parallel_safe=True | SATISFIED | Class attribute verified True |
| COAST-PY-04: No regressions | SATISFIED | 1916 tests pass, 0 failures |
| COAST-PY-05: Backend selection --coast-backend python/r/auto | SATISFIED | All 4 files wired; auto probes AllelicSeries specifically |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `tests/unit/test_coast_python_comparison.py` | 256 | `R golden values (PLACEHOLDER - update after running...)` | Warning | The R-vs-Python comparison is not yet implemented; tolerance tests are self-consistency only |

### Human Verification Required

#### 1. R Golden Value Comparison (COAST-PY-01)

**Test:** With R and AllelicSeries installed, run:
```bash
Rscript scripts/generate_coast_golden.R
```
Then copy the `_R_GOLDEN_*` output dicts into `tests/unit/test_coast_python_comparison.py` as constants and replace the `test_tier1_tolerance_relative_check` and `test_tier2_tolerance_log10_check` tests with actual Python-vs-R assertions using the tiered tolerance formula.

**Expected:** For each of the 5 scenarios:
- p > 1e-4: relative difference < 1e-4
- p <= 1e-4: |log10(p_python) - log10(p_R)| < 0.5

**Why human:** Requires R with AllelicSeries installed. AllelicSeries is not installed on the CI machine (verified: `PackageNotInstalledError: The R package "AllelicSeries" is not installed.`). The R script is fully written and ready to run; the bottleneck is R environment availability.

### Gaps Summary

The implementation is 4/5 complete. The infrastructure for comparing Python COAST output to R AllelicSeries::COAST() exists — the R script is written (238 lines covering 5 scenarios), the Python test file has the tolerance scaffold, and the tolerance math (tier1: rel < 1e-4; tier2: log10 < 0.5) is correctly specified. However, the actual R reference values were never generated and hardcoded. The SUMMARY explicitly documents this as intentional: "Placeholder pattern: comparison tests use statistical ranges (no hardcoded R values) until user runs R script offline."

All other requirements are fully met:
- Pure Python backend with no R/rpy2 dependency: VERIFIED
- parallel_safe=True as class attribute: VERIFIED
- --coast-backend python/r/auto wired through all layers: VERIFIED
- Zero regressions (1916 tests, 0 failures): VERIFIED

The gap is COAST-PY-01 only: the tiered tolerance comparison against R has not been executed and the R golden values have not been hardcoded into the test suite.

---

_Verified: 2026-02-22T03:34:48Z_
_Verifier: Claude (gsd-verifier)_
