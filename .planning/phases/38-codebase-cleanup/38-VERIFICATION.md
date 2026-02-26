---
phase: 38-codebase-cleanup
verified: 2026-02-26T18:40:22Z
status: passed
score: 5/5 must-haves verified
---

# Phase 38: Codebase Cleanup Verification Report

**Phase Goal:** The codebase contains no dead code, stale documentation, or unresolved minor tech debt — any developer reading the source sees accurate, current state.
**Verified:** 2026-02-26T18:40:22Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `stage_info.py` is deleted and nothing imports it | VERIFIED | File absent; `from variantcentrifuge.stage_info` import search returns zero results. Variable-name uses of `stage_info` in `display_utils.py` and `interactive_resume.py` reference `stage_registry.StageInfo` — a different module. |
| 2 | `--coast-backend r` is not present in CLI choices; passing it raises an argument error | VERIFIED | `cli.py` line 439: `choices=["auto", "python"]`. Runtime check: `python -c "create_parser().parse_args(['--coast-backend', 'r', ...])"` produces `error: argument --coast-backend: invalid choice: 'r' (choose from 'auto', 'python')`. |
| 3 | Dead functions in `prioritizer.py` (6 functions) and `analyzer.py` (4 functions) are removed with no broken call sites | VERIFIED | `prioritizer.py` contains only 3 functions: `prioritize_patterns`, `calculate_confidence`, `get_pattern_description` (191 lines). `analyzer.py` has no `get_inheritance_summary`, `filter_by_inheritance_pattern`, `export_inheritance_report`, or public `create_inheritance_details`; only `_create_inheritance_details` (private) exists. All 48 inheritance tests pass. |
| 4 | `docs/faq.md`, `docs/association_testing.md`, and `docs/changelog.md` accurately reflect the current CLI — no removed flags or outdated column names | VERIFIED | `faq.md`: no `--bcftools-filter`, `--chunks`, `--gzip-intermediates` (store_true variant), `--interval-expansion`, or `--no-filtering` references. `association_testing.md`: `skat_pvalue` changed to `skat_o_pvalue` at both locations (lines 334 and 439). `changelog.md`: "classic and stage-based pipeline modes" updated to "pipeline modes". Fisher lambda_GC diagnostic-only note present at lines 548-552. |
| 5 | `make ci-check` passes with zero lint errors and zero test failures | VERIFIED | `make ci-check` completes: 2066 passed, 3 skipped (missing fixture data), 0 failed. Integration tests are excluded from ci-check by design (`not integration` marker). |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|---------|--------|---------|
| `variantcentrifuge/stage_info.py` | DELETED — must not exist | VERIFIED (absent) | File deleted; zero module-level import references anywhere |
| `variantcentrifuge/cli.py` | `coast_backend` choices `["auto", "python"]`, no `--gzip-intermediates store_true` | VERIFIED | Line 439: `choices=["auto", "python"]`. Lines 938-944: only `--no-gzip-intermediates` with `default=True` remains. |
| `variantcentrifuge/inheritance/prioritizer.py` | Only `prioritize_patterns`, `calculate_confidence`, `get_pattern_description` | VERIFIED | 191 lines, 3 public functions. All 6 dead functions absent. |
| `variantcentrifuge/inheritance/analyzer.py` | No `get_inheritance_summary`, `filter_by_inheritance_pattern`, `export_inheritance_report`; `create_inheritance_details` renamed to `_create_inheritance_details` | VERIFIED | Only 4 functions: `analyze_inheritance`, `_finalize_inheritance_patterns`, `_create_inheritance_details`, `process_inheritance_output`. Public dead functions absent. |
| `docs/source/faq.md` | Accurate CLI flag documentation | VERIFIED | No references to removed flags |
| `docs/source/guides/association_testing.md` | `skat_o_pvalue` everywhere, lambda_GC note | VERIFIED | `skat_o_pvalue` at lines 334 and 439. Fisher lambda_GC diagnostic-only note at lines 548-552. |
| `docs/source/changelog.md` | No "classic pipeline" reference | VERIFIED | "pipeline modes" at lines 46 and 51 (no "classic" qualifier) |
| `variantcentrifuge/stages/__init__.py` | Complete `__all__` including 6 previously missing stages | VERIFIED | `AssociationAnalysisStage`, `ClinVarPM5Stage`, `VariantAnalysisStage`, `DataSortingStage`, `ParallelCompleteProcessingStage`, `PCAComputationStage` all present in imports and `__all__`. |
| `variantcentrifuge/pipeline_core/runner.py` | No `TODO: Implement intelligent batching` | VERIFIED | Line 724: `# NOTE: Intelligent batching deferred — current fixed-size batching is sufficient` |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `analyzer.py` | `_create_inheritance_details` | Internal call from `_finalize_inheritance_patterns` | VERIFIED | Line 296: `details = _create_inheritance_details(...)` |
| `cli.py` | `coast_backend` validation | `analysis_stages.py` line 2195 | VERIFIED | Validation checks `not in ("auto", "python")` — no `r` in allowed set |
| `test_performance_optimizations.py` | `_create_inheritance_details` | Import at line 11 | VERIFIED | `from variantcentrifuge.inheritance.analyzer import _create_inheritance_details` — uses private name correctly |

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| DEAD-01: Delete `stage_info.py` | SATISFIED | File absent, zero import references |
| DEAD-02: Remove `--coast-backend r` choice | SATISFIED | CLI and analysis_stages.py validation both updated |
| DEAD-03: Remove 6 dead prioritizer functions | SATISFIED | All 6 absent from `prioritizer.py` |
| DEAD-04: Remove 3 dead analyzer functions; rename `create_inheritance_details` | SATISFIED | 3 dead functions absent; `_create_inheritance_details` is the only form |
| DEAD-05: Remove `--gzip-intermediates store_true` | SATISFIED | Only `--no-gzip-intermediates` with `action="store_false", default=True` |
| DOCS-01: Fix removed CLI flags in `faq.md` | SATISFIED | All removed flags replaced with current equivalents |
| DOCS-02: Fix `skat_pvalue` -> `skat_o_pvalue` in association guide | SATISFIED | Both occurrences updated |
| DOCS-03: Remove "classic pipeline" from changelog | SATISFIED | Entry reads "pipeline modes" |
| CLEAN-01: Remove "refactored pipeline" from test docstrings | SATISFIED | Zero matches in `variantcentrifuge/` and `tests/` |
| CLEAN-02: Remove "original pipeline" from source comments | SATISFIED | Zero matches in `variantcentrifuge/` and `tests/` |
| CLEAN-03: Add missing stage exports to `stages/__init__.py` | SATISFIED | All 6 stages exported |
| CLEAN-04: Remove dead `test_parallel_variant_extraction` xfail | SATISFIED | `TestParallelProcessing` class absent from test file |
| CLEAN-05: Remove dead `chunks` config TODO in `analysis_stages.py` | SATISFIED | Verified absent |
| CLEAN-06: Resolve TODO intelligent batching in `runner.py` | SATISFIED | Replaced with deferred NOTE at line 724 |
| CLEAN-07: Add Fisher lambda_GC diagnostic-only note | SATISFIED | Note present at lines 548-552 of association_testing.md |

### Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `tests/integration/test_pipeline_with_mocked_tools.py` | Stale mock `patch("variantcentrifuge.stages.processing_stages.run_command")` — attribute does not exist in that module | Warning | 5 integration tests ERROR on setup. Pre-existing before phase 38 (confirmed via git history). NOT in `make ci-check` scope (integration marker excluded). |

**Note on integration test failures:** The `test_pipeline_with_mocked_tools.py` failures (`AttributeError: processing_stages has no attribute run_command`) were present before phase 38 began. This file's mock targets the wrong module namespace. This is a pre-existing tech debt item not in the phase-38 requirements scope, and it does not affect `make ci-check` which excludes integration tests.

### Human Verification Required

None. All success criteria are mechanically verifiable.

### Gaps Summary

No gaps. All 5 success criteria from the ROADMAP are satisfied:

1. `stage_info.py` is deleted with zero import references remaining — confirmed absent.
2. `--coast-backend r` is absent from CLI choices — runtime check confirms argument error on invalid value.
3. All 9 dead functions removed (6 from `prioritizer.py`, 3 from `analyzer.py`) — grep returns zero matches for all deleted names.
4. All three docs (`faq.md`, `association_testing.md`, `changelog.md`) are accurate — verified line by line.
5. `make ci-check` passes — 2066 passed, 0 failed.

**Notable execution complexity:** Plan-02 ran concurrently with plan-01 and briefly restored the deleted functions (commit `7a34e70`). Plan-01 then issued a supplemental commit (`b9dc03d`) to re-remove them after plan-02 cleaned the test imports. The final state of the codebase is correct.

---

_Verified: 2026-02-26T18:40:22Z_
_Verifier: Claude (gsd-verifier)_
