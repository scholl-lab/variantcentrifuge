---
phase: 29-classic-pipeline-deprecation-and-removal
verified: 2026-02-23T06:52:11Z
status: passed
score: 5/5 must-haves verified
gaps: []
---

# Phase 29: Classic Pipeline Deprecation and Removal Verification Report

**Phase Goal:** Clean up vestigial "classic pipeline" and "refactored" naming from the codebase. Research confirmed pipeline.py IS the stage-based pipeline (no classic code path remains). This phase removes the phantom `--use-new-pipeline` flag from tests/docs, renames `run_refactored_pipeline` to `run_pipeline`, and updates documentation to describe a single pipeline architecture.
**Verified:** 2026-02-23T06:52:11Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `run_refactored_pipeline` renamed to `run_pipeline` in pipeline.py; cli.py imports and calls `run_pipeline` | VERIFIED | `pipeline.py:369: def run_pipeline(args: argparse.Namespace) -> None:`; `cli.py:12: from .pipeline import run_pipeline`; `cli.py:1665: run_pipeline(refactored_args)` |
| 2 | Zero references to `--use-new-pipeline` in any test file, fixture, or documentation | VERIFIED | `grep -rn "use.new.pipeline" tests/` returns nothing; `grep -in "use.new.pipeline" README.md CLAUDE.md` returns nothing |
| 3 | pipeline.py docstring describes stage-based pipeline (no "refactored" or "classic" language) | VERIFIED | Module docstring is "Stage-based pipeline orchestration. This module builds and runs the stage-based pipeline using PipelineRunner, PipelineContext, and registered Stage subclasses from pipeline_core/." |
| 4 | All existing tests pass after renaming; no behavioral changes | VERIFIED | 2001 passed, 3 skipped (missing fixture data — expected), 0 failures; lint clean |
| 5 | CLAUDE.md, README updated to reflect single pipeline architecture (no "Two Pipeline Modes") | VERIFIED | CLAUDE.md "Two Pipeline Modes" section replaced with "Pipeline Architecture"; README.md `(--use-new-pipeline)` parenthetical removed; `grep -in "two pipeline\|classic pipeline\|use.new.pipeline" CLAUDE.md README.md` returns nothing |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/pipeline.py` | Stage-based pipeline with clean naming; `def run_pipeline` | VERIFIED | 566 lines; `def run_pipeline` at line 369; no `def main`; clean module docstring; no "refactored" in docstring |
| `variantcentrifuge/cli.py` | CLI entry point importing `run_pipeline` | VERIFIED | 1520 lines; `from .pipeline import run_pipeline` at line 12; `run_pipeline(refactored_args)` at line 1665 |
| `tests/integration/test_pipeline_with_mocked_tools.py` | Updated import and call sites | VERIFIED | Line 15: `from variantcentrifuge.pipeline import build_pipeline_stages, run_pipeline`; 6 call sites use `run_pipeline(args)` |
| `tests/test_gene_list_integration.py` | Mock function and monkeypatch target updated | VERIFIED | Lines 140-150: `mock_run_pipeline`, `"variantcentrifuge.cli.run_pipeline"` |
| `tests/unit/test_cli_debug_logging.py` | @patch decorators updated | VERIFIED | Lines 13, 76: `@patch("variantcentrifuge.cli.run_pipeline")` |
| `tests/conftest.py` | Gene burden fixture without `--use-new-pipeline` | VERIFIED | `grep "use.new.pipeline" tests/conftest.py` returns nothing |
| `tests/test_gene_burden_comprehensive.py` | 8 test commands without phantom flag | VERIFIED | `grep "use.new.pipeline" tests/test_gene_burden_comprehensive.py` returns nothing |
| `tests/test_cli.py` | Docstrings without `--use-new-pipeline` | VERIFIED | `grep -in "use.new.pipeline" tests/test_cli.py` returns nothing |
| `tests/fixtures/geneburden/README.md` | Example commands without phantom flag | VERIFIED | `grep -in "use.new.pipeline" tests/fixtures/geneburden/README.md` returns nothing |
| `CLAUDE.md` | Single pipeline architecture description | VERIFIED | "Pipeline Architecture" section (not "Two Pipeline Modes"); describes `pipeline.py` and `pipeline_core/` as the single architecture |
| `README.md` | No legacy pipeline references | VERIFIED | `grep -in "use.new.pipeline\|classic pipeline\|two pipeline" README.md` returns nothing |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `variantcentrifuge/cli.py` | `variantcentrifuge/pipeline.py` | `from .pipeline import run_pipeline` | WIRED | Import at line 12; call at line 1665 |
| `tests/integration/test_pipeline_with_mocked_tools.py` | `variantcentrifuge/pipeline.py` | `from variantcentrifuge.pipeline import run_pipeline` | WIRED | Import line 15; 6 call sites |
| `CLAUDE.md` | `variantcentrifuge/pipeline.py` | Architecture description referencing `pipeline.py` as stage-based entry point | WIRED | "pipeline.py — Builds and runs the stage-based pipeline" at line 36 |

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| `run_refactored_pipeline` renamed to `run_pipeline` in pipeline.py; cli.py imports and calls `run_pipeline` | SATISFIED | Verified in both files |
| Zero references to `--use-new-pipeline` in any test file, fixture, or documentation | SATISFIED | Full codebase scan returns nothing outside `.planning/` dir |
| pipeline.py docstring describes stage-based pipeline (no "refactored" or "classic" language) | SATISFIED | Module docstring at lines 1-6 is clean |
| All existing tests pass after renaming; no behavioral changes | SATISFIED | 2001/2001 fast tests pass |
| CLAUDE.md, README updated to reflect single pipeline architecture (no "Two Pipeline Modes") | SATISFIED | Both files updated |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `variantcentrifuge/pipeline.py` | 509 | `"refactored_pipeline"` as fallback default string for checkpoint version | Info | Not user-facing; never reached in practice — `cli.py` always sets `cfg["pipeline_version"] = __version__` before calling `run_pipeline`; `setup_stages.py` also sets it from `__version__` as a secondary guard |
| `variantcentrifuge/pipeline_core/runner.py` | 119 | `"refactored_pipeline"` as fallback default string for checkpoint version | Info | Same as above — internal checkpoint compatibility identifier, not pipeline naming |

These two `"refactored_pipeline"` strings are internal checkpoint fallback defaults written before the phase, not user-facing pipeline mode naming. They are unreachable in normal operation (cli.py always injects `__version__` first) and exist solely as compatibility guards for hypothetical malformed state files. They do not contradict the phase goal of removing "refactored" naming from the pipeline's public interface, docstrings, or test/documentation references.

### Human Verification Required

None. All success criteria are mechanically verifiable.

### Gaps Summary

No gaps. All five success criteria from the phase goal are fully achieved in the codebase:

1. `run_pipeline` is the canonical function name throughout source and tests — no instances of `run_refactored_pipeline` remain anywhere in the codebase.
2. `--use-new-pipeline` has been completely eliminated from tests, fixtures, and documentation.
3. `pipeline.py` module docstring is clean — describes stage-based pipeline orchestration.
4. The test suite passes without modification — 2001 tests, 0 failures.
5. CLAUDE.md and README.md now describe a single pipeline architecture.

---

_Verified: 2026-02-23T06:52:11Z_
_Verifier: Claude (gsd-verifier)_
