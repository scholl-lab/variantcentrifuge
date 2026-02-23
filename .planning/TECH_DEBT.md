# Tech Debt Backlog

Carried forward from v0.15.0 milestone audit (2026-02-23).

## Pipeline Wiring

### TD-01: PCAComputationStage not wired into pipeline.py

**Priority:** High
**Source:** v0.15.0 audit, cross-phase integration check
**Files:** `variantcentrifuge/pipeline.py`, `variantcentrifuge/stages/processing_stages.py:2075`, `variantcentrifuge/stages/stage_registry.py:481`

`PCAComputationStage` is defined and registered in the stage registry but never imported or added in `pipeline.py:build_pipeline_stages()`. The `--pca-tool akt` CLI flag silently does nothing — AKT never executes, no PCA file is produced, and PCA merge is skipped.

**Workaround:** Use `--pca-file` with a pre-computed PCA eigenvec file.

**Fix:** Import `PCAComputationStage` in `pipeline.py` and conditionally append it when `pca_tool` is set.

---

### TD-02: create_stages_from_config() missing association mapping

**Priority:** Medium
**Source:** v0.15.0 audit, integration checker
**Files:** `variantcentrifuge/pipeline.py:330-366`

`create_stages_from_config()` does not set `args.perform_association` or `args.perform_gene_burden` from the config dict. Programmatic callers (e.g., Snakemake workflows) cannot trigger association analysis via config dict. The primary CLI path is unaffected.

**Fix:** Map `perform_association` and `perform_gene_burden` keys in `create_stages_from_config()`.

---

## Validation Gaps

### TD-03: COAST R golden values not generated

**Priority:** Medium
**Source:** Phase 24 verification (COAST-PY-01)
**Files:** `tests/unit/test_coast_python_comparison.py`, `scripts/generate_coast_golden.R`, `tests/fixtures/coast_golden/`

Tolerance tests in `test_coast_python_comparison.py` compare Python-vs-Python (self-consistency), not Python-vs-R. The R script is written (238 lines, 5 scenarios) and ready to run. Requires R + AllelicSeries package installed to generate reference values and hardcode them as `_R_GOLDEN_*` constants in the test file.

**Fix:** Run `Rscript scripts/generate_coast_golden.R`, copy output constants into test file, replace self-consistency assertions with Python-vs-R comparisons.

---

## Specification Mismatches

### TD-04: skat_o_p_value column naming

**Priority:** Low
**Source:** Phase 22 verification (DIAG-01)
**Files:** `variantcentrifuge/association/engine.py`, `variantcentrifuge/association/tests/skat_python.py`

DIAG-01 specifies `skat_o_p` as a distinct output column. Implementation uses `skat_p_value` for all SKAT methods (SKAT, SKATO, Burden) since `test.name` returns `"skat"` regardless of method. The p-value is computed correctly; only the column name differs from spec. The `skat_method` extra column indicates which method was used.

**Options:** (a) Rename column to `skat_o_p_value` when `skat_method=SKATO`, (b) update DIAG-01 spec to match implementation, (c) accept as-is with documentation note.

---

### TD-05: Fisher lambda_GC criterion

**Priority:** Low
**Source:** Phase 22 verification (criterion #5)
**Files:** `.planning/ROADMAP.md:170`

ROADMAP criterion #5 states lambda_GC should fall in [0.95, 1.05] on permuted null. Fisher's exact test is inherently conservative (sub-uniform p-values under null), yielding lambda_GC ~0.85. The criterion is valid for score-based tests (burden, SKAT) but inapplicable to Fisher.

**Fix:** Clarify ROADMAP criterion to exclude Fisher, or add a permuted-null test specifically for score-based tests.

---

## Minor Cleanup

### TD-06: Internal "refactored_pipeline" checkpoint strings

**Priority:** Low
**Source:** Phase 29 verification
**Files:** `variantcentrifuge/pipeline.py:509`, `variantcentrifuge/pipeline_core/runner.py:119`

Two internal `"refactored_pipeline"` strings remain as checkpoint fallback defaults. They are unreachable in normal operation — `cli.py` always sets `cfg["pipeline_version"] = __version__` before calling `run_pipeline()`. They exist as compatibility guards for malformed state files.

**Fix:** Replace with `"pipeline"` or `"unknown"`.

---

*Created: 2026-02-23 from v0.15.0 milestone audit*
