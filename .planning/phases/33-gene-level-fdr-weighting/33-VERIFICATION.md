---
phase: 33-gene-level-fdr-weighting
verified: 2026-02-24T08:06:17Z
status: passed
score: 8/8 must-haves verified
---

# Phase 33: Gene-Level FDR Weighting Verification Report

**Phase Goal:** Users with biological prior knowledge (pLI scores, GWAS signals) can increase statistical power by providing per-gene weights for weighted Benjamini-Hochberg FDR correction, with correct weight normalization and transparent diagnostics.
**Verified:** 2026-02-24T08:06:17Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Passing --gene-prior-weights applies weighted BH with weights renormalized to mean=1.0; q-values differ from unweighted BH in proportion to weight variation | VERIFIED | `apply_weighted_correction()` at correction.py:141-264 performs `w_norm = raw_weights * m / weight_sum` (mean=1.0 renorm); `test_non_uniform_weights_differ_from_unweighted` and `test_weights_renormalized_to_mean_one` both PASS |
| 2 | Genes absent from the weight file receive weight=1.0 (neutral) | VERIFIED | correction.py:198-204: `raw_weights[i] = 1.0` default in loop for genes not in `weight_map`; `test_missing_genes_get_weight_one` PASSES |
| 3 | If more than 50% of tested genes are missing from the weight file, a warning is emitted | VERIFIED | correction.py:207-218: `missing_frac > 0.50` branch emits `logger.warning`; `test_coverage_warning_50_percent` PASSES |
| 4 | If more than 80% missing, a stronger warning is emitted | VERIFIED | correction.py:208-213: `missing_frac > 0.80` branch emits `logger.warning` with "STRONG WARNING" prefix; `test_coverage_warning_80_percent` PASSES |
| 5 | Effective number of tests (sum(w)^2 / sum(w^2)) appears in diagnostics output | VERIFIED | correction.py:243-248: `n_eff` computed and logged in `apply_weighted_correction()`; diagnostics.py:552-554: also computed and logged via `write_fdr_weight_diagnostics()`; `test_effective_number_of_tests_logged` PASSES |
| 6 | Running without --gene-prior-weights produces identical results to current plain BH | VERIFIED | engine.py:504-505: `else` branch calls `apply_correction(raw_pvals, self._config.correction_method)` — unchanged plain BH path; `test_uniform_weights_match_apply_correction[fdr]` and `[bonferroni]` both PASS; `TestApplyCorrectionBackwardCompat` class (4 tests) all PASS |
| 7 | Zero or negative weights cause an immediate error naming offending genes | VERIFIED | correction.py:122-135: `load_gene_weights()` collects `bad_genes` list and raises `ValueError` with offending gene names; `test_zero_weight_raises`, `test_negative_weight_raises`, `test_multiple_bad_weights_listed` all PASS |
| 8 | Single gene tested skips weighted BH with info log | VERIFIED | correction.py:193-195: `if m == 1: logger.info("Single gene tested — FDR weighting has no effect, skipping"); return pvals_array.copy(), default_norm`; `test_single_gene_returned_unchanged` and `test_single_gene_logs_info` PASS |

**Score:** 8/8 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/association/correction.py` | `load_gene_weights()` and `apply_weighted_correction()` | VERIFIED | File: 265 lines. Both functions present with full implementation. `def load_gene_weights` at line 73, `def apply_weighted_correction` at line 141 |
| `variantcentrifuge/association/base.py` | `gene_prior_weights` and `gene_prior_weight_column` fields on AssociationConfig | VERIFIED | base.py:182-188: `gene_prior_weights: str | None = None` and `gene_prior_weight_column: str = "weight"` under "# Phase 33" comment |
| `variantcentrifuge/association/diagnostics.py` | `write_fdr_weight_diagnostics()` function | VERIFIED | diagnostics.py:438-556: full 119-line implementation writing `fdr_weight_diagnostics.tsv`, computing effective N, logging coverage stats |
| `variantcentrifuge/association/engine.py` | Weighted correction path in run_all FDR pass | VERIFIED | engine.py:443-508: `fdr_weights_by_gene` dict, `use_weighted` branch, `load_gene_weights()` + `apply_weighted_correction()` calls, diagnostics wired. `fdr_weight` column at line 569 |
| `variantcentrifuge/cli.py` | --gene-prior-weights and --gene-prior-weight-column CLI flags | VERIFIED | cli.py:603-617: both `--gene-prior-weights` and `--gene-prior-weight-column` added to stats_group with file-existence validation at cli.py:1304-1308 |
| `variantcentrifuge/stages/analysis_stages.py` | Weight loading, fdr_weight column, diagnostics wiring | VERIFIED | analysis_stages.py: `gene_prior_weights` in `VALID_ASSOCIATION_KEYS` (line 1976), in str validation set (line 2022), and passed to `AssociationConfig` at line 2207-2208 |
| `tests/unit/test_weighted_correction.py` | Unit tests for weighted BH correction and weight loading (min 80 lines) | VERIFIED | 470 lines, 34 tests across 4 classes: `TestLoadGeneWeights`, `TestApplyWeightedCorrection`, `TestApplyCorrectionBackwardCompat`, `TestWriteFdrWeightDiagnostics`. All 34 PASS |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `variantcentrifuge/cli.py` | `variantcentrifuge/stages/analysis_stages.py` | `cfg['gene_prior_weights']` set by CLI, read by `_build_assoc_config_from_context` | WIRED | cli.py:1307 sets `cfg["gene_prior_weights"]`; analysis_stages.py:2207 reads it via `_get("gene_prior_weights", ...)` into AssociationConfig |
| `variantcentrifuge/stages/analysis_stages.py` | `variantcentrifuge/association/correction.py` | `load_gene_weights()` called in AssociationAnalysisStage.run() | WIRED | engine.py:457-466 imports and calls `load_gene_weights()` when `use_weighted` is True; the engine is invoked by `AssociationAnalysisStage` |
| `variantcentrifuge/association/engine.py` | `variantcentrifuge/association/correction.py` | `apply_weighted_correction()` called in FDR pass when weights present | WIRED | engine.py:457-459 imports `apply_weighted_correction`, engine.py:473-475 calls it with `(raw_pvals, testable_genes, weight_map, method)` |
| `variantcentrifuge/association/engine.py` | `variantcentrifuge/association/diagnostics.py` | `fdr_weight` column added to results_df rows | WIRED | engine.py:567-569: `if fdr_weights_by_gene: row["fdr_weight"] = fdr_weights_by_gene.get(gene, 1.0)`; engine.py:485-503 imports and calls `write_fdr_weight_diagnostics()` |

### Requirements Coverage (from ROADMAP.md success criteria)

| Requirement | Status | Notes |
|-------------|--------|-------|
| Weighted BH with mean=1.0 renormalization, q-values differ from plain BH | SATISFIED | Verified in correction.py + 3 unit tests |
| Genes absent get weight=1.0; >50% missing warning | SATISFIED | Both coverage thresholds implemented and tested |
| Effective number of tests in diagnostics | SATISFIED | Appears in both correction.py logger.info and diagnostics.py TSV log |
| Without --gene-prior-weights: identical to current plain BH | SATISFIED | Engine else-branch is unchanged; backward-compat test class passes |
| IHW is NOT implemented — no CLI flag, no stub, no error message | SATISFIED | No "ihw" or "IHW" string anywhere in variantcentrifuge/ source |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | — | — | — | No TODO/FIXME, no placeholder content, no empty handlers found in phase artifacts |

### Human Verification Required

None. All truths can be verified structurally. The mathematical correctness of the weighted BH implementation (Genovese et al. 2006 method) is validated by unit tests, including the `test_uniform_weights_match_apply_correction` test confirming backward compatibility and `test_non_uniform_weights_differ_from_unweighted` confirming the weighting has measurable effect.

### Gaps Summary

No gaps. All 8 truths are verified against actual code:

- `correction.py` contains substantive, fully-implemented `load_gene_weights()` and `apply_weighted_correction()` with correct renormalization math, coverage warnings at 50%/80% thresholds, effective-N computation, and single-gene edge case handling.
- `base.py` carries the two new `AssociationConfig` fields.
- `diagnostics.py` has `write_fdr_weight_diagnostics()` producing a structured TSV with per-gene raw weight, normalized weight, unweighted vs weighted q-values, and significance-change classification.
- `engine.py` has a clean `use_weighted` branch that loads weights, runs both weighted and unweighted correction (for diagnostics comparison), writes diagnostics, and adds the `fdr_weight` output column.
- `cli.py` adds both CLI flags with file-existence validation.
- `analysis_stages.py` passes the new fields through to AssociationConfig.
- 34 unit tests all pass, covering every must-have truth.
- IHW is cleanly absent from the codebase as specified.

---

_Verified: 2026-02-24T08:06:17Z_
_Verifier: Claude (gsd-verifier)_
