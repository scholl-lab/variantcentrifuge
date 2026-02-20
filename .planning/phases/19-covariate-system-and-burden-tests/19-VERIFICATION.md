---
phase: 19-covariate-system-and-burden-tests
verified: 2026-02-20T07:37:35Z
status: passed
score: 5/5 must-haves verified
---

# Phase 19: Covariate System and Burden Tests Verification Report

**Phase Goal:** Users can run logistic or linear regression burden tests with covariate adjustment, and the genotype matrix builder handles all real-world genotype encodings correctly.
**Verified:** 2026-02-20T07:37:35Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Covariate file with rows in different order than VCF produces identical aligned matrix as a file in VCF order | VERIFIED | `load_covariates()` uses `df.reindex(vcf_samples)` + NaN assertion (covariates.py:113-117). Manual test confirmed identical numpy arrays. `test_load_covariates_alignment_produces_identical_results` asserts `np.array_equal`. |
| 2 | Logistic burden test reports OR + 95% CI matching `statsmodels.Logit` on same inputs | VERIFIED | `LogisticBurdenTest.run()` extracts beta/CI at index 1, converts to OR via `exp(beta)`. Manual verification: p-value 0.44392944 matches exactly. `test_logistic_burden_matches_manual_statsmodels` asserts `pytest.approx(rel=1e-6)`. |
| 3 | Linear burden test reports beta + SE matching `statsmodels.OLS` on same inputs | VERIFIED | `LinearBurdenTest.run()` returns `effect_size=beta` (not OR). Manual verification: beta 2.45245672 matches exactly. `test_linear_burden_matches_manual_statsmodels` asserts `pytest.approx(rel=1e-6)`. |
| 4 | Beta(MAF;1,25) weights up-weight rare variants; uniform weights produce all-ones | VERIFIED | `beta_maf_weights()` uses `scipy.stats.beta.pdf` with `np.clip`. Verified: weights[MAF=0.001] > weights[MAF=0.01] > weights[MAF=0.1]. `uniform_weights()` returns `np.ones(n)`. `test_beta_weights_monotonic_decreasing` and `test_uniform_weights_all_ones` confirm. |
| 5 | Genotype strings `1/2`, `./.`, `0\|1` parse correctly — multi-allelic counted as 1 dosage, missing to None, phased summed to 0/1/2 | VERIFIED | `parse_gt_to_dosage()` handles all cases: `1/2`→`(1, True)`, `./.`→`(None, False)`, `0\|1`→`(1, False)`. Missing values imputed to `round(2*MAF)` for binary, `2*MAF` for quantitative. 33 tests confirm all edge cases. |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/association/covariates.py` | `load_covariates()` with sample alignment, one-hot, multicollinearity | VERIFIED | 148 lines, exports `load_covariates`, used in analysis_stages.py |
| `variantcentrifuge/association/genotype_matrix.py` | `parse_gt_to_dosage()`, `build_genotype_matrix()` with imputation | VERIFIED | 257 lines, both functions exported and wired to stage |
| `variantcentrifuge/association/weights.py` | `beta_maf_weights()`, `uniform_weights()`, `get_weights()` | VERIFIED | 139 lines, all three functions exported, used in burden tests |
| `variantcentrifuge/association/base.py` | Extended `AssociationConfig` with 8 new fields | VERIFIED | 198 lines; all 8 fields present with defaults (trait_type, variant_weights, covariate_file, etc.) |
| `variantcentrifuge/association/tests/logistic_burden.py` | `LogisticBurdenTest` with Firth fallback, OR + CI, warning codes | VERIFIED | 436 lines, exports `LogisticBurdenTest`, registered in engine |
| `variantcentrifuge/association/tests/linear_burden.py` | `LinearBurdenTest` with beta + SE output | VERIFIED | 212 lines, exports `LinearBurdenTest`, registered in engine |
| `variantcentrifuge/cli.py` | 5 new CLI args: `--covariate-file`, `--covariates`, `--categorical-covariates`, `--trait-type`, `--variant-weights` | VERIFIED | All 5 args present at lines 430-470; config mapping at lines 1092-1102; validation at lines 1242-1251 |
| `variantcentrifuge/stages/analysis_stages.py` | Covariate loading and genotype matrix building in `AssociationAnalysisStage` | VERIFIED | `load_covariates` called at line 2194; `build_genotype_matrix` called at line 2282; both behind `covariate_file` and `needs_genotype_matrix` guards |
| `tests/unit/test_association_covariates.py` | 17 tests for COV-01..04 | VERIFIED | 366 lines, 17 tests, all passing |
| `tests/unit/test_association_genotype_matrix.py` | 33 tests for BURDEN-03 | VERIFIED | 504 lines, 33 tests, all passing |
| `tests/unit/test_association_weights.py` | 23 tests for WEIGHT-01..02 | VERIFIED | 204 lines, 23 tests, all passing |
| `tests/unit/test_association_logistic_burden.py` | 15 tests for BURDEN-01 + Firth | VERIFIED | 422 lines, 15 tests, all passing |
| `tests/unit/test_association_linear_burden.py` | 11 tests for BURDEN-02 + CI | VERIFIED | 330 lines, 11 tests, all passing |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `covariates.py` | `pandas reindex` | `df.reindex(vcf_samples)` + NaN assertion | WIRED | Lines 113-117 of covariates.py |
| `genotype_matrix.py` | `weights.py` | `variant_mafs` passed to `get_weights()` in burden tests | WIRED | logistic_burden.py:319, linear_burden.py:149 |
| `base.py` | `covariates.py` | `AssociationConfig.covariate_file` drives loading in stage | WIRED | analysis_stages.py:2117, 2191 |
| `engine.py` | `logistic_burden.py` | `_build_registry()` imports and registers `LogisticBurdenTest` | WIRED | engine.py:42-47 |
| `engine.py` | `linear_burden.py` | `_build_registry()` imports and registers `LinearBurdenTest` | WIRED | engine.py:42-48 |
| `analysis_stages.py` | `genotype_matrix.py` | `build_genotype_matrix` called per gene when regression tests requested | WIRED | analysis_stages.py:2269, 2282 |
| `analysis_stages.py` | `covariates.py` | `load_covariates` called once before per-gene loop | WIRED | analysis_stages.py:2192-2194 |
| `cli.py` | `base.py` | CLI args map to `AssociationConfig` fields via `cfg[...]` | WIRED | cli.py:1092-1102 |
| `test_association_logistic_burden.py` | `statsmodels.api.Logit` | Manual reference calculation vs `LogisticBurdenTest` output | WIRED | `sm.Logit` called at line ~150 in test file |
| `test_association_linear_burden.py` | `statsmodels.api.OLS` | Manual reference calculation vs `LinearBurdenTest` output | WIRED | `sm.OLS` called at line ~60 in test file |

### Requirements Coverage

| Requirement | Status | Evidence |
|-------------|--------|----------|
| COV-01: Sample alignment by reindex | SATISFIED | `df.reindex(vcf_samples)`, ValueError for missing samples, test confirms shuffled == ordered |
| COV-02: Categorical one-hot encoding | SATISFIED | `pd.get_dummies(drop_first=True, dtype=float)`, auto-detect + explicit modes |
| COV-03: Multicollinearity warning | SATISFIED | `np.linalg.cond(X) > 1000` triggers `logger.warning()`, test confirms warning via caplog |
| COV-04: Column selection | SATISFIED | `covariate_columns` parameter filters columns before alignment |
| BURDEN-01: Logistic OR + CI matching statsmodels | SATISFIED | p-value, beta, OR all match to rel=1e-6 in `test_logistic_burden_matches_manual_statsmodels` |
| BURDEN-02: Linear beta + SE matching statsmodels | SATISFIED | beta, p-value match to rel=1e-6 in `test_linear_burden_matches_manual_statsmodels` |
| BURDEN-03: Correct genotype encoding (1/2, ./., 0\|1) | SATISFIED | `parse_gt_to_dosage()` handles all edge cases; 13 direct parse tests + matrix tests |
| WEIGHT-01: Beta(MAF;1,25) upweights rare variants | SATISFIED | Monotonically decreasing over MAF grid; 8 weight tests confirm |
| WEIGHT-02: Uniform weights = all 1.0 | SATISFIED | `np.ones(n)` returned; test confirms `gmat @ uniform == gmat.sum(axis=1)` |

### Anti-Patterns Found

No stub patterns, TODOs, FIXMEs, or placeholder content found in any Phase 19 source files. No empty handlers. No console.log-only implementations.

### Human Verification Required

None. All key behaviors are verified programmatically:
- Numerical correctness proven against reference `statsmodels` calculations
- Alignment correctness proven by `np.array_equal` on shuffled vs ordered inputs
- All 99 Phase 19 unit tests pass
- All 72 Phase 18 regression tests pass (no regressions)

## Test Summary

- Phase 19 tests: **99 passed** (covariates: 17, genotype_matrix: 33, weights: 23, logistic_burden: 15, linear_burden: 11)
- Phase 18 regression: **72 passed** (base: 9, engine: 12, fisher: 22, stage: 21 + 8)
- Total at phase completion: 1,441 tests (per 19-03-SUMMARY.md)

---

_Verified: 2026-02-20T07:37:35Z_
_Verifier: Claude (gsd-verifier)_
