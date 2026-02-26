---
phase: 31-coast-fix
verified: 2026-02-23T18:32:11Z
status: passed
score: 5/5 must-haves verified
gaps: []
post_verification_fix: "Added vcf_fields list to variable_assignment_config.json, updated cli.py injection to read from vcf_fields, made analysis_stages.py validation model-aware. Commit 5641221."
real_data_validation: "Tested on GCKD cohort (22GB VCF, ~5000 samples). 3 runs: 7-gene sift_polyphen (4 genes, 0 None), 7-gene CADD (7 genes, 0 None), 500-gene sift_polyphen (444 genes, 0 None). Found and fixed auto-field injection delimiter mismatch (space vs comma in fields_to_extract)."
---

# Phase 31: COAST Fix Verification Report

**Phase Goal:** COAST produces valid p-values for real-world cohorts where most genes lack all three variant categories (BMV/DMV/PTV), genotype matrices are reliably available to the test, multi-transcript SnpEff effect strings are handled correctly, and classification scoring is configurable via the existing formula engine.
**Verified:** 2026-02-23T18:32:11Z
**Status:** passed (after orchestrator fix for auto-field injection gap)
**Re-verification:** No — initial verification + inline fix

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|---------|
| 1 | COAST returns a numeric p-value for genes with only 1 or 2 of 3 variant categories | VERIFIED | Both Python and R backends use `len(present) == 0` guard; 7 tests in `test_coast_partial_category.py` pass |
| 2 | COAST classification correctly handles SnpEff "&"-concatenated effect strings | VERIFIED | `_resolve_effect()` at allelic_series.py:106 implements PTV > missense > first priority; 5 dedicated tests + `test_classify_variants_with_resolve_effect` all pass |
| 3 | SIFT and PolyPhen fields auto-injected when --association-tests coast specified | FAILED | Auto-injection code exists but produces empty field list for built-in models; no actual VCF fields are injected |
| 4 | User can specify --coast-classification cadd and pipeline uses CADD thresholds | VERIFIED | CLI option exists (cli.py:436), resolves path, propagates to AssociationConfig; CADD model tested in test_classify_variants_cadd_config |
| 5 | scoring/coast_classification/ config exists as authoritative classification source | VERIFIED | Three models exist (sift_polyphen, cadd, custom_template), each with valid formula_config.json and variable_assignment_config.json |

**Score:** 4/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/stages/analysis_stages.py` | GT matrix fallback (COAST-01) | VERIFIED | Lines 2376-2397: fallback_df -> gt_cols_df from df; 3875 lines total |
| `variantcentrifuge/association/tests/allelic_series_python.py` | Partial-category fallback, coast_status, finalize | VERIFIED | 508 lines; _n_complete/_n_partial/_n_skipped counters; finalize() logs INFO summary |
| `variantcentrifuge/association/tests/allelic_series.py` | _resolve_effect, partial-category fallback, formula engine path | VERIFIED | 1026 lines; _resolve_effect at line 106; _classify_via_formula_engine; classify_variants with model_dir |
| `variantcentrifuge/association/base.py` | coast_classification field on AssociationConfig | VERIFIED | Line 177: `coast_classification: str \| None = None` |
| `variantcentrifuge/cli.py` | --coast-classification CLI option and auto-field injection | PARTIAL | CLI option verified (line 436); auto-injection block exists (lines 1321-1378) but injects zero fields for built-in models |
| `scoring/coast_classification/sift_polyphen/formula_config.json` | SIFT+PolyPhen classification model | VERIFIED | 32 lines; outputs coast_category; replicates hardcoded logic |
| `scoring/coast_classification/cadd/formula_config.json` | CADD threshold classification model | VERIFIED | 20 lines; CADD>=15 -> DMV, 0<CADD<15 -> BMV |
| `scoring/coast_classification/custom_template/` | Template for custom models | VERIFIED | Both JSON files exist with annotated placeholder structure |
| `tests/unit/test_coast_partial_category.py` | Tests for partial-category behavior | VERIFIED | 246 lines; 7 tests all passing |
| `tests/unit/test_genotype_matrix_fallback.py` | Tests for GT matrix recovery | VERIFIED | 238 lines; 11 tests all passing |
| `tests/unit/test_coast_classification.py` | Tests for effect resolution, formula engine, auto-injection | PARTIAL | 325 lines; 15 tests pass, but test_auto_field_injection_sift_polyphen only asserts isinstance(vcf_fields, list) — does not assert that any fields are actually injected |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `allelic_series_python.py` | `allelic_series.py:classify_variants` | model_dir from config.coast_classification | VERIFIED | Line 308: `coast_model_dir = getattr(config, "coast_classification", None)` → passed to classify_variants |
| `allelic_series.py:classify_variants` | `variantcentrifuge/scoring.py` | model_dir != None → _classify_via_formula_engine → read_scoring_config + apply_scoring | VERIFIED | Lines 324-332 delegate to _classify_via_formula_engine when model_dir provided |
| `variantcentrifuge/cli.py` | `variantcentrifuge/association/base.py` | coast_classification path in cfg | VERIFIED | Line 1377: cfg["coast_classification"] = _coast_model_dir; analysis_stages.py:2201 propagates to AssociationConfig |
| `variantcentrifuge/cli.py` | `scoring/coast_classification/*/variable_assignment_config.json` | auto-field injection reads config keys | PARTIAL | Reads config (line 1346-1350) but filters all COAST_* keys → injects zero VCF fields |
| `analysis_stages.py:GT-fallback` | `allelic_series_python.py` via genotype_matrix | df_with_per_sample_gt → gene_data["genotype_matrix"] | VERIFIED | Lines 2539-2541 build genotype matrix from gt_source_df |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| COAST-01: GT matrix available when variants_df is None | SATISFIED | Fixed in analysis_stages.py; 11 tests verify |
| COAST-02: SIFT/PolyPhen auto-injected when COAST selected | NOT SATISFIED | Auto-injection produces empty field list for built-in models; users must still add fields manually |
| COAST-03: Valid p-value with 1 or 2 of 3 categories | SATISFIED | Both backends fixed; 7 tests verify |
| COAST-04: "&"-concatenated effects handled | SATISFIED | _resolve_effect() implemented; 5 tests verify |
| COAST-05: --coast-classification CLI option | SATISFIED | CLI option, path resolution, AssociationConfig propagation verified |
| COAST-06: scoring/coast_classification/ as authoritative config | SATISFIED | Three model configs exist and parse correctly |
| COAST-07: Multi-transcript priority-order resolution | SATISFIED | _resolve_effect implements PTV > missense > first; 5+ tests verify (NOTE: ROADMAP said majority-vote but PLAN/RESEARCH refined to priority-order) |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `variantcentrifuge/cli.py` | 1354-1355 | `_vcf_fields = [f for f in required_fields if not f.startswith("COAST_")]` returns empty list for built-in models | Blocker | COAST-02 success criterion unmet: no VCF fields injected for default sift_polyphen model |
| `tests/unit/test_coast_classification.py` | 257-259 | `test_auto_field_injection_sift_polyphen` asserts only `isinstance(vcf_fields, list)`, not that any fields are injected | Warning | Test passes when producing empty list, which is the broken behavior |
| `variantcentrifuge/stages/analysis_stages.py` | 2323-2352 | Hardcoded COAST column check still instructs users to add fields manually; not updated to use model-derived field list | Blocker | When SIFT/PolyPhen missing, COAST is skipped with error message "Add them to --fields" |

### Human Verification Required

None identified — all items were verifiable programmatically.

### Gaps Summary

One gap blocks full goal achievement:

**COAST-02 — Auto-field injection does not inject any fields for built-in models.** The implementation made an intentional design choice to use `COAST_*` normalized column names in `variable_assignment_config.json` files (rather than raw VCF field names like `dbNSFP_SIFT_pred`). This correctly avoids sanitization ambiguity inside `classify_variants()`, but it means the auto-injection block in `cli.py` (which filters out `COAST_*` keys) produces an empty `_vcf_fields` list and injects nothing.

The consequence: when a user runs `--association-tests coast` with the default `--coast-classification sift_polyphen` and has NOT added `dbNSFP_SIFT_pred` and `dbNSFP_Polyphen2_HDIV_pred` to `--fields`, the pipeline reaches `analysis_stages.py` line 2342 and silently drops COAST from the test list while printing an error instructing the user to add those fields manually. This is the exact UX problem COAST-02 was supposed to fix.

The four items that ARE working (partial-category fallback, effect string resolution, formula engine classification, CLI option with CADD model) are all fully implemented and tested (142 COAST tests pass). The gap is isolated to the SIFT/PolyPhen auto-injection path.

**Root cause of gap:** The variable_assignment_config design deviation documented in the SUMMARY was architecturally correct for the formula engine but broke the auto-injection mechanism. A supplemental fix is needed: either (a) add a separate `vcf_fields` section to variable_assignment_config specifying which raw VCF fields to inject, or (b) hardcode the per-model injection in cli.py, or (c) update the hardcoded analysis_stages.py validation to not require SIFT/PolyPhen when a formula-engine model is configured.

---

_Verified: 2026-02-23T18:32:11Z_
_Verifier: Claude (gsd-verifier)_
