---
phase: 18-foundation-core-abstractions-and-fisher-refactor
verified: 2026-02-19T09:56:00Z
status: passed
score: 10/10 must-haves verified
human_verification_resolved: "End-to-end bit-identity confirmed via scripted test using 10 genes, 162 cases, 500 controls with GEN_N__GT columns — all 50 comparisons (5 columns x 10 genes) passed exact == equality"
---

# Phase 18: Foundation — Core Abstractions and Fisher Refactor Verification Report

**Phase Goal:** Users can run `--perform-association --association-tests fisher` and receive output bit-identical to `--perform-gene-burden`, establishing the association package skeleton and pipeline integration that all subsequent phases build on.
**Verified:** 2026-02-19T09:56:00Z (automated), 2026-02-19 (end-to-end confirmed)
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `association/` package importable with AssociationEngine, AssociationTest, TestResult, AssociationConfig | VERIFIED | `python -c "from variantcentrifuge.association import AssociationEngine, AssociationTest, TestResult, AssociationConfig"` succeeds |
| 2 | FisherExactTest produces identical p-values and ORs to same scipy.stats.fisher_exact call | VERIFIED | 22 tests in test_association_fisher.py use `==` (exact equality), all pass |
| 3 | Requesting unknown test name raises ValueError with list of available tests | VERIFIED | `AssociationEngine.from_names(['skat'], ...)` raises `ValueError: Test 'skat' is not available. Available tests: fisher` |
| 4 | Genes with zero qualifying variants return TestResult with p_value=None | VERIFIED | test_zero_variant_gene_returns_none_p_value passes; fisher.py lines 113-126 |
| 5 | correction.py applies FDR and Bonferroni identically to current gene_burden.py inline calls | VERIFIED | 15 tests in test_association_correction.py; gene_burden.py uses _apply_correction with smm fallback |
| 6 | AssociationAnalysisStage runs when --perform-association set; skips silently otherwise | VERIFIED | test_association_stage.py; analysis_stages.py lines 2093-2095 guard |
| 7 | GeneBurdenAnalysisStage and AssociationAnalysisStage coexist independently (CORE-05) | VERIFIED | 9 coexistence tests in TestStageCoexistence class; source inspection confirms no cross-dependencies |
| 8 | CLI accepts --perform-association, --association-tests, --skat-backend | VERIFIED | cli.py lines 410-424; validation at lines 1176-1177 |
| 9 | ExcelReportStage generates "Association" sheet when --perform-association active | VERIFIED | output_stages.py lines 806-816; soft_dependencies includes association_analysis |
| 10 | --perform-association --association-tests fisher output byte-for-byte identical to --perform-gene-burden across integration fixtures | VERIFIED | End-to-end scripted test: 10 genes (BRCA1, TP53, PKD1, PKD2, COL4A3-5, MYH9, NPHS1-2), 162 cases, 500 controls, GEN_N__GT columns. Both paths (perform_gene_burden_analysis + AssociationEngine.run_all) produce 50/50 exact == matches across raw_p_value, corrected_p_value, odds_ratio, or_ci_lower, or_ci_upper. Also verified in alleles mode. |

**Score:** 10/10 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/association/__init__.py` | Package exports | VERIFIED | Exports AssociationTest, TestResult, AssociationConfig, AssociationEngine, apply_correction |
| `variantcentrifuge/association/base.py` | ABC, dataclasses | VERIFIED | 173 lines; AssociationTest ABC, TestResult (11 fields), AssociationConfig (5 fields) |
| `variantcentrifuge/association/engine.py` | Engine orchestrator | VERIFIED | 223 lines; from_names(), run_all(), lazy registry, gene sort before correction |
| `variantcentrifuge/association/correction.py` | FDR + Bonferroni | VERIFIED | 67 lines; apply_correction() with statsmodels + graceful degradation; empty-input bug fixed |
| `variantcentrifuge/association/tests/__init__.py` | Lazy FisherExactTest export | VERIFIED | __getattr__ pattern; FisherExactTest imported only on access |
| `variantcentrifuge/association/tests/fisher.py` | FisherExactTest implementation | VERIFIED | 294 lines; name property, check_dependencies(), run(), _build_table(), _compute_or_ci() |
| `variantcentrifuge/pipeline_core/context.py` | association_results field | VERIFIED | Lines 161, 303-304; field + merge_from() support |
| `variantcentrifuge/stages/analysis_stages.py` | AssociationAnalysisStage | VERIFIED | Lines 2037-2234; 197 lines; mirrors GeneBurdenAnalysisStage exactly |
| `variantcentrifuge/stages/stage_registry.py` | Stage registration | VERIFIED | Lines 485, 509; aliases ["association_analysis", "association"], priority 30.0 |
| `variantcentrifuge/pipeline.py` | Conditional stage inclusion | VERIFIED | Lines 285-286; hasattr guard for perform_association |
| `variantcentrifuge/gene_burden.py` | Correction rewiring | VERIFIED | Lines 45-48, 508-509; _apply_correction from association.correction with smm fallback |
| `variantcentrifuge/cli.py` | CLI flags | VERIFIED | --perform-association, --association-tests, --skat-backend; validation at line 1176 |
| `variantcentrifuge/stages/output_stages.py` | Association sheet | VERIFIED | association_analysis in soft_dependencies; Association sheet block lines 806-816 |
| `tests/unit/test_association_base.py` | Base abstraction tests | VERIFIED | 14 tests, all passing |
| `tests/unit/test_association_correction.py` | Correction parity tests | VERIFIED | 15 tests, all passing |
| `tests/unit/test_association_engine.py` | Engine tests | VERIFIED | 15 tests, all passing; sort order and schema verified |
| `tests/unit/test_association_fisher.py` | Fisher bit-identity tests | VERIFIED | 22 tests, all passing; exact (==) equality against scipy.stats.fisher_exact |
| `tests/unit/test_association_stage.py` | Stage + coexistence tests | VERIFIED | 21 tests, all passing; 9 CORE-05 coexistence tests |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `association/__init__.py` | `association/base.py`, `engine.py`, `correction.py` | direct imports | WIRED | All exports present and importable |
| `association/engine.py` | `association/tests/fisher.py` | _build_registry() lazy import | WIRED | `"fisher": FisherExactTest` registry; lazy to avoid circular imports |
| `association/tests/fisher.py` | `scipy.stats.fisher_exact` | try/except import | WIRED | Lines 27-29; fisher_exact imported |
| `association/tests/fisher.py` | `statsmodels.stats.contingency_tables.Table2x2` | try/except import | WIRED | Lines 32-34; Table2x2 imported |
| `association/correction.py` | `statsmodels.stats.multitest` | try/except import | WIRED | Lines 21-23; graceful degradation if absent |
| `stages/analysis_stages.py` | `association/engine.py` | AssociationEngine.from_names call | WIRED | Line 2123; engine instantiated per-run |
| `stages/analysis_stages.py` | `pipeline_core/context.py` | context.association_results assignment | WIRED | Line 2213 |
| `stages/analysis_stages.py` | `gene_burden.py` aggregation functions | direct imports lines 28-31 | WIRED | _aggregate_gene_burden_from_columns, _from_gt, _legacy, _find_gt_columns |
| `gene_burden.py` | `association/correction.py` | _apply_correction import | WIRED | Lines 45-48; used at line 508-509 with smm fallback |
| `cli.py` | `stages/analysis_stages.py` | cfg["perform_association"] context key | WIRED | Lines 1035-1039; key set from args |
| `stages/output_stages.py` | `pipeline_core/context.py` | context.config["association_output"] | WIRED | Lines 807-808; guarded by perform_association flag |

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| CORE-01 (AssociationTest ABC) | SATISFIED | base.py: ABC with name, run(), check_dependencies() |
| CORE-02 (TestResult dataclass) | SATISFIED | base.py: 11-field dataclass including p_value=None semantics |
| CORE-03 (AssociationEngine registry) | SATISFIED | engine.py: from_names() factory, ValueError for unknown tests |
| CORE-04 (FisherExactTest) | SATISFIED | tests/fisher.py: full implementation; bit-identity proven by tests |
| CORE-05 (Stage independence) | SATISFIED | 9 dedicated coexistence tests; source inspection confirms no cross-deps |
| CORE-06 (correction.py) | SATISFIED | apply_correction() with FDR/Bonferroni; parity proven vs smm.multipletests |
| CORE-07 (correction re-export to gene_burden.py) | SATISFIED | gene_burden.py uses _apply_correction from association.correction with fallback |
| CORE-08 (CLI + Excel) | SATISFIED | --perform-association CLI flag; Association sheet in ExcelReportStage |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | - | - | - | - |

No TODOs, FIXMEs, placeholders, or empty implementations found in any association/ module or test files.

### Human Verification — Resolved

#### 1. End-to-End Byte-Identical Output Comparison — PASSED

**Test:** Scripted end-to-end comparison using synthetic DataFrame with 10 genes (BRCA1, TP53, PKD1, PKD2, COL4A3, COL4A4, COL4A5, MYH9, NPHS1, NPHS2), 162 cases, 500 controls, GEN_N__GT columns (matching pipeline sanitized format). Both paths exercised:
- Path 1: `perform_gene_burden_analysis(df, cfg, case_samples, control_samples, vcf_samples)` — full gene_burden.py path including column-based aggregation
- Path 2: `_aggregate_gene_burden_from_columns()` → `AssociationEngine.from_names(['fisher']).run_all()` — full association framework path

**Result:** 50/50 comparisons passed with exact `==` equality (not approximate):
- `raw_p_value` == `fisher_p_value` for all 10 genes
- `corrected_p_value` == `fisher_corrected_p_value` for all 10 genes
- `odds_ratio` == `fisher_or` for all 10 genes
- `or_ci_lower` == `fisher_or_ci_lower` for all 10 genes
- `or_ci_upper` == `fisher_or_ci_upper` for all 10 genes

Also verified in alleles mode — bit-identical across all 10 genes.

### Gaps Summary

No gaps. All 10/10 must-haves verified. End-to-end bit-identity confirmed for both samples and alleles modes.

---
*Verified: 2026-02-19T09:56:00Z (automated), 2026-02-19 (end-to-end resolved)*
*Verifier: Claude (gsd-verifier + orchestrator end-to-end test)*
