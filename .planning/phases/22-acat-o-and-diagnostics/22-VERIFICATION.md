---
phase: 22-acat-o-and-diagnostics
verified: 2026-02-21T18:27:35Z
status: gaps_found
score: 3/5 must-haves verified
gaps:
  - truth: "The per-gene association TSV contains all standard columns including skat_o_p, and explicit carrier counts"
    status: failed
    reason: "There is no distinct skat_o_p_value column. When SKAT runs with skat_method=SKATO the p-value appears in skat_p_value. DIAG-01 lists skat_o_p as a required column name separate from skat_p. Additionally, carrier counts (proband_carrier_count, control_carrier_count) are not explicit engine output columns — they appear only as embedded contingency table in Fisher's 'table' extra column."
    artifacts:
      - path: "variantcentrifuge/association/engine.py"
        issue: "No skat_o_p_value column produced. Engine uses {test_name}_p_value pattern so SKAT-O method produces skat_p_value, not skat_o_p_value. No proband_carrier_count or control_carrier_count in output."
      - path: "variantcentrifuge/association/tests/skat_python.py"
        issue: "test.name returns 'skat' regardless of skat_method (SKAT vs SKATO vs Burden). So skat_method=SKATO produces skat_p_value, not skat_o_p_value."
    missing:
      - "Either a dedicated skat_o_p_value column when skat_method=SKATO, or documentation that skat_p_value serves as skat_o_p when SKATO method is configured"
      - "Explicit proband_carrier_count and control_carrier_count columns in engine output (or confirmation that 'table' extra satisfies DIAG-01)"

  - truth: "The --diagnostics-output directory contains lambda_gc.txt (criterion specifies .txt extension)"
    status: failed
    reason: "ROADMAP success criterion #3 specifies lambda_gc.txt but the implementation creates lambda_gc.tsv. The CONTEXT.md decision specifies .tsv, so this may be a criterion typo, but as written the criterion is not met."
    artifacts:
      - path: "variantcentrifuge/association/diagnostics.py"
        issue: "write_diagnostics() creates lambda_gc.tsv at line 299, not lambda_gc.txt as criterion #3 requires"
    missing:
      - "Either rename lambda_gc.tsv to lambda_gc.txt, or update ROADMAP criterion #3 to say lambda_gc.tsv"

  - truth: "Lambda_GC computed on a permuted null phenotype falls within [0.95, 1.05] for Fisher, burden, and SKAT tests"
    status: failed
    reason: "Fisher's exact test is inherently conservative — p-values are sub-uniform under the null. Actual lambda_GC on a permuted null with Fisher test is ~0.84-0.85, not within [0.95, 1.05]. The test suite uses [0.85, 1.15] tolerance for uniform null p-values, not the [0.95, 1.05] criterion. No test validates the [0.95, 1.05] criterion on a permuted null phenotype with the actual association test machinery."
    artifacts:
      - path: "tests/unit/test_association_diagnostics.py"
        issue: "test_lambda_gc_uniform uses [0.85, 1.15] tolerance on 1000 uniform p-values, not [0.95, 1.05] on permuted null. No end-to-end permuted null test exists."
      - path: "variantcentrifuge/association/tests/fisher.py"
        issue: "Fisher's exact test is inherently conservative; its p-values are not uniformly distributed under null, so lambda_GC < 1.0 is expected behavior, not a calibration bug."
    missing:
      - "Clarification in ROADMAP whether criterion #5 applies to Fisher's exact test (which is inherently conservative and cannot have lambda_GC near 1.0)"
      - "OR: a permuted null test specifically for score-based tests (burden/SKAT) that CAN achieve lambda_GC in [0.95, 1.05]"
---

# Phase 22: ACAT-O and Diagnostics Verification Report

**Phase Goal:** Users receive omnibus ACAT-O p-values combining burden and SKAT results per gene, a single FDR-corrected set of ACAT-O p-values across all genes, and a diagnostics directory containing lambda_GC per test and QQ plot data TSV.
**Verified:** 2026-02-21T18:27:35Z
**Status:** gaps_found
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Per-gene TSV has all standard columns (fisher_p, burden_p, skat_p, skat_o_p, acat_o_p, effect sizes, CIs, variant/carrier counts) | FAILED | No skat_o_p_value column; carrier counts not explicit; see gap detail |
| 2 | ACAT-O computed via single Cauchy combination per gene, single FDR pass across all genes | VERIFIED | engine._compute_acat_o() post-loop confirmed; apply_correction on acat_o only; 91 tests pass |
| 3 | --diagnostics-output creates lambda_gc.txt (criterion) / lambda_gc.tsv (implementation), qq_data.tsv, summary.txt | PARTIAL | lambda_gc.tsv exists (not .txt as criterion requires), qq_data.tsv and summary.txt confirmed |
| 4 | Pipeline warns on case_carriers < 10, n_cases < 200, ratio > 1:20, flags genes in TSV | VERIFIED | compute_per_gene_warnings(), emit_sample_size_warnings() confirmed; warnings column in results_df |
| 5 | Lambda_GC on permuted null within [0.95, 1.05] for Fisher, burden, SKAT | FAILED | Fisher's exact test is conservative; lambda_GC ~0.84 on null, not within [0.95, 1.05] |

**Score:** 2/5 truths fully verified (truth 2 and 4); 1/5 partially verified (truth 3)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/association/tests/acat.py` | cauchy_combination(), compute_acat_o() | VERIFIED | 167 lines, substantive, Liu & Xie formula with overflow guard |
| `variantcentrifuge/association/diagnostics.py` | compute_lambda_gc(), compute_qq_data(), emit_sample_size_warnings(), compute_per_gene_warnings(), write_diagnostics() | VERIFIED | 357 lines, all 5 functions present and substantive |
| `variantcentrifuge/association/engine.py` | _compute_acat_o() post-loop, ARCH-03 FDR | VERIFIED | 406 lines; _compute_acat_o() at line 165; single apply_correction() on acat_o |
| `variantcentrifuge/association/base.py` | AssociationConfig with min_cases, max_case_control_ratio, min_case_carriers, diagnostics_output | VERIFIED | All 4 fields present at lines 143-153 |
| `variantcentrifuge/cli.py` | --diagnostics-output CLI arg | VERIFIED | Lines 465-473; config mapping at line 1112; validation at line 1256 |
| `variantcentrifuge/stages/analysis_stages.py` | diagnostics integration in AssociationAnalysisStage._process() | VERIFIED | diagnostics_output config at line 2134; per-gene warnings at line 2382; write_diagnostics at line 2422 |
| `tests/unit/test_association_acat.py` | ACAT formula tests (27 tests) | VERIFIED | 240 lines, 27 tests pass |
| `tests/unit/test_association_diagnostics.py` | Diagnostics tests (40 tests) | VERIFIED | 603 lines, 40 tests pass |
| `tests/unit/test_association_engine.py` | Engine ACAT-O integration tests (7 new tests) | VERIFIED | 501 lines, 7 ACAT-O tests pass |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `engine.py:_compute_acat_o()` | `tests/acat.py:compute_acat_o()` | lazy import at line 190 | VERIFIED | `from variantcentrifuge.association.tests.acat import compute_acat_o` confirmed |
| `engine.py:run_all()` | `correction.py:apply_correction()` | single call on acat_o p-values | VERIFIED | apply_correction() called once at line 332, on testable_genes from acat_o_results only |
| `analysis_stages.py` | `diagnostics.py` | lazy import after engine.run_all() | VERIFIED | Lines 2382 and 2417 lazy import `from ..association.diagnostics import ...` |
| `cli.py:--diagnostics-output` | `analysis_stages.py:assoc_config` | config key `diagnostics_output` | VERIFIED | cfg["diagnostics_output"] at line 1112; AssociationConfig(diagnostics_output=...) at line 2134 |

### Requirements Coverage

| Requirement | Description | Status | Blocking Issue |
|-------------|-------------|--------|----------------|
| OMNI-01 | ACAT-V per-variant Cauchy combination | OUT OF SCOPE | Explicitly deferred to Phase 23 per CONTEXT.md; cauchy_combination() infrastructure created |
| OMNI-02 | ACAT-O omnibus combining burden + SKAT (+ ACAT-V) per gene | PARTIAL | ACAT-O implemented for burden+SKAT; ACAT-V input deferred to Phase 23 |
| OMNI-03 | Single FDR on ACAT-O across genes | VERIFIED | ARCH-03 implemented and tested |
| DIAG-01 | Per-gene TSV with skat_o_p, all standard columns | PARTIAL | No distinct skat_o_p_value column; carrier counts not explicit |
| DIAG-02 | Lambda_GC per test | VERIFIED | compute_lambda_gc() per test in write_diagnostics() |
| DIAG-03 | QQ plot data TSV | VERIFIED | qq_data.tsv with test, expected_neg_log10_p, observed_neg_log10_p |
| DIAG-05 | --diagnostics-output CLI arg | VERIFIED | Present, validated, maps to config |
| DIAG-06 | Warn n_cases < 200, ratio > 1:20, flag genes with case_carriers < 10 | VERIFIED | emit_sample_size_warnings() and compute_per_gene_warnings() confirmed |

### Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `variantcentrifuge/association/diagnostics.py` | None found | - | - |
| `variantcentrifuge/association/tests/acat.py` | None found | - | - |
| `variantcentrifuge/association/engine.py` | None found | - | - |

No blocker anti-patterns (placeholder implementations, empty handlers, TODO stubs) found in any Phase 22 artifacts.

### ACAT-O Numerical Verification

Direct code execution confirmed:
- `cauchy_combination([0.05, 0.05, 0.05])` = 0.050000 (equal-input identity, expected ~0.05)
- `cauchy_combination([0.001, 0.01, 0.05])` = 0.002679 (matches Liu & Xie 2020 Table 1 reference value ~0.00268)
- `compute_acat_o({"burden": 0.05, "skat": None})` = 0.05 (pass-through)
- `compute_acat_o({"burden": None, "skat": None})` = None
- `cauchy_combination([1e-20, 0.5])` = valid float (overflow guard working)
- Lambda_GC on 1000 uniform(0,1) null p-values = 1.0150 (within [0.9, 1.1])

### Gaps Summary

Three gaps block full criterion compliance:

**Gap 1 (DIAG-01 / Criterion 1): No distinct skat_o_p_value column.** DIAG-01 and criterion #1 list `skat_o_p` as a required column name. The implementation uses `skat_p_value` for the SKAT test p-value regardless of whether `skat_method` is SKAT, SKATO, or Burden. This is a naming gap: the column that serves as `skat_o_p` when `skat_method=SKATO` is named `skat_p_value`, not `skat_o_p_value`. This gap is primarily a specification vs implementation mismatch — the p-value is computed and stored, just under a different column name than DIAG-01 specifies. Additionally, carrier counts (proband/control) are not explicit output columns.

**Gap 2 (Criterion 3): lambda_gc.tsv vs lambda_gc.txt.** The ROADMAP criterion says `lambda_gc.txt` but the implementation creates `lambda_gc.tsv`. The CONTEXT.md decision explicitly chose `.tsv`. This is likely a typo in the ROADMAP criterion, but as written it fails the criterion. Resolution: update ROADMAP criterion #3 to say `lambda_gc.tsv`.

**Gap 3 (Criterion 5): lambda_GC [0.95, 1.05] on permuted null.** Fisher's exact test is inherently conservative — its p-values are sub-uniform under the null, and lambda_GC consistently falls around 0.84-0.85, not within [0.95, 1.05]. The criterion as written cannot be met for Fisher's exact test. Score-based tests (logistic burden, SKAT) may meet this criterion but no end-to-end test validates it. Resolution options: (a) restrict criterion #5 to score-based tests only, excluding Fisher; (b) add end-to-end permuted null test for burden/SKAT tests; (c) clarify that criterion #5 is a calibration property of score tests, not Fisher.

**Note on OMNI-01/OMNI-02 scope:** REQUIREMENTS.md still marks OMNI-01 and OMNI-02 as "Pending". The CONTEXT.md scopes ACAT-V out of Phase 22. The cauchy_combination() infrastructure for ACAT-V is ready; the per-variant application is Phase 23 scope.

---

*Verified: 2026-02-21T18:27:35Z*
*Verifier: Claude (gsd-verifier)*
