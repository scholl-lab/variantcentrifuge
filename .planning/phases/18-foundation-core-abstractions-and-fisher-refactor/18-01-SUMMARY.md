---
phase: 18-foundation-core-abstractions-and-fisher-refactor
plan: "01"
subsystem: association
tags: [scipy, statsmodels, fisher-exact, abc, dataclass, multiple-testing, fdr, bonferroni]

# Dependency graph
requires: []
provides:
  - "variantcentrifuge.association package importable with AssociationEngine, AssociationTest, TestResult, AssociationConfig"
  - "FisherExactTest: bit-identical reimplementation of gene_burden.py Fisher pipeline"
  - "correction.py: standalone apply_correction(pvals, method) for FDR and Bonferroni"
  - "AssociationEngine.from_names() with test registry; raises ValueError for unknown tests"
  - "Wide-format DataFrame output: gene, n_cases, n_controls, n_variants, {test}_p_value, {test}_corrected_p_value, {test}_or, {test}_or_ci_lower, {test}_or_ci_upper"
affects:
  - 18-02 (parity tests compare FisherExactTest output vs gene_burden.py)
  - 18-03 (AssociationAnalysisStage uses AssociationEngine)
  - 19 (burden regression implements AssociationTest ABC)
  - 20 (SKAT R backend implements AssociationTest ABC)
  - 21 (pure Python SKAT implements AssociationTest ABC)
  - 22 (ACAT-O implements AssociationTest ABC; uses apply_correction)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "AssociationTest ABC: name property + run(gene, contingency_data, config) + check_dependencies()"
    - "Lazy test registry: _build_registry() defers imports to avoid circular imports"
    - "Wide-format output: {test_name}_{field} column naming for multi-test DataFrames"
    - "Graceful degradation: apply_correction returns raw pvals when statsmodels absent"
    - "Continuity correction: 0.5 added to zero cells before CI computation (Haldane-Anscombe)"

key-files:
  created:
    - variantcentrifuge/association/__init__.py
    - variantcentrifuge/association/base.py
    - variantcentrifuge/association/correction.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/association/tests/__init__.py
    - variantcentrifuge/association/tests/fisher.py
  modified: []

key-decisions:
  - "Clean reimplementation (not delegation): fisher.py does NOT import from gene_burden.py — correct coupling direction for future deprecation"
  - "Lazy test registry via _build_registry() function to avoid circular imports at module load"
  - "Gene sort before correction in run_all() matches gene_burden.py line 411 sort — required for bit-identical FDR ordering"
  - "p_value=None semantics: skipped genes (zero variants) are excluded from output, not coded as 1.0"
  - "check_dependencies() decorated with noqa: B027 — intentional non-abstract method with default no-op for extensibility"

patterns-established:
  - "AssociationTest ABC: all future tests implement name, run(), optionally check_dependencies()"
  - "TestResult dataclass: standard result container for all tests; corrected_p_value filled post-hoc by engine"
  - "AssociationConfig: single config object passed through engine to all tests; mirrors gene_burden.py cfg dict keys"

# Metrics
duration: 6min
completed: 2026-02-19
---

# Phase 18 Plan 01: Foundation — Core Abstractions and Fisher Refactor Summary

**AssociationTest ABC, TestResult/AssociationConfig dataclasses, AssociationEngine orchestrator,
correction.py (FDR+Bonferroni), and FisherExactTest — a clean reimplementation of gene_burden.py's
Fisher pipeline producing bit-identical p-values and ORs via shared scipy.stats.fisher_exact calls.**

## Performance

- **Duration:** 6 min 17 sec
- **Started:** 2026-02-19T07:29:28Z
- **Completed:** 2026-02-19T07:35:45Z
- **Tasks:** 2
- **Files created:** 6

## Accomplishments

- Created `variantcentrifuge/association/` package with 5 modules fully importable
- FisherExactTest produces valid p-values (smoke test BRCA1: p=0.0065, OR=3.35) via identical
  scipy.stats.fisher_exact + statsmodels.Table2x2.oddsratio_confint call chain as gene_burden.py
- AssociationEngine.from_names(["skat"], ...) raises ValueError listing available tests
- Genes with zero qualifying variants return TestResult with p_value=None and are excluded from output
- apply_correction([0.01, 0.05, 0.1], "fdr") returns [0.03, 0.075, 0.1] — correct Benjamini-Hochberg
- All code passes ruff lint and format; no new typecheck errors in association/ module

## Task Commits

Each task was committed atomically:

1. **Task 1: Create association/ package with base abstractions and correction module** - `0dba65a` (feat)
2. **Task 2: Implement FisherExactTest and AssociationEngine** - `b02a1b4` (feat)

## Files Created/Modified

- `variantcentrifuge/association/__init__.py` — Package exports (AssociationTest, TestResult, AssociationConfig, AssociationEngine, apply_correction)
- `variantcentrifuge/association/base.py` — AssociationTest ABC, TestResult dataclass (11 fields), AssociationConfig dataclass (5 fields)
- `variantcentrifuge/association/correction.py` — apply_correction() with FDR/Bonferroni via statsmodels; graceful degradation if absent
- `variantcentrifuge/association/engine.py` — AssociationEngine: from_names() factory with registry, run_all() producing wide-format DataFrame
- `variantcentrifuge/association/tests/__init__.py` — Lazy __getattr__ so fisher.py imported only when FisherExactTest accessed
- `variantcentrifuge/association/tests/fisher.py` — FisherExactTest: table construction, fisher_exact, CI computation (score->normal->logit fallback)

## Decisions Made

- **Clean reimplementation**: fisher.py does NOT import from gene_burden.py. Bit-identity is
  guaranteed by sharing the same scipy/statsmodels functions with the same parameters, not by code
  sharing. This keeps the coupling direction correct for future gene_burden.py deprecation.
- **Lazy test registry**: _build_registry() function defers FisherExactTest import to avoid
  circular imports at package load time.
- **Gene sort order**: run_all() sorts by GENE key before running tests, matching gene_burden.py
  line 411. This is critical for bit-identical FDR correction ordering.
- **p_value=None semantics**: Zero-variant genes return TestResult with p_value=None and are
  excluded from the output DataFrame. This is distinct from p_value=1.0 (test ran, no signal).

## Deviations from Plan

None — plan executed exactly as written.

## Issues Encountered

- Ruff lint found 4 issues in initial implementation: RUF022 (__all__ not sorted), E501 (line too
  long in abstract method), B027 (empty non-abstract method), B905 (zip without strict=). All
  fixed before committing — no separate fix commit needed.
- Mypy flagged type narrowing issue in engine.py raw_pvals list comprehension (list[float | None]
  vs list[float]). Resolved with explicit type annotation and type: ignore comment.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- Phase 18 Plan 02 (parity tests) can now be written: FisherExactTest.run() is implemented and
  will be compared against gene_burden.perform_gene_burden_analysis() for identical contingency data
- Phase 18 Plan 03 (AssociationAnalysisStage) can reference AssociationEngine directly
- The ABC interface (name, run, check_dependencies) is the stable contract for Phases 19-22

---
*Phase: 18-foundation-core-abstractions-and-fisher-refactor*
*Completed: 2026-02-19*
