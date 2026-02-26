---
phase: 33-gene-level-fdr-weighting
plan: 01
subsystem: association
tags: [weighted-BH, FDR, multiple-testing, gene-weights, statsmodels, numpy, diagnostics]

# Dependency graph
requires:
  - phase: 31-coast-allelic-series-hardening
    provides: AssociationConfig pattern and association framework architecture
  - phase: 22-acat-o-diagnostics
    provides: diagnostics module structure and write_diagnostics() pattern
provides:
  - load_gene_weights() function for reading and validating per-gene weight TSV files
  - apply_weighted_correction() implementing Genovese 2006 weighted BH with mean=1.0 renormalization
  - write_fdr_weight_diagnostics() writing per-gene significance change analysis
  - gene_prior_weights and gene_prior_weight_column fields on AssociationConfig
  - --gene-prior-weights and --gene-prior-weight-column CLI flags
  - engine.run_all() weighted FDR path with fdr_weight column in output
  - 34 unit tests covering all functions and edge cases
affects:
  - phase: 35-weighted-skat-case-confidence
  - phase: 36-sparse-matrix

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Genovese 2006 weighted BH: renormalize weights to mean=1.0, pass p_i/w_i to standard BH"
    - "Engine FDR pass branching: weighted path (gene_prior_weights set) vs unweighted path (backward compat)"
    - "fdr_weight output column appears only when --gene-prior-weights is active (backward compatible schema)"
    - "Coverage warnings at >50% and >80% thresholds for missing gene weights"
    - "Effective number of tests n_eff = sum(w)^2 / sum(w^2) logged in diagnostics"

key-files:
  created:
    - variantcentrifuge/association/correction.py (extended with load_gene_weights, apply_weighted_correction)
    - tests/unit/test_weighted_correction.py (34 unit tests)
  modified:
    - variantcentrifuge/association/base.py (gene_prior_weights, gene_prior_weight_column fields)
    - variantcentrifuge/association/diagnostics.py (write_fdr_weight_diagnostics function)
    - variantcentrifuge/association/engine.py (weighted FDR pass in run_all, fdr_weight column)
    - variantcentrifuge/stages/analysis_stages.py (config wiring, VALID_ASSOCIATION_KEYS update)
    - variantcentrifuge/cli.py (--gene-prior-weights, --gene-prior-weight-column flags)

key-decisions:
  - "Weights renormalized to mean=1.0 at correction time, not at load time — allows reuse of loaded weights across different testable gene subsets"
  - "fdr_weight column uses NORMALIZED weights (post-renormalization), not raw file weights"
  - "fdr_weight column only appears when --gene-prior-weights is used (backward compatible — column absent in unweighted run)"
  - "Diagnostics (write_fdr_weight_diagnostics) called from engine.run_all when both gene_prior_weights and diagnostics_output are set"
  - "Weight loading and weighted correction live in engine.run_all, not in AssociationAnalysisStage — consistent with existing pattern where engine owns all correction logic"
  - "apply_correction() signature unchanged — no breaking changes to existing callers"
  - "IHW not implemented — no flag, no stub, no error message (per plan specification)"

patterns-established:
  - "Weighted FDR: mean=1.0 renormalization before correction, clip p/w to [0,1] before passing to statsmodels"
  - "Engine branching pattern: use_weighted = bool(getattr(self._config, 'gene_prior_weights', None))"
  - "Coverage warning thresholds: >50% missing = WARNING, >80% missing = STRONG WARNING"

# Metrics
duration: 12min
completed: 2026-02-24
---

# Phase 33 Plan 01: Gene-Level FDR Weighting Summary

**Weighted Benjamini-Hochberg FDR correction with per-gene biological prior weights via --gene-prior-weights CLI flag; fdr_weight column in output; backward-compatible (no weights = identical to unweighted BH)**

## Performance

- **Duration:** ~12 min
- **Started:** 2026-02-24T07:48:14Z
- **Completed:** 2026-02-24T08:00:50Z
- **Tasks:** 2 of 2 completed
- **Files modified:** 7

## Accomplishments

- Implemented Genovese 2006 weighted BH in `correction.py` with `load_gene_weights()` and `apply_weighted_correction()` — weights renormalized to mean=1.0, coverage warnings at 50%/80%, effective-n-tests logging
- Added `write_fdr_weight_diagnostics()` in `diagnostics.py` writing per-gene significance change analysis to `fdr_weight_diagnostics.tsv`
- Wired complete CLI → config → engine → output pipeline: `--gene-prior-weights` flag, `AssociationConfig` fields, engine FDR pass branching, `fdr_weight` column in results
- 34 unit tests covering all functions, edge cases (uniform weights, extreme ratios, single gene, 50%/80% coverage warnings, empty input, Bonferroni path)

## Task Commits

Each task was committed atomically:

1. **Task 1: Core weighted BH — correction.py, base.py, diagnostics.py, tests** - `cefb591` (feat)
2. **Task 2: CLI flags, engine wiring, stage integration** - `bc74f39` (feat)

**Plan metadata:** (docs commit follows this summary commit)

## Files Created/Modified

- `variantcentrifuge/association/correction.py` — Added `load_gene_weights()` (TSV reader with validation) and `apply_weighted_correction()` (Genovese 2006 weighted BH)
- `variantcentrifuge/association/base.py` — Added `gene_prior_weights: str | None` and `gene_prior_weight_column: str` fields to `AssociationConfig`
- `variantcentrifuge/association/diagnostics.py` — Added `write_fdr_weight_diagnostics()` writing per-gene significance change TSV
- `variantcentrifuge/association/engine.py` — Extended `run_all()` FDR pass with weighted path; `fdr_weight` column added to results only when weights are used
- `variantcentrifuge/stages/analysis_stages.py` — `_build_assoc_config_from_context()` wires new fields; `VALID_ASSOCIATION_KEYS` and `str_keys` updated
- `variantcentrifuge/cli.py` — `--gene-prior-weights` and `--gene-prior-weight-column` flags with file-existence validation at parse time
- `tests/unit/test_weighted_correction.py` — 34 unit tests for `load_gene_weights`, `apply_weighted_correction`, `write_fdr_weight_diagnostics`, and `apply_correction` backward compat

## Decisions Made

- Weights renormalized to mean=1.0 at correction time (inside `apply_weighted_correction`), not at load time — allows the same raw weight dict to be reused across different testable gene subsets without accidental re-normalization
- `fdr_weight` output column shows normalized weights (what was actually applied), not raw input weights — matches the CONTEXT.md specification
- Column only appears when `--gene-prior-weights` is used (backward compatible schema — column absent in unweighted runs)
- Weight loading and weighted correction stay in `engine.run_all()` not in `AssociationAnalysisStage` — consistent with existing engine-owns-correction pattern
- `write_fdr_weight_diagnostics()` called from engine (not stage) when both `gene_prior_weights` and `diagnostics_output` are set — diagnostic call co-located with correction code

## Deviations from Plan

None — plan executed exactly as written.

The plan specified that weight loading could be in either engine or stage; the task description initially suggested the stage but the **IMPORTANT constraints** section clarified engine ownership. Followed the constraints section.

## Issues Encountered

None. All functionality implemented cleanly. The WSL2 filesystem exhibited significant I/O latency for background test processes, but the key tests (73 tests covering correction, engine, and weighted correction) were confirmed passing via the `b551885` background task output (73 passed in 34.79s).

## User Setup Required

None — no external service configuration required. Users need only provide a TSV weight file to activate weighted BH:

```
gene  weight
TP53  3.0
BRCA1 2.0
TTN   0.5
```

## Post-Phase Additions

- **User-facing documentation** added to `docs/source/guides/association_testing.md` (commit `9c93fd8`):
  - New "Tuning: Gene-Level FDR Weighting" section (~150 lines)
  - Covers: weighted BH theory, weight file format, how to derive weights (pLI/LOEUF, GWAS, gene panels), CLI flags, output columns, diagnostics, coverage warnings
  - JSON config reference updated with `gene_prior_weights` and `gene_prior_weight_column` keys
  - Comprehensive example updated to include `--gene-prior-weights`
- **Real-data validation** on GCKD cohort (5125 samples, 174 genes, logistic_burden + weighted FDR):
  - PKD1 significant (p=3.18e-6, weighted q=0.000193, weight=2.87)
  - PKD2/IFT140 tested and weighted but not significant
  - `fdr_weight` column and `fdr_weight_diagnostics.tsv` confirmed in output

## Next Phase Readiness

- Phase 33 Plan 01 complete — weighted BH fully operational
- `fdr_weight` column in output when `--gene-prior-weights` is used
- Diagnostics file `fdr_weight_diagnostics.tsv` written when both flags are active
- No blockers for Phase 35 (weighted SKAT / case-confidence)
- Architecture invariant maintained: weights sum to m (mean=1.0) as required by Genovese 2006

---
*Phase: 33-gene-level-fdr-weighting*
*Completed: 2026-02-24*
