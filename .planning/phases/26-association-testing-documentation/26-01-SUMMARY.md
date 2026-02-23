---
phase: 26-association-testing-documentation
plan: 01
subsystem: documentation
tags: [association-testing, skat, coast, fisher, burden-test, acat-o, sphinx, myst-markdown]

# Dependency graph
requires:
  - phase: 25-python-default-backends
    provides: Python default backends, R deprecation, ACAT-V — all documented in this guide
  - phase: 22-acat-diagnostics
    provides: ACAT-O omnibus, diagnostics (lambda_GC, QQ) described in guide
  - phase: 23-pca-weights-coast-config
    provides: PCA integration, variant weights, COAST test, JSON config described in guide
  - phase: 24-pure-python-coast
    provides: Pure Python COAST backend documented
provides:
  - Comprehensive association testing user guide (docs/source/guides/association_testing.md, 843 lines)
  - Updated guides index toctree (docs/source/guides/index.md) linking to the new guide
affects:
  - 26-02 (remaining documentation tasks in this phase)
  - 27-association-performance-optimizations (guide accuracy must be maintained as performance changes are made)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "MyST Markdown for Sphinx documentation with admonitions (:::{ note/warning/tip })"
    - "Command + output snippet style: CLI invocation with representative output lines"

key-files:
  created:
    - docs/source/guides/association_testing.md
  modified:
    - docs/source/guides/index.md

key-decisions:
  - "Guide placed after cohort_analysis in toctree (association testing extends cohort workflows)"
  - "All 13 sections written from source code inspection — no guessed CLI flags or output columns"
  - "Test comparison table covers all 6 test types with when-to-use guidance"
  - "Troubleshooting section covers 8 common failure modes from RESEARCH.md pitfalls"

patterns-established:
  - "Source-verified accuracy: all CLI flags from cli.py, output columns from engine.py, config keys from analysis_stages.py"
  - "Dual-audience writing: bioinformatician and clinical geneticist level simultaneously"

# Metrics
duration: 4min
completed: 2026-02-22
---

# Phase 26 Plan 01: Association Testing Guide Summary

**843-line MyST Markdown guide covering Fisher, logistic/linear burden, SKAT-O, COAST, and ACAT-O with covariates, PCA, weights, diagnostics, JSON config, test comparison table, and 8-section troubleshooting**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-22T19:58:12Z
- **Completed:** 2026-02-22T20:02:06Z
- **Tasks:** 2/2
- **Files modified:** 2

## Accomplishments

- Created `docs/source/guides/association_testing.md` (843 lines, all 13 sections from CONTEXT.md)
- All CLI flags verified against `cli.py` (lines 408-530): `--perform-association`, `--association-tests`, `--skat-backend`, `--coast-backend`, `--covariate-file`, `--covariates`, `--categorical-covariates`, `--trait-type`, `--variant-weights`, `--variant-weight-params`, `--coast-weights`, `--pca-file`, `--pca-tool`, `--pca-components`, `--diagnostics-output`
- All output columns verified against `engine.py`: `fisher_p_value`, `fisher_or`, `logistic_burden_beta`, `skat_p_value`, `skat_o_rho`, `acat_v_p`, `coast_p_value`, `coast_n_bmv/dmv/ptv`, `acat_o_corrected_p_value`
- All 28 JSON config keys documented with defaults, types, and descriptions from `VALID_ASSOCIATION_KEYS`
- Added `association_testing` to `docs/source/guides/index.md` toctree after `cohort_analysis`

## Task Commits

1. **Task 1: Write association testing guide** - `957f149` (docs)
2. **Task 2: Add guide to guides index toctree** - `ee9c4cc` (docs)

## Files Created/Modified

- `docs/source/guides/association_testing.md` — 843-line comprehensive user guide
- `docs/source/guides/index.md` — added `association_testing` toctree entry after `cohort_analysis`

## Decisions Made

- Placed guide after `cohort_analysis` in the guides toctree (association testing is a natural extension of cohort analysis workflows)
- All content written by reading source code directly as recommended by RESEARCH.md — no content guessed from training data
- Used MyST admonitions (`:::{note}`, `:::{warning}`, `:::{tip}`) for important callouts per CONTEXT.md discretion
- Used "Command + output snippet" style throughout for practical, copy-paste-ready examples

## Deviations from Plan

None — plan executed exactly as written. All 13 sections present, 843 lines (>400 minimum), CLI flags accurate, guides index updated.

## Issues Encountered

None.

## User Setup Required

None — documentation-only plan, no external service configuration required.

## Next Phase Readiness

- Plan 26-01 deliverables complete: association testing guide written and linked
- Plan 26-02 can proceed: existing docs updates (usage.md, cohort_analysis.md, faq.md, README.md, API stubs, changelog)
- No blockers

---
*Phase: 26-association-testing-documentation*
*Completed: 2026-02-22*
