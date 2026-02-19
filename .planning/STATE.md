# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-19)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.15.0 Modular Rare Variant Association Framework

## Current Position

Phase: Not started (defining requirements)
Plan: —
Status: Defining requirements
Last activity: 2026-02-19 — Milestone v0.15.0 started

Progress: ░░░░░░░░░░░░░░░░░░░░░ 0%

## Milestone Overview

(Phases not yet defined — roadmap creation in progress)

## Accumulated Context

### Decisions

| ID | Decision | Plan | Rationale |
|----|----------|------|-----------|
| MILE-01 | Dual SKAT backend (R via rpy2 + pure Python) | Milestone | R is gold standard for validation; Python for portability |
| MILE-02 | Compiled Davies method via ctypes (qfc.c) | Milestone | Exact p-values without R dependency; Liu fallback for safety |
| MILE-03 | Support both AKT and PLINK for PCA | Milestone | AKT for VCF-native workflows, PLINK for pre-existing BED files |
| MILE-04 | Full scope (Steps 1-7 from design doc) | Milestone | Include allelic series, functional weights, quantitative traits, JSON config |

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-19
Stopped at: Milestone v0.15.0 initialization — defining requirements
Resume file: None
Next: Define requirements and create roadmap
