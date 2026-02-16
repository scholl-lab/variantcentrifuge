# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.14.0 Report UX Overhaul — Phase 13 (JS Stack Modernization)

## Current Position

Phase: 13 of 5 (JS Stack Modernization)
Plan: 01 of 2 completed
Status: In progress
Last activity: 2026-02-16 — Completed 13-01-PLAN.md

Progress: [██░░░░░░░░░░░░░░░░░░] 10% (1/10 total plans across all phases)

## Milestone Overview

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 13 | JS Stack Modernization | 5 | In Progress (1/2 plans) |
| 14 | Information Hierarchy and Semantic Color Coding | 7 | Not Started |
| 15 | Table Redesign | 10 | Not Started |
| 16 | Column-Level Filtering and Visualization | 8 | Not Started |
| 17 | Accessibility and Print/PDF | 8 | Not Started |

## Accumulated Context

### Decisions

| ID | Decision | Plan | Rationale |
|----|----------|------|-----------|
| MILE-01 | Individual HTML report only (cohort report out of scope for v0.14.0) | Milestone | Focus on single-report UX first |
| MILE-02 | Expandable rows only (no side panel) | Milestone | Simpler implementation, better mobile support |
| MILE-03 | Modernize JS stack: DataTables v2, drop jQuery, add Tippy.js | Milestone | Modern libraries, better performance |
| MILE-04 | Replace Plotly with Chart.js (~65KB vs ~3.5MB) | Milestone | 98% size reduction, faster load times |
| MILE-05 | Phase 13 is foundation -- everything depends on the stack modernization | Milestone | Dependency planning |
| MILE-06 | COLOR and HIER grouped into Phase 14 (tightly coupled) | Milestone | Colored summary cards + colored badges belong together |
| MILE-07 | FILTER and VIZ grouped into Phase 16 (both above-table) | Milestone | Interactive features that update together |
| MILE-08 | A11Y and PRINT grouped into Phase 17 (polish phase) | Milestone | Needs all components to exist first |
| 13-01-01 | Vendor all JS/CSS libraries rather than using CDN | 13-01 | Enables offline HTML reports, improves reproducibility |
| 13-01-02 | Chart.js 4.4.8 replaces Plotly (65KB vs 3.5MB) | 13-01 | 98% size reduction, sufficient functionality |
| 13-01-03 | Use jQuery slim build (no AJAX/effects) | 13-01 | DataTables v2 requires jQuery, slim build saves 30KB |
| 13-01-04 | Inline asset embedding pattern | 13-01 | Enables true single-file HTML reports |

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-16 16:05 UTC
Stopped at: Completed 13-01-PLAN.md (Vendor JS/CSS Libraries)
Resume file: None
Next: Execute 13-02-PLAN.md (Template Rewrite)
