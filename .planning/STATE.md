# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.14.0 Report UX Overhaul — Phase 13 (JS Stack Modernization)

## Current Position

Phase: 13 — JS Stack Modernization
Plan: Not started
Status: Roadmap created, awaiting phase planning
Last activity: 2026-02-16 — Roadmap created for v0.14.0

Progress: [░░░░░░░░░░░░░░░░░░░░] 0% (0/5 phases)

## Milestone Overview

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 13 | JS Stack Modernization | 5 | Not Started |
| 14 | Information Hierarchy and Semantic Color Coding | 7 | Not Started |
| 15 | Table Redesign | 10 | Not Started |
| 16 | Column-Level Filtering and Visualization | 8 | Not Started |
| 17 | Accessibility and Print/PDF | 8 | Not Started |

## Accumulated Context

### Decisions

- Individual HTML report only (cohort report out of scope for v0.14.0)
- Expandable rows only (no side panel)
- Modernize JS stack: DataTables v2, drop jQuery, add Tippy.js
- Replace Plotly with Chart.js (~65KB vs ~3.5MB)
- Phase 13 is foundation -- everything depends on the stack modernization
- COLOR and HIER grouped into Phase 14 (tightly coupled: colored summary cards + colored badges)
- FILTER and VIZ grouped into Phase 16 (both are above-table interactive features, VIZ charts update with filters)
- A11Y and PRINT grouped into Phase 17 (polish phase: needs all components to exist first)

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-16
Stopped at: Roadmap created for v0.14.0
Resume file: None
Next: Plan Phase 13 (JS Stack Modernization)
