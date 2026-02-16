# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.14.0 Report UX Overhaul — Phase 13 complete, Phase 14 next

## Current Position

Phase: 13 of 5 (JS Stack Modernization)
Plan: 3 of 3 completed
Status: Phase complete, verified ✓
Last activity: 2026-02-16 — Phase 13 complete (3/3 plans, 5/5 must-haves verified)

Progress: [████░░░░░░░░░░░░░░░░] 20% (1/5 phases)

## Milestone Overview

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 13 | JS Stack Modernization | 5 | Complete ✓ (verified) |
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
| 13-02-01 | Fixed asset key collision by namespacing with subdirectory prefix | 13-02 | Prevents CSS overwriting JS for same-stem filenames |
| 13-02-02 | Semantic color scheme for impact levels (HIGH=red, MODERATE=orange, LOW=amber, MODIFIER=gray) | 13-02 | Establishes consistent color language for Phase 14 |
| 13-02-03 | Conditional Tippy.js tooltip attachment | 13-02 | Only tooltips on actually truncated content, improves performance |
| 13-02-04 | Loading skeleton with shimmer during DataTable init | 13-02 | Better perceived performance, modern UX pattern |
| 13-03-01 | Class-scoped pytest fixtures for rendered HTML | 13-03 | Efficient test execution by reusing expensive rendering |
| 13-03-02 | Tests tolerate parallel execution (13-02 and 13-03 in parallel) | 13-03 | Graceful skips if template not yet modernized |
| 13-03-03 | Dedicated test for asset key namespacing | 13-03 | Prevents regression of 13-02 blocking issue |
| 13-FIX-01 | Added Popper.js as separate asset (tippy-bundle UMD requires window.Popper) | Post-verify | Tippy.js UMD won't initialize without Popper loaded first |
| 13-FIX-02 | Table must be visible during DataTables init (scrollX header sizing) | Post-verify | Hidden tables produce 0-height headers; skeleton hides on initComplete |
| 13-FIX-03 | Hide dt-scroll-body thead via CSS (DataTables scrollX duplicate) | Post-verify | Prevents ghost header row with duplicate sort arrows |
| 13-FIX-04 | Force left-align on all columns including dt-type-numeric | Post-verify | DataTables auto-right-aligns numeric columns, causing sort icon overlap |
| 13-FIX-05 | Franklin, SpliceAI, autopvs1 hidden by default; link order: ClinVar, gnomAD_2, Varsome first | Post-verify | Reduce column clutter; most-used links visible by default |

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-16
Stopped at: Phase 13 complete, verified, post-verification UI fixes applied
Resume file: None
Next: Plan Phase 14 (Information Hierarchy and Semantic Color Coding)
