# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.14.0 Report UX Overhaul — Phase 17 complete, verified ✓

## Current Position

Phase: 16 of 6 (Column-Level Filtering and Visualization)
Plan: 1 of 2 completed
Status: In progress
Last activity: 2026-02-18 — Completed 16-01-PLAN.md (noUiSlider assets and HTML/CSS scaffolding)

Progress: [███████████████████░] ~92% (Phase 16 Plan 1/2 in progress)

## Milestone Overview

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 13 | JS Stack Modernization | 5 | Complete ✓ (verified) |
| 14 | Information Hierarchy and Semantic Color Coding | 7 | Complete ✓ (verified) |
| 15 | Table Redesign | 10 | Complete ✓ (10/10 verified) |
| 16 | Column-Level Filtering and Visualization | 8 | In progress (Plan 1/2 complete) |
| 17 | Accessibility and Print/PDF | 8 | Complete ✓ (8/8 verified) |

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
| 14-01-01 | Dashboard positioned above variant table | 14-01 | Overview-first UX pattern shows summary before details |
| 14-01-02 | Two-row dashboard grid: metric cards (row 1) + dual charts (row 2) | 14-01 | Separates quick metrics from visualizations |
| 14-01-03 | Chart height reduced to 220px (from 300px) | 14-01 | Dashboard must fit above-the-fold on 1080p (~600px vertical budget) |
| 14-01-04 | Inheritance chart placeholder when no pedigree | 14-01 | Better UX than empty chart or error |
| 14-01-05 | Sample count extraction from GT column | 14-01 | Useful metric for multi-sample VCF reports |
| 14-01-06 | Filter "none" and "reference" from inheritance_distribution | 14-01 | These aren't real patterns, would clutter visualization |
| 14-02-01 | Badge colors applied via inline styles (not CSS classes) | 14-02 | Simplifies implementation and improves performance |
| 14-02-02 | DataTables render functions check type === 'display' | 14-02 | Preserves sorting/filtering on raw data values |
| 14-02-03 | ClinVar missing data shows empty cell, Inheritance shows 'Unknown' badge | 14-02 | ClinVar is optional metadata, Inheritance is core analysis field |
| 14-02-04 | LOW impact color updated to #f59e0b (WCAG-compliant amber) | 14-02 | Accessibility improvement, passes WCAG AA contrast on white |
| 14-02-05 | Metadata footer uses flexbox with gap and wrap | 14-02 | Modern responsive layout that adapts to narrow viewports |
| 14-02-06 | Filter expression in monospace code tag | 14-02 | Technical syntax should use monospace font for readability |
| 14-03-01 | Class-scoped fixtures for HTML rendering | 14-03 | Reuse expensive template rendering across tests for performance |
| 14-03-02 | String position comparison for layout testing | 14-03 | Verify ordering with html.find() positions instead of DOM parsing |
| 14-03-03 | Pattern-based badge testing | 14-03 | Test render function code patterns rather than runtime DOM output |
| 15-01-01 | Dark header gradient from #2c3e50 to #1a252f with white text | 15-01 | Charcoal gradient provides strong visual hierarchy while aligning with Phase 14 color palette |
| 15-01-02 | Zebra striping with #f8f9fa (odd) and white (even) | 15-01 | Subtle contrast improves row scannability without being visually distracting |
| 15-01-03 | Tooltip max-width reduced from 400px to 300px | 15-01 | TABLE-01 requirement for narrower tooltips, better for longer HGVS notation |
| 15-01-04 | Truncated-cell max-width increased from 150px to 180px | 15-01 | Better fit for HGVS notation while maintaining table compactness |
| 15-01-05 | FixedColumns JS loaded after datatables.min.js, before buttons | 15-01 | Extension dependencies: FixedColumns requires DataTables core, should load before optional buttons |
| 15-02-01 | FixedColumns left: 2 freezes both chevron and GENE columns | 15-02 | Control column (chevron) is index 0, GENE is index 1 - freezing 2 columns keeps both visible during horizontal scroll |
| 15-02-02 | formatChildRow auto-categorizes fields into Identifiers/Annotations/Scores/Links | 15-02 | Reduces cognitive load by grouping related fields rather than flat field list |
| 15-02-03 | Density defaults to compact mode with 25 rows per page | 15-02 | Maximizes data density for analysis-focused users, easily toggled to regular/relaxed |
| 15-02-04 | Column intelligence based on field name patterns | 15-02 | CHROM/POS fixed width, GENE grows, HGVS monospace, scores right-aligned for type-appropriate rendering |
| 15-02-05 | Enhanced tooltip triggers (mouseenter/focus/touch/immediate) | 15-02 | Enables keyboard navigation, mobile interaction, and reduces latency for accessibility |
| 17-01-01 | Use role="table" not role="grid" for DataTables | 17-01 | ARIA grid implies editable cells; DataTables is read-only so role="table" is semantically correct |
| 17-01-02 | SVG icons with aria-hidden plus sr-only text | 17-01 | Dual accessibility: SVG hidden from screen readers, adjacent span provides descriptive text |
| 17-01-03 | All badge colors meet WCAG AA 4.5:1 contrast ratio | 17-01 | MODERATE #c05000 (4.6:1), LOW #b45309 (5.2:1), plus ClinVar and Inheritance badge fixes |
| 17-01-04 | Truncated cells keyboard-accessible via tabindex="0" | 17-01 | Enables keyboard focus to trigger Tippy.js tooltips configured with focus trigger |
| 17-02-01 | Canvas elements hidden from screen readers (aria-hidden="true") with data table fallbacks | 17-02 | Charts inaccessible to screen readers; data tables provide equivalent text alternative |
| 17-02-02 | Data tables use sr-only class for visual hiding but screen reader accessibility | 17-02 | Prevents duplication for sighted users while providing data access to AT |
| 17-02-03 | Print mode shows chart data tables instead of canvas (canvas display:none) | 17-02 | Canvas doesn't print well; data tables provide printable chart information |
| 17-02-04 | FixedColumns collapse to static positioning in print to prevent duplicate columns | 17-02 | DataTables FixedColumns creates clones that would print as duplicates |
| 17-02-05 | Detail panels hidden by default in print mode | 17-02 | Users can expand specific rows before printing; default hide prevents clutter |
| 17-02-06 | Download PDF button uses window.print() for browser-native print dialog | 17-02 | Leverages browser PDF export capabilities without server-side rendering |
| 17-03-01 | Two test class organization: TestPhase17Accessibility and TestPhase17ChartAndPrint | 17-03 | Logical separation by feature domain (15 + 19 tests = 34 total) |
| 17-03-02 | WCAG contrast helper function with luminance calculation | 17-03 | Automated validation prevents future color regressions, reusable for other components |
| 17-03-03 | Print stylesheet section testing extracts @media print block | 17-03 | Validates comprehensive print behavior (hide controls, collapse columns, prevent breaks) |
| 17-03-04 | Two-tiered testing: template source patterns + rendered HTML checks | 17-03 | Comprehensive coverage of both static patterns and dynamic Jinja rendering |
| 16-01-01 | noUiSlider JS loaded after jQuery slim and before DataTables | 16-01 | No deps on other libs, should load before DataTables custom code |
| 16-01-02 | chevron-icon CSS class for collapsible section (not chevron used for row expansion) | 16-01 | Avoids CSS conflict with Phase 15 table row expansion toggle |
| 16-01-03 | Filter row starts hidden with empty mount point; Plan 02 populates dynamically | 16-01 | Separates structure (Plan 01) from behavior (Plan 02) for context budget |
| 16-01-04 | Visualization section collapsed by default (aria-expanded=false, hidden attribute) | 16-01 | Keeps above-the-fold clean; user clicks to expand |

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-18
Stopped at: Completed 16-01-PLAN.md (noUiSlider assets vendored, HTML/CSS scaffolding added)
Resume file: None
Next: Execute 16-02-PLAN.md (JavaScript behavior for filter controls, chip strip, and visualization charts)
