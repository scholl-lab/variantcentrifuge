# Phase 16: Column-Level Filtering and Visualization - Context

**Gathered:** 2026-02-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Add per-column filtering controls (numeric sliders, categorical dropdowns, text search) and expanded visualizations (variant type breakdown, chromosome distribution, allele frequency histogram) to the individual HTML report. Filters and charts interact — charts update when filters change. The cohort report's filter system is the reference pattern.

</domain>

<decisions>
## Implementation Decisions

### Filter control placement
- Header row filters: each filterable column gets an inline control directly in/below the column header
- Toggle button to show/hide the filter row (not always visible)
- Keep DataTables global search alongside column filters (AND logic between them)
- Curated subset of filterable columns only (IMPACT, ClinVar, Inheritance, scores/AF, GENE — not every visible column)

### Filter chip behavior
- Active filter chips displayed in a dedicated strip above the table, below the toolbar
- Strip shows "N filters active" count alongside individual removable chips
- "Reset All" button only (no keyboard shortcut — Escape is already used for tooltips)
- No filter persistence — filters reset on page reload (clean slate each time)

### Chart layout and interaction
- New charts (variant type breakdown, chromosome distribution, allele frequency histogram) go in a separate collapsible "Visualizations" section between the Phase 14 dashboard and the table
- Collapsed by default — users click to expand
- No animation on filter change — snap to new values instantly
- Phase 14 dashboard charts (impact distribution, inheritance) also update when column filters change — all charts reflect the filtered view

### Missing values handling
- Global "Include missing values" toggle (not per-column) defaults to checked (include missing)
- Toggle is its own UI element, not represented as a chip
- Missing/empty cells left blank — no dash or N/A marker

### Claude's Discretion
- Specific filter control widgets per column type (exact slider vs input range, dropdown style)
- Chart sizing and arrangement within the collapsible section
- Filter row toggle button placement and icon
- How to wire DataTables column search API with custom filter controls
- Allele frequency histogram bin sizing and log-scale axis implementation

</decisions>

<specifics>
## Specific Ideas

- Cohort report's filter system is the reference pattern for column-level filtering
- Charts should use the semantic color palette established in Phase 14 (HIGH=red, MODERATE=orange, LOW=amber, MODIFIER=gray)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 16-column-level-filtering-and-visualization*
*Context gathered: 2026-02-17*
