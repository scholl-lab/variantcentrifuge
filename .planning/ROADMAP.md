# Roadmap: VariantCentrifuge

## Milestones

- SHIPPED **v0.12.1 Baseline** — Phases 1-5 (shipped 2026-02-14, pre-GSD)
- SHIPPED **v0.13.0 Performance Optimization** — Phases 6-12 (shipped 2026-02-16) — [archive](milestones/v0.13.0-ROADMAP.md)
- ACTIVE **v0.14.0 Report UX Overhaul** — Phases 13-17

## Phases

<details>
<summary>v0.12.1 Baseline (Phases 1-5) — SHIPPED 2026-02-14</summary>

**Milestone Goal:** Full-featured variant analysis pipeline with two architectures, 40+ stages, inheritance analysis, gene burden, scoring, multi-format output. All 30 historical issues resolved. 1035 tests passing. CI/CD with Docker, docs, and multi-platform testing.

**Note:** Phases 1-5 represent completed work before GSD tracking. Details available in git history.

</details>

<details>
<summary>v0.13.0 Performance Optimization (Phases 6-12) — SHIPPED 2026-02-16</summary>

- [x] Phase 6: Benchmark Framework (4/4 plans) — completed 2026-02-14
- [x] Phase 7: Quick Wins - Tier 1 (3/3 plans) — completed 2026-02-14
- [x] Phase 8: DataFrame Optimization (4/4 plans) — completed 2026-02-14
- [x] Phase 9: Inheritance Analysis Optimization (5/5 plans) — completed 2026-02-14
- [x] Phase 10: Output Optimization (3/3 plans) — completed 2026-02-15
- [x] Phase 11: Pipeline I/O Elimination (3/3 plans) — completed 2026-02-15
- [x] Phase 12: Parallelization & Chunking (4/4 plans) — completed 2026-02-16

Full details: [milestones/v0.13.0-ROADMAP.md](milestones/v0.13.0-ROADMAP.md)

</details>

### v0.14.0 Report UX Overhaul (Phases 13-17)

Transform the individual HTML report from a functional-but-dated developer tool (scored 4.8/10) into a modern, accessible clinical genomics report with semantic color coding, enterprise-grade table design, column-level filtering, expanded visualizations, and print support. All changes target the individual HTML report only. Implementation involves Python template code and embedded HTML/CSS/JS.

---

#### Phase 13: JS Stack Modernization

**Goal:** The individual report runs on a modern, lightweight JS stack with no jQuery dependency, enabling all subsequent UX work.

**Dependencies:** None (foundation phase)

**Requirements:** STACK-01, STACK-02, STACK-03, STACK-04, STACK-05

**Plans:** 3 plans

Plans:
- [x] 13-01-PLAN.md — Vendor JS/CSS assets and update report generator to load them
- [x] 13-02-PLAN.md — Rewrite HTML template for modern stack (DataTables v2, Chart.js, Tippy.js)
- [x] 13-03-PLAN.md — Add tests for asset loading and template rendering

**Success Criteria:**

1. The individual HTML report loads with DataTables v2 and all existing functionality (sorting, pagination, search, column visibility, horizontal scroll) works identically to the current report
2. No jQuery is loaded or referenced anywhere in the generated HTML -- all JS is vanilla or library-specific
3. Tippy.js is loaded and available for tooltip use (integration with table cells happens in Phase 15)
4. The impact distribution chart renders using Chart.js instead of Plotly, reducing JS bundle from ~3.5MB to ~65KB
5. A loading skeleton or spinner is visible during DataTable initialization, replaced by the table once ready

---

#### Phase 14: Information Hierarchy and Semantic Color Coding

**Goal:** Users see a clinically meaningful overview (summary dashboard with colored cards, charts) before the variant table, with semantic color coding applied throughout.

**Dependencies:** Phase 13 (Chart.js for colored charts, DataTables v2 for badge rendering)

**Requirements:** HIER-01, HIER-02, HIER-03, COLOR-01, COLOR-02, COLOR-03, COLOR-04

**Plans:** 3 plans

Plans:
- [x] 14-01-PLAN.md — Expand summary data and restructure template with dashboard above table
- [x] 14-02-PLAN.md — Add semantic color badges to IMPACT/ClinVar/Inheritance columns and metadata footer
- [x] 14-03-PLAN.md — Tests for dashboard, badges, and metadata footer

**Success Criteria:**

1. Opening the report shows summary cards and charts above the variant table without scrolling -- the table is below the fold on a standard 1080p display
2. Summary cards display impact breakdown, inheritance distribution, top genes, sample count, and filter criteria -- each card color-coded to match severity palette (red for HIGH, orange for MODERATE, etc.)
3. IMPACT column values in the table render as colored badges (HIGH=red, MODERATE=orange, LOW=amber, MODIFIER=gray) instead of plain text
4. ClinVar classification values render as colored badges (Pathogenic=red through Benign=green) and inheritance patterns render as colored badges (de novo=red, compound het=purple, AD=blue, AR=green, X-linked=teal)
5. A report metadata footer displays filter criteria, VCF source, reference genome, pipeline version, and run date

---

#### Phase 15: Table Redesign

**Goal:** Users can efficiently scan, explore, and read variant data in a modern enterprise-grade table with tooltips, sticky columns, expandable row details, and configurable density.

**Dependencies:** Phase 13 (DataTables v2, Tippy.js), Phase 14 (color badges in cells)

**Requirements:** TABLE-01, TABLE-02, TABLE-03, TABLE-04, TABLE-05, TABLE-06, TABLE-07, TABLE-08, TABLE-09, TABLE-10

**Plans:** 3 plans

Plans:
- [x] 15-01-PLAN.md — Vendor FixedColumns extension and add CSS styles for table redesign
- [x] 15-02-PLAN.md — Wire up JS behavior (column intelligence, expandable rows, density toggle, sticky columns)
- [x] 15-03-PLAN.md — Tests for Phase 15 table redesign features

**Success Criteria:**

1. Hovering over or focusing (keyboard) a truncated cell shows a Tippy.js tooltip with the full content, positioned viewport-aware (no clipping at edges), dismissable with Escape -- the old yellow hover-expand is completely removed
2. The GENE column stays visible (sticky) when scrolling horizontally, and the header row has a dark background with white text and a colored bottom border
3. Clicking a chevron on any row expands an inline detail panel showing all variant fields grouped by category (Identifiers, Annotations, Scores, Links) in a key-value layout
4. Users can switch between Compact, Regular, and Relaxed density modes via a toolbar control, with their preference persisted across sessions (localStorage)
5. The table shows "Showing X of Y variants" prominently, defaults to 25 rows per page, uses zebra striping, applies intelligent column widths (fixed for CHROM/POS, grow for GENE, truncation with monospace for HGVS, right-aligned numbers), and uses middle truncation for variant IDs

---

#### Phase 16: Column-Level Filtering and Visualization

**Goal:** Users can filter variants by individual columns (numeric ranges, categorical dropdowns, text search) and see expanded visualizations that update with filters.

**Dependencies:** Phase 13 (Chart.js), Phase 15 (stable table structure)

**Requirements:** FILTER-01, FILTER-02, FILTER-03, FILTER-04, VIZ-01, VIZ-02, VIZ-03, VIZ-04

**Success Criteria:**

1. Each filterable column has an appropriate control: numeric sliders for scores/frequencies, categorical dropdowns for IMPACT/ClinVar/inheritance, and text search for gene names -- matching the cohort report's filter system
2. Active filters are displayed as removable chips/tags above the table, and a "Reset all filters" button clears all filters in one click
3. Each filter includes an "Include missing values" checkbox so users can control whether variants with missing data are shown or hidden
4. Impact distribution, variant type breakdown, chromosome distribution, and allele frequency histogram charts are displayed above the table using semantic colors, and they visually update when filters change
5. The allele frequency histogram uses log-scale for its axis to properly visualize the distribution of rare variants

---

#### Phase 17: Accessibility and Print/PDF

**Goal:** The report meets WCAG 2.1 AA accessibility standards and can be printed or exported to PDF with a clean, usable layout.

**Dependencies:** Phases 14-16 (all visual components and interactive features must exist before accessibility audit and print optimization)

**Requirements:** A11Y-01, A11Y-02, A11Y-03, A11Y-04, A11Y-05, A11Y-06, PRINT-01, PRINT-02

**Success Criteria:**

1. The data table has `role="grid"`, all filter controls have `aria-label` attributes, summary cards have semantic meaning for screen readers, and skip-to-content links allow keyboard users to jump past the header
2. All tooltips (Tippy.js) are keyboard-accessible via focus trigger with `tabindex="0"`, dismissable with Escape, and work on touch devices
3. Every chart has a text alternative (data table fallback or aria-label summary) so screen reader users can access the visualized data
4. Emoji link icons are replaced with SVG icons that have `aria-hidden="true"` plus screen-reader-only text, and all color-coded badges maintain WCAG AA contrast ratios (4.5:1 minimum)
5. A print stylesheet hides interactive controls (filters, pagination, column toggles, density selector) and optimizes the table layout for paper, and a "Download PDF" button triggers browser-based PDF export

---

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1-5. Baseline | v0.12.1 | N/A | Complete | 2026-02-14 |
| 6. Benchmark Framework | v0.13.0 | 4/4 | Complete | 2026-02-14 |
| 7. Quick Wins - Tier 1 | v0.13.0 | 3/3 | Complete | 2026-02-14 |
| 8. DataFrame Optimization | v0.13.0 | 4/4 | Complete | 2026-02-14 |
| 9. Inheritance Analysis Optimization | v0.13.0 | 5/5 | Complete | 2026-02-14 |
| 10. Output Optimization | v0.13.0 | 3/3 | Complete | 2026-02-15 |
| 11. Pipeline I/O Elimination | v0.13.0 | 3/3 | Complete | 2026-02-15 |
| 12. Parallelization & Chunking | v0.13.0 | 4/4 | Complete | 2026-02-16 |
| 13. JS Stack Modernization | v0.14.0 | 3/3 | Complete | 2026-02-16 |
| 14. Information Hierarchy and Semantic Color Coding | v0.14.0 | 3/3 | Complete | 2026-02-16 |
| 15. Table Redesign | v0.14.0 | 3/3 | Complete | 2026-02-17 |
| 16. Column-Level Filtering and Visualization | v0.14.0 | 0/? | Not Started | — |
| 17. Accessibility and Print/PDF | v0.14.0 | 0/? | Not Started | — |
