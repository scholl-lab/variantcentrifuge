# Requirements: v0.14.0 Report UX Overhaul

## Milestone Requirements

### JS Stack Modernization

- [ ] **STACK-01**: Upgrade DataTables from v1.10.20 to v2 (jQuery-optional)
- [ ] **STACK-02**: Remove jQuery dependency, replace with vanilla JS
- [ ] **STACK-03**: Add Tippy.js for viewport-aware, accessible tooltips
- [ ] **STACK-04**: Replace Plotly with Chart.js for lighter chart bundle (~65KB vs ~3.5MB)
- [ ] **STACK-05**: Add loading skeleton/spinner during DataTable initialization

### Information Hierarchy

- [ ] **HIER-01**: Move summary dashboard (cards + charts) above the variant table
- [ ] **HIER-02**: Expand summary cards: impact breakdown, inheritance distribution, top genes, sample count, filter criteria
- [ ] **HIER-03**: Add report metadata footer (filter criteria, VCF source, reference genome, pipeline version, run date)

### Semantic Color Coding

- [ ] **COLOR-01**: IMPACT column as colored badges (HIGH=red, MODERATE=orange, LOW=amber, MODIFIER=gray)
- [ ] **COLOR-02**: ClinVar classification badges (Pathogenic=red, Likely Path=orange, VUS=yellow, Likely Benign=light green, Benign=green)
- [ ] **COLOR-03**: Inheritance pattern badges (de novo=red, compound het=purple, AD=blue, AR=green, X-linked=teal)
- [ ] **COLOR-04**: Color-coded summary cards matching severity palette

### Data Table Design

- [ ] **TABLE-01**: Replace hover-expand with Tippy.js tooltips (viewport-aware, keyboard+hover trigger, Escape dismiss)
- [ ] **TABLE-02**: Sticky first column (GENE) using DataTables FixedColumns
- [ ] **TABLE-03**: Expandable row detail with chevron (click to show all fields in grouped key-value layout)
- [ ] **TABLE-04**: Zebra striping (alternating row backgrounds)
- [ ] **TABLE-05**: Enhanced header row (darker background, white text, bottom border)
- [ ] **TABLE-06**: Intelligent column widths per type (fixed for CHROM/POS, grow for GENE, truncation for HGVS, right-align for numbers)
- [ ] **TABLE-07**: Content density toggle (Compact/Regular/Relaxed) with localStorage persistence
- [ ] **TABLE-08**: Middle truncation for variant IDs, end truncation for HGVS, monospace for genetic notation
- [ ] **TABLE-09**: Prominent record count ("Showing X of Y variants")
- [ ] **TABLE-10**: Default page length 25 (up from 10)

### Column-Level Filtering

- [ ] **FILTER-01**: Port cohort report's column-level filtering to individual report (numeric sliders, categorical dropdowns, text search)
- [ ] **FILTER-02**: Active filter indicator (removable filter chips/tags)
- [ ] **FILTER-03**: "Include missing values" checkboxes for each filter
- [ ] **FILTER-04**: Reset all filters button

### Data Visualization

- [ ] **VIZ-01**: Impact distribution chart with semantic colors (above table)
- [ ] **VIZ-02**: Variant type breakdown chart (SNV/indel/MNV)
- [ ] **VIZ-03**: Chromosome distribution chart (horizontal bar)
- [ ] **VIZ-04**: Allele frequency histogram (log-scale)

### Print & PDF

- [ ] **PRINT-01**: Print stylesheet (@media print) hiding interactive controls, optimizing table layout
- [ ] **PRINT-02**: PDF export button (browser-based via html2pdf.js or window.print)

### Accessibility

- [ ] **A11Y-01**: ARIA roles and labels on table (role="grid"), filter controls, summary cards
- [ ] **A11Y-02**: Keyboard-accessible tooltips (tabindex="0", focus trigger, Escape dismiss)
- [ ] **A11Y-03**: Chart text alternatives (data table fallback or aria-label summary)
- [ ] **A11Y-04**: Skip-to-content links for keyboard navigation
- [ ] **A11Y-05**: Replace emoji link icons with SVG icons + screen-reader text
- [ ] **A11Y-06**: Sufficient color contrast (WCAG AA minimum)

## Future Requirements

- Cohort report improvements (port individual report enhancements back)
- Dark mode support (CSS custom properties)
- Filter presets (pathogenic rare, de novo candidates, compound het)
- URL-based filter state for shareable views
- Variant bookmarking/flagging with localStorage
- Gene info tooltips (OMIM, constraint scores)
- Keyboard shortcuts (j/k navigation, Enter/Esc for details)
- Column resize handles
- Virtual scrolling for 10,000+ variant datasets
- Pedigree visualization

## Out of Scope

- Cohort HTML report changes — focus on individual report gap
- Side panel variant detail — using expandable rows instead (simpler)
- Server-side data loading — reports must remain self-contained HTML files
- Custom user-defined filter presets — deferred to future milestone
- Excel report UX changes — separate scope

## Traceability

| Requirement | Phase |
|-------------|-------|
| STACK-01 | — |
| STACK-02 | — |
| STACK-03 | — |
| STACK-04 | — |
| STACK-05 | — |
| HIER-01 | — |
| HIER-02 | — |
| HIER-03 | — |
| COLOR-01 | — |
| COLOR-02 | — |
| COLOR-03 | — |
| COLOR-04 | — |
| TABLE-01 | — |
| TABLE-02 | — |
| TABLE-03 | — |
| TABLE-04 | — |
| TABLE-05 | — |
| TABLE-06 | — |
| TABLE-07 | — |
| TABLE-08 | — |
| TABLE-09 | — |
| TABLE-10 | — |
| FILTER-01 | — |
| FILTER-02 | — |
| FILTER-03 | — |
| FILTER-04 | — |
| VIZ-01 | — |
| VIZ-02 | — |
| VIZ-03 | — |
| VIZ-04 | — |
| PRINT-01 | — |
| PRINT-02 | — |
| A11Y-01 | — |
| A11Y-02 | — |
| A11Y-03 | — |
| A11Y-04 | — |
| A11Y-05 | — |
| A11Y-06 | — |
