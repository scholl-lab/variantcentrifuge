# Requirements: v0.14.0 Report UX Overhaul

## Milestone Requirements

### JS Stack Modernization

- [x] **STACK-01**: Upgrade DataTables from v1.10.20 to v2 (jQuery-optional)
- [x] **STACK-02**: Remove jQuery dependency, replace with vanilla JS
- [x] **STACK-03**: Add Tippy.js for viewport-aware, accessible tooltips
- [x] **STACK-04**: Replace Plotly with Chart.js for lighter chart bundle (~65KB vs ~3.5MB)
- [x] **STACK-05**: Add loading skeleton/spinner during DataTable initialization

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

- [x] **TABLE-01**: Replace hover-expand with Tippy.js tooltips (viewport-aware, keyboard+hover trigger, Escape dismiss)
- [x] **TABLE-02**: Sticky first column (GENE) using DataTables FixedColumns
- [x] **TABLE-03**: Expandable row detail with chevron (click to show all fields in grouped key-value layout)
- [x] **TABLE-04**: Zebra striping (alternating row backgrounds)
- [x] **TABLE-05**: Enhanced header row (darker background, white text, bottom border)
- [x] **TABLE-06**: Intelligent column widths per type (fixed for CHROM/POS, grow for GENE, truncation for HGVS, right-align for numbers)
- [x] **TABLE-07**: Content density toggle (Compact/Regular/Relaxed) with localStorage persistence
- [x] **TABLE-08**: Middle truncation for variant IDs, end truncation for HGVS, monospace for genetic notation
- [x] **TABLE-09**: Prominent record count ("Showing X of Y variants")
- [x] **TABLE-10**: Default page length 25 (up from 10)

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

- [x] **PRINT-01**: Print stylesheet (@media print) hiding interactive controls, optimizing table layout
- [x] **PRINT-02**: PDF export button (browser-based via html2pdf.js or window.print)

### Accessibility

- [x] **A11Y-01**: ARIA roles and labels on table (role="table"), filter controls, summary cards
- [x] **A11Y-02**: Keyboard-accessible tooltips (tabindex="0", focus trigger, Escape dismiss)
- [x] **A11Y-03**: Chart text alternatives (data table fallback or aria-label summary)
- [x] **A11Y-04**: Skip-to-content links for keyboard navigation
- [x] **A11Y-05**: Replace emoji link icons with SVG icons + screen-reader text
- [x] **A11Y-06**: Sufficient color contrast (WCAG AA minimum)

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

| Requirement | Phase | Status |
|-------------|-------|--------|
| STACK-01 | Phase 13 | Complete |
| STACK-02 | Phase 13 | Complete |
| STACK-03 | Phase 13 | Complete |
| STACK-04 | Phase 13 | Complete |
| STACK-05 | Phase 13 | Complete |
| HIER-01 | Phase 14 | Complete |
| HIER-02 | Phase 14 | Complete |
| HIER-03 | Phase 14 | Complete |
| COLOR-01 | Phase 14 | Complete |
| COLOR-02 | Phase 14 | Complete |
| COLOR-03 | Phase 14 | Complete |
| COLOR-04 | Phase 14 | Complete |
| TABLE-01 | Phase 15 | Complete |
| TABLE-02 | Phase 15 | Complete |
| TABLE-03 | Phase 15 | Complete |
| TABLE-04 | Phase 15 | Complete |
| TABLE-05 | Phase 15 | Complete |
| TABLE-06 | Phase 15 | Complete |
| TABLE-07 | Phase 15 | Complete |
| TABLE-08 | Phase 15 | Complete |
| TABLE-09 | Phase 15 | Complete |
| TABLE-10 | Phase 15 | Complete |
| FILTER-01 | Phase 16 | Pending |
| FILTER-02 | Phase 16 | Pending |
| FILTER-03 | Phase 16 | Pending |
| FILTER-04 | Phase 16 | Pending |
| VIZ-01 | Phase 16 | Pending |
| VIZ-02 | Phase 16 | Pending |
| VIZ-03 | Phase 16 | Pending |
| VIZ-04 | Phase 16 | Pending |
| A11Y-01 | Phase 17 | Complete |
| A11Y-02 | Phase 17 | Complete |
| A11Y-03 | Phase 17 | Complete |
| A11Y-04 | Phase 17 | Complete |
| A11Y-05 | Phase 17 | Complete |
| A11Y-06 | Phase 17 | Complete |
| PRINT-01 | Phase 17 | Complete |
| PRINT-02 | Phase 17 | Complete |
