# VariantCentrifuge Report UX/UI Assessment

**Assessor:** Senior UI/UX Designer — Scientific & Clinical Genomics Reports
**Date:** 2026-02-16
**Scope:** HTML Individual Report, HTML Cohort Report, Excel Report

---

## Executive Summary

VariantCentrifuge's report system is **functionally solid** — it delivers the right data with interactive tables, configurable filters, and external database links. However, it falls short of modern scientific reporting standards in **information hierarchy, visual design maturity, accessibility, and user workflow support**. The reports feel like "developer-built tools" rather than "designer-crafted experiences," which is common in bioinformatics but increasingly unacceptable as clinicians and non-bioinformaticians become end users.

**Overall Score: 4.8 / 10** — Functional but with significant room for improvement, especially in table design.

---

## Detailed Category Ratings

### 1. Information Architecture & Hierarchy — 4/10

**Current state:**
- The individual report places the **variant table first**, summary cards second, and the impact chart last — this is **inverted** from best practice. Users need context (summary, key metrics) before diving into row-level data.
- The cohort report does better (dashboard → chart → filters → table) but buries filters in a separate section instead of integrating them near the table.
- No clear narrative flow from "overview → insight → detail" (the inverted pyramid pattern recommended by [Nielsen Norman Group](https://www.nngroup.com/articles/dashboards-preattentive/)).

**Problems:**
- Summary cards show only 2 metrics (Total Variants, Unique Genes) — far too sparse for a genomic variant report. Missing: impact breakdown, inheritance pattern distribution, pathogenicity summary, sample count, filter criteria applied.
- No "at a glance" clinical summary — a geneticist opening this report must scroll past hundreds of table rows to understand the result landscape.
- The impact distribution chart is at the **bottom** of the page, after the data table — most users will never scroll to it.

**Best practice:** Place high-level summary → key visualizations → interactive filters → detailed table. The most critical insight should be visible without scrolling ([Smashing Magazine dashboard design](https://www.smashingmagazine.com/2021/11/dashboard-design-research-decluttering-data-viz/)).

---

### 2. Visual Design & Aesthetics — 5/10

**Current state:**
- Clean, modern system font stack (Apple/Segoe UI/Roboto) — good.
- Card-based layout with subtle shadows — acceptable.
- Blue gradient header is visually appealing.
- 13px table font size is appropriate for data-dense tables.

**Problems:**
- **Monotone color palette.** Everything is the same blue (#007bff). No semantic color coding — HIGH impact variants look identical to MODIFIER variants. No visual distinction between pathogenic and benign.
- **Stats cards are generic.** All cards have the same blue left-border. No color coding (e.g., red for HIGH impact count, green for benign).
- **The yellow hover-expand tooltip (#FFFFE0)** looks dated and clashes with the blue theme. Modern tooltips use subtle grays or theme-consistent colors.
- **No dark mode support.** Scientific users often work in low-light environments.
- **Footer is minimal** — "Created with VariantCentrifuge" provides no value. Should contain run parameters, filter criteria, or version info.
- **No favicon or branding** — the report looks anonymous in browser tabs.

**Best practice:** Use semantic color coding for clinical severity (red = pathogenic/HIGH, orange = VUS, green = benign/LOW). Consistent, purposeful use of color improves scanning speed by up to 55% ([dashboard visualization research](https://hopara.io/blog/data-visualization-dashboard/)).

---

### 3. Data Table Design — 4/10

> **Score revised downward** from initial 6/10 after deep-dive analysis of truncation, hover-expand, and content overflow behavior. The current hover-expand implementation is fundamentally broken in multiple scenarios and the table lacks critical enterprise-grade patterns.

**Current state:**
- DataTables integration with sorting, pagination, search, column visibility — good fundamentals.
- Horizontal scrolling for many columns — necessary and implemented.
- Hover-to-expand for truncated cells — *attempted* but problematic (see deep-dive below).
- Compact styling (4px/6px padding) — appropriate for genomic data density.

#### 3A. Truncation & Hover-Expand — Critical Issues

The current hover-expand system uses `position: absolute` with a yellow (#FFFFE0) tooltip overlay. This approach has **multiple failure modes:**

**Bug 1 — Hover tooltip clips at table boundaries.**
The expanded `.hover-expand-content.expanded` uses `position: absolute` relative to the table cell. When a cell is near the right edge of the visible table area, the expanded tooltip overflows off-screen. There is no viewport-aware repositioning — the tooltip blindly expands rightward regardless of available space. This is a [known anti-pattern](https://github.com/schrodinger/fixed-data-table-2/issues/48) in table tooltip design.

**Bug 2 — Hover doesn't work with horizontal scroll.**
When `scrollX: true` is enabled, the DataTables wrapper has `overflow: hidden` on the outer container. The `position: absolute` tooltip is clipped by this overflow boundary. Users hovering over truncated cells in scrolled-right columns may see **partially visible or completely invisible** expanded content.

**Bug 3 — Hover is mouse-only — inaccessible.**
The expand is triggered solely by `mouseenter`/`mouseleave` events. There is no `focus`/`blur` handler, no `tabindex`, no keyboard activation. This means:
- Keyboard-only users **cannot access truncated content at all**
- Screen reader users get only the truncated text
- Touch device users (iPad in clinical settings) have **no way to expand**
- This violates [WCAG 2.1 SC 1.4.13 Content on Hover or Focus](https://www.w3.org/WAI/WCAG21/Understanding/content-on-hover-or-focus.html)

**Bug 4 — Tooltip obscures adjacent cells.**
The expanded tooltip (`z-index: 1000`) covers neighboring cells. If a user needs to compare the expanded value with the cell next to it, the tooltip blocks the view. There's no way to "pin" or dismiss the tooltip — it vanishes the moment the cursor moves away.

**Bug 5 — No indication of truncation before hover.**
The `can-expand` class (cursor: help) is only applied *after* DataTables draw callback checks `offsetWidth < scrollWidth`. But this check can give false negatives when content fits the `max-width` but is still meaningfully shortened. There is no consistent visual indicator (e.g., ellipsis with a subtle icon) that a cell contains more content.

**Bug 6 — Yellow tooltip clashes with theme and feels dated.**
The `#FFFFE0` background with `#666` border looks like a Windows XP-era tooltip. Modern design systems (Carbon, PatternFly, Material) use theme-consistent, subtle surfaces with proper shadows.

**What best practices say:**

Per [PatternFly's truncation guidelines](https://www.patternfly.org/ux-writing/truncation/):
- Truncated items **must** always include a tooltip showing the full string
- Preserve **no fewer than 4 characters** when truncating
- Never truncate column headings
- Provide **at least one method** for users to view the entire string
- For strings over ~150 characters, use expandable rows or overlays instead of tooltips

Per [Carbon Design System's overflow content pattern](https://carbondesignsystem.com/patterns/overflow-content/):
- Do **not** include interactive elements within tooltips
- Tooltips should be triggered by **both hover and focus**
- Use popovers (not absolute-positioned spans) for content that needs to persist

Per [Tippy.js accessibility docs](https://atomiks.github.io/tippyjs/v5/accessibility/):
- Tooltips should be applied to focusable elements (`tabindex="0"`)
- Must be dismissible with Escape key
- Must work with keyboard navigation and touch

#### 3B. Column & Row Structure Issues

- **No frozen/sticky first column.** When scrolling horizontally across 30+ columns, users lose track of which variant they're looking at. The GENE or VAR_ID column should be pinned. Per [Pencil & Paper](https://www.pencilandpaper.io/articles/ux-pattern-analysis-enterprise-data-tables): *"In a horizontal scroll situation, having the leftmost column sticky is just as important as the fixed header."*
- **No row highlighting/striping.** Alternating row colors (zebra striping) are absent, making it hard to track across wide tables ([UX Booth data tables](https://uxbooth.com/articles/designing-user-friendly-data-tables/)).
- **No conditional formatting.** IMPACT column values (HIGH, MODERATE, LOW, MODIFIER) are plain text — should be color-coded badges/chips. Same for inheritance patterns.
- **Header row doesn't stand out enough.** Light gray (#e9ecef) on white is too subtle for a 30+ column table.
- **Column order is data-driven, not user-driven.** Columns appear in TSV order, but clinically relevant columns (GENE, IMPACT, HGVS, inheritance) should be grouped and prioritized.
- **No column grouping.** Related columns (all dbNSFP scores, all allele frequencies) should be visually grouped with spanning headers.
- **PageLength defaults to 10.** For genomic data review, 25 or 50 is more practical — geneticists typically review dozens of variants at once.
- **No row selection or bulk actions.** Users can't select rows to export a subset, flag variants for review, or compare selected variants.
- **No expandable row detail.** Clicking a row should expand inline to show all variant details in a readable key-value format — the most scalable pattern for handling tables with 30+ columns where most are hidden by default ([Modern Enterprise UI Tables](https://medium.com/pulsar/modern-enterprise-ui-design-part-1-tables-ad8ee1b9feb)).
- **No content density toggle.** Users should be able to switch between condensed (40px rows), regular (48px), and relaxed (56px) modes depending on their task — scanning vs. reading ([Pencil & Paper density patterns](https://www.pencilandpaper.io/articles/ux-pattern-analysis-enterprise-data-tables)).
- **No column resize handles.** Users cannot adjust column widths to fit their data. Especially important for columns like HGVS_C where content length varies dramatically.

**Best practice:** Sticky columns, zebra striping, semantic color coding, expandable rows, and accessible tooltips are enterprise table design essentials per [Nielsen Norman Group](https://www.nngroup.com/articles/data-tables/), [Stéphanie Walter's enterprise table guide](https://stephaniewalter.design/blog/essential-resources-design-complex-data-tables/), and [Denovers enterprise UX](https://www.denovers.com/blog/enterprise-table-ux-design).

---

### 4. Data Visualization — 4/10

**Current state:**
- Single bar chart (Impact Distribution) using Plotly.
- Cohort report has a horizontal bar chart (Top 20 Genes by Sample Count) with click-to-filter.

**Problems:**
- **Only 1 chart in the individual report** — a variant report should visualize: impact distribution, chromosome distribution, variant type breakdown (SNV/indel/MNV), allele frequency spectrum, inheritance pattern distribution.
- **The bar chart is basic.** No color coding by impact severity (all bars are the same blue). HIGH should be red, MODERATE orange, LOW yellow, MODIFIER gray.
- **No interactive tooltips on chart bars** showing exact counts + percentages.
- **No chromosome ideogram or genomic position visualization.** This is a genomic tool — users expect to see where variants fall across the genome.
- **No allele frequency distribution plot.** AF is one of the most critical variant filters; visualizing the distribution helps users calibrate filters.
- **The chart container is fixed at 400px height** which may be too much for a simple 4-bar chart but too little for a more complex visualization.
- **Cohort chart lacks color coding** — all genes are the same blue regardless of the number of distinct variants.

**Best practice:** Scientific reports should include 3-5 key visualizations above the detail table. Each chart should use color purposefully and support interaction (filter, drill-down). See [Plotly scientific visualization examples](https://plotly.com/python/scientific-charts/).

---

### 5. Filtering & Search — 7/10

**Current state (cohort report — individual report has only DataTables search):**
- Smart column type detection (numeric, categorical, text) — impressive.
- noUiSlider for numeric range filtering — good UX pattern.
- Categorical dropdowns for low-cardinality columns.
- "Include missing values" checkboxes — excellent for genomic data where missingness is meaningful.
- Reset all filters button.
- Dynamic chart/dashboard updates on filter changes.
- Click-to-filter from gene chart.

**Problems:**
- **Individual report has NO column-level filtering.** Only a global search box. This is a major gap — users need to filter by GENE, IMPACT, AF, inheritance pattern independently.
- **No saved/preset filter profiles.** Common use cases (e.g., "pathogenic rare variants," "de novo candidates," "compound het") should be one-click presets.
- **No active filter indicator.** When filters are applied, there's no clear visual indicator of what's currently filtered (no filter chips/tags).
- **Numeric filter labels show raw column names** (e.g., "dbNSFP_CADD_phred") — should show human-readable labels.
- **No filter count indicator** — users can't see "Showing 23 of 1,247 variants" prominently.
- **Text filters require exact typing** — no autocomplete or fuzzy matching for gene names.

**Best practice:** Filters should be discoverable, show active state clearly, and allow quick reset. Preset filter combinations dramatically improve workflow speed ([Justinmind data table filtering](https://www.justinmind.com/ui-design/data-table)).

---

### 6. Accessibility — 3/10

**Current state:**
- Semantic HTML structure (header, table, footer).
- Alt text absent for all visual elements.
- Title attributes on hover-expand cells ("Hover to see full content").

**Problems:**
- **No ARIA roles or labels.** The table lacks `role="grid"`, filter controls lack `aria-label`, stats cards lack semantic meaning for screen readers ([Yale accessibility guide](https://usability.yale.edu/web-accessibility/articles/tables)).
- **Color is the only differentiator** in several places (link icons, IGV vs. standard links). Fails WCAG 2.1 "Use of Color" criterion.
- **Hover-expand tooltips are mouse-only.** No keyboard activation, no focus management. Keyboard users and screen reader users cannot access truncated content.
- **Contrast issues.** The gray text (#6c757d) on white cards has a contrast ratio of ~4.6:1 — passes AA but fails AAA. The header subtitle (white at 0.9 opacity on blue) is borderline.
- **No skip-to-content links.** For reports with 1000+ variants, keyboard navigation is impractical without skip links.
- **Charts have no text alternative.** The Plotly chart is purely visual — no data table fallback or `aria-label` summary.
- **Emoji link icons (🔗)** may not render consistently across platforms and are not meaningful to screen readers.

**Best practice:** Scientific reports must be accessible to all researchers. WCAG 2.1 AA compliance is the minimum standard. See [W3C WAI table tutorials](https://www.w3.org/WAI/tutorials/tables/).

---

### 7. Responsiveness & Print — 3/10

**Current state:**
- max-width: 95% provides some margin.
- Plotly charts are responsive.
- DataTables scrollX handles horizontal overflow.

**Problems:**
- **No mobile/tablet layout.** The 95% width with 30+ column tables is unusable on anything smaller than a desktop.
- **No print stylesheet.** Printing this report produces a broken layout with cut-off tables and invisible interactive elements. Geneticists frequently print or PDF reports for clinical review boards.
- **No PDF export option.** Clinical workflows require PDF attachments for patient records and multidisciplinary team meetings.
- **Stats cards don't stack gracefully** on narrow viewports — the grid breaks at small sizes.
- **No breakpoint-aware column hiding.** On smaller screens, less important columns should auto-hide.

**Best practice:** Scientific reports need print-first design or PDF export. Clinical reports are legal documents in many contexts. A `@media print` stylesheet is essential ([HTML for Bioinformatics guide](https://omicstutorials.com/html-for-bioinformatics-a-comprehensive-guide-for-structuring-and-presenting-scientific-data/)).

---

### 8. Performance & Loading — 6/10

**Current state:**
- CDN-loaded dependencies (DataTables, Plotly, jQuery).
- Deferred rendering in DataTables (`deferRender: true` in cohort).
- Individual report uses server-side rendering (Jinja2 templates) — all HTML pre-generated.

**Problems:**
- **All variant data is embedded in the HTML file.** For large datasets (10,000+ variants), this creates multi-megabyte HTML files that are slow to open and parse.
- **No lazy loading for the chart.** Plotly loads even if the user never scrolls to the chart section.
- **plotly-latest.min.js is ~3.5MB.** Loading the full Plotly bundle for a single bar chart is wasteful. Should use a minimal bundle or a lightweight chart library.
- **jQuery dependency (97KB)** for what could be vanilla JS — adds unnecessary weight.
- **No loading indicator.** Users see a blank page while DataTables initializes large datasets.
- **DataTables 1.10.20 (individual) vs 1.13.6 (cohort)** — version inconsistency; 1.10.20 is outdated (2019).

**Best practice:** For data-heavy scientific reports, consider virtual scrolling, paginated data loading, or Web Worker-based parsing. Keep initial bundle under 500KB for fast load ([web performance best practices](https://www.nngroup.com/articles/complex-application-design/)).

---

### 9. Interactivity & Workflow Support — 5/10

**Current state:**
- Column show/hide toggle.
- Sortable columns.
- Hover-expand for truncated cells.
- External database links (SpliceAI, Franklin, Varsome, gnomAD, ClinVar) — excellent for clinical workflow.
- IGV report links per sample — very useful.

**Problems:**
- **No variant bookmarking/flagging.** Geneticists review hundreds of variants and need to mark candidates for follow-up. No way to flag, star, or annotate variants in the report.
- **No variant comparison mode.** Can't select 2-3 variants to compare side-by-side.
- **No keyboard shortcuts.** j/k for next/previous row, Enter to expand, Esc to close — standard in data-heavy tools.
- **No copy-to-clipboard for individual cells.** Users frequently need to copy a variant ID, HGVS notation, or gene name.
- **No shareable filtered views.** Filter state isn't reflected in the URL — can't share a link to "all HIGH impact variants in TP53."
- **External links open in new tabs** (good) but there's no visual grouping of related links per variant.
- **No variant detail panel.** Clicking a row should expand an inline detail view or side panel showing all fields for that variant in a readable format.

**Best practice:** Interactive reports should support the user's actual workflow (review → flag → discuss → report). URL-based state enables collaboration ([UX Planet data tables](https://uxplanet.org/best-practices-for-usable-and-efficient-data-table-in-applications-4a1d1fb29550)).

---

### 10. Clinical/Scientific Context — 5/10

**Current state:**
- ACMG-relevant external links (ClinVar, Varsome).
- Impact categorization (HIGH/MODERATE/LOW/MODIFIER).
- Inheritance pattern information available.
- Gene burden analysis results in Excel.

**Problems:**
- **No ACMG classification badges.** If ClinVar_CLNSIG data is present, it should be displayed as colored badges (Pathogenic = red, VUS = orange, Benign = green) matching [ACMG/AMP standards](https://www.nature.com/articles/gim201530).
- **No variant interpretation aids.** No inline ACMG criteria evidence, no strength indicators, no classification decision support.
- **No phenotype/HPO term display.** If phenotype data is loaded in the pipeline, it's not surfaced in the HTML report.
- **No sample pedigree visualization.** For family studies, a mini pedigree diagram would provide essential context.
- **No gene-level summary panel.** Clicking a gene name should show: OMIM link, known disease associations, constraint scores (pLI, LOEUF), previous hits in this cohort.
- **Inheritance pattern not visually prominent.** This is one of the most clinically important fields but is treated like any other text column.
- **No quality metrics display.** No variant quality scores, read depth indicators, or confidence metrics visualized.
- **Report metadata (filter criteria, VCF source, reference genome) is buried** in a separate Metadata sheet in Excel, not shown in HTML.

**Best practice:** Clinical genomic reports should follow [ACMG reporting guidelines](https://pubmed.ncbi.nlm.nih.gov/25741868/) with clear pathogenicity indicators and interpretive context.

---

## Score Summary

| Category | Score | Weight | Weighted |
|---|---|---|---|
| 1. Information Architecture & Hierarchy | 4/10 | High | 4.0 |
| 2. Visual Design & Aesthetics | 5/10 | Medium | 5.0 |
| 3. Data Table Design | 4/10 | High | 4.0 |
| 4. Data Visualization | 4/10 | Medium | 4.0 |
| 5. Filtering & Search | 7/10 | High | 7.0 |
| 6. Accessibility | 3/10 | High | 3.0 |
| 7. Responsiveness & Print | 3/10 | Medium | 3.0 |
| 8. Performance & Loading | 6/10 | Medium | 6.0 |
| 9. Interactivity & Workflow | 5/10 | Medium | 5.0 |
| 10. Clinical/Scientific Context | 5/10 | High | 5.0 |
| **Overall (weighted average)** | | | **4.8/10** |

---

## Priority Improvement Recommendations

### Tier 1 — High Impact, Moderate Effort

These changes would most dramatically improve the user experience:

#### 1.1 Restructure Information Hierarchy (Individual Report)
Move summary cards and impact chart **above** the variant table. Add more summary metrics:
- Impact breakdown with colored badges
- Inheritance pattern distribution
- Top affected genes (mini list)
- Filter criteria used
- Sample count and VCF source

#### 1.2 Add Semantic Color Coding
Implement a consistent color system:
- **IMPACT:** HIGH = `#dc3545` (red), MODERATE = `#fd7e14` (orange), LOW = `#ffc107` (amber), MODIFIER = `#6c757d` (gray)
- **ACMG/ClinVar:** Pathogenic = red badge, Likely Pathogenic = orange, VUS = yellow, Likely Benign = light green, Benign = green
- **Inheritance:** De novo = red, Compound het = purple, AD = blue, AR = green, X-linked = teal
- Apply these as colored chips/badges in table cells and summary cards.

#### 1.3 Comprehensive Table Redesign

This is the **single most impactful improvement** — the table is where users spend 90% of their time.

**A. Replace hover-expand with proper tooltip system:**
- Remove the custom `position: absolute` hover-expand entirely
- Implement [Tippy.js](https://atomiks.github.io/tippyjs/) (7KB gzipped) or [Floating UI](https://floating-ui.com/) for viewport-aware tooltips
- Tooltips must be triggered by **both hover and focus** (keyboard accessible)
- Add `tabindex="0"` to truncated cells so they're keyboard-focusable
- Tooltips auto-position (flip to top/left when near edges) — no more off-screen clipping
- Dismiss with Escape key per WCAG 2.1 SC 1.4.13
- Style with theme-consistent colors (white bg, subtle shadow, `#333` text) instead of yellow

**B. Sticky first column + enhanced header:**
- Pin GENE (or first identifier column) using DataTables [FixedColumns](https://datatables.net/extensions/fixedcolumns/) extension
- Strengthen header: `#343a40` background, white text, bottom border `2px solid #007bff`
- Add alternating row backgrounds: `#f8f9fa` / `#ffffff` (zebra striping)
- Active row highlight on hover: `#e8f0fe` (light blue)

**C. Expandable row detail:**
- Add a chevron (▶) in the first column — clicking expands the row inline
- Expanded view shows ALL fields for that variant in a 2-column key-value layout
- Groups related fields: Identifiers | Annotations | Scores | Links | Quality
- This replaces the need to scroll right through 30+ columns for most review tasks
- Use DataTables [RowDetails](https://datatables.net/examples/api/row_details.html) pattern

**D. Content density toggle:**
- Add a toolbar control: Compact | Regular | Relaxed
- Compact: 32px row height, 11px font — for scanning large result sets
- Regular: 44px row height, 13px font — default working mode
- Relaxed: 56px row height, 14px font — for detailed reading
- Persist preference in localStorage

**E. Column width strategy:**
- Set intelligent defaults per column type:
  - CHROM: 80px fixed
  - POS: 100px fixed
  - GENE: 100px min, grow
  - IMPACT: 90px (badge width)
  - HGVS_C/P: 150px with truncation + tooltip
  - AF columns: 80px right-aligned
  - Score columns: 70px right-aligned
  - Link columns: 50px (icon only)
- Add resize handles on column borders (DataTables ColReorder)
- Remember user column width adjustments in localStorage

#### 1.4 Add Column-Level Filtering to Individual Report
Port the cohort report's filter system (sliders, dropdowns, text search) to the individual report. This is the single biggest feature gap.

---

### Tier 2 — High Impact, Higher Effort

#### 2.1 Add Print Stylesheet + PDF Export
```css
@media print {
    header { background: #333 !important; -webkit-print-color-adjust: exact; }
    .chart-section { page-break-inside: avoid; }
    #variants_table { font-size: 10px; }
    .dataTables_wrapper .dataTables_filter,
    .dataTables_wrapper .dataTables_length,
    .dataTables_wrapper .dataTables_paginate { display: none; }
}
```
Add a "Download PDF" button using browser print or a library like html2pdf.js.

#### 2.2 Expand Visualizations
Add to the individual report:
- **Chromosome distribution** (horizontal bar or ideogram)
- **Variant type pie chart** (SNV vs. indel vs. MNV)
- **Allele frequency histogram** (log-scale)
- Use color-coded bars matching the semantic color system.

#### 2.3 Variant Detail Side Panel (Advanced)
Beyond the expandable row (Tier 1.3C), add a slide-out side panel for deep variant review:
- Triggered by double-click or dedicated "detail" button per row
- Shows all fields in organized sections with copy-to-clipboard
- All external links grouped as icon buttons (SpliceAI, Franklin, Varsome, gnomAD, ClinVar)
- IGV links prominently displayed with thumbnail preview
- Stays open while user navigates table rows (updates content on row selection)
- Width: ~400px, slides from right edge — most scalable pattern per [Pencil & Paper](https://www.pencilandpaper.io/articles/ux-pattern-analysis-enterprise-data-tables)

#### 2.4 ARIA Accessibility Pass
- Add `role="grid"` to the data table
- Add `aria-label` to all interactive controls
- Make hover-expand work with keyboard (Enter/Space to toggle)
- Add chart data table alternatives
- Replace emoji 🔗 with SVG icons with `aria-hidden="true"` + screen-reader text

---

### Tier 3 — Nice-to-Have Enhancements

#### 3.1 Filter Presets
Add one-click filter presets:
- "Pathogenic Rare Variants" (AF < 0.01, IMPACT = HIGH)
- "De Novo Candidates" (inheritance = de_novo)
- "Compound Het Candidates" (inheritance = compound_het)
- Custom user-defined presets

#### 3.2 URL-Based Filter State
Encode filter state in URL query parameters so filtered views can be shared:
`report.html?impact=HIGH&gene=TP53&af_max=0.01`

#### 3.3 Performance Optimization
- Replace Plotly with a lighter charting library (Chart.js ~65KB vs Plotly ~3.5MB) for simple charts
- Replace jQuery with vanilla JS (or upgrade DataTables to v2 which is jQuery-optional)
- Add a loading skeleton/spinner during DataTable initialization
- For large datasets, consider virtual scrolling (e.g., DataTables Scroller extension)

#### 3.4 Dark Mode
Add a theme toggle with CSS custom properties:
```css
:root { --bg-primary: #ffffff; --text-primary: #333333; }
[data-theme="dark"] { --bg-primary: #1a1a2e; --text-primary: #e0e0e0; }
```

#### 3.5 Variant Flagging/Notes
Add localStorage-based variant flagging:
- Star/flag button per row
- Notes field per variant
- Export flagged variants list
- Persistent across page reloads (localStorage)

#### 3.6 Gene Info Tooltips
On hover over gene names, show a mini-panel with:
- OMIM disease association
- pLI / LOEUF constraint scores
- Link to GeneCards, OMIM, UniProt

---

## Comparison with Industry Tools

| Feature | VariantCentrifuge | Franklin | Varsome | VarSeq | IGV.js |
|---|---|---|---|---|---|
| Color-coded pathogenicity | No | Yes | Yes | Yes | N/A |
| Sticky columns | No | Yes | Yes | Yes | N/A |
| Filter presets | No | Yes | Yes | Yes | N/A |
| Print/PDF support | No | Yes | Yes | Yes | Limited |
| Variant detail panel | No | Yes | Yes | Yes | Yes |
| Expandable row detail | No | Yes | Yes | Yes | N/A |
| ACMG badges | No | Yes | Yes | Yes | N/A |
| Keyboard shortcuts | No | Limited | Limited | Yes | Yes |
| Accessible tooltips | No (mouse-only) | Yes | Yes | Yes | Fair |
| Viewport-aware tooltips | No (clips at edges) | Yes | Yes | Yes | Yes |
| Content density toggle | No | No | No | Yes | No |
| Column resizing | No | Yes | Yes | Yes | N/A |
| Accessibility (WCAG) | Poor | Good | Good | Good | Fair |
| Pedigree visualization | No | Yes | No | Yes | No |

---

## Conclusion

VariantCentrifuge's report system has a strong **data foundation** — the right information is captured and made available. The main gaps are in **presentation, hierarchy, table design, and clinical context**. The Tier 1 improvements (restructure hierarchy, add color coding, comprehensive table redesign, individual report filtering) would elevate the report from a 4.8 to approximately a **7.5/10** with moderate development effort and no architectural changes.

The cohort report is notably more advanced than the individual report — porting cohort features (advanced filtering, dashboard, dynamic updates) back to the individual report should be a priority.

---

## Appendix A: Data Table Design Deep-Dive

This appendix provides detailed design specifications for the table redesign recommended in Tier 1.3. It covers every aspect of the table from truncation to tooltips to layout strategy.

### A1. The Truncation Problem — Current vs. Proposed

#### Current Implementation (Broken)

```
index.html lines 121-153:
- .hover-expand-content uses display: inline-block; overflow: hidden; text-overflow: ellipsis
- .hover-expand-content.expanded uses position: absolute with z-index: 1000
- Triggered by mouseenter/mouseleave only
- Background: #FFFFE0 (light yellow)
- Max-width: 400px hardcoded
```

**Failure scenarios documented:**

| Scenario | Expected | Actual |
|---|---|---|
| Cell near right table edge | Tooltip visible | Tooltip overflows off-screen, clipped by scrollX container |
| User scrolls table right, hovers | Tooltip in correct position | Tooltip may render at wrong position (absolute coords don't account for scroll offset) |
| Keyboard navigation to cell | Content accessible | No focus handler — content invisible to keyboard users |
| Touch device (iPad) | Tap to expand | No touch handler — content inaccessible |
| Two adjacent truncated cells | Both expandable | First tooltip covers second cell; no way to see both |
| Cell with 200+ char content | Readable tooltip | 400px max-width forces tiny text or excessive wrapping |
| Bottom-row cell | Tooltip below | Tooltip may extend below viewport with no flip logic |

#### Proposed Implementation

**Technology:** [Tippy.js v6](https://atomiks.github.io/tippyjs/) (7KB gzipped, built on Floating UI)

**Why Tippy.js:**
- Viewport-aware auto-positioning (flips, shifts to stay visible)
- Keyboard accessible out of the box (focus trigger)
- Touch device support (tap to show, tap elsewhere to dismiss)
- Accessible: adds `aria-describedby` automatically
- Dismissible with Escape key (WCAG 2.1 SC 1.4.13)
- Theming system for consistent styling
- Animation support (fade, scale, shift)
- Can handle interactive content (links inside tooltips)
- Tiny bundle size — smaller than the current jQuery overhead

**Implementation sketch:**

```javascript
// In DataTables drawCallback or createdCell:
import tippy from 'tippy.js';

// For each truncated cell:
columnDefs.push({
    targets: truncatedColumnIndices,
    createdCell: function(td, cellData) {
        if (cellData && String(cellData).length > threshold) {
            const span = td.querySelector('.truncated');
            span.setAttribute('tabindex', '0'); // Keyboard focusable

            tippy(span, {
                content: cellData,           // Full untruncated text
                placement: 'bottom-start',   // Default position
                flip: true,                  // Auto-flip when near edges
                maxWidth: 500,               // Wider than current 400px
                theme: 'variant-detail',     // Custom theme
                trigger: 'mouseenter focus', // Both mouse AND keyboard
                interactive: false,          // No clickable content inside
                appendTo: document.body,     // Escape overflow:hidden containers
                delay: [200, 0],             // 200ms show delay, instant hide
                aria: {
                    content: 'describedby',  // Proper ARIA association
                },
            });
        }
    }
});
```

**CSS for truncated cells:**
```css
/* Truncated cell indicator */
.truncated {
    display: inline-block;
    max-width: var(--col-max-width, 150px);
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    vertical-align: bottom;
}

/* Focus ring for keyboard users */
.truncated:focus {
    outline: 2px solid #007bff;
    outline-offset: 1px;
    border-radius: 2px;
}

/* Tippy theme - consistent with report design */
.tippy-box[data-theme~='variant-detail'] {
    background-color: #ffffff;
    color: #333333;
    border: 1px solid #dee2e6;
    border-radius: 6px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
    font-size: 13px;
    line-height: 1.5;
    padding: 8px 12px;
    max-width: 500px;
    word-break: break-word;
}

.tippy-box[data-theme~='variant-detail'] .tippy-arrow {
    color: #ffffff;
}

/* For HGVS notation — use monospace */
.tippy-box[data-theme~='variant-detail'] .hgvs-content {
    font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace;
    font-size: 12px;
}
```

### A2. Truncation Strategy per Column Type

Not all columns should truncate the same way. Per [PatternFly guidelines](https://www.patternfly.org/ux-writing/truncation/):

| Column | Strategy | Max Width | Truncation Position | Tooltip |
|---|---|---|---|---|
| GENE | No truncation (short names) | auto | N/A | No |
| CHROM | No truncation | 80px | N/A | No |
| POS | No truncation | 100px | N/A | No |
| IMPACT | Badge — never truncate | 90px | N/A | No |
| HGVS_C | End truncation with ellipsis | 180px | End | Yes — monospace |
| HGVS_P | End truncation with ellipsis | 150px | End | Yes — monospace |
| VAR_ID | **Middle** truncation | 120px | Middle | Yes |
| GT | End truncation | 250px (configurable) | End | Yes — structured (sample list) |
| REF/ALT | End truncation | 80px | End | Yes — monospace |
| EFFECT | End truncation | 150px | End | Yes |
| ClinVar_CLNSIG | Badge — never truncate | auto | N/A | No (use badge) |
| dbNSFP scores | No truncation (numbers) | 70px right-aligned | N/A | No |
| AF columns | No truncation (numbers) | 80px right-aligned | N/A | No |
| Link columns | Icon only — no text | 40px | N/A | Yes (full URL) |

**Middle truncation** for VAR_ID is important because identifiers like `chr1:12345678:ATCGATCG:G` have the most unique part at the end. Display as `chr1:123...TCG:G` so users can differentiate variants. Implementation:

```javascript
function middleTruncate(text, maxLen) {
    if (text.length <= maxLen) return text;
    const keep = Math.floor((maxLen - 3) / 2);
    return text.slice(0, keep) + '...' + text.slice(-keep);
}
```

### A3. Expandable Row Detail Pattern

For tables with 30+ columns where most are hidden by default, the **expandable row** is the industry-standard solution. It replaces the need to scroll horizontally or use hidden column tooltips.

**Design specification:**

```
┌──────────────────────────────────────────────────────────┐
│ ▶ BRCA1  chr17:41234567  c.5266dupC  HIGH  De novo  ... │  <- Collapsed row
├──────────────────────────────────────────────────────────┤
│ ▼ TP53   chr17:7577120   c.818G>A    HIGH  AD       ... │  <- Expanded row
│ ┌──────────────────────────────────────────────────────┐ │
│ │  Identifiers          │  Annotations                │ │
│ │  CHROM: chr17         │  EFFECT: missense_variant   │ │
│ │  POS: 7577120         │  IMPACT: HIGH               │ │
│ │  REF: G               │  GENE: TP53                 │ │
│ │  ALT: A               │  HGVS_C: c.818G>A          │ │
│ │  VAR_ID: chr17:757... │  HGVS_P: p.Arg273His       │ │
│ │                        │                             │ │
│ │  Scores               │  Population Freq            │ │
│ │  CADD: 35.0           │  gnomAD AF: 0.00001        │ │
│ │  REVEL: 0.932         │  ExAC AF: 0.00002          │ │
│ │  SIFT: deleterious    │                             │ │
│ │                        │                             │ │
│ │  External Links                                     │ │
│ │  [ClinVar] [Varsome] [gnomAD] [Franklin] [SpliceAI]│ │
│ │                                                      │ │
│ │  IGV Reports                                        │ │
│ │  [Sample1 🔬] [Sample2 🔬] [Sample3 🔬]            │ │
│ └──────────────────────────────────────────────────────┘ │
│ ▶ CFTR   chr7:117199646  c.1521_...  HIGH  AR       ... │  <- Collapsed row
└──────────────────────────────────────────────────────────┘
```

**Implementation:** DataTables child row API:
```javascript
$('#variants_table tbody').on('click', 'td.details-control', function () {
    const tr = $(this).closest('tr');
    const row = table.row(tr);

    if (row.child.isShown()) {
        row.child.hide();
        tr.removeClass('shown');
    } else {
        row.child(formatDetailPanel(row.data())).show();
        tr.addClass('shown');
    }
});

function formatDetailPanel(data) {
    // Returns HTML with grouped key-value pairs
    // Uses a 2-column CSS grid for compact display
    return `<div class="variant-detail-panel">
        <div class="detail-section">
            <h4>Identifiers</h4>
            <dl>
                <dt>CHROM</dt><dd>${data.CHROM}</dd>
                <dt>POS</dt><dd>${data.POS}</dd>
                ...
            </dl>
        </div>
        ...
    </div>`;
}
```

### A4. Sticky Column Implementation

**Problem:** With `scrollX: true` and 30+ columns, horizontal scrolling causes users to lose track of which gene/variant they're viewing.

**Solution:** Use DataTables [FixedColumns](https://datatables.net/extensions/fixedcolumns/) extension:

```javascript
$('#variants_table').DataTable({
    scrollX: true,
    fixedColumns: {
        left: 2,  // Pin GENE + expand chevron
        right: 0
    }
});
```

**CSS for frozen column visual separation:**
```css
/* Visual divider between frozen and scrollable columns */
.dtfc-fixed-left {
    border-right: 2px solid #007bff !important;
    box-shadow: 4px 0 8px rgba(0, 0, 0, 0.08);
}

/* Ensure frozen column has opaque background */
table.dataTable.dtfc-has-left tbody tr td.dtfc-fixed-left {
    background-color: inherit; /* Respects zebra striping */
}
```

### A5. Conditional Formatting — Color-Coded Badges

Replace plain text values with colored badges for categorical clinical columns:

```javascript
// In DataTables column render:
render: function(data, type, row) {
    if (type !== 'display') return data; // Preserve raw value for sort/filter

    const badges = {
        // IMPACT
        'HIGH':     '<span class="badge badge-high">HIGH</span>',
        'MODERATE': '<span class="badge badge-moderate">MODERATE</span>',
        'LOW':      '<span class="badge badge-low">LOW</span>',
        'MODIFIER': '<span class="badge badge-modifier">MODIFIER</span>',

        // ClinVar
        'Pathogenic':        '<span class="badge badge-pathogenic">Pathogenic</span>',
        'Likely_pathogenic': '<span class="badge badge-likely-path">Likely Path.</span>',
        'Uncertain_significance': '<span class="badge badge-vus">VUS</span>',
        'Likely_benign':     '<span class="badge badge-likely-benign">Likely Benign</span>',
        'Benign':            '<span class="badge badge-benign">Benign</span>',
    };
    return badges[data] || data;
}
```

```css
.badge {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 12px;
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 0.3px;
    white-space: nowrap;
}

/* Impact severity */
.badge-high      { background: #fce4e4; color: #c0392b; border: 1px solid #e6a0a0; }
.badge-moderate  { background: #fef3e2; color: #d68910; border: 1px solid #f0c987; }
.badge-low       { background: #fef9e7; color: #b7950b; border: 1px solid #f0e68c; }
.badge-modifier  { background: #f0f0f0; color: #666666; border: 1px solid #d0d0d0; }

/* ClinVar classification */
.badge-pathogenic   { background: #fce4e4; color: #c0392b; }
.badge-likely-path  { background: #fdebd0; color: #ca6f1e; }
.badge-vus          { background: #fef9e7; color: #b7950b; }
.badge-likely-benign{ background: #eafaf1; color: #1e8449; }
.badge-benign       { background: #d5f5e3; color: #196f3d; }

/* Inheritance patterns */
.badge-denovo       { background: #fce4e4; color: #c0392b; }
.badge-comp-het     { background: #f0e0f5; color: #7d3c98; }
.badge-ad           { background: #d6eaf8; color: #2471a3; }
.badge-ar           { background: #d5f5e3; color: #1e8449; }
.badge-xlinked      { background: #d0ece7; color: #148f77; }
```

### A6. Table Toolbar Design

The current toolbar (DataTables default DOM `'lBfrtip'`) is generic. Redesign:

```
┌─────────────────────────────────────────────────────────────────┐
│  🔍 Search: [____________]   Showing 23 of 1,247 variants      │
│                                                                  │
│  [Compact] [Regular] [Relaxed]  │  [Show/Hide Cols ▾]  │  [⬇ Export ▾]  │
│  ───────────────────────────────────────────────────────────────  │
│  Active filters: [IMPACT: HIGH ✕] [AF < 0.01 ✕] [Clear all]   │
└─────────────────────────────────────────────────────────────────┘
```

Key elements:
- **Prominent record count** ("Showing X of Y variants") — always visible
- **Density toggle** — three buttons for compact/regular/relaxed
- **Active filter chips** — removable tags showing current filter state
- **Export dropdown** — CSV, Excel, Copy to Clipboard
- **Column manager** — grouped by category (Identifiers, Annotations, Scores, etc.)

### A7. Do's and Don'ts Summary

| DO | DON'T |
|---|---|
| Use viewport-aware tooltips (Tippy.js/Floating UI) | Use `position: absolute` in overflow:hidden containers |
| Make tooltips keyboard-accessible (`tabindex="0"` + focus trigger) | Rely on mouseenter/mouseleave only |
| Show ellipsis (`...`) as visual truncation indicator | Silently clip content with `overflow: hidden` alone |
| Preserve at least 4 characters when truncating ([PatternFly](https://www.patternfly.org/ux-writing/truncation/)) | Truncate to 1-2 characters or empty string |
| Use middle truncation for identifiers (unique part at end) | Always truncate from the end |
| Pin identifier columns when enabling horizontal scroll | Let all columns scroll away together |
| Use zebra striping for wide tables (30+ columns) | Use same background for all rows |
| Add expandable row detail for hidden column content | Force users to scroll right through 30+ columns |
| Use colored badges for categorical clinical values | Show classification as plain text |
| Provide content density toggle (compact/regular/relaxed) | Lock users into one density |
| Allow column resizing with drag handles | Use fixed column widths that can't adapt to data |
| Show active filter state as removable chips | Apply filters silently with no visible indicator |
| Right-align numeric columns (scores, frequencies) | Left-align everything |
| Use monospace font for genetic notation (HGVS, GT) | Use proportional font for sequence data |
| Dismiss tooltips with Escape key (WCAG 2.1 SC 1.4.13) | Trap tooltip open with no keyboard dismiss |
| Append tooltips to `document.body` to escape containers | Append tooltips inside `overflow: hidden` parents |
| Add 150-200ms delay before showing tooltip | Show instantly (causes tooltip flicker during scanning) |
| Never truncate table column headers ([PatternFly](https://www.patternfly.org/ux-writing/truncation/)) | Apply truncation rules to header row |

---

## Sources

- [Nielsen Norman Group — Dashboards: Making Charts Easier to Understand](https://www.nngroup.com/articles/dashboards-preattentive/)
- [Nielsen Norman Group — Data Tables: Four Major User Tasks](https://www.nngroup.com/articles/data-tables/)
- [Nielsen Norman Group — 8 Design Guidelines for Complex Applications](https://www.nngroup.com/articles/complex-application-design/)
- [Pencil & Paper — Data Table Design UX Patterns & Best Practices](https://www.pencilandpaper.io/articles/ux-pattern-analysis-enterprise-data-tables)
- [Stéphanie Walter — Enterprise UX: Essential Resources for Complex Data Tables](https://stephaniewalter.design/blog/essential-resources-design-complex-data-tables/)
- [Smashing Magazine — Dashboard Design: Research, Decluttering, Data Viz](https://www.smashingmagazine.com/2021/11/dashboard-design-research-decluttering-data-viz/)
- [UX Booth — Designing User-Friendly Data Tables](https://uxbooth.com/articles/designing-user-friendly-data-tables/)
- [UX Planet — Best Practices for Data Tables](https://uxplanet.org/best-practices-for-usable-and-efficient-data-table-in-applications-4a1d1fb29550)
- [Justinmind — Designing Effective Data Table UI](https://www.justinmind.com/ui-design/data-table)
- [Yale Usability — Tables & Web Accessibility](https://usability.yale.edu/web-accessibility/articles/tables)
- [HTML for Bioinformatics — Structuring Scientific Data](https://omicstutorials.com/html-for-bioinformatics-a-comprehensive-guide-for-structuring-and-presenting-scientific-data/)
- [ACMG/AMP Standards for Variant Interpretation](https://www.nature.com/articles/gim201530)
- [ClinGen — Sequence Variant Interpretation](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/)
- [UX Design Institute — Top UX Trends 2026](https://www.uxdesigninstitute.com/blog/the-top-ux-design-trends-in-2026/)
- [Denovers — 6 Best Practices for Enterprise Table UX](https://www.denovers.com/blog/enterprise-table-ux-design)
- [PatternFly — Truncation Design Guidelines](https://www.patternfly.org/ux-writing/truncation/)
- [PatternFly — Truncate Component](https://www.patternfly.org/components/truncate/design-guidelines/)
- [Carbon Design System — Overflow Content Patterns](https://carbondesignsystem.com/patterns/overflow-content/)
- [Carbon Design System — Data Table Usage](https://v10.carbondesignsystem.com/components/data-table/usage/)
- [Carbon Design System — Tooltip Usage](https://carbondesignsystem.com/components/tooltip/usage/)
- [Tippy.js — Tooltip and Popover Library](https://atomiks.github.io/tippyjs/)
- [Tippy.js — Accessibility Guide](https://atomiks.github.io/tippyjs/v5/accessibility/)
- [Floating UI — Tooltip Documentation](https://floating-ui.com/docs/tooltip)
- [Rick Strahl — HTML Table Cell Overflow Handling](https://weblog.west-wind.com/posts/2023/Jan/26/HTML-Table-Cell-Overflow-Handling)
- [codegenes — CSS text-overflow in Table Cells](https://www.codegenes.net/blog/css-text-overflow-in-a-table-cell/)
- [Modern Enterprise UI Design — Tables (Pulsar/Medium)](https://medium.com/pulsar/modern-enterprise-ui-design-part-1-tables-ad8ee1b9feb)
- [Data Table Design Patterns (UX Design/Bootcamp)](https://bootcamp.uxdesign.cc/data-table-design-patterns-4e38188a0981)
- [MobileSpoon — 20 Rules for Data Table Design](https://www.mobilespoon.net/2019/11/design-ui-tables-20-rules-guide.html)
- [MindK — 17 UX Tips for Data Table Design](https://www.mindk.com/blog/better-data-table-design/)
- [WCAG 2.1 SC 1.4.13 — Content on Hover or Focus](https://www.w3.org/WAI/WCAG21/Understanding/content-on-hover-or-focus.html)
