# Phase 17: Accessibility and Print/PDF - Research

**Researched:** 2026-02-17
**Domain:** Web accessibility (WCAG 2.1 AA) and print/PDF optimization for single-file HTML reports
**Confidence:** HIGH

## Summary

Phase 17 implements WCAG 2.1 Level AA accessibility compliance and print/PDF optimization for the VariantCentrifuge HTML report. The research reveals that accessibility is primarily a matter of adding appropriate ARIA attributes, ensuring keyboard navigation, providing text alternatives, and maintaining sufficient color contrast. Print optimization requires CSS @media print stylesheets that hide interactive elements and optimize table layouts. PDF export can be achieved through browser-native window.print() functionality with proper print stylesheets, avoiding the complexity and limitations of JavaScript PDF libraries.

**Key findings:**
- WCAG 2.1 AA requires 4.5:1 contrast ratio for normal text, keyboard accessibility, text alternatives for non-text content, and skip navigation mechanisms
- DataTables v2 has built-in ARIA label support for pagination but requires manual ARIA implementation for the table itself (use role="table" not role="grid" for static tables)
- Chart.js requires manual accessibility implementation via aria-label on canvas elements or fallback data tables
- Tippy.js already supports keyboard accessibility with default "mouseenter focus" trigger but requires tabindex="0" on non-focusable elements
- Print stylesheets should use display: none for interactive elements and page-break-inside: avoid for table rows
- Browser-native window.print() is superior to JavaScript PDF libraries (html2pdf.js, jsPDF) for this use case

**Primary recommendation:** Implement accessibility through progressive enhancement - add ARIA attributes to existing components, verify color contrast ratios, implement skip-to-content link, replace emoji icons with SVG + sr-only text, and create a comprehensive @media print stylesheet. Avoid JavaScript PDF libraries in favor of browser-native printing.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| WCAG 2.1 | Level AA | Accessibility standard | International standard, U.S. government requirement as of April 2026 |
| @media print | CSS3 | Print stylesheets | Native CSS, no dependencies, universal browser support |
| window.print() | Browser API | PDF export trigger | Native browser capability, maintains text selectability |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| WebAIM Contrast Checker | Web tool | Verify color contrast | During development to check badge colors meet 4.5:1 ratio |
| axe DevTools | Browser extension | Automated a11y testing | Initial accessibility audit, catches ~30-40% of issues |
| WAVE | Browser extension | Visual accessibility assessment | Secondary validation, visual presentation of issues |
| Google Lighthouse | Browser DevTools | Accessibility scoring | Ongoing monitoring, CI integration possible |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| window.print() | html2pdf.js | html2pdf renders to image (non-selectable text, large files, mobile limits) |
| window.print() | jsPDF | Requires manual layout coding, doesn't preserve HTML semantics |
| role="table" | role="grid" | Grid requires 2D keyboard navigation implementation, overkill for static table with sortable columns |

**Installation:**
No new dependencies required. All accessibility features use standard HTML5/ARIA attributes and CSS. Testing tools are browser extensions or web-based.

## Architecture Patterns

### Recommended Project Structure
```
variantcentrifuge/templates/
├── index.html                    # Single template file (inline assets)
│   ├── <head>                    # Add skip-link styles, print stylesheet
│   ├── <body>                    # Add skip-to-content link as first element
│   ├── <style>                   # Add @media print styles, sr-only class
│   └── <script>                  # Add ARIA attributes via JavaScript after DataTables init

variantcentrifuge/assets/
├── icons/                        # New directory for SVG icons
│   ├── external-link.svg         # Replace emoji link icons
│   └── (other icons as needed)
```

### Pattern 1: Skip-to-Content Link
**What:** A hidden link as the first focusable element that jumps keyboard users to main content
**When to use:** Required for WCAG 2.4.1 (Bypass Blocks - Level A) on any page with repeated navigation
**Example:**
```html
<!-- First element inside <body> -->
<a href="#main-content" class="skip-link">Skip to main content</a>

<!-- Target element -->
<main id="main-content" tabindex="-1">
  <!-- Main content here -->
</main>

<style>
/* Hide off-screen by default */
.skip-link {
    position: absolute;
    left: -9999px;
    top: auto;
    width: 1px;
    height: 1px;
    overflow: hidden;
}

/* Visible when focused */
.skip-link:focus {
    position: fixed;
    top: 10px;
    left: 10px;
    z-index: 9999;
    padding: 10px 20px;
    background-color: #007bff;
    color: white;
    text-decoration: none;
    border-radius: 4px;
    width: auto;
    height: auto;
    overflow: visible;
}
</style>
```

### Pattern 2: ARIA Roles for DataTables
**What:** Add semantic ARIA attributes to the DataTables table after initialization
**When to use:** For any data table that uses DataTables but needs screen reader support
**Example:**
```javascript
// After DataTables initialization
var table = $('#variants_table').DataTable({ /* config */ });

// Add ARIA attributes
$('#variants_table').attr('role', 'table');
$('#variants_table thead').attr('role', 'rowgroup');
$('#variants_table tbody').attr('role', 'rowgroup');
$('#variants_table th').attr('role', 'columnheader');
$('#variants_table td').attr('role', 'cell');

// Add aria-labels to filter controls
$('.dataTables_filter input').attr('aria-label', 'Search variants');
$('.dataTables_length select').attr('aria-label', 'Number of variants per page');

// Add aria-label to column visibility button if present
$('.buttons-colvis').attr('aria-label', 'Show or hide columns');
```
**Important:** Use `role="table"` NOT `role="grid"`. Grid pattern requires 2D arrow key navigation and is for interactive widgets. DataTables with sortable columns is still a table, not a grid.

### Pattern 3: Chart.js Accessibility with Data Table Fallback
**What:** Hide canvas from screen readers, provide accessible data table as alternative
**When to use:** Every Chart.js chart in the report (dashboard section)
**Example:**
```html
<div class="chart-wrapper">
    <h3>Variants by Impact</h3>
    <div class="chart-container">
        <canvas id="impactChart" aria-hidden="true"></canvas>
    </div>
    <!-- Accessible data table fallback -->
    <table class="chart-data-table sr-only" role="table" aria-label="Variants by Impact data">
        <thead>
            <tr>
                <th>Impact</th>
                <th>Count</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>HIGH</td>
                <td>{{ summary.impact_counts.HIGH }}</td>
            </tr>
            <!-- ... other rows ... -->
        </tbody>
    </table>
</div>

<style>
/* Screen reader only class */
.sr-only {
    position: absolute;
    left: -10000px;
    width: 1px;
    height: 1px;
    overflow: hidden;
}
</style>
```

### Pattern 4: SVG Icons with Screen Reader Text
**What:** Replace emoji link icons with SVG icons hidden from screen readers, plus descriptive text
**When to use:** For all external link icons in the report (SpliceAI, Franklin, Varsome, gnomAD, ClinVar)
**Example:**
```html
<!-- Old emoji approach -->
<a href="{{ url }}">🔗</a>

<!-- New accessible approach -->
<a href="{{ url }}" class="external-link">
    <svg aria-hidden="true" class="icon-external" width="16" height="16" viewBox="0 0 16 16">
        <path d="M13 3h-2V1h4v4h-2V3.5L9.5 7 8 5.5 11.5 2H13z"/>
        <path d="M12 10v3H3V4h3V2H3a2 2 0 00-2 2v9a2 2 0 002 2h9a2 2 0 002-2v-3h-2z"/>
    </svg>
    <span class="sr-only">View in {{ tool_name }} (opens in new tab)</span>
</a>

<style>
.external-link {
    display: inline-flex;
    align-items: center;
    gap: 4px;
}
.icon-external {
    flex-shrink: 0;
}
</style>
```

### Pattern 5: Print Stylesheet Structure
**What:** @media print CSS that hides interactive elements and optimizes table layout
**When to use:** Every interactive HTML report that may be printed or exported to PDF
**Example:**
```css
@media print {
    /* Hide interactive controls */
    .dataTables_filter,
    .dataTables_length,
    .dataTables_paginate,
    .dataTables_info,
    .dt-buttons,
    .buttons-colvis,
    .density-toggle,
    .skip-link,
    footer {
        display: none !important;
    }

    /* Prevent page breaks inside rows */
    #variants_table tbody tr {
        page-break-inside: avoid;
        break-inside: avoid;
    }

    /* Repeat headers on each page */
    #variants_table thead {
        display: table-header-group;
    }

    /* Expand fixed columns */
    .dtfc-fixed-left {
        position: static !important;
        box-shadow: none !important;
    }

    /* Collapse detail panels or remove if too verbose */
    .detail-panel {
        page-break-inside: avoid;
        break-inside: avoid;
        font-size: 10pt;
    }

    /* Optimize badge colors for grayscale printing */
    .badge {
        border: 1px solid #333;
        color: #000 !important;
        background-color: #fff !important;
    }

    /* Ensure links show URLs */
    a[href^="http"]:after {
        content: " (" attr(href) ")";
        font-size: 9pt;
        color: #666;
    }
}
```

### Pattern 6: Tippy.js Keyboard Accessibility (Already Implemented)
**What:** Ensure tooltips work with keyboard navigation via tabindex and focus trigger
**When to use:** All Tippy.js tooltips on non-focusable elements
**Example:**
```javascript
// For truncated cells (which are inside <td>, inherently not focusable)
tippy('.truncated-cell', {
    trigger: 'mouseenter focus',  // Responds to both mouse and keyboard
    interactive: false,
    allowHTML: false,
    content(reference) {
        return reference.getAttribute('data-tippy-content');
    }
});

// For standalone tooltips on non-focusable elements
// Add tabindex="0" to the element
<span class="info-icon" tabindex="0" data-tippy-content="Explanation">ℹ️</span>
```
**Note:** VariantCentrifuge already uses `trigger: 'mouseenter focus'` (from Phase 15-02-05), so tooltips are keyboard-accessible. Verify non-focusable trigger elements have `tabindex="0"`.

### Anti-Patterns to Avoid
- **role="grid" on DataTables:** Don't use grid role for tables with sortable columns. Grid requires full 2D arrow-key navigation. Tables can be sortable, filterable, and still be tables.
- **html2pdf.js for reports:** Renders content as images (non-selectable text, accessibility lost, large file sizes, mobile memory limits). Use window.print() instead.
- **Visible-only skip links:** Skip links must be hidden by default but visible on focus. Never use `display: none` (removes from tab order).
- **Generic chart alt text:** Don't use aria-label="Chart" - be specific: "Variants by Impact bar chart showing 45 HIGH, 120 MODERATE"
- **Color-only information:** Don't use color alone to convey meaning. Badges must have text labels, not just background colors.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| PDF generation from HTML | Custom canvas rendering, server-side PDF generation | Browser window.print() with @media print stylesheet | Native browser feature preserves text selectability, accessibility, works offline, no library bloat |
| Screen reader only text | Custom JavaScript to hide/show content | CSS `.sr-only` class with position absolute | Standard pattern, works across all screen readers, no JavaScript needed |
| Contrast ratio validation | Visual inspection, manual calculation | WebAIM Contrast Checker, browser DevTools | Automated tools provide precise ratios, suggest alternative colors, free and fast |
| Accessibility testing | Manual testing only | axe DevTools + manual testing | Automated tools catch ~40% of issues reliably, but manual testing is still required for the other ~60% |
| Skip navigation implementation | Custom scrolling JavaScript | Standard `<a href="#main-content">` link with CSS | Simple, reliable, screen reader compatible, no JavaScript failure modes |
| ARIA attribute management | Custom attribute toggling logic | Declarative ARIA in HTML, set once after DataTables init | Less brittle, easier to audit, follows ARIA authoring practices |
| Keyboard trap prevention | Custom focus management | Ensure all interactive elements are in logical tab order, test with Tab key | Browser handles most cases correctly if DOM order is logical |
| SVG icon accessibility | Custom icon fonts, iframes | Inline SVG with aria-hidden="true" + span.sr-only | No HTTP requests, works offline, cacheable, screen reader compatible |

**Key insight:** Accessibility is primarily additive, not transformative. Existing HTML/CSS/JavaScript needs ARIA attributes and semantic structure added, not wholesale rewriting.

## Common Pitfalls

### Pitfall 1: Using aria-label on non-labelable elements
**What goes wrong:** Adding aria-label to `<div>`, `<span>`, or `<table>` without a role doesn't work - screen readers ignore it
**Why it happens:** Misunderstanding which elements accept aria-label natively
**How to avoid:** Only use aria-label on: interactive elements (button, input, select), elements with appropriate roles (role="img", role="button"), or form controls. For tables, use `<caption>` or `aria-labelledby` pointing to a heading.
**Warning signs:** Screen reader announces element but not the label you added

### Pitfall 2: Insufficient color contrast on colored badges
**What goes wrong:** Semantic color badges (IMPACT: HIGH = red, MODERATE = orange, etc.) fail WCAG AA if background color doesn't have 4.5:1 contrast with white text
**Why it happens:** Colors chosen for visual appeal without contrast verification
**How to avoid:** Test every badge color combination with WebAIM Contrast Checker. For Phase 17, audit all existing badge colors in template (lines 669-747 in index.html). Darken background colors if needed.
**Warning signs:**
- Red (#dc3545) with white text: 4.51:1 (PASS, barely)
- Orange (#fd7e14) with white text: 2.55:1 (FAIL - needs darkening)
- Amber (#f59e0b) with white text: 2.17:1 (FAIL - needs darkening)
**Fix:** Darken failed colors or use black text on light backgrounds

### Pitfall 3: Print stylesheets not hiding FixedColumns duplicates
**What goes wrong:** DataTables FixedColumns creates cloned column DOM elements that appear as duplicates when printed
**Why it happens:** FixedColumns uses absolute positioning, which isn't removed by standard print stylesheets
**How to avoid:** In @media print, set `.dtfc-fixed-left { position: static !important; }` to collapse the fixed column back into normal flow
**Warning signs:** Printed table shows first column twice

### Pitfall 4: Forgetting tabindex="-1" on skip link targets
**What goes wrong:** Skip link href="#main-content" jumps to element, but focus doesn't follow (confusing for keyboard users)
**Why it happens:** Non-interactive elements can't receive focus by default
**How to avoid:** Add `tabindex="-1"` to the target element (e.g., `<main id="main-content" tabindex="-1">`). This allows programmatic focus without adding element to normal tab order.
**Warning signs:** Skip link scrolls page but next Tab press doesn't continue from the target

### Pitfall 5: Chart canvas elements exposed to screen readers
**What goes wrong:** Screen readers try to interpret canvas pixel data, resulting in "unlabeled graphic" or gibberish
**Why it happens:** Forgetting to add aria-hidden="true" to canvas elements
**How to avoid:** Always add `aria-hidden="true"` to Chart.js canvas elements and provide alternative content (data table or detailed aria-label on wrapper)
**Warning signs:** Screen reader users report "unlabeled" or "image" announcements without context

### Pitfall 6: Using display: none or visibility: hidden on skip links
**What goes wrong:** Skip link is removed from keyboard tab order entirely, defeating the purpose
**Why it happens:** Misunderstanding CSS hiding techniques
**How to avoid:** Use off-screen positioning (position: absolute; left: -9999px;) instead of display: none. Bring back on-screen with position: fixed when focused.
**Warning signs:** Tab key doesn't focus the skip link even though it's in the DOM

### Pitfall 7: Over-testing with automated tools, under-testing manually
**What goes wrong:** Achieving 100% Lighthouse accessibility score but still having major usability issues for keyboard/screen reader users
**Why it happens:** Automated tools catch only ~40% of WCAG issues
**How to avoid:** Use axe/WAVE/Lighthouse for initial audit, then manually test: navigate entire page with keyboard only (Tab, Enter, Esc), test with actual screen reader (NVDA/JAWS on Windows, VoiceOver on macOS), verify all interactive elements reachable and usable
**Warning signs:** Perfect automated scores but real users report accessibility problems

## Code Examples

Verified patterns from official sources:

### WCAG AA Minimum Contrast
```css
/* Source: https://www.w3.org/TR/WCAG21/#contrast-minimum
 * Success Criterion 1.4.3 (Level AA)
 * Normal text: 4.5:1 minimum
 * Large text (18pt+ or 14pt+ bold): 3:1 minimum
 */

/* PASS: White on dark blue */
.header {
    background-color: #0056b3;  /* Contrast: 8.6:1 */
    color: white;
}

/* FAIL: White on orange */
.badge-moderate {
    background-color: #fd7e14;  /* Contrast: 2.55:1 - TOO LOW */
    color: white;
}

/* FIX: Darken orange or use black text */
.badge-moderate {
    background-color: #c05000;  /* Contrast: 4.6:1 - PASS */
    color: white;
}
```

### Skip-to-Content Link (WebAIM Pattern)
```html
<!-- Source: https://webaim.org/techniques/skipnav/ -->
<body>
    <a href="#main-content" class="skip-link">Skip to main content</a>

    <header>
        <!-- Logo, navigation, etc. -->
    </header>

    <main id="main-content" tabindex="-1">
        <!-- Main content starts here -->
    </main>
</body>

<style>
/* Off-screen by default */
.skip-link {
    position: absolute;
    left: -10000px;
    top: auto;
    width: 1px;
    height: 1px;
    overflow: hidden;
}

/* Visible when focused */
.skip-link:focus {
    position: fixed;
    top: 10px;
    left: 10px;
    z-index: 9999;
    padding: 12px 20px;
    background-color: #007bff;
    color: white;
    text-decoration: none;
    border-radius: 4px;
    box-shadow: 0 0 10px rgba(0,0,0,0.3);
    width: auto;
    height: auto;
    overflow: visible;
}
</style>
```

### Chart.js with Data Table Fallback
```html
<!-- Source: https://www.chartjs.org/docs/latest/general/accessibility.html -->
<div class="chart-wrapper">
    <h3 id="impact-chart-heading">Variants by Impact</h3>

    <!-- Canvas hidden from screen readers -->
    <canvas id="impactChart"
            aria-hidden="true"
            aria-labelledby="impact-chart-heading">
    </canvas>

    <!-- Alternative content for screen readers -->
    <table class="sr-only" aria-labelledby="impact-chart-heading">
        <thead>
            <tr>
                <th scope="col">Impact Level</th>
                <th scope="col">Variant Count</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>HIGH</td>
                <td>{{ summary.impact_counts.HIGH }}</td>
            </tr>
            <tr>
                <td>MODERATE</td>
                <td>{{ summary.impact_counts.MODERATE }}</td>
            </tr>
            <!-- etc. -->
        </tbody>
    </table>
</div>

<style>
.sr-only {
    position: absolute;
    left: -10000px;
    width: 1px;
    height: 1px;
    overflow: hidden;
}
</style>
```

### SVG Icon with Screen Reader Text
```html
<!-- Source: https://www.scottohara.me/blog/2019/05/22/contextual-images-svgs-and-a11y.html -->

<!-- Decorative icon (hidden from screen readers) -->
<a href="https://spliceailookup.broadinstitute.org/" target="_blank" class="external-link">
    <svg aria-hidden="true" focusable="false" class="icon" width="16" height="16">
        <use xlink:href="#icon-external-link"/>
    </svg>
    <span class="sr-only">View in SpliceAI Lookup (opens in new tab)</span>
    SpliceAI
</a>

<!-- SVG sprite definition (place once in document) -->
<svg xmlns="http://www.w3.org/2000/svg" style="display: none;">
    <symbol id="icon-external-link" viewBox="0 0 16 16">
        <path d="M13.5 1A1.5 1.5 0 0015 2.5v11a1.5 1.5 0 01-1.5 1.5h-11A1.5 1.5 0 011 13.5v-11A1.5 1.5 0 012.5 1h11zm-1 1h-10v12h10V2z"/>
        <path d="M6.5 5.5L10 5.5 10 9M10 5.5L5.5 10"/>
    </symbol>
</svg>

<style>
.external-link {
    display: inline-flex;
    align-items: center;
    gap: 4px;
}
.icon {
    flex-shrink: 0;
    fill: currentColor;
}
</style>
```

### DataTables ARIA Attributes
```javascript
// Source: https://datatables.net/forums/discussion/23920/datatables-and-aria
// Applied after DataTables initialization

var table = $('#variants_table').DataTable({
    // ... config ...
});

// Add table role structure (use 'table' not 'grid' for static tables)
$('#variants_table').attr('role', 'table');
$('#variants_table thead').attr('role', 'rowgroup');
$('#variants_table tbody').attr('role', 'rowgroup');
$('#variants_table thead tr').attr('role', 'row');
$('#variants_table tbody tr').attr('role', 'row');
$('#variants_table th').attr('role', 'columnheader');
$('#variants_table td').attr('role', 'cell');

// Add labels to controls
$('.dataTables_filter input').attr({
    'aria-label': 'Search variants',
    'type': 'search'  // Semantic input type
});
$('.dataTables_length select').attr('aria-label', 'Number of variants per page');
$('.buttons-colvis').attr('aria-label', 'Show or hide columns');

// Add live region for table updates (helps screen readers)
$('.dataTables_info').attr({
    'role': 'status',
    'aria-live': 'polite'
});
```

### Print Stylesheet for DataTables
```css
/* Source: Best practices from DataTables forums and CSS-Tricks */
@media print {
    /* Hide all interactive controls */
    .dataTables_wrapper .dataTables_filter,
    .dataTables_wrapper .dataTables_length,
    .dataTables_wrapper .dataTables_paginate,
    .dataTables_wrapper .dt-buttons,
    .density-toggle,
    .buttons-colvis,
    footer,
    .skip-link {
        display: none !important;
    }

    /* Collapse FixedColumns back to normal flow */
    .dtfc-fixed-left,
    .dtfc-fixed-right {
        position: static !important;
        box-shadow: none !important;
    }

    /* Prevent page breaks inside table rows */
    #variants_table tbody tr {
        page-break-inside: avoid;      /* Legacy */
        break-inside: avoid;            /* Modern */
    }

    /* Repeat table header on each page */
    #variants_table thead {
        display: table-header-group;
    }

    #variants_table tfoot {
        display: table-footer-group;
    }

    /* Optimize for grayscale printing */
    .badge {
        border: 1px solid #333 !important;
        color: #000 !important;
        background-color: #fff !important;
    }

    /* Show URLs for external links */
    a[href^="http"]:not(.skip-link)::after {
        content: " (" attr(href) ")";
        font-size: 9pt;
        color: #555;
    }

    /* Reduce zebra striping to light gray only */
    #variants_table tbody tr:nth-child(odd) {
        background-color: #f5f5f5 !important;
    }

    /* Keep chart data tables visible, hide canvas */
    canvas {
        display: none !important;
    }
    .chart-data-table {
        position: static !important;
        width: auto !important;
        height: auto !important;
        overflow: visible !important;
        clip: auto !important;
    }
}
```

### Window.print() PDF Export Button
```html
<!-- Simple browser-based PDF export -->
<button id="export-pdf-btn" class="btn btn-primary">
    <svg aria-hidden="true" class="icon" width="16" height="16">
        <use xlink:href="#icon-download"/>
    </svg>
    Download PDF
</button>

<script>
document.getElementById('export-pdf-btn').addEventListener('click', function() {
    // Browser's native print dialog includes "Save as PDF" option
    window.print();
});
</script>

<style>
/* Hide the PDF export button when printing */
@media print {
    #export-pdf-btn {
        display: none !important;
    }
}
</style>
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| WCAG 2.0 | WCAG 2.1 Level AA | 2018 (updated requirement 2026) | Added mobile accessibility, low vision requirements; U.S. government enforcement starts April 2026 |
| role="grid" for all interactive tables | role="table" for sortable/filterable tables | ~2020 (ARIA 1.2 clarification) | Grid pattern only for 2D navigation/editing; tables can be interactive without being grids |
| Emoji icons (🔗) | SVG icons with aria-hidden + sr-only text | 2018-2020 | Emoji inconsistent across platforms, not reliably hidden from screen readers, sizing issues |
| html2pdf.js / jsPDF | window.print() with @media print | 2022+ (mobile memory limits) | JavaScript PDF libs create raster images, fail on mobile, large files; browser print preserves accessibility |
| Manual ARIA testing only | Automated (axe/WAVE/Lighthouse) + manual | 2019+ (tooling maturity) | Automated tools now catch 40% of issues reliably, but manual testing still required for other 60% |
| display: none for skip links | Off-screen positioning | 2015+ (WCAG 2.0 compliance) | display: none removes from tab order; off-screen keeps in tab order but hidden |

**Deprecated/outdated:**
- **WCAG 2.0 as sole standard:** WCAG 2.1 is now the requirement (Level AA for government sites as of April 2026)
- **Using role="application":** Disables all screen reader shortcuts; almost never appropriate for web content
- **Font icon libraries (Font Awesome 4):** SVG sprites are preferred (inline, cacheable, accessible, no HTTP request)
- **JavaScript-only tooltips without aria-describedby:** Tippy.js now handles this automatically for non-interactive tooltips

## Open Questions

Things that couldn't be fully resolved:

1. **Exact badge color values for Phase 17**
   - What we know: Template uses inline styles for badges (lines 669-747), some colors likely fail WCAG AA
   - What's unclear: Which specific badge colors fail contrast requirements (need to test each)
   - Recommendation: Task 1 should audit all badge colors with WebAIM Contrast Checker, create color constant map with compliant values

2. **Detail panel accessibility in print mode**
   - What we know: DataTables child rows contain expandable detail panels with grouped variant fields
   - What's unclear: Should print mode expand all detail panels, collapse all, or allow user choice? Could make printed reports very long.
   - Recommendation: Default to collapsed (main table only) in print, add optional "Expand All for Print" button that expands all rows before calling window.print()

3. **Screen reader verbosity for metric cards**
   - What we know: Dashboard has metric cards with statistics and gene lists
   - What's unclear: Best ARIA structure for card grids - should each card be role="region" with aria-label, or is semantic HTML sufficient?
   - Recommendation: Use semantic HTML (`<section>` with `<h3>` headings) first, test with screen reader, add role="region" + aria-labelledby only if needed

4. **DataTables FixedColumns ARIA implications**
   - What we know: FixedColumns clones DOM elements, creating two copies of the first column
   - What's unclear: Do screen readers announce the cloned column twice? Should cloned column have aria-hidden="true"?
   - Recommendation: Test with screen reader (NVDA/JAWS). If duplicate announcement occurs, add aria-hidden="true" to `.dtfc-fixed-left` clone after initialization

## Sources

### Primary (HIGH confidence)
- [WCAG 2.1 Official Specification](https://www.w3.org/TR/WCAG21/) - Contrast ratios, keyboard accessibility, bypass blocks
- [Chart.js Accessibility Documentation](https://www.chartjs.org/docs/latest/general/accessibility.html) - Canvas aria-label and fallback content
- [Tippy.js Accessibility Documentation](https://atomiks.github.io/tippyjs/v5/accessibility/) - Keyboard triggers, aria-describedby, tabindex recommendations
- [WebAIM: Skip Navigation Links](https://webaim.org/techniques/skipnav/) - Implementation pattern and CSS
- [WebAIM: Contrast Checker](https://webaim.org/resources/contrastchecker/) - Tool for verifying WCAG AA compliance
- [MDN: ARIA table role](https://developer.mozilla.org/en-US/docs/Web/Accessibility/ARIA/Reference/Roles/table_role) - When to use table vs grid
- [W3C ARIA Authoring Practices: Table Pattern](https://www.w3.org/WAI/ARIA/apg/patterns/table/) - ARIA structure for tables
- [W3C ARIA Authoring Practices: Grid Pattern](https://www.w3.org/WAI/ARIA/apg/patterns/grid/) - When grid is appropriate

### Secondary (MEDIUM confidence)
- [WCAG 2.1 AA Compliance Checklist 2026](https://www.webability.io/blog/wcag-2-1-aa-the-standard-for-accessible-web-design) - Implementation guide
- [DataTables ARIA Discussion](https://datatables.net/forums/discussion/23920/datatables-and-aria) - Community patterns for ARIA
- [CSS-Tricks: Accessible SVG Icons](https://css-tricks.com/accessible-svg-icons/) - SVG + sr-only pattern
- [Scott O'Hara: Contextual Images and SVGs](https://www.scottohara.me/blog/2019/05/22/contextual-images-svgs-and-a11y.html) - SVG accessibility patterns
- [MDN: Printing Guide](https://developer.mozilla.org/en-US/docs/Web/CSS/Guides/Media_queries/Printing) - @media print best practices
- [CSS-Tricks: Print Stylesheet Approaches](https://css-tricks.com/print-stylesheet-approaches-blacklist-vs-whitelist/) - Hiding interactive elements
- [Adrian Roselli: ARIA Grid as Anti-Pattern](https://adrianroselli.com/2020/07/aria-grid-as-an-anti-pattern.html) - Why not to use grid for tables

### Secondary (MEDIUM confidence - Testing Tools)
- [axe DevTools vs WAVE vs Lighthouse Comparison](https://inclly.com/resources/accessibility-testing-tools-comparison) - Tool capabilities and limitations
- [Automated Accessibility Testing Guide](https://a11y.coffee/start-testing/) - Testing workflow
- [WebAIM: Color Contrast and Accessibility](https://webaim.org/articles/contrast/) - Understanding contrast requirements

### Tertiary (LOW confidence - Community best practices)
- [Print Page Breaks for Tables - 2026 Best Practices](https://copyprogramming.com/howto/avoid-page-break-inside-row-of-table) - CSS fragmentation patterns
- [html2pdf.js Documentation](https://ekoopmans.github.io/html2pdf.js/) - JavaScript PDF library (NOT recommended for this use case)
- [How to Generate PDF with jsPDF](https://pdfnoodle.com/blog/generating-pdfs-from-html-with-jspdf) - Alternative PDF library (NOT recommended)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - WCAG 2.1 is official W3C standard, browser APIs are stable, no library dependencies
- Architecture: HIGH - Patterns from official documentation (WCAG, Chart.js, Tippy.js, WebAIM), widely used
- Pitfalls: MEDIUM - Based on community forums and blog posts, validated against official docs where possible
- Color contrast specifics: LOW - Need to audit actual badge colors in template with contrast checker tool

**Research date:** 2026-02-17
**Valid until:** 30 days (WCAG 2.1 stable, browser APIs stable, tools mature)

**Note on DataTables v2 ARIA:** DataTables v2 has improved ARIA for pagination controls but does NOT automatically add table structure roles. These must be added manually after initialization. Source: [DataTables forums](https://datatables.net/forums/discussion/23920/datatables-and-aria)
