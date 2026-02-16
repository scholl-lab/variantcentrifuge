# Phase 13: JS Stack Modernization - Research

**Researched:** 2026-02-16
**Domain:** Frontend JavaScript modernization for single-file HTML reports
**Confidence:** HIGH

## Summary

Phase 13 replaces the current jQuery + DataTables v1.10 + Plotly stack with vanilla JS + DataTables v2 + Chart.js + Tippy.js for the individual HTML report. The migration must preserve all existing functionality (sorting, pagination, search, column visibility, horizontal scroll, hover-expand cells) while modernizing the stack and reducing bundle size from ~3.5MB (Plotly) to ~65KB (Chart.js).

The report architecture requires **complete self-containment**: all JavaScript, CSS, and data must be embedded directly in a single HTML file for offline use in clinical/HPC environments where internet access is restricted. This constraint eliminates CDN links and requires version-locked minified library files shipped with the Python package.

**Primary recommendation:** Use DataTables v2 with vanilla JS API (jQuery still required as dependency but no jQuery code written), Chart.js v4 with chartjs-plugin-datalabels for on-bar labels, and Tippy.js v6 for accessible tooltips. Inline all assets using Jinja2 template with `<script>` and `<style>` tags. Package minified files in `variantcentrifuge/assets/` using hatchling's default file inclusion.

## Standard Stack

The established libraries for this migration:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| DataTables | 2.x (latest: 2.1+) | Interactive table with search/sort/pagination | Industry standard for feature-rich HTML tables; v2 provides vanilla JS API while maintaining jQuery under the hood |
| Chart.js | 4.x (latest: 4.5.1) | Lightweight charting library | Dominant in JavaScript charting (smaller bundle than Plotly, canvas-based for performance) |
| Tippy.js | 6.x (latest: 6.3.7) | Accessible tooltips with viewport awareness | Powered by Popper.js, handles positioning edge cases automatically |
| chartjs-plugin-datalabels | 2.x | Display value labels on Chart.js bars | Official Chart.js plugin for on-bar text labels |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Jinja2 | (already in deps) | Template rendering | Generate single-file HTML with inlined assets |
| hatchling | (already build system) | Package static assets | Include minified JS/CSS files in Python wheel |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| DataTables v2 | simple-datatables (vanilla JS) | Would avoid jQuery dependency entirely but lose feature parity (buttons extension, mature ecosystem) and require rewriting all table logic |
| Chart.js | Plotly.js | Current implementation; 54x larger bundle (~3.5MB vs ~65KB), overkill for simple horizontal bar chart |
| Tippy.js | Custom vanilla JS tooltip | User's CONTEXT.md explicitly chose Tippy.js to avoid temporary solutions |

**Installation:**
```bash
# Download minified files manually and place in variantcentrifuge/assets/
# DataTables v2: dataTables.js, dataTables.dataTables.css
# DataTables Buttons: dataTables.buttons.js, buttons.dataTables.css, buttons.colVis.js
# Chart.js: chart.umd.js (or chart.min.js)
# chartjs-plugin-datalabels: chartjs-plugin-datalabels.min.js
# Tippy.js: tippy-bundle.umd.min.js (includes Popper)
```

**Bundle Size Reference:**
- Chart.js v4: ~65-70KB minified (estimate based on ecosystem comparisons)
- Tippy.js v6.3.7: ~15-20KB minified with Popper bundled
- DataTables v2: ~90-100KB core + ~20KB Buttons extension
- **Total new stack:** ~200KB vs. **Plotly alone:** ~3.5MB (94% reduction)

## Architecture Patterns

### Recommended Project Structure
```
variantcentrifuge/
├── assets/                          # NEW: Static library files
│   ├── js/
│   │   ├── datatables.min.js
│   │   ├── datatables.buttons.min.js
│   │   ├── buttons.colVis.min.js
│   │   ├── chart.umd.min.js
│   │   ├── chartjs-plugin-datalabels.min.js
│   │   └── tippy-bundle.umd.min.js
│   └── css/
│       ├── datatables.min.css
│       └── buttons.dataTables.min.css
├── templates/
│   └── index.html                   # MODIFIED: Inline assets, modern stack
└── generate_html_report.py          # MODIFIED: Load assets, pass to template
```

### Pattern 1: Inline Asset Embedding
**What:** Embed minified JS/CSS directly in HTML template using Jinja2 variables
**When to use:** Single-file HTML reports for offline environments

**Example:**
```python
# generate_html_report.py
def generate_html_report(...):
    assets_dir = Path(__file__).parent / "assets"

    # Load all library files
    with open(assets_dir / "js" / "datatables.min.js") as f:
        datatables_js = f.read()
    with open(assets_dir / "js" / "chart.umd.min.js") as f:
        chartjs_js = f.read()
    # ... load other assets

    html_content = template.render(
        variants=variants_data,
        summary=summary,
        # Pass assets to template
        datatables_js=datatables_js,
        chartjs_js=chartjs_js,
        # ... other assets
    )
```

```html
<!-- templates/index.html -->
<head>
    <style>
        {{ datatables_css }}
        {{ buttons_css }}
        /* Custom styles below libraries */
    </style>
</head>
<body>
    <!-- Content -->
    <script>{{ datatables_js }}</script>
    <script>{{ datatables_buttons_js }}</script>
    <script>{{ chartjs_js }}</script>
    <script>{{ chartjs_datalabels_js }}</script>
    <script>{{ tippy_js }}</script>
    <script>
        // Custom initialization code
    </script>
</body>
```

### Pattern 2: DataTables v2 Vanilla JS Initialization
**What:** Use vanilla JS API for DataTables (no jQuery code) while jQuery runs under the hood
**When to use:** All DataTable initialization in v2

**Example:**
```javascript
// Source: https://datatables.net/blog/2024/datatables-2
// Vanilla JS initialization (DataTables v2)
let table = new DataTable('#variants_table', {
    autoWidth: false,
    pageLength: 10,
    scrollX: true,
    layout: {
        topStart: 'pageLength',
        topEnd: 'buttons',
        bottomStart: 'info',
        bottomEnd: 'paging'
    },
    lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
    buttons: [
        {
            extend: 'colvis',
            text: 'Show/Hide Columns'
        }
    ],
    columnDefs: [
        { targets: hiddenColumnIndices, visible: false }
    ],
    drawCallback: function(settings) {
        this.api().columns.adjust();
        // Check for truncated cells, add Tippy tooltips
    }
});
```

### Pattern 3: Chart.js Horizontal Bar with On-Bar Labels
**What:** Chart.js v4 horizontal bar chart with datalabels plugin
**When to use:** Impact distribution chart

**Example:**
```javascript
// Source: https://www.chartjs.org/docs/latest/charts/bar.html
// Horizontal bar chart with on-bar labels
const ctx = document.getElementById('impact_chart').getContext('2d');
new Chart(ctx, {
    type: 'bar',
    data: {
        labels: ['HIGH', 'MODERATE', 'LOW', 'MODIFIER'],
        datasets: [{
            data: [10, 45, 23, 5],
            backgroundColor: [
                'rgba(220, 53, 69, 0.8)',   // HIGH: red
                'rgba(255, 159, 64, 0.8)',  // MODERATE: orange
                'rgba(255, 205, 86, 0.8)',  // LOW: amber
                'rgba(156, 163, 175, 0.8)'  // MODIFIER: gray
            ]
        }]
    },
    options: {
        indexAxis: 'y',  // Horizontal bars
        responsive: true,
        plugins: {
            legend: { display: false },
            tooltip: {
                callbacks: {
                    label: function(context) {
                        return 'Count: ' + context.parsed.x;
                    }
                }
            },
            datalabels: {
                anchor: 'end',
                align: 'end',
                color: '#333',
                font: { weight: 'bold' },
                formatter: (value) => value
            }
        }
    },
    plugins: [ChartDataLabels]
});
```

### Pattern 4: Tippy.js Tooltip Initialization
**What:** Replace hover-expand cells with Tippy.js tooltips
**When to use:** Truncated table cells (HGVS_C, HGVS_P, GT, etc.)

**Example:**
```javascript
// Source: https://atomiks.github.io/tippyjs/
// Initialize Tippy.js for truncated cells
table.on('draw', function() {
    document.querySelectorAll('.truncated-cell').forEach(function(cell) {
        if (cell.scrollWidth > cell.offsetWidth) {
            tippy(cell, {
                content: cell.textContent,
                placement: 'top',
                theme: 'dark',
                maxWidth: 400,
                interactive: false
            });
        }
    });
});
```

### Pattern 5: Skeleton Screen with CSS Animation
**What:** Show loading skeleton while DataTables initializes
**When to use:** Table-heavy pages with noticeable initialization delay

**Example:**
```html
<!-- Source: https://www.freecodecamp.org/news/how-to-build-skeleton-screens-using-css-for-better-user-experience/ -->
<div id="skeleton" class="skeleton-wrapper">
    <div class="skeleton-row"></div>
    <div class="skeleton-row"></div>
    <div class="skeleton-row"></div>
    <div class="skeleton-row"></div>
    <div class="skeleton-row"></div>
</div>
<table id="variants_table" style="display:none;">...</table>

<style>
.skeleton-row {
    height: 30px;
    background: linear-gradient(90deg, #f0f0f0 25%, #e0e0e0 50%, #f0f0f0 75%);
    background-size: 200% 100%;
    animation: shimmer 1.5s infinite;
    margin-bottom: 8px;
    border-radius: 4px;
}

@keyframes shimmer {
    0% { background-position: -200% 0; }
    100% { background-position: 200% 0; }
}
</style>

<script>
let table = new DataTable('#variants_table', {
    initComplete: function() {
        document.getElementById('skeleton').style.display = 'none';
        document.getElementById('variants_table').style.display = 'table';
    }
});
</script>
```

### Anti-Patterns to Avoid

- **Using jQuery syntax in custom code:** DataTables v2 still requires jQuery as a dependency, but you should write vanilla JS (`document.querySelector()` not `$()`) to avoid spreading jQuery usage
- **CDN links in template:** Reports must work offline; CDN links break in HPC/clinical environments
- **Mixing old/new DataTables APIs:** Don't mix `$('#table').DataTable()` with `new DataTable()` — choose vanilla JS API consistently
- **Inlining unminified libraries:** Bundle size matters; always use minified versions
- **Custom tooltip positioning logic:** Tippy.js handles viewport awareness automatically; don't reimplement

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Tooltip viewport positioning | Custom absolute positioning with boundary detection | Tippy.js | Edge cases are complex (viewport edges, scrolling, parent transforms, RTL languages). Tippy + Popper handle all of this. |
| Table column visibility UI | Custom dropdown with checkboxes | DataTables Buttons colvis extension | Already built into DataTables v2, handles show/hide state, accessibility, and responsive behavior |
| Skeleton loading animation | Custom div manipulation and CSS | CSS keyframe animation with `background-position` | Modern browsers optimize CSS animations better than JS manipulation; shimmer effect is now a CSS pattern |
| Chart responsiveness | Manual canvas resize listeners | Chart.js `responsive: true` option | Chart.js handles container resize observers, debouncing, and redraw optimization |
| XSS protection in tooltips | Manual HTML escaping | Tippy.js `allowHTML: false` (default) | Tippy.js escapes content by default; only enable HTML when content is trusted |

**Key insight:** Single-file HTML reports seem simple but have hard constraints (offline, no external deps, security). Using battle-tested libraries prevents edge case bugs and security issues.

## Common Pitfalls

### Pitfall 1: DataTables v2 File Naming
**What goes wrong:** Importing `jquery.dataTables.js` (v1 filename) fails silently or loads wrong version
**Why it happens:** DataTables v2 renamed core file from `jquery.dataTables.js` to `dataTables.js`
**How to avoid:** Always use `dataTables.js` (v2 filename) and `dataTables.dataTables.css` for default styling
**Warning signs:** DataTables API missing methods, layout option not recognized, console errors about undefined properties

### Pitfall 2: deferRender Default Change
**What goes wrong:** `table.row().node()` returns null for non-visible rows in v2
**Why it happens:** DataTables v2 changed `deferRender` default from `false` to `true` (rows created only when visible)
**How to avoid:** If code relies on accessing non-visible row DOM nodes, set `deferRender: false` explicitly in config
**Warning signs:** TypeError on `.node()` calls, Tippy.js initialization failing for off-page rows

### Pitfall 3: Chart.js Plugin Registration
**What goes wrong:** chartjs-plugin-datalabels doesn't display labels on bars
**Why it happens:** Plugin must be explicitly registered in Chart.js v4 (not auto-registered)
**How to avoid:** Pass plugin in chart config: `plugins: [ChartDataLabels]` and ensure plugin script loads before chart initialization
**Warning signs:** Chart renders but no labels appear, no console errors (silent failure)

### Pitfall 4: Tippy.js on Dynamically Created Content
**What goes wrong:** Tooltips don't appear on table cells after pagination/search
**Why it happens:** Tippy.js initialized once on page load, doesn't know about new DOM elements created by DataTables
**How to avoid:** Initialize Tippy in DataTables `drawCallback` (fires after every table redraw)
**Warning signs:** Tooltips work on page 1 but not on page 2+, tooltips disappear after filtering

### Pitfall 5: Inline Script CSP Violations
**What goes wrong:** Browser blocks inline scripts with strict Content-Security-Policy headers
**Why it happens:** Single-file HTML requires inline scripts, conflicts with CSP `script-src` directives
**How to avoid:** Reports are static files (not served by web server), but if CSP becomes an issue, use CSP nonce or hash in template
**Warning signs:** Console errors "Refused to execute inline script", blank page in high-security environments

### Pitfall 6: jQuery Dependency Confusion
**What goes wrong:** Developer assumes "vanilla JS API" means no jQuery needed, removes jQuery script tag
**Why it happens:** DataTables v2 marketing emphasizes vanilla JS API but still requires jQuery internally
**How to avoid:** DataTables v2 **requires jQuery as a dependency** even when using vanilla JS syntax. jQuery must be loaded but you don't write jQuery code.
**Warning signs:** ReferenceError: $ is not defined, DataTable is not defined, table initialization fails

### Pitfall 7: Asset Loading Order
**What goes wrong:** TypeError: DataTable.Buttons is undefined, ChartDataLabels is not defined
**Why it happens:** Extensions/plugins loaded before their base libraries
**How to avoid:** Strict loading order: jQuery → DataTables core → DataTables Buttons → Chart.js → chartjs-plugin-datalabels → Tippy.js
**Warning signs:** Initialization errors referencing undefined constructors, console shows script loaded but object undefined

## Code Examples

Verified patterns from official sources:

### jQuery Removal Common Patterns
```javascript
// Source: https://tobiasahlin.com/blog/move-from-jquery-to-vanilla-javascript/

// jQuery → Vanilla JS equivalents
// DOM Selection
$('.class')          → document.querySelectorAll('.class')
$('#id')             → document.getElementById('id')

// Event Handling
$('.btn').click(fn)  → document.querySelectorAll('.btn').forEach(btn =>
                          btn.addEventListener('click', fn))

// Class Manipulation
$el.addClass('x')    → el.classList.add('x')
$el.removeClass('x') → el.classList.remove('x')
$el.toggleClass('x') → el.classList.toggle('x')

// Element Manipulation
$el.show()           → el.style.display = 'block'
$el.hide()           → el.style.display = 'none'

// Attribute Access
$el.attr('data-x')   → el.getAttribute('data-x')
$el.data('x')        → el.dataset.x
```

### DataTables v2 Layout vs. DOM
```javascript
// Source: https://datatables.net/upgrade/2

// Old (v1): dom string
$('#table').DataTable({
    dom: 'lBfrtip'  // l=length, B=buttons, f=filter, r=processing, t=table, i=info, p=pagination
});

// New (v2): layout object (recommended)
new DataTable('#table', {
    layout: {
        topStart: 'pageLength',
        topEnd: 'buttons',
        bottomStart: 'info',
        bottomEnd: 'paging'
    }
});
```

### Packaging Static Assets with Hatchling
```python
# Source: https://packaging.python.org/tutorials/packaging-projects/

# pyproject.toml (current config uses hatchling)
[tool.hatch.build.targets.wheel]
packages = ["variantcentrifuge"]
# By default, hatchling includes all files in packages
# assets/ directory will be included if inside variantcentrifuge/

# generate_html_report.py
from pathlib import Path

def generate_html_report(...):
    # Load assets from package
    assets_dir = Path(__file__).parent / "assets"

    # Read minified library files
    js_files = {
        'datatables': (assets_dir / "js" / "datatables.min.js").read_text(encoding='utf-8'),
        'chartjs': (assets_dir / "js" / "chart.umd.min.js").read_text(encoding='utf-8'),
        'tippy': (assets_dir / "js" / "tippy-bundle.umd.min.js").read_text(encoding='utf-8'),
    }

    css_files = {
        'datatables': (assets_dir / "css" / "datatables.min.css").read_text(encoding='utf-8'),
    }

    # Pass to template
    html_content = template.render(
        variants=variants_data,
        summary=summary,
        js_libraries=js_files,
        css_libraries=css_files,
    )
```

### HTML Template Structure
```html
<!-- templates/index.html -->
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>VariantCentrifuge Report</title>
    <style>
        /* Inline library CSS */
        {{ css_libraries.datatables }}

        /* Custom styles */
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }
        .skeleton-row { /* ... */ }
    </style>
</head>
<body>
    <div id="skeleton"><!-- Skeleton rows --></div>
    <table id="variants_table" style="display:none;"><!-- Table content --></table>

    <!-- Inline library JS (order matters!) -->
    <script>{{ js_libraries.datatables }}</script>
    <script>{{ js_libraries.chartjs }}</script>
    <script>{{ js_libraries.tippy }}</script>

    <!-- Custom initialization -->
    <script>
        // Initialize DataTables
        let table = new DataTable('#variants_table', {
            // config
            initComplete: function() {
                document.getElementById('skeleton').style.display = 'none';
                document.getElementById('variants_table').style.display = 'table';
            }
        });

        // Initialize Chart.js
        const ctx = document.getElementById('impact_chart').getContext('2d');
        new Chart(ctx, { /* config */ });
    </script>
</body>
</html>
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| DataTables v1.10 (2015) | DataTables v2 (2024) | January 2024 | Vanilla JS API added, layout system replaces dom, class naming standardized with dt- prefix |
| jQuery `.DataTable()` syntax | `new DataTable()` vanilla JS syntax | DataTables v2 (2024) | Can write vanilla JS code but jQuery still required under the hood |
| Plotly.js for charts | Chart.js for simple charts | Ongoing trend (2020+) | Chart.js dominant for basic charts (smaller bundle, faster), Plotly used for complex scientific visualizations |
| Chart.js horizontal bar type | Chart.js `indexAxis: 'y'` | Chart.js v3 (2021) | Removed separate horizontalBar type, unified with bar + indexAxis option |
| Custom tooltip positioning | Tippy.js/Popper-based solutions | 2018+ | Viewport-aware tooltips now standard, manual positioning considered legacy |
| CDN-linked reports | Self-contained single-file HTML | Always for clinical/HPC | Offline requirement never changed, but CDN pattern common in web dev |

**Deprecated/outdated:**
- **DataTables ColVis extension (separate):** Merged into Buttons extension as colvis button type (v1.6+, 2019)
- **jQuery dependency for simple DOM tasks:** Vanilla JS now has feature parity (querySelector, classList, fetch) since IE11 retirement
- **Chart.js v2 horizontalBar type:** Removed in v3 (2021), use `type: 'bar', options: { indexAxis: 'y' }`
- **Inline event handlers (onclick=""):** Considered anti-pattern for CSP and maintainability; use addEventListener

## Open Questions

Things that couldn't be fully resolved:

1. **Chart.js exact bundle size**
   - What we know: Chart.js v4 is significantly lighter than Plotly (~65KB estimate vs ~3.5MB confirmed), multiple sources confirm Chart.js as "lightweight"
   - What's unclear: Exact minified+gzip size for Chart.js v4.5.1 not found in search results (Bundlephobia page didn't display numbers)
   - Recommendation: Download chart.umd.min.js from official release and measure actual file size before finalizing. Likely 60-80KB minified.

2. **DataTables v2 jQuery removal timeline**
   - What we know: DataTables developers stated they lack resources to build a jQuery-free version (confirmed in forum discussions 2024)
   - What's unclear: Whether v3 or future versions will drop jQuery dependency
   - Recommendation: Accept jQuery dependency for now. If it becomes a blocker, consider migrating to simple-datatables (vanilla JS fork) but that's Phase 14+ scope.

3. **CSP nonce/hash for inline scripts**
   - What we know: Single-file HTML requires inline scripts which conflict with strict CSP policies
   - What's unclear: Whether clinical/HPC environments serve reports through web servers with CSP headers (likely not — files opened directly in browser)
   - Recommendation: Implement basic version without CSP nonces. If environment requires CSP, add nonce generation in Phase 14+ as enhancement.

4. **Tippy.js bundle variants**
   - What we know: Tippy.js v6 offers multiple bundles (with/without Popper, UMD vs ESM)
   - What's unclear: Which bundle variant is smallest for single-file inlining (tippy-bundle.umd.min.js vs alternatives)
   - Recommendation: Use tippy-bundle.umd.min.js (includes Popper, UMD format for browser globals). Verify bundle includes all dependencies by testing tooltip positioning edge cases.

## Sources

### Primary (HIGH confidence)
- [DataTables v2 Upgrade Guide](https://datatables.net/upgrade/2) - Migration guide, breaking changes, file naming
- [DataTables v2 Blog Announcement](https://datatables.net/blog/2024/datatables-2) - New features, layout system, vanilla JS API
- [Chart.js Official Documentation](https://www.chartjs.org/docs/latest/) - Getting started, configuration, chart types
- [Chart.js Bar Chart Documentation](https://www.chartjs.org/docs/latest/charts/bar.html) - Horizontal bar configuration, indexAxis option
- [Tippy.js Official Documentation](https://atomiks.github.io/tippyjs/) - Installation, configuration, accessibility features
- [Tippy.js Accessibility Guide](https://atomiks.github.io/tippyjs/v5/accessibility/) - ARIA support, focus handling
- [chartjs-plugin-datalabels Documentation](https://chartjs-plugin-datalabels.netlify.app/) - On-bar label configuration

### Secondary (MEDIUM confidence)
- [jQuery vs Vanilla JS Migration Cheat Sheet](https://tobiasahlin.com/blog/move-from-jquery-to-vanilla-javascript/) - Common pattern conversions
- [DataTables Column Visibility Buttons](https://datatables.net/extensions/buttons/examples/column_visibility/index.html) - colvis button examples
- [Skeleton Loading CSS Tutorial](https://www.freecodecamp.org/news/how-to-build-skeleton-screens-using-css-for-better-user-experience/) - Shimmer animation patterns
- [Python Packaging Static Assets](https://packaging.python.org/tutorials/packaging-projects/) - Include non-code files in packages
- [OWASP XSS Prevention Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Cross_Site_Scripting_Prevention_Cheat_Sheet.html) - Security best practices for HTML reports

### Tertiary (LOW confidence)
- [Luzmo JavaScript Chart Library Comparison](https://www.luzmo.com/blog/best-javascript-chart-libraries) - Chart.js vs Plotly bundle size discussion
- [LogRocket Charting Libraries Comparison](https://blog.logrocket.com/comparing-most-popular-javascript-charting-libraries/) - Feature comparison, marked LOW due to lack of specific version details
- [DataTables Forum: jQuery 2.0 Discussion](https://datatables.net/forums/discussion/67258/will-datatables-2-0-require-jquery) - Community discussion on jQuery dependency (not official roadmap)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official documentation confirmed DataTables v2, Chart.js v4, Tippy.js v6 current versions and capabilities
- Architecture: HIGH - Jinja2 template pattern verified, hatchling packaging confirmed in pyproject.toml
- Pitfalls: HIGH - Breaking changes documented in official upgrade guides, common mistakes catalogued in forums

**Research date:** 2026-02-16
**Valid until:** 2026-03-16 (30 days - stack is stable, DataTables v2 released 2024, Chart.js v4 released 2022)
