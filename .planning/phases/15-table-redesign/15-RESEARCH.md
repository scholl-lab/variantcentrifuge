# Phase 15: Table Redesign - Research

**Researched:** 2026-02-17
**Domain:** DataTables v2 table UX enhancement with tooltips, sticky columns, expandable details, and density controls
**Confidence:** HIGH

## Summary

Phase 15 transforms the variant data table into a modern enterprise-grade interface. The current implementation (post-Phase 13 and 14) uses DataTables v2.2.2 with vanilla JS, Tippy.js tooltips on truncated cells, and semantic color badges. This phase adds:

1. **Sticky first column (GENE)** using DataTables FixedColumns extension
2. **Expandable row details** via DataTables child row API with custom renderer
3. **Content density toggle** (Compact/Regular/Relaxed) with localStorage persistence
4. **Enhanced styling** including zebra striping, dark header, intelligent column widths
5. **Improved truncation** with middle truncation for variant IDs, monospace fonts, right-aligned numbers

The standard approach is to extend the existing DataTables initialization with additional extensions (FixedColumns), custom render functions for column widths, child row handlers, and a density control UI.

**Primary recommendation:** Vendor DataTables FixedColumns extension (JS + CSS), use DataTables child row API (`row.child()`) for expandable details, implement density modes via CSS classes controlled by localStorage, and enhance column configuration with `columnDefs` width/alignment/truncation rules.

## Standard Stack

### Core Extensions (Required)

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| datatables.net-fixedcolumns | 5.0.x | Sticky columns with horizontal scroll | Official DataTables extension, maintained by SpryMedia, battle-tested for freeze columns |
| datatables.net-fixedcolumns-dt | 5.0.x | FixedColumns styling for DataTables default theme | Matching CSS for the dt (default) theme used in current template |

### Already Vendored (Phase 13)

| Library | Version | Purpose | Status |
|---------|---------|---------|--------|
| datatables.net | 2.2.2 | Core DataTables v2 | ✓ Already vendored (95KB) |
| tippy.js | 6.3.7 | Tooltip library | ✓ Already vendored (26KB bundle with Popper) |
| Chart.js | 4.4.8 | Dashboard charts | ✓ Already vendored (202KB) |
| jQuery slim | 3.7.1 | DataTables dependency | ✓ Already vendored (70KB) |

### Installation

```bash
# Download FixedColumns from jsDelivr CDN
curl -o variantcentrifuge/assets/js/fixedcolumns.min.js \
  https://cdn.jsdelivr.net/npm/datatables.net-fixedcolumns@5.0.5/js/dataTables.fixedColumns.min.js

curl -o variantcentrifuge/assets/css/fixedcolumns.dataTables.min.css \
  https://cdn.jsdelivr.net/npm/datatables.net-fixedcolumns-dt@5.0.5/css/fixedColumns.dataTables.min.css
```

**Note:** Latest version is 5.0.5 (as of Feb 2026). Verify compatibility with DataTables 2.2.2 before vendoring.

## Architecture Patterns

### Current Template Structure (Post-Phase 13/14)

```
variantcentrifuge/templates/index.html
├── <head>
│   ├── Inlined CSS assets (datatables.min.css, buttons, tippy.css)
│   └── Custom <style> block (body, header, table, badges, tooltips)
├── <body>
│   ├── <header> — Report title and metadata
│   ├── <div class="dashboard"> — Metric cards + Chart.js charts
│   ├── <div class="table-section">
│   │   ├── <div id="table-skeleton"> — Loading shimmer
│   │   └── <table id="variants_table"> — Main variant table
│   └── <footer class="metadata-footer"> — Pipeline metadata
└── <script>
    ├── Inlined JS assets (jquery, datatables, buttons, chart, tippy)
    └── DOMContentLoaded handler
        ├── DataTable initialization (new DataTable('#variants_table', {...}))
        ├── Badge render functions (IMPACT, ClinVar, Inheritance_Pattern)
        ├── Tippy.js tooltip attachment (initTippyTooltips())
        └── Chart.js initialization (impact chart, inheritance chart)
```

### Current DataTables Configuration (Line 606-634)

```javascript
var table = new DataTable('#variants_table', {
    autoWidth: false,
    pageLength: 25,
    scrollX: true,
    layout: {
        topStart: {
            pageLength: { menu: [10, 25, 50, 100, -1] }
        },
        topEnd: {
            search: true,
            buttons: [{ extend: 'colvis', text: 'Show/Hide Columns' }]
        }
    },
    columnDefs: dtColumnDefs,  // Hidden columns + badge render functions
    initComplete: function() {
        skeleton.style.display = 'none';
        this.api().columns.adjust();
    },
    drawCallback: function(settings) {
        initTippyTooltips();  // Re-attach tooltips after pagination/sort
    }
});
```

**Key observations:**
- `scrollX: true` enables horizontal scrolling (required for FixedColumns)
- `layout` API is DataTables v2 syntax (replaces v1's `dom` option)
- `columnDefs` array already used for hidden columns and badge rendering
- `drawCallback` re-initializes Tippy tooltips on table redraw

### Pattern 1: Sticky Column with FixedColumns

**What:** Freeze the GENE column (or any left column) so it remains visible during horizontal scroll

**When to use:** Tables with many columns (>10) where horizontal scroll is necessary, and key identifier column (GENE) should always be visible

**Example:**
```javascript
// Source: https://datatables.net/extensions/fixedcolumns/examples/initialisation/left_right_columns.html
var table = new DataTable('#variants_table', {
    scrollX: true,           // Required for FixedColumns
    fixedColumns: {
        left: 1,             // Freeze first column (GENE)
        // Optional shadow separator via CSS (not JS config)
    },
    // ... other config
});
```

**CSS for shadow separator (UX best practice):**
```css
/* Phase 15 addition: Drop shadow on frozen column for "floating" effect */
.dtfc-fixed-left {
    box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);  /* Right shadow on frozen column */
}
```

**Note:** DataTables FixedColumns works by cloning the fixed columns into separate `<div>` containers positioned absolutely. The original columns remain in the scrolling table but are visually hidden via CSS.

### Pattern 2: Expandable Row Details (Child Rows)

**What:** Click a chevron/icon to expand a row and show all variant fields in a grouped layout

**When to use:** Tables with many columns where not all data fits on screen, providing "progressive disclosure" UX pattern

**Example:**
```javascript
// Source: https://datatables.net/examples/api/row_details.html
// 1. Add control column to columnDefs
columnDefs: [
    {
        targets: 0,
        className: 'dt-control',
        orderable: false,
        data: null,
        defaultContent: '<span class="chevron">›</span>'  // Chevron icon
    },
    // ... other column defs
]

// 2. Define format function for child row content
function formatChildRow(rowData) {
    // Return HTML with grouped fields
    return `
        <div class="detail-panel">
            <div class="detail-section">
                <h4>Identifiers</h4>
                <dl>
                    <dt>VAR_ID</dt><dd>${rowData.VAR_ID}</dd>
                    <dt>ID</dt><dd>${rowData.ID || 'N/A'}</dd>
                </dl>
            </div>
            <div class="detail-section">
                <h4>Annotations</h4>
                <dl>
                    <dt>HGVS_C</dt><dd>${rowData.HGVS_C}</dd>
                    <dt>HGVS_P</dt><dd>${rowData.HGVS_P}</dd>
                </dl>
            </div>
            <!-- More sections: Scores, Links -->
        </div>
    `;
}

// 3. Attach click handler for child row toggle
table.on('click', 'tbody td.dt-control', function(e) {
    let tr = e.target.closest('tr');
    let row = table.row(tr);

    if (row.child.isShown()) {
        row.child.hide();
        tr.classList.remove('shown');
    } else {
        row.child(formatChildRow(row.data())).show();
        tr.classList.add('shown');
    }
});
```

**CSS for child row styling:**
```css
/* Chevron rotation animation */
tr.shown td.dt-control .chevron {
    transform: rotate(90deg);
    transition: transform 200ms ease-in-out;
}

/* Detail panel layout */
.detail-panel {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
    gap: 16px;
    padding: 16px;
    background-color: #f8f9fa;
    border-left: 3px solid #007bff;
}

.detail-section {
    background: white;
    padding: 12px;
    border-radius: 4px;
    border: 1px solid #dee2e6;
}

.detail-section h4 {
    margin: 0 0 8px 0;
    font-size: 13px;
    text-transform: uppercase;
    color: #6c757d;
}

.detail-section dl {
    display: grid;
    grid-template-columns: auto 1fr;
    gap: 4px 12px;
    margin: 0;
    font-size: 13px;
}

.detail-section dt {
    font-weight: 600;
    color: #495057;
}

.detail-section dd {
    margin: 0;
    word-break: break-word;
}
```

### Pattern 3: Content Density Modes

**What:** Toggle between Compact/Regular/Relaxed padding modes, persisted via localStorage

**When to use:** Tables used by power users who need to scan many rows (compact) vs. casual users who prefer whitespace (relaxed)

**Example:**
```javascript
// 1. Define density CSS classes
/* In <style> block */
#variants_table.density-compact td { padding: 4px 6px; }
#variants_table.density-regular td { padding: 8px 10px; }
#variants_table.density-relaxed td { padding: 12px 14px; }

// 2. Load density from localStorage on init
var savedDensity = localStorage.getItem('variantTableDensity') || 'compact';
document.getElementById('variants_table').classList.add('density-' + savedDensity);

// 3. Add density toggle to layout.topEnd
layout: {
    topEnd: {
        search: true,
        buttons: [
            {
                extend: 'colvis',
                text: 'Show/Hide Columns'
            },
            {
                text: 'Density: Compact',
                action: function(e, dt, node, config) {
                    cycleDensity(dt);
                }
            }
        ]
    }
}

// 4. Density cycle function
var densityModes = ['compact', 'regular', 'relaxed'];
var currentDensityIndex = densityModes.indexOf(savedDensity);

function cycleDensity(dt) {
    var tableNode = dt.table().node();

    // Remove old class
    tableNode.classList.remove('density-' + densityModes[currentDensityIndex]);

    // Cycle to next
    currentDensityIndex = (currentDensityIndex + 1) % densityModes.length;
    var newDensity = densityModes[currentDensityIndex];

    // Apply new class
    tableNode.classList.add('density-' + newDensity);

    // Save to localStorage
    localStorage.setItem('variantTableDensity', newDensity);

    // Update button text
    this.text('Density: ' + newDensity.charAt(0).toUpperCase() + newDensity.slice(1));

    // Redraw to adjust row heights
    dt.columns.adjust().draw();
}
```

### Pattern 4: Intelligent Column Widths

**What:** Apply type-specific column widths and text handling (fixed, grow, truncate, monospace, right-align)

**When to use:** Always — improves scannability by matching visual presentation to data semantics

**Example:**
```javascript
// In columnDefs array
columnDefs: [
    // Fixed width for genomic coordinates
    {
        targets: ['CHROM', 'POS'],  // Can use column names if columns.name is defined
        width: '80px',
        className: 'monospace'
    },
    // Grow for gene names
    {
        targets: 'GENE',
        width: 'auto',
        className: 'gene-column'
    },
    // Truncate with ellipsis for HGVS
    {
        targets: ['HGVS_C', 'HGVS_P'],
        width: '180px',
        className: 'monospace truncate',
        render: function(data, type, row) {
            if (type === 'display' && data && data.length > 30) {
                var truncated = data.substring(0, 30) + '...';
                return '<span class="truncated-cell" data-tippy-content="' + data + '">' + truncated + '</span>';
            }
            return data;
        }
    },
    // Right-align numeric scores
    {
        targets: ['dbNSFP_CADD_phred', 'dbNSFP_REVEL_score'],
        className: 'dt-type-numeric dt-right',
        render: function(data, type, row) {
            if (type === 'display' && data) {
                return parseFloat(data).toFixed(2);
            }
            return data;
        }
    },
    // Middle truncation for variant IDs
    {
        targets: 'VAR_ID',
        width: '120px',
        render: function(data, type, row) {
            if (type === 'display' && data && data.length > 20) {
                var start = data.substring(0, 10);
                var end = data.substring(data.length - 7);
                return '<span class="middle-truncate" data-tippy-content="' + data + '">' + start + '...' + end + '</span>';
            }
            return data;
        }
    }
]
```

**CSS for truncation classes:**
```css
/* Standard end truncation */
.truncate {
    max-width: 180px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}

/* Monospace for genomic notation */
.monospace {
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
}

/* Right-align numbers */
.dt-right {
    text-align: right !important;
}

/* Middle truncation (CSS can't do this, requires JS render) */
.middle-truncate {
    display: inline-block;
    max-width: 120px;
}
```

### Pattern 5: Enhanced Header Styling

**What:** Dark header with white text, colored bottom border, zebra striping on rows

**When to use:** Always — improves visual hierarchy and reduces eye strain on large tables

**Example:**
```css
/* Dark header row */
#variants_table thead th,
.dt-scroll-head thead th {
    background: linear-gradient(180deg, #2c3e50 0%, #1a252f 100%);
    color: white;
    border-bottom: 3px solid #007bff;  /* Accent color from phase 14 */
    font-weight: 600;
    padding: 12px 8px;
}

/* Zebra striping */
#variants_table tbody tr:nth-child(odd) {
    background-color: #f8f9fa;
}

#variants_table tbody tr:nth-child(even) {
    background-color: white;
}

/* Hover row highlight */
#variants_table tbody tr:hover {
    background-color: #e3f2fd;
}
```

### Anti-Patterns to Avoid

- **Don't use FixedColumns without scrollX** — FixedColumns requires horizontal scrolling to function
- **Don't modify fixed column clone manually** — DataTables manages cloned columns, manual DOM changes will be overwritten
- **Don't put interactive elements in child rows without event delegation** — Child rows are created dynamically, use delegated event handlers
- **Don't forget to call `columns.adjust()` after density change** — DataTables needs to recalculate column widths
- **Don't use CSS `position: sticky` instead of FixedColumns** — Breaks with DataTables' scrolling implementation

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Sticky columns with horizontal scroll | Custom position:sticky + scroll sync | DataTables FixedColumns extension | Handles edge cases: column visibility toggle, dynamic column widths, sorting, filtering, pagination, responsive behavior |
| Child row expand/collapse | Custom row insertion + animation | DataTables `row.child()` API | Maintains proper DOM structure, handles pagination edge cases, cleans up child rows automatically on redraw |
| Column width management | Manual CSS width + JS calculation | DataTables `columnDefs.width` | Auto-adjusts for content, handles hidden columns, responsive width on window resize |
| Tooltip positioning at viewport edges | Custom tooltip library or manual calculation | Tippy.js (already vendored) | Handles viewport clipping, flipping, arrow positioning, touch devices, keyboard accessibility |
| Table state persistence | Custom localStorage + JSON serialization | DataTables `stateSave` option | Saves page, sort, search, column visibility — but for density only, manual localStorage is simpler |

**Key insight:** DataTables has a mature extension ecosystem (FixedColumns, Responsive, RowGroup, etc.) with edge case handling. Don't reinvent sticky columns or child rows — use the official APIs.

## Common Pitfalls

### Pitfall 1: FixedColumns Not Appearing

**What goes wrong:** FixedColumns extension loaded but frozen column doesn't stick during scroll

**Why it happens:**
- Extension CSS not loaded (only JS vendored)
- `scrollX: true` not enabled in DataTables config
- Table width doesn't exceed container width (no horizontal scroll needed)

**How to avoid:**
- Vendor BOTH `fixedcolumns.min.js` AND `fixedcolumns.dataTables.min.css`
- Inline both in template like other assets
- Ensure `scrollX: true` is in DataTables config
- Test with many columns visible to force horizontal scroll

**Warning signs:**
- Console error: "FixedColumns requires scrollX to be enabled"
- Frozen column appears but doesn't stay fixed during scroll
- Layout breaks on pagination/sort

### Pitfall 2: Child Row Content Not Showing Full Data

**What goes wrong:** Child row formatter function shows `undefined` or missing fields

**Why it happens:**
- `row.data()` returns raw data object, but column names may be sanitized (Phase 8 sanitization)
- Column names in JSON use original names, but template expects sanitized names
- Hidden columns may not be in `row.data()` if not rendered

**How to avoid:**
- Use `columnData` array from Jinja2 to map `original_name` to display
- Access row data using original field names (e.g., `rowData['HGVS_C']` not `rowData['HGVS.C']`)
- Test with all hidden columns visible to ensure data availability

**Warning signs:**
- Child row shows "undefined" or blank values
- Only visible columns appear in child row
- Link columns missing href values

### Pitfall 3: Tippy Tooltips Not Re-Appearing After Density Change

**What goes wrong:** Tooltips work initially but break after switching density mode

**Why it happens:**
- Tippy instances attached to specific DOM elements
- Density change triggers DataTables `draw()` which recreates `<td>` elements
- Old Tippy instances still reference destroyed DOM nodes

**How to avoid:**
- Call `initTippyTooltips()` in DataTables `drawCallback`
- Check `if (!cell._tippy)` before creating new instances (prevents duplicates)
- Use `tippy.destroyAll()` if recreating entire table

**Warning signs:**
- Tooltips work on page load but not after pagination/sort/density change
- Console error: "Cannot read property 'getBoundingClientRect' of null"

### Pitfall 4: Column Widths Not Respected with FixedColumns

**What goes wrong:** Column widths set via `columnDefs.width` are ignored when FixedColumns is enabled

**Why it happens:**
- FixedColumns clones columns into separate containers
- Width calculations happen before FixedColumns initializes
- Fixed and scrolling tables have different width constraints

**How to avoid:**
- Set `autoWidth: false` in DataTables config (already present)
- Use explicit pixel widths (e.g., `'80px'`) not percentages
- Call `columns.adjust()` in `initComplete` callback (already present)
- For frozen column, set min-width in CSS as well as columnDefs

**Warning signs:**
- Frozen column narrower than expected
- Column widths change after horizontal scroll
- Fixed and scrolling column headers misaligned

### Pitfall 5: localStorage Density Persisting Across Different Reports

**What goes wrong:** User sets density to "relaxed" in one report, all future reports open in relaxed mode even if they prefer compact

**Why it happens:**
- Single `variantTableDensity` key in localStorage
- No namespacing per VCF file or report instance

**How to avoid:**
- Use namespaced key: `variantTableDensity_${vcfFilename}` or simple global default
- For Phase 15, use global key (simpler) — per-report persistence is Phase 18+ enhancement
- Document that density is a user preference, not per-report setting

**Warning signs:**
- User confused why new report opens in unexpected density
- Different projects have different optimal densities

## Code Examples

Verified patterns from official sources and current codebase:

### Current Column Configuration (Phase 13/14 Baseline)

```javascript
// Source: variantcentrifuge/templates/index.html lines 492-602
var allColumnOriginalNames = {{ column_data | map(attribute='original_name') | list | tojson }};
var defaultHiddenColumnOriginalNames = {{ default_hidden_columns | tojson }};
var columnData = {{ column_data | tojson }};
var hiddenColumnIndices = [];
var dtColumnDefs = [];

// Find indices of columns to hide by default
defaultHiddenColumnOriginalNames.forEach(function(hiddenColName) {
    var index = allColumnOriginalNames.indexOf(hiddenColName);
    if (index !== -1) {
        hiddenColumnIndices.push(index);
    }
});

// Setup default hidden columns
if (hiddenColumnIndices.length > 0) {
    dtColumnDefs.push({
        targets: hiddenColumnIndices,
        visible: false
    });
}

// Badge render functions for semantic color coding (Phase 14)
columnData.forEach(function(col, index) {
    if (col.original_name === 'IMPACT') {
        dtColumnDefs.push({
            targets: index,
            render: function(data, type, row) {
                if (type === 'display' && data) {
                    var colors = {
                        'HIGH': '#dc3545',
                        'MODERATE': '#fd7e14',
                        'LOW': '#f59e0b',
                        'MODIFIER': '#6c757d'
                    };
                    var color = colors[data] || '#6c757d';
                    return '<span class="badge" style="background-color:' + color + ';">' + data + '</span>';
                }
                return data;
            }
        });
    }
    // ... ClinVar, Inheritance_Pattern badges
});
```

### Enhanced Column Configuration (Phase 15 Addition)

```javascript
// Source: Phase 15 requirements + DataTables columnDefs API
// Add to existing columnData.forEach loop:

// Fixed width for CHROM/POS (TABLE-06)
if (col.original_name === 'CHROM') {
    dtColumnDefs.push({
        targets: index,
        width: '70px',
        className: 'monospace'
    });
}

if (col.original_name === 'POS') {
    dtColumnDefs.push({
        targets: index,
        width: '90px',
        className: 'monospace dt-right'
    });
}

// Grow for GENE (TABLE-06)
if (col.original_name === 'GENE') {
    dtColumnDefs.push({
        targets: index,
        className: 'gene-column'  // No fixed width, let it grow
    });
}

// Monospace for REF/ALT (TABLE-08)
if (col.original_name === 'REF' || col.original_name === 'ALT') {
    dtColumnDefs.push({
        targets: index,
        className: 'monospace'
    });
}

// End truncation for HGVS (TABLE-08)
if (col.original_name === 'HGVS_C' || col.original_name === 'HGVS_P') {
    dtColumnDefs.push({
        targets: index,
        width: '180px',
        className: 'monospace',
        render: function(data, type, row) {
            if (type === 'display' && data) {
                var maxLen = 30;
                if (data.length > maxLen) {
                    var truncated = data.substring(0, maxLen) + '...';
                    return '<span class="truncated-cell" style="max-width:180px;" data-tippy-content="' + data + '">' + truncated + '</span>';
                }
            }
            return data;
        }
    });
}

// Right-align numeric scores (TABLE-06)
if (col.original_name.includes('_score') || col.original_name.includes('_phred') || col.original_name.includes('_AF')) {
    dtColumnDefs.push({
        targets: index,
        className: 'dt-right',
        render: function(data, type, row) {
            if (type === 'display' && data && !isNaN(data)) {
                return parseFloat(data).toFixed(3);
            }
            return data;
        }
    });
}

// Middle truncation for VAR_ID (TABLE-08)
if (col.original_name === 'VAR_ID') {
    dtColumnDefs.push({
        targets: index,
        width: '120px',
        render: function(data, type, row) {
            if (type === 'display' && data && data.length > 18) {
                var start = data.substring(0, 10);
                var end = data.substring(data.length - 6);
                return '<span class="middle-truncate" data-tippy-content="' + data + '">' + start + '...' + end + '</span>';
            }
            return data;
        }
    });
}
```

### FixedColumns Integration (TABLE-02)

```javascript
// Source: https://datatables.net/extensions/fixedcolumns/
// Modify DataTable initialization:

var table = new DataTable('#variants_table', {
    autoWidth: false,
    pageLength: 25,
    scrollX: true,
    fixedColumns: {
        left: 1  // Freeze first column (GENE)
    },
    layout: {
        // ... existing layout config
    },
    columnDefs: dtColumnDefs,
    initComplete: function() {
        if (skeleton) { skeleton.style.display = 'none'; }
        this.api().columns.adjust();
    },
    drawCallback: function(settings) {
        initTippyTooltips();
    }
});
```

**CSS for shadow separator:**
```css
/* Source: UX research from Phase 15 CONTEXT.md */
.dtfc-fixed-left {
    box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
    z-index: 2;
}
```

### Child Row Expandable Details (TABLE-03)

```javascript
// Source: https://datatables.net/examples/api/row_details.html
// 1. Add control column (chevron) to columnDefs
dtColumnDefs.unshift({
    targets: 0,
    className: 'dt-control',
    orderable: false,
    data: null,
    defaultContent: '<span class="chevron">›</span>'
});

// 2. Adjust all other target indices by +1 (control column is now first)
// ... existing dtColumnDefs need targets updated

// 3. Format function for child row
function formatChildRow(rowData) {
    var html = '<div class="detail-panel">';

    // Identifiers section
    html += '<div class="detail-section">';
    html += '<h4>Identifiers</h4>';
    html += '<dl>';
    html += '<dt>Variant ID</dt><dd>' + (rowData.VAR_ID || 'N/A') + '</dd>';
    html += '<dt>rsID</dt><dd>' + (rowData.ID || 'N/A') + '</dd>';
    html += '<dt>Position</dt><dd>' + rowData.CHROM + ':' + rowData.POS + '</dd>';
    html += '</dl></div>';

    // Annotations section
    html += '<div class="detail-section">';
    html += '<h4>Annotations</h4>';
    html += '<dl>';
    html += '<dt>Gene</dt><dd>' + (rowData.GENE || 'N/A') + '</dd>';
    html += '<dt>Effect</dt><dd>' + (rowData.EFFECT || 'N/A') + '</dd>';
    html += '<dt>Impact</dt><dd>' + (rowData.IMPACT || 'N/A') + '</dd>';
    html += '<dt>HGVS.c</dt><dd>' + (rowData.HGVS_C || 'N/A') + '</dd>';
    html += '<dt>HGVS.p</dt><dd>' + (rowData.HGVS_P || 'N/A') + '</dd>';
    html += '</dl></div>';

    // Scores section
    html += '<div class="detail-section">';
    html += '<h4>Prediction Scores</h4>';
    html += '<dl>';
    html += '<dt>CADD</dt><dd>' + (rowData.dbNSFP_CADD_phred || 'N/A') + '</dd>';
    html += '<dt>REVEL</dt><dd>' + (rowData.dbNSFP_REVEL_score || 'N/A') + '</dd>';
    html += '<dt>ClinVar</dt><dd>' + (rowData.ClinVar_CLNSIG || rowData.dbNSFP_clinvar_clnsig || 'N/A') + '</dd>';
    html += '</dl></div>';

    // Links section
    html += '<div class="detail-section">';
    html += '<h4>External Resources</h4>';
    html += '<dl>';
    if (rowData.ClinVar) {
        html += '<dt>ClinVar</dt><dd><a href="' + rowData.ClinVar + '" target="_blank">View in ClinVar ↗</a></dd>';
    }
    if (rowData.gnomAD_2) {
        html += '<dt>gnomAD</dt><dd><a href="' + rowData.gnomAD_2 + '" target="_blank">View in gnomAD ↗</a></dd>';
    }
    if (rowData.Varsome) {
        html += '<dt>Varsome</dt><dd><a href="' + rowData.Varsome + '" target="_blank">View in Varsome ↗</a></dd>';
    }
    html += '</dl></div>';

    html += '</div>';
    return html;
}

// 4. Click handler for toggle
table.on('click', 'tbody td.dt-control', function(e) {
    var tr = e.target.closest('tr');
    var row = table.row(tr);

    if (row.child.isShown()) {
        row.child.hide();
        tr.classList.remove('shown');
    } else {
        row.child(formatChildRow(row.data())).show();
        tr.classList.add('shown');
    }
});
```

**CSS for child row (add to <style> block):**
```css
/* Chevron in control column */
td.dt-control {
    cursor: pointer;
    text-align: center;
    width: 30px;
    user-select: none;
}

td.dt-control .chevron {
    display: inline-block;
    font-size: 18px;
    font-weight: bold;
    color: #007bff;
    transition: transform 200ms ease-in-out;
}

tr.shown td.dt-control .chevron {
    transform: rotate(90deg);
}

/* Detail panel multi-column grid */
.detail-panel {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
    gap: 16px;
    padding: 16px;
    background-color: #f8f9fa;
    border-left: 3px solid #007bff;
    animation: slideDown 200ms ease-out;
}

@keyframes slideDown {
    from {
        opacity: 0;
        transform: translateY(-10px);
    }
    to {
        opacity: 1;
        transform: translateY(0);
    }
}

/* Detail section cards */
.detail-section {
    background: white;
    padding: 12px 14px;
    border-radius: 6px;
    border: 1px solid #dee2e6;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
}

.detail-section h4 {
    margin: 0 0 10px 0;
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: #6c757d;
    font-weight: 700;
}

.detail-section dl {
    display: grid;
    grid-template-columns: auto 1fr;
    gap: 6px 12px;
    margin: 0;
    font-size: 13px;
}

.detail-section dt {
    font-weight: 600;
    color: #495057;
}

.detail-section dd {
    margin: 0;
    word-break: break-word;
    color: #212529;
}

.detail-section a {
    color: #007bff;
    text-decoration: none;
}

.detail-section a:hover {
    text-decoration: underline;
}
```

### Content Density Toggle (TABLE-07)

```javascript
// Source: localStorage best practices + DataTables buttons API
// 1. Load saved density on page load (before DataTable init)
var densityModes = ['compact', 'regular', 'relaxed'];
var savedDensity = localStorage.getItem('variantTableDensity') || 'compact';
var currentDensityIndex = densityModes.indexOf(savedDensity);

// 2. Apply density class to table
document.addEventListener('DOMContentLoaded', function() {
    var tableElement = document.getElementById('variants_table');
    if (tableElement) {
        tableElement.classList.add('density-' + savedDensity);
    }
});

// 3. Add density button to layout
layout: {
    topStart: {
        pageLength: { menu: [10, 25, 50, 100, -1] }
    },
    topEnd: {
        search: true,
        buttons: [
            {
                text: 'Density: ' + savedDensity.charAt(0).toUpperCase() + savedDensity.slice(1),
                action: function(e, dt, node, config) {
                    var tableNode = dt.table().node();

                    // Remove current density class
                    tableNode.classList.remove('density-' + densityModes[currentDensityIndex]);

                    // Cycle to next mode
                    currentDensityIndex = (currentDensityIndex + 1) % densityModes.length;
                    var newDensity = densityModes[currentDensityIndex];

                    // Apply new density class
                    tableNode.classList.add('density-' + newDensity);

                    // Save to localStorage
                    localStorage.setItem('variantTableDensity', newDensity);

                    // Update button text
                    node.text('Density: ' + newDensity.charAt(0).toUpperCase() + newDensity.slice(1));

                    // Redraw table to adjust row heights
                    dt.columns.adjust().draw(false);  // false = stay on same page
                }
            },
            {
                extend: 'colvis',
                text: 'Show/Hide Columns'
            }
        ]
    }
}
```

**CSS for density modes:**
```css
/* TABLE-07: Compact/Regular/Relaxed density modes */
#variants_table.density-compact td {
    padding: 4px 6px;
    line-height: 1.3;
}

#variants_table.density-regular td {
    padding: 8px 10px;
    line-height: 1.5;
}

#variants_table.density-relaxed td {
    padding: 12px 14px;
    line-height: 1.7;
}

/* Header padding scales with density */
#variants_table.density-compact th {
    padding: 8px 6px;
}

#variants_table.density-regular th {
    padding: 10px 10px;
}

#variants_table.density-relaxed th {
    padding: 12px 14px;
}
```

### Enhanced Header and Zebra Striping (TABLE-04, TABLE-05)

```css
/* Source: Phase 15 requirements + Phase 14 color palette */
/* TABLE-05: Dark header with white text and bottom border */
#variants_table thead th,
.dt-scroll-head thead th,
.dataTables_scrollHead thead th {
    background: linear-gradient(180deg, #2c3e50 0%, #1a252f 100%);
    color: white;
    border-bottom: 3px solid #007bff;  /* Accent color */
    font-weight: 600;
    padding: 10px 8px;
    text-align: left !important;
}

/* TABLE-04: Zebra striping */
#variants_table tbody tr:nth-child(odd) {
    background-color: #f8f9fa;  /* Light gray */
}

#variants_table tbody tr:nth-child(even) {
    background-color: white;
}

/* Hover highlight */
#variants_table tbody tr:hover {
    background-color: #e3f2fd;  /* Light blue */
}

/* Preserve zebra striping when row expanded */
#variants_table tbody tr.shown + tr {
    background-color: transparent !important;  /* Child row inherits parent bg */
}
```

### Prominent Record Count (TABLE-09)

```css
/* Source: DataTables default info element styling enhancement */
/* Make the "Showing X to Y of Z entries" more prominent */
.dt-info,
.dataTables_info {
    font-size: 14px;
    font-weight: 600;
    color: #212529;
    padding: 8px 0;
}

.dt-info::before,
.dataTables_info::before {
    content: "📊 ";  /* Optional icon for visual prominence */
}
```

**Note:** DataTables automatically generates the info element with text like "Showing 1 to 25 of 137 entries" based on `pageLength` and total records. No JS needed, just style enhancement.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Plotly.js for charts (~3.5MB) | Chart.js (~200KB) | Phase 13 (Feb 2026) | 94% bundle size reduction, faster page load |
| jQuery-based DataTables v1.10 | Vanilla JS DataTables v2.2 | Phase 13 (Feb 2026) | Modern API, no jQuery user code, layout options replace dom strings |
| Hover-expand with yellow background | Tippy.js tooltips on truncated cells | Phase 13 (Feb 2026) | Viewport-aware, keyboard accessible, dismissable with Escape |
| Plain IMPACT text | Color badges (Phase 14) | Phase 14 (Feb 2026) | Semantic color coding, better scannability |
| 10 rows per page default | 25 rows per page | TABLE-10 (Phase 15) | Reduces pagination clicks for power users |
| No sticky columns | FixedColumns extension | Phase 15 (planned) | GENE column visible during horizontal scroll |
| No row details | Child rows with grouped fields | Phase 15 (planned) | Progressive disclosure, cleaner main table |
| Fixed padding | Density modes (Compact/Regular/Relaxed) | Phase 15 (planned) | User preference for information density |

**Deprecated/outdated:**
- DataTables v1 `dom` option syntax → Replaced by `layout` API in v2
- jQuery plugin syntax `$('#table').DataTable()` → Vanilla `new DataTable('#table')` in v2 (both still work)
- Plotly.js for simple charts → Chart.js for lighter bundle
- Manual tooltip positioning → Tippy.js handles viewport clipping

## Open Questions

### Question 1: FixedColumns Compatibility with DataTables 2.2.2

**What we know:**
- DataTables v2.2.2 released Nov 2024
- FixedColumns 5.0.5 released Oct 2024
- Both are latest stable versions

**What's unclear:**
- Exact compatibility matrix not published on DataTables site
- Some forum posts mention breaking changes between v1 and v2

**Recommendation:**
- Test FixedColumns 5.0.5 with DataTables 2.2.2 in dev environment before vendoring
- If incompatible, check DataTables releases page for compatible FixedColumns version
- Fallback: Use CSS `position: sticky` on first column (less robust but no dependency)

### Question 2: Column Index Shift with Control Column

**What we know:**
- Adding control column (chevron) to columnDefs as `targets: 0` shifts all indices
- Current columnDefs use `forEach(function(col, index))` to assign targets

**What's unclear:**
- Should control column be inserted before or after loop?
- Do existing badge render functions need index adjustment?

**Recommendation:**
- Insert control column definition BEFORE existing columnDefs loop
- Adjust all subsequent `targets: index` to `targets: index + 1`
- Alternative: Use column names instead of indices if DataTables supports (check API)

### Question 3: Middle Truncation Browser Support

**What we know:**
- CSS `text-overflow: ellipsis` only supports end truncation
- Middle truncation requires JavaScript to split string

**What's unclear:**
- Performance impact of JS render function on 1000+ rows
- Whether Tippy tooltip still works with custom truncated HTML

**Recommendation:**
- Implement middle truncation in `render` function (as shown in code examples)
- Benchmark with 1000+ row test dataset
- If too slow, fallback to end truncation for VAR_ID

### Question 4: Density Button Text Length

**What we know:**
- Button text "Density: Compact" is ~17 characters
- DataTables toolbar may have limited space with search + colvis button

**What's unclear:**
- Does button text wrap or truncate on narrow screens?
- Should button use icon instead of text label?

**Recommendation:**
- Start with text label "Density: Compact" as specified
- Test on 1024px width (common laptop resolution)
- If text wraps, consider shorter label "🔲 Compact" with icon

## Sources

### Primary (HIGH confidence)

- DataTables v2.2.2 Documentation - https://datatables.net/
- DataTables FixedColumns Extension - https://datatables.net/extensions/fixedcolumns/
- DataTables Child Rows Example - https://datatables.net/examples/api/row_details.html
- DataTables ColumnDefs API - https://datatables.net/reference/option/columnDefs
- Tippy.js v6.3.7 Documentation - https://atomiks.github.io/tippyjs/ (already vendored in Phase 13)
- Current template implementation - variantcentrifuge/templates/index.html (lines 1-840)

### Secondary (MEDIUM confidence)

- jsDelivr CDN for FixedColumns - https://www.jsdelivr.com/package/npm/datatables.net-fixedcolumns-dt
- DataTables forums on column widths - https://datatables.net/forums/discussion/71994/set-column-max-width
- DataTables ellipsis plugin - https://datatables.net/plug-ins/dataRender/ellipsis
- Phase 13 verification report - .planning/phases/13-js-stack-modernization/13-VERIFICATION.md

### Tertiary (LOW confidence)

- Web search results on DataTables v2 responsive details - https://datatables.net/reference/option/responsive.details
  (Note: Phase 15 uses manual child rows, not Responsive extension — this is for reference only)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official DataTables extensions, verified versions, clear installation path
- Architecture: HIGH - Current template structure documented, DataTables API well-documented
- Pitfalls: MEDIUM - Based on forum discussions and extrapolation from Phase 13 fixes, not direct experience

**Research date:** 2026-02-17
**Valid until:** 60 days (DataTables v2 is stable, minimal churn expected)

**Key dependencies:**
- Phase 13 completion (DataTables v2.2.2, Tippy.js vendored) ✓ Complete
- Phase 14 completion (Badge render functions, semantic colors) ✓ Complete
- FixedColumns extension compatibility testing (BLOCKER if incompatible)
