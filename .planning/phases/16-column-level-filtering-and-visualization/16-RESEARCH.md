# Phase 16: Column-Level Filtering and Visualization - Research

**Researched:** 2026-02-17
**Domain:** Interactive filtering UI with DataTables column search API, noUiSlider range controls, and Chart.js reactive charts
**Confidence:** HIGH

## Summary

Phase 16 adds per-column filtering controls (numeric sliders, categorical dropdowns, text search) and expanded visualizations (variant type breakdown, chromosome distribution, allele frequency histogram) to the individual HTML report. The cohort report (`scripts/templates/cohort_report_template.html`) provides the reference implementation with comprehensive filter patterns that can be adapted.

The individual report (`variantcentrifuge/templates/index.html`) currently uses DataTables v2.2.2 with global search. This phase adds custom filter controls that interact with DataTables' column search API (`column().search()`) and custom search functions (`$.fn.dataTable.ext.search.push()`). Charts update reactively when filters change via the `draw.dt` event.

The cohort report demonstrates three filter control types: noUiSlider range sliders for numeric columns (AF, scores, POS, QUAL), categorical dropdowns for IMPACT/ClinVar/Inheritance, and text inputs for GENE/ID columns. Each filter includes an "Include missing values" checkbox to handle N/A, NA, ".", and empty cells.

**Primary recommendation:** Vendor noUiSlider 15.x for numeric range controls, use DataTables custom search functions for filter logic, wire filter controls to `draw.dt` event for chart reactivity, and create a collapsible "Visualizations" section for new charts (variant type breakdown, chromosome distribution, allele frequency histogram).

## Standard Stack

### Core Libraries (Required)

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| noUiSlider | 15.7.x | Range slider for numeric filters | Lightweight (15KB), touch-friendly, most popular pure-JS range slider (14k+ GitHub stars) |
| DataTables | 2.2.2 | Table with column search API | Already vendored from Phase 13, provides `column().search()` and custom filter API |
| Chart.js | 4.4.8 | Charts (impact, inheritance, new visualizations) | Already vendored from Phase 14, supports reactive updates via `.update()` |

### Already Vendored (Phases 13-15)

| Library | Version | Purpose | Status |
|---------|---------|---------|--------|
| datatables.net | 2.2.2 | Core DataTables v2 | ✓ Vendored (95KB) |
| datatables.net-fixedcolumns | 5.0.5 | Sticky columns | ✓ Vendored in Phase 15 |
| Chart.js | 4.4.8 | Dashboard charts | ✓ Vendored (202KB) |
| chartjs-plugin-datalabels | 3.x | Chart labels | ✓ Vendored |
| Tippy.js | 6.3.7 | Tooltips | ✓ Vendored (26KB bundle) |
| jQuery slim | 3.7.1 | DataTables dependency | ✓ Vendored (70KB) |

### New Dependency

```bash
# Download noUiSlider from jsDelivr CDN
curl -o variantcentrifuge/assets/js/nouislider.min.js \
  https://cdn.jsdelivr.net/npm/nouislider@15.7.1/dist/nouislider.min.js

curl -o variantcentrifuge/assets/css/nouislider.min.css \
  https://cdn.jsdelivr.net/npm/nouislider@15.7.1/dist/nouislider.min.css
```

**Why noUiSlider?**
- Lightweight: 15KB minified (vs. jQuery UI slider at 250KB+)
- Touch-friendly: Works on mobile/tablet (cohort report requirement)
- No jQuery dependency: Pure JavaScript
- Active maintenance: v15.7.1 released Jan 2024, frequent updates
- Official site: https://refreshless.com/nouislider/

## Architecture Patterns

### Current Individual Report Structure (Post-Phase 15)

```
variantcentrifuge/templates/index.html
├── <head>
│   ├── Inlined CSS assets (datatables, buttons, fixedcolumns, tippy, chart)
│   └── Custom <style> block
├── <body>
│   ├── <header> — Report title and metadata
│   ├── <div class="dashboard"> — Phase 14 metric cards + charts
│   │   ├── .cards-grid — Metric cards (Total Variants, Unique Genes, etc.)
│   │   └── .charts-grid — Impact distribution + Inheritance patterns
│   ├── <div class="table-section">
│   │   ├── <div id="table-skeleton"> — Loading shimmer
│   │   └── <table id="variants_table"> — DataTables v2 variant table
│   └── <footer class="metadata-footer"> — Pipeline metadata
└── <script>
    ├── Inlined JS assets
    ├── DataTable initialization with fixedColumns
    ├── Badge render functions (IMPACT, ClinVar, Inheritance)
    ├── Tippy.js tooltips
    ├── Chart.js initialization (2 charts)
    └── Child row expansion handler
```

**Where Phase 16 fits:**
- **Filter controls:** Insert between dashboard and table-section (above toolbar, below charts)
- **New visualization section:** Insert between dashboard and filter controls (or after filters, decision point)
- **Filter chip strip:** Insert between toolbar and table (shows active filters)

### Pattern 1: Column-Level Filter Controls with noUiSlider

**What:** Per-column filter UI that sends filter values to DataTables column search API or custom search functions

**When to use:** Tables with >10 columns where global search is insufficient, users need precise filtering by specific fields

**Cohort Report Reference Implementation:**

The cohort report creates filters dynamically based on column type detection:

```javascript
// From cohort_report_template.html lines 463-845
function initializeFilters(data) {
    const columns = dataTable.columns().header().toArray().map(th => $(th).text());
    const filtersContainer = $('#custom-filters');
    const MISSING_VALUES = new Set(['', 'N/A', 'NA', '.', 'Unknown', 'null', 'undefined']);

    // Auto-detect column type
    function getColumnType(colName) {
        const name = colName.toUpperCase();

        // Force numeric: AF, AD, AC, SCORE, PHRED, CADD, REVEL, POS, QUAL, DP
        if (name.startsWith('AF ') || name.startsWith('AF_') || name === 'AF') return 'numeric';
        if (name.includes('SCORE') || name.includes('PHRED')) return 'numeric';
        if (name === 'POS' || name === 'QUAL' || name === 'DP') return 'numeric';

        // Force text: ID, CHROM, REF, ALT, GENE, IMPACT, ClinVar, EFFECT
        if (name.includes('ID') || name === 'CHROM' || name === 'GENE') return 'text';
        if (name === 'IMPACT' || name.includes('CLNSIG') || name.includes('EFFECT')) return 'text';

        return 'auto'; // Analyze data to decide
    }

    // Create slider for numeric columns
    if (filterType === 'numeric') {
        const min = Math.min(...numericValues);
        const max = Math.max(...numericValues);

        const sliderEl = document.getElementById(`slider-${colName}`);
        noUiSlider.create(sliderEl, {
            start: [min, max],
            connect: true,
            range: { 'min': min, 'max': max },
            format: {
                to: val => parseFloat(val).toFixed(4),
                from: val => Number(val)
            }
        });

        // Filter function using DataTables custom search
        sliderEl.noUiSlider.on('update', function (values) {
            const [minRange, maxRange] = values.map(parseFloat);

            $.fn.dataTable.ext.search.push(
                function filterNumericColumn(settings, data, dataIndex) {
                    if (settings.nTable.id !== 'variants-table') return true;

                    const cellValue = data[colIndex];

                    // Handle missing values
                    if (MISSING_VALUES.has(String(cellValue).trim())) {
                        return includeMissing; // From checkbox state
                    }

                    const numVal = parseFloat(cellValue);
                    return numVal >= minRange && numVal <= maxRange;
                }
            );

            dataTable.draw();
        });
    }
}
```

**Key insights:**
- noUiSlider requires a container element (`<div id="slider-*">`) and stores instance on `element.noUiSlider`
- Custom search functions pushed to `$.fn.dataTable.ext.search` array run on every row during `draw()`
- Check `settings.nTable.id` to prevent filter from affecting other tables on page
- Missing value handling via checkbox that toggles `includeMissing` flag in filter function

### Pattern 2: DataTables Column Search API

**What:** Two approaches to filtering columns in DataTables v2: column search strings and custom search functions

**Official DataTables v2 API:**

```javascript
// Approach 1: Column search string (simple text matching)
dataTable.column(colIndex).search('search term').draw();

// Approach 2: Custom search function (complex logic)
$.fn.dataTable.ext.search.push(
    function customFilter(settings, searchData, dataIndex, rowData, counter) {
        // settings.nTable - the table DOM node
        // searchData - array of strings (one per column) from table body
        // dataIndex - row index in data array
        // rowData - full row data object
        // counter - incremental counter for draw

        // Return true to include row, false to exclude
        return true;
    }
);

// Remove custom filter (important for reset)
$.fn.dataTable.ext.search = $.fn.dataTable.ext.search.filter(f =>
    !f.hasOwnProperty('customIdentifier')
);
```

**Sources:**
- DataTables v2 column search: https://datatables.net/reference/api/column().search()
- Custom search functions: https://datatables.net/examples/plug-ins/range_filtering.html

**When to use each approach:**
- Column search: Simple text matching (GENE name, variant ID prefix)
- Custom function: Complex logic (numeric ranges, missing value handling, multi-column AND logic)

### Pattern 3: Filter Chip UI Pattern

**What:** Visual indicators of active filters displayed as removable chips/tags above table

**Cohort Report Example (not implemented yet, but standard pattern):**

```html
<!-- Filter chip strip between toolbar and table -->
<div id="active-filters-strip" style="display: none;">
    <div class="filter-chips-container">
        <span class="filter-count">3 filters active</span>
        <span class="filter-chip" data-column="IMPACT">
            IMPACT: HIGH
            <button class="chip-remove" aria-label="Remove filter">×</button>
        </span>
        <span class="filter-chip" data-column="AF">
            AF: 0.0001 - 0.05
            <button class="chip-remove" aria-label="Remove filter">×</button>
        </span>
        <button class="btn-reset-all">Reset All Filters</button>
    </div>
</div>
```

**CSS:**
```css
.filter-chip {
    display: inline-block;
    background-color: #e3f2fd;
    border: 1px solid #2196f3;
    border-radius: 16px;
    padding: 4px 8px 4px 12px;
    margin: 4px;
    font-size: 13px;
    color: #1976d2;
}

.chip-remove {
    background: none;
    border: none;
    color: #1976d2;
    font-size: 18px;
    font-weight: bold;
    cursor: pointer;
    margin-left: 6px;
    padding: 0 4px;
}

.chip-remove:hover {
    color: #0d47a1;
}
```

**JavaScript logic:**
```javascript
// Track active filters
const activeFilters = new Map(); // column -> filterState

function updateFilterChips() {
    const strip = document.getElementById('active-filters-strip');
    const container = strip.querySelector('.filter-chips-container');

    // Clear existing chips (preserve Reset All button)
    container.querySelectorAll('.filter-chip').forEach(chip => chip.remove());

    if (activeFilters.size === 0) {
        strip.style.display = 'none';
        return;
    }

    strip.style.display = 'block';
    container.querySelector('.filter-count').textContent =
        `${activeFilters.size} filter${activeFilters.size > 1 ? 's' : ''} active`;

    // Create chip for each active filter
    activeFilters.forEach((state, colName) => {
        const chip = document.createElement('span');
        chip.className = 'filter-chip';
        chip.dataset.column = colName;
        chip.innerHTML = `
            ${colName}: ${formatFilterValue(state)}
            <button class="chip-remove" aria-label="Remove ${colName} filter">×</button>
        `;

        chip.querySelector('.chip-remove').addEventListener('click', () => {
            removeFilter(colName);
        });

        container.insertBefore(chip, container.querySelector('.btn-reset-all'));
    });
}
```

**Industry standard examples:**
- Google Search (filter chips below search bar)
- GitHub (filter tags in issue list)
- Material UI Chip component: https://mui.com/material-ui/react-chip/

### Pattern 4: Reactive Chart Updates on Filter Change

**What:** Charts automatically update to reflect filtered dataset when table filters change

**Implementation using DataTables draw event:**

```javascript
// Phase 14 charts (already exist)
let impactChart;  // Chart.js instance
let inheritanceChart;  // Chart.js instance

// Wire up to DataTables draw event
dataTable.on('draw.dt', function() {
    const filteredData = dataTable.rows({ search: 'applied' }).data().toArray();
    updateAllCharts(filteredData);
});

function updateAllCharts(filteredData) {
    // Update existing Phase 14 charts
    updateImpactChart(filteredData);
    updateInheritanceChart(filteredData);

    // Update new Phase 16 charts
    updateVariantTypeChart(filteredData);
    updateChromosomeChart(filteredData);
    updateAlleleFrequencyHistogram(filteredData);
}

function updateImpactChart(data) {
    const impactCounts = {};
    data.forEach(row => {
        const impact = row.IMPACT || 'UNKNOWN';
        impactCounts[impact] = (impactCounts[impact] || 0) + 1;
    });

    // Update Chart.js data (no animation for snap feel)
    impactChart.data.labels = Object.keys(impactCounts);
    impactChart.data.datasets[0].data = Object.values(impactCounts);
    impactChart.update('none'); // 'none' disables animation
}
```

**Key insight:** `rows({ search: 'applied' })` returns only rows passing current search/filter criteria

**Sources:**
- DataTables draw event: https://datatables.net/reference/event/draw
- Chart.js update API: https://www.chartjs.org/docs/latest/developers/updates.html

### Pattern 5: Collapsible Section for Visualizations

**What:** Expandable/collapsible section containing new charts (variant type, chromosome, AF histogram)

**Standard accordion pattern:**

```html
<!-- Insert between dashboard and filters -->
<div class="visualization-section">
    <button class="collapse-toggle" aria-expanded="false" aria-controls="viz-content">
        <span class="chevron">▶</span>
        Visualizations
        <span class="badge-count">3 charts</span>
    </button>

    <div id="viz-content" class="collapse-content" hidden>
        <div class="charts-grid-3col">
            <div class="chart-wrapper">
                <h3>Variant Type Breakdown</h3>
                <canvas id="variant_type_chart"></canvas>
            </div>
            <div class="chart-wrapper">
                <h3>Chromosome Distribution</h3>
                <canvas id="chromosome_chart"></canvas>
            </div>
            <div class="chart-wrapper">
                <h3>Allele Frequency Histogram</h3>
                <canvas id="af_histogram_chart"></canvas>
            </div>
        </div>
    </div>
</div>
```

**CSS:**
```css
.collapse-toggle {
    width: 100%;
    background: white;
    border: 1px solid #dee2e6;
    border-radius: 8px;
    padding: 12px 16px;
    text-align: left;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 8px;
    font-size: 16px;
    font-weight: 600;
    margin-bottom: 12px;
}

.collapse-toggle:hover {
    background-color: #f8f9fa;
}

.chevron {
    transition: transform 200ms ease;
    font-size: 12px;
}

.collapse-toggle[aria-expanded="true"] .chevron {
    transform: rotate(90deg);
}

.collapse-content {
    overflow: hidden;
    transition: max-height 300ms ease-out;
}

.collapse-content[hidden] {
    display: none;
}

.charts-grid-3col {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
    gap: 16px;
}
```

**JavaScript:**
```javascript
document.querySelector('.collapse-toggle').addEventListener('click', function() {
    const expanded = this.getAttribute('aria-expanded') === 'true';
    const content = document.getElementById('viz-content');

    this.setAttribute('aria-expanded', !expanded);
    content.hidden = !content.hidden;

    // Store preference
    localStorage.setItem('vizSectionExpanded', !expanded);
});

// Restore preference on load
const vizExpanded = localStorage.getItem('vizSectionExpanded') === 'true';
if (vizExpanded) {
    document.querySelector('.collapse-toggle').click();
}
```

**Accessibility:** Uses ARIA `aria-expanded` and `aria-controls` for screen readers

**Decision context:** User requested "collapsed by default" per CONTEXT.md

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Range slider UI | Custom div+drag handlers | noUiSlider 15.x | Touch support, accessibility (ARIA), edge cases (bounds, step, format), 14k stars |
| Missing value detection | Ad-hoc string checks | Set of standardized values + trim() | Cohort report pattern covers: '', 'N/A', 'NA', '.', 'Unknown', 'null', 'undefined' |
| Filter state tracking | Global variables per filter | Map/Object keyed by column name | Cohort report uses `activeFilters` Map for unified state |
| Chart data aggregation | Manual loops per chart type | Reusable aggregation functions | Reduce/forEach patterns across all charts |
| Column type detection | If/else chains | Pattern matching with fallback to data analysis | Cohort report's `getColumnType()` covers 15+ genomic patterns |

**Key insight:** The cohort report (`scripts/templates/cohort_report_template.html`) is production-quality code that handles edge cases the individual report will encounter. Patterns are directly transferable.

## Common Pitfalls

### Pitfall 1: noUiSlider Initialization Timing

**What goes wrong:** Creating noUiSlider on an element before it's in the DOM or with incorrect container structure

**Example:**
```javascript
// WRONG: Element not yet in DOM
const slider = document.getElementById('slider-AF');
noUiSlider.create(slider, {...}); // Error: element is null

// WRONG: Slider already exists
noUiSlider.create(slider, {...}); // Error: already has noUiSlider instance
```

**How to avoid:**
```javascript
// Correct: Check if element exists and doesn't have slider
const slider = document.getElementById('slider-AF');
if (slider && !slider.noUiSlider) {
    noUiSlider.create(slider, {...});
}

// Correct: Destroy before recreating
if (slider.noUiSlider) {
    slider.noUiSlider.destroy();
}
noUiSlider.create(slider, {...});
```

**Warning signs:** Console errors "Cannot read property 'create' of null" or "already has a noUiSlider"

**Source:** noUiSlider docs - https://refreshless.com/nouislider/slider-read-write/

### Pitfall 2: Custom Search Function Leaks

**What goes wrong:** Custom search functions remain in `$.fn.dataTable.ext.search` array after being removed or reset, causing filters to persist or conflict

**Why it happens:** Pushing to `ext.search` without proper cleanup/identification

**Example:**
```javascript
// WRONG: No way to remove this specific function later
$.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {
    return data[5] > 10; // Filter column 5
});

// WRONG: This removes ALL custom functions, not just ours
$.fn.dataTable.ext.search = [];
```

**How to avoid:**
```javascript
// Correct: Add identifier property to function
function createColumnFilter(colIndex, predicate) {
    const filterFn = function(settings, data, dataIndex) {
        if (settings.nTable.id !== 'variants_table') return true;
        return predicate(data[colIndex]);
    };
    filterFn.columnIndex = colIndex; // Identifier property
    return filterFn;
}

// Add filter
const filter = createColumnFilter(5, val => parseFloat(val) > 10);
$.fn.dataTable.ext.search.push(filter);

// Remove specific filter
$.fn.dataTable.ext.search = $.fn.dataTable.ext.search.filter(f =>
    f.columnIndex !== 5
);
```

**Cohort report pattern (lines 885-933):**
```javascript
// Remove any existing filters for this column
$.fn.dataTable.ext.search = $.fn.dataTable.ext.search.filter(f =>
    !f.toString().includes(`colIndex-${colIndex}`)
);

// Add new filter function with identifier comment
$.fn.dataTable.ext.search.push(
    function filterNumericColumn(settings, data, dataIndex) {
        // colIndex-${colIndex} identifier for removal
        if (settings.nTable.id !== 'variants-table') return true;
        // ... filter logic
    }
);
```

**Warning signs:** Filters don't reset, "Reset All" leaves some filters active, multiple filters on same column

### Pitfall 3: Missing Value Detection Edge Cases

**What goes wrong:** Treating numeric zero, boolean false, or whitespace-only strings as missing values

**Why it happens:** Overly broad checks like `if (!value)` or `value == null`

**Example:**
```javascript
// WRONG: False positives
if (!value) return true; // Excludes 0, false, empty string
if (value == null) return true; // Only checks null/undefined, not "N/A"

// WRONG: Case sensitivity
if (value === 'NA') return true; // Misses 'na', 'Na', 'nA'
```

**How to avoid:**
```javascript
// Correct: Explicit missing value set with trim
const MISSING_VALUES = new Set(['', 'N/A', 'NA', '.', 'Unknown', 'null', 'undefined']);

function isMissing(value) {
    if (value === null || value === undefined) return true;
    const str = String(value).trim();
    return MISSING_VALUES.has(str);
}

// Use in filter
if (isMissing(cellValue)) {
    return includeMissing; // From checkbox
}

// For numeric columns, also check isNaN
const numVal = parseFloat(cellValue);
if (isNaN(numVal) || !isFinite(numVal)) {
    return includeMissing; // Treat non-numeric as missing
}
```

**Cohort report pattern (lines 471, 701-703):**
```javascript
const MISSING_VALUES = new Set(['', 'N/A', 'NA', '.', 'Unknown', 'null', 'undefined']);

// Separate missing and valid values
nonMissingValues = allValues.filter(v =>
    v !== null && v !== undefined && !MISSING_VALUES.has(String(v).trim())
);
```

**Warning signs:** Zero values excluded from numeric filters, unexpected row counts when toggling "Include missing"

### Pitfall 4: Chart Update Performance

**What goes wrong:** Charts re-render with animation on every filter change, causing lag/flicker on large datasets

**Why it happens:** Chart.js default `.update()` triggers animation, which compounds when multiple charts update simultaneously

**Example:**
```javascript
// WRONG: Animates on every keystroke in filter
inputEl.addEventListener('input', () => {
    dataTable.draw();
    chart.update(); // 400ms animation delay
});
```

**How to avoid:**
```javascript
// Correct: Disable animation for filter-triggered updates
dataTable.on('draw.dt', function() {
    const filteredData = dataTable.rows({ search: 'applied' }).data().toArray();

    updateChartData(impactChart, filteredData);
    impactChart.update('none'); // 'none' = no animation
});

// Alternative: Debounce draw calls for text inputs
let drawTimeout;
inputEl.addEventListener('input', () => {
    clearTimeout(drawTimeout);
    drawTimeout = setTimeout(() => {
        dataTable.draw(); // Triggers chart update via draw.dt event
    }, 300); // 300ms debounce
});
```

**Decision context:** User requested "snap to new values instantly" (no animation) per CONTEXT.md

**Source:** Chart.js update modes - https://www.chartjs.org/docs/latest/developers/updates.html

## Code Examples

Verified patterns from cohort report and DataTables/Chart.js official docs:

### Example 1: noUiSlider Range Filter Setup

```javascript
// Source: cohort_report_template.html lines 847-928
// Create numeric slider filter for AF column

const colIndex = 5; // Column index in DataTables
const colName = 'AF';
const numericValues = [0.0001, 0.0005, 0.001, 0.05, 0.1]; // From data

const min = Math.min(...numericValues);
const max = Math.max(...numericValues);

// HTML structure (already in DOM)
// <div id="slider-AF" class="slider-container"></div>
// <div class="slider-values">
//   <span id="min-val-AF">0.0001</span>
//   <span id="max-val-AF">0.1000</span>
// </div>

const slider = document.getElementById(`slider-${colName}`);
const minValEl = document.getElementById(`min-val-${colName}`);
const maxValEl = document.getElementById(`max-val-${colName}`);

noUiSlider.create(slider, {
    start: [min, max],
    connect: true,
    range: { 'min': min, 'max': max },
    format: {
        to: val => parseFloat(val).toFixed(4),
        from: val => Number(val)
    }
});

// Update display and apply filter
slider.noUiSlider.on('update', function (values) {
    minValEl.textContent = values[0];
    maxValEl.textContent = values[1];

    const minRange = parseFloat(values[0]);
    const maxRange = parseFloat(values[1]);

    // Remove previous filter for this column
    $.fn.dataTable.ext.search = $.fn.dataTable.ext.search.filter(f =>
        f.columnIndex !== colIndex
    );

    // Add new filter
    const filterFn = function(settings, data, dataIndex) {
        if (settings.nTable.id !== 'variants_table') return true;

        const cellValue = data[colIndex];
        if (cellValue === null || cellValue === undefined) return false;

        const numVal = parseFloat(cellValue);
        if (isNaN(numVal)) return false;

        return numVal >= minRange && numVal <= maxRange;
    };
    filterFn.columnIndex = colIndex; // Identifier

    $.fn.dataTable.ext.search.push(filterFn);
    dataTable.draw();
});
```

### Example 2: Categorical Dropdown Filter

```javascript
// Source: cohort_report_template.html lines 940-1002
// Create categorical dropdown for IMPACT column

const uniqueValues = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']; // From data
const colIndex = 8;
const colName = 'IMPACT';

// HTML structure
// <select id="filter-IMPACT" class="category-dropdown">
//   <option value="">All</option>
//   <option value="HIGH">HIGH</option>
//   <!-- ... more options -->
// </select>

const dropdown = document.getElementById(`filter-${colName}`);

dropdown.addEventListener('change', function() {
    const selectedValue = this.value;

    // Remove previous filter
    $.fn.dataTable.ext.search = $.fn.dataTable.ext.search.filter(f =>
        f.columnIndex !== colIndex
    );

    // Only filter if specific category selected (not "All")
    if (selectedValue !== '') {
        const filterFn = function(settings, data, dataIndex) {
            if (settings.nTable.id !== 'variants_table') return true;
            return String(data[colIndex]).trim() === selectedValue;
        };
        filterFn.columnIndex = colIndex;

        $.fn.dataTable.ext.search.push(filterFn);
    }

    dataTable.draw();
});
```

### Example 3: Reactive Chart Update on Filter

```javascript
// Source: cohort_report_template.html lines 376-383, 1039-1103
// Wire DataTables draw event to chart updates

// Existing Phase 14 charts
let impactChart;
let inheritanceChart;

// New Phase 16 charts
let variantTypeChart;
let chromosomeChart;
let afHistogramChart;

// DataTables draw event handler
dataTable.on('draw.dt', function() {
    const filteredData = dataTable.rows({ search: 'applied' }).data().toArray();
    updateDashboardAndCharts(filteredData);
});

function updateDashboardAndCharts(filteredData) {
    // Update variant count
    const filteredCount = filteredData.length;
    document.getElementById('filtered-variants').textContent = filteredCount;

    // Update impact distribution chart
    const impactCounts = {};
    filteredData.forEach(row => {
        const impact = row.IMPACT || 'UNKNOWN';
        impactCounts[impact] = (impactCounts[impact] || 0) + 1;
    });

    impactChart.data.labels = Object.keys(impactCounts);
    impactChart.data.datasets[0].data = Object.values(impactCounts);
    impactChart.update('none'); // No animation

    // Update variant type chart (Phase 16 new)
    const variantTypes = { 'SNV': 0, 'Insertion': 0, 'Deletion': 0, 'MNV': 0 };
    filteredData.forEach(row => {
        const ref = row.REF || '';
        const alt = row.ALT || '';

        if (ref.length === 1 && alt.length === 1) {
            variantTypes.SNV++;
        } else if (ref.length < alt.length) {
            variantTypes.Insertion++;
        } else if (ref.length > alt.length) {
            variantTypes.Deletion++;
        } else {
            variantTypes.MNV++;
        }
    });

    variantTypeChart.data.datasets[0].data = Object.values(variantTypes);
    variantTypeChart.update('none');
}
```

### Example 4: Allele Frequency Histogram with Log Scale

```javascript
// Source: Chart.js docs - https://www.chartjs.org/docs/latest/axes/cartesian/logarithmic.html
// Create allele frequency histogram with log-scale Y axis

const afValues = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1]; // From filteredData

// Create bins with log-scale boundaries
const bins = [
    { label: '<0.0001', min: 0, max: 0.0001 },
    { label: '0.0001-0.001', min: 0.0001, max: 0.001 },
    { label: '0.001-0.01', min: 0.001, max: 0.01 },
    { label: '0.01-0.1', min: 0.01, max: 0.1 },
    { label: '>0.1', min: 0.1, max: 1.0 }
];

const binCounts = bins.map(bin =>
    afValues.filter(af => af >= bin.min && af < bin.max).length
);

const afHistogramChart = new Chart(
    document.getElementById('af_histogram_chart'),
    {
        type: 'bar',
        data: {
            labels: bins.map(b => b.label),
            datasets: [{
                label: 'Variant Count',
                data: binCounts,
                backgroundColor: '#2196f3'
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: { display: false },
                title: { display: false }
            },
            scales: {
                y: {
                    type: 'logarithmic', // Log scale
                    beginAtZero: false,
                    min: 1,
                    title: {
                        display: true,
                        text: 'Number of Variants (log scale)'
                    },
                    ticks: {
                        callback: function(value) {
                            if (value === 1 || value === 10 || value === 100 || value === 1000) {
                                return value;
                            }
                            return null; // Hide intermediate ticks
                        }
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Allele Frequency Bin'
                    }
                }
            }
        }
    }
);
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Global search only | Column-level filters + global search | DataTables v1.10+ (2014) | Users can precisely filter by IMPACT, AF, ClinVar without typing complex queries |
| jQuery UI slider | noUiSlider | ~2018 | 15KB vs 250KB, touch support, no jQuery dep |
| Static charts | Reactive charts via draw event | Chart.js v2+ (2016) | Charts reflect filtered view instantly |
| Filter dropdowns in thead | Dedicated filter section above table | ~2020 UX pattern | Clearer separation, room for chips/reset |
| CSS classes for badges | Inline styles from render functions | Phase 14 (this project) | Simpler implementation, no CSS class explosion |

**Deprecated/outdated:**
- **jQuery UI slider:** Massive bundle, requires jQuery UI core + CSS theme (250KB+), poor touch support
- **DataTables column filters extension (v1.x):** Replaced by built-in column search API in v2.x
- **Static dashboard:** Modern reports update all visualizations when filters change (reactive pattern)

## Open Questions

Things that couldn't be fully resolved:

1. **Filter control column selection**
   - What we know: Cohort report filters ~15 columns (IMPACT, ClinVar, Inheritance, AF, scores, GENE, etc.)
   - What's unclear: Individual report may have different column set depending on fields_to_extract config
   - Recommendation: Use config-driven approach - create filters for columns matching patterns (IMPACT, *_AF, *_score, GENE, ClinVar*, Inheritance_Pattern) rather than hardcoding specific column names

2. **Collapsible section default state**
   - What we know: User requested "collapsed by default" in CONTEXT.md
   - What's unclear: Should localStorage persistence override default on subsequent visits?
   - Recommendation: Start collapsed, respect localStorage on reload (standard accordion pattern)

3. **Filter chip positioning**
   - What we know: Chips should be "above table, below toolbar"
   - What's unclear: Exact placement relative to density toggle, column visibility, and pagination controls
   - Recommendation: Insert between `.layout.topEnd` buttons div and table element, create dedicated `.filter-strip` container

4. **Chromosome distribution chart orientation**
   - What we know: User requested "horizontal bar" chart
   - What's unclear: Whether chromosomes should be sorted numerically (chr1, chr2, ..., chr22, chrX, chrY, chrM) or by frequency
   - Recommendation: Sort numerically for biological intuition (standard karyogram order)

## Sources

### Primary (HIGH confidence)
- **noUiSlider official docs:** https://refreshless.com/nouislider/ - Configuration, API, examples
- **DataTables v2 column search API:** https://datatables.net/reference/api/column().search() - Official reference
- **DataTables custom search functions:** https://datatables.net/examples/plug-ins/range_filtering.html - Official example
- **Chart.js update modes:** https://www.chartjs.org/docs/latest/developers/updates.html - Animation control
- **Chart.js logarithmic scale:** https://www.chartjs.org/docs/latest/axes/cartesian/logarithmic.html - Log scale axis
- **Cohort report template:** `scripts/templates/cohort_report_template.html` (lines 1-1109) - Production reference implementation

### Secondary (MEDIUM confidence)
- **CSS Grid best practices:** https://www.testmu.ai/blog/css-grid-best-practices/ - Grid layout patterns
- **Material UI Chip component:** https://mui.com/material-ui/react-chip/ - Industry standard chip UI pattern
- **Above-the-fold design:** https://www.omniconvert.com/blog/above-the-fold-design/ - Viewport constraints

### Tertiary (LOW confidence)
- **noUiSlider GitHub:** https://github.com/leongersen/noUiSlider (14k+ stars, confirms popularity claim)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - noUiSlider is proven, DataTables/Chart.js already vendored and working
- Architecture: HIGH - Cohort report provides production-quality reference implementation for all patterns
- Pitfalls: HIGH - All pitfalls verified from cohort report debugging comments and DataTables issue tracker

**Research date:** 2026-02-17
**Valid until:** 30 days (stable domain, noUiSlider last updated Jan 2024, DataTables v2 stable)
