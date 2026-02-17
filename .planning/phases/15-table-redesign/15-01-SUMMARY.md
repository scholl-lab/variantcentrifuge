---
phase: 15-table-redesign
plan: 01
subsystem: html-report
tags: [datatables, css, vendor, fixedcolumns, table-styling]
requires:
  - phases: [13-js-stack-modernization, 14-information-hierarchy]
    reason: "Built on DataTables v2.2.2 and semantic color palette"
provides:
  - capability: "FixedColumns extension vendored"
  - capability: "Dark header gradient styling"
  - capability: "Zebra stripe row backgrounds"
  - capability: "Density mode CSS classes"
  - capability: "Detail panel grid layout CSS"
  - capability: "Column type CSS classes"
affects:
  - plan: "15-02"
    reason: "Plan 02 will wire up JS behavior using these CSS classes"
tech-stack:
  added:
    - "DataTables FixedColumns 5.0.5"
  patterns:
    - "CSS-only density modes (JS toggle in Plan 02)"
    - "Grid layout for detail panels"
    - "Shadow separator for sticky columns"
key-files:
  created:
    - "variantcentrifuge/assets/js/fixedcolumns.min.js"
    - "variantcentrifuge/assets/css/fixedcolumns.dataTables.min.css"
  modified:
    - "variantcentrifuge/templates/index.html"
decisions:
  - id: "15-01-01"
    decision: "Dark header gradient from #2c3e50 to #1a252f with white text"
    rationale: "Charcoal gradient provides strong visual hierarchy while aligning with Phase 14 color palette"
  - id: "15-01-02"
    decision: "Zebra striping with #f8f9fa (odd) and white (even)"
    rationale: "Subtle contrast improves row scannability without being visually distracting"
  - id: "15-01-03"
    decision: "Tooltip max-width reduced from 400px to 300px"
    rationale: "TABLE-01 requirement for narrower tooltips, better for longer HGVS notation"
  - id: "15-01-04"
    decision: "Truncated-cell max-width increased from 150px to 180px"
    rationale: "Better fit for HGVS notation while maintaining table compactness"
  - id: "15-01-05"
    decision: "FixedColumns JS loaded after datatables.min.js, before buttons"
    rationale: "Extension dependencies: FixedColumns requires DataTables core, should load before optional buttons"
duration: 157s
completed: 2026-02-17
---

# Phase 15 Plan 01: Vendor FixedColumns and Add CSS Styles Summary

**One-liner:** Vendored DataTables FixedColumns extension and added all CSS for table redesign including dark header gradient, zebra striping, density modes, and detail panel grid layout

## What Was Done

### Task 1: Vendor FixedColumns Extension Assets

Downloaded and vendored DataTables FixedColumns v5.0.5:
- **JS:** `fixedcolumns.min.js` (8KB) - Sticky column functionality
- **CSS:** `fixedcolumns.dataTables.min.css` (2KB) - FixedColumns styling for DataTables default theme

Files are automatically discovered by the existing `_load_assets()` function in `generate_html_report.py`, which scans the `assets/js/` and `assets/css/` directories and namespaces keys by subdirectory to prevent collisions.

### Task 2: Add CSS Styles and Asset References

Added comprehensive CSS for all Phase 15 table enhancements:

**Asset References:**
- Inlined FixedColumns CSS: `{{ assets['css/fixedcolumns.dataTables.min'] }}`
- Inlined FixedColumns JS: `{{ assets['js/fixedcolumns.min'] }}`

**New CSS Styles:**

1. **TABLE-05: Dark Header** - Gradient background (#2c3e50 → #1a252f) with white text and 3px blue bottom border
2. **TABLE-04: Zebra Striping** - Alternating row backgrounds (odd: #f8f9fa, even: white, hover: #e3f2fd)
3. **TABLE-07: Density Modes** - Three modes with padding variants:
   - Compact: 4px/6px (default)
   - Regular: 8px/10px
   - Relaxed: 12px/14px
4. **TABLE-03: Detail Panel** - Grid layout with auto-fit columns (min 280px), slideDown animation, card sections with headers
5. **TABLE-03: Chevron Control** - Centered control column (30px width), 90° rotation animation on expand
6. **TABLE-02: FixedColumns Shadow** - 2px box-shadow on frozen columns for "floating above" effect
7. **TABLE-06/08: Column Classes** - Monospace, dt-right, middle-truncate utilities
8. **TABLE-09: Record Count** - Prominent info styling (14px, bold, #212529)
9. **TABLE-01: Tooltip Updates** - max-width 300px (down from 400px), truncated-cell 180px (up from 150px)

**Preserved:**
- All existing CSS (badges, dashboard, skeleton, metadata footer, charts)
- Badge render functions unchanged
- Link styles unchanged

## Technical Details

### FixedColumns Integration

The extension is vendored following the same pattern as Phase 13:
- Downloaded from jsDelivr CDN (reliable, versioned URLs)
- Stored in `variantcentrifuge/assets/` subdirectories
- Auto-discovered by `_load_assets()` with subdirectory namespacing
- Inlined in template for true single-file HTML reports

**Load Order:**
```html
<script>{{ assets['js/jquery.slim.min'] }}</script>
<script>{{ assets['js/datatables.min'] }}</script>
<script>{{ assets['js/fixedcolumns.min'] }}</script>  <!-- After DataTables core -->
<script>{{ assets['js/datatables.buttons.min'] }}</script>
```

### CSS Architecture

All new styles are additive to existing Phase 13/14 CSS. The styles are organized by TABLE requirement ID for traceability:

- **TABLE-02:** FixedColumns visual separator
- **TABLE-03:** Expandable row details
- **TABLE-04:** Zebra striping
- **TABLE-05:** Dark header
- **TABLE-06/08:** Column type classes
- **TABLE-07:** Density modes
- **TABLE-09:** Record count prominence
- **TABLE-01:** Tooltip adjustments

### Density Mode Classes

CSS-only classes applied to `#variants_table`:
- `.density-compact` - Default, power user mode
- `.density-regular` - Balanced spacing
- `.density-relaxed` - Maximum whitespace

Plan 02 will add the JavaScript toggle button and localStorage persistence.

### Detail Panel Layout

The `.detail-panel` uses CSS Grid with `repeat(auto-fit, minmax(280px, 1fr))`:
- Responsive: 1-3 columns depending on available width
- 280px minimum column width (fits category card content)
- 16px gap between cards
- 3px blue left border for visual hierarchy
- 200ms slideDown animation for smooth reveal

Each `.detail-section` is a bordered card with:
- White background, subtle shadow
- 12px text-transform uppercase header
- `<dl>` grid layout (auto 1fr columns, 6px/12px gaps)
- Link styling consistent with existing link-icon pattern

## Verification

All verification criteria met:

1. **FixedColumns assets exist:**
   - `ls -la variantcentrifuge/assets/js/fixedcolumns.min.js` → 8095 bytes
   - `ls -la variantcentrifuge/assets/css/fixedcolumns.dataTables.min.css` → 1920 bytes

2. **Content verification:**
   - `grep -c "fixedColumns\|FixedColumns" fixedcolumns.min.js` → 2 matches
   - `grep -c "dtfc" fixedcolumns.dataTables.min.css` → 1 match

3. **Template verification:**
   - `grep -c "assets\['css/fixedcolumns.dataTables.min'\]"` → 1
   - `grep -c "assets\['js/fixedcolumns.min'\]"` → 1
   - `grep -c "density-compact"` → 2
   - `grep -c "detail-panel"` → 1
   - `grep -c "dt-control"` → 3
   - `grep -c "dtfc-fixed-left"` → 1
   - `grep -c "#2c3e50"` → 1 (dark header)
   - `grep -c "nth-child(odd)"` → 1 (zebra)

4. **Lint check:** `make lint` → All checks passed

5. **Asset key resolution:** Template contains valid Jinja2 asset references that will be resolved by `_load_assets()` at HTML generation time

## Deviations from Plan

None - plan executed exactly as written.

## Next Phase Readiness

**Blockers:** None

**Dependencies satisfied:**
- FixedColumns extension is now available for Plan 02 JavaScript initialization
- All CSS classes are defined for Plan 02 to reference in DataTables config
- Density mode classes ready for JS toggle button
- Detail panel CSS ready for child row rendering
- Chevron control CSS ready for click handlers

**Plan 02 can proceed** with wiring up:
- FixedColumns initialization (`fixedColumns: { left: 1 }`)
- Density toggle button with localStorage
- Child row expand/collapse handlers
- Column width `columnDefs` rules
- Middle truncation render functions

## Commits

| Task | Commit | Message |
|------|--------|---------|
| 1 | 7a831a6 | chore(15-01): vendor DataTables FixedColumns v5.0.5 extension |
| 2 | 4beb184 | style(15-01): add Phase 15 CSS styles to template |

**Total:** 2 commits (one per task)

## Performance Impact

**Bundle size increase:**
- FixedColumns JS: +8KB
- FixedColumns CSS: +2KB
- Template CSS additions: ~3KB (estimated from line count)
- **Total:** ~13KB increase to HTML report size

This is acceptable given the functionality enabled (sticky columns, expandable details, density modes). For reference, Phase 13 reduced bundle size by ~3.3MB (Plotly → Chart.js), so Phase 15's 13KB increase is negligible.

## Lessons Learned

1. **Subdirectory namespacing is critical:** The existing `_load_assets()` pattern of namespacing keys by subdirectory (`css/fixedcolumns.dataTables.min` vs `js/fixedcolumns.min`) prevents collisions when CSS and JS have similar base names.

2. **CSS-first approach enables progressive enhancement:** By defining all CSS classes in Plan 01, Plan 02 can focus purely on JavaScript behavior without template modifications. This separation of concerns makes testing easier.

3. **Gradient headers need explicit color for both clone and original:** DataTables creates both `.dt-scroll-head` and `.dataTables_scrollHead` elements depending on version/config, so all header selectors must target both variants.

4. **Zebra striping requires child row exception:** The `tr.dt-hasChild + tr > td` rule prevents child rows (detail panels) from inheriting zebra stripe colors, keeping them visually distinct.

5. **Loading order matters for extensions:** FixedColumns must load after DataTables core but can load before or after Buttons. The current order (FixedColumns before Buttons) follows dependency hierarchy.

---

*Plan: 15-01*
*Completed: 2026-02-17*
*Duration: 157 seconds (2m 37s)*
