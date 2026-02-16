---
phase: 13-js-stack-modernization
plan: 02
subsystem: ui
tags: [html, jinja2, datatables-v2, chartjs, tippyjs, vanilla-js]

# Dependency graph
requires:
  - phase: 13-01
    provides: Vendored JS/CSS assets and asset loading infrastructure
provides:
  - Fully modernized HTML report template with inlined assets
  - DataTables v2 vanilla JS initialization
  - Chart.js horizontal bar chart with semantic colors
  - Tippy.js tooltips for truncated content
  - Loading skeleton with shimmer animation
affects: [14-information-hierarchy, 15-table-redesign, 16-column-filtering]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Jinja2 asset inlining via {{ assets['js/library'] }}"
    - "DataTables v2 vanilla JS API (new DataTable, not $.DataTable)"
    - "Chart.js horizontal bar with semantic color mapping"
    - "Tippy.js tooltip attachment in drawCallback"
    - "Loading skeleton with CSS shimmer animation"

key-files:
  created:
    - tests/unit/test_html_report_assets.py
  modified:
    - variantcentrifuge/templates/index.html
    - variantcentrifuge/generate_html_report.py

key-decisions:
  - "Fixed asset key collision by namespacing with subdirectory prefix (js/, css/)"
  - "Semantic color scheme: HIGH=red, MODERATE=orange, LOW=amber, MODIFIER=gray"
  - "Tippy.js tooltips only attach when content is actually truncated"
  - "Loading skeleton shows 7 shimmer rows during DataTable initialization"

patterns-established:
  - "Asset inlining: All JS/CSS embedded in single HTML file via Jinja2 {{ assets['...'] }}"
  - "No jQuery user code: jQuery present only for DataTables internal use"
  - "Vanilla JS event handling: document.addEventListener, not $().ready()"
  - "DataTables v2 layout API for button/pagination placement"

# Metrics
duration: 35min
completed: 2026-02-16
---

# Phase 13 Plan 02: Template Rewrite Summary

**HTML report template fully modernized: DataTables v2 with vanilla JS, Chart.js horizontal bar with semantic colors, Tippy.js tooltips, loading skeleton, zero CDN dependencies**

## Performance

- **Duration:** 35 min
- **Started:** 2026-02-16T17:15:00Z
- **Completed:** 2026-02-16T17:50:00Z
- **Tasks:** 1 (plus 1 blocking fix)
- **Files modified:** 3

## Accomplishments

- Complete HTML template rewrite with all assets inlined via Jinja2
- DataTables v2 initialization using vanilla JS API (new DataTable)
- Chart.js horizontal bar chart with semantic impact colors replacing Plotly
- Tippy.js tooltips replacing old yellow hover-expand CSS popup
- Loading skeleton with shimmer animation for better perceived performance
- Zero jQuery user code (jQuery only present for DataTables internal dependency)
- Zero external CDN links (all libraries embedded in single HTML file)

## Task Commits

Each task was committed atomically:

1. **Blocking fix: Namespace asset keys** - `431c3d5` (fix)
2. **Task 1: Rewrite template** - `488ef8b` (feat)

## Files Created/Modified

- `variantcentrifuge/generate_html_report.py` - Fixed _load_assets() to namespace keys by subdirectory (js/, css/)
- `variantcentrifuge/templates/index.html` - Complete rewrite with modern JS stack
- `tests/unit/test_html_report_assets.py` - New test file for asset loading and template rendering

## Decisions Made

1. **Asset key namespacing:** Changed from `f.stem` to `f"{subdir}/{f.stem}"` to prevent key collisions between same-stem files (datatables.min.js vs datatables.min.css)

2. **Semantic color scheme for impact levels:**
   - HIGH: #dc3545 (red)
   - MODERATE: #fd7e14 (orange)
   - LOW: #ffc107 (amber)
   - MODIFIER: #6c757d (gray)

3. **Conditional Tippy.js attachment:** Only attach tooltips to cells where content is actually truncated (offsetWidth < scrollWidth), not all cells with `apply_hover_expand` flag

4. **Loading skeleton timing:** Show 7 shimmer rows during table initialization, hide in DataTable initComplete callback

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Asset key collision between JS and CSS files**
- **Found during:** Task 1 (Template rewrite planning)
- **Issue:** _load_assets() used f.stem as keys, causing datatables.min.js and datatables.min.css to collide (CSS overwrites JS since css/ iterates second)
- **Fix:** Changed asset keys to namespace by subdirectory: `js/datatables.min` and `css/datatables.min`
- **Files modified:** variantcentrifuge/generate_html_report.py
- **Verification:** Ran _load_assets() and confirmed all 10 assets have unique keys with proper namespacing
- **Committed in:** 431c3d5 (separate blocking fix commit before template work)

---

**Total deviations:** 1 auto-fixed (Rule 3 blocking)
**Impact on plan:** Essential fix to enable template to load correct assets. Without namespacing, template would embed CSS content where JS was expected.

## Issues Encountered

None - template rewrite proceeded smoothly after blocking fix.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 14 (Information Hierarchy and Semantic Color Coding):**
- Modern JS stack foundation in place
- All assets inlined and working
- DataTables v2, Chart.js, Tippy.js functional
- Template renders valid HTML with ~498KB single-file output

**Blockers/Concerns:** None

**Note for Phase 14:** The semantic color scheme established here (HIGH=red, MODERATE=orange, etc.) should be used consistently for impact badges in Phase 14.

---
*Phase: 13-js-stack-modernization*
*Completed: 2026-02-16*
