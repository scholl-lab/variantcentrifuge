# Phase 13: JS Stack Modernization - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace the current jQuery + DataTables v1.10 + Plotly stack with vanilla JS + DataTables v2 + Chart.js + Tippy.js. The individual HTML report must retain all existing functionality (sorting, pagination, search, column visibility, horizontal scroll) on the new stack. This is the foundation phase — all Phase 14-17 UX work depends on this migration.

</domain>

<decisions>
## Implementation Decisions

### Library loading strategy
- Inline everything: all JS and CSS embedded directly in the HTML file (no CDN links)
- Reports must be fully self-contained single-file HTML — works offline in clinical/HPC environments
- Minified library files shipped in the Python package (`variantcentrifuge/assets/` or similar), version-locked
- CSS also inlined — true zero-external-dependency single file
- Data (variants, summary) also embedded in the HTML — no separate JSON files needed

### Chart appearance
- Horizontal bar chart (not vertical like current Plotly chart)
- Semantic colors applied immediately: HIGH=red, MODERATE=orange, LOW=amber, MODIFIER=gray
- Hover tooltips showing count values (no click-to-filter behavior)
- On-bar value labels — count numbers displayed directly on each bar for immediate readability
- No click interactivity (filtering is Phase 16 scope)

### Loading experience
- Skeleton screen for the table area only (header and chart render immediately)
- 5-8 skeleton rows mimicking table layout
- Shimmer/pulse animation on skeleton elements
- Skeleton replaced by actual table once DataTables initialization completes

### Tooltip migration
- Migrate hover-expand cells to Tippy.js immediately (not temporary vanilla JS)
- Standard tooltip appearance: dark background, white text (replacing current yellow hover-expand)
- Tippy.js handles all truncated cell hover behavior from Phase 13 onward

### Column visibility
- Use DataTables v2 built-in column visibility buttons (direct migration from v1 Buttons extension)
- No custom dropdown — leverage DT v2 native support

### Header styling
- Light refresh of the page header: modern font, cleaner spacing, keep the blue theme
- Not a full redesign — just modernized to match the new component feel

### Claude's Discretion
- Exact asset directory structure and embedding mechanism
- DataTables v2 configuration details and API migration
- Tippy.js configuration (placement, delay, max-width)
- Skeleton screen CSS implementation
- Chart.js responsive behavior and sizing
- How to structure the single-file HTML template (script tag placement, etc.)

</decisions>

<specifics>
## Specific Ideas

- Single-file HTML is key — clinicians share reports via email and USB drives, can't rely on internet or folder structures
- The horizontal bar chart with on-bar labels should feel like a clean dashboard element, not a legacy data viz
- Loading skeleton should feel modern (think GitHub's loading states)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 13-js-stack-modernization*
*Context gathered: 2026-02-16*
