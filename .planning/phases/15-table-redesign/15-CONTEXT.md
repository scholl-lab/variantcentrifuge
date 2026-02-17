# Phase 15: Table Redesign - Context

**Gathered:** 2026-02-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Transform the variant data table into a modern enterprise-grade table with Tippy.js tooltips on truncated cells, sticky GENE column, expandable row detail panels, configurable density modes, and intelligent column widths. All changes target the individual HTML report's embedded DataTables v2 table. Phase 13 (modern JS stack) and Phase 14 (semantic color badges) are prerequisites.

</domain>

<decisions>
## Implementation Decisions

### Expandable row details
- Multi-column grid layout for the detail panel (2-3 columns depending on space)
- Fields grouped by category (Identifiers, Annotations, Scores, Links), each in a bordered card with a header
- Smooth slide animation (~200ms) for expand/collapse transitions
- Links (ClinVar, gnomAD, Varsome, etc.) render as clickable URLs with external-link icons

### Column behavior & widths
- End truncation for HGVS notation (standard ellipsis at end, e.g., `NM_000350.3:c.5882G...`)
- Monospace font for technical columns only: CHROM, POS, REF, ALT, HGVS
- Right-align numeric score columns (CADD, SpliceAI, allele frequency, etc.)
- Fixed widths for CHROM/POS, grow for GENE, truncation for HGVS
- "Showing X of Y variants" displayed prominently, 25 rows per page default, zebra striping

### Header & density styling
- Dark header row with white text and a static accent-colored bottom border (purely decorative)
- Default density mode: Compact (power users scanning many variants)
- Three density modes: Compact, Regular, Relaxed — persisted via localStorage
- Toolbar positioned above table, right-aligned — contains density toggle, column visibility, row count

### Tooltip behavior
- Immediate trigger on hover (no delay) for truncated cells
- Content formatted with line breaks at logical separators (semicolons, pipes) for readability
- Maximum tooltip width: 300px (narrow)
- Touch devices: tap to show, tap elsewhere to dismiss
- Keyboard accessible via focus, dismissable with Escape
- Old yellow hover-expand completely removed

### Claude's Discretion
- Dark header exact color (slate/charcoal vs dark blue — based on Phase 14 palette)
- Sticky GENE column visual separator (shadow preferred based on UX research — validate against DataTables v2)
- Detail panel grid column count (2 vs 3 based on available space)
- Exact density mode padding values
- Zebra stripe color intensity
- Column visibility defaults beyond what Phase 13 established

</decisions>

<specifics>
## Specific Ideas

- Sticky column separator: UX best practice research indicates drop shadow creates "floating above" illusion for spatial depth — used by Google Sheets, Airtable, Linear. Preferred over border unless DataTables v2 constraints dictate otherwise.
- Category cards in expandable detail should feel visually consistent with Phase 14's dashboard cards
- Compact mode is the default because clinical geneticists typically scan many variants quickly

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 15-table-redesign*
*Context gathered: 2026-02-17*
