# Phase 14: Information Hierarchy and Semantic Color Coding - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Users see a clinically meaningful overview (summary dashboard with colored cards, charts) before the variant table, with semantic color coding applied throughout. The dashboard includes summary cards, impact/inheritance charts, top genes, and a metadata footer. Color badges are added to IMPACT, ClinVar, and inheritance columns in the table.

This phase does NOT include: column-level filtering (Phase 16), table redesign/tooltips/expandable rows (Phase 15), or accessibility/print (Phase 17).

</domain>

<decisions>
## Implementation Decisions

### Summary dashboard layout
- Single row of compact metric cards at top, charts below in a second row
- All summary content (cards + charts) must fit strictly above the fold on a 1080p display (~600px vertical budget)
- Cards use neutral style (white/gray background) with colored numbers matching severity palette
- Top genes displayed as a ranked list in a card (top 5-10 genes with variant counts), not a chart
- Chart row: impact distribution chart (existing, enhanced) + inheritance distribution chart side by side

### Badge design in table
- Filled pill badges with white text on colored background (like GitHub labels)
- All badge types (IMPACT, ClinVar, inheritance) share the same visual style — same pill shape/size, differentiated only by color palette
- IMPACT colors: HIGH=#dc3545 (red), MODERATE=#fd7e14 (orange), LOW=#ffc107 (amber), MODIFIER=#6c757d (gray)
- ClinVar colors: Pathogenic=red, Likely Pathogenic=orange-red, VUS=yellow/amber, Likely Benign=light green, Benign=green
- Inheritance colors: de novo=red, compound het=purple, AD=blue, AR=green, X-linked=teal
- Unknown/unresolved inheritance: gray badge labeled "Unknown"

### Metadata footer
- Always-visible compact bar at bottom of report (not collapsible)
- Subtle gray background, small text — reference material, not a focal point
- Content limited to success criteria essentials: filter criteria, VCF source, reference genome, pipeline version, run date
- Filter criteria shown as actual filter expression (the SnpSift filter string), not a human-readable summary
- No analysis parameters (gene list, scoring model) — keep it minimal

### Data completeness handling
- When inheritance analysis wasn't run: show placeholder card saying "Inheritance analysis not available (no pedigree)" — signals the capability exists
- When ClinVar data is missing for a variant: empty cell (no badge, no text)
- Zero counts: always show (e.g., "HIGH: 0") — explicitly confirms categories were checked
- Trivial charts (single category): still render the chart — consistency over cleverness

### Claude's Discretion
- Exact spacing, padding, and typography within the dashboard grid
- Chart dimensions and aspect ratios (within the above-fold constraint)
- Card border/shadow styling details
- Exact amber shade for VUS badge (as long as it's distinct from pathogenic red and benign green)
- How to handle very long filter expressions in the footer (truncation/wrapping)

</decisions>

<specifics>
## Specific Ideas

- Current report already has the semantic color map in Chart.js (HIGH=red, MODERATE=orange, LOW=amber, MODIFIER=gray) — extend this to badges
- Current summary has 2 stat cards (Total Variants, Unique Genes) — expand to include sample count and filter criteria cards
- The impact distribution chart already exists — enhance it and add inheritance distribution alongside

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 14-information-hierarchy-and-semantic-color-coding*
*Context gathered: 2026-02-16*
