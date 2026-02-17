---
phase: 17-accessibility-and-print-pdf
plan: 01
subsystem: ui
tags: [wcag, accessibility, aria, svg, contrast, a11y, html]

# Dependency graph
requires:
  - phase: 15-table-redesign
    provides: DataTables v2 initialization, Tippy.js tooltips, badge rendering, table structure
provides:
  - Skip-to-content link for keyboard navigation
  - ARIA roles and labels on all interactive elements
  - Accessible SVG link icons with screen reader text
  - WCAG AA 4.5:1 contrast-compliant badge colors
  - Keyboard-accessible tooltip triggers
affects: [17-02, 17-03, print-stylesheet, pdf-generation]

# Tech tracking
tech-stack:
  added: []
  patterns: [skip-link pattern, sr-only utility class, SVG sprite pattern, WCAG color compliance]

key-files:
  created: []
  modified: [variantcentrifuge/templates/index.html]

key-decisions:
  - "Use role='table' not 'grid' for DataTables (research finding: grid implies editable)"
  - "SVG icons with aria-hidden='true' plus sr-only text for screen readers"
  - "All badge colors verified at 4.5:1 contrast ratio using WebAIM contrast checker"
  - "Truncated cells get tabindex='0' for keyboard focus to trigger Tippy.js tooltips"

patterns-established:
  - "sr-only utility class for screen-reader-only text (position: absolute off-screen)"
  - "skip-link pattern: hidden by default, visible on :focus with high z-index"
  - "SVG sprite with defs/symbol pattern for reusable icons"
  - "ARIA attributes added in DataTables initComplete callback after DOM ready"

# Metrics
duration: 5min
completed: 2026-02-17
---

# Phase 17 Plan 01: Core Accessibility Features Summary

**WCAG 2.1 AA foundation: skip-link, ARIA roles, SVG icons, 4.5:1 contrast badges, keyboard tooltips**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-17T14:22:45Z
- **Completed:** 2026-02-17T14:27:34Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments
- Skip-to-content link as first focusable element enables keyboard users to bypass header
- All DataTables controls (search, page length, colvis, density, info) have proper ARIA labels
- Emoji link icons replaced with accessible SVG icons plus screen-reader-only text
- All badge colors (IMPACT, ClinVar, Inheritance) pass WCAG AA 4.5:1 contrast with white text
- Truncated cells keyboard-accessible via tabindex="0" for tooltip interaction

## Task Commits

Each task was committed atomically:

1. **Task 1: Add skip-link, ARIA roles, and keyboard tooltip support** - `d319435` (feat)
2. **Task 2: Replace emoji icons with SVG and fix badge contrast** - `a456d86` (feat)

## Files Created/Modified
- `variantcentrifuge/templates/index.html` - Added skip-link, ARIA roles, SVG sprite, contrast-fixed badge colors, keyboard-accessible tooltips

## Decisions Made

**1. Use role="table" not role="grid" for DataTables**
- Research finding from 17-RESEARCH.md: ARIA grid implies editable cells
- DataTables is read-only, so role="table" is semantically correct
- Applied in initComplete callback after table initialization

**2. SVG icons with dual accessibility approach**
- SVG marked `aria-hidden="true"` so screen readers ignore visual element
- Adjacent `.sr-only` span provides descriptive text for screen readers
- Pattern: `<svg aria-hidden>...</svg><span class="sr-only">IGV report (opens in new tab)</span>`

**3. Badge color contrast fixes validated**
- Used WebAIM Contrast Checker to verify all colors meet WCAG AA 4.5:1 ratio
- MODERATE: #fd7e14 (2.55:1 FAIL) → #c05000 (4.6:1 PASS)
- LOW: #f59e0b (2.17:1 FAIL) → #b45309 (5.2:1 PASS)
- Similar fixes for ClinVar and Inheritance badges (see commit details)
- Colors remain visually distinct while meeting accessibility requirements

**4. Truncated cells keyboard-accessible via tabindex="0"**
- Tippy.js already configured with `trigger: 'mouseenter focus'` from Phase 15
- Adding `tabindex="0"` makes truncated cells focusable by keyboard
- Users can Tab to truncated content and trigger tooltip on focus
- Applied to Jinja template and JavaScript render functions

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all accessibility features integrated cleanly with existing Phase 15 infrastructure.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

Core accessibility foundation complete. Ready for Phase 17-02 (Focus Management) and 17-03 (Print Stylesheet).

**Dependencies satisfied:**
- Skip-link present and functional
- ARIA roles established on table and controls
- SVG icons ready for print stylesheet optimization
- Badge colors won't need adjustment for print (already high contrast)

**No blockers for:**
- Focus management enhancements (17-02)
- Print/PDF stylesheet work (17-03)
- Screen reader testing validation

---
*Phase: 17-accessibility-and-print-pdf*
*Completed: 2026-02-17*
