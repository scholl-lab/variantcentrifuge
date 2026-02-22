---
phase: 26-association-testing-documentation
plan: 02
subsystem: documentation
tags: [docs, association-testing, api-reference, changelog, sphinx, myst]

# Dependency graph
requires:
  - phase: 26-01
    provides: association_testing.md guide that cross-references point to
  - phase: 18-25
    provides: implemented association testing framework (AssociationConfig, AssociationEngine, all tests)
provides:
  - API reference stub for variantcentrifuge.association package
  - v0.15.0 changelog entry following Keep a Changelog format
  - Association testing cross-references in 5 existing doc files
  - Association testing feature bullet in README
affects:
  - 27-performance-optimizations (will need to update changelog/docs when adding optimizations)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "API stub follows gene_burden.md pattern: title + description + automodule directive"
    - "Cross-references use brief mention + link, never duplicate content from target guide"

key-files:
  created:
    - docs/source/api/association.md
    - .planning/phases/26-association-testing-documentation/26-02-SUMMARY.md
  modified:
    - docs/source/api/index.md
    - docs/source/usage.md
    - docs/source/guides/cohort_analysis.md
    - docs/source/faq.md
    - docs/source/index.md
    - docs/source/changelog.md
    - README.md

key-decisions:
  - "API stub covers public exports from __init__.py via single automodule directive (not per-submodule stubs)"
  - "FAQ entries explain ACAT-O semantics and R deprecation in self-contained answers with links to guide"
  - "index.md toctree updated to include guides/association_testing (discoverability from main nav)"

patterns-established:
  - "New features get: API stub + feature bullet in index.md + CLI flags table in usage.md + FAQ entries + changelog entry"

# Metrics
duration: 3min
completed: 2026-02-22
---

# Phase 26 Plan 02: Association Testing Docs - Cross-references and API Stub Summary

**API reference stub, v0.15.0 changelog, and association testing cross-references added to all 7 doc touchpoints (usage.md, cohort_analysis.md, faq.md, index.md, api/index.md, changelog.md, README.md)**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-22T19:58:21Z
- **Completed:** 2026-02-22T20:01:00Z
- **Tasks:** 3
- **Files modified:** 8

## Accomplishments
- Created `docs/source/api/association.md` API stub following gene_burden.md pattern with automodule directive
- Updated API index toctree and Module Overview to include association module
- Added 15-row CLI flag table for association testing flags to usage.md with link to guide
- Added association testing pointer paragraph to cohort_analysis.md after Gene Burden section
- Added 3 FAQ entries covering covariate adjustment, ACAT-O semantics, and R deprecation
- Updated index.md feature description bullet and toctree to include association testing
- Added v0.15.0 changelog entry with all 8 feature groups in Keep a Changelog format
- Added modular association testing feature bullet to README.md

## Task Commits

Each task was committed atomically:

1. **Task 1: Create API reference stub and update API index** - `eda4abd` (feat)
2. **Task 2: Update existing docs with cross-references** - `782e894` (docs)
3. **Task 3: Add changelog entry and update README** - `515e2a8` (docs)

## Files Created/Modified
- `docs/source/api/association.md` - New API stub with automodule:: variantcentrifuge.association
- `docs/source/api/index.md` - Added association to toctree and Module Overview
- `docs/source/usage.md` - Added Association Testing subsection with 15-row CLI flag table
- `docs/source/guides/cohort_analysis.md` - Added advanced association testing pointer paragraph
- `docs/source/faq.md` - Added Association Testing section with 3 FAQ entries
- `docs/source/index.md` - Updated feature bullet; added association_testing to toctree
- `docs/source/changelog.md` - Added [0.15.0] entry with 8 feature groups
- `README.md` - Added modular association testing feature bullet

## Decisions Made
- API stub uses single `automodule` directive (not per-submodule stubs) because `__init__.py` re-exports all public classes (AssociationConfig, AssociationEngine, AssociationTest, TestResult, apply_correction)
- FAQ entries are self-contained answers with links to the guide, not forwarding stubs, so they serve users who land on FAQ directly
- `index.md` main toctree updated to include `guides/association_testing` (in addition to guides/index.md) because the main nav is built from index.md's toctree list

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Phase 26 documentation is now complete (Plan 01: guide; Plan 02: cross-references and stubs)
- Plan 01 (association_testing.md guide) was found already committed before this plan ran
- All cross-references in this plan point to docs/source/guides/association_testing.md which exists
- Ready for Phase 27: Association Performance Optimizations

---
*Phase: 26-association-testing-documentation*
*Completed: 2026-02-22*
