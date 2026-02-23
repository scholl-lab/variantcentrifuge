# Phase 32: Region Restriction and PCA Wiring - Context

**Gathered:** 2026-02-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Add BED-based region prefilter restricting variant extraction to callable regions, and wire PCA computation into the pipeline so covariates flow automatically to association tests. Both features are standalone and independent of each other. Existing CLI and output formats remain unchanged.

</domain>

<decisions>
## Implementation Decisions

### BED mismatch handling
- Hard fail with clear error on chromosome naming mismatches (e.g., 'chr1' vs '1') — no auto-correction
- Strict BED validation — fail on malformed BED files (overlapping intervals, unsorted, missing columns)
- No coverage check — trust that the user knows their capture kit; just apply regions as given
- Zero-variant genes after restriction: log at DEBUG level (silent unless verbose)

### Region restriction scope
- Global prefilter — variants outside BED regions never enter the pipeline at all; every report reflects restricted regions
- Intersection happens once upstream, before chunk splitting — all gene chunks use the same restricted BED
- Log restriction summary at INFO level (e.g., "Region restriction: 1,234 of 5,678 gene regions retained")
- No metadata in output files — output format unchanged for backward compatibility

### PCA covariate flow
- Default 10 PCs as covariates (standard GWAS convention)
- Single `--pca` flag: auto-detect whether value is a tool name ("akt") or a file path (file exists → load pre-computed eigenvectors)
- Accept external eigenvector files from both akt and PLINK — infer format from filename/content
- If PCA computation fails (too few samples, singular matrix), fail the pipeline — don't silently run tests without covariates
- Cache computed eigenvectors using the existing pipeline cache system for reuse across runs

### CLI flag design
- `--regions-bed` for region restriction (takes BED file path)
- `--pca` for PCA covariates (takes tool name or eigenvector file path)
- `--pca-components N` to override default 10 PCs
- Both flags fully independent — can use either or both, no coupling
- Unify with any existing PCA-related flags (check `--pca-tool` from current CLI)

### Claude's Discretion
- BED intersection implementation (bedtools vs bcftools -R vs in-pipeline)
- PCA tool detection heuristic (how to distinguish "akt" from a file path)
- Eigenvector file format auto-detection logic
- Cache invalidation strategy for PCA eigenvectors
- Exact error messages and log formatting

</decisions>

<specifics>
## Specific Ideas

- User wants `--pca` to be smart about its argument: if it looks like a file, load it; if it's a tool name like "akt", compute. Infer format from filename where possible.
- Cache PCA results using the existing cache infrastructure — no new caching mechanism.
- Keep the experience simple: one flag per feature, sensible defaults, minimal required configuration.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 32-region-restriction-and-pca-wiring*
*Context gathered: 2026-02-23*
