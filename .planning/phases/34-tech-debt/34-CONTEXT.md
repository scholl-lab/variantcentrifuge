# Phase 34: Tech Debt - Context

**Gathered:** 2026-02-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Fix association subsystem config mapping, validate COAST against R golden values, standardize association column naming, and document Fisher lambda_GC exclusion scope. All changes are internal correctness/consistency fixes — no new user-facing features.

</domain>

<decisions>
## Implementation Decisions

### COAST golden values
- One-time validation exercise, NOT a permanent CI regression test
- Run R AllelicSeries package on GCKD testing data (testing/ directory) and compare against Python COAST output
- Cover core + edge case scenarios: all 3 categories present, 1-2 categories missing, single-variant gene, all variants in one category, tied scores
- Fix any discrepancies found between R and Python implementations
- R script does not need to be committed — this is an alignment pass, not ongoing infrastructure

### Config mapping fix
- `create_stages_from_config()` should warn and skip on unrecognized config keys (not fail)
- Validate that required stage dependencies are met when activating stages from config — fail early if not
- Integration test only (no unit tests): set config dict keys and verify correct stages activate end-to-end
- Add docstring to the function with a mapping table showing config key → stage(s)

### Column naming standardization
- Audit ALL association output column names, not just SKAT-O
- Convention: `test_name_pvalue` (no underscore in "pvalue") — e.g., `burden_pvalue`, `fisher_pvalue`, `skat_o_pvalue`, `coast_pvalue`
- Same pattern for derived columns: `acat_o_qvalue`, `fdr_weight`, etc.
- Clean break — no backward compatibility aliases or deprecation period
- Q-value and all derived columns also follow the standardized pattern

### Lambda_GC documentation
- Code comments only (maintainer-facing, not user docs)
- Key points to document:
  - Fisher's exact test is exempt from lambda_GC correction (exact test, not asymptotic — applying GC would be statistically invalid)
  - Score-based tests (SKAT, burden) get lambda_GC correction
  - lambda_GC > 1.2 = significant inflation, investigate population stratification
  - Over-correction risk: GC correction in large studies can reduce power (cite 2025 literature)
- Do NOT include full threshold table — keep comments concise

</decisions>

<specifics>
## Specific Ideas

- COAST validation should use real GCKD VCF data for broader, realistic input — not just synthetic datasets
- Lambda_GC comment should reference: Fisher exempt (exact test), score-tests corrected, >1.2 investigate, over-correction loses real associations
- Column naming follows `test_pvalue` / `test_qvalue` pattern consistently across all association test types

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 34-tech-debt*
*Context gathered: 2026-02-24*
