# Phase 18: Foundation — Core Abstractions and Fisher Refactor - Context

**Gathered:** 2026-02-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Establish the `association/` package skeleton with core abstractions (AssociationTest ABC, TestResult, AssociationEngine, AssociationConfig, correction.py) and refactor the existing Fisher's exact gene burden test into this new framework. Users can run `--perform-association --association-tests fisher` and receive output bit-identical to `--perform-gene-burden`. The old gene burden path stays fully functional and untouched. No new statistical tests are added — only Fisher moves to the new framework.

</domain>

<decisions>
## Implementation Decisions

### CLI Interface Design
- Test selection via single flag with comma-separated values: `--association-tests fisher,burden,skat`
- `--perform-association` is required as a master switch (mirrors `--perform-gene-burden` pattern)
- `--association-tests` without `--perform-association` is an error
- `--perform-association` without `--association-tests` defaults to Fisher only
- `--perform-association` and `--perform-gene-burden` can coexist in the same run — both stages run independently, producing their own output (useful for validation during transition)

### Association Output Format
- Separate "Association" sheet in Excel alongside existing Gene Burden sheet (both present when both flags active)
- New column schema from the start, designed for multi-test output (not cloning Gene Burden columns)
- Wide format: one row per gene with test-specific columns (fisher_p, burden_p, skat_p, acat_o_p, effect sizes, CIs, counts) — consistent with Phase 22 per-gene TSV spec
- Phase 18 fills Fisher columns; other test columns added as phases deliver them

### Deprecation Strategy
- No deprecation warning for `--perform-gene-burden` in Phase 18 — wait until Phase 22+ when association framework is feature-complete
- When deprecation happens (Phase 22+), `--perform-gene-burden` becomes a thin wrapper that internally calls the association framework with Fisher test — one code path to maintain
- For Phase 18, both paths remain fully independent

### Fisher Code Approach
- **Clean reimplementation** in `association/tests/fisher.py` — NOT delegation to gene_burden.py internals
- The ~80 lines of statistical core (scipy.stats.fisher_exact, statsmodels.Table2x2.oddsratio_confint, smm.multipletests) are ported fresh to conform to the AssociationTest ABC interface
- Bit-identity guaranteed by test suite (cross-validation against gene_burden.py output), not by code sharing
- Correct coupling direction: `gene_burden.py` imports from `association/correction.py` (CORE-08), never the reverse
- Look for optimization/speedup opportunities during reimplementation (e.g., vectorized table construction, batch fisher_exact calls)
- Aggregation functions stay in `gene_burden.py` for now — FisherExactTest receives pre-built contingency data via the standard genotype matrix interface

### Error Messaging
- Requesting an unregistered test name (e.g., `--association-tests skat` before Phase 20) is a hard error with exit — "Test 'skat' is not available. Available tests: fisher"
- Missing dependencies (scipy, statsmodels) checked eagerly at startup when `--perform-association` is specified — fail before any processing begins
- Genes with zero variants passing filters: skip silently, report p_value=NA in output (no log noise)

### Logging
- Minimal at INFO level: summary only (e.g., "Association analysis: 342 genes tested, 12 significant (FDR < 0.05)")
- Per-gene detail available at DEBUG level for troubleshooting

### Claude's Discretion
- Exact AssociationTest ABC method signatures (as long as genotype matrix interface supports all planned tests)
- Internal structure of AssociationEngine (iteration strategy, result collection)
- TestResult dataclass field names and types
- AssociationConfig field organization
- Optimization strategies for Fisher reimplementation
- How AssociationAnalysisStage integrates with PipelineRunner dependency graph

</decisions>

<specifics>
## Specific Ideas

- The clean reimplementation approach was chosen after researching strangler fig / adapter patterns — delegation would invert the coupling direction (association/ -> gene_burden.py), blocking clean deprecation later
- The parity test suite (Phase 18 plan 18-04) serves as a permanent regression guard through all 5 subsequent phases
- Aggregation layer (three strategies in gene_burden.py) stays put for now — the ABC's genotype matrix interface separates aggregation from testing cleanly

</specifics>

<deferred>
## Deferred Ideas

- Deprecation warning for `--perform-gene-burden` — Phase 22+ (after association framework is feature-complete)
- Converting `--perform-gene-burden` to a thin wrapper — Phase 22+
- Moving aggregation functions out of `gene_burden.py` — future cleanup phase

</deferred>

---

*Phase: 18-foundation-core-abstractions-and-fisher-refactor*
*Context gathered: 2026-02-19*
