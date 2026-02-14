# Phase 9: Inheritance Analysis Optimization - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

10-100x speedup on inheritance analysis through vectorization of the three-pass analysis (deduction, compound het, prioritization), which currently consumes 40-60% of total pipeline time. Also includes making comp_het_vectorized.py the sole implementation, refactoring parallel_analyzer.py to use new vectorized code, and removing the original non-vectorized implementations.

</domain>

<decisions>
## Implementation Decisions

### Correctness validation
- Clinically equivalent output required (not byte-identical): same inheritance patterns and prioritization, row order may differ
- When sorted by a stable key, all values must match exactly — row order is the only acceptable difference
- Golden file comparison as the validation strategy: generate reference outputs from current code, diff against vectorized output
- Golden files from both synthetic data (Phase 6 generators) and real anonymized clinical data
- User will provide a large anonymized VCF file (~5000 cases) with phenotype and gene set groupings for validation
- Golden file comparison runs as a separate validation script, not part of the normal pytest suite
- Any difference in vectorized output is a build failure — fix the vectorized code, no tolerance for "acceptable" discrepancies
- Edge cases (missing genotypes, unusual ploidy) must produce identical results or the build fails

### Rollback strategy
- Replace original implementations in-place once validated — no dual code paths, no feature flags
- Vectorize pass-by-pass: Pass 1 (deduction) → validate → commit, then Pass 2 (comp het) → validate → commit, then Pass 3 (prioritization) → validate → commit
- Each pass is independently revertable via git history
- Partial vectorization is acceptable if a specific pass has sections too complex to vectorize fully — accept partial speedup rather than blocking
- Old code preserved in git history only — no archive files, no legacy docs

### Vectorization approach
- Go directly to full NumPy vectorization for Pass 1 — skip itertuples intermediate step
- Research best practices for vectorizing complex conditional logic (de novo, AD, AR, X-linked, mitochondrial pattern branches) — balance optimization with maintainability as a senior systems engineer would
- Readability first: clear variable names, comments explaining genetic logic, step-by-step vectorized operations — even if slightly slower than maximally dense NumPy
- Apply KISS, DRY, SOLID principles — modular, maintainable code that works correctly
- Edge case handling strategy to be determined during research — research best practices rather than prescribing scalar fallback vs full vectorization upfront

### Compound het handling
- Make comp_het_vectorized.py the sole implementation, remove original comp_het.py after validation
- Gene-level grouping via groupby('gene').apply() with vectorized inner logic — safer, easier to debug
- Refactor comp_het_vectorized.py to use Phase 7-8 patterns (categorical dtypes, observed=True, optimized groupby)
- Update parallel_analyzer.py in this phase to use the new vectorized implementations — it's naturally coupled to the vectorization work

### Claude's Discretion
- Specific NumPy vectorization patterns for each inheritance rule (boolean masks, np.select, matrix ops — research and decide)
- How to structure the validation script (CLI interface, output format, reporting)
- Optimal ordering of boolean mask evaluations for performance
- How to decompose deduce_patterns_for_variant into vectorized sub-operations
- Whether to use pandas or pure NumPy for specific operations

</decisions>

<specifics>
## Specific Ideas

- User wants a systems engineering approach: research best practices for vectorizing genomic analysis pipelines, not ad-hoc optimization
- Large-scale validation with ~5000 real clinical cases provides high confidence before removing fallbacks
- Readability is valued over micro-optimizations — this is clinical software where correctness and auditability matter
- The pass-by-pass approach allows incremental validation and reduces risk of large-scale breakage

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 09-inheritance-analysis-optimization*
*Context gathered: 2026-02-14*
