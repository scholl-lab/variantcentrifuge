# Phase 7: Quick Wins - Tier 1 - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Immediate 30-40% speedup through zero-risk standard optimizations: dead code removal, groupby parameter fixes, temp file cleanup, and inter-stage memory management. No new features, no architectural changes. Phase 6 benchmarks prove improvement.

</domain>

<decisions>
## Implementation Decisions

### Dead code handling
- Verify the GT parsing loop (gene_burden.py:220-249) is truly dead before removing — add logging/assertion to confirm it never executes, then remove in the same PR
- Regression test validates semantically identical output (same data values, formatting differences acceptable)
- Clean removal — no comments about what was removed, git history is the record
- If verification reveals the loop IS reached in some edge case, refactor it inline as part of this phase rather than deferring

### groupby migration
- Update all 12+ groupby call sites to use `observed=True` in a single atomic commit
- Correctness first: if `observed=True` changes behavior at any call site (drops empty groups that affect output), skip that site and mark it for Phase 8 when categorical dtypes are added
- Add a Ruff custom lint rule that flags `groupby()` without `observed=True` — configured to block CI (treat as error)

### Memory management
- Place `gc.collect()` after every pipeline stage, not just heavy ones — maximum memory recovery with predictable profile
- Debug-level logging of memory before/after `gc.collect()` calls — visible with `-v` flag, silent otherwise
- Audit all temp files across the codebase for similar leaks, not just the known gene_bed.py leak — fix all found in this phase
- Refactor temp file handling to use context managers (`tempfile.NamedTemporaryFile` with `with` statements) instead of explicit cleanup

### Validation strategy
- Run full Phase 6 benchmark suite before/after — not just gene burden, ensure no regressions anywhere
- If combined optimizations yield less than 30% speedup, investigate further for additional quick wins before moving to Phase 8
- Update `.planning/performance-analysis-report.md` with before/after benchmark numbers for the performance record
- Benchmark on local dev machine — consistent same-machine comparison, absolute numbers may vary

### Claude's Discretion
- Exact Ruff rule configuration for groupby enforcement
- gc.collect() logging format and memory measurement approach
- Order of applying optimizations within the phase
- Which temp file patterns qualify as "leaks" during the audit

</decisions>

<specifics>
## Specific Ideas

- Verification-first approach for dead code: prove it's dead, then remove — don't trust analysis alone
- Context managers are the preferred cleanup pattern (Pythonic, automatic)
- Performance report should be a living document tracking gains across phases

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 07-quick-wins-tier-1*
*Context gathered: 2026-02-14*
