# Phase 8: DataFrame Optimization - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Optimize pandas DataFrame loading, dtypes, and iteration patterns to achieve 50-70% memory reduction and 2-3x I/O speedup. Changes apply to stage-based pipeline only (classic pipeline is slated for deprecation — no changes). Must preserve byte-identical output for all pipeline modes.

</domain>

<decisions>
## Implementation Decisions

### Categorical dtype scope
- Auto-detect ALL low-cardinality columns (<50 unique values) and convert to categorical, not just the 5 listed in roadmap
- Convert at CSV load time (pass dtype dict to pd.read_csv), not post-load
- All-or-nothing rollback strategy: if any categorical column causes test failures, revert all categoricals until root cause is fixed

### Categorical pd.NA handling
- **Research needed:** Researcher should investigate best practices for handling pd.NA comparison breakage with categorical dtypes (pd.NA == 'value' returns pd.NA, not False). Options include fixing comparison sites, pre-populating categories, or late conversion. Pick the approach a senior software engineer would recommend for maintainability.

### Iteration replacement strategy
- Replace iterrows with itertuples on ALL hot paths, not just inheritance Pass 2-3 (includes gene_burden, scoring, statistics)
- Diff-based validation: run before/after on test data, byte-compare output files to catch any divergence
- Permanently rename columns with invalid Python identifiers (e.g., 'GEN[0].GT') at load time — all downstream code uses new clean names
- No temporary rename-restore patterns

### In-memory pass-through design
- Use shared memory reference for cross-process DataFrame access between stages
- **Research needed:** Investigate multiprocessing.shared_memory vs mmap-backed files for pandas DataFrames — pick best approach for performance and development ergonomics as a senior systems engineer would
- Always write TSV to disk as standard output artifact (users depend on it, useful for debugging)
- Auto-fallback to disk-based pass-through when DataFrame exceeds size threshold
- **Research needed:** Determine optimal memory threshold for auto-fallback targeting normal desktops with 8-16GB RAM
- Context managers for shared memory cleanup (guaranteed cleanup even on exceptions)
- Stage-based pipeline only — no changes to classic pipeline

### PyArrow migration scope
- Main variant DataFrame load only (not config reads, gene lists, BED files, etc.)
- **Research needed:** Investigate ArrowDtype string maturity in pandas 2.x — whether to convert back to object strings after loading or keep ArrowDtype throughout. Prioritize maintainability and performance.
- PyArrow becomes a required dependency (not optional with fallback)
- Pin to latest stable PyArrow version (>= 14.0)
- Trust existing benchmarks and unit tests for validation — no separate engine comparison tests

### Claude's Discretion
- Exact column rename mapping for invalid Python identifiers
- Implementation details of shared memory lifecycle management
- Which specific groupby/comparison sites need pd.NA-safe updates
- Ordering of optimizations within the phase plans

</decisions>

<specifics>
## Specific Ideas

- Target systems are normal desktops with 8-16GB RAM — optimizations should work well on these, not just beefy servers
- Classic pipeline (pipeline.py) gets zero changes — it's headed for deprecation (DEPR-01 in backlog)
- Context managers preferred for resource cleanup patterns

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 08-dataframe-optimization*
*Context gathered: 2026-02-14*
