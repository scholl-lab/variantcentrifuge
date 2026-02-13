# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** Phase 6 - Benchmark Framework

## Current Position

Phase: 6 of 12 (Benchmark Framework)
Plan: None yet (ready to plan)
Status: Ready to plan
Last activity: 2026-02-14 — Roadmap created for milestone v0.13.0 Performance Optimization

Progress: [█████░░░░░░░░░░░░░░░] 25% (Phases 1-5 complete, pre-GSD baseline)

## Performance Metrics

**Velocity:**
- Total plans completed: 0 (milestone just started)
- Average duration: N/A
- Total execution time: 0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1-5. Baseline | N/A | N/A | N/A (pre-GSD) |

**Recent Trend:**
- Last 5 plans: None yet
- Trend: N/A (milestone initialization)

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Benchmarks before optimization: Can't prove improvement without baseline measurements
- Accept risk on Tier 3 (vectorization, Cython): User wants maximum performance gains, will revert if tests break
- Local benchmarks only, no CI workflow: Get benchmarks working first, CI integration deferred
- Target v0.13.0: Performance improvements warrant minor version bump

### Pending Todos

None yet.

### Blockers/Concerns

**Phase 6 (Benchmark Framework):**
- Must use synthetic or public datasets only, never reference private cohort names in committed code
- Ratio assertions needed for zero-flakiness comparisons within same run
- Memory profiling via tracemalloc required for budget enforcement

**Phase 8 (DataFrame Optimization):**
- PyArrow engine + pandas 3.0 string dtype changes require careful testing to avoid `pd.NA` comparison breakage
- Categorical dtypes depend on observed=True from Phase 7 to prevent 3500x groupby slowdown

**Phase 9 (Inheritance Analysis Optimization):**
- Full vectorization (INHER-03) is high-risk Tier 3 work, requires extensive validation to preserve clinical correctness
- Must preserve byte-identical output for all inheritance patterns

**Phase 11 (Pipeline & Cython Optimization):**
- Subprocess piping error handling needs platform-specific validation (Windows vs Linux)
- Cython extension requires hatch-cython integration with existing hatchling build

## Session Continuity

Last session: 2026-02-14 (roadmap creation)
Stopped at: Roadmap and STATE.md initialized for v0.13.0
Resume file: None
