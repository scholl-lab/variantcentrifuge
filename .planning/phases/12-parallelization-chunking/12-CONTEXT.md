# Phase 12: Parallelization & Chunking - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Dynamic work distribution and memory-efficient processing at scale. Improve chunking and parallelism for inheritance analysis on large cohorts (5000+ samples). This phase must first investigate what optimizations are ACTUALLY needed given the dramatic pipeline changes in Phases 7-11 (82-84% memory reduction, 7-hour genotype replacement eliminated, 2.7-hour SnpSift replaced with bcftools).

</domain>

<decisions>
## Implementation Decisions

### Realistic scope — investigate first, implement only what's needed
- Phase MAY shrink significantly if investigation shows bottlenecks are gone
- Drop work stealing (adaptive redistribution) — smarter chunking is sufficient
- Drop memory pools (pre-allocated reusable buffers) — investigate whether needed given Phase 8 memory gains
- Drop async I/O + memory mapping — investigate whether I/O bottlenecks remain after Phase 11 eliminated SnpSift and genotype replacement
- Research-driven: act as senior systems engineer, use web search for best practices, thoroughly investigate codebase state
- Principles: KISS, DRY, SOLID, modularity, no antipatterns, clean code, no dead code

### Chunking strategy
- Auto-detect environment: detect available resources (CPU, RAM, SLURM/PBS/cgroups) and choose strategy accordingly
- Remove --chunks CLI flag: chunking becomes fully automatic, no manual override
- Log auto-detected chunk size at INFO level so users can see the system's decision
- Technical details (static vs dynamic, calculation method) to be determined by research

### Parallelism scope
- Focus on improving existing inheritance analysis parallelism only — do NOT expand to gene burden or scoring stages
- Leave runner.py stage-level parallelism (ThreadPoolExecutor) unchanged
- Auto-detect optimal worker count based on CPU cores and available memory — no --workers CLI flag
- Auto-disable parallel processing for small datasets where overhead isn't worth it (threshold determined by benchmarking)

### Memory management approach
- Replace inheritance_memory_manager.py with new pipeline-wide memory management (no dead code left behind)
- New module usable by any pipeline stage that needs memory awareness, not just inheritance
- Low memory behavior: warn and continue (log warning, let OS handle pressure)
- Report memory usage statistics (total peak, per-stage breakdown) at INFO level after pipeline completes

### Claude's Discretion
- Whether memory pools, async I/O, or mmap are actually needed (research-driven decision)
- Chunking calculation algorithm (static at start vs runtime adaptive)
- Small-dataset parallelism threshold
- Internal architecture of new memory management module
- Whether to consolidate or restructure existing parallel_analyzer.py

</decisions>

<specifics>
## Specific Ideas

- User emphasized: "reevaluate based on Phase 11 and recent changes" — the roadmap success criteria were written before Phases 7-11 transformed the pipeline
- User emphasized: "deeply investigate phases and code as a senior developer" — research must be thorough, not superficial
- User emphasized: "no dead code" — if inheritance_memory_manager.py is replaced, remove it completely
- User emphasized: "no regressions, no antipatterns" — all changes must be tested

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 12-parallelization-chunking*
*Context gathered: 2026-02-16*
