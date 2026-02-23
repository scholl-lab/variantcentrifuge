# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** Planning next milestone

## Current Position

Phase: 29 of 29 (v0.15.0 complete)
Plan: N/A
Status: Milestone v0.15.0 shipped — ready for next milestone
Last activity: 2026-02-23 — v0.15.0 milestone complete

Progress: ████████████████████████ 100% (v0.15.0 shipped)

## Accumulated Context

### Decisions

(Cleared at milestone — full log in PROJECT.md Key Decisions table and milestones/v0.15.0-ROADMAP.md)

### Architecture Invariants

- R backend: parallel_safe=False; rpy2 calls only from main thread
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (memory constraint)
- P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-23 — v0.15.0 milestone completion
Stopped at: Milestone archived
Resume file: None
Next: `/gsd:new-milestone` for next milestone (requires `/clear` first for fresh context)
