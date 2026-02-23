# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.16.0 — Association Hardening & Multi-Cohort Features

## Current Position

Phase: Not started (defining requirements)
Plan: —
Status: Defining requirements
Last activity: 2026-02-23 — Milestone v0.16.0 started

Progress: ░░░░░░░░░░ 0%

## Accumulated Context

### Decisions

(New milestone — decisions logged as they're made)

### Architecture Invariants

- R backend: parallel_safe=False; rpy2 calls only from main thread
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (memory constraint)
- P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-23 — v0.16.0 milestone initialization
Stopped at: Requirements definition
Resume file: None
Next: Define requirements → create roadmap
