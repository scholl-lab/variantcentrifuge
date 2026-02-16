# VariantCentrifuge

## What This Is

A Python CLI tool for filtering, annotating, and analyzing variants from multi-sample VCF files in clinical genomics. It processes large cohort VCFs (up to 22 GB, 5,000+ samples) through a pipeline: gene BED creation, variant extraction, filtering, field extraction, analysis (inheritance, scoring, statistics, gene burden), and output generation (TSV/Excel/HTML/IGV reports). Used by the Scholl Lab for nephrology and rare disease variant analysis.

## Core Value

Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats.

## Requirements

### Validated

- CLI entry point with comprehensive argument parsing and config mapping
- Two pipeline architectures: classic (default) and stage-based (40+ modular stages)
- Gene BED creation with caching (snpEff + bedtools)
- Variant extraction and filtering via bcftools and SnpSift
- Field extraction via bcftools query (replaced SnpSift extractFields — 19x faster) — v0.13.0
- Genotype replacement deferred to output time (eliminated 7-hour bottleneck) — v0.13.0
- Three-pass inheritance analysis: deduction, compound het, prioritization
- Inheritance analysis vectorized with NumPy (40-47% faster) — v0.13.0
- Variant scoring with configurable formula models
- Gene burden analysis with Fisher's exact test and multiple testing correction
- Statistics generation (gene-level, impact summary, variant type summary)
- Output: TSV, Excel (xlsxwriter + openpyxl two-pass), HTML (interactive), IGV reports — v0.13.0
- Pseudonymization and archiving
- Checkpoint/resume system with atomic file operations
- Pipeline-wide ResourceManager for memory/CPU management (replaced InheritanceMemoryManager) — v0.13.0
- Parallel processing via ThreadPoolExecutor/ProcessPoolExecutor with gene sorting load balancing — v0.13.0
- Per-stage memory reporting at INFO level — v0.13.0
- DataFrame optimization: PyArrow engine, categorical dtypes, itertuples (82-84% memory reduction) — v0.13.0
- Performance benchmark suite: 60+ tests, regression detection, ratio assertions — v0.13.0
- Snakemake workflow for batch VCF processing
- CI: pytest + ruff + mypy on push/PR, Docker build with GHCR
- 1100+ tests passing, 0 failures, cross-platform (Windows + Linux)

### Active

(No active requirements — next milestone not yet defined)

### Out of Scope

- CI benchmark workflow (github-action-benchmark) — defer to after local benchmarks are stable
- Real-world test datasets (#60) — separate milestone
- Report generation validation (#61) — separate milestone
- Polars migration — unacceptable risk for clinical tool with deep pandas integration
- Free-threaded Python 3.13+ — numpy/pandas lack support
- Classic pipeline deprecation (DEPR-01) — tracked in backlog, needs deprecation warning first

## Context

- Current version: v0.13.0
- Shipped v0.13.0 with 28,378 LOC Python
- Tech stack: Python 3.10+, pandas, NumPy, openpyxl, xlsxwriter, bcftools, SnpSift
- Pipeline time reduced from 10+ hours to under 1 hour on large cohorts (>10x improvement)
- 4 open issues remain (#58, #60, #61, #62)
- 2 minor tech debt items from v0.13.0 (aspirational speedup target, dead code in GenotypeReplacementStage)

## Constraints

- **Correctness**: All optimizations must preserve identical output — clinical tool, zero tolerance for behavioral changes
- **Test coverage**: Every optimization requires benchmark measurements before and after
- **Backwards compatibility**: CLI interface and output formats must not change
- **Dependencies**: No new heavy dependencies except pytest-benchmark (dev)
- **Platform**: Must pass on both Windows and Linux (CI runs Ubuntu)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Benchmarks before optimization | Can't prove improvement without baseline measurements | Good |
| Local benchmarks only, no CI workflow | Get benchmarks working first, CI integration deferred | Good |
| Target v0.13.0 | Performance improvements warrant minor version bump | Good |
| Column sanitization at load time | GEN[0].GT → GEN_0__GT for itertuples compatibility, restored at output | Good |
| Conservative compound het | No pairing when phase ambiguous, avoids false positives | Good |
| bcftools query replaces SnpSift extractFields | 19x faster, C-based, per-sample GT columns | Good |
| Genotype replacement eliminated | Deferred to output time, saves 7 hours on large cohorts | Good |
| Pipeline-wide ResourceManager | Unified memory/CPU management, replaces inheritance-specific manager | Good |
| Vectorized comp het as sole implementation | Removed original comp_het.py, single implementation simplifies maintenance | Good |
| PyArrow + categorical dtypes | 82-84% memory reduction, 3x I/O speedup | Good |
| itertuples migration | 30.9x iteration speedup over iterrows | Good |
| xlsxwriter + openpyxl two-pass | Fast bulk write + rich formatting finalization | Good |

---
*Last updated: 2026-02-16 after v0.13.0 milestone*
