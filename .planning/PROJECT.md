# VariantCentrifuge

## What This Is

A Python CLI tool for filtering, annotating, and analyzing variants from multi-sample VCF files in clinical genomics. It processes large cohort VCFs (up to 22 GB, 5,000+ samples) through a 7-phase pipeline: gene BED creation, variant extraction, filtering, field extraction, genotype replacement, analysis (inheritance, scoring, statistics, gene burden), and output generation (TSV/Excel/HTML/IGV reports). Used by the Scholl Lab for nephrology and rare disease variant analysis.

## Core Value

Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats.

## Requirements

### Validated

- CLI entry point with comprehensive argument parsing and config mapping
- Two pipeline architectures: classic (default) and stage-based (40+ modular stages)
- Gene BED creation with caching (snpEff + bedtools)
- Variant extraction and filtering via bcftools and SnpSift
- Field extraction and genotype replacement (vectorized + sequential)
- Three-pass inheritance analysis: deduction, compound het, prioritization
- Variant scoring with configurable formula models
- Gene burden analysis with Fisher's exact test and multiple testing correction
- Statistics generation (gene-level, impact summary, variant type summary)
- Output: TSV, Excel (with hyperlinks/formatting), HTML (interactive), IGV reports
- Pseudonymization and archiving
- Checkpoint/resume system with atomic file operations
- Memory management for HPC environments (SLURM, PBS, cgroups)
- Parallel processing via ThreadPoolExecutor/ProcessPoolExecutor
- Snakemake workflow for batch VCF processing
- CI: pytest + ruff + mypy on push/PR, Docker build with GHCR
- 1035 tests passing, 0 failures, cross-platform (Windows + Linux)

### Active

- [ ] Performance benchmark framework (pytest-benchmark, local)
- [ ] Performance optimization: Tier 1 quick wins (observed=True, dead code removal, pyarrow engine, temp file leak)
- [ ] Performance optimization: Tier 2 medium effort (categoricals, itertuples, vectorize comp_het, DataFrame pass-through)
- [ ] Performance optimization: Tier 3 strategic (full inheritance vectorization, pipe fusion, xlsxwriter, Cython/Rust genotype kernel)

### Out of Scope

- CI benchmark workflow (github-action-benchmark) — defer to after local benchmarks are stable
- Real-world test datasets (#60) — separate milestone
- Report generation validation (#61) — separate milestone
- Polars migration — unacceptable risk for clinical tool with deep pandas integration
- Free-threaded Python 3.13+ — numpy/pandas lack support
- SIMD/cache-aware layouts — diminishing returns, high complexity
- SharedMemory for parallel processing — defer until profiling shows serialization as bottleneck
- GPU acceleration — insufficient numeric computation to justify

## Context

- Current version: v0.12.1
- All 30 historical issues resolved; 4 open issues remain (#58, #60, #61, #62)
- Performance analysis report (2026-02-13) identifies critical bottlenecks:
  - `df.apply(axis=1)` in inheritance deduction: 40-60% of analysis time
  - Nested `iterrows()` in gene burden GT parsing: 20-30% of analysis time (dead code)
  - No `observed=True` on any groupby: risk of 3500x slowdown with categoricals
  - `dtype=str` everywhere prevents numeric ops and PyArrow optimizations
  - Write-then-read temp files in analysis and output stages
- Expected pipeline speedup: ~3-4x for large cohorts (22 GB multi-sample VCF, 5000+ samples)
- Existing test infrastructure: 1035 tests, pytest markers (unit/integration/slow), ~3 min fast suite

## Constraints

- **Correctness**: All optimizations must preserve identical output — clinical tool, zero tolerance for behavioral changes
- **Test coverage**: Every optimization requires benchmark measurements before and after
- **Backwards compatibility**: CLI interface and output formats must not change
- **Dependencies**: No new heavy dependencies except pytest-benchmark (dev) and optionally Cython
- **Platform**: Must pass on both Windows and Linux (CI runs Ubuntu)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Benchmarks before optimization | Can't prove improvement without baseline measurements | — Pending |
| Accept risk on Tier 3 (vectorization, Cython) | User wants maximum performance gains, will revert if tests break | — Pending |
| Local benchmarks only, no CI workflow | Get benchmarks working first, CI integration deferred | — Pending |
| Target v0.13.0 | Performance improvements warrant minor version bump | — Pending |

---
*Last updated: 2026-02-14 after milestone v0.13.0 initialization*
