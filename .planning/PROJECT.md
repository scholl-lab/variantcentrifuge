# VariantCentrifuge

## What This Is

A Python CLI tool for filtering, annotating, and analyzing variants from multi-sample VCF files in clinical genomics. It processes large cohort VCFs (up to 22 GB, 5,000+ samples) through a stage-based pipeline: gene BED creation, variant extraction, filtering, field extraction, analysis (inheritance, scoring, statistics, association testing), and output generation (TSV/Excel/HTML/IGV reports). Used by the Scholl Lab for nephrology and rare disease variant analysis.

## Core Value

Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats.

## Requirements

### Validated

- CLI entry point with comprehensive argument parsing and config mapping
- Stage-based pipeline architecture (40+ modular stages) with dependency graph and parallel execution
- Gene BED creation with caching (snpEff + bedtools)
- Variant extraction and filtering via bcftools and SnpSift
- Field extraction via bcftools query (replaced SnpSift extractFields — 19x faster) — v0.13.0
- Genotype replacement deferred to output time (eliminated 7-hour bottleneck) — v0.13.0
- Three-pass inheritance analysis: deduction, compound het, prioritization
- Inheritance analysis vectorized with NumPy (40-47% faster) — v0.13.0
- Variant scoring with configurable formula models
- Gene burden analysis with Fisher's exact test and multiple testing correction
- Individual HTML report: modern JS stack, semantic color coding, table redesign, column-level filtering, accessibility, print/PDF — v0.14.0
- Statistics generation (gene-level, impact summary, variant type summary)
- Output: TSV, Excel (xlsxwriter + openpyxl two-pass), HTML (interactive), IGV reports — v0.13.0
- Pseudonymization and archiving
- Checkpoint/resume system with atomic file operations
- Pipeline-wide ResourceManager for memory/CPU management — v0.13.0
- Parallel processing via ThreadPoolExecutor/ProcessPoolExecutor with gene sorting load balancing — v0.13.0
- Per-stage memory reporting at INFO level — v0.13.0
- DataFrame optimization: PyArrow engine, categorical dtypes, itertuples (82-84% memory reduction) — v0.13.0
- Performance benchmark suite: 60+ tests, regression detection, ratio assertions — v0.13.0
- Snakemake workflow for batch VCF processing
- CI: pytest + ruff + mypy on push/PR, Docker build with GHCR
- Modular association engine with Fisher, burden (logistic/linear), SKAT-O, COAST, ACAT-O omnibus tests — v0.15.0
- Pure Python SKAT and COAST backends (R deprecated with opt-in) — v0.15.0
- Covariate adjustment with sample alignment, PCA integration (PLINK/AKT file formats) — v0.15.0
- Functional variant weights: Beta(MAF), CADD/REVEL-based, uniform, combined — v0.15.0
- ACAT-V per-variant and ACAT-O per-gene omnibus combination with single FDR — v0.15.0
- Diagnostics: lambda_GC, QQ plot data/visualization, sample size warnings — v0.15.0
- Gene-level parallelization via ProcessPoolExecutor with BLAS thread control — v0.15.0
- JSON config mode for association analysis — v0.15.0
- Compiled Davies method (qfc.c via ctypes) with saddlepoint/Liu fallback chain — v0.15.0
- 2001 tests passing, 0 failures, cross-platform (Windows + Linux)

### Active

(No active requirements — next milestone not yet started)

### Out of Scope

- CI benchmark workflow (github-action-benchmark) — defer to after local benchmarks are stable
- Real-world test datasets (#60) — separate milestone
- Report generation validation (#61) — separate milestone
- Polars migration — unacceptable risk for clinical tool with deep pandas integration
- Free-threaded Python 3.13+ — numpy/pandas lack support
- Meta-analysis (RAREMETAL, Meta-SAIGE) — different problem: combining results across cohorts
- Mixed model / GRM (SAIGE-GENE approach) — biobank-scale; PCs + kinship exclusion sufficient for GCKD
- Adaptive test selection per gene — post-hoc best-test selection inflates type I error

## Context

- Current version: v0.15.0
- Tech stack: Python 3.10+, pandas, NumPy, openpyxl, xlsxwriter, bcftools, SnpSift, scipy, statsmodels
- Source: 38,948 LOC Python; Tests: 53,970 LOC Python
- Pipeline time reduced from 10+ hours to under 1 hour on large cohorts (>10x improvement) — v0.13.0
- Individual HTML report overhauled with modern JS stack, semantic colors, filtering, accessibility — v0.14.0
- Association framework validated on GCKD cohort (5,125 samples): PKD1 Fisher p=4.7e-39, SKAT p=3.47e-17 — v0.15.0
- GCKD cohort (5,125 samples) is the primary validation dataset

## Constraints

- **Correctness**: Clinical tool, zero tolerance for behavioral changes in existing functionality
- **Backwards compatibility**: Existing CLI interface and output formats must not change; new features are additive
- **Test coverage**: Statistical tests validated against R SKAT reference implementation
- **Dependencies**: Optional deps (rpy2, matplotlib) degrade gracefully; core tests use scipy/statsmodels only
- **Platform**: Must pass on both Windows and Linux (CI runs Ubuntu); compiled qfc.c needs cross-platform build
- **Numerical precision**: Python SKAT p-values within 10% relative difference of R; Davies ctypes within 1%

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Benchmarks before optimization | Can't prove improvement without baseline measurements | ✓ Good |
| Column sanitization at load time | GEN[0].GT → GEN_0__GT for itertuples compatibility, restored at output | ✓ Good |
| Conservative compound het | No pairing when phase ambiguous, avoids false positives | ✓ Good |
| bcftools query replaces SnpSift extractFields | 19x faster, C-based, per-sample GT columns | ✓ Good |
| Genotype replacement eliminated | Deferred to output time, saves 7 hours on large cohorts | ✓ Good |
| Pipeline-wide ResourceManager | Unified memory/CPU management, replaces inheritance-specific manager | ✓ Good |
| PyArrow + categorical dtypes | 82-84% memory reduction, 3x I/O speedup | ✓ Good |
| itertuples migration | 30.9x iteration speedup over iterrows | ✓ Good |
| xlsxwriter + openpyxl two-pass | Fast bulk write + rich formatting finalization | ✓ Good |
| Expandable rows over side panel | Simpler, DataTables-native pattern, avoids complexity | ✓ Good |
| Modernize JS stack (DataTables v2) | jQuery-optional; reduce bundle size, future-proof | ✓ Good |
| Dual SKAT backend (R + Python) | R is gold standard for validation; Python for portability | ✓ Good |
| Compiled Davies via ctypes | Exact p-values without R dependency; Liu fallback for safety | ✓ Good |
| Single FDR on ACAT-O across genes | Separate per-test FDR is statistically incorrect | ✓ Good |
| Python default backends (R deprecated) | R probe expensive/error-prone; Python always available | ✓ Good |
| GL quadrature for SKAT-O | 46x speedup; 128 nodes sufficient for smooth integrands | ✓ Good |
| ProcessPoolExecutor for gene parallelization | Near-linear speedup; BLAS thread control prevents oversubscription | ✓ Good |
| Clean reimplementation for FisherExactTest | Correct coupling direction for future gene_burden deprecation | ✓ Good |

---
*Last updated: 2026-02-23 after v0.15.0 milestone*
