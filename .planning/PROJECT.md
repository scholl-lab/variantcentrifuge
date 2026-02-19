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
- Individual HTML report: modern JS stack, semantic color coding, table redesign, column-level filtering, accessibility, print/PDF — v0.14.0
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

**Milestone v0.15.0 — Modular Rare Variant Association Framework**

- Modular association engine with pluggable statistical tests (Fisher, burden, SKAT, SKAT-O, ACAT-O)
- Dual SKAT backend: R SKAT via rpy2 (gold standard) + pure Python fallback with iterative validation
- Covariate adjustment system (age, sex, custom covariates + PCA components)
- PCA integration: pre-computed file loading + AKT/PLINK computation wrappers
- Variant weighting schemes: Beta(MAF), CADD/REVEL-based, uniform
- Logistic/linear regression burden tests with covariates (statsmodels)
- ACAT-V (variant-level) and ACAT-O (omnibus p-value combination)
- Compiled Davies method (qfc.c via ctypes) for mixture-of-χ² p-values
- Diagnostics: lambda_GC, QQ plot data, Manhattan plot data
- Allelic series test (LoF vs damaging vs benign ordered alternatives)
- JSON config mode for power users
- Backward compatible: existing --perform-gene-burden unchanged

### Out of Scope

- CI benchmark workflow (github-action-benchmark) — defer to after local benchmarks are stable
- Real-world test datasets (#60) — separate milestone
- Report generation validation (#61) — separate milestone
- Polars migration — unacceptable risk for clinical tool with deep pandas integration
- Free-threaded Python 3.13+ — numpy/pandas lack support
- Classic pipeline deprecation (DEPR-01) — tracked in backlog, needs deprecation warning first

## Context

- Current version: v0.14.0
- Tech stack: Python 3.10+, pandas, NumPy, openpyxl, xlsxwriter, bcftools, SnpSift, scipy, statsmodels
- Pipeline time reduced from 10+ hours to under 1 hour on large cohorts (>10x improvement) — v0.13.0
- Individual HTML report overhauled with modern JS stack, semantic colors, filtering, accessibility — v0.14.0
- Current gene burden is Fisher's exact only — adequate for strong Mendelian signals (PKD1 p=2e-80) but underpowered for novel associations
- Design document: `.planning/ASSOCIATION-FRAMEWORK-DESIGN.md` (detailed architecture, math specs, validation strategy)
- R SKAT package (Lee Lab) is the gold standard for rare variant association — we port from it
- Davies method (Algorithm AS 155) is the critical piece for mixture-of-χ² p-values
- GCKD cohort (5,125 samples) is the primary validation dataset

## Constraints

- **Correctness**: Clinical tool, zero tolerance for behavioral changes in existing functionality
- **Backwards compatibility**: Existing CLI interface and output formats must not change; new features are additive
- **Test coverage**: Statistical tests validated against R SKAT reference implementation
- **Dependencies**: New optional deps (rpy2, matplotlib) must degrade gracefully; core tests use scipy/statsmodels only
- **Platform**: Must pass on both Windows and Linux (CI runs Ubuntu); compiled qfc.c needs cross-platform build
- **Numerical precision**: Python SKAT p-values within 10% relative difference of R; Davies ctypes within 1%

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

| Individual report only for v0.14.0 | Cohort report is already better; focus on biggest gap first | ✓ Good |
| Expandable rows over side panel | Simpler, DataTables-native pattern, avoids complexity | ✓ Good |
| Modernize JS stack | DataTables v2 is jQuery-optional; reduce bundle size, future-proof | ✓ Good |
| Dual SKAT backend (R + Python) | R is gold standard for validation; Python for portability | — Pending |
| Compiled Davies via ctypes | Exact p-values without R dependency; Liu fallback for safety | — Pending |
| AKT + PLINK for PCA | Both supported — AKT for VCF-native, PLINK for pre-existing BED files | — Pending |

---
*Last updated: 2026-02-19 after v0.15.0 milestone start*
