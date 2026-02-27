# Project Milestones: VariantCentrifuge

## v0.17.0 Tech Debt Cleanup & Compound Het Parallelization (Shipped: 2026-02-27)

**Delivered:** Cleaned dead code (509-line module, 11 dead functions), fixed stale documentation, resolved minor tech debt. Compound het parallelization attempted and reverted after real-data regression — lesson: always benchmark with production-scale data.

**Phases completed:** 38-39 (5 plans total; Phase 39 optimization reverted)

**Key accomplishments:**

- Deleted `stage_info.py` (509 lines dead code) and 11 dead inheritance functions from prioritizer/analyzer
- Removed `--coast-backend r` CLI choice and fixed `--gzip-intermediates` flag confusion
- Corrected stale docs: faq.md removed flags, association_testing.md column name, changelog.md classic pipeline ref
- Added 6 missing `__all__` exports to `stages/__init__.py`, resolved stale TODOs
- Created compound het parallelization benchmark suite (retained for future optimization)
- Post-mortem: numpy-only workers 2x slower on real data (5125 samples) despite synthetic benchmark gains

**Stats:**

- 45 files created/modified
- 93,390 lines of Python (38,407 source + 54,983 tests)
- 2 phases, 5 plans (1 phase reverted)
- 2 days from start to ship (2026-02-26 → 2026-02-27)

**Git range:** `docs(38)` → `docs(39)`

**What's next:** TBD — candidates include case-confidence weights (#85), real-world test datasets (#60), report validation (#61)

---

## v0.16.0 Association Hardening & Multi-Cohort Features (Shipped: 2026-02-26)

**Delivered:** Hardened the v0.15.0 association framework for real-world use — fixed COAST p=None bugs, added BED region restriction, PCA pipeline wiring, weighted FDR correction, and streaming genotype matrices to prevent OOM on large gene panels.

**Phases completed:** 30-37 (12 plans total; Phases 35+36 deferred)

**Key accomplishments:**

- Fixed COAST producing p=None on partial-category genes with configurable classification scoring (3 models) and multi-transcript effect resolution
- Added --regions-bed BED-based region restriction with chromosome mismatch detection
- Wired PCAComputationStage with unified --pca flag (file or AKT autodetect)
- Implemented Genovese 2006 weighted Benjamini-Hochberg FDR correction with per-gene biological priors
- Standardized association column naming (_pvalue/_qvalue), added COAST golden value regression tests, fixed config mapping bugs
- Streaming genotype matrices (build-test-discard) with shared ResourceManager — O(1 gene) peak memory

**Stats:**

- 62 files created/modified
- 94,812 lines of Python (39,208 source + 55,604 tests)
- 6 phases, 12 plans (2 phases deferred)
- 3 days from start to ship (2026-02-23 → 2026-02-26)

**Git range:** `docs(30)` → `docs(37)`

**What's next:** TBD — candidates include case-confidence weights (#85), sparse matrices, real-world test datasets (#60), report validation (#61)

---

## v0.15.0 Modular Rare Variant Association Framework (Shipped: 2026-02-23)

**Delivered:** Multi-test rare variant association engine with SKAT-O, COAST, burden tests, ACAT-O omnibus, covariate/PCA adjustment, and pure Python backends — backward compatible with existing gene burden pipeline.

**Phases completed:** 18-29 (35 plans total)

**Key accomplishments:**

- Built modular AssociationTest framework with Fisher, logistic/linear burden, SKAT-O, COAST, and ACAT-O omnibus tests running simultaneously with unified output schema
- Implemented pure Python SKAT and COAST backends matching R reference output, eliminating R dependency as default; achieved 46x SKAT-O speedup via Gauss-Legendre quadrature
- Added covariate adjustment with sample alignment validation, PCA integration (PLINK/AKT file formats), and functional variant weights (CADD/REVEL/Beta-MAF)
- Created ACAT-O/ACAT-V omnibus combination with Cauchy method and single FDR correction strategy across genes
- Built diagnostics module with lambda-GC, QQ plot data/visualization, sample size warnings, and per-gene carrier count warnings
- Added ProcessPoolExecutor gene-level parallelization, JSON config mode, and 20 CLI args for reproducible HPC workflows

**Stats:**

- 210 files created/modified
- 92,918 lines of Python (38,948 source + 53,970 tests)
- 12 phases, 35 plans, 214 commits
- 5 days from start to ship (2026-02-19 → 2026-02-23)

**Git range:** `feat(18-01)` → `docs(29)`

**What's next:** TBD — candidates include real-world test datasets, report generation validation, CI benchmark integration

---

## v0.14.0 Report UX Overhaul (Shipped: 2026-02-19)

**Delivered:** Complete HTML report overhaul with modern JS stack, semantic color coding, table redesign, column-level filtering, accessibility, and print/PDF support.

**Phases completed:** 13-17 (15 plans total)

**Key accomplishments:**

- Modernized JS stack with DataTables v2 (jQuery-optional)
- Implemented semantic color coding for variant significance
- Redesigned table layout with expandable rows
- Added column-level filtering and visualization
- Added accessibility features and print/PDF optimization

**Stats:**

- 5 phases, 15 plans
- 4 days (2026-02-16 → 2026-02-19)

**Git range:** Phases 13-17

---

## v0.13.0 Performance Optimization (Shipped: 2026-02-16)

**Delivered:** Reduced large cohort pipeline time from 10+ hours to under 1 hour through systematic optimization of DataFrame operations, inheritance analysis vectorization, I/O elimination, and pipeline-wide resource management.

**Phases completed:** 6-12 (26 plans total)

**Key accomplishments:**

- Built 60-test benchmark framework with regression detection and cross-phase comparison
- Achieved 48-98% gene burden speedup through dead code removal and memory management
- Delivered 82-84% memory reduction and 30.9x iteration speedup via categorical dtypes and itertuples
- Achieved 40-47% inheritance analysis speedup through three-pass NumPy vectorization
- Eliminated genotype replacement stage (7 hrs) and replaced SnpSift extractFields with bcftools query (19x faster)
- Created pipeline-wide ResourceManager with auto-tuned parallelism and per-stage memory reporting
- Overall: >10x pipeline speedup (exceeding 3-4x target)

**Stats:**

- 183 files created/modified
- 28,378 lines of Python (+30,654 net)
- 7 phases, 26 plans, 129 commits
- 3 days from start to ship (2026-02-14 → 2026-02-16)

**Git range:** `docs(06)` → `docs(12)`

**What's next:** Next milestone TBD — candidates include CI benchmark integration, classic pipeline deprecation, real-world test datasets

---

## v0.12.1 Baseline (Shipped: 2026-02-14)

**Delivered:** Full-featured variant analysis pipeline with two architectures, 40+ stages, inheritance analysis, gene burden, scoring, multi-format output.

**Phases completed:** 1-5 (pre-GSD tracking)

**Key accomplishments:**

- Two pipeline architectures: classic (default) and stage-based (40+ modular stages)
- Three-pass inheritance analysis: deduction, compound het, prioritization
- Gene burden analysis with Fisher's exact test and multiple testing correction
- Multi-format output: TSV, Excel (with hyperlinks/formatting), HTML, IGV reports
- All 30 historical issues resolved, 1035 tests passing
- CI/CD with Docker, docs, and multi-platform testing

**Stats:**

- 1035 tests passing
- Cross-platform (Windows + Linux)
- Full CI/CD pipeline

---
