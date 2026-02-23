# Project Milestones: VariantCentrifuge

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
