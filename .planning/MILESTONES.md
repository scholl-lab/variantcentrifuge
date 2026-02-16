# Project Milestones: VariantCentrifuge

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
