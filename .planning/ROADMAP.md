# Roadmap: VariantCentrifuge

## Milestones

- [x] **v0.12.1 Baseline** - Phases 1-5 (shipped 2026-02-14)
- [x] **v0.13.0 Performance Optimization** - Phases 6-12 (shipped 2026-02-16)
- [x] **v0.14.0 Report UX Overhaul** - Phases 13-17 (shipped 2026-02-19)
- [x] **v0.15.0 Modular Rare Variant Association Framework** - Phases 18-29 (shipped 2026-02-23)
- [x] **v0.16.0 Association Hardening & Multi-Cohort Features** - Phases 30-37 (shipped 2026-02-26)
- [ ] **v0.17.0 Tech Debt Cleanup & Compound Het Parallelization** - Phases 38-39 (in progress)

## Phases

<details>
<summary>Phases 1-37 (v0.12.1 through v0.16.0) — SHIPPED</summary>

Phases 1-5: v0.12.1 Baseline (pre-GSD tracking)
Phases 6-12: v0.13.0 Performance Optimization
Phases 13-17: v0.14.0 Report UX Overhaul
Phases 18-29: v0.15.0 Modular Rare Variant Association Framework
Phases 30-37: v0.16.0 Association Hardening & Multi-Cohort Features

See MILESTONES.md for details.

</details>

---

## v0.17.0 — Tech Debt Cleanup & Compound Het Parallelization

**Goal:** Eliminate dead code, correct stale documentation, resolve minor tech debt, and parallelize the GIL-bound compound het pass for measurable speedup on multi-core systems.

**Scope:** 16 requirements across 4 categories (Performance, Dead Code, Documentation, Cleanup). The cleanup work is entirely low-risk and behaviorally inert. The compound het parallelization is the single complex item and is isolated to its own phase.

---

### Phase 38: Codebase Cleanup

**Goal:** The codebase contains no dead code, stale documentation, or unresolved minor tech debt — any developer reading the source sees accurate, current state.

**Dependencies:** None (independent of Phase 39)

**Requirements:** DEAD-01, DEAD-02, DEAD-03, DEAD-04, DEAD-05, DOCS-01, DOCS-02, DOCS-03, CLEAN-01, CLEAN-02, CLEAN-03, CLEAN-04, CLEAN-05, CLEAN-06, CLEAN-07

**Plans:** 3 plans

Plans:
- [x] 38-01-PLAN.md — Dead code removal (DEAD-01 through DEAD-05, CLEAN-04, CLEAN-05)
- [x] 38-02-PLAN.md — Documentation and comment fixes (DOCS-01, DOCS-02, DOCS-03, CLEAN-01, CLEAN-02, CLEAN-07)
- [x] 38-03-PLAN.md — Code cleanup (CLEAN-03, CLEAN-06)

**Success Criteria:**

1. `stage_info.py` is deleted and no import of it exists anywhere in the codebase; all tests still pass.
2. `--coast-backend r` is not present in CLI help output or argument choices; attempting to pass it raises an argument error.
3. The dead functions identified in `prioritizer.py` (5 functions) and `analyzer.py` (4 functions) are removed and no call sites reference them.
4. `docs/faq.md`, `docs/association_testing.md`, and `docs/changelog.md` accurately reflect the current CLI interface and column naming — no references to removed flags or outdated column names remain.
5. `make ci-check` passes with zero lint errors and zero test failures after all cleanup is applied.

---

### Phase 39: Compound Het Parallelization

**Goal:** Users running large cohorts on multi-core machines get measurably faster compound het analysis — the GIL-bound Pass 2 loop replaced by true parallelism.

**Dependencies:** Phase 38 (clean codebase as starting point, no functional dependency)

**Requirements:** PERF-01, PERF-02

**Plans:** 2 plans

Plans:
- [ ] 39-01-PLAN.md — Synthetic benchmark script and baseline capture
- [ ] 39-02-PLAN.md — Core optimization (pre-dispatch dedup, pedigree arrays, numpy-only workers, batch size env var)

**Success Criteria:**

1. `parallel_analyzer.py` eliminates GIL contention in compound het Pass 2 workers by pre-dispatching DataFrame operations and pre-computing pedigree arrays as NumPy integer arrays; workers receive only NumPy arrays in the hot path.
2. A benchmark or test demonstrates measurable wall-time speedup when compound het analysis runs on a multi-gene dataset with 2+ available CPU cores.
3. All existing compound het tests pass without modification — no behavioral change in pairing logic or output.
4. The implementation handles edge cases (single gene, gene with no het variants, 1 CPU) without error or regression.

---

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 30. Dead Code Cleanup | v0.16.0 | 1/1 | Complete | 2026-02-23 |
| 31. COAST Fix | v0.16.0 | 2/2 | Complete | 2026-02-23 |
| 32. Region Restriction and PCA Wiring | v0.16.0 | 2/2 | Complete | 2026-02-23 |
| 33. Gene-Level FDR Weighting | v0.16.0 | 1/1 | Complete | 2026-02-24 |
| 34. Tech Debt | v0.16.0 | 3/3 | Complete | 2026-02-24 |
| 35. Case-Confidence Weights | v0.16.0 | 0/2 | Deferred | - |
| 36. Performance — Sparse Genotype Matrices | v0.16.0 | 0/1 | Deferred | - |
| 37. Association Resource Management & Memory Streaming | v0.16.0 | 3/3 | Complete | 2026-02-25 |
| 38. Codebase Cleanup | v0.17.0 | 3/3 | Complete | 2026-02-26 |
| 39. Compound Het Parallelization | v0.17.0 | 0/2 | Pending | - |
