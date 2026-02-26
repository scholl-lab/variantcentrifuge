# Requirements: VariantCentrifuge v0.17.0

**Defined:** 2026-02-26
**Core Value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats

## v0.17.0 Requirements

### Performance

- [x] **PERF-01**: Compound het Pass 2 eliminates GIL contention via pre-dispatch dedup, pedigree arrays, and numpy-only workers
- [x] **PERF-02**: Compound het parallelization achieves measurable speedup on multi-core systems

### Dead Code

- [x] **DEAD-01**: `stage_info.py` module deleted (509 lines, never imported)
- [x] **DEAD-02**: `--coast-backend r` choice removed from CLI (no implementation exists)
- [x] **DEAD-03**: Dead inheritance functions removed from `prioritizer.py` (5 functions)
- [x] **DEAD-04**: Dead inheritance functions removed from `analyzer.py` (4 functions)
- [x] **DEAD-05**: Redundant `--gzip-intermediates` flag fixed or removed

### Documentation

- [x] **DOCS-01**: `faq.md` references to 4 removed CLI flags corrected
- [x] **DOCS-02**: `association_testing.md` column name corrected (skat_pvalue → skat_o_pvalue)
- [x] **DOCS-03**: `changelog.md` classic pipeline reference updated

### Cleanup

- [x] **CLEAN-01**: Stale "refactored pipeline" docstrings updated in 3 integration tests
- [x] **CLEAN-02**: Stale "original pipeline" comments updated in 3 source files
- [x] **CLEAN-03**: Missing `__all__` exports added to `stages/__init__.py`
- [x] **CLEAN-04**: Dead integration test methods removed (TestBCFToolsPrefilter, test_parallel_variant_extraction)
- [x] **CLEAN-05**: TODO(12-02) legacy chunking check resolved
- [x] **CLEAN-06**: TODO intelligent batching resolved or documented as deferred
- [x] **CLEAN-07**: TD-05 Fisher lambda_GC doc clarification addressed

## Future Requirements

### Deferred from v0.16.0

- **CONF-01**: Case-confidence weighted sample status (#85) — needs HPO infrastructure
- **SPARSE-01**: Sparse genotype matrices — streaming solves OOM; centering destroys sparsity

### Future Milestones

- Remove R backends for SKAT/COAST entirely
- Single eigendecomposition optimization for SKAT-O
- CI benchmark workflow (github-action-benchmark)
- Real-world test datasets (#60)
- Report generation validation (#61)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Polars migration | Unacceptable risk for clinical tool with deep pandas integration |
| Free-threaded Python 3.13+ | numpy/pandas lack support |
| IHW (FDR-04) | No Python implementation exists |
| Kinship matrix / mixed models | Requires SAIGE-style sparse GRM; overkill for GCKD |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| PERF-01 | Phase 39 | Complete |
| PERF-02 | Phase 39 | Complete |
| DEAD-01 | Phase 38 | Complete |
| DEAD-02 | Phase 38 | Complete |
| DEAD-03 | Phase 38 | Complete |
| DEAD-04 | Phase 38 | Complete |
| DEAD-05 | Phase 38 | Complete |
| DOCS-01 | Phase 38 | Complete |
| DOCS-02 | Phase 38 | Complete |
| DOCS-03 | Phase 38 | Complete |
| CLEAN-01 | Phase 38 | Complete |
| CLEAN-02 | Phase 38 | Complete |
| CLEAN-03 | Phase 38 | Complete |
| CLEAN-04 | Phase 38 | Complete |
| CLEAN-05 | Phase 38 | Complete |
| CLEAN-06 | Phase 38 | Complete |
| CLEAN-07 | Phase 38 | Complete |

**Coverage:**
- v0.17.0 requirements: 16 total
- Mapped to phases: 16
- Unmapped: 0

---
*Requirements defined: 2026-02-26*
*Last updated: 2026-02-26 after Phase 39 completion*
