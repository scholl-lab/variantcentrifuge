# Requirements: VariantCentrifuge v0.17.0

**Defined:** 2026-02-26
**Core Value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats

## v0.17.0 Requirements

### Performance

- [ ] **PERF-01**: Compound het Pass 2 uses ProcessPoolExecutor or Numba instead of GIL-bound ThreadPoolExecutor
- [ ] **PERF-02**: Compound het parallelization achieves measurable speedup on multi-core systems

### Dead Code

- [ ] **DEAD-01**: `stage_info.py` module deleted (509 lines, never imported)
- [ ] **DEAD-02**: `--coast-backend r` choice removed from CLI (no implementation exists)
- [ ] **DEAD-03**: Dead inheritance functions removed from `prioritizer.py` (5 functions)
- [ ] **DEAD-04**: Dead inheritance functions removed from `analyzer.py` (4 functions)
- [ ] **DEAD-05**: Redundant `--gzip-intermediates` flag fixed or removed

### Documentation

- [ ] **DOCS-01**: `faq.md` references to 4 removed CLI flags corrected
- [ ] **DOCS-02**: `association_testing.md` column name corrected (skat_pvalue → skat_o_pvalue)
- [ ] **DOCS-03**: `changelog.md` classic pipeline reference updated

### Cleanup

- [ ] **CLEAN-01**: Stale "refactored pipeline" docstrings updated in 3 integration tests
- [ ] **CLEAN-02**: Stale "original pipeline" comments updated in 3 source files
- [ ] **CLEAN-03**: Missing `__all__` exports added to `stages/__init__.py`
- [ ] **CLEAN-04**: Dead integration test methods removed (TestBCFToolsPrefilter, test_parallel_variant_extraction)
- [ ] **CLEAN-05**: TODO(12-02) legacy chunking check resolved
- [ ] **CLEAN-06**: TODO intelligent batching resolved or documented as deferred
- [ ] **CLEAN-07**: TD-05 Fisher lambda_GC doc clarification addressed

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
| PERF-01 | TBD | Pending |
| PERF-02 | TBD | Pending |
| DEAD-01 | TBD | Pending |
| DEAD-02 | TBD | Pending |
| DEAD-03 | TBD | Pending |
| DEAD-04 | TBD | Pending |
| DEAD-05 | TBD | Pending |
| DOCS-01 | TBD | Pending |
| DOCS-02 | TBD | Pending |
| DOCS-03 | TBD | Pending |
| CLEAN-01 | TBD | Pending |
| CLEAN-02 | TBD | Pending |
| CLEAN-03 | TBD | Pending |
| CLEAN-04 | TBD | Pending |
| CLEAN-05 | TBD | Pending |
| CLEAN-06 | TBD | Pending |
| CLEAN-07 | TBD | Pending |

**Coverage:**
- v0.17.0 requirements: 16 total
- Mapped to phases: 0
- Unmapped: 16 ⚠️

---
*Requirements defined: 2026-02-26*
*Last updated: 2026-02-26 after initial definition*
