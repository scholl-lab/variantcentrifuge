# Requirements: VariantCentrifuge v0.16.0

**Defined:** 2026-02-23
**Core Value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats

## v1 Requirements

### COAST Fix (#87)

- [x] **COAST-01**: Genotype matrix is available to all regression-based tests (SKAT, burden, COAST) when per-sample GT columns exist in the DataFrame
- [x] **COAST-02**: When `--association-tests coast` is specified, SIFT and PolyPhen annotation fields are auto-injected into field extraction
- [x] **COAST-03**: COAST test produces valid p-values when 1 or 2 of 3 variant categories (BMV/DMV/PTV) are missing, with a warning logged (matching R AllelicSeries behavior)
- [x] **COAST-04**: COAST classification handles SnpEff "&"-concatenated effect strings (e.g., `stop_gained&splice_region_variant`)
- [x] **COAST-05**: User can select a COAST classification model via `--coast-classification <model-name>` (default: `sift_polyphen`)
- [x] **COAST-06**: A `scoring/coast_classification/` config defines BMV/DMV/PTV thresholds as configuration, supporting SIFT+PolyPhen, CADD-based, and REVEL-based classification strategies
- [x] **COAST-07**: Multi-transcript annotation conflicts use configurable resolution strategy (default: majority-vote) in the classification scoring model

### Gene-Level FDR Weighting (#86)

- [x] **FDR-01**: User can provide a gene-to-weight TSV file via `--gene-prior-weights <file>` for weighted FDR correction
- [x] **FDR-02**: Weighted BH (Genovese 2006) is applied when weights are provided: p-values adjusted by gene weight, weights renormalized to average=1.0
- [x] **FDR-03**: Genes absent from the weight file receive weight=1.0 (neutral); warning emitted if >50% of tested genes are missing
- [ ] **FDR-04**: IHW (Independent Hypothesis Weighting) is available via `--gene-prior-method ihw` when rpy2 and Bioconductor IHW package are installed
- [x] **FDR-05**: Effective number of tests (`sum(w)^2 / sum(w^2)`) is reported in diagnostics output
- [x] **FDR-06**: When `--gene-prior-weights` is not provided, all genes receive equal weight (identical to current plain BH/Bonferroni)

### Case-Confidence Weights (#85)

- [ ] **CONF-01**: User can provide per-sample confidence weights from phenotype file via `--case-confidence file --case-confidence-column <col>`
- [ ] **CONF-02**: User can auto-compute confidence weights from HPO term overlap via `--case-confidence count`
- [ ] **CONF-03**: Logistic burden test uses `freq_weights` parameter from statsmodels for weighted regression
- [ ] **CONF-04**: Linear burden test uses WLS (weighted least squares) with sample weights
- [ ] **CONF-05**: Fisher's exact test ignores sample weights (integer counts only) with info log message
- [ ] **CONF-06**: Effective sample size (`sum(w)^2 / sum(w^2)`) is reported in diagnostics output
- [ ] **CONF-07**: `--case-confidence equal` (default) produces uniform weights, preserving exact backward compatibility

### Region Restriction (#84)

- [x] **REGION-01**: User can restrict variant extraction to a BED file via `--regions-bed <bed_file>`
- [x] **REGION-02**: The gene BED is intersected with the restriction BED via `bedtools intersect` before `bcftools view -R`
- [x] **REGION-03**: Chromosome naming mismatches between BED files are detected and reported as an error
- [x] **REGION-04**: The intersection happens once before chunk splitting (all chunks respect the restriction)

### PCA Pipeline Integration (TD-01)

- [x] **PCA-01**: `PCAComputationStage` is wired into `pipeline.py` and executes when `--pca` is set
- [x] **PCA-02**: PCA eigenvectors computed by the stage are available to `AssociationAnalysisStage` via `context.stage_results`
- [x] **PCA-03**: `--pca akt` runs AKT via subprocess and produces eigenvec output

### Performance Optimization

- [ ] **PERF-01**: Sparse genotype matrix representation (scipy.sparse.csr_array) is used when genotype density < threshold (opt-in via `--sparse-genotypes`)
- [ ] **PERF-02**: SKAT kernel computation works correctly with sparse input matrices
- [ ] **PERF-03**: SKAT-O eigendecomposition is profiled; if redundant eigh calls are confirmed, single precomputed eigendecomposition with rank-1 updates is implemented

### Dead Code Cleanup (#88)

- [x] **CLEAN-01**: 8 dead stage classes removed from processing_stages.py, analysis_stages.py, output_stages.py
- [x] **CLEAN-02**: Corresponding `register_stage()` calls removed from stage_registry.py
- [x] **CLEAN-03**: Corresponding imports and `__all__` entries removed from stages/__init__.py
- [x] **CLEAN-04**: `"refactored_pipeline"` default strings replaced with `"pipeline"` in pipeline.py and runner.py
- [x] **CLEAN-05**: "Old pipeline" comments (~15 occurrences) reworded or removed
- [x] **CLEAN-06**: "Refactored" docstrings updated in error_handling.py and analyze_variants.py

### Tech Debt

- [ ] **TD-02**: `create_stages_from_config()` maps `perform_association` and `perform_gene_burden` keys from config dict
- [ ] **TD-03**: COAST R golden values generated and hardcoded in test file for Python-vs-R comparison
- [ ] **TD-04**: SKAT-O output uses `skat_o_p_value` column name when `skat_method=SKATO`
- [ ] **TD-05**: Fisher lambda_GC criterion clarified in ROADMAP (excluded from Fisher, applied to score-based tests)
- [ ] **TD-06**: Internal `"refactored_pipeline"` checkpoint strings replaced (covered by CLEAN-04)

## v2 Requirements

### Deferred from v0.16.0

- **CONF-SKAT**: Weighted SKAT with sample weights in null model — no reference implementation exists; needs published methodology
- **R-REMOVE**: Remove R backends entirely (rpy2 bridge) — deferred to v0.17.0
- **KINSHIP**: Kinship matrix / mixed models (SAIGE-GENE approach) — biobank-scale
- **AKT-PLINK-AUTO**: Automated PCA tool version detection and output format validation

## Out of Scope

| Feature | Reason |
|---------|--------|
| Weighted SKAT (sample weights in null model) | No reference implementation in any published tool; type I error control unvalidated |
| IHW without rpy2 | No Python IHW implementation exists; requires R dependency |
| Adaptive permutation p-values | High complexity; flag low-MAC genes in output instead |
| Meta-analysis (RAREMETAL) | Different problem: combining results across cohorts |
| Mixed model / GRM | Biobank-scale; PCs + kinship exclusion sufficient for GCKD |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| COAST-01 | Phase 31 | Complete |
| COAST-02 | Phase 31 | Complete |
| COAST-03 | Phase 31 | Complete |
| COAST-04 | Phase 31 | Complete |
| COAST-05 | Phase 31 | Complete |
| COAST-06 | Phase 31 | Complete |
| COAST-07 | Phase 31 | Complete |
| FDR-01 | Phase 33 | Complete |
| FDR-02 | Phase 33 | Complete |
| FDR-03 | Phase 33 | Complete |
| FDR-04 | Phase 33 | Excluded (IHW deferred) |
| FDR-05 | Phase 33 | Complete |
| FDR-06 | Phase 33 | Complete |
| CONF-01 | Phase 35 | Pending |
| CONF-02 | Phase 35 | Pending |
| CONF-03 | Phase 35 | Pending |
| CONF-04 | Phase 35 | Pending |
| CONF-05 | Phase 35 | Pending |
| CONF-06 | Phase 35 | Pending |
| CONF-07 | Phase 35 | Pending |
| REGION-01 | Phase 32 | Complete |
| REGION-02 | Phase 32 | Complete |
| REGION-03 | Phase 32 | Complete |
| REGION-04 | Phase 32 | Complete |
| PCA-01 | Phase 32 | Complete |
| PCA-02 | Phase 32 | Complete |
| PCA-03 | Phase 32 | Complete |
| PERF-01 | Phase 36 | Pending |
| PERF-02 | Phase 36 | Pending |
| PERF-03 | Phase 36 | Pending |
| CLEAN-01 | Phase 30 | Complete |
| CLEAN-02 | Phase 30 | Complete |
| CLEAN-03 | Phase 30 | Complete |
| CLEAN-04 | Phase 30 | Complete |
| CLEAN-05 | Phase 30 | Complete |
| CLEAN-06 | Phase 30 | Complete |
| TD-02 | Phase 34 | Pending |
| TD-03 | Phase 34 | Pending |
| TD-04 | Phase 34 | Pending |
| TD-05 | Phase 34 | Pending |
| TD-06 | Phase 30 | Complete (same change as CLEAN-04) |

**Coverage:**
- v1 requirements: 41 total
- Mapped to phases: 41
- Unmapped: 0

---
*Requirements defined: 2026-02-23*
*Last updated: 2026-02-24 — Phase 33 FDR requirements marked Complete (FDR-04 IHW excluded)*
