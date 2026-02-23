# Requirements: VariantCentrifuge v0.16.0

**Defined:** 2026-02-23
**Core Value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats

## v1 Requirements

### COAST Fix (#87)

- [ ] **COAST-01**: Genotype matrix is available to all regression-based tests (SKAT, burden, COAST) when per-sample GT columns exist in the DataFrame
- [ ] **COAST-02**: When `--association-tests coast` is specified, SIFT and PolyPhen annotation fields are auto-injected into field extraction
- [ ] **COAST-03**: COAST test produces valid p-values when 1 or 2 of 3 variant categories (BMV/DMV/PTV) are missing, with a warning logged (matching R AllelicSeries behavior)
- [ ] **COAST-04**: COAST classification handles SnpEff "&"-concatenated effect strings (e.g., `stop_gained&splice_region_variant`)
- [ ] **COAST-05**: User can select a COAST classification model via `--coast-classification <model-name>` (default: `sift_polyphen`)
- [ ] **COAST-06**: A `scoring/coast_classification/` config defines BMV/DMV/PTV thresholds as configuration, supporting SIFT+PolyPhen, CADD-based, and REVEL-based classification strategies
- [ ] **COAST-07**: Multi-transcript annotation conflicts use configurable resolution strategy (default: majority-vote) in the classification scoring model

### Gene-Level FDR Weighting (#86)

- [ ] **FDR-01**: User can provide a gene-to-weight TSV file via `--gene-prior-weights <file>` for weighted FDR correction
- [ ] **FDR-02**: Weighted BH (Genovese 2006) is applied when weights are provided: p-values adjusted by gene weight, weights renormalized to average=1.0
- [ ] **FDR-03**: Genes absent from the weight file receive weight=1.0 (neutral); warning emitted if >50% of tested genes are missing
- [ ] **FDR-04**: IHW (Independent Hypothesis Weighting) is available via `--gene-prior-method ihw` when rpy2 and Bioconductor IHW package are installed
- [ ] **FDR-05**: Effective number of tests (`sum(w)^2 / sum(w^2)`) is reported in diagnostics output
- [ ] **FDR-06**: When `--gene-prior-weights` is not provided, all genes receive equal weight (identical to current plain BH/Bonferroni)

### Case-Confidence Weights (#85)

- [ ] **CONF-01**: User can provide per-sample confidence weights from phenotype file via `--case-confidence file --case-confidence-column <col>`
- [ ] **CONF-02**: User can auto-compute confidence weights from HPO term overlap via `--case-confidence count`
- [ ] **CONF-03**: Logistic burden test uses `freq_weights` parameter from statsmodels for weighted regression
- [ ] **CONF-04**: Linear burden test uses WLS (weighted least squares) with sample weights
- [ ] **CONF-05**: Fisher's exact test ignores sample weights (integer counts only) with info log message
- [ ] **CONF-06**: Effective sample size (`sum(w)^2 / sum(w^2)`) is reported in diagnostics output
- [ ] **CONF-07**: `--case-confidence equal` (default) produces uniform weights, preserving exact backward compatibility

### Region Restriction (#84)

- [ ] **REGION-01**: User can restrict variant extraction to a BED file via `--restrict-regions <bed_file>`
- [ ] **REGION-02**: The gene BED is intersected with the restriction BED via `bedtools intersect` before `bcftools view -R`
- [ ] **REGION-03**: Chromosome naming mismatches between BED files are detected and reported as an error
- [ ] **REGION-04**: The intersection happens once before chunk splitting (all chunks respect the restriction)

### PCA Pipeline Integration (TD-01)

- [ ] **PCA-01**: `PCAComputationStage` is wired into `pipeline.py` and executes when `--pca-tool` is set
- [ ] **PCA-02**: PCA eigenvectors computed by the stage are available to `AssociationAnalysisStage` via `context.stage_results`
- [ ] **PCA-03**: `--pca-tool akt` runs AKT via subprocess and produces eigenvec output

### Performance Optimization

- [ ] **PERF-01**: Sparse genotype matrix representation (scipy.sparse.csr_array) is used when genotype density < threshold (opt-in via `--sparse-genotypes`)
- [ ] **PERF-02**: SKAT kernel computation works correctly with sparse input matrices
- [ ] **PERF-03**: SKAT-O eigendecomposition is profiled; if redundant eigh calls are confirmed, single precomputed eigendecomposition with rank-1 updates is implemented

### Dead Code Cleanup (#88)

- [ ] **CLEAN-01**: 8 dead stage classes removed from processing_stages.py, analysis_stages.py, output_stages.py
- [ ] **CLEAN-02**: Corresponding `register_stage()` calls removed from stage_registry.py
- [ ] **CLEAN-03**: Corresponding imports and `__all__` entries removed from stages/__init__.py
- [ ] **CLEAN-04**: `"refactored_pipeline"` default strings replaced with `"pipeline"` in pipeline.py and runner.py
- [ ] **CLEAN-05**: "Old pipeline" comments (~15 occurrences) reworded or removed
- [ ] **CLEAN-06**: "Refactored" docstrings updated in error_handling.py and analyze_variants.py

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
| COAST-01 | — | Pending |
| COAST-02 | — | Pending |
| COAST-03 | — | Pending |
| COAST-04 | — | Pending |
| COAST-05 | — | Pending |
| COAST-06 | — | Pending |
| COAST-07 | — | Pending |
| FDR-01 | — | Pending |
| FDR-02 | — | Pending |
| FDR-03 | — | Pending |
| FDR-04 | — | Pending |
| FDR-05 | — | Pending |
| FDR-06 | — | Pending |
| CONF-01 | — | Pending |
| CONF-02 | — | Pending |
| CONF-03 | — | Pending |
| CONF-04 | — | Pending |
| CONF-05 | — | Pending |
| CONF-06 | — | Pending |
| CONF-07 | — | Pending |
| REGION-01 | — | Pending |
| REGION-02 | — | Pending |
| REGION-03 | — | Pending |
| REGION-04 | — | Pending |
| PCA-01 | — | Pending |
| PCA-02 | — | Pending |
| PCA-03 | — | Pending |
| PERF-01 | — | Pending |
| PERF-02 | — | Pending |
| PERF-03 | — | Pending |
| CLEAN-01 | — | Pending |
| CLEAN-02 | — | Pending |
| CLEAN-03 | — | Pending |
| CLEAN-04 | — | Pending |
| CLEAN-05 | — | Pending |
| CLEAN-06 | — | Pending |
| TD-02 | — | Pending |
| TD-03 | — | Pending |
| TD-04 | — | Pending |
| TD-05 | — | Pending |
| TD-06 | — | Pending |

**Coverage:**
- v1 requirements: 41 total
- Mapped to phases: 0 (pending roadmap)
- Unmapped: 41

---
*Requirements defined: 2026-02-23*
*Last updated: 2026-02-23 after initial definition*
