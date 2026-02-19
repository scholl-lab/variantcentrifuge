# Requirements: v0.15.0 Modular Rare Variant Association Framework

**Defined:** 2026-02-19
**Core Value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats

## v1 Requirements

### Core Framework

- [ ] **CORE-01**: AssociationEngine orchestrates test selection, execution per gene, and result collection
- [ ] **CORE-02**: Abstract AssociationTest base class with standardized TestResult output (p-value, effect size, CI, variant counts)
- [ ] **CORE-03**: Fisher's exact test refactored into association/tests/fisher.py with bit-identical output to current --perform-gene-burden
- [ ] **CORE-04**: AssociationAnalysisStage registered in pipeline with dependency on dataframe_loading and sample_config_loading
- [ ] **CORE-05**: GeneBurdenAnalysisStage shim that skips when --perform-association is active (backward compatible)
- [ ] **CORE-06**: CLI args: --perform-association, --association-tests, --skat-backend (auto/r/python)
- [ ] **CORE-07**: PipelineContext gains association_results field; ExcelReportStage gains Association sheet
- [ ] **CORE-08**: Multiple testing correction refactored into association/correction.py (FDR + Bonferroni), re-exported in gene_burden.py

### Covariate System

- [ ] **COV-01**: Covariate file loading with explicit sample ID alignment to VCF sample order (assert no NaN after reindex)
- [ ] **COV-02**: Automatic one-hot encoding for categorical covariate columns
- [ ] **COV-03**: Multicollinearity check (condition number warning)
- [ ] **COV-04**: CLI args: --covariate-file, --covariates (column selection), --trait-type (binary/quantitative)

### Burden Tests

- [ ] **BURDEN-01**: Logistic regression burden test via statsmodels.Logit with Wald test, OR + 95% CI
- [ ] **BURDEN-02**: Linear regression burden test via statsmodels.OLS with beta + SE for quantitative traits
- [ ] **BURDEN-03**: Genotype matrix builder with correct handling of multi-allelic (1/2), missing (./.), phased (0|1) genotypes; mean imputation for missing

### Variant Weights

- [ ] **WEIGHT-01**: Beta(MAF; 1, 25) weights as default (scipy.stats.beta.pdf)
- [ ] **WEIGHT-02**: Uniform weights option (backward compatible with current behavior)
- [ ] **WEIGHT-03**: Functional weights from annotation columns (CADD-normalized, REVEL-based)
- [ ] **WEIGHT-04**: Combined weights: Beta(MAF) x functional score
- [ ] **WEIGHT-05**: CLI args: --variant-weights (beta/uniform/cadd/revel/combined), --variant-weight-params

### SKAT/SKAT-O

- [ ] **SKAT-01**: R SKAT backend via rpy2 with automatic R/SKAT package detection and graceful fallback
- [ ] **SKAT-02**: SKATBinary used by default for binary traits (never continuous-trait SKAT on binary phenotypes)
- [ ] **SKAT-03**: Small-sample moment adjustment (SKAT_Null_Model_MomentAdjust) as R backend default
- [ ] **SKAT-04**: SKAT-O with rho grid search and method="optimal.adj" correction
- [ ] **SKAT-05**: Pure Python SKAT backend validated against R within 10% relative difference on log10(p)
- [ ] **SKAT-06**: Davies method via ctypes-compiled qfc.c with corrected defaults (acc=1e-9, lim=10^6)
- [ ] **SKAT-07**: Fallback chain: Davies -> saddlepoint -> Liu moment-matching; p_method recorded in output
- [ ] **SKAT-08**: R backend declares parallel_safe=False; rpy2 calls only from main thread
- [ ] **SKAT-09**: R memory management: explicit del + rpy2 gc() every 100 genes
- [ ] **SKAT-10**: Eigenvalue stability: scipy.linalg.eigh, threshold max(eigenvalues, 0), skip if matrix_rank < 2

### Omnibus Tests

- [ ] **OMNI-01**: ACAT-V per-variant Cauchy combination of marginal score test p-values
- [ ] **OMNI-02**: ACAT-O omnibus combining burden + SKAT + ACAT-V p-values per gene
- [ ] **OMNI-03**: Single FDR correction applied to ACAT-O p-values across genes (not separate per test)

### Diagnostics and Output

- [ ] **DIAG-01**: Per-gene TSV output with standard columns (fisher_p, burden_p, skat_p, skat_o_p, acat_o_p, effect sizes, CIs, variant counts)
- [ ] **DIAG-02**: Lambda_GC genomic inflation factor computed per test
- [ ] **DIAG-03**: QQ plot data TSV (observed vs expected -log10(p))
- [ ] **DIAG-04**: Optional matplotlib QQ plot PNG/SVG (lazy import, graceful if matplotlib absent)
- [ ] **DIAG-05**: CLI arg: --diagnostics-output (path for diagnostics directory)
- [ ] **DIAG-06**: Warn when n_cases < 200 or case:control ratio > 1:20; flag genes with case_carriers < 10

### PCA Integration

- [ ] **PCA-01**: PCA file loading supporting PLINK .eigenvec, AKT output, and generic TSV formats
- [ ] **PCA-02**: PCA components merged as covariates (default 10 PCs; warn if >20)
- [ ] **PCA-03**: AKT wrapper as PCAComputationStage in pipeline (optional, requires akt in PATH)
- [ ] **PCA-04**: CLI args: --pca-file, --pca-tool akt, --pca-components

### Allelic Series

- [ ] **SERIES-01**: COAST allelic series test with BMV/DMV/PTV variant classification from PolyPhen/SIFT annotations
- [ ] **SERIES-02**: Configurable variant category weights (default w=1,2,3 for BMV/DMV/PTV)

### Configuration

- [ ] **CONFIG-01**: JSON config mode for association analysis (--association-config file.json)
- [ ] **CONFIG-02**: Config supports all CLI association options as JSON fields

## Future Requirements

- Permutation-based p-values for low-MAC situations (adaptive permutation)
- Kinship matrix support / mixed-model SKAT (GENESIS package)
- Conditional analysis (testing after conditioning on known signals)
- Manhattan plot data generation (gene-level)
- PLINK 2.0 wrapper for PCA computation
- Score statistics export for meta-analysis (RAREMETAL-compatible format)
- Dark mode for diagnostic plots

## Out of Scope

| Feature | Reason |
|---------|--------|
| Meta-analysis (RAREMETAL, Meta-SAIGE) | Different problem: combining results across cohorts; out of scope for single-cohort tool |
| Mixed model / GRM (SAIGE-GENE approach) | Biobank-scale; PCs + kinship exclusion sufficient for GCKD (<10K samples) |
| DeepRVAT neural network weights | Requires massive training data; no nephrology model exists |
| Phased haplotype tests | Existing comp_het.py already handles compound heterozygous detection |
| REGENIE-style whole-genome split | Designed for >100K samples; not needed for lab cohorts |
| Adaptive test selection per gene | Post-hoc best-test selection inflates type I error |
| PLINK 2.0 wrapper | AKT is sufficient for VCF-native PCA; PLINK file loading (pre-computed) supported |

## Traceability

(Populated during roadmap creation)

| Requirement | Phase | Status |
|-------------|-------|--------|

**Coverage:**
- v1 requirements: 44 total
- Mapped to phases: 0
- Unmapped: 44

---
*Requirements defined: 2026-02-19*
*Last updated: 2026-02-19 after initial definition*
