# Feature Landscape — Modular Rare Variant Association Framework

**Domain:** Rare variant association testing for genomic sequencing studies
**Researched:** 2026-02-19 (v0.15.0), updated 2026-02-23 (v0.16.0)
**Target:** Adding association testing to variantcentrifuge (nephrology lab, ~5,000 samples, multi-sample VCF)
**Existing baseline:** Fisher's exact test with FDR/Bonferroni correction in `gene_burden.py`

---

## Context: What This Milestone Adds

The existing `gene_burden.py` implements a CMC/CAST collapsing test (Li & Leal 2008) with
Fisher's exact test. That works for strong Mendelian signals (e.g., PKD1 p=2e-80) but is
underpowered for novel associations, cannot adjust for covariates, and cannot handle
quantitative traits. This feature landscape covers what a credible rare variant association
framework must include versus what is genuinely differentiating versus what to avoid.

---

## Table Stakes

Features that users and reviewers expect from any credible rare variant association framework.
Missing these = framework is not publishable or scientifically credible.

| Feature | Why Expected | Complexity | Dependencies |
|---------|--------------|------------|--------------|
| **Fisher's exact test (refactored)** | Already exists; backward compatibility is non-negotiable | Low | Existing `gene_burden.py` — refactor to `association/tests/fisher.py` |
| **Logistic regression burden test** | Standard for binary traits with covariates; OR + 95% CI expected by reviewers | Medium | `statsmodels.api.Logit` (already a dependency); covariate system |
| **SKAT (variance-component score test)** | Gold standard for detecting heterogeneous effects; Wu et al. 2011 | High | numpy eigenvalues, Davies method (mixture of chi-squared), null model fitting |
| **SKAT-O (optimal rho combination)** | Combines burden + SKAT; robust across genetic architectures; Lee et al. 2012 | High | SKAT implementation + rho grid search + optimal.adj correction |
| **Beta(MAF; 1, 25) variant weights** | The dominant weighting scheme in all major tools (R SKAT, SAIGE-GENE, REGENIE); strongly upweights rare variants | Low | MAF computation from genotype matrix |
| **Covariate adjustment (age, sex, PCs)** | Any analysis without covariate adjustment is methodologically incomplete; standard since 2015 | Medium | Covariate file loading + null model fitting shared across tests |
| **FDR multiple testing correction** | Benjamini-Hochberg is the default in field; already in `gene_burden.py` | Low | Already implemented — carry over from existing code |
| **Bonferroni correction** | Conservative alternative required for reporting | Trivial | Already implemented |
| **Per-gene results output (TSV)** | Every published tool outputs TSV per gene with p-values, effect sizes, variant counts | Low | Results consolidation module |
| **Effect size + confidence intervals** | OR + 95% CI for logistic burden; beta + SE for linear; SKAT outputs test statistic only | Medium | statsmodels results objects |
| **Binary trait support** | Case-control is the primary design in nephrology cohorts | Medium | Logistic null model |
| **Variant count per gene in output** | Users need n_variants, n_cases_carriers, n_controls_carriers to assess results | Low | Aggregate from genotype matrix |
| **MAF filter for included variants** | Standard: only include variants below MAF threshold (e.g., 1%, 5%) | Low | Applied during genotype matrix construction |

---

## Standard Output Columns (Table Stakes)

Based on SAIGE-GENE+ and seqMeta output conventions (HIGH confidence — from official documentation):

```
GENE  n_variants  n_cases  n_controls  case_carriers  control_carriers
fisher_p  fisher_or  fisher_ci_lower  fisher_ci_upper
burden_p  burden_beta  burden_se  burden_or  burden_or_ci_lower  burden_or_ci_upper
skat_p  skat_stat
skat_o_p  skat_o_rho_opt  skat_o_stat
acat_o_p
corrected_fisher_p  corrected_skat_o_p  corrected_acat_o_p
```

The seqMeta convention of `test_p`, `test_beta`, `test_se` as prefixed column groups is the
de facto standard. SAIGE-GENE+ uses `Pvalue_Burden`, `Pvalue_SKAT`, `Pvalue` (SKAT-O overall)
plus `markerIDs`, `markerAFs`.

---

## Differentiators

Features that set this framework apart from a basic implementation. Not universally expected
but provide genuine scientific and usability value.

| Feature | Value Proposition | Complexity | Dependencies |
|---------|-------------------|------------|--------------|
| **Dual SKAT backend (R via rpy2 + pure Python)** | R SKAT is the published gold standard; Python fallback removes R requirement for users who don't have R; iterative validation ensures correctness | High | rpy2 (optional), numpy, ctypes wrapper for Davies C code (`qfc.c` from SKAT source) |
| **ACAT-O omnibus test** | Combines burden + SKAT + ACAT-V without correlation estimation; robust across genetic architectures; Liu et al. 2019 | Medium | Requires p-values from other tests; Cauchy combination formula is simple analytically |
| **ACAT-V (per-variant Cauchy combination)** | Particularly powerful when only a few causal variants are present; combines marginal score test p-values | Medium | Per-variant score test p-values from SKAT null model |
| **Quantitative trait support (linear burden + linear SKAT)** | Many nephrology phenotypes are continuous (eGFR, proteinuria); without this, framework is half-functional | Medium | `statsmodels.api.OLS` (already dependency); linear SKAT score test |
| **Functional variant weights (CADD, REVEL)** | DeepRVAT (Nature Genetics 2024) shows functional annotations improve power; REVEL already in variantcentrifuge pipeline columns | Medium | Score columns already extracted by SnpSift; weight normalization module |
| **PCA file integration** | Population stratification control; 10 PCs standard practice; support PLINK .eigenvec, AKT, and generic TSV formats | Medium | PCA loading + sample ID matching + covariate merging |
| **Lambda_GC diagnostic** | Genomic inflation factor; expected in any analysis submission; immediately reveals miscalibration | Low | Compute from observed p-values: median(chi2) / 0.4549 |
| **QQ plot data output** | Expected in any methodology paper; output TSV of observed vs expected -log10(p) — plotting is optional | Low | Sort p-values, compute expected quantiles |
| **JSON config mode** | Power users want reproducible analysis configs as files; lab SOP-able | Low | Config dataclass serialization; no new dependencies |
| **Allelic series test (COAST)** | 29-82% more significant associations than SKAT-O in UK Biobank (AJHG 2023); uses BMV/DMV/PTV variant categories already present in variantcentrifuge annotation | High | Requires variant classification into benign missense / damaging missense / PTV; PolyPhen + SIFT annotations must be available; `AllelicSeries` R package or re-implementation of COAST |
| **Uniform weight option** | Needed for comparison/sensitivity analysis; also backward-compatible with current behavior | Trivial | Already implied by current Fisher implementation |

---

## Test Suite Detail (What's Expected and Why)

### Standard Burden Tests

| Test | What It Does | When Best | Covariate Support | Output |
|------|-------------|-----------|-------------------|--------|
| **Fisher's exact** (CMC/CAST) | 2x2 contingency table; collapse variants per gene | Strong Mendelian signals, small N, no covariates | No | OR, 95% CI, p |
| **Logistic regression burden** | Weighted sum → binary outcome via logistic regression | Case-control with covariates (age, sex, PCs) | Yes | OR, 95% CI, p (Wald or LRT) |
| **Linear regression burden** | Weighted sum → continuous outcome via OLS | Quantitative traits with covariates | Yes | beta, SE, p (t-test) |
| **Madsen-Browning / WSS** | Wilcoxon rank-sum on weighted burden score; permutation p-values | Protective variants present; non-parametric | Partial | p (permutation) |

Madsen-Browning is historically important but rarely used standalone in 2024 — most tools
default to logistic regression burden. Including WSS as an option adds completeness without
complexity (scipy has Wilcoxon). LOW confidence on current adoption: single WebSearch finding.

### Variance-Component Tests

| Test | Key Parameters | Output |
|------|---------------|--------|
| **SKAT** | kernel (linear.weighted default), weights.beta (1,25 default), missing_cutoff (0.15), max_maf (1.0) | Q statistic, p-value via Davies/Liu |
| **SKAT-O** | rho grid (0 to 1, 11 values), optimal.adj correction, method="SKATO" | p-value, optimal rho |

### Omnibus Tests

| Test | What It Combines | When Best |
|------|-----------------|-----------|
| **ACAT-O** | ACAT-V + SKAT (two weight configs) + Burden (two weight configs) = 6 tests | Robust across unknown genetic architectures; recommended default by REGENIE |
| **SKATO-ACAT** | SKAT-O via Cauchy combination (faster than integration) | When SKAT-O integration is slow |

REGENIE (v3+) supports both SKATO and ACATO as named tests — the REGENIE documentation
confirms ACATO is now a peer option alongside SKATO, not a replacement (MEDIUM confidence
from official REGENIE documentation).

### Allelic Series Tests

| Test | Variant Categories | Weighting | Status |
|------|-------------------|-----------|--------|
| **COAST** | BMV (benign missense), DMV (damaging missense), PTV (protein-truncating) | Default w=(1,2,3); customizable | Published AJHG 2023; R package on CRAN |
| **COAST-SS** | Summary statistics version | Same | Extends COAST; AJHG 2025 |

COAST is genuinely new (2023) and not yet universally available in tools — this would be
differentiated. However, it requires PolyPhen + SIFT variant classification, which
variantcentrifuge already extracts via SnpSift annotations. HIGH confidence on publication
findings from PMC full-text.

---

## Variant Weighting Schemes

| Scheme | Formula | When to Use | Standard In |
|--------|---------|-------------|-------------|
| **Beta(MAF; 1, 25)** | `scipy.stats.beta.pdf(maf, 1, 25)` | Default for rare variant analyses | R SKAT, SAIGE-GENE, REGENIE, all major tools |
| **Beta(MAF; 1, 1)** | Uniform; identical weight for all variants | Sensitivity/comparison | Same tools, as "flat" option |
| **CADD-normalized** | `CADD_phred / 40` (normalized to ~[0,1]) | When CADD scores available in annotations | RAVA-FIRST (PLOS Genetics 2022) |
| **REVEL-based** | REVEL score directly (already 0-1) | Missense-rich gene sets; REVEL already in pipeline | Research tool; growing adoption |
| **Combined** | `Beta(MAF; 1, 25) × CADD_normalized` | When both functional and frequency signal matter | DeepRVAT (Nature Genetics 2024) framework |

HIGH confidence on Beta(1,25) as dominant standard — confirmed in original SKAT paper (PMC),
SKAT R package documentation (RDocumentation), and multiple tool manuals.

---

## PCA / Population Stratification

| Aspect | Standard Practice | Notes |
|--------|------------------|-------|
| **Number of PCs** | 10 PCs as default; 3 PCs often sufficient for homogeneous cohorts | Including 10 PCs reduces spurious associations to nominal levels (PMC6283567) |
| **Pre-computed vs on-the-fly** | Pre-computed is strongly preferred; PCA should use common variants not the rare variants being tested | Rare variants alone are insufficient for PCA — must use common SNPs (MAF > 5%) |
| **File format** | PLINK `.eigenvec` (tab-separated, FID/IID/PC1..PCn) is most common; AKT output is similar | Simple TSV parser handles both |
| **PCA tools** | PLINK 2.0 (most common); AKT (Illumina, VCF-native); EIGENSOFT | Design doc recommends supporting both PLINK and AKT |
| **When to compute on-the-fly** | Only when user does not have pre-computed PCA and has AKT/PLINK available | Should be optional pipeline stage, not required |

MEDIUM confidence — cross-referenced PMC articles on population stratification, consistent
across multiple sources.

---

## Diagnostics

| Diagnostic | What It Shows | Implementation | Expected by |
|-----------|--------------|----------------|-------------|
| **Lambda_GC** | Genomic inflation factor; lambda ≈ 1.0 expected; >1.1 suggests stratification | `median(observed_chi2) / 0.4549` | Any reviewer; any paper |
| **QQ plot data (TSV)** | Observed vs expected -log10(p); inflation visible | Sort p-values, compute expected quantiles | Standard in any GWAS/RVAS paper |
| **QQ plots (PNG/SVG)** | Visual output of QQ data | matplotlib (optional dependency) | Nice-to-have; TSV data mandatory |
| **Per-test lambda_GC** | Separate lambda per test (Fisher, SKAT-O, ACAT-O) | Compute from each p-value vector | Allows diagnosis of test-specific inflation |
| **Permutation p-values** | Non-parametric p-values when sample size too small for asymptotic assumptions | Computationally expensive; adaptive permutation possible | Small-sample situations (<50 cases) |

**On permutation-based p-values:** When minor allele count (MAC) is very low, analytical
p-values from SKAT/burden tests have inflated type I error (PMC, biostatistics literature).
Permutation is the correct solution but expensive. For lab-scale data (hundreds to thousands
of samples), this is relevant. However, implementing adaptive permutation is high complexity
and should be deferred — flag low-MAC genes in output as a mitigation.

---

## Anti-Features

Features to explicitly NOT build. Common mistakes or scope creep in this domain.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Meta-analysis (RAREMETAL, Meta-SAIGE)** | Different problem: combining results across cohorts; requires score statistics export format; out of scope for a single-cohort tool | Output score statistics in standard format if users want to feed into RAREMETAL — but do not implement meta-analysis itself |
| **Mixed model / GRM (SAIGE-GENE approach)** | Correct for large relatedness in biobank cohorts; requires genetic relationship matrix computation, REML fitting; adds weeks of implementation complexity | For the GCKD nephrology cohort, covariate PCs + kinship exclusion is sufficient; recommend excluding close relatives as a preprocessing step |
| **Phased haplotype tests** | Compound heterozygous detection is already handled by inheritance analysis in the existing pipeline; don't duplicate | Use existing `comp_het.py` output as a covariate or stratification variable if needed |
| **Deep learning weights (DeepRVAT)** | DeepRVAT (Nature Genetics 2024) trains neural network weights per phenotype; requires massive training data; no off-the-shelf model for nephrology | Support CADD/REVEL as functional weights; do not build DeepRVAT |
| **Conditional analysis** | SAIGE-GENE+ implements conditional analysis (testing after conditioning on known signals); complex, requires additional iteration | Not needed for discovery analyses; add to future roadmap if a reviewer requests it |
| **Manhattan plots (full genome)** | This is a gene-burden tool, not a GWAS tool; gene-level Manhattan (one dot per gene) is less standard | Output QQ plot data; let users make their own plots from the TSV |
| **LD computation within variant sets** | LD between rare variants within a gene is negligible and computing it adds complexity for no gain | Proceed without LD computation; SKAT's variance-component approach implicitly handles correlation |
| **Adaptive test selection per gene** | Automatically choosing best test per gene post-hoc inflates type I error | Pre-specify tests; report all; let user interpret |
| **Real-time progress visualization** | Nice for interactive use but adds dependency (tqdm already present) and complicates pipeline logging | Use existing logging; progress visible via log output |
| **REGENIE-style whole-genome split** | REGENIE splits analysis into step 1 (null model on common variants) + step 2 (rare variant tests); designed for biobank scale (>100K samples) | Not needed for lab cohorts (<10K samples) where standard logistic/SKAT null model is fast |

---

## Feature Dependencies

```
Existing Fisher's exact (gene_burden.py)
    |
    +--> Refactor into association/tests/fisher.py   [Step 1: Foundation]
    |    Verify bit-identical output
    |
    +--> AssociationEngine + abstract base classes    [Step 1: Foundation]
    |    + CLI extensions (backward compatible)
    |
    +--> Covariate system (loading, validation)       [Step 2]
    |    + Beta weight implementation
    |    + Uniform weight implementation
    |    |
    |    +--> Logistic regression burden test         [Step 2]
    |    |    (statsmodels.api.Logit — already a dep)
    |    |
    |    +--> Linear regression burden test           [Step 2]
    |         (statsmodels.api.OLS — already a dep)
    |
    +--> SKAT null model fitting                      [Step 3]
    |    (shared across SKAT, SKAT-O, ACAT-O)
    |    |
    |    +--> SKAT score test + Davies/Liu p-value    [Step 3 or 4]
    |    |    (R backend via rpy2 preferred;
    |    |     Python fallback: numpy eigenvalues
    |    |     + Davies ctypes from qfc.c)
    |    |
    |    +--> SKAT-O (rho grid search)                [Step 3 or 4]
    |    |    (requires SKAT; optimal.adj correction)
    |    |
    |    +--> ACAT-V (per-variant Cauchy)             [Step 5]
    |    |    (requires per-variant score test p-values from null model)
    |    |
    |    +--> ACAT-O (omnibus combination)            [Step 5]
    |         (requires SKAT p-value + burden p-value + ACAT-V p-value)
    |
    +--> Diagnostics (lambda_GC, QQ data)             [Step 5]
    |    (requires all p-values collected)
    |
    +--> PCA file loading                             [Step 6]
    |    (pre-computed PCA → covariates)
    |
    +--> Functional weights (CADD, REVEL)             [Step 7 / Future]
    |    (requires score columns in variant data)
    |
    +--> Allelic series test (COAST)                  [Step 7 / Future]
         (requires BMV/DMV/PTV classification;
          PolyPhen + SIFT already in annotations)
```

---

## MVP Recommendation

For a credible first release that substantially improves on Fisher-only:

### MVP (Must Have — Steps 1-3)
1. **Refactor Fisher** into association framework (backward compatible)
2. **Logistic regression burden** with covariate adjustment (age, sex, PCs from file)
3. **SKAT + SKAT-O** via R backend (rpy2) with pure Python fallback
4. **Beta(1,25) weights** as default
5. **FDR + Bonferroni** (already done — carry over)
6. **Lambda_GC** diagnostic output
7. **QQ data TSV** output
8. **Standard TSV output** (per-gene, all tests, effect sizes, CIs)

### Post-MVP (High Value, Steps 4-6)
- ACAT-O omnibus test (analytically simple once SKAT/burden p-values exist)
- ACAT-V (per-variant Cauchy combination)
- PCA computation stage (AKT/PLINK wrappers)
- Linear burden test for quantitative traits
- Functional weights (CADD, REVEL)
- JSON config mode

### Future (Steps 7+)
- Allelic series test (COAST) — requires BMV/DMV/PTV classification
- Permutation-based p-values for low-MAC situations
- Kinship matrix support / mixed model correction

---

## Small-Sample Considerations

The GCKD cohort has ~5,000 samples. Most nephrology candidate-gene studies have 50-500
cases. Key implications:

| Sample Size | Risk | Mitigation |
|-------------|------|-----------|
| N_cases < 50 | SKAT analytical p-values may be inflated for binary traits; Liu/Davies approximations less accurate | Flag genes with case_carriers < 10; consider exact/permutation p-values |
| N_cases 50-200 | SKAT-O "optimal.adj" method handles this better than "optimal" | Use method="optimal.adj" as default in R SKAT |
| N_cases > 500 | Standard asymptotic approximations valid | Default settings appropriate |

The SKAT R package `SKAT_Null_Model_MomentAdjust()` function provides small-sample
adjustment for binary traits — this should be the default when R backend is used.

---

## Complexity Ratings

| Rating | Definition |
|--------|------------|
| **Trivial** | < 1 hour; single formula or lookup |
| **Low** | 2-8 hours; well-understood algorithm, clear reference implementation |
| **Medium** | 1-3 days; multiple components, testing required |
| **High** | 1-2 weeks; numerical methods, validation required, multiple failure modes |

---

## Confidence Assessment

| Area | Confidence | Basis |
|------|------------|-------|
| SKAT/SKAT-O parameters and behavior | HIGH | Official R package documentation (RDocumentation), original SKAT paper (PMC3135811) |
| ACAT-O mechanism | HIGH | Full-text PMC article (PMC6407498) |
| Beta(1,25) as dominant weight | HIGH | SKAT paper, SKAT R package docs, confirmed in multiple tools |
| SAIGE-GENE+ output format | HIGH | Official SAIGE documentation (saigegit.github.io) |
| REGENIE supported tests | HIGH | Official REGENIE documentation (rgcgithub.github.io) |
| Allelic series (COAST) | HIGH | PMC full-text (PMC10432147) |
| PCA standard practice (10 PCs) | MEDIUM | Multiple PMC articles, consistent finding |
| Madsen-Browning current adoption rate | LOW | Single WebSearch; not verified in current tool documentation |
| ACAT-O replacing SKAT-O | LOW | No evidence found; they appear complementary; REGENIE offers both |
| Functional weight adoption in 2025 publications | MEDIUM | DeepRVAT (Nature Genetics 2024) verified; general CADD/REVEL usage confirmed |

---

# v0.16.0 Addendum: Association Hardening and Multi-Cohort Features

**Researched:** 2026-02-23
**Scope:** Features NEW to v0.16.0. The sections above document the v0.15.0 foundation (already built). This section documents what v0.16.0 adds on top.

---

## v0.16.0 Feature 1: Weighted FDR Correction (Issue #86)

### How Comparable Tools Handle It

The standard genomic approach has two distinct layers:

**Layer A — Static weighted BH (Genovese et al. 2006)**
- Assigns per-hypothesis non-negative weights `w_i` such that `sum(w_i) == n_hypotheses`.
- Computes weighted p-values: `p_i_weighted = p_i / w_i`.
- Applies standard BH to weighted p-values.
- FDR is controlled at nominal alpha when `sum(w) == n`.
- Used when: weights are pre-specified from biological priors (gene-level GWAS signal, pathway membership, clinical gene list).
- Reference: Genovese, Roeder & Wasserman (2006); PMC3458887

**Layer B — IHW (Ignatiadis et al. Nature Methods 2016)**
- Data-driven weight assignment using a covariate independent of p-values under H0 but informative of power.
- Splits hypotheses into strata by covariate, subdivides into cross-validation folds to prevent overfitting.
- Learns per-stratum weights; applies weighted BH.
- Used when: a good covariate exists (MAF, gene expression, pathway score) and cohort is large enough for cross-validation (typically >1000 tests).
- Available as Bioconductor `IHW` R package. No mature Python equivalent.
- Reference: PMC4930141; Nature Methods 10.1038/nmeth.3885

**Applied to rare variant association (Huang et al. Genes 2022, PMC8872452):**
- Gene-level ML prediction score used as IHW covariate. Concentrates FDR budget on high-prior genes.
- Reported 109% more discoveries at FDR alpha=0.3 for schizophrenia data.
- Recommended covariate types: gene expression in disease-relevant tissue, GWAS locus membership, pathway risk scores, gene-level constraint (pLI/LOEUF).

### Current State in VariantCentrifuge

`correction.py` implements flat BH and Bonferroni only. No weighting mechanism exists.

### Table Stakes for v0.16.0

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Static weighted BH (Genovese) | Standard approach; `statsmodels.stats.multitest.multipletests` accepts `weights` kwarg | Low | ~10-20 lines added to `correction.py`. Verified: statsmodels >= 0.14 supports `weights` parameter. |
| Per-gene weight file loader | Mechanism to supply biological prior weights | Low-Med | Add `--fdr-weight-file` flag; TSV format (gene, weight); normalize to sum=n internally before passing to `multipletests` |

### Differentiators for v0.16.0

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| IHW integration (optional, R) | Up to 2x discovery improvement when good covariate exists; state-of-art in biobank studies | Medium | Wrap IHW R package via rpy2. Add `--fdr-method ihw --ihw-covariate <col>` where covariate comes from the output table (e.g., n_variants, mean_maf, or a user-supplied column). |
| Built-in covariate diagnostic | Show covariate-to-pvalue correlation so users know whether IHW will help | Low | Spearman correlation of each numeric output column vs p-values; flagged in diagnostics output |

### Anti-Features for v0.16.0

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Mandatory biological priors | Most users won't have them; blocks basic usage | Make weights optional; default is flat BH (unchanged behavior) |
| IHW from scratch in Python | Active research area; subtle cross-validation correctness requirements | Wrap IHW R package via rpy2; graceful fallback to weighted BH if R/IHW absent |
| Separate FDR per test type | Running BH separately for SKAT vs burden pools inflates discovery | Apply correction jointly across all genes for each test type, or use omnibus p-value as the correction target |

**Dependencies:** `correction.py` `apply_correction()` — add `weights` parameter. `gene_burden.py` — thread optional weight vector. No new external tools required for weighted BH.

---

## v0.16.0 Feature 2: Sample-Level Case-Confidence Weights (Issue #85)

### How Comparable Tools Handle It

**REGENIE and SAIGE:** Neither natively supports case-confidence weights for uncertain binary phenotypes. They treat binary phenotype as hard 0/1. Standard workarounds:
1. **Liability threshold model** — convert binary phenotype to continuous liability score (e.g., from EHR probability), then run as quantitative trait. Recent application: Nature Genetics 2025 (10.1038/s41588-025-02370-4).
2. **Super-normal control design** — tighten control inclusion criteria to remove borderline cases. Shown to recover power under moderate-to-high misclassification (medRxiv 2025.12.14).

**Weighted logistic regression:**
- Technically sound: down-weight samples with uncertain diagnosis using `sample_weight` in the logistic regression.
- `statsmodels.Logit` does NOT support per-sample weights. `sklearn.LogisticRegression` supports `sample_weight` in `fit()`.
- For SKAT the weighting modifies the score vector: `S = G^T * (w * residuals)` — requires modifying null model residuals and eigenvalue recomputation. Research-grade modification, not off-the-shelf.
- No established reference implementation in a published genomics tool as of research date.

**Current state in VariantCentrifuge:** Binary phenotype only, hard 0/1. No weighting path.

### Table Stakes for v0.16.0

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Weighted logistic burden (binary path) | Enables case-confidence in burden tests; sklearn has clean support | Medium | Use `sklearn.LogisticRegression(solver='lbfgs', sample_weight=...)` for weighted path. Add `--sample-weight-file` flag (TSV: sample_id, weight). |
| Liability/continuous proxy mode documentation | Best-practice workaround when SKAT weighting is not available | Low | Already possible via `--trait-type quantitative`; needs documentation. |

### Differentiators for v0.16.0

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Weighted SKAT (per-sample weights in score) | Statistically correct handling of uncertain phenotypes in kernel tests | High | Modify `S = G^T * (sqrt(w) * residuals)` and recompute null model with weighted GLM. Flag as experimental. Requires numerical validation. |
| Sample-weight QC report | ESS (effective sample size = (sum w)^2 / sum w^2); warns when ESS << N | Low-Med | Standard in survey statistics; prevents silent power loss |

### Anti-Features for v0.16.0

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Conflating sample weights with FDR weights | These are different layers (data vs. correction) | Keep sample-weight file and FDR-weight file as separate CLI flags |
| Frequency-based sample weighting | No theoretical basis in phenotype uncertainty | Use case-confidence directly; MAF weights belong to variant weighting |
| Modifying phenotype vector in-place | Destroys auditability | Keep original phenotype array; pass weights as separate vector |

**Dependencies:** `logistic_burden.py` — add `sample_weight` parameter; `python_backend.py` `fit_null_model()` — add weighted GLM path (deferred for SKAT). New CLI flag `--sample-weight-file`.

---

## v0.16.0 Feature 3: COAST Classification Fix and Scoring Model (Issue #87)

### How the Reference Implementation Classifies Variants

COAST (McCaw et al. AJHG 2023, PMC10432147) uses three canonical categories:

**PTV (Protein Truncating Variant):**
- Effect types: `stop_gained`, `frameshift_variant`, `splice_acceptor_variant`, `splice_donor_variant`
- Filtered by LOFTEE to remove low-confidence LoF

**DMV (Damaging Missense Variant):**
- Effect: `missense_variant`
- PolyPhen-2 "probably damaging" OR "possibly damaging"
- AND SIFT "deleterious" OR "deleterious low confidence"
- Note: the original paper uses OR logic for at least one tool predicting damaging

**BMV (Benign Missense Variant):**
- Effect: `missense_variant`
- PolyPhen-2 "benign" AND SIFT "tolerated" OR "tolerated low confidence"
- Both tools must agree on benign

**CADD/REVEL:** NOT used in the reference COAST classification scheme. The `AllelicSeries` CRAN vignette confirms annotation codes are user-supplied; COAST does not dictate classification methodology. CADD and REVEL are discussed as potential extensions for continuous severity.

### Current State in VariantCentrifuge

`allelic_series.py` `classify_variants()` is **correctly implemented** matching the reference:
- PTV: HIGH impact + PTV effects (correct)
- DMV: missense + OR(SIFT-damaging, PolyPhen-damaging) (correct)
- BMV: missense + AND(SIFT-benign, PolyPhen-benign) (correct)

**The real COAST bug** is NOT classification. Root causes of `p=None` for all genes:
1. **Strict 3-category gate** (lines 551-581 of `allelic_series.py`): genes missing any one of BMV/DMV/PTV return `None`. In real cohorts most genes lack one category.
2. **gene_df alignment**: `contingency_data` may not carry `gene_df` correctly through the Python COAST path, triggering `coast_skip_reason: NO_GENE_DF` or `ANNOTATION_GENOTYPE_MISMATCH`.

### Table Stakes for v0.16.0

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Configurable BMV/DMV/PTV classification criteria | COAST annotation is user-supplied by design; different annotation pipelines differ | Medium | Add `--coast-classification-config` JSON: define which effect types map to which category, which SIFT/PolyPhen values are damaging/benign |
| Partial-category COAST fallback | Most real genes have only 2 of 3 categories; returning None for all is not useful | Medium | When only 2 categories present: run 2-category allelic series test or drop the missing category and test with available ones |
| CADD score as optional DMV classifier | Single-source classification when PolyPhen/SIFT unavailable; CADD>20 is standard DMV threshold in the field | Low-Med | Add `--coast-cadd-threshold`; classify missense as DMV if CADD>=threshold, BMV if below |
| REVEL as optional DMV classifier | REVEL outperforms SIFT+PolyPhen for missense pathogenicity; REVEL>=0.5 is standard threshold | Low-Med | Add `--coast-revel-threshold`; same structure as CADD |
| gene_df wiring audit and fix | Core bug: genotype matrix build must preserve gene_df row alignment | High | Audit why `ANNOTATION_GENOTYPE_MISMATCH` occurs; fix row alignment between `gene_df` and genotype matrix columns in `engine.py` + `genotype_matrix.py` |

### Differentiators for v0.16.0

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| N-category COAST support | AllelicSeries accepts arbitrary integer annotation levels; enable >3 categories for fine-grained severity | Medium | Allow user to map annotation values to integer codes 1..k; python backend handles arbitrary k |
| LOFTEE flag integration | Removes low-confidence LoF from PTVs; reduces false positives in PTV category | Low | Add LOFTEE column detection; exclude LOFTEE_HC=LC variants from PTV |

### Anti-Features for v0.16.0

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Hard 3-category gate as only COAST mode | Returns p=None for most real-world genes; defeats the purpose | Return None only when zero classifiable variants; use partial-category test for 1-2 categories |
| Recomputing SIFT/PolyPhen from sequence | Tool is post-annotation analysis framework | Accept pre-computed annotations from SnpEff/dbNSFP columns |
| CADD-only classification replacing SIFT/PolyPhen | Changes results vs. reference; poor for reproducibility | Treat CADD as supplementary fallback when SIFT/PolyPhen absent, not primary classifier |

**Dependencies:** `allelic_series.py` `classify_variants()` + `coast_python.py` `PythonCOASTBackend` + `engine.py` `contingency_data` assembly + `genotype_matrix.py` row ordering.

---

## v0.16.0 Feature 4: Region Restriction in Prefilter (Issue #84)

### How Comparable Tools Handle It

BED-file region restriction is a standard prefilter step in all VCF pipelines:

**bcftools native approach:**
- `bcftools view --regions-file <bed> input.vcf.gz` requires bgzipped + tabix-indexed VCF
- BED coordinates are 0-based half-open; VCF coordinates are 1-based — bcftools handles this conversion
- Known issue: `--regions-file` can produce duplicates for overlapping BED intervals (bcftools GitHub issue #221); use `--regions-file` not `--targets-file` for exact intersection without memory issues
- Streaming: if BED is bgzipped + tabix-indexed (.tbi), bcftools streams efficiently; uncompressed BED is loaded into memory entirely

**Cross-cohort use case:**
- Restricting to shared genomic regions captured in ALL cohorts reduces false negatives from differential capture
- Standard in multi-cohort analysis: compute intersection of capture kits, apply as prefilter BED

### Current State in VariantCentrifuge

`BCFtoolsPrefilterStage` applies expression-based filtering (`bcftools view -i "..."`) but does NOT support `--regions-file`. Gene extraction uses a gene BED for bcftools view (coordinate-based) but this is the gene-selection BED, not a user-supplied region restriction.

### Table Stakes for v0.16.0

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| `--regions-bed` CLI flag | Standard cross-cohort analysis requirement; bcftools native, trivial to add | Low | Prepend `bcftools view -R <bed>` (or `--regions-file <bed>`) before expression filter; can combine with existing `--bcftools-prefilter` |
| BED format validation | bcftools is finicky about 0-based coords; silent errors are common | Low | Validate BED has 3 columns; warn if chromosomes don't match VCF contigs |
| Tabix index detection | Streaming via tabix is faster for large BEDs | Low | Check `<bed>.tbi` exists; log whether streaming or in-memory mode |

### Differentiators for v0.16.0

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Capture kit intersection utility | Compute BED intersection of N cohort capture kits before analysis; enables reproducible multi-cohort runs | Medium | Wrapper around `bedtools multiinter -i *.bed | awk '$4==N'`; document as separate utility script |

### Anti-Features for v0.16.0

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Merging region BED with gene BED internally | Creates confusing behavior; output genes depend on intersection ordering | Apply region BED as pre-step; gene BED as subsequent step; keep logically separate |
| Requiring tabix-indexed BED | Not all users have tabix installed | Support both; index just speeds things up |
| Region restriction on output TSV (post-extraction) | Very wasteful; must be applied at VCF-level | Apply at bcftools view level, before SnpSift |

**Dependencies:** `BCFtoolsPrefilterStage` — add `--regions-file` flag to bcftools invocation. `extract_variants()` in `filters.py` — thread regions-file through. No new external tools; bcftools already required.

---

## v0.16.0 Feature 5: Single Eigendecomposition for SKAT-O Optimization

### How Reference Implementations Handle It

SKAT-O requires eigendecomposition of the kernel matrix once per gene. The dominant cost is kernel matrix formation and projection-adjusted Z1 computation.

**SAIGE-GENE+ optimization (Nature Genetics 2022):** 1,407x speedup achieved by:
1. Sparse genotype matrix (see Feature 6)
2. Reusing eigenvalues across Burden/SKAT/SKAT-O — all tests share the same kernel
3. Precomputing projection-adjusted Z1_adj once, reusing across the rho grid

**Current VariantCentrifuge Python backend:**
- `_get_lambda()` in `python_backend.py` called per rho value in SKAT-O via `_skato_optimal_param()`.
- Z1_half decomposition computed at SKAT-O call time; not shared with standard SKAT call.

### Table Stakes for v0.16.0

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Profiling to confirm bottleneck | Must verify eigendecomposition is actually the bottleneck before optimizing | Low | Add `--profile-association` timing flag; check whether kernel formation or eigendecomposition dominates at n=5125 |

### Differentiators for v0.16.0

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Shared Z1_adj across rho values | Avoids recomputing projection within the rho grid loop | Low-Med | Refactor `_skato_optimal_param()` to accept pre-projected Z1_adj; verify numerical equivalence |
| Shared null model residuals across SKAT/SKAT-O/COAST | Null model fit is O(n^3) once; already shared — confirm wiring is correct | Low | Audit whether `fit_null_model()` is called once per cohort vs. once per gene |

### Anti-Features for v0.16.0

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Caching eigenvalues across different genes | Gene matrices have different shapes; would be wrong | Cache only within a single gene's multi-test run |
| Pre-computing all gene kernels before testing | Memory O(n_genes * n_variants^2); infeasible genome-wide | Compute per-gene, free after test |

**Dependencies:** `python_backend.py` `test_gene()` and `_skato_optimal_param()` — internal refactor. `coast_python.py` `_run_allelic_skat()` — shares `_skato_get_pvalue` infrastructure.

---

## v0.16.0 Feature 6: Sparse Genotype Matrices

### When Sparsity Pays Off

Rare variant genotype matrices are inherently sparse: MAF < 1-5% means >95% of entries are 0 (hom-ref). After imputation, missing values are imputed to `round(2*MAF)` = 0 when MAF < 0.25 (always true for rare variants), preserving sparsity.

**Measured density figures from literature:**
- At n=50,000 WES samples, 96% of loci have ALT allele frequency < 0.1% (spVCF paper, Bioinformatics 2021)
- SKAT with MAF < 1%: expected density approximately 0.01 (1% non-zero entries)
- scipy sparse matrix breaks even around 20-25% density; clearly beneficial below 10%

**SAIGE-GENE+ (Nature Genetics 2022):** Sparse genotype matrix reduced memory from 65 GB to 2.1 GB for TTN gene in UK Biobank.

**R SKAT package:** Supports sparse matrix format (confirmed in SKAT 2.0+ CRAN PDF, July 2025).

**Practical threshold for GCKD (n=5,125):**
- Typical gene: 5-50 rare variants at MAF < 1%; density ~0.01-0.02
- scipy CSR sparse IS beneficial for memory at this scale
- For a 5125x50 matrix at 1% density: ~2MB dense vs ~0.3MB sparse
- Sparse format breakeven: `n_samples * n_variants > ~50,000` cells AND density < 20%
- At GCKD scale, sparse is a memory nicety; at UK Biobank scale it is essential

### Table Stakes for v0.16.0

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| scipy.sparse CSR option in `build_genotype_matrix` | Memory efficiency for large cohorts; R SKAT already supports it | Medium | Add `sparse=True` parameter; return `scipy.sparse.csr_matrix`. Most numpy ops need `.toarray()` — profile overhead before full commitment |
| Sparse-aware SKAT kernel computation | `K = G.T @ G` works identically for sparse or dense scipy matrices | Low-Med | Verify `scipy.sparse.csr_matrix @ G` numerical equivalence with dense path via tests |

### Differentiators for v0.16.0

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Automatic sparsity detection | Avoid sparse overhead for small genes; use dense for <50 variants, sparse for >50 at large n | Low | `if n_samples * n_variants > THRESHOLD and density < 0.20: use_sparse = True` |
| Memory usage diagnostics | Show peak genotype matrix memory per gene; guides whether sparse optimization is warranted | Low | Log at DEBUG: `n_samples x n_variants, density=X%, format=sparse/dense` |

### Anti-Features for v0.16.0

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Forcing sparse for all genes | For small genes (5-10 variants), sparse overhead exceeds benefit; adds conversion cost | Threshold-based auto-selection based on density * size |
| Persisting sparse matrices to disk between stages | The entire point of sparse is fast in-memory ops; serialization adds overhead | Keep as in-memory objects within the gene loop |
| Sparse format in COAST burden tests | COAST aggregates by category using small dense ops; sparse adds complexity without benefit | Dense in COAST; sparse only for SKAT kernel computation |

**Dependencies:** `genotype_matrix.py` `build_genotype_matrix()` — add optional sparse return. `python_backend.py` — SKAT kernel `Z1 @ Z1.T` and score `G.T @ residuals` can be sparse-aware. `coast_python.py` — do NOT convert to sparse.

---

## v0.16.0 Feature 7: PCAComputationStage Wiring (AKT/PLINK Integration)

### Current State

`pca.py` `load_pca_file()` reads pre-computed PCA files in three formats (PLINK `.eigenvec` with/without header, AKT stdout). Must be provided via `--pca-file`. No stage computes PCA.

### What Comparable Tools Do

- **REGENIE:** Requires pre-computed PCA as covariates; users run PLINK2 `--pca` separately
- **SAIGE-GENE+:** Provides a PCA step as a separate tool (`createSparseGRM.R`)
- **STAARpipeline:** Accepts pre-computed GRM/PCA
- **Convention:** PCA is computed once per cohort because it requires common variants (MAF > 5%), not the rare variants being tested

### Table Stakes for v0.16.0

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Documentation of PLINK2 PCA workflow | Most users will use pre-computed PCA; clear instructions needed | Low | Document in CLI help: which `plink2 --pca` command to run, what `--pca-file` format to use |
| Optional PCAComputationStage (plink2 wrapper) | Removes need for manual PLINK2 invocation when a VCF + sample list is available | High | Requires separate common-variant VCF input (rare variant VCF is wrong for PCA). Add `--pca-compute-vcf` flag. Calls `plink2 --pca`. |

### Anti-Features for v0.16.0

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| PCA from rare variant VCF | Rare variants are poor for PCA (sparse, no LD structure) | Require separate common-variant VCF or pre-computed PC file |
| Implementing PCA computation in Python | Well-solved by PLINK2 | Wrap `plink2 --pca` via subprocess |

**Dependencies:** `pca.py` `load_pca_file()` — already correct. New external tool dependency: `plink2` (optional; fails gracefully). New CLI flag `--pca-compute-vcf`.

---

## v0.16.0 MVP Recommendation

### Must Ship (Table Stakes — fix real bugs and enable primary use cases)

1. **COAST gene_df alignment fix** — Core bug causing p=None; targeted fix to engine/contingency_data wiring
2. **COAST partial-category fallback** — Eliminates the "all None" problem for real cohorts
3. **COAST configurable classification** — Allow CADD/REVEL when SIFT/PolyPhen unavailable
4. **Weighted BH (static, Genovese)** — Low-complexity, high-value; statsmodels already supports it
5. **Per-gene weight file** — Input mechanism for #4; simple TSV loader
6. **`--regions-bed` prefilter** — Essential for cross-cohort analysis; bcftools-native, trivial

### Defer to Post-MVP (Differentiators — harder or need validation)

1. **IHW integration** — R/rpy2 dependency, cross-validation complexity; separate sub-milestone
2. **Weighted SKAT (per-sample phenotype weights)** — Research-grade; needs numerical validation
3. **Sparse genotype matrices** — Memory nicety at GCKD scale; defer until cohort grows or performance becomes a bottleneck
4. **PCAComputationStage** — Documentation gap, not a bug; users can run PLINK2 manually
5. **Single eigendecomposition optimization** — Profile first; current code may already be fine

---

## v0.16.0 Confidence Assessment

| Feature Area | Confidence | Basis |
|--------------|------------|-------|
| Weighted BH (Genovese) | HIGH | PMC3458887, statsmodels docs verified |
| IHW covariate approach | HIGH | PMC4930141, PMC8872452, Bioconductor IHW vignette |
| COAST classification reference | HIGH | PMC10432147, CRAN vignette confirming user-supplied annotations |
| COAST p=None root cause | HIGH | Direct code reading of allelic_series.py lines 551-581 |
| bcftools --regions-file behavior | HIGH | Official bcftools docs; known issue #221 verified |
| Sparse matrix breakeven density | MEDIUM | spVCF paper; SKAT R package docs; general CSR literature |
| Case-confidence weighted SKAT | LOW | No reference implementation found; theoretical extension only |
| SKAT-O eigendecomposition sharing | MEDIUM | SAIGE-GENE+ paper; direct code reading of python_backend.py |

---

## Sources (v0.16.0 Addendum)

- Genovese, Roeder & Wasserman (2006) weighted BH: [PMC3458887](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458887/)
- Ignatiadis et al. IHW (Nature Methods 2016): [PMC4930141](https://pmc.ncbi.nlm.nih.gov/articles/PMC4930141/) | [Nature Methods](https://www.nature.com/articles/nmeth.3885)
- Huang et al. gene-level IHW for rare variants (Genes 2022): [PMC8872452](https://pmc.ncbi.nlm.nih.gov/articles/PMC8872452/)
- Bioconductor IHW package vignette: [bioconductor.org/packages/IHW](https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html)
- McCaw et al. COAST (AJHG 2023): [PMC10432147](https://pmc.ncbi.nlm.nih.gov/articles/PMC10432147/)
- AllelicSeries CRAN vignette (annotation codes user-supplied): [cran.r-project.org](https://cran.r-project.org/web/packages/AllelicSeries/vignettes/coast.html)
- SAIGE-GENE+ sparse matrix performance (Nature Genetics 2022): [Nature Genetics](https://www.nature.com/articles/s41588-022-01178-w)
- spVCF sparse genotype matrix paper (Bioinformatics 2021): [Oxford Academic](https://academic.oup.com/bioinformatics/article/36/22-23/5537/6029516)
- bcftools regions-file documentation: [samtools.github.io](https://samtools.github.io/bcftools/bcftools.html)
- REGENIE overview and binary trait handling: [rgcgithub.github.io/regenie](https://rgcgithub.github.io/regenie/overview/)
- Meta-SAIGE multi-cohort analysis (Nature Genetics 2025): [Nature Genetics](https://www.nature.com/articles/s41588-025-02403-y)
- Phenotype misclassification SuperControl (medRxiv 2025): [medRxiv](https://www.medrxiv.org/content/10.64898/2025.12.14.25342213v1.full)
- Liability threshold model GWAS (Nature Genetics 2025): [Nature Genetics](https://www.nature.com/articles/s41588-025-02370-4)
- STAARpipeline functional annotation categories (Nature Methods 2022): [Nature Methods](https://www.nature.com/articles/s41592-022-01641-w)
- Weighted multiple testing in genomics (BioData Mining 2012): [PMC3458887](https://biodatamining.biomedcentral.com/articles/10.1186/1756-0381-5-4)

## Sources (v0.15.0 baseline, carried forward)

- [SKAT R Package Documentation — RDocumentation](https://www.rdocumentation.org/packages/SKAT/versions/2.0.1/topics/SKAT) — HIGH
- [SKAT original paper: Wu et al. 2011 — PMC3135811](https://pmc.ncbi.nlm.nih.gov/articles/PMC3135811/) — HIGH
- [SKAT-O paper: Lee et al. 2012 — PMC3415556](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/) — HIGH
- [ACAT paper: Liu et al. 2019 — PMC6407498](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/) — HIGH
- [COAST allelic series paper — PMC10432147](https://pmc.ncbi.nlm.nih.gov/articles/PMC10432147/) — HIGH
- [SAIGE-GENE+ Nature Genetics 2022](https://www.nature.com/articles/s41588-022-01178-w) — HIGH
- [REGENIE gene-based test options — official docs](https://rgcgithub.github.io/regenie/options/) — HIGH
- [Rare-Variant Association Analysis review — PMC4085641](https://pmc.ncbi.nlm.nih.gov/articles/PMC4085641/) — HIGH
- [RVTESTS GitHub and documentation](https://github.com/zhanxw/rvtests) — MEDIUM
- [Population stratification in rare variant analysis — PMC6283567](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283567/) — MEDIUM
- [RAVA-FIRST (CADD-based weighting) — PLOS Genetics](https://pmc.ncbi.nlm.nih.gov/articles/PMC9518893/) — MEDIUM
- [DeepRVAT — Nature Genetics 2024](https://www.nature.com/articles/s41588-024-01919-z) — MEDIUM
- [AllelicSeries R package — CRAN](https://cran.r-project.org/package=AllelicSeries) — HIGH
- [seqMeta output format conventions — GitHub genepi-freiburg/seqmeta](https://github.com/genepi-freiburg/seqmeta) — MEDIUM
