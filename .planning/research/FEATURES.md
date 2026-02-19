# Feature Landscape — Modular Rare Variant Association Framework

**Domain:** Rare variant association testing for genomic sequencing studies
**Researched:** 2026-02-19
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

## Sources

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
