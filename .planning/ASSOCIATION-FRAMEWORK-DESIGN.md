# Modular Rare Variant Association Framework — Design Document

## Motivation

The current gene burden pipeline uses Fisher's exact test only — adequate for strong
Mendelian signals (PKD1 p=2e-80) but underpowered for discovering novel associations.
Modern rare variant association studies require SKAT-O, covariate adjustment, variant
weighting, population stratification control, and omnibus combination tests.

**Goal**: Build a modular, all-Python statistical testing framework that lets users
choose what data they have and what tests to run, from simple Fisher's exact to
full SKAT-O with covariates and PCA.

---

## AKT vs PLINK for PCA

| Feature | AKT (Illumina) | PLINK 2.0 |
|---------|----------------|-----------|
| **Input format** | VCF/BCF directly | Requires conversion to BED/BGEN |
| **PCA quality** | Identical to PLINK for top PCs (AKT paper, Bioinformatics 2017) |  Reference implementation |
| **Speed (2,504 samples)** | 47 seconds | ~similar, but add conversion time |
| **Kinship** | Built-in (`akt kin`) | `--make-king` |
| **WGS workflow** | Native — uses tabix on VCF | Requires VCF→PLINK conversion |
| **Site selection** | Ships with curated common SNP VCFs | User must filter |
| **Dependencies** | C++, htslib | C/C++ |

**Recommendation**: Support both. AKT is ideal for users already working with
multi-sample VCFs (like GCKD). PLINK is the standard when users have pre-existing
BED files or want KING kinship matrices. The framework should accept pre-computed
PCA files from either tool, plus optionally run AKT or PLINK internally.

**References**:
- AKT paper: https://academic.oup.com/bioinformatics/article/33/1/142/2525685
- PLINK 2.0: https://www.cog-genomics.org/plink/2.0/

---

## Architecture: Modular Association Engine

### Design Principles

1. **User chooses what they have** — covariates? PCA? pedigree? phenotype type?
2. **User chooses what tests to run** — Fisher, logistic regression burden, SKAT-O, ACAT-O
3. **Dual backend** — R SKAT via rpy2 (gold standard) + pure Python fallback (no R needed)
4. **Backward compatible** — existing `--perform-gene-burden` still works identically
5. **Pipeline-integrated** — PCA computation is a pipeline stage, not a separate tool
6. **Iterative validation** — Python implementation validated against R results, improved until convergence

### Module Layout

```
variantcentrifuge/
  association/                    # NEW package
    __init__.py
    engine.py                     # AssociationEngine: orchestrates test selection + execution
    config.py                     # AssociationConfig dataclass, CLI parsing helpers
    tests/                        # Statistical test implementations
      __init__.py
      base.py                     # Abstract AssociationTest class
      fisher.py                   # Fisher's exact (current, refactored)
      logistic_burden.py          # Logistic regression burden test
      skat.py                     # SKAT (variance-component score test)
      skat_o.py                   # SKAT-O (optimal rho combination)
      acat.py                     # ACAT-V and ACAT-O omnibus
      allelic_series.py           # Allelic series test (LoF vs damaging vs benign)
    weights/                      # Variant weighting schemes
      __init__.py
      base.py                     # Abstract VariantWeighter
      beta_weights.py             # Beta(MAF; a1, a2) — Madsen-Browning / Wu weights
      functional_weights.py       # CADD, REVEL, PolyPhen-based
      uniform.py                  # Current behavior (all weight = 1)
    covariates.py                 # Covariate loading, validation, PCA integration
    pca.py                        # PCA computation (AKT wrapper, PLINK wrapper, or load from file)
    kinship.py                    # Kinship matrix (AKT kin, KING from PLINK)
    correction.py                 # Multiple testing (FDR, Bonferroni, Holm) — refactored from gene_burden.py
    diagnostics.py                # QQ plot, lambda_GC, Manhattan plot data
    results.py                    # Unified result format, effect sizes, CI computation
    backends/                     # SKAT/SKAT-O computation backends
      __init__.py
      base.py                     # Abstract SKATBackend class
      r_backend.py                # R SKAT via rpy2 (gold standard, preferred when R available)
      python_backend.py           # Pure Python implementation (fallback, iteratively validated)
      davies.py                   # Davies method: ctypes wrapper for qfc.c + Liu moment-matching fallback
```

### Core Abstractions

```python
# association/tests/base.py

class AssociationTest(ABC):
    """Base class for all gene-level association tests."""

    name: str                    # e.g., "fisher", "skat_o", "logistic_burden"
    supports_covariates: bool    # Can this test adjust for covariates?
    supports_weights: bool       # Can this test use variant weights?
    supports_quantitative: bool  # Can this test handle continuous traits?
    requires_genotype_matrix: bool  # Needs full G matrix (not just counts)?

    @abstractmethod
    def test_gene(
        self,
        genotype_matrix: np.ndarray,    # (n_samples, n_variants) — 0/1/2 dosages
        phenotype: np.ndarray,           # (n_samples,) — binary or continuous
        covariates: np.ndarray | None,   # (n_samples, n_covariates) or None
        weights: np.ndarray | None,      # (n_variants,) or None
        **kwargs,
    ) -> TestResult:
        """Run association test for one gene. Return TestResult."""
        ...

class TestResult:
    """Standardized result from any association test."""
    gene: str
    test_name: str
    p_value: float
    effect_size: float          # OR for binary, beta for continuous
    effect_ci_lower: float
    effect_ci_upper: float
    n_variants: int
    n_cases: int
    n_controls: int
    case_carriers: int
    control_carriers: int
    statistic: float            # Test statistic (chi2, Q, T, etc.)
    extra: dict                 # Test-specific extras (rho_opt for SKAT-O, etc.)
```

### AssociationEngine Flow

```
User config (CLI/JSON)
    │
    ▼
┌──────────────────────────────────────────┐
│ AssociationEngine.configure()            │
│                                          │
│  1. Load phenotype (binary/continuous)   │
│  2. Load covariates (age, sex, ...)      │
│  3. Load PCA (from file or compute)      │
│  4. Merge covariates + PCA              │
│  5. Select tests based on user config    │
│  6. Select weighting scheme              │
│  7. Validate compatibility               │
└──────────────────────────────────────────┘
    │
    ▼
┌──────────────────────────────────────────┐
│ For each gene:                           │
│                                          │
│  1. Extract genotype matrix (n × m)      │
│  2. Apply variant weights                │
│  3. Run each selected test               │
│  4. Collect TestResults                  │
│  5. If ACAT-O: combine p-values          │
└──────────────────────────────────────────┘
    │
    ▼
┌──────────────────────────────────────────┐
│ Post-processing:                         │
│                                          │
│  1. Multiple testing correction (FDR)    │
│  2. Compute diagnostics (λ_GC, QQ)      │
│  3. Format output (TSV + optional plots) │
└──────────────────────────────────────────┘
```

---

## Statistical Tests: Mathematical Specifications

### 1. Fisher's Exact Test (existing — refactor into module)

2×2 contingency table, two-tailed. No covariates, no weights.
Already implemented in `gene_burden.py`. Move to `association/tests/fisher.py`.

### 2. Logistic Regression Burden Test

For binary traits with covariates:

```
logit(P(Y=1)) = α + β·B + γ₁·C₁ + γ₂·C₂ + ... + γ_k·C_k

where:
  B = Σⱼ wⱼ · Gⱼ     (weighted burden score, summed over variants j in gene)
  Cₖ = covariates (age, sex, PC1, PC2, ..., PC_n)
  wⱼ = variant weight (Beta(MAF; 1, 25) or CADD-based or uniform)

Test: Wald test or LRT on β
P-value: from χ²(1) distribution
Effect size: OR = exp(β), CI from profile likelihood or Wald
```

Implementation: `statsmodels.api.Logit` — already a dependency.

### 3. Linear Regression Burden Test

For quantitative traits with covariates:

```
Y = α + β·B + γ₁·C₁ + ... + γ_k·C_k + ε

Test: t-test on β
P-value: from t(n-p-1) distribution
Effect size: β, CI from OLS
```

Implementation: `statsmodels.api.OLS` — already a dependency.

### 4. SKAT (Variance-Component Score Test)

**Dual-backend**: R SKAT via rpy2 (gold standard) with pure Python fallback.

**Null model** (fit once, reuse across all genes):

```
Binary:   logit(P(Y=1)) = α + Γ·C        (logistic)
Continuous:         Y    = α + Γ·C + ε    (linear)
```

**Score statistic** for gene with genotype matrix G (n × m) and weight matrix W:

```
Q_SKAT = (Y - μ̂)ᵀ · G · W · Gᵀ · (Y - μ̂)

where:
  μ̂ = fitted values from null model
  W = diag(w₁², w₂², ..., wₘ²)     (variant weights squared)
  G = n × m genotype matrix (0/1/2)
```

**P-value**: Q follows a mixture of χ² distributions under H₀:

```
Q ~ Σᵢ λᵢ · χ²(1)

where λᵢ are eigenvalues of:
  P₀^(1/2) · G · W · Gᵀ · P₀^(1/2)
  P₀ = V - V·X·(Xᵀ·V·X)⁻¹·Xᵀ·V    (projection matrix)
  V = diag(μ̂ᵢ(1 - μ̂ᵢ)) for binary, σ²·I for continuous
```

P-value via **Davies' method** (exact) or **moment-matching** to a scaled χ²
(faster, sufficiently accurate):

```
Method: Match first 3 moments of Σ λᵢ·χ²(1) to a·χ²(d) + b
  a = Var(Q) / (2·E[Q])
  d = 2·E[Q]² / Var(Q)
  b = E[Q] - a·d
  p = P(a·χ²(d) + b > Q_observed)
```

Dependencies: `numpy` (eigenvalues), `scipy.stats.chi2` (p-value). Already available.

### 5. SKAT-O (Optimal ρ Combination)

Combines burden and SKAT by searching over ρ ∈ {0, 0.1, 0.2, ..., 1.0}:

```
Q(ρ) = (1-ρ)·Q_SKAT + ρ·Q_burden

where:
  Q_burden = (Y - μ̂)ᵀ · G · W · 1 · 1ᵀ · W · Gᵀ · (Y - μ̂)
           = [Σⱼ wⱼ(Y-μ̂)ᵀGⱼ]²
```

For each ρ, compute p(ρ) via the mixture-of-χ² method above.
Then:

```
p_SKAT-O = min_ρ p(ρ), corrected for searching over ρ
```

Correction uses the correlation structure between Q(ρ) values
(Lee et al. 2012, "optimal.adj" method from SKAT R package).

### SKAT/SKAT-O Implementation Strategy: Dual Backend

**Why dual backend?** The SKAT-O algorithm is mathematically non-trivial —
especially the Davies method for mixture-of-χ² p-values and the ρ-correction
in SKAT-O. Rather than build from scratch and hope for numerical correctness,
we use R as the oracle and iteratively improve the Python implementation.

#### Backend Architecture

```python
class SKATBackend(ABC):
    """Abstract backend for SKAT/SKAT-O computation."""

    @abstractmethod
    def fit_null_model(self, phenotype, covariates, trait_type): ...

    @abstractmethod
    def test_gene(self, genotype_matrix, weights, method="SKAT-O") -> TestResult: ...

class RSKATBackend(SKATBackend):
    """Gold standard: calls R SKAT package via rpy2."""
    # Uses: SKAT::SKAT_Null_Model(), SKAT::SKAT()
    # Exact same algorithm as published papers

class PythonSKATBackend(SKATBackend):
    """Pure Python fallback, iteratively validated against R."""
    # Uses: numpy, scipy, compiled Davies via ctypes
```

#### Backend Selection (Automatic)

```python
def get_skat_backend() -> SKATBackend:
    """Auto-detect best available backend."""
    try:
        import rpy2
        # Check R and SKAT package available
        return RSKATBackend()
    except ImportError:
        logger.info("R/rpy2 not available, using Python SKAT backend")
        return PythonSKATBackend()
```

Users can also force a backend: `--skat-backend r` or `--skat-backend python`.

#### R Backend Details

The R SKAT package ([leelabsg/SKAT](https://github.com/leelabsg/SKAT)) contains:

| R file | What it does | Python equivalent |
|--------|-------------|-------------------|
| `R/SKAT_Optimal.R` | SKAT-O rho search + p-value combination | `skat_o.py` |
| `R/SKAT_CompQuadForm.R` | Davies method wrapper + Liu fallback | `backends/davies.py` |
| `R/Function.R` | `Get_Lambda`, `Beta.Weights`, `Get_Liu_Params_Mod` | `backends/python_backend.py` |
| `R/MAIN.R` | Entry point, method dispatch | `skat.py` |
| `R/SKAT_Linear.R` | Score test for continuous traits | `backends/python_backend.py` |
| `R/SKAT_Logistic.R` | Score test for binary traits | `backends/python_backend.py` |
| `R/Null_Model.R` | Null model fitting | `backends/python_backend.py` |
| `src/qfc.cpp` | **C++ Davies algorithm** (Algorithm AS 155) | `backends/davies.py` (ctypes) |

The R package is ~54% R + ~38% C++ (performance-critical Davies method in compiled code).

#### Davies Method: The Critical Piece

The Davies method computes exact p-values for Q ~ Σλᵢχ²(1). The R SKAT
package uses `qfc.cpp` (~300 lines of C++, Algorithm AS 155). Our strategy:

1. **Primary**: Compile `qfc.c` (standalone, no R dependency) and call via `ctypes`
2. **Fallback**: Liu moment-matching approximation (pure Python, less accurate for
   very small p-values but adequate for most genes)
3. **Validation**: Compare ctypes Davies vs R Davies on identical eigenvalue sets

The `qfc.c` file is self-contained C code implementing Davies (1980) "Algorithm AS 155:
The distribution of a linear combination of χ² random variables". It can be compiled
as a shared library without any R infrastructure.

```python
# backends/davies.py — simplified concept

def davies_pvalue(q_stat: float, eigenvalues: np.ndarray) -> float:
    """Compute P(Q > q) via Davies method."""
    try:
        # Try compiled C implementation first
        return _davies_ctypes(q_stat, eigenvalues)
    except (OSError, RuntimeError):
        # Fall back to Liu moment-matching approximation
        return _liu_moment_matching(q_stat, eigenvalues)

def _liu_moment_matching(q_stat: float, eigenvalues: np.ndarray) -> float:
    """Liu et al. approximation — pure Python."""
    c1 = np.sum(eigenvalues)           # E[Q]
    c2 = np.sum(eigenvalues**2)        # Var(Q)/2
    c3 = np.sum(eigenvalues**3)
    c4 = np.sum(eigenvalues**4)
    # ... moment matching to scaled chi-squared
    # Ported from R/Function.R: Get_Liu_Params_Mod()
```

#### Iterative Validation Workflow

```
Step 1: Implement R backend via rpy2
        └── Run on GCKD data → reference p-values for all genes

Step 2: Implement Python backend (score test + Liu approximation)
        └── Compare against R p-values
        └── Identify genes where |log10(p_python) - log10(p_R)| > threshold

Step 3: Add compiled Davies method (ctypes)
        └── Re-compare — should close gap on small p-values
        └── Target: p-values within 10% relative difference

Step 4: Port SKAT-O rho correction from R/SKAT_Optimal.R
        └── Compare SKAT-O p-values
        └── Use R source code as reference for "optimal.adj" correction

Step 5: Once Python matches R within tolerance → Python becomes default
        └── R backend kept as validation/comparison option
```

#### R Source Code Reference

The R SKAT package source is fully available for study:

- **Repository**: https://github.com/leelabsg/SKAT (Lee Lab, primary)
- **Mirror**: https://github.com/lin-lab/SKAT (Lin Lab)
- **CRAN**: https://cran.r-project.org/package=SKAT
- **Key paper**: Lee S et al. (2012) "Optimal tests for rare variant effects
  in sequencing association studies." Biostatistics 13(4):762-775

The C++ Davies implementation (`src/qfc.cpp`) can be extracted and compiled
independently — it has no R-specific dependencies beyond the `.C()` call interface,
which we replace with a standard C function signature for ctypes.

### 6. ACAT-O (Omnibus P-Value Combination)

Combines p-values from different tests using Cauchy distribution:

```
T_ACAT = Σᵢ wᵢ · tan((0.5 - pᵢ)·π)

p_ACAT = 0.5 - arctan(T_ACAT / Σwᵢ) / π
```

Typically combines: p_burden, p_SKAT, p_ACAT-V (sparse signal test).
No permutations needed. Robust to correlated p-values.

**ACAT-V** (variant-level): Combines per-variant score test p-values:

```
T_ACAT-V = Σⱼ wⱼ · tan((0.5 - pⱼ)·π)
```

where pⱼ is the marginal score test p-value for variant j.

### 7. Allelic Series Test (Future)

Tests whether variant classes (LoF, damaging missense, benign missense) have
monotonically ordered effects:

```
H₀: β_LoF = β_dmis = β_bmis = 0
H₁: |β_LoF| ≥ |β_dmis| ≥ |β_bmis|   (ordered alternatives)
```

---

## Variant Weighting Schemes

### Beta Weights (Madsen-Browning / Wu et al.)

```
w(MAF) = Beta(MAF; a₁, a₂)

Standard:   a₁=1, a₂=25   → strongly up-weights rare variants
Flat:       a₁=1, a₂=1    → uniform weights (= current behavior)
SKAT default: a₁=1, a₂=25
```

### Functional Weights

```
w(v) = f(CADD(v), REVEL(v), PolyPhen(v), ...)

Options:
  - CADD-based:  w = CADD_phred / 40  (normalized to ~[0,1])
  - REVEL-based:  w = REVEL_score
  - Combined:  w = Beta(MAF) × CADD_normalized
```

---

## PCA Integration

### Supported Modes

```
--pca-file eigenvec.txt             # Pre-computed (any tool)
--pca-tool akt --pca-n-components 20   # Compute with AKT
--pca-tool plink --pca-n-components 20 # Compute with PLINK 2.0
```

### AKT PCA Workflow (Pipeline Stage)

```bash
# 1. Extract common biallelic SNPs (MAF > 5%, LD-pruned sites shipped with AKT)
akt pca -R akt_sites_GRCh37.vcf.gz input.vcf.gz > pca_results.txt

# 2. Optionally compute kinship
akt kin -R akt_sites_GRCh37.vcf.gz input.vcf.gz > kinship.txt
```

AKT outputs eigenvectors + eigenvalues. The framework parses these and merges
with sample phenotype/covariate data.

### PLINK PCA Workflow

```bash
# 1. Convert VCF to PLINK format (if not already available)
plink2 --vcf input.vcf.gz --make-bed --out temp_plink

# 2. LD-prune
plink2 --bfile temp_plink --indep-pairwise 50 5 0.2 --out pruned

# 3. Compute PCA
plink2 --bfile temp_plink --extract pruned.prune.in --pca 20 --out pca_results
```

### PCA File Format (Accepted Input)

Tab-separated, first column = sample ID, remaining columns = PC1, PC2, ...:

```
#SAMPLE_ID  PC1         PC2         PC3        ...
Sample001   0.00234     -0.01456    0.00089    ...
Sample002   -0.00567    0.00234     0.00345    ...
```

Supports: PLINK `.eigenvec`, AKT output, EIGENSOFT `.evec`, custom TSV.

---

## Covariate System

### CLI Interface

```
--covariate-file covariates.tsv      # TSV with sample ID + covariate columns
--covariates age,sex,bmi             # Which columns to use (default: all)
--pca-file pca.eigenvec              # PCA components (merged as covariates)
--pca-components 10                  # How many PCs to include (default: 10)
```

### Covariate File Format

```
#SAMPLE_ID  age     sex     bmi     batch    center
Sample001   52      1       24.5    batch1   center_A
Sample002   67      0       31.2    batch2   center_B
```

- Numeric columns used directly
- Categorical columns (string) automatically one-hot encoded
- Missing values: warn and exclude sample, or impute to median

---

## User-Facing CLI Design

### Simple Mode (Backward Compatible)

```bash
# Existing behavior — works exactly as before
variantcentrifuge ... --perform-gene-burden --gene-burden-mode alleles
```

### Extended Mode

```bash
variantcentrifuge ... \
  --perform-gene-burden \
  --association-tests fisher,skat_o,acat_o \
  --gene-burden-mode alleles \
  --covariate-file covariates.tsv \
  --covariates age,sex \
  --pca-file pca.eigenvec \
  --pca-components 10 \
  --variant-weights beta \
  --variant-weight-params "a1=1,a2=25" \
  --correction-method fdr \
  --association-output gene_association_results.tsv \
  --diagnostics-output diagnostics/
```

### JSON Config Mode (Power Users)

```bash
variantcentrifuge ... --association-config association_config.json
```

```json
{
  "tests": ["fisher", "logistic_burden", "skat_o", "acat_o"],
  "trait_type": "binary",
  "covariates": {
    "file": "covariates.tsv",
    "columns": ["age", "sex", "bmi"],
    "pca_file": "pca.eigenvec",
    "pca_components": 10
  },
  "weights": {
    "scheme": "beta",
    "params": {"a1": 1, "a2": 25}
  },
  "correction": "fdr",
  "diagnostics": true,
  "output": {
    "per_test_results": true,
    "combined_results": true,
    "qq_plot": true,
    "lambda_gc": true
  }
}
```

---

## Output Format

### Combined Results TSV (one row per gene, all tests)

```
GENE  n_variants  n_cases  n_controls  case_carriers  control_carriers  \
  fisher_p  fisher_or  fisher_or_ci_lower  fisher_or_ci_upper  \
  burden_p  burden_beta  burden_se  \
  skat_p  skat_stat  \
  skat_o_p  skat_o_rho_opt  skat_o_stat  \
  acat_o_p  \
  corrected_fisher_p  corrected_skat_o_p  corrected_acat_o_p
```

### Diagnostics Output

- `lambda_gc.txt` — genomic inflation factor per test
- `qq_data.tsv` — observed vs expected -log10(p) for QQ plots
- Optionally: PNG/SVG plots via matplotlib (if available)

---

## Implementation Steps

### Step 1: Core Framework + Fisher Refactor (Foundation)
- Abstract test/weight/result classes
- Refactor Fisher's exact into `association/tests/fisher.py`
- AssociationEngine orchestrator
- CLI extensions (backward compatible)
- Tests: verify Fisher refactor produces identical output

### Step 2: Covariate System + Logistic Regression Burden
- Covariate loading system
- Logistic regression burden test (statsmodels)
- Linear regression burden test (statsmodels)
- Variant weighting (Beta weights, uniform)
- Tests: validate against manual logistic regression

### Step 3: R SKAT Backend (Gold Standard)
- rpy2 integration for R SKAT package
- Null model fitting via `SKAT::SKAT_Null_Model()`
- SKAT/SKAT-O/Burden via `SKAT::SKAT()` with method dispatch
- Automatic R/rpy2 availability detection
- Tests: run on synthetic data with known signals, store reference p-values

### Step 4: Pure Python SKAT Backend (Iterative)
- Port score test from `R/SKAT_Linear.R` and `R/SKAT_Logistic.R`
- Compile `qfc.c` as shared library, wrap with ctypes for Davies method
- Liu moment-matching fallback (from `R/Function.R: Get_Liu_Params_Mod`)
- Port SKAT-O rho search from `R/SKAT_Optimal.R`
- Port "optimal.adj" correction from `R/SKAT_Optimal.R`
- **Validate**: compare all p-values against R backend on same data
- **Iterate**: fix discrepancies using R source as reference until convergence

### Step 5: ACAT-O + Diagnostics
- ACAT-V (per-variant score combination)
- ACAT-O (omnibus: Fisher + burden + SKAT + ACAT-V)
- Lambda_GC computation
- QQ plot data generation
- Multiple testing correction across all tests

### Step 6: PCA Integration
- PCA file loader (PLINK, AKT, generic formats)
- AKT wrapper (optional pipeline stage)
- PLINK wrapper (optional pipeline stage)
- PCA as covariates in regression/SKAT
- Kinship matrix support (future: mixed models)

### Step 7: Advanced Features (Future)
- Functional variant weights (CADD, REVEL)
- Quantitative trait support (linear regression, linear SKAT)
- Allelic series test
- JSON config mode
- Documentation and examples

---

## Dependencies

### Already Available
- `scipy.stats` — Fisher's exact, chi2 distribution
- `statsmodels` — Logit, OLS, Table2x2, multipletests
- `numpy` — linear algebra, eigenvalues
- `pandas` — data handling

### New (Optional)
- `rpy2` — R SKAT backend (preferred for SKAT/SKAT-O when R is available)
- `matplotlib` — diagnostic plots (optional, graceful fallback)

### External Tools (Optional)
- `R` + `SKAT` package — gold standard SKAT/SKAT-O (used via rpy2)
- `akt` — PCA/kinship from VCF (user's existing tool)
- `plink2` — PCA/kinship (alternative)

### Compiled (Bundled)
- `qfc.c` — Davies method for mixture-of-χ² p-values (extracted from SKAT,
  compiled as shared library, called via ctypes). Fallback to Liu approximation
  if compilation not available.

---

## Validation Strategy

### Reference Datasets
1. **Synthetic data**: Known ground truth (genes with injected burden signal)
2. **R SKAT backend as oracle**: Run R backend on same data → store reference p-values
3. **Python vs R comparison**: Automated test suite comparing Python backend against R on
   identical inputs across hundreds of genes
4. **Published results**: Reproduce known gene-disease associations from public data
5. **GCKD PKD1/PKD2**: Verify new tests still detect PKD1 (p<1e-50) and PKD2

### Acceptance Criteria
- Fisher refactor: **bit-identical** output to current implementation
- R SKAT backend: **exact match** to running R SKAT directly (same code path)
- Python SKAT vs R SKAT: p-values within **10% relative difference** (numerical precision)
- Python SKAT-O vs R SKAT-O: p-values within **10% relative difference**
- Davies (ctypes) vs R Davies: p-values within **1% relative difference** (same algorithm)
- ACAT-O: exact match to analytical formula (deterministic)
- Lambda_GC on null data: within **[0.95, 1.05]** (well-calibrated)
- All existing gene burden tests unchanged (backward compatibility)

### Iterative Improvement Process
1. Run both backends on GCKD cohort (all genes)
2. Identify genes where `|log10(p_python) - log10(p_R)| > 0.1`
3. Diagnose: score test, eigenvalues, Davies, or rho-correction?
4. Fix Python implementation using R source as reference
5. Re-run comparison → iterate until all genes converge

---

## References

### Papers
1. Wu MC et al. (2011) SKAT. Am J Hum Genet 89(1):82-93
2. Lee S et al. (2012) SKAT-O. Biostatistics 13(4):762-775
3. Liu Y et al. (2019) ACAT. Am J Hum Genet 104(3):410-421
4. Li X et al. (2020) STAAR. Nat Genet 52(4):437-445
5. Li B, Leal SM (2008) CMC. Am J Hum Genet 83(3):311-321
6. Madsen BE, Browning SR (2009) Variable threshold. PLoS Genet 5(3):e1000384
7. Arthur R et al. (2017) AKT. Bioinformatics 33(1):142-144
8. Lavallee-Adam M et al. (2023) Allelic series. Am J Hum Genet 90(7)
9. Davies RB (1980) Algorithm AS 155. J R Stat Soc C (Applied Statistics) 29(3):323-333

### Source Code (for Python port)
- R SKAT package (Lee Lab): https://github.com/leelabsg/SKAT
- R SKAT package (Lin Lab): https://github.com/lin-lab/SKAT
- CRAN documentation: https://cran.r-project.org/package=SKAT
- Key files for porting:
  - `R/SKAT_Optimal.R` — SKAT-O rho search + optimal.adj correction
  - `R/SKAT_CompQuadForm.R` — Davies wrapper + Liu fallback
  - `R/Function.R` — Get_Lambda, Get_Liu_Params_Mod, Beta.Weights
  - `R/SKAT_Logistic.R` — binary trait score test
  - `R/SKAT_Linear.R` — continuous trait score test
  - `R/Null_Model.R` — null model fitting
  - `src/qfc.cpp` — C++ Davies algorithm (AS 155, ~300 lines, standalone)
- Python reference: https://github.com/cozygene/RL-SKAT (RL-SKAT, uses FaST-LMM Davies)
- Moment matching: https://github.com/deanbodenham/momentchi2py (4 approximation methods)
