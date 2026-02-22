# Phase 26: Association Testing Documentation - Research

**Researched:** 2026-02-22
**Domain:** Technical documentation (Sphinx/MyST Markdown) for genomics software
**Confidence:** HIGH

## Summary

Phase 26 is a documentation-only phase requiring no code changes. The association testing framework (Phases 18-25) is fully implemented and the documentation task is to write comprehensive user-facing content in the existing Sphinx + MyST Markdown system.

The existing docs infrastructure is mature: Sphinx with the Furo theme, MyST parser for Markdown, autodoc for API reference, sphinxcontrib.mermaid for diagrams. All existing guides in `docs/source/guides/` are plain Markdown files. The new guide belongs at `docs/source/guides/association_testing.md`.

The association framework is production-ready with 8 test variants (Fisher, LogisticBurden, LinearBurden, SKAT/SKAT-O via Python or R backends, COAST via Python or R backends, ACAT-V internal to SKAT, ACAT-O omnibus), a full covariate system, PCA integration, functional variant weights, and three diagnostics files. Every CLI flag, config key, and output column is verifiable directly from the source code.

**Primary recommendation:** Write all content by reading source code directly — no "stale from training" risk. The source is the ground truth for flag names, defaults, output columns, and behavior.

## Documentation Infrastructure

### Format and Location (HIGH confidence — directly verified)

| Item | Value |
|------|-------|
| Format | MyST Markdown (`.md` files) throughout `docs/source/` |
| Build system | Sphinx with Furo theme |
| Extensions | `myst_parser`, `sphinx.ext.autodoc`, `sphinx.ext.napoleon`, `sphinxcontrib.mermaid`, `sphinx_autodoc_typehints` |
| New guide location | `docs/source/guides/association_testing.md` |
| API stubs location | `docs/source/api/association.md` (new) |
| Index update | `docs/source/index.md` toctree under "Practical Guides" |
| Guides index update | `docs/source/guides/index.md` toctree |
| API index update | `docs/source/api/index.md` toctree under new "Association Testing" section |

### Existing API Stub Pattern (HIGH confidence — verified from `docs/source/api/gene_burden.md`)

```markdown
# Gene Burden Module

Statistical methods for gene burden analysis

```{eval-rst}
.. automodule:: variantcentrifuge.gene_burden
   :members:
   :undoc-members:
   :show-inheritance:
```
```

The association module API stub should follow this exact pattern, once per submodule.

## Association Module: Public API (HIGH confidence — verified from source)

### Classes exported from `variantcentrifuge.association` (`__init__.py`)

```python
from variantcentrifuge.association import (
    AssociationConfig,   # Configuration dataclass
    AssociationEngine,   # Orchestrator
    AssociationTest,     # Abstract base class
    TestResult,          # Per-gene result dataclass
    apply_correction,    # Standalone FDR/Bonferroni function
)
```

### `AssociationConfig` Fields (all with defaults)

| Field | Default | Description |
|-------|---------|-------------|
| `correction_method` | `"fdr"` | `"fdr"` (Benjamini-Hochberg) or `"bonferroni"` |
| `gene_burden_mode` | `"samples"` | `"samples"` (carrier-based CMC/CAST) or `"alleles"` |
| `confidence_interval_method` | `"normal_approx"` | CI computation method |
| `confidence_interval_alpha` | `0.05` | CI significance level (0.05 = 95% CI) |
| `continuity_correction` | `0.5` | Haldane-Anscombe correction for zero cells |
| `covariate_file` | `None` | Path to covariate TSV/CSV (first col = sample ID) |
| `covariate_columns` | `None` | Subset of columns; None = all |
| `categorical_covariates` | `None` | Columns to one-hot encode; None = auto-detect |
| `trait_type` | `"binary"` | `"binary"` (logistic) or `"quantitative"` (OLS) |
| `variant_weights` | `"beta:1,25"` | Weight scheme string |
| `variant_weight_params` | `None` | Extra weight parameters (dict) |
| `missing_site_threshold` | `0.10` | Exclude variants with >10% missing site-wide |
| `missing_sample_threshold` | `0.80` | Exclude samples with >80% missing |
| `firth_max_iter` | `25` | Newton-Raphson iterations for Firth fallback |
| `skat_backend` | `"python"` | `"python"`, `"r"` (deprecated), or `"auto"` |
| `skat_method` | `"SKAT"` | `"SKAT"`, `"Burden"`, or `"SKATO"` |
| `coast_backend` | `"python"` | `"python"`, `"r"` (deprecated), or `"auto"` |
| `min_cases` | `200` | Warn if n_cases < this |
| `max_case_control_ratio` | `20.0` | Warn if controls/cases > this |
| `min_case_carriers` | `10` | Per-gene warning if case_carriers < this |
| `diagnostics_output` | `None` | Path to diagnostics directory |
| `pca_file` | `None` | Pre-computed PCA file path |
| `pca_tool` | `None` | `"akt"` to invoke AKT as subprocess |
| `pca_components` | `10` | Number of PCs to use; warn if >20 |
| `coast_weights` | `None` | Category weights for COAST [BMV, DMV, PTV]; default [1,2,3] |

### `AssociationEngine` Methods

```python
# Construction (preferred)
engine = AssociationEngine.from_names(["fisher", "skat"], config)

# Run all tests, returns wide DataFrame
result_df = engine.run_all(gene_burden_data)  # list of dicts
```

### `AssociationTest` Abstract Interface

```python
class AssociationTest(ABC):
    @property
    @abstractmethod
    def name(self) -> str: ...          # "fisher", "logistic_burden", etc.

    @abstractmethod
    def run(self, gene, contingency_data, config) -> TestResult: ...

    def effect_column_names(self) -> dict: ...   # column suffix mapping
    def check_dependencies(self) -> None: ...    # raise ImportError if missing
    def prepare(self, gene_count: int) -> None:  # pre-loop hook
    def finalize(self) -> None: ...              # post-loop hook
```

### Concrete Test Classes (for API stubs)

| Class | Module | Registry Name |
|-------|--------|---------------|
| `FisherExactTest` | `variantcentrifuge.association.tests.fisher` | `"fisher"` |
| `LogisticBurdenTest` | `variantcentrifuge.association.tests.logistic_burden` | `"logistic_burden"` |
| `LinearBurdenTest` | `variantcentrifuge.association.tests.linear_burden` | `"linear_burden"` |
| `PurePythonSKATTest` | `variantcentrifuge.association.tests.skat_python` | `"skat"` / `"skat_python"` |
| `RSKATTest` | `variantcentrifuge.association.tests.skat_r` | `"skat"` (when R backend) |
| `PurePythonCOASTTest` | `variantcentrifuge.association.tests.allelic_series_python` | `"coast"` (Python) |
| `COASTTest` | `variantcentrifuge.association.tests.allelic_series` | `"coast"` (R backend) |

ACAT-O is NOT a test class — it is a post-loop meta-test computed by `AssociationEngine`. The functions `compute_acat_o` and `cauchy_combination` in `variantcentrifuge.association.tests.acat` are internal.

## CLI Arguments (HIGH confidence — verified from `variantcentrifuge/cli.py` lines 408-530)

All arguments are in the `stats_group` argument group (`--stats-config` group):

| Flag | Type/Choices | Default | Help |
|------|-------------|---------|------|
| `--perform-association` | `store_true` | False | Enable association analysis |
| `--association-tests` | `str` (comma-sep) | None (→ fisher) | Tests to run: fisher, logistic_burden, linear_burden, skat, skat_python, coast |
| `--skat-backend` | auto/r/python | `python` | SKAT backend |
| `--coast-backend` | auto/r/python | `python` | COAST backend |
| `--covariate-file` | `str` | None | TSV/CSV covariate file (first col = sample ID) |
| `--covariates` | `str` (comma-sep) | None | Subset of covariate columns |
| `--categorical-covariates` | `str` (comma-sep) | None | Force categorical encoding |
| `--trait-type` | binary/quantitative | `binary` | Phenotype scale |
| `--diagnostics-output` | `str` (path) | None | Directory for lambda_GC, QQ data, summary |
| `--variant-weights` | `str` | `beta:1,25` | Weight scheme: beta:a,b, uniform, cadd, revel, combined |
| `--variant-weight-params` | `str` (JSON) | None | Extra weight params, e.g. `'{"cadd_cap": 30}'` |
| `--coast-weights` | `str` (comma-sep) | None | BMV,DMV,PTV weights, e.g. `1,2,3` |
| `--pca-file` | `str` (path) | None | Pre-computed PCA file (PLINK .eigenvec, AKT, or generic TSV) |
| `--pca-tool` | akt | None | PCA tool to invoke |
| `--pca-components` | `int` | 10 | Number of PCs (warn if >20) |

**Validation rules (from `cli.py` lines 1340-1365):**
- `--association-tests` requires `--perform-association`
- `--covariate-file` requires `--perform-association`
- `--diagnostics-output` requires `--perform-association`
- `--pca-file` and `--pca-tool` require `--perform-association`

**Note:** There is NO `--association-config` CLI flag. JSON config is passed via the existing `-c`/`--config` flag as the `"association"` sub-section in the JSON config file.

## JSON Config Format (HIGH confidence — verified from `analysis_stages.py` lines 2044-2073)

The association config lives under an `"association"` key in the main config JSON:

```json
{
  "association": {
    "association_tests": ["fisher", "logistic_burden"],
    "correction_method": "fdr",
    "gene_burden_mode": "samples",
    "trait_type": "binary",
    "covariate_file": "covariates.tsv",
    "covariate_columns": ["age", "sex", "batch"],
    "categorical_covariates": ["sex", "batch"],
    "pca_file": "pca.eigenvec",
    "pca_components": 10,
    "variant_weights": "beta:1,25",
    "variant_weight_params": {"cadd_cap": 40.0},
    "coast_weights": [1, 2, 3],
    "skat_backend": "python",
    "coast_backend": "python",
    "diagnostics_output": "diagnostics/",
    "min_cases": 200,
    "max_case_control_ratio": 20.0,
    "min_case_carriers": 10,
    "confidence_interval_alpha": 0.05,
    "continuity_correction": 0.5,
    "missing_site_threshold": 0.10,
    "missing_sample_threshold": 0.80,
    "firth_max_iter": 25
  }
}
```

**Valid keys (complete list):** All fields in `VALID_ASSOCIATION_KEYS` (from `analysis_stages.py` lines 2044-2073). Unknown keys raise `ValueError` at startup.

**Precedence:** CLI args > JSON `"association"` section > `AssociationConfig` defaults.

## Output Format (HIGH confidence — verified from engine.py and analysis_stages.py)

### Association TSV File

Filename: `{base_name}.association.tsv` (or `.tsv.gz` if compression on) in `--output-dir`.

**Shared columns (always present):**
| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `n_cases` | Total case samples |
| `n_controls` | Total control samples |
| `n_variants` | Qualifying variants for this gene |
| `warnings` | Semicolon-separated warning codes |

**Fisher test columns (when fisher is active):**
| Column | Description |
|--------|-------------|
| `fisher_p_value` | Raw (uncorrected) Fisher p-value |
| `fisher_or` | Odds ratio |
| `fisher_or_ci_lower` | OR 95% CI lower bound |
| `fisher_or_ci_upper` | OR 95% CI upper bound |

**Logistic burden columns (when logistic_burden is active):**
| Column | Description |
|--------|-------------|
| `logistic_burden_p_value` | Raw p-value |
| `logistic_burden_beta` | Beta coefficient |
| `logistic_burden_se` | Standard error |
| `logistic_burden_beta_ci_lower` | Beta 95% CI lower |
| `logistic_burden_beta_ci_upper` | Beta 95% CI upper |

**Linear burden columns (when linear_burden is active):**
| Column | Description |
|--------|-------------|
| `linear_burden_p_value` | Raw p-value |
| `linear_burden_beta` | Beta coefficient |
| `linear_burden_se` | Standard error |
| `linear_burden_beta_ci_lower` | Beta 95% CI lower |
| `linear_burden_beta_ci_upper` | Beta 95% CI upper |

**SKAT columns (when skat is active):**
| Column | Description |
|--------|-------------|
| `skat_p_value` | Raw SKAT p-value |
| `skat_o_rho` | Optimal rho for SKAT-O (in extra) |
| `skat_warnings` | SKAT-specific warning codes |
| `acat_v_p` | ACAT-V per-variant score (in extra, feeds ACAT-O) |

**COAST columns (when coast is active):**
| Column | Description |
|--------|-------------|
| `coast_p_value` | Raw COAST p-value |
| `coast_burden_p_value` | Cauchy of 6 burden components (extra) |
| `coast_skat_p_value` | Allelic SKAT p-value (extra) |
| `coast_n_bmv` | BMV variant count (extra) |
| `coast_n_dmv` | DMV variant count (extra) |
| `coast_n_ptv` | PTV variant count (extra) |

**ACAT-O columns (always present when any test runs):**
| Column | Description |
|--------|-------------|
| `acat_o_p_value` | Raw ACAT-O omnibus p-value |
| `acat_o_corrected_p_value` | FDR/Bonferroni corrected ACAT-O p-value |

**Critical design note:** Individual test p-values (`fisher_p_value`, `skat_p_value`, etc.) are NOT corrected for multiple testing. Only `acat_o_corrected_p_value` is the primary significance measure (ARCH-03 design decision).

### Diagnostics Directory (when `--diagnostics-output` is set)

Three files written by `write_diagnostics()`:

| File | Format | Description |
|------|--------|-------------|
| `lambda_gc.tsv` | TSV: test_name, lambda_gc, n_tests | Genomic inflation factor per test |
| `qq_data.tsv` | TSV: test, expected_neg_log10_p, observed_neg_log10_p | QQ plot data for all tests |
| `summary.txt` | Plain text | Human-readable summary: sample sizes, warnings, lambda_GC values |
| `qq_plot.png` | PNG (optional) | QQ plot — only written if matplotlib is installed |

**lambda_GC interpretation:**
- ~1.0: well-calibrated
- Warn logged if n < 100 (unreliable)
- Values >1.05 suggest inflation (population stratification, cryptic relatedness, or true signal)

## Statistical Tests Reference (HIGH confidence — verified from source code)

### Test Comparison Table

| Test | CLI Name | Trait | Sample Size Need | What It Detects | Computational Cost | When to Use |
|------|---------|-------|-----------------|-----------------|-------------------|-------------|
| Fisher's Exact | `fisher` | Binary | Small OK | Carrier enrichment (collapsed) | Negligible | First-pass, always suitable |
| Logistic Burden | `logistic_burden` | Binary | ≥50-100 cases | Cumulative rare variant burden with covariate adjustment | Low | Moderate cohort, need covariate adjustment |
| Linear Burden | `linear_burden` | Quantitative | ≥50 | Continuous phenotype burden | Low | Quantitative traits |
| SKAT/SKAT-O | `skat` | Binary or Quantitative | ≥100 cases recommended | Heterogeneous effects (not all variants causal) | Moderate | Large cohorts, uncertain effect direction |
| COAST | `coast` | Binary or Quantitative | ≥100 cases recommended | Ordered allelic series (PTV > DMV > BMV) | Moderate | When monotone effect size trend expected |
| ACAT-O | (automatic) | Both | Depends on components | Omnibus combining all tests | Negligible (post-hoc) | Always computed, primary significance measure |

### Fisher's Exact Test
- **Model:** 2×2 contingency table (carrier vs non-carrier in cases vs controls)
- **When appropriate:** Any cohort size; no covariate adjustment
- **Assumptions:** Random sampling; independence of carriers
- **Output:** Odds ratio with 95% CI, raw p-value
- **Collapsing modes:** `samples` (unique carrier count) or `alleles` (max dosage)

### Logistic Burden Test
- **Model:** Weighted burden score → Logistic regression with Firth fallback for separation
- **When appropriate:** Binary trait, moderate-to-large cohort (≥50 cases), covariate adjustment needed
- **Assumptions:** Logistic model; variants have same direction of effect
- **Firth fallback:** Auto-triggered when convergence fails or BSE >100 (separation detection)
- **Warning codes:** `PERFECT_SEPARATION`, `QUASI_SEPARATION`, `ZERO_CARRIERS_ONE_GROUP`, `LOW_CARRIER_COUNT`, `FIRTH_CONVERGE_FAIL`, `NO_GENOTYPE_MATRIX`

### Linear Burden Test
- **Model:** Weighted burden score → OLS regression
- **When appropriate:** Quantitative phenotype (biomarker, measurement), covariate adjustment needed
- **Assumptions:** OLS; continuous, approximately normal phenotype
- **Output:** Beta coefficient + SE + 95% CI

### SKAT / SKAT-O
- **Model:** Score-based quadratic kernel test (SKAT); optimal-weight omnibus combining SKAT and burden (SKAT-O)
- **When appropriate:** Large cohorts; heterogeneous effects (some variants protective, some harmful)
- **Python backend:** numpy/scipy (default, thread-safe)
- **R backend:** rpy2 (deprecated, not thread-safe, requires R + SKAT package)
- **Also computes:** ACAT-V per-variant score test (stored in extra, feeds ACAT-O)
- **Output:** p-value only (no effect size for SKAT)

### COAST (Combined Omnibus Association Test)
- **Model:** Allelic series test with BMV/DMV/PTV categories; Cauchy combination of 6 burden components + allelic SKAT
- **Reference:** McCaw et al. AJHG 2023
- **Variant classification:** PTV (HIGH impact + stop_gained/frameshift/splice), DMV (missense + SIFT D/PolyPhen P/D), BMV (missense + SIFT T + PolyPhen B)
- **Requires:** SIFT and PolyPhen predictions in VCF annotations (columns tried: dbNSFP_SIFT_pred, SIFT_pred, dbNSFP_Polyphen2_HDIV_pred, PolyPhen_pred)
- **When appropriate:** When monotone effect size trend expected (PTVs most damaging, BMVs least)
- **Output:** p-value + per-category variant counts

### ACAT-O (Aggregated Cauchy Association Test — Omnibus)
- **Model:** Cauchy combination of all primary test p-values per gene (including ACAT-V when SKAT ran)
- **Reference:** Liu et al. 2019 AJHG; Liu and Xie 2020 JASA
- **FDR correction:** Applied ONLY to ACAT-O p-values across all genes (single correction pass, ARCH-03)
- **When appropriate:** Always use `acat_o_corrected_p_value` as primary significance measure

## Variant Weights (HIGH confidence — verified from `weights.py`)

| Scheme | CLI Spec | Description |
|--------|---------|-------------|
| Beta(MAF) | `beta:1,25` | Default; rare variants upweighted; matches SKAT R package convention |
| Uniform | `uniform` | All variants equal weight |
| CADD | `cadd` | Beta(MAF) × min(CADD_phred/40, 1.0); requires CADD annotations |
| REVEL | `revel` | Beta(MAF) × REVEL_score [0,1]; requires REVEL annotations |
| Combined | `combined` | Beta(MAF) × functional score; uses CADD, falls back to REVEL |

**Custom params via `--variant-weight-params`:** JSON string, e.g. `'{"cadd_cap": 30}'` to cap CADD at 30 instead of 40.

**Fallback for missing scores:** LoF variants → functional weight 1.0 (conservative). Missense without predictions → 1.0 (no up/down weighting). Warning logged with per-category missing counts.

## PCA Integration (HIGH confidence — verified from `pca.py`)

Supported PCA file formats:
1. **PLINK .eigenvec with header** (line starts with `#FID` or `FID`)
2. **PLINK .eigenvec without header** (two non-numeric ID cols then numeric)
3. **AKT stdout / generic TSV** (one sample ID col + numeric columns)

When `--pca-tool akt` is set, AKT is invoked as a subprocess (AKT must be in PATH).

PCA columns are appended to the covariate matrix. Use `--pca-components N` to limit PCs (warn if N > 20).

## Covariate System (HIGH confidence — verified from `covariates.py`)

Covariate file requirements:
- First column: sample IDs matching VCF sample IDs
- First row: column headers
- Delimiter: auto-detected (.tsv/.tab → tab, .csv → comma, otherwise sniff)

Auto-detection for categorical columns: non-numeric columns with ≤5 unique values are one-hot encoded. Override with `--categorical-covariates`.

## Existing Docs to Update (HIGH confidence — verified by reading files)

### `docs/source/usage.md`
- Current state: Has "Statistical Analysis" table (lines ~189-198) with `--perform-gene-burden`, `--correction-method`, etc.
- Missing: All `--perform-association` and related flags
- Action: Add "Association Testing" subsection after "Statistical Analysis" with table of association flags + brief note + link to guide

### `docs/source/guides/cohort_analysis.md`
- Current state: Covers cohort reports, gene burden testing (`--perform-gene-burden`), Snakemake
- Missing: Association testing mention
- Action: Add brief section under "Advanced Cohort Analysis" pointing to association_testing.md guide

### `docs/source/faq.md`
- Current state: General questions, installation, basic filtering
- Missing: No association testing content
- Action: Add 2-3 FAQ entries about association testing (e.g., "How do I run a gene burden test with covariates?", "What is ACAT-O?")

### `docs/source/index.md`
- Current state: Mentions "Gene Burden Testing: Statistical analysis for case-control studies using Fisher's exact test" in Comprehensive Analysis Tools
- Missing: New association framework
- Action: Update bullet to mention full framework; add `association_testing` to the toctree under "Practical Guides"

### `README.md`
- Current state: Feature list includes "Gene burden analysis with Fisher's exact test"
- Action: Add feature bullet: "Modular association testing framework: SKAT-O, COAST, burden tests, ACAT-O omnibus, covariate/PCA support"

### `docs/source/api/index.md`
- Current state: No `association` submodule listed anywhere
- Action: Add "Association Testing" section in toctree with `api/association` entry

### `docs/source/changelog.md`
- Current state: Uses Keep a Changelog format; latest is `[0.14.0] - 2026-02-18`
- Action: Add `[0.15.0] - 2026-02-22` section under `## [Unreleased]` or as new section

## Changelog Entry Format (HIGH confidence — verified from existing changelog)

Format follows Keep a Changelog:

```markdown
## [0.15.0] - 2026-02-22

### Added
- **Modular association testing framework** (v0.15.0) across 8 phases:
  - **Core abstractions**: AssociationTest ABC, TestResult, AssociationConfig, AssociationEngine orchestrator
  - **Simple tests**: Fisher's exact test (carrier/allele modes, OR + 95% CI)
  - **Burden regression**: Logistic burden (binary, Firth fallback for separation) and linear burden (quantitative, OLS)
  - **Kernel tests**: SKAT and SKAT-O with pure Python backend (numpy/scipy); R backend available but deprecated
  - **Allelic series**: COAST test (BMV/DMV/PTV categories, ordered effect hypothesis); pure Python backend default
  - **Omnibus combining**: ACAT-O per-gene omnibus via Cauchy combination; ACAT-V per-variant score within SKAT
  - **Covariate system**: TSV/CSV covariate files, auto-detected delimiter, one-hot encoding for categorical variables
  - **PCA integration**: PLINK .eigenvec, AKT output, and generic TSV formats; optional AKT subprocess invocation
  - **Variant weights**: Beta(MAF), uniform, CADD, REVEL, and combined functional weight schemes
  - **Diagnostics**: lambda_GC, QQ data TSV, optional matplotlib QQ plot PNG
  - **JSON config**: `"association"` section in config.json with full validation and CLI override precedence
  - **FDR strategy**: Single Benjamini-Hochberg correction pass on ACAT-O p-values only (ARCH-03)
```

## Architecture Patterns for the Guide

### Section Structure (from CONTEXT.md decisions)

```
1. Quick Start (5-10 lines, copy-paste-ready)
2. Setup
   - Covariates
   - PCA
3. Simple Tests
   - Fisher's Exact Test
   - Burden Tests (Logistic and Linear)
4. Advanced Tests
   - SKAT / SKAT-O
   - COAST
5. Combining Tests: ACAT-O
6. Tuning
   - Variant Weights
   - Diagnostics (lambda_GC, QQ plot)
   - JSON Config
7. Test Selection Reference (comparison table)
8. Troubleshooting
```

### Quick Start Example (from CLI flags)

```bash
# Minimal: Fisher's exact test
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --output-file results.tsv
```

Output: `results.association.tsv` with `gene`, `n_cases`, `n_controls`, `n_variants`, `fisher_p_value`, `fisher_or`, `fisher_or_ci_lower`, `fisher_or_ci_upper`, `acat_o_p_value`, `acat_o_corrected_p_value`, `warnings`.

### Comprehensive Example

```bash
# SKAT-O + logistic burden + COAST with covariates and diagnostics
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests skat,logistic_burden,coast \
  --trait-type binary \
  --covariate-file covariates.tsv \
  --covariates age,sex,batch \
  --categorical-covariates sex,batch \
  --pca-file ancestry.eigenvec \
  --pca-components 10 \
  --variant-weights beta:1,25 \
  --diagnostics-output diagnostics/ \
  --output-file results.tsv
```

## Common Pitfalls

### Pitfall 1: COAST Fails with "NO_CLASSIFIABLE_VARIANTS"
**What goes wrong:** COAST requires SIFT and PolyPhen predictions to classify BMV/DMV. If the VCF lacks these annotations, all missense variants get code 0 (unclassified) and COAST is skipped.
**How to avoid:** Annotate with dbNSFP (provides `dbNSFP_SIFT_pred` and `dbNSFP_Polyphen2_HDIV_pred`). COAST also requires at least one variant in each category (BMV, DMV, PTV).
**Warning signs:** `coast_skip_reason: NO_CLASSIFIABLE_VARIANTS` in output.

### Pitfall 2: Separation in Logistic Burden
**What goes wrong:** When all carriers are cases (or controls), logistic regression cannot converge. Firth fallback is triggered automatically but may also fail.
**Warning codes:** `PERFECT_SEPARATION`, `QUASI_SEPARATION`, `FIRTH_CONVERGE_FAIL`
**How to avoid:** Ensure adequate controls; use Fisher's exact test for small cohorts.

### Pitfall 3: Primary Test P-values are Not Corrected
**What goes wrong:** Users report `fisher_p_value` as significant without understanding it is uncorrected.
**Correct approach:** Use `acat_o_corrected_p_value` as the primary significance measure. Individual test p-values are for diagnostic signal decomposition only.

### Pitfall 4: Sample ID Mismatch Between VCF and Covariate File
**What goes wrong:** `ValueError` if covariate file sample IDs don't match VCF sample IDs (case-sensitive).
**How to avoid:** Run `bcftools query -l input.vcf.gz` to get exact sample IDs, use them in the covariate file first column.

### Pitfall 5: lambda_GC Unreliable on Small Panels
**What goes wrong:** lambda_GC is flagged `[UNRELIABLE: n < 100]` in `diagnostics/lambda_gc.tsv`.
**How to avoid:** Use lambda_GC interpretation only with ≥100 tested genes. For small panels, rely on individual gene results.

### Pitfall 6: R Backend Not Found
**What goes wrong:** `ImportError: rpy2 not installed` when `--skat-backend r` is specified.
**How to avoid:** Use the default Python backend (`--skat-backend python`). R backend requires `pip install rpy2` and R with the SKAT package.

### Pitfall 7: ACAT-O Column Absent from Output
**What goes wrong:** User runs only `fisher` and doesn't see ACAT-O columns.
**Clarification:** ACAT-O columns are always present — when only one test runs, `acat_o_p_value` equals the single test's p-value (pass-through behavior per Cauchy combination spec).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Multiple testing correction | Custom BH implementation | `apply_correction()` from `variantcentrifuge.association` | Already handles both FDR and Bonferroni; tested |
| P-value combination | Custom Cauchy formula | `cauchy_combination()` from acat.py (internal) | Handles numerical stability for p < 1e-16 |
| Covariate loading | Manual CSV parsing | `load_covariates()` from covariates.py | Handles delimiter detection, alignment, encoding |
| PCA file loading | Manual format detection | `load_pca_file()` from pca.py | Handles 3 PLINK/AKT formats |
| lambda_GC computation | Manual chi2 conversion | `compute_lambda_gc()` from diagnostics.py | Already handles None/NaN filtering |

## Existing API Stub Patterns to Follow

The existing stub pattern from `docs/source/api/gene_burden.md`:

```markdown
# [Module Name]

[One-line description]

```{eval-rst}
.. automodule:: variantcentrifuge.[module_path]
   :members:
   :undoc-members:
   :show-inheritance:
```
```

For the association module, one stub file covering the entire subpackage, or individual files per submodule. Given the size of the association package, individual stubs for the major submodules are cleaner:

Recommended stubs:
- `docs/source/api/association.md` — the main `__init__.py` exports (AssociationConfig, AssociationEngine, AssociationTest, TestResult, apply_correction)
- Could expand to submodule stubs if desired, but the CONTEXT.md says "public classes only"

## Sources

### Primary (HIGH confidence)
- `variantcentrifuge/association/__init__.py` — public exports verified
- `variantcentrifuge/association/base.py` — AssociationConfig, AssociationTest, TestResult verified
- `variantcentrifuge/association/engine.py` — AssociationEngine, output columns verified
- `variantcentrifuge/association/diagnostics.py` — diagnostics output files verified
- `variantcentrifuge/association/weights.py` — weight schemes verified
- `variantcentrifuge/association/pca.py` — PCA formats verified
- `variantcentrifuge/association/covariates.py` — covariate loading verified
- `variantcentrifuge/association/tests/fisher.py` — FisherExactTest verified
- `variantcentrifuge/association/tests/logistic_burden.py` — LogisticBurdenTest verified
- `variantcentrifuge/association/tests/linear_burden.py` — LinearBurdenTest verified
- `variantcentrifuge/association/tests/skat_python.py` — PurePythonSKATTest verified
- `variantcentrifuge/association/tests/acat.py` — ACAT-V, ACAT-O, cauchy_combination verified
- `variantcentrifuge/association/tests/allelic_series.py` — COASTTest, classify_variants verified
- `variantcentrifuge/association/tests/allelic_series_python.py` — PurePythonCOASTTest verified
- `variantcentrifuge/cli.py` lines 408-530 — all CLI argument definitions verified
- `variantcentrifuge/stages/analysis_stages.py` lines 2044-2291 — JSON config keys, VALID_ASSOCIATION_KEYS, _build_assoc_config_from_context verified
- `docs/source/conf.py` — Sphinx extensions and format verified
- `docs/source/index.md` — toctree structure verified
- `docs/source/guides/index.md` — guides toctree verified
- `docs/source/api/index.md` — API index structure verified
- `docs/source/api/gene_burden.md` — API stub pattern verified
- `docs/source/changelog.md` — changelog format and version history verified
- `docs/source/usage.md` — existing CLI flag documentation structure verified
- `docs/source/guides/cohort_analysis.md` — existing cohort guide content verified
- `docs/source/faq.md` — FAQ content and gaps verified
- `README.md` — existing feature list and structure verified

## Metadata

**Confidence breakdown:**
- Documentation infrastructure: HIGH — all paths and formats verified from source
- CLI arguments: HIGH — verified from `cli.py` argparse definitions
- Association module API: HIGH — verified from `__init__.py`, `base.py`, `engine.py`
- Output columns: HIGH — verified from `engine.py` docstring and `run_all()` implementation
- Test statistics descriptions: HIGH — verified from docstrings in test files
- JSON config format: HIGH — verified from `VALID_ASSOCIATION_KEYS` and `_build_assoc_config_from_context()`
- Diagnostics output: HIGH — verified from `write_diagnostics()`

**Research date:** 2026-02-22
**Valid until:** 2026-03-22 (stable — documentation-only phase, code not changing)
