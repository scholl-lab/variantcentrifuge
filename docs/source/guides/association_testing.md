# Association Testing Guide

This guide covers all association testing functionality in VariantCentrifuge v0.16.0, including
Fisher's exact test, logistic and linear burden tests, SKAT/SKAT-O, COAST (allelic series),
the ACAT-O omnibus, and gene-level FDR weighting. It assumes VariantCentrifuge is already installed and you have a multi-sample
VCF annotated with functional predictions and population allele frequencies.

---

## Quick Start

Run a Fisher's exact test on all genes in your gene file:

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --output-file results.tsv
```

Output: `results.association.tsv` with columns:

```
gene  n_cases  n_controls  n_variants  fisher_pvalue  fisher_or  fisher_or_ci_lower  fisher_or_ci_upper  acat_o_pvalue  acat_o_qvalue  warnings
GENE1 312      488         7           0.00031         3.82       1.87                7.78                0.00031         0.0047
GENE2 312      488         3           0.41            1.23       0.74                2.05                0.41            1.0
...
```

**Primary significance measure:** `acat_o_qvalue` — this is the FDR-corrected omnibus
p-value. Use it for all significance decisions.

---

## Setup: Covariates

### Covariate File Format

Covariate files are TSV or CSV with sample IDs in the first column and a header row:

```
sample_id   age   sex   batch   gfr
SAMPLE_001  45    F     batch1  62.3
SAMPLE_002  51    M     batch1  88.1
SAMPLE_003  38    F     batch2  45.7
...
```

- First column: sample IDs matching VCF sample IDs exactly (case-sensitive)
- Delimiter: auto-detected (`.tsv`/`.tab` → tab, `.csv` → comma, otherwise sniffed)
- Header row required

### CLI Flags

```bash
--covariate-file covariates.tsv     # Path to covariate file
--covariates age,sex,batch          # Use only these columns (default: all)
--categorical-covariates sex,batch  # Force one-hot encoding for these columns
```

**Auto-detection of categorical columns:** Non-numeric columns with 5 or fewer unique values are
automatically one-hot encoded. Override with `--categorical-covariates` if needed (e.g., a batch
column with many unique integer codes that is still categorical).

### Example with Covariates

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests logistic_burden \
  --covariate-file covariates.tsv \
  --covariates age,sex,batch \
  --categorical-covariates sex,batch \
  --output-file results.tsv
```

:::{note}
`--covariate-file` requires `--perform-association`. Specifying it without enabling association
analysis raises a validation error at startup.
:::

---

## Setup: PCA

Population stratification is a major confounder in case-control association studies. Include
principal components as covariates to correct for it.

### Supported PCA File Formats

VariantCentrifuge auto-detects three PCA file formats:

**PLINK .eigenvec with header** (line starts with `#FID` or `FID`):
```
#FID  IID         PC1       PC2       PC3
FAM1  SAMPLE_001  0.0142    -0.0231   0.0089
FAM1  SAMPLE_002  -0.0305   0.0118    -0.0042
```

**PLINK .eigenvec without header** (two non-numeric ID columns then numeric):
```
FAM1  SAMPLE_001  0.0142    -0.0231   0.0089
FAM1  SAMPLE_002  -0.0305   0.0118    -0.0042
```

**AKT stdout / generic TSV** (one sample ID column then numeric):
```
sample_id   PC1       PC2       PC3
SAMPLE_001  0.0142    -0.0231   0.0089
SAMPLE_002  -0.0305   0.0118    -0.0042
```

PLINK format uses IID (column 2) as the sample identifier — this matches VCF sample names.

### CLI Flags

```bash
--pca-file ancestry.eigenvec    # Pre-computed PCA file
--pca-components 10             # Number of PCs to include (default: 10, warn if >20)
--pca-tool akt                  # Invoke AKT as subprocess to compute PCA on-the-fly
```

PCA columns are automatically appended to the covariate matrix before any regression. If you
specify both `--covariate-file` and `--pca-file`, PCs are added as additional columns.

:::{warning}
Using `--pca-tool akt` requires AKT to be in your PATH. A missing AKT binary raises a hard error
(not a silent skip) to prevent undetected misconfiguration.
:::

:::{tip}
Using more than 20 PCs (`--pca-components >20`) triggers a warning. Over-conditioning on PCs can
absorb true signal from rare variant tests.
:::

---

## Simple Tests: Fisher's Exact Test

Fisher's exact test uses a 2×2 contingency table (carrier vs non-carrier in cases vs controls)
to detect carrier enrichment. It requires no distributional assumptions and works at any sample
size, making it an appropriate first-pass screen for all cohorts. It does not adjust for
covariates.

**When to use:** Any cohort size; first-pass screening; when covariate adjustment is not required.

**CLI name:** `fisher` (default when `--perform-association` is specified without
`--association-tests`)

### CLI Example

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests fisher \
  --output-file results.tsv
```

### Output Columns

| Column | Description |
|--------|-------------|
| `fisher_pvalue` | Raw (uncorrected) Fisher p-value |
| `fisher_or` | Odds ratio |
| `fisher_or_ci_lower` | OR 95% CI lower bound |
| `fisher_or_ci_upper` | OR 95% CI upper bound |

### Collapsing Modes

The `gene_burden_mode` setting controls how variants are collapsed per sample:

- `samples` (default): Binary carrier indicator — a sample is a carrier if it has ≥1 qualifying
  variant allele. This is the CMC/CAST approach.
- `alleles`: Maximum dosage across variants per sample (0, 1, or 2). More sensitive to dosage
  effects.

Set via JSON config (see [JSON Config](#tuning-json-config)):
```json
{ "association": { "gene_burden_mode": "alleles" } }
```

---

## Simple Tests: Burden Tests

Burden tests collapse all qualifying variants in a gene into a single weighted burden score per
sample, then regress phenotype on that score. This approach maximises power when variants in a
gene all act in the same direction.

### Logistic Burden Test (Binary Traits)

The logistic burden test fits a logistic regression model with the weighted burden score (plus
optional covariates) as predictors and case/control status as outcome. It automatically falls
back to Firth penalised likelihood when standard logistic regression fails to converge or when
quasi-complete separation is detected (BSE > 100).

**When to use:** Binary trait (case/control), moderate-to-large cohort (≥50 cases), covariate
adjustment required.

**CLI name:** `logistic_burden`

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests logistic_burden \
  --trait-type binary \
  --covariate-file covariates.tsv \
  --covariates age,sex,batch \
  --categorical-covariates sex,batch \
  --output-file results.tsv
```

**Output columns:**

| Column | Description |
|--------|-------------|
| `logistic_burden_pvalue` | Raw p-value |
| `logistic_burden_beta` | Beta coefficient (log-odds per unit burden) |
| `logistic_burden_se` | Standard error of beta |
| `logistic_burden_beta_ci_lower` | Beta 95% CI lower bound |
| `logistic_burden_beta_ci_upper` | Beta 95% CI upper bound |

:::{note}
Burden tests report **beta + SE**, not odds ratios. The per-unit burden OR would be numerically
close to 1.0 due to beta-distributed variant weights and is not interpretable. Beta matches the
SKAT/SAIGE-GENE convention.
:::

**Separation warning codes:** `PERFECT_SEPARATION`, `QUASI_SEPARATION`, `FIRTH_CONVERGE_FAIL`,
`ZERO_CARRIERS_ONE_GROUP`, `LOW_CARRIER_COUNT`. These appear in the `warnings` column.

### Linear Burden Test (Quantitative Traits)

The linear burden test fits an OLS regression with the weighted burden score (plus optional
covariates) as predictors and a continuous phenotype as outcome.

**When to use:** Quantitative trait (biomarker, laboratory measurement), covariate adjustment
required.

**CLI name:** `linear_burden`

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column gfr \
  --perform-association \
  --association-tests linear_burden \
  --trait-type quantitative \
  --covariate-file covariates.tsv \
  --covariates age,sex,batch \
  --output-file results.tsv
```

**Output columns:**

| Column | Description |
|--------|-------------|
| `linear_burden_pvalue` | Raw p-value |
| `linear_burden_beta` | Beta coefficient |
| `linear_burden_se` | Standard error |
| `linear_burden_beta_ci_lower` | Beta 95% CI lower bound |
| `linear_burden_beta_ci_upper` | Beta 95% CI upper bound |

---

## Advanced Tests: SKAT / SKAT-O

The Sequence Kernel Association Test (SKAT) is a score-based quadratic kernel test that does not
assume all variants in a gene have the same direction of effect. It is most powerful when effects
are heterogeneous — some variants protective, others harmful, or only a subset truly causal.
SKAT-O (optimal SKAT) mixes SKAT and burden test statistics to maximise power over a range of
effect architectures, estimating the optimal mixing parameter rho.

**When to use:** Large cohorts (≥100 cases recommended); when effect direction is uncertain or
heterogeneous.

**CLI name:** `skat` (activates SKAT-O by default)

### Backends

| Backend | CLI Flag | Notes |
|---------|----------|-------|
| Pure Python (default) | `--skat-backend python` | numpy/scipy; thread-safe; no R required |
| Auto | `--skat-backend auto` | Selects Python (same as `python` in v0.15.0) |
| R (deprecated) | `--skat-backend r` | rpy2 + R SKAT package; not thread-safe; raises DeprecationWarning |

The Python backend uses the compiled Davies method (qfc.c via ctypes) for exact p-values with a
saddlepoint approximation fallback for extreme values, and Liu fallback as a last resort. Results
are validated to match the R SKAT package within 10%.

### CLI Example

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests skat \
  --trait-type binary \
  --covariate-file covariates.tsv \
  --pca-file ancestry.eigenvec \
  --output-file results.tsv
```

### Output Columns

| Column | Description |
|--------|-------------|
| `skat_o_pvalue` | Raw SKAT-O p-value |
| `skat_o_rho` | Optimal rho (mixing parameter between SKAT and burden) |
| `acat_v_p` | ACAT-V per-variant score (internal; feeds ACAT-O) |

SKAT produces no effect size estimate. The p-value is the primary output.

### ACAT-V

When SKAT runs, ACAT-V (Aggregated Cauchy Association Test — Variant) is automatically computed
as a per-variant score combination. It is not a separate test to activate; its p-value is stored
in `acat_v_p` and fed into the ACAT-O omnibus alongside the SKAT p-value.

---

## Advanced Tests: COAST

COAST (Combined Omnibus Association Test for allelic Series) tests the hypothesis that variants
with greater predicted functional impact have larger effect sizes, following the ordered
alternative PTV > DMV > BMV. It combines six burden components (BMV sum/max, DMV sum/max, PTV
sum/max) and an allelic SKAT via Cauchy combination. This test is appropriate when you have a
biological prior that PTVs are most damaging, with missense variants spanning a range of
functional severity.

**Reference:** McCaw et al. AJHG 2023.

**When to use:** When a monotone effect size trend (PTVs most damaging) is expected; requires
SIFT and PolyPhen predictions in VCF annotations.

**CLI name:** `coast`

### Variant Classification

COAST classifies qualifying variants into three categories:

| Category | Code | Criteria |
|----------|------|---------|
| **PTV** (Protein-Truncating) | 3 | HIGH impact AND effect in {stop_gained, frameshift_variant, splice_acceptor_variant, splice_donor_variant} |
| **DMV** (Deleterious Missense) | 2 | Missense + (SIFT: D OR PolyPhen: P/D) — either prediction suffices |
| **BMV** (Benign Missense) | 1 | Missense + (SIFT: T) + (PolyPhen: B) — both required |
| Unclassified | 0 | Missense without SIFT/PolyPhen predictions — excluded from COAST only |

Annotation columns tried (first found wins):
- SIFT: `dbNSFP_SIFT_pred`, `SIFT_pred`, `sift_pred`
- PolyPhen: `dbNSFP_Polyphen2_HDIV_pred`, `dbNSFP_Polyphen2_HVAR_pred`, `Polyphen2_HDIV_pred`, `Polyphen2_HVAR_pred`, `polyphen2_hdiv_pred`

Unclassified variants (code 0) are excluded from COAST but remain in SKAT and burden tests.

### Backends

| Backend | CLI Flag | Notes |
|---------|----------|-------|
| Pure Python (default) | `--coast-backend python` | No R required; thread-safe; matches R AllelicSeries |
| Auto | `--coast-backend auto` | Probes rpy2 AND AllelicSeries R package; falls back to Python |
| R (deprecated) | `--coast-backend r` | Requires rpy2 + R AllelicSeries package; raises DeprecationWarning |

### CLI Example

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests coast \
  --trait-type binary \
  --covariate-file covariates.tsv \
  --coast-weights 1,2,3 \
  --output-file results.tsv
```

### COAST Weights

`--coast-weights BMV,DMV,PTV` sets the ordered alternative weights. Default is `1,2,3`
(PTVs three times the weight of BMVs). Adjust based on biological expectation.

### Output Columns

| Column | Description |
|--------|-------------|
| `coast_pvalue` | Raw COAST omnibus p-value |
| `coast_burden_pvalue` | Cauchy of the 6 burden components |
| `coast_skat_pvalue` | Allelic SKAT p-value |
| `coast_n_bmv` | Number of BMV variants |
| `coast_n_dmv` | Number of DMV variants |
| `coast_n_ptv` | Number of PTV variants |

---

## Combining Tests: ACAT-O

ACAT-O (Aggregated Cauchy Association Test — Omnibus) combines all primary test p-values for a
gene into a single omnibus p-value using the Cauchy combination. This captures signal regardless
of which test is most powerful for a given gene. The Cauchy combination is particularly robust to
correlation between test statistics.

**References:** Liu et al. 2019 AJHG; Liu and Xie 2020 JASA.

ACAT-O is **not** a test you activate separately — it is computed automatically after all
selected tests have run. When SKAT is active, ACAT-V is also included in the combination.

### FDR Correction Strategy

A **single** Benjamini-Hochberg FDR correction pass is applied to the ACAT-O p-values across all
genes. Individual test p-values (`fisher_pvalue`, `skat_o_pvalue`, etc.) are **not** corrected.

:::{warning}
Do not apply additional correction to individual test p-values, and do not apply FDR separately
per test. Both approaches are statistically incorrect given this design. Use
`acat_o_qvalue` as your primary significance measure.
:::

### Output Columns

| Column | Description |
|--------|-------------|
| `acat_o_pvalue` | Raw (uncorrected) ACAT-O omnibus p-value |
| `acat_o_qvalue` | FDR-corrected (Benjamini-Hochberg) ACAT-O p-value |

### Pass-Through Behaviour

When only one test is active, `acat_o_pvalue` equals that test's raw p-value (pass-through).
FDR correction is still applied across genes. This is by design — the omnibus is always a safe
primary measure.

---

## Tuning: Variant Weights

Variant weights control the contribution of each qualifying variant to the burden score in
logistic/linear burden tests and to the SKAT kernel matrix. Well-chosen weights improve power
by upweighting likely pathogenic variants.

### Weight Schemes

| Scheme | CLI Spec | Description |
|--------|---------|-------------|
| Beta(MAF) | `beta:1,25` | Default. Beta(MAF; a=1, b=25) — rare variants upweighted; matches R SKAT convention |
| Beta custom | `beta:a,b` | Custom Beta distribution parameters, e.g. `beta:1,1` for flat weighting |
| Uniform | `uniform` | All qualifying variants receive weight 1.0 |
| CADD | `cadd` | Beta(MAF) × min(CADD_phred / 40, 1.0); requires `CADD_phred` or `dbNSFP_CADD_phred` |
| REVEL | `revel` | Beta(MAF) × REVEL_score [0,1]; requires `REVEL_score` or `dbNSFP_REVEL_score` |
| Combined | `combined` | Beta(MAF) × functional score; prefers CADD, falls back to REVEL |

### CLI Flags

```bash
--variant-weights beta:1,25              # Weight scheme (default: beta:1,25)
--variant-weight-params '{"cadd_cap":30}'  # JSON string with extra parameters
```

`--variant-weight-params` accepts a JSON string for scheme-specific tuning. Current supported
parameter: `cadd_cap` (float, default 40.0) — cap CADD phred score before normalisation.

### Fallback Behaviour for Missing Annotations

When functional annotations are absent:
- LoF (HIGH impact) variants: functional weight 1.0 (conservative, no down-weighting)
- Missense without SIFT/PolyPhen/CADD/REVEL: functional weight 1.0
- A warning is logged with per-category missing annotation counts

:::{tip}
Start with the default `beta:1,25` weights. Switch to `combined` or `cadd` only when CADD
annotations are present in your VCF and you want to differentiate variant functional severity
beyond MAF alone.
:::

---

## Tuning: Diagnostics

Enable diagnostics to assess test calibration:

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests fisher,skat,logistic_burden \
  --covariate-file covariates.tsv \
  --diagnostics-output diagnostics/ \
  --output-file results.tsv
```

### Output Files

| File | Format | Description |
|------|--------|-------------|
| `lambda_gc.tsv` | TSV: `test_name`, `lambda_gc`, `n_tests` | Genomic inflation factor per test |
| `qq_data.tsv` | TSV: `test`, `expected_neg_log10_p`, `observed_neg_log10_p` | QQ plot data (sorted ascending by expected) |
| `summary.txt` | Plain text | Human-readable summary: sample sizes, per-gene warnings, lambda_GC values |
| `qq_plot.png` | PNG (optional) | QQ plot — written only if `matplotlib` is installed |

### Interpreting lambda_GC

The genomic inflation factor lambda_GC is the ratio of the median observed chi-squared statistic
to the expected median (0.455 for chi2 with 1 df):

| Value | Interpretation |
|-------|---------------|
| ~1.0 | Well-calibrated — test statistics match the null distribution |
| >1.05 | Inflation — investigate population stratification, cryptic relatedness, or true polygenic signal |
| <0.95 | Deflation — may indicate overly conservative test or sparse variants |
| Flagged `[UNRELIABLE: n < 100]` | Fewer than 100 tested genes — lambda_GC is unreliable at this scale |

:::{note}
lambda_GC is unreliable when fewer than 100 genes were tested. For small gene panels, skip
lambda_GC and rely on individual gene results directly.
:::

:::{note}
For Fisher's exact test, lambda_GC is computed as a diagnostic indicator only — it is **not**
used to correct or adjust Fisher p-values. Use it to detect gross calibration issues, not as a
correction factor.
:::

---

## Tuning: Gene-Level FDR Weighting

By default, Benjamini-Hochberg FDR correction treats all genes equally. When you have prior
biological knowledge about which genes are more likely to be associated with your phenotype, you
can provide per-gene weights to increase statistical power for those genes while maintaining
overall FDR control.

This implements the **weighted Benjamini-Hochberg** procedure (Genovese et al. 2006), which
divides each gene's p-value by its weight before BH ranking. Genes with higher weights are
effectively given a larger share of the FDR budget, making them easier to detect. The overall
FDR guarantee is preserved as long as weights average to 1.0, which VariantCentrifuge enforces
automatically through renormalization.

### How It Works

1. You provide a TSV file mapping gene names to positive weights (higher = more likely relevant)
2. At correction time, weights are renormalized so their mean across *tested* genes equals 1.0
3. Each gene's ACAT-O p-value is divided by its normalized weight before BH ranking
4. Genes absent from the weight file receive weight 1.0 (neutral)

The net effect: high-weight genes become easier to call significant, low-weight genes become
harder, and unweighted genes are slightly penalized to compensate. If your priors are correct,
you gain power. If wrong, you lose some — but FDR is still controlled at the stated level.

### Weight File Format

A two-column TSV with a header row. The first column is the gene name (must match gene names in
the association output exactly), and the second column is the weight (positive number):

```
gene	weight
PKD1	3.0
PKD2	2.5
IFT140	2.0
BRCA1	1.8
BRCA2	1.8
TP53	1.5
TTN	0.5
```

Weights do not need to sum or average to any particular value — renormalization is automatic. A
weight of 2.0 means "I believe this gene is twice as likely to matter as a typical gene." A
weight below 1.0 means "this gene is less likely to matter" (e.g., TTN, which produces many
rare variants due to its size but is often a nuisance signal).

### Deriving Weights

There is no single correct way to derive weights. The choice depends on your study and what
prior information is available. Here are practical approaches:

**Gene constraint scores (recommended starting point):**
Use gnomAD pLI or LOEUF scores, which quantify how intolerant a gene is to loss-of-function
variation. Highly constrained genes are more likely to be disease-relevant:

```python
# Example: convert pLI to weights
# pLI ranges from 0 to 1; higher = more constrained
import pandas as pd
df = pd.read_csv("gnomad_constraint.tsv", sep="\t")
df["weight"] = 1.0 + 2.0 * df["pLI"]  # range: 1.0 to 3.0
df[["gene", "weight"]].to_csv("gene_weights.tsv", sep="\t", index=False)
```

**GWAS proximity:**
Genes near genome-wide significant GWAS loci for your phenotype get higher weights:

```
gene	weight
PKD1	3.0    # GWAS hit for kidney disease
PKD2	2.5    # GWAS hit for kidney disease
UMOD	2.0    # GWAS hit for kidney function
APOL1	2.0    # Known risk gene in specific populations
```

**Clinical gene panels:**
Genes on an established diagnostic panel for your phenotype get higher weights, non-panel genes
get neutral or lower weights:

```python
panel_genes = set(open("panel_genes.txt").read().split())
# Panel genes: weight 2.0, others: weight 1.0
weights = {g: 2.0 if g in panel_genes else 1.0 for g in all_genes}
```

**Literature-based:**
Manually assign weights based on published evidence. This is subjective but often the most
practical approach for small targeted studies.

**Combining sources:**
Multiply independent evidence sources (each normalized to mean ~1.0):

```python
weight = pli_weight * panel_weight * gwas_weight
```

:::{tip}
Keep weights moderate (0.5–5.0 range). Extreme weights (e.g., 100) concentrate nearly all FDR
budget on one gene and leave almost nothing for discoveries in other genes. If you are certain
a gene is causal, you probably do not need a statistical test for it.
:::

:::{warning}
A weight of 0 is not allowed (division by zero). Use a small positive value like 0.1 if you
want to strongly down-weight a gene.
:::

### CLI Flags

```bash
--gene-prior-weights weights.tsv         # Path to gene-weight TSV file
--gene-prior-weight-column weight        # Column name for weights (default: "weight")
```

### CLI Example

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests logistic_burden \
  --correction-method fdr \
  --gene-prior-weights gene_weights.tsv \
  --diagnostics-output diagnostics/ \
  --output-file results.tsv
```

### Output

When `--gene-prior-weights` is provided, the association output gains one additional column:

| Column | Description |
|--------|-------------|
| `fdr_weight` | Normalized weight used for this gene (after renormalization to mean=1.0) |

When `--gene-prior-weights` is **not** provided, the `fdr_weight` column is absent (backward
compatible).

### Diagnostics

When both `--gene-prior-weights` and `--diagnostics-output` are set, an additional diagnostics
file is written:

**`fdr_weight_diagnostics.tsv`** — per-gene breakdown:

| Column | Description |
|--------|-------------|
| `gene` | Gene name |
| `raw_weight` | Weight from the input file (1.0 if absent) |
| `normalized_weight` | Weight after renormalization to mean=1.0 |
| `unweighted_q` | q-value without weighting (standard BH) |
| `weighted_q` | q-value with weighting (weighted BH) |
| `significance_change` | `gained` if newly significant, `lost` if no longer significant, blank if unchanged |

The `significance_change` column highlights genes where weighting changed the significance
call at FDR < 0.05. This is the primary way to assess the impact of your weight choices.

### Coverage Warnings

If more than 50% of tested genes are absent from the weight file, a warning is logged:

```
WARNING: 85% of tested genes (148/174) are missing from gene weight file — weights may not be informative
```

This is not an error — missing genes receive weight 1.0 (neutral) — but it suggests the weight
file may be too narrow to meaningfully redistribute FDR budget.

### Effective Number of Tests

The diagnostics output includes the effective number of tests, computed as
`sum(w)^2 / sum(w^2)`. When all weights are equal, this equals the actual number of tests. When
weights vary, it is lower, reflecting the reduced multiple testing burden from concentrating
power on fewer genes. This metric helps assess how much the weights deviate from uniform.

---

## Tuning: JSON Config

All association parameters can be set via the main JSON config file (passed with `-c`/`--config`).
Use this for reproducible analyses or when parameters are too numerous for the command line.

### Annotated Example

```json
{
  "association": {
    "association_tests": ["fisher", "logistic_burden", "skat", "coast"],

    "correction_method": "fdr",
    "gene_burden_mode": "samples",
    "trait_type": "binary",

    "covariate_file": "covariates.tsv",
    "covariate_columns": ["age", "sex", "batch"],
    "categorical_covariates": ["sex", "batch"],

    "pca_file": "ancestry.eigenvec",
    "pca_components": 10,

    "variant_weights": "beta:1,25",
    "variant_weight_params": {"cadd_cap": 40.0},

    "coast_weights": [1, 2, 3],

    "gene_prior_weights": "gene_weights.tsv",
    "gene_prior_weight_column": "weight",

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

### All 29 Valid Keys

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `association_tests` | list of str | `["fisher"]` | Tests to run: `fisher`, `logistic_burden`, `linear_burden`, `skat`, `coast` |
| `correction_method` | str | `"fdr"` | `"fdr"` (Benjamini-Hochberg) or `"bonferroni"` |
| `gene_burden_mode` | str | `"samples"` | `"samples"` (CMC/CAST carrier) or `"alleles"` (max dosage) |
| `trait_type` | str | `"binary"` | `"binary"` (case/control) or `"quantitative"` (continuous) |
| `covariate_file` | str | `null` | Path to covariate TSV/CSV |
| `covariate_columns` | list of str | `null` | Subset of covariate columns (null = all) |
| `categorical_covariates` | list of str | `null` | Columns to one-hot encode (null = auto-detect) |
| `pca_file` | str | `null` | Pre-computed PCA file (PLINK .eigenvec, AKT, or generic TSV) |
| `pca_tool` | str | `null` | `"akt"` to invoke AKT as subprocess |
| `pca_components` | int | `10` | Number of PCs to include (warn if >20) |
| `variant_weights` | str | `"beta:1,25"` | Weight scheme: `beta:a,b`, `uniform`, `cadd`, `revel`, `combined` |
| `variant_weight_params` | dict | `null` | Extra weight params, e.g. `{"cadd_cap": 30}` |
| `coast_weights` | list of float | `[1,2,3]` | COAST category weights: `[BMV, DMV, PTV]` |
| `gene_prior_weights` | str | `null` | Path to gene-weight TSV for weighted FDR correction |
| `gene_prior_weight_column` | str | `"weight"` | Column name in gene-weight file containing weights |
| `skat_backend` | str | `"python"` | `"python"`, `"auto"`, or `"r"` (deprecated) |
| `coast_backend` | str | `"python"` | `"python"`, `"auto"`, or `"r"` (deprecated) |
| `diagnostics_output` | str | `null` | Path to diagnostics directory |
| `min_cases` | int | `200` | Warn if case count below this |
| `max_case_control_ratio` | float | `20.0` | Warn if control-to-case ratio exceeds this |
| `min_case_carriers` | int | `10` | Per-gene warn if case carriers below this |
| `confidence_interval_alpha` | float | `0.05` | CI significance level (0.05 = 95% CI) |
| `continuity_correction` | float | `0.5` | Haldane-Anscombe continuity correction for zero cells in Fisher |
| `missing_site_threshold` | float | `0.10` | Exclude variants with >10% missing genotypes site-wide |
| `missing_sample_threshold` | float | `0.80` | Exclude samples with >80% missing genotypes |
| `skat_method` | str | `"SKAT"` | SKAT variant: `"SKAT"` (default), `"Burden"`, or `"SKATO"` (omnibus) |
| `confidence_interval_method` | str | `"normal_approx"` | Odds ratio CI method: tries score → normal → logit in sequence |
| `association_workers` | int | `1` | Number of parallel workers for association testing |
| `firth_max_iter` | int | `25` | Newton-Raphson iterations for Firth penalised likelihood |

### Precedence

CLI flags override JSON config values, which override `AssociationConfig` defaults:

```
CLI args > JSON "association" section > AssociationConfig defaults
```

Pass the config file via the existing `-c`/`--config` flag — there is no separate
`--association-config` flag.

:::{note}
Unknown keys in the `"association"` section raise a `ValueError` at startup with a list of
valid keys. This prevents silent misconfiguration from typos.
:::

---

## Test Selection Reference

Use this table to choose the right test(s) for your study:

| Test | CLI Name | Trait | Sample Size | What It Detects | Computational Cost | When to Use |
|------|---------|-------|-------------|-----------------|-------------------|-------------|
| Fisher's Exact | `fisher` | Binary | Any | Carrier enrichment (collapsed) | Negligible | First-pass screen; always appropriate |
| Logistic Burden | `logistic_burden` | Binary | ≥50 cases | Cumulative rare variant burden with covariate adjustment | Low | Binary trait with covariates |
| Linear Burden | `linear_burden` | Quantitative | ≥50 | Continuous phenotype burden | Low | Quantitative phenotype |
| SKAT/SKAT-O | `skat` | Binary or Quantitative | ≥100 cases | Heterogeneous effects (not all variants causal or same direction) | Moderate | Large cohorts; uncertain effect direction |
| COAST | `coast` | Binary or Quantitative | ≥100 cases | Ordered allelic series (PTV > DMV > BMV) | Moderate | Biological prior of monotone effect sizes; requires SIFT/PolyPhen |
| ACAT-O | (automatic) | Both | Depends on active tests | Omnibus combining all active tests | Negligible | Always present; primary significance measure |

**Recommended starting point for most studies:**

```bash
--association-tests fisher,logistic_burden,skat
```

This covers carrier enrichment (Fisher), covariate-adjusted regression (logistic burden), and
heterogeneous effects (SKAT-O), with ACAT-O combining all three.

Add `coast` when you have:
- SIFT and PolyPhen annotations in your VCF
- A biological prior that functional severity predicts effect magnitude
- Sufficient sample size (≥100 cases)

Use `linear_burden` instead of `logistic_burden` when:
- Your phenotype is a continuous measurement (e.g., GFR, blood pressure, biomarker level)
- You set `--trait-type quantitative`

---

## Troubleshooting

### R not found / rpy2 not installed

**Error:** `ImportError: rpy2 not installed` or `ImportError: Cannot find R installation`

**Solution:** Use the default Python backend — it requires no R installation:

```bash
--skat-backend python    # explicit
--coast-backend python   # explicit
# (These are the defaults in v0.15.0; you don't need to specify them)
```

The R backend is deprecated as of v0.15.0. If you specifically need R validation, install rpy2
and the R SKAT/AllelicSeries packages, then use `--skat-backend r`. A DeprecationWarning is
logged.

### Low Sample Size Warnings

**Warning in `warnings` column:** `LOW_CASE_COUNT`, `LOW_CARRIER_COUNT`, `IMBALANCED_COHORT`

These are triggered by thresholds in `AssociationConfig`:

| Warning | Default Threshold | Config Key |
|---------|------------------|------------|
| Total cases below threshold | 200 | `min_cases` |
| Per-gene case carriers below threshold | 10 | `min_case_carriers` |
| Control-to-case ratio exceeds threshold | 20.0 | `max_case_control_ratio` |

Adjust thresholds in JSON config if your study design intentionally has fewer cases. The analysis
still runs — these are warnings, not failures.

### Logistic Burden Convergence Failure

**Warning:** `FIRTH_CONVERGE_FAIL` in the `warnings` column, `logistic_burden_pvalue` = `None`

**What happened:** Both standard logistic regression and Firth penalised likelihood failed to
converge for this gene. Common cause: extreme imbalance (all carriers are cases), tiny carrier
count, or collinear covariates.

**Solutions:**
- Check `coast_n_ptv`/`coast_n_dmv`/`coast_n_bmv` (or `n_variants`) for the gene — very few
  carriers cannot support regression
- Remove highly correlated covariates
- Use Fisher's exact test for genes with very few carriers (it does not converge-fail)

### Convergence warnings `PERFECT_SEPARATION` or `QUASI_SEPARATION`

These appear when all (or nearly all) carriers in a gene are cases or controls. The Firth
fallback handles this automatically. If Firth also fails (`FIRTH_CONVERGE_FAIL`), the p-value
is set to `None` and the gene is excluded from ACAT-O.

### lambda_GC Inflation (lambda_GC > 1.05)

**In `diagnostics/lambda_gc.tsv`:** `lambda_gc` > 1.05

**Possible causes and actions:**

1. **Population stratification** (most common): Add PCA components via `--pca-file`
2. **Cryptic relatedness**: Verify your cohort does not include related individuals in
   case and control groups
3. **True polygenic signal**: If the trait is highly polygenic, mild inflation is expected
4. **Batch effects**: Add batch indicator to covariates via `--covariates` and
   `--categorical-covariates`

If `n_tests < 100`, the lambda_GC estimate is flagged as unreliable — do not over-interpret it.

### COAST Fails: `NO_CLASSIFIABLE_VARIANTS`

**In output:** `coast_pvalue` = `None`, and `coast_skip_reason` extra key contains `NO_CLASSIFIABLE_VARIANTS`

**What happened:** All variants in this gene received classification code 0 — they are missense
variants without SIFT or PolyPhen predictions. COAST requires at least one BMV, DMV, and PTV.

**Solutions:**
- Annotate your VCF with dbNSFP (provides `dbNSFP_SIFT_pred` and `dbNSFP_Polyphen2_HDIV_pred`)
- COAST requires PTVs — if a gene has only missense variants, COAST is not applicable
- Use SKAT or logistic burden instead

Check which annotation columns are present:
```bash
bcftools view -h input.vcf.gz | grep "##INFO" | grep -E "SIFT|Polyphen|PolyPhen"
```

### Sample ID Mismatch

**Error:** `ValueError` during covariate or PCA loading, mentioning sample IDs

**What happened:** Sample IDs in the covariate/PCA file do not match VCF sample IDs exactly
(comparison is case-sensitive).

**Fix:** Get the exact VCF sample IDs and use them in the first column of your covariate file:

```bash
bcftools query -l input.vcf.gz > vcf_samples.txt
head vcf_samples.txt
```

Ensure the first column of `covariates.tsv` matches these names exactly.

### ACAT-O Equals Single Test P-value

**What you see:** `acat_o_pvalue` is identical to `fisher_pvalue` for all genes.

**This is correct behavior.** When only one test is active, the Cauchy combination of a single
p-value returns that p-value unchanged (pass-through). FDR correction still runs across genes.

To get a true omnibus combination, enable multiple tests:
```bash
--association-tests fisher,logistic_burden,skat
```

### Missing ACAT-O Columns

ACAT-O columns are always present in the output when `--perform-association` is active. If you
do not see them, check that the output file is `results.association.tsv` (not `results.tsv`
itself — the association results go to a separate `.association.tsv` file).

---

## Comprehensive Example

The following command runs Fisher, logistic burden, SKAT-O, and COAST together with
covariate adjustment, PCA, functional weights, gene-level FDR weighting, and diagnostics:

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --phenotype-file phenotypes.tsv \
  --phenotype-sample-column sample_id \
  --phenotype-value-column case_control \
  --perform-association \
  --association-tests fisher,logistic_burden,skat,coast \
  --trait-type binary \
  --covariate-file covariates.tsv \
  --covariates age,sex,batch \
  --categorical-covariates sex,batch \
  --pca-file ancestry.eigenvec \
  --pca-components 10 \
  --variant-weights beta:1,25 \
  --coast-weights 1,2,3 \
  --gene-prior-weights gene_weights.tsv \
  --diagnostics-output diagnostics/ \
  --output-file results.tsv
```

Expected output files:

```
results.association.tsv       # Gene-level association results (includes fdr_weight column)
diagnostics/
  lambda_gc.tsv               # Genomic inflation factors
  qq_data.tsv                 # QQ plot data
  summary.txt                 # Human-readable summary
  qq_plot.png                 # QQ plot (if matplotlib installed)
  fdr_weight_diagnostics.tsv  # Per-gene weight impact analysis (when --gene-prior-weights set)
```

Representative output (first 3 columns omitted for brevity):

```
gene    n_cases  n_controls  n_variants  fisher_pvalue  ...  acat_o_qvalue  warnings
GENE1   312      488         7           0.00031             0.0047
GENE2   312      488         3           0.41                1.0
GENE3   312      488         12          8.2e-06             1.2e-04             LOW_CARRIER_COUNT
```
