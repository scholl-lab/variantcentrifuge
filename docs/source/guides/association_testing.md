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
gene  n_cases  n_controls  n_variants  fisher_pvalue  fisher_qvalue  fisher_or  fisher_or_ci_lower  fisher_or_ci_upper  primary_test  primary_pvalue  primary_qvalue  warnings
GENE1 312      488         7           0.00031         0.0047         3.82       1.87                7.78                fisher        0.00031         0.0047
GENE2 312      488         3           0.41            1.0            1.23       0.74                2.05                fisher        0.41            1.0
...
```

`n_cases` and `n_controls` count only case/control samples that are present in the VCF
(listed samples without genotypes are ignored); for Fisher they are the samples with a call
for the gene.

**Primary significance measure:** `primary_qvalue`. For Fisher-only runs this is the same
value as `fisher_qvalue`; for multi-test runs it is usually the same value as
`acat_o_qvalue`.

## Input Prerequisites: Callability and Genotype QC

Association tests count every `0/0` genotype as an observed reference call. They cannot tell
whether that call was supported by sequencing data. Check the points below before you interpret
burden or association results, especially when the input VCF was merged from variant-only
per-sample VCFs.

### Reference Callability

A site-level `PASS` and a stored `0/0` genotype do not show that each sample was callable at
that position.

- **Missing-to-reference merging.** `bcftools merge -0` (`--missing-to-ref`) writes `0/0` for
  samples whose input has no record at a site. A variant-only VCF has no record both where the
  sample is confidently reference and where it has no coverage, so the filled `0/0` carries no
  callability evidence. Its per-sample fields such as `FORMAT/DP` and `FORMAT/GQ` stay missing.
- **Merged FILTER values.** With `bcftools merge -F x`, a merged record is `PASS` if any
  contributing input passes. `PASS` therefore does not mean that every retained genotype meets
  depth or quality criteria.
- **Allele counts do not fix this.** Recomputing `AC`/`AN` from the stored genotypes makes them
  consistent with those calls, but it cannot recover coverage evidence that the merge never had.

Synthetic two-sample example. `S1.vcf.gz` has a PASS record `chr1:1000 A>G` with `GT=0/1` and
`DP=35`. The two S2 rows show two different ways the cohort VCF can get a genotype for S2 at that
position: from joint genotyping, or from merging a variant-only `S2.vcf.gz` that has no record there.

| Sample | Source of the `chr1:1000` genotype | Merged `GT` | `FORMAT/DP` | Evidence |
|--------|------------------------------------|-------------|-------------|----------|
| S1 | Called in `S1.vcf.gz` | `0/1` | `35` | Observed ALT call |
| S2 | Joint genotyping or a gVCF reference block | `0/0` | `28` | Observed reference call |
| S2 | `bcftools merge -0 -F x` of variant-only VCFs | `0/0` | `.` | Filled reference; S2 may have had no coverage |

The site is `PASS` either way, and both S2 rows look the same to an association test. Only the
S2 row with observed reference evidence (`DP=28`) supports a reference call. The filled `0/0`
does not.

Before you interpret a burden result, establish reference evidence for every sample at the tested
sites. For example:

- Joint-call the cohort, or merge gVCFs, so that reference genotypes come with their own
  `DP`/`GQ`.
- Use per-sample callable-region or coverage files (for example from mosdepth or GATK
  CallableLoci). Either restrict the analysis to regions callable in all samples, or set
  genotypes outside a sample's callable regions to missing.
- Check that callability does not differ between cases and controls, for example by sequencing
  batch or capture kit. Differential callability can create or hide burden signals.

These steps are upstream assumptions; VariantCentrifuge does not introduce the merge behavior.
How to treat reference calls with missing `DP` depends on how the VCF was produced, so no single
blanket conversion rule applies.

### Genotype QC Applied by VariantCentrifuge

VariantCentrifuge does not filter individual genotypes on `DP`, `GQ` or allele balance before
association testing. (`--genotype-filter` selects carriers by genotype class, such as `het` or
`comp_het`; it is not a QC filter.) Apply genotype-level QC upstream, for example with
`bcftools filter` or `bcftools +setGT` on `FMT/DP` and `FMT/GQ`.

The following genotypes are treated as missing, which makes them subject to the missingness
thresholds and imputation below:

- fully missing calls (`./.`, `.|.`, `.`, or an empty/`NA` value);
- partial biallelic calls (`0/.`, `./1`, `1|.`);
- unparseable strings (for example `N/A` or `malformed`).

Filled `0/0` calls are not missing, so the checks below cannot detect them:

- **Fisher's exact test** uses the complete-case policy described in
  [Missing Genotypes and Complete-Case Policy](#missing-genotypes-and-complete-case-policy),
  unless `--missing-as-ref` is set (then missing calls count as reference, as in every other
  test of the run).
- **Matrix-based tests** (logistic/linear burden, SKAT/SKAT-O, COAST) and `--score-artifact`
  apply two per-gene thresholds. Set them with the CLI flags below or in the `"association"`
  section of the JSON config (the CLI wins). Values must be between 0 and 1.

| CLI flag / JSON key | Default | Effect |
|-----|---------|--------|
| `--missing-site-threshold` / `missing_site_threshold` | `0.10` | Drops a variant when more than this fraction of all samples have a missing genotype (call rate below 90%). R SKAT's default `missing_cutoff` is 0.15. |
| `--missing-sample-threshold` / `missing_sample_threshold` | `0.80` | Drops a sample from that gene's test (genotypes, phenotype and covariates) when more than this fraction of its genotypes at the kept variants are missing. |
| `--missing-as-ref` / `missing_as_ref` | off | Sets every missing call to reference (dosage 0) *before* both thresholds, so nothing is dropped for missingness or imputed. Applies to all association tests of the run, Fisher included, and to `--score-artifact`. See below. |

Allele frequencies, the imputation of the remaining missing genotypes, and minor-allele recoding
all use the retained samples only, so excluded samples cannot shift them. Missing genotypes are
imputed as `round(ploidy * alt_freq)` for binary traits and `ploidy * alt_freq` for quantitative
traits. For binary traits, a variant whose missing rate differs between cases and controls by
more than 0.05 is logged as a differential-missingness warning. Genes with fewer than 5 minor
allele copies after these steps report `NA` for the matrix-based tests.

The default imputation follows R SKAT: binary traits use `impute.method="bestguess"` (the
`SKATBinary` default) and quantitative traits use `"fixed"` (the `SKAT` default). Both are
checked against SKAT 2.2.5 on a synthetic matrix with missing calls
(`tests/unit/test_missing_policy_parity.py`).

Dropped variants are reported. When site-missingness QC removes variants, the run logs one
WARNING with the totals (`--score-artifact` logs its own). When matrix-based tests run, the
results also gain an `n_variants_qc` column. It counts the variants that passed site-missingness
QC, before the per-gene MAC check and any per-test exclusions. Per-test columns such as
`skat_n_marker_test` give the number each test actually used. `n_variants` is the count that the
first test of the run reports.

#### Merged VCFs and `--missing-as-ref`

Best-guess imputation assumes that a missing call is missing at random. That does not hold for
VCFs merged from single-sample VCFs (`bcftools merge` without `-0`): there, `./.` usually means
"no variant called here", and the observed calls are mostly the carriers. A site then has a high
missing rate and an inflated observed allele frequency. The default site threshold drops such
sites, and raising the threshold to keep them makes things worse: at an observed AF above 0.25,
`round(2·AF)` imputes every missing sample as a heterozygous carrier.

For these inputs, use `--missing-as-ref`. It gives the same result as
`bcftools +setGT -- -t . -n 0` before the run. Partial calls such as `./1` also become `0/0`, and
in `--score-artifact` the allele number counts every sample. Where `./.` is a genuine no-call
(low coverage), carriers among the no-calls are counted as reference. That biases carrier
counts toward reference. If missingness differs between cases and controls, it can also make
p-values too small, not only too large. Prefer genotypes merged with `-0` or gVCF-based joint
calling when you have them.

The packed `GT` column interacts with this as follows. When per-sample genotype columns are
rebuilt from the packed column, a sample that is not listed for a variant is filled as `0/0`. A
sample listed with an explicit `./.` stays missing and is subject to the policy above.

### Sample Subsetting at Extraction

Like PLINK `--keep`, REGENIE and SAIGE, an association run (`--perform-association`) whose
analysed samples are a strict subset of the VCF samples extracts only those samples'
genotypes (`bcftools query -S` at field extraction). The analysed set is the union of cases
and controls (`--case-samples[-file]`, `--control-samples[-file]`, or `--phenotype-file` with
`--case-phenotypes`/`--control-phenotypes`), or, for `--trait-type quantitative`, the samples
with a value in `--phenotype-value-column`. Unassigned samples then never reach the
intermediate TSV, the chunks or the genotype matrices. Association results and
`--score-artifact` scores are unchanged:

- Variant filters (`--bcftools-prefilter`, `--filters`) still run on all samples, and INFO
  fields such as `AC`/`AN`/`AF` are read as stored, not recomputed for the subset.
- A sample listed in a sample file but absent from the VCF is handled as before (it is never
  passed to bcftools).
- Per-sample columns keep the VCF sample order.

What changes: the intermediate TSV has fewer per-sample columns, and the output `GT` column
(and anything derived from it, such as the phenotypes column or HTML/Excel carrier lists)
names carriers among the analysed samples only. The SHA-256 of the sorted subset sample IDs
is recorded as `sample_subset_sha256` in the run metadata and in the score-artifact
provenance (`"none"` when every sample was kept). Edited sample/phenotype files, `--no-sample-subset`
or a change to any flag listed below make `--resume` refuse the checkpoint (configuration
changed), since the extracted samples would differ.

All samples are kept, with one log line naming the reason, when a feature reads samples
outside the analysed set: inheritance/pedigree analysis (`--ped`, `--sample-sheet`,
`--inheritance-mode`, or a scoring model that needs inheritance),
`--append-extra-sample-fields`, `--no-replacement`, `--igv`, `--pseudonymize`, a tumor/normal
sample index, a custom `--stats-config` (its expressions may read per-sample columns),
`--late-filtering` or `--final-filter`. The default statistics do not block subsetting: their
`total_samples` still reports the VCF sample count. `--no-sample-subset` turns
subsetting off.

### Multi-Allelic Sites Must Be Split

The genotype matrix encodes one ALT allele per variant. Any genotype that contains a non-missing
allele index above 1 stops the analysis with a `ValueError` that names the gene, variant and
sample column. It is not recoded. This includes fully called genotypes (`0/2`, `1/2`, `2/2`,
haploid `2`, and their phased forms), partial calls (`./2`, `2/.`, `.|2`) and polyploid calls
(`0/1/2`). Split and normalize multi-allelic records before
running VariantCentrifuge:

```bash
bcftools norm -m- -f reference.fa input.vcf.gz -Oz -o input.split.vcf.gz
```

---

## Transcript-Selected Consequence Filtering

When association or gene burden analysis uses a transcript list and consequence
predicates from SnpEff, make the consequence filter row-level by splitting SnpEff
annotations before SnpSift:

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --transcript-file mane_select_transcripts.txt \
  --preset high_or_lof_or_nmd \
  --split-snpeff-lines before_filters \
  --perform-association \
  --output-file results.tsv
```

This matters for filters that inspect the whole annotation array, such as
`ANN[ANY].IMPACT`, `LOF[*].PERC`, or `NMD[*].PERC`. In no-split mode, or with
`--split-snpeff-lines after_filters`, SnpSift evaluates those predicates on the
original multi-annotation VCF record before transcript selection. A record can
therefore pass because one gene has a qualifying annotation while transcript
filtering later retains a different gene. Use `before_filters` when the retained
row itself must satisfy the consequence predicate before entering burden or
association analysis.

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

### Covariate Pre-flight Checks

Before the gene loop, the covariates (including `--pca-file` PCs) are checked once, for every
test:

- A covariate with fewer than 5 values that differ from its most common value (for example an
  indicator equal to 1 in a single sample) is **dropped** with a warning naming it and its count.
  Such a column carries almost no information and can separate the phenotype by itself, which
  used to make every logistic burden fit fail.
- Rank-deficient or collinear covariates (condition number above 1e8 after centring and scaling)
  **stop the run** with an error naming the columns involved.
- For binary traits, a covariate whose case and control values do not overlap or touch (it
  completely separates the phenotype, so no logistic model has a finite estimate) **stops the run**
  with an error naming it. Quasi-separation, such as a batch with controls only, is left to the
  Firth fallback (and to the batch warning of `--diagnostics-output`).

Well-conditioned covariates pass through unchanged. For the logistic burden fit only, covariates
are centred and scaled when `[1, covariates]` has a condition number above 1e6 (for example
values in very different units); the burden beta, SE and p-value do not depend on this, so beta
stays on the dosage scale.

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

**AKT stdout** (`akt pca ... > pca.txt`; no header, a tab after the sample ID, then
space-padded values):
```
SAMPLE_001	  1.26366 -0.243411    1.4646
SAMPLE_002	-0.512342   2.03114  -1.10982
```

**Generic TSV/CSV** (one sample ID column then numeric; header optional):
```
sample_id   PC1       PC2       PC3
SAMPLE_001  0.0142    -0.0231   0.0089
SAMPLE_002  -0.0305   0.0118    -0.0042
```

Files ending in `.csv` are comma-separated. Other files are split on tabs when every line
has the same number of tabs and no value is space-padded (so IDs and header names may contain
spaces); otherwise on any mix of tabs and spaces. Headerless files whose first two columns could
both match some VCF samples need a header line; column 1 equal to column 2 on every row is read
as PLINK `FID IID`. Duplicate sample rows are
an error. A header is assumed when the first value
after the sample ID column(s) is not a number. Any NaN or non-numeric value in a used PC column
stops the run with an error naming the file, the columns and example samples.

PLINK format uses IID (column 2) as the sample identifier — this matches VCF sample names.

### CLI Flags

```bash
--pca ancestry.eigenvec         # Pre-computed PCA file (or --pca akt for on-the-fly)
--pca-components 10             # Number of PCs to include (default: 10, warn if >20)
--pca-sites sites.vcf.gz        # Optional VCF/BCF file of high-quality variant sites for PCA
```

PCA columns are automatically appended to the covariate matrix before any regression. If you
specify both `--covariate-file` and `--pca ancestry.eigenvec`, PCs are added as additional columns.

:::{warning}
Using `--pca akt` requires AKT to be in your PATH. A missing AKT binary raises a hard error
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

### Zero-Cell Handling and Continuity Correction

When any cell in the 2×2 contingency table is zero (common for rare variants where no controls or no cases carry a variant), the sample odds ratio would be zero or undefined (infinite). VariantCentrifuge applies the Haldane-Anscombe continuity correction (adding `continuity_correction`, default `0.5`, to all four cells) to evaluate both the reported odds ratio (`fisher_or`) and its normal-approximation confidence interval (`fisher_or_ci_lower`, `fisher_or_ci_upper`) from the same corrected table via `statsmodels.stats.contingency_tables.Table2x2`. If marginal totals are zero (structural zero), the odds ratio and confidence intervals return `NaN`.

### Missing Genotypes and Complete-Case Policy

VariantCentrifuge applies principled missingness handling across association tests:
(see [Input Prerequisites](#input-prerequisites-callability-and-genotype-qc) for what counts
as missing and why filled reference calls are not):

- **Fisher's Exact Test (Complete-Case Policy):** Samples with missing genotype calls (`./.`, `.`, `.|.`) at any qualifying variant in a gene are excluded from the test's called counts (`proband_called_count`, `control_called_count`) and reference denominator (`proband_called_chromosomes`, `control_called_chromosomes`). Samples with missing calls are not assumed to be reference (`0/0`) unless `--missing-as-ref` is set. If all cases or all controls are uncalled for a gene, the test is safely skipped.
- **Matrix-Based Tests (Burden Regression, SKAT, SKAT-O):** Genotype matrices use ploidy-scaled mean imputation (`ploidy * alt_allele_freq`, i.e. `2 * alt_freq` for diploid and `1 * alt_freq` for haploid; or `round(ploidy * alt_allele_freq)` bestguess for binary traits). Dosage units are `{0, 1, 2}` for diploid variants and `{0, 1}` for haploid calls.

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

**Convergence failures:** a gene where both Logit and Firth fail is reported as NA with
`FIRTH_CONVERGE_FAIL` plus the reason Logit failed: `CONVERGENCE_SEPARATION` (BSE > 100),
`CONVERGENCE_MAX_ITER`, `CONVERGENCE_SINGULAR_HESSIAN` or `CONVERGENCE_FIT_ERROR`. The stage logs
the count per reason, and **fails** (non-zero exit, no results table) when one reason hits more
than 50 % of the tested genes, listing the reason counts and the covariates whose null-model
standard error exceeds 100.

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

**SKAT p-values.** The Python backend computes the mixture-of-chi-square tail with the compiled
Davies method (`qfc`), starting at R's settings (absolute accuracy 1e-6, 10,000 terms) and
refining the accuracy so the error stays within about 1e-4 of the p-value. Below p = 1e-5,
where Davies' absolute error is too coarse, the tail is computed with Davies on the
exponentially tilted mixture (`skat_p_method = davies_tilted`), which keeps relative accuracy
in the deep tail. The Kuonen saddlepoint and the Liu moment-matching approximation are used
only as fallbacks; they are reported with `skat_p_converged = False`.

**SKAT-O p-values.** The omnibus p-value integrates over the chi-square(1) mixing variable in
survival form (no "1 − ..." cancellation and no truncation at 40) with an adaptive outer
integral and an error estimate. If Davies or the integral cannot meet the accuracy budget, the
Liu approximation is used instead and the result is reported with `skat_p_converged = False`.

**Agreement with R.** SKAT, SKAT-O and Burden match R SKAT 2.2.5 for moderate p-values. In
the deep tail and in some edge cases (Davies at R's 1e-6 absolute accuracy, R's
`integrate(0, 40)` outer integral, or R falling back to Liu after Davies' 10,000-term limit)
R itself is less accurate, and results can differ from R by more than 10%. These discrepancies
were checked against independent references (importance sampling, high-accuracy Davies, a fine
grid of R's own integrand), which agree with VariantCentrifuge.

### Imbalance-Robust SKAT (`--skat-robust`)

With few cases or few controls and ultra-rare variants, the asymptotic SKAT/SKAT-O null is
anti-conservative for binary traits (Zhao et al. 2020, AJHG). `--skat-robust` (or the
config-file key `skat_robust: true`) switches binary SKAT/SKAT-O to a port of R
`SKAT::SKATBinary_Robust`:

- each variant's score is tested with the normal approximation; where |z| > 2 the p-value is
  replaced by a saddlepoint p-value (MAC ≥ 10) or an exact p-value from enumerating the
  carriers' case/control labels (MAC < 10, as R — not ≤ 10; mid-p);
- each variant's score variance is recalibrated to `V~ = S² / qchisq(1 − p~, 1)`, the score
  covariance is rescaled to these variances, and the burden direction is recalibrated with a
  saddlepoint burden p-value;
- SKAT or SKAT-O is then computed from the calibrated score and covariance, with the same
  Davies / SKAT-O integration as the default path.

Results are reported in the usual `skat_*` columns with `skat_p_method = spa_er`. Genes with a
single informative variant get a p-value (R behaviour) instead of being skipped. The option
needs `--trait-type binary`, `skat` in `--association-tests`, the Python backend and
`--skat-method SKAT` or `SKATO`; other combinations fail before extraction. It is off by
default; the default output is unchanged.

**Agreement with R.** On synthetic 1:9 and 8:1 cohorts, the per-variant saddlepoint/exact
p-values and the calibrated score/covariance match `SKATBinary_Robust` (SKAT 2.2.5, SPAtest
3.1.2) to 1e-11 relative. Gene SKAT p-values agree with R's Davies at 1e-12 accuracy
on R's own calibrated statistics (R's printed p differs only where its 1e-6 absolute Davies
falls back to Liu). Gene SKAT-O p-values use the same SKAT-O omnibus as the default path;
they differ from R's printed p by up to 4.7e-4 absolute (relative up to 6e-3 observed, at
p ≈ 1e-4, where the difference is below 1e-6 absolute),
which comes from R's SKAT-O integration (Davies at 1e-6 absolute inside the integrand),
not from the calibration (`scripts/generate_skat_robust_parity.R`).

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

:::{warning}
ACAT-V is anti-conservative for ultra-rare variants when cases and controls are imbalanced,
unless it is saddlepoint (SPA) corrected; VariantCentrifuge has no SPA correction. A single
ultra-rare variant can then dominate ACAT-O and `primary_qvalue`. In imbalanced cohorts use
`--acat-o-exclude acat_v` (the case:control imbalance warning suggests it), and consider
`--primary-test skato` or `--primary-test burden`. See [Choosing the Primary Test](#choosing-the-primary-test).
:::

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
| Auto | `--coast-backend auto` | Probes environment; defaults to pure Python |

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

COAST tests the variants kept by genotype site QC (`--missing-site-threshold`), the same
variants burden and SKAT test.

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

ACAT-O (Aggregated Cauchy Association Test — Omnibus) combines the gene's set-test p-values
into a single omnibus p-value using the Cauchy combination: the burden test (`logistic_burden` or
`linear_burden`), SKAT/SKAT-O, ACAT-V, and COAST when it runs. This captures signal regardless
of which test is most powerful for a given gene. The Cauchy combination is particularly robust to
correlation between test statistics.

**Fisher's exact test is not part of ACAT-O.** ACAT-O is defined over continuous set tests
(burden, SKAT, ACAT-V; Liu et al. 2019 AJHG, and STAAR). Fisher is a discrete test: its p-value
has a point mass at exactly 1 (for example a gene with equal case and control carrier rates), and
in the Cauchy combination one p = 1 vetoes the whole omnibus. Fisher is still computed and
reported in its own `fisher_*` columns. In a run of `fisher` plus one other test, ACAT-O is that
test's p-value (single-input pass-through) and stays the primary result; `fisher` alone has no
ACAT-O.

The combination follows the reference R `ACAT()`: a p-value of exactly 1 is kept, and it
dominates (tan((0.5 - 1)·π) ≈ -1.6e16), so a gene where any combined test has p = 1 gets an ACAT-O
of about 1, even if another test is small. The same holds inside ACAT-V: every marginal variant
p-value in [0, 1] is passed on, so a variant whose score is exactly 0 (p = 1) gives ACAT-V = 1,
as in STAAR. Tests that did not run are left out. When one input is exactly 0 and another exactly
1 the combination is undefined and ACAT-O is empty (R `ACAT()` stops with an error).

**References:** Liu et al. 2019 AJHG; Liu and Xie 2020 JASA; Li et al. 2020 Nat Genet (STAAR).

ACAT-O is **not** a test you activate separately — it is computed automatically after
multiple selected tests have run. When SKAT is active, ACAT-V is also included in the
combination unless `--acat-o-exclude acat_v` is given.

### FDR Correction Strategy

A **single** Benjamini-Hochberg FDR correction pass is applied to the primary p-values across
all genes. By default, in multi-test and SKAT runs the primary p-value is ACAT-O
(`--primary-test` changes this). In single-test runs
such as Fisher-only, the primary p-value is that test's own p-value, and the corrected value is
also exposed under a test-specific column such as `fisher_qvalue`.

:::{warning}
Do not apply additional correction to the output p-values, and do not apply FDR separately per
test. Both approaches are statistically incorrect given this design. Use `primary_qvalue` as
your primary significance measure.
:::

### Output Columns

| Column | Description |
|--------|-------------|
| `primary_test` | Name of the result used for correction (`fisher`, `acat_o`, etc.) |
| `primary_pvalue` | Raw p-value used for the single correction pass |
| `primary_qvalue` | Corrected primary p-value; use this for significance decisions |
| `<test>_qvalue` | Corrected p-value for a single-test run, e.g. `fisher_qvalue` |
| `acat_o_pvalue` | Raw ACAT-O omnibus p-value; present whenever ACAT-O combines tests |
| `acat_o_qvalue` | Corrected ACAT-O p-value; present when ACAT-O is the primary result |

### Choosing the Primary Test

`--primary-test` sets which result drives `primary_*` and the single FDR pass:

| Value | Primary result | `primary_test` |
|-------|----------------|----------------|
| `acat_o` (default) | ACAT-O for multi-test or SKAT runs; the test itself in single-test runs | `acat_o`, or the test name |
| `skato` | SKAT-O; needs `skat` in `--association-tests` and `--skat-method SKATO` | `skat_o` |
| `burden` | `logistic_burden` (binary trait) or `linear_burden` (quantitative); must be in `--association-tests` | that test name |

With `skato` or `burden`, the chosen test's corrected value is also written as
`<test>_qvalue` (`skat_o_qvalue`, `logistic_burden_qvalue`). ACAT-O is still computed and
`acat_o_pvalue` is kept, but it is not corrected, so `acat_o_qvalue` is absent. If the chosen
test is not run, the run stops before extraction with an error naming the missing test.

`--acat-o-exclude acat_v` leaves ACAT-V out of ACAT-O (see the ACAT-V warning above). In a
SKAT-only run ACAT-O would then be a pass-through, so SKAT itself becomes the primary result
(`skat_qvalue` or `skat_o_qvalue`) and no `acat_o_*` columns are written.

Both options also apply to gene-set units, and can be set in the config file's `association`
section as `"primary_test"` and `"acat_o_exclude"` (the CLI wins). A config-file
`"acat_o_exclude"` can only be removed in the config file; there is no CLI value that unsets it. They do not affect
`--score-artifact` files or meta-analysis, which has its own primary (SKAT-O, else burden).
Defaults reproduce earlier output exactly.

### Pass-Through Behaviour

For single-test runs without a real omnibus combination, VariantCentrifuge no longer writes
ACAT-O pass-through columns. For example, Fisher-only output uses `fisher_pvalue`,
`fisher_qvalue`, and the generic `primary_*` columns.

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
| Score column | `column:<name>` or `score_column` | Beta(MAF) × any numeric per-variant score column, or raw score-only weights |

### CLI Flags

```bash
--variant-weights beta:1,25              # Weight scheme (default: beta:1,25)
--variant-weight-params '{"cadd_cap":30}'  # JSON string with extra parameters
```

`--variant-weight-params` accepts a JSON object for scheme-specific tuning.
Supported keys depend on the selected weight scheme. For CADD-based schemes,
`cadd_cap` (float, default 40.0) caps the CADD phred score before
normalisation. For score-column weights, use the parameters in the table below.

### Score-Column Variant Weights

Variant-level weights affect each variant's contribution to burden scores and
SKAT kernels. They are separate from gene-level FDR prior weights, which affect
multiple-testing correction only.

Use `column:<name>` to reference any numeric per-variant column already present
in the association input DataFrame:

```bash
--perform-association \
--association-tests logistic_burden,skat \
--variant-weights column:nephro_candidate_score \
--variant-weight-params '{"score_min":0,"score_max":10,"floor":0.1,"combine_with_beta":true}'
```

The equivalent friendly form is:

```bash
--perform-association \
--association-tests logistic_burden,skat \
--variant-weights score_column \
--variant-weight-column nephro_candidate_score \
--variant-weight-params '{"score_min":0,"score_max":10,"floor":0.1,"combine_with_beta":true}'
```

Score-column parameters:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `score_min` | `null` | Lower bound for explicit linear normalization. |
| `score_max` | `null` | Upper bound for explicit linear normalization. |
| `floor` | `0.0` | Lower bound after clipping finite normalized scores. |
| `ceiling` | `1.0` | Upper bound after clipping finite normalized scores. |
| `combine_with_beta` | `true` | Multiply the functional score by `Beta(MAF; beta_a,beta_b)`. |
| `missing` | `null` | `neutral` in Beta-combined mode, `floor` in raw score-only mode. |
| `beta_a` | `1.0` | Beta alpha used when `combine_with_beta=true`. |
| `beta_b` | `25.0` | Beta beta used when `combine_with_beta=true`. |

When `combine_with_beta=false`, the normalized score itself is the final variant
weight. In that raw score-only mode, `missing="neutral"` is invalid because
there is no multiplicative Beta baseline.

Python SKAT, SKAT-O, and ACAT-V use the same resolved weight vector as burden
tests. For `beta:*`, Python SKAT preserves its historical MAF source:
`geno.mean(axis=0) / 2.0` after genotype imputation. Python SKAT now honors
`cadd`, `revel`, and `combined`; earlier versions silently reduced those
schemes to Beta-only weights in SKAT.

Fisher exact tests ignore variant weights because they use collapsed contingency
tables. COAST uses its own ordered allelic-series category weights. The
deprecated R SKAT backend supports only `beta:a,b` and `uniform`; use the
default Python backend for score-column or functional weights.

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
| `covariate_case_control.tsv` | TSV: `covariate`, `level`, `n_cases`, `n_controls` | Cases/controls per level of each categorical covariate (raw labels, pre-encoding); written when `--categorical-covariates` is used, to expose batch confounding |
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
When weighted FDR is active, the main association output also includes explicit weighted and
unweighted aliases for the primary result, such as `fisher_weighted_qvalue`,
`fisher_unweighted_qvalue`, `primary_weighted_qvalue`, and `primary_unweighted_qvalue`.

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

## Gene Sets as Testing Units

`--gene-set-file sets.tsv` tests predefined gene sets (for example a positive-control set of
known disease genes) in addition to the single genes. The file is a TSV with the columns
`set_id` and `gene`, one row per member gene:

```text
set_id	gene
known_disease	PKD1
known_disease	PKD2
```

Each set is one extra testing unit: the union of its member genes' qualifying variants, with a
variant annotated to several member genes counted once. It runs through the same tests,
weights and covariates as the genes, so a one-gene set gives exactly that gene's result and a
two-gene set gives the result of one gene holding both genes' variants.

- **Shared variants.** A variant annotated to several member genes enters the set once and
  keeps the annotation of its first row in input order, which can change its COAST class or
  functional weight within the set compared with the other gene.
- **Run time.** Cost grows with sets × variants: each set rebuilds its genotype matrix and
  reruns every test (for example about 20 s extra for 200 sets on 20,000 variants × 2,000
  samples).
- **Output.** Set rows carry the set id in the `gene` column. A `unit_type` column (`gene` or
  `set`) is added after `gene` only when `--gene-set-file` is given.
- **Multiple testing.** Genes and sets are separate correction families: the gene rows and
  their q-values are identical to a run without sets, and sets are corrected among themselves.
  `--gene-prior-weights` applies to genes only; sets are corrected unweighted. The
  `--diagnostics-output` QQ/λ_GC files describe genes only.
- **Score artifact.** With `--score-artifact`, each set is stored as an ordinary unit under its
  set id, so `variantcentrifuge meta` combines it like a gene. `meta` has no unit type, so it
  corrects sets and genes together.
- **Validation.** Missing `set_id`/`gene` columns, an empty or `::`-containing `set_id`,
  duplicate `set_id`/`gene` rows and a `set_id` equal to a gene symbol in the data are errors.
  Member genes without qualifying variants (or unknown) are warned with their count; a set with
  no such gene is skipped with a warning.

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

    "variant_weights": "score_column",
    "variant_weight_column": "nephro_candidate_score",
    "variant_weight_params": {
      "score_min": 0,
      "score_max": 10,
      "floor": 0.1,
      "combine_with_beta": true
    },

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
| `variant_weights` | str | `"beta:1,25"` | Weight scheme: `beta:a,b`, `uniform`, `cadd`, `revel`, `combined`, `column:<name>`, `score_column` |
| `variant_weight_column` | str or null | `null` | Score column used when `variant_weights` is `score_column` |
| `variant_weight_params` | dict | `null` | Extra weight params, e.g. `{"cadd_cap": 30}` |
| `coast_weights` | list of float | `[1,2,3]` | COAST category weights: `[BMV, DMV, PTV]` |
| `gene_prior_weights` | str | `null` | Path to gene-weight TSV for weighted FDR correction |
| `gene_prior_weight_column` | str | `"weight"` | Column name in gene-weight file containing weights |
| `skat_backend` | str | `"python"` | `"python"`, `"auto"`, or `"r"` (deprecated) |
| `coast_backend` | str | `"python"` | `"python"` or `"auto"` |
| `diagnostics_output` | str | `null` | Path to diagnostics directory |
| `min_cases` | int | `200` | Warn if case count below this |
| `max_case_control_ratio` | float | `20.0` | Warn if the larger/smaller of cases and controls exceeds this (either direction; the message names which group is larger) |
| `min_case_carriers` | int | `10` | Per-gene warn if case carriers below this |
| `confidence_interval_alpha` | float | `0.05` | CI significance level (0.05 = 95% CI) |
| `continuity_correction` | float | `0.5` | Haldane-Anscombe continuity correction for zero cells in Fisher |
| `missing_site_threshold` | float | `0.10` | Exclude variants with >10% missing genotypes site-wide |
| `missing_sample_threshold` | float | `0.80` | Exclude samples with >80% missing genotypes |
| `missing_as_ref` | bool | `false` | Treat missing genotypes as reference before missingness QC (`--missing-as-ref`) |
| `skat_method` | str | `"SKAT"` | SKAT variant: `"SKAT"` (default), `"Burden"`, or `"SKATO"` (omnibus) |
| `confidence_interval_method` | str | `"normal_approx"` | Odds ratio CI method: tries score → normal → logit in sequence |
| `association_workers` | int | `0` | Worker processes for the per-gene tests (`--association-workers`). `0` = one per `--threads`, `-1` = all cores, `1` = sequential. Each worker uses one BLAS thread; panels under 50 genes and R backends run sequentially. Parallel runs need Linux (fork); elsewhere genes run sequentially. Results are identical for any worker count |
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
| ACAT-O | (automatic) | Both | Depends on active tests | Omnibus of burden, SKAT/SKAT-O, ACAT-V (and COAST); not Fisher | Negligible | Always present; primary significance measure |

**Recommended starting point for most studies:**

```bash
--association-tests fisher,logistic_burden,skat
```

This covers carrier enrichment (Fisher), covariate-adjusted regression (logistic burden), and
heterogeneous effects (SKAT-O), with ACAT-O combining logistic burden, SKAT-O and ACAT-V
(Fisher is reported on its own).

Add `coast` when you have:
- SIFT and PolyPhen annotations in your VCF
- A biological prior that functional severity predicts effect magnitude
- Sufficient sample size (≥100 cases)

**Quantitative traits** (a continuous measurement such as GFR, blood pressure or a biomarker
level) use a different combination:

```bash
--trait-type quantitative --association-tests linear_burden,skat
```

`fisher` and `logistic_burden` are binary-trait tests and are rejected at startup with
`--trait-type quantitative`. Because `fisher` is the default, `--association-tests` must be
given explicitly. `coast` may be added.

- **Outcome.** The outcome is the numeric value in `--phenotype-value-column` of
  `--phenotype-file`, matched to VCF samples through `--phenotype-sample-column` (after
  `--remove-sample-substring`, if set).
- **Missing outcome.** Samples whose value is empty, `NA`, `NaN` or `.` are excluded from the
  genotypes, the outcome and the covariates, and the number excluded is logged. Any other
  non-numeric or non-finite value, or conflicting values for one sample, stops the run.
- **No case/control labels.** Case/control definitions (`--case-phenotypes`,
  `--case-samples` and so on) are neither required nor used, and `n_cases`/`n_controls` are
  reported as 0.
- **No score artifact.** `--score-artifact` (cohort meta-analysis) supports binary traits only.

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
| Case/control imbalance (larger:smaller, either direction) exceeds threshold | 20.0 | `max_case_control_ratio` |

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

### Untestable Genes: Skip and Warning Codes

When a gene cannot be tested reliably, the affected test reports no p-value (empty
`{test}_pvalue`) and a code that says why, rather than a spurious p-value. The gene is left
out of FDR correction and ACAT-O for that test. Codes appear in these columns:

- `skat_skip_reason` and `skat_warnings`: SKAT / SKAT-O. Codes from the null-model fit (upper
  case, e.g. `CONSTANT_PHENOTYPE`) appear in both columns; per-gene kernel codes (lower case,
  e.g. `rank_deficient`) appear in `skat_skip_reason`.
- `warnings`: linear and logistic burden (and the per-gene carrier warnings)
- `coast_skip_reason`: COAST

| Code | Column | Meaning | What to do |
|------|--------|---------|------------|
| `NONFINITE_PHENOTYPE` | `skat_skip_reason`, `skat_warnings`, `warnings`, `coast_skip_reason` | The outcome contains NaN or infinite values. | Check the phenotype values passed to the analysis. |
| `CONSTANT_PHENOTYPE` | `skat_skip_reason`, `skat_warnings`, `warnings`, `coast_skip_reason` | The outcome has no variance among the analysed samples. | Check the phenotype column and how many samples have a value. |
| `NONFINITE_COVARIATES` | `warnings` (linear burden) | The covariate matrix contains NaN or infinite values. | Fix or drop the affected covariate values. |
| `ill_conditioned_covariates` | `skat_skip_reason` | Covariates are too ill-conditioned to decide which variants carry information beyond them. | Remove redundant or near-collinear covariates (e.g. duplicated columns, nested categories, too many PCs). |
| `ILL_CONDITIONED_COVARIATES` | `warnings` (linear burden), `coast_skip_reason` | Same as above for the burden or the allelic SKAT-O projection. | As above. |
| `no_informative_variants` | `skat_skip_reason` | Every qualifying variant is fully explained by the covariates. | Nothing to fix; the gene has no testable signal under this model. |
| `BURDEN_NOT_ESTIMABLE` | `warnings` (linear burden) | The burden score is constant or lies in the span of the intercept and covariates. | Typical for genes with very few carriers; use SKAT or relax variant filters. |
| `no_positive_eigenvalues` | `skat_skip_reason` | Binary trait: the covariate-adjusted SKAT kernel has no eigenvalue above round-off, so no p-value exists (R SKAT stops with "No Eigenvalue is bigger than 0!!"). | Nothing to fix; the variants carry no information beyond the covariates. |
| `rank_deficient` | `skat_skip_reason` | The genotype matrix has rank below 2 (for example a single variant, or identical carrier patterns). | Expected for single-variant genes; use burden or Fisher. |
| `zero_variance` | `skat_skip_reason` | Burden-method SKAT: the collapsed burden is in the covariate span. | As for `BURDEN_NOT_ESTIMABLE`. |
| `OLS_FIT_FAILED` / `NONFINITE_STATISTICS` | `warnings` (linear burden) | The regression raised, or returned a non-finite estimate, SE, p-value or CI. | Check for extreme outcome or covariate values. |
| `OUTCOME_EXPLAINED_BY_COVARIATES` | `skat_skip_reason`, `skat_warnings`, `warnings` (linear burden), `coast_skip_reason` | The covariates explain the outcome exactly, so no residual variation is left to test. | Remove covariates derived from the outcome. |
| `INSUFFICIENT_RESIDUAL_DF` | `skat_skip_reason`, `skat_warnings`, `coast_skip_reason` (linear burden reports `OUTCOME_EXPLAINED_BY_COVARIATES` when the covariates alone use up the samples, or `NONFINITE_STATISTICS` when adding the burden does) | Too few samples remain (after per-gene missingness QC) for the number of covariates. | Use fewer covariates or PCs, or more samples. |

COAST also reports structural skips in `coast_skip_reason`: `NO_GENOTYPE_MATRIX`,
`NO_PHENOTYPE_OR_EMPTY_MATRIX`, `NO_GENE_DF`, `MISSING_EFFECT_IMPACT_COLUMNS`,
`ANNOTATION_GENOTYPE_MISMATCH`, `NO_CLASSIFIABLE_VARIANTS` (see below), `ALL_CATEGORIES_EMPTY`
and `BACKEND_ERROR:<message>`.

A p-value that comes from an approximation rather than converged Davies is still reported, but
with `skat_p_converged = False` (see [SKAT p-values](#backends)).

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

### Missing ACAT-O Columns In Fisher-Only Output

**What you see:** Fisher-only output contains `fisher_qvalue` and `primary_qvalue`, but no
`acat_o_pvalue` or `acat_o_qvalue`.

**This is correct behavior.** ACAT-O columns are written only when ACAT-O is the primary result.
For Fisher-only runs, the corrected value is the Fisher result, so it is written as
`fisher_qvalue` and mirrored in `primary_qvalue`.

To get a true omnibus combination, enable multiple tests:
```bash
--association-tests fisher,logistic_burden,skat
```

### Missing ACAT-O Columns

ACAT-O columns are present for multi-test and SKAT runs. If you expect them but do not see
them, first check that the output file is `results.association.tsv` (not `results.tsv` itself
— the association results go to a separate `.association.tsv` file), then check the selected
tests in `--association-tests`.

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
  --pca ancestry.eigenvec \
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
