# Phase 19: Covariate System and Burden Tests - Context

**Gathered:** 2026-02-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Enable logistic and linear regression burden tests with covariate adjustment, built on a genotype matrix builder that correctly handles all real-world genotype encodings. Users run `--perform-association --association-tests logistic_burden` (or `linear_burden`) with a covariate file and get regression-based p-values with effect sizes. Weight schemes (Beta(MAF;1,25), uniform) and trait type specification are included.

</domain>

<decisions>
## Implementation Decisions

### Covariate file format & handling
- Auto-detect TSV vs CSV from file content or extension (support both)
- First column = sample ID, remaining columns = covariates, header row required
- Samples in covariate file NOT in VCF: warn (log WARNING listing extra sample IDs) but continue with intersection
- VCF samples MISSING from covariate file: error and abort — every VCF sample must have covariates
- Categorical covariates: auto-detect (column with <=5 unique non-numeric values), with `--categorical-covariates` CLI flag to override (force categorical or numeric)
- One-hot encode categorical covariates automatically; drop one level to avoid multicollinearity

### Genotype edge cases & missing data
- **Two-layer missing genotype strategy:**
  1. Pre-filter: exclude variants with >10% missing calls site-wide; exclude samples missing >80% of variants in a gene
  2. Imputation: mean imputation (2x MAF computed from ALL samples, not stratified by phenotype) for continuous traits; bestguess (round to nearest 0/1/2) for binary traits — matches SKAT R package convention
- Differential missingness warning: flag variants where case vs control missing rates differ by >5% absolute
- Never silently treat ./. as dosage 0 — always apply the imputation strategy
- **Multi-allelic genotypes (1/2):** count as dosage=1 (het-equivalent) + emit WARNING telling user to run `bcftools norm -m-`. This is a defensive fallback — the recommended pipeline splits multi-allelics upstream
- **Phased genotypes (0|1 vs 1|0):** treat identically to unphased — phase information irrelevant for burden/SKAT (sum alleles to 0/1/2 dosage)

### Burden test output & effect sizes
- **Logistic burden (binary traits):** p_value, odds_ratio, ci_lower, ci_upper, beta, se, n_variants, n_carriers, n_carriers_cases, n_carriers_controls
- **Linear burden (quantitative traits):** p_value, beta, se, ci_lower, ci_upper, n_variants, n_carriers
- **Convergence warnings:** Add a `warnings` column to output with structured codes:
  - `PERFECT_SEPARATION` — all carriers are cases or all are controls
  - `QUASI_SEPARATION` — near-perfect separation (>=90% carriers in one group)
  - `FIRTH_CONVERGE_FAIL` — Firth correction attempted and failed
  - `LOW_CARRIER_COUNT` — total carriers below threshold
  - `ZERO_CARRIERS_ONE_GROUP` — no carriers in cases or controls
- **Firth fallback:** Always attempt Firth's penalized likelihood regression when standard logistic regression has separation/convergence issues. Report Firth p-value if it converges; report p_value=NA with warning code if Firth also fails
- Never report a naive non-converged p-value without a warning flag
- Never silently omit genes — always include the row with NA + warning code

### Default behaviors & guard rails
- **Default weighting:** Beta(MAF; 1, 25) for all test types (burden, SKAT, SKAT-O) — field standard used by SKAT R, regenie, SAIGE-GENE+
  - `--variant-weights uniform` for equal weights
  - `--variant-weights beta:a1,a2` for custom Beta parameters
- **Trait type:** Default to binary (HPO-based workflow). Require explicit `--trait-type quantitative` for continuous traits. Log a confirmation of which type is being used
- **Phenotype encoding:** 0/1 (control/case) — matches SKAT, regenie, SAIGE conventions
- **Full tiered warning system:**

  | Condition | Level | Action |
  |-----------|-------|--------|
  | Cases < 10 | ERROR | Refuse to run; results invalid |
  | Cases < 50 | WARN (high) | No practical power; note prominently |
  | Cases < 200 | WARN | Underpowered for SKAT; interpret with caution |
  | Case:control > 1:20 | WARN | Type I error inflation risk without SPA/Firth |
  | Per-gene MAC < 5 | SKIP | Insufficient data; gene excluded |
  | Per-gene case carriers < 10 | NOTE | Fall back to Fisher's exact instead of regression |
  | Per-gene case carriers < 3 | WARN | Extremely sparse; single-observation signal |
  | Total N < 1000 (binary) | WARN | Must use small-sample SKAT correction (MA/ER) |

### Claude's Discretion
- Exact implementation of delimiter auto-detection heuristic
- Categorical auto-detection threshold tuning (<=5 unique values is starting point)
- Multicollinearity check implementation details
- Firth regression library choice (statsmodels vs custom)
- Exact log message wording and formatting
- Internal genotype matrix storage format (dense vs sparse)
- Mean imputation implementation details (per-variant vs per-gene MAF calculation)

</decisions>

<specifics>
## Specific Ideas

- Covariate file format should feel familiar to users of PLINK .cov files and regenie phenotype files
- Warning column in output should follow PLINK2's ERRCODE pattern — structured codes, not free text
- Beta(MAF;1,25) weighting formula: `w = 25 * (1-p)^24` — should be documented in output headers or help text
- Two-layer missing genotype strategy matches SKAT R conventions: `impute.method="fixed"` for continuous, `"bestguess"` for binary traits
- The pipeline should work end-to-end for small exome studies (50-200 cases) common in rare disease research, with appropriate warnings but not blocking execution

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 19-covariate-system-and-burden-tests*
*Context gathered: 2026-02-19*
