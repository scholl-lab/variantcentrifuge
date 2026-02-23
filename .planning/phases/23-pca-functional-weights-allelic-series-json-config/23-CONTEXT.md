# Phase 23: PCA Integration, Functional Weights, Allelic Series, and JSON Config - Context

**Gathered:** 2026-02-21
**Status:** Ready for planning

<domain>
## Phase Boundary

Users can supply a pre-computed PCA file or trigger AKT-based PCA computation as a pipeline stage, apply CADD/REVEL functional variant weights, run the COAST allelic series test, and specify all association options via a JSON config file — with optional matplotlib QQ plot generation for users who have it installed.

</domain>

<decisions>
## Implementation Decisions

### PCA Integration

- Support three PCA file formats: PLINK `.eigenvec`, AKT kinship output, and generic TSV with sample ID + numeric PC columns (auto-detect format)
- `--pca-tool` supports AKT only (not PLINK). PLINK users pre-compute and provide via `--pca-file`
- **Hard error** when `--pca-tool akt` is specified but AKT is not in PATH — user explicitly requested PCA, so missing tool is a failure (not a skip-with-warning)
- Default to **10 PCs** when `--pca-components` not specified (standard GWAS convention matching REGENIE/SAIGE). Warn at >20 PCs
- PCs merge as additional covariate columns, sample alignment verified against VCF sample order

### Functional Weight Schemes

- CADD and REVEL scores sourced from **VCF INFO annotations** only (CADD_PHRED, REVEL_score fields). No external TSV lookup
- **CADD normalization:** Phred / 40, capped at 1.0. Combined weight = `Beta(MAF) × min(CADD_Phred / 40, 1.0)`. Cap value (40) configurable via `--variant-weight-params`
- **REVEL normalization:** Used directly (already 0-1 range). Combined weight = `Beta(MAF) × REVEL_score`
- **Variant-type-aware fallback for missing scores:**
  - LoF variants (frameshift, nonsense, splice) missing CADD/REVEL → **maximum weight (1.0)** for functional component. Rationale: LoF variants are most likely damaging; excluding them would be statistically disastrous
  - Missense variants missing scores → **uniform weight (functional component = 1.0)**, falling back to Beta(MAF)-only. No penalty, no boost
  - Other variants (synonymous, intronic) missing scores → Beta(MAF)-only
- Log warning with count of variants per category that had missing functional scores

### COAST Allelic Series

- **Classification logic:**
  - PTV (code 3): SnpEff HIGH impact — `stop_gained`, `frameshift_variant`, `splice_acceptor_variant`, `splice_donor_variant`
  - DMV (code 2): `missense_variant` (MODERATE impact) AND predicted damaging by SIFT (`deleterious`) OR PolyPhen (`probably_damaging` / `possibly_damaging`) from SnpSift dbNSFP annotations
  - BMV (code 1): `missense_variant` AND predicted benign by BOTH SIFT (`tolerated`) AND PolyPhen (`benign`)
- **Missing SIFT/PolyPhen predictions:** Missense variants without pathogenicity predictions are **excluded from COAST only** (still included in SKAT/burden tests). COAST requires clean BMV/DMV/PTV classification; adding unclassified variants would break the allelic series monotonicity assumption. Log warning per gene with count of excluded variants
- **Default category weights:** w=(1, 2, 3) for BMV, DMV, PTV matching the COAST paper (Mccaw et al., AJHG 2023). Configurable in JSON config, overridable via CLI flag
- **Output columns:** Report component + omnibus p-values matching R AllelicSeries package default:
  - `coast_burden_p_value` — burden component
  - `coast_skat_p_value` — SKAT component
  - `coast_p_value` — omnibus (feeds into ACAT-O and FDR)
  - `coast_n_bmv`, `coast_n_dmv`, `coast_n_ptv` — variant counts per category

### JSON Config

- **No new file:** Add an `"association"` section to the existing `config.json` (loaded via `-c/--config`). One config file for the whole pipeline
- **CLI overrides JSON:** JSON config is the base; CLI flags override individual settings. Standard convention
- **Fail fast with clear error** on invalid config values (unknown test names, invalid backends, bad thresholds). Validate all association config at startup, list all invalid fields in error message. No partial runs
- **Documentation only** — no `--generate-association-config` template generator. Document the schema in README/docs

### Claude's Discretion

- Exact AKT subprocess invocation and parameter passing
- PCA file format auto-detection heuristics
- CADD/REVEL VCF INFO field name variations (e.g., `CADD_PHRED` vs `CADD_phred`)
- COAST internal implementation details (score computation, null model)
- JSON config schema validation library choice
- Matplotlib QQ plot styling and layout
- Exact error message wording for config validation failures

</decisions>

<specifics>
## Specific Ideas

- CADD normalization cap of 40 chosen because CADD Phred ≥ 40 = top 0.01% most deleterious (research-backed threshold)
- CADD documentation recommends raw scores for statistical comparisons, but Phred/40 is more interpretable and practical for multiplicative weights
- REVEL only scores missense variants — LoF variants inherently have no REVEL score, so type-aware fallback is essential
- UK Biobank exome study (Backman et al., Nature 2021) used 10 common-variant PCs + 20 exome-derived PCs; 10 PCs is the conservative standard default
- COAST's allelic series assumption (BMV < DMV < PTV effect sizes) breaks if unclassified variants are shoehorned into arbitrary categories

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 23-pca-functional-weights-allelic-series-json-config*
*Context gathered: 2026-02-21*
