# Phase 22: ACAT-O and Diagnostics - Context

**Gathered:** 2026-02-21
**Status:** Ready for planning

<domain>
## Phase Boundary

Users receive omnibus ACAT-O p-values combining burden and SKAT results per gene, a single FDR-corrected set of ACAT-O p-values across all genes, and a diagnostics directory containing lambda_GC per test, QQ plot data TSV, and a human-readable summary. ACAT-V (variant-level Cauchy) and matplotlib QQ plots are out of scope (Phase 23).

</domain>

<decisions>
## Implementation Decisions

### Output columns & TSV layout
- **Fixed schema with NA for inactive tests** — all test columns always present regardless of which tests were run. Downstream scripts can hardcode column names. Matches SAIGE-GENE and REGENIE conventions.
- **Wide format (one row per gene)** — columns like fisher_p_value, burden_p_value, skat_p_value, skat_o_p_value, acat_o_p_value all on one row. Matches SAIGE-GENE style.
- **Prefix naming pattern** — `{test_name}_{metric}` (e.g., `fisher_p_value`, `burden_beta`, `skat_p_value`). Consistent with existing engine.py column assembly at line 251.
- **Default sort by acat_o_corrected_p_value** — genes sorted by omnibus corrected p-value to surface strongest signals. Excel conditional formatting keyed on this column.

### Diagnostics directory structure
- **Single lambda_gc.tsv** with a `test_name` column — one file, rows like: fisher | 1.02, burden | 0.98, skat | 1.01.
- **Single qq_data.tsv** with a `test_name` column — filter by test for plotting. Consistent with lambda_gc.tsv layout.
- **Main output only for association results** — diagnostics directory contains only diagnostics files (lambda_gc.tsv, qq_data.tsv, summary.txt). Association results TSV lives in the standard output location.
- **Include summary.txt** — human-readable overview with sample sizes, active tests, lambda_GC values, cohort-level warnings, and flagged gene count.

### Sample size warnings & flagging
- **Warnings column in TSV** — a text column per gene with semicolon-separated flags (e.g., `low_carrier_count;sparse_genotype`). Filterable in Excel/R.
- **Thresholds configurable via AssociationConfig** — fields like `min_cases`, `max_case_control_ratio`, `min_case_carriers` with sensible defaults (200, 20, 10). Exposed as CLI args, mapped in analysis_stages.py. Never hardcoded.
- **Two severity levels** — `critical` for cohort-level issues (n_cases < 200, case:control > 1:20) logged prominently and in summary.txt header; `warn` for per-gene flags (case_carriers < 10) in TSV warnings column only.
- **Keep p-values for flagged genes** — results present, warning column lets users filter. Matches REGENIE behavior — interpretation is the user's responsibility.

### ACAT-O input flexibility
- **Combine whatever is available (2+ p-values)** — ACAT-O runs on any 2+ p-values from active tests.
- **Pass-through for k=1** — if only one test p-value exists, copy it to `acat_o_p_value` with a logged warning. NA for k=0.
- **No ACAT-V** — ACAT-O = Cauchy(burden_p, skat_p). ACAT-V is out of scope for Phase 22.
- **Single FDR on ACAT-O only** — per ARCH-03. FDR/Bonferroni applied to `acat_o_p_value` across all genes. Individual test p-values (fisher, burden, skat) remain uncorrected in output — they exist for diagnostic signal decomposition, not independent hypothesis testing.

### Claude's Discretion
- ACAT-O Cauchy combination implementation details (equal weights vs. adaptive)
- QQ data TSV column format (observed vs expected -log10(p), ranks)
- summary.txt exact formatting and sections
- lambda_gc.tsv column ordering beyond test_name and lambda_gc
- Warning flag string naming conventions

</decisions>

<specifics>
## Specific Ideas

- Column naming must follow existing `{test_name}_{metric}` prefix pattern from engine.py line 251
- Warning thresholds go as AssociationConfig dataclass fields (same pattern as Phase 19 covariate fields)
- ACAT-O FDR is the single correction pass — aligns with Liu et al. 2019 and SAIGE-GENE+ conventions
- REGENIE and SAIGE-GENE both use fixed-schema output — our approach matches industry standard

</specifics>

<deferred>
## Deferred Ideas

- ACAT-V (variant-level Cauchy combination) — Phase 23
- matplotlib QQ plot PNG/SVG generation — Phase 23 (DIAG-04)
- Per-test FDR columns — deliberately excluded; statistically incorrect per ARCH-03

</deferred>

---

*Phase: 22-acat-o-and-diagnostics*
*Context gathered: 2026-02-21*
