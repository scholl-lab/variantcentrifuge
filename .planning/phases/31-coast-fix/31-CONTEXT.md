# Phase 31: COAST Fix - Context

**Gathered:** 2026-02-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Fix the COAST rare-variant association test so it produces valid p-values for real-world cohorts where most genes lack all three variant categories (BMV/DMV/PTV), genotype matrices are reliably available to the test, multi-transcript SnpEff effect strings are handled correctly, and classification scoring is configurable via the existing formula engine. Creating new association tests or adding new statistical methods are separate phases.

</domain>

<decisions>
## Implementation Decisions

### Partial-category fallback
- Match the R reference implementation (insitro/AllelicSeries): drop empty categories, proceed with present ones, return None only when ALL categories are empty
- Log WHICH specific categories (BMV/DMV/PTV) are missing per gene at DEBUG level
- Add `coast_status` to the `extra` dict: `"complete"` when all 3 categories present, `"partial"` when 1-2 present (p-value is valid but reduced power)
- Add `coast_missing_categories` to `extra` when status is partial (e.g., `"BMV,PTV"`)
- When all categories empty: existing `coast_skip_reason` pattern with `p_value=None` — unchanged
- Category counts (`coast_n_bmv`, `coast_n_dmv`, `coast_n_ptv`) remain in `extra` for all outcomes

### Classification scoring config
- Reuse the existing scoring formula engine (`scoring/` directory with `variable_assignment_config.json` + `formula_config.json`)
- Classification configs live in `scoring/coast_classification/` (and subdirectories per model)
- Three configs ship: SnpEff-based default, CADD-threshold model, and an empty custom template documenting the format
- ALL classification models use the formula engine — including the current SnpEff-based classification (migrate from hardcoded Python `classify_variants()` to a formula config)
- Formula outputs category codes: 0 (unclassified/excluded), 1 (BMV), 2 (DMV), 3 (PTV)
- CADD thresholds are configurable in the JSON config, shipped with literature defaults
- `--coast-classification <name>` selects the active classification model by directory name

### Auto-field injection
- Pipeline reads the active classification config's `variable_assignment_config.json` to determine required annotation fields
- Missing fields are auto-added to the SnpSift extraction field list at INFO log level
- Auto-injected fields are kept in the final output (TSV/Excel) — full transparency
- If a required annotation field doesn't exist in the VCF at all (e.g., no CADD annotations but CADD classification requested), fail early with a clear error message

### Warning & diagnostics
- Per-gene partial-category messages at DEBUG level (quiet by default, visible with -v/--verbose)
- INFO-level summary after COAST completes: "COAST: 312 genes tested (265 complete, 47 partial, 14 skipped)"
- SnpEff '&'-concatenated effect string resolution logged at DEBUG only
- When `--diagnostics-output` is set, write a `coast_classification.tsv` audit file with per-variant rows: gene, variant_id, original_effect, resolved_effect, assigned_category, score_used

### Claude's Discretion
- Exact formula expressions for SnpEff and CADD classification configs
- Internal implementation of the category-dropping logic (how Aggregator equivalent works in Python)
- coast_classification.tsv column ordering and formatting details
- How the formula engine resolves '&'-concatenated strings (split-and-prioritize vs. regex match)

</decisions>

<specifics>
## Specific Ideas

- R reference (insitro/AllelicSeries) is the authoritative behavioral spec for partial-category handling: drop empty columns from aggregated genotype matrix, run COAST on remaining categories, Cauchy omnibus combines whatever components survive
- The existing scoring formula engine is the correct infrastructure — classification is just a scoring formula that outputs discrete codes instead of continuous scores
- "I want it like pg_dump" spirit: if a classification model requires fields the VCF doesn't have, fail fast and tell the user exactly what's missing — no silent degradation

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 31-coast-fix*
*Context gathered: 2026-02-23*
