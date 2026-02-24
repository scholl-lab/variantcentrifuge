# Phase 33: Gene-Level FDR Weighting - Context

**Gathered:** 2026-02-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Implement weighted Benjamini-Hochberg FDR correction with per-gene biological prior weights. Users provide a weight file via CLI; weights are renormalized to mean=1.0 (Genovese 2006 guarantee). Backward-compatible: no weights = identical to current plain BH. IHW is explicitly excluded from this phase.

</domain>

<decisions>
## Implementation Decisions

### Weight file format
- TSV with header row; first two columns are gene identifier and weight
- Additional columns allowed (e.g. source, annotation) — ignored by pipeline, preserved for user reference
- Default weight column name: `weight`; overridable via `--gene-prior-weight-column <name>` CLI flag
- Gene matching: primary match on HGNC gene symbol; optional second column with ENSG IDs used as fallback for unmatched symbols
- Genes absent from file receive weight=1.0 (neutral)

### Diagnostics & transparency
- Per-gene `fdr_weight` column added to association results output — shows the applied (renormalized) weight for each gene
- Detailed summary diagnostics logged when weighted BH runs:
  - Effective number of tests: `sum(w)^2 / sum(w^2)`
  - Weight distribution: min, max, median, mean
  - Coverage: number/percentage of tested genes with weights from file vs defaulting to 1.0
  - Number of genes at default weight=1.0
- Impact comparison: log count of genes that gained/lost significance (crossed FDR threshold) compared to unweighted BH
- Output to BOTH logger (stderr/log file) AND a dedicated diagnostics file (e.g. `fdr_diagnostics.tsv`) with per-gene detail

### IHW decision
- IHW (Independent Hypothesis Weighting) is **not implemented** in this phase — no CLI flag, no stub, no error message
- Rationale: project policy is Python-first clean reimplementation; IHW is Bioconductor-only and would require substantial Python rewrite
- IHW is deferred to backlog as a potential future phase

### Method flag
- No `--gene-prior-method` flag — weighted BH is the only method
- Providing `--gene-prior-weights <file>` activates weighted BH; omitting it produces standard (unweighted) BH

### Edge cases and validation
- **Zero/negative weights:** Error and stop — weights must be strictly positive. Exit with clear error naming offending genes
- **Missing weight coverage:**
  - \>50% of tested genes missing from weight file → WARNING: "Most tested genes have no prior weight — weighting may not improve power"
  - \>80% missing → STRONG WARNING: "Only X% of tested genes have prior weights — consider whether weighting is appropriate"
  - Neither is a blocking error — math remains correct
- **Extreme weights:** Warn only, no cap
  - Warn if max/min weight ratio exceeds 100
  - Warn if any single gene's normalized weight exceeds 10
  - No hard cap — FDR control is mathematically guaranteed regardless of magnitude (Genovese 2006)
- **Single gene (m=1):** Skip weighted BH with info log: "Single gene tested — FDR weighting has no effect, skipping." Mathematically a no-op since renormalization forces weight to 1.0

### Claude's Discretion
- Weight file loader implementation details (pandas vs manual parsing)
- Diagnostics file format (column layout, ordering)
- Exact log message wording and log levels
- How to structure the comparison between weighted and unweighted BH internally (run both vs compute delta)
- Weight renormalization implementation (in-place vs copy)

</decisions>

<specifics>
## Specific Ideas

- Weight renormalization must enforce mean=1.0 (sum to m) — this is the Genovese 2006 FDR guarantee, already noted as an architecture invariant
- The `fdr_weight` output column should show the renormalized weight (what was actually applied), not the raw input weight
- Research references for downstream agents:
  - Genovese, Roeder & Wasserman (2006) — weighted BH theory and sum-to-m constraint
  - Roeder & Wasserman (2009) — robustness to moderate weight misspecification
  - Roquain & van de Wiel (2009) — optimal weighting for FDR control

</specifics>

<deferred>
## Deferred Ideas

- **Python IHW implementation** — Independent Hypothesis Weighting as a future phase. Would require understanding Ignatiadis et al. 2016 algorithm, covariate binning, weight learning, and golden-file validation against R/Bioconductor reference. Substantial scope.

</deferred>

---

*Phase: 33-gene-level-fdr-weighting*
*Context gathered: 2026-02-24*
