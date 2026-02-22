# Phase 26: Association Testing Documentation - Context

**Gathered:** 2026-02-22
**Status:** Ready for planning

<domain>
## Phase Boundary

Comprehensive documentation for the association testing framework added in Phases 18-24. Includes a new user guide, updates to existing docs, API reference stubs, and a v0.15.0 changelog entry. No code changes — documentation only.

</domain>

<decisions>
## Implementation Decisions

### Guide structure & flow
- Linear tutorial progression with a minimal quick-start section at the top (5-10 lines, copy-paste-ready)
- Workflow-based ordering: Setup (covariates/PCA) → Simple tests (Fisher/Burden) → Advanced tests (SKAT/COAST) → Combining (ACAT-O) → Tuning (weights/diagnostics/config)
- Moderate theory per test: one paragraph explaining the statistical model, when it's appropriate, and its assumptions — enough to choose the right test, but still usage-focused
- Include a comparison table for test selection (trait type, sample size needs, what it detects, computational cost, when to use)

### Example depth & data
- Command + output snippet style: CLI invocation plus 2-5 representative lines of output
- Generic placeholders for data files (`input.vcf`, `covariates.tsv`, `pca.eigenvec`) — user substitutes their own
- Full annotated JSON config example (with inline comments explaining each key, values, and defaults)

### Audience & prerequisites
- Dual audience: bioinformaticians and clinical geneticists. Write at the geneticist level but don't over-explain CLI basics — assume basic command-line competence
- Practical & direct tone: "Run this. You'll see this. If X, do Y." Minimal formality
- Assume VariantCentrifuge is already installed and the reader has a multi-sample VCF — don't cover installation or basic pipeline setup
- Dedicated troubleshooting section at the end (R not found, low sample size warnings, convergence failures, lambda_GC interpretation)

### Existing doc updates
- Cross-reference only in usage.md, cohort_analysis.md, FAQ — brief mention + link, no duplicated content
- README.md: feature bullet + link to the guide. Keep README concise
- API reference stubs: public classes only (AssociationTest, AssociationEngine, AssociationConfig, test classes). Skip internal helpers
- v0.15.0 changelog: feature summary grouped by capability (new tests, covariates/PCA, weights, diagnostics, JSON config). One bullet each

### Claude's Discretion
- Exact section heading wording and markdown formatting
- How much of each output snippet to show
- Whether to use admonitions/callouts for warnings and notes
- Ordering within sections when multiple approaches are equally valid
- FAQ additions beyond what's explicitly discussed

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches for technical genomics documentation.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 26-association-testing-documentation*
*Context gathered: 2026-02-22*
