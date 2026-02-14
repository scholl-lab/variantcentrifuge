# Phase 6: Benchmark Framework - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Performance benchmarking infrastructure with synthetic data generators, benchmark tests at micro/meso/macro levels, regression detection with hard-fail canary assertions, and memory profiling. Local-only tooling for the v0.13.0 optimization sprint — not a permanent CI fixture.

</domain>

<decisions>
## Implementation Decisions

### Synthetic data design
- Derive realistic distributions from real-world VCF data, then fully anonymize: obfuscate all identifiers, randomize/mix genotypes, remove all traceable information
- Generation process runs in a gitignored working directory; generated output is checked for sensitive/traceable content before use
- Never reference or mention source data files in committed code
- Research best practices for synthetic genomic data generation and document the approach
- Support trios for inheritance-specific benchmarks (simple family structure)
- Primary focus is large cohort benchmarking — variantcentrifuge is built for rapid filtering of huge multi-sample VCFs and association testing
- Variant count axis: 100, 1K, 10K, 50K variants
- Sample count axis: 10, 100, 500, 1000 samples
- Reproducible generators (seeded randomness)

### Benchmark coverage
- Three granularity levels: micro (individual hot functions), meso (module-level operations), macro (full pipeline runs)
- Components benchmarked: inheritance analysis, genotype replacement, gene burden, scoring, DataFrame I/O
- Compound het detection gets its own dedicated benchmark suite (separate from general inheritance)
- End-to-end benchmarks cover BOTH classic and stage-based pipelines
- Python code paths only — no external tool (bcftools, SnpSift) benchmarks

### Regression policy
- Uniform 20% regression threshold across all benchmarks
- Hard fail when regression detected — test failure, blocks the change
- Memory budget violations (tracemalloc) are warnings only, not hard failures
- Benchmarks run explicitly only: `pytest -m performance` or `pytest tests/performance/` — not part of normal test suite

### Results & workflow
- Console table output for quick review during development
- JSON file saved to `benchmarks/` at project root (gitignored)
- Simple diff/comparison helper: load two JSON result files, show what got faster/slower
- Ephemeral results — no historical tracking, no timestamped storage
- Local only — no CI integration for benchmarks
- Keep it simple — this is optimization sprint tooling, not permanent infrastructure

### Claude's Discretion
- Exact pytest fixture and parameterization design
- tracemalloc API usage and memory budget thresholds
- Console table formatting
- JSON schema for result files
- Benchmark file organization within tests/performance/
- Specific anonymization techniques for synthetic data generation

</decisions>

<specifics>
## Specific Ideas

- "variantcentrifuge is best for large cohorts, trios is just a side case, its for rapid filtering of huge multisample VCF files" — benchmarks must reflect this primary use case
- Synthetic data derived from real-world distributions ensures benchmarks are meaningful, not just testing random noise
- The comparison helper should be dead simple — load two files, print a diff table

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 06-benchmark-framework*
*Context gathered: 2026-02-14*
