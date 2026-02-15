# Phase 11: Pipeline I/O Elimination - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Eliminate two major pipeline bottlenecks: genotype replacement stage (7 hrs) and SnpSift extractFields (2.7 hrs). Target: reduce total pipeline time from 10+ hours to under 1 hour on large cohorts. SnpSift filter stays unchanged (28 min, needs `na`/`has`/`ANN[ANY]` operators that bcftools lacks). Issues: #77 (genotype replacement elimination), #76 (parent bottleneck issue).

</domain>

<decisions>
## Implementation Decisions

### Genotype deferral strategy
- **Output stages only**: TSV and Excel stages format genotypes just before writing. Analysis stages work with raw GT codes throughout
- GenotypeReplacementStage is eliminated entirely (not skipped, removed)
- Raw per-sample GT columns (GEN[0].GT, GEN[1].GT, etc.) flow through analysis pipeline as individual DataFrame columns
- Output stages reconstruct human-readable `Sample1(0/1);Sample2(1/1)` format at write time using `context.vcf_samples` (ordered sample list) + per-sample column values
- No new context fields needed — sample-to-column mapping is deterministic from sample index to column name (VCF specification order)

### GT cache removal
- Remove `_GT_PARSED` cache column entirely — it was created in Phase 10 but never consumed by any production code
- Delete `parse_gt_column()` function from `dataframe_optimizer.py`
- Remove `GT_PATTERN` regex from `dataframe_optimizer.py` (keep it in `converter.py` where it's actually used)
- Cache cleanup code in output stages (`startswith("_")` drop) can remain for safety

### bcftools query migration
- Replace SnpSift extractFields with `bcftools query` — measured **19.6x faster** on test VCF (3,747ms vs 191ms for 7,738 variants)
- Use **bcftools query + Python ANN parse** approach (not split-vep, which requires fragile VCF header rewriting for SnpEff compatibility)
- bcftools query extracts: fixed VCF fields (CHROM, POS, REF, ALT, ID, FILTER, QUAL), INFO fields (AC, dbNSFP_*, ClinVar, hgmd, splice scores), raw ANN string, raw NMD string, and per-sample GT columns
- Python post-processing parses ANN sub-fields by pipe-splitting first annotation entry (`split('|')`) — 18ms for 7,738 variants (trivial)
- NMD[0].PERC extracted by pipe-splitting NMD field, taking index 3
- **bcftools is the only supported extraction backend** — no SnpSift fallback (bcftools already a hard requirement for prefiltering)

### Field configuration format
- **New bcftools-native format** in config.json — no legacy SnpSift-style field names (ANN[0].GENE, GEN[*].GT)
- Update config.json `fields_to_extract` to use bcftools query format (%INFO/ANN, [\t%GT], etc.)
- Must be properly documented and tested
- Existing user configs will need migration (breaking change, documented in release notes)

### Missing value handling
- **Update downstream to handle both '.' and 'NA'** as missing values
- bcftools outputs '.' for missing; SnpSift used 'NA'
- All downstream comparison code updated to treat both as missing — more robust long-term

### SnpSift filter retention
- SnpSift filter stage kept as-is (28 min, uses `na`/`has`/`ANN[ANY]` operators)
- Document that SnpSift/Java is still required for the filtering stage
- bcftools replaces extraction only

### Claude's Discretion
- Exact bcftools query format string construction
- Python ANN/NMD parsing implementation details
- Order and approach for downstream missing-value updates
- How to handle edge cases (multi-allelic variants, missing ANN field, phased genotypes)

</decisions>

<specifics>
## Specific Ideas

- Benchmarked on actual test VCF (`testing/gckd_all.GRCh37.annotated.vcf.gz`): 5,125 samples, SnpEff-annotated
- bcftools +split-vep incompatible with SnpEff ANN header format (expects "Format: " prefix, SnpEff uses "Functional annotations: '") — confirmed limitation, drove decision toward bcftools query approach
- Per-sample GT columns from bcftools output align perfectly with genotype deferral — each sample gets its own tab-separated column, no packed format parsing needed
- Estimated full cohort speedup for field extraction: 2.7 hours -> ~8-16 minutes (conservative 10-20x)

</specifics>

<deferred>
## Deferred Ideas

- Investigate replacing SnpSift filter with bcftools filter expressions — could eliminate Java dependency entirely (future phase)
- Auto-detection of old vs new checkpoint format — decided to invalidate old checkpoints instead

</deferred>

## Output & Transition Decisions

### Output compatibility
- **Semantically equivalent** output — same data, columns, values; minor formatting differences acceptable (column order, '.' vs 'NA')
- Not byte-identical — relaxed from original roadmap criterion given format changes
- GT output format **may be improved** (e.g., consistent ordering, phased genotypes preserved) — document any changes
- **Golden output files** for CI validation — generated from current pipeline using existing test fixtures (not full GCKD)

### Transition approach
- **Remove old code entirely** — delete SnpSift extraction code and GenotypeReplacementStage (clean break, like Phase 9's comp_het.py removal)
- **Invalidate old checkpoints** — pipeline detects incompatible checkpoint files and warns user to restart
- **Update Snakemake workflow** — review and update `snakemake/variantcentrifuge.smk` for new extraction approach and removed stages

---

*Phase: 11-pipeline-io-elimination*
*Context gathered: 2026-02-15*
