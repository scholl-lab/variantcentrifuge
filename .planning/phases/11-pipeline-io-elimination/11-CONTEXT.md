# Phase 11: Pipeline I/O Elimination - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Eliminate two major pipeline bottlenecks: (1) genotype replacement stage (7 hrs on large cohorts) by deferring GT formatting to output time, and (2) SnpSift extractFields (2.7 hrs) by replacing with bcftools query (19x faster, measured). Target: reduce total pipeline time from 10+ hours to under 1 hour. SnpSift filter stays (28 min, needs `na`/`has`/`ANN[ANY]` operators that bcftools lacks).

</domain>

<decisions>
## Implementation Decisions

### Genotype deferral strategy
- **Output stages only**: TSV and Excel stages format genotypes at write time. Analysis stages work with raw per-sample GT columns throughout — no intermediate formatting step
- **Remove _GT_PARSED cache entirely**: The Phase 10 cache column is unused dead code (created but never consumed by any production code) and structurally incompatible with raw SnpSift format (no embedded sample names). Delete `parse_gt_column()` and all `_GT_PARSED` references
- **Reconstruct GT from existing context**: Use `context.vcf_samples` (ordered sample list) + per-sample GT columns (`GEN[i].GT`) to build `"Sample1(0/1);Sample2(1/1)"` format at output time. No new context fields needed — mapping is deterministic from sample index to column name per VCF specification
- **Skip reference genotypes**: Output only non-reference genotypes (skip 0/0, ./.) in the reconstructed GT column, matching current behavior exactly

### bcftools query migration
- **bcftools query + Python post-parse**: Use `bcftools query` for field extraction (19x faster than SnpSift, measured on 100K variants / 5,125 samples). Extract raw ANN as single pipe-separated string, parse sub-fields in Python
- **NOT bcftools +split-vep**: Tested and rejected — drops variants without ANN annotations (233/100K lost), requires fragile header rewrite hack (SnpEff uses "Functional annotations:" vs required "Format:" prefix), and "Annotation" must be renamed to "Consequence"
- **Dynamic field mapping from config**: Build bcftools query format string dynamically from `config.json` `fields_to_extract`. Map SnpSift field names (`ANN[0].GENE`, `GEN[*].GT`) to bcftools equivalents automatically
- **Missing value normalization**: bcftools outputs "." for missing values; normalize to "NA" in the Python post-parse step to match current pipeline expectations
- **Require bcftools >= 1.10**: Stable query format strings and good performance. bcftools is already a hard dependency (BCFToolsPrefilterStage)

### Output compatibility
- **Semantically equivalent output**: Same data, same columns, same values — but allow "." vs "NA" normalization, column reordering, and minor whitespace differences. Not byte-identical
- **Extend golden file infrastructure**: Add full TSV/Excel golden files alongside existing inheritance golden files (Phase 9). Compare before/after refactor on test VCF
- **Update phenotype extraction**: phenotype.py currently parses replaced GT format — must be updated in this phase to work with raw per-sample columns

### Transition & rollback
- **Clean switch, no feature flags**: Replace old extraction/replacement path entirely. Tests + golden files prove equivalence. No --legacy-extraction flags or dead code paths
- **SnpSift stays required**: Still needed for SnpSift filter stage. No conditional dependency logic
- **Old checkpoints incompatible**: Users must re-run pipelines from scratch. Old genotype-replaced TSV checkpoints are not supported
- **Snakemake works as-is**: CLI interface unchanged, Snakemake just invokes the CLI

### Claude's Discretion
- bcftools query format string construction details
- ANN field parsing implementation (string split vs regex)
- Order of operations for the Python post-parse step
- Golden file test scenario selection (which variants, how many)
- How to handle edge cases in ANN parsing (multiple annotations, empty fields)
- Phenotype extraction refactoring approach

</decisions>

<specifics>
## Specific Ideas

### Benchmark results (measured on test VCF)
- **Test dataset**: `testing/gckd_all.GRCh37.annotated.vcf.gz` — 5,125 samples
- **100K variants extraction timing**:
  - SnpSift extractFields: 29m 50s
  - bcftools query: 1m 32s (**19.4x faster**)
  - bcftools +split-vep: 1m 35s (18.9x faster, but drops 233 variants)
- **Column structure difference**:
  - SnpSift: packed GT in 1 column (colon-separated, all 5125 samples)
  - bcftools: per-sample GT in separate tab columns (5125 columns)
- **Per-sample columns are exactly what Phase 11 needs** — eliminates need for genotype replacement entirely

### ANN field parsing
- SnpEff ANN format: `Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos/cDNA.length|CDS.pos/CDS.length|AA.pos/AA.length|Distance|ERRORS/WARNINGS/INFO`
- Multiple annotations comma-separated; take first (index [0]) to match current `ANN[0].*` behavior
- Simple: `ann_string.split(",")[0].split("|")` per row

### Existing infrastructure that supports this change
- `create_sample_columns_from_gt_intelligent()` already handles both replaced and raw SnpSift GT formats
- `context.vcf_samples` already populated by SampleConfigLoadingStage
- `context.column_rename_map` already handles column name sanitization/restoration
- Output stages already restore original column names before writing

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 11-pipeline-io-elimination*
*Context gathered: 2026-02-15*
