# Row-Level Consequence Qualification For Burden Design

## Problem

GitHub issue #106 reports gene burden rows that do not themselves satisfy the selected
consequence predicate when VariantCentrifuge is run with SnpEff splitting, transcript
selection, and gene burden or association.

The problematic run shape is:

```text
--preset high_or_lof_or_nmd
--preset rare
--preset not_benign
--split-snpeff-lines before_filters
--transcript-file ...
--perform-gene-burden
```

Observed examples include records where a multi-annotation VCF record has a qualifying
`HIGH` annotation for one gene, but the retained row after transcript filtering belongs
to another gene with `ANN[0].IMPACT=MODIFIER` or `LOW`. That retained row can then be
counted in gene burden or association even though the row-level equivalent of the preset
is false:

```text
ANN[0].IMPACT == HIGH
OR LOF[0].PERC > 0.9
OR NMD[0].PERC > 0.9
```

This inflates gene-level counts for neighboring or co-annotated genes.

## Debug Findings

There are two likely contributing defects in the current code.

### Parallel Processing Ignores SnpEff Splitting

The single-threaded pipeline wires `--split-snpeff-lines before_filters` as:

```text
VariantExtractionStage
-> MultiAllelicSplitStage
-> SnpSiftFilterStage
-> TranscriptFilterStage
-> FieldExtractionStage
```

This is encoded in `variantcentrifuge/pipeline.py`, where `MultiAllelicSplitStage` and
`SnpSiftFilterStage(split_before_filtering=True)` are only added in the `threads <= 1`
branch.

For `threads > 1`, the pipeline uses `ParallelCompleteProcessingStage`. Its chunk worker
currently runs:

```text
extract_variants()
-> apply_snpsift_filter()
-> filter_vcf_to_transcripts()
-> extract_fields_bcftools()
```

It does not pass `snpeff_splitting_mode` to worker config and does not call the splitter
at all. Therefore a run using `--threads > 1` can still apply SnpSift to an unsplit
multi-annotation record even when the user requested `before_filters`.

This matches the issue mechanism exactly: record-level `ANN[ANY]` passes because of
gene A, then transcript filtering retains gene B.

The same dataflow bug can survive a partial fix if `before_filters` splitting writes a
split VCF but the `late_filtering` or no-filter fallback still assigns
`filtered_vcf = chunk_vcf`. Those fallback branches must preserve the current working
VCF, not reset to the unsplit extraction output.

### Splitting ANN Does Not Preserve LOF/NMD Provenance

`variantcentrifuge.vcf_eff_one_per_line.split_vcf_effects()` splits `ANN` or `EFF`, but
all non-ANN INFO fields are copied verbatim onto every split line.

For `high_or_lof_or_nmd`, this leaves a second leak path even when `ANN` splitting is
honored: a split row for gene B can retain `LOF=(GENE_A|...|0.95)` or
`NMD=(GENE_A|...|0.95)` from the original record and pass the LOF/NMD branch.

`variantcentrifuge.transcript_filter.filter_vcf_line_by_transcripts()` also preserves
all non-ANN INFO fields verbatim, so transcript filtering does not repair LOF/NMD
provenance after the fact. The fix still belongs in the splitter: in both supported
split modes, splitting should happen before transcript filtering, so transcript
filtering receives already provenance-correct split rows.

## Goals

- Make `--split-snpeff-lines before_filters` mean the same thing in sequential and
  parallel/chunked processing.
- Ensure split rows used for SnpSift filtering carry only LOF/NMD entries that belong to
  the same gene as the retained ANN row.
- Prevent non-qualifying transcript-selected rows from reaching gene burden or
  association when the selected preset qualified only another annotation on the original
  VCF record.
- Add focused regression coverage for the issue #106 cross-gene `HIGH` and LOF/NMD
  cases.
- Preserve existing valid split, transcript-filter, SnpSift, burden, and association
  behavior.

## Non-Goals

- Do not rewrite gene burden or association aggregation.
- Do not build a Python parser for arbitrary SnpSift expressions.
- Do not change the biological meaning of `high_or_lof_or_nmd`, `rare`, or
  `not_benign`.
- Do not make `after_filters` silently claim row-level semantics. Filtering before
  splitting is inherently record-level for `ANN[ANY]`.

## Design

### 1. Honor SnpEff Splitting In Parallel Chunk Workers

Extend `ParallelCompleteProcessingStage._process_chunks_parallel()` worker config with:

```python
"snpeff_splitting_mode": context.config.get("snpeff_splitting_mode"),
```

Then update `_process_single_chunk()` to mirror the sequential dataflow:

```text
extract variants -> chunk_vcf

if snpeff_splitting_mode == "before_filters":
    split chunk_vcf -> chunk_base.split_annotations.vcf.gz
    working_vcf = split_vcf
else:
    working_vcf = chunk_vcf

if late_filtering or no filter expression:
    filtered_vcf = working_vcf
else:
    apply SnpSift filter to working_vcf -> filtered_vcf

if snpeff_splitting_mode == "after_filters":
    split filtered_vcf -> chunk_base.split_annotations.vcf.gz
    filtered_vcf = split_vcf

transcript filter filtered_vcf, if requested
extract fields
```

This fixes the direct `ANN[ANY]` cross-gene leak for parallel runs using
`before_filters`.

### 2. Make Split Rows LOF/NMD-Gene-Aware

Add helper functions in `variantcentrifuge.vcf_eff_one_per_line`:

- parse the gene name and gene ID from the retained `ANN` entry
- parse comma-separated SnpEff LOF/NMD entries of the form
  `(Gene|Gene_ID|num_transcripts|percent_affected)`
- keep only LOF/NMD entries where the LOF/NMD gene or gene ID matches the retained
  ANN gene name or gene ID
- drop the `LOF` or `NMD` INFO key from that split row when no entries match

Scope this pruning to `ANN` rows only. The splitter also supports legacy `EFF`, but EFF
has a different subfield layout; ANN indexes must not be applied to EFF records. Leave
EFF behavior unchanged unless a separate EFF-specific parser is implemented.

The splitter should continue to preserve all other INFO fields. For ANN split rows,
malformed LOF/NMD entries should be dropped from split rows and logged as warnings. That
is fail-closed and avoids reintroducing cross-gene qualification through an unparseable
gene-scoped annotation.

Separating LOF/NMD and re-appending matching values may move those INFO keys later in the
INFO field. This is functionally harmless because INFO fields are keyed by name and VCF
does not require a specific INFO key order.

This makes SnpSift evaluate `LOF[*].PERC` and `NMD[*].PERC` against gene-matched
provenance on split rows.

### 3. Treat Missing Or Late Splitting As Record-Level

No splitting and `after_filters` both let SnpSift evaluate the original multi-annotation
record. In those modes, `ANN[ANY]`, `LOF[*]`, and `NMD[*]` remain record-level
predicates. For gene burden or association plus transcript filtering, these modes cannot
guarantee row-level consequence qualification.

Recommended policy:

- Keep no-split and `after_filters` available for workflows that accept record-level
  filtering.
- When `perform_gene_burden` or `perform_association` is true, transcript filtering is
  requested, `snpeff_splitting_mode != "before_filters"`, and the filter expression
  contains `ANN[ANY]`, `LOF[*]`, or `NMD[*]`, log a clear warning.
- Prefer a hard error only if maintainers want strict fail-closed behavior for
  burden/association correctness. The safer long-term default is probably to require
  `before_filters` for transcript-selected burden/association with consequence presets,
  but that is a compatibility decision.

### 4. Do Not Filter Inside Gene Burden

`perform_gene_burden_analysis()` and `AssociationAnalysisStage` currently aggregate the
DataFrame they receive. That is the right boundary: filtering should be complete before
rows enter the analysis DataFrame.

Adding a burden-local interpretation of `high_or_lof_or_nmd` would duplicate SnpSift
semantics, would not generalize to custom presets, and would leave association with a
parallel implementation burden. Fix the VCF processing contract instead.

## Verification Strategy

- Unit test `split_vcf_effects()` with two ANN entries:
  - `GENE_A` `HIGH` with matching `LOF=(GENE_A|GENEA_ID|1|0.95)`
  - `GENE_B` `MODIFIER`
  - expected: the `GENE_B` split row has no `LOF`/`NMD` from `GENE_A`
- Unit test `split_vcf_effects()` with gene ID matching when gene symbols differ.
- Unit test `ParallelCompleteProcessingStage._process_single_chunk()` for
  `before_filters`:
  - splitter is called after extraction and before `apply_snpsift_filter`
  - `apply_snpsift_filter()` receives the split VCF path
  - transcript filtering receives the filtered VCF path
  - `late_filtering` and no-filter fallbacks keep the split working VCF
- Unit test `ParallelCompleteProcessingStage._process_single_chunk()` for
  `after_filters`:
  - `apply_snpsift_filter()` receives the extracted chunk VCF
  - splitter is called on the filtered VCF
  - transcript filtering/extraction consume the split VCF
- Real-tool integration should run in a fresh temporary directory so stale chunk TSV
  reuse cannot mask the fix.
- Add one real-tool test that calls `_process_single_chunk()` directly. Mocked tests
  verify ordering, but a direct chunk-worker test catches wiring regressions across the
  actual external-tool flow.
- Integration regression with real SnpSift/bcftools/bgzip when available:
  - tiny VCF with one record containing `GENE_A HIGH` and `GENE_B MODIFIER`
  - transcript file selects `GENE_B`
  - run the chunk worker or full pipeline shape with `high_or_lof_or_nmd`
  - expected: `GENE_B` is absent from extracted rows and burden output
- Existing SnpSift filter guard tests, transcript filter tests, chunked processing tests,
  gene burden tests, and association stage tests must continue to pass.

## Risks

- LOF/NMD gene matching uses SnpEff field conventions. If an input VCF uses a
  non-standard LOF/NMD format, the helper must not crash; for ANN split rows it should
  drop malformed LOF/NMD entries and warn.
- No splitting and `after_filters` remain record-level by design. If strict row-level
  semantics are required for all transcript-filtered burden/association consequence runs,
  the pipeline must reject any `snpeff_splitting_mode` other than `before_filters`.
- Parallel worker changes affect performance-sensitive paths. The regression tests
  should mock external tools for unit coverage and keep real-tool coverage tiny.
- `_process_chunks_parallel()` can reuse existing chunk TSVs by filename after checkpoint
  or rerun, and `_validate_chunk_tsv()` currently checks only basic parseability. A user
  who reruns in the same output directory after upgrading could reuse TSV chunks produced
  by the old buggy dataflow. Config-hash validation for chunk outputs is a separate
  hardening task; issue #106 tests should run in fresh directories.
