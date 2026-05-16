# Row-Level Consequence Qualification For Burden Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix issue #106 so transcript-selected split SnpEff rows cannot be counted by
gene burden or association unless the retained row itself satisfies the selected
consequence predicate.

**Architecture:** Correct the VCF processing contract before analysis. First, make
parallel chunk processing honor `snpeff_splitting_mode` exactly like sequential
processing. Second, make the SnpEff splitter preserve LOF/NMD provenance by keeping only
LOF/NMD entries that match the retained ANN gene on each split row. Then add unit and
integration regressions that prove cross-gene `ANN[ANY]` and LOF/NMD leaks are removed.

**Tech Stack:** Python, pytest, pandas, SnpSift, bcftools, bgzip, existing
VariantCentrifuge stages.

---

## Files

- Modify `variantcentrifuge/stages/processing_stages.py`
  - Pass `snpeff_splitting_mode` into parallel worker config.
  - Call `split_snpeff_annotations()` in `_process_single_chunk()` before or after
    SnpSift according to the configured mode.
- Modify `variantcentrifuge/vcf_eff_one_per_line.py`
  - Add LOF/NMD parsing and per-ANN gene matching during split.
- Modify `tests/unit/stages/test_parallel_transcript_filter.py`
  - Add before/after filter split ordering tests for the parallel chunk worker.
- Modify `tests/unit/stages/test_processing_stages_critical.py` or create focused
  splitter tests in `tests/unit/test_vcf_eff_one_per_line.py`
  - Add LOF/NMD provenance tests for split rows.
- Modify or add `tests/integration/test_snpeff_row_level_burden.py`
  - Add a tiny real-tool regression covering issue #106.

---

## Task 1: Prove Parallel Chunk Split Ordering Bugs

**Files:**
- Modify `tests/unit/stages/test_parallel_transcript_filter.py`

- [ ] Add a failing test named
  `test_parallel_chunk_splits_before_snpsift_when_requested`.

Use mocks for `extract_variants`, `split_snpeff_annotations`, `apply_snpsift_filter`,
`filter_vcf_to_transcripts`, and `extract_fields_bcftools`.

Scenario:

```python
config={
    "threads_per_chunk": 1,
    "snpeff_splitting_mode": "before_filters",
    "filters": "(ANN[ANY].IMPACT has 'HIGH')",
    "fields_to_extract": "CHROM POS ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
    "transcript_ids": {"NM_gene_b.1"},
    "vcf_samples": ["S1"],
}
```

Expected assertions:

```python
split_vcf = tmp_path / "sample.chunk_0.split_annotations.vcf.gz"
filtered_vcf = tmp_path / "sample.chunk_0.filtered.vcf.gz"

mock_split.assert_called_once_with(str(chunk_vcf), str(split_vcf))
assert mock_apply_filter.call_args.args[0] == str(split_vcf)
assert mock_filter_transcripts.call_args.args[0] == str(filtered_vcf)
```

- [ ] Add a failing test named
  `test_parallel_chunk_splits_after_snpsift_when_requested`.

Scenario:

```python
config={
    "threads_per_chunk": 1,
    "snpeff_splitting_mode": "after_filters",
    "filters": "(ANN[ANY].IMPACT has 'HIGH')",
    "fields_to_extract": "CHROM POS ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
    "transcript_ids": {"NM_gene_b.1"},
    "vcf_samples": ["S1"],
}
```

Expected assertions:

```python
chunk_vcf = tmp_path / "sample.chunk_0.variants.vcf.gz"
filtered_vcf = tmp_path / "sample.chunk_0.filtered.vcf.gz"
split_vcf = tmp_path / "sample.chunk_0.split_annotations.vcf.gz"

assert mock_apply_filter.call_args.args[0] == str(chunk_vcf)
mock_split.assert_called_once_with(str(filtered_vcf), str(split_vcf))
assert mock_filter_transcripts.call_args.args[0] == str(split_vcf)
```

- [ ] Add a parameterized failing test for the `before_filters` fallback paths:
  `late_filtering=True` and no filter expression. Both should split first and then
  feed the split VCF to transcript filtering and field extraction instead of falling
  back to the unsplit `chunk_vcf`.

- [ ] Run the new test and confirm it fails because the splitter is not called and
  `apply_snpsift_filter()` receives `chunk_vcf`.

Command:

```bash
pytest tests/unit/stages/test_parallel_transcript_filter.py -q
```

---

## Task 2: Implement Parallel Split Ordering

**Files:**
- Modify `variantcentrifuge/stages/processing_stages.py`

- [ ] Add `snpeff_splitting_mode` to `worker_config` in `_process_chunks_parallel()`.

```python
"snpeff_splitting_mode": context.config.get("snpeff_splitting_mode"),
```

- [ ] In `_process_single_chunk()`, create a local `working_vcf` after extraction.

Pseudo-code:

```python
working_vcf = chunk_vcf
split_mode = config.get("snpeff_splitting_mode")

if split_mode == "before_filters":
    split_vcf = intermediate_dir / f"{chunk_base}.split_annotations.vcf.gz"
    split_snpeff_annotations(str(working_vcf), str(split_vcf))
    working_vcf = split_vcf
```

- [ ] Preserve `working_vcf` through fallback paths. In both late-filtering mode and
  the no-filter-expression branch, set `filtered_vcf = working_vcf`, not
  `filtered_vcf = chunk_vcf`.

- [ ] Use `working_vcf` as the input to `apply_snpsift_filter()` when a SnpSift filter
  is present.

- [ ] After SnpSift filtering, split on `after_filters`:

```python
if split_mode == "after_filters":
    split_vcf = intermediate_dir / f"{chunk_base}.split_annotations.vcf.gz"
    split_snpeff_annotations(str(filtered_vcf), str(split_vcf))
    filtered_vcf = split_vcf
```

- [ ] Keep transcript filtering after filtering/splitting, and keep field extraction
  pointed at the most recent `filtered_vcf`.

- [ ] Run:

```bash
pytest tests/unit/stages/test_parallel_transcript_filter.py -q
```

Expected: the new before-filter test and existing transcript-filter tests pass.

---

## Task 3: Prove LOF/NMD Cross-Gene Leakage In Split Rows

**Files:**
- Add or modify `tests/unit/test_vcf_eff_one_per_line.py`

- [ ] Add a failing unit test named
  `test_split_ann_prunes_lof_nmd_to_matching_gene`.

Input VCF line:

```text
chr1  100  .  A  T  .  PASS
ANN=T|splice_donor_variant|HIGH|GENE_A|ENSGA|transcript|NM_a|protein_coding||||||||,
    T|downstream_gene_variant|MODIFIER|GENE_B|ENSGB|transcript|NM_b|protein_coding||||||||
LOF=(GENE_A|ENSGA|1|0.95)
NMD=(GENE_A|ENSGA|1|0.91)
```

Expected:

- the `GENE_A` split row keeps `LOF` and `NMD`
- the `GENE_B` split row does not contain those `LOF` or `NMD` entries
- both rows preserve unrelated INFO fields like `AC=1`

- [ ] Add a second test for gene-ID matching when the LOF/NMD gene symbol differs but
  gene ID matches the retained ANN gene ID.

- [ ] Run the tests and confirm they fail because `other_info_parts` are currently
  copied unchanged to every split row.

Command:

```bash
pytest tests/unit/test_vcf_eff_one_per_line.py -q
```

---

## Task 4: Implement LOF/NMD Provenance-Aware Splitting

**Files:**
- Modify `variantcentrifuge/vcf_eff_one_per_line.py`

- [ ] Add constants for ANN gene positions:

```python
ANN_GENE_NAME_INDEX = 3
ANN_GENE_ID_INDEX = 4
```

- [ ] Add helper `_ann_gene_keys(ann_entry: str) -> set[str]` that returns non-empty
  gene name and gene ID from the retained ANN entry. Guard with `len(parts) > 4`
  before indexing.

- [ ] Add helper `_filter_gene_scoped_info_value(value: str, gene_keys: set[str]) -> str | None`
  for LOF/NMD values:

```python
entries = value.split(",")
for entry in entries:
    stripped = entry.strip()
    parts = stripped.strip("()").split("|")
    if len(parts) >= 2 and (parts[0] in gene_keys or parts[1] in gene_keys):
        keep original entry text
```

Return comma-joined kept entries, or `None` when no well-formed entry matches.

- [ ] In `split_vcf_effects()`, separate `LOF=` and `NMD=` from other INFO fields,
  then, for each retained ANN entry, append only matching LOF/NMD values.

- [ ] Scope LOF/NMD pruning to `field_name == "ANN"`. Legacy `EFF` has a different
  subfield layout, so leave EFF split behavior unchanged unless an EFF-specific parser
  is added later.

- [ ] For malformed LOF/NMD entries on ANN split rows, drop the malformed entry and
  log a warning. This is fail-closed and avoids cross-gene qualification through an
  unparseable gene-scoped annotation.

- [ ] Note in code review that LOF/NMD INFO keys may be re-appended after `ANN`; INFO
  key order is not semantically significant.

- [ ] Preserve current behavior for records without multiple ANN/EFF entries.

- [ ] Run:

```bash
pytest tests/unit/test_vcf_eff_one_per_line.py tests/unit/test_transcript_filter.py -q
```

Expected: all pass.

---

## Task 5: Add Real-Tool Regression For Issue #106

**Files:**
- Add `tests/integration/test_snpeff_row_level_burden.py`

- [ ] Mark the module skipped unless `SnpSift`, `bcftools`, and `bgzip` are available,
  matching `tests/integration/test_snpsift_filter_guards.py`.

- [ ] Use `tmp_path`/fresh temporary directories only. Do not reuse existing output
  directories, because `_process_chunks_parallel()` can reuse stale chunk TSVs by
  filename and config-hash validation is out of scope for this fix.

- [ ] Create a tiny VCF with:

```text
##INFO=<ID=ANN,Number=.,Type=String,...>
##INFO=<ID=LOF,Number=.,Type=String,...>
##INFO=<ID=NMD,Number=.,Type=String,...>
##FORMAT=<ID=GT,Number=1,Type=String,...>
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT CASE CTRL
chr1 100 . A T 500 PASS ANN=...GENE_A HIGH NM_a..., ...GENE_B MODIFIER NM_b...;LOF=(GENE_A|ENSGA|1|0.95);NMD=(GENE_A|ENSGA|1|0.91) GT 0/1 0/0
```

- [ ] Exercise the same core flow as the pipeline:

```python
split_snpeff_annotations(input_vcf, split_vcf)
apply_snpsift_filter(split_vcf, high_or_lof_or_nmd_expr, {"threads": 1}, filtered_vcf)
filter_vcf_to_transcripts(filtered_vcf, transcript_vcf, {"NM_b"})
extract_fields_bcftools(... fields include ANN[0].GENE ANN[0].IMPACT NMD[0].PERC GEN[*].GT ...)
```

Expected: the extracted TSV is header-only or contains no `GENE_B` row.

- [ ] Add a positive-control record where `GENE_B` itself is `HIGH` or has matching
  `NMD=(GENE_B|ENSGB|1|0.95)`. Expected: `GENE_B` remains.

- [ ] Add one real-tool test that calls
  `ParallelCompleteProcessingStage()._process_single_chunk(...)` directly with
  `snpeff_splitting_mode="before_filters"`, the tiny VCF, a one-line BED, transcript
  IDs selecting the non-qualifying transcript, and the `high_or_lof_or_nmd` filter.
  Expected: the chunk TSV has no negative-control `GENE_B` row. This complements the
  mocked ordering tests with actual worker wiring.

- [ ] Run:

```bash
pytest tests/integration/test_snpeff_row_level_burden.py -q
```

---

## Task 6: Add Gene Burden/Association Contract Tests

**Files:**
- Add to `tests/integration/test_snpeff_row_level_burden.py` or a focused stage test.

- [ ] Build a DataFrame from the issue #106 extracted TSV and pass it to
  `perform_gene_burden_analysis()` with one case and one control.

- [ ] Confirm `perform_gene_burden_analysis()` has no minimum-sample guard before
  relying on the 1-case/1-control fixture. If that changes later, increase the fixture
  sample count rather than weakening the assertion.

Expected:

- no burden result for `GENE_B` in the negative-control case
- one burden result for `GENE_B` in the positive-control case

- [ ] For association, keep this lightweight: test that the pre-analysis DataFrame used
  by `AssociationAnalysisStage` excludes `GENE_B` for the negative-control case. Avoid
  invoking regression tests that require large sample sizes unless the test uses Fisher
  only and an existing small-sample path.

- [ ] Run:

```bash
pytest tests/test_gene_burden.py tests/unit/test_association_stage.py -q
```

---

## Task 7: Optional Guard For Ambiguous Record-Level Modes

**Files:**
- Modify `variantcentrifuge/cli.py` or `variantcentrifuge/stages/setup_stages.py`
- Add unit tests in `tests/test_cli.py` or setup-stage tests

- [ ] Decide policy:
  - warning only, for compatibility
  - hard error, for strict burden/association correctness

- [ ] If warning:
  log when all are true:
  - `snpeff_splitting_mode != "before_filters"` (covers both no splitting and
    `after_filters`)
  - transcript filtering is requested
  - `perform_gene_burden` or `perform_association` is true
  - filter expression contains `ANN[ANY]`, `LOF[*]`, or `NMD[*]`

- [ ] If hard error:
  raise `ValueError` with a message telling the user to use
  `--split-snpeff-lines before_filters`.

This can be deferred if maintainers want to keep issue #106 narrowly scoped to
`before_filters`.

---

## Task 8: Verification Sweep

- [ ] Run focused tests:

```bash
pytest \
  tests/unit/test_vcf_eff_one_per_line.py \
  tests/unit/test_transcript_filter.py \
  tests/unit/stages/test_parallel_transcript_filter.py \
  tests/unit/stages/test_processing_stages_critical.py \
  tests/integration/test_snpeff_row_level_burden.py \
  -q
```

- [ ] Run affected analysis tests:

```bash
pytest tests/test_gene_burden.py tests/unit/test_association_stage.py -q
```

- [ ] Run existing SnpSift guard tests:

```bash
pytest tests/unit/test_filters.py tests/integration/test_snpsift_filter_guards.py -q
```

- [ ] If external tools are available, run a small CLI end-to-end command with:

```text
--threads 2
--split-snpeff-lines before_filters
--preset high_or_lof_or_nmd
--transcript-file selecting the non-qualifying transcript
--perform-gene-burden
```

Expected:

- no final row for the non-qualifying gene
- no burden result for the non-qualifying gene
- positive-control qualifying gene remains

---

## Acceptance Criteria

- `--split-snpeff-lines before_filters` is honored in both sequential and parallel
  processing.
- Split rows no longer carry LOF/NMD entries from other genes when SnpEff ANN is split.
- The issue #106 negative-control gene does not appear in extracted TSV, burden output,
  or association input.
- Positive controls still pass when the retained row itself has `HIGH`, matching LOF,
  or matching NMD evidence.
- Existing issue #104 SnpSift fail-closed coverage remains green.
