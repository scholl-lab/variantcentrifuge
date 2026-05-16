# SnpSift Filter Guards Design

## Problem

GitHub issue #104 reports common and/or benign variants surviving a run that used:

```text
--field-profile dbnsfp5
--preset high_or_lof_or_nmd
--preset rare
--preset not_benign
--split-snpeff-lines before_filters
--force-chunked-processing
```

The concrete BMP2K example, `chr4:78912028 C>T`, has a HIGH SnpEff consequence, gnomAD4.1 joint AF `0.00221454`, and benign ClinVar text. It should pass the first `high_or_lof_or_nmd` term but fail both `rare` and `not_benign`.

The linked issue comment adds more examples from the same output, including HELLS, DEFB1, SMARCE1, PPP5C, COL9A2, HAP1, EFCAB3, and CPLANE1. These are useful regression cases because several are HIGH-impact stop/start-loss records that should still be removed by frequency and benign clinical-significance filters.

## Root Cause

There are two root causes:

1. Three built-in presets are syntactically malformed after field-profile resolution in both `dbnsfp4` and `dbnsfp5`:
   - `high_or_lof_or_nmd`: one extra closing parenthesis.
   - `high_or_pathogenic`: one missing closing parenthesis.
   - `moderate_and_high_prediction`: one missing closing parenthesis.
2. `SnpSift filter` can print parser diagnostics to stderr and still exit with status `0`. VariantCentrifuge currently treats status `0` as success, compresses the temporary VCF, indexes it, and continues.

The failure is most visible when presets are combined through the CLI. The CLI wraps each preset in parentheses and joins them with `&`. With `high_or_lof_or_nmd + rare + not_benign`, SnpSift reports:

```text
line 1:126 mismatched input ')' expecting ...
```

and still emits the BMP2K-like record. With `high_or_pathogenic` and `moderate_and_high_prediction` combinations, SnpSift reports:

```text
missing ')' at '<EOF>'
```

and still exits `0`.

## Affected Surfaces

### Early SnpSift Filtering

`variantcentrifuge.filters.apply_snpsift_filter()` calls `run_command(snpsift_cmd, output_file=tmp_vcf)`. `run_command()` captures stderr but only raises on nonzero return codes. This allows SnpSift parser diagnostics to be ignored.

Affected user paths:

- normal non-chunked SnpSift filtering
- chunked processing subtasks that call `apply_snpsift_filter()`
- `--split-snpeff-lines before_filters` because the split VCF becomes the SnpSift input
- any combination of malformed presets with other presets or custom `--filters`

### Preset Construction

`variantcentrifuge.config.json` contains known malformed presets, and `tests/unit/test_field_profile.py::test_all_presets_have_balanced_parens` explicitly excludes them. This converts a correctness bug into accepted behavior.

### Late and Final TSV Filtering

`filter_tsv_with_expression()` and `filter_dataframe_with_query()` also fail open: on invalid pandas-query expressions they log an error and return or copy unfiltered data. This is not the same SnpSift parser issue, but it is the same safety class: a requested filter fails and unfiltered rows survive. The current test suite expects this behavior for `FinalFilteringStage`.

## Goals

- Fix all built-in malformed presets.
- Fail closed when SnpSift reports parser diagnostics even if it exits `0`.
- Fail closed for invalid late/final pandas filters.
- Add regression coverage using issue #104 examples so common/benign HIGH-impact variants do not leak through `high_or_lof_or_nmd + rare + not_benign`.
- Preserve valid filter behavior and existing preset semantics.

## Non-Goals

- Do not build a full SnpSift expression parser in Python.
- Do not change the biological meaning of `rare`, `not_benign`, `high_or_lof_or_nmd`, `high_or_pathogenic`, or `moderate_and_high_prediction`.
- Do not rewrite the filtering pipeline or change chunking architecture.
- Do not add post-hoc QC as a substitute for failing invalid filters. QC can be a later enhancement, but the primary fix is fail-closed filtering.

## Design

### Preset Syntax Repair

Change only the minimum parentheses needed:

```json
"moderate_and_high_prediction": "((ANN[ANY].IMPACT has 'MODERATE') & ((dbNSFP_REVEL_score >= 0.9) | (dbNSFP_CADD_phred >= 30)))"
"high_or_lof_or_nmd": "(((exists LOF[*].PERC) & (LOF[*].PERC > 0.9)) | ((exists NMD[*].PERC) & (NMD[*].PERC > 0.9)) | (ANN[ANY].IMPACT has 'HIGH'))"
"high_or_pathogenic": "((ANN[ANY].IMPACT has 'HIGH') | (((dbNSFP_clinvar_clnsig =~ '[Pp]athogenic') & !(dbNSFP_clinvar_clnsig =~ '[Cc]onflicting')) | ((ClinVar_CLNSIG =~ '[Pp]athogenic') & !(ClinVar_CLNSIG =~ '[Cc]onflicting'))))"
```

Then remove those preset names from the known-unbalanced exception list. The balanced-parentheses test should cover every built-in preset after field-profile resolution.

### SnpSift stderr Guard

Add a SnpSift-specific stderr classifier in `variantcentrifuge.filters` and use it only for `SnpSift filter` calls.

Fatal diagnostics should include:

- `mismatched input`
- `missing ')'`
- `token recognition error`
- `Error parsing`
- `Cannot parse`
- Java exception headers from SnpSift filter evaluation, such as `Exception in thread "main"`

The guard should raise `RuntimeError` with the filter expression and the first diagnostic lines. This is intentionally not global in `run_command()` because some tools write benign warnings to stderr.

### Command Execution Support

`run_command()` currently hides stderr from callers when the command succeeds. Add an optional `return_result: bool = False` parameter that returns `subprocess.CompletedProcess[str]` instead of stdout/path when callers need stderr. Keep the default behavior unchanged for compatibility.

`apply_snpsift_filter()` should call:

```python
result = run_command(snpsift_cmd, output_file=tmp_vcf, return_result=True)
_raise_for_snpsift_filter_stderr(result.stderr, filter_string)
```

If the guard raises, do not bgzip or index the temporary VCF.

### Late and Final Filter Fail-Closed Behavior

Change `filter_dataframe_with_query()` to raise `ValueError` on invalid query syntax instead of returning the unfiltered DataFrame. `FinalFilteringStage` should allow that exception to fail the pipeline.

Change `filter_tsv_with_expression()` to raise `ValueError` for query failures and file-level filtering errors instead of writing/copying unfiltered input to output. This prevents silent leakage in standalone TSV filtering paths.

### Regression VCF

Add unit tests that create a minimal realistic SnpEff-style VCF with two records:

1. A BMP2K-like common/benign HIGH-impact record:
   - `ANN=T|stop_gained|HIGH|BMP2K|...`
   - `dbNSFP_gnomAD4.1_joint_AF=0.00221454`
   - `dbNSFP_gnomAD4.1_joint_AC=3546`
   - `dbNSFP_clinvar_clnsig=Benign/Likely_benign`
   - `ClinVar_CLNSIG=Benign/Likely_benign`
2. A rare/non-benign HIGH-impact positive-control record that should remain.

Filtering with `high_or_lof_or_nmd + rare + not_benign` under `dbnsfp5` should retain only the positive-control record.

## Alternatives Considered

### Only Fix Preset Parentheses

This fixes the known bad config, but it leaves VariantCentrifuge vulnerable to user-provided malformed `--filters`, future preset regressions, and SnpSift's exit-code behavior. This is insufficient.

### Validate Every Filter by Running SnpSift Before Processing

A preflight command can catch errors before expensive chunk processing, but it still needs the same stderr guard and a representative input VCF. It is useful later, but it is more moving parts than needed for the first robust fix.

### Parse SnpSift Expressions in Python

This would be brittle and hard to keep aligned with SnpSift. Parentheses balance plus real SnpSift stderr validation gives better coverage with less risk.

## Verification Strategy

- Unit tests for malformed SnpSift stderr classification.
- Unit tests for `apply_snpsift_filter()` raising when SnpSift exits `0` with parser diagnostics.
- Config tests that all resolved presets have balanced parentheses for both `dbnsfp4` and `dbnsfp5`.
- Integration-style test with a minimal VCF proving issue #104 records are filtered out.
- Unit tests proving invalid final and TSV filters raise instead of returning unfiltered rows.
- Existing filtering, field-profile, output-stage, and chunked-processing tests must pass.

## Risks

- Some users may rely on invalid final filters failing open. That behavior is unsafe and should be treated as a bug fix.
- SnpSift can emit non-fatal warnings to stderr. The guard must match parser/evaluation diagnostics narrowly instead of failing on any stderr.
- Tests that invoke real SnpSift should be skipped when external tools are unavailable, following existing external-tool test patterns.
