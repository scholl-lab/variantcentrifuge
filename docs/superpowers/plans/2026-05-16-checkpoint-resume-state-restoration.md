# Checkpoint Resume State Restoration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans`.
> For parallel work, use `superpowers:subagent-driven-development` with disjoint file
> ownership. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Fix issue #101 so checkpoint resume cannot report success unless required
runtime state is restored or required outputs are produced.

**Architecture:** Make checkpoint resume policy explicit. Completed stages may be
skipped only when downstream state is durable or restorable. Volatile in-memory
DataFrame stages must have their persisted checkpoint completion cleared and then be
re-executed. Requested analysis/output stages fail on missing required runtime state,
but valid zero-row DataFrames write header-only outputs.

**Tech Stack:** Python, pytest, pandas, VariantCentrifuge pipeline runner and checkpoint
system.

---

## Files

- Modify `variantcentrifuge/pipeline_core/stage.py`
  - Record checkpoint failures in `Stage.__call__()`.
  - Add a conservative checkpoint resume policy hook.
- Modify `variantcentrifuge/checkpoint.py`
  - Make `fail_step()` thread-safe.
  - Allow `--resume-from` to target failed/stale/incomplete stages when safe.
  - Add small status/downstream clearing helpers if useful.
- Modify `variantcentrifuge/pipeline_core/runner.py`
  - Replace blind resume pre-marking with policy-aware restore/recompute planning.
  - Clear recompute-stage checkpoint completion persistently so `Stage.__call__()` does
    not skip through its own checkpoint gate.
  - Replace the hard-coded selective-resume prerequisite stage set with policy logic.
- Modify `variantcentrifuge/stages/processing_stages.py`
  - Harden `ParallelCompleteProcessingStage._handle_checkpoint_skip()`.
  - Split parallel VCF completion state from DataFrame chunked-analysis completion.
  - Make `GenotypeReplacementStage` resume skip a clean Phase-11 no-op.
  - Keep `PhenotypeIntegrationStage` artifact checks conditional on real output.
- Modify `variantcentrifuge/stages/analysis_stages.py`
  - Harden `DataFrameLoadingStage._handle_checkpoint_skip()`.
  - Add `ChunkedAnalysisStage._handle_checkpoint_skip()`.
  - Split `dataframe_chunked_analysis_complete` from parallel VCF processing.
  - Fail requested association/gene-burden/statistics stages on missing DataFrame.
  - Write header-only association/gene-burden outputs for valid empty DataFrames.
- Modify `variantcentrifuge/stages/output_stages.py`
  - Fail `TSVOutputStage` on missing DataFrame and missing chunked output.
  - Fail requested Excel generation on missing final TSV.
  - Add checkpoint output-file tracking for TSV/Excel/metadata reports.
- Modify `variantcentrifuge/pipeline.py`
  - Add final output validation before logging successful completion.
- Modify `docs/source/resume_system.md`
  - Document durable/restorable/recompute resume policies.
  - Document safe failed-stage `--resume-from` behavior.
- Add or modify tests under:
  - `tests/unit/pipeline/test_stage.py`
  - `tests/unit/pipeline/test_runner.py`
  - `tests/unit/test_parallel_stage_resume.py`
  - `tests/unit/stages/test_analysis_stages.py`
  - `tests/unit/stages/test_output_stages.py`
  - `tests/unit/test_association_stage.py`
  - `tests/unit/test_gene_burden_regression.py` or a focused gene-burden stage test
  - `tests/test_checkpoint.py`
  - `tests/integration/test_checkpoint_resume_state_restoration.py`

---

## Task 1: Prove Stage Failures Are Not Recorded

**Files:**
- Modify `tests/unit/pipeline/test_stage.py`

- [x] Add `test_stage_failure_marks_checkpoint_failed`.

Test sketch:

```python
class FailingStage(Stage):
    @property
    def name(self):
        return "failing_stage"

    def _process(self, context):
        raise ValueError("boom")


context.checkpoint_state = Mock()

with pytest.raises(ValueError, match="boom"):
    FailingStage()(context)

context.checkpoint_state.start_step.assert_called_once()
context.checkpoint_state.fail_step.assert_called_once_with("failing_stage", "boom")
context.checkpoint_state.complete_step.assert_not_called()
```

- [x] Run and confirm the test fails:

```bash
pytest tests/unit/pipeline/test_stage.py -q
```

---

## Task 2: Record Failures And Lock `fail_step()`

**Files:**
- Modify `variantcentrifuge/pipeline_core/stage.py`
- Modify `variantcentrifuge/checkpoint.py`

- [x] In `Stage.__call__()` exception handling, record failure without hiding the
  original exception:

```python
except Exception as e:
    elapsed = time.time() - start_time
    if context.checkpoint_state:
        try:
            context.checkpoint_state.fail_step(self.name, str(e))
        except Exception as checkpoint_error:
            logger.warning(
                "Failed to record checkpoint failure for stage '%s': %s",
                self.name,
                checkpoint_error,
            )
    logger.error(f"Stage '{self.name}' failed after {elapsed:.1f}s: {e}")
    raise
```

- [x] Wrap `PipelineState.fail_step()` body in `with self._state_lock:` and call
  `self.save()` after releasing the lock, matching the `start_step()` and
  `complete_step()` pattern.
- [x] Audit current process-pool behavior. If checkpoint-relevant stages can execute in a
  child process, add a TODO or guard so parent-side checkpoint state remains authoritative.
- [x] Run:

```bash
pytest tests/unit/pipeline/test_stage.py tests/test_checkpoint.py -q
```

---

## Task 3: Prove Flag Split And Parallel Skip Requirements

**Files:**
- Modify `tests/unit/test_parallel_stage_resume.py`
- Modify `tests/unit/stages/test_analysis_stages.py`

- [x] Add `test_handle_checkpoint_skip_restores_parallel_vcf_completion_only`.

Expected:

```python
result.extracted_tsv == merged_tsv
result.data == merged_tsv
result.config["parallel_vcf_processing_complete"] is True
assert "dataframe_chunked_analysis_complete" not in result.config
```

- [x] Add `test_parallel_skip_does_not_suppress_dataframe_chunking`.

Set up `DataFrameLoadingStage._should_use_chunks()` to return `True` after
`ParallelCompleteProcessingStage._handle_checkpoint_skip()`. Expected:

```python
DataFrameLoadingStage()._handle_checkpoint_skip(context)
assert context.config["use_chunked_processing"] is True
assert context.current_dataframe is None
```

- [x] Add `test_handle_checkpoint_skip_raises_for_missing_merged_tsv`.

Expected:

```python
with pytest.raises(RuntimeError, match="Cannot restore parallel_complete_processing"):
    stage._handle_checkpoint_skip(context)
```

- [x] Run and confirm failures:

```bash
pytest tests/unit/test_parallel_stage_resume.py tests/unit/stages/test_analysis_stages.py -q
```

---

## Task 4: Implement Flag Split And Harden Processing Skips

**Files:**
- Modify `variantcentrifuge/stages/processing_stages.py`
- Modify `variantcentrifuge/stages/analysis_stages.py`

- [x] Replace the parallel VCF completion write:

```python
# ParallelCompleteProcessingStage._process() and _handle_checkpoint_skip()
context.config["parallel_vcf_processing_complete"] = True
context.config.pop("chunked_processing_complete", None)  # only if safe in local tests
```

If removing the legacy key immediately breaks unrelated tests, keep a compatibility read
elsewhere but do not let this key suppress DataFrame chunking.

- [x] In `ParallelCompleteProcessingStage._handle_checkpoint_skip()`, validate before
  mutating context:

```python
if not merged_tsv.exists() or not self._validate_chunk_tsv(merged_tsv):
    raise RuntimeError(
        "Cannot restore parallel_complete_processing from checkpoint: "
        f"merged TSV missing or invalid: {merged_tsv}"
    )

context.extracted_tsv = merged_tsv
context.data = merged_tsv
context.config["parallel_vcf_processing_complete"] = True
context.mark_complete("variant_extraction")
context.mark_complete("snpsift_filtering")
context.mark_complete("field_extraction")
```

- [x] Update `DataFrameLoadingStage` chunking checks to use:

```python
dataframe_chunks_done = context.config.get("dataframe_chunked_analysis_complete", False)
if self._should_use_chunks(context, input_file) and not dataframe_chunks_done:
    context.config["use_chunked_processing"] = True
    return context
```

- [x] Make `GenotypeReplacementStage._handle_checkpoint_skip()` a clean no-op because
  `_process()` is a Phase-11 no-op:

```python
logger.info("Genotype replacement is a Phase-11 no-op; nothing to restore")
return context
```

- [x] For `PhenotypeIntegrationStage._handle_checkpoint_skip()`, raise only when
  phenotype data/config means the stage should have produced a TSV. If no phenotype data
  is configured, return cleanly.
- [x] Run:

```bash
pytest tests/unit/test_parallel_stage_resume.py tests/unit/stages/test_analysis_stages.py -q
```

---

## Task 5: Prove DataFrame And Chunked Analysis Restoration

**Files:**
- Modify `tests/unit/stages/test_analysis_stages.py`

- [x] Add `test_dataframe_loading_skip_raises_when_no_tsv_available`.
- [x] Add `test_dataframe_loading_skip_defers_to_chunked_analysis_for_large_tsv`.
- [x] Add `test_chunked_analysis_skip_restores_dataframe_from_artifact`.
- [x] Add `test_chunked_analysis_skip_restores_header_only_empty_dataframe`.
- [x] Add `test_chunked_analysis_skip_raises_when_artifact_missing`.

Expected for valid chunked restore:

```python
result = ChunkedAnalysisStage()._handle_checkpoint_skip(context)
assert result.current_dataframe is not None
assert result.chunked_analysis_tsv == chunked_tsv
assert result.config["dataframe_chunked_analysis_complete"] is True
```

- [x] Run and confirm failures:

```bash
pytest tests/unit/stages/test_analysis_stages.py -q
```

---

## Task 6: Implement DataFrame And Chunked Restore

**Files:**
- Modify `variantcentrifuge/stages/analysis_stages.py`

- [x] Harden `DataFrameLoadingStage._handle_checkpoint_skip()`:

```python
input_file = self._find_input_file(context)
if not input_file or not Path(input_file).exists():
    raise RuntimeError(
        "Cannot restore dataframe_loading from checkpoint: no TSV input file found"
    )

if self._should_use_chunks(context, input_file) and not context.config.get(
    "dataframe_chunked_analysis_complete", False
):
    context.config["use_chunked_processing"] = True
    return context

df, rename_map = load_optimized_dataframe(...)
context.current_dataframe = df
context.column_rename_map = rename_map
```

- [x] Add `ChunkedAnalysisStage._handle_checkpoint_skip()`:

```python
def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
    chunked_tsv = context.chunked_analysis_tsv
    if not chunked_tsv:
        base = context.workspace.get_intermediate_path("chunked_analysis_results.tsv")
        gz = Path(str(base) + ".gz")
        chunked_tsv = gz if gz.exists() else base

    if not chunked_tsv or not Path(chunked_tsv).exists():
        raise RuntimeError(
            "Cannot restore chunked_analysis from checkpoint: "
            f"chunked output missing: {chunked_tsv}"
        )

    compression = "gzip" if str(chunked_tsv).endswith(".gz") else None
    df, rename_map = load_optimized_dataframe(
        str(chunked_tsv),
        sep="\t",
        compression=compression,
        on_bad_lines="warn",
    )
    context.current_dataframe = df
    context.column_rename_map = rename_map
    context.chunked_analysis_tsv = Path(chunked_tsv)
    context.config["dataframe_chunked_analysis_complete"] = True
    context.config["use_chunked_processing"] = True
    return context
```

- [x] Ensure `ChunkedAnalysisStage._process()` writes
  `chunked_analysis_results.tsv(.gz)` even when all chunks are empty but a header is
  known. If no chunks are processed at all due to missing input, raise instead of warning
  and returning.
- [x] Run:

```bash
pytest tests/unit/stages/test_analysis_stages.py -q
```

---

## Task 7: Prove Missing-State Failures And Empty Outputs

**Files:**
- Modify `tests/unit/test_association_stage.py`
- Modify `tests/unit/test_gene_burden_regression.py` or a focused stage test
- Modify `tests/unit/stages/test_output_stages.py`

- [x] Add `test_association_requested_without_dataframe_raises`.

Make this assert the error occurs before `AssociationEngine.from_names()` is called.

- [x] Add `test_association_empty_dataframe_writes_header_only_output`.
- [x] Add `test_gene_burden_requested_without_dataframe_raises`.
- [x] Add `test_gene_burden_empty_dataframe_writes_header_only_output`.
- [x] Add `test_statistics_enabled_without_dataframe_raises`.
- [x] Add `test_tsv_output_without_dataframe_or_chunked_output_raises`.
- [x] Add `test_excel_requested_without_tsv_raises`.
- [x] Run and confirm failures:

```bash
pytest tests/unit/test_association_stage.py tests/unit/test_gene_burden_regression.py tests/unit/stages/test_output_stages.py -q
```

---

## Task 8: Implement Required Input Contracts And Empty Outputs

**Files:**
- Modify `variantcentrifuge/stages/analysis_stages.py`
- Modify `variantcentrifuge/stages/output_stages.py`

- [x] Move the association DataFrame lookup before engine construction:

```python
df = context.current_dataframe if context.current_dataframe is not None else context.variants_df
if df is None:
    raise RuntimeError("No DataFrame loaded for association analysis")
```

- [x] In association, keep `df.empty` valid and write an empty output:

```python
if df.empty:
    assoc_output, compression = self._resolve_association_output(context)
    empty = pd.DataFrame(columns=self._association_output_columns(test_names))
    empty.to_csv(assoc_output, sep="\t", index=False, compression=compression)
    context.config["association_output"] = str(assoc_output)
    context.association_results = empty
    return context
```

Use a small helper for `_association_output_columns(test_names)` so tests can pin the
schema. At minimum include `gene`, `n_cases`, `n_controls`, `n_variants`,
`primary_test`, `primary_pvalue`, `primary_qvalue`, and `warnings`, plus configured
test p-value/q-value columns where known.

- [x] In gene burden, distinguish `None` from empty:

```python
df = context.current_dataframe
if df is None:
    raise RuntimeError("No DataFrame loaded for gene burden analysis")
if df.empty:
    burden_output, compression = self._resolve_gene_burden_output(context)
    empty = pd.DataFrame(columns=self._gene_burden_output_columns(context))
    empty.to_csv(burden_output, sep="\t", index=False, compression=compression)
    context.config["gene_burden_output"] = str(burden_output)
    context.gene_burden_results = empty
    return context
```

- [x] In statistics, raise only when stats are enabled and `current_dataframe is None`.
  Empty DataFrames should still compute or write an empty/default statistics file if the
  stats engine supports that.
- [x] In `TSVOutputStage`, raise when no DataFrame and no existing
  `context.chunked_analysis_tsv` can be copied:

```python
if df is None and not (context.chunked_analysis_tsv and context.chunked_analysis_tsv.exists()):
    raise RuntimeError("No data to write")
```

- [x] In `ExcelReportStage`, raise when Excel is requested and `context.final_output_path`
  is missing or nonexistent.
- [x] Run focused tests from Task 7.

---

## Task 9: Prove Policy-Aware Plain Resume Clears Persisted State

**Files:**
- Modify `tests/unit/pipeline/test_runner.py`

- [x] Add `test_plain_resume_reruns_volatile_dataframe_stage`.

Expected:

- restorable `dataframe_loading` skip handler is called
- volatile `custom_annotation` is not marked complete before execution
- `checkpoint_state.clear_step_completion("custom_annotation")` is called
- `custom_annotation` executes
- downstream `association_analysis` executes after it

- [x] Add `test_completed_volatile_stage_clears_downstream_completion`.

Checkpoint has completed `custom_annotation`, `inheritance_analysis`, and `tsv_output`.
Expected calls include all three stages. "Ignore in memory" is not sufficient because
`Stage.__call__()` has its own checkpoint skip gate.

- [x] Add `test_restore_failure_aborts_resume_before_downstream_execution`.
- [x] Add `test_stage_call_would_skip_if_runner_does_not_clear_checkpoint` as a guard for
  the dual-skip-gate regression.
- [x] Run and confirm failures:

```bash
pytest tests/unit/pipeline/test_runner.py -q
```

---

## Task 10: Implement Policy-Aware Plain Resume

**Files:**
- Modify `variantcentrifuge/pipeline_core/stage.py`
- Modify `variantcentrifuge/pipeline_core/runner.py`
- Modify stage classes that need explicit policies.

- [x] Add a conservative policy hook:

```python
class Stage(ABC):
    @property
    def checkpoint_resume_policy(self) -> str:
        return "recompute"
```

- [x] Mark restorable stages only after handlers exist:

```python
class DataFrameLoadingStage(Stage):
    @property
    def checkpoint_resume_policy(self) -> str:
        return "restore"
```

Do the same for setup loaders with valid handlers, `ParallelCompleteProcessingStage`,
and `ChunkedAnalysisStage`. Leave volatile DataFrame mutators as the default
`"recompute"`.

- [x] Add helpers in `PipelineRunner`:

```python
def _ordered_stages(self, stages: list[Stage]) -> list[Stage]:
    ordered = []
    for level in self._create_execution_plan(stages):
        ordered.extend(level)
    return ordered

def _clear_stage_and_downstream(
    self, stage_name: str, ordered: list[Stage], checkpoint_state
) -> None:
    names = [stage.name for stage in ordered]
    start = names.index(stage_name)
    for name in names[start:]:
        checkpoint_state.clear_step_completion(name)
```

- [x] Replace the plain `--resume` pre-loop with execution-order logic:

```python
ordered = self._ordered_stages(stages)
for stage in ordered:
    if not context.checkpoint_state.should_skip_step(stage.name):
        continue

    policy = stage.checkpoint_resume_policy
    if policy == "restore":
        if not hasattr(stage, "_handle_checkpoint_skip"):
            raise RuntimeError(f"Stage '{stage.name}' is restorable but has no skip handler")
        context = stage._handle_checkpoint_skip(context)
        context.mark_complete(stage.name)
    elif policy == "skip":
        context.mark_complete(stage.name)
    else:
        self._clear_stage_and_downstream(stage.name, ordered, context.checkpoint_state)
        break
```

- [x] Avoid repeated stale-stage side effects from `should_skip_step()` where possible.
  Prefer status helpers for already-cleaned checkpoint state if you add them.
- [x] Run:

```bash
pytest tests/unit/pipeline/test_runner.py -q
```

---

## Task 11: Fix `--resume-from` Restart Semantics

**Files:**
- Modify `variantcentrifuge/checkpoint.py`
- Modify `variantcentrifuge/pipeline_core/runner.py`
- Modify `tests/test_checkpoint.py`
- Modify `tests/unit/pipeline/test_runner.py`

- [x] Add tests:
  - `validate_resume_from_stage()` accepts an available stage whose checkpoint status is
    `failed`
  - unknown stages still fail
  - unsafe `--resume-from association_analysis` suggests `--resume-from dataframe_loading`
  - the hard-coded prerequisite stage set no longer marks policy-recompute stages
    complete

- [x] Change validation to restart semantics:

```python
if stage_name not in available_stages:
    return False, f"Stage '{stage_name}' is not available in current configuration"
if not self._loaded_from_file:
    return False, "No checkpoint file loaded"
return True, ""
```

Detailed safety belongs in `PipelineRunner`, where the stage graph and policies are
available.

- [x] In `_handle_selective_resume()`, replace the hard-coded `prerequisite_stages` set
  with policy-based restoration of ordered prerequisites. If a prerequisite is volatile
  and is not being re-executed, raise with the nearest safe earlier stage.
- [x] Run:

```bash
pytest tests/test_checkpoint.py tests/unit/pipeline/test_runner.py -q
```

---

## Task 12: Add Output Tracking And Final Validation

**Files:**
- Modify `variantcentrifuge/stages/output_stages.py`
- Modify `variantcentrifuge/pipeline.py`
- Add focused tests in existing output/pipeline test files.

- [x] Add `get_output_files()` to output stages that currently lack it:

```python
class TSVOutputStage(Stage):
    def get_output_files(self, context):
        return [context.final_output_path] if context.final_output_path else []

class ExcelReportStage(Stage):
    def get_output_files(self, context):
        path = context.report_paths.get("excel")
        return [path] if path else []

class MetadataGenerationStage(Stage):
    def get_output_files(self, context):
        path = context.report_paths.get("metadata")
        return [path] if path else []
```

- [x] Add a final validator in `pipeline.py`:

```python
def _validate_requested_outputs(context: PipelineContext) -> None:
    output_file = context.config.get("output_file")
    if output_file not in (None, "stdout", "-"):
        if not context.final_output_path or not context.final_output_path.exists():
            raise RuntimeError("Requested TSV output was not written")

    if context.config.get("perform_association"):
        assoc = context.config.get("association_output")
        if not assoc or not Path(assoc).exists():
            raise RuntimeError("Requested association output was not written")

    if context.config.get("perform_gene_burden"):
        burden = context.config.get("gene_burden_output")
        if not burden or not Path(burden).exists():
            raise RuntimeError("Requested gene burden output was not written")

    if context.config.get("xlsx") or context.config.get("excel"):
        excel = context.report_paths.get("excel")
        if not excel or not Path(excel).exists():
            raise RuntimeError("Requested Excel output was not written")
```

- [x] Call the validator after `runner.run(stages, context)` and before
  `logger.info("Pipeline completed successfully!")`.
- [x] Add tests for missing TSV, association, gene burden, and Excel outputs.

---

## Task 13: Update Resume Documentation

**Files:**
- Modify `docs/source/resume_system.md`

- [x] Document durable file, restorable memory, and volatile memory resume policies.
- [x] State that volatile memory stages are re-run after process restart.
- [x] Document that recompute stages must have checkpoint completion cleared because
  both the runner and `Stage.__call__()` can skip completed checkpoint stages.
- [x] Document that missing required checkpoint artifacts abort resume.
- [x] Update `--resume-from` docs to explain failed-stage recovery and suggested earlier
  safe stages.
- [x] Note that `finalizing` stale cleanup should be treated like `running` or explain
  why it is not.
- [x] State clearly that after this fix, the AGDE reproduction may correctly fail again
  in `association_analysis` until the separate categorical genotype bug is fixed.

---

## Task 14: Add Issue #101 Chunked Regression

**Files:**
- Add `tests/integration/test_checkpoint_resume_state_restoration.py`

- [x] Build a small mock-stage integration that matches the issue shape:
  - durable TSV-producing stage completed
  - `dataframe_loading` completed but deferred to chunked analysis
  - `chunked_analysis` completed and wrote a chunked TSV
  - volatile DataFrame stage completed
  - `association_analysis` failed or stale

- [x] On plain resume:
  - assert the chunked artifact is loaded
  - assert `context.current_dataframe is not None` before association executes
  - assert volatile DataFrame stages re-run
  - assert final TSV and association sidecar exist for success cases

- [x] Add a missing-artifact variant:
  - remove the chunked TSV
  - resume must raise
  - log must not contain `Pipeline completed successfully!`

- [x] Add a zero-row variant:
  - chunked TSV contains only a header
  - resume restores an empty DataFrame
  - TSV and requested sidecars are written as header-only files

---

## Task 15: Run Verification

Focused tests:

```bash
pytest tests/unit/pipeline/test_stage.py -q
pytest tests/unit/test_parallel_stage_resume.py -q
pytest tests/unit/stages/test_analysis_stages.py -q
pytest tests/unit/pipeline/test_runner.py -q
pytest tests/unit/test_association_stage.py -q
pytest tests/unit/test_gene_burden_regression.py -q
pytest tests/unit/stages/test_output_stages.py -q
pytest tests/test_checkpoint.py -q
```

Broader affected suites:

```bash
pytest tests/test_checkpoint.py tests/test_checkpoint_resume.py tests/test_checkpoint_parallel.py tests/test_checkpoint_parallel_resume.py -q
pytest tests/unit/stages tests/unit/pipeline -q
```

Full suite if runtime allows:

```bash
pytest -q
```

---

## Completion Criteria

- [x] Issue #101 reproduction no longer produces a false successful pipeline.
- [x] Plain resume restores chunked DataFrame state or aborts before downstream stages.
- [x] Recompute stages have persisted checkpoint completion cleared before execution.
- [x] Requested analysis/output stages cannot complete with `current_dataframe=None`.
- [x] Valid empty DataFrames produce header-only requested outputs.
- [x] `--resume-from` gives a safe recovery path for failed stages.
- [x] Resume documentation matches implementation.
- [x] Focused and broader checkpoint tests pass.
