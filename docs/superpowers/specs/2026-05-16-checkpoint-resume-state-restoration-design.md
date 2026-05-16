# Checkpoint Resume State Restoration Design

## Problem

GitHub issue #101 reports that `--resume` can report a successful pipeline after a
previous failure even though the resumed process did not restore the in-memory
DataFrame or produce the requested result files.

Observed run shape:

```text
initial run:
  completed through inheritance_analysis
  association_analysis failed

resume:
  stale association_analysis marked for re-execution
  completed upstream stages skipped
  association_analysis completed in ~1s
  tsv_output completed in 0s
  excel_report warned about missing input
  Pipeline completed successfully
```

The output directory then contained metadata, logs, and checkpoint files, but no valid
final `variants.tsv` or association result table.

This is a data integrity bug. A successful exit code must mean that requested analyses
and outputs were actually produced, or that the pipeline intentionally produced valid
empty outputs.

## Debug Findings

### 1. Resume Skip Happens Before Runtime State Is Restored

`PipelineRunner.run()` handles plain `--resume` by iterating over all stages before
execution planning:

```python
for stage in stages:
    if context.checkpoint_state.should_skip_step(stage.name):
        context.mark_complete(stage.name)
        if hasattr(stage, "_handle_checkpoint_skip"):
            context = stage._handle_checkpoint_skip(context)
```

This is only safe for stages whose state can be reconstructed from durable artifacts.
Many completed stages mutate `context.current_dataframe` only in memory and have no
checkpoint skip handler. This illustrative list includes:

- `custom_annotation`
- `inheritance_analysis`
- `variant_analysis`
- `variant_identifier`
- `variant_scoring`
- `final_filtering`
- `pseudonymization`
- `clinvar_pm5`

Skipping these stages after a process restart marks them complete without restoring
their DataFrame mutations. The runner pre-loop is not the only skip gate:
`Stage.__call__()` independently calls `checkpoint_state.should_skip_step()` before
executing a stage. Therefore a policy-aware runner must also clear persisted checkpoint
completion for any stage that must be recomputed; leaving the checkpoint entry completed
will still cause `Stage.__call__()` to skip it.

### 2. The Resume Documentation Promises State Restoration

`docs/source/resume_system.md` says memory-based stages should restore state from disk
when skipped and that `_handle_checkpoint_skip()` is the stage integration point. The
current implementation only partially follows that contract.

The practical rule should be:

- a completed stage may be skipped only if all context state needed by downstream stages
  is restored, or
- the stage must be re-executed, or
- resume must fail with a clear instruction.

### 3. Some Skip Handlers Are Permissive On Missing Artifacts

`ParallelCompleteProcessingStage._handle_checkpoint_skip()` reconstructs the expected
merged TSV path and logs a warning if it is missing or invalid, but still assigns that
path to `context.extracted_tsv`, marks constituent stages complete, and returns.

It also does not restore `context.config["chunked_processing_complete"] = True`, even
though normal `_process()` sets that flag. Without it, `DataFrameLoadingStage` can enter
chunked mode during checkpoint skip and avoid loading `context.current_dataframe`.

The issue #101 log also shows `chunked_analysis` being skipped. That matters because
`DataFrameLoadingStage._handle_checkpoint_skip()` deliberately returns without loading a
DataFrame when the input should be processed by `ChunkedAnalysisStage`. A correct fix for
the observed path must restore `chunked_analysis` output or re-run it; restoring only
`dataframe_loading` is insufficient.

The current `chunked_processing_complete` flag is overloaded:
`ParallelCompleteProcessingStage` sets it after VCF chunk extraction/merge, and
`ChunkedAnalysisStage` sets it after DataFrame chunk analysis. These are different
milestones. Reusing one flag can mask missing `chunked_analysis` state by forcing
`DataFrameLoadingStage` down a full-load path.

`GenotypeReplacementStage` is currently a universal no-op after Phase 11. Its skip
handler warning about a missing `genotype_replaced.tsv` is misleading and should become
a clean no-op skip, not a fatal artifact check. `PhenotypeIntegrationStage` is
conditionally no-op when no phenotype data is loaded; it should only require an output
artifact when it actually produced phenotype-integrated TSV state.

### 4. Requested Analysis Stages Can Complete With No DataFrame

Several stages log an error or warning and return `context` instead of failing when their
required runtime state is missing:

- `AssociationAnalysisStage`: `No DataFrame loaded for association analysis`
- `GeneBurdenAnalysisStage`: `No DataFrame loaded for gene burden analysis`
- `StatisticsGenerationStage`: `No DataFrame loaded for statistics`
- `TSVOutputStage`: `No data to write`
- `ExcelReportStage`: `No input file for Excel generation`

The stage wrapper marks any returned context as successful and completes the checkpoint.
That turns a broken resume into completed `association_analysis`, `tsv_output`, and
`excel_report` entries.

Important distinction:

- `context.current_dataframe is None` means runtime state is missing and should fail for
  stages that require a DataFrame.
- `context.current_dataframe.empty` can be a valid zero-variant result and should produce
  valid header-only outputs where applicable.

Current empty-DataFrame handling is inconsistent with that goal. `AssociationAnalysisStage`
and `GeneBurdenAnalysisStage` return without writing sidecar output when `df.empty`, and
`ChunkedAnalysisStage` only writes `chunked_analysis_results.tsv(.gz)` when at least one
chunk object is returned. Top-level output validation would turn those valid empty runs
into false failures unless the empty-output paths write header-only artifacts.

### 5. Stage Failures Are Not Recorded By The Stage Wrapper

The `checkpoint()` decorator and `CheckpointContext` context manager call
`PipelineState.fail_step()` on exceptions. `Stage.__call__()` starts and completes
checkpoint entries but its exception path only logs and re-raises. The failed stage is
left as `running`, later detected as stale.

Stale cleanup prevents some partial-output reuse, but it loses the original failure
status and error message. The wrapper should record failures directly. Because stages can
run under `ThreadPoolExecutor`, `PipelineState.fail_step()` must be protected by the same
state lock used by `start_step()` and `complete_step()` before `Stage.__call__()` starts
calling it from concurrent stage wrappers.

### 6. `--resume-from` Semantics Contradict The Docs

The docs describe `--resume-from STAGE` as restart behavior:

```text
Restart from specific stage (re-execute stage and all subsequent stages)
```

`PipelineState.validate_resume_from_stage()` currently rejects a target stage unless
that target was completed. That blocks the natural recovery command for a failed stage:

```text
Cannot resume from 'association_analysis':
Stage 'association_analysis' was not completed in previous run.
Cannot resume from an incomplete stage.
```

The direct target may still be unsafe if prerequisite in-memory state cannot be restored.
The CLI should either rewind to a safe durable stage or refuse with the exact safe
`--resume-from` stage to use.

## Goals

- Never report pipeline success when requested final outputs are missing.
- Never mark a requested analysis or output stage complete when its required input state
  is absent.
- Make checkpoint skip behavior explicit: skip only durable or restorable stages.
- Re-execute volatile in-memory stages after resume unless their post-stage state has a
  durable checkpoint artifact.
- Restore chunked DataFrame analysis state when a previous run used chunked analysis.
- Preserve the ability to skip expensive VCF extraction/filtering/chunk processing when
  their outputs are valid.
- Improve `--resume-from` so failed-stage recovery has a safe, understandable path.
- Update tests and docs so future stages follow the same checkpoint contract.

## Non-Goals

- Do not fix the original association categorical error from the first run in this
  issue. That is a separate analysis bug.
- Do not serialize the entire `PipelineContext` object as a quick fix.
- Do not change biological analysis semantics.
- Do not require every no-op disabled stage to write an output file.
- Do not make `GenotypeReplacementStage` resume fail for its intentionally missing
  Phase-11 output file.

## Resume Model

Introduce an explicit checkpoint resume policy for stages.

### Durable File Stage

A stage is safe to skip when all downstream-needed state can be restored from files.

Examples:

- `parallel_complete_processing`
- `variant_extraction`
- `field_extraction`
- `phenotype_integration` when it produced a TSV

Requirements:

- checkpoint output files must exist and validate
- skip handler must restore context paths and flags
- missing expected artifacts must raise a resume error, not warn and continue

### Restorable Memory Stage

A memory stage is safe to skip only if it can rebuild the runtime object from a durable
artifact.

Examples:

- `dataframe_loading` can reload `context.current_dataframe` from the latest TSV
- `chunked_analysis` can reload from `chunked_analysis_tsv` if that file exists

Requirements:

- skip handler must set the exact runtime state downstream stages need
- `None` DataFrame after the handler is a hard failure unless the next stages do not need
  a DataFrame

### Volatile Memory Stage

A stage that mutates in-memory data without a durable post-stage artifact is not safe to
skip after process restart.

Examples:

- `custom_annotation`
- `inheritance_analysis`
- `variant_analysis`
- `variant_identifier`
- `variant_scoring`
- `final_filtering`
- `pseudonymization`
- `clinvar_pm5`

Default behavior:

- plain `--resume` should clear completion for the earliest completed volatile stage
  after the last restored durable state, then re-run that stage and all affected
  downstream stages
- selective `--resume-from` should reject or rewind when the target depends on volatile
  state that cannot be restored

This may re-run cheap DataFrame stages, but it avoids silent data loss.

## Proposed Design

### 1. Add Stage Resume Policy

Add a small API to `Stage`:

```python
class Stage(ABC):
    checkpoint_resume_policy = "recompute"
```

Allowed values:

- `"skip"`: safe to skip without a custom handler because it has no runtime state needed
  downstream
- `"restore"`: safe to skip only by calling `_handle_checkpoint_skip()`
- `"recompute"`: do not skip after restart; re-run even if checkpoint says completed

Prefer conservative defaults. `recompute` is the safest default for stages that do not
declare a policy.

For the first fix, mark known setup and durable/restorable stages explicitly:

- setup loaders with working skip handlers: `"restore"`
- `parallel_complete_processing`: `"restore"`
- `dataframe_loading`: `"restore"`
- `chunked_analysis`: `"restore"` after its handler exists
- disabled/no-op report stages can use `"skip"` only when the feature is disabled
- volatile DataFrame mutators: `"recompute"`

The dual-state invariant is mandatory:

- restored stages: call the skip handler, restore context, mark the stage complete in
  the in-memory `PipelineContext`, and leave the checkpoint entry completed
- recompute stages: do not mark the stage complete in memory, and call
  `clear_step_completion()` for that stage and every downstream stage so
  `Stage.__call__()` cannot skip them through its own checkpoint gate

### 2. Change Plain Resume Planning

Instead of blindly marking every completed checkpoint stage as complete:

1. Load and validate checkpoint state.
2. Build the execution order.
3. Walk completed stages in execution order.
4. For `"restore"` stages:
   - call `_handle_checkpoint_skip()`
   - require it to leave the needed context state valid
   - mark complete only after restoration succeeds
5. For `"skip"` stages:
   - mark complete only when the stage is truly disabled or has no downstream state
6. For `"recompute"` stages:
   - clear that stage and every downstream stage from the checkpoint and save the change
   - stop pre-marking later completed stages
   - let normal execution re-run from that point

This creates the desired behavior for issue #101:

```text
restore parallel_complete_processing
-> restore dataframe_loading if it loaded a DataFrame
-> restore or re-run chunked_analysis if dataframe_loading deferred to chunking
-> re-run custom_annotation/inheritance/association/tsv/excel as needed
```

Because the original AGDE run failed inside `association_analysis`, the first corrected
resume may fail there again with the original categorical genotype bug. That is the
right outcome for issue #101: a real failure must stay visible instead of becoming a
false successful pipeline.

### 2a. Split Chunked Completion State

Replace the overloaded `chunked_processing_complete` meaning with two explicit flags:

```python
context.config["parallel_vcf_processing_complete"] = True
context.config["dataframe_chunked_analysis_complete"] = True
```

`ParallelCompleteProcessingStage` owns `parallel_vcf_processing_complete`.
`ChunkedAnalysisStage` owns `dataframe_chunked_analysis_complete`. For backward
compatibility during migration, read the legacy `chunked_processing_complete` only as a
fallback when neither explicit flag is present; do not let the parallel VCF flag suppress
DataFrame chunking. `DataFrameLoadingStage` should decide whether to defer to
`ChunkedAnalysisStage` based on file size and `dataframe_chunked_analysis_complete`, not
on the parallel VCF processing flag.

### 3. Harden Skip Handlers

`ParallelCompleteProcessingStage._handle_checkpoint_skip()`:

- validate merged TSV before assigning it to context
- raise if the expected merged TSV is absent or invalid
- restore:

```python
context.extracted_tsv = merged_tsv
context.data = merged_tsv
context.config["parallel_vcf_processing_complete"] = True
```

- mark constituent stages complete only after validation succeeds

`DataFrameLoadingStage._handle_checkpoint_skip()`:

- raise if no TSV can be found
- raise if a TSV path exists but cannot be loaded
- if it chooses full loading, require `context.current_dataframe is not None`
- if it intentionally defers to chunked analysis, set `context.config["use_chunked_processing"]`
  and require a later `ChunkedAnalysisStage` restore or recompute before analysis/output
  stages run

`ChunkedAnalysisStage._handle_checkpoint_skip()`:

- find `context.chunked_analysis_tsv` or the default
  `chunked_analysis_results.tsv(.gz)`
- load it into `context.current_dataframe`
- set `context.config["dataframe_chunked_analysis_complete"] = True`
- set `context.chunked_analysis_tsv`
- if the prior chunked run legitimately produced zero rows, load or write a valid
  header-only artifact rather than treating missing data rows as corruption

`GenotypeReplacementStage` and `PhenotypeIntegrationStage`:

- `GenotypeReplacementStage` is a Phase-11 no-op and should skip cleanly without looking
  for `genotype_replaced.tsv`
- if `PhenotypeIntegrationStage` was a true no-op in this configuration, skip cleanly
- if the checkpoint says it produced a file and that file is missing, raise

### 4. Make Missing Required Inputs Fatal

Replace warning-and-return behavior with explicit errors when the feature is requested
and the required input state is missing.

Required hard failures:

- `AssociationAnalysisStage`: `perform_association=True` and no DataFrame
- `GeneBurdenAnalysisStage`: `perform_gene_burden=True` and no DataFrame
- `StatisticsGenerationStage`: stats enabled and no DataFrame
- `TSVOutputStage`: no DataFrame and no valid chunked output file
- `ExcelReportStage`: Excel requested and final TSV path is missing

Keep zero-row DataFrames valid. They should result in valid empty/header-only output
files instead of `None` state.

For `GeneBurdenAnalysisStage`, an empty-but-present DataFrame should still resolve and
store `gene_burden_output`, set `context.gene_burden_results` to an empty DataFrame with
the normal output schema, and write that header. For `AssociationAnalysisStage`, an
empty-but-present DataFrame should resolve and store `association_output`, set
`context.association_results` to an empty DataFrame with the configured association
schema, and write that header. The missing-DataFrame check in association must happen
before `AssociationEngine.from_names()` so the error explains the resume-state problem
instead of surfacing an unrelated backend/dependency error.

### 5. Record Stage Failures In Checkpoint State

Update `Stage.__call__()` exception handling:

```python
except Exception as e:
    if context.checkpoint_state:
        context.checkpoint_state.fail_step(self.name, str(e))
    logger.error(...)
    raise
```

This aligns stage execution with the existing decorator and context-manager checkpoint
behavior.

`PipelineState.fail_step()` should also take `_state_lock`, matching `start_step()` and
`complete_step()`. Process-pool execution remains a separate limitation: checkpoint
mutations made in a child process do not automatically update the parent state, so
checkpoint-relevant stages should either run in the parent/thread executor or return
their checkpoint updates for parent-side persistence.

### 6. Fix `--resume-from`

Make the validation match documented restart semantics:

- the stage must exist in the current pipeline
- the target may be completed, failed, stale, finalizing, or pending
- at least one valid prior durable/restorable checkpoint must exist when dependencies
  are not being re-executed
- the new policy system must replace the hard-coded prerequisite stage name set in
  `_handle_selective_resume()`; two independent mechanisms will drift

When a requested target cannot safely run because prerequisite in-memory state is not
restorable, fail with an actionable message:

```text
Cannot resume from 'association_analysis' safely because required DataFrame state from
'inheritance_analysis' is not durable. Resume from 'dataframe_loading' to reload the TSV
and re-run in-memory analysis stages.
```

An alternative acceptable implementation is automatic rewind:

```text
Requested --resume-from association_analysis, but the nearest safe restart point is
dataframe_loading. Re-running from dataframe_loading.
```

Failing with a suggestion is less surprising for a first fix.

### 7. Add Final Output Validation

After `runner.run(stages, context)` and before logging successful completion:

- if `output_file` is not stdout, require `context.final_output_path` exists
- if `perform_association`, require `context.config["association_output"]` exists after
  association completed
- if `perform_gene_burden`, require `context.config["gene_burden_output"]` exists after
  gene burden completed
- if Excel requested, require the Excel report path exists

This is a last-resort guard. The stage contracts should catch the issue earlier, but the
top-level success log should not rely only on per-stage optimism.

## Acceptance Criteria

- Reproducing issue #101 no longer logs `Pipeline completed successfully` unless the
  final TSV and requested association output exist.
- Plain `--resume` after a failed `association_analysis` restores or re-runs the
  DataFrame-producing path before association executes.
- `association_analysis`, `tsv_output`, and `excel_report` cannot be checkpoint-completed
  when their required inputs are missing.
- Stale or failed stages are recorded as failed in checkpoint state with the original
  error message.
- `--resume-from association_analysis` either works safely or fails with the exact earlier
  stage to use.
- Existing valid checkpoint resumes for expensive extraction/filtering stages still skip
  expensive work when outputs validate.
- Valid zero-variant runs produce header-only final, association, and gene-burden
  outputs where those outputs were requested.
- The overloaded `chunked_processing_complete` flag no longer controls both parallel VCF
  processing and DataFrame chunked analysis.

## Test Strategy

- Unit tests for `Stage.__call__()` failure checkpoint recording.
- Unit tests for `ParallelCompleteProcessingStage._handle_checkpoint_skip()`:
  - valid merged TSV restores path and `parallel_vcf_processing_complete`
  - missing or invalid merged TSV raises
- Unit tests for `DataFrameLoadingStage._handle_checkpoint_skip()`:
  - valid TSV reloads DataFrame
  - missing TSV raises
  - large/deferred TSV sets `use_chunked_processing` without pretending the DataFrame is
    restored
- Unit tests for `ChunkedAnalysisStage._handle_checkpoint_skip()`:
  - valid chunked output reloads DataFrame and flags
  - header-only chunked output restores an empty DataFrame
- Unit tests for required input guards:
  - association requested plus `current_dataframe=None` raises
  - association requested plus an empty DataFrame writes a header-only output
  - gene burden requested plus an empty DataFrame writes a header-only output
  - TSV output with no DataFrame and no chunked result raises
  - Excel requested plus missing final TSV raises
- Runner resume tests:
  - plain resume does not skip volatile completed DataFrame stages
  - completed volatile stage causes its own and downstream checkpoint completion to be
    cleared persistently before `Stage.__call__()` can run
  - restore failures abort resume before downstream stages can run
- Flag-split tests:
  - parallel VCF resume restoration does not suppress DataFrame chunking by itself
  - DataFrame chunked analysis completion is represented by
    `dataframe_chunked_analysis_complete`
- Selective resume tests:
  - incomplete failed target is accepted only when prerequisites are safely restored
  - unsafe target emits the suggested earlier stage
- Integration regression:
  - create a checkpoint with upstream stages completed and `association_analysis`
    stale/failed
  - include the actual issue #101 shape where `dataframe_loading` deferred to
    `chunked_analysis`
  - resume
  - assert `context.current_dataframe` is restored before association executes
  - assert the pipeline exits non-zero if DataFrame restoration is impossible, or
    produces final TSV plus association output if restoration succeeds

## Risks

- Re-running volatile DataFrame stages may add minutes to some resumes, but it is safer
  than trusting stale in-memory state that no longer exists.
- A more scalable future design would make expensive volatile stages such as
  `inheritance_analysis` persist their post-stage DataFrame and become restorable. The
  recompute policy is the minimal correctness fix, not the final performance design.
- Some tests currently expect warning-and-return behavior for missing data. Those tests
  should be updated to distinguish disabled/no-op stages from requested stages with
  broken inputs.
- Splitting `chunked_processing_complete` touches several guard branches in analysis and
  output stages. Keep a legacy fallback during migration, but make tests assert the new
  explicit flag semantics.
- File-producing stages without skip handlers may need follow-up hardening. The first
  fix should cover the issue #101 path and add the policy framework so gaps are visible.
- `cleanup_stale_stages()` currently considers `running` stale but not `finalizing`.
  Include `finalizing` in the hardening pass or document why it remains separate.
