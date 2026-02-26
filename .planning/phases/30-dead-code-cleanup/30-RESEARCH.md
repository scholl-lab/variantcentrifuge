# Phase 30: Dead Code Cleanup - Research

**Researched:** 2026-02-23
**Domain:** Python codebase cleanup — dead stage classes, naming artifacts
**Confidence:** HIGH

## Summary

Phase 30 removes 8 stage classes that exist in source files but are never added to the live
pipeline in `pipeline.py`. These stages are registered in `stage_registry.py` (and some are in
`stages/__init__.py`) but `build_pipeline_stages()` never instantiates them. The pipeline
continues to function correctly after removal because none of these stages produce data that
downstream live stages consume.

Two additional discrepancies exist outside the "8 dead stages" list: `DataSortingStage` and
`ClinVarPM5Stage` are used in `pipeline.py` but are not registered in the registry and not in
`stages/__init__.py`. These are live stages with a registration gap — they should be added to
the registry, not removed.

Beyond stage removal, three naming artifacts need cleaning: the string `"refactored_pipeline"`
(default value in pipeline.py:509 and runner.py:119), the docstring "Enhanced error handling
utilities for the refactored pipeline." in `error_handling.py:2`, and the line "Refactored to
be more modular:" in `analyze_variants.py:7`. Eleven inline comments containing "old pipeline"
also need rewording.

**Primary recommendation:** Remove all 8 dead stage classes, their tests, their registry
entries, and their `__init__.py` references in a single atomic plan. Fix cross-references in live
stages. Then address naming artifacts. Keep `DataSortingStage` and `ClinVarPM5Stage` (live stages
that just lack registration — they are out of scope for removal).

## The 8 Dead Stage Classes

All 8 stages are registered in `stage_registry.py` and several appear in `stages/__init__.py`,
but none are ever passed to `stages.append()` in `pipeline.py`'s `build_pipeline_stages()`.

| Stage Class | File | Registered Name | In `__init__.py` |
|---|---|---|---|
| `ParallelVariantExtractionStage` | `stages/processing_stages.py` | `parallel_variant_extraction` | YES (import + `__all__`) |
| `BCFToolsPrefilterStage` | `stages/processing_stages.py` | `bcftools_prefilter` | YES (import + `__all__`) |
| `TranscriptFilterStage` | `stages/processing_stages.py` | `transcript_filtering` | NO |
| `StreamingDataProcessingStage` | `stages/processing_stages.py` | `streaming_processing` | YES (import ONLY, missing from `__all__`) |
| `PCAComputationStage` | `stages/processing_stages.py` | `pca_computation` | NO |
| `GenotypeFilterStage` | `stages/analysis_stages.py` | `genotype_filtering` | NO |
| `ParallelAnalysisOrchestrator` | `stages/analysis_stages.py` | `analysis_orchestrator` | YES (import + `__all__`) |
| `ParallelReportGenerationStage` | `stages/output_stages.py` | `parallel_reports` | YES (import + `__all__`) |

**Note:** `DataSortingStage` and `ClinVarPM5Stage` appear in `pipeline.py` but NOT in the
registry and NOT in `stages/__init__.py`. They are live but unregistered. Do NOT remove them;
the scope of this phase is removing the 8 dead stages listed above.

## Where register_stage() Calls Live (CLEAN-02)

All `register_stage()` calls for dead stages are in
`variantcentrifuge/stages/stage_registry.py` inside the respective `_register_*_stages()`
functions:

```python
# In _register_processing_stages() — lines ~443-481:
from .processing_stages import (
    BCFToolsPrefilterStage,          # REMOVE import
    ...
    ParallelVariantExtractionStage,  # REMOVE import
    PCAComputationStage,             # REMOVE import
    StreamingDataProcessingStage,    # REMOVE import
    TranscriptFilterStage,           # REMOVE import
    ...
)
register_stage(ParallelVariantExtractionStage, "processing", ["parallel_extraction"], 20.0)  # REMOVE
register_stage(BCFToolsPrefilterStage, "processing", ["bcftools_prefilter", "prefilter"], 15.0)  # REMOVE
register_stage(TranscriptFilterStage, "processing", ["transcript_filter"], 10.0)  # REMOVE
register_stage(StreamingDataProcessingStage, "processing", ["streaming_processing"], 30.0)  # REMOVE
register_stage(PCAComputationStage, "processing", ["pca_computation", "pca"], 15.0)  # REMOVE

# In _register_analysis_stages() — lines ~484-514:
from .analysis_stages import (
    ...
    GenotypeFilterStage,           # REMOVE import
    ParallelAnalysisOrchestrator,  # REMOVE import
    ...
)
register_stage(GenotypeFilterStage, "analysis", ["genotype_filter"], 5.0)  # REMOVE
register_stage(ParallelAnalysisOrchestrator, "analysis", ["analysis_orchestrator"], 45.0)  # REMOVE

# In _register_output_stages() — lines ~517-541:
from .output_stages import (
    ...
    ParallelReportGenerationStage,  # REMOVE import
    ...
)
register_stage(ParallelReportGenerationStage, "output", ["parallel_reports"], 25.0)  # REMOVE
```

## Where imports and `__all__` Entries Live (CLEAN-03)

File: `variantcentrifuge/stages/__init__.py`

Dead stages currently in `__init__.py`:

```python
# These imports must be removed:
from .analysis_stages import (
    ...
    ParallelAnalysisOrchestrator,   # line 17 — REMOVE
    ...
)
from .output_stages import (
    ...
    ParallelReportGenerationStage,  # line 28 — REMOVE
    ...
)
from .processing_stages import (
    BCFToolsPrefilterStage,         # line 34 — REMOVE
    ...
    ParallelVariantExtractionStage, # line 37 — REMOVE
    ...
    StreamingDataProcessingStage,   # line 43 — REMOVE
    ...
)

# These __all__ entries must be removed:
"BCFToolsPrefilterStage",          # line 63
"ParallelAnalysisOrchestrator",    # line 81
"ParallelReportGenerationStage",   # line 82
"ParallelVariantExtractionStage",  # line 83
# NOTE: "StreamingDataProcessingStage" is imported but NOT in __all__ — just remove import
```

Dead stages NOT in `__init__.py` (no change needed there):
- `TranscriptFilterStage` — not imported, not in `__all__`
- `PCAComputationStage` — not imported, not in `__all__`
- `GenotypeFilterStage` — not imported, not in `__all__`

## Where "refactored_pipeline" Appears (CLEAN-04, TD-06)

Exactly 2 occurrences in production code:

| File | Line | Content |
|------|------|---------|
| `variantcentrifuge/pipeline.py` | 509 | `pipeline_version = initial_config.get("pipeline_version", "refactored_pipeline")` |
| `variantcentrifuge/pipeline_core/runner.py` | 119 | `pipeline_version = context.config.get("pipeline_version", "refactored_pipeline")` |

**Change:** Replace the default string `"refactored_pipeline"` with `"pipeline"` in both
locations. The `pipeline_version` key is used for checkpoint compatibility checks
(`checkpoint_state.can_resume(config, pipeline_version)`). Changing the default means new
runs write `"pipeline"` to checkpoint state. Old checkpoints written with `"refactored_pipeline"`
will not resume cleanly — this is acceptable since the requirement explicitly calls for this change.

Note: `variantcentrifuge/cli.py` contains the variable name `refactored_args` (lines 1660-1665).
This is a local variable name, NOT a naming artifact — do NOT rename it.

## Where "old pipeline" and "Refactored" Comments/Docstrings Appear (CLEAN-05, CLEAN-06)

### "Refactored" docstrings (CLEAN-06) — 2 occurrences:

| File | Line | Content |
|------|------|---------|
| `variantcentrifuge/pipeline_core/error_handling.py` | 2 | `Enhanced error handling utilities for the refactored pipeline.` |
| `variantcentrifuge/analyze_variants.py` | 7 | `Refactored to be more modular:` |

**Change for error_handling.py:** Replace `"for the refactored pipeline"` with `"for the pipeline"`.

**Change for analyze_variants.py:** Replace `"Refactored to be more modular:"` with
`"Modular design:"` or simply remove the word "Refactored". The lines below it still accurately
describe the design.

### "old pipeline" inline comments (CLEAN-05) — 11 occurrences in production code:

| File | Line | Comment Text |
|------|------|---|
| `variantcentrifuge/pipeline.py` | 492 | `# Map igv argument to igv_enabled for consistency with old pipeline` |
| `variantcentrifuge/stages/analysis_stages.py` | 1375 | `# Create default stats file path like the old pipeline does` |
| `variantcentrifuge/stages/output_stages.py` | 186 | `# Use VAR_ID to match old pipeline` |
| `variantcentrifuge/stages/output_stages.py` | 194 | `# Generate IDs based on key fields with hash like old pipeline` |
| `variantcentrifuge/stages/output_stages.py` | 244 | `id_column = "VAR_ID"  # Match the old pipeline` |
| `variantcentrifuge/stages/output_stages.py` | 741 | `# Add additional sheets like the old pipeline does` |
| `variantcentrifuge/stages/output_stages.py` | 753 | docstring: `Add additional sheets to the Excel file like the old pipeline does.` |
| `variantcentrifuge/stages/output_stages.py` | 1028 | `# Create TSV metadata file like the old pipeline` |
| `variantcentrifuge/stages/output_stages.py` | 1033 | `# Write metadata in TSV format like the old pipeline` |
| `variantcentrifuge/stages/processing_stages.py` | 1327 | docstring line: `This stage replicates the old pipeline's behavior where each chunk` |
| `variantcentrifuge/stages/processing_stages.py` | 1927 | `# Use the sort function from old pipeline` |

**Strategy:** Remove the phrase "like the old pipeline does/did" or replace with a neutral
description of what the code actually does. For docstrings that describe historical behavior
(e.g., "replicates the old pipeline's behavior"), rewrite to describe current behavior.

Note: The requirements say "~15 occurrences" but the actual count in production code is 11.
There are no "old pipeline" occurrences in test files.

## Cross-References in Live Stages That Must Be Cleaned

When removing dead stages, live stage dependency declarations that reference dead stage names
must also be updated to avoid stale dependency strings:

| Live Stage | File | Dependency Type | Dead Stage Name to Remove |
|---|---|---|---|
| `SnpSiftFilterStage` | `processing_stages.py:722` | `soft_dependencies` | `"parallel_variant_extraction"` |
| `FieldExtractionStage` | `processing_stages.py:790` | `soft_dependencies` | `"transcript_filtering"` |
| `VariantAnalysisStage` | `analysis_stages.py:1511` | `soft_dependencies` | `"genotype_filtering"` |
| `ParallelCompleteProcessingStage` | `processing_stages.py:1447, 1491` | `context.mark_complete()` calls | `"parallel_variant_extraction"` |

The `context.mark_complete("parallel_variant_extraction")` calls in `ParallelCompleteProcessingStage`
exist to satisfy dependency graph compatibility when running parallel mode. When
`ParallelVariantExtractionStage` is removed from the registry, these calls become no-ops
(they mark a stage name that no longer exists in the registry). They can safely be removed or
left in place — removing them is cleaner.

## Are Any "Dead" Stages Actually Referenced in Tests That Would Break?

YES. Six of the 8 dead stages have dedicated test classes that will fail with `ImportError`
after the class is removed. These test classes and their containing test sections must be removed
alongside the stage classes:

| Dead Stage | Test File | Test Class Name |
|---|---|---|
| `BCFToolsPrefilterStage` | `tests/unit/stages/test_processing_stages_critical.py` | `TestBCFToolsPrefilterStage` (~8 test refs) |
| `ParallelVariantExtractionStage` | `tests/unit/stages/test_processing_stages.py` | `TestParallelVariantExtractionStage` |
| `ParallelVariantExtractionStage` | `tests/unit/test_parallel_stage_resume.py` | `TestParallelVariantExtractionStage` |
| `StreamingDataProcessingStage` | `tests/unit/stages/test_processing_stages_critical.py` | `TestStreamingDataProcessingStage` (~7 test refs) |
| `PCAComputationStage` | `tests/unit/test_pca.py` | entire file (~13 test refs) |
| `ParallelAnalysisOrchestrator` | `tests/unit/stages/test_analysis_stages.py` | `TestParallelAnalysisOrchestrator` (~9 test methods) |
| `ParallelReportGenerationStage` | `tests/unit/stages/test_output_stages_simple.py` | `TestParallelReportGenerationStage` (~5 test methods) |
| `TranscriptFilterStage` | (none) | no tests |
| `GenotypeFilterStage` | (none) | no tests |

All these test files run as part of `make test-fast` (none have `@pytest.mark.slow` or
`@pytest.mark.integration` markers). The plan MUST remove the test classes alongside the
stage classes, otherwise `make test-fast` will fail with `ImportError` even though the
stage removal itself is correct.

The file `tests/unit/test_pca.py` appears dedicated entirely to `PCAComputationStage` — the
whole file should be deleted.

## Safe Removal Order

Removing all 8 stages in a single commit is safe because they have no dependencies between them
(removing A doesn't break B). However, the test files need to be updated atomically with the
source changes to avoid `ImportError` at test collection time.

**Recommended approach: single atomic commit** covering:
1. Remove 8 dead stage class definitions from processing/analysis/output stage files
2. Remove their test classes from test files (or delete entire test file if dedicated)
3. Remove their registry entries from `stage_registry.py`
4. Remove their imports and `__all__` entries from `stages/__init__.py`
5. Remove stale dependency references in live stages (cross-references above)
6. Run `make test-fast` to verify

If doing sequentially (for bisect purposes), the order that minimizes cascading failures:
1. `TranscriptFilterStage` and `GenotypeFilterStage` first (no tests, no upstream deps)
2. `BCFToolsPrefilterStage` (has tests, has cross-refs in now-gone TranscriptFilterStage)
3. `ParallelVariantExtractionStage` (has tests, has cross-refs in live SnpSiftFilterStage)
4. `StreamingDataProcessingStage` (has tests, no cross-refs in live stages)
5. `PCAComputationStage` (has tests, its only dependency "bcftools_prefilter" is also dead)
6. `GenotypeFilterStage` (already done — has a cross-ref in VariantAnalysisStage)
7. `ParallelAnalysisOrchestrator` (has tests)
8. `ParallelReportGenerationStage` (has tests)

## Naming Artifact Strategy Summary

| Artifact | Location | Change |
|---|---|---|
| `"refactored_pipeline"` default | `pipeline.py:509` | `"pipeline"` |
| `"refactored_pipeline"` default | `runner.py:119` | `"pipeline"` |
| `"refactored pipeline"` docstring | `error_handling.py:2` | Remove "refactored " |
| `"Refactored to be more modular:"` | `analyze_variants.py:7` | `"Modular design:"` |
| 11x `"old pipeline"` comments | various stage files | Reword or remove |

## Architecture Patterns

### Removal Pattern

When removing a stage class `FooStage`:

```python
# Step 1: Search for ALL string references to the stage name
# Stage name = return value of FooStage.name property
# E.g., if FooStage.name returns "foo_filtering":
grep -rn "FooStage\|foo_filtering" variantcentrifuge/ tests/

# Step 2: Remove from stage file
# Delete the class definition and ALL its methods

# Step 3: Remove from stage_registry.py
# - Remove from the import block in _register_*_stages()
# - Remove the register_stage() call

# Step 4: Remove from stages/__init__.py
# - Remove from import block
# - Remove from __all__ list

# Step 5: Remove from test files
# - Remove test class and all its methods
# - Remove the import from the test file's import section

# Step 6: Remove cross-references in live stages
# - Remove stage name from soft_dependencies / dependencies sets
# - Remove context.mark_complete("stage_name") calls
```

### Anti-Patterns to Avoid

- **Removing class but leaving registry entry**: Registry tries to instantiate non-existent class
  at startup → `ImportError` or `NameError`.
- **Removing class but leaving test import**: Test file fails to collect → `ImportError` in
  pytest.
- **Removing stage name from registry but leaving it in `__init__.py` `__all__`**: Ruff will
  flag unused imports → CI failure.
- **Leaving stale soft_dependency strings**: Does not cause test failure (PipelineRunner ignores
  soft deps that aren't in active stages), but is misleading dead code.
- **Renaming `refactored_args` variable in cli.py**: This is a local variable name, NOT a naming
  artifact. The requirements target comments, docstrings, and configuration strings.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead |
|---|---|---|
| Finding all references to a stage | Manual search | `grep -rn "ClassName\|stage_name" variantcentrifuge/ tests/` |
| Verifying no stale references remain | Trust your memory | Grep the exact class name AND stage name string after removal |

## Common Pitfalls

### Pitfall 1: Removing class without removing test import

**What goes wrong:** `from variantcentrifuge.stages.processing_stages import BCFToolsPrefilterStage`
in a test file causes `ImportError` when pytest tries to collect the file.

**How to avoid:** Always grep for the class name in `tests/` before and after removal.

### Pitfall 2: "bcftools_prefilter" as config key vs stage name

**What goes wrong:** The string `"bcftools_prefilter"` appears in `checkpoint.py:820`,
`cli.py:1440`, and `filters.py:141` — but as a **config parameter key**, not a stage name.
These must NOT be removed when removing `BCFToolsPrefilterStage`.

**How to avoid:** Distinguish between `BCFToolsPrefilterStage` (class name, remove) and
`"bcftools_prefilter"` as a config dict key (keep). The config key controls the bcftools
filter expression applied in variant extraction; the feature still works via `filters.py`
directly, independent of the stage class.

### Pitfall 3: `refactored_args` variable name in cli.py

**What goes wrong:** CLEAN-05/06 requirements target naming artifacts in comments/docstrings.
The variable `refactored_args` in `cli.py:1660` is a local variable name (technical noun),
not a naming artifact. Renaming it would be unrelated churn.

**How to avoid:** Only change strings that appear in comments (`# ... old pipeline ...`),
docstrings (`"""..."""`), and config string defaults (`"refactored_pipeline"`).

### Pitfall 4: Checkpoint incompatibility after CLEAN-04

**What goes wrong:** Users with existing checkpoint files written before this change will have
`pipeline_version: "refactored_pipeline"` stored. After the change, the pipeline defaults to
`"pipeline"`, causing resume checks to fail (version mismatch).

**How to avoid:** This is expected behavior per requirements — the change is intentional. No
compatibility shim is needed. Document in commit message that existing checkpoints from before
this phase will not auto-resume (users must restart).

### Pitfall 5: PCAComputationStage has real tests

**What goes wrong:** `PCAComputationStage` has 13 test references in
`tests/unit/test_pca.py` covering real functionality. Simply deleting the class will break
these tests.

**How to avoid:** Delete `tests/unit/test_pca.py` entirely when removing `PCAComputationStage`.
The stage and its tests are dead code — the underlying `akt pca` tool invocation is not wired
into the pipeline.

## Open Questions

1. **Should DataSortingStage and ClinVarPM5Stage be added to registry in this phase?**
   - They are used in `pipeline.py` but not registered, creating inconsistency
   - The requirements (CLEAN-01 through CLEAN-06) only cover removal of dead stages
   - MEDIUM confidence: Recommend leaving this as a separate small fix, not blocking this phase

2. **Do existing user checkpoint files matter?**
   - The `"refactored_pipeline"` → `"pipeline"` change breaks resume for existing checkpoints
   - Confidence: HIGH — this breakage is intentional per requirements

## Sources

### Primary (HIGH confidence)

- Direct code inspection of all referenced files at the paths listed in this document
- `variantcentrifuge/stages/processing_stages.py` — all Stage class definitions read
- `variantcentrifuge/stages/analysis_stages.py` — Stage class definitions identified via grep
- `variantcentrifuge/stages/output_stages.py` — Stage class definitions identified via grep
- `variantcentrifuge/stages/stage_registry.py` — all `register_stage()` calls read
- `variantcentrifuge/stages/__init__.py` — all imports and `__all__` entries read
- `variantcentrifuge/pipeline.py` — all `stages.append()` calls read (lines 144-325)
- `.planning/REQUIREMENTS.md` — CLEAN-01 through CLEAN-06 requirements
- `.planning/ROADMAP.md` — Phase 30 goals and success criteria
- `.planning/research/PITFALLS.md` — dead stage pitfalls section

## Metadata

**Confidence breakdown:**
- Which 8 stages are dead: HIGH — verified by cross-referencing registry vs pipeline.py
- Register/import locations: HIGH — direct file reading
- Test file coverage: HIGH — grep confirmed
- Cross-references in live stages: HIGH — direct file reading
- "refactored_pipeline" locations: HIGH — grep confirmed (2 occurrences)
- "old pipeline" comment count: HIGH — grep confirmed (11 occurrences in production code)
- Checkpoint breakage consequence: HIGH — code logic read directly

**Research date:** 2026-02-23
**Valid until:** Until next code change in stages/ — 90 days if code is stable
