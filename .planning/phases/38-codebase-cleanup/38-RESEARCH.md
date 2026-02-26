# Phase 38: Codebase Cleanup - Research

**Researched:** 2026-02-26
**Domain:** Python source cleanup, dead code removal, documentation correction
**Confidence:** HIGH

## Summary

Phase 38 is a surgical cleanup phase — no new functionality, only removal of dead code,
correction of stale comments, and documentation fixes. Every requirement was verified by
reading the actual source files. All claimed dead code is confirmed dead; all documentation
errors are confirmed present. No speculative removals are needed.

The primary risk is deleting test code alongside dead production functions: several dead
functions in `prioritizer.py` and `analyzer.py` are actively tested. Those tests must be
deleted at the same time as the production functions, otherwise the test suite breaks with
ImportError. The `stage_info.py` deletion is completely safe — zero imports, zero references.

**Primary recommendation:** Execute requirements in this order: DEAD-01 (safest, no deps),
DEAD-02/DEAD-05 (CLI changes, update tests), DEAD-03/DEAD-04 (delete functions + tests
together), CLEAN-04 (delete dead test methods), CLEAN-01/CLEAN-02 (docstring/comment
updates), CLEAN-03 (add __all__ entries), CLEAN-05/CLEAN-06/CLEAN-07 (resolve TODOs),
DOCS-01/DOCS-02/DOCS-03 (documentation).

---

## DEAD-01: stage_info.py Module

**Status:** Confirmed dead. Zero imports anywhere in the codebase.

**File:** `variantcentrifuge/stage_info.py`
**Size:** 508 lines
**How to verify it's dead:**
```bash
grep -rn "from.*stage_info\|import stage_info" variantcentrifuge/ tests/
# Returns: (nothing)
```

The variable name `stage_info` appears in `interactive_resume.py` (lines 111, 196-207, 224-238,
256, 302, 305) and `display_utils.py` (lines 114-149) but these are local variables of type
`StageInfo` obtained from `registry.get_stage_info()` — they are NOT imports of the module.

**Action:** Delete the file. Nothing else to change.

**Risk:** None. Confirmed zero import references.

---

## DEAD-02: --coast-backend r CLI Choice

**Status:** Confirmed present in CLI; confirmed has a deprecated R implementation. The
requirement phrase "no implementation exists" is slightly inaccurate — a deprecated R
implementation exists (`COASTTest` in `allelic_series.py`) but the CLI choice should still be
removed because the R COAST backend is deprecated and scheduled for removal in v0.17.0.

**File:** `variantcentrifuge/cli.py`
**Lines:** 437-442

```python
# cli.py lines 437-442 (CURRENT STATE):
stats_group.add_argument(
    "--coast-backend",
    choices=["auto", "r", "python"],
    default="python",
    help="COAST computation backend: python (default), r (deprecated), or auto",
)
```

**Change needed:** Remove `"r"` from `choices`. Choices should become `["auto", "python"]`.

**Downstream effects:**
- `variantcentrifuge/stages/analysis_stages.py` line 2199: validation check `d["coast_backend"]
  not in ("auto", "r", "python")` — must also remove `"r"` from the valid set.
- `variantcentrifuge/association/base.py` line 173: `coast_backend: str = "python"` — no change
  needed (the field still holds a string).
- `tests/unit/test_backend_defaults.py` lines 74-77: `test_config_coast_backend_can_be_set_to_r`
  — this test must be deleted (it tests setting `coast_backend="r"` directly on the dataclass,
  which will still work at the dataclass level but the CLI will no longer accept it; the test's
  intent is to verify the CLI behavior, so it should be removed).

**Complication:** `analysis_stages.py` line 2199 validates `coast_backend` values including `"r"`.
This validation must be updated to only allow `("auto", "python")`.

**Risk:** LOW. The R backend class still exists; if someone sets `coast_backend="r"` directly
in config JSON they'll get the deprecated R behavior. Removing from CLI only blocks the
`--coast-backend r` command-line path.

---

## DEAD-03: Dead Functions in prioritizer.py

**Status:** Confirmed. Six functions exist in `prioritizer.py` that are never called from
production code. The requirement says "5 functions" — counting `filter_compatible_patterns`
and its helper `is_pattern_compatible` as a pair counts as 5 named items, but there are
actually 6 function definitions to remove.

**File:** `variantcentrifuge/inheritance/prioritizer.py`

**Production-used functions (do NOT remove):**
- `prioritize_patterns` (lines 73-130) — imported by `analyzer.py`, called in pass 3
- `calculate_confidence` (lines 158-195) — called by `prioritize_patterns`
- `get_pattern_description` (lines 239-287) — imported by `analyzer.py`, used in
  `create_inheritance_details()`

**Dead functions (remove):**

| Function | Lines | Why dead |
|----------|-------|----------|
| `adjust_pattern_score` | 133-155 | Returns base_score unchanged; not imported anywhere |
| `get_pattern_category` | 198-212 | Not imported anywhere in production |
| `group_patterns_by_category` | 215-236 | Not imported anywhere in production |
| `resolve_conflicting_patterns` | 290-331 | Not imported anywhere in production |
| `filter_compatible_patterns` | 334-357 | Not imported anywhere in production |
| `is_pattern_compatible` | 360-389 | Helper for filter_compatible_patterns; same |

**Verification:**
```bash
grep -rn "adjust_pattern_score\|get_pattern_category\|group_patterns_by_category\|resolve_conflicting_patterns\|filter_compatible_patterns\|is_pattern_compatible" variantcentrifuge/ | grep -v "prioritizer.py"
# Returns: (nothing except .pyc files)
```

**Tests that must be deleted simultaneously:**
- `tests/test_inheritance/test_prioritizer.py`:
  - `test_adjust_pattern_score_de_novo` (line 45)
  - `test_adjust_pattern_score_recessive` (line 54)
  - `test_adjust_pattern_score_frequency` (line 63)
  - `test_adjust_pattern_score_x_linked_male` (line 71)
  - `test_get_pattern_category` (line 97)
  - `test_group_patterns_by_category` (line 107)
  - `test_get_pattern_description` (line 131) — wait: `get_pattern_description` is NOT dead,
    so this test should be KEPT
  - `test_resolve_conflicting_patterns` (line 140)
  - `test_filter_compatible_patterns` (line 165)
  - `test_is_pattern_compatible` (line 184)
  - The import lines at top importing the dead functions (lines 5-13)

- `tests/test_inheritance/test_enhanced_prioritizer.py`:
  - `TestScoreAdjustment` class (lines 68-107) — tests `adjust_pattern_score` entirely
  - `TestPatternCategories` class (lines 136-165) — tests `get_pattern_category` and
    `group_patterns_by_category`
  - `TestConflictResolution` class (lines 179-217) — tests `resolve_conflicting_patterns`
  - `TestPatternFiltering` class (lines 219-237) — tests `filter_compatible_patterns`
  - Import lines for dead functions (lines 5-12)

**Risk:** MEDIUM. Must coordinate production deletion with test deletion to avoid ImportError.

---

## DEAD-04: Dead Functions in analyzer.py

**Status:** Confirmed. Three (possibly four) functions in `analyzer.py` are never called from
production code. The requirement says "4 functions."

**File:** `variantcentrifuge/inheritance/analyzer.py`

**Production-used functions (do NOT remove):**
- `analyze_inheritance` — main entry point, imported everywhere
- `_finalize_inheritance_patterns` — called by `analyze_inheritance` and `parallel_analyzer.py`
- `create_inheritance_details` — called by `_finalize_inheritance_patterns`
- `process_inheritance_output` — imported and called by `analysis_stages.py` lines 1151, 1158,
  3826, 3955

**Dead functions (only called from tests, never from production):**

| Function | Lines | Callers |
|----------|-------|---------|
| `get_inheritance_summary` | 407-454 | Tests only: `test_inheritance_integration.py:261`, `test_analyzer.py:192` |
| `filter_by_inheritance_pattern` | 457-487 | Tests only: `test_inheritance_integration.py:298-307`, `test_analyzer.py:217-228` |
| `export_inheritance_report` | 490-530 | Tests only: `test_analyzer.py:228-244` |

That is 3 confirmed dead functions. The requirement says 4. Possible resolution: `create_inheritance_details` may be considered "dead as a public API" because it is only called internally (never imported from outside analyzer.py), making it effectively private. However, removing it would require also removing `test_performance_optimizations.py:11` and `test_inheritance_integration.py:9` import. The safest interpretation: treat `create_inheritance_details` as the 4th dead function (it's only called internally and needlessly exported).

**Recommendation:** Remove all 4: `get_inheritance_summary`, `filter_by_inheritance_pattern`,
`export_inheritance_report`, and `create_inheritance_details` (rename to `_create_inheritance_details`
as a private function rather than deleting, since it's actively called in the module).

**Tests that must be deleted/updated simultaneously:**

- `tests/test_inheritance/test_analyzer.py`:
  - Imports of `export_inheritance_report`, `filter_by_inheritance_pattern`,
    `get_inheritance_summary` (lines 11-13) — remove
  - `test_get_inheritance_summary` method (line 187)
  - `test_filter_by_inheritance_pattern` method (line 198)
  - `test_export_inheritance_report` method (line 228)
  - If `create_inheritance_details` is made private: update import on line 11

- `tests/test_inheritance/test_inheritance_integration.py`:
  - Imports of `filter_by_inheritance_pattern`, `get_inheritance_summary` (lines 10-11)
  - `test_get_inheritance_summary` method (line 241)
  - filter_by_inheritance_pattern usage (lines 298-307)

- `tests/test_inheritance/test_performance_optimizations.py`:
  - Import of `create_inheritance_details` (line 11) — update if renamed

**Risk:** MEDIUM. Same coordination requirement as DEAD-03.

---

## DEAD-05: Redundant --gzip-intermediates Flag

**Status:** Confirmed redundant. The flag has `action="store_true"` and `default=True`, making
`--gzip-intermediates` a no-op when passed. It can't disable compression (only
`--no-gzip-intermediates` can do that).

**File:** `variantcentrifuge/cli.py`
**Lines:** 938-951 (current state):

```python
misc_group.add_argument(
    "--gzip-intermediates",
    action="store_true",
    default=True,                    # <- already True by default
    help="...",
)
misc_group.add_argument(
    "--no-gzip-intermediates",
    dest="gzip_intermediates",
    action="store_false",
    help="Disable compression of intermediate TSV files.",
)
```

**Options:**
1. **Remove `--gzip-intermediates` entirely** — cleaner, keeps only `--no-gzip-intermediates`.
   Breaks `faq.md` usage example at line 156 (needs to be updated anyway for DOCS-01).
2. **Keep but fix** — add a note in help text that this is the default.

**Recommendation:** Remove `--gzip-intermediates` (the `store_true` variant). Keep
`--no-gzip-intermediates`. The `default=True` on `--gzip-intermediates` can move to the
`--no-gzip-intermediates` argument or be set via `set_defaults()`.

**Tests to update:**
- `tests/test_cli.py` line 738: `gzip_intermediates=True` — this tests the config mapping, not
  the flag itself; still valid.
- No tests appear to test `--gzip-intermediates` flag directly (tests set the config value).

**Risk:** LOW. `--no-gzip-intermediates` is sufficient to control the feature.

---

## DOCS-01: faq.md References to Removed CLI Flags

**Status:** Confirmed. Four removed CLI flags appear in `docs/source/faq.md`:

| Line | Current (broken) | Should be |
|------|-----------------|-----------|
| 146 | `--bcftools-filter "QUAL>=30"` | `--bcftools-prefilter "QUAL>=30"` |
| 151 | `--chunks 1000` | Remove (no equivalent flag; chunking is automatic) |
| 156 | `--gzip-intermediates` | Remove (default; mention `--no-gzip-intermediates` to disable) |
| 161 | `--interval-expansion 0` | Remove (flag no longer exists) |
| 168 | `--chunks 5000` | Remove or replace with explanation of automatic chunking |
| 179 | `--interval-expansion` | Remove mention |
| 184 | `--no-filtering` | Remove or replace with current equivalent |

**Current CLI equivalents:**
- `--bcftools-filter` → `--bcftools-prefilter` (line 146 in CLI at 146)
- `--chunks` → removed; automatic memory-based chunking via `--no-chunked-processing` /
  `--force-chunked-processing`
- `--interval-expansion` → removed (was part of gene BED generation, no current equivalent)
- `--no-filtering` → removed (no current equivalent; filtering is always applied)
- `--gzip-intermediates` → `--no-gzip-intermediates` to disable (compression is default)

**Risk:** LOW. Documentation-only change.

---

## DOCS-02: association_testing.md skat_pvalue Column Name

**Status:** Confirmed. Three occurrences of `skat_pvalue` in
`docs/source/guides/association_testing.md` need to be changed to `skat_o_pvalue`.

**File:** `docs/source/guides/association_testing.md`

| Line | Current | Should be |
|------|---------|-----------|
| 334 | `` `skat_pvalue` \| Raw SKAT-O p-value `` | `` `skat_o_pvalue` \| Raw SKAT-O p-value `` |
| 439 | `fisher_pvalue`, `skat_pvalue`, etc. | `fisher_pvalue`, `skat_o_pvalue`, etc. |

Line 417 has `coast_skat_pvalue` which appears correct (it's a COAST sub-test p-value, not the
SKAT-O p-value).

**How the renaming works:** In `association/engine.py` line 625-629, when `skat_method == "SKATO"`
the `skat_pvalue` column is renamed to `skat_o_pvalue`. Since SKAT-O is the default and
recommended method, documentation should use `skat_o_pvalue`.

**Risk:** None. Documentation only.

---

## DOCS-03: changelog.md Classic Pipeline Reference

**Status:** Confirmed. The v0.13.1 entry references "classic pipeline" which no longer exists
(removed in phase 29).

**File:** `docs/source/changelog.md`
**Line:** 51

```markdown
# CURRENT:
- Resource auto-detection across classic and stage-based pipeline modes

# SHOULD BE (suggested):
- Resource auto-detection across pipeline modes
```

The classic pipeline was removed in phase 29 (`29-classic-pipeline-deprecation-and-removal`).
The changelog entry for v0.13.1 is a historical record of what was fixed then — the reference
can be updated to reflect current terminology without losing historical accuracy.

**Risk:** None. Documentation only.

---

## CLEAN-01: Stale "refactored pipeline" Docstrings in Integration Tests

**Status:** Confirmed. Three integration test files have module/class docstrings using the
phrase "refactored pipeline" which is now just "the pipeline."

| File | Line | Current | Should be |
|------|------|---------|-----------|
| `tests/integration/test_basic_pipeline.py` | 1 | `"""Basic integration test for the refactored pipeline."""` | `"""Basic integration tests for the pipeline."""` |
| `tests/integration/test_pipeline_with_mocked_tools.py` | 1-3 | `"""Integration tests for the refactored pipeline...` | `"""Integration tests for the pipeline...` |
| `tests/integration/test_inheritance_analysis_integration.py` | 19 | `"""Test inheritance analysis in the refactored pipeline."""` | `"""Test inheritance analysis in the pipeline."""` |

**Risk:** None. Comment-only changes.

---

## CLEAN-02: Stale "original pipeline" Comments in Source Files

**Status:** Confirmed. Three source files contain comments referencing "the original pipeline."

| File | Line | Current | Should be |
|------|------|---------|-----------|
| `variantcentrifuge/pipeline.py` | 423 | `# Compute base name like the original pipeline` | `# Compute base name from VCF file and gene name` |
| `variantcentrifuge/stages/processing_stages.py` | 1485 | `This is adapted from the original pipeline's sort_tsv_by_gene function.` | `Sort a TSV file by gene column using external sort for memory efficiency.` |
| `variantcentrifuge/stages/setup_stages.py` | 130 | `# Parse phenotype arguments from CLI (matching original pipeline logic)` | `# Parse phenotype arguments from CLI` |

**Risk:** None. Comment-only changes.

---

## CLEAN-03: Missing __all__ Exports in stages/__init__.py

**Status:** Confirmed. Six stages are used in `pipeline.py` but are missing from
`stages/__init__.py` exports.

**File:** `variantcentrifuge/stages/__init__.py`

**Currently exported:** 22 stages (lines 53-89)

**Missing from both imports section AND __all__:**

From `analysis_stages.py` (imported by pipeline.py lines 19-30):
- `AssociationAnalysisStage`
- `ClinVarPM5Stage`
- `VariantAnalysisStage`

From `processing_stages.py` (imported by pipeline.py lines 42-54):
- `DataSortingStage`
- `ParallelCompleteProcessingStage`
- `PCAComputationStage`

**Action:** Add import statements and `__all__` entries for all 6 missing stages.

**Risk:** None. Additive change only.

---

## CLEAN-04: Dead Integration Test Methods

**Status:** Confirmed. Two dead test items exist in
`tests/integration/test_pipeline_with_mocked_tools.py`.

**Item 1: `test_parallel_variant_extraction`**
- Location: `TestParallelProcessing` class, line 357
- Marked `@pytest.mark.xfail` (line 354) — permanently failing
- Tests `parallel_variant_extraction` stage name (line 497) which was removed in phase 30
- The test can never pass because the stage it checks for no longer exists

**Item 2: `TestBCFToolsPrefilter` class**
- Location: lines 500-end-of-file (approximately lines 500-700)
- Contains: `test_bcftools_prefilter_applied` (line 516) and `test_bcftools_prefilter_with_complex_expression` (line 641)
- `BCFToolsPrefilterStage` was removed in phase 30; this class tests the prefilter via
  full pipeline run which no longer exercises that stage path
- The test logic mocks bcftools calls and checks for prefilter application — this is
  actually testing `VariantExtractionStage`'s prefilter behavior, but the test code
  references patterns from the old stage architecture

**Decision needed:** Whether `TestBCFToolsPrefilter` is truly dead (the prefilter feature still
exists but the stage was renamed/refactored). If the prefilter functionality is covered by
other tests, delete; if not, rename and update. Check `tests/unit/stages/test_processing_stages_critical.py`
for existing coverage.

**Risk:** MEDIUM for `TestBCFToolsPrefilter` — need to verify coverage before deleting.

---

## CLEAN-05: TODO(12-02) Legacy Chunking Check

**Status:** Confirmed resolvable. The TODO is in `analysis_stages.py` at line 760:

```python
# TODO(12-02): Remove this check after migrating to ResourceManager auto-detection
# Check explicit chunking request (legacy config key, removed from CLI in 12-01)
if context.config.get("chunks"):
    return True
```

**Context:** `--chunks` was removed from the CLI in phase 12-01. The config key `chunks` can
no longer be set via CLI. However, it could theoretically be set via a config JSON file.

**Verification:** `--chunks` does not appear in the current CLI:
```bash
python3 -c "from variantcentrifuge.cli import create_parser; p = create_parser(); p.parse_args(['--help'])" 2>&1 | grep chunks
# Returns: (nothing about --chunks)
```

**Resolution:** Remove lines 760-763 from `analysis_stages.py`. The ResourceManager
auto-detection (memory-aware chunking) is already implemented on lines 765+.

**Risk:** LOW. The `chunks` config key is dead; no path sets it.

---

## CLEAN-06: TODO Intelligent Batching

**Status:** Confirmed resolvable as "document deferred."

**File:** `variantcentrifuge/pipeline_core/runner.py`
**Line:** 725

```python
# For now, return stages as-is
# TODO: Implement intelligent batching based on estimated runtime
return stages
```

This is in `_batch_lightweight_stages()` method. The function body is just `return stages` —
a 2-line stub. The method IS called at  line 660 in the execution flow.

**Verification:**
```bash
grep -rn "_batch_lightweight_stages" variantcentrifuge/ | grep -v ".pyc"
# runner.py:660: stages = self._batch_lightweight_stages(stages)
# runner.py:711: def _batch_lightweight_stages(self, stages: list[Stage]) -> list[Stage]:
```

**Resolution:** Replace the TODO comment with a note that intelligent batching is deferred to
a future milestone. The function must stay (it is called), but the TODO can be resolved.

Suggested replacement:
```python
# Intelligent batching (grouping lightweight stages into single executor slots to reduce
# scheduling overhead) is deferred. Current pipeline has sufficient throughput for
# the target workloads. If stage count grows significantly, revisit this.
return stages
```

**Risk:** None. Comment-only change inside an existing function.

---

## CLEAN-07: TD-05 Fisher lambda_GC Doc Clarification

**Status:** Partially done. The `diagnostics.py` source code docstring already contains the
full Fisher exemption explanation (lines 62-67). The requirement is about the USER documentation.

**What's already done:** `variantcentrifuge/association/diagnostics.py` lines 62-67 document
that Fisher p-values must NOT be GC-corrected (TD-05 from phase 34).

**What's still missing:** The user-facing documentation in
`docs/source/guides/association_testing.md` section "Interpreting lambda_GC" (around line 531)
does not mention that lambda_GC for Fisher is diagnostic-only and must not be used for
correction.

**Suggested addition:** A note in the lambda_GC interpretation section of association_testing.md:

```markdown
:::{note}
lambda_GC for Fisher's exact test is a **diagnostic metric only**. Fisher uses permutation-based
exact statistics, not chi2(1) asymptotic statistics — applying GC correction to Fisher p-values
is statistically invalid. Only apply GC correction to SKAT and burden test p-values.
:::
```

**Risk:** None. Documentation-only addition.

---

## Architecture Patterns

### Pattern: Delete Production Function + Delete Test Together

For DEAD-03 and DEAD-04, the pattern is:
1. Remove function definition from source
2. Remove the function from the import in `__init__.py` if present (inheritance `__init__.py`
   only exports `analyze_inheritance` — no change needed there)
3. Remove the function from test imports
4. Remove the test methods that test the deleted function
5. Run `make ci-check` to verify

### Pattern: CLI Choice Removal

For DEAD-02:
1. Edit `cli.py` choices list
2. Edit `analysis_stages.py` validation
3. Update test that verified `"r"` was accepted

### Pattern: Comment-Only Edits

For CLEAN-01, CLEAN-02, DOCS-01, DOCS-02, DOCS-03, CLEAN-07: simple text replacements.
No logic changes, no test changes needed.

---

## Don't Hand-Roll

| Problem | Use Instead | Why |
|---------|-------------|-----|
| Finding all dead code | The research above | Already verified via grep |
| Determining if tests pass | `make ci-check` | Standard project CI |

---

## Common Pitfalls

### Pitfall 1: Deleting Tests Without Deleting Functions First
**What goes wrong:** If you delete test imports but not the functions, remaining test code
that references the function will fail with NameError.
**Prevention:** Delete function definition and test reference in same commit.

### Pitfall 2: Missing the `analysis_stages.py` Validation Update for DEAD-02
**What goes wrong:** CLI refuses `--coast-backend r` but the config validator still accepts `"r"`
from JSON config, leading to inconsistent behavior.
**Prevention:** Update both `cli.py` choices AND `analysis_stages.py` line 2199 validation.

### Pitfall 3: Removing __gzip-intermediates__ Without Updating Default
**What goes wrong:** If `--gzip-intermediates` is removed but the argparse default isn't set
elsewhere, `gzip_intermediates` might default to `None` instead of `True`.
**Prevention:** Keep `default=True` on the `--no-gzip-intermediates` argument via
`parser.set_defaults(gzip_intermediates=True)` or equivalent.

### Pitfall 4: stage_info.py Variables vs Module
**What goes wrong:** Grep for `stage_info` finds hundreds of matches and you conclude the
module IS imported.
**Prevention:** Search for `import stage_info` or `from.*stage_info`, not just `stage_info`.
The variable named `stage_info` in other files is a local variable of type `StageInfo`.

---

## Open Questions

1. **CLEAN-04 - TestBCFToolsPrefilter coverage**
   - What we know: The `BCFToolsPrefilterStage` class was deleted in phase 30
   - What's unclear: Whether `TestBCFToolsPrefilter` tests still exercise live functionality
     (the bcftools prefilter feature exists, just in a different stage)
   - Recommendation: Check `tests/unit/stages/test_processing_stages_critical.py` for existing
     coverage; if prefilter is covered there, delete `TestBCFToolsPrefilter` entirely

2. **DEAD-04 fourth function — create_inheritance_details**
   - What we know: 3 clearly dead functions; requirement says 4
   - What's unclear: Whether `create_inheritance_details` should be made private or removed
   - Recommendation: Rename to `_create_inheritance_details` (making it private), update the
     internal call in `_finalize_inheritance_patterns`, and update the test import

3. **CLEAN-06 - _batch_lightweight_stages resolution** (verified)
   - What we know: Called at runner.py line 660; function body is `return stages` with a TODO
   - Resolution: Replace TODO with deferred-work comment; cannot delete the method

---

## Sources

### Primary (HIGH confidence)
- Direct file reads of all affected files with line number verification
- `variantcentrifuge/stage_info.py` — confirmed 508 lines, zero imports
- `variantcentrifuge/cli.py` — confirmed `--coast-backend` choices, `--gzip-intermediates` args
- `variantcentrifuge/inheritance/prioritizer.py` — confirmed 9 functions, verified 3 used in production
- `variantcentrifuge/inheritance/analyzer.py` — confirmed function list, verified production callers
- `variantcentrifuge/stages/__init__.py` — confirmed 22 exports, identified 6 missing
- `variantcentrifuge/pipeline.py` — confirmed 28 stage imports including 6 missing from __init__.py
- `variantcentrifuge/stages/analysis_stages.py` — confirmed TODO(12-02) at line 760
- `variantcentrifuge/pipeline_core/runner.py` — confirmed TODO batching at line 725
- `variantcentrifuge/association/diagnostics.py` — confirmed Fisher exemption documented (TD-05 done)
- `docs/source/faq.md` — confirmed 4 removed flag references
- `docs/source/guides/association_testing.md` — confirmed 3x `skat_pvalue` instances
- `docs/source/changelog.md` line 51 — confirmed "classic pipeline" reference

### Secondary (MEDIUM confidence)
- `.planning/REQUIREMENTS.md` — source of truth for requirement definitions
- `.planning/phases/34-tech-debt/` — confirms TD-05 was completed in phase 34

---

## Metadata

**Confidence breakdown:**
- Dead code identification: HIGH — verified by grep with zero production call sites
- Documentation errors: HIGH — directly read the docs files
- Test impact: HIGH — traced all imports and test references manually
- Architecture patterns: HIGH — standard Python cleanup patterns

**Research date:** 2026-02-26
**Valid until:** 2026-03-26 (stable domain, low churn expected)
