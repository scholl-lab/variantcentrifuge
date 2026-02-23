# Phase 29 Research: Classic Pipeline Deprecation and Removal

## Key Finding

**The classic pipeline code path no longer exists.** The current `pipeline.py` IS the stage-based pipeline — it imports from `pipeline_core/` and uses `PipelineRunner`, `Stage`, `PipelineContext`, etc. The `--use-new-pipeline` flag is **not defined in argparse** and has no effect. Phase 29 is primarily a **cleanup and documentation task**, not a code removal task.

## Current State

### cli.py Routing (No Branching)
- Line 12: `from .pipeline import run_refactored_pipeline`
- Line 1623: Calls `run_refactored_pipeline(refactored_args)` unconditionally
- **No `if use_new_pipeline` check exists**
- Two additional imports of `create_stages_from_config` (lines 1462, 1501)

### pipeline.py IS Stage-Based (590 lines)
- Docstring: "Refactored pipeline using the new Stage architecture"
- Imports: PipelineContext, PipelineRunner, Stage, Workspace from pipeline_core/
- Imports all stages from stages/ directory
- Functions: `run_refactored_pipeline()`, `build_pipeline_stages()`, `create_stages_from_config()`, `check_scoring_requires_inheritance()`
- **No monolithic/classic code path remains**

### --use-new-pipeline Flag
- **NOT defined in argparse** (searched entire cli.py, 1634 lines)
- Never parsed, never checked — it's a documentation phantom
- Still referenced in tests and docs as if it exists

## Files Requiring Changes

### Source Code
| File | Change | Lines |
|------|--------|-------|
| `variantcentrifuge/pipeline.py` | Rename docstring; remove "refactored" naming; clean up `main()` demo function | Docstring, function names |
| `variantcentrifuge/cli.py` | Rename `run_refactored_pipeline` → `run_pipeline` (optional clarity) | L12, L1623 |

### Tests
| File | Change | Lines |
|------|--------|-------|
| `tests/conftest.py` | Remove `--use-new-pipeline` from gene_burden_cmd_template | L176 |
| `tests/test_gene_burden_comprehensive.py` | Remove `--use-new-pipeline` from 8 test commands | L82,120,164,202,242,277,325,360 |
| `tests/test_cli.py` | Update docstrings mentioning `--use-new-pipeline` | L60, L176 |
| `tests/integration/test_pipeline_with_mocked_tools.py` | Update import if renamed | L15 |
| `tests/integration/test_inheritance_analysis_integration.py` | Update import if renamed | L15 |

### Documentation
| File | Change | Lines |
|------|--------|-------|
| `CLAUDE.md` | Remove "Two Pipeline Modes" section; describe single stage-based architecture | L34-39 |
| `README.md` | Remove `(--use-new-pipeline)` reference | L26 |
| `tests/fixtures/geneburden/README.md` | Remove `--use-new-pipeline` from 5 examples | L73,87,99,114,129 |
| `docs/source/api/scoring.md` | Update `run_pipeline` import reference | L153,168 |
| `docs/source/development.md` | Update `run_pipeline` reference | L351 |

## Risk Assessment

**Risk: LOW** — This is cleanup, not architectural change.

1. **No functional code paths change** — pipeline.py already IS stage-based
2. **Flag removal is safe** — flag is not parsed; removing references won't break behavior
3. **Rename is optional** — `run_refactored_pipeline` → `run_pipeline` improves clarity but isn't required
4. **Test changes are string removal** — removing vestigial flag from command strings
5. **Snakemake unaffected** — workflow uses CLI entry point, no pipeline.py imports, no flag

## Recommended Plan Structure

Given the scope is smaller than originally anticipated (cleanup vs removal):

- **Plan 29-01**: Clean up pipeline.py naming + cli.py imports (rename refactored→standard)
- **Plan 29-02**: Remove all --use-new-pipeline references from tests and fixtures
- **Plan 29-03**: Update CLAUDE.md, README.md, and docs to reflect single pipeline architecture

All three plans are Wave 1 (independent, can execute in parallel).
