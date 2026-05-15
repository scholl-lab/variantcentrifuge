# Annotation Config Normalization Design

## Goal

Fix the silent custom annotation failure reported in GitHub issue #95 by making BED and gene-list annotation configuration use one normalized internal contract across CLI, config-driven pipeline construction, setup stages, output stages, checkpoint hashing, and the low-level annotator.

This spec deliberately does not implement the larger issue #94 feature request for first-class technical-region boolean columns. It restores and hardens the existing `Custom_Annotation` behavior so issue #94 can be built on a reliable base.

## Background

The current codebase has two incompatible key families for the same annotation inputs:

- BED annotations:
  - Legacy/public key: `annotate_bed`
  - Current low-level annotator key: `annotate_bed_files`
- Gene-list annotations:
  - Legacy/public key: `annotate_gene_list`
  - Older helper key: `annotate_gene_list_files`
  - Current low-level annotator key: `annotate_gene_lists`

The CLI parses `--annotate-bed` and `--annotate-gene-list` correctly. During config assembly, it writes BED inputs to `annotate_bed_files`, gene-list inputs to `annotate_gene_list_files`, and gene-list inputs again to `annotate_gene_lists`. Later pipeline code still reads only `annotate_bed` and `annotate_gene_list` in important places.

The result is a silent analysis bug: annotation stages may be present because they were selected from argparse attributes, but `AnnotationConfigLoadingStage` sees no BED files in `context.config`, so `context.annotation_configs` remains empty and `CustomAnnotationStage` does nothing.

## Evidence

Local reproduction:

- `AnnotationConfigLoadingStage` with `{"annotate_bed": [bed]}` loads `annotation_configs["bed_files"]`.
- `AnnotationConfigLoadingStage` with `{"annotate_bed_files": [bed]}` leaves `annotation_configs` empty.
- `load_custom_features({"annotate_bed_files": [bed]})` loads BED intervals successfully.
- `create_stages_from_config({"annotate_bed_files": ["x.bed"]})` omits annotation stages entirely.

This means the BED file contents and `intervaltree` dependency are not the root cause. The root cause is config key drift.

## Scope

In scope:

- Normalize BED and gene-list annotation config keys.
- Preserve backward compatibility with existing CLI flags and config files.
- Make stage selection, setup-stage loading, checkpoint skip restore, output-column insertion, validation, and checkpoint hashing use normalized annotation config.
- Add focused regression tests covering BED and gene-list aliases.
- Add a small DataFrame-level BED overlap regression test proving a CLI-style config reaches `Custom_Annotation`.

Out of scope:

- New named BED syntax such as `name=path.bed`.
- Dedicated technical-region yes/no columns from issue #94.
- Changes to BED interval semantics.
- Changes to existing JSON gene annotation behavior, except preserving it while normalizing BED and gene-list keys.
- Snakemake workflow feature expansion for BED lists. The workflow currently only exposes `gene_list_files`; this bug fix should not broaden workflow schema unless a later issue asks for it.

## Design Options

### Option A: Central Normalizer, Keep Compatibility Aliases

Create a small helper in `variantcentrifuge/config.py`:

- `normalize_annotation_config(config: dict[str, Any]) -> dict[str, Any]`
- It mutates and returns `config`.
- It normalizes scalar strings and list-like inputs into lists.
- It chooses canonical internal keys:
  - `annotate_bed_files`
  - `annotate_gene_lists`
- It also populates compatibility aliases:
  - `annotate_bed`
  - `annotate_gene_list`
  - `annotate_gene_list_files`

All pipeline code can then consume either the canonical helper result or the canonical keys. This is the recommended approach because it removes the root cause at the boundary while keeping old config files and tests valid.

### Option B: Patch Every Read Site With OR Logic

Change every read from `config.get("annotate_bed")` to `config.get("annotate_bed") or config.get("annotate_bed_files")`, and similar for gene lists.

This is smaller initially but spreads alias logic across the codebase. It makes the next annotation option more likely to repeat the same drift.

### Option C: Hard Rename Everything To Canonical Keys

Remove legacy keys and update all code to use only `annotate_bed_files` and `annotate_gene_lists`.

This is the cleanest internal model but breaks existing config files and tests. It is not appropriate for a bug fix on a user-facing CLI option.

## Selected Design

Use Option A.

The normalization helper is the single source of truth for alias handling. Public CLI flags remain unchanged. Existing JSON config files using old keys continue to work. Internal annotator code continues to receive the canonical keys it already expects.

## Normalization Rules

BED aliases are merged in this precedence order:

1. `annotate_bed_files`
2. `annotate_bed`

Gene-list aliases are merged in this precedence order:

1. `annotate_gene_lists`
2. `annotate_gene_list_files`
3. `annotate_gene_list`

The first non-empty value wins. This prevents duplicate annotation if CLI config assembly has already written the same inputs under several aliases. Empty values do not mask later non-empty aliases.

Accepted input shapes:

- `None` -> `[]`
- `""` -> `[]`
- `"path.bed"` -> `["path.bed"]`
- `["a.bed", "b.bed"]` -> unchanged list without empty entries
- Tuple/set inputs may be converted to lists defensively.

Output keys after normalization:

```python
cfg["annotate_bed_files"] == bed_files
cfg["annotate_bed"] == bed_files
cfg["annotate_gene_lists"] == gene_lists
cfg["annotate_gene_list_files"] == gene_lists
cfg["annotate_gene_list"] == gene_lists
```

## Data Flow

CLI run:

1. `create_parser()` parses `--annotate-bed` into `args.annotate_bed` and `--annotate-gene-list` into `args.annotate_gene_list`.
2. `cli.main()` writes parsed values into the config dict.
3. `normalize_annotation_config(cfg)` fills canonical and compatibility keys.
4. `run_pipeline()` receives `args.config = cfg`.
5. `build_pipeline_stages()` sees annotation values and adds `AnnotationConfigLoadingStage` and `CustomAnnotationStage`.
6. `AnnotationConfigLoadingStage` reads normalized canonical keys and writes `context.annotation_configs`.
7. `CustomAnnotationStage` passes canonical keys to `load_custom_features()`.
8. `annotate_dataframe_with_features()` writes `Custom_Annotation`.

Config-driven run:

1. `create_stages_from_config(config)` normalizes the config before creating its compatibility namespace.
2. Annotation stages are added when any BED, gene-list, or JSON gene annotation is configured.

Resume/checkpoint run:

1. Checkpoint config hashing includes canonical annotation keys.
2. Changing BED or gene-list inputs changes the hash.
3. `_handle_checkpoint_skip()` restores `context.annotation_configs` from normalized canonical keys.

## Files To Change

- `variantcentrifuge/config.py`
  - Add `_as_list()` or similar private list coercion helper.
  - Add `normalize_annotation_config()`.

- `variantcentrifuge/cli.py`
  - Import and call `normalize_annotation_config()` after annotation args are copied into `cfg`.

- `variantcentrifuge/pipeline.py`
  - Normalize config in `create_stages_from_config()`.
  - Use normalized config values when building stage lists.

- `variantcentrifuge/stages/setup_stages.py`
  - Read canonical normalized keys in `_process()`.
  - Read canonical normalized keys in `_handle_checkpoint_skip()`.

- `variantcentrifuge/stages/output_stages.py`
  - Treat either normalized BED or gene-list inputs as a custom annotation request, preferably via canonical keys after normalization.

- `variantcentrifuge/checkpoint.py`
  - Normalize config before hashing, or add canonical annotation keys to the hash after normalizing.

- Tests:
  - `tests/unit/test_annotation_config_normalization.py`
  - `tests/unit/stages/test_setup_stages.py`
  - `tests/integration/test_create_stages_from_config.py`
  - `tests/unit/stages/test_output_stages.py` or `tests/unit/stages/test_output_stages_simple.py`
  - `tests/test_checkpoint.py`
  - `tests/test_gene_list_integration.py`
  - Optional: `tests/unit/stages/test_analysis_stages.py` for a DataFrame-level BED overlap regression.

## Error Handling

This fix should not add new fatal validation. Existing behavior logs warnings for missing BED or gene-list files in the low-level loaders and validation helpers. The normalizer should not check file existence. It should only normalize shape and key aliases.

If conflicting non-empty aliases are supplied with different values, the first non-empty alias in the precedence order wins. This is deterministic and avoids accidental duplication. A debug log can be added later, but it is not necessary for this bug fix.

## Testing Requirements

The tests must prove the bug before fixing it and protect the fixed config contract afterward:

- Canonical BED config reaches `AnnotationConfigLoadingStage`.
- Legacy BED config still works.
- Canonical gene-list config reaches `AnnotationConfigLoadingStage`.
- Legacy gene-list config still works.
- `annotate_gene_list_files` still works for older integration code.
- Empty legacy aliases do not mask non-empty canonical aliases.
- `create_stages_from_config()` adds annotation stages for canonical BED and gene-list keys.
- CLI parsing normalizes both BED and gene-list aliases into `args.config`.
- Checkpoint hash changes when canonical BED or gene-list lists change.
- A tiny BED interval overlapping a tiny DataFrame variant produces a non-empty `Custom_Annotation`.

## Acceptance Criteria

- A command with one or more `--annotate-bed` inputs results in non-empty BED features being loaded when files are non-empty.
- A variant overlapping an input BED interval receives a `Region=...` entry in `Custom_Annotation`.
- A command with `--annotate-gene-list` still supports existing behavior and reaches both the unified annotator and legacy-compatible config keys.
- Config-only runs using `annotate_bed_files` or `annotate_gene_lists` activate custom annotation stages.
- Config-only runs using `annotate_bed`, `annotate_gene_list`, or `annotate_gene_list_files` remain supported.
- Changing BED or gene-list annotation inputs invalidates checkpoint resume compatibility.
- Existing JSON gene annotation behavior remains unchanged.

## Self-Review

- Placeholder scan: no placeholders remain.
- Scope check: the design is focused on issue #95 and gene-list alias drift. It intentionally excludes issue #94 boolean technical-region columns.
- Ambiguity check: canonical keys and alias precedence are explicit.
- Consistency check: all consuming code is directed to canonical keys after normalization.
