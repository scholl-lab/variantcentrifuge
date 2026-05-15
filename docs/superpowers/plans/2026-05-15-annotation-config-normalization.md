# Annotation Config Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix custom BED and gene-list annotation key drift so CLI and config-driven runs reliably populate `Custom_Annotation` and checkpoint state changes when annotation inputs change.

**Architecture:** Add one annotation config normalizer that canonicalizes BED inputs to `annotate_bed_files` and gene-list inputs to `annotate_gene_lists`, while preserving legacy aliases. Call it at CLI/config boundaries and use canonical keys inside setup, pipeline construction, output, and checkpoint code.

**Tech Stack:** Python, argparse config dictionaries, pandas DataFrame stages, pytest, existing VariantCentrifuge stage pipeline.

---

## File Structure

- Modify `variantcentrifuge/config.py`
  - Owns annotation key normalization helper.
- Modify `variantcentrifuge/cli.py`
  - Calls the helper after annotation CLI values are copied into `cfg`.
- Modify `variantcentrifuge/pipeline.py`
  - Normalizes config-driven stage creation and uses normalized values for annotation stage decisions.
- Modify `variantcentrifuge/stages/setup_stages.py`
  - Reads canonical normalized annotation keys during normal execution and checkpoint skip restore.
- Modify `variantcentrifuge/stages/output_stages.py`
  - Detects requested custom annotations through canonical normalized keys.
- Modify `variantcentrifuge/checkpoint.py`
  - Includes canonical annotation keys in configuration hashes.
- Create `tests/unit/test_annotation_config_normalization.py`
  - Unit tests for all alias and list-shape rules.
- Modify `tests/test_gene_list_integration.py`
  - Regression test that CLI populates canonical and compatibility gene-list keys.
- Modify `tests/test_cli_argument_parsing.py`
  - Regression test that CLI populates canonical and compatibility BED keys.
- Modify `tests/unit/stages/test_setup_stages.py`
  - Stage-loading tests for canonical BED and gene-list keys.
- Modify `tests/integration/test_create_stages_from_config.py`
  - Stage selection tests for canonical annotation keys.
- Modify `tests/unit/stages/test_output_stages.py`
  - Output-column insertion test for canonical annotation keys.
- Modify `tests/test_checkpoint.py`
  - Hash compatibility tests for canonical annotation keys.
- Modify `tests/unit/stages/test_analysis_stages.py`
  - Small DataFrame plus BED overlap regression proving canonical config reaches `Custom_Annotation`.

---

### Task 1: Add Annotation Config Normalizer

**Files:**
- Modify: `variantcentrifuge/config.py`
- Create: `tests/unit/test_annotation_config_normalization.py`

- [ ] **Step 1: Write failing normalizer tests**

Create `tests/unit/test_annotation_config_normalization.py`:

```python
"""Tests for annotation config key normalization."""

from variantcentrifuge.config import normalize_annotation_config


def test_normalizes_canonical_bed_key_to_all_bed_aliases():
    cfg = {"annotate_bed_files": ["regions.bed"]}

    result = normalize_annotation_config(cfg)

    assert result is cfg
    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_bed"] == ["regions.bed"]


def test_normalizes_legacy_bed_key_to_canonical_key():
    cfg = {"annotate_bed": "regions.bed"}

    result = normalize_annotation_config(cfg)

    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_bed"] == ["regions.bed"]


def test_empty_legacy_bed_key_does_not_mask_canonical_bed_key():
    cfg = {"annotate_bed": [], "annotate_bed_files": ["regions.bed"]}

    result = normalize_annotation_config(cfg)

    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_bed"] == ["regions.bed"]


def test_normalizes_canonical_gene_list_key_to_all_gene_list_aliases():
    cfg = {"annotate_gene_lists": ["genes.txt"]}

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_normalizes_legacy_gene_list_key_to_canonical_key():
    cfg = {"annotate_gene_list": "genes.txt"}

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_normalizes_gene_list_files_key_to_canonical_key():
    cfg = {"annotate_gene_list_files": ["genes.txt"]}

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_empty_legacy_gene_list_key_does_not_mask_canonical_gene_list_key():
    cfg = {
        "annotate_gene_list": [],
        "annotate_gene_list_files": [],
        "annotate_gene_lists": ["genes.txt"],
    }

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_missing_annotation_keys_are_normalized_to_empty_lists():
    result = normalize_annotation_config({})

    assert result["annotate_bed_files"] == []
    assert result["annotate_bed"] == []
    assert result["annotate_gene_lists"] == []
    assert result["annotate_gene_list_files"] == []
    assert result["annotate_gene_list"] == []


def test_filters_empty_items_from_annotation_lists():
    cfg = {
        "annotate_bed_files": ["regions.bed", "", None],
        "annotate_gene_lists": ["genes.txt", "", None],
    }

    result = normalize_annotation_config(cfg)

    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_gene_lists"] == ["genes.txt"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_annotation_config_normalization.py -q
```

Expected: FAIL with `ImportError` or `AttributeError` because `normalize_annotation_config` does not exist yet.

- [ ] **Step 3: Implement the normalizer**

Edit `variantcentrifuge/config.py`:

```python
import json
import os
from collections.abc import Iterable
from typing import Any
```

Add these helpers after `load_config()`:

```python
def _annotation_value_as_list(value: Any) -> list[str]:
    """Normalize annotation config values to a list of non-empty strings."""
    if value is None:
        return []
    if isinstance(value, str):
        stripped = value.strip()
        return [stripped] if stripped else []
    if isinstance(value, Iterable):
        normalized = []
        for item in value:
            if item is None:
                continue
            item_str = str(item).strip()
            if item_str:
                normalized.append(item_str)
        return normalized
    item_str = str(value).strip()
    return [item_str] if item_str else []


def _first_non_empty_annotation_value(config: dict[str, Any], keys: list[str]) -> list[str]:
    """Return the first non-empty normalized annotation value for key aliases."""
    for key in keys:
        value = _annotation_value_as_list(config.get(key))
        if value:
            return value
    return []


def normalize_annotation_config(config: dict[str, Any]) -> dict[str, Any]:
    """Normalize custom annotation config aliases in-place.

    Canonical internal keys:
    - annotate_bed_files
    - annotate_gene_lists

    Compatibility aliases:
    - annotate_bed
    - annotate_gene_list
    - annotate_gene_list_files
    """
    bed_files = _first_non_empty_annotation_value(
        config,
        ["annotate_bed_files", "annotate_bed"],
    )
    gene_lists = _first_non_empty_annotation_value(
        config,
        ["annotate_gene_lists", "annotate_gene_list_files", "annotate_gene_list"],
    )

    config["annotate_bed_files"] = bed_files
    config["annotate_bed"] = bed_files
    config["annotate_gene_lists"] = gene_lists
    config["annotate_gene_list_files"] = gene_lists
    config["annotate_gene_list"] = gene_lists
    return config
```

- [ ] **Step 4: Run normalizer tests to verify they pass**

Run:

```bash
pytest tests/unit/test_annotation_config_normalization.py -q
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add variantcentrifuge/config.py tests/unit/test_annotation_config_normalization.py
git commit -m "fix: normalize annotation config aliases"
```

---

### Task 2: Normalize CLI Annotation Config

**Files:**
- Modify: `variantcentrifuge/cli.py`
- Modify: `tests/test_gene_list_integration.py`
- Modify: `tests/test_cli_argument_parsing.py`

- [ ] **Step 1: Add failing CLI BED config regression test**

Append to `tests/test_cli_argument_parsing.py`:

```python
def test_cli_annotate_bed_populates_canonical_and_legacy_config_keys(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    captured_configs = []

    def fake_run_pipeline(args):
        captured_configs.append(args.config.copy())

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t99\t101\tregion\n")

    monkeypatch.setattr(cli_module, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--annotate-bed",
            str(bed),
            "--output-dir",
            str(tmp_path),
        ],
    )

    exit_code = cli_module.main()

    assert exit_code == 0
    assert len(captured_configs) == 1
    cfg = captured_configs[0]
    assert cfg["annotate_bed_files"] == [str(bed)]
    assert cfg["annotate_bed"] == [str(bed)]
```

- [ ] **Step 2: Extend existing gene-list CLI test**

In `tests/test_gene_list_integration.py`, extend the assertions around the captured `cfg`:

```python
        assert "annotate_gene_list_files" in cfg
        assert "annotate_gene_lists" in cfg
        assert "annotate_gene_list" in cfg
        assert len(cfg["annotate_gene_list_files"]) == 2
        assert cfg["annotate_gene_lists"] == cfg["annotate_gene_list_files"]
        assert cfg["annotate_gene_list"] == cfg["annotate_gene_list_files"]
        assert gene_lists["cancer_genes"] in cfg["annotate_gene_list_files"]
        assert gene_lists["apc_genes"] in cfg["annotate_gene_list_files"]
```

- [ ] **Step 3: Run CLI tests to verify they fail**

Run:

```bash
pytest tests/test_cli_argument_parsing.py::test_cli_annotate_bed_populates_canonical_and_legacy_config_keys tests/test_gene_list_integration.py::test_cli_gene_list_integration -q
```

Expected: BED test FAILS because `annotate_bed` is absent or empty in `args.config`. Gene-list test may fail because `annotate_gene_list` is absent or not normalized.

- [ ] **Step 4: Call normalizer from CLI**

Edit `variantcentrifuge/cli.py`.

Add import near existing local imports from `variantcentrifuge.config` if present, or near other top-level imports:

```python
from .config import normalize_annotation_config
```

After these existing config assignments:

```python
    cfg["annotate_bed_files"] = args.annotate_bed
    cfg["annotate_json_genes"] = args.annotate_json_genes
    cfg["json_gene_mapping"] = args.json_gene_mapping
    cfg["json_genes_as_columns"] = args.json_genes_as_columns
    # Note: args.annotate_gene_list is already handled above as "annotate_gene_list_files"
    # We'll map it to the new unified system
    cfg["annotate_gene_lists"] = args.annotate_gene_list
```

Add:

```python
    normalize_annotation_config(cfg)
```

- [ ] **Step 5: Run CLI tests to verify they pass**

Run:

```bash
pytest tests/test_cli_argument_parsing.py::test_cli_annotate_bed_populates_canonical_and_legacy_config_keys tests/test_gene_list_integration.py::test_cli_gene_list_integration -q
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add variantcentrifuge/cli.py tests/test_cli_argument_parsing.py tests/test_gene_list_integration.py
git commit -m "fix: normalize CLI annotation config"
```

---

### Task 3: Normalize Stage Selection And Setup Loading

**Files:**
- Modify: `variantcentrifuge/pipeline.py`
- Modify: `variantcentrifuge/stages/setup_stages.py`
- Modify: `tests/unit/stages/test_setup_stages.py`
- Modify: `tests/integration/test_create_stages_from_config.py`

- [ ] **Step 1: Add failing setup-stage tests for canonical keys**

Append to `TestAnnotationConfigLoadingStage` in `tests/unit/stages/test_setup_stages.py`:

```python
    def test_annotation_loading_from_canonical_keys(self, context):
        """Canonical annotation keys must load into annotation_configs."""
        bed_file = context.config["annotate_bed"][0]
        gene_list = context.config["annotate_gene_list"][0]
        context.config["annotate_bed"] = []
        context.config["annotate_gene_list"] = []
        context.config["annotate_bed_files"] = [bed_file]
        context.config["annotate_gene_lists"] = [gene_list]

        stage = AnnotationConfigLoadingStage()
        result = stage(context)

        assert result.annotation_configs["bed_files"] == [bed_file]
        assert result.annotation_configs["gene_lists"] == [gene_list]

    def test_checkpoint_skip_restores_canonical_annotation_keys(self, context):
        """Checkpoint skip restore must use canonical annotation keys."""
        bed_file = context.config["annotate_bed"][0]
        gene_list = context.config["annotate_gene_list"][0]
        context.config["annotate_bed"] = []
        context.config["annotate_gene_list"] = []
        context.config["annotate_bed_files"] = [bed_file]
        context.config["annotate_gene_lists"] = [gene_list]

        stage = AnnotationConfigLoadingStage()
        result = stage._handle_checkpoint_skip(context)

        assert result.annotation_configs["bed_files"] == [bed_file]
        assert result.annotation_configs["gene_lists"] == [gene_list]
```

- [ ] **Step 2: Add failing config-driven stage selection tests**

Append to `tests/integration/test_create_stages_from_config.py`:

```python
from variantcentrifuge.stages.analysis_stages import CustomAnnotationStage
from variantcentrifuge.stages.setup_stages import AnnotationConfigLoadingStage


@pytest.mark.integration
def test_canonical_annotate_bed_files_activates_custom_annotation_stages():
    stage_types = _stage_types({"annotate_bed_files": ["regions.bed"]})

    assert AnnotationConfigLoadingStage in stage_types
    assert CustomAnnotationStage in stage_types


@pytest.mark.integration
def test_canonical_annotate_gene_lists_activates_custom_annotation_stages():
    stage_types = _stage_types({"annotate_gene_lists": ["genes.txt"]})

    assert AnnotationConfigLoadingStage in stage_types
    assert CustomAnnotationStage in stage_types


@pytest.mark.integration
def test_gene_list_files_alias_activates_custom_annotation_stages():
    stage_types = _stage_types({"annotate_gene_list_files": ["genes.txt"]})

    assert AnnotationConfigLoadingStage in stage_types
    assert CustomAnnotationStage in stage_types
```

- [ ] **Step 3: Run stage tests to verify they fail**

Run:

```bash
pytest tests/unit/stages/test_setup_stages.py::TestAnnotationConfigLoadingStage tests/integration/test_create_stages_from_config.py::test_canonical_annotate_bed_files_activates_custom_annotation_stages tests/integration/test_create_stages_from_config.py::test_canonical_annotate_gene_lists_activates_custom_annotation_stages tests/integration/test_create_stages_from_config.py::test_gene_list_files_alias_activates_custom_annotation_stages -q
```

Expected: canonical-key setup tests and config-driven stage selection tests FAIL.

- [ ] **Step 4: Normalize config in pipeline stage creation**

Edit `variantcentrifuge/pipeline.py`.

Add import:

```python
from .config import load_config, normalize_annotation_config
```

Replace the existing `from .config import load_config` import with the combined import above.

In `build_pipeline_stages(args)`, after:

```python
    config = args.config if hasattr(args, "config") and isinstance(args.config, dict) else {}
```

Add:

```python
    if config:
        normalize_annotation_config(config)
```

Change both annotation-stage condition blocks from argparse-only checks to config-aware checks:

```python
    annotations_requested = any(
        [
            config.get("annotate_bed_files"),
            config.get("annotate_gene_lists"),
            config.get("annotate_json_genes"),
            hasattr(args, "annotate_bed") and args.annotate_bed,
            hasattr(args, "annotate_gene_list") and args.annotate_gene_list,
            hasattr(args, "annotate_json_genes") and args.annotate_json_genes,
        ]
    )
```

Use `annotations_requested` for both the setup-stage insertion and the `CustomAnnotationStage` insertion.

In `create_stages_from_config(config)`, add at the start:

```python
    config = normalize_annotation_config(config.copy())
```

Then change:

```python
    args.annotate_bed = config.get("annotate_bed", [])
    args.annotate_gene_list = config.get("annotate_gene_list", [])
```

to:

```python
    args.annotate_bed = config.get("annotate_bed_files", [])
    args.annotate_gene_list = config.get("annotate_gene_lists", [])
```

- [ ] **Step 5: Read canonical keys in setup stage**

Edit `variantcentrifuge/stages/setup_stages.py`.

Add import:

```python
from ..config import load_config, normalize_annotation_config
```

Replace the existing `from ..config import load_config` import with the combined import above.

At the start of `AnnotationConfigLoadingStage._process()` after `annotations = {}`:

```python
        normalize_annotation_config(context.config)
```

Change:

```python
        bed_files = context.config.get("annotate_bed", [])
```

to:

```python
        bed_files = context.config.get("annotate_bed_files", [])
```

Change:

```python
        gene_lists = context.config.get("annotate_gene_list", [])
```

to:

```python
        gene_lists = context.config.get("annotate_gene_lists", [])
```

Repeat the same normalization and key changes inside `_handle_checkpoint_skip()`.

- [ ] **Step 6: Run stage tests to verify they pass**

Run:

```bash
pytest tests/unit/stages/test_setup_stages.py::TestAnnotationConfigLoadingStage tests/integration/test_create_stages_from_config.py -q
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add variantcentrifuge/pipeline.py variantcentrifuge/stages/setup_stages.py tests/unit/stages/test_setup_stages.py tests/integration/test_create_stages_from_config.py
git commit -m "fix: route annotation aliases through stages"
```

---

### Task 4: Fix Output Column Detection And Checkpoint Hashing

**Files:**
- Modify: `variantcentrifuge/stages/output_stages.py`
- Modify: `variantcentrifuge/checkpoint.py`
- Modify: `tests/unit/stages/test_output_stages.py`
- Modify: `tests/test_checkpoint.py`

- [ ] **Step 1: Add failing output-stage test for canonical annotation keys**

Append to `TestVariantIdentifierStage` in `tests/unit/stages/test_output_stages.py`:

```python
    def test_adds_custom_annotation_column_for_canonical_annotation_keys(self, context):
        """Canonical annotation config keys should request Custom_Annotation output."""
        context.config["annotate_bed_files"] = ["regions.bed"]
        context.config["annotate_gene_lists"] = ["genes.txt"]
        context.config["annotate_bed"] = []
        context.config["annotate_gene_list"] = []

        stage = VariantIdentifierStage()
        result = stage(context)

        assert "Custom_Annotation" in result.current_dataframe.columns
```

- [ ] **Step 2: Add failing checkpoint hash tests**

Append to `TestPipelineState` in `tests/test_checkpoint.py`:

```python
    def test_pipeline_state_resume_detects_canonical_bed_annotation_change(self, tmp_path):
        """Changing annotate_bed_files must invalidate resume compatibility."""
        state = PipelineState(str(tmp_path))
        config = {
            "gene_name": "BRCA1",
            "vcf_file": "input.vcf",
            "annotate_bed_files": ["regions_a.bed"],
        }
        state.initialize(config, "1.0.0")
        state.save()

        new_state = PipelineState(str(tmp_path))
        new_state.load()

        changed_config = {
            "gene_name": "BRCA1",
            "vcf_file": "input.vcf",
            "annotate_bed_files": ["regions_b.bed"],
        }
        assert new_state.can_resume(changed_config, "1.0.0") is False

    def test_pipeline_state_resume_detects_canonical_gene_list_annotation_change(self, tmp_path):
        """Changing annotate_gene_lists must invalidate resume compatibility."""
        state = PipelineState(str(tmp_path))
        config = {
            "gene_name": "BRCA1",
            "vcf_file": "input.vcf",
            "annotate_gene_lists": ["genes_a.txt"],
        }
        state.initialize(config, "1.0.0")
        state.save()

        new_state = PipelineState(str(tmp_path))
        new_state.load()

        changed_config = {
            "gene_name": "BRCA1",
            "vcf_file": "input.vcf",
            "annotate_gene_lists": ["genes_b.txt"],
        }
        assert new_state.can_resume(changed_config, "1.0.0") is False
```

- [ ] **Step 3: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/stages/test_output_stages.py::TestVariantIdentifierStage::test_adds_custom_annotation_column_for_canonical_annotation_keys tests/test_checkpoint.py::TestPipelineState::test_pipeline_state_resume_detects_canonical_bed_annotation_change tests/test_checkpoint.py::TestPipelineState::test_pipeline_state_resume_detects_canonical_gene_list_annotation_change -q
```

Expected: output-stage test FAILS because it checks only legacy keys. Checkpoint tests FAIL because the hash ignores canonical keys.

- [ ] **Step 4: Normalize output-stage annotation detection**

Edit `variantcentrifuge/stages/output_stages.py`.

Add import:

```python
from ..config import normalize_annotation_config
```

In `VariantIdentifierStage._process()`, before `custom_annotations_requested = any(...)`, add:

```python
        normalize_annotation_config(context.config)
```

Change the request check to canonical keys:

```python
        custom_annotations_requested = any(
            [
                context.config.get("annotate_bed_files", []),
                context.config.get("annotate_gene_lists", []),
                context.config.get("annotate_json_genes", []),
            ]
        )
```

- [ ] **Step 5: Normalize checkpoint hash input**

Edit `variantcentrifuge/checkpoint.py`.

Add import near other local imports:

```python
from .config import normalize_annotation_config
```

At the start of `_hash_configuration()`:

```python
        config = normalize_annotation_config(config.copy())
```

Add canonical and compatibility annotation keys to `relevant_keys`:

```python
            "annotate_bed",
            "annotate_bed_files",
            "annotate_gene_list",
            "annotate_gene_list_files",
            "annotate_gene_lists",
            "annotate_json_genes",
```

Ensure the final list contains these keys once each. Keeping both canonical and compatibility keys is acceptable after normalization because they carry equal values and preserve compatibility with existing stored state behavior.

- [ ] **Step 6: Run output and checkpoint tests to verify they pass**

Run:

```bash
pytest tests/unit/stages/test_output_stages.py::TestVariantIdentifierStage::test_adds_custom_annotation_column_for_canonical_annotation_keys tests/test_checkpoint.py::TestPipelineState -q
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add variantcentrifuge/stages/output_stages.py variantcentrifuge/checkpoint.py tests/unit/stages/test_output_stages.py tests/test_checkpoint.py
git commit -m "fix: include annotations in output and checkpoint config"
```

---

### Task 5: Add End-To-End Stage Regression For BED And Gene Lists

**Files:**
- Modify: `tests/unit/stages/test_analysis_stages.py`

- [ ] **Step 1: Add failing DataFrame BED overlap regression**

Append to `TestCustomAnnotationStage` in `tests/unit/stages/test_analysis_stages.py`:

```python
    def test_canonical_bed_config_reaches_custom_annotation(self, base_context, tmp_path):
        """A canonical BED config should annotate overlapping variants."""
        from variantcentrifuge.stages.setup_stages import AnnotationConfigLoadingStage

        bed_file = tmp_path / "GIAB_TANDEM_REPEAT.bed"
        bed_file.write_text("chr1\t99\t101\tTR\n")
        base_context.config["annotate_bed_files"] = [str(bed_file)]
        base_context.config["annotate_bed"] = []
        base_context.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 500],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["HTT", "PKD1"],
            }
        )

        AnnotationConfigLoadingStage()._process(base_context)
        result = CustomAnnotationStage()._process(base_context)

        assert "Custom_Annotation" in result.current_dataframe.columns
        assert "Region=TR_GIAB_TANDEM_REPEAT" in result.current_dataframe.iloc[0][
            "Custom_Annotation"
        ]
        assert result.current_dataframe.iloc[1]["Custom_Annotation"] == ""
```

- [ ] **Step 2: Add failing canonical gene-list regression**

Append to `TestCustomAnnotationStage` in `tests/unit/stages/test_analysis_stages.py`:

```python
    def test_canonical_gene_list_config_reaches_custom_annotation(self, base_context, tmp_path):
        """A canonical gene-list config should annotate matching genes."""
        from variantcentrifuge.stages.setup_stages import AnnotationConfigLoadingStage

        gene_file = tmp_path / "repeat_review_genes.txt"
        gene_file.write_text("HTT\n")
        base_context.config["annotate_gene_lists"] = [str(gene_file)]
        base_context.config["annotate_gene_list"] = []
        base_context.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr16"],
                "POS": [100, 200],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["HTT", "PKD1"],
            }
        )

        AnnotationConfigLoadingStage()._process(base_context)
        result = CustomAnnotationStage()._process(base_context)

        assert "Custom_Annotation" in result.current_dataframe.columns
        assert "InGeneList=repeat_review_genes" in result.current_dataframe.iloc[0][
            "Custom_Annotation"
        ]
        assert result.current_dataframe.iloc[1]["Custom_Annotation"] == ""
```

- [ ] **Step 3: Run analysis-stage regressions**

Run:

```bash
pytest tests/unit/stages/test_analysis_stages.py::TestCustomAnnotationStage::test_canonical_bed_config_reaches_custom_annotation tests/unit/stages/test_analysis_stages.py::TestCustomAnnotationStage::test_canonical_gene_list_config_reaches_custom_annotation -q
```

Expected before Tasks 1-4 are complete: FAIL. Expected after Tasks 1-4: PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/unit/stages/test_analysis_stages.py
git commit -m "test: cover canonical custom annotation flow"
```

---

### Task 6: Full Verification

**Files:**
- No code files changed in this task.

- [ ] **Step 1: Run focused annotation regression suite**

Run:

```bash
pytest tests/unit/test_annotation_config_normalization.py tests/test_cli_argument_parsing.py::test_cli_annotate_bed_populates_canonical_and_legacy_config_keys tests/test_gene_list_integration.py::test_cli_gene_list_integration tests/unit/stages/test_setup_stages.py::TestAnnotationConfigLoadingStage tests/integration/test_create_stages_from_config.py tests/unit/stages/test_output_stages.py::TestVariantIdentifierStage tests/test_checkpoint.py::TestPipelineState tests/unit/stages/test_analysis_stages.py::TestCustomAnnotationStage -q
```

Expected: PASS.

- [ ] **Step 2: Run broader existing annotation tests**

Run:

```bash
pytest tests/test_annotator.py tests/test_gene_list_annotation.py tests/test_gene_list_integration.py -q
```

Expected: PASS.

- [ ] **Step 3: Run pipeline-stage smoke tests**

Run:

```bash
pytest tests/unit/stages/test_setup_stages.py tests/unit/stages/test_analysis_stages.py tests/unit/stages/test_output_stages.py tests/integration/test_create_stages_from_config.py -q
```

Expected: PASS.

- [ ] **Step 4: Inspect changed files**

Run:

```bash
git diff --stat HEAD
git diff --check
```

Expected: `git diff --check` exits 0 with no whitespace errors.

- [ ] **Step 5: Final commit if previous task commits were skipped**

If implementation was done without per-task commits, create one scoped commit:

```bash
git add variantcentrifuge/config.py variantcentrifuge/cli.py variantcentrifuge/pipeline.py variantcentrifuge/stages/setup_stages.py variantcentrifuge/stages/output_stages.py variantcentrifuge/checkpoint.py tests/unit/test_annotation_config_normalization.py tests/test_cli_argument_parsing.py tests/test_gene_list_integration.py tests/unit/stages/test_setup_stages.py tests/integration/test_create_stages_from_config.py tests/unit/stages/test_output_stages.py tests/test_checkpoint.py tests/unit/stages/test_analysis_stages.py
git commit -m "fix: normalize custom annotation config keys"
```

---

## Plan Self-Review

Spec coverage:

- BED alias drift is covered in Tasks 1, 2, 3, 4, and 5.
- Gene-list alias drift is covered in Tasks 1, 2, 3, 4, and 5.
- CLI config is covered in Task 2.
- Config-driven stage creation is covered in Task 3.
- Setup-stage normal and checkpoint-skip paths are covered in Task 3.
- Output-column insertion is covered in Task 4.
- Checkpoint invalidation is covered in Task 4.
- DataFrame annotation behavior is covered in Task 5.

Placeholder scan:

- No placeholder steps remain.
- Every code-changing step includes concrete code.
- Every test step has a command and expected result.

Type consistency:

- Canonical BED key is `annotate_bed_files` throughout.
- Canonical gene-list key is `annotate_gene_lists` throughout.
- Compatibility gene-list aliases are `annotate_gene_list` and `annotate_gene_list_files` throughout.
