# Transcript Filtering and MANE Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `--transcript-list` and `--transcript-file` actually restrict extracted SnpEff annotations to requested transcripts, and document the supported MANE workflow.

**Architecture:** Implement transcript filtering as a VCF INFO/ANN transform before field extraction, so the selected transcript becomes `ANN[0]` and the existing extraction pipeline continues to work. Add the stage to both single-thread and parallel chunk processing, fix the `--split-snpeff-lines` config drift, and document MANE as a transcript-list driven workflow unless the input VCF carries explicit MANE metadata.

**Tech Stack:** Python 3.12, pandas, pytest, bgzip/tabix/bcftools via Conda, existing stage pipeline in `variantcentrifuge/stages/processing_stages.py`.

---

## Environment Check

Local Conda environments currently installed:

```text
base                     C:\Users\bernt\anaconda3
muconeup-manuscript      C:\Users\bernt\anaconda3\envs\muconeup-manuscript
snakemake                C:\Users\bernt\anaconda3\envs\snakemake
vntyper-manuscript       C:\Users\bernt\anaconda3\envs\vntyper-manuscript
```

Project-declared environments:

- `environment_local.yml` declares `name: annotation` and includes `bcftools=1.21`, `bedtools=2.31.1`, `snpeff`, `snpsift`, pytest, and runtime dependencies.
- `conda/environment.yml` declares `name: variantcentrifuge` and includes `bcftools=1.21`, `bedtools=2.31.1`, `snpeff=5.2`, `snpsift=5.2`, Python 3.12, pytest, ruff, mypy, and editable install.

Observed problem:

- Neither `annotation` nor `variantcentrifuge` exists locally.
- Existing `snakemake` env has `Python 3.13.12` and `pytest 9.0.2`, but does not have `bcftools`, `SnpSift`, or `snpEff`.

Use `conda/environment.yml` as the proper development env because it is concise, repo-local, and named for this project.

## File Structure

- Modify `variantcentrifuge/transcript_filter.py`: new focused module for reading transcript IDs, filtering SnpEff `ANN`, writing filtered VCFs, and composing the transform for chunk and stage use.
- Modify `variantcentrifuge/stages/processing_stages.py`: add `TranscriptFilterStage`, wire it before `FieldExtractionStage`, and call the same helper inside `ParallelCompleteProcessingStage._process_single_chunk`.
- Modify `variantcentrifuge/pipeline.py`: add transcript filter stage when transcript options are present, and fix `--split-snpeff-lines` wiring to use `snpeff_splitting_mode`.
- Modify `variantcentrifuge/cli.py`: validate transcript input files early, normalize transcript IDs into config, and update help text for MANE usage.
- Modify `docs/source/usage.md` and `docs/source/guides/annotation_strategies.md`: document defined-transcript extraction and MANE transcript-list workflow.
- Test `tests/unit/test_transcript_filter.py`: pure unit tests for ANN filtering and transcript ID loading.
- Test `tests/unit/stages/test_transcript_filter_stage.py`: stage behavior with mocked helper.
- Test `tests/unit/stages/test_processing_stages_critical.py`: update split-stage config expectations.
- Test `tests/integration/test_transcript_filtering_pipeline.py`: small mocked pipeline-level behavior if existing integration mocks make this cheap.

---

### Task 1: Create and Verify the Proper Conda Environment

**Files:**
- Read: `conda/environment.yml`
- No source modifications.

- [ ] **Step 1: Create the environment if absent**

Run:

```powershell
conda env list
conda env create -f conda\environment.yml
```

Expected if missing:

```text
Preparing transaction: done
Executing transaction: done
```

Expected if already created by another worker:

```text
variantcentrifuge        C:\Users\bernt\anaconda3\envs\variantcentrifuge
```

- [ ] **Step 2: Verify external tools**

Run:

```powershell
conda run -n variantcentrifuge python --version
conda run -n variantcentrifuge bcftools --version
conda run -n variantcentrifuge SnpSift -h
conda run -n variantcentrifuge snpEff -h
conda run -n variantcentrifuge pytest --version
```

Expected:

```text
Python 3.12.x
bcftools 1.21
SnpSift version output/help
snpEff version output/help
pytest 8.x or newer
```

- [ ] **Step 3: Verify a narrow existing test before changes**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\test_bcftools_extractor.py tests\unit\stages\test_processing_stages_critical.py -q
```

Expected:

```text
passed
```

---

### Task 2: Add Pure ANN Transcript Filtering Helpers

**Files:**
- Create: `variantcentrifuge/transcript_filter.py`
- Test: `tests/unit/test_transcript_filter.py`

- [ ] **Step 1: Write failing unit tests**

Create `tests/unit/test_transcript_filter.py`:

```python
from pathlib import Path

import pytest

from variantcentrifuge.transcript_filter import (
    filter_ann_value_by_transcripts,
    load_transcript_ids,
)


def test_load_transcript_ids_from_list_and_file(tmp_path: Path) -> None:
    transcript_file = tmp_path / "transcripts.txt"
    transcript_file.write_text(
        "# comment\nNM_000059.4\n\nENST00000357654.9\n",
        encoding="utf-8",
    )

    result = load_transcript_ids("NM_007294.4, NM_000546.6", str(transcript_file))

    assert result == {
        "NM_007294.4",
        "NM_000546.6",
        "NM_000059.4",
        "ENST00000357654.9",
    }


def test_load_transcript_ids_rejects_empty_file(tmp_path: Path) -> None:
    transcript_file = tmp_path / "empty.txt"
    transcript_file.write_text("# only comments\n\n", encoding="utf-8")

    with pytest.raises(ValueError, match="No transcript IDs"):
        load_transcript_ids(None, str(transcript_file))


def test_filter_ann_keeps_matching_feature_id_and_reorders_to_single_ann() -> None:
    ann = (
        "A|synonymous_variant|LOW|BRCA2|BRCA2|transcript|NM_other.1|protein_coding||||||||,"
        "A|missense_variant|MODERATE|BRCA2|BRCA2|transcript|NM_000059.4|protein_coding|"
        "10/27|c.1234G>A|p.Arg412His|||||,"
        "A|intron_variant|MODIFIER|BRCA2|BRCA2|transcript|ENST00000380152.8|protein_coding||||||||"
    )

    result = filter_ann_value_by_transcripts(ann, {"NM_000059.4"})

    assert result == (
        "A|missense_variant|MODERATE|BRCA2|BRCA2|transcript|NM_000059.4|protein_coding|"
        "10/27|c.1234G>A|p.Arg412His|||||"
    )


def test_filter_ann_keeps_all_matching_feature_ids_in_original_order() -> None:
    ann = (
        "A|missense_variant|MODERATE|GENE|GENE|transcript|NM_1|protein_coding||||||||,"
        "A|splice_region_variant|LOW|GENE|GENE|transcript|NM_2|protein_coding||||||||,"
        "A|intron_variant|MODIFIER|GENE|GENE|transcript|NM_3|protein_coding||||||||"
    )

    result = filter_ann_value_by_transcripts(ann, {"NM_1", "NM_3"})

    assert result == (
        "A|missense_variant|MODERATE|GENE|GENE|transcript|NM_1|protein_coding||||||||,"
        "A|intron_variant|MODIFIER|GENE|GENE|transcript|NM_3|protein_coding||||||||"
    )


def test_filter_ann_returns_none_when_no_transcript_matches() -> None:
    ann = "A|missense_variant|MODERATE|BRCA2|BRCA2|transcript|NM_other.1|protein_coding||||||||"

    assert filter_ann_value_by_transcripts(ann, {"NM_000059.4"}) is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\test_transcript_filter.py -q
```

Expected:

```text
ModuleNotFoundError: No module named 'variantcentrifuge.transcript_filter'
```

- [ ] **Step 3: Implement pure helpers**

Create `variantcentrifuge/transcript_filter.py`:

```python
"""Transcript-level filtering for SnpEff ANN annotations."""

from __future__ import annotations

import gzip
import io
import logging
import os
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

ANN_FEATURE_ID_INDEX = 6


def load_transcript_ids(
    transcript_list: str | None = None,
    transcript_file: str | None = None,
) -> set[str]:
    """Load transcript IDs from a comma-separated list and/or one-per-line file."""
    transcript_ids: set[str] = set()

    if transcript_list:
        transcript_ids.update(
            item.strip() for item in transcript_list.split(",") if item.strip()
        )

    if transcript_file:
        path = Path(transcript_file)
        if not path.exists():
            raise FileNotFoundError(f"Transcript file not found: {transcript_file}")
        with path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line and not line.startswith("#"):
                    transcript_ids.add(line)

    if (transcript_list or transcript_file) and not transcript_ids:
        raise ValueError("No transcript IDs were provided")

    return transcript_ids


def filter_ann_value_by_transcripts(ann_value: str, transcript_ids: set[str]) -> str | None:
    """Return ANN entries whose Feature_ID matches the requested transcripts."""
    if not ann_value or ann_value == ".":
        return None

    kept_entries: list[str] = []
    for entry in ann_value.split(","):
        parts = entry.split("|")
        if len(parts) > ANN_FEATURE_ID_INDEX and parts[ANN_FEATURE_ID_INDEX] in transcript_ids:
            kept_entries.append(entry)

    if not kept_entries:
        return None

    return ",".join(kept_entries)
```

- [ ] **Step 4: Run tests to verify helper behavior**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\test_transcript_filter.py -q
```

Expected:

```text
5 passed
```

---

### Task 3: Add VCF Record Transform and Writer

**Files:**
- Modify: `variantcentrifuge/transcript_filter.py`
- Test: `tests/unit/test_transcript_filter.py`

- [ ] **Step 1: Add failing tests for INFO/ANN and VCF filtering**

Append to `tests/unit/test_transcript_filter.py`:

```python
import gzip

from variantcentrifuge.transcript_filter import (
    filter_vcf_line_by_transcripts,
    filter_vcf_to_transcripts,
)


def test_filter_vcf_line_rewrites_ann_and_preserves_other_info() -> None:
    line = (
        "1\t100\t.\tA\tG\t.\tPASS\t"
        "AC=1;ANN=G|synonymous_variant|LOW|GENE|GENE|transcript|NM_other.1|protein_coding||||||||,"
        "G|missense_variant|MODERATE|GENE|GENE|transcript|NM_keep.1|protein_coding||||||||;DP=20"
        "\tGT\t0/1\n"
    )

    result = filter_vcf_line_by_transcripts(line, {"NM_keep.1"})

    assert result is not None
    assert "AC=1;" in result
    assert "DP=20" in result
    assert "NM_keep.1" in result
    assert "NM_other.1" not in result


def test_filter_vcf_line_drops_no_match_record() -> None:
    line = (
        "1\t100\t.\tA\tG\t.\tPASS\t"
        "ANN=G|synonymous_variant|LOW|GENE|GENE|transcript|NM_other.1|protein_coding||||||||"
        "\tGT\t0/1\n"
    )

    assert filter_vcf_line_by_transcripts(line, {"NM_keep.1"}) is None


def test_filter_vcf_to_transcripts_writes_only_matching_records(tmp_path: Path) -> None:
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf.gz"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t100\t.\tA\tG\t.\tPASS\tANN=G|missense_variant|MODERATE|GENE|GENE|transcript|NM_keep.1|protein_coding||||||||\tGT\t0/1\n"
        "1\t200\t.\tC\tT\t.\tPASS\tANN=T|synonymous_variant|LOW|GENE|GENE|transcript|NM_other.1|protein_coding||||||||\tGT\t0/1\n",
        encoding="utf-8",
    )

    count = filter_vcf_to_transcripts(str(input_vcf), str(output_vcf), {"NM_keep.1"})

    assert count == 1
    with gzip.open(output_vcf, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert "#CHROM" in text
    assert "NM_keep.1" in text
    assert "NM_other.1" not in text
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\test_transcript_filter.py -q
```

Expected:

```text
ImportError: cannot import name 'filter_vcf_line_by_transcripts'
```

- [ ] **Step 3: Implement VCF transform**

Append to `variantcentrifuge/transcript_filter.py`:

```python
def filter_vcf_line_by_transcripts(line: str, transcript_ids: set[str]) -> str | None:
    """Filter one VCF record line by ANN Feature_ID."""
    if line.startswith("#"):
        return line

    columns = line.rstrip("\n").split("\t")
    if len(columns) < 8:
        return line

    info_parts = columns[7].split(";") if columns[7] else []
    rewritten_parts: list[str] = []
    found_ann = False

    for info_part in info_parts:
        if info_part.startswith("ANN="):
            found_ann = True
            filtered_ann = filter_ann_value_by_transcripts(
                info_part.removeprefix("ANN="),
                transcript_ids,
            )
            if filtered_ann is None:
                return None
            rewritten_parts.append(f"ANN={filtered_ann}")
        else:
            rewritten_parts.append(info_part)

    if not found_ann:
        return None

    columns[7] = ";".join(rewritten_parts)
    return "\t".join(columns) + "\n"


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def filter_vcf_to_transcripts(
    input_vcf: str,
    output_vcf: str,
    transcript_ids: set[str],
    *,
    index_output: bool = False,
) -> int:
    """Write a VCF containing only records with matching ANN transcript IDs."""
    kept_records = 0

    if output_vcf.endswith(".gz"):
        output_handle = gzip.open(output_vcf, "wt", encoding="utf-8", compresslevel=1)
    else:
        output_handle = open(output_vcf, "w", encoding="utf-8")

    with _open_text(input_vcf) as input_handle, output_handle:
        for line in input_handle:
            filtered = filter_vcf_line_by_transcripts(line, transcript_ids)
            if filtered is None:
                continue
            if not filtered.startswith("#"):
                kept_records += 1
            output_handle.write(filtered)

    if index_output and output_vcf.endswith(".gz"):
        subprocess.run(["bcftools", "index", "-f", output_vcf], check=True)

    logger.info(
        "Transcript filtering retained %d records for %d transcript IDs",
        kept_records,
        len(transcript_ids),
    )
    return kept_records
```

- [ ] **Step 4: Run pure tests**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\test_transcript_filter.py -q
```

Expected:

```text
8 passed
```

---

### Task 4: Wire Transcript Filtering Into Single-Thread Pipeline

**Files:**
- Modify: `variantcentrifuge/stages/processing_stages.py`
- Modify: `variantcentrifuge/pipeline.py`
- Test: `tests/unit/stages/test_transcript_filter_stage.py`

- [ ] **Step 1: Write failing stage tests**

Create `tests/unit/stages/test_transcript_filter_stage.py`:

```python
from pathlib import Path
from unittest.mock import patch

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.processing_stages import TranscriptFilterStage


@pytest.fixture
def context(tmp_path: Path) -> PipelineContext:
    workspace = Workspace(output_dir=tmp_path, base_name="sample")
    return PipelineContext(
        args=None,
        config={},
        workspace=workspace,
        data=None,
    )


def test_transcript_filter_skips_without_requested_transcripts(context: PipelineContext) -> None:
    input_vcf = Path("input.vcf.gz")
    context.extracted_vcf = input_vcf
    context.data = input_vcf
    context.config = {}

    result = TranscriptFilterStage()._process(context)

    assert result.data == input_vcf
    assert not hasattr(result, "transcript_filtered_vcf")


@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
def test_transcript_filter_uses_extracted_vcf_and_sets_context(
    mock_filter,
    context: PipelineContext,
    tmp_path: Path,
) -> None:
    input_vcf = tmp_path / "input.vcf.gz"
    input_vcf.write_text("dummy", encoding="utf-8")
    context.extracted_vcf = input_vcf
    context.data = input_vcf
    context.config = {"transcript_list": "NM_000059.4,NM_007294.4"}

    result = TranscriptFilterStage()._process(context)

    mock_filter.assert_called_once()
    assert mock_filter.call_args.args[0] == str(input_vcf)
    assert mock_filter.call_args.args[2] == {"NM_000059.4", "NM_007294.4"}
    assert result.transcript_filtered_vcf.name == "sample.transcript_filtered.vcf.gz"
    assert result.data == result.transcript_filtered_vcf
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\stages\test_transcript_filter_stage.py -q
```

Expected:

```text
ImportError: cannot import name 'TranscriptFilterStage'
```

- [ ] **Step 3: Implement `TranscriptFilterStage`**

Modify imports near the top of `variantcentrifuge/stages/processing_stages.py`:

```python
from ..transcript_filter import filter_vcf_to_transcripts, load_transcript_ids
```

Add this class after `SnpSiftFilterStage`:

```python
class TranscriptFilterStage(Stage):
    """Filter SnpEff ANN annotations to requested transcript IDs."""

    @property
    def name(self) -> str:
        return "transcript_filtering"

    @property
    def description(self) -> str:
        return "Filter ANN annotations by transcript ID"

    @property
    def dependencies(self) -> set[str]:
        return {"variant_extraction"}

    @property
    def soft_dependencies(self) -> set[str]:
        return {"snpsift_filtering", "multiallelic_split"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        transcript_ids = load_transcript_ids(
            context.config.get("transcript_list"),
            context.config.get("transcript_file"),
        )
        if not transcript_ids:
            logger.debug("No transcript filtering requested")
            return context

        input_vcf = (
            getattr(context, "split_annotations_vcf", None)
            or getattr(context, "filtered_vcf", None)
            or getattr(context, "extracted_vcf", None)
            or context.data
        )
        if not input_vcf:
            raise ValueError("No input VCF available for transcript filtering")

        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.transcript_filtered.vcf.gz"
        )
        filter_vcf_to_transcripts(
            str(input_vcf),
            str(output_vcf),
            transcript_ids,
            index_output=True,
        )

        context.transcript_filtered_vcf = output_vcf  # type: ignore[attr-defined]
        if not getattr(context, "extracted_tsv", None):
            context.data = output_vcf
        return context
```

- [ ] **Step 4: Wire the stage in `pipeline.py`**

In `variantcentrifuge/pipeline.py`, import `TranscriptFilterStage` with the other processing stages if needed.

Replace the current split gate:

```python
if hasattr(args, "snpeff_split_by_transcript") and args.snpeff_split_by_transcript:
    stages.append(MultiAllelicSplitStage())
```

with:

```python
if config.get("snpeff_splitting_mode"):
    stages.append(MultiAllelicSplitStage())

if config.get("transcript_list") or config.get("transcript_file"):
    stages.append(TranscriptFilterStage())
```

- [ ] **Step 5: Run stage tests**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\stages\test_transcript_filter_stage.py -q
```

Expected:

```text
2 passed
```

---

### Task 5: Fix `--split-snpeff-lines` Config Drift

**Files:**
- Modify: `variantcentrifuge/stages/processing_stages.py`
- Modify: `tests/unit/stages/test_processing_stages_critical.py`

- [ ] **Step 1: Update failing tests to current CLI config**

In `tests/unit/stages/test_processing_stages_critical.py`, change split-stage configs from:

```python
base_context.config = {
    "snpeff_split_by_transcript": True,
    "snpeff_split_before_filter": True,
}
```

to:

```python
base_context.config = {
    "snpeff_splitting_mode": "before_filters",
}
```

For the after-filter test, use:

```python
base_context.config = {
    "snpeff_splitting_mode": "after_filters",
}
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\stages\test_processing_stages_critical.py::TestMultiAllelicSplitStage -q
```

Expected:

```text
FAILED ... split not called
```

- [ ] **Step 3: Fix `MultiAllelicSplitStage`**

In `variantcentrifuge/stages/processing_stages.py`, replace:

```python
if not context.config.get("snpeff_split_by_transcript"):
    logger.debug("SNPeff transcript splitting not requested")
    return context

# Determine when to split
split_before = context.config.get("snpeff_split_before_filter", True)
if split_before and context.is_complete("snpsift_filtering"):
    # Already filtered, don't split again
    return context
```

with:

```python
splitting_mode = context.config.get("snpeff_splitting_mode")
if not splitting_mode:
    logger.debug("SNPeff transcript splitting not requested")
    return context

if splitting_mode == "before_filters" and context.is_complete("snpsift_filtering"):
    return context
```

- [ ] **Step 4: Run split-stage tests**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\stages\test_processing_stages_critical.py::TestMultiAllelicSplitStage -q
```

Expected:

```text
passed
```

---

### Task 6: Wire Transcript Filtering Into Parallel Chunk Processing

**Files:**
- Modify: `variantcentrifuge/stages/processing_stages.py`
- Test: `tests/unit/stages/test_parallel_processing_subtasks.py` or create `tests/unit/stages/test_parallel_transcript_filter.py`

- [ ] **Step 1: Add a focused test with mocks**

Create `tests/unit/stages/test_parallel_transcript_filter.py`:

```python
from pathlib import Path
from unittest.mock import patch

from variantcentrifuge.stages.processing_stages import ParallelCompleteProcessingStage


@patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
@patch("variantcentrifuge.stages.processing_stages.extract_variants")
def test_parallel_chunk_applies_transcript_filter_before_field_extraction(
    mock_extract_variants,
    mock_filter_transcripts,
    mock_extract_fields,
    tmp_path: Path,
) -> None:
    chunk_bed = tmp_path / "chunk.bed"
    chunk_bed.write_text("1\t0\t100\n", encoding="utf-8")

    stage = ParallelCompleteProcessingStage()
    result = stage._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file="input.vcf.gz",
        base_name="sample",
        intermediate_dir=tmp_path,
        config={
            "threads_per_chunk": 1,
            "fields_to_extract": "CHROM POS ANN[0].FEATUREID GEN[*].GT",
            "transcript_list": "NM_keep.1",
            "vcf_samples": ["S1"],
        },
    )

    transcript_vcf = tmp_path / "sample.chunk_0.transcript_filtered.vcf.gz"
    mock_filter_transcripts.assert_called_once()
    assert mock_filter_transcripts.call_args.args[2] == {"NM_keep.1"}
    assert mock_extract_fields.call_args.kwargs["variant_file"] == str(transcript_vcf)
    assert result.name == "sample.chunk_0.extracted.tsv.gz"
```

- [ ] **Step 2: Run test to verify failure**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\stages\test_parallel_transcript_filter.py -q
```

Expected:

```text
FAILED ... Expected 'filter_vcf_to_transcripts' to have been called once
```

- [ ] **Step 3: Implement parallel hook**

In `_process_single_chunk`, after filtering and before field extraction, insert:

```python
transcript_ids = load_transcript_ids(
    config.get("transcript_list"),
    config.get("transcript_file"),
)
if transcript_ids:
    transcript_filtered_vcf = intermediate_dir / f"{chunk_base}.transcript_filtered.vcf.gz"
    filter_vcf_to_transcripts(
        str(filtered_vcf),
        str(transcript_filtered_vcf),
        transcript_ids,
        index_output=True,
    )
    filtered_vcf = transcript_filtered_vcf
```

- [ ] **Step 4: Run parallel transcript test**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\stages\test_parallel_transcript_filter.py -q
```

Expected:

```text
1 passed
```

---

### Task 7: CLI Validation and Help Text

**Files:**
- Modify: `variantcentrifuge/cli.py`
- Test: `tests/test_cli_argument_parsing.py`

- [ ] **Step 1: Add parser/help tests**

Append to `tests/test_cli_argument_parsing.py`:

```python
def test_transcript_filtering_parameters_are_documented():
    from variantcentrifuge.cli import create_parser

    parser = create_parser()
    help_text = parser.format_help()

    assert "--transcript-list" in help_text
    assert "--transcript-file" in help_text
    assert "MANE" in help_text
```

- [ ] **Step 2: Run test to verify failure**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\test_cli_argument_parsing.py::test_transcript_filtering_parameters_are_documented -q
```

Expected:

```text
FAILED ... assert 'MANE' in help_text
```

- [ ] **Step 3: Update help text**

In `variantcentrifuge/cli.py`, replace the `--transcript-file` help with:

```python
help=(
    "Path to a file containing transcript IDs, one per line. "
    "For GRCh38, use MANE Select transcript lists matching the VCF namespace "
    "(RefSeq NM_* or Ensembl ENST*). For hg19/GRCh37 RefSeq annotations, use "
    "a GRCh37 RefSeq Select transcript list because MANE is defined on GRCh38."
),
```

- [ ] **Step 4: Add early validation after config load**

After `validate_vcf_file(args.vcf_file, logger)` in `main()`, add:

```python
if args.transcript_file and not Path(args.transcript_file).is_file():
    logger.error(f"Transcript file not found: {args.transcript_file}")
    sys.exit(1)
```

- [ ] **Step 5: Run parser tests**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\test_cli_argument_parsing.py -q
```

Expected:

```text
passed
```

---

### Task 8: Documentation for Defined Transcripts and MANE

**Files:**
- Modify: `docs/source/usage.md`
- Modify: `docs/source/guides/annotation_strategies.md`

- [ ] **Step 1: Update usage table**

In `docs/source/usage.md`, replace the transcript rows with:

```markdown
| `--transcript-list IDS` | — | Comma-separated transcript IDs. Matching `ANN` entries are retained and become `ANN[0]` for extraction |
| `--transcript-file PATH` | — | File with transcript IDs, one per line. Recommended for MANE Select workflows when the file uses the same transcript ID namespace as the VCF |
```

- [ ] **Step 2: Add MANE workflow section**

In `docs/source/guides/annotation_strategies.md`, add before "Custom Gene Annotations with VariantCentrifuge":

```markdown
## Transcript-Specific Extraction and MANE

VariantCentrifuge can restrict SnpEff `ANN` extraction to defined transcript IDs:

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file annotated.vcf.gz \
  --transcript-file mane_refseq_transcripts.txt \
  --fields "CHROM POS REF ALT ANN[0].GENE ANN[0].FEATUREID ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P GEN[*].GT" \
  --output-file mane_only.tsv
```

The transcript file must match the transcript namespace present in `ANN[0].FEATUREID`.
For SnpEff RefSeq annotations this is typically `NM_*`; for Ensembl/GENCODE annotations
this is typically `ENST*`. Transcript versions matter: `NM_000059.4` and
`NM_000059.3` are different IDs.

### MANE transcript lists

For GRCh38, use the official NCBI MANE summary file. It contains both RefSeq and
Ensembl transcript accessions:

```bash
curl -L -o MANE.GRCh38.v1.5.summary.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz

# RefSeq IDs for SnpEff RefSeq annotations
zcat MANE.GRCh38.v1.5.summary.txt.gz \
  | awk -F'\t' 'NR > 1 && $10 == "MANE Select" {print $6}' \
  > mane_select_refseq_transcripts.txt

# Ensembl IDs for Ensembl/GENCODE annotations
zcat MANE.GRCh38.v1.5.summary.txt.gz \
  | awk -F'\t' 'NR > 1 && $10 == "MANE Select" {print $8}' \
  > mane_select_ensembl_transcripts.txt
```

For GRCh37/hg19, there is no separate MANE release equivalent to the GRCh38 MANE
summary. MANE is defined on GRCh38, although NCBI provides mappings to GRCh37 for
clinical interpretation. For GRCh37 RefSeq-annotated VCFs, use the GRCh37 RefSeq
Select set instead:

```bash
curl -L -o GCF_000001405.25_GRCh37.p13_genomic.gff.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_000001405.25-RS_2024_09/GCF_000001405.25_GRCh37.p13_genomic.gff.gz

zcat GCF_000001405.25_GRCh37.p13_genomic.gff.gz \
  | awk -F'\t' '$3 == "mRNA" && $9 ~ /tag=RefSeq Select/ {
      match($9, /transcript_id=([^;]+)/, a);
      if (a[1] != "") print a[1]
    }' \
  | sort -u \
  > grch37_refseq_select_transcripts.txt
```

SnpEff also supports MANE-aware upstream annotation for GRCh38 with MANE-specific
databases such as `GRCh38.mane.1.2.refseq` and `GRCh38.mane.1.2.ensembl`, and can
restrict annotation with `-onlyTr` or `-tag MANE_Select`. VariantCentrifuge does not
infer MANE status from gene symbols alone; provide a transcript list that matches the
VCF's annotation source.
```

- [ ] **Step 3: Run docs-sensitive smoke checks**

Run:

```powershell
conda run -n variantcentrifuge pytest tests\unit\test_vcf_header_parser.py tests\test_cli_argument_parsing.py -q
```

Expected:

```text
passed
```

---

### Task 9: End-to-End Narrow Verification

**Files:**
- No new files unless a fixture is needed.

- [ ] **Step 1: Run unit tests for changed areas**

Run:

```powershell
conda run -n variantcentrifuge pytest `
  tests\unit\test_transcript_filter.py `
  tests\unit\stages\test_transcript_filter_stage.py `
  tests\unit\stages\test_parallel_transcript_filter.py `
  tests\unit\stages\test_processing_stages_critical.py `
  tests\unit\test_bcftools_extractor.py `
  tests\test_cli_argument_parsing.py `
  -q
```

Expected:

```text
passed
```

- [ ] **Step 2: Run formatting/linting on touched files**

Run:

```powershell
conda run -n variantcentrifuge ruff format variantcentrifuge\transcript_filter.py variantcentrifuge\stages\processing_stages.py variantcentrifuge\pipeline.py variantcentrifuge\cli.py tests\unit\test_transcript_filter.py tests\unit\stages\test_transcript_filter_stage.py tests\unit\stages\test_parallel_transcript_filter.py
conda run -n variantcentrifuge ruff check variantcentrifuge\transcript_filter.py variantcentrifuge\stages\processing_stages.py variantcentrifuge\pipeline.py variantcentrifuge\cli.py tests\unit\test_transcript_filter.py tests\unit\stages\test_transcript_filter_stage.py tests\unit\stages\test_parallel_transcript_filter.py
```

Expected:

```text
All checks passed!
```

- [ ] **Step 3: Run broader fast tests if time allows**

Run:

```powershell
conda run -n variantcentrifuge pytest -m "unit and not slow" -q
```

Expected:

```text
passed
```

---

## Implementation Notes

- Do not implement automatic MANE inference in this first patch. SnpEff supports MANE-aware GRCh38 annotation upstream, but VariantCentrifuge should consume explicit transcript lists rather than infer MANE status from gene symbols.
- Filtering `ANN` before extraction is preferable to filtering only VCF records, because existing `ANN[0].*` extraction takes the first comma-separated ANN entry.
- Keep all transcript filtering logic in `variantcentrifuge/transcript_filter.py`; pipeline stages should only load IDs, call the transform, and update context paths.
- Use `ANN` Feature_ID index `6`, matching existing `ANN_SUBFIELD_MAP["ANN[0].FEATUREID"]`.
- Preserve non-ANN INFO fields and all FORMAT/sample columns exactly.
- If a record has no `ANN`, drop it when transcript filtering is requested. That is consistent with "extract only defined transcripts."

## Self-Review

- Spec coverage: the plan fixes existing `--transcript-list` / `--transcript-file`, preserves downstream `ANN[0]` extraction behavior, covers single-thread and parallel processing, fixes split config drift, and documents GRCh38 MANE plus GRCh37/hg19 RefSeq Select guidance.
- Placeholder scan: no `TBD`, vague "add tests", or undefined helper names remain.
- Type consistency: helper names are consistent across tasks: `load_transcript_ids`, `filter_ann_value_by_transcripts`, `filter_vcf_line_by_transcripts`, `filter_vcf_to_transcripts`, and `TranscriptFilterStage`.
