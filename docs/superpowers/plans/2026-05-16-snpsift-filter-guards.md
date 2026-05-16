# SnpSift Filter Guards Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix issue #104 by correcting malformed built-in SnpSift presets and making all filtering layers fail closed when requested filters cannot be applied.

**Architecture:** Repair the three malformed built-in preset expressions and remove their test exemption. Add a SnpSift-specific stderr diagnostic guard around `apply_snpsift_filter()` by extending `run_command()` to optionally return `CompletedProcess`. Change late/final TSV filtering from fail-open to fail-closed, then add regression tests using BMP2K-like issue #104 records.

**Tech Stack:** Python, pytest, pandas, SnpSift, bgzip, bcftools, existing VariantCentrifuge pipeline stages.

---

## Files

- Modify `variantcentrifuge/config.json`
  - Correct `moderate_and_high_prediction`, `high_or_lof_or_nmd`, and `high_or_pathogenic`.
- Modify `variantcentrifuge/utils.py`
  - Add optional `return_result` support to `run_command()`.
- Modify `variantcentrifuge/filters.py`
  - Add SnpSift stderr diagnostic detection.
  - Use it in `apply_snpsift_filter()`.
  - Make TSV and DataFrame filters raise on invalid expressions.
- Modify `tests/unit/test_field_profile.py`
  - Remove malformed-preset exceptions.
  - Check both `dbnsfp4` and `dbnsfp5`.
- Modify `tests/unit/test_filters.py` or create it if absent.
  - Add unit tests for SnpSift stderr classification and `run_command()` result return.
- Modify `tests/unit/stages/test_output_stages.py`
  - Update invalid final-filter behavior from fail-open to fail-closed.
- Modify `tests/test_final_filter_simple.py`
  - Add direct invalid-query fail-closed test.
- Create `tests/integration/test_snpsift_filter_guards.py`
  - Add real-SnpSift regression coverage for issue #104.

---

### Task 1: Prove Built-In Presets Are Malformed

**Files:**
- Modify: `tests/unit/test_field_profile.py`

- [ ] **Step 1: Write the failing tests**

Replace `test_all_presets_have_balanced_parens` with this stricter version:

```python
def test_all_resolved_presets_have_balanced_parens_for_each_profile():
    """All built-in presets must remain syntactically composable."""
    from variantcentrifuge.config import load_config

    for profile in ("dbnsfp4", "dbnsfp5"):
        cfg = load_config()
        cfg["field_profile"] = profile
        resolve_profile(cfg)

        for name, expr in cfg["presets"].items():
            open_count = expr.count("(")
            close_count = expr.count(")")
            assert open_count == close_count, (
                f"Preset '{name}' has unbalanced parentheses under {profile}: "
                f"{open_count} opening, {close_count} closing. Expression: {expr}"
            )
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
pytest tests/unit/test_field_profile.py::TestFieldProfileResolution::test_all_resolved_presets_have_balanced_parens_for_each_profile -q
```

Expected: FAIL naming `moderate_and_high_prediction`, `high_or_lof_or_nmd`, or `high_or_pathogenic`.

- [ ] **Step 3: Fix the preset expressions**

In `variantcentrifuge/config.json`, replace these three values:

```json
"moderate_and_high_prediction": "((ANN[ANY].IMPACT has 'MODERATE') & ((dbNSFP_REVEL_score >= 0.9) | (dbNSFP_CADD_phred >= 30)))",
"high_or_lof_or_nmd": "(((exists LOF[*].PERC) & (LOF[*].PERC > 0.9)) | ((exists NMD[*].PERC) & (NMD[*].PERC > 0.9)) | (ANN[ANY].IMPACT has 'HIGH'))",
"high_or_pathogenic": "((ANN[ANY].IMPACT has 'HIGH') | (((dbNSFP_clinvar_clnsig =~ '[Pp]athogenic') & !(dbNSFP_clinvar_clnsig =~ '[Cc]onflicting')) | ((ClinVar_CLNSIG =~ '[Pp]athogenic') & !(ClinVar_CLNSIG =~ '[Cc]onflicting'))))",
```

- [ ] **Step 4: Run test to verify it passes**

Run:

```bash
pytest tests/unit/test_field_profile.py::TestFieldProfileResolution::test_all_resolved_presets_have_balanced_parens_for_each_profile -q
```

Expected: PASS.

- [ ] **Step 5: Verify related field-profile tests**

Run:

```bash
pytest tests/unit/test_field_profile.py -q
```

Expected: PASS.

- [ ] **Step 6: Commit checkpoint**

```bash
git add variantcentrifuge/config.json tests/unit/test_field_profile.py
git commit -m "fix: correct malformed SnpSift presets"
```

---

### Task 2: Fail Closed On SnpSift Parser Diagnostics

**Files:**
- Modify: `variantcentrifuge/utils.py`
- Modify: `variantcentrifuge/filters.py`
- Test: `tests/unit/test_filters.py`

- [ ] **Step 1: Write failing unit tests for command result access and stderr classification**

Add these tests to `tests/unit/test_filters.py`:

```python
import subprocess
from unittest.mock import patch

import pytest

from variantcentrifuge.filters import (
    _raise_for_snpsift_filter_stderr,
    apply_snpsift_filter,
)
from variantcentrifuge.utils import run_command


def test_run_command_can_return_completed_process():
    result = run_command(["python", "-c", "print('ok')"], return_result=True)

    assert isinstance(result, subprocess.CompletedProcess)
    assert result.returncode == 0
    assert result.stdout.strip() == "ok"


@pytest.mark.parametrize(
    "stderr",
    [
        "line 1:126 mismatched input ')' expecting {'*'}",
        "line 1:293 missing ')' at '<EOF>'",
        "line 1:10 token recognition error at: '@'",
        "Error parsing expression",
        "Exception in thread \"main\" java.lang.RuntimeException: Cannot parse EffectType 'HIGH'",
    ],
)
def test_snpsift_filter_stderr_parser_diagnostics_raise(stderr):
    with pytest.raises(RuntimeError, match="SnpSift filter reported parser diagnostics"):
        _raise_for_snpsift_filter_stderr(stderr, "BAD_EXPR")


def test_snpsift_filter_stderr_allows_empty_or_nonfatal_warnings():
    _raise_for_snpsift_filter_stderr("", "EXPR")
    _raise_for_snpsift_filter_stderr("Some non-fatal note", "EXPR")


@patch("variantcentrifuge.filters.run_command")
@patch("variantcentrifuge.filters.os.remove")
@patch("variantcentrifuge.filters.os.path.exists", return_value=True)
@patch("variantcentrifuge.filters._snpsift_memory_flag", return_value="-Xmx1g")
def test_apply_snpsift_filter_raises_before_compressing_when_snpsift_reports_parser_error(
    mock_memory,
    mock_exists,
    mock_remove,
    mock_run_command,
    tmp_path,
):
    input_vcf = tmp_path / "input.vcf.gz"
    output_vcf = tmp_path / "output.vcf.gz"
    input_vcf.write_text("placeholder")

    mock_run_command.return_value = subprocess.CompletedProcess(
        args=["SnpSift"],
        returncode=0,
        stdout="",
        stderr="line 1:126 mismatched input ')' expecting {'*'}",
    )

    with pytest.raises(RuntimeError, match="SnpSift filter reported parser diagnostics"):
        apply_snpsift_filter(str(input_vcf), "BAD_EXPR", {"threads": 1}, str(output_vcf))

    assert mock_run_command.call_count == 1
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_filters.py::test_run_command_can_return_completed_process tests/unit/test_filters.py::test_snpsift_filter_stderr_parser_diagnostics_raise tests/unit/test_filters.py::test_snpsift_filter_stderr_allows_empty_or_nonfatal_warnings tests/unit/test_filters.py::test_apply_snpsift_filter_raises_before_compressing_when_snpsift_reports_parser_error -q
```

Expected: FAIL because `run_command()` has no `return_result` argument and the SnpSift stderr helper does not exist.

- [ ] **Step 3: Add optional result return to `run_command()`**

In `variantcentrifuge/utils.py`, update the signature and return logic:

```python
def run_command(
    cmd: list, output_file: str | None = None, return_result: bool = False
) -> str | subprocess.CompletedProcess[str]:
```

After the nonzero return-code check, return the process when requested:

```python
    logger.debug("Command completed successfully.")
    if return_result:
        return result
    if output_file:
        return output_file
    return result.stdout
```

Keep existing default behavior unchanged.

- [ ] **Step 4: Add SnpSift stderr guard**

In `variantcentrifuge/filters.py`, add near `_snpsift_memory_flag()`:

```python
SNPSIFT_FATAL_STDERR_PATTERNS = (
    "mismatched input",
    "missing ')'",
    "token recognition error",
    "Error parsing",
    "Cannot parse",
    'Exception in thread "main"',
)


def _raise_for_snpsift_filter_stderr(stderr: str | None, filter_string: str) -> None:
    """Raise when SnpSift reports parser/evaluation diagnostics despite exit code 0."""
    if not stderr:
        return

    diagnostic_lines = [
        line.strip()
        for line in stderr.splitlines()
        if any(pattern in line for pattern in SNPSIFT_FATAL_STDERR_PATTERNS)
    ]
    if not diagnostic_lines:
        return

    preview = "\n".join(diagnostic_lines[:5])
    raise RuntimeError(
        "SnpSift filter reported parser diagnostics while exiting successfully. "
        f"Filter expression: {filter_string}\n{preview}"
    )
```

- [ ] **Step 5: Use the guard in `apply_snpsift_filter()`**

Replace:

```python
run_command(snpsift_cmd, output_file=tmp_vcf)
```

with:

```python
result = run_command(snpsift_cmd, output_file=tmp_vcf, return_result=True)
_raise_for_snpsift_filter_stderr(result.stderr, filter_string)
```

- [ ] **Step 6: Run focused tests to verify they pass**

Run:

```bash
pytest tests/unit/test_filters.py::test_run_command_can_return_completed_process tests/unit/test_filters.py::test_snpsift_filter_stderr_parser_diagnostics_raise tests/unit/test_filters.py::test_snpsift_filter_stderr_allows_empty_or_nonfatal_warnings tests/unit/test_filters.py::test_apply_snpsift_filter_raises_before_compressing_when_snpsift_reports_parser_error -q
```

Expected: PASS.

- [ ] **Step 7: Run filter tests**

Run:

```bash
pytest tests/unit/test_filters.py tests/unit/stages/test_processing_stages_critical.py -q
```

Expected: PASS.

- [ ] **Step 8: Commit checkpoint**

```bash
git add variantcentrifuge/utils.py variantcentrifuge/filters.py tests/unit/test_filters.py
git commit -m "fix: fail closed on SnpSift parser diagnostics"
```

---

### Task 3: Add Real SnpSift Regression For Issue #104

**Files:**
- Create: `tests/integration/test_snpsift_filter_guards.py`

- [ ] **Step 1: Write the failing integration test**

Create `tests/integration/test_snpsift_filter_guards.py`:

```python
"""Regression tests for SnpSift filter fail-closed behavior."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from variantcentrifuge.config import load_config
from variantcentrifuge.field_profile import resolve_profile
from variantcentrifuge.filters import apply_snpsift_filter


pytestmark = pytest.mark.skipif(
    not all(shutil.which(tool) for tool in ("SnpSift", "bgzip", "bcftools")),
    reason="SnpSift, bgzip, and bcftools are required",
)


def _write_issue_104_vcf(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO">',
                '##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects">',
                '##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects">',
                '##INFO=<ID=dbNSFP_gnomAD4.1_joint_AF,Number=A,Type=Float,Description="AF">',
                '##INFO=<ID=dbNSFP_gnomAD4.1_joint_AC,Number=A,Type=Integer,Description="AC">',
                '##INFO=<ID=dbNSFP_clinvar_clnsig,Number=.,Type=String,Description="ClinVar">',
                '##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="ClinVar">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr4\t78912028\t.\tC\tT\t500\tPASS\tANN=T|stop_gained|HIGH|BMP2K|ENSG00000138756|transcript|ENST00000264889|protein_coding|1/10|c.3481C>T|p.Gln1161Ter|1/100|1/100|1/33||;dbNSFP_gnomAD4.1_joint_AF=0.00221454;dbNSFP_gnomAD4.1_joint_AC=3546;dbNSFP_clinvar_clnsig=Benign/Likely_benign;ClinVar_CLNSIG=Benign/Likely_benign",
                "chr4\t78913000\t.\tG\tA\t500\tPASS\tANN=A|stop_gained|HIGH|KEEPER|ENSGKEEP|transcript|ENSTKEEP|protein_coding|1/10|c.1G>A|p.Trp1Ter|1/100|1/100|1/33||;dbNSFP_gnomAD4.1_joint_AF=0.00001;dbNSFP_gnomAD4.1_joint_AC=1;dbNSFP_clinvar_clnsig=Pathogenic;ClinVar_CLNSIG=Pathogenic",
                "",
            ]
        )
    )


def _dbnsfp5_filter(*presets: str) -> str:
    cfg = load_config()
    cfg["field_profile"] = "dbnsfp5"
    resolve_profile(cfg)
    return " & ".join(f"({cfg['presets'][name]})" for name in presets)


def test_issue_104_high_common_benign_record_is_filtered_out(tmp_path):
    input_vcf = tmp_path / "issue104.vcf"
    output_vcf = tmp_path / "filtered.vcf.gz"
    _write_issue_104_vcf(input_vcf)

    filter_expr = _dbnsfp5_filter("high_or_lof_or_nmd", "rare", "not_benign")

    apply_snpsift_filter(str(input_vcf), filter_expr, {"threads": 1}, str(output_vcf))

    result = subprocess.run(
        ["bcftools", "view", "-H", str(output_vcf)],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = [line for line in result.stdout.splitlines() if line]
    assert len(rows) == 1
    assert "KEEPER" in rows[0]
    assert "BMP2K" not in rows[0]
```

- [ ] **Step 2: Run test to verify it fails before Task 1 or Task 2 fixes**

Run:

```bash
pytest tests/integration/test_snpsift_filter_guards.py::test_issue_104_high_common_benign_record_is_filtered_out -q
```

Expected before fixes: FAIL because BMP2K leaks or because SnpSift parser diagnostics now raise before preset repair. Expected after Tasks 1 and 2: PASS.

- [ ] **Step 3: Run integration test after fixes**

Run:

```bash
pytest tests/integration/test_snpsift_filter_guards.py -q
```

Expected: PASS, or SKIP if required external tools are unavailable.

- [ ] **Step 4: Commit checkpoint**

```bash
git add tests/integration/test_snpsift_filter_guards.py
git commit -m "test: cover issue 104 SnpSift filter leakage"
```

---

### Task 4: Make Late And Final Filters Fail Closed

**Files:**
- Modify: `variantcentrifuge/filters.py`
- Modify: `tests/test_final_filter_simple.py`
- Modify: `tests/unit/stages/test_output_stages.py`

- [ ] **Step 1: Write failing direct filter tests**

Add to `tests/test_final_filter_simple.py`:

```python
import pytest
```

and add:

```python
def test_final_filter_invalid_expression_raises():
    df = pd.DataFrame({"CHROM": ["chr1"], "score": ["0.1"]})

    with pytest.raises(ValueError, match="Invalid final filter expression"):
        filter_dataframe_with_query(df, 'CHROM === "chr1"')
```

- [ ] **Step 2: Update output-stage invalid filter test to expect fail-closed**

In `tests/unit/stages/test_output_stages.py`, replace `test_invalid_filter_expression` with:

```python
def test_invalid_filter_expression_raises(self, context):
    """Invalid final filters must fail closed instead of returning unfiltered data."""
    context.config["final_filter"] = 'CHROM === "chr1"'

    stage = FinalFilteringStage()

    with pytest.raises(ValueError, match="Invalid final filter expression"):
        stage(context)
```

- [ ] **Step 3: Run tests to verify they fail**

Run:

```bash
pytest tests/test_final_filter_simple.py::test_final_filter_invalid_expression_raises tests/unit/stages/test_output_stages.py::TestFinalFilteringStage::test_invalid_filter_expression_raises -q
```

Expected: FAIL because invalid filters currently return unfiltered data.

- [ ] **Step 4: Change DataFrame filtering to raise**

In `variantcentrifuge/filters.py`, replace the `except Exception as e` block in `filter_dataframe_with_query()` with:

```python
    except Exception as e:
        message = f"Invalid final filter expression '{filter_expression}': {e}"
        logger.error(message)
        raise ValueError(message) from e
```

- [ ] **Step 5: Change TSV filtering to raise**

In `filter_tsv_with_expression()`, replace the inner query exception handler with:

```python
        except Exception as e:
            message = f"Invalid TSV filter expression '{filter_expression}': {e}"
            logger.error(message)
            raise ValueError(message) from e
```

Replace the outer exception fallback that copies input to output with:

```python
    except Exception as e:
        if isinstance(e, ValueError):
            raise
        message = f"Error in TSV filtering: {e}"
        logger.error(message)
        raise RuntimeError(message) from e
```

- [ ] **Step 6: Run focused tests to verify they pass**

Run:

```bash
pytest tests/test_final_filter_simple.py::test_final_filter_invalid_expression_raises tests/unit/stages/test_output_stages.py::TestFinalFilteringStage::test_invalid_filter_expression_raises -q
```

Expected: PASS.

- [ ] **Step 7: Run output and final-filter tests**

Run:

```bash
pytest tests/test_final_filter_simple.py tests/unit/stages/test_output_stages.py -q
```

Expected: PASS.

- [ ] **Step 8: Commit checkpoint**

```bash
git add variantcentrifuge/filters.py tests/test_final_filter_simple.py tests/unit/stages/test_output_stages.py
git commit -m "fix: fail closed on invalid final filters"
```

---

### Task 5: Full Verification And Cleanup

**Files:**
- Review all modified files.

- [ ] **Step 1: Run focused filter suite**

Run:

```bash
pytest tests/unit/test_field_profile.py tests/unit/test_filters.py tests/test_final_filter_simple.py tests/unit/stages/test_processing_stages_critical.py tests/unit/stages/test_output_stages.py tests/integration/test_snpsift_filter_guards.py -q
```

Expected: PASS, with `test_snpsift_filter_guards.py` skipped only if external tools are unavailable.

- [ ] **Step 2: Run formatting**

Run:

```bash
make format
```

Expected: command completes successfully.

- [ ] **Step 3: Run lint**

Run:

```bash
make lint
```

Expected: PASS.

- [ ] **Step 4: Run typecheck**

Run:

```bash
make typecheck
```

Expected: PASS.

- [ ] **Step 5: Run format check**

Run:

```bash
make format-check
```

Expected: PASS.

- [ ] **Step 6: Run CI check**

Run:

```bash
make ci-check
```

Expected: PASS.

- [ ] **Step 7: Commit any formatting-only changes**

If `make format` changed files:

```bash
git status --short
git add <formatted-files>
git commit -m "style: format SnpSift filter guard changes"
```

If no files changed, skip this step.

- [ ] **Step 8: Final branch status**

Run:

```bash
git status --short --branch
git log --oneline --decorate -5
```

Expected: clean working tree on the issue #104 branch, with checkpoint commits present.

---

## Self-Review

- Issue #104 examples are covered by Task 3 with a BMP2K-like common/benign HIGH-impact record and a rare/pathogenic positive control.
- The three known malformed presets are covered by Task 1.
- SnpSift exit-code `0` with parser diagnostics is covered by Task 2.
- Broader fail-open filter behavior in late/final filtering is covered by Task 4.
- Each task has red/green TDD steps and checkpoint commits.
- No production implementation step appears before a failing test step in its task.
