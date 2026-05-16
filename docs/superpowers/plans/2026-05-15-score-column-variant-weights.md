# Score Column Variant Weights Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement GitHub issue #97 by allowing arbitrary numeric variant-level score columns to drive burden and Python SKAT variant weights while preserving existing Beta, uniform, Fisher, and COAST behavior.

**Architecture:** Extend `variantcentrifuge.association.weights` so one shared resolver returns concrete per-variant weight vectors for burden tests and Python SKAT. Move annotation alignment to the lazy genotype matrix builder by propagating the real genotype keep mask, then pass aligned score/CADD/REVEL/effect arrays through the engine to tests. Python SKAT resolves weights in `PurePythonSKATTest.run()` with its current post-imputation MAF source and passes the concrete vector to the backend.

**Tech Stack:** Python, numpy, pandas, scipy, statsmodels, argparse, pytest, existing VariantCentrifuge association stage and SKAT backend.

---

## Verified Current Code Paths

- `variantcentrifuge/association/weights.py`
  - `get_weights()` starts at line 305 and currently supports only `beta:*`, `uniform`, `cadd`, `revel`, and `combined`.
  - CADD/REVEL parsing already uses `_parse_scores_to_float()` and missing-score logging helpers.
- `variantcentrifuge/association/tests/logistic_burden.py`
  - `LogisticBurdenTest.run()` calls `get_weights()` at lines 332-340 and already passes CADD, REVEL, effects, and params.
- `variantcentrifuge/association/tests/linear_burden.py`
  - `LinearBurdenTest.run()` calls `get_weights()` at lines 154-162 with the same functional-weight inputs as logistic burden.
- `variantcentrifuge/association/tests/skat_python.py`
  - `PurePythonSKATTest.run()` currently calls `parse_weights_beta()` at lines 266-274 and passes only `weights_beta` to `PythonSKATBackend.test_gene()`.
  - ACAT-V currently recomputes Beta weights at lines 277-289 from `geno.mean(axis=0) / 2.0`.
- `variantcentrifuge/association/backends/python_backend.py`
  - `PythonSKATBackend.test_gene()` takes `weights_beta` at lines 655-662.
  - `_test_skat()`, `_test_burden()`, and `_test_skato()` each recompute Beta weights from `geno.mean(axis=0) / 2.0`.
- `variantcentrifuge/association/genotype_matrix.py`
  - `build_genotype_matrix()` computes `keep_variants_mask` at line 249 but currently returns only `(geno, mafs, sample_mask, warnings_list)`.
- `variantcentrifuge/stages/analysis_stages.py`
  - `_GenotypeMatrixBuilder.__call__()` builds the matrix lazily at lines 80-88.
  - The association stage currently extracts CADD/REVEL/effect arrays at lines 2721-2785 by duplicating the genotype missingness mask with a regex.
- `variantcentrifuge/association/engine.py`
  - The engine invokes lazy builders in the first-gene parallel path, remaining parallel path, and sequential path, then copies only genotype, MAF, phenotype, and covariate keys.
- `variantcentrifuge/association/tests/skat_r.py`
  - R SKAT still parses only Beta parameters at lines 288-296 and delegates to an R backend that uses `weights.beta`.

---

## File Structure

- Modify `variantcentrifuge/association/weights.py`
  - Own score-column weight semantics, config spec parsing, and `get_weights(..., score_values=...)` dispatch.
- Modify `variantcentrifuge/association/base.py`
  - Add `AssociationConfig.variant_weight_column`.
- Modify `variantcentrifuge/cli.py`
  - Add `--variant-weight-column`, parse only dict-valued `--variant-weight-params`, and validate simple CLI-only score-column errors.
- Modify `variantcentrifuge/stages/analysis_stages.py`
  - Accept `variant_weight_column` in JSON config, validate dict-valued params, resolve score columns, validate DataFrame columns, and pass annotation column names into `_GenotypeMatrixBuilder`.
- Modify `variantcentrifuge/association/genotype_matrix.py`
  - Add an opt-in return of the real `keep_variants_mask` without breaking existing four-tuple callers.
- Modify `variantcentrifuge/association/engine.py`
  - Copy all lazy builder annotation arrays into `gene_data` before tests and discard them after tests.
- Modify `variantcentrifuge/association/tests/logistic_burden.py`
  - Pass aligned `score_values` to `get_weights()`.
- Modify `variantcentrifuge/association/tests/linear_burden.py`
  - Pass aligned `score_values` to `get_weights()`.
- Modify `variantcentrifuge/association/backends/base.py`
  - Update the SKAT backend contract to accept concrete weight vectors.
- Modify `variantcentrifuge/association/backends/python_backend.py`
  - Use provided concrete weights in SKAT, Burden, and SKAT-O backend methods.
- Modify `variantcentrifuge/association/tests/skat_python.py`
  - Resolve concrete weights with `get_weights()` using `geno.mean(axis=0) / 2.0`, then pass the same vector to SKAT and ACAT-V.
- Modify `variantcentrifuge/association/tests/skat_r.py`
  - Fail clearly for non-Beta and non-uniform R SKAT weights.
- Modify `variantcentrifuge/association/engine.py`
  - Fail before R dependency checks when R SKAT is selected with unsupported weight specs.
- Create `tests/unit/test_score_column_weights.py`
  - Unit tests for score-column semantics and dispatch.
- Modify `tests/unit/test_functional_weights.py`
  - Regression coverage proving existing specs still dispatch unchanged after the supported-spec message changes.
- Modify `tests/unit/test_json_config.py`
  - JSON association validation and `AssociationConfig` propagation tests.
- Modify `tests/unit/test_cli_association_args.py`
  - Direct `_build_assoc_config_from_context()` propagation tests for `variant_weight_column`.
- Modify `tests/test_cli_argument_parsing.py`
  - Parser/main tests for CLI propagation and invalid weight params.
- Modify `tests/unit/test_association_genotype_matrix.py`
  - Keep-mask opt-in tests.
- Modify `tests/unit/test_streaming_matrix.py`
  - Lazy builder annotation alignment and engine handoff tests.
- Create `tests/unit/test_score_column_association_stage.py`
  - Stage-level score-column validation and builder wiring tests.
- Modify `tests/unit/test_association_logistic_burden.py`
  - Logistic burden score-column regression tests.
- Modify `tests/unit/test_association_linear_burden.py`
  - Linear burden raw score-only regression tests.
- Modify `tests/unit/test_skat_python_backend.py`
  - Backend tests that explicit weights are used by SKAT, Burden, and SKAT-O paths.
- Create `tests/unit/test_skat_python_weights.py`
  - Pure Python SKAT wrapper tests for shared resolver use, CADD/REVEL/combined behavior, ACAT-V weights, and post-imputation MAF compatibility.
- Modify `tests/unit/test_skat_r_test.py`
  - R SKAT unsupported-weight guard tests.
- Modify `tests/unit/test_association_engine.py`
  - Engine-level fail-fast test before R dependency checks.
- Modify `tests/unit/test_association_fisher.py`
  - Fisher ignores variant-weight config regression.
- Modify `tests/unit/test_coast_python.py`
  - COAST remains independent from variant-weight config regression.
- Create `tests/unit/test_score_column_association_integration.py`
  - Small DataFrame regression for `logistic_burden,skat` with score-column weights.
- Modify `docs/source/guides/association_testing.md`
  - User docs for score-column weights, params, examples, Fisher/COAST notes, Python SKAT behavior change, and R SKAT limitation.
- Modify `docs/source/configuration.md`
  - Top-level cross-reference to the association testing guide for score-column variant weights.

---

### Task 1: Add Config And CLI Plumbing

**Files:**
- Modify: `variantcentrifuge/association/base.py`
- Modify: `variantcentrifuge/cli.py`
- Modify: `variantcentrifuge/stages/analysis_stages.py`
- Test: `tests/unit/test_json_config.py`
- Test: `tests/unit/test_cli_association_args.py`
- Test: `tests/test_cli_argument_parsing.py`

- [ ] **Step 1: Write failing JSON config tests**

Add these tests to `tests/unit/test_json_config.py`:

```python
def test_variant_weight_column_validated_and_propagated():
    ctx = _make_context(
        {
            "association": {
                "variant_weights": "score_column",
                "variant_weight_column": "nephro_candidate_score",
                "variant_weight_params": {
                    "score_min": 0,
                    "score_max": 10,
                    "floor": 0.1,
                    "combine_with_beta": True,
                },
            }
        }
    )

    config = _build_assoc_config_from_context(ctx)

    assert config.variant_weights == "score_column"
    assert config.variant_weight_column == "nephro_candidate_score"
    assert config.variant_weight_params == {
        "score_min": 0,
        "score_max": 10,
        "floor": 0.1,
        "combine_with_beta": True,
    }


def test_variant_weight_column_must_be_string_when_present():
    with pytest.raises(ValueError, match="'variant_weight_column' must be a string"):
        _validate_association_config_dict({"variant_weight_column": ["ncs"]})


def test_variant_weight_params_must_be_dict_when_present():
    with pytest.raises(ValueError, match="'variant_weight_params' must be an object"):
        _validate_association_config_dict({"variant_weight_params": ["score_min", 0]})
```

- [ ] **Step 2: Write failing direct config propagation tests**

Add these tests to `tests/unit/test_cli_association_args.py`:

```python
@pytest.mark.unit
def test_variant_weight_column_from_cli_config():
    ctx = _make_context(
        {
            "variant_weights": "score_column",
            "variant_weight_column": "nephro_candidate_score",
        }
    )

    assoc_config = _build_assoc_config_from_context(ctx)

    assert assoc_config.variant_weights == "score_column"
    assert assoc_config.variant_weight_column == "nephro_candidate_score"


@pytest.mark.unit
def test_variant_weight_column_json_fallback():
    ctx = _make_context(
        {
            "association": {
                "variant_weights": "score_column",
                "variant_weight_column": "json_score",
            }
        }
    )

    assoc_config = _build_assoc_config_from_context(ctx)

    assert assoc_config.variant_weights == "score_column"
    assert assoc_config.variant_weight_column == "json_score"


@pytest.mark.unit
def test_cli_variant_weight_column_overrides_json():
    ctx = _make_context(
        {
            "variant_weight_column": "cli_score",
            "association": {
                "variant_weights": "score_column",
                "variant_weight_column": "json_score",
            },
        }
    )

    assoc_config = _build_assoc_config_from_context(ctx)

    assert assoc_config.variant_weight_column == "cli_score"
```

- [ ] **Step 3: Write failing CLI parser tests**

Add `import pytest` near the top of `tests/test_cli_argument_parsing.py`, then add these tests:

```python
def test_variant_weight_column_cli_populates_config(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    captured_configs = []

    def fake_run_pipeline(args):
        captured_configs.append(args.config.copy())

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    monkeypatch.setattr(cli_module, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--output-dir",
            str(tmp_path),
            "--perform-association",
            "--association-tests",
            "logistic_burden",
            "--variant-weights",
            "score_column",
            "--variant-weight-column",
            "nephro_candidate_score",
            "--variant-weight-params",
            '{"score_min":0,"score_max":10,"floor":0.1}',
        ],
    )

    exit_code = cli_module.main()

    assert exit_code == 0
    assert captured_configs[0]["variant_weights"] == "score_column"
    assert captured_configs[0]["variant_weight_column"] == "nephro_candidate_score"
    assert captured_configs[0]["variant_weight_params"] == {
        "score_min": 0,
        "score_max": 10,
        "floor": 0.1,
    }


def test_variant_weight_params_cli_requires_json_object(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--output-dir",
            str(tmp_path),
            "--perform-association",
            "--association-tests",
            "skat",
            "--variant-weight-params",
            '["score_min",0]',
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli_module.main()

    assert exc_info.value.code == 2


def test_score_column_cli_requires_variant_weight_column(monkeypatch, tmp_path):
    from variantcentrifuge import cli as cli_module

    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    monkeypatch.setattr(
        "sys.argv",
        [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf),
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--output-dir",
            str(tmp_path),
            "--perform-association",
            "--association-tests",
            "skat",
            "--variant-weights",
            "score_column",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli_module.main()

    assert exc_info.value.code == 2
```

- [ ] **Step 4: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_json_config.py::test_variant_weight_column_validated_and_propagated tests/unit/test_json_config.py::test_variant_weight_column_must_be_string_when_present tests/unit/test_json_config.py::test_variant_weight_params_must_be_dict_when_present -q
pytest tests/unit/test_cli_association_args.py::test_variant_weight_column_from_cli_config tests/unit/test_cli_association_args.py::test_variant_weight_column_json_fallback tests/unit/test_cli_association_args.py::test_cli_variant_weight_column_overrides_json -q
pytest tests/test_cli_argument_parsing.py::test_variant_weight_column_cli_populates_config tests/test_cli_argument_parsing.py::test_variant_weight_params_cli_requires_json_object tests/test_cli_argument_parsing.py::test_score_column_cli_requires_variant_weight_column -q
```

Expected:

- JSON tests fail because `variant_weight_column` is an unknown key or absent from `AssociationConfig`.
- CLI propagation tests fail because `AssociationConfig` has no `variant_weight_column`.
- CLI parser tests fail because `--variant-weight-column` is unrecognized and JSON arrays are accepted as `variant_weight_params`.

- [ ] **Step 5: Implement config and CLI plumbing**

In `variantcentrifuge/association/base.py`, add the field after `variant_weight_params`:

```python
variant_weight_column: str | None = None
"""Column name used when variant_weights='score_column'. None for inline column:<name> and non-score schemes."""
```

In `variantcentrifuge/cli.py`, add this parser argument immediately after `--variant-weights`:

```python
stats_group.add_argument(
    "--variant-weight-column",
    type=str,
    default=None,
    help=(
        "Variant-level numeric score column to use when --variant-weights score_column is set. "
        "Ignored when --variant-weights uses inline column:<name>."
    ),
)
```

In `variantcentrifuge/cli.py`, set the config key immediately after `cfg["variant_weights"]`:

```python
cfg["variant_weights"] = getattr(args, "variant_weights", "beta:1,25")
cfg["variant_weight_column"] = getattr(args, "variant_weight_column", None)
```

In the `--variant-weight-params` JSON parsing block in `variantcentrifuge/cli.py`, replace the assignment with:

```python
parsed_params = _json.loads(_vwp_raw)
if not isinstance(parsed_params, dict):
    raise TypeError(
        f"expected a JSON object, got {type(parsed_params).__name__}"
    )
cfg["variant_weight_params"] = parsed_params
```

In `variantcentrifuge/cli.py`, add these association argument validations near the existing association validation block:

```python
if getattr(args, "variant_weight_column", None) and not args.perform_association:
    parser.error("--variant-weight-column requires --perform-association to be set")

if args.perform_association:
    variant_weights_arg = getattr(args, "variant_weights", "beta:1,25")
    if variant_weights_arg == "score_column" and not getattr(args, "variant_weight_column", None):
        parser.error("--variant-weight-column is required when --variant-weights score_column")
    if variant_weights_arg.startswith("column:") and not variant_weights_arg[len("column:") :].strip():
        parser.error("--variant-weights column:<name> requires a non-empty column name")
```

In `variantcentrifuge/stages/analysis_stages.py`, add `"variant_weight_column"` to `VALID_ASSOCIATION_KEYS`.

In `_validate_association_config_dict()`, add `"variant_weight_column"` to `str_keys` and add a dict-validation branch:

```python
dict_keys = {"variant_weight_params"}

for key in dict_keys & set(d):
    if d[key] is not None and not isinstance(d[key], dict):
        errors.append(f"'{key}' must be an object, got {type(d[key]).__name__}")
```

In `_build_assoc_config_from_context()`, pass the new field into `AssociationConfig` immediately after `variant_weight_params`:

```python
variant_weight_column=_get("variant_weight_column", default=None, nullable=True),
```

- [ ] **Step 6: Run tests to verify they pass**

Run:

```bash
pytest tests/unit/test_json_config.py::test_variant_weight_column_validated_and_propagated tests/unit/test_json_config.py::test_variant_weight_column_must_be_string_when_present tests/unit/test_json_config.py::test_variant_weight_params_must_be_dict_when_present -q
pytest tests/unit/test_cli_association_args.py::test_variant_weight_column_from_cli_config tests/unit/test_cli_association_args.py::test_variant_weight_column_json_fallback tests/unit/test_cli_association_args.py::test_cli_variant_weight_column_overrides_json -q
pytest tests/test_cli_argument_parsing.py::test_variant_weight_column_cli_populates_config tests/test_cli_argument_parsing.py::test_variant_weight_params_cli_requires_json_object tests/test_cli_argument_parsing.py::test_score_column_cli_requires_variant_weight_column -q
```

Expected: all selected tests pass.

- [ ] **Step 7: Commit config and CLI plumbing**

```bash
git add variantcentrifuge/association/base.py variantcentrifuge/cli.py variantcentrifuge/stages/analysis_stages.py tests/unit/test_json_config.py tests/unit/test_cli_association_args.py tests/test_cli_argument_parsing.py
git commit -m "feat: add score-column weight config plumbing"
```

---

### Task 2: Add Score-Column Weight Resolver

**Files:**
- Modify: `variantcentrifuge/association/weights.py`
- Create: `tests/unit/test_score_column_weights.py`
- Modify: `tests/unit/test_functional_weights.py`

- [ ] **Step 1: Write failing score-column weight tests**

Create `tests/unit/test_score_column_weights.py`:

```python
"""Tests for arbitrary score-column variant weights."""

from __future__ import annotations

import logging

import numpy as np
import pytest

from variantcentrifuge.association.weights import (
    beta_maf_weights,
    get_weights,
    resolve_score_weight_column,
    score_column_weights,
)


@pytest.fixture
def mafs() -> np.ndarray:
    return np.array([0.01, 0.02, 0.05], dtype=np.float64)


def test_resolve_score_weight_column_from_inline_column():
    assert resolve_score_weight_column("column:nephro_candidate_score", None) == (
        "nephro_candidate_score"
    )


def test_resolve_score_weight_column_inline_wins_over_separate(caplog):
    caplog.set_level(logging.DEBUG, logger="variantcentrifuge")

    result = resolve_score_weight_column("column:inline_score", "separate_score")

    assert result == "inline_score"
    assert "separate_score" in caplog.text


def test_resolve_score_weight_column_from_friendly_form():
    assert resolve_score_weight_column("score_column", "nephro_candidate_score") == (
        "nephro_candidate_score"
    )


def test_resolve_score_weight_column_requires_column_for_friendly_form():
    with pytest.raises(ValueError, match="score_column requires variant_weight_column"):
        resolve_score_weight_column("score_column", None)


def test_resolve_score_weight_column_rejects_empty_inline_column():
    with pytest.raises(ValueError, match="column:<name> requires a non-empty column name"):
        resolve_score_weight_column("column:", None)


def test_score_column_default_combines_beta_with_unit_scale_scores(mafs):
    scores = np.array([0.2, 0.5, 1.0], dtype=np.float64)

    result = score_column_weights(mafs, scores)

    expected = beta_maf_weights(mafs, a=1.0, b=25.0) * scores
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_score_column_explicit_range_floor_and_beta_params(mafs):
    scores = np.array([0.0, 5.0, 10.0], dtype=np.float64)

    result = score_column_weights(
        mafs,
        scores,
        weight_params={
            "score_min": 0,
            "score_max": 10,
            "floor": 0.1,
            "beta_a": 0.5,
            "beta_b": 0.5,
        },
    )

    functional = np.array([0.1, 0.5, 1.0], dtype=np.float64)
    expected = beta_maf_weights(mafs, a=0.5, b=0.5) * functional
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_score_column_raw_score_only_mode(mafs):
    scores = np.array([0.2, 0.5, 2.0], dtype=np.float64)

    result = score_column_weights(
        mafs,
        scores,
        weight_params={"combine_with_beta": False, "ceiling": 1.0},
    )

    np.testing.assert_allclose(result, np.array([0.2, 0.5, 1.0]), rtol=1e-12)


def test_score_column_missing_default_is_neutral_when_combined_with_beta(mafs):
    scores = np.array([0.2, np.nan, "."], dtype=object)

    result = score_column_weights(mafs, scores)

    expected = beta_maf_weights(mafs) * np.array([0.2, 1.0, 1.0])
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_score_column_missing_default_is_floor_in_raw_mode(mafs):
    scores = np.array([0.2, np.nan, "."], dtype=object)

    result = score_column_weights(
        mafs,
        scores,
        weight_params={"combine_with_beta": False, "floor": 0.1},
    )

    np.testing.assert_allclose(result, np.array([0.2, 0.1, 0.1]), rtol=1e-12)


def test_score_column_missing_neutral_invalid_in_raw_mode(mafs):
    scores = np.array([0.2, np.nan, 0.5], dtype=object)

    with pytest.raises(ValueError, match="missing='neutral'.*combine_with_beta=false"):
        score_column_weights(
            mafs,
            scores,
            weight_params={"combine_with_beta": False, "missing": "neutral"},
        )


@pytest.mark.parametrize(
    "params, message",
    [
        ({"score_min": 10, "score_max": 0}, "score_min must be less than score_max"),
        ({"floor": -0.1}, "floor must be >= 0"),
        ({"ceiling": 0}, "ceiling must be > 0"),
        ({"floor": 0.8, "ceiling": 0.5}, "floor must be <= ceiling"),
        ({"beta_a": 0}, "beta_a must be > 0"),
        ({"beta_b": -1}, "beta_b must be > 0"),
        ({"missing": "mean"}, "missing must be one of"),
    ],
)
def test_score_column_invalid_params_raise(mafs, params, message):
    with pytest.raises(ValueError, match=message):
        score_column_weights(mafs, np.array([0.1, 0.2, 0.3]), weight_params=params)


def test_score_column_length_mismatch_raises(mafs):
    with pytest.raises(ValueError, match="same length"):
        score_column_weights(mafs, np.array([0.1, 0.2]))


def test_get_weights_dispatches_inline_column(mafs):
    scores = np.array([0.2, 0.5, 1.0], dtype=np.float64)

    result = get_weights(mafs, "column:nephro_candidate_score", score_values=scores)

    expected = score_column_weights(mafs, scores)
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_get_weights_dispatches_score_column_alias(mafs):
    scores = np.array([0.2, 0.5, 1.0], dtype=np.float64)

    result = get_weights(mafs, "score_column", score_values=scores)

    expected = score_column_weights(mafs, scores)
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_get_weights_score_column_requires_score_values(mafs):
    with pytest.raises(ValueError, match="requires score_values"):
        get_weights(mafs, "score_column")


def test_get_weights_column_empty_name_raises(mafs):
    with pytest.raises(ValueError, match="column:<name> requires a non-empty column name"):
        get_weights(mafs, "column:", score_values=np.array([0.1, 0.2, 0.3]))
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_score_column_weights.py -q
```

Expected: import errors for `resolve_score_weight_column` and `score_column_weights`, followed by dispatch failures after those names exist.

- [ ] **Step 3: Implement score-column helper functions**

In `variantcentrifuge/association/weights.py`, update the module docstring supported list to include:

```python
- ``score_column_weights``: Beta(MAF) x arbitrary numeric score-column weights,
  or raw normalized score weights when combine_with_beta is false.
```

Add this helper after `_parse_scores_to_float()`:

```python
def resolve_score_weight_column(
    weight_spec: str,
    variant_weight_column: str | None = None,
) -> str | None:
    """Resolve a score-column weight spec to a DataFrame column name."""
    if weight_spec.startswith("column:"):
        inline_column = weight_spec[len("column:") :].strip()
        if not inline_column:
            raise ValueError("column:<name> requires a non-empty column name")
        if variant_weight_column:
            logger.debug(
                "variant_weights=%r includes inline column %r; ignoring variant_weight_column=%r",
                weight_spec,
                inline_column,
                variant_weight_column,
            )
        return inline_column

    if weight_spec == "score_column":
        if not variant_weight_column:
            raise ValueError("score_column requires variant_weight_column")
        return variant_weight_column

    return None
```

Add parameter normalization helpers before `score_column_weights()`:

```python
_SCORE_WEIGHT_DEFAULTS = {
    "score_min": None,
    "score_max": None,
    "floor": 0.0,
    "ceiling": 1.0,
    "combine_with_beta": True,
    "missing": None,
    "beta_a": 1.0,
    "beta_b": 25.0,
}


def _normalize_score_weight_params(weight_params: dict | None) -> dict:
    if weight_params is not None and not isinstance(weight_params, dict):
        raise ValueError(
            f"variant_weight_params must be a dict, got {type(weight_params).__name__}"
        )

    params = dict(_SCORE_WEIGHT_DEFAULTS)
    params.update(weight_params or {})

    score_min = params["score_min"]
    score_max = params["score_max"]
    if (score_min is None) != (score_max is None):
        raise ValueError("score_min and score_max must be provided together")
    if score_min is not None and float(score_min) >= float(score_max):
        raise ValueError("score_min must be less than score_max")

    floor = float(params["floor"])
    ceiling = float(params["ceiling"])
    if floor < 0:
        raise ValueError("floor must be >= 0")
    if ceiling <= 0:
        raise ValueError("ceiling must be > 0")
    if floor > ceiling:
        raise ValueError("floor must be <= ceiling")

    beta_a = float(params["beta_a"])
    beta_b = float(params["beta_b"])
    if beta_a <= 0:
        raise ValueError("beta_a must be > 0")
    if beta_b <= 0:
        raise ValueError("beta_b must be > 0")

    missing = params["missing"]
    if missing not in (None, "neutral", "floor"):
        raise ValueError("missing must be one of None, 'neutral', or 'floor'")

    combine_with_beta = bool(params["combine_with_beta"])
    if missing == "neutral" and not combine_with_beta:
        raise ValueError("missing='neutral' is invalid when combine_with_beta=false")

    params["floor"] = floor
    params["ceiling"] = ceiling
    params["beta_a"] = beta_a
    params["beta_b"] = beta_b
    params["combine_with_beta"] = combine_with_beta
    if score_min is not None:
        params["score_min"] = float(score_min)
        params["score_max"] = float(score_max)
    return params
```

Add the public resolver before `combined_weights()`:

```python
def score_column_weights(
    mafs: np.ndarray,
    score_values: np.ndarray,
    *,
    variant_effects: np.ndarray | None = None,
    weight_params: dict | None = None,
) -> np.ndarray:
    """Compute weights from an arbitrary numeric variant score column."""
    mafs_arr = np.asarray(mafs, dtype=np.float64)
    scores_f = _parse_scores_to_float(np.asarray(score_values, dtype=object))
    assert scores_f is not None

    if len(mafs_arr) != len(scores_f):
        raise ValueError(
            f"score_values and mafs must have the same length "
            f"(got {len(scores_f)} and {len(mafs_arr)})"
        )

    params = _normalize_score_weight_params(weight_params)
    finite_mask = np.isfinite(scores_f)
    missing_mask = ~finite_mask

    normalized = np.empty(len(scores_f), dtype=np.float64)
    normalized[:] = np.nan

    if params["score_min"] is not None:
        score_min = float(params["score_min"])
        score_max = float(params["score_max"])
        normalized[finite_mask] = (scores_f[finite_mask] - score_min) / (score_max - score_min)
    else:
        normalized[finite_mask] = scores_f[finite_mask]
        out_of_unit = finite_mask & ((scores_f < 0.0) | (scores_f > 1.0))
        if bool(out_of_unit.any()):
            logger.warning(
                "score_column weights: %d finite score(s) outside [0, 1] without explicit range",
                int(out_of_unit.sum()),
            )

    functional = np.clip(normalized, params["floor"], params["ceiling"])

    missing_mode = params["missing"]
    if missing_mode is None:
        missing_mode = "neutral" if params["combine_with_beta"] else "floor"

    if missing_mode == "neutral":
        functional[missing_mask] = 1.0
    else:
        functional[missing_mask] = params["floor"]

    _log_missing_score_counts(missing_mask, variant_effects, "score_column")

    if bool(missing_mask.all()) and len(missing_mask) > 0:
        logger.warning("score_column weights: all score values for this gene are missing or invalid")

    if params["combine_with_beta"]:
        maf_w = beta_maf_weights(mafs_arr, a=params["beta_a"], b=params["beta_b"])
        return np.asarray(maf_w * functional, dtype=np.float64)

    return np.asarray(functional, dtype=np.float64)
```

- [ ] **Step 4: Extend `get_weights()` dispatch**

In `variantcentrifuge/association/weights.py`, add the new keyword to the signature:

```python
score_values: np.ndarray | None = None,
```

Add this dispatch block after the `beta:*` branch and before functional CADD/REVEL branches:

```python
if weight_spec == "score_column" or weight_spec.startswith("column:"):
    if weight_spec.startswith("column:") and not weight_spec[len("column:") :].strip():
        raise ValueError("column:<name> requires a non-empty column name")
    if score_values is None:
        raise ValueError(
            f"weight_spec='{weight_spec}' requires score_values to be provided"
        )
    return score_column_weights(
        mafs_arr,
        score_values,
        variant_effects=variant_effects,
        weight_params=weight_params,
    )
```

Update the unknown-spec error text:

```python
"Supported specs: 'beta:a,b', 'uniform', 'cadd', 'revel', 'combined', "
"'column:<name>', 'score_column'."
```

- [ ] **Step 5: Run tests to verify score-column behavior passes**

Run:

```bash
pytest tests/unit/test_score_column_weights.py tests/unit/test_functional_weights.py -q
```

Expected: all tests pass, including existing functional-weight regression tests.

- [ ] **Step 6: Commit weight resolver**

```bash
git add variantcentrifuge/association/weights.py tests/unit/test_score_column_weights.py tests/unit/test_functional_weights.py
git commit -m "feat: add score-column variant weight resolver"
```

---

### Task 3: Expose Genotype Keep Mask Without Breaking Existing Callers

**Files:**
- Modify: `variantcentrifuge/association/genotype_matrix.py`
- Test: `tests/unit/test_association_genotype_matrix.py`

- [ ] **Step 1: Write failing keep-mask tests**

Add these tests to `tests/unit/test_association_genotype_matrix.py`:

```python
def test_build_genotype_matrix_can_return_keep_variants_mask() -> None:
    gt_values = [
        ["0/1", "0/1", "0/0", "0/0"],
        ["./.", "./.", "0/0", "0/0"],
        ["1/1", "0/1", "0/0", "0/0"],
    ]
    gene_df, vcf_samples, gt_cols = _make_gene_df(4, 3, gt_values)

    geno, mafs, sample_mask, warnings_list, keep_mask = build_genotype_matrix(
        gene_df,
        vcf_samples,
        gt_cols,
        missing_site_threshold=0.25,
        return_keep_mask=True,
    )

    assert geno.shape == (4, 2)
    assert len(mafs) == 2
    assert sample_mask == [True, True, True, True]
    assert warnings_list == []
    np.testing.assert_array_equal(keep_mask, np.array([True, False, True]))


def test_build_genotype_matrix_default_return_shape_stays_four_tuple() -> None:
    gt_values = [["0/1", "0/0"]]
    gene_df, vcf_samples, gt_cols = _make_gene_df(2, 1, gt_values)

    result = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

    assert len(result) == 4
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_association_genotype_matrix.py::test_build_genotype_matrix_can_return_keep_variants_mask tests/unit/test_association_genotype_matrix.py::test_build_genotype_matrix_default_return_shape_stays_four_tuple -q
```

Expected: the first test fails with `TypeError` because `return_keep_mask` is not accepted; the second test passes before the implementation and must still pass after it.

- [ ] **Step 3: Implement opt-in keep-mask return**

In `variantcentrifuge/association/genotype_matrix.py`, update the signature:

```python
def build_genotype_matrix(
    gene_df: pd.DataFrame,
    vcf_samples: Sequence[str],
    gt_columns: Sequence[str],
    is_binary: bool = True,
    missing_site_threshold: float = 0.10,
    missing_sample_threshold: float = 0.80,
    phenotype_vector: np.ndarray | None = None,
    return_keep_mask: bool = False,
) -> tuple[np.ndarray, np.ndarray, list[bool], list[str]] | tuple[
    np.ndarray, np.ndarray, list[bool], list[str], np.ndarray
]:
```

At the end of `build_genotype_matrix()`, replace the return statement with:

```python
if return_keep_mask:
    return geno, mafs, sample_mask, warnings_list, keep_variants_mask
return geno, mafs, sample_mask, warnings_list
```

Update the docstring parameter and returns sections to state that the fifth return value is present only when `return_keep_mask=True`.

- [ ] **Step 4: Run tests to verify they pass**

Run:

```bash
pytest tests/unit/test_association_genotype_matrix.py::test_build_genotype_matrix_can_return_keep_variants_mask tests/unit/test_association_genotype_matrix.py::test_build_genotype_matrix_default_return_shape_stays_four_tuple -q
```

Expected: both tests pass.

- [ ] **Step 5: Run genotype matrix regression tests**

Run:

```bash
pytest tests/unit/test_association_genotype_matrix.py -q
```

Expected: all genotype matrix tests pass, proving the default four-tuple behavior is preserved.

- [ ] **Step 6: Commit keep-mask support**

```bash
git add variantcentrifuge/association/genotype_matrix.py tests/unit/test_association_genotype_matrix.py
git commit -m "feat: expose genotype keep mask for lazy alignment"
```

---

### Task 4: Align Annotation Arrays Inside The Lazy Builder

**Files:**
- Modify: `variantcentrifuge/stages/analysis_stages.py`
- Test: `tests/unit/test_streaming_matrix.py`

- [ ] **Step 1: Write failing lazy builder alignment tests**

Add these tests to `tests/unit/test_streaming_matrix.py`:

```python
def test_builder_aligns_score_values_with_real_keep_mask():
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    df = pd.DataFrame(
        {
            "GENE": ["A", "A", "A"],
            "nephro_candidate_score": [1.0, 99.0, 7.0],
            "GEN_0__GT": ["0/1", "./.", "1/1"],
            "GEN_1__GT": ["0/1", "./.", "0/1"],
            "GEN_2__GT": ["0/1", "0/0", "0/1"],
            "GEN_3__GT": ["0/1", "0/0", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        score_column="nephro_candidate_score",
    )

    result = builder()

    assert result["genotype_matrix"].shape[1] == 2
    np.testing.assert_array_equal(result["score_values"], np.array([1.0, 7.0], dtype=object))
    assert result["variant_weight_column"] == "nephro_candidate_score"


def test_builder_aligns_functional_annotation_arrays_with_real_keep_mask():
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    df = pd.DataFrame(
        {
            "GENE": ["A", "A", "A"],
            "dbNSFP_CADD_phred": [30.0, 99.0, 12.0],
            "dbNSFP_REVEL_score": [0.9, 0.1, 0.4],
            "ANN_0__EFFECT": ["missense_variant", "stop_gained", "frameshift_variant"],
            "GEN_0__GT": ["0/1", "./.", "1/1"],
            "GEN_1__GT": ["0/1", "./.", "0/1"],
            "GEN_2__GT": ["0/1", "0/0", "0/1"],
            "GEN_3__GT": ["0/1", "0/0", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        cadd_column="dbNSFP_CADD_phred",
        revel_column="dbNSFP_REVEL_score",
        effect_column="ANN_0__EFFECT",
    )

    result = builder()

    np.testing.assert_array_equal(result["cadd_scores"], np.array([30.0, 12.0], dtype=object))
    np.testing.assert_array_equal(result["revel_scores"], np.array([0.9, 0.4], dtype=object))
    np.testing.assert_array_equal(
        result["variant_effects"],
        np.array(["missense_variant", "frameshift_variant"], dtype=object),
    )


def test_builder_uses_parser_keep_mask_for_malformed_genotypes():
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    df = pd.DataFrame(
        {
            "GENE": ["A", "A", "A"],
            "nephro_candidate_score": [1.0, 99.0, 7.0],
            "GEN_0__GT": ["0/1", "0/x", "1/1"],
            "GEN_1__GT": ["0/1", "0/x", "0/1"],
            "GEN_2__GT": ["0/1", "0/x", "0/1"],
            "GEN_3__GT": ["0/1", "0/x", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        score_column="nephro_candidate_score",
    )

    result = builder()

    assert result["genotype_matrix"].shape[1] == 2
    np.testing.assert_array_equal(result["score_values"], np.array([1.0, 7.0], dtype=object))
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_streaming_matrix.py::TestGenotypeMatrixBuilder::test_builder_result_keys tests/unit/test_streaming_matrix.py::test_builder_aligns_score_values_with_real_keep_mask tests/unit/test_streaming_matrix.py::test_builder_aligns_functional_annotation_arrays_with_real_keep_mask tests/unit/test_streaming_matrix.py::test_builder_uses_parser_keep_mask_for_malformed_genotypes -q
```

Expected:

- New tests fail because `_GenotypeMatrixBuilder` does not accept annotation column fields or return annotation arrays.
- Existing `test_builder_result_keys` continues to pass after implementation because optional annotation keys are only present when requested.

- [ ] **Step 3: Add annotation fields to `_GenotypeMatrixBuilder`**

In `variantcentrifuge/stages/analysis_stages.py`, add these fields at the end of `_GenotypeMatrixBuilder` so current constructor calls remain valid:

```python
score_column: str | None = None
cadd_column: str | None = None
revel_column: str | None = None
effect_column: str | None = None
```

Add this helper method inside `_GenotypeMatrixBuilder`:

```python
def _aligned_annotation_payload(self, keep_variants_mask: np.ndarray) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    if self.score_column is not None:
        payload["score_values"] = self.gene_df[self.score_column].to_numpy(dtype=object)[
            keep_variants_mask
        ]
        payload["variant_weight_column"] = self.score_column
    if self.cadd_column is not None:
        payload["cadd_scores"] = self.gene_df[self.cadd_column].to_numpy(dtype=object)[
            keep_variants_mask
        ]
    if self.revel_column is not None:
        payload["revel_scores"] = self.gene_df[self.revel_column].to_numpy(dtype=object)[
            keep_variants_mask
        ]
    if self.effect_column is not None:
        payload["variant_effects"] = self.gene_df[self.effect_column].to_numpy(dtype=object)[
            keep_variants_mask
        ]
    return payload
```

Add this helper method inside `_GenotypeMatrixBuilder`:

```python
def _empty_annotation_payload(self) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    if self.score_column is not None:
        payload["score_values"] = np.asarray([], dtype=object)
        payload["variant_weight_column"] = self.score_column
    if self.cadd_column is not None:
        payload["cadd_scores"] = np.asarray([], dtype=object)
    if self.revel_column is not None:
        payload["revel_scores"] = np.asarray([], dtype=object)
    if self.effect_column is not None:
        payload["variant_effects"] = np.asarray([], dtype=object)
    return payload
```

- [ ] **Step 4: Use real keep mask inside `_GenotypeMatrixBuilder.__call__()`**

In the empty `gene_df` return branch, build `result` and merge optional empty payload:

```python
result = {
    "genotype_matrix": np.zeros((n_samples, 0), dtype=float),
    "variant_mafs": np.zeros(0, dtype=float),
    "phenotype_vector": self.phenotype_vector,
    "covariate_matrix": self.covariate_matrix,
    "gt_warnings": [],
    "mac_filtered": False,
}
result.update(self._empty_annotation_payload())
return result
```

Change the genotype matrix call to request the keep mask:

```python
geno, mafs, sample_mask, gt_warnings, keep_variants_mask = build_genotype_matrix(
    self.gene_df,
    self.vcf_samples,
    self.gt_columns,
    is_binary=self.is_binary,
    missing_site_threshold=self.missing_site_threshold,
    missing_sample_threshold=self.missing_sample_threshold,
    phenotype_vector=self.phenotype_vector,
    return_keep_mask=True,
)
annotation_payload = self._aligned_annotation_payload(keep_variants_mask)
```

In the MAC filter branch, clear aligned annotation arrays:

```python
if total_mac < 5:
    geno = np.zeros((geno.shape[0], 0), dtype=float)
    mafs = np.zeros(0, dtype=float)
    annotation_payload = self._empty_annotation_payload()
    mac_filtered = True
```

Merge the annotation payload into the returned dict:

```python
result = {
    "genotype_matrix": geno,
    "variant_mafs": mafs,
    "phenotype_vector": pv,
    "covariate_matrix": cm,
    "gt_warnings": gt_warnings,
    "mac_filtered": mac_filtered,
}
result.update(annotation_payload)
return result
```

- [ ] **Step 5: Run tests to verify they pass**

Run:

```bash
pytest tests/unit/test_streaming_matrix.py::TestGenotypeMatrixBuilder::test_builder_result_keys tests/unit/test_streaming_matrix.py::test_builder_aligns_score_values_with_real_keep_mask tests/unit/test_streaming_matrix.py::test_builder_aligns_functional_annotation_arrays_with_real_keep_mask tests/unit/test_streaming_matrix.py::test_builder_uses_parser_keep_mask_for_malformed_genotypes -q
```

Expected: all selected tests pass.

- [ ] **Step 6: Commit lazy builder alignment**

```bash
git add variantcentrifuge/stages/analysis_stages.py tests/unit/test_streaming_matrix.py
git commit -m "feat: align variant annotations in lazy genotype builder"
```

---

### Task 5: Wire Stage-Level Score Column Resolution

**Files:**
- Modify: `variantcentrifuge/stages/analysis_stages.py`
- Create: `tests/unit/test_score_column_association_stage.py`

- [ ] **Step 1: Write failing stage tests**

Create `tests/unit/test_score_column_association_stage.py`:

```python
"""Stage-level tests for score-column variant weights."""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage


def _workspace(tmp_path):
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.intermediate_dir.mkdir()
    workspace.get_intermediate_path = lambda name: workspace.intermediate_dir / name
    workspace.get_output_path = lambda name, ext=".tsv": workspace.output_dir / f"{name}{ext}"
    workspace.base_name = "score_column_stage"
    return workspace


def _context(tmp_path, df, extra_config):
    cases = [f"CASE{i}" for i in range(10)]
    controls = [f"CTRL{i}" for i in range(10)]
    config = {
        "perform_association": True,
        "case_samples": cases,
        "control_samples": controls,
        "association_tests": ["logistic_burden"],
        "variant_weights": "score_column",
        "variant_weight_column": "nephro_candidate_score",
        "output_dir": str(tmp_path),
        "output_file_base": "score_column_stage",
        "gzip_intermediates": False,
    }
    config.update(extra_config)
    ctx = PipelineContext(args=Namespace(), config=config, workspace=_workspace(tmp_path))
    ctx.current_dataframe = df
    ctx.vcf_samples = cases + controls
    return ctx


def _df_with_score():
    sample_cols = {}
    for i in range(20):
        sample_cols[f"GEN_{i}__GT"] = ["0/1", "0/1", "0/0"]
    return pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1", "GENE1"],
            "nephro_candidate_score": [2.0, 5.0, 9.0],
            "dbNSFP_CADD_phred": [10.0, 20.0, 30.0],
            "dbNSFP_REVEL_score": [0.1, 0.2, 0.3],
            "ANN_0__EFFECT": ["missense_variant", "stop_gained", "frameshift_variant"],
            **sample_cols,
        }
    )


def test_stage_fails_when_score_column_missing(tmp_path):
    df = _df_with_score().drop(columns=["nephro_candidate_score"])
    context = _context(tmp_path, df, {})

    with pytest.raises(ValueError, match="variant weight column 'nephro_candidate_score'"):
        AssociationAnalysisStage()._process(context)


def test_stage_passes_score_column_to_lazy_builder(monkeypatch, tmp_path):
    from variantcentrifuge.association.engine import AssociationEngine

    captured_gene_data = []

    def fake_run_all(self, gene_burden_data):
        captured_gene_data.extend(gene_burden_data)
        return pd.DataFrame(
            {
                "gene": ["GENE1"],
                "n_cases": [10],
                "n_controls": [10],
                "n_variants": [3],
                "logistic_burden_pvalue": [1.0],
            }
        )

    monkeypatch.setattr(AssociationEngine, "run_all", fake_run_all)
    context = _context(tmp_path, _df_with_score(), {})

    AssociationAnalysisStage()._process(context)

    builder = captured_gene_data[0]["_genotype_matrix_builder"]
    assert builder.score_column == "nephro_candidate_score"
    assert builder.effect_column == "ANN_0__EFFECT"


def test_stage_inline_column_spec_wins_over_variant_weight_column(monkeypatch, tmp_path):
    from variantcentrifuge.association.engine import AssociationEngine

    captured_gene_data = []

    def fake_run_all(self, gene_burden_data):
        captured_gene_data.extend(gene_burden_data)
        return pd.DataFrame(
            {
                "gene": ["GENE1"],
                "n_cases": [10],
                "n_controls": [10],
                "n_variants": [3],
                "logistic_burden_pvalue": [1.0],
            }
        )

    monkeypatch.setattr(AssociationEngine, "run_all", fake_run_all)
    df = _df_with_score()
    df["inline_score"] = [1.0, 2.0, 3.0]
    context = _context(
        tmp_path,
        df,
        {
            "variant_weights": "column:inline_score",
            "variant_weight_column": "nephro_candidate_score",
        },
    )

    AssociationAnalysisStage()._process(context)

    builder = captured_gene_data[0]["_genotype_matrix_builder"]
    assert builder.score_column == "inline_score"
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_score_column_association_stage.py -q
```

Expected: tests fail because the stage does not resolve score columns or pass annotation column names to `_GenotypeMatrixBuilder`.

- [ ] **Step 3: Implement stage column detection helpers**

In `variantcentrifuge/stages/analysis_stages.py`, import the resolver near the existing association imports:

```python
from ..association.weights import resolve_score_weight_column
```

Add these module-level helpers near `_GenotypeMatrixBuilder`:

```python
def _find_cadd_column(df: pd.DataFrame) -> str | None:
    return next(
        (c for c in df.columns if c.lower() in ("dbnsfp_cadd_phred", "cadd_phred")),
        None,
    )


def _find_revel_column(df: pd.DataFrame) -> str | None:
    return next(
        (c for c in df.columns if c.lower() in ("dbnsfp_revel_score", "revel_score")),
        None,
    )


def _find_effect_column(df: pd.DataFrame) -> str | None:
    return next((c for c in df.columns if c.upper() in ("EFFECT", "ANN_0__EFFECT")), None)
```

- [ ] **Step 4: Resolve and validate score column before grouping**

In `AssociationAnalysisStage._process()`, after `assoc_config` is built and before gene grouping, add:

```python
score_weight_column = resolve_score_weight_column(
    assoc_config.variant_weights,
    assoc_config.variant_weight_column,
)
if score_weight_column is not None and score_weight_column not in df.columns:
    raise ValueError(
        f"Requested variant weight column '{score_weight_column}' is missing from "
        "the association input DataFrame."
    )

needs_functional_annotations = assoc_config.variant_weights in {
    "cadd",
    "revel",
    "combined",
}
```

- [ ] **Step 5: Pass annotation column names into the builder and remove duplicated regex mask**

Inside the per-gene loop where `_GenotypeMatrixBuilder` is created, compute columns before instantiating the builder:

```python
cadd_column = _find_cadd_column(gene_df) if needs_functional_annotations else None
revel_column = _find_revel_column(gene_df) if needs_functional_annotations else None
effect_column = (
    _find_effect_column(gene_df)
    if needs_functional_annotations or score_weight_column is not None
    else None
)
```

Pass them to `_GenotypeMatrixBuilder`:

```python
builder = _GenotypeMatrixBuilder(
    gene_df=gene_df,
    vcf_samples=vcf_samples_list,
    gt_columns=gt_columns_for_matrix,
    is_binary=is_binary,
    missing_site_threshold=assoc_config.missing_site_threshold,
    missing_sample_threshold=assoc_config.missing_sample_threshold,
    phenotype_vector=phenotype_vector,
    covariate_matrix=covariate_matrix,
    score_column=score_weight_column,
    cadd_column=cadd_column,
    revel_column=revel_column,
    effect_column=effect_column,
)
```

Delete the current CADD/REVEL extraction block that starts with:

```python
# Phase 23: Extract functional annotation columns for CADD/REVEL weight schemes
if assoc_config.variant_weights in ("cadd", "revel", "combined"):
```

and ends after the `variant_effects` assignment.

- [ ] **Step 6: Run tests to verify stage wiring passes**

Run:

```bash
pytest tests/unit/test_score_column_association_stage.py tests/unit/test_streaming_matrix.py -q
```

Expected: all selected tests pass.

- [ ] **Step 7: Commit stage wiring**

```bash
git add variantcentrifuge/stages/analysis_stages.py tests/unit/test_score_column_association_stage.py tests/unit/test_streaming_matrix.py
git commit -m "feat: wire score-column annotations through association stage"
```

---

### Task 6: Pass Lazy Builder Annotation Results Through The Engine

**Files:**
- Modify: `variantcentrifuge/association/engine.py`
- Test: `tests/unit/test_streaming_matrix.py`

- [ ] **Step 1: Write failing engine handoff test**

Add this test to `tests/unit/test_streaming_matrix.py`:

```python
def test_engine_passes_builder_score_values_to_tests_and_discards_payload():
    from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
    from variantcentrifuge.association.engine import AssociationEngine
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    seen: dict[str, object] = {}

    class RecordingTest(AssociationTest):
        parallel_safe = True

        @property
        def name(self) -> str:
            return "recording"

        def run(self, gene, contingency_data, config):
            seen["score_values"] = contingency_data.get("score_values")
            seen["variant_weight_column"] = contingency_data.get("variant_weight_column")
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=1.0,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=2,
                n_controls=2,
                n_variants=2,
            )

    df = pd.DataFrame(
        {
            "GENE": ["A", "A"],
            "nephro_candidate_score": [2.0, 8.0],
            "GEN_0__GT": ["0/1", "1/1"],
            "GEN_1__GT": ["0/1", "0/1"],
            "GEN_2__GT": ["0/0", "0/1"],
            "GEN_3__GT": ["0/0", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        score_column="nephro_candidate_score",
    )
    gene_data = [
        {
            "GENE": "A",
            "proband_count": 2,
            "control_count": 2,
            "n_qualifying_variants": 2,
            "_genotype_matrix_builder": builder,
        }
    ]

    engine = AssociationEngine([RecordingTest()], AssociationConfig())
    result_df = engine.run_all(gene_data)

    assert result_df is not None
    np.testing.assert_array_equal(seen["score_values"], np.array([2.0, 8.0], dtype=object))
    assert seen["variant_weight_column"] == "nephro_candidate_score"
    assert "score_values" not in gene_data[0]
    assert "variant_weight_column" not in gene_data[0]
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
pytest tests/unit/test_streaming_matrix.py::test_engine_passes_builder_score_values_to_tests_and_discards_payload -q
```

Expected: the test fails because `AssociationEngine` does not copy annotation keys from builder results into `gene_data`.

- [ ] **Step 3: Implement engine helper functions**

In `variantcentrifuge/association/engine.py`, add these module-level constants and helpers near `_run_gene_worker()`:

```python
_BUILDER_RESULT_KEYS = (
    "genotype_matrix",
    "variant_mafs",
    "phenotype_vector",
    "covariate_matrix",
    "score_values",
    "variant_weight_column",
    "cadd_scores",
    "revel_scores",
    "variant_effects",
)

_PER_GENE_MATRIX_PAYLOAD_KEYS = (
    "genotype_matrix",
    "variant_mafs",
    "score_values",
    "variant_weight_column",
    "cadd_scores",
    "revel_scores",
    "variant_effects",
)


def _apply_builder_result_to_gene_data(
    gene_data: dict[str, Any],
    result: dict[str, Any],
    gene: str,
) -> None:
    for key in _BUILDER_RESULT_KEYS:
        if key in result:
            gene_data[key] = result[key]
    for warning in result.get("gt_warnings", []):
        logger.warning(f"Gene {gene}: {warning}")
    if result.get("mac_filtered"):
        logger.debug(f"Gene {gene}: MAC < 5 — regression will report NA")


def _discard_per_gene_matrix_payload(gene_data: dict[str, Any]) -> None:
    for key in _PER_GENE_MATRIX_PAYLOAD_KEYS:
        gene_data.pop(key, None)
```

- [ ] **Step 4: Replace repeated builder-result assignment blocks**

In the first-gene parallel path, replace the manual assignments and warning loop with:

```python
_apply_builder_result_to_gene_data(first_gene_data, _result, first_gene)
```

Replace:

```python
first_gene_data.pop("genotype_matrix", None)
first_gene_data.pop("variant_mafs", None)
```

with:

```python
_discard_per_gene_matrix_payload(first_gene_data)
```

In the remaining parallel path, replace the manual assignments and warning loop with:

```python
_apply_builder_result_to_gene_data(gd, _result, _gene)
```

Replace the cleanup loop body with:

```python
_discard_per_gene_matrix_payload(gd)
```

In the sequential path, replace the manual assignments and warning loop with:

```python
_apply_builder_result_to_gene_data(gene_data, result, gene)
```

Replace:

```python
gene_data.pop("genotype_matrix", None)
gene_data.pop("variant_mafs", None)
```

with:

```python
_discard_per_gene_matrix_payload(gene_data)
```

- [ ] **Step 5: Run tests to verify they pass**

Run:

```bash
pytest tests/unit/test_streaming_matrix.py::test_engine_passes_builder_score_values_to_tests_and_discards_payload tests/unit/test_streaming_matrix.py::TestEngineBuilderConsumption::test_matrix_discarded_after_sequential_run -q
```

Expected: all selected tests pass.

- [ ] **Step 6: Commit engine handoff**

```bash
git add variantcentrifuge/association/engine.py tests/unit/test_streaming_matrix.py
git commit -m "feat: pass aligned annotation arrays from lazy builder"
```

---

### Task 7: Use Score Values In Burden Tests

**Files:**
- Modify: `variantcentrifuge/association/tests/logistic_burden.py`
- Modify: `variantcentrifuge/association/tests/linear_burden.py`
- Test: `tests/unit/test_association_logistic_burden.py`
- Test: `tests/unit/test_association_linear_burden.py`

- [ ] **Step 1: Write failing logistic burden test**

Add this test to `tests/unit/test_association_logistic_burden.py`:

```python
def test_logistic_burden_passes_score_values_to_weight_resolver(
    monkeypatch, synthetic_binary_data
) -> None:
    from variantcentrifuge.association import weights as weights_module

    geno, phenotype, mafs = synthetic_binary_data
    score_values = np.array([0.2, 0.4, 0.6, 0.8, 1.0], dtype=np.float64)
    seen_kwargs: dict[str, object] = {}

    def fake_get_weights(mafs_arg, weight_spec, **kwargs):
        seen_kwargs.update(kwargs)
        return np.ones(len(mafs_arg), dtype=np.float64)

    monkeypatch.setattr(weights_module, "get_weights", fake_get_weights)
    config = AssociationConfig(variant_weights="column:nephro_candidate_score")
    contingency = _make_contingency(geno, phenotype, mafs)
    contingency["score_values"] = score_values

    result = LogisticBurdenTest().run("GENE1", contingency, config)

    assert result.p_value is not None
    assert seen_kwargs["score_values"] is score_values
```

- [ ] **Step 2: Write failing linear burden raw score-only test**

Add this test to `tests/unit/test_association_linear_burden.py`:

```python
def test_linear_burden_raw_score_only_matches_manual_ols(
    synthetic_quantitative_data,
) -> None:
    geno, phenotype, mafs = synthetic_quantitative_data
    score_values = np.array([0.2, 0.4, 0.6, 0.8, 1.0], dtype=np.float64)
    config = AssociationConfig(
        variant_weights="column:nephro_candidate_score",
        variant_weight_params={"combine_with_beta": False},
    )

    burden = geno @ score_values
    design = sm.add_constant(burden.reshape(-1, 1))
    fit = sm.OLS(phenotype, design).fit()

    contingency = _make_contingency(geno, phenotype, mafs)
    contingency["score_values"] = score_values
    result = LinearBurdenTest().run("GENE1", contingency, config)

    assert result.p_value is not None
    assert result.p_value == pytest.approx(float(fit.pvalues[1]), rel=1e-6)
    assert result.effect_size == pytest.approx(float(fit.params[1]), rel=1e-6)
```

- [ ] **Step 3: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_association_logistic_burden.py::test_logistic_burden_passes_score_values_to_weight_resolver tests/unit/test_association_linear_burden.py::test_linear_burden_raw_score_only_matches_manual_ols -q
```

Expected:

- Logistic test fails because `score_values` is not passed into `get_weights()`.
- Linear test fails because `get_weights("column:...")` requires `score_values`, but `LinearBurdenTest` does not pass it yet.

- [ ] **Step 4: Pass score values in burden tests**

In both `variantcentrifuge/association/tests/logistic_burden.py` and `variantcentrifuge/association/tests/linear_burden.py`, update the `get_weights()` call:

```python
weights = get_weights(
    mafs,
    config.variant_weights,
    cadd_scores=contingency_data.get("cadd_scores"),
    revel_scores=contingency_data.get("revel_scores"),
    score_values=contingency_data.get("score_values"),
    variant_effects=contingency_data.get("variant_effects"),
    weight_params=config.variant_weight_params,
)
```

- [ ] **Step 5: Run tests to verify they pass**

Run:

```bash
pytest tests/unit/test_association_logistic_burden.py::test_logistic_burden_passes_score_values_to_weight_resolver tests/unit/test_association_linear_burden.py::test_linear_burden_raw_score_only_matches_manual_ols -q
```

Expected: both tests pass.

- [ ] **Step 6: Run burden regression suites**

Run:

```bash
pytest tests/unit/test_association_logistic_burden.py tests/unit/test_association_linear_burden.py -q
```

Expected: all logistic and linear burden tests pass.

- [ ] **Step 7: Commit burden integration**

```bash
git add variantcentrifuge/association/tests/logistic_burden.py variantcentrifuge/association/tests/linear_burden.py tests/unit/test_association_logistic_burden.py tests/unit/test_association_linear_burden.py
git commit -m "feat: use score-column weights in burden tests"
```

---

### Task 8: Make Python SKAT Backend Consume Concrete Weights

**Files:**
- Modify: `variantcentrifuge/association/backends/base.py`
- Modify: `variantcentrifuge/association/backends/python_backend.py`
- Test: `tests/unit/test_skat_python_backend.py`

- [ ] **Step 1: Write failing backend weight tests**

Add these tests to `tests/unit/test_skat_python_backend.py`:

```python
def test_backend_skat_uses_explicit_weight_vector(monkeypatch):
    from variantcentrifuge.association.backends.base import NullModelResult

    backend = PythonSKATBackend()
    geno = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    weights = np.array([2.0, 3.0, 5.0], dtype=np.float64)
    null_model = NullModelResult(
        model=None,
        trait_type="quantitative",
        n_samples=4,
        adjustment=False,
        extra={
            "residuals": np.array([0.5, -0.5, 0.25, -0.25]),
            "sigma2": 1.0,
            "design_matrix": np.ones((4, 1), dtype=np.float64),
        },
    )
    captured: dict[str, np.ndarray] = {}

    def fake_eigenvalues(geno_weighted, null_model_arg):
        captured["geno_weighted"] = geno_weighted.copy()
        return np.array([1.0], dtype=np.float64)

    monkeypatch.setattr(backend, "_compute_eigenvalues_filtered", fake_eigenvalues)
    monkeypatch.setattr(
        "variantcentrifuge.association.backends.python_backend.compute_pvalue",
        lambda q_stat, lambdas: (0.5, "liu", True),
    )

    result = backend._test_skat("GENE1", geno, null_model, weights)

    assert result["p_value"] == 0.5
    np.testing.assert_allclose(captured["geno_weighted"], geno * weights[np.newaxis, :])


def test_backend_burden_uses_explicit_weight_vector():
    from variantcentrifuge.association.backends.base import NullModelResult

    backend = PythonSKATBackend()
    geno = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    weights = np.array([1.0, 10.0, 100.0], dtype=np.float64)
    residuals = np.array([1.0, -1.0, 0.5, -0.5], dtype=np.float64)
    null_model = NullModelResult(
        model=None,
        trait_type="quantitative",
        n_samples=4,
        adjustment=False,
        extra={"residuals": residuals, "sigma2": 1.0},
    )

    result = backend._test_burden("GENE1", geno, null_model, weights)

    burden = geno @ weights
    score = float(residuals @ burden)
    variance = float(burden @ burden)
    expected_p = float(2.0 * scipy.stats.norm.sf(abs(score / np.sqrt(variance))))
    assert result["p_value"] == pytest.approx(expected_p)


def test_backend_rejects_weight_length_mismatch():
    from variantcentrifuge.association.backends.base import NullModelResult

    backend = PythonSKATBackend()
    geno = np.ones((4, 3), dtype=np.float64)
    null_model = NullModelResult(
        model=None,
        trait_type="quantitative",
        n_samples=4,
        adjustment=False,
        extra={"residuals": np.ones(4), "sigma2": 1.0},
    )

    with pytest.raises(ValueError, match="weights length"):
        backend.test_gene(
            gene="GENE1",
            genotype_matrix=geno,
            null_model=null_model,
            method="SKAT",
            weights=np.array([1.0, 2.0]),
        )
```

If `scipy` is not imported in `tests/unit/test_skat_python_backend.py`, add:

```python
import scipy.stats
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_skat_python_backend.py::test_backend_skat_uses_explicit_weight_vector tests/unit/test_skat_python_backend.py::test_backend_burden_uses_explicit_weight_vector tests/unit/test_skat_python_backend.py::test_backend_rejects_weight_length_mismatch -q
```

Expected: tests fail because backend methods still accept `weights_beta` and recompute Beta weights internally.

- [ ] **Step 3: Update backend abstract contract**

In `variantcentrifuge/association/backends/base.py`, change `test_gene()` signature to:

```python
def test_gene(
    self,
    gene: str,
    genotype_matrix: np.ndarray,
    null_model: NullModelResult,
    method: str,
    weights: np.ndarray | None = None,
    weights_beta: tuple[float, float] | None = None,
) -> dict[str, Any]:
```

Update the docstring:

```python
weights : np.ndarray or None
    Concrete per-variant weights aligned to genotype_matrix columns. Python SKAT
    requires this path.
weights_beta : tuple of (float, float) or None
    Legacy Beta parameters used by the deprecated R backend.
```

- [ ] **Step 4: Update Python backend to use explicit weights**

In `variantcentrifuge/association/backends/python_backend.py`, remove the `beta_maf_weights` import if no longer used in the backend.

Add this helper method inside `PythonSKATBackend`:

```python
def _validate_variant_weights(self, genotype_matrix: np.ndarray, weights: np.ndarray | None) -> np.ndarray:
    if weights is None:
        raise ValueError("PythonSKATBackend requires concrete variant weights")
    weights_arr = np.asarray(weights, dtype=np.float64)
    n_variants = genotype_matrix.shape[1]
    if len(weights_arr) != n_variants:
        raise ValueError(
            f"weights length must match genotype_matrix columns "
            f"(got {len(weights_arr)} and {n_variants})"
        )
    if not np.isfinite(weights_arr).all():
        raise ValueError("weights must be finite")
    return weights_arr
```

Change `PythonSKATBackend.test_gene()` signature to:

```python
def test_gene(
    self,
    gene: str,
    genotype_matrix: np.ndarray,
    null_model: NullModelResult,
    method: str,
    weights: np.ndarray | None = None,
    weights_beta: tuple[float, float] | None = None,
) -> dict[str, Any]:
```

At the start of `test_gene()` after `n_variants`, add:

```python
weights_arr = self._validate_variant_weights(genotype_matrix, weights)
```

Update method dispatch to pass `weights_arr`:

```python
if method_upper == "SKAT":
    result = self._test_skat(gene, genotype_matrix, null_model, weights_arr)
elif method_upper == "BURDEN":
    result = self._test_burden(gene, genotype_matrix, null_model, weights_arr)
elif method_upper in ("SKATO", "OPTIMAL.ADJ", "OPTIMAL"):
    result = self._test_skato(gene, genotype_matrix, null_model, weights_arr)
else:
    logger.warning(f"Unknown SKAT method '{method}' for gene {gene}; using SKAT")
    result = self._test_skat(gene, genotype_matrix, null_model, weights_arr)
```

Change `_test_skat()`, `_test_burden()`, and `_test_skato()` signatures from `weights_beta` to:

```python
weights: np.ndarray
```

Delete the internal `mafs = geno.mean(axis=0) / 2.0` and `beta_maf_weights(...)` lines in each method. Use the provided `weights` when forming `geno_weighted` or burden:

```python
geno_weighted = geno * weights[np.newaxis, :]
```

and:

```python
burden = geno @ weights
```

- [ ] **Step 5: Update existing backend tests that call `test_gene()`**

In `tests/unit/test_skat_python_backend.py`, replace calls that pass a Beta tuple
either by keyword or positionally:

```python
weights_beta=(1.0, 25.0)
```

or:

```python
backend.test_gene("GENE1", genotype_matrix, null_model, "SKAT", (1.0, 25.0))
```

with:

```python
weights=np.ones(genotype_matrix.shape[1], dtype=np.float64)
```

For tests that specifically compare Beta behavior, compute:

```python
mafs = genotype_matrix.mean(axis=0) / 2.0
weights = beta_maf_weights(mafs, a=1.0, b=25.0)
```

and pass `weights=weights`.

- [ ] **Step 6: Run backend tests to verify they pass**

Run:

```bash
pytest tests/unit/test_skat_python_backend.py -q
```

Expected: all Python backend tests pass.

- [ ] **Step 7: Commit Python backend explicit weights**

```bash
git add variantcentrifuge/association/backends/base.py variantcentrifuge/association/backends/python_backend.py tests/unit/test_skat_python_backend.py
git commit -m "feat: pass explicit weights to Python SKAT backend"
```

---

### Task 9: Resolve Shared Weights In Pure Python SKAT

**Files:**
- Modify: `variantcentrifuge/association/tests/skat_python.py`
- Create: `tests/unit/test_skat_python_weights.py`
- Modify: `tests/unit/test_acat_v.py`
- Modify: `tests/unit/test_skat_python_comparison.py`

- [ ] **Step 1: Write failing PurePythonSKAT shared resolver tests**

Create `tests/unit/test_skat_python_weights.py`:

```python
"""Tests for PurePythonSKATTest variant weight resolution."""

from __future__ import annotations

from types import MethodType

import numpy as np
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest
from variantcentrifuge.association.weights import beta_maf_weights, get_weights


def _contingency(geno, phenotype, **extra):
    data = {
        "genotype_matrix": geno,
        "phenotype_vector": phenotype,
        "covariate_matrix": None,
        "proband_count": int(phenotype.sum()),
        "control_count": int((phenotype == 0).sum()),
        "n_qualifying_variants": geno.shape[1],
    }
    data.update(extra)
    return data


def _make_test_with_capture(monkeypatch):
    test = PurePythonSKATTest()
    test.check_dependencies()
    captured: dict[str, np.ndarray] = {}

    def fake_test_gene(self, gene, genotype_matrix, null_model, method, weights=None, weights_beta=None):
        captured["backend_weights"] = np.asarray(weights, dtype=np.float64).copy()
        return {
            "p_value": 0.5,
            "rho": None,
            "n_variants": genotype_matrix.shape[1],
            "n_marker_test": genotype_matrix.shape[1],
            "warnings": [],
            "p_method": "test",
            "p_converged": True,
            "skip_reason": None,
        }

    test._backend.test_gene = MethodType(fake_test_gene, test._backend)

    def fake_acat_v(**kwargs):
        captured["acat_v_weights"] = np.asarray(kwargs["weights"], dtype=np.float64).copy()
        return 0.75

    monkeypatch.setattr(
        "variantcentrifuge.association.tests.skat_python.compute_acat_v",
        fake_acat_v,
    )
    return test, captured


def test_python_skat_resolves_score_column_weights_for_backend_and_acat_v(monkeypatch):
    geno = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    phenotype = np.array([1, 1, 1, 0, 0, 0], dtype=float)
    score_values = np.array([0.2, 0.5, 1.0], dtype=np.float64)
    config = AssociationConfig(
        trait_type="binary",
        variant_weights="column:nephro_candidate_score",
        variant_weight_params={"combine_with_beta": False},
    )
    test, captured = _make_test_with_capture(monkeypatch)

    result = test.run(
        "GENE1",
        _contingency(geno, phenotype, score_values=score_values),
        config,
    )

    expected = get_weights(
        geno.mean(axis=0) / 2.0,
        "column:nephro_candidate_score",
        score_values=score_values,
        weight_params={"combine_with_beta": False},
    )
    assert result.extra["acat_v_p"] == 0.75
    np.testing.assert_allclose(captured["backend_weights"], expected)
    np.testing.assert_allclose(captured["acat_v_weights"], expected)


@pytest.mark.parametrize(
    "variant_weights, extra_key, extra_values",
    [
        ("cadd", "cadd_scores", np.array([30.0, 20.0, 10.0], dtype=np.float64)),
        ("revel", "revel_scores", np.array([0.9, 0.5, 0.1], dtype=np.float64)),
        ("combined", "cadd_scores", np.array([30.0, 20.0, 10.0], dtype=np.float64)),
    ],
)
def test_python_skat_honors_functional_weight_schemes(
    monkeypatch,
    variant_weights,
    extra_key,
    extra_values,
):
    geno = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    phenotype = np.array([1, 1, 1, 0, 0, 0], dtype=float)
    config = AssociationConfig(trait_type="binary", variant_weights=variant_weights)
    test, captured = _make_test_with_capture(monkeypatch)
    contingency = _contingency(geno, phenotype, **{extra_key: extra_values})

    test.run("GENE1", contingency, config)

    expected = get_weights(
        geno.mean(axis=0) / 2.0,
        variant_weights,
        cadd_scores=contingency.get("cadd_scores"),
        revel_scores=contingency.get("revel_scores"),
    )
    np.testing.assert_allclose(captured["backend_weights"], expected)


def test_python_skat_beta_uses_post_imputation_genotype_maf(monkeypatch):
    geno = np.array(
        [
            [0.0, 1.0, 2.0],
            [0.0, 1.0, 2.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    phenotype = np.array([1, 1, 1, 0, 0, 0], dtype=float)
    pre_imputation_mafs = np.array([0.01, 0.01, 0.01], dtype=np.float64)
    config = AssociationConfig(trait_type="binary", variant_weights="beta:1,25")
    test, captured = _make_test_with_capture(monkeypatch)

    test.run(
        "GENE1",
        _contingency(geno, phenotype, variant_mafs=pre_imputation_mafs),
        config,
    )

    expected = beta_maf_weights(geno.mean(axis=0) / 2.0, a=1.0, b=25.0)
    not_expected = beta_maf_weights(pre_imputation_mafs, a=1.0, b=25.0)
    np.testing.assert_allclose(captured["backend_weights"], expected)
    assert not np.allclose(captured["backend_weights"], not_expected)
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_skat_python_weights.py -q
```

Expected:

- Score-column test fails because `PurePythonSKATTest.run()` still calls `parse_weights_beta()`.
- Functional scheme tests fail because `cadd`, `revel`, and `combined` are still reduced to Beta-only in Python SKAT.
- Post-imputation MAF test may pass for Beta only before the implementation, and it must still pass after the implementation.

- [ ] **Step 3: Resolve weights in `PurePythonSKATTest.run()`**

In `variantcentrifuge/association/tests/skat_python.py`, replace imports:

```python
from variantcentrifuge.association.tests._utils import parse_weights_beta
from variantcentrifuge.association.weights import beta_maf_weights
```

with:

```python
from variantcentrifuge.association.weights import get_weights
```

Replace the `weights_beta = parse_weights_beta(...)` block and backend call with:

```python
skat_mafs = geno.mean(axis=0) / 2.0
weights = get_weights(
    skat_mafs,
    config.variant_weights,
    cadd_scores=contingency_data.get("cadd_scores"),
    revel_scores=contingency_data.get("revel_scores"),
    score_values=contingency_data.get("score_values"),
    variant_effects=contingency_data.get("variant_effects"),
    weight_params=config.variant_weight_params,
)

result = self._backend.test_gene(
    gene=gene,
    genotype_matrix=geno,
    null_model=self._null_model,
    method=config.skat_method,
    weights=weights,
)
```

Replace ACAT-V weight computation:

```python
mafs = geno.mean(axis=0) / 2.0
a1, a2 = weights_beta
acat_v_weights = beta_maf_weights(mafs, a=a1, b=a2)
```

with:

```python
acat_v_weights = weights
```

- [ ] **Step 4: Update Pure Python SKAT comparison tests for backend signature**

In `tests/unit/test_skat_python_comparison.py`, update direct `PythonSKATBackend.test_gene()` calls the same way as Task 8:

```python
weights=np.ones(genotype_matrix.shape[1], dtype=np.float64)
```

or Beta vectors computed with:

```python
weights = beta_maf_weights(genotype_matrix.mean(axis=0) / 2.0, a=1.0, b=25.0)
```

- [ ] **Step 5: Update ACAT-V tests that indirectly depend on `PurePythonSKATTest.run()`**

In `tests/unit/test_acat_v.py`, update expected comments so they state ACAT-V uses the same resolved variant weights as SKAT. No assertion should expect Beta-only weights after this change.

If any test monkeypatches `PythonSKATBackend.test_gene()` with a `weights_beta`-only fake, update it to accept:

```python
def fake_test_gene(self, gene, genotype_matrix, null_model, method, weights=None, weights_beta=None):
    ...
```

- [ ] **Step 6: Run SKAT weight tests and comparison tests**

Run:

```bash
pytest tests/unit/test_skat_python_weights.py tests/unit/test_acat_v.py tests/unit/test_skat_python_comparison.py -q
```

Expected: all selected tests pass.

- [ ] **Step 7: Commit Pure Python SKAT shared resolver use**

```bash
git add variantcentrifuge/association/tests/skat_python.py tests/unit/test_skat_python_weights.py tests/unit/test_acat_v.py tests/unit/test_skat_python_comparison.py
git commit -m "feat: resolve shared weights for Python SKAT"
```

---

### Task 10: Guard Deprecated R SKAT And Preserve Fisher/COAST Semantics

**Files:**
- Modify: `variantcentrifuge/association/engine.py`
- Modify: `variantcentrifuge/association/tests/skat_r.py`
- Test: `tests/unit/test_skat_r_test.py`
- Test: `tests/unit/test_association_engine.py`
- Test: `tests/unit/test_association_fisher.py`
- Test: `tests/unit/test_coast_python.py`

- [ ] **Step 1: Write failing R SKAT guard tests**

Add this test to `tests/unit/test_association_engine.py`:

```python
def test_engine_rejects_r_skat_with_score_column_before_dependency_check():
    from variantcentrifuge.association.base import AssociationConfig
    from variantcentrifuge.association.engine import AssociationEngine

    config = AssociationConfig(
        skat_backend="r",
        variant_weights="column:nephro_candidate_score",
    )

    with pytest.raises(ValueError, match="R SKAT backend supports only beta:a,b and uniform"):
        AssociationEngine.from_names(["skat"], config)
```

Add this test to `tests/unit/test_skat_r_test.py`:

```python
def test_r_skat_run_rejects_functional_weights_without_r_backend():
    from variantcentrifuge.association.base import AssociationConfig
    from variantcentrifuge.association.tests.skat_r import RSKATTest

    config = AssociationConfig(skat_backend="r", variant_weights="cadd")
    contingency_data = {
        "genotype_matrix": np.ones((4, 2), dtype=np.float64),
        "phenotype_vector": np.array([1.0, 1.0, 0.0, 0.0]),
        "covariate_matrix": None,
        "proband_count": 2,
        "control_count": 2,
        "n_qualifying_variants": 2,
    }

    test = RSKATTest()

    with pytest.raises(ValueError, match="Use --skat-backend python"):
        test.run("GENE1", contingency_data, config)
```

- [ ] **Step 2: Write Fisher and COAST semantic regression tests**

Add this test to `tests/unit/test_association_fisher.py`:

```python
def test_fisher_ignores_score_column_variant_weights():
    from variantcentrifuge.association.base import AssociationConfig
    from variantcentrifuge.association.tests.fisher import FisherExactTest

    contingency_data = {
        "proband_carrier_count": 4,
        "control_carrier_count": 1,
        "proband_count": 10,
        "control_count": 10,
        "n_qualifying_variants": 3,
        "score_values": np.array([0.1, 0.5, 1.0]),
    }
    config = AssociationConfig(
        variant_weights="column:nephro_candidate_score",
        variant_weight_column="nephro_candidate_score",
    )

    result = FisherExactTest().run("GENE1", contingency_data, config)

    assert result.test_name == "fisher"
    assert result.p_value is not None
    assert result.n_variants == 3
```

Add this test to `tests/unit/test_coast_python.py` in an existing `PurePythonCOASTTest` integration section:

```python
def test_python_coast_does_not_require_score_values_for_score_column_variant_weights():
    from variantcentrifuge.association.base import AssociationConfig
    from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest

    config = AssociationConfig(
        trait_type="binary",
        variant_weights="column:nephro_candidate_score",
        variant_weight_column="nephro_candidate_score",
        coast_backend="python",
    )
    test = PurePythonCOASTTest()
    contingency_data = {
        "genotype_matrix": np.zeros((4, 0), dtype=np.float64),
        "phenotype_vector": np.array([1.0, 1.0, 0.0, 0.0]),
        "covariate_matrix": None,
        "gene_df": pd.DataFrame({"GENE": pd.Series([], dtype=str)}),
        "proband_count": 2,
        "control_count": 2,
        "n_qualifying_variants": 0,
    }

    result = test.run("GENE1", contingency_data, config)

    assert result.test_name == "coast"
    assert result.p_value is None
```

- [ ] **Step 3: Run tests to verify they fail**

Run:

```bash
pytest tests/unit/test_association_engine.py::test_engine_rejects_r_skat_with_score_column_before_dependency_check tests/unit/test_skat_r_test.py::test_r_skat_run_rejects_functional_weights_without_r_backend tests/unit/test_association_fisher.py::test_fisher_ignores_score_column_variant_weights tests/unit/test_coast_python.py::test_python_coast_does_not_require_score_values_for_score_column_variant_weights -q
```

Expected:

- R SKAT tests fail because unsupported specs are not rejected before Beta parsing.
- Fisher should pass before implementation if the test data matches existing Fisher keys, and it must still pass after implementation.
- COAST should pass before implementation if the chosen fixture reaches an early return, and it must still pass after implementation.

- [ ] **Step 4: Implement R SKAT support guard**

In `variantcentrifuge/association/tests/skat_r.py`, add this helper near `_parse_weights_beta()`:

```python
def _is_r_skat_weight_supported(variant_weights: str) -> bool:
    return variant_weights == "uniform" or variant_weights.startswith("beta:")


def _raise_unsupported_r_skat_weight(variant_weights: str) -> None:
    raise ValueError(
        "R SKAT backend supports only beta:a,b and uniform variant weights. "
        f"Got '{variant_weights}'. Use --skat-backend python for cadd, revel, "
        "combined, score_column, or column:<name> weights."
    )
```

At the start of `RSKATTest.run()` after early matrix/phenotype checks and before backend initialization checks, add:

```python
if not _is_r_skat_weight_supported(config.variant_weights):
    _raise_unsupported_r_skat_weight(config.variant_weights)
```

In `variantcentrifuge/association/engine.py`, before the dependency-check loop in `AssociationEngine.from_names()`, add:

```python
if skat_backend == "r" and "skat" in test_names:
    from variantcentrifuge.association.tests.skat_r import _is_r_skat_weight_supported

    if not _is_r_skat_weight_supported(config.variant_weights):
        raise ValueError(
            "R SKAT backend supports only beta:a,b and uniform variant weights. "
            f"Got '{config.variant_weights}'. Use --skat-backend python for cadd, revel, "
            "combined, score_column, or column:<name> weights."
        )
```

- [ ] **Step 5: Run tests to verify they pass**

Run:

```bash
pytest tests/unit/test_association_engine.py::test_engine_rejects_r_skat_with_score_column_before_dependency_check tests/unit/test_skat_r_test.py::test_r_skat_run_rejects_functional_weights_without_r_backend tests/unit/test_association_fisher.py::test_fisher_ignores_score_column_variant_weights tests/unit/test_coast_python.py::test_python_coast_does_not_require_score_values_for_score_column_variant_weights -q
```

Expected: all selected tests pass.

- [ ] **Step 6: Commit R guard and semantic regressions**

```bash
git add variantcentrifuge/association/engine.py variantcentrifuge/association/tests/skat_r.py tests/unit/test_skat_r_test.py tests/unit/test_association_engine.py tests/unit/test_association_fisher.py tests/unit/test_coast_python.py
git commit -m "feat: guard unsupported R SKAT functional weights"
```

---

### Task 11: Add End-To-End Score-Column Association Regression

**Files:**
- Create: `tests/unit/test_score_column_association_integration.py`

- [ ] **Step 1: Write failing integration regression**

Create `tests/unit/test_score_column_association_integration.py`:

```python
"""Small association regression for score-column weights across burden and Python SKAT."""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage


def test_score_column_weights_run_logistic_burden_and_python_skat(tmp_path):
    cases = [f"CASE{i}" for i in range(10)]
    controls = [f"CTRL{i}" for i in range(10)]
    samples = cases + controls

    rows = []
    for variant_idx, score in enumerate([2.0, 5.0, 9.0]):
        row = {
            "GENE": "GENE1",
            "nephro_candidate_score": score,
            "ANN_0__EFFECT": "missense_variant",
        }
        for sample_idx, _sample in enumerate(samples):
            if variant_idx == 0:
                gt = "0/1" if sample_idx < 8 else "0/0"
            elif variant_idx == 1:
                gt = "0/1" if 4 <= sample_idx < 14 else "0/0"
            else:
                gt = "1/1" if sample_idx in (0, 1, 10) else "0/0"
            row[f"GEN_{sample_idx}__GT"] = gt
        rows.append(row)

    df = pd.DataFrame(rows)
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.intermediate_dir.mkdir()
    workspace.get_intermediate_path = lambda name: workspace.intermediate_dir / name
    workspace.get_output_path = lambda name, ext=".tsv": workspace.output_dir / f"{name}{ext}"
    workspace.base_name = "score_column_integration"

    context = PipelineContext(
        args=Namespace(),
        config={
            "perform_association": True,
            "case_samples": cases,
            "control_samples": controls,
            "association_tests": ["logistic_burden", "skat"],
            "skat_backend": "python",
            "variant_weights": "score_column",
            "variant_weight_column": "nephro_candidate_score",
            "variant_weight_params": {
                "score_min": 0,
                "score_max": 10,
                "floor": 0.1,
                "combine_with_beta": True,
            },
            "output_dir": str(tmp_path),
            "output_file_base": "score_column_integration",
            "gzip_intermediates": False,
        },
        workspace=workspace,
    )
    context.current_dataframe = df
    context.vcf_samples = samples

    result_context = AssociationAnalysisStage()._process(context)

    result_df = result_context.association_results
    assert result_df is not None
    assert not result_df.empty
    assert "logistic_burden_pvalue" in result_df.columns
    assert "skat_pvalue" in result_df.columns
    assert "acat_o_pvalue" in result_df.columns
```

- [ ] **Step 2: Run test to verify it fails before all integration work is complete**

Run:

```bash
pytest tests/unit/test_score_column_association_integration.py -q
```

Expected: this fails until config, stage, builder, burden, and Python SKAT tasks are complete.

- [ ] **Step 3: Run test to verify it passes after prior tasks**

Run:

```bash
pytest tests/unit/test_score_column_association_integration.py -q
```

Expected: the integration regression passes and produces association results with burden, SKAT, and ACAT-O columns.

- [ ] **Step 4: Commit integration regression**

```bash
git add tests/unit/test_score_column_association_integration.py
git commit -m "test: cover score-column association workflow"
```

---

### Task 12: Update Association Documentation

**Files:**
- Modify: `docs/source/guides/association_testing.md`
- Modify: `docs/source/configuration.md`

- [ ] **Step 1: Add score-column docs**

In `docs/source/guides/association_testing.md`, update the variant weighting section around the existing `--variant-weights` documentation to include:

```markdown
### Score-Column Variant Weights

Variant-level weights affect each variant's contribution to burden scores and
SKAT kernels. They are separate from gene-level FDR prior weights, which affect
multiple-testing correction only.

Use `column:<name>` to reference any numeric per-variant column already present
in the association input DataFrame:

```bash
--perform-association \
--association-tests logistic_burden,skat \
--variant-weights column:nephro_candidate_score \
--variant-weight-params '{"score_min":0,"score_max":10,"floor":0.1,"combine_with_beta":true}'
```

The equivalent friendly form is:

```bash
--perform-association \
--association-tests logistic_burden,skat \
--variant-weights score_column \
--variant-weight-column nephro_candidate_score \
--variant-weight-params '{"score_min":0,"score_max":10,"floor":0.1,"combine_with_beta":true}'
```

Score-column parameters:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `score_min` | `null` | Lower bound for explicit linear normalization. |
| `score_max` | `null` | Upper bound for explicit linear normalization. |
| `floor` | `0.0` | Lower bound after clipping finite normalized scores. |
| `ceiling` | `1.0` | Upper bound after clipping finite normalized scores. |
| `combine_with_beta` | `true` | Multiply the functional score by `Beta(MAF; beta_a,beta_b)`. |
| `missing` | `null` | `neutral` in Beta-combined mode, `floor` in raw score-only mode. |
| `beta_a` | `1.0` | Beta alpha used when `combine_with_beta=true`. |
| `beta_b` | `25.0` | Beta beta used when `combine_with_beta=true`. |

When `combine_with_beta=false`, the normalized score itself is the final variant
weight. In that raw score-only mode, `missing="neutral"` is invalid because
there is no multiplicative Beta baseline.

Python SKAT, SKAT-O, and ACAT-V use the same resolved weight vector as burden
tests. For `beta:*`, Python SKAT preserves its historical MAF source:
`geno.mean(axis=0) / 2.0` after genotype imputation. Python SKAT now honors
`cadd`, `revel`, and `combined`; earlier versions silently reduced those
schemes to Beta-only weights in SKAT.

Fisher exact tests ignore variant weights because they use collapsed contingency
tables. COAST uses its own ordered allelic-series category weights. The
deprecated R SKAT backend supports only `beta:a,b` and `uniform`; use the
default Python backend for score-column or functional weights.
```
```

Update the existing config example to include:

```json
{
  "association": {
    "variant_weights": "score_column",
    "variant_weight_column": "nephro_candidate_score",
    "variant_weight_params": {
      "score_min": 0,
      "score_max": 10,
      "floor": 0.1,
      "combine_with_beta": true
    }
  }
}
```

Update the association config table row for `variant_weights` to list:

```markdown
`beta:a,b`, `uniform`, `cadd`, `revel`, `combined`, `column:<name>`, `score_column`
```

Add a row:

```markdown
| `variant_weight_column` | str or null | `null` | Score column used when `variant_weights` is `score_column` |
```

- [ ] **Step 2: Add top-level configuration cross-reference**

Add this short cross-reference to `docs/source/configuration.md`:

```markdown
### Association Score-Column Weights

For score-column variant weights, set `association.variant_weights` to
`score_column` and set `association.variant_weight_column` to the numeric
per-variant column name. See the association testing guide for normalization
parameters.
```

- [ ] **Step 3: Verify docs contain the new user-facing forms**

Run:

```bash
rg -n 'column:<name>|score_column|variant_weight_column|nephro_candidate_score|missing="neutral"' docs/source/guides/association_testing.md docs/source/configuration.md
```

Expected: output includes the score-column forms, the NCS example, `variant_weight_column`, and the raw-mode missing-value warning.

- [ ] **Step 4: Commit docs**

```bash
git add docs/source/guides/association_testing.md docs/source/configuration.md
git commit -m "docs: document score-column association weights"
```

---

### Task 13: Final Verification And Feature Commit Review

**Files:**
- Verify all modified implementation, test, and docs files from prior tasks.

- [ ] **Step 1: Run focused weight and config tests**

Run:

```bash
pytest tests/unit/test_score_column_weights.py tests/unit/test_functional_weights.py tests/unit/test_json_config.py tests/unit/test_cli_association_args.py tests/test_cli_argument_parsing.py -q
```

Expected: all selected tests pass.

- [ ] **Step 2: Run alignment and stage tests**

Run:

```bash
pytest tests/unit/test_association_genotype_matrix.py tests/unit/test_streaming_matrix.py tests/unit/test_score_column_association_stage.py tests/unit/test_score_column_association_integration.py -q
```

Expected: all selected tests pass.

- [ ] **Step 3: Run burden and SKAT tests**

Run:

```bash
pytest tests/unit/test_association_logistic_burden.py tests/unit/test_association_linear_burden.py tests/unit/test_skat_python_backend.py tests/unit/test_skat_python_weights.py tests/unit/test_acat_v.py tests/unit/test_skat_python_comparison.py -q
```

Expected: all selected tests pass.

- [ ] **Step 4: Run R guard and unaffected-test regressions**

Run:

```bash
pytest tests/unit/test_association_engine.py tests/unit/test_skat_r_test.py tests/unit/test_association_fisher.py tests/unit/test_coast_python.py -q
```

Expected: all selected tests pass. R-dependent tests that are already skipped in the local environment should continue to skip for the same dependency reasons; the new R guard tests must not require rpy2.

- [ ] **Step 5: Run compile check**

Run:

```bash
python -m compileall variantcentrifuge
```

Expected: command exits with status 0 and no syntax errors.

- [ ] **Step 6: Run full unit suite**

Run:

```bash
pytest tests/unit -q
```

Expected: all unit tests pass or retain pre-existing environment skips.

- [ ] **Step 7: Inspect diff for accidental unrelated edits**

Run:

```bash
git status --short
git diff --stat
git diff --check
```

Expected:

- `git status --short` shows only files listed in this plan.
- `git diff --check` reports no whitespace errors.

- [ ] **Step 8: Commit verification fixes only when verification changed files**

Run `git status --short`. When verification fixes changed files, commit them:

```bash
git add variantcentrifuge tests docs
git commit -m "test: verify score-column association weights"
```

Expected: the branch contains focused commits for config, resolver, alignment, burden, SKAT, docs, and verification fixes. When `git status --short` is empty after Step 7, skip this commit.

---

## Self-Review Against The Approved Spec

- Shared resolver: Tasks 2, 7, and 9 route burden tests and Python SKAT through `get_weights()` plus `score_column_weights()`.
- Existing `beta:*` and `uniform` behavior: Tasks 2, 8, 9, and 13 preserve existing unit tests and explicitly test Python SKAT Beta MAF source from `geno.mean(axis=0) / 2.0`.
- Python SKAT honors `cadd`, `revel`, and `combined`: Task 9 adds parametrized tests proving those specs pass concrete functional weights into the backend.
- `column:<name>` and `score_column + variant_weight_column`: Tasks 1, 2, and 5 cover CLI/config parsing, resolution, and stage validation.
- `variant_weight_column` in CLI, JSON validation, and `AssociationConfig`: Task 1 covers all three.
- Exact annotation/genotype alignment: Tasks 3, 4, 5, and 6 propagate and consume the real `keep_variants_mask`.
- Alignment inside lazy builder: Task 4 moves score/CADD/REVEL/effect extraction to `_GenotypeMatrixBuilder`; Task 5 deletes the stage-level duplicated regex mask.
- Raw score-only missing values: Task 2 rejects `missing="neutral"` when `combine_with_beta=false` and tests floor fallback.
- `beta_a` and `beta_b`: Task 2 implements and tests custom Beta shape parameters for score-column weights.
- R SKAT behavior: Task 10 fails clearly for non-Beta functional and score-column weights before R dependency checks.
- Fisher and COAST semantics: Task 10 adds regressions showing Fisher ignores variant weights and COAST does not require score values.
- Documentation: Task 12 documents score-column usage, parameters, defaults, Python SKAT behavior change, Fisher/COAST notes, and R SKAT limitation.
- Placeholder scan: This plan contains no unresolved placeholders and every code-changing task includes concrete test code, implementation snippets, commands, expected outcomes, and commit steps.
