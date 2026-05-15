# Score Column Variant Weights Design

## Goal

Implement GitHub issue #97: allow users to use an arbitrary numeric variant-level
score column as the per-variant weight source for association tests.

The primary user story is to run burden and SKAT-style association tests with
VariantCentrifuge-computed scores such as `nephro_candidate_score`:

```bash
--perform-association \
--association-tests logistic_burden,skat \
--variant-weights column:nephro_candidate_score
```

The feature must support both Beta(MAF)-combined score weights and raw score-only
weights without changing existing behavior for `beta:*`, `uniform`, `cadd`,
`revel`, or `combined`.

## Background

VariantCentrifuge already has two distinct prior-weight concepts:

- Variant-level association weights: affect each variant's contribution to the
  burden score or SKAT kernel.
- Gene-level FDR prior weights: affect multiple-testing correction only.

Issue #97 is about the first concept. Scores such as `nephro_candidate_score`
are computed per variant before association analysis and are present in the final
variant DataFrame. They should be usable directly as prior information during
association testing instead of only as hard filters or post-hoc report columns.

## Current State

Existing weight schemes are implemented in
`variantcentrifuge/association/weights.py`:

- `beta:a,b`
- `uniform`
- `cadd`
- `revel`
- `combined`

`logistic_burden` and `linear_burden` call `get_weights(...)` directly and can
receive annotation arrays through `contingency_data`.

Python SKAT follows a different path. It parses `variant_weights` into
`weights_beta`, passes only the two Beta parameters to the backend, and the
backend recomputes Beta(MAF) weights internally. As a result, arbitrary score
columns cannot be added safely by changing only `get_weights(...)`; SKAT would
ignore them or fall back to Beta(1,25).

The association stage already has a partial precedent for aligned per-variant
annotation arrays. For `cadd`, `revel`, and `combined`, it extracts annotation
columns from each per-gene DataFrame and replicates the site missingness filter
used by `build_genotype_matrix()`. That same alignment concern is central to
score-column weights.

## Scope

In scope:

- Add a generic score-column variant weight scheme.
- Support both `column:<name>` and a friendly `score_column` plus
  `--variant-weight-column <name>` form.
- Add `variant_weight_column` to CLI, JSON association config, and
  `AssociationConfig`.
- Extend weight resolution so score columns work in:
  - `logistic_burden`
  - `linear_burden`
  - Python `skat` / `skat_python`, including SKAT-O and ACAT-V
- Preserve all existing weight schemes and defaults.
- Validate requested columns and weight parameters with helpful errors.
- Add focused unit, stage, and regression tests.
- Update association documentation.

Out of scope:

- Changing Fisher exact behavior. Fisher uses collapsed contingency tables and
  does not consume per-variant weights.
- Changing COAST category weighting. COAST has its own ordered allelic-series
  weights; score-column support for COAST can be designed later.
- Supporting arbitrary Python expressions in `--variant-weights`. Users should
  create score columns through the existing scoring system, then reference the
  resulting column.
- Changing gene-level FDR weighting behavior.
- Supporting custom score-column weights in the deprecated R SKAT backend unless
  that can be done with a small explicit `weights=` vector path. The Python SKAT
  backend is the supported default and should be the implementation target.

## Design Options

### Option A: Add Score Handling Only To Burden Tests

This would extend `get_weights(...)` and the burden tests, leaving SKAT unchanged.

Pros:

- Smallest initial patch.
- Matches the existing burden-test flow.

Cons:

- Fails one of the issue's explicit requirements.
- Creates a dangerous behavior gap where `--association-tests logistic_burden,skat`
  uses score weights for one test and Beta weights for the other.

### Option B: Resolve A Shared Weight Vector Per Gene

Create one weight-resolution path that produces a concrete per-variant weight
array aligned to the genotype matrix. Burden tests use it for `geno @ weights`;
Python SKAT uses it for its weighted genotype matrix and ACAT-V weights.

Pros:

- One statistical interpretation across burden, SKAT, SKAT-O, and ACAT-V.
- Keeps weight-spec parsing centralized.
- Avoids silent SKAT divergence.
- Fits the existing `get_weights(...)` abstraction with a moderate extension.

Cons:

- Requires backend signature changes for Python SKAT.
- Needs careful regression tests for alignment and backward compatibility.

### Option C: Encode Score Weights As Synthetic Beta Parameters

Try to express score-column weights through `weights_beta` or special sentinel
values passed to existing SKAT backend code.

Pros:

- Minimal signature changes.

Cons:

- Conceptually wrong. Score columns are not Beta distribution parameters.
- Increases hidden branching and fallback risk.
- Harder to test and document clearly.

## Selected Design

Use Option B: resolve a shared per-gene weight vector and pass that vector to all
association tests that consume variant weights.

`column:<name>` is the canonical user-facing spec. `score_column` plus
`variant_weight_column` is an alias for readability and JSON configuration.

Recommended example for `nephro_candidate_score`:

```bash
--variant-weights column:nephro_candidate_score \
--variant-weight-params '{"score_min":0,"score_max":10,"floor":0.1,"combine_with_beta":true}'
```

Equivalent friendly form:

```bash
--variant-weights score_column \
--variant-weight-column nephro_candidate_score \
--variant-weight-params '{"score_min":0,"score_max":10,"floor":0.1,"combine_with_beta":true}'
```

## CLI And Config Contract

Add:

```bash
--variant-weight-column <column-name>
```

`AssociationConfig` gains:

```python
variant_weight_column: str | None = None
```

The JSON `association` section accepts:

```json
{
  "variant_weights": "score_column",
  "variant_weight_column": "nephro_candidate_score",
  "variant_weight_params": {
    "score_min": 0,
    "score_max": 10,
    "floor": 0.1,
    "combine_with_beta": true
  }
}
```

Supported specs after the change:

- `beta:a,b`
- `uniform`
- `cadd`
- `revel`
- `combined`
- `column:<column_name>`
- `score_column` with `variant_weight_column`

Validation rules:

- `score_column` requires `variant_weight_column`.
- `column:<name>` must have a non-empty `<name>`.
- If both `column:<name>` and `variant_weight_column` are set, the inline
  `column:<name>` value wins. A debug log may note the ignored separate column.
- Unknown specs continue to raise `ValueError`.
- `variant_weight_column` must be a string in JSON config when present.
- `variant_weight_params` must be a dict after CLI JSON parsing.

## Score Weight Semantics

Score-column weights use these params:

```json
{
  "score_min": null,
  "score_max": null,
  "floor": 0.0,
  "ceiling": 1.0,
  "combine_with_beta": true,
  "missing": "neutral"
}
```

Meaning:

- `score_min` / `score_max`: explicit range for linear normalization. When both
  are provided, compute `(score - score_min) / (score_max - score_min)`.
- If no range is provided, treat scores as already on a 0-1 scale.
- `floor`: lower bound for finite normalized score weights after clipping.
- `ceiling`: upper bound for finite normalized score weights after clipping.
- `combine_with_beta`: if true, final weight is `Beta(MAF; a,b) * functional`.
  If false, final weight is the normalized score weight only.
- `missing`: initially supports `neutral` and `floor`.
  - `neutral`: missing or invalid scores get functional weight 1.0 when
    `combine_with_beta=true`, or 1.0 when raw score-only mode is used.
  - `floor`: missing or invalid scores get the configured floor.

Default behavior:

- `combine_with_beta=true`
- no implicit per-gene min-max normalization
- missing scores are neutral

No implicit per-gene normalization should be added. Per-gene scaling would make
weights depend on which variants survive filtering in each gene and would reduce
comparability across genes.

For NCS-like scores, users should pass explicit range parameters:

```json
{"score_min": 0, "score_max": 10, "floor": 0.1, "combine_with_beta": true}
```

## Weight Resolver API

Extend `get_weights(...)` to accept:

```python
score_values: np.ndarray | None = None
weight_column: str | None = None
```

Add a public helper:

```python
score_column_weights(
    mafs: np.ndarray,
    score_values: np.ndarray,
    *,
    variant_effects: np.ndarray | None = None,
    weight_params: dict | None = None,
) -> np.ndarray
```

`score_column_weights()` is responsible for:

- numeric parsing using the same missing-value conventions as CADD/REVEL
- explicit score range normalization
- clipping to floor/ceiling
- missing/invalid fallback
- optional Beta(MAF) multiplication
- logging missing/invalid counts
- returning a finite `float64` array with the same length as `mafs`

Length mismatch between `mafs` and `score_values` is a hard `ValueError` because
it means annotation arrays are not aligned to the genotype matrix.

## Association Data Flow

The association stage should detect score-column weighting before building
per-gene `gene_data`:

1. Resolve the requested score column name from `variant_weights` or
   `variant_weight_column`.
2. Validate the column exists in the association input DataFrame before the gene
   loop starts.
3. For each gene, extract the score column from `gene_df`.
4. Apply the same site missingness keep mask used for CADD/REVEL alignment.
5. Store the aligned array as `gene_data["score_values"]`.
6. Store the chosen column name as `gene_data["variant_weight_column"]` for
   diagnostics and test error messages.

The existing CADD/REVEL extraction block should be generalized rather than
duplicated again. A small helper can compute the site keep mask once per gene
and apply it to any requested annotation arrays.

## Alignment Strategy

Correct alignment means:

- `genotype_matrix.shape[1]`
- `variant_mafs.shape[0]`
- `score_values.shape[0]`
- `cadd_scores.shape[0]` when present
- `revel_scores.shape[0]` when present
- `variant_effects.shape[0]` when present

must all describe the same post-site-filter variants in the same order.

The first implementation can preserve the current missingness-mask replication,
but it should isolate it in a helper to reduce future drift. A better follow-up
would be for `build_genotype_matrix()` to return `keep_variants_mask`, letting
annotation alignment reuse the exact mask instead of recomputing it.

The implementation plan should consider returning `keep_variants_mask` from
`build_genotype_matrix()` if that can be done compatibly. If not, the helper must
mirror the existing vectorized missingness logic and tests must guard it.

## Burden Tests

`logistic_burden` and `linear_burden` continue to call `get_weights(...)`.

They should pass:

- `score_values=contingency_data.get("score_values")`
- `variant_effects=contingency_data.get("variant_effects")`
- `weight_params=config.variant_weight_params`

No additional burden-specific branch is needed.

## Python SKAT Integration

Python SKAT must receive the same resolved variant weights as burden tests.

Change the Python backend contract from Beta-only parameters to either:

```python
variant_weights: np.ndarray | None = None
weights_beta: tuple[float, float] | None = None
```

or, preferably after resolving at the test layer:

```python
weights: np.ndarray
```

`PurePythonSKATTest.run()` should resolve weights with `get_weights(...)` using
`variant_mafs` from `contingency_data`, then pass the concrete vector to
`PythonSKATBackend.test_gene()`.

The backend should use that vector in:

- `_test_skat()`
- `_test_burden()`
- `_test_skato()`

ACAT-V should also use the same resolved vector instead of recomputing Beta(MAF)
from `weights_beta`. This keeps ACAT-O's included ACAT-V component consistent
with the selected variant weighting scheme.

Backward compatibility:

- For `beta:*`, resolved vectors should numerically match current backend Beta
  weights.
- For `uniform`, resolved vectors should be all ones.
- Existing SKAT tests should continue to pass after reference calculations are
  updated to use the resolved vector where appropriate.

## R SKAT Backend

The R backend is deprecated and is not the primary target.

Two acceptable implementation choices:

1. If low risk, pass explicit weights to R SKAT through its supported `weights`
   argument when a non-Beta scheme is selected.
2. If not low risk, fail fast for score-column weights with `skat_backend="r"`
   and tell users to use the default Python backend.

The second option is acceptable for this issue because the Python backend is the
default and avoids adding fragile rpy2 behavior.

## Error Handling And Diagnostics

Hard errors:

- `score_column` without a resolved column name.
- Requested score column missing from the association DataFrame.
- Invalid `score_min >= score_max`.
- Negative `floor`.
- `ceiling <= 0`.
- `floor > ceiling`.
- Aligned score array length differs from `variant_mafs`.

Warnings:

- Score-column values outside 0-1 when no explicit range is provided.
- Missing or invalid scores are present.
- All score values for a gene are missing or invalid.
- R SKAT requested with score-column weights, if explicit R weights are not
  implemented.

Logging should include the score column name and counts. Per-gene warnings can
be log messages initially; adding output columns for weight diagnostics is not
required for the first implementation.

## Tests

Unit tests for `weights.py`:

- `column:<name>` dispatches to score-column weights.
- `score_column` dispatches when `weight_column` is provided.
- explicit score range normalization works.
- raw score-only mode works with `combine_with_beta=false`.
- missing scores use neutral fallback.
- missing scores use floor fallback.
- invalid params raise `ValueError`.
- score array length mismatch raises `ValueError`.
- existing specs remain unchanged.

Association stage tests:

- `variant_weight_column` is accepted by JSON config validation.
- CLI config propagation sets `variant_weight_column`.
- requested score column missing from DataFrame fails clearly.
- aligned `score_values` are attached to `gene_data`.
- one variant filtered by site missingness also filters the corresponding score.

Burden regression tests:

- `logistic_burden` with score-column weights changes the burden vector relative
  to `uniform`.
- `linear_burden` uses the same resolved weight vector.
- `combine_with_beta=false` produces burden equal to `geno @ normalized_scores`.

Python SKAT tests:

- SKAT uses explicit score weights instead of Beta-only weights.
- SKAT-O uses explicit score weights.
- ACAT-V uses explicit score weights.
- `beta:1,25` and `uniform` remain backward compatible.

Integration or stage-level regression:

- A small DataFrame containing `nephro_candidate_score` and per-sample GT columns
  can run `logistic_burden,skat` with `--variant-weights score_column` and
  produce results without falling back to Beta silently.

Documentation tests or docs review:

- Association guide lists `column:<name>` and `score_column`.
- Example command shows `nephro_candidate_score`.
- Config table includes `variant_weight_column`.

## Documentation

Update `docs/source/guides/association_testing.md`:

- Add a row for score-column weights.
- Explain the difference between variant-level weights and gene-level FDR
  weights.
- Document normalization parameters and defaults.
- Add the NCS example.
- Note that Fisher ignores variant weights and COAST currently uses its own
  category weights.

Update API docs if needed so `AssociationConfig.variant_weight_column` appears
in generated documentation.

## Acceptance Criteria

- Existing association workflows produce identical results with default
  `beta:1,25`.
- `--variant-weights column:nephro_candidate_score` runs for
  `logistic_burden`, `linear_burden`, and Python `skat`.
- `--variant-weights score_column --variant-weight-column nephro_candidate_score`
  is equivalent to `column:nephro_candidate_score`.
- Score arrays remain aligned after site missingness filtering.
- Missing or invalid score values are handled deterministically and logged.
- Python SKAT, SKAT-O, and ACAT-V all use the same resolved weight vector.
- Invalid configuration fails before silently running with the wrong weights.
- Documentation includes a working AGDE/NCS-style example.

## Self-Review

- Placeholder scan: no placeholder sections or unresolved markers remain.
- Scope check: the design is limited to score-column variant weights and does not
  change gene-level FDR priors, Fisher, or COAST semantics.
- Ambiguity check: canonical CLI forms, normalization defaults, missing-value
  behavior, and SKAT backend expectations are explicit.
- Consistency check: the design uses one resolved weight vector for every test
  that consumes variant-level weights, avoiding burden/SKAT divergence.
