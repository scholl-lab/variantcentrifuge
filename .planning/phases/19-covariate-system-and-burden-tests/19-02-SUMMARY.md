---
phase: 19-covariate-system-burden-tests
plan: "02"
title: "LogisticBurdenTest + LinearBurdenTest + Stage + CLI"
subsystem: association
tags: [logistic-regression, linear-regression, firth, burden-test, covariates, genotype-matrix, cli, statsmodels, numpy]

dependency-graph:
  requires:
    - 19-01  # covariates.py, genotype_matrix.py, weights.py, AssociationConfig phase 19 fields
  provides:
    - LogisticBurdenTest: "Logit with Firth N-R fallback; OR + 95% CI; warning codes"
    - LinearBurdenTest: "OLS regression; beta + SE + CI"
    - AssociationEngine-extended: "Registry now includes logistic_burden and linear_burden"
    - AssociationAnalysisStage-phase19: "Covariate loading + genotype matrix building + tiered sample size warnings"
    - CLI-phase19: "5 new args: --covariate-file, --covariates, --categorical-covariates, --trait-type, --variant-weights"
  affects:
    - 20-01  # SKAT backend will call same genotype matrix path
    - 22-01  # Column renaming (linear_burden_or -> linear_burden_beta) deferred to Phase 22

tech-stack:
  added: []
  patterns:
    - "firth-nr: self-contained Newton-Raphson (~130 lines) with Jeffreys prior + step-halving"
    - "separation-detection: mle_retvals['converged'] AND bse.max()>100 (both needed)"
    - "ci-ndarray: conf_int() returns ndarray (n_params,2) not DataFrame; use ci[1,0]/ci[1,1]"
    - "contingency-data-extension: regression keys added to per-gene dict; Fisher ignores them"
    - "tiered-sample-size-warnings: <10=error abort, <50=warn, <200=warn, ratio>1:20=warn"
    - "mac-check: per-gene MAC<5 returns empty matrix, regression reports NA"

key-files:
  created:
    - variantcentrifuge/association/tests/logistic_burden.py
    - variantcentrifuge/association/tests/linear_burden.py
  modified:
    - variantcentrifuge/association/tests/__init__.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/cli.py
    - tests/unit/test_association_stage.py

decisions:
  - id: IMPL-12
    decision: "LogisticBurdenTest builds design matrix inline (not in genotype_matrix.py)"
    rationale: "Firth fallback and separation checks are logistic-specific; keeping in logistic_burden.py avoids coupling"
  - id: IMPL-13
    decision: "Tiered sample size check at n_cases<10 aborts with logger.error() not exception"
    rationale: "Stage returns context cleanly; exception would crash pipeline; logger.error signals severity"
  - id: IMPL-14
    decision: "linear_burden effect_size = beta (not OR); engine column named *_or contains beta"
    rationale: "Column renaming deferred to Phase 22 per ARCH decision; no semantic change needed now"

metrics:
  duration: "~9 minutes"
  completed: "2026-02-20"
  tasks: 2
  tests-passing: 72
---

# Phase 19 Plan 02: LogisticBurdenTest + LinearBurdenTest + Stage + CLI Summary

**One-liner:** Logit burden test with Firth N-R fallback (OR + CI + warning codes) + OLS linear burden test + stage genotype matrix builder + 5 new CLI args for regression-based association.

## What Was Built

### variantcentrifuge/association/tests/logistic_burden.py (new)

`LogisticBurdenTest(AssociationTest)`:

- `name` = `"logistic_burden"`, `check_dependencies()` lazily imports statsmodels
- `run()` flow:
  1. Check for `genotype_matrix` in `contingency_data`; return `NO_GENOTYPE_MATRIX` if absent
  2. Call `get_weights(mafs, config.variant_weights)` for Beta(MAF;1,25) or uniform weights
  3. Compute `burden = geno @ weights` (n_samples,)
  4. Build design matrix `[intercept, burden, covariates...]` via `sm.add_constant` + `np.column_stack`
  5. Pre-flight checks: `PERFECT_SEPARATION`, `QUASI_SEPARATION`, `ZERO_CARRIERS_ONE_GROUP`, `LOW_CARRIER_COUNT`
  6. Fit `sm.Logit(...).fit(disp=False, maxiter=100)` inside `warnings.catch_warnings(record=True)`
  7. Detect separation: `not mle_retvals['converged']` OR `bse.max() > 100.0`
  8. Firth fallback via `_firth_logistic()`; if Firth returns None → `FIRTH_CONVERGE_FAIL` + p_value=None
  9. Extract index 1 (not 0): `beta=params[1]`, `se=bse[1]`, `p_value=pvalues[1]`, `ci=conf_int()[1,:]`
  10. Return `effect_size=exp(beta)`, `ci_lower=exp(ci_lower)`, `ci_upper=exp(ci_upper)`

Self-contained `_firth_logistic(y, x, max_iter, tol)` Newton-Raphson implementation:
- Jeffreys prior penalty: penalized score `U = X^T(y - pi + h_diag*(0.5 - pi))`
- Hat matrix diagonal computed memory-efficiently: `diag(W^0.5 X V X^T W^0.5)`
- Step-halving: tries 10 halvings to ensure penalized log-likelihood increases
- Returns `FirthResult` namespace with `.params`, `.bse`, `.pvalues`, `.conf_int()`
- Returns None on `LinAlgError` (singular Fisher information)

`_firth_loglik(beta, y, x)`: penalized log-likelihood = log-lik + 0.5 * log(det(I))

### variantcentrifuge/association/tests/linear_burden.py (new)

`LinearBurdenTest(AssociationTest)`:

- `name` = `"linear_burden"`, `check_dependencies()` lazily imports statsmodels
- `run()` flow (simpler — no Firth needed):
  1. Check `genotype_matrix` presence
  2. Compute weighted burden score
  3. Build design matrix
  4. Fit `sm.OLS(phenotype, design).fit()`
  5. Extract index 1: `beta`, `se`, `p_value`, `ci = conf_int()[1, :]` (ndarray, not DataFrame)
  6. Return `effect_size=beta`, `ci_lower`, `ci_upper`

### variantcentrifuge/association/tests/__init__.py (updated)

Lazy loading via `__getattr__` for `LogisticBurdenTest` and `LinearBurdenTest`. Updated `__all__`.

### variantcentrifuge/association/engine.py (updated)

`_build_registry()` now imports and registers:
- `"logistic_burden": LogisticBurdenTest`
- `"linear_burden": LinearBurdenTest`

Registry grows from 1 to 3 entries. `AssociationEngine.from_names(["fisher", "logistic_burden"], config)` works immediately.

### variantcentrifuge/cli.py (updated)

Five new arguments added to the stats argument group after `--skat-backend`:

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--covariate-file` | str | None | TSV/CSV covariate file path |
| `--covariates` | str | None | Comma-separated column names to use |
| `--categorical-covariates` | str | None | Force columns to categorical |
| `--trait-type` | choice | "binary" | "binary" or "quantitative" |
| `--variant-weights` | str | "beta:1,25" | Weight scheme spec |

Updated `--association-tests` help text to include `logistic_burden, linear_burden`.

New config mapping: all 5 args mapped to `cfg[...]` in the config section.

New validations:
- `--covariate-file` requires `--perform-association`
- `--trait-type quantitative` + `logistic_burden` → `parser.error()`

### variantcentrifuge/stages/analysis_stages.py (updated)

`AssociationAnalysisStage._process()` extended:

1. **Pass Phase 19 fields to AssociationConfig**: `covariate_file`, `covariate_columns`, `categorical_covariates`, `trait_type`, `variant_weights`

2. **Tiered sample size warnings** (CONTEXT.md spec):
   - `n_cases < 10` → `logger.error()` + `return context` (abort)
   - `n_cases < 50` → `logger.warning()` (no practical power)
   - `n_cases < 200` → `logger.warning()` (underpowered for SKAT)
   - `n_controls / n_cases > 20` → `logger.warning()` (Type I error risk)

3. **Covariate loading once**: `load_covariates()` called before per-gene loop when `covariate_file` set

4. **Phenotype vector once**: `np.array([1 if s in case_set else 0 for s in vcf_samples_list])`

5. **Per-gene genotype matrix augmentation**: when `logistic_burden`, `linear_burden`, `skat`, or `skat_python` in test_names:
   - `build_genotype_matrix()` called per gene
   - Sample mask applied to phenotype and covariate matrices
   - MAC < 5 → empty matrix (regression reports NA, not error)
   - Keys added to gene_data dict: `genotype_matrix`, `variant_mafs`, `phenotype_vector`, `covariate_matrix`, `vcf_samples`
   - Backward compatible: FisherExactTest ignores these keys

## Decisions Made

| ID | Decision | Rationale |
|----|----------|-----------|
| IMPL-12 | Design matrix built inline in LogisticBurdenTest, not in genotype_matrix.py | Firth + separation checks are logistic-specific; keeps coupling clean |
| IMPL-13 | n_cases<10 aborts with logger.error() not exception | Stage returns context cleanly; exception would crash pipeline |
| IMPL-14 | linear_burden effect_size=beta; engine column `linear_burden_or` contains beta | Column renaming deferred to Phase 22; no behavioral change needed now |

## Verification Results

All plan-specified verifications passed:

```
python -c "from variantcentrifuge.association.engine import _build_registry; r = _build_registry(); assert len(r) == 3"  # OK
python -c "from variantcentrifuge.association.engine import AssociationEngine; from variantcentrifuge.association.base import AssociationConfig; engine = AssociationEngine.from_names(['fisher', 'logistic_burden'], AssociationConfig())"  # OK
python -c "from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest; ..."  # OK
make lint  # 0 errors
pytest tests/unit/test_association_{base,engine,fisher,stage}.py  # 72 passed
```

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed N806, E501, SIM102, RUF059 lint violations in new files**

- **Found during:** Task 1 (lint check)
- **Issue:** Variable names `n, k = x.shape` → `_n, k` (RUF059); line length violations (E501); nested if that could be combined (SIM102)
- **Fix:** Renamed unused variable to `_n`; refactored covariate check to multi-line `has_covariates` variable; combined nested if with `and`
- **Files modified:** `logistic_burden.py`, `linear_burden.py`

**2. [Rule 1 - Bug] Updated test fixtures from 2 cases to 10 cases**

- **Found during:** Task 2 verification (test run)
- **Issue:** `case_control_config` fixture used `["CASE1", "CASE2"]` (2 cases). New tiered sample size check in `AssociationAnalysisStage._process()` aborts when `n_cases < 10` with `logger.error()`. This caused 1 test to fail and 5 more to be at risk.
- **Fix:** Updated `case_control_config` fixture and all inline test configs with `perform_association=True` to use `[f"CASE{i}" for i in range(1, 11)]` (10 cases, 10 controls). The sample size abort is correct per CONTEXT.md spec — tests needed to be representative of real usage.
- **Files modified:** `tests/unit/test_association_stage.py`

## Task Commit Table

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Implement LogisticBurdenTest and LinearBurdenTest | 90cfc98 | logistic_burden.py, linear_burden.py, __init__.py, engine.py |
| 2 | Wire covariate loading and genotype matrix into stage + CLI args | 1e86202 | analysis_stages.py, cli.py, test_association_stage.py |

## Next Phase Readiness

Plan 20-01 (R SKAT Backend) can proceed:
- `build_genotype_matrix()` already called per-gene in `AssociationAnalysisStage._process()` when `skat` is in test_names
- `phenotype_vector` and `covariate_matrix` already available in `contingency_data`
- `AssociationConfig.trait_type` already gates binary vs quantitative path
- All infrastructure from Plans 19-01 and 19-02 is in place
