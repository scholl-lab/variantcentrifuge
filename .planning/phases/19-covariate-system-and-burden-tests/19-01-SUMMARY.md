---
phase: 19-covariate-system-burden-tests
plan: "01"
title: "Covariate System and Data Infrastructure"
subsystem: association
tags: [covariates, genotype-matrix, weights, burden-tests, numpy, pandas, scipy]

dependency-graph:
  requires:
    - 18-04  # AssociationConfig ABC, TestResult, AssociationEngine from Phase 18
  provides:
    - load_covariates: "Covariate loading with sample alignment, one-hot encoding, multicollinearity check"
    - parse_gt_to_dosage: "VCF GT string to (dosage, is_multi_allelic) tuple"
    - build_genotype_matrix: "Two-layer missing strategy, imputed (n_samples, n_variants) matrix"
    - beta_maf_weights: "Beta(MAF;1,25) weights upweighting rare variants (SKAT convention)"
    - uniform_weights: "All-ones weight vector for unweighted burden tests"
    - AssociationConfig-extended: "8 new fields for covariate/trait/weight configuration"
  affects:
    - 19-02  # LogisticBurdenTest, LinearBurdenTest will call these functions
    - 20-01  # SKAT backend will reuse genotype_matrix.py

tech-stack:
  added: []
  patterns:
    - "reindex-align: df.reindex(vcf_samples) + NaN assertion for covariate alignment"
    - "two-layer-missing: pre-filter variants/samples then mean imputation"
    - "maf-from-all: MAF from all samples combined (never stratified) per SKAT R convention"
    - "multi-allelic-flag: parse_gt_to_dosage returns (dosage, bool) tuple"

key-files:
  created:
    - variantcentrifuge/association/covariates.py
    - variantcentrifuge/association/genotype_matrix.py
    - variantcentrifuge/association/weights.py
  modified:
    - variantcentrifuge/association/base.py

decisions:
  - id: IMPL-09
    decision: "parse_gt_to_dosage returns (int|None, bool) not int|None"
    rationale: "Multi-allelic flag needed to emit 'run bcftools norm' warning without second parse pass"
  - id: IMPL-10
    decision: "load_covariates returns (np.ndarray, list[str]) tuple"
    rationale: "Column names returned alongside matrix for diagnostics; callers can ignore second element"
  - id: IMPL-11
    decision: "build_genotype_matrix: sample_mask is list[bool], all samples remain in G"
    rationale: "Callers (logistic burden test) decide whether to exclude high-missing samples; matrix returned for all"

metrics:
  duration: "~4 minutes"
  completed: "2026-02-20"
  tasks: 2
  tests-passing: 51
---

# Phase 19 Plan 01: Covariate System and Data Infrastructure Summary

**One-liner:** Covariate file loading with reindex alignment + VCF GT matrix builder with two-layer missing strategy + Beta(MAF;1,25) weights + 8 new AssociationConfig fields.

## What Was Built

Three new modules and one extended module providing the data infrastructure
that Phase 19's burden test implementations (Plan 02) will consume:

### variantcentrifuge/association/base.py (extended)

Added 8 new fields to `AssociationConfig` — all optional with backward-compatible
defaults. No existing fields or docstrings changed.

New fields:
- `covariate_file: str | None = None`
- `covariate_columns: list[str] | None = None`
- `categorical_covariates: list[str] | None = None`
- `trait_type: str = "binary"`
- `variant_weights: str = "beta:1,25"`
- `missing_site_threshold: float = 0.10`
- `missing_sample_threshold: float = 0.80`
- `firth_max_iter: int = 25`

### variantcentrifuge/association/covariates.py (new)

`load_covariates(filepath, vcf_samples, covariate_columns, categorical_columns) -> (np.ndarray, list[str])`

- Auto-detects delimiter from extension (`.tsv`/`.tab` -> tab, `.csv` -> comma,
  else `csv.Sniffer` on first 2KB, fallback tab)
- `pd.read_csv(filepath, sep=sep, index_col=0)` — first column = sample ID
- Column selection if `covariate_columns` not None
- VCF samples missing from covariate file -> `ValueError` listing all missing IDs
- Extra covariate samples not in VCF -> logged WARNING (up to 5 IDs shown)
- `df.reindex(vcf_samples)` + `assert not df_aligned.isnull().any().any()`
- Auto-detect categorical: non-numeric AND `nunique() <= 5`
- `pd.get_dummies(..., drop_first=True, dtype=float)` for one-hot encoding
- `np.linalg.cond(X)` — WARNING if condition number > 1000
- Returns `(float64 ndarray, column_names)` — shape `(n_samples, k)`

### variantcentrifuge/association/genotype_matrix.py (new)

`parse_gt_to_dosage(gt) -> tuple[int | None, bool]`

- Handles all GT edge cases: `0/0` -> `(0, False)`, `0/1` -> `(1, False)`,
  `1/1` -> `(2, False)`, `./.` -> `(None, False)`, `1/2` -> `(1, True)`
- Phased genotypes (`|`) treated identically to unphased (`/`)
- Partial missing (`./1`) -> `(None, False)`
- Multi-allelic (any allele > 1) -> dosage 1, flag True
- NOT a delegation to `gene_burden._gt_to_dosage()` (which returns 0 for `1/2`)

`build_genotype_matrix(gene_df, vcf_samples, gt_columns, ...) -> (G, mafs, sample_mask, warnings_list)`

- Parses all GT values into `(n_variants, n_samples)` raw dosage matrix
- **Layer 1a:** Remove variants with > `missing_site_threshold` (10%) fraction missing
- **Layer 1b:** Flag samples with > `missing_sample_threshold` (80%) fraction missing
- **Layer 2:** Compute MAF from ALL samples (never stratified), impute: `round(2*MAF)` for binary, `2*MAF` for quantitative
- Returns fully imputed matrix, `np.isnan(G).sum() == 0` guaranteed
- Emits structured warning if multi-allelic genotypes detected

### variantcentrifuge/association/weights.py (new)

- `beta_maf_weights(mafs, a=1.0, b=25.0)`: `scipy.stats.beta.pdf` with `np.clip(mafs, 1e-8, 1-1e-8)`
- `uniform_weights(n_variants)`: `np.ones(n_variants)`
- `get_weights(mafs, weight_spec)`: parses `"beta:a,b"` or `"uniform"`, raises `ValueError` for unknown

## Decisions Made

| ID | Decision | Rationale |
|----|----------|-----------|
| IMPL-09 | `parse_gt_to_dosage` returns `(int\|None, bool)` | Multi-allelic flag needed for warning without second parse |
| IMPL-10 | `load_covariates` returns `(ndarray, list[str])` | Column names for diagnostics; callers can ignore |
| IMPL-11 | `sample_mask` is `list[bool]`, all samples remain in `G` | Callers decide exclusion; matrix covers all samples |

## Verification Results

All plan-specified verifications passed:

```
python -c "from variantcentrifuge.association.base import AssociationConfig; c = AssociationConfig(); assert c.trait_type == 'binary'"  # OK
python -c "from variantcentrifuge.association.covariates import load_covariates"  # OK
python -c "from variantcentrifuge.association.genotype_matrix import parse_gt_to_dosage, build_genotype_matrix"  # OK
python -c "from variantcentrifuge.association.weights import beta_maf_weights, uniform_weights, get_weights"  # OK
make lint  # 0 errors
pytest tests/unit/test_association_{base,engine,fisher}.py  # 51 passed
```

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Renamed uppercase matrix variables to avoid N806 lint violation**

- **Found during:** Task 2 (lint check)
- **Issue:** Variable names `X` and `G` trigger `N806` (variable in function should be lowercase)
- **Fix:** Renamed `X` to `covariate_matrix` in `covariates.py`; renamed `G` to `geno` in `genotype_matrix.py`
- **Files modified:** `covariates.py`, `genotype_matrix.py`

**2. [Rule 1 - Bug] Fixed B904 lint error in weights.py**

- **Found during:** Task 1 (lint check)
- **Issue:** `raise ValueError` inside `except` clause without `from err`
- **Fix:** Changed to `raise ValueError(...) from err`
- **Files modified:** `weights.py`

**3. [Rule 1 - Bug] Changed `from typing import Sequence` to `from collections.abc import Sequence`**

- **Found during:** Task 1 and Task 2 (lint check UP035)
- **Issue:** Python 3.10+ prefers `collections.abc` imports
- **Fix:** Updated both `covariates.py` and `genotype_matrix.py`
- **Files modified:** `covariates.py`, `genotype_matrix.py`

## Task Commit Table

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Extend AssociationConfig + create covariates.py + weights.py | ac79a68 | base.py, covariates.py, weights.py |
| 2 | Create genotype_matrix.py | 5e3761b | genotype_matrix.py |

## Next Phase Readiness

Plan 19-02 (LogisticBurdenTest + LinearBurdenTest) can proceed immediately:
- `load_covariates()` is ready to be called from `AssociationAnalysisStage._process()`
- `build_genotype_matrix()` is ready to be called per-gene from `LogisticBurdenTest.run()`
- `get_weights()` is ready to be called with `AssociationConfig.variant_weights` spec
- `AssociationConfig.trait_type` gates logistic vs linear path
