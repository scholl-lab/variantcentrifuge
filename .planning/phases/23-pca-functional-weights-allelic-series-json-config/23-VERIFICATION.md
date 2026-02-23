---
phase: 23-pca-functional-weights-allelic-series-json-config
verified: 2026-02-22T00:20:14Z
status: passed
score: 18/18 must-haves verified
---

# Phase 23: PCA Integration, Functional Weights, Allelic Series, and JSON Config Verification Report

**Phase Goal:** Users can supply a pre-computed PCA file or trigger AKT-based PCA computation as a pipeline stage, apply CADD/REVEL functional variant weights, run the COAST allelic series test, and specify all association options via a JSON config file — with optional matplotlib QQ plot generation for users who have it installed.
**Verified:** 2026-02-22T00:20:14Z
**Status:** passed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| #  | Truth                                                                                 | Status     | Evidence                                                                                             |
|----|---------------------------------------------------------------------------------------|------------|------------------------------------------------------------------------------------------------------|
| 1  | User can load a PLINK .eigenvec file and have PCs merged as covariates                | VERIFIED   | `load_pca_file()` handles `plink_header`/`plink_nohdr` formats; `merge_pca_covariates()` exists and is called in `analysis_stages.py:2496-2499` |
| 2  | User can load an AKT output file and have PCs merged as covariates                    | VERIFIED   | `_detect_pca_format()` identifies `akt_or_generic` format; `load_pca_file()` parses it with no header |
| 3  | User can specify `--pca-tool akt` and AKT runs as a pipeline stage                   | VERIFIED   | `PCAComputationStage` in `processing_stages.py:2074`; registered in `stage_registry.py:481`; `--pca-tool` CLI arg wired at `cli.py:514` |
| 4  | Missing AKT binary raises ToolNotFoundError (hard error, not skip)                   | VERIFIED   | `processing_stages.py:2106` raises `ToolNotFoundError("akt", self.name)`; test `test_pca.py:380-387` verifies |
| 5  | Requesting >20 PCA components triggers a logged warning                               | VERIFIED   | `pca.py:130-135` logs warning when `n_components > 20`; `test_pca.py:255-270` verifies |
| 6  | PCA sample IDs are aligned to VCF sample order with validation                        | VERIFIED   | `pca.py:236-243` raises `ValueError` listing missing samples; `pca.py:257` reindexes to `vcf_samples` order |
| 7  | User can specify `--variant-weights cadd` and get CADD-normalized weights applied     | VERIFIED   | `weights.py:183-220` implements `cadd_weights()`; `get_weights()` dispatches at line 381; wired in `logistic_burden.py:331-337` and `linear_burden.py:153-159` |
| 8  | User can specify `--variant-weights revel` and get REVEL-based weights applied        | VERIFIED   | `weights.py:223-257` implements `revel_weights()`; dispatched by `get_weights()` at line 390 |
| 9  | User can specify `--variant-weights combined` and get Beta(MAF) x CADD weights       | VERIFIED   | `weights.py:260-300` implements `combined_weights()` preferring CADD over REVEL; dispatched at line 398 |
| 10 | LoF/missense variants missing scores receive functional weight 1.0                    | VERIFIED   | `weights.py:216` uses `np.where(nan_mask, 1.0, ...)` for CADD; `_log_missing_score_counts()` logs per-category counts |
| 11 | Functional weight schemes produce different p-values than uniform on same data        | VERIFIED   | `test_functional_weights.py:407-430` explicitly tests this (TestNumericalDifference class) |
| 12 | COAST test classifies variants into BMV/DMV/PTV categories                            | VERIFIED   | `allelic_series.py:101-238` implements `classify_variants()` with codes 3=PTV, 2=DMV, 1=BMV, 0=unclassified |
| 13 | COAST wraps R AllelicSeries::COAST() via rpy2                                         | VERIFIED   | `allelic_series.py:243` class `COASTTest`; `check_dependencies()` at line 289-317 imports `rpy2.robjects` and `AllelicSeries` R package |
| 14 | Engine registers "coast" test and ACAT-O integrates coast p-values                    | VERIFIED   | `engine.py:71` maps `"coast": COASTTest`; `analysis_stages.py:2406` adds `"coast"` to `needs_regression` tuple |
| 15 | User can add an 'association' section to config.json and have all fields applied       | VERIFIED   | `_build_assoc_config_from_context()` at `analysis_stages.py:2161` reads `context.config["association"]` dict and applies all fields |
| 16 | CLI flags override JSON config values                                                 | VERIFIED   | `_get()` closure in `_build_assoc_config_from_context()` implements CLI > JSON > default precedence; `test_json_config.py:214-238` verifies override |
| 17 | Invalid JSON config keys produce a clear error listing all unknown keys               | VERIFIED   | `_validate_association_config_dict()` at line 2095-2097 collects all unknown keys and raises with complete list; `test_json_config.py:86-93` verifies |
| 18 | When matplotlib is absent, pipeline completes without error and logs INFO              | VERIFIED   | `diagnostics.py:255-256` catches `ImportError` and returns `False`; `test_qq_plot.py:140-166` verifies no exception and no file written |

**Score:** 18/18 truths verified

---

### Required Artifacts

| Artifact                                                                   | Expected                                              | Status     | Details                                              |
|---------------------------------------------------------------------------|-------------------------------------------------------|------------|------------------------------------------------------|
| `variantcentrifuge/association/pca.py`                                     | `load_pca_file`, `merge_pca_covariates`               | VERIFIED   | 306 lines, both functions exported, no stubs         |
| `tests/unit/test_pca.py`                                                   | min 100 lines, 35 test methods                        | VERIFIED   | 453 lines, 35 test methods, all pass                 |
| `variantcentrifuge/association/weights.py`                                 | `cadd_weights`, `revel_weights`, `combined_weights`, `get_weights` | VERIFIED | 411 lines, all 4 functions present and exported |
| `tests/unit/test_functional_weights.py`                                    | min 80 lines, 45 test methods                         | VERIFIED   | 492 lines, 45 test methods, all pass                 |
| `variantcentrifuge/association/tests/allelic_series.py`                    | `COASTTest`, `classify_variants`                      | VERIFIED   | 775 lines, both symbols present                      |
| `tests/unit/test_coast.py`                                                 | 43 test methods                                       | VERIFIED   | 1123 lines, 43 test methods, all pass                |
| `variantcentrifuge/stages/analysis_stages.py`                              | `_build_assoc_config_from_context`, `VALID_ASSOCIATION_KEYS` | VERIFIED | Both symbols present; called in `_process()` at line 2371 |
| `variantcentrifuge/association/diagnostics.py`                             | `write_qq_plot` with lazy matplotlib + Agg backend    | VERIFIED   | `write_qq_plot()` at line 225; lazy import + `matplotlib.use("Agg")` at lines 251-253; called from `write_diagnostics()` at line 433 |
| `tests/unit/test_json_config.py`                                           | min 80 lines, 24 test methods                         | VERIFIED   | 412 lines, 24 test methods, all pass                 |
| `tests/unit/test_qq_plot.py`                                               | min 40 lines, 11 test methods                         | VERIFIED   | 286 lines, 11 test methods, all pass                 |

---

### Key Link Verification

| From                                           | To                                          | Via                                                                   | Status | Details                                                        |
|------------------------------------------------|---------------------------------------------|-----------------------------------------------------------------------|--------|----------------------------------------------------------------|
| `pca.py:load_pca_file`                         | `analysis_stages.py:AssociationAnalysisStage` | `from ..association.pca import load_pca_file, merge_pca_covariates` | WIRED  | `analysis_stages.py:2493-2499` imports and calls both functions |
| `processing_stages.py:PCAComputationStage`     | `stage_registry.py`                          | `register_stage(PCAComputationStage, ...)` at line 481                | WIRED  | Stage registered as `pca_computation`                          |
| `weights.py:get_weights`                       | `logistic_burden.py` / `linear_burden.py`    | `get_weights()` call passing `cadd_scores`, `revel_scores`, `variant_effects` from `contingency_data` | WIRED | `logistic_burden.py:331-337`, `linear_burden.py:153-159` |
| `analysis_stages.py`                           | CADD/REVEL annotation columns                | Site-filter-aligned extraction in genotype matrix loop `2675-2685`   | WIRED  | `gene_data["cadd_scores"]`, `gene_data["revel_scores"]`, `gene_data["variant_effects"]` populated |
| `allelic_series.py:COASTTest`                  | `association/engine.py`                      | `_build_registry()` maps `"coast": COASTTest` at line 71             | WIRED  | Engine lazy-loads via `tests/__init__.py:34-37`                |
| `allelic_series.py:COASTTest`                  | `analysis_stages.py:needs_regression`        | `"coast"` in `needs_regression` tuple at line 2406                   | WIRED  | Triggers regression code path when coast test selected         |
| `analysis_stages.py:gene_data["gene_df"]`      | `allelic_series.py:COASTTest.run()`           | `gene_data["gene_df"] = gene_df.reset_index(drop=True)` at line 2697 | WIRED  | COASTTest accesses per-variant annotations via `gene_data["gene_df"]` |
| `_build_assoc_config_from_context()`           | `context.config["association"]`              | Reads `cfg.get("association", {})` at line 2188                      | WIRED  | JSON section applied with CLI > JSON > default precedence      |
| `_validate_association_config_dict()`          | `_build_assoc_config_from_context()`         | Called at line 2190 when `json_assoc` is non-empty                   | WIRED  | Validation runs before `AssociationConfig` construction        |
| `diagnostics.py:write_qq_plot`                 | `write_diagnostics()`                        | Called at line 433 with `qq_combined` data and `qq_plot_path`        | WIRED  | `write_diagnostics()` calls `write_qq_plot()` after TSV output  |

---

### Anti-Patterns Found

No anti-patterns detected in any of the four implementation files:
- `variantcentrifuge/association/pca.py` — no TODO/FIXME/placeholder, no empty returns, real implementation throughout
- `variantcentrifuge/association/weights.py` — no TODO/FIXME/placeholder, all weight functions fully implemented
- `variantcentrifuge/association/tests/allelic_series.py` — no TODO/FIXME, real rpy2 wrapper with proper BMV/DMV/PTV classification
- `variantcentrifuge/association/diagnostics.py` — no TODO/FIXME, `write_qq_plot()` fully implemented with graceful matplotlib-absent path

---

### Test Execution Results

All 158 Phase 23 unit tests pass:
- `tests/unit/test_pca.py` — 35 passed
- `tests/unit/test_functional_weights.py` — 45 passed
- `tests/unit/test_coast.py` — 43 passed
- `tests/unit/test_json_config.py` — 24 passed
- `tests/unit/test_qq_plot.py` — 11 passed

---

### Human Verification Required

| Test                               | What to do                                                                                   | Expected                                                                 | Why human                                                       |
|------------------------------------|----------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|-----------------------------------------------------------------|
| AKT PCA end-to-end                 | Run `variantcentrifuge --pca-tool akt --perform-association` with `akt` installed            | AKT subprocess executes, PCA file written, loaded and merged as covariates in output | Requires live AKT binary and VCF input; can't mock subprocess meaningfully |
| COAST with real R environment      | Run with `--association-tests coast` with R + AllelicSeries installed                       | COAST p-values appear in per-gene output; ACAT-O omnibus p includes COAST component | Requires live R + AllelicSeries package installation            |
| QQ plot PNG visual inspection      | Run with `--diagnostics-output /tmp/diag` with matplotlib installed                          | `qq_plot.png` generated; axes labelled; diagonal line visible; points follow expected pattern | Visual appearance cannot be verified programmatically           |
| CADD/REVEL from real annotation    | Run with `--variant-weights cadd` on VCF with dbNSFP CADD scores                            | Weights differ from beta default; per-gene p-values change; log shows score extraction | Requires annotated VCF with real CADD/REVEL columns            |

---

### Validation Timing Note

The `_validate_association_config_dict()` function runs within `AssociationAnalysisStage._process()` — not at CLI parsing or pipeline initialization time. The plan's claim of "validated at startup before any processing" is imprecise: validation occurs before any association computation begins but after VCF preprocessing stages complete. This is architecturally correct (the stage fails fast and clearly) but is not true pipeline-startup validation. This is an informational note, not a gap.

---

## Summary

All 18 must-have truths for Phase 23 are verified. All 10 required artifacts exist, are substantive (no stubs), and are correctly wired into the pipeline:

- **Plan 23-01 (PCA):** `pca.py` implements three-format auto-detection, sample alignment with ValueError on missing samples, >20 component warning, and `PCAComputationStage` with hard `ToolNotFoundError` for missing AKT. Fully wired in `analysis_stages.py`.
- **Plan 23-02 (Functional Weights):** `weights.py` implements `cadd_weights`, `revel_weights`, `combined_weights` with NaN-safe parsing, type-aware missing-score fallback logging, and backward-compatible `get_weights()` extension. Wired end-to-end through stage annotation extraction to both burden test implementations.
- **Plan 23-03 (COAST):** `allelic_series.py` implements `classify_variants()` with 6-candidate SIFT/PolyPhen column detection, `COASTTest` wrapping R AllelicSeries via rpy2, and returns `p_value=None` on missing categories. Registered in engine and in `needs_regression` tuple. `gene_df` stored in `gene_data` for annotation access.
- **Plan 23-04 (JSON Config + QQ Plot):** `VALID_ASSOCIATION_KEYS` frozenset, fail-fast `_validate_association_config_dict()` collecting all errors, `_build_assoc_config_from_context()` with CLI > JSON > default precedence. `write_qq_plot()` uses lazy matplotlib import with Agg backend; headless-safe; returns `False` gracefully when matplotlib absent.

---

_Verified: 2026-02-22T00:20:14Z_
_Verifier: Claude (gsd-verifier)_
