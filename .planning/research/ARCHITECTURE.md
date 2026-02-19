# Architecture Patterns: Modular Rare Variant Association Framework

**Project:** variantcentrifuge — Association Framework Integration
**Milestone:** v0.15.0
**Researched:** 2026-02-19
**Confidence:** HIGH (based on direct codebase reading)

---

## Executive Summary

The association framework integrates into the existing stage-based pipeline through four key decisions:

1. **Two new stages, not one** — `AssociationAnalysisStage` replaces `GeneBurdenAnalysisStage` for users requesting association tests. The old stage becomes a shim that delegates to `AssociationEngine.run_fisher_only()` for backward compatibility.

2. **PCA as an optional pre-processing stage** — `PCAComputationStage` runs early, after VCF extraction but before DataFrame analysis. Its output (eigenvectors) flows through `PipelineContext` to `AssociationAnalysisStage`.

3. **rpy2 isolation via factory pattern** — The R backend is never imported at module level. `get_skat_backend()` is called lazily inside `AssociationAnalysisStage._process()`, keeping the import boundary explicit and the package importable without R.

4. **Association results extend the existing output pattern** — `context.gene_burden_results` is the established field for gene-level DataFrames; association results augment it. Output stages check `context.config.get("perform_association")` alongside `perform_gene_burden` to decide what sheets to write.

---

## Current Architecture Baseline

### Pipeline Flow (as-is)

```
CLI (cli.py) → argparse.Namespace → config dict
                                         |
                              run_refactored_pipeline()
                                         |
                              PipelineRunner.run(stages, context)
                                         |
  SETUP (7 stages)
    ConfigurationLoadingStage
    PhenotypeLoadingStage
    ScoringConfigLoadingStage
    PedigreeLoadingStage
    AnnotationConfigLoadingStage
    SampleConfigLoadingStage
    PhenotypeCaseControlAssignmentStage
                                         |
  PROCESSING (5-7 stages)
    GeneBedCreationStage
    VariantExtractionStage / ParallelVariantExtractionStage
    BCFToolsPrefilterStage (optional)
    SnpSiftFilterStage
    FieldExtractionStage
    GenotypeReplacementStage (no-op since v0.13, kept for dep graph)
    PhenotypeIntegrationStage
    ExtraColumnRemovalStage
                                         |
  ANALYSIS (8-10 stages)
    DataFrameLoadingStage           ← context.current_dataframe set here
    CustomAnnotationStage
    InheritanceAnalysisStage
    VariantScoringStage
    GenotypeFilterStage
    StatisticsGenerationStage
    VariantAnalysisStage
    GeneBurdenAnalysisStage         ← context.gene_burden_results set here
    ChunkedAnalysisStage (alt path)
                                         |
  OUTPUT (5-7 stages)
    VariantIdentifierStage
    FinalFilteringStage
    PseudonymizationStage
    TSVOutputStage
    ExcelReportStage                ← reads context.gene_burden_results
    HTMLReportStage
    IGVReportStage
```

### PipelineContext Fields Relevant to Association

```python
# From context.py (actual fields)
vcf_samples: list[str]           # Sample IDs from VCF — ORDER MATCHES GT COLUMNS
pedigree_data: dict | None       # {sample_id: {father, mother, sex, affection}}
phenotype_data: dict | None      # {sample_id: {hpo_terms, ...}}
current_dataframe: pd.DataFrame  # Variant × feature matrix (post-analysis)
variants_df: pd.DataFrame        # Optimized variant DataFrame (alt pass-through)
gene_burden_results: pd.DataFrame | None  # Gene-level results table
config: dict[str, Any]          # Merged CLI + file config
  config["case_samples"]: list[str]       # Set by PhenotypeCaseControlAssignmentStage
  config["control_samples"]: list[str]    # Set by PhenotypeCaseControlAssignmentStage
  config["perform_gene_burden"]: bool
stage_results: dict[str, Any]    # Arbitrary per-stage result storage
```

### GeneBurdenAnalysisStage: What Exists

Reading the actual implementation (`analysis_stages.py:1810-2050`), the current stage:

- **Hard dependencies:** `{"dataframe_loading", "sample_config_loading"}`
- **Soft dependencies:** `{"custom_annotation"}`
- **Guard:** Checks `context.config.get("perform_gene_burden")` first; returns early if False
- **Input:** `context.current_dataframe` with `GENE` column and either packed `GT` or per-sample `GEN_N__GT` columns
- **Core call:** `perform_gene_burden_analysis(df, burden_config, case_samples, control_samples, vcf_samples)`
- **Output:** Sets `context.gene_burden_results` and writes to file path in `context.config["gene_burden_output"]`
- **Chunked processing awareness:** Defers if `use_chunked_processing=True` and chunks not complete

The stage has a `_handle_checkpoint_skip` method that restores the gene burden output path from disk — new AssociationAnalysisStage must implement the same pattern.

### gene_burden.py: Reusable Components

The Fisher-specific logic in `gene_burden.py` is structured as pure functions:

```python
perform_gene_burden_analysis(df, cfg, case_samples, control_samples, vcf_samples) -> pd.DataFrame
  → _aggregate_gene_burden_from_columns(df, case_samples, control_samples, vcf_samples, gt_columns)
  → _aggregate_gene_burden_from_gt(df, case_samples, control_samples)
  → _aggregate_gene_burden_legacy(df)
  → _compute_or_confidence_interval(...)
  → _gt_to_dosage(gt) -> int
  → _find_gt_columns(df) -> list[str]
```

These aggregation functions extract the **genotype matrix** per gene — exactly what the association engine needs. They are shared utilities, not `GeneBurdenAnalysisStage`-specific.

---

## Recommended Architecture

### New Components Overview

```
variantcentrifuge/
  association/                          # NEW package
    __init__.py
    engine.py                           # AssociationEngine orchestrator
    config.py                           # AssociationConfig dataclass
    results.py                          # TestResult, AssociationResults container
    correction.py                       # Multiple testing (refactored from gene_burden.py)
    covariates.py                       # Covariate loading + PCA merging
    pca.py                              # PCA file parser + AKT/PLINK wrappers
    diagnostics.py                      # Lambda_GC, QQ plot data
    tests/                              # Statistical test implementations
      __init__.py
      base.py                           # Abstract AssociationTest + TestResult
      fisher.py                         # Fisher's exact (refactored from gene_burden.py)
      logistic_burden.py                # Logistic regression burden
      linear_burden.py                  # Linear regression burden
      skat.py                           # SKAT + SKAT-O dispatch
      acat.py                           # ACAT-V, ACAT-O
    weights/
      __init__.py
      base.py                           # Abstract VariantWeighter
      beta_weights.py                   # Beta(MAF) weights
      uniform.py                        # Uniform weights (current behavior)
    backends/                           # SKAT/SKAT-O computation only
      __init__.py
      base.py                           # Abstract SKATBackend
      r_backend.py                      # rpy2 bridge to R SKAT package
      python_backend.py                 # Pure Python SKAT implementation
      davies.py                         # Davies method: ctypes(qfc.c) + Liu fallback

  stages/
    analysis_stages.py                  # MODIFY: add AssociationAnalysisStage
                                        # MODIFY: GeneBurdenAnalysisStage becomes shim
    processing_stages.py                # MODIFY: add PCAComputationStage

  pipeline_core/
    context.py                          # MODIFY: add association_results field
```

**What is NOT changed:**
- `gene_burden.py` — aggregation functions reused by association engine
- `stage_registry.py` — new stages registered in new `_register_association_stages()` function
- `pipeline_core/stage.py` — abstract Stage base unchanged
- `pipeline_core/runner.py` — no changes needed
- All existing stages — unchanged

---

## Integration Points in Detail

### 1. Stage Integration: AssociationAnalysisStage vs GeneBurdenAnalysisStage

**Decision: Wrap, do not replace.**

`GeneBurdenAnalysisStage` continues to exist and run for `--perform-gene-burden` without `--perform-association`. This preserves the existing user experience and avoids breaking changes in the classic pipeline path.

When `--perform-association` is specified (with or without `--perform-gene-burden`), `AssociationAnalysisStage` runs instead. It:
- Runs `fisher` test always (equivalent to gene burden)
- Runs any additional requested tests (logistic burden, SKAT, SKAT-O, ACAT-O)
- Writes `context.gene_burden_results` with the Fisher results (backward compat)
- Writes `context.association_results` with the full multi-test results

```python
class AssociationAnalysisStage(Stage):
    """Multi-test rare variant association analysis."""

    @property
    def name(self) -> str:
        return "association_analysis"

    @property
    def dependencies(self) -> set[str]:
        # Same as GeneBurdenAnalysisStage — needs loaded DataFrame and case/control samples
        return {"dataframe_loading", "sample_config_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        # Run after annotation, and after PCA if present
        return {"custom_annotation", "pca_computation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.config.get("perform_association"):
            return context

        # Chunked processing awareness (same pattern as GeneBurdenAnalysisStage)
        if context.config.get("use_chunked_processing") and not context.config.get(
            "chunked_processing_complete"
        ):
            logger.info("Association analysis deferred — chunked processing not complete")
            return context

        from ..association.engine import AssociationEngine
        from ..association.config import AssociationConfig

        assoc_config = AssociationConfig.from_pipeline_context(context)
        engine = AssociationEngine(assoc_config)

        df = context.current_dataframe
        results = engine.run(df, context)

        # Backward compat: put Fisher results into gene_burden_results
        context.gene_burden_results = results.fisher_results
        # New field for full multi-test results
        context.association_results = results

        # Write outputs
        self._write_outputs(results, context)
        return context
```

**GeneBurdenAnalysisStage becomes a conditional shim:**

```python
class GeneBurdenAnalysisStage(Stage):
    def _process(self, context: PipelineContext) -> PipelineContext:
        # If association is requested, AssociationAnalysisStage handles it
        if context.config.get("perform_association"):
            logger.debug("Skipping GeneBurdenAnalysisStage — AssociationAnalysisStage active")
            return context

        if not context.config.get("perform_gene_burden"):
            return context

        # Original implementation unchanged...
        burden_results = perform_gene_burden_analysis(...)
        context.gene_burden_results = burden_results
        return context
```

**Stage ordering in registry:** Both stages are registered. The dependency graph ensures they don't conflict — `AssociationAnalysisStage` has a soft dependency on `gene_burden_analysis` so if both are present, association runs after (though in practice the shim check prevents double-running).

### 2. Data Flow: Genotype Matrix Extraction

The existing `_aggregate_gene_burden_from_columns()` and `_aggregate_gene_burden_from_gt()` in `gene_burden.py` already solve the genotype matrix extraction problem. The association engine reuses them:

```python
# association/engine.py
from ..gene_burden import (
    _find_gt_columns,
    _gt_to_dosage,
    _aggregate_gene_burden_from_columns,  # For Fisher test only
)

def _extract_genotype_matrix(
    self,
    gene_df: pd.DataFrame,
    vcf_samples: list[str],
    gt_columns: list[str],
) -> np.ndarray:
    """Extract (n_samples, n_variants) dosage matrix for association tests.

    Returns a float64 array with values in {0, 1, 2} per sample per variant.
    Ordered by vcf_samples row order.
    """
    n_samples = len(vcf_samples)
    n_variants = len(gene_df)
    G = np.zeros((n_samples, n_variants), dtype=np.float64)

    for var_idx, (_, row) in enumerate(gene_df.iterrows()):
        for samp_idx, col in enumerate(gt_columns[:n_samples]):
            G[samp_idx, var_idx] = _gt_to_dosage(str(row.get(col, "0/0")))

    return G
```

**Data availability in PipelineContext:**

| Data | Context Field | Available From Stage |
|------|---------------|---------------------|
| Sample IDs (ordered) | `context.vcf_samples` | `sample_config_loading` |
| Phenotype (binary) | `context.config["case_samples"]`, `["control_samples"]` | `phenotype_case_control_assignment` |
| Genotype matrix | Reconstructed from `current_dataframe` GT columns | `dataframe_loading` |
| Covariates | Loaded from file (new) | `association_config_loading` (new) |
| PCA eigenvectors | `context.stage_results["pca_computation"]` | `pca_computation` (new) |
| Pedigree | `context.pedigree_data` | `pedigree_loading` |

**Critical:** The per-sample GT columns (`GEN_0__GT`, `GEN_1__GT`, ...) are the native format since v0.13. The `_find_gt_columns()` function already handles both the sanitized `GEN_N__GT` and original `GEN[N].GT` formats.

### 3. Module Layout and Shared Utilities

**Relationships between packages:**

```
association/engine.py
    uses gene_burden._find_gt_columns()       # locate per-sample GT columns
    uses gene_burden._gt_to_dosage()          # genotype string → dosage
    uses gene_burden._aggregate_gene_burden_from_columns()  # for Fisher test
    uses association/tests/fisher.py          # runs Fisher test (wraps gene_burden logic)
    uses association/tests/logistic_burden.py
    uses association/tests/skat.py
    uses association/correction.py            # FDR/Bonferroni (extracted from gene_burden.py)
    uses association/covariates.py            # covariate matrix assembly
    uses association/weights/                 # variant weights

association/correction.py
    refactors from gene_burden.py             # _apply_multiple_testing_correction()
    gene_burden.py imports from correction.py # backward compat re-export

association/tests/fisher.py
    thin wrapper around gene_burden.perform_gene_burden_analysis()
    ensures identical output to current behavior

inheritance/ package
    no interaction with association/          # inheritance runs before, adds columns
    association reads Inheritance_Pattern column as potential stratification variable

scoring/
    no interaction                            # scoring runs before association
    scored variants available in current_dataframe
```

**What stays in `gene_burden.py` indefinitely:**
- `perform_gene_burden_analysis()` — public API, backward compat
- `_aggregate_gene_burden_from_columns()` — reused by association engine
- `_aggregate_gene_burden_from_gt()` — fallback aggregation
- `_gt_to_dosage()` — utility used everywhere
- `_find_gt_columns()` — utility used everywhere

**What moves to `association/`:**
- `_apply_multiple_testing_correction()` (new shared function, referenced from `gene_burden.py`)
- `_compute_or_confidence_interval()` (moved to `association/results.py`, re-exported from gene_burden)

The `gene_burden.py` module will have a `# TODO: In v0.16, deprecate these in favor of association/` comment but is not removed during this milestone.

### 4. R Backend Isolation

**Never import rpy2 at module level.** The `backends/r_backend.py` module itself is imported lazily:

```python
# association/backends/__init__.py
from .base import SKATBackend

def get_skat_backend(force: str | None = None) -> "SKATBackend":
    """Factory: auto-detect best available SKAT backend.

    Parameters
    ----------
    force : str | None
        "r" to force R backend, "python" to force Python backend, None to auto-detect.
    """
    if force == "python":
        from .python_backend import PythonSKATBackend
        return PythonSKATBackend()

    if force == "r" or force is None:
        try:
            import rpy2  # noqa: F401 — just test importability
            from .r_backend import RSKATBackend
            backend = RSKATBackend()
            # Test R SKAT package availability
            backend.check_skat_available()
            return backend
        except (ImportError, RuntimeError) as e:
            if force == "r":
                raise RuntimeError(
                    f"R SKAT backend requested but not available: {e}"
                ) from e
            logger.info(f"R/rpy2 not available ({e}), using Python SKAT backend")

    from .python_backend import PythonSKATBackend
    return PythonSKATBackend()
```

This factory is called inside `AssociationEngine.__init__()` or `AssociationEngine.run()`, never at import time. The `association/` package can always be imported regardless of R availability.

**The `r_backend.py` file structure:**

```python
# association/backends/r_backend.py
# This module is ONLY imported when rpy2 is available.
# Do not import rpy2 at function level inside this module
# (it's already confirmed available by the caller).

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr

from .base import SKATBackend, NullModel

class RSKATBackend(SKATBackend):
    """Calls R SKAT package via rpy2. Gold standard implementation."""

    def __init__(self):
        self._skat = None  # Lazy R package import

    def check_skat_available(self) -> None:
        """Verify SKAT R package is installed. Raises RuntimeError if not."""
        try:
            self._skat = importr("SKAT")
        except Exception as e:
            raise RuntimeError(f"R SKAT package not available: {e}") from e

    def fit_null_model(self, phenotype, covariates, trait_type="binary") -> NullModel:
        numpy2ri.activate()
        y_r = ro.FloatVector(phenotype)
        if covariates is not None:
            x_r = ro.r.matrix(covariates, nrow=len(phenotype))
            null_model = self._skat.SKAT_Null_Model(y_r, x_r, out_type=trait_type)
        else:
            null_model = self._skat.SKAT_Null_Model(y_r, out_type=trait_type)
        return NullModel(r_obj=null_model)

    def test_gene(self, genotype_matrix, weights, null_model, method="SKAT-O") -> dict:
        numpy2ri.activate()
        G_r = ro.r.matrix(genotype_matrix, nrow=genotype_matrix.shape[0])
        w_r = ro.FloatVector(weights)
        result = self._skat.SKAT(G_r, null_model.r_obj, weights=w_r, method=method)
        return {"p_value": result.rx2("p.value")[0], "statistic": result.rx2("Q")[0]}
```

### 5. PCA as a Pipeline Stage

**New stage: `PCAComputationStage`**

PCA computation needs the multi-sample VCF, which is available as `context.filtered_vcf` (or `context.extracted_vcf`) after the processing stages. It runs before `DataFrameLoadingStage` so the eigenvectors are available when `AssociationAnalysisStage` runs.

```python
class PCAComputationStage(Stage):
    """Compute PCA for population stratification control.

    Supports: loading pre-computed eigenvec file, running AKT, running PLINK 2.0.
    PCA eigenvectors are stored in stage_results and merged as covariates
    by AssociationAnalysisStage.
    """

    @property
    def name(self) -> str:
        return "pca_computation"

    @property
    def dependencies(self) -> set[str]:
        # Needs VCF to exist; runs before DataFrame analysis
        return {"variant_extraction"}  # or "snpsift_filter" for filtered VCF

    @property
    def soft_dependencies(self) -> set[str]:
        return {"bcftools_prefilter", "snpsift_filter"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        pca_file = context.config.get("pca_file")
        pca_tool = context.config.get("pca_tool")  # "akt", "plink", or None

        if not pca_file and not pca_tool:
            logger.debug("PCA not requested")
            return context

        from ..association.pca import load_pca_file, run_akt_pca, run_plink_pca

        if pca_file:
            eigenvec_df = load_pca_file(pca_file, context.vcf_samples)
        elif pca_tool == "akt":
            vcf_path = context.filtered_vcf or context.extracted_vcf
            eigenvec_df = run_akt_pca(vcf_path, context.config)
        elif pca_tool == "plink":
            vcf_path = context.filtered_vcf or context.extracted_vcf
            eigenvec_df = run_plink_pca(vcf_path, context.config)
        else:
            raise ValueError(f"Unknown pca_tool: {pca_tool}")

        # Store eigenvectors keyed by sample ID, accessible to association stage
        context.stage_results["pca_computation"] = eigenvec_df
        logger.info(f"PCA computed: {eigenvec_df.shape[1]-1} components, {len(eigenvec_df)} samples")
        return context
```

**Pipeline position:** `PCAComputationStage` belongs in `processing_stages.py` because it wraps external tools (akt, plink) that consume the VCF. It is registered in the `"processing"` category. The dependency graph ensures it runs after VCF is available but before DataFrame analysis stages need it.

**Why not in setup?** Setup stages load config and metadata. PCA requires actual VCF data from the extraction stages. Processing is the right category.

**VCF availability timeline:**

```
GeneBedCreationStage            ← gene_bed_file available
VariantExtractionStage          ← extracted_vcf available
BCFToolsPrefilterStage          ← (optional)
SnpSiftFilterStage              ← filtered_vcf available   ← PCAComputationStage runs here
FieldExtractionStage            ← extracted_tsv available
...
DataFrameLoadingStage           ← current_dataframe available
...
AssociationAnalysisStage        ← reads stage_results["pca_computation"]
```

### 6. Configuration

**Approach: Extend existing config dict, with optional JSON override.**

The existing pattern uses `context.config` (a `dict[str, Any]`) populated from both `config.json` and CLI args. New association keys follow the same pattern:

```python
# In cli.py — new argument group
assoc_group = parser.add_argument_group("Association Analysis")
assoc_group.add_argument("--perform-association", action="store_true", ...)
assoc_group.add_argument("--association-tests", type=str, default="fisher",
                          help="Comma-separated: fisher,logistic_burden,skat,skat_o,acat_o")
assoc_group.add_argument("--trait-type", choices=["binary", "quantitative"], default="binary")
assoc_group.add_argument("--covariate-file", type=str)
assoc_group.add_argument("--covariates", type=str, help="Comma-separated column names")
assoc_group.add_argument("--pca-file", type=str)
assoc_group.add_argument("--pca-tool", choices=["akt", "plink"])
assoc_group.add_argument("--pca-components", type=int, default=10)
assoc_group.add_argument("--variant-weights", choices=["uniform", "beta", "cadd"], default="uniform")
assoc_group.add_argument("--skat-backend", choices=["auto", "r", "python"], default="auto")
assoc_group.add_argument("--association-config", type=str, help="Path to association JSON config")
assoc_group.add_argument("--association-output", type=str)
```

**JSON config mode** (power user override):

```python
# In AssociationAnalysisStage._process() or a new AssociationConfigLoadingStage:
assoc_config_file = context.config.get("association_config")
if assoc_config_file:
    import json
    with open(assoc_config_file) as f:
        assoc_overrides = json.load(f)
    context.config.update({f"assoc_{k}": v for k, v in assoc_overrides.items()})
```

**No separate JSON file needed for the default case.** The `variantcentrifuge/config.json` does not need association-specific defaults — those live in `AssociationConfig.defaults()`. A separate `association_config.json` is only for power users overriding everything.

**AssociationConfig dataclass:**

```python
# association/config.py
from dataclasses import dataclass, field
from typing import Any

@dataclass
class AssociationConfig:
    tests: list[str] = field(default_factory=lambda: ["fisher"])
    trait_type: str = "binary"
    covariate_file: str | None = None
    covariate_columns: list[str] = field(default_factory=list)
    pca_components: int = 10
    variant_weights: str = "uniform"
    weight_params: dict[str, float] = field(default_factory=dict)
    skat_backend: str = "auto"
    correction_method: str = "fdr"
    output_file: str | None = None
    diagnostics: bool = False

    @classmethod
    def from_pipeline_context(cls, context: "PipelineContext") -> "AssociationConfig":
        cfg = context.config
        tests_str = cfg.get("association_tests", "fisher")
        return cls(
            tests=[t.strip() for t in tests_str.split(",")],
            trait_type=cfg.get("trait_type", "binary"),
            covariate_file=cfg.get("covariate_file"),
            covariate_columns=[c.strip() for c in cfg.get("covariates", "").split(",") if c],
            pca_components=int(cfg.get("pca_components", 10)),
            variant_weights=cfg.get("variant_weights", "uniform"),
            skat_backend=cfg.get("skat_backend", "auto"),
            correction_method=cfg.get("correction_method", "fdr"),
            output_file=cfg.get("association_output"),
            diagnostics=bool(cfg.get("association_diagnostics", False)),
        )
```

### 7. Output Integration

**ExcelReportStage already handles gene burden via `context.config["gene_burden_output"]`.** The same pattern extends to association:

```python
# In ExcelReportStage._process() — MODIFICATION (add after gene burden block):
if context.config.get("perform_association") and hasattr(context, "association_results"):
    assoc_file = context.config.get("association_output")
    if assoc_file and Path(assoc_file).exists():
        try:
            append_tsv_as_sheet(xlsx_file, assoc_file, sheet_name="Association")
            # If diagnostics enabled, add QQ data sheet
            qq_file = context.config.get("qq_data_output")
            if qq_file and Path(qq_file).exists():
                append_tsv_as_sheet(xlsx_file, qq_file, sheet_name="QQ Data")
        except Exception as e:
            logger.error(f"Failed to add Association sheet: {e}")
```

**TSVOutputStage:** Unchanged. The primary variant-level TSV is unchanged; association results are a separate file.

**New output files written by AssociationAnalysisStage:**
- `{base_name}.association.tsv` — combined results, one row per gene, all tests
- `{base_name}.diagnostics/lambda_gc.txt` — inflation factors (if diagnostics=True)
- `{base_name}.diagnostics/qq_data.tsv` — QQ plot data (if diagnostics=True)

**Backward compatibility for `--perform-gene-burden` alone:** Zero changes. `GeneBurdenAnalysisStage` runs, sets `context.gene_burden_results`, writes to `context.config["gene_burden_output"]`, and `ExcelReportStage` picks up the "Gene Burden" sheet exactly as before.

### 8. PipelineContext Extension

Add one field to `context.py`:

```python
@dataclass
class PipelineContext:
    # ... existing fields ...

    # Association analysis results (NEW — v0.15.0)
    association_results: Any | None = None  # AssociationResults dataclass

    # ... in merge_from():
    def merge_from(self, other: "PipelineContext") -> None:
        # ... existing merges ...
        if other.association_results is not None and self.association_results is None:
            self.association_results = other.association_results
```

This is the minimal context change. `association_results` stores the full multi-test results object. `gene_burden_results` continues to hold the Fisher-equivalent results for backward compat.

---

## Component Boundaries

| Component | Responsibility | Communicates With |
|-----------|---------------|-------------------|
| `association/engine.py` | Orchestrate test selection, per-gene execution loop | `tests/`, `weights/`, `covariates.py`, `gene_burden.py` |
| `association/config.py` | Parse and validate association configuration | `pipeline_core/context.py` |
| `association/covariates.py` | Load covariate file, merge PCA, validate sample alignment | `context.stage_results["pca_computation"]` |
| `association/pca.py` | Load eigenvec files or wrap akt/plink | External tools, context.filtered_vcf |
| `association/tests/fisher.py` | Thin wrapper around existing gene_burden logic | `gene_burden.py` |
| `association/tests/logistic_burden.py` | Logistic regression via statsmodels | `statsmodels.api.Logit` |
| `association/tests/skat.py` | SKAT/SKAT-O via backend factory | `backends/get_skat_backend()` |
| `association/backends/r_backend.py` | rpy2 bridge to R SKAT package | R runtime (lazy import) |
| `association/backends/python_backend.py` | Pure Python SKAT | `backends/davies.py` |
| `association/backends/davies.py` | Davies method: ctypes(qfc.c) + Liu fallback | `qfc.so` (compiled C) |
| `AssociationAnalysisStage` | Pipeline integration, context I/O | `association/engine.py`, PipelineContext |
| `PCAComputationStage` | Pipeline integration for PCA | `association/pca.py`, context.filtered_vcf |
| `GeneBurdenAnalysisStage` | Fisher-only backward compat | `gene_burden.py` (unchanged) |
| `ExcelReportStage` | Write Association sheet | context.config["association_output"] |

---

## Data Flow Changes

### Before (v0.14)

```
current_dataframe (GENE, GEN_N__GT, ...)
         |
GeneBurdenAnalysisStage
  _find_gt_columns() → column-based aggregation
  Fisher's exact → FDR correction → gene_burden_results DataFrame
         |
ExcelReportStage → "Gene Burden" sheet
```

### After (v0.15) — with association

```
current_dataframe (GENE, GEN_N__GT, ...)
         |
         |              context.stage_results["pca_computation"]
         |                    |
         ↓                    ↓
AssociationAnalysisStage
  AssociationEngine.run(df, context)
    |
    ├─ covariates.py: merge covariates + PCA eigenvectors → covariate matrix (n × k)
    |
    ├─ For each gene:
    |   _find_gt_columns(df) → gt_cols
    |   _extract_genotype_matrix(gene_df, vcf_samples, gt_cols) → G (n × m)
    |   weights.compute(G, gene_df) → w (m,)
    |   |
    |   ├─ fisher.test_gene(G, phenotype) → TestResult
    |   ├─ logistic_burden.test_gene(G, phenotype, covariates, w) → TestResult
    |   ├─ skat.test_gene(G, phenotype, covariates, w) → TestResult
    |   └─ acat.combine([fisher_p, burden_p, skat_p]) → TestResult
    |
    ├─ correction.apply_fdr([all_p_values])
    |
    └─ results.AssociationResults
         |
         ├─ context.gene_burden_results = results.fisher_as_dataframe()  # backward compat
         ├─ context.association_results = results
         └─ write to {base}.association.tsv

ExcelReportStage
  → "Gene Burden" sheet (from context.config["gene_burden_output"])  ← unchanged
  → "Association" sheet (from context.config["association_output"])  ← NEW
```

---

## Suggested Build Order

The dependency graph dictates this sequence. Each phase is independently testable.

### Phase 1: Foundation — Core Abstractions + Fisher Refactor

**What to build:**
- `association/__init__.py`, `config.py`, `results.py`
- `association/tests/base.py` (AssociationTest ABC, TestResult)
- `association/tests/fisher.py` (wrap `gene_burden.perform_gene_burden_analysis`)
- `association/engine.py` (run single test, produce AssociationResults)
- `AssociationAnalysisStage` (minimal — Fisher only, no covariates)
- `PipelineContext.association_results` field
- `GeneBurdenAnalysisStage` shim check
- CLI args: `--perform-association`, `--association-tests`
- `ExcelReportStage` extension for association sheet
- Register `AssociationAnalysisStage` in `stage_registry.py`

**Why first:** All subsequent phases depend on the stage being in the pipeline and the engine orchestration pattern. Fisher test provides a validatable baseline — output must match current `--perform-gene-burden` exactly.

**Validation gate:** `--perform-association --association-tests fisher` produces bit-identical output to `--perform-gene-burden`.

### Phase 2: Covariate System + Logistic/Linear Burden

**What to build:**
- `association/covariates.py` (load file, validate sample alignment, handle missing)
- `association/weights/` (beta_weights.py, uniform.py, base.py)
- `association/tests/logistic_burden.py` (statsmodels.Logit wrapper)
- `association/tests/linear_burden.py` (statsmodels.OLS wrapper)
- `association/correction.py` (refactor from gene_burden.py)
- CLI args: `--covariate-file`, `--covariates`, `--variant-weights`, `--trait-type`

**Dependencies from Phase 1:** Engine loop, TestResult structure, AssociationAnalysisStage.

**Why second:** Logistic burden uses existing statsmodels (no new dependencies) and provides covariate-adjusted results before the complex SKAT infrastructure is needed. Also establishes the covariate matrix assembly pattern that SKAT reuses.

### Phase 3: R SKAT Backend

**What to build:**
- `association/backends/base.py` (SKATBackend ABC, NullModel)
- `association/backends/__init__.py` (get_skat_backend factory)
- `association/backends/r_backend.py` (RSKATBackend via rpy2)
- `association/tests/skat.py` (dispatch to backend, both SKAT and SKAT-O)
- CLI args: `--skat-backend`
- Integration into AssociationEngine loop

**Dependencies from Phase 2:** AssociationEngine structure, covariate matrix.

**Why third:** R backend is the oracle for validating Python backend. Build and validate R backend first, then use it as reference in Phase 4.

### Phase 4: Pure Python SKAT Backend

**What to build:**
- `association/backends/davies.py` (Liu fallback first, then ctypes qfc.c)
- `association/backends/python_backend.py` (PythonSKATBackend, iterative validation against R)
- `qfc.c` compiled as `qfc.so` (or bundled .pyd for Windows)
- Validation test suite comparing Python vs R p-values

**Dependencies from Phase 3:** R backend output as reference ground truth.

**Why fourth:** Requires R backend for correctness validation. The iterative development loop (implement → compare → fix) needs both backends running simultaneously.

### Phase 5: ACAT-O + Diagnostics

**What to build:**
- `association/tests/acat.py` (ACAT-V and ACAT-O)
- `association/diagnostics.py` (lambda_GC, QQ plot data)
- Diagnostics output in AssociationAnalysisStage
- CLI args: `--diagnostics-output`

**Dependencies from Phase 4:** Requires p-values from all tests; runs as post-processing step.

### Phase 6: PCA Integration

**What to build:**
- `association/pca.py` (file loader for PLINK/AKT/generic formats)
- `PCAComputationStage` (AKT wrapper, PLINK wrapper)
- Register `PCAComputationStage` in `stage_registry.py` under "processing"
- CLI args: `--pca-file`, `--pca-tool`, `--pca-components`
- Integration into covariates.py (merge PCA eigenvectors as covariates)

**Dependencies from Phase 2:** Covariate system must be in place to accept PCA columns.

**Why last:** PCA is optional (pre-computed files can be used from Phase 2 onward via `--pca-file`). The computation stage adds external tool dependencies (akt, plink2) and is architecturally independent of all statistical test phases.

---

## Anti-Patterns to Avoid

### Anti-Pattern 1: Importing rpy2 at Module Level

**What goes wrong:** `import variantcentrifuge` fails on systems without R. Any CI environment without R breaks.

**Prevention:** The `r_backend.py` module is only imported inside the factory function after confirming `rpy2` is importable. The `association/__init__.py` never re-exports from `r_backend`.

### Anti-Pattern 2: Replacing GeneBurdenAnalysisStage

**What goes wrong:** Existing users lose `--perform-gene-burden` behavior or get subtle output differences. Classic pipeline (non-stage-based) also calls gene burden directly.

**Prevention:** `GeneBurdenAnalysisStage` stays. It becomes a conditional shim. `gene_burden.py` is not modified (only extended with imports from `association/correction.py`).

### Anti-Pattern 3: Building Genotype Matrix Before DataFrame is Loaded

**What goes wrong:** `PCAComputationStage` or `AssociationAnalysisStage` depends on `current_dataframe` but runs before `DataFrameLoadingStage`.

**Prevention:** `AssociationAnalysisStage` hard-depends on `"dataframe_loading"`. `PCAComputationStage` hard-depends on `"variant_extraction"` (not DataFrame loading) since it operates on the VCF, not the DataFrame.

### Anti-Pattern 4: Separate Association Config JSON in Default Path

**What goes wrong:** Users must create a second config file for basic association tests; confusing alongside the existing `config.json`.

**Prevention:** All association options available as CLI flags (Phase 1 target). JSON config mode is a power-user override, not the primary interface. The `--association-config` JSON file is optional.

### Anti-Pattern 5: Storing Full Genotype Matrix in PipelineContext

**What goes wrong:** For 5,000 samples × 50,000 variants, the full genotype matrix is 200M float64 values = 1.6 GB. Storing it in context makes memory management impossible.

**Prevention:** The genotype matrix is **never stored in context**. It is extracted per-gene inside the engine loop, used, and discarded. Only the eigenvector DataFrame (n_samples × n_components, typically 5,000 × 20 = tiny) is stored in `context.stage_results`.

### Anti-Pattern 6: Making SKAT Mandatory for Fisher Users

**What goes wrong:** `AssociationAnalysisStage` tries to initialize `PythonSKATBackend` even when only Fisher test is requested, causing import errors or slow startup.

**Prevention:** `get_skat_backend()` is only called when `"skat"` or `"skat_o"` is in the test list. The engine checks `assoc_config.tests` before any backend initialization.

---

## Scalability Considerations

| Concern | Small cohort (<100 samples) | Medium cohort (100-1000 samples) | Large cohort (5000+ samples) |
|---------|---------------------------|----------------------------------|------------------------------|
| Genotype matrix memory | Negligible | ~10 MB per gene batch | Extract per-gene, never full matrix |
| SKAT runtime | <1s per gene | 1-5s per gene | May need parallelization (Phase 6+) |
| PCA computation | Load from file preferred | AKT: <5min | AKT: <5min (optimized for VCF) |
| Covariate alignment | Warn on missing | Warn on missing | Error on >5% missing |
| R rpy2 overhead | Startup cost ~2s amortized | Negligible | Negligible |

---

## Testing Strategy

### Correctness Tests (by phase)

**Phase 1 — Fisher parity:**
```python
# tests/association/test_fisher_parity.py
def test_association_fisher_matches_gene_burden(tmp_path, variant_df, case_samples, control_samples):
    """AssociationAnalysisStage Fisher output must match GeneBurdenAnalysisStage exactly."""
    # Run old stage
    old_results = run_gene_burden_stage(variant_df, case_samples, control_samples)

    # Run new stage with fisher only
    new_results = run_association_stage(variant_df, case_samples, control_samples,
                                        tests=["fisher"])

    pd.testing.assert_frame_equal(
        old_results.sort_values("GENE").reset_index(drop=True),
        new_results.fisher_as_dataframe().sort_values("GENE").reset_index(drop=True),
        check_exact=False, rtol=1e-10
    )
```

**Phase 3/4 — SKAT backend validation:**
```python
# tests/association/test_skat_parity.py
@pytest.mark.skipif(not r_available(), reason="R not available")
def test_python_skat_matches_r_skat(synthetic_cohort):
    """Python SKAT p-values within 10% relative difference of R SKAT."""
    r_backend = RSKATBackend()
    py_backend = PythonSKATBackend()

    for gene, G, pheno, cov in synthetic_cohort.iter_genes():
        r_result = r_backend.test_gene(G, pheno, cov, method="SKAT-O")
        py_result = py_backend.test_gene(G, pheno, cov, method="SKAT-O")

        # Relative difference on log10(p) scale
        rel_diff = abs(np.log10(r_result["p_value"]) - np.log10(py_result["p_value"]))
        assert rel_diff < 0.1 * abs(np.log10(r_result["p_value"])), (
            f"Gene {gene}: R p={r_result['p_value']:.2e}, "
            f"Python p={py_result['p_value']:.2e}"
        )
```

**Synthetic data generation:**
```python
# tests/association/fixtures/synthetic_cohort.py
def make_synthetic_cohort(
    n_cases: int = 100,
    n_controls: int = 100,
    n_genes: int = 50,
    signal_genes: int = 5,  # genes with injected burden signal
    n_variants_per_gene: int = 10,
    seed: int = 42,
) -> SyntheticCohort:
    """Generate synthetic case-control cohort with known signals."""
    rng = np.random.default_rng(seed)
    # ... generate genotype matrix, phenotypes, covariates
```

**Backward compatibility (regression):**
```python
# tests/integration/test_association_backward_compat.py
def test_perform_gene_burden_unchanged(tmp_path, giab_vcf):
    """--perform-gene-burden output unchanged after association framework added."""
    # Run with old flag only
    result = run_pipeline(giab_vcf, flags=["--perform-gene-burden"])
    assert result.gene_burden_tsv.exists()
    assert_gene_burden_matches_golden(result.gene_burden_tsv)
```

---

## Sources

All findings from direct codebase reading:

- `/mnt/c/development/scholl-lab/variantcentrifuge/variantcentrifuge/pipeline_core/context.py` — PipelineContext fields
- `/mnt/c/development/scholl-lab/variantcentrifuge/variantcentrifuge/pipeline_core/stage.py` — Stage ABC, dependency system
- `/mnt/c/development/scholl-lab/variantcentrifuge/variantcentrifuge/stages/analysis_stages.py` — GeneBurdenAnalysisStage (line 1810), InheritanceAnalysisStage, data flow patterns
- `/mnt/c/development/scholl-lab/variantcentrifuge/variantcentrifuge/gene_burden.py` — perform_gene_burden_analysis, _aggregate_gene_burden_from_columns, _gt_to_dosage, _find_gt_columns
- `/mnt/c/development/scholl-lab/variantcentrifuge/variantcentrifuge/stages/stage_registry.py` — registration patterns, category conventions
- `/mnt/c/development/scholl-lab/variantcentrifuge/variantcentrifuge/stages/output_stages.py:791` — ExcelReportStage gene burden integration pattern
- `/mnt/c/development/scholl-lab/variantcentrifuge/.planning/ASSOCIATION-FRAMEWORK-DESIGN.md` — design intent and mathematical specs
