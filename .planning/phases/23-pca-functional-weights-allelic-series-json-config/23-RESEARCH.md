# Phase 23: PCA Integration, Functional Weights, Allelic Series, and JSON Config - Research

**Researched:** 2026-02-21
**Domain:** PCA file parsing, functional variant weights, COAST allelic series test, JSON config
**Confidence:** HIGH

## Summary

Phase 23 extends the association framework with four orthogonal capabilities: PCA integration (for population stratification correction), functional variant weighting (CADD/REVEL), the COAST allelic series test (via the R AllelicSeries package), and JSON config mode (association section in existing config.json). All four capabilities build on the established architecture from Phases 18-22.

PCA support involves parsing three file formats (PLINK .eigenvec, AKT stdout, generic TSV), optionally invoking AKT as a subprocess via a new `PCAComputationStage`, and merging the resulting PC columns into the covariate matrix already managed by `covariates.py`. The existing `load_covariates()` function is file-based and needs a companion function to merge an in-memory numpy array of PCs.

Functional weights extend `weights.py` with two new normalized score functions (CADD: `min(CADD_phred/cap, 1.0)`, REVEL: value as-is) and a combined multiplier scheme. The scores are read from the DataFrame columns `dbNSFP_CADD_phred` and `dbNSFP_REVEL_score` (confirmed from `config.json`). The existing `get_weights()` dispatch function needs new spec strings: `"cadd"`, `"revel"`, `"combined"`.

COAST (Coding-variant Allelic Series Test) is implemented in the R AllelicSeries package (CRAN, v0.1.1.5). It expects an integer annotation vector (1=BMV, 2=DMV, 3=PTV), a genotype matrix, phenotype, and covariates. AllelicSeries is NOT currently installed in the R environment. The COAST test must be wrapped as an `AssociationTest` subclass following the same pattern as `RSKATTest` (rpy2-based, main-thread-only, null model cached, `parallel_safe=False`).

JSON config mode adds an `"association"` section to the existing `config.json` read by `-c/--config`. No new file or flag is introduced for config loading; the config.json loading path already exists. The `AssociationConfig` dataclass fields need 1:1 JSON key mappings. Config validation should use plain dict-key checking with explicit error messages (no `jsonschema` dep needed; codebase currently has no schema validation dependency).

**Primary recommendation:** Implement the four sub-plans in order: 23-01 (PCA), 23-02 (functional weights), 23-03 (COAST), 23-04 (JSON config + QQ plot). Each is independently testable. No sub-plan requires any other sub-plan to be merged first.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| numpy | >=1.24 (in deps) | PCA array operations, weight multiplication | Already used everywhere |
| scipy.stats | >=1.10 (in deps) | beta.pdf for existing weights (no new scipy needed) | Already in deps |
| pandas | >=2.0 (in deps) | DataFrame column lookup for CADD/REVEL/EFFECT/IMPACT | Already used |
| rpy2 | >=3.5.0 (in deps) | R AllelicSeries package calls for COAST | Already used for SKAT |
| subprocess | stdlib | AKT invocation in PCAComputationStage | Already used in processing_stages.py |
| shutil.which | stdlib | AKT PATH detection | Pattern already used in pipeline |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| matplotlib | optional (lazy import) | QQ plot PNG/SVG output | Only if user requests --diagnostics-output and matplotlib is installed |
| json | stdlib | Association config section parsing | Config loading path |
| logging | stdlib | Per-category missing-score warnings | Always |
| pathlib.Path | stdlib | PCA file detection, output paths | Always |
| csv.Sniffer | stdlib | Generic TSV format auto-detection | PCA parser for unknown delimiters |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Plain dict checking for JSON validation | jsonschema library | jsonschema is not in deps; adding it for config validation is over-engineering; explicit error messages are clearer |
| rpy2 for COAST | Pure Python re-implementation | COAST math is non-trivial (ACATV omnibus of burden+SKAT); using R AllelicSeries is authoritative and matches the paper |
| New `--association-config` CLI flag | Association section in existing config.json | Decision locked: add "association" section to existing config.json |

**Installation:** No new Python dependencies. AllelicSeries R package must be installable: `install.packages('AllelicSeries')`. Matplotlib is already in the optional deps universe (check pyproject.toml before adding).

## Architecture Patterns

### Recommended File Structure
```
variantcentrifuge/
├── association/
│   ├── pca.py                    # NEW: PCA file parsers + merge_pca_covariates()
│   ├── weights/                  # NEW: sub-package (rename weights.py -> weights/)
│   │   ├── __init__.py           # Re-export get_weights, beta_maf_weights, uniform_weights
│   │   ├── maf_weights.py        # Existing beta/uniform logic (moved)
│   │   └── functional_weights.py # NEW: cadd_weights, revel_weights, combined_weights
│   ├── tests/
│   │   └── allelic_series.py     # NEW: COASTTest (AssociationTest subclass)
│   └── [existing files unchanged]
└── stages/
    └── processing_stages.py      # EXTEND: add PCAComputationStage at end
```

**Alternative (simpler):** Keep `weights.py` as a single file, add new functions directly. Avoid the sub-package refactor. This is preferred — refactoring to a sub-package risks import churn and test breakage with no benefit.

### Pattern 1: PCA File Parsing (pca.py)

**What:** Auto-detect format from file content, parse to `(sample_ids: list[str], pcs: np.ndarray)`, then `merge_pca_covariates()` handles alignment and column creation for injection into the existing `covariate_matrix`.

**Auto-detection heuristics (verified against PLINK 1.9 and AKT docs):**
- PLINK .eigenvec (no header): first line has FID + IID + numeric columns; 2 sample ID columns before PCs
- PLINK .eigenvec (with header): first line starts with `#FID` or `FID`; columns named PC1, PC2, ...
- AKT stdout format: header-less, first column sample ID, next N columns numeric (single ID column, no FID)
- Generic TSV: first column is sample ID (string), remaining columns are all numeric

**Merge strategy:** PCs become additional covariate columns. The `merge_pca_covariates()` function takes the existing `covariate_matrix` (or None) and horizontally concatenates the PC columns after aligning to VCF sample order. Returns an augmented covariate matrix.

```python
# Source: codebase covariates.py pattern (load_covariates)
def load_pca_file(
    filepath: str,
    vcf_samples: Sequence[str],
    n_components: int = 10,
) -> tuple[np.ndarray, list[str]]:
    """Returns (pc_matrix, column_names) aligned to vcf_samples order."""
    ...

def merge_pca_covariates(
    pca_matrix: np.ndarray,              # (n_samples, n_pcs)
    pca_col_names: list[str],             # ["PC1", "PC2", ...]
    covariate_matrix: np.ndarray | None,  # existing covariates or None
    covariate_col_names: list[str],       # existing names
) -> tuple[np.ndarray, list[str]]:
    """Appends PCA columns to covariate matrix. Returns (merged, names)."""
    ...
```

### Pattern 2: PCAComputationStage (processing_stages.py)

**What:** New Stage subclass that invokes `akt pca` as a subprocess when `--pca-tool akt` is specified. Writes stdout to a temp file, then stores the path as `context.pca_file`. The `DataFrameLoadingStage` → `AssociationAnalysisStage` chain reads this file later.

**AKT invocation:**
```python
# Source: illumina.github.io/akt — akt pca output is space-delimited
# Sample ID in first column, P0 P1 P2 ... columns following
cmd = ["akt", "pca", vcf_file, "-N", str(n_components), "-R", sites_file]
# stdout captured to temp file
```

**Hard error pattern (matches existing ToolNotFoundError):**
```python
if not shutil.which("akt"):
    raise ToolNotFoundError("akt", self.name)
```

**Stage dependencies:** `{"bcftools_prefilter"}` (soft: need a VCF). Add to `_register_processing_stages()` in stage_registry.py.

### Pattern 3: Functional Weights (weights.py extension)

**What:** Two new weight functions added directly to `weights.py`. The `get_weights()` dispatcher is extended to accept new spec strings. Importantly, functional weights require the DataFrame to compute (they read per-variant scores), while the existing `get_weights()` only needs MAF values. This means functional weights have a different calling signature.

**New function signatures:**
```python
def cadd_weights(
    mafs: np.ndarray,
    cadd_scores: np.ndarray,   # from df["dbNSFP_CADD_phred"]
    cap: float = 40.0,
    fallback_lof: float = 1.0,
    fallback_missense: float = 1.0,  # combined = beta(maf) only
    fallback_other: float = 1.0,     # combined = beta(maf) only
    variant_types: np.ndarray | None = None,  # EFFECT column values
) -> np.ndarray:
    """Beta(MAF) * min(CADD_phred / cap, 1.0). Missing scores use type-aware fallback."""
    ...

def revel_weights(
    mafs: np.ndarray,
    revel_scores: np.ndarray,  # from df["dbNSFP_REVEL_score"]
    variant_types: np.ndarray | None = None,
) -> np.ndarray:
    """Beta(MAF) * REVEL_score. REVEL is missense-only; LoF get max weight."""
    ...
```

**get_weights() extension pattern:**
```python
# New spec strings: "cadd", "revel", "combined"
# These require additional kwargs because they need the DataFrame scores
def get_weights(
    mafs: np.ndarray,
    weight_spec: str,
    # NEW optional kwargs for functional weights:
    cadd_scores: np.ndarray | None = None,
    revel_scores: np.ndarray | None = None,
    variant_types: np.ndarray | None = None,
    weight_params: dict | None = None,
) -> np.ndarray:
```

**Column names confirmed from config.json (HIGH confidence):**
- CADD: `dbNSFP_CADD_phred` (in fields_to_extract)
- REVEL: `dbNSFP_REVEL_score` (in fields_to_extract)
- EFFECT: `EFFECT` (aliased from `ANN_0__EFFECT` by `_create_standard_column_aliases`)
- IMPACT: `IMPACT` (aliased from `ANN_0__IMPACT`)

**Type-aware fallback logic (from CONTEXT.md):**
```python
LOF_EFFECTS = frozenset({
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant",
})
LOF_IMPACTS = frozenset({"HIGH"})
MISSENSE_EFFECTS = frozenset({"missense_variant"})

# LoF → functional component = 1.0 (max weight)
# Missense missing score → functional component = 1.0 (no penalty, no boost)
# Other missing score → functional component = 1.0 (Beta(MAF) only)
```

### Pattern 4: COAST Test (allelic_series.py)

**What:** `COASTTest` implements `AssociationTest`. Classification of variants into BMV/DMV/PTV uses the `EFFECT` and `IMPACT` columns plus dbNSFP PolyPhen/SIFT prediction columns from the gene-level DataFrame. The R `AllelicSeries::COAST()` function is called via rpy2.

**Variant classification (from CONTEXT.md, HIGH confidence):**
```python
PTV_EFFECTS = frozenset({
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant",
})
# code=3 (PTV): IMPACT == "HIGH" AND effect in PTV_EFFECTS

# code=2 (DMV): EFFECT == "missense_variant" AND
#   (SIFT == "deleterious" OR PolyPhen in {"probably_damaging", "possibly_damaging"})

# code=1 (BMV): EFFECT == "missense_variant" AND
#   SIFT == "tolerated" AND PolyPhen == "benign"

# Missing SIFT/PolyPhen -> exclude from COAST, include in SKAT/burden
```

**PolyPhen/SIFT column names from dbNSFP:** These come from SnpSift dbNSFP annotations. Typical SnpSift field names are `dbNSFP_SIFT_pred` and `dbNSFP_Polyphen2_HDIV_pred` (or `dbNSFP_Polyphen2_HVAR_pred`). After sanitization `[` → `_` and `.` → `__`, names remain as-is (no brackets in these names). The exact column names depend on user's fields_to_extract config — COAST classification must check for multiple possible column name patterns.

**COAST function signature (verified from rdrr.io/cran/AllelicSeries):**
```r
AllelicSeries::COAST(
  anno,          # integer vector: 1=BMV, 2=DMV, 3=PTV
  geno,          # n_samples x n_variants genotype matrix
  pheno,         # phenotype vector
  covar = NULL,  # covariate matrix (PCs, age, sex)
  weights = c(1, 2, 3),  # default configurable
  is_pheno_binary = FALSE,
  min_mac = 0
)
```

**Return value slots (verified from CRAN vignette):**
- `@Pvals`: DataFrame with burden, SKAT, and omnibus p-values
- `@Counts`: DataFrame with n_bmv, n_dmv, n_ptv counts
- `@Betas`: Effect sizes

**P-value column extraction pattern (rpy2):**
```python
# Source: r_backend.py pattern adapted for AllelicSeries
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

allelic_series = importr("AllelicSeries")
result = allelic_series.COAST(anno=r_anno, geno=r_geno, pheno=r_pheno, ...)
pvals = result.slots["Pvals"]
# Column names in Pvals: "Burden", "SKAT", "Omnibus" (or similar — verify at runtime)
```

**Parallel safety:** `parallel_safe=False` (rpy2 constraint, same as RSKATTest).

**Engine integration:** Add `"coast"` to `_build_registry()` in `engine.py`. The COAST test must also be excluded from ACAT-O computation (it IS the allelic series omnibus — feeding coast into ACAT-O would double-count). Decision needed: coast_p_value feeds into ACAT-O or not. From requirements: `coast_p_value` (omnibus) feeds into ACAT-O. Add it to the ACAT-O combination in `_compute_acat_o()`.

### Pattern 5: JSON Config Mode (association section in config.json)

**What:** The existing `-c/--config` flag already loads `config.json` into `context.config`. Add a new code path in `AssociationAnalysisStage._process()` (or a helper) that reads `context.config.get("association", {})` and populates `AssociationConfig` fields before applying CLI overrides.

**Precedence:** JSON config → CLI overrides. The stage already reads from `context.config`; CLI args are already parsed into `context.config` by cli.py. Since CLI args are set before the stage runs, CLI values are already in `context.config`. The JSON "association" section provides defaults, CLI provides overrides.

**Implementation approach:**
```python
def _build_assoc_config_from_context(context: PipelineContext) -> AssociationConfig:
    """Build AssociationConfig respecting JSON config < CLI override precedence."""
    # Start from JSON "association" section (if present)
    json_assoc = context.config.get("association", {})  # from config.json

    # Validate all JSON keys upfront — fail fast with all errors listed
    _validate_association_config_dict(json_assoc)  # raises ValueError with all bad keys

    # Merge: JSON provides defaults, context.config (CLI) overrides
    # Since CLI values are already parsed into context.config, just read from there
    # JSON only provides values for keys NOT already set by CLI
    return AssociationConfig(
        correction_method=context.config.get("correction_method")
            or json_assoc.get("correction_method", "fdr"),
        ...
    )
```

**Validation approach (no jsonschema dep):**
```python
VALID_ASSOCIATION_KEYS = frozenset({
    "correction_method", "gene_burden_mode", "trait_type",
    "variant_weights", "variant_weight_params", "skat_backend",
    "covariate_file", "covariate_columns", "categorical_covariates",
    "pca_file", "pca_components", "coast_weights",
    "association_tests", "min_cases", "max_case_control_ratio",
    "min_case_carriers", "diagnostics_output",
})

def _validate_association_config_dict(d: dict) -> None:
    unknown = sorted(set(d) - VALID_ASSOCIATION_KEYS)
    if unknown:
        raise ValueError(
            f"Invalid keys in config.json 'association' section: {unknown}. "
            f"Valid keys: {sorted(VALID_ASSOCIATION_KEYS)}"
        )
    # Type checks on known keys...
```

### Pattern 6: Matplotlib QQ Plot (diagnostics.py extension)

**What:** Add `write_qq_plot()` function to `diagnostics.py`. Called from `write_diagnostics()` if matplotlib is importable. Uses lazy import + `matplotlib.use("Agg")` before pyplot import for headless HPC.

```python
def write_qq_plot(
    qq_data: pd.DataFrame,
    output_path: str | Path,
) -> bool:
    """Write QQ plot PNG/SVG. Returns True if successful, False if matplotlib absent."""
    try:
        import matplotlib
        matplotlib.use("Agg")   # MUST be before pyplot import; headless HPC
        import matplotlib.pyplot as plt
    except ImportError:
        logger.debug("matplotlib not available — skipping QQ plot")
        return False

    # Plot logic using qq_data (expected/observed neg log10 p columns)
    ...
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return True
```

**`matplotlib.use("Agg")` timing:** Must be called before `import matplotlib.pyplot`. If pyplot was already imported elsewhere (unlikely in HPC batch mode), use of `Agg` must happen at the first matplotlib import in the process. The lazy-import approach here is correct — call `matplotlib.use()` in the function body before importing pyplot.

### Anti-Patterns to Avoid

- **Do not** store the PCA matrix as a PipelineContext attribute. Context invariant: no large matrices in PipelineContext. Store the file path, read it in AssociationAnalysisStage.
- **Do not** add a new config file for JSON config mode. Decision is locked: use "association" section in existing config.json.
- **Do not** modify `get_weights()` to require a DataFrame parameter. Keep MAF-only weights with the current signature; add functional weights as separate functions with different signatures, called explicitly in AssociationAnalysisStage.
- **Do not** import matplotlib at module level in diagnostics.py (breaks HPC imports).
- **Do not** run COAST for genes where any classification category is absent. The allelic series assumption breaks with empty categories. Log and skip with `None` p-value.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| COAST statistical test | Custom SKAT + burden combination | R AllelicSeries::COAST() | Allelic series math involves ordered alternative hypotheses and ACATV combination; re-implementing is error-prone |
| PCA computation | Custom SVD on VCF | AKT subprocess | AKT implements best practices (RedSVD, LD-thinning recommendations); custom PCA on VCF is a research project |
| CADD normalization | Custom mapping | min(CADD_phred/40, 1.0) | Direct formula from CONTEXT.md; no external lookup |
| JSON schema validation | jsonschema library | Plain dict checking with explicit error messages | No new dependency needed; existing codebase pattern |
| QQ plot | Custom plotting | matplotlib with Agg backend | matplotlib is the standard; graceful degradation pattern already established |

**Key insight:** COAST is the one component that cannot be Python-only without significant research effort. The R AllelicSeries package is the reference implementation.

## Common Pitfalls

### Pitfall 1: PLINK .eigenvec FID vs IID Mismatch
**What goes wrong:** PLINK .eigenvec has two ID columns (FID + IID). Using FID when VCF uses IID (or vice versa) causes 100% sample alignment failure.
**Why it happens:** PLINK distinguishes family ID and individual ID; VCF uses IID only.
**How to avoid:** Always use the second column (IID) when parsing PLINK eigenvec files. Document this explicitly.
**Warning signs:** All samples report as missing from covariate file.

### Pitfall 2: AKT Output Column Order
**What goes wrong:** AKT stdout has no column headers. Column order is sample ID, then P0, P1, ... (zero-indexed). Not PC1, PC2.
**Why it happens:** AKT uses zero-based indexing (P0, P1, ...) in its documentation.
**How to avoid:** Rename to PC1, PC2, ... internally after parsing. Verify this against actual AKT output.
**Warning signs:** Covariate column names in log show "P0" instead of "PC1".

### Pitfall 3: dbNSFP Column Names for PolyPhen/SIFT
**What goes wrong:** PolyPhen/SIFT prediction column names vary by dbNSFP version and which fields the user extracted (e.g., `dbNSFP_Polyphen2_HDIV_pred` vs `dbNSFP_Polyphen2_HVAR_pred` vs `dbNSFP_SIFT_pred`).
**Why it happens:** dbNSFP has multiple PolyPhen predictors; the codebase config.json does not extract PolyPhen columns by default — users must add them.
**How to avoid:** COAST classification should try multiple column name patterns and log clearly which columns it found/used. If no SIFT or PolyPhen columns are present, log a clear error and skip COAST for all genes (not a hard crash).
**Warning signs:** All missense variants end up unclassified (excluded from COAST).

### Pitfall 4: AllelicSeries R Package Not Installed
**What goes wrong:** `COASTTest.check_dependencies()` fails at engine construction time.
**Why it happens:** AllelicSeries is not in the base R installation and is not installed on the test environment (verified: `AllelicSeries NOT installed` in R check).
**How to avoid:** Add install instructions to docs. `check_dependencies()` should raise `ImportError` with exact install command. Tests for COAST must be marked with a custom marker and skipped when AllelicSeries is absent (same pattern as `skat_r` tests).
**Warning signs:** `ImportError` at engine construction with `"--association-tests coast"`.

### Pitfall 5: `matplotlib.use("Agg")` After pyplot Import
**What goes wrong:** `UserWarning: matplotlib is currently using <backend>` or backend switch fails silently.
**Why it happens:** Calling `matplotlib.use()` after pyplot is imported has no effect.
**How to avoid:** Call `matplotlib.use("Agg")` before `import matplotlib.pyplot`. The lazy-import pattern (import inside function body) guarantees this if the function is only called in batch mode.

### Pitfall 6: COAST Not Added to ACAT-O Combination
**What goes wrong:** `coast_p_value` is a primary test result but `_compute_acat_o()` only reads from `self._tests`. ACAT-O would silently omit COAST unless it is registered in the engine.
**Why it happens:** ACAT-O reads registered tests from `self._tests`; COAST must be in the registry.
**How to avoid:** Add `"coast": COASTTest` to `_build_registry()`. The ACAT-O loop will automatically include it.

### Pitfall 7: AssociationConfig Missing New Fields
**What goes wrong:** New config fields (pca_file, pca_components, coast_weights, variant_weight_params) are not in `AssociationConfig` dataclass, causing `AttributeError` in the stage.
**Why it happens:** AssociationConfig must be updated in tandem with new CLI args.
**How to avoid:** Update `AssociationConfig` in `base.py` with new fields and their defaults before implementing the stage logic.

## Code Examples

### PLINK .eigenvec Parsing
```python
# Source: PLINK 1.9 format docs (cog-genomics.org/plink/1.9/formats)
# Format: FID IID PC1 PC2 ... (space-delimited, no header by default)
def _parse_plink_eigenvec(lines: list[str]) -> tuple[list[str], np.ndarray]:
    """Parse PLINK .eigenvec. Uses IID (column 1) as sample ID."""
    sample_ids, rows = [], []
    for line in lines:
        parts = line.split()
        if not parts:
            continue
        # FID=parts[0], IID=parts[1], PCs=parts[2:]
        sample_ids.append(parts[1])
        rows.append([float(x) for x in parts[2:]])
    return sample_ids, np.array(rows, dtype=float)
```

### CADD Weight Computation
```python
# Source: CONTEXT.md — CADD normalization: min(CADD_phred / 40, 1.0)
def cadd_weights(
    mafs: np.ndarray,
    cadd_scores: np.ndarray,
    variant_types: np.ndarray | None = None,
    cap: float = 40.0,
) -> np.ndarray:
    maf_w = beta_maf_weights(mafs)
    nan_mask = np.isnan(cadd_scores)
    functional = np.where(nan_mask, 1.0, np.minimum(cadd_scores / cap, 1.0))
    # LoF with missing score -> functional=1.0 (already default), log warning count
    return maf_w * functional
```

### COAST Classification
```python
# Source: CONTEXT.md classification logic
PTV_IMPACTS = frozenset({"HIGH"})
PTV_EFFECTS = frozenset({
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant",
})
MISSENSE_EFFECT = "missense_variant"
SIFT_DAMAGING = frozenset({"deleterious"})
POLYPHEN_DAMAGING = frozenset({"probably_damaging", "possibly_damaging"})
SIFT_BENIGN = frozenset({"tolerated"})
POLYPHEN_BENIGN = frozenset({"benign"})

def classify_variants(
    effect: pd.Series,
    impact: pd.Series,
    sift_pred: pd.Series | None,
    polyphen_pred: pd.Series | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (annotation_codes, include_mask).
    annotation_codes: int array, 1=BMV, 2=DMV, 3=PTV, 0=unclassified
    include_mask: bool array, True if variant should be included in COAST
    """
    ...
```

### Matplotlib QQ Plot (lazy import)
```python
# Source: matplotlib.org/stable/users/explain/figure/backends.html
def write_qq_plot(qq_data: pd.DataFrame, output_path: str | Path) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")       # headless HPC — MUST be before pyplot
        import matplotlib.pyplot as plt
    except ImportError:
        logger.debug("matplotlib not available — QQ plot skipped")
        return False

    fig, ax = plt.subplots(figsize=(6, 6))
    for test_name, group in qq_data.groupby("test"):
        ax.scatter(group["expected_neg_log10_p"], group["observed_neg_log10_p"],
                   s=4, alpha=0.6, label=test_name)
    max_val = max(qq_data["expected_neg_log10_p"].max(),
                  qq_data["observed_neg_log10_p"].max()) + 0.5
    ax.plot([0, max_val], [0, max_val], "k--", linewidth=0.8)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True
```

### rpy2 AllelicSeries Call Pattern
```python
# Source: r_backend.py pattern + AllelicSeries CRAN vignette
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr

numpy2ri.activate()
allelic_series_pkg = importr("AllelicSeries")

result = allelic_series_pkg.COAST(
    anno=ro.IntVector(anno_codes),
    geno=ro.r.matrix(geno.T.ravel(), nrow=n_samples, ncol=n_variants),
    pheno=ro.FloatVector(phenotype),
    covar=covar_r if covar is not None else ro.NULL,
    weights=ro.FloatVector([1.0, 2.0, 3.0]),
    is_pheno_binary=ro.BoolVector([True]),
)
pvals = result.slots["Pvals"]
# Column access: pvals.rx2("Omnibus")[0], pvals.rx2("Burden")[0], pvals.rx2("SKAT")[0]
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Population structure control via family-based studies only | PCA-based stratification control (10-20 PCs) | GWAS era ~2006; standard since | Required for unrelated case-control studies |
| Uniform variant weighting in burden tests | Beta(MAF) + functional score combined weights | SKAT paper 2011; CADD/REVEL integration 2014+ | Increases power for pathogenic variants |
| Single-test burden/SKAT for gene discovery | Allelic series tests (COAST) | McCaw et al. AJHG 2023 | Better power when effect magnitude increases with damage severity |
| CLI-only configuration | JSON config file support | Common pipeline pattern | Reproducibility for batch HPC workflows |

**Deprecated/outdated:**
- Using FID as sample ID for PLINK files: always use IID for VCF matching
- CADD raw scores for weighting: Phred normalization (divide by 40) is the standard for multiplicative weights

## Open Questions

1. **COAST Pvals column names**
   - What we know: AllelicSeries COAST returns a `@Pvals` slot with burden, SKAT, and omnibus components
   - What's unclear: Exact R column names in the `@Pvals` data frame (rdrr.io documentation was not specific; the CRAN vignette showed "Burden", "SKAT", "Omnibus" but these need runtime verification)
   - Recommendation: In `COASTTest.check_dependencies()`, add a smoke-test that calls COAST on a tiny synthetic dataset and inspects the Pvals column names; log them at DEBUG level. This verifies the exact column names at runtime.

2. **PolyPhen/SIFT column names in user data**
   - What we know: dbNSFP provides `dbNSFP_SIFT_pred` and `dbNSFP_Polyphen2_HDIV_pred` (among others); default config.json does NOT extract these
   - What's unclear: Whether users commonly extract these; what the exact column names are after SnpSift extraction and sanitization
   - Recommendation: Support a configurable list of PolyPhen/SIFT column names via `--variant-weight-params` or hardcode a list of common patterns to try (in priority order). Log which were found.

3. **n_components slice in load_pca_file**
   - What we know: Default is 10 PCs; warn at >20
   - What's unclear: Whether `n_components` should slice the file (take first N PCs) or validate that the file has exactly N PCs
   - Recommendation: Slice to first `n_components` columns. If file has fewer PCs than requested, log a warning and use what's available (do not error).

4. **COAST + ACAT-O interaction**
   - What we know: `coast_p_value` feeds into ACAT-O (from CONTEXT.md)
   - What's unclear: Whether the COAST omnibus p-value should be the one fed to ACAT-O, or whether COAST's burden and SKAT components should be fed separately
   - Recommendation: Feed `coast_p_value` (omnibus, `@Pvals$Omnibus`) to ACAT-O as a single entry — this avoids double-counting COAST's internal combination.

## Sources

### Primary (HIGH confidence)
- CRAN AllelicSeries vignette (rdrr.io/cran/AllelicSeries) — COAST function signature, return value structure
- PLINK 1.9 format reference (cog-genomics.org/plink/1.9/formats) — .eigenvec file format with FID+IID columns
- AKT documentation (illumina.github.io/akt) — PCA output format (space-delimited, SAMPLE_ID P0 P1 ...)
- variantcentrifuge/config.json — column names: dbNSFP_CADD_phred, dbNSFP_REVEL_score
- variantcentrifuge/association/base.py — AssociationConfig dataclass fields
- variantcentrifuge/association/covariates.py — load_covariates() pattern for merge_pca_covariates()
- variantcentrifuge/stages/analysis_stages.py — _STANDARD_COLUMN_ALIASES (EFFECT, IMPACT column names), AssociationAnalysisStage._process() pattern
- variantcentrifuge/association/backends/r_backend.py — rpy2 call pattern, detect_environment() for AllelicSeries check_dependencies() template

### Secondary (MEDIUM confidence)
- McCaw et al. AJHG 2023 (pmc.ncbi.nlm.nih.gov/articles/PMC10432147) — COAST allelic series test methodology, default weights w=(1,2,3)
- matplotlib backends documentation (matplotlib.org/stable/users/explain/figure/backends.html) — Agg backend for headless HPC, use() timing requirement
- AllelicSeries CRAN package description (insitro.r-universe.dev/AllelicSeries) — package exists, version 0.1.1.5

### Tertiary (LOW confidence)
- Exact column names in AllelicSeries COAST `@Pvals` slot — inferred from vignette text, needs runtime verification

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all dependencies already in codebase or standard stdlib
- Architecture (PCA, weights, JSON config): HIGH — follows established codebase patterns exactly
- Architecture (COAST): MEDIUM — rpy2 call pattern is HIGH confidence; exact Pvals slot column names are LOW (needs runtime verification)
- Column names (CADD, REVEL, EFFECT, IMPACT): HIGH — confirmed from config.json and analysis_stages.py
- Pitfalls: HIGH — derived from existing codebase patterns and verified external tool docs

**Research date:** 2026-02-21
**Valid until:** 2026-03-23 (30 days — stable domain; AllelicSeries is mature, PLINK/AKT formats are frozen)
