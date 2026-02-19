# Phase 18: Foundation — Core Abstractions and Fisher Refactor - Research

**Researched:** 2026-02-19
**Domain:** Python ABC pattern, association test framework, pipeline stage integration, Excel report extension
**Confidence:** HIGH (all findings from direct codebase inspection)

## Summary

Phase 18 creates the `association/` package skeleton and refactors the Fisher's exact test
from `gene_burden.py` into the new framework. The codebase has a well-established Stage
ABC pattern (`pipeline_core/stage.py`) and a clear mechanism for extending `PipelineContext`
with new result fields. The existing `GeneBurdenAnalysisStage` (in `stages/analysis_stages.py`)
provides the exact dependency graph pattern and pipeline integration pattern to follow.

The key risk area is guaranteeing bit-identical output between the refactored Fisher path
and the existing `gene_burden.py` Fisher path. The existing code uses a three-strategy
aggregation system — the correct approach is to reuse all three aggregation strategies
unchanged in `gene_burden.py` and feed their output through the new `FisherExactTest`
ABC implementation.

**Primary recommendation:** Model `AssociationAnalysisStage` exactly after `GeneBurdenAnalysisStage`
— same dependency declaration, same config-flag guard in `_process`, same output-to-context
pattern. Add `association_results: pd.DataFrame | None = None` to `PipelineContext` as a
plain dataclass field (no other pattern exists for result fields). Add "Association" sheet
in `ExcelReportStage._add_additional_sheets` using the existing `append_tsv_as_sheet` call
pattern established for "Gene Burden".

## Standard Stack

### Core (all already present in the project)
| Library | Version | Purpose | Source |
|---------|---------|---------|--------|
| `scipy.stats.fisher_exact` | project dep | Fisher 2x2 test | `gene_burden.py` line 32 |
| `statsmodels.stats.multitest` | project dep | FDR/Bonferroni correction | `gene_burden.py` line 39 |
| `statsmodels.stats.contingency_tables.Table2x2` | project dep | OR confidence intervals | `gene_burden.py` line 40 |
| `pandas` | project dep | Gene-level result DataFrames | used throughout |
| `numpy` | project dep | NaN handling, array ops | used throughout |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `openpyxl` | project dep | Appending Excel sheets | `converter.append_tsv_as_sheet` |
| `abc.ABC`, `abc.abstractmethod` | stdlib | AssociationTest base class | New `association/` package |

### Alternatives Considered
None — the correct libraries are already used in `gene_burden.py`. The refactor uses
the identical imports to guarantee identical numerics.

**No new dependencies needed for Phase 18.**

## Architecture Patterns

### Recommended Project Structure
```
variantcentrifuge/
├── association/                   # New package
│   ├── __init__.py                # Exports: AssociationEngine, AssociationTest, TestResult, AssociationConfig
│   ├── base.py                    # AssociationTest ABC + TestResult dataclass + AssociationConfig
│   ├── engine.py                  # AssociationEngine: orchestrates test dispatch per gene
│   ├── correction.py              # apply_correction(pvals, method) — FDR + Bonferroni
│   └── tests/
│       ├── __init__.py            # Exports: FisherExactTest
│       └── fisher.py              # FisherExactTest(AssociationTest) — clean reimplementation
```

### Pattern 1: AssociationTest ABC
**What:** Abstract base class with a single required method; TestResult is a dataclass.
**When to use:** Every association test (Fisher, Burden, SKAT) inherits from this.

```python
# Source: modelled on pipeline_core/stage.py ABC pattern
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

@dataclass
class TestResult:
    gene: str
    test_name: str                   # "fisher", "burden", "skat"
    p_value: float | None            # None = not run / zero-variant skip
    corrected_p_value: float | None
    effect_size: float | None        # OR for Fisher/burden; beta/tau for regression tests
    ci_lower: float | None
    ci_upper: float | None
    n_cases: int
    n_controls: int
    n_variants: int
    extra: dict = field(default_factory=dict)  # test-specific fields

class AssociationTest(ABC):
    @property
    @abstractmethod
    def name(self) -> str:
        """Unique short identifier, e.g. 'fisher'."""

    @abstractmethod
    def run(
        self,
        gene: str,
        contingency_data: dict,           # pre-aggregated counts from gene_burden.py aggregators
        config: "AssociationConfig",
    ) -> TestResult:
        """Run test on one gene. Return TestResult with p_value=None to skip."""
```

**Key design decision from CONTEXT.md:** `FisherExactTest.run()` receives
pre-aggregated contingency data (the same dicts that `_aggregate_gene_burden_from_columns`
and `_aggregate_gene_burden_from_gt` already produce). This keeps the aggregation
logic in `gene_burden.py` untouched for Phase 18.

### Pattern 2: AssociationEngine
**What:** Orchestrates test registration, per-gene dispatch, and result collection.
**When to use:** Called from `AssociationAnalysisStage._process`.

```python
# Source: modelled on existing pipeline orchestration patterns
class AssociationEngine:
    def __init__(self, tests: list[AssociationTest], config: AssociationConfig):
        self._tests = {t.name: t for t in tests}
        self._config = config

    @classmethod
    def from_names(cls, test_names: list[str], config: AssociationConfig) -> "AssociationEngine":
        """Resolve names to registered test instances. Hard error on unknown name."""
        available = {"fisher": FisherExactTest}   # grows with phases
        tests = []
        for name in test_names:
            if name not in available:
                available_str = ", ".join(sorted(available))
                raise ValueError(
                    f"Test '{name}' is not available. Available tests: {available_str}"
                )
            tests.append(available[name]())
        return cls(tests, config)

    def run_all(self, gene_burden_data: list[dict]) -> pd.DataFrame:
        """Run all registered tests per gene. Returns wide-format DataFrame."""
```

### Pattern 3: AssociationAnalysisStage
**What:** Stage that calls AssociationEngine and stores results in PipelineContext.
**When to use:** Registered in stage_registry.py; added to pipeline.py when `--perform-association` active.

```python
# Source: mirrors GeneBurdenAnalysisStage in stages/analysis_stages.py
class AssociationAnalysisStage(Stage):
    @property
    def name(self) -> str:
        return "association_analysis"

    @property
    def dependencies(self) -> set[str]:
        return {"dataframe_loading", "sample_config_loading"}  # same as GeneBurdenAnalysisStage

    @property
    def soft_dependencies(self) -> set[str]:
        return {"custom_annotation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.config.get("perform_association"):
            logger.debug("Association analysis not requested")
            return context
        # ... validate eager dependencies (scipy, statsmodels) ...
        # ... call AssociationEngine.run_all(gene_burden_data) ...
        context.association_results = results_df
        return context
```

### Pattern 4: Extending PipelineContext
**What:** Add `association_results` field as a plain dataclass field.
**When to use:** Any new named result that stages produce for downstream consumers.

```python
# Source: context.py — follow existing pattern at line 160
# In PipelineContext dataclass, after gene_burden_results:
association_results: pd.DataFrame | None = None

# Also update merge_from() to handle the new field:
if other.association_results is not None and self.association_results is None:
    self.association_results = other.association_results
```

### Pattern 5: Adding an Excel Sheet
**What:** Extend `ExcelReportStage._add_additional_sheets` with "Association" sheet.
**When to use:** When `context.config.get("perform_association")` is True.

```python
# Source: output_stages.py lines 790-804 (Gene Burden sheet pattern)
if context.config.get("perform_association"):
    assoc_file = context.config.get("association_output")
    if (
        assoc_file
        and Path(assoc_file).exists()
        and Path(assoc_file).stat().st_size > 0
    ):
        try:
            append_tsv_as_sheet(xlsx_file, assoc_file, sheet_name="Association")
        except Exception as e:
            logger.error(f"Failed to add Association sheet: {e}")
```

### Pattern 6: Registering a New Stage
**What:** Add new stage to both `stage_registry.py` and `pipeline.py`.
**When to use:** For `AssociationAnalysisStage`.

```python
# In stage_registry.py _register_analysis_stages():
from .analysis_stages import AssociationAnalysisStage
register_stage(AssociationAnalysisStage, "analysis", ["association_analysis", "association"], 30.0)

# In pipeline.py build_stage_list (after GeneBurdenAnalysisStage block):
if hasattr(args, "perform_association") and args.perform_association:
    stages.append(AssociationAnalysisStage())
```

### Pattern 7: Soft Dependency in ExcelReportStage
**What:** ExcelReportStage already has `soft_dependencies = {"gene_burden_analysis"}`.
Must also add `"association_analysis"` so the sheet is present.

```python
# In ExcelReportStage.soft_dependencies (output_stages.py line 676-678):
@property
def soft_dependencies(self) -> set[str]:
    return {"gene_burden_analysis", "association_analysis"}  # add new stage
```

### Pattern 8: correction.py shared module
**What:** Extract `smm.multipletests` calls into `association/correction.py`; re-export from `gene_burden.py`.

```python
# association/correction.py
def apply_correction(pvals: list[float], method: str) -> list[float]:
    """Apply FDR (fdr_bh) or Bonferroni correction. Returns corrected p-values."""
    if smm is None:
        return pvals
    if method == "bonferroni":
        return list(smm.multipletests(pvals, method="bonferroni")[1])
    return list(smm.multipletests(pvals, method="fdr_bh")[1])

# gene_burden.py — replace inline multipletests with:
from .association.correction import apply_correction
corrected_pvals = apply_correction(list(pvals), correction_method)
```

### Anti-Patterns to Avoid
- **Delegation from association to gene_burden:** Importing `gene_burden.perform_gene_burden_analysis`
  inside `association/tests/fisher.py` inverts the dependency direction and blocks Phase 22 deprecation.
- **Storing genotype matrix in PipelineContext:** Architecture invariant. Pass pre-aggregated
  contingency dicts into `AssociationTest.run()`, not raw DataFrames.
- **Importing association package from gene_burden.py (except correction.py):**
  The only coupling allowed is `gene_burden.py` importing `association.correction`.
- **Discovering CLI args as `getattr(args, ...)` inside a Stage:** All CLI args flow through
  `context.config` (see line 1011 in cli.py: `cfg["perform_gene_burden"] = args.perform_gene_burden`).
  Do the same for `perform_association` and `association_tests`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Fisher 2x2 test | custom p-value formula | `scipy.stats.fisher_exact` | numerical precision, handles edge cases |
| FDR correction | manual BH step-up | `statsmodels.stats.multitest.multipletests(method="fdr_bh")` | already in project, exact same numerics |
| OR confidence intervals | manual logit CI | `statsmodels.stats.contingency_tables.Table2x2.oddsratio_confint` | robust fallback chain already implemented |
| Excel sheet append | xlsxwriter direct write | `converter.append_tsv_as_sheet` | handles openpyxl mode="a", tested |
| Stage dependency graph | custom ordering | Stage ABC `dependencies` + `soft_dependencies` | PipelineRunner already does topo-sort |
| Test registry | custom dict | `AssociationEngine.from_names` dict lookup | straightforward, phase-by-phase additions |

## Common Pitfalls

### Pitfall 1: Correction Module Circular Import
**What goes wrong:** `association/correction.py` imports from `statsmodels`; `gene_burden.py`
imports from `association.correction`; if `association/__init__.py` also imports from
`gene_burden`, you get a circular import.
**Why it happens:** Package initialization order.
**How to avoid:** `association/__init__.py` must not import anything from `gene_burden.py`.
Only `gene_burden.py` imports from `association.correction`. Keep `correction.py` leaf-level
(only imports stdlib + statsmodels).
**Warning signs:** `ImportError: cannot import name 'X' from partially initialized module`.

### Pitfall 2: Bit-Identity Broken by Floating Point Evaluation Order
**What goes wrong:** `fisher_exact(table)` returns `(odds_ratio, pval)` — the tuple order
is `(statistic, pvalue)`. The current `gene_burden.py` at line 468 does:
`odds_ratio, pval = fisher_exact(table)`.
**Why it happens:** Reversing the tuple assignment gives wrong values silently.
**How to avoid:** Use exact same unpacking order in `FisherExactTest.run()`. The test suite
(Plan 18-04) will catch this, but get it right from the start.

### Pitfall 3: Gene Sort Order Breaks Correction
**What goes wrong:** Multiple testing correction depends on the order of p-values passed
to `multipletests`. If `AssociationEngine` iterates genes in a different order than
`gene_burden.py`, corrected p-values differ even when raw p-values match.
**Why it happens:** `gene_burden.py` sorts by GENE at line 411:
`grouped = grouped.sort_values("GENE").reset_index(drop=True)`.
**How to avoid:** `AssociationEngine.run_all` must sort gene_burden_data by gene name
before iterating, exactly mirroring this sort. The parity test in Plan 18-04 will
catch any mismatch.

### Pitfall 4: AssociationAnalysisStage Runs Before Case/Control Resolved
**What goes wrong:** `context.config.get("case_samples")` returns empty list if
`SampleConfigLoadingStage` hasn't resolved file-based sample lists yet.
**Why it happens:** Same issue exists in `GeneBurdenAnalysisStage` (lines 1882-1904).
**How to avoid:** Declare `"sample_config_loading"` as a hard dependency (same as
`GeneBurdenAnalysisStage`). Also add the same warning log for empty case/control
and return early without error — the exact same guard pattern.

### Pitfall 5: `--association-tests` Without `--perform-association`
**What goes wrong:** If the CLI accepts `--association-tests` but the validation is in
the Stage's `_process` rather than in `cli.py`, the user gets a confusing silent no-op.
**Why it happens:** Stage guards check `context.config.get("perform_association")`.
**How to avoid:** Add validation in `cli.py` after arg parsing: if `args.association_tests`
is set and `args.perform_association` is False, call `parser.error(...)` immediately.

### Pitfall 6: Missing scipy/statsmodels Causes AttributeError Deep in Pipeline
**What goes wrong:** `fisher_exact = None` when scipy not installed; calling `fisher_exact(table)`
raises `TypeError: 'NoneType' object is not callable` at gene-iteration time.
**Why it happens:** `gene_burden.py` has the same try/except pattern at import time.
**How to avoid:** Eager validation at `AssociationAnalysisStage._process` entry:
```python
from .association.tests.fisher import _check_dependencies
_check_dependencies()  # raises ImportError with clear message if missing
```
This matches the CONTEXT.md requirement: "Missing dependencies checked eagerly at startup."

### Pitfall 7: Zero-Variant Genes Not Skipped Silently
**What goes wrong:** Gene with zero qualifying variants produces p_value=NaN from
fisher_exact([[0,0],[n,m]]) and logs confusing warnings.
**Why it happens:** `gene_burden.py` skips genes where `p_count == 0 and c_count == 0`
(line 425) but doesn't explicitly handle zero-variant case.
**How to avoid:** In `FisherExactTest.run()`, check `contingency_data["n_qualifying_variants"] == 0`
and return `TestResult(..., p_value=None, ...)`. The CONTEXT.md specifies: "skip silently,
report p_value=NA."

## Code Examples

### Existing Fisher Core to Port (from gene_burden.py)

```python
# Source: variantcentrifuge/gene_burden.py lines 467-475
# This is the ~80 lines of statistical core to port cleanly

if fisher_exact is not None:
    odds_ratio, pval = fisher_exact(table)
else:
    odds_ratio = float("nan")
    pval = 1.0

ci_lower, ci_upper = _compute_or_confidence_interval(
    table, odds_ratio, ci_method, ci_alpha, continuity_correction
)
```

### Existing Correction Call to Extract

```python
# Source: variantcentrifuge/gene_burden.py lines 507-510
# This block becomes association/correction.py
if correction_method == "bonferroni":
    corrected_pvals = smm.multipletests(pvals, method="bonferroni")[1]
else:
    corrected_pvals = smm.multipletests(pvals, method="fdr_bh")[1]
```

### Existing GeneBurdenAnalysisStage dependency pattern (exact model to follow)

```python
# Source: variantcentrifuge/stages/analysis_stages.py lines 1823-1833
@property
def dependencies(self) -> set[str]:
    return {"dataframe_loading", "sample_config_loading"}

@property
def soft_dependencies(self) -> set[str]:
    return {"custom_annotation"}
```

### Existing CLI arg mapping pattern

```python
# Source: variantcentrifuge/cli.py lines 380-396, 1011-1014
# Statistical Analysis group in create_parser():
stats_group.add_argument("--perform-gene-burden", action="store_true", ...)
stats_group.add_argument("--gene-burden-mode", choices=["samples", "alleles"], ...)
stats_group.add_argument("--correction-method", choices=["fdr", "bonferroni"], ...)

# In main() after arg parsing:
cfg["perform_gene_burden"] = args.perform_gene_burden
cfg["gene_burden_mode"] = args.gene_burden_mode
cfg["correction_method"] = args.correction_method

# Follow same pattern for:
# cfg["perform_association"] = args.perform_association
# cfg["association_tests"] = [t.strip() for t in (args.association_tests or "fisher").split(",")]
# cfg["skat_backend"] = args.skat_backend
```

### Existing Excel sheet append pattern

```python
# Source: variantcentrifuge/stages/output_stages.py lines 791-804
# Gene Burden sheet — Association sheet follows identical pattern:
if context.config.get("perform_gene_burden"):
    gene_burden_file = context.config.get("gene_burden_output")
    if (
        gene_burden_file
        and Path(gene_burden_file).exists()
        and Path(gene_burden_file).stat().st_size > 0
    ):
        try:
            append_tsv_as_sheet(xlsx_file, gene_burden_file, sheet_name="Gene Burden")
        except Exception as e:
            logger.error(f"Failed to add Gene Burden sheet: {e}")
```

### Existing PipelineContext result field pattern

```python
# Source: variantcentrifuge/pipeline_core/context.py lines 160-161
# Analysis results — follow exactly for association_results:
gene_burden_results: pd.DataFrame | None = None

# And in merge_from() at line 300-301:
if other.gene_burden_results is not None and self.gene_burden_results is None:
    self.gene_burden_results = other.gene_burden_results
```

## State of the Art

| Old Approach | Current Approach | Impact for Phase 18 |
|--------------|------------------|---------------------|
| Monolithic `pipeline.py` | Stage-based `pipeline_core/` | AssociationAnalysisStage follows Stage ABC |
| Direct Fisher in gene_burden.py | Refactored to FisherExactTest ABC | Fisher code extracted, interface standardized |
| FDR inline in gene_burden.py | Shared `association/correction.py` | gene_burden.py re-imports from there |
| Single output sheet (Gene Burden) | Multiple sheets per analysis mode | "Association" sheet added alongside |

## Open Questions

1. **Wide-format Association sheet column schema**
   - What we know: One row per gene; fisher_p, fisher_corrected_p, fisher_or, fisher_or_ci_lower,
     fisher_or_ci_upper, n_cases, n_controls, n_variants columns specified in context.
   - What's unclear: Exact column naming convention for test-prefixed columns.
   - Recommendation: Define column names in `TestResult` extra dict or in `AssociationEngine`
     output columns. Establish the naming pattern in Plan 18-01 (the ABC/TestResult definition)
     so it flows consistently to the Excel sheet.

2. **Association output file path convention**
   - What we know: Gene burden writes to `{base_name}.gene_burden.tsv` (line 1995 of
     analysis_stages.py). The context stores path in `context.config["gene_burden_output"]`.
   - What's unclear: Whether association output should be `{base_name}.association.tsv` (one file
     with wide format) or one file per test.
   - Recommendation: Single wide-format TSV `{base_name}.association.tsv` following the same
     convention, stored in `context.config["association_output"]`.

3. **AssociationConfig vs reusing gene_burden config keys**
   - What we know: `GeneBurdenAnalysisStage` reads `context.config.get("gene_burden_mode")`,
     `context.config.get("correction_method")`, etc.
   - What's unclear: Should `AssociationConfig` read the same `correction_method` key
     (shared with gene burden) or use a new `association_correction_method` key?
   - Recommendation: Reuse `correction_method` for Phase 18 since Fisher uses the same
     correction. Phase 22+ can add per-test correction config if needed.

## Sources

### Primary (HIGH confidence)
All findings derived from direct codebase inspection:

- `variantcentrifuge/gene_burden.py` — full Fisher implementation, aggregation strategies,
  FDR/Bonferroni correction, CI calculation
- `variantcentrifuge/pipeline_core/stage.py` — Stage ABC with `name`, `dependencies`,
  `soft_dependencies`, `_process` abstractmethod
- `variantcentrifuge/pipeline_core/context.py` — PipelineContext dataclass field patterns,
  `merge_from` method, `gene_burden_results` field
- `variantcentrifuge/stages/analysis_stages.py` — `GeneBurdenAnalysisStage` full implementation
  (dependency declarations, config guards, result storage)
- `variantcentrifuge/stages/output_stages.py` — `ExcelReportStage._add_additional_sheets`,
  Gene Burden sheet append pattern, `soft_dependencies`
- `variantcentrifuge/stages/stage_registry.py` — `register_stage`, `_register_analysis_stages`,
  `initialize_registry` pattern
- `variantcentrifuge/pipeline.py` — stage list construction, `GeneBurdenAnalysisStage` gating
- `variantcentrifuge/cli.py` — `create_parser`, arg-to-config mapping for gene burden args
- `variantcentrifuge/converter.py` — `append_tsv_as_sheet` signature and openpyxl mode

## Metadata

**Confidence breakdown:**
- Package structure: HIGH — follows existing pipeline_core ABC pattern exactly
- AssociationTest ABC design: HIGH — derived from Stage ABC; Claude's discretion per CONTEXT.md
- Fisher refactor correctness: HIGH — gene_burden.py is the reference; parity test guards
- Pipeline integration: HIGH — mirrors GeneBurdenAnalysisStage line-for-line patterns
- Excel sheet addition: HIGH — identical to Gene Burden sheet append
- Correction module extraction: HIGH — straightforward refactor of two smm calls
- CLI args: HIGH — mirrors existing --perform-gene-burden group exactly

**Research date:** 2026-02-19
**Valid until:** Stable — this is internal codebase knowledge, not external library APIs.
No expiry.
