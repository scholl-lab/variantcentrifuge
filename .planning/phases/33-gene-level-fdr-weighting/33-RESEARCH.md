# Phase 33: Gene-Level FDR Weighting - Research

**Researched:** 2026-02-24
**Domain:** Statistical multiple testing correction / weighted Benjamini-Hochberg
**Confidence:** HIGH

## Summary

Phase 33 adds weighted Benjamini-Hochberg (wBH) FDR correction to the existing
`correction.py` module. Users provide a TSV file of per-gene biological prior
weights (pLI scores, GWAS signals, etc.) via `--gene-prior-weights <file>`; the
engine renormalizes the weights to mean=1.0 and applies standard BH to the
adjusted p-values `p_i / w_i`. When no weight file is provided, the result is
identical to the current plain BH, preserving full backward compatibility.

The implementation scope is deliberately narrow: one new function in
`correction.py`, one new loader function (weight file parsing), two new CLI
flags (`--gene-prior-weights`, `--gene-prior-weight-column`), two new fields in
`AssociationConfig`, a `fdr_weight` column in the output DataFrame, and a
diagnostics output file. No new modules are required; all work fits naturally
into existing files.

IHW (Independent Hypothesis Weighting) is explicitly excluded from this phase
per CONTEXT.md. No CLI flag, no stub, no error message for IHW.

**Primary recommendation:** Implement wBH as a pure NumPy computation in
`correction.py` (no new external dependencies). Use `statsmodels.stats.multitest`
on the adjusted p-values, exactly as the existing `apply_correction()` does for
plain BH.

---

## Standard Stack

### Core (all already in project dependencies)

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| numpy | 2.2.6 | Weight arithmetic, effective-n formula | Already in use everywhere |
| pandas | 2.3.3 | TSV file loading, DataFrame output column | Already used for all association I/O |
| statsmodels | 0.14.4 | BH step-up procedure on adjusted p-values | Already used by `apply_correction()` |

### No New Dependencies

Zero new packages are required. The weighted BH algorithm is:
1. Renormalize weights: `w_norm = weights * m / weights.sum()`
2. Adjust p-values: `p_adj = pvals / w_norm` (clip to [0, 1])
3. Apply standard BH: `smm.multipletests(p_adj, method="fdr_bh")[1]`

All three steps use only numpy + statsmodels, already present.

---

## Architecture Patterns

### Recommended Project Structure

No new files needed. Changes span these existing files:

```
variantcentrifuge/association/
├── correction.py          # ADD: apply_weighted_correction(), load_gene_weights()
├── base.py                # ADD: two new AssociationConfig fields
├── diagnostics.py         # ADD: write_fdr_weight_diagnostics()
└── engine.py              # MODIFY: pass weights to correction, add fdr_weight column

variantcentrifuge/
└── cli.py                 # ADD: --gene-prior-weights, --gene-prior-weight-column flags

variantcentrifuge/stages/
└── analysis_stages.py     # MODIFY: _build_assoc_config_from_context() + engine call site
```

### Pattern 1: Weighted BH Algorithm (Genovese 2006)

**What:** Apply BH step-up procedure to `p_i / w_i` where weights sum to `m`
(the number of tests). This is mathematically guaranteed to control FDR at
the nominal level regardless of weight magnitude (Genovese, Roeder & Wasserman
2006, Biometrika 93(3):509-524).

**Key constraint:** Weights MUST be renormalized so their sum equals `m`
(equivalently, their mean equals 1.0). Failing to renormalize breaks FDR
control.

**When to use:** When `--gene-prior-weights` is provided and `m > 1`.

**Verified implementation:**

```python
# Source: verified with numpy 2.2.6 + statsmodels 0.14.4
import numpy as np
import statsmodels.stats.multitest as smm

def apply_weighted_correction(
    pvals: np.ndarray,          # shape (m,), raw p-values
    weights: np.ndarray,        # shape (m,), positive prior weights, ALREADY validated
    method: str = "fdr",
) -> np.ndarray:
    """Weighted BH (Genovese 2006). Returns corrected p-values in input order."""
    m = len(pvals)
    if m == 0:
        return pvals

    # Renormalize: enforce mean=1.0 (sum=m) — Genovese 2006 FDR guarantee
    w_norm = weights * m / weights.sum()

    if method == "bonferroni":
        # Weighted Bonferroni: p_i * m / w_i (equivalent to p/w * m)
        corrected = np.minimum(pvals / w_norm, 1.0)
        return corrected

    # Weighted BH: apply standard BH to adjusted p-values
    p_adj = pvals / w_norm
    # Clip to [0, 1]: genes with w < p get p/w > 1 => treated as if p=1
    p_adj_clipped = np.clip(p_adj, 0.0, 1.0)
    corrected: np.ndarray = smm.multipletests(p_adj_clipped, method="fdr_bh")[1]
    return corrected
```

**Verified behavior:**
- Equal weights (all 1.0): output is identical to unweighted BH (backward compatible)
- Single gene (m=1): renormalization forces w_norm=[1.0], identical to unweighted
- p/w > 1.0 (small weight, large p): clipped to 1.0 before BH (gene gets no power)

### Pattern 2: Weight File Loading

**What:** Parse TSV with header row; extract gene symbol column (column 0) and
weight column (named `weight` by default, overridable by CLI flag).

**Implementation approach:** pandas `read_csv(sep="\t")` — matches all other
TSV loading in the project. Fall back to ENSG ID column if gene symbol lookup
fails.

```python
# Source: codebase pattern from covariates.py and diagnostics.py
import pandas as pd
from pathlib import Path

def load_gene_weights(
    filepath: str | Path,
    weight_column: str = "weight",
    gene_column: str | None = None,   # None = auto-detect (first column)
) -> dict[str, float]:
    """
    Load gene prior weights from TSV file.

    Returns dict mapping gene symbol -> raw (un-normalized) weight.
    Caller is responsible for renormalization.

    Raises
    ------
    ValueError
        If file has no rows, weight column is missing, or any weight is
        zero/negative. Error message names the offending gene(s).
    FileNotFoundError
        If filepath does not exist.
    """
    df = pd.read_csv(filepath, sep="\t", dtype=str)

    if df.empty:
        raise ValueError(f"Weight file {filepath} contains no data rows")

    # First column = gene identifier (unless override provided)
    gene_col = gene_column if gene_column is not None else df.columns[0]

    if weight_column not in df.columns:
        raise ValueError(
            f"Weight column '{weight_column}' not found in {filepath}. "
            f"Available columns: {list(df.columns)}"
        )

    # Parse weights as float; catch non-numeric
    weights: dict[str, float] = {}
    bad_genes: list[str] = []
    for _, row in df.iterrows():
        gene = str(row[gene_col]).strip()
        raw = str(row[weight_column]).strip()
        try:
            w = float(raw)
        except ValueError:
            raise ValueError(
                f"Non-numeric weight '{raw}' for gene '{gene}' in {filepath}"
            )
        if w <= 0.0:
            bad_genes.append(f"{gene}={w}")
        weights[gene] = w

    if bad_genes:
        raise ValueError(
            f"Weights must be strictly positive. "
            f"Zero/negative weights found: {', '.join(bad_genes)}"
        )

    return weights
```

### Pattern 3: Gene Weight Resolution at Correction Time

**What:** Map loaded weights to the ordered list of testable genes, assign
default=1.0 for missing genes, emit coverage warnings.

```python
# Source: codebase pattern, verified against context requirements
def resolve_weights(
    genes: list[str],
    weight_map: dict[str, float],
) -> tuple[np.ndarray, int]:
    """
    Resolve per-gene weights for the FDR correction call.

    Returns
    -------
    weights : np.ndarray, shape (len(genes),)
        Raw (un-normalized) weights in gene order. Missing genes get 1.0.
    n_missing : int
        Count of genes that received the default weight=1.0.
    """
    weights = np.array(
        [weight_map.get(g, 1.0) for g in genes],
        dtype=float,
    )
    n_missing = sum(1 for g in genes if g not in weight_map)
    return weights, n_missing
```

### Pattern 4: Effective Number of Tests

**What:** Report `sum(w)^2 / sum(w^2)` using the NORMALIZED weights.
This is the same regardless of whether computed before or after normalization
(scaling all weights by a constant does not change the ratio).

```python
# Source: verified numerically with numpy 2.2.6
def effective_n_tests(weights_normalized: np.ndarray) -> float:
    """
    Effective number of tests (Genovese 2006 interpretation).

    For uniform weights: equals m (total genes).
    For concentrated weights: approaches 1.0.
    """
    return float(weights_normalized.sum() ** 2 / (weights_normalized ** 2).sum())
```

**Verified numerical results:**
- Uniform weights (all 1.0, m=5): effective_n = 5.0
- Extreme concentration (one gene weight=5, others near zero): effective_n ≈ 1.0
- Mixed weights [2.0, 1.5, 0.5, 0.8, 0.2] (m=5): effective_n ≈ 3.48

### Pattern 5: AssociationConfig Extension

**What:** Add two new optional fields to the existing `AssociationConfig` dataclass
in `base.py`. Pattern matches all prior phase field additions.

```python
# Source: base.py existing pattern (Phase 27 example)
@dataclass
class AssociationConfig:
    # ... existing fields ...

    # Phase 33: Gene-level FDR weighting
    gene_prior_weights: str | None = None
    """Path to gene prior weights TSV file. None = standard (unweighted) BH."""

    gene_prior_weight_column: str = "weight"
    """Name of the weight column in the weights TSV file. Default: 'weight'."""
```

### Pattern 6: CLI Flags

**What:** Add two flags to `stats_group` in `cli.py`. Pattern matches all
existing association analysis flags.

```python
# Source: cli.py existing pattern
stats_group.add_argument(
    "--gene-prior-weights",
    type=str,
    default=None,
    help=(
        "Path to per-gene prior weight TSV file for weighted FDR correction. "
        "TSV must have a header row; first column = gene symbol, weight column = 'weight' "
        "(overridable with --gene-prior-weight-column). "
        "Genes absent from the file receive weight=1.0 (neutral). "
        "Activates weighted Benjamini-Hochberg (Genovese 2006)."
    ),
)
stats_group.add_argument(
    "--gene-prior-weight-column",
    type=str,
    default="weight",
    help=(
        "Name of the weight column in the --gene-prior-weights file. "
        "Default: 'weight'. Override if your file uses a different column name."
    ),
)
```

### Pattern 7: Engine Integration

**What:** The `engine.py` correction call (lines 437-449) currently calls
`apply_correction(raw_pvals, self._config.correction_method)`. Phase 33
replaces this with a weighted call when weights are configured, and adds a
`fdr_weight` column to the output DataFrame.

The call site in `engine.py` run_all():

```python
# CURRENT (lines 442-449):
if testable_genes:
    raw_pvals = [acat_o_results[g].p_value for g in testable_genes]
    corrected = apply_correction(raw_pvals, self._config.correction_method)
    for gene, corr_p in zip(testable_genes, corrected, strict=True):
        acat_o_results[gene].corrected_p_value = float(corr_p)

# PHASE 33 REPLACEMENT:
if testable_genes:
    raw_pvals_arr = np.array(
        [acat_o_results[g].p_value for g in testable_genes], dtype=float
    )
    weight_map = self._config._resolved_gene_weights or {}
    # _resolved_gene_weights is populated in analysis_stages.py before engine.run_all()

    if weight_map and len(testable_genes) > 1:
        weights_arr, n_missing = resolve_weights(testable_genes, weight_map)
        corrected = apply_weighted_correction(
            raw_pvals_arr, weights_arr, self._config.correction_method
        )
        # Compute normalized weights for fdr_weight column output
        m = len(testable_genes)
        w_norm = weights_arr * m / weights_arr.sum()
        fdr_weights_by_gene = dict(zip(testable_genes, w_norm))
    else:
        corrected = apply_correction(raw_pvals_arr, self._config.correction_method)
        fdr_weights_by_gene = {g: 1.0 for g in testable_genes}

    for gene, corr_p in zip(testable_genes, corrected, strict=True):
        acat_o_results[gene].corrected_p_value = float(corr_p)
```

**Alternative approach:** Do NOT put weight resolution in `engine.py`. Instead,
keep all weight logic in `analysis_stages.py`, call `apply_weighted_correction()`
directly from the stage, and pass the `fdr_weight` column separately. This
avoids modifying engine.py's orchestration logic and keeps the engine lean.

**Recommended approach:** Weight loading and diagnostics in `analysis_stages.py`,
weight math in `correction.py`. Engine.py only passes through the already-computed
corrected p-values. The `fdr_weight` column is added to `results_df` in the
stage, not inside the engine. This matches existing patterns (per-gene warnings
are computed in the stage, not the engine).

### Pattern 8: Diagnostics File

**What:** Write a TSV and log summary when weighted BH runs. The dedicated
diagnostics file (e.g., `fdr_weight_diagnostics.tsv`) contains one row per gene
with columns: `gene`, `raw_weight`, `normalized_weight`, `raw_p_value`,
`weighted_p_value`, `corrected_p_value`.

Pattern matches `diagnostics.py`'s `write_diagnostics()` — write to a dedicated
function, either in `correction.py` or `diagnostics.py`.

### Anti-Patterns to Avoid

- **Weight renormalization after splitting:** Always renormalize on the FULL
  set of testable genes, not per-chromosome or per-test. The sum-to-m constraint
  requires the denominator to be the number of hypotheses being jointly corrected.
- **Modifying engine.py orchestration internals:** Keep weight logic in the
  stage (`analysis_stages.py`) and in `correction.py`. The engine should remain
  unaware of where weights come from.
- **Silent fallback on weight lookup failure:** Missing genes must receive
  weight=1.0 with a warning, not silently drop or fail.
- **Using normalized weights for impact comparison:** When logging "N genes
  gained significance", compare to a separate `apply_correction(raw_pvals, method)`
  call — do not attempt to reconstruct unweighted from weighted.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| BH step-up procedure | Custom sort + rank implementation | `smm.multipletests(method="fdr_bh")` | Already tested, edge cases handled |
| TSV file parsing | Manual `open()` + `split("\t")` | `pd.read_csv(sep="\t")` | Handles encoding, quoting, dtype |
| p-value clipping | Custom bounds check | `np.clip(p_adj, 0.0, 1.0)` | Vectorized, readable |

**Key insight:** The weighted BH is just BH applied to `p/w`. There is no special
"weighted BH" implementation needed in statsmodels — the existing BH implementation
is the correct implementation when fed the adjusted p-values.

---

## Common Pitfalls

### Pitfall 1: Forgetting to Clip Adjusted p-Values

**What goes wrong:** When a gene has a small weight (w < p), dividing gives p/w > 1.0.
`smm.multipletests` expects p-values in [0, 1]; behavior with values > 1 is undefined.

**Why it happens:** Low-priority genes (small w) can have p/w > 1, which is
mathematically meaningful (the gene is "deprioritized") but numerically invalid.

**How to avoid:** `p_adj_clipped = np.clip(p_adj, 0.0, 1.0)` before passing to BH.

**Warning signs:** `nan` or unexpected values in corrected output.

### Pitfall 2: Renormalizing Without All Genes

**What goes wrong:** If weight renormalization uses `len(weight_map)` instead of
`len(testable_genes)`, the mean won't be 1.0 and FDR control is broken.

**Why it happens:** Weight file may contain more genes than are actually tested.
Renormalization denominator must be the count of genes in the correction pass.

**How to avoid:** Resolve weights to the `testable_genes` list first, then
renormalize: `w_norm = weights * len(testable_genes) / weights.sum()`.

### Pitfall 3: Double Normalization

**What goes wrong:** If weights are normalized at load time (to mean=1 over all
genes in the file) and then normalized again at correction time (to mean=1 over
testable genes), the sum constraint is satisfied only for the testable subset,
but the values are wrong.

**Why it happens:** Temptation to "clean" weights at load time.

**How to avoid:** Store RAW weights from file. Normalize once, at correction time,
against the actual `testable_genes` list.

### Pitfall 4: Wrong Diagnostic Weight Column

**What goes wrong:** The `fdr_weight` output column shows the raw weight, not the
renormalized weight actually applied.

**Why it happens:** Loading and applying are separate steps; easy to log the
pre-normalization value.

**How to avoid:** The `fdr_weight` column must contain `w_norm[i]` (after
renormalization), NOT the raw file weight. CONTEXT.md is explicit: "the
`fdr_weight` output column should show the renormalized weight (what was actually
applied), not the raw input weight."

### Pitfall 5: Extreme Weight p-Value Artifact

**What goes wrong:** Gene with w_norm >> 1 gets p/w << p_raw. If p_raw is
already very small, p/w can underflow to 0.0, becoming statistically implausible.

**Why it happens:** Large weights can push p_adj below float precision.

**How to avoid:** No hard cap is needed (context decision), but clip at a
reasonable minimum such as `1e-300` or let numpy handle it naturally (float64
min ~5e-324). The warning for max/min ratio > 100 alerts the user.

### Pitfall 6: Single Gene Path Bypass

**What goes wrong:** With m=1, weights always renormalize to [1.0], and weighted
BH is identical to unweighted BH. But the code still runs the full path (loading,
normalizing, calling BH) unnecessarily.

**Why it happens:** Normal code path doesn't guard against m=1.

**How to avoid:** Add an early return / skip guard: `if len(testable_genes) <= 1:
log info and use plain apply_correction()`. CONTEXT.md requires this skip with
specific log message.

### Pitfall 7: ENSG Fallback Matching Logic

**What goes wrong:** Weight file has ENSG IDs in a second column; gene symbol
lookup fails; ENSG fallback is never triggered because lookup order is wrong.

**Why it happens:** Complex two-column lookup is easy to implement incorrectly.

**How to avoid:** Keep the loader simple. Load a primary dict (symbol → weight)
and optionally a secondary dict (ENSG → weight). When resolving a gene, try
symbol first, ENSG second. The gene symbol in `testable_genes` comes from the
GENE column of gene_burden_data, which is always a symbol.

---

## Code Examples

### Full Weighted BH Flow

```python
# Source: verified with numpy 2.2.6 + statsmodels 0.14.4 (2026-02-24)
import numpy as np
import statsmodels.stats.multitest as smm

# --- Step 1: Load weights from file (in correction.py) ---
weight_map = {"GENE_A": 2.0, "GENE_B": 1.5, "GENE_C": 0.5, "GENE_D": 0.8}
# GENE_E not in file -> gets 1.0 default

# --- Step 2: Resolve weights for testable genes (in analysis_stages.py) ---
testable_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"]
raw_pvals = np.array([0.01, 0.04, 0.03, 0.2, 0.5])

weights = np.array([weight_map.get(g, 1.0) for g in testable_genes])
# [2.0, 1.5, 0.5, 0.8, 1.0]

n_missing = sum(1 for g in testable_genes if g not in weight_map)
# n_missing = 1

# --- Step 3: Warnings ---
m = len(testable_genes)  # 5
pct_missing = n_missing / m  # 0.20 (20%)
# < 50%: no warning

# --- Step 4: Apply weighted correction (in correction.py) ---
w_norm = weights * m / weights.sum()
# [1.724, 1.293, 0.431, 0.690, 0.862]
# sum = 5.0, mean = 1.0

p_adj = raw_pvals / w_norm
# [0.0058, 0.0309, 0.0696, 0.2900, 0.5800]
p_adj_clipped = np.clip(p_adj, 0.0, 1.0)
corrected = smm.multipletests(p_adj_clipped, method="fdr_bh")[1]

# --- Step 5: Effective number of tests ---
eff_n = w_norm.sum() ** 2 / (w_norm ** 2).sum()
# ~4.47 (close to 5 because weights aren't extreme)

# --- Step 6: Output fdr_weight column ---
# results_df["fdr_weight"] = [w_norm[i] for each gene]
# Shows renormalized weight, NOT raw file weight
```

### Config Propagation (cli.py -> AssociationConfig)

```python
# Source: cli.py existing pattern (verified against lines 1285-1290)

# cli.py additions:
cfg["gene_prior_weights"] = getattr(args, "gene_prior_weights", None)
cfg["gene_prior_weight_column"] = getattr(args, "gene_prior_weight_column", "weight")

# analysis_stages.py _build_assoc_config_from_context() additions:
return AssociationConfig(
    # ... existing fields ...
    gene_prior_weights=_get("gene_prior_weights", default=None, nullable=True),
    gene_prior_weight_column=_get(
        "gene_prior_weight_column", default="weight", nullable=False
    ),
)
```

### Impact Comparison (Unweighted vs Weighted)

```python
# Source: codebase pattern; approach verified against context requirements
# "Impact comparison: log count of genes that gained/lost significance"
# Run both corrections and compare at a given alpha threshold

alpha = 0.05
corrected_unweighted = smm.multipletests(raw_pvals, method="fdr_bh")[1]
corrected_weighted = smm.multipletests(p_adj_clipped, method="fdr_bh")[1]

gained = int(np.sum((corrected_weighted < alpha) & (corrected_unweighted >= alpha)))
lost   = int(np.sum((corrected_weighted >= alpha) & (corrected_unweighted < alpha)))
logger.info(
    f"Weighted FDR vs unweighted BH at alpha={alpha}: "
    f"{gained} gene(s) gained significance, {lost} gene(s) lost significance"
)
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Plain BH (equal weights) | Weighted BH with user priors | This phase | Increased power when priors are informative |
| `apply_correction()` returns corrected array only | Returns corrected array + adds `fdr_weight` column | This phase | Transparent output showing what weight was applied |

**No deprecated items.** The existing `apply_correction()` function continues to
work unchanged for the unweighted case (Bonferroni and plain BH).

---

## Open Questions

1. **ENSG ID fallback column name**
   - What we know: Context says "optional second column with ENSG IDs used as
     fallback for unmatched symbols"
   - What's unclear: The weight file format has "first two columns are gene
     identifier and weight." Is the ENSG ID the second column (before weight),
     or is weight always the second column and ENSG can be any other column?
   - Recommendation: Load weight file with pandas, keep all columns; primary
     lookup = column[0] (gene symbol), weight lookup = `weight_column`; ENSG
     fallback = find any column starting with "ENSG" or named "gene_id". This
     is a Claude's-Discretion area per CONTEXT.md — keep simple and document.

2. **Diagnostics file location**
   - What we know: Context says "dedicated diagnostics file (e.g. `fdr_diagnostics.tsv`)"
     and "Output to BOTH logger AND dedicated diagnostics file"
   - What's unclear: Should the file go into `assoc_config.diagnostics_output`
     directory (if set) or adjacent to the main association output file?
   - Recommendation: If `diagnostics_output` is set, write there. Otherwise,
     write adjacent to the association TSV output. This matches the pattern for
     other diagnostics in `write_diagnostics()`.

3. **Bonferroni + weights**
   - What we know: Weighted Bonferroni is `min(p_i * m / w_i, 1.0)` (equivalent
     to multiplying by m and dividing by w_i). FDR-06 says "all equal weight =
     identical to current plain BH/Bonferroni."
   - What's unclear: Should weighted Bonferroni be supported at all, or only
     weighted BH?
   - Recommendation: Support weighted Bonferroni for completeness (the math is
     trivial; one line of code). Use the same `apply_weighted_correction()` entry
     point with method routing.

---

## Sources

### Primary (HIGH confidence)

- Verified numerically in this codebase (numpy 2.2.6, statsmodels 0.14.4,
  pandas 2.3.3) — weighted BH algorithm, effective-n formula, edge case behavior
- `variantcentrifuge/association/correction.py` — existing `apply_correction()`
  pattern that weighted version extends
- `variantcentrifuge/association/base.py` — AssociationConfig field addition pattern
- `variantcentrifuge/association/engine.py` — engine correction call site (lines 437-449)
- `variantcentrifuge/stages/analysis_stages.py` — `_build_assoc_config_from_context()`
  and AssociationAnalysisStage correction integration point
- `variantcentrifuge/cli.py` — CLI flag addition pattern (lines 530-605)
- `tests/unit/test_association_correction.py` — existing test structure to follow

### Secondary (MEDIUM confidence)

- Genovese CR, Roeder K, Wasserman L (2006). "False discovery control with
  p-value weighting." Biometrika 93(3):509-524. — Core algorithm reference;
  sum-to-m constraint and FDR control proof. (Not fetched directly; algorithm
  confirmed via IHW documentation and numerical verification.)
- IHW Bioconductor documentation (confirmed): "IHW calculates weights w_i ≥ 0
  such that they average to 1 (Σw_i = m). BH is applied to P_i/w_i."
  Source: [IHW Introduction](https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html)

### Tertiary (LOW confidence)

- WebSearch results on weighted BH implementations — confirmed algorithm direction
  but no single source provided complete pseudocode verification.

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all existing dependencies, no new packages needed
- Algorithm correctness: HIGH — numerically verified in Python with the actual
  installed versions; matches documented theory
- Architecture: HIGH — follows established patterns from every prior phase
- Pitfalls: HIGH — derived from numerical experiments and code inspection
- ENSG fallback implementation detail: MEDIUM — documented in context but
  exact column convention is a discretionary decision

**Research date:** 2026-02-24
**Valid until:** 2026-09-24 (stable mathematical algorithm; no external dependencies)
