# Phase 25: Python Default Backends and Quick Wins — Research

**Researched:** 2026-02-22
**Domain:** Association testing backend configuration, p-value fallback chain, ACAT omnibus test
**Confidence:** HIGH (all findings verified directly from codebase)

---

## Summary

Phase 25 consists of four targeted changes to the association framework, all in Wave 1
(defaults + quick wins). Research was performed entirely by reading the existing codebase
— no external library research needed because all work involves modifying already-written
Python files.

The four changes are:
1. **Backend default swap** — `auto` mode in `get_skat_backend()` and `engine.py`
   currently tries R first. Change to return Python directly.
2. **Deprecation warnings** — Add `warnings.warn(DeprecationWarning)` in the R backend
   `__init__()` methods of `RSKATTest` and `COASTTest`.
3. **Saddlepoint-before-Liu fallback** — In `compute_pvalue()` (davies.py), when Davies
   returns out-of-range (p > 1 or p <= 0), try `_kuonen_pvalue()` before `_liu_pvalue()`.
   This is a 3-line change.
4. **ACAT-V per-variant score test** — New function `compute_acat_v()` in `acat.py`, plus
   wiring into `compute_acat_o()` and the engine. Needs the SKAT null model residuals and
   genotype matrix, which are available in the engine's per-gene context via
   `contingency_data["genotype_matrix"]` and `null_model.extra["residuals"]`.

**Primary recommendation:** All four changes are independent and can be implemented in
exactly the files identified below. No new modules, no new dependencies, no API changes
outside the listed files.

---

## Standard Stack

No new libraries are introduced in this phase. All changes use existing imports.

### Existing imports already in scope

| Module | Purpose | Used In |
|--------|---------|---------|
| `warnings` (stdlib) | `warnings.warn()` for deprecation | `r_backend.py`, `skat_r.py`, `allelic_series.py` |
| `scipy.stats.norm` | `norm.sf()` for per-variant score test p-values | `python_backend.py` (already imported) |
| `numpy` | Array operations for ACAT-V | `acat.py` (already imported) |
| `variantcentrifuge.association.backends.davies` | `_kuonen_pvalue()` | `davies.py` (internal) |
| `variantcentrifuge.association.weights` | `beta_maf_weights()` | `python_backend.py` (already imported) |

---

## Architecture Patterns

### Part 1: Backend Default Swap

**Files:** `backends/__init__.py`, `engine.py`, `base.py`, `cli.py`

**Current `auto` behavior in `get_skat_backend()` (`backends/__init__.py` lines 73-95):**
```python
if backend_name == "auto":
    try:
        from variantcentrifuge.association.backends.r_backend import RSKATBackend
        backend = RSKATBackend()
        try:
            import rpy2.robjects
            _ = rpy2.robjects
            return backend        # R-first: returns R backend if rpy2 is importable
        except Exception:
            pass
    except Exception:
        pass
    from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
    return PythonSKATBackend()    # Only reaches here if R fails
```

**New `auto` behavior:**
```python
if backend_name == "auto":
    # Python is now the default; no R probe needed
    from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
    return PythonSKATBackend()
```

**Current `auto` behavior in `engine.py` (`from_names()` lines 130-173):**

SKAT swap block (lines 130-148):
```python
skat_backend = getattr(config, "skat_backend", "r")
if skat_backend == "python":
    registry["skat"] = PurePythonSKATTest
elif skat_backend == "auto":
    try:
        get_skat_backend("r")   # probe R
        # keep RSKATTest if R is available
    except ...:
        registry["skat"] = PurePythonSKATTest  # only swaps if R fails
```

COAST swap block (lines 150-173):
```python
coast_backend = getattr(config, "coast_backend", "auto")
if coast_backend == "python":
    registry["coast"] = PurePythonCOASTTest
elif coast_backend == "auto":
    try:
        import rpy2.robjects
        importr("AllelicSeries")
        # keep COASTTest if R+AllelicSeries available
    except ...:
        registry["coast"] = PurePythonCOASTTest  # only swaps if R fails
```

**New `auto` behavior in engine.py — swap immediately without R probe:**
```python
skat_backend = getattr(config, "skat_backend", "python")  # default changes too
if skat_backend in ("python", "auto"):          # auto now means Python
    registry["skat"] = PurePythonSKATTest
# "r" branch: unchanged (explicit opt-in)

coast_backend = getattr(config, "coast_backend", "python")
if coast_backend in ("python", "auto"):
    registry["coast"] = PurePythonCOASTTest
# "r" branch: unchanged
```

**`AssociationConfig` defaults in `base.py` (lines 139, 173):**
```python
skat_backend: str = "auto"    # Change to "python"
coast_backend: str = "auto"   # Change to "python"
```

**CLI defaults in `cli.py` (lines 424-434):**
```python
"--skat-backend",
choices=["auto", "r", "python"],
default="auto",   # Change to "python"
help="...",

"--coast-backend",
choices=["auto", "r", "python"],
default="auto",   # Change to "python"
help="...",
```

**Note on `"auto"` semantics after the change:** The string `"auto"` no longer tries R
first — it becomes an alias for `"python"`. The docstring in `get_skat_backend()` and the
`AssociationConfig` docstrings must be updated to reflect this. The `"r"` value remains
fully functional as an explicit opt-in.

### Part 2: Deprecation Warnings

**Files:** `backends/r_backend.py` (`RSKATBackend.__init__`), `tests/skat_r.py`
(`RSKATTest.__init__`), `tests/allelic_series.py` (`COASTTest.__init__`)

**Pattern:**
```python
import warnings

class RSKATBackend(SKATBackend):
    def __init__(self) -> None:
        warnings.warn(
            "RSKATBackend is deprecated and will be removed in v0.17.0. "
            "Use PythonSKATBackend (--skat-backend python, now the default).",
            DeprecationWarning,
            stacklevel=2,
        )
        # ... existing init code ...
```

The same pattern applies to `RSKATTest.__init__()` and `COASTTest.__init__()`.

**Important:** `stacklevel=2` is correct here because `__init__` is called by the caller
of `RSKATBackend()` / `RSKATTest()` — typically `from_names()` in the engine or
`check_dependencies()`. Using `stacklevel=2` points the warning to the call site.

**Warning message guidance:**
- Include the removal version: `"will be removed in v0.17.0"`
- Include the alternative: `"Use --skat-backend python (now the default)"`
- Keep concise (one sentence per warning)

### Part 3: Saddlepoint-Before-Liu Fallback

**File:** `backends/davies.py`, function `compute_pvalue()` (lines 391-439)

**Current fallback at lines 401-409 (when Davies returns out-of-range):**
```python
# R: if p > 1 or p <= 0, fall back to Liu
if p_val > 1.0 or p_val <= 0.0:
    logger.info(
        "Davies p=%s out of range (ifault=%d); using Liu fallback (p=%s)",
        p_val, ifault, p_liu,
    )
    return p_liu, "liu", False
```

**New fallback (insert saddlepoint attempt before Liu):**
```python
# If Davies returns out-of-range, try saddlepoint before Liu
if p_val > 1.0 or p_val <= 0.0:
    p_sp = _kuonen_pvalue(q, lambdas)
    if p_sp is not None and 0.0 < p_sp < 1.0:
        logger.info(
            "Davies p=%s out of range; saddlepoint p=%s used.", p_val, p_sp
        )
        return float(p_sp), "saddlepoint", False
    logger.info(
        "Davies p=%s out of range (ifault=%d); saddlepoint also failed; "
        "using Liu fallback (p=%s)", p_val, ifault, p_liu,
    )
    return p_liu, "liu", False
```

**Important implementation note:** The `_PROACTIVE_THRESHOLD` constant (line 53,
`1e-5`) is **defined but not currently used** in `compute_pvalue()`. The module docstring
(line 19) and IMPL-26 specify "proactive saddlepoint at p<=1e-5 even when Davies
ifault=0." This was declared as a decision but never implemented. Phase 25 scope is
ONLY the out-of-range fallback (Part 3.1), not the proactive threshold — but planners
should be aware that `_PROACTIVE_THRESHOLD` exists unused and the proactive trigger
is separate from the out-of-range fallback being added here.

**No API changes:** `compute_pvalue()` signature and return type are unchanged. The
`p_method` return value now includes `"saddlepoint"` for a new case (out-of-range
Davies recovery), but this was already a documented valid return value.

### Part 4: ACAT-V Per-Variant Score Test

**Files:** `tests/acat.py` (new function `compute_acat_v()`), `engine.py` (wire in
ACAT-V result to `_compute_acat_o()`), `tests/skat_python.py` or directly in `acat.py`.

**Mathematical definition (from ASSOCIATION_REMAINING_WORK.md, Part 3.6):**

For each variant j in a gene with n_variants variants:
- Score: `S_j = G_j' @ residuals` (dot product of genotype column j with null model residuals)
- Variance: `V_j = sigma2 * G_j' @ G_j` (for quantitative) or `G_j' @ diag(mu*(1-mu)) @ G_j` (for binary)
- Per-variant p-value: `p_j = 2 * norm.sf(abs(S_j) / sqrt(V_j))`
- ACAT-V: `p_acat_v = cauchy_combination(p_1, ..., p_M, weights=beta_maf_weights(mafs))`

**Key question: Where is ACAT-V computed?**

ACAT-V requires:
1. Genotype matrix `G` (n_samples, n_variants) — available in engine via `contingency_data["genotype_matrix"]`
2. Null model residuals — available from the fitted null model in `PurePythonSKATTest`
3. `sigma2` / `mu_hat` for variance computation — in `null_model.extra`

ACAT-V is NOT a separate AssociationTest registered in the engine — it is a component of
ACAT-O. Per decision IMPL-34, ACAT-O is a post-loop meta-test only, not in `_TEST_REGISTRY`.

**Two implementation options:**

Option A (simpler, preferred per ASSOCIATION_REMAINING_WORK.md): Implement `compute_acat_v()`
as a pure function in `acat.py` that takes `geno`, `residuals`, `sigma2_or_variance`, `mafs`,
and `trait_type`. Wire it into `_compute_acat_o()` in `engine.py` by passing genotype matrix
data from `gene_data` (which already contains `"genotype_matrix"` when SKAT is registered).

Option B (alternative): Add a new `"acat_v"` registered test that runs per-gene alongside
other tests, then feeds into ACAT-O. This contradicts IMPL-34 (ACAT-O components are not
registered tests).

**Recommended approach: Option A.** Implement `compute_acat_v()` in `acat.py` as a utility
function. Call it from `_compute_acat_o()` in `engine.py`.

**Complication:** The current `_compute_acat_o()` only receives `results_by_test` and
`sorted_data`. The genotype matrix is in `sorted_data[i]["genotype_matrix"]` (available when
SKAT/COAST tests ran). The null model `extra` (residuals, sigma2) is NOT currently threaded
into `_compute_acat_o()`. Options:

1. Pass `null_model` reference to `_compute_acat_o()` — requires storing it on the engine
   from the SKAT test's `prepare()`/`finalize()` lifecycle.
2. Compute per-variant score stats in the SKAT test's `run()` method and store in
   `result.extra["acat_v_p"]` — then `_compute_acat_o()` can read it from `results_by_test`.
3. Have `compute_acat_v()` be called in the SKAT test's `run()` and store the result in
   `extra`, so it flows through naturally.

**Option 3 is cleanest:** `PurePythonSKATTest.run()` already has access to both the
genotype matrix and null model residuals. It can call `compute_acat_v()` and store
`result.extra["acat_v_p"] = p_acat_v`. Then `_compute_acat_o()` reads it from
`results_by_test["skat"]["gene_name"].extra["acat_v_p"]` (or `"skat_python"`).

**Variance computation for ACAT-V:**

For quantitative traits: `V_j = sigma2 * (G_j' @ G_j)` where `G_j` is the raw (unweighted)
column. This is the simplest form.

For binary traits: `V_j = G_j' @ diag(variance) @ G_j` where `variance = mu_hat * (1 - mu_hat)`.
Both `mu_hat` and `sigma2` are in `null_model.extra` (set in `fit_null_model()` lines 621-630).

**Weight for ACAT-V:** `beta_maf_weights(mafs, a1, a2)` — same Beta(MAF; 1, 25) weights
already used for SKAT. This is already imported in `python_backend.py`.

**The `compute_acat_v()` function signature:**
```python
def compute_acat_v(
    geno: np.ndarray,           # (n_samples, n_variants) — raw unweighted
    residuals: np.ndarray,      # (n_samples,)
    trait_type: str,            # "binary" or "quantitative"
    sigma2: float,              # from null_model.extra["sigma2"]
    mu_hat: np.ndarray | None,  # from null_model.extra["mu_hat"], None for quantitative
    weights: np.ndarray,        # (n_variants,) — Beta(MAF) weights
) -> float | None:
```

**Example implementation:**
```python
def compute_acat_v(geno, residuals, trait_type, sigma2, mu_hat, weights):
    # Score vector: S_j = G_j' r
    score_vec = geno.T @ residuals          # (n_variants,)

    # Variance vector: V_j
    if trait_type == "binary" and mu_hat is not None:
        variance_obs = mu_hat * (1.0 - mu_hat)   # (n_samples,)
        var_vec = np.einsum('ij,i,ij->j', geno, variance_obs, geno)  # G_j' diag(v) G_j
    else:
        var_vec = sigma2 * np.einsum('ij,ij->j', geno, geno)         # sigma2 * ||G_j||^2

    # Guard: zero-variance variants excluded
    valid = var_vec > 0.0
    if not np.any(valid):
        return None

    z_stats = score_vec[valid] / np.sqrt(var_vec[valid])
    p_vals = 2.0 * scipy.stats.norm.sf(np.abs(z_stats))

    return cauchy_combination(p_vals, weights=weights[valid])
```

**Note:** `scipy.stats.norm` is already imported in `python_backend.py`. `cauchy_combination`
is already in `acat.py` (same file where `compute_acat_v` will live).

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Per-variant score test | Custom chi-squared implementation | `scipy.stats.norm.sf(abs(z))` via existing pattern | Already implemented in `_test_burden()` lines 937-938, exact same formula |
| Beta(MAF) weights | Custom beta density | `beta_maf_weights()` from `weights.py` | Already imported in `python_backend.py` |
| Cauchy combination | Custom p-value combiner | `cauchy_combination()` from `acat.py` | Already implemented, handles tiny p-values |
| Saddlepoint approximation | New implementation | `_kuonen_pvalue()` already in `davies.py` | Already implemented and validated |

---

## Common Pitfalls

### Pitfall 1: DeprecationWarning stacklevel

**What goes wrong:** Using `stacklevel=1` in `warnings.warn()` inside `__init__()` points
the warning to the `__init__` method line, not the call site. Users see the internal file
path instead of their code.

**How to avoid:** Use `stacklevel=2` in `__init__()`. This makes the warning point to
where `RSKATBackend()` / `RSKATTest()` / `COASTTest()` is constructed, which is in
`from_names()` or `check_dependencies()`.

**Note:** If the deprecation warning is added to both `RSKATBackend.__init__()` AND
`RSKATTest.__init__()` (or its `check_dependencies()`), it will fire twice because
`RSKATTest.check_dependencies()` calls `get_skat_backend("r")` which constructs
`RSKATBackend()`. Consider adding the warning only at the test-wrapper level (RSKATTest,
COASTTest) to avoid double-firing.

### Pitfall 2: `auto` semantics must be consistent across 4 locations

**What goes wrong:** Changing the `auto` default in `base.py` but not in `cli.py`, or
changing `engine.py` but not `backends/__init__.py`, leaves inconsistent behavior.

**All 4 locations that handle `auto`:**
1. `backends/__init__.py` — `get_skat_backend("auto")` (line 73)
2. `engine.py` — SKAT swap block (line 137)
3. `engine.py` — COAST swap block (line 159)
4. `base.py` — `AssociationConfig.skat_backend` default (line 139) and `coast_backend` default (line 173)
5. `cli.py` — `--skat-backend` default (line 426) and `--coast-backend` default (line 432)

### Pitfall 3: ACAT-V with no genotype matrix

**What goes wrong:** ACAT-V requires a genotype matrix in `contingency_data`. If only
Fisher or burden tests are registered (no SKAT), there is no null model and potentially
no genotype matrix.

**How to avoid:** In `_compute_acat_o()`, make ACAT-V computation conditional on the SKAT
test being registered AND having produced a result with genotype data. If SKAT is not
registered, ACAT-V is silently skipped (returns `None`).

**Detection:** Check `"skat" in self._tests or "skat_python" in self._tests` before
attempting ACAT-V computation.

### Pitfall 4: ACAT-V stored in extra vs. returned as new result

**What goes wrong:** If `acat_v_p` is stored in `results_by_test["skat"][gene].extra`,
but the engine reads from `results_by_test` key `"acat_v"` (which doesn't exist),
ACAT-V is silently skipped.

**How to avoid:** Define the contract clearly. If Option 3 (store in extra) is used:
- In `PurePythonSKATTest.run()`: `result.extra["acat_v_p"] = compute_acat_v(...)`
- In `_compute_acat_o()`: read `results_by_test.get("skat", {}).get(gene, TestResult()).extra.get("acat_v_p")`

### Pitfall 5: `_PROACTIVE_THRESHOLD` confusion

**What goes wrong:** The constant `_PROACTIVE_THRESHOLD = 1e-5` (davies.py line 53) is
defined but NOT currently used in `compute_pvalue()`. Phase 25 Part 3.1 only adds the
out-of-range fallback — it does NOT implement the proactive threshold trigger (IMPL-26
is out of scope).

**How to avoid:** Do not implement the proactive threshold as part of this phase. Only
modify the `if p_val > 1.0 or p_val <= 0.0:` branch. Leave `_PROACTIVE_THRESHOLD`
unused (it can be implemented in a future performance phase).

### Pitfall 6: Double deprecation warning

**What goes wrong:** `RSKATTest.check_dependencies()` calls `get_skat_backend("r")` which
constructs `RSKATBackend()`. If both classes have `DeprecationWarning` in `__init__()`,
every `from_names(["skat"])` call fires the warning twice.

**How to avoid:** Add `DeprecationWarning` only to `RSKATTest.__init__()` and
`COASTTest.__init__()` (the user-facing wrappers), NOT to `RSKATBackend.__init__()`.
Keep `RSKATBackend` as an internal implementation detail without user-facing deprecation.

---

## Code Examples

### Saddlepoint-before-Liu insertion (verified from codebase)

Current code in `compute_pvalue()` at lines 401-409:
```python
# Source: variantcentrifuge/association/backends/davies.py, lines 401-409
if p_val > 1.0 or p_val <= 0.0:
    logger.info(
        "Davies p=%s out of range (ifault=%d); using Liu fallback (p=%s)",
        p_val,
        ifault,
        p_liu,
    )
    return p_liu, "liu", False
```

New code (replace those 8 lines with):
```python
# Source: Phase 25 change — saddlepoint-before-Liu
if p_val > 1.0 or p_val <= 0.0:
    p_sp = _kuonen_pvalue(q, lambdas)
    if p_sp is not None and 0.0 < p_sp < 1.0:
        logger.info(
            "Davies p=%s out of range (ifault=%d); using saddlepoint p=%s",
            p_val, ifault, p_sp,
        )
        return float(p_sp), "saddlepoint", False
    logger.info(
        "Davies p=%s out of range (ifault=%d); saddlepoint failed; "
        "using Liu fallback (p=%s)",
        p_val, ifault, p_liu,
    )
    return p_liu, "liu", False
```

### Backend default swap in `get_skat_backend()` (verified from codebase)

```python
# Source: variantcentrifuge/association/backends/__init__.py
if backend_name == "auto":
    # Python is now the default backend — no R probe needed.
    # Use --skat-backend r to explicitly opt in to R backend.
    from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
    return PythonSKATBackend()
```

### Backend default swap in `engine.py` (verified from codebase)

```python
# Source: variantcentrifuge/association/engine.py — from_names() SKAT block
skat_backend = getattr(config, "skat_backend", "python")
if skat_backend in ("python", "auto"):
    from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest
    registry["skat"] = PurePythonSKATTest
    logger.debug("SKAT: using Python backend (skat_backend=%s)", skat_backend)
# elif skat_backend == "r": keep RSKATTest (explicit R opt-in, no change needed)

coast_backend = getattr(config, "coast_backend", "python")
if coast_backend in ("python", "auto"):
    from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest
    registry["coast"] = PurePythonCOASTTest
    logger.debug("COAST: using Python backend (coast_backend=%s)", coast_backend)
```

### Deprecation warning pattern (verified against stdlib warnings module)

```python
# Source: stdlib warnings module — standard pattern
import warnings

class RSKATTest(AssociationTest):
    def __init__(self) -> None:
        warnings.warn(
            "RSKATTest is deprecated and will be removed in v0.17.0. "
            "The Python SKAT backend is now the default. "
            "Remove --skat-backend r to use the default Python backend.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._backend: RSKATBackend | None = None
        # ... rest of existing __init__ ...
```

### ACAT-V function in `acat.py` (new code, uses existing imports)

```python
import scipy.stats  # already imported in acat.py via cauchy.sf
import numpy as np  # already imported

def compute_acat_v(
    geno: np.ndarray,
    residuals: np.ndarray,
    trait_type: str,
    sigma2: float,
    mu_hat: np.ndarray | None,
    weights: np.ndarray,
) -> float | None:
    """
    Compute ACAT-V p-value: per-variant score test p-values combined via Cauchy.

    For each variant j: p_j = 2 * Phi(-|S_j| / sqrt(V_j))
    where S_j = G_j' r (score) and V_j is the marginal variance.
    Combined via cauchy_combination() with Beta(MAF) weights.
    """
    if geno.shape[1] == 0:
        return None

    score_vec = geno.T @ residuals  # (n_variants,)

    if trait_type == "binary" and mu_hat is not None:
        variance_obs = mu_hat * (1.0 - mu_hat)   # (n_samples,)
        var_vec = np.array([
            float(geno[:, j] @ (variance_obs * geno[:, j]))
            for j in range(geno.shape[1])
        ])
    else:
        var_vec = sigma2 * np.einsum("ij,ij->j", geno, geno)

    valid = var_vec > 0.0
    if not np.any(valid):
        return None

    z_stats = score_vec[valid] / np.sqrt(var_vec[valid])
    p_vals = 2.0 * scipy.stats.norm.sf(np.abs(z_stats))

    return cauchy_combination(p_vals, weights=weights[valid])
```

Note: The binary-trait loop can be vectorized with `np.einsum('ij,i,ij->j', geno, variance_obs, geno)` if `n_variants` is large. The loop form is correct and readable for initial implementation.

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| R SKAT backend default | Python backend default | Phase 25 | No R dependency required by default |
| Davies → Liu (out-of-range) | Davies → Saddlepoint → Liu | Phase 25 | Improved accuracy for extreme tails (p < 1e-8) |
| ACAT-O without ACAT-V | ACAT-O includes ACAT-V component | Phase 25 | Better power for sparse signal (few causal variants) |

**Why saddlepoint over Liu for extreme tails (HIGH confidence from ASSOCIATION_REMAINING_WORK.md):**
Liu moment-matching systematically overestimates p-values when eigenvalues are non-uniform
(non-spherical kernel). Saddlepoint (Kuonen 1999, Lugannani-Rice) is accurate to
exponential precision in the tail. SAIGE-GENE+ uses SPA as its primary approximation.

---

## Open Questions

1. **ACAT-V location decision (see Pitfall 3 / Pitfall 4)**
   - What we know: ACAT-V can be computed in `PurePythonSKATTest.run()` and stored in
     `extra`, or computed on-the-fly in `_compute_acat_o()` if the engine has access to
     gene data.
   - What's unclear: Does the planner want ACAT-V wired only through the Python SKAT test,
     or should it work even when SKAT is not registered?
   - Recommendation: Wire through SKAT test `extra` for simplicity. If SKAT is not
     registered, ACAT-V is absent and ACAT-O combines whatever tests are registered.

2. **`"auto"` string deprecation**
   - Should `"auto"` become a hard alias for `"python"` (current plan) or emit a separate
     `FutureWarning` that `"auto"` previously meant "try R first"?
   - The current plan treats `"auto"` as silently equivalent to `"python"` post-change.
     This is the simplest approach and is consistent with the requirement spec.

3. **Test updates for backend default change**
   - Existing tests that check `skat_backend="auto"` behavior (R-first probe) will need
     updating. Search for `skat_backend.*auto` or `coast_backend.*auto` in tests and
     update expectations.

---

## Sources

### Primary (HIGH confidence)

- Direct codebase reading — `variantcentrifuge/association/backends/__init__.py` (factory)
- Direct codebase reading — `variantcentrifuge/association/backends/davies.py` (fallback chain)
- Direct codebase reading — `variantcentrifuge/association/tests/acat.py` (ACAT-O, cauchy_combination)
- Direct codebase reading — `variantcentrifuge/association/engine.py` (from_names, _compute_acat_o)
- Direct codebase reading — `variantcentrifuge/association/base.py` (AssociationConfig defaults)
- Direct codebase reading — `variantcentrifuge/association/backends/python_backend.py` (score test patterns)
- Direct codebase reading — `variantcentrifuge/cli.py` (CLI backend arguments)
- `.planning/ASSOCIATION_REMAINING_WORK.md` — Part 1, 3.1, 3.6 specs (authoritative scope document)

### Secondary (MEDIUM confidence)

- ASSOCIATION_REMAINING_WORK.md references: Liu & Xie (2020) ACAT — Cauchy combination
- ASSOCIATION_REMAINING_WORK.md references: SAIGE-GENE+ (SPA preferred over Liu for extreme tails)

---

## Metadata

**Confidence breakdown:**
- Backend default swap: HIGH — exact code locations verified, pattern is mechanical
- Deprecation warnings: HIGH — stdlib `warnings.warn` API is stable, stacklevel behavior verified
- Saddlepoint fallback: HIGH — `_kuonen_pvalue()` already exists, insertion point is exact
- ACAT-V implementation: MEDIUM-HIGH — math is clear, wiring location requires a design choice (see Open Questions)

**Research date:** 2026-02-22
**Valid until:** Stable — changes are mechanical refactors of existing code with no external dependencies
