# Cohort Meta-Analysis Implementation Plan

> **For agentic workers:** implement task-by-task with TDD. Steps use checkbox (`- [ ]`) syntax.
> Commits are DEFERRED — AGENTS.md forbids commits/push/PR without explicit user approval. Do the
> work on branch `feat/cohort-meta-analysis`; leave it uncommitted for review.

**Goal:** Add a score-statistic meta-analysis engine that combines per-stratum rare-variant
association evidence (burden, SKAT, SKAT-O) from persisted score/covariance artifacts.

**Architecture:** New `variantcentrifuge/association/meta/` package. Per-stratum artifacts hold
unweighted score vector `U` and projection-adjusted covariance `V` (aggregates only). A
`variantcentrifuge-meta` CLI aligns strata per gene, harmonises weights from pooled MAF, and
reuses the existing Davies chain and SKAT-O omnibus to produce meta p-values plus heterogeneity.

**Tech Stack:** Python 3.10+, numpy, scipy, statsmodels (all already deps). Reuses
`association/backends/{python_backend,davies}.py`.

## Global Constraints

- Every new production module ≤ 600 physical lines (`make lint-loc`).
- Do NOT add code to `python_backend.py` (allowlisted at 1132; must not grow) — import from it.
- Never run `variantcentrifuge` on real data. Tests use synthetic fixtures only.
- Artifacts/reports contain only aggregate quantities — never sample IDs or phenotype rows.
- Single-stratum meta must reproduce `_test_skat`/`_test_skato`/burden (golden gate).
- ruff N-series: scientific matrices use descriptive names (matches existing backend style).

---

### Task 1: Artifact schema + gene-level artifact computation

**Files:**
- Create: `variantcentrifuge/association/meta/__init__.py`
- Create: `variantcentrifuge/association/meta/artifact.py`
- Test: `tests/unit/test_meta_artifact.py`

**Interfaces — Produces:**
- `@dataclass ScoreCovarianceArtifact` with fields: `schema_version:int`, `stratum_id:str`,
  `genome_build:str`, `annotation_version:str`, `mask_id:str`, `mask_hash:str`,
  `trait_type:str`, `n_cases:int`, `n_controls:int`, `sample_id_hashes:list[str]`,
  `genes:dict[str, GeneScore]`.
- `@dataclass GeneScore` with: `variant_keys:list[str]`, `u:np.ndarray (p,)`,
  `v:np.ndarray (p,p)`, `ac:np.ndarray (p,)`, `an:np.ndarray (p,)`, `null_converged:bool`.
- `compute_gene_artifact(genotype_matrix, phenotype, covariates, variant_keys, trait_type) -> GeneScore`
- `save_artifact(artifact, path) -> None` (`.npz` arrays + `<path>.meta.json`)
- `load_artifact(path) -> ScoreCovarianceArtifact`

**Key implementation (`compute_gene_artifact`):**
```python
from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

def compute_gene_artifact(genotype_matrix, phenotype, covariates, variant_keys, trait_type):
    backend = PythonSKATBackend(); backend.detect_environment()
    nm = backend.fit_null_model(phenotype, covariates, trait_type)
    r = nm.extra["residuals"]
    g = np.asarray(genotype_matrix, dtype=np.float64)
    u = g.T @ r                                   # unweighted score  (p,)
    x = nm.extra["design_matrix"]
    if trait_type == "binary":
        pi = nm.extra["mu_hat"] * (1.0 - nm.extra["mu_hat"]); phi = np.sqrt(pi)
        z_phi = phi[:, None] * g
        x_phi = phi[:, None] * x
        x_pi = pi[:, None] * x
        z_adj = z_phi - x_phi @ np.linalg.solve(x_phi.T @ x_phi, x_pi.T @ g)
    else:
        s2 = nm.extra["sigma2"]
        z_adj = (g - x @ np.linalg.solve(x.T @ x, x.T @ g)) / np.sqrt(s2)
    v = z_adj.T @ z_adj                            # projection-adjusted cov (p,p)
    an = 2.0 * np.sum(~np.isnan(g), axis=0)        # (built before imputation in caller)
    ac = np.nansum(g, axis=0)
    return GeneScore(variant_keys, u, v, ac, an, bool(nm.extra))
```
Note: `an`/`ac` should be computed from the pre-imputation dosage; caller passes the raw matrix's
column sums. For synthetic tests the imputed matrix is acceptable (documented).

**Steps:** write failing test asserting `compute_gene_artifact` returns `u.shape==(p,)`,
`v.shape==(p,p)`, `v` symmetric PSD → run (fail) → implement → run (pass) → round-trip
save/load equality test → (commit deferred).

---

### Task 2: Golden consistency — single-stratum SKAT/SKAT-O/burden

**Files:**
- Create: `variantcentrifuge/association/meta/combine.py`
- Test: `tests/unit/test_meta_combine_golden.py`

**Interfaces — Produces:**
- `meta_skat(u_list, v_list, weights) -> tuple[float|None, str]`  (p_value, p_method)
- `meta_burden(u_list, v_list, weights) -> tuple[float|None, float, int]`  (p, beta, direction)
- `meta_skato(u_list, v_list, weights) -> tuple[float|None, float|None]`  (p, rho)

where `u_list[c]`, `v_list[c]` are ALREADY aligned to the union order (Task 3 provides alignment),
`weights` is the harmonised per-variant weight vector.

**Key implementation:**
```python
from variantcentrifuge.association.backends.davies import compute_pvalue
from variantcentrifuge.association.backends.python_backend import _skato_get_pvalue, _SKATO_RHO_GRID
import scipy.linalg, scipy.stats, numpy as np

def _weighted(u_list, v_list, w):
    s = sum(w * u for u in u_list)
    sigma = sum((w[:, None] * v) * w[None, :] for v in v_list)
    return s, sigma

def meta_skat(u_list, v_list, w):
    s, sigma = _weighted(u_list, v_list, w)
    q = float(s @ s) / 2.0
    lam = scipy.linalg.eigh(sigma / 2.0, eigvals_only=True, driver="evr")
    pos = lam[lam >= 0]
    if pos.size == 0: return None, "skip"
    lam = lam[lam > pos.mean() / 100_000.0]
    if lam.size == 0: return 1.0, "liu"
    p, method, _ = compute_pvalue(q, lam)
    return (None if not np.isfinite(p) else float(p)), method

def meta_burden(u_list, v_list, w):
    s, sigma = _weighted(u_list, v_list, w)
    score = float(s.sum()); var = float(sigma.sum())
    if var <= 0.0: return None, 0.0, 0
    z = score / np.sqrt(var)
    return float(2.0 * scipy.stats.norm.sf(abs(z))), score / var, int(np.sign(score))

def meta_skato(u_list, v_list, w):
    s, sigma = _weighted(u_list, v_list, w)
    a = sigma / 2.0
    lam, vec = scipy.linalg.eigh(a)
    keep = lam > (lam[lam >= 0].mean() / 100_000.0 if np.any(lam >= 0) else 0)
    if not np.any(keep): return None, None
    b = (vec[:, keep] * np.sqrt(lam[keep])).T          # (k,p): b.T@b == a (kept subspace)
    q_rho = np.array([((1 - r) * float(s @ s) + r * float(s.sum())**2) / 2.0
                      for r in _SKATO_RHO_GRID])
    p, per_rho = _skato_get_pvalue(q_rho, b, _SKATO_RHO_GRID)
    rho = _SKATO_RHO_GRID[int(np.argmin(per_rho))]
    return float(p), rho
```

**Golden test:** synthetic `(G, y, X)`; compute artifact `U`,`V`; weights = `Beta(maf;1,25)` from
`G.mean(0)/2`; assert `meta_skat([U],[V],w)` p ≈ `PythonSKATBackend._test_skat` p (rtol 1e-6),
`meta_skato` p ≈ `_test_skato` p, and `meta_burden` p ≈ projection-adjusted burden
(cross-check `_test_skato` per-rho at rho=1).

---

### Task 3: Variant alignment + pooled-MAF weights

**Files:**
- Create: `variantcentrifuge/association/meta/align.py`
- Test: `tests/unit/test_meta_align.py`

**Interfaces — Produces:**
- `canonical_key(chrom, pos, ref, alt) -> str`
- `align_strata(gene_scores: list[GeneScore]) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray, np.ndarray, list[str]]`
  returns aligned `u_list`, `v_list` (union order, zero-padded), pooled `ac`, `an`, `keys`.
- `pooled_maf_weights(ac, an, scheme="beta:1,25") -> np.ndarray`

Union = sorted set of keys across strata; each stratum's `U`/`V` embedded via index map
(absent variant → zero row/col). Test: two strata with overlapping+disjoint keys align correctly;
a variant present in one only contributes `U=0` in the other; pooled MAF `= ΣAC/ΣAN`.

---

### Task 4: Manifest + comparability/safety validation

**Files:**
- Create: `variantcentrifuge/association/meta/manifest.py`
- Test: `tests/unit/test_meta_manifest.py`

**Interfaces — Produces:**
- `read_study_manifest(path) -> list[StratumRef]` (columns: `stratum_id, artifact_path`)
- `validate_comparability(artifacts) -> None` — raises `MetaCompatibilityError` on mismatched
  `genome_build` / `annotation_version` / `mask_hash`, or on sample-hash overlap across strata.

Tests (acceptance #4, #5): build/annotation/mask mismatch → raises; overlapping `sample_id_hashes`
→ raises; well-formed manifest → passes; missing artifact file → clear error.

---

### Task 5: Heterogeneity + leave-one-out

**Files:**
- Create: `variantcentrifuge/association/meta/heterogeneity.py`
- Test: `tests/unit/test_meta_heterogeneity.py`

**Interfaces — Produces:**
- `cochran_q(betas, ses) -> tuple[float, float, float]`  (Q, p, I²)
- `per_stratum_burden(u_list, v_list, w) -> tuple[np.ndarray, np.ndarray]`  (betas, ses)
- `leave_one_out(u_list, v_list, w, combine_fn) -> list[tuple[str, float]]`

Test (acceptance #1, #2): concordant strata → `meta_burden_p < min(stratum_p)`; opposite-direction
strata → `het_p` significant while common-effect burden attenuated.

---

### Task 6: Runner + exome-wide multiplicity + output

**Files:**
- Create: `variantcentrifuge/association/meta/runner.py`
- Test: `tests/unit/test_meta_runner.py`

**Interfaces — Produces:**
- `run_meta(manifest_path, output_path, weight_scheme, correction) -> pd.DataFrame`

Per gene×mask: align → weights → burden/skat/skato → heterogeneity → LOO; per-gene omnibus across
masks via ACAT (`cauchy_combination`); BH/Bonferroni across genes on the omnibus; write `meta.tsv`
with the §8 schema. Test: 3-gene synthetic run produces the expected columns and a monotone qvalue.

---

### Task 7: CLI console_script

**Files:**
- Create: `variantcentrifuge/association/meta/cli.py`
- Modify: `pyproject.toml` `[project.scripts]` → add `variantcentrifuge-meta = "variantcentrifuge.association.meta.cli:main"`
- Test: `tests/unit/test_meta_cli.py`

`main()` argparse: `--manifest`, `--output`, `--weights` (default `beta:1,25`),
`--correction` (default `fdr`), `--log-level`. Reuses `run_meta`. Test invokes `main([...])` on a
synthetic manifest and asserts the output TSV exists with correct header; safe-logging test
(acceptance #6) asserts no sample IDs appear in output/logs.

---

### Task 8: Docs

**Files:**
- Create: `docs/source/guides/meta_analysis.md`
- Modify: `docs/source/guides/association_testing.md` (cross-link; ACAT-O vs meta distinction)

Explain the study manifest, artifact/output schema, weight harmonisation, the within-gene ACAT-O
vs across-stratum meta distinction, the calibration workflow, and the documented severe-imbalance
limitation (SPA deferred).

## Self-Review

- Spec §4 stats → Tasks 1,2,3. §5 safety → Tasks 3,4. §6 modules → Tasks 1–7. §8 output → Task 6.
  §9 acceptance → Tasks 2,4,5,7. §10 deferred → docs (Task 8). All covered.
- No placeholders; all load-bearing code shown. Signatures consistent across tasks
  (`u_list`,`v_list`,`w`; `GeneScore`; `_skato_get_pvalue`).
- Golden gate (Task 2) is the correctness anchor.
