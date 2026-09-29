# Design — Cohort-stratified rare-variant score-statistic meta-analysis

Issue: scholl-lab/variantcentrifuge#112
Date: 2026-07-20
Status: Design (autonomous execution under session goal directive)

## 0. Provenance and scope decision

This design was produced under an explicit session directive to brainstorm a spec
and plan, adversarially review them, and implement end-to-end. The brainstorming
skill's interactive per-section approval gate was executed autonomously (the user
is the issue author and pre-authorised the workflow); every material decision is
recorded here for audit.

**Governance note (must be ratified by a maintainer).** `.planning/PROJECT.md`
currently lists *"Meta-analysis (RAREMETAL, Meta-SAIGE)"* and *"Mixed model / GRM
(SAIGE-GENE approach)"* as **Out of Scope**. This feature intentionally moves the
first item in-scope and adds an *interface seam* (not an implementation) for the
second. `PROJECT.md` must be updated to reflect the new scope. This is a
deliberate scope change, not an oversight.

## 1. Problem

VariantCentrifuge is a single-cohort tool: one VCF, one case/control assignment,
one phenotype vector, one null model, one output sidecar. It has a mature
association engine (Fisher, logistic/linear burden with Firth fallback,
pure-Python SKAT/SKAT-O, COAST, covariates/PCA, within-gene ACAT-O) but no way to
combine two clinically distinct, potentially severely imbalanced case-control
cohorts in one calibrated exome-wide rare-variant discovery analysis without
either (a) pooling with a cohort covariate (ignores ascertainment/outcome/batch
differences) or (b) combining per-gene p-values (loses signed score and covariance
information needed for burden and kernel meta-analysis).

## 2. Key enabling observation

The sufficient statistics for score-statistic meta-analysis already exist
transiently in the SKAT backend and are discarded:

- `PythonSKATBackend._test_skat` computes `score_vec = (G∘W)ᵀ r` and eigenvalues
  of the projection-adjusted covariance, then collapses to a p-value.
- `NullModelResult.extra` carries `{residuals, sigma2, mu_hat, design_matrix}` —
  everything needed to regenerate the per-gene score vector `U` and covariance `V`
  for any variant set.

`U` (length p) and `V` (p×p) are **gene-level aggregates** — no per-sample rows —
so persisting them satisfies the issue's privacy non-goal directly. This feature
is fundamentally: *stop discarding `U` and `V`, persist them per stratum, then
combine*, reusing the existing Davies→saddlepoint→Liu chain and the existing
SKAT-O omnibus integration unchanged.

## 3. Terminology

The word "cohort" already means "one case/control study population / an aggregated
single-sample report" throughout the codebase. To avoid collision, this feature
uses **stratum** (plural: strata) for a meta-analysis contributing study, and
**study manifest** for the multi-stratum configuration.

## 4. Statistical design (verified against existing code)

For stratum `c`, with null residuals `r_c = y_c − μ̂_c` and the GLS projection
`P_c = W_c − W_c X_c (X_cᵀ W_c X_c)⁻¹ X_cᵀ W_c` (binary: `W_c = diag(μ̂(1−μ̂))`;
quantitative: `W_c = I/σ²`), and a per-gene genotype dosage matrix `G_c`
(columns = union variants, aligned by canonical key, absent → zero column):

- **Unweighted score** `U_c = G_cᵀ r_c` (length p).
- **Unweighted projection-adjusted covariance** `V_c = Z̃_cᵀ Z̃_c` (p×p), where
  `Z̃_c = φ⊙G_c − φ⊙X_c (X_cᵀ diag(φ²) X_c)⁻¹ X_cᵀ diag(φ²) G_c`, `φ = √diag(W_c)`.
  This is exactly the unweighted form of the `z_adj` matrix in
  `_compute_eigenvalues_filtered` / `_test_skato`.

Artifacts store `U_c`, `V_c` **unweighted** plus per-variant `AC_c`, `AN_c`.
Weights are harmonised at meta time from the **pooled** MAF
`maf_j = (Σ_c AC_cj)/(Σ_c AN_cj)` via `Beta(maf; 1, 25)` (or the configured
scheme), applied identically to every stratum. Rationale: cohort-specific MAF
weights put `U_c`/`V_c` on incompatible scales; pooled-MAF harmonisation is the
MetaSKAT convention.

At meta time, with harmonised weights `w` over the union:
`S = Σ_c diag(w) U_c` (weighted meta score), `Σ = Σ_c diag(w) V_c diag(w)`
(weighted meta covariance). Then:

- **Meta-SKAT:** `Q = ‖S‖²/2`, null `Q ~ Σ λ_i χ²₁` with `λ = eig(Σ/2)` →
  `davies.compute_pvalue(Q, λ)`.
- **Meta-burden (common effect):** `S_b = wᵀ Σ_c U_c`, `Var_b = wᵀ Σ_c V_c w`
  (equivalently `1ᵀS`, `1ᵀΣ1`), `Z = S_b/√Var_b`, `p = 2·Φ̄(|Z|)`,
  direction `= sign(S_b)`. This is the projection-adjusted score burden
  (more rigorous than the legacy `_test_burden`, which omits covariate projection).
- **Meta-SKAT-O:** factor `A = Σ/2 = BᵀB` (via `eigh`, dropping tiny/negative
  eigenvalues); compute `q_rho[i] = ((1−ρ)‖S‖² + ρ(1ᵀS)²)/2`; call the existing
  `_skato_get_pvalue(q_rho, B, RHO_GRID)`. This is exact because
  `_skato_optimal_param` and `_skato_get_pvalue` depend on the per-sample matrix
  **only through `A = z1ᵀz1`**; any `B` with `BᵀB = A` yields identical results.

**Golden correctness gate:** a single-stratum meta MUST reproduce
`_test_skat`/`_test_skato`/(projection-adjusted burden) p-values bit-for-bit
(within numerical tolerance) for the same data and weights. This ties all new
statistics to already-trusted code.

**Heterogeneity:** per-stratum burden effect `β_c = (wᵀU_c)/(wᵀV_c w)`,
`SE_c = 1/√(wᵀV_c w)`; Cochran's `Q = Σ (β_c−β̄)²/SE_c²`, `df = C−1`,
`I² = max(0,(Q−df)/Q)`, `p_het = χ²_df.sf(Q)`. Leave-one-stratum-out re-runs the
meta on each subset. An explicit, opt-in ACAT combination of per-stratum p-values
is offered as a heterogeneity-robust **companion** (never a silent substitute for
the score-based common-effect test).

## 5. Comparability and safety (the "never silently pool" requirement)

Each artifact records `genome_build`, an `annotation_version`, and a
`mask_definition_hash`. The meta step **refuses** (hard error, no fallback) when
strata disagree on build/annotation/mask hash. Additional checks:

- Variant keys are canonicalised (`chrom:pos:ref:alt`, normalised, left-aligned).
  Alignment is by exact key; allele orientation fixes burden sign.
- MAC filtering is applied on the **pooled** count at meta time, never per stratum
  (a variant monomorphic in one stratum correctly contributes `U=0`, zero `V`
  block — it is kept, not dropped).
- Cross-stratum **sample-ID overlap** is checked; overlap violates the additive
  covariance assumption and is refused (issue non-goal: no shared controls). Only
  hashed sample IDs are compared; raw IDs never enter artifacts or reports.
- Artifacts and reports contain only aggregate quantities (`U`, `V`, `AC`, `AN`,
  counts, convergence flags). No sample IDs, no phenotype rows. Enforced by test.

## 6. Architecture (module boundaries, all ≤600 lines per LOC policy)

New package `variantcentrifuge/association/meta/`:

- `artifact.py` — `ScoreCovarianceArtifact` dataclass; `compute_gene_artifact()`
  (builds `U`,`V`,`AC`,`AN` from a genotype matrix + `NullModelResult`, reusing
  `PythonSKATBackend.fit_null_model`); `save_artifact()` / `load_artifact()`
  (`.npz` arrays + JSON sidecar metadata; schema-versioned).
- `align.py` — canonical variant keys; union/alignment of per-stratum `U`/`V` into
  a shared variant space; pooled-MAF computation.
- `combine.py` — `meta_burden()`, `meta_skat()`, `meta_skato()` on aligned
  artifacts; imports `davies.compute_pvalue`, `_skato_get_pvalue`, `_SKATO_RHO_GRID`.
- `heterogeneity.py` — Cochran Q, I², per-stratum direction, leave-one-out,
  optional ACAT companion.
- `manifest.py` — study-manifest reader/validator (extends the dormant
  `helpers.read_sequencing_manifest`); comparability checks from §5.
- `runner.py` — orchestrate: load artifacts → align per gene×mask → combine →
  heterogeneity → wide-format results; exome-wide multiplicity (per-gene omnibus
  across masks via ACAT, then BH/Bonferroni across genes).
- `cli.py` — `variantcentrifuge-meta` console_script (new `[project.scripts]`
  entry). No VCF, no external tools; consumes artifacts only.

**Emission from the real pipeline** is a thin, later addition: an
`--emit-score-artifact` flag whose stage handler calls
`meta.artifact.compute_gene_artifact` on the per-gene matrices the association
stage already builds. To respect the 600-line ceiling on the allowlisted
`python_backend.py` (1132 lines), no code is added there; the meta package imports
its public/module-level functions.

## 7. Data flow

Per stratum (existing pipeline + emit flag, or the `compute_gene_artifact` API):
VCF → association stage → `{gene: (variant_keys, U, V, AC, AN, n_cases,
n_controls, null_convergence)}` → `stratum.artifact.npz` + `.meta.json`.

Meta (`variantcentrifuge-meta --manifest study.tsv --output meta.tsv`):
manifest → load N artifacts → validate comparability → per gene: align union,
pooled-MAF weights, `meta_burden`/`meta_skat`/`meta_skato`, heterogeneity, LOO →
exome-wide correction → `meta.tsv` (+ optional per-stratum diagnostics).

## 8. Output schema (`meta.tsv`)

`gene, mask, n_strata, n_cases_total, n_controls_total, n_variants_union,
meta_burden_p, meta_burden_beta, meta_burden_direction, meta_skat_p,
meta_skato_p, meta_skato_rho, het_q, het_p, het_i2, per_stratum_directions,
loo_min_p, primary_p, primary_q, status, warnings`. Machine-readable and
reproducible from artifacts alone.

## 9. Testing → acceptance-criteria mapping

1. **Concordant effects gain evidence** — 2 strata, same-direction burden signal;
   assert `meta_burden_p < min(stratum_p)`.
2. **Opposite directions distinguishable** — common-effect burden attenuated while
   `het_p` significant.
3. **Highly imbalanced stratum** — null simulation (modest scale) checks calibration
   with documented tolerance; stable failure reporting. The severe-imbalance
   anti-conservatism of the normal/Davies approximation is **documented as a known
   limitation** requiring the deferred SPA backend; the test asserts the honest
   behaviour, not an over-claim.
4. **Incompatible callability** — build/annotation/mask-hash mismatch → explicit
   error, never silent pooling.
5. **Config/alignment validation failures** — manifest errors and sample overlap → error.
6. **Safe aggregate outputs** — artifact/report contain no sample IDs / phenotype rows.
7. **Golden consistency** — single-stratum meta == `_test_skat`/`_test_skato`/burden.

Independent cross-validation against MetaSKAT/seqMeta on shared synthetic data is
recommended as a follow-up (R not guaranteed in CI; captured values are Python-native).

## 10. Deferred (interface-ready, not implemented)

- **Genotype-level SPA** (SAIGE-GENE style) for calibrated *severe*-imbalance
  burden/ACAT-V tails. The artifact carries the raw `U`/`V`; an SPA meta path can
  be added without schema change.
- **GLMM/GRM null backend.** `V` is stored as a full matrix (not assumed diagonal),
  so a relatedness-aware null (non-diagonal residual covariance) drops in without
  reworking the schema.

## 11. Non-goals (unchanged from issue)

No phenotype harmonisation, no control borrowing, no external controls, no raw
genotype/clinical upload, no claim that a pooled cohort covariate solves
ascertainment/outcome differences.

## 12. Revision 1 — corrections from adversarial (Codex, high-effort) review

Each item below overrides earlier text where they conflict. All were verified
against the math and the backend before acceptance.

1. **Binary-only scope.** For Gaussian residuals `Var(G'r) = σ²·G'PG`, so the
   §4 pairing `U=G'r`, `V=G'PG/σ²` is inconsistent for **quantitative** traits.
   For **binary** traits `σ²=1` and `V=G'PG=Var(U)` exactly, so the design is
   correct. This implementation is therefore **binary-only**;
   `compute_gene_artifact`/`combine` raise `NotImplementedError` for quantitative,
   which is deferred until a coherent likelihood-score convention is added.
2. **Exact SKAT-O factorisation.** Factor `A=Σ/2` keeping **all non-negative**
   eigenmodes (clip only round-off negatives to 0), so `BᵀB=A` exactly. The
   per-rho `_get_lambda` filtering stays *inside* `_skato_get_pvalue`. (Earlier
   plan Task 2 truncated positive modes — wrong.)
3. **Raw allele counts.** `compute_gene_artifact` takes the **pre-imputation**
   dosage (NaN=missing); `AN = 2·#non-missing`, `AC = nansum(raw)`. Never derive
   `AC`/`AN` from the imputed matrix (imputation forces `AN=2N` and pollutes `AC`).
4. **Real convergence.** Store `null_converged = bool(fit_result.converged)` from
   statsmodels (not `bool(nm.extra)`, which is always true). The meta step
   **excludes** non-converged strata with an explicit status and propagates the
   Davies p-method/convergence flag into results.
5. **Weight surface = Beta/uniform only.** Only pooled-frequency Beta and uniform
   weights are reproducible from artifacts; CADD/REVEL/score-column schemes are
   rejected with a clear error (would require persisting per-variant annotation
   provenance — deferred).
6. **ALT vs minor allele.** `ΣAC/ΣAN` is pooled **ALT** frequency; the Beta weight
   uses `MAF = min(af, 1−af)`. Callers must apply rare-AF filtering upstream.
7. **Heterogeneity df** = (strata with non-zero information) − 1. Manifest
   validation additionally rejects **mixed trait types** across strata.
8. Sample-hash non-overlap is checked but does **not** by itself establish
   independence for related/shared-control populations — that remains an explicit
   caller assumption (documented).
