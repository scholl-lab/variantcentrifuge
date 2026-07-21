# Cohort-Stratified Meta-Analysis

VariantCentrifuge can combine two or more clinically distinct case-control
**strata** into one calibrated rare-variant discovery analysis using
**score-statistic meta-analysis** (issue #112). Unlike pooling with a cohort
covariate, each stratum keeps its own null model, ascertainment, and case:control
ratio; unlike combining per-gene p-values, the signed score and covariance are
preserved so burden direction and kernel (SKAT) evidence are retained.

> **Scope.** This is an analysis-engine feature for **binary traits**. Study-specific
> phenotype eligibility, covariate choice, and QC remain caller responsibilities.
> Quantitative traits, a mixed-model/GRM null, and a genotype-level saddlepoint
> (SPA) for *severe* imbalance are deferred (see *Limitations*).

## How it works

The workflow mirrors RAREMETAL / MetaSKAT / seqMeta: each stratum exports a
compact **score/covariance artifact**, then a separate step combines them.

1. **Per-stratum artifact.** For every gene, VariantCentrifuge stores the
   unweighted score vector `U = Gᵀr` (r = null residuals) and the
   projection-adjusted covariance `V = GᵀPG`, plus per-variant observed ALT-allele
   counts (`AC`) and allele numbers (`AN`). These are gene-level **aggregates** —
   no per-sample genotypes or phenotype values — so artifacts are safe to move
   between sites.
2. **Meta step (`variantcentrifuge-meta`).** For each gene it aligns the strata on
   a canonical `chrom:pos:ref:alt` key (union of variants; a variant absent in a
   stratum contributes zero information there), harmonises weights from the
   **pooled** frequency, and combines:
   - **meta-burden** — common-effect score burden (`Z = 1ᵀS / √(1ᵀΣ1)`), with a
     signed direction;
   - **meta-SKAT** — kernel test via the Davies mixture chain;
   - **meta-SKAT-O** — the optimal omnibus (the reported *primary* result);
   - **heterogeneity** — Cochran's Q, I², per-stratum directions, and
     leave-one-stratum-out.
   Exome-wide multiplicity correction (BH/Bonferroni) is applied to the primary
   p-value across genes.

Weights are harmonised from the pooled ALT frequency `ΣAC/ΣAN`, using
`MAF = min(af, 1-af)` for the `Beta(MAF; 1, 25)` weight. This is essential: a
stratum's own MAF would put its `U`/`V` on an incompatible scale.

### ACAT-O vs meta-analysis

These are **different** operations and are kept distinct:

- **Within-gene ACAT-O** (see the Association Testing guide) combines *tests/masks
  within one cohort* via Cauchy combination.
- **Cross-stratum meta-analysis** (this guide) combines the *score and covariance*
  across cohorts. Per-cohort p-value combination is **never** used as a silent
  substitute for the score-based common-effect test; an explicit ACAT companion is
  available as a heterogeneity-robust complement only.

## Study manifest

A TSV, one row per stratum:

```
stratum_id	artifact_path
adult_cohort	/data/meta/adult.npz
pediatric_cohort	/data/meta/pediatric.npz
```

The meta step **refuses** (hard error, no silent fallback) when strata disagree on
genome build, annotation version, mask definition, or trait type, or when their
(hashed) sample sets overlap — overlap would violate the additive-covariance
independence assumption. Sample non-overlap does not by itself establish
independence for related or shared-control populations; that remains a caller
assumption.

## Running it

```bash
variantcentrifuge-meta \
  --manifest study.tsv \
  --output meta_results.tsv \
  --weights beta:1,25 \
  --correction fdr
```

Only `beta:a,b` and `uniform` weights are supported (the only schemes reproducible
from the artifact). Output columns include `meta_burden_p`,
`meta_burden_direction`, `meta_skat_p`, `meta_skato_p`, `meta_skato_rho`, `het_q`,
`het_p`, `het_i2`, `per_stratum_directions`, `loo_min_p`, `primary_p`,
`primary_q`, and a `status`/`warnings` pair (e.g. a non-converged stratum is
excluded with an explicit warning).

## API

```python
from variantcentrifuge.association.meta import (
    compute_gene_score, meta_burden, meta_skat, meta_skato,
)
```

`compute_gene_score(imputed_matrix, raw_dosage, null_model, variant_keys)` builds a
`GeneScore`; the `meta_*` functions combine aligned strata. A single-stratum call
reproduces the single-cohort SKAT/SKAT-O/burden result exactly — the correctness
gate that ties the meta engine to the validated backend.

## Calibration and limitations

The common-effect score meta uses a normal/Davies approximation that is well
calibrated under **moderate** imbalance. Under **severe** case:control imbalance the
per-variant/burden normal approximation is anti-conservative in the tail; a
calibrated result there requires a genotype-level saddlepoint (SPA) backend, which
is **deferred**. The artifact stores the raw `U`/`V`, so an SPA meta path can be
added without changing the schema. Likewise, `V` is stored as a full matrix, so a
future relatedness-aware (GLMM/GRM) null drops in without a schema change.

Independent cross-validation against MetaSKAT/seqMeta on shared synthetic data is
recommended before publishing exome-wide claims (R is not guaranteed in CI, so the
in-repo reference values are Python-native regression anchors).
