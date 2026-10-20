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
2. **Meta step (`variantcentrifuge meta`).** For each gene it aligns the strata on
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
  across cohorts. Per-cohort p-value combination is **never** used as a
  substitute for the score-based common-effect test.

## Creating a cohort artifact

Add `--score-artifact <path>` to a normal association run of each cohort
(binary trait, case/control samples set):

```bash
variantcentrifuge ... --perform-association --association-tests skat \
  --score-artifact /data/meta/adult.npz
```

This writes `adult.npz` and `adult.npz.meta.json` next to the usual results. It
uses the same site/sample missingness QC as the association tests but keeps genes
below the per-gene MAC filter, because a gene too rare to test in one cohort still
adds information to the meta. The artifact records the reference (`--reference`)
as build and annotation version, a hash of the variant selection (`--filters`,
`--final-filter`, `--bcftools-prefilter`) as the mask, and SHA-256 hashes of the
sample IDs for the overlap check. Run every cohort with the same reference and
filters, or the meta step will refuse to pool them.

The artifact (schema version 2) also records the missing-genotype policy as
`missing_policy`: the imputation mode (`bestguess`, or `missing_as_ref` with
`--missing-as-ref`) and the site and sample missingness thresholds. Version 1 artifacts
still load. They are read as `bestguess` with unknown thresholds, which is what they used.
With `--missing-as-ref`, the stored `AN` counts every sample.

`--mask-label NAME` stores a readable name for the mask (e.g. `pLoF`) next to its
hash. The artifact also records annotation provenance: the snpEff version
(`##SnpEffVersion`) and genome database (from `##SnpEffCmd`), the dbNSFP and gnomAD
database files used by annotation command lines (`##SnpSiftCmd` and other `##...Cmd`/
`##...Command` lines: the basename of the `-db` value, else of the first file argument,
never the input VCF), and SHA-256 hashes of the `--regions-bed` and
`--gene-file` contents (`none` when not used). A header line the input VCF does not
carry is recorded as unknown.

## Study manifest

A TSV, one row per stratum:

```
stratum_id	artifact_path
adult_cohort	/data/meta/adult.npz
pediatric_cohort	/data/meta/pediatric.npz
```

Relative paths are read relative to the manifest file. The `stratum_id` names the
cohort in the output (e.g. `per_stratum_directions` = `adult_cohort:1,pediatric_cohort:-1`).
It must not contain `:`, `;` or `,`, which separate the per-stratum output fields.

### Several masks

An optional `mask_id` column runs several masks in one `meta` call, one row per mask
and stratum (each row points at that mask's artifact for the stratum):

```
stratum_id	artifact_path	mask_id
adult_cohort	adult.plof.npz	pLoF
pediatric_cohort	pediatric.plof.npz	pLoF
adult_cohort	adult.missense.npz	missense
pediatric_cohort	pediatric.missense.npz	missense
```

Each mask needs at least two strata and is validated and combined on its own, exactly
as a single-mask run (identical results); the output stacks one row per gene and mask,
with `mask` set to the `mask_id`. The manifest `mask_id` is your own label (letters,
digits, `.`, `_`, `-`; it also names the diagnostics subdirectory), unrelated to the
artifact's internal mask id, which is a prefix of the mask hash. Two columns are added, repeated on each of a gene's
rows: `acat_masks_p`, the Cauchy (ACAT, Liu et al. 2019) combination with equal
weights of the gene's per-mask `primary_p`, and `acat_masks_q`, corrected with
`--correction` across genes (one ACAT test per gene). Only masks inside their own
correction family (`primary_q` set; not `low_mac`, `not_in_all_strata` or degenerate)
enter the combination; a gene with none has no ACAT value. Report `acat_masks_p`
rather than the best mask: picking the smallest per-mask p-value inflates type I error.
Without a `mask_id` column the output is unchanged, and `mask` is the artifacts'
`--mask-label` (or the mask hash prefix).

The meta step **refuses** (hard error, no silent fallback) when strata disagree on
genome build, annotation version, mask definition, trait type, or missing-genotype
imputation mode, or when their (hashed) sample sets overlap — overlap would violate the additive-covariance
independence assumption. Sample non-overlap does not by itself establish
independence for related or shared-control populations; that remains a caller
assumption. Strata whose missingness thresholds differ are pooled, with a warning.

It also refuses strata whose annotation provenance differs, naming the field and each
stratum's value; `--allow-annotation-mismatch` pools them anyway with a warning.
Unknown provenance (an older artifact, or a VCF without those header lines) is only
warned about, never refused.

## Running it

```bash
variantcentrifuge meta \
  --manifest study.tsv \
  --output meta_results.tsv \
  --weights beta:1,25 \
  --correction fdr \
  --min-mac 5 \
  --diagnostics-output meta_diagnostics/
```

Only `beta:a,b` and `uniform` weights are supported (the only schemes reproducible
from the artifact). `--weights` takes several schemes, e.g.
`--weights beta:1,25 uniform`: every column that depends on the weights (p-values,
`meta_burden_beta`/`_direction`, `meta_skato_rho`, heterogeneity, leave-one-out,
`primary_p`/`primary_q`, `status`, and with several masks
`acat_masks_p`/`_q`) then
appears once per scheme, suffixed by the scheme name with non-alphanumerics replaced
by `_` (`primary_p_beta_1_25`, `primary_p_uniform`). Each scheme is corrected on its
own, and its diagnostics go to `DIR/[<mask_id>/]<scheme>/`. Gene, mask, sample
counts, MAC, carrier and `warnings` columns are shared. With one scheme the column
names are unchanged. Output columns include `meta_burden_p`,
`meta_burden_direction`, `meta_skat_p`, `meta_skato_p`, `meta_skato_rho`, `het_q`,
`het_p`, `het_i2`, `per_stratum_directions`, `loo_min_p`, `primary_p`,
`primary_q`, and a `status`/`warnings` pair (e.g. a non-converged stratum is
excluded with an explicit warning).

### Leave-one-stratum-out

Each leave-one-out result is a full meta-analysis of the remaining strata (variants
re-aligned, weights re-derived from their pooled frequency), so it equals running
`meta` on the manifest without that stratum (checked against R MetaSKAT). Fixed
columns, independent of the number of strata:

- `loo_burden_by_stratum` — `stratum:p:direction;...`, the burden p-value and
  direction (`+1`/`-1`/`+0`) without each stratum, e.g. `adult_cohort:0.0123:+1`;
- `loo_skato_by_stratum` — `stratum:p;...` for SKAT-O (`NA` when not computable);
- `loo_max_burden_p` / `loo_max_burden_stratum` and `loo_max_skato_p` /
  `loo_max_skato_stratum` — the largest LOO p-value and the stratum whose removal
  gives it. A gene that is significant only while one stratum is included has a
  large `loo_max_*_p` naming that stratum (a "needs replication" flag);
- `loo_min_p` — the smallest burden LOO p-value (kept for compatibility).

With a single contributing stratum the LOO columns are empty.

### Gene sets

Gene sets from `--gene-set-file` are stored in the artifact as ordinary units under their
set id, so `meta` combines them like genes and corrects them in the same family as the
genes (one larger family, which is conservative for both).

### Minor-allele counts and `--min-mac`

`cmac_total` is the cumulative minor-allele count over the gene's union of variants,
the sum of the per-stratum counts in `cmac_by_stratum` (`stratum:mac;...`, integers),
each `Σ min(AC, AN − AC)` over that stratum's variants from the artifact's
per-variant `AC`/`AN`. Both count only variants that carry score
information (non-zero variance after the covariate projection): a site that is
monomorphic within a stratum, or explained by its covariates, adds nothing to the
test and so nothing to the MAC used by `--min-mac`.

`--min-mac N` (a finite number ≥ 0) labels genes with `cmac_total < N` as `status=low_mac`. Their p-values
are still reported, but `primary_q` is empty: like untestable genes in the main
engine, they are left out of the multiple-testing family (and of the diagnostics).
`--require-variants-in-all-strata` labels genes as `status=not_in_all_strata` unless
every stratum in the manifest contributes at least one informative minor allele
(per-stratum MAC > 0; a gene absent from a stratum, or whose null model did not
converge there, fails). Such genes are reported but handled like `low_mac` genes:
no `primary_q`, no diagnostics.

`carriers_by_stratum` (`stratum:n;...`) gives, per stratum, the number of ALT-allele
carriers: samples with at least one observed ALT allele among the gene's variants (the
allele the artifact's scores are coded on; at a site where ALT is the major allele
nearly every sample counts), after the artifact's
missingness QC (a missing call is never a carrier; with `--missing-as-ref` it is REF);
`n_carriers_total` is their sum. Artifacts written before this field existed show
`NA` for their stratum and leave the total empty. Like `cmac_*`, both cover the
strata whose null model converged.

### Diagnostics (λ_GC and QQ)

`--diagnostics-output DIR` writes the same files as a single-cohort run
(`lambda_gc.tsv`, `qq_data.tsv`, `summary.txt`, and `qq_plot.png` when matplotlib is
installed) for `meta_burden`, `meta_skat`, `meta_skato` and `primary`, per mask:
with several masks each gets its own subdirectory `DIR/<mask_id>/`. λ and QQ use exactly the genes in the correction family
(those with a `primary_q`; `low_mac` and `degenerate` genes are left out), while the
warning count in `summary.txt` covers every gene. λ from fewer than 100 genes is
flagged as unreliable.

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

The meta engine is cross-validated against R MetaSKAT 0.90 (`MetaSKAT_wZ`,
`combined.weight=TRUE`, `weights.beta=c(1,25)`): burden and SKAT agree to ~1e-15 and
SKAT-O to ~1e-4 (`tests/unit/test_r_parity_skato.py` stores the R values, so CI needs
no R). Below p≈1e-5 the per-rho tail probabilities use a saddlepoint approximation,
whereas R falls back to the Liu approximation there; Liu is anti-conservative in that
range, so very small p-values can be an order of magnitude or more larger (more
conservative) than R's.
