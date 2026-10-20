# Changelog

All notable changes to VariantCentrifuge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Developer benchmark harness `scripts/bench_association.py` for the association stage (seeded synthetic S/M/L data, packed or per-sample GT layout, bit-exact `--compare` result gate, slow wall-clock regression test) and per-sub-step timings recorded in `context.stage_results["association_analysis"]`.
- `--filters` is checked before extraction: SnpSift runs once on the VCF's real header plus one synthetic record; unquoted string literals, unknown INFO fields and syntax errors now fail immediately with SnpSift's message and the expression (#211). Skipped with `--late-filtering`.
- `--temp-dir DIR` moves the workspace temp dir, bulk Python temp files, `sort -T` and SnpSift's JVM tmp dir off `/tmp`/output (passed explicitly; process `TMPDIR` is not changed); default unchanged (#213).
- Start-up free-space warning on the temp and output filesystems when free space is below 16x the compressed input VCF (advisory only; gene-panel runs need far less) (#213).

### Changed

- Association runs whose case/control or phenotyped samples are a strict subset of the VCF now extract only those samples' genotypes (`bcftools query -S`), as PLINK `--keep`, REGENIE and SAIGE do. Results and score artifacts are unchanged; the intermediate TSV and the output `GT` carrier list cover the analysed samples only. Samples are kept when inheritance, `--append-extra-sample-fields`, `--no-replacement`, `--igv`, `--pseudonymize`, a tumor/normal sample index, statistics, `--late-filtering` or `--final-filter` need them; `--no-sample-subset` opts out. The subset hash is recorded as `sample_subset_sha256` in the metadata and score-artifact provenance, edited sample/phenotype files invalidate resume checkpoints, and so does a change to `--no-sample-subset` or any of these flags (#212).
- Faster per-gene association kernels with bit-identical results: the covariate (nuisance) QR basis and its conditioning are computed once per null model instead of on every projection, and the Liu fallback p-value is computed only when Davies fails (Python SKAT/SKAT-O, ACAT-V and score-artifact U/V; ~20% less kernel time per gene for SKAT-O + ACAT-V on 5,000 samples with 13 covariates).

### Deprecated

- `--perform-gene-burden` now warns and will be removed in the next minor release; use `--perform-association` (e.g. `--association-tests fisher`). The legacy path ignores `--missing-as-ref` (#217).

### Fixed

- COAST no longer skips a gene (`ANNOTATION_GENOTYPE_MISMATCH`) when genotype site QC drops a variant; it tests the kept variants, like burden and SKAT (#215).
- Single-cohort ACAT (ACAT-O, ACAT-V, COAST omnibus) follows R `ACAT()`, which changes results: a gene where any combined test has p exactly 1 (common, e.g. Fisher p = 1.0 for a gene with equal case and control carrier rates) now gets ACAT-O ≈ 1, as in R; before, that p = 1 was silently dropped and ACAT-O was driven by the other tests. A gene with one test at p exactly 0 and another at exactly 1 gets an empty ACAT-O (R stops with an error there). p outside [0, 1] is an error. An allelic SKAT-O that COAST cannot compute (rank < 2) is now reported as undefined (`coast_skat_pvalue` empty) and left out of the COAST omnibus instead of entering as p = 1 (#216).
- Python SKAT (binary trait): a gene whose covariate-adjusted kernel has no eigenvalue above round-off is now skipped with `skat_skip_reason=no_positive_eigenvalues` (R SKAT stops there) instead of reporting a placeholder p = 1 that forced ACAT-O towards 1. Quantitative traits keep p = 1, as R computes.

## [0.22.1] - 2026-10-20

### Changed

- Test fixtures and R reference scripts use neutral synthetic scenario names; no reference values changed (#218).

## [0.22.0] - 2026-10-20

### Added

- `--annotate-bed NAME=path.bed[.gz]` adds a yes/no column `NAME` (variant REF span overlaps a region); plain `--annotate-bed path.bed` is unchanged (#94).
- `.bed.gz` input for `--annotate-bed`.
- Association diagnostics: with `--diagnostics-output` and categorical covariates, `covariate_case_control.tsv` lists cases/controls per level of each categorical covariate (raw labels); levels with no cases or no controls are warned (#202).
- New AF presets `popmax_1percent` and `popmax_rare` filter on the dbNSFP gnomAD popmax allele frequency (`gnomAD_exomes/genomes_POPMAX_AF` for dbnsfp4, `gnomAD4.1_joint_POPMAX_AF` for dbnsfp5) (#202).
- `meta --diagnostics-output DIR` writes λ_GC, QQ data and a summary for the meta burden, SKAT, SKAT-O and primary p-values, one set per run (= per mask) (#192).
- `meta` reports leave-one-stratum-out p-values named by stratum for burden (with direction) and SKAT-O, plus the largest LOO p-value and the stratum whose removal gives it (`loo_burden_by_stratum`, `loo_skato_by_stratum`, `loo_max_{burden,skato}_{p,stratum}`) (#193).
- `meta` reports cumulative minor-allele counts (`cmac_total`, `cmac_by_stratum`); `--min-mac N` marks genes below N as `status=low_mac`, and `--require-variants-in-all-strata` marks genes without an informative minor allele in every stratum as `status=not_in_all_strata`; both are left out of the multiple-testing correction and the diagnostics (#194).
- `--mask-label NAME` stores a readable mask name in the `--score-artifact`, which now also records annotation provenance (snpEff version and genome, dbNSFP and gnomAD files from the VCF header, SHA-256 of `--regions-bed`/`--gene-file`) as optional schema-v2 fields. `meta` refuses strata whose provenance differs (naming field and strata) unless `--allow-annotation-mismatch` is given; unknown provenance only warns (#197).
- Score-artifact provenance also hashes `--gene-set-file`, so `meta` refuses strata whose gene sets share ids but differ in members (unless `--allow-annotation-mismatch`), and takes dbNSFP/gnomAD sources only from the database operands of `SnpSift annotate/dbnsfp` and `bcftools annotate -a` header lines (#197, #201).
- `--score-artifact` records per-gene carrier counts (samples with ≥1 ALT allele after QC; optional schema-v2 field) and `meta` reports them as `carriers_by_stratum` and `n_carriers_total`, completing #194.
- `meta` runs several masks at once from an optional manifest `mask_id` column: each mask is combined exactly as in a single-mask run, and `acat_masks_p` (Cauchy combination of the masks' `primary_p`) plus `acat_masks_q` (corrected across genes) are added; diagnostics go to one subdirectory per mask. Without the column the output is unchanged, except that `mask` shows the `--mask-label` when the artifacts have one (#196).
- `meta --weights` accepts several schemes (e.g. `--weights beta:1,25 uniform`), giving one column set per scheme suffixed by the scheme name (`primary_p_beta_1_25`, `primary_p_uniform`), each corrected separately; diagnostics go to one subdirectory per scheme (`DIR/[<mask_id>/]<scheme>/`); a single scheme keeps today's column names (#202).
- `--gene-set-file` (TSV `set_id`, `gene`) tests each gene set as an extra unit (union of its genes' qualifying variants, shared variants once) with the same tests; set rows get `unit_type=set`, are corrected separately from genes, and appear in `--score-artifact` as ordinary units; `--resume` reruns the analysis when the set file (path or content) or `--mask-label` changes (#201).
- `--gene-set-file` (TSV `set_id`, `gene`) tests each gene set as an extra unit (union of its genes' qualifying variants, shared variants once) with the same tests; set rows get `unit_type=set`, are corrected separately from genes, and appear in `--score-artifact` as ordinary units (#201).

### Changed

- Python API: `build_stratum_artifact()` and `ScoreCovarianceArtifact` require an explicit keyword `missing_policy`; the `bestguess` default now applies only to version 1 artifacts on load, so pre-imputed matrices cannot be mislabelled (#200).
- `meta` gives zero weight to loci monomorphic in the pooled strata, so `--weights beta:a,b` with `a < 1` no longer crashes leave-one-out when a locus varies only in the dropped stratum (#193).
- `meta` leave-one-out (including `loo_min_p`) re-derives the variant weights from the remaining strata, so each LOO result equals a meta-analysis without that stratum (#193).
- **Genotype missing data (#200).** Added `--missing-as-ref`, which sets missing genotypes to
  reference in every association test (Fisher, the matrix-based tests and `--score-artifact`)
  before missingness QC. Use it for VCFs merged from single-sample VCFs, where `./.`
  means "not called as variant". Before, raising the site threshold to keep such sites imputed
  every missing sample as a carrier. `--missing-site-threshold` and
  `--missing-sample-threshold` expose the existing thresholds (defaults unchanged, range 0–1).
  Variants dropped by site-missingness QC are now logged as one WARNING with totals, and the
  results gain an `n_variants_qc` column when matrix-based tests run. The default imputation is
  unchanged, and it is now checked against R SKAT (bestguess for binary traits, fixed for
  quantitative traits).
- **Score artifact schema version 2 (#200).** `--score-artifact` records `missing_policy`
  (imputation mode and both thresholds). Version 1 artifacts still load as `bestguess` with
  unknown thresholds. `variantcentrifuge meta` refuses to pool strata with different
  imputation modes and warns when their thresholds differ. Artifacts written by this version
  cannot be read by older versions.

### Fixed

- Per-sample genotype columns built by the iterative path were misaligned after transcript de-duplication, so real carriers became missing; and positional comma-separated GT turned explicit no-calls (`./.`) into reference without any flag. Both affected association input (#200).
- `--annotate-bed` treated 1-based VCF POS as 0-based, so hits at BED interval edges were wrong (first base missed ... last base excluded); overlap now uses the REF span [POS-1, POS-1+len(REF)), so deletions reaching into a region are hit. Legacy `Region=` output changes at interval edges (#94).
- An unreadable or missing `--annotate-bed` file now raises instead of silently annotating nothing (#94).
- Association diagnostics: `IMBALANCED_COHORT` is now two-sided (max/min of cases and controls) and the warning names the direction, e.g. "cases outnumber controls 100.0:1" (#199).
- Firth logistic regression (the `logistic_burden` separation fallback) no longer fails with `FIRTH_CONVERGE_FAIL` when it is already at the optimum but no step can change the penalised log-likelihood beyond rounding noise (common with uncentred covariates such as birth year); such fits are now accepted when the Newton decrement is below 1e-8 (#198).

## [0.21.0] - 2026-10-20

Fixes the findings of the v0.19.0 audit revalidation (#173–#187).

### Fixed

- Unsplit multi-allelic genotypes (`0/2`, `2/2`, haploid `2`, `./2`, `0/1/2`) are rejected at the genotype-matrix boundary with an error naming gene, site, sample and GT, instead of being miscoded; split with `bcftools norm -m-` (#175).
- Allele frequencies, minor-allele recoding, imputation and missing-call ploidy use only the samples retained after QC, in the matrix builder and the score-artifact export; a gene left without case/control contrast is omitted from the artifact with a warning (#174).
- Quantitative SKAT residual variance uses RSS/(n − rank X), with Q and the kernel eigenvalues sharing one nuisance projection; results are invariant under y → a·y + b (#179).
- SKAT-backend Burden score variance is b′Pb, with the Bernoulli-weighted projection for binary traits (#180).
- Redundant covariates (duplicated, constant, translated) no longer crash with `Singular matrix`: a full-rank orthonormal design drives the fit, residual df and every projection; variants explained by covariates are dropped, and strongly ill-conditioned covariates give an explicit `ill_conditioned_covariates` skip (#182).
- Non-finite or constant outcomes and constant or covariate-collinear burdens return no p-value with a warning code instead of spurious or NaN p-values; the linear-burden intercept is added explicitly (#186).
- SKAT-O omnibus integral: survival form without "1 −" cancellation or truncation at 40, support-edge breakpoints, exact chi-square for rank one, adaptive outer integral with an error estimate; failures fall back to Liu and report `p_converged=False`. Davies calls are serialized for thread safety (#181).
- SKAT tail p-values below 1e-5 use Davies on the exponentially tilted mixture (`p_method="davies_tilted"`) instead of the saddlepoint (#181).
- Forced chunked processing resolves the gene column like the non-chunked loader, sorts with a locale-independent single-field key and restores input row order; results match non-chunked runs (#183).
- Checkpoint resume re-normalizes CLI sample lists (comma strings no longer become sets of characters), and restoring a no-op chunk stage no longer requires an artifact it never wrote (#184).
- Quantitative association fits the numeric `--phenotype-value-column` outcome aligned to VCF samples (also after `--remove-sample-substring`) instead of failing for lack of case/control cohorts (#187).
- `python -m variantcentrifuge.cli` propagates the exit code of `main()`; CLI tests run with `sys.executable` (#177).
- HTML report toolbar, filters, chart grid and footer filter expression fit narrow screens; the table still scrolls in its container (#185).
- Faster GT parsing: each ~1M-cell block's distinct GT strings are factorized and classified once instead of string-processing every cell, with identical output (#176; the reconstruction speed-up for #173 shipped in 0.20.1).
- Association guide: new section "Input Prerequisites: Callability and Genotype QC" (#178).

### Changed

- Unsupported `ANN[*]`/`EFF[*]` extraction fields are rejected at startup; use `ANN[0].<FIELD>` (#183).
- Chunked processing requires a rectangular TSV and fails with the line number on a malformed row (non-chunked mode skips or pads it); use `--no-chunked-processing` for such input (#183).
- `--trait-type quantitative` rejects `fisher`, `logistic_burden` and `--score-artifact` at startup; case/control definitions are not used, samples without an outcome are excluded, and `n_cases`/`n_controls` are reported as 0 (#187).
- New skip and warning codes: `ill_conditioned_covariates`, `no_informative_variants`, `CONSTANT_PHENOTYPE`, `NONFINITE_PHENOTYPE`, `NONFINITE_COVARIATES`, `BURDEN_NOT_ESTIMABLE`, `ILL_CONDITIONED_COVARIATES`, `NONFINITE_STATISTICS`, `OUTCOME_EXPLAINED_BY_COVARIATES`, `INSUFFICIENT_RESIDUAL_DF` (see the association guide's troubleshooting section).
- In the deep tail and some edge cases, SKAT/SKAT-O p-values deliberately differ from R SKAT 2.2.5 where R's Davies accuracy or `integrate()` is the less accurate side; moderate p-values are unchanged.
- Dependencies: smart-open 8 (`smart_open.open` replaces the removed `smart_open.smart_open`), ruff 0.17, sphinx>=9.1, github/codeql-action 4.39.

### Added

- Permanent stage-level invariance test for repeated identical transcript annotations, for binary and quantitative traits (#188).

## [0.20.1] - 2026-10-20

### Changed

- Faster association and variant analysis on wide cohorts (identical results): genotype-matrix parsing runs once per gene block instead of once per sample column, GT-column reconstruction normalises all cells at once and drops per-sample columns in one step (deleting them one by one was quadratic), and case/control counting no longer parses an implicit `0/0` for every absent sample. Panels under 50 genes run association in-process instead of starting a worker pool. One 1000 Genomes cohort (3,391 variants x 1,039 samples): 21.5 s -> 6.2 s. Peak memory stays bounded: work runs in blocks of about 1M cells.

### Fixed

- BED chunks for parallel extraction were written to a `chunks/` directory next to the input BED (for `--regions-bed`, the user's input directory). They now go to the run's workspace, so read-only input directories work and concurrent runs on one BED cannot overwrite each other's chunks.
- `make test-e2e` (and the Docker CI e2e step) now tests the package installed in the image instead of the mounted, uncompiled checkout, so SKAT runs with the compiled Davies extension as shipped rather than silently on the Liu/saddlepoint fallback. A guard test fails if that regresses.

## [0.20.0] - 2026-10-20

Python 3.12+ and a faster CI.

### Changed

- **Breaking:** requires Python 3.12+ (was 3.10+). Python 3.10 is end-of-life, pandas 3 needs 3.11+, and the Docker/conda builds already ship 3.12. Users on 3.10/3.11 keep getting the last 0.19.x release from pip.
- CI: PRs test Python 3.12 only; `main` and releases test 3.12 and 3.14. Tests run in parallel (`-n auto --dist worksteal`), coverage runs once (in the Coverage workflow, with the `sysmon` tracer), feature branches no longer run CI twice (push + PR), and superseded PR runs are cancelled. Actions moved to Node 24 releases, every job has a timeout, and Dependabot now updates actions and `uv.lock` weekly.
- Tests: the Excel hyperlink-cap regression tests (#159) reach xlsxwriter's 65,530-link cap by pre-filling its counter instead of building 10,000-row workbooks; the fast suite drops from about 60 s to 13 s locally.

## [0.19.2] - 2026-10-20

- Re-release of 0.19.1 with version metadata synchronised across `version.py`, `CITATION.cff` and `.zenodo.json`. The v0.19.1 tag failed release validation because only `version.py` had been bumped.

## [0.19.1] - 2026-10-20

### Fixed

- Snakemake workflow environment pins `smart-open<8` like `pyproject.toml`; 8.x removed `smart_open.smart_open` (#189).

## [0.19.0] - 2026-10-20

Cohort-stratified meta-analysis (#112, #114).

### Added

- Cohort-stratified rare-variant score-statistic meta-analysis (#112): `--score-artifact` writes a cohort's per-gene score/covariance artifact during association analysis, and the new `variantcentrifuge meta` subcommand combines two or more cohorts (burden, SKAT, SKAT-O, Cochran Q/I², leave-one-out, exome-wide FDR). Binary traits only.
- Real-data end-to-end test: two 1000 Genomes cohorts (EUR, SAS) for 10 nephrology genes, annotated like the AGDE cohort (snpEff 5.2 GRCh38.p14, ClinVar), run through association, `--score-artifact` and `variantcentrifuge meta` inside the Docker image (`make test-e2e`, also in the Docker CI workflow). `make fixture-1000g` rebuilds the fixture byte-identically from public sources (#60).

### Fixed

- SKAT-O (and the COAST allelic SKAT-O component) p-values: the omnibus integral used one 128-node Gauss-Legendre rule over the chi-square(1) density, whose x^-1/2 singularity at 0 and rho-switching kinks made p-values 1–8% too large on typical genes and up to ~2x too large near p=0.01. The integral is now split at the kinks and taken in t=sqrt(x), matching R SKAT 2.2.5 to ~1e-4 at the same speed. Burden and SKAT p-values were not affected.
- Burden/association with `--split-snpeff-lines`: a variant annotated on several transcripts of the same gene entered the regression genotype matrix (logistic/linear burden, SKAT, SKAT-O, COAST) once per transcript and inflated `n_qualifying_variants`. Rows are now collapsed to one per variant per gene, keeping snpEff's most severe annotation. Fisher carrier/allele counts were already per-sample and unaffected.
- `variantcentrifuge meta` crashed (ZeroDivisionError) on genes with a single qualifying variant; SKAT-O now returns no value there, as in single-cohort analysis, and the primary p-value falls back to burden.
- Provenance recorded SnpSift 5.2 as `N/A`: `SnpSift annotate -version` prints only the number.

## [0.18.0] - 2026-10-20

Remediates all 45 M1–M6 audit issues (#127–#171) and the follow-up six-part review (#172).

### Added

- Optional subcommands (`run`, `report`, `validate-config`, `build-pm5`, `download-references`) alongside the legacy CLI (#168).
- `--from-tsv` mode to rescore, recompute inheritance and regenerate reports without VCF extraction (#169).
- HTML report variant pinning, reviewer notes, "Pinned Only" filter and TSV/JSON shortlist export (#170).
- Analysis profiles (`--profile`, `--save-profile`, `--list-profiles`) and a validated `--sample-sheet` adapter (#171).
- Explicit `|type:numeric` / `|type:categorical` scoring variable types with a fixed type order; ambiguous untyped variables now fail with a clear error instead of being guessed.
- Coverage workflow and README badge; `make build-ext` builds the Davies C extension, and `make test*` run in parallel.

### Fixed

- Security: sandboxed stats expressions, HTML autoescaping and URL checks, Excel formula-injection off, private PM5 downloads with SHA-256 checks and timeouts, pseudonym map and sample IDs kept out of `run_summary.json`, HTML, metadata and `Inheritance_Details` under `--pseudonymize` (#127–#131).
- Docker: Jackson 2.22.4 in the rebuilt snpEff/SnpSift jars (CVE-2026-89407, -89425, -68497, -91776, -91777, -19032, -83557).
- Statistics: missing genotypes excluded from burden/Fisher denominators, haploid dosages, ALT-major dosage flip, covariate-adjusted ACAT-V, Haldane-corrected Fisher OR/CI, Firth step-halving, Davies fallbacks, cases removed from controls when both definitions match (#137–#145).
- Inheritance: X PAR1/PAR2 treated as autosomal using the configured GRCh37/GRCh38 build, PLINK chr26 as MT, stricter de novo calls, symmetric compound-het pairing, scalar/vectorized X-linked parity (#139).
- Pipeline: final filtering kept in the main process, inheritance computed for large in-memory data, chunked output fails on schema mismatch instead of misaligning columns, stage soft dependencies ordered, checkpoint hardening, stage failures propagated (#147–#153, #157).
- Wheels ship the HTML report templates and assets; shipped scoring configs reduce multi-value cells (`reduce:max`); the COAST custom scoring template loads.

### Changed

- Production modules split to ≤ 600 lines each with an empty `.loc-allowlist` (#154, #157); strict mypy, Python 3.10–3.13 CI matrix, `uv.lock` and PyPI Trusted Publishing (#164–#166).

## [0.17.16] - 2026-10-16

### Fixed

- Parse `snpEff -version` and `SnpSift annotate -version` from combined stdout/stderr and record active `Filter`, `Preset`, `Final Filter`, `Tumor Samples`, `Normal Samples`, `Git SHA`, and `Transcript File` (with MD5) in `*_metadata.tsv` (#123).

## [0.17.15] - 2026-10-16

### Fixed

- Use explicit minute time units (`m`) for `runtime` in Snakemake profile `default-resources` and `set-resources` so `humanfriendly` does not interpret bare integers as seconds (#125).

## [0.17.14] - 2026-10-15

### Added

- Support position-specific `FILTER` exceptions (`--filter-exception-file`) for somatic hotspot rescue while preserving original filters in `INFO/RAW_FILTER` (#121).
- Add `--fail-on-empty` to exit with code 2 when zero variants remain after filtering (#119).

### Fixed

- Resolve tumor and normal sample indices from `##tumor_sample`/`##normal_sample` VCF headers, unify `mutect2_*` and `somatic*` presets, and fail fast when `bcftools query -l` fails (#117).
- Run `TranscriptFilterStage` before `SnpSiftFilterStage` in sequential and parallel modes, match base transcript IDs across RefSeq version skew, emit all kept transcripts per gene, and include `FEATUREID` in merge keys (#118).
- Log per-stage input/output variant counts and accept valid 1-line header-only parallel chunk TSVs (#119).
- Preserve `config.json` options when CLI flags are omitted, forward P0 flags in Snakemake helpers, and append `extra_sample_fields` (retaining `0/0` normal samples) in reconstructed `GT` columns (#120).
- Decompose multiallelic VCF records via `bcftools norm -m -any` before allele-indexed filtering and restrict `ANN`/`LOF`/`NMD` annotations to the decomposed `ALT` allele (#122).

## [0.17.13] - 2026-07-19

### Changed

- Add shared multi-LLM contributor guidance and a minimal Claude Code entrypoint.
- Enforce the 600-line production Python module budget locally and in CI, with
  a fixed baseline for existing oversized modules and a required refactor when
  touching them.

## [0.17.12] - 2026-07-19

### Security

- Refresh the digest-pinned container base, rebuild the pinned SnpEff and SnpSift 5.2 JARs with
  fixed Java runtime dependencies, and verify compatibility with deterministic behavioral smoke
  tests. The container workflow now gates publication on an all-severity actionable Trivy scan,
  retains the complete fixed and vendor-unfixed JSON audit, and signs the gated image by immutable
  digest after publication.

## [0.17.11] - 2026-05-16

### Fixed
- Restore checkpoint resume state for chunked and downstream analysis stages so
  continuation runs recompute the necessary DataFrame context before
  association and output stages.
- Clear persisted checkpoint completion for recomputed stages, including the
  `Stage.__call__` checkpoint skip path.
- Preserve valid empty outputs by writing header-only requested files.
- Normalize categorical genotype values before association genotype matrix
  encoding.

## [0.17.10] - 2026-05-16

### Fixed
- Expose Fisher-only weighted FDR results under `fisher_qvalue`,
  `fisher_weighted_qvalue`, and `primary_qvalue` instead of ACAT-O pass-through
  column names.
- Count association-stage significance summaries from `primary_qvalue` so
  weighted and unweighted q-value aliases do not affect the reported number of
  significant genes.

## [0.17.9] - 2026-05-16

### Fixed
- Honor `--split-snpeff-lines before_filters` in parallel chunk processing so SnpSift consequence predicates are evaluated on split SnpEff rows before transcript selection.
- Prune SnpEff `LOF` and `NMD` INFO entries on ANN split rows to the retained row's matching gene symbol or gene ID, preventing cross-gene burden and association leakage.
- Warn and document when transcript-filtered gene burden or association uses record-level `ANN[ANY]`, `LOF[*]`, or `NMD[*]` predicates without `before_filters`.

## [0.17.8] - 2026-05-16

### Fixed
- Correct malformed built-in SnpSift presets that could leak common or benign variants when combined with other filters.
- Fail closed on SnpSift parser/evaluation diagnostics even when SnpSift exits successfully.
- Fail closed on invalid final DataFrame and TSV filter expressions instead of returning unfiltered results.
- Preserve stderr on command failures for downstream diagnostics.

## [0.17.7] - 2026-05-16

### Fixed
- Parse SnpEff `NMD[0].PERC` from the documented percent-affected field instead of returning the NMD gene symbol.

## [0.17.5] - 2026-05-16

### Fixed
- Order phenotype integration after threaded complete processing so auto-threaded runs do not start before field extraction has completed.
- Allow chunked analysis scheduling when optional genotype replacement is omitted, including `--no-replacement` runs.
- Handle categorical per-sample GT columns in column-based gene burden aggregation.
- Write gene burden sidecar output next to the requested final output using the expected `.gene_burden.tsv` name.
- Keep temporary scoring formula variables out of user-facing integration output.

## [0.17.4] - 2026-05-15

### Fixed
- Build `snpEff genes2bed` commands with expansion options before genome and gene positional arguments.
- Capture and propagate stdout/stderr diagnostics from gene BED subprocess failures.
- Allow `--regions-bed` with `--gene-name all` to use the supplied BED directly.

## [0.17.3] - 2026-05-15

### Fixed
- Restore and harden custom BED and gene-list annotation config handling so CLI and config-driven runs preserve `Custom_Annotation` behavior across legacy and canonical keys.
- Preserve checkpoint resume compatibility across annotation config key migration while still invalidating checkpoints when annotation inputs change.

## [0.17.2] - 2026-05-14

### Fixed
- Prevent raw per-sample `GEN[N].GT` columns from leaking into final TSV and Excel outputs by default after variant analysis re-attaches them for downstream burden/association stages.
- Preserve explicit `--no-replacement` behavior so users can still request raw per-sample genotype columns in final outputs.

## [0.17.1] - 2026-05-13

### Fixed
- Apply `--transcript-list` and `--transcript-file` to extracted SnpEff `ANN` transcript annotations, including MANE transcript filtering.
- Preserve multiallelic split ordering across transcript filtering and SnpSift filtering pipeline stages.
- Recompress transcript-filtered VCF outputs with `bcftools` when indexing is requested.
- Restore the Docker scan workflow by using a valid Trivy action tag and retrying transient GHCR logins.

## [0.15.0] - 2026-02-22

### Added
- **Modular association testing framework** with pluggable test architecture:
  - Fisher's exact test with carrier/allele modes, odds ratio, and 95% CI
  - Logistic burden regression for binary traits with Firth fallback for separation
  - Linear burden regression for quantitative traits (OLS)
  - SKAT and SKAT-O with pure Python backend (numpy/scipy, thread-safe); R backend deprecated
  - COAST allelic series test (BMV/DMV/PTV categories); pure Python backend default
  - ACAT-O per-gene omnibus via Cauchy combination; ACAT-V per-variant score within SKAT
- **Covariate system**: TSV/CSV covariate files with auto-detected delimiter and categorical encoding
- **PCA integration**: PLINK `.eigenvec`, AKT output, and generic TSV formats; optional AKT subprocess
- **Variant weights**: Beta(MAF), uniform, CADD, REVEL, and combined functional weight schemes
- **Diagnostics**: Genomic inflation factor (lambda_GC), QQ data TSV, optional matplotlib QQ plot
- **JSON config**: `"association"` section in config.json with validation and CLI override precedence
- **Single FDR strategy**: Benjamini-Hochberg correction on ACAT-O p-values only (not per-test)
- **Association testing guide**: Comprehensive user documentation with examples and troubleshooting

## [0.14.0] - 2026-02-18

### Added
- **Modern HTML report** — complete UX overhaul of the individual variant report across 5 phases:
  - **JS Stack Modernization**: DataTables v2, Chart.js (65KB vs Plotly 3.5MB, 98% reduction), Tippy.js tooltips, all assets vendored for offline single-file reports
  - **Summary dashboard**: metric cards (total variants, genes, samples, impact breakdown, top genes), impact distribution and inheritance pattern charts
  - **Semantic color badges**: IMPACT (red/orange/amber/gray), ClinVar (red-to-green severity), inheritance patterns (de novo=red, compound het=purple, AD=blue, AR=green, X-linked=teal)
  - **Table redesign**: sticky GENE column (FixedColumns), dark header, expandable row detail panels, content density toggle (Compact/Regular/Relaxed with localStorage), intelligent column widths, zebra striping
  - **Column-level filtering**: noUiSlider range sliders for numeric columns (POS, gnomAD AF, CADD), categorical dropdowns (IMPACT, ClinVar, Inheritance), text search (GENE), removable filter chips, "Include missing values" toggle, reactive chart updates
  - **Unified toolbar**: all controls in one 28px row — entries/page, filters, missing, search, density, show/hide columns, PDF export
  - **Accessibility (WCAG 2.1 AA)**: skip-link, ARIA roles/labels, keyboard-accessible tooltips, SVG icons with screen-reader text, chart data table fallbacks, 4.5:1 contrast ratios on all badges
  - **Print/PDF support**: @media print stylesheet hiding interactive controls, PDF export via browser print dialog
  - Report metadata footer with filter criteria, VCF source, reference genome, version, and date
  - Expanded summary.json with inheritance distribution, top genes, and sample count
  - Loading skeleton with shimmer animation during DataTable initialization
- 107 HTML report-specific tests (structure, behavior, assets, accessibility, print)
- Unified resource auto-detection across pipeline modes

## [0.13.1] - 2026-02-16

### Fixed
- Resource auto-detection across pipeline modes

## [0.13.0] - 2026-02-16

### Added
- **Performance optimization** across 7 phases:
  - Benchmark framework with timing instrumentation
  - Vectorized genotype replacement (Pandas-native operations)
  - DataFrame optimization with sanitized column names for itertuples compatibility
  - Inheritance analysis optimization with vectorized deduction
  - Output stage optimization
  - Pipeline I/O elimination (reduced intermediate file writes)
  - Parallelization and chunking with memory-aware chunk sizing
- Memory management system with SLURM/PBS/cgroup detection
- Golden file infrastructure for inheritance validation
- Performance benchmark test suite

## [0.12.0] - 2025-12-01

### Added
- **Stage-based pipeline architecture** (`pipeline_core/`) with modular stages, dependency graph, topological sort, and parallel execution via ThreadPoolExecutor/ProcessPoolExecutor
- **Restructured Snakemake workflow** matching lab pipeline conventions (Issue #68):
  - Snakemake 8+ with `min_version("8.0")` and executor plugins
  - Schema-validated config (`config/config_vc.yaml`) and sample sheet (`config/samples.tsv`)
  - Profile layering: `profiles/default/` (resources) + `profiles/{bih,charite,local}/` (executor)
  - Auto-detecting launcher script (`scripts/run_snakemake.sh`) for BIH, Charite, and local
  - Singularity/Apptainer container support via `container:` directive with conda fallback
- **Docker image** on GHCR (`ghcr.io/scholl-lab/variantcentrifuge`) with all bioinformatics tools pre-installed:
  - Multi-stage build using micromamba for minimal image size
  - CI/CD pipeline with automated builds, Trivy security scanning, and cosign image signing
  - `docker-compose.yml` with volume mount patterns for data, snpEff databases, and custom scoring configs
  - Non-root container for security best practices
- **Field profiles** for annotation database version compatibility (`--field-profile`):
  - Built-in `dbnsfp4` and `dbnsfp5` profiles for seamless gnomAD field switching
  - Template syntax `{{fragment:param}}` for parameterized filter presets
  - `--list-field-profiles` to show available profiles
  - Custom profiles configurable in `config.json` without code changes
- **Inheritance analysis** with three-pass pipeline (deduction, compound het, prioritization):
  - Supported patterns: de novo, AD, AR, X-linked (XLR/XLD), compound heterozygous, mitochondrial
  - PED file integration for family-based analysis (`--ped`, `--inheritance-mode`)
  - Vectorized compound het detection (10-50x faster than original)
  - Segregation analysis with Fisher's exact test
  - Pattern prioritization by clinical significance
- **bcftools pre-filtering** (`--bcftools-prefilter`) for early variant filtering during extraction
- **Final filtering** (`--final-filter`) using pandas query syntax on any column including computed scores
- **Sample pseudonymization** for privacy-preserving data sharing (Issue #34):
  - Multiple naming schemas: sequential, categorical, anonymous, and custom patterns
  - Consistent pseudonym mapping across all output formats (TSV, Excel, HTML)
  - PED file pseudonymization support (`--pseudonymize-ped`)
  - Secure mapping table storage in parent directory
- **Checkpoint and resume system** for robust pipeline execution:
  - Automatic pipeline state tracking with `.variantcentrifuge_state.json`
  - Resume capability after interruptions (`--enable-checkpoint` and `--resume`)
  - Optional file checksum validation (`--checkpoint-checksum`)
  - Interactive resume point selection (`--interactive-resume`)
  - Thread-safe state updates for parallel chunk processing
- Unified annotation system supporting BED files, gene lists, and JSON gene data
- JSON gene annotation feature with flexible field mapping (`--annotate-json-genes` and `--json-gene-mapping`)
- Comprehensive Sphinx documentation with modern Furo theme
- GitHub Actions workflow for automated documentation deployment
- Tumor-normal filtering presets (`somatic`, `loh`, `tumor_only`) with configurable sample indices and thresholds
- VCF annotation inspection (`--show-vcf-annotations`) for field discovery
- Genotype filtering (`--genotype-filter`) with per-gene override support
- Transcript-level filtering (`--transcript-list`, `--transcript-file`)
- ClinVar PM5 annotation support (`--clinvar-pm5-lookup`)

### Changed
- Documentation migrated from README to structured Sphinx documentation
- Enhanced filtering with three-stage approach (bcftools pre-filter, SnpSift filter, final filter)

### Fixed
- Numeric type conversion in final filtering to handle mixed data types correctly
- **Gene burden analysis edge cases** (Issue #31): Improved handling of infinite and zero odds ratios

## [0.5.0] - 2024-08-01

### Added
- Interactive HTML report generation with sortable tables
- IGV.js integration for genomic visualization
- Cohort analysis and reporting functionality
- Gene burden analysis with Fisher's exact test
- Phenotype data integration and filtering
- Preset filter system for common analysis workflows
- Excel output format support
- External database links (SpliceAI, Franklin, Varsome, gnomAD, ClinVar)

### Changed
- Improved command-line interface with better argument organization
- Enhanced configuration system with JSON-based presets
- Optimized variant filtering pipeline
- Better error handling and user feedback

### Fixed
- VCF header normalization for indexed fields
- Sample identification in cohort reports
- IGV report generation with proper FASTA handling

## [0.4.0] - 2024-06-01

### Added
- Comprehensive test suite with pytest
- Pre-commit hooks for code quality
- Gene list annotation functionality
- Variant statistics and metadata generation

### Changed
- Modular code architecture with clear separation of concerns
- Improved logging and debugging capabilities
- Enhanced VCF processing pipeline

### Fixed
- Gene BED file generation edge cases
- Genotype replacement functionality

## [0.3.0] - 2024-04-01

### Added
- Gene-centric variant filtering
- SnpSift integration for field extraction
- Basic HTML report generation
- Phenotype integration capabilities

### Changed
- Migrated from Bash/R to Python-based pipeline
- Improved error handling and validation

## [0.2.0] - 2024-02-01

### Added
- Initial Python CLI implementation
- VCF processing with bcftools integration
- Basic filtering capabilities
- Configuration file support

### Changed
- Complete rewrite from shell scripts to Python

## [0.1.0] - 2024-01-01

### Added
- Initial release
- Basic variant filtering functionality
- Shell script-based pipeline
- Simple output generation

---

## Legend

- **Added** for new features
- **Changed** for changes in existing functionality
- **Deprecated** for soon-to-be removed features
- **Removed** for now removed features
- **Fixed** for any bug fixes
- **Security** for vulnerability fixes
