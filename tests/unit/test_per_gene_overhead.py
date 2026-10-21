"""Per-gene overhead cuts keep every result bit-identical to v0.25.0.

The genotype-matrix build, the score-artifact record and ACAT-V skip work that cannot
change the result (no missing calls, all samples retained, one tail call per gene); the
references below are the v0.25.0 code, verbatim where it matters.
"""

from __future__ import annotations

import gc
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from variantcentrifuge.association.backends.null_model import _project, classify_projected
from variantcentrifuge.association.genotype_matrix import build_genotype_matrix
from variantcentrifuge.association.gt_codes import GenotypeCodes
from variantcentrifuge.association.meta.align import canonical_key
from variantcentrifuge.association.meta.artifact import (
    _build_for_artifact,
    fit_stratum_null,
)
from variantcentrifuge.association.tests import acat
from variantcentrifuge.association.tests.acat import cauchy_combination, compute_acat_v
from variantcentrifuge.stages.association_stage_helpers import _score_artifact_hook
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

N = 300
# Diploid, phased, partial, haploid (0/1/.), generic missing ("", NA): every call kind.
VOCAB = ["0/0", "0/1", "1/1", "./.", "0|1", "./1", "0", "1", ".", "", "NA"]


def _ref_infer(ploidy, kind, retained):
    kept = kind[:, retained]
    haploid_rows = (kept == 1).any(axis=1) & ~(kept == 2).any(axis=1)
    ploidy[haploid_rows[:, None] & (kind == 3)] = 1.0


def _ref_build(parsed, n_samples, is_binary, site, sample, pheno, missing_as_ref):
    """v0.25.0 build_genotype_matrix after parsing (steps 1-6)."""
    raw, ploidy, kind = (a.copy() for a in parsed)
    warnings_list = []
    if missing_as_ref:
        raw[np.isnan(raw)] = 0.0
    keep = np.isnan(raw).mean(axis=1) <= site
    raw, ploidy, kind = raw[keep], ploidy[keep], kind[keep]
    if raw.shape[0] > 0:
        missing_per_sample = np.isnan(raw).mean(axis=0)
    else:
        missing_per_sample = np.zeros(n_samples, dtype=np.float64)
    sample_mask = [bool(missing_per_sample[i] <= sample) for i in range(n_samples)]
    if is_binary and pheno is not None and raw.shape[0] > 0:
        case_mask, ctrl_mask = pheno == 1, pheno == 0
        if case_mask.sum() > 0 and ctrl_mask.sum() > 0:
            miss = np.isnan(raw)
            case_rates = miss[:, case_mask].mean(axis=1)
            ctrl_rates = miss[:, ctrl_mask].mean(axis=1)
            for v in np.flatnonzero(np.abs(case_rates - ctrl_rates) > 0.05):
                warnings_list.append(
                    f"Variant index {v}: differential missingness "
                    f"(cases={case_rates[v]:.2f}, controls={ctrl_rates[v]:.2f})."
                )
    retained = np.asarray(sample_mask, dtype=bool)
    _ref_infer(ploidy, kind, retained)
    missing = np.isnan(raw)
    obs = ~missing & retained
    total_ac = np.where(obs, raw, 0.0).sum(axis=1)
    total_an = np.where(obs, ploidy, 0.0).sum(axis=1)
    zero = ~obs.any(axis=1) | (total_an <= 0.0) | (total_ac >= total_an)
    alt_freq = np.divide(total_ac, total_an, out=np.zeros_like(total_ac), where=~zero)
    mafs = np.where(zero, 0.0, np.minimum(alt_freq, 1.0 - alt_freq))
    imputed = ploidy * alt_freq[:, None]
    geno = np.where(missing, np.round(imputed) if is_binary else imputed, raw)
    flip = ~zero & (alt_freq > 0.5)
    geno[flip] = ploidy[flip] - geno[flip]
    geno[zero] = 0.0
    return geno.T, mafs, sample_mask, warnings_list, keep


def _ref_record(b, built, null):
    """v0.25.0 gene_score_record + compute_gene_score arithmetic (one null model)."""
    keep, mask = built["variant_keep_mask"], built["sample_mask"]
    raw, plo, kind = (a.copy() for a in built["parsed"])
    raw, plo = raw[keep], plo[keep]
    _ref_infer(plo, kind[keep], mask)
    if b.missing_as_ref:
        raw[np.isnan(raw)] = 0.0
    geno = np.array(built["qc_genotype_matrix"], dtype=np.float64, order="C")
    raw, plo = raw.T[mask], plo.T[mask]
    observed = ~np.isnan(raw)
    ac = np.where(observed, raw, 0.0).sum(axis=0)
    an = np.where(observed, plo, 0.0).sum(axis=0)
    flipped = (ac > an / 2.0) & (ac < an)
    geno[:, flipped] = plo[:, flipped] - geno[:, flipped]
    rows = b.gene_df.loc[keep]
    sites = (rows[c] for c in ("CHROM", "POS", "REF", "ALT"))
    keys = [canonical_key(*s) for s in zip(*sites, strict=True)]
    z_adj, z, bound = _project(geno, null)
    explained = classify_projected(z_adj, z, bound)
    u = np.where(explained, 0.0, geno.T @ null.extra["residuals"])
    z_adj[:, explained] = 0.0
    v = z_adj.T @ z_adj
    v = 0.5 * (v + v.T)
    observed = ~np.isnan(raw)
    return (
        keys,
        u,
        v,
        np.where(observed, raw, 0.0).sum(axis=0),
        np.where(observed, plo, 0.0).sum(axis=0),
        int((raw > 0).any(axis=1).sum()),
    )


def _cases(seed: int):
    """(name, codes) per gene: edge shapes over one shared vocabulary."""
    rng = np.random.default_rng(seed)
    common = rng.choice(4, size=(6, N), p=[0.9, 0.07, 0.02, 0.01])
    many_missing = common.copy()
    many_missing[:, :20] = 3  # 20 samples fail sample QC
    many_missing[2, : N // 2] = 3  # one site fails site QC
    haploid = rng.choice([6, 7, 8, 9, 10], size=(4, N), p=[0.85, 0.08, 0.03, 0.02, 0.02])
    flip = rng.choice([1, 2, 4, 5, 9], size=(3, N), p=[0.2, 0.6, 0.1, 0.05, 0.05])
    no_carrier = np.zeros((5, N), dtype=int)
    one_variant = common[:1].copy()
    no_missing = rng.choice(3, size=(5, N), p=[0.9, 0.08, 0.02])
    mixed = rng.choice(len(VOCAB), size=(8, N))
    return [
        ("common", common),
        ("many_missing", many_missing),
        ("haploid", haploid),
        ("flip", flip),
        ("no_carrier", no_carrier),
        ("one_variant", one_variant),
        ("no_missing", no_missing),
        ("mixed", mixed),
    ]


def _builders(seed: int, missing_as_ref: bool, is_binary: bool = True):
    cases = _cases(seed)
    codes = np.concatenate([c for _, c in cases]).astype(np.int8)
    gtc = GenotypeCodes(codes, VOCAB, [f"S{i}" for i in range(N)])
    df = pd.DataFrame(
        {
            "CHROM": "chr1",
            "POS": np.arange(len(codes)) * 7 + 11,
            "REF": "a",
            "ALT": "G",
        }
    )
    pheno = np.r_[np.ones(60), np.zeros(N - 60)]
    cov = np.random.default_rng(seed + 1).normal(size=(N, 2))
    builders, start = {}, 0
    for name, c in cases:
        rows = np.arange(start, start + len(c))
        builders[name] = _GenotypeMatrixBuilder(
            gene_df=df.iloc[start : start + len(c)],
            vcf_samples=[f"S{i}" for i in range(N)],
            gt_columns=[],
            is_binary=is_binary,
            missing_site_threshold=0.2,
            missing_sample_threshold=0.5,
            phenotype_vector=pheno,
            covariate_matrix=cov,
            missing_as_ref=missing_as_ref,
            genotypes=gtc,
            rows=rows,
        )
        start += len(c)
    return df, builders


def _same(a, b) -> None:
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape
    assert a.tobytes() == b.tobytes()


@pytest.mark.unit
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("missing_as_ref", [False, True])
@pytest.mark.parametrize("is_binary", [True, False])
def test_build_matches_v025(seed: int, missing_as_ref: bool, is_binary: bool) -> None:
    _df, builders = _builders(seed, missing_as_ref, is_binary)
    for name, b in builders.items():
        parsed = b.gene_arrays()
        ref = _ref_build(parsed, N, is_binary, 0.2, 0.5, b.phenotype_vector, missing_as_ref)
        got = build_genotype_matrix(
            b.gene_df,
            b.vcf_samples,
            [],
            is_binary=is_binary,
            missing_site_threshold=0.2,
            missing_sample_threshold=0.5,
            phenotype_vector=b.phenotype_vector,
            return_keep_mask=True,
            missing_as_ref=missing_as_ref,
            parsed=parsed,
        )
        _same(got[0], ref[0])
        assert got[0].flags.f_contiguous == ref[0].flags.f_contiguous, name
        _same(got[1], ref[1])
        assert got[2] == ref[2] and all(type(x) is bool for x in got[2]), name
        assert got[3] == ref[3], name
        _same(got[4], ref[4])


@pytest.mark.unit
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("missing_as_ref", [False, True])
def test_artifact_record_matches_v025(seed: int, missing_as_ref: bool) -> None:
    df, builders = _builders(seed, missing_as_ref)
    hook = _score_artifact_hook(builders, df=df)
    for name, b in builders.items():
        built = b()
        mask = built["sample_mask"]
        null = fit_stratum_null(b.phenotype_vector[mask], b.covariate_matrix[mask], "binary")
        ref = _ref_record(b, built, null)
        for record in (hook(name, b, built), hook(name, b, _build_for_artifact(b))):
            n_sites, n_kept, score = record
            assert (n_sites, n_kept) == (len(b.rows), int(built["variant_keep_mask"].sum()))
            if not ref[0]:  # every site failed QC
                assert score is None, name
                continue
            assert score.variant_keys == ref[0], name
            for got, want in zip((score.u, score.v, score.ac, score.an), ref[1:5], strict=True):
                _same(got, want)
            assert score.n_carriers == ref[5], name
    # the sample-masked path is covered (missing-as-ref keeps every sample)
    assert builders["many_missing"]()["sample_mask"].all() == missing_as_ref


@pytest.mark.unit
def test_gene_arrays_multiallelic_error_unchanged() -> None:
    """A multi-allelic code still raises naming the first offending sample, then row."""
    codes = np.zeros((3, 5), dtype=np.int8)
    codes[2, 1] = codes[1, 3] = 1
    gtc = GenotypeCodes(codes, ["0/0", "1/2"], [f"col{i}" for i in range(5)])
    df = pd.DataFrame({"CHROM": "1", "POS": [1, 2, 3], "REF": "A", "ALT": "G", "GENE": "X"})
    with pytest.raises(ValueError, match=r"'1/2' \(gene X, variant 1:3:A:G, column 'col1'\)"):
        gtc.gene_arrays(np.arange(3), df)
    ok = GenotypeCodes(codes, ["0/0", "0/1"], [f"col{i}" for i in range(5)])
    dosage, ploidy, kind = ok.gene_arrays(np.array([2, 1]), df)
    _same(dosage, ok.dosage[codes[[2, 1]]])
    _same(ploidy, ok.ploidy[codes[[2, 1]]])
    _same(kind, ok.kind[codes[[2, 1]]])


def _acat_loop(geno, residuals, trait_type, sigma2, mu_hat, weights) -> float | None:
    """v0.25.0 ACAT-V without a null model: per-variant loop."""
    score_vec = geno.T @ residuals
    if trait_type == "binary" and mu_hat is not None:
        g_ss = np.sum(geno**2 * (mu_hat * (1.0 - mu_hat))[:, np.newaxis], axis=0)
        var_vec = g_ss
    else:
        g_ss = np.sum(geno**2, axis=0)
        var_vec = sigma2 * g_ss
    p, w = [], []
    for j in range(geno.shape[1]):
        if g_ss[j] <= 1e-12 or var_vec[j] <= 0.0:
            continue
        p.append(2.0 * float(norm.sf(abs(score_vec[j] / np.sqrt(var_vec[j])))))
        w.append(float(weights[j]))
    return cauchy_combination(p, w) if p else None


@pytest.mark.unit
@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_acat_v_without_null_model_matches_loop(trait: str) -> None:
    rng = np.random.default_rng(4)
    for gene in range(60):
        p = int(rng.integers(1, 10))
        geno = rng.binomial(2, rng.uniform(0.0, 0.2, p), size=(N, p)).astype(float)
        if gene % 4 == 0:
            geno[:, 0] = 0.0  # zero variance: skipped
        residuals = rng.normal(size=N) * (50 if gene % 7 == 0 else 1)  # tiny p too
        mu = rng.uniform(0.1, 0.9, N) if trait == "binary" else None
        sigma2 = 0.0 if gene % 9 == 0 else 1.3
        w = rng.uniform(0.0, 3.0, p)
        got = compute_acat_v(geno, residuals, trait, sigma2, mu, w)
        assert got == _acat_loop(geno, residuals, trait, sigma2, mu, w), gene


@pytest.mark.unit
def test_acat_v_evaluates_normal_tail_once_per_gene() -> None:
    """One vectorised tail call per gene, not one per variant (fails on v0.25.0)."""
    rng = np.random.default_rng(2)
    geno = rng.binomial(2, 0.1, size=(N, 12)).astype(float)
    calls = 0
    real = acat.norm.sf

    def counting(x):
        nonlocal calls
        calls += 1
        return real(x)

    with patch.object(acat.norm, "sf", counting):
        compute_acat_v(geno, rng.normal(size=N), "quantitative", 1.0, None, np.ones(12))
    assert calls == 1


@pytest.mark.unit
def test_stage_matrix_builders_run_with_gc_paused(tmp_path) -> None:
    """The per-gene builder loop pauses the cyclic GC (fails on v0.25.0); it is restored."""
    from scripts.bench_association import make_df, parse_args, run_stage

    df, samples = make_df(6, 80, "packed", seed=3, miss=0.03)
    args = parse_args(["--seed", "3", "--case-frac", "0.3", "--tests", "logistic_burden"])
    states = []
    real = _GenotypeMatrixBuilder.__init__

    def recording(self, *args, **kwargs):
        states.append(gc.isenabled())
        real(self, *args, **kwargs)

    assert gc.isenabled()
    with patch.object(_GenotypeMatrixBuilder, "__init__", recording):
        run_stage(df, samples, args, tmp_path)
    assert states and not any(states)
    assert gc.isenabled()


@pytest.mark.unit
def test_code_matrix_prefaulted_after_parse_workers_fork(monkeypatch) -> None:
    """Pages are touched by threads once the workers exist (fails on v0.25.0: no _prefault);
    touching them before the fork would make the merge copy them on write."""
    from variantcentrifuge.association import engine_workers, gt_codes

    if not engine_workers.fork_available():
        pytest.skip("fork only")
    events = []
    real_map, real_prefault = gt_codes.fork_map, gt_codes._prefault

    def fork_map(*args):
        events.append("fork")
        yield from real_map(*args)

    def prefault(codes, workers):
        events.append(("prefault", codes.shape, workers))
        real_prefault(codes, workers)

    monkeypatch.setattr(gt_codes, "fork_map", fork_map)
    monkeypatch.setattr(gt_codes, "_prefault", prefault)
    df = pd.DataFrame({"GT": ["S0(0/1);S2(1/1)", "", "S1(./.)", "S3(0/1)"]})
    samples = [f"S{i}" for i in range(4)]
    got = gt_codes.parse_genotype_codes(df, samples, samples, [], workers=3)
    want = gt_codes.parse_genotype_codes(df, samples, samples, [], workers=1)
    assert events == ["fork", ("prefault", (4, 4), 3)]
    _same(got.codes, want.codes)
    assert got.vocab == want.vocab
