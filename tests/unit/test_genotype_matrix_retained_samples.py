"""#174: frequencies, imputation and recoding use only the retained (sample-QC) samples."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.genotype_matrix import build_genotype_matrix
from variantcentrifuge.association.meta.artifact import emit_stratum_artifact
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

GT = {0: "0/0", 1: "0/1", 2: "1/1"}


def _issue_frame() -> pd.DataFrame:
    """Issue #174 example: 10 retained samples, one high-missingness 1/1 carrier."""
    df = pd.DataFrame(
        {
            f"s{i}": [("0/1" if i < 2 else "0/0")] + [("0/1" if i == 0 else "0/0")] * 9
            for i in range(10)
        }
    )
    df["excluded"] = ["1/1"] + ["./."] * 9
    return df


def _builder(df: pd.DataFrame, samples: list[str], pheno: np.ndarray, site: float):
    return _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=samples,
        gt_columns=samples,
        is_binary=True,
        missing_site_threshold=site,
        missing_sample_threshold=0.8,
        phenotype_vector=pheno,
        covariate_matrix=None,
    )


@pytest.mark.unit
def test_issue_frequency_uses_retained_samples() -> None:
    """Retained samples carry 2 ALT alleles of 20 -> 0.10 (baseline reported 2/11)."""
    df = _issue_frame()
    r = _GenotypeMatrixBuilder(df, list(df), list(df), True, 0.15, 0.8, None, None)()
    assert r["genotype_matrix"].shape == (10, 10)
    assert r["variant_mafs"][0] == pytest.approx(0.10)
    assert r["genotype_matrix"][:, 0].sum() / 20 == pytest.approx(0.10)


@pytest.mark.unit
def test_issue_excluded_carriers_do_not_change_imputation() -> None:
    """A retained missing call is imputed from retained samples (1/9 -> round(2/9)=0)."""
    base = _issue_frame().drop(columns="excluded")
    base.loc[0, "s9"] = "./."
    full = base.copy()
    for i in range(8):
        full[f"excluded{i}"] = ["1/1"] + ["./."] * 9
    rf = _GenotypeMatrixBuilder(full, list(full), list(full), True, 0.8, 0.8, None, None)()
    rr = _GenotypeMatrixBuilder(base, list(base), list(base), True, 0.8, 0.8, None, None)()
    assert rf["genotype_matrix"][9, 0] == 0.0
    assert rf["variant_mafs"][0] == pytest.approx(1 / 9)
    np.testing.assert_array_equal(rf["genotype_matrix"], rr["genotype_matrix"])
    np.testing.assert_array_equal(rf["variant_mafs"], rr["variant_mafs"])


def _cohort(n_ret: int = 40, n_exc: int = 40, seed: int = 7):
    """Retained AF(v0)=15/76 with 2 missing calls; excluded samples are 1/1 at v0 only.

    All-sample AF(v0) = (15 + 2*n_exc) / (76 + 2*n_exc) > 0.5, so the baseline both
    imputed and minor-allele recoded v0 from a population the tests never see.
    """
    rng = np.random.default_rng(seed)
    ret = [f"R{i}" for i in range(n_ret)]
    exc = [f"X{i}" for i in range(n_exc)]
    v0 = ["0/0"] * n_ret
    for i in range(13):
        v0[i] = "0/1"
    v0[13] = "1/1"
    v0[38], v0[39] = "./.", "./."
    rows = {s: [v0[i]] for i, s in enumerate(ret)}
    for s in ret:
        rows[s] += [GT[int(d)] for d in rng.binomial(2, 0.1, 9)]
    for s in exc:
        rows[s] = ["1/1"] + ["./."] * 9
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"] * 10,
            "POS": list(range(100, 1100, 100)),
            "REF": ["A"] * 10,
            "ALT": ["G"] * 10,
            **rows,
        }
    )
    pheno = (np.arange(n_ret + n_exc) % 2).astype(float)
    return df, ret, exc, pheno


@pytest.mark.unit
@pytest.mark.parametrize("is_binary", [True, False])
def test_full_cohort_equals_explicit_retained_input(is_binary: bool) -> None:
    """Excluded samples cannot alter retained frequencies, imputation or recoding."""
    df, ret, exc, _ = _cohort()
    samples = ret + exc
    geno, mafs, mask, _w, keep = build_genotype_matrix(
        df,
        samples,
        samples,
        is_binary=is_binary,
        missing_site_threshold=0.6,
        return_keep_mask=True,
    )
    assert mask == [True] * len(ret) + [False] * len(exc)
    ref_geno, ref_mafs, ref_mask, _ = build_genotype_matrix(
        df, ret, ret, is_binary=is_binary, missing_site_threshold=0.6
    )
    assert all(ref_mask) and keep.all()
    np.testing.assert_array_equal(mafs, ref_mafs)
    np.testing.assert_array_equal(geno[np.asarray(mask)], ref_geno)
    # Independent check on v0: retained AF = 15/76, observed ALT-coded (not flipped).
    assert mafs[0] == pytest.approx(15 / 76)
    expected_missing = round(2 * 15 / 76) if is_binary else 2 * 15 / 76
    np.testing.assert_allclose(geno[38:40, 0], expected_missing)
    assert not np.isnan(geno).any()


@pytest.mark.unit
def test_score_artifact_matches_explicit_retained_input(tmp_path) -> None:
    """The score-artifact recoding/AC/AN use the same retained population (#174)."""
    df, ret, exc, pheno = _cohort()
    kwargs = {
        "stratum_id": "s",
        "genome_build": "GRCh38",
        "annotation_version": "v",
        "mask_id": "m",
        "mask_hash": "h",
    }
    full = emit_stratum_artifact(
        str(tmp_path / "full"), {"G": _builder(df, ret + exc, pheno, 0.6)}, **kwargs
    )
    ref = emit_stratum_artifact(
        str(tmp_path / "ref"), {"G": _builder(df, ret, pheno[: len(ret)], 0.6)}, **kwargs
    )
    a, b = full.genes["G"], ref.genes["G"]
    assert a.variant_keys == b.variant_keys
    np.testing.assert_allclose(a.u, b.u, rtol=0, atol=1e-10)
    np.testing.assert_allclose(a.v, b.v, rtol=0, atol=1e-10)
    np.testing.assert_array_equal(a.ac, b.ac)
    np.testing.assert_array_equal(a.an, b.an)
    # Independent v0 values: retained observed AC=15, AN=76.
    assert (a.ac[0], a.an[0]) == (15.0, 76.0)
