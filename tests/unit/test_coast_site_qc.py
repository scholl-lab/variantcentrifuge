"""COAST tests the variants kept by genotype site QC instead of skipping the gene (#215)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.allelic_series import COASTTest
from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

N = 200
SAMPLES = [f"S{i}" for i in range(N)]
CONFIG = AssociationConfig(trait_type="binary", coast_weights=[1.0, 2.0, 3.0])


def _gene_df(sparse_row: int | None) -> pd.DataFrame:
    """Synthetic 6-variant gene (2 BMV, 2 DMV, 2 PTV); ``sparse_row`` is ~50% uncalled."""
    rng = np.random.default_rng(7)
    rows = []
    for j in range(6):
        gts = np.where(rng.random(N) < 0.15, "0/1", "0/0").astype(object)
        if j == sparse_row:
            gts[::2] = "./."
        rows.append(dict(zip(SAMPLES, gts, strict=True)))
    df = pd.DataFrame(rows)
    df.insert(0, "EFFECT", ["missense_variant"] * 4 + ["frameshift_variant"] * 2)
    df.insert(1, "IMPACT", ["MODERATE"] * 4 + ["HIGH"] * 2)
    df.insert(2, "dbNSFP_SIFT_pred", ["tolerated"] * 2 + ["deleterious"] * 2 + ["."] * 2)
    polyphen = ["benign"] * 2 + ["probably_damaging"] * 2 + ["."] * 2
    df.insert(3, "dbNSFP_Polyphen2_HDIV_pred", polyphen)
    return df


def _gene_data(df: pd.DataFrame) -> dict:
    pheno = (np.arange(N) % 2).astype(float)
    built = _GenotypeMatrixBuilder(df, SAMPLES, SAMPLES, True, 0.10, 0.80, pheno, None)()
    return {
        **built,
        "gene_df": df.reset_index(drop=True),
        "proband_count": 100,
        "control_count": 100,
        "n_qualifying_variants": len(df),
    }


def _python_coast(data: dict):
    test = PurePythonCOASTTest()
    test.check_dependencies()
    return test.run("G", data, CONFIG)


def test_python_coast_tests_variants_kept_by_site_qc():
    data = _gene_data(_gene_df(sparse_row=0))
    assert data["variant_keep_mask"].tolist() == [False] + [True] * 5
    result = _python_coast(data)
    assert result.p_value is not None, result.extra
    # Identical to running COAST on the gene without the dropped site.
    expected = _python_coast(_gene_data(_gene_df(sparse_row=0).iloc[1:]))
    assert result.p_value == expected.p_value
    assert (result.extra["coast_n_bmv"], result.extra["coast_n_dmv"]) == (1, 2)
    assert result.extra["coast_n_ptv"] == 2


def test_mismatch_guard_kept_without_mask():
    data = _gene_data(_gene_df(sparse_row=0))
    del data["variant_keep_mask"]
    assert _python_coast(data).extra["coast_skip_reason"] == "ANNOTATION_GENOTYPE_MISMATCH"


def test_mac_filtered_gene_is_skipped_as_empty_matrix():
    df = _gene_df(sparse_row=None)
    df.loc[:, SAMPLES] = "0/0"
    data = _gene_data(df)
    assert data["mac_filtered"] and data["variant_keep_mask"].all()
    assert _python_coast(data).extra["coast_skip_reason"] == "NO_PHENOTYPE_OR_EMPTY_MATRIX"


def test_r_coast_receives_kept_variants(monkeypatch):
    test = COASTTest()
    seen = {}

    def fake_run(*, geno_filtered, anno_codes_filtered, **_):
        seen["geno"], seen["anno"] = geno_filtered, anno_codes_filtered
        return 0.5, 0.5, 0.5

    monkeypatch.setattr(test, "_run_r_coast", fake_run)
    data = _gene_data(_gene_df(sparse_row=2))  # drop the first DMV
    result = test.run("G", data, CONFIG)
    assert result.p_value == pytest.approx(0.5)
    assert seen["anno"].tolist() == [1, 1, 2, 3, 3]
    assert seen["geno"].shape[1] == 5
