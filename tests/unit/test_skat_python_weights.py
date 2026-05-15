"""Tests for PurePythonSKATTest variant weight resolution."""

from __future__ import annotations

from types import MethodType

import numpy as np
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest
from variantcentrifuge.association.weights import beta_maf_weights, get_weights


def _contingency(geno, phenotype, **extra):
    data = {
        "genotype_matrix": geno,
        "phenotype_vector": phenotype,
        "covariate_matrix": None,
        "proband_count": int(phenotype.sum()),
        "control_count": int((phenotype == 0).sum()),
        "n_qualifying_variants": geno.shape[1],
    }
    data.update(extra)
    return data


def _make_test_with_capture(monkeypatch):
    test = PurePythonSKATTest()
    test.check_dependencies()
    captured: dict[str, np.ndarray] = {}

    def fake_test_gene(self, gene, genotype_matrix, null_model, method, weights=None, weights_beta=None):
        captured["backend_weights"] = np.asarray(weights, dtype=np.float64).copy()
        return {
            "p_value": 0.5,
            "rho": None,
            "n_variants": genotype_matrix.shape[1],
            "n_marker_test": genotype_matrix.shape[1],
            "warnings": [],
            "p_method": "test",
            "p_converged": True,
            "skip_reason": None,
        }

    test._backend.test_gene = MethodType(fake_test_gene, test._backend)

    def fake_acat_v(**kwargs):
        captured["acat_v_weights"] = np.asarray(kwargs["weights"], dtype=np.float64).copy()
        return 0.75

    monkeypatch.setattr(
        "variantcentrifuge.association.tests.skat_python.compute_acat_v",
        fake_acat_v,
    )
    return test, captured


def test_python_skat_resolves_score_column_weights_for_backend_and_acat_v(monkeypatch):
    geno = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    phenotype = np.array([1, 1, 1, 0, 0, 0], dtype=float)
    score_values = np.array([0.2, 0.5, 1.0], dtype=np.float64)
    config = AssociationConfig(
        trait_type="binary",
        variant_weights="column:nephro_candidate_score",
        variant_weight_params={"combine_with_beta": False},
    )
    test, captured = _make_test_with_capture(monkeypatch)

    result = test.run(
        "GENE1",
        _contingency(geno, phenotype, score_values=score_values),
        config,
    )

    expected = get_weights(
        geno.mean(axis=0) / 2.0,
        "column:nephro_candidate_score",
        score_values=score_values,
        weight_params={"combine_with_beta": False},
    )
    assert result.extra["acat_v_p"] == 0.75
    np.testing.assert_allclose(captured["backend_weights"], expected)
    np.testing.assert_allclose(captured["acat_v_weights"], expected)


@pytest.mark.parametrize(
    "variant_weights, extra_key, extra_values",
    [
        ("cadd", "cadd_scores", np.array([30.0, 20.0, 10.0], dtype=np.float64)),
        ("revel", "revel_scores", np.array([0.9, 0.5, 0.1], dtype=np.float64)),
        ("combined", "cadd_scores", np.array([30.0, 20.0, 10.0], dtype=np.float64)),
    ],
)
def test_python_skat_honors_functional_weight_schemes(
    monkeypatch,
    variant_weights,
    extra_key,
    extra_values,
):
    geno = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    phenotype = np.array([1, 1, 1, 0, 0, 0], dtype=float)
    config = AssociationConfig(trait_type="binary", variant_weights=variant_weights)
    test, captured = _make_test_with_capture(monkeypatch)
    contingency = _contingency(geno, phenotype, **{extra_key: extra_values})

    test.run("GENE1", contingency, config)

    expected = get_weights(
        geno.mean(axis=0) / 2.0,
        variant_weights,
        cadd_scores=contingency.get("cadd_scores"),
        revel_scores=contingency.get("revel_scores"),
    )
    np.testing.assert_allclose(captured["backend_weights"], expected)


def test_python_skat_beta_uses_post_imputation_genotype_maf(monkeypatch):
    geno = np.array(
        [
            [0.0, 1.0, 2.0],
            [0.0, 1.0, 2.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    phenotype = np.array([1, 1, 1, 0, 0, 0], dtype=float)
    pre_imputation_mafs = np.array([0.01, 0.01, 0.01], dtype=np.float64)
    config = AssociationConfig(trait_type="binary", variant_weights="beta:1,25")
    test, captured = _make_test_with_capture(monkeypatch)

    test.run(
        "GENE1",
        _contingency(geno, phenotype, variant_mafs=pre_imputation_mafs),
        config,
    )

    expected = beta_maf_weights(geno.mean(axis=0) / 2.0, a=1.0, b=25.0)
    not_expected = beta_maf_weights(pre_imputation_mafs, a=1.0, b=25.0)
    np.testing.assert_allclose(captured["backend_weights"], expected)
    assert not np.allclose(captured["backend_weights"], not_expected)
