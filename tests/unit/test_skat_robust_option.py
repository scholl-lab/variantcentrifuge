"""--skat-robust (#223): option validation and the end-to-end Python SKAT path.

R references: ``tests/fixtures/skat_robust_parity`` (see test_skat_robust_parity.py).
Here the null model is fitted by the Python backend (statsmodels), not taken from R,
so p-values agree with R to the GLM fit tolerance (both stop at ~1e-8 deviance change).
"""

from __future__ import annotations

from unittest.mock import MagicMock

import numpy as np
import pytest

from tests.unit.test_association_primary_test import _parse, _validate
from tests.unit.test_skat_robust_parity import REF, _load
from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest


@pytest.mark.unit
class TestValidation:
    def test_flag_unset_by_default_and_mapped(self):
        from variantcentrifuge.cli_association_config import populate_association_config

        assert _parse([]).skat_robust is None
        cfg: dict = {"perform_association": True}
        populate_association_config(_parse(["--skat-robust"]), cfg, pytest.fail)
        assert cfg["skat_robust"] is True

    def test_valid_combinations(self):
        assert _validate(["--skat-robust"], ["skat"]) is None
        assert _validate(["--skat-robust", "--skat-method", "SKATO"], "fisher,skat") is None
        assert _validate([], ["skat"], {"skat_robust": True}) is None

    @pytest.mark.parametrize(
        ("argv", "tests", "fragment"),
        [
            (["--skat-robust", "--trait-type", "quantitative"], ["skat"], "binary traits only"),
            (["--skat-robust"], ["fisher"], "'skat' in --association-tests"),
            (["--skat-robust", "--skat-backend", "r"], ["skat"], "python SKAT backend"),
            (["--skat-robust", "--skat-method", "Burden"], ["skat"], "SKAT or SKATO"),
        ],
    )
    def test_invalid_combinations_fail_before_extraction(self, argv, tests, fragment):
        msg = _validate(argv, tests)
        assert msg is not None and fragment in msg

    def test_config_file_flag_is_validated_too(self):
        msg = _validate(["--trait-type", "quantitative"], ["skat"], {"skat_robust": True})
        assert msg is not None and "binary traits only" in msg

    def test_requires_perform_association(self):
        from variantcentrifuge.cli_association_config import validate_association_args

        args = _parse(["--skat-robust"])
        args.perform_association = False
        args.association_tests = None
        fail = MagicMock(side_effect=SystemExit)
        with pytest.raises(SystemExit):
            validate_association_args(args, {}, fail)
        assert "--skat-robust requires --perform-association" in fail.call_args[0][0]

    def test_config_file_key(self):
        from variantcentrifuge.stages.analysis_stage_helpers import (
            _build_assoc_config_from_context,
            _validate_association_config_dict,
        )

        ctx = MagicMock()
        ctx.config = {"skat_robust": None, "association": {"skat_robust": True}}
        assert _build_assoc_config_from_context(ctx).skat_robust is True
        ctx.config = {}
        assert _build_assoc_config_from_context(ctx).skat_robust is False
        with pytest.raises(ValueError, match="'skat_robust' must be a boolean"):
            _validate_association_config_dict({"skat_robust": "yes"})


def _run(cohort: str, gene: str, **config):
    geno, w, y, _mu, design = _load(cohort, gene)
    test = PurePythonSKATTest()
    test.check_dependencies()
    data = {
        "genotype_matrix": geno,
        "phenotype_vector": y,
        "covariate_matrix": design[:, 1:],
        "variant_mafs": geno.mean(axis=0) / 2.0,
        "proband_count": int(y.sum()),
        "control_count": int(y.size - y.sum()),
        "n_qualifying_variants": geno.shape[1],
    }
    return test.run(gene, data, AssociationConfig(variant_weights="beta:1,25", **config)), w


@pytest.mark.unit
@pytest.mark.parametrize("method", ["SKAT", "SKATO"])
@pytest.mark.parametrize(("cohort", "gene"), [("r19", "le10"), ("r81", "le10"), ("r19", "mixed")])
def test_skat_robust_end_to_end_matches_r(cohort, gene, method):
    res, _w = _run(cohort, gene, skat_method=method, skat_robust=True)
    assert res.extra["skat_p_method"] == "spa_er"
    ref = REF[(cohort, gene, method)]
    if method == "SKAT":
        assert res.p_value == pytest.approx(float(ref["p_davies_1e12"]), rel=1e-4, abs=2e-12)
    else:
        assert res.p_value == pytest.approx(float(ref["p_value"]), rel=1e-3, abs=1e-6)


@pytest.mark.unit
def test_default_path_is_unchanged_and_less_calibrated():
    # 1:9 imbalance, ultra-rare variants: the asymptotic p is smaller than the robust one.
    default, _w = _run("r19", "le10", skat_method="SKAT")
    robust, _w = _run("r19", "le10", skat_method="SKAT", skat_robust=True)
    assert default.extra["skat_p_method"] != "spa_er"
    assert default.p_value < robust.p_value


@pytest.mark.unit
def test_single_variant_gene_gets_a_robust_p_value():
    res, _w = _run("r19", "single", skat_method="SKATO", skat_robust=True)
    assert res.p_value == pytest.approx(float(REF[("r19", "single", "SKAT")]["p_value"]), rel=1e-6)
    assert res.extra["skat_o_rho"] == 0.0


@pytest.mark.unit
def test_robust_acat_v_uses_calibrated_variant_pvalues():
    """10:990, a case-only singleton and a variant in four controls (Codex review case).

    v0.24 dev: robust SKAT-O 0.0076 but ACAT-V (normal score p) 4.7e-23 drove the
    default ACAT-O primary to 9.4e-23. ACAT-V now combines the SPA/ER p-values.
    """
    from variantcentrifuge.association.engine import AssociationEngine
    from variantcentrifuge.association.spa import variant_pvalues
    from variantcentrifuge.association.tests.acat import cauchy_combination
    from variantcentrifuge.association.weights import get_weights

    n = 1000
    y = np.zeros(n)
    y[:10] = 1.0
    geno = np.zeros((n, 2))
    geno[0, 0] = 1.0
    geno[[100, 200, 300, 400], 1] = 1.0
    built = {
        "genotype_matrix": geno,
        "phenotype_vector": y,
        "covariate_matrix": None,
        "variant_mafs": geno.mean(axis=0) / 2.0,
    }
    gene = {
        "GENE": "G",
        "proband_count": 10,
        "control_count": 990,
        "n_qualifying_variants": 2,
        "_genotype_matrix_builder": lambda: dict(built),
    }

    def run(robust):
        config = AssociationConfig(skat_method="SKATO", skat_robust=robust)
        return AssociationEngine.from_names(["skat"], config).run_all([dict(gene)]).iloc[0]

    default, robust = run(False), run(True)
    assert default["acat_v_p"] < 1e-20  # default path unchanged
    _s, p, _g, _v = variant_pvalues(geno, y, np.full(n, 0.01), np.ones((n, 1)))
    w = get_weights(geno.mean(axis=0) / 2.0, "beta:1,25")
    assert robust["acat_v_p"] == pytest.approx(cauchy_combination(p, w), rel=1e-9)
    assert robust["primary_pvalue"] > 1e-3
