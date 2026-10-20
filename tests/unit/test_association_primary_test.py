"""--primary-test and --acat-o-exclude (#210): which result drives primary_* and FDR."""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig, TestResult
from variantcentrifuge.association.engine import AssociationEngine
from variantcentrifuge.association.tests.acat import compute_acat_o as cauchy

GENES = ["G1", "G2", "G3"]
SKAT_P = {"G1": 0.2, "G2": 0.04, "G3": 0.5}
ACAT_V_P = {"G1": 1e-5, "G2": 0.3, "G3": 0.6}
BURDEN_P = {"G1": 0.08, "G2": 0.01, "G3": None}


class _FakeTest:
    """Minimal AssociationTest stand-in returning fixed p-values per gene."""

    def __init__(self, name: str, pvals: dict[str, float | None], extra: dict | None = None):
        self.name = name
        self._p = pvals
        self._extra = extra or {}

    def effect_column_names(self) -> dict[str, Any]:
        return {"effect": None, "se": None, "ci_lower": None, "ci_upper": None}

    def check_dependencies(self) -> None:
        pass

    def prepare(self, n: int) -> None:
        pass

    def finalize(self) -> None:
        pass

    def run(self, gene: str, data: dict, config: Any) -> TestResult:
        extra = {k: (v[gene] if isinstance(v, dict) else v) for k, v in self._extra.items()}
        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=self._p[gene],
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=10,
            n_controls=10,
            n_variants=2,
            extra=extra,
        )


def _skat(method: str = "SKATO") -> _FakeTest:
    return _FakeTest("skat", SKAT_P, {"skat_method": method, "acat_v_p": ACAT_V_P})


def _burden(name: str = "logistic_burden") -> _FakeTest:
    return _FakeTest(name, BURDEN_P)


def _run(tests: list[_FakeTest], **cfg: Any):
    config = AssociationConfig(**cfg)
    data = [{"GENE": g, "n_qualifying_variants": 2} for g in GENES]
    return AssociationEngine(tests, config).run_all(data).set_index("gene")


@pytest.mark.unit
class TestPrimaryTestDefault:
    def test_default_is_acat_o_with_acat_v(self):
        df = _run([_burden(), _skat()])
        assert set(df["primary_test"]) == {"acat_o"}
        exp = cauchy({"logistic_burden": 0.08, "skat": 0.2, "acat_v": 1e-5})
        assert df.loc["G1", "primary_pvalue"] == pytest.approx(exp, rel=1e-12)
        assert "acat_o_qvalue" in df.columns

    def test_explicit_acat_o_equals_default(self):
        a = _run([_burden(), _skat()])
        b = _run([_burden(), _skat()], primary_test="acat_o")
        assert a.equals(b)


@pytest.mark.unit
class TestPrimarySkato:
    def test_skato_drives_primary_and_fdr(self):
        df = _run([_burden(), _skat()], primary_test="skato", skat_method="SKATO")
        assert set(df["primary_test"]) == {"skat_o"}
        assert df["primary_pvalue"].to_dict() == SKAT_P
        assert df["primary_qvalue"].tolist() == df["skat_o_qvalue"].tolist()
        # BH on 3 SKAT-O p-values: 0.04 -> 0.12, 0.2 -> 0.3, 0.5 -> 0.5
        assert df["primary_qvalue"].to_dict() == pytest.approx({"G1": 0.3, "G2": 0.12, "G3": 0.5})
        # ACAT-O is still computed and reported, but is not FDR-corrected.
        assert "acat_o_pvalue" in df.columns
        assert "acat_o_qvalue" not in df.columns

    def test_skato_requires_skato_method(self):
        with pytest.raises(ValueError, match="--skat-method SKATO"):
            AssociationEngine([_skat("SKAT")], AssociationConfig(primary_test="skato"))

    def test_skato_requires_skat_test(self):
        with pytest.raises(ValueError, match="skat"):
            AssociationEngine(
                [_burden()], AssociationConfig(primary_test="skato", skat_method="SKATO")
            )

    def test_skato_alone(self):
        df = _run([_skat()], primary_test="skato", skat_method="SKATO")
        assert set(df["primary_test"]) == {"skat_o"}
        assert df["primary_pvalue"].to_dict() == SKAT_P


@pytest.mark.unit
class TestPrimaryBurden:
    def test_binary_burden_is_logistic(self):
        df = _run([_burden(), _skat()], primary_test="burden")
        assert set(df["primary_test"]) == {"logistic_burden"}
        assert df.loc["G1", "primary_pvalue"] == 0.08
        # Gene with no burden p is excluded from the FDR family (BH on 2 values).
        q = df["logistic_burden_qvalue"]
        assert q[["G1", "G2"]].tolist() == pytest.approx([0.08, 0.02])
        assert q.isna()["G3"]

    def test_quantitative_burden_is_linear(self):
        df = _run(
            [_burden("linear_burden"), _skat()], primary_test="burden", trait_type="quantitative"
        )
        assert set(df["primary_test"]) == {"linear_burden"}

    def test_burden_missing_raises(self):
        with pytest.raises(ValueError, match="logistic_burden"):
            AssociationEngine([_skat()], AssociationConfig(primary_test="burden"))

    def test_unknown_primary_raises(self):
        with pytest.raises(ValueError, match="primary_test"):
            AssociationEngine([_skat()], AssociationConfig(primary_test="fisher"))


@pytest.mark.unit
class TestAcatOExclude:
    def test_exclude_acat_v_from_omnibus(self):
        df = _run([_burden(), _skat()], acat_o_exclude="acat_v")
        exp = cauchy({"logistic_burden": 0.08, "skat": 0.2})
        assert df.loc["G1", "acat_o_pvalue"] == pytest.approx(exp, rel=1e-12)
        assert set(df["primary_test"]) == {"acat_o"}

    def test_skat_only_without_acat_v_is_single_test(self):
        # ACAT-O of SKAT alone would be a pass-through, so SKAT is the primary result.
        df = _run([_skat()], acat_o_exclude="acat_v", skat_method="SKATO")
        assert set(df["primary_test"]) == {"skat_o"}
        assert "acat_o_pvalue" not in df.columns
        assert df["primary_pvalue"].to_dict() == SKAT_P

    def test_unknown_exclude_raises(self):
        with pytest.raises(ValueError, match="acat_o_exclude"):
            AssociationEngine([_skat()], AssociationConfig(acat_o_exclude="skat"))


def _parse(argv: list[str]):
    import argparse

    from variantcentrifuge.cli_association_groups import add_association_argument_group

    parser = argparse.ArgumentParser()
    add_association_argument_group(parser)
    return parser.parse_args(argv)


def _validate(argv: list[str], tests: list[str] | str, file_assoc: dict | None = None):
    from unittest.mock import Mock

    from variantcentrifuge.cli_association_config import validate_association_args

    args = _parse(argv)
    args.perform_association = True
    args.association_tests = tests
    fail = Mock(side_effect=SystemExit)
    try:
        cfg = {"association_tests": tests, "association": file_assoc or {}}
        validate_association_args(args, cfg, fail)
    except SystemExit:
        return fail.call_args[0][0]
    return None


@pytest.mark.unit
class TestCliAndConfig:
    def test_cli_defaults_unset_and_mapped(self):
        from variantcentrifuge.cli_association_config import populate_association_config

        args = _parse([])
        assert args.primary_test is None and args.acat_o_exclude is None
        cfg: dict = {"perform_association": True}
        populate_association_config(_parse(["--primary-test", "burden"]), cfg, pytest.fail)
        assert cfg["primary_test"] == "burden" and cfg["acat_o_exclude"] is None

    def test_cli_rejects_primary_not_run(self):
        msg = _validate(["--primary-test", "burden"], ["skat"])
        assert msg and "logistic_burden" in msg
        assert _validate(["--primary-test", "burden"], "logistic_burden,skat") is None
        assert _validate(["--primary-test", "skato"], ["skat"]) is not None
        assert _validate(["--primary-test", "skato", "--skat-method", "SKATO"], ["skat"]) is None
        assert _validate([], ["skat"]) is None

    def test_config_file_primary_fails_before_extraction(self):
        msg = _validate([], ["skat"], {"primary_test": "burden"})
        assert msg and "logistic_burden" in msg
        assert _validate([], ["skat"], {"primary_test": "skato", "skat_method": "SKATO"})
        assert _validate(["--skat-method", "SKATO"], ["skat"], {"primary_test": "skato"}) is None
        # The CLI value wins over the config file.
        assert _validate(["--primary-test", "acat_o"], ["skat"], {"primary_test": "burden"}) is None
        msg = _validate([], ["skat"], {"acat_o_exclude": "skat"})
        assert msg and "acat_o_exclude" in msg

    def test_cli_rejects_unknown_choice(self):
        with pytest.raises(SystemExit):
            _parse(["--acat-o-exclude", "skat"])

    def test_config_file_keys(self):
        from unittest.mock import MagicMock

        from variantcentrifuge.stages.analysis_stage_helpers import (
            _build_assoc_config_from_context,
            _validate_association_config_dict,
        )

        ctx = MagicMock()
        ctx.config = {
            "primary_test": None,
            "acat_o_exclude": None,
            "association": {"primary_test": "skato", "acat_o_exclude": "acat_v"},
        }
        cfg = _build_assoc_config_from_context(ctx)
        assert (cfg.primary_test, cfg.acat_o_exclude) == ("skato", "acat_v")
        ctx.config = {"primary_test": "burden", "association": {"primary_test": "skato"}}
        assert _build_assoc_config_from_context(ctx).primary_test == "burden"
        ctx.config = {}
        cfg = _build_assoc_config_from_context(ctx)
        assert (cfg.primary_test, cfg.acat_o_exclude) == ("acat_o", None)
        with pytest.raises(ValueError, match=r"(?s)primary_test.*acat_o_exclude"):
            _validate_association_config_dict({"primary_test": "x", "acat_o_exclude": "y"})


@pytest.mark.unit
def test_imbalance_warning_suggests_acat_v_exclusion():
    from variantcentrifuge.association.diagnostics import imbalance_message

    msg = imbalance_message(900, 30)
    assert msg is not None and "--acat-o-exclude acat_v" in msg


FISHER_P = {"G1": 1.0, "G2": 1e-6, "G3": 0.3}


def _fisher() -> _FakeTest:
    return _FakeTest("fisher", FISHER_P)


@pytest.mark.unit
class TestFisherNotInAcatO:
    """ACAT-O combines only burden, SKAT/SKAT-O and ACAT-V (Liu 2019; STAAR) (#216, #210)."""

    def test_fisher_with_burden_is_pass_through_of_burden(self):
        df = _run([_fisher(), _burden()])
        assert set(df["primary_test"]) == {"acat_o"}  # schema as any multi-test run
        assert df.loc["G1", "acat_o_pvalue"] == 0.08  # fisher p = 1 no longer vetoes
        assert df.loc["G2", "primary_pvalue"] == 0.01  # fisher p = 1e-6 no longer drives
        assert pd.isna(df.loc["G3", "acat_o_pvalue"])  # burden missing: no fallback to fisher
        assert df.loc["G1", "fisher_pvalue"] == 1.0  # Fisher still reported on its own

    def test_fisher_with_skat_combines_skat_and_acat_v(self):
        df = _run([_fisher(), _skat()])
        for g in GENES:
            exp = cauchy({"skat": SKAT_P[g], "acat_v": ACAT_V_P[g]})
            assert df.loc[g, "acat_o_pvalue"] == pytest.approx(exp, rel=1e-12)

    def test_fisher_burden_skat_equals_run_without_fisher(self):
        with_f = _run([_fisher(), _burden(), _skat()])
        without = _run([_burden(), _skat()])
        cols = ["acat_o_pvalue", "acat_o_qvalue", "primary_pvalue", "primary_qvalue"]
        pd.testing.assert_frame_equal(with_f[cols], without[cols])

    def test_fisher_only_has_no_acat_o(self):
        df = _run([_fisher()])
        assert set(df["primary_test"]) == {"fisher"}
        assert "acat_o_pvalue" not in df.columns

    def test_primary_burden_override_with_fisher(self):
        df = _run([_fisher(), _burden()], primary_test="burden")
        assert set(df["primary_test"]) == {"logistic_burden"}
        assert df.loc["G1", "acat_o_pvalue"] == 0.08
