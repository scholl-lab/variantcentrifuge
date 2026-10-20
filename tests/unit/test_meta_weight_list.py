"""Several weight schemes in one meta run: one suffixed column set per scheme (#202)."""

from __future__ import annotations

import pandas as pd
import pytest

from tests.unit.test_meta_multimask import _two_mask_study, _write
from tests.unit.test_meta_runner_cli import _write_study
from variantcentrifuge.association.meta.cli import main
from variantcentrifuge.association.meta.runner import OUTPUT_COLUMNS, run_meta

SCHEME_FREE = [
    "gene",
    "mask",
    "n_strata",
    "n_cases_total",
    "n_controls_total",
    "n_variants_union",
    "cmac_total",
    "cmac_by_stratum",
    "n_carriers_total",
    "carriers_by_stratum",
    "warnings",
]
WEIGHTED = [c for c in OUTPUT_COLUMNS if c not in SCHEME_FREE]


def test_each_scheme_block_equals_its_single_scheme_run(tmp_path):
    manifest = _write_study(tmp_path)
    multi = run_meta(manifest, str(tmp_path / "m.tsv"), weight_scheme=["beta:1,25", "uniform"])
    assert list(multi.columns) == SCHEME_FREE + [
        f"{c}_{tag}" for tag in ("beta_1_25", "uniform") for c in WEIGHTED
    ]
    for scheme, tag in (("beta:1,25", "beta_1_25"), ("uniform", "uniform")):
        single = run_meta(manifest, str(tmp_path / f"{tag}.tsv"), weight_scheme=scheme)
        got = multi[SCHEME_FREE + [f"{c}_{tag}" for c in WEIGHTED]]
        got.columns = SCHEME_FREE + WEIGHTED
        pd.testing.assert_frame_equal(got[OUTPUT_COLUMNS], single, check_dtype=False)
    # The schemes really differ.
    assert (multi["meta_burden_p_beta_1_25"] != multi["meta_burden_p_uniform"]).any()


def test_single_scheme_list_keeps_todays_names(tmp_path):
    df = run_meta(_write_study(tmp_path), str(tmp_path / "m.tsv"), weight_scheme=["uniform"])
    assert list(df.columns) == OUTPUT_COLUMNS


def test_duplicate_schemes_rejected(tmp_path):
    with pytest.raises(ValueError, match="distinct"):
        run_meta(_write_study(tmp_path), str(tmp_path / "m"), weight_scheme=["uniform"] * 2)


def test_cli_weight_list_with_masks_and_diagnostics(tmp_path):
    manifest = _write(tmp_path, _two_mask_study(tmp_path))
    out, diag = tmp_path / "o.tsv", tmp_path / "diag"
    argv = ["--manifest", manifest, "--output", str(out), "--diagnostics-output", str(diag)]
    assert main([*argv, "--weights", "beta:1,25", "uniform"]) == 0
    df = pd.read_csv(out, sep="\t")
    for tag in ("beta_1_25", "uniform"):
        assert {f"acat_masks_p_{tag}", f"acat_masks_q_{tag}", f"primary_q_{tag}"} <= set(df)
        for mask in ("pLoF", "missense"):
            assert (diag / mask / tag / "lambda_gc.tsv").exists()
    assert "primary_p" not in df.columns


def test_firth_columns_follow_scheme_suffix(tmp_path):
    """#198 x #202: Firth IVW columns are weight-dependent, so they get the scheme suffix."""
    from variantcentrifuge.association.meta.firth_meta import FIRTH_COLUMNS

    df = run_meta(
        _write_study(tmp_path), str(tmp_path / "m.tsv"), weight_scheme=["beta:1,25", "uniform"]
    )
    for col in FIRTH_COLUMNS:
        assert {f"{col}_beta_1_25", f"{col}_uniform"} <= set(df.columns)
        assert col not in df.columns
