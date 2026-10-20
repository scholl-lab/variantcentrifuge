"""End-to-end runner + CLI: manifest → meta TSV, and refusal paths."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.meta.artifact import save_artifact
from variantcentrifuge.association.meta.cli import main
from variantcentrifuge.association.meta.runner import (
    OUTPUT_COLUMNS,
    meta_analyze_artifacts,
    run_meta,
)

GENE_KEYS = {
    "SIG": [f"1:{100 + j}:A:G" for j in range(14)],
    "NULLG": [f"2:{100 + j}:A:G" for j in range(10)],
}


def _write_study(tmp_path, effect_signs=(+1.0, +1.0), **compat):
    a = make_stratum(
        "A", seed=1, gene_keys=GENE_KEYS, effect_gene="SIG", effect_sign=effect_signs[0], **compat
    )
    b = make_stratum(
        "B", seed=2, gene_keys=GENE_KEYS, effect_gene="SIG", effect_sign=effect_signs[1]
    )
    pa, pb = str(tmp_path / "A"), str(tmp_path / "B")
    save_artifact(a, pa)
    save_artifact(b, pb)
    manifest = tmp_path / "study.tsv"
    manifest.write_text(f"stratum_id\tartifact_path\nA\t{pa}\nB\t{pb}\n")
    return str(manifest)


def test_run_meta_produces_expected_schema_and_correction(tmp_path):
    manifest = _write_study(tmp_path)
    out = str(tmp_path / "meta.tsv")
    df = run_meta(manifest, out, weight_scheme="beta:1,25", correction="fdr")

    assert list(df.columns) == OUTPUT_COLUMNS
    assert set(df["gene"]) == {"SIG", "NULLG"}
    on_disk = pd.read_csv(out, sep="\t")
    assert len(on_disk) == 2
    # Every tested gene carries n_strata=2 and an "ok" status.
    assert (df.loc[df["primary_p"].notna(), "n_strata"] == 2).all()
    assert set(df.loc[df["primary_p"].notna(), "status"]) <= {"ok"}
    # q-values are valid probabilities.
    q = df["primary_q"].dropna()
    assert ((q >= 0) & (q <= 1)).all()


def test_run_meta_detects_concordant_signal(tmp_path):
    manifest = _write_study(tmp_path, effect_signs=(+1.0, +1.0))
    out = str(tmp_path / "meta.tsv")
    df = run_meta(manifest, out).set_index("gene")
    # The signal gene is more significant than the null gene.
    assert df.loc["SIG", "primary_p"] < df.loc["NULLG", "primary_p"]


def test_meta_subcommand_success(tmp_path):
    from variantcentrifuge.cli import main as vc_main

    manifest = _write_study(tmp_path)
    out = str(tmp_path / "meta.tsv")
    rc = vc_main(["meta", "--manifest", manifest, "--output", out, "--weights", "beta:1,25"])
    assert rc == 0
    assert (tmp_path / "meta.tsv").exists()


def test_cli_refuses_incomparable_strata(tmp_path):
    # Stratum B built on a different genome build -> comparability refusal (exit 2).
    a = make_stratum("A", seed=1, gene_keys=GENE_KEYS)
    b = make_stratum("B", seed=2, gene_keys=GENE_KEYS, genome_build="GRCh37")
    pa, pb = str(tmp_path / "A"), str(tmp_path / "B")
    save_artifact(a, pa)
    save_artifact(b, pb)
    manifest = tmp_path / "study.tsv"
    manifest.write_text(f"stratum_id\tartifact_path\nA\t{pa}\nB\t{pb}\n")
    rc = main(["--manifest", str(manifest), "--output", str(tmp_path / "out.tsv")])
    assert rc == 2


def test_single_variant_gene_falls_back_to_burden(tmp_path):
    """Rank-1 genes crashed SKAT-O (ZeroDivisionError) on real 1000 Genomes data."""
    keys = {"ONE": ["3:500:C:T"], "SIG": GENE_KEYS["SIG"]}
    strata = [make_stratum(s, seed=i, gene_keys=keys) for i, s in enumerate(("A", "B"))]
    for art in strata:
        save_artifact(art, str(tmp_path / art.stratum_id))
    manifest = tmp_path / "study.tsv"
    manifest.write_text("stratum_id\tartifact_path\nA\tA\nB\tB\n")
    row = run_meta(str(manifest), str(tmp_path / "out.tsv")).set_index("gene").loc["ONE"]
    assert pd.isna(row["meta_skato_p"])
    assert row["primary_p"] == row["meta_burden_p"]
    assert row["status"] == "ok"


# --- #192 / #193 / #194: diagnostics, named LOO, cumulative MAC -------------------


def _three_strata():
    return [
        make_stratum(s, seed=i + 1, gene_keys=GENE_KEYS, effect_gene="SIG", effect_sign=1.0)
        for i, s in enumerate(("A", "B", "C"))
    ]


def _parse(field):
    """'A:0.1:+1;B:NA' -> {'A': ['0.1', '+1'], 'B': ['NA']}"""
    return {p.split(":")[0]: p.split(":")[1:] for p in field.split(";")}


def test_loo_equals_direct_meta_on_remaining_strata():
    arts = _three_strata()
    full = meta_analyze_artifacts(arts).set_index("gene")
    for gene in ("SIG", "NULLG"):
        row = full.loc[gene]
        burden = _parse(row["loo_burden_by_stratum"])
        skato = _parse(row["loo_skato_by_stratum"])
        assert set(burden) == set(skato) == {"A", "B", "C"}
        for drop in ("A", "B", "C"):
            sub = meta_analyze_artifacts([a for a in arts if a.stratum_id != drop])
            ref = sub.set_index("gene").loc[gene]
            assert float(burden[drop][0]) == pytest.approx(ref["meta_burden_p"], rel=1e-9)
            assert int(burden[drop][1]) == ref["meta_burden_direction"]
            assert float(skato[drop][0]) == pytest.approx(ref["meta_skato_p"], rel=1e-9)
        bmax = max(burden, key=lambda s: float(burden[s][0]))
        smax = max(skato, key=lambda s: float(skato[s][0]))
        assert row["loo_max_burden_stratum"] == bmax
        assert row["loo_max_burden_p"] == pytest.approx(float(burden[bmax][0]), rel=1e-9)
        assert row["loo_max_skato_stratum"] == smax
        assert row["loo_max_skato_p"] == pytest.approx(float(skato[smax][0]), rel=1e-9)
        assert row["loo_min_p"] == pytest.approx(min(float(b[0]) for b in burden.values()))


def test_single_stratum_has_no_loo():
    row = meta_analyze_artifacts(_three_strata()[:1]).iloc[0]
    assert pd.isna(row["loo_burden_by_stratum"]) and pd.isna(row["loo_max_skato_p"])


def test_cumulative_mac_from_artifact_ac_an():
    arts = _three_strata()
    df = meta_analyze_artifacts(arts).set_index("gene")
    for gene in ("SIG", "NULLG"):
        assert all((np.diag(a.genes[gene].v) > 0).all() for a in arts)  # all informative
        per = {
            a.stratum_id: np.minimum(a.genes[gene].ac, a.genes[gene].an - a.genes[gene].ac).sum()
            for a in arts
        }
        assert df.loc[gene, "cmac_total"] == pytest.approx(sum(per.values()))
        assert {k: float(v[0]) for k, v in _parse(df.loc[gene, "cmac_by_stratum"]).items()} == (
            pytest.approx(per)
        )


def test_min_mac_flags_low_mac_and_excludes_from_correction():
    arts = _three_strata()
    base = meta_analyze_artifacts(arts, correction="bonferroni").set_index("gene")
    cut = (base["cmac_total"].min() + base["cmac_total"].max()) / 2
    low = base["cmac_total"].idxmin()
    high = base["cmac_total"].idxmax()
    df = meta_analyze_artifacts(arts, correction="bonferroni", min_mac=cut).set_index("gene")
    assert df.loc[low, "status"] == "low_mac"
    assert pd.notna(df.loc[low, "primary_p"])  # still reported, not dropped
    assert pd.isna(df.loc[low, "primary_q"])  # but outside the correction family
    assert df.loc[high, "status"] == "ok"
    assert df.loc[high, "primary_q"] == pytest.approx(min(1.0, df.loc[high, "primary_p"]))
    # Without the filter both genes share the Bonferroni family of two.
    assert base.loc[high, "primary_q"] == pytest.approx(min(1.0, 2 * base.loc[high, "primary_p"]))


def test_cli_min_mac_and_diagnostics_output(tmp_path):
    manifest = _write_study(tmp_path)
    out, diag = tmp_path / "meta.tsv", tmp_path / "diag"
    rc = main(
        [
            "--manifest",
            manifest,
            "--output",
            str(out),
            "--diagnostics-output",
            str(diag),
            "--min-mac",
            "1e9",
        ]
    )
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert set(df["status"]) == {"low_mac"}
    # Every gene is low-MAC, so no p-value enters lambda/QQ.
    assert pd.read_csv(diag / "lambda_gc.tsv", sep="\t").empty


def test_diagnostics_output_lambda_per_meta_test(tmp_path):
    manifest = _write_study(tmp_path)
    diag = tmp_path / "diag"
    df = run_meta(manifest, str(tmp_path / "meta.tsv"), diagnostics_dir=str(diag))
    lam = pd.read_csv(diag / "lambda_gc.tsv", sep="\t").set_index("test_name")
    assert set(lam.index) == {"meta_burden", "meta_skat", "meta_skato", "primary"}
    assert (lam["n_tests"] == 2).all()
    qq = pd.read_csv(diag / "qq_data.tsv", sep="\t")
    assert set(qq["test"]) == set(lam.index)
    summary = (diag / "summary.txt").read_text()
    assert f"n_cases    = {int(df['n_cases_total'].max())}" in summary
    assert "Genes with per-gene warnings: 0" in summary


def _save_study(tmp_path, arts):
    lines = ["stratum_id\tartifact_path"]
    for art in arts:
        save_artifact(art, str(tmp_path / art.stratum_id))
        lines.append(f"{art.stratum_id}\t{art.stratum_id}")
    manifest = tmp_path / "study.tsv"
    manifest.write_text("\n".join(lines) + "\n")
    return str(manifest)


def test_diagnostics_family_equals_correction_family(tmp_path):
    arts = _three_strata()[:2]
    for art in arts:  # NULLG carries no information -> degenerate, but meta_skat_p == 1
        gs = art.genes["NULLG"]
        gs.u, gs.v = np.zeros_like(gs.u), np.zeros_like(gs.v)
    arts[0].genes["SIG"].null_converged = False  # SIG keeps a per-gene warning
    manifest = _save_study(tmp_path, arts)
    diag = tmp_path / "diag"
    df = run_meta(manifest, str(tmp_path / "m.tsv"), diagnostics_dir=str(diag)).set_index("gene")
    assert df.loc["NULLG", "status"] == "degenerate" and df.loc["NULLG", "meta_skat_p"] == 1.0
    lam = pd.read_csv(diag / "lambda_gc.tsv", sep="\t").set_index("test_name")
    assert (lam["n_tests"] == 1).all()  # only SIG is in the correction family
    # Low-MAC genes leave lambda/QQ but their warnings are still counted.
    run_meta(manifest, str(tmp_path / "m2.tsv"), min_mac=1e9, diagnostics_dir=str(diag))
    assert pd.read_csv(diag / "lambda_gc.tsv", sep="\t").empty
    assert "Genes with per-gene warnings: 1" in (diag / "summary.txt").read_text()


def test_mac_ignores_zero_information_variants():
    """A locus fixed REF in one stratum and fixed ALT in the other adds no score
    information, so it must not count toward cumulative MAC / --min-mac."""
    arts = _three_strata()[:2]
    before = meta_analyze_artifacts(arts).set_index("gene").loc["NULLG"]
    for art, fixed_alt in zip(arts, (False, True), strict=True):
        gs = art.genes["NULLG"]
        gs.u[0], gs.v[0, :], gs.v[:, 0] = 0.0, 0.0, 0.0
        gs.ac[0] = gs.an[0] if fixed_alt else 0.0
    row = meta_analyze_artifacts(arts, min_mac=1.0).set_index("gene").loc["NULLG"]
    per = {
        a.stratum_id: int(np.minimum(g.ac, g.an - g.ac)[1:].sum())
        for a in arts
        for g in [a.genes["NULLG"]]
    }
    expected = sum(per.values())
    assert row["cmac_total"] == pytest.approx(expected)
    assert row["cmac_total"] < before["cmac_total"]
    assert row["cmac_by_stratum"] == ";".join(f"{s}:{m}" for s, m in per.items())


@pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf")])
def test_min_mac_must_be_finite_and_nonnegative(tmp_path, bad):
    with pytest.raises(ValueError, match="min_mac"):
        meta_analyze_artifacts(_three_strata()[:2], min_mac=bad)
    manifest = _write_study(tmp_path)
    assert (
        main(["--manifest", manifest, "--output", str(tmp_path / "o.tsv"), "--min-mac", str(bad)])
        == 1
    )


@pytest.mark.parametrize("name", ["A;B", "A:B", "A,B"])
def test_stratum_ids_with_output_separators_are_rejected(name):
    arts = _three_strata()[:2]
    arts[0].stratum_id = name
    with pytest.raises(ValueError, match="separator"):
        meta_analyze_artifacts(arts)


def test_require_variants_in_all_strata(tmp_path):
    keys = {**GENE_KEYS, "ONLYA": [f"4:{100 + j}:A:G" for j in range(5)]}
    a = make_stratum("A", seed=1, gene_keys=keys, effect_gene="SIG", effect_sign=1.0)
    b = make_stratum("B", seed=2, gene_keys=GENE_KEYS, effect_gene="SIG", effect_sign=1.0)
    gs = b.genes["NULLG"]  # present in B but without any informative variant there
    gs.u, gs.v = np.zeros_like(gs.u), np.zeros_like(gs.v)
    loose = meta_analyze_artifacts([a, b], correction="bonferroni").set_index("gene")
    assert loose.loc["ONLYA", "status"] == "single_stratum"
    assert loose["primary_q"].notna().all()

    df = meta_analyze_artifacts(
        [a, b], correction="bonferroni", require_variants_in_all_strata=True
    ).set_index("gene")
    assert list(df["status"]) == ["not_in_all_strata", "not_in_all_strata", "ok"]
    assert df.loc[["NULLG", "ONLYA"], "primary_p"].notna().all()  # reported, not dropped
    assert df.loc[["NULLG", "ONLYA"], "primary_q"].isna().all()  # outside the family
    assert df.loc["SIG", "primary_q"] == pytest.approx(min(1.0, df.loc["SIG", "primary_p"]))

    manifest = _save_study(tmp_path, [a, b])
    out = tmp_path / "o.tsv"
    assert (
        main(["--manifest", manifest, "--output", str(out), "--require-variants-in-all-strata"])
        == 0
    )
    assert (pd.read_csv(out, sep="\t")["status"] == "not_in_all_strata").sum() == 2


def test_cmac_total_is_sum_of_stratum_macs_with_fixed_alt_site():
    arts = _three_strata()[:2]
    gs = arts[0].genes["NULLG"]  # site 0 fixed ALT in A (no information there)
    gs.u[0], gs.v[0, :], gs.v[:, 0], gs.ac[0] = 0.0, 0.0, 0.0, gs.an[0]
    row = meta_analyze_artifacts(arts).set_index("gene").loc["NULLG"]
    per = {k: float(v[0]) for k, v in _parse(row["cmac_by_stratum"]).items()}
    assert row["cmac_total"] == sum(per.values())
