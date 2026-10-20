"""End-to-end runner + CLI: manifest → meta TSV, and refusal paths."""

from __future__ import annotations

import pandas as pd

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.meta.artifact import save_artifact
from variantcentrifuge.association.meta.cli import main
from variantcentrifuge.association.meta.runner import OUTPUT_COLUMNS, run_meta

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
