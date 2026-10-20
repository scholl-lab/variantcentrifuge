"""--gene-set-file: gene sets as extra association testing units (#201)."""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.gene_sets import add_gene_set_units, load_gene_sets
from variantcentrifuge.association.meta.artifact import load_artifact
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage

N = 80
CASES = [f"CASE{i}" for i in range(N // 2)]
CONTROLS = [f"CTRL{i}" for i in range(N // 2)]
GT = {0: "0/0", 1: "0/1", 2: "1/1"}


def _variants() -> pd.DataFrame:
    """Genes A (2 variants), B (1), C (1); variant at POS 40 is annotated to both A and B."""
    rng = np.random.default_rng(201)
    rows = [("A", 10), ("A", 20), ("B", 30), ("C", 50), ("A", 40), ("B", 40)]
    dosage = {pos: rng.binomial(2, 0.08, N) for pos in (10, 20, 30, 40, 50)}
    dosage[10][: N // 4] = 1  # case-enriched signal so tests are not degenerate
    return pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": pos,
                "REF": "A",
                "ALT": "G",
                "GENE": gene,
                **{f"GEN_{i}__GT": GT[int(dosage[pos][i])] for i in range(N)},
            }
            for gene, pos in rows
        ]
    )


def _run(tmp_path, df, name, **extra):
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.get_output_path = lambda n, ext=".tsv": tmp_path / f"{name}{ext}"
    workspace.base_name = name
    config = {
        "perform_association": True,
        "case_samples": CASES,
        "control_samples": CONTROLS,
        "association_tests": ["fisher", "logistic_burden"],
        "association_min_cases": 1,
        "gzip_intermediates": False,
        **extra,
    }
    ctx = PipelineContext(args=Namespace(), config=config, workspace=workspace)
    ctx.current_dataframe = df
    ctx.vcf_samples = CASES + CONTROLS
    AssociationAnalysisStage()._process(ctx)
    return ctx.association_results.set_index("gene")


def _set_file(tmp_path, text):
    path = tmp_path / "sets.tsv"
    path.write_text(text)
    return str(path)


def _stats(row: pd.Series) -> pd.Series:
    return row.drop([c for c in row.index if "qvalue" in c or c == "unit_type"])


def test_set_units_match_gene_and_merged_gene_runs(tmp_path):
    df = _variants()
    sets = _set_file(tmp_path, "set_id\tgene\nONE\tA\nTWO\tA\nTWO\tB\n")
    res = _run(tmp_path, df, "with_sets", gene_set_file=sets)
    base = _run(tmp_path, df, "genes_only")

    assert res.loc[["A", "B", "C"], "unit_type"].eq("gene").all()
    assert res.loc[["ONE", "TWO"], "unit_type"].eq("set").all()
    # Genes keep their own correction family: identical to a run without sets.
    pd.testing.assert_frame_equal(res.loc[base.index, base.columns], base, check_dtype=False)
    # A one-gene set is exactly the gene.
    pd.testing.assert_series_equal(_stats(res.loc["ONE"]), _stats(res.loc["A"]), check_names=False)
    # A two-gene set equals one merged gene, with the shared variant counted once.
    merged = df.assign(GENE=df["GENE"].replace({"A": "AB", "B": "AB"}))
    ref = _run(tmp_path, merged, "merged").loc["AB"]
    assert res.loc["TWO", "n_variants"] == 4
    pd.testing.assert_series_equal(_stats(res.loc["TWO"]), _stats(ref), check_names=False)


def test_sets_are_ordinary_units_in_score_artifact(tmp_path):
    sets = _set_file(tmp_path, "set_id\tgene\nTWO\tA\nTWO\tB\n")
    artifact = str(tmp_path / "cohort.npz")
    _run(
        tmp_path,
        _variants(),
        "art",
        gene_set_file=sets,
        score_artifact=artifact,
        reference="GRCh38.99",
        diagnostics_output=str(tmp_path / "diag"),
    )
    assert (tmp_path / "diag" / "lambda_gc.tsv").exists()
    art = load_artifact(artifact)
    assert art.genes["TWO"].variant_keys == ["1:10:A:G", "1:20:A:G", "1:30:A:G", "1:40:A:G"]
    np.testing.assert_allclose(art.genes["TWO"].u[:2], art.genes["A"].u[:2])


@pytest.mark.parametrize(
    ("text", "match"),
    [
        ("set\tgene\nS\tA\n", "set_id"),
        ("set_id\tgene\nS\tA\nS\tA\n", "duplicate"),
        ("set_id\tgene\n\tA\n", "empty set_id"),
        ("set_id\tgene\nS::x\tA\n", "::"),
    ],
)
def test_load_gene_sets_rejects_bad_files(tmp_path, text, match):
    with pytest.raises(ValueError, match=match):
        load_gene_sets(_set_file(tmp_path, text))


def test_set_id_equal_to_gene_symbol_is_an_error():
    with pytest.raises(ValueError, match="gene symbols"):
        add_gene_set_units(_variants(), {"A": ["B"]})


def test_unknown_genes_warn_with_count_and_empty_sets_are_skipped(tmp_path, caplog):
    sets = load_gene_sets(_set_file(tmp_path, "set_id\tgene\nS\tA\nS\tX1\nS\tX2\nE\tX3\nN\t\n"))
    assert sets == {"S": ["A", "X1", "X2"], "E": ["X3"], "N": []}
    out, added = add_gene_set_units(_variants(), sets)
    assert added == ["S"]
    assert (out["GENE"] == "S").sum() == 3
    assert "3 gene(s)" in caplog.text
    assert "'E'" in caplog.text and "'N'" in caplog.text
    df = _variants()
    out, added = add_gene_set_units(df, {"E": ["X3"]})
    assert out is df and added == []


def _parse(argv: list[str]) -> Namespace:
    import argparse

    from variantcentrifuge.cli_association_groups import add_association_argument_group

    parser = argparse.ArgumentParser()
    add_association_argument_group(parser)
    return parser.parse_args(argv)


def test_cli_flag_maps_to_config_and_checks_file(tmp_path):
    from variantcentrifuge.cli_association_config import populate_association_config

    path = _set_file(tmp_path, "set_id\tgene\nS\tA\n")
    args = _parse(["--gene-set-file", path])
    cfg: dict = {}
    populate_association_config(args, cfg, pytest.fail)
    assert cfg["gene_set_file"] == path

    fail = Mock(side_effect=SystemExit)
    with pytest.raises(SystemExit):
        populate_association_config(_parse(["--gene-set-file", "nope.tsv"]), {}, fail)
    assert "--gene-set-file file not found" in fail.call_args[0][0]


def test_cli_flag_requires_perform_association(tmp_path):
    from variantcentrifuge.cli_association_config import validate_association_args

    args = _parse(["--gene-set-file", "x.tsv"])
    args.perform_association = False
    args.association_tests = None
    fail = Mock(side_effect=SystemExit)
    with pytest.raises(SystemExit):
        validate_association_args(args, {}, fail)
    assert "--gene-set-file requires --perform-association" in fail.call_args[0][0]
