"""Meta-analysis gene loop in forked workers == sequential loop, with progress lines."""

from __future__ import annotations

import logging

import pandas as pd
import pytest

from tests.unit.meta_helpers import make_stratum
from tests.unit.test_meta_multimask import _two_mask_study, _write
from variantcentrifuge.association.engine_workers import fork_available
from variantcentrifuge.association.meta.cli import main
from variantcentrifuge.association.meta.runner import meta_analyze_artifacts, run_meta

# Many genes, so fork_map splits them into several chunks per worker.
KEYS = {f"G{i:02d}": [f"{i + 1}:{100 + j}:A:G" for j in range(4 + i % 7)] for i in range(12)}


def _three_strata():
    return [
        make_stratum(s, seed=i + 1, gene_keys=KEYS, effect_gene="G03", effect_sign=1.0)
        for i, s in enumerate(("A", "B", "C"))
    ]


@pytest.mark.skipif(not fork_available(), reason="fork_map needs fork (Linux)")
@pytest.mark.parametrize("workers", [3, -1])
def test_parallel_equals_sequential_with_loo(workers):
    arts = _three_strata()
    seq = meta_analyze_artifacts(arts, workers=1)
    par = meta_analyze_artifacts(arts, workers=workers)
    assert seq["loo_skato_by_stratum"].notna().all()  # LOO ran for every gene
    pd.testing.assert_frame_equal(seq, par, check_exact=True)


@pytest.mark.skipif(not fork_available(), reason="fork_map needs fork (Linux)")
def test_run_meta_parallel_equals_sequential_masks_and_schemes(tmp_path):
    manifest = _write(tmp_path, _two_mask_study(tmp_path))
    schemes = ["beta:1,25", "uniform"]
    seq = run_meta(manifest, str(tmp_path / "seq.tsv"), weight_scheme=schemes)
    par = run_meta(manifest, str(tmp_path / "par.tsv"), weight_scheme=schemes, workers=4)
    pd.testing.assert_frame_equal(seq, par, check_exact=True)
    assert (tmp_path / "seq.tsv").read_bytes() == (tmp_path / "par.tsv").read_bytes()


@pytest.mark.parametrize("workers", [1, 2])
def test_progress_line_logged(workers, caplog):
    caplog.set_level(logging.INFO, logger="variantcentrifuge")
    meta_analyze_artifacts(_three_strata(), workers=workers)
    lines = [r.getMessage() for r in caplog.records if "Meta-analysis genes" in r.getMessage()]
    assert lines and lines[-1].startswith("Meta-analysis genes: 12/12 (100.0%)")


def test_cli_workers_option(tmp_path):
    manifest = _write(tmp_path, _two_mask_study(tmp_path))
    out1, out2 = tmp_path / "w1.tsv", tmp_path / "w2.tsv"
    assert main(["--manifest", manifest, "--output", str(out1)]) == 0
    assert main(["--manifest", manifest, "--output", str(out2), "--workers", "2"]) == 0
    assert out1.read_bytes() == out2.read_bytes()


@pytest.mark.parametrize("bad", ["0", "-2"])
def test_cli_rejects_bad_workers(tmp_path, bad):
    with pytest.raises(SystemExit) as exc:
        main(["--manifest", "m.tsv", "--output", str(tmp_path / "o.tsv"), "--workers", bad])
    assert exc.value.code == 2


@pytest.mark.parametrize("workers", [1, 2])
def test_gene_failure_names_the_gene(workers, monkeypatch):
    import variantcentrifuge.association.meta.runner as runner

    real = runner._gene_row

    def flaky(gene, *a, **k):
        if gene == "G05":
            raise ValueError("boom")
        return real(gene, *a, **k)

    monkeypatch.setattr(runner, "_gene_row", flaky)
    with pytest.raises(ValueError, match="boom") as exc:
        meta_analyze_artifacts(_three_strata(), workers=workers)
    assert "meta-analysis gene G05" in getattr(exc.value, "__notes__", [])


def test_sequential_fallback_without_fork_is_logged(monkeypatch, caplog):
    import variantcentrifuge.association.meta.runner as runner

    monkeypatch.setattr(runner, "fork_available", lambda: False)
    caplog.set_level(logging.INFO, logger="variantcentrifuge")
    meta_analyze_artifacts(_three_strata(), workers=2)
    assert "Parallel meta-analysis needs fork (Linux); running genes sequentially" in caplog.text


def test_cli_error_log_names_the_gene(tmp_path, monkeypatch, caplog):
    import variantcentrifuge.association.meta.runner as runner

    def broken(gene, *a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(runner, "_gene_row", broken)
    manifest = _write(tmp_path, _two_mask_study(tmp_path))
    caplog.set_level(logging.ERROR, logger="variantcentrifuge")
    assert main(["--manifest", manifest, "--output", str(tmp_path / "o.tsv")]) == 1
    assert "Meta-analysis failed: boom; meta-analysis gene NULLG" in caplog.text
