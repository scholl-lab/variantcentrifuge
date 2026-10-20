"""Association progress and stage-summary logging (#220)."""

from __future__ import annotations

import logging
import re

import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine
from variantcentrifuge.association.timing import Progress, log_stage_summary


def _genes(n: int) -> list[dict]:
    return [
        {
            "GENE": f"G{i}",
            "proband_count": 100,
            "control_count": 100,
            "proband_carrier_count": 10 + i,
            "control_carrier_count": 5 + i,
            "proband_allele_count": 10 + i,
            "control_allele_count": 5 + i,
            "n_qualifying_variants": 2,
        }
        for i in range(n)
    ]


def test_progress_throttled_by_clock_and_logs_completion(caplog: pytest.LogCaptureFixture) -> None:
    now = [0.0]
    progress = Progress("genes", 10, interval=30.0, clock=lambda: now[0])
    with caplog.at_level(logging.INFO):
        progress.update(1)  # 0 s: throttled
        now[0] = 31.0
        progress.update(4)  # logged
        now[0] = 40.0
        progress.update(5)  # throttled
        progress.update(10)  # completion always logged
    lines = [r.getMessage() for r in caplog.records]
    assert len(lines) == 2
    assert lines[0].startswith("genes: 4/10 (40.0%)") and "ETA 0:00:46" in lines[0]
    assert lines[1].startswith("genes: 10/10 (100.0%)")


def test_stage_summary_lines(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO):
        log_stage_summary({"aggregation": 1.5, "engine": 2.0}, ["a fact"])
    text = "\n".join(r.getMessage() for r in caplog.records)
    assert "sub-step aggregation: 1.5 s" in text
    assert "total: 3.5 s, peak RSS" in text
    assert "Association: a fact" in text


@pytest.mark.parametrize("workers", [1, 2])
def test_engine_progress_and_decisions(
    workers: int, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("variantcentrifuge.association.engine.MIN_GENES_FOR_PARALLEL", 2)
    engine = AssociationEngine.from_names(
        ["fisher"], AssociationConfig(association_workers=workers)
    )
    with caplog.at_level(logging.INFO):
        engine.run_all(_genes(8))
    progress = [r for r in caplog.records if "Association genes:" in r.getMessage()]
    assert progress and "(100.0%)" in progress[-1].getMessage()  # completion line
    assert "BLAS pinned" in engine.run_info[0]
    assert ("parallel" in engine.run_info[0]) == (workers == 2)


def test_rss_units_by_platform(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    import resource
    from types import SimpleNamespace

    monkeypatch.setattr(
        resource, "getrusage", lambda _who: SimpleNamespace(ru_maxrss=2 * 1024 * 1024)
    )
    for platform, mib in (("linux", 2048), ("darwin", 2)):
        monkeypatch.setattr("sys.platform", platform)
        caplog.clear()
        with caplog.at_level(logging.INFO):
            log_stage_summary({}, [])
        assert f"peak RSS {mib} MiB" in caplog.text


def test_subset_summary_variants() -> None:
    from types import SimpleNamespace

    from variantcentrifuge.sample_subset import subset_summary

    def ctx(**cfg):
        return SimpleNamespace(config=cfg, vcf_samples=["a", "b", "c"])

    assert subset_summary(ctx(sample_subset_sha256="0123456789abcdef")) == (
        "sample subset: 3 VCF samples extracted (sha 01234567)"
    )
    assert "--no-sample-subset" in subset_summary(ctx(no_sample_subset=True))
    assert "--igv" in subset_summary(ctx(igv_enabled=True, no_stats=True))
    assert "no sample file" in subset_summary(ctx(no_stats=True))


@pytest.mark.parametrize("workers", [1, 4])
def test_progress_total_is_all_genes(
    workers: int, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("variantcentrifuge.association.engine.MIN_GENES_FOR_PARALLEL", 2)
    engine = AssociationEngine.from_names(
        ["fisher"], AssociationConfig(association_workers=workers)
    )
    with caplog.at_level(logging.INFO):
        engine.run_all(_genes(100))
    last = [r.getMessage() for r in caplog.records if "Association genes:" in r.getMessage()][-1]
    assert "100/100 (100.0%)" in last
    given, with_p, without = map(int, re.findall(r"\d+", engine.run_info[1]))
    assert (given, with_p + without) == (100, 100) and without > 0  # synthetic extremes untestable


def test_gene_sets_get_own_label_and_run_info(
    caplog: pytest.LogCaptureFixture,
) -> None:
    from variantcentrifuge.association.gene_sets import run_gene_and_set_units

    cfg = AssociationConfig()
    engine = AssociationEngine.from_names(["fisher"], cfg)
    units = _genes(3) + [dict(g, GENE=f"SET{i}") for i, g in enumerate(_genes(2))]
    with caplog.at_level(logging.INFO):
        run_gene_and_set_units(engine, ["fisher"], cfg, units, ["SET0", "SET1"], None)
    text = caplog.text
    assert "Association genes: 3/3" in text and "Association gene sets: 2/2" in text
    assert any(i.startswith("gene sets: 2 given") for i in engine.run_info)
    assert any(i.startswith("genes: 3 given") for i in engine.run_info)


def test_summary_omits_zero_steps_and_child_rss(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO):
        log_stage_summary({"sample_setup": 0.0, "engine": 1.0}, [])
    assert "sample_setup" not in caplog.text and "worker" not in caplog.text
