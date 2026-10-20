"""Association progress and stage-summary logging (#220)."""

from __future__ import annotations

import logging

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
