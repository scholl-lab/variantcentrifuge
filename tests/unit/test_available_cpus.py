"""-1 / "all cores" worker counts mean the CPUs this process may use (SLURM/cgroup affinity)."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import variantcentrifuge
from tests.unit.test_meta_parallel import _three_strata
from tests.unit.test_parallel_association import _make_gene_data
from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine, dispatch_parallel_genes
from variantcentrifuge.association.meta.runner import meta_analyze_artifacts
from variantcentrifuge.stages.analysis_stage_helpers import _build_assoc_config_from_context


@pytest.fixture
def four_of_64(monkeypatch):
    """A 4-CPU job on a 64-core node."""
    monkeypatch.setattr(os, "cpu_count", lambda: 64)
    monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: set(range(4)), raising=False)


def test_available_cpus_uses_affinity(four_of_64):
    from variantcentrifuge.association.engine_workers import available_cpus

    assert available_cpus() == 4


@pytest.mark.parametrize(("count", "expected"), [(8, 8), (None, 1)])
def test_available_cpus_without_affinity(monkeypatch, count, expected):
    from variantcentrifuge.association.engine_workers import available_cpus

    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: count)
    assert available_cpus() == expected


def test_association_workers_minus_one(four_of_64):
    context = MagicMock()
    context.config = {"association_workers": -1}
    assert _build_assoc_config_from_context(context).association_workers == 4


def test_engine_minus_one(four_of_64, monkeypatch):
    monkeypatch.setattr("variantcentrifuge.association.engine.MIN_GENES_FOR_PARALLEL", 2)
    engine = AssociationEngine.from_names(["fisher"], AssociationConfig(association_workers=-1))
    target = "variantcentrifuge.association.engine.dispatch_parallel_genes"
    with patch(target, wraps=dispatch_parallel_genes) as dispatch:
        engine.run_all(_make_gene_data(20))
    assert dispatch.call_args.kwargs["actual_workers"] == 4


def test_meta_cap(four_of_64, caplog):
    caplog.set_level(logging.INFO, logger="variantcentrifuge")
    with patch("variantcentrifuge.association.meta.runner.fork_available", lambda: False):
        meta_analyze_artifacts(_three_strata(), workers=64)
    assert "--workers 64 capped at 4 CPU cores" in caplog.text


def test_no_raw_cpu_count_in_worker_resolution():
    """Every association/stage worker count resolves through available_cpus()."""
    root = Path(variantcentrifuge.__file__).parent
    hits = [
        f"{p.relative_to(root)}:{i}"
        for sub in ("association", "stages")
        for p in sorted((root / sub).rglob("*.py"))
        if p.name != "engine_workers.py"
        for i, line in enumerate(p.read_text().splitlines(), 1)
        if re.search(r"os\.cpu_count\(\)", line)
    ]
    assert hits == []
