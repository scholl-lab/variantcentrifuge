"""Stage execution order must not depend on PYTHONHASHSEED (#219)."""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

from variantcentrifuge.pipeline_core.runner_planner import create_execution_plan
from variantcentrifuge.pipeline_core.stage import Stage

_PLAN_SCRIPT = """
import argparse, json
from variantcentrifuge.pipeline import build_pipeline_stages
from variantcentrifuge.pipeline_core.runner_planner import create_execution_plan
args = argparse.Namespace(
    config={"ped_file": "x.ped", "inheritance_mode": "simple"},
    perform_association=True, perform_gene_burden=True, ped_file="x.ped",
    inheritance_mode="simple", threads=1,
)
plan = create_execution_plan(build_pipeline_stages(args))
print(json.dumps([[s.name for s in level] for level in plan]))
"""


def _plan_under_seed(seed: str) -> list[list[str]]:
    env = {**os.environ, "PYTHONHASHSEED": seed}
    out = subprocess.run(
        [sys.executable, "-c", _PLAN_SCRIPT], env=env, capture_output=True, text=True, check=True
    ).stdout
    return json.loads(out.strip().splitlines()[-1])


@pytest.mark.unit
def test_full_cli_plan_is_identical_across_hash_seeds():
    plans = {seed: _plan_under_seed(seed) for seed in ("0", "1", "2", "3", "8")}
    flat = [name for level in plans["0"] for name in level]
    assert {"association_analysis", "variant_analysis", "inheritance_analysis"} <= set(flat)
    assert all(plan == plans["0"] for plan in plans.values()), plans


class _S(Stage):
    def __init__(self, name: str, deps: set[str], runtime: float = 1.0):
        self._n, self._d, self._r = name, deps, runtime

    @property
    def name(self) -> str:
        return self._n

    @property
    def dependencies(self) -> set[str]:
        return self._d

    @property
    def estimated_runtime(self) -> float:
        return self._r

    def _process(self, context):  # pragma: no cover - never executed
        return context


@pytest.mark.unit
def test_level_ties_sorted_by_name_after_runtime():
    stages = [
        _S("root", set()),
        _S("zeta", {"root"}),
        _S("alpha", {"root"}),
        _S("slow", {"root"}, runtime=9.0),
        _S("mid", {"root"}),
    ]
    plan = create_execution_plan(stages)
    assert [[s.name for s in level] for level in plan] == [
        ["root"],
        ["slow", "alpha", "mid", "zeta"],
    ]


@pytest.mark.unit
def test_circular_dependency_still_raises():
    with pytest.raises(ValueError, match="Circular"):
        create_execution_plan([_S("a", {"b"}), _S("b", {"a"})])
