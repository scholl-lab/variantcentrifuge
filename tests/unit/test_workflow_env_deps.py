"""The Snakemake workflow conda env must provide every runtime dependency in pyproject.toml."""

import re
import tomllib
from pathlib import Path

from packaging import requirements

ROOT = Path(__file__).resolve().parents[2]


def _norm(name: str) -> str:
    return name.lower().replace("_", "-")


def test_workflow_env_has_all_pyproject_dependencies():
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["dependencies"]
    env = (ROOT / "workflow" / "envs" / "variantcentrifuge.yml").read_text()
    entries = [ln.strip()[2:].strip() for ln in env.splitlines() if ln.strip().startswith("- ")]
    listed = {_norm(re.split(r"[=<>\s#]", e)[0]) for e in entries}
    missing = [d for d in deps if _norm(requirements.Requirement(d).name) not in listed]
    assert not missing, f"workflow env lacks pyproject dependencies: {missing}"
