"""Wall-clock regression guard for the association stage (size S, synthetic data).

Refresh the baseline after an intended speed change, on the machine that runs this test:
    python scripts/bench_association.py --size S --repeat 3 \
        --write-baseline tests/performance/association_bench_baseline.json
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
BASELINE = Path(__file__).with_name("association_bench_baseline.json")


@pytest.mark.slow
@pytest.mark.performance
def test_association_stage_size_s_not_slower_than_baseline():
    base = json.loads(BASELINE.read_text())
    # Subprocess so BLAS thread variables are set before numpy is first imported.
    cmd = [sys.executable, str(ROOT / "scripts" / "bench_association.py"), "--size", "S"]
    cmd += ["--repeat", "3", "--workers", "1", "--blas-threads", "1"]
    out = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True).stdout
    line = next(ln for ln in out.splitlines() if ln.startswith("RESULT "))
    result = dict(kv.split("=", 1) for kv in line.split()[1:])
    assert int(result["rows"]) == base["rows"]
    assert float(result["stage_s"]) <= 1.2 * base["stage_s"], (result, base)


def _bench(*args: str) -> dict[str, str]:
    cmd = [sys.executable, str(ROOT / "scripts" / "bench_association.py"), *args]
    out = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True).stdout
    line = next(ln for ln in out.splitlines() if ln.startswith("RESULT "))
    return dict(kv.split("=", 1) for kv in line.split()[1:])


@pytest.mark.slow
@pytest.mark.performance
def test_meta_mode_runs_on_two_disjoint_artifacts(tmp_path):
    common = ["--size", "S", "--genes", "8", "--samples", "300", "--artifact"]
    _bench(*common, "--out", str(tmp_path / "a.tsv"))
    _bench(*common, "--seed", "43", "--sample-prefix", "T", "--out", str(tmp_path / "b.tsv"))
    npz = [str(tmp_path / f"{x}.tsv.scores.npz") for x in "ab"]
    result = _bench("--meta", *npz, "--out", str(tmp_path / "meta.tsv"))
    assert result["rc"] == "0" and int(result["rows"]) == 8, result
    assert (tmp_path / "meta.tsv.pkl").exists()


def test_pss_peak_samples_short_runs():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "bench", ROOT / "scripts" / "bench_association.py"
    )
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    with bench.PssPeak() as pss:
        pass
    assert pss.peak_mb > 0
