#!/usr/bin/env python
"""Seeded synthetic benchmark for ``AssociationAnalysisStage._process``.

Builds a synthetic variant table (no real data), runs the association stage in-process
with a mocked Workspace and prints one machine-readable ``RESULT`` line plus a table of
the stage's sub-step timings (``context.stage_results["association_analysis"]``).

Examples
--------
    python scripts/bench_association.py --size S
    python scripts/bench_association.py --size L --layout columns --tests fisher --out a.tsv
    python scripts/bench_association.py --compare before.tsv after.tsv
    python scripts/bench_association.py --size S --repeat 3 \
        --write-baseline tests/performance/association_bench_baseline.json
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import sys
import tempfile
import time
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

PRESETS = {"S": (50, 1000), "M": (500, 3000), "L": (2000, 6000)}  # (genes, samples)
BLAS_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
GT_CODES = ["0/0", "0/1", "1/1", "./."]
CHUNK_CELLS = 4_000_000  # bounds the transient float64 draw to ~32 MB per chunk


def make_df(n_genes, n_samples, layout, seed=42, vmin=5, vmax=45, miss=0.005):
    """Return (variants DataFrame, sample IDs); same genotypes for both layouts."""
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    n_var = rng.integers(vmin, vmax + 1, n_genes)
    genes = np.repeat([f"G{i:05d}" for i in range(n_genes)], n_var)
    n = len(genes)
    af = np.clip(rng.beta(0.5, 150, n), 1e-4, 0.05)
    effect = rng.choice(["missense_variant", "stop_gained", "synonymous_variant"], n)
    samples = [f"S{i:05d}" for i in range(n_samples)]
    sids = np.array(samples, dtype=object)
    codes = np.empty((n, n_samples), dtype=np.int8) if layout == "columns" else None
    packed: list[str] = []
    step = max(1, CHUNK_CELLS // n_samples)
    for lo in range(0, n, step):
        hi = min(n, lo + step)
        u = rng.random((hi - lo, n_samples))
        a = af[lo:hi, None]
        c = (u < a**2).astype(np.int8) + (u < 2 * a).astype(np.int8)  # binomial(2, af)
        c[rng.random((hi - lo, n_samples)) < miss] = 3
        for row in c:
            nz = np.flatnonzero(row)
            packed.append(";".join(f"{sids[j]}({GT_CODES[row[j]]})" for j in nz))
        if codes is not None:
            codes[lo:hi] = c
    df = pd.DataFrame(
        {
            "CHROM": "1",
            "POS": np.arange(n) * 10 + 1000,
            "REF": "A",
            "ALT": "G",
            "GENE": genes,
            "EFFECT": effect,
            "IMPACT": "MODERATE",
            "GT": packed,
        }
    )
    if codes is not None:  # per-sample categorical GT columns, as after dataframe_optimizer
        cols = {
            f"GEN_{j}__GT": pd.Categorical.from_codes(codes[:, j], categories=GT_CODES)
            for j in range(n_samples)
        }
        del codes
        df = pd.concat([df, pd.DataFrame(cols)], axis=1)
    return df, samples


def run_stage(df, samples, args, tmp):
    """Run the association stage once; return (results, wall seconds, sub-step seconds)."""
    import numpy as np

    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace
    from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage

    rng = np.random.default_rng(args.seed + 1)
    perm = list(rng.permutation(samples))
    n_cases = int(len(samples) * args.case_frac)
    n_assigned = int(len(samples) * (1 - args.unassigned))
    ws = Mock(spec=Workspace)
    ws.output_dir = ws.intermediate_dir = tmp
    cfg = {
        "perform_association": True,
        "case_samples": perm[:n_cases],
        "control_samples": perm[n_cases:n_assigned],
        "association_tests": args.tests.split(","),
        "output_dir": str(tmp),
        "output_file_base": "bench",
        "gzip_intermediates": False,
        "association_workers": args.workers,
        "threads": args.workers,
        "skat_method": args.skat_method,
    }
    if args.artifact:
        cfg["score_artifact"] = str(tmp / "bench.scores.npz")
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=ws)
    ctx.current_dataframe = df
    ctx.vcf_samples = samples
    stage = AssociationAnalysisStage()
    t0 = time.perf_counter()
    out = stage._process(ctx)
    wall = time.perf_counter() - t0
    timings = out.stage_results.get(stage.name, {}).get("substep_seconds", {})
    return out.association_results, wall, dict(timings)


def _as_float(text: str) -> float:
    """Parse exactly like Python's float() (pandas' fast parser is not round-trip)."""
    try:
        return float(text)
    except ValueError:
        return float("nan")


def compare(path_a, path_b) -> int:
    """Return 0 if both result tables are identical, else print per-column diffs and return 1."""
    import numpy as np
    import pandas as pd

    # Strings are compared exactly; to_csv writes float repr, so equal text == equal bits.
    a = pd.read_csv(path_a, sep="\t", dtype=str, keep_default_na=False)
    b = pd.read_csv(path_b, sep="\t", dtype=str, keep_default_na=False)
    if list(a.columns) != list(b.columns) or a.shape != b.shape:
        print(f"DIFF shape/columns: {a.shape} {list(a.columns)} vs {b.shape} {list(b.columns)}")
        return 1
    n_diff = 0
    for col in a.columns:
        neq = (a[col] != b[col]).to_numpy()
        if not neq.any():
            continue
        n_diff += 1
        x = np.array([_as_float(v) for v in a[col][neq]])
        y = np.array([_as_float(v) for v in b[col][neq]])
        both = ~(np.isnan(x) | np.isnan(y))
        msg = f"DIFF {col}: {int(neq.sum())} rows"
        if both.any():
            d = np.abs(x[both] - y[both])
            rel = d / np.maximum(np.abs(x[both]), np.finfo(float).tiny)
            msg += f", max_abs={d.max():.3g}, max_rel={rel.max():.3g}"
        if not both.all():
            msg += f", {int((~both).sum())} non-numeric/NaN mismatches"
        print(msg)
    if n_diff:
        return 1
    print(f"IDENTICAL {a.shape[0]} rows x {a.shape[1]} columns")
    return 0


def parse_args(argv=None):
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--size", choices=sorted(PRESETS), default="S")
    p.add_argument("--genes", type=int, help="override the preset gene count")
    p.add_argument("--samples", type=int, help="override the preset sample count")
    p.add_argument("--layout", choices=["packed", "columns"], default="packed")
    p.add_argument("--tests", default="fisher,logistic_burden,skat")
    p.add_argument("--skat-method", default="SKATO")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--blas-threads", type=int, default=1, help="sets OMP/OPENBLAS/MKL threads")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--case-frac", type=float, default=0.1)
    p.add_argument("--unassigned", type=float, default=0.0, help="fraction of samples unassigned")
    p.add_argument("--artifact", action="store_true", help="also write a score artifact")
    p.add_argument("--repeat", type=int, default=1, help="run N times; report the fastest")
    p.add_argument("--out", help="write the association result table (TSV)")
    p.add_argument("--write-baseline", metavar="JSON", help="store this run as a baseline")
    p.add_argument("--tmp", help="scratch directory (default: system temp)")
    p.add_argument("--compare", nargs=2, metavar=("A", "B"), help="compare two result tables")
    return p.parse_args(argv)


def main(argv=None) -> int:
    """Run the benchmark (or ``--compare``) and return the process exit code."""
    args = parse_args(argv)
    if args.compare:
        return compare(*args.compare)
    for var in BLAS_VARS:  # must precede the first numpy import to take effect
        os.environ[var] = str(args.blas_threads)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import logging

    logging.basicConfig(level=logging.ERROR)  # per-gene warnings would drown the table

    genes, samples_n = PRESETS[args.size]
    genes = args.genes or genes
    samples_n = args.samples or samples_n
    t0 = time.perf_counter()
    df, samples = make_df(genes, samples_n, args.layout, args.seed)
    gen_s = time.perf_counter() - t0
    gen_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    best = None
    with tempfile.TemporaryDirectory(prefix="bench_assoc_", dir=args.tmp) as tmp:
        for _ in range(max(1, args.repeat)):
            run = run_stage(df, samples, args, Path(tmp))
            if best is None or run[1] < best[1]:
                best = run
    assert best is not None
    res, wall, timings = best
    rss_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    rss_children = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
    record = {
        "size": args.size,
        "layout": args.layout,
        "genes": genes,
        "samples": samples_n,
        "variants": len(df),
        "tests": args.tests,
        "skat_method": args.skat_method,
        "workers": args.workers,
        "blas_threads": args.blas_threads,
        "repeat": args.repeat,
        "gen_s": round(gen_s, 3),
        "stage_s": round(wall, 3),
        "gen_peak_rss_mb": round(gen_rss),
        "peak_rss_mb": round(rss_self),
        "children_peak_rss_mb": round(rss_children),
        "rows": 0 if res is None else len(res),
        "substeps": {k: round(v, 3) for k, v in timings.items()},
    }
    flat = {k: v for k, v in record.items() if k != "substeps"}
    flat.update({f"sub_{k}": v for k, v in record["substeps"].items()})
    print("RESULT " + " ".join(f"{k}={v}" for k, v in flat.items()))
    print(f"\n{'sub-step':<22}{'seconds':>10}{'% stage':>9}")
    for name, sec in timings.items():
        print(f"{name:<22}{sec:>10.3f}{100 * sec / wall if wall else 0:>8.1f}%")
    print(f"{'stage total':<22}{wall:>10.3f}")
    if args.out and res is not None:
        res.to_csv(args.out, sep="\t", index=False)
    if args.write_baseline:
        record["machine"] = f"{platform.machine()} {platform.processor()} cpus={os.cpu_count()}"
        record["python"] = platform.python_version()
        Path(args.write_baseline).write_text(json.dumps(record, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
