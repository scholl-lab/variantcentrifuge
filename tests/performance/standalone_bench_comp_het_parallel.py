#!/usr/bin/env python3
"""
Standalone synthetic benchmark for compound het parallelization.

Measures wall-time for the full analyze_inheritance_parallel() call
across multiple min_variants_for_parallel thresholds and worker counts.

Usage:
    python tests/performance/standalone_bench_comp_het_parallel.py

No assertions — exit 0 always. Human reviews the timing table.
"""

import os
import statistics
import sys
import time

# Allow running from project root without installation
_project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

import numpy as np  # noqa: I001
import pandas as pd


# ---------------------------------------------------------------------------
# Synthetic data generation
# ---------------------------------------------------------------------------


def _generate_skewed_genes(
    n_genes: int,
    n_variants: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Generate a gene assignment array with a realistic skewed distribution.

    A few "large" genes (like PKD1) get 50-200 variants each.
    Most genes get 2-5 variants each.
    """
    gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]

    # Assign variant counts: 2% of genes get 50-200 variants (the "large" genes)
    n_large = max(1, int(n_genes * 0.02))
    n_small = n_genes - n_large

    # Large gene counts: uniform 50-200
    large_counts = rng.integers(50, 201, size=n_large)
    # Small gene counts: uniform 2-5
    small_counts = rng.integers(2, 6, size=n_small)

    all_counts = np.concatenate([large_counts, small_counts])
    # Scale to match n_variants
    total = all_counts.sum()
    scaled = np.round(all_counts * n_variants / total).astype(int)
    # Adjust last gene to hit target exactly
    diff = n_variants - scaled.sum()
    scaled[-1] = max(2, scaled[-1] + diff)

    # Build gene assignment array
    gene_array = np.repeat(gene_names, scaled)
    # Shuffle genes (not variants within gene — keep grouped for realistic groupby)
    rng.shuffle(gene_array)
    return gene_array[:n_variants]


def generate_benchmark_data(
    n_genes: int = 1000,
    n_variants: int = 10_000,
    n_samples: int = 6,
    seed: int = 42,
) -> tuple[pd.DataFrame, dict, list[str]]:
    """
    Generate synthetic variant DataFrame + pedigree for compound het benchmarking.

    Returns
    -------
    (df, pedigree_data, sample_list)
        df has columns: CHROM, POS, REF, ALT, GENE, and per-sample GT columns
        pedigree_data has sample_id -> {sample_id, father_id, mother_id, sex, affected_status}
        sample_list is a list of sample IDs
    """
    rng = np.random.default_rng(seed)

    # Build pedigree: n_trios = n_samples // 3 trios
    n_trios = max(1, n_samples // 3)
    pedigree_data = {}
    sample_list = []

    for t in range(n_trios):
        child_id = f"CHILD_{t + 1:03d}"
        father_id = f"FATHER_{t + 1:03d}"
        mother_id = f"MOTHER_{t + 1:03d}"

        pedigree_data[child_id] = {
            "sample_id": child_id,
            "father_id": father_id,
            "mother_id": mother_id,
            "sex": "1" if rng.random() < 0.5 else "2",
            "affected_status": "2",  # affected
        }
        pedigree_data[father_id] = {
            "sample_id": father_id,
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",  # unaffected
        }
        pedigree_data[mother_id] = {
            "sample_id": mother_id,
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",  # unaffected
        }
        sample_list.extend([child_id, father_id, mother_id])

    # Use only n_samples samples (trim if necessary)
    sample_list = sample_list[:n_samples]

    # Generate chromosomes with realistic distribution
    chrom_choices = [str(i) for i in range(1, 23)] + ["X"]
    chrom_weights = np.array([0.14] + [0.075] * 9 + [0.015] * 12 + [0.005])
    chrom_weights = chrom_weights / chrom_weights.sum()
    chroms = rng.choice(chrom_choices, size=n_variants, p=chrom_weights)

    # Generate positions
    positions = rng.integers(1, 250_000_000, size=n_variants)

    # Generate REF/ALT
    bases = ["A", "C", "G", "T"]
    refs = rng.choice(bases, size=n_variants)
    alts = np.array([rng.choice([b for b in bases if b != r]) for r in refs])

    # Generate genes with skewed distribution
    genes = _generate_skewed_genes(n_genes, n_variants, rng)

    # Generate per-sample genotypes
    # Realistic mix: ~60% 0/0, ~25% 0/1, ~10% 1/1, ~5% ./.
    gt_choices = ["0/0", "0/1", "1/1", "./."]
    gt_weights = [0.60, 0.25, 0.10, 0.05]
    sample_gts = {
        sample_id: rng.choice(gt_choices, size=n_variants, p=gt_weights)
        for sample_id in sample_list
    }

    # Build DataFrame
    data = {
        "CHROM": chroms,
        "POS": positions,
        "REF": refs,
        "ALT": alts,
        "GENE": genes,
    }
    for sample_id, gts in sample_gts.items():
        data[sample_id] = gts

    df = pd.DataFrame(data)

    return df, pedigree_data, sample_list


# ---------------------------------------------------------------------------
# Benchmark runner
# ---------------------------------------------------------------------------


def _time_call(fn, n_runs: int = 3) -> tuple[float, float, float]:
    """Run fn n_runs times, return (median, min, max) wall-times in seconds."""
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), min(times), max(times)


def run_benchmark(
    df: pd.DataFrame,
    pedigree_data: dict,
    sample_list: list[str],
    n_runs: int = 3,
) -> list[dict]:
    """
    Run benchmark sweeping min_variants_for_parallel thresholds and worker counts.

    Returns list of result dicts with keys:
        label, min_variants_for_parallel, n_workers, median_s, min_s, max_s
    """
    from variantcentrifuge.inheritance.parallel_analyzer import analyze_inheritance_parallel

    results = []

    # Sequential baseline: force sequential by setting threshold very high
    print("  Measuring sequential baseline ...", flush=True)

    def sequential():
        return analyze_inheritance_parallel(
            df.copy(),
            dict(pedigree_data),
            list(sample_list),
            min_variants_for_parallel=999_999_999,
            n_workers=1,
        )

    med, mn, mx = _time_call(sequential, n_runs)
    results.append(
        {
            "label": "sequential",
            "min_variants_for_parallel": "N/A",
            "n_workers": 1,
            "median_s": med,
            "min_s": mn,
            "max_s": mx,
        }
    )

    # Determine available CPUs
    cpu_count = os.cpu_count() or 2
    worker_counts = sorted({1, 2, 4} | ({cpu_count} if cpu_count > 4 else set()))

    thresholds = [25, 50, 100, 200, 500]

    for threshold in thresholds:
        for n_workers in worker_counts:
            label = f"parallel (workers={n_workers}, threshold={threshold})"
            print(f"  Measuring {label} ...", flush=True)

            def parallel(t=threshold, w=n_workers):
                return analyze_inheritance_parallel(
                    df.copy(),
                    dict(pedigree_data),
                    list(sample_list),
                    min_variants_for_parallel=t,
                    n_workers=w,
                )

            med, mn, mx = _time_call(parallel, n_runs)
            results.append(
                {
                    "label": label,
                    "min_variants_for_parallel": threshold,
                    "n_workers": n_workers,
                    "median_s": med,
                    "min_s": mn,
                    "max_s": mx,
                }
            )

    return results


def _print_table(
    results: list[dict],
    n_genes: int,
    n_variants: int,
    n_samples: int,
    n_trios: int,
) -> None:
    """Print benchmark results as a formatted table."""
    print()
    print("Compound Het Pass 2 Benchmark")
    print("=" * 60)
    print(f"Dataset: {n_genes} genes, {n_variants} variants, {n_samples} samples ({n_trios} trios)")
    print()

    col_label = 44
    col_num = 10

    header = (
        f"{'Configuration':<{col_label}} | "
        f"{'Median (s)':>{col_num}} | "
        f"{'Min (s)':>{col_num}} | "
        f"{'Max (s)':>{col_num}}"
    )
    separator = "-" * len(header)
    print(header)
    print(separator)

    sequential_median = None
    for r in results:
        if r["label"] == "sequential":
            sequential_median = r["median_s"]
        speedup = ""
        if sequential_median and r["label"] != "sequential" and r["median_s"] > 0:
            ratio = sequential_median / r["median_s"]
            speedup = f"  ({ratio:.2f}x)"

        label = r["label"] + speedup
        print(
            f"{label:<{col_label + 12}} | "
            f"{r['median_s']:>{col_num}.3f} | "
            f"{r['min_s']:>{col_num}.3f} | "
            f"{r['max_s']:>{col_num}.3f}"
        )

    print()
    print("Note: speedup relative to sequential baseline.")
    print()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    n_genes = 1000
    n_variants = 10_000
    n_samples = 6
    n_runs = 3
    seed = 42

    n_trios = max(1, n_samples // 3)

    print(
        f"Generating synthetic dataset: {n_genes} genes, {n_variants} variants, "
        f"{n_samples} samples ({n_trios} trios) ...",
        flush=True,
    )
    t0 = time.perf_counter()
    df, pedigree_data, sample_list = generate_benchmark_data(
        n_genes=n_genes,
        n_variants=n_variants,
        n_samples=n_samples,
        seed=seed,
    )
    print(f"Data generation: {time.perf_counter() - t0:.2f}s", flush=True)
    print(f"Actual genes: {df['GENE'].nunique()}, actual variants: {len(df)}", flush=True)
    print(f"Sample list: {sample_list}", flush=True)
    print(f"Running each config {n_runs} times ...", flush=True)
    print()

    results = run_benchmark(df, pedigree_data, sample_list, n_runs=n_runs)

    _print_table(results, n_genes, n_variants, n_samples, n_trios)


if __name__ == "__main__":
    main()
