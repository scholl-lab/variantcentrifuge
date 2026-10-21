# File: variantcentrifuge/association/engine_workers.py
# Location: variantcentrifuge/association/engine_workers.py
"""Worker execution routines and bounded dispatch for association tests."""

from __future__ import annotations

import concurrent.futures
import contextlib
import ctypes
import gc
import logging
import multiprocessing
import os
import sys
from collections.abc import Callable, Iterator, Sequence
from typing import TYPE_CHECKING, Any

from variantcentrifuge.association.timing import Progress

if TYPE_CHECKING:
    from variantcentrifuge.association.base import (
        AssociationConfig,
        AssociationTest,
        TestResult,
    )

logger = logging.getLogger("variantcentrifuge")


# (task, read-only state) set by the parent just before forking: workers inherit the state
# copy-on-write, so only index ranges and small results cross the process boundary (#212).
_FORK_STATE: tuple[Callable[[Any, int, int], Any], Any] | None = None

# (setter, getter) of OpenBLAS (plain, scipy-openblas32/64 wheels) and MKL.
_BLAS_THREAD_API = (
    ("openblas_set_num_threads", "openblas_get_num_threads"),
    ("openblas_set_num_threads64_", "openblas_get_num_threads64_"),
    ("scipy_openblas_set_num_threads", "scipy_openblas_get_num_threads"),
    ("scipy_openblas_set_num_threads64_", "scipy_openblas_get_num_threads64_"),
    ("MKL_Set_Num_Threads", "MKL_Get_Max_Threads"),
)


def _loaded_blas_libraries() -> list[str]:
    """Paths of the BLAS shared libraries mapped into this process (Linux only).

    ponytail: matches OpenBLAS and the single-library MKL (mkl_rt) only; split MKL
    builds (libmkl_core + threading layer) stay unpinned. Add their names if one shows up.
    """
    try:
        with open("/proc/self/maps") as fh:
            paths = {line.split()[-1] for line in fh if "/" in line}
    except OSError:
        return []
    return sorted(
        p for p in paths if ".so" in p and ("openblas" in p.lower() or "mkl_rt" in p.lower())
    )


def _blas_thread_controls() -> list[tuple[Any, Any]]:
    """(setter, getter) for every loaded BLAS library that exports a known pair."""
    controls = []
    for path in _loaded_blas_libraries():
        lib = ctypes.CDLL(path)  # already mapped: returns the loaded library
        for set_name, get_name in _BLAS_THREAD_API:
            if hasattr(lib, set_name) and hasattr(lib, get_name):
                controls.append((getattr(lib, set_name), getattr(lib, get_name)))
    return controls


@contextlib.contextmanager
def single_blas_thread() -> Iterator[None]:
    """Run per-gene tests with one BLAS thread, restoring the previous counts afterwards.

    Per-gene matrices are small, so multithreaded BLAS is several times slower there;
    parallelism comes from worker processes instead (#212).
    """
    controls = _blas_thread_controls()
    before = [get() for _set, get in controls]
    for set_threads, _get in controls:
        set_threads(1)
    try:
        yield
    finally:
        for (set_threads, _get), n in zip(controls, before, strict=True):
            set_threads(n)


def _worker_initializer() -> None:
    """Limit BLAS to one thread in a worker process.

    A forked worker inherits BLAS libraries already initialised by the parent, so the
    environment variables alone have no effect there; call each library's setter.
    """
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    for set_threads, _get in _blas_thread_controls():
        set_threads(1)


_BUILDER_RESULT_KEYS = (
    "genotype_matrix",
    "variant_mafs",
    "phenotype_vector",
    "covariate_matrix",
    "score_values",
    "variant_weight_column",
    "cadd_scores",
    "revel_scores",
    "variant_effects",
    "variant_keep_mask",
    "n_variants_qc",  # kept after the matrix is discarded: the stage reports it (#200)
)

_PER_GENE_MATRIX_PAYLOAD_KEYS = (
    "_genotype_matrix_builder",
    "genotype_matrix",
    "variant_mafs",
    "score_values",
    "variant_weight_column",
    "cadd_scores",
    "revel_scores",
    "variant_effects",
    "variant_keep_mask",
)


def _apply_builder_result_to_gene_data(
    gene_data: dict[str, Any],
    result: dict[str, Any],
    gene: str,
) -> None:
    for key in _BUILDER_RESULT_KEYS:
        if key in result:
            gene_data[key] = result[key]
    for warning in result.get("gt_warnings", []):
        logger.warning(f"Gene {gene}: {warning}")
    if result.get("mac_filtered"):
        logger.debug(f"Gene {gene}: MAC < 5 — regression will report NA")


def _discard_per_gene_matrix_payload(gene_data: dict[str, Any]) -> None:
    for key in _PER_GENE_MATRIX_PAYLOAD_KEYS:
        gene_data.pop(key, None)


def _single_test_output_id(test_name: str, result: TestResult | None) -> str:
    """Return the user-facing column prefix for a single-test primary result."""
    if test_name == "skat" and result is not None and result.extra.get("skat_method") == "SKATO":
        return "skat_o"
    return test_name


# Called right after a gene's matrix is built, before its tests: (gene, builder, build result).
GeneHook = Callable[[str, Any, dict[str, Any]], Any]


def run_sequential_genes(
    gene_data_list: list[dict[str, Any]],
    tests: dict[str, AssociationTest],
    config: AssociationConfig,
    results_by_test: dict[str, dict[str, TestResult]],
    gene_hook: GeneHook | None = None,
    hook_results: dict[str, Any] | None = None,
    progress: Progress | None = None,
) -> None:
    """Run tests sequentially across a list of genes with per-gene matrix lifecycle.

    ``gene_hook``'s return value for each gene whose matrix is built here goes to
    ``hook_results[gene]`` (the score artifact uses it to reuse the tests' matrix, #212).
    """
    for i, gene_data in enumerate(gene_data_list, 1):
        gene = gene_data.get("GENE", "")
        try:
            # PERF-06: Invoke lazy builder if present (per-gene, O(1) peak memory)
            if "_genotype_matrix_builder" in gene_data and "genotype_matrix" not in gene_data:
                builder = gene_data.pop("_genotype_matrix_builder")
                result = builder()
                if gene_hook is not None and hook_results is not None:
                    hook_results[gene] = gene_hook(gene, builder, result)
                _apply_builder_result_to_gene_data(gene_data, result, gene)

            for test_name, test in tests.items():
                result = test.run(gene, gene_data, config)
                results_by_test[test_name][gene] = result
                logger.debug(
                    f"Gene {gene} | {test_name}: p={result.p_value}, OR={result.effect_size}"
                )
        finally:
            # PERF-06: Discard matrix after all tests for this gene complete
            _discard_per_gene_matrix_payload(gene_data)
            if progress is not None:
                progress.update(i)


def _fork_task(bounds: tuple[int, int]) -> Any:
    assert _FORK_STATE is not None, "worker started without inherited fork state"
    task, state = _FORK_STATE
    return task(state, *bounds)


def _task_bounds(
    n_items: int, workers: int, costs: Sequence[float] | None = None
) -> list[tuple[int, int]]:
    """``fork_map``'s ``(lo, hi)`` tasks in submission order.

    Without ``costs``: ~8 chunks per worker in item order. With ``costs`` (e.g. variants
    per gene): tasks of at most 8 items; tasks costing over twice the mean go first (costliest
    first), so a huge gene starts early instead of leaving one worker busy at the end
    (audit 0.25 fix 5).
    """
    chunk = max(1, n_items // (workers * 8))
    if costs is not None:
        chunk = min(chunk, 8)
    bounds = [(lo, min(lo + chunk, n_items)) for lo in range(0, n_items, chunk)]
    if costs is not None and bounds:
        task_costs = {b: sum(costs[b[0] : b[1]]) for b in bounds}
        big = 2 * sum(task_costs.values()) / len(bounds)
        # Only costly tasks jump the queue; the rest keep item order so the in-order
        # results (progress, early errors) keep flowing during the run (review P3-r1).
        bounds.sort(key=lambda b: -task_costs[b] if task_costs[b] > big else 0)
    return bounds


def fork_map(
    task: Callable[[Any, int, int], Any],
    state: Any,
    n_items: int,
    workers: int,
    costs: Sequence[float] | None = None,
) -> Iterator[tuple[int, int, Any]]:
    """Yield ``(lo, hi, task(state, lo, hi))`` over ``range(n_items)`` in item order.

    Runs in forked worker processes (one BLAS thread each) that inherit ``state``; ``task``
    must be a module-level function. Tasks are planned by ``_task_bounds`` (``costs``:
    small tasks, costly ones first); results are still yielded in item order.
    Linux only: callers check ``fork_available()`` first.
    """
    global _FORK_STATE
    bounds = _task_bounds(n_items, workers, costs)
    _FORK_STATE = (task, state)
    # Frozen objects are skipped by the children's GC, so its passes don't dirty (copy) the
    # parent's heap pages that the workers share copy-on-write. Skipped when the caller
    # already froze objects: unfreeze() would undo the caller's freeze too.
    freeze = gc.get_freeze_count() == 0
    if freeze:
        gc.collect()
        gc.freeze()
    try:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers,
            mp_context=multiprocessing.get_context("fork"),
            initializer=_worker_initializer,
        ) as executor:
            futures = {b: executor.submit(_fork_task, b) for b in bounds}
            try:
                for lo, hi in sorted(futures):
                    yield lo, hi, futures[lo, hi].result()
            finally:  # on error or early close, like executor.map: drop tasks not started
                for future in futures.values():
                    future.cancel()
    finally:
        _FORK_STATE = None
        if freeze:
            gc.unfreeze()


def available_cpus() -> int:
    """CPUs this process may use (SLURM/cgroup affinity), what every -1 worker count means."""
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def fork_available() -> bool:
    """Whether fork_map can run (fork start method with copy-on-write: Linux)."""
    return sys.platform.startswith("linux")


def _run_gene_range(state: tuple, lo: int, hi: int) -> tuple[dict, dict, dict]:
    """Worker task: build matrices and run all tests (and the gene hook) for genes[lo:hi]."""
    genes, tests, config, gene_hook = state
    chunk = genes[lo:hi]
    results: dict[str, dict[str, TestResult]] = {name: {} for name in tests}
    hook_results: dict[str, Any] = {}
    run_sequential_genes(chunk, tests, config, results, gene_hook, hook_results)
    qc = {gd.get("GENE", ""): gd["n_variants_qc"] for gd in chunk if "n_variants_qc" in gd}
    return results, qc, hook_results


def dispatch_parallel_genes(
    remaining: list[dict[str, Any]],
    tests: dict[str, AssociationTest],
    config: AssociationConfig,
    actual_workers: int,
    results_by_test: dict[str, dict[str, TestResult]],
    gene_hook: GeneHook | None = None,
    hook_results: dict[str, Any] | None = None,
    progress: Progress | None = None,
) -> None:
    """Run genes in forked workers (see ``fork_map``); results match the sequential loop.

    Workers inherit the genes, the tests (null models already fitted by the parent),
    the config and ``gene_hook``, call the lazy matrix builders themselves and return
    only small results. Without fork (non-Linux) genes run sequentially.
    """
    if not fork_available():
        logger.info("Parallel association needs fork (Linux); running genes sequentially")
        run_sequential_genes(
            remaining, tests, config, results_by_test, gene_hook, hook_results,
            progress or Progress("Association genes", len(remaining)),
        )  # fmt: skip
        return
    logger.info(f"Parallel association: {actual_workers} workers for {len(remaining)} genes")
    state = (remaining, tests, config, gene_hook)
    progress = progress or Progress("Association genes", len(remaining))
    for lo, hi, (chunk_results, qc, chunk_hooks) in fork_map(
        _run_gene_range,
        state,
        len(remaining),
        actual_workers,
        [gd.get("n_qualifying_variants", 1) for gd in remaining],
    ):
        for test_name, by_gene in chunk_results.items():
            results_by_test[test_name].update(by_gene)
        if hook_results is not None:
            hook_results.update(chunk_hooks)
        for gd in remaining[lo:hi]:
            _discard_per_gene_matrix_payload(gd)  # built in the worker, not here
            gene = gd.get("GENE", "")
            if gene in qc:
                gd["n_variants_qc"] = qc[gene]
        progress.update(hi)
