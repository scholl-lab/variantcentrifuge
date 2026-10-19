# File: variantcentrifuge/association/engine_workers.py
# Location: variantcentrifuge/association/engine_workers.py
"""Worker execution routines and bounded dispatch for association tests."""

from __future__ import annotations

import concurrent.futures
import logging
import os
import pickle
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from variantcentrifuge.association.base import (
        AssociationConfig,
        AssociationTest,
        TestResult,
    )

logger = logging.getLogger("variantcentrifuge")


def _worker_initializer() -> None:
    """Set BLAS thread counts to 1 in worker processes to prevent oversubscription."""
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"


def _run_gene_worker(
    args: tuple[str, dict[str, Any], bytes, AssociationConfig],
) -> tuple[str, dict[str, TestResult]]:
    """Process a single gene in a subprocess worker.

    Called by ProcessPoolExecutor. Receives pre-pickled test instances
    (with null models already fitted) to avoid redundant null model fitting
    in each worker.

    Parameters
    ----------
    args : tuple
        (gene_name, gene_data_dict, pickled_tests_bytes, config)

    Returns
    -------
    tuple
        (gene_name, dict mapping test_name -> TestResult)
    """
    gene, gene_data, pickled_tests, config = args
    tests: dict[str, AssociationTest] = pickle.loads(pickled_tests)
    results: dict[str, TestResult] = {}
    for test_name, test in tests.items():
        result = test.run(gene, gene_data, config)
        results[test_name] = result
    return gene, results


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
)

_PER_GENE_MATRIX_PAYLOAD_KEYS = (
    "genotype_matrix",
    "variant_mafs",
    "score_values",
    "variant_weight_column",
    "cadd_scores",
    "revel_scores",
    "variant_effects",
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


def run_sequential_genes(
    gene_data_list: list[dict[str, Any]],
    tests: dict[str, AssociationTest],
    config: AssociationConfig,
    results_by_test: dict[str, dict[str, TestResult]],
) -> None:
    """Run tests sequentially across a list of genes with per-gene matrix lifecycle."""
    for gene_data in gene_data_list:
        gene = gene_data.get("GENE", "")
        try:
            # PERF-06: Invoke lazy builder if present (per-gene, O(1) peak memory)
            if "_genotype_matrix_builder" in gene_data and "genotype_matrix" not in gene_data:
                builder = gene_data.pop("_genotype_matrix_builder")
                result = builder()
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


def dispatch_parallel_genes(
    remaining: list[dict[str, Any]],
    tests: dict[str, AssociationTest],
    config: AssociationConfig,
    actual_workers: int,
    results_by_test: dict[str, dict[str, TestResult]],
) -> None:
    """Dispatch genes across worker processes in bounded batches.

    To prevent memory exhaustion on large panels, matrices are built
    in bounded batches proportional to worker count (2 * actual_workers),
    and each batch's matrix payloads are discarded immediately after its
    results are collected.
    """
    logger.info(
        f"Parallel association: {actual_workers} workers for {len(remaining)} remaining genes"
    )

    pickled_tests = pickle.dumps(tests)
    batch_size = max(1, actual_workers * 2)

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=actual_workers,
        initializer=_worker_initializer,
    ) as executor:
        for batch_start in range(0, len(remaining), batch_size):
            batch = remaining[batch_start : batch_start + batch_size]
            try:
                # Build matrices lazily only for this batch
                for gd in batch:
                    if "_genotype_matrix_builder" in gd and "genotype_matrix" not in gd:
                        _builder = gd.pop("_genotype_matrix_builder")
                        _result = _builder()
                        _gene = gd.get("GENE", "")
                        _apply_builder_result_to_gene_data(gd, _result, _gene)

                args_list = [(gd.get("GENE", ""), gd, pickled_tests, config) for gd in batch]

                for gene, gene_results in executor.map(_run_gene_worker, args_list):
                    for test_name, result in gene_results.items():
                        results_by_test[test_name][gene] = result
                        logger.debug(
                            f"Gene {gene} | {test_name}: "
                            f"p={result.p_value}, OR={result.effect_size}"
                        )
            finally:
                # Discard matrix payloads immediately for this batch
                for gd in batch:
                    _discard_per_gene_matrix_payload(gd)
