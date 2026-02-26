# File: variantcentrifuge/association/engine.py
# Location: variantcentrifuge/variantcentrifuge/association/engine.py
"""
AssociationEngine — orchestrator for multi-test rare variant association.

The engine accepts a list of AssociationTest instances, runs each test on
every gene in the input data, computes ACAT-O omnibus p-values per gene,
applies a single round of multiple testing correction to ACAT-O p-values
(ARCH-03), and returns a wide-format DataFrame with one row per gene.

Column naming convention: {test_name}_{field}, e.g.:
  fisher_pvalue, fisher_or,
  fisher_or_ci_lower, fisher_or_ci_upper

Effect size column names are test-aware (via AssociationTest.effect_column_names()):
  - Fisher: fisher_or, fisher_or_ci_lower, fisher_or_ci_upper
  - Burden tests: logistic_burden_beta, logistic_burden_se,
    logistic_burden_beta_ci_lower, logistic_burden_beta_ci_upper

ACAT-O columns (Phase 22):
  acat_o_pvalue  — raw omnibus Cauchy combination of per-test p-values
  acat_o_qvalue  — FDR/Bonferroni corrected across all genes

SKAT-O note: when skat_method="SKATO", the skat_pvalue column is renamed to
  skat_o_pvalue to distinguish SKAT-O from the plain SKAT result.

FDR strategy (ARCH-03):
  Single correction pass on ACAT-O p-values only. Individual test p-values
  (fisher, burden, skat) remain uncorrected — they exist for diagnostic signal
  decomposition, not independent hypothesis testing. corrected_p_value on
  primary TestResult objects is always None.

Shared columns (gene-level metadata, not test-specific):
  gene, n_cases, n_controls, n_variants
"""

from __future__ import annotations

import logging
import os
from typing import Any

import pandas as pd

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.correction import apply_correction

logger = logging.getLogger("variantcentrifuge")

# Test registry: grows one entry per phase as new tests are implemented.
# Phase 18: Fisher only.
# Phase 19: + "burden" (logistic/linear burden regression)
# Phase 20: + "skat" (R backend)
# Phase 21: + "skat_python" (pure Python backend)
# Phase 22: ACAT-O is a meta-test; NOT in registry (computed post-loop).
# Phase 23: + "coast" (R AllelicSeries COAST allelic series test)
_TEST_REGISTRY: dict[str, type[AssociationTest]] = {}


def _build_registry() -> dict[str, type[AssociationTest]]:
    """Build the test registry lazily to avoid circular imports at module load."""
    from variantcentrifuge.association.tests.allelic_series import COASTTest
    from variantcentrifuge.association.tests.fisher import FisherExactTest
    from variantcentrifuge.association.tests.linear_burden import LinearBurdenTest
    from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest
    from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest
    from variantcentrifuge.association.tests.skat_r import RSKATTest

    return {
        "fisher": FisherExactTest,
        "logistic_burden": LogisticBurdenTest,
        "linear_burden": LinearBurdenTest,
        "skat": RSKATTest,
        "skat_python": PurePythonSKATTest,
        "coast": COASTTest,
    }


def _worker_initializer() -> None:
    """Set BLAS thread counts to 1 in worker processes to prevent oversubscription."""
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"


def _run_gene_worker(
    args: tuple[str, dict, bytes, AssociationConfig],
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
    import pickle

    gene, gene_data, pickled_tests, config = args
    tests: dict[str, AssociationTest] = pickle.loads(pickled_tests)
    results: dict[str, TestResult] = {}
    for test_name, test in tests.items():
        result = test.run(gene, gene_data, config)
        results[test_name] = result
    return gene, results


class AssociationEngine:
    """
    Orchestrates association testing across multiple genes and multiple tests.

    Usage
    -----
    >>> config = AssociationConfig(gene_burden_mode="samples")
    >>> engine = AssociationEngine.from_names(["fisher"], config)
    >>> result_df = engine.run_all(gene_burden_data)

    Parameters
    ----------
    tests : list of AssociationTest
        Instantiated test objects to run. Use from_names() for the common case.
    config : AssociationConfig
        Configuration shared across all tests.
    """

    def __init__(self, tests: list[AssociationTest], config: AssociationConfig) -> None:
        self._tests: dict[str, AssociationTest] = {t.name: t for t in tests}
        self._config = config

    @classmethod
    def from_names(
        cls,
        test_names: list[str],
        config: AssociationConfig,
    ) -> AssociationEngine:
        """
        Construct engine from a list of test name strings.

        Parameters
        ----------
        test_names : list of str
            Names of tests to run (e.g. ["fisher"]).
        config : AssociationConfig
            Runtime configuration.

        Returns
        -------
        AssociationEngine

        Raises
        ------
        ValueError
            If any test_name is not in the registry. The error message lists
            all available tests so users know what's valid. Note: "acat_o" is
            NOT a valid test name here — ACAT-O is computed post-loop as a
            meta-test and does not run per-gene.
        ImportError
            If a test's required dependencies are not installed. Raised eagerly
            (before any data processing) by check_dependencies().
        """
        registry = _build_registry()

        # Backend-aware swap: --skat-backend python/auto routes "skat" to PurePythonSKATTest.
        # This runs BEFORE the unknown-name check so that "skat" resolves correctly.
        # "auto" and "python" follow identical code paths (unified pattern).
        skat_backend = getattr(config, "skat_backend", "python")
        if skat_backend in ("python", "auto"):
            from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

            registry["skat"] = PurePythonSKATTest

        # Backend-aware swap: --coast-backend python/auto routes "coast" to PurePythonCOASTTest.
        # This runs BEFORE the unknown-name check so that "coast" resolves correctly.
        # "auto" and "python" follow identical code paths (unified pattern).
        coast_backend = getattr(config, "coast_backend", "python")
        if coast_backend in ("python", "auto"):
            from variantcentrifuge.association.tests.allelic_series_python import (
                PurePythonCOASTTest,
            )

            registry["coast"] = PurePythonCOASTTest

        available = sorted(registry.keys())

        unknown = [name for name in test_names if name not in registry]
        if unknown:
            for name in unknown:
                raise ValueError(
                    f"Test '{name}' is not available. Available tests: {', '.join(available)}"
                )

        tests: list[AssociationTest] = []
        for name in test_names:
            test_cls = registry[name]
            test = test_cls()
            test.check_dependencies()  # eager dependency check
            tests.append(test)

        return cls(tests, config)

    def _compute_acat_o(
        self,
        results_by_test: dict[str, dict[str, TestResult]],
        sorted_data: list[dict[str, Any]],
    ) -> dict[str, TestResult]:
        """
        Compute ACAT-O omnibus p-value per gene by combining all primary test p-values.

        This is a post-loop meta-test: it runs after all primary tests have completed
        for all genes and combines their p-values via Cauchy combination. ACAT-O is
        NOT in the test registry and does NOT run in the per-gene loop.

        Parameters
        ----------
        results_by_test : dict
            results_by_test[test_name][gene] = TestResult from the per-gene loop.
        sorted_data : list of dict
            Gene burden data in sorted (alphabetical by gene) order.

        Returns
        -------
        dict mapping gene -> TestResult for acat_o
            Empty dict if no primary tests produced any p-values.
        """
        # Lazy import to avoid circular imports at module load time
        from variantcentrifuge.association.tests.acat import compute_acat_o

        if not self._tests:
            logger.warning("_compute_acat_o: no primary tests registered; ACAT-O skipped.")
            return {}

        acat_results: dict[str, TestResult] = {}
        n_valid = 0
        n_total = len(sorted_data)

        for gene_data in sorted_data:
            gene = gene_data.get("GENE", "")

            # Collect p-values from all primary tests for this gene
            test_pvals: dict[str, float | None] = {}
            for test_name in self._tests:
                res = results_by_test[test_name].get(gene)
                test_pvals[test_name] = res.p_value if res is not None else None

            # Include ACAT-V if available from SKAT test (Phase 25)
            # ACAT-V is NOT a primary test — it's a per-variant score test stored
            # in the SKAT result's extra dict. Adding it to the ACAT-O omnibus
            # improves power for sparse signals.
            for test_name in self._tests:
                res = results_by_test[test_name].get(gene)
                if res is not None and "acat_v_p" in res.extra:
                    acat_v_p = res.extra["acat_v_p"]
                    if acat_v_p is not None:
                        test_pvals["acat_v"] = acat_v_p
                    break  # only one SKAT test can produce ACAT-V

            # Combine p-values via Cauchy formula
            acat_p = compute_acat_o(test_pvals)

            # Find the first primary test result with data for metadata
            first_result = next(
                (
                    results_by_test[tn][gene]
                    for tn in self._tests
                    if results_by_test[tn].get(gene) is not None
                ),
                None,
            )
            if first_result is None:
                continue

            acat_result = TestResult(
                gene=gene,
                test_name="acat_o",
                p_value=acat_p,
                corrected_p_value=None,  # populated later by FDR pass
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=first_result.n_cases,
                n_controls=first_result.n_controls,
                n_variants=first_result.n_variants,
            )
            acat_results[gene] = acat_result

            if acat_p is not None:
                n_valid += 1

        logger.info(
            f"ACAT-O: {n_valid}/{n_total} genes with valid omnibus p-values "
            f"({n_total - n_valid} skipped — no valid primary test p-values)"
        )

        if n_valid == 0:
            logger.warning(
                "ACAT-O: no genes had valid p-values. Check that primary tests produced results."
            )

        return acat_results

    def run_all(self, gene_burden_data: list[dict[str, Any]]) -> pd.DataFrame:
        """
        Run all registered tests across all genes and return wide-format results.

        FDR correction (ARCH-03): applied only to ACAT-O p-values across all
        genes. Primary test columns (fisher_pvalue, burden_pvalue, etc.) are
        uncorrected — they are diagnostic signal decomposition, not independent
        hypotheses. corrected_p_value on primary TestResults is always None.

        Parameters
        ----------
        gene_burden_data : list of dict
            One dict per gene with keys: GENE, proband_count, control_count,
            proband_carrier_count, control_carrier_count, proband_allele_count,
            control_allele_count, n_qualifying_variants.

        Returns
        -------
        pd.DataFrame
            Wide-format results. One row per gene that has at least one test
            result (genes where all tests returned p_value=None are excluded).
            Columns: gene, n_cases, n_controls, n_variants, then per-test
            columns: {test}_pvalue (uncorrected), plus test-aware effect
            columns (e.g. fisher_or or logistic_burden_beta), and finally
            acat_o_pvalue and acat_o_qvalue.
        """
        if not gene_burden_data:
            logger.warning("No gene burden data provided to AssociationEngine.")
            return pd.DataFrame()

        # Sort by gene name for deterministic correction order (matches gene_burden.py line 411)
        sorted_data = sorted(gene_burden_data, key=lambda d: d.get("GENE", ""))

        logger.info(
            f"Association analysis: running {list(self._tests.keys())} "
            f"tests on {len(sorted_data)} genes"
        )

        # Collect results: results_by_test[test_name][gene] = TestResult
        results_by_test: dict[str, dict[str, TestResult]] = {name: {} for name in self._tests}

        # Lifecycle hook: prepare() before gene loop (allows progress setup, panel warnings)
        for test in self._tests.values():
            test.prepare(len(sorted_data))

        # Determine whether to use parallel execution (Phase 27)
        n_workers = self._config.association_workers
        all_parallel_safe = all(
            getattr(test, "parallel_safe", False) for test in self._tests.values()
        )
        use_parallel = n_workers != 1 and all_parallel_safe and len(sorted_data) > 1

        if n_workers != 1 and not all_parallel_safe:
            logger.warning(
                "association_workers=%d requested but not all tests are parallel_safe "
                "(falling back to sequential). Non-parallel tests: %s",
                n_workers,
                [n for n, t in self._tests.items() if not getattr(t, "parallel_safe", False)],
            )

        if use_parallel:
            # Run first gene sequentially to trigger lazy null model fitting
            # (SKAT/COAST tests fit null model on first call and cache it).
            first_gene_data = sorted_data[0]
            first_gene = first_gene_data.get("GENE", "")

            # PERF-06: Invoke lazy builder for first gene before sequential run
            if (
                "_genotype_matrix_builder" in first_gene_data
                and "genotype_matrix" not in first_gene_data
            ):
                _builder = first_gene_data.pop("_genotype_matrix_builder")
                _result = _builder()
                first_gene_data["genotype_matrix"] = _result["genotype_matrix"]
                first_gene_data["variant_mafs"] = _result["variant_mafs"]
                first_gene_data["phenotype_vector"] = _result["phenotype_vector"]
                first_gene_data["covariate_matrix"] = _result["covariate_matrix"]
                for _w in _result.get("gt_warnings", []):
                    logger.warning(f"Gene {first_gene}: {_w}")
                if _result.get("mac_filtered"):
                    logger.debug(f"Gene {first_gene}: MAC < 5 — regression will report NA")

            for test_name, test in self._tests.items():
                result = test.run(first_gene, first_gene_data, self._config)
                results_by_test[test_name][first_gene] = result
                logger.debug(
                    f"Gene {first_gene} | {test_name}: p={result.p_value}, OR={result.effect_size}"
                )

            # PERF-06: Discard first gene matrix after tests complete
            first_gene_data.pop("genotype_matrix", None)
            first_gene_data.pop("variant_mafs", None)

            # Pickle test instances now that null models are fitted
            import concurrent.futures
            import pickle

            pickled_tests = pickle.dumps(self._tests)

            # Dispatch remaining genes to worker processes
            remaining = sorted_data[1:]
            if remaining:
                actual_workers = (os.cpu_count() or 1) if n_workers == -1 else n_workers
                # Don't over-provision workers for small remaining panels
                if len(remaining) < actual_workers * 2:
                    actual_workers = max(1, len(remaining) // 2)

                logger.info(
                    f"Parallel association: {actual_workers} workers for "
                    f"{len(remaining)} remaining genes"
                )

                # PERF-06: Build matrices eagerly before pickling for worker dispatch.
                # Builders capture gene_df which would be large to pickle alongside matrix;
                # building here avoids sending both builder (with gene_df) and matrix to workers.
                for gd in remaining:
                    if "_genotype_matrix_builder" in gd and "genotype_matrix" not in gd:
                        _builder = gd.pop("_genotype_matrix_builder")
                        _result = _builder()
                        gd["genotype_matrix"] = _result["genotype_matrix"]
                        gd["variant_mafs"] = _result["variant_mafs"]
                        gd["phenotype_vector"] = _result["phenotype_vector"]
                        gd["covariate_matrix"] = _result["covariate_matrix"]
                        _gene = gd.get("GENE", "")
                        for _w in _result.get("gt_warnings", []):
                            logger.warning(f"Gene {_gene}: {_w}")
                        if _result.get("mac_filtered"):
                            logger.debug(f"Gene {_gene}: MAC < 5 — regression will report NA")

                args_list = [
                    (gd.get("GENE", ""), gd, pickled_tests, self._config) for gd in remaining
                ]
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=actual_workers,
                    initializer=_worker_initializer,
                ) as executor:
                    for gene, gene_results in executor.map(_run_gene_worker, args_list):
                        for test_name, result in gene_results.items():
                            results_by_test[test_name][gene] = result
                            logger.debug(
                                f"Gene {gene} | {test_name}: p={result.p_value}, "
                                f"OR={result.effect_size}"
                            )

                # PERF-06: Discard matrices after parallel batch completes
                for gd in remaining:
                    gd.pop("genotype_matrix", None)
                    gd.pop("variant_mafs", None)
        else:
            # Sequential gene loop (default path or fallback)
            for gene_data in sorted_data:
                gene = gene_data.get("GENE", "")

                # PERF-06: Invoke lazy builder if present (per-gene, O(1) peak memory)
                if "_genotype_matrix_builder" in gene_data and "genotype_matrix" not in gene_data:
                    builder = gene_data.pop("_genotype_matrix_builder")
                    result = builder()
                    gene_data["genotype_matrix"] = result["genotype_matrix"]
                    gene_data["variant_mafs"] = result["variant_mafs"]
                    gene_data["phenotype_vector"] = result["phenotype_vector"]
                    gene_data["covariate_matrix"] = result["covariate_matrix"]
                    for w in result.get("gt_warnings", []):
                        logger.warning(f"Gene {gene}: {w}")
                    if result.get("mac_filtered"):
                        logger.debug(f"Gene {gene}: MAC < 5 — regression will report NA")

                for test_name, test in self._tests.items():
                    result = test.run(gene, gene_data, self._config)
                    results_by_test[test_name][gene] = result
                    logger.debug(
                        f"Gene {gene} | {test_name}: p={result.p_value}, OR={result.effect_size}"
                    )

                # PERF-06: Discard matrix after all tests for this gene complete
                gene_data.pop("genotype_matrix", None)
                gene_data.pop("variant_mafs", None)
                # Keep phenotype_vector and covariate_matrix — shared references, not per-gene

        # Lifecycle hook: finalize() after gene loop (allows timing summary, cleanup)
        for test in self._tests.values():
            test.finalize()

        # --- ARCH-03: Single FDR pass on ACAT-O only ---
        # Primary tests do NOT get corrected_p_value — they are diagnostic decomposition,
        # not independent hypotheses. ACAT-O is the single omnibus significance measure.

        # Step 1: Compute ACAT-O p-values (post-loop meta-test)
        acat_o_results = self._compute_acat_o(results_by_test, sorted_data)
        results_by_test["acat_o"] = acat_o_results

        # Step 2: Apply FDR correction to ACAT-O p-values only
        gene_order = [d.get("GENE", "") for d in sorted_data]
        testable_genes = [
            g
            for g in gene_order
            if acat_o_results.get(g) is not None and acat_o_results[g].p_value is not None
        ]

        # Phase 33: weighted or unweighted correction path
        fdr_weights_by_gene: dict[str, float] = {}
        unweighted_corrected_by_gene: dict[str, float] = {}
        weighted_norm_arr: list[float] = []

        if testable_genes:
            raw_pvals: list[float] = [
                acat_o_results[g].p_value  # type: ignore[misc]
                for g in testable_genes
            ]

            use_weighted = bool(getattr(self._config, "gene_prior_weights", None))

            if use_weighted:
                from variantcentrifuge.association.correction import (
                    apply_weighted_correction,
                    load_gene_weights,
                )

                weight_col = getattr(self._config, "gene_prior_weight_column", "weight")
                weight_map = load_gene_weights(
                    self._config.gene_prior_weights,  # type: ignore[arg-type]
                    weight_col,
                )

                # Unweighted for diagnostics comparison
                unweighted = apply_correction(raw_pvals, self._config.correction_method)
                for g, uw in zip(testable_genes, unweighted, strict=True):
                    unweighted_corrected_by_gene[g] = float(uw)

                corrected, norm_w = apply_weighted_correction(
                    raw_pvals, testable_genes, weight_map, self._config.correction_method
                )
                weighted_norm_arr = norm_w.tolist()
                for g, nw in zip(testable_genes, norm_w, strict=True):
                    fdr_weights_by_gene[g] = float(nw)

                # Write diagnostics if requested (Phase 33)
                diag_out = getattr(self._config, "diagnostics_output", None)
                if diag_out and len(testable_genes) > 1:
                    import numpy as _np

                    from variantcentrifuge.association.diagnostics import (
                        write_fdr_weight_diagnostics,
                    )

                    write_fdr_weight_diagnostics(
                        genes=testable_genes,
                        raw_weights=weight_map,
                        normalized_weights=_np.asarray(weighted_norm_arr, dtype=float),
                        unweighted_corrected=_np.asarray(
                            [unweighted_corrected_by_gene[g] for g in testable_genes],
                            dtype=float,
                        ),
                        weighted_corrected=_np.asarray(
                            [float(corrected[i]) for i in range(len(testable_genes))],
                            dtype=float,
                        ),
                        fdr_threshold=0.05,
                        diagnostics_dir=diag_out,
                    )
            else:
                corrected = apply_correction(raw_pvals, self._config.correction_method)

            for gene, corr_p in zip(testable_genes, corrected, strict=True):
                acat_o_results[gene].corrected_p_value = float(corr_p)

        # Build wide-format DataFrame
        rows = []
        for gene_data in sorted_data:
            gene = gene_data.get("GENE", "")
            # Include gene if it has a result or a skip reason (rank-deficient)
            any_result = any(results_by_test[tn][gene].p_value is not None for tn in self._tests)
            any_skip = any(
                results_by_test[tn][gene].extra.get("skat_skip_reason") for tn in self._tests
            )
            if not any_result and not any_skip:
                continue

            # Use metadata from first result with data, or first for skipped
            first_result = next(
                (r for tn in self._tests if (r := results_by_test[tn][gene]).p_value is not None),
                next(results_by_test[tn][gene] for tn in self._tests),
            )

            row: dict[str, Any] = {
                "gene": gene,
                "n_cases": first_result.n_cases,
                "n_controls": first_result.n_controls,
                "n_variants": first_result.n_variants,
            }

            for test_name, test in self._tests.items():
                res = results_by_test[test_name][gene]
                col_names = test.effect_column_names()
                row[f"{test_name}_pvalue"] = res.p_value
                # Primary test corrected_p_value is always None (ARCH-03)
                # It is written for schema completeness but remains None.
                # None-effect guard: skip column creation when effect/CI names are None.
                # Required for SKAT which has no effect size, SE, or confidence interval.
                # Without this guard, col_names['effect'] = None would produce
                # "skat_None" as a column name.
                if col_names.get("effect") is not None:
                    row[f"{test_name}_{col_names['effect']}"] = res.effect_size
                if col_names.get("se") is not None:
                    row[f"{test_name}_{col_names['se']}"] = res.se
                if col_names.get("ci_lower") is not None:
                    row[f"{test_name}_{col_names['ci_lower']}"] = res.ci_lower
                if col_names.get("ci_upper") is not None:
                    row[f"{test_name}_{col_names['ci_upper']}"] = res.ci_upper
                # Write test-specific extra columns (e.g. skat_o_rho, skat_warnings).
                # Keys in res.extra are already namespaced by the test (e.g. "skat_o_rho"),
                # so we write them without an additional test_name prefix.
                # Fisher and burden tests have empty extra dicts, so this is a no-op for them.
                for extra_key, extra_val in res.extra.items():
                    row[extra_key] = extra_val

                # SKAT-O post-processing: rename skat_pvalue -> skat_o_pvalue when
                # skat_method is SKATO. The skat_method key is stored in res.extra
                # and already written to the row above.
                if res.extra.get("skat_method") == "SKATO" and f"{test_name}_pvalue" in row:
                    row[f"{test_name}_o_pvalue"] = row.pop(f"{test_name}_pvalue")

            # ACAT-O columns (Phase 22) — omnibus significance measure
            acat_res = acat_o_results.get(gene)
            row["acat_o_pvalue"] = acat_res.p_value if acat_res is not None else None
            row["acat_o_qvalue"] = acat_res.corrected_p_value if acat_res is not None else None

            # Phase 33: fdr_weight column — only added when weighted BH is active
            if fdr_weights_by_gene:
                row["fdr_weight"] = fdr_weights_by_gene.get(gene, 1.0)

            rows.append(row)

        if not rows:
            logger.warning("Association analysis: no genes with testable variants.")
            return pd.DataFrame()

        result_df = pd.DataFrame(rows)
        n_sig = (result_df.get("acat_o_qvalue", pd.Series(dtype=float)) < 0.05).sum()
        logger.info(
            f"Association analysis complete: {len(result_df)} genes tested, "
            f"{n_sig} significant (ACAT-O corrected p < 0.05)"
        )
        return result_df
