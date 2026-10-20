"""
AssociationEngine — orchestrator for multi-test rare variant association.

Runs registered association tests per gene, determines primary p-values,
applies multiple testing correction, and returns wide-format DataFrames.
"""

from __future__ import annotations

import logging
import os
from typing import Any

import pandas as pd

from variantcentrifuge.association.acat_o import compute_acat_o
from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.correction import apply_correction
from variantcentrifuge.association.engine_workers import (
    _BUILDER_RESULT_KEYS,
    _PER_GENE_MATRIX_PAYLOAD_KEYS,
    _apply_builder_result_to_gene_data,
    _discard_per_gene_matrix_payload,
    _single_test_output_id,
    _worker_initializer,
    dispatch_parallel_genes,
    run_sequential_genes,
    single_blas_thread,
)

__all__ = [
    "_BUILDER_RESULT_KEYS",
    "_PER_GENE_MATRIX_PAYLOAD_KEYS",
    "AssociationEngine",
    "_apply_builder_result_to_gene_data",
    "_build_registry",
    "_discard_per_gene_matrix_payload",
    "_single_test_output_id",
    "_worker_initializer",
    "dispatch_parallel_genes",
    "run_sequential_genes",
]

logger = logging.getLogger("variantcentrifuge")

# Below this many genes, worker start-up costs more than the tests themselves.
MIN_GENES_FOR_PARALLEL = 50


def _null_model_pending(tests: dict[str, AssociationTest]) -> bool:
    """True while a test with a lazy null model has neither fitted nor failed to fit it."""
    return any(
        getattr(t, "_null_model", False) is None and getattr(t, "_null_failed_key", None) is None
        for t in tests.values()
    )


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

        if skat_backend == "r" and "skat" in test_names:
            from variantcentrifuge.association.tests.skat_r import _is_r_skat_weight_supported

            if not _is_r_skat_weight_supported(config.variant_weights):
                raise ValueError(
                    "R SKAT backend supports only beta:a,b and uniform variant weights. "
                    f"Got '{config.variant_weights}'. Use --skat-backend python for cadd, revel, "
                    "combined, score_column, or column:<name> weights."
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
        """Compute ACAT-O per gene; see ``association.acat_o.compute_acat_o``."""
        return compute_acat_o(results_by_test, sorted_data, list(self._tests))

    def run_all(self, gene_burden_data: list[dict[str, Any]]) -> pd.DataFrame:
        """
        Run all registered tests across all genes and return wide-format results.

        FDR correction is applied once to the primary p-value across all genes.
        Multi-test and SKAT runs use ACAT-O as the primary result. Single
        non-SKAT runs use the test's own p-value as the primary result.

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
            columns: {test}_pvalue, plus test-aware effect columns (e.g.
            fisher_or or logistic_burden_beta), and primary_pvalue/primary_qvalue.
            ACAT-O columns are present only when ACAT-O is the primary result.
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
        use_parallel = (
            n_workers != 1 and all_parallel_safe and len(sorted_data) >= MIN_GENES_FOR_PARALLEL
        )

        if n_workers != 1 and not all_parallel_safe:
            logger.warning(
                "association_workers=%d requested but not all tests are parallel_safe "
                "(falling back to sequential). Non-parallel tests: %s",
                n_workers,
                [n for n, t in self._tests.items() if not getattr(t, "parallel_safe", False)],
            )

        with single_blas_thread():  # gene-level parallelism only, never BLAS (#212)
            if use_parallel:
                # Run genes in the parent until every lazy null model is fitted (or failed), so
                # forked workers inherit it instead of each refitting it (#212).
                done = 0
                while done < len(sorted_data) and (
                    done == 0
                    or (done < MIN_GENES_FOR_PARALLEL and _null_model_pending(self._tests))
                ):
                    run_sequential_genes(
                        sorted_data[done : done + 1], self._tests, self._config, results_by_test
                    )
                    done += 1

                # Dispatch remaining genes to worker processes
                remaining = sorted_data[done:]
                if remaining:
                    actual_workers = (os.cpu_count() or 1) if n_workers == -1 else n_workers
                    # Don't over-provision workers for small remaining panels
                    if len(remaining) < actual_workers * 2:
                        actual_workers = max(1, len(remaining) // 2)

                    dispatch_parallel_genes(
                        remaining=remaining,
                        tests=self._tests,
                        config=self._config,
                        actual_workers=actual_workers,
                        results_by_test=results_by_test,
                    )
            else:
                # Sequential gene loop (default path or fallback)
                run_sequential_genes(
                    sorted_data,
                    self._tests,
                    self._config,
                    results_by_test,
                )

        # Lifecycle hook: finalize() after gene loop (allows timing summary, cleanup)
        for test in self._tests.values():
            test.finalize()

        if not self._tests:
            logger.warning("Association analysis: no association tests registered.")
            return pd.DataFrame()

        # Single correction pass on the run's primary p-value.
        #
        # ACAT-O is the right primary result for multi-test runs. It also stays
        # primary for SKAT because the engine can include ACAT-V in that omnibus
        # combination. For single non-SKAT tests, ACAT-O would be a pure
        # pass-through, so correction is attached to the actual test result
        # (for example fisher_qvalue) instead of misleading acat_o_* columns.
        primary_result_id = next(iter(self._tests))
        write_acat_o = len(self._tests) > 1 or "skat" in self._tests
        acat_o_results: dict[str, TestResult] = {}
        correction_results: dict[str, TestResult]
        if write_acat_o:
            acat_o_results = self._compute_acat_o(results_by_test, sorted_data)
            results_by_test["acat_o"] = acat_o_results
            correction_results = acat_o_results
            primary_result_id = "acat_o"
        else:
            correction_results = results_by_test[primary_result_id]

        gene_order = [d.get("GENE", "") for d in sorted_data]
        testable_genes = [
            g
            for g in gene_order
            if correction_results.get(g) is not None and correction_results[g].p_value is not None
        ]

        # Phase 33: weighted or unweighted correction path
        fdr_weights_by_gene: dict[str, float] = {}
        unweighted_corrected_by_gene: dict[str, float] = {}
        weighted_corrected_by_gene: dict[str, float] = {}
        weighted_norm_arr: list[float] = []

        if testable_genes:
            raw_pvals: list[float] = [
                correction_results[g].p_value  # type: ignore[misc]
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
                for g, wt in zip(testable_genes, corrected, strict=True):
                    weighted_corrected_by_gene[g] = float(wt)

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
                correction_results[gene].corrected_p_value = float(corr_p)

        # Build wide-format DataFrame
        rows = []
        for gene_data in sorted_data:
            gene = gene_data.get("GENE", "")
            # Include gene if it has a result, a skip reason or a warning code
            any_result = any(results_by_test[tn][gene].p_value is not None for tn in self._tests)
            any_skip = any(
                (extra := results_by_test[tn][gene].extra).get("skat_skip_reason")
                or extra.get("coast_skip_reason")
                or extra.get("warnings")
                for tn in self._tests
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
                # Corrected significance is emitted through primary_* columns
                # and, for single-test runs, through {test}_qvalue aliases.
                # Do not reintroduce legacy {test}_corrected_pvalue columns.
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

            # Primary significance columns. In single-test pass-through mode,
            # also expose a test-specific q-value such as fisher_qvalue.
            primary_res = correction_results.get(gene)
            primary_output_id = primary_result_id
            if not write_acat_o:
                primary_output_id = _single_test_output_id(primary_result_id, primary_res)
                row[f"{primary_output_id}_qvalue"] = (
                    primary_res.corrected_p_value if primary_res is not None else None
                )
                if weighted_corrected_by_gene:
                    row[f"{primary_output_id}_weighted_qvalue"] = weighted_corrected_by_gene.get(
                        gene
                    )
                    row[f"{primary_output_id}_unweighted_qvalue"] = (
                        unweighted_corrected_by_gene.get(gene)
                    )

            row["primary_test"] = primary_output_id
            row["primary_pvalue"] = primary_res.p_value if primary_res is not None else None
            row["primary_qvalue"] = (
                primary_res.corrected_p_value if primary_res is not None else None
            )
            if weighted_corrected_by_gene:
                row["primary_weighted_qvalue"] = weighted_corrected_by_gene.get(gene)
                row["primary_unweighted_qvalue"] = unweighted_corrected_by_gene.get(gene)

            # ACAT-O columns are only emitted when ACAT-O is a real primary
            # result, not for single-test pass-throughs such as Fisher-only.
            if write_acat_o:
                acat_res = acat_o_results.get(gene)
                row["acat_o_pvalue"] = acat_res.p_value if acat_res is not None else None
                row["acat_o_qvalue"] = acat_res.corrected_p_value if acat_res is not None else None
                if weighted_corrected_by_gene:
                    row["acat_o_weighted_qvalue"] = weighted_corrected_by_gene.get(gene)
                    row["acat_o_unweighted_qvalue"] = unweighted_corrected_by_gene.get(gene)

            # Phase 33: fdr_weight column — only added when weighted BH is active
            if fdr_weights_by_gene:
                row["fdr_weight"] = fdr_weights_by_gene.get(gene, 1.0)

            rows.append(row)

        if not rows:
            logger.warning("Association analysis: no genes with testable variants.")
            return pd.DataFrame()

        result_df = pd.DataFrame(rows)
        n_sig = (result_df.get("primary_qvalue", pd.Series(dtype=float)) < 0.05).sum()
        logger.info(
            f"Association analysis complete: {len(result_df)} genes tested, "
            f"{n_sig} significant (primary corrected p < 0.05)"
        )
        return result_df
