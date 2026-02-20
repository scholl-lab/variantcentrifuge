# File: variantcentrifuge/association/engine.py
# Location: variantcentrifuge/variantcentrifuge/association/engine.py
"""
AssociationEngine — orchestrator for multi-test rare variant association.

The engine accepts a list of AssociationTest instances, runs each test on
every gene in the input data, applies multiple testing correction, and
returns a wide-format DataFrame with one row per gene.

Column naming convention: {test_name}_{field}, e.g.:
  fisher_p_value, fisher_corrected_p_value, fisher_or,
  fisher_or_ci_lower, fisher_or_ci_upper

Effect size column names are test-aware (via AssociationTest.effect_column_names()):
  - Fisher: fisher_or, fisher_or_ci_lower, fisher_or_ci_upper
  - Burden tests: logistic_burden_beta, logistic_burden_se,
    logistic_burden_beta_ci_lower, logistic_burden_beta_ci_upper

Shared columns (gene-level metadata, not test-specific):
  gene, n_cases, n_controls, n_variants
"""

from __future__ import annotations

import logging
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
# Phase 22: + "acat_o" (omnibus)
_TEST_REGISTRY: dict[str, type[AssociationTest]] = {}


def _build_registry() -> dict[str, type[AssociationTest]]:
    """Build the test registry lazily to avoid circular imports at module load."""
    from variantcentrifuge.association.tests.fisher import FisherExactTest
    from variantcentrifuge.association.tests.linear_burden import LinearBurdenTest
    from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest
    from variantcentrifuge.association.tests.skat_r import RSKATTest

    return {
        "fisher": FisherExactTest,
        "logistic_burden": LogisticBurdenTest,
        "linear_burden": LinearBurdenTest,
        "skat": RSKATTest,
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
            all available tests so users know what's valid.
        ImportError
            If a test's required dependencies are not installed. Raised eagerly
            (before any data processing) by check_dependencies().
        """
        registry = _build_registry()
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

    def run_all(self, gene_burden_data: list[dict[str, Any]]) -> pd.DataFrame:
        """
        Run all registered tests across all genes and return wide-format results.

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
            columns: {test}_p_value, {test}_corrected_p_value, plus
            test-aware effect columns (e.g. fisher_or or logistic_burden_beta).
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

        for gene_data in sorted_data:
            gene = gene_data.get("GENE", "")
            for test_name, test in self._tests.items():
                result = test.run(gene, gene_data, self._config)
                results_by_test[test_name][gene] = result
                logger.debug(
                    f"Gene {gene} | {test_name}: p={result.p_value}, OR={result.effect_size}"
                )

        # Lifecycle hook: finalize() after gene loop (allows timing summary, cleanup)
        for test in self._tests.values():
            test.finalize()

        # Apply multiple testing correction per test (on genes with non-None p-values)
        for test_name in self._tests:
            gene_results = results_by_test[test_name]
            # Collect testable genes in sorted order (same order as correction in gene_burden.py)
            testable_genes = [
                g
                for g in (d.get("GENE", "") for d in sorted_data)
                if gene_results[g].p_value is not None
            ]
            if testable_genes:
                # p_value is not None for all testable_genes (filtered above)
                raw_pvals: list[float] = [
                    gene_results[g].p_value  # type: ignore[misc]
                    for g in testable_genes
                ]
                corrected = apply_correction(raw_pvals, self._config.correction_method)
                for gene, corr_p in zip(testable_genes, corrected, strict=True):
                    gene_results[gene].corrected_p_value = float(corr_p)

        # Build wide-format DataFrame
        rows = []
        for gene_data in sorted_data:
            gene = gene_data.get("GENE", "")
            # Include gene if at least one test has a real result
            any_result = any(results_by_test[tn][gene].p_value is not None for tn in self._tests)
            if not any_result:
                continue

            # Use metadata from the first test result that has data
            first_result = next(
                r for tn in self._tests if (r := results_by_test[tn][gene]).p_value is not None
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
                row[f"{test_name}_p_value"] = res.p_value
                row[f"{test_name}_corrected_p_value"] = res.corrected_p_value
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

            rows.append(row)

        if not rows:
            logger.warning("Association analysis: no genes with testable variants.")
            return pd.DataFrame()

        result_df = pd.DataFrame(rows)
        n_sig = (result_df.get("fisher_corrected_p_value", pd.Series(dtype=float)) < 0.05).sum()
        logger.info(
            f"Association analysis complete: {len(result_df)} genes tested, "
            f"{n_sig} significant (corrected p < 0.05)"
        )
        return result_df
