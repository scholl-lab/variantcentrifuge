"""ACAT-O omnibus combination of primary association test p-values per gene."""

from __future__ import annotations

import logging
from typing import Any

from variantcentrifuge.association.base import TestResult

logger = logging.getLogger("variantcentrifuge")


PRIMARY_TESTS = ("acat_o", "skato", "burden")
ACAT_O_EXCLUDABLE = ("acat_v",)
# Discrete tests reported on their own but never combined: ACAT-O is defined over the
# continuous burden / SKAT(-O) / ACAT-V tests (Liu 2019 AJHG; STAAR), and Fisher's p = 1
# mass would veto the Cauchy combination (#216, #210).
NOT_ACAT_O_COMPONENTS = ("fisher",)


def resolve_primary_test(test_names: list[str], config: Any) -> str | None:
    """Return the test id that drives primary_* and FDR, or None for ACAT-O (#210).

    Raises ValueError when the choice is unknown or its test is not being run, so a
    misconfigured run fails before the gene loop.
    """
    primary = getattr(config, "primary_test", "acat_o")
    exclude = getattr(config, "acat_o_exclude", None)
    if exclude not in (None, *ACAT_O_EXCLUDABLE):
        raise ValueError(f"acat_o_exclude must be one of {ACAT_O_EXCLUDABLE}, got {exclude!r}")
    if primary == "acat_o":
        return None
    if primary == "skato":
        if "skat" not in test_names or getattr(config, "skat_method", "SKAT") != "SKATO":
            raise ValueError(
                "--primary-test skato requires 'skat' in --association-tests "
                "and --skat-method SKATO"
            )
        return "skat"
    if primary == "burden":
        trait = getattr(config, "trait_type", "binary")
        name = "linear_burden" if trait == "quantitative" else "logistic_burden"
        if name not in test_names:
            raise ValueError(
                f"--primary-test burden requires '{name}' in --association-tests "
                f"(trait type {trait})"
            )
        return name
    raise ValueError(f"primary_test must be one of {PRIMARY_TESTS}, got {primary!r}")


def compute_acat_o(
    results_by_test: dict[str, dict[str, TestResult]],
    sorted_data: list[dict[str, Any]],
    test_names: list[str],
    exclude: str | None = None,
) -> dict[str, TestResult]:
    """
    Compute ACAT-O omnibus p-value per gene by combining the burden, SKAT/SKAT-O and
    ACAT-V p-values (Fisher excluded, see ``NOT_ACAT_O_COMPONENTS``).

    This is a post-loop meta-test: it runs after all primary tests have completed
    for all genes and combines their p-values via Cauchy combination. ACAT-O is
    NOT in the test registry and does NOT run in the per-gene loop.

    Parameters
    ----------
    results_by_test : dict
        results_by_test[test_name][gene] = TestResult from the per-gene loop.
    sorted_data : list of dict
        Gene burden data in sorted (alphabetical by gene) order.
    test_names : list of str
        Names of the registered primary tests, in engine order.
    exclude : str or None
        ``"acat_v"`` leaves the SKAT test's ACAT-V p-value out of the combination.

    Returns
    -------
    dict mapping gene -> TestResult for acat_o
        Empty dict if no primary tests produced any p-values.
    """
    # Lazy import to avoid circular imports at module load time
    from variantcentrifuge.association.tests.acat import compute_acat_o as _cauchy_combine

    if not test_names:
        logger.warning("_compute_acat_o: no primary tests registered; ACAT-O skipped.")
        return {}

    acat_results: dict[str, TestResult] = {}
    n_valid = 0
    n_total = len(sorted_data)

    for gene_data in sorted_data:
        gene = gene_data.get("GENE", "")

        # Collect p-values from all primary tests for this gene
        test_pvals: dict[str, float | None] = {}
        for test_name in test_names:
            if test_name in NOT_ACAT_O_COMPONENTS:
                continue
            res = results_by_test[test_name].get(gene)
            test_pvals[test_name] = res.p_value if res is not None else None

        # Include ACAT-V if available from SKAT test (Phase 25)
        # ACAT-V is NOT a primary test — it's a per-variant score test stored
        # in the SKAT result's extra dict. Adding it to the ACAT-O omnibus
        # improves power for sparse signals.
        # --acat-o-exclude acat_v: ACAT-V is anti-conservative for ultra-rare variants
        # under case:control imbalance without SPA (#210).
        if exclude != "acat_v":
            for test_name in test_names:
                res = results_by_test[test_name].get(gene)
                if res is not None and "acat_v_p" in res.extra:
                    acat_v_p = res.extra["acat_v_p"]
                    if acat_v_p is not None:
                        test_pvals["acat_v"] = acat_v_p
                    break  # only one SKAT test can produce ACAT-V

        # Combine p-values via Cauchy formula
        acat_p = _cauchy_combine(test_pvals)

        # Find the first primary test result with data for metadata
        first_result = next(
            (
                results_by_test[tn][gene]
                for tn in test_names
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
