# File: variantcentrifuge/association/tests/fisher.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/fisher.py
"""
Fisher's exact test implementation for the association framework.

Clean reimplementation of the Fisher's exact test from gene_burden.py,
conforming to the AssociationTest ABC. Produces bit-identical p-values and
odds ratios to gene_burden.perform_gene_burden_analysis() for the same
contingency data — verified by the Phase 18 parity test suite.

Statistical approach: 2x2 contingency table Fisher's exact test with:
- Carrier-based (CMC/CAST) or allele-based table construction
- Confidence intervals via statsmodels Table2x2.oddsratio_confint
  (score -> normal -> logit fallback, same as gene_burden.py lines 97-111)
- Haldane-Anscombe continuity correction for zero cells
"""

from __future__ import annotations

import logging
from math import isnan
from typing import Any

import numpy as np

try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None  # type: ignore[assignment]

try:
    from statsmodels.stats.contingency_tables import Table2x2
except ImportError:
    Table2x2 = None  # type: ignore[assignment]

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult

logger = logging.getLogger("variantcentrifuge")


class FisherExactTest(AssociationTest):
    """
    Two-sided Fisher's exact test on a 2x2 contingency table.

    Implements the same statistical pipeline as gene_burden.py:
    1. Table construction (samples mode: carriers vs non-carriers;
       alleles mode: alt alleles vs ref alleles)
    2. scipy.stats.fisher_exact for p-value and odds ratio
    3. statsmodels Table2x2.oddsratio_confint for CIs (score -> normal ->
       logit fallback with continuity correction for zero cells)

    Coupling direction: this module imports from variantcentrifuge.association.
    It does NOT import from gene_burden.py. Bit-identity is guaranteed by
    sharing the same scipy/statsmodels calls with the same parameters.
    """

    @property
    def name(self) -> str:
        """Short identifier used in test registry and output column prefixes."""
        return "fisher"

    def check_dependencies(self) -> None:
        """
        Raise ImportError if scipy or statsmodels are not available.

        Called eagerly by AssociationEngine.from_names() so users get a
        clear error before any data processing begins.
        """
        missing = []
        if fisher_exact is None:
            missing.append("scipy")
        if Table2x2 is None:
            missing.append("statsmodels")
        if missing:
            raise ImportError(
                f"FisherExactTest requires {' and '.join(missing)}. "
                f"Install with: pip install {' '.join(missing)}"
            )

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """
        Run Fisher's exact test for one gene.

        Parameters
        ----------
        gene : str
            Gene symbol being tested.
        contingency_data : dict
            Gene-level aggregated counts with keys:
              proband_count, control_count,
              proband_carrier_count, control_carrier_count,
              proband_allele_count, control_allele_count,
              n_qualifying_variants.
        config : AssociationConfig
            Runtime configuration.

        Returns
        -------
        TestResult
            p_value=None when the gene is skipped (zero qualifying variants
            or both group counts are zero). Otherwise contains p-value, OR,
            and CIs computed via scipy/statsmodels matching gene_burden.py.
        """
        n_variants = int(contingency_data.get("n_qualifying_variants", 0))
        p_count = int(contingency_data.get("proband_count", 0))
        c_count = int(contingency_data.get("control_count", 0))

        # Skip genes with zero qualifying variants (silent, per spec)
        if n_variants == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                n_cases=p_count,
                n_controls=c_count,
                n_variants=n_variants,
            )

        # Skip when both groups are empty (matches gene_burden.py line 425-426)
        if p_count == 0 and c_count == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                n_cases=p_count,
                n_controls=c_count,
                n_variants=n_variants,
            )

        # Build 2x2 table — mirrors gene_burden.py lines 428-465 exactly
        table, skip = self._build_table(gene, contingency_data, config)
        if skip:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                n_cases=p_count,
                n_controls=c_count,
                n_variants=n_variants,
            )

        # Run Fisher's exact test — same call as gene_burden.py line 468
        if fisher_exact is not None:
            odds_ratio, pval = fisher_exact(table)
        else:
            odds_ratio = float("nan")
            pval = 1.0

        # Compute confidence intervals — same logic as gene_burden.py lines 97-114
        ci_lower, ci_upper = self._compute_or_ci(table, config)

        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=float(pval),
            corrected_p_value=None,  # filled in by engine after correction
            effect_size=float(odds_ratio),
            ci_lower=float(ci_lower) if not isnan(float(ci_lower)) else None,
            ci_upper=float(ci_upper) if not isnan(float(ci_upper)) else None,
            n_cases=p_count,
            n_controls=c_count,
            n_variants=n_variants,
            extra={"table": table},
        )

    def _build_table(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> tuple[list[list[int]], bool]:
        """
        Build 2x2 contingency table from gene-level counts.

        Mirrors gene_burden.py lines 428-465. Returns (table, skip) where
        skip=True means the gene should be excluded from results.
        """
        p_count = int(contingency_data.get("proband_count", 0))
        c_count = int(contingency_data.get("control_count", 0))

        if config.gene_burden_mode == "samples":
            # Collapsing test (CMC/CAST): carrier vs non-carrier
            # Uses carrier counts when available, falls back to variant counts
            if "proband_carrier_count" in contingency_data:
                p_var = int(contingency_data["proband_carrier_count"])
                c_var = int(contingency_data["control_carrier_count"])
            else:
                p_var = int(contingency_data.get("proband_variant_count", 0))
                c_var = int(contingency_data.get("control_variant_count", 0))

            p_ref = p_count - p_var
            c_ref = c_count - c_var
            table: list[list[int]] = [[p_var, c_var], [p_ref, c_ref]]

            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference counts: p_count={p_count}, "
                    f"p_var={p_var}, p_ref={p_ref}, c_count={c_count}, "
                    f"c_var={c_var}, c_ref={c_ref}"
                )
                return table, True

        else:
            # Allele-based test: alt alleles vs ref alleles
            p_all = int(contingency_data.get("proband_allele_count", 0))
            c_all = int(contingency_data.get("control_allele_count", 0))
            p_ref = p_count * 2 - p_all
            c_ref = c_count * 2 - c_all
            table = [[p_all, c_all], [p_ref, c_ref]]

            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference allele counts: "
                    f"p_count={p_count}, p_all={p_all}, p_ref={p_ref}, "
                    f"c_count={c_count}, c_all={c_all}, c_ref={c_ref}"
                )
                return table, True

        return table, False

    def _compute_or_ci(
        self,
        table: list[list[int]],
        config: AssociationConfig,
    ) -> tuple[float, float]:
        """
        Compute confidence intervals for the odds ratio.

        Ports _compute_or_confidence_interval() from gene_burden.py
        (lines 48-114) with identical logic: score -> normal -> logit fallback,
        continuity correction for zero cells, structural-zero detection.
        """
        a = table[0][0]
        b = table[0][1]
        c = table[1][0]
        d = table[1][1]

        # Structural zeros: marginal totals are zero (OR undefined)
        if (a + b == 0) or (c + d == 0) or (a + c == 0) or (b + d == 0):
            logger.debug(f"Structural zero in table {table}, cannot compute OR or CI.")
            return np.nan, np.nan

        if Table2x2 is None:
            logger.warning("statsmodels not available. Cannot compute confidence intervals.")
            return np.nan, np.nan

        alpha = config.confidence_interval_alpha
        continuity_correction = config.continuity_correction

        # Apply continuity correction if any cell is zero
        table_np = np.array(table, dtype=float)
        if 0 in table_np.flatten():
            logger.debug(
                f"Applying continuity correction ({continuity_correction}) to table {table}"
            )
            table_for_ci = table_np + continuity_correction
        else:
            table_for_ci = table_np

        try:
            cont_table = Table2x2(table_for_ci)
            # Try score method first — robust for sparse data
            ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="score")

            if isnan(ci_lower) or isnan(ci_upper):
                logger.debug("Score method failed, trying normal approximation.")
                ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="normal")

                if isnan(ci_lower) or isnan(ci_upper):
                    logger.debug("Normal approximation failed, trying logit method.")
                    ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="logit")

            return float(ci_lower), float(ci_upper)
        except Exception as e:
            logger.warning(f"Failed to compute CI for table {table_for_ci.tolist()}: {e}")
            return np.nan, np.nan
