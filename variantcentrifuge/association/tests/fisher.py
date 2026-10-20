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
  (normal approximation, same as gene_burden.py)
- Haldane-Anscombe continuity correction for zero cells
"""

from __future__ import annotations

import logging
from math import isnan
from typing import Any

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
    3. statsmodels Table2x2.oddsratio_confint for CIs (normal approximation
       with continuity correction for zero cells)

    Coupling direction: this module imports from variantcentrifuge.association.
    It does NOT import from gene_burden.py. Bit-identity is guaranteed by
    sharing the same scipy/statsmodels calls with the same parameters.
    """

    parallel_safe: bool = True  # Thread-safe -- pure scipy/statsmodels, no R

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
        p_count = int(
            contingency_data.get("proband_called_count", contingency_data.get("proband_count", 0))
        )
        c_count = int(
            contingency_data.get("control_called_count", contingency_data.get("control_count", 0))
        )

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
                se=None,
                n_cases=p_count,
                n_controls=c_count,
                n_variants=n_variants,
            )

        # Skip when either group is empty (matches gene_burden.py line 425-426)
        if p_count == 0 or c_count == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
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
                se=None,
                n_cases=p_count,
                n_controls=c_count,
                n_variants=n_variants,
            )

        # Run Fisher's exact test on original integer counts for p-value
        if fisher_exact is not None:
            _, pval = fisher_exact(table)
        else:
            pval = 1.0

        # Compute OR estimate and confidence intervals (Haldane-Anscombe corrected if any cell is 0)
        # Local import: gene_burden imports the association package at module load.
        from variantcentrifuge.gene_burden import _compute_or_confidence_interval

        odds_ratio, ci_lower, ci_upper = _compute_or_confidence_interval(
            table, config.confidence_interval_alpha, config.continuity_correction
        )

        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=float(pval),
            corrected_p_value=None,  # filled in by engine after correction
            effect_size=float(odds_ratio) if not isnan(float(odds_ratio)) else None,
            ci_lower=float(ci_lower) if not isnan(float(ci_lower)) else None,
            ci_upper=float(ci_upper) if not isnan(float(ci_upper)) else None,
            se=None,
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
        p_count = int(
            contingency_data.get("proband_called_count", contingency_data.get("proband_count", 0))
        )
        c_count = int(
            contingency_data.get("control_called_count", contingency_data.get("control_count", 0))
        )

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
            p_tot_chrom = int(contingency_data.get("proband_called_chromosomes", p_count * 2))
            c_tot_chrom = int(contingency_data.get("control_called_chromosomes", c_count * 2))
            p_ref = p_tot_chrom - p_all
            c_ref = c_tot_chrom - c_all
            table = [[p_all, c_all], [p_ref, c_ref]]

            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference allele counts: "
                    f"p_count={p_count}, p_all={p_all}, p_ref={p_ref}, "
                    f"c_count={c_count}, c_all={c_all}, c_ref={c_ref}"
                )
                return table, True

        return table, False
