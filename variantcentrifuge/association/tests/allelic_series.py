# File: variantcentrifuge/association/tests/allelic_series.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/allelic_series.py
"""
COASTTest — AssociationTest wrapper for the R AllelicSeries::COAST() test.

Implements the Combined Omnibus Association Test (COAST) for allelic series
testing. COAST has higher power than standard burden/SKAT when effect sizes
are ordered by variant damage severity: BMV (benign missense) < DMV (damaging
missense) < PTV (protein-truncating variants). This ordered alternative
provides substantial power gains for variants with monotone effect size trends.

Reference: McCaw et al. AJHG 2023. "Operating within the context of an allelic
series, COAST maintains type 1 error control and demonstrates superior power
compared to burden and SKAT tests."

Classification logic is implemented in ``allelic_series_classify.py`` and
re-exported here for backward compatibility.
R bridge functions are implemented in ``allelic_series_r.py``.

Thread safety
--------------
COASTTest is not thread-safe (rpy2 restriction). The stage that uses COASTTest
must declare parallel_safe=False.

Lifecycle hooks
----------------
prepare(n_genes): logs start message, sets up progress tracking
finalize(): timing summary, R gc() cleanup
"""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any

import numpy as np

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.tests.allelic_series_classify import (
    CADD_COLUMN_CANDIDATES,
    EFFECT_COLUMN_CANDIDATES,
    IMPACT_COLUMN_CANDIDATES,
    MISSENSE_EFFECT,
    POLYPHEN_BENIGN,
    POLYPHEN_COLUMN_CANDIDATES,
    POLYPHEN_DAMAGING,
    PTV_EFFECTS,
    SIFT_BENIGN,
    SIFT_COLUMN_CANDIDATES,
    SIFT_DAMAGING,
    _resolve_effect,
    classify_variants,
)
from variantcentrifuge.association.tests.allelic_series_r import (
    _GC_INTERVAL,
    _extract_pval,
    check_r_dependencies,
    run_r_coast,
)

if TYPE_CHECKING:
    pass

logger = logging.getLogger("variantcentrifuge")

__all__ = [
    "CADD_COLUMN_CANDIDATES",
    "EFFECT_COLUMN_CANDIDATES",
    "IMPACT_COLUMN_CANDIDATES",
    "MISSENSE_EFFECT",
    "POLYPHEN_BENIGN",
    "POLYPHEN_COLUMN_CANDIDATES",
    "POLYPHEN_DAMAGING",
    "PTV_EFFECTS",
    "SIFT_BENIGN",
    "SIFT_COLUMN_CANDIDATES",
    "SIFT_DAMAGING",
    "_GC_INTERVAL",
    "COASTTest",
    "_extract_pval",
    "_resolve_effect",
    "check_r_dependencies",
    "classify_variants",
    "run_r_coast",
]


class COASTTest(AssociationTest):
    """
    # DEPRECATED
    COAST allelic series test via R AllelicSeries::COAST().

    Registered in the engine as ``"coast"``. Classifies variants into
    BMV/DMV/PTV categories using SnpEff EFFECT/IMPACT and dbNSFP SIFT/PolyPhen
    annotations, then invokes the R AllelicSeries COAST omnibus test.

    COAST combines a burden and a SKAT component that both test the allelic
    series hypothesis. The omnibus p-value from COAST feeds into ACAT-O.
    The burden and SKAT component p-values are stored in extra for diagnostics.

    Parameters
    ----------
    None — instantiated by AssociationEngine.from_names() via the registry.

    Raises
    ------
    ImportError
        From check_dependencies() if rpy2 is not importable or the
        AllelicSeries R package is not installed.
    """

    parallel_safe: bool = False  # rpy2 restriction: main thread only

    def __init__(self) -> None:
        import warnings

        warnings.warn(
            "COASTTest (R COAST backend) is deprecated and will be removed in v0.17.0. "
            "Use --coast-backend python (default).",
            DeprecationWarning,
            stacklevel=2,
        )
        # Lifecycle tracking (set by prepare())
        self._total_genes: int = 0
        self._log_interval: int = 50
        self._genes_processed: int = 0
        self._start_time: float = 0.0

        # Category-level counters (set by prepare(), incremented in run())
        self._n_complete: int = 0
        self._n_partial: int = 0
        self._n_skipped: int = 0

    @property
    def name(self) -> str:
        """Short identifier for this test used in registry lookup and column prefixes."""
        return "coast"

    def check_dependencies(self) -> None:
        """Verify that rpy2 is importable and the R AllelicSeries package is installed."""
        check_r_dependencies()

    def effect_column_names(self) -> dict[str, str | None]:
        """Column name suffixes for COAST output."""
        return {
            "effect": None,
            "se": None,
            "ci_lower": None,
            "ci_upper": None,
        }

    def prepare(self, gene_count: int) -> None:
        """Called by the engine before the per-gene loop."""
        self._total_genes = gene_count
        self._genes_processed = 0
        self._log_interval = max(10, min(50, gene_count // 10)) if gene_count > 0 else 50
        self._start_time = time.time()
        self._n_complete = 0
        self._n_partial = 0
        self._n_skipped = 0

        logger.info(f"COAST: beginning allelic series analysis of {gene_count} genes")

    def finalize(self) -> None:
        """Called by the engine after the per-gene loop completes."""
        elapsed = time.time() - self._start_time
        n_tested = self._n_complete + self._n_partial
        logger.info(f"COAST complete: {self._genes_processed} genes in {elapsed:.1f}s")
        logger.info(
            f"COAST: {n_tested} genes tested "
            f"({self._n_complete} complete, {self._n_partial} partial, "
            f"{self._n_skipped} skipped)"
        )

        # Final R gc() to free heap
        try:
            import rpy2.robjects as ro

            ro.r("gc()")
        except Exception:
            pass  # R cleanup is best-effort

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """Run COAST for a single gene."""
        n_cases = int(contingency_data.get("proband_count", 0))
        n_controls = int(contingency_data.get("control_count", 0))
        n_variants = int(contingency_data.get("n_qualifying_variants", 0))

        # Require genotype matrix and phenotype
        if "genotype_matrix" not in contingency_data:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_GENOTYPE_MATRIX"},
            )

        geno: np.ndarray = contingency_data["genotype_matrix"]
        phenotype: np.ndarray | None = contingency_data.get("phenotype_vector")
        covariate_matrix: np.ndarray | None = contingency_data.get("covariate_matrix")
        gene_df = contingency_data.get("gene_df")

        if phenotype is None or geno.shape[0] == 0 or geno.shape[1] == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_PHENOTYPE_OR_EMPTY_MATRIX"},
            )

        if gene_df is None or len(gene_df) == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_GENE_DF"},
            )

        # Find EFFECT and IMPACT columns in gene_df
        effect_col = self._find_column(gene_df, EFFECT_COLUMN_CANDIDATES, "EFFECT")
        impact_col = self._find_column(gene_df, IMPACT_COLUMN_CANDIDATES, "IMPACT")

        if effect_col is None or impact_col is None:
            logger.warning(
                f"COAST [{gene}]: could not find EFFECT ({effect_col}) or "
                f"IMPACT ({impact_col}) columns. Skipping."
            )
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "MISSING_EFFECT_IMPACT_COLUMNS"},
            )

        # Classify variants into COAST categories
        coast_model_dir = getattr(config, "coast_classification", None)
        anno_codes, include_mask = classify_variants(
            gene_df, effect_col, impact_col, model_dir=coast_model_dir
        )

        keep = contingency_data.get("variant_keep_mask")  # matrix builder site QC (#215)
        if keep is not None and len(keep) == len(anno_codes):
            anno_codes, include_mask = anno_codes[keep], include_mask[keep]
        if len(anno_codes) != geno.shape[1]:
            logger.warning(
                f"COAST [{gene}]: annotation codes length ({len(anno_codes)}) "
                f"does not match genotype matrix columns ({geno.shape[1]}). "
                "This may happen when site filtering removed variants. Skipping."
            )
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "ANNOTATION_GENOTYPE_MISMATCH"},
            )

        # Filter to COAST-eligible variants only
        anno_filtered = anno_codes[include_mask]
        geno_filtered = geno[:, include_mask]

        if geno_filtered.shape[1] == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_CLASSIFIABLE_VARIANTS"},
            )

        # Count variants per category
        n_bmv = int(np.sum(anno_filtered == 1))
        n_dmv = int(np.sum(anno_filtered == 2))
        n_ptv = int(np.sum(anno_filtered == 3))

        missing = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n == 0]
        present = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n > 0]

        if len(present) == 0:
            self._n_skipped += 1
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={
                    "coast_skip_reason": "ALL_CATEGORIES_EMPTY",
                    "coast_n_bmv": 0,
                    "coast_n_dmv": 0,
                    "coast_n_ptv": 0,
                    "coast_status": "skipped",
                },
            )

        coast_status = "complete" if not missing else "partial"
        if missing:
            logger.debug(
                f"COAST [{gene}]: partial -- missing {', '.join(missing)} "
                f"(BMV={n_bmv}, DMV={n_dmv}, PTV={n_ptv})"
            )

        # Run R AllelicSeries::COAST()
        try:
            burden_p, skat_p, omnibus_p = self._run_r_coast(
                geno_filtered=geno_filtered,
                anno_codes_filtered=anno_filtered,
                phenotype=phenotype,
                covariate_matrix=covariate_matrix,
                config=config,
            )
        except Exception as exc:
            logger.error(
                f"COAST [{gene}]: R AllelicSeries::COAST() failed: {type(exc).__name__}: {exc!r}"
            )
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={
                    "coast_skip_reason": f"R_ERROR:{exc!s}",
                    "coast_n_bmv": n_bmv,
                    "coast_n_dmv": n_dmv,
                    "coast_n_ptv": n_ptv,
                },
            )

        # Increment counter and periodic GC
        self._genes_processed += 1
        if coast_status == "complete":
            self._n_complete += 1
        else:
            self._n_partial += 1
        if self._genes_processed % _GC_INTERVAL == 0:
            try:
                import rpy2.robjects as ro

                ro.r("gc()")
            except Exception:
                pass

        # Progress logging
        if self._total_genes > 0 and self._genes_processed % self._log_interval == 0:
            pct = 100.0 * self._genes_processed / self._total_genes
            logger.info(
                f"COAST progress: {self._genes_processed}/{self._total_genes} genes ({pct:.0f}%)"
            )

        extra_base: dict[str, Any] = {
            "coast_burden_pvalue": burden_p,
            "coast_skat_pvalue": skat_p,
            "coast_n_bmv": n_bmv,
            "coast_n_dmv": n_dmv,
            "coast_n_ptv": n_ptv,
            "coast_status": coast_status,
        }
        if missing:
            extra_base["coast_missing_categories"] = ",".join(missing)

        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=omnibus_p,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=n_cases,
            n_controls=n_controls,
            n_variants=n_variants,
            extra=extra_base,
        )

    def _find_column(self, gene_df: Any, candidates: list[str], label: str) -> str | None:
        """Find the first available column name from a list of candidates."""
        cols = set(gene_df.columns)
        for candidate in candidates:
            if candidate in cols:
                return candidate
        logger.debug(
            f"COASTTest: could not find {label} column. "
            f"Tried: {candidates}. Available: {sorted(cols)}"
        )
        return None

    def _run_r_coast(
        self,
        geno_filtered: np.ndarray,
        anno_codes_filtered: np.ndarray,
        phenotype: np.ndarray,
        covariate_matrix: np.ndarray | None,
        config: AssociationConfig,
    ) -> tuple[float | None, float | None, float | None]:
        """Delegate R AllelicSeries::COAST() execution to allelic_series_r module."""
        return run_r_coast(
            geno_filtered=geno_filtered,
            anno_codes_filtered=anno_codes_filtered,
            phenotype=phenotype,
            covariate_matrix=covariate_matrix,
            config=config,
        )
