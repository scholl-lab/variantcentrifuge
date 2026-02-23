# File: variantcentrifuge/association/tests/allelic_series_python.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/allelic_series_python.py
"""
PurePythonCOASTTest -- AssociationTest wrapper for the pure Python COAST backend.

Wraps PythonCOASTBackend to integrate pure-Python COAST into the AssociationEngine
framework. Registered in the engine registry as ``"coast"`` (same name as COASTTest)
so that the backend can be swapped in at engine construction time without changing
downstream consumers.

This class follows the EXACT structural pattern of PurePythonSKATTest (skat_python.py):
- Uses PythonCOASTBackend instead of RSKATBackend/PythonSKATBackend
- Is thread-safe (parallel_safe=True) -- no rpy2 restrictions
- Fits null model lazily on first gene call
- Reuses classify_variants() from allelic_series.py for BMV/DMV/PTV annotation

Variant Classification
-----------------------
PurePythonCOASTTest reuses classify_variants() from allelic_series.py to classify
each gene's variants into:
  - Code 3 (PTV): HIGH impact AND effect in {stop_gained, frameshift_variant, ...}
  - Code 2 (DMV): missense_variant AND (SIFT deleterious OR PolyPhen probably/possibly_damaging)
  - Code 1 (BMV): missense_variant AND SIFT tolerated AND PolyPhen benign
  - Code 0 (unclassified): excluded from COAST

Skip conditions
---------------
The run() method replicates ALL skip-condition guards from COASTTest.run():
  - NO_GENOTYPE_MATRIX: genotype_matrix key absent from contingency_data
  - NO_PHENOTYPE_OR_EMPTY_MATRIX: phenotype_vector is None or geno shape is 0
  - NO_GENE_DF: gene_df key absent from contingency_data
  - MISSING_EFFECT_IMPACT_COLUMNS: EFFECT/IMPACT columns not found in gene_df
  - ANNOTATION_GENOTYPE_MISMATCH: gene_df rows != genotype matrix columns
  - NO_CLASSIFIABLE_VARIANTS: no BMV/DMV/PTV variants after filtering
  - MISSING_CATEGORIES:<cats>: one or more of BMV/DMV/PTV has zero variants

Output extra dict
-----------------
Matches COASTTest output keys for downstream consumers (engine, diagnostics):
  - coast_burden_p_value: Cauchy combination of 6 burden component p-values
  - coast_skat_p_value: allelic SKAT p-value
  - coast_n_bmv, coast_n_dmv, coast_n_ptv: variant counts per category

Thread safety
-------------
PurePythonCOASTTest is thread-safe. The stage using this test can declare
parallel_safe=True (unlike COASTTest which requires parallel_safe=False).

Lifecycle hooks
---------------
prepare(n_genes): logs start message, sets up progress tracking
finalize(): timing summary (no R gc() -- pure Python backend)
"""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any

import numpy as np

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.tests.allelic_series import (
    EFFECT_COLUMN_CANDIDATES,
    IMPACT_COLUMN_CANDIDATES,
    classify_variants,
)

if TYPE_CHECKING:
    from variantcentrifuge.association.backends.base import NullModelResult
    from variantcentrifuge.association.backends.coast_python import PythonCOASTBackend

logger = logging.getLogger("variantcentrifuge")


class PurePythonCOASTTest(AssociationTest):
    """
    COAST allelic series test using the pure Python backend (no R required).

    Registered in the engine as ``"coast"`` -- same name as COASTTest -- so that
    the backend can be swapped in at engine construction time. When R is unavailable
    or ``--coast-backend python`` is specified, the engine uses this class instead
    of COASTTest.

    Parameters
    ----------
    None -- instantiated by AssociationEngine.from_names() via the registry.

    Notes
    -----
    parallel_safe=True: PythonCOASTBackend is thread-safe. No rpy2 restrictions.
    The null model is fit lazily on the first gene call and cached.
    """

    parallel_safe: bool = True  # Thread-safe -- no rpy2

    def __init__(self) -> None:
        self._backend: PythonCOASTBackend | None = None
        self._null_model: NullModelResult | None = None

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
        """Short identifier for this test -- same as COASTTest for registry swap."""
        return "coast"

    def check_dependencies(self) -> None:
        """
        Initialize the Python COAST backend.

        Python COAST always succeeds (numpy/scipy/statsmodels are required
        variantcentrifuge dependencies). This method initializes internal state
        and logs environment info via PythonSKATBackend.log_environment().
        """
        from variantcentrifuge.association.backends.coast_python import PythonCOASTBackend

        self._backend = PythonCOASTBackend()
        # Ensure internal SKAT backend is ready and log environment
        self._backend._ensure_skat_backend()
        if self._backend._skat_backend is not None:
            self._backend._skat_backend.log_environment()

    def effect_column_names(self) -> dict[str, str | None]:
        """
        Column name suffixes for COAST output.

        COAST produces only p-values (omnibus + burden/SKAT components in extra).
        There is no single effect size. All four effect slots are None.

        Returns
        -------
        dict with all values None.
        """
        return {
            "effect": None,
            "se": None,
            "ci_lower": None,
            "ci_upper": None,
        }

    def prepare(self, gene_count: int) -> None:
        """
        Called by the engine before the per-gene loop.

        Sets up progress tracking and logs the start message.

        Parameters
        ----------
        gene_count : int
            Total number of genes that will be processed.
        """
        self._total_genes = gene_count
        self._genes_processed = 0
        self._log_interval = max(10, min(50, gene_count // 10)) if gene_count > 0 else 50
        self._start_time = time.time()
        self._n_complete = 0
        self._n_partial = 0
        self._n_skipped = 0

        logger.info(f"Python COAST: beginning allelic series analysis of {gene_count} genes")

    def finalize(self) -> None:
        """
        Called by the engine after the per-gene loop completes.

        Logs aggregate timing summary. No R gc() (pure Python backend).
        """
        elapsed = time.time() - self._start_time
        n_tested = self._n_complete + self._n_partial
        logger.info(f"Python COAST complete: {self._genes_processed} genes in {elapsed:.1f}s")
        logger.info(
            f"COAST: {n_tested} genes tested "
            f"({self._n_complete} complete, {self._n_partial} partial, "
            f"{self._n_skipped} skipped)"
        )

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """
        Run COAST for a single gene using the pure Python backend.

        Replicates all skip-condition guards from COASTTest.run(), then delegates
        to PythonCOASTBackend.test_gene() for the statistical computation.

        Parameters
        ----------
        gene : str
            Gene symbol being tested.
        contingency_data : dict
            Must contain:
              - ``genotype_matrix``  : np.ndarray (n_samples, n_variants)
              - ``phenotype_vector`` : np.ndarray (n_samples,)
              - ``gene_df``          : pd.DataFrame with per-variant annotation rows
            Optional:
              - ``covariate_matrix`` : np.ndarray (n_samples, k) or None
        config : AssociationConfig
            Runtime configuration including coast_weights and trait_type.

        Returns
        -------
        TestResult
            p_value=None when test is skipped (missing variant categories,
            no genotype matrix, or insufficient data).
            extra contains: coast_burden_p_value, coast_skat_p_value,
            coast_n_bmv, coast_n_dmv, coast_n_ptv.
        """
        n_cases = int(contingency_data.get("proband_count", 0))
        n_controls = int(contingency_data.get("control_count", 0))
        n_variants = int(contingency_data.get("n_qualifying_variants", 0))

        # ── Guard: genotype matrix required ──────────────────────────────────────
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

        # ── Guard: phenotype and non-empty genotype matrix ────────────────────────
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

        # ── Guard: gene_df required for annotation ────────────────────────────────
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

        # ── Find EFFECT and IMPACT columns ────────────────────────────────────────
        effect_col = _find_column(gene_df, EFFECT_COLUMN_CANDIDATES, "EFFECT")
        impact_col = _find_column(gene_df, IMPACT_COLUMN_CANDIDATES, "IMPACT")

        if effect_col is None or impact_col is None:
            logger.warning(
                f"Python COAST [{gene}]: could not find EFFECT ({effect_col}) or "
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

        # ── Classify variants into COAST categories ───────────────────────────────
        anno_codes, include_mask = classify_variants(gene_df, effect_col, impact_col)

        # ── Guard: annotation/genotype dimension alignment ────────────────────────
        if len(anno_codes) != geno.shape[1]:
            logger.warning(
                f"Python COAST [{gene}]: annotation codes length ({len(anno_codes)}) "
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

        # ── Filter to COAST-eligible variants ─────────────────────────────────────
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

        # ── Count variants per category ───────────────────────────────────────────
        n_bmv = int(np.sum(anno_filtered == 1))
        n_dmv = int(np.sum(anno_filtered == 2))
        n_ptv = int(np.sum(anno_filtered == 3))

        # ── Partial-category logic: skip only when ALL categories are empty ───────
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
                f"Python COAST [{gene}]: partial -- missing {', '.join(missing)} "
                f"(BMV={n_bmv}, DMV={n_dmv}, PTV={n_ptv})"
            )

        # Backend must have been initialized by check_dependencies()
        if self._backend is None:
            raise RuntimeError(
                "PurePythonCOASTTest.run() called without check_dependencies(). "
                "Use AssociationEngine.from_names() which calls check_dependencies() "
                "at construction time."
            )

        # ── Lazy null model: fit once, reuse for all genes ────────────────────────
        if self._null_model is None:
            self._backend._ensure_skat_backend()
            assert self._backend._skat_backend is not None
            self._null_model = self._backend._skat_backend.fit_null_model(
                phenotype=phenotype,
                covariates=covariate_matrix,
                trait_type=config.trait_type,
            )

        # ── COAST weights from config or default ──────────────────────────────────
        coast_weights = getattr(config, "coast_weights", None) or [1.0, 2.0, 3.0]

        # ── Run full COAST via backend ─────────────────────────────────────────────
        try:
            result = self._backend.test_gene(
                gene=gene,
                geno_filtered=geno_filtered,
                anno_codes_filtered=anno_filtered,
                phenotype=phenotype,
                covariates=covariate_matrix,
                coast_weights=coast_weights,
                trait_type=config.trait_type,
                null_model=self._null_model,
            )
        except Exception as exc:
            logger.error(f"Python COAST [{gene}]: backend test_gene() failed: {exc}")
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
                    "coast_skip_reason": f"BACKEND_ERROR:{exc!s}",
                    "coast_n_bmv": n_bmv,
                    "coast_n_dmv": n_dmv,
                    "coast_n_ptv": n_ptv,
                },
            )

        # ── Increment counter and progress logging ────────────────────────────────
        self._genes_processed += 1
        if coast_status == "complete":
            self._n_complete += 1
        else:
            self._n_partial += 1
        if self._total_genes > 0 and self._genes_processed % self._log_interval == 0:
            pct = 100.0 * self._genes_processed / self._total_genes
            logger.info(
                f"Python COAST progress: "
                f"{self._genes_processed}/{self._total_genes} genes ({pct:.0f}%)"
            )

        # ── Compute coast_burden_p_value as Cauchy of 6 burden components ─────────
        # This mirrors the convention expected by downstream consumers (engine, diagnostics).
        # The omnibus p_value already combines all 7 (6 burden + 1 SKAT); the extra
        # coast_burden_p_value provides a standalone summary of the burden sub-test.
        from variantcentrifuge.association.tests.acat import cauchy_combination

        burden_p_values: list[float | None] = result.get("burden_p_values", [])
        coast_burden_p = cauchy_combination(burden_p_values) if burden_p_values else None

        extra_base: dict[str, Any] = {
            "coast_burden_p_value": coast_burden_p,
            "coast_skat_p_value": result.get("skat_p_value"),
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
            p_value=result.get("p_value"),
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


def _find_column(gene_df: Any, candidates: list[str], label: str) -> str | None:
    """Find the first available column name from a list of candidates."""
    cols = set(gene_df.columns)
    for candidate in candidates:
        if candidate in cols:
            return candidate
    logger.debug(
        f"PurePythonCOASTTest: could not find {label} column. "
        f"Tried: {candidates}. Available: {sorted(cols)}"
    )
    return None
