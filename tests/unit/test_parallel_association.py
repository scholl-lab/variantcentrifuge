# File: tests/unit/test_parallel_association.py
# Location: tests/unit/test_parallel_association.py
"""
Tests for parallel gene-level association execution (Phase 27 Plan 03).

Verifies that:
- Parallel execution (association_workers > 1) produces identical results
  to sequential execution (association_workers = 1) when all tests are
  parallel_safe.
- Non-parallel-safe tests trigger a sequential fallback with a warning.
- Single-gene data stays sequential even if workers > 1.
- _worker_initializer sets BLAS environment variables in worker processes.
"""

from __future__ import annotations

import os
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.engine import AssociationEngine, _worker_initializer

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_gene_data(n: int = 5) -> list[dict[str, Any]]:
    """Generate synthetic gene burden data for n genes."""
    genes = []
    for i in range(n):
        genes.append(
            {
                "GENE": f"GENE{i}",
                "proband_count": 100,
                "control_count": 100,
                "proband_carrier_count": 10 + i,
                "control_carrier_count": 5 + i,
                "proband_allele_count": 12 + i,
                "control_allele_count": 6 + i,
                "n_qualifying_variants": 3 + i,
            }
        )
    return genes


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_parallel_results_match_sequential() -> None:
    """Parallel execution with Fisher test produces identical DataFrame to sequential."""
    genes = _make_gene_data(5)

    config_seq = AssociationConfig(association_workers=1)
    engine_seq = AssociationEngine.from_names(["fisher"], config_seq)
    df_seq = engine_seq.run_all(genes)

    config_par = AssociationConfig(association_workers=2)
    engine_par = AssociationEngine.from_names(["fisher"], config_par)
    df_par = engine_par.run_all(genes)

    # Both DataFrames should have the same shape
    assert df_seq.shape == df_par.shape, (
        f"Sequential shape {df_seq.shape} != parallel shape {df_par.shape}"
    )

    # Sort both by gene to ensure consistent comparison order
    df_seq_sorted = df_seq.sort_values("gene").reset_index(drop=True)
    df_par_sorted = df_par.sort_values("gene").reset_index(drop=True)

    # Assert DataFrames are equal (same genes, same p-values, same effect sizes)
    pd.testing.assert_frame_equal(
        df_seq_sorted,
        df_par_sorted,
        check_exact=False,
        rtol=1e-10,  # Allow tiny floating-point differences from process boundary
    )


@pytest.mark.unit
def test_parallel_fallback_on_unsafe_test() -> None:
    """Non-parallel-safe test triggers sequential fallback with warning."""

    class UnsafeTest(AssociationTest):
        """Mock test with parallel_safe=False."""

        parallel_safe: bool = False

        @property
        def name(self) -> str:
            return "unsafe_test"

        def run(
            self, gene: str, contingency_data: dict[str, Any], config: AssociationConfig
        ) -> TestResult:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=0.5,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=contingency_data.get("proband_count", 0),
                n_controls=contingency_data.get("control_count", 0),
                n_variants=contingency_data.get("n_qualifying_variants", 0),
            )

        def effect_column_names(self) -> dict[str, str | None]:
            return {"effect": None, "se": None, "ci_lower": None, "ci_upper": None}

    genes = _make_gene_data(3)
    config = AssociationConfig(association_workers=2)
    engine = AssociationEngine([UnsafeTest()], config)

    # Capture the warning log message
    with patch("variantcentrifuge.association.engine.logger") as mock_logger:
        df = engine.run_all(genes)
        # Verify warning was emitted
        warning_calls = [str(call) for call in mock_logger.warning.call_args_list]
        fallback_warned = any("falling back to sequential" in call for call in warning_calls)
        assert fallback_warned, (
            f"Expected 'falling back to sequential' warning, got: {warning_calls}"
        )

    # Still produces valid output (sequential fallback works)
    assert len(df) == 3
    assert "unsafe_test_p_value" in df.columns


@pytest.mark.unit
def test_parallel_single_gene_stays_sequential() -> None:
    """With only 1 gene, parallel path is not taken even if workers > 1."""
    genes = _make_gene_data(1)

    # workers=2 but only 1 gene — use_parallel = (workers!=1 and parallel_safe and len>1)
    # len(sorted_data) > 1 is False, so use_parallel = False
    config = AssociationConfig(association_workers=2)
    engine = AssociationEngine.from_names(["fisher"], config)

    # Should complete without error and produce a single-row result
    df = engine.run_all(genes)
    assert len(df) == 1
    assert df["gene"].iloc[0] == "GENE0"
    assert "fisher_p_value" in df.columns
    # p_value should be a valid float (not None)
    assert df["fisher_p_value"].iloc[0] is not None


@pytest.mark.unit
def test_worker_initializer_sets_env() -> None:
    """_worker_initializer sets OPENBLAS_NUM_THREADS and MKL_NUM_THREADS to 1."""
    # Save original values
    orig_openblas = os.environ.get("OPENBLAS_NUM_THREADS")
    orig_mkl = os.environ.get("MKL_NUM_THREADS")
    orig_omp = os.environ.get("OMP_NUM_THREADS")

    try:
        _worker_initializer()
        assert os.environ["OPENBLAS_NUM_THREADS"] == "1"
        assert os.environ["MKL_NUM_THREADS"] == "1"
        assert os.environ["OMP_NUM_THREADS"] == "1"
    finally:
        # Restore original values to avoid contaminating other tests
        if orig_openblas is None:
            os.environ.pop("OPENBLAS_NUM_THREADS", None)
        else:
            os.environ["OPENBLAS_NUM_THREADS"] = orig_openblas

        if orig_mkl is None:
            os.environ.pop("MKL_NUM_THREADS", None)
        else:
            os.environ["MKL_NUM_THREADS"] = orig_mkl

        if orig_omp is None:
            os.environ.pop("OMP_NUM_THREADS", None)
        else:
            os.environ["OMP_NUM_THREADS"] = orig_omp
