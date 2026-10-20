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
- _worker_initializer limits BLAS threads in worker processes.
- Workers build the genotype matrices; the parent fits null models once (#212).
"""

from __future__ import annotations

import os
from typing import Any, ClassVar
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.engine import AssociationEngine, _worker_initializer

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _parallel_for_small_panels(monkeypatch: pytest.MonkeyPatch) -> None:
    """Exercise the parallel path with the small synthetic panels used here."""
    monkeypatch.setattr("variantcentrifuge.association.engine.MIN_GENES_FOR_PARALLEL", 2)


def test_small_panels_stay_sequential_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    """Panels below MIN_GENES_FOR_PARALLEL never start a worker pool."""
    monkeypatch.setattr("variantcentrifuge.association.engine.MIN_GENES_FOR_PARALLEL", 50)
    engine = AssociationEngine.from_names(["fisher"], AssociationConfig(association_workers=4))
    with patch("variantcentrifuge.association.engine.dispatch_parallel_genes") as dispatch:
        df = engine.run_all(_make_gene_data(10))
    dispatch.assert_not_called()
    assert len(df) == 10


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
    assert "unsafe_test_pvalue" in df.columns


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
    assert "fisher_pvalue" in df.columns
    # p_value should be a valid float (not None)
    assert df["fisher_pvalue"].iloc[0] is not None


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


@pytest.mark.unit
def test_parallel_association_payload_cleanup_on_batch_error() -> None:
    """Matrix payloads in batch are discarded when a builder or worker raises an error."""

    def faulty_builder() -> dict[str, Any]:
        raise RuntimeError("Simulated builder failure in batch")

    def normal_builder() -> dict[str, Any]:
        return {
            "genotype_matrix": np.zeros((4, 2), dtype=float),
            "variant_mafs": np.array([0.1, 0.2]),
        }

    genes = [
        {
            "GENE": "G0",
            "proband_count": 10,
            "control_count": 10,
            "_genotype_matrix_builder": normal_builder,
        },
        {
            "GENE": "G1",
            "proband_count": 10,
            "control_count": 10,
            "_genotype_matrix_builder": faulty_builder,
        },
    ]

    config = AssociationConfig(association_workers=2)
    engine = AssociationEngine.from_names(["fisher"], config)
    with pytest.raises(RuntimeError, match="Simulated builder failure"):
        engine.run_all(genes)

    # All materialized matrix payloads must be cleaned up via finally
    for g in genes:
        assert "genotype_matrix" not in g
        assert "variant_mafs" not in g


@pytest.mark.unit
def test_parallel_vs_sequential_parity_with_real_builder() -> None:
    """Parallel execution with real _GenotypeMatrixBuilder matches sequential parity."""
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    def make_test_gene(name: str, n_samples: int = 30) -> dict[str, Any]:
        vcf_samples = [f"sample_{i}" for i in range(n_samples)]
        gt_cols = [f"GEN_{i}__GT" for i in range(n_samples)]
        gene_df = pd.DataFrame(
            {
                "GENE": [name] * 3,
                "CHROM": ["chr1"] * 3,
                "POS": [100, 200, 300],
                "REF": ["A", "C", "G"],
                "ALT": ["T", "T", "T"],
                **{col: ["0/0", "0/0", "0/0"] for col in gt_cols},
            }
        )
        # Sufficient carriers so MAC >= 5
        gene_df.loc[0, gt_cols[0]] = "0/1"
        gene_df.loc[0, gt_cols[1]] = "0/1"
        gene_df.loc[1, gt_cols[2]] = "0/1"
        gene_df.loc[1, gt_cols[3]] = "0/1"
        gene_df.loc[2, gt_cols[16]] = "0/1"
        gene_df.loc[2, gt_cols[17]] = "0/1"
        phenotype = np.array([1.0] * 15 + [0.0] * 15)
        builder = _GenotypeMatrixBuilder(
            gene_df=gene_df,
            vcf_samples=vcf_samples,
            gt_columns=gt_cols,
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.1,
            phenotype_vector=phenotype,
            covariate_matrix=None,
        )
        return {
            "GENE": name,
            "proband_count": 15,
            "control_count": 15,
            "n_qualifying_variants": 3,
            "phenotype_vector": phenotype,
            "_genotype_matrix_builder": builder,
        }

    genes_seq = [make_test_gene(f"GENE_{i}") for i in range(6)]
    genes_par = [make_test_gene(f"GENE_{i}") for i in range(6)]

    cfg_seq = AssociationConfig(association_workers=1)
    cfg_par = AssociationConfig(association_workers=2)

    eng_seq = AssociationEngine.from_names(["logistic_burden"], cfg_seq)
    eng_par = AssociationEngine.from_names(["logistic_burden"], cfg_par)

    df_seq = eng_seq.run_all(genes_seq).sort_values("gene").reset_index(drop=True)
    df_par = eng_par.run_all(genes_par).sort_values("gene").reset_index(drop=True)

    assert len(df_seq) == 6
    pd.testing.assert_frame_equal(df_seq, df_par, check_exact=False, rtol=1e-8)


def _real_gene(name: str, seed: int, n_samples: int = 60) -> dict[str, Any]:
    """A gene with a real lazy builder, binary phenotype and two covariates."""
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    rng = np.random.default_rng(seed)
    gt_cols = [f"GEN_{i}__GT" for i in range(n_samples)]
    codes = rng.choice(["0/0", "0/1", "1/1"], size=(6, n_samples), p=[0.85, 0.12, 0.03])
    gene_df = pd.DataFrame(
        {"GENE": name, "CHROM": "1", "POS": np.arange(6) * 10 + 100, "REF": "A", "ALT": "T"}
    )
    gene_df = pd.concat([gene_df, pd.DataFrame(codes, columns=gt_cols)], axis=1)
    phenotype = np.array([1.0] * 20 + [0.0] * (n_samples - 20))
    covariates = np.random.default_rng(0).normal(size=(n_samples, 2))
    builder = _GenotypeMatrixBuilder(
        gene_df=gene_df,
        vcf_samples=[f"s{i}" for i in range(n_samples)],
        gt_columns=gt_cols,
        is_binary=True,
        missing_site_threshold=0.1,
        missing_sample_threshold=0.8,
        phenotype_vector=phenotype,
        covariate_matrix=covariates,
    )
    return {
        "GENE": name,
        "proband_count": 20,
        "control_count": n_samples - 20,
        "n_qualifying_variants": 6,
        "_genotype_matrix_builder": builder,
    }


@pytest.mark.unit
def test_worker_count_invariance_exact() -> None:
    """1 and 3 workers give bit-identical tables and the same n_variants_qc (#212)."""
    cfg = {"skat_method": "SKATO", "trait_type": "binary"}
    tables, qc = [], []
    for workers in (1, 3):
        genes = [_real_gene(f"G{i:02d}", i) for i in range(14)]
        engine = AssociationEngine.from_names(
            ["fisher", "logistic_burden", "skat_python"],
            AssociationConfig(association_workers=workers, **cfg),
        )
        tables.append(engine.run_all(genes))
        qc.append([g.get("n_variants_qc") for g in genes])
        assert not any("_genotype_matrix_builder" in g or "genotype_matrix" in g for g in genes)
    assert len(tables[0]) == 14
    pd.testing.assert_frame_equal(tables[0], tables[1], check_exact=True)
    assert qc[0] == qc[1] and None not in qc[0]


@pytest.mark.unit
def test_null_model_fitted_once_across_workers() -> None:
    """The SKAT null model is fitted in the parent only, never again in a worker (#212)."""
    import multiprocessing

    from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

    calls = multiprocessing.get_context("fork").Value("i", 0)  # shared with forked workers
    real_fit = PythonSKATBackend.fit_null_model

    def counting_fit(self, *args, **kwargs):
        with calls.get_lock():
            calls.value += 1
        return real_fit(self, *args, **kwargs)

    genes = [_real_gene(f"G{i:02d}", i) for i in range(12)]
    genes.insert(0, _real_gene("A_EMPTY", 0))  # first gene sorts first and fits no null model
    genes[0]["_genotype_matrix_builder"].gene_df = genes[1][
        "_genotype_matrix_builder"
    ].gene_df.iloc[:0]
    engine = AssociationEngine.from_names(
        ["skat_python"], AssociationConfig(association_workers=3, skat_method="SKATO")
    )
    with patch.object(PythonSKATBackend, "fit_null_model", counting_fit):
        df = engine.run_all(genes)
    assert len(df) == 13
    assert calls.value == 1


_BLAS_GETTER: Any = None  # inherited by the forked worker below


def _blas_threads_in_worker() -> int:
    return int(_BLAS_GETTER())


def _openblas_controls() -> tuple[Any, Any]:
    """(getter, setter) of a loaded OpenBLAS, found independently of the code under test."""
    import ctypes

    with open("/proc/self/maps") as fh:
        paths = sorted({ln.split()[-1] for ln in fh if "openblas" in ln and ".so" in ln})
    pairs = [
        (getattr(lib, g), getattr(lib, g.replace("get", "set")))
        for lib in map(ctypes.CDLL, paths)
        for g in ("openblas_get_num_threads", "scipy_openblas_get_num_threads64_")
        if hasattr(lib, g)
    ]
    if not pairs:
        pytest.skip("no OpenBLAS with a known thread getter loaded")
    return pairs[0]


@pytest.mark.unit
def test_worker_initializer_limits_loaded_blas() -> None:
    """In a forked worker the already-loaded BLAS runs one thread, whatever the parent uses."""
    global _BLAS_GETTER
    import concurrent.futures
    import multiprocessing

    _BLAS_GETTER, set_threads = _openblas_controls()
    before = _BLAS_GETTER()
    set_threads(2)
    try:
        with concurrent.futures.ProcessPoolExecutor(
            1, mp_context=multiprocessing.get_context("fork"), initializer=_worker_initializer
        ) as ex:
            assert ex.submit(_blas_threads_in_worker).result() == 1
        assert _BLAS_GETTER() == 2  # the parent keeps its own setting
    finally:
        set_threads(before)


class _BlasProbeTest(AssociationTest):
    """Records the BLAS thread count seen while the engine runs a gene."""

    parallel_safe = True
    seen: ClassVar[list[int]] = []

    @property
    def name(self) -> str:
        return "probe"

    def run(self, gene: str, contingency_data: dict[str, Any], config: AssociationConfig):
        _BlasProbeTest.seen.append(int(_BLAS_GETTER()))
        return TestResult(gene, self.name, 0.5, None, None, None, None, None, 1, 1, 1)

    def effect_column_names(self) -> dict[str, str | None]:
        return {"effect": None, "se": None, "ci_lower": None, "ci_upper": None}


@pytest.mark.unit
def test_sequential_engine_runs_genes_with_one_blas_thread() -> None:
    """workers=1 also pins BLAS to one thread during the gene loop and restores it after."""
    global _BLAS_GETTER
    _BLAS_GETTER, set_threads = _openblas_controls()
    before = _BLAS_GETTER()
    set_threads(2)
    try:
        _BlasProbeTest.seen = []
        AssociationEngine([_BlasProbeTest()], AssociationConfig(association_workers=1)).run_all(
            _make_gene_data(3)
        )
        assert _BlasProbeTest.seen == [1, 1, 1]
        assert _BLAS_GETTER() == 2
    finally:
        set_threads(before)


class _PidRecordingBuilder:
    """Wraps a builder and counts calls made in the parent process."""

    def __init__(self, inner: Any, parent_calls: Any) -> None:
        self.inner, self.parent_calls, self.parent = inner, parent_calls, os.getpid()

    def __call__(self) -> dict[str, Any]:
        if os.getpid() == self.parent:
            with self.parent_calls.get_lock():
                self.parent_calls.value += 1
        return self.inner()


@pytest.mark.unit
def test_matrices_are_built_in_workers() -> None:
    """Only the parent's null-model prefix builds a matrix there; workers build the rest."""
    import multiprocessing

    parent_calls = multiprocessing.get_context("fork").Value("i", 0)
    genes = [_real_gene(f"G{i:02d}", i) for i in range(12)]
    for g in genes:
        g["_genotype_matrix_builder"] = _PidRecordingBuilder(
            g["_genotype_matrix_builder"], parent_calls
        )
    engine = AssociationEngine.from_names(
        ["logistic_burden", "skat_python"], AssociationConfig(association_workers=3)
    )
    assert len(engine.run_all(genes)) == 12
    assert parent_calls.value == 1


@pytest.mark.unit
def test_without_fork_genes_run_sequentially(monkeypatch: pytest.MonkeyPatch) -> None:
    """Off Linux (no fork) the parallel request falls back to the sequential loop."""
    monkeypatch.setattr(
        "variantcentrifuge.association.engine_workers.fork_available", lambda: False
    )
    tables = []
    for workers in (1, 3):
        genes = [_real_gene(f"G{i:02d}", i) for i in range(8)]
        engine = AssociationEngine.from_names(
            ["logistic_burden"], AssociationConfig(association_workers=workers)
        )
        with patch("variantcentrifuge.association.engine_workers.fork_map") as fork_map:
            tables.append(engine.run_all(genes))
        fork_map.assert_not_called()
    pd.testing.assert_frame_equal(tables[0], tables[1], check_exact=True)


@pytest.mark.unit
def test_fork_task_runs_gene_range_in_process(monkeypatch: pytest.MonkeyPatch) -> None:
    """The worker entry point (normally only run in a child) returns results and QC counts."""
    from variantcentrifuge.association import engine_workers as ew

    genes = [_real_gene(f"G{i:02d}", i) for i in range(3)]
    tests = {
        "logistic_burden": AssociationEngine.from_names(
            ["logistic_burden"], AssociationConfig()
        )._tests["logistic_burden"]
    }
    monkeypatch.setattr(
        ew, "_FORK_STATE", (ew._run_gene_range, (genes, tests, AssociationConfig(), None))
    )
    results, qc, hooks = ew._fork_task((1, 3))
    assert hooks == {}
    assert list(results["logistic_burden"]) == ["G01", "G02"]
    assert qc == {"G01": 6, "G02": 6}


@pytest.mark.unit
def test_gc_unfrozen_after_parallel_run_and_on_error() -> None:
    """fork_map freezes the GC for the fork and always unfreezes it (#212)."""
    import gc

    frozen = gc.get_freeze_count()  # objects frozen elsewhere stay frozen (may die meanwhile)
    genes = [_real_gene(f"G{i:02d}", i) for i in range(6)]
    engine = AssociationEngine.from_names(
        ["logistic_burden"], AssociationConfig(association_workers=2)
    )
    engine.run_all(genes)
    assert gc.get_freeze_count() <= frozen

    def boom() -> dict[str, Any]:
        raise RuntimeError("worker failure")

    genes = [_real_gene(f"G{i:02d}", i) for i in range(6)]
    genes[-1]["_genotype_matrix_builder"] = boom
    with pytest.raises(RuntimeError, match="worker failure"):
        engine.run_all(genes)
    assert gc.get_freeze_count() <= frozen


def _range_len(state: Any, lo: int, hi: int) -> int:
    return hi - lo


@pytest.mark.unit
def test_fork_map_keeps_a_callers_gc_freeze() -> None:
    """Codex FX-4: fork_map must not unfreeze objects the caller froze."""
    import gc

    from variantcentrifuge.association.engine_workers import fork_map

    gc.freeze()
    try:
        before = gc.get_freeze_count()
        assert before > 0
        assert sum(r for _, _, r in fork_map(_range_len, None, 10, 2)) == 10
        assert before - 100 < gc.get_freeze_count() <= before  # frozen objects may be freed
    finally:
        gc.unfreeze()
