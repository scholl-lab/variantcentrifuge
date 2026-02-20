"""
Unit tests for RSKATBackend with mocked rpy2.

All tests mock rpy2 so that no actual R installation is required.
The tests verify correct R function dispatch, NA handling, memory
management, and thread safety.

Covers requirements: SKAT-01 through SKAT-09

Mocking strategy
----------------
Since r_backend.py uses deferred imports (``import rpy2.robjects as ro``
inside method bodies), we patch via ``sys.modules`` injection AND link
parent module attributes to their sub-modules. The critical pattern is:

    mock_rpy2.robjects = mock_ro
    mock_ro.packages = mock_rpacks

This ensures that ``import rpy2.robjects.packages as rpacks`` inside the
method body resolves to our mock (Python walks mock_rpy2.robjects.packages).
"""

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Core helper: build a properly-linked rpy2 mock hierarchy
# ---------------------------------------------------------------------------


def _build_rpy2_mocks(
    *,
    isinstalled_return: bool = True,
    rpy2_version: str = "3.5.14",
    na_real_sentinel: object | None = None,
) -> tuple[MagicMock, MagicMock, MagicMock, object]:
    """
    Build a properly-linked rpy2 mock hierarchy.

    Returns
    -------
    (mock_rpy2, mock_ro, mock_rpacks, na_sentinel)
        Linked mock objects where mock_rpy2.robjects = mock_ro and
        mock_ro.packages = mock_rpacks.
    """
    if na_real_sentinel is None:
        na_real_sentinel = object()

    # rpy2.robjects.packages
    mock_rpacks = MagicMock(name="rpy2.robjects.packages")
    mock_rpacks.isinstalled = MagicMock(return_value=isinstalled_return)
    mock_skat_pkg = MagicMock(name="SKAT_pkg")
    mock_base_pkg = MagicMock(name="base_pkg")
    mock_rpacks.importr = MagicMock(
        side_effect=lambda name: mock_skat_pkg if name == "SKAT" else mock_base_pkg
    )

    # rpy2.robjects
    mock_ro = MagicMock(name="rpy2.robjects")
    mock_ro.packages = mock_rpacks  # parent attribute link
    mock_ro.globalenv = {}
    mock_ro.FloatVector = MagicMock(side_effect=lambda x: list(x))
    mock_ro.Formula = MagicMock(side_effect=lambda x: x)

    def default_r(code: str) -> Any:
        if "R.Version()" in code:
            ver_list = MagicMock()
            ver_list.names = ["major", "minor"]
            ver_list.__getitem__ = lambda self, idx: [["4"], ["3.0"]][idx]
            return ver_list
        if "packageVersion" in code:
            ver = MagicMock()
            ver.__getitem__ = lambda self, idx: "2.2.5"
            ver.__len__ = lambda self: 1
            return ver
        result = MagicMock()
        result.__len__ = lambda self: 0
        return result

    mock_ro.r = MagicMock(side_effect=default_r)
    mock_ro.r.__getitem__ = MagicMock(return_value=MagicMock(return_value=MagicMock()))

    # rpy2.rinterface
    mock_rinterface = MagicMock(name="rpy2.rinterface")
    mock_rinterface.NA_Real = na_real_sentinel

    # rpy2 root — link children
    mock_rpy2 = MagicMock(name="rpy2")
    mock_rpy2.__version__ = rpy2_version
    mock_rpy2.robjects = mock_ro
    mock_rpy2.rinterface = mock_rinterface

    return mock_rpy2, mock_ro, mock_rpacks, na_real_sentinel


def _patch_rpy2(
    *,
    isinstalled_return: bool = True,
    rpy2_version: str = "3.5.14",
    na_real_sentinel: object | None = None,
) -> tuple[dict[str, Any], MagicMock, MagicMock, MagicMock, object]:
    """
    Build sys.modules patch dict and return (mods, mock_rpy2, mock_ro, mock_rpacks, na_sentinel).

    Use with ``patch.dict("sys.modules", mods)``.
    """
    mock_rpy2, mock_ro, mock_rpacks, na_sentinel = _build_rpy2_mocks(
        isinstalled_return=isinstalled_return,
        rpy2_version=rpy2_version,
        na_real_sentinel=na_real_sentinel,
    )

    mock_embedded = MagicMock(name="rpy2.rinterface_lib.embedded")

    class RRuntimeError(Exception):
        pass

    mock_embedded.RRuntimeError = RRuntimeError

    mock_rinterface_lib = MagicMock(name="rpy2.rinterface_lib")
    mock_rinterface_lib.embedded = mock_embedded

    mods = {
        "rpy2": mock_rpy2,
        "rpy2.robjects": mock_ro,
        "rpy2.robjects.packages": mock_rpacks,
        "rpy2.rinterface": mock_rpy2.rinterface,
        "rpy2.rinterface_lib": mock_rinterface_lib,
        "rpy2.rinterface_lib.embedded": mock_embedded,
    }
    return mods, mock_rpy2, mock_ro, mock_rpacks, na_sentinel


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def _make_backend_with_mocked_env() -> Any:
    """Create an RSKATBackend instance that skips detect_environment()."""
    from variantcentrifuge.association.backends.r_backend import RSKATBackend

    backend = RSKATBackend()
    backend._skat_pkg = MagicMock(name="skat_pkg")
    backend._base_pkg = MagicMock(name="base_pkg")
    backend._r_gc_func = MagicMock(name="r_gc_func")
    backend._rpy2_version = "3.5.14"
    backend._r_version = "4.3"
    backend._skat_version = "2.2.5"
    backend._r_home = "/usr/lib/R"
    return backend


def _make_null_model(trait_type: str = "binary") -> Any:
    """Create a NullModelResult with a mock R model object."""
    from variantcentrifuge.association.backends.base import NullModelResult

    return NullModelResult(
        model=MagicMock(name="r_null_model"),
        trait_type=trait_type,
        n_samples=100,
        adjustment=True,
    )


def _make_r_outer_result(
    p_value_scalar: Any,
    rho_scalar: Any = None,
    n_marker_test_scalar: Any = None,
    warnings_list: list[str] | None = None,
    na_real_sentinel: object | None = None,
) -> Any:
    """
    Build a mock for the R outer list returned by ro.r(r_code).

    Structure: outer.rx2("result") -> inner result
               outer.rx2("warnings") -> char vector
    """
    if warnings_list is None:
        warnings_list = []

    # Inner result mock
    inner_result = MagicMock(name="inner_result")

    def inner_rx2_dispatch(key: str) -> Any:
        if key == "p.value":
            v = MagicMock()
            v.__getitem__ = lambda self, idx: p_value_scalar
            return v
        if key == "param":
            p = MagicMock(name="param")

            def param_rx2(k: str) -> Any:
                v = MagicMock()
                if k == "rho_est":
                    v.__getitem__ = lambda self, idx: rho_scalar
                elif k == "n.marker.test":
                    v.__getitem__ = lambda self, idx: n_marker_test_scalar
                else:
                    v.__getitem__ = lambda self, idx: na_real_sentinel
                return v

            p.rx2 = param_rx2
            return p
        return MagicMock()

    inner_result.rx2 = inner_rx2_dispatch

    # Warnings vector
    warns_mock = MagicMock()
    warns_mock.__len__ = lambda self: len(warnings_list)
    warns_mock.__iter__ = lambda self: iter(warnings_list)

    # Outer mock
    outer = MagicMock(name="outer_result")

    def outer_rx2_dispatch(key: str) -> Any:
        if key == "result":
            return inner_result
        if key == "warnings":
            return warns_mock
        return MagicMock()

    outer.rx2 = outer_rx2_dispatch
    return outer


def _build_ro_r_with_outer(
    outer_result: Any,
    na_real_sentinel: object,
    extra_rm_tracker: list[str] | None = None,
) -> MagicMock:
    """
    Build a mock for ``mock_ro.r`` that:
    - Returns outer_result for SKAT test code (contains 'local({')
    - Tracks rm() calls in extra_rm_tracker if provided
    - Returns a no-op mock for other calls
    """

    def r_side_effect(code: str) -> Any:
        if extra_rm_tracker is not None and "rm(" in code:
            extra_rm_tracker.append(code)
        if "local({" in code:
            return outer_result
        result = MagicMock()
        result.__len__ = lambda self: 0
        return result

    mock_r = MagicMock(side_effect=r_side_effect)
    mock_r.__getitem__ = MagicMock(return_value=MagicMock(return_value=MagicMock()))
    return mock_r


# ---------------------------------------------------------------------------
# detect_environment() tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestDetectEnvironment:
    """Tests for RSKATBackend.detect_environment() error paths and success."""

    def test_detect_environment_no_rpy2(self):
        """detect_environment raises ImportError with R_HOME hint when rpy2 unavailable."""
        from variantcentrifuge.association.backends.r_backend import RSKATBackend

        backend = RSKATBackend()

        # Block all rpy2 imports
        original_import = __import__

        def blocking_import(name: str, *args: Any, **kwargs: Any) -> Any:
            if name.startswith("rpy2"):
                raise ImportError(f"Mocked: no module '{name}'")
            return original_import(name, *args, **kwargs)

        with (
            patch("builtins.__import__", side_effect=blocking_import),
            pytest.raises(ImportError) as exc_info,
        ):
            backend.detect_environment()

        assert "R_HOME" in str(exc_info.value)

    def test_detect_environment_no_skat_package(self):
        """detect_environment raises ImportError with install.packages hint when SKAT missing."""
        from variantcentrifuge.association.backends.r_backend import RSKATBackend

        backend = RSKATBackend()
        mods, _, _, _, _ = _patch_rpy2(isinstalled_return=False)

        with patch.dict("sys.modules", mods), pytest.raises(ImportError) as exc_info:
            backend.detect_environment()

        assert "install.packages('SKAT')" in str(exc_info.value)

    def test_detect_environment_success(self):
        """detect_environment succeeds and caches _skat_pkg and _base_pkg."""
        from variantcentrifuge.association.backends.r_backend import RSKATBackend

        backend = RSKATBackend()
        mods, _, _, _, _ = _patch_rpy2(isinstalled_return=True)

        with patch.dict("sys.modules", mods):
            # Should not raise
            backend.detect_environment()

        assert backend._skat_pkg is not None
        assert backend._base_pkg is not None


# ---------------------------------------------------------------------------
# log_environment() tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLogEnvironment:
    """Tests for RSKATBackend.log_environment()."""

    def test_log_environment_versions(self, caplog):
        """log_environment emits INFO with all versions and R_HOME."""
        import logging

        backend = _make_backend_with_mocked_env()
        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            backend.log_environment()

        assert any("3.5.14" in msg for msg in caplog.messages)
        assert any("4.3" in msg for msg in caplog.messages)
        assert any("2.2.5" in msg for msg in caplog.messages)
        assert any("/usr/lib/R" in msg for msg in caplog.messages)

    def test_log_environment_warns_old_r_version(self, caplog):
        """log_environment emits WARNING when R version < 4.0."""
        import logging

        backend = _make_backend_with_mocked_env()
        backend._r_version = "3.6"

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            backend.log_environment()

        warning_msgs = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warning_msgs) >= 1
        assert any("3.6" in m.message for m in warning_msgs)

    def test_log_environment_no_warnings_for_good_versions(self, caplog):
        """log_environment emits no WARNING when all versions meet requirements."""
        import logging

        backend = _make_backend_with_mocked_env()
        backend._r_version = "4.3"
        backend._skat_version = "2.2.5"
        backend._rpy2_version = "3.5.14"

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            backend.log_environment()

        warning_msgs = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warning_msgs) == 0


# ---------------------------------------------------------------------------
# fit_null_model() tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFitNullModel:
    """Tests for RSKATBackend.fit_null_model()."""

    def test_fit_null_model_binary_with_covariates(self):
        """fit_null_model calls SKAT_Null_Model with out_type='D' and Adjustment=True."""
        backend = _make_backend_with_mocked_env()
        phenotype = np.array([1.0, 0.0, 1.0, 0.0, 1.0] * 20)
        covariates = np.random.randn(100, 2)

        mods, _, mock_ro, _, _ = _patch_rpy2()
        mock_ro.globalenv = {}

        with patch.dict("sys.modules", mods):
            result = backend.fit_null_model(phenotype, covariates, "binary")

        backend._skat_pkg.SKAT_Null_Model.assert_called_once()
        call_kwargs = backend._skat_pkg.SKAT_Null_Model.call_args
        assert call_kwargs.kwargs.get("out_type") == "D"
        assert call_kwargs.kwargs.get("Adjustment") is True
        assert result.trait_type == "binary"
        assert result.n_samples == 100

    def test_fit_null_model_continuous_no_covariates(self):
        """fit_null_model calls SKAT_Null_Model with out_type='C' and intercept formula."""
        backend = _make_backend_with_mocked_env()
        phenotype = np.random.randn(100)

        mods, _, mock_ro, _, _ = _patch_rpy2()
        mock_ro.globalenv = {}
        captured_formula: dict[str, Any] = {}

        mock_ro.Formula = MagicMock(
            side_effect=lambda fs: captured_formula.update({"formula": fs}) or fs
        )

        with patch.dict("sys.modules", mods):
            result = backend.fit_null_model(phenotype, None, "quantitative")

        assert captured_formula.get("formula") == "._vc_y ~ 1"
        call_kwargs = backend._skat_pkg.SKAT_Null_Model.call_args
        assert call_kwargs.kwargs.get("out_type") == "C"
        assert result.trait_type == "quantitative"

    def test_fit_null_model_asserts_main_thread(self):
        """fit_null_model raises RuntimeError when called from a non-main thread."""
        backend = _make_backend_with_mocked_env()
        phenotype = np.array([1.0, 0.0] * 50)
        error_holder: list[Exception] = []

        def run_in_thread() -> None:
            try:
                backend.fit_null_model(phenotype, None, "binary")
            except RuntimeError as e:
                error_holder.append(e)

        t = threading.Thread(target=run_in_thread)
        t.start()
        t.join()

        assert len(error_holder) == 1
        assert "main thread" in str(error_holder[0]).lower()

    def test_fit_null_model_cleans_r_globals(self):
        """fit_null_model calls ro.r('rm(...)') to clean up ._vc_ globals."""
        backend = _make_backend_with_mocked_env()
        phenotype = np.array([1.0, 0.0] * 50)

        mods, _, mock_ro, _, _ = _patch_rpy2()
        mock_ro.globalenv = {}
        rm_calls: list[str] = []

        original_side = mock_ro.r.side_effect

        def capturing_r(code: str) -> Any:
            if "rm(" in code:
                rm_calls.append(code)
            return original_side(code)

        mock_ro.r = MagicMock(side_effect=capturing_r)
        mock_ro.r.__getitem__ = MagicMock(return_value=MagicMock(return_value=MagicMock()))

        with patch.dict("sys.modules", mods):
            backend.fit_null_model(phenotype, None, "binary")

        assert any("rm(" in c and "._vc_" in c for c in rm_calls), (
            f"Expected rm(list=ls(pattern='._vc_')) to be called, got: {rm_calls}"
        )

    def test_fit_null_model_resets_gene_counter(self):
        """fit_null_model resets _genes_processed to 0 after fitting."""
        backend = _make_backend_with_mocked_env()
        backend._genes_processed = 999
        phenotype = np.array([1.0, 0.0] * 50)

        mods, _, mock_ro, _, _ = _patch_rpy2()
        mock_ro.globalenv = {}

        with patch.dict("sys.modules", mods):
            backend.fit_null_model(phenotype, None, "binary")

        assert backend._genes_processed == 0


# ---------------------------------------------------------------------------
# test_gene() tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestTestGene:
    """Tests for RSKATBackend.test_gene()."""

    def _make_mods_with_result(
        self,
        outer_result: Any,
        na_sentinel: object,
        rm_tracker: list[str] | None = None,
        code_tracker: list[str] | None = None,
        geno_3x2: bool = False,
    ) -> dict[str, Any]:
        """Build patched sys.modules where ro.r returns outer_result for SKAT calls."""
        mods, _, mock_ro, _, _ = _patch_rpy2(na_real_sentinel=na_sentinel)
        mock_ro.globalenv = {}

        def r_side_effect(code: str) -> Any:
            if code_tracker is not None:
                code_tracker.append(code)
            if rm_tracker is not None and "rm(" in code:
                rm_tracker.append(code)
            if "local({" in code:
                return outer_result
            result = MagicMock()
            result.__len__ = lambda self: 0
            return result

        mock_ro.r = MagicMock(side_effect=r_side_effect)
        mock_ro.r.__getitem__ = MagicMock(return_value=MagicMock(return_value=MagicMock()))

        if geno_3x2:
            captured_fv: list[list] = []
            captured_mx: dict[str, Any] = {}

            def capture_fv(data: Any) -> Any:
                captured_fv.append(list(data))
                return data

            matrix_func = MagicMock()

            def capture_matrix(data: Any, nrow: int, ncol: int) -> Any:
                captured_mx.update({"data": data, "nrow": nrow, "ncol": ncol})
                return MagicMock()

            matrix_func.side_effect = capture_matrix
            mock_ro.FloatVector = MagicMock(side_effect=capture_fv)
            mock_ro.r.__getitem__ = MagicMock(return_value=matrix_func)
            # store on mods for access in test
            mods["_captured_fv"] = captured_fv
            mods["_captured_mx"] = captured_mx

        return mods

    def test_test_gene_binary_calls_skatbinary(self):
        """test_gene dispatches to SKATBinary for binary trait."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))
        na_sentinel = object()

        outer_result = _make_r_outer_result(
            p_value_scalar=0.05, n_marker_test_scalar=3, na_real_sentinel=na_sentinel
        )
        code_tracker: list[str] = []
        mods = self._make_mods_with_result(outer_result, na_sentinel, code_tracker=code_tracker)

        with patch.dict("sys.modules", mods):
            result = backend.test_gene(
                gene="BRCA1",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKAT",
                weights_beta=(1.0, 25.0),
            )

        assert any("SKATBinary" in code for code in code_tracker), (
            f"Expected 'SKATBinary' in R code, got: {code_tracker}"
        )
        assert result["p_value"] == 0.05

    def test_test_gene_continuous_calls_skat(self):
        """test_gene dispatches to SKAT (not SKATBinary) for continuous trait."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="quantitative")
        geno = np.zeros((10, 3))
        na_sentinel = object()

        outer_result = _make_r_outer_result(
            p_value_scalar=0.12, n_marker_test_scalar=3, na_real_sentinel=na_sentinel
        )
        code_tracker: list[str] = []
        mods = self._make_mods_with_result(outer_result, na_sentinel, code_tracker=code_tracker)

        with patch.dict("sys.modules", mods):
            backend.test_gene(
                gene="MYH7",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKAT",
                weights_beta=(1.0, 25.0),
            )

        skat_binary_calls = [c for c in code_tracker if "SKATBinary" in c]
        skat_calls = [c for c in code_tracker if "SKAT(" in c and "SKATBinary" not in c]
        assert len(skat_binary_calls) == 0, "SKATBinary should NOT be called for continuous"
        assert len(skat_calls) >= 1, f"SKAT should be called for continuous, got: {code_tracker}"

    def test_test_gene_skat_o_extracts_rho(self):
        """SKAT-O method extracts optimal rho from result.rx2('param').rx2('rho_est')."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))
        na_sentinel = object()

        outer_result = _make_r_outer_result(
            p_value_scalar=0.03,
            rho_scalar=0.5,
            n_marker_test_scalar=3,
            na_real_sentinel=na_sentinel,
        )
        mods = self._make_mods_with_result(outer_result, na_sentinel)

        with patch.dict("sys.modules", mods):
            result = backend.test_gene(
                gene="TP53",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKATO",
                weights_beta=(1.0, 25.0),
            )

        assert result["rho"] == 0.5

    def test_test_gene_na_real_returns_none(self):
        """NA_Real from R maps to p_value=None in the result dict."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))

        # Use a unique object as NA_Real sentinel
        na_sentinel = object()
        # p.value returns the sentinel (== NA_Real)
        outer_result = _make_r_outer_result(
            p_value_scalar=na_sentinel,
            n_marker_test_scalar=3,
            na_real_sentinel=na_sentinel,
        )
        mods = self._make_mods_with_result(outer_result, na_sentinel)

        with patch.dict("sys.modules", mods):
            result = backend.test_gene(
                gene="GENE1",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKAT",
                weights_beta=(1.0, 25.0),
            )

        assert result["p_value"] is None

    def test_test_gene_captures_r_warnings(self):
        """Warnings returned from R code are included in result['warnings']."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))
        na_sentinel = object()

        outer_result = _make_r_outer_result(
            p_value_scalar=0.04,
            n_marker_test_scalar=3,
            warnings_list=["moment adjustment applied"],
            na_real_sentinel=na_sentinel,
        )
        mods = self._make_mods_with_result(outer_result, na_sentinel)

        with patch.dict("sys.modules", mods):
            result = backend.test_gene(
                gene="GENE2",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKAT",
                weights_beta=(1.0, 25.0),
            )

        assert "moment adjustment applied" in result["warnings"]

    def test_test_gene_matrix_column_major(self):
        """Genotype matrix is passed to R in column-major (Fortran) order."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.array([[1, 2], [3, 4], [5, 6]])  # 3 samples x 2 variants
        na_sentinel = object()

        outer_result = _make_r_outer_result(
            p_value_scalar=0.1, n_marker_test_scalar=2, na_real_sentinel=na_sentinel
        )
        mods = self._make_mods_with_result(outer_result, na_sentinel, geno_3x2=True)

        with patch.dict("sys.modules", mods):
            backend.test_gene(
                gene="GTEST",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKAT",
                weights_beta=(1.0, 25.0),
            )

        captured_fv: list[list] = mods["_captured_fv"]
        captured_mx: dict[str, Any] = mods["_captured_mx"]

        # Column-major of [[1,2],[3,4],[5,6]] = [1,3,5,2,4,6]
        expected_col_major = list(geno.ravel(order="F").tolist())
        assert any(fv == expected_col_major for fv in captured_fv), (
            f"Expected column-major order {expected_col_major}, got: {captured_fv}"
        )
        assert captured_mx.get("nrow") == 3
        assert captured_mx.get("ncol") == 2

    def test_test_gene_r_exception_propagates(self):
        """R runtime errors propagate directly (not caught) — abort on infrastructure failure."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))
        na_sentinel = object()

        mods, _, mock_ro, _, _ = _patch_rpy2(na_real_sentinel=na_sentinel)
        mock_ro.globalenv = {}

        class RRuntimeError(Exception):
            pass

        mock_ro.r = MagicMock(side_effect=RRuntimeError("R computation failed"))
        mock_ro.r.__getitem__ = MagicMock(return_value=MagicMock(return_value=MagicMock()))

        with patch.dict("sys.modules", mods), pytest.raises(RRuntimeError):
            backend.test_gene(
                gene="GENE3",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKAT",
                weights_beta=(1.0, 25.0),
            )

    def test_test_gene_cleans_r_globals_per_gene(self):
        """test_gene calls ro.r('rm(...)') to clean up ._vc_ globals after each call."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))
        na_sentinel = object()

        outer_result = _make_r_outer_result(
            p_value_scalar=0.05, n_marker_test_scalar=3, na_real_sentinel=na_sentinel
        )
        rm_tracker: list[str] = []
        mods = self._make_mods_with_result(outer_result, na_sentinel, rm_tracker=rm_tracker)

        with patch.dict("sys.modules", mods):
            backend.test_gene(
                gene="GENEX",
                genotype_matrix=geno,
                null_model=null_model,
                method="SKAT",
                weights_beta=(1.0, 25.0),
            )

        assert any("rm(" in c and "._vc_" in c for c in rm_tracker), (
            f"Expected R globals cleanup, got rm_tracker: {rm_tracker}"
        )

    def test_test_gene_increments_genes_processed(self):
        """test_gene increments _genes_processed on the backend after each call."""
        backend = _make_backend_with_mocked_env()
        backend._genes_processed = 0
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))
        na_sentinel = object()

        outer_result = _make_r_outer_result(
            p_value_scalar=0.05, n_marker_test_scalar=3, na_real_sentinel=na_sentinel
        )
        mods = self._make_mods_with_result(outer_result, na_sentinel)

        with patch.dict("sys.modules", mods):
            backend.test_gene("GENE_A", geno, null_model, "SKAT", (1.0, 25.0))
            backend.test_gene("GENE_B", geno, null_model, "SKAT", (1.0, 25.0))

        assert backend._genes_processed == 2


# ---------------------------------------------------------------------------
# Memory management tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestMemoryManagement:
    """Tests for RSKATBackend GC interval, heap check, and cleanup."""

    def test_gc_called_at_interval(self):
        """GC is triggered every GC_INTERVAL (100) genes — tested via spy on _run_r_gc."""
        backend = _make_backend_with_mocked_env()

        gc_call_count = [0]

        def counting_gc() -> None:
            gc_call_count[0] += 1

        backend._run_r_gc = counting_gc

        from variantcentrifuge.association.backends.r_backend import RSKATBackend

        # Simulate the GC check logic from RSKATTest.run()
        for i in range(1, 251):  # 250 gene calls
            backend._genes_processed = i
            if backend._genes_processed % RSKATBackend.GC_INTERVAL == 0:
                backend._run_r_gc()

        # At genes 100 and 200 -> exactly 2 calls
        assert gc_call_count[0] == 2

    def test_check_r_heap_returns_mb(self):
        """_check_r_heap() returns heap usage in MB from R gc() Vcells matrix."""
        backend = _make_backend_with_mocked_env()

        # vcells_used = 1048576 => heap_mb = 1048576 * 8 / 1024^2 = 8.0 MB
        vcells_used = 1048576
        expected_mb = vcells_used * 8.0 / (1024.0 * 1024.0)

        gc_matrix = MagicMock()
        rx_result = MagicMock()
        rx_result.__getitem__ = lambda self, idx: float(vcells_used)
        gc_matrix.rx = MagicMock(return_value=rx_result)

        mods, _, mock_ro, _, _ = _patch_rpy2()
        mock_ro.r = MagicMock(return_value=gc_matrix)

        with patch.dict("sys.modules", mods):
            result = backend._check_r_heap()

        assert result == pytest.approx(expected_mb, rel=1e-6)

    def test_check_r_heap_returns_none_on_error(self):
        """_check_r_heap() returns None if R gc() call fails."""
        backend = _make_backend_with_mocked_env()

        mods, _, mock_ro, _, _ = _patch_rpy2()
        mock_ro.r = MagicMock(side_effect=RuntimeError("R gc() failed"))

        with patch.dict("sys.modules", mods):
            result = backend._check_r_heap()

        assert result is None

    def test_cleanup_calls_gc(self):
        """cleanup() calls _run_r_gc."""
        backend = _make_backend_with_mocked_env()

        gc_called = [False]

        def spy_gc() -> None:
            gc_called[0] = True

        backend._run_r_gc = spy_gc
        backend.cleanup()

        assert gc_called[0] is True

    def test_run_r_gc_uses_cached_func(self):
        """_run_r_gc calls the cached _r_gc_func when available."""
        backend = _make_backend_with_mocked_env()

        gc_called = [False]
        backend._r_gc_func = MagicMock(side_effect=lambda: gc_called.__setitem__(0, True))

        backend._run_r_gc()

        assert gc_called[0] is True

    def test_run_r_gc_fallback_when_no_cache(self):
        """_run_r_gc falls back to ro.r['gc']() when _r_gc_func is None."""
        backend = _make_backend_with_mocked_env()
        backend._r_gc_func = None

        gc_called = [False]
        gc_func = MagicMock(side_effect=lambda: gc_called.__setitem__(0, True))

        mods, _, mock_ro, _, _ = _patch_rpy2()
        mock_ro.r.__getitem__ = MagicMock(return_value=gc_func)

        with patch.dict("sys.modules", mods):
            backend._run_r_gc()

        # No exception means fallback ran successfully


# ---------------------------------------------------------------------------
# Thread safety tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestThreadSafety:
    """Tests for RSKATBackend thread safety guard."""

    def test_assert_main_thread_from_worker(self):
        """_assert_main_thread raises RuntimeError from a ThreadPoolExecutor worker."""
        from variantcentrifuge.association.backends.r_backend import RSKATBackend

        backend = RSKATBackend()
        error_holder: list[Exception] = []

        def worker_task() -> None:
            try:
                backend._assert_main_thread()
            except RuntimeError as e:
                error_holder.append(e)

        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(worker_task)
            future.result()

        assert len(error_holder) == 1
        assert "main thread" in str(error_holder[0]).lower()

    def test_assert_main_thread_from_main_thread(self):
        """_assert_main_thread does not raise when called from the main thread."""
        from variantcentrifuge.association.backends.r_backend import RSKATBackend

        backend = RSKATBackend()
        # Should not raise
        backend._assert_main_thread()

    def test_fit_null_model_thread_safety(self):
        """fit_null_model raises RuntimeError from non-main thread."""
        backend = _make_backend_with_mocked_env()
        phenotype = np.array([1.0, 0.0] * 50)
        errors: list[str] = []

        def run_in_thread() -> None:
            try:
                backend.fit_null_model(phenotype, None, "binary")
            except RuntimeError as e:
                errors.append(str(e))

        t = threading.Thread(target=run_in_thread)
        t.start()
        t.join()

        assert len(errors) == 1
        assert "main thread" in errors[0].lower()

    def test_test_gene_thread_safety(self):
        """test_gene raises RuntimeError from non-main thread."""
        backend = _make_backend_with_mocked_env()
        null_model = _make_null_model(trait_type="binary")
        geno = np.zeros((10, 3))
        errors: list[str] = []

        def run_in_thread() -> None:
            try:
                backend.test_gene("GENE", geno, null_model, "SKAT", (1.0, 25.0))
            except RuntimeError as e:
                errors.append(str(e))

        t = threading.Thread(target=run_in_thread)
        t.start()
        t.join()

        assert len(errors) == 1
        assert "main thread" in errors[0].lower()
