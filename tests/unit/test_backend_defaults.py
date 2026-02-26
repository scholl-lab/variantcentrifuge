# File: tests/unit/test_backend_defaults.py
# Location: tests/unit/test_backend_defaults.py
"""
Unit tests for Python-default backend behavior (Phase 25).

Verifies:
- get_skat_backend("auto") returns PythonSKATBackend
- get_skat_backend("python") returns PythonSKATBackend
- AssociationConfig() default fields are "python"
- RSKATTest() emits DeprecationWarning
- COASTTest() emits DeprecationWarning
"""

from __future__ import annotations

import warnings

import pytest

from variantcentrifuge.association.backends import get_skat_backend
from variantcentrifuge.association.base import AssociationConfig


@pytest.mark.unit
class TestGetSKATBackendDefaults:
    """Tests for get_skat_backend() auto/python behavior."""

    def test_get_skat_backend_auto_returns_python(self):
        """auto backend selector must return PythonSKATBackend (not RSKATBackend)."""
        backend = get_skat_backend("auto")
        assert type(backend).__name__ == "PythonSKATBackend"

    def test_get_skat_backend_python_explicit(self):
        """Explicit 'python' selector returns PythonSKATBackend."""
        backend = get_skat_backend("python")
        assert type(backend).__name__ == "PythonSKATBackend"

    def test_get_skat_backend_auto_and_python_same_type(self):
        """'auto' and 'python' return the same backend type."""
        auto_backend = get_skat_backend("auto")
        python_backend = get_skat_backend("python")
        assert type(auto_backend) is type(python_backend)

    def test_get_skat_backend_r_returns_r_backend(self):
        """Explicit 'r' selector returns RSKATBackend."""
        backend = get_skat_backend("r")
        assert type(backend).__name__ == "RSKATBackend"

    def test_get_skat_backend_invalid_raises_value_error(self):
        """Unknown backend name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown SKAT backend"):
            get_skat_backend("julia")


@pytest.mark.unit
class TestAssociationConfigDefaults:
    """Tests for AssociationConfig default field values."""

    def test_config_skat_backend_default_is_python(self):
        """AssociationConfig().skat_backend defaults to 'python'."""
        config = AssociationConfig()
        assert config.skat_backend == "python"

    def test_config_coast_backend_default_is_python(self):
        """AssociationConfig().coast_backend defaults to 'python'."""
        config = AssociationConfig()
        assert config.coast_backend == "python"

    def test_config_skat_backend_can_be_set_to_r(self):
        """skat_backend can be set to 'r' for explicit R backend use."""
        config = AssociationConfig(skat_backend="r")
        assert config.skat_backend == "r"

    def test_config_skat_backend_auto_still_accepted(self):
        """'auto' is still a valid value for skat_backend (backward compat)."""
        config = AssociationConfig(skat_backend="auto")
        assert config.skat_backend == "auto"


@pytest.mark.unit
class TestRSKATTestDeprecationWarning:
    """Tests for RSKATTest deprecation warning."""

    def test_rskat_test_emits_deprecation_warning(self):
        """RSKATTest() must emit DeprecationWarning on instantiation."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            RSKATTest()
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "RSKATTest" in str(w[0].message)
            assert "v0.17.0" in str(w[0].message)
            assert "python" in str(w[0].message)

    def test_rskat_test_deprecation_warning_stacklevel(self):
        """RSKATTest deprecation warning points to caller (stacklevel=2)."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            RSKATTest()
            assert len(w) == 1
            # The filename in the warning should be this test file, not skat_r.py
            assert "test_backend_defaults" in w[0].filename


@pytest.mark.unit
class TestCOASTTestDeprecationWarning:
    """Tests for COASTTest deprecation warning."""

    def test_coast_test_emits_deprecation_warning(self):
        """COASTTest() must emit DeprecationWarning on instantiation."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            COASTTest()
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "COASTTest" in str(w[0].message)
            assert "v0.17.0" in str(w[0].message)
            assert "python" in str(w[0].message)

    def test_coast_test_deprecation_warning_stacklevel(self):
        """COASTTest deprecation warning points to caller (stacklevel=2)."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            COASTTest()
            assert len(w) == 1
            assert "test_backend_defaults" in w[0].filename
