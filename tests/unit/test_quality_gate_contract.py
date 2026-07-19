"""Static contract tests for repository-wide local quality gates."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_makefile_quality_targets_cover_the_repository() -> None:
    """Lint and formatting targets must inspect every tracked Python surface."""
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert "$(PYTHON) -m ruff check ." in makefile
    assert "$(PYTHON) -m ruff format ." in makefile
    assert "$(PYTHON) -m ruff format --check --diff ." in makefile
    assert "lint: ## Run ruff linter (same as CI)" not in makefile
    assert "format-check: ## Check code formatting (same as CI)" not in makefile
    assert "mirroring GitHub Actions" not in makefile


def test_makefile_typecheck_is_blocking() -> None:
    """A mypy failure must propagate through ``make typecheck`` and ``ci-check``."""
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert "\t$(PYTHON) -m mypy variantcentrifuge/" in makefile
    assert "\t-$(PYTHON) -m mypy variantcentrifuge/" not in makefile
    assert "non-blocking" not in makefile


def test_dev_dependencies_keep_numpy_stubs_compatible_with_python_310() -> None:
    """Dev installs must use NumPy stubs parsable by the configured mypy target."""
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert '"numpy<2.5",' in pyproject


def test_runtime_dependencies_exclude_incompatible_smart_open_8() -> None:
    """Fresh installs must retain the ``smart_open`` callable used by the application."""
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert '"smart-open<8",' in pyproject
