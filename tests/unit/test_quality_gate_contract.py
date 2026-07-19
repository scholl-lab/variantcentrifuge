"""Static contract tests for repository-wide local quality gates."""

from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.10
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[2]


def _make_recipe(makefile: str, target: str) -> list[str]:
    """Return recipe commands belonging only to one Make target."""
    lines = makefile.splitlines()
    for _index, line in enumerate(lines):
        if line.startswith((" ", "\t", "#")):
            continue
        target_list, separator, _ = line.partition(":")
        if separator and target in target_list.split():
            break
    else:
        raise AssertionError(f"Make target not found: {target}")

    recipe = []
    for line in lines[_index + 1 :]:
        if line.startswith("\t"):
            recipe.append(line.removeprefix("\t").rstrip())
        elif line.strip() and not line.startswith("#"):
            break

    if not recipe:
        raise AssertionError(f"Make target has no recipe: {target}")
    return recipe


def _load_pyproject() -> dict[str, Any]:
    """Load project configuration using the runtime's TOML implementation."""
    with (ROOT / "pyproject.toml").open("rb") as pyproject_file:
        return tomllib.load(pyproject_file)


def test_makefile_quality_targets_cover_the_repository() -> None:
    """Lint and formatting targets must inspect every tracked Python surface."""
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert "$(PYTHON) -m ruff check ." in _make_recipe(makefile, "lint")
    assert "$(PYTHON) -m ruff format ." in _make_recipe(makefile, "format")
    assert "$(PYTHON) -m ruff format --check --diff ." in _make_recipe(makefile, "format-check")
    assert "lint: ## Run ruff linter (same as CI)" not in makefile
    assert "format-check: ## Check code formatting (same as CI)" not in makefile
    assert "mirroring GitHub Actions" not in makefile


def test_makefile_typecheck_is_blocking() -> None:
    """A mypy failure must propagate through ``make typecheck`` and ``ci-check``."""
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    typecheck_recipe = _make_recipe(makefile, "typecheck")
    mypy_commands = [command for command in typecheck_recipe if "-m mypy" in command]

    assert mypy_commands == ["$(PYTHON) -m mypy variantcentrifuge/"]
    assert not any(command.lstrip("@").startswith("-") for command in typecheck_recipe)
    assert not any("|| true" in command for command in typecheck_recipe)


def test_ci_check_invokes_every_local_quality_gate() -> None:
    """The aggregate target must delegate to each documented local gate."""
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    ci_recipe = {command.lstrip("@") for command in _make_recipe(makefile, "ci-check")}

    assert {
        "$(MAKE) lint",
        "$(MAKE) format-check",
        "$(MAKE) typecheck",
        "$(MAKE) test-fast",
    } <= ci_recipe


def test_dev_dependencies_keep_numpy_stubs_compatible_with_python_310() -> None:
    """Dev installs must use NumPy stubs parsable by the configured mypy target."""
    pyproject = _load_pyproject()

    assert "numpy<2.5" in pyproject["project"]["optional-dependencies"]["dev"]


def test_runtime_dependencies_exclude_incompatible_smart_open_8() -> None:
    """Fresh installs must retain the ``smart_open`` callable used by the application."""
    pyproject = _load_pyproject()

    assert "smart-open<8" in pyproject["project"]["dependencies"]
