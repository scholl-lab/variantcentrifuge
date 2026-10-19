"""Static contract tests for repository-wide local quality gates."""

import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

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


def test_pytest_disables_benchmark_by_default() -> None:
    """Benchmark calibration must be disabled in default pytest options."""
    pyproject = _load_pyproject()
    addopts = pyproject["tool"]["pytest"]["ini_options"]["addopts"]

    assert "--benchmark-disable" in addopts


def test_test_result_class_is_exempt_from_pytest_collection() -> None:
    """TestResult is a data class and must not trigger PytestCollectionWarning."""
    from variantcentrifuge.association.base import TestResult

    assert getattr(TestResult, "__test__", True) is False


def test_comp_het_original_benchmark_is_removed() -> None:
    """Benchmark suite must not contain misleading non-existent original implementation test."""
    benchmark_code = (ROOT / "tests/performance/benchmark_comp_het.py").read_text(encoding="utf-8")

    assert "test_comp_het_original_scaling" not in benchmark_code
    assert "original implementations" not in benchmark_code


def test_makefile_installs_use_locked_dependencies() -> None:
    """Makefile install targets must enforce committed lockfile."""
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert any("uv sync --locked" in cmd for cmd in _make_recipe(makefile, "install"))
    assert any(
        "uv sync --locked --extra dev" in cmd for cmd in _make_recipe(makefile, "install-dev")
    )


def test_ci_workflow_uses_locked_installations() -> None:
    """GitHub Actions test workflow must install using uv sync --locked."""
    test_yml = (ROOT / ".github/workflows/test.yml").read_text(encoding="utf-8")

    assert "uv sync --locked --extra dev" in test_yml
    assert "uv pip install" not in test_yml


def test_release_workflow_contract() -> None:
    """Release workflow must publish only on release events using PyPI Trusted Publishing."""
    release_yml_path = ROOT / ".github/workflows/release.yml"
    assert release_yml_path.exists(), "release.yml workflow must exist"

    release_yml = release_yml_path.read_text(encoding="utf-8")
    assert "id-token: write" in release_yml
    assert "pypa/gh-action-pypi-publish" in release_yml
    assert "check_version_sync.py" in release_yml
    assert "pull_request" not in release_yml


def test_ruff_per_file_ignores_has_no_stale_exemptions() -> None:
    """Ruff per-file-ignores must not exempt F401/F811/F841 outside __init__.py."""
    pyproject = _load_pyproject()
    per_file_ignores = pyproject["tool"]["ruff"]["lint"]["per-file-ignores"]

    for pattern, codes in per_file_ignores.items():
        if pattern == "__init__.py":
            assert codes == ["F401"]
            continue
        assert "F401" not in codes, f"{pattern} has stale F401 exemption"
        assert "F811" not in codes, f"{pattern} has stale F811 exemption"
        assert "F841" not in codes, f"{pattern} has stale F841 exemption"


def test_mypy_configuration_disallows_untyped_defs() -> None:
    """Mypy configuration must disallow untyped defs and restrict missing import ignores."""
    pyproject = _load_pyproject()
    mypy_conf = pyproject["tool"]["mypy"]

    assert mypy_conf.get("disallow_untyped_defs") is True
    assert mypy_conf.get("ignore_missing_imports") is not True


def test_coverage_report_defines_fail_under() -> None:
    """Coverage configuration must enforce a positive coverage floor of 75%."""
    pyproject = _load_pyproject()
    fail_under = pyproject["tool"]["coverage"]["report"]["fail_under"]

    assert fail_under == 75


def test_pre_commit_contains_mypy_hook() -> None:
    """Pre-commit configuration must include local mypy quality hook."""
    pre_commit = (ROOT / ".pre-commit-config.yaml").read_text(encoding="utf-8")
    assert "id: mypy" in pre_commit


def test_pre_commit_mypy_uses_locked_virtualenv() -> None:
    """Pre-commit mypy hook must route through locked virtualenv or uv run."""
    pre_commit = (ROOT / ".pre-commit-config.yaml").read_text(encoding="utf-8")
    data = yaml.safe_load(pre_commit)
    local_repo = next(r for r in data["repos"] if r.get("repo") == "local")
    mypy_hook = next(h for h in local_repo["hooks"] if h.get("id") == "mypy")
    entry = mypy_hook["entry"]
    assert ".venv/bin/python" in entry
    assert "uv run --locked --extra dev" in entry


def test_makefile_defaults_to_virtualenv_python() -> None:
    """Makefile must default PYTHON to locked .venv if present."""
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    assert "PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python)" in makefile


def test_ci_workflow_mypy_is_blocking() -> None:
    """GitHub Actions mypy step must be blocking without continue-on-error."""
    test_yml = (ROOT / ".github/workflows/test.yml").read_text(encoding="utf-8")
    data = yaml.safe_load(test_yml)
    type_check_steps = data["jobs"]["type-check"]["steps"]
    mypy_step = next(s for s in type_check_steps if "python -m mypy" in s.get("run", ""))
    assert mypy_step.get("continue-on-error") is not True


def test_ci_workflow_matrix_covers_python_versions() -> None:
    """GitHub Actions test matrix must cover Python 3.10, 3.11, 3.12, and 3.13."""
    test_yml = (ROOT / ".github/workflows/test.yml").read_text(encoding="utf-8")
    data = yaml.safe_load(test_yml)
    matrix_versions = data["jobs"]["test"]["strategy"]["matrix"]["python-version"]
    for py_ver in ["3.10", "3.11", "3.12", "3.13"]:
        assert py_ver in matrix_versions, f"Missing Python {py_ver} in test matrix"


def test_ci_workflow_runs_mocked_integration_tests() -> None:
    """CI workflow must define and strictly require a mocked integration test job."""
    test_yml = (ROOT / ".github/workflows/test.yml").read_text(encoding="utf-8")
    data = yaml.safe_load(test_yml)

    # Verify integration-test job
    integration_job = data["jobs"]["integration-test"]
    assert any("integration and not slow" in s.get("run", "") for s in integration_job["steps"])

    # Verify ci-success gate requirements
    ci_success = data["jobs"]["ci-success"]
    assert set(ci_success["needs"]) == {"lint", "type-check", "test", "integration-test"}

    # Simulate ci-success condition script in bash
    check_step = next(
        s for s in ci_success["steps"] if "All required jobs succeeded!" in s.get("run", "")
    )
    raw_script = check_step["run"]

    # When all jobs succeed
    all_success_script = (
        raw_script.replace("${{ needs.lint.result }}", "success")
        .replace("${{ needs.type-check.result }}", "success")
        .replace("${{ needs.test.result }}", "success")
        .replace("${{ needs.integration-test.result }}", "success")
    )
    res = subprocess.run(["bash", "-c", all_success_script], capture_output=True, text=True)
    assert res.returncode == 0, f"Expected 0 when all succeed: {res.stderr}"

    # When any job fails, cancels, or is skipped, it must exit non-zero
    for bad_job in ["lint", "type-check", "test", "integration-test"]:
        for bad_status in ["failure", "cancelled", "skipped"]:
            job_script = (
                raw_script.replace(f"${{{{ needs.{bad_job}.result }}}}", bad_status)
                .replace("${{ needs.lint.result }}", "success")
                .replace("${{ needs.type-check.result }}", "success")
                .replace("${{ needs.test.result }}", "success")
                .replace("${{ needs.integration-test.result }}", "success")
            )
            res = subprocess.run(["bash", "-c", job_script], capture_output=True, text=True)
            assert res.returncode != 0, f"Expected non-zero when {bad_job} is {bad_status}"


def test_installation_docs_document_apptainer() -> None:
    """Installation docs must cover Apptainer/Singularity execution.

    Must include binds, SnpEff storage, and entrypoint activation.
    """
    install_md = (ROOT / "docs" / "source" / "installation.md").read_text(encoding="utf-8")
    assert "Apptainer" in install_md
    assert "Singularity" in install_md
    assert "docker://ghcr.io/scholl-lab/variantcentrifuge" in install_md
    assert "-B " in install_md
    assert "/data/snpeff" in install_md
    assert "sha256:" in install_md
    assert "apptainer run" in install_md
    assert "/usr/local/bin/_entrypoint.sh" in install_md
    # Ensure bare apptainer exec ... variantcentrifuge without entrypoint is not instructed
    assert "apptainer exec docker://" not in install_md


def test_release_workflow_isolated_wheel_installation() -> None:
    """Release workflow must verify wheel in isolated environment outside source checkout."""
    release_yml = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    data = yaml.safe_load(release_yml)
    steps = data["jobs"]["build"]["steps"]
    verify_step = next(s for s in steps if s.get("name") == "Verify isolated installation of wheel")
    run_cmd = verify_step["run"]
    assert "RUNNER_TEMP" in run_cmd
    assert "-I -c" in run_cmd
    assert "site-packages" in run_cmd
    assert "m.version('variantcentrifuge')" in run_cmd
    assert "vc.__version__ == exp" in run_cmd


def test_psutil_types_enforced_in_dev_and_mypy() -> None:
    """types-psutil must be in dev dependencies, excluded from mypy overrides, and enforced."""
    pyproject = _load_pyproject()
    dev_deps = pyproject["project"]["optional-dependencies"]["dev"]
    assert "types-psutil" in dev_deps

    overrides = pyproject["tool"]["mypy"]["overrides"]
    modules = [m for o in overrides for m in o.get("module", [])]
    assert "psutil" not in modules
    assert "psutil.*" not in modules

    # Verify comments exist explaining third-party overrides in pyproject.toml
    pyproject_text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "Third-party libraries without official PEP 561" in pyproject_text

    # Verify mypy type checks psutil attributes
    res = subprocess.run(
        [sys.executable, "-m", "mypy", "-c", "import psutil; psutil.nonexistent_attr"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert res.returncode != 0
    assert 'Module has no attribute "nonexistent_attr"' in res.stdout
