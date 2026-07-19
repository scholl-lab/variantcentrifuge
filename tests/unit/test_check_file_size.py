"""Behavior tests for the production Python file-size guard."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CHECKER = REPOSITORY_ROOT / "scripts" / "check_file_size.py"


def _write_python_file(path: Path, line_count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join("pass" for _ in range(line_count)) + "\n")


def _run_checker(
    source_root: Path, allowlist: Path | None = None
) -> subprocess.CompletedProcess[str]:
    command = [sys.executable, str(CHECKER), str(source_root)]
    if allowlist is not None:
        command.extend(["--allowlist", str(allowlist)])
    return subprocess.run(command, check=False, capture_output=True, text=True)


def test_accepts_production_module_at_600_lines(tmp_path: Path) -> None:
    source_root = tmp_path / "package"
    _write_python_file(source_root / "module.py", 600)

    result = _run_checker(source_root)

    assert result.returncode == 0, result.stdout + result.stderr


def test_rejects_unallowlisted_production_module_over_600_lines(tmp_path: Path) -> None:
    source_root = tmp_path / "package"
    _write_python_file(source_root / "module.py", 601)

    result = _run_checker(source_root)

    assert result.returncode == 1
    assert "module.py: 601 lines exceeds the 600-line limit" in result.stdout


def test_accepts_allowlisted_module_at_its_baseline_ceiling(tmp_path: Path) -> None:
    source_root = tmp_path / "package"
    source_file = source_root / "module.py"
    _write_python_file(source_file, 601)
    allowlist = tmp_path / "allowlist"
    allowlist.write_text(f"{source_file}:601\n")

    result = _run_checker(source_root, allowlist)

    assert result.returncode == 0, result.stdout + result.stderr


def test_ignores_test_modules(tmp_path: Path) -> None:
    source_root = tmp_path / "package"
    _write_python_file(source_root / "tests" / "test_large.py", 601)

    result = _run_checker(source_root)

    assert result.returncode == 0, result.stdout + result.stderr
