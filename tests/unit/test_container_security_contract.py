from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def _contains_dependency_line(environment: str, expected: str) -> bool:
    return expected in environment.splitlines()


@pytest.mark.parametrize(
    "environment",
    [
        "#  - snpeff=5.2=hdfd78af_3\n",
        "    - snpeff=5.2=hdfd78af_3\n",
    ],
)
def test_dependency_matching_rejects_inactive_lines(environment: str) -> None:
    assert not _contains_dependency_line(
        environment,
        "  - snpeff=5.2=hdfd78af_3",
    )


def test_docker_environment_pins_release_era_java_tool_builds() -> None:
    environment = _text("conda/environment-docker.yml")
    assert _contains_dependency_line(environment, "  - snpeff=5.2=hdfd78af_3")
    assert _contains_dependency_line(environment, "  - snpsift=5.2=hdfd78af_0")


def test_docker_environment_has_secure_pip_floor() -> None:
    environment = _text("conda/environment-docker.yml")
    assert _contains_dependency_line(environment, "  - pip>=26.1.2")
