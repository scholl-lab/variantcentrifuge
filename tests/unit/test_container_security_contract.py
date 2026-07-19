from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_docker_environment_pins_release_era_java_tool_builds() -> None:
    environment = _text("conda/environment-docker.yml")
    assert "  - snpeff=5.2=hdfd78af_3\n" in environment
    assert "  - snpsift=5.2=hdfd78af_0\n" in environment


def test_docker_environment_has_secure_pip_floor() -> None:
    environment = _text("conda/environment-docker.yml")
    assert "  - pip>=26.1.2\n" in environment
