"""Tests for deterministic Trivy report summarization."""

import json
from pathlib import Path

import pytest

from scripts.summarize_trivy import (
    has_fix,
    main,
    render_markdown,
    severity_counts,
    vulnerabilities,
)


def _document(*items: dict[str, object]) -> dict[str, object]:
    return {"Results": [{"Vulnerabilities": list(items)}]}


def test_vulnerabilities_flattens_multiple_results_and_skips_absent_list() -> None:
    first = {"VulnerabilityID": "CVE-1"}
    second = {"VulnerabilityID": "CVE-2"}
    document = {
        "Results": [
            {"Vulnerabilities": [first]},
            {"Target": "empty-result"},
            {"Vulnerabilities": [second]},
        ]
    }

    assert vulnerabilities(document) == [first, second]


@pytest.mark.parametrize(
    ("fixed_version", "expected"),
    [(None, False), ("", False), ("1.2.3", True)],
)
def test_has_fix_requires_nonempty_string(fixed_version: object, expected: bool) -> None:
    vulnerability = {}
    if fixed_version is not None:
        vulnerability["FixedVersion"] = fixed_version

    assert has_fix(vulnerability) is expected


def test_severity_counts_returns_every_bucket_and_maps_unknown_values() -> None:
    items = [
        {"Severity": "LOW"},
        {"Severity": "HIGH"},
        {"Severity": "not-a-trivy-severity"},
        {},
    ]

    assert severity_counts(items) == {
        "UNKNOWN": 2,
        "LOW": 1,
        "MEDIUM": 0,
        "HIGH": 1,
        "CRITICAL": 0,
    }


def test_render_markdown_partitions_complete_findings_by_fixed_version() -> None:
    complete = vulnerabilities(
        _document(
            {"Severity": "LOW", "FixedVersion": ""},
            {"Severity": "MEDIUM", "FixedVersion": "2.0"},
            {"Severity": "HIGH"},
        )
    )

    report = render_markdown(complete, [])

    assert report == (
        "## Trivy vulnerability summary\n"
        "\n"
        "- Complete findings: 3\n"
        "- Vendor-fixed findings: 1\n"
        "- Vendor-unfixed findings: 2\n"
        "- Actionable findings: 0\n"
        "\n"
        "### Complete findings by severity\n"
        "\n"
        "| Severity | Count |\n"
        "| --- | ---: |\n"
        "| UNKNOWN | 0 |\n"
        "| LOW | 1 |\n"
        "| MEDIUM | 1 |\n"
        "| HIGH | 1 |\n"
        "| CRITICAL | 0 |\n"
        "\n"
        "### Actionable findings by severity\n"
        "\n"
        "| Severity | Count |\n"
        "| --- | ---: |\n"
        "| UNKNOWN | 0 |\n"
        "| LOW | 0 |\n"
        "| MEDIUM | 0 |\n"
        "| HIGH | 0 |\n"
        "| CRITICAL | 0 |\n"
    )


def test_render_markdown_excludes_actionable_entries_without_nonempty_string_fix() -> None:
    actionable = [
        {"Severity": "CRITICAL", "FixedVersion": "9.9"},
        {"Severity": "HIGH", "FixedVersion": ""},
        {"Severity": "MEDIUM", "FixedVersion": 7},
        {"Severity": "LOW"},
    ]

    report = render_markdown([], actionable)

    assert "- Actionable findings: 1" in report
    actionable_section = report.split("### Actionable findings by severity", maxsplit=1)[1]
    assert "| UNKNOWN | 0 |" in actionable_section
    assert "| LOW | 0 |" in actionable_section
    assert "| MEDIUM | 0 |" in actionable_section
    assert "| HIGH | 0 |" in actionable_section
    assert "| CRITICAL | 1 |" in actionable_section


def test_main_counts_actionable_severity_but_never_becomes_the_gate(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    complete_path = tmp_path / "complete.json"
    actionable_path = tmp_path / "actionable.json"
    complete_path.write_text(
        json.dumps(_document({"Severity": "CRITICAL", "FixedVersion": "9.9"})),
        encoding="utf-8",
    )
    actionable_path.write_text(
        json.dumps(_document({"Severity": "CRITICAL", "FixedVersion": "9.9"})),
        encoding="utf-8",
    )

    result = main([str(complete_path), str(actionable_path)])

    assert result == 0
    output = capsys.readouterr().out
    assert "Actionable findings: 1" in output
    actionable_section = output.split("### Actionable findings by severity", maxsplit=1)[1]
    assert "| CRITICAL | 1 |" in actionable_section


@pytest.mark.parametrize(
    ("payload", "expected_error"),
    [
        ("[]", "top-level JSON value must be an object"),
        ("{not-json", "invalid JSON"),
        ('{"Results": NaN}', "invalid JSON"),
        ('{"Results": Infinity}', "invalid JSON"),
        ('{"Results": -Infinity}', "invalid JSON"),
    ],
)
def test_main_rejects_malformed_documents(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    payload: str,
    expected_error: str,
) -> None:
    complete_path = tmp_path / "complete.json"
    actionable_path = tmp_path / "actionable.json"
    complete_path.write_text(payload, encoding="utf-8")
    actionable_path.write_text(json.dumps(_document()), encoding="utf-8")

    result = main([str(complete_path), str(actionable_path)])

    assert result == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "error:" in captured.err
    assert expected_error in captured.err
