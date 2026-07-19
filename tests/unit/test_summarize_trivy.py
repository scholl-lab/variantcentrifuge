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

    assert "Complete findings: 3" in report
    assert "Vendor-fixed findings: 1" in report
    assert "Vendor-unfixed findings: 2" in report
    assert "Actionable findings: 0" in report


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
