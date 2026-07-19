#!/usr/bin/env python3
"""Render deterministic Markdown counts from complete and actionable Trivy JSON."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import cast

SEVERITIES = ("UNKNOWN", "LOW", "MEDIUM", "HIGH", "CRITICAL")


def vulnerabilities(document: dict[str, object]) -> list[dict[str, object]]:
    """Flatten vulnerability dictionaries from every Trivy result."""
    results = document.get("Results", [])
    if results is None:
        return []
    if not isinstance(results, list):
        raise ValueError("'Results' must be a list")

    flattened: list[dict[str, object]] = []
    for result in results:
        if not isinstance(result, dict):
            raise ValueError("each entry in 'Results' must be an object")
        result_vulnerabilities = result.get("Vulnerabilities", [])
        if result_vulnerabilities is None:
            continue
        if not isinstance(result_vulnerabilities, list):
            raise ValueError("'Vulnerabilities' must be a list")
        for vulnerability in result_vulnerabilities:
            if not isinstance(vulnerability, dict):
                raise ValueError("each vulnerability must be an object")
            flattened.append(cast(dict[str, object], vulnerability))
    return flattened


def has_fix(vulnerability: dict[str, object]) -> bool:
    """Return whether Trivy reports a nonempty vendor-provided fixed version."""
    fixed_version = vulnerability.get("FixedVersion")
    return isinstance(fixed_version, str) and bool(fixed_version)


def severity_counts(items: list[dict[str, object]]) -> dict[str, int]:
    """Count vulnerabilities in deterministic Trivy severity buckets."""
    counts = dict.fromkeys(SEVERITIES, 0)
    for item in items:
        severity = item.get("Severity")
        normalized = severity.upper() if isinstance(severity, str) else "UNKNOWN"
        bucket = normalized if normalized in counts else "UNKNOWN"
        counts[bucket] += 1
    return counts


def render_markdown(complete: list[dict[str, object]], actionable: list[dict[str, object]]) -> str:
    """Render complete, vendor-fixed, vendor-unfixed, and actionable finding counts."""
    fixed = [vulnerability for vulnerability in complete if has_fix(vulnerability)]
    unfixed = [vulnerability for vulnerability in complete if not has_fix(vulnerability)]
    complete_by_severity = severity_counts(complete)
    actionable_by_severity = severity_counts(actionable)

    lines = [
        "## Trivy vulnerability summary",
        "",
        f"- Complete findings: {len(complete)}",
        f"- Vendor-fixed findings: {len(fixed)}",
        f"- Vendor-unfixed findings: {len(unfixed)}",
        f"- Actionable findings: {len(actionable)}",
        "",
        "### Complete findings by severity",
        "",
        "| Severity | Count |",
        "| --- | ---: |",
    ]
    lines.extend(f"| {severity} | {complete_by_severity[severity]} |" for severity in SEVERITIES)
    lines.extend(
        [
            "",
            "### Actionable findings by severity",
            "",
            "| Severity | Count |",
            "| --- | ---: |",
        ]
    )
    lines.extend(f"| {severity} | {actionable_by_severity[severity]} |" for severity in SEVERITIES)
    return "\n".join(lines) + "\n"


def _load_document(path: str) -> dict[str, object]:
    try:
        parsed = json.loads(Path(path).read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in {path}: {error.msg}") from error
    except OSError as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    if not isinstance(parsed, dict):
        raise ValueError(f"top-level JSON value must be an object in {path}")
    return cast(dict[str, object], parsed)


def main(argv: list[str] | None = None) -> int:
    """Read two Trivy JSON documents and print their Markdown summary."""
    arguments = sys.argv[1:] if argv is None else argv
    if len(arguments) != 2:
        print("error: expected COMPLETE_JSON ACTIONABLE_JSON", file=sys.stderr)
        return 2

    try:
        complete = vulnerabilities(_load_document(arguments[0]))
        actionable = vulnerabilities(_load_document(arguments[1]))
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(render_markdown(complete, actionable), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
