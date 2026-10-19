"""Regression tests for Issue #167:
- Organize variantcentrifuge --help into task-oriented argument groups
- Demote deprecated/internal flags to argparse.SUPPRESS or Advanced group
- Add curated examples epilog
"""

from __future__ import annotations

import pytest

from variantcentrifuge.cli import create_parser


@pytest.mark.unit
class TestIssue167HelpOrganization:
    """Verify task-oriented argument groups and examples epilog in --help."""

    def test_help_contains_task_oriented_groups(self) -> None:
        """--help output must organize arguments into task-oriented groups."""
        parser = create_parser()
        help_text = parser.format_help()

        expected_group_keywords = [
            "Core Input/Output",
            "Gene",
            "Filtering",
            "Inheritance",
            "Reporting & Visualization",
            "Checkpoint & Resume",
            "Privacy",
            "Performance",
            "Advanced & Tuning Options",
        ]

        for keyword in expected_group_keywords:
            assert any(
                keyword.lower() in group.title.lower()
                for group in parser._action_groups
                if group.title
            ), f"Expected group matching '{keyword}' in argument groups"
            assert keyword.lower() in help_text.lower(), f"Expected '{keyword}' in help text"

    def test_help_contains_curated_examples_epilog(self) -> None:
        """--help output must include a curated examples epilog."""
        parser = create_parser()
        help_text = parser.format_help()

        assert "Examples:" in help_text
        assert "variantcentrifuge -g" in help_text
        assert "--preset rare,coding" in help_text

    def test_internal_tuning_flags_demoted_to_advanced_group(self) -> None:
        """Internal tuning flags must be located in the Advanced/Internal group."""
        parser = create_parser()
        advanced_groups = [g for g in parser._action_groups if g.title and "Advanced" in g.title]
        assert len(advanced_groups) == 1, "Expected exactly one Advanced/Internal options group"

        advanced_group = advanced_groups[0]
        advanced_flags = {
            opt for action in advanced_group._group_actions for opt in action.option_strings
        }

        # Check internal knobs are in the advanced group
        internal_knobs = [
            "--force-chunked-processing",
            "--force-inheritance-processing",
            "--memory-safety-factor",
            "--inheritance-memory-fraction",
            "--sample-column-creation-method",
        ]
        for knob in internal_knobs:
            assert knob in advanced_flags, f"Expected internal flag '{knob}' in Advanced group"

    def test_deprecated_flags_are_suppressed_from_help(self) -> None:
        """Deprecated flags must be hidden from --help but remain functional."""
        parser = create_parser()
        help_text = parser.format_help()

        # Deprecated aliases should NOT appear in help
        assert "--pca-file" not in help_text
        assert "--pca-tool" not in help_text
        assert "--ped-file" not in help_text

        # But they must still parse without error (backward compatibility)
        parsed_pca = parser.parse_args(["-v", "test.vcf", "--pca-file", "my.eigenvec"])
        assert parsed_pca.pca == "my.eigenvec"

        parsed_ped = parser.parse_args(["-v", "test.vcf", "--ped-file", "family.ped"])
        assert parsed_ped.ped == "family.ped"
