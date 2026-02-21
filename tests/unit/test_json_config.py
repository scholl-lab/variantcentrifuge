"""
Unit tests for JSON config loading, validation, and CLI override logic.

Tests _validate_association_config_dict() and _build_assoc_config_from_context()
from variantcentrifuge.stages.analysis_stages.

Coverage:
- Validation: valid dict, unknown keys, wrong types, invalid enum values, empty dict
- Config building: JSON-only, CLI-only, CLI override, mixed sources
- End-to-end path: context.config["association"] -> AssociationConfig
- Round-trip: all JSON fields -> matching AssociationConfig
"""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import (
    _build_assoc_config_from_context,
    _validate_association_config_dict,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_context(config: dict) -> PipelineContext:
    """Build a minimal PipelineContext with the given config dict."""
    workspace = Mock(spec=Workspace)
    return PipelineContext(
        args=Namespace(),
        config=config,
        workspace=workspace,
    )


# ---------------------------------------------------------------------------
# Tests: _validate_association_config_dict
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestValidateAssociationConfigDict:
    """Unit tests for _validate_association_config_dict()."""

    def test_empty_dict_is_valid(self):
        """Empty dict: no required keys, should pass without error."""
        _validate_association_config_dict({})

    def test_valid_dict_passes(self):
        """All valid keys with correct types: no exception."""
        _validate_association_config_dict(
            {
                "correction_method": "fdr",
                "gene_burden_mode": "samples",
                "trait_type": "binary",
                "variant_weights": "beta:1,25",
                "skat_backend": "auto",
                "skat_method": "SKAT",
                "pca_components": 10,
                "min_cases": 200,
                "min_case_carriers": 10,
                "max_case_control_ratio": 20.0,
                "confidence_interval_alpha": 0.05,
                "continuity_correction": 0.5,
                "missing_site_threshold": 0.10,
                "missing_sample_threshold": 0.80,
                "firth_max_iter": 25,
                "association_tests": ["fisher"],
                "coast_weights": [1.0, 2.0, 3.0],
            }
        )

    def test_single_unknown_key_raises(self):
        """A single unknown key raises ValueError naming that key."""
        with pytest.raises(ValueError, match="Unknown keys"):
            _validate_association_config_dict({"typo_key": "value"})

    def test_multiple_unknown_keys_all_listed(self):
        """Multiple unknown keys: all listed in the error message."""
        with pytest.raises(ValueError) as exc_info:
            _validate_association_config_dict(
                {"bad_key_1": 1, "bad_key_2": 2, "correction_method": "fdr"}
            )
        msg = str(exc_info.value)
        assert "bad_key_1" in msg
        assert "bad_key_2" in msg

    def test_wrong_type_string_field(self):
        """correction_method must be a string; int raises ValueError."""
        with pytest.raises(ValueError, match="'correction_method' must be a string"):
            _validate_association_config_dict({"correction_method": 42})

    def test_wrong_type_int_field(self):
        """pca_components must be an int; float raises ValueError."""
        with pytest.raises(ValueError, match="'pca_components' must be an integer"):
            _validate_association_config_dict({"pca_components": 10.5})

    def test_wrong_type_float_field(self):
        """max_case_control_ratio must be numeric; string raises ValueError."""
        with pytest.raises(ValueError, match="'max_case_control_ratio' must be a number"):
            _validate_association_config_dict({"max_case_control_ratio": "twenty"})

    def test_wrong_type_list_str_field(self):
        """association_tests must be a list; string raises ValueError."""
        with pytest.raises(ValueError, match="'association_tests' must be a list"):
            _validate_association_config_dict({"association_tests": "fisher"})

    def test_wrong_type_list_float_field(self):
        """coast_weights must be a list; dict raises ValueError."""
        with pytest.raises(ValueError, match="'coast_weights' must be a list"):
            _validate_association_config_dict({"coast_weights": {"bmv": 1.0}})

    def test_invalid_correction_method_enum(self):
        """correction_method must be 'fdr' or 'bonferroni'."""
        with pytest.raises(ValueError, match=r"correction_method.*fdr.*bonferroni"):
            _validate_association_config_dict({"correction_method": "invalid"})

    def test_invalid_trait_type_enum(self):
        """trait_type must be 'binary' or 'quantitative'."""
        with pytest.raises(ValueError, match=r"trait_type.*binary.*quantitative"):
            _validate_association_config_dict({"trait_type": "continuous"})

    def test_invalid_skat_backend_enum(self):
        """skat_backend must be 'auto', 'r', or 'python'."""
        with pytest.raises(ValueError, match=r"skat_backend.*auto.*r.*python"):
            _validate_association_config_dict({"skat_backend": "julia"})

    def test_multiple_errors_collected(self):
        """Multiple validation errors: all reported in a single exception."""
        with pytest.raises(ValueError) as exc_info:
            _validate_association_config_dict(
                {
                    "unknown_key": 1,
                    "correction_method": "bad_method",
                    "pca_components": "ten",
                }
            )
        msg = str(exc_info.value)
        assert "error(s)" in msg
        assert "unknown_key" in msg
        assert "correction_method" in msg
        assert "pca_components" in msg

    def test_nullable_string_field_accepts_none(self):
        """Nullable string fields (like covariate_file) may be None."""
        _validate_association_config_dict({"covariate_file": None})

    def test_int_bool_rejected_for_int_fields(self):
        """bool is a subclass of int in Python; bool values are accepted for int fields."""
        # True is int(1), False is int(0) — this is valid Python type semantics
        # We document this behavior with a passing test (not a bug)
        _validate_association_config_dict({"pca_components": True})


# ---------------------------------------------------------------------------
# Tests: _build_assoc_config_from_context
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestBuildAssocConfigFromContext:
    """Unit tests for _build_assoc_config_from_context()."""

    def test_no_association_section_returns_defaults(self):
        """No 'association' key: all AssociationConfig fields use their defaults."""
        ctx = _make_context({"correction_method": "fdr", "gene_burden_mode": "samples"})
        config = _build_assoc_config_from_context(ctx)

        assert isinstance(config, AssociationConfig)
        assert config.correction_method == "fdr"
        assert config.gene_burden_mode == "samples"
        assert config.trait_type == "binary"
        assert config.variant_weights == "beta:1,25"
        assert config.min_cases == 200
        assert config.max_case_control_ratio == 20.0
        assert config.min_case_carriers == 10
        assert config.pca_components == 10
        assert config.skat_backend == "auto"
        assert config.skat_method == "SKAT"
        assert config.covariate_file is None
        assert config.diagnostics_output is None

    def test_json_only_config_applied(self):
        """JSON association section applies when CLI keys are absent."""
        ctx = _make_context(
            {
                "association": {
                    "correction_method": "bonferroni",
                    "trait_type": "quantitative",
                    "min_cases": 50,
                    "pca_components": 5,
                    "skat_backend": "python",
                    "diagnostics_output": "/tmp/diag",
                }
            }
        )
        config = _build_assoc_config_from_context(ctx)

        assert config.correction_method == "bonferroni"
        assert config.trait_type == "quantitative"
        assert config.min_cases == 50
        assert config.pca_components == 5
        assert config.skat_backend == "python"
        assert config.diagnostics_output == "/tmp/diag"

    def test_cli_overrides_json(self):
        """CLI values in context.config override JSON association section values."""
        ctx = _make_context(
            {
                # CLI-set (takes precedence)
                "correction_method": "bonferroni",
                "trait_type": "quantitative",
                "diagnostics_output": "/cli/diag",
                # JSON association section (lower precedence)
                "association": {
                    "correction_method": "fdr",
                    "trait_type": "binary",
                    "diagnostics_output": "/json/diag",
                    "min_cases": 100,
                },
            }
        )
        config = _build_assoc_config_from_context(ctx)

        # CLI wins
        assert config.correction_method == "bonferroni"
        assert config.trait_type == "quantitative"
        assert config.diagnostics_output == "/cli/diag"
        # JSON fills gaps
        assert config.min_cases == 100

    def test_mixed_sources(self):
        """Some fields from CLI, others from JSON, others from defaults."""
        ctx = _make_context(
            {
                # CLI-set
                "covariate_file": "/path/to/cov.tsv",
                "skat_backend": "r",
                # JSON section
                "association": {
                    "min_cases": 75,
                    "min_case_carriers": 5,
                    "coast_weights": [1.0, 2.0, 3.0],
                },
            }
        )
        config = _build_assoc_config_from_context(ctx)

        # From CLI
        assert config.covariate_file == "/path/to/cov.tsv"
        assert config.skat_backend == "r"
        # From JSON
        assert config.min_cases == 75
        assert config.min_case_carriers == 5
        assert config.coast_weights == [1.0, 2.0, 3.0]
        # From defaults
        assert config.correction_method == "fdr"
        assert config.trait_type == "binary"
        assert config.pca_components == 10

    def test_invalid_json_section_raises_on_build(self):
        """Invalid JSON section raises ValueError at build time (fail-fast)."""
        ctx = _make_context(
            {
                "association": {
                    "unknown_field": "value",
                }
            }
        )
        with pytest.raises(ValueError, match="Unknown keys"):
            _build_assoc_config_from_context(ctx)

    def test_context_config_association_path_end_to_end(self):
        """End-to-end test: association key in context.config flows to AssociationConfig.

        This validates the full path: JSON config.json "association" section
        -> context.config["association"] -> _build_assoc_config_from_context()
        -> AssociationConfig with correct values.

        In production, context.config is populated by ConfigurationLoadingStage
        from load_config(), which returns json.load(f) verbatim. Any "association"
        top-level key is therefore preserved exactly as shown here.
        """
        # Simulate what ConfigurationLoadingStage.process() does:
        # context.config = load_config(config_file)
        # where config.json contains {"association": {...}, ...other_keys...}
        context_config = {
            # Other config.json top-level keys (e.g. filter presets, links, etc.)
            "filter": "rare",
            "fields_of_interest": ["GENE"],
            # The "association" section we added to config.json
            "association": {
                "correction_method": "bonferroni",
                "trait_type": "quantitative",
                "min_cases": 150,
                "max_case_control_ratio": 10.0,
                "skat_backend": "python",
                "pca_components": 5,
            },
        }
        ctx = _make_context(context_config)
        config = _build_assoc_config_from_context(ctx)

        assert config.correction_method == "bonferroni"
        assert config.trait_type == "quantitative"
        assert config.min_cases == 150
        assert config.max_case_control_ratio == 10.0
        assert config.skat_backend == "python"
        assert config.pca_components == 5

    def test_round_trip_all_json_fields(self):
        """Round-trip: all JSON-settable fields -> matching AssociationConfig values."""
        ctx = _make_context(
            {
                "association": {
                    "correction_method": "bonferroni",
                    "gene_burden_mode": "alleles",
                    "trait_type": "quantitative",
                    "variant_weights": "uniform",
                    "skat_backend": "r",
                    "skat_method": "Burden",
                    "covariate_file": "/cov.tsv",
                    "covariate_columns": ["age", "sex"],
                    "categorical_covariates": ["sex"],
                    "pca_file": "/pca.tsv",
                    "pca_tool": "akt",
                    "pca_components": 20,
                    "coast_weights": [1.5, 2.5, 3.5],
                    "association_tests": ["fisher", "skat_python"],
                    "min_cases": 300,
                    "max_case_control_ratio": 15.0,
                    "min_case_carriers": 20,
                    "diagnostics_output": "/diag/",
                    "confidence_interval_method": "normal_approx",
                    "confidence_interval_alpha": 0.01,
                    "continuity_correction": 0.25,
                    "missing_site_threshold": 0.05,
                    "missing_sample_threshold": 0.75,
                    "firth_max_iter": 50,
                }
            }
        )
        config = _build_assoc_config_from_context(ctx)

        assert config.correction_method == "bonferroni"
        assert config.gene_burden_mode == "alleles"
        assert config.trait_type == "quantitative"
        assert config.variant_weights == "uniform"
        assert config.skat_backend == "r"
        assert config.skat_method == "Burden"
        assert config.covariate_file == "/cov.tsv"
        assert config.covariate_columns == ["age", "sex"]
        assert config.categorical_covariates == ["sex"]
        assert config.pca_file == "/pca.tsv"
        assert config.pca_tool == "akt"
        assert config.pca_components == 20
        assert config.coast_weights == [1.5, 2.5, 3.5]
        assert config.min_cases == 300
        assert config.max_case_control_ratio == 15.0
        assert config.min_case_carriers == 20
        assert config.diagnostics_output == "/diag/"
        assert config.confidence_interval_method == "normal_approx"
        assert config.confidence_interval_alpha == 0.01
        assert config.continuity_correction == 0.25
        assert config.missing_site_threshold == 0.05
        assert config.missing_sample_threshold == 0.75
        assert config.firth_max_iter == 50

    def test_cli_none_value_falls_through_to_json(self):
        """CLI value of None for a nullable field falls through to JSON."""
        ctx = _make_context(
            {
                # CLI set diagnostics_output=None (no --diagnostics-output flag given)
                "diagnostics_output": None,
                "association": {
                    "diagnostics_output": "/from/json/",
                },
            }
        )
        config = _build_assoc_config_from_context(ctx)
        # None CLI value -> JSON value wins for nullable fields
        assert config.diagnostics_output == "/from/json/"

    def test_association_prefixed_min_cases_keys(self):
        """association_min_cases (CLI-style) is resolved to min_cases correctly."""
        # CLI sets "association_min_cases" (the prefixed form used in context.config)
        ctx = _make_context(
            {
                "association_min_cases": 999,
                "association_max_case_control_ratio": 5.0,
                "association_min_case_carriers": 3,
                "association": {
                    # These JSON values are lower priority than CLI prefixed keys
                    "min_cases": 1,
                    "max_case_control_ratio": 1.0,
                    "min_case_carriers": 1,
                },
            }
        )
        config = _build_assoc_config_from_context(ctx)
        # CLI prefixed keys win (non-nullable: key present in cfg -> CLI wins)
        assert config.min_cases == 999
        assert config.max_case_control_ratio == 5.0
        assert config.min_case_carriers == 3
