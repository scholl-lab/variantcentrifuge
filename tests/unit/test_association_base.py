"""
Unit tests for association framework base abstractions.

Tests TestResult dataclass, AssociationConfig defaults, and
AssociationTest ABC enforcement.
"""

from __future__ import annotations

import pytest

from variantcentrifuge.association.base import (
    AssociationConfig,
    AssociationTest,
    TestResult,
)

# ---------------------------------------------------------------------------
# TestResult tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestTestResult:
    """Tests for the TestResult dataclass."""

    def test_creation_all_fields(self):
        """TestResult can be created with all fields populated."""
        result = TestResult(
            gene="BRCA1",
            test_name="fisher",
            p_value=0.01,
            corrected_p_value=0.05,
            effect_size=2.5,
            ci_lower=1.1,
            ci_upper=5.7,
            n_cases=100,
            n_controls=200,
            n_variants=5,
        )
        assert result.gene == "BRCA1"
        assert result.test_name == "fisher"
        assert result.p_value == 0.01
        assert result.corrected_p_value == 0.05
        assert result.effect_size == 2.5
        assert result.ci_lower == 1.1
        assert result.ci_upper == 5.7
        assert result.n_cases == 100
        assert result.n_controls == 200
        assert result.n_variants == 5
        assert result.extra == {}

    def test_creation_p_value_none_skipped_gene(self):
        """TestResult with p_value=None represents a skipped gene (not failed)."""
        result = TestResult(
            gene="PCSK9",
            test_name="fisher",
            p_value=None,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            n_cases=50,
            n_controls=50,
            n_variants=0,
        )
        assert result.p_value is None
        assert result.corrected_p_value is None
        assert result.effect_size is None
        assert result.n_variants == 0

    def test_extra_dict_usage(self):
        """TestResult.extra stores test-specific ancillary data."""
        result = TestResult(
            gene="TP53",
            test_name="fisher",
            p_value=0.001,
            corrected_p_value=0.01,
            effect_size=3.0,
            ci_lower=1.5,
            ci_upper=6.0,
            n_cases=75,
            n_controls=75,
            n_variants=3,
            extra={"table": [[5, 1], [70, 74]], "method": "score"},
        )
        assert result.extra["table"] == [[5, 1], [70, 74]]
        assert result.extra["method"] == "score"

    def test_extra_default_is_empty_dict(self):
        """TestResult.extra defaults to empty dict (not shared mutable)."""
        r1 = TestResult(
            gene="A",
            test_name="fisher",
            p_value=0.5,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            n_cases=10,
            n_controls=10,
            n_variants=1,
        )
        r2 = TestResult(
            gene="B",
            test_name="fisher",
            p_value=0.5,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            n_cases=10,
            n_controls=10,
            n_variants=1,
        )
        # Mutating one does not affect the other
        r1.extra["key"] = "val"
        assert "key" not in r2.extra

    def test_corrected_p_value_none_until_engine_fills(self):
        """corrected_p_value=None on construction; engine fills it post-correction."""
        result = TestResult(
            gene="MYH7",
            test_name="fisher",
            p_value=0.03,
            corrected_p_value=None,
            effect_size=1.8,
            ci_lower=0.9,
            ci_upper=3.6,
            n_cases=120,
            n_controls=120,
            n_variants=2,
        )
        assert result.corrected_p_value is None
        # Engine sets it later
        result.corrected_p_value = 0.15
        assert result.corrected_p_value == 0.15


# ---------------------------------------------------------------------------
# AssociationConfig tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationConfig:
    """Tests for AssociationConfig defaults."""

    def test_defaults_match_gene_burden_expectations(self):
        """Default values mirror gene_burden.py cfg keys (lines 373-377)."""
        cfg = AssociationConfig()
        assert cfg.correction_method == "fdr"
        assert cfg.gene_burden_mode == "samples"
        assert cfg.confidence_interval_method == "normal_approx"
        assert cfg.confidence_interval_alpha == 0.05
        assert cfg.continuity_correction == 0.5

    def test_custom_values_accepted(self):
        """AssociationConfig accepts non-default values."""
        cfg = AssociationConfig(
            correction_method="bonferroni",
            gene_burden_mode="alleles",
            confidence_interval_alpha=0.01,
            continuity_correction=0.0,
        )
        assert cfg.correction_method == "bonferroni"
        assert cfg.gene_burden_mode == "alleles"
        assert cfg.confidence_interval_alpha == 0.01
        assert cfg.continuity_correction == 0.0

    def test_confidence_interval_alpha_zero_point_zero_five(self):
        """Default alpha=0.05 gives 95% CIs."""
        cfg = AssociationConfig()
        assert cfg.confidence_interval_alpha == pytest.approx(0.05)

    def test_continuity_correction_haldane_anscombe(self):
        """Default continuity_correction=0.5 is the Haldane-Anscombe value."""
        cfg = AssociationConfig()
        assert cfg.continuity_correction == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# AssociationTest ABC tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationTestABC:
    """Tests for the AssociationTest abstract base class."""

    def test_cannot_instantiate_directly(self):
        """AssociationTest is abstract and cannot be instantiated directly."""
        with pytest.raises(TypeError):
            AssociationTest()  # type: ignore[abstract]

    def test_subclass_must_implement_name(self):
        """Subclass missing name property raises TypeError on instantiation."""

        class MissingName(AssociationTest):
            def run(self, gene, contingency_data, config):  # type: ignore[override]
                return None

        with pytest.raises(TypeError):
            MissingName()

    def test_subclass_must_implement_run(self):
        """Subclass missing run method raises TypeError on instantiation."""

        class MissingRun(AssociationTest):
            @property
            def name(self) -> str:
                return "missing_run"

        with pytest.raises(TypeError):
            MissingRun()

    def test_concrete_subclass_instantiates_successfully(self):
        """A complete subclass with name and run instantiates successfully."""

        class ConcreteTest(AssociationTest):
            @property
            def name(self) -> str:
                return "concrete"

            def run(self, gene, contingency_data, config):  # type: ignore[override]
                return TestResult(
                    gene=gene,
                    test_name=self.name,
                    p_value=1.0,
                    corrected_p_value=None,
                    effect_size=None,
                    ci_lower=None,
                    ci_upper=None,
                    n_cases=0,
                    n_controls=0,
                    n_variants=0,
                )

        test = ConcreteTest()
        assert test.name == "concrete"

    def test_check_dependencies_default_noop(self):
        """Default check_dependencies is a no-op (no ImportError)."""

        class NoDepTest(AssociationTest):
            @property
            def name(self) -> str:
                return "nodep"

            def run(self, gene, contingency_data, config):  # type: ignore[override]
                pass

        test = NoDepTest()
        # Should not raise
        test.check_dependencies()
