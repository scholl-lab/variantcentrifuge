"""Unit tests for helpers module NoneType fixes."""

import pandas as pd

from variantcentrifuge.helpers import determine_case_control_sets


class TestDetermineCaseControlSetsNoneTypeFix:
    """Test the fixes for NoneType errors in determine_case_control_sets function."""

    def test_none_case_samples_handling(self):
        """Test that None case_samples is handled correctly without NoneType errors."""
        all_samples = {"SAMPLE1", "SAMPLE2", "SAMPLE3"}

        # Config with None case_samples (this previously caused NoneType iteration error)
        cfg = {
            "case_samples": None,
            "control_samples": ["SAMPLE1"],
            "case_phenotypes": None,
            "control_phenotypes": None,
            "phenotypes": {},
        }

        df = pd.DataFrame()

        # This should not raise NoneType iteration error
        case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)

        # When case_samples is None but control_samples is provided,
        # case_samples should be derived as all_samples - control_samples
        assert case_samples == {"SAMPLE2", "SAMPLE3"}
        assert control_samples == {"SAMPLE1"}

    def test_none_control_samples_handling(self):
        """Test that None control_samples is handled correctly without NoneType errors."""
        all_samples = {"SAMPLE1", "SAMPLE2", "SAMPLE3"}

        # Config with None control_samples
        cfg = {
            "case_samples": ["SAMPLE1"],
            "control_samples": None,
            "case_phenotypes": None,
            "control_phenotypes": None,
            "phenotypes": {},
        }

        df = pd.DataFrame()

        # This should not raise NoneType iteration error
        case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)

        # When control_samples is None but case_samples is provided,
        # control_samples should be derived as all_samples - case_samples
        assert case_samples == {"SAMPLE1"}
        assert control_samples == {"SAMPLE2", "SAMPLE3"}

    def test_none_case_phenotypes_handling(self):
        """Test that None case_phenotypes is handled correctly without NoneType errors."""
        all_samples = {"SAMPLE1", "SAMPLE2", "SAMPLE3"}

        # Config with None case_phenotypes (this previously caused NoneType iteration error)
        cfg = {
            "case_samples": None,
            "control_samples": None,
            "case_phenotypes": None,  # This was causing 'NoneType' object is not iterable
            "control_phenotypes": ["healthy"],
            "phenotypes": {"SAMPLE1": {"healthy"}, "SAMPLE2": {"affected"}, "SAMPLE3": {"healthy"}},
        }

        df = pd.DataFrame()

        # This should not raise NoneType iteration error
        case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)

        # With control_phenotypes only, samples matching control terms become controls
        # and samples not matching become cases
        assert case_samples == {"SAMPLE2"}  # SAMPLE2 doesn't match "healthy"
        assert control_samples == {"SAMPLE1", "SAMPLE3"}  # These match "healthy"

    def test_none_control_phenotypes_handling(self):
        """Test that None control_phenotypes is handled correctly without NoneType errors."""
        all_samples = {"SAMPLE1", "SAMPLE2", "SAMPLE3"}

        # Config with None control_phenotypes
        cfg = {
            "case_samples": None,
            "control_samples": None,
            "case_phenotypes": ["affected"],
            "control_phenotypes": None,  # This was causing 'NoneType' object is not iterable
            "phenotypes": {"SAMPLE1": {"healthy"}, "SAMPLE2": {"affected"}, "SAMPLE3": {"healthy"}},
        }

        df = pd.DataFrame()

        # This should not raise NoneType iteration error
        case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)

        # With case_phenotypes only, samples matching case terms become cases
        # and samples not matching become controls
        assert case_samples == {"SAMPLE2"}  # SAMPLE2 matches "affected"
        assert control_samples == {"SAMPLE1", "SAMPLE3"}  # These don't match "affected"

    def test_all_none_values_default_behavior(self):
        """Test behavior when all relevant config values are None."""
        all_samples = {"SAMPLE1", "SAMPLE2", "SAMPLE3"}

        # Config with all None values
        cfg = {
            "case_samples": None,
            "control_samples": None,
            "case_phenotypes": None,
            "control_phenotypes": None,
            "phenotypes": {},
        }

        df = pd.DataFrame()

        # This should not raise NoneType iteration error
        case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)

        # When no criteria are provided, all samples should become controls
        assert case_samples == set()
        assert control_samples == all_samples

    def test_mixed_none_and_empty_list_handling(self):
        """Test correct handling of mixed None and empty list values."""
        all_samples = {"SAMPLE1", "SAMPLE2", "SAMPLE3"}

        # Config mixing None and empty lists
        cfg = {
            "case_samples": [],  # Empty list
            "control_samples": None,  # None
            "case_phenotypes": None,  # None
            "control_phenotypes": [],  # Empty list
            "phenotypes": {},
        }

        df = pd.DataFrame()

        # This should not raise NoneType iteration error
        case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)

        # With no valid criteria, all samples should become controls
        assert case_samples == set()
        assert control_samples == all_samples

    def test_phenotype_none_values_in_compute_function(self):
        """Test compute_phenotype_based_case_control_assignment handles None."""
        from variantcentrifuge.helpers import (
            compute_phenotype_based_case_control_assignment,
        )

        vcf_samples = ["SAMPLE1", "SAMPLE2", "SAMPLE3"]
        phenotype_map = {"SAMPLE1": {"healthy"}, "SAMPLE2": {"affected"}, "SAMPLE3": {"healthy"}}

        # Test with None case_phenotypes (this was causing the error)
        case_samples, control_samples = compute_phenotype_based_case_control_assignment(
            vcf_samples=vcf_samples,
            phenotype_map=phenotype_map,
            case_phenotypes=None,  # This was causing 'NoneType' object is not iterable
            control_phenotypes=["healthy"],
        )

        # Should work without NoneType error
        assert isinstance(case_samples, set)
        assert isinstance(control_samples, set)
        assert control_samples == {"SAMPLE1", "SAMPLE3"}  # Match "healthy"
        assert case_samples == {"SAMPLE2"}  # Doesn't match "healthy"

    def test_edge_case_empty_phenotype_map_with_none_values(self):
        """Test edge case with empty phenotype map and None values."""
        all_samples = {"SAMPLE1", "SAMPLE2"}

        cfg = {
            "case_samples": None,
            "control_samples": None,
            "case_phenotypes": None,
            "control_phenotypes": None,
            "phenotypes": {},  # Empty phenotype map
        }

        df = pd.DataFrame()

        # Should not raise any errors
        case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)

        # Default behavior: all samples become controls
        assert case_samples == set()
        assert control_samples == all_samples
