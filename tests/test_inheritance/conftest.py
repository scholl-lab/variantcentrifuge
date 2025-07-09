"""Pytest configuration and fixtures for inheritance tests."""

from typing import Any, Dict

import pandas as pd
import pytest


@pytest.fixture
def sample_pedigree_data() -> Dict[str, Dict[str, Any]]:
    """Provide sample pedigree data for testing."""
    return {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",  # Male
            "affected_status": "2",  # Affected
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",  # Male
            "affected_status": "1",  # Unaffected
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",  # Female
            "affected_status": "1",  # Unaffected
        },
    }


@pytest.fixture
def single_sample_pedigree() -> Dict[str, Dict[str, Any]]:
    """Provide single sample pedigree data."""
    return {
        "Sample1": {
            "sample_id": "Sample1",
            "father_id": "0",
            "mother_id": "0",
            "sex": "0",  # Unknown
            "affected_status": "2",  # Affected
        }
    }


@pytest.fixture
def multi_generation_pedigree() -> Dict[str, Dict[str, Any]]:
    """Provide multi-generation pedigree data."""
    return {
        "grandparent": {
            "sample_id": "grandparent",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "2",
        },
        "parent": {
            "sample_id": "parent",
            "father_id": "grandparent",
            "mother_id": "grandparent_spouse",
            "sex": "2",
            "affected_status": "2",
        },
        "grandparent_spouse": {
            "sample_id": "grandparent_spouse",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
        "child": {
            "sample_id": "child",
            "father_id": "spouse",
            "mother_id": "parent",
            "sex": "1",
            "affected_status": "2",
        },
        "spouse": {
            "sample_id": "spouse",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
    }


@pytest.fixture
def sample_variant_df() -> pd.DataFrame:
    """Provide sample variant DataFrame."""
    return pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "GENE1",
                "IMPACT": "HIGH",
                "AF": "0.001",
                "Sample1": "0/1",
                "Sample2": "0/0",
                "Sample3": "1/1",
            },
            {
                "CHROM": "2",
                "POS": "2000",
                "REF": "C",
                "ALT": "G",
                "GENE": "GENE2",
                "IMPACT": "MODERATE",
                "AF": "0.01",
                "Sample1": "0/0",
                "Sample2": "0/1",
                "Sample3": "0/1",
            },
        ]
    )


@pytest.fixture
def compound_het_df() -> pd.DataFrame:
    """Provide DataFrame with potential compound heterozygous variants."""
    return pd.DataFrame(
        [
            # GENE1 - multiple heterozygous variants
            {
                "CHROM": "1",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "GENE1",
                "Sample1": "0/1",
                "Sample2": "0/0",
            },
            {
                "CHROM": "1",
                "POS": "2000",
                "REF": "C",
                "ALT": "G",
                "GENE": "GENE1",
                "Sample1": "0/1",
                "Sample2": "0/1",
            },
            {
                "CHROM": "1",
                "POS": "3000",
                "REF": "G",
                "ALT": "A",
                "GENE": "GENE1",
                "Sample1": "0/1",
                "Sample2": "0/0",
            },
            # GENE2 - single variant
            {
                "CHROM": "2",
                "POS": "4000",
                "REF": "T",
                "ALT": "C",
                "GENE": "GENE2",
                "Sample1": "0/1",
                "Sample2": "0/1",
            },
        ]
    )


@pytest.fixture
def x_linked_variants() -> pd.DataFrame:
    """Provide X-linked variant data."""
    return pd.DataFrame(
        [
            {
                "CHROM": "X",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "DMD",
                "affected_male": "0/1",  # Hemizygous
                "carrier_female": "0/1",  # Heterozygous
                "unaffected_male": "0/0",
                "unaffected_female": "0/0",
            }
        ]
    )


@pytest.fixture
def mitochondrial_variants() -> pd.DataFrame:
    """Provide mitochondrial variant data."""
    return pd.DataFrame(
        [
            {
                "CHROM": "MT",
                "POS": "1000",
                "REF": "A",
                "ALT": "T",
                "GENE": "MT-ND1",
                "affected_mother": "1/1",
                "affected_child1": "1/1",
                "affected_child2": "1/1",
                "unaffected_father": "0/0",
            }
        ]
    )
