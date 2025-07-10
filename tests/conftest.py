"""Shared pytest fixtures for all test modules."""

from pathlib import Path
from typing import Any, Dict

import pandas as pd
import pytest


@pytest.fixture
def gene_burden_test_config() -> Dict[str, Any]:
    """Standard gene burden analysis configuration."""
    return {
        "gene_burden_mode": "samples",
        "correction_method": "fdr",
        "confidence_interval_method": "normal_approx",
        "confidence_interval_alpha": 0.05,
        "continuity_correction": 0.5,
    }


@pytest.fixture
def gene_burden_alleles_config() -> Dict[str, Any]:
    """Gene burden configuration for alleles mode."""
    return {
        "gene_burden_mode": "alleles",
        "correction_method": "fdr",
        "confidence_interval_method": "normal_approx",
        "confidence_interval_alpha": 0.05,
    }


@pytest.fixture
def gene_burden_edge_case_data() -> pd.DataFrame:
    """Standard edge case data for gene burden testing."""
    data = [
        # Gene with zero variants in controls (infinite OR)
        {
            "GENE": "GENE1",
            "proband_count": 100,
            "control_count": 100,
            "proband_variant_count": 10,
            "control_variant_count": 0,
            "proband_allele_count": 15,
            "control_allele_count": 0,
        },
        # Gene with zero variants in cases (zero OR)
        {
            "GENE": "GENE2",
            "proband_count": 100,
            "control_count": 100,
            "proband_variant_count": 0,
            "control_variant_count": 10,
            "proband_allele_count": 0,
            "control_allele_count": 15,
        },
        # Gene with zero variants in both groups
        {
            "GENE": "GENE3",
            "proband_count": 100,
            "control_count": 100,
            "proband_variant_count": 0,
            "control_variant_count": 0,
            "proband_allele_count": 0,
            "control_allele_count": 0,
        },
        # Gene with normal distribution
        {
            "GENE": "GENE4",
            "proband_count": 100,
            "control_count": 100,
            "proband_variant_count": 10,
            "control_variant_count": 5,
            "proband_allele_count": 12,
            "control_allele_count": 6,
        },
        # Gene with single variant
        {
            "GENE": "GENE5",
            "proband_count": 100,
            "control_count": 100,
            "proband_variant_count": 1,
            "control_variant_count": 1,
            "proband_allele_count": 1,
            "control_allele_count": 1,
        },
    ]
    return pd.DataFrame(data)


@pytest.fixture
def gene_burden_simple_test_data() -> pd.DataFrame:
    """Simple gene burden test data for unit testing."""
    test_data = {
        "GENE": ["GENE1", "GENE1", "GENE2", "GENE2"],
        "proband_count": [10, 10, 10, 10],
        "control_count": [10, 10, 10, 10],
        "proband_variant_count": [2, 1, 3, 2],
        "control_variant_count": [1, 0, 1, 1],
        "proband_allele_count": [3, 1, 4, 3],
        "control_allele_count": [1, 0, 2, 1],
    }
    return pd.DataFrame(test_data)


@pytest.fixture
def gene_burden_test_data_paths():
    """Standard test data file paths for gene burden integration tests."""
    fixtures_dir = Path(__file__).parent / "fixtures" / "geneburden" / "output"
    return {
        "vcf": fixtures_dir / "enhanced_test_data.vcf.gz",
        "case_samples": fixtures_dir / "case_samples.txt",
        "control_samples": fixtures_dir / "control_samples.txt",
        "genes": fixtures_dir / "test_genes.txt",
    }


@pytest.fixture
def gene_burden_comprehensive_test_data_paths():
    """Test data paths for comprehensive gene burden testing."""
    fixtures_dir = Path(__file__).parent / "fixtures" / "geneburden" / "comprehensive_test_data"
    return {
        "vcf": fixtures_dir / "test_data.vcf.gz",
        "genes": fixtures_dir / "test_genes.txt",
        "case_samples": fixtures_dir / "case_samples.txt",
        "control_samples": fixtures_dir / "control_samples.txt",
        "phenotypes": fixtures_dir / "phenotypes.tsv",
    }


@pytest.fixture
def gene_burden_cmd_template():
    """Template for creating gene burden analysis commands."""

    def _create_cmd(
        vcf_file,
        gene_file,
        case_samples=None,
        control_samples=None,
        case_samples_file=None,
        control_samples_file=None,
        output_dir="output",
        output_file="results.tsv",
        preset="rare,coding",
        **kwargs,
    ):
        """Create a command list for gene burden analysis."""
        cmd = [
            "variantcentrifuge",
            "--vcf-file",
            str(vcf_file),
            "--gene-file",
            str(gene_file),
            "--perform-gene-burden",
            "--preset",
            preset,
            "--output-dir",
            str(output_dir),
            "--output-file",
            output_file,
            "--use-new-pipeline",
        ]

        if case_samples_file and control_samples_file:
            cmd.extend(["--case-samples-file", str(case_samples_file)])
            cmd.extend(["--control-samples-file", str(control_samples_file)])
        elif case_samples and control_samples:
            cmd.extend(["--case-samples", ",".join(case_samples)])
            cmd.extend(["--control-samples", ",".join(control_samples)])

        for key, value in kwargs.items():
            cmd.extend([f"--{key.replace('_', '-')}", str(value)])

        return cmd

    return _create_cmd


@pytest.fixture
def gene_burden_required_columns():
    """Required columns for gene burden results validation."""
    return [
        "GENE",
        "proband_count",
        "control_count",
        "proband_allele_count",
        "control_allele_count",
        "raw_p_value",
        "corrected_p_value",
        "odds_ratio",
    ]


@pytest.fixture
def case_control_gt_sample_data():
    """Sample data with GT column for case/control testing."""
    return {
        "GENE": ["GENE1", "GENE1", "GENE2"],
        "GT": [
            "CASE_001(1/1);CASE_002(0/1);CTRL_001(0/0);CTRL_002(0/1)",
            "CASE_001(0/1);CASE_003(0/0);CTRL_001(1/1);CTRL_003(0/1)",
            "CASE_002(1/1);CTRL_001(0/1);CTRL_002(0/0)",
        ],
    }


@pytest.fixture
def case_control_sample_sets():
    """Standard case and control sample sets."""
    return {
        "case_samples": {"CASE_001", "CASE_002", "CASE_003"},
        "control_samples": {"CTRL_001", "CTRL_002", "CTRL_003"},
    }
