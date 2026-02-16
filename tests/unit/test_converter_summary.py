"""Unit tests for expanded summary.json generation in converter.py."""

import json

import pandas as pd
import pytest

from variantcentrifuge.converter import produce_report_json


def _create_test_tsv(tmp_path, data_dict):
    """Create a minimal TSV file from a dict of column->values."""
    df = pd.DataFrame(data_dict)
    tsv_path = tmp_path / "test_variants.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return str(tsv_path)


@pytest.mark.unit
def test_summary_has_inheritance_distribution(tmp_path):
    """Test that summary.json contains inheritance_distribution with correct counts."""
    # Create TSV with inheritance patterns
    tsv_path = _create_test_tsv(
        tmp_path,
        {
            "CHROM": ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6"],
            "POS": [100, 200, 300, 400, 500, 600],
            "REF": ["A", "C", "G", "T", "A", "C"],
            "ALT": ["G", "T", "A", "C", "T", "G"],
            "GENE": ["GENE1", "GENE2", "GENE3", "GENE1", "GENE2", "GENE3"],
            "IMPACT": ["HIGH", "MODERATE", "LOW", "HIGH", "MODERATE", "LOW"],
            "Inheritance_Pattern": [
                "de_novo",
                "autosomal_dominant",
                "autosomal_recessive",
                "de_novo",
                "none",  # Should be excluded
                "reference",  # Should be excluded
            ],
        },
    )

    # Call produce_report_json
    produce_report_json(tsv_path, str(tmp_path))

    # Load the resulting summary.json
    summary_path = tmp_path / "report" / "summary.json"
    assert summary_path.exists(), "summary.json was not created"

    with open(summary_path, encoding="utf-8") as f:
        summary = json.load(f)

    # Assert inheritance_distribution exists
    assert "inheritance_distribution" in summary

    # Assert correct counts (excluding "none" and "reference")
    expected = {
        "de_novo": 2,
        "autosomal_dominant": 1,
        "autosomal_recessive": 1,
    }
    assert summary["inheritance_distribution"] == expected


@pytest.mark.unit
def test_summary_has_top_genes(tmp_path):
    """Test that summary.json contains top_genes sorted by count."""
    # Create TSV with repeated gene names
    tsv_path = _create_test_tsv(
        tmp_path,
        {
            "CHROM": ["chr1"] * 15,
            "POS": list(range(100, 1600, 100)),
            "REF": ["A"] * 15,
            "ALT": ["G"] * 15,
            "GENE": [
                "GENE1",
                "GENE1",
                "GENE1",
                "GENE1",
                "GENE1",  # 5 occurrences
                "GENE2",
                "GENE2",
                "GENE2",  # 3 occurrences
                "GENE3",
                "GENE3",  # 2 occurrences
                "GENE4",  # 1 occurrence
                "GENE5",
                "GENE5",
                "GENE5",
                "GENE5",  # 4 occurrences
            ],
            "IMPACT": ["HIGH"] * 15,
        },
    )

    # Call produce_report_json
    produce_report_json(tsv_path, str(tmp_path))

    # Load the resulting summary.json
    summary_path = tmp_path / "report" / "summary.json"
    with open(summary_path, encoding="utf-8") as f:
        summary = json.load(f)

    # Assert top_genes exists
    assert "top_genes" in summary
    assert isinstance(summary["top_genes"], list)

    # Assert genes are sorted by count descending
    top_genes = summary["top_genes"]
    assert len(top_genes) > 0
    assert all("gene" in item and "count" in item for item in top_genes)

    # Expected order: GENE1(5), GENE5(4), GENE2(3), GENE3(2), GENE4(1)
    assert top_genes[0] == {"gene": "GENE1", "count": 5}
    assert top_genes[1] == {"gene": "GENE5", "count": 4}
    assert top_genes[2] == {"gene": "GENE2", "count": 3}
    assert top_genes[3] == {"gene": "GENE3", "count": 2}
    assert top_genes[4] == {"gene": "GENE4", "count": 1}


@pytest.mark.unit
def test_summary_top_genes_max_10(tmp_path):
    """Test that top_genes contains at most 10 entries."""
    # Create TSV with 15 unique genes
    genes = [f"GENE{i}" for i in range(1, 16)]
    tsv_path = _create_test_tsv(
        tmp_path,
        {
            "CHROM": ["chr1"] * 15,
            "POS": list(range(100, 1600, 100)),
            "REF": ["A"] * 15,
            "ALT": ["G"] * 15,
            "GENE": genes,
            "IMPACT": ["HIGH"] * 15,
        },
    )

    # Call produce_report_json
    produce_report_json(tsv_path, str(tmp_path))

    # Load the resulting summary.json
    summary_path = tmp_path / "report" / "summary.json"
    with open(summary_path, encoding="utf-8") as f:
        summary = json.load(f)

    # Assert max 10 entries
    assert "top_genes" in summary
    assert len(summary["top_genes"]) <= 10


@pytest.mark.unit
def test_summary_has_num_samples(tmp_path):
    """Test that summary.json contains num_samples from GT column."""
    # Create TSV with GT column containing multiple samples
    tsv_path = _create_test_tsv(
        tmp_path,
        {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [100, 200, 300],
            "REF": ["A", "C", "G"],
            "ALT": ["G", "T", "A"],
            "GENE": ["GENE1", "GENE2", "GENE3"],
            "IMPACT": ["HIGH", "MODERATE", "LOW"],
            "GT": [
                "Sample1(0/1);Sample2(0/0);Sample3(1/1)",
                "Sample1(0/1);Sample4(1/1)",
                "Sample2(0/1);Sample3(0/0);Sample5(1/1)",
            ],
        },
    )

    # Call produce_report_json
    produce_report_json(tsv_path, str(tmp_path))

    # Load the resulting summary.json
    summary_path = tmp_path / "report" / "summary.json"
    with open(summary_path, encoding="utf-8") as f:
        summary = json.load(f)

    # Assert num_samples equals unique sample count
    # Unique samples: Sample1, Sample2, Sample3, Sample4, Sample5 = 5
    assert "num_samples" in summary
    assert summary["num_samples"] == 5


@pytest.mark.unit
def test_summary_empty_dataframe(tmp_path):
    """Test that summary.json handles empty DataFrame correctly."""
    # Create empty TSV (header only)
    tsv_path = _create_test_tsv(
        tmp_path,
        {
            "CHROM": [],
            "POS": [],
            "REF": [],
            "ALT": [],
            "GENE": [],
            "IMPACT": [],
            "Inheritance_Pattern": [],
            "GT": [],
        },
    )

    # Call produce_report_json
    produce_report_json(tsv_path, str(tmp_path))

    # Load the resulting summary.json
    summary_path = tmp_path / "report" / "summary.json"
    with open(summary_path, encoding="utf-8") as f:
        summary = json.load(f)

    # Assert correct empty values
    assert summary["inheritance_distribution"] == {}
    assert summary["top_genes"] == []
    assert summary["num_samples"] == 0
    assert summary["num_variants"] == 0
    assert summary["num_genes"] == 0


@pytest.mark.unit
def test_summary_no_inheritance_column(tmp_path):
    """Test that summary.json handles missing Inheritance_Pattern column."""
    # Create TSV without Inheritance_Pattern column
    tsv_path = _create_test_tsv(
        tmp_path,
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "GENE": ["GENE1", "GENE2"],
            "IMPACT": ["HIGH", "MODERATE"],
        },
    )

    # Call produce_report_json
    produce_report_json(tsv_path, str(tmp_path))

    # Load the resulting summary.json
    summary_path = tmp_path / "report" / "summary.json"
    with open(summary_path, encoding="utf-8") as f:
        summary = json.load(f)

    # Assert inheritance_distribution is empty
    assert "inheritance_distribution" in summary
    assert summary["inheritance_distribution"] == {}
