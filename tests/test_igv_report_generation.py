#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Integration tests for IGV report generation in variantcentrifuge."""

import json
import os
import shutil
import tempfile

import pytest

from variantcentrifuge.generate_igv_report import generate_igv_report
from variantcentrifuge.utils import generate_igv_safe_filename_base


@pytest.fixture
def create_test_data():
    """Create temporary test files for IGV report generation."""
    # Create temporary directory
    temp_dir = tempfile.mkdtemp()
    report_dir = os.path.join(temp_dir, "report")
    os.makedirs(report_dir, exist_ok=True)

    # Create a simple TSV with a variant with long alleles
    variants_tsv = os.path.join(temp_dir, "variants.tsv")
    with open(variants_tsv, "w", encoding="utf-8") as f:
        f.write("CHROM\tPOS\tREF\tALT\tGT\n")
        f.write("chr1\t12345\tACGTACGT\tTGCATGCA\tSAMPLE1(1/1);SAMPLE2(0/1)\n")

    # Create a BAM mapping file
    bam_mapping = os.path.join(temp_dir, "bam_mapping.tsv")
    with open(bam_mapping, "w", encoding="utf-8") as f:
        f.write("sample_id,bam_path\n")
        f.write("SAMPLE1,/path/to/sample1.bam\n")
        f.write("SAMPLE2,/path/to/sample2.bam\n")

    yield {
        "temp_dir": temp_dir,
        "report_dir": report_dir,
        "variants_tsv": variants_tsv,
        "bam_mapping": bam_mapping,
    }

    # Cleanup
    shutil.rmtree(temp_dir)


def test_igv_report_generation_with_long_alleles(monkeypatch, create_test_data):
    """Test IGV report generation with long alleles that would cause filename issues."""
    test_data = create_test_data

    # Mock subprocess.run to avoid actually running igv-reports
    def mock_run(cmd, *args, **kwargs):
        # Just create the report file that would have been created
        output_path = None
        for i, arg in enumerate(cmd):
            if arg == "--output" and i + 1 < len(cmd):
                output_path = cmd[i + 1]

        if output_path:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write("<html><body>Mock IGV Report</body></html>")

        # Fixed MockCompletedProcess with stdout and stderr attributes
        class MockCompletedProcess:
            def __init__(self):
                self.returncode = 0
                self.stdout = ""
                self.stderr = ""

        return MockCompletedProcess()

    monkeypatch.setattr("subprocess.run", mock_run)

    # Run the report generation function
    generate_igv_report(
        variants_tsv=test_data["variants_tsv"],
        output_dir=test_data["report_dir"],
        bam_mapping_file=test_data["bam_mapping"],
        igv_reference="hg38",  # Use a reference genome name
        igv_max_allele_len_filename=8,  # Use smaller values to force truncation
        igv_hash_len_filename=5,
        igv_max_variant_part_filename=40,
    )

    # Check that the igv_reports directory was created
    igv_dir = os.path.join(test_data["report_dir"], "igv")
    assert os.path.exists(igv_dir)

    # Check that the mapping file was created
    mapping_file = os.path.join(igv_dir, "igv_reports_map.json")
    assert os.path.exists(mapping_file)

    # Read the mapping file to verify structure
    with open(mapping_file, "r", encoding="utf-8") as f:
        mapping = json.load(f)

    # Verify that the mapping contains the expected keys
    assert "variants" in mapping
    assert len(mapping["variants"]) > 0

    # Check that the first variant has the expected structure
    first_variant = mapping["variants"][0]
    assert "chrom" in first_variant
    assert "pos" in first_variant
    assert "ref" in first_variant
    assert "alt" in first_variant
    assert "sample_reports" in first_variant

    # Verify that the full original alleles are preserved in the mapping
    assert first_variant["ref"] == "ACGTACGT"
    assert first_variant["alt"] == "TGCATGCA"

    # Check that at least one report file was generated (for SAMPLE2 which has 0/1 genotype)
    assert len(first_variant["sample_reports"]) > 0

    # For each sample with a report, verify the file exists and uses the shortened name
    for sample_id, report_path in first_variant["sample_reports"].items():
        # The report path is relative to the report directory
        full_path = os.path.join(test_data["report_dir"], report_path)
        assert os.path.exists(full_path)

        # Extract the filename part
        filename = os.path.basename(full_path)

        # Verify that the filename contains the truncated REF/ALT alleles
        # The safe filename should match our utility function result
        safe_base = generate_igv_safe_filename_base(
            sample_id=sample_id.upper(),  # Sample IDs are uppercase in BAM mapping
            chrom="chr1",
            pos="12345",
            ref="ACGTACGT",  # Changed to match the actual input values
            alt="TGCATGCA",  # Changed to match the actual input values
            max_allele_len=8,
            hash_len=5,
            max_variant_part_len=40,
        )
        # Verify the filename uses the safe base
        expected_filename = f"{safe_base}_igv_report.html"
        assert filename == expected_filename

        # Ensure that the full allele lengths are not in the filename
        assert len(filename) < 100, "Filename is still too long"
        assert "_" in filename, "Special characters not properly handled in filename"


def test_igv_report_with_local_fasta(monkeypatch, create_test_data):
    """Test IGV report generation with local FASTA file instead of reference genome."""
    test_data = create_test_data

    # Create a mock local FASTA file and index
    fasta_file = os.path.join(test_data["temp_dir"], "reference.fa")
    with open(fasta_file, "w") as f:
        f.write(">chr1\nACGTACGT\n>chr2\nTGCATGCA\n")

    # Create a mock FASTA index file with standard naming convention (FASTA filename + .fai)
    # This follows the convention where the index exists alongside the FASTA with .fai extension
    with open(f"{fasta_file}.fai", "w") as f:
        f.write("chr1\t8\t6\t8\t9\n")
        f.write("chr2\t8\t21\t8\t9\n")

    # Mock subprocess.run to avoid actually running igv-reports
    def mock_run(cmd, *args, **kwargs):
        # Extract the output path
        output_path = None
        for i, arg in enumerate(cmd):
            if arg == "--output" and i + 1 < len(cmd):
                output_path = cmd[i + 1]

        if output_path:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write("<html><body>Mock IGV Report With Local FASTA</body></html>")

        # Fixed MockCompletedProcess with stdout and stderr attributes
        class MockCompletedProcess:
            def __init__(self):
                self.returncode = 0
                self.stdout = ""
                self.stderr = ""

        return MockCompletedProcess()

    monkeypatch.setattr("subprocess.run", mock_run)

    # Run the report generation function with local FASTA
    generate_igv_report(
        variants_tsv=test_data["variants_tsv"],
        output_dir=test_data["report_dir"],
        bam_mapping_file=test_data["bam_mapping"],
        igv_fasta=fasta_file,  # Use local FASTA instead of reference genome
        # FASTA index is expected to be at fasta_file + '.fai' following standard convention
        igv_max_allele_len_filename=8,
        igv_hash_len_filename=5,
        igv_max_variant_part_filename=40,
    )

    # Check that the igv_reports directory was created
    igv_dir = os.path.join(test_data["report_dir"], "igv")
    assert os.path.exists(igv_dir)

    # Check that the mapping file was created
    mapping_file = os.path.join(igv_dir, "igv_reports_map.json")
    assert os.path.exists(mapping_file)

    # Verify some report files were created
    report_files = [f for f in os.listdir(igv_dir) if f.endswith("_igv_report.html")]
    assert len(report_files) > 0
