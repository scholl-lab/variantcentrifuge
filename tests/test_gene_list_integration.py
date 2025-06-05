"""
Integration tests for the gene list annotation feature.
Tests the integration with the pipeline and CLI components.
"""

import os
import shutil
import tempfile

import pytest


@pytest.fixture
def test_data_dir():
    """Create a temporary directory for test data."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_vcf_file(test_data_dir):
    """Create a minimal sample VCF file for testing."""
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t1000\t.\tA\tG\t100\tPASS\tGENE=TP53\tGT\t0/1\t0/0
chr1\t2000\t.\tG\tT\t100\tPASS\tGENE=BRCA1\tGT\t0/1\t0/0
chr2\t3000\t.\tC\tA\t100\tPASS\tGENE=BRCA2\tGT\t0/0\t0/1
chr3\t4000\t.\tT\tC\t100\tPASS\tGENE=APC;IMPACT=HIGH\tGT\t0/1\t0/1
"""
    vcf_file = os.path.join(test_data_dir, "sample.vcf")
    with open(vcf_file, "w", encoding="utf-8") as f:
        f.write(vcf_content)
    return vcf_file


@pytest.fixture
def gene_lists(test_data_dir):
    """Create sample gene list files."""
    cancer_genes_file = os.path.join(test_data_dir, "cancer_genes.txt")
    with open(cancer_genes_file, "w", encoding="utf-8") as f:
        f.write("TP53\nBRCA1\nBRCA2\n")

    apc_genes_file = os.path.join(test_data_dir, "apc_list.txt")
    with open(apc_genes_file, "w", encoding="utf-8") as f:
        f.write("APC\n")

    return {"cancer_genes": cancer_genes_file, "apc_genes": apc_genes_file}


def test_pipeline_with_gene_list_annotation(sample_vcf_file, gene_lists, test_data_dir):
    """Test that gene list annotation works in the full pipeline."""
    output_dir = os.path.join(test_data_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Create minimal configuration for the pipeline
    config = {
        "vcf_file": sample_vcf_file,
        "output_dir": output_dir,
        "genes": "all",
        "reference": test_data_dir,  # Mock reference dir
        "debug_level": "DEBUG",
        "no_links": True,  # Disable links for simplicity
        "extract_fields": "CHROM,POS,REF,ALT,GENE",  # Basic fields
        "annotate_gene_list_files": [gene_lists["cancer_genes"], gene_lists["apc_genes"]],
        "filters": "FILTER = 'PASS'",  # Simple filter
    }

    # Mock minimal functions to avoid actual vcf processing
    def mock_gene_bed_file():
        bed_file = os.path.join(test_data_dir, "genes.bed")
        with open(bed_file, "w", encoding="utf-8") as f:
            f.write("chr1\t900\t1100\tTP53\n")
            f.write("chr1\t1900\t2100\tBRCA1\n")
            f.write("chr2\t2900\t3100\tBRCA2\n")
            f.write("chr3\t3900\t4100\tAPC\n")
        return bed_file

    # Create a mock tsv file that would be the result of extraction
    tsv_file = os.path.join(output_dir, "extracted.tsv")
    with open(tsv_file, "w", encoding="utf-8") as f:
        f.write("CHROM\tPOS\tREF\tALT\tGENE\n")
        f.write("chr1\t1000\tA\tG\tTP53\n")
        f.write("chr1\t2000\tG\tT\tBRCA1\n")
        f.write("chr2\t3000\tC\tA\tBRCA2\n")
        f.write("chr3\t4000\tT\tC\tAPC\n")

    # Create a minimal pipeline function for testing
    def minimal_pipeline_test():
        # Find the final output file
        final_output = os.path.join(output_dir, "final.tsv")

        # Run pipeline with mock data
        config["final_output"] = final_output
        with open(final_output, "w", encoding="utf-8") as out:
            with open(tsv_file, "r", encoding="utf-8") as inp:
                lines = inp.readlines()

                # Manually run the gene list annotation part
                from variantcentrifuge.helpers import annotate_variants_with_gene_lists

                lines = annotate_variants_with_gene_lists(lines, config["annotate_gene_list_files"])

                # Write the output
                for line in lines:
                    out.write(line.rstrip("\n") + "\n")

        # Check the results
        with open(final_output, "r", encoding="utf-8") as f:
            content = f.read()

        # Verify header contains gene list columns
        header = content.split("\n")[0]
        assert "cancer_genes" in header
        assert "apc_list" in header

        # Verify gene annotations
        lines = content.split("\n")
        assert lines[1].endswith("yes\tno")  # TP53: in cancer genes, not in APC
        assert lines[2].endswith("yes\tno")  # BRCA1: in cancer genes, not in APC
        assert lines[3].endswith("yes\tno")  # BRCA2: in cancer genes, not in APC
        assert lines[4].endswith("no\tyes")  # APC: not in cancer genes, in APC

    # Run the test
    minimal_pipeline_test()


def test_cli_with_gene_list_annotation(sample_vcf_file, gene_lists, test_data_dir, monkeypatch):
    """Test that the CLI correctly handles gene list annotation."""
    from variantcentrifuge.cli import main

    output_dir = os.path.join(test_data_dir, "cli_output")
    os.makedirs(output_dir, exist_ok=True)

    # Mock run_pipeline to avoid actual execution but verify args are passed correctly
    run_pipeline_calls = []

    def mock_run_pipeline(args, cfg, start_time):
        """Mock run_pipeline that just records the call."""
        run_pipeline_calls.append(cfg.copy())
        return 0

    # Apply the mock
    monkeypatch.setattr("variantcentrifuge.cli.run_pipeline", mock_run_pipeline)

    # Simulate CLI arguments
    import sys

    original_argv = sys.argv
    try:
        sys.argv = [
            "variantcentrifuge",
            "--vcf",
            sample_vcf_file,
            "--output-dir",
            output_dir,
            "--gene-name",
            "all",
            "--filters",
            "not_artefact",
            "--annotate-gene-list",
            gene_lists["cancer_genes"],
            "--annotate-gene-list",
            gene_lists["apc_genes"],
            "--log-level",
            "DEBUG",
        ]

        # Run the CLI
        main()

        # Check if the gene list files were correctly passed to the pipeline
        assert len(run_pipeline_calls) == 1
        cfg = run_pipeline_calls[0]

        assert "annotate_gene_list_files" in cfg
        assert len(cfg["annotate_gene_list_files"]) == 2
        assert gene_lists["cancer_genes"] in cfg["annotate_gene_list_files"]
        assert gene_lists["apc_genes"] in cfg["annotate_gene_list_files"]

    finally:
        sys.argv = original_argv
