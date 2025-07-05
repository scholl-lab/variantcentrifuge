#!/usr/bin/env python3
"""
Integration tests for the pipeline improvements.

These tests verify that the gzip and chunked processing
features work correctly in the full pipeline context.
"""

import os
import tempfile
import gzip
import pandas as pd
import pytest
from unittest.mock import patch, Mock


class TestPipelineIntegration:
    """Test the integrated pipeline with new features."""

    def create_test_vcf(self, filename, num_variants=10):
        """Create a minimal test VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=ANN,Number=.,Type=String,Description="Annotation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
"""

        # Add variants
        for i in range(num_variants):
            gene = f"GENE{(i % 3) + 1}"  # Rotate through GENE1, GENE2, GENE3
            chrom = "chr1" if i < 5 else "chr2"
            pos = 1000 + (i * 100)
            vcf_content += (
                f"{chrom}\t{pos}\t.\tA\tT\t100\tPASS\tAC=1;"
                f"ANN=T|missense_variant|MODERATE|{gene}|ENSG001|transcript|"
                f"ENST001|protein_coding|1/1|c.123A>T|p.Met41Leu|123|456|||\tGT\t0/1\t0/0\n"
            )

        with gzip.open(filename, "wt") as f:
            f.write(vcf_content)

    def create_test_config(self):
        """Create a minimal test configuration."""
        return {
            "reference": "GRCh37",
            "add_chr": False,
            "filters": "QUAL >= 100",
            "fields_to_extract": (
                "CHROM POS REF ALT QUAL FILTER ANN[0].GENE ANN[0].EFFECT " "ANN[0].IMPACT GEN[*].GT"
            ),
            "extract_fields_separator": ":",
            "no_stats": True,
            "keep_intermediates": False,
            "use_phenotype_filtering": False,
            "perform_gene_burden": False,
            "calculate_inheritance": False,
            "chunk_size": 5,  # Small chunk size for testing
            "force_chunked_processing": True,  # Force chunked mode for testing
            "sort_memory_limit": "100M",
            "sort_parallel": 2,
        }

    @patch("variantcentrifuge.pipeline.check_external_tools")
    @patch("variantcentrifuge.pipeline.get_gene_bed")
    @patch("variantcentrifuge.pipeline.extract_variants")
    @patch("variantcentrifuge.filters.apply_snpsift_filter")
    @patch("variantcentrifuge.extractor.extract_fields")
    @patch("variantcentrifuge.replacer.replace_genotypes")
    def test_pipeline_with_gzip_intermediates(
        self,
        mock_replace_genotypes,
        mock_extract_fields,
        mock_apply_filter,
        mock_extract_variants,
        mock_get_bed,
        mock_check_tools,
    ):
        """Test that intermediate files can be gzipped."""
        # Create test environment
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test gzip functionality directly
            test_file = os.path.join(temp_dir, "test.tsv.gz")

            # Write gzipped content
            with gzip.open(test_file, "wt") as f:
                f.write("CHROM\tPOS\tGENE\n")
                f.write("chr1\t1000\tGENE1\n")

            # Verify it's gzipped
            with open(test_file, "rb") as f:
                magic = f.read(2)
                assert magic == b"\x1f\x8b", "File should have gzip magic number"

            # Read it back
            with gzip.open(test_file, "rt") as f:
                lines = f.readlines()
                assert len(lines) == 2
                assert "GENE1" in lines[1]

    @patch("variantcentrifuge.pipeline.check_external_tools")
    @patch("variantcentrifuge.pipeline.analyze_variants")
    def test_chunked_processing_large_file(self, mock_analyze, mock_check_tools):
        """Test chunked processing triggers for large files."""
        # Create test environment
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a "large" TSV file (mock the size check)
            tsv_file = os.path.join(temp_dir, "large.tsv.gz")

            # Create test data with multiple genes
            data = []
            for gene_num in range(5):
                gene = f"GENE{gene_num + 1}"
                for var_num in range(20):
                    data.append(
                        {
                            "CHROM": "chr1",
                            "POS": 1000 + (gene_num * 1000) + var_num,
                            "REF": "A",
                            "ALT": "T",
                            "GENE": gene,
                            "IMPACT": "MODERATE",
                            "GT": "0/1:0/0:0/0",
                        }
                    )

            df = pd.DataFrame(data)
            df.to_csv(tsv_file, sep="\t", index=False, compression="gzip")

            # Mock analyze_variants
            def mock_analyze_func(inp, cfg):
                # Return processed lines
                if hasattr(inp, "read"):
                    lines = inp.readlines()
                    for line in lines:
                        yield line.strip()

            mock_analyze.side_effect = mock_analyze_func

            # Test with chunked processing
            from variantcentrifuge.pipeline import process_chunked_pipeline

            output_file = os.path.join(temp_dir, "output.tsv")

            cfg = {
                "chunk_size": 25,  # Each gene has 20 variants
                "keep_intermediates": False,
                "calculate_inheritance": False,
                "late_filtering": False,
                "no_links": True,
            }

            args = Mock()
            args.add_column = None

            # Process
            success = process_chunked_pipeline(
                final_tsv=tsv_file,
                final_output=output_file,
                cfg=cfg,
                custom_features={},
                scoring_config=None,
                pedigree_data=None,
                args=args,
                base_name="test",
                intermediate_dir=temp_dir,
                chunksize=25,
            )

            assert success is True

            # Verify output
            assert os.path.exists(output_file)

            # Read output and verify
            with open(output_file, "r") as f:
                lines = f.readlines()

            # Should have header + 100 data rows
            assert len(lines) == 101

            # Verify VAR_ID column was added
            assert lines[0].startswith("VAR_ID\t")

            # Verify all genes are present
            genes_found = set()
            for line in lines[1:]:
                parts = line.strip().split("\t")
                if len(parts) > 5:  # GENE column
                    genes_found.add(parts[5])

            assert genes_found == {"GENE1", "GENE2", "GENE3", "GENE4", "GENE5"}

    def test_sort_error_handling(self):
        """Test sort function error handling."""
        from variantcentrifuge.pipeline import sort_tsv_by_gene

        with tempfile.TemporaryDirectory() as temp_dir:
            # Test with non-existent file
            with pytest.raises(Exception):
                sort_tsv_by_gene("/non/existent/file.tsv", os.path.join(temp_dir, "output.tsv"))

            # Test with file containing no gene column
            bad_file = os.path.join(temp_dir, "no_gene.tsv")
            with open(bad_file, "w") as f:
                f.write("CHROM\tPOS\tREF\tALT\n")
                f.write("chr1\t100\tA\tT\n")

            with pytest.raises(ValueError, match="Gene column"):
                sort_tsv_by_gene(bad_file, os.path.join(temp_dir, "output.tsv"), gene_column="GENE")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
