"""
Integration tests for --final-filter functionality.

Tests the complete pipeline with final filtering applied to ensure
filters work correctly on computed columns like scores and inheritance patterns.
"""

import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock
from variantcentrifuge.cli import main


class TestFinalFilterIntegration:
    """Integration tests for final filter feature."""
    
    def test_final_filter_with_score(self, monkeypatch, tmp_path):
        """Test that final filter works with computed score columns."""
        # Create minimal test files
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text("""##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t1000\t.\tA\tT\t100\tPASS\tAF=0.01\tGT\t0/1
chr1\t2000\t.\tG\tC\t100\tPASS\tAF=0.5\tGT\t1/1
""")
        
        output_dir = tmp_path / "output"
        output_file = output_dir / "test.final.tsv"
        
        # Mock external tools
        def mock_check_tools(*args, **kwargs):
            return True
            
        def mock_run_command(cmd, output_file=None):
            # Simulate tool outputs
            if "genes2bed" in cmd:
                bed_file = cmd[-1]
                with open(bed_file, "w") as f:
                    f.write("chr1\t500\t3000\tTEST_GENE\n")
            elif "bcftools" in cmd and "view" in cmd:
                # Simulate variant extraction
                output = cmd[cmd.index("-o") + 1]
                vcf_file.rename(output)
            elif "SnpSift" in cmd and "filter" in cmd:
                # Pass through
                pass
            elif "SnpSift" in cmd and "extractFields" in cmd:
                # Create a simple TSV output
                output = output_file or cmd[-1]
                with open(output, "w") as f:
                    f.write("CHROM\tPOS\tREF\tALT\tGENE\tIMPACT\tScore\n")
                    f.write("chr1\t1000\tA\tT\tTEST_GENE\tHIGH\t0.9\n")
                    f.write("chr1\t2000\tG\tC\tTEST_GENE\tMODERATE\t0.3\n")
            elif "bgzip" in cmd:
                # Mock compression
                pass
            elif "bcftools" in cmd and "index" in cmd:
                # Mock indexing
                pass
            return 0
            
        # Mock subprocess.run for gene_bed module
        def mock_subprocess_run(cmd, **kwargs):
            if "genes2bed" in cmd:
                # Write to the output file specified in stdout
                if 'stdout' in kwargs and hasattr(kwargs['stdout'], 'write'):
                    kwargs['stdout'].write("chr1\t500\t3000\tTEST_GENE\n")
            elif "sort" in cmd:
                # For sort command, just pass through
                pass
            return MagicMock(returncode=0)
        
        monkeypatch.setattr("subprocess.run", mock_subprocess_run)
        monkeypatch.setattr("variantcentrifuge.pipeline.check_external_tools", mock_check_tools)
        monkeypatch.setattr("variantcentrifuge.utils.run_command", mock_run_command)
        monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)
        monkeypatch.setattr("variantcentrifuge.extractor.run_command", mock_run_command)
        
        # Run variantcentrifuge with final filter
        test_args = [
            "variantcentrifuge",
            "--vcf-file", str(vcf_file),
            "--gene-name", "TEST_GENE",
            "--output-dir", str(output_dir),
            "--output-file", "test.final.tsv",
            "--filters", "FILTER='PASS'",
            "--fields", "CHROM POS REF ALT GENE IMPACT Score",
            "--reference", "GRCh38",
            "--no-replacement",
            "--no-stats",
            "--final-filter", "Score > 0.5"  # This is the key test
        ]
        
        with patch("sys.argv", test_args):
            main()
        
        # Check that output was created and filtered correctly
        assert output_file.exists()
        
        # Read the output and verify filtering worked
        with open(output_file) as f:
            lines = f.readlines()
            
        # Should have header + 1 data row (Score > 0.5 filters out the second variant)
        assert len(lines) == 2
        assert "VAR_ID\tCHROM\tPOS\tREF\tALT\tGENE\tIMPACT\tScore" in lines[0]
        assert "chr1\t1000\tA\tT\tTEST_GENE\tHIGH\t0.9" in lines[1]
        assert "chr1\t2000" not in lines[1]  # This variant should be filtered out
        
    def test_final_filter_with_string_columns(self, monkeypatch, tmp_path):
        """Test that final filter works with string column comparisons."""
        # Create minimal test files
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t1000\t.\tA\tT\t100\tPASS\t.\tGT\t0/1
chr2\t2000\t.\tG\tC\t100\tPASS\t.\tGT\t1/1
chr3\t3000\t.\tC\tG\t100\tPASS\t.\tGT\t0/1
""")
        
        output_dir = tmp_path / "output"
        output_file = output_dir / "test.final.tsv"
        
        # Mock external tools
        def mock_check_tools(*args, **kwargs):
            return True
            
        def mock_run_command(cmd, output_file=None):
            if "genes2bed" in cmd:
                bed_file = cmd[-1]
                with open(bed_file, "w") as f:
                    f.write("chr1\t500\t1500\tGENE1\n")
                    f.write("chr2\t1500\t2500\tGENE2\n")
                    f.write("chr3\t2500\t3500\tGENE3\n")
            elif "bcftools" in cmd and "view" in cmd:
                output = cmd[cmd.index("-o") + 1]
                vcf_file.rename(output)
            elif "SnpSift" in cmd and "extractFields" in cmd:
                output = output_file or cmd[-1]
                with open(output, "w") as f:
                    f.write("CHROM\tPOS\tREF\tALT\tGENE\tIMPACT\tInheritance_Pattern\n")
                    f.write("chr1\t1000\tA\tT\tGENE1\tHIGH\tde_novo\n")
                    f.write("chr2\t2000\tG\tC\tGENE2\tMODERATE\tautosomal_recessive\n")
                    f.write("chr3\t3000\tC\tG\tGENE3\tHIGH\tcompound_heterozygous\n")
            elif "bgzip" in cmd or "bcftools" in cmd and "index" in cmd:
                pass
            return 0
            
        # Mock subprocess.run for gene_bed module
        def mock_subprocess_run(cmd, **kwargs):
            if "genes2bed" in cmd:
                # Write to the output file specified in stdout
                if 'stdout' in kwargs and hasattr(kwargs['stdout'], 'write'):
                    kwargs['stdout'].write("chr1\t500\t3000\tTEST_GENE\n")
            elif "sort" in cmd:
                # For sort command, just pass through
                pass
            return MagicMock(returncode=0)
        
        monkeypatch.setattr("subprocess.run", mock_subprocess_run)
        monkeypatch.setattr("variantcentrifuge.pipeline.check_external_tools", mock_check_tools)
        monkeypatch.setattr("variantcentrifuge.utils.run_command", mock_run_command)
        monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)
        monkeypatch.setattr("variantcentrifuge.extractor.run_command", mock_run_command)
        
        # Run with filter for specific inheritance patterns
        test_args = [
            "variantcentrifuge",
            "--vcf-file", str(vcf_file),
            "--gene-name", "GENE1 GENE2 GENE3",
            "--output-dir", str(output_dir),
            "--output-file", "test.final.tsv",
            "--filters", "FILTER='PASS'",
            "--fields", "CHROM POS REF ALT GENE IMPACT Inheritance_Pattern",
            "--reference", "GRCh38",
            "--no-replacement",
            "--no-stats",
            "--final-filter", 'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]'
        ]
        
        with patch("sys.argv", test_args):
            main()
        
        # Check results
        assert output_file.exists()
        
        with open(output_file) as f:
            lines = f.readlines()
            
        # Should have header + 2 data rows (filtering out autosomal_recessive)
        assert len(lines) == 3
        assert "de_novo" in "".join(lines)
        assert "compound_heterozygous" in "".join(lines)
        assert "autosomal_recessive" not in "".join(lines)
        
    def test_final_filter_complex_expression(self, monkeypatch, tmp_path):
        """Test complex filter expressions with AND/OR logic."""
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t1000\t.\tA\tT\t100\tPASS\t.\tGT\t0/1
""")
        
        output_dir = tmp_path / "output"
        output_file = output_dir / "test.final.tsv"
        
        def mock_check_tools(*args, **kwargs):
            return True
            
        def mock_run_command(cmd, output_file=None):
            if "genes2bed" in cmd:
                bed_file = cmd[-1]
                with open(bed_file, "w") as f:
                    f.write("chr1\t500\t5000\tBRCA1\n")
            elif "bcftools" in cmd and "view" in cmd:
                output = cmd[cmd.index("-o") + 1]
                vcf_file.rename(output)
            elif "SnpSift" in cmd and "extractFields" in cmd:
                output = output_file or cmd[-1]
                with open(output, "w") as f:
                    f.write("CHROM\tPOS\tGENE\tIMPACT\tinheritance_score\tCADD_score\n")
                    f.write("chr1\t1000\tBRCA1\tHIGH\t0.95\t25.5\n")
                    f.write("chr1\t2000\tBRCA1\tMODERATE\t0.4\t15.0\n")
                    f.write("chr1\t3000\tBRCA1\tHIGH\t0.3\t30.0\n")
                    f.write("chr1\t4000\tBRCA1\tLOW\t0.9\t10.0\n")
            elif "bgzip" in cmd or "bcftools" in cmd and "index" in cmd:
                pass
            return 0
            
        # Mock subprocess.run for gene_bed module
        def mock_subprocess_run(cmd, **kwargs):
            if "genes2bed" in cmd:
                # Write to the output file specified in stdout
                if 'stdout' in kwargs and hasattr(kwargs['stdout'], 'write'):
                    kwargs['stdout'].write("chr1\t500\t3000\tTEST_GENE\n")
            elif "sort" in cmd:
                # For sort command, just pass through
                pass
            return MagicMock(returncode=0)
        
        monkeypatch.setattr("subprocess.run", mock_subprocess_run)
        monkeypatch.setattr("variantcentrifuge.pipeline.check_external_tools", mock_check_tools)
        monkeypatch.setattr("variantcentrifuge.utils.run_command", mock_run_command)
        monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)
        monkeypatch.setattr("variantcentrifuge.extractor.run_command", mock_run_command)
        
        # Complex filter: (HIGH impact AND high score) OR (very high CADD)
        test_args = [
            "variantcentrifuge",
            "--vcf-file", str(vcf_file),
            "--gene-name", "BRCA1",
            "--output-dir", str(output_dir),
            "--output-file", "test.final.tsv",
            "--filters", "FILTER='PASS'",
            "--fields", "CHROM POS GENE IMPACT inheritance_score CADD_score",
            "--reference", "GRCh38",
            "--no-replacement",
            "--no-stats",
            "--final-filter", '(IMPACT == "HIGH" and inheritance_score > 0.5) or CADD_score > 28'
        ]
        
        with patch("sys.argv", test_args):
            main()
        
        assert output_file.exists()
        
        with open(output_file) as f:
            lines = f.readlines()
            
        # Should retain:
        # - Row 1: HIGH impact AND score 0.95 > 0.5
        # - Row 3: CADD 30.0 > 28
        # Should filter out:
        # - Row 2: MODERATE impact
        # - Row 4: LOW impact and CADD 10.0 < 28
        assert len(lines) == 3  # header + 2 matching rows
        data = "".join(lines)
        assert "chr1\t1000" in data  # HIGH + high score
        assert "chr1\t3000" in data  # High CADD
        assert "chr1\t2000" not in data  # Filtered out
        assert "chr1\t4000" not in data  # Filtered out