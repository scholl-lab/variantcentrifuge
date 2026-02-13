"""Unit tests for gene list annotation functionality."""

import tempfile
from pathlib import Path

from variantcentrifuge.helpers import (
    _sanitize_column_name,
    annotate_variants_with_gene_lists,
    load_gene_list,
)


class TestGeneListLoading:
    """Test gene list loading functionality."""

    def test_load_empty_file(self, tmp_path):
        """Test loading an empty gene list file."""
        # Create empty file
        empty_file = tmp_path / "empty.txt"
        empty_file.write_text("")

        gene_set = load_gene_list(str(empty_file))
        assert len(gene_set) == 0
        assert isinstance(gene_set, set)

    def test_load_basic_gene_list(self, tmp_path):
        """Test loading a basic gene list with one gene per line."""
        genes = ["TP53", "BRCA1", "BRCA2", "APC"]
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("\n".join(genes))

        gene_set = load_gene_list(str(gene_file))
        assert len(gene_set) == 4
        for gene in genes:
            assert gene.upper() in gene_set

    def test_load_mixed_case_genes(self, tmp_path):
        """Test that genes are loaded case-insensitively."""
        genes = ["tp53", "Brca1", "BRCA2", "aPc"]
        expected = {"TP53", "BRCA1", "BRCA2", "APC"}  # All uppercase
        gene_file = tmp_path / "mixed_case.txt"
        gene_file.write_text("\n".join(genes))

        gene_set = load_gene_list(str(gene_file))
        assert gene_set == expected

    def test_load_with_empty_lines(self, tmp_path):
        """Test loading a gene list with empty lines."""
        content = "TP53\n\nBRCA1\n\n\nBRCA2\nAPC\n"
        gene_file = tmp_path / "with_empty.txt"
        gene_file.write_text(content)

        gene_set = load_gene_list(str(gene_file))
        assert len(gene_set) == 4
        assert gene_set == {"TP53", "BRCA1", "BRCA2", "APC"}

    def test_load_with_whitespace(self, tmp_path):
        """Test loading a gene list with whitespace."""
        content = "  TP53  \n BRCA1 \nBRCA2\t\nAPC"
        gene_file = tmp_path / "with_whitespace.txt"
        gene_file.write_text(content)

        gene_set = load_gene_list(str(gene_file))
        assert len(gene_set) == 4
        assert gene_set == {"TP53", "BRCA1", "BRCA2", "APC"}

    def test_file_not_found(self):
        """Test handling of a non-existent file."""
        gene_set = load_gene_list("/non/existent/path.txt")
        assert len(gene_set) == 0
        assert isinstance(gene_set, set)


class TestColumnNameSanitization:
    """Test column name sanitization functionality."""

    def test_basic_filename(self):
        """Test basic filename sanitization."""
        path = "/path/to/my_gene_list.txt"
        assert _sanitize_column_name(path) == "my_gene_list"

    def test_strip_multiple_extensions(self):
        """Test stripping of common gene list extensions."""
        assert _sanitize_column_name("/path/cancer_genes.txt") == "cancer_genes"
        assert _sanitize_column_name("/path/cancer_genes.genes") == "cancer_genes"
        assert _sanitize_column_name("/path/cancer_genes.gene.list") == "cancer_genes"

    def test_invalid_characters(self):
        """Test replacement of invalid characters."""
        assert _sanitize_column_name("/path/gene@list#1.txt") == "gene_list_1"
        assert _sanitize_column_name("/path/gene list 1.txt") == "gene_list_1"

    def test_purely_numeric(self):
        """Test handling of purely numeric filenames."""
        assert _sanitize_column_name("/path/123.txt") == "list_123"

    def test_empty_name(self):
        """Test handling of filenames that sanitize to empty."""
        assert _sanitize_column_name("/path/-----.txt") == "unnamed_list"

    def test_starting_with_number(self):
        """Test handling of filenames that start with a number."""
        assert _sanitize_column_name("/path/123genes.txt") == "list_123genes"


class TestVariantAnnotation:
    """Test variant TSV annotation with gene lists."""

    def setup_method(self):
        """Set up test data."""
        self.header = "CHROM\tPOS\tREF\tALT\tGENE\tIMPACT"
        self.variants = [
            "chr1\t1000\tA\tG\tTP53\tHIGH",
            "chr1\t2000\tG\tT\tBRCA1\tMODERATE",
            "chr2\t3000\tC\tA\tBRCA2\tLOW",
            "chr3\t4000\tT\tC\tAPC,TP53\tMODERATE",  # Multiple genes
            "chr4\t5000\tA\tG\tunknown\tLOW",
            "chr5\t6000\tG\tC\t\tLOW",  # Empty GENE
        ]
        self.tsv_lines = [self.header, *self.variants]

    def create_gene_list_file(self, tmp_path, name, genes):
        """Create a gene list file."""
        file_path = tmp_path / f"{name}.txt"
        file_path.write_text("\n".join(genes))
        return str(file_path)

    def test_no_gene_lists(self):
        """Test annotation with no gene lists."""
        result = annotate_variants_with_gene_lists(self.tsv_lines, [])
        assert result == self.tsv_lines

    def test_empty_tsv(self):
        """Test annotation with empty TSV."""
        gene_list_file = self.create_gene_list_file(Path(tempfile.gettempdir()), "test", ["TP53"])
        result = annotate_variants_with_gene_lists([], [gene_list_file])
        assert result == []

    def test_single_gene_list(self, tmp_path):
        """Test annotation with a single gene list."""
        cancer_genes = ["TP53", "BRCA1"]
        cancer_gene_file = self.create_gene_list_file(tmp_path, "cancer_genes", cancer_genes)

        result = annotate_variants_with_gene_lists(self.tsv_lines, [cancer_gene_file])

        # Check header has new column
        assert result[0] == f"{self.header}\tcancer_genes"

        # Check annotations
        assert result[1].endswith("\tyes")  # TP53
        assert result[2].endswith("\tyes")  # BRCA1
        assert result[3].endswith("\tno")  # BRCA2
        assert result[4].endswith("\tyes")  # APC,TP53
        assert result[5].endswith("\tno")  # unknown
        assert result[6].endswith("\tno")  # empty

    def test_multiple_gene_lists(self, tmp_path):
        """Test annotation with multiple gene lists."""
        cancer_genes = ["TP53", "BRCA1"]
        cancer_gene_file = self.create_gene_list_file(tmp_path, "cancer_genes", cancer_genes)

        apc_genes = ["APC"]
        apc_gene_file = self.create_gene_list_file(tmp_path, "apc_genes", apc_genes)

        result = annotate_variants_with_gene_lists(
            self.tsv_lines, [cancer_gene_file, apc_gene_file]
        )

        # Check header has new columns
        assert result[0] == f"{self.header}\tcancer_genes\tapc_genes"

        # Check cancer genes column
        assert result[1].split("\t")[-2] == "yes"  # TP53
        assert result[2].split("\t")[-2] == "yes"  # BRCA1
        assert result[3].split("\t")[-2] == "no"  # BRCA2
        assert result[4].split("\t")[-2] == "yes"  # APC,TP53

        # Check APC genes column
        assert result[1].split("\t")[-1] == "no"  # TP53
        assert result[2].split("\t")[-1] == "no"  # BRCA1
        assert result[3].split("\t")[-1] == "no"  # BRCA2
        assert result[4].split("\t")[-1] == "yes"  # APC,TP53

    def test_case_insensitive_matching(self, tmp_path):
        """Test case-insensitive gene matching."""
        # Create variants with differently cased genes
        variants_with_case = [
            "chr1\t1000\tA\tG\ttp53\tHIGH",  # lowercase
            "chr1\t2000\tG\tT\tBrCa1\tMODERATE",  # mixed case
        ]
        tsv_lines = [self.header, *variants_with_case]

        # Create gene list with uppercase genes
        cancer_genes = ["TP53", "BRCA1"]
        cancer_gene_file = self.create_gene_list_file(tmp_path, "cancer_genes", cancer_genes)

        result = annotate_variants_with_gene_lists(tsv_lines, [cancer_gene_file])

        # Both should match despite case differences
        assert result[1].endswith("\tyes")  # tp53 matches TP53
        assert result[2].endswith("\tyes")  # BrCa1 matches BRCA1

    def test_multple_gene_separators(self, tmp_path):
        """Test handling of different gene separators."""
        # Create variants with differently separated genes
        variants_with_separators = [
            "chr1\t1000\tA\tG\tTP53,BRCA1\tHIGH",  # comma
            "chr1\t2000\tG\tT\tBRCA2;APC\tMODERATE",  # semicolon
            "chr2\t3000\tC\tA\tTP53 BRCA2\tLOW",  # space
            "chr3\t4000\tT\tC\tAPC, TP53; BRCA1\tLOW",  # mixed
        ]
        tsv_lines = [self.header, *variants_with_separators]

        # Create gene lists
        cancer_genes = ["TP53", "BRCA1"]
        cancer_gene_file = self.create_gene_list_file(tmp_path, "cancer_genes", cancer_genes)

        result = annotate_variants_with_gene_lists(tsv_lines, [cancer_gene_file])

        # All should have matches
        assert result[1].endswith("\tyes")  # TP53,BRCA1
        assert result[2].endswith("\tno")  # BRCA2;APC
        assert result[3].endswith("\tyes")  # TP53 BRCA2
        assert result[4].endswith("\tyes")  # APC, TP53; BRCA1

    def test_duplicate_column_names(self, tmp_path):
        """Test handling of duplicate column names."""
        # Create two files that would sanitize to the same column name
        genes1 = ["TP53"]
        genes2 = ["BRCA1"]

        # Both would sanitize to "genes"
        file1 = self.create_gene_list_file(tmp_path, "genes", genes1)
        file2 = self.create_gene_list_file(tmp_path, "genes.txt", genes2)

        result = annotate_variants_with_gene_lists(self.tsv_lines, [file1, file2])

        # Header should have uniquified column names
        assert result[0] == f"{self.header}\tgenes\tgenes_1"

        # TP53 variant should match first list but not second
        assert result[1].split("\t")[-2] == "yes"
        assert result[1].split("\t")[-1] == "no"

        # BRCA1 variant should match second list but not first
        assert result[2].split("\t")[-2] == "no"
        assert result[2].split("\t")[-1] == "yes"

    def test_missing_gene_column(self, tmp_path):
        """Test behavior when GENE column is missing."""
        # Create TSV without GENE column
        header = "CHROM\tPOS\tREF\tALT\tIMPACT"
        variants = [
            "chr1\t1000\tA\tG\tHIGH",
            "chr1\t2000\tG\tT\tMODERATE",
        ]
        tsv_lines = [header, *variants]

        cancer_genes = ["TP53", "BRCA1"]
        cancer_gene_file = self.create_gene_list_file(tmp_path, "cancer_genes", cancer_genes)

        # Should return original lines unchanged with a warning
        result = annotate_variants_with_gene_lists(tsv_lines, [cancer_gene_file])
        assert result == tsv_lines


# Integration test will be added to test the pipeline integration
