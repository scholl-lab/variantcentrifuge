"""
Regression tests comparing old and new pipeline outputs.

This module ensures that the refactored pipeline produces identical
results to the original monolithic pipeline.
"""

import tempfile
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import pytest

from variantcentrifuge.cli import parse_args
from variantcentrifuge.config import load_config
from variantcentrifuge.pipeline import run_pipeline


def compare_tsv_files(file1: Path, file2: Path, ignore_columns: List[str] = None) -> bool:
    """Compare two TSV files, optionally ignoring certain columns."""
    df1 = pd.read_csv(file1, sep="\t", dtype=str)
    df2 = pd.read_csv(file2, sep="\t", dtype=str)

    # Check if DataFrames have same shape
    if df1.shape != df2.shape:
        print(f"Shape mismatch: {file1} has {df1.shape}, {file2} has {df2.shape}")
        return False

    # Check if columns match
    if set(df1.columns) != set(df2.columns):
        print(f"Column mismatch: {file1} has {set(df1.columns)}, {file2} has {set(df2.columns)}")
        return False

    # Ignore specified columns
    if ignore_columns:
        cols_to_compare = [col for col in df1.columns if col not in ignore_columns]
        df1 = df1[cols_to_compare]
        df2 = df2[cols_to_compare]

    # Sort both DataFrames to ensure consistent ordering
    sort_cols = (
        ["CHROM", "POS", "REF", "ALT"]
        if all(col in df1.columns for col in ["CHROM", "POS", "REF", "ALT"])
        else df1.columns.tolist()
    )
    df1_sorted = df1.sort_values(by=sort_cols).reset_index(drop=True)
    df2_sorted = df2.sort_values(by=sort_cols).reset_index(drop=True)

    # Compare DataFrames
    try:
        pd.testing.assert_frame_equal(df1_sorted, df2_sorted, check_dtype=False)
        return True
    except AssertionError as e:
        print(f"DataFrame comparison failed: {e}")
        return False


class TestPipelineRegression:
    """Test that new pipeline produces identical output to old pipeline."""

    @pytest.fixture
    def test_data_dir(self):
        """Path to test data directory."""
        return Path(__file__).parent / "data"

    @pytest.fixture
    def small_vcf(self, test_data_dir):
        """Path to small test VCF file."""
        vcf_path = test_data_dir / "small_test.vcf.gz"
        if not vcf_path.exists():
            pytest.skip(f"Test VCF not found: {vcf_path}")
        return vcf_path

    def run_pipeline_comparison(
        self, args_list: List[str], compare_func=compare_tsv_files, ignore_columns: List[str] = None
    ) -> Tuple[bool, str]:
        """Run both pipelines and compare outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run old pipeline
            old_output_dir = Path(tmpdir) / "old_pipeline"
            old_output_dir.mkdir()

            old_args = args_list + [
                "--output-dir",
                str(old_output_dir),
            ]

            # Parse args and load config
            args = parse_args(old_args)
            cfg = load_config(args.config) if args.config else {}

            # Merge CLI args into config
            cfg.update(vars(args))
            cfg["use_new_pipeline_architecture"] = False

            # Run old pipeline
            import datetime

            start_time = datetime.datetime.now()
            run_pipeline(args, cfg, start_time)

            # Run new pipeline
            new_output_dir = Path(tmpdir) / "new_pipeline"
            new_output_dir.mkdir()

            new_args = args_list + [
                "--output-dir",
                str(new_output_dir),
            ]

            # Parse args and load config
            args = parse_args(new_args)
            cfg = load_config(args.config) if args.config else {}

            # Merge CLI args into config
            cfg.update(vars(args))
            cfg["use_new_pipeline_architecture"] = True

            # Run new pipeline
            start_time = datetime.datetime.now()
            run_pipeline(args, cfg, start_time)

            # Compare outputs
            old_files = sorted(old_output_dir.glob("*.tsv"))
            new_files = sorted(new_output_dir.glob("*.tsv"))

            if len(old_files) != len(new_files):
                return (
                    False,
                    f"Different number of output files: old={len(old_files)}, new={len(new_files)}",
                )

            for old_file, new_file in zip(old_files, new_files):
                if old_file.name != new_file.name:
                    return False, f"File name mismatch: {old_file.name} vs {new_file.name}"

                if not compare_func(old_file, new_file, ignore_columns):
                    return False, f"File content mismatch: {old_file.name}"

            return True, "All outputs match"

    @pytest.mark.slow
    def test_basic_gene_extraction(self, small_vcf):
        """Test basic gene extraction produces identical output."""
        args = [
            "--vcf-file",
            str(small_vcf),
            "--gene-name",
            "BRCA1",
            "--fields",
            "CHROM,POS,REF,ALT,QUAL",
        ]

        success, message = self.run_pipeline_comparison(args)
        assert success, message

    @pytest.mark.slow
    def test_with_filtering(self, small_vcf):
        """Test pipeline with filtering produces identical output."""
        args = [
            "--vcf-file",
            str(small_vcf),
            "--gene-name",
            "BRCA1",
            "--fields",
            "CHROM,POS,REF,ALT,QUAL,ANN[*].GENE",
            "--filter",
            "QUAL >= 30",
        ]

        success, message = self.run_pipeline_comparison(args)
        assert success, message

    @pytest.mark.slow
    def test_with_genotype_replacement(self, small_vcf):
        """Test pipeline with genotype replacement produces identical output."""
        args = [
            "--vcf-file",
            str(small_vcf),
            "--gene-name",
            "BRCA1",
            "--fields",
            "CHROM,POS,REF,ALT,QUAL,GEN[*].GT",
            "--replace-genotypes",
        ]

        success, message = self.run_pipeline_comparison(args)
        assert success, message

    @pytest.mark.slow
    def test_parallel_extraction(self, small_vcf):
        """Test parallel variant extraction produces identical output."""
        args = [
            "--vcf-file",
            str(small_vcf),
            "--gene-name",
            "BRCA1,BRCA2,TP53",
            "--fields",
            "CHROM,POS,REF,ALT,QUAL",
            "--threads",
            "4",
        ]

        success, message = self.run_pipeline_comparison(args)
        assert success, message

    @pytest.mark.slow
    @pytest.mark.parametrize("preset", ["rare", "coding", "pathogenic"])
    def test_with_presets(self, small_vcf, preset):
        """Test pipeline with different presets produces identical output."""
        args = [
            "--vcf-file",
            str(small_vcf),
            "--gene-name",
            "BRCA1",
            "--preset",
            preset,
            "--fields",
            "CHROM,POS,REF,ALT,QUAL,ANN[*].IMPACT",
        ]

        success, message = self.run_pipeline_comparison(args)
        assert success, message


if __name__ == "__main__":
    # Allow running directly for debugging
    import sys

    pytest.main([__file__] + sys.argv[1:])
