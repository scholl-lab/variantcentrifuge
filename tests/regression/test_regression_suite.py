"""Regression test suite for validating the new pipeline against the old pipeline.

This test suite runs both pipelines with identical inputs and validates that outputs
are identical, ensuring no regression in functionality.
"""

import hashlib
import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pytest

logger = logging.getLogger(__name__)


class RegressionTestConfig:
    """Configuration for a regression test case."""

    def __init__(
        self,
        name: str,
        gene_name: Optional[str] = None,
        gene_file: Optional[str] = None,
        preset: Optional[str] = None,
        extra_args: Optional[List[str]] = None,
        description: str = "",
    ):
        self.name = name
        self.gene_name = gene_name
        self.gene_file = gene_file
        self.preset = preset
        self.extra_args = extra_args or []
        self.description = description


# Define regression test cases
REGRESSION_TESTS = [
    RegressionTestConfig(
        name="single_gene_basic",
        gene_name="BRCA1",
        description="Basic single gene extraction",
    ),
    RegressionTestConfig(
        name="single_gene_rare_coding",
        gene_name="TP53",
        preset="rare,coding",
        description="Single gene with rare coding filter preset",
    ),
    RegressionTestConfig(
        name="multi_gene_file",
        gene_file="test_genes.txt",
        description="Multiple genes from file",
    ),
    RegressionTestConfig(
        name="gene_with_scoring",
        gene_name="CFTR",
        extra_args=["--scoring-config-path", "scoring/inheritance_score"],
        description="Gene with variant scoring",
    ),
    RegressionTestConfig(
        name="gene_with_inheritance",
        gene_name="PKD1",
        extra_args=["--ped", "test_family.ped", "--inheritance-mode", "columns"],
        description="Gene with inheritance analysis",
    ),
    RegressionTestConfig(
        name="gene_with_annotations",
        gene_name="LDLR",
        extra_args=["--annotate-bed", "test_regions.bed"],
        description="Gene with custom BED annotations",
    ),
    RegressionTestConfig(
        name="complex_filtering",
        gene_name="APOB",
        preset="pathogenic",
        extra_args=["--final-filter", "AF < 0.001"],
        description="Complex filtering with preset and final filter",
    ),
]


class PipelineRunner:
    """Runs the pipeline and captures outputs."""

    def __init__(self, use_new_pipeline: bool = False):
        self.use_new_pipeline = use_new_pipeline

    def run(
        self,
        vcf_file: str,
        output_dir: Path,
        config: RegressionTestConfig,
    ) -> Dict[str, Path]:
        """Run the pipeline with given configuration.

        Returns dict of output file paths.
        """
        output_file = output_dir / f"{config.name}_output.tsv"

        # Build command
        cmd = [
            "variantcentrifuge",
            "--vcf-file",
            vcf_file,
            "--output-file",
            str(output_file),
        ]

        if self.use_new_pipeline:
            cmd.append("--use-new-pipeline")

        if config.gene_name:
            cmd.extend(["--gene-name", config.gene_name])
        elif config.gene_file:
            cmd.extend(["--gene-file", config.gene_file])

        if config.preset:
            cmd.extend(["--preset", config.preset])

        cmd.extend(config.extra_args)

        # Add common args for reproducibility
        cmd.extend(
            [
                "--threads",
                "1",  # Single thread for deterministic results
                "--no-metadata",  # Skip metadata to avoid timestamp differences
            ]
        )

        # Run pipeline
        logger.info(
            f"Running {'new' if self.use_new_pipeline else 'old'} pipeline: {' '.join(cmd)}"
        )

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(output_dir),
                check=True,
            )
            logger.debug(f"Pipeline stdout: {result.stdout}")

            if result.stderr:
                logger.warning(f"Pipeline stderr: {result.stderr}")

        except subprocess.CalledProcessError as e:
            logger.error(f"Pipeline failed: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise

        # Collect output files
        outputs = {
            "tsv": output_file,
        }

        # Check for additional output files
        stats_file = output_dir / f"{config.name}_output_stats.json"
        if stats_file.exists():
            outputs["stats"] = stats_file

        return outputs


class OutputValidator:
    """Validates that outputs from both pipelines match."""

    def compare_tsv_files(
        self,
        old_file: Path,
        new_file: Path,
        ignore_columns: Optional[List[str]] = None,
    ) -> Tuple[bool, List[str]]:
        """Compare TSV files, returning (match, differences)."""
        ignore_columns = ignore_columns or []
        differences = []

        # Read files
        try:
            old_df = pd.read_csv(old_file, sep="\t", dtype=str, na_values=[])
            new_df = pd.read_csv(new_file, sep="\t", dtype=str, na_values=[])
        except Exception as e:
            return False, [f"Failed to read files: {e}"]

        # Check shape
        if old_df.shape != new_df.shape:
            differences.append(f"Shape mismatch: old={old_df.shape}, new={new_df.shape}")
            return False, differences

        # Check columns
        old_cols = set(old_df.columns) - set(ignore_columns)
        new_cols = set(new_df.columns) - set(ignore_columns)

        if old_cols != new_cols:
            missing = old_cols - new_cols
            extra = new_cols - old_cols
            if missing:
                differences.append(f"Missing columns: {missing}")
            if extra:
                differences.append(f"Extra columns: {extra}")
            return False, differences

        # Sort both dataframes by key columns for consistent comparison
        key_cols = ["CHROM", "POS", "REF", "ALT"]
        if all(col in old_df.columns for col in key_cols):
            old_df = old_df.sort_values(key_cols).reset_index(drop=True)
            new_df = new_df.sort_values(key_cols).reset_index(drop=True)

        # Compare values column by column
        for col in sorted(old_cols):
            if col in ignore_columns:
                continue

            # Handle NaN comparison
            old_vals = old_df[col].fillna("")
            new_vals = new_df[col].fillna("")

            if not old_vals.equals(new_vals):
                # Find first difference
                for idx in range(len(old_vals)):
                    if old_vals.iloc[idx] != new_vals.iloc[idx]:
                        differences.append(
                            f"Column '{col}' differs at row {idx}: "
                            f"old='{old_vals.iloc[idx]}', new='{new_vals.iloc[idx]}'"
                        )
                        break

                # Count total differences
                n_diff = (old_vals != new_vals).sum()
                differences.append(f"Column '{col}' has {n_diff} differences")

        return len(differences) == 0, differences

    def compare_stats_files(
        self,
        old_file: Path,
        new_file: Path,
    ) -> Tuple[bool, List[str]]:
        """Compare statistics JSON files."""
        differences = []

        try:
            with open(old_file) as f:
                old_stats = json.load(f)
            with open(new_file) as f:
                new_stats = json.load(f)
        except Exception as e:
            return False, [f"Failed to read stats files: {e}"]

        # Compare keys
        if set(old_stats.keys()) != set(new_stats.keys()):
            differences.append(
                f"Stats keys differ: old={set(old_stats.keys())}, " f"new={set(new_stats.keys())}"
            )
            return False, differences

        # Compare values (with tolerance for floating point)
        for key in old_stats:
            old_val = old_stats[key]
            new_val = new_stats[key]

            if isinstance(old_val, (int, float)) and isinstance(new_val, (int, float)):
                # Numeric comparison with tolerance
                if abs(old_val - new_val) > 1e-6:
                    differences.append(f"Stats '{key}' differs: old={old_val}, new={new_val}")
            elif old_val != new_val:
                differences.append(f"Stats '{key}' differs: old={old_val}, new={new_val}")

        return len(differences) == 0, differences

    def calculate_file_hash(self, file_path: Path) -> str:
        """Calculate SHA256 hash of file contents."""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()


@pytest.fixture(scope="module")
def test_data_dir():
    """Directory containing test data files."""
    # Look for test data in several locations
    possible_dirs = [
        Path(__file__).parent / "data",
        Path(__file__).parent.parent / "data",
        Path(__file__).parent.parent.parent / "test_data",
    ]

    for dir_path in possible_dirs:
        if dir_path.exists():
            return dir_path

    pytest.skip("Test data directory not found")


@pytest.fixture(scope="module")
def test_vcf_file(test_data_dir):
    """Path to test VCF file."""
    vcf_file = test_data_dir / "test_cohort.vcf.gz"
    if not vcf_file.exists():
        pytest.skip(f"Test VCF file not found: {vcf_file}")
    return str(vcf_file)


@pytest.mark.regression
class TestPipelineRegression:
    """Regression tests comparing old and new pipeline outputs."""

    def setup_method(self):
        """Set up test environment."""
        self.old_runner = PipelineRunner(use_new_pipeline=False)
        self.new_runner = PipelineRunner(use_new_pipeline=True)
        self.validator = OutputValidator()

    @pytest.mark.parametrize("config", REGRESSION_TESTS, ids=lambda c: c.name)
    def test_pipeline_output_match(self, config, test_vcf_file, tmp_path):
        """Test that new pipeline produces identical output to old pipeline."""
        # Skip if dependencies missing
        self._check_test_dependencies(config, tmp_path)

        # Create output directories
        old_output_dir = tmp_path / "old_pipeline"
        new_output_dir = tmp_path / "new_pipeline"
        old_output_dir.mkdir()
        new_output_dir.mkdir()

        # Run both pipelines
        try:
            old_outputs = self.old_runner.run(test_vcf_file, old_output_dir, config)
            new_outputs = self.new_runner.run(test_vcf_file, new_output_dir, config)
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Pipeline execution failed: {e}")

        # Compare outputs
        assert set(old_outputs.keys()) == set(
            new_outputs.keys()
        ), f"Output files differ: old={set(old_outputs.keys())}, new={set(new_outputs.keys())}"

        # Compare TSV files
        tsv_match, tsv_diffs = self.validator.compare_tsv_files(
            old_outputs["tsv"],
            new_outputs["tsv"],
            ignore_columns=["Variant_ID"],  # May differ due to processing order
        )

        if not tsv_match:
            # Save files for debugging
            debug_dir = tmp_path / "debug"
            debug_dir.mkdir()

            import shutil

            shutil.copy(old_outputs["tsv"], debug_dir / "old_output.tsv")
            shutil.copy(new_outputs["tsv"], debug_dir / "new_output.tsv")

            pytest.fail(
                f"TSV outputs differ for {config.name}:\n"
                + "\n".join(tsv_diffs)
                + f"\nDebug files saved to: {debug_dir}"
            )

        # Compare stats if present
        if "stats" in old_outputs:
            stats_match, stats_diffs = self.validator.compare_stats_files(
                old_outputs["stats"],
                new_outputs["stats"],
            )

            if not stats_match:
                pytest.fail(f"Statistics differ for {config.name}:\n" + "\n".join(stats_diffs))

        # Log success
        logger.info(f"✓ Regression test passed: {config.name}")

    def test_deterministic_output(self, test_vcf_file, tmp_path):
        """Test that pipeline produces deterministic output across runs."""
        config = RegressionTestConfig(
            name="deterministic_test",
            gene_name="BRCA2",
            preset="rare",
        )

        # Run new pipeline twice
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        run1_dir.mkdir()
        run2_dir.mkdir()

        outputs1 = self.new_runner.run(test_vcf_file, run1_dir, config)
        outputs2 = self.new_runner.run(test_vcf_file, run2_dir, config)

        # Compare outputs
        tsv_match, diffs = self.validator.compare_tsv_files(
            outputs1["tsv"],
            outputs2["tsv"],
        )

        assert tsv_match, f"Pipeline output is not deterministic:\n" + "\n".join(diffs)

        # Compare file hashes
        hash1 = self.validator.calculate_file_hash(outputs1["tsv"])
        hash2 = self.validator.calculate_file_hash(outputs2["tsv"])

        assert hash1 == hash2, f"File hashes differ: {hash1} != {hash2}"

    def _check_test_dependencies(self, config: RegressionTestConfig, tmp_path: Path):
        """Check if test dependencies are available."""
        # Create dummy files if needed for testing
        if config.gene_file:
            gene_file = tmp_path / config.gene_file
            if not gene_file.exists():
                gene_file.write_text("BRCA1\nTP53\nEGFR\n")

        if "--ped" in config.extra_args:
            ped_idx = config.extra_args.index("--ped") + 1
            ped_file = tmp_path / config.extra_args[ped_idx]
            if not ped_file.exists():
                # Create minimal PED file
                ped_file.write_text(
                    "FAM001\tCHILD\tFATHER\tMOTHER\t1\t2\n"
                    "FAM001\tFATHER\t0\t0\t1\t1\n"
                    "FAM001\tMOTHER\t0\t0\t2\t1\n"
                )

        if "--annotate-bed" in config.extra_args:
            bed_idx = config.extra_args.index("--annotate-bed") + 1
            bed_file = tmp_path / config.extra_args[bed_idx]
            if not bed_file.exists():
                # Create minimal BED file
                bed_file.write_text(
                    "chr1\t1000000\t2000000\tTestRegion1\n" "chr2\t5000000\t6000000\tTestRegion2\n"
                )


if __name__ == "__main__":
    # Allow running as script for debugging
    pytest.main([__file__, "-v", "-s", "-m", "regression"])
