"""Regression test suite for validating the new pipeline against the old pipeline.

This test suite runs both pipelines with identical inputs and validates that outputs
are identical, ensuring no regression in functionality.
"""

import hashlib
import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pytest

logger = logging.getLogger(__name__)

# Check for required external tools
REQUIRED_TOOLS = ["bcftools", "snpEff", "SnpSift", "bedtools"]
TOOLS_AVAILABLE = all(shutil.which(tool) is not None for tool in REQUIRED_TOOLS)

if not TOOLS_AVAILABLE:
    pytest.skip(
        (
            "Skipping regression tests: Missing required tools: "
            f"{[t for t in REQUIRED_TOOLS if not shutil.which(t)]}"
        ),
        allow_module_level=True,
    )


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
        extra_args=["--filters", "FILTER = 'PASS'"],
        description="Basic single gene extraction",
    ),
    RegressionTestConfig(
        name="single_gene_rare_coding",
        gene_name="TP53",
        extra_args=["--filters", "(ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE')"],
        description="Single gene with coding filter",
    ),
    RegressionTestConfig(
        name="multi_gene_file",
        gene_file="test_genes.txt",
        extra_args=["--filters", "FILTER = 'PASS'"],
        description="Multiple genes from file",
    ),
    RegressionTestConfig(
        name="gene_with_scoring",
        gene_name="CFTR",
        extra_args=[
            "--filters",
            "FILTER = 'PASS'",
            "--scoring-config-path",
            "scoring/inheritance_score",
        ],
        description="Gene with variant scoring",
    ),
    RegressionTestConfig(
        name="gene_with_inheritance",
        gene_name="PKD1",
        extra_args=[
            "--filters",
            "FILTER = 'PASS'",
            "--ped",
            "test_family.ped",
            "--inheritance-mode",
            "columns",
        ],
        description="Gene with inheritance analysis",
    ),
    RegressionTestConfig(
        name="gene_with_annotations",
        gene_name="LDLR",
        extra_args=["--filters", "FILTER = 'PASS'", "--annotate-bed", "test_regions.bed"],
        description="Gene with custom BED annotations",
    ),
    RegressionTestConfig(
        name="complex_filtering",
        gene_name="APOB",
        extra_args=["--filters", "ClinVar_CLNSIG =~ 'Pathogenic'", "--final-filter", "AF < 0.001"],
        description="Complex filtering with ClinVar and final filter",
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
        # Just use the filename, not the full path, since pipeline will put it in output_dir
        output_filename = f"{config.name}_output.tsv"
        output_file = output_dir / output_filename

        # Build command
        cmd = [
            "variantcentrifuge",
            "--vcf-file",
            vcf_file,
            "--output-file",
            output_filename,  # Just the filename, not full path
            "--output-dir",
            str(output_dir),
            "--config",
            str(Path(__file__).parent.parent / "fixtures" / "test_config.json"),
        ]

        if self.use_new_pipeline:
            cmd.append("--use-new-pipeline")

        if config.gene_name:
            cmd.extend(["--gene-name", config.gene_name])
        elif config.gene_file:
            gene_file_path = config.gene_file
            if not Path(gene_file_path).is_absolute():
                # Resolve relative path
                abs_path = Path(__file__).parent.parent.parent / gene_file_path
                if not abs_path.exists():
                    abs_path = Path(__file__).parent.parent / "fixtures" / gene_file_path
                gene_file_path = str(abs_path)
            cmd.extend(["--gene-file", gene_file_path])

        if config.preset:
            cmd.extend(["--preset", config.preset])

        # Fix relative paths in extra_args
        fixed_extra_args = []
        i = 0
        while i < len(config.extra_args):
            arg = config.extra_args[i]
            fixed_extra_args.append(arg)

            # Check if this is a file argument that needs path resolution
            if arg in ["--ped", "--annotate-bed", "--scoring-config-path", "--gene-file"]:
                if i + 1 < len(config.extra_args):
                    file_path = config.extra_args[i + 1]
                    if not Path(file_path).is_absolute():
                        # Try relative to project root first
                        abs_path = Path(__file__).parent.parent.parent / file_path
                        if not abs_path.exists():
                            # Try relative to tests/fixtures
                            abs_path = Path(__file__).parent.parent / "fixtures" / file_path
                        fixed_extra_args.append(str(abs_path))
                        i += 1  # Skip the file path we just processed
                else:
                    # This shouldn't happen with well-formed args
                    pass
            i += 1

        cmd.extend(fixed_extra_args)

        # Add common args for reproducibility
        cmd.extend(
            [
                "--threads",
                "1",  # Single thread for deterministic results
                "--fields",
                (
                    "CHROM POS REF ALT ID FILTER QUAL AC AF ANN[0].GENE ANN[0].EFFECT "
                    "ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P CADD_phred "
                    "dbNSFP_gnomAD_exomes_AF ClinVar_CLNSIG GEN[*].GT"
                ),
                "--keep-intermediates",  # Keep intermediate files for debugging
            ]
        )

        # Run pipeline
        # Filter out any None values before logging
        cmd = [str(arg) for arg in cmd if arg is not None]
        logger.info(
            f"Running {'new' if self.use_new_pipeline else 'old'} pipeline: {' '.join(cmd)}"
        )

        # Convert relative paths to absolute
        project_root = Path(__file__).parent.parent.parent

        # Update VCF path if relative
        if not Path(vcf_file).is_absolute():
            vcf_file = str(project_root / vcf_file)
            cmd[cmd.index("--vcf-file") + 1] = vcf_file

        logger.info(f"Running command from {project_root}: {' '.join(cmd)}")
        logger.info(f"Expected output file: {output_file}")
        logger.info(f"Output directory: {output_dir}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(project_root),  # Run from project root
                check=True,
            )
            logger.debug(f"Pipeline stdout: {result.stdout}")

            if result.stderr:
                logger.warning(f"Pipeline stderr: {result.stderr}")

            # Log directory contents after run
            if output_dir.exists():
                logger.info(f"Output directory contents: {list(output_dir.iterdir())}")
                for subdir in output_dir.iterdir():
                    if subdir.is_dir():
                        logger.info(f"  {subdir}: {list(subdir.iterdir())[:5]}...")  # First 5 files

        except subprocess.CalledProcessError as e:
            logger.error(f"Pipeline failed: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise

        # Collect output files
        # Check if the output file exists at the specified location
        if not output_file.exists():
            # Check if it's in the output directory with a different name
            possible_files = list(output_dir.glob("*.tsv"))
            if possible_files:
                logger.warning(
                    f"Output file not at expected location {output_file}, found: {possible_files}"
                )
                output_file = possible_files[0]  # Use the first TSV found
            else:
                # Check if it's in a subdirectory
                possible_files = list(output_dir.rglob("*.tsv"))
                if possible_files:
                    logger.warning(
                        (
                            f"Output file not at expected location {output_file}, "
                            f"found in subdirs: {possible_files}"
                        )
                    )
                    output_file = possible_files[0]
                else:
                    raise FileNotFoundError(f"No TSV output found in {output_dir}")

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
    # First try the regression test VCF in fixtures
    regression_vcf = Path("tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz")
    if regression_vcf.exists():
        return str(regression_vcf)

    # Fall back to test_cohort.vcf.gz in test data dir
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

        assert tsv_match, "Pipeline output is not deterministic:\n" + "\n".join(diffs)

        # Compare file hashes
        hash1 = self.validator.calculate_file_hash(outputs1["tsv"])
        hash2 = self.validator.calculate_file_hash(outputs2["tsv"])

        assert hash1 == hash2, f"File hashes differ: {hash1} != {hash2}"

    def _check_test_dependencies(self, config: RegressionTestConfig, tmp_path: Path):
        """Check if test dependencies are available."""
        # Update config to use test fixtures
        if config.gene_file == "test_genes.txt":
            config.gene_file = "tests/fixtures/test_genes.txt"

        # Update extra args to use test fixtures
        new_extra_args = []
        for i, arg in enumerate(config.extra_args):
            new_extra_args.append(arg)
            if arg == "--ped" and i + 1 < len(config.extra_args):
                if config.extra_args[i + 1] == "test_family.ped":
                    new_extra_args.append("tests/fixtures/test_family.ped")
                    i += 1  # Skip the next item
                else:
                    new_extra_args.append(config.extra_args[i + 1])
            elif arg == "--annotate-bed" and i + 1 < len(config.extra_args):
                if config.extra_args[i + 1] == "test_regions.bed":
                    new_extra_args.append("tests/fixtures/test_regions.bed")
                    i += 1  # Skip the next item
                else:
                    new_extra_args.append(config.extra_args[i + 1])
        config.extra_args = new_extra_args


if __name__ == "__main__":
    # Allow running as script for debugging
    pytest.main([__file__, "-v", "-s", "-m", "regression"])
