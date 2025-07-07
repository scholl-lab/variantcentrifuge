#!/usr/bin/env python
"""Generate baseline outputs using the old pipeline for regression testing.

This script runs the old pipeline with various configurations to create
baseline outputs that the new pipeline should match exactly.
"""

import argparse
import json
import logging
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from tests.regression.test_regression_suite import REGRESSION_TESTS, RegressionTestConfig

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


class BaselineGenerator:
    """Generates baseline outputs for regression testing."""

    def __init__(self, output_dir: Path, vcf_file: str):
        self.output_dir = output_dir
        self.vcf_file = vcf_file
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        self.configs_dir = self.output_dir / "configs"
        self.outputs_dir = self.output_dir / "outputs"
        self.logs_dir = self.output_dir / "logs"

        for dir_path in [self.configs_dir, self.outputs_dir, self.logs_dir]:
            dir_path.mkdir(exist_ok=True)

    def generate_all_baselines(self) -> Dict[str, bool]:
        """Generate baseline outputs for all test configurations.

        Returns dict mapping test name to success status.
        """
        results = {}

        logger.info(f"Generating baselines for {len(REGRESSION_TESTS)} test cases")
        logger.info(f"Using VCF: {self.vcf_file}")
        logger.info(f"Output directory: {self.output_dir}")

        for config in REGRESSION_TESTS:
            logger.info(f"\nProcessing: {config.name}")
            logger.info(f"Description: {config.description}")

            try:
                self._generate_baseline(config)
                results[config.name] = True
                logger.info(f"✓ Successfully generated baseline for {config.name}")
            except Exception as e:
                results[config.name] = False
                logger.error(f"✗ Failed to generate baseline for {config.name}: {e}")

        # Write summary
        self._write_summary(results)

        return results

    def _generate_baseline(self, config: RegressionTestConfig) -> None:
        """Generate baseline for a single test configuration."""
        # Save configuration
        config_file = self.configs_dir / f"{config.name}.json"
        with open(config_file, "w") as f:
            json.dump(
                {
                    "name": config.name,
                    "description": config.description,
                    "gene_name": config.gene_name,
                    "gene_file": config.gene_file,
                    "preset": config.preset,
                    "extra_args": config.extra_args,
                },
                f,
                indent=2,
            )

        # Prepare output files
        output_filename = f"{config.name}.tsv"
        log_file = self.logs_dir / f"{config.name}.log"
        
        # Ensure directories exist
        self.outputs_dir.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        # Build command
        cmd = [
            "variantcentrifuge",
            "--vcf-file",
            self.vcf_file,
            "--output-file",
            output_filename,
            "--output-dir",
            str(self.outputs_dir),
        ]

        if config.gene_name:
            cmd.extend(["--gene-name", config.gene_name])
        elif config.gene_file:
            # Create gene file if needed
            gene_file_path = self.output_dir / config.gene_file
            if not gene_file_path.exists():
                self._create_test_gene_file(gene_file_path)
            cmd.extend(["--gene-file", str(gene_file_path)])

        # Skip preset handling since we're using explicit filters in extra_args

        # Process extra args
        processed_args = []
        for i, arg in enumerate(config.extra_args):
            if arg in ["--ped", "--annotate-bed", "--scoring-config-path"]:
                # Next arg is a file path
                if i + 1 < len(config.extra_args):
                    file_arg = config.extra_args[i + 1]
                    file_path = self._prepare_test_file(file_arg)
                    processed_args.append(arg)
                    processed_args.append(str(file_path))
                else:
                    processed_args.append(arg)
            elif config.extra_args[i - 1] not in [
                "--ped",
                "--annotate-bed",
                "--scoring-config-path",
            ]:
                # Not a file path argument
                processed_args.append(arg)

        cmd.extend(processed_args)

        # Add standard args
        cmd.extend(
            [
                "--threads",
                "1",
                "--log-file",
                str(log_file),
                "--fields",
                "CHROM POS REF ALT ID FILTER QUAL AC AF ANN[0].GENE ANN[0].EFFECT ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P CADD_phred dbNSFP_gnomAD_exomes_AF ClinVar_CLNSIG GEN[*].GT",
            ]
        )

        # Run command
        logger.debug(f"Running: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            logger.error(f"Command failed: {' '.join(cmd)}")
            logger.error(f"stdout: {result.stdout}")
            logger.error(f"stderr: {result.stderr}")
            raise RuntimeError(f"Pipeline failed with return code {result.returncode}")

        # Save command for reproducibility
        cmd_file = self.configs_dir / f"{config.name}_command.txt"
        cmd_file.write_text(" ".join(cmd))

        # Check outputs
        expected_output = self.outputs_dir / output_filename
        if not expected_output.exists():
            raise RuntimeError(f"Expected output file not created: {expected_output}")

        # Generate statistics if not disabled
        if "--no-stats" not in cmd:
            stats_file = self.outputs_dir / f"{config.name}_stats.json"
            if stats_file.exists():
                logger.debug(f"Statistics file created: {stats_file}")

    def _create_test_gene_file(self, path: Path) -> None:
        """Create a test gene file."""
        genes = [
            "BRCA1",
            "BRCA2",
            "TP53",
            "EGFR",
            "KRAS",
            "BRAF",
            "PIK3CA",
            "PTEN",
            "APC",
            "MLH1",
        ]
        path.write_text("\n".join(genes))
        logger.debug(f"Created gene file: {path}")

    def _prepare_test_file(self, file_arg: str) -> Path:
        """Prepare test files referenced in arguments."""
        file_path = self.output_dir / file_arg

        if file_arg.endswith(".ped"):
            # Create test PED file
            if not file_path.exists():
                file_path.write_text(
                    "FAM001\tCHILD\tFATHER\tMOTHER\t1\t2\n"
                    "FAM001\tFATHER\t0\t0\t1\t1\n"
                    "FAM001\tMOTHER\t0\t0\t2\t1\n"
                )
                logger.debug(f"Created PED file: {file_path}")

        elif file_arg.endswith(".bed"):
            # Create test BED file
            if not file_path.exists():
                file_path.write_text(
                    "chr1\t1000000\t2000000\tHotspot1\n"
                    "chr2\t5000000\t6000000\tHotspot2\n"
                    "chr17\t41196312\t41277500\tBRCA1_region\n"
                )
                logger.debug(f"Created BED file: {file_path}")

        elif file_arg.startswith("scoring/"):
            # Handle scoring config path
            scoring_dir = Path("scoring") / file_arg.split("/")[1]
            if scoring_dir.exists():
                return scoring_dir
            else:
                logger.warning(f"Scoring config not found: {scoring_dir}")
                # Return dummy path
                return Path(file_arg)

        return file_path

    def _write_summary(self, results: Dict[str, bool]) -> None:
        """Write summary of baseline generation."""
        summary = {
            "generated_at": datetime.now().isoformat(),
            "vcf_file": self.vcf_file,
            "output_dir": str(self.output_dir),
            "total_tests": len(results),
            "successful": sum(1 for v in results.values() if v),
            "failed": sum(1 for v in results.values() if not v),
            "results": results,
        }

        summary_file = self.output_dir / "baseline_summary.json"
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2)

        logger.info(f"\nSummary written to: {summary_file}")
        logger.info(f"Total tests: {summary['total_tests']}")
        logger.info(f"Successful: {summary['successful']}")
        logger.info(f"Failed: {summary['failed']}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate baseline outputs for regression testing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "vcf_file",
        nargs="?",
        default="tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz",
        help="VCF file to use for testing (default: tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz)",
    )

    parser.add_argument(
        "--output-dir",
        default="baseline_outputs",
        help="Directory for baseline outputs (default: baseline_outputs)",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Check VCF file exists
    vcf_path = Path(args.vcf_file)
    if not vcf_path.exists():
        logger.error(f"VCF file not found: {vcf_path}")
        sys.exit(1)

    # Check for required tools
    required_tools = ["variantcentrifuge", "bcftools", "snpEff", "SnpSift"]
    missing_tools = []

    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)

    if missing_tools:
        logger.error(f"Required tools not found in PATH: {', '.join(missing_tools)}")
        logger.error("Please install missing tools and ensure they are in PATH")
        sys.exit(1)

    # Generate baselines
    output_dir = Path(args.output_dir)
    generator = BaselineGenerator(output_dir, str(vcf_path))

    results = generator.generate_all_baselines()

    # Exit with error if any tests failed
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
