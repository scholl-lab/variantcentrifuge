"""Test that compares new pipeline outputs against baseline outputs."""

import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pytest

from variantcentrifuge.cli import parse_args
from variantcentrifuge.config import load_config
from variantcentrifuge.pipeline import run_pipeline

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


class TestBaselineComparison:
    """Compare new pipeline outputs against pre-generated baseline outputs."""

    @pytest.fixture
    def baseline_dir(self):
        """Path to baseline outputs directory."""
        return Path(__file__).parent / "baseline_outputs"

    @pytest.fixture
    def baseline_configs(self, baseline_dir):
        """Load all baseline configurations."""
        configs_dir = baseline_dir / "configs"
        if not configs_dir.exists():
            pytest.skip(f"Baseline configs directory not found: {configs_dir}")

        configs = {}
        for config_file in configs_dir.glob("*.json"):
            with open(config_file) as f:
                config = json.load(f)
                configs[config["name"]] = config

        return configs

    def compare_tsv_outputs(
        self, baseline_file: Path, new_file: Path, ignore_columns: Optional[List[str]] = None
    ) -> Tuple[bool, str]:
        """Compare two TSV files for equality."""
        # Read files
        baseline_df = pd.read_csv(baseline_file, sep="\t", dtype=str)
        new_df = pd.read_csv(new_file, sep="\t", dtype=str)

        # Fill NaN with empty string for comparison
        baseline_df = baseline_df.fillna("")
        new_df = new_df.fillna("")

        # Check shape
        if baseline_df.shape != new_df.shape:
            return False, f"Shape mismatch: baseline={baseline_df.shape}, new={new_df.shape}"

        # Check columns
        if set(baseline_df.columns) != set(new_df.columns):
            missing_in_new = set(baseline_df.columns) - set(new_df.columns)
            missing_in_baseline = set(new_df.columns) - set(baseline_df.columns)
            return (
                False,
                (
                    f"Column mismatch: missing_in_new={missing_in_new}, "
                    f"missing_in_baseline={missing_in_baseline}"
                ),
            )

        # Ignore specified columns
        if ignore_columns:
            cols_to_compare = [col for col in baseline_df.columns if col not in ignore_columns]
            baseline_df = baseline_df[cols_to_compare]
            new_df = new_df[cols_to_compare]

        # Sort both DataFrames by key columns for consistent comparison
        sort_cols = []
        if "CHROM" in baseline_df.columns and "POS" in baseline_df.columns:
            sort_cols = ["CHROM", "POS"]
            if "REF" in baseline_df.columns and "ALT" in baseline_df.columns:
                sort_cols.extend(["REF", "ALT"])

        if sort_cols:
            baseline_df = baseline_df.sort_values(by=sort_cols).reset_index(drop=True)
            new_df = new_df.sort_values(by=sort_cols).reset_index(drop=True)

        # Compare DataFrames
        try:
            pd.testing.assert_frame_equal(baseline_df, new_df, check_dtype=False)
            return True, "Files match"
        except AssertionError as e:
            return False, str(e)

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "config_name",
        [
            "single_gene_basic",
            "single_gene_rare_coding",
            "multi_gene_file",
            "gene_with_scoring",
            "gene_with_inheritance",
            "gene_with_annotations",
        ],
    )
    def test_baseline_comparison(self, config_name, baseline_dir, baseline_configs, tmp_path):
        """Test that new pipeline produces same output as baseline."""
        if config_name not in baseline_configs:
            pytest.skip(f"Baseline config not found: {config_name}")

        config = baseline_configs[config_name]

        # Get baseline output file
        baseline_output = baseline_dir / "outputs" / f"{config_name}.tsv"
        if not baseline_output.exists():
            pytest.skip(f"Baseline output not found: {baseline_output}")

        # Read command file to get exact arguments
        command_file = baseline_dir / "configs" / f"{config_name}_command.txt"
        if command_file.exists():
            command = command_file.read_text().strip().split()
            # Remove 'variantcentrifuge' from beginning
            if command[0] == "variantcentrifuge":
                command = command[1:]
        else:
            # Reconstruct command from config
            command = self._reconstruct_command(config, baseline_dir)

        # Update output directory to temp directory
        new_command = []
        skip_next = False
        for i, arg in enumerate(command):
            if skip_next:
                skip_next = False
                continue
            if arg == "--output-dir":
                new_command.extend(["--output-dir", str(tmp_path)])
                skip_next = True
            elif arg == "--log-file":
                new_command.extend(["--log-file", str(tmp_path / "test.log")])
                skip_next = True
            else:
                new_command.append(arg)

        # Add flag to use new pipeline
        new_command.append("--use-new-pipeline")

        # Parse args and run pipeline
        args = parse_args(new_command)
        cfg = load_config(args.config) if args.config else {}
        cfg.update(vars(args))
        cfg["use_new_pipeline_architecture"] = True

        # Run new pipeline
        import datetime

        start_time = datetime.datetime.now()
        run_pipeline(args, cfg, start_time)

        # Find output file
        output_filename = f"{config_name}.tsv"
        new_output = tmp_path / output_filename

        if not new_output.exists():
            # Try to find it in intermediate directory
            intermediate_dir = tmp_path / "intermediate"
            if intermediate_dir.exists():
                for file in intermediate_dir.glob("*.tsv"):
                    if config_name in file.name:
                        new_output = file
                        break

        assert new_output.exists(), f"New pipeline output not found: {new_output}"

        # Compare outputs
        # Ignore columns that might differ between runs (e.g., timestamps, temp file paths)
        ignore_columns = ["Pipeline_Run_ID", "Analysis_Timestamp", "Temp_File_Path"]

        match, message = self.compare_tsv_outputs(
            baseline_output, new_output, ignore_columns=ignore_columns
        )

        assert match, f"Output mismatch for {config_name}: {message}"

    def _reconstruct_command(self, config: Dict, baseline_dir: Path) -> List[str]:
        """Reconstruct command from config."""
        command = []

        # VCF file
        vcf_file = (
            baseline_dir.parent.parent
            / "fixtures"
            / "test_regression_variants.GRCh37.annotated.vcf.gz"
        )
        command.extend(["--vcf-file", str(vcf_file)])

        # Gene selection
        if config.get("gene_name"):
            command.extend(["--gene-name", config["gene_name"]])
        elif config.get("gene_file"):
            gene_file = baseline_dir / config["gene_file"]
            command.extend(["--gene-file", str(gene_file)])

        # Extra args
        if config.get("extra_args"):
            # Process extra args to fix file paths
            extra_args = []
            skip_next = False
            for i, arg in enumerate(config["extra_args"]):
                if skip_next:
                    skip_next = False
                    continue

                if arg in ["--ped", "--annotate-bed", "--scoring-config-path"]:
                    extra_args.append(arg)
                    if i + 1 < len(config["extra_args"]):
                        file_arg = config["extra_args"][i + 1]
                        if file_arg.startswith("scoring/"):
                            # Use actual scoring path
                            extra_args.append(file_arg)
                        else:
                            # Use baseline directory file
                            extra_args.append(str(baseline_dir / file_arg))
                        skip_next = True
                else:
                    extra_args.append(arg)

            command.extend(extra_args)

        # Standard fields
        command.extend(
            [
                "--fields",
                (
                    "CHROM POS REF ALT ID FILTER QUAL AC AF ANN[0].GENE ANN[0].EFFECT "
                    "ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P CADD_phred "
                    "dbNSFP_gnomAD_exomes_AF ClinVar_CLNSIG GEN[*].GT"
                ),
            ]
        )

        return command
