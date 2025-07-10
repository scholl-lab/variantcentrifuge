"""
Setup stages for configuration and data loading.

This module contains stages that load configuration files, phenotype data,
scoring configurations, and other setup tasks that can run in parallel.
"""

import logging
import os
from typing import List, Optional

from ..config import load_config
from ..helpers import get_vcf_samples
from ..ped_reader import read_pedigree
from ..phenotype import load_phenotypes
from ..pipeline_core import PipelineContext, Stage
from ..scoring import read_scoring_config

logger = logging.getLogger(__name__)


def load_terms_from_file(file_path: Optional[str], logger: logging.Logger) -> List[str]:
    """
    Load terms (HPO terms, sample IDs, etc.) from a file, one per line.

    Parameters
    ----------
    file_path : str or None
        Path to a file containing one term per line.
    logger : logging.Logger
        Logger instance for error logging.

    Returns
    -------
    list of str
        A list of terms loaded from the file.

    Raises
    ------
    SystemExit
        If the file is missing or empty and a file_path was specified.
    """
    terms: List[str] = []
    if file_path:
        if not os.path.exists(file_path):
            logger.error(f"Required file not found: {file_path}")
            raise FileNotFoundError(f"Required file not found: {file_path}")
        with open(file_path, "r", encoding="utf-8") as f:
            found_any = False
            for line in f:
                t = line.strip()
                if t:
                    found_any = True
                    terms.append(t)
            if not found_any:
                logger.error(f"File {file_path} is empty or invalid.")
                raise ValueError(f"File {file_path} is empty or invalid.")
    return terms


class ConfigurationLoadingStage(Stage):
    """Load and merge configuration from file and CLI arguments."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "configuration_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load configuration and merge with CLI arguments"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load configuration and merge with CLI arguments."""
        logger.info("Loading configuration...")

        # Check if config was already passed as a dict (from main pipeline)
        if hasattr(context.args, "config") and isinstance(context.args.config, dict):
            # Config already fully merged by cli.py, just use it
            config = context.args.config
            logger.info("Using pre-merged configuration from main pipeline")
        else:
            # Load base configuration
            if hasattr(context.args, "config") and context.args.config:
                config = load_config(context.args.config)
                logger.info(f"Loaded configuration from {context.args.config}")
            else:
                config = load_config()  # Load default config
                logger.info("Loaded default configuration")

            # Merge CLI arguments (CLI takes precedence)
            cli_args = vars(context.args)
            for key, value in cli_args.items():
                if value is not None and key not in ["config", "fields"]:
                    # Skip internal arguments and fields (handled as fields_to_extract)
                    if key.startswith("_"):
                        continue
                    config[key] = value

            # Ensure phenotype arguments are properly set even if they come as None
            # This fixes the issue where --phenotype-file args don't make it to config
            phenotype_args = ["phenotype_file", "phenotype_sample_column", "phenotype_value_column"]
            for arg in phenotype_args:
                if hasattr(context.args, arg):
                    config[arg] = getattr(context.args, arg)

        # Merge any existing config from context (e.g., from tests)
        if context.config:
            for key, value in context.config.items():
                if key not in config:
                    config[key] = value

        # Set pipeline version
        config["pipeline_version"] = config.get("pipeline_version", "1.0.0")

        # Convert fields_to_extract to extract for downstream stages
        if "extract" not in config and "fields_to_extract" in config:
            config["extract"] = config["fields_to_extract"].split()
            logger.debug(
                f"Converted fields_to_extract to extract list with {len(config['extract'])} fields"
            )

        # Parse phenotype arguments from CLI (matching original pipeline logic)
        # This fixes the bug where phenotype-based case/control assignment doesn't work
        case_hpo_terms = []
        control_hpo_terms = []

        # Parse case phenotypes
        if config.get("case_phenotypes"):
            if isinstance(config["case_phenotypes"], str):
                # Parse comma-separated string
                case_hpo_terms = [
                    t.strip() for t in config["case_phenotypes"].split(",") if t.strip()
                ]
                logger.debug(f"Parsed {len(case_hpo_terms)} case phenotype terms from CLI")
            else:
                # Already a list
                case_hpo_terms = config["case_phenotypes"]

        # Load additional case phenotypes from file
        if config.get("case_phenotypes_file"):
            file_terms = load_terms_from_file(config["case_phenotypes_file"], logger)
            case_hpo_terms.extend(file_terms)
            logger.debug(f"Loaded {len(file_terms)} additional case phenotype terms from file")

        # Parse control phenotypes
        if config.get("control_phenotypes"):
            if isinstance(config["control_phenotypes"], str):
                # Parse comma-separated string
                control_hpo_terms = [
                    t.strip() for t in config["control_phenotypes"].split(",") if t.strip()
                ]
                logger.debug(f"Parsed {len(control_hpo_terms)} control phenotype terms from CLI")
            else:
                # Already a list
                control_hpo_terms = config["control_phenotypes"]

        # Load additional control phenotypes from file
        if config.get("control_phenotypes_file"):
            file_terms = load_terms_from_file(config["control_phenotypes_file"], logger)
            control_hpo_terms.extend(file_terms)
            logger.debug(f"Loaded {len(file_terms)} additional control phenotype terms from file")

        # Update config with parsed lists
        config["case_phenotypes"] = case_hpo_terms
        config["control_phenotypes"] = control_hpo_terms

        # Also parse case/control samples in the same way
        case_samples = []
        control_samples = []

        # Parse case samples
        if config.get("case_samples"):
            if isinstance(config["case_samples"], str):
                # Parse comma-separated string
                case_samples = [s.strip() for s in config["case_samples"].split(",") if s.strip()]
                logger.debug(f"Parsed {len(case_samples)} case samples from CLI")
            else:
                # Already a list
                case_samples = config["case_samples"]

        # Parse control samples
        if config.get("control_samples"):
            if isinstance(config["control_samples"], str):
                # Parse comma-separated string
                control_samples = [
                    s.strip() for s in config["control_samples"].split(",") if s.strip()
                ]
                logger.debug(f"Parsed {len(control_samples)} control samples from CLI")
            else:
                # Already a list
                control_samples = config["control_samples"]

        # Update config with parsed lists
        config["case_samples"] = case_samples
        config["control_samples"] = control_samples

        # Update context
        context.config = config

        logger.debug(f"Configuration loaded with {len(config)} settings")
        return context


class PhenotypeLoadingStage(Stage):
    """Load phenotype data from file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "phenotype_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load phenotype data"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load phenotype data if provided."""
        phenotype_file = context.config.get("phenotype_file")
        phenotype_sample_column = context.config.get("phenotype_sample_column")
        phenotype_value_column = context.config.get("phenotype_value_column")

        logger.debug(f"Phenotype file: {phenotype_file}")
        logger.debug(f"Phenotype sample column: {phenotype_sample_column}")
        logger.debug(f"Phenotype value column: {phenotype_value_column}")
        logger.debug(f"Config keys: {list(context.config.keys())}")

        if phenotype_file and phenotype_sample_column and phenotype_value_column:
            logger.info(f"Loading phenotypes from {phenotype_file}")

            phenotypes = load_phenotypes(
                phenotype_file, phenotype_sample_column, phenotype_value_column
            )

            if not phenotypes:
                raise ValueError(
                    f"No phenotype data loaded from {phenotype_file}. "
                    "Check file formatting and column names."
                )

            context.phenotype_data = phenotypes
            context.config["phenotypes"] = phenotypes  # For compatibility
            logger.info(f"Loaded phenotypes for {len(phenotypes)} samples")
        else:
            logger.debug("No phenotype file specified, skipping phenotype loading")

        return context


class ScoringConfigLoadingStage(Stage):
    """Load variant scoring configuration."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "scoring_config_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load scoring configuration"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load scoring configuration if provided."""
        scoring_config_path = context.config.get("scoring_config_path")

        if scoring_config_path:
            logger.info(f"Loading scoring configuration from {scoring_config_path}")

            try:
                scoring_config = read_scoring_config(scoring_config_path)
                context.scoring_config = scoring_config
                logger.info(
                    f"Successfully loaded scoring configuration with "
                    f"{len(scoring_config.get('formulas', {}))} formulas"
                )
            except Exception as e:
                raise ValueError(f"Failed to load scoring configuration: {e}")
        else:
            logger.debug("No scoring configuration specified")

        return context


class PedigreeLoadingStage(Stage):
    """Load pedigree data from PED file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "pedigree_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load pedigree data"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load pedigree data if provided."""
        ped_file = context.config.get("ped_file") or context.config.get("ped")

        if ped_file:
            logger.info(f"Loading pedigree from {ped_file}")

            try:
                pedigree_data = read_pedigree(ped_file)
                context.pedigree_data = pedigree_data

                # Log pedigree summary
                n_samples = len(pedigree_data)
                n_affected = sum(
                    1 for s in pedigree_data.values() if s.get("affected_status") == "2"
                )
                logger.info(f"Loaded pedigree with {n_samples} samples ({n_affected} affected)")
            except Exception as e:
                raise ValueError(f"Failed to load pedigree file: {e}")
        else:
            logger.debug("No pedigree file specified")

        return context


class AnnotationConfigLoadingStage(Stage):
    """Load custom annotation configurations (BED files, gene lists, JSON)."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "annotation_config_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load annotation configurations"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load all annotation configurations."""
        annotations = {}

        # Load BED file annotations
        bed_files = context.config.get("annotate_bed", [])
        if bed_files:
            if isinstance(bed_files, str):
                bed_files = [bed_files]
            annotations["bed_files"] = bed_files
            logger.info(f"Configured {len(bed_files)} BED file annotations")

        # Load gene list annotations
        gene_lists = context.config.get("annotate_gene_list", [])
        if gene_lists:
            if isinstance(gene_lists, str):
                gene_lists = [gene_lists]
            annotations["gene_lists"] = gene_lists
            logger.info(f"Configured {len(gene_lists)} gene list annotations")

        # Load JSON gene annotations
        json_genes = context.config.get("annotate_json_genes")
        json_mapping = context.config.get("json_gene_mapping")
        if json_genes and json_mapping:
            annotations["json_genes"] = json_genes
            annotations["json_mapping"] = json_mapping
            logger.info(f"Configured JSON gene annotations from {json_genes}")

        if annotations:
            context.annotation_configs = annotations

        return context


class SampleConfigLoadingStage(Stage):
    """Load sample configurations and VCF sample names."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "sample_config_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load sample configurations"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    @property
    def dependencies(self) -> set:
        """Return the set of stage names this stage depends on."""
        # Needs config to know VCF file path
        return {"configuration_loading"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load VCF samples and any sample-specific configurations."""
        vcf_file = context.config.get("vcf_file")

        if not vcf_file:
            raise ValueError("No VCF file specified in configuration")

        # Get VCF samples
        logger.info(f"Reading sample names from {vcf_file}")
        vcf_samples = get_vcf_samples(vcf_file)

        # Apply sample substring removal if configured
        remove_substring = context.config.get("remove_sample_substring")
        if remove_substring and remove_substring.strip():
            logger.info(f"Removing substring '{remove_substring}' from sample names")
            original_samples = vcf_samples[:]  # Create a copy for logging
            vcf_samples = [s.replace(remove_substring, "") for s in vcf_samples]

            # Log changes for debugging
            changes = []
            for orig, new in zip(original_samples, vcf_samples):
                if orig != new:
                    changes.append(f"{orig} -> {new}")

            if changes:
                logger.info(f"Applied substring removal to {len(changes)} samples")
                logger.debug(f"Sample name changes: {changes[:5]}...")  # Log first 5 changes
            else:
                logger.info(f"No samples contained substring '{remove_substring}'")

        context.vcf_samples = vcf_samples
        logger.info(f"Found {len(vcf_samples)} samples in VCF")

        # CRITICAL: Set sample_list in config for downstream analysis
        # This is required by analyze_variants.py for case/control assignment
        context.config["sample_list"] = ",".join(vcf_samples)
        logger.debug(f"Set sample_list config with {len(vcf_samples)} processed sample names")

        # Load case/control sample lists if provided via files
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        logger.debug(f"Case samples file from config: {case_samples_file}")
        logger.debug(f"Control samples file from config: {control_samples_file}")

        if case_samples_file:
            logger.debug(f"Loading case samples from file: {case_samples_file}")
            with open(case_samples_file, "r") as f:
                case_samples = [line.strip() for line in f if line.strip()]
            context.config["case_samples"] = case_samples
            logger.info(f"Loaded {len(case_samples)} case samples from file")
            logger.debug(f"Case samples: {case_samples[:10]}...")  # Log first 10 for debugging
        else:
            logger.debug("No case samples file provided")

        if control_samples_file:
            logger.debug(f"Loading control samples from file: {control_samples_file}")
            with open(control_samples_file, "r") as f:
                control_samples = [line.strip() for line in f if line.strip()]
            context.config["control_samples"] = control_samples
            logger.info(f"Loaded {len(control_samples)} control samples from file")
            logger.debug(
                f"Control samples: {control_samples[:10]}..."
            )  # Log first 10 for debugging
        else:
            logger.debug("No control samples file provided")

        # CRUCIAL: Handle automatic control assignment when only case samples are provided
        case_samples = context.config.get("case_samples", [])
        control_samples = context.config.get("control_samples", [])

        if case_samples and not control_samples:
            # Automatically assign all non-case samples as controls
            case_samples_set = set(case_samples)
            vcf_samples_set = set(vcf_samples)
            auto_control_samples = list(vcf_samples_set - case_samples_set)
            context.config["control_samples"] = auto_control_samples
            logger.info(
                f"Automatically assigned {len(auto_control_samples)} control samples (all non-case samples)"
            )
            logger.debug(f"Auto control samples: {auto_control_samples[:10]}...")
        elif control_samples and not case_samples:
            # Automatically assign all non-control samples as cases
            control_samples_set = set(control_samples)
            vcf_samples_set = set(vcf_samples)
            auto_case_samples = list(vcf_samples_set - control_samples_set)
            context.config["case_samples"] = auto_case_samples
            logger.info(
                f"Automatically assigned {len(auto_case_samples)} case samples (all non-control samples)"
            )
            logger.debug(f"Auto case samples: {auto_case_samples[:10]}...")

        # If no explicit sample files provided but phenotype information exists,
        # we need to derive case/control assignments from phenotypes after phenotype loading
        # This will be handled in a separate stage after phenotype data is loaded
        if not case_samples_file and not control_samples_file:
            # Check if we have phenotype-based assignment criteria
            case_phenotypes = context.config.get("case_phenotypes", [])
            control_phenotypes = context.config.get("control_phenotypes", [])
            phenotype_file = context.config.get("phenotype_file")

            if (case_phenotypes or control_phenotypes) and phenotype_file:
                logger.info(
                    "Phenotype-based case/control assignment will be performed "
                    "after phenotype loading"
                )
            elif not (context.config.get("case_samples") or context.config.get("control_samples")):
                logger.debug(
                    "No explicit case/control assignment specified - "
                    "will use all samples as controls"
                )

        return context


class PhenotypeCaseControlAssignmentStage(Stage):
    """Assign case/control samples based on phenotype data when no explicit files provided."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "phenotype_case_control_assignment"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Assign case/control samples based on phenotype data"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    @property
    def dependencies(self) -> set:
        """Return the set of stage names this stage depends on."""
        return {"phenotype_loading", "sample_config_loading"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Assign case/control samples based on phenotype data if no explicit assignments exist."""

        # Check if case/control samples are already explicitly assigned
        existing_case_samples = context.config.get("case_samples", [])
        existing_control_samples = context.config.get("control_samples", [])
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        # Skip if explicit assignments exist
        if (
            existing_case_samples
            or existing_control_samples
            or case_samples_file
            or control_samples_file
        ):
            logger.debug(
                "Explicit case/control assignments exist, skipping phenotype-based assignment"
            )
            return context

        # Get phenotype criteria - these should already be parsed as lists
        case_phenotypes = context.config.get("case_phenotypes", [])
        control_phenotypes = context.config.get("control_phenotypes", [])

        # If they're still strings, parse them here
        if isinstance(case_phenotypes, str):
            case_phenotypes = [t.strip() for t in case_phenotypes.split(",") if t.strip()]
        if isinstance(control_phenotypes, str):
            control_phenotypes = [t.strip() for t in control_phenotypes.split(",") if t.strip()]

        # Skip if no phenotype criteria
        if not case_phenotypes and not control_phenotypes:
            logger.debug("No phenotype criteria specified, skipping phenotype-based assignment")
            return context

        # Get VCF samples
        vcf_samples = getattr(context, "vcf_samples", [])
        if not vcf_samples:
            logger.warning("No VCF samples found, cannot perform phenotype-based assignment")
            return context

        # Get phenotype file parameters
        phenotype_file = context.config.get("phenotype_file")
        phenotype_sample_column = context.config.get("phenotype_sample_column")
        phenotype_value_column = context.config.get("phenotype_value_column")

        logger.debug(f"All config keys: {sorted(context.config.keys())}")
        logger.debug(f"Phenotype file: {phenotype_file}")
        logger.debug(f"Sample column: {phenotype_sample_column}")
        logger.debug(f"Value column: {phenotype_value_column}")
        logger.debug(f"Case phenotypes: {case_phenotypes}")
        logger.debug(f"Control phenotypes: {control_phenotypes}")

        if not phenotype_file or not phenotype_sample_column or not phenotype_value_column:
            logger.warning(
                "Missing phenotype file parameters, cannot perform phenotype-based assignment"
            )
            logger.warning(
                f"phenotype_file={phenotype_file}, sample_column={phenotype_sample_column}, value_column={phenotype_value_column}"
            )
            return context

        # Simple logic: load CSV and filter by phenotype terms to get sample lists
        try:
            remove_substring = context.config.get("remove_sample_substring", "")
            case_samples, control_samples = self._compute_samples_from_phenotype_file(
                phenotype_file,
                phenotype_sample_column,
                phenotype_value_column,
                case_phenotypes,
                control_phenotypes,
                vcf_samples,
                remove_substring,
            )

            # Update configuration with classified samples
            context.config["case_samples"] = case_samples
            context.config["control_samples"] = control_samples

            logger.info(
                f"Phenotype-based assignment complete: {len(case_samples)} cases, "
                f"{len(control_samples)} controls"
            )

            if case_phenotypes and len(case_samples) == 0:
                logger.warning("No samples matched the specified case phenotypes")
            if control_phenotypes and len(control_samples) == 0:
                logger.warning("No samples matched the specified control phenotypes")

            # Log sample assignment summary
            logger.debug(f"Case samples (first 10): {case_samples[:10]}")
            logger.debug(f"Control samples (first 10): {control_samples[:10]}")

        except Exception as e:
            logger.error(f"Failed to perform phenotype-based assignment: {e}")
            return context

        return context

    def _compute_samples_from_phenotype_file(
        self,
        phenotype_file,
        sample_column,
        value_column,
        case_phenotypes,
        control_phenotypes,
        vcf_samples,
        remove_substring="",
    ):
        """
        Simple function: load CSV, filter by phenotype terms, get sample lists.

        As user requested: filter the phenotype-value-column by the phenotype list,
        get the phenotype-sample-column by group and deduplicate -> sample list.
        """
        import pandas as pd

        logger.info(f"Loading phenotype file: {phenotype_file}")

        # Determine delimiter
        if phenotype_file.endswith(".csv"):
            df = pd.read_csv(phenotype_file)
        else:
            df = pd.read_csv(phenotype_file, sep="\t")

        logger.info(f"Loaded {len(df)} rows from phenotype file")
        logger.debug(f"Columns: {list(df.columns)}")

        # Check required columns exist
        if sample_column not in df.columns:
            raise ValueError(f"Sample column '{sample_column}' not found in phenotype file")
        if value_column not in df.columns:
            raise ValueError(f"Value column '{value_column}' not found in phenotype file")

        # Convert to sets for efficient matching
        case_phenotype_set = set(case_phenotypes) if case_phenotypes else set()
        control_phenotype_set = set(control_phenotypes) if control_phenotypes else set()
        # Convert VCF samples to strings for consistent comparison
        vcf_samples_set = set(str(s) for s in vcf_samples)

        logger.info(f"Looking for case phenotypes: {case_phenotype_set}")
        logger.info(f"Looking for control phenotypes: {control_phenotype_set}")

        # Debug: show sample matching process
        logger.debug(f"First 5 phenotype file samples: {df[sample_column].unique()[:5].tolist()}")
        logger.debug(f"First 5 VCF samples: {vcf_samples[:5]}")

        # Convert sample IDs to strings for consistent comparison
        # VCF samples are strings, but phenotype file samples might be integers
        vcf_samples_str = set(str(s) for s in vcf_samples)
        pheno_samples_raw = set(str(s) for s in df[sample_column].unique())

        # Apply sample substring removal to phenotype samples if configured
        # This is necessary because VCF samples have already been processed
        if remove_substring and remove_substring.strip():
            logger.debug(
                f"Applying substring removal '{remove_substring}' to phenotype samples for matching"
            )
            pheno_samples_str = set(s.replace(remove_substring, "") for s in pheno_samples_raw)
        else:
            pheno_samples_str = pheno_samples_raw

        # Check if there are any VCF samples in the phenotype file
        vcf_in_pheno = vcf_samples_str & pheno_samples_str
        logger.debug(f"{len(vcf_in_pheno)} VCF samples found in phenotype file")
        if vcf_in_pheno:
            logger.debug(f"Example matching samples: {list(vcf_in_pheno)[:5]}")

        # Show some phenotype values being searched for
        unique_phenotypes = set(df[value_column].unique())
        logger.debug(f"Found {len(unique_phenotypes)} unique phenotypes in file")
        logger.debug(f"Example phenotypes: {list(unique_phenotypes)[:10]}")

        # Check for exact matches with case phenotypes
        case_matches = case_phenotype_set & unique_phenotypes
        logger.debug(f"{len(case_matches)} case phenotypes found in file: {case_matches}")

        # Find samples that match case phenotypes
        case_samples = set()
        if case_phenotype_set:
            case_mask = df[value_column].isin(case_phenotype_set)
            case_rows = df[case_mask]
            # Convert to strings for consistent comparison
            case_samples_raw = set(str(s) for s in case_rows[sample_column].unique())
            # Apply substring removal if configured
            if remove_substring and remove_substring.strip():
                case_samples = set(s.replace(remove_substring, "") for s in case_samples_raw)
            else:
                case_samples = case_samples_raw
            logger.info(f"Found {len(case_samples)} samples matching case phenotypes")

        # Find samples that match control phenotypes
        control_samples = set()
        if control_phenotype_set:
            control_mask = df[value_column].isin(control_phenotype_set)
            control_rows = df[control_mask]
            # Convert to strings for consistent comparison
            control_samples_raw = set(str(s) for s in control_rows[sample_column].unique())
            # Apply substring removal if configured
            if remove_substring and remove_substring.strip():
                control_samples = set(s.replace(remove_substring, "") for s in control_samples_raw)
            else:
                control_samples = control_samples_raw
            logger.info(f"Found {len(control_samples)} samples matching control phenotypes")

        # Apply assignment logic
        if case_phenotype_set and not control_phenotype_set:
            # Only case phenotypes specified - all non-matching VCF samples are controls
            final_case_samples = case_samples & vcf_samples_set
            final_control_samples = vcf_samples_set - final_case_samples
        elif control_phenotype_set and not case_phenotype_set:
            # Only control phenotypes specified - all non-matching VCF samples are cases
            final_control_samples = control_samples & vcf_samples_set
            final_case_samples = vcf_samples_set - final_control_samples
        else:
            # Both specified - exact matching only
            final_case_samples = case_samples & vcf_samples_set
            final_control_samples = control_samples & vcf_samples_set

        logger.info(
            f"Final assignment: {len(final_case_samples)} cases, {len(final_control_samples)} controls"
        )

        return list(final_case_samples), list(final_control_samples)

    def _classify_samples_by_phenotypes(
        self, vcf_samples, phenotype_data, case_phenotypes, control_phenotypes
    ):
        """Simple modular function to classify samples based on phenotype data."""
        logger.info(
            f"Performing phenotype-based case/control assignment for {len(vcf_samples)} VCF samples"
        )
        logger.info(
            f"Using {len(case_phenotypes)} case phenotypes and "
            f"{len(control_phenotypes)} control phenotypes"
        )
        logger.debug(f"Case phenotypes: {case_phenotypes}")
        logger.debug(f"Control phenotypes: {control_phenotypes}")

        classified_cases = set()
        classified_controls = set()

        for sample in vcf_samples:
            sample_phenotypes = phenotype_data.get(sample, set())

            # Check if sample matches case or control phenotypes
            matches_case = (
                any(pheno in sample_phenotypes for pheno in case_phenotypes)
                if case_phenotypes
                else False
            )
            matches_control = (
                any(pheno in sample_phenotypes for pheno in control_phenotypes)
                if control_phenotypes
                else False
            )

            if case_phenotypes and not control_phenotypes:
                # Only case phenotypes specified - all non-matching samples are controls
                if matches_case:
                    classified_cases.add(sample)
                else:
                    classified_controls.add(sample)
            elif control_phenotypes and not case_phenotypes:
                # Only control phenotypes specified - all non-matching samples are cases
                if matches_control:
                    classified_controls.add(sample)
                else:
                    classified_cases.add(sample)
            else:
                # Both case and control phenotypes specified
                if matches_case and not matches_control:
                    classified_cases.add(sample)
                elif matches_control and not matches_case:
                    classified_controls.add(sample)
                # Samples matching both or neither are not classified

        return classified_cases, classified_controls
