"""Sample and pedigree loading stages for VariantCentrifuge."""

import logging
import os
import sys
from collections import Counter
from typing import Any

from ..pipeline_core import PipelineContext, Stage

logger = logging.getLogger(__name__)


def get_vcf_samples(vcf_file: str) -> list[str]:
    """Get sample names from VCF, respecting setup_stages patches."""
    setup_mod = sys.modules.get("variantcentrifuge.stages.setup_stages")
    if (
        setup_mod
        and hasattr(setup_mod, "get_vcf_samples")
        and setup_mod.get_vcf_samples is not get_vcf_samples
    ):
        return setup_mod.get_vcf_samples(vcf_file)  # type: ignore[no-any-return]
    from ..helpers import get_vcf_samples as _real_get_vcf_samples

    return _real_get_vcf_samples(vcf_file)


def read_pedigree(ped_file: str) -> dict[str, Any]:
    """Read pedigree file, respecting setup_stages patches."""
    setup_mod = sys.modules.get("variantcentrifuge.stages.setup_stages")
    if (
        setup_mod
        and hasattr(setup_mod, "read_pedigree")
        and setup_mod.read_pedigree is not read_pedigree
    ):
        return setup_mod.read_pedigree(ped_file)  # type: ignore[no-any-return]
    from ..ped_reader import read_pedigree as _real_read_pedigree

    return _real_read_pedigree(ped_file)


class PedigreeLoadingStage(Stage):
    """Load pedigree data from PED file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "pedigree_loading"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Pedigree data can be restored by the checkpoint skip handler."""
        return "restore"

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
        sample_sheet = context.config.get("sample_sheet") or getattr(
            context.args, "sample_sheet", None
        )

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
                raise ValueError(f"Failed to load pedigree file: {e}") from e
        elif sample_sheet:
            from ..sample_sheet import load_sample_sheet

            logger.info(f"Loading pedigree from sample sheet {sample_sheet}")
            sheet_data = load_sample_sheet(str(sample_sheet))
            context.pedigree_data = sheet_data.pedigree_data
        else:
            logger.debug("No pedigree file specified")

        return context

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        Since this is the pedigree loading stage, we need to restore pedigree data
        that would have been loaded during normal execution.
        """
        ped_file = context.config.get("ped_file") or context.config.get("ped")
        sample_sheet = context.config.get("sample_sheet") or getattr(
            context.args, "sample_sheet", None
        )

        if ped_file:
            logger.debug(f"Restoring pedigree from {ped_file} (checkpoint skip)")

            try:
                pedigree_data = read_pedigree(ped_file)
                context.pedigree_data = pedigree_data

                # Log pedigree summary
                n_samples = len(pedigree_data)
                n_affected = sum(
                    1 for s in pedigree_data.values() if s.get("affected_status") == "2"
                )
                logger.debug(
                    f"Restored pedigree with {n_samples} samples ({n_affected} affected) "
                    f"(checkpoint skip)"
                )
            except Exception as e:
                raise ValueError(f"Failed to restore pedigree file: {e}") from e
        elif sample_sheet:
            from ..sample_sheet import load_sample_sheet

            sheet_data = load_sample_sheet(str(sample_sheet))
            context.pedigree_data = sheet_data.pedigree_data
        else:
            logger.debug("No pedigree file specified during checkpoint skip")

        return context


def _require_unique_samples(samples: list[str], remove_substring: str | None) -> None:
    """Reject duplicate sample IDs (e.g. "A.bam" and "A" after removing ".bam").

    Downstream genotype, phenotype and cohort lookups are keyed by sample ID; duplicates
    would silently give one sample's genotypes to both.
    """
    dups = sorted(s for s, n in Counter(samples).items() if n > 1)
    if dups:
        how = f" after removing substring '{remove_substring}'" if remove_substring else ""
        raise ValueError(
            f"Sample IDs are not unique{how}: {', '.join(dups[:10])}. Choose a "
            "--remove-sample-substring that keeps every sample ID distinct."
        )


class SampleConfigLoadingStage(Stage):
    """Load sample configurations and VCF sample names."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "sample_config_loading"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Sample configuration can be restored by the checkpoint skip handler."""
        return "restore"

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

    @property
    def soft_dependencies(self) -> set[str]:
        """Run after pedigree, phenotype, and scoring config loading when present."""
        return {"pedigree_loading", "phenotype_loading", "scoring_config_loading"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load VCF samples and any sample-specific configurations."""
        vcf_file = context.config.get("vcf_file")
        from_tsv = context.config.get("from_tsv") or getattr(context.args, "from_tsv", None)

        if from_tsv and not vcf_file:
            from ..from_tsv import prepare_from_tsv_context

            prepare_from_tsv_context(context)
            vcf_samples = list(context.vcf_samples or [])
            logger.info(f"Loaded {len(vcf_samples)} sample(s) from TSV {from_tsv}")
        else:
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
                for orig, new in zip(original_samples, vcf_samples, strict=False):
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

        _require_unique_samples(vcf_samples, context.config.get("remove_sample_substring"))

        # Load case/control sample lists if provided via files
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        logger.debug(f"Case samples file from config: {case_samples_file}")
        logger.debug(f"Control samples file from config: {control_samples_file}")

        if case_samples_file:
            logger.debug(f"Loading case samples from file: {case_samples_file}")
            with open(case_samples_file) as f:
                case_samples = [line.strip() for line in f if line.strip()]
            context.config["case_samples"] = case_samples
            logger.info(f"Loaded {len(case_samples)} case samples from file")
            logger.debug(f"Case samples: {case_samples[:10]}...")  # Log first 10 for debugging
        else:
            logger.debug("No case samples file provided")

        if control_samples_file:
            logger.debug(f"Loading control samples from file: {control_samples_file}")
            with open(control_samples_file) as f:
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
                f"Automatically assigned {len(auto_control_samples)} control samples "
                f"(all non-case samples)"
            )
            logger.debug(f"Auto control samples: {auto_control_samples[:10]}...")
        elif control_samples and not case_samples:
            # Automatically assign all non-control samples as cases
            control_samples_set = set(control_samples)
            vcf_samples_set = set(vcf_samples)
            auto_case_samples = list(vcf_samples_set - control_samples_set)
            context.config["case_samples"] = auto_case_samples
            logger.info(
                f"Automatically assigned {len(auto_case_samples)} case samples "
                f"(all non-control samples)"
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

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore VCF samples and case/control
        assignments that would have been loaded during normal execution.
        """
        vcf_file = context.config.get("vcf_file")
        from_tsv = context.config.get("from_tsv") or getattr(context.args, "from_tsv", None)

        if from_tsv and not vcf_file:
            from ..from_tsv import prepare_from_tsv_context

            prepare_from_tsv_context(context)
            vcf_samples = list(context.vcf_samples or [])
        else:
            if not vcf_file:
                raise ValueError("No VCF file specified in configuration during checkpoint skip")

            # Restore VCF samples
            logger.debug(f"Restoring VCF samples from {vcf_file} (checkpoint skip)")
            vcf_samples = get_vcf_samples(vcf_file)

            # Apply sample substring removal if configured
            remove_substring = context.config.get("remove_sample_substring")
            if remove_substring and remove_substring.strip():
                logger.debug(
                    f"Applying substring removal '{remove_substring}' during checkpoint skip"
                )
                vcf_samples = [s.replace(remove_substring, "") for s in vcf_samples]

            context.vcf_samples = vcf_samples

            # CRITICAL: Set sample_list in config for downstream analysis
            context.config["sample_list"] = ",".join(vcf_samples)

        _require_unique_samples(vcf_samples, context.config.get("remove_sample_substring"))

        # Restore case/control sample lists if provided via files
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        if case_samples_file and os.path.exists(case_samples_file):
            with open(case_samples_file) as f:
                case_samples = [line.strip() for line in f if line.strip()]
            context.config["case_samples"] = case_samples
            logger.debug(f"Restored {len(case_samples)} case samples from file (checkpoint skip)")

        if control_samples_file and os.path.exists(control_samples_file):
            with open(control_samples_file) as f:
                control_samples = [line.strip() for line in f if line.strip()]
            context.config["control_samples"] = control_samples
            logger.debug(
                f"Restored {len(control_samples)} control samples from file (checkpoint skip)"
            )

        # Handle automatic control assignment
        case_samples = context.config.get("case_samples", [])
        control_samples = context.config.get("control_samples", [])

        if case_samples and not control_samples:
            case_samples_set = set(case_samples)
            vcf_samples_set = set(vcf_samples)
            auto_control_samples = list(vcf_samples_set - case_samples_set)
            context.config["control_samples"] = auto_control_samples
            logger.debug(
                f"Auto-assigned {len(auto_control_samples)} control samples (checkpoint skip)"
            )
        elif control_samples and not case_samples:
            control_samples_set = set(control_samples)
            vcf_samples_set = set(vcf_samples)
            auto_case_samples = list(vcf_samples_set - control_samples_set)
            context.config["case_samples"] = auto_case_samples
            logger.debug(f"Auto-assigned {len(auto_case_samples)} case samples (checkpoint skip)")

        logger.info(
            f"Restored {len(vcf_samples)} VCF samples and case/control assignments "
            f"(checkpoint skip)"
        )
        return context


__all__ = [
    "PedigreeLoadingStage",
    "SampleConfigLoadingStage",
    "get_vcf_samples",
    "read_pedigree",
]
