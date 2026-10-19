"""Phenotype-based case/control assignment stage for VariantCentrifuge."""

import logging

import pandas as pd

from ..pipeline_core import PipelineContext, Stage

logger = logging.getLogger(__name__)


class PhenotypeCaseControlAssignmentStage(Stage):
    """Assign case/control samples based on phenotype data when no explicit files provided."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "phenotype_case_control_assignment"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Case/control assignment can be restored by the checkpoint skip handler."""
        return "restore"

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
        # Check if case/control samples are already explicitly assigned.
        # Use `or []` to handle None values (CLI sets case_samples=None when
        # --case-samples is not used, but --case-samples-file may still load them
        # via SampleConfigLoadingStage).
        existing_case_samples = context.config.get("case_samples") or []
        existing_control_samples = context.config.get("control_samples") or []
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        # Skip if case/control samples are already assigned (by SampleConfigLoadingStage
        # or by explicit --case-samples CLI argument)
        if existing_case_samples or existing_control_samples:
            logger.info(
                "Case/control samples already assigned (%d cases, %d controls), "
                "skipping re-assignment",
                len(existing_case_samples),
                len(existing_control_samples),
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

        # Check for file-based case/control assignment first
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        if case_samples_file or control_samples_file:
            logger.info("File-based case/control assignment detected")
            result = self._handle_file_based_assignment(
                context, case_samples_file, control_samples_file
            )
            return result  # type: ignore[no-any-return]

        # Skip if no phenotype criteria
        if not case_phenotypes and not control_phenotypes:
            logger.debug("No phenotype criteria specified, skipping case/control assignment")
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
                f"phenotype_file={phenotype_file}, "
                f"sample_column={phenotype_sample_column}, "
                f"value_column={phenotype_value_column}"
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
                context,
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
            raise

        return context

    def _compute_samples_from_phenotype_file(
        self,
        phenotype_file: str,
        sample_column: str,
        value_column: str,
        case_phenotypes: list[str],
        control_phenotypes: list[str],
        vcf_samples: list[str],
        remove_substring: str = "",
        context: PipelineContext | None = None,
    ) -> tuple[list[str], list[str]]:
        """Load CSV, filter by phenotype terms, get sample lists."""
        logger.info(f"Loading phenotype file: {phenotype_file}")

        # Determine delimiter
        if phenotype_file.endswith(".csv"):
            df = pd.read_csv(phenotype_file, low_memory=False)
        else:
            df = pd.read_csv(phenotype_file, sep="\t", low_memory=False)

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

        logger.info(f"Looking for case phenotypes: {case_phenotype_set}")
        logger.info(f"Looking for control phenotypes: {control_phenotype_set}")

        # Debug: show sample matching process
        logger.debug(f"First 5 phenotype file samples: {df[sample_column].unique()[:5].tolist()}")
        logger.debug(f"First 5 VCF samples: {vcf_samples[:5]}")

        # Convert sample IDs to strings for consistent comparison
        # VCF samples are strings, but phenotype file samples might be integers
        vcf_samples_str = {str(s) for s in vcf_samples}

        # CRITICAL FIX: Apply same substring removal to VCF samples as case/control samples
        if remove_substring and remove_substring.strip():
            vcf_samples_str = {s.replace(remove_substring, "") for s in vcf_samples_str}
            logger.debug(f"Applied substring removal '{remove_substring}' to VCF samples")

        pheno_samples_raw = {str(s) for s in df[sample_column].unique()}

        # Apply sample substring removal to phenotype samples if configured
        # This is necessary because VCF samples have already been processed
        if remove_substring and remove_substring.strip():
            logger.debug(
                f"Applying substring removal '{remove_substring}' to phenotype samples for matching"
            )
            pheno_samples_str = {s.replace(remove_substring, "") for s in pheno_samples_raw}
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
            case_samples_raw = {str(s) for s in case_rows[sample_column].unique()}
            # Apply substring removal if configured
            if remove_substring and remove_substring.strip():
                case_samples = {s.replace(remove_substring, "") for s in case_samples_raw}
            else:
                case_samples = case_samples_raw
            logger.info(f"Found {len(case_samples)} samples matching case phenotypes")

        # Find samples that match control phenotypes
        control_samples = set()
        if control_phenotype_set:
            control_mask = df[value_column].isin(control_phenotype_set)
            control_rows = df[control_mask]
            # Convert to strings for consistent comparison
            control_samples_raw = {str(s) for s in control_rows[sample_column].unique()}
            # Apply substring removal if configured
            if remove_substring and remove_substring.strip():
                control_samples = {s.replace(remove_substring, "") for s in control_samples_raw}
            else:
                control_samples = control_samples_raw
            logger.info(f"Found {len(control_samples)} samples matching control phenotypes")

        # DEBUG: Show the actual sample formats to diagnose mismatch
        logger.debug("DEBUGGING SAMPLE ASSIGNMENT:")
        logger.debug(f"Case samples format: {list(case_samples)[:3] if case_samples else 'None'}")
        logger.debug(f"VCF samples format: {list(vcf_samples_str)[:3]}")
        logger.debug(f"Configured remove_substring: '{remove_substring}'")
        if context:
            direct_config = context.config.get("remove_sample_substring")
            logger.debug(f"Direct config remove_sample_substring: '{direct_config}'")
        else:
            logger.debug("No context available to check direct config")

        # Debug: show what we're intersecting
        logger.debug(f"VCF samples (processed): {len(vcf_samples_str)} samples")
        logger.debug(f"Case samples from phenotypes: {len(case_samples)} samples")
        logger.debug(f"First 5 VCF samples: {list(vcf_samples_str)[:5]}")
        logger.debug(f"First 5 case samples: {list(case_samples)[:5]}")

        # Debug intersection issue
        intersection = case_samples & vcf_samples_str
        logger.debug(f"Intersection of case samples with VCF samples: {len(intersection)} samples")
        if len(intersection) > 0:
            logger.debug(f"First 5 intersecting samples: {list(intersection)[:5]}")
        else:
            logger.warning("NO INTERSECTION between case samples and VCF samples!")
            # Show sample name comparison
            case_sample_list = list(case_samples)[:10]
            vcf_sample_list = list(vcf_samples_str)[:10]
            logger.warning(f"Example case samples: {case_sample_list}")
            logger.warning(f"Example VCF samples: {vcf_sample_list}")

        # Apply assignment logic using the correct VCF samples set (with substring removal applied)
        if case_phenotype_set and not control_phenotype_set:
            # Only case phenotypes specified - all non-matching VCF samples are controls
            final_case_samples = case_samples & vcf_samples_str
            final_control_samples = vcf_samples_str - final_case_samples
        elif control_phenotype_set and not case_phenotype_set:
            # Only control phenotypes specified - all non-matching VCF samples are cases
            final_control_samples = control_samples & vcf_samples_str
            final_case_samples = vcf_samples_str - final_control_samples
        else:
            # Both specified - exact matching only
            final_case_samples = case_samples & vcf_samples_str
            final_control_samples = control_samples & vcf_samples_str
        logger.info(
            f"Final assignment: {len(final_case_samples)} cases, "
            f"{len(final_control_samples)} controls"
        )
        return list(final_case_samples), list(final_control_samples)

    def _classify_samples_by_phenotypes(
        self,
        vcf_samples: list[str],
        phenotype_data: dict[str, set[str]],
        case_phenotypes: list[str],
        control_phenotypes: list[str],
    ) -> tuple[set[str], set[str]]:
        """Classify samples based on phenotype data."""
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

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        Since this stage performs phenotype-based case/control assignment, we need to
        restore or re-compute the assignments that would have been made during normal execution.
        """
        # Check if case/control samples are already explicitly assigned in config
        existing_case_samples = context.config.get("case_samples", [])
        existing_control_samples = context.config.get("control_samples", [])
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        # If explicit assignments exist, we don't need to do anything
        if (
            existing_case_samples
            or existing_control_samples
            or case_samples_file
            or control_samples_file
        ):
            logger.debug(
                "Explicit case/control assignments exist, no restoration needed (checkpoint skip)"
            )
            return context

        # Get phenotype criteria - these should already be parsed as lists
        case_phenotypes = context.config.get("case_phenotypes", [])
        control_phenotypes = context.config.get("control_phenotypes", [])

        # If no phenotype criteria, no assignment was done
        if not case_phenotypes and not control_phenotypes:
            logger.debug("No phenotype criteria specified, no restoration needed (checkpoint skip)")
            return context

        # Get VCF samples
        vcf_samples = getattr(context, "vcf_samples", [])
        if not vcf_samples:
            raise ValueError(
                "No VCF samples found, cannot restore phenotype-based assignment (checkpoint skip)"
            )

        # Get phenotype file parameters
        phenotype_file = context.config.get("phenotype_file")
        phenotype_sample_column = context.config.get("phenotype_sample_column")
        phenotype_value_column = context.config.get("phenotype_value_column")

        if not phenotype_file or not phenotype_sample_column or not phenotype_value_column:
            raise ValueError(
                "Missing phenotype file parameters, cannot restore phenotype-based assignment "
                "(checkpoint skip)"
            )

        # Re-compute phenotype-based assignments
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
                context,
            )

            # Update configuration with classified samples
            context.config["case_samples"] = case_samples
            context.config["control_samples"] = control_samples

            logger.debug(
                f"Restored phenotype-based assignment: {len(case_samples)} cases, "
                f"{len(control_samples)} controls (checkpoint skip)"
            )

        except Exception as e:
            raise ValueError(
                f"Failed to restore phenotype-based assignment during checkpoint skip: {e}"
            ) from e

        return context

    def _handle_file_based_assignment(
        self,
        context: PipelineContext,
        case_samples_file: str | None,
        control_samples_file: str | None,
    ) -> PipelineContext:
        """Handle file-based case/control assignment."""
        logger.info("Processing file-based case/control assignment")
        # Get VCF samples
        vcf_samples = getattr(context, "vcf_samples", [])
        if not vcf_samples:
            logger.warning("No VCF samples found, cannot perform file-based assignment")
            return context
        vcf_samples_set = set(vcf_samples)
        case_samples = set()
        control_samples = set()

        # Load case samples from file
        if case_samples_file:
            logger.info(f"Loading case samples from: {case_samples_file}")
            try:
                with open(case_samples_file) as f:
                    for line in f:
                        sample = line.strip()
                        if sample:
                            case_samples.add(sample)
                logger.info(f"Loaded {len(case_samples)} case samples from file")
            except Exception as e:
                logger.error(f"Failed to load case samples from {case_samples_file}: {e}")
                raise

        # Load control samples from file
        if control_samples_file:
            logger.info(f"Loading control samples from: {control_samples_file}")
            try:
                with open(control_samples_file) as f:
                    for line in f:
                        sample = line.strip()
                        if sample:
                            control_samples.add(sample)
                logger.info(f"Loaded {len(control_samples)} control samples from file")
            except Exception as e:
                logger.error(f"Failed to load control samples from {control_samples_file}: {e}")
                raise

        # Apply substring removal if configured (same as VCF samples)
        remove_substring = context.config.get("remove_sample_substring", "")
        if remove_substring and remove_substring.strip():
            logger.info(f"Applying substring removal '{remove_substring}' to case/control samples")
            case_samples = {s.replace(remove_substring, "") for s in case_samples}
            control_samples = {s.replace(remove_substring, "") for s in control_samples}

        # Find intersection with VCF samples
        case_in_vcf = case_samples & vcf_samples_set
        control_in_vcf = control_samples & vcf_samples_set

        logger.info(f"Case samples in VCF: {len(case_in_vcf)} out of {len(case_samples)}")
        logger.info(f"Control samples in VCF: {len(control_in_vcf)} out of {len(control_samples)}")

        # Warn about samples not found in VCF
        missing_cases = case_samples - vcf_samples_set
        missing_controls = control_samples - vcf_samples_set

        if missing_cases:
            logger.warning(
                f"{len(missing_cases)} case samples not found in VCF: {list(missing_cases)[:10]}"
            )
        if missing_controls:
            logger.warning(
                f"{len(missing_controls)} control samples not found in VCF: "
                f"{list(missing_controls)[:10]}"
            )

        # Assignment logic: if only one type is specified, others become the opposite
        if case_samples_file and not control_samples_file:
            # Only case file provided - all non-case VCF samples become controls
            final_case_samples = case_in_vcf
            final_control_samples = vcf_samples_set - final_case_samples
        elif control_samples_file and not case_samples_file:
            # Only control file provided - all non-control VCF samples become cases
            final_control_samples = control_in_vcf
            final_case_samples = vcf_samples_set - final_control_samples
        else:
            # Both files provided - use intersection
            final_case_samples = case_in_vcf
            final_control_samples = control_in_vcf

        # Guard: do not overwrite existing assignments with empty results.
        # This can happen when VCF samples haven't been fully parsed yet
        # or when sample name formats don't match between file and VCF.
        if not final_case_samples and not final_control_samples:
            logger.warning(
                "File-based assignment produced 0 cases and 0 controls "
                "(no sample names matched VCF samples). "
                "Preserving any previously loaded case/control assignments."
            )
            return context

        # Update configuration with classified samples
        context.config["case_samples"] = list(final_case_samples)
        context.config["control_samples"] = list(final_control_samples)

        logger.info(
            f"File-based assignment complete: {len(final_case_samples)} cases, "
            f"{len(final_control_samples)} controls"
        )

        # Log sample assignment summary
        logger.debug(f"Case samples (first 10): {list(final_case_samples)[:10]}")
        logger.debug(f"Control samples (first 10): {list(final_control_samples)[:10]}")

        return context


__all__ = ["PhenotypeCaseControlAssignmentStage"]
