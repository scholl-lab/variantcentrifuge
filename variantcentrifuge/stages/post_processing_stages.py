"""
Post-extraction processing stages for VariantCentrifuge.

Contains:
- GenotypeReplacementStage
- PhenotypeIntegrationStage
- ExtraColumnRemovalStage
- DataSortingStage
- PCAComputationStage
"""

import contextlib
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, cast

import pandas as pd

from ..phenotype import extract_phenotypes_for_gt_row, extract_phenotypes_from_sample_columns
from ..pipeline_core import PipelineContext, Stage

logger = logging.getLogger(__name__)


class GenotypeReplacementStage(Stage):
    """Replace genotype values with sample IDs."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "genotype_replacement"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Replace genotype values with sample IDs"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"field_extraction", "sample_config_loading"}

    @property
    def checkpoint_resume_policy(self) -> str:
        """Phase-11 no-op can be skipped without restoring a file."""
        return "skip"

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Clean no-op restore for Phase-11 genotype replacement."""
        logger.info("Genotype replacement is a Phase-11 no-op; nothing to restore")
        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """No-op: genotype replacement eliminated in Phase 11.

        GT reconstruction is deferred to output time (TSV/Excel stages).
        Raw per-sample GT columns from bcftools flow directly to analysis.
        """
        logger.info(
            "Genotype replacement skipped — raw per-sample GT columns flow to analysis (Phase 11)"
        )
        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        if context.data:
            return [context.data]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        return []


class PhenotypeIntegrationStage(Stage):
    """Integrate phenotype data into the variant table."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "phenotype_integration"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Phenotype integration can restore its output path or no-op cleanly."""
        return "restore"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Add phenotype data to variant table"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Field extraction may be performed by FieldExtractionStage in sequential mode
        # or internally by ParallelCompleteProcessingStage in threaded mode.
        return set()

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # These stages should run before if they exist in the pipeline
        return {
            "field_extraction",
            "parallel_complete_processing",
            "data_sorting",
            "genotype_replacement",
            "phenotype_loading",
        }

    def _has_genotype_replacement(self, context: PipelineContext) -> bool:
        """Check if genotype replacement is enabled."""
        # Check if genotype replacement is not disabled and GT field is present
        if context.config.get("no_replacement", False):
            return False

        # Check if genotype_replaced_tsv exists in context
        return (
            hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv is not None
        )

    def _has_phenotype_data(self, context: PipelineContext) -> bool:
        """Check if phenotype data is available."""
        # Check if phenotypes were loaded
        if not hasattr(context, "phenotype_data") or context.phenotype_data is None:
            return False

        # Check if phenotypes dict is populated
        return bool(context.phenotype_data)

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the context paths
        so that subsequent stages can find the phenotypes_added_tsv file.
        """
        # Reconstruct the expected output path
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.phenotypes_added.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        # Check if the file exists and restore context
        if output_tsv.exists():
            context.phenotypes_added_tsv = output_tsv
            context.data = output_tsv
            logger.info(f"Restored phenotypes_added_tsv path to {output_tsv} after checkpoint skip")
        elif not self._has_phenotype_data(context):
            logger.info("No phenotype data configured; nothing to restore")
        else:
            raise RuntimeError(
                "Cannot restore phenotype_integration from checkpoint: "
                f"expected TSV missing: {output_tsv}"
            )

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Add phenotype data to the table."""
        if not context.phenotype_data:
            logger.debug("No phenotype data to integrate")
            return context

        # Determine input file - use genotype_replaced if available, otherwise extracted
        if hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            input_tsv = context.genotype_replaced_tsv
        elif hasattr(context, "extracted_tsv") and context.extracted_tsv:
            input_tsv = context.extracted_tsv
        else:
            input_tsv = context.data
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.phenotypes_added.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        logger.info(f"Adding phenotype data for {len(context.vcf_samples)} samples")

        # Read input TSV
        compression = "gzip" if str(input_tsv).endswith(".gz") else None
        df = pd.read_csv(input_tsv, sep="\t", dtype=str, compression=cast(Any, compression))

        # Check for per-sample columns (sample names, GEN[i].GT, or GEN_i__GT)
        has_sample_columns = bool(context.vcf_samples) and (
            all(sample in df.columns for sample in context.vcf_samples)
            or any(
                f"GEN[{i}].GT" in df.columns or f"GEN_{i}__GT" in df.columns
                for i in range(len(context.vcf_samples))
            )
        )

        if has_sample_columns:
            # Phase 11: Use per-sample columns directly
            logger.debug("Using per-sample GT columns for phenotype extraction")
            df["Phenotypes"] = df.apply(  # type: ignore[call-overload]
                lambda row: extract_phenotypes_from_sample_columns(
                    row, context.vcf_samples, context.phenotype_data
                ),
                axis=1,
            )
        elif "GT" in df.columns:
            # Backwards compatibility: Use packed GT column (from genotype replacement)
            logger.debug("Using packed GT column for phenotype extraction (legacy mode)")
            phenotype_data: dict[str, set[str]] = context.phenotype_data or {}
            df["Phenotypes"] = df["GT"].apply(
                lambda gt_val: extract_phenotypes_for_gt_row(gt_val, phenotype_data)
            )
        else:
            logger.warning("No GT column or per-sample columns found, cannot add phenotype data")
            return context

        # Write output
        compression = "gzip" if str(output_tsv).endswith(".gz") else None
        df.to_csv(output_tsv, sep="\t", index=False, compression=cast(Any, compression))

        context.phenotypes_added_tsv = output_tsv
        context.data = output_tsv
        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        if hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            return [context.genotype_replaced_tsv]
        elif hasattr(context, "extracted_tsv") and context.extracted_tsv:
            return [context.extracted_tsv]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        if hasattr(context, "phenotypes_added_tsv") and context.phenotypes_added_tsv:
            return [context.phenotypes_added_tsv]
        return []


class ExtraColumnRemovalStage(Stage):
    """Remove extra columns from the TSV file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "extra_column_removal"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Remove specified columns from TSV"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"phenotype_integration", "genotype_replacement", "field_extraction"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Remove extra columns if specified."""
        columns_to_remove = context.config.get("extra_columns_to_remove", [])

        if not columns_to_remove:
            logger.debug("No columns to remove")
            return context

        input_tsv = context.data
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.columns_removed.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        logger.info(f"Removing {len(columns_to_remove)} columns")

        # Read, remove columns, and write
        compression = "gzip" if str(input_tsv).endswith(".gz") else None
        df = pd.read_csv(input_tsv, sep="\t", dtype=str, compression=cast(Any, compression))

        # Remove columns that exist
        existing_cols = [col for col in columns_to_remove if col in df.columns]
        if existing_cols:
            df = df.drop(columns=existing_cols)
            logger.info(f"Removed columns: {existing_cols}")

        # Write output
        compression = "gzip" if str(output_tsv).endswith(".gz") else None
        df.to_csv(output_tsv, sep="\t", index=False, compression=cast(Any, compression))

        context.extra_columns_removed_tsv = output_tsv
        context.data = output_tsv
        return context


class DataSortingStage(Stage):
    """Sort extracted TSV data by gene column for efficient chunked processing."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "data_sorting"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Sort TSV data by gene column for efficient processing"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"field_extraction"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Sort the extracted TSV by gene column if gene column is present."""
        from . import processing_stages as _ps_mod

        if not context.extracted_tsv:
            logger.warning("No extracted TSV found for sorting")
            return context

        # Check if sorting is needed and resolve gene column
        gene_col_cfg = context.config.get("gene_column_name")
        has_gene_column = False
        gene_column = "GENE"

        try:
            with _ps_mod.smart_open(str(context.extracted_tsv), "r") as f:
                header = f.readline().strip()
                if header:
                    columns = header.split("\t")
                    if gene_col_cfg and gene_col_cfg != "GENE":
                        gene_column = gene_col_cfg
                        has_gene_column = gene_column in columns
                    elif "GENE" in columns:
                        gene_column = "GENE"
                        has_gene_column = True
                    elif "ANN[0].GENE" in columns:
                        gene_column = "ANN[0].GENE"
                        has_gene_column = True
                    elif "ANN_0__GENE" in columns:
                        gene_column = "ANN_0__GENE"
                        has_gene_column = True
        except Exception as e:
            logger.warning(f"Could not read header from {context.extracted_tsv}: {e}")
            return context

        if not has_gene_column:
            logger.info(f"Gene column '{gene_column}' not found in TSV, skipping sorting")
            return context

        # Create sorted output file
        sorted_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.extracted.sorted.tsv"
        )

        # Handle gzip compression if requested
        if context.config.get("gzip_intermediates"):
            sorted_tsv = Path(str(sorted_tsv) + ".gz")

        logger.info(f"Sorting TSV by gene column '{gene_column}' for efficient processing")

        # Use the standard sort function
        self._sort_tsv_by_gene(
            input_file=str(context.extracted_tsv),
            output_file=str(sorted_tsv),
            gene_column=gene_column,
            temp_dir=str(context.workspace.intermediate_dir),
            memory_limit=context.config.get("sort_memory_limit", "2G"),
            parallel=context.config.get("sort_parallel", 4),
        )

        # Clean up unsorted file if not keeping intermediates
        if not context.config.get("keep_intermediates", False):
            with contextlib.suppress(OSError):
                os.remove(str(context.extracted_tsv))

        # Update context with sorted file
        context.extracted_tsv = sorted_tsv
        context.data = sorted_tsv

        return context

    def _sort_tsv_by_gene(
        self,
        input_file: str,
        output_file: str,
        gene_column: str = "GENE",
        temp_dir: str | None = None,
        memory_limit: str = "2G",
        parallel: int = 4,
    ) -> str:
        """
        Sort a TSV file by gene column using external sort command for memory efficiency.

        This is adapted from the legacy sort_tsv_by_gene function.
        """
        from . import processing_stages as _ps_mod

        logger.info(f"Sorting TSV file by gene column: {input_file} -> {output_file}")
        logger.info(f"Using memory limit: {memory_limit}, parallel threads: {parallel}")

        # First, find the gene column index
        with _ps_mod.smart_open(input_file, "r") as f:
            header = f.readline().strip()
            if not header:
                raise ValueError(f"Input file '{input_file}' is empty or has no header")
            columns = header.split("\t")

        try:
            gene_col_idx = columns.index(gene_column) + 1  # sort uses 1-based indexing
        except ValueError as e:
            raise ValueError(f"Gene column '{gene_column}' not found in TSV file") from e

        logger.info(f"Found gene column '{gene_column}' at position {gene_col_idx}")

        # Resolve full paths for shell utilities to avoid PATH issues in
        # subprocess shells (e.g., WSL, conda envs, non-interactive bash).
        import shlex

        gzip_bin = shutil.which("gzip") or "gzip"
        sort_bin = shutil.which("sort") or "sort"
        tail_bin = shutil.which("tail") or "tail"

        # Build sort arguments with memory efficiency options
        sort_args = [
            f"-k{gene_col_idx},{gene_col_idx}",
            f"--buffer-size={memory_limit}",
            f"--parallel={parallel}",
            "--stable",
            f"--compress-program={gzip_bin}",
        ]
        if temp_dir:
            sort_args.extend(["-T", temp_dir])

        sort_cmd = shlex.join(sort_args)
        safe_input, safe_output = shlex.quote(input_file), shlex.quote(output_file)
        safe_gzip, safe_sort, safe_tail = (
            shlex.quote(gzip_bin),
            shlex.quote(sort_bin),
            shlex.quote(tail_bin),
        )

        # Build the complete command based on compression
        if input_file.endswith(".gz"):
            pipe_in = (
                f"{safe_gzip} -cd {safe_input} | "
                f"{{ IFS= read -r header; echo \"$header\"; {safe_sort} {sort_cmd} -t$'\\t'; }}"
            )
        else:
            pipe_in = (
                f'{{ IFS= read -r header < {safe_input}; echo "$header"; '
                f"{safe_tail} -n +2 {safe_input} | {safe_sort} {sort_cmd} -t$'\\t'; }}"
            )

        cmd = (
            f"{pipe_in} | {safe_gzip} -c > {safe_output}"
            if output_file.endswith(".gz")
            else f"{pipe_in} > {safe_output}"
        )

        # Execute the command
        logger.debug(f"Running sort command: {cmd}")

        # Find bash executable
        bash_candidate = shutil.which("bash") or "/bin/bash"
        bash_path: str | None = bash_candidate if os.path.exists(bash_candidate) else None
        if bash_path is None:
            logger.warning(f"bash not found at {bash_candidate}, falling back to shell default")

        # Execute command
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, executable=bash_path
        )

        if result.returncode != 0:
            logger.error(f"Sort command failed: {result.stderr}")
            raise RuntimeError(f"Failed to sort TSV file: {result.stderr}")

        # Verify output file was created
        if not os.path.exists(output_file):
            raise RuntimeError(f"Output file '{output_file}' was not created")

        # Check if output file has content
        with _ps_mod.smart_open(output_file, "r") as f:
            first_line = f.readline()
            if not first_line:
                raise RuntimeError(f"Output file '{output_file}' is empty")

        logger.info("Successfully sorted TSV file by gene column")
        return output_file

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the sorted TSV file."""
        if context.extracted_tsv:
            return [context.extracted_tsv]
        return []


class PCAComputationStage(Stage):
    """Compute or load PCA eigenvectors for association covariate adjustment."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "pca_computation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Compute or load PCA eigenvectors"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"configuration_loading"}

    @property
    def parallel_safe(self) -> bool:
        """PCA writes to its own file — safe to run alongside other stages."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Compute or load PCA eigenvectors."""
        pca_arg = context.config.get("pca")
        if not pca_arg:
            return context

        n_components = context.config.get("pca_components", 10)

        if os.path.isfile(pca_arg):
            pca_file = pca_arg
            logger.info(f"PCA: using pre-computed eigenvectors from {pca_file}")
        elif pca_arg == "akt":
            pca_file = self._run_akt(context, n_components)
        else:
            raise ValueError(
                f"--pca argument '{pca_arg}' is neither a valid file path "
                "nor a recognized tool name ('akt')."
            )

        # Store pca_file in config so _build_assoc_config_from_context picks it up
        context.config["pca_file"] = pca_file
        context.mark_complete(self.name, result={"pca_file": pca_file})
        return context

    def _run_akt(self, context: PipelineContext, n_components: int) -> str:
        """Run AKT PCA computation and return path to eigenvector file."""
        vcf_file = context.config.get("vcf_file")
        if not vcf_file:
            raise ValueError("--pca akt requires a VCF file (--vcf-file)")

        output_path = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.pca.eigenvec"
        )

        # Simple cache: reuse if output exists and is non-empty
        if output_path.exists() and output_path.stat().st_size > 0:
            logger.info(f"PCA: reusing cached eigenvectors from {output_path}")
            return str(output_path)

        pca_sites = context.config.get("pca_sites")
        cmd = ["akt", "pca", vcf_file, "-N", str(n_components)]
        if pca_sites:
            cmd.extend(["-R", pca_sites])
        else:
            cmd.append("--force")
        logger.info(f"Running AKT PCA: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output_path.write_text(result.stdout)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"AKT PCA failed (exit code {e.returncode}): {e.stderr}") from e
        except FileNotFoundError as e:
            raise RuntimeError(
                "AKT not found in PATH. Install akt or use --pca <eigenvec_file> "
                "to supply pre-computed eigenvectors."
            ) from e

        if not output_path.exists() or output_path.stat().st_size == 0:
            raise RuntimeError(
                "AKT PCA completed but produced no output. "
                "Check that the VCF has sufficient samples and variants."
            )

        return str(output_path)
