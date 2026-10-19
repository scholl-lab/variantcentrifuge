"""BED generation and restriction stages."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

from ..gene_bed import get_gene_bed, normalize_genes
from ..pipeline_core import PipelineContext, Stage
from ..pipeline_core.error_handling import (
    FileFormatError,
    ToolNotFoundError,
    graceful_error_handling,
    retry_on_failure,
    validate_file_exists,
)

logger = logging.getLogger(__name__)

__all__ = ["GeneBedCreationStage"]


def _ps_attr(name: str, fallback: Any) -> Any:
    """Get attribute from this module or processing_stages (supports test patching)."""
    this_mod = sys.modules.get(__name__)
    if this_mod is not None:
        val = getattr(this_mod, name, None)
        if val is not None and val is not fallback:
            return val
    ps_mod = sys.modules.get("variantcentrifuge.stages.processing_stages")
    if ps_mod is not None:
        val = getattr(ps_mod, name, None)
        if val is not None:
            return val
    return fallback


class GeneBedCreationStage(Stage):
    """Convert gene names to BED file with genomic intervals."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "gene_bed_creation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Create BED file from gene names"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"configuration_loading"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - creates independent BED file

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Create BED file for specified genes."""
        with graceful_error_handling(self.name):
            validate_fn = _ps_attr("validate_annotation_consequence_filtering", None)
            if validate_fn is None:
                from .processing_stages import (
                    validate_annotation_consequence_filtering as validate_fn,
                )
            validate_fn(context.config)

            # Get normalized gene names
            gene_name = context.config.get("gene_name", "")
            gene_file = context.config.get("gene_file")

            # Validate gene file if provided
            if gene_file:
                gene_file = str(validate_file_exists(gene_file, self.name))

            # Normalize genes
            normalize_fn = _ps_attr("normalize_genes", normalize_genes)
            try:
                normalized_genes = normalize_fn(gene_name, gene_file, logger)
                context.config["normalized_genes"] = normalized_genes
            except SystemExit as e:
                # Convert SystemExit to proper exception
                # Check the specific error condition
                if not gene_name and not gene_file:
                    raise ValueError(
                        "No gene name provided. Use --gene-name or --gene-file to specify genes."
                    ) from e
                elif gene_file and not os.path.exists(gene_file):
                    raise FileNotFoundError(f"Gene file not found: {gene_file}") from e
                elif gene_name and os.path.exists(gene_name):
                    raise ValueError(
                        f"File path '{gene_name}' provided to --gene-name. "
                        f"Use --gene-file for gene list files."
                    ) from e
                else:
                    # Generic error for other SystemExit cases
                    raise ValueError("Invalid gene specification") from e

            # Generate BED file
            reference = context.config.get("reference", "GRCh37.75")
            interval_expand = context.config.get("interval_expand", 0)
            add_chr = context.config.get("add_chr", True)
            regions_bed = context.config.get("regions_bed")

            if regions_bed and normalized_genes.lower().strip() == "all":
                context.gene_bed_file = validate_file_exists(regions_bed, self.name)
                logger.info(
                    "Using regions BED directly because gene selection is 'all': %s",
                    context.gene_bed_file,
                )
                return context

            logger.info(f"Creating BED file for genes: {normalized_genes}")

            # Retry logic for external tool failures
            @retry_on_failure(
                max_attempts=3, delay=2.0, exceptions=(subprocess.CalledProcessError,)
            )
            def create_bed() -> str:
                get_bed_fn = _ps_attr("get_gene_bed", get_gene_bed)
                res = get_bed_fn(
                    reference=reference,
                    gene_name=normalized_genes,
                    interval_expand=interval_expand,
                    add_chr=add_chr,
                    output_dir=str(context.workspace.output_dir),
                )
                return str(res)

            try:
                bed_file = create_bed()
            except subprocess.CalledProcessError as e:
                # Check if it's a gene not found error before the broader tool-not-found case.
                if "not found in database" in str(e):
                    snpeff_reason = getattr(e, "stderr", None) or str(e)
                    raise FileFormatError(
                        normalized_genes,
                        (
                            f"valid gene names for {reference}; "
                            f"snpEff genes2bed reported: {snpeff_reason.strip()}"
                        ),
                        self.name,
                    ) from e
                # Check if it's a tool not found error
                if "snpEff" in str(e) and ("not found" in str(e) or "No such file" in str(e)):
                    raise ToolNotFoundError("snpEff", self.name) from e
                raise

            context.gene_bed_file = Path(bed_file)
            logger.info(f"Created BED file: {bed_file}")

            # Region restriction: intersect gene BED with restriction BED
            if regions_bed:
                context.gene_bed_file = self._intersect_with_restriction_bed(
                    context.gene_bed_file, regions_bed, context
                )
                logger.info(f"Gene BED restricted to regions: {context.gene_bed_file}")

            return context

    def _read_chromosomes(self, bed_file: Path) -> set[str]:
        """Read chromosome names from the first column of a BED file."""
        chromosomes: set[str] = set()
        with open(bed_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                chrom = line.split("\t")[0]
                if chrom:
                    chromosomes.add(chrom)
        return chromosomes

    def _count_regions(self, bed_file: Path) -> int:
        """Count non-comment, non-empty lines in a BED file."""
        count = 0
        with open(bed_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                count += 1
        return count

    def _intersect_with_restriction_bed(
        self, gene_bed: Path, restriction_bed: str, context: PipelineContext
    ) -> Path:
        """Intersect gene BED with a restriction BED file (e.g., capture kit)."""
        restriction_path = Path(restriction_bed)
        if not restriction_path.is_file():
            raise FileNotFoundError(f"Restriction BED file not found: {restriction_bed}")

        # Detect chr-prefix mismatch
        gene_chroms = self._read_chromosomes(gene_bed)
        restriction_chroms = self._read_chromosomes(restriction_path)

        if gene_chroms and restriction_chroms:
            gene_has_chr = any(c.startswith("chr") for c in gene_chroms)
            restriction_has_chr = any(c.startswith("chr") for c in restriction_chroms)
            if gene_has_chr != restriction_has_chr:
                raise ValueError(
                    "Chromosome naming mismatch between gene BED and restriction BED: "
                    f"gene BED {'has' if gene_has_chr else 'lacks'} 'chr' prefix "
                    f"(e.g., {next(iter(gene_chroms))}), "
                    f"restriction BED {'has' if restriction_has_chr else 'lacks'} 'chr' prefix "
                    f"(e.g., {next(iter(restriction_chroms))}). "
                    "Ensure both files use the same chromosome naming convention."
                )

        # Count original regions for logging
        orig_count = self._count_regions(gene_bed)

        # Run bedtools intersect
        out_path = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.restricted.bed"
        )
        with open(out_path, "w") as out_f:
            subprocess.run(
                ["bedtools", "intersect", "-a", str(gene_bed), "-b", str(restriction_bed)],
                stdout=out_f,
                check=True,
            )

        # Check result is non-empty
        if out_path.stat().st_size == 0:
            raise ValueError(
                f"Intersection of gene BED with restriction BED produced no regions. "
                f"Gene BED: {gene_bed}, Restriction BED: {restriction_bed}. "
                "Check that both files cover overlapping genomic regions."
            )

        new_count = self._count_regions(out_path)
        logger.info(
            f"Region restriction: {new_count} of {orig_count} gene regions retained "
            f"after intersection with {restriction_bed}"
        )

        return out_path

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the generated BED file."""
        if context.gene_bed_file:
            return [context.gene_bed_file]
        return []
