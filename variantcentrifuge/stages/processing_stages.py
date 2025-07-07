"""
Processing stages for variant extraction and data transformation.

This module contains stages that handle the core data processing tasks:
- Gene BED file creation
- Variant extraction from VCF files
- Filtering with bcftools and SnpSift
- Field extraction to TSV format
- Genotype replacement
- Phenotype integration
"""

import logging
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd

from ..extractor import extract_fields
from ..filters import apply_snpsift_filter, extract_variants
from ..gene_bed import get_gene_bed, normalize_genes
from ..phenotype import aggregate_phenotypes_for_samples
from ..pipeline_core import PipelineContext, Stage
from ..replacer import replace_genotypes
from ..utils import run_command, split_bed_file
from ..vcf_eff_one_per_line import process_vcf_file as split_snpeff_annotations

logger = logging.getLogger(__name__)


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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"configuration_loading"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - creates independent BED file

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Create BED file for specified genes."""
        # Get normalized gene names
        gene_name = context.config.get("gene_name", "")
        gene_file = context.config.get("gene_file")

        # Normalize genes
        normalized_genes = normalize_genes(gene_name, gene_file, logger)
        context.config["normalized_genes"] = normalized_genes

        # Generate BED file
        reference = context.config.get("reference", "GRCh37.75")
        interval_expand = context.config.get("interval_expand", 0)
        add_chr = context.config.get("add_chr", True)

        logger.info(f"Creating BED file for genes: {normalized_genes}")
        bed_file = get_gene_bed(
            reference=reference,
            gene_name=normalized_genes,
            interval_expand=interval_expand,
            add_chr=add_chr,
            output_dir=str(context.workspace.output_dir),
        )

        context.gene_bed_file = Path(bed_file)
        logger.info(f"Created BED file: {bed_file}")

        return context

    def get_output_files(self, context: PipelineContext) -> List[Path]:
        """Return the generated BED file."""
        if context.gene_bed_file:
            return [context.gene_bed_file]
        return []


class VariantExtractionStage(Stage):
    """Extract variants from VCF file using BED regions (single-threaded)."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "variant_extraction"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Extract variants from VCF by genomic regions"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"gene_bed_creation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract variants using the BED file."""
        vcf_file = context.config["vcf_file"]
        bed_file = str(context.gene_bed_file)

        # Output path
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.variants.vcf.gz"
        )

        # Prepare config for extract_variants
        extract_config = {
            "threads": context.config.get("threads", 1),
            "bcftools_prefilter": context.config.get("bcftools_filter"),
        }

        logger.info(f"Extracting variants from {vcf_file} using {bed_file}")
        extract_variants(
            vcf_file=vcf_file,
            bed_file=bed_file,
            cfg=extract_config,
            output_file=str(output_vcf),
        )

        context.extracted_vcf = output_vcf
        context.data = output_vcf

        return context

    def get_output_files(self, context: PipelineContext) -> List[Path]:
        """Return the extracted VCF file."""
        if context.extracted_vcf:
            return [context.extracted_vcf]
        return []


class ParallelVariantExtractionStage(Stage):
    """Extract variants in parallel by splitting BED file into chunks."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "parallel_variant_extraction"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Extract variants in parallel using multiple genomic regions"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"gene_bed_creation"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        # This stage manages its own parallelism internally
        return False

    @property
    def estimated_runtime(self) -> float:
        """Return the estimated runtime in seconds."""
        return 30.0  # Typically takes longer

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract variants in parallel by processing BED chunks."""
        threads = context.config.get("threads", 1)

        if threads <= 1:
            # Fall back to single-threaded extraction
            logger.info("Thread count <= 1, using single-threaded extraction")
            extractor = VariantExtractionStage()
            return extractor._process(context)

        # Split BED file into chunks
        bed_chunks = self._split_bed_file(context.gene_bed_file, threads)
        logger.info(f"Split BED file into {len(bed_chunks)} chunks for parallel processing")

        # Process chunks in parallel
        chunk_outputs = self._process_chunks_parallel(context, bed_chunks)

        # Merge results
        merged_output = self._merge_chunk_outputs(context, chunk_outputs)

        context.extracted_vcf = merged_output
        context.data = merged_output

        # Cleanup chunks
        self._cleanup_chunks(bed_chunks, chunk_outputs)

        return context

    def _split_bed_file(self, bed_file: Path, n_chunks: int) -> List[Path]:
        """Split BED file into roughly equal chunks."""
        chunk_dir = bed_file.parent / "chunks"
        chunk_dir.mkdir(exist_ok=True)

        chunks = split_bed_file(str(bed_file), n_chunks, str(chunk_dir))
        return [Path(chunk) for chunk in chunks]

    def _process_chunks_parallel(
        self, context: PipelineContext, bed_chunks: List[Path]
    ) -> List[Path]:
        """Process BED chunks in parallel."""
        vcf_file = context.config["vcf_file"]
        threads = context.config.get("threads", 1)

        # Each worker gets limited threads
        threads_per_worker = max(1, threads // len(bed_chunks))

        # Prepare config for extract_variants
        extract_config = {
            "threads": threads_per_worker,
            "bcftools_prefilter": context.config.get("bcftools_filter"),
        }

        chunk_outputs = []
        with ProcessPoolExecutor(max_workers=len(bed_chunks)) as executor:
            # Submit all jobs
            future_to_chunk = {}
            for i, chunk_bed in enumerate(bed_chunks):
                output_vcf = context.workspace.get_temp_path(f"chunk_{i}.variants.vcf.gz")

                future = executor.submit(
                    extract_variants,
                    vcf_file=vcf_file,
                    bed_file=str(chunk_bed),
                    cfg=extract_config,
                    output_file=str(output_vcf),
                )
                future_to_chunk[future] = (chunk_bed, output_vcf)

            # Collect results
            for future in as_completed(future_to_chunk):
                chunk_bed, output_vcf = future_to_chunk[future]
                try:
                    future.result()  # Raises any exceptions
                    chunk_outputs.append(output_vcf)
                    logger.debug(f"Completed extraction for chunk: {chunk_bed.name}")
                except Exception as e:
                    logger.error(f"Failed to process chunk {chunk_bed}: {e}")
                    raise

        return chunk_outputs

    def _merge_chunk_outputs(self, context: PipelineContext, chunk_outputs: List[Path]) -> Path:
        """Merge VCF chunks into a single file."""
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.variants.vcf.gz"
        )

        if len(chunk_outputs) == 1:
            # Just move the single file
            shutil.move(str(chunk_outputs[0]), str(output_vcf))
        else:
            # Merge with bcftools
            cmd = [
                "bcftools",
                "concat",
                "-a",  # Allow overlaps
                "-Oz",  # Output compressed
                "-o",
                str(output_vcf),
            ] + [str(f) for f in chunk_outputs]

            run_command(cmd)

            # Index the merged file
            run_command(["bcftools", "index", str(output_vcf)])

        logger.info(f"Merged {len(chunk_outputs)} chunks into {output_vcf}")
        return output_vcf

    def _cleanup_chunks(self, bed_chunks: List[Path], vcf_chunks: List[Path]) -> None:
        """Clean up temporary chunk files."""
        for chunk in bed_chunks + vcf_chunks:
            if chunk.exists():
                chunk.unlink()

        # Remove chunks directory if empty
        if bed_chunks and bed_chunks[0].parent.name == "chunks":
            chunk_dir = bed_chunks[0].parent
            if chunk_dir.exists() and not any(chunk_dir.iterdir()):
                chunk_dir.rmdir()


class BCFToolsPrefilterStage(Stage):
    """Apply bcftools pre-filtering for performance optimization."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "bcftools_prefilter"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply bcftools filter expression"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"variant_extraction", "parallel_variant_extraction"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply bcftools filter if specified."""
        bcftools_filter = context.config.get("bcftools_filter")

        if not bcftools_filter:
            logger.debug("No bcftools filter specified, skipping")
            return context

        # Note: This is typically applied during extraction for efficiency
        # This stage is here for cases where we need a separate filtering step
        logger.info(f"Applying bcftools filter: {bcftools_filter}")

        input_vcf = context.extracted_vcf or context.data
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.bcftools_filtered.vcf.gz"
        )

        cmd = [
            "bcftools",
            "view",
            "-i",
            bcftools_filter,
            "-Oz",
            "-o",
            str(output_vcf),
            str(input_vcf),
        ]

        run_command(cmd)
        run_command(["bcftools", "index", str(output_vcf)])

        context.filtered_vcf = output_vcf
        context.data = output_vcf

        return context


class MultiAllelicSplitStage(Stage):
    """Split multi-allelic SNPeff annotations into one per line."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "multiallelic_split"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Split multi-allelic SNPeff annotations"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Must run after gene bed creation and variant data exists
        return {"gene_bed_creation"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - independent transformation, creates new file

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Split SNPeff annotations if requested."""
        if not context.config.get("snpeff_split_by_transcript"):
            logger.debug("SNPeff transcript splitting not requested")
            return context

        # Determine when to split
        split_before = context.config.get("snpeff_split_before_filter", True)
        if split_before and context.is_complete("snpsift_filtering"):
            # Already filtered, don't split again
            return context

        input_vcf = context.filtered_vcf or context.extracted_vcf or context.data
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.split_annotations.vcf.gz"
        )

        logger.info("Splitting SNPeff annotations by transcript")
        split_snpeff_annotations(str(input_vcf), str(output_vcf))

        context.data = output_vcf
        return context


class SnpSiftFilterStage(Stage):
    """Apply SnpSift filter expressions."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "snpsift_filtering"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply SnpSift filter expressions"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Must run after some form of variant extraction
        return {"gene_bed_creation"}

    def _split_before_filter(self) -> bool:
        """Check if we should split before filtering."""
        # This is a simplification - in real implementation would check context
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply SnpSift filters."""
        filter_expr = context.config.get("filter")

        if not filter_expr or context.config.get("late_filtering"):
            logger.debug("No SnpSift filter or using late filtering")
            return context

        input_vcf = context.data
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.filtered.vcf.gz"
        )

        logger.info(f"Applying SnpSift filter: {filter_expr}")
        filter_config = {"threads": context.config.get("threads", 1)}
        apply_snpsift_filter(str(input_vcf), filter_expr, filter_config, str(output_vcf))

        context.filtered_vcf = output_vcf
        context.data = output_vcf

        return context


class FieldExtractionStage(Stage):
    """Extract specified fields from VCF to TSV format."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "field_extraction"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Extract fields from VCF to TSV format"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Must have gene bed and some form of extraction
        return {"gene_bed_creation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract fields to TSV format."""
        input_vcf = context.data
        fields = context.config.get("extract", [])

        if not fields:
            raise ValueError("No fields specified for extraction")

        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.extracted.tsv"
        )

        # Handle gzip compression if requested
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        logger.info(f"Extracting {len(fields)} fields from VCF to TSV")

        # Prepare config for extract_fields
        extract_config = {
            "extract_fields_separator": context.config.get("extract_fields_separator", ","),
            "debug_level": context.config.get("log_level", "INFO"),
        }

        # Convert fields list to space-separated string
        fields_str = " ".join(fields)

        extract_fields(
            variant_file=str(input_vcf),
            fields=fields_str,
            cfg=extract_config,
            output_file=str(output_tsv),
        )

        context.extracted_tsv = output_tsv
        context.data = output_tsv

        return context

    def get_output_files(self, context: PipelineContext) -> List[Path]:
        """Return the extracted TSV file."""
        if context.extracted_tsv:
            return [context.extracted_tsv]
        return []


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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"field_extraction", "sample_config_loading"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Replace genotypes with sample IDs."""
        if not context.config.get("replace_genotypes"):
            logger.debug("Genotype replacement not requested")
            return context

        input_tsv = context.data
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.genotype_replaced.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        # Get sample names
        samples = context.vcf_samples
        if not samples:
            logger.warning("No VCF samples loaded, skipping genotype replacement")
            return context

        logger.info(f"Replacing genotypes for {len(samples)} samples")
        replace_genotypes(
            str(input_tsv),
            samples,
            str(output_tsv),
            missing_string=context.config.get("missing_string", ""),
            output_to_stdout=False,
        )

        context.data = output_tsv
        return context


class PhenotypeIntegrationStage(Stage):
    """Integrate phenotype data into the variant table."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "phenotype_integration"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Add phenotype data to variant table"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        deps = {"field_extraction"}
        if self._has_genotype_replacement():
            deps.add("genotype_replacement")
        if self._has_phenotype_data():
            deps.add("phenotype_loading")
        return deps

    def _has_genotype_replacement(self) -> bool:
        """Check if genotype replacement is enabled."""
        # Simplified - would check context in real implementation
        return True

    def _has_phenotype_data(self) -> bool:
        """Check if phenotype data is available."""
        # Simplified - would check context in real implementation
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Add phenotype data to the table."""
        if not context.phenotype_data:
            logger.debug("No phenotype data to integrate")
            return context

        input_tsv = context.data
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.phenotypes_added.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        # Aggregate phenotypes
        samples = context.vcf_samples
        aggregated = aggregate_phenotypes_for_samples(
            phenotypes=context.phenotype_data,
            samples=samples,
            missing_string=context.config.get("missing_string_phenotype", ""),
        )

        logger.info(f"Adding phenotype data for {len(samples)} samples")

        # Read input TSV
        compression = "gzip" if str(input_tsv).endswith(".gz") else None
        df = pd.read_csv(input_tsv, sep="\t", dtype=str, compression=compression)

        # Add phenotype column
        df["Phenotypes"] = aggregated

        # Write output
        compression = "gzip" if str(output_tsv).endswith(".gz") else None
        df.to_csv(output_tsv, sep="\t", index=False, compression=compression)

        context.data = output_tsv
        return context


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
    def dependencies(self) -> Set[str]:
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
        df = pd.read_csv(input_tsv, sep="\t", dtype=str, compression=compression)

        # Remove columns that exist
        existing_cols = [col for col in columns_to_remove if col in df.columns]
        if existing_cols:
            df = df.drop(columns=existing_cols)
            logger.info(f"Removed columns: {existing_cols}")

        # Write output
        compression = "gzip" if str(output_tsv).endswith(".gz") else None
        df.to_csv(output_tsv, sep="\t", index=False, compression=compression)

        context.data = output_tsv
        return context


class StreamingDataProcessingStage(Stage):
    """Single-pass streaming processing for memory efficiency."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "streaming_processing"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Stream process large files in a single pass"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Can replace the separate genotype/phenotype/column stages
        return {"field_extraction"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return False  # Streaming is inherently sequential

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process file in streaming fashion for memory efficiency."""
        # Check if we should use streaming
        if not context.config.get("use_streaming_processing"):
            logger.debug("Streaming processing not enabled")
            return context

        input_file = context.data
        output_file = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.processed.tsv"
        )

        # Handle compression
        if context.config.get("gzip_intermediates"):
            output_file = Path(str(output_file) + ".gz")

        logger.info("Processing file in streaming mode for memory efficiency")

        # Get processing parameters
        replace_genotypes = context.config.get("replace_genotypes", False)
        samples = context.vcf_samples if replace_genotypes else []
        phenotypes = context.phenotype_data if context.phenotype_data else {}
        columns_to_remove = context.config.get("extra_columns_to_remove", [])

        # Stream process the file
        self._stream_process_file(
            input_file=input_file,
            output_file=output_file,
            samples=samples,
            phenotypes=phenotypes,
            columns_to_remove=columns_to_remove,
            missing_string=context.config.get("missing_string", "./."),
        )

        context.data = output_file
        return context

    def _stream_process_file(
        self,
        input_file: Path,
        output_file: Path,
        samples: List[str],
        phenotypes: Dict[str, str],
        columns_to_remove: List[str],
        missing_string: str,
    ) -> None:
        """Stream process file line by line."""
        import gzip

        # Determine compression
        in_compressed = str(input_file).endswith(".gz")
        out_compressed = str(output_file).endswith(".gz")

        # Open files
        if in_compressed:
            in_fh = gzip.open(input_file, "rt")
        else:
            in_fh = open(input_file, "r")

        if out_compressed:
            out_fh = gzip.open(output_file, "wt")
        else:
            out_fh = open(output_file, "w")

        try:
            # Process header
            header = in_fh.readline().strip()
            columns = header.split("\t")

            # Find columns to process
            gt_columns = [i for i, col in enumerate(columns) if col.endswith(".GT")]
            columns_to_remove_idx = [i for i, col in enumerate(columns) if col in columns_to_remove]

            # Update header
            new_columns = []
            for i, col in enumerate(columns):
                if i in columns_to_remove_idx:
                    continue
                if samples and i < len(columns) and columns[i].endswith(".GT"):
                    # Replace with sample name
                    sample_idx = gt_columns.index(i)
                    if sample_idx < len(samples):
                        new_columns.append(samples[sample_idx])
                    else:
                        new_columns.append(col)
                else:
                    new_columns.append(col)

            # Add phenotype column if needed
            if phenotypes:
                new_columns.append("Phenotypes")

            # Write header
            out_fh.write("\t".join(new_columns) + "\n")

            # Process data lines
            for line in in_fh:
                fields = line.strip().split("\t")

                # Build new row
                new_fields = []
                for i, field in enumerate(fields):
                    if i in columns_to_remove_idx:
                        continue

                    if samples and i in gt_columns:
                        # Replace genotype encoding
                        if field in ["0/1", "1/0", "0 < /dev/null | 1", "1|0"]:
                            new_fields.append("het")
                        elif field in ["1/1", "1|1"]:
                            new_fields.append("hom")
                        elif field in ["0/0", "0|0"]:
                            new_fields.append("ref")
                        elif field == missing_string:
                            new_fields.append("missing")
                        else:
                            new_fields.append(field)
                    else:
                        new_fields.append(field)

                # Add phenotypes if needed
                if phenotypes:
                    pheno_values = []
                    for sample in samples:
                        pheno_values.append(phenotypes.get(sample, ""))
                    new_fields.append(";".join(pheno_values))

                # Write line
                out_fh.write("\t".join(new_fields) + "\n")

        finally:
            in_fh.close()
            out_fh.close()

        logger.info(f"Streaming processing complete: {output_file}")
