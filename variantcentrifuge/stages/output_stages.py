"""
Output stages for generating final results and reports.

This module contains stages that produce the final outputs:
- Variant identifier generation
- Final filtering
- Pseudonymization
- TSV output
- Excel reports
- HTML reports
- IGV reports
- Metadata generation
- Archive creation
"""

import json
import logging
import sys
import tarfile
from pathlib import Path
from typing import Set
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from ..converter import convert_to_excel, produce_report_json
from ..filters import filter_dataframe_with_query
from ..generate_html_report import generate_html_report
from ..generate_igv_report import generate_igv_report
from ..helpers import match_IGV_link_columns
from ..links import add_links_to_table
from ..pipeline_core import PipelineContext, Stage
from ..pseudonymizer import create_pseudonymizer
from ..utils import get_tool_version, sanitize_metadata_field

logger = logging.getLogger(__name__)


class VariantIdentifierStage(Stage):
    """Generate unique variant identifiers."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "variant_identifier"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate variant identifiers"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Always depends on dataframe_loading
        deps = {"dataframe_loading"}
        # Note: Should run after variant_analysis but we handle this through
        # pipeline ordering rather than hard dependency, so the stage can be
        # used independently in tests or other contexts
        return deps

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Should run after variant_analysis because it creates a new DataFrame
        return {"variant_analysis"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return False  # Must run in main process - modifies DataFrame in-place

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Add variant identifiers."""
        # Handle streaming mode for large files
        if (
            context.config.get("stream_variant_ids")
            and context.data
            and not context.current_dataframe
        ):
            return self._process_streaming(context)

        # Standard DataFrame mode
        df = context.current_dataframe
        if df is None:
            if context.config.get("use_chunked_processing"):
                logger.debug("Variant IDs added during chunked processing")
                return context
            logger.warning("No DataFrame for variant identifier generation")
            return context

        # Use VAR_ID to match old pipeline
        id_column = "VAR_ID"
        if id_column in df.columns:
            logger.debug(f"Variant ID column '{id_column}' already exists")
            return context

        logger.info("Generating variant identifiers")

        # Generate IDs based on key fields with hash like old pipeline
        key_fields = ["CHROM", "POS", "REF", "ALT"]
        if all(field in df.columns for field in key_fields):
            import hashlib

            var_ids = []
            for idx, row in enumerate(df.itertuples(index=False), 1):
                chrom = str(getattr(row, "CHROM", ""))
                pos = str(getattr(row, "POS", ""))
                ref = str(getattr(row, "REF", ""))
                alt = str(getattr(row, "ALT", ""))

                combined = f"{chrom}{pos}{ref}{alt}"
                short_hash = hashlib.md5(combined.encode("utf-8")).hexdigest()[:4]
                var_id = f"var_{idx:04d}_{short_hash}"
                var_ids.append(var_id)

            # Insert VAR_ID as the first column
            df.insert(0, id_column, var_ids)
        else:
            # Fallback to simple index
            df.insert(0, id_column, [f"var_{i:04d}_0000" for i in range(1, len(df) + 1)])

        # Also ensure Custom_Annotation column exists (even if empty)
        if "Custom_Annotation" not in df.columns:
            # Find a good position for it - after GT column if it exists
            if "GT" in df.columns:
                gt_pos = df.columns.get_loc("GT") + 1
                df.insert(gt_pos, "Custom_Annotation", "")
            else:
                df["Custom_Annotation"] = ""

        context.current_dataframe = df
        return context

    def _process_streaming(self, context: PipelineContext) -> PipelineContext:
        """Process large files in streaming mode without loading into memory."""
        import gzip
        from pathlib import Path

        input_file = Path(context.data)
        output_file = context.workspace.get_intermediate_path("with_variant_ids.tsv.gz")
        id_column = "VAR_ID"  # Match the old pipeline

        logger.info(f"Adding variant IDs in streaming mode: {input_file} -> {output_file}")

        # Determine if input is compressed
        open_func = gzip.open if str(input_file).endswith(".gz") else open

        with open_func(input_file, "rt") as inp, gzip.open(output_file, "wt") as out:
            # Process header
            header_line = inp.readline().rstrip("\n")
            header_fields = header_line.split("\t")

            # Check if ID column already exists
            if id_column in header_fields:
                logger.debug(
                    f"Variant ID column '{id_column}' already exists, copying file unchanged"
                )
                out.write(header_line + "\n")
                for line in inp:
                    out.write(line)
                context.data = output_file
                return context

            # Find key field indices
            key_indices = {}
            for field in ["CHROM", "POS", "REF", "ALT"]:
                if field in header_fields:
                    key_indices[field] = header_fields.index(field)

            # Write new header with ID column first
            out.write(f"{id_column}\t{header_line}\n")

            # Process data lines
            line_count = 0
            use_key_fields = len(key_indices) == 4  # All key fields present
            import hashlib

            for line in inp:
                line = line.rstrip("\n")
                if not line:
                    continue

                fields = line.split("\t")
                line_count += 1

                # Generate variant ID matching regular mode format
                if use_key_fields:
                    chrom = fields[key_indices["CHROM"]]
                    pos = fields[key_indices["POS"]]
                    ref = fields[key_indices["REF"]]
                    alt = fields[key_indices["ALT"]]

                    combined = f"{chrom}{pos}{ref}{alt}"
                    short_hash = hashlib.md5(combined.encode("utf-8")).hexdigest()[:4]
                    var_id = f"var_{line_count:04d}_{short_hash}"
                else:
                    var_id = f"var_{line_count:04d}_0000"

                # Write line with ID
                out.write(f"{var_id}\t{line}\n")

                if line_count % 100000 == 0:
                    logger.debug(f"Processed {line_count:,} variants")

        logger.info(f"Added variant IDs to {line_count:,} variants")
        context.data = output_file
        return context


class FinalFilteringStage(Stage):
    """Apply final filters on computed columns."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "final_filtering"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply final filters on all columns"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depends on dataframe_loading since filtering can work on any columns
        return {"dataframe_loading"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - only filters DataFrame, no external I/O

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Final filtering often filters on computed columns from these stages
        return {"variant_scoring", "inheritance_analysis", "custom_annotation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply final filtering."""
        # Check for late filtering or final filter
        late_filtering = context.config.get("late_filtering")
        final_filter = context.config.get("final_filter")
        original_filter = context.config.get("filter") if late_filtering else None

        logger.debug(f"FinalFilteringStage - late_filtering: {late_filtering}")
        logger.debug(f"FinalFilteringStage - final_filter: {final_filter}")
        logger.debug(f"FinalFilteringStage - original_filter: {original_filter}")

        if not any([late_filtering, final_filter]):
            logger.debug("No final filtering configured")
            return context

        df = context.current_dataframe
        if df is None:
            if context.config.get("use_chunked_processing"):
                logger.debug("Final filtering applied during chunked processing")
                return context
            logger.warning("No DataFrame for final filtering")
            return context

        initial_count = len(df)
        # Debug: Check what columns are available
        logger.debug(f"FinalFilteringStage - Available columns: {list(df.columns)}")
        if "nephro_candidate_score" in df.columns:
            logger.debug("FinalFilteringStage - nephro_candidate_score column exists")
            # Show some sample values
            sample_values = df["nephro_candidate_score"].head(10).tolist()
            logger.debug(
                f"FinalFilteringStage - Sample nephro_candidate_score values: " f"{sample_values}"
            )
        else:
            logger.debug("FinalFilteringStage - nephro_candidate_score column NOT found")

        # Apply late filtering if configured
        if late_filtering and original_filter:
            logger.info(f"Applying late filter: {original_filter}")
            df = filter_dataframe_with_query(df, original_filter)
            logger.info(f"Late filter retained {len(df)}/{initial_count} variants")

        # Apply final filter if configured
        if final_filter:
            logger.info(f"Applying final filter: {final_filter}")
            pre_filter_count = len(df)
            df = filter_dataframe_with_query(df, final_filter)
            post_filter_count = len(df)
            logger.info(
                f"Final filter retained {post_filter_count}/{pre_filter_count} "
                f"variants (total: {post_filter_count}/{initial_count})"
            )

        context.current_dataframe = df
        return context


class PseudonymizationStage(Stage):
    """Apply sample pseudonymization for privacy."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "pseudonymization"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply sample pseudonymization"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depend on dataframe_loading as minimum requirement
        # Other stages (filtering, variant_id) are optional
        return {"dataframe_loading"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply pseudonymization if requested."""
        if not context.config.get("pseudonymize"):
            logger.debug("Pseudonymization not requested")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame for pseudonymization")
            return context

        # Get samples
        samples = context.vcf_samples or []
        if not samples:
            logger.warning("No samples to pseudonymize")
            return context

        # Create pseudonymizer
        schema = context.config.get("pseudonymize_schema", "sequential")
        metadata_file = context.config.get("pseudonymize_metadata")

        logger.info(f"Creating pseudonymization ({schema} schema) for {len(samples)} samples")

        # Create appropriate pseudonymizer
        if schema == "sequential":
            pseudonymizer = create_pseudonymizer(
                "sequential", prefix=context.config.get("pseudonymize_prefix", "SAMPLE")
            )
        elif schema == "categorical":
            pseudonymizer = create_pseudonymizer(
                "categorical",
                category_field=context.config.get("pseudonymize_category_field", "phenotype"),
            )
        else:
            pseudonymizer = create_pseudonymizer(schema)

        # Create mapping for samples
        metadata_dict = {}
        if metadata_file:
            with open(metadata_file, "r") as f:
                metadata_dict = json.load(f)

        mapping = pseudonymizer.create_mapping(samples, metadata_dict)

        # Apply to DataFrame
        df = pseudonymizer.pseudonymize_dataframe(df)

        # Save mapping
        mapping_file = (
            context.workspace.output_dir.parent
            / f"pseudonymization_mapping_{context.workspace.timestamp}.json"
        )
        with open(mapping_file, "w") as f:
            json.dump(mapping, f, indent=2)
        logger.info(f"Saved pseudonymization mapping to {mapping_file}")

        context.current_dataframe = df
        context.config["pseudonymization_mapping"] = mapping
        context.config["pseudonymizer"] = pseudonymizer

        return context


class TSVOutputStage(Stage):
    """Write final TSV output file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "tsv_output"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Write final TSV output"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # TSV output should run after all processing and analysis stages
        # We depend on dataframe_loading as the base requirement
        # All analysis stages that modify the dataframe should run before this
        deps = {"dataframe_loading"}

        # IMPORTANT: We must explicitly depend on variant_identifier to ensure
        # VAR_ID is added before we write the output. The runner cannot infer
        # this dependency automatically since both stages modify the same DataFrame.
        # This prevents them from running in parallel and ensures VAR_ID is present.
        deps.add("variant_identifier")

        # Note: We don't add other optional dependencies here like final_filtering
        # or pseudonymization because they may not exist in all pipelines.
        # The runner will still ensure proper ordering based on the dependency
        # graph of stages that are actually present.

        return deps

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Soft dependencies are optional - if these stages exist, we should run after them
        # This includes all stages that might modify the DataFrame
        return {"variant_scoring", "final_filtering", "pseudonymization"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Write final TSV output."""
        df = context.current_dataframe

        # Handle chunked processing
        if context.config.get("use_chunked_processing") and df is None:
            # Output already written during chunked processing
            output_file = context.config.get("chunked_output_file")
            if output_file:
                # Ensure the path is within the output directory
                output_filename = Path(output_file).name
                context.final_output_path = context.workspace.output_dir / output_filename
                logger.info(f"Chunked processing output: {context.final_output_path}")
            return context

        if df is None:
            logger.error("No data to write")
            return context

        # Determine output path
        output_file = context.config.get("output_file")
        if output_file in [None, "stdout", "-"]:
            # Write to stdout
            logger.info("Writing output to stdout")
            df.to_csv(sys.stdout, sep="\t", index=False, na_rep="")
            return context

        # Write to file
        output_file_config = context.config.get("output_file")
        logger.debug(f"TSVOutputStage - output_file from config: {output_file_config}")

        if output_file_config:
            # Use specified filename but place it in the output directory
            output_filename = Path(output_file_config).name
            output_path = context.workspace.output_dir / output_filename
            logger.debug(f"Using specified filename in output dir: {output_path}")
        else:
            # Generate default output path
            output_path = context.workspace.get_output_path("", ".tsv")
            logger.debug(f"Using default output path: {output_path}")

        logger.info(f"Writing {len(df)} variants to {output_path}")

        # Add external links if configured (default is to add links unless no_links is True)
        if not context.config.get("no_links", False):
            # Get link configurations from config
            link_configs = context.config.get("links", {})
            if link_configs:
                # Convert DataFrame to lines, add links, then convert back
                # First, save to lines format
                lines = []
                # Header
                lines.append("\t".join(df.columns))
                # Data rows
                for _, row in df.iterrows():
                    lines.append("\t".join(str(val) for val in row))

                # Add links
                lines_with_links = add_links_to_table(lines, link_configs)

                # Convert back to DataFrame
                import io

                df = pd.read_csv(io.StringIO("\n".join(lines_with_links)), sep="\t", dtype=str)

        # Write file
        compression = None
        if context.config.get("gzip_output") or str(output_path).endswith(".gz"):
            compression = "gzip"
            if not str(output_path).endswith(".gz"):
                output_path = Path(str(output_path) + ".gz")

        df.to_csv(output_path, sep="\t", index=False, na_rep="", compression=compression)

        logger.info(f"Successfully wrote output to: {output_path}")

        context.final_output_path = output_path
        context.report_paths["tsv"] = output_path

        return context


class ExcelReportStage(Stage):
    """Generate Excel report with multiple sheets."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "excel_report"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate Excel report"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"tsv_output"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate Excel report."""
        if not context.config.get("xlsx") and not context.config.get("excel"):
            logger.debug("Excel report not requested")
            return context

        input_file = context.final_output_path
        if not input_file or not input_file.exists():
            logger.warning("No input file for Excel generation")
            return context

        output_path = context.workspace.get_output_path("", ".xlsx")

        logger.info(f"Generating Excel report: {output_path}")

        # Convert TSV to Excel
        convert_to_excel(
            str(input_file),
            str(output_path),
            add_stats_sheet=not context.config.get("no_stats"),
            stats_file=context.config.get("stats_output_file"),
        )

        context.report_paths["excel"] = output_path
        return context


class HTMLReportStage(Stage):
    """Generate interactive HTML report."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "html_report"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate HTML report"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"tsv_output"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate HTML report."""
        if not context.config.get("html_report"):
            logger.debug("HTML report not requested")
            return context

        input_file = context.final_output_path
        if not input_file or not input_file.exists():
            logger.warning("No input file for HTML generation")
            return context

        # Create report directory
        report_dir = context.workspace.output_dir / "report"
        report_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Generating JSON data for HTML report")

        # Call produce_report_json with correct signature (TSV file, output_dir)
        # This will create variants.json and summary.json in output_dir/report/
        produce_report_json(str(input_file), str(context.workspace.output_dir))

        # Path to the generated JSON files
        variants_json = report_dir / "variants.json"
        summary_json = report_dir / "summary.json"

        if not variants_json.exists() or not summary_json.exists():
            logger.error(f"Failed to generate JSON files: {variants_json}, {summary_json}")
            return context

        logger.info(f"Generating HTML report in: {report_dir}")

        # Call generate_html_report with correct signature
        generate_html_report(
            variants_json=str(variants_json),
            summary_json=str(summary_json),
            output_dir=str(report_dir),
            cfg=context.config,  # Pass the full configuration
        )

        # Check if HTML was generated successfully
        html_file = report_dir / "index.html"
        if html_file.exists():
            context.report_paths["html"] = html_file
            logger.info("HTML report generated successfully.")
        else:
            logger.error("HTML report generation failed - index.html not found")

        context.report_paths["json"] = variants_json

        return context


class IGVReportStage(Stage):
    """Generate IGV.js visualization report."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "igv_report"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate IGV.js report"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"tsv_output"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate IGV report."""
        if not context.config.get("igv_enabled"):
            logger.debug("IGV report not requested")
            return context

        # Check required inputs
        input_file = context.final_output_path
        if not input_file or not input_file.exists():
            logger.warning("No input file for IGV generation")
            return context

        bam_mapping_file = context.config.get("bam_mapping_file")
        if not bam_mapping_file:
            logger.warning("No BAM mapping file specified for IGV report")
            return context

        # Create IGV reports directory
        report_dir = context.workspace.output_dir / "report"
        report_dir.mkdir(parents=True, exist_ok=True)

        # Get IGV parameters
        igv_reference = context.config.get("igv_reference")
        igv_fasta = context.config.get("igv_fasta")
        igv_ideogram = context.config.get("igv_ideogram")

        logger.info("Starting IGV report generation.")

        try:
            # Use threads config for parallel IGV generation, defaulting to 4 workers
            max_workers = min(context.config.get("threads", 4), 8)  # Cap at 8 to avoid overwhelming system
            
            generate_igv_report(
                variants_tsv=str(input_file),
                output_dir=str(report_dir),
                bam_mapping_file=bam_mapping_file,
                igv_reference=igv_reference,
                integrate_into_main=True,  # Always integrate if IGV is enabled
                igv_fasta=igv_fasta,
                igv_ideogram=igv_ideogram,
                igv_max_allele_len_filename=context.config.get("igv_max_allele_len_filename", 10),
                igv_hash_len_filename=context.config.get("igv_hash_len_filename", 6),
                igv_max_variant_part_filename=context.config.get("igv_max_variant_part_filename", 50),
                igv_flanking=context.config.get("igv_flanking", 50),
                max_workers=max_workers,
            )
            
            logger.info("IGV reports and mapping file generated.")
            
            # Store the IGV reports path
            igv_reports_dir = report_dir / "igv"
            context.report_paths["igv"] = igv_reports_dir
            
        except Exception as e:
            logger.error(f"Failed to generate IGV reports: {e}")
            # Don't fail the pipeline, just skip IGV reports
            
        return context


class MetadataGenerationStage(Stage):
    """Generate metadata file with run information."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "metadata_generation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate metadata file"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"tsv_output"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate metadata file."""
        if context.config.get("no_metadata"):
            logger.debug("Metadata generation disabled")
            return context

        metadata_path = context.workspace.get_output_path("_metadata", ".json")

        logger.info(f"Generating metadata: {metadata_path}")

        # Collect metadata
        metadata = {
            "pipeline_version": context.config.get("pipeline_version", "unknown"),
            "run_date": context.start_time.isoformat(),
            "execution_time": f"{context.get_execution_time():.1f} seconds",
            "input_files": {
                "vcf": context.config.get("vcf_file"),
                "genes": context.config.get("normalized_genes"),
                "phenotypes": context.config.get("phenotype_file"),
                "pedigree": context.config.get("ped_file"),
            },
            "parameters": {
                "reference": context.config.get("reference"),
                "filter": context.config.get("filter"),
                "threads": context.config.get("threads"),
                "inheritance_mode": context.config.get("inheritance_mode"),
                "scoring_config": context.config.get("scoring_config_path"),
            },
            "tool_versions": {
                "bcftools": get_tool_version("bcftools"),
                "snpeff": get_tool_version("snpEff"),
                "snpsift": get_tool_version("SnpSift"),
                "bedtools": get_tool_version("bedtools"),
            },
            "output_files": {name: str(path) for name, path in context.report_paths.items()},
            "statistics": context.statistics if context.statistics else None,
        }

        # Sanitize metadata
        metadata = sanitize_metadata_field(metadata)

        # Write metadata
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2, default=str)

        context.report_paths["metadata"] = metadata_path
        return context


class ArchiveCreationStage(Stage):
    """Create compressed archive of results."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "archive_creation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Create results archive"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Archive should run after all outputs
        # Just depend on TSV output as minimum
        return {"tsv_output"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - creates file in parent directory, no conflicts

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Create archive if requested."""
        if not context.config.get("archive_results"):
            logger.debug("Archive creation not requested")
            return context

        archive_path = context.workspace.get_archive_path()

        logger.info(f"Creating archive: {archive_path}")

        # Create tarball
        with tarfile.open(archive_path, "w:gz") as tar:
            # Add output directory
            tar.add(context.workspace.output_dir, arcname=context.workspace.output_dir.name)

        logger.info(
            f"Archive created: {archive_path} ({archive_path.stat().st_size / 1024 / 1024:.1f} MB)"
        )

        context.report_paths["archive"] = archive_path
        return context


class ParallelReportGenerationStage(Stage):
    """Generate all reports in parallel."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "parallel_report_generation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate all reports in parallel"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        return {"tsv_output"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return False  # Manages its own parallelism

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate all requested reports in parallel."""
        # Check which reports are requested
        reports_to_generate = []

        if context.config.get("xlsx") or context.config.get("excel"):
            reports_to_generate.append(("excel", ExcelReportStage()))

        if context.config.get("html_report"):
            reports_to_generate.append(("html", HTMLReportStage()))

        if context.config.get("igv_enabled"):
            reports_to_generate.append(("igv", IGVReportStage()))

        if not context.config.get("no_metadata"):
            reports_to_generate.append(("metadata", MetadataGenerationStage()))

        if not reports_to_generate:
            logger.debug("No reports to generate")
            return context

        if len(reports_to_generate) == 1:
            # Just run the single report
            _, stage = reports_to_generate[0]
            return stage._process(context)

        logger.info(f"Generating {len(reports_to_generate)} reports in parallel")

        # Run reports in parallel
        with ThreadPoolExecutor(max_workers=len(reports_to_generate)) as executor:
            futures = {}
            for report_type, stage in reports_to_generate:
                future = executor.submit(stage._process, context)
                futures[future] = report_type

            # Wait for completion
            for future in as_completed(futures):
                report_type = futures[future]
                try:
                    future.result()
                    logger.debug(f"Completed {report_type} report generation")
                except Exception as e:
                    logger.error(f"Failed to generate {report_type} report: {e}")

        return context
