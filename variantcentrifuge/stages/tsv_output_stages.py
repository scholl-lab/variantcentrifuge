"""
TSV and DataFrame output stages for VariantCentrifuge.

Contains:
- VariantIdentifierStage
- FinalFilteringStage
- TSVOutputStage
- ArchiveCreationStage
"""

import gzip
import logging
import shutil
import sys
import tarfile
from pathlib import Path
from typing import Any, cast

import pandas as pd

from ..checkpoint_models import AtomicFileOperation
from ..config import normalize_annotation_config
from ..gt_reconstruction import _collapse_per_sample_gt_columns_for_output
from ..links import add_links_to_table, resolve_links
from ..pipeline_core import PipelineContext, Stage

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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Always depends on dataframe_loading
        deps = {"dataframe_loading"}
        # Note: Should run after variant_analysis but we handle this through
        # pipeline ordering rather than hard dependency, so the stage can be
        # used independently in tests or other contexts
        return deps

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Should run after variant_analysis because it creates a new DataFrame
        return {"variant_analysis", "chunked_analysis"}

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

        # Use VAR_ID as the variant identifier column
        id_column = "VAR_ID"
        if id_column in df.columns:
            logger.debug(f"Variant ID column '{id_column}' already exists")
            return context

        logger.info("Generating variant identifiers")

        # Generate IDs based on key fields with hash
        key_fields = ["CHROM", "POS", "REF", "ALT"]
        if all(field in df.columns for field in key_fields):
            import hashlib

            var_ids = []
            for idx, r in enumerate(df.itertuples(index=False), 1):
                key = (
                    f"{getattr(r, 'CHROM', '')}{getattr(r, 'POS', '')}"
                    f"{getattr(r, 'REF', '')}{getattr(r, 'ALT', '')}"
                )
                h = hashlib.md5(key.encode("utf-8")).hexdigest()[:4]
                var_ids.append(f"var_{idx:04d}_{h}")
            df.insert(0, id_column, var_ids)
        else:
            # Fallback to simple index
            df.insert(0, id_column, [f"var_{i:04d}_0000" for i in range(1, len(df) + 1)])

        annotation_config = normalize_annotation_config(context.config.copy())

        # Only add Custom_Annotation column if custom annotations are actually requested
        custom_annotations_requested = any(
            [
                annotation_config.get("annotate_bed_files", []),
                annotation_config.get("annotate_gene_lists", []),
                annotation_config.get("annotate_json_genes", []),
            ]
        )

        if custom_annotations_requested and "Custom_Annotation" not in df.columns:
            # Find a good position for it - after GT column if it exists
            if "GT" in df.columns:
                gt_pos = int(cast(Any, df.columns.get_loc("GT"))) + 1
                df.insert(gt_pos, "Custom_Annotation", "")
            else:
                df["Custom_Annotation"] = ""

        context.current_dataframe = df
        return context

    def _process_streaming(self, context: PipelineContext) -> PipelineContext:
        """Process large files in streaming mode without loading into memory."""
        input_file = Path(context.data)
        output_file = context.workspace.get_intermediate_path("with_variant_ids.tsv.gz")
        id_column = "VAR_ID"  # Default variant identifier column

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
                    combined = "".join(
                        fields[key_indices[k]] for k in ["CHROM", "POS", "REF", "ALT"]
                    )
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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depends on dataframe_loading since filtering can work on any columns
        return {"dataframe_loading"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - only filters DataFrame, no external I/O

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Final filtering often filters on computed columns from these stages
        return {
            "variant_analysis",
            "variant_identifier",
            "variant_scoring",
            "inheritance_analysis",
            "custom_annotation",
            "clinvar_pm5",
        }

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply final filtering."""
        from . import output_stages as _os_mod

        # Check for late filtering or final filter
        late_filtering = context.config.get("late_filtering")
        final_filter = context.config.get("final_filter")
        original_filter = (
            (context.config.get("filters") or context.config.get("filter"))
            if late_filtering
            else None
        )

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
                f"FinalFilteringStage - Sample nephro_candidate_score values: {sample_values}"
            )
        else:
            logger.debug("FinalFilteringStage - nephro_candidate_score column NOT found")

        # Apply late filtering if configured
        if late_filtering and original_filter:
            logger.info(f"Applying late filter: {original_filter}")
            df = _os_mod.filter_dataframe_with_query(df, original_filter)
            logger.info(f"Late filter retained {len(df)}/{initial_count} variants")

        # Apply final filter if configured
        if final_filter:
            logger.info(f"Applying final filter: {final_filter}")
            pre_filter_count = len(df)
            context.pre_final_filter_variants = pre_filter_count
            context.config["pre_final_filter_variants"] = pre_filter_count
            df = _os_mod.filter_dataframe_with_query(df, final_filter)
            post_filter_count = len(df)
            logger.info(
                f"Final filter retained {post_filter_count}/{pre_filter_count} "
                f"variants (total: {post_filter_count}/{initial_count})"
            )

        context.current_dataframe = df
        return context


def _tsv_has_data_rows(tsv_path: Path | str) -> bool:
    """Check if TSV file contains at least one data row beyond the header."""
    open_fn = gzip.open if str(tsv_path).endswith(".gz") else open
    try:
        with open_fn(tsv_path, "rt", encoding="utf-8", errors="replace") as f:
            header = f.readline()
            if not header:
                return False
            return any(line.strip() for line in f)
    except Exception:
        return False


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
    def dependencies(self) -> set[str]:
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
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Soft dependencies are optional - if these stages exist, we should run after them
        # This includes all stages that might modify the DataFrame
        return {"variant_scoring", "final_filtering", "pseudonymization", "chunked_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Write final TSV output."""
        df = context.current_dataframe

        # Handle chunked processing
        if context.config.get("use_chunked_processing") and df is None:
            # Check if chunked analysis results are available
            if hasattr(context, "chunked_analysis_tsv") and context.chunked_analysis_tsv:
                if not Path(context.chunked_analysis_tsv).exists():
                    raise RuntimeError("No data to write")
                if context.config.get("fail_on_empty", False) and not _tsv_has_data_rows(
                    context.chunked_analysis_tsv
                ):
                    logger.error("Zero variants remaining in final output (--fail-on-empty)")
                    sys.exit(2)
                logger.info(f"Using chunked analysis results from: {context.chunked_analysis_tsv}")
                # Copy chunked results to final output location
                output_file = context.config.get("output_file")
                if not output_file or output_file in ["stdout", "-"]:
                    open_fn = (
                        gzip.open if str(context.chunked_analysis_tsv).endswith(".gz") else open
                    )
                    with open_fn(context.chunked_analysis_tsv, "rt") as f:
                        for line in f:
                            sys.stdout.write(line)
                    return context
                final_output_path = Path(output_file)
                if not final_output_path.is_absolute():
                    final_output_path = context.workspace.output_dir / final_output_path
                src_gz = str(context.chunked_analysis_tsv).endswith(".gz")
                dst_gz = str(final_output_path).endswith(".gz")
                with AtomicFileOperation(final_output_path) as temp_out:
                    if src_gz == dst_gz:
                        shutil.copy2(context.chunked_analysis_tsv, temp_out)
                    elif src_gz:
                        with (
                            gzip.open(context.chunked_analysis_tsv, "rb") as fi,
                            open(temp_out, "wb") as fo,
                        ):
                            shutil.copyfileobj(fi, fo)
                    else:
                        with (
                            open(context.chunked_analysis_tsv, "rb") as fi,
                            gzip.open(temp_out, "wb") as fo,
                        ):
                            shutil.copyfileobj(fi, fo)
                context.final_output_path = final_output_path
                logger.info(f"Chunked analysis results copied to: {final_output_path}")
                return context

            # Fallback: Output already written during chunked processing
            output_file = context.config.get("chunked_output_file")
            if output_file:
                # Ensure the path is within the output directory
                output_filename = Path(output_file).name
                context.final_output_path = context.workspace.output_dir / output_filename
                logger.info(f"Chunked processing output: {context.final_output_path}")
                if context.config.get("fail_on_empty", False) and not _tsv_has_data_rows(
                    context.final_output_path
                ):
                    logger.error("Zero variants remaining in final output (--fail-on-empty)")
                    sys.exit(2)
                return context
            raise RuntimeError("No data to write")

        if df is None:
            raise RuntimeError("No data to write")

        if len(df) == 0 and context.config.get("fail_on_empty", False):
            logger.error("Zero variants remaining in final output (--fail-on-empty)")
            sys.exit(2)

        # Phase 11: Keep per-sample GT columns available for analysis, but collapse
        # them out of final materialized outputs unless users requested raw output.
        extra_sample_fields = (
            context.config.get("extra_sample_fields")
            if context.config.get(
                "append_extra_sample_fields", bool(context.config.get("extra_sample_fields"))
            )
            else None
        )
        df = _collapse_per_sample_gt_columns_for_output(
            df,
            context.vcf_samples,
            no_replacement=context.config.get("no_replacement", False),
            extra_sample_fields=extra_sample_fields,
            delimiter=context.config.get("extra_sample_field_delimiter") or ":",
            sample_mapping=context.config.get("pseudonymization_mapping"),
        )
        context.current_dataframe = df

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
            link_configs = resolve_links(context.config)
            if link_configs:
                import io

                lines = ["\t".join(df.columns)]
                lines.extend(
                    "\t".join(str(val) for val in row) for row in df.itertuples(index=False)
                )
                df = pd.read_csv(
                    io.StringIO("\n".join(add_links_to_table(lines, link_configs))),
                    sep="\t",
                    dtype=str,
                )

        # Drop internal cache columns before output
        cache_cols = [c for c in df.columns if c.startswith("_")]
        if cache_cols:
            df = df.drop(columns=cache_cols)
            logger.debug(f"Dropped {len(cache_cols)} cache columns: {cache_cols}")

        # Restore original column names before writing output
        if context.column_rename_map:
            reverse_map = {v: k for k, v in context.column_rename_map.items()}
            df = df.rename(columns=reverse_map)
            logger.debug(f"Restored {len(reverse_map)} column names for output")

        # Write file
        compression = None
        if context.config.get("gzip_output") or str(output_path).endswith(".gz"):
            compression = "gzip"
            if not str(output_path).endswith(".gz"):
                output_path = Path(str(output_path) + ".gz")

        with AtomicFileOperation(output_path) as temp_out:
            df.to_csv(
                temp_out, sep="\t", index=False, na_rep="", compression=cast(Any, compression)
            )

        logger.info(f"Successfully wrote output to: {output_path}")

        context.final_output_path = output_path
        context.report_paths["tsv"] = output_path

        return context

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the final TSV output for checkpoint tracking."""
        return [context.final_output_path] if context.final_output_path else []


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
    def dependencies(self) -> set[str]:
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


from .metadata_helpers import collect_metadata_rows as _collect_metadata_rows  # noqa: F401
