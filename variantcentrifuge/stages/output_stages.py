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
from pathlib import Path
from typing import Any

import pandas as pd

from ..checkpoint_models import AtomicFileOperation
from ..converter import (
    convert_to_excel,
    produce_report_json,
)
from ..filters import filter_dataframe_with_query
from ..generate_html_report import generate_html_report
from ..generate_igv_report import generate_igv_report
from ..gt_reconstruction import (
    _collapse_per_sample_gt_columns_for_output,
    _find_per_sample_gt_columns,
    reconstruct_gt_column,
)
from ..pipeline_core import PipelineContext, Stage
from ..pseudonymizer import create_pseudonymizer
from ..run_summary import _collect_metadata_rows, write_run_summary
from ..utils import get_tool_version, sanitize_metadata_field
from .tsv_output_stages import (
    ArchiveCreationStage,
    FinalFilteringStage,
    TSVOutputStage,
    VariantIdentifierStage,
)

logger = logging.getLogger(__name__)

__all__ = [
    "ArchiveCreationStage",
    "ExcelReportStage",
    "FinalFilteringStage",
    "HTMLReportStage",
    "IGVReportStage",
    "MetadataGenerationStage",
    "PseudonymizationStage",
    "TSVOutputStage",
    "VariantIdentifierStage",
    "_collapse_per_sample_gt_columns_for_output",
    "_find_per_sample_gt_columns",
    "filter_dataframe_with_query",
    "reconstruct_gt_column",
]


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
    def dependencies(self) -> set[str]:
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
            with open(metadata_file) as f:
                metadata_dict = json.load(f)

        mapping = pseudonymizer.create_mapping(samples, metadata_dict)

        # Apply to DataFrame
        df = pseudonymizer.pseudonymize_dataframe(df)

        # Save mapping to user-specified table path
        mapping_file = context.config.get("pseudonymize_table")
        if mapping_file:
            # Use the pseudonymizer's save_mapping method for proper TSV format
            pseudonymizer.save_mapping(mapping_file, include_metadata=True)
            logger.info(f"Saved pseudonymization mapping to {mapping_file}")
        else:
            # Fallback (shouldn't happen due to CLI validation)
            mapping_file = (
                context.workspace.output_dir.parent
                / f"pseudonymization_mapping_{context.workspace.timestamp}.tsv"
            )
            pseudonymizer.save_mapping(str(mapping_file), include_metadata=True)
            logger.warning(f"No pseudonymize_table specified, saved to {mapping_file}")

        context.current_dataframe = df
        context.vcf_samples = [mapping.get(s, s) for s in samples]
        context.config["pseudonymization_mapping"] = mapping
        context.config["pseudonymizer"] = pseudonymizer

        # Create pseudonymized PED file if requested
        if context.config.get("pseudonymize_ped") and context.config.get("ped_file"):
            # Create PED filename based on mapping file path
            mapping_path = Path(mapping_file)
            ped_output = mapping_path.parent / (mapping_path.stem + "_pedigree.ped")
            pseudonymizer.pseudonymize_ped_file(context.config["ped_file"], str(ped_output))
            logger.info(f"Pseudonymized PED file saved to: {ped_output}")

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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"tsv_output", "metadata_generation", "statistics_generation"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of optional stage names this stage depends on."""
        # This ensures Excel generation waits for gene burden/association analysis if running.
        return {"gene_burden_analysis", "association_analysis"}

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
            raise RuntimeError("No input file for Excel generation")

        output_path = context.workspace.get_output_path("", ".xlsx")

        logger.info(f"Generating Excel report: {output_path}")

        # Use in-memory DataFrame if available (avoids redundant disk read)
        # Prefer current_dataframe (has all columns from scoring/analysis stages)
        # over variants_df (stale snapshot from dataframe_loading — see #80)
        excel_df = None
        source_df = (
            context.current_dataframe
            if context.current_dataframe is not None
            else context.variants_df
        )
        if source_df is not None:
            excel_df = source_df.copy()
            extra_sample_fields = (
                context.config.get("extra_sample_fields")
                if context.config.get(
                    "append_extra_sample_fields", bool(context.config.get("extra_sample_fields"))
                )
                else None
            )
            excel_df = _collapse_per_sample_gt_columns_for_output(
                excel_df,
                context.vcf_samples,
                no_replacement=context.config.get("no_replacement", False),
                extra_sample_fields=extra_sample_fields,
                delimiter=context.config.get("extra_sample_field_delimiter") or ":",
                sample_mapping=context.config.get("pseudonymization_mapping"),
            )
            # Drop internal cache columns before Excel output
            cache_cols = [c for c in excel_df.columns if c.startswith("_")]
            if cache_cols:
                excel_df = excel_df.drop(columns=cache_cols)
                logger.debug(f"Dropped {len(cache_cols)} cache columns: {cache_cols}")
            # Restore original column names for Excel output
            if context.column_rename_map:
                reverse_map = {v: k for k, v in context.column_rename_map.items()}
                excel_df = excel_df.rename(columns=reverse_map)
                logger.debug(f"Restored {len(reverse_map)} column names for Excel generation")
            logger.info("Using in-memory DataFrame for Excel generation (skipping disk read)")

        # Convert TSV to Excel using the correct API
        # The convert_to_excel function automatically creates the output path
        # by changing the .tsv extension to .xlsx
        # Collect additional sheets (Metadata, Statistics, Gene Burden, Association)
        extra_sheets = self._collect_additional_sheets(context)

        # Convert to Excel in a single XlsxWriter session (all sheets + formatting)
        xlsx_file = convert_to_excel(
            str(input_file),
            context.config,
            df=excel_df,
            extra_sheets=extra_sheets,
            output_file=str(output_path),
        )

        # Move the generated file to the desired location if needed
        if xlsx_file != str(output_path):
            import shutil

            shutil.move(xlsx_file, str(output_path))
            xlsx_file = str(output_path)

        logger.info(f"Excel report generated with all sheets: {xlsx_file}")

        context.report_paths["excel"] = output_path
        return context

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the Excel report for checkpoint tracking."""
        path = context.report_paths.get("excel")
        return [path] if path else []

    def _collect_additional_sheets(self, context: PipelineContext) -> dict[str, pd.DataFrame]:
        """Collect available optional sheets as DataFrames for single-session Excel write."""
        sheets: dict[str, pd.DataFrame] = {}

        # Add Metadata sheet
        metadata_file = context.report_paths.get("metadata")
        if (
            metadata_file
            and Path(metadata_file).exists()
            and Path(metadata_file).stat().st_size > 0
        ):
            try:
                sheets["Metadata"] = pd.read_csv(
                    metadata_file, sep="\t", header=0, low_memory=False
                )
                logger.info("Successfully collected Metadata sheet")
            except Exception as e:
                logger.error(f"Failed to read Metadata sheet: {e}")

        # Add Statistics sheet
        if not context.config.get("no_stats"):
            stats_file = context.config.get("stats_output_file")
            if stats_file and Path(stats_file).exists() and Path(stats_file).stat().st_size > 0:
                try:
                    sheets["Statistics"] = pd.read_csv(
                        stats_file, sep="\t", header=0, low_memory=False
                    )
                    logger.info("Successfully collected Statistics sheet")
                except Exception as e:
                    logger.error(f"Failed to read Statistics sheet: {e}")

        # Add Gene Burden sheet
        if context.config.get("perform_gene_burden"):
            gb_file = context.config.get("gene_burden_output")
            if gb_file and Path(gb_file).exists() and Path(gb_file).stat().st_size > 0:
                try:
                    sheets["Gene Burden"] = pd.read_csv(
                        gb_file, sep="\t", header=0, low_memory=False
                    )
                    logger.info("Successfully collected Gene Burden sheet")
                except Exception as e:
                    logger.error(f"Failed to read Gene Burden sheet: {e}")

        # Add Association sheet
        if context.config.get("perform_association"):
            assoc_file = context.config.get("association_output")
            if assoc_file and Path(assoc_file).exists() and Path(assoc_file).stat().st_size > 0:
                try:
                    sheets["Association"] = pd.read_csv(
                        assoc_file, sep="\t", header=0, low_memory=False
                    )
                    logger.info("Successfully collected Association sheet")
                except Exception as e:
                    logger.error(f"Failed to read Association sheet: {e}")

        return sheets

    def _add_additional_sheets(self, xlsx_file: str, context: PipelineContext) -> None:
        """Add additional sheets to existing Excel file (kept for backward compatibility)."""
        extra_sheets = self._collect_additional_sheets(context)
        import shutil

        with AtomicFileOperation(xlsx_file) as temp_xlsx:
            shutil.copy2(xlsx_file, temp_xlsx)
            with pd.ExcelWriter(
                str(temp_xlsx), engine="openpyxl", mode="a", if_sheet_exists="new"
            ) as writer:
                for sheet_name, sheet_df in extra_sheets.items():
                    sheet_df.to_excel(writer, index=False, sheet_name=sheet_name)


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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"tsv_output"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # HTML report should wait for IGV reports if they're being generated
        # so that IGV links can be included in the HTML report
        return {"igv_report"}

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

        try:
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
        except Exception as e:
            logger.error(f"HTML report generation failed: {e}")

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
    def dependencies(self) -> set[str]:
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
            max_workers = min(
                context.config.get("threads", 4), 8
            )  # Cap at 8 to avoid overwhelming system

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
                igv_max_variant_part_filename=context.config.get(
                    "igv_max_variant_part_filename", 50
                ),
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
    def dependencies(self) -> set[str]:
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

        self._generate_metadata(context)
        return context

    def _generate_metadata(
        self, context: PipelineContext, metadata_path: Path | None = None
    ) -> Path:
        """Generate metadata file with run provenance."""
        # Create TSV metadata file
        metadata_path = metadata_path or context.workspace.get_output_path("_metadata", ".tsv")

        logger.info(f"Generating metadata: {metadata_path}")

        # Write metadata in TSV format atomically
        with (
            AtomicFileOperation(metadata_path) as temp_mf,
            open(temp_mf, "w", encoding="utf-8") as mf,
        ):
            mf.write("Parameter\tValue\n")

            def meta_write(param: Any, val: Any) -> None:
                p = sanitize_metadata_field(param)
                v = sanitize_metadata_field(val)
                mf.write(f"{p}\t{v}\n")

            for param, val in _collect_metadata_rows(context, get_tool_version):
                meta_write(param, val)

        context.report_paths["metadata"] = metadata_path

        # Consolidate provenance manifest in run_summary.json
        summary_path = write_run_summary(context)
        context.report_paths["run_summary"] = summary_path

        return metadata_path

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the metadata report and run summary for checkpoint tracking."""
        paths = []
        if path := context.report_paths.get("metadata"):
            paths.append(path)
        if summary := context.report_paths.get("run_summary"):
            paths.append(summary)
        return paths
