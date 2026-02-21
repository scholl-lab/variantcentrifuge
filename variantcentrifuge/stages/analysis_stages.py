"""
Analysis stages for variant processing and computation.

This module contains stages that perform analysis on the extracted data:
- DataFrame loading and chunked processing
- Custom annotations (BED, gene lists, JSON)
- Inheritance pattern analysis
- Variant scoring
- Statistics generation
- Gene burden analysis
"""

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ..analyze_variants import analyze_variants
from ..annotator import annotate_dataframe_with_features, load_custom_features
from ..association.base import AssociationConfig
from ..association.engine import AssociationEngine
from ..dataframe_optimizer import load_optimized_dataframe, should_use_memory_passthrough
from ..filters import filter_final_tsv_by_genotype
from ..gene_burden import (
    _aggregate_gene_burden_from_columns,
    _aggregate_gene_burden_from_gt,
    _aggregate_gene_burden_legacy,
    _find_gt_columns,
    perform_gene_burden_analysis,
)
from ..inheritance import analyze_inheritance
from ..pipeline_core import PipelineContext, Stage
from ..scoring import apply_scoring
from ..stats_engine import StatsEngine

logger = logging.getLogger(__name__)

# Standard column aliases: sanitized name -> canonical short name.
# When ANN[0].GENE is sanitized to ANN_0__GENE, downstream stages still expect "GENE".
_STANDARD_COLUMN_ALIASES = {
    "ANN_0__GENE": "GENE",
    "ANN_0__EFFECT": "EFFECT",
    "ANN_0__IMPACT": "IMPACT",
    "ANN_0__HGVS_C": "HGVS_C",
    "ANN_0__HGVS_P": "HGVS_P",
}


def _create_standard_column_aliases(df: pd.DataFrame) -> None:
    """Create standard short-name aliases for sanitized ANN columns.

    After column sanitization (ANN[0].GENE -> ANN_0__GENE), many downstream
    stages expect plain column names like GENE, EFFECT, IMPACT. This function
    creates those aliases in-place if the sanitized source column exists and
    the short name doesn't already exist.
    """
    for sanitized_name, short_name in _STANDARD_COLUMN_ALIASES.items():
        if short_name not in df.columns and sanitized_name in df.columns:
            df[short_name] = df[sanitized_name]
            logger.debug(f"Created column alias: {sanitized_name} -> {short_name}")


def create_sample_columns_from_gt_vectorized(
    df: pd.DataFrame, vcf_samples: list[str], separator: str = ";", snpsift_sep: str = ","
) -> pd.DataFrame:
    """Create individual sample columns from GT column data using vectorized operations.

    This function is a high-performance version of create_sample_columns_from_gt that uses
    pandas vectorized string operations instead of row-by-row regex parsing. Provides
    5-10x performance improvement for large datasets.

    Args
    ----
        df: DataFrame containing GT column
        vcf_samples: List of sample IDs to create columns for
        separator: Separator used in replaced genotype format (default: ";")
        snpsift_sep: Separator used in SnpSift format (default: ",")

    Returns
    -------
        DataFrame with individual sample columns added

    Raises
    ------
        ValueError: If GT column is missing or empty
    """
    # Phase 11: Check if sample columns already exist FIRST (before checking GT)
    # This allows bcftools query output (per-sample columns, no GT) to skip creation
    sample_columns_exist = any(sample_id in df.columns for sample_id in vcf_samples)

    if sample_columns_exist:
        logger.debug("Sample columns already exist, skipping creation")
        return df

    # If sample columns don't exist, we need GT column to create them
    if "GT" not in df.columns:
        raise ValueError("GT column not found in DataFrame")

    if len(df) == 0:
        logger.warning("Empty DataFrame provided to create_sample_columns_from_gt_vectorized")
        return df

    # Get the first GT value to determine format
    first_gt = str(df.iloc[0]["GT"]) if len(df) > 0 else ""

    # Check if genotype replacement has already been done
    # (format: "Sample1(0/1);Sample2(0/0)")
    if "(" in first_gt and ")" in first_gt:
        logger.debug(
            "GT column contains replaced genotypes, extracting sample genotypes using "
            "vectorized operations"
        )

        # Create a copy to avoid modifying original
        df_copy = df.copy()

        # Initialize all sample columns with missing genotypes efficiently
        # Create all sample columns at once to avoid DataFrame fragmentation
        sample_data = dict.fromkeys(vcf_samples, "./.")
        sample_columns = pd.DataFrame(sample_data, index=df_copy.index)
        df_copy = pd.concat([df_copy, sample_columns], axis=1)

        # Handle missing/null GT values
        gt_series = df_copy["GT"].fillna("").astype(str)

        # Process only non-empty GT entries
        valid_mask = (gt_series != "") & (gt_series != "NA") & (gt_series != "nan")
        valid_gt = gt_series[valid_mask]

        if len(valid_gt) > 0:
            # Split each GT entry by separator and explode to individual sample entries
            # This creates a long-format DataFrame with one row per sample entry
            gt_split = valid_gt.str.split(separator).explode().reset_index()
            gt_split.columns = ["original_index", "sample_entry"]

            # Filter out empty entries
            gt_split = gt_split[gt_split["sample_entry"].str.strip() != ""]

            if len(gt_split) > 0:
                # Extract sample names and genotypes using vectorized regex
                # Pattern: "SampleName(genotype:extra:fields)" -> extract sample name and genotype
                extraction = gt_split["sample_entry"].str.extract(r"^([^(]+)\(([^:)]+)")
                gt_split["sample_name"] = extraction[0].str.strip()
                gt_split["genotype"] = extraction[1]

                # Remove entries where extraction failed
                valid_extractions = gt_split.dropna(subset=["sample_name", "genotype"])

                # Process each sample
                for sample_id in vcf_samples:
                    # Find entries for this sample
                    sample_mask = valid_extractions["sample_name"] == sample_id
                    sample_entries = valid_extractions[sample_mask]

                    if len(sample_entries) > 0:
                        # Map genotypes back to original dataframe indices
                        for row in sample_entries.itertuples(index=False):
                            original_idx = getattr(row, "original_index", None)
                            genotype = getattr(row, "genotype", "")
                            df_copy.loc[original_idx, sample_id] = genotype

        logger.debug(
            f"Created {len(vcf_samples)} sample columns from replaced genotypes using "
            f"vectorized operations"
        )
        return df_copy

    elif snpsift_sep in first_gt:
        logger.debug(
            f"Extracting sample genotypes from GT column using separator '{snpsift_sep}' "
            f"with vectorized operations"
        )

        # Create a copy to avoid modifying original
        df_copy = df.copy()

        # Initialize all sample columns with missing genotypes efficiently
        # Create all sample columns at once to avoid DataFrame fragmentation
        sample_data = dict.fromkeys(vcf_samples, "./.")
        sample_columns = pd.DataFrame(sample_data, index=df_copy.index)
        df_copy = pd.concat([df_copy, sample_columns], axis=1)

        # Handle missing/null GT values
        gt_series = df_copy["GT"].fillna("").astype(str)

        # Process only non-empty GT entries
        valid_mask = (gt_series != "") & (gt_series != "NA") & (gt_series != "nan")
        valid_gt = gt_series[valid_mask]

        if len(valid_gt) > 0:
            # Split GT values by SnpSift separator
            gt_split = valid_gt.str.split(snpsift_sep, expand=True)

            # Map each column to corresponding sample
            num_samples_in_data = min(len(vcf_samples), gt_split.shape[1])

            for i in range(num_samples_in_data):
                sample_id = vcf_samples[i]
                genotypes = gt_split.iloc[:, i].fillna("./.")

                # Map back to original dataframe indices
                df_copy.loc[valid_gt.index, sample_id] = genotypes.values

        logger.debug(
            f"Created {len(vcf_samples)} sample columns for inheritance analysis using "
            f"vectorized operations"
        )
        return df_copy

    else:
        logger.warning(
            f"GT column format not recognized: {first_gt[:100]}. "
            f"Inheritance analysis may not work correctly."
        )
        return df


def create_sample_columns_from_gt_intelligent(
    df: pd.DataFrame,
    vcf_samples: list[str],
    separator: str = ";",
    snpsift_sep: str = ",",
    method: str = "auto",
) -> pd.DataFrame:
    """Intelligent sample column creation with automatic method selection.

    This function automatically selects between iterative and vectorized methods
    based on dataset characteristics and user preferences.

    Args
    ----
        df: DataFrame containing GT column
        vcf_samples: List of sample IDs to create columns for
        separator: Separator used in replaced genotype format (default: ";")
        snpsift_sep: Separator used in SnpSift format (default: ",")
        method: Creation method ("auto", "iterative", "vectorized")

    Returns
    -------
        DataFrame with individual sample columns added
    """
    import time

    # Auto-select method based on dataset characteristics
    if method == "auto":
        num_variants = len(df)
        num_samples = len(vcf_samples)
        total_operations = num_variants * num_samples

        # Use vectorized method for large datasets (significant performance gain)
        # Threshold: > 1M operations (e.g., 1000 variants x 1000 samples)
        if total_operations > 1_000_000:
            selected_method = "vectorized"
            logger.info(
                f"Auto-selected vectorized sample column creation for {num_variants:,} variants x "
                f"{num_samples:,} samples ({total_operations:,} operations)"
            )
        else:
            selected_method = "iterative"
            logger.debug(
                f"Auto-selected iterative sample column creation for {num_variants:,} variants x "
                f"{num_samples:,} samples ({total_operations:,} operations)"
            )
    else:
        selected_method = method
        logger.debug(f"Using user-specified sample column creation method: {selected_method}")

    # Execute selected method with timing
    start_time = time.time()

    if selected_method == "vectorized":
        result_df = create_sample_columns_from_gt_vectorized(
            df, vcf_samples, separator, snpsift_sep
        )
    else:  # iterative
        result_df = create_sample_columns_from_gt(df, vcf_samples, separator, snpsift_sep)

    end_time = time.time()
    processing_time = end_time - start_time

    # Log performance information
    if processing_time > 1.0:  # Only log if significant
        num_variants = len(df)
        num_samples = len(vcf_samples)
        ops_per_second = (
            (num_variants * num_samples) / processing_time if processing_time > 0 else 0
        )
        logger.info(
            f"Sample column creation ({selected_method}): {processing_time:.2f}s "
            f"for {num_variants:,} variants x {num_samples:,} samples "
            f"({ops_per_second:,.0f} ops/sec)"
        )

    return result_df


def create_sample_columns_from_gt(
    df: pd.DataFrame, vcf_samples: list[str], separator: str = ";", snpsift_sep: str = ","
) -> pd.DataFrame:
    """Create individual sample columns from GT column data.

    This function extracts genotype information from the GT column and creates
    individual sample columns for inheritance analysis. It handles both replaced
    genotype format (e.g., "Sample1(0/1);Sample2(0/0)") and SnpSift format
    (e.g., "0/1,0/0,1/1").

    Args
    ----
        df: DataFrame containing GT column
        vcf_samples: List of sample IDs to create columns for
        separator: Separator used in replaced genotype format (default: ";")
        snpsift_sep: Separator used in SnpSift format (default: ",")

    Returns
    -------
        DataFrame with individual sample columns added

    Raises
    ------
        ValueError: If GT column is missing or empty
    """
    # Phase 11: Check if sample columns already exist FIRST (before checking GT)
    # This allows bcftools query output (per-sample columns, no GT) to skip creation
    sample_columns_exist = any(sample_id in df.columns for sample_id in vcf_samples)

    if sample_columns_exist:
        logger.debug("Sample columns already exist, skipping creation")
        return df

    # If sample columns don't exist, we need GT column to create them
    if "GT" not in df.columns:
        raise ValueError("GT column not found in DataFrame")

    if len(df) == 0:
        logger.warning("Empty DataFrame provided to create_sample_columns_from_gt")
        return df

    # Get the first GT value to determine format
    first_gt = str(df.iloc[0]["GT"]) if len(df) > 0 else ""

    # Check if genotype replacement has already been done
    # (format: "Sample1(0/1);Sample2(0/0)")
    if "(" in first_gt and ")" in first_gt:
        logger.debug("GT column contains replaced genotypes, extracting sample genotypes")

        # Pre-create all sample columns data
        sample_data: dict[str, list[str]] = {sample_id: [] for sample_id in vcf_samples}

        # Parse replaced genotype format
        for row in df.itertuples(index=False):
            gt_value = str(getattr(row, "GT", ""))
            # Initialize all samples with missing genotype
            row_genotypes = dict.fromkeys(vcf_samples, "./.")

            if gt_value and gt_value != "NA" and gt_value != "nan":
                # Split by separator (usually semicolon) for different samples
                sample_entries = gt_value.split(separator)
                for entry in sample_entries:
                    # Parse "SampleName(genotype:extra:fields)"
                    match = re.match(r"^([^(]+)\(([^)]+)\)", entry)
                    if match:
                        sample_name = match.group(1)
                        genotype_info = match.group(2)
                        # Extract just the genotype part (before the first colon)
                        genotype = genotype_info.split(":")[0]
                        if sample_name in row_genotypes:
                            row_genotypes[sample_name] = genotype

            # Add genotypes for all samples
            for sample_id in vcf_samples:
                sample_data[sample_id].append(row_genotypes[sample_id])

        # Create DataFrame from sample data and concatenate
        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        logger.debug(f"Created {len(sample_data)} sample columns from replaced genotypes")

    elif snpsift_sep in first_gt:
        logger.debug(f"Extracting sample genotypes from GT column using separator '{snpsift_sep}'")

        # Pre-create all sample columns data
        sample_data = {sample_id: [] for sample_id in vcf_samples}

        # Extract genotypes for each row
        for row in df.itertuples(index=True):
            gt_value = str(getattr(row, "GT", ""))
            if gt_value and gt_value != "NA" and gt_value != "nan":
                genotypes = gt_value.split(snpsift_sep)
                if len(genotypes) != len(vcf_samples):
                    logger.warning(
                        f"Row {row.Index}: Expected {len(vcf_samples)} genotypes "
                        f"but found {len(genotypes)}"
                    )
                for i, sample_id in enumerate(vcf_samples):
                    if i < len(genotypes):
                        sample_data[sample_id].append(genotypes[i])
                    else:
                        sample_data[sample_id].append("./.")
            else:
                # Missing GT data
                for sample_id in vcf_samples:
                    sample_data[sample_id].append("./.")

        # Create DataFrame from sample data and concatenate
        sample_df = pd.DataFrame(sample_data)
        df = pd.concat([df, sample_df], axis=1)

        logger.debug(f"Created {len(sample_data)} sample columns for inheritance analysis")

    else:
        logger.warning(
            f"GT column format not recognized: {first_gt[:100]}. "
            f"Inheritance analysis may not work correctly."
        )

    return df


def handle_inheritance_analysis_error(
    df: pd.DataFrame,
    error: Exception,
    preserve_details_for_scoring: bool = False,
    context_description: str = "inheritance analysis",
) -> pd.DataFrame:
    """Handle errors in inheritance analysis with consistent error recovery.

    This function provides unified error handling for both chunked and non-chunked
    processing. It ensures the DataFrame schema is maintained even when errors occur.

    Args
    ----
        df: DataFrame that encountered an error
        error: The exception that occurred
        preserve_details_for_scoring: Whether to add Inheritance_Details column
        context_description: Description of the processing context for logging

    Returns
    -------
        DataFrame with error-state inheritance columns added

    Note
    ----
        - Logs the error for debugging
        - Adds Inheritance_Pattern column with "error" value
        - Optionally adds Inheritance_Details column if needed for scoring
        - Maintains DataFrame schema consistency
    """
    logger.error(f"Error in {context_description}: {error}")

    # Add error-state inheritance columns to maintain schema
    df = df.copy()
    df["Inheritance_Pattern"] = "error"

    if preserve_details_for_scoring:
        df["Inheritance_Details"] = "{}"

    return df


def cleanup_sample_columns(
    df: pd.DataFrame, vcf_samples: list[str], preserve_columns: list[str] | None = None
) -> pd.DataFrame:
    """Clean up sample columns after inheritance analysis.

    This function removes temporary sample columns that were created for inheritance
    analysis while preserving essential columns. It provides a unified approach
    for both chunked and non-chunked processing.

    Args
    ----
        df: DataFrame containing sample columns to clean up
        vcf_samples: List of sample IDs that should be removed
        preserve_columns: List of column names to preserve even if they match sample names
                         (default: ["GT", "Inheritance_Pattern", "Inheritance_Details"])

    Returns
    -------
        DataFrame with sample columns cleaned up

    Note
    ----
        - Removes columns that match VCF sample names
        - Preserves specified columns even if they match sample names
        - Logs cleanup statistics for debugging
    """
    if not vcf_samples:
        return df

    if preserve_columns is None:
        preserve_columns = ["GT", "Inheritance_Pattern", "Inheritance_Details"]

    # Find sample columns to remove
    sample_columns_to_remove = [
        col for col in df.columns if col in vcf_samples and col not in preserve_columns
    ]

    if sample_columns_to_remove:
        logger.debug(
            f"Removing {len(sample_columns_to_remove)} sample columns after inheritance analysis"
        )
        df = df.drop(columns=sample_columns_to_remove)
        logger.debug(f"Columns after cleanup: {len(df.columns)}")
    else:
        logger.debug("No sample columns to remove")

    return df


class DataFrameLoadingStage(Stage):
    """Load TSV data into DataFrame or prepare for chunked processing."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "dataframe_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load data into DataFrame for analysis"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # In sequential mode, depends on field_extraction
        # In parallel mode, ParallelCompleteProcessingStage handles field extraction internally
        # We return an empty set here and check for data availability in _process()
        # This allows the stage to work in both sequential and parallel modes
        return set()

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after all data transformation stages
        # In sequential mode: after field_extraction and processing stages
        # In parallel mode: after parallel_complete_processing
        return {
            "field_extraction",
            "parallel_complete_processing",
            "genotype_replacement",
            "phenotype_integration",
            "extra_column_removal",
        }

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the DataFrame
        from the most recent TSV file so that subsequent stages can proceed.
        """
        # Find the most recent TSV file using the same logic as _process
        input_file = self._find_input_file(context)

        if not input_file:
            logger.warning("No TSV file found to restore DataFrame during checkpoint skip")
            return context

        # Check if we should use chunked processing - skip if parallel processing done
        if self._should_use_chunks(context, input_file) and not context.config.get(
            "chunked_processing_complete"
        ):
            logger.info("File too large, will use chunked processing (checkpoint skip)")
            context.config["use_chunked_processing"] = True
            # Don't load full DataFrame - chunked stages will handle it
            return context

        # Load full DataFrame with optimizations
        logger.info(f"Restoring DataFrame from {input_file} (checkpoint skip)")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        df, rename_map = load_optimized_dataframe(
            str(input_file),
            sep="\t",
            compression=compression,
            on_bad_lines="warn",  # Warn about problematic lines instead of failing
        )

        context.current_dataframe = df
        context.column_rename_map = rename_map
        logger.info(f"Restored {len(df)} variants into DataFrame after checkpoint skip")
        _create_standard_column_aliases(df)

        # Check if DataFrame is small enough for in-memory pass-through
        if should_use_memory_passthrough(df):
            context.variants_df = df
            logger.info("DataFrame stored for in-memory pass-through (checkpoint skip)")

        return context

    def _find_input_file(self, context: PipelineContext) -> Path:
        """Find the most recent TSV file in the pipeline."""
        # Priority: extra_columns_removed > phenotypes_added > genotype_replaced > extracted
        if hasattr(context, "extra_columns_removed_tsv") and context.extra_columns_removed_tsv:
            return context.extra_columns_removed_tsv
        elif hasattr(context, "phenotypes_added_tsv") and context.phenotypes_added_tsv:
            return context.phenotypes_added_tsv
        elif hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            return context.genotype_replaced_tsv
        elif hasattr(context, "extracted_tsv") and context.extracted_tsv:
            return context.extracted_tsv
        else:
            # Fallback to context.data, but verify it's a TSV file
            input_file = context.data
            if input_file and (
                str(input_file).endswith(".vcf.gz") or str(input_file).endswith(".vcf")
            ):
                logger.warning(
                    f"DataFrameLoadingStage received VCF file instead of TSV: {input_file}"
                )
                # Try to find any TSV file in order of preference
                for attr in [
                    "extracted_tsv",
                    "genotype_replaced_tsv",
                    "phenotypes_added_tsv",
                    "extra_columns_removed_tsv",
                ]:
                    if hasattr(context, attr):
                        tsv_file = getattr(context, attr)
                        if tsv_file and tsv_file.exists():
                            logger.info(f"Using {attr} instead: {tsv_file}")
                            return Path(tsv_file)
                else:
                    raise ValueError(
                        f"DataFrameLoadingStage requires a TSV file, but got VCF: {input_file}. "
                        f"No TSV files found in context."
                    )
            return Path(input_file)

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load data into DataFrame or prepare chunked processing."""
        # Find the input file
        input_file = self._find_input_file(context)

        # Check if we should use chunked processing - skip if parallel processing done
        if self._should_use_chunks(context, input_file) and not context.config.get(
            "chunked_processing_complete"
        ):
            logger.info("File too large, will use chunked processing")
            context.config["use_chunked_processing"] = True
            # Don't load full DataFrame - chunked stages will handle it
            return context

        # Load full DataFrame with optimizations
        logger.info(f"Loading data from {input_file}")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        df, rename_map = load_optimized_dataframe(
            str(input_file),
            sep="\t",
            compression=compression,
            on_bad_lines="warn",  # Warn about problematic lines instead of failing
        )

        context.current_dataframe = df
        context.column_rename_map = rename_map
        logger.info(f"Loaded {len(df)} variants into DataFrame")

        # Create standard column aliases from sanitized ANN/NMD column names.
        # After sanitization, ANN[0].GENE becomes ANN_0__GENE, but many downstream
        # stages expect a plain "GENE" column. Create aliases for common fields.
        _create_standard_column_aliases(df)

        # Check if DataFrame is small enough for in-memory pass-through
        if should_use_memory_passthrough(df):
            context.variants_df = df
            logger.info("DataFrame stored for in-memory pass-through")

        return context

    def _should_use_chunks(self, context: PipelineContext, file_path: Path) -> bool:
        """Determine if chunked processing should be used."""
        # Check if chunked processing is explicitly disabled
        if context.config.get("no_chunked_processing"):
            return False

        # Check if chunked processing is explicitly forced
        if context.config.get("force_chunked_processing"):
            return True

        # TODO(12-02): Remove this check after migrating to ResourceManager auto-detection
        # Check explicit chunking request (legacy config key, removed from CLI in 12-01)
        if context.config.get("chunks"):
            return True

        # Memory-aware chunking decision
        try:
            # Get sample count for memory estimation
            sample_count = len(context.vcf_samples) if context.vcf_samples else 0

            # Estimate variant count from file size (rough approximation)
            file_size_mb = file_path.stat().st_size / (1024 * 1024)
            estimated_variants = max(1, int(file_size_mb * 1000))  # ~1000 variants per MB

            # Calculate memory footprint estimation
            if sample_count > 0:
                # Each genotype cell ~= 8 bytes in memory after processing
                estimated_memory_mb = (estimated_variants * sample_count * 8) / (1024 * 1024)

                # Use adaptive thresholds based on sample count
                if sample_count > 1000:
                    # Large multi-sample VCFs: use chunking for files > 50MB
                    threshold_mb = context.config.get("chunk_threshold_mb", 50)
                    memory_threshold_mb = 500  # 500MB memory limit
                elif sample_count > 100:
                    # Medium multi-sample VCFs: use chunking for files > 100MB
                    threshold_mb = context.config.get("chunk_threshold_mb", 100)
                    memory_threshold_mb = 1000  # 1GB memory limit
                else:
                    # Small sample sets: use original threshold
                    threshold_mb = context.config.get("chunk_threshold_mb", 500)
                    memory_threshold_mb = 2000  # 2GB memory limit

                # Use chunking if file size OR estimated memory exceeds threshold
                if file_size_mb > threshold_mb:
                    logger.info(
                        f"Using chunked processing: file size {file_size_mb:.1f}MB > "
                        f"{threshold_mb}MB threshold"
                    )
                    return True

                if estimated_memory_mb > memory_threshold_mb:
                    logger.info(
                        f"Using chunked processing: estimated memory {estimated_memory_mb:.1f}MB > "
                        f"{memory_threshold_mb}MB threshold"
                    )
                    logger.info(
                        f"Memory estimation: {estimated_variants} variants x {sample_count} samples"
                    )
                    return True
            else:
                # Fallback to file size only
                threshold_mb = context.config.get("chunk_threshold_mb", 500)
                return bool(file_size_mb > threshold_mb)

        except Exception as e:
            logger.debug(f"Error in chunking decision: {e}")
            return False

        return False

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        # Check in priority order for available TSV files
        for attr in [
            "extra_columns_removed_tsv",
            "phenotypes_added_tsv",
            "genotype_replaced_tsv",
            "extracted_tsv",
        ]:
            if hasattr(context, attr):
                tsv_file = getattr(context, attr)
                if tsv_file and tsv_file.exists():
                    return [tsv_file]
        # Fallback to context.data
        if hasattr(context, "data") and context.data:
            return [context.data]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        # DataFrameLoadingStage doesn't produce file outputs, just loads data into memory
        return []


class CustomAnnotationStage(Stage):
    """Apply custom annotations from BED files, gene lists, and JSON."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "custom_annotation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply custom annotations to variants"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Always depends on dataframe being loaded
        # Annotation config is handled in setup stages if needed
        return {"dataframe_loading"}

    def _has_annotation_config(self) -> bool:
        """Check if annotation config stage exists."""
        return True  # Simplified

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply custom annotations."""
        if not context.annotation_configs:
            logger.debug("No custom annotations configured")
            return context

        # Check if using chunked processing
        if context.config.get("use_chunked_processing"):
            logger.info("Annotations will be applied during chunked processing")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for annotation")
            return context

        # Load annotation features
        # Create config dict for load_custom_features
        annotation_cfg = {
            "annotate_bed_files": context.annotation_configs.get("bed_files", []),
            "annotate_gene_lists": context.annotation_configs.get("gene_lists", []),
            "annotate_json_genes": context.annotation_configs.get("json_genes", []),
            "json_gene_mapping": context.annotation_configs.get("json_mapping", ""),
            "json_genes_as_columns": context.config.get("json_genes_as_columns", False),
        }
        features = load_custom_features(annotation_cfg)

        # Apply annotations
        logger.info(f"Applying {len(features)} custom annotation features")
        df = annotate_dataframe_with_features(df, features)

        context.current_dataframe = df
        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        input_files = []
        # Add BED files if available
        if hasattr(context, "annotation_configs") and context.annotation_configs:
            for bed_file in context.annotation_configs.get("bed_files", []):
                if Path(bed_file).exists():
                    input_files.append(Path(bed_file))
            for gene_list in context.annotation_configs.get("gene_lists", []):
                if Path(gene_list).exists():
                    input_files.append(Path(gene_list))
            for json_file in context.annotation_configs.get("json_genes", []):
                if Path(json_file).exists():
                    input_files.append(Path(json_file))
        return input_files

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        # CustomAnnotationStage modifies DataFrame in memory, no file outputs
        return []


class InheritanceAnalysisStage(Stage):
    """Calculate Mendelian inheritance patterns."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "inheritance_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Calculate inheritance patterns"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depend on dataframe_loading as a hard requirement
        # custom_annotation is optional and will run before if present
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after custom_annotation if it exists
        return {"custom_annotation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Calculate inheritance patterns."""
        # Check if we have necessary data
        df = context.current_dataframe
        if df is None and not context.config.get("use_chunked_processing"):
            logger.warning("No DataFrame loaded for inheritance analysis")
            return context

        # Determine inheritance mode
        inheritance_mode = context.config.get("inheritance_mode", "simple")

        # Check if using chunked processing - only defer if chunked processing is not complete
        if context.config.get("use_chunked_processing") and not context.config.get(
            "chunked_processing_complete"
        ):
            logger.info(
                f"Inheritance analysis ({inheritance_mode} mode) "
                "will be applied during chunked processing"
            )
            # Check if scoring configuration requires Inheritance_Details
            preserve_details = self._check_if_scoring_needs_details(context)
            if preserve_details:
                logger.debug("Preserving Inheritance_Details for chunked processing scoring")

            # Store inheritance config for chunked processing
            context.config["inheritance_analysis_config"] = {
                "inheritance_mode": inheritance_mode,
                "vcf_samples": context.vcf_samples,
                "pedigree_data": context.pedigree_data,
                "use_vectorized_comp_het": not context.config.get("no_vectorized_comp_het", False),
                "preserve_details_for_scoring": preserve_details,
            }
            return context

        # Get samples info
        vcf_samples = context.vcf_samples or []
        pedigree_data: dict[str, Any] | None = context.pedigree_data

        if not vcf_samples:
            logger.warning("No VCF samples available for inheritance analysis")
            return context

        # Ensure vcf_samples is a deterministically ordered list
        if isinstance(vcf_samples, set):
            # If somehow still a set, convert to sorted list for deterministic order
            vcf_samples = sorted(vcf_samples)
            logger.warning(
                "VCF samples was a set - converted to sorted list for deterministic order"
            )

        # At this point df must be a DataFrame (None case handled above)
        assert df is not None, "DataFrame must not be None for inheritance analysis"

        # Use ResourceManager to decide processing strategy
        from ..memory import ResourceManager

        rm = ResourceManager(config=context.config)

        # Calculate optimal chunk size and worker count
        chunk_size = rm.auto_chunk_size(len(df), len(vcf_samples))

        # Estimate memory per gene for worker calculation
        # Use average variants per gene as rough estimate (conservative)
        num_genes = df["GENE"].nunique() if "GENE" in df.columns else 1
        avg_variants_per_gene = len(df) / max(1, num_genes)
        memory_per_gene_gb = rm.estimate_memory(int(avg_variants_per_gene), len(vcf_samples))
        max_workers = rm.auto_workers(
            task_count=df["GENE"].nunique() if "GENE" in df.columns else 1,
            memory_per_task_gb=memory_per_gene_gb,
        )

        # Determine if dataset fits in memory or needs chunking
        should_chunk = len(df) > chunk_size

        if should_chunk:
            logger.info(
                f"Using chunked processing: chunk_size={chunk_size:,}, "
                f"max_workers={max_workers}, ~{memory_per_gene_gb:.1f}GB per gene"
            )
            # Force chunked processing and delegate to chunked analysis
            context.config["use_chunked_processing"] = True
            context.config["force_chunked_processing"] = True

            # Check if scoring configuration requires Inheritance_Details
            preserve_details = self._check_if_scoring_needs_details(context)
            if preserve_details:
                logger.info("Preserving Inheritance_Details for chunked processing scoring")

            context.config["inheritance_analysis_config"] = {
                "inheritance_mode": inheritance_mode,
                "vcf_samples": vcf_samples,
                "pedigree_data": pedigree_data,
                "use_vectorized_comp_het": not context.config.get("no_vectorized_comp_het", False),
                "preserve_details_for_scoring": preserve_details,
            }
            return context

        logger.info(
            f"Dataset fits in memory ({len(df):,} variants < {chunk_size:,} chunk size), "
            f"processing as single batch"
        )

        logger.info(
            f"Calculating inheritance patterns ({inheritance_mode} mode) "
            f"for {len(vcf_samples)} samples"
        )

        # Prepare sample columns for inheritance analysis.
        # Two paths: (1) packed GT column exists -> parse it into sample columns
        #            (2) Phase 11 per-sample GT columns (GEN_N__GT) -> rename to sample IDs
        sample_columns_exist = any(s in df.columns for s in vcf_samples)

        if not sample_columns_exist and "GT" in df.columns:
            logger.debug("Preparing sample columns from GT field for inheritance analysis")
            prep_start = self._start_subtask("sample_column_preparation")

            # Use the intelligent sample column creation function
            sample_column_method = context.config.get("sample_column_creation_method", "auto")
            df = create_sample_columns_from_gt_intelligent(
                df=df,
                vcf_samples=vcf_samples,
                separator=context.config.get("separator", ";"),
                snpsift_sep=context.config.get("extract_fields_separator", ","),
                method=sample_column_method,
            )

            self._end_subtask("sample_column_preparation", prep_start)

        elif not sample_columns_exist:
            # Phase 11: per-sample GT columns (GEN_N__GT) from bcftools query
            from ..stages.output_stages import _find_per_sample_gt_columns

            gt_cols = _find_per_sample_gt_columns(df)
            if gt_cols:
                n = min(len(gt_cols), len(vcf_samples))
                logger.info(
                    f"Creating sample columns from {n} per-sample GT columns "
                    f"for inheritance analysis"
                )
                prep_start = self._start_subtask("sample_column_preparation")
                for i in range(n):
                    df[vcf_samples[i]] = df[gt_cols[i]]
                self._end_subtask("sample_column_preparation", prep_start)

        # Apply inheritance analysis with error handling
        analysis_start = self._start_subtask("inheritance_calculation")

        try:
            # Check if we should use parallel processing
            threads = context.config.get("threads", 1)
            min_variants_for_parallel = context.config.get(
                "min_variants_for_parallel_inheritance", 100
            )

            if threads > 1 and len(df) >= min_variants_for_parallel:
                # Use parallel analyzer for better performance
                logger.info(f"Using parallel inheritance analyzer with {threads} workers")
                # Track detailed timing
                import time

                from ..inheritance.parallel_analyzer import analyze_inheritance_parallel

                comp_het_start = time.time()

                df = analyze_inheritance_parallel(
                    df=df,
                    sample_list=vcf_samples,
                    pedigree_data=pedigree_data or {},
                    use_vectorized_comp_het=not context.config.get("no_vectorized_comp_het", False),
                    n_workers=threads,
                    min_variants_for_parallel=min_variants_for_parallel,
                )

                # Record compound het timing (this is a major component)
                comp_het_time = time.time() - comp_het_start
                if comp_het_time > 1.0:  # Only record if significant
                    self._subtask_times["compound_het_analysis"] = comp_het_time
            else:
                # Use sequential analyzer for small datasets or single-threaded
                df = analyze_inheritance(
                    df=df,
                    sample_list=vcf_samples,
                    pedigree_data=pedigree_data or {},
                    use_vectorized_comp_het=not context.config.get("no_vectorized_comp_het", False),
                )

        except Exception as e:
            # Check if scoring configuration requires Inheritance_Details
            preserve_details = self._check_if_scoring_needs_details(context)
            df = handle_inheritance_analysis_error(
                df,
                e,
                preserve_details_for_scoring=preserve_details,
                context_description="non-chunked inheritance analysis",
            )

        self._end_subtask("inheritance_calculation", analysis_start)

        # Process output based on inheritance mode
        output_start = self._start_subtask("output_processing")
        from ..inheritance.analyzer import process_inheritance_output

        # Check if scoring configuration requires Inheritance_Details
        preserve_details = self._check_if_scoring_needs_details(context)
        if preserve_details:
            logger.debug("Preserving Inheritance_Details for scoring configuration")

        df = process_inheritance_output(
            df, inheritance_mode, preserve_details_for_scoring=preserve_details
        )
        self._end_subtask("output_processing", output_start)

        # Remove individual sample columns that were added for inheritance analysis
        # Use unified cleanup function for consistency with chunked processing
        if vcf_samples:
            cleanup_start = self._start_subtask("column_cleanup")
            # For non-chunked processing, remove ALL sample columns (don't preserve GT)
            df = cleanup_sample_columns(df, vcf_samples, preserve_columns=[])
            self._end_subtask("column_cleanup", cleanup_start)

        context.current_dataframe = df
        return context

    def _check_if_scoring_needs_details(self, context: PipelineContext) -> bool:
        """Check if scoring configuration requires Inheritance_Details column."""
        if not context.scoring_config:
            return False

        # Check if any scoring formulas reference details or segregation_p_value
        if "formulas" in context.scoring_config:
            for formula_dict in context.scoring_config["formulas"]:
                for formula_name, formula_expr in formula_dict.items():
                    if isinstance(formula_expr, str) and any(
                        var in formula_expr for var in ["details", "segregation_p_value"]
                    ):
                        logger.info(
                            f"Scoring formula '{formula_name}' requires Inheritance_Details"
                        )
                        return True

        # Check variable assignments for Inheritance_Details column dependency
        if "variables" in context.scoring_config:
            variables = context.scoring_config["variables"]
            if "Inheritance_Details" in variables:
                logger.info("Scoring configuration maps Inheritance_Details column")
                return True

        return False

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        input_files = []
        # Include PED file if available for inheritance analysis
        if hasattr(context, "ped_file") and context.ped_file:
            input_files.append(Path(context.ped_file))
        return input_files

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        # InheritanceAnalysisStage modifies DataFrame in memory, no file outputs
        return []


class VariantScoringStage(Stage):
    """Apply variant scoring formulas."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "variant_scoring"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Calculate variant scores"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depend on dataframe_loading as a hard requirement
        # Other stages (custom_annotation, inheritance_analysis) are optional
        # The pipeline runner will ensure correct ordering based on what stages are present
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after annotations and inheritance if they exist
        return {"custom_annotation", "inheritance_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply variant scoring."""
        if not context.scoring_config:
            logger.debug("No scoring configuration loaded")
            return context

        # Check if using chunked processing - only defer if chunked processing is not complete
        if context.config.get("use_chunked_processing") and not context.config.get(
            "chunked_processing_complete"
        ):
            logger.info("Scoring will be applied during chunked processing")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for scoring")
            return context

        logger.info("Applying variant scoring formulas")

        # Apply scoring
        df = apply_scoring(df, context.scoring_config)

        context.current_dataframe = df
        return context


class GenotypeFilterStage(Stage):
    """Filter variants by genotype patterns."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "genotype_filtering"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Filter variants by genotype patterns"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Must run after dataframe is loaded and scored
        return {"dataframe_loading", "variant_scoring"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply genotype filtering if requested."""
        # Check if genotype filtering is requested
        genotype_filter = context.config.get("genotype_filter")
        gene_genotype_file = context.config.get("gene_genotype_file")

        if not genotype_filter and not gene_genotype_file:
            logger.debug("No genotype filtering requested")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for genotype filtering")
            return context

        # Parse genotype modes
        genotype_modes = set()
        if genotype_filter:
            genotype_modes = {g.strip() for g in genotype_filter.split(",") if g.strip()}

        logger.info(f"Applying genotype filtering with modes: {genotype_modes}")
        if gene_genotype_file:
            logger.info(f"Using gene-specific genotype file: {gene_genotype_file}")

        # Save DataFrame to temporary file for filtering
        import tempfile

        # Use compressed temporary files for genotype filtering to save space
        use_compression = context.config.get("gzip_intermediates", True)
        suffix = ".tsv.gz" if use_compression else ".tsv"

        with tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False) as tmp_input:
            if use_compression:
                # For compressed files, write without pandas compression first
                import gzip

                tmp_input.close()  # Close the file handle
                with gzip.open(tmp_input.name, "wt", compresslevel=1) as gz_file:
                    df.to_csv(gz_file, sep="\t", index=False)
                tmp_input_path = tmp_input.name
            else:
                df.to_csv(tmp_input.name, sep="\t", index=False)
                tmp_input_path = tmp_input.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False) as tmp_output:
            tmp_output_path = tmp_output.name

        try:
            # Apply genotype filtering
            filter_final_tsv_by_genotype(
                input_tsv=tmp_input_path,
                output_tsv=tmp_output_path,
                global_genotypes=genotype_modes,
                gene_genotype_file=gene_genotype_file,
            )

            # Load filtered results back into DataFrame
            filtered_df = pd.read_csv(tmp_output_path, sep="\t", low_memory=False)

            logger.info(f"Genotype filtering: {len(df)} → {len(filtered_df)} variants")

            # Update context with filtered DataFrame
            context.current_dataframe = filtered_df

        finally:
            # Clean up temporary files
            import os

            try:
                os.unlink(tmp_input_path)
                os.unlink(tmp_output_path)
            except OSError:
                pass

        return context


class StatisticsGenerationStage(Stage):
    """Generate summary statistics."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "statistics_generation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate summary statistics"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depends on dataframe being loaded
        # Statistics can be run on whatever columns are available
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Should run after chunked analysis if it's present
        return {"chunked_analysis"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - read-only computation on DataFrame

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the statistics file path
        so that subsequent stages (like Excel generation) can find the statistics.
        """
        if context.config.get("no_stats"):
            return context

        # Try to find existing statistics file
        stats_output = context.config.get("stats_output_file")
        if not stats_output:
            # Check for default statistics file path
            default_stats_path = context.workspace.get_intermediate_path("statistics.tsv")
            if default_stats_path.exists():
                stats_output = str(default_stats_path)
                context.config["stats_output_file"] = stats_output
                logger.info(f"Restored statistics file path: {stats_output}")
            else:
                logger.warning("Statistics file not found during checkpoint skip")

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate statistics."""
        if context.config.get("no_stats"):
            logger.debug("Statistics generation disabled")
            return context

        # Check if using chunked processing - generate statistics after chunks are combined
        if context.config.get("use_chunked_processing") and not context.config.get(
            "chunked_processing_complete"
        ):
            logger.info("Statistics will be generated after chunked processing completes")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for statistics")
            return context

        # Get stats config
        stats_config_path = context.config.get("stats_config")
        config = stats_config_path if stats_config_path else self._get_default_stats_config()

        logger.info("Generating summary statistics")

        # Create stats engine and compute statistics
        engine = StatsEngine(config)
        stats = engine.compute(df)

        context.statistics = stats

        # Write stats if output file specified or create default path
        stats_output = context.config.get("stats_output_file")
        if not stats_output and not context.config.get("no_stats"):
            # Create default stats file path like the old pipeline does
            stats_output = context.workspace.get_intermediate_path("statistics.tsv")
            logger.debug(f"Created default statistics output path: {stats_output}")

        if stats_output:
            # Always save the stats output path to context for other stages to use
            context.config["stats_output_file"] = str(stats_output)
            self._write_statistics(stats, stats_output)
            logger.debug(f"Statistics written to: {stats_output}")

        return context

    def _write_statistics(self, stats: dict[str, Any], output_file: str) -> None:
        """Write statistics to file in a consistent format."""
        all_data = []

        # Collect dataset statistics
        if "dataset" in stats and not stats["dataset"].empty:
            dataset_df = stats["dataset"].copy()
            dataset_df["category"] = "dataset"
            all_data.append(dataset_df)

        # Collect gene statistics
        if "genes" in stats and not stats["genes"].empty:
            genes_df = stats["genes"].copy()
            genes_df["category"] = "gene"
            all_data.append(genes_df)

        # Collect any other grouped statistics
        for key, value in stats.items():
            if key not in ["dataset", "genes"] and hasattr(value, "empty") and not value.empty:
                other_df = value.copy()
                other_df["category"] = key
                all_data.append(other_df)

        # Combine all data into a single consistent DataFrame
        if all_data:
            combined_df = pd.concat(all_data, ignore_index=True, sort=False)
            combined_df.to_csv(output_file, sep="\t", index=False)
        else:
            # Write a minimal file if no data
            with open(output_file, "w") as f:
                f.write("statistic\tvalue\tcategory\n")
                f.write("total_variants\t0\tdataset\n")

    def _get_default_stats_config(self) -> dict[str, Any]:
        """Get default statistics configuration using basic columns."""
        return {
            "stats_version": "1.0",
            "description": "Default statistics configuration using basic columns",
            "dataset_stats": [
                {
                    "name": "total_variants",
                    "expression": "len(df)",
                    "description": "Total number of variants",
                },
                {
                    "name": "total_samples",
                    "expression": "len([col for col in df.columns if '_GT' in col])",
                    "description": "Total number of samples",
                },
                {
                    "name": "unique_genes",
                    "expression": "df['GENE'].nunique() if 'GENE' in df.columns else 0",
                    "description": "Number of unique genes",
                },
                {
                    "name": "high_impact_variants",
                    "expression": "(df['IMPACT'] == 'HIGH').sum() if 'IMPACT' in df.columns else 0",
                    "description": "Number of HIGH impact variants",
                },
                {
                    "name": "moderate_impact_variants",
                    "expression": (
                        "(df['IMPACT'] == 'MODERATE').sum() if 'IMPACT' in df.columns else 0"
                    ),
                    "description": "Number of MODERATE impact variants",
                },
            ],
            "gene_stats": [
                {
                    "name": "variant_count",
                    "expression": "len(group_df)",
                    "groupby": "GENE",
                    "description": "Number of variants per gene",
                },
                {
                    "name": "high_impact_count",
                    "expression": (
                        "(group_df['IMPACT'] == 'HIGH').sum() "
                        "if 'IMPACT' in group_df.columns else 0"
                    ),
                    "groupby": "GENE",
                    "description": "Number of HIGH impact variants per gene",
                },
                {
                    "name": "moderate_impact_count",
                    "expression": (
                        "(group_df['IMPACT'] == 'MODERATE').sum() "
                        "if 'IMPACT' in group_df.columns else 0"
                    ),
                    "groupby": "GENE",
                    "description": "Number of MODERATE impact variants per gene",
                },
            ],
        }


class VariantAnalysisStage(Stage):
    """Run variant-level analysis to add analysis columns."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "variant_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Perform variant-level analysis"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on having a DataFrame with variants
        deps = {"dataframe_loading"}
        return deps

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Note: variant_identifier should run AFTER this stage, not before,
        # because analyze_variants creates a new DataFrame
        # Run after custom_annotation to ensure deterministic column ordering
        # Run after genotype filtering if present
        # CRITICAL: Run after inheritance_analysis to preserve inheritance columns
        return {"custom_annotation", "genotype_filtering", "inheritance_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Run variant analysis."""
        if context.current_dataframe is None:
            logger.warning("No DataFrame available for variant analysis")
            return context

        df = context.current_dataframe

        # Check if inheritance analysis was requested
        should_have_inheritance = self._should_calculate_inheritance(context)

        # Check if inheritance columns exist in the input DataFrame
        inheritance_cols_present = [
            col for col in df.columns if col in ["Inheritance_Pattern", "Inheritance_Details"]
        ]
        if inheritance_cols_present:
            logger.debug(f"Input DataFrame has inheritance columns: {inheritance_cols_present}")
            for col in inheritance_cols_present:
                non_null_count = df[col].notna().sum()
                logger.debug(f"Column '{col}' has {non_null_count} non-null values")
        else:
            if should_have_inheritance:
                logger.warning(
                    "No inheritance columns found in input DataFrame - "
                    "they were lost before VariantAnalysisStage"
                )
            else:
                logger.debug(
                    "No inheritance columns found - inheritance analysis was not requested"
                )

        # Prepare config for analyze_variants
        analysis_config = {
            "vcf_sample_names": context.vcf_samples or [],
            "sample_list": ",".join(context.vcf_samples or []),
            "no_db_links": context.config.get("no_db_links", False),
            "base_name": context.workspace.base_name,
            "reference": context.config.get("reference", "GRCh37"),
            # CRITICAL: Pass through case/control samples from phenotype assignment
            "case_samples": context.config.get("case_samples") or [],
            "control_samples": context.config.get("control_samples") or [],
            # Pass through phenotype terms for legacy compatibility
            "case_phenotypes": context.config.get("case_phenotypes") or [],
            "control_phenotypes": context.config.get("control_phenotypes") or [],
        }

        # Phase 11: Reconstruct packed GT column for analyze_variants if needed.
        # analyze_variants expects a packed GT column ("Sample1(0/1);Sample2(1/1)").
        # With Phase 11, we have per-sample GEN_N__GT columns instead.
        if "GT" not in df.columns and context.vcf_samples:
            from ..stages.output_stages import _find_per_sample_gt_columns, reconstruct_gt_column

            gt_cols = _find_per_sample_gt_columns(df)
            if gt_cols:
                logger.info("Reconstructing GT column for variant analysis")
                # Work on a copy to avoid modifying the original DataFrame
                df = reconstruct_gt_column(df.copy(), context.vcf_samples)

        # Restore original column names before writing temp TSV.
        # analyze_variants expects unsanitized names (GENE, not ANN_0__GENE).
        if context.column_rename_map:
            reverse_map = {v: k for k, v in context.column_rename_map.items()}
            df = df.rename(columns=reverse_map)

        # Create temporary TSV file for analyze_variants
        temp_tsv = context.workspace.get_intermediate_path("temp_for_analysis.tsv")

        # Use compression for temporary analysis files to save space
        use_compression = context.config.get("gzip_intermediates", True)
        if use_compression:
            temp_tsv = Path(str(temp_tsv) + ".gz")
            compression = "gzip"
        else:
            compression = None

        df.to_csv(temp_tsv, sep="\t", index=False, compression=compression)

        logger.info("Running variant-level analysis")

        # Run analyze_variants and collect results
        analysis_results = []
        # Use appropriate file opener based on compression
        if use_compression:
            import gzip

            def open_func(f):
                return gzip.open(f, "rt")

        else:

            def open_func(f):
                return open(f)

        with open_func(temp_tsv) as inp:
            for line in analyze_variants(inp, analysis_config):
                analysis_results.append(line)

        # Parse results back into DataFrame
        logger.debug(f"analyze_variants produced {len(analysis_results)} result lines")
        if analysis_results:
            # First line is header
            header = analysis_results[0].split("\t")
            data = [line.split("\t") for line in analysis_results[1:]]

            # Create new DataFrame with analysis results
            analysis_df = pd.DataFrame(data, columns=header)
            logger.debug(
                f"Analysis DataFrame created with {len(analysis_df)} rows and "
                f"columns: {list(analysis_df.columns)}"
            )

            # The analyze_variants function may filter or reorder rows, so we need to match
            # them properly. We'll use key columns (CHROM, POS, REF, ALT) to merge back missing
            # columns
            key_cols = ["CHROM", "POS", "REF", "ALT"]

            # Check key columns for merging
            orig_key_cols_present = [col for col in key_cols if col in df.columns]
            analysis_key_cols_present = [col for col in key_cols if col in analysis_df.columns]
            logger.debug(f"Original DataFrame key columns: {orig_key_cols_present}")
            logger.debug(f"Analysis DataFrame key columns: {analysis_key_cols_present}")

            # Check if we have the key columns in both dataframes
            if all(col in df.columns for col in key_cols) and all(
                col in analysis_df.columns for col in key_cols
            ):
                # Get columns from original df that are not in analysis_df
                original_only_cols = [col for col in df.columns if col not in analysis_df.columns]

                if original_only_cols:
                    # Prioritize inheritance columns and other critical columns
                    inheritance_cols = [
                        col
                        for col in original_only_cols
                        if col in ["Inheritance_Pattern", "Inheritance_Details"]
                    ]
                    if inheritance_cols:
                        logger.debug(f"Preserving inheritance columns: {inheritance_cols}")

                    # Create a subset of original df with key columns and missing columns
                    merge_df = df[key_cols + original_only_cols].copy()

                    # Merge to get the missing columns
                    try:
                        analysis_df = pd.merge(analysis_df, merge_df, on=key_cols, how="left")
                        logger.debug(
                            f"Successfully preserved columns from original dataframe: "
                            f"{original_only_cols}"
                        )

                        # Verify inheritance columns were preserved
                        if inheritance_cols:
                            preserved_inheritance = [
                                col for col in inheritance_cols if col in analysis_df.columns
                            ]
                            if preserved_inheritance:
                                logger.debug(
                                    f"Inheritance columns successfully preserved: "
                                    f"{preserved_inheritance}"
                                )
                                # Check if inheritance columns have actual data
                                for col in preserved_inheritance:
                                    non_null_count = analysis_df[col].notna().sum()
                                    logger.debug(
                                        f"Column '{col}' has {non_null_count} non-null values"
                                    )
                            else:
                                if should_have_inheritance:
                                    logger.error(
                                        f"Failed to preserve inheritance columns: "
                                        f"{inheritance_cols}"
                                    )
                                else:
                                    logger.debug(
                                        f"No inheritance columns to preserve: {inheritance_cols}"
                                    )
                    except Exception as e:
                        logger.error(f"Error merging columns: {e}")
                        logger.warning("Inheritance columns may be lost due to merge failure")

                    # Reorder columns to put VAR_ID first if it exists
                    if "VAR_ID" in analysis_df.columns:
                        cols = analysis_df.columns.tolist()
                        cols.remove("VAR_ID")
                        cols = ["VAR_ID", *cols]
                        analysis_df = analysis_df[cols]

                    # Put Custom_Annotation after GT if both exist
                    if "Custom_Annotation" in analysis_df.columns and "GT" in analysis_df.columns:
                        cols = analysis_df.columns.tolist()
                        cols.remove("Custom_Annotation")
                        gt_idx = cols.index("GT")
                        cols.insert(gt_idx + 1, "Custom_Annotation")
                        analysis_df = analysis_df[cols]

            else:
                logger.warning(
                    "Could not preserve columns from original dataframe - missing key columns "
                    "for merging"
                )
                # Check if inheritance columns exist in original but not in analysis result
                inheritance_cols = [
                    col
                    for col in df.columns
                    if col in ["Inheritance_Pattern", "Inheritance_Details"]
                ]
                if inheritance_cols and should_have_inheritance:
                    logger.error(
                        f"CRITICAL: Inheritance columns {inheritance_cols} will be lost "
                        f"due to missing key columns for merge"
                    )

            context.current_dataframe = analysis_df
            logger.info(f"Variant analysis complete: {len(analysis_df)} variants analyzed")

            # Final check: Verify if inheritance columns survived
            final_inheritance_cols = [
                col
                for col in analysis_df.columns
                if col in ["Inheritance_Pattern", "Inheritance_Details"]
            ]
            if final_inheritance_cols:
                logger.debug(f"Final DataFrame has inheritance columns: {final_inheritance_cols}")
                for col in final_inheritance_cols:
                    non_null_count = analysis_df[col].notna().sum()
                    logger.debug(f"Final column '{col}' has {non_null_count} non-null values")
            else:
                if should_have_inheritance:
                    logger.error(
                        "INHERITANCE COLUMNS LOST - Final DataFrame missing inheritance columns!"
                    )
                else:
                    logger.debug("Final DataFrame has no inheritance columns - none were expected")
        else:
            logger.warning("No analysis results generated - DataFrame unchanged")
            # Check if inheritance columns exist in the unchanged DataFrame
            inheritance_cols_present = [
                col for col in df.columns if col in ["Inheritance_Pattern", "Inheritance_Details"]
            ]
            if inheritance_cols_present:
                logger.debug(
                    f"Unchanged DataFrame still has inheritance columns: {inheritance_cols_present}"
                )
            else:
                if should_have_inheritance:
                    logger.warning("Unchanged DataFrame also missing inheritance columns")
                else:
                    logger.debug(
                        "Unchanged DataFrame has no inheritance columns - none were expected"
                    )

        # Clean up temp file
        if temp_tsv.exists():
            temp_tsv.unlink()

        return context

    def _should_calculate_inheritance(self, context: PipelineContext) -> bool:
        """Check if inheritance analysis should be calculated."""
        # Use the same logic as pipeline.py
        has_ped_file = context.pedigree_data is not None and len(context.pedigree_data) > 0
        has_inheritance_mode = context.config.get("inheritance_mode")
        has_calculate_inheritance_config = context.config.get("calculate_inheritance", False)
        requires_inheritance_for_scoring = self._check_if_scoring_needs_details(context)

        should_calculate = (
            has_ped_file
            or has_inheritance_mode
            or has_calculate_inheritance_config
            or requires_inheritance_for_scoring
        )

        logger.debug(
            f"Should calculate inheritance: {should_calculate} (ped={has_ped_file}, "
            f"mode={has_inheritance_mode}, config={has_calculate_inheritance_config}, "
            f"scoring={requires_inheritance_for_scoring})"
        )
        return bool(should_calculate)

    def _check_if_scoring_needs_details(self, context: PipelineContext) -> bool:
        """Check if scoring configuration requires Inheritance_Details."""
        if not context.scoring_config:
            return False

        # Check if any scoring formula references Inheritance_Details
        scoring_formulas = context.scoring_config.get("scoring_formulas", {})
        for formula_name, formula_data in scoring_formulas.items():
            formula_str = formula_data.get("formula", "")
            if "Inheritance_Details" in formula_str:
                logger.debug(f"Scoring formula '{formula_name}' requires Inheritance_Details")
                return True

        return False

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        # VariantAnalysisStage works with DataFrame in memory, no specific input files
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        # VariantAnalysisStage modifies DataFrame in memory, no file outputs
        return []


class GeneBurdenAnalysisStage(Stage):
    """Perform gene burden analysis for case-control studies."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "gene_burden_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Perform gene burden analysis"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on having a DataFrame with variants and case/control samples
        return {"dataframe_loading", "sample_config_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after custom_annotation if it exists
        return {"custom_annotation"}

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the gene burden output file path
        so that subsequent stages (like Excel generation) can find the results.
        """
        if not context.config.get("perform_gene_burden"):
            return context

        # Try to find existing gene burden output file
        burden_output = context.config.get("gene_burden_output")
        if not burden_output:
            # Check for default gene burden output file path
            output_dir = context.config.get("output_dir", "output")
            base_name = context.config.get("output_file_base", "gene_burden_results")

            # Try both compressed and uncompressed versions
            default_path = Path(output_dir) / f"{base_name}.gene_burden.tsv"
            default_path_gz = Path(str(default_path) + ".gz")

            if default_path_gz.exists():
                burden_output = str(default_path_gz)
                context.config["gene_burden_output"] = burden_output
                logger.info(f"Restored gene burden file path: {burden_output}")
            elif default_path.exists():
                burden_output = str(default_path)
                context.config["gene_burden_output"] = burden_output
                logger.info(f"Restored gene burden file path: {burden_output}")
            else:
                logger.warning("Gene burden output file not found during checkpoint skip")

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Perform gene burden analysis."""
        if not context.config.get("perform_gene_burden"):
            logger.debug("Gene burden analysis not requested")
            return context

        # Check if using chunked processing - gene burden should be done after chunks are combined
        if context.config.get("use_chunked_processing") and not context.config.get(
            "chunked_processing_complete"
        ):
            logger.info("Gene burden analysis will be performed after chunked processing completes")
            return context

        # Check required inputs
        case_samples = context.config.get("case_samples", [])
        control_samples = context.config.get("control_samples", [])

        # Debug logging to understand what's available
        logger.debug(
            f"Gene burden analysis context check: "
            f"case_samples={len(case_samples) if case_samples else 0}, "
            f"control_samples={len(control_samples) if control_samples else 0}"
        )

        # Debug all context config keys to understand what's available
        config_keys = list(context.config.keys())
        sample_related_keys = [k for k in config_keys if "sample" in k.lower()]
        logger.debug(f"All sample-related config keys: {sample_related_keys}")

        if not case_samples or not control_samples:
            logger.warning(
                f"Case/control samples not defined for gene burden analysis: "
                f"case_samples={len(case_samples) if case_samples else 0}, "
                f"control_samples={len(control_samples) if control_samples else 0}"
            )
            logger.debug(f"Available config keys with 'sample': {sample_related_keys}")
            return context

        df = context.current_dataframe
        if df is None or df.empty:
            logger.warning("No DataFrame loaded for gene burden analysis")
            return context

        # Check if DataFrame has required GENE column
        if "GENE" not in df.columns:
            logger.error("DataFrame missing required 'GENE' column for gene burden analysis")
            logger.debug(f"Available columns: {list(df.columns)[:20]}...")
            return context

        # Phase 11: Reconstruct packed GT column if missing
        if "GT" not in df.columns and context.vcf_samples:
            from ..stages.output_stages import _find_per_sample_gt_columns, reconstruct_gt_column

            gt_cols = _find_per_sample_gt_columns(df)
            if gt_cols:
                logger.info("Reconstructing GT column for gene burden analysis")
                df = reconstruct_gt_column(df.copy(), context.vcf_samples)
                context.current_dataframe = df

        if "GT" not in df.columns:
            logger.error("DataFrame missing required 'GT' column for gene burden analysis")
            return context

        logger.debug(
            f"Performing gene burden analysis: {len(case_samples)} cases, "
            f"{len(control_samples)} controls"
        )

        # Run analysis with proper config
        burden_config = {
            "gene_burden_mode": context.config.get("gene_burden_mode", "alleles"),
            "correction_method": context.config.get("correction_method", "fdr"),
            "confidence_interval_method": context.config.get(
                "confidence_interval_method", "normal_approx"
            ),
            "confidence_interval_alpha": context.config.get("confidence_interval_alpha", 0.05),
        }

        # Add case/control count columns to DataFrame first
        from ..helpers import assign_case_control_counts

        # Use VCF samples as all_samples - this is critical!
        # The assign_case_control_counts function needs to iterate over ALL samples
        # that appear in the VCF, not just those in case+control lists
        all_vcf_samples = set(context.vcf_samples) if context.vcf_samples else set()

        logger.info(f"Using {len(all_vcf_samples)} VCF samples for gene burden analysis")

        import time as _time

        _t_cc = _time.monotonic()
        df_with_counts = assign_case_control_counts(
            df=df,
            case_samples=set(case_samples),
            control_samples=set(control_samples),
            all_samples=all_vcf_samples,  # This should be ALL samples in VCF
        )
        _t_cc_elapsed = _time.monotonic() - _t_cc
        logger.info(
            f"Case/control count assignment completed in {_t_cc_elapsed:.2f}s "
            f"({len(df_with_counts)} variants, {len(all_vcf_samples)} samples)"
        )

        # Perform gene burden analysis with proper per-sample collapsing
        # Pass vcf_samples to enable fast column-based aggregation (Phase 11)
        _t_burden = _time.monotonic()
        burden_results = perform_gene_burden_analysis(
            df=df_with_counts,
            cfg=burden_config,
            case_samples=set(case_samples),
            control_samples=set(control_samples),
            vcf_samples=list(context.vcf_samples) if context.vcf_samples else None,
        )
        _t_burden_elapsed = _time.monotonic() - _t_burden
        n_genes = len(burden_results) if burden_results is not None else 0
        logger.info(
            f"Gene burden analysis completed in {_t_burden_elapsed:.2f}s ({n_genes} genes tested)"
        )

        context.gene_burden_results = burden_results

        # Write results to output file
        burden_output = context.config.get("gene_burden_output")
        if not burden_output:
            # Create default gene burden output path
            output_dir = context.config.get("output_dir", "output")
            base_name = context.config.get("output_file_base", "gene_burden_results")
            burden_output = str(Path(output_dir) / f"{base_name}.gene_burden.tsv")

            # Save the generated path back into the context for other stages to use.
            context.config["gene_burden_output"] = burden_output

        # Apply compression to gene burden results based on configuration
        use_compression = context.config.get("gzip_intermediates", True)
        if use_compression and not str(burden_output).endswith(".gz"):
            burden_output = str(burden_output) + ".gz"
            context.config["gene_burden_output"] = burden_output
            compression = "gzip"
        else:
            compression = None

        burden_results.to_csv(burden_output, sep="\t", index=False, compression=compression)
        logger.info(f"Wrote gene burden results to {burden_output}")

        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        # GeneBurdenAnalysisStage works with DataFrame in memory, no specific input files
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        # Return gene burden output file if configured
        if hasattr(context, "config") and context.config.get("gene_burden_output"):
            burden_output = context.config["gene_burden_output"]
            if Path(burden_output).exists():
                return [Path(burden_output)]
        return []


class AssociationAnalysisStage(Stage):
    """Perform association analysis using the modular AssociationEngine framework."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "association_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Perform association analysis"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading", "sample_config_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {"custom_annotation"}

    @property
    def parallel_safe(self) -> bool:
        """Return False — rpy2/R SKAT backend is not thread-safe (SKAT-08).

        rpy2 is not thread-safe. Calling rpy2 functions from a ThreadPoolExecutor
        worker thread causes segfaults with no Python traceback. This explicit
        property documents the SKAT-08 requirement and prevents parallel execution
        regardless of which tests are active.
        """
        return False

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the association output file path
        so that subsequent stages can find the results.
        """
        if not context.config.get("perform_association"):
            return context

        # Try to find existing association output file
        assoc_output = context.config.get("association_output")
        if not assoc_output:
            output_dir = context.config.get("output_dir", "output")
            base_name = context.config.get("output_file_base", "association_results")

            default_path = Path(output_dir) / f"{base_name}.association.tsv"
            default_path_gz = Path(str(default_path) + ".gz")

            if default_path_gz.exists():
                assoc_output = str(default_path_gz)
                context.config["association_output"] = assoc_output
                logger.info(f"Restored association output file path: {assoc_output}")
            elif default_path.exists():
                assoc_output = str(default_path)
                context.config["association_output"] = assoc_output
                logger.info(f"Restored association output file path: {assoc_output}")
            else:
                logger.warning("Association output file not found during checkpoint skip")

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Perform association analysis."""
        if not context.config.get("perform_association"):
            logger.debug("Association analysis not requested")
            return context

        # Check required inputs
        case_samples = context.config.get("case_samples", [])
        control_samples = context.config.get("control_samples", [])

        if not case_samples or not control_samples:
            logger.warning(
                f"Case/control samples not defined for association analysis: "
                f"case_samples={len(case_samples) if case_samples else 0}, "
                f"control_samples={len(control_samples) if control_samples else 0}"
            )
            return context

        # Build AssociationConfig from context (Phase 19: includes covariate/trait/weight fields)
        # Phase 22: includes diagnostics_output and sample size warning thresholds
        assoc_config = AssociationConfig(
            correction_method=context.config.get("correction_method", "fdr"),
            gene_burden_mode=context.config.get("gene_burden_mode", "samples"),
            confidence_interval_method=context.config.get(
                "confidence_interval_method", "normal_approx"
            ),
            confidence_interval_alpha=context.config.get("confidence_interval_alpha", 0.05),
            covariate_file=context.config.get("covariate_file"),
            covariate_columns=context.config.get("covariate_columns"),
            categorical_covariates=context.config.get("categorical_covariates"),
            trait_type=context.config.get("trait_type", "binary"),
            variant_weights=context.config.get("variant_weights", "beta:1,25"),
            variant_weight_params=context.config.get("variant_weight_params"),
            diagnostics_output=context.config.get("diagnostics_output"),
            min_cases=context.config.get("association_min_cases", 200),
            max_case_control_ratio=context.config.get("association_max_case_control_ratio", 20.0),
            min_case_carriers=context.config.get("association_min_case_carriers", 10),
            # Phase 23: PCA fields
            pca_file=context.config.get("pca_file"),
            pca_tool=context.config.get("pca_tool"),
            pca_components=context.config.get("pca_components", 10),
        )
        logger.info(f"Association analysis: trait type = {assoc_config.trait_type}")

        # Get test names; default to fisher
        test_names: list[str] = context.config.get("association_tests", ["fisher"])

        # Eager dependency check: instantiate engine (also calls check_dependencies per test)
        engine = AssociationEngine.from_names(test_names, assoc_config)

        # Get DataFrame
        df = context.current_dataframe
        if df is None:
            df = context.variants_df
        if df is None or df.empty:
            logger.warning("No DataFrame loaded for association analysis")
            return context

        if "GENE" not in df.columns:
            logger.error("DataFrame missing required 'GENE' column for association analysis")
            return context

        # Phase 11: Reconstruct packed GT column if missing.
        # Phase 19: When regression tests need per-sample GT columns, keep a
        # reference to the pre-reconstruction DataFrame for genotype matrix
        # building. reconstruct_gt_column drops per-sample columns, so we
        # must save them first.
        needs_regression = any(
            t in test_names for t in ("logistic_burden", "linear_burden", "skat", "skat_python")
        )
        df_with_per_sample_gt: pd.DataFrame | None = None
        if "GT" not in df.columns and context.vcf_samples:
            from ..stages.output_stages import _find_per_sample_gt_columns, reconstruct_gt_column

            gt_cols = _find_per_sample_gt_columns(df)
            if gt_cols:
                # Save DataFrame with per-sample GT columns for genotype matrix
                if needs_regression:
                    df_with_per_sample_gt = df
                logger.info("Reconstructing GT column for association analysis")
                df = reconstruct_gt_column(df.copy(), context.vcf_samples)
                context.current_dataframe = df
        elif "GT" in df.columns and needs_regression and context.vcf_samples:
            # GT already reconstructed (e.g. by gene_burden_analysis at same level).
            # Try variants_df as fallback source for per-sample GT columns.
            from ..stages.output_stages import _find_per_sample_gt_columns

            fallback_df = context.variants_df
            if fallback_df is not None:
                gt_cols_fb = _find_per_sample_gt_columns(fallback_df)
                if gt_cols_fb:
                    logger.info(
                        f"Association analysis: recovered {len(gt_cols_fb)} "
                        "per-sample GT columns from variants_df for genotype matrix"
                    )
                    df_with_per_sample_gt = fallback_df

        # Determine aggregation strategy (same priority as perform_gene_burden_analysis)
        case_set = set(case_samples)
        control_set = set(control_samples)
        vcf_samples_list = list(context.vcf_samples) if context.vcf_samples else None

        # ------------------------------------------------------------------
        # Phase 19: Tiered sample size warnings (CONTEXT.md)
        # ------------------------------------------------------------------
        n_cases_total = len(case_samples)
        n_controls_total = len(control_samples)
        if n_cases_total < 10:
            logger.error(
                f"Association analysis: only {n_cases_total} case(s) found — "
                "fewer than 10 cases produces invalid results. Aborting."
            )
            return context
        if n_cases_total < 50:
            logger.warning(
                f"Association analysis: only {n_cases_total} case(s) — "
                "fewer than 50 cases provides no practical power for regression tests"
            )
        elif n_cases_total < 200:
            logger.warning(
                f"Association analysis: {n_cases_total} case(s) — "
                "fewer than 200 cases is underpowered for SKAT; interpret results with caution"
            )
        if n_controls_total > 0 and n_cases_total > 0:
            ratio = n_controls_total / n_cases_total
            if ratio > 20:
                logger.warning(
                    f"Association analysis: case:control ratio is 1:{ratio:.0f} — "
                    "exceeds 1:20; Type I error inflation risk without SPA/Firth correction"
                )

        # ------------------------------------------------------------------
        # Phase 19: Load covariates once (if covariate_file provided)
        # ------------------------------------------------------------------
        covariate_matrix = None
        covariate_col_names: list[str] = []
        if assoc_config.covariate_file and vcf_samples_list:
            from ..association.covariates import load_covariates

            covariate_matrix, covariate_col_names = load_covariates(
                assoc_config.covariate_file,
                vcf_samples_list,
                assoc_config.covariate_columns,
                assoc_config.categorical_covariates,
            )
            logger.info(
                f"Association analysis: loaded {covariate_matrix.shape[1]} "
                f"covariate(s): {covariate_col_names}"
            )

        # ------------------------------------------------------------------
        # Phase 23: Load PCA file and merge with covariates inline
        # ------------------------------------------------------------------
        pca_file = assoc_config.pca_file
        if pca_file and vcf_samples_list:
            from ..association.pca import load_pca_file, merge_pca_covariates

            n_pcs = assoc_config.pca_components
            pca_matrix, pca_col_names = load_pca_file(
                pca_file, vcf_samples_list, n_components=n_pcs
            )
            covariate_matrix, covariate_col_names = merge_pca_covariates(
                pca_matrix,
                pca_col_names,
                covariate_matrix,
                covariate_col_names,
            )
            logger.info(
                f"PCA: merged {len(pca_col_names)} PCs with existing covariates "
                f"-> {len(covariate_col_names)} total columns"
            )

        # ------------------------------------------------------------------
        # Phase 19: Build phenotype vector once (0=control, 1=case)
        # ------------------------------------------------------------------
        import numpy as np

        phenotype_vector = None
        if vcf_samples_list:
            case_set_all = set(case_samples)
            phenotype_vector = np.array(
                [1.0 if s in case_set_all else 0.0 for s in vcf_samples_list],
                dtype=float,
            )
            n_pv_cases = int(phenotype_vector.sum())
            n_pv_controls = int((1.0 - phenotype_vector).sum())
            logger.info(
                f"Association analysis: phenotype vector built — "
                f"{n_pv_cases} cases, {n_pv_controls} controls"
            )

        # ------------------------------------------------------------------
        # Standard aggregation (existing paths — unchanged)
        # ------------------------------------------------------------------
        has_case_ctrl = True  # already checked above
        gt_columns = _find_gt_columns(df)
        use_column_aggregation = bool(
            has_case_ctrl
            and gt_columns
            and vcf_samples_list
            and len(gt_columns) <= len(vcf_samples_list)
        )
        use_gt_aggregation = has_case_ctrl and not use_column_aggregation and "GT" in df.columns

        if use_column_aggregation:
            logger.info(
                f"Association analysis: using column-based aggregation "
                f"({len(gt_columns)} GT columns, "
                f"{len(case_set)} cases, {len(control_set)} controls)"
            )
            gene_burden_data = _aggregate_gene_burden_from_columns(
                df, case_set, control_set, vcf_samples_list, gt_columns
            )
        elif use_gt_aggregation:
            logger.info(
                "Association analysis: using packed GT string collapsing "
                f"({len(case_set)} cases, {len(control_set)} controls)"
            )
            gene_burden_data = _aggregate_gene_burden_from_gt(df, case_set, control_set)
        else:
            logger.info("Association analysis: using pre-computed per-variant counts (legacy mode)")
            gene_burden_data = _aggregate_gene_burden_legacy(df)

        if not gene_burden_data:
            logger.warning("No genes found with variant data for association analysis.")
            return context

        # ------------------------------------------------------------------
        # Phase 19: Augment gene_burden_data with genotype matrix for
        # regression tests (logistic_burden, linear_burden, skat, skat_python)
        # Backward compatible: FisherExactTest ignores the new keys.
        #
        # Use df_with_per_sample_gt when available (per-sample GT columns
        # were dropped by reconstruct_gt_column for the aggregation step).
        # Fall back to gt_columns from current df if columns still present.
        # ------------------------------------------------------------------
        gt_source_df = df_with_per_sample_gt if df_with_per_sample_gt is not None else df
        gt_columns_for_matrix = _find_gt_columns(gt_source_df)
        if needs_regression and gt_columns_for_matrix and vcf_samples_list:
            from ..association.genotype_matrix import build_genotype_matrix

            is_binary = assoc_config.trait_type == "binary"
            for gene_data in gene_burden_data:
                gene_name = gene_data.get("GENE", "")
                gene_df = gt_source_df[gt_source_df["GENE"] == gene_name]
                if gene_df.empty:
                    gene_data["genotype_matrix"] = np.zeros((len(vcf_samples_list), 0), dtype=float)
                    gene_data["variant_mafs"] = np.zeros(0, dtype=float)
                    gene_data["phenotype_vector"] = phenotype_vector
                    gene_data["covariate_matrix"] = covariate_matrix
                    continue

                geno, mafs, sample_mask, gt_warnings = build_genotype_matrix(
                    gene_df,
                    vcf_samples_list,
                    gt_columns_for_matrix,
                    is_binary=is_binary,
                    missing_site_threshold=assoc_config.missing_site_threshold,
                    missing_sample_threshold=assoc_config.missing_sample_threshold,
                    phenotype_vector=phenotype_vector,
                )
                for w in gt_warnings:
                    logger.warning(f"Gene {gene_name}: {w}")

                # Apply sample mask to phenotype and covariates if any high-missing samples
                if not all(sample_mask):
                    mask_arr = np.array(sample_mask, dtype=bool)
                    pv = phenotype_vector[mask_arr] if phenotype_vector is not None else None
                    cm = covariate_matrix[mask_arr] if covariate_matrix is not None else None
                    geno = geno[mask_arr]
                else:
                    pv = phenotype_vector
                    cm = covariate_matrix

                # Per-gene MAC check: skip regression if < 5 minor allele copies
                total_mac = int(geno.sum()) if geno.size > 0 else 0
                if total_mac < 5:
                    logger.debug(
                        f"Gene {gene_name}: MAC={total_mac} < 5 — "
                        "regression will report NA (insufficient data)"
                    )
                    gene_data["genotype_matrix"] = np.zeros((geno.shape[0], 0), dtype=float)
                    gene_data["variant_mafs"] = np.zeros(0, dtype=float)
                else:
                    gene_data["genotype_matrix"] = geno
                    gene_data["variant_mafs"] = mafs

                # Phase 23: Extract functional annotation columns for CADD/REVEL weight schemes
                # (WEIGHT-05). Arrays must align with variant_mafs (post site-filter length).
                if assoc_config.variant_weights in ("cadd", "revel", "combined"):
                    _cadd_col = next(
                        (
                            c
                            for c in gene_df.columns
                            if c.lower() in ("dbnsfp_cadd_phred", "cadd_phred")
                        ),
                        None,
                    )
                    _revel_col = next(
                        (
                            c
                            for c in gene_df.columns
                            if c.lower() in ("dbnsfp_revel_score", "revel_score")
                        ),
                        None,
                    )
                    _effect_col = next(
                        (c for c in gene_df.columns if c.upper() in ("EFFECT", "ANN_0__EFFECT")),
                        None,
                    )
                    # Align annotations with variant_mafs (build_genotype_matrix applies
                    # a site missing-rate filter; replicate it here for correct alignment).
                    # The filter marks variants with >missing_site_threshold missing GTs.
                    # Missing GTs are identified by parse_gt_to_dosage returning None,
                    # i.e. strings like "./.", ".|.", ".", or empty.
                    _n_df = len(gene_df)
                    _n_kept = len(gene_data["variant_mafs"])
                    if _n_kept < _n_df:
                        from ..association.genotype_matrix import parse_gt_to_dosage as _pgd

                        _gt_cols_list = list(gt_columns_for_matrix)
                        _n_samples_gt = len(_gt_cols_list)
                        _keep_mask_ann = np.ones(_n_df, dtype=bool)
                        for _vi, (_, _row) in enumerate(gene_df.iterrows()):
                            _n_miss = sum(
                                1
                                for _col in _gt_cols_list
                                if _pgd(str(_row.get(_col) or ""))[0] is None
                            )
                            _miss_frac = _n_miss / _n_samples_gt if _n_samples_gt > 0 else 0.0
                            if _miss_frac > assoc_config.missing_site_threshold:
                                _keep_mask_ann[_vi] = False
                    else:
                        _keep_mask_ann = None

                    if _cadd_col:
                        _vals = gene_df[_cadd_col].values
                        gene_data["cadd_scores"] = (
                            _vals[_keep_mask_ann] if _keep_mask_ann is not None else _vals
                        )
                    if _revel_col:
                        _vals = gene_df[_revel_col].values
                        gene_data["revel_scores"] = (
                            _vals[_keep_mask_ann] if _keep_mask_ann is not None else _vals
                        )
                    if _effect_col:
                        _vals = gene_df[_effect_col].values
                        gene_data["variant_effects"] = (
                            _vals[_keep_mask_ann] if _keep_mask_ann is not None else _vals
                        )

                gene_data["phenotype_vector"] = pv
                gene_data["covariate_matrix"] = cm
                gene_data["vcf_samples"] = vcf_samples_list

            # Release the per-sample GT DataFrame (can be large with many samples)
            del df_with_per_sample_gt, gt_source_df

        # Phase 22: Build gene lookup for per-gene warnings BEFORE engine.run_all()
        # gene_burden_data dicts use "GENE" (uppercase) as the gene key
        gene_data_by_gene = {d.get("GENE", d.get("gene", "")): d for d in gene_burden_data}

        # Run association tests
        results_df = engine.run_all(gene_burden_data)

        if results_df.empty:
            logger.warning("Association analysis produced no results.")
            return context

        # ------------------------------------------------------------------
        # Phase 22: Per-gene warnings column
        # ------------------------------------------------------------------
        from ..association.diagnostics import compute_per_gene_warnings

        warnings_by_gene: dict[str, str] = {}
        for gene in results_df["gene"].unique():
            gdata = gene_data_by_gene.get(gene, {})
            case_carriers = gdata.get("proband_carrier_count", 0)
            gene_warnings = compute_per_gene_warnings(gene, case_carriers, assoc_config)
            if gene_warnings:
                warnings_by_gene[gene] = ";".join(gene_warnings)
        results_df["warnings"] = results_df["gene"].map(warnings_by_gene).fillna("")

        # Write results to output file
        assoc_output = context.config.get("association_output")
        if not assoc_output:
            output_dir = context.config.get("output_dir", "output")
            base_name = context.config.get("output_file_base", "association_results")
            assoc_output = str(Path(output_dir) / f"{base_name}.association.tsv")
            context.config["association_output"] = assoc_output

        # Apply compression based on configuration
        use_compression = context.config.get("gzip_intermediates", True)
        if use_compression and not str(assoc_output).endswith(".gz"):
            assoc_output = str(assoc_output) + ".gz"
            context.config["association_output"] = assoc_output
            compression = "gzip"
        else:
            compression = None

        results_df.to_csv(assoc_output, sep="\t", index=False, compression=compression)
        logger.info(f"Wrote association results to {assoc_output}")

        # ------------------------------------------------------------------
        # Phase 22: Write diagnostics if requested
        # ------------------------------------------------------------------
        if assoc_config.diagnostics_output:
            from ..association.diagnostics import emit_sample_size_warnings, write_diagnostics

            cohort_warnings = emit_sample_size_warnings(
                n_cases_total, n_controls_total, assoc_config
            )
            write_diagnostics(
                results_df=results_df,
                diagnostics_dir=assoc_config.diagnostics_output,
                test_names=test_names,
                n_cases=n_cases_total,
                n_controls=n_controls_total,
                cohort_warnings=cohort_warnings,
            )

        # Store results in context
        context.association_results = results_df

        # Log summary
        n_genes = len(results_df)
        # Count significant genes across any corrected p-value column
        corr_cols = [c for c in results_df.columns if c.endswith("_corrected_p_value")]
        n_sig = int((results_df[corr_cols].min(axis=1) < 0.05).sum()) if corr_cols else 0
        logger.info(
            f"Association analysis: {n_genes} genes tested, {n_sig} significant (FDR < 0.05)"
        )

        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        if hasattr(context, "config") and context.config.get("association_output"):
            assoc_output = context.config["association_output"]
            if Path(assoc_output).exists():
                return [Path(assoc_output)]
        return []


class ClinVarPM5Stage(Stage):
    """Annotate variants with ACMG PM5 criterion from ClinVar lookup."""

    @property
    def name(self) -> str:
        return "clinvar_pm5"

    @property
    def description(self) -> str:
        return "Annotate variants with ClinVar PM5 criterion"

    @property
    def dependencies(self) -> set[str]:
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        return {"variant_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        lookup_path = context.config.get("clinvar_pm5_lookup")
        if not lookup_path:
            logger.debug("No --clinvar-pm5-lookup provided, skipping PM5 annotation")
            return context

        if context.current_dataframe is None:
            logger.warning("No DataFrame available for PM5 annotation")
            return context

        from ..clinvar_pm5 import annotate_pm5, load_pm5_lookup

        lookup = load_pm5_lookup(lookup_path)

        # Resolve gene and HGVS_P column names
        df = context.current_dataframe
        gene_col = "GENE" if "GENE" in df.columns else "ANN[0].GENE"
        hgvs_col = "HGVS_P" if "HGVS_P" in df.columns else "ANN[0].HGVS_P"

        if gene_col not in df.columns:
            logger.warning("Gene column not found in DataFrame, skipping PM5")
            return context

        if hgvs_col not in df.columns:
            logger.warning("HGVS_P column not found in DataFrame, skipping PM5")
            return context

        logger.info("Running PM5 annotation on %d variants", len(df))
        context.current_dataframe = annotate_pm5(df, lookup, gene_col=gene_col, hgvs_col=hgvs_col)
        logger.info("PM5 annotation complete")
        return context


class ChunkedAnalysisStage(Stage):
    """Process large files in chunks with all analysis steps."""

    def __init__(self, max_variants_per_gene: int = 10000):
        """Initialize chunked analysis stage with configurable gene size limits.

        Args
        ----
            max_variants_per_gene: Maximum number of variants allowed per gene
                                 before forcing chunk split (default: 10,000)
        """
        super().__init__()
        self.max_variants_per_gene = max_variants_per_gene

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "chunked_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Process large files in memory-efficient chunks"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on configuration and processing stages but not full DataFrame loading
        # This stage only runs when use_chunked_processing is True
        deps = set()

        # For sequential processing, we need individual stages
        deps.add("field_extraction")  # Fallback for sequential mode
        deps.add("phenotype_integration")  # Always needed
        deps.add("genotype_replacement")  # Usually needed

        return deps

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Soft dependencies - these stages will run before if they exist, but are not required
        return {
            "parallel_complete_processing",  # For parallel processing mode
            "extra_column_removal",
            "dataframe_loading",  # Should run before regular DataFrame analysis
            "data_sorting",  # Ensure data is sorted for gene-aware chunking
            "annotation_config_loading",  # Optional annotation config
            "scoring_config_loading",  # Optional scoring config
            "pedigree_loading",  # Optional pedigree data
        }

    def _has_configs(self) -> bool:
        """Check if config stages exist."""
        return True  # Simplified

    def _setup_inheritance_config(self, context: PipelineContext) -> None:
        """Set up inheritance analysis configuration for chunked processing."""
        # Check if inheritance analysis should be performed
        should_calculate_inheritance = self._should_calculate_inheritance(context)

        if should_calculate_inheritance:
            logger.info("Setting up inheritance analysis configuration for chunked processing")

            # Determine inheritance mode
            inheritance_mode = context.config.get("inheritance_mode", "simple")

            # Check if scoring configuration requires Inheritance_Details
            preserve_details = self._check_if_scoring_needs_details(context)
            if preserve_details:
                logger.debug("Preserving Inheritance_Details for chunked processing scoring")

            # Store inheritance config for chunked processing
            context.config["inheritance_analysis_config"] = {
                "inheritance_mode": inheritance_mode,
                "vcf_samples": context.vcf_samples,
                "pedigree_data": context.pedigree_data,
                "use_vectorized_comp_het": not context.config.get("no_vectorized_comp_het", False),
                "preserve_details_for_scoring": preserve_details,
            }

            logger.debug(
                f"Inheritance config set up for chunked processing: mode={inheritance_mode}, "
                f"samples={len(context.vcf_samples or [])}, "
                f"pedigree={len(context.pedigree_data or {})}"
            )
        else:
            logger.debug("Inheritance analysis not needed for chunked processing")

    def _should_calculate_inheritance(self, context: PipelineContext) -> bool:
        """Check if inheritance analysis should be calculated."""
        # Use the same logic as pipeline.py
        has_ped_file = context.pedigree_data is not None and len(context.pedigree_data) > 0
        has_inheritance_mode = context.config.get("inheritance_mode")
        has_calculate_inheritance_config = context.config.get("calculate_inheritance", False)
        requires_inheritance_for_scoring = self._check_if_scoring_needs_details(context)

        should_calculate = (
            has_ped_file
            or has_inheritance_mode
            or has_calculate_inheritance_config
            or requires_inheritance_for_scoring
        )

        logger.debug(
            f"Should calculate inheritance: {should_calculate} (ped={has_ped_file}, "
            f"mode={has_inheritance_mode}, config={has_calculate_inheritance_config}, "
            f"scoring={requires_inheritance_for_scoring})"
        )
        return bool(should_calculate)

    def _check_if_scoring_needs_details(self, context: PipelineContext) -> bool:
        """Check if scoring configuration requires Inheritance_Details."""
        if not context.scoring_config:
            return False

        # Check if any scoring formulas reference inheritance variables
        if "formulas" in context.scoring_config:
            for formula_dict in context.scoring_config["formulas"]:
                for formula_name, formula_expr in formula_dict.items():
                    if isinstance(formula_expr, str) and any(
                        var in formula_expr
                        for var in [
                            "pattern",
                            "details",
                            "Inheritance_Pattern",
                            "Inheritance_Details",
                        ]
                    ):
                        logger.debug(
                            f"Scoring formula '{formula_name}' requires inheritance analysis"
                        )
                        return True

        # Check variable assignments for inheritance column dependencies
        if "variables" in context.scoring_config:
            variables = context.scoring_config["variables"]
            if any(col in variables for col in ["Inheritance_Pattern", "Inheritance_Details"]):
                logger.debug("Scoring configuration requires inheritance columns")
                return True

        return False

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process file in chunks."""
        if not context.config.get("use_chunked_processing"):
            logger.debug("Chunked processing not enabled")
            return context

        # Set up inheritance analysis configuration if needed
        self._setup_inheritance_config(context)

        # Auto-detect chunk size using ResourceManager
        from ..memory import ResourceManager

        # Get sample count from inheritance config or estimate
        inheritance_config = context.config.get("inheritance_analysis_config", {})
        vcf_samples = inheritance_config.get("vcf_samples", [])
        num_samples = len(vcf_samples) if vcf_samples else 1

        # Estimate total variants from input file (for proper chunk size calculation)
        # Use a conservative estimate based on file size if exact count unavailable
        total_variants = 100000  # Conservative default for chunk sizing

        rm = ResourceManager(config=context.config)
        chunk_size = rm.auto_chunk_size(total_variants, num_samples)

        # Determine input file
        if context.extra_columns_removed_tsv:
            input_file = context.extra_columns_removed_tsv
        elif context.phenotypes_added_tsv:
            input_file = context.phenotypes_added_tsv
        elif context.genotype_replaced_tsv:
            input_file = context.genotype_replaced_tsv
        elif context.extracted_tsv:
            input_file = context.extracted_tsv
        else:
            logger.error("No TSV file available for chunked processing")
            return context

        logger.info(f"Processing {input_file} in chunks of {chunk_size} variants")

        # Import required modules for chunk processing - done where needed

        # Find gene column
        gene_col_name = None
        # Use appropriate file opener based on compression
        if str(input_file).endswith(".gz"):
            import gzip

            def open_func(f):
                return gzip.open(f, "rt")

        else:

            def open_func(f):
                return open(f)

        with open_func(input_file) as f:
            header = f.readline().strip().split("\t")
            for gene_col in ["GENE", "Gene", "gene", "SYMBOL", "ANN[*].GENE", "EFF[*].GENE"]:
                if gene_col in header:
                    gene_col_name = gene_col
                    break

        if not gene_col_name:
            logger.warning("No GENE column found, falling back to non-chunked processing")
            return context

        # Sort file by gene column if not already sorted
        sorted_file = self._ensure_sorted_by_gene(input_file, gene_col_name, context)

        # Process chunks - use parallel processing if enabled
        # cfg["threads"] is already resolved (auto-detected or user-specified)
        threads_config = context.config.get("threads", 1)
        # Cap workers by memory constraint to avoid OOM
        memory_per_chunk_gb = rm.estimate_memory(chunk_size, num_samples)
        memory_capped_workers = rm.auto_workers(
            task_count=max(1, threads_config),
            memory_per_task_gb=memory_per_chunk_gb,
        )
        max_workers = min(threads_config, memory_capped_workers)

        parallel_chunks = context.config.get("parallel_chunks", False) or max_workers > 1

        if parallel_chunks and max_workers > 1:
            logger.info(f"Processing chunks in parallel with {max_workers} workers")
            logger.info(
                "This will significantly improve performance for large datasets "
                "with inheritance analysis"
            )
            output_chunks = self._process_chunks_parallel(
                sorted_file, gene_col_name, chunk_size, context, max_workers
            )
        else:
            logger.info("Processing chunks sequentially")
            if max_workers > 1:
                logger.info(
                    "To enable parallel chunk processing, set --parallel-chunks "
                    "or increase --threads"
                )
            output_chunks = self._process_chunks_sequential(
                sorted_file, gene_col_name, chunk_size, context
            )

        # Combine all chunks
        if output_chunks:
            logger.info(f"Combining {len(output_chunks)} processed chunks")
            # Debug: Log chunk sizes
            for i, chunk in enumerate(output_chunks):
                logger.debug(f"Chunk {i + 1} size: {len(chunk) if chunk is not None else 'None'}")

            # Filter out empty chunks
            non_empty_chunks = [
                chunk for chunk in output_chunks if chunk is not None and not chunk.empty
            ]
            if non_empty_chunks:
                context.current_dataframe = pd.concat(non_empty_chunks, ignore_index=True)
                logger.info(
                    f"Combined {len(non_empty_chunks)} non-empty chunks into "
                    f"{len(context.current_dataframe)} variants"
                )
            else:
                logger.warning("All chunks are empty - creating empty DataFrame")
                logger.warning(
                    f"Output chunks: "
                    f"{[len(chunk) if chunk is not None else 'None' for chunk in output_chunks]}"
                )
                context.current_dataframe = pd.DataFrame()
            context.config["chunked_processing_complete"] = True

            # Write chunked analysis results to a new TSV file for other stages
            chunked_output_path = context.workspace.get_intermediate_path(
                "chunked_analysis_results.tsv"
            )

            # Always compress chunked analysis results for space efficiency
            use_compression = context.config.get("gzip_intermediates", True)
            if use_compression:
                chunked_output_path = Path(str(chunked_output_path) + ".gz")
                compression = "gzip"
            else:
                compression = None

            context.current_dataframe.to_csv(
                chunked_output_path, sep="\t", index=False, compression=compression
            )

            # Update context paths for downstream stages
            context.chunked_analysis_tsv = chunked_output_path

            logger.info(
                f"Chunked processing completed: {len(context.current_dataframe)} total variants"
            )
            logger.info(f"Chunked analysis results written to: {chunked_output_path}")
        else:
            logger.warning("No chunks were processed")

        return context

    def _process_chunks_sequential(
        self, sorted_file: Path, gene_col_name: str, chunk_size: int, context: PipelineContext
    ) -> list[pd.DataFrame]:
        """Process chunks sequentially (original implementation)."""
        output_chunks = []

        for chunk_num, chunk_df in enumerate(
            self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size, context)
        ):
            logger.debug(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")
            processed_chunk = self._process_single_chunk(chunk_df, context, chunk_num)
            output_chunks.append(processed_chunk)

        return output_chunks

    def _process_chunks_parallel(
        self,
        sorted_file: Path,
        gene_col_name: str,
        chunk_size: int,
        context: PipelineContext,
        max_workers: int,
    ) -> list[pd.DataFrame]:
        """Process chunks in parallel using ThreadPoolExecutor."""
        # Pre-load all chunks into memory for parallel processing
        chunks = list(
            self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size, context)
        )
        logger.info(f"Loaded {len(chunks)} chunks for parallel processing")

        output_chunks = [None] * len(chunks)  # Pre-allocate to maintain order

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all chunks for processing
            future_to_index = {}
            for i, chunk_df in enumerate(chunks):
                future = executor.submit(self._process_single_chunk, chunk_df, context, i)
                future_to_index[future] = i

            # Collect results in order
            for future in as_completed(future_to_index):
                chunk_index = future_to_index[future]
                try:
                    processed_chunk = future.result()
                    output_chunks[chunk_index] = processed_chunk
                    logger.info(f"Completed processing chunk {chunk_index + 1}")
                except Exception as e:
                    logger.error(f"Error processing chunk {chunk_index + 1}: {e}")
                    raise

        return output_chunks

    def _process_single_chunk(
        self, chunk_df: pd.DataFrame, context: PipelineContext, chunk_num: int
    ) -> pd.DataFrame:
        """Process a single chunk with all analysis steps."""
        logger.debug(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")

        # Ensure chunk_df is not empty
        if chunk_df.empty:
            logger.warning(f"Chunk {chunk_num + 1} is empty, skipping processing")
            return chunk_df

        # Apply custom annotations if configured
        if context.annotation_configs:
            chunk_df = self._apply_chunk_annotations(chunk_df, context)

        # Apply inheritance analysis if configured
        inheritance_config = context.config.get("inheritance_analysis_config")
        if inheritance_config and not context.config.get("skip_inheritance", False):
            logger.debug(f"Applying inheritance analysis to chunk {chunk_num + 1}")
            logger.debug(f"Inheritance config: {inheritance_config}")
            chunk_df = self._apply_chunk_inheritance_analysis(chunk_df, context)
            logger.debug(f"Chunk {chunk_num + 1} after inheritance: {list(chunk_df.columns)}")
        else:
            logger.warning(
                f"Inheritance analysis not applied to chunk {chunk_num + 1}: "
                f"config={inheritance_config is not None}, "
                f"skip={context.config.get('skip_inheritance', False)}"
            )

        # Apply scoring if configured
        if context.scoring_config:
            logger.debug(f"Applying scoring to chunk {chunk_num + 1}")
            from variantcentrifuge.scoring import apply_scoring

            chunk_df = apply_scoring(chunk_df, context.scoring_config)

        # Apply analysis using the line-based interface - only if needed
        # Since inheritance and scoring are already done, we may skip this step
        # analyze_variants is causing issues with the chunk data format, so skip it

        logger.debug(
            f"Chunk {chunk_num + 1} skipping analyze_variants - inheritance and scoring complete"
        )

        if False:  # Never execute this block
            # Prepare analysis config - inheritance is handled separately above
            analysis_cfg = {
                "add_variant_id": context.config.get("add_variant_id", True),
                "perform_gene_burden": False,  # Gene burden analysis should be done after chunking
                "custom_gene_list": context.config.get("custom_gene_list"),
                "case_samples_file": context.config.get("case_samples_file"),
                "control_samples_file": context.config.get("control_samples_file"),
                "vcf_samples": context.vcf_samples,
                "aggregate_operation": context.config.get("aggregate_operation", "max"),
                "aggregate_column": context.config.get("aggregate_column"),
                "sample_values": context.config.get("sample_values", {}),
                "skip_inheritance": True,  # Already handled above
                "calculate_inheritance": False,  # Explicitly disable since handled above
                "pedigree_data": {},  # Empty dict since inheritance handled above
                "sample_list": ",".join(context.vcf_samples) if context.vcf_samples else "",
            }

            # Convert DataFrame to lines for analyze_variants function
            import io

            try:
                logger.info(
                    f"Chunk {chunk_num + 1} before analyze_variants: {len(chunk_df)} variants"
                )
                tsv_lines = chunk_df.to_csv(sep="\t", index=False)
                lines_iter = iter(tsv_lines.strip().split("\n"))

                # Call analyze_variants and collect results
                result_lines = list(analyze_variants(lines_iter, analysis_cfg))
                logger.info(
                    f"Chunk {chunk_num + 1} analyze_variants returned {len(result_lines)} lines"
                )

                # Convert back to DataFrame
                if result_lines:
                    result_text = "\n".join(result_lines)
                    chunk_df = pd.read_csv(io.StringIO(result_text), sep="\t", dtype=str)
                    logger.info(f"Chunk {chunk_num + 1} processed: {len(chunk_df)} variants")
                else:
                    logger.warning(
                        f"No results from analyze_variants for chunk {chunk_num + 1}, "
                        f"keeping original chunk with {len(chunk_df)} variants"
                    )
                    # Keep original chunk instead of losing data - inheritance and scoring
                    # already applied
            except Exception as e:
                logger.error(f"Error processing chunk {chunk_num + 1}: {e}")
                # Return original chunk if processing fails
                return chunk_df
        else:
            logger.debug(
                f"Skipping analyze_variants for chunk {chunk_num + 1} - "
                f"inheritance and scoring already complete"
            )

        logger.debug(f"Chunk {chunk_num + 1} final size: {len(chunk_df)} variants")
        logger.debug(f"Chunk {chunk_num + 1} columns: {list(chunk_df.columns)}")
        return chunk_df

    def _ensure_sorted_by_gene(
        self, input_file: Path, gene_col: str, context: PipelineContext
    ) -> Path:
        """Ensure the TSV file is sorted by gene column for gene-aware chunking.

        This method ensures that the input TSV file is sorted by the specified gene
        column, which is essential for gene-aware chunking. Gene-aware chunking keeps
        variants from the same gene together to enable proper compound heterozygous
        detection during inheritance analysis.

        Args
        ----
            input_file: Path to the input TSV file
            gene_col: Name of the gene column to sort by
            context: Pipeline context containing configuration

        Returns
        -------
            Path to the sorted TSV file (may be the same as input if already sorted)

        Note
        ----
            - Uses system sort for memory efficiency with large files
            - Handles both compressed (.gz) and uncompressed files
            - Creates cached sorted file to avoid repeated sorting
            - Preserves header line during sorting
        """
        import subprocess
        import tempfile

        # Check if already sorted
        if context.config.get("assume_sorted", False):
            return input_file

        sorted_file = context.workspace.intermediate_dir / f"{input_file.stem}.sorted.tsv"

        if sorted_file.exists():
            logger.info(f"Using existing sorted file: {sorted_file}")
            return sorted_file

        logger.info(f"Sorting {input_file} by {gene_col} column")

        # Find column index (1-based for sort command)
        # Use appropriate file opener based on compression
        if str(input_file).endswith(".gz"):
            import gzip

            def open_func(f):
                return gzip.open(f, "rt")

        else:

            def open_func(f):
                return open(f)

        with open_func(input_file) as f:
            header = f.readline().strip().split("\t")
            try:
                col_idx = header.index(gene_col) + 1  # 1-based index
            except ValueError:
                logger.error(f"Column {gene_col} not found in header")
                return input_file

        # Use system sort for memory efficiency
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
            # Write header to temp file
            with open_func(input_file) as f:
                tmp.write(f.readline())
            tmp_path = tmp.name

        try:
            # Sort data (skip header)
            sort_cmd = [
                "sort",
                "-t",
                "\t",  # Tab delimiter
                "-k",
                str(col_idx),  # Sort by gene column
                "--stable",  # Stable sort
                "-T",
                str(context.workspace.temp_dir),  # Use temp directory
            ]

            # Adaptive memory limit based on file size
            file_size_mb = input_file.stat().st_size / (1024 * 1024)
            if file_size_mb > 1000:  # Files > 1GB
                memory_limit = context.config.get("sort_memory_limit", "2G")
            elif file_size_mb > 100:  # Files > 100MB
                memory_limit = context.config.get("sort_memory_limit", "1G")
            else:
                memory_limit = context.config.get("sort_memory_limit", "512M")

            if memory_limit:
                sort_cmd.extend(["-S", memory_limit])

            logger.info(
                f"Sorting {input_file} ({file_size_mb:.1f}MB) with {memory_limit} memory limit"
            )

            # Sort command: handle both compressed and uncompressed files
            with open(tmp_path, "a") as out_f:
                # Use appropriate command for compressed vs uncompressed files
                if str(input_file).endswith(".gz"):
                    # For gzipped files: zcat input | tail -n +2 | sort ...
                    zcat_proc = subprocess.Popen(["zcat", str(input_file)], stdout=subprocess.PIPE)
                    tail_proc = subprocess.Popen(
                        ["tail", "-n", "+2"], stdin=zcat_proc.stdout, stdout=subprocess.PIPE
                    )
                    assert zcat_proc.stdout is not None
                    zcat_proc.stdout.close()
                else:
                    # For uncompressed files: tail -n +2 input | sort ...
                    tail_proc = subprocess.Popen(
                        ["tail", "-n", "+2", str(input_file)], stdout=subprocess.PIPE
                    )

                assert tail_proc.stdout is not None
                sort_proc = subprocess.Popen(sort_cmd, stdin=tail_proc.stdout, stdout=out_f)
                tail_proc.stdout.close()
                sort_proc.wait()

                if sort_proc.returncode != 0:
                    raise subprocess.CalledProcessError(sort_proc.returncode, sort_cmd)

            # Move temp file to final location
            import shutil

            shutil.move(tmp_path, sorted_file)

            logger.info(f"Sorting complete: {sorted_file}")
            return sorted_file

        except Exception as e:
            logger.error(f"Sorting failed: {e}")
            # Clean up temp file
            try:
                import os

                os.unlink(tmp_path)
            except Exception:
                pass
            # Return original file if sorting fails
            return input_file

    def _calculate_inheritance_safe_chunk_size(
        self, sample_count: int, context: Optional["PipelineContext"] = None
    ) -> int:
        """Calculate safe chunk size based on available memory to prevent inheritance skipping.

        Uses intelligent memory management to determine optimal chunk size based on system memory.
        Ensures chunks can be processed by inheritance analysis without being skipped.

        Args
        ----
            sample_count: Number of samples in the dataset
            context: Pipeline context containing memory configuration (optional)

        Returns
        -------
            Safe chunk size in number of variants
        """
        logger.debug(f"Calculating inheritance-safe chunk size for {sample_count} samples")

        try:
            # Use ResourceManager if context available
            if context and hasattr(context, "config"):
                from ..memory import ResourceManager

                rm = ResourceManager(config=context.config)
                # Use conservative estimate for total variants
                # (actual count determined during processing)
                total_variants = 100000  # Conservative estimate for chunking
                chunk_size = rm.auto_chunk_size(total_variants, sample_count)
                logger.info(
                    f"Using memory-aware chunk size: {chunk_size} variants "
                    f"for {sample_count} samples"
                )
                return chunk_size
        except Exception as e:
            logger.debug(f"Could not use ResourceManager: {e}. Falling back to improved defaults.")

        # Improved fallback defaults - much more aggressive than before
        if sample_count <= 100:
            chunk_size = 50000  # 50K variants (was 10K) - high-memory systems can handle this
        elif sample_count <= 500:
            chunk_size = 25000  # 25K variants (was 5K) - 12.5M cells, reasonable for modern systems
        elif sample_count <= 1000:
            chunk_size = 15000  # 15K variants (was 3K) - 15M cells, good balance
        elif sample_count <= 2500:
            chunk_size = 8000  # 8K variants (was 2.5K) - 20M cells, within modern limits
        elif sample_count <= 5000:
            chunk_size = 6000  # 6K variants - 30M cells, optimized for large cohorts
        else:
            # For very large cohorts (>5000 samples), use memory-based calculation
            # Target ~50M cells per chunk (reasonable for 256GB+ systems)
            target_cells = 50_000_000
            chunk_size = max(1000, target_cells // sample_count)

        # Apply reasonable bounds
        chunk_size = max(1000, min(chunk_size, 100000))  # Between 1K and 100K variants

        logger.info(
            f"Using inheritance-safe chunk size: {chunk_size} variants for {sample_count} samples "
            f"(~{chunk_size * sample_count / 1_000_000:.1f}M cells)"
        )
        return chunk_size

    def _read_tsv_in_gene_chunks(
        self,
        filepath: Path,
        gene_column: str,
        chunksize: int,
        context: Optional["PipelineContext"] = None,
    ):
        """Read TSV file in gene-aware chunks with configurable gene size limits."""
        import pandas as pd

        # Calculate inheritance-safe chunk size if context and VCF samples are available
        if context and hasattr(context, "vcf_samples") and context.vcf_samples:
            sample_count = len(context.vcf_samples)
            safe_chunk_size = self._calculate_inheritance_safe_chunk_size(sample_count, context)

            # Use the smaller of user-specified chunk size or safe chunk size
            if chunksize > safe_chunk_size:
                logger.warning(
                    f"User-specified chunk size ({chunksize}) exceeds inheritance-safe limit "
                    f"({safe_chunk_size}) for {sample_count} samples. Using safe chunk size."
                )
                chunksize = safe_chunk_size
            else:
                logger.info(
                    f"Using user-specified chunk size ({chunksize}) which is within "
                    f"inheritance-safe limit ({safe_chunk_size}) for {sample_count} samples"
                )
        else:
            logger.debug("No VCF samples available, using original chunk size")

        logger.info(
            f"Reading {filepath} in gene-aware chunks (column: {gene_column}, size: {chunksize})"
        )

        # Create reader for chunked reading
        reader = pd.read_csv(
            filepath,
            sep="\t",
            chunksize=chunksize,
            low_memory=False,
            dtype=str,  # Read all as strings to preserve data
        )

        gene_buffer = pd.DataFrame()
        chunks_yielded = 0

        # Configurable gene size limits to prevent memory explosion
        # Get from context config if available, otherwise use instance attribute or default
        max_variants_per_gene = 10000  # Default
        if context and hasattr(context, "config"):
            max_variants_per_gene = context.config.get(
                "max_variants_per_gene", max_variants_per_gene
            )
        elif hasattr(self, "max_variants_per_gene"):
            max_variants_per_gene = self.max_variants_per_gene

        # Adjust buffer sizes based on inheritance-safe chunk size
        # For inheritance-safe chunking, we want tighter limits to avoid memory issues
        max_buffer_size = min(
            chunksize * 5, max_variants_per_gene
        )  # Reduced from 20x to 5x for better memory control
        emergency_threshold = min(chunksize * 10, max_variants_per_gene)  # Reduced from 50x to 10x

        for chunk in reader:
            # Combine with buffer
            if gene_buffer.empty:
                gene_buffer = chunk
            else:
                gene_buffer = pd.concat([gene_buffer, chunk], ignore_index=True)

            # Check buffer size and warn about large genes
            if len(gene_buffer) > max_buffer_size:
                logger.warning(
                    f"Gene buffer has grown to {len(gene_buffer)} rows. "
                    f"This may indicate a gene with too many variants. "
                    f"Consider increasing chunk size or checking data sorting."
                )
                # Force yield to prevent memory issues
                if len(gene_buffer) > emergency_threshold:  # Configurable emergency threshold
                    logger.error(
                        f"Gene buffer exceeded emergency threshold ({emergency_threshold} rows, "
                        f"max_variants_per_gene={max_variants_per_gene}). "
                        f"Forcing chunk yield to prevent memory exhaustion."
                    )
                    # Find any gene boundary to split at
                    gene_values = gene_buffer[gene_column].values
                    for i in range(chunksize, len(gene_values)):
                        if gene_values[i] != gene_values[i - 1]:
                            # Found a gene boundary, split here
                            yield_df = gene_buffer.iloc[:i].copy()
                            yield yield_df
                            chunks_yielded += 1
                            gene_buffer = gene_buffer.iloc[i:].reset_index(drop=True)
                            break
                    else:
                        # No gene boundary found, yield the entire buffer
                        logger.warning("No gene boundaries found, yielding entire buffer")
                        yield gene_buffer.copy()
                        chunks_yielded += 1
                        gene_buffer = pd.DataFrame()

            # Find gene boundaries for normal processing
            if len(gene_buffer) >= chunksize:
                # Find where genes change
                gene_values = gene_buffer[gene_column].values
                gene_change_indices = [0]
                current_gene = gene_values[0]

                for i in range(1, len(gene_values)):
                    if gene_values[i] != current_gene:
                        gene_change_indices.append(i)
                        current_gene = gene_values[i]

                # Find split point after chunksize
                split_index = None
                for idx in gene_change_indices:
                    if idx >= chunksize:
                        split_index = idx
                        break

                if split_index is None:
                    # All rows belong to same gene or last gene extends beyond chunksize
                    # Keep accumulating unless buffer is getting too large
                    continue
                else:
                    # Yield complete genes
                    yield_df = gene_buffer.iloc[:split_index].copy()
                    # Log gene distribution in chunk
                    if logger.isEnabledFor(logging.DEBUG):
                        genes_in_chunk = yield_df[gene_column].nunique()
                        logger.debug(
                            f"Yielding chunk with {len(yield_df)} variants across "
                            f"{genes_in_chunk} genes"
                        )

                    yield yield_df
                    chunks_yielded += 1

                    # Keep remaining for next iteration
                    gene_buffer = gene_buffer.iloc[split_index:].reset_index(drop=True)

        # Yield any remaining data
        if not gene_buffer.empty:
            if logger.isEnabledFor(logging.DEBUG):
                genes_in_final_chunk = gene_buffer[gene_column].nunique()
                logger.debug(
                    f"Yielding final chunk with {len(gene_buffer)} variants across "
                    f"{genes_in_final_chunk} genes"
                )
            yield gene_buffer
            chunks_yielded += 1

        logger.info(f"Completed reading {chunks_yielded} gene-aware chunks")

    def _apply_chunk_annotations(
        self, chunk_df: pd.DataFrame, context: PipelineContext
    ) -> pd.DataFrame:
        """Apply custom annotations to a chunk."""
        try:
            from ..annotator import (
                annotate_dataframe_with_features,
                load_custom_features,
            )

            # Create config dict for load_custom_features
            annotation_cfg = {
                "annotate_bed_files": context.annotation_configs.get("bed_files", []),
                "annotate_gene_lists": context.annotation_configs.get("gene_lists", []),
                "annotate_json_genes": context.annotation_configs.get("json_genes", []),
                "json_gene_mapping": context.annotation_configs.get("json_mapping", ""),
                "json_genes_as_columns": context.config.get("json_genes_as_columns", False),
            }
            features = load_custom_features(annotation_cfg)

            if features:
                chunk_df = annotate_dataframe_with_features(chunk_df, features)
                logger.debug(f"Applied {len(features)} annotation features to chunk")

            return chunk_df
        except Exception as e:
            logger.error(f"Error applying annotations to chunk: {e}")
            return chunk_df

    def _apply_chunk_inheritance_analysis(
        self, chunk_df: pd.DataFrame, context: PipelineContext
    ) -> pd.DataFrame:
        """Apply inheritance analysis to a chunk with gene-aware processing.

        This method processes a chunk of variant data and applies inheritance analysis
        to determine inheritance patterns. It handles sample column creation from GT
        data, performs inheritance analysis, and cleans up temporary columns.

        Args
        ----
            chunk_df: DataFrame containing variant data for a chunk
            context: Pipeline context containing configuration and samples

        Returns
        -------
            DataFrame with inheritance analysis results and cleaned columns

        Note
        ----
            - Creates individual sample columns from GT column if needed
            - Applies inheritance analysis using the configured mode
            - Removes temporary sample columns after analysis to prevent huge files
            - Handles memory safety checks for large datasets
        """
        try:
            from ..inheritance.analyzer import (
                analyze_inheritance,
                process_inheritance_output,
            )

            # Extract config parameters
            inheritance_config = context.config.get("inheritance_analysis_config", {})
            vcf_samples = inheritance_config.get("vcf_samples", [])
            pedigree_data = inheritance_config.get("pedigree_data", {})
            use_vectorized_comp_het = inheritance_config.get("use_vectorized_comp_het", True)
            inheritance_mode = inheritance_config.get("inheritance_mode", "simple")
            preserve_details_for_scoring = inheritance_config.get(
                "preserve_details_for_scoring", False
            )

            if not vcf_samples:
                logger.warning("No VCF samples available for chunk inheritance analysis")
                return chunk_df

            # Use ResourceManager for memory safety check
            from ..memory import ResourceManager

            rm = ResourceManager(config=context.config)
            estimated_memory_gb = rm.estimate_memory(len(chunk_df), len(vcf_samples))
            safe_memory_gb = rm.memory_gb * rm.memory_safety_factor

            logger.debug(
                f"Chunk memory: {estimated_memory_gb:.2f}GB estimated, "
                f"{safe_memory_gb:.2f}GB available"
            )

            # Check if this chunk should be processed
            force_processing = context.config.get("force_inheritance_processing", False)
            should_process = estimated_memory_gb <= safe_memory_gb or force_processing

            if not should_process:
                logger.error(
                    f"MEMORY LIMIT: Skipping inheritance analysis for chunk with "
                    f"{len(chunk_df)} variants x {len(vcf_samples)} samples "
                    f"({estimated_memory_gb:.2f}GB exceeds {safe_memory_gb:.2f}GB limit). "
                    f"Consider increasing --max-memory-gb or using --force-inheritance-processing."
                )
                # Add placeholder inheritance columns to maintain schema
                chunk_df["Inheritance_Pattern"] = "memory_limit_exceeded"
                if preserve_details_for_scoring:
                    chunk_df["Inheritance_Details"] = '{"error": "memory_limit_exceeded"}'
                return chunk_df

            logger.info(
                f"Starting inheritance analysis for chunk: {len(chunk_df)} variants, "
                f"{len(vcf_samples)} samples, inheritance mode: {inheritance_mode}"
            )

            # Track timing for sample column preparation
            import time

            # CRITICAL: Create individual sample columns from GT column
            # This is required for inheritance analysis to work correctly
            if "GT" in chunk_df.columns and len(chunk_df) > 0:
                logger.debug("Creating sample columns from GT column for inheritance analysis")

                prep_start = time.time()

                # Use the intelligent sample column creation function
                sample_column_method = context.config.get("sample_column_creation_method", "auto")
                chunk_df = create_sample_columns_from_gt_intelligent(
                    df=chunk_df,
                    vcf_samples=vcf_samples,
                    separator=context.config.get("separator", ";"),
                    snpsift_sep=context.config.get("extract_fields_separator", ","),
                    method=sample_column_method,
                )

                prep_time = time.time() - prep_start
                if prep_time > 0.1:  # Only log if significant
                    logger.debug(f"Sample column preparation took {prep_time:.2f}s")
            else:
                logger.warning(
                    "No GT column found in chunk, inheritance analysis may not work correctly"
                )

            # Apply inheritance analysis with parallel processing support
            analysis_start = time.time()
            threads = context.config.get("threads", 1)
            min_variants_for_parallel = context.config.get(
                "min_variants_for_parallel_inheritance", 100
            )

            if threads > 1 and len(chunk_df) >= min_variants_for_parallel:
                # Use parallel analyzer for better performance
                logger.info(f"Using parallel inheritance analyzer for chunk with {threads} workers")

                from ..inheritance.parallel_analyzer import analyze_inheritance_parallel

                comp_het_start = time.time()

                chunk_df = analyze_inheritance_parallel(
                    df=chunk_df,
                    sample_list=vcf_samples,
                    pedigree_data=pedigree_data,
                    use_vectorized_comp_het=use_vectorized_comp_het,
                    n_workers=threads,
                    min_variants_for_parallel=min_variants_for_parallel,
                )

                # Record compound het timing (this is a major component)
                comp_het_time = time.time() - comp_het_start
                if comp_het_time > 1.0:  # Only record if significant
                    logger.debug(f"Parallel compound het analysis took {comp_het_time:.2f}s")
            else:
                # Use sequential analyzer for small chunks or single-threaded
                chunk_df = analyze_inheritance(
                    df=chunk_df,
                    pedigree_data=pedigree_data,
                    sample_list=vcf_samples,
                    use_vectorized_comp_het=use_vectorized_comp_het,
                )

            analysis_time = time.time() - analysis_start
            logger.info(
                f"Inheritance analysis completed in {analysis_time:.2f}s "
                f"for {len(chunk_df)} variants"
            )

            # Process inheritance output based on mode and scoring requirements
            output_start = time.time()
            if preserve_details_for_scoring:
                logger.debug("Preserving Inheritance_Details for chunk scoring")

            chunk_df = process_inheritance_output(
                chunk_df,
                inheritance_mode,
                preserve_details_for_scoring=preserve_details_for_scoring,
            )

            output_time = time.time() - output_start
            if output_time > 0.1:  # Only log if significant
                logger.debug(f"Output processing took {output_time:.2f}s")

            # CRITICAL: Clean up sample columns after inheritance analysis
            # The sample columns were created from GT for inheritance analysis
            # but should not be in the final output as they make the file huge
            # Use unified cleanup function for consistency with non-chunked processing
            cleanup_start = time.time()
            chunk_df = cleanup_sample_columns(
                chunk_df,
                vcf_samples,
                preserve_columns=["GT", "Inheritance_Pattern", "Inheritance_Details"],
            )

            cleanup_time = time.time() - cleanup_start
            if cleanup_time > 0.1:  # Only log if significant
                logger.debug(f"Column cleanup took {cleanup_time:.2f}s")

            # Log successful completion with inheritance pattern summary
            inheritance_patterns = chunk_df["Inheritance_Pattern"].value_counts()
            logger.info(
                f"Chunk inheritance analysis SUCCESS: {len(chunk_df)} variants processed. "
                f"Patterns: {dict(inheritance_patterns)}"
            )

            return chunk_df

        except Exception as e:
            return handle_inheritance_analysis_error(
                chunk_df,
                e,
                preserve_details_for_scoring=inheritance_config.get(
                    "preserve_details_for_scoring", False
                ),
                context_description="chunk inheritance analysis",
            )

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        # Check in priority order for available TSV files
        for attr in [
            "extra_columns_removed_tsv",
            "phenotypes_added_tsv",
            "genotype_replaced_tsv",
            "extracted_tsv",
        ]:
            if hasattr(context, attr):
                tsv_file = getattr(context, attr)
                if tsv_file and tsv_file.exists():
                    return [tsv_file]
        # Fallback to context.data
        if hasattr(context, "data") and context.data:
            return [context.data]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        # ChunkedAnalysisStage processes data in memory, no direct file outputs
        return []


class ParallelAnalysisOrchestrator(Stage):
    """Orchestrate parallel analysis of independent genes."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "parallel_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Run analysis in parallel by gene"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Can run after basic data loading
        return {"dataframe_loading", "custom_annotation"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        # This stage manages its own parallelism
        return False

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Run analysis in parallel by gene."""
        if not context.config.get("parallel_analysis"):
            logger.debug("Parallel analysis not enabled")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for parallel analysis")
            return context

        threads = context.config.get("threads", 1)
        if threads <= 1:
            logger.debug("Thread count <= 1, skipping parallel analysis")
            return context

        gene_column = context.config.get("gene_column", "GENE")

        # Find the gene column if not specified
        if gene_column not in df.columns:
            for col in ["GENE", "Gene", "gene", "SYMBOL", "ANN[*].GENE", "EFF[*].GENE"]:
                if col in df.columns:
                    gene_column = col
                    break
            else:
                logger.warning("No gene column found, cannot perform parallel analysis by gene")
                return context

        unique_genes = df[gene_column].unique()
        num_genes = len(unique_genes)

        logger.info(f"Running parallel analysis for {num_genes} genes using {threads} threads")

        # Import required modules
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Split DataFrame by gene
        gene_groups = df.groupby(gene_column, observed=True)

        # Prepare analysis configuration
        analysis_cfg = {
            "add_variant_id": context.config.get("add_variant_id", True),
            "perform_gene_burden": context.config.get("perform_gene_burden", False),
            "custom_gene_list": context.config.get("custom_gene_list"),
            "case_samples_file": context.config.get("case_samples_file"),
            "control_samples_file": context.config.get("control_samples_file"),
            "vcf_samples": context.vcf_samples,
            "aggregate_operation": context.config.get("aggregate_operation", "max"),
            "aggregate_column": context.config.get("aggregate_column"),
            "sample_values": context.config.get("sample_values", {}),
        }

        # Track successful and failed genes
        successful_results = []
        failed_genes = []

        # Process genes in parallel
        with ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit tasks
            futures = {}
            for gene_name, gene_df in gene_groups:
                future = executor.submit(
                    self._process_gene_group,
                    gene_name,
                    gene_df,
                    context.pedigree_data,
                    context.scoring_config,
                    analysis_cfg,
                    context.config.get("skip_analysis", False),
                )
                futures[future] = gene_name

            # Process results as they complete
            try:
                for future in as_completed(futures):
                    gene_name = futures[future]
                    try:
                        result_df = future.result()
                        if result_df is not None:
                            successful_results.append(result_df)
                            logger.debug(f"Successfully processed gene: {gene_name}")
                    except Exception as e:
                        logger.error(f"Failed to process gene {gene_name}: {e!s}")
                        failed_genes.append(gene_name)
                        # Don't cancel other tasks, continue processing

            except KeyboardInterrupt:
                logger.warning("Parallel processing interrupted by user")
                # Cancel remaining futures
                for future in futures:
                    future.cancel()
                raise

        # Merge results
        if successful_results:
            logger.info(
                f"Merging results from {len(successful_results)} successfully processed genes"
            )
            context.current_dataframe = pd.concat(successful_results, ignore_index=True)

            # Restore original column order
            original_columns = df.columns.tolist()
            current_columns = context.current_dataframe.columns.tolist()
            # Add any new columns at the end
            new_columns = [col for col in current_columns if col not in original_columns]
            ordered_columns = original_columns + new_columns
            # Only reorder columns that exist in the result
            ordered_columns = [col for col in ordered_columns if col in current_columns]
            context.current_dataframe = context.current_dataframe[ordered_columns]

            if failed_genes:
                logger.warning(
                    f"Failed to process {len(failed_genes)} genes: {', '.join(failed_genes[:10])}"
                )
                if len(failed_genes) > 10:
                    logger.warning(f"... and {len(failed_genes) - 10} more")
        else:
            logger.error("No genes were successfully processed")
            context.current_dataframe = df  # Keep original data

        logger.info(
            f"Parallel analysis completed: {len(successful_results)}/{num_genes} genes processed"
        )
        return context

    @staticmethod
    def _process_gene_group(
        gene_name: str,
        gene_df: pd.DataFrame,
        pedigree_data: dict | None,
        scoring_config: dict | None,
        analysis_cfg: dict,
        skip_analysis: bool,
    ) -> pd.DataFrame:
        """Process a single gene's variants (runs in worker process)."""
        import io

        from variantcentrifuge.analyze_variants import analyze_variants
        from variantcentrifuge.scoring import apply_scoring

        # Make a copy to avoid modifying the original
        result_df = gene_df.copy()

        try:
            # Apply scoring if configured
            if scoring_config:
                result_df = apply_scoring(result_df, scoring_config)

            # Apply analysis if not skipped
            if not skip_analysis:
                # Convert DataFrame to lines for analyze_variants function
                tsv_lines = result_df.to_csv(sep="\t", index=False)
                lines_iter = iter(tsv_lines.strip().split("\n"))

                # Call analyze_variants and collect results
                result_lines = list(analyze_variants(lines_iter, analysis_cfg))

                # Convert back to DataFrame
                if result_lines:
                    result_text = "\n".join(result_lines)
                    result_df = pd.read_csv(io.StringIO(result_text), sep="\t", dtype=str)

            return result_df

        except Exception as e:
            # Log error and re-raise
            import logging

            logger = logging.getLogger(__name__)
            logger.error(f"Error processing gene {gene_name}: {e!s}")
            raise
