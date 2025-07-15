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
from typing import Any, Dict, List, Optional, Set

import pandas as pd

from ..analyze_variants import analyze_variants
from ..annotator import annotate_dataframe_with_features, load_custom_features
from ..filters import filter_final_tsv_by_genotype
from ..gene_burden import perform_gene_burden_analysis
from ..inheritance import analyze_inheritance
from ..pipeline_core import PipelineContext, Stage
from ..scoring import apply_scoring
from ..stats_engine import StatsEngine

logger = logging.getLogger(__name__)


def create_sample_columns_from_gt(
    df: pd.DataFrame, 
    vcf_samples: List[str], 
    separator: str = ";", 
    snpsift_sep: str = ","
) -> pd.DataFrame:
    """Create individual sample columns from GT column data.
    
    This function extracts genotype information from the GT column and creates
    individual sample columns for inheritance analysis. It handles both replaced
    genotype format (e.g., "Sample1(0/1);Sample2(0/0)") and SnpSift format
    (e.g., "0/1,0/0,1/1").
    
    Args:
        df: DataFrame containing GT column
        vcf_samples: List of sample IDs to create columns for
        separator: Separator used in replaced genotype format (default: ";")
        snpsift_sep: Separator used in SnpSift format (default: ",")
        
    Returns:
        DataFrame with individual sample columns added
        
    Raises:
        ValueError: If GT column is missing or empty
    """
    if "GT" not in df.columns:
        raise ValueError("GT column not found in DataFrame")
    
    if len(df) == 0:
        logger.warning("Empty DataFrame provided to create_sample_columns_from_gt")
        return df
    
    # Check if sample columns already exist
    sample_columns_exist = any(
        sample_id in df.columns for sample_id in vcf_samples
    )
    
    if sample_columns_exist:
        logger.debug("Sample columns already exist, skipping creation")
        return df
    
    # Get the first GT value to determine format
    first_gt = str(df.iloc[0]["GT"]) if len(df) > 0 else ""
    
    # Check if genotype replacement has already been done
    # (format: "Sample1(0/1);Sample2(0/0)")
    if "(" in first_gt and ")" in first_gt:
        logger.debug("GT column contains replaced genotypes, extracting sample genotypes")
        
        # Pre-create all sample columns data
        sample_data = {sample_id: [] for sample_id in vcf_samples}
        
        # Parse replaced genotype format
        for idx, row in df.iterrows():
            gt_value = str(row["GT"])
            # Initialize all samples with missing genotype
            row_genotypes = {sample_id: "./." for sample_id in vcf_samples}
            
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
        logger.debug(
            f"Extracting sample genotypes from GT column using separator '{snpsift_sep}'"
        )
        
        # Pre-create all sample columns data
        sample_data = {sample_id: [] for sample_id in vcf_samples}
        
        # Extract genotypes for each row
        for idx, row in df.iterrows():
            gt_value = str(row["GT"])
            if gt_value and gt_value != "NA" and gt_value != "nan":
                genotypes = gt_value.split(snpsift_sep)
                if len(genotypes) != len(vcf_samples):
                    logger.warning(
                        f"Row {idx}: Expected {len(vcf_samples)} genotypes "
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
    context_description: str = "inheritance analysis"
) -> pd.DataFrame:
    """Handle errors in inheritance analysis with consistent error recovery.
    
    This function provides unified error handling for both chunked and non-chunked
    processing. It ensures the DataFrame schema is maintained even when errors occur.
    
    Args:
        df: DataFrame that encountered an error
        error: The exception that occurred
        preserve_details_for_scoring: Whether to add Inheritance_Details column
        context_description: Description of the processing context for logging
        
    Returns:
        DataFrame with error-state inheritance columns added
        
    Note:
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
    df: pd.DataFrame, 
    vcf_samples: List[str], 
    preserve_columns: List[str] = None
) -> pd.DataFrame:
    """Clean up sample columns after inheritance analysis.
    
    This function removes temporary sample columns that were created for inheritance
    analysis while preserving essential columns. It provides a unified approach
    for both chunked and non-chunked processing.
    
    Args:
        df: DataFrame containing sample columns to clean up
        vcf_samples: List of sample IDs that should be removed
        preserve_columns: List of column names to preserve even if they match sample names
                         (default: ["GT", "Inheritance_Pattern", "Inheritance_Details"])
        
    Returns:
        DataFrame with sample columns cleaned up
        
    Note:
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
        col for col in df.columns 
        if col in vcf_samples and col not in preserve_columns
    ]
    
    if sample_columns_to_remove:
        logger.debug(f"Removing {len(sample_columns_to_remove)} sample columns after inheritance analysis")
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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # In sequential mode, depends on field_extraction
        # In parallel mode, ParallelCompleteProcessingStage handles field extraction internally
        # We return an empty set here and check for data availability in _process()
        # This allows the stage to work in both sequential and parallel modes
        return set()

    @property
    def soft_dependencies(self) -> Set[str]:
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

        # Check if we should use chunked processing
        if self._should_use_chunks(context, input_file):
            logger.info("File too large, will use chunked processing (checkpoint skip)")
            context.config["use_chunked_processing"] = True
            # Don't load full DataFrame - chunked stages will handle it
            return context

        # Load full DataFrame
        logger.info(f"Restoring DataFrame from {input_file} (checkpoint skip)")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        # Use quoting=3 (QUOTE_NONE) to handle data with quotes
        # Use low_memory=False to avoid dtype inference issues
        df = pd.read_csv(
            input_file,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_values=[""],
            compression=compression,
            quoting=3,  # QUOTE_NONE - don't treat quotes specially
            low_memory=False,
            on_bad_lines="warn",  # Warn about problematic lines instead of failing
        )

        context.current_dataframe = df
        logger.info(f"Restored {len(df)} variants into DataFrame after checkpoint skip")

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
                            return tsv_file
                else:
                    raise ValueError(
                        f"DataFrameLoadingStage requires a TSV file, but got VCF: {input_file}. "
                        f"No TSV files found in context."
                    )
            return input_file

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load data into DataFrame or prepare chunked processing."""
        # Find the input file
        input_file = self._find_input_file(context)

        # Check if we should use chunked processing
        if self._should_use_chunks(context, input_file):
            logger.info("File too large, will use chunked processing")
            context.config["use_chunked_processing"] = True
            # Don't load full DataFrame - chunked stages will handle it
            return context

        # Load full DataFrame
        logger.info(f"Loading data from {input_file}")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        # Use quoting=3 (QUOTE_NONE) to handle data with quotes
        # Use low_memory=False to avoid dtype inference issues
        df = pd.read_csv(
            input_file,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_values=[""],
            compression=compression,
            quoting=3,  # QUOTE_NONE - don't treat quotes specially
            low_memory=False,
            on_bad_lines="warn",  # Warn about problematic lines instead of failing
        )

        context.current_dataframe = df
        logger.info(f"Loaded {len(df)} variants into DataFrame")

        return context

    def _should_use_chunks(self, context: PipelineContext, file_path: Path) -> bool:
        """Determine if chunked processing should be used."""
        # Check if chunked processing is explicitly disabled
        if context.config.get("no_chunked_processing"):
            return False

        # Check if chunked processing is explicitly forced
        if context.config.get("force_chunked_processing"):
            return True

        # Check explicit chunking request
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
                        f"Memory estimation: {estimated_variants} variants × {sample_count} samples"
                    )
                    return True
            else:
                # Fallback to file size only
                threshold_mb = context.config.get("chunk_threshold_mb", 500)
                return file_size_mb > threshold_mb

        except Exception as e:
            logger.debug(f"Error in chunking decision: {e}")
            return False

        return False

    def get_input_files(self, context: PipelineContext) -> List[Path]:
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

    def get_output_files(self, context: PipelineContext) -> List[Path]:
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
    def dependencies(self) -> Set[str]:
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

    def get_input_files(self, context: PipelineContext) -> List[Path]:
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

    def get_output_files(self, context: PipelineContext) -> List[Path]:
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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depend on dataframe_loading as a hard requirement
        # custom_annotation is optional and will run before if present
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> Set[str]:
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

        # Check if using chunked processing
        if context.config.get("use_chunked_processing"):
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
        pedigree_data = context.pedigree_data

        if not vcf_samples:
            logger.warning("No VCF samples available for inheritance analysis")
            return context

        # Ensure vcf_samples is a deterministically ordered list
        if isinstance(vcf_samples, set):
            # If somehow still a set, convert to sorted list for deterministic order
            vcf_samples = sorted(list(vcf_samples))
            logger.warning(
                "VCF samples was a set - converted to sorted list for deterministic order"
            )

        # Memory safety check - avoid creating massive matrices
        estimated_memory_cells = len(df) * len(vcf_samples)
        max_memory_cells = context.config.get(
            "max_inheritance_memory_cells", 100_000_000
        )  # 100M cells

        if estimated_memory_cells > max_memory_cells:
            logger.warning(
                f"Dataset too large for in-memory inheritance analysis: "
                f"{len(df)} variants × {len(vcf_samples)} samples = "
                f"{estimated_memory_cells:,} cells. "
                f"Forcing chunked processing (limit: {max_memory_cells:,} cells)"
            )
            # Force chunked processing and delegate to chunked analysis
            context.config["use_chunked_processing"] = True
            context.config["force_chunked_processing"] = True
            # Check if scoring configuration requires Inheritance_Details
            preserve_details = self._check_if_scoring_needs_details(context)
            if preserve_details:
                logger.info(
                    "Preserving Inheritance_Details for large dataset chunked processing scoring"
                )

            context.config["inheritance_analysis_config"] = {
                "inheritance_mode": inheritance_mode,
                "vcf_samples": vcf_samples,
                "pedigree_data": pedigree_data,
                "use_vectorized_comp_het": not context.config.get("no_vectorized_comp_het", False),
                "preserve_details_for_scoring": preserve_details,
            }
            return context

        logger.info(
            f"Calculating inheritance patterns ({inheritance_mode} mode) "
            f"for {len(vcf_samples)} samples"
        )

        # Prepare sample columns if GT field exists and contains multiple samples
        if "GT" in df.columns:
            logger.debug("Preparing sample columns from GT field for inheritance analysis")
            prep_start = self._start_subtask("sample_column_preparation")

            # Use the unified sample column creation function
            df = create_sample_columns_from_gt(
                df=df,
                vcf_samples=vcf_samples,
                separator=context.config.get("separator", ";"),
                snpsift_sep=context.config.get("extract_fields_separator", ",")
            )

            self._end_subtask("sample_column_preparation", prep_start)

        # Apply inheritance analysis with error handling
        analysis_start = self._start_subtask("inheritance_calculation")

        try:
            # Check if we should use parallel processing
            threads = context.config.get("threads", 1)
            min_variants_for_parallel = context.config.get("min_variants_for_parallel_inheritance", 100)

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
                    pedigree_data=pedigree_data,
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
                    pedigree_data=pedigree_data,
                    use_vectorized_comp_het=not context.config.get("no_vectorized_comp_het", False),
                )

        except Exception as e:
            # Check if scoring configuration requires Inheritance_Details
            preserve_details = self._check_if_scoring_needs_details(context)
            df = handle_inheritance_analysis_error(
                df, 
                e, 
                preserve_details_for_scoring=preserve_details,
                context_description="non-chunked inheritance analysis"
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
                    if isinstance(formula_expr, str):
                        # Check if formula references details or extracted fields from details
                        if any(var in formula_expr for var in ["details", "segregation_p_value"]):
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

    def get_input_files(self, context: PipelineContext) -> List[Path]:
        """Return input files for checkpoint tracking."""
        input_files = []
        # Include PED file if available for inheritance analysis
        if hasattr(context, "ped_file") and context.ped_file:
            input_files.append(Path(context.ped_file))
        return input_files

    def get_output_files(self, context: PipelineContext) -> List[Path]:
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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depend on dataframe_loading as a hard requirement
        # Other stages (custom_annotation, inheritance_analysis) are optional
        # The pipeline runner will ensure correct ordering based on what stages are present
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after annotations and inheritance if they exist
        return {"custom_annotation", "inheritance_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply variant scoring."""
        if not context.scoring_config:
            logger.debug("No scoring configuration loaded")
            return context

        # Check if using chunked processing
        if context.config.get("use_chunked_processing"):
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
    def dependencies(self) -> Set[str]:
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
            genotype_modes = set(g.strip() for g in genotype_filter.split(",") if g.strip())

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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depends on dataframe being loaded
        # Statistics can be run on whatever columns are available
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> Set[str]:
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

    def _write_statistics(self, stats: Dict[str, Any], output_file: str) -> None:
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

    def _get_default_stats_config(self) -> Dict[str, Any]:
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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on having a DataFrame with variants
        deps = {"dataframe_loading"}
        return deps

    @property
    def soft_dependencies(self) -> Set[str]:
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
            logger.warning(
                "No inheritance columns found in input DataFrame - "
                "they were lost before VariantAnalysisStage"
            )

        # Prepare config for analyze_variants
        analysis_config = {
            "vcf_sample_names": context.vcf_samples or [],
            "sample_list": ",".join(context.vcf_samples or []),
            "no_db_links": context.config.get("no_db_links", False),
            "base_name": context.workspace.base_name,
            "reference": context.config.get("reference", "GRCh37"),
        }

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
                return open(f, "r")

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
                                logger.error(
                                    f"Failed to preserve inheritance columns: {inheritance_cols}"
                                )
                    except Exception as e:
                        logger.error(f"Error merging columns: {e}")
                        logger.warning("Inheritance columns may be lost due to merge failure")

                    # Reorder columns to put VAR_ID first if it exists
                    if "VAR_ID" in analysis_df.columns:
                        cols = analysis_df.columns.tolist()
                        cols.remove("VAR_ID")
                        cols = ["VAR_ID"] + cols
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
                if inheritance_cols:
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
                logger.error(
                    "INHERITANCE COLUMNS LOST - Final DataFrame missing inheritance columns!"
                )
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
                logger.warning("Unchanged DataFrame also missing inheritance columns")

        # Clean up temp file
        if temp_tsv.exists():
            temp_tsv.unlink()

        return context

    def get_input_files(self, context: PipelineContext) -> List[Path]:
        """Return input files for checkpoint tracking."""
        # VariantAnalysisStage works with DataFrame in memory, no specific input files
        return []

    def get_output_files(self, context: PipelineContext) -> List[Path]:
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
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on having a DataFrame with variants and case/control samples
        return {"dataframe_loading", "sample_config_loading"}

    @property
    def soft_dependencies(self) -> Set[str]:
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

        if not case_samples or not control_samples:
            logger.warning("Case/control samples not defined for gene burden analysis")
            return context

        df = context.current_dataframe
        if df is None or df.empty:
            logger.warning("No DataFrame loaded for gene burden analysis")
            return context

        # Check if DataFrame has required GENE column
        if "GENE" not in df.columns:
            logger.error("DataFrame missing required 'GENE' column for gene burden analysis")
            logger.debug(f"Available columns: {list(df.columns)}")
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

        df_with_counts = assign_case_control_counts(
            df=df,
            case_samples=set(case_samples),
            control_samples=set(control_samples),
            all_samples=all_vcf_samples,  # This should be ALL samples in VCF
        )

        # Perform gene burden analysis on prepared DataFrame
        burden_results = perform_gene_burden_analysis(df=df_with_counts, cfg=burden_config)

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

    def get_input_files(self, context: PipelineContext) -> List[Path]:
        """Return input files for checkpoint tracking."""
        # GeneBurdenAnalysisStage works with DataFrame in memory, no specific input files
        return []

    def get_output_files(self, context: PipelineContext) -> List[Path]:
        """Return output files for checkpoint tracking."""
        # Return gene burden output file if configured
        if hasattr(context, "config") and context.config.get("gene_burden_output"):
            burden_output = context.config["gene_burden_output"]
            if Path(burden_output).exists():
                return [Path(burden_output)]
        return []


class ChunkedAnalysisStage(Stage):
    """Process large files in chunks with all analysis steps."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "chunked_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Process large files in memory-efficient chunks"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on configuration and processing stages but not full DataFrame loading
        # This stage only runs when use_chunked_processing is True
        deps = set()

        # Need basic processing to be complete
        # For parallel processing, we need the complete processing stage
        # For sequential processing, we need individual stages
        deps.add("parallel_complete_processing")  # Will be ignored if not present
        deps.add("field_extraction")  # Fallback for sequential mode
        deps.add("phenotype_integration")  # Always needed
        deps.add("genotype_replacement")  # Usually needed

        return deps

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Soft dependencies - these stages will run before if they exist, but are not required
        return {
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
        # Use the same logic as pipeline_refactored.py
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
        return should_calculate

    def _check_if_scoring_needs_details(self, context: PipelineContext) -> bool:
        """Check if scoring configuration requires Inheritance_Details."""
        if not context.scoring_config:
            return False

        # Check if any scoring formulas reference inheritance variables
        if "formulas" in context.scoring_config:
            for formula_dict in context.scoring_config["formulas"]:
                for formula_name, formula_expr in formula_dict.items():
                    if isinstance(formula_expr, str):
                        # Check if formula references common inheritance variables
                        if any(
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

        chunk_size = context.config.get("chunks", 10000)

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
                return open(f, "r")

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
        max_workers = context.config.get("threads", 1)
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
                logger.debug(f"Chunk {i+1} size: {len(chunk) if chunk is not None else 'None'}")

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
    ) -> List[pd.DataFrame]:
        """Process chunks sequentially (original implementation)."""
        output_chunks = []
        chunk_num = 0

        for chunk_df in self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size):
            logger.debug(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")
            processed_chunk = self._process_single_chunk(chunk_df, context, chunk_num)
            output_chunks.append(processed_chunk)
            chunk_num += 1

        return output_chunks

    def _process_chunks_parallel(
        self,
        sorted_file: Path,
        gene_col_name: str,
        chunk_size: int,
        context: PipelineContext,
        max_workers: int,
    ) -> List[pd.DataFrame]:
        """Process chunks in parallel using ThreadPoolExecutor."""
        # Pre-load all chunks into memory for parallel processing
        chunks = list(self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size))
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
        
        Args:
            input_file: Path to the input TSV file
            gene_col: Name of the gene column to sort by
            context: Pipeline context containing configuration
            
        Returns:
            Path to the sorted TSV file (may be the same as input if already sorted)
            
        Note:
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
                return open(f, "r")

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
                    zcat_proc.stdout.close()
                else:
                    # For uncompressed files: tail -n +2 input | sort ...
                    tail_proc = subprocess.Popen(
                        ["tail", "-n", "+2", str(input_file)], stdout=subprocess.PIPE
                    )

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

    def _read_tsv_in_gene_chunks(self, filepath: Path, gene_column: str, chunksize: int):
        """Read TSV file in gene-aware chunks."""
        import pandas as pd

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
        max_buffer_size = chunksize * 20  # Increased buffer size for very large genes

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
                if len(gene_buffer) > chunksize * 50:  # Emergency threshold
                    logger.error(
                        f"Gene buffer exceeded emergency threshold ({chunksize * 50} rows). "
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
            from ..annotator import annotate_dataframe_with_features, load_custom_features

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
        
        Args:
            chunk_df: DataFrame containing variant data for a chunk
            context: Pipeline context containing configuration and samples
            
        Returns:
            DataFrame with inheritance analysis results and cleaned columns
            
        Note:
            - Creates individual sample columns from GT column if needed
            - Applies inheritance analysis using the configured mode
            - Removes temporary sample columns after analysis to prevent huge files
            - Handles memory safety checks for large datasets
        """
        try:
            import re
            from ..inheritance.analyzer import analyze_inheritance, process_inheritance_output

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

            # Memory safety check for chunk - use more generous limit for chunks
            estimated_memory_cells = len(chunk_df) * len(vcf_samples)
            max_chunk_memory_cells = 100_000_000  # 100M cells per chunk (increased from 10M)

            if estimated_memory_cells > max_chunk_memory_cells:
                logger.warning(
                    f"Chunk too large for inheritance analysis: "
                    f"{len(chunk_df)} variants × {len(vcf_samples)} samples = "
                    f"{estimated_memory_cells:,} cells. "
                    f"Skipping inheritance analysis for this chunk."
                )
                # Add empty inheritance columns to maintain schema
                chunk_df["Inheritance_Pattern"] = "skipped"
                if preserve_details_for_scoring:
                    chunk_df["Inheritance_Details"] = "{}"
                return chunk_df

            logger.debug(
                f"Applying inheritance analysis to chunk: {len(chunk_df)} variants, "
                f"{len(vcf_samples)} samples"
            )

            # Track timing for sample column preparation
            import time
            
            # CRITICAL: Create individual sample columns from GT column
            # This is required for inheritance analysis to work correctly
            if "GT" in chunk_df.columns and len(chunk_df) > 0:
                logger.debug("Creating sample columns from GT column for inheritance analysis")
                
                prep_start = time.time()
                
                # Use the unified sample column creation function
                chunk_df = create_sample_columns_from_gt(
                    df=chunk_df,
                    vcf_samples=vcf_samples,
                    separator=context.config.get("separator", ";"),
                    snpsift_sep=context.config.get("extract_fields_separator", ",")
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
            min_variants_for_parallel = context.config.get("min_variants_for_parallel_inheritance", 100)
            
            if threads > 1 and len(chunk_df) >= min_variants_for_parallel:
                # Use parallel analyzer for better performance
                logger.debug(f"Using parallel inheritance analyzer for chunk with {threads} workers")
                
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
            if analysis_time > 0.5:  # Only log if significant
                logger.debug(f"Inheritance analysis took {analysis_time:.2f}s")

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
                chunk_df, vcf_samples, 
                preserve_columns=["GT", "Inheritance_Pattern", "Inheritance_Details"]
            )
            
            cleanup_time = time.time() - cleanup_start
            if cleanup_time > 0.1:  # Only log if significant
                logger.debug(f"Column cleanup took {cleanup_time:.2f}s")

            return chunk_df

        except Exception as e:
            return handle_inheritance_analysis_error(
                chunk_df, 
                e, 
                preserve_details_for_scoring=inheritance_config.get("preserve_details_for_scoring", False),
                context_description="chunk inheritance analysis"
            )

    def get_input_files(self, context: PipelineContext) -> List[Path]:
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

    def get_output_files(self, context: PipelineContext) -> List[Path]:
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
    def dependencies(self) -> Set[str]:
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
        gene_groups = df.groupby(gene_column)

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
                        logger.error(f"Failed to process gene {gene_name}: {str(e)}")
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
        pedigree_data: Optional[Dict],
        scoring_config: Optional[Dict],
        analysis_cfg: Dict,
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
            logger.error(f"Error processing gene {gene_name}: {str(e)}")
            raise
