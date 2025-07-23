"""
Vectorized genotype replacement module using pandas.

This module provides a high-performance alternative to the streaming genotype replacement
in replacer.py by using pandas vectorized string operations for processing TSV files.

Key optimizations:
- Uses pandas.str accessor for vectorized string operations
- Leverages DataFrame operations instead of row-by-row processing
- Implements regex-based genotype transformations with Series.str.replace
- Handles missing allele normalization with vectorized operations
- Memory-efficient chunked processing for large datasets

Performance considerations:
- Trades memory usage for processing speed
- Best suited for datasets that fit in memory or can be processed in chunks
- Significantly faster than line-by-line processing for large sample counts
- Scales well with number of samples due to vectorization
"""

import logging
import re
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Union

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def _normalize_snpeff_column_name(col_name: str) -> str:
    """
    Normalize known SnpEff or SnpSift prefixes from a single column name.

    For example, "GEN[*].DP" or "ANN[*].Effect" might become "DP" or "Effect".
    """
    return (
        col_name.replace("ANN[*].", "")
        .replace("ANN[0].", "")
        .replace("GEN[*].", "")
        .replace("GEN[0].", "")
        .replace("NMD[*].", "NMD_")
        .replace("NMD[0].", "NMD_")
        .replace("LOF[*].", "LOF_")
        .replace("LOF[0].", "LOF_")
    )


class VectorizedGenotypeReplacer:
    """
    High-performance vectorized genotype replacement using pandas operations.

    This class provides a pandas-based alternative to the streaming genotype replacement
    approach, optimizing for speed through vectorized string operations.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the vectorized genotype replacer.

        Parameters
        ----------
        config : Dict[str, Any]
            Configuration dictionary containing replacement parameters
        """
        self.config = config
        self.separator = config.get("separator", ";")
        self.genotype_replacement_map = config.get("genotype_replacement_map", {r"[2-9]": "1"})
        self.snpsift_sep = config.get("extract_fields_separator", ":")

        # Extra fields configuration
        self.append_extra_fields = config.get("append_extra_sample_fields", False)
        self.extra_sample_fields = config.get("extra_sample_fields", [])
        self.extra_field_delimiter = config.get("extra_sample_field_delimiter", ":")

        # Parse and validate samples
        sample_list_str = config.get("sample_list", "")
        self.samples = [s.strip() for s in sample_list_str.split(",") if s.strip()]

        if not self.samples:
            raise ValueError("No samples found in config['sample_list']")

        # Compile regex patterns for performance
        self.compiled_patterns = []
        for pattern, replacement in self.genotype_replacement_map.items():
            compiled = re.compile(pattern)
            self.compiled_patterns.append((compiled, replacement))

        logger.debug(f"Initialized vectorized replacer for {len(self.samples)} samples")

    def process_file(self, input_path: Union[str, Path], output_path: Union[str, Path]) -> None:
        """
        Process an entire TSV file with vectorized genotype replacement.

        Parameters
        ----------
        input_path : str or Path
            Path to input TSV file
        output_path : str or Path
            Path to output TSV file
        """
        input_path = Path(input_path)
        output_path = Path(output_path)

        logger.info(f"Processing genotype replacement: {input_path} -> {output_path}")

        # Read the TSV file
        compression = "gzip" if str(input_path).endswith(".gz") else None
        df = pd.read_csv(input_path, sep="\t", dtype=str, compression=compression)

        # Process genotypes
        df = self._process_dataframe(df)

        # Write output
        output_compression = "gzip" if str(output_path).endswith(".gz") else None
        df.to_csv(output_path, sep="\t", index=False, compression=output_compression)

        logger.info(f"Completed vectorized genotype replacement for {len(df)} variants")

    def process_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply vectorized genotype replacement to a DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame with GT column

        Returns
        -------
        pd.DataFrame
            DataFrame with processed GT column
        """
        return self._process_dataframe(df.copy())

    def _process_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Process DataFrame with vectorized operations.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame to process

        Returns
        -------
        pd.DataFrame
            Processed DataFrame
        """
        if "GT" not in df.columns:
            logger.warning("No GT column found in DataFrame. Returning unchanged.")
            return df

        # Create column mapping for extra fields if needed
        extra_field_indices = {}
        if self.append_extra_fields and self.extra_sample_fields:
            # Build normalized column mapping
            normalized_to_col = {}
            for col in df.columns:
                norm_col = _normalize_snpeff_column_name(col)
                normalized_to_col[norm_col] = col

            # Map requested extra fields to actual columns
            for raw_field in self.extra_sample_fields:
                norm_field = _normalize_snpeff_column_name(raw_field)
                if norm_field in normalized_to_col:
                    actual_col = normalized_to_col[norm_field]
                    extra_field_indices[raw_field] = actual_col
                    logger.debug(f"Mapped extra field '{raw_field}' to column '{actual_col}'")
                else:
                    logger.warning(f"Extra field '{raw_field}' not found in DataFrame")

        # Process GT column with vectorized operations
        logger.debug("Starting vectorized genotype processing")

        # Split GT column into individual genotype series per sample
        gt_split = df["GT"].str.split(self.snpsift_sep, expand=True)

        # Ensure we have enough columns for all samples
        num_samples = len(self.samples)
        if gt_split.shape[1] < num_samples:
            # Pad with NaN columns if needed
            missing_cols = num_samples - gt_split.shape[1]
            for i in range(missing_cols):
                gt_split[gt_split.shape[1] + i] = np.nan

        # Process each sample's genotypes
        processed_genotypes = []

        for sample_idx, sample_name in enumerate(self.samples):
            if sample_idx >= gt_split.shape[1]:
                logger.warning(f"Sample {sample_name} has no corresponding genotype column")
                continue

            # Get genotype series for this sample
            genotype_series = gt_split[sample_idx].fillna("")

            # Apply vectorized genotype transformations
            processed_series = self._process_genotype_series(genotype_series)

            # Filter out non-variant genotypes (0/0, ./.)
            variant_mask = ~processed_series.isin(["0/0", "./.", ""])

            # Build sample output strings
            sample_strings = self._build_sample_strings(
                processed_series, sample_name, sample_idx, df, extra_field_indices, variant_mask
            )

            processed_genotypes.append(sample_strings)

        # Combine all samples into final GT column
        # Stack all sample strings and join with separator
        final_gt = self._combine_sample_genotypes(processed_genotypes)
        df["GT"] = final_gt

        logger.debug("Completed vectorized genotype processing")
        return df

    def _process_genotype_series(self, genotype_series: pd.Series) -> pd.Series:
        """
        Apply vectorized genotype transformations to a pandas Series.

        Parameters
        ----------
        genotype_series : pd.Series
            Series of genotype strings for one sample

        Returns
        -------
        pd.Series
            Processed genotype series
        """
        # Step 1: Replace phased genotypes (| -> /)
        processed = genotype_series.str.replace(r"\|", "/", regex=True)

        # Step 2: Normalize missing alleles (./N -> 0/N, N/. -> N/0)
        # Handle ./N case
        processed = processed.str.replace(r"^\.\/(\d+)$", r"0/\1", regex=True)
        # Handle N/. case
        processed = processed.str.replace(r"^(\d+)\/\.$", r"\1/0", regex=True)

        # Step 3: Apply configured genotype replacements
        for compiled_pattern, replacement in self.compiled_patterns:
            # Apply regex replacement to each genotype
            processed = processed.str.replace(compiled_pattern, replacement, regex=True)

        # Step 4: Normalize multi-digit alleles to binary (any non-0 digit -> 1)
        # Split diploid genotypes and process each allele
        diploid_mask = processed.str.contains("/", na=False)

        # Process diploid genotypes
        if diploid_mask.any():
            diploid_genotypes = processed[diploid_mask]

            # Split into alleles
            allele_split = diploid_genotypes.str.split("/", expand=True)
            if allele_split.shape[1] >= 2:
                # Process each allele: convert any digit > 0 to 1
                allele1 = allele_split[0].str.replace(r"^[1-9]\d*$", "1", regex=True)
                allele2 = allele_split[1].str.replace(r"^[1-9]\d*$", "1", regex=True)

                # Rejoin alleles
                processed_diploid = allele1 + "/" + allele2
                processed.loc[diploid_mask] = processed_diploid

        # Process haploid genotypes (single digits)
        haploid_mask = processed.str.match(r"^\d+$", na=False)
        if haploid_mask.any():
            processed.loc[haploid_mask] = processed.loc[haploid_mask].str.replace(
                r"^[1-9]\d*$", "1", regex=True
            )

        return processed

    def _build_sample_strings(
        self,
        genotype_series: pd.Series,
        sample_name: str,
        sample_idx: int,
        df: pd.DataFrame,
        extra_field_indices: Dict[str, str],
        variant_mask: pd.Series,
    ) -> pd.Series:
        """
        Build sample output strings with optional extra fields.

        Parameters
        ----------
        genotype_series : pd.Series
            Processed genotype series
        sample_name : str
            Name of the sample
        sample_idx : int
            Index of the sample
        df : pd.DataFrame
            Original DataFrame
        extra_field_indices : Dict[str, str]
            Mapping of extra fields to column names
        variant_mask : pd.Series
            Boolean mask for variant genotypes

        Returns
        -------
        pd.Series
            Series of sample strings
        """
        # Start with empty strings
        result = pd.Series("", index=genotype_series.index, dtype=str)

        if not variant_mask.any():
            return result

        # Get variant genotypes
        variant_genotypes = genotype_series[variant_mask]

        if self.append_extra_fields and extra_field_indices:
            # Build extra field strings for each variant
            extra_parts = []

            for raw_field, col_name in extra_field_indices.items():
                if col_name in df.columns:
                    # Split extra field column by separator and get sample_idx
                    extra_split = df[col_name].str.split(self.snpsift_sep, expand=True)
                    if sample_idx < extra_split.shape[1]:
                        extra_values = extra_split[sample_idx].fillna("")
                        extra_parts.append(extra_values[variant_mask])
                    else:
                        extra_parts.append(pd.Series("", index=variant_genotypes.index))
                else:
                    extra_parts.append(pd.Series("", index=variant_genotypes.index))

            if extra_parts:
                # Combine extra fields with delimiter
                extra_combined = extra_parts[0].astype(str)
                for extra_part in extra_parts[1:]:
                    extra_combined = (
                        extra_combined + self.extra_field_delimiter + extra_part.astype(str)
                    )

                # Remove any empty extra fields (replace empty strings with empty)
                extra_combined = extra_combined.str.replace(r"^:+$", "", regex=True)
                extra_combined = extra_combined.str.replace(r"^:+", "", regex=True)
                extra_combined = extra_combined.str.replace(r":+$", "", regex=True)

                # Combine genotype with extra fields (only if extra fields are not empty)
                mask_with_extra = extra_combined != ""
                genotype_with_extra = variant_genotypes.copy()
                genotype_with_extra.loc[mask_with_extra] = (
                    variant_genotypes.loc[mask_with_extra]
                    + self.extra_field_delimiter
                    + extra_combined.loc[mask_with_extra]
                )
            else:
                genotype_with_extra = variant_genotypes
        else:
            genotype_with_extra = variant_genotypes

        # Format as sample(genotype)
        sample_strings = f"{sample_name}(" + genotype_with_extra + ")"

        # Place sample strings in result at correct positions
        result.loc[variant_mask] = sample_strings

        return result

    def _combine_sample_genotypes(self, sample_genotype_lists: List[pd.Series]) -> pd.Series:
        """
        Combine sample genotype strings into final GT column.

        Parameters
        ----------
        sample_genotype_lists : List[pd.Series]
            List of sample genotype series

        Returns
        -------
        pd.Series
            Combined GT series
        """
        if not sample_genotype_lists:
            return pd.Series("", index=range(0))

        # Get index from first series
        index = sample_genotype_lists[0].index

        # Combine all non-empty strings with separator
        result = pd.Series("", index=index, dtype=str)

        for row_idx in index:
            row_parts = []
            for sample_series in sample_genotype_lists:
                if row_idx in sample_series.index:
                    value = sample_series.loc[row_idx]
                    if value and value != "":
                        row_parts.append(value)

            if row_parts:
                result.loc[row_idx] = self.separator.join(row_parts)
            else:
                result.loc[row_idx] = ""

        return result


def replace_genotypes_vectorized(
    input_path: Union[str, Path], output_path: Union[str, Path], config: Dict[str, Any]
) -> None:
    """
    High-level function for vectorized genotype replacement.

    Parameters
    ----------
    input_path : str or Path
        Path to input TSV file
    output_path : str or Path
        Path to output TSV file
    config : Dict[str, Any]
        Configuration dictionary
    """
    replacer = VectorizedGenotypeReplacer(config)
    replacer.process_file(input_path, output_path)


def process_chunked_vectorized(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    config: Dict[str, Any],
    chunk_size: int = 10000,
) -> None:
    """
    Process large files using chunked vectorized replacement.

    Parameters
    ----------
    input_path : str or Path
        Path to input TSV file
    output_path : str or Path
        Path to output TSV file
    config : Dict[str, Any]
        Configuration dictionary
    chunk_size : int
        Number of rows per chunk
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    replacer = VectorizedGenotypeReplacer(config)

    # Read and process in chunks
    compression = "gzip" if str(input_path).endswith(".gz") else None
    output_compression = "gzip" if str(output_path).endswith(".gz") else None

    first_chunk = True

    for chunk_df in pd.read_csv(
        input_path, sep="\t", dtype=str, compression=compression, chunksize=chunk_size
    ):
        processed_chunk = replacer.process_dataframe(chunk_df)

        # Write chunk to output
        mode = "w" if first_chunk else "a"
        header = first_chunk

        processed_chunk.to_csv(
            output_path,
            sep="\t",
            index=False,
            compression=output_compression,
            mode=mode,
            header=header,
        )

        first_chunk = False

        logger.debug(f"Processed chunk of {len(chunk_df)} rows")

    logger.info(f"Completed chunked vectorized genotype replacement: {input_path} -> {output_path}")


def _process_chunk_worker(
    chunk_data: pd.DataFrame, config: Dict[str, Any], chunk_id: int, temp_dir: str
) -> str:
    """
    Worker function to process a single chunk in parallel.

    Parameters
    ----------
    chunk_data : pd.DataFrame
        Chunk of data to process
    config : Dict[str, Any]
        Configuration dictionary
    chunk_id : int
        Unique identifier for this chunk
    temp_dir : str
        Temporary directory for chunk outputs

    Returns
    -------
    str
        Path to the processed chunk file
    """
    try:
        replacer = VectorizedGenotypeReplacer(config)
        processed_chunk = replacer.process_dataframe(chunk_data)

        # Write chunk to temporary file
        chunk_output = Path(temp_dir) / f"chunk_{chunk_id:06d}.tsv"
        processed_chunk.to_csv(chunk_output, sep="\t", index=False)

        logger.debug(f"Worker completed chunk {chunk_id} with {len(processed_chunk)} rows")
        return str(chunk_output)

    except Exception as e:
        logger.error(f"Worker failed processing chunk {chunk_id}: {e}")
        raise


def process_parallel_chunked_vectorized(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    config: Dict[str, Any],
    chunk_size: int = 10000,
    max_workers: int = None,
    available_memory_gb: float = None,
) -> None:
    """
    Process large files using parallel chunked vectorized replacement.

    This function splits the input file into chunks and processes multiple chunks
    concurrently using separate processes to maximize CPU utilization while
    maintaining memory safety.

    Parameters
    ----------
    input_path : str or Path
        Path to input TSV file
    output_path : str or Path
        Path to output TSV file
    config : Dict[str, Any]
        Configuration dictionary
    chunk_size : int
        Number of rows per chunk
    max_workers : int, optional
        Maximum number of parallel workers. If None, calculated based on memory.
    available_memory_gb : float, optional
        Available system memory in GB. If None, will be detected.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    # Import psutil for memory detection if not provided
    if available_memory_gb is None:
        try:
            import psutil

            memory = psutil.virtual_memory()
            available_memory_gb = memory.available / (1024**3)
        except ImportError:
            logger.warning("psutil not available, using default memory estimate")
            available_memory_gb = 8.0

    # Estimate memory per chunk and calculate safe worker count
    sample_list_str = config.get("sample_list", "")
    num_samples = len([s.strip() for s in sample_list_str.split(",") if s.strip()])

    # Conservative estimate: base + (samples * bytes_per_sample) * pandas_overhead * chunk_rows
    base_bytes_per_row = 200
    bytes_per_sample = 60 if num_samples > 1000 else 50
    pandas_overhead = 4 if num_samples > 3000 else 3
    estimated_bytes_per_row = (
        base_bytes_per_row + (bytes_per_sample * num_samples)
    ) * pandas_overhead
    estimated_chunk_memory_gb = (estimated_bytes_per_row * chunk_size) / (1024**3)

    # Calculate safe number of workers with more aggressive scaling
    # Since we're streaming output (not loading chunks back), we can use more memory per worker
    safe_memory_gb = available_memory_gb * 0.8  # Use 80% of available memory
    if max_workers is None:
        # Use more threads and higher memory per worker since scaling is favorable
        memory_based_workers = int(safe_memory_gb / estimated_chunk_memory_gb)
        # Don't cap at 8 - use actual thread count provided, up to memory limit
        max_workers = max(1, min(memory_based_workers, 16))  # Cap at 16 for system stability

    logger.info(
        f"Using parallel chunked processing with {max_workers} workers "
        f"(estimated {estimated_chunk_memory_gb:.2f}GB per chunk, "
        f"safe memory limit {safe_memory_gb:.1f}GB)"
    )

    compression = "gzip" if str(input_path).endswith(".gz") else None
    output_compression = "gzip" if str(output_path).endswith(".gz") else None

    # Create temporary directory for chunk outputs
    with tempfile.TemporaryDirectory(prefix="variantcentrifuge_chunks_") as temp_dir:
        logger.debug(f"Using temporary directory: {temp_dir}")

        # Read file and split into chunks for parallel processing
        chunk_files = []
        header = None

        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit chunk processing jobs
                future_to_chunk = {}
                chunk_id = 0

                for chunk_df in pd.read_csv(
                    input_path, sep="\t", dtype=str, compression=compression, chunksize=chunk_size
                ):
                    if header is None:
                        header = chunk_df.columns.tolist()

                    # Submit chunk for processing
                    future = executor.submit(
                        _process_chunk_worker, chunk_df, config, chunk_id, temp_dir
                    )
                    future_to_chunk[future] = chunk_id
                    chunk_id += 1

                total_chunks = chunk_id
                logger.info(f"Submitted {total_chunks} chunks for parallel processing")

                # Collect results in order with progress logging
                chunk_results = {}
                completed_count = 0
                for future in as_completed(future_to_chunk):
                    chunk_id = future_to_chunk[future]
                    try:
                        chunk_file = future.result()
                        chunk_results[chunk_id] = chunk_file
                        completed_count += 1

                        # Progress logging every 10% or minimum every 5 chunks
                        progress_interval = max(5, total_chunks // 10)
                        if (
                            completed_count % progress_interval == 0
                            or completed_count == total_chunks
                        ):
                            progress_pct = (completed_count / total_chunks) * 100
                            logger.info(
                                f"Chunk progress: {completed_count}/{total_chunks} "
                                f"({progress_pct:.1f}%) completed"
                            )
                        else:
                            logger.debug(f"Completed chunk {chunk_id}")
                    except Exception as e:
                        logger.error(f"Failed to process chunk {chunk_id}: {e}")
                        raise

                # Stream combine results in correct order (memory-efficient)
                logger.info(f"Stream combining {len(chunk_results)} processed chunks")

                # Open output file for writing
                if output_compression == "gzip":
                    import gzip

                    output_file = gzip.open(output_path, "wt", encoding="utf-8")
                else:
                    output_file = open(output_path, "w", encoding="utf-8")

                try:
                    first_chunk = True
                    combined_count = 0
                    for chunk_id in sorted(chunk_results.keys()):
                        chunk_file = chunk_results[chunk_id]

                        # Stream copy chunk file to output
                        with open(chunk_file, "r", encoding="utf-8") as chunk_input:
                            if first_chunk:
                                # First chunk: copy everything including header
                                for line in chunk_input:
                                    output_file.write(line)
                                first_chunk = False
                            else:
                                # Subsequent chunks: skip header line
                                next(chunk_input)  # Skip header
                                for line in chunk_input:
                                    output_file.write(line)

                        # Clean up chunk file immediately after processing
                        Path(chunk_file).unlink()
                        combined_count += 1

                        # Progress logging for combination phase
                        combine_interval = max(5, total_chunks // 10)
                        if combined_count % combine_interval == 0 or combined_count == total_chunks:
                            combine_pct = (combined_count / total_chunks) * 100
                            logger.info(
                                f"Combination progress: {combined_count}/{total_chunks} "
                                f"({combine_pct:.1f}%) combined"
                            )
                        else:
                            logger.debug(f"Stream copied and cleaned chunk {chunk_id}")

                finally:
                    output_file.close()

        except Exception as e:
            logger.error(f"Parallel chunked processing failed: {e}")
            raise

    final_chunk_count = len(chunk_results) if "chunk_results" in locals() else 0
    logger.info(
        f"✓ Completed parallel chunked vectorized genotype replacement: "
        f"{input_path} -> {output_path} "
        f"({final_chunk_count} chunks processed with {max_workers} workers)"
    )
