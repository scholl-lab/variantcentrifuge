"""Helper functions for creating and cleaning sample genotype columns."""

from __future__ import annotations

import logging
import re
import time

import pandas as pd

logger = logging.getLogger(__name__)


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

                # Process each sample using vectorized assignment
                for sample_id in vcf_samples:
                    # Find entries for this sample
                    sample_mask = valid_extractions["sample_name"] == sample_id
                    sample_entries = valid_extractions[sample_mask]

                    if len(sample_entries) > 0:
                        valid_rows = sample_entries.dropna(subset=["original_index"])
                        df_copy.loc[valid_rows["original_index"].values, sample_id] = (  # type: ignore[index]
                            valid_rows["genotype"].values
                        )

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
