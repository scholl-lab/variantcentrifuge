"""
VCF genotype extractor for inheritance analysis.

This module extracts genotype information directly from VCF files
to enable inheritance analysis without relying on TSV extraction.
"""

import subprocess
import logging
from typing import Dict, List, Tuple
import tempfile
import os

logger = logging.getLogger(__name__)


def extract_genotypes_from_vcf(
    vcf_file: str, variant_positions: List[Tuple[str, str, str, str]], sample_list: List[str]
) -> Dict[str, Dict[str, str]]:
    """
    Extract genotypes for specific variants from a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file
    variant_positions : List[Tuple[str, str, str, str]]
        List of (chrom, pos, ref, alt) tuples
    sample_list : List[str]
        List of sample IDs

    Returns
    -------
    Dict[str, Dict[str, str]]
        Dictionary mapping variant_key -> {sample_id: genotype}
    """
    if not variant_positions:
        return {}

    genotype_data = {}

    # Create a temporary regions file for bcftools
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
        for chrom, pos, ref, alt in variant_positions:
            tmp.write(f"{chrom}\t{pos}\n")
        regions_file = tmp.name

    try:
        # Use bcftools query to extract genotypes
        cmd = [
            "bcftools",
            "query",
            "-R",
            regions_file,
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
            vcf_file,
        ]

        logger.debug(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Parse the output
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                continue

            chrom, pos, ref, alt = parts[:4]
            genotypes = parts[4:] if len(parts) > 4 else []

            # Create variant key
            variant_key = f"{chrom}:{pos}:{ref}>{alt}"

            # Map genotypes to samples
            sample_genotypes = {}
            for i, gt in enumerate(genotypes):
                if i < len(sample_list):
                    sample_genotypes[sample_list[i]] = gt

            genotype_data[variant_key] = sample_genotypes

    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to extract genotypes from VCF: {e}")
        logger.error(f"stderr: {e.stderr}")
    finally:
        # Clean up temporary file
        if os.path.exists(regions_file):
            os.unlink(regions_file)

    return genotype_data


def create_variant_genotype_map(
    df, vcf_file: str, sample_list: List[str]
) -> Dict[int, Dict[str, str]]:
    """
    Create a mapping of DataFrame index to sample genotypes.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with variant data
    vcf_file : str
        Path to the VCF file
    sample_list : List[str]
        List of sample IDs

    Returns
    -------
    Dict[int, Dict[str, str]]
        Dictionary mapping DataFrame index -> {sample_id: genotype}
    """
    # Extract variant positions from DataFrame
    variant_positions = []
    index_to_variant = {}

    for idx, row in df.iterrows():
        chrom = str(row.get("CHROM", ""))
        pos = str(row.get("POS", ""))
        ref = str(row.get("REF", ""))
        alt = str(row.get("ALT", ""))

        if chrom and pos and ref and alt:
            variant_positions.append((chrom, pos, ref, alt))
            variant_key = f"{chrom}:{pos}:{ref}>{alt}"
            index_to_variant[idx] = variant_key

    # Extract genotypes from VCF
    genotype_data = extract_genotypes_from_vcf(vcf_file, variant_positions, sample_list)

    # Map back to DataFrame indices
    index_genotype_map = {}

    # Include all DataFrame indices, even if no variant data
    for idx in df.index:
        if idx in index_to_variant:
            variant_key = index_to_variant[idx]
            if variant_key in genotype_data:
                index_genotype_map[idx] = genotype_data[variant_key]
            else:
                index_genotype_map[idx] = {}
        else:
            # No variant data for this row
            index_genotype_map[idx] = {}

    return index_genotype_map


def augment_dataframe_with_vcf_genotypes(df, vcf_file: str, sample_list: List[str]):
    """
    Augment DataFrame with genotype columns from VCF.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with variant data
    vcf_file : str
        Path to the VCF file
    sample_list : List[str]
        List of sample IDs

    Returns
    -------
    pandas.DataFrame
        DataFrame with added sample genotype columns
    """
    import pandas as pd
    
    # Get genotype mapping
    genotype_map = create_variant_genotype_map(df, vcf_file, sample_list)

    # Build all genotype columns at once to avoid fragmentation
    genotype_columns = {}
    for sample_id in sample_list:
        genotype_columns[sample_id] = df.index.map(
            lambda idx: genotype_map.get(idx, {}).get(sample_id, "./.")
        ).tolist()
    
    # Create a new DataFrame with genotype columns and concatenate
    genotype_df = pd.DataFrame(genotype_columns, index=df.index)
    
    # Concatenate original DataFrame with genotype columns
    result_df = pd.concat([df, genotype_df], axis=1)
    
    return result_df
