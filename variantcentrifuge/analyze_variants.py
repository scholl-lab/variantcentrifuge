# File: variantcentrifuge/analyze_variants.py
# Location: variantcentrifuge/variantcentrifuge/analyze_variants.py

"""
Variant analysis module for gene burden and other statistics.

Modular design:
- Basic and comprehensive stats computations moved to stats.py.
- Gene burden analysis moved to gene_burden.py.
- Helper functions moved to helpers.py.

This module now:
- Reads a TSV of variants and their annotations (including a GT column listing sample genotypes).
- Classifies samples into case/control sets based on user input (sample lists or phenotype terms).
- If no case/control criteria are provided, defaults to making all samples controls.
- Invokes external modules (stats.py, gene_burden.py) for computations.
- Orchestrates flow and writes or yields output lines.

Maintains previous functionality, CLI interface, and output format.
"""

import io
import logging
from collections.abc import Iterator
from typing import Any, cast

import pandas as pd

from . import gene_burden, scoring, stats
from .helpers import (
    assign_case_control_counts,
    build_sample_phenotype_map,
    determine_case_control_sets,
    extract_sample_and_genotype,
    genotype_to_allele_count,
)
from .inheritance.analyzer import analyze_inheritance

logger = logging.getLogger("variantcentrifuge")


def analyze_variants(lines: Iterator[str], cfg: dict[str, Any]) -> Iterator[str]:
    """
    Analyze variants and optionally perform gene burden analysis.

    Steps
    -----
    1. Parse input TSV into a DataFrame.
    2. Retrieve the full sample list from `cfg["sample_list"]`.
    3. Determine case/control sets based on `cfg` (samples or phenotypes).
    4. Compute per-variant case/control allele counts.
    5. Compute basic and optionally comprehensive gene-level stats.
    6. If requested, perform gene burden (Fisher's exact test + correction).

    Parameters
    ----------
    lines : Iterator[str]
        Input lines representing a TSV with variant data.
    cfg : Dict[str, Any]
        Configuration dictionary. Keys include:
        - sample_list (str): comma-separated full sample list from VCF
        - case_samples, control_samples (List[str]): optional lists of samples
        - case_phenotypes, control_phenotypes (List[str]): optional lists of phenotypes
        - perform_gene_burden (bool): Whether to perform gene burden analysis
        - gene_burden_mode (str): "samples" or "alleles"
        - correction_method (str): "fdr" or "bonferroni"
        - no_stats (bool): Skip stats if True
        - stats_output_file (str): Path to stats output
        - gene_burden_output_file (str, optional): Path to gene burden output
        - xlsx (bool): If True, might append to Excel after analysis

    Yields
    ------
    str
        Processed lines of output TSV or gene-level burden results.
    """
    perform_gene_burden_flag = cfg.get("perform_gene_burden", False)
    stats_output_file = cfg.get("stats_output_file")
    no_stats = cfg.get("no_stats", False)

    logger.debug(
        "analyze_variants: perform_gene_burden=%s, stats_output_file=%s, no_stats=%s",
        perform_gene_burden_flag,
        stats_output_file,
        no_stats,
    )

    # Read all input lines into a text block
    text_data = "".join(line for line in lines)
    if not text_data.strip():
        logger.debug("No input data provided to analyze_variants.")
        # MODIFIED: Start of empty report generation - Return empty buffer
        # This ensures that even if no input text data is provided,
        # we still return something (for cases where empty input is acceptable)
        return
        # MODIFIED: End of empty report generation

    df = pd.read_csv(io.StringIO(text_data), sep="\t", dtype=str)
    if len(df) == 0:
        logger.debug("Empty DataFrame after reading input in analyze_variants.")
        # MODIFIED: Start of empty report generation
        # Extract header from the input text and return it
        header_line = text_data.strip().split("\n")[0] if text_data.strip() else ""
        if header_line:
            logger.debug("Returning only header row from analyze_variants.")
            yield header_line
        return
        # MODIFIED: End of empty report generation

    required_columns = ["CHROM", "POS", "REF", "ALT", "GENE", "GT"]
    missing_columns = [c for c in required_columns if c not in df.columns]
    if missing_columns:
        logger.error(
            "Missing required columns: %s. Returning unchanged lines.",
            ", ".join(missing_columns),
        )
        for line in text_data.strip().split("\n"):
            yield line
        return

    if "sample_list" not in cfg or not cfg["sample_list"].strip():
        logger.error("No sample_list found in cfg. Unable to determine the full sample set.")
        for line in text_data.strip().split("\n"):
            yield line
        return

    # all_samples: full set of samples
    all_samples = set(cfg["sample_list"].split(","))
    logger.debug("Total samples from VCF header: %d", len(all_samples))

    # Determine case/control sets
    case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)
    logger.debug(
        "Number of case samples: %d; Number of control samples: %d",
        len(case_samples),
        len(control_samples),
    )

    # Assign case/control counts per variant
    df = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    # Apply inheritance analysis first, as scores might depend on it
    logger.debug(
        f"calculate_inheritance: {cfg.get('calculate_inheritance')}, "
        f"pedigree_data type: {type(cfg.get('pedigree_data'))}"
    )
    # Check if inheritance analysis has already been performed (columns exist)
    inheritance_already_done = (
        "Inheritance_Pattern" in df.columns and "Inheritance_Details" in df.columns
    )

    if (
        cfg.get("calculate_inheritance")
        and cfg.get("pedigree_data") is not None
        and not inheritance_already_done
        and not cfg.get("skip_inheritance", False)
    ):
        logger.info("Performing inheritance pattern analysis...")
        try:
            sample_list = cfg.get("sample_list", "").split(",") if cfg.get("sample_list") else []
            logger.debug(f"Sample list for inheritance: {sample_list}")
            # Use vectorized compound het implementation by default for better performance
            use_vectorized = cfg.get("use_vectorized_comp_het", True)
            df = analyze_inheritance(
                df, cfg["pedigree_data"], sample_list, use_vectorized_comp_het=use_vectorized
            )
            logger.info("Inheritance analysis complete.")
        except Exception as e:
            logger.error(f"Inheritance analysis failed: {e}", exc_info=True)
            # Add empty columns to prevent downstream errors
            df["Inheritance_Pattern"] = "error"
            df["Inheritance_Details"] = "{}"
            # Clean up any temporary columns that might have been created
            for col in ["_inheritance_patterns", "_comp_het_info"]:
                if col in df.columns:
                    df = df.drop(columns=[col])
    elif inheritance_already_done:
        logger.info("Inheritance analysis already performed, skipping re-calculation.")

    # Apply scoring if configuration is provided (can now use inheritance results)
    scoring_config = cfg.get("scoring_config")
    if scoring_config:
        logger.info("Applying custom scoring model to variants.")
        try:
            df = scoring.apply_scoring(df, scoring_config)
            logger.debug(
                "Scoring applied successfully. New columns: %s",
                [col for col in df.columns if col not in required_columns],
            )
        except Exception as e:
            logger.error(f"Failed to apply scoring: {e}")
            # Decide whether to exit or continue without scores. Continuing is more robust.

    # Compute basic stats
    logger.debug("Computing basic statistics...")
    basic_stats_df = stats.compute_basic_stats(df, all_samples)
    logger.debug("Basic statistics computed.")

    # Write basic stats if requested
    if stats_output_file:
        logger.debug("Writing basic stats to %s", stats_output_file)

        # Apply compression to statistics output for space efficiency
        use_compression = "gzip" if str(stats_output_file).endswith(".gz") else None
        if use_compression:
            basic_stats_df.to_csv(
                stats_output_file, sep="\t", index=False, header=True, mode="w", compression="gzip"
            )
        else:
            basic_stats_df.to_csv(stats_output_file, sep="\t", index=False, header=True, mode="w")

    # Compute comprehensive stats if not skipped
    if not no_stats:
        logger.debug("Computing comprehensive gene-level statistics...")
        comprehensive_df = stats.compute_gene_stats(df)
        impact_summary = stats.compute_impact_summary(df)
        variant_type_summary = stats.compute_variant_type_summary(df)

        # Compute custom statistics if configured
        custom_stats = None
        stats_config_path = cfg.get("stats_config")
        if stats_config_path:
            logger.info(f"Computing custom statistics using config: {stats_config_path}")
            custom_stats = stats.compute_custom_stats(df, stats_config_path)

        combined_stats = stats.merge_and_format_stats(
            comprehensive_df, impact_summary, variant_type_summary, custom_stats
        )
        logger.debug("Comprehensive statistics computed.")

        comp_stats_list = []
        for _, row in combined_stats.iterrows():
            gene = row["GENE"]
            for col in combined_stats.columns:
                if col == "GENE":
                    continue
                val = row[col]
                comp_stats_list.append([f"{gene}_{col}", str(val)])
        if comp_stats_list and stats_output_file:
            logger.debug("Appending comprehensive stats to the stats file.")
            comp_stats_df = pd.DataFrame(comp_stats_list, columns=["metric", "value"])
            # Apply compression consistently with the main stats file
            use_compression = "gzip" if str(stats_output_file).endswith(".gz") else None
            if use_compression:
                comp_stats_df.to_csv(
                    stats_output_file,
                    sep="\t",
                    index=False,
                    header=False,
                    mode="a",
                    compression="gzip",
                )
            else:
                comp_stats_df.to_csv(
                    stats_output_file, sep="\t", index=False, header=False, mode="a"
                )

            # Write custom stats if available
            if custom_stats:
                # Dataset-level custom stats
                if "dataset" in custom_stats and not custom_stats["dataset"].empty:
                    logger.debug("Appending custom dataset stats to the stats file.")
                    custom_stats["dataset"].to_csv(
                        stats_output_file, sep="\t", index=False, header=False, mode="a"
                    )

                # Grouped custom stats (formatted as metric/value pairs)
                if "groups" in custom_stats:
                    for stat_name, stat_df in custom_stats["groups"].items():
                        logger.debug(
                            f"Appending custom grouped stats '{stat_name}' to the stats file."
                        )
                        # Convert to metric/value format
                        for idx, row in stat_df.iterrows():
                            if isinstance(row, pd.Series):
                                for col in stat_df.columns:
                                    if col not in stat_df.index.names:
                                        metric_name = f"{stat_name}_{idx}_{col}"
                                        comp_stats_list.append([metric_name, str(row[col])])

                    if comp_stats_list:
                        additional_stats_df = pd.DataFrame(
                            comp_stats_list, columns=["metric", "value"]
                        )
                        additional_stats_df.to_csv(
                            stats_output_file, sep="\t", index=False, header=False, mode="a"
                        )
    else:
        logger.debug("No comprehensive stats requested (no_stats=True).")

    # Perform gene burden analysis if requested
    if perform_gene_burden_flag:
        logger.debug("Performing configurable gene burden analysis...")
        gene_burden_results = gene_burden.perform_gene_burden_analysis(df, cfg)

        if gene_burden_results.empty:
            logger.warning("Gene burden results are empty. Returning variant-level data.")
            out_str = df.to_csv(sep="\t", index=False)
            for line in out_str.strip().split("\n"):
                yield line
            return

        burden_output_file = cfg.get("gene_burden_output_file", "gene_burden_results.tsv")
        logger.info("Writing gene burden results to %s", burden_output_file)
        # Apply compression to gene burden output
        use_compression = "gzip" if str(burden_output_file).endswith(".gz") else None
        gene_burden_results.to_csv(
            burden_output_file, sep="\t", index=False, compression=cast(Any, use_compression)
        )

        # Yield burden results as output
        burden_text = gene_burden_results.to_csv(sep="\t", index=False)
        for line in burden_text.strip().split("\n"):
            yield line
    else:
        logger.debug("Gene burden not requested, returning processed data.")
        out_str = df.to_csv(sep="\t", index=False)
        for line in out_str.strip().split("\n"):
            yield line
