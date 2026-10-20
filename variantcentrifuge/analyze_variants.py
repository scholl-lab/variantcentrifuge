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

from . import gene_burden, stats
from .helpers import (
    assign_case_control_counts,
    determine_case_control_sets,
)

logger = logging.getLogger("variantcentrifuge")


__all__ = ["analyze_variants", "analyze_variants_df"]


def analyze_variants_df(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """
    Analyze variants DataFrame and compute allele frequencies, counts, and statistics.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with variant data.
    cfg : dict[str, Any]
        Configuration dictionary with sample lists, case/control assignments, stats settings.

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with variant counts and statistics columns.
    """
    if len(df) == 0:
        logger.debug("Empty DataFrame in analyze_variants_df.")
        return df.copy()

    required_columns = ["CHROM", "POS", "REF", "ALT", "GENE", "GT"]
    missing_columns = [c for c in required_columns if c not in df.columns]
    if missing_columns:
        logger.error(
            "Missing required columns: %s. Returning unchanged DataFrame.",
            ", ".join(missing_columns),
        )
        return df

    if "sample_list" not in cfg or not cfg["sample_list"].strip():
        logger.error("No sample_list found in cfg. Unable to determine the full sample set.")
        return df

    df = df.copy(deep=False)
    all_samples = set(cfg["sample_list"].split(","))
    logger.debug("Total samples from VCF header: %d", len(all_samples))

    case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)
    logger.debug(
        "Number of case samples: %d; Number of control samples: %d",
        len(case_samples),
        len(control_samples),
    )

    df = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    logger.debug("Computing basic statistics...")
    basic_stats_df = stats.compute_basic_stats(df, all_samples)
    logger.debug("Basic statistics computed.")

    stats_output_file = cfg.get("stats_output_file")
    no_stats = cfg.get("no_stats", False)

    if stats_output_file:
        logger.debug("Writing basic stats to %s", stats_output_file)
        use_compression = "gzip" if str(stats_output_file).endswith(".gz") else None
        if use_compression:
            basic_stats_df.to_csv(
                stats_output_file, sep="\t", index=False, header=True, mode="w", compression="gzip"
            )
        else:
            basic_stats_df.to_csv(stats_output_file, sep="\t", index=False, header=True, mode="w")

    if not no_stats:
        logger.debug("Computing comprehensive gene-level statistics...")
        comprehensive_df = stats.compute_gene_stats(df)
        impact_summary = stats.compute_impact_summary(df)
        variant_type_summary = stats.compute_variant_type_summary(df)

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

            if custom_stats:
                if "dataset" in custom_stats and not custom_stats["dataset"].empty:
                    logger.debug("Appending custom dataset stats to the stats file.")
                    custom_stats["dataset"].to_csv(
                        stats_output_file, sep="\t", index=False, header=False, mode="a"
                    )

                if "groups" in custom_stats:
                    for stat_name, stat_df in custom_stats["groups"].items():
                        logger.debug(
                            f"Appending custom grouped stats '{stat_name}' to the stats file."
                        )
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

    return df


def analyze_variants(lines: Iterator[str], cfg: dict[str, Any]) -> Iterator[str]:
    """
    Analyze variants and optionally perform gene burden analysis.

    Steps
    -----
    1. Parse input TSV into a DataFrame.
    2. Compute counts and stats via analyze_variants_df.
    3. If requested, perform gene burden (Fisher's exact test + correction).

    Parameters
    ----------
    lines : Iterator[str]
        Input lines representing a TSV with variant data.
    cfg : Dict[str, Any]
        Configuration dictionary.

    Yields
    ------
    str
        Processed lines of output TSV or gene-level burden results.
    """
    perform_gene_burden_flag = cfg.get("perform_gene_burden", False)

    # Read all input lines into a text block
    text_data = "".join(line for line in lines)
    if not text_data.strip():
        logger.debug("No input data provided to analyze_variants.")
        return

    df = pd.read_csv(io.StringIO(text_data), sep="\t", dtype=str)
    if len(df) == 0:
        logger.debug("Empty DataFrame after reading input in analyze_variants.")
        header_line = text_data.strip().split("\n")[0] if text_data.strip() else ""
        if header_line:
            logger.debug("Returning only header row from analyze_variants.")
            yield header_line
        return

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

    df = analyze_variants_df(df, cfg)

    # Apply inheritance analysis if requested (for public analyze_variants callers)
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
        sample_list = cfg.get("sample_list", "").split(",") if cfg.get("sample_list") else []
        logger.debug(f"Sample list for inheritance: {sample_list}")
        use_vectorized = cfg.get("use_vectorized_comp_het", True)
        from .inheritance import analyze_inheritance

        df = analyze_inheritance(
            df,
            cfg["pedigree_data"],
            sample_list,
            use_vectorized_comp_het=use_vectorized,
            reference=cfg.get("reference"),
        )
        logger.info("Inheritance analysis complete.")
    elif inheritance_already_done:
        logger.info("Inheritance analysis already performed, skipping re-calculation.")

    # Apply scoring if configuration is provided (can use inheritance results)
    scoring_config = cfg.get("scoring_config")
    if scoring_config:
        logger.info("Applying custom scoring model to variants.")
        from . import scoring

        df = scoring.apply_scoring(df, scoring_config)
        logger.debug("Scoring applied successfully.")

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
