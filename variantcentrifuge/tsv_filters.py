# File: variantcentrifuge/tsv_filters.py
# Location: variantcentrifuge/variantcentrifuge/tsv_filters.py

"""
TSV and DataFrame filtering utilities.

Provides genotype-based TSV row filtering (`filter_final_tsv_by_genotype`),
expression-based TSV filtering (`filter_tsv_with_expression`), and in-memory
DataFrame query filtering (`filter_dataframe_with_query`).
"""

import logging
import os
import re
from io import StringIO
from typing import Any, cast

import pandas as pd

from .utils import smart_open

logger = logging.getLogger("variantcentrifuge")


def filter_final_tsv_by_genotype(
    input_tsv: str,
    output_tsv: str,
    global_genotypes: set[str] | None = None,
    gene_genotype_file: str | None = None,
    gene_column_name: str = "GENE",
    gt_column_name: str = "GT",
) -> None:
    """
    Filter the final TSV rows by genotype.

    This can be done globally (using a single set of requested genotypes like {"het"}, {"hom"},
    or {"comp_het"}) or on a per-gene basis if gene_genotype_file is provided.
    The gene_genotype_file must contain at least two columns:
      1) GENE
      2) GENOTYPES (one or more of het, hom, comp_het, comma-separated)

    The logic is:
      - 'het'       => keep samples with genotype 0/1 or 1/0
      - 'hom'       => keep samples with genotype 1/1
      - 'comp_het'  => keep samples that have at least two distinct variants in the same gene
                       with genotype 0/1 or 1/0

    If a gene is defined in gene_genotype_file, the union of those genotype rules is applied.
    If a gene is not defined in that file (or if none is provided), global_genotypes is used.

    The resulting TSV keeps only lines that have at least one sample fulfilling the chosen
    genotype filters. If no samples remain on a line, that line is discarded.

    Additionally, if a sample passes because of 'het' or 'hom' or 'comp_het',
    we append a reason marker. For example:

    325879(0/1:53,55:108) => 325879(0/1:53,55:108)(het) or 325879(0/1:53,55:108)(het,comphet)

    Parameters
    ----------
    input_tsv : str
        Path to the input TSV, e.g. the final genotype_replaced TSV.
    output_tsv : str
        Path to the filtered output TSV.
    global_genotypes : set of str, optional
        A set of genotype filters to apply globally if no gene-specific rule is found.
        E.g. {"het"}, {"hom"}, or {"comp_het"}, or any combination thereof.
    gene_genotype_file : str, optional
        Path to a file with columns 'GENE' and 'GENOTYPES'. Each row can specify one or more
        genotype filters for a given gene, e.g. "BRCA2    het,comp_het".
    gene_column_name : str
        The column name in the TSV that specifies the gene.
    gt_column_name : str
        The column name in the TSV that specifies the genotype(s) for each sample.

    Returns
    -------
    None
        A new TSV is written to output_tsv.

    Notes
    -----
    For 'comp_het', we group rows by gene, identify samples that appear at least twice
    (with 0/1 or 1/0), then keep those rows for those samples only. We also annotate
    each genotype with the reason(s) it passed (het, hom, comphet).

    Examples
    --------
    >>> filter_final_tsv_by_genotype(
    ...     "input.genotype_replaced.tsv",
    ...     "output.genotype_filtered.tsv",
    ...     global_genotypes={"comp_het"}
    ... )
    """
    if global_genotypes is None:
        global_genotypes = set()

    logger.debug("filter_final_tsv_by_genotype called with:")
    logger.debug("  input_tsv=%s", input_tsv)
    logger.debug("  output_tsv=%s", output_tsv)
    logger.debug("  global_genotypes=%s", global_genotypes)
    logger.debug("  gene_genotype_file=%s", gene_genotype_file)

    # Show the first few lines of the input_tsv (if available)
    try:
        with open(input_tsv, encoding="utf-8") as inp_debug:
            logger.debug("First lines from input_tsv:")
            for i in range(3):
                line = inp_debug.readline()
                if not line:
                    break
                logger.debug("Line %d: %s", i + 1, line.rstrip("\n"))
    except FileNotFoundError:
        logger.error("input_tsv file not found: %s", input_tsv)
        raise

    # Read the gene -> genotype(s) mapping if provided
    gene_to_genotypes: dict[str, set[str]] = {}
    if gene_genotype_file and os.path.exists(gene_genotype_file):
        logger.debug("Attempting to read gene -> genotype rules from: %s", gene_genotype_file)
        with open(gene_genotype_file, encoding="utf-8") as gfile:
            all_lines = gfile.readlines()

        logger.debug("First lines from gene_genotype_file:")
        for i, line in enumerate(all_lines[:3]):
            logger.debug("Line %d: %s", i + 1, line.rstrip("\n"))

        if not all_lines:
            logger.debug("Gene genotype file is empty, skipping parsing.")
        else:
            gfile_replay = StringIO("".join(all_lines))
            gg_header_cols = next(gfile_replay).strip().split("\t")

            # Expect at least columns: GENE, GENOTYPES
            gg_gene_idx: int | None = None
            geno_idx: int | None = None
            for i, col in enumerate(gg_header_cols):
                if col.upper() == "GENE":
                    gg_gene_idx = i
                elif col.upper() == "GENOTYPES":
                    geno_idx = i
            if gg_gene_idx is None or geno_idx is None:
                raise ValueError(
                    "gene_genotype_file must have columns named 'GENE' and 'GENOTYPES'."
                )
            line_count = 0
            for line in gfile_replay:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) <= max(gg_gene_idx, geno_idx):
                    continue
                gname = parts[gg_gene_idx].strip()
                genos = parts[geno_idx].replace(" ", "").split(",")
                if gname:
                    if gname not in gene_to_genotypes:
                        gene_to_genotypes[gname] = set()
                    gene_to_genotypes[gname].update(genos)
                line_count += 1
            logger.debug("Finished reading %d data lines from gene_genotype_file", line_count)
    else:
        logger.debug(
            "No valid gene_genotype_file found, or file does not exist at: %s",
            gene_genotype_file,
        )

    # Helpers to detect genotype as 'het' (0/1 or 1/0) or 'hom' (1/1)
    def is_het(gt_string: str) -> bool:
        return gt_string in ("0/1", "1/0", "0|1", "1|0")

    def is_hom(gt_string: str) -> bool:
        return gt_string in ("1/1", "1|1")

    # We'll parse the lines grouped by gene, so we can do 'comp_het' logic.
    lines_by_gene: dict[str, list[tuple[str, dict[str, str]]]] = {}

    with smart_open(input_tsv, "r") as inp:
        header: str = next(inp).rstrip("\n")
        header_cols = header.split("\t")
        # Identify gene and GT columns
        try:
            gene_idx: int = header_cols.index(gene_column_name)
        except ValueError as e:
            raise ValueError(
                f"Could not find gene column '{gene_column_name}' in TSV header."
            ) from e
        try:
            gt_idx: int = header_cols.index(gt_column_name)
        except ValueError as e:
            raise ValueError(
                f"Could not find genotype column '{gt_column_name}' in TSV header."
            ) from e

        for line in inp:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) <= max(gene_idx, gt_idx):
                continue

            gene_val = parts[gene_idx].strip()
            gt_val = parts[gt_idx].strip()

            # Parse each sample( genotype ) entry in GT
            sample_genotypes = {}
            if gt_val:
                entries = gt_val.split(";")
                for entry in entries:
                    entry = entry.strip()
                    if not entry:
                        continue
                    # example: "325879(0/1:53,55:108)"
                    # We'll store the entire substring "0/1:53,55:108"
                    # but also parse out the "main genotype" before a colon if present
                    if "(" in entry and ")" in entry:
                        sample_name = entry.split("(")[0].strip()
                        inside_paren = entry.split("(")[1].strip(")")
                        sample_genotypes[sample_name] = inside_paren

            if gene_val not in lines_by_gene:
                lines_by_gene[gene_val] = []
            lines_by_gene[gene_val].append((line, sample_genotypes))

    # Identify which genes require 'comp_het'
    def gene_filters(g: str) -> set[str]:
        if g in gene_to_genotypes:
            return gene_to_genotypes[g]
        return global_genotypes

    comp_het_qualified: dict[str, set[str]] = {}
    all_genes = set(lines_by_gene.keys())
    relevant_for_comp_het = {g for g in all_genes if "comp_het" in gene_filters(g)}

    # For each gene that uses comp_het, gather samples that have >=2 het variants
    for g in relevant_for_comp_het:
        sample_het_count: dict[str, int] = {}
        for _line_str, sample_gt_dict in lines_by_gene[g]:
            for sample_name, genotype_substring in sample_gt_dict.items():
                # e.g. genotype_substring = "0/1:53,55:108"
                # parse the main genotype by splitting on the first colon
                main_gt = (
                    genotype_substring.split(":")[0]
                    if ":" in genotype_substring
                    else genotype_substring
                )
                if is_het(main_gt):
                    sample_het_count[sample_name] = sample_het_count.get(sample_name, 0) + 1
        qualified_samples = {s for s, c in sample_het_count.items() if c >= 2}
        comp_het_qualified[g] = qualified_samples
        logger.debug("Gene '%s' => comp_het_qualified samples: %s", g, qualified_samples)

    # We'll build our filtered lines
    filtered_lines: list[str] = [header]  # keep original header as is

    # Now go through each gene's lines
    for g, items in lines_by_gene.items():
        filters_for_gene = gene_filters(g)
        do_het = "het" in filters_for_gene
        do_hom = "hom" in filters_for_gene
        do_comp_het = "comp_het" in filters_for_gene

        for line_str, sample_gt_dict in items:
            parts = line_str.split("\t")
            new_sample_entries = []

            # For each sample, check if it meets the filter(s)
            for sample_name, genotype_substring in sample_gt_dict.items():
                main_gt = (
                    genotype_substring.split(":")[0]
                    if ":" in genotype_substring
                    else genotype_substring
                )

                reasons = []
                if do_het and is_het(main_gt):
                    reasons.append("het")
                if do_hom and is_hom(main_gt):
                    reasons.append("hom")
                if (
                    do_comp_het
                    and is_het(main_gt)
                    and sample_name in comp_het_qualified.get(g, set())
                ):
                    reasons.append("comphet")

                if reasons:
                    # This sample passes the filter => we append the reason(s)
                    reason_str = f"({','.join(reasons)})"
                    # e.g. "325879(0/1:53,55:108)(het,comphet)"
                    new_sample_entries.append(f"{sample_name}({genotype_substring}){reason_str}")

            # If we have at least one sample that passes => keep the line
            if new_sample_entries:
                parts[gt_idx] = "; ".join(new_sample_entries)
                filtered_lines.append("\t".join(parts))

    # Write out the filtered lines
    with smart_open(output_tsv, "w") as out:
        for line in filtered_lines:
            out.write(line + "\n")

    logger.info("Genotype filtering complete. Input: %s, Output: %s", input_tsv, output_tsv)


def filter_tsv_with_expression(
    input_tsv: str, output_tsv: str, filter_expression: str, pandas_query: bool = True
) -> None:
    """
    Filter a TSV file using a filter expression.

    This function enables filtering on any column in the TSV, including computed
    columns like scores, inheritance patterns, and custom annotations that are
    added during the analysis pipeline.

    Parameters
    ----------
    input_tsv : str
        Path to the input TSV file
    output_tsv : str
        Path to the output filtered TSV file
    filter_expression : str
        Filter expression. If pandas_query is True, this should be a pandas query
        string (e.g., "Score > 0.5 & Impact == 'HIGH'"). If False, it should be
        a SnpSift-style expression that will be translated.
    pandas_query : bool
        If True, use pandas query syntax. If False, translate from SnpSift syntax.
    """
    logger.info(f"Applying TSV filter: {filter_expression}")

    try:
        # Read the TSV file
        compression = "gzip" if input_tsv.endswith(".gz") else None
        df = pd.read_csv(
            input_tsv,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            compression=cast(Any, compression),
        )
        initial_count = len(df)

        if initial_count == 0:
            logger.warning("Input TSV is empty, writing empty output")
            compression_out = "gzip" if output_tsv.endswith(".gz") else None
            df.to_csv(output_tsv, sep="\t", index=False, compression=cast(Any, compression_out))
            return

        if not pandas_query:
            # Convert SnpSift-style expression to pandas query
            # This is a simplified converter - extend as needed
            expression = filter_expression

            # Handle basic conversions
            expression = expression.replace("=", "==")
            expression = expression.replace("!==", "!=")
            expression = re.sub(r"\b&\b", " & ", expression)
            expression = re.sub(r"\b\|\b", " | ", expression)

            # Convert numeric comparisons for score columns
            # Assume columns ending with _Score or Score are numeric
            score_columns = [
                col for col in df.columns if col.endswith("Score") or col.endswith("_score")
            ]
            for col in score_columns:
                # Convert string column to numeric for comparison
                df[col] = pd.to_numeric(df[col], errors="coerce")

            filter_expression = expression

        # Apply the filter
        try:
            filtered_df = df.query(filter_expression)
        except Exception as e:
            message = f"Invalid TSV filter expression '{filter_expression}': {e}"
            logger.error(message)
            raise ValueError(message) from e

        final_count = len(filtered_df)
        pct = final_count / initial_count * 100
        logger.info(f"Filter retained {final_count}/{initial_count} variants ({pct:.1f}%)")

        # Write the filtered data
        compression_out = "gzip" if output_tsv.endswith(".gz") else None
        filtered_df.to_csv(
            output_tsv, sep="\t", index=False, compression=cast(Any, compression_out)
        )

    except Exception as e:
        if isinstance(e, ValueError):
            raise
        message = f"Error in TSV filtering: {e}"
        logger.error(message)
        raise RuntimeError(message) from e


def filter_dataframe_with_query(df: pd.DataFrame, filter_expression: str) -> pd.DataFrame:
    """
    Filter a pandas DataFrame using a query expression.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame.
    filter_expression : str
        The query string to apply.

    Returns
    -------
    pd.DataFrame
        The filtered DataFrame.
    """
    if not filter_expression:
        return df

    logger.info(f"Applying final filter expression: {filter_expression}")
    original_count = len(df)

    try:
        # Create a copy to avoid modifying the original DataFrame
        df_copy = df.copy()

        # Convert string/object columns to numeric types for comparisons like > or <.
        # In pandas 3.0+, string columns use "str" dtype instead of "object".
        for col in df_copy.columns:
            if pd.api.types.is_string_dtype(df_copy[col]) or df_copy[col].dtype == object:
                converted = pd.to_numeric(df_copy[col], errors="coerce")
                if converted.notna().any():
                    df_copy[col] = converted

        filtered_df = df_copy.query(filter_expression)
        final_count = len(filtered_df)
        logger.info(f"Final filter retained {final_count} of {original_count} variants.")
        return filtered_df

    except Exception as e:
        message = f"Invalid final filter expression '{filter_expression}': {e}"
        logger.error(message)
        raise ValueError(message) from e
