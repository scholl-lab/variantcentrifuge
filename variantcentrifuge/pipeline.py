# File: variantcentrifuge/pipeline.py
# Location: variantcentrifuge/variantcentrifuge/pipeline.py

"""
Pipeline orchestration module for variantcentrifuge.

This module provides high-level orchestration of the analysis steps:
- Checks external tools
- Normalizes genes
- Loads phenotypes and phenotype terms
- Extracts variants and applies filters
- Performs optional genotype replacement
- Integrates phenotype data
- Runs variant-level and gene-burden analyses using analyze_variants
- Generates metadata and optionally converts results to Excel
- Cleans up intermediates if requested

All steps are coordinated within the run_pipeline function.
"""

import argparse
import datetime
import hashlib
import logging
import locale
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from .analyze_variants import analyze_variants
from .annotator import (
    load_custom_features,
    annotate_dataframe_with_features,
    validate_annotation_config,
)
from .converter import (
    append_tsv_as_sheet,
    convert_to_excel,
    finalize_excel_file,
    produce_report_json,
)
from .extractor import extract_fields
from .filters import (
    apply_snpsift_filter,
    extract_variants,
    filter_final_tsv_by_genotype,
    filter_tsv_with_expression,
)
from .gene_bed import get_gene_bed, normalize_genes
from .links import add_links_to_table
from .phenotype import aggregate_phenotypes_for_samples, load_phenotypes
from .phenotype_filter import filter_phenotypes
from .ped_reader import read_pedigree
from .replacer import replace_genotypes
from .scoring import read_scoring_config
from .utils import (
    check_external_tools,
    ensure_fields_in_extract,
    get_tool_version,
    normalize_snpeff_headers,
    run_command,
    sanitize_metadata_field,
)
from variantcentrifuge.helpers import (
    annotate_variants_with_gene_lists,
    check_file,
    determine_case_control_sets,
    dump_df_to_xlsx,
    extract_gencode_id,
    get_vcf_names,
    get_vcf_regions,
    get_vcf_samples,
    get_vcf_size,
    match_IGV_link_columns,
    read_sequencing_manifest,
)

# Import the SNPeff annotation splitting function
from .vcf_eff_one_per_line import process_vcf_file

logger = logging.getLogger("variantcentrifuge")


def remove_vcf_extensions(filename: str) -> str:
    """
    Remove common VCF-related extensions from a filename.

    Parameters
    ----------
    filename : str
        The input filename, possibly ending in .vcf, .vcf.gz, or .gz.

    Returns
    -------
    str
        The filename base without VCF-related extensions.
    """
    if filename.endswith(".vcf.gz"):
        return filename[:-7]
    elif filename.endswith(".vcf"):
        return filename[:-4]
    elif filename.endswith(".gz"):
        return filename[:-3]
    return filename


def compute_base_name(vcf_path: str, gene_name: str) -> str:
    """
    Compute a base name for output files based on the VCF filename and genes.

    If multiple genes are specified, create a hash to represent them.
    If 'all' is specified, append '.all'.
    Otherwise, append the gene name if it's not already in the VCF base name.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.
    gene_name : str
        The normalized gene name string.

    Returns
    -------
    str
        A base name for output files.
    """
    genes = gene_name.strip()
    vcf_base = os.path.basename(vcf_path)
    vcf_base = remove_vcf_extensions(vcf_base)

    if genes.lower() == "all":
        return f"{vcf_base}.all"
    split_genes = genes.split()
    if len(split_genes) > 1:
        gene_hash = hashlib.md5(genes.encode("utf-8")).hexdigest()[:8]
        return f"{vcf_base}.multiple-genes-{gene_hash}"
    else:
        if split_genes and split_genes[0].lower() in vcf_base.lower():
            return vcf_base
        else:
            return f"{vcf_base}.{split_genes[0]}" if split_genes else vcf_base


def load_terms_from_file(file_path: Optional[str], logger: logging.Logger) -> List[str]:
    """
    Load terms (HPO terms, sample IDs, etc.) from a file, one per line.

    Parameters
    ----------
    file_path : str or None
        Path to a file containing one term per line.
    logger : logging.Logger
        Logger instance for error logging.

    Returns
    -------
    list of str
        A list of terms loaded from the file.

    Raises
    ------
    SystemExit
        If the file is missing or empty and a file_path was specified.
    """
    terms: List[str] = []
    if file_path:
        if not os.path.exists(file_path):
            logger.error(f"Required file not found: {file_path}")
            sys.exit(1)
        with open(file_path, "r", encoding="utf-8") as f:
            found_any = False
            for line in f:
                t = line.strip()
                if t:
                    found_any = True
                    terms.append(t)
            if not found_any:
                logger.error(f"File {file_path} is empty or invalid.")
                sys.exit(1)
    return terms


def parse_samples_from_vcf(vcf_file: str) -> List[str]:
    """
    Parse sample names from a VCF file by reading its header line.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file.

    Returns
    -------
    list of str
        A list of sample names extracted from the VCF header.
    """
    output = run_command(["bcftools", "view", "-h", vcf_file], output_file=None)
    chrom_line = None
    for line in output.splitlines():
        if line.startswith("#CHROM"):
            chrom_line = line
            break
    if not chrom_line:
        return []
    fields = chrom_line.strip().split("\t")
    if len(fields) <= 9:
        return []
    samples = fields[9:]
    return samples


def run_pipeline(
    args: argparse.Namespace, cfg: Dict[str, Any], start_time: datetime.datetime
) -> None:
    """
    High-level orchestration of the pipeline steps.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    cfg : dict
        Configuration dictionary merged from CLI and config file.
    start_time : datetime.datetime
        The start time of the run.

    Returns
    -------
    None
        Writes output files and may print results to stdout.
    """
    # Check external tools
    check_external_tools()

    # Normalize genes
    gene_name = normalize_genes(args.gene_name if args.gene_name else "", args.gene_file, logger)
    logger.debug(f"Normalized gene list: {gene_name}")

    # Determine base name and output paths
    if args.output_file is not None:
        if args.output_file in ["stdout", "-"]:
            final_to_stdout = True
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = None
        else:
            final_to_stdout = False
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = os.path.join(args.output_dir, os.path.basename(args.output_file))
    else:
        final_to_stdout = False
        base_name = compute_base_name(args.vcf_file, gene_name)
        final_output = os.path.join(args.output_dir, f"{base_name}.final.tsv")

    os.makedirs(args.output_dir, exist_ok=True)
    intermediate_dir = os.path.join(args.output_dir, "intermediate")
    os.makedirs(intermediate_dir, exist_ok=True)

    if not cfg["no_stats"] and not args.stats_output_file:
        cfg["stats_output_file"] = os.path.join(intermediate_dir, f"{base_name}.statistics.tsv")
    else:
        cfg["stats_output_file"] = args.stats_output_file

    # Load phenotypes if provided
    phenotypes = {}
    use_phenotypes = False
    if args.phenotype_file and args.phenotype_sample_column and args.phenotype_value_column:
        phenotypes = load_phenotypes(
            args.phenotype_file,
            args.phenotype_sample_column,
            args.phenotype_value_column,
        )
        if not phenotypes:
            logger.error(
                f"No phenotype data loaded from {args.phenotype_file}. "
                "Check file formatting and column names."
            )
            sys.exit(1)
        use_phenotypes = True

    # Store phenotypes in cfg for determine_case_control_sets usage
    cfg["phenotypes"] = phenotypes

    # Load scoring configuration if path is provided
    scoring_config = None
    if cfg.get("scoring_config_path"):
        try:
            scoring_config = read_scoring_config(cfg["scoring_config_path"])
            logger.info(
                f"Successfully loaded scoring configuration from {cfg['scoring_config_path']}"
            )
        except Exception as e:
            logger.error(f"Failed to load scoring configuration: {e}")
            sys.exit(1)

    # Load pedigree data if provided
    pedigree_data = None
    if cfg.get("calculate_inheritance"):
        if cfg.get("ped_file"):
            try:
                pedigree_data = read_pedigree(cfg["ped_file"])
                cfg["pedigree_data"] = pedigree_data
                logger.info(f"Loaded pedigree data for {len(pedigree_data)} individuals.")
            except Exception as e:
                logger.error(f"Failed to load PED file: {e}")
                sys.exit(1)
        else:
            # Create single-sample pedigree data for inheritance analysis
            logger.info("No PED file provided - treating all samples as affected individuals.")
            # When no PED file is provided, create empty pedigree_data
            # Will be populated later when we have sample information
            pedigree_data = {}
            cfg["pedigree_data"] = pedigree_data

    # Validate and load custom annotation features
    annotation_errors = validate_annotation_config(cfg)
    if annotation_errors:
        logger.error("Custom annotation configuration errors:")
        for error in annotation_errors:
            logger.error(f"  - {error}")
        sys.exit(1)

    custom_features = load_custom_features(cfg)

    # Load phenotype terms
    case_hpo_terms = []
    control_hpo_terms = []
    if args.case_phenotypes:
        case_hpo_terms = [t.strip() for t in args.case_phenotypes.split(",") if t.strip()]
    case_hpo_terms += load_terms_from_file(args.case_phenotypes_file, logger)

    if args.control_phenotypes:
        control_hpo_terms = [t.strip() for t in args.control_phenotypes.split(",") if t.strip()]
    control_hpo_terms += load_terms_from_file(args.control_phenotypes_file, logger)

    cfg["case_phenotypes"] = case_hpo_terms
    cfg["control_phenotypes"] = control_hpo_terms

    # Load explicit sample sets
    case_samples = []
    control_samples = []
    if args.case_samples:
        case_samples = [s.strip() for s in args.case_samples.split(",") if s.strip()]
    case_samples += load_terms_from_file(args.case_samples_file, logger)

    if args.control_samples:
        control_samples = [s.strip() for s in args.control_samples.split(",") if s.strip()]
    control_samples += load_terms_from_file(args.control_samples_file, logger)

    cfg["case_samples"] = case_samples
    cfg["control_samples"] = control_samples

    # Extract gene BED
    bed_file = get_gene_bed(
        cfg["reference"],
        gene_name,
        interval_expand=cfg.get("interval_expand", 0),
        add_chr=cfg.get("add_chr", True),
        output_dir=args.output_dir,
    )
    logger.debug(f"Gene BED created at: {bed_file}")

    if not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
        logger.error("Gene BED file could not be created or is empty. Check genes or reference.")
        sys.exit(1)

    # Check if requested genes found in reference
    if gene_name.lower() != "all":
        requested_genes = gene_name.split()
        found_genes = set()
        with open(bed_file, "r", encoding="utf-8") as bf:
            for line in bf:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    g_name = parts[3].split(";")[0].strip()
                    found_genes.add(g_name)
        missing = [g for g in requested_genes if g not in found_genes]
        if missing:
            logger.warning(
                "The following gene(s) were not found in the reference: " f"{', '.join(missing)}"
            )
            if cfg.get("debug_level", "INFO") == "ERROR":
                sys.exit(1)

    # Filenames for intermediate steps
    variants_file = os.path.join(intermediate_dir, f"{base_name}.variants.vcf.gz")
    splitted_before_file = os.path.join(
        intermediate_dir, f"{base_name}.splitted_before_filters.vcf.gz"
    )
    splitted_after_file = os.path.join(
        intermediate_dir, f"{base_name}.splitted_after_filters.vcf.gz"
    )
    filtered_file = os.path.join(intermediate_dir, f"{base_name}.filtered.vcf.gz")
    transcript_filtered_file = os.path.join(
        intermediate_dir, f"{base_name}.transcript_filtered.vcf.gz"
    )
    extracted_tsv = os.path.join(intermediate_dir, f"{base_name}.extracted.tsv")
    genotype_replaced_tsv = os.path.join(intermediate_dir, f"{base_name}.genotype_replaced.tsv")
    phenotype_added_tsv = os.path.join(intermediate_dir, f"{base_name}.phenotypes_added.tsv")
    gene_burden_tsv = os.path.join(args.output_dir, f"{base_name}.gene_burden.tsv")

    # Parse samples from VCF
    original_samples = parse_samples_from_vcf(args.vcf_file)
    if cfg.get("remove_sample_substring"):
        substring_to_remove = cfg["remove_sample_substring"]
        if substring_to_remove and substring_to_remove.strip():
            original_samples = [s.replace(substring_to_remove, "") for s in original_samples]

    # If no PED file was provided, populate pedigree_data with samples as affected
    if pedigree_data is not None and not pedigree_data and cfg.get("calculate_inheritance", False):
        logger.debug("Populating pedigree data for single sample analysis")
        for sample in original_samples:
            pedigree_data[sample] = {
                "sample_id": sample,
                "father_id": "0",
                "mother_id": "0",
                "sex": "0",  # Unknown sex
                "affected_status": "2",  # Affected
            }

    # -----------------------------------------------------------------------
    # Step 1: Extract variants => variants_file
    # -----------------------------------------------------------------------
    extract_variants(args.vcf_file, bed_file, cfg, variants_file)

    # Handle snpeff splitting mode (None, 'before_filters', or 'after_filters')
    splitting_mode = cfg.get("snpeff_splitting_mode", None)

    if splitting_mode == "before_filters":
        logger.info("Splitting multiple SNPeff (EFF/ANN) annotations before main filtering.")
        process_vcf_file(variants_file, splitted_before_file)
        prefiltered_for_snpsift = splitted_before_file
    else:
        # either None or 'after_filters'
        prefiltered_for_snpsift = variants_file

    # -----------------------------------------------------------------------
    # Step 3: Apply SnpSift filter => filtered_file (skip if late filtering is enabled)
    # -----------------------------------------------------------------------
    if cfg.get("late_filtering", False):
        logger.info("Late filtering enabled - skipping early SnpSift filter step")
        # Just copy the file without filtering
        import shutil

        shutil.copy2(prefiltered_for_snpsift, filtered_file)
    else:
        apply_snpsift_filter(prefiltered_for_snpsift, cfg["filters"], cfg, filtered_file)

    # If splitting after filters
    if splitting_mode == "after_filters":
        logger.info(
            "Splitting multiple SNPeff (EFF/ANN) annotations after main filters, before transcript filter."
        )
        process_vcf_file(filtered_file, splitted_after_file)
        post_filter_for_transcripts = splitted_after_file
    else:
        post_filter_for_transcripts = filtered_file

    # -----------------------------------------------------------------------
    # Step 4: If transcripts are provided, filter by transcripts => transcript_filtered_file
    # -----------------------------------------------------------------------
    transcripts = []
    if args.transcript_list:
        transcripts.extend([t.strip() for t in args.transcript_list.split(",") if t.strip()])
    if args.transcript_file:
        if not os.path.exists(args.transcript_file):
            logger.error(f"Transcript file not found: {args.transcript_file}")
            sys.exit(1)
        with open(args.transcript_file, "r", encoding="utf-8") as tf:
            for line in tf:
                tr = line.strip()
                if tr:
                    transcripts.append(tr)

    transcripts = list(set(transcripts))  # remove duplicates

    if transcripts:
        logger.info("Filtering for transcripts using SnpSift.")
        or_clauses = [f"(EFF[*].TRID = '{tr}')" for tr in transcripts]
        transcript_filter_expr = " | ".join(or_clauses)

        apply_snpsift_filter(
            post_filter_for_transcripts,
            transcript_filter_expr,
            cfg,
            transcript_filtered_file,
        )
        final_filtered_for_extraction = transcript_filtered_file
    else:
        final_filtered_for_extraction = post_filter_for_transcripts

    # -----------------------------------------------------------------------
    # Step 5: Extract fields => extracted_tsv
    # -----------------------------------------------------------------------
    # If user wants to append extra sample fields, ensure they're in the main fields
    # (Unify them into cfg["fields_to_extract"])
    if cfg.get("append_extra_sample_fields", False) and cfg.get("extra_sample_fields", []):
        updated_fields = ensure_fields_in_extract(
            cfg["fields_to_extract"], cfg["extra_sample_fields"]
        )
        cfg["fields_to_extract"] = updated_fields
        logger.debug(
            f"Updated fields_to_extract with extra sample fields: {cfg['fields_to_extract']}"
        )

    field_list = " ".join((cfg["fields_to_extract"] or args.fields).strip().split())
    logger.debug(f"Extracting fields: {field_list} -> {extracted_tsv}")

    snpsift_sep = cfg.get("extract_fields_separator", ":")
    logger.debug(f"Using SnpSift subfield separator for extractFields: '{snpsift_sep}'")

    temp_cfg = cfg.copy()
    temp_cfg["extract_fields_separator"] = snpsift_sep
    extract_fields(final_filtered_for_extraction, field_list, temp_cfg, extracted_tsv)

    # Genotype replacement logic
    with open(extracted_tsv, "r", encoding="utf-8") as f:
        header_line = f.readline().strip()
    columns = header_line.split("\t")
    gt_present = "GT" in columns

    cfg["sample_list"] = ",".join(original_samples) if original_samples else ""
    if not args.no_replacement and gt_present:
        with open(extracted_tsv, "r", encoding="utf-8") as inp, open(
            genotype_replaced_tsv, "w", encoding="utf-8"
        ) as out:
            for line in replace_genotypes(inp, cfg):
                out.write(line + "\n")
        replaced_tsv = genotype_replaced_tsv
    else:
        replaced_tsv = extracted_tsv

    # If user has appended extra fields, they might want them removed after genotype assembly
    if cfg.get("append_extra_sample_fields", False) and cfg.get("extra_sample_fields"):
        logger.debug(
            "User config => removing columns after replacement: %s",
            cfg["extra_sample_fields"],
        )
        stripped_tsv = os.path.join(intermediate_dir, f"{base_name}.stripped_extras.tsv")

        with open(replaced_tsv, "r", encoding="utf-8") as inp, open(
            stripped_tsv, "w", encoding="utf-8"
        ) as out:
            # Read the header line from replaced_tsv
            raw_header_line = next(inp).rstrip("\n")
            original_header_cols = raw_header_line.split("\t")

            # Create a normalized copy of the header to see how columns might have been renamed
            # by SnpSift or other tools:
            # Keep normalized header code for potential future use
            # normed_header_line = normalize_snpeff_headers([raw_header_line])[0]

            # Build a map from normalized_col -> original_index
            normalized_to_index = {}
            for i, col in enumerate(original_header_cols):
                # Single-column approach to re-using normalize_snpeff_headers:
                normed_col = normalize_snpeff_headers([col])[0]
                normalized_to_index[normed_col] = i

            remove_indices = []
            # For each user-supplied field, normalize it and see if it matches a normalized header col
            for raw_col_name in cfg["extra_sample_fields"]:
                single_line = normalize_snpeff_headers([raw_col_name])[0]
                if single_line in normalized_to_index:
                    idx = normalized_to_index[single_line]
                    remove_indices.append(idx)
                    logger.debug(
                        "Removing normalized column '%s' => real index %d",
                        single_line,
                        idx,
                    )
                else:
                    logger.warning(
                        "Column '%s' was requested for removal but not found in header!",
                        raw_col_name,
                    )

            remove_indices.sort(reverse=True)
            new_header = [h for i, h in enumerate(original_header_cols) if i not in remove_indices]
            out.write("\t".join(new_header) + "\n")

            for line_num, line in enumerate(inp, start=2):
                line = line.rstrip("\n")
                if not line.strip():
                    out.write(line + "\n")
                    continue
                parts = line.split("\t")
                for idx in remove_indices:
                    if idx < len(parts):
                        parts.pop(idx)
                out.write("\t".join(parts) + "\n")

        os.rename(stripped_tsv, replaced_tsv)
    logger.debug("Extra column removal (if requested) complete.")

    # Integrate phenotypes if provided
    if use_phenotypes:
        pattern = re.compile(r"^([^()]+)(?:\([^)]+\))?$")

        with open(replaced_tsv, "r", encoding="utf-8") as inp, open(
            phenotype_added_tsv, "w", encoding="utf-8"
        ) as out:
            header = next(inp).rstrip("\n")
            header_fields = header.split("\t")
            header_fields.append("phenotypes")
            out.write("\t".join(header_fields) + "\n")

            wrote_data = False
            gt_idx = header_fields.index("GT") if "GT" in header_fields else None
            for line in inp:
                line = line.rstrip("\n")
                if not line.strip():
                    out.write(line + "\n")
                    continue
                fields_line = line.split("\t")

                if gt_idx is not None and gt_idx < len(fields_line):
                    gt_value = fields_line[gt_idx]
                    samples_in_line = []
                    if gt_value.strip():
                        sample_entries = gt_value.split(";")
                        for entry in sample_entries:
                            entry = entry.strip()
                            if entry:
                                m = pattern.match(entry)
                                if m:
                                    s = m.group(1).strip()
                                    if s:
                                        samples_in_line.append(s)
                    pheno_str = ""
                    if samples_in_line:
                        pheno_str = aggregate_phenotypes_for_samples(samples_in_line, phenotypes)
                    fields_line.append(pheno_str)
                else:
                    fields_line.append("")
                out.write("\t".join(fields_line) + "\n")
                wrote_data = True

            if not wrote_data:
                logger.warning(
                    "Phenotype integration found no data rows to process. Only header row will be written."
                )
        final_tsv = phenotype_added_tsv
    else:
        final_tsv = replaced_tsv

    # Genotype filtering if requested
    if getattr(args, "genotype_filter", None) or getattr(args, "gene_genotype_file", None):
        genotype_filtered_tsv = os.path.join(args.output_dir, f"{base_name}.genotype_filtered.tsv")
        genotype_modes = set()
        if getattr(args, "genotype_filter", None):
            genotype_modes = set(g.strip() for g in args.genotype_filter.split(",") if g.strip())
        filter_final_tsv_by_genotype(
            input_tsv=final_tsv,
            output_tsv=genotype_filtered_tsv,
            global_genotypes=genotype_modes,
            gene_genotype_file=args.gene_genotype_file,
        )
        final_tsv = genotype_filtered_tsv

    # If Excel requested
    if args.xlsx:
        excel_file = os.path.splitext(final_output or f"{base_name}.final.tsv")[0] + ".xlsx"
        cfg["final_excel_file"] = excel_file
    else:
        cfg["final_excel_file"] = None

    # Load TSV into DataFrame for unified annotation and analysis
    import pandas as pd

    try:
        df = pd.read_csv(final_tsv, sep="\t", dtype=str, keep_default_na=False)
        logger.info(f"Loaded {len(df)} variants for annotation and analysis")
    except Exception as e:
        logger.error(f"Failed to load TSV file {final_tsv}: {e}")
        sys.exit(1)

    # Apply unified custom annotations
    has_custom_features = (
        bool(custom_features.get("regions_by_chrom"))
        or bool(custom_features.get("gene_lists"))
        or bool(custom_features.get("json_gene_data"))
    )

    if has_custom_features:
        logger.info("Applying unified custom annotations...")
        df = annotate_dataframe_with_features(df, custom_features)
    else:
        # Add empty Custom_Annotation column for consistency
        df["Custom_Annotation"] = ""

    # If inheritance analysis is requested, augment DataFrame with VCF genotypes
    if cfg.get("calculate_inheritance", False) and pedigree_data is not None:
        logger.info("Augmenting DataFrame with VCF genotypes for inheritance analysis...")
        from .vcf_genotype_extractor import augment_dataframe_with_vcf_genotypes

        # Use the final filtered VCF file for genotype extraction
        vcf_for_genotypes = final_filtered_for_extraction
        df = augment_dataframe_with_vcf_genotypes(df, vcf_for_genotypes, original_samples)

    # Convert DataFrame back to file format for analyze_variants
    annotated_tsv = os.path.join(intermediate_dir, f"{base_name}.annotated.tsv")
    df.to_csv(annotated_tsv, sep="\t", index=False, na_rep="")

    # Run analyze_variants for variant-level results
    temp_cfg = cfg.copy()
    temp_cfg["perform_gene_burden"] = False
    temp_cfg["scoring_config"] = scoring_config  # Pass the loaded config
    temp_cfg["pedigree_data"] = pedigree_data  # Pass the loaded pedigree data
    temp_cfg["calculate_inheritance"] = cfg.get("calculate_inheritance", False)
    temp_cfg["sample_list"] = cfg.get("sample_list", "")  # Ensure sample_list is passed
    buffer = []
    with open(annotated_tsv, "r", encoding="utf-8") as inp:
        line_count = 0
        for line in analyze_variants(inp, temp_cfg):
            buffer.append(line)
            line_count += 1
        if line_count <= 1:  # Only header or nothing
            logger.warning(
                "No variant-level results produced after filtering. Generating empty output files with headers."
            )

    # Note: Gene list annotation is now handled by the unified annotation system above
    # This ensures --annotate-gene-list files are processed through the new Custom_Annotation column
    # The old separate column system is deprecated but kept here for backward compatibility checking
    gene_list_annotation_files = cfg.get("annotate_gene_list_files", [])
    if gene_list_annotation_files and len(buffer) > 1:
        logger.warning(
            "Gene list annotation is now integrated into the unified Custom_Annotation system. "
            "Individual gene list columns are no longer added separately."
        )

    # Apply late filtering if enabled and filters are specified
    if cfg.get("late_filtering", False) and cfg.get("filters") and len(buffer) > 1:
        logger.info("Applying late filtering on scored and annotated data...")

        # Write buffer to a temporary file
        temp_scored_tsv = os.path.join(intermediate_dir, f"{base_name}.scored.tsv")
        with open(temp_scored_tsv, "w", encoding="utf-8") as f:
            for line in buffer:
                f.write(line + "\n")

        # Apply the filter
        temp_filtered_tsv = os.path.join(intermediate_dir, f"{base_name}.late_filtered.tsv")
        filter_tsv_with_expression(
            temp_scored_tsv,
            temp_filtered_tsv,
            cfg["filters"],
            pandas_query=False,  # Use SnpSift-style syntax
        )

        # Read the filtered results back into buffer
        with open(temp_filtered_tsv, "r", encoding="utf-8") as f:
            buffer = [line.rstrip("\n") for line in f]

        if len(buffer) <= 1:
            logger.warning("No variants passed late filtering")

    # Add links if not disabled and we have data rows
    if not cfg.get("no_links", False) and len(buffer) > 1:
        links_config = cfg.get("links", {})
        if links_config:
            logger.debug("Adding link columns to final output.")
            buffer = add_links_to_table(buffer, links_config)
    else:
        logger.debug("Link columns are disabled by configuration or no variants to link.")

    # Insert a leading VAR_ID column
    def add_variant_identifier(lines: List[str]) -> List[str]:
        """
        Insert a leading VAR_ID column into each line.

        The VAR_ID is generated from a zero-padded incrementing counter plus
        a short 4-character hash based on CHROM, POS, REF, and ALT.

        Parameters
        ----------
        lines : list of str
            The lines of the final variant table.

        Returns
        -------
        list of str
            A new list of lines with a VAR_ID column inserted at the beginning.
        """
        new_lines = []
        header_indices = {}
        row_counter = 1

        # Identify header and relevant columns
        header_line = lines[0].rstrip("\n")
        header_parts = header_line.split("\t")
        for idx, col in enumerate(header_parts):
            header_indices[col] = idx

        # Prepend header
        new_header = ["VAR_ID"] + header_parts
        new_lines.append("\t".join(new_header))

        chrom_idx = header_indices.get("CHROM", None)
        pos_idx = header_indices.get("POS", None)
        ref_idx = header_indices.get("REF", None)
        alt_idx = header_indices.get("ALT", None)

        for line in lines[1:]:
            line = line.rstrip("\n")
            if not line.strip():
                new_lines.append(line)
                continue

            fields = line.split("\t")
            chrom_val = fields[chrom_idx] if chrom_idx is not None else ""
            pos_val = fields[pos_idx] if pos_idx is not None else ""
            ref_val = fields[ref_idx] if ref_idx is not None else ""
            alt_val = fields[alt_idx] if alt_idx is not None else ""

            combined = f"{chrom_val}{pos_val}{ref_val}{alt_val}"
            short_hash = hashlib.md5(combined.encode("utf-8")).hexdigest()[:4]

            var_id = f"var_{row_counter:04d}_{short_hash}"
            row_counter += 1
            new_line = [var_id] + fields
            new_lines.append("\t".join(new_line))

        return new_lines

    # Add variant identifiers only if we have data rows
    if len(buffer) > 1:
        buffer = add_variant_identifier(buffer)
    else:
        # For empty results, just add the VAR_ID column to the header
        if buffer and "\t" in buffer[0]:
            header_parts = buffer[0].split("\t")
            new_header = ["VAR_ID"] + header_parts
            buffer[0] = "\t".join(new_header)

    # Append named blank columns if requested
    def add_named_columns(lines: List[str], col_names: List[str]) -> List[str]:
        """
        Append columns named in col_names to the table, with blank data for each row.

        Parameters
        ----------
        lines : list of str
            Lines (header + data) of the final table.
        col_names : list of str
            The column names to add.

        Returns
        -------
        list of str
            The updated table lines with extra named blank columns.
        """
        if not lines or not col_names:
            return lines

        new_lines = []
        header_line = lines[0].rstrip("\n")
        header_parts = header_line.split("\t")
        header_parts.extend(col_names)
        new_header = "\t".join(header_parts)
        new_lines.append(new_header)

        for line in lines[1:]:
            line = line.rstrip("\n")
            if not line.strip():
                new_lines.append(line)
                continue
            fields = line.split("\t")
            fields.extend([""] * len(col_names))
            new_lines.append("\t".join(fields))

        return new_lines

    if args.add_column:
        buffer = add_named_columns(buffer, args.add_column)

    if final_to_stdout:
        if buffer:
            sys.stdout.write("\n".join(buffer) + "\n")
        final_out_path = None
    else:
        final_file = final_output
        with open(final_file, "w", encoding="utf-8") as out:
            for line in buffer:
                out.write(line + "\n")
        final_out_path = final_file

    # Gene burden if requested
    if cfg.get("perform_gene_burden", False):
        gene_burden_cfg = cfg.copy()
        gene_burden_cfg["perform_gene_burden"] = True
        line_count = 0
        with open(final_tsv, "r", encoding="utf-8") as inp, open(
            gene_burden_tsv, "w", encoding="utf-8"
        ) as out:
            for line in analyze_variants(inp, gene_burden_cfg):
                out.write(line + "\n")
                line_count += 1
        if line_count <= 1:  # Only header or nothing
            logger.warning("Gene burden requested but no variant data produced.")
        else:
            logger.debug(f"Gene burden analysis complete: {gene_burden_tsv}")

    end_time = datetime.datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"Run ended at {end_time.isoformat()}, duration: {duration} seconds")

    # Metadata
    metadata_file = os.path.join(args.output_dir, f"{base_name}.metadata.tsv")
    with open(metadata_file, "w", encoding="utf-8") as mf:
        mf.write("Parameter\tValue\n")

        def meta_write(param, val):
            p = sanitize_metadata_field(param)
            v = sanitize_metadata_field(val)
            mf.write(f"{p}\t{v}\n")

        meta_write("Tool", "variantcentrifuge")
        meta_write("Version", cfg.get("version", "N/A"))
        meta_write("Run_start_time", start_time.isoformat())
        meta_write("Run_end_time", end_time.isoformat())
        meta_write("Run_duration_seconds", str(duration))
        meta_write("Date", datetime.datetime.now().isoformat())
        meta_write("Command_line", " ".join([sanitize_metadata_field(x) for x in sys.argv]))

        for k, v in cfg.items():
            meta_write(f"config.{k}", str(v))

        # Versions
        snpeff_ver = get_tool_version("snpEff")
        bcftools_ver = get_tool_version("bcftools")
        snpsift_ver = get_tool_version("SnpSift")
        bedtools_ver = get_tool_version("bedtools")

        meta_write("tool.snpeff_version", snpeff_ver)
        meta_write("tool.bcftools_version", bcftools_ver)
        meta_write("tool.snpsift_version", snpsift_ver)
        meta_write("tool.bedtools_version", bedtools_ver)

    # Define report directory consistently
    report_dir = os.path.join(args.output_dir, "report")
    os.makedirs(report_dir, exist_ok=True)  # Ensure report_dir exists

    # Phase 1: Excel conversion (without finalization) if requested
    xlsx_file = None
    if args.xlsx and final_out_path and not final_to_stdout:
        if not os.path.exists(final_out_path) or os.path.getsize(final_out_path) == 0:
            logger.warning("Final output file is empty. Cannot convert to Excel.")
        else:
            xlsx_file = convert_to_excel(final_out_path, cfg)
            append_tsv_as_sheet(xlsx_file, metadata_file, sheet_name="Metadata")
            if (
                not cfg["no_stats"]
                and cfg.get("stats_output_file")
                and os.path.exists(cfg["stats_output_file"])
            ):
                if os.path.getsize(cfg["stats_output_file"]) > 0:
                    append_tsv_as_sheet(
                        xlsx_file, cfg["stats_output_file"], sheet_name="Statistics"
                    )
                else:
                    logger.warning("Stats file is empty, skipping Statistics sheet.")
            if (
                cfg.get("perform_gene_burden", False)
                and os.path.exists(gene_burden_tsv)
                and os.path.getsize(gene_burden_tsv) > 0
            ):
                append_tsv_as_sheet(xlsx_file, gene_burden_tsv, sheet_name="Gene Burden")

    # Phase 2: Generate IGV reports and map if IGV is enabled
    # This must happen before finalize_excel_file and produce_report_json
    igv_enabled = cfg.get("igv_enabled", False)
    # Check if the file exists and has content (at least header)
    if igv_enabled and final_out_path and os.path.exists(final_out_path):
        # MODIFIED: Start of local IGV FASTA feature
        bam_map_file = cfg.get("bam_mapping_file")
        igv_reference_genome = cfg.get("igv_reference")
        igv_fasta_file = cfg.get("igv_fasta")
        igv_ideogram_file = cfg.get("igv_ideogram")

        # Validate the required parameters and files
        if not bam_map_file:
            logger.error("For IGV integration, --bam-mapping-file must be provided.")
            sys.exit(1)

        if not igv_fasta_file and not igv_reference_genome:
            logger.error(
                "For IGV integration, either --igv-reference or --igv-fasta must be provided."
            )
            sys.exit(1)

        # Validate local FASTA files if provided
        if igv_fasta_file or igv_ideogram_file:
            from .validators import validate_igv_files

            # Get the FASTA index file if provided
            # FASTA index now follows standard naming convention (FASTA file + .fai)
            validate_igv_files(igv_fasta_file, igv_ideogram_file)

        if igv_fasta_file and igv_reference_genome:
            logger.warning(
                "Both local FASTA file and genome reference ID provided. Local FASTA file will take precedence."
            )
        # MODIFIED: End of local IGV FASTA feature

        from .generate_igv_report import generate_igv_report  # Ensure import is present

        # The output_dir for generate_igv_report should be where igv_reports_map.json is expected
        # which is output_dir/report, and report_path in map are relative to output_dir/report
        # MODIFIED: Start of local IGV FASTA feature
        generate_igv_report(
            variants_tsv=final_out_path,  # final_out_path is the path to the main results TSV
            output_dir=report_dir,  # This ensures map is in report/igv/
            bam_mapping_file=bam_map_file,
            igv_reference=igv_reference_genome,
            integrate_into_main=True,  # Always integrate if IGV is enabled
            igv_fasta=igv_fasta_file,
            igv_ideogram=igv_ideogram_file,
            # MODIFIED: Pass filename shortening parameters
            igv_max_allele_len_filename=cfg.get("igv_max_allele_len_filename", 10),
            igv_hash_len_filename=cfg.get("igv_hash_len_filename", 6),
            igv_max_variant_part_filename=cfg.get("igv_max_variant_part_filename", 50),
            # MODIFIED: Start of IGV flanking feature - pass configurable flanking value
            igv_flanking=cfg.get("igv_flanking", 50),
            # MODIFIED: End of IGV flanking feature
        )
        # MODIFIED: End of local IGV FASTA feature
        logger.info("IGV reports and mapping file generated.")

    # Phase 3: Produce HTML report if requested - after IGV reports are generated
    if args.html_report and final_out_path and os.path.exists(final_out_path):
        # produce_report_json will look for the IGV map created in Phase 2
        produce_report_json(final_out_path, args.output_dir)

        from .generate_html_report import generate_html_report  # Ensure import

        variants_json = os.path.join(report_dir, "variants.json")
        summary_json = os.path.join(report_dir, "summary.json")
        generate_html_report(
            variants_json=variants_json,
            summary_json=summary_json,
            output_dir=report_dir,  # HTML report itself is in report_dir
            cfg=cfg,  # Pass the main configuration dictionary
        )
        logger.info("HTML report generated successfully.")

    # Phase 4: Finalize Excel file after IGV reports are generated
    if xlsx_file:
        # Now finalize the Excel file with awareness of IGV integration
        finalize_excel_file(xlsx_file, cfg)  # cfg contains igv_enabled
        logger.info("Excel file finalized with IGV links (if enabled).")

    # MODIFIED: Start of intermediate cleanup feature
    # Check if intermediate files should be kept (from CLI args or config)
    keep_intermediates = args.keep_intermediates or cfg.get("keep_intermediates", False)

    # Clean up intermediate files if not keeping them
    if not keep_intermediates:
        logger.info("Cleaning up intermediate files...")
        intermediates = [
            variants_file,
            splitted_before_file,
            splitted_after_file,
            filtered_file,
            extracted_tsv,
        ]
        if not args.no_replacement and gt_present:
            intermediates.append(genotype_replaced_tsv)
        if use_phenotypes:
            intermediates.append(phenotype_added_tsv)

        # Safety check: never delete the final output
        if final_out_path and final_out_path in intermediates:
            intermediates.remove(final_out_path)

        # Clean up intermediate directory safely
        intermediate_dir = os.path.join(args.output_dir, "intermediate")
        files_deleted = 0

        for f in intermediates:
            if f and os.path.exists(f):
                try:
                    os.remove(f)
                    files_deleted += 1
                except (IOError, OSError) as e:
                    logger.warning(f"Could not remove intermediate file {f}: {str(e)}")

        logger.info(f"Deleted {files_deleted} intermediate files.")
    else:
        logger.info("Keeping intermediate files (use --keep-intermediates=False to delete them)")
    # MODIFIED: End of intermediate cleanup feature

    logger.info(
        f"Processing completed. Output saved to " f"{'stdout' if final_to_stdout else final_output}"
    )
