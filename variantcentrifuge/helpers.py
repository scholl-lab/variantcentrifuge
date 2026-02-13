# File: variantcentrifuge/helpers.py
# Location: variantcentrifuge/variantcentrifuge/helpers.py

"""
Helper functions for analyze_variants and related processes.

Provides:
- Case/control assignment
- Phenotype map building
- Genotype parsing and allele count conversion
- Sample and phenotype classification logic
"""

import logging
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def check_file(file_path: str, exit_on_error: bool = True) -> bool:
    """
    Check if a file exists and is readable.

    Parameters
    ----------
    file_path : str
        Path to the file to check
    exit_on_error : bool, optional
        Whether to exit the program if the file does not exist, by default True

    Returns
    -------
    bool
        True if the file exists and is readable, False otherwise

    Raises
    ------
    SystemExit
        If exit_on_error is True and the file does not exist
    """
    if not os.path.exists(file_path):
        if exit_on_error:
            logger.error(f"File not found: {file_path}")
            sys.exit(1)
        return False

    if not os.path.isfile(file_path):
        if exit_on_error:
            logger.error(f"Not a file: {file_path}")
            sys.exit(1)
        return False

    return True


def compute_phenotype_based_case_control_assignment(
    vcf_samples: list[str],
    phenotype_map: dict[str, set[str]],
    case_phenotypes: list[str],
    control_phenotypes: list[str] | None = None,
    remove_substring: str = "",
) -> tuple[set[str], set[str]]:
    """
    Unified function for phenotype-based case/control assignment.

    This function handles the core logic for matching phenotype data to VCF samples
    and assigning case/control status based on phenotype criteria.

    Parameters
    ----------
    vcf_samples : List[str]
        List of sample names from VCF file
    phenotype_map : Dict[str, Set[str]]
        Map of sample_name -> set_of_phenotypes
    case_phenotypes : List[str]
        List of phenotype terms that define cases
    control_phenotypes : List[str], optional
        List of phenotype terms that define controls
    remove_substring : str, optional
        Substring to remove from sample names for matching

    Returns
    -------
    Tuple[Set[str], Set[str]]
        Tuple of (case_samples, control_samples)
    """
    if not case_phenotypes and not control_phenotypes:
        logger.debug("No phenotype criteria provided")
        return set(), set(vcf_samples)

    case_terms = set(case_phenotypes) if case_phenotypes else set()
    control_terms = set(control_phenotypes) if control_phenotypes else set()

    # Apply substring removal to VCF samples for consistent matching
    vcf_samples_processed = set(vcf_samples)
    if remove_substring and remove_substring.strip():
        vcf_samples_processed = {s.replace(remove_substring, "") for s in vcf_samples_processed}
        logger.debug(f"Applied substring removal '{remove_substring}' to VCF samples")

    classified_cases = set()
    classified_controls = set()

    for s in vcf_samples_processed:
        # Try to get phenotypes using current sample name
        phenos = phenotype_map.get(s, set())

        # If not found and substring removal was applied, try original name
        if not phenos and remove_substring:
            original_sample_name = s + remove_substring
            phenos = phenotype_map.get(original_sample_name, set())
            if phenos:
                logger.debug(f"Found phenotypes for {s} using original name {original_sample_name}")

        match_case = any(p in phenos for p in case_terms) if case_terms else False
        match_control = any(p in phenos for p in control_terms) if control_terms else False

        if case_terms and not control_terms:
            # Only case terms
            if match_case:
                classified_cases.add(s)
            else:
                classified_controls.add(s)
        elif control_terms and not case_terms:
            # Only control terms
            if match_control:
                classified_controls.add(s)
            else:
                classified_cases.add(s)
        else:
            # Both sets of terms
            if match_case and not match_control:
                classified_cases.add(s)
            elif match_control and not match_case:
                classified_controls.add(s)
            # If matches both or none, sample is not classified

    logger.debug(
        f"Phenotype-based assignment: {len(classified_cases)} cases, {len(classified_controls)} controls"
    )
    return classified_cases, classified_controls


def determine_case_control_sets(
    all_samples: set[str], cfg: dict[str, Any], df: pd.DataFrame
) -> tuple[set[str], set[str]]:
    """
    Determine case/control sample sets based on configuration.

    Logic:
    - If explicit case/control samples are provided, use them directly.
    - Else if phenotype terms are given, classify samples based on those using phenotypes from the phenotype file.
    - Else, all samples become controls.

    Parameters
    ----------
    all_samples : set of str
        The full set of sample names.
    cfg : dict
        Configuration dictionary that may include "case_samples", "control_samples",
        "case_phenotypes", "control_phenotypes", and "phenotypes" from the phenotype file.
    df : pd.DataFrame
        DataFrame with variant information (not used for phenotype assignment anymore).

    Returns
    -------
    Tuple[set, set]
        (set_of_case_samples, set_of_control_samples)

    Raises
    ------
    SystemExit
        If no valid configuration for assigning case/control sets can be determined.
    """
    logger.debug("Determining case/control sets...")

    case_samples = set(cfg.get("case_samples") or [])
    control_samples = set(cfg.get("control_samples") or [])

    # Step 1: If explicit sets are provided
    if case_samples or control_samples:
        logger.debug(
            "Explicit sample sets provided: %d cases, %d controls",
            len(case_samples),
            len(control_samples),
        )
        if case_samples and not control_samples:
            control_samples = all_samples - case_samples
        elif control_samples and not case_samples:
            case_samples = all_samples - control_samples
        return case_samples, control_samples

    # Step 2: Phenotype-based logic
    case_terms = cfg.get("case_phenotypes") or []
    control_terms = cfg.get("control_phenotypes") or []

    # If no phenotype terms are given, default to all controls
    if not case_terms and not control_terms:
        logger.debug("No phenotype terms or sets provided. All samples = controls.")
        return set(), all_samples

    # Use phenotypes from the phenotype input file (if provided)
    # If no phenotype file was provided, this will be an empty dict.
    sample_phenotype_map = cfg.get("phenotypes", {})
    logger.debug("Phenotype map has %d samples (from phenotype file).", len(sample_phenotype_map))

    # Use the unified phenotype-based assignment function
    remove_substring = cfg.get("remove_sample_substring", "")
    classified_cases, classified_controls = compute_phenotype_based_case_control_assignment(
        vcf_samples=list(all_samples),
        phenotype_map=sample_phenotype_map,
        case_phenotypes=case_terms,
        control_phenotypes=control_terms,
        remove_substring=remove_substring,
    )

    if case_terms and len(classified_cases) == 0:
        logger.warning("No samples match the case phenotype terms.")
    if control_terms and len(classified_controls) == 0:
        logger.warning("No samples match the control phenotype terms.")

    logger.debug(
        "Classified %d cases, %d controls.",
        len(classified_cases),
        len(classified_controls),
    )
    return classified_cases, classified_controls


def build_sample_phenotype_map(df: pd.DataFrame) -> dict[str, set[str]]:
    """
    Build a map of sample -> set_of_phenotypes from the 'phenotypes' column if present.

    If multiple phenotype groups appear and multiple samples in GT are present:
    - If counts match, assign phenotypes group-wise.
    - Otherwise, skip phenotype assignment for that row to prevent incorrect inflation.
    """
    if "phenotypes" not in df.columns:
        return {}

    from collections import defaultdict

    sample_phenos = defaultdict(set)

    for idx, row in df.iterrows():
        pheno_str = row.get("phenotypes", "")
        if not isinstance(pheno_str, str) or not pheno_str.strip():
            continue

        gt_val = row.get("GT", "")
        sample_entries = gt_val.split(";") if gt_val and isinstance(gt_val, str) else []
        sample_names = []
        for entry in sample_entries:
            entry = entry.strip()
            if not entry:
                continue
            sname, _ = extract_sample_and_genotype(entry)
            if sname:
                sample_names.append(sname)

        pheno_groups = [pg.strip() for pg in pheno_str.split(";") if pg.strip()]

        if len(pheno_groups) > 1:
            # Multiple phenotype groups
            if len(pheno_groups) == len(sample_names):
                # Perfect alignment
                for sname, pgroup in zip(sample_names, pheno_groups, strict=False):
                    phenos = {p.strip() for p in pgroup.split(",") if p.strip()}
                    sample_phenos[sname].update(phenos)
            else:
                # Mismatch in counts, skip assignment
                logger.warning(
                    "Phenotype groups (%d) != sample count (%d) at row %d. "
                    "Skipping phenotype assignment for this row.",
                    len(pheno_groups),
                    len(sample_names),
                    idx,
                )
                # No assignment here to prevent incorrect inflation
        else:
            # Only one phenotype group (or none)
            phenos = {p.strip() for p in pheno_str.split(",") if p.strip()}
            for sname in sample_names:
                sample_phenos[sname].update(phenos)

    return dict(sample_phenos)


def assign_case_control_counts(
    df: pd.DataFrame,
    case_samples: set[str],
    control_samples: set[str],
    all_samples: set[str],
) -> pd.DataFrame:
    """
    Assign case/control counts, allele counts, and homozygous variant counts per variant.

    Creates columns:
    - proband_count/control_count: total number of case/control samples
    - proband_variant_count/control_variant_count: number of case/control samples with a variant allele
    - proband_allele_count/control_allele_count: sum of variant alleles in case/control samples
    - proband_homozygous_count/control_homozygous_count: how many case/control samples have a homozygous variant (1/1)

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of variants with a "GT" column listing variants per sample.
    case_samples : set of str
        Set of samples classified as cases.
    control_samples : set of str
        Set of samples classified as controls.
    all_samples : set of str
        All samples present in the VCF.

    Returns
    -------
    pd.DataFrame
        DataFrame with assigned case/control counts and alleles, including homozygous counts.
    """
    logger.debug("Assigning case/control counts to variants...")

    total_proband = len(case_samples)
    total_control = len(control_samples)
    logger.debug(
        "assign_case_control_counts: Total proband samples: %d, Total control samples: %d",
        total_proband,
        total_control,
    )

    proband_variant_count_list = []
    control_variant_count_list = []
    proband_allele_count_list = []
    control_allele_count_list = []
    p_hom_count_list = []
    c_hom_count_list = []

    for idx, val in enumerate(df["GT"]):
        # Log progress every 1000 variants
        if idx % 1000 == 0:
            logger.debug("Processing variant row %d for allele counts...", idx)

        samples_with_variant = {}
        if isinstance(val, str) and val.strip():
            for s in val.split(";"):
                s = s.strip()
                if not s:
                    continue
                sample_name, genotype = extract_sample_and_genotype(s)
                if sample_name is not None:
                    samples_with_variant[sample_name] = genotype

        p_variant_count = 0
        c_variant_count = 0
        p_allele_count = 0
        c_allele_count = 0
        p_hom_count = 0
        c_hom_count = 0

        # For each sample, if not in samples_with_variant, assume genotype=0/0
        for sample_name in all_samples:
            genotype = samples_with_variant.get(sample_name, "0/0")
            allele_count = genotype_to_allele_count(genotype)

            # Check if homozygous variant
            is_hom_variant = genotype == "1/1"

            if sample_name in case_samples:
                if allele_count > 0:
                    p_variant_count += 1
                    p_allele_count += allele_count
                    if is_hom_variant:
                        p_hom_count += 1
            elif sample_name in control_samples and allele_count > 0:
                c_variant_count += 1
                c_allele_count += allele_count
                if is_hom_variant:
                    c_hom_count += 1

        proband_variant_count_list.append(p_variant_count)
        control_variant_count_list.append(c_variant_count)
        proband_allele_count_list.append(p_allele_count)
        control_allele_count_list.append(c_allele_count)
        p_hom_count_list.append(p_hom_count)
        c_hom_count_list.append(c_hom_count)

    df["proband_count"] = total_proband
    df["control_count"] = total_control
    df["proband_variant_count"] = proband_variant_count_list
    df["control_variant_count"] = control_variant_count_list
    df["proband_allele_count"] = proband_allele_count_list
    df["control_allele_count"] = control_allele_count_list
    df["proband_homozygous_count"] = p_hom_count_list
    df["control_homozygous_count"] = c_hom_count_list

    logger.debug("Case/control counts assigned.")
    if len(df) > 0:
        logger.debug(
            "Example for first row: %s",
            df.iloc[0][
                [
                    "proband_count",
                    "proband_variant_count",
                    "proband_allele_count",
                    "proband_homozygous_count",
                    "control_count",
                    "control_variant_count",
                    "control_allele_count",
                    "control_homozygous_count",
                ]
            ].to_dict(),
        )

    return df


def extract_sample_and_genotype(sample_field: str) -> tuple[str, str]:
    """
    Extract sample name and genotype from a field like 'sample(0/1:172:110,62)' or 'sample(0/1)'.

    If parentheses are missing, assume no genotype is specified -> no variant (0/0).

    The genotype portion may include extra fields after a colon (e.g. coverage),
    so we split on the first colon to isolate the actual genotype string (e.g., '0/1').

    Parameters
    ----------
    sample_field : str
        A string like 'sample(0/1)', 'sample(0/1:172:110,62)', or 'sample'.

    Returns
    -------
    (str, str)
        (sample_name, genotype) where genotype is typically '0/1', '1/1', '0/0', etc.
    """
    start_idx = sample_field.find("(")
    end_idx = sample_field.find(")")
    # Default is the entire string as the sample name, and genotype empty
    sample_name = sample_field.strip()
    genotype = ""

    if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
        sample_name = sample_field[:start_idx].strip()
        genotype_str = sample_field[start_idx + 1 : end_idx].strip()

        # Split on the first colon to separate e.g. '0/1' from '172:110,62'
        genotype = genotype_str.split(":", 1)[0] if ":" in genotype_str else genotype_str

    return sample_name, genotype


def genotype_to_allele_count(genotype: str) -> int:
    """
    Convert genotype string to allele count.

    - '1/1' -> 2
    - '0/1' or '1/0' -> 1
    - '0/0' or '' -> 0

    Parameters
    ----------
    genotype : str
        Genotype string, expected to be one of '0/0', '0/1', '1/0', '1/1', or ''.

    Returns
    -------
    int
        The allele count for the given genotype.
    """
    if genotype == "1/1":
        return 2
    elif genotype in ["0/1", "1/0"]:
        return 1
    return 0


def load_gene_list(file_path: str) -> set[str]:
    """
    Load genes from a file (one gene per line) into a set of uppercase gene names.

    Parameters
    ----------
    file_path : str
        Path to a file containing gene names, one per line

    Returns
    -------
    Set[str]
        A set of gene names in uppercase for case-insensitive matching
    """
    genes: set[str] = set()
    if not os.path.exists(file_path):
        logger.warning(f"Gene list file not found: {file_path}. Skipping this list.")
        return genes

    try:
        with open(file_path, encoding="utf-8") as f:
            for line_num, line_content in enumerate(f, 1):
                gene = line_content.strip()
                if gene:  # Avoid empty lines
                    if (
                        not gene.isalnum() and "_" not in gene and "-" not in gene
                    ):  # Basic check for valid gene name characters
                        logger.debug(
                            f"Potentially non-standard character in gene '{gene}' from file {file_path} at line {line_num}. Including anyway."
                        )
                    genes.add(gene.upper())  # Store as uppercase for case-insensitive matching
    except Exception as e:
        logger.error(f"Error reading gene list file {file_path}: {e}")

    if not genes:
        logger.warning(f"Gene list file {file_path} was empty or no valid genes found.")
    else:
        logger.debug(f"Loaded {len(genes)} genes from {file_path}")

    return genes


def _sanitize_column_name(file_path: str) -> str:
    """
    Convert a file path to a valid TSV column name.

    Extract the filename and remove invalid characters.

    Parameters
    ----------
    file_path : str
        Path to a file

    Returns
    -------
    str
        A sanitized column name derived from the file's basename
    """
    base_name = os.path.basename(file_path)
    name_part, ext = os.path.splitext(base_name)

    # Iteratively remove common extensions
    common_extensions = [".txt", ".list", ".tsv", ".csv", ".genes", ".gene"]
    original_name_part = name_part
    while True:
        found_ext = False
        for common_ext in common_extensions:
            if name_part.lower().endswith(common_ext):
                name_part = name_part[: -len(common_ext)]
                found_ext = True
                break
        if not found_ext:
            break

    if not name_part:  # If all was extension
        name_part = original_name_part  # Fallback to name before extension stripping

    # Replace non-alphanumeric characters (except underscore) with underscore
    sanitized_name = re.sub(r"[^\w_]", "_", name_part)
    # Remove leading/trailing underscores and collapse multiple underscores
    sanitized_name = re.sub(r"_+", "_", sanitized_name).strip("_")

    if not sanitized_name:  # If name becomes empty after sanitization (e.g. "---.txt")
        sanitized_name = "unnamed_list"

    # Ensure it doesn't start with a number or is purely numeric
    if not sanitized_name or sanitized_name[0].isdigit() or sanitized_name.isdigit():
        sanitized_name = "list_" + sanitized_name

    return sanitized_name


def annotate_variants_with_gene_lists(lines: list[str], gene_list_files: list[str]) -> list[str]:
    """
    Add new columns to variant data lines indicating membership in provided gene lists.

    For each gene list file, a new column is added to the TSV file. The column will contain 'yes' if
    any of the genes in the GENE column (which may be comma-separated) is present in the gene list,
    and 'no' otherwise.

    Parameters
    ----------
    lines : List[str]
        Lines of the TSV file to annotate
    gene_list_files : List[str]
        List of file paths containing gene names (one per line)

    Returns
    -------
    List[str]
        The input lines with additional columns for each gene list

    Notes
    -----
    - If the GENE column is missing, an error is logged and the input lines are returned unchanged
    - Multiple genes in the GENE column can be separated by commas, semicolons, or spaces
    - If a gene list file fails to load, it is skipped with an error message
    - Column names are sanitized from the file basename, duplicates are suffixed with numbers
    """
    if not gene_list_files or not lines or len(lines) < 2:
        return lines

    loaded_gene_sets: dict[str, set[str]] = {}
    ordered_new_column_names: list[str] = []

    for file_path in gene_list_files:
        genes = load_gene_list(file_path)
        if genes:  # Only process if genes were successfully loaded
            base_col_name = _sanitize_column_name(file_path)
            final_col_name = base_col_name
            suffix = 1
            # Ensure unique column names if multiple files sanitize to the same name
            while final_col_name in loaded_gene_sets:
                final_col_name = f"{base_col_name}_{suffix}"
                suffix += 1

            loaded_gene_sets[final_col_name] = genes
            ordered_new_column_names.append(final_col_name)
            logger.info(
                f"Loaded {len(genes)} genes for annotation column '{final_col_name}' from {file_path}"
            )

    if not loaded_gene_sets:
        logger.warning("No valid gene lists were loaded; no annotation columns will be added.")
        return lines

    # Process header
    header_line = lines[0].rstrip("\n")
    header_parts = header_line.split("\t")

    try:
        gene_col_idx = header_parts.index("GENE")
    except ValueError:
        logger.error(
            "CRITICAL: 'GENE' column not found in TSV header. Cannot perform gene list annotation."
        )
        # Return lines unmodified to prevent pipeline breakage if GENE column is unexpectedly missing
        return lines

    new_header_string = "\t".join(header_parts + ordered_new_column_names)
    modified_lines = [new_header_string]

    # Process data lines
    for data_line_str in lines[1:]:
        data_line_str = data_line_str.rstrip("\n")
        if not data_line_str.strip():  # Preserve empty lines if any
            modified_lines.append(data_line_str)
            continue

        fields = data_line_str.split("\t")

        # Ensure fields list is long enough for gene_col_idx
        if gene_col_idx >= len(fields):
            logger.warning(
                f"Line has too few fields to access GENE column (index {gene_col_idx}): {data_line_str}"
            )
            # Append "no" for all new columns for this malformed line
            fields.extend(["" for _ in range(gene_col_idx - len(fields) + 1)])
            fields.extend(["no" for _ in ordered_new_column_names])
            modified_lines.append("\t".join(fields))
            continue

        variant_gene_field_value = fields[gene_col_idx].strip().upper()  # Match in uppercase

        # Handle if GENE column contains multiple genes (e.g., comma-separated or other delimiters)
        if variant_gene_field_value:
            variant_genes_in_record = {
                g.strip() for g in re.split(r"[,; ]+", variant_gene_field_value) if g.strip()
            }
        else:
            variant_genes_in_record = set()

        for col_name in ordered_new_column_names:
            gene_set_for_current_list = loaded_gene_sets[col_name]
            is_member = "no"
            # Check for intersection between variant genes and gene list
            if variant_genes_in_record and not variant_genes_in_record.isdisjoint(
                gene_set_for_current_list
            ):
                is_member = "yes"
            fields.append(is_member)

        modified_lines.append("\t".join(fields))

    logger.info(
        f"Annotation with gene lists complete. Added {len(ordered_new_column_names)} columns."
    )
    return modified_lines


def dump_df_to_xlsx(
    df: pd.DataFrame, output_path: str | Path, sheet_name: str = "Sheet1", index: bool = False
) -> None:
    """
    Export a pandas DataFrame to an Excel (xlsx) file.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to export
    output_path : str or Path
        Path where the Excel file will be saved
    sheet_name : str, optional
        Name of the worksheet in the Excel file, by default "Sheet1"
    index : bool, optional
        Whether to include the DataFrame's index in the Excel file, by default False

    Returns
    -------
    None
        The function writes the DataFrame to an Excel file and logs the action
    """
    try:
        output_path = str(output_path)  # Convert Path to string if needed
        df.to_excel(output_path, sheet_name=sheet_name, index=index)
        logger.info(f"DataFrame exported to Excel file: {output_path}")
    except Exception as e:
        logger.error(f"Error exporting DataFrame to Excel: {e}")
        raise


def extract_gencode_id(gene_name: str) -> str:
    """
    Extract the GENCODE ID (ENSG) from a gene name string if present.

    Parameters
    ----------
    gene_name : str
        Gene name, possibly containing an ENSG ID

    Returns
    -------
    str
        The ENSG ID if found, otherwise an empty string
    """
    match = re.search(r"(ENSG\d+)", gene_name)
    return match.group(1) if match else ""


def get_vcf_names(vcf_file: str) -> list[str]:
    """
    Get sample names from a VCF file header.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    List[str]
        List of sample names from the VCF file
    """
    try:
        from subprocess import PIPE, Popen

        cmd = ["bcftools", "query", "-l", vcf_file]
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.error(f"Error getting sample names from VCF: {stderr}")
            return []

        return [line.strip() for line in stdout.split("\n") if line.strip()]
    except Exception as e:
        logger.error(f"Error running bcftools query: {e}")
        return []


def get_vcf_regions(vcf_file: str) -> list[str]:
    """
    Get the chromosomal regions covered in a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    List[str]
        List of regions in format "chr:start-end"
    """
    try:
        from subprocess import PIPE, Popen

        cmd = ["bcftools", "query", "-f", "%CHROM:%POS\n", vcf_file]
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.error(f"Error getting regions from VCF: {stderr}")
            return []

        positions = [line.strip() for line in stdout.split("\n") if line.strip()]
        if not positions:
            return []

        # Group by chromosome and find min/max positions
        by_chrom = defaultdict(list)
        for pos in positions:
            if ":" in pos:
                chrom, position = pos.split(":", 1)
                by_chrom[chrom].append(int(position))

        return [
            f"{chrom}:{min(positions)}-{max(positions)}" for chrom, positions in by_chrom.items()
        ]
    except Exception as e:
        logger.error(f"Error analyzing VCF regions: {e}")
        return []


def get_vcf_samples(vcf_file: str) -> list[str]:
    """
    Get an ordered list of sample names from a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    List[str]
        Ordered list of sample names (preserves VCF header order)
    """
    return get_vcf_names(vcf_file)


def get_vcf_size(vcf_file: str) -> int:
    """
    Get the number of variants in a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    int
        Number of variants in the VCF file
    """
    try:
        from subprocess import PIPE, Popen

        cmd = ["bcftools", "view", "-H", vcf_file, "|wc", "-l"]
        process = Popen(
            " ".join(cmd), shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.error(f"Error getting VCF size: {stderr}")
            return 0

        try:
            return int(stdout.strip())
        except ValueError:
            logger.error(f"Could not parse VCF size output: {stdout}")
            return 0
    except Exception as e:
        logger.error(f"Error counting variants in VCF: {e}")
        return 0


def match_igv_link_columns(header_fields: list[str]) -> dict[str, int]:
    """
    Find columns in a TSV header that should contain IGV links.

    Parameters
    ----------
    header_fields : List[str]
        List of header field names

    Returns
    -------
    Dict[str, int]
        Mapping of column name to column index for IGV link columns
    """
    # Columns that might contain coordinates for IGV links
    igv_column_patterns = ["CHROM", "POS", "START", "END", "POSITION", "LOCATION"]

    result = {}
    for i, field in enumerate(header_fields):
        if any(pattern in field.upper() for pattern in igv_column_patterns):
            result[field] = i

    return result


def read_sequencing_manifest(file_path: str) -> dict[str, dict[str, str]]:
    """
    Read a sequencing manifest file containing sample metadata.

    Expected format: TSV with header row and sample IDs in first column.

    Parameters
    ----------
    file_path : str
        Path to the manifest file

    Returns
    -------
    Dict[str, Dict[str, str]]
        Mapping of sample IDs to metadata dictionaries
    """
    result: dict[str, dict[str, str]] = {}
    if not os.path.exists(file_path):
        logger.warning(f"Manifest file not found: {file_path}")
        return result

    try:
        df = pd.read_csv(file_path, sep="\t", dtype=str)
        if df.empty:
            logger.warning(f"Empty manifest file: {file_path}")
            return result

        # Assuming first column contains sample IDs
        id_col = df.columns[0]

        # Convert each row to a dictionary
        for _, row in df.iterrows():
            sample_id = row[id_col]
            if sample_id and not pd.isna(sample_id):
                result[str(sample_id)] = row.to_dict()

        logger.info(f"Loaded metadata for {len(result)} samples from {file_path}")
        return result
    except Exception as e:
        logger.error(f"Error reading sequencing manifest: {e}")
        return result
