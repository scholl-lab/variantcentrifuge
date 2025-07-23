"""IGV report generation for VariantCentrifuge."""

import json
import logging
import os
import re
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Tuple

from .utils import generate_igv_safe_filename_base

logger = logging.getLogger("variantcentrifuge")


def _generate_single_igv_report(
    task_data: Tuple[str, str, str, str, str, str, str, str, str, int, str, str],
) -> Tuple[bool, str, str, str, str, str, str]:
    """
    Generate a single IGV report for one variant/sample combination.

    Parameters
    ----------
    task_data : tuple
        Tuple containing all parameters needed for one IGV report generation:
        (sample_id, chrom, pos, ref_allele, alt_allele, bam_path, variant_tsv_path,
         sample_report_path, igv_flanking, igv_reference, igv_fasta, igv_ideogram)

    Returns
    -------
    tuple
        (success, sample_id, chrom, pos, ref_allele, alt_allele, sample_report_path)
    """
    (
        sample_id,
        chrom,
        pos,
        ref_allele,
        alt_allele,
        bam_path,
        variant_tsv_path,
        sample_report_path,
        igv_flanking,
        igv_reference,
        igv_fasta,
        igv_ideogram,
    ) = task_data

    try:
        # Create variant TSV file for this specific report
        with open(variant_tsv_path, "w", encoding="utf-8") as sf:
            sf.write("CHROM\tPOS\tREF\tALT\n")
            sf.write(f"{chrom}\t{pos}\t{ref_allele}\t{alt_allele}\n")

        # Build create_report command
        cmd = [
            "create_report",
            variant_tsv_path,
            "--sequence",
            "1",
            "--begin",
            "2",
            "--end",
            "2",
            "--flanking",
            str(igv_flanking),
            "--info-columns",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "--tracks",
            bam_path,
            "--output",
            sample_report_path,
        ]

        # Add genome reference or local FASTA
        if igv_fasta:
            cmd.extend(["--fasta", igv_fasta])
            if igv_ideogram:
                cmd.extend(["--ideogram", igv_ideogram])
        elif igv_reference:
            cmd.extend(["--genome", igv_reference])

        logger.debug(
            f"Running create_report for sample {sample_id}, variant {chrom}:{pos} {ref_allele}>{alt_allele}"
        )

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"create_report failed for {sample_report_path}: {result.stderr}")
            return False, sample_id, chrom, pos, ref_allele, alt_allele, sample_report_path
        else:
            logger.info(f"IGV report generated: {sample_report_path}")
            return True, sample_id, chrom, pos, ref_allele, alt_allele, sample_report_path

    except Exception as e:
        logger.error(f"Error generating IGV report for {sample_id} {chrom}:{pos}: {e}")
        return False, sample_id, chrom, pos, ref_allele, alt_allele, sample_report_path


def parse_bam_mapping(bam_mapping_file: str) -> Dict[str, str]:
    """
    Parse a BAM mapping file (TSV or CSV) that maps sample_id to BAM file paths.

    Parameters
    ----------
    bam_mapping_file : str
        Path to the BAM mapping file.

    Returns
    -------
    dict
        Dictionary mapping sample IDs to BAM file paths.
    """
    mapping = {}
    with open(bam_mapping_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", "\t").split("\t")
            if len(parts) < 2:
                continue
            sample_id = parts[0].strip()
            bam_path = parts[1].strip()
            if sample_id and bam_path:
                mapping[sample_id] = bam_path
    return mapping


# MODIFIED: Start of local IGV FASTA feature
def generate_igv_report(
    variants_tsv: str,
    output_dir: str,
    bam_mapping_file: str,
    igv_reference: str = None,
    integrate_into_main: bool = False,
    igv_fasta: str = None,
    igv_ideogram: str = None,
    igv_max_allele_len_filename: int = 10,
    igv_hash_len_filename: int = 6,
    igv_max_variant_part_filename: int = 50,
    # MODIFIED: Start of IGV flanking feature
    igv_flanking: int = 50,
    max_workers: int = 4,
) -> None:
    # MODIFIED: End of local IGV FASTA feature
    """
    Generate per-variant per-sample IGV reports using igv-reports 'create_report'.

    For each variant present in one or more samples (genotype not 0/0 or ./.), we:

    - Create a one-line TSV describing that variant with columns: CHROM, POS, REF, ALT.
    - For each sample carrying that variant, run create_report:
      create_report single_variant.tsv
        --genome igv_reference (if using an online genome reference)
        OR
        --fasta local_fasta_file (if using a local FASTA file)
        --sequence 1 --begin 2 --end 2
        --flanking 50
        --info-columns CHROM POS REF ALT
        --tracks corresponding_sample.bam
        --output sample_report.html
    - Produce an HTML report named {sample_id}_{CHROM}_{POS}_{REF}_{ALT}_igv_report.html.

    Changes:

    1. Place IGV reports into a subfolder report/igv/ within output_dir.
    2. Create a JSON mapping file (igv_reports_map.json) that maps each generated IGV report
       back to its sample and variant for integration in the interactive HTML report.
    3. Support for local FASTA files to avoid network access requirements.
    4. Implement filename shortening for long REF/ALT alleles to prevent "File name too long" errors.
    5. Support for configurable flanking region size for IGV reports.

    Parameters
    ----------
    variants_tsv : str
        Path to the final variants TSV file.
    output_dir : str
        Output directory where the main reports are stored.
        IGV reports will be placed in output_dir/report/igv/.
    bam_mapping_file : str
        Path to BAM mapping file.
    igv_reference : str, optional
        Genome reference (e.g., hg19 or hg38) for IGV.
        Not used if igv_fasta is provided.
    integrate_into_main : bool
        If True, integrate into main report (placeholder).
    igv_fasta : str, optional
        Path to a local FASTA file for IGV reports. The index file (.fai)
        must exist in the same directory with the same name + '.fai' extension
        following the standard convention (e.g., genome.fa.fai).
    igv_ideogram : str, optional
        Path to an ideogram file for IGV visualization.
    igv_max_allele_len_filename : int, optional
        Maximum length for REF/ALT alleles in filenames.
        Default is 10.
    igv_hash_len_filename : int, optional
        Length of hash to append when truncating alleles.
        Default is 6.
    igv_max_variant_part_filename : int, optional
        Maximum length for the variant part of filenames.
        Default is 50.
    igv_flanking : int, optional
        Flanking region size in base pairs for IGV reports.
        Default is 50.
    max_workers : int, optional
        Maximum number of parallel workers for IGV report generation.
        Default is 4. Higher values can speed up generation but may overwhelm the system.

    Returns
    -------
    None
    """
    logger.info("Starting IGV report generation.")
    # Place IGV reports into report/igv/ subfolder
    igv_dir = os.path.join(output_dir, "igv")
    os.makedirs(igv_dir, exist_ok=True)

    bam_mapping = parse_bam_mapping(bam_mapping_file)
    logger.debug(f"BAM mapping loaded: {len(bam_mapping)} entries.")

    # Identify columns
    with open(variants_tsv, "r", encoding="utf-8") as f:
        header_line = f.readline().strip()
    header = header_line.split("\t")

    # Required columns
    required_cols = ["CHROM", "POS", "REF", "ALT", "GT"]
    for rc in required_cols:
        if rc not in header:
            logger.error(f"Missing required column '{rc}' in {variants_tsv}. Cannot proceed.")
            return

    chrom_idx = header.index("CHROM")
    pos_idx = header.index("POS")
    ref_idx = header.index("REF")
    alt_idx = header.index("ALT")
    gt_idx = header.index("GT")

    pattern = re.compile(r"([^()]+)\(([^)]+)\)")

    # Count total variants that will be processed
    total_variants = 0
    with open(variants_tsv, "r", encoding="utf-8") as f:
        f.readline()  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            gt_value = parts[gt_idx].strip()
            if gt_value:
                sample_entries = gt_value.split(";")
                for entry in sample_entries:
                    entry = entry.strip()
                    if not entry:
                        continue
                    m = pattern.match(entry)
                    if m:
                        genotype = m.group(2).strip()
                        if genotype not in ["0/0", "./."]:
                            total_variants += 1
    logger.info(f"Total variants to process for IGV reports: {total_variants}")

    # Initialize list to store mapping of variants to reports
    # The structure will be: [{ chrom, pos, ref, alt, sample_reports: {sample_id: report_path, ...} }, ...]
    variant_report_map = []
    # Dictionary to map variant identifiers (chrom_pos_ref_alt) to their index in variant_report_map
    variant_index_map = {}

    # Collect all tasks for parallel processing
    tasks = []

    # Collect all variant/sample combinations first
    with open(variants_tsv, "r", encoding="utf-8") as f:
        f.readline()  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            chrom = parts[chrom_idx]
            pos = parts[pos_idx]
            ref_allele = parts[ref_idx]
            alt_allele = parts[alt_idx]
            gt_value = parts[gt_idx].strip()

            if gt_value:
                sample_entries = gt_value.split(";")
                for entry in sample_entries:
                    entry = entry.strip()
                    if not entry:
                        continue
                    m = pattern.match(entry)
                    if m:
                        sample_id = m.group(1).strip()
                        genotype = m.group(2).strip()
                        # If the sample carries the variant (not 0/0 or ./.)
                        if genotype not in ["0/0", "./."]:
                            if sample_id not in bam_mapping:
                                logger.warning(f"No BAM found for sample {sample_id}, skipping.")
                                continue

                            # Generate safe filename
                            safe_filename_base = generate_igv_safe_filename_base(
                                sample_id,
                                chrom,
                                pos,
                                ref_allele,
                                alt_allele,
                                max_allele_len=igv_max_allele_len_filename,
                                hash_len=igv_hash_len_filename,
                                max_variant_part_len=igv_max_variant_part_filename,
                            )

                            variant_tsv_path = os.path.join(
                                igv_dir,
                                f"{safe_filename_base}_variant.tsv",
                            )

                            sample_report_path = os.path.join(
                                igv_dir,
                                f"{safe_filename_base}_igv_report.html",
                            )

                            # Create task tuple for parallel processing
                            task = (
                                sample_id,
                                chrom,
                                pos,
                                ref_allele,
                                alt_allele,
                                bam_mapping[sample_id],
                                variant_tsv_path,
                                sample_report_path,
                                igv_flanking,
                                igv_reference,
                                igv_fasta,
                                igv_ideogram,
                            )
                            tasks.append(task)

    logger.info(f"Starting parallel IGV report generation with {max_workers} workers")
    logger.info(f"Total tasks to process: {len(tasks)}")

    # Process tasks in parallel
    processed_variants = 0
    map_lock = threading.Lock()  # Lock for thread-safe map updates

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(_generate_single_igv_report, task): task for task in tasks
        }

        # Process completed tasks
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                success, sample_id, chrom, pos, ref_allele, alt_allele, sample_report_path = (
                    future.result()
                )

                if success:
                    # Add entry to variant report map
                    variant_key = f"{chrom}_{pos}_{ref_allele}_{alt_allele}"

                    # Get the relative path to the report
                    rel_report_path = os.path.relpath(sample_report_path, output_dir)

                    # Thread-safe update of variant map
                    with map_lock:
                        if variant_key in variant_index_map:
                            variant_idx = variant_index_map[variant_key]
                            variant_report_map[variant_idx]["sample_reports"][
                                sample_id
                            ] = rel_report_path
                        else:
                            # First time seeing this variant, create a new entry
                            new_variant = {
                                "chrom": chrom,
                                "pos": pos,
                                "ref": ref_allele,
                                "alt": alt_allele,
                                "sample_reports": {sample_id: rel_report_path},
                            }
                            variant_report_map.append(new_variant)
                            variant_index_map[variant_key] = len(variant_report_map) - 1

                processed_variants += 1
                if processed_variants % 10 == 0 or processed_variants == len(tasks):
                    logger.info(
                        f"Progress: {processed_variants}/{len(tasks)} IGV reports completed"
                    )

            except Exception as e:
                logger.error(f"Task failed: {e}")
                processed_variants += 1

    # Write JSON mapping
    igv_map_file = os.path.join(igv_dir, "igv_reports_map.json")
    with open(igv_map_file, "w", encoding="utf-8") as jf:
        # Wrap the variant report map in a dictionary with 'variants' key
        json.dump({"variants": variant_report_map}, jf, indent=2)
    logger.info(f"IGV report mapping JSON created at {igv_map_file}")

    logger.info("IGV report generation completed.")
