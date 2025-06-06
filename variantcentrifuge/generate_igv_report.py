# File: variantcentrifuge/reporting/generate_igv_report.py

import json
import logging
import os
import re
import subprocess
from typing import Dict

from .utils import generate_igv_safe_filename_base

logger = logging.getLogger("variantcentrifuge")


def parse_bam_mapping(bam_mapping_file: str) -> Dict[str, str]:
    """Parse a BAM mapping file (TSV or CSV) that maps sample_id to BAM file paths."""
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

    Args:
        variants_tsv (str): Path to the final variants TSV file.
        output_dir (str): Output directory where the main reports are stored.
                          IGV reports will be placed in output_dir/report/igv/.
        bam_mapping_file (str): Path to BAM mapping file.
        igv_reference (str, optional): Genome reference (e.g., hg19 or hg38) for IGV.
                                      Not used if igv_fasta is provided.
        integrate_into_main (bool): If True, integrate into main report (placeholder).
        igv_fasta (str, optional): Path to a local FASTA file for IGV reports. The index file (.fai)
                                   must exist in the same directory with the same name + '.fai' extension
                                   following the standard convention (e.g., genome.fa.fai).
        igv_ideogram (str, optional): Path to an ideogram file for IGV visualization.
        igv_max_allele_len_filename (int, optional): Maximum length for REF/ALT alleles in filenames.
                                                   Default is 10.
        igv_hash_len_filename (int, optional): Length of hash to append when truncating alleles.
                                             Default is 6.
        igv_max_variant_part_filename (int, optional): Maximum length for the variant part of filenames.
                                                     Default is 50.
        igv_flanking (int, optional): Flanking region size in base pairs for IGV reports.
                                    Default is 50.
    # MODIFIED: End of IGV flanking feature
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
    processed_variants = 0  # Initialize counter for processed variants

    # Now process variants again
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

                            # MODIFIED: Use filename shortening for variant TSV files
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
                            with open(variant_tsv_path, "w", encoding="utf-8") as sf:
                                # Columns: CHROM POS REF ALT
                                sf.write("CHROM\tPOS\tREF\tALT\n")
                                # Write the original full REF/ALT alleles to the TSV for proper visualization
                                sf.write(f"{chrom}\t{pos}\t{ref_allele}\t{alt_allele}\n")

                            sample_report_path = os.path.join(
                                igv_dir,
                                f"{safe_filename_base}_igv_report.html",
                            )

                            # MODIFIED: Start of local IGV FASTA feature
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
                                # MODIFIED: Start of IGV flanking feature - use configurable value
                                str(igv_flanking),
                                # MODIFIED: End of IGV flanking feature
                                "--info-columns",
                                "CHROM",
                                "POS",
                                "REF",
                                "ALT",
                                "--tracks",
                                bam_mapping[sample_id],
                                "--output",
                                sample_report_path,
                            ]
                            # If local FASTA is provided, use it instead of genome reference
                            if igv_fasta:
                                cmd.extend(["--fasta", igv_fasta])
                                # Optional ideogram file
                                if igv_ideogram:
                                    cmd.extend(["--ideogram", igv_ideogram])
                            else:
                                # Use genome reference if no local FASTA is provided
                                cmd.extend(["--genome", igv_reference])
                            # MODIFIED: End of local IGV FASTA feature

                            logger.debug(
                                f"Running create_report for sample {sample_id}, variant {chrom}:{pos} {ref_allele}>{alt_allele}"
                            )
                            result = subprocess.run(cmd, capture_output=True, text=True)

                            # Log output from create_report
                            if result.stdout:
                                logger.info(f"create_report stdout: {result.stdout.strip()}")
                            if result.stderr:
                                logger.warning(f"create_report stderr: {result.stderr.strip()}")

                            if result.returncode != 0:
                                logger.error(f"create_report failed for {sample_report_path}")
                            else:
                                logger.info(f"IGV report generated: {sample_report_path}")
                                # Add entry to variant report map
                                # Create a unique key for this variant
                                variant_key = f"{chrom}_{pos}_{ref_allele}_{alt_allele}"

                                # Get the relative path to the report
                                rel_report_path = os.path.relpath(sample_report_path, output_dir)

                                # If this variant is already in our map, add this sample's report
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
                            logger.info(
                                f"Progress: {processed_variants}/{total_variants} variants processed."
                            )

    # Write JSON mapping
    igv_map_file = os.path.join(igv_dir, "igv_reports_map.json")
    with open(igv_map_file, "w", encoding="utf-8") as jf:
        # Wrap the variant report map in a dictionary with 'variants' key
        json.dump({"variants": variant_report_map}, jf, indent=2)
    logger.info(f"IGV report mapping JSON created at {igv_map_file}")

    logger.info("IGV report generation completed.")
