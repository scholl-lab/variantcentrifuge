# File: variantcentrifuge/reporting/generate_igv_report.py

import os
import subprocess
import re
import json
import logging
from typing import Dict

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


def generate_igv_report(
    variants_tsv: str,
    output_dir: str,
    bam_mapping_file: str,
    igv_reference: str,
    integrate_into_main: bool = False
) -> None:
    """
    Generate per-variant per-sample IGV reports using igv-reports 'create_report'.

    For each variant present in one or more samples (genotype not 0/0 or ./.), we:
    - Create a one-line TSV describing that variant with columns: CHROM, POS, REF, ALT.
    - For each sample carrying that variant, run create_report:
      create_report single_variant.tsv
        --genome igv_reference
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

    Args:
        variants_tsv (str): Path to the final variants TSV file.
        output_dir (str): Output directory where the main reports are stored.
                          IGV reports will be placed in output_dir/report/igv/.
        bam_mapping_file (str): Path to BAM mapping file.
        igv_reference (str): Genome reference (e.g., hg19 or hg38) for IGV.
        integrate_into_main (bool): If True, integrate into main report (placeholder).
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

    processed_variants = 0
    variant_report_map = []  # will store JSON mapping entries

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

                            variant_tsv_path = os.path.join(
                                igv_dir,
                                f"{sample_id}_{chrom}_{pos}_{ref_allele}_{alt_allele}_variant.tsv"
                            )
                            with open(variant_tsv_path, "w", encoding="utf-8") as sf:
                                # Columns: CHROM POS REF ALT
                                sf.write("CHROM\tPOS\tREF\tALT\n")
                                sf.write(f"{chrom}\t{pos}\t{ref_allele}\t{alt_allele}\n")

                            sample_report_path = os.path.join(
                                igv_dir,
                                f"{sample_id}_{chrom}_{pos}_{ref_allele}_{alt_allele}_igv_report.html"
                            )

                            cmd = [
                                "create_report",
                                variant_tsv_path,
                                "--genome", igv_reference,
                                "--sequence", "1", "--begin", "2", "--end", "2",
                                "--flanking", "50",
                                "--info-columns", "CHROM", "POS", "REF", "ALT",
                                "--tracks", bam_mapping[sample_id],
                                "--output", sample_report_path
                            ]

                            logger.debug(f"Running create_report for sample {sample_id}, variant {chrom}:{pos} {ref_allele}>{alt_allele}")
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
                                variant_report_map.append({
                                    "sample_id": sample_id,
                                    "chrom": chrom,
                                    "pos": pos,
                                    "ref": ref_allele,
                                    "alt": alt_allele,
                                    "report_path": os.path.relpath(sample_report_path, output_dir)
                                })

                            processed_variants += 1
                            logger.info(f"Progress: {processed_variants}/{total_variants} variants processed.")

    # Write JSON mapping
    igv_map_file = os.path.join(igv_dir, "igv_reports_map.json")
    with open(igv_map_file, "w", encoding="utf-8") as jf:
        json.dump(variant_report_map, jf, indent=2)
    logger.info(f"IGV report mapping JSON created at {igv_map_file}")

    logger.info("IGV report generation completed.")
